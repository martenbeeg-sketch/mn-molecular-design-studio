"""Families, clustering and diversity statistics (docs/families.md).

A family is one BindCraft backbone (one trajectory, the file stem up to _s{seed}). Families are counted,
sequences are not. Per run the module builds the family records, clusters the families at two levels
(binder fold by TM-align, binding site by contacted-residue overlap) and gives the Chao-Jost
"do I need more designs" readout and a replay of it. Estimates are per run; the pooled block is descriptive.

The statistics are a faithful port of Family.cluster, Family.chaoJost and Family.advice of the mockup
(mockups/build/web/family.js), because the viewer recomputes them in the browser and the two must agree at
the defaults (tests/test_families_js.py runs the real family.js under Node against this code). Where the
mockup breaks (coverage is NaN for f1 = 0 and f2 = 0 at n = 3) this module returns the correct limit and
the divergence is listed in docs/families.md.

Numbering: footprints are UniProt residue sets (TargetConstruct.to_uniprot); when a run has no
uniprot_offset yet the module falls back to construct numbers and says so in `warnings`.
"""

from __future__ import annotations

import base64
import hashlib
import json
import math
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Optional, Sequence

import numpy as np

from . import paratope as _para
from . import structure as S
from .config import filters_for_run
from .model import CONTACT_CUTOFF, DesignRecord, RunData, binder_chain_count

SCHEMA = 1
BASES = ("fold", "site", "both")
#: Average-linkage cut distance (1 - similarity) per basis: TM 0.5 is "same fold" (Foldseek convention used
#: by the plan), Jaccard 0.5 is "same site", the combined distance uses the same cut.
DEFAULT_DISTANCE = {"fold": 0.5, "site": 0.5, "both": 0.5}
FOLD_WEIGHT = 0.5
BASIS_DEFAULT = "both"
#: Cuts reported next to the defaults so a threshold artefact is visible (0.2 = TM 0.8 "same backbone").
SWEEP_DISTANCES = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
#: Forecast horizon of the mockup's text ("600 more trajectories").
FORECAST_HORIZON = 600
WILSON_Z80 = 1.2816
#: Gates (plan, "Families and diversity").
MIN_FAMILIES = 10
LOW_CONFIDENCE_BELOW = 30
COVERAGE_KEEP = 0.7
COVERAGE_ENOUGH = 0.9
#: Scaffolded runs: two loop lengths are the same length class when no more than this many residues apart (a class is a chain of
#: lengths each within this many residues of the next). A convention of this module, written next to the result.
LOOP_CLASS_GAP = 3
FOLD_OFF_WHY = "every design is built on the same scaffold, so fold similarity measures the scaffold, not the designs"
#: The bases that include fold distance, and the one that does not.
FOLD_BASES = ("fold", "both")
SITE_BASIS = "site"

#: Stems summarised per family (mean over the accepted siblings, CSV scale).
SUMMARY_STEMS = (
    "pLDDT", "pTM", "i_pTM", "i_pAE", "i_pLDDT", "ShapeComplementarity", "dG", "dSASA",
    "Surface_Hydrophobicity", "n_InterfaceHbonds", "n_InterfaceUnsatHbonds", "Binder_RMSD", "Hotspot_RMSD",
    "Binder_pLDDT",
)


# ------------------------------------------------------------------------------------------------ statistics
def _clamp(x: float, a: float, b: float) -> float:
    return a if x < a else b if x > b else x


def counts(sizes: Sequence[int]) -> tuple[int, int]:
    """(f1, f2): classes seen once and twice."""
    return sum(1 for s in sizes if s == 1), sum(1 for s in sizes if s == 2)


def hill(sizes: Sequence[int], q: float) -> float:
    n = sum(sizes)
    if not n:
        return 0.0
    p = [s / n for s in sizes]
    if q == 0:
        return float(len(sizes))
    if q == 1:
        return math.exp(-sum(x * math.log(x) for x in p if x))
    return sum(x ** q for x in p) ** (1.0 / (1.0 - q))


def coverage(sizes: Sequence[int]) -> float:
    """Chao-Jost sample coverage. Classes are clusters, individuals are families.

    C = 1 - (f1/n) (n-1) f1 / ((n-1) f1 + 2 f2); for f2 = 0 the bias-corrected form (n-1)(f1-1) and 2.
    f1 = 0 gives 1 (the mockup formula is NaN there for n = 3)."""
    n = sum(sizes)
    if not n:
        return 0.0
    f1, f2 = counts(sizes)
    if f1 == 0:
        return 1.0
    if f2 > 0:
        c = 1 - (f1 / n) * (((n - 1) * f1) / (((n - 1) * f1) + 2 * f2))
    else:
        c = 1 - (f1 / n) * (((n - 1) * (f1 - 1)) / (((n - 1) * (f1 - 1)) + 2))
    return _clamp(c, 0.0, 1.0)


@dataclass
class ChaoJost:
    n: int
    S: int
    f1: int
    f2: int
    C: float
    f0: float
    chao1: float

    def snew(self, m: float) -> float:
        """Expected new clusters after m more families (valid for m <= n)."""
        if self.f0 <= 0:
            return 0.0
        return self.f0 * (1 - (1 - self.f1 / (self.n * self.f0 + self.f1)) ** max(0.0, m))


def chao_jost(sizes: Sequence[int]) -> ChaoJost:
    n, S_ = sum(sizes), len(sizes)
    if not n:
        return ChaoJost(0, 0, 0, 0, 0.0, 0.0, 0.0)
    f1, f2 = counts(sizes)
    C = coverage(sizes)
    if n > 1:
        f0 = ((n - 1) / n) * ((f1 * f1) / (2 * f2) if f2 > 0 else (f1 * (f1 - 1)) / 2)
    else:
        f0 = 0.0
    return ChaoJost(n, S_, f1, f2, C, f0, S_ + f0)


def wilson(k: float, N: float, z: float = WILSON_Z80) -> tuple[float, float]:
    """Wilson score interval of k/N (z = 1.2816 gives an 80 percent range)."""
    if not N:
        return 0.0, 0.0
    p = k / N
    d = 1 + (z * z) / N
    c = p + (z * z) / (2 * N)
    r = z * math.sqrt((p * (1 - p)) / N + (z * z) / (4 * N * N))
    return (c - r) / d, (c + r) / d


def level_for(n_families: int, cov: float) -> str:
    """early (fewer than 10 families, yield only), keep (C < 0.7), maybe (C < 0.9), enough."""
    if n_families < MIN_FAMILIES:
        return "early"
    if cov < COVERAGE_KEEP:
        return "keep"
    if cov < COVERAGE_ENOUGH:
        return "maybe"
    return "enough"


def curve(order: Sequence[tuple[int, Any]]) -> list[dict[str, Any]]:
    """Accumulation curve in generation order. `order` is [(traj, cluster id)] sorted by traj."""
    pts: list[dict[str, Any]] = [{"t": 0, "fams": 0, "clusters": 0, "new": False}]
    seen: set = set()
    for i, (t, c) in enumerate(order):
        seen.add(c)
        pts.append({"t": t, "fams": i + 1, "clusters": len(seen), "new": pts[-1]["clusters"] < len(seen)})
    return pts


@dataclass
class Advice:
    """Port of Family.advice (family.js). `n_traj` is the number of trajectories in the units of the run's
    traj_basis (here: passed trajectories)."""

    coverage: float
    chao1: float
    f0: float
    f1: int
    f2: int
    S: int
    n_families: int
    n_traj: int
    p_new_family: float
    recent_new: int
    early_new: int
    expected_more: float
    fam_more: float
    yield_p: float
    yield_ci: tuple[float, float]
    low_conf: bool
    hill1: float
    hill2: float
    dominance: float
    level: str
    snew_list: list[float]
    cj: ChaoJost
    pts: list[dict[str, Any]]

    def at(self, t: float) -> float:
        """Expected clusters at trajectory t (t >= n_traj) at the observed yield."""
        return self.S + self.cj.snew(min(self.n_families, max(0.0, (t - self.n_traj) * self.yield_p)))

    def traj_to_families(self, k: float) -> float:
        return max(0.0, (k - self.n_families) / (self.yield_p or 1e-9))


def advice(order: Sequence[tuple[int, Any]], sizes: Sequence[int], n_traj: int, horizon: int = FORECAST_HORIZON) -> Advice:
    """Plain-language readout inputs (Family.advice): families in generation order with their cluster ids,
    the cluster sizes and the number of trajectories so far."""
    sizes = list(sizes)
    n, S_ = len(order), len(sizes)
    cj = chao_jost(sizes)
    pts = curve(order)
    yield_p = n / n_traj if n_traj else 0.0
    ci = wilson(n, n_traj)
    last_win = n_traj * 0.3
    recent = sum(1 for p in pts if p["t"] > n_traj - last_win and p["new"])
    early = sum(1 for p in pts if p["t"] <= last_win and p["new"])
    adv = Advice(
        coverage=cj.C, chao1=cj.chao1, f0=cj.f0, f1=cj.f1, f2=cj.f2, S=S_, n_families=n, n_traj=n_traj,
        p_new_family=_clamp(1 - cj.C, 0.0, 1.0), recent_new=recent, early_new=early, expected_more=0.0,
        fam_more=min(n, horizon * yield_p), yield_p=yield_p, yield_ci=ci, low_conf=n < LOW_CONFIDENCE_BELOW,
        hill1=hill(sizes, 1) if sizes else 0.0, hill2=hill(sizes, 2) if sizes else 0.0,
        dominance=(max(sizes) / n) if n and sizes else 0.0, level=level_for(n, cj.C),
        snew_list=[cj.snew(m) for m in range(1, n + 1)], cj=cj, pts=pts,
    )
    adv.expected_more = adv.at(n_traj + horizon) - S_
    return adv


# ------------------------------------------------------------------------------------------------ clustering
def average_linkage(D: np.ndarray, theta: float) -> list[list[int]]:
    """Average-linkage agglomerative clustering cut at distance theta (Family.cluster).

    Merges while the smallest average distance is <= theta (inclusive). Ties are broken like the mockup: the
    first pair in row-major order i < j. Clusters are returned largest first (stable), members in merge order."""
    D = np.array(D, dtype=np.float64, copy=True)
    n = len(D)
    if not n:
        return []
    members = [[i] for i in range(n)]
    alive = np.ones(n, dtype=bool)
    size = np.ones(n)
    upper = np.triu(np.ones((n, n), dtype=bool), 1)
    live = n
    while live > 1:
        M = np.where(upper & alive[:, None] & alive[None, :], D, np.inf)
        flat = int(np.argmin(M))
        bi, bj = divmod(flat, n)
        best = M[bi, bj]
        if not np.isfinite(best) or best > theta:
            break
        ks = np.flatnonzero(alive)
        ks = ks[(ks != bi) & (ks != bj)]
        d = (size[bi] * D[bi, ks] + size[bj] * D[bj, ks]) / (size[bi] + size[bj])
        D[bi, ks] = d
        D[ks, bi] = d
        members[bi] = members[bi] + members[bj]
        size[bi] += size[bj]
        alive[bj] = False
        live -= 1
    out = [members[i] for i in range(n) if alive[i]]
    return sorted(out, key=lambda m: -len(m))


def _jaccard_distance(a: frozenset, b: frozenset) -> float:
    i = len(a & b)
    u = len(a) + len(b) - i
    return 1 - i / u if u else 1.0


# ------------------------------------------------------------------------------------------------ geometry
@dataclass
class DesignGeometry:
    """What the family code needs from one complex: the binder CA trace (for TM-align), its sequence (only
    used to seed the alignment, never written out) and the contacted target residues in construct numbers."""

    run_id: str
    design_id: str
    ca: np.ndarray
    seq: str
    contacts: tuple[int, ...]
    structure_model: Optional[int] = None

    @property
    def chain(self) -> tuple[np.ndarray, str]:
        return self.ca, self.seq


def _ca_trace(binder: S.Atoms) -> tuple[np.ndarray, str]:
    ca = binder.ca_xyz()
    keep = ~np.isnan(ca[:, 0])
    seq = "".join(np.array(list(binder.sequence))[keep])
    return np.ascontiguousarray(ca[keep]), seq


def design_geometry(design: DesignRecord, cutoff: float = CONTACT_CUTOFF) -> DesignGeometry:
    """Load one design's structure. Contacts are heavy-atom pairs within `cutoff` (4.5 A, model.CONTACT_CUTOFF)."""
    if not design.structure_path:
        raise ValueError(f"{design.run_id}/{design.id}: no structure file")
    cx = S.load_complex(design.structure_path, n_binder_chains=binder_chain_count(design))
    ca, seq = _ca_trace(cx.binder)
    cont = S.residue_contacts(cx.target, cx.binder, cutoff)
    return DesignGeometry(design.run_id, design.id, ca, seq, tuple(int(n) for n in cont.a_residues), design.structure_model)


def load_geometry(designs: Iterable[DesignRecord], cutoff: float = CONTACT_CUTOFF) -> dict[tuple[str, str], DesignGeometry]:
    """{(run, design id): DesignGeometry} for designs that have a structure (about 20 ms each)."""
    return {d.key: design_geometry(d, cutoff) for d in designs if d.structure_path}


# ------------------------------------------------------------------------------------------------ numbering and epoch
def _numberer(run: RunData, uniprot_offset: Optional[int] = None) -> tuple[Callable[[int], int], str]:
    off = run.target.uniprot_offset if run.target.uniprot_offset is not None else uniprot_offset
    if off is None:
        return (lambda n: int(n)), "construct"
    return (lambda n, o=off: int(n) + o), "uniprot"


def _runs_ranges(numbers: Iterable[int]) -> list[list[int]]:
    out: list[list[int]] = []
    for n in sorted(numbers):
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return out


def epoch_components(run: RunData, uniprot_offset: Optional[int] = None) -> dict[str, Any]:
    """What defines a settings epoch: hotspots, target construct, length range and the filter set.
    A change in any of them starts a new epoch (estimates must not carry across)."""
    num, numbering = _numberer(run, uniprot_offset)
    st = run.settings
    lo, hi = st.length_min, st.length_max
    if lo is None and run.backbones:
        lo, hi = min(b.length for b in run.backbones), max(b.length for b in run.backbones)
    return {
        "numbering": numbering,
        "hotspots": sorted(num(h) for h in st.hotspots),
        "construct": _runs_ranges(num(n) for n in run.target.numbers),
        "length": [lo, hi],
        "filters": {
            "name": st.filters_name,
            "rules": [[r.stem, r.threshold, r.op, r.per_model] for r in filters_for_run(run).rules],   # the run's own filter set
        },
    }


def epoch_hash(run: RunData, uniprot_offset: Optional[int] = None) -> str:
    """Twelve hex characters over `epoch_components`. Equal hashes: estimates of the two runs may be pooled."""
    blob = json.dumps(epoch_components(run, uniprot_offset), sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(blob.encode()).hexdigest()[:12]


def epoch_differences(runs: Sequence[RunData], uniprot_offset: Optional[int] = None) -> list[str]:
    """Names of the epoch components in which the runs differ."""
    comps = [epoch_components(r, uniprot_offset) for r in runs]
    return [k for k in comps[0] if k != "numbering" and any(c[k] != comps[0][k] for c in comps[1:])] if comps else []


def _sequence_cap(run: RunData) -> tuple[int, str]:
    """Accepted sequences per backbone BindCraft allows (max_mpnn_sequences). Not recorded in many runs:
    the BindCraft default 2 is used when no backbone exceeds it, otherwise the observed maximum."""
    raw = run.settings.settings_raw or {}
    if "max_mpnn_sequences" in raw:
        return int(raw["max_mpnn_sequences"]), "settings"
    most = max((len(b.accepted_ids) for b in run.backbones), default=0)
    if most <= 2:
        return 2, "BindCraft default (not recorded in the run folder; no backbone has more than 2 accepted sequences)"
    return most, "inferred: highest number of accepted sequences on one backbone"


def traj_info(run: RunData) -> tuple[str, Optional[int]]:
    """(traj_basis, number of trajectories in that unit). Runs often only have the row order of
    trajectory_stats.csv, so the unit is 'passed' trajectories (rows) and the attempted total stays unknown."""
    n = run.funnel.trajectory_rows if run.funnel.trajectory_rows is not None else len(run.backbones)
    if not n:
        return "none", None
    return "passed", int(n)


# ------------------------------------------------------------------------------------------------ family records
def design_sort_key(d: DesignRecord) -> tuple:
    """Best first: BindCraft rank, then Average i_pTM, then id."""
    return (d.rank is None, d.rank if d.rank is not None else 0, -(d.metrics.get("i_pTM") or 0.0), d.id)


def quality_scores(run: RunData) -> dict[str, float]:
    """design id -> q in (0, 1] with the best design (by `design_sort_key`) highest within the run, ties
    impossible. The exporter writes this as design.q so that the viewer's representative of a family
    (largest q) is the same design as FamilyRecord.best."""
    acc = sorted(run.accepted(), key=design_sort_key)
    n = len(acc)
    return {d.id: round(1 - i / (n + 1), 6) for i, d in enumerate(acc)}


@dataclass
class FamilyRecord:
    """One backbone and what BindCraft kept from it."""

    id: str
    run_id: str
    traj: int                                  # 1-based row of trajectory_stats.csv (0 = unknown)
    length: int
    seed: int
    members: list[str]                         # every design with a row, accepted first (by rank), then rejected
    accepted: list[str]
    rejected: list[str]
    n_rejected_with_structure: int
    best: Optional[str]                        # best accepted design: the representative for fold distance
    ranks: list[int]
    footprint: list[int]                       # target residues touched by at least half of the accepted siblings
    footprint_numbering: str                   # "uniprot" or "construct"
    summary: dict[str, Optional[float]]        # mean of the accepted designs per stem
    best_metrics: dict[str, Optional[float]]
    trajectory_metrics: dict[str, Optional[float]]
    cum_time_h: Optional[float]                # sum of TrajectoryTime up to and including this trajectory
    at_sequence_cap: bool

    @property
    def n_accepted(self) -> int:
        return len(self.accepted)

    @property
    def n_rejected(self) -> int:
        return len(self.rejected)

    def to_dict(self) -> dict[str, Any]:
        d = dict(self.__dict__)
        d["n_accepted"], d["n_rejected"] = self.n_accepted, self.n_rejected
        return d


def build_families(
    run: RunData,
    geometry: Optional[dict[tuple[str, str], DesignGeometry]] = None,
    *,
    include_rejected_only: bool = False,
    uniprot_offset: Optional[int] = None,
    warnings: Optional[list[str]] = None,
) -> list[FamilyRecord]:
    """Family records of one run in generation order (trajectory_stats.csv row order).

    By default only backbones with at least one accepted design are families (rejected-only backbones are left out); rejected
    siblings are listed as members and counted. `include_rejected_only` also returns backbones whose every
    scored sequence was rejected. `geometry` (from `load_geometry`) avoids reloading structures; without it the
    accepted structures are read here."""
    warn = warnings if warnings is not None else []
    num, numbering = _numberer(run, uniprot_offset)
    if numbering == "construct":
        warn.append(f"{run.run_id}: uniprot_offset is not set, footprints use construct numbers")
    by_bb: dict[str, list[DesignRecord]] = {}
    for d in run.designs:
        by_bb.setdefault(d.backbone_id, []).append(d)
    bb_by_id = {b.id: b for b in run.backbones}
    unknown = sorted(set(by_bb) - set(bb_by_id))
    if unknown:
        warn.append(f"{run.run_id}: {len(unknown)} backbones have design rows but no trajectory row (traj 0)")
    cap, _ = _sequence_cap(run)
    cum, running = {}, 0.0
    for b in run.backbones:
        running += (b.trajectory_time_s or 0.0)
        cum[b.id] = running / 3600.0
    ids = sorted(by_bb, key=lambda k: (bb_by_id[k].order if k in bb_by_id else 10**9, k))
    out: list[FamilyRecord] = []
    for bid in ids:
        ds = by_bb[bid]
        acc = sorted((d for d in ds if d.status == "accepted"), key=design_sort_key)
        rej = sorted((d for d in ds if d.status != "accepted"), key=lambda d: (d.mpnn_index is None, d.mpnn_index or 0, d.id))
        if not acc and not include_rejected_only:
            continue
        b = bb_by_id.get(bid)
        # footprint: target residues touched by at least half of the accepted siblings (Family.group)
        cnt: dict[int, int] = {}
        usable = 0
        for d in acc:
            g = (geometry or {}).get(d.key)
            if g is None and d.structure_path:
                g = design_geometry(d)
            if g is None:
                continue
            usable += 1
            for r in g.contacts:
                cnt[r] = cnt.get(r, 0) + 1
        fp = sorted(num(r) for r, n in cnt.items() if n >= usable / 2) if usable else []
        summary: dict[str, Optional[float]] = {}
        for stem in SUMMARY_STEMS:
            vals = [d.metrics.get(stem) for d in acc if d.metrics.get(stem) is not None]
            summary[stem] = float(np.mean(vals)) if vals else None
        best = acc[0] if acc else None
        out.append(FamilyRecord(
            id=bid, run_id=run.run_id, traj=b.traj if b else 0,
            length=(b.length if b else ds[0].length), seed=(b.seed if b else ds[0].seed),
            members=[d.id for d in acc] + [d.id for d in rej], accepted=[d.id for d in acc], rejected=[d.id for d in rej],
            n_rejected_with_structure=sum(1 for d in rej if d.has_structure),
            best=best.id if best else None, ranks=[d.rank for d in acc if d.rank is not None],
            footprint=fp, footprint_numbering=numbering, summary=summary,
            best_metrics={s: best.metrics.get(s) for s in SUMMARY_STEMS} if best else {},
            trajectory_metrics={k: v for k, v in (b.trajectory_metrics.items() if b else []) if v is not None},
            cum_time_h=round(cum[bid], 3) if bid in cum else None,
            at_sequence_cap=len(acc) >= cap,
        ))
    return out


# ------------------------------------------------------------------------------------------------ fold matrix
def project_design_ids(runs: Sequence[RunData]) -> dict[tuple[str, str], str]:
    """{(run, design id): project design id}. The id itself when it is unique across the project (all runs, all
    designs), otherwise '<run>/<id>' (data_contract.md, design.id)."""
    seen: dict[str, int] = {}
    for r in runs:
        for d in r.designs:
            seen[d.id] = seen.get(d.id, 0) + 1
    return {d.key: (d.id if seen[d.id] == 1 else f"{d.run_id}/{d.id}") for r in runs for d in r.designs}


def fold_matrix(
    designs: Iterable[DesignRecord],
    geometry: Optional[dict[tuple[str, str], DesignGeometry]] = None,
    *,
    id_of: Optional[Callable[[DesignRecord], str]] = None,
    n_jobs: Optional[int] = None,
    progress: Optional[Callable[[int, int], None]] = None,
) -> dict[str, Any]:
    """All-against-all TM-align of the binder chains: the campaign blob's `fold`.

    `tm` is the base64 of an n x n uint8 matrix (row-major), round(100 x TM-score) normalised by the shorter
    chain, symmetric, diagonal 100. Rows follow the order of `designs`. `runs` and `names` are extra
    parallel lists (run id and the design id without any run prefix) so a row can be found by (run, design)."""
    ds = list(designs)
    geo = geometry or {}
    chains = []
    for d in ds:
        g = geo.get(d.key) or design_geometry(d)
        chains.append(g.chain)
    n = len(ds)
    M = np.zeros((n, n), dtype=np.uint8)
    if n:
        tm = S.tm_matrix(chains, n_jobs=n_jobs, progress=progress).tm
        M = np.array([[int(round(100 * float(v))) for v in row] for row in tm], dtype=np.int64)
        M = np.clip(M, 0, 100).astype(np.uint8)
        M = np.maximum(M, M.T)
        np.fill_diagonal(M, 100)
    ids = [(id_of(d) if id_of else d.id) for d in ds]
    if len(set(ids)) != len(ids):
        raise ValueError("fold_matrix: design ids are not unique; pass id_of=project_design_ids(...) lookups")
    return {
        "n": n, "ids": ids, "tm": base64.b64encode(M.tobytes()).decode("ascii"),
        "runs": [d.run_id for d in ds], "names": [d.id for d in ds],
        "note": ("TM-align (tmtools) on the binder CA trace, TM-score x 100 normalised by the shorter chain, over "
                 "the accepted designs of all runs. Backbones of different runs are compared too; clusters and "
                 "estimates are formed per run."),
    }


@dataclass
class FoldLookup:
    """Decoded fold matrix: similarity (TM x 100) by (run, design id)."""

    ids: list[str]
    M: np.ndarray                              # (n, n) uint8
    index: dict[tuple[Optional[str], str], int]

    def row(self, run_id: str, design_id: str) -> int:
        for key in ((run_id, design_id), (None, design_id), (None, f"{run_id}/{design_id}")):
            if key in self.index:
                return self.index[key]
        raise KeyError(f"design {run_id}/{design_id} is not in the fold matrix")

    def distance(self, a: tuple[str, str], b: tuple[str, str]) -> float:
        return 1 - float(self.M[self.row(*a), self.row(*b)]) / 100


def decode_fold(fold: dict[str, Any]) -> FoldLookup:
    n = int(fold["n"])
    M = np.frombuffer(base64.b64decode(fold["tm"]), dtype=np.uint8).reshape(n, n)
    idx: dict[tuple[Optional[str], str], int] = {}
    runs, names = fold.get("runs"), fold.get("names")
    for i, fid in enumerate(fold["ids"]):
        idx[(None, fid)] = i
        if runs and names:
            idx[(runs[i], names[i])] = i
    return FoldLookup(list(fold["ids"]), M, idx)


def fold_labels(fold: dict[str, Any], cut: float = DEFAULT_DISTANCE["fold"]) -> dict[str, str]:
    """{design id: 'F<n>'}: average-linkage clusters of all designs in the fold matrix cut at distance `cut`
    (TM 0.5 by default), numbered largest first. Pooled over runs and descriptive: design.sc of the pages."""
    fl = decode_fold(fold)
    D = 1 - fl.M.astype(np.float64) / 100
    labels: dict[str, str] = {}
    for ci, members in enumerate(average_linkage(D, cut), 1):
        for i in members:
            labels[fl.ids[i]] = f"F{ci}"
    return labels


@dataclass
class FoldReps:
    """Accepted designs with their labels and CA traces, for labelling designs that are not in the matrix
    (rejected designs get the label of their nearest accepted design)."""

    ids: list[str]
    labels: list[str]
    chains: list[tuple[np.ndarray, str]]
    cut: float = DEFAULT_DISTANCE["fold"]


def build_fold_reps(
    designs: Sequence[DesignRecord], fold: dict[str, Any], geometry: Optional[dict[tuple[str, str], DesignGeometry]] = None,
    cut: float = DEFAULT_DISTANCE["fold"],
) -> FoldReps:
    lab = fold_labels(fold, cut)
    geo = geometry or {}
    fl = decode_fold(fold)
    ids, labels, chains = [], [], []
    for d in designs:
        i = fl.row(d.run_id, d.id)
        ids.append(fl.ids[i])
        labels.append(lab[fl.ids[i]])
        chains.append((geo.get(d.key) or design_geometry(d)).chain)
    return FoldReps(ids, labels, chains, cut)


def _trace_of(x: Any) -> tuple[np.ndarray, str]:
    if isinstance(x, S.Atoms):
        return _ca_trace(x)
    if isinstance(x, DesignGeometry):
        return x.chain
    return np.ascontiguousarray(x[0], dtype=np.float64), x[1]


def fold_label_details(binder: Any, reps: FoldReps) -> tuple[str, float, Optional[str]]:
    """(label, TM x 100 of the nearest accepted design, its id). The label is that design's fold label when
    TM >= 1 - cut (0.5 by default), otherwise 'unclustered'. TM is normalised by the shorter chain."""
    xyz, seq = _trace_of(binder)
    best, bi = -1.0, -1
    for i, (cxyz, cseq) in enumerate(reps.chains):
        t = S.tm_align(xyz, seq, cxyz, cseq).tm_shorter
        if t > best:
            best, bi = t, i
    if bi < 0:
        return "unclustered", 0.0, None
    label = reps.labels[bi] if best >= (1 - reps.cut) - 1e-9 else "unclustered"
    return label, round(100 * best, 1), reps.ids[bi]


def fold_label_for(binder_atoms: Any, matrix_reps: FoldReps) -> str:
    """Fold label of a design that is not in the matrix (a rejected design): the label of its nearest accepted
    design by TM-align, or 'unclustered' when no accepted design reaches TM 1 - cut."""
    return fold_label_details(binder_atoms, matrix_reps)[0]


_REF_STATE: dict[str, Any] = {}


def _ref_init(reps: FoldReps) -> None:
    _REF_STATE["reps"] = reps


def _ref_one(item: tuple[np.ndarray, str]) -> tuple[str, float, Optional[str]]:
    return fold_label_details(item, _REF_STATE["reps"])


def fold_labels_for(chains: Sequence[Any], reps: FoldReps, n_jobs: Optional[int] = None) -> list[tuple[str, float, Optional[str]]]:
    """`fold_label_details` for many designs, in a spawn pool (scripts need an `if __name__ == '__main__'` guard)."""
    import multiprocessing as mp
    import os

    items = [_trace_of(c) for c in chains]
    jobs = max(1, (os.cpu_count() or 2) - 1) if n_jobs is None else max(1, int(n_jobs))
    jobs = min(jobs, max(1, len(items)))
    if jobs == 1 or len(items) < 8:
        return [fold_label_details(it, reps) for it in items]
    with mp.get_context("spawn").Pool(jobs, initializer=_ref_init, initargs=(reps,)) as pool:
        return pool.map(_ref_one, items, chunksize=4)


# ------------------------------------------------------------------------------------------------ diversity blocks
def _distance_matrix(fams: Sequence[FamilyRecord], fold: FoldLookup, w: float) -> np.ndarray:
    n = len(fams)
    D = np.zeros((n, n))
    fp = [frozenset(f.footprint) for f in fams]
    rows = [fold.row(f.run_id, f.best) if f.best is not None else None for f in fams]
    for i in range(n):
        for j in range(i + 1, n):
            fd = 1 - float(fold.M[rows[i], rows[j]]) / 100 if rows[i] is not None and rows[j] is not None else 1.0
            sd = _jaccard_distance(fp[i], fp[j])
            D[i, j] = D[j, i] = w * fd + (1 - w) * sd
    return D


def _cluster_basis(fams: Sequence[FamilyRecord], D: np.ndarray, theta: float) -> list[dict[str, Any]]:
    out = []
    for ci, m in enumerate(average_linkage(D, theta), 1):
        out.append({"id": ci, "size": len(m), "fams": [fams[i].id for i in m]})
    return out


def _count_block(clusters: list[dict[str, Any]], n: int) -> dict[str, int]:
    sizes = [c["size"] for c in clusters]
    f1, f2 = counts(sizes)
    return {"n": n, "S": len(sizes), "f1": f1, "f2": f2}


def _advice_block(adv: Advice, basis: str, traj_basis: str) -> dict[str, Any]:
    return {
        "coverage": round(adv.coverage, 6), "f0": round(adv.f0, 6), "chao1": round(adv.chao1, 6),
        "hill1": round(adv.hill1, 6), "hill2": round(adv.hill2, 6), "dominance": round(adv.dominance, 6),
        "level": adv.level, "low_confidence": adv.low_conf, "n_traj": adv.n_traj, "traj_basis": traj_basis,
        "yield": {"p": round(adv.yield_p, 6), "lo": round(adv.yield_ci[0], 6), "hi": round(adv.yield_ci[1], 6)},
        "snew": [round(v, 6) for v in adv.snew_list],
        # extras the board text uses
        "p_new_family": round(adv.p_new_family, 6), "recent_new": adv.recent_new, "early_new": adv.early_new,
        "expected_more_600": round(adv.expected_more, 6), "families_more_600": round(adv.fam_more, 6),
        "n_families": adv.n_families, "S": adv.S, "f1": adv.f1, "f2": adv.f2,
    }


def _replay(fams: Sequence[FamilyRecord], clusters: list[dict[str, Any]], n_traj: int) -> list[dict[str, Any]]:
    """Backtest after each family in generation order: what the readout would have said at that trajectory,
    with the clusters of the full clustering restricted to the families found so far (as the board does), and
    what the run went on to find. Coverage is null below 10 families (yield only)."""
    cl_of = {fid: c["id"] for c in clusters for fid in c["fams"]}
    order = sorted((f for f in fams if f.traj > 0), key=lambda f: (f.traj, f.id))
    S_full = len(clusters)
    out = []
    seen: dict[int, int] = {}
    for k, f in enumerate(order, 1):
        seen[cl_of[f.id]] = seen.get(cl_of[f.id], 0) + 1
        sizes = list(seen.values())
        cj = chao_jost(sizes)
        lvl = level_for(k, cj.C)
        yield_p = k / f.traj
        end = len(sizes) + cj.snew(min(k, max(0.0, (n_traj - f.traj) * yield_p)))
        out.append({
            "upto": f.traj, "n_families": k, "S": len(sizes),
            "coverage": None if lvl == "early" else round(cj.C, 6), "level": lvl,
            "found_later": S_full - len(sizes), "chao1": round(cj.chao1, 6), "forecast_at_end": round(end, 6),
        })
    return out


def _sweep(fams: Sequence[FamilyRecord], D: np.ndarray, distances: Sequence[float]) -> list[dict[str, Any]]:
    out = []
    for th in distances:
        sizes = [len(m) for m in average_linkage(D, th)]
        f1, f2 = counts(sizes)
        cj = chao_jost(sizes)
        out.append({"distance": th, "S": len(sizes), "f1": f1, "f2": f2, "coverage": round(cj.C, 6) if sizes else None,
                    "chao1": round(cj.chao1, 6)})
    return out


def _diversity_block(
    run_id: Optional[str], fams: Sequence[FamilyRecord], n_designs: int, fold: FoldLookup, cuts: dict[str, float],
    fold_weight: float, *, estimates_valid: bool, traj: Optional[tuple[str, Optional[int]]], with_by_run: bool = False,
) -> dict[str, Any]:
    weights = {"fold": 1.0, "site": 0.0, "both": fold_weight}
    block: dict[str, Any] = {
        "run": run_id, "n_designs": n_designs, "n_families": len(fams), "estimates_valid": estimates_valid,
        "distance": dict(cuts), "clusters": {}, "counts": {}, "sweep": {},
    }
    advice_d: dict[str, Any] = {}
    replay_by: dict[str, Any] = {}
    for basis in BASES:
        D = _distance_matrix(fams, fold, weights[basis])
        cl = _cluster_basis(fams, D, cuts[basis])
        if with_by_run:
            run_of = {f.id: f.run_id for f in fams}
            for c in cl:
                by: dict[str, int] = {}
                for fid in c["fams"]:
                    by[run_of[fid]] = by.get(run_of[fid], 0) + 1
                c["by_run"] = by
        block["clusters"][basis] = cl
        block["counts"][basis] = _count_block(cl, len(fams))
        block["sweep"][basis] = _sweep(fams, D, SWEEP_DISTANCES)
        if estimates_valid and traj and traj[0] != "none" and traj[1] and fams and all(f.traj > 0 for f in fams):
            cl_of = {fid: c["id"] for c in cl for fid in c["fams"]}
            order = [(f.traj, cl_of[f.id]) for f in sorted(fams, key=lambda f: (f.traj, f.id))]
            adv = advice(order, [c["size"] for c in cl], traj[1])
            advice_d[basis] = _advice_block(adv, basis, traj[0])
            replay_by[basis] = _replay(fams, cl, traj[1])
    if advice_d:
        block["advice"] = advice_d
        block["replay"] = replay_by[BASIS_DEFAULT]
        block["replay_by_basis"] = replay_by
        block["order"] = {"traj_basis": traj[0], "n_traj": traj[1]}
    return block


def diversity(
    run: RunData, fold: dict[str, Any], basis_cuts: Optional[dict[str, float]] = None, *,
    families: Optional[list[FamilyRecord]] = None, fold_weight: float = FOLD_WEIGHT,
    geometry: Optional[dict[tuple[str, str], DesignGeometry]] = None, uniprot_offset: Optional[int] = None,
) -> dict[str, Any]:
    """The diversity block of one run (data_contract.md 3.9): clusters per basis, counts, Chao-Jost advice, replay,
    a distance sweep and the pitfall flags. `fold` is the dict of `fold_matrix` (it must contain the run's
    accepted designs)."""
    cuts = {**DEFAULT_DISTANCE, **(basis_cuts or {})}
    fams = families if families is not None else build_families(run, geometry, uniprot_offset=uniprot_offset)
    block = _diversity_block(
        run.run_id, fams, len(run.accepted()), decode_fold(fold), cuts, fold_weight, estimates_valid=True, traj=traj_info(run),
    )
    block["epoch"] = {"hash": epoch_hash(run, uniprot_offset), "components": epoch_components(run, uniprot_offset)}
    block["epitope"] = epitope_diversity(run, fams, geometry, uniprot_offset=uniprot_offset)
    # One scaffold makes fold distance measure the scaffold; a short or cyclic
    # binder puts it outside what TM-score can say at all. Either way the figure
    # would be read as diversity, so it is not reported.
    why_off = FOLD_OFF_WHY if _para.has_scaffold(run) else fold_unreliable(run)
    if why_off:
        _switch_fold_off(block, why_off)
    block["flags"] = pitfall_flags(run, fams, block)
    if why_off:
        block["flags"].insert(0, _fold_off_flag(why_off))
    return block


def diversity_all(
    runs: Sequence[RunData], fold: dict[str, Any], basis_cuts: Optional[dict[str, float]] = None, *,
    families: Optional[dict[str, list[FamilyRecord]]] = None, fold_weight: float = FOLD_WEIGHT,
    geometry: Optional[dict[tuple[str, str], DesignGeometry]] = None, uniprot_offset: Optional[int] = None,
) -> dict[str, Any]:
    """campaign.diversity: per-run blocks plus a pooled block. The pooled block clusters the families of all runs
    together (site overlap compares UniProt numbers, fold compares TM-align across runs) and is descriptive: it
    has `estimates_valid` false, and no advice or replay, when the runs' epochs differ."""
    cuts = {**DEFAULT_DISTANCE, **(basis_cuts or {})}
    fams_by = {r.run_id: (families or {}).get(r.run_id) or build_families(r, geometry, uniprot_offset=uniprot_offset) for r in runs}
    runs_out = {r.run_id: diversity(r, fold, cuts, families=fams_by[r.run_id], fold_weight=fold_weight, geometry=geometry,
                                    uniprot_offset=uniprot_offset) for r in runs}
    hashes = {r.run_id: epoch_hash(r, uniprot_offset) for r in runs}
    valid = len(set(hashes.values())) <= 1
    all_fams = [f for r in runs for f in fams_by[r.run_id]]
    pooled = _diversity_block(
        None, all_fams, sum(len(r.accepted()) for r in runs), decode_fold(fold), cuts, fold_weight,
        estimates_valid=valid, traj=traj_info(runs[0]) if valid and len(runs) == 1 else None, with_by_run=True,
    )
    pooled["epoch"] = {"hashes": hashes, "differences": epoch_differences(runs, uniprot_offset)}
    off_why = {r.run_id: (FOLD_OFF_WHY if _para.has_scaffold(r) else fold_unreliable(r)) for r in runs}
    scaffolded_runs = [rid for rid, why in off_why.items() if why]
    if scaffolded_runs and len(scaffolded_runs) == len(runs):
        # Every run has a reason; the pooled figure would be the same mistake.
        reasons = sorted(set(off_why.values()))
        _switch_fold_off(pooled, reasons[0] if len(reasons) == 1 else "; and ".join(reasons))
        pooled["flags"] = [_fold_off_flag(reasons[0] if len(reasons) == 1 else "; and ".join(reasons))]
    elif scaffolded_runs:
        pooled["flags"] = [{
            "code": "fold_scaffolded", "severity": "info",
            "text": f"{', '.join(scaffolded_runs)} {'is' if len(scaffolded_runs) == 1 else 'are'} built on a fixed scaffold, so within "
                    f"{'it' if len(scaffolded_runs) == 1 else 'them'} fold distance measures the scaffold; the pooled fold and combined counts mix that "
                    "with the other runs and are descriptive only.", "data": {"runs": scaffolded_runs}}]
    if not valid:
        pooled["flags"] = list(pooled.get("flags", [])) + [{
            "code": "epoch", "severity": "warn",
            "text": "The runs used different settings (" + ", ".join(pooled["epoch"]["differences"]) + "). The pooled view is "
                    "descriptive; coverage and forecasts are not given across epochs.",
            "data": {"hashes": hashes, "differences": pooled["epoch"]["differences"]},
        }]
    return {
        "schema": SCHEMA, "basis_default": BASIS_DEFAULT, "distance_default": dict(cuts), "fold_weight": fold_weight,
        "runs": runs_out, "pooled": pooled,
    }


# ------------------------------------------------------------------------------------------------ scaffolded runs and epitope spread

#: Below this many residues, fold distance is not reported. TM-align normalises
#: by d0 = 1.24*(L-15)^(1/3) - 1.8: at L <= 15 that is negative and the
#: implementation clamps d0 to 0.5, and it only rises above the clamp past
#: L ~ 21, so a score on a short chain is not length-calibrated. The "TM > 0.5
#: means the same fold" rule it would be read with was established on protein
#: domains, not on peptides: two unrelated 15-mers, both extended or both
#: helical, clear 0.5 routinely. This code also normalises by the SHORTER chain,
#: the most generous of the choices, so a 12-mer against a 24-mer is scored with
#: d0 taken from 12.
FOLD_MIN_LENGTH = 30
FOLD_SHORT_WHY = (
    "the binders are too short for fold distance to mean anything: TM-align's normalisation is clamped below about "
    "21 residues and undefined at 15, and the 0.5 cut that would be read as \"the same fold\" was calibrated on protein "
    "domains, so short peptides of unrelated sequence score above it as a matter of course"
)
FOLD_MULTICHAIN_WHY = (
    "the binder is written as several chains, and fold distance compares one CA trace: run end to end, the trace would "
    "jump from the last residue of one chain to the first of the next and align a connection the molecule does not have"
)
FOLD_CYCLIC_WHY = (
    "a cyclic binder has no chain ends, and TM-align compares it as an open CA trace, so the alignment it reports is "
    "not of the molecule that was designed"
)


def fold_unreliable(run: RunData) -> Optional[str]:
    """Why fold distance should not be reported for this run, or None when it can be.

    Two reasons, both about the measure rather than the designs. A campaign of
    short binders (peptides, and most cyclic peptides) is scored outside the
    range TM-score is calibrated on; a cyclic campaign is scored as a molecule
    it is not. In both cases the site axis -- where on the target a design binds,
    from its contact footprint -- is unaffected and carries the diversity, which
    is what `_switch_fold_off` leaves in place.
    """
    modality = str((run.annotations or {}).get("modality") or "").lower()
    if "cyclic" in modality:
        return FOLD_CYCLIC_WHY
    chains = max((binder_chain_count(d) for d in run.accepted()), default=1)
    if chains > 1:
        return FOLD_MULTICHAIN_WHY
    lengths = sorted(d.length for d in run.accepted() if d.length)
    if lengths and lengths[len(lengths) // 2] < FOLD_MIN_LENGTH:
        shortest, longest = lengths[0], lengths[-1]
        span = f"{shortest}" if shortest == longest else f"{shortest} to {longest}"
        return f"{FOLD_SHORT_WHY} (the accepted binders here are {span} residues)"
    return None


def _switch_fold_off(block: dict[str, Any], why: str = FOLD_OFF_WHY) -> None:
    """Switch the fold-based figures off in a diversity block (in place): every basis that includes fold distance is removed, and the
    block says so and names the basis that is left. A scaffolded run's designs share one fold by construction, so clustering by fold (or
    by fold and site together) reports a saturated coverage that means "these share a scaffold", not "the space is explored"."""
    for key in ("clusters", "counts", "sweep", "advice", "replay_by_basis"):
        sub = block.get(key)
        if isinstance(sub, dict):
            for b in FOLD_BASES:
                sub.pop(b, None)
            if not sub:
                del block[key]
    rb = block.get("replay_by_basis")
    if isinstance(rb, dict) and SITE_BASIS in rb:
        block["replay"] = rb[SITE_BASIS]
    else:
        block.pop("replay", None)
    block["fold"] = {"status": "off", "why": why}
    block["basis_off"] = list(FOLD_BASES)
    block["basis_default"] = SITE_BASIS


def _fold_off_flag(why: str = FOLD_OFF_WHY) -> dict[str, Any]:
    return {"code": "fold_off", "severity": "info",
            "text": f"Fold diversity is switched off for this run: {why}. Clusters are counted on the binding site alone; see the paratope and "
                    "epitope figures for what varies.", "data": {"bases_off": list(FOLD_BASES)}}


_RESIDUE_TOKEN = re.compile(r"^[A-Za-z]?(\d+)$")


def _target_residues(d: DesignRecord) -> Optional[set[int]]:
    """The target residues BindCraft 2 lists for a design (``Interface_Target_Residues``, construct numbers), or None when the column is
    absent, empty or in a layout this does not read (several target chains)."""
    cell = d.extras.get("Interface_Target_Residues") if d.extras else None
    if not isinstance(cell, str):
        return None
    text = cell.strip().strip("[]()")
    if not text or "/" in text:
        return None
    out: set[int] = set()
    for tok in re.split(r"[,;\s|]+", text):
        tok = tok.strip("'\"")
        if not tok:
            continue
        m = _RESIDUE_TOKEN.match(tok)
        if not m:
            return None
        out.add(int(m.group(1)))
    return out or None


def epitope_diversity(
    run: RunData, fams: Sequence[FamilyRecord], geometry: Optional[dict[tuple[str, str], DesignGeometry]] = None, *,
    uniprot_offset: Optional[int] = None,
) -> dict[str, Any]:
    """How far apart the families' epitopes are, for any run (``diversity.runs[run].epitope``).

    The epitope of a family is the set of target residues contacted by any of its accepted designs. The block gives the range of the
    pairwise Jaccard overlap between families and the residues every family contacts. It does not depend on the fold or on the modality,
    which is why it is the axis that still means something when fold diversity is switched off.

    The residues come from BindCraft 2's own ``Interface_Target_Residues`` when every accepted design carries it, else from the contacts
    recomputed from the structures (4.5 angstrom, as everywhere on the pages); the block says which (``source``) and never mixes the two.
    Status ``unknown`` (with a reason) when fewer than two families have an epitope: nothing is invented.
    """
    num, numbering = _numberer(run, uniprot_offset)
    members = [(f, [run.design(i) for i in f.accepted]) for f in fams if f.accepted]
    base: dict[str, Any] = {"n_families": len(members)}
    if len(members) < 2:
        return {"status": "unknown", "why": "fewer than two families, so there is no overlap between epitopes to measure", **base}
    cols = [[_target_residues(d) for d in ds] for _f, ds in members]
    if all(c is not None for cs in cols for c in cs):
        source = "BindCraft Interface_Target_Residues"
        sets = [set().union(*cs) for cs in cols]
    else:
        geo_sets: list[set[int]] = []
        missing = False
        for _f, ds in members:
            u: set[int] = set()
            for d in ds:
                g = (geometry or {}).get(d.key)
                if g is None and d.structure_path:
                    g = design_geometry(d)
                if g is None:
                    missing = True
                    break
                u.update(g.contacts)
            if missing:
                break
            geo_sets.append(u)
        if missing:
            return {"status": "unknown", "why": "the accepted designs list no target residues and have no structure to take contacts from", **base}
        source = "contacts recomputed from the structures (4.5 A)"
        sets = geo_sets
    kept = [(f, {num(r) for r in st}) for (f, _d), st in zip(members, sets) if st]
    if len(kept) < 2:
        return {"status": "unknown", "why": "fewer than two families contact any target residue", **base}
    js = []
    for i in range(len(kept)):
        for j in range(i + 1, len(kept)):
            a, b = kept[i][1], kept[j][1]
            js.append(len(a & b) / len(a | b))
    shared = sorted(set.intersection(*(st for _f, st in kept)))
    n = len(kept)
    return {
        "status": "observed", "n_families": n, "pairwise_jaccard": [round(min(js), 2), round(max(js), 2)], "shared_by_all": shared,
        "sizes": {f.id: len(st) for f, st in kept}, "source": source, "numbering": numbering,
        "rule": "the epitope of a family is every target residue contacted by any of its accepted designs",
        "low_confidence": n < LOW_CONFIDENCE_BELOW,
        "note": f"{n} families" + (f": too few to call a trend (estimates are low-confidence below {LOW_CONFIDENCE_BELOW} families)." if n < LOW_CONFIDENCE_BELOW else "."),
    }


def loop_length_classes(lengths: Sequence[int], gap: int = LOOP_CLASS_GAP) -> list[list[int]]:
    """Lengths grouped into classes: sorted, a new class starts where the next length is more than ``gap`` residues past the last one."""
    out: list[list[int]] = []
    for v in sorted(lengths):
        if out and v - out[-1][-1] <= gap:
            out[-1].append(v)
        else:
            out.append([v])
    return out


def diversity_paratope(para: "_para.RunParatope", families: dict[str, Sequence[str]], epitope: Optional[dict[str, Any]]) -> dict[str, Any]:
    """``runs[].diversity_paratope``: the axes that replace fold diversity for a scaffolded run. Each axis is ``observed`` or ``unknown``
    (or ``off`` for fold, with its reason); none is invented. ``families`` maps a family id to its accepted design ids."""
    fam_of = {d: f for f, ids in families.items() for d in ids}
    acc = [para.designs[i] for i in para.accepted if i in para.designs]
    fam_ids = sorted({fam_of[d.id] for d in acc if d.id in fam_of})
    n_fam = len(fam_ids)
    low = n_fam < LOW_CONFIDENCE_BELOW
    out: dict[str, Any] = {"fold": {"status": "off", "why": FOLD_OFF_WHY}}
    if acc:
        profiles = sorted({tuple(d.loops) for d in acc}, key=lambda t: (-len(t), t))
        out["paratope"] = {"status": "observed", "profiles": len(profiles), "n": len(acc), "n_families": n_fam,
                           "profile_list": [list(p) for p in profiles],
                           "note": "distinct region-usage profiles among accepted designs: which designed segments touch the target"}
    else:
        out["paratope"] = {"status": "unknown", "why": "no accepted design could be placed on the scaffold"}
    # loop length: the segment whose length the spec lets vary most
    varying = [s for s in para.layout.segments if s.variable and s.length[1] > s.length[0]]
    if not varying:
        out["loop_length"] = {"status": "off", "why": "the campaign's edit spec fixes the length of every designed segment, so there is no loop length to compare"}
    else:
        seg = max(varying, key=lambda s: (s.length[1] - s.length[0], -s.index))
        by_fam = para.family_lengths(seg.name, fam_of)
        mixed = sorted(f for f, v in by_fam.items() if len(v) > 1)
        if len(by_fam) < 2:
            out["loop_length"] = {"status": "unknown", "segment": seg.name, "why": "fewer than two families with a placed design, so there are no lengths to compare"}
        elif mixed:
            out["loop_length"] = {"status": "unknown", "segment": seg.name,
                                  "why": f"the accepted designs of {len(mixed)} famil{'y' if len(mixed) == 1 else 'ies'} differ in {seg.name} length, which one backbone cannot do"}
        else:
            values = sorted(v[0] for v in by_fam.values())
            cls = loop_length_classes(values)
            out["loop_length"] = {
                "status": "observed", "segment": seg.name, "values": values, "distinct": len(set(values)), "classes": len(cls),
                "class_ranges": [[c[0], c[-1]] for c in cls], "class_families": [len(c) for c in cls], "class_gap": LOOP_CLASS_GAP,
                "note": f"length of the {seg.name} stretch the campaign let vary (scaffold {', '.join(f'{a}-{b}' for a, b, _l, _h in seg.variable)}), one value "
                        f"per family; lengths within {LOOP_CLASS_GAP} residues of the next belong to one class"}
    out["epitope"] = epitope if epitope is not None else {"status": "unknown", "why": "not computed"}
    out["n_families"] = n_fam
    out["low_confidence"] = low
    if low:
        out["note"] = f"{n_fam} families: too few to call a trend (estimates are low-confidence below {LOW_CONFIDENCE_BELOW} families)."
    return out


# ------------------------------------------------------------------------------------------------ pitfalls
def pitfall_flags(run: RunData, fams: Sequence[FamilyRecord], block: dict[str, Any]) -> list[dict[str, Any]]:
    """Things that make the estimates weaker than they look. Each flag: code, severity (info or warn), text, data."""
    flags: list[dict[str, Any]] = []
    n = len(fams)
    basis, n_traj = traj_info(run)
    if run.funnel.trajectories_attempted is None:
        est = run.funnel.trajectories_attempted_estimate
        flags.append({
            "code": "unknown_attempts", "severity": "info",
            "text": "The number of attempted trajectories is unknown (no Trajectory folders), so yield is per trajectory that "
                    "completed and passed, not per attempt." + (f" A derived estimate of attempts is {est}." if est else ""),
            "data": {"traj_basis": basis, "n_traj": n_traj, "attempted_estimate": est},
        })
    if n < MIN_FAMILIES:
        flags.append({"code": "small_sample", "severity": "warn", "text": f"{n} families: too few for coverage estimates (yield only below {MIN_FAMILIES}).", "data": {"n_families": n}})
    elif n < LOW_CONFIDENCE_BELOW:
        flags.append({"code": "small_sample", "severity": "info", "text": f"{n} families: estimates are shown with low confidence below {LOW_CONFIDENCE_BELOW}.", "data": {"n_families": n}})
    cap, cap_src = _sequence_cap(run)
    at_cap = sum(1 for f in fams if f.at_sequence_cap)
    if at_cap:
        flags.append({
            "code": "sibling_cap", "severity": "info",
            "text": f"BindCraft keeps at most {cap} accepted sequences per backbone. {at_cap} of {n} families are at that limit, so "
                    "their sibling counts are censored (counts of sequences per family are lower bounds); counting families is not affected.",
            "data": {"cap": cap, "cap_source": cap_src, "families_at_cap": at_cap, "n_families": n},
        })
    fc = run.failure_counters or {}
    af2 = {k: fc[k] for k in ("pLDDT", "i_pTM", "i_pAE") if fc.get(k)}
    if af2:
        flags.append({
            "code": "unlogged_af2_failures", "severity": "info",
            "text": "Sequences that fail the AlphaFold2 pre-filters leave no row (only a counter, and a sequence counts once per filter it fails), "
                    "so the sequence-level funnel cannot be reconstructed from the rows.",
            "data": {"counters": af2},
        })
    fold_sweep = {r["distance"]: r for r in block.get("sweep", {}).get("fold", [])}
    tight = fold_sweep.get(0.2)
    if tight and n >= 3 and tight["S"] >= 0.9 * n:
        flags.append({
            "code": "threshold_artefact", "severity": "warn",
            "text": f"At TM 0.8 (distance 0.2) {tight['S']} of {n} families are separate fold clusters, so the fold count cannot saturate at that "
                    "setting. Use TM 0.5 for 'same fold'; 0.8 is the 'same backbone' cut.",
            "data": {"S": tight["S"], "n_families": n, "distance": 0.2},
        })
    flags.append({
        "code": "restart_counters", "severity": "info",
        "text": "A restarted BindCraft run resets its own counters (failure_csv.csv) and cannot be detected from the output files; the row "
                "order of trajectory_stats.csv is the only generation order available.",
        "data": {},
    })
    flags.append({
        "code": "extrapolation_limit", "severity": "info",
        "text": "Forecasts of new clusters are valid for at most as many further families as already seen (m <= n).",
        "data": {"n_families": n},
    })
    flags.append({
        "code": "predicted_only", "severity": "info",
        "text": "A passing design is a predicted binder, not a measured one.", "data": {},
    })
    return flags


def chronology(run: RunData) -> dict[str, Any]:
    """What can and cannot be said about the order of trajectories in this run."""
    n = len(run.backbones)
    times = [b.trajectory_time_s for b in run.backbones if b.trajectory_time_s is not None]
    return {
        "order_source": "row order of trajectory_stats.csv (BindCraft appends a row when a trajectory finishes and passes)",
        "n_rows": n,
        "time_column": "TrajectoryTime",
        "time_is": "duration of each trajectory, not a timestamp; the order cannot be verified from it",
        "total_trajectory_hours": round(sum(times) / 3600.0, 2) if times else None,
        "verified": False,
    }


# ------------------------------------------------------------------------------------------------ orchestration
@dataclass
class FamilyAnalysis:
    runs: list[str]
    families: dict[str, list[FamilyRecord]]
    fold: dict[str, Any]
    fold_labels: dict[str, str]
    diversity: dict[str, Any]
    design_ids: dict[tuple[str, str], str]
    chronology: dict[str, dict[str, Any]]
    rejected_fold: dict[str, dict[str, Any]] = field(default_factory=dict)   # project id -> {label, tm, nearest}
    warnings: list[str] = field(default_factory=list)
    seconds: float = 0.0

    def run_result(self, run_id: str) -> dict[str, Any]:
        """JSON-serialisable result for one run (work/<project>/families/<run>.json)."""
        return {
            "schema": SCHEMA, "run": run_id,
            "families": [f.to_dict() for f in self.families[run_id]],
            "diversity": self.diversity["runs"][run_id],
            "chronology": self.chronology[run_id],
            "fold_labels": {pid: self.fold_labels[pid] for (r, _), pid in sorted(self.design_ids.items())
                            if r == run_id and pid in self.fold_labels},
            "warnings": [w for w in self.warnings if w.startswith(run_id)],
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA, "runs": {r: self.run_result(r) for r in self.runs}, "fold": self.fold,
            "fold_labels": self.fold_labels, "diversity": self.diversity,
            "rejected_fold": self.rejected_fold, "warnings": self.warnings,
        }


def analyze_families(
    runs: Sequence[RunData], *, uniprot_offset: Optional[int] = None, n_jobs: Optional[int] = None,
    basis_cuts: Optional[dict[str, float]] = None, fold_weight: float = FOLD_WEIGHT, rejected_labels: bool = False,
) -> FamilyAnalysis:
    """Everything this module produces for a project: structures are read once, families per run, the fold matrix
    over the accepted designs of all runs, pooled fold labels, the diversity block and (optionally) the fold label of
    every rejected design with a structure."""
    t0 = time.time()
    warnings: list[str] = []
    accepted = [d for r in runs for d in r.accepted() if d.structure_path]
    geo = load_geometry(accepted)
    pids = project_design_ids(runs)
    fams = {r.run_id: build_families(r, geo, uniprot_offset=uniprot_offset, warnings=warnings) for r in runs}
    fold = fold_matrix(accepted, geo, id_of=lambda d: pids[d.key], n_jobs=n_jobs)
    div = diversity_all(runs, fold, basis_cuts, families=fams, fold_weight=fold_weight, geometry=geo, uniprot_offset=uniprot_offset)
    labels = fold_labels(fold, {**DEFAULT_DISTANCE, **(basis_cuts or {})}["fold"])
    rejected: dict[str, dict[str, Any]] = {}
    if rejected_labels:
        reps = build_fold_reps(accepted, fold, geo, {**DEFAULT_DISTANCE, **(basis_cuts or {})}["fold"])
        rej = [d for r in runs for d in r.rejected_with_structure()]
        if rej:
            chains = [design_geometry(d).chain for d in rej]
            for d, (lab, tm, near) in zip(rej, fold_labels_for(chains, reps, n_jobs)):
                rejected[pids[d.key]] = {"label": lab, "tm": tm, "nearest": near, "run": d.run_id}
    return FamilyAnalysis(
        runs=[r.run_id for r in runs], families=fams, fold=fold, fold_labels=labels, diversity=div, design_ids=pids,
        chronology={r.run_id: chronology(r) for r in runs}, rejected_fold=rejected, warnings=warnings,
        seconds=round(time.time() - t0, 2),
    )


def write_artefacts(analysis: FamilyAnalysis, out_dir: str | Path) -> list[Path]:
    """Write <run>.json for every run and fold.json (the fold matrix, pooled labels, rejected labels) into out_dir."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for rid in analysis.runs:
        p = out / f"{rid}.json"
        p.write_text(json.dumps(analysis.run_result(rid), indent=1, sort_keys=True, default=_json_default) + "\n")
        paths.append(p)
    p = out / "fold.json"
    p.write_text(json.dumps({
        "schema": SCHEMA, "fold": analysis.fold, "labels": analysis.fold_labels, "rejected": analysis.rejected_fold,
        "diversity_pooled": analysis.diversity["pooled"], "distance_default": analysis.diversity["distance_default"],
    }, indent=1, sort_keys=True, default=_json_default) + "\n")
    paths.append(p)
    return paths


def _json_default(o: Any) -> Any:
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """python -m binder_review.families --run run01=PATH[,SETTINGS.json] --run run02=PATH --construct FILE --offset N --out DIR"""
    import argparse

    from .readers import bindcraft1

    ap = argparse.ArgumentParser(prog="binder_review.families", description="Families, fold clusters and diversity per run.")
    ap.add_argument("--run", action="append", required=True, help="id=run folder[,settings json]")
    ap.add_argument("--construct", help="construct PDB handed to BindCraft")
    ap.add_argument("--offset", type=int, help="uniprot = construct + offset when the runs carry none")
    ap.add_argument("--out", required=True)
    ap.add_argument("--rejected", action="store_true", help="also label rejected designs with a structure")
    a = ap.parse_args(argv)
    runs = []
    for spec in a.run:
        rid, rest = spec.split("=", 1)
        path, _, settings = rest.partition(",")
        runs.append(bindcraft1.read_run(path, rid, settings_path=settings or None, construct_path=a.construct))
    res = analyze_families(runs, uniprot_offset=a.offset, rejected_labels=a.rejected)
    for p in write_artefacts(res, a.out):
        print(p)
    print(f"{len(runs)} runs in {res.seconds} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
