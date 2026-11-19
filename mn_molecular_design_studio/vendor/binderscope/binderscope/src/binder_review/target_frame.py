"""Target frame: judge a binder designed against a truncated construct as if it sat on the full-length protein.

BindCraft designs a binder against a fragment of the target (a construct: a stretch of the chain, possibly with loops
deleted). The frame places every design on the full-length reference chain by a number-matched CA superposition,
weighted around the binder (so the epitope is where it should be even when the rest of the fragment moved) or over the whole
target, and then asks what the fragment hid: the other subunits of an oligomer, the rest of the chain, the glycans, the cut
ends and the membrane. See docs/target_frame.md for the API, the definitions and the verification.

Conventions (docs/CONTRACTS.md)
  * Three numberings: construct (what BindCraft wrote, 1 to N, with gaps where loops were deleted), the author numbering
    of the reference structure and, when an accession is given, UniProt. The reference numbering must equal UniProt
    numbering then (checked by sequence); without an accession the reference numbering is the common numbering. The
    construct offset is found by sequence, never assumed. Everything this module returns for the blobs is numbered in
    that common numbering (called "UniProt" in the field names); construct numbers appear only where noted.
  * Every coordinate returned for a blob is in the reference frame (the reference chain). A transform is (R, t) with
    ``x' = x @ R.T + t``, the same as ``x' = R x + t`` for a column vector.
  * A frame is data: ``TargetFrame.meta`` (JSON-able) and ``TargetFrame.arrays`` (numpy). ``save`` and ``load`` write
    and read exactly those two, so later stages never recompute it.
  * Metrics only, no verdicts, in what goes to the blobs. ``verdicts`` exists for QA and reads the same thresholds as
    profile.checks.context (``DEFAULT_CONFIG["checks"]``, kept equal to the built-in default profile by a test).
  * UniProt, SIFTS and AlphaFold records are optional: without an accession the frame is derived from the reference
    structure and the construct alone (no domain strip, no named topology, no membrane plane) and says so in ``notes``.
"""
from __future__ import annotations

import copy
import datetime as _dt
import hashlib
import json
import math
from collections import Counter, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.spatial import cKDTree

from . import model, structure as S

__all__ = [
    "FRAME_SCHEMA", "DEFAULT_CONFIG", "FrameError", "RunSpec", "Glycan", "Placement", "ContextResult", "DesignCheck",
    "TargetFrame", "build_frame", "classify_reference", "auto_partner_chains", "reference_label", "ReferenceChains", "analyse_reference",
    "describe_reference", "format_reference_chains", "place_design", "check_design", "check_run", "verdicts", "membrane_excess",
    "plane_from_points", "trimmed_fit", "interface_weights", "interface_fit", "PLACEMENT_FITS", "cut_points_for", "ss_records", "report_tables",
    "METRIC_KEYS",
]

FRAME_SCHEMA = "binderscope.target_frame/1"
FRAME_JSON = "frame.json"
FRAME_NPZ = "frame.npz"

#: Thresholds and geometry settings. "checks" mirrors the default profile (profile.checks.context);
#: tests/test_target_frame.py asserts the two are equal. "geometry" holds what the profile does not carry.
DEFAULT_CONFIG: dict[str, Any] = {
    "checks": {
        "partner": {"fail_n3": 10, "warn_min": 4.0, "warn_buried": 3},
        "rest": {"fail_n3": 10, "warn_min": 4.0},
        "cut": {"k": 4, "warn_n": 2},
        "cover": {"min_area": 10.0, "warn_n": 2},
        "near": {"bad_lt": 3.0, "warn_lt": 4.0, "info_lt": 8.0},
    },
    "geometry": {
        "contact_cutoff": model.CONTACT_CUTOFF,   # a contact on the pages (tc, hot)
        "clash": 2.2,                             # binder heavy atom this close to partner or rest atoms is a steric clash
        "n_close": [3.0, 4.0, 6.0],               # atom-count radii (olig_n3, olig_n4, ctx.n6)
        "epitope_radius": 5.0,                    # epitope = hotspots + residues within this of a hotspot residue
        "olig_min_area": 5.0,                     # target.olig.min_area: buried at the partner above this (A^2)
        "membrane_offset": 6.0,                   # the plane lies this far short of the chain end (A)
        "place_warn_rmsd": 1.5,                   # exporter warning threshold for the CA RMSD near the interface (A)
        "fit": "interface",                       # how a design is placed: "interface" (weighted around the binder), "all" (every shared CA) or "core" (trimmed)
        "iface_full": 8.0,                        # interface fit: a CA this close to a binder heavy atom has weight 1 (A)
        "iface_decay": 6.0,                       # beyond that the weight falls as exp(-(d / this)^2): 1/e this far beyond iface_full (A)
        "iface_floor": 0.01,                      # the smallest weight of any shared CA
        "iface_cutoff": 2.5,                      # a residue near the binder that deviates more than this is left out of the refit (A)
        "iface_min_n": 12,                        # fewer residues near the binder than this: the whole-target fit is used instead
        "core_cutoff": 2.0,                       # trimmed refit drops residues deviating more than this (A)
        "core_cycles": 5,
        "sasa_probe": 1.4,
        "sasa_points": 960,
    },
}

AA3 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS",
       "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
       "Y": "TYR", "V": "VAL"}

#: Colour classes of the frame strip (profile.frame.segments[].kind): x outside, v first domain, c second, l linker. Domains
#: are named by their CATH homology name (a profile may name its own segments).
DOMAIN_KINDS = ("v", "c", "b", "s", "h")
#: Chain ids given to the sugars of in-construct glycan sites in target_pdb and glycan_atoms (A and B are the polymers).
GLYCAN_CHAINS = "CDEFGHIJKLMNOPQRSTUVWXYZ"


class FrameError(ValueError):
    """The inputs do not fit together (numbering, chains, runs) or a frame file is unusable."""


# ----------------------------------------------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------------------------------------------
def _ints(xs: Iterable) -> list[int]:
    return [int(x) for x in xs]


def _ranges(numbers: Iterable[int]) -> list[list[int]]:
    """Contiguous [first, last] stretches of a set of residue numbers, ascending."""
    out: list[list[int]] = []
    for n in sorted(set(int(x) for x in numbers)):
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return out


def _expand(ranges: Iterable[Sequence[int]]) -> set[int]:
    s: set[int] = set()
    for a, b in ranges:
        s.update(range(int(a), int(b) + 1))
    return s


def _complement(ranges: Sequence[Sequence[int]], lo: int, hi: int) -> list[list[int]]:
    """Stretches of lo..hi that the given ranges do not cover."""
    have = _expand(ranges)
    return _ranges(n for n in range(lo, hi + 1) if n not in have)


def cut_points_for(resi_ranges: Sequence[Sequence[int]], native_range: Sequence[int]) -> list[int]:
    """Residues next to a cut: the first and last residue unless the construct reaches the native end (first within one
    residue of the native start, or last at the native end) and both flanks of every removed stretch."""
    rr = [list(map(int, r)) for r in resi_ranges]
    if not rr:
        return []
    pts: set[int] = set()
    first, last = rr[0][0], rr[-1][1]
    if first > native_range[0] + 1:
        pts.add(first)
    if last < native_range[1]:
        pts.add(last)
    for a, b in zip(rr, rr[1:]):
        pts.update((a[1], b[0]))
    return sorted(pts)


def _round(x: Any, nd: int) -> Any:
    if x is None:
        return None
    if isinstance(x, np.ndarray):
        return np.round(x.astype(float), nd).tolist()
    if isinstance(x, (list, tuple)):
        return [_round(v, nd) for v in x]
    return round(float(x), nd)


def _flat(xyz: np.ndarray, nd: int = 2) -> list[float]:
    return np.round(np.asarray(xyz, dtype=float).reshape(-1), nd).tolist()


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n < 1e-9:
        raise FrameError("zero-length vector where a direction is needed")
    return np.asarray(v, dtype=float) / n


def _sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _json_default(o: Any):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    raise TypeError(f"not JSON serialisable: {type(o)}")


def ss_records(ss: str, numbers: Sequence[int], names: Sequence[str], chain: str = "A") -> list[str]:
    """HELIX records, then SHEET records, from an H/E/L string (one letter per entry of ``numbers``). A helix needs 4
    residues, a strand 3, and a stretch must be contiguous in residue number. The record layout is the one the mockup
    pages read."""
    helices: list[str] = []
    strands: list[str] = []
    n = len(ss)
    i = 0
    while i < n:
        if ss[i] in "HE":
            j = i
            while j + 1 < n and ss[j + 1] == ss[i] and numbers[j + 1] == numbers[j] + 1:
                j += 1
            if j - i + 1 >= (4 if ss[i] == "H" else 3):
                a, b = numbers[i], numbers[j]
                if ss[i] == "H":
                    k = len(helices) + 1
                    helices.append("HELIX  %3d %3d %3s %1s %4d  %3s %1s %4d  1                               %5d"
                                   % (k, k, names[i], chain, a, names[j], chain, b, j - i + 1))
                else:
                    k = len(strands) + 1
                    strands.append("SHEET  %3d %3s 1 %3s %1s%4d  %3s %1s%4d  0"
                                   % (k, "S%02d" % k, names[i], chain, a, names[j], chain, b))
            i = j + 1
        else:
            i += 1
    return helices + strands


def default_view_axis(ca: np.ndarray) -> np.ndarray:
    """A unit vector to orient a view by when no site gives one: the axis along which the points (CA atoms) vary least, so that a view
    along it shows the flat face; the sign is fixed by making its largest component positive. [0, 0, 1] for fewer than 3 points or
    points on a line."""
    pts = np.asarray(ca, dtype=float).reshape(-1, 3)
    if len(pts) >= 3:
        _, sv, vt = np.linalg.svd(pts - pts.mean(axis=0), full_matrices=False)
        if sv[0] > 1e-6:
            axis = vt[-1] / np.linalg.norm(vt[-1])
            return axis if axis[int(np.argmax(np.abs(axis)))] >= 0 else -axis
    return np.array([0.0, 0.0, 1.0])


def plane_from_points(centroid: np.ndarray, end: np.ndarray, short_by: float) -> dict[str, Any]:
    """Exclusion plane at ``short_by`` angstrom before ``end``, seen from ``centroid``: axis is the unit vector from the
    centroid to the end, offset the distance between them minus ``short_by``. A point x is past the plane when
    (x - centroid) . axis > offset."""
    d = np.asarray(end, dtype=float) - np.asarray(centroid, dtype=float)
    axis = _unit(d)
    return {"centroid": np.asarray(centroid, dtype=float), "axis": axis, "offset": float(d @ axis) - float(short_by)}


def membrane_excess(xyz: np.ndarray, plane: Mapping[str, Any]) -> float:
    """Largest signed distance of the points beyond the plane (negative: all points are clear by that much)."""
    xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
    if not len(xyz):
        return float("nan")
    centroid, axis = np.asarray(plane["centroid"], dtype=float), np.asarray(plane["axis"], dtype=float)
    return float(((xyz - centroid) @ axis - float(plane["offset"])).max())


def trimmed_fit(P: np.ndarray, Q: np.ndarray, cutoff: float = 2.0, cycles: int = 5, min_frac: float = 0.5
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Kabsch fit with outlier rejection: refit on the pairs that deviate by at most ``cutoff`` until the set stops
    changing (at most ``cycles`` rounds). Returns (R, t, keep). Falls back to all pairs when fewer than ``min_frac``
    would stay."""
    keep = np.ones(len(P), dtype=bool)
    R, t, _ = S.kabsch(P, Q)
    for _ in range(max(1, cycles)):
        dev = np.linalg.norm(P @ R.T + t - Q, axis=1)
        new = dev <= cutoff
        if new.sum() < max(3, math.ceil(min_frac * len(P))) or (new == keep).all():
            break
        keep = new
        R, t, _ = S.kabsch(P[keep], Q[keep])
    return R, t, keep


#: The ways a design can be placed on the reference chain (DEFAULT_CONFIG["geometry"]["fit"]).
PLACEMENT_FITS = ("interface", "all", "core")
#: A shared CA with at least this weight is "near" the binder: only such residues are left out of the interface refit as outliers, and
#: at least ``iface_min_n`` of them are needed for an interface fit.
IFACE_NEAR_WEIGHT = 0.3
#: ``Placement.iface_rmsd`` and ``iface_n`` cover the shared CA with at least this weight (within about 13 A of a binder atom).
IFACE_REPORT_WEIGHT = 0.5
#: Rounds of the interface refit (as many as the trimmed fit).
IFACE_CYCLES = 5
#: A residue that is not near the binder is left out of the interface refit when it deviates by more than this many times ``iface_cutoff``
#: (a part of the target that hinges against the part the binder sits on: far from the binder it would act as a long lever on the rotation).
IFACE_FAR_FACTOR = 3.0


def interface_weights(dist: np.ndarray, full: float = 8.0, decay: float = 6.0, floor: float = 0.01) -> np.ndarray:
    """Weight of every shared CA in the interface fit, from its distance to the nearest binder heavy atom (in the design's own frame):
    1 up to ``full``, then a Gaussian ``exp(-((d - full) / decay)^2)`` (``decay`` is the distance beyond ``full`` at which the weight is
    1/e), never below ``floor`` (the far part of the target still keeps the transform from drifting, but counts for little)."""
    d = np.asarray(dist, dtype=float)
    w = np.where(d <= float(full), 1.0, np.exp(-(((d - float(full)) / max(float(decay), 1e-9)) ** 2)))
    return np.maximum(w, float(floor))


def interface_fit(P: np.ndarray, Q: np.ndarray, dist: np.ndarray, full: float = 8.0, decay: float = 6.0, floor: float = 0.01,
                  cutoff: float = 2.5, min_n: int = 12, cycles: int = IFACE_CYCLES
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, bool]:
    """Weighted Kabsch fit of the design's target CA ``P`` on the reference CA ``Q``, weighted around the binder.

    ``dist`` is the distance of every P to the nearest binder heavy atom. The first fit uses ``interface_weights``; residues that are
    near the binder (weight at least ``IFACE_NEAR_WEIGHT``) and deviate by more than ``cutoff`` under the transform are then left out and the
    rest refitted, as are residues farther away that deviate by more than ``IFACE_FAR_FACTOR`` times ``cutoff`` (a far part of the target that
    moved against the part around the binder); the rounds go on until the set stops changing (at most ``cycles``; a round that would leave
    fewer than ``min_n`` near residues is not taken). With fewer than ``min_n`` near residues to begin with there is nothing to weight around
    and the plain fit over every pair is returned. Returns (R, t, keep, weights, interface): ``interface`` is False for that fallback."""
    P, Q = np.asarray(P, dtype=float), np.asarray(Q, dtype=float)
    min_n = max(3, int(min_n))
    w = interface_weights(dist, full, decay, floor)
    near = w >= IFACE_NEAR_WEIGHT
    keep = np.ones(len(P), dtype=bool)
    if len(P) < 3 or int(near.sum()) < min_n:
        R, t, _ = S.kabsch(P, Q)
        return R, t, keep, w, False
    R, t, _ = S.kabsch(P, Q, w)
    for _ in range(max(1, int(cycles))):
        dev = np.linalg.norm(P @ R.T + t - Q, axis=1)
        new = ~(dev > np.where(near, float(cutoff), IFACE_FAR_FACTOR * float(cutoff)))
        if (new == keep).all() or int((new & near).sum()) < min_n:
            break
        keep = new
        R, t, _ = S.kabsch(P[keep], Q[keep], w[keep])
    return R, t, keep, w, True


# ----------------------------------------------------------------------------------------------------------------
# dataclasses
# ----------------------------------------------------------------------------------------------------------------
@dataclass
class RunSpec:
    """What the frame needs to know about one run: its hotspots and which construct residues its target carries.

    Numbers are construct numbers (the numbering BindCraft wrote). ``numbers`` None means the whole construct file.
    ``offset`` overrides the frame's construct-to-UniProt offset for a run built from a different construct file.
    """

    run_id: str
    hotspots: list[int]
    numbers: list[int] | None = None
    label: str = ""
    offset: int | None = None
    sequence: str | None = None          # one letter per entry of numbers, to check against the reference
    consistent: bool | None = None       # every structure of the run carries the same residue set (from the reader)


def _as_run_specs(runs: Any) -> list[RunSpec]:
    """Accept {run id: RunData | RunSpec | {"hotspots": [...], "numbers": [...]} | [hotspot, ...]} or a list of those."""
    if runs is None:
        raise FrameError("runs are required (run id -> RunData, RunSpec or hotspot list)")
    items = list(runs.items()) if isinstance(runs, Mapping) else [(getattr(r, "run_id", None), r) for r in runs]
    out: list[RunSpec] = []
    for rid, r in items:
        if isinstance(r, RunSpec):
            out.append(r)
        elif hasattr(r, "settings") and hasattr(r, "target"):   # a model.RunData
            t = r.target
            out.append(RunSpec(run_id=str(rid or r.run_id), hotspots=_ints(r.settings.hotspots),
                               numbers=_ints(t.numbers) or None, label=r.label or str(rid),
                               offset=t.uniprot_offset, sequence=t.sequence or None, consistent=t.consistent))
        elif isinstance(r, Mapping):
            out.append(RunSpec(run_id=str(rid), hotspots=_ints(r.get("hotspots", [])),
                               numbers=_ints(r["numbers"]) if r.get("numbers") else None,
                               label=str(r.get("label", rid)), offset=r.get("offset")))
        else:
            out.append(RunSpec(run_id=str(rid), hotspots=_ints(r), label=str(rid)))
    ids = [s.run_id for s in out]
    if len(set(ids)) != len(ids) or not all(ids):
        raise FrameError(f"run ids must be unique and non-empty: {ids}")
    return out


@dataclass
class Glycan:
    """One N-glycan of the reference structure: the sugar atoms attached to one Asn of one polymer chain."""

    owner: str                 # polymer chain the Asn belongs to
    site: int                  # Asn residue number (author numbering of the reference)
    resnames: list[str]        # sugar residue names in order, e.g. ["NAG", "NAG"]
    src: str                   # where the atoms are in the file: chain:residue per sugar, e.g. "D:1+D:2" or "A:3141"
    link_distance: float | None
    anchor: np.ndarray         # Asn ND2 position
    atoms: S.Atoms             # sugar heavy atoms; resnum is the 1-based sugar residue index inside this glycan

    def meta(self) -> dict[str, Any]:
        return {"owner": self.owner, "site": self.site, "resnames": list(self.resnames), "src": self.src,
                "link_distance": None if self.link_distance is None else round(float(self.link_distance), 3),
                "n_atoms": len(self.atoms)}


@dataclass
class Placement:
    """A design placed on the reference chain: the transform, how well the target fits, residue by residue."""

    run_id: str
    R: np.ndarray                  # 3x3, x' = x @ R.T + t
    t: np.ndarray
    rmsd: float                    # CA RMSD of all shared residues under this transform (place_rmsd)
    n: int                         # shared residues
    numbers: np.ndarray            # construct numbers of the shared residues, ascending
    offset: int                    # construct -> UniProt
    deviation: np.ndarray          # CA distance to the reference per shared residue after the transform (A)
    core_rmsd: float               # RMSD after a trimmed refit (residues deviating more than core_cutoff dropped)
    core_n: int
    core_numbers: np.ndarray       # construct numbers kept by the trimmed refit
    fit: str                       # the fit the transform comes from: "interface" (weighted around the binder), "all" (every shared CA) or "core" (trimmed)
    n_mismatch: int = 0            # shared residues whose type differs from the reference
    warnings: list[str] = field(default_factory=list)
    iface_rmsd: float = float("nan")    # CA RMSD under this transform over the residues near the binder (weight >= 0.5, within about 13 A); NaN without any
    iface_n: int = 0                    # how many residues that is
    global_rmsd: float = float("nan")   # CA RMSD the whole-target fit gives (the smallest RMSD over every shared residue any rigid transform can reach)
    fit_note: str = ""                  # why the fit is not the one asked for (the interface fit needs residues near the binder), else empty

    def apply(self, obj: Any) -> Any:
        """Transform an (n, 3) array or a structure.Atoms table into the reference frame."""
        if isinstance(obj, S.Atoms):
            return obj.transformed(self.R, self.t)
        return S.apply_transform(obj, self.R, self.t)

    @property
    def uniprot_numbers(self) -> np.ndarray:
        return self.numbers + self.offset

    def deviation_by_uniprot(self) -> dict[int, float]:
        return {int(u): float(d) for u, d in zip(self.uniprot_numbers, self.deviation)}

    @property
    def per_residue(self) -> dict[int, float]:
        """{UniProt number: CA deviation from the reference chain in A} for every shared residue (the construct number is
        ``number - offset``)."""
        return self.deviation_by_uniprot()

    @property
    def max_deviation(self) -> float:
        return float(self.deviation.max()) if len(self.deviation) else float("nan")

    def dropped_numbers(self) -> list[int]:
        """UniProt numbers the trimmed refit rejected (the mobile part of the target)."""
        kept = set(int(x) for x in self.core_numbers)
        return [int(n + self.offset) for n in self.numbers if int(n) not in kept]

    def rf(self) -> list[list[float]]:
        """[R(9), t(3)] row-major, the blobs' transform layout (x' = R x + t)."""
        return [self.R.reshape(-1).tolist(), self.t.tolist()]

    def fit_fields(self) -> dict[str, Any]:
        """How the design was placed and how well that fits where it matters: the fit used, the CA RMSD over the residues near the binder
        and over the whole target (what the whole-target fit gives), and the reason when the fit is not the one asked for. NaN becomes None."""
        return {"fit": self.fit,
                "iface_rmsd": round(self.iface_rmsd, 3) if math.isfinite(self.iface_rmsd) else None, "iface_n": int(self.iface_n),
                "global_rmsd": round(self.global_rmsd, 3) if math.isfinite(self.global_rmsd) else None, "fit_note": self.fit_note}

    def to_dict(self, residues: bool = False) -> dict[str, Any]:
        d = {"run": self.run_id, "rmsd": round(self.rmsd, 3), "n": int(self.n),
             "core_rmsd": round(self.core_rmsd, 3), "core_n": int(self.core_n), "max_dev": round(self.max_deviation, 2),
             "dropped": self.dropped_numbers(), "n_mismatch": int(self.n_mismatch), "R": _round(self.R.reshape(-1), 6),
             "t": _round(self.t, 4), **self.fit_fields()}
        if residues:
            d["deviation"] = {str(u): round(v, 2) for u, v in self.deviation_by_uniprot().items()}
        if self.warnings:
            d["warnings"] = list(self.warnings)
        return d


@dataclass
class ContextResult:
    """What the fragment hid, for one placed binder: metrics (the data contract names), the inspector's ctx and octx
    objects and the diagnostics behind them. ``metrics`` has exactly the keys of data_contract.md 8.5 that A2 owns."""

    run_id: str
    metrics: dict[str, Any]
    glycan: dict[str, float]
    glycan_min: dict[str, Any]
    ctx: dict[str, Any]
    octx: dict[str, Any]
    diag: dict[str, Any]
    config: dict[str, Any] = field(repr=False, default_factory=dict)

    def verdicts(self, config: Mapping[str, Any] | None = None) -> dict[str, dict[str, str]]:
        return verdicts(self.metrics, self.diag, config or self.config.get("checks"))

    def to_dict(self) -> dict[str, Any]:
        return {"metrics": {k: (None if v is None else (int(v) if isinstance(v, (int, np.integer)) else round(float(v), 3)))
                            for k, v in self.metrics.items()},
                "glycan": {k: round(v, 2) for k, v in self.glycan.items()}, "glycan_min": self.glycan_min,
                "ctx": self.ctx, "octx": self.octx, "diag": self.diag}


METRIC_KEYS = ("olig_min", "olig_n3", "olig_n4", "olig_n", "rest_min", "rest_n3", "rest_n4", "cut_n", "cov_n",
               "mem_excess")


def _partner_words(diag: Mapping[str, Any] | None) -> tuple[str, str, str]:
    """How a finding names what the binder runs into among the partner chains: ``(what, consequence, where)``. A clash with another subunit
    of the assembly means the assembled protein has no room for the binder there; a clash with a different molecule that was counted as a
    partner means the binder competes with it, which matters only if that molecule belongs to the complex one cares about; a clash with
    the rest of the target's own chain is a different claim again (the ``rest`` part). Without ``diag`` (or from an older frame) the
    older wording stands."""
    roles = list((diag or {}).get("partner_roles") or [])
    if not roles and (diag or {}).get("partner_role"):
        roles = [str(diag["partner_role"])]
    pair = (diag or {}).get("partner_pair") or {}
    chain = str((pair.get("partner") or ["", ])[0] or "")
    near = f", chain {chain}" if chain and roles else ""
    if not roles:
        return "the other subunit", "", ""
    if set(roles) == {"subunit"}:
        return "another subunit of the assembly", ": the binder would not fit on the assembled protein", near
    if "subunit" in roles:
        return "another subunit of the assembly and a different molecule of the complex", "", near
    if set(roles) == {"other_molecule"}:
        return "a different molecule of the complex", ": the binder would compete with it for this site", near
    return "a chain outside the assembly that was counted as a partner", "", near


def verdicts(metrics: Mapping[str, Any], diag: Mapping[str, Any] | None = None,
             checks: Mapping[str, Any] | None = None) -> dict[str, dict[str, str]]:
    """Pass, warn or fail of the Context check parts and the membrane flag, with a plain-language reason.

    The rules are profile.checks.context (docs/CONTRACTS.md 2.5): partner fails at olig_n3 >= fail_n3 and warns below
    warn_min, at any atom within 3 A or at warn_buried buried contacted residues; rest likewise without the buried
    part; cut and cover only warn. A missing metric takes no part (a missing metric never fails a design). ``overall``
    is the worst of partner, rest, cut and cover; ``membrane`` is a separate flag (fail above the plane)."""
    c = copy.deepcopy(DEFAULT_CONFIG["checks"])
    for k, v in (checks or {}).items():
        c.setdefault(k, {}).update(v)
    diag = diag or {}
    out: dict[str, dict[str, str]] = {}

    def has(k: str) -> bool:
        return metrics.get(k) is not None

    if has("olig_n3") and has("olig_min"):
        n3, dm, n4, nb = metrics["olig_n3"], metrics["olig_min"], metrics.get("olig_n4") or 0, metrics.get("olig_n") or 0
        p = c["partner"]
        what, claim, near = _partner_words(diag)
        if n3 >= p["fail_n3"]:
            out["partner"] = {"level": "fail", "reason": f"{n3} binder atoms within 3 A of {what} (closest {dm:.1f} A{near}){claim}"}
        elif n3 > 0 or dm < p["warn_min"] or nb >= p["warn_buried"]:
            why = []
            if n3 > 0 or dm < p["warn_min"]:
                why.append(f"closest binder atom {dm:.1f} A from {what} ({n4} within 4 A{near})")
            if nb >= p["warn_buried"]:
                why.append(f"{nb} contacted residues are buried at the interface with {what}")
            out["partner"] = {"level": "warn", "reason": "; ".join(why)}
        else:
            out["partner"] = {"level": "pass", "reason": f"clear of {what} ({dm:.1f} A)"}
    if has("rest_n3") and has("rest_min"):
        n3, dm = metrics["rest_n3"], metrics["rest_min"]
        r = c["rest"]
        if n3 >= r["fail_n3"]:
            out["rest"] = {"level": "fail", "reason": f"{n3} binder atoms within 3 A of the rest of the target chain, the part the construct left "
                                                      f"out (closest {dm:.1f} A): the full-length protein has material where this binder sits"}
        elif n3 > 0 or dm < r["warn_min"]:
            out["rest"] = {"level": "warn", "reason": f"closest binder atom {dm:.1f} A from the rest of the target chain, the part the construct left out"}
        else:
            out["rest"] = {"level": "pass", "reason": f"clear of the rest of the target chain ({dm:.1f} A)"}
    if has("cut_n"):
        n, k = metrics["cut_n"], c["cut"]
        near = list(diag.get("cut_contacts") or [])
        if n >= k["warn_n"]:
            out["cut"] = {"level": "warn", "reason": f"{n} contacted residues within {k['k']} of a cut end ({', '.join(map(str, near))})"}
        else:
            out["cut"] = {"level": "pass", "reason": f"{n} contacted residues within {k['k']} of a cut end"}
    if has("cov_n"):
        n, k = metrics["cov_n"], c["cover"]
        cov = [r[0] for r in diag.get("covered") or []]
        if n >= k["warn_n"]:
            out["cover"] = {"level": "warn", "reason": f"{n} contacted residues are buried by the removed regions ({', '.join(map(str, cov))})"}
        else:
            out["cover"] = {"level": "pass", "reason": f"{n} contacted residues are buried by the removed regions"}
    order = {"pass": 0, "warn": 1, "fail": 2}
    parts = [out[k]["level"] for k in ("partner", "rest", "cut", "cover") if k in out]
    if parts:
        worst = max(parts, key=lambda lv: order[lv])
        why = [f"{k}: {out[k]['reason']}" for k in ("partner", "rest", "cut", "cover") if k in out and out[k]["level"] == worst]
        out["overall"] = {"level": worst, "reason": "; ".join(why)}
    if has("mem_excess"):
        m = metrics["mem_excess"]
        out["membrane"] = ({"level": "fail", "reason": f"reaches {m:.1f} A past the membrane-side plane"} if m > 0 else
                           {"level": "pass", "reason": f"clear of the membrane-side plane by {-m:.1f} A"})
    return out


# ----------------------------------------------------------------------------------------------------------------
# reading the reference structure
# ----------------------------------------------------------------------------------------------------------------
def _first_residue_per_number(atoms: S.Atoms) -> tuple[np.ndarray, list[int], list[str]]:
    """The atoms of the first residue (in file order) of every residue number, as a mask over ``atoms``; the numbers that more than one residue
    carries (an insertion code, alternate residue types); and the labels of the residues that are not the first ("150A")."""
    res = atoms.residues
    seen: set[int] = set()
    shared: set[int] = set()
    keep = np.zeros(len(res), dtype=bool)
    later: list[str] = []
    for i in range(len(res)):
        n = int(res.number[i])
        if n in seen:
            shared.add(n)
            later.append(f"{n}{res.icode[i]}" if str(res.icode[i]).strip() else f"{n} {res.name[i]}")
        else:
            seen.add(n)
            keep[i] = True
    return (keep[atoms.res_index] if len(atoms) else np.zeros(0, dtype=bool)), sorted(shared), later


def _find_glycans(sf: S.StructureFile) -> list[Glycan]:
    """N-glycans from the covalent links of the file: Asn ND2 to a sugar C1, then sugar to sugar, so a branched chain
    (NAG-NAG) is followed to its end. The sugar residues may sit in a branched chain (C, D, E) or be single ligands of
    the polymer chain (A:ligands); both are found through the links, never by chain name."""
    sugar: dict[tuple[str, int, str], tuple[str, int]] = {}
    for gname, kind in sf.kinds.items():
        if kind not in ("branched", "ligand"):
            continue
        r = sf.groups[gname].residues
        for i in range(len(r)):
            sugar[(str(r.chain[i]), int(r.number[i]), str(r.name[i]))] = (gname, i)
    adj: dict[tuple[str, int, str], set] = {k: set() for k in sugar}
    anchors: dict[tuple[str, int], dict[tuple[str, int, str], float]] = {}
    for ln in sf.links:
        if ln.kind != "covale":
            continue
        ka, kb = (ln.a.chain, ln.a.resnum, ln.a.resname), (ln.b.chain, ln.b.resnum, ln.b.resname)
        sa, sb = ka in sugar, kb in sugar
        if sa and sb:
            adj[ka].add(kb)
            adj[kb].add(ka)
        elif sa != sb:
            prot, sug = (ln.b, ka) if sa else (ln.a, kb)
            if prot.resname == "ASN" and prot.atom == "ND2":
                site = (str(prot.chain), int(prot.resnum))
                best = anchors.setdefault(site, {})
                d = ln.distance if ln.distance else 9.9
                best[sug] = min(d, best.get(sug, 99.0))
    out: list[Glycan] = []
    for (chain, site), linked in sorted(anchors.items()):
        if chain not in sf.groups or sf.kinds.get(chain) != "polymer":
            continue
        order: list[tuple[str, int, str]] = []
        queue = deque(sorted(linked))
        seen: set = set()
        while queue:
            k = queue.popleft()
            if k in seen:
                continue
            seen.add(k)
            order.append(k)
            queue.extend(sorted(adj[k] - seen))
        parts, names, src_bits = [], [], []
        for idx, k in enumerate(order, start=1):
            gname, ordinal = sugar[k]
            at = sf.groups[gname].select_residues(ordinals=[ordinal]).heavy()
            parts.append(S.Atoms(at.xyz, at.name, at.element, at.resname, np.full(len(at), idx), np.full(len(at), chain),
                                 at.bfactor, hetero=np.ones(len(at), dtype=bool)))
            names.append(k[2])
            src_bits.append(f"{k[0]}:{k[1]}")
        poly = sf.groups[chain]
        try:
            anchor = poly.xyz[poly.atom(site, "ND2")].copy()
        except KeyError:
            continue
        out.append(Glycan(chain, site, names, "+".join(src_bits), min(linked.values()), anchor, S.Atoms.concat(parts)))
    return out


def _parse_assembly(path: str) -> dict[str, Any]:
    """Assembly records of the mmCIF: _pdbx_struct_assembly, _pdbx_struct_assembly_gen, _pdbx_struct_oper_list."""
    import gemmi

    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except (ValueError, RuntimeError, OSError):         # a PDB-format reference has no assembly records
        block = gemmi.cif.Block("none")
    asm = []
    for row in block.find("_pdbx_struct_assembly.", ["id", "details", "method_details", "oligomeric_details", "oligomeric_count"]):
        asm.append({"id": gemmi.cif.as_string(row[0]), "details": gemmi.cif.as_string(row[1]),
                    "method": gemmi.cif.as_string(row[2]), "oligomeric_details": gemmi.cif.as_string(row[3]),
                    "oligomeric_count": gemmi.cif.as_string(row[4])})
    gen = []
    for row in block.find("_pdbx_struct_assembly_gen.", ["assembly_id", "oper_expression", "asym_id_list"]):
        gen.append({"assembly_id": gemmi.cif.as_string(row[0]), "oper_expression": gemmi.cif.as_string(row[1]),
                    "asym_ids": [a.strip() for a in gemmi.cif.as_string(row[2]).split(",") if a.strip()]})
    opers = {}
    for row in block.find("_pdbx_struct_oper_list.", ["id", "type"]):
        opers[gemmi.cif.as_string(row[0])] = gemmi.cif.as_string(row[1])
    st = gemmi.read_structure(str(path))
    label_asym = {}
    for ch in st[0]:
        for res in ch:
            if res.entity_id and res.subchain:
                label_asym.setdefault(ch.name, res.subchain)
                break
    return {"assemblies": asm, "generators": gen, "operators": opers, "label_asym": label_asym}


# ----------------------------------------------------------------------------------------------------------------
# packing structure tables into numpy arrays (frame.npz)
# ----------------------------------------------------------------------------------------------------------------
_ATOM_FIELDS = ("xyz", "name", "element", "resname", "resnum", "chain", "bfactor")


def _pack(prefix: str, atoms: S.Atoms, out: dict[str, np.ndarray]) -> None:
    for f in _ATOM_FIELDS:
        out[f"{prefix}.{f}"] = np.asarray(getattr(atoms, f))
    out[f"{prefix}.hetero"] = np.asarray(atoms.hetero)


def _unpack(prefix: str, arrays: Mapping[str, np.ndarray]) -> S.Atoms:
    g = {f: arrays[f"{prefix}.{f}"] for f in _ATOM_FIELDS}
    return S.Atoms(g["xyz"], g["name"], g["element"], g["resname"], g["resnum"], g["chain"], g["bfactor"],
                   hetero=arrays[f"{prefix}.hetero"])


# ----------------------------------------------------------------------------------------------------------------
# environment atoms (partner, rest of the chain) as flat arrays with a lazy KD-tree
# ----------------------------------------------------------------------------------------------------------------
class _Env:
    """Atoms a binder must stay clear of. res is the residue number a hit is reported against (for a sugar the Asn it
    hangs on); chain is the polymer chain it belongs to."""

    def __init__(self, xyz, res, resname, name, sugar, chain):
        self.xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
        self.res = np.asarray(res, dtype=np.int64)
        self.resname = np.asarray(resname, dtype="<U3")
        self.name = np.asarray(name, dtype="<U4")
        self.sugar = np.asarray(sugar, dtype=bool)
        self.chain = np.asarray(chain, dtype="<U4")
        self._tree: cKDTree | None = None

    def __len__(self) -> int:
        return len(self.xyz)

    @property
    def tree(self) -> cKDTree:
        if self._tree is None:
            self._tree = cKDTree(self.xyz)
        return self._tree

    @classmethod
    def build(cls, polymer: S.Atoms, mask: np.ndarray, glycans: Sequence[Glycan]) -> "_Env":
        sel = polymer.select(mask)
        xyz = [sel.xyz]
        res = [sel.resnum]
        rn = [sel.resname]
        nm = [sel.name]
        sg = [np.zeros(len(sel), dtype=bool)]
        ch = [sel.chain]
        for g in glycans:
            n = len(g.atoms)
            xyz.append(g.atoms.xyz)
            res.append(np.full(n, g.site))
            rn.append(g.atoms.resname)
            nm.append(g.atoms.name)
            sg.append(np.ones(n, dtype=bool))
            ch.append(np.full(n, g.owner))
        return cls(np.vstack(xyz), np.concatenate(res), np.concatenate(rn), np.concatenate(nm), np.concatenate(sg),
                   np.concatenate(ch))


# ----------------------------------------------------------------------------------------------------------------
# the frame
# ----------------------------------------------------------------------------------------------------------------
class TargetFrame:
    """The reference chain, its dimer partner, glycans, annotations and per-run constructs (see docs/target_frame.md).

    ``meta`` is JSON-able, ``arrays`` numpy; nothing else is state. Build one with :func:`build_frame` or read one with
    :meth:`load`."""

    def __init__(self, meta: dict[str, Any], arrays: dict[str, np.ndarray]):
        self.meta = meta
        self.arrays = arrays
        self.config: dict[str, Any] = meta["config"]
        self._cache: dict[str, Any] = {}
        n = meta["numbering"]
        self.accession: str = n.get("accession") or ""      # "" when no UniProt accession was given
        self.length: int = int(n["length"])
        self.offset: int = int(n["offset"])                 # construct -> reference (UniProt) numbering
        self.native_range: list[int] = list(n.get("native_range") or [1, self.length])
        self.run_ids: list[str] = [r["run_id"] for r in meta["runs"]]
        self._runs = {r["run_id"]: r for r in meta["runs"]}

    @property
    def reference_given(self) -> bool:
        """False when the frame was built without a full-length reference: the construct is then the frame, nothing outside it is known and
        every check that needs the rest of the protein, the other subunits, sugars or a membrane is off. A frame saved before this field
        existed always had a reference."""
        return bool(self.meta["reference"].get("given", True))

    @property
    def reference_class(self) -> dict[str, Any]:
        """``classify_reference`` of the reference this frame was built from (``kind`` none without one). A frame saved before the
        classification existed reads as ``unknown``, which is treated as predicted."""
        c = self.meta["reference"].get("class")
        if c:
            return c
        return ({"kind": "unknown", "treated_as": "predicted", "method": None, "resolution": None,
                 "evidence": ["the frame predates the classification"]} if self.reference_given else classify_reference(None))

    # ---- persistence -----------------------------------------------------------------------------------------
    def save(self, directory: str | Path) -> Path:
        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        (d / FRAME_JSON).write_text(json.dumps(self.meta, indent=1, sort_keys=True, default=_json_default) + "\n")
        np.savez_compressed(d / FRAME_NPZ, **self.arrays)
        return d

    @classmethod
    def load(cls, directory: str | Path) -> "TargetFrame":
        d = Path(directory)
        try:
            meta = json.loads((d / FRAME_JSON).read_text())
            with np.load(d / FRAME_NPZ, allow_pickle=False) as z:
                arrays = {k: z[k] for k in z.files}
        except (OSError, ValueError) as ex:
            raise FrameError(f"cannot read a target frame from {d}: {ex}") from ex
        if meta.get("schema") != FRAME_SCHEMA:
            raise FrameError(f"{d}: schema {meta.get('schema')!r}, expected {FRAME_SCHEMA!r}")
        return cls(meta, arrays)

    # ---- reference tables ------------------------------------------------------------------------------------
    @property
    def chain_a(self) -> S.Atoms:
        if "chain_a" not in self._cache:
            self._cache["chain_a"] = _unpack("A", self.arrays)
        return self._cache["chain_a"]

    @property
    def numbers(self) -> np.ndarray:
        """Residue numbers (UniProt) of the reference chain, file order."""
        return self.arrays["res.number"]

    def _index(self, u: int) -> int:
        idx = self._cache.get("index")
        if idx is None:
            idx = {int(n): i for i, n in enumerate(self.numbers)}
            self._cache["index"] = idx
        return idx[int(u)]

    def has_residue(self, u: int) -> bool:
        try:
            self._index(u)
            return True
        except KeyError:
            return False

    @property
    def ambiguous_numbers(self) -> set[int]:
        """Residue numbers that more than one residue of the reference chain carries (insertion codes, alternate residue types). Designs are
        matched to the reference by number, which cannot tell such residues apart. The frame keeps the first residue of each such number in
        its tables (``numbers``, the per-residue arrays, ``chain_a``) and so in the target structure and every per-residue list of the
        blobs; the others are in ``dropped_atoms`` and count only as part of the rest of the chain."""
        if "ambiguous" not in self._cache:
            u, c = np.unique(self.numbers, return_counts=True)
            self._cache["ambiguous"] = {int(x) for x in u[c > 1]} | {int(x) for x in self.meta["numbering"].get("ambiguous", [])}
        return self._cache["ambiguous"]

    @property
    def dropped_atoms(self) -> S.Atoms:
        """Heavy atoms of the reference residues that share their number with an earlier residue (150A next to 150): not part of the residue
        tables, but part of the chain."""
        if "dropped" not in self._cache:
            self._cache["dropped"] = _unpack("X", self.arrays) if "X.xyz" in self.arrays else S.Atoms.empty()
        return self._cache["dropped"]

    @property
    def ref_ca(self) -> dict[int, np.ndarray]:
        """{UniProt number: CA position} of the reference chain. A number that two residues share (``ambiguous_numbers``) is left out, so that
        the placement never pairs a design residue with the wrong one."""
        if "ref_ca" not in self._cache:
            amb = self.ambiguous_numbers
            self._cache["ref_ca"] = {int(n): c for n, c in zip(self.numbers, self.arrays["res.ca"])
                                     if not np.isnan(c[0]) and int(n) not in amb}
        return self._cache["ref_ca"]

    @property
    def ref_one(self) -> dict[int, str]:
        """{UniProt number: one-letter code} of the reference chain."""
        if "ref_one" not in self._cache:
            self._cache["ref_one"] = {int(n): str(a) for n, a in zip(self.numbers, self.arrays["res.one"])}
        return self._cache["ref_one"]

    @property
    def uniprot_sequence(self) -> str:
        """The UniProt sequence of the target; empty when no accession was given or its record was not read."""
        return self.meta["numbering"].get("uniprot_sequence") or ""

    @property
    def partner_chains(self) -> list[str]:
        return list(self.meta["partner_chains"])

    @property
    def composition(self) -> dict[str, Any] | None:
        """What the reference file holds and what was done with each of its chains (``meta["composition"]``, see
        :func:`analyse_reference`); None for a frame without a reference or one saved before the composition was kept."""
        return self.meta.get("composition")

    def partner_role(self, chain: str) -> str:
        """What a partner is: ``subunit`` (the same protein as the target, in its assembly), ``other_molecule`` (a different molecule
        that was named as a partner) or ``outside_assembly`` (a chain of the file that the assembly does not hold, named as a partner).
        A frame saved before roles existed counted only subunits, so its partners are subunits."""
        comp = self.composition
        if comp:
            for row in comp["chains"]:
                if row["chain"] == chain:
                    return str(row["role"])
        return "subunit" if chain in self.partner_chains else ""

    def partner_label(self, chain: str) -> str:
        """The words that name a partner chain on hover: which kind of partner it is and where it comes from."""
        ref = self.meta["reference"]
        where = f"{ref['pdb']} assembly {ref['assembly']}" if ref.get("assembly") else str(ref["pdb"])
        row = next((r for r in (self.composition or {}).get("chains", []) if r["chain"] == chain), None) or {}
        role = self.partner_role(chain)
        if role == "other_molecule":
            return f"different molecule (chain {chain} of {where}" + (f": {row['name']}" if row.get("name") else "") + ")"
        if role == "outside_assembly":
            return f"neighbour outside the assembly (chain {chain} of {ref['pdb']}), counted because it was named"
        if row.get("generated"):
            return f"other subunit (symmetry copy {chain} of chain {row.get('copy_of')}, {where})"
        return f"other subunit (chain {chain} of {where})"

    @property
    def glycans(self) -> list[Glycan]:
        if "glycans" not in self._cache:
            out = []
            for k, gm in enumerate(self.meta["glycans"]):
                atoms = _unpack(f"G{k}", self.arrays)
                out.append(Glycan(gm["owner"], int(gm["site"]), list(gm["resnames"]), gm["src"], gm["link_distance"],
                                  self.arrays[f"G{k}.anchor"], atoms))
            self._cache["glycans"] = out
        return self._cache["glycans"]

    # ---- numbering -------------------------------------------------------------------------------------------
    def _offset_of(self, run_id: str | None) -> int:
        return self.offset if run_id is None else int(self._runs[run_id]["offset"])

    def construct_to_uniprot(self, i: int, run_id: str | None = None) -> int:
        return int(i) + self._offset_of(run_id)

    def uniprot_to_construct(self, u: int, run_id: str | None = None) -> int | None:
        """Construct number of a UniProt residue, or None when it is outside the construct file (or, with a run id,
        outside that run's target)."""
        i = int(u) - self._offset_of(run_id)
        if run_id is not None:
            return i if int(u) in self._run_present(run_id) else None
        c = self.meta["numbering"]
        return i if c["construct_first"] <= i <= c["construct_last"] else None

    def residue(self, u: int) -> str:
        """One-letter code of a UniProt residue (from the UniProt sequence, else the reference chain)."""
        seq = self.uniprot_sequence
        if seq and 1 <= int(u) <= len(seq):
            return seq[int(u) - 1]
        return str(self.arrays["res.one"][self._index(u)])

    def label(self, u: int) -> str:
        """"R12": one-letter code and residue number (UniProt numbering)."""
        return f"{self.residue(u)}{int(u)}"

    def labels(self, us: Iterable[int]) -> list[str]:
        return [self.label(u) for u in us]

    def assign(self, runs: Any) -> dict[str, Any]:
        """Set uniprot_offset and uniprot_accession on every run's TargetConstruct (what the readers leave empty).
        Returns {run id: {"offset", "n_present", "sequence_matches_reference"}}."""
        out = {}
        items = runs.items() if isinstance(runs, Mapping) else [(r.run_id, r) for r in runs]
        for rid, run in items:
            if rid not in self._runs:
                raise FrameError(f"run {rid!r} is not in this frame (runs: {self.run_ids})")
            t = run.target
            t.uniprot_offset = self._offset_of(rid)
            t.uniprot_accession = self.accession or None
            match = None
            if t.numbers and t.sequence:
                match = all(self.residue(n + t.uniprot_offset) == aa for n, aa in zip(t.numbers, t.sequence))
            out[rid] = {"offset": t.uniprot_offset, "n_present": len(t.numbers), "sequence_matches_reference": match}
        return out

    # ---- runs --------------------------------------------------------------------------------------------------
    def run_info(self, run_id: str) -> dict[str, Any]:
        try:
            return self._runs[run_id]
        except KeyError:
            raise FrameError(f"run {run_id!r} is not in this frame (runs: {self.run_ids})") from None

    def _run_present(self, run_id: str) -> set[int]:
        key = f"present.{run_id}"
        if key not in self._cache:
            self._cache[key] = _expand(self.run_info(run_id)["resi_ranges"])
        return self._cache[key]

    def hotspots(self, run_id: str | None = None) -> list[int]:
        if run_id is not None:
            return list(self.run_info(run_id)["hotspots"])
        return sorted({h for r in self.meta["runs"] for h in r["hotspots"]})

    def epitope(self, run_id: str | None = None) -> list[int]:
        if run_id is not None:
            return list(self.run_info(run_id)["epitope"])
        return sorted({h for r in self.meta["runs"] for h in r["epitope"]})

    def union_residues(self) -> list[int]:
        """UniProt numbers of the union of the runs' constructs (target.resi)."""
        return sorted(set().union(*(self._run_present(r) for r in self.run_ids)))

    def union_array(self, values: Mapping[int, Any], fill: Any = None, numbering: str = "construct") -> list[Any]:
        """Per-entry list over target.resi from {residue number: value}. numbering "construct" (default) takes construct
        numbers, as per-design analysis produces them; "uniprot" takes UniProt numbers. Residues without a value get
        ``fill`` (a design's run may lack part of the union)."""
        out = []
        for u in self.union_residues():
            key = u - self.offset if numbering == "construct" else u
            out.append(values.get(key, fill))
        return out

    def hotspots_touched(self, run_id: str, contacted_uniprot: Iterable[int]) -> list[int]:
        """design.hot: contacted residues that are hotspots of the design's own run."""
        return sorted(set(int(c) for c in contacted_uniprot) & set(self.hotspots(run_id)))

    # ---- environments ------------------------------------------------------------------------------------------
    def _partner_env(self) -> _Env | None:
        if "partner" not in self._cache:
            parts = [_unpack(f"P_{c}", self.arrays) for c in self.partner_chains]
            if not parts:
                self._cache["partner"] = None
            else:
                poly = S.Atoms.concat(parts)
                self._cache["partner"] = _Env.build(poly, np.ones(len(poly), dtype=bool),
                                                    [g for g in self.glycans if g.owner in self.partner_chains])
        return self._cache["partner"]

    def _rest_env(self, run_id: str) -> _Env | None:
        key = f"rest.{run_id}"
        if not self.reference_given:
            return None             # the stretches a run deleted from the construct are not "the rest of the protein": nothing outside it is known
        if key not in self._cache:
            present = self._run_present(run_id)
            A, X = self.chain_a, self.dropped_atoms
            poly = S.Atoms.concat([A, X]) if len(X) else A                  # a residue with an insertion code is part of the rest of the chain too
            mask = ~np.isin(poly.resnum, np.fromiter(present, dtype=np.int64))
            gl = [g for g in self.glycans if g.owner == self.meta["chain"] and g.site not in present]
            env = _Env.build(poly, mask, gl)
            self._cache[key] = env if len(env) else None
        return self._cache[key]

    # ---- membrane ------------------------------------------------------------------------------------------------
    @property
    def membrane(self) -> dict[str, Any]:
        m = self.meta["membrane"]
        return {"centroid": np.asarray(m["centroid"], dtype=float), "axis": np.asarray(m["axis"], dtype=float),
                "offset": float(m["offset"])}

    # ---- placement -------------------------------------------------------------------------------------------------
    def place(self, design: Any, run_id: str, fit: str | None = None, n_binder_chains: int | None = None) -> Placement:
        """``place_design`` with this frame; ``fit`` defaults to the frame's ``geometry.fit``."""
        return place_design(self, design, run_id, fit=fit, n_binder_chains=n_binder_chains)

    # ---- context checks ---------------------------------------------------------------------------------------------
    def context(self, binder: Any, run_id: str, contacted: Iterable[int] = (), *, config: Mapping[str, Any] | None = None
                ) -> ContextResult:
        """Context of a binder already placed in the reference frame.

        ``binder``: structure.Atoms (heavy atoms; hydrogens are dropped) or an (n, 3) array. ``contacted``: UniProt numbers
        of the target residues the binder contacts (the interface module's tc, converted). Residues outside the run's
        target are ignored. Returns the data contract metrics (olig_*, rest_*, cut_n, cov_n, mem_excess), the glycan
        distances, the inspector objects ctx and octx and diagnostics."""
        cfg = copy.deepcopy(self.config)
        for sect, vals in (config or {}).items():
            cfg.setdefault(sect, {}).update(vals)
        geo, chk = cfg["geometry"], cfg["checks"]
        run = self.run_info(run_id)
        present = self._run_present(run_id)
        if isinstance(binder, S.Atoms):
            b = binder.heavy() if binder.is_hydrogen.any() else binder
            xyz, bres, bname, bresn = b.xyz, b.resnum, b.name, b.resname
        else:
            xyz = np.asarray(binder, dtype=float).reshape(-1, 3)
            bres = np.zeros(len(xyz), dtype=np.int64)
            bname = np.array([""] * len(xyz))
            bresn = np.array([""] * len(xyz))
        if not len(xyz):
            raise FrameError("no binder atoms to check")
        contacted = sorted(int(c) for c in set(contacted) if int(c) in present)
        n_close = [float(x) for x in geo["n_close"]]
        r_clash, (r3, r4, r6) = float(geo["clash"]), n_close

        metrics: dict[str, Any] = {k: None for k in METRIC_KEYS}
        diag: dict[str, Any] = {"contacted": contacted}
        if not self.reference_given:
            # Without a reference the construct is the whole frame. Every metric here compares the binder with something the construct does not
            # contain (the other subunit, the rest of the chain, sugars, a membrane), and cut_n and cov_n need to know which ends were cut and
            # what the rest would bury. They are absent, not zero: a zero would read as "checked, clear".
            diag["glycan_min"] = {"site": None, "d": None}
            return ContextResult(run_id, metrics, {}, {"site": None, "d": None}, {}, {}, diag, cfg)

        def nearest(env: _Env | None):
            if env is None:
                return None
            d, ix = env.tree.query(xyz)
            j = int(np.argmin(d))
            return d, ix, j

        # partner: the other subunit and its sugars
        pe = self._partner_env()
        hit = nearest(pe)
        octx: dict[str, Any] = {}
        if hit is not None:
            d, ix, j = hit
            buried = self.partner_buried()
            olig_buried = [c for c in contacted if buried.get(c, 0.0) >= float(geo["olig_min_area"])]
            metrics.update(olig_min=float(d[j]), olig_n3=int((d < r3).sum()), olig_n4=int((d < r4).sum()),
                           olig_n=len(olig_buried))
            diag["partner_pair"] = {"binder": [int(bres[j]), str(bresn[j]), str(bname[j])],
                                    "partner": [str(pe.chain[ix[j]]), int(pe.res[ix[j]]), str(pe.resname[ix[j]]),
                                                str(pe.name[ix[j]])], "d": round(float(d[j]), 2)}
            diag["partner_role"] = self.partner_role(str(pe.chain[ix[j]]))
            diag["partner_roles"] = sorted({self.partner_role(str(c)) for c in set(pe.chain[ix[d < r4]].tolist())}) or [diag["partner_role"]]
            diag["partner_clash"] = int((d < r_clash).sum())
            diag["dimer_buried_contacts"] = olig_buried
            octx = {"olig_min": round(float(d[j]), 2), "olig_n3": metrics["olig_n3"], "olig_n4": metrics["olig_n4"],
                    "olig_buried": olig_buried}
        # a homodimer carries two copies of the binder: the copy on the other protomer is the binder moved by the
        # operator that maps chain A onto chain B; the two copies must not overlap (information only, no verdict)
        tf = (self.meta.get("dimer") or {}).get("two_fold") or {}
        if tf.get("R") and tf.get("t"):
            twin = xyz @ np.asarray(tf["R"], dtype=float).reshape(3, 3).T + np.asarray(tf["t"], dtype=float)
            ds = cKDTree(xyz).query(twin)[0]
            diag["sym_min"] = round(float(ds.min()), 2)
            octx["sym_min"], octx["sym_n3"] = diag["sym_min"], int((ds < r3).sum())
        # rest of the chain: residues this run's target lacks, with the sugars attached to them
        re_env = self._rest_env(run_id)
        hit = nearest(re_env)
        ctx: dict[str, Any] = {}
        areas_all: list[list[float]] = []
        covered: list[list[float]] = []
        if hit is not None:
            d, ix, j = hit
            metrics.update(rest_min=float(d[j]), rest_n3=int((d < r3).sum()), rest_n4=int((d < r4).sum()))
            sugar = bool(re_env.sugar[ix[j]])
            pair_res = int(re_env.res[ix[j]])
            pair_resn = "ASN" if sugar else str(re_env.resname[ix[j]])
            diag["rest_pair"] = {"binder": [int(bres[j]), str(bresn[j]), str(bname[j])],
                                 "rest": [int(re_env.res[ix[j]]), str(re_env.resname[ix[j]]), str(re_env.name[ix[j]])],
                                 "glycan": sugar, "d": round(float(d[j]), 2)}
            diag["rest_clash"] = int((d < r_clash).sum())
            diag["rest_n_by_region"] = self._region_counts(run_id, re_env, ix[d < r4])
            ctx = {"ref": self.meta["reference"]["pdb"], "removed": [list(r) for r in run["removed_native"]],
                   "nearest": round(float(d[j]), 2), "n3": metrics["rest_n3"], "n4": metrics["rest_n4"],
                   "n6": int((d < r6).sum()), "pair": [int(bres[j]), pair_res], "pair_resn": pair_resn}
            if sugar:
                ctx["pair_glycan"] = str(re_env.resname[ix[j]])
        # residues the removed regions bury (their area with and without the rest of the chain present)
        area = self.covered_area(run_id)
        min_area = float(chk["cover"]["min_area"])
        areas_all = [[int(u), round(float(a), 1)] for u, a in sorted(area.items()) if a > min_area]
        covered = [[u, a] for u, a in areas_all if u in set(contacted)]
        metrics["cov_n"] = len(covered)
        diag["covered"] = covered
        # in the whole assembly (other subunit and all sugars included) rather than the rest of the chain alone
        asm_loss = self.assembly_loss(run_id)
        diag["buried_in_assembly"] = [[c, round(asm_loss[c], 1)] for c in contacted if asm_loss.get(c, 0.0) > min_area]
        # cut ends
        cuts = [int(p) for p in run["cut_points"]]
        k = int(chk["cut"]["k"])
        near_cut = [c for c in contacted if any(abs(c - p) <= k for p in cuts)]
        metrics["cut_n"] = len(near_cut)
        diag["cut_contacts"] = near_cut
        ctx.setdefault("ref", self.meta["reference"]["pdb"])
        ctx.setdefault("removed", [list(r) for r in run["removed_native"]])
        ctx["covered"] = covered
        ctx["covered_all"] = areas_all
        ctx["cut"] = {"points": cuts, "near": near_cut}
        # membrane (only when a transmembrane region is known: without one the plane is a placeholder and no distance is reported)
        metrics["mem_excess"] = membrane_excess(xyz, self.membrane) if self.meta["membrane"].get("known", True) else None
        # glycans of the construct (sites inside this run's target)
        gly: dict[str, float] = {}
        for site in self.meta["glycan_sites_in_construct"]:
            if site not in present:
                continue
            gs = [g for g in self.glycans if g.owner == self.meta["chain"] and g.site == site]
            if gs:
                pts = np.vstack([g.atoms.xyz for g in gs])
            else:
                pts = self.arrays["glycan.anchor_xyz"][list(self.meta["glycan_sites_in_construct"]).index(site)][None, :]
            gly[str(site)] = float(cKDTree(pts).query(xyz)[0].min())
        gmin = {"site": None, "d": None}
        if gly:
            s = min(gly, key=gly.get)
            gmin = {"site": int(s), "d": round(gly[s], 2)}
        diag["glycan_min"] = gmin
        return ContextResult(run_id, metrics, gly, gmin, ctx, octx, diag, cfg)

    def _region_counts(self, run_id: str, env: _Env, ix: np.ndarray) -> dict[str, int]:
        """Binder atoms within 4 A of the rest, by the removed region they touch: N-terminal side of the construct,
        C-terminal side, or a loop deleted inside it; sugars are counted as 'glycan'."""
        run = self.run_info(run_id)
        first, last = run["first"], run["last"]
        deleted = _expand(run["removed"])
        out: Counter = Counter()
        for r, sugar in zip(env.res[ix], env.sugar[ix]):
            if sugar:
                out["glycan"] += 1
            elif int(r) < first:
                out["before"] += 1
            elif int(r) > last:
                out["after"] += 1
            elif int(r) in deleted:
                out["deleted_loop"] += 1
            else:
                out["other"] += 1
        return dict(out)

    def partner_buried(self) -> dict[int, float]:
        """{residue: area A^2} buried by the other subunit (and its sugars) for every reference-chain residue above
        olig_min_area."""
        if "buried" not in self._cache:
            thr = float(self.config["geometry"]["olig_min_area"])
            b = self.arrays["res.buried_partner"]
            self._cache["buried"] = {int(n): float(v) for n, v in zip(self.numbers, b) if v >= thr}
        return self._cache["buried"]

    def covered_area(self, run_id: str) -> dict[int, float]:
        """{residue: area A^2} a run's residues lose to the removed regions (SASA of the run's target alone minus with
        the rest of the chain and its sugars), all residues with a positive loss."""
        key = f"cov.{run_id}"
        if key not in self._cache:
            alone = self.arrays[f"run.{run_id}.sasa_alone"]
            rest = self.arrays[f"run.{run_id}.sasa_rest"]
            loss = alone - rest
            self._cache[key] = {int(n): float(v) for n, v, a in zip(self.numbers, loss, alone)
                                if not np.isnan(a) and v > 1e-6}
        return self._cache[key]

    def assembly_loss(self, run_id: str) -> dict[int, float]:
        """{residue: area A^2} a run's residues lose in the whole assembly (the other subunit, the rest of the chain and every
        sugar present) compared with the run's target alone, all residues with a positive loss. A superset of
        ``covered_area`` and ``partner_buried`` together."""
        key = f"asm.{run_id}"
        if key not in self._cache:
            alone = self.arrays[f"run.{run_id}.sasa_alone"]
            loss = alone - self.arrays["res.sasa_dimer"]
            self._cache[key] = {int(n): float(v) for n, v, a in zip(self.numbers, loss, alone) if not np.isnan(a) and v > 1e-6}
        return self._cache[key]

    # ---- blocks for the blobs ---------------------------------------------------------------------------------------
    def _selected(self, runs: Any) -> list[str]:
        if runs is None:
            return list(self.run_ids)
        if isinstance(runs, str):
            return [runs]
        ids = [r if isinstance(r, str) else getattr(r, "run_id", None) for r in (runs.keys() if isinstance(runs, Mapping) else runs)]
        bad = [r for r in ids if r not in self._runs]
        if bad:
            raise FrameError(f"runs not in this frame: {bad} (frame runs: {self.run_ids})")
        return ids

    def run_blocks(self, runs: Any = None) -> dict[str, dict[str, Any]]:
        """{run id: runs[].target} (data_contract.md 3.3), UniProt numbering."""
        out = {}
        for rid in self._selected(runs):
            r = self._runs[rid]
            out[rid] = {"resi_ranges": [list(x) for x in r["resi_ranges"]], "removed": [list(x) for x in r["removed"]],
                        "first": r["first"], "last": r["last"], "n_resi": r["n_resi"], "cut_points": list(r["cut_points"]),
                        "input_offset": -int(r["offset"]), "hotspots": list(r["hotspots"]), "epitope": list(r["epitope"]),
                        "sequence_matches_reference": r["sequence_matches_reference"], "consistent": r["consistent"],
                        "reference_given": self.reference_given}
        return out

    def profile_frame(self) -> dict[str, Any]:
        """profile.frame as derived from UniProt, SIFTS and the reference structure (native_range, segments, glycan_sites,
        disulfides). Without an accession the segments are empty and the sites are those the reference structure shows."""
        a = self.meta["annotations"]
        return {"native_range": list(self.native_range), "segments": copy.deepcopy(a["segments"]),
                "glycan_sites": [g["site"] for g in a["glycosylation"]],
                "disulfides": [list(d["pair"]) for d in a["disulfides"]]}

    def _glycan_letter(self, site: int) -> str:
        """Chain id of a site's sugars in target_pdb and glycan_atoms: C, D, E ... in site order (A and B belong to the
        polymer chains), independent of which runs are selected."""
        return GLYCAN_CHAINS[list(self.meta["glycan_sites_in_construct"]).index(int(site))]

    def _glycan_block(self) -> tuple[list[int], dict[str, list[float]], dict[str, dict[str, Any]]]:
        """Glycan sites of the union construct, their ND2 anchors and the sugar atoms of the reference chain."""
        sites = sorted(set(self.union_residues()) & set(self.meta["glycan_sites_in_construct"]))
        xyz, atoms = {}, {}
        for site in sites:
            idx = list(self.meta["glycan_sites_in_construct"]).index(site)
            xyz[str(site)] = _round(self.arrays["glycan.anchor_xyz"][idx], 2)
            gs = [g for g in self.glycans if g.owner == self.meta["chain"] and g.site == site]
            if gs:
                g = gs[0]
                atoms[str(site)] = {"chain": self._glycan_letter(site), "resn": list(g.resnames),
                                    "el": g.atoms.element.tolist(), "xyz": _flat(g.atoms.xyz), "name": g.atoms.name.tolist(),
                                    "res": (g.atoms.resnum - 1).astype(int).tolist(),
                                    "src": f"{self.meta['reference']['pdb']} {g.src}"}
        return sites, xyz, atoms

    def target_block(self, runs: Any = None, name: str | None = None, site_residues: Iterable[int] | None = None) -> dict[str, Any]:
        """The campaign blob's ``target`` (data_contract.md 3.2) for the selected runs (default all), UniProt numbering,
        reference frame, coordinates rounded to 0.01 A and areas to 0.1 A^2. ``name`` defaults to the UniProt name.

        Site A1 is the intended epitope: the centre of the hotspots, its residues the hotspots and their neighbours. A BindCraft run
        may have no hotspot residue at all; when none of the selected runs records one, A1 is where the binders sit instead:
        ``site_residues`` (the residues the designs contact most, from the pose records) give the centre and the residues and the
        site says so in ``source``; without them the centre is the construct's and the normal is the axis the construct is thinnest
        along (``default_view_axis``). The normal is always a unit vector, since the pages orient their views by it."""
        ids = self._selected(runs)
        present = sorted(set().union(*(self._run_present(r) for r in ids)))
        idx = np.array([self._index(u) for u in present], dtype=np.int64)
        ca = self.arrays["res.ca"][idx]
        atoms, ares, asasa = self._atoms_union_for(present)
        gl_sites, gl_xyz, gl_atoms = self._glycan_block()
        gl_sites = [s for s in gl_sites if s in set(present)]
        hs = {r: list(self._runs[r]["hotspots"]) for r in ids}
        ep = {r: list(self._runs[r]["epitope"]) for r in ids}
        hot_all = sorted({h for v in hs.values() for h in v})
        epi_all = sorted({h for v in ep.values() for h in v})
        uni_ca = {int(n): c for n, c in zip(self.numbers, self.arrays["res.ca"])}
        centroid = ca.mean(axis=0)
        site: dict[str, Any] = {}
        if hot_all:
            center, site_res = np.mean([uni_ca[h] for h in hot_all], axis=0), epi_all
        else:
            site_res = sorted(int(r) for r in (site_residues or ()) if int(r) in set(present) and not np.isnan(uni_ca[int(r)][0]))
            center = np.mean([uni_ca[r] for r in site_res], axis=0) if site_res else centroid
            site = {"source": ("no hotspots recorded: the residues the designs contact most" if site_res else
                               "no hotspots recorded and no design contacts: the whole construct, viewed along its thinnest axis")}
        direction = center - centroid
        normal = direction / np.linalg.norm(direction) if float(np.linalg.norm(direction)) > 1e-6 else default_view_axis(ca)
        m = self.membrane
        t = {
            "name": name or self.meta["annotations"]["name"], "uniprot": self.accession, "construct": [present[0], present[-1]],
            "length": self.length,
            "resi": _ints(present), "resn": [self.residue(u) for u in present], "ss": "".join(self.arrays["res.ss"][idx].tolist()),
            "rsa": _round(self.arrays["res.rsa"][idx], 2), "ca": _flat(ca),
            "atoms": {"x": _flat(atoms.xyz), "res": ares.tolist(), "el": atoms.element.tolist(), "name": atoms.name.tolist(),
                      "sasa": _round(asasa, 1)},
            "hotspots": hot_all, "epitope": epi_all, "glycan": gl_sites, "glycan_xyz": {k: v for k, v in gl_xyz.items() if int(k) in gl_sites},
            "membrane": {"axis": _round(m["axis"], 4), "centroid": _round(m["centroid"], 2), "offset": round(m["offset"], 2),
                         "label": self.meta["membrane"]["label"], "source": self.meta["membrane"]["source"]},
            "sites": {"A1": {"center": _round(center, 2), "normal": _round(normal, 4), "residues": site_res, **site}},
            "olig": self._olig_block(present),
            "orig_resi": [u - self.offset for u in present],
            "numbering": {"scheme": f"UniProt {self.accession}" if self.accession else "reference structure numbering",
                          "input_offset": -self.offset, "input_label": "BindCraft input (construct numbering)"},
            "reference": {"pdb": self.meta["reference"]["pdb"] if self.reference_given else "", "chain": self.meta["chain"],
                          "assembly": self.meta["reference"]["assembly"], "uniprot": self.accession,
                          "label": self.meta["reference"]["label"], "given": self.reference_given},
            "unresolved": _complement(self.meta["annotations"]["observed"], present[0], present[-1]) if self.meta["annotations"]["observed"] else [],
            "hotspots_by_run": hs, "epitope_by_run": ep,
            "glycan_atoms": {k: v for k, v in gl_atoms.items() if int(k) in gl_sites},
            "rest": self._rest_block(ids),
            "dom": [self.domain_of(u) for u in present],
            "frame": self.profile_frame(),
        }
        return t

    def _atoms_union_for(self, present: list[int]) -> tuple[S.Atoms, np.ndarray, np.ndarray]:
        A = self.chain_a
        mask = np.isin(A.resnum, np.array(present, dtype=np.int64))
        sel = A.select(mask)
        pos = {int(n): i for i, n in enumerate(present)}
        return sel, np.array([pos[int(n)] for n in sel.resnum], dtype=np.int64), self.arrays["A.sasa"][mask]

    def domain_of(self, u: int) -> str:
        for s in self.meta["annotations"]["segments"]:
            if s["start"] <= u <= s["end"]:
                return s["label"]
        return ""

    def _olig_block(self, present: list[int]) -> dict[str, Any]:
        if not self.partner_chains:
            return {}
        pset = set(present)
        neigh = {}
        labels = {}
        for c in self.partner_chains:
            at = _unpack(f"P_{c}", self.arrays)
            ca = at.select((at.name == "CA") & (at.element == "C"))
            neigh[c] = _flat(ca.xyz)
            labels[c] = self.partner_label(c)
        ref = self.meta["reference"]
        return {"buried": {str(u): round(a, 1) for u, a in sorted(self.partner_buried().items()) if u in pset},
                "neigh": neigh, "chains": list(self.partner_chains), "labels": labels,
                "roles": {c: self.partner_role(c) for c in self.partner_chains},
                "source": f"{ref['pdb']} assembly {ref['assembly']}" if ref.get("assembly") else str(ref["pdb"]),
                "min_area": float(self.config["geometry"]["olig_min_area"])}

    def _rest_block(self, ids: list[str]) -> dict[str, Any]:
        """Rest-of-chain context: one trace per stretch of the reference chain that at least one run's target lacks."""
        if not self.reference_given:
            return {"ref": self.meta["reference"]["pdb"], "ranges": [], "n_res": 0, "traces": [], "not_in_reference": [],
                    "source": "no reference structure was given: nothing outside the construct is known"}
        sets = [self._run_present(r) for r in ids]
        union = set().union(*sets)
        lacking = set()
        for u in self.numbers:
            if any(int(u) not in s for s in sets):
                lacking.add(int(u))
        traces = []
        ca = {int(n): c for n, c in zip(self.numbers, self.arrays["res.ca"])}
        for a, b in _ranges(lacking):
            where = "before the construct" if b < min(union) else "after the construct" if a > max(union) else "inside the construct"
            who = [r for r, s in zip(ids, sets) if any(n not in s for n in range(a, b + 1))]
            label = f"{where} ({a}-{b})" if where != "inside the construct" else f"removed in {', '.join(who)} ({a}-{b})"
            traces.append({"label": label, "chain": self.meta["chain"], "ranges": [[a, b]],
                           "ca": _flat(np.array([ca[n] for n in range(a, b + 1)]))})
        outside = [r for r in _ranges(lacking) if r[1] < min(union) or r[0] > max(union)]
        return {"ref": f"{self.meta['reference']['pdb']} chain {self.meta['chain']}", "ranges": outside,
                "n_res": sum(b - a + 1 for a, b in outside), "traces": traces,
                "not_in_reference": copy.deepcopy(self.meta["annotations"]["unresolved"]),
                "source": "reference chain residues no run's target contains, placed from the crystal structure"}

    def target_pdb(self, runs: Any = None) -> str:
        """PDB text of the union construct of the reference chain (chain A, UniProt numbering, heavy atoms) with HELIX and
        SHEET records, the in-construct sugars as HETATM (chains C, D, ...) and END."""
        ids = self._selected(runs)
        present = sorted(set().union(*(self._run_present(r) for r in ids)))
        idx = np.array([self._index(u) for u in present], dtype=np.int64)
        atoms, _, _ = self._atoms_union_for(present)
        ss = "".join(self.arrays["res.ss"][idx].tolist())
        names = [AA3.get(self.arrays["res.one"][i], "UNK") for i in idx]
        lines = ss_records(ss, present, names, "A")                   # always chain A: the pages tell the binder (chain B) from the target by the chain id,
        lines += S.to_pdb_lines(atoms, 1, chain="A")                  # whatever the chain of the reference structure is called
        serial = len(atoms) + 1
        for site in self._glycan_block()[0]:
            if site not in set(present):
                continue
            gs = [x for x in self.glycans if x.owner == self.meta["chain"] and x.site == site]
            if gs:
                lines += S.to_pdb_lines(gs[0].atoms, serial, chain=self._glycan_letter(site))
                serial += len(gs[0].atoms)
        lines.append("END")
        return "\n".join(lines) + "\n"

    def inspector_block(self, runs: Any = None, site_residues: Iterable[int] | None = None) -> dict[str, Any]:
        """The target-side part of the inspector blob (data_contract.md 4): target_pdb, target_seq, target_resi, t_ss,
        hotspots, epitope, glycan, glycan_xyz, glycan_atoms, membrane, runs, olig, rest."""
        ids = self._selected(runs)
        t = self.target_block(ids, site_residues=site_residues)
        rb = self.run_blocks(ids)
        return {"target_pdb": self.target_pdb(ids), "target_seq": "".join(t["resn"]), "target_resi": t["resi"], "t_ss": t["ss"],
                "hotspots": t["hotspots"], "epitope": t["epitope"], "glycan": t["glycan"], "glycan_xyz": t["glycan_xyz"],
                "glycan_atoms": t["glycan_atoms"], "membrane": t["membrane"],
                "runs": [{"id": r, "label": self._runs[r]["label"], "hotspots": list(self._runs[r]["hotspots"]),
                          "target": rb[r]} for r in ids],
                "olig": t["olig"], "rest": t["rest"]}

    # ---- summaries ---------------------------------------------------------------------------------------------------
    def residue_table(self, run_id: str | None = None) -> list[dict[str, Any]]:
        """One row per residue of a run's target (default: the union construct): UniProt and construct number, residue,
        secondary structure and the SASA of the residue in four settings (A^2): the run's target alone (what BindCraft
        saw), the full reference chain, the chain with its sugars, and the dimer with all sugars."""
        present = sorted(self._run_present(run_id)) if run_id else self.union_residues()
        alone = self.arrays[f"run.{run_id}.sasa_alone"] if run_id else None
        rows = []
        for u in present:
            i = self._index(u)
            rows.append({"uniprot": u, "construct": u - self.offset, "aa": self.residue(u), "ss": str(self.arrays["res.ss"][i]),
                         "rsa": round(float(self.arrays["res.rsa"][i]), 3),
                         "sasa_target_alone": None if alone is None else round(float(alone[i]), 1),
                         "sasa_chain": round(float(self.arrays["res.sasa"][i]), 1),
                         "sasa_chain_glycans": round(float(self.arrays["res.sasa_chain_gly"][i]), 1),
                         "sasa_dimer": round(float(self.arrays["res.sasa_dimer"][i]), 1),
                         "buried_by_partner": round(float(self.arrays["res.buried_partner"][i]), 1)})
        return rows

    def summary(self) -> dict[str, Any]:
        n = self.meta["numbering"]
        return {"uniprot": self.accession or None, "reference": f"{self.meta['reference']['pdb']} chain {self.meta['chain']}",
                "offset": self.offset, "construct": [n["construct_first"], n["construct_last"]],
                "runs": {r: {"n_resi": v["n_resi"], "ranges": v["resi_ranges"], "hotspots": self.labels(v["hotspots"]),
                             "cut_points": v["cut_points"]} for r, v in self._runs.items()},
                "dimer": {k: self.meta["dimer"].get(k) for k in ("oligomeric_details", "bsa_total", "n_interface")},
                "warnings": list(self.meta.get("warnings", []))}


# ----------------------------------------------------------------------------------------------------------------
# building the frame
# ----------------------------------------------------------------------------------------------------------------
def _sasa_parts(parts: Sequence[S.Atoms], probe: float, points: int) -> list[np.ndarray]:
    """Per-atom SASA of several atom tables when all of them are present together."""
    parts = [p for p in parts]
    cat = S.Atoms.concat(parts)
    s = S.sasa_atoms(cat.xyz, cat.element, probe, points)
    out, pos = [], 0
    for p in parts:
        out.append(s[pos:pos + len(p)])
        pos += len(p)
    return out


def _glycan_atoms(gl: Sequence[Glycan]) -> S.Atoms:
    return S.Atoms.concat([g.atoms for g in gl]) if gl else S.Atoms.empty()


def _read_construct(path: str) -> S.Atoms:
    cx = S.load_complex(path, binder_chain=None)
    if not len(cx.target):
        raise FrameError(f"{path}: no atoms in chain A")
    return cx.target


def _vote_offset(items: Iterable[tuple[int, str]], target: Mapping[int, str]) -> tuple[int, int]:
    """Offset d for which ``target[number + d]`` equals the residue type at the most (number, type) pairs. The numbering
    of the two sequences is compared residue by residue, so gaps in either need no special handling. Returns (d, matches);
    an offset tied with another one is an error rather than a guess."""
    by_aa: dict[str, list[int]] = {}
    for n, aa in target.items():
        by_aa.setdefault(aa, []).append(int(n))
    votes: Counter = Counter()
    for n, aa in items:
        for tn in by_aa.get(str(aa), ()):
            votes[tn - int(n)] += 1
    if not votes:
        raise FrameError("the two sequences have no residue in common")
    best = votes.most_common(2)
    if len(best) > 1 and best[0][1] == best[1][1]:
        raise FrameError(f"numbering offset is ambiguous ({best[0]} and {best[1]})")
    return int(best[0][0]), int(best[0][1])


def _safe(fn, warnings: list[str], what: str, default=None):
    """Run a reference-data lookup; a missing cache entry or network failure becomes a warning, not an error."""
    from . import refdata as rd

    try:
        return fn()
    except (rd.RefDataError, KeyError, ValueError, OSError) as ex:
        warnings.append(f"{what}: {type(ex).__name__}: {ex}")
        return default


def _annotations(chain_a: S.Atoms, glycans: Sequence[Glycan], links: Sequence[S.Link], chain: str,
                 partner_chains: Sequence[str], uni_seq: str, acc: str | None, pdb: str, warnings: list[str],
                 construct_range: tuple[int, int], native_range: tuple[int, int], reference_given: bool = True) -> dict[str, Any]:
    """UniProt, SIFTS and reference-structure annotations: topology, domains, segments, glycosylation, disulfides,
    isoforms, unresolved stretches and functional sites. Everything cited with its source; a source that cannot be read
    is left out and noted in ``warnings``. Without an accession (``acc`` None) no UniProt or SIFTS record is read: the
    topology, domains, isoforms and functional sites stay empty and the glycosylation sites and disulfide bonds are the
    ones the reference structure itself shows. ``reference_given`` false (the construct is the frame): no SIFTS record is read, because there is
    no PDB entry to ask about, and nothing is reported as unresolved, since the construct file says nothing about what a full-length structure
    would lack."""
    from . import refdata as rd

    entry = (_safe(lambda: rd.get_uniprot(acc), warnings, f"UniProt {acc}", {}) or {}) if acc else {}
    feats = entry.get("features", [])

    def loc(f):
        return int(f["location"]["start"]["value"]), int(f["location"]["end"]["value"])

    ref_nums = chain_a.residues.number
    ref_last = int(ref_nums.max())
    desc = (entry.get("proteinDescription") or {})
    rec = ((desc.get("recommendedName") or {}).get("fullName") or {}).get("value")
    name = rec or acc or ""
    genes = entry.get("genes") or []
    gene = (genes[0].get("geneName") or {}).get("value") if genes else None

    topology = []
    for f in feats:
        if f["type"] == "Topological domain":
            a, b = loc(f)
            topology.append({"label": f.get("description", ""), "start": a, "end": b, "source": "UniProt topological domain"})
        elif f["type"] == "Transmembrane":
            a, b = loc(f)
            topology.append({"label": "Transmembrane", "start": a, "end": b, "note": f.get("description", ""),
                             "source": "UniProt transmembrane region"})
    topology.sort(key=lambda x: x["start"])

    # domains: CATH from SIFTS for the reference chain (the folded domains), Pfam and InterPro as finer annotation
    domains: list[dict[str, Any]] = []
    # SIFTS maps PDB entries. A predicted model has no entry, so there is nothing
    # to ask for and asking would only leave a warning on a healthy review.
    sifts_possible = bool(acc and reference_given and rd.is_pdb_id(pdb))
    cath = (_safe(lambda: rd.get_sifts(pdb, "cath"), warnings, "SIFTS CATH", {}) or {}) if sifts_possible else {}
    for cid, d in (cath.get("CATH") or {}).items():
        for m in d.get("mappings", []):
            if m["chain_id"] == chain:
                domains.append({"source": "CATH", "id": cid, "name": d.get("homology") or d.get("name", ""),
                                "architecture": d.get("architecture", ""), "domain": m.get("domain"),
                                "start": int(m["start"]["author_residue_number"]), "end": int(m["end"]["author_residue_number"])})
    fine = []
    for kind, key in (("pfam", "Pfam"), ("interpro", "InterPro")) if sifts_possible else ():
        data = _safe(lambda: rd.get_sifts(pdb, kind), warnings, f"SIFTS {key}", {}) or {}
        for fid, d in (data.get(key) or {}).items():
            spans = [m for m in d.get("mappings", []) if m["chain_id"] == chain]
            if spans:
                m = spans[0]
                fine.append({"source": key, "id": fid, "name": d.get("name", ""), "description": d.get("description", d.get("name", "")),
                             "start": int(m["start"]["author_residue_number"]), "end": int(m["end"]["author_residue_number"])})
    domains.sort(key=lambda x: x["start"])
    fine.sort(key=lambda x: (x["source"], x["start"]))

    # segments for the frame strip: UniProt topology up to the first folded domain (the stretch between the membrane
    # anchor and the first domain is the linker), then the CATH domains; a domain runs on to the next one, and the last
    # one to the end of the extracellular region (the residues after the last observed one are not in the structure)
    segments: list[dict[str, Any]] = []
    nat_lo, nat_hi = native_range
    for t in topology:
        if t["label"] == "Transmembrane":
            segments.append({"label": "TM", "start": t["start"], "end": t["end"], "kind": "x"})
        elif "ytoplasmic" in t["label"]:
            segments.append({"label": "Cytoplasmic", "start": t["start"], "end": t["end"], "kind": "x"})
    ecto = next((t for t in topology if t["label"].lower().startswith("extracellular")), None)
    if ecto is not None and domains:
        if domains[0]["start"] > ecto["start"]:
            segments.append({"label": "Linker", "start": ecto["start"], "end": domains[0]["start"] - 1, "kind": "l"})
        for i, d in enumerate(domains):
            label, phrase = d["name"], "the " + d["name"]
            if i + 1 < len(domains):
                end = max(d["end"], domains[i + 1]["start"] - 1)
            else:
                end = ecto["end"] if d["end"] >= ref_last else d["end"]
            segments.append({"label": label, "start": d["start"], "end": end,
                             "kind": DOMAIN_KINDS[min(i, len(DOMAIN_KINDS) - 1)], "phrase": phrase})
    segments.sort(key=lambda s: s["start"])
    if not segments and acc and reference_given:
        warnings.append("no frame segments: UniProt topology or SIFTS domains were not available")

    # N-glycosylation: UniProt sites with the sequon and what the reference structure shows for each chain (without an
    # accession: the sites the reference structure shows on the reference chain, with the sequon read from its residues)
    gly_ref = {(g.owner, g.site): g for g in glycans}
    ref_letters = {int(n): str(a) for n, a in zip(chain_a.residues.number, chain_a.residues.one)}
    glycosylation = []
    sites = ([loc(f)[0] for f in feats if f["type"] == "Glycosylation" and "N-linked" in f.get("description", "")] if acc
             else sorted({g.site for g in glycans if g.owner == chain}))
    for site in sorted(sites):
        sequon = (uni_seq[site - 1:site + 2] if uni_seq else
                  "".join(ref_letters.get(site + i, "") for i in range(3)) if not acc else "")
        ev: dict[str, Any] = {}
        for c in (chain, *partner_chains):
            g = gly_ref.get((c, site))
            ev[c] = None if g is None else {"sugars": list(g.resnames), "src": g.src}
        glycosylation.append({"site": site, "sequon": sequon, "sequon_ok": bool(sequon) and sequon[0] == "N" and len(sequon) == 3
                              and sequon[1] != "P" and sequon[2] in "ST",
                              "in_construct": construct_range[0] <= site <= construct_range[1], "reference": ev,
                              "source": ("UniProt carbohydrate feature (N-linked GlcNAc)" if acc else
                                         "N-glycan attached to this Asn in the reference structure")})
    # disulfides: UniProt bonds and the SG-SG distances the reference file reports for each chain
    ds_links = {}
    for ln in links:
        if ln.kind == "disulf":
            ds_links[(ln.a.chain, tuple(sorted((ln.a.resnum, ln.b.resnum))))] = ln.distance
    disulfides = []
    pairs = ([loc(f) for f in feats if f["type"] == "Disulfide bond"] if acc
             else sorted({pair for (c, pair) in ds_links if c == chain and len(pair) == 2}))
    for a, b in pairs:
        disulfides.append({"pair": [a, b], "in_construct": construct_range[0] <= a and b <= construct_range[1],
                           "reference": {c: ds_links.get((c, (a, b))) for c in (chain, *partner_chains)},
                           "source": ("UniProt disulfide feature; distances from the reference links" if acc else
                                      "disulfide link of the reference structure")})
    disulfides.sort(key=lambda d: d["pair"])

    # unresolved: native range minus the observed reference residues
    observed = _ranges(int(n) for n in ref_nums)
    unresolved = _complement(observed, nat_lo, nat_hi) if reference_given else []

    isoforms = []
    for c in entry.get("comments", []):
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in c.get("isoforms", []):
                isoforms.append({"id": (iso.get("isoformIds") or [None])[0], "name": (iso.get("name") or {}).get("value"),
                                 "status": iso.get("isoformSequenceStatus"), "synonyms": [s.get("value") for s in iso.get("synonyms", [])]})
    for f in feats:
        if f["type"] == "Alternative sequence" and "isoform" in f.get("description", ""):
            a, b = loc(f)
            for iso in isoforms:
                if f["description"].endswith("isoform " + str(iso["name"])):
                    iso.setdefault("removed", []).append([a, b])
    for iso in isoforms:
        rem = _expand(iso.get("removed", []))
        iso["contains_construct"] = not (rem & set(range(construct_range[0], construct_range[1] + 1)))
    chains_u = [{"label": f.get("description", ""), "start": loc(f)[0], "end": loc(f)[1], "id": f.get("featureId")}
                for f in feats if f["type"] == "Chain"]
    signal = [{"start": loc(f)[0], "end": loc(f)[1]} for f in feats if f["type"] == "Signal"]
    functional = []
    for f in feats:
        if f["type"] in ("Active site", "Binding site", "Site"):
            a, b = loc(f)
            what = f.get("description") or (f.get("ligand") or {}).get("name", "")
            functional.append({"type": f["type"], "start": a, "end": b, "description": what,
                               "residues": ([a] if a == b else list(range(a, b + 1))),
                               "in_construct": construct_range[0] <= a and b <= construct_range[1]})
    signal_note = ""
    if topology and not signal:
        signal_note = "no cleaved signal peptide: the transmembrane helix acts as a signal-anchor (type II membrane protein)"
    return {"name": name, "recommended_name": rec, "gene": gene, "organism": (entry.get("organism") or {}).get("scientificName"),
            "topology": topology, "domains": domains, "domains_fine": fine, "segments": segments, "glycosylation": glycosylation,
            "disulfides": disulfides, "signal": signal, "signal_note": signal_note, "chains": chains_u, "isoforms": isoforms,
            "unresolved": unresolved, "observed": observed, "functional_sites": functional,
            "sources": ({"uniprot": f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                         **({"sifts": f"PDBe SIFTS mappings for {pdb.lower()}"} if reference_given else {})} if acc else {})}


def _dimer_report(path: str, chain: str, partners: Sequence[str], A: S.Atoms, P: dict[str, S.Atoms],
                  glycans: Sequence[Glycan], probe: float, points: int, olig_min_area: float,
                  assembly_id: str = "1", subunit_first: bool = True) -> dict[str, Any]:
    """Assembly records, buried surface, chain A interface residues, the no-overlap check and the two-fold axis. The two-fold axis
    relates the target chain to the first partner, so it is only read when that partner is the same protein (``subunit_first``): the
    residue numbers of another molecule have nothing to do with the target's."""
    asm = _parse_assembly(path)
    out: dict[str, Any] = {"assemblies": asm["assemblies"], "generators": asm["generators"]}
    a1 = next((a for a in asm["assemblies"] if a["id"] == assembly_id), {})
    gen = [g for g in asm["generators"] if g["assembly_id"] == a1.get("id")]
    labels = [asm["label_asym"].get(c) for c in (chain, *partners)]
    ops = [asm["operators"].get(e.strip(), "") for g in gen for e in g["oper_expression"].split(",")]
    out.update(oligomeric_details=a1.get("oligomeric_details", ""), oligomeric_count=a1.get("oligomeric_count", ""),
               method=a1.get("method", ""), identity_only=bool(ops) and all("identity" in o.lower() for o in ops),
               chains_in_assembly=bool(gen) and all(l in gen[0]["asym_ids"] for l in labels if l))
    if not partners:
        return out
    B = S.Atoms.concat([P[c] for c in partners])
    sA = S.sasa_atoms(A.xyz, A.element, probe, points)
    sB = S.sasa_atoms(B.xyz, B.element, probe, points)
    cat = S.Atoms.concat([A, B])
    sAB = S.sasa_atoms(cat.xyz, cat.element, probe, points)
    dA, dB = sA - sAB[:len(A)], sB - sAB[len(A):]
    res_d = A.per_residue(dA)
    out["bsa_chain_a"] = round(float(dA.sum()), 1)
    out["bsa_partner"] = round(float(dB.sum()), 1)
    out["bsa_total"] = round(float(dA.sum() + dB.sum()), 1)
    out["interface_residues_any"] = [int(n) for n, v in zip(A.residues.number, res_d) if v > 0.5]
    out["interface_residues"] = [int(n) for n, v in zip(A.residues.number, res_d) if v >= olig_min_area]
    out["n_interface"] = len(out["interface_residues"])
    # overlap
    tb = cKDTree(B.xyz)
    d, i = cKDTree(A.xyz).query(B.xyz)
    j = int(np.argmin(d))
    pairs = cKDTree(A.xyz).sparse_distance_matrix(tb, 3.0, output_type="coo_matrix")
    out["min_distance"] = round(float(d[j]), 3)
    out["closest_pair"] = {"partner": [str(B.chain[j]), int(B.resnum[j]), str(B.resname[j]), str(B.name[j])],
                           "chain": [chain, int(A.resnum[i[j]]), str(A.resname[i[j]]), str(A.name[i[j]])]}
    pd = np.sort(pairs.data)
    out["pairs_below"] = {"2.2": int((pd < 2.2).sum()), "2.5": int((pd < 2.5).sum()), "3.0": int((pd < 3.0).sum())}
    heavy_c = [(a, b) for a, b in zip(A.element[pairs.row], B.element[pairs.col])]
    out["pairs_below_3_carbon_involved"] = int(sum(1 for a, b in heavy_c if "C" in (a, b)))
    allB = S.Atoms.concat([B, _glycan_atoms([g for g in glycans if g.owner in partners])])
    allA = S.Atoms.concat([A, _glycan_atoms([g for g in glycans if g.owner == chain])])
    out["min_distance_with_sugars"] = round(float(cKDTree(allB.xyz).query(allA.xyz)[0].min()), 3)
    # two-fold axis relating chain A to the partner (needs the same residues in both)
    pa = P[partners[0]]
    keep_p, shared_p, _ = _first_residue_per_number(pa)
    if shared_p:                                                          # the same residue of a shared number on both sides of the pair
        pa = pa.select(keep_p)
    ca_a, ca_b = A.ca_by_number(), pa.ca_by_number()
    letter_a = {int(n): str(o) for n, o in zip(A.residues.number, A.residues.one)}
    letter_b = {int(n): str(o) for n, o in zip(pa.residues.number, pa.residues.one)}
    common = sorted(set(ca_a) & set(ca_b))
    agree = sum(1 for n in common if letter_a.get(n) == letter_b.get(n)) / len(common) if common else 0.0
    Pm = Qm = None
    if subunit_first and len(common) > 10 and agree >= 0.9:
        Pm = np.array([ca_a[n] for n in common])
        Qm = np.array([ca_b[n] for n in common])
    elif subunit_first:
        # the numbers of the other subunit do not pair the same residues (the second chain of many entries starts at 301): pair them by sequence
        ia, ib = _matched_pairs("".join(letter_a.values()), "".join(letter_b.values()))
        cA, cB = A.ca_xyz(), pa.ca_xyz()
        ok = [(i, j) for i, j in zip(ia, ib) if not np.isnan(cA[i, 0]) and not np.isnan(cB[j, 0])]
        if len(ok) > 10:
            Pm, Qm = cA[[i for i, _ in ok]], cB[[j for _, j in ok]]
            common = [i for i, _ in ok]
    if Pm is not None and Qm is not None:
        R, t, rms = S.kabsch(Pm, Qm)
        if rms > _TWO_FOLD_MAX_RMSD:
            # an operator that does not superpose the two subunits is no operator: the twin check is left out, not read from it
            out["two_fold_skipped"] = (f"chain {chain} and its first partner superpose at {rms:.1f} A over {len(common)} residues, "
                                       f"more than {_TWO_FOLD_MAX_RMSD:g} A, so they are not copies of one another")
            return out
        angle = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
        w, v = np.linalg.eig(R)
        axis = _unit(np.real(v[:, int(np.argmin(abs(w - 1)))]))
        tperp = t - (t @ axis) * axis
        point = np.linalg.lstsq(np.eye(3) - R, tperp, rcond=None)[0]
        out["two_fold"] = {"rmsd": round(rms, 3), "n": len(common), "angle": round(angle, 2), "axis": _round(axis, 4),
                           "point": _round(point, 2), "screw": round(float(t @ axis), 3),
                           "R": [round(float(x), 6) for x in R.ravel()], "t": [round(float(x), 3) for x in t]}
    return out


def auto_partner_chains(sf: S.StructureFile, chain: str, min_shared: int = 20, min_identity: float = 0.9) -> list[str]:
    """The other polymer chains of the reference that carry the same protein as ``chain`` (the other subunits of an
    oligomer): at least ``min_shared`` residue numbers in common, ``min_identity`` of them with the same residue type."""
    def letters(c: str) -> dict[int, str]:
        r = sf.groups[c].residues
        return {int(n): str(a) for n, a in zip(r.number, r.one)}

    ref = letters(chain)
    out = []
    for c in sf.polymer_chains():
        if c == chain:
            continue
        other = letters(c)
        shared = [n for n in other if n in ref]
        if len(shared) >= min_shared and sum(1 for n in shared if ref[n] == other[n]) >= min_identity * len(shared):
            out.append(c)
    return out


_EXPERIMENTAL_METHODS = ("X-RAY", "ELECTRON", "NMR", "NEUTRON", "FIBER", "SOLUTION SCATTERING", "EPR", "INFRARED", "FLUORESCENCE", "POWDER")


def classify_reference(path: str | Path | None) -> dict[str, Any]:
    """Is the reference an experimental structure or a prediction? Decided here, from the file, so it holds for a command-line user with no
    web service in sight.

    The distinction matters because the context checks measure the binder against the reference as one static model: against an experimental
    structure that is a measurement, against a prediction it is only as good as the prediction (AlphaFold places a region that has no fixed
    position, such as a flexible linker, essentially at random). ``kind`` is ``experimental``, ``predicted``, ``unknown`` or ``none`` (no
    reference); ``treated_as`` is what the pages act on, and ``unknown`` is treated as predicted, the safer reading. ``evidence`` lists what
    was found, in words, so a reader can disagree with the call.

    Experimental needs a positive record (the deposited method, ``_exptl.method`` or ``EXPDTA``) and no prediction record. Prediction needs
    a model record (the ``_ma_`` categories of a model-archive file, which AlphaFold DB writes) or a theoretical-model method; per-residue
    constant B-factors between 0 and 100 (how a pLDDT is stored) count only when the file names no method at all."""
    if path is None:
        return {"kind": "none", "treated_as": "none", "method": None, "resolution": None, "evidence": ["no reference structure was given"]}
    text_path = str(path)
    method: str | None = None
    resolution: float | None = None
    pred: list[str] = []
    gemmi = S._gemmi()
    try:
        low = text_path.lower()
        if low.endswith((".cif", ".mmcif", ".cif.gz")):
            doc = gemmi.cif.read(text_path)
            b = doc[0]
            m = b.find_value("_exptl.method")
            method = m.strip("'\" ") if m and m not in ("?", ".") else None
            for tag in ("_refine.ls_d_res_high", "_reflns.d_resolution_high", "_em_3d_reconstruction.resolution"):
                v = b.find_value(tag)
                if v and v not in ("?", "."):
                    try:
                        resolution = float(v)
                        break
                    except ValueError:
                        pass
            mt = b.find_value("_ma_model_list.model_type")
            if mt:
                pred.append(f"model record _ma_model_list.model_type = {mt.strip(chr(34) + chr(39))}")
            names = b.find_values("_ma_qa_metric.name")
            if names and any("plddt" in str(n).lower() for n in names):
                pred.append("model-quality record _ma_qa_metric: pLDDT")
            grp = b.find_value("_ma_model_list.model_group_name")
            if grp and not mt:
                pred.append(f"model record _ma_model_list.model_group_name = {grp.strip(chr(34) + chr(39))}")
        else:
            with open(text_path, "r", errors="replace") as fh:
                for line in fh:
                    if line.startswith("EXPDTA"):
                        method = line[6:].strip().strip(";") or None
                    elif line.startswith("REMARK   2 RESOLUTION."):
                        try:
                            resolution = float(line.split("RESOLUTION.")[1].split()[0])
                        except (IndexError, ValueError):
                            pass
                    elif line.startswith(("ATOM", "HETATM")):
                        break
    except (OSError, ValueError, RuntimeError):
        return {"kind": "unknown", "treated_as": "predicted", "method": None, "resolution": None,
                "evidence": ["the file could not be read for a method or a model record"]}
    ev: list[str] = []
    if method:
        ev.append(f"deposited method: {method}" + (f", resolution {resolution:g} A" if resolution else ""))
    ev += pred
    theoretical = bool(method) and "THEORETICAL" in method.upper()
    experimental_method = bool(method) and not theoretical and any(k in method.upper() for k in _EXPERIMENTAL_METHODS)
    if pred or theoretical:
        if experimental_method:
            ev.append("the file records both a method and a model: read as a prediction")
        return {"kind": "predicted", "treated_as": "predicted", "method": method, "resolution": resolution, "evidence": ev}
    if experimental_method:
        return {"kind": "experimental", "treated_as": "experimental", "method": method, "resolution": resolution, "evidence": ev}
    # no method and no model record: look at the B-factor column, which a prediction fills with pLDDT (one value per residue, 0 to 100)
    try:
        st = S._read_gemmi(text_path)
        per_res_const, total, bmax = 0, 0, 0.0
        for ch in st[0]:
            for r in ch:
                bs = {round(a.b_iso, 2) for a in r}
                total += 1
                per_res_const += len(bs) == 1
                bmax = max(bmax, max(bs))
        if total >= 20 and per_res_const >= 0.98 * total and bmax <= 100.0:
            ev.append("no method is recorded, and every residue has one B-factor between 0 and 100, which is how a pLDDT is stored")
            return {"kind": "predicted", "treated_as": "predicted", "method": None, "resolution": None, "evidence": ev}
    except Exception:  # noqa: BLE001 - the B-factor look is a last resort; failing it leaves "unknown"
        pass
    ev.append("no method and no model record in the file: cannot tell, treated as predicted")
    return {"kind": "unknown", "treated_as": "predicted", "method": None, "resolution": None, "evidence": ev}


def reference_label(reference_id: str) -> str:
    """The id of the reference structure as it is written in labels: a PDB id is upper case, any other id is kept as given."""
    rid = str(reference_id)
    return rid.upper() if len(rid) == 4 and rid[0].isdigit() and rid.isalnum() else rid


# ----------------------------------------------------------------------------------------------------------------
# what the reference file holds: its chains, its assembly, and what the frame does with each chain
# ----------------------------------------------------------------------------------------------------------------
#: The operator that superposes the target chain on its first partner is only read when it does so within this many angstrom (CA RMSD).
_TWO_FOLD_MAX_RMSD = 5.0
#: When chains carry the construct equally, the target is the first in the file; if the construct's coordinates follow another more closely by this much (A, CA RMSD) the reason says so.
_COORD_CLEAR = 0.25
#: Chains that carry the target's own sequence (at least this share of the shorter chain, and at least ``_MIN_SAME`` residues).
_SAME_IDENTITY = 0.9
_MIN_SAME = 20
#: A copy of the assembly is only kept when it lies within this many angstrom of the target chain, if the assembly is large.
_COPY_REACH = 30.0
_MAX_COPIES = 40
#: Contact radius between a chain and the target chain (heavy atoms), as for the contacts everywhere else.
_TOUCH = 4.5
_TOUCH_LIGAND = 4.0

# Residue names that are, in nearly every entry, something that was put into the crystal or the buffer and not a part of the
# biology. The name is a hint shown to the reader; nothing is decided by it, because none of these is ever measured.
_ADDITIVES = frozenset({"GOL", "EDO", "PEG", "PGE", "PG4", "PE4", "1PE", "2PE", "P6G", "PGO", "MPD", "DMS", "SO4", "PO4", "ACT",
                        "ACY", "FMT", "IMD", "TRS", "MES", "EPE", "BME", "CIT", "TAR", "TLA", "NO3", "SCN", "EOH", "IPA", "BU3",
                        "CAC", "HED", "DTT", "MLI", "MLA", "OXL", "FLC", "BTB", "BCT", "NH4", "PEG", "SPD", "SPM", "UNX", "UNL"})
_SUGAR_NAMES = frozenset({"NAG", "NDG", "BMA", "MAN", "FUC", "FUL", "GAL", "GLC", "BGC", "SIA", "NGA", "A2G", "XYP", "GLA", "FRU", "RIB"})
_NUCLEIC = frozenset({"DA", "DC", "DG", "DT", "DU", "DI", "A", "C", "G", "U", "I", "N"})
# A hint from the name of an entity that a chain is an antibody-like binder or a stabilising partner, which crystallographers
# add to a target in order to get it to crystallise and which is not a part of the target's biology.
_AID_WORDS = ("antibody", "nanobody", "sybody", "megabody", "monobody", "affibody", "darpin", "immunoglobulin", " fab ", "fab fragment",
              "heavy chain", "light chain", "scfv", "vhh", "single-chain variable", " fv ", "chaperone")


def _wlist(chains: Sequence[str]) -> str:
    """'chain A', 'chains A and B', 'chains A, B and C'."""
    names = [str(c) for c in chains]
    if len(names) == 1:
        return f"chain {names[0]}"
    return f"chains {', '.join(names[:-1])} and {names[-1]}"


def _moved(at: S.Atoms, R: np.ndarray, t: np.ndarray, name: str) -> S.Atoms:
    """The atoms placed by ``x' = x @ R.T + t`` and given the chain name ``name``."""
    return S.Atoms(at.xyz @ np.asarray(R, dtype=float).T + np.asarray(t, dtype=float), at.name, at.element, at.resname, at.resnum,
                   np.full(len(at), name), at.bfactor, at.icode, at.hetero)


def _moved_glycan(g: "Glycan", R: np.ndarray, t: np.ndarray, name: str) -> "Glycan":
    return Glycan(name, g.site, list(g.resnames), g.src, g.link_distance, g.anchor @ np.asarray(R, dtype=float).T + np.asarray(t, dtype=float),
                  _moved(g.atoms, R, t, name))


def _matched_positions(a: str, b: str) -> tuple[set[int], set[int]]:
    """The positions of ``a`` and of ``b`` that lie in a common stretch of the two sequences (difflib blocks, junk heuristic off,
    so a 150-residue sequence is not quietly compared on its rarer letters only)."""
    import difflib

    ia: set[int] = set()
    ib: set[int] = set()
    if a and b:
        for blk in difflib.SequenceMatcher(None, a, b, autojunk=False).get_matching_blocks():
            ia.update(range(blk.a, blk.a + blk.size))
            ib.update(range(blk.b, blk.b + blk.size))
    return ia, ib


def _matched_pairs(a: str, b: str) -> tuple[list[int], list[int]]:
    """Paired positions of two sequences that lie in a common stretch, in order: ``(positions in a, positions in b)``."""
    import difflib

    ia: list[int] = []
    ib: list[int] = []
    if a and b:
        for blk in difflib.SequenceMatcher(None, a, b, autojunk=False).get_matching_blocks():
            ia.extend(range(blk.a, blk.a + blk.size))
            ib.extend(range(blk.b, blk.b + blk.size))
    return ia, ib


def _entity_names(path: str) -> tuple[dict[str, str], dict[str, str]]:
    """What the file calls its molecules: ({entity id: description}, {chain: description}). An mmCIF has the first; a PDB file's COMPND records
    give the second. Either may be empty, and a name is only ever shown, never used to decide anything."""
    gemmi = S._gemmi()
    by_entity: dict[str, str] = {}
    by_chain: dict[str, str] = {}
    try:
        block = gemmi.cif.read(str(path)).sole_block()
        for row in block.find("_entity.", ["id", "pdbx_description"]):
            d = gemmi.cif.as_string(row[1]).strip()
            if d and d not in ("?", "."):
                by_entity[gemmi.cif.as_string(row[0])] = d
        return by_entity, by_chain
    except (ValueError, RuntimeError, OSError):
        pass
    try:
        text = Path(path).read_text(errors="ignore")
    except OSError:
        return by_entity, by_chain
    mol = ""
    pending = ""
    for line in text.splitlines():
        if not line.startswith("COMPND"):
            continue
        body = line[10:].strip()
        pending = (pending + " " + body).strip() if pending else body
        if not pending.endswith(";"):
            continue
        key, _, value = pending.rstrip(";").partition(":")
        pending = ""
        key, value = key.strip().upper(), value.strip()
        if key == "MOLECULE":
            mol = value
        elif key == "CHAIN" and mol:
            for c in value.replace(" ", "").split(","):
                if c:
                    by_chain[c] = mol
    return by_entity, by_chain


def _tidy_name(text: str) -> str:
    """Older entries write every description in capitals."""
    text = " ".join(str(text).split())
    return text.capitalize() if text.isupper() and len(text) > 3 else text


class ReferenceChains:
    """What a reference file holds and what the frame does with each chain of it.

    ``data`` is plain JSON (the table the service shows and ``frame.meta["composition"]`` keeps): every chain with the role it was given
    and the reason, the target and why it is the target, the partners and why, the assembly, and the non-polymer molecules the file
    holds. ``target`` is the chosen chain, ``partners`` the chain names the frame will measure the binder against (chains of the file and
    copies that the assembly records make from them), and ``problem`` is a sentence when no target could be chosen. :meth:`atoms` and
    :meth:`glycans` give the atoms of a partner, whether it is a chain of the file or a generated copy."""

    def __init__(self, data: dict[str, Any], sf: S.StructureFile, copies: dict[str, tuple[str, np.ndarray, np.ndarray]],
                 glycans: Sequence["Glycan"]):
        self.data = data
        self.target: str | None = data["target"]["chain"]
        self.partners: list[str] = list(data["partners"]["chains"])
        self.problem: str | None = data.get("problem")
        self._sf = sf
        self._copies = copies
        self._glycans = list(glycans)

    def atoms(self, name: str) -> S.Atoms:
        """Heavy atoms of a chain of the file, or of the copy that the assembly makes from one."""
        if name in self._copies:
            base, R, t = self._copies[name]
            return _moved(self._sf.groups[base].heavy(), R, t, name)
        return self._sf.groups[name].heavy()

    def glycans(self, name: str) -> list["Glycan"]:
        """The N-glycans of a chain of the file, or of a generated copy (moved with it)."""
        if name in self._copies:
            base, R, t = self._copies[name]
            return [_moved_glycan(g, R, t, name) for g in self._glycans if g.owner == base]
        return [g for g in self._glycans if g.owner == name]

    def role(self, name: str) -> str:
        for row in self.data["chains"]:
            if row["chain"] == name:
                return row["role"]
        return ""


def _parse_partner_request(partner_chains: Any) -> tuple[str, list[str]]:
    """``(mode, names)``: ``auto`` (None, "auto", ""), ``none`` ((), [], "none"), ``all`` ("all") or ``given`` (a list, or a comma list)."""
    if partner_chains is None:
        return "auto", []
    if isinstance(partner_chains, str):
        text = partner_chains.strip()
        low = text.lower()
        if low in ("", "auto"):
            return "auto", []
        if low == "none":
            return "none", []
        if low == "all":
            return "all", []
        names = [c.strip() for c in text.split(",") if c.strip()]
    else:
        names = [str(c).strip() for c in partner_chains if str(c).strip()]
    return ("given", names) if names else ("none", [])


def analyse_reference(reference: str | Path, construct_sequence: str | None = None, *, chain: str | None = None,
                      partner_chains: Sequence[str] | str | None = None, assembly_id: str | None = None,
                      sf: S.StructureFile | None = None, glycans: Sequence["Glycan"] | None = None,
                      construct_atoms: S.Atoms | None = None) -> ReferenceChains:
    """Read what a reference structure holds and decide, explainably, which chain is the target and which are its partners.

    ``construct_sequence``: the one-letter sequence of the construct (chain A of the construct file); with it the target chain is
    found by sequence when ``chain`` is None, and every protein chain is reported with how much of the construct it holds.
    ``construct_atoms``: the construct's atoms (chain A) instead of its sequence; when several chains carry the construct equally
    (the copies of a homo-oligomer) the first in the file is the target, and the reason says whether the construct's coordinates
    follow another chain more closely, which is for the person to act on.
    ``chain``: the chain of the target, or None to let the sequence decide. ``partner_chains``: None or ``"auto"`` (the default: the other
    subunits, that is chains that carry the same protein as the target and belong to its assembly, plus the copies the assembly makes
    from them by symmetry), ``"none"`` or ``()``, ``"all"`` (every other polymer chain of the assembly) or the names of chains (also the
    generated copies' names, which :func:`describe_reference` lists). ``assembly_id``: the assembly of the file to take the subunits from
    (default: the first one that holds the target chain). Never raises for a reference that holds no usable target: ``problem`` says why.
    Nothing here guesses a role the file cannot support: chains of another molecule are never counted unless asked for, and ligands, ions,
    sugars that are not N-glycans and waters are listed but are never partners (docs/target_frame.md, "Chains of the reference")."""
    gemmi = S._gemmi()
    path = str(reference)
    sf = sf or S.read_structure(path)
    st = S._read_gemmi(path)
    pdb = reference_label(Path(path).stem)
    mode, asked = _parse_partner_request(partner_chains)
    notes: list[str] = []

    # ---- the polymer chains of the file
    by_entity, by_chain = _entity_names(path)
    entity_of: dict[str, str] = {}
    for ch in st[0]:
        for res in ch:
            if res.entity_id and ch.name in sf.groups and sf.kinds.get(ch.name) == "polymer" and ch.name not in entity_of:
                ent = st.get_entity(res.entity_id)
                if ent is not None and ent.entity_type == gemmi.EntityType.Polymer:
                    entity_of[ch.name] = res.entity_id
                    break
    info: dict[str, dict[str, Any]] = {}
    for c in sf.polymer_chains():
        at = sf.groups[c]
        r = at.residues
        names = [str(n) for n in r.name]
        n_aa = sum(1 for n in names if n in S.THREE_TO_ONE or n in S._MODIFIED)
        n_nuc = sum(1 for n in names if n in _NUCLEIC)
        kind = "protein" if n_aa >= max(1, 0.5 * len(names)) else "nucleic acid" if n_nuc >= max(1, 0.5 * len(names)) else "other polymer"
        numbers = [int(x) for x in r.number]
        info[c] = {"kind": kind, "seq": "".join(str(x) for x in r.one) if kind == "protein" else "", "n": len(names),
                   "first": numbers[0] if numbers else None, "last": numbers[-1] if numbers else None,
                   "entity": entity_of.get(c, ""),
                   "name": _tidy_name(by_entity.get(entity_of.get(c, ""), "") or by_chain.get(c, ""))}
    protein = [c for c, v in info.items() if v["kind"] == "protein"]
    c_ca: np.ndarray | None = None
    if construct_atoms is not None and len(construct_atoms):
        construct_sequence = "".join(str(a) for a in construct_atoms.residues.one)
        c_ca = construct_atoms.ca_xyz()
    cseq = "".join(ch for ch in (construct_sequence or "") if ch.isalpha()).upper()

    # ---- which chain is the target, and why
    problem: str | None = None
    name_it = False                      # the problem is solved by naming the chain
    how, why, warning = "", "", ""
    alternatives: list[str] = []
    construct_hold: dict[str, int] = {}
    if cseq:
        for c in protein:
            construct_hold[c] = len(_matched_positions(cseq, info[c]["seq"])[0])
    construct_rmsd: dict[str, float] = {}
    if c_ca is not None and cseq:
        for c in protein:
            if construct_hold.get(c, 0) < _SAME_IDENTITY * len(cseq):
                continue
            ia, ib = _matched_pairs(cseq, info[c]["seq"])
            cc = sf.groups[c].ca_xyz()
            ok = [(i, j) for i, j in zip(ia, ib) if not np.isnan(c_ca[i, 0]) and not np.isnan(cc[j, 0])]
            if len(ok) >= 10:
                construct_rmsd[c] = float(S.kabsch(c_ca[[i for i, _ in ok]], cc[[j for _, j in ok]])[2])
    target: str | None = None
    if chain:
        if chain not in info:
            problem, name_it = f"{Path(path).name}: chain {chain!r} is not a polymer chain of the file (its polymer chains: {', '.join(info) or 'none'}).", True
        elif info[chain]["kind"] != "protein":
            problem, name_it = f"chain {chain} of {Path(path).name} is {info[chain]['kind']}, not a protein: the target must be a protein chain.", True
        else:
            target, how = chain, "given"
            why = "as asked" + (f"; it holds {construct_hold[chain]} of the construct's {len(cseq)} residues" if cseq else "")
            if cseq and construct_hold[chain] < _SAME_IDENTITY * len(cseq):
                better = max(protein, key=lambda c: construct_hold[c])
                warning = (f"chain {chain} holds only {construct_hold[chain]} of the construct's {len(cseq)} residues"
                           + (f"; chain {better} holds {construct_hold[better]}" if construct_hold[better] > construct_hold[chain] else ""))
    elif not protein:
        problem = f"{Path(path).name} holds no protein chain, so there is no target in it."
    elif not cseq:
        if len(protein) == 1:
            target, how, why = protein[0], "only_chain", "the only protein chain of the file"
        else:
            problem, name_it = (f"{Path(path).name} has several protein chains ({', '.join(protein)}) and no construct was given to choose "
                                "between them: name the chain of the target."), True
    else:
        need = _SAME_IDENTITY * len(cseq)
        fits = sorted((c for c in protein if construct_hold[c] >= need), key=lambda c: (-construct_hold[c], protein.index(c)))
        if fits:
            target = fits[0]
            top = construct_hold[target]
            ties = [c for c in fits if construct_hold[c] == top]
            alternatives = [c for c in fits if c != target]
            if not alternatives:
                how, why = "only_match", f"the only chain that carries the construct ({top} of {len(cseq)} residues)"
            elif len(ties) > 1:
                how = "first_of_equal"
                fit = ""
                if all(c in construct_rmsd for c in ties):
                    ranked = sorted(ties, key=lambda c: construct_rmsd[c])
                    if construct_rmsd[ranked[1]] - construct_rmsd[ranked[0]] >= _COORD_CLEAR and ranked[0] != target:
                        fit = (f"; the construct's own coordinates follow chain {ranked[0]} more closely (CA RMSD {construct_rmsd[ranked[0]]:.2f} A, "
                               f"against {construct_rmsd[target]:.2f} A for {target}), so choose {ranked[0]} if the construct was cut from it")
                    else:
                        fit = " and follow the construct's coordinates about equally closely (CA RMSD " + ", ".join(f"{construct_rmsd[c]:.2f}" for c in ties) + " A)"
                why = (f"{_wlist(ties)} carry the construct equally ({top} of {len(cseq)} residues each){fit}; {target} was taken because it "
                       "comes first in the file. Choose another chain to use it instead")
            else:
                how = "best_match"
                held = ", ".join(f"{c} {construct_hold[c]}" for c in fits)
                why = f"it carries the most of the construct ({held} of {len(cseq)} residues)"
        else:
            best = max(protein, key=lambda c: construct_hold[c])
            union = set()
            for c in protein:
                union |= _matched_positions(cseq, info[c]["seq"])[0]
            split = [c for c in protein if construct_hold[c] >= 0.1 * len(cseq)]
            if len(split) > 1 and len(union) >= need:
                problem = (f"the construct is split over {_wlist(split)} of {Path(path).name} "
                           f"({', '.join(f'{c} {construct_hold[c]}' for c in split)} of {len(cseq)} residues) and no single chain holds it; "
                           "the review reads the target from one chain, so a reference that holds the construct in one chain is needed.")
            else:
                problem = (f"no chain of {Path(path).name} carries the construct: the closest, chain {best}, holds "
                           f"{construct_hold[best]} of its {len(cseq)} residues. Is this the right protein?")

    # ---- the assembly the target belongs to
    asms = list(st.assemblies)
    sub_of: dict[str, str] = {}
    for ch in st[0]:
        poly = ch.get_polymer()
        if len(poly) and ch.name in info:
            sub_of[ch.name] = poly.subchain_id()
    chosen = None
    if asms and target:
        def holds(a, c: str) -> bool:
            return any(c in list(g.chains) or sub_of.get(c) in list(g.subchains) for g in a.generators)
        if assembly_id:
            chosen = next((a for a in asms if a.name == str(assembly_id)), None)
            if chosen is None:
                problem = problem or f"{Path(path).name} has no assembly {assembly_id} (it has {', '.join(a.name for a in asms)})."
        else:
            chosen = next((a for a in asms if holds(a, target)), None)
        if chosen is not None and not holds(chosen, target):
            notes.append(f"chain {target} is not part of assembly {chosen.name} of the file")
    elif assembly_id and not asms:
        notes.append(f"the file has no assembly records, so assembly {assembly_id} could not be used")

    # ---- the copies the assembly makes (instances beyond the chains of the file) and the table of chains
    first_of = {c: i for i, c in enumerate(info)}
    in_assembly: dict[str, bool | None] = {c: None for c in info}
    generated: list[tuple[str, str, str, np.ndarray, np.ndarray]] = []       # (name, base chain, operator name, R, t)
    asm_info: dict[str, Any] = {"records": bool(asms), "id": chosen.name if chosen is not None else None,
                                "oligomeric": (chosen.oligomeric_details if chosen is not None else "") or "",
                                "available": [a.name for a in asms], "generated": []}
    expanded_ok = False
    if chosen is not None and target:
        members: set[str] = set()
        instances: list[tuple[str, Any]] = []                                 # (chain, operator)
        for g in chosen.generators:
            cs = [c for c in info if c in list(g.chains) or sub_of.get(c) in list(g.subchains)]
            members.update(cs)
            for op in g.operators:
                instances.extend((c, op) for c in cs)
        for c in info:
            in_assembly[c] = c in members
        identity_target = any(c == target and op.transform.is_identity() for c, op in instances)
        copies = [(c, op) for c, op in instances if not op.transform.is_identity()]
        if copies and not identity_target:
            notes.append(f"assembly {chosen.name} does not place chain {target} where the file has it, so its symmetry copies were not generated")
        elif copies:
            expanded_ok = True
            tgt_ca = sf.groups[target].heavy()
            tgt_tree = cKDTree(tgt_ca.xyz)
            kept: list[tuple[str, Any, np.ndarray, np.ndarray]] = []
            dropped = 0
            for c, op in copies:
                R = np.array(op.transform.mat.tolist(), dtype=float)
                t = np.array(op.transform.vec.tolist(), dtype=float)
                if len(copies) > _MAX_COPIES:
                    ca = sf.groups[c].heavy().xyz[::4] @ R.T + t
                    if float(tgt_tree.query(ca)[0].min()) > _COPY_REACH:
                        dropped += 1
                        continue
                kept.append((c, op, R, t))
            taken = set(info)
            for c, op, R, t in kept:
                nm = f"{c}-{op.name}"
                k = 1
                while len(nm) > 4 or nm in taken:
                    nm = f"{c[:1]}~{k}"
                    k += 1
                taken.add(nm)
                generated.append((nm, c, str(op.name), R, t))
            if dropped:
                asm_info["left_out"] = dropped
                notes.append(f"{dropped} further copies made by assembly {chosen.name} lie more than {_COPY_REACH:g} A from the target chain and were left out")
    asm_info["generated"] = [g[0] for g in generated]
    asm_info["expanded"] = expanded_ok

    # ---- every chain and copy as a row: what it is, how it relates to the target, what the frame does with it
    cinfo = dict(info)
    inst_atoms: dict[str, S.Atoms] = {}
    for c in info:
        inst_atoms[c] = sf.groups[c].heavy()
    for nm, base, opn, R, t in generated:
        inst_atoms[nm] = _moved(sf.groups[base].heavy(), R, t, nm)
        cinfo[nm] = {**info[base], "base": base, "operator": opn}
    same: dict[str, bool] = {}
    if target:
        t_info = info[target]
        for nm, v in info.items():
            if nm == target:
                same[nm] = True
            elif v["kind"] != "protein":
                same[nm] = False
            elif t_info["entity"] and v["entity"] and t_info["entity"] == v["entity"]:
                same[nm] = True
            else:
                short = min(len(v["seq"]), len(t_info["seq"]))
                ident = len(_matched_positions(v["seq"], t_info["seq"])[0]) if short else 0
                same[nm] = short >= _MIN_SAME and ident >= _SAME_IDENTITY * short
        for nm, v in cinfo.items():
            if nm not in same:
                same[nm] = same[v["base"]]
    tgt_tree = cKDTree(inst_atoms[target].xyz) if target else None
    tgt_resnum = inst_atoms[target].resnum if target else None
    rows: list[dict[str, Any]] = []
    for nm, v in cinfo.items():
        base = v.get("base", nm)
        row: dict[str, Any] = {"chain": nm, "kind": v["kind"], "entity": v["entity"], "name": v["name"], "n_residues": v["n"],
                               "first": v["first"], "last": v["last"], "generated": nm != base, "role": "", "used": False,
                               "why": "", "touches": None, "closest": None, "hint": ""}
        if nm != base:
            row["copy_of"], row["operator"] = base, v["operator"]
        if cseq and v["kind"] == "protein":
            row["construct"] = {"found": construct_hold.get(base, 0), "of": len(cseq)}
            if base in construct_rmsd:
                row["construct"]["rmsd"] = round(construct_rmsd[base], 2)
        if tgt_tree is not None and nm != target:
            d, ix = tgt_tree.query(inst_atoms[nm].xyz, distance_upper_bound=_TOUCH)
            near = np.isfinite(d)
            row["touches"] = int(len(set(int(x) for x in tgt_resnum[ix[near]])))
            row["closest"] = round(float(tgt_tree.query(inst_atoms[nm].xyz)[0].min()), 1)
        low = (v["name"] or "").lower()
        if v["kind"] == "protein" and any(w in f" {low} " for w in _AID_WORDS):
            row["hint"] = "named like an antibody-like binder or a crystallisation aid"
        rows.append(row)
    # roles
    for row in rows:
        nm = row["chain"]
        base = row.get("copy_of", nm)
        if target is None:
            row["role"] = "unassigned"                       # no target, so no chain has a relation to it
        elif nm == target:
            row["role"] = "target"
        elif same.get(nm):
            row["role"] = "outside_assembly" if in_assembly.get(base) is False else "subunit"
        else:
            row["role"] = "outside_assembly" if in_assembly.get(base) is False else "other_molecule"

    # ---- the partners
    by_name = {r["chain"]: r for r in rows}
    asm_word = f"assembly {chosen.name}" if chosen is not None else "the file"
    chosen_partners: list[str] = []
    partners_why = ""
    if target is None:
        partners_why = "no target chain was found, so none were chosen"
    elif mode == "none":
        partners_why = "as asked: the binder is not checked against any other chain"
    elif mode == "given":
        bad = [c for c in asked if c not in by_name]
        if bad:
            problem = problem or (f"{Path(path).name}: partner chain {', '.join(map(repr, bad))} is not a polymer chain of the file "
                                  f"(chains: {', '.join(by_name)}).")
        elif target in asked:
            problem = problem or f"chain {target} is the target and cannot also be its partner."
        else:
            chosen_partners = list(dict.fromkeys(asked))
            partners_why = "as asked"
            for c in chosen_partners:
                if by_name[c]["role"] == "outside_assembly":
                    notes.append(f"chain {c} is not part of {asm_word}, but it was asked for and is counted")
                if by_name[c]["kind"] != "protein":
                    notes.append(f"chain {c} is {by_name[c]['kind']}: it is counted as a partner as asked")
    elif mode == "all":
        chosen_partners = [r["chain"] for r in rows if r["role"] in ("subunit", "other_molecule")]
        partners_why = f"every other polymer chain of {asm_word}" if chosen_partners else f"{asm_word} holds no other polymer chain"
    else:
        chosen_partners = [r["chain"] for r in rows if r["role"] == "subunit"]
        if chosen_partners:
            copy_n = sum(1 for c in chosen_partners if by_name[c]["generated"])
            partners_why = (f"the same protein as the target, in {asm_word}"
                            + (f"; {copy_n} made from the assembly's symmetry operators, because the file holds only part of it" if copy_n else ""))
        elif len(rows) == 1:
            partners_why = "the file holds no other polymer chain"
        else:
            partners_why = f"the file holds no other chain with the target's sequence in {asm_word}"
    # subunits first: the first partner is the one the oligomer's two-fold operator is read from
    chosen_partners.sort(key=lambda c: (by_name[c]["role"] != "subunit", list(by_name).index(c)))
    for row in rows:
        row["used"] = row["chain"] in chosen_partners
        nm, role = row["chain"], row["role"]
        base = row.get("copy_of", nm)
        if role == "unassigned":
            continue
        if row["touches"] is None:
            touch = ""
        elif row["touches"]:
            touch = f" It touches the target at {row['touches']} residues."
        else:
            touch = f" It does not touch the target (closest {row['closest']} A)."
        if role == "target":
            row["why"] = why
        elif role == "subunit":
            same_in = f"the same protein as the target and part of {asm_word}" if in_assembly.get(base) else \
                "the same protein as the target (the file has no assembly record to say more)"
            if row["generated"]:
                row["why"] = (f"A copy of chain {base} that {asm_word} makes with operator {row['operator']}. The file holds only part of the "
                              f"assembly, so the copy was generated from the records.{touch}")
            else:
                row["why"] = f"{same_in[0].upper() + same_in[1:]}.{touch}"
        elif role == "other_molecule":
            row["why"] = (f"A different molecule{(': ' + row['name']) if row['name'] else ''}. The file cannot say whether it belongs to the real "
                          f"complex or was added to crystallise the target, so it is not counted as a partner unless you include it.{touch}")
        else:
            row["why"] = (f"{'The same protein as the target' if same.get(nm) else 'A different molecule'}, but not part of {asm_word}: "
                          f"a neighbour in the crystal, not a partner.{touch}")
        if row["used"] and role in ("other_molecule", "outside_assembly"):
            row["why"] = "Counted as a partner because it was asked for. " + row["why"]
    for r in rows:
        r["generated"] = bool(r["generated"])

    # ---- everything that is not a polymer chain
    gl = list(glycans) if glycans is not None else _find_glycans(sf)
    non_poly, n_water, ligand_notes = _non_polymers(st, tgt_tree, inst_atoms.get(target) if target else None)
    measured = {target, *chosen_partners} if target else set()
    glycan_rows = [{"chain": g.owner, "site": g.site, "sugars": list(g.resnames),
                    "measured": g.owner in measured} for g in gl]
    notes.extend(ligand_notes)

    # ---- the summary the log and the service show first
    parts = []
    if target:
        parts.append(f"Target: chain {target} ({why})")
        parts.append("Partners: " + (_wlist(chosen_partners) if chosen_partners else "none") + f" ({partners_why})")
        left = [r for r in rows if r["role"] in ("other_molecule", "outside_assembly") and not r["used"]]
        if left:
            parts.append("Not counted: " + "; ".join(
                f"chain {r['chain']} ({'another molecule' if r['role'] == 'other_molecule' else 'not in the assembly'}"
                + (f", {r['name']}" if r["name"] else "") + (", touches the target" if r["touches"] else "") + ")" for r in left))
    cases = []
    if target:
        if not any(r["role"] != "target" for r in rows):
            cases.append("single_chain")
        if any(r["role"] == "subunit" and not r["generated"] for r in rows):
            cases.append("subunits_in_file")
        if any(r["role"] == "subunit" and r["generated"] for r in rows):
            cases.append("subunits_by_symmetry")
        if any(r["role"] == "other_molecule" for r in rows):
            cases.append("other_molecules")
        if any(r["role"] == "outside_assembly" for r in rows):
            cases.append("crystal_neighbours")
        if alternatives:
            cases.append("target_ambiguous")
    data = {
        "file": Path(path).name, "entry": pdb, "resolution": sf.resolution,
        "construct": {"length": len(cseq)} if cseq else None,
        "target": {"chain": target, "how": how, "why": why, "alternatives": alternatives, "warning": warning},
        "chains": rows,
        "partners": {"mode": mode, "chains": chosen_partners, "why": partners_why},
        "assembly": asm_info,
        "non_polymers": non_poly, "waters": n_water, "glycans": glycan_rows,
        "cases": cases, "notes": notes, "summary": ". ".join(parts) if parts else (problem or ""),
        "problem": problem, "name_the_chain": name_it,
    }
    return ReferenceChains(data, sf, {nm: (base, R, t) for nm, base, _op, R, t in generated}, gl)


def _non_polymers(st, tgt_tree, tgt_atoms) -> tuple[list[dict[str, Any]], int, list[str]]:
    """The molecules of the file that are not polymer chains: ``(rows, number of waters, notes)``. Ions, ligands and additives are
    grouped by name with a count and the target residues they touch (within 4 A); N-glycans are reported separately. None of them is
    ever a partner: the clash check would measure a detergent molecule or a crystallisation additive as if it were part of the protein."""
    gemmi = S._gemmi()
    groups: dict[str, dict[str, Any]] = {}
    waters = 0
    for ch in st[0]:
        for res in ch:
            ent = st.get_entity(res.entity_id) if res.entity_id else None
            etype = ent.entity_type if ent is not None else None
            if res.name in S._WATER or etype == gemmi.EntityType.Water:
                waters += 1
                continue
            if etype == gemmi.EntityType.Branched:
                continue                                                      # an N-glycan (listed under glycans)
            if etype == gemmi.EntityType.Polymer or (etype is None and res.het_flag != "H") or \
                    res.name in S.THREE_TO_ONE or (res.name in S._MODIFIED and res.het_flag != "H"):
                continue
            xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name.upper() != "H"], dtype=float).reshape(-1, 3)
            if not len(xyz):
                continue
            g = groups.setdefault(res.name, {"count": 0, "chains": set(), "xyz": [], "heavy": len(xyz)})
            g["count"] += 1
            g["chains"].add(ch.name)
            g["xyz"].append(xyz)
    rows: list[dict[str, Any]] = []
    notes: list[str] = []
    for name, g in sorted(groups.items(), key=lambda kv: (-kv[1]["count"], kv[0])):
        kind = "ion" if g["heavy"] == 1 else "sugar" if name in _SUGAR_NAMES else "ligand"
        touches: list[int] = []
        closest = None
        if tgt_tree is not None and tgt_atoms is not None:
            for xyz in g["xyz"]:
                d, ix = tgt_tree.query(xyz, distance_upper_bound=_TOUCH_LIGAND)
                near = np.isfinite(d)
                if near.any():
                    touches.extend(int(x) for x in tgt_atoms.resnum[ix[near]])
            closest = round(float(min(tgt_tree.query(xyz)[0].min() for xyz in g["xyz"])), 1)
        touches = sorted(set(touches))
        rows.append({"name": name, "kind": kind, "count": g["count"], "chains": sorted(g["chains"]),
                     "touches_target": touches[:12], "n_touch": len(touches), "closest": closest,
                     "hint": "named like a buffer, salt or cryoprotectant" if name in _ADDITIVES else ""})
    held = [r for r in rows if r["n_touch"] and r["kind"] in ("ion", "ligand", "sugar") and not r["hint"]]
    if held:
        notes.append("bound to the target and not measured: " + ", ".join(
            f"{r['name']} x{r['count']} (residues {', '.join(map(str, r['touches_target'][:4]))}{' ...' if r['n_touch'] > 4 else ''})" for r in held[:4]))
    return rows, waters, notes


def describe_reference(reference: str | Path, construct_sequence: str | None = None, **kw: Any) -> dict[str, Any]:
    """``analyse_reference(...).data``: what the reference holds and what would be done with each chain, as plain JSON. This is the
    table the web service shows before a build; ``binder-review chains`` prints it."""
    return analyse_reference(reference, construct_sequence, **kw).data


def format_reference_chains(data: Mapping[str, Any]) -> list[str]:
    """The composition of a reference as lines of text: one for the target, one for the partners, one per chain that is not counted,
    and the non-polymer molecules. What the build log and ``binder-review chains`` print."""
    out: list[str] = []
    if data.get("problem"):
        out.append(f"PROBLEM: {data['problem']}")
    t = data["target"]
    if not t["chain"]:
        for r in data["chains"]:
            out.append(f"  chain {r['chain']}: {r['kind']}, {r['n_residues']} residues" + (f" ({r['name']})" if r["name"] else ""))
        return out
    out.append(f"target: chain {t['chain']}: {t['why']}")
    if t.get("warning"):
        out.append(f"WARNING: {t['warning']}")
    got = data["partners"]["chains"]
    out.append(f"partners: {_wlist(got) if got else 'none'} ({data['partners']['why']})")
    asm = data["assembly"]
    if asm["id"]:
        extra = f" ({asm['oligomeric']})" if asm["oligomeric"] else ""
        n = len(asm["generated"])
        out.append(f"assembly {asm['id']}{extra}" + (f"; made {n} cop{'y' if n == 1 else 'ies'} from its symmetry operators: {', '.join(asm['generated'])}" if n else ""))
    elif not asm["records"]:
        out.append("the file has no assembly records: every chain with the target's sequence counts as a subunit")
    for r in data["chains"]:
        if r["role"] == "target":
            continue
        used = "counted" if r["used"] else "not counted"
        what = {"subunit": "subunit", "other_molecule": "other molecule", "outside_assembly": "not in the assembly"}[r["role"]]
        out.append(f"  chain {r['chain']}: {what}, {used}" + (f" ({r['name']})" if r["name"] else "") + (f"; {r['hint']}" if r["hint"] else ""))
    if data["non_polymers"] or data["waters"]:
        bits = [f"{r['name']} x{r['count']}" + (" (touches the target)" if r["n_touch"] else "") for r in data["non_polymers"][:10]]
        if data["waters"]:
            bits.append(f"{data['waters']} waters")
        out.append("never counted as partners: " + ", ".join(bits))
    for n in data.get("notes") or []:
        out.append(f"note: {n}")
    return out


def build_frame(reference_cif: str | Path | None, construct_pdb: str | Path, runs: Any, uniprot_acc: str | None = None, *,
                chain: str | None = None, partner_chains: Sequence[str] | str | None = None, assembly_id: str | None = None,
                config: Mapping[str, Any] | None = None, use_alphafold: bool = True,
                reference_id: str | None = None) -> TargetFrame:
    """Build the target frame.

    ``reference_cif``: the full-length reference structure. ``construct_pdb``: the file BindCraft was given (numbering 1
    to N). ``runs``: {run id: RunData | RunSpec | {"hotspots": [...], "numbers": [...]} | [hotspot, ...]} with hotspots
    and residues in construct numbers; a run given by hotspots only is assumed to carry the whole construct.
    ``uniprot_acc``: the UniProt accession of the target, optional. With it the numbering is verified (construct against
    the reference chain and UniProt, reference against UniProt and SIFTS) and the topology, domains, glycosylation,
    disulfides and a membrane plane are annotated from the refdata cache (a missing record becomes a warning and the
    part that needs it is left out); without it the frame comes from the reference structure and the construct alone.
    The construct offset is found by sequence, not assumed. ``chain``: the chain of the target in the reference (None: the
    chain that carries the construct's sequence; the frame records which one it took and why, in ``meta["composition"]``).
    ``partner_chains``: what the binder is measured against besides the target chain. None or ``"auto"``: the other subunits, that
    is chains of the same protein that belong to the target's assembly, and the copies that the assembly's symmetry operators make
    from them when the file holds only part of the assembly; ``"none"`` or ``()``: nothing; ``"all"``: every other polymer chain of
    the assembly; a list: exactly those chains. Another molecule is never counted unless it is named, and ligands, ions, additives
    and waters never are (see :func:`analyse_reference`). ``assembly_id``: the assembly of the file the subunits come from (default:
    the first one that holds the target chain). ``reference_id``: the id the reference is named by in labels (default: the file stem).
    Raises FrameError when the files do not fit.

    ``reference_cif=None`` builds the frame without a full-length reference: the construct is the frame (chain A of the construct file),
    designs are superposed on it so that they share one frame, numbering stays the construct's, and nothing outside the construct is
    known. ``frame.reference_given`` is then false and ``context()`` returns no metric at all (see docs/target_frame.md, "Without a
    reference"). There are no other subunits, no sugars and no membrane plane, and the cut ends of the construct itself are unknown:
    an accession can still be given, for the numbering check, the topology and the specificity panel."""
    from . import refdata as rd

    acc = (uniprot_acc or "").strip() or None
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    for sect, vals in (config or {}).items():
        cfg.setdefault(sect, {}).update(vals)
    geo = cfg["geometry"]
    probe, points = float(geo["sasa_probe"]), int(geo["sasa_points"])
    specs = _as_run_specs(runs)
    warnings: list[str] = []
    notes: list[str] = []
    # -- reference structure (or, without one, the construct: the only structure that is known)
    given = reference_cif is not None
    if given:
        frame_file = reference_cif
    else:
        frame_file = construct_pdb
        if chain not in (None, "A", "auto") or _parse_partner_request(partner_chains)[0] not in ("auto", "none") or reference_id:
            raise FrameError("without a reference structure the construct is the frame: its target chain is A, there is no other subunit "
                             "and no reference id (chain, partner chains and reference id name parts of a reference structure)")
        partner_chains = []
        chain = "A"
    sf = S.read_structure(str(frame_file))
    construct = _read_construct(str(construct_pdb))
    construct_seq = "".join(str(a) for a in construct.residues.one)
    glycans_in_file = _find_glycans(sf) if given else []
    comp: ReferenceChains | None = None
    if given:
        # which chain is the target and which chains are partners is decided in one place, with its reasons, so that the command
        # line, the service and the frame agree and a person can read what was done (meta["composition"])
        comp = analyse_reference(reference_cif, construct_seq, chain=None if chain in (None, "auto") else chain, partner_chains=partner_chains,
                                 assembly_id=assembly_id, sf=sf, glycans=glycans_in_file, construct_atoms=construct)
        if comp.problem:
            raise FrameError(comp.problem.rstrip(".") + ("; name the chain of the target with --reference-chain"
                                                         if comp.data["name_the_chain"] and "name the chain" not in comp.problem else ""))
        chain, partner_chains = str(comp.target), list(comp.partners)
        if comp.data["target"]["warning"]:
            # a chain that was named but does not carry the construct: say which one does, before the numbering vote can fail on its own terms
            held = (next(r for r in comp.data["chains"] if r["chain"] == chain).get("construct") or {}).get("found", 0)
            better = [r for r in comp.data["chains"] if r["chain"] != chain and not r["generated"]
                      and (r.get("construct") or {}).get("found", 0) >= 0.9 * len(construct_seq)]
            raise FrameError(f"the construct matches the reference chain at only {held} of {len(construct_seq)} residues"
                             + (f"; chain {better[0]['chain']} carries it ({better[0]['construct']['found']} residues): name it with --reference-chain"
                                if better else ""))
        notes.append("reference chains: " + comp.data["summary"] + ".")
        notes.extend(comp.data["notes"])
    elif sf.kinds.get(chain) != "polymer":
        raise FrameError(f"{frame_file}: chain {chain!r} is not a polymer chain (chains: {sf.polymer_chains()}); the construct's target must be chain A")
    A_all = sf.groups[chain].heavy()
    keep_a, shared, later = _first_residue_per_number(A_all)
    A = A_all.select(keep_a) if shared else A_all                          # one residue per number: the tables and the arrays are read by number
    X = A_all.select(~keep_a) if shared else S.Atoms.empty()               # the later ones (150A) are kept apart: part of the chain, not of its tables
    P = {c: comp.atoms(c) for c in partner_chains} if comp else {}
    glycans = list(glycans_in_file)                      # sugars of the construct file are never a measurement of the real protein
    for c in partner_chains:
        if comp and c in comp._copies:
            glycans.extend(comp.glycans(c))              # the sugars of a copy that the assembly made move with it
    pdb = reference_label(reference_id or Path(reference_cif).stem) if given else Path(construct_pdb).stem
    if not given:
        notes.append("no reference structure was given: the construct is the frame. Designs are superposed on it so that they share one frame, "
                     "and nothing outside the construct is known, so the checks against the other subunits, the rest of the protein, cut ends, "
                     "sugars and the membrane were not run (their values are absent, not zero)")
    if not acc:
        notes.append("no UniProt accession was given: the numbering is that of the reference structure and no UniProt, SIFTS or "
                     "AlphaFold record was read (no domain strip, no named topology, no membrane plane)")
    # -- sequences and numbering
    ref_one = {int(n): str(a) for n, a in zip(A.residues.number, A.residues.one)}
    if shared:
        warnings.append(f"reference chain {chain}: {len(shared)} residue numbers are carried by more than one residue (insertion codes?), "
                        f"first {shared[:5]}; designs are matched to the reference by residue number, so these residues are left out of the placement")
        warnings.append(f"reference chain {chain}: of the residues that share a number only the first is kept in the frame's residue tables and in the "
                        f"target structure of the Inspector; the later ones ({', '.join(later[:5])}{' ...' if len(later) > 5 else ''}) are not drawn and "
                        "have no values of their own, but their atoms count as part of the rest of the chain wherever their number is outside a run's target")
    uni_seq = ""
    rec = _safe(lambda: rd.fasta_sequence(rd.get_uniprot_fasta(acc)), warnings, f"UniProt FASTA {acc}", "") if acc else ""
    if rec:
        uni_seq = rec
    c_numbers = [int(n) for n in construct.residues.number]
    c_one = [str(a) for a in construct.residues.one]
    d_ref, _ = _vote_offset(zip(c_numbers, c_one), ref_one)     # reference number = construct number + d_ref
    identical = sum(1 for n, a in zip(c_numbers, c_one) if ref_one.get(n + d_ref) == a)
    if identical < 0.9 * len(c_numbers):
        raise FrameError(f"the construct matches the reference chain at only {identical} of {len(c_numbers)} residues")
    # reference author numbering against UniProt: the sequence decides the offset, SIFTS is recorded as a second witness
    ref_off, ref_vs_uni = 0, None
    if uni_seq:
        ref_off, same = _vote_offset(ref_one.items(), {i + 1: a for i, a in enumerate(uni_seq)})
        ref_vs_uni = {"n": len(ref_one), "identical": same, "offset": ref_off}
        if same < 0.9 * len(ref_one):
            raise FrameError(f"reference chain {chain} matches UniProt {acc} at only {same} of {len(ref_one)} residues")
    sifts = (_safe(lambda: rd.get_sifts(pdb, "uniprot_segments"), warnings, "SIFTS uniprot_segments", None)
             if acc and given and rd.is_pdb_id(pdb) else None)
    sifts_info: dict[str, Any] | None = None
    if sifts:
        segs = [m for m in sifts["UniProt"].get(acc, {}).get("mappings", []) if m["chain_id"] == chain]
        if segs:
            sifts_info = {"segments": [[int(m["start"]["author_residue_number"]), int(m["end"]["author_residue_number"]),
                                        int(m["unp_start"]), int(m["unp_end"])] for m in segs],
                          "author_equals_uniprot": all(m["unp_start"] - m["start"]["author_residue_number"] == 0 and
                                                       m["unp_end"] - m["end"]["author_residue_number"] == 0 for m in segs)}
    if ref_off != 0:
        raise FrameError(f"{'reference author' if given else 'construct'} numbering differs from UniProt by {ref_off}: not supported (all outputs would "
                         "need renumbering); " + ("give a reference structure whose residue numbers follow UniProt, or leave the accession out" if given else
                         "give a construct whose residue numbers follow UniProt, or leave the accession out"))
    offset = d_ref + ref_off                                    # construct -> UniProt
    construct_range = (c_numbers[0] + offset, c_numbers[-1] + offset)
    uni_identical = None
    if uni_seq:
        uni_identical = sum(1 for n, a in zip(c_numbers, c_one) if 1 <= n + offset <= len(uni_seq) and uni_seq[n + offset - 1] == a)
    length = len(uni_seq) if uni_seq else int(max(ref_one))
    # Without a reference and without an accession nothing says where the protein ends, so the native range is the construct's own span:
    # the construct's ends are then not reported as cuts (we cannot tell), only the stretches a run deleted inside it are.
    native_all = [1, length] if given or uni_seq else [construct_range[0], construct_range[1]]
    numbering = {"accession": acc or "", "length": length, "native_range": native_all, "offset": offset, "reference_offset": ref_off,
                 "construct_first": c_numbers[0], "construct_last": c_numbers[-1], "n_construct": len(c_numbers),
                 "construct_numbers": c_numbers, "uniprot_first": construct_range[0], "uniprot_last": construct_range[1],
                 "reference_first": int(min(ref_one)), "reference_last": int(max(ref_one)), "uniprot_sequence": uni_seq,
                 "checks": {"construct_vs_reference": {"n": len(c_numbers), "identical": identical},
                            "construct_vs_uniprot": None if uni_identical is None else {"n": len(c_numbers), "identical": uni_identical},
                            "reference_vs_uniprot": ref_vs_uni, "sifts": sifts_info}}
    if shared:
        numbering["ambiguous"] = shared                                     # numbers the placement leaves out (carried by more than one residue)
    # -- per run
    run_rows: list[dict[str, Any]] = []
    A_by_num = {int(n): i for i, n in enumerate(A.residues.number)}
    for sp in specs:
        off = offset if sp.offset is None else int(sp.offset)
        nums = sorted(set(sp.numbers)) if sp.numbers else list(c_numbers)
        outside = [n for n in nums if n not in set(c_numbers)]
        if outside and sp.offset is None:
            raise FrameError(f"run {sp.run_id}: target residues {outside[:5]} are not in the construct file")
        present = [n + off for n in nums]
        lost = [u for u in present if u not in A_by_num]
        if lost:
            raise FrameError(f"run {sp.run_id}: residues {lost[:5]} are not in the reference chain")
        bad_h = [h for h in sp.hotspots if (h + off) not in set(present)]
        if bad_h:
            raise FrameError(f"run {sp.run_id}: hotspots {bad_h} are not in the run's target")
        seq_ok = None
        if sp.sequence is not None and len(sp.sequence) == len(nums):
            seq_ok = all(ref_one.get(n + off) == a for n, a in zip(nums, sp.sequence))
        rr = _ranges(present)
        removed = [[a[1] + 1, b[0] - 1] for a, b in zip(rr, rr[1:])]
        native = list(native_all)
        run_rows.append({
            "run_id": sp.run_id, "label": sp.label or sp.run_id, "offset": off, "hotspots_construct": sorted(sp.hotspots),
            "hotspots": sorted(h + off for h in sp.hotspots), "resi_ranges": rr, "removed": removed, "first": present[0],
            "last": present[-1], "n_resi": len(present), "cut_points": cut_points_for(rr, native),
            "removed_native": _complement(rr, native[0], native[1]),   # everything of the native chain the run lacks
            "sequence_matches_reference": seq_ok, "consistent": sp.consistent, "epitope": []})
    present_union = sorted(set().union(*(_expand(r["resi_ranges"]) for r in run_rows)))
    union_range = (present_union[0], present_union[-1])
    # -- annotations
    ann = _annotations(A, glycans, sf.links if given else [], chain, partner_chains, uni_seq, acc, pdb, warnings, union_range,
                       tuple(native_all), reference_given=given)
    sites_in = sorted(g.site for g in glycans if g.owner == chain and union_range[0] <= g.site <= union_range[1])
    uni_sites = [g["site"] for g in ann["glycosylation"] if g["in_construct"]]
    glyc_sites = sorted(set(sites_in) | set(uni_sites)) if given else []   # the sites UniProt lists stay on the strip, but no distance to one is measured
    anchors = []
    for s in glyc_sites:
        try:
            anchors.append(A.xyz[A.atom(s, "ND2")])
        except KeyError:
            raise FrameError(f"glycan site {s} has no Asn ND2 in the reference chain") from None
    # -- secondary structure and SASA
    ss = S.assign_ss(A)
    one = A.residues.one
    sasa_A = S.sasa_atoms(A.xyz, A.element, probe, points)
    sasa_A_res = A.per_residue(sasa_A)
    rsa = np.clip(S.relative_sasa(sasa_A_res, one), 0.0, 1.0)
    ga = [g for g in glycans if g.owner == chain]
    gp = [g for g in glycans if g.owner in partner_chains]
    Bpoly = S.Atoms.concat([P[c] for c in partner_chains]) if partner_chains else S.Atoms.empty()
    s_chain_gly = _sasa_parts([A, _glycan_atoms(ga)], probe, points)[0]
    s_dimer = _sasa_parts([A, Bpoly, _glycan_atoms(ga), _glycan_atoms(gp)], probe, points)[0]
    buried = A.per_residue(s_chain_gly) - A.per_residue(s_dimer)
    arrays: dict[str, np.ndarray] = {}
    _pack("A", A, arrays)
    if len(X):
        _pack("X", X, arrays)
    arrays["A.sasa"] = sasa_A
    for c in partner_chains:
        _pack(f"P_{c}", P[c], arrays)
    for k, g in enumerate(glycans):
        _pack(f"G{k}", g.atoms, arrays)
        arrays[f"G{k}.anchor"] = np.asarray(g.anchor, dtype=float)
    arrays["glycan.anchor_xyz"] = np.array(anchors, dtype=float).reshape(-1, 3)
    arrays["res.number"] = A.residues.number.astype(np.int64)
    arrays["res.one"] = np.asarray(one, dtype="<U1")
    arrays["res.ss"] = np.asarray(list(ss), dtype="<U1")
    arrays["res.ca"] = A.ca_xyz()
    arrays["res.sasa"] = sasa_A_res
    arrays["res.rsa"] = rsa
    arrays["res.sasa_chain_gly"] = A.per_residue(s_chain_gly)
    arrays["res.sasa_dimer"] = A.per_residue(s_dimer)
    arrays["res.buried_partner"] = buried
    # per run: SASA of the run's own target alone, and with the rest of the chain (and its sugars) present
    for row in run_rows:
        present = _expand(row["resi_ranges"])
        mask = np.isin(A.resnum, np.fromiter(present, dtype=np.int64))
        own = A.select(mask)
        rest_gl = [g for g in ga if g.site not in present]
        rest = S.Atoms.concat([A.select(~mask), _glycan_atoms(rest_gl)])
        s_alone = S.sasa_atoms(own.xyz, own.element, probe, points)
        s_with = _sasa_parts([own, rest], probe, points)[0]
        full = np.full(len(one), np.nan)
        full_rest = np.full(len(one), np.nan)
        sel_idx = np.nonzero(np.isin(A.residues.number, np.fromiter(present, dtype=np.int64)))[0]
        full[sel_idx] = own.per_residue(s_alone)
        full_rest[sel_idx] = own.per_residue(s_with)
        arrays[f"run.{row['run_id']}.sasa_alone"] = full
        arrays[f"run.{row['run_id']}.sasa_rest"] = full_rest
        # epitope: hotspots plus residues within epitope_radius of a hotspot residue, inside this run's target
        hs = set(row["hotspots"])
        tree = cKDTree(own.xyz)
        hit = set(hs)
        hot_idx = np.nonzero(np.isin(own.resnum, np.fromiter(hs, dtype=np.int64)))[0]
        for ball in tree.query_ball_point(own.xyz[hot_idx], float(geo["epitope_radius"])):
            hit.update(int(own.resnum[a]) for a in ball)
        row["epitope"] = sorted(hit & present)
    # -- membrane: plane 6 A short of the chain end that leads to the transmembrane anchor (only when UniProt names one;
    #    otherwise the plane is a placeholder that no check reads: ``known`` is false and no membrane distance is reported)
    ca_by = {int(n): c for n, c in zip(A.residues.number, arrays["res.ca"])}
    centroid = np.mean([ca_by[u] for u in present_union], axis=0)
    tm = [t for t in ann["topology"] if t["label"] == "Transmembrane"] if given else []   # the plane is built from a chain end of the reference
    first_ref, last_ref = int(A.residues.number.min()), int(A.residues.number.max())
    end_ref = last_ref if tm and tm[0]["start"] > (first_ref + last_ref) / 2 else first_ref   # the end nearer the anchor in sequence
    short = float(geo["membrane_offset"])
    plane = plane_from_points(centroid, ca_by[end_ref], short)
    anchor = f" (UniProt {tm[0]['start']}-{tm[0]['end']})" if tm else ""
    membrane = {"centroid": plane["centroid"].tolist(), "axis": plane["axis"].tolist(), "offset": plane["offset"],
                "end_residue": end_ref, "short_by": short, "known": bool(tm),
                "label": "membrane-side plane (approximate)" if tm else "no membrane plane (no transmembrane region known)",
                "source": (f"centroid of the construct CA atoms towards CA of residue {end_ref} of the reference chain, the end "
                           f"that leads to the transmembrane anchor{anchor}; the plane is {short:g} Å short of that end") if tm else
                          ("no UniProt topology was read, so no transmembrane region is known: the plane is a placeholder and no membrane distance is reported"
                           if given else "no reference structure was given, so the chain end that leads to a membrane is unknown: the plane is a placeholder and no membrane distance is reported")}
    if not given:
        membrane["label"] = "no membrane plane (no reference structure was given)"
    asm_used = comp.data["assembly"]["id"] if comp else None
    dimer = _dimer_report(str(reference_cif), chain, partner_chains, A, P, glycans, probe, points, float(geo["olig_min_area"]),
                          asm_used or "", subunit_first=bool(comp and partner_chains and comp.role(partner_chains[0]) == "subunit")) if given else {}
    if "two_fold" in dimer:
        ang = float(np.degrees(np.arccos(abs(np.clip(plane["axis"] @ np.asarray(dimer["two_fold"]["axis"]), -1, 1)))))
        pt = np.asarray(dimer["two_fold"]["point"])
        ax2 = np.asarray(dimer["two_fold"]["axis"])
        sign = 1.0 if (ca_by[end_ref] - pt) @ ax2 >= 0 else -1.0
        membrane["diagnostics"] = {
            "angle_to_two_fold_axis": round(ang, 1),
            "chain_end_height_on_two_fold_axis": round(float(sign * (ca_by[end_ref] - pt) @ ax2), 1),
            "construct_height_on_two_fold_axis": round(float(sign * (centroid - pt) @ ax2), 1)}
    if use_alphafold and tm and acc and given:
        af = _safe(lambda: rd.get_alphafold(acc, with_pae=False), warnings, "AlphaFold model", None)
        if af:
            try:
                m = S.read_structure(str(af["model_path"])).groups["A"]
                sup = S.superpose_residues(m, A)
                cam = m.ca_xyz() @ sup.R.T + sup.t
                nums = m.residues.number
                tmask = (nums >= tm[0]["start"]) & (nums <= tm[0]["end"])
                beyond = ((cam - plane["centroid"]) @ plane["axis"]) - plane["offset"]
                membrane.setdefault("diagnostics", {}).update({
                    "alphafold_fit_rmsd": round(sup.rmsd, 2), "alphafold_fit_n": int(sup.n),
                    "alphafold_tm_beyond_plane": round(float(beyond[tmask].mean()), 1),
                    "alphafold_tm_plddt": round(float(m.plddt[tmask].mean()), 1),
                    "alphafold_note": "AlphaFold DB model placed by CA superposition over the reference residues; the "
                                      "anchor and linker positions are low-confidence and only a plausibility check"})
            except (S.StructureError, KeyError, ValueError) as ex:
                warnings.append(f"AlphaFold membrane check skipped: {ex}")
    meta = {
        "schema": FRAME_SCHEMA, "generated": {"tool": "binder_review.target_frame", "date": _dt.date.today().isoformat()},
        "inputs": {"reference": ({"file": Path(reference_cif).name, "sha256": _sha256(reference_cif)} if given else {"file": None, "sha256": None}),
                   "construct": {"file": Path(construct_pdb).name, "sha256": _sha256(construct_pdb)}},
        "config": cfg, "chain": chain, "partner_chains": list(partner_chains),
        "composition": comp.data if comp else None,
        "reference": ({"pdb": pdb, "assembly": asm_used or None, "resolution": sf.resolution, "given": True, "class": classify_reference(reference_cif),
                       "label": f"{pdb} chain {chain}" + (f", assembly {asm_used}" if asm_used else "")
                       + (" (dimer)" if dimer.get("oligomeric_count") == "2" else "")}
                      if given else
                      {"pdb": pdb, "assembly": None, "resolution": None, "given": False, "class": classify_reference(None),
                       "label": f"{pdb} (the construct: no full-length reference was given)"}),
        "numbering": numbering, "runs": run_rows, "annotations": ann, "dimer": dimer, "membrane": membrane,
        "glycans": [g.meta() for g in glycans], "glycan_sites_in_construct": glyc_sites,
        "warnings": warnings, "notes": notes,
    }
    return TargetFrame(meta, arrays)


# ----------------------------------------------------------------------------------------------------------------
# placing a design
# ----------------------------------------------------------------------------------------------------------------
def place_design(frame: TargetFrame, design: Any, run_id: str, *, fit: str | None = None,
                 n_binder_chains: int | None = None) -> Placement:
    """Place a design on the reference chain.

    ``design``: a structure.Complex or the path of a BindCraft PDB. The design's chain A is superposed on the reference
    chain over the residues they share, matched by number (construct number + offset = reference number) with CA atoms;
    deleted loops and any other absent residue simply do not take part. ``fit`` (default: the frame's ``geometry.fit``,
    "interface"):

    * "interface" weights every shared CA by its distance to the nearest binder heavy atom (``interface_weights``: 1 up to
      ``iface_full``, Gaussian beyond, never below ``iface_floor``) and drops near residues that deviate by more than
      ``iface_cutoff`` before refitting (``interface_fit``), so the epitope sits where it does on the reference even when the rest of
      the target moved against it (a hinge, a flexible loop). With fewer than ``iface_min_n`` residues near the binder, or no binder, the
      whole-target fit is used and ``fit_note`` says why.
    * "all" uses every shared CA with the same weight (docs/CONTRACTS.md).
    * "core" refits on the residues a trimmed fit keeps.

    ``Placement.fit`` is the fit actually used. ``Placement.rmsd`` is the RMSD of all shared CA under the transform used,
    ``iface_rmsd`` the one over the residues near the binder (weight at least 0.5), ``global_rmsd`` what the whole-target fit gives and
    ``core_rmsd`` the trimmed fit's own RMSD."""
    geo = {**DEFAULT_CONFIG["geometry"], **frame.config.get("geometry", {})}
    mode = geo["fit"] if fit is None else fit
    if mode not in PLACEMENT_FITS:
        raise ValueError(f"fit must be one of {', '.join(PLACEMENT_FITS)}, not {mode!r}")
    cx = design if isinstance(design, S.Complex) else S.load_complex(str(design), n_binder_chains=n_binder_chains)
    tgt = cx.target
    off = int(frame.run_info(run_id)["offset"])
    ca_t = tgt.ca_by_number()
    ref_ca, ref_one = frame.ref_ca, frame.ref_one
    one_t = {int(n): str(a) for n, a in zip(tgt.residues.number, tgt.residues.one)}
    common = [n for n in sorted(ca_t) if (n + off) in ref_ca]
    if len(common) < 10:
        raise FrameError(f"{getattr(cx, 'path', 'design')}: only {len(common)} target residues match the reference by number")
    P = np.array([ca_t[n] for n in common])
    Q = np.array([ref_ca[n + off] for n in common])
    Rg, tg, global_rmsd = S.kabsch(P, Q)
    Rc, tc, keep = trimmed_fit(P, Q, float(geo["core_cutoff"]), int(geo["core_cycles"]))
    core_rmsd = float(np.sqrt(((P[keep] @ Rc.T + tc - Q[keep]) ** 2).sum(1).mean()))
    # where the binder is, in the design's own frame: the distance of every shared CA to its nearest heavy atom
    binder = cx.binder.heavy() if len(cx.binder) else cx.binder
    dist = cKDTree(binder.xyz).query(P)[0] if len(binder) else np.full(len(P), np.inf)
    weights = interface_weights(dist, float(geo["iface_full"]), float(geo["iface_decay"]), float(geo["iface_floor"]))
    used, note = mode, ""
    if mode == "interface":
        R, t, _, _, ok = interface_fit(P, Q, dist, float(geo["iface_full"]), float(geo["iface_decay"]), float(geo["iface_floor"]),
                                       float(geo["iface_cutoff"]), int(geo["iface_min_n"]))
        if not ok:
            used, R, t = "all", Rg, tg
            n_near = int((weights >= IFACE_NEAR_WEIGHT).sum())
            note = ("the structure has no binder atoms" if not len(binder) else
                    f"only {n_near} target residues lie near the binder (at least {int(geo['iface_min_n'])} are needed)") + \
                   ": placed by the whole-target fit"
    elif mode == "all":
        R, t = Rg, tg
    else:
        R, t = Rc, tc
    dev = np.linalg.norm(P @ R.T + t - Q, axis=1)
    rmsd = float(np.sqrt((dev ** 2).mean()))
    near = weights >= IFACE_REPORT_WEIGHT
    iface_n = int(near.sum())
    iface_rmsd = float(np.sqrt((dev[near] ** 2).mean())) if iface_n else float("nan")
    nums = np.array(common, dtype=np.int64)
    mism = sum(1 for n in common if one_t.get(n) != ref_one.get(n + off))
    warn: list[str] = []
    present = frame._run_present(run_id)
    extra = [n for n in ca_t if (n + off) not in present]
    if extra:
        warn.append(f"design target has {len(extra)} residues outside run {run_id}'s target (first: {extra[0] + off})")
    missing = len(present) - len([n for n in common if (n + off) in present])
    if missing:
        warn.append(f"design target lacks {missing} residues of run {run_id}'s target")
    if mism:
        warn.append(f"{mism} target residues differ in type from the reference")
    if iface_n and iface_rmsd > float(geo["place_warn_rmsd"]):
        warn.append(f"placement RMSD {iface_rmsd:.2f} A near the binder (over {iface_n} residues) is above {geo['place_warn_rmsd']} A: the design's "
                    "target differs from the crystal chain at the interface, so context distances are less certain")
    return Placement(run_id, R, t, rmsd, len(common), nums, off, dev, core_rmsd, int(keep.sum()), nums[keep], used, mism, warn,
                     iface_rmsd, iface_n, float(global_rmsd), note)


# ----------------------------------------------------------------------------------------------------------------
# checking designs
# ----------------------------------------------------------------------------------------------------------------
@dataclass
class DesignCheck:
    """Placement and context of one design: one row of work/<project>/target_checks.json."""

    run_id: str
    design_id: str
    status: str | None
    rank: int | None
    structure: str
    model: int | None
    placement: Placement
    context: ContextResult
    contacted: list[int]
    hotspots_touched: list[int]
    hotspot_distance: dict[int, float]
    contact_deviation_max: float | None
    verdicts: dict[str, dict[str, str]]
    termini: dict[str, float] = field(default_factory=dict)   # nterm_d, cterm_d: binder end CA to the nearest target atom (A)

    def to_dict(self) -> dict[str, Any]:
        c = self.context
        pl = self.placement
        return {
            "run": self.run_id, "id": self.design_id, "status": self.status, "rank": self.rank, "structure": self.structure,
            "model": self.model,
            "placement": {"rmsd": round(pl.rmsd, 3), "n": pl.n, "core_rmsd": round(pl.core_rmsd, 3), "core_n": pl.core_n,
                          "max_dev": round(pl.max_deviation, 2), "dropped": pl.dropped_numbers(), "R": _round(pl.R.reshape(-1), 6),
                          "t": _round(pl.t, 4), "n_mismatch": pl.n_mismatch, "warnings": pl.warnings, **pl.fit_fields()},
            "metrics": {k: (None if v is None else (int(v) if isinstance(v, (int, np.integer)) else round(float(v), 3)))
                        for k, v in c.metrics.items()},
            "place_rmsd": round(pl.rmsd, 3), "glycan": {k: round(v, 2) for k, v in c.glycan.items()}, "glycan_min": c.glycan_min,
            "nterm_d": None if "nterm_d" not in self.termini else round(self.termini["nterm_d"], 3),
            "cterm_d": None if "cterm_d" not in self.termini else round(self.termini["cterm_d"], 3),
            "tc": self.contacted, "hot": self.hotspots_touched, "hot_distance": {str(k): round(v, 2) for k, v in self.hotspot_distance.items()},
            "contact_dev_max": None if self.contact_deviation_max is None else round(self.contact_deviation_max, 2),
            "cut_contacts": c.diag.get("cut_contacts", []), "dimer_buried_contacts": c.diag.get("dimer_buried_contacts", []),
            "covered": c.diag.get("covered", []), "buried_in_assembly": c.diag.get("buried_in_assembly", []),
            "clash": {"partner": c.diag.get("partner_clash"), "rest": c.diag.get("rest_clash")},
            "rest_by_region": c.diag.get("rest_n_by_region", {}), "rest_pair": c.diag.get("rest_pair"),
            "partner_pair": c.diag.get("partner_pair"), "ctx": c.ctx, "octx": c.octx, "verdicts": self.verdicts,
        }


def check_design(frame: TargetFrame, design: Any, run_id: str, *, design_id: str | None = None, status: str | None = None,
                 rank: int | None = None, model_k: int | None = None, fit: str | None = None,
                 config: Mapping[str, Any] | None = None, n_binder_chains: int | None = None) -> DesignCheck:
    """Place one design (``fit`` as in ``place_design``, default the frame's ``geometry.fit``) and run every context check on it (contacts at
    the page cut-off, 4.5 A heavy atom to heavy atom, taken in the design's own frame, then converted to UniProt)."""
    cx = design if isinstance(design, S.Complex) else S.load_complex(str(design), n_binder_chains=n_binder_chains)
    pl = place_design(frame, cx, run_id, fit=fit)
    cutoff = float(frame.config["geometry"]["contact_cutoff"])
    con = S.residue_contacts(cx.target, cx.binder, cutoff)
    contacted = sorted(int(n) + pl.offset for n in con.a_residues)
    placed = pl.apply(cx.binder)
    res = frame.context(placed, run_id, contacted, config=config)
    hs = frame.hotspots(run_id)
    hot_d: dict[int, float] = {}
    tgt_p = pl.apply(cx.target).heavy()
    tree = cKDTree(placed.heavy().xyz)
    for h in hs:
        m = tgt_p.resnum == h - pl.offset
        if m.any():
            hot_d[h] = float(tree.query(tgt_p.xyz[m])[0].min())
    dev = pl.deviation_by_uniprot()
    cdev = max((dev[c] for c in contacted if c in dev), default=None)
    sid = design_id or Path(getattr(cx, "path", "design")).stem
    res.octx["rms_head"] = round(pl.core_rmsd, 2)    # rigid part of the target on the reference (schema: octx.rms_head)
    # the ends of the binder: CA of its first and last residue to the nearest heavy atom of the design's own target (the
    # interface module computes the same two numbers; repeated here so a QA row is complete)
    ca_b = cx.binder.ca_xyz()
    ttree = cKDTree(cx.target.heavy().xyz)
    ends = {"nterm_d": float(ttree.query(ca_b[0])[0]), "cterm_d": float(ttree.query(ca_b[-1])[0])}
    return DesignCheck(run_id, sid, status, rank, Path(getattr(cx, "path", "")).name, model_k, pl, res, contacted,
                       frame.hotspots_touched(run_id, contacted), hot_d, cdev, res.verdicts(), ends)


def check_run(frame: TargetFrame, run: Any, *, run_id: str | None = None, include_rejected: bool = True, fit: str | None = None,
              progress: Any = None) -> list[DesignCheck]:
    """check_design for every design of a model.RunData that has a structure (accepted first, in rank order, then the
    rejected ones that have a structure)."""
    rid = run_id or run.run_id
    designs = list(run.accepted()) + (list(run.rejected_with_structure()) if include_rejected else [])
    out = []
    for i, d in enumerate(designs):
        out.append(check_design(frame, d.structure_path, rid, design_id=d.id, status=d.status, rank=d.rank,
                                model_k=d.structure_model, fit=fit))
        if progress:
            progress(i + 1, len(designs))
    return out


# ----------------------------------------------------------------------------------------------------------------
# report: the tables of docs/target_frame.md, from a frame and the rows of target_checks.json
# ----------------------------------------------------------------------------------------------------------------
def report_tables(frame: TargetFrame, rows: Sequence[Mapping[str, Any]]) -> str:
    """Markdown tables of the placement and context results (what docs/target_frame.md shows). ``rows`` are
    ``DesignCheck.to_dict()`` rows. Designs are named by id; no sequence appears."""
    out: list[str] = []
    runs = [r for r in frame.run_ids if any(x["run"] == r for x in rows)]
    acc = {r: [x for x in rows if x["run"] == r and x["status"] == "accepted"] for r in runs}
    lv = ("pass", "warn", "fail")

    def dist(label: str, rs: Sequence[Mapping[str, Any]]) -> str:
        r = np.array([x["placement"]["rmsd"] for x in rs])
        c = np.array([x["placement"]["core_rmsd"] for x in rs])
        cn = np.array([x["placement"]["core_n"] for x in rs])
        i = np.array([x["placement"]["iface_rmsd"] for x in rs if x["placement"].get("iface_rmsd") is not None] or [float("nan")])
        g = np.array([x["placement"].get("global_rmsd") if x["placement"].get("global_rmsd") is not None else x["placement"]["rmsd"] for x in rs])
        fits = Counter(x["placement"].get("fit", "all") for x in rs)
        warn = float(frame.config["geometry"]["place_warn_rmsd"])
        return (f"| {label} | {len(rs)} | {'/'.join(f'{k} {v}' for k, v in sorted(fits.items()))} | {np.median(r):.2f} | {np.percentile(r, 90):.2f} | {r.max():.2f} | "
                f"{np.median(i):.2f} | {np.percentile(i, 90):.2f} | {np.nanmax(i):.2f} | {int((i > warn).sum())} | {np.median(g):.2f} | {np.median(c):.2f} | {int(np.median(cn))} |")

    out += ["Placement (A, CA of the design's target on the reference chain: over all shared residues, and over the residues near the binder, "
            "which is what the context checks depend on; the whole-target fit is the smallest RMSD over all shared residues any rigid transform reaches)", "",
            "| set | n | fit | all residues median | p90 | max | near the binder median | p90 | max | above warn | whole-target fit median | core RMSD median | core residues median |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in runs:
        if acc[r]:
            out.append(dist(f"{r} accepted", acc[r]))
        rej = [x for x in rows if x["run"] == r and x["status"] == "rejected"]
        if rej:
            out.append(dist(f"{r} rejected", rej))
    out.append(dist("all", list(rows)))
    out += ["", "Context verdicts of the accepted designs (profile.checks.context thresholds)", "",
            "| run | n | check | pass | warn | fail |", "|---|---|---|---|---|---|"]
    for r in runs:
        for part in ("partner", "rest", "cut", "cover", "overall", "membrane"):
            c = Counter(x["verdicts"][part]["level"] for x in acc[r] if part in x["verdicts"])
            out.append(f"| {r} | {len(acc[r])} | {part} | " + " | ".join(str(c.get(k, 0)) for k in lv) + " |")
    out += ["", "Hotspots touched (binder heavy atom within 4.5 A of the hotspot residue, accepted designs)", "",
            "| run | hotspot | designs touching | of |", "|---|---|---|---|"]
    for r in runs:
        hs = frame.hotspots(r)
        for h in hs:
            out.append(f"| {r} | {frame.label(h)} | {sum(h in x['hot'] for x in acc[r])} | {len(acc[r])} |")
        out.append(f"| {r} | all {len(hs)} hotspots | {sum(len(x['hot']) == len(hs) for x in acc[r])} | {len(acc[r])} |")
    # where the hotspots sit: exposure, the area the rest of the chain takes, distance to the sugars and the other subunit
    pe = frame._partner_env()
    chain_a = frame.chain_a
    sites = [x for x in frame.glycans if x.owner == frame.meta["chain"] and x.site in frame.meta["glycan_sites_in_construct"]]
    out += ["", "Hotspot facts (reference chain; areas in A^2; distances from the hotspot residue to the nearest atom, in A)", "",
            "| Hotspot | Runs | Construct no. | Relative SASA | Area lost to the rest | " + " | ".join(f"To N{g.site} sugar" for g in sites)
            + " | To the other subunit |", "|---|---|" + "---|" * (4 + len(sites))]
    for h in frame.hotspots():
        xyz = chain_a.xyz[chain_a.resnum == h]
        ruls = [r for r in frame.run_ids if h in frame.hotspots(r)]
        lost = " / ".join(f"{frame.covered_area(r).get(h, 0.0):.0f}" for r in ruls)
        dg = " | ".join(f"{float(cKDTree(g.atoms.xyz).query(xyz)[0].min()):.1f}" for g in sites)
        dp = "-" if pe is None else f"{float(pe.tree.query(xyz)[0].min()):.1f}"
        out.append(f"| {frame.label(h)} | {', '.join(ruls)} | {h - frame.offset} | {frame.arrays['res.rsa'][frame._index(h)]:.2f} | {lost} | {dg} | {dp} |")
    # the area a residue loses to the removed regions: how the cover warning depends on its threshold
    out += ["", "Cover check sensitivity: accepted designs with at least 2 contacted residues that lose more than the area to the removed regions", "",
            "| run | n | >10 A^2 (profile) | >20 | >30 | >50 | target residues above 10 / 30 A^2 |", "|---|---|---|---|---|---|---|"]
    for r in runs:
        cov = frame.covered_area(r)
        cells = [str(sum(sum(1 for u in x["tc"] if cov.get(u, 0) > a) >= 2 for x in acc[r])) for a in (10, 20, 30, 50)]
        out.append(f"| {r} | {len(acc[r])} | " + " | ".join(cells) +
                   f" | {sum(1 for a in cov.values() if a > 10)} / {sum(1 for a in cov.values() if a > 30)} of {frame.run_info(r)['n_resi']} |")
    out += ["", "Per-design table, accepted designs (distances in A; cut = contacted residues within the cut window; cover = contacted residues buried by the removed regions; "
            "membrane = signed distance past the plane, negative is clear)", "",
            "| design | run | rank | RMSD | partner min | partner n3 | rest min | rest n3 | glycan min (site) | cut | cover | membrane | partner | rest | cut | cover | overall |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    shout = {"pass": "pass", "warn": "warn", "fail": "FAIL"}
    for r in runs:
        for x in sorted(acc[r], key=lambda y: (y["rank"] is None, y["rank"] or 0)):
            m, v, g = x["metrics"], x["verdicts"], x["glycan_min"]
            gl = "-" if g["d"] is None else f"{g['d']:.1f} ({g['site']})"
            out.append(f"| {x['id']} | {r} | {x['rank']} | {x['place_rmsd']:.2f} | {m['olig_min']:.1f} | {m['olig_n3']} | {m['rest_min']:.1f} | {m['rest_n3']} | {gl} | "
                       f"{m['cut_n']} | {m['cov_n']} | {'-' if m['mem_excess'] is None else format(m['mem_excess'], '.0f')} | " + " | ".join(shout[v[k]["level"]] for k in ("partner", "rest", "cut", "cover", "overall")) + " |")
    fails = [x for r in runs for x in acc[r] if x["verdicts"]["overall"]["level"] == "fail"]
    out += ["", f"Accepted designs that fail the Context check ({len(fails)}): the other subunit or the rest of the chain", "",
            "| design | failing check | binder atoms within 3 A / 4 A | where | nearest rest residue |", "|---|---|---|---|---|"]
    for x in fails:
        failing = [k for k in ("partner", "rest") if x["verdicts"][k]["level"] == "fail"]
        where = ", ".join(f"{k} {v}" for k, v in (x.get("rest_by_region") or {}).items())
        pr = (x.get("rest_pair") or {}).get("rest")
        out.append(f"| {x['id']} | {', '.join(failing)} | {x['metrics']['rest_n3']} / {x['metrics']['rest_n4']} | {where} | "
                   f"{'' if not pr else pr[1] + str(pr[0])} |")
    return "\n".join(out) + "\n"


# ----------------------------------------------------------------------------------------------------------------
# command line: build the frame and the per-design table for BindCraft 1 runs
# ----------------------------------------------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(prog="python -m binder_review.target_frame", description=__doc__.split("\n")[0])
    ap.add_argument("--reference", required=True, help="reference structure of the full-length target (mmCIF or PDB)")
    ap.add_argument("--construct", required=True, help="the target file BindCraft was given")
    ap.add_argument("--run", action="append", default=[], metavar="ID=PATH", help="a BindCraft 1 run folder (repeat)")
    ap.add_argument("--uniprot", help="UniProt accession of the target (optional: without it no UniProt, SIFTS or AlphaFold record is read)")
    ap.add_argument("--reference-id", help="name of the reference structure in labels (default: the file stem)")
    ap.add_argument("--reference-chain", default="A", help="chain of the target in the reference (default A)")
    ap.add_argument("--partner-chains", default="auto",
                    help="chains of the other subunits: comma separated, none, or auto (the other chains with the same sequence)")
    ap.add_argument("--out", required=True, help="work directory (writes target/frame.* and target_checks.json)")
    ap.add_argument("--placement", choices=("interface", "whole", "core"), default="interface",
                    help="how designs are placed on the reference: interface (default: weighted around the binder), whole (every shared CA) or core (trimmed)")
    ap.add_argument("--no-checks", action="store_true", help="build and save the frame only")
    ap.add_argument("--report", metavar="FILE", help="also write the markdown tables of docs/target_frame.md to FILE")
    args = ap.parse_args(argv)
    from .readers import bindcraft1

    runs = {}
    for item in args.run:
        rid, _, path = item.partition("=")
        runs[rid] = bindcraft1.read_run(path, rid, construct_path=args.construct)
    partners = None if args.partner_chains == "auto" else () if args.partner_chains == "none" else tuple(args.partner_chains.split(","))
    frame = build_frame(args.reference, args.construct, runs, args.uniprot, chain=args.reference_chain, partner_chains=partners,
                        reference_id=args.reference_id, config={"geometry": {"fit": {"whole": "all"}.get(args.placement, args.placement)}})
    frame.assign(runs)
    out = Path(args.out)
    frame.save(out / "target")
    print(json.dumps(frame.summary(), indent=1, default=_json_default))
    if not args.no_checks:
        rows = []
        for rid, run in runs.items():
            rows += [c.to_dict() for c in check_run(frame, run)]
        doc = {"schema": "binderscope.target_checks/1", "generated": {"date": _dt.date.today().isoformat()},
               "config": frame.config, "n": len(rows), "designs": rows}
        (out / "target_checks.json").write_text(json.dumps(doc, indent=1, default=_json_default) + "\n")
        print(f"wrote {len(rows)} design rows to {out / 'target_checks.json'}")
        if args.report:
            Path(args.report).write_text(report_tables(frame, rows))
            print(f"wrote the report tables to {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
