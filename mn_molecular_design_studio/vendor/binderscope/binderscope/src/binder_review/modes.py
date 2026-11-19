"""Binding modes, the 2-D map embedding, hotspot coverage and near-miss rescue (task B2b).

Everything here is keyed by a design key (a string the caller chooses; `project_ids` gives the
blob ids) and works on UniProt-numbered contact sets, so it needs no knowledge of the target:
mode names are generated from the data and nothing is hard-coded to a target.

Pipeline (`build_modes`)
  1. footprint vectors (contacted target residues, UniProt numbers, 4.5 A heavy atom contacts) and
     contact centroids (mean CA of the contacted residues in one reference table);
  2. distance = 0.5 * Jaccard distance of the footprints + 1.0 * centroid distance (scaled by its
     maximum over the accepted designs), the mockup's recipe (real_export.py);
  3. Ward clustering of the accepted designs on classical-MDS coordinates of that distance, k chosen by
     silhouette among 3 to 7 under a balance rule (see `choose_k`); modes are ordered by size M1..Mk;
  4. rejected designs are assigned to the nearest mode and never influence clusters or the embedding;
  5. t-SNE of the accepted designs (fixed seed, MDS initialisation), scaled to 0..1; rejected designs
     are placed among their nearest accepted designs; classical MDS replaces t-SNE for tiny sets;
  6. hotspot coverage per design and per mode, auto names and one-sentence notes.

`rescue_records` is the near-miss part: failed filters, the numeric gap and the `nm` summary for
rejected designs (rule of docs/CONTRACTS.md 2.4, `model.near_miss_record`).
"""

from __future__ import annotations

import collections
import math
import os
import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Optional, Sequence

import numpy as np

from . import config, model, structure as S

MAX_MODES = 7  # the pages have seven colours and shapes
DEFAULT_SEED = 5  # t-SNE seed of the mockup
JACCARD_WEIGHT = 0.5
CENTROID_WEIGHT = 1.0
#: a hotspot counts as covered by a mode when this share of its designs touch it
HOT_SHARE_NAMED = 0.5
#: enrichment (share inside the mode minus share outside) needed to list a residue as enriched
ENRICH_MIN = 0.2
NAME_MAX_CHARS = 40
NEAR_HOTSPOT_MAX_A = 12.0


# ---------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------


@dataclass
class RunInfo:
    """What the mode builder needs to know about one run.

    hotspots   UniProt numbers of the run's hotspots (RunData.hotspots_uniprot()).
    designs    {design key: "accepted" | "rejected"} of the designs of this run that have a footprint.
    residues   UniProt numbers of the run's target (optional; only used to drop a footprint residue
               the run's target cannot contain, which would be a numbering error).
    """

    run_id: str
    hotspots: list[int] = field(default_factory=list)
    designs: dict[str, str] = field(default_factory=dict)
    residues: Optional[set[int]] = None
    label: str = ""


def project_ids(runs: Mapping[str, model.RunData]) -> dict[tuple[str, str], str]:
    """Blob id of every design: the id itself, prefixed with the run id ("run01/id") only when two
    runs would otherwise share it (data_contract.md, design.id)."""
    count = collections.Counter(d.id for r in runs.values() for d in r.designs)
    return {d.key: (d.id if count[d.id] == 1 else f"{d.run_id}/{d.id}") for r in runs.values() for d in r.designs}


def residue_names_from_runs(runs: Mapping[str, model.RunData]) -> dict[int, str]:
    """{UniProt number: one-letter code} over the targets of all runs (needs uniprot_offset)."""
    out: dict[int, str] = {}
    for r in runs.values():
        t = r.target
        for n, a in zip(t.numbers, t.sequence):
            out[t.to_uniprot(n)] = a
    return out


def reference_ca_from_construct(path: str | os.PathLike, uniprot_offset: int) -> dict[int, np.ndarray]:
    """{UniProt number: CA xyz} from a single-chain construct file (construct number + offset).
    Any complete CA table by UniProt number works, for instance the reference chain in the reference frame."""
    cx = S.load_complex(str(path), binder_chain=None)
    return {int(n) + uniprot_offset: np.asarray(x, float) for n, x in cx.target.ca_by_number().items()}


def design_footprint(structure_path: str, uniprot_offset: int, cutoff: float = model.CONTACT_CUTOFF,
                     n_binder_chains: int = 1) -> list[int]:
    """Contacted target residues of one design, UniProt numbers (heavy atoms within `cutoff`)."""
    cx = S.load_complex(str(structure_path), n_binder_chains=n_binder_chains)
    c = S.residue_contacts(cx.target, cx.binder, cutoff)
    return [int(n) + uniprot_offset for n in c.a_residues]


def collect_footprints(runs: Mapping[str, model.RunData], cutoff: float = model.CONTACT_CUTOFF) -> tuple[dict[str, list[int]], dict[str, RunInfo]]:
    """Footprints of every accepted design and every rejected design that has a structure, plus the
    RunInfo of each run. Runs need target.uniprot_offset. Costs about 20 ms per design."""
    ids = project_ids(runs)
    foot: dict[str, list[int]] = {}
    infos: dict[str, RunInfo] = {}
    for rid, run in runs.items():
        off = run.target.uniprot_offset
        if off is None:
            raise ValueError(f"run {rid}: target.uniprot_offset is not set (apply the target frame first)")
        info = RunInfo(rid, hotspots=run.hotspots_uniprot(), residues=set(run.target.uniprot_numbers()), label=run.label)
        for d in run.accepted() + run.rejected_with_structure():
            key = ids[d.key]
            foot[key] = design_footprint(d.structure_path, off, cutoff, model.binder_chain_count(d))
            info.designs[key] = d.status
        infos[rid] = info
    return foot, infos


def contact_centroids(footprints: Mapping[str, Iterable[int]], ref_ca: Mapping[int, Any]) -> dict[str, Optional[np.ndarray]]:
    """Mean CA of the contacted residues in the reference table; None for an empty footprint."""
    out: dict[str, Optional[np.ndarray]] = {}
    for k, fp in footprints.items():
        pts = [np.asarray(ref_ca[r], float) for r in fp if r in ref_ca]
        out[k] = np.mean(pts, axis=0) if pts else None
    return out


def default_sites(ref_ca: Mapping[int, Any], runs_info: Mapping[str, RunInfo]) -> dict[str, dict[str, Any]]:
    """Fallback `target.sites` when the target frame supplies none: A1 = the pooled hotspots."""
    hot = sorted({h for i in runs_info.values() for h in i.hotspots})
    return {"A1": {"residues": hot}} if hot else {}


# ---------------------------------------------------------------------------------------------
# Distance, clustering, embedding
# ---------------------------------------------------------------------------------------------


def footprint_matrix(footprints: Sequence[Iterable[int]], residues: Sequence[int]) -> np.ndarray:
    pos = {r: i for i, r in enumerate(residues)}
    m = np.zeros((len(footprints), len(residues)))
    for i, fp in enumerate(footprints):
        for r in fp:
            if r in pos:
                m[i, pos[r]] = 1.0
    return m


def jaccard_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """1 - |A and B| / |A or B| between the rows of a and the rows of b; two empty sets are 1 apart."""
    inter = a @ b.T
    union = a.sum(1)[:, None] + b.sum(1)[None, :] - inter
    return 1.0 - inter / np.maximum(union, 1.0)


def _euclid(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))


def classical_mds(d: np.ndarray, dims: int = 2) -> np.ndarray:
    """Classical (Torgerson) MDS coordinates of a distance matrix; columns with a non-positive
    eigenvalue are zero. Deterministic up to the sign of the eigenvectors, which is fixed here."""
    n = d.shape[0]
    if n == 0:
        return np.zeros((0, dims))
    j = np.eye(n) - 1.0 / n
    b = -0.5 * j @ (d ** 2) @ j
    w, v = np.linalg.eigh((b + b.T) / 2)
    idx = np.argsort(-w)[:dims]
    x = np.zeros((n, dims))
    for c, i in enumerate(idx):
        if w[i] > 1e-12:
            vec = v[:, i]
            vec = vec * (1.0 if vec[np.argmax(np.abs(vec))] > 0 else -1.0)
            x[:, c] = vec * math.sqrt(w[i])
    return x


def _silhouette(d: np.ndarray, labels: np.ndarray) -> float:
    from sklearn.metrics import silhouette_score

    if len(set(labels)) < 2 or len(set(labels)) >= len(labels):
        return float("nan")
    return float(silhouette_score(d, labels, metric="precomputed"))


def _ward_labels(coords: np.ndarray, k: int) -> np.ndarray:
    from scipy.cluster.hierarchy import fcluster, linkage

    z = linkage(coords, "ward")
    return fcluster(z, k, "maxclust")


@dataclass
class KChoice:
    k: int
    silhouette: Optional[float]
    rule: str
    table: list[dict[str, Any]]


def choose_k(d: np.ndarray, coords: np.ndarray, k_range: tuple[int, int] = (3, MAX_MODES), max_share: float = 0.4,
             min_size: int = 3) -> tuple[np.ndarray, KChoice]:
    """Ward labels for the k with the best silhouette among k_range under a balance rule.

    Silhouette alone favours k = 3 with one mode holding most designs, which does not triage anything,
    so a candidate must also have no cluster above `max_share` of the designs and none below `min_size`
    designs. If no k satisfies both, the size limit is dropped, then both; the rule that applied is
    reported. Small sets: k is capped at n - 1, and a set of fewer than 3 designs is one mode.
    """
    n = d.shape[0]
    if n < 3:
        return np.ones(n, dtype=int), KChoice(1, None, "fewer than 3 designs: one mode", [])
    lo, hi = k_range
    hi = max(2, min(hi, MAX_MODES, n - 1))
    lo = max(2, min(lo, hi))
    min_size = min(min_size, max(1, n // 8))
    rows = []
    for k in range(lo, hi + 1):
        lab = _ward_labels(coords, k)
        sizes = sorted(collections.Counter(lab).values(), reverse=True)
        rows.append({"k": k, "silhouette": round(_silhouette(d, lab), 4), "sizes": sizes, "labels": lab})

    def ok(r, share, size):
        return (not share or r["sizes"][0] / n <= max_share) and (not size or r["sizes"][-1] >= min_size)

    for share, size, rule in ((True, True, "best silhouette with no mode above %.0f%% of the designs and none under %d designs" % (100 * max_share, min_size)),
                              (False, True, "best silhouette with no mode under %d designs (no balanced k found)" % min_size),
                              (False, False, "best silhouette (no k met the size rules)")):
        cand = [r for r in rows if ok(r, share, size) and r["silhouette"] == r["silhouette"]]
        if cand:
            best = max(cand, key=lambda r: (r["silhouette"], -r["k"]))
            return best["labels"], KChoice(best["k"], best["silhouette"], rule, [{k: v for k, v in r.items() if k != "labels"} for r in rows])
    best = rows[0]
    return best["labels"], KChoice(best["k"], None, "no silhouette available", [{k: v for k, v in r.items() if k != "labels"} for r in rows])


def embed_2d(d: np.ndarray, seed: int = DEFAULT_SEED) -> tuple[np.ndarray, str]:
    """2-D coordinates in 0..1 for a distance matrix: t-SNE (precomputed metric, classical MDS start,
    fixed seed) from 10 designs on, classical MDS below that."""
    n = d.shape[0]
    if n == 0:
        return np.zeros((0, 2)), "none"
    if n < 3:
        pts = np.array([[0.0, 0.5], [1.0, 0.5]])[:n] if n == 2 else np.array([[0.5, 0.5]])
        return pts, "trivial"
    init = classical_mds(d, 2)
    method = "mds"
    xy = init
    if n >= 10:
        try:
            from sklearn.manifold import TSNE

            scale = init.std() or 1.0
            xy = TSNE(n_components=2, metric="precomputed", init=init / scale * 1e-4, perplexity=float(min(30.0, max(5.0, (n - 1) / 3.0))),
                      random_state=seed, learning_rate="auto").fit_transform(d)
            method = "tsne"
        except Exception:  # pragma: no cover - sklearn failure falls back to the MDS start
            xy, method = init, "mds"
    return _unit_square(xy), method


def _unit_square(xy: np.ndarray) -> np.ndarray:
    lo, hi = xy.min(0), xy.max(0)
    span = np.where(hi - lo > 1e-12, hi - lo, 1.0)
    out = (xy - lo) / span
    out[:, (hi - lo) <= 1e-12] = 0.5
    return out


# ---------------------------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------------------------


@dataclass
class ModeResult:
    """Output of build_modes.

    modes             campaign `modes` block (M1..Mk by size) with the extra keys hot_share, n_wrong_site,
                      mean_hot_frac, nearest_hot and centroid_dist (additional properties are allowed).
    assign            {design key: {mode, site, emb: [u, v], status, run, mode_dist}} for every design given.
    hotspots_touched  {design key: {run, hot, hot_of, frac, wrong_site}} (UniProt numbers; frac None
                      when the run has no hotspots).
    k_selection       how k was chosen (silhouette per k, rule).
    """

    modes: dict[str, dict[str, Any]]
    assign: dict[str, dict[str, Any]]
    hotspots_touched: dict[str, dict[str, Any]]
    k_selection: dict[str, Any]
    embedding: dict[str, Any]
    params: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {"modes": self.modes, "assign": self.assign, "hotspots_touched": self.hotspots_touched,
                "k_selection": self.k_selection, "embedding": self.embedding, "params": self.params}


# ---------------------------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------------------------


def _label(r: int, names: Mapping[int, str]) -> str:
    return f"{names.get(r, '')}{r}"


def _site_overlap(core: Iterable[int], residues: Iterable[int]) -> float:
    core, residues = set(core), set(residues)
    return len(core & residues) / len(core) if core else 0.0


def _pick_site(core: Iterable[int], centroid: Optional[np.ndarray], sites: Mapping[str, Mapping[str, Any]],
               site_centroids: Mapping[str, Optional[np.ndarray]]) -> Optional[str]:
    if not sites:
        return None
    keys = sorted(sites)
    ov = {k: _site_overlap(core, sites[k].get("residues", [])) for k in keys}
    best = max(ov.values())
    if best > 0:
        return next(k for k in keys if ov[k] == best)
    if centroid is not None:
        dist = {k: float(np.linalg.norm(centroid - site_centroids[k])) for k in keys if site_centroids.get(k) is not None}
        if dist:
            return min(dist, key=lambda k: (dist[k], k))
    return keys[0]


def _mode_name(enr: Sequence[int], all_hot: set[int], near: Optional[tuple[int, float]], names: Mapping[int, str], used: set[str]) -> str:
    """Generated name: up to three enriched residues ("S183/L276/K278 site"); when none of them is a
    hotspot and the nearest hotspot is close, "near <hotspot>" is added if it fits in 40 characters."""
    parts = [_label(r, names) for r in enr[:3]]
    base = ("/".join(parts) + " site") if parts else "Mixed contacts"
    name = base
    if parts and not (set(enr[:3]) & all_hot) and near is not None and near[1] <= NEAR_HOTSPOT_MAX_A:
        cand = f"{base} · near {_label(near[0], names)}"
        if len(cand) <= NAME_MAX_CHARS:
            name = cand
    if name in used:
        for extra in enr[3:]:
            cand = name.replace(" site", "/" + _label(extra, names) + " site", 1)
            if cand not in used and len(cand) <= NAME_MAX_CHARS:
                name = cand
                break
    i = 2
    unique = name
    while unique in used:
        unique = f"{name} ({i})"
        i += 1
    used.add(unique)
    return unique


def _pct(x: float) -> str:
    return f"{round(100 * x)}%"


def _join(items: Sequence[str]) -> str:
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


# ---------------------------------------------------------------------------------------------
# build_modes
# ---------------------------------------------------------------------------------------------


def build_modes(footprints: Mapping[str, Iterable[int]], centroids: Optional[Mapping[str, Any]], runs_info: Mapping[str, RunInfo], *,
                residue_names: Optional[Mapping[int, str]] = None, ref_ca: Optional[Mapping[int, Any]] = None,
                sites: Optional[Mapping[str, Mapping[str, Any]]] = None, k_range: tuple[int, int] = (3, MAX_MODES),
                max_share: float = 0.4, min_size: int = 3, seed: int = DEFAULT_SEED) -> ModeResult:
    """Cluster the accepted designs into binding modes and place every design.

    footprints  {design key: contacted target residues (UniProt)} for every accepted design and every
                rejected design that has a structure.
    centroids   {design key: xyz of the contact centroid} in any one frame; None computes them from
                `ref_ca` (CA by UniProt number). A design whose footprint is empty gets the mean centroid.
    runs_info   {run id: RunInfo} (hotspots, design keys with status).
    sites       target.sites (key -> {"residues": [...]}); None builds A1 from the pooled hotspots.
    The sites only name the site a design or mode overlaps most; they never move a design between modes.
    """
    names = dict(residue_names or {})
    foot = {k: sorted({int(r) for r in v}) for k, v in footprints.items()}
    status: dict[str, str] = {}
    run_of: dict[str, str] = {}
    for rid, info in runs_info.items():
        for k, st in info.designs.items():
            status[k], run_of[k] = st, rid
    missing = [k for k in foot if k not in status]
    if missing:
        raise ValueError(f"{len(missing)} footprints belong to no run in runs_info, e.g. {missing[0]}")
    for rid, info in runs_info.items():
        if info.residues is not None:
            bad = sorted({r for k, st in info.designs.items() for r in foot.get(k, []) if r not in info.residues})
            if bad:
                raise ValueError(f"run {rid}: footprint residues {bad[:5]} are not in the run's target (numbering error?)")

    acc = sorted(k for k in foot if status[k] == "accepted")
    rej = sorted(k for k in foot if status[k] != "accepted")
    if not acc:
        raise ValueError("no accepted design with a footprint")

    if centroids is None:
        if ref_ca is None:
            raise ValueError("pass centroids or ref_ca")
        cen_raw: dict[str, Any] = contact_centroids(foot, ref_ca)
    else:
        cen_raw = {k: (None if centroids.get(k) is None else np.asarray(centroids[k], float)) for k in foot}
    have = [c for c in cen_raw.values() if c is not None]
    fallback = np.mean(have, axis=0) if have else np.zeros(3)
    cen = {k: (cen_raw[k] if cen_raw[k] is not None else fallback) for k in foot}
    no_contacts = sorted(k for k in foot if not foot[k])

    residues = sorted({r for fp in foot.values() for r in fp} | {h for i in runs_info.values() for h in i.hotspots})
    fa = footprint_matrix([foot[k] for k in acc], residues)
    fr = footprint_matrix([foot[k] for k in rej], residues)
    ca_ = np.array([cen[k] for k in acc])
    cr_ = np.array([cen[k] for k in rej]) if rej else np.zeros((0, 3))

    dc_aa = _euclid(ca_, ca_)
    dc_scale = float(dc_aa.max()) or 1.0
    d_aa = JACCARD_WEIGHT * jaccard_distance(fa, fa) + CENTROID_WEIGHT * dc_aa / dc_scale
    np.fill_diagonal(d_aa, 0.0)
    d_ra = (JACCARD_WEIGHT * jaccard_distance(fr, fa) + CENTROID_WEIGHT * _euclid(cr_, ca_) / dc_scale) if rej else np.zeros((0, len(acc)))

    # -- clusters ----------------------------------------------------------------------------
    coords = classical_mds(d_aa, dims=min(10, max(2, len(acc) - 1)))
    labels, kc = choose_k(d_aa, coords, k_range, max_share, min_size)
    clusters = sorted(set(labels), key=lambda c: (-int((labels == c).sum()), float(np.mean([np.mean(foot[acc[i]] or [0]) for i in np.where(labels == c)[0]])), int(c)))
    if len(clusters) > MAX_MODES:  # not reachable with k <= 7; kept so the 7-colour rule is enforced in one place
        keep, rest = clusters[:MAX_MODES - 1], clusters[MAX_MODES - 1:]
        for c in rest:
            labels[labels == c] = rest[0]
        clusters = keep + [rest[0]]
    mode_of_cluster = {c: f"M{i + 1}" for i, c in enumerate(clusters)}
    acc_mode = {k: mode_of_cluster[labels[i]] for i, k in enumerate(acc)}

    # rejected: mean of the (up to three) smallest distances to each mode's members
    rej_mode: dict[str, str] = {}
    rej_dist: dict[str, float] = {}
    for j, k in enumerate(rej):
        best_m, best_v = None, None
        for c in clusters:
            idx = np.where(labels == c)[0]
            v = float(np.sort(d_ra[j, idx])[:3].mean())
            if best_v is None or v < best_v - 1e-12:
                best_m, best_v = mode_of_cluster[c], v
        rej_mode[k], rej_dist[k] = best_m, best_v

    # -- embedding ---------------------------------------------------------------------------
    emb_acc, emb_method = embed_2d(d_aa, seed)
    emb: dict[str, list[float]] = {k: [float(emb_acc[i, 0]), float(emb_acc[i, 1])] for i, k in enumerate(acc)}
    if rej:
        rng = np.random.default_rng(seed)
        kk = min(3, len(acc))
        nn_scale = np.sort(d_aa + np.eye(len(acc)) * 9, axis=1)[:, 0]
        typical = float(np.median(nn_scale)) if len(acc) > 1 else 0.1
        for j, k in enumerate(rej):
            idx = np.argsort(d_ra[j])[:kk]
            w = 1.0 / (d_ra[j, idx] + 0.05)
            p = (emb_acc[idx] * w[:, None]).sum(0) / w.sum()
            # jitter grows with the distance to the nearest accepted design so that coincident points stay readable
            amp = min(0.04, 0.01 + 0.1 * max(0.0, float(d_ra[j, idx[0]]) - typical))
            p = p + rng.normal(0, 1, 2) * amp
            emb[k] = [float(min(1.0, max(0.0, p[0]))), float(min(1.0, max(0.0, p[1])))]

    # -- hotspot coverage --------------------------------------------------------------------
    touched: dict[str, dict[str, Any]] = {}
    for k in foot:
        info = runs_info[run_of[k]]
        hs = list(info.hotspots)
        hot = sorted(set(foot[k]) & set(hs))
        touched[k] = {"run": run_of[k], "hot": hot, "hot_of": len(hs), "frac": (len(hot) / len(hs)) if hs else None,
                      "wrong_site": bool(hs) and not hot}

    # -- sites -------------------------------------------------------------------------------
    if sites is None:
        sites = default_sites(ref_ca or {}, runs_info)
    site_cen: dict[str, Optional[np.ndarray]] = {}
    for sk, s in sites.items():
        pts = [np.asarray(ref_ca[r], float) for r in s.get("residues", []) if ref_ca and r in ref_ca]
        site_cen[sk] = np.mean(pts, axis=0) if pts else None

    # -- per-mode summaries ------------------------------------------------------------------
    all_hot = {h for i in runs_info.values() for h in i.hotspots}
    modes: dict[str, dict[str, Any]] = {}
    used_names: set[str] = set()
    n_acc = len(acc)
    for c in clusters:
        mid = mode_of_cluster[c]
        m_idx = np.where(labels == c)[0]
        m_keys = [acc[i] for i in m_idx]
        freq_in = fa[m_idx].mean(0)
        freq_out = fa[labels != c].mean(0) if (labels != c).any() else np.zeros(len(residues))
        enrich = freq_in - freq_out
        order_enr = [i for i in np.lexsort((np.array(residues), -enrich)) if enrich[i] >= ENRICH_MIN][:6]
        enr_ranked = [residues[i] for i in order_enr]
        top_ranked = [residues[i] for i in np.lexsort((np.array(residues), -freq_in)) if freq_in[i] > 0][:8]
        core = [residues[i] for i in range(len(residues)) if freq_in[i] >= 0.5] or top_ranked[:4]
        mean_cen = np.mean([cen[k] for k in m_keys], axis=0)
        rej_keys = [k for k in rej if rej_mode[k] == mid]
        by_run = {rid: sum(1 for k in m_keys if run_of[k] == rid) for rid in runs_info if any(run_of[k] == rid for k in m_keys)}

        hot_share: dict[str, float] = {}
        for h in sorted(all_hot):
            app = [k for k in m_keys if h in runs_info[run_of[k]].hotspots]
            if app:
                hot_share[str(h)] = round(sum(1 for k in app if h in touched[k]["hot"]) / len(app), 3)
        fr_vals = [touched[k]["frac"] for k in m_keys if touched[k]["frac"] is not None]
        n_wrong = sum(1 for k in m_keys if touched[k]["wrong_site"])
        nearest_hot: dict[str, Any] = {}
        for rid, info in runs_info.items():
            if not ref_ca or not any(run_of[k] == rid for k in m_keys):
                continue
            dist = {h: float(np.linalg.norm(mean_cen - np.asarray(ref_ca[h], float))) for h in info.hotspots if h in ref_ca}
            if dist:
                h = min(dist, key=lambda x: (dist[x], x))
                nearest_hot[rid] = {"residue": h, "label": _label(h, names), "dist": round(dist[h], 1)}
        near_pooled = None
        if nearest_hot:
            v = min(nearest_hot.values(), key=lambda x: (x["dist"], x["residue"]))
            near_pooled = (v["residue"], v["dist"])

        name = _mode_name(sorted(enr_ranked[:3]) if enr_ranked else sorted(top_ranked[:3]), all_hot, near_pooled, names, used_names)
        site = _pick_site(core, mean_cen, sites, site_cen)
        modes[mid] = {
            "name": name, "note": "", "n": len(m_keys), "n_rejected": len(rej_keys), "by_run": by_run,
            "top": sorted(top_ranked), "enr": sorted(enr_ranked), "auto_name": True,
            "hot_share": hot_share, "n_wrong_site": n_wrong,
            "mean_hot_frac": round(float(np.mean(fr_vals)), 3) if fr_vals else None,
            "nearest_hot": nearest_hot,
        }
        if site is not None:
            modes[mid]["site"] = site
        modes[mid]["note"] = _mode_note(modes[mid], n_acc, names, enr_ranked)

    assign: dict[str, dict[str, Any]] = {}
    for k in foot:
        is_acc = status[k] == "accepted"
        md = acc_mode[k] if is_acc else rej_mode[k]
        core_k = foot[k]
        sk = _pick_site(core_k, cen[k] if k not in no_contacts else None, sites, site_cen)
        a = {"mode": md, "emb": [round(emb[k][0], 6), round(emb[k][1], 6)], "status": status[k], "run": run_of[k]}
        if sk is not None:
            a["site"] = sk
        if not is_acc:
            a["mode_dist"] = round(rej_dist[k], 4)
        if not foot[k]:
            a["no_contacts"] = True
        assign[k] = a

    return ModeResult(
        modes=modes, assign=assign, hotspots_touched=touched,
        k_selection={"k": kc.k, "silhouette": kc.silhouette, "rule": kc.rule, "table": kc.table, "method": "ward on classical MDS coordinates of the combined distance"},
        embedding={"method": emb_method, "seed": seed, "n_accepted": len(acc), "n_rejected": len(rej), "rejected_placement": "weighted mean of the 3 nearest accepted designs plus a small deterministic jitter"},
        params={"jaccard_weight": JACCARD_WEIGHT, "centroid_weight": CENTROID_WEIGHT, "centroid_scale": round(dc_scale, 3), "k_range": list(k_range),
                "max_share": max_share, "min_size": min_size, "contact_cutoff": model.CONTACT_CUTOFF, "no_contacts": no_contacts},
    )


def _mode_note(mode: Mapping[str, Any], n_total: int, names: Mapping[int, str], enr_ranked: Sequence[int]) -> str:
    """One to three sentences, built from the numbers only."""
    n = mode["n"]
    runs = ", ".join(f"{r} {c}" for r, c in mode["by_run"].items())
    s1 = f"{n} accepted design{'s' if n != 1 else ''} ({runs}), {_pct(n / n_total)} of the pool."
    hs = sorted(mode["hot_share"].items(), key=lambda kv: (-kv[1], int(kv[0])))
    parts = [f"{_label(int(h), names)} {_pct(v)}" for h, v in hs]
    s2 = ""
    if parts:
        s2 = " Touch rate per hotspot: " + _join(parts) + "."
        if mode["n_wrong_site"]:
            w = mode["n_wrong_site"]
            s2 += f" {w} design{'s' if w != 1 else ''} touch{'' if w != 1 else 'es'} none of {'their' if w != 1 else 'its'} run's hotspots."
    s3 = ""
    if enr_ranked:
        s3 = " Enriched against the other modes: " + _join([_label(r, names) for r in sorted(enr_ranked[:5])]) + "."
    if mode.get("site"):
        s3 += f" Overlaps site {mode['site']} most."
    if mode["n_rejected"]:
        s3 += f" {mode['n_rejected']} rejected design{'s' if mode['n_rejected'] != 1 else ''} {'falls' if mode['n_rejected'] == 1 else 'fall'} nearest to it."
    return (s1 + s2 + s3).strip()


# ---------------------------------------------------------------------------------------------
# End to end for RunData
# ---------------------------------------------------------------------------------------------


def build_from_runs(runs: Mapping[str, model.RunData], ref_ca: Mapping[int, Any], sites: Optional[Mapping[str, Mapping[str, Any]]] = None,
                    **kwargs: Any) -> tuple[ModeResult, dict[tuple[str, str], str]]:
    """Footprints, centroids and modes for RunData objects whose targets carry uniprot_offset.
    Returns the result and the {(run, design id): key} map used for the keys."""
    foot, infos = collect_footprints(runs)
    names = residue_names_from_runs(runs)
    res = build_modes(foot, None, infos, residue_names=names, ref_ca=ref_ca, sites=sites, **kwargs)
    return res, project_ids(runs)


# ---------------------------------------------------------------------------------------------
# Near-miss rescue
# ---------------------------------------------------------------------------------------------


def _round_fail(rec: dict[str, Any]) -> dict[str, Any]:
    for k in ("value", "threshold", "gap"):
        rec[k] = round(float(rec[k]), 6)
    return rec


def _bare_name(name: str) -> str:
    """A name from BindCraft 2's failed_filters without the "(not measured)" tag or an ".<state>" suffix."""
    return re.sub(r"\s*\(not measured\)\s*$", "", name).partition(".")[0].strip()


def design_failures(d: model.DesignRecord, rules: "Iterable[config.FilterRule] | config.FilterSet" = config.DEFAULT_FILTERS) -> list[config.FilterFailure]:
    """Filters a design fails under `rules`.

    Plain rules (and a FilterSet that is not a campaign's): BindCraft's Average and per-model checks against the metrics,
    nothing else, which is exactly what this function always did.

    A campaign FilterSet (the thresholds the run recorded itself): BindCraft's own verdict first. The reader recorded in
    `d.failed_filters` every failure BindCraft named, with a measured gap where the CSV shows the miss (scope "Average")
    and scope "Model", gap 0.0, where it does not (BindCraft judges confidence filters per validation model and the CSV
    holds the rounded mean). Those stand: where BindCraft said a design failed, it failed. A re-evaluation against the metrics
    then adds only a stem BindCraft did not name, because it stops at the first validation model that fails and so never
    names everything the design misses (an Average failure, measured). A stem is never listed twice.
    """
    fs = config.as_filter_set(rules)
    recomputed = config.check_filters(d.metrics, d.models, d.interface_aas, fs.rules)
    if not fs.from_campaign:
        return recomputed
    recorded = list(d.failed_filters)
    named = {f.stem for f in recorded}
    return recorded + [f for f in recomputed if f.stem not in named]


def unscored_failures(d: model.DesignRecord, failures: Sequence[config.FilterFailure]) -> list[str]:
    """Names BindCraft reported for a design that no failure represents (a threshold or a value the reader did not have,
    or a mandatory metric that was not measured). Order and spelling are BindCraft's. Only BindCraft 2 records names."""
    have = {f.stem.lower() for f in failures}
    return [n for n in d.extras.get("failed_filter_names", []) if _bare_name(n).lower() not in have]


def _fail_dict(f: config.FilterFailure, d: model.DesignRecord, fs: config.FilterSet) -> dict[str, Any]:
    """design.fail[] entry. For a campaign set it also says who found the failure and whether its size is known."""
    rec = _round_fail(model.failure_to_dict(f))
    if fs.from_campaign:
        reported = any(f.stem == r.stem for r in d.failed_filters)
        unmeasured = reported and f.scope != "Average" and f.gap == 0.0
        rec["reported"] = reported            # True: BindCraft named it. False: measured here from the CSV mean, BindCraft stopped before it
        rec["measured"] = not unmeasured      # False: BindCraft's verdict; the CSV cannot say by how much (gap is a lower bound, rel is null)
        if unmeasured:
            rec["rel"] = None
    return rec


def near_miss(d: model.DesignRecord, failures: Sequence[config.FilterFailure], fs: config.FilterSet) -> Optional[dict[str, Any]]:
    """design.nm of a rejected design under `fs`, None when it fails nothing.

    The rule is model.near_miss_record's: exactly one failed stem, missed by at most NEAR_MISS_MAX_REL of its threshold.
    For a campaign set it needs the miss to be a measured one. A failure BindCraft reported but whose size the CSV does not
    show (scope "Model", gap 0.0) is a lower bound, not a measured gap, and reading it as a miss of 0 percent would put
    it at the top of the rescue list; a name BindCraft reported that no failure represents has no size at all. Either way
    the design is not called near (`near` False, `score` 0, the names in `unmeasured`).
    """
    nm = model.near_miss_record(failures)
    if not fs.from_campaign:
        return nm
    unmeasured = sorted({f.stem for f in failures if _fail_dict(f, d, fs)["measured"] is False})
    unscored = unscored_failures(d, failures)
    if not unmeasured and not unscored:
        return nm
    extra = sorted({_bare_name(n) for n in unscored} - {""})
    stems = sorted(set(nm["stems"] if nm else []) | set(extra))
    return {"n_failed": len(stems), "stems": stems, "worst_rel": nm["worst_rel"] if nm else 0.0, "score": 0.0, "near": False,
            "unmeasured": sorted(set(unmeasured) | set(extra))}


def _filter_set(run: model.RunData, rules: "Optional[Iterable[config.FilterRule] | config.FilterSet]") -> config.FilterSet:
    """`rules` when the caller gave some, else the filters the run was judged by (config.filters_for_run)."""
    return config.filters_for_run(run) if rules is None else config.as_filter_set(rules)


def rescue_records(run: model.RunData, rules: "Optional[Iterable[config.FilterRule] | config.FilterSet]" = None, *,
                   include_metric_only: bool = False) -> dict[str, dict[str, Any]]:
    """{design id: {"fail": [...], "nm": {...}}} for the rejected designs of a run.

    `rules` default to the filters the run was judged by: the thresholds it recorded itself (BindCraft 2) or, when it recorded
    none, BindCraft 1's reconstructed defaults. Pass rules to judge by others. By default only designs that have a structure
    are listed, because only those enter designs[]. include_metric_only adds the rows without a structure for the
    Query statistics. Ids are the run-local design ids. A design that fails nothing under `rules`
    (BindCraft rejected it for a reason the rules do not cover) has fail [] and nm None.
    """
    fs = _filter_set(run, rules)
    out: dict[str, dict[str, Any]] = {}
    pool = run.rejected() if include_metric_only else run.rejected_with_structure()
    for d in pool:
        fl = design_failures(d, fs)
        out[d.id] = {"fail": [_fail_dict(f, d, fs) for f in fl], "nm": near_miss(d, fl, fs)}
        if out[d.id]["nm"]:
            out[d.id]["nm"]["worst_rel"] = round(out[d.id]["nm"]["worst_rel"], 3)
    return out


def rescue_order(records: Mapping[str, Mapping[str, Any]]) -> list[str]:
    """Design ids from the closest miss to the furthest: near-misses first, then by nm.score, then id."""
    def key(i: str):
        nm = records[i].get("nm")
        if not nm:
            return (2, 0.0, i)
        return (0 if nm["near"] else 1, -nm["score"], i)
    return sorted(records, key=key)


def accepted_failures(run: model.RunData, rules: "Optional[Iterable[config.FilterRule] | config.FilterSet]" = None) -> dict[str, list[dict[str, Any]]]:
    """Accepted designs that fail a rule (should be empty; non-empty means the filter set is wrong).

    Under the run's own recorded thresholds this is a real check: BindCraft passed these designs, so a failure here means a
    threshold was misread or mapped to the wrong metric. Under BindCraft 1's reconstructed defaults it is the check
    that the reconstruction is right, and it is what showed that a BindCraft 2 campaign was being judged by a rule it
    never applied (Surface_Hydrophobicity)."""
    fs = _filter_set(run, rules)
    out = {}
    for d in run.accepted():
        fl = design_failures(d, fs)
        if fl:
            out[d.id] = [_fail_dict(f, d, fs) for f in fl]
    return out


def rescue_summary(run: model.RunData, rules: "Optional[Iterable[config.FilterRule] | config.FilterSet]" = None) -> dict[str, Any]:
    """Counts for the doc and for qualitative checks against failure_csv.

    per_stem[stem] is the number of not-accepted rows that fail the stem under `rules` (a row counts once
    per stem), next to run.failure_counters[stem]. They agree for the counter-verified rules; for
    rules BindCraft applies before a sequence gets a row (pLDDT, i_pTM, i_pAE) the rows hold no
    failure and the counter is far larger. `by_n_failed` is over the rows with a structure. For a run that recorded its own
    thresholds the stems are the campaign's filter names, the same names as the counters, and the summary adds `filters`
    (where the set came from and what could not be evaluated).
    """
    fs = _filter_set(run, rules)
    rows = run.rejected()
    stem_counts: collections.Counter = collections.Counter()
    nfail_struct: collections.Counter = collections.Counter()
    n_near = 0
    n_none = 0
    for d in rows:
        fl = design_failures(d, fs)
        stems = list(dict.fromkeys(f.stem for f in fl))
        stem_counts.update(stems)
        if d.has_structure:
            nfail_struct[len(stems)] += 1
            nm = near_miss(d, fl, fs)
            n_near += 1 if nm and nm["near"] else 0
        if not fl and not (fs.from_campaign and unscored_failures(d, fl)):
            n_none += 1
    per_stem = []
    for r in fs.rules:
        per_stem.append({"stem": r.stem, "provenance": r.provenance, "rows_failing": stem_counts.get(r.stem, 0),
                         "failure_csv": run.failure_counters.get(r.stem)})
    out = {"run": run.run_id, "n_rejected_rows": len(rows), "n_with_structure": sum(1 for d in rows if d.has_structure),
           "n_near": n_near, "by_n_failed": dict(sorted(nfail_struct.items())), "per_stem": per_stem, "no_failure": n_none}
    if fs.from_campaign:
        out["filters"] = {"source": fs.source, "unevaluated": [u.name for u in fs.unevaluated], "unlimited": list(fs.unlimited)}
    return out
