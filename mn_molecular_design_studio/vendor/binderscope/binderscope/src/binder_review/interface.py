"""Recomputed interface analysis of one binder-target complex, and its per-residue tracks.

``analyze_design`` takes a complex loaded through :mod:`binder_review.structure` and returns a
:class:`DesignAnalysis`: contacts, hydrogen bonds, salt bridges, buried unsatisfied polar atoms, buried area,
hydrophobic patches, secondary structure, termini, hotspots touched, per-residue tracks and a comparison with
BindCraft's own columns. ``analyze_run`` does this for many designs with a process pool and a disk cache.

Conventions (docs/CONTRACTS.md): target residue numbers are CONSTRUCT numbers (the numbering of the file BindCraft
was given, gaps allowed: a run may have deleted loops); binder residues are 1-based; analysis happens
in the design's own frame (geometry does not depend on the frame). A contact is a heavy-atom pair within 4.5 A
(``model.CONTACT_CUTOFF``). BindCraft's own interface rule (4.0 A, hydrogens counted) is only used in
``bindcraft_bc`` and ``vs_bindcraft`` to compare with its columns. Nothing here needs ipSAE, a PAE matrix, a
second predictor or Rosetta energies; those quantities are not computed.
"""
from __future__ import annotations

import concurrent.futures as cf
import dataclasses
import hashlib
import json
import multiprocessing as mp
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Optional, Sequence

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import disulfides as D
from . import liabilities as L
from . import structure as S
from .model import (BINDCRAFT_INTERFACE_CUTOFF, CONTACT_CUTOFF, HBOND_MAX_DIST, HOTSPOT_TOUCH_CUTOFF,
                    SALT_BRIDGE_MAX_DIST, DesignRecord, binder_chain_count)

ANALYSIS_VERSION = 2          # bump when the numbers change; part of every cache key
DEFAULT_N_POINTS = "accurate"

CLOSE_CUTOFF = 3.5            # the closer contact set
# Buried unsatisfied polar atoms: a polar heavy atom counts when binding buries it (its SASA drops by at
# least BUNS_MIN_BURIAL), it is closed off in the complex (no accessible surface: SASA at most BUNS_MAX_SASA with
# the 1.4 A probe and 960 sphere points) and it has no hydrogen bond in the complex. Rosetta's BuriedUnsatHbonds
# (BindCraft) has the same idea with a 1.1 A probe and a delta between complex and separated chains; see
# docs/interface_validation.md for the sensitivity.
BUNS_MIN_BURIAL = 0.5         # A^2
BUNS_MAX_SASA = 0.1           # A^2
PATCH_MIN_ATOM_SASA = 1.0     # A^2, exposed hydrophobic atom
PATCH_MAX_INTERFACE = 0.5     # A^2, an atom with more burial than this belongs to the interface, not a patch
PATCH_LINK = 5.0              # A, hydrophobic atoms closer than this are one patch
MAX_PATCHES = 3
HYDROPHOBIC_AA = frozenset("ACFILMPVWY")   # BindCraft's hydrophobic set (Interface_Hydrophobicity)
PATCH_AA = frozenset("AVLIMFWYC")          # the mockup's patch residues
AROMATIC_AA = frozenset("FWY")
POSITIVE_AA = frozenset("KR")
NEGATIVE_AA = frozenset("DE")
POLAR_AA = frozenset("STNQH")
NEAR_LOOP_DIST = 6.0          # A, a loop closer than this to the target counts as near the interface
SS_BINDCRAFT = dict(helix="HGI", strand="E", fill_bulges=False)   # BindCraft counts H, G, I as helix and E as sheet

_PATCH_BACKBONE = frozenset(("C", "CA", "N", "O"))


# --------------------------------------------------------------------------------------------- result object
@dataclass
class DesignAnalysis:
    """Everything ``analyze_design`` computes. JSON-serialisable through :meth:`to_dict`.

    Target residue numbers are construct numbers, binder residue numbers are 1-based indexes into the sequence
    (``bc`` is 1-based here; the campaign blob's ``design.bc`` is 0-based and the exporter converts). Per-residue
    lists follow the binder (``*_res``, ``ss``) in residue order. ``None`` means "not computable".
    """

    id: Optional[str] = None
    run_id: Optional[str] = None
    n: int = 0                                  # binder length
    n_target: int = 0                           # target residues present
    # contacts
    tc: list = field(default_factory=list)      # target residues in contact (4.5 A heavy atoms)
    bc: list = field(default_factory=list)      # binder residues in contact (4.5 A), 1-based
    nc: int = 0                                 # atom pairs within 4.5 A
    close_tc: list = field(default_factory=list)  # same at 3.5 A
    close_bc: list = field(default_factory=list)
    bindcraft_bc: list = field(default_factory=list)  # BindCraft's rule: 4.0 A, hydrogens counted (for validation)
    bindcraft_bc_heavy: list = field(default_factory=list)  # the same with heavy atoms only (the unrelaxed prediction has no H)
    contacts: list = field(default_factory=list)      # [{t, b, d, n}] residue pairs at 4.5 A
    contact_map: dict = field(default_factory=dict)   # {t: [...], b: [...], d: [[...]]} min distances, None = no contact
    # polar interactions
    hbonds: list = field(default_factory=list)  # [{t, tn, ta, b, bn, ba, d, k}]  k: D>A target donates, A<D target accepts
    salt: list = field(default_factory=list)    # [{t, tn, b, bn, d}] one per residue pair
    buns: list = field(default_factory=list)    # [{side, resi, resn, atom}] side t (construct number) or b (1-based)
    n_hbonds: int = 0
    n_salt: int = 0
    n_unsat: int = 0
    n_unsat_binder: int = 0
    n_unsat_target: int = 0
    unsat_density: float = 0.0                  # unsatisfied atoms per 1000 A^2 of buried area
    hbond_basis: str = "explicit"               # hydrogens from the file on both chains ('explicit'), one chain ('partial') or none ('computed': amide H placed, side chains by heavy-atom geometry)
    # surface
    dsasa: float = 0.0                          # buried area, both sides
    dsasa_binder: float = 0.0
    dsasa_target_total: float = 0.0
    dsasa_res: list = field(default_factory=list)   # binder side, per residue
    t_dsasa: dict = field(default_factory=dict)     # target side, {construct number: A^2}
    binder_sasa_alone: float = 0.0
    hyd_frac: float = 0.0                       # share of the buried area that is carbon or sulfur
    hphob_patch: float = 0.0                    # largest exposed hydrophobic patch, A^2
    patches: list = field(default_factory=list)  # at most 3, [{area, res}]
    comp: dict = field(default_factory=dict)    # interface composition, binder and target sides
    # secondary structure
    ss: str = ""                                # H, E, L per binder residue (DSSP of the complex)
    loop_frac: float = 0.0                      # share of interface residues in loops
    helix_frac: float = 0.0                     # share of all binder residues
    strand_frac: float = 0.0
    iface_helix: float = 0.0                    # share of interface residues in helices
    iface_strand: float = 0.0
    segments: list = field(default_factory=list)       # [{ss, start, end, n_iface, min_d}]
    loops_near: list = field(default_factory=list)     # loop segments within NEAR_LOOP_DIST of the target
    # termini
    nterm_d: float = 0.0                        # N-terminal CA to the nearest target heavy atom
    cterm_d: float = 0.0
    termini: dict = field(default_factory=dict)  # {n: {d, contact}, c: {d, contact}}
    center: list = field(default_factory=list)  # centre of the interface, design frame
    # per-residue
    rsa_res: list = field(default_factory=list)     # binder relative SASA in the complex
    plddt_res: list = field(default_factory=list)   # B-factor, 0 to 100
    plddt_t: dict = field(default_factory=dict)     # {construct number: pLDDT 0 to 100}
    tracks: dict = field(default_factory=dict)
    ssbonds: list = field(default_factory=list)     # [[i, j]] the binder's disulfides, 1-based
    ssbond_basis: str = "none"                      # "recorded" (the file lists them), "inferred" (SG-SG geometry) or "none"
    disulfides: list = field(default_factory=list)  # [{a, b, d (A), band ("recorded", "close" or "loose")}]
    # hotspots (the run's own)
    hotspots: list = field(default_factory=list)    # hotspot construct numbers present in the target
    hotspot_dist: dict = field(default_factory=dict)  # {construct number: min heavy-atom distance to the binder}
    hot: list = field(default_factory=list)         # hotspots touched: tc and hotspots
    # sequence
    seq_props: dict = field(default_factory=dict)   # pI, charge, gravy, n_cys, n_met, n_trp, n_lys, length
    liab: list = field(default_factory=list)        # liabilities.find_liabilities
    global_flags: list = field(default_factory=list)
    interface_aas: dict = field(default_factory=dict)  # amino-acid counts of bindcraft_bc
    # comparison and bookkeeping
    vs_bindcraft: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return _jsonable(dataclasses.asdict(self))

    @classmethod
    def from_dict(cls, d: Mapping) -> "DesignAnalysis":
        d = dict(d)
        for key in ("t_dsasa", "plddt_t", "hotspot_dist"):
            if key in d and d[key] is not None:
                d[key] = {int(k): v for k, v in d[key].items()}
        names = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in names})


def _jsonable(o: Any) -> Any:
    if isinstance(o, dict):
        return {(int(k) if isinstance(k, (np.integer,)) else k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return None if f != f or f in (float("inf"), float("-inf")) else f
    if isinstance(o, np.bool_):
        return bool(o)
    return o


# --------------------------------------------------------------------------------------------- loading
def load_for_analysis(path: str | os.PathLike, target_chain: str = "A", binder_chain: str = "B",
                      n_binder_chains: int | None = None) -> S.Complex:
    """Load a BindCraft complex for analysis: hydrogens kept (H-bonds use them), heavy-atom tables otherwise."""
    return S.load_complex(str(path), target_chain=target_chain, binder_chain=binder_chain, strip_h=False,
                          n_binder_chains=n_binder_chains)


# --------------------------------------------------------------------------------------------- helpers
def _frac(count: float, total: float) -> float:
    return float(count) / float(total) if total else 0.0


def _runs(s: str) -> list[tuple[str, int, int]]:
    """Runs of equal letters as (letter, first, last) with 1-based inclusive positions."""
    out, i = [], 0
    while i < len(s):
        j = i
        while j < len(s) and s[j] == s[i]:
            j += 1
        out.append((s[i], i + 1, j))
        i = j
    return out


def composition(letters: Sequence[str]) -> dict:
    """Counts and shares of amino-acid classes in a list of one-letter codes."""
    n = len(letters)
    sets = {"hydrophobic": HYDROPHOBIC_AA, "aromatic": AROMATIC_AA, "positive": POSITIVE_AA,
            "negative": NEGATIVE_AA, "polar": POLAR_AA}
    out: dict = {"n": n, "aa": {a: 0 for a in S.AA_ORDER}}
    for c in letters:
        if c in out["aa"]:
            out["aa"][c] += 1
    for name, grp in sets.items():
        k = sum(1 for c in letters if c in grp)
        out[name] = k
        out[f"{name}_frac"] = _frac(k, n)
    k = sum(1 for c in letters if c in POSITIVE_AA or c in NEGATIVE_AA)
    out["charged"], out["charged_frac"] = k, _frac(k, n)
    return out


def _dedupe_hbonds(hb: "S.HBonds") -> np.ndarray:
    """Row indices of one entry per atom pair (a Ser/Thr/Tyr pair that passes both ways is listed twice)."""
    if not len(hb):
        return np.zeros(0, dtype=np.int64)
    order = np.lexsort((np.where(np.isnan(hb.h_distance), 99.0, hb.h_distance), hb.b_index, hb.a_index))
    seen, keep = set(), []
    for r in order:
        key = (int(hb.a_index[r]), int(hb.b_index[r]))
        if key not in seen:
            seen.add(key)
            keep.append(int(r))
    return np.array(sorted(keep, key=lambda r: (hb.a_index[r], hb.b_index[r])), dtype=np.int64)


def _satisfied_atoms(comb_all: "S.Atoms", heavy_pos: np.ndarray) -> np.ndarray:
    """Boolean per heavy atom of the complex: takes part in at least one hydrogen bond."""
    hb = S.hbonds(comb_all, max_dist=HBOND_MAX_DIST)
    sat = np.zeros(int((heavy_pos >= 0).sum()), dtype=bool)
    for idx in np.concatenate([hb.a_index, hb.b_index]) if len(hb) else []:
        p = heavy_pos[int(idx)]
        if p >= 0:
            sat[p] = True
    return sat


def _label_components(xyz: np.ndarray, link: float) -> np.ndarray:
    n = len(xyz)
    if n == 0:
        return np.zeros(0, dtype=int)
    pairs = cKDTree(xyz).query_pairs(link, output_type="ndarray")
    if not len(pairs):
        return np.arange(n)
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(g, directed=False)[1]


# --------------------------------------------------------------------------------------------- main
def analyze_design(cx: S.Complex, design: Optional[DesignRecord] = None, hotspots_construct: Iterable[int] = (), *,
                   n_points: int | str = DEFAULT_N_POINTS, buns_min_burial: float = BUNS_MIN_BURIAL,
                   buns_max_sasa: float = BUNS_MAX_SASA, run_id: Optional[str] = None,
                   design_id: Optional[str] = None) -> DesignAnalysis:
    """Analyse one complex.

    cx                  from :func:`load_for_analysis` (hydrogens kept). A complex loaded with ``strip_h=True``
                        is reloaded from ``cx.path`` so the H-bond rules can use the hydrogens in the file.
    design              the design's :class:`DesignRecord`: its BindCraft columns fill ``vs_bindcraft`` and its
                        id and run are copied. Optional.
    hotspots_construct  the hotspots of the design's own run, construct numbers (e.g. 26, 122, 125).
                        Hotspots absent from the design's target (a deleted loop) are
                        listed in ``meta.hotspots_missing`` and never count as touched.
    n_points            SASA sphere points ('accurate' 960, the default; 'default' is 200 and 0.05 percent off in
                        dSASA but too coarse for the burial test of the unsatisfied polar atoms).

    Returns a :class:`DesignAnalysis`. Target numbers are construct numbers; the exporter converts them to
    UniProt numbers through ``TargetConstruct.to_uniprot``.
    """
    notes: list[str] = []
    if cx.hydrogens_stripped and cx.has_hydrogens and cx.path:
        # the same chains as before: reloading with the default would drop all but the first binder chain
        cx = S.load_complex(cx.path, target_chain=cx.target_chain,
                            binder_chain=list(cx.binder_chains) or cx.binder_chain, strip_h=False)
    tgt_all, bnd_all = cx.target, cx.binder          # may carry hydrogens
    tgt, bnd = tgt_all.heavy(), bnd_all.heavy()
    if not len(tgt) or not len(bnd):
        raise S.StructureError(f"{cx.path}: target or binder has no heavy atoms")
    t_h, b_h = bool(tgt_all.is_hydrogen.any()), bool(bnd_all.is_hydrogen.any())
    has_h = t_h and b_h
    tr, br = tgt.residues, bnd.residues
    n, n_t = len(br), len(tr)
    bnum = br.number
    if not np.array_equal(bnum, np.arange(1, n + 1)):
        notes.append("binder residue numbers are not 1..n; per-binder values are indexed by residue order (1-based)")
    b_of = {int(num): i + 1 for i, num in enumerate(bnum)}          # file number -> 1-based index
    seq = br.sequence
    tletters = list(tr.one)

    a = DesignAnalysis(id=design_id or (design.id if design else None), run_id=run_id or (design.run_id if design else None),
                       n=n, n_target=n_t)

    # ---- contacts (page definition, 4.5 A heavy atoms) and the closer set
    c45 = S.residue_contacts(tgt, bnd, CONTACT_CUTOFF)
    c35 = S.residue_contacts(tgt, bnd, CLOSE_CUTOFF)
    a.tc = c45.a_residues
    a.bc = sorted(b_of[x] for x in c45.b_residues)
    a.nc = int(c45.n_atom_contacts)
    a.close_tc = c35.a_residues
    a.close_bc = sorted(b_of[x] for x in c35.b_residues)
    a.bindcraft_bc = sorted(b_of[x] for x in S.bindcraft_interface_residues(cx, BINDCRAFT_INTERFACE_CUTOFF, True))
    a.bindcraft_bc_heavy = sorted(b_of[x] for x in S.bindcraft_interface_residues(cx, BINDCRAFT_INTERFACE_CUTOFF, False))
    a.contacts = [{"t": int(t), "b": b_of[int(b)], "d": round(float(d), 2), "n": int(k)}
                  for t, b, d, k in zip(c45.num_a, c45.num_b, c45.min_dist, c45.n_atom_pairs)]
    a.contacts.sort(key=lambda r: (r["t"], r["b"]))
    ti = {t: i for i, t in enumerate(a.tc)}
    bi = {b: i for i, b in enumerate(a.bc)}
    dm = [[None] * len(a.tc) for _ in a.bc]
    for r in a.contacts:
        dm[bi[r["b"]]][ti[r["t"]]] = round(r["d"], 1)
    a.contact_map = {"t": a.tc, "b": a.bc, "d": dm}
    bc_set = set(a.bc)

    # ---- buried area
    ids = S.interface_dsasa(cx, n_points=n_points)
    a.dsasa = float(ids.dsasa_total)
    a.dsasa_binder, a.dsasa_target_total = float(ids.dsasa_binder), float(ids.dsasa_target)
    a.binder_sasa_alone = float(ids.binder_sasa_alone.sum())
    a.dsasa_res = [float(x) for x in ids.binder_residue]
    a.t_dsasa = {int(num): float(v) for num, v in zip(tr.number, ids.target_residue)}
    b_res_sasa = bnd.per_residue(ids.binder_sasa_complex)
    a.rsa_res = [float(x) if x == x else None for x in S.relative_sasa(b_res_sasa, br.one)]
    dsa_c = np.concatenate([ids.target_atom, ids.binder_atom])           # burial per heavy atom, target then binder
    sasa_c = np.concatenate([ids.target_sasa_complex, ids.binder_sasa_complex])
    el_all = np.concatenate([tgt.element, bnd.element])
    carbon = np.isin(el_all, ("C", "S"))
    a.hyd_frac = _frac(dsa_c[carbon].sum(), a.dsasa) if a.dsasa > 0 else 0.0

    # ---- secondary structure (DSSP of the complex, as BindCraft does, taken on the binder)
    d8 = S.assign_dssp(S.Atoms.concat([tgt, bnd]))[n_t:]
    a.ss = S.ss3_from_dssp(d8, **SS_BINDCRAFT)
    ssl = list(a.ss)
    a.helix_frac = _frac(ssl.count("H"), n)
    a.strand_frac = _frac(ssl.count("E"), n)
    ifs = [ssl[i - 1] for i in a.bc]
    a.loop_frac = _frac(ifs.count("L"), len(ifs))
    a.iface_helix = _frac(ifs.count("H"), len(ifs))
    a.iface_strand = _frac(ifs.count("E"), len(ifs))

    # ---- distances of binder residues to the target
    ttree = cKDTree(tgt.xyz)
    d_atom, _ = ttree.query(bnd.xyz)
    res_min_d = np.minimum.reduceat(d_atom, br.start)
    ca_i = br.ca
    nterm_d = float(ttree.query(bnd.xyz[ca_i[0]])[0]) if ca_i[0] >= 0 else float(res_min_d[0])
    cterm_d = float(ttree.query(bnd.xyz[ca_i[-1]])[0]) if ca_i[-1] >= 0 else float(res_min_d[-1])
    a.nterm_d, a.cterm_d = nterm_d, cterm_d
    a.termini = {"n": {"d": nterm_d, "contact": 1 in bc_set}, "c": {"d": cterm_d, "contact": n in bc_set}}
    for ss_letter, first, last in _runs(a.ss):
        seg = {"ss": ss_letter, "start": first, "end": last,
               "n_iface": int(sum(1 for r in range(first, last + 1) if r in bc_set)),
               "min_d": float(res_min_d[first - 1:last].min())}
        a.segments.append(seg)
        if ss_letter == "L" and seg["min_d"] <= NEAR_LOOP_DIST:
            a.loops_near.append(seg)

    # ---- interface centre: mean midpoint of the atom pairs in contact
    ai, bj, dd = S.atom_contacts(tgt, bnd, CONTACT_CUTOFF, heavy_only=True)
    if len(ai):
        a.center = [float(x) for x in (0.5 * (tgt.xyz[ai] + bnd.xyz[bj])).mean(0)]
    else:
        a.center = [float(x) for x in (0.5 * (tgt.xyz.mean(0) + bnd.xyz.mean(0)))]

    # ---- hydrogen bonds across the interface (explicit hydrogens when present)
    hb = S.hbonds(tgt_all, bnd_all, max_dist=HBOND_MAX_DIST)
    keep = _dedupe_hbonds(hb)
    a.hbond_basis = "explicit" if has_h else ("partial" if t_h or b_h else "computed")
    for r in keep:
        ia_, ib_ = int(hb.a_index[r]), int(hb.b_index[r])
        a.hbonds.append({"t": int(tgt_all.resnum[ia_]), "tn": str(tgt_all.resname[ia_]), "ta": str(tgt_all.name[ia_]),
                         "b": b_of[int(bnd_all.resnum[ib_])], "bn": str(bnd_all.resname[ib_]), "ba": str(bnd_all.name[ib_]),
                         "d": round(float(hb.distance[r]), 2), "k": "D>A" if hb.donor_is_a[r] else "A<D"})
    a.n_hbonds = len(a.hbonds)

    # ---- salt bridges, one per residue pair
    sb: dict[tuple[int, int], dict] = {}
    for s_ in S.salt_bridges(tgt, bnd, cutoff=SALT_BRIDGE_MAX_DIST):
        key = (int(s_.a_resnum), b_of[int(s_.b_resnum)])
        if key not in sb or s_.distance < sb[key]["d"]:
            sb[key] = {"t": key[0], "tn": s_.a_resname, "b": key[1], "bn": s_.b_resname, "d": round(float(s_.distance), 2)}
    a.salt = [sb[k] for k in sorted(sb)]
    a.n_salt = len(a.salt)

    # ---- buried unsatisfied polar atoms
    comb_all = S.Atoms.concat([tgt_all, bnd_all])
    heavy_mask = ~comb_all.is_hydrogen
    heavy_pos = np.full(len(comb_all), -1, dtype=np.int64)
    heavy_pos[heavy_mask] = np.arange(int(heavy_mask.sum()))
    sat = _satisfied_atoms(comb_all, heavy_pos)
    comb_heavy = comb_all.select(heavy_mask)
    don, acc, _ = S._polar_roles(comb_heavy)
    polar = don | acc
    cand = polar & (dsa_c >= buns_min_burial) & (sasa_c <= buns_max_sasa) & ~sat
    n_tgt_heavy = len(tgt)
    for i in np.nonzero(cand)[0]:
        if i < n_tgt_heavy:
            a.buns.append({"side": "t", "resi": int(tgt.resnum[i]), "resn": str(tgt.resname[i]), "atom": str(tgt.name[i])})
        else:
            j = i - n_tgt_heavy
            a.buns.append({"side": "b", "resi": b_of[int(bnd.resnum[j])], "resn": str(bnd.resname[j]), "atom": str(bnd.name[j])})
    a.n_unsat = len(a.buns)
    a.n_unsat_binder = sum(1 for b_ in a.buns if b_["side"] == "b")
    a.n_unsat_target = a.n_unsat - a.n_unsat_binder
    a.unsat_density = a.n_unsat / (a.dsasa / 1000.0) if a.dsasa > 0 else 0.0

    # ---- exposed hydrophobic patches on the binder
    b_sasa_c, b_d = ids.binder_sasa_complex, ids.binder_atom
    side_chain = ~np.isin(bnd.name, list(_PATCH_BACKBONE))
    hyd_atoms = np.nonzero(np.isin(bnd.element, ("C", "S")) & np.isin(bnd.resname, [S.ONE_TO_THREE[c] for c in PATCH_AA])
                           & side_chain & (b_sasa_c > PATCH_MIN_ATOM_SASA) & (b_d < PATCH_MAX_INTERFACE))[0]
    patches = []
    if len(hyd_atoms):
        lab = _label_components(bnd.xyz[hyd_atoms], PATCH_LINK)
        for c in np.unique(lab):
            members = hyd_atoms[lab == c]
            patches.append({"area": float(b_sasa_c[members].sum()),
                            "res": sorted({b_of[int(bnd.resnum[m])] for m in members})})
        patches.sort(key=lambda p: -p["area"])
    a.patches = patches[:MAX_PATCHES]
    a.hphob_patch = float(patches[0]["area"]) if patches else 0.0

    # ---- interface composition
    a.comp = {"binder": composition([seq[i - 1] for i in a.bc]),
              "target": composition([tletters[tr.index_of(t)] for t in a.tc])}
    bcl = [seq[i - 1] for i in a.bindcraft_bc]
    a.interface_aas = composition(bcl)["aa"]

    # ---- disulfides: the file's record when it has one, else the SG-SG geometry (see disulfides.py for why)
    sg_idx = np.nonzero((bnd.name == "SG") & (bnd.element == "S"))[0]
    sg_xyz = {b_of[int(bnd.resnum[k])]: bnd.xyz[k] for k in sg_idx}
    recorded = [(b_of[x], b_of[y]) for x, y in cx.recorded_disulfides() if x in b_of and y in b_of]
    ds = D.assess(sg_xyz, recorded)
    a.ssbond_basis, a.disulfides = ds["basis"], ds["pairs"]
    a.ssbonds = [[p["a"], p["b"]] for p in ds["pairs"]]
    for p in ds["pairs"]:
        if p["band"] == "recorded" and p["d"] is not None and p["d"] > D.LOOSE_MAX:
            notes.append(f"the file records a disulfide between binder residues {p['a']} and {p['b']} but their sulfurs are "
                         f"{p['d']:.1f} A apart; the record is kept")

    # ---- hotspots of the design's run
    hs_all = [int(h) for h in hotspots_construct]
    present = [h for h in hs_all if tr.has(h)]
    missing = [h for h in hs_all if not tr.has(h)]
    if missing:
        notes.append(f"hotspots absent from this design's target: {missing}")
    a.hotspots = present
    btree = cKDTree(bnd.xyz)
    for h in present:
        ri = tr.index_of(h)
        pts = tgt.xyz[tr.start[ri]:tr.stop[ri]]
        a.hotspot_dist[h] = float(btree.query(pts)[0].min())
    tc_set = set(a.tc)
    a.hot = [h for h in present if h in tc_set]      # tc intersected with the hotspots, as the exporter asserts
    by_dist = [h for h in present if a.hotspot_dist[h] <= HOTSPOT_TOUCH_CUTOFF]
    if by_dist != a.hot:
        notes.append(f"hotspots by distance {by_dist} differ from tc intersection {a.hot} (boundary rounding)")

    # ---- sequence properties and liabilities
    a.seq_props = L.sequence_properties(seq)
    a.liab = L.find_liabilities(seq, a.bc, [x if x is not None else float("nan") for x in a.rsa_res],
                                patches=a.patches, disulfides=ds)
    a.global_flags = L.global_flags(seq)

    # ---- per-residue tracks
    pl = br.plddt
    a.plddt_res = [float(x) if x == x else None for x in pl]
    a.plddt_t = {int(num): (float(v) if v == v else None) for num, v in zip(tr.number, tr.plddt)}
    pairs_per_b = np.zeros(n)
    tgt_res_per_b = np.zeros(n)
    for r in a.contacts:
        pairs_per_b[r["b"] - 1] += r["n"]
        tgt_res_per_b[r["b"] - 1] += 1
    kd = [S.KYTE_DOOLITTLE.get(c, 0.0) for c in seq]
    sev_rank = {"bad": 3, "warn": 2, "info": 1}
    worst: dict[int, str] = {}
    for e in a.liab:
        if sev_rank[e["sev"]] > sev_rank.get(worst.get(e["pos"], ""), 0):
            worst[e["pos"]] = e["sev"]
    a.tracks = {
        "plddt": a.plddt_res, "rsa": a.rsa_res, "dsasa": a.dsasa_res,
        "iface": [1 if i + 1 in bc_set else 0 for i in range(n)],
        "n_target_res": [int(x) for x in tgt_res_per_b], "n_atom_pairs": [int(x) for x in pairs_per_b],
        "dist": [float(x) for x in res_min_d], "kd": kd, "ss": a.ss,
        "liab": [{"pos": p, "sev": s_} for p, s_ in sorted(worst.items())],
    }

    # ---- comparison with BindCraft's columns
    a.meta = {"version": ANALYSIS_VERSION, "hydrogens": has_h, "structure_model": design.structure_model if design else None,
              "buns": {"min_burial": buns_min_burial, "max_sasa": buns_max_sasa}, "n_points": n_points,
              "notes": notes, "hotspots_missing": missing}
    if design is not None:
        a.vs_bindcraft = compare_with_bindcraft(a, design, cx, bnd)
    return a


# --------------------------------------------------------------------------------------------- BindCraft comparison
#: stem -> (how the recomputed number is defined, same definition as BindCraft's?)
VS_NOTES = {
    "n_InterfaceResidues": ("binder residues within 4.0 A of any target atom, hydrogens counted", True),
    "Interface_Hydrophobicity": ("percent of those residues in A C F I L M P V W Y", True),
    "n_InterfaceHbonds": ("target-binder hydrogen bonds from the file's hydrogens (geometry, no energy term)", False),
    "n_InterfaceUnsatHbonds": ("buried polar atoms without an H-bond in the complex, both sides (not Rosetta's delta filter)", False),
    "dSASA": ("buried area, Shrake-Rupley with Bondi radii, both sides", False),
    "Interface_SASA_%": ("buried area of both sides over the binder's SASA in the complex (reading identified by comparison with the CSV)", False),
    "Interface_Helix%": ("share of the 4.0 A heavy-atom interface residues in DSSP H, G, I", False),
    "Interface_BetaSheet%": ("share of the 4.0 A heavy-atom interface residues in DSSP E", False),
    "Interface_Loop%": ("share of the 4.0 A heavy-atom interface residues in other DSSP states", False),
    "Binder_Helix%": ("share of binder residues in DSSP H, G, I (DSSP of the relaxed complex)", False),
    "Binder_BetaSheet%": ("share of binder residues in DSSP E", False),
    "Binder_Loop%": ("share of binder residues in other DSSP states", False),
    "Surface_Hydrophobicity": ("approximation of Rosetta's surface-layer fraction, binder alone", False),
    "InterfaceHbondsPercentage": ("n_InterfaceHbonds per 4.0 A interface residue, percent", False),
    "InterfaceUnsatHbondsPercentage": ("n_InterfaceUnsatHbonds per 4.0 A interface residue, percent", False),
}
#: stems BindCraft reports that cannot be recomputed from the relaxed complex alone
NOT_RECOMPUTABLE = {
    "Hotspot_RMSD": "needs BindCraft's trajectory target; the file holds only the predicted complex",
    "dG": "Rosetta interface energy",
    "ShapeComplementarity": "Rosetta shape complementarity",
    "PackStat": "Rosetta packing statistic",
    "i_pAE": "needs the predicted aligned error matrix, which BindCraft 1 does not write",
}


def _percent_ss(ss: Sequence[str], idx: Optional[Sequence[int]] = None) -> tuple[float, float, float]:
    """(helix, sheet, loop) percent over the residues at 1-based ``idx`` (all when None)."""
    sel = [ss[i - 1] for i in idx] if idx is not None else list(ss)
    if not sel:
        return (float("nan"),) * 3
    return tuple(100.0 * sel.count(c) / len(sel) for c in "HEL")  # type: ignore[return-value]


def compare_with_bindcraft(a: DesignAnalysis, design: DesignRecord, cx: S.Complex, bnd: S.Atoms) -> dict:
    """Recomputed values next to the copied BindCraft values, in the CSV scale.

    The BindCraft value is the one of the AlphaFold model whose pose is in the file (``_model{k}``;
    ``design.pose_metric``). ``InterfaceResidues`` in the CSV is always model 2's, so its comparison is exact only
    when the file is model 2 (``same_model`` False otherwise).
    """
    k = design.structure_model
    bcset = a.bindcraft_bc
    ss_b = _percent_ss(a.ss)
    ss_i = _percent_ss(a.ss, a.bindcraft_bc_heavy)
    seq = bnd.residues.sequence
    hyd = 100.0 * sum(1 for i in bcset if seq[i - 1] in HYDROPHOBIC_AA) / len(bcset) if bcset else 0.0
    in_complex = a.binder_sasa_alone - a.dsasa_binder     # the binder's own SASA in the complex
    mine = {
        "n_InterfaceResidues": float(len(bcset)),
        "Interface_Hydrophobicity": hyd,
        "n_InterfaceHbonds": float(a.n_hbonds),
        "n_InterfaceUnsatHbonds": float(a.n_unsat),
        "dSASA": a.dsasa,
        "Interface_SASA_%": 100.0 * a.dsasa / in_complex if in_complex else float("nan"),
        "Interface_Helix%": ss_i[0], "Interface_BetaSheet%": ss_i[1], "Interface_Loop%": ss_i[2],
        "Binder_Helix%": ss_b[0], "Binder_BetaSheet%": ss_b[1], "Binder_Loop%": ss_b[2],
        "Surface_Hydrophobicity": S.surface_hydrophobicity(bnd),
        "InterfaceHbondsPercentage": 100.0 * a.n_hbonds / len(bcset) if bcset else float("nan"),
        "InterfaceUnsatHbondsPercentage": 100.0 * a.n_unsat / len(bcset) if bcset else float("nan"),
    }
    out: dict = {}
    for stem, val in mine.items():
        ref = design.pose_metric(stem)
        note, same = VS_NOTES[stem]
        out[stem] = {"recomputed": val, "bindcraft": ref, "average": design.metrics.get(stem), "model": k,
                     "diff": None if ref is None or val != val else float(val - ref),
                     "same_definition": same, "note": note}
    csv_set = sorted(int(x[1:]) for x in design.interface_residues)
    if csv_set:
        inter = len(set(csv_set) & set(bcset))
        union = len(set(csv_set) | set(bcset))
        out["InterfaceResidues"] = {"recomputed": bcset, "bindcraft": csv_set, "model": 2, "same_model": k == 2,
                                    "jaccard": inter / union if union else 1.0, "exact": csv_set == bcset,
                                    "same_definition": True,
                                    "note": "CSV column is model 2's interface; exact only for _model2 files"}
    for stem, why in NOT_RECOMPUTABLE.items():
        out[stem] = {"recomputed": None, "bindcraft": design.pose_metric(stem), "model": k, "diff": None,
                     "same_definition": False, "note": "not recomputable: " + why}
    return out


# --------------------------------------------------------------------------------------------- run level
def file_sha256(path: str | os.PathLike) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


_CODE_FINGERPRINT: Optional[str] = None


def code_fingerprint() -> str:
    """Hash of the source files the numbers depend on, so a code change can never serve a stale cache."""
    global _CODE_FINGERPRINT
    if _CODE_FINGERPRINT is None:
        here = Path(__file__).resolve().parent
        h = hashlib.sha256()
        for name in ("interface.py", "liabilities.py", "disulfides.py", "structure.py", "model.py"):
            try:
                h.update((here / name).read_bytes())
            except OSError:
                h.update(name.encode())
        _CODE_FINGERPRINT = h.hexdigest()[:16]
    return _CODE_FINGERPRINT


def _cache_key(path: str, hotspots: Sequence[int], design: DesignRecord) -> str:
    """Structure hash plus everything else the analysis depends on (hotspots, the BindCraft columns it is compared with)."""
    # "nb": the same file analysed with one binder chain and with two gives
    # different numbers, and the file's hash cannot tell them apart.
    meta = json.dumps({"v": ANALYSIS_VERSION, "code": code_fingerprint(), "hot": sorted(int(h) for h in hotspots), "model": design.structure_model,
                       "nb": binder_chain_count(design),
                       "buns": [BUNS_MIN_BURIAL, BUNS_MAX_SASA], "points": DEFAULT_N_POINTS, "metrics": design.metrics,
                       "models": {str(k): v for k, v in design.models.items()},
                       "iface": design.interface_residues}, sort_keys=True, default=str)
    return hashlib.sha256((file_sha256(path) + meta).encode()).hexdigest()[:24]


def _work(args: tuple) -> tuple[str, dict]:
    key, design, hotspots, loader_name = args
    cx = load_for_analysis(design.structure_path, n_binder_chains=binder_chain_count(design))
    a = analyze_design(cx, design, hotspots)
    return key, a.to_dict()


def _hotspots_for(design: DesignRecord, hotspots: Any) -> list[int]:
    if hotspots is None:
        return []
    if callable(hotspots):
        return list(hotspots(design))
    if isinstance(hotspots, Mapping):
        return list(hotspots.get(design.run_id, []))
    return list(hotspots)


def analyze_run(designs: Iterable[DesignRecord], loader: Optional[Callable[[str], S.Complex]] = None,
                n_jobs: Optional[int] = None, cache_dir: Optional[str | os.PathLike] = None,
                hotspots: Any = None, progress: Optional[Callable[[int, int], None]] = None) -> dict[str, DesignAnalysis]:
    """Analyse many designs. Returns ``{design id: DesignAnalysis}``.

    designs    designs with a structure (``DesignRecord.structure_path``); others are skipped.
    loader     ``path -> Complex``. Default :func:`load_for_analysis`. A custom loader forces ``n_jobs=1``
               (a function cannot always be sent to a worker).
    n_jobs     worker processes (spawn); default CPU count minus one; 1 runs in this process.
    cache_dir  directory for the disk cache (``<cache_dir>/<run>/<id>.json``), keyed by the structure file's
               SHA-256 plus the hotspots and the analysis version. A hit is read back, a miss is computed and
               written. None disables the cache.
    hotspots   the hotspots (construct numbers) of each design's own run: a mapping ``{run_id: [numbers]}``,
               a single list used for all designs, or a callable ``design -> [numbers]``.
    progress   called as ``progress(done, total)``.

    Design ids are the keys; if two designs of different runs share an id every key becomes ``"<run>/<id>"``.
    """
    todo = [d for d in designs if d.structure_path]
    ids = [d.id for d in todo]
    clash = len(set(ids)) != len(ids)

    def key_of(d: DesignRecord) -> str:
        return f"{d.run_id}/{d.id}" if clash else d.id

    results: dict[str, DesignAnalysis] = {}
    cache = Path(cache_dir) if cache_dir else None
    pending: list[tuple[str, DesignRecord, list[int], str]] = []
    paths: dict[str, Path] = {}
    keys: dict[str, str] = {}
    for d in todo:
        hs = _hotspots_for(d, hotspots)
        k = key_of(d)
        if cache is not None:
            ck = _cache_key(d.structure_path, hs, d)
            p = cache / d.run_id / f"{d.id}.json"
            paths[k], keys[k] = p, ck
            if p.is_file():
                try:
                    blob = json.loads(p.read_text())
                    if blob.get("key") == ck:
                        results[k] = DesignAnalysis.from_dict(blob["analysis"])
                        continue
                except (OSError, ValueError):
                    pass
        pending.append((k, d, hs, ""))
    done = len(results)
    total = len(todo)
    if progress and done:
        progress(done, total)

    def store(k: str, ad: dict) -> None:
        results[k] = DesignAnalysis.from_dict(ad)
        if cache is not None:
            p = paths[k]
            p.parent.mkdir(parents=True, exist_ok=True)
            tmp = p.with_suffix(".tmp")
            tmp.write_text(json.dumps({"key": keys[k], "analysis": ad}))
            tmp.replace(p)

    workers = n_jobs if n_jobs is not None else max(1, (os.cpu_count() or 2) - 1)
    if loader is not None or workers <= 1 or len(pending) <= 2:
        for k, d, hs, _ in pending:
            cx = (loader or load_for_analysis)(d.structure_path)
            store(k, analyze_design(cx, d, hs).to_dict())
            done += 1
            if progress:
                progress(done, total)
    else:
        ctx = mp.get_context("spawn")
        with cf.ProcessPoolExecutor(max_workers=min(workers, len(pending)), mp_context=ctx) as ex:
            futs = {ex.submit(_work, (k, d, hs, "")): k for k, d, hs, _ in pending}
            for fut in cf.as_completed(futs):
                k = futs[fut]
                _, ad = fut.result()
                store(k, ad)
                done += 1
                if progress:
                    progress(done, total)
    return results


# --------------------------------------------------------------------------------------------- validation table
#: tolerance per stem: (value, "abs" in the CSV unit or "rel" as a share of BindCraft's value)
VALIDATION_TOLERANCES: dict[str, tuple[float, str]] = {
    "n_InterfaceResidues": (0.0, "abs"),            # a count with the same rule: must be exact
    "Interface_Hydrophobicity": (0.05, "abs"),      # percent; the CSV rounds to 2 decimals
    "n_InterfaceHbonds": (1.0, "abs"),              # one bond
    "n_InterfaceUnsatHbonds": (1.0, "abs"),         # one atom
    "dSASA": (0.05, "rel"),                         # different radii: about 2.5 percent is expected
    "Interface_SASA_%": (5.0, "abs"),               # percentage points
    "Interface_Helix%": (5.0, "abs"), "Interface_BetaSheet%": (5.0, "abs"), "Interface_Loop%": (5.0, "abs"),
    "Binder_Helix%": (5.0, "abs"), "Binder_BetaSheet%": (5.0, "abs"), "Binder_Loop%": (5.0, "abs"),
    "Surface_Hydrophobicity": (0.02, "abs"),        # fraction; the CSV rounds to 2 decimals
    "InterfaceHbondsPercentage": (5.0, "abs"), "InterfaceUnsatHbondsPercentage": (5.0, "abs"),
}


def validation_table(analyses: Iterable[DesignAnalysis], stems: Optional[Sequence[str]] = None,
                     reference: str = "bindcraft") -> list[dict]:
    """Agreement between recomputed values and BindCraft's columns over many designs.

    ``reference`` is "bindcraft" (the column of the model whose pose is in the file, the default) or "average"
    (the ``Average_`` column over the validated models).

    One row per stem: ``n``, Pearson ``r``, Spearman ``rho`` (None when a side is constant), ``mean_diff``
    (recomputed minus BindCraft), ``mad`` (mean absolute difference), ``mean_bindcraft``, ``tol`` with
    ``tol_kind`` ("abs" or "rel") and ``within`` (share of designs inside the tolerance).
    """
    from scipy.stats import spearmanr

    analyses = list(analyses)
    out = []
    for stem in (stems or list(VALIDATION_TOLERANCES)):
        x, y = [], []
        for a in analyses:
            v = a.vs_bindcraft.get(stem)
            ref = v.get(reference) if v else None
            if ref is not None and v.get("recomputed") is not None and v["recomputed"] == v["recomputed"]:
                x.append(float(v["recomputed"]))
                y.append(float(ref))
        if not x:
            continue
        x_, y_ = np.array(x), np.array(y)
        const = x_.std() == 0 or y_.std() == 0
        tol, kind = VALIDATION_TOLERANCES.get(stem, (0.0, "abs"))
        err = np.abs(x_ - y_)
        ok = err <= (tol * np.abs(y_) if kind == "rel" else tol) + 1e-9
        out.append({"stem": stem, "n": len(x_), "r": None if const else float(np.corrcoef(x_, y_)[0, 1]),
                    "rho": None if const else float(spearmanr(x_, y_)[0]), "mean_diff": float((x_ - y_).mean()),
                    "mad": float(err.mean()), "mean_bindcraft": float(y_.mean()), "tol": tol, "tol_kind": kind,
                    "within": float(ok.mean())})
    return out
