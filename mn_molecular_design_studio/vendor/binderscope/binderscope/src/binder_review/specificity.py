"""Specificity panel: orthologs and paralogs of the target aligned to the target's reference sequence.

Everything here works from public sequences in the offline reference cache (refdata). No binder
sequence or design id is involved, and nothing leaves the machine. The output is the HOM blob of
docs/data_contract.md section 5 (docs/schema/homolog.schema.json).

The panel is data: the target's UniProt accession (the reference sequence), the species of the orthologs and
the paralogs come from the profile (``specificity.species`` as ``acc:<accession>`` matches, ``specificity.paralogs.members``
as accessions) or from a panel file in the reference cache (``panel.json``, written by scripts/fetch_reference.py), or
are passed to ``build_panel``. Without an accession, or with neither species nor paralogs, there is no panel:
``empty_homolog_json`` writes the blob with no orthologs and no paralogs, and the pages hide the specificity features.

Residue numbers are UniProt numbers throughout (the reference structure's numbering must equal them). Construct
numbers never appear here: the caller hands in `resis` and per-run hotspots already in UniProt numbers.

Method in short (details and caveats in docs/specificity.md)
  * every member is aligned globally to the human reference with Bio.Align.PairwiseAligner (BLOSUM62,
    gap open -10, extend -0.5, as the mockup did), progressive pairwise to the reference, not a
    multiple alignment: no MSA tool is installed and the panel is a star around the reference;
  * paralogs can be too distant (below about 30 percent identity) for a sequence alignment to place the loops, so
    when AlphaFold DB models of the reference and the paralog are in the cache the per-residue map comes
    from a TM-align structural alignment (a CA pair further apart than STRUCT_MAX_DIST after
    superposition is low confidence); the whole-protein identity percentages stay on the sequence
    alignment. Without models the paralog falls back to the sequence alignment;
  * in a sequence alignment a human residue is "low confidence" when the identity in the 15-residue
    window around it is below LOWCONF_MIN_IDENT or it lies within that window of an indel of
    INDEL_MIN residues or more (a screen for loops that aligned by chance, not a verdict);
  * grades are rank based (1 variable to 9 conserved) and computed from the orthologs only; the
    paralog track is a separate set of arrays and means "conserved across the family".
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import re
import sys
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from . import __version__
from . import refdata as rd

SCHEMA_ID = "binderscope.homolog/1"

#: Bio.Align parameters of every alignment (mockup parity: homolog_conservation.py).
ALIGN_OPEN = -10.0
ALIGN_EXTEND = -0.5
#: A reference position is low confidence in a member when the share of identical positions in the window
#: of 2 * LOWCONF_HALF_WINDOW + 1 residues around it is below this. A close paralog has few flagged
#: positions, distant ones many, and a shuffled sequence aligned instead 43 to 61 percent (docs/specificity.md):
#: the flag marks windows that are not much better than chance. It is a screen, not a verdict.
LOWCONF_HALF_WINDOW = 7
LOWCONF_MIN_IDENT = 0.25
#: A sequence alignment is also unreliable next to a long gap (an insertion in either sequence).
INDEL_MIN = 8
#: Structural alignment of paralogs: an aligned CA pair further apart than this (after the TM-align
#: superposition of the whole chain) does not superpose, so the pair is low confidence.
STRUCT_MAX_DIST = 4.0
#: Ortholog identity to the reference below this suggests a wrong accession.
MIN_ORTHOLOG_IDENTITY = rd.MIN_ORTHOLOG_IDENTITY
#: Paralogs whose whole-protein identity to the reference is below this are called distant in the method note.
DISTANT_IDENTITY = 0.40
#: Tie handling and cover floor of the grade formula (data_contract.md 5).
COVER_FLOOR = 0.3
N_GRADES = 9

AA = "ACDEFGHIKLMNPQRSTVWY"
#: The viewer's notion of a similar (not identical) residue (common.js S.similar).
SIMILAR_GROUPS = ("ILVM", "FYW", "KRH", "DE", "ST", "NQ", "AG")


def similar(a: str, b: str) -> bool:
    """True for two different residues of the same group (the viewer's definition)."""
    return a != b and any(a in g and b in g for g in SIMILAR_GROUPS)


# ---------------------------------------------------------------------------------------------
# alignment
# ---------------------------------------------------------------------------------------------

def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices

    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = ALIGN_OPEN
    al.extend_gap_score = ALIGN_EXTEND
    al.mode = "global"
    return al


def align_map(reference: str, other: str, aligner=None, *, with_indels: bool = False):
    """The residue of `other` aligned to each residue of `reference` ('-' for a gap), a string of
    len(reference). Insertions in `other` are dropped; the first optimal alignment is used.

    With with_indels the result is (map, near) where near is a boolean array: the position lies within
    LOWCONF_HALF_WINDOW of a gap in either sequence that is INDEL_MIN residues or longer (a long gap
    means the neighbouring register is a guess)."""
    al = aligner or _aligner()
    best = al.align(reference, other)[0]
    out = ["-"] * len(reference)
    blocks = list(zip(*best.aligned))
    for (r0, r1), (o0, o1) in blocks:
        out[r0:r1] = other[o0:o1]
    mp = "".join(out)
    if not with_indels:
        return mp
    near = np.zeros(len(reference), dtype=bool)
    for ((_, r1), (_, o1)), ((r2, _), (o2, _)) in zip(blocks, blocks[1:]):
        if max(r2 - r1, o2 - o1) >= INDEL_MIN:
            near[max(0, r1 - LOWCONF_HALF_WINDOW):min(len(reference), r2 + LOWCONF_HALF_WINDOW)] = True
    return mp, near


def window_identity(reference: str, aligned: str, half_window: int = LOWCONF_HALF_WINDOW) -> np.ndarray:
    """Share of identical aligned positions in the window around every reference position. The window
    is clipped at the ends of the sequence (the divisor is the clipped length)."""
    same = np.fromiter((a == b for a, b in zip(aligned, reference)), dtype=float, count=len(reference))
    csum = np.concatenate([[0.0], np.cumsum(same)])
    n = len(reference)
    idx = np.arange(n)
    lo = np.maximum(0, idx - half_window)
    hi = np.minimum(n, idx + half_window + 1)
    return (csum[hi] - csum[lo]) / (hi - lo)


def low_confidence(reference: str, aligned: str, min_ident: float = LOWCONF_MIN_IDENT,
                   half_window: int = LOWCONF_HALF_WINDOW) -> np.ndarray:
    """Boolean per reference position: aligned (not a gap) but inside a poorly aligned window."""
    aligned_arr = np.fromiter((c != "-" for c in aligned), dtype=bool, count=len(aligned))
    return aligned_arr & (window_identity(reference, aligned, half_window) < min_ident)


def _model_ca(path: Path):
    """(residue numbers, one-letter sequence, CA coordinates) of the first chain of a model file."""
    import gemmi

    ch = gemmi.read_structure(str(path))[0][0]
    nums, seq, xyz = [], [], []
    for r in ch:
        ca = r.find_atom("CA", "*")
        if ca is None:
            continue
        nums.append(r.seqid.num)
        info = gemmi.find_tabulated_residue(r.name)
        seq.append(info.one_letter_code.upper() if info else "X")
        xyz.append([ca.pos.x, ca.pos.y, ca.pos.z])
    return nums, "".join(seq), np.array(xyz)


@lru_cache(maxsize=32)
def structure_map(reference_model: str, other_model: str, max_dist: float = STRUCT_MAX_DIST) -> tuple[str, np.ndarray, dict[str, float]]:
    """Structure-guided map of a paralog onto the reference from two model files (AlphaFold DB, UniProt
    numbering).

    TM-align (tmtools) superposes the chains and gives a sequential residue pairing. Returns the paralog
    residue paired with every reference position ('-' where TM-align leaves the position unpaired), the
    CA distance of every pair after the superposition (nan for a gap) and {tm, rmsd, aligned}. Pairs
    further apart than max_dist are kept in the map but are the caller's low-confidence flags."""
    from tmtools import tm_align

    rn, rs, rx = _model_ca(Path(reference_model))
    on, os_, ox = _model_ca(Path(other_model))
    res = tm_align(rx, ox, rs, os_)
    xt = rx @ np.asarray(res.u).T + np.asarray(res.t)
    out = ["-"] * len(rs)
    dist = np.full(len(rs), np.nan)
    i = j = 0
    for a, b in zip(res.seqxA, res.seqyA):
        if a != "-" and b != "-":
            out[i] = os_[j]
            dist[i] = float(np.linalg.norm(xt[i] - ox[j]))
        i += a != "-"
        j += b != "-"
    # numbering: AlphaFold DB models number residues 1..n, so index i is residue i + 1; models that
    # skip residues are shifted onto the reference numbering through rn
    full = ["-"] * (max(rn) if rn else 0)
    fdist = np.full(len(full), np.nan)
    for k, num in enumerate(rn):
        full[num - 1] = out[k]
        fdist[num - 1] = dist[k]
    info = {"tm": float(res.tm_norm_chain1), "rmsd": float(res.rmsd), "aligned": int(np.isfinite(dist).sum())}
    return "".join(full), fdist, info


# ---------------------------------------------------------------------------------------------
# per-column statistics and grades
# ---------------------------------------------------------------------------------------------

def column_stats(reference: str, resis: Sequence[int], maps: Sequence[str],
                 include_reference: bool = False, usable: Sequence[Sequence[bool]] | None = None):
    """cons, ident, cover per entry of resis from aligned members.

    maps: one full-length aligned string per member (index = reference position - 1). usable: optional
    boolean arrays (same shape) marking residues that count; others are treated as gaps. With
    include_reference the reference residue is added to every column's counts (family track).
    cons = 1 - H / log2(20) with H the Shannon entropy of the residue counts (+1e-9, as the mockup);
    ident = share of counting members that carry the reference residue; cover = share of members that
    count. A column without a counting member scores cons 0, ident 0, cover 0.
    """
    n = len(resis)
    cons = np.zeros(n)
    ident = np.zeros(n)
    cover = np.zeros(n)
    n_members = len(maps)
    for i, r in enumerate(resis):
        col = []
        for k, m in enumerate(maps):
            c = m[r - 1]
            if c == "-" or (usable is not None and not usable[k][r - 1]):
                continue
            col.append(c)
        cover[i] = len(col) / n_members if n_members else 0.0
        if not col:
            continue
        ref_res = reference[r - 1]
        ident[i] = float(np.mean([c == ref_res for c in col]))
        counted = col + [ref_res] if include_reference else col
        cnt = np.array([counted.count(a) for a in AA], float) + 1e-9
        p = cnt / cnt.sum()
        cons[i] = 1.0 - float(-(p * np.log2(p)).sum()) / np.log2(len(AA))
    return cons, ident, cover


def grades_from(cons: np.ndarray, cover: np.ndarray) -> np.ndarray:
    """ConSurf-style grade 1 (variable) to 9 (conserved): rank of cons * max(cover, COVER_FLOOR) over
    the positions in nine bins.

    Tied scores keep one grade (a tie is never split by chance). A tie at the top that is larger than
    its share of a bin (a small panel has many columns that are invariant) takes the highest grade still
    free, and the remaining positions are ranked again into the remaining grades; otherwise the bins
    would stay empty or, with average ranks as in the mockup, put the invariant majority into the middle
    grades. Monotone: a higher score never gets a lower grade, and a column that is invariant in every
    member while all members cover it grades 9 whenever invariant columns outnumber one bin."""
    from scipy.stats import rankdata

    score = np.asarray(cons, float) * np.maximum(np.asarray(cover, float), COVER_FLOOR)
    n = score.size
    grade = np.ones(n, dtype=int)
    rem = np.arange(n)
    top = N_GRADES
    while rem.size:
        sc = score[rem]
        if top > 1:
            block = sc == sc.max()
            if block.sum() > rem.size / top:
                grade[rem[block]] = top
                rem = rem[~block]
                top -= 1
                continue
        rk = rankdata(sc, method="max") / rem.size
        grade[rem] = np.clip(np.ceil(rk * top), 1, top).astype(int)
        break
    return grade


# ---------------------------------------------------------------------------------------------
# data classes
# ---------------------------------------------------------------------------------------------

@dataclass
class Member:
    """One aligned ortholog or paralog."""

    kind: str                      # 'ortholog' | 'paralog'
    key: str                       # profile species key (mouse) or gene symbol of a paralog
    acc: str
    label: str
    name: str                      # display name, same as label for paralogs
    gene: str
    organism: str
    taxon: int | None
    reviewed: bool
    length: int
    role: str                      # 'default' | 'optional' (refdata panel)
    full_map: str                  # aligned residue per human position 1..len(reference), '-' for a gap
    low: np.ndarray                # low-confidence flag per human position (0-based index)
    extra: bool = False
    match: str = ""
    seq_map: str | None = None     # sequence alignment, kept when full_map comes from a structural alignment
    aligner: str = "sequence"      # 'sequence' | 'structure': what full_map was derived from
    struct: dict[str, float] | None = None   # TM-align tm, rmsd, aligned (structure only)

    @property
    def optional(self) -> bool:
        return self.role == "optional"

    def identical_count(self, reference: str, residues: Iterable[int] | None = None) -> int:
        """Identical residues over the given residues (the per-residue map), or over the whole protein
        from the sequence alignment (the figure UniProt and refdata quote)."""
        mp = self.full_map if residues is not None else (self.seq_map or self.full_map)
        rs = range(1, len(reference) + 1) if residues is None else residues
        return sum(1 for r in rs if mp[r - 1] == reference[r - 1])


@dataclass
class SpecificityPanel:
    reference: dict[str, Any]
    ref_seq: str
    resis: list[int]
    members: list[Member]
    tracks: dict[str, Any]
    spec: dict[str, Any]           # mirrors of profile.specificity: species, show, subject, kind
    consurf: dict[str, Any] | None = None
    warnings: list[str] = field(default_factory=list)
    params: dict[str, Any] = field(default_factory=dict)
    observed: tuple[int, int] | None = None   # residues of the reference chain that the reference structure observes (first, last)

    @property
    def orthologs(self) -> list[Member]:
        return [m for m in self.members if m.kind == "ortholog"]

    @property
    def paralogs(self) -> list[Member]:
        return [m for m in self.members if m.kind == "paralog"]

    @property
    def seq(self) -> str:
        return "".join(self.ref_seq[r - 1] for r in self.resis)

    def member(self, key: str) -> Member:
        for m in self.members:
            if m.key == key or m.acc == key:
                return m
        raise KeyError(f"no member {key!r}; have {[m.key for m in self.members]}")

    def label_of(self, residue: int) -> str:
        return f"{self.ref_seq[residue - 1]}{residue}"


# ---------------------------------------------------------------------------------------------
# identity at a residue set
# ---------------------------------------------------------------------------------------------

def _overlap(reference: str, full_map: str, residues: Sequence[int]) -> dict[str, Any]:
    ident = sim = gap = 0
    differing: list[str] = []
    gaps: list[str] = []
    for r in residues:
        if not 1 <= r <= len(reference):
            raise ValueError(f"residue {r} is outside the reference (1 to {len(reference)})")
        a, b = full_map[r - 1], reference[r - 1]
        if a == "-":
            gap += 1
            gaps.append(f"{b}{r}")
        elif a == b:
            ident += 1
        else:
            differing.append(f"{b}{r}{a}")
            if similar(a, b):
                sim += 1
    n = len(residues)
    return {
        "n": n, "identical": ident, "similar": sim, "gap": gap,
        "identity": ident / n if n else 0.0,
        "similarity": (ident + sim) / n if n else 0.0,
        "differing": differing, "gaps": gaps,
    }


def identity_at(panel: SpecificityPanel, key: str, residues_uniprot: Sequence[int]) -> dict[str, Any]:
    """Identity and similarity of one member at a set of human residues (UniProt numbers).

    Returns n, identical, similar, gap (counts), identity = identical / n, similarity = (identical +
    similar) / n with the viewer's residue groups, differing (['A57K', ...] reference residue, number,
    member residue; positions without an aligned residue are in gaps, not in differing) and lowconf
    (labels of residues in poorly aligned windows). The pages compute the same numbers at run time
    from the aligned maps."""
    m = panel.member(key)
    out = _overlap(panel.ref_seq, m.full_map, list(residues_uniprot))
    out["lowconf"] = [panel.label_of(r) for r in residues_uniprot if m.low[r - 1]]
    return out


def region_identity(panel: SpecificityPanel, key: str, first: int, last: int) -> dict[str, Any]:
    """identity_at over every residue from first to last (inclusive)."""
    return identity_at(panel, key, list(range(first, last + 1)))


def hotspot_identity_table(panel: SpecificityPanel, hotspots_by_run: Mapping[str, Sequence[int]],
                           keys: Sequence[str] | None = None) -> list[dict[str, Any]]:
    """One row per (run, member): identity at the run's hotspots (UniProt numbers). Rows hold run,
    member key, kind, label, hotspots (labels), aligned (the member's residue at each hotspot, '-'
    for a gap) and the identity_at fields."""
    rows = []
    members = [panel.member(k) for k in keys] if keys else panel.members
    for run, hs in hotspots_by_run.items():
        hs = list(hs)
        for m in members:
            r = identity_at(panel, m.key, hs)
            rows.append({
                "run": run, "key": m.key, "kind": m.kind, "label": m.label,
                "hotspots": [panel.label_of(h) for h in hs],
                "aligned": [m.full_map[h - 1] for h in hs],
                **r,
            })
    return rows


# ---------------------------------------------------------------------------------------------
# building
# ---------------------------------------------------------------------------------------------

def _species_key(common_name: str) -> str:
    return re.sub(r"[^a-z0-9_]+", "_", (common_name or "").lower()).strip("_")


def default_species(panel: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Species entries of a panel file in profile form ({match, label, key, optional}): the key is the lower-case common name."""
    out = []
    for e in panel.get("orthologs") or []:
        key = _species_key(e["common_name"])
        ent = {"match": f"acc:{e['accession']}", "label": e["common_name"].capitalize(), "key": key}
        if e.get("role") == "optional":
            ent["optional"] = True
        out.append(ent)
    return out


def _load_profile_spec(profile: Mapping[str, Any] | None, rpanel: Mapping[str, Any], target: Mapping[str, Any]) -> dict[str, Any]:
    spec = dict((profile or {}).get("specificity") or {})
    spec.setdefault("subject", target.get("gene") or "target")
    spec.setdefault("kind", "host")
    if not spec.get("species"):
        spec["species"] = default_species(rpanel)
    spec.setdefault("show", [x["key"] for x in spec["species"][:2]])
    return spec


def _sequence_of(acc: str) -> str:
    return rd.fasta_sequence(rd.get_uniprot_fasta(acc))


def _consurf_overlay(ref_seq: str, resis: Sequence[int], warnings: list[str], pdb: str | None, chain: str) -> dict[str, Any] | None:
    """The ConSurf-DB grades of a PDB chain (optional overlay); None when the reference is not a PDB entry or has no record."""
    try:
        rows = rd.get_consurfdb(pdb, chain) if pdb else None
    except (rd.RefDataError, ValueError):
        rows = None
    if not rows:
        return None
    by_num = {}
    for r in rows:
        try:
            by_num[int(r["pdb_resnum"])] = r
        except (TypeError, ValueError):
            continue
    grades: list[int | None] = []
    insufficient: list[int] = []
    bad = 0
    for r in resis:
        row = by_num.get(r)
        if row is None:
            grades.append(None)
            continue
        if row["aa"] != ref_seq[r - 1]:
            bad += 1
            grades.append(None)
            continue
        grades.append(int(row["grade"]))
        if row["insufficient"]:
            insufficient.append(r)
    if bad:
        warnings.append(f"ConSurf-DB residue letters disagree with the reference at {bad} positions; those grades are dropped")
    entry = rd.manifest().get("entries", {}).get(f"consurfdb/{pdb.upper()}_{chain}_consurf_grades.txt", {})
    return {
        "source": "ConSurf-DB (consurfdb.tau.ac.il)",
        "pdb": pdb.upper(),
        "chain": chain,
        "retrieved": entry.get("retrieved") or _dt.date.today().isoformat(),
        "grade": grades,
        "insufficient": insufficient,
        "licence": "ConSurf-DB data is for academic, non-commercial use. Off by default; its alignment mixes "
                   "orthologs and paralogs, unlike the ortholog-only grades of this panel.",
    }


def _reference_model(reference: str, path: str | Path | None, warnings: list[str]) -> str | None:
    """The reference model file when it can serve as the structural reference (UniProt numbering and the
    reference sequence), else None with a warning."""
    if path is None:
        return None
    try:
        nums, seq, _ = _model_ca(Path(path))
    except Exception as exc:  # unreadable model: fall back to sequence alignment
        warnings.append(f"reference model {path} unreadable ({exc}); paralogs use the sequence alignment")
        return None
    if any(n < 1 or n > len(reference) or reference[n - 1] != c for n, c in zip(nums, seq)):
        warnings.append("reference model does not match the reference sequence; paralogs use the sequence alignment")
        return None
    return str(path)


def _apply_structure(m: Member, ref_model: str, model: str | Path, sequence: str, reference: str,
                     warnings: list[str]) -> None:
    """Replace the per-residue map of a paralog by its structural alignment to the reference. The sequence
    map stays in seq_map (whole-protein identity percentages). A pair that does not superpose within
    STRUCT_MAX_DIST is low confidence; an unpaired position is a gap."""
    try:
        sm, dist, info = structure_map(ref_model, str(model))
    except Exception as exc:
        warnings.append(f"paralog {m.key} ({m.acc}): structural alignment failed ({exc}); using the sequence alignment")
        return
    n = len(reference)
    fm = (sm + "-" * n)[:n]
    d = np.full(n, np.nan)
    d[:min(n, dist.size)] = dist[:n]
    m.seq_map, m.full_map, m.aligner, m.struct = m.full_map, fm, "structure", info
    m.low = np.isfinite(d) & (d > STRUCT_MAX_DIST)


def build_panel_from_sequences(reference: str, orthologs: Sequence[Mapping[str, Any]],
                               paralogs: Sequence[Mapping[str, Any]], resis: Sequence[int], *,
                               spec: Mapping[str, Any] | None = None,
                               reference_info: Mapping[str, Any] | None = None,
                               track_paralogs: Sequence[str] | None = None,
                               reference_model: str | Path | None = None,
                               observed: Sequence[int] | None = None) -> SpecificityPanel:
    """Build a panel from sequences without touching refdata (tests, toy alignments).

    orthologs and paralogs are dicts with `sequence`, `key`, `acc` and optional label, name, gene,
    organism, taxon, reviewed, role, extra and, for paralogs, `model` (a structure file). A paralog with
    a model is compared with the reference by structural alignment when reference_model is given.
    track_paralogs: gene keys that enter the paralog track (default: those whose role is not 'optional').
    observed: [first, last] residue of the reference that the reference structure observes (a statistics scope)."""
    resis = sorted(set(int(r) for r in resis))
    if resis and (resis[0] < 1 or resis[-1] > len(reference)):
        raise ValueError(f"resis must lie within 1 to {len(reference)}")
    al = _aligner()
    members: list[Member] = []
    warnings: list[str] = []
    ref_model = _reference_model(reference, reference_model, warnings)
    for kind, items in (("ortholog", orthologs), ("paralog", paralogs)):
        for it in items:
            if kind == "paralog":
                fm, near = align_map(reference, it["sequence"], al, with_indels=True)
                low = low_confidence(reference, fm) | (near & np.array([c != "-" for c in fm]))
            else:
                fm = align_map(reference, it["sequence"], al)
                low = low_confidence(reference, fm)
            m = Member(
                kind=kind, key=it["key"], acc=it["acc"], label=it.get("label") or it["key"],
                name=it.get("name") or it.get("label") or it["key"], gene=it.get("gene") or "",
                organism=it.get("organism") or "", taxon=it.get("taxon"), reviewed=bool(it.get("reviewed", False)),
                length=len(it["sequence"]), role=it.get("role", "default"), full_map=fm,
                low=low, extra=bool(it.get("extra", False)), match=it.get("match", ""),
            )
            if kind == "paralog" and ref_model is not None and it.get("model"):
                _apply_structure(m, ref_model, it["model"], it["sequence"], reference, warnings)
            elif kind == "paralog" and ref_model is not None:
                warnings.append(f"paralog {m.key} ({m.acc}): no AlphaFold DB model in the reference cache, per-residue "
                                "comparison uses the sequence alignment and is unreliable in loops")
            members.append(m)
    panel = SpecificityPanel(
        reference=dict(reference_info or {}), ref_seq=reference, resis=resis, members=members, tracks={},
        spec=dict(spec or {}), warnings=warnings, observed=(int(observed[0]), int(observed[1])) if observed else None,
        params={"matrix": "BLOSUM62", "gap_open": ALIGN_OPEN, "gap_extend": ALIGN_EXTEND, "mode": "global",
                "lowconf_half_window": LOWCONF_HALF_WINDOW, "lowconf_min_identity": LOWCONF_MIN_IDENT,
                "indel_min": INDEL_MIN},
    )
    structural = [m for m in panel.paralogs if m.aligner == "structure"]
    if structural:
        panel.params["paralog_alignment"] = {
            "method": "TM-align of AlphaFold DB models onto the reference model", "max_pair_distance": STRUCT_MAX_DIST,
            "members": {m.key: {k: _r(v, 3) for k, v in (m.struct or {}).items()} for m in structural}}
    for m in panel.orthologs:
        idn = m.identical_count(reference) / len(reference)
        if idn < MIN_ORTHOLOG_IDENTITY and not m.extra:
            warnings.append(f"ortholog {m.key} ({m.acc}) is only {idn:.1%} identical to the reference; check the accession")
    # ortholog track: reference excluded, every non-extra ortholog counts, gaps do not
    core = [m for m in panel.orthologs if not m.extra]
    cons, ident, cover = column_stats(reference, resis, [m.full_map for m in core])
    grade = grades_from(cons, cover)
    tracks: dict[str, Any] = {"n_orth": len(core), "cons": cons, "ident": ident, "cover": cover, "grade": grade}
    # paralog track: reference included (is the residue shared across the family), low-confidence
    # residues count as gaps
    pk = set(track_paralogs) if track_paralogs is not None else {m.key for m in panel.paralogs if not m.optional}
    tp = [m for m in panel.paralogs if m.key in pk]
    pc, pi, pv = column_stats(reference, resis, [m.full_map for m in tp], include_reference=True,
                              usable=[~m.low for m in tp])
    tracks.update({"n_par": len(tp), "par_keys": [m.key for m in tp], "par_cons": pc, "par_ident": pi,
                   "par_cover": pv, "par_grade": grades_from(pc, pv) if tp else np.ones(len(resis), dtype=int)})
    panel.tracks = tracks
    return panel


def _alphafold_model(acc: str, warnings: list[str]) -> str | None:
    """Path of the cached AlphaFold DB model of an accession, or None (offline miss, no model)."""
    try:
        return str(rd.get_alphafold(acc, with_pae=False)["model_path"])
    except (rd.RefDataError, KeyError, IndexError) as exc:
        warnings.append(f"no AlphaFold DB model for {acc} in the reference cache ({exc})")
        return None


def _cache_panel() -> dict[str, Any]:
    """The panel file of the reference cache (data/reference/panel.json, written by scripts/fetch_reference.py), or an empty panel."""
    try:
        return rd.get_panel()
    except rd.RefDataError:
        return {}


def panel_requested(profile: Mapping[str, Any] | None) -> bool:
    """True when the profile asks for a panel: it names species or has a paralogs block. A profile without either has no specificity features."""
    sp = (profile or {}).get("specificity") or {}
    return bool(sp.get("species")) or bool(sp.get("paralogs"))


def build_panel(resis: Sequence[int], profile: Mapping[str, Any] | None = None, *, include_consurf: bool = True,
                accession: str | None = None, panel: Mapping[str, Any] | None = None, observed: Sequence[int] | None = None,
                pdb: str | None = None, chain: str = "A") -> SpecificityPanel:
    """Specificity panel of the target from the cached public sequences (refdata, offline-first).

    resis: target residues in UniProt numbers (campaign target.resi: the union of the runs' constructs).
    accession: UniProt accession of the target, the reference sequence (default: the profile's target.uniprot, else the target of
    ``panel``). profile: a profile dict; its specificity.species names the orthologs ({match: 'acc:<accession>' or an organism name,
    label, key, optional}) and specificity.paralogs.members the paralogs ([{acc, gene, optional}]); a paralogs block without members
    takes the paralogs of ``panel``, ordered by specificity.paralogs.order. panel: {"orthologs": [...], "paralogs": [...]} as written by
    scripts/fetch_reference.py (default: data/reference/panel.json when it exists); it supplies the details of an ortholog or paralog the
    profile names by accession and the species of a profile that names none. observed: [first, last] residue the reference structure
    observes. pdb, chain: the reference structure, for the optional ConSurf-DB overlay (a PDB id).
    """
    if resis is None:
        raise ValueError("build_panel needs the target residues (resis)")
    rpanel = dict(panel) if panel is not None else _cache_panel()
    acc = accession or ((profile or {}).get("target") or {}).get("uniprot") or (rpanel.get("target") or {}).get("accession")
    if not acc:
        raise ValueError("a specificity panel needs the UniProt accession of the target (it is the reference sequence)")
    target = rd.uniprot_summary(rd.get_uniprot(acc))
    reference = _sequence_of(acc)
    spec = _load_profile_spec(profile, rpanel, target)
    by_acc = {e["accession"]: e for e in rpanel.get("orthologs") or []}
    by_name = {}
    for e in rpanel.get("orthologs") or []:
        by_name.setdefault(e["species"], e)
        by_name.setdefault(e["organism"], e)
    warnings: list[str] = []
    orthologs = []
    for s in spec["species"]:
        match = s["match"]
        if match.startswith("acc:"):
            oacc = match[4:]
            e = by_acc.get(oacc)
            if e is None:
                try:
                    entry = rd.uniprot_summary(rd.get_uniprot(oacc))
                except (rd.RefDataError, ValueError) as exc:
                    warnings.append(f"species {s['key']} ({oacc}) skipped: not in the reference cache ({exc})")
                    continue
                e = {"accession": oacc, "gene": entry.get("gene"), "species": entry["organism"], "organism": entry["organism"],
                     "taxon_id": entry.get("taxon_id"), "reviewed": entry.get("reviewed"),
                     "role": "optional" if s.get("optional") else "default"}
        else:
            e = by_name.get(match)
            if e is None:
                warnings.append(f"species {s['key']} ({match}) skipped: no ortholog of that name in the panel")
                continue
        role = "optional" if s.get("optional") else e.get("role", "default")
        orthologs.append({
            "sequence": _sequence_of(e["accession"]), "key": s["key"], "acc": e["accession"], "label": s.get("label") or s["key"],
            "gene": e.get("gene") or "", "organism": e.get("species") or e.get("organism") or "", "taxon": e.get("taxon_id"),
            "reviewed": bool(e.get("reviewed")), "role": role, "match": match,
        })
    par_spec = spec.get("paralogs") or {}
    if par_spec.get("members"):
        par_entries = []
        for m in par_spec["members"]:
            try:
                ps = rd.uniprot_summary(rd.get_uniprot(m["acc"]))
            except (rd.RefDataError, ValueError) as exc:
                warnings.append(f"paralog {m.get('gene') or m['acc']} skipped: not in the reference cache ({exc})")
                continue
            par_entries.append({"accession": m["acc"], "gene": m.get("gene") or ps.get("gene") or m["acc"], "organism": ps.get("organism"),
                                "taxon_id": ps.get("taxon_id"), "reviewed": ps.get("reviewed"),
                                "role": "optional" if m.get("optional") else "default"})
    elif par_spec:
        par_entries = list(rpanel.get("paralogs") or [])
    else:
        par_entries = []
    order = par_spec.get("order") or []
    rank = {g: i for i, g in enumerate(order)}
    ref_model = _alphafold_model(acc, warnings)
    paralogs = []
    for e in sorted(par_entries, key=lambda e: rank.get(e["gene"], len(rank))):
        paralogs.append({
            "model": _alphafold_model(e["accession"], warnings) if ref_model else None,
            "sequence": _sequence_of(e["accession"]), "key": e["gene"], "acc": e["accession"], "label": e["gene"], "name": e["gene"],
            "gene": e["gene"], "organism": e.get("organism") or "", "taxon": e.get("taxon_id"),
            "reviewed": bool(e.get("reviewed")), "role": e.get("role", "default"),
        })
    panel_out = build_panel_from_sequences(
        reference, orthologs, paralogs, resis, spec=spec, reference_model=ref_model, observed=observed,
        reference_info={"acc": acc, "gene": target.get("gene"), "organism": target.get("organism"),
                        "length": len(reference), "reviewed": target.get("reviewed")},
    )
    panel_out.warnings[:0] = warnings
    if include_consurf:
        panel_out.consurf = _consurf_overlay(reference, panel_out.resis, panel_out.warnings, pdb, chain)
    return panel_out


def empty_homolog_json(resis: Sequence[int], letters: Mapping[int, str], *, subject: str = "", kind: str = "host",
                       note: str = "", date: str | None = None) -> dict[str, Any]:
    """The HOM blob of a project without a specificity panel: the reference residue at every entry of resis (``letters``: residue number
    to one-letter code), no orthologs, no paralogs and flat tracks. The profile's capability flags ``species_panel`` and ``paralog_panel``
    are then false and the pages hide the specificity features."""
    resis = [int(r) for r in resis]
    n = len(resis)
    return {
        "schema": SCHEMA_ID, "resis": resis, "seq": "".join(letters.get(r, "X") for r in resis),
        "cons": [0.0] * n, "ident": [0.0] * n, "grade": [1] * n, "cover": [0.0] * n, "n_orth": 0, "orth": [], "par": [],
        "par_cons": [0.0] * n, "par_ident": [0.0] * n, "par_cover": [0.0] * n, "par_grade": [1] * n, "n_par": 0, "par_keys": [],
        "species": [], "show": [], "subject": subject, "kind": kind,
        "note": note or "No specificity panel: no orthologs or paralogs were given.",
        "generated": {"tool": "binderscope specificity.py", "version": __version__, "date": date or _dt.date.today().isoformat(),
                      "params": {}, "warnings": []},
    }


# ---------------------------------------------------------------------------------------------
# HOM blob
# ---------------------------------------------------------------------------------------------

def _r(x: float, nd: int = 3) -> float:
    return round(float(x), nd)


def _seq_map(m: Member) -> str:
    """The sequence alignment of a member (equal to full_map unless a structural alignment replaced it)."""
    return m.seq_map or m.full_map


def _stat_nd(m: Member) -> int:
    """Decimals of reported identities: paralog identity moves by up to 5 points with the aligner
    settings (22 to 30 percent identical), so it is reported to whole percent."""
    return 2 if m.kind == "paralog" else 4


def member_stats(panel: SpecificityPanel, m: Member) -> dict[str, Any]:
    """Identity and similarity over the whole protein, the region the reference structure observes (when known) and the
    construct (the panel's residues), plus identity over the member's own length (the figure refdata quotes).

    full, observed and identity_of_member_length always come from the sequence alignment. construct
    comes from the map used for the per-residue comparison (the structural alignment for a paralog that
    has one); construct_sequence repeats it for the sequence alignment in that case."""
    ref = panel.ref_seq
    nd = _stat_nd(m)
    scopes = {"full": (1, len(ref))}
    if panel.observed:
        scopes["observed"] = panel.observed
    scopes["construct"] = (panel.resis[0], panel.resis[-1])
    out: dict[str, Any] = {}
    for name, (a, b) in scopes.items():
        rs = list(range(a, min(b, len(ref)) + 1)) if name != "construct" else list(panel.resis)
        o = _overlap(ref, m.full_map if name == "construct" else _seq_map(m), rs)
        out[name] = {"n": o["n"], "identity": _r(o["identity"], nd), "similarity": _r(o["similarity"], nd), "gap": o["gap"]}
    if m.aligner == "structure":
        o = _overlap(ref, _seq_map(m), list(panel.resis))
        out["construct_sequence"] = {"n": o["n"], "identity": _r(o["identity"], nd), "similarity": _r(o["similarity"], nd), "gap": o["gap"]}
    out["identity_of_member_length"] = _r(m.identical_count(ref) / m.length, nd)
    out["aligner"] = m.aligner
    return out


def _orth_json(panel: SpecificityPanel, m: Member) -> dict[str, Any]:
    ref = panel.ref_seq
    d: dict[str, Any] = {
        "acc": m.acc, "org": m.organism, "ident": _r(m.identical_count(ref) / len(ref), 4),
        "map": "".join(m.full_map[r - 1] for r in panel.resis),
        "gene": m.gene, "label": m.label, "key": m.key,
    }
    if m.taxon is not None:
        d["taxon"] = int(m.taxon)
    d["reviewed"] = m.reviewed
    d["length"] = m.length
    d["optional"] = m.optional
    if m.extra:
        d["extra"] = True
    d["lowconf"] = [r for r in panel.resis if m.low[r - 1]]
    d["stats"] = member_stats(panel, m)
    return d


def _par_json(panel: SpecificityPanel, m: Member) -> dict[str, Any]:
    ref = panel.ref_seq
    mp = "".join(m.full_map[r - 1] for r in panel.resis)
    aligned = [r for r in panel.resis if m.full_map[r - 1] != "-"]
    same = sum(1 for r in aligned if m.full_map[r - 1] == ref[r - 1])
    d: dict[str, Any] = {
        "acc": m.acc, "name": m.name, "gene": m.gene or m.key,
        "cov": _r(len(aligned) / len(panel.resis), 4) if panel.resis else 0.0,
        "ident_ref": _r(same / len(aligned), 2) if aligned else 0.0,
        "ident_full": _r(m.identical_count(ref) / len(ref), 2),
        "map": mp, "label": m.label, "key": m.key, "length": m.length,
        "optional": m.optional,
        "distant": m.optional,
        "lowconf": [r for r in panel.resis if m.low[r - 1]],
        "stats": member_stats(panel, m),
        "aligner": m.aligner,
    }
    if m.struct:
        d["struct"] = {k: _r(v, 3) for k, v in m.struct.items()}
    if m.taxon is not None:
        d["taxon"] = int(m.taxon)
    d["reviewed"] = m.reviewed
    return d


def par_method_note(panel: SpecificityPanel) -> str:
    """One honest sentence on how the paralogs were compared (shown under the paralog bars)."""
    st = [m for m in panel.paralogs if m.aligner == "structure" and not m.optional]
    sq = [m for m in panel.paralogs if m.aligner != "structure" and not m.optional]
    parts = []
    if st:
        fit = ", ".join(f"{m.key} {m.struct['rmsd']:.1f} \u00c5 over {int(m.struct['aligned'])} residues" for m in st)
        organism = panel.reference.get("organism")
        onto = ("human " if organism == "Homo sapiens" else f"{organism} " if organism else "") + str(panel.reference.get("gene", "reference"))
        far = [m.key for m in st if m.identical_count(panel.ref_seq) / len(panel.ref_seq) < DISTANT_IDENTITY]
        parts.append(
            f"Paralog residues are matched by structure (TM-align of AlphaFold DB models onto {onto}: {fit}). "
            f"A * marks a residue whose partner is more than {STRUCT_MAX_DIST:g} \u00c5 away after superposition, so the match is a guess and the residue "
            "is left out of the family conservation track."
            + (f" Loops of the distant members ({', '.join(far)}) match least well." if far else ""))
    if sq:
        parts.append(f"{', '.join(m.key for m in sq)} {'is' if len(sq) == 1 else 'are'} compared by sequence alignment only, which can misplace loops; "
                     "positions next to a long insertion or in a poorly matching window are marked *.")
    parts.append("Whole-protein identities come from the sequence alignment and move by several points with the aligner settings, so read them as approximate.")
    return " ".join(parts)


def to_homolog_json(panel: SpecificityPanel, *, date: str | None = None) -> dict[str, Any]:
    """The HOM blob as a plain dict (data_contract.md section 5), valid against homolog.schema.json.

    Additions beyond the contract, all additive: orth[].optional, .lowconf, .stats; par[].optional,
    .distant, .ident_full, .lowconf, .stats; the paralog track (par_cons, par_ident, par_cover,
    par_grade, n_par, par_keys); consurf.insufficient; generated.params and generated.warnings."""
    t = panel.tracks
    spec = panel.spec
    species = []
    for s in spec.get("species", []):
        species.append([s["match"], s.get("label") or s["key"], s["key"]])
    present = {m.key for m in panel.orthologs}
    species = [s for s in species if s[2] in present]
    out: dict[str, Any] = {
        "schema": SCHEMA_ID,
        "resis": list(panel.resis),
        "seq": panel.seq,
        "cons": [_r(x) for x in t["cons"]],
        "ident": [_r(x) for x in t["ident"]],
        "grade": [int(x) for x in t["grade"]],
        "cover": [_r(x) for x in t["cover"]],
        "n_orth": int(t["n_orth"]),
        "orth": [_orth_json(panel, m) for m in panel.orthologs],
        "par": [_par_json(panel, m) for m in panel.paralogs],
        "par_cons": [_r(x) for x in t["par_cons"]],
        "par_ident": [_r(x) for x in t["par_ident"]],
        "par_cover": [_r(x) for x in t["par_cover"]],
        "par_grade": [int(x) for x in t["par_grade"]],
        "n_par": int(t["n_par"]),
        "par_keys": list(t["par_keys"]),
        "par_note": par_method_note(panel),
        "species": species,
        "show": [k for k in spec.get("show", []) if k in present],
        "subject": spec.get("subject", ""),
        "kind": spec.get("kind", "host"),
        "note": (
            f"Public UniProt sequences from the offline reference cache ({t['n_orth']} orthologs, "
            f"{len(panel.paralogs)} paralogs of {panel.reference.get('gene', 'the target')}), each aligned globally to the "
            f"{'human ' if panel.reference.get('organism') == 'Homo sapiens' else ''}reference (BLOSUM62, gap open {-ALIGN_OPEN:g}, extend {-ALIGN_EXTEND:g}) and cut to the target residues. "
            "cons, ident, cover and grade come from the orthologs only; the par_* arrays describe the default paralogs "
            f"({', '.join(t['par_keys'])}) with low-confidence positions counted as gaps. Paralogs are matched by structure when "
            "AlphaFold DB models are cached, else by sequence."
        ),
        "generated": {
            "tool": "binderscope specificity.py",
            "version": __version__,
            "date": date or _dt.date.today().isoformat(),
            "params": panel.params,
            "warnings": list(panel.warnings),
        },
    }
    if panel.consurf:
        out["consurf"] = panel.consurf
    return out


# ---------------------------------------------------------------------------------------------
# documentation tables and command line
# ---------------------------------------------------------------------------------------------

def _pct(x: float, whole: bool = False) -> str:
    return f"{100 * x:.0f}" if whole else f"{100 * x:.1f}"


def render_doc_tables(panel: SpecificityPanel, hotspots_by_run: Mapping[str, Sequence[int]] | None = None) -> str:
    """Markdown tables of a panel: identities per member and, for each run in ``hotspots_by_run`` (UniProt numbers), the aligned residues at
    its hotspots. A project's documentation can embed them (``update_doc``)."""
    hs = hotspots_by_run or {}
    lo, hi = panel.resis[0], panel.resis[-1]
    obs = f"; observed = reference structure residues {panel.observed[0]} to {panel.observed[1]}" if panel.observed else ""
    cols = ["Member", "Accession", "Length", "Reviewed", "Alignment", "Full %", "Member length %"] + (["Observed %"] if panel.observed else []) + \
           ["Construct %", "Construct similarity %", "Construct low-confidence positions"]
    lines = [
        f"Identity is the share of reference residues with the same residue in the member (gaps count as different). "
        f"Construct = residues {lo} to {hi}{obs}. "
        "Member length is the identity over the member's own length, the figure refdata and UniProt report. Full, observed "
        "and member-length figures come from the sequence alignment; the construct figures and every per-residue cell below "
        "come from the Alignment column (paralogs are matched by structure when AlphaFold DB models are cached). Paralog "
        "identities are given to whole percent because they move by up to 5 points with the aligner settings.",
        "",
        "| " + " | ".join(cols) + " |",
        "|" + "---|" * len(cols),
    ]
    for m in panel.members:
        s = member_stats(panel, m)
        low = int(sum(1 for r in panel.resis if m.low[r - 1]))
        tag = m.key + (" (optional)" if m.optional else "")
        w = m.kind == "paralog"
        cells = [tag, m.acc, str(m.length), "yes" if m.reviewed else "no", m.aligner, _pct(s["full"]["identity"], w),
                 _pct(s["identity_of_member_length"], w)] + ([_pct(s["observed"]["identity"], w)] if panel.observed else []) + \
                [_pct(s["construct"]["identity"], w), _pct(s["construct"]["similarity"], w), f"{low} of {len(panel.resis)}"]
        lines.append("| " + " | ".join(cells) + " |")
    for run, nums in hs.items():
        nums = list(nums)
        labels = [panel.label_of(h) for h in nums]
        lines += ["", f"Hotspots of {run}: {', '.join(labels)}.", "",
                  "| Member | " + " | ".join(labels) + " | Identical | Similar | Differing |",
                  "|---|" + "---|" * len(labels) + "---|---|---|"]
        for m in panel.members:
            r = identity_at(panel, m.key, nums)
            cells = []
            for h in nums:
                a = m.full_map[h - 1]
                cells.append(("**" if a != panel.ref_seq[h - 1] else "") + a + ("**" if a != panel.ref_seq[h - 1] else "")
                             + ("*" if m.low[h - 1] else ""))
            lines.append(
                f"| {m.key}{' (optional)' if m.optional else ''} | " + " | ".join(cells)
                + f" | {r['identical']}/{r['n']} | {r['similar']} | {', '.join(r['differing']) or 'none'}"
                + (f"; gap {', '.join(r['gaps'])}" if r["gaps"] else "") + " |"
            )
    lines += ["", "Bold: differs from the reference. A trailing * marks low confidence: for a structural alignment a partner more than "
              f"{STRUCT_MAX_DIST:g} \u00c5 from the reference residue after superposition, for a sequence alignment a poorly "
              "matching window or a position next to a long insertion."]
    return "\n".join(lines)


_BEGIN = "<!-- BEGIN GENERATED: specificity.tables -->"
_END = "<!-- END GENERATED: specificity.tables -->"


def update_doc(path: Path, text: str) -> bool:
    """Replace the generated block of the doc. Returns True when the file changed."""
    src = path.read_text()
    new = re.sub(re.escape(_BEGIN) + r".*?" + re.escape(_END), lambda _m: f"{_BEGIN}\n{text}\n{_END}", src, flags=re.S)
    if new != src:
        path.write_text(new)
    return new != src


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m binder_review.specificity",
                                 description="Build the specificity panel of a target from the offline reference cache.")
    ap.add_argument("--uniprot", help="UniProt accession of the target (default: target.uniprot of the profile)")
    ap.add_argument("--residues", required=True, metavar="FIRST-LAST", help="target residues of the panel, UniProt numbers (the construct)")
    ap.add_argument("--observed", metavar="FIRST-LAST", help="residues the reference structure observes")
    ap.add_argument("--reference-id", help="PDB id of the reference structure (optional ConSurf-DB overlay)")
    ap.add_argument("--out", type=Path, help="write the HOM blob (homolog.json) here")
    ap.add_argument("--profile", type=Path, help="profile JSON whose specificity block names the panel")
    ap.add_argument("--hotspots", action="append", default=[], metavar="RUN=R1,R2", help="hotspots of a run for the tables (UniProt numbers)")
    ap.add_argument("--update-doc", type=Path, help="refresh the generated tables of this markdown file")
    ap.add_argument("--print-tables", action="store_true")
    args = ap.parse_args(argv)
    profile = json.loads(args.profile.read_text()) if args.profile else None
    lo, hi = (int(x) for x in args.residues.split("-"))
    obs = [int(x) for x in args.observed.split("-")] if args.observed else None
    hs = {k: [int(x) for x in v.split(",")] for k, v in (h.split("=", 1) for h in args.hotspots)}
    panel = build_panel(range(lo, hi + 1), profile, accession=args.uniprot, observed=obs, pdb=args.reference_id)
    for w in panel.warnings:
        print("warning:", w, file=sys.stderr)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(to_homolog_json(panel), separators=(",", ":")))
        print("wrote", args.out)
    if args.update_doc:
        print("doc updated" if update_doc(args.update_doc, render_doc_tables(panel, hs)) else "doc unchanged")
    if args.print_tables:
        print(render_doc_tables(panel, hs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
