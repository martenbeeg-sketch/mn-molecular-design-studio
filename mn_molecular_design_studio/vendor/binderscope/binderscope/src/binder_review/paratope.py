"""Paratope analysis for campaigns built on a fixed scaffold (docs/paratope.md).

A nanobody, scFv, Fab or ankyrin-repeat campaign designs a few segments of one scaffold and leaves the rest alone.
Every design is then the same fold, so fold similarity measures the scaffold and not the designs. What varies is
which segments of the scaffold touch the target, how long the loops that were allowed to change length came out
and where on the target the binders land. This module answers the first two; ``families.py`` combines them with
the epitope spread into ``runs[].diversity_paratope`` (docs/families.md).

Where the regions come from. BindCraft 2 records which scaffold positions it was allowed to change
(``settings.filters.Scaffold_Sequence_Retained_Fraction.params.scaffold_edits``, surfaced by the reader as
``RunData.annotations["scaffold"]``) and which designed stretches may change length. That record is exact, so
there is no numbering by inference (no IMGT, no ANARCI, no HMMER): the designed segments are the stretches the
spec lists, and everything else is scaffold framework. Because the retained-sequence filter was mandatory at 1.0,
every position the spec does not list keeps the scaffold residue, so the invariant stretches between the designed
segments are exact anchors that fix each design's length offsets. A design is aligned by finding the one way to
choose the lengths of its variable stretches under which every invariant stretch matches; a design that matches
in no way, or in more than one, is reported as unaligned and never guessed at.

Grammar of the spec: ``<chain><start>[-<end>][*][(<lo>-<hi>)]``, comma separated. A bare or starred entry is a
fixed-length substitution; ``(lo-hi)`` marks a stretch whose designed length varies in that range. A token that
does not parse switches the analysis off for the run, with the reason in ``off_reason``.

Segment names come from the modality: CDR1 to CDR3 for a VHH, H1 to H3 and L1 to L3 for an scFv or Fab, and
"designed segment 1..n" for everything else (an ankyrin-repeat protein has no CDRs). Where the scaffold does not
fit the modality's expected pattern the module falls back to the numbered names and says so in ``notes``.

Everything here is computed from sequences and BindCraft's own ``Interface_Binder_Residues`` column (read by the
reader into ``DesignRecord.interface_residues``), so the contact criterion is BindCraft's, not Binderscope's
recomputed one. Region boundaries are the campaign's own edit spec, not IMGT or Kabat.
"""
from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

SCHEMA = 1

EDIT_TOKEN = re.compile(r"^([A-Za-z])(\d+)(?:-(\d+))?(\*)?(?:\((\d+)-(\d+)\))?$")
#: Designed stretches separated by at most this many undesigned residues are one segment (27-33 and 35 of a VHH
#: are one loop; 59 and 97 are not). A convention of this module, not a property of the scaffold.
MERGE_GAP = 1
#: VHH only: the former VL face of framework 2, in the scaffold's own numbering (Kabat-like 37/44/45/47).
FACE_POSITIONS = (37, 44, 45, 47)
FACE_REGION = "FR2-face"
#: How many framework positions are listed in a run block.
MAX_FRAMEWORK_POSITIONS = 20
#: How many unaligned designs are named in a run block.
MAX_UNALIGNED_LISTED = 10

NOT_SCAFFOLDED = None


class ParatopeError(ValueError):
    """The analysis cannot run for this run; the message says why in words a person can read."""


# --------------------------------------------------------------------------------------------- the spec
@dataclass(frozen=True)
class Edit:
    """One entry of the edit spec: scaffold positions start..end of a chain, substituted (``lo`` is None) or of variable length."""

    chain: str
    start: int
    end: int
    lo: Optional[int] = None
    hi: Optional[int] = None

    @property
    def variable(self) -> bool:
        return self.lo is not None

    @property
    def n(self) -> int:
        return self.end - self.start + 1


def parse_edits(text: Any) -> list[Edit]:
    """The edit spec as sorted ``Edit`` entries. Raises ``ParatopeError`` on any token that does not parse, on a range that
    runs backwards, on a length range with lo > hi and on entries that overlap."""
    if not isinstance(text, str) or not text.strip():
        raise ParatopeError("the scaffold edit spec is empty")
    out: list[Edit] = []
    for raw in re.sub(r"\s+", "", text).split(","):
        if not raw:
            continue
        m = EDIT_TOKEN.match(raw)
        if not m:
            raise ParatopeError(f"the scaffold edit spec has a token that is not in the form A12, A12-18, A12* or A12-18(5-7): {raw!r}")
        chain, a, b, _star, lo, hi = m.groups()
        start = int(a)
        end = int(b) if b is not None else start
        if end < start:
            raise ParatopeError(f"the scaffold edit spec has a range that runs backwards: {raw!r}")
        if lo is not None and int(lo) > int(hi):
            raise ParatopeError(f"the scaffold edit spec has a length range with its limits the wrong way round: {raw!r}")
        out.append(Edit(chain.upper(), start, end, int(lo) if lo is not None else None, int(hi) if hi is not None else None))
    if not out:
        raise ParatopeError("the scaffold edit spec lists no positions")
    out.sort(key=lambda e: (e.chain, e.start))
    for a, b in zip(out, out[1:]):
        if a.chain == b.chain and b.start <= a.end:
            raise ParatopeError(f"the scaffold edit spec lists overlapping entries: {a.chain}{a.start}-{a.end} and {b.chain}{b.start}-{b.end}")
    return out


# --------------------------------------------------------------------------------------------- the scaffold
def read_scaffold_file(path: str | Path) -> dict[str, dict[int, str]]:
    """{chain: {residue number: one-letter code}} of the polymer chains of a scaffold structure file (gemmi)."""
    import gemmi

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    if not len(st):
        raise ParatopeError(f"the scaffold file {Path(path).name} holds no model")
    out: dict[str, dict[int, str]] = {}
    for ch in st[0]:
        seq: dict[int, str] = {}
        for res in ch.get_polymer():
            info = gemmi.find_tabulated_residue(res.name)
            code = (info.one_letter_code or "").upper() if info else ""
            seq[res.seqid.num] = code if code.strip() else "X"
        if seq:
            out[ch.name.upper()] = seq
    if not out:
        raise ParatopeError(f"the scaffold file {Path(path).name} has no polymer chain")
    return out


def find_scaffold_file(recorded: Optional[str], run_dir: Optional[str | Path] = None, extra: Iterable[str | Path] = ()) -> Optional[Path]:
    """The scaffold structure file: the path the campaign recorded, else a file of the same name next to the campaign folder, in a
    ``scaffolds`` folder beside it or inside it, else in one of the ``extra`` folders. None when nothing is found: a campaign moved to
    another machine has usually left its scaffold behind."""
    if not recorded:
        return None
    p = Path(recorded)
    if p.is_file():
        return p
    name = p.name
    roots: list[Path] = [Path(x) for x in extra]
    if run_dir:
        rd = Path(run_dir)
        roots += [rd, rd / "scaffolds", rd.parent, rd.parent / "scaffolds"]
    for r in roots:
        if (r / name).is_file():
            return r / name
    return None


# --------------------------------------------------------------------------------------------- the layout
@dataclass
class Piece:
    """One stretch of the scaffold in order: ``inv`` (invariant, must match), ``fix`` (designed, fixed length) or ``var`` (designed,
    length lo..hi). ``seg`` is the index of the segment a designed piece belongs to (None for a designed position outside every
    segment, such as the VHH framework-2 face)."""

    kind: str
    chain: str
    start: int
    end: int
    seq: str = ""
    lo: int = 0
    hi: int = 0
    seg: Optional[int] = None


@dataclass
class Segment:
    """A designed segment: the span from its first to its last designed position, in scaffold numbering."""

    index: int
    chain: str
    first: int
    last: int
    designed: list[tuple[int, int]]
    variable: list[tuple[int, int, int, int]]          # (start, end, lo, hi) of its variable pieces
    name: str = ""

    @property
    def length(self) -> tuple[int, int]:
        """What a design's segment length is allowed to be: the variable stretch's range when it has one, else the span."""
        if self.variable:
            return (sum(v[2] for v in self.variable), sum(v[3] for v in self.variable))
        n = self.last - self.first + 1
        return (n, n)


@dataclass
class Layout:
    """The scaffold as the spec divides it: pieces per chain, segments, and the region of every scaffold position."""

    modality: str
    chains: list[str]
    scaffold: dict[str, dict[int, str]]
    pieces: dict[str, list[Piece]]
    segments: list[Segment]
    region_of: dict[tuple[str, int], str]               # (chain, scaffold position) -> region name
    framework_designed: list[tuple[str, int]]
    face: tuple[int, ...]
    naming: str                                         # "VHH" | "antibody" | "generic"
    scaffold_file: str = ""
    notes: list[str] = field(default_factory=list)

    @property
    def region_order(self) -> list[str]:
        """Segments in order, then the framework regions in the order of the scaffold."""
        names = [s.name for s in self.segments]
        for key in sorted(self.region_of):
            if self.region_of[key] not in names:
                names.append(self.region_of[key])
        return names


def _group_segments(edits: Sequence[Edit], chain: str, skip: set[int]) -> list[tuple[list[tuple[int, int]], list[tuple[int, int, int, int]]]]:
    """Designed stretches of one chain grouped into segments: [(designed ranges, variable pieces)], positions in ``skip`` left out."""
    units: list[tuple[int, int, Optional[tuple[int, int]]]] = []
    for e in edits:
        if e.chain != chain:
            continue
        if e.variable:
            units.append((e.start, e.end, (e.lo, e.hi)))
            continue
        run_start: Optional[int] = None
        for p in range(e.start, e.end + 2):
            ok = p <= e.end and p not in skip
            if ok and run_start is None:
                run_start = p
            elif not ok and run_start is not None:
                units.append((run_start, p - 1, None))
                run_start = None
    units.sort()
    groups: list[list[tuple[int, int, Optional[tuple[int, int]]]]] = []
    for u in units:
        if groups and u[0] - groups[-1][-1][1] - 1 <= MERGE_GAP:
            groups[-1].append(u)
        else:
            groups.append([u])
    out = []
    for g in groups:
        designed: list[tuple[int, int, Optional[tuple[int, int]]]] = []
        for a, b, var in g:
            if designed and var is None and designed[-1][2] is None and designed[-1][1] + 1 == a:
                designed[-1] = (designed[-1][0], b, None)
            else:
                designed.append((a, b, var))
        out.append(([(a, b) for a, b, _v in designed], [(u[0], u[1], u[2][0], u[2][1]) for u in g if u[2] is not None]))
    return out


def _heavy_or_light(seq: dict[int, str], after: int) -> Optional[str]:
    """H if the scaffold residues just after position ``after`` read like a heavy-chain J region (WGxG), L if like a light-chain one (FGxG)."""
    text = "".join(seq.get(p, "X") for p in range(after + 1, after + 9))
    if re.search(r"WG.G", text):
        return "H"
    if re.search(r"FG.G", text):
        return "L"
    return None


def _build(edits: Sequence[Edit], scaffold: dict[str, dict[int, str]], modality: str, scaffold_file: str, face_rule: bool) -> Layout:
    chains = sorted({e.chain for e in edits})
    for c in chains:
        if c not in scaffold:
            raise ParatopeError(f"the edit spec names scaffold chain {c}, which the scaffold file does not have")
        nums = sorted(scaffold[c])
        if nums != list(range(nums[0], nums[-1] + 1)):
            raise ParatopeError(f"scaffold chain {c} is not numbered consecutively, so the edit spec cannot be placed on it")
        for e in edits:
            if e.chain == c and (e.start < nums[0] or e.end > nums[-1]):
                raise ParatopeError(f"the edit spec lists {c}{e.start}-{e.end}, outside the scaffold chain ({c}{nums[0]}-{nums[-1]})")
    face: tuple[int, ...] = ()
    skip: set[int] = set()
    if face_rule:
        face = tuple(p for p in FACE_POSITIONS)
        skip = set(face)
    # pieces, in scaffold order, with the segment of each designed piece
    segs: list[Segment] = []
    pieces: dict[str, list[Piece]] = {}
    for c in chains:
        groups = _group_segments(edits, c, skip)
        first, last = min(scaffold[c]), max(scaffold[c])
        seg_of_pos: dict[int, int] = {}
        for ranges, var in groups:
            idx = len(segs)
            segs.append(Segment(idx, c, ranges[0][0], ranges[-1][1], ranges, var))
            for a, b in ranges:
                for p in range(a, b + 1):
                    seg_of_pos[p] = idx
        # a face position inside a segment span would make the face rule ambiguous
        lst: list[Piece] = []
        pos = first
        for e in (x for x in edits if x.chain == c):
            if e.start > pos:
                lst.append(Piece("inv", c, pos, e.start - 1, "".join(scaffold[c][p] for p in range(pos, e.start))))
            if e.variable:
                lst.append(Piece("var", c, e.start, e.end, lo=e.lo, hi=e.hi, seg=seg_of_pos.get(e.start)))
            else:
                # a fixed entry may straddle the face rule's skipped positions: one piece per maximal run of the same segment
                a = e.start
                while a <= e.end:
                    s = seg_of_pos.get(a)
                    b = a
                    while b + 1 <= e.end and seg_of_pos.get(b + 1) == s:
                        b += 1
                    lst.append(Piece("fix", c, a, b, lo=b - a + 1, hi=b - a + 1, seg=s))
                    a = b + 1
            pos = e.end + 1
        if pos <= last:
            lst.append(Piece("inv", c, pos, last, "".join(scaffold[c][p] for p in range(pos, last + 1))))
        pieces[c] = lst
    # a segment span must not swallow a face position
    for s in segs:
        for p in face:
            if s.first <= p <= s.last and s.chain == chains[0]:
                raise ParatopeError("face positions fall inside a designed segment")
    region_of, naming, notes = _name_regions(segs, scaffold, chains, modality, face)
    fw_designed = []
    for c in chains:
        for pc in pieces[c]:
            if pc.kind != "inv" and pc.seg is None:
                fw_designed += [(c, p) for p in range(pc.start, pc.end + 1)]
    return Layout(modality, chains, scaffold, pieces, segs, region_of, fw_designed, face if naming == "VHH" else (), naming, scaffold_file, notes)


def _name_regions(segs: list[Segment], scaffold: dict[str, dict[int, str]], chains: list[str], modality: str,
                  face: tuple[int, ...]) -> tuple[dict[tuple[str, int], str], str, list[str]]:
    """Name the segments and the framework between them from the modality. Returns (region of every scaffold position, how it was
    named, notes); the names are written into ``segs``."""
    notes: list[str] = []
    region: dict[tuple[str, int], str] = {}
    mod = (modality or "").lower()

    def generic() -> tuple[dict[tuple[str, int], str], str, list[str]]:
        for i, s in enumerate(segs, 1):
            s.name = f"designed segment {i}"
        for c in chains:
            for p in scaffold[c]:
                inside = next((s for s in segs if s.chain == c and s.first <= p <= s.last), None)
                region[(c, p)] = inside.name if inside else "framework"
        return region, "generic", notes

    if mod == "vhh" and len(chains) == 1 and len(segs) == 3:
        c = chains[0]
        sc = scaffold[c]
        for i, s in enumerate(segs, 1):
            s.name = f"CDR{i}"
        for p in sc:
            inside = next((s for s in segs if s.first <= p <= s.last), None)
            if inside:
                region[(c, p)] = inside.name
            elif p in face:
                region[(c, p)] = FACE_REGION
            elif p < segs[0].first:
                region[(c, p)] = "FR1"
            elif p < segs[1].first:
                region[(c, p)] = "FR2"
            elif p < segs[2].first:
                region[(c, p)] = "FR3"
            else:
                region[(c, p)] = "FR4"
        return region, "VHH", notes
    if mod in ("scfv", "fab") and segs:
        per_chain: dict[str, list[Segment]] = defaultdict(list)
        for s in segs:
            per_chain[s.chain].append(s)
        if all(len(v) % 3 == 0 for v in per_chain.values()):
            domains: list[tuple[str, list[Segment], str]] = []
            ok = True
            for c in chains:
                for k in range(0, len(per_chain[c]), 3):
                    trio = per_chain[c][k:k + 3]
                    hl = _heavy_or_light(scaffold[c], trio[2].last)
                    if hl is None:
                        ok = False
                    domains.append((c, trio, hl or "?"))
            if ok and len({d[2] for d in domains}) == len(domains):
                for c, trio, hl in domains:
                    for i, s in enumerate(trio, 1):
                        s.name = f"{hl}{i}"
                for c in chains:
                    doms = sorted((d for d in domains if d[0] == c), key=lambda d: d[1][0].first)
                    for p in scaffold[c]:
                        inside = next((s for s in segs if s.chain == c and s.first <= p <= s.last), None)
                        if inside:
                            region[(c, p)] = inside.name
                            continue
                        name = "inter-domain"
                        for di, (_c, trio, hl) in enumerate(doms):
                            if p < trio[0].first:
                                name = f"{hl}-FR1" if di == 0 else "inter-domain"
                            elif p < trio[1].first:
                                name = f"{hl}-FR2"
                            elif p < trio[2].first:
                                name = f"{hl}-FR3"
                            elif di == len(doms) - 1:
                                name = f"{hl}-FR4"
                            else:
                                continue
                            break
                        region[(c, p)] = name
                return region, "antibody", notes
        notes.append(f"The {modality} scaffold's designed segments do not fall into heavy and light variable domains of three segments each "
                     "(or the domains cannot be told apart from the scaffold sequence), so the segments are numbered, not named H1 to L3.")
        return generic()
    if mod == "vhh":
        notes.append(f"A VHH has three CDRs but the edit spec lists {len(segs)} designed segments, so the segments are numbered, not named "
                     "CDR1 to CDR3.")
    return generic()


def build_layout(edits: Sequence[Edit], scaffold: dict[str, dict[int, str]], modality: str, scaffold_file: str = "") -> Layout:
    """The scaffold divided as the edit spec says. For a VHH the framework-2 face (37, 44, 45, 47) is set apart from the CDRs when the scaffold
    carries the landmarks that make those numbers mean what they say (Cys22 and Trp36)."""
    chains = sorted({e.chain for e in edits})
    can_face = (modality or "").lower() == "vhh" and len(chains) == 1 and scaffold.get(chains[0], {}).get(36) == "W" \
        and scaffold.get(chains[0], {}).get(22) == "C"
    if can_face:
        try:
            lay = _build(edits, scaffold, modality, scaffold_file, True)
            if lay.naming == "VHH":
                return lay
        except ParatopeError:
            pass
    lay = _build(edits, scaffold, modality, scaffold_file, False)
    if (modality or "").lower() == "vhh" and len(chains) == 1 and not can_face:
        lay.notes.append("The scaffold lacks the landmarks (Cys22, Trp36) that identify the framework-2 face, so no face is reported.")
    return lay


# --------------------------------------------------------------------------------------------- alignment
def align_chain(pieces: Sequence[Piece], seq: str) -> list[tuple[Piece, int, int]]:
    """Place the pieces of one chain on a design sequence: [(piece, start index in seq, length)]. Raises ``ParatopeError`` when the
    invariant stretches match at no choice of the variable lengths, or at more than one."""
    n = len(pieces)
    lo_suffix = [0] * (n + 1)
    hi_suffix = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        p = pieces[i]
        a, b = (len(p.seq), len(p.seq)) if p.kind == "inv" else (p.lo, p.hi)
        lo_suffix[i] = lo_suffix[i + 1] + a
        hi_suffix[i] = hi_suffix[i + 1] + b
    if not (lo_suffix[0] <= len(seq) <= hi_suffix[0]):
        raise ParatopeError(f"length {len(seq)} is outside the {lo_suffix[0]} to {hi_suffix[0]} the edit spec allows")
    sols: list[list[tuple[Piece, int, int]]] = []
    acc: list[tuple[Piece, int, int]] = []

    def rec(i: int, pos: int) -> None:
        if len(sols) > 1:
            return
        if i == n:
            if pos == len(seq):
                sols.append(list(acc))
            return
        if pos + lo_suffix[i] > len(seq) or pos + hi_suffix[i] < len(seq):
            return
        p = pieces[i]
        if p.kind == "inv":
            if seq.startswith(p.seq, pos):
                acc.append((p, pos, len(p.seq)))
                rec(i + 1, pos + len(p.seq))
                acc.pop()
        else:
            for length in range(p.lo, p.hi + 1):
                if pos + length > len(seq):
                    break
                acc.append((p, pos, length))
                rec(i + 1, pos + length)
                acc.pop()

    rec(0, 0)
    if not sols:
        raise ParatopeError("the scaffold's invariant stretches do not match at any length the edit spec allows")
    if len(sols) > 1:
        raise ParatopeError("more than one placement of the variable-length segments fits the scaffold's invariant stretches")
    return sols[0]


# --------------------------------------------------------------------------------------------- per design
@dataclass
class DesignPara:
    """One design placed on the scaffold, with where its interface residues fall."""

    id: str
    regions: dict[str, int]                      # region -> interface residues (zero entries left out)
    loops: list[str]                             # designed segments with at least one interface residue, in scaffold order
    framework: int                               # interface residues outside every designed segment
    lengths: dict[str, int]                      # segment -> length of its variable stretch in this design
    face: dict[str, str]                         # designed framework-2 face positions -> the residue in this design
    face_mutated: list[str]                      # those positions whose residue differs from the scaffold's
    framework_at: list[tuple[str, int]]          # (region, scaffold position) of each framework interface residue, as a list
    n_interface: int

    def to_blob(self, layout: Layout) -> dict[str, Any]:
        out: dict[str, Any] = {"segments": dict(self.regions), "loops": list(self.loops), "framework": self.framework}
        if self.lengths:
            out["lengths"] = dict(self.lengths)
        if layout.face:
            out["face"] = dict(self.face)
            out["face_mutated"] = list(self.face_mutated)
        if self.framework_at:
            out["framework_at"] = [[r, p] for r, p in self.framework_at]
        return out


def _chain_letters(runs_residues: Iterable[Sequence[str]]) -> list[str]:
    return sorted({t[0].upper() for res in runs_residues for t in res if t})


def place_design(layout: Layout, design_id: str, sequence: str, interface: Sequence[str],
                 chain_map: Mapping[str, int]) -> DesignPara:
    """Align one design (``sequence`` is the binder chains joined by '/') and tally where its interface residues fall.
    ``interface`` is ["B13", ...] (chain letter and 1-based residue number in the chain), ``chain_map`` maps the letters to the
    index of the chain in ``sequence``. Raises ``ParatopeError`` for a design that does not align."""
    parts = sequence.split("/")
    if len(parts) != len(layout.chains):
        raise ParatopeError(f"the design has {len(parts)} binder chain(s) but the edit spec covers {len(layout.chains)}")
    # residue index (1-based, per chain) -> (region, scaffold position or None); designed-variable residues have no scaffold position
    maps: list[list[tuple[str, Optional[int]]]] = []
    lengths: dict[str, int] = {}
    scaf_at: list[dict[int, int]] = []                  # scaffold position -> 0-based index in the design chain
    for ci, c in enumerate(layout.chains):
        placed = align_chain(layout.pieces[c], parts[ci])
        m: list[tuple[str, Optional[int]]] = [("", None)] * len(parts[ci])
        at: dict[int, int] = {}
        for piece, start, length in placed:
            if piece.kind == "var":
                seg = layout.segments[piece.seg] if piece.seg is not None else None
                name = seg.name if seg else "framework"
                for k in range(length):
                    m[start + k] = (name, None)
                if seg:
                    lengths[seg.name] = lengths.get(seg.name, 0) + length
            else:
                for k, sp in enumerate(range(piece.start, piece.end + 1)):
                    m[start + k] = (layout.region_of.get((c, sp), "framework"), sp)
                    at[sp] = start + k
        maps.append(m)
        scaf_at.append(at)
    regions: Counter[str] = Counter()
    fw_at: list[tuple[str, int]] = []
    seg_names = {s.name for s in layout.segments}
    n_if = 0
    for tok in interface:
        tok = tok.strip()
        if not tok:
            continue
        letter, num = tok[0].upper(), tok[1:]
        if not num.isdigit() or letter not in chain_map:
            raise ParatopeError(f"interface residue {tok!r} does not name a binder chain and a number")
        ci = chain_map[letter]
        k = int(num) - 1
        if not (0 <= k < len(maps[ci])):
            raise ParatopeError(f"interface residue {tok!r} lies outside the sequence")
        region, sp = maps[ci][k]
        regions[region] += 1
        n_if += 1
        if region not in seg_names and sp is not None:
            fw_at.append((region, sp))
    loops = [s.name for s in layout.segments if regions.get(s.name)]
    framework = sum(v for r, v in regions.items() if r not in seg_names)
    face: dict[str, str] = {}
    mutated: list[str] = []
    if layout.face:
        c0 = layout.chains[0]
        designed = {p for (c, p) in layout.framework_designed if c == c0}
        for p in layout.face:
            if p in designed and p in scaf_at[0]:
                aa = parts[0][scaf_at[0][p]]
                face[str(p)] = aa
                if aa != layout.scaffold[c0].get(p):
                    mutated.append(str(p))
    order = {r: i for i, r in enumerate(layout.region_order)}
    return DesignPara(design_id, {r: regions[r] for r in sorted(regions, key=lambda r: order.get(r, 10 ** 6))}, loops, framework, lengths,
                      face, mutated, sorted(fw_at, key=lambda t: (t[1], t[0])), n_if)


# --------------------------------------------------------------------------------------------- the run
@dataclass
class RunParatope:
    """The paratope analysis of one scaffolded run."""

    run_id: str
    layout: Layout
    designs: dict[str, DesignPara]                    # aligned designs
    unaligned: dict[str, str]                         # design id -> why
    accepted: list[str]                               # ids of accepted designs (aligned or not), in run order
    rejected: list[str]
    spec_text: str
    retained: dict[str, Any]
    notes: list[str] = field(default_factory=list)

    # -- groups
    def group(self, status: str) -> list[DesignPara]:
        ids = self.accepted if status == "accepted" else self.rejected
        return [self.designs[i] for i in ids if i in self.designs]

    def segments_blob(self) -> list[dict[str, Any]]:
        out = []
        for s in self.layout.segments:
            lo, hi = s.length
            d: dict[str, Any] = {"name": s.name, "chain": s.chain, "scaffold": [s.first, s.last], "designed": [[a, b] for a, b in s.designed],
                                 "length": [lo, hi]}
            if s.variable:
                d["variable"] = [[a, b, l, h] for a, b, l, h in s.variable]
            out.append(d)
        return out

    def block(self, status: str = "accepted") -> Optional[dict[str, Any]]:
        """``runs[].paratope`` (status "accepted") or ``runs[].paratope_rejected`` ("rejected"): None when no design of that group aligned."""
        grp = self.group(status)
        if not grp:
            return None
        lay = self.layout
        seg_names = [s.name for s in lay.segments]
        total = sum(d.n_interface for d in grp)
        tally: Counter[str] = Counter()
        for d in grp:
            tally.update(d.regions)
        order = {r: i for i, r in enumerate(lay.region_order)}
        ranked = sorted(tally, key=lambda r: (-tally[r], order.get(r, 10 ** 6)))
        share = {r: round(100.0 * tally[r] / total, 1) for r in ranked} if total else {}
        fw_total = sum(v for r, v in tally.items() if r not in seg_names)
        uses: Counter[tuple[str, ...]] = Counter(tuple(d.loops) for d in grp)
        loop_use = [{"loops": list(k), "n": n, "share": round(100.0 * n / len(grp), 1)}
                    for k, n in sorted(uses.items(), key=lambda kv: (-kv[1], -len(kv[0]), kv[0]))]
        hits: Counter[tuple[str, int]] = Counter()
        for d in grp:
            hits.update(set(d.framework_at))
        fw_pos = []
        for (reg, sp), n in sorted(hits.items(), key=lambda kv: (-kv[1], kv[0][1], kv[0][0]))[:MAX_FRAMEWORK_POSITIONS]:
            aa = lay.scaffold[lay.chains[0]].get(sp, "X") if len(lay.chains) == 1 else "X"
            fw_pos.append({"region": reg, "scaffold": sp, "aa": aa, "n": n})
        group_word = "accepted" if status == "accepted" else "rejected"
        block: dict[str, Any] = {
            "source": "scaffold_edits", "scaffold": lay.scaffold_file, "modality": lay.modality, "naming": lay.naming,
            "segments": self.segments_blob(), "framework_designed": [p for _c, p in lay.framework_designed],
            "aligned": len(self.designs), "unaligned": len(self.unaligned),
            "n": len(grp), "group": group_word,
            "share": share, "framework_share": round(100.0 * fw_total / total, 1) if total else 0.0,
            "loop_use": loop_use, "framework_positions": fw_pos,
            "note": f"Shares are of interface residues on {group_word} designs. Regions are this campaign's own edit spec, not IMGT or Kabat; "
                    "the contact criterion is BindCraft's Interface_Binder_Residues.",
        }
        if len(lay.chains) > 1:
            block["chains"] = list(lay.chains)
        if self.unaligned:
            block["unaligned_designs"] = [{"id": i, "why": w} for i, w in list(self.unaligned.items())[:MAX_UNALIGNED_LISTED]]
        no_if = sum(1 for d in grp if d.n_interface == 0)
        if no_if:
            block["no_interface"] = no_if
        if lay.notes or self.notes:
            block["notes"] = list(dict.fromkeys([*lay.notes, *self.notes]))
        return block

    def design_blobs(self) -> dict[str, dict[str, Any]]:
        return {i: d.to_blob(self.layout) for i, d in self.designs.items()}

    # -- what the diversity axes need
    def family_lengths(self, segment: str, family_of: Mapping[str, str]) -> dict[str, list[int]]:
        """{family id: sorted distinct lengths of ``segment``'s variable stretch among its accepted designs}."""
        out: dict[str, set[int]] = defaultdict(set)
        for i in self.accepted:
            d = self.designs.get(i)
            if d is not None and segment in d.lengths and i in family_of:
                out[family_of[i]].add(d.lengths[segment])
        return {f: sorted(v) for f, v in out.items()}


def has_scaffold(run: Any) -> bool:
    """True when the run recorded a scaffold edit spec (``RunData.annotations["scaffold"]``): a scaffolded run, whatever else is true of it."""
    sc = (getattr(run, "annotations", None) or {}).get("scaffold")
    return bool(isinstance(sc, Mapping) and sc.get("scaffold_edits"))


def off_reason(run: Any) -> Optional[str]:
    """Why the paratope analysis cannot run on a scaffolded run, from what the record says alone (None when it may be tried)."""
    sc = (getattr(run, "annotations", None) or {}).get("scaffold") or {}
    try:
        parse_edits(sc.get("scaffold_edits"))
    except ParatopeError as ex:
        return str(ex)
    thr = sc.get("threshold")
    try:
        short = thr is None or float(thr) < 1.0
    except (TypeError, ValueError):
        short = True
    if sc.get("mandatory") is False or short or sc.get("higher") is False:
        return ("the campaign did not require Scaffold_Sequence_Retained_Fraction to be 1.0 as a mandatory filter, so positions outside the edit "
                "spec may have changed and the scaffold's invariant stretches cannot anchor an exact alignment")
    return None


def analyse_run(run: Any, *, scaffold_dirs: Iterable[str | Path] = ()) -> tuple[Optional[RunParatope], Optional[str]]:
    """The paratope analysis of a run: (result, None) when it ran, (None, reason) when the run is scaffolded but the analysis is off, and
    (None, None) for a run that is not scaffolded. Never raises for a problem with the data: that is what the reason is for."""
    if not has_scaffold(run):
        return None, None
    sc = run.annotations["scaffold"]
    why = off_reason(run)
    if why:
        return None, why
    try:
        edits = parse_edits(sc["scaffold_edits"])
        path = find_scaffold_file(sc.get("scaffold"), run.source_dir, scaffold_dirs)
        label = Path(sc["scaffold"]).name if sc.get("scaffold") else ""
        if path is None:
            raise ParatopeError(f"the scaffold file {label or '(not recorded)'} is not on this machine, and its sequence is what anchors the "
                                "alignment of every design")
        modality = str((run.annotations or {}).get("modality") or "")
        layout = build_layout(edits, read_scaffold_file(path), modality, label)
        designs = list(run.designs)
        letters = _chain_letters(d.interface_residues for d in designs)
        if len(layout.chains) == 1:
            chain_map = {c: 0 for c in letters}
        else:
            if len(letters) != len(layout.chains):
                raise ParatopeError(f"the designs' interface residues name {len(letters)} binder chain(s) but the edit spec covers {len(layout.chains)}, "
                                    "so the chains cannot be matched to the scaffold's")
            chain_map = {c: i for i, c in enumerate(letters)}
        placed: dict[str, DesignPara] = {}
        unaligned: dict[str, str] = {}
        for d in designs:
            try:
                placed[d.id] = place_design(layout, d.id, d.sequence, d.interface_residues, chain_map)
            except ParatopeError as ex:
                unaligned[d.id] = str(ex)
        if not placed:
            raise ParatopeError("no design aligns to the scaffold" + (f" (for example {next(iter(unaligned.values()))})" if unaligned else ""))
    except ParatopeError as ex:
        return None, str(ex)
    except Exception as ex:                                    # a damaged scaffold file or similar: off, with the reason
        return None, f"the scaffold could not be read ({type(ex).__name__}: {ex})"
    notes = []
    if len(layout.chains) > 1:
        notes.append("Several binder chains: the scaffold's chains are matched to the design's chains in order, which no real campaign has checked.")
    res = RunParatope(run.run_id, layout, placed, unaligned, [d.id for d in run.accepted()], [d.id for d in run.rejected()],
                      sc["scaffold_edits"], {"threshold": sc.get("threshold"), "mandatory": sc.get("mandatory")}, notes)
    return res, None
