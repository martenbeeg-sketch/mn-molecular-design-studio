"""Which cysteines of a binder form a disulfide, when the structure file does not say.

BindCraft 2 writes mmCIF files with no ``_struct_conn`` block, and the PDB files Binderscope places carry no
SSBOND records, so the file is silent about every bond. A predictor does not model covalent bonds either: it
places the two sulfur atoms of a disulfide close together without closing the bond to its exact length. Reading
"no record" as "no bond" calls the conserved framework disulfide of every nanobody an unpaired cysteine, so
the pairing is read from the geometry instead, and the result says that it was.

The two cut-offs, both on the SG-SG distance:

* A closed disulfide is about 2.03 A (the S-S bond; Pauling, and the crystal surveys of Richardson 1981 and
  Petersen 1999 put it at 2.0 to 2.1 A).
* Two sulfur atoms that do not interact sit at the sum of their van der Waals radii, 3.6 A (Bondi 1964), or
  farther. Contacts well inside that figure are rare without a bond holding them.

Checked against the predicted complexes of a scaffolded VHH campaign, where each VHH carries exactly the two
framework cysteines of the fold: every bonded sulfur pair fell between the bond length and van der Waals
contact, the great majority at 3.0 A or less and a thin tail above it. The 2.6 A limit the code used before
this module threw away bonded pairs, which is what prompted this. The physical limit and the observed
distribution agree on where a pair ends, so the cut-offs sit where each argument ends:

``CLOSE_MAX`` 3.0 A   0.6 A inside van der Waals contact. A pair this close is one the model drew together
                      on purpose, whatever bond length it settled on. Called "close".
``LOOSE_MAX`` 3.4 A   0.2 A inside van der Waals contact, so the thin tail of looser bonded pairs stays
                      inside the band while sulfurs that are only near each other stay out. Between 3.0 and
                      3.4 A the sulfurs are still compressed past contact, but a model that ends here is
                      less sure of itself, so the pair is kept and called "loose". Beyond it two cysteines
                      are merely near each other and are not paired.

A single cut-off would either throw away the loose pairs (turning a real disulfide into a false liability)
or hide how far a pair is from a closed bond. Two bands keep the pair and tell the reader which kind it is.

Pairing is not proximity. A sulfur takes part in one disulfide at most, so among the close pairs the
combination that uses each sulfur once is chosen, maximising the summed margin inside the loose cut-off
(LOOSE_MAX minus distance). That prefers more pairs and, between rival pairings, the tighter bonds: a 2.1 A
pair beats two 3.3 A pairs that would need the same sulfurs.

A disulfide recorded in the file always wins: inference is only the fallback for a file that records none.
"""
from __future__ import annotations

import itertools
import math
from typing import Iterable, Mapping, Optional, Sequence

CLOSE_MAX = 3.0
LOOSE_MAX = 3.4

BASIS_RECORDED = "recorded"      # the file lists the bonds (mmCIF _struct_conn, PDB SSBOND or LINK)
BASIS_INFERRED = "inferred"      # read from SG-SG geometry because the file records none
BASIS_NONE = "none"              # no record and no sulfur coordinates: nothing can be said

_EXACT_LIMIT = 14                # group size up to which every pairing is tried


def _dist(p: Sequence[float], q: Sequence[float]) -> float:
    return math.sqrt(sum((float(x) - float(y)) ** 2 for x, y in zip(p, q)))


def band_of(d: float) -> str:
    """"close" up to CLOSE_MAX, "loose" up to LOOSE_MAX; farther than that is not a pair."""
    return "close" if d <= CLOSE_MAX else "loose"


def nearest_sg(sg: Mapping[int, Sequence[float]]) -> dict[int, Optional[float]]:
    """Distance from each sulfur to the nearest other sulfur (None when it is the only one)."""
    out: dict[int, Optional[float]] = {}
    for p, xyz in sg.items():
        ds = [_dist(xyz, q) for r, q in sg.items() if r != p]
        out[int(p)] = min(ds) if ds else None
    return out


def _components(nodes: list[int], edges: dict[tuple[int, int], float]) -> list[list[int]]:
    """Groups of sulfurs joined by candidate pairs; a sulfur with no candidate is in no group."""
    adj: dict[int, set[int]] = {n: set() for n in nodes}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen: set[int] = set()
    out = []
    for n in nodes:
        if n in seen or not adj[n]:
            continue
        stack, comp = [n], []
        seen.add(n)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        out.append(sorted(comp))
    return out


def _best_matching(comp: list[int], edges: dict[tuple[int, int], float]) -> list[tuple[int, int]]:
    """The disjoint pairs inside one group with the largest summed margin (LOOSE_MAX minus distance)."""
    if len(comp) > _EXACT_LIMIT:
        # More sulfurs in one cluster than a protein has; take the closest pairs first.
        members = set(comp)
        used: set[int] = set()
        out = []
        for (a, b), _ in sorted(((e, v) for e, v in edges.items() if e[0] in members), key=lambda x: x[1]):
            if a not in used and b not in used:
                used.update((a, b))
                out.append((a, b))
        return out

    memo: dict[tuple[int, ...], tuple[float, list[tuple[int, int]]]] = {}

    def solve(nodes: tuple[int, ...]) -> tuple[float, list[tuple[int, int]]]:
        if len(nodes) < 2:
            return 0.0, []
        if nodes in memo:
            return memo[nodes]
        first, rest = nodes[0], nodes[1:]
        best: Optional[tuple[float, list[tuple[int, int]]]] = None
        for other in rest:                       # the pairing options come first so that a tie keeps the pair
            key = (first, other)
            if key not in edges:
                continue
            w, pairs = solve(tuple(x for x in rest if x != other))
            w += LOOSE_MAX - edges[key]
            if best is None or w > best[0]:
                best = (w, [key] + pairs)
        w, pairs = solve(rest)                   # or leave `first` unpaired
        if best is None or w > best[0]:
            best = (w, pairs)
        memo[nodes] = best
        return best

    return solve(tuple(comp))[1]


def infer_pairs(sg: Mapping[int, Sequence[float]]) -> list[dict]:
    """Disulfides read from the geometry of ``sg`` {residue: (x, y, z)} of SG atoms.

    Returns ``[{"a", "b", "d", "band"}]`` sorted by ``a`` with ``a < b``. A residue is in at most one pair.
    ``d`` is the SG-SG distance in A.
    """
    res = sorted(int(p) for p in sg)
    xyz = {int(p): v for p, v in sg.items()}
    edges: dict[tuple[int, int], float] = {}
    for a, b in itertools.combinations(res, 2):
        dd = _dist(xyz[a], xyz[b])
        if dd <= LOOSE_MAX:
            edges[(a, b)] = dd
    pairs: list[tuple[int, int]] = []
    for comp in _components(res, edges):
        pairs += _best_matching(comp, edges)
    return [{"a": a, "b": b, "d": round(edges[(a, b)], 2), "band": band_of(edges[(a, b)])} for a, b in sorted(pairs)]


def assess(sg: Mapping[int, Sequence[float]], recorded: Optional[Iterable[Sequence[int]]] = None) -> dict:
    """Disulfides of one binder and how they are known.

    ``sg``        {1-based residue: (x, y, z)} of the binder's cysteine SG atoms.
    ``recorded``  pairs the structure file lists as disulfides, 1-based. When there are any they are the answer,
                  untouched: a record is not overridden by geometry, even where the geometry disagrees.

    Returns ``{"basis", "pairs", "nearest"}``. ``basis`` is "recorded", "inferred" or "none" (no record and no
    sulfur coordinates, so nothing is judged). ``pairs`` is ``[{"a", "b", "d", "band"}]``; ``d`` is the SG-SG
    distance when both sulfurs are in ``sg`` (otherwise None) and ``band`` is "recorded", "close" or "loose".
    ``nearest`` is :func:`nearest_sg` of ``sg``.
    """
    sg = {int(p): xyz for p, xyz in sg.items()}
    nearest = nearest_sg(sg)
    rec = sorted({(min(int(a), int(b)), max(int(a), int(b))) for a, b in (recorded or []) if int(a) != int(b)})
    if rec:
        pairs = [{"a": a, "b": b, "d": round(_dist(sg[a], sg[b]), 2) if a in sg and b in sg else None, "band": "recorded"}
                 for a, b in rec]
        return {"basis": BASIS_RECORDED, "pairs": pairs, "nearest": nearest}
    if not sg:
        return {"basis": BASIS_NONE, "pairs": [], "nearest": nearest}
    return {"basis": BASIS_INFERRED, "pairs": infer_pairs(sg), "nearest": nearest}
