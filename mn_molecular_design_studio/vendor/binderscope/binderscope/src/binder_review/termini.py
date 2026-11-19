"""Terminus analysis: which end of a designed binder can carry a tag or a fusion, and what it would cost (docs/termini.md).

A tag decision is not one number, and not one tag. The two people usually choose between are His6 (6 residues, unstructured, at either end) and
Twin-Strep (28 residues, also unstructured but several times the bulk, captured by the tetrameric Strep-Tactin), and a terminus that carries
the first may not carry the second. Distance to the target is necessary but not sufficient: an end that is far in angstrom can still point into
the binder's own fold, into the target or into a partner subunit. So every end is judged for four tags. For the two unstructured ones the
question is whether a flexible chain of their length can route from the end into open solvent without crossing the target or the binder, which
is a much weaker condition than an empty cone; for the two folded fusions (a GFP-like barrel, an Fc) it is whether a sphere of their size fits
against the end in the cone the chain can leave in, with the linker each tag size needs. The evidence behind the answer (distance, direction,
secondary structure, burial, confidence) is kept next to it so that a person can disagree.

What it measures, per end (N and C):

* ``d``       terminal CA to the nearest target heavy atom (the number the interface module has always exported as nterm_d and cterm_d).
* ``aim``     the angle between the direction the chain leaves the fold in and the direction of the nearest target atom (0 = pointing at it).
* ``away``    where the end sits relative to the paratope: +1 on the far side of the binder, -1 on the paratope side. Immobilisation by a tag
              at an end with a negative value holds the binder by its binding face.
* ``contact`` / ``iface_n``  whether the terminal residue (or how many of the three nearest the end) is itself in the interface.
* ``rsa``, ``plddt``, ``ss``, ``elem``, ``tail``  burial, confidence and secondary structure of the terminal residues.
* ``room``    the radius of the largest sphere (a folded fusion) that fits against the end without a linker, and ``limit``, what stops the first tag that is not free.
* ``g``       per tag: ``free`` (fits as it is), ``linker`` (fits after a linker of ``l`` residues) or ``blocked`` (no way out or no room within the
              longest linker tried; for a folded fusion also when the end is itself in the interface). A flexible tag carries ``esc``, its path to open solvent.

Everything is a prediction made from one predicted conformation. The probes are a reach search and a sphere in a cone, not a simulation of a tag; the sizes
are round figures for the tag families, not measurements (see ``TAGS``). A flexible tail's position is the least reliable part of a model,
so a disordered tail is judged from the last ordered residue and says so (``conf``).

The geometry is rigid-body invariant, so it is computed in the design's own frame. The placed pose enters only through ``env``: a function
that answers, for a set of points of the design's frame, how close they come to what lies outside the construct (the other subunit, the rest
of the chain, the glycans), as the target frame sees them (``target_frame.TargetFrame.context``). Without it that part of the check is
off and says so; it is never filled with a default.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Mapping, Optional, Sequence

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

SCHEMA = 1

#: Tag families as (key, what it stands for, size in angstrom). Two kinds are tested differently, because they behave differently:
#:
#: * ``flexible`` tags (His6, Twin-Strep) are unstructured peptides. The question is not whether a corridor of their full length is empty, but
#:   whether a flexible chain of that many residues can get from the end into open solvent without crossing the target or the binder, which is a
#:   much weaker condition: the chain routes around obstacles. ``FLEXIBLE`` gives, per tag, its length in residues and how many of them (the
#:   motif: the histidines, one Strep-tag II unit) have to lie in open solvent to be bound by Ni-NTA or Strep-Tactin; the rest is reach.
#:   The size is the clearance in angstrom, from the nearest heavy atom, that counts as open solvent for that tag: a Twin-Strep is captured by
#:   Strep-Tactin XT, a tetramer, so it needs more open space around the motif than a His6 does.
#: * ``globular`` tags (a GFP, SNAP or mCherry barrel, an Fc or MBP) are folded and cannot route around anything: the size is the radius of
#:   the sphere that stands for them, tested in a cone against the surface (the probe below).
#:
#: Round figures and a judgement, not measurements: His6 is 6 residues (about 0.8 kDa); Twin-Strep is WSHPQFEK-(GGGS)3-WSHPQFEK, 28 residues
#: (about 3 kDa), commonly used at the C-terminus and at both ends; GFP-like barrels are about 24 x 42 A and can turn, so 16 A is an effective
#: radius; an Fc is the largest thing people fuse. Changing them changes verdicts, so they live here and are written into the data. The tag
#: sequences are not carried anywhere: a tag is not part of the designed sequence, and its residues (two tryptophans in a Twin-Strep) never
#: enter the design's liabilities.
TAGS = (
    ("his6", "His6 (6 residues, IMAC)", 5.0),
    ("strep", "Twin-Strep (28 residues, Strep-Tactin)", 10.0),
    ("domain", "GFP, SNAP, mCherry", 16.0),
    ("large", "Fc, MBP", 24.0),
)
FLEXIBLE = {"his6": (6, 3), "strep": (28, 8)}          # key -> (residues, residues of motif that must lie in open solvent)
TAG_KEYS = tuple(t[0] for t in TAGS)
#: how far the chain's exit direction may be bent by a flexible linker, in degrees
CONE_DEG = 40.0
CONE_RINGS = 2           # rings of directions between the axis and the cone edge
CONE_AZIMUTHS = 12
#: gap between the tip and the tag's surface, in angstrom (about one linker residue)
STANDOFF = 3.0
#: reach of one linker residue, angstrom (an extended chain is 3.3-3.8 A per residue; 3.3 is the honest, shorter figure)
REACH = 3.3
#: linker lengths tried, in residues; 15 is (G4S)3. An end that needs more is called blocked.
LINKER_STEPS = (0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15)
#: heavy-atom radius and the overlap a tag's envelope may have with the surface (a tag is not a rigid ball: it packs and gives)
ATOM_R = 1.7
SLACK = 2.0
#: a linker's own thickness: no heavy atom nearer than this to the path from the tip to the tag
PATH_CLEAR = 3.0
PATH_STEP = 1.0
#: the exit direction is the vector from the mean of this many CA positions along the chain to the end
EXIT_WINDOW = 6
#: a terminal stretch outside helix and strand is a tail from this length, and is read as disordered below this pLDDT
TAIL_MIN = 3
DISORDER_PLDDT = 70.0
TAIL_CREDIT_MAX = 10
#: the largest sphere radius the room search tries, and so the largest room it reports
ROOM_CAP = 30.0
#: terminal-residue summaries
RSA_N = 3
PLDDT_N = 5
#: "an end within this distance of the target is near it": the figure BindCraft's own termini_together uses for the two ends
NEAR = 10.0
#: which side of the binder an end is on: the paratope side when it is within NEAR of the target or lies toward the interface (``away`` below
#: PARA_AWAY), the far side when it lies away from it (``away`` at least FAR_AWAY) and clear of the target, else lateral
PARA_AWAY = -0.25
FAR_AWAY = 0.25
#: environment candidates tried per linker length before the next length is tried
ENV_TRIES = 4

#: the flexible-chain search: open solvent is looked for on a grid of this spacing within this distance of the end, through voxels at least
#: PATH_CLEAR from every atom, starting from the voxels within FLEX_START of the end that are not inside an atom (the chain leaves the end itself, wherever the model has
#: put it). The environment (other subunits, the rest of the protein) is asked about the best path found and blocks the voxels it crosses.
FLEX_SPACING = 3.0
FLEX_HALF = 60.0
FLEX_START = 4.0
#: within FLEX_START a voxel is open from this clearance (the end sits closer to its own neighbours than a chain may pass anywhere else), but never inside an atom
FLEX_START_CLEAR = 1.5
FLEX_ENV_ROUNDS = 6

STATUS_RANK = {"free": 0, "linker": 1, "blocked": 2}


# ---------------------------------------------------------------------------------------------------------------------------------
# format: what kind of molecule this is, and whether the question has an answer for it
# ---------------------------------------------------------------------------------------------------------------------------------
@dataclass
class Format:
    """What the run's molecule is, read from its own settings (never assumed to be a de novo binder)."""

    kind: str = "de novo"             # de novo, VHH, scFv, Fab, ARP, cyclic peptide or the raw word
    known: bool = False               # the run recorded its modality (False: read as de novo because nothing said otherwise)
    copies: int = 1                   # protomers of a homo-oligomer binder
    chains: int = 1                   # binder chains the run designed
    applies: bool = True              # False: the question has no answer for this format
    why: str = ""                     # the reason, when it does not
    convention: Optional[str] = None  # the end the format's usual practice tags ("C" for a VHH), with the reason in ``convention_why``
    convention_why: str = ""
    scaffold: bool = False
    peptide: bool = False             # short enough that a tag is a large part of the molecule
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind, "known": self.known, "copies": self.copies, "chains": self.chains, "applies": self.applies,
                "why": self.why, "convention": self.convention, "convention_why": self.convention_why, "scaffold": self.scaffold,
                "peptide": self.peptide, "notes": list(self.notes)}


_KIND = (("cyclic peptide", re.compile(r"cyc", re.I)), ("VHH", re.compile(r"vhh|nanobod|sdab", re.I)), ("scFv", re.compile(r"scfv", re.I)),
         ("Fab", re.compile(r"\bfab\b", re.I)), ("ARP", re.compile(r"\barp\b|ankyrin|darpin", re.I)),
         ("de novo", re.compile(r"de ?novo|^binder$|minibinder|mini-binder|^protein$", re.I)))
#: lengths up to this are called peptides for the one sentence that says a tag is then a large share of the molecule
PEPTIDE_MAX = 30


def _dig(raw: Any, *keys: str) -> Any:
    """First value found under any of ``keys`` in a settings dict, looking at the top level and under "settings", "resolved" and "metadata"."""
    if not isinstance(raw, Mapping):
        return None
    for scope in (raw, raw.get("settings"), raw.get("resolved"), raw.get("metadata")):
        if isinstance(scope, Mapping):
            for k in keys:
                if k in scope and scope[k] is not None:
                    return scope[k]
    return None


def _num(x: Any) -> Optional[float]:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def format_of(run: Any) -> Format:
    """The format of a run (a ``model.RunData``): its modality, copies and chains, and whether a terminus question exists for it.

    A cyclic peptide has no free ends, so it says so and carries no numbers. A run with several binder chains (an Fab, say) is switched off
    too: such a molecule has two ends per chain and this readout answers for one pair. A homo-oligomer binder still gets an answer, with
    the note that a tag appears on every copy.

    The chain count comes from the designs when the settings do not record it. A BindCraft 2 campaign often does not write
    ``binder_chains``, and reading 1 from a silent settings file is how a two-chain binder came to be given a one-chain answer."""
    from .model import binder_chain_count        # local: this module takes duck-typed objects and imports nothing else of the package

    raw = getattr(getattr(run, "settings", None), "settings_raw", None) or {}
    ann = getattr(run, "annotations", None) or {}
    word = ann.get("modality") or _dig(raw, "modality", "settings_modality") or getattr(getattr(run, "settings", None), "protocol", None)
    if isinstance(word, str) and re.fullmatch(r"\d+stage|4stage", word.strip(), re.I):
        word = None                                                   # BindCraft 1's protocol name is not a modality
    kind, known = "de novo", False
    if isinstance(word, str) and word.strip():
        known = True
        kind = next((k for k, rx in _KIND if rx.search(word)), word.strip())
    copies = int(_num(_dig(raw, "copies")) or 1)
    cyc = _dig(raw, "cyclize_peptide")
    chains_v = _dig(raw, "binder_chains")
    chains = len(chains_v) if isinstance(chains_v, (list, tuple)) else 0
    if chains < 2:
        # The designs are the better witness: the reader counted the chains from
        # Binder_Sequence, which the campaign always writes.
        chains = max((binder_chain_count(d) for d in getattr(run, "designs", []) or []), default=1)
    f = Format(kind=kind, known=known, copies=max(1, copies), chains=max(1, chains))
    f.scaffold = bool(ann.get("scaffold")) or kind in ("VHH", "scFv", "Fab", "ARP")
    lens = [d.length for d in getattr(run, "designs", []) or [] if getattr(d, "length", None)]
    f.peptide = bool(lens) and max(lens) <= PEPTIDE_MAX
    if kind == "cyclic peptide" or (cyc and str(cyc).lower() not in ("false", "0", "none")):
        f.kind, f.applies = "cyclic peptide", False
        f.why = "A head-to-tail cyclic peptide has no free ends, so there is no end to put a tag on. A tag has to go on a side chain instead."
    elif f.copies > 1:
        f.notes.append(f"A homo-oligomer binder of {f.copies} copies: a tag at one end is on every copy at once, and the copies' tags sit near each other.")
    if f.applies and f.chains > 1 and f.copies == 1:
        # The whole binder is loaded now, so this is no longer about what can be
        # read: it is that the readout is one pair of ends and such a binder has
        # two per chain, each with its own exposure and its own distance from the
        # interface. Choosing between them is a judgement this does not make yet.
        f.applies = False
        f.why = (f"The run designed {f.chains} binder chains, so the molecule has {2 * f.chains} ends. This readout answers for one "
                 "pair of ends, and which chain's end to tag is a choice it does not make; look at the chains in the inspector.")
    if kind == "VHH":
        f.convention = "C"
        f.convention_why = ("A VHH's N-terminus lies at the same end as the CDRs and its C-terminus at the opposite end of the domain, "
                            "which is why nanobodies are conventionally tagged at the C-terminus.")
    elif kind in ("scFv", "Fab"):
        f.notes.append(f"A {kind} is built from framework domains whose ends are read from the structure like any other binder's; "
                       "no end is assumed here from convention.")
    if f.peptide and f.applies:
        f.notes.append(f"At {max(lens)} residues or fewer a tag is a large share of the molecule and can change how it folds and binds.")
    return f


def designed_for(run: Any) -> dict[str, Any]:
    """Which terminus properties BindCraft 2 was asked to optimise in this run. The analysis confirms what a run asked for and discovers
    what it did not, and a reader has to know which. ``known`` is False when the run recorded no settings (BindCraft 1, a missing
    file): then nothing is claimed either way."""
    raw = getattr(getattr(run, "settings", None), "settings_raw", None) or {}
    inner = raw.get("settings") if isinstance(raw.get("settings"), Mapping) else raw
    if not isinstance(inner, Mapping) or not inner or not any(str(k).startswith(("weights_", "min_", "max_")) for k in inner):
        return {"known": False, "accessible": False, "together": False, "buried": False, "items": []}
    w = lambda k: (_num(inner.get(k)) or 0.0)  # noqa: E731 - a missing or zero weight is "not asked for"
    acc = any(w(k) > 0 for k in ("weights_termini_angle", "weights_n_terminus_away", "weights_c_terminus_away"))
    tog = w("weights_termini_distance") > 0
    bur = w("weights_exposed_termini") > 0
    items = []
    if acc:
        items.append("both ends angled away from the target (termini_accessible)")
    if tog:
        items.append("the two ends within about 10 A of each other (termini_together)")
    if bur:
        items.append("the ends buried against the fold (protease_stable), the opposite of what a tag wants")
    return {"known": True, "accessible": acc, "together": tog, "buried": bur, "items": items}


# ---------------------------------------------------------------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------------------------------------------------------------
def _unit(v: np.ndarray) -> Optional[np.ndarray]:
    n = float(np.linalg.norm(v))
    return None if not math.isfinite(n) or n < 1e-6 else v / n


def exit_direction(ca: np.ndarray, idx: int, end: str, window: int = EXIT_WINDOW) -> Optional[np.ndarray]:
    """The direction the chain leaves the fold in at residue ``idx`` of one end ('n' or 'c'): from the mean of the next ``window`` CA
    positions inward to the CA itself. The mean sits near the axis of a helix and along a strand, so the vector follows the element and not
    the wobble of one residue. None when fewer than three inward residues exist or a CA is missing."""
    inner = ca[idx + 1: idx + 1 + window] if end == "n" else ca[max(0, idx - window): idx]
    if len(inner) < 3 or np.isnan(inner).any() or np.isnan(ca[idx]).any():
        return None
    return _unit(ca[idx] - inner.mean(axis=0))


def cone_directions(u: np.ndarray, half: float = CONE_DEG, rings: int = CONE_RINGS, azimuths: int = CONE_AZIMUTHS) -> np.ndarray:
    """Unit vectors in a cone of half-angle ``half`` degrees around ``u``: the axis and ``rings`` rings of ``azimuths`` directions."""
    a = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(u, a)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    out = [u]
    for r in range(1, rings + 1):
        th = math.radians(half * r / rings)
        for k in range(azimuths):
            ph = 2.0 * math.pi * k / azimuths
            out.append(math.cos(th) * u + math.sin(th) * (math.cos(ph) * e1 + math.sin(ph) * e2))
    return np.asarray(out)


def _dist(tree: Optional[cKDTree], pts: np.ndarray) -> np.ndarray:
    return np.full(len(pts), np.inf) if tree is None else tree.query(pts)[0]


#: what can stand in the way of a tag, as the codes written into the data ("fold" is the binder's own; the others come from ``env``)
LIMITS = ("target", "fold", "olig", "rest", "glycan")
Env = Callable[[np.ndarray], Mapping[str, Optional[float]]]


def _env_min(env: Optional[Env], pts: np.ndarray) -> tuple[Optional[float], Optional[str]]:
    """(smallest distance, source) from a set of points to what lies outside the construct, or (None, None) when it cannot be checked."""
    if env is None:
        return None, None
    got = {k: v for k, v in (env(pts) or {}).items() if v is not None and math.isfinite(v)}
    if not got:
        return None, None
    k = min(got, key=got.get)
    return float(got[k]), k


class _Flex:
    """Can a flexible chain of a given length get from one end into open solvent without crossing the target or the binder?

    This is the question for an unstructured tag (His6, Twin-Strep), and it is deliberately the weak condition. A flexible chain is not a rod:
    it routes around obstacles, so the end does not need an empty cone or corridor of the tag's full extended length, only some path through
    free space (at least ``PATH_CLEAR`` from every atom) that is no longer than the chain, to a point with the clearance that counts as open
    solvent for the tag. The search is the shortest such path on a grid around the end (``FLEX_SPACING``, ``FLEX_HALF``); what it ignores is
    self-avoidance of the chain and the entropy cost of a tight path, so it can only err toward saying a tag fits. The other subunits and the
    rest of the protein are asked about the best path (``env``); the voxels it crosses that are too close are blocked and the search repeated."""

    def __init__(self, origin: np.ndarray, target: Optional[cKDTree], body: Optional[cKDTree], env: Optional[Env]):
        self.env = env
        n = int(round(FLEX_HALF / FLEX_SPACING))
        self.N = 2 * n + 1
        ax = np.arange(-n, n + 1) * FLEX_SPACING
        self.rel = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1).reshape(-1, 3)
        self.pts = origin + self.rel
        self.near = np.linalg.norm(self.rel, axis=1) <= FLEX_START
        self.c = np.minimum(_dist(target, self.pts), _dist(body, self.pts))
        self.free = (self.c >= PATH_CLEAR) | (self.near & (self.c >= FLEX_START_CLEAR))
        self.free[n * self.N * self.N + n * self.N + n] = True                  # the end itself
        self.blocked = np.zeros(len(self.rel), dtype=bool)
        self.start = n * self.N * self.N + n * self.N + n
        self._dist: Optional[np.ndarray] = None
        self._pred: Optional[np.ndarray] = None
        self._limit = 0.0

    def _graph(self):
        N = self.N
        ok = (self.free & ~self.blocked).reshape(N, N, N)
        idx = np.arange(N ** 3).reshape(N, N, N)
        a_s, b_s, w_s = [], [], []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    if (dx, dy, dz) <= (0, 0, 0):                      # one of each pair of opposite offsets
                        continue
                    sa = tuple(slice(max(0, -d), N - max(0, d)) for d in (dx, dy, dz))
                    sb = tuple(slice(max(0, d), N - max(0, -d)) for d in (dx, dy, dz))
                    m = ok[sa] & ok[sb]
                    a, b = idx[sa][m], idx[sb][m]
                    a_s.append(a)
                    b_s.append(b)
                    w_s.append(np.full(len(a), FLEX_SPACING * math.sqrt(dx * dx + dy * dy + dz * dz)))
        a, b, w = np.concatenate(a_s), np.concatenate(b_s), np.concatenate(w_s)
        return coo_matrix((w, (a, b)), shape=(N ** 3, N ** 3)).tocsr()

    def _solve(self, limit: float) -> None:
        if self._dist is None or limit > self._limit:
            self._dist, self._pred = dijkstra(self._graph(), directed=False, indices=self.start, return_predecessors=True, limit=limit)
            self._limit = limit

    def _path(self, j: int) -> list[int]:
        out = [j]
        while out[-1] != self.start and self._pred[out[-1]] >= 0:
            out.append(int(self._pred[out[-1]]))
        return out

    def escape(self, bulk: float, limit: float) -> Optional[float]:
        """The length in angstrom of the shortest path from the end to a point whose clearance is at least ``bulk``, or None when there is
        none within ``limit`` (or the environment closes every one that is)."""
        excluded = np.zeros(len(self.rel), dtype=bool)
        for _ in range(FLEX_ENV_ROUNDS + 1):
            self._solve(limit)
            cand = (self.c >= bulk) & np.isfinite(self._dist) & ~excluded & (self._dist <= limit)
            if not cand.any():
                return None
            j = int(np.argmin(np.where(cand, self._dist, np.inf)))
            if self.env is None:
                return float(self._dist[j])
            bad_path, bad_end = [], False
            for v in self._path(j):
                if self.near[v] and v != j:
                    continue                                          # the chain leaves the end wherever the model has put it
                d, _src = _env_min(self.env, self.pts[v][None, :])
                if d is not None and d < (bulk if v == j else PATH_CLEAR):
                    if v == j:
                        bad_end = True
                    else:
                        bad_path.append(v)
            if not bad_path and not bad_end:
                return float(self._dist[j])
            if bad_end:
                excluded[j] = True
            for v in bad_path:
                self.blocked |= (np.linalg.norm(self.pts - self.pts[v], axis=1) <= PATH_CLEAR) & ~self.near
            if bad_path:
                self._dist = None
        return None


class _Probe:
    """The free-space profile along a cone of directions from one origin, against the target and the binder's own fold, with the
    environment asked only about the few candidates that pass those two."""

    def __init__(self, origin: np.ndarray, u: np.ndarray, target: Optional[cKDTree], body: Optional[cKDTree], env: Optional[Env]):
        self.o, self.env, self.tt, self.bt = origin, env, target, body
        self.dirs = cone_directions(u)
        max_reach = TAGS[-1][2] + STANDOFF + REACH * LINKER_STEPS[-1] + 12 * REACH          # room for a tail credit as well
        self.s = np.arange(0.0, max_reach + PATH_STEP, PATH_STEP)
        pts = origin[None, None, :] + self.dirs[:, None, :] * self.s[None, :, None]            # (K, S, 3)
        flat = pts.reshape(-1, 3)
        self.dt = _dist(target, flat).reshape(len(self.dirs), len(self.s))
        self.db = _dist(body, flat).reshape(len(self.dirs), len(self.s))
        comb = np.minimum(self.dt, self.db)
        # the first PATH_STEP*2 of the path is the end's own surroundings: the path is judged from there outward
        start = int(round(2.0 / PATH_STEP))
        masked = comb.copy()
        masked[:, :start] = np.inf
        self.path_min = np.minimum.accumulate(masked, axis=1)                                  # smallest distance on the path up to s

    def _at(self, k: int, reach: float) -> np.ndarray:
        return self.o + self.dirs[k] * reach

    def _idx(self, reach: float) -> int:
        return min(len(self.s) - 1, int(round(reach / PATH_STEP)))

    def fits(self, r: float, reach: float) -> list[int]:
        """Directions whose sphere of radius ``r`` at ``reach`` along them clears the target and the fold, by a path that clears too."""
        i = self._idx(reach)
        need = r + ATOM_R - SLACK
        path_to = self._idx(max(0.0, reach - r))
        ok = np.minimum(self.dt[:, i], self.db[:, i]) >= need
        ok &= self.path_min[:, path_to] >= PATH_CLEAR
        return [int(k) for k in np.nonzero(ok)[0]]

    def env_ok(self, k: int, r: float, reach: float) -> bool:
        """Does the candidate also clear what lies outside the construct (the sphere and the path to it)? True when that cannot be checked."""
        if self.env is None:
            return True
        dc, _ = _env_min(self.env, self._at(k, reach)[None, :])
        if dc is not None and dc < r + ATOM_R - SLACK:
            return False
        if reach - r > 2.0:
            path = np.array([self._at(k, 2.0 + j * 4.0) for j in range(max(1, int((reach - r - 2.0) // 4.0) + 1))])
            dp, _ = _env_min(self.env, path)
            if dp is not None and dp < PATH_CLEAR:
                return False
        return True

    def clear(self, r: float, reach: float) -> bool:
        """Is there a direction with room for a sphere of radius ``r`` at ``reach``, that the environment leaves clear as well?"""
        cand = self.fits(r, reach)
        if not cand:
            return False
        i = self._idx(reach)
        cand.sort(key=lambda k: -min(self.dt[k, i], self.db[k, i]))
        return any(self.env_ok(k, r, reach) for k in cand[:ENV_TRIES])

    def room(self, credit: int) -> float:
        """The radius of the largest sphere that fits at its standoff with no linker (0 when none does). The fold and the target are
        asked first, for every size; the environment only for the largest few that pass, from the top."""
        sizes = [float(r) for r in np.arange(2.0, ROOM_CAP + 1.0, 1.0)]
        reach = lambda r: r + STANDOFF + REACH * credit  # noqa: E731
        own = [r for r in sizes if self.fits(r, reach(r))]
        for r in sorted(own, reverse=True):
            if self.clear(r, reach(r)):
                return r
        return 0.0

    def limit(self, r: float, reach: float) -> str:
        """What limits a sphere of radius ``r`` at ``reach``: of the directions, the one with the most room, and the source nearest it."""
        i = self._idx(reach)
        k = int(np.argmax(np.minimum(self.dt[:, i], self.db[:, i])))
        got = {"target": float(self.dt[k, i]), "fold": float(self.db[k, i])}
        if self.env is not None:
            d, src = _env_min(self.env, self._at(k, reach)[None, :])
            if d is not None and src in LIMITS:
                got[src] = d
        return min(got, key=got.get)


# ---------------------------------------------------------------------------------------------------------------------------------
# one design
# ---------------------------------------------------------------------------------------------------------------------------------
def _terminal_stats(a: Any, end: str) -> dict[str, Any]:
    """Secondary structure, burial and confidence of the residues at one end (a ``DesignAnalysis``'s per-residue lists)."""
    n = int(a.n)
    ss = str(a.ss or "")
    order = range(n) if end == "n" else range(n - 1, -1, -1)
    tail = 0
    for i in order:
        if i < len(ss) and ss[i] == "L":
            tail += 1
        else:
            break
    anchor = None if tail >= n or not ss else (tail if end == "n" else n - 1 - tail)       # 0-based index of the last ordered residue
    elem_ss, elem_n = None, 0
    if anchor is not None and anchor < len(ss):
        elem_ss = ss[anchor]
        j = anchor
        step = 1 if end == "n" else -1
        while 0 <= j < len(ss) and ss[j] == elem_ss:
            elem_n += 1
            j += step
    idx = list(order)[:max(RSA_N, PLDDT_N)]
    rsa = [a.rsa_res[i] for i in idx[:RSA_N] if i < len(a.rsa_res) and a.rsa_res[i] is not None]
    pl = [a.plddt_res[i] for i in idx[:PLDDT_N] if i < len(a.plddt_res) and a.plddt_res[i] is not None]
    tail_pl = [a.plddt_res[i] for i in list(order)[:tail] if i < len(a.plddt_res) and a.plddt_res[i] is not None]
    bc = set(int(b) for b in (a.bc or []))
    res1 = [i + 1 for i in idx[:3]]                                                           # 1-based binder numbers of the three terminal residues
    return {"tail": tail, "anchor": anchor, "elem_ss": elem_ss, "elem_n": elem_n,
            "rsa": float(np.mean(rsa)) if rsa else None, "plddt": float(np.mean(pl)) if pl else None,
            "tail_plddt": float(np.mean(tail_pl)) if tail_pl else None,
            "contact": (res1[0] in bc) if res1 else False, "iface_n": sum(1 for r in res1 if r in bc)}


def _round(x: Optional[float], nd: int = 1) -> Optional[float]:
    return None if x is None or not math.isfinite(x) else round(float(x), nd)


def _side(d: Optional[float], away: Optional[float], contact: bool) -> str:
    """'paratope', 'lateral' or 'far': where an end sits relative to the binding face. A tag there holds or displays the binder by the face
    it binds with, whatever room the end has."""
    if contact or (d is not None and d < NEAR) or (away is not None and away < PARA_AWAY):
        return "paratope"
    if away is not None and away >= FAR_AWAY:
        return "far"
    return "lateral"


def _end(a: Any, bnd: Any, ca: np.ndarray, end: str, target_xyz: np.ndarray, t_tree: cKDTree, env: Optional[Env],
         centroid: np.ndarray, centre: Optional[np.ndarray]) -> dict[str, Any]:
    n = len(ca)
    i0 = 0 if end == "n" else n - 1
    tip = ca[i0]
    out: dict[str, Any] = {"d": None, "d_env": None, "env_what": None, "at": None, "dir": None, "aim": None, "away": None, "side": None, "contact": False, "iface_n": 0, "rsa": None,
                           "plddt": None, "ss": None, "elem": 0, "tail": 0, "room": None, "limit": None, "g": None, "conf": "ok", "via_tail": False}
    if np.isnan(tip).any():
        out["conf"] = "off"
        out["why"] = "the terminal residue has no CA atom"
        return out
    st = _terminal_stats(a, end)
    out.update(rsa=_round(st["rsa"], 2), plddt=_round(st["plddt"], 1), ss=st["elem_ss"], elem=st["elem_n"], tail=st["tail"],
               contact=bool(st["contact"]), iface_n=int(st["iface_n"]))
    dd, _ = t_tree.query(tip)
    out["d"] = _round(float(dd), 1)
    if env is not None:
        d_env, src = _env_min(env, tip[None, :])
        out["d_env"], out["env_what"] = _round(d_env, 1), src
    # A long disordered tail is judged from its last ordered residue, with the tail counted as linker: where the tip of a floppy tail
    # lies is a guess, and the model is no better a guide to where it would lie in a tagged protein. Everything that decides the
    # verdict (the side of the binder, whether the end touches the target, the direction) is then read at that residue.
    disordered = st["tail"] >= TAIL_MIN and st["anchor"] is not None and (st["tail_plddt"] is None or st["tail_plddt"] < DISORDER_PLDDT)
    low = disordered or (st["plddt"] is not None and st["plddt"] < DISORDER_PLDDT)
    oi = st["anchor"] if disordered else i0
    origin = ca[oi]
    credit = min(st["tail"], TAIL_CREDIT_MAX) if disordered else 0
    d_origin = float(t_tree.query(origin)[0])
    touch = bool(st["contact"]) if not disordered else (oi + 1) in set(int(b) for b in (a.bc or []))
    # where the end sits relative to the paratope: +1 on the far side of the binder, -1 on the interface side
    if centre is not None:
        v1, v2 = _unit(origin - centroid), _unit(centroid - centre)
        if v1 is not None and v2 is not None:
            out["away"] = _round(float(np.dot(v1, v2)), 2)
    out["side"] = _side(d_origin, out["away"], touch)
    u = exit_direction(ca, oi, end)
    if u is None:
        out["conf"] = "off"
        out["why"] = "the end is too short to give a direction"
        return out
    out["via_tail"] = bool(disordered)
    out["conf"] = "low" if low else "ok"
    out["at"] = int(oi) + 1                                                 # the residue the probe starts from (1-based): the end itself, or the last ordered one
    out["dir"] = [round(float(x), 3) for x in u]                            # in the design's frame; the exporter turns it into the placed frame
    # the direction relative to the nearest target atom, from where the probe starts (0 = pointing at it)
    _, j = t_tree.query(origin)
    vu = _unit(target_xyz[j] - origin)
    if vu is not None:
        out["aim"] = _round(math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(u, vu)))))), 0)
    # the binder's fold, without the residues a tag grows out of: the terminal one, or the whole disordered tail and the residue it hangs from
    gone = {i0} if not disordered else set(range(0, oi + 1)) if end == "n" else set(range(oi, n))
    keep = ~np.isin(bnd.res_index, list(gone))
    body_tree = cKDTree(bnd.xyz[keep]) if keep.any() else None
    pr = _Probe(origin, u, t_tree, body_tree, env)
    flex = _Flex(origin, t_tree, body_tree, env)
    g: dict[str, Any] = {}
    for key, _what, r in TAGS:
        if key in FLEXIBLE:
            # an unstructured chain: can it route into open solvent? Its own length, less the motif that has to lie in the open, is its reach.
            res, motif = FLEXIBLE[key]
            budget = (res - motif + credit) * REACH
            esc = flex.escape(r, budget + LINKER_STEPS[-1] * REACH)
            found = None if esc is None else next((s for s in LINKER_STEPS if esc <= budget + s * REACH), None)
            if found is None:
                g[key] = {"s": "blocked", "l": None, "why": "room", "esc": None}
            else:
                g[key] = {"s": "free" if found == 0 else "linker", "l": int(found), "esc": _round(esc, 0)}
            continue
        if touch:
            g[key] = {"s": "blocked", "l": None, "why": "contact"}            # the end itself touches the target: a folded tag grows out of the interface
            continue
        found = next((s for s in LINKER_STEPS if pr.clear(r, r + STANDOFF + REACH * (credit + s))), None)
        g[key] = {"s": "blocked", "l": None, "why": "room"} if found is None else {"s": "free" if found == 0 else "linker", "l": int(found)}
    out["g"] = g
    out["room"] = 0.0 if touch else pr.room(credit)
    top = next((t for t in TAGS if g[t[0]]["s"] != "free"), None)
    if top is not None:
        out["limit"] = "target" if touch else pr.limit(top[2], top[2] + STANDOFF + REACH * credit)
    return out


def _better(gn: Mapping[str, Any], gc: Mapping[str, Any]) -> str:
    """The better of two usable ends by fit alone: 'N', 'C' or 'both'. A linker three or more residues shorter decides between two ends that
    both need one; below that the difference is within what a model can tell."""
    rn, rc = STATUS_RANK[gn["s"]], STATUS_RANK[gc["s"]]
    if rn != rc:
        return "N" if rn < rc else "C"
    if rn == 0:
        return "both"
    d = gn["l"] - gc["l"]
    return "both" if abs(d) < 3 else ("N" if d < 0 else "C")


def _pick(n_e: Mapping[str, Any], c_e: Mapping[str, Any], key: str) -> str:
    """Which end to tag for one tag size: 'N', 'C', 'both' (either end is as good) or 'none'.

    Room comes first: an end with none for this tag is out. Then the side: of two ends with room, one on the paratope side loses to one that
    is not, because a tag there would hold or display the binder by its binding face. Only when both are on that side does fit decide."""
    gn, gc = n_e["g"][key], c_e["g"][key]
    un, uc = gn["s"] != "blocked", gc["s"] != "blocked"
    if not un and not uc:
        return "none"
    if un != uc:
        return "N" if un else "C"
    pn, pc = n_e["side"] == "paratope", c_e["side"] == "paratope"
    if pn != pc:
        return "C" if pn else "N"
    return _better(gn, gc)


def assess(a: Any, binder: Any, target: Any, fmt: Format, env: Optional[Env] = None) -> dict[str, Any]:
    """The terminus readout of one design.

    ``a`` is its ``interface.DesignAnalysis`` (secondary structure, burial, pLDDT, interface residues and centre); ``binder`` and ``target`` are
    the design's own chains as ``structure.Atoms`` in the design's frame; ``fmt`` is the run's ``Format``. ``env`` answers, for points of
    the design's frame, how close they come to what lies outside the construct (None: that check is off). Returns a plain dict, ready for
    the blobs (``docs/termini.md``)."""
    base: dict[str, Any] = {"v": SCHEMA, "applies": bool(fmt.applies), "format": fmt.kind, "copies": fmt.copies}
    if not fmt.applies:
        base["why"] = fmt.why
        return base
    bnd = binder.heavy() if binder.is_hydrogen.any() else binder
    tgt = target.heavy() if target.is_hydrogen.any() else target
    ca = bnd.ca_xyz()
    if len(ca) < 6 or len(tgt) == 0:
        base.update(applies=False, why="the structure has too few residues or no target to judge an end against")
        return base
    t_xyz = tgt.xyz
    t_tree = cKDTree(t_xyz)
    centroid = np.nanmean(ca, axis=0)
    centre = np.asarray(a.center, dtype=float) if getattr(a, "center", None) else None
    n_e = _end(a, bnd, ca, "n", t_xyz, t_tree, env, centroid, centre)
    c_e = _end(a, bnd, ca, "c", t_xyz, t_tree, env, centroid, centre)
    base.update(n=n_e, c=c_e)
    nc = None if np.isnan(ca[0]).any() or np.isnan(ca[-1]).any() else float(np.linalg.norm(ca[0] - ca[-1]))
    base["nc_d"] = _round(nc, 1)
    base["env"] = env is not None
    if n_e.get("g") and c_e.get("g"):
        base["pick"] = {k: _pick(n_e, c_e, k) for k in TAG_KEYS}
        base["immob"] = _immobilise(n_e, c_e)
    else:
        base["pick"] = None
        base["immob"] = None
    return base


def _immobilise(n_e: Mapping[str, Any], c_e: Mapping[str, Any]) -> Optional[str]:
    """The end to hold the binder by (SPR, BLI, a coated plate) so that the paratope faces away from the surface: of the ends that can
    take a His6, those on the far side of the binder from the interface; the better one, or 'either' when they are alike. None when
    no such end exists, which is the case a person has to know about."""
    good = [(e, d["away"]) for e, d in (("N", n_e), ("C", c_e)) if d["g"]["his6"]["s"] != "blocked" and d["away"] is not None and d["away"] > 0]
    if not good:
        return None
    if len(good) == 1:
        return good[0][0]
    return "either" if abs(good[0][1] - good[1][1]) < 0.3 else max(good, key=lambda x: x[1])[0]


def family_pick(picks: Sequence[str]) -> Optional[str]:
    """The end to tag for a family: what most of its sequences say, and C, N, both, none in that order when they are tied. Sequences on one
    backbone refold with their ends in somewhat different places, so one sequence is not the family; a family whose sequences do not all agree
    is marked as split by the pages. The rule does not look at any quality score, so the pages and the exporter judge a family alike
    (``Common.term.famPick``)."""
    if not picks:
        return None
    cnt = {k: sum(1 for p in picks if p == k) for k in ("C", "N", "both", "none")}
    top = max(cnt.values())
    return next(k for k, v in cnt.items() if v == top)


def tag_table() -> list[dict[str, Any]]:
    """The tag sizes as written into the run block, so the pages and the data agree on them. ``kind`` says how a tag is tested; ``r`` is the
    radius of a globular tag's sphere or the clearance that counts as open solvent for a flexible one; ``res``, ``motif`` and ``budget`` (the
    reach in angstrom a flexible chain has beyond its motif) are given for the flexible tags only."""
    out = []
    for k, w, r in TAGS:
        row: dict[str, Any] = {"key": k, "what": w, "r": r, "kind": "flexible" if k in FLEXIBLE else "globular"}
        if k in FLEXIBLE:
            res, motif = FLEXIBLE[k]
            row.update(res=res, motif=motif, budget=round((res - motif) * REACH, 1))
        out.append(row)
    return out


def constants() -> dict[str, Any]:
    """Every threshold of the analysis, for the manifest and the documentation."""
    return {"cone_deg": CONE_DEG, "standoff": STANDOFF, "reach": REACH, "linker_steps": list(LINKER_STEPS), "atom_r": ATOM_R, "slack": SLACK,
            "path_clear": PATH_CLEAR, "flex_spacing": FLEX_SPACING, "flex_half": FLEX_HALF, "flex_start": FLEX_START, "room_cap": ROOM_CAP, "para_away": PARA_AWAY, "far_away": FAR_AWAY, "exit_window": EXIT_WINDOW, "tail_min": TAIL_MIN, "disorder_plddt": DISORDER_PLDDT, "near": NEAR}
