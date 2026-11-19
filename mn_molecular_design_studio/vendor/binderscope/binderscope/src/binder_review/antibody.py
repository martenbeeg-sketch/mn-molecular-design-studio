"""Antibody-style profiles: CDR numbering, framework and CDR spans, hallmark residues, cysteine pairing, CDR
liabilities and the CDR share of an interface, for VH, VHH, VL and scFv sequences.

Everything here is a screen, not a verdict. Numbering is IMGT and is stored beside the design; coordinates are
never renumbered. Positions are 1-based indices into the sequence that was passed in.

Two backends give the same result shape:

* ``abnumber`` (with ``anarci`` and ``hmmscan``) when importable: real IMGT numbering, the ground truth.
* a heuristic fallback (``fallback``) that needs nothing but the standard library. It finds the four anchors
  every well-formed V domain has (Cys23, Trp41, Cys104 and the J-region motif WGxG or FGxG) and places the
  IMGT CDR limits at fixed distances from them. It was calibrated against abnumber on 11 germline-like VH,
  VHH, kappa and lambda sequences and reproduces all CDR limits there. It does not number insertions, so it
  returns no IMGT labels beyond the hallmark positions, and it gives up (``ok`` False) on sequences without
  the anchors. Use it for well-formed domains only.

Entry points:

* :func:`profile_sequence` - one sequence (VHH, VH, VL, scFv or Fv with a linker): the full profile dict.
* :func:`number_chain` - one V domain: regions, spans, IMGT labels.
* :func:`cdr_interface_share` - CDR share of a contact list.
* :func:`cdr_liabilities` - liabilities of :mod:`liabilities` tagged with their region, plus CDR rules.
* :func:`cysteine_report` - pair-aware Cys check.

Exporter hook: ``profile_sequence(seq, contact_residues=..., sg_xyz=...)`` returns a dict whose ``regions``,
``region_counts`` and ``cdrs`` entries are exactly the inspector blob fields of the same names, and ``nH``
is set for scFv (residues in the VH domain). See docs/antibody.md.
"""
from __future__ import annotations

import itertools
import re
import warnings
from typing import Iterable, Mapping, Optional, Sequence

from . import liabilities as _liab

REGION_ORDER = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
CDR_NAMES = ("CDR1", "CDR2", "CDR3")
FLANK = "tag"                     # residues outside the V domain (leader, tag, linker)
CANONICAL_CYS = (23, 104)         # IMGT numbers of the conserved disulfide
SG_SG_MAX = 2.5                   # A, a disulfide S-S bond is 2.03 A; allow slack for models

# Hallmark positions of camelid VHH (IMGT 42, 49, 50, 52 = Kabat 37, 44, 45, 47) and the usual human VH.
HALLMARK_POSITIONS = (42, 49, 50, 52)
_VH_HALLMARK = "VGLW"
_VHH_ALLOWED = ({"F", "Y"}, {"E", "Q"}, {"R", "C"}, {"G", "L", "F"})

try:                                                  # pragma: no cover - depends on the environment
    from abnumber import Chain as _AbChain
    from abnumber.exceptions import ChainParseError as _ChainParseError
    HAVE_ABNUMBER = True
except Exception:                                     # ImportError, or hmmscan/anarci missing at import
    _AbChain = None
    _ChainParseError = Exception
    HAVE_ABNUMBER = False


def _clean(seq: str) -> str:
    return "".join(c for c in str(seq).upper() if c.isalpha())


# --------------------------------------------------------------------------------------------- fallback

_FR4_H = re.compile(r"W[GA].G[A-Z]{2}[VLI]T[VLI][SA]")
_FR4_L = re.compile(r"F[GA].G[A-Z][KRT][VLT][A-Z]{2}[KRL]")
_J_TAIL = set("ACDEFIKLMNPQRSTVWY")   # allowed inside FR4 after the J motif; G and H stop it (linker, His tag)
_FR4_H_SHORT = re.compile(r"W[GA].G[TQSAK]")
_FR4_L_SHORT = re.compile(r"F[GA].G[TQSAG]")
_FR2_ANY = re.compile(r"W[YFVIAL][RQ][QRK]")


def fallback(seq: str) -> Optional[dict]:
    """Heuristic IMGT region limits of one V domain; None when the anchors are not found.

    Returns ``{"chain_type", "start", "end", "spans", "imgt"}`` with 0-based inclusive indices. Anchors:
    Cys23 (``[LVIM]S?C`` followed by four residues and a CDR1 that ends two residues before Trp41),
    Trp41 (FR2 ``W[YFVIA][RQ][QRK]``), Cys104 (``Y[YFHW]C`` before the J motif) and the J motif.
    CDR1 = Cys23+4 to Trp41-3; CDR2 starts 15 residues after Trp41 and ends 38 (heavy) or 36 (light)
    residues before Cys104; CDR3 runs from Cys104+1 to the residue before the J motif.
    """
    s = _clean(seq)
    if len(s) < 85:
        return None
    best = None
    for wm in _FR2_ANY.finditer(s):
        w = wm.start()
        if not 28 <= w <= 55:
            continue
        c23 = s.rfind("C", 0, w - 9)
        if c23 < 10 or not 10 <= w - c23 <= 22:
            continue
        # Cys104 sits 53 (light) to 60 (heavy) residues after Trp41 and follows an aromatic residue; later Cys
        # in CDR3 (extra disulfides) must not be taken for it.
        c104 = next((c for c in range(w + 50, min(len(s) - 1, w + 64)) if s[c] == "C" and s[c - 1] in "YFHWV"), None)
        if c104 is None:
            continue
        # Full J motif first; the short WGxG / FGxG form covers domains whose FR4 is truncated (crystal Fv).
        for fr4s in ((_FR4_H, _FR4_L), (_FR4_H_SHORT, _FR4_L_SHORT)):
            for kind, fr4 in zip("HL", fr4s):
                m = fr4.search(s, c104 + 3)
                if m and m.start() - c104 <= 40:
                    best = (kind, m.start(), c104, w, c23)
                    break
            if best:
                break
        if best:
            break
    if best is None:
        return None
    kind, j, c104, w, c23 = best
    chain_type = "H"
    if kind == "L":
        chain_type = "K" if s[j + 9:j + 10] in ("K", "R") or s[j + 10:j + 11] in ("K", "R") else "L"
    drop = 38 if kind == "H" else 36
    # IMGT 1 lies 22 residues before Cys23 (heavy) or 22 to 23 (light); anything before it is a leader.
    start = max(0, c23 - (21 if kind == "H" else 22))
    cdr2_end = c104 - drop
    spans = {
        "FR1": (start, c23 + 3),
        "CDR1": (c23 + 4, w - 3),
        "FR2": (w - 2, w + 14),
        "CDR2": (w + 15, cdr2_end),
        "FR3": (cdr2_end + 1, c104),
        "CDR3": (c104 + 1, j - 1),
    }
    fr4_end = j + 3
    while fr4_end + 1 < min(len(s), j + (11 if kind == "H" else 10)) and s[fr4_end + 1] in _J_TAIL:
        fr4_end += 1
    spans["FR4"] = (j, fr4_end)
    for name, (a, b) in spans.items():
        if name != "FR1" and b < a and name in CDR_NAMES:
            return None
    if not all(spans[c][1] >= spans[c][0] for c in ("CDR1", "CDR2", "CDR3")):
        return None
    imgt = {23: c23, 41: w, 42: w + 1, 49: w + 8, 50: w + 9, 52: w + 11, 104: c104}
    return {"chain_type": chain_type, "start": start, "end": fr4_end, "spans": spans, "imgt": imgt}


# --------------------------------------------------------------------------------------------- numbering

def _abnumber_chain(seq: str) -> Optional[dict]:
    if not HAVE_ABNUMBER:
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ch = _AbChain(seq, scheme="imgt")
    except Exception:
        return None
    start = seq.find(ch.seq)
    if start < 0:
        return None
    regions, labels, imgt = [], [], {}
    for i, pos in enumerate(ch.positions):
        regions.append(pos.get_region())
        lab = f"{pos.number}{pos.letter}".strip()
        labels.append(lab)
        if not pos.letter:
            imgt[int(pos.number)] = start + i
    # ANARCI stops at IMGT 126 or so; the last residues of the J segment still belong to FR4.
    fr4_len = 11 if ch.chain_type == "H" else 10
    end = start + len(ch.seq) - 1
    while regions[-1] == "FR4" and regions.count("FR4") < fr4_len and end + 1 < len(seq) and seq[end + 1] in _J_TAIL:
        regions.append("FR4")
        labels.append(None)
        end += 1
    return {"chain_type": ch.chain_type, "start": start, "end": end,
            "regions": regions, "labels": labels, "imgt": imgt}


def number_chain(seq: str, backend: str = "auto") -> dict:
    """Number one V domain (VHH, VH, VL). ``backend`` is "auto", "abnumber" or "fallback".

    Result: ``ok``, ``backend``, ``chain_type`` (H, K or L), ``regions`` (per residue of ``seq``: FR1 ... FR4,
    or "tag" outside the domain), ``labels`` (IMGT label per residue, None outside the domain or with the
    fallback), ``spans`` (``{"CDR1": (first, last), ...}``, 1-based inclusive, also FR1 to FR4), ``cdrs``
    (CDR sequences) and ``imgt`` (``{IMGT number: 1-based sequence position}`` of the unlettered positions).
    On failure ``ok`` is False and ``error`` says why.
    """
    s = _clean(seq)
    n = len(s)
    res: Optional[dict] = None
    used = None
    if backend in ("auto", "abnumber"):
        res = _abnumber_chain(s)
        used = "abnumber" if res else None
    if res is None and backend in ("auto", "fallback"):
        fb = fallback(s)
        if fb:
            regions = [""] * n
            for name, (a, b) in fb["spans"].items():
                for i in range(a, b + 1):
                    regions[i] = name
            res = {"chain_type": fb["chain_type"], "start": fb["start"], "end": fb["end"], "regions": regions[fb["start"]:fb["end"] + 1],
                   "labels": None, "imgt": dict(fb["imgt"])}
            used = "fallback"
    if res is None:
        return {"ok": False, "backend": None, "error": "no V domain found", "n": n}
    regions = [FLANK] * n
    labels: list = [None] * n
    for i, r in enumerate(res["regions"]):
        regions[res["start"] + i] = r
        if res["labels"]:
            labels[res["start"] + i] = res["labels"][i]
    spans: dict = {}
    for r in REGION_ORDER:
        idx = [i + 1 for i, x in enumerate(regions) if x == r]
        if idx:
            spans[r] = (idx[0], idx[-1])
    cdrs = {c: (s[spans[c][0] - 1:spans[c][1]] if c in spans else "") for c in CDR_NAMES}
    imgt = {k: v + 1 for k, v in res["imgt"].items()}
    return {"ok": True, "backend": used, "chain_type": res["chain_type"], "n": n, "regions": regions,
            "labels": labels, "spans": spans, "cdrs": cdrs, "imgt": imgt,
            "domain": (res["start"] + 1, res["end"] + 1)}


# --------------------------------------------------------------------------------------------- hallmarks

def hallmarks(seq: str, numbering: Mapping) -> dict:
    """Hallmark residues of a heavy-chain domain at IMGT 42, 49, 50 and 52 and a VHH-versus-VH call.

    Camelid VHH replace the VH/VL interface residues V42, G49, L50, W52 by hydrophilic ones (typically
    F/Y42, E/Q49, R/C50, G/L/F52). ``call`` is "VHH-like" when at least three of the four fit that set,
    "VH-like" when at least three match V, G, L, W, else "mixed". Light chains return ``call`` None.
    """
    if not numbering.get("ok") or numbering.get("chain_type") != "H":
        return {"residues": {}, "string": "", "call": None}
    s = _clean(seq)
    imgt = numbering["imgt"]
    res = {}
    for p in HALLMARK_POSITIONS:
        i = imgt.get(p)
        if i:
            res[p] = s[i - 1]
    string = "".join(res.get(p, "-") for p in HALLMARK_POSITIONS)
    vhh = sum(1 for p, ok in zip(HALLMARK_POSITIONS, _VHH_ALLOWED) if res.get(p) in ok)
    vh = sum(1 for p, a in zip(HALLMARK_POSITIONS, _VH_HALLMARK) if res.get(p) == a)
    call = "VHH-like" if vhh >= 3 else "VH-like" if vh >= 3 else "mixed"
    return {"residues": res, "string": string, "call": call}


# --------------------------------------------------------------------------------------------- cysteines

def cysteine_report(seq: str, numbering: Mapping, sg_xyz: Optional[Mapping[int, Sequence[float]]] = None) -> dict:
    """Pair-aware cysteine check of one V domain.

    The conserved pair is Cys23 and Cys104. Every other Cys is "extra". Extra Cys that can pair up (an even
    number, or pairs confirmed by SG-SG distance when ``sg_xyz`` {position: (x, y, z)} is given) are not free
    thiols: VHH commonly carry a second disulfide between CDR3 and CDR1, CDR2 or FR2 (position 50). Result keys:

    ``canonical``  {23: pos or None, 104: pos or None, "present": bool}
    ``extra``      [{pos, region}]
    ``pairs``      [{a, b, regions, source ("distance" or "parity"), canonical}] including the canonical pair
    ``free``       positions of Cys that stay unpaired (a free thiol)
    ``status``     "ok" (only the canonical pair), "extra disulfide" (extra Cys all paired), "free cysteine",
                   "canonical pair broken" (Cys23 or Cys104 missing)
    ``note``       one line for the page
    """
    s = _clean(seq)
    regions = numbering.get("regions") or [""] * len(s)
    imgt = numbering.get("imgt") or {}
    c23, c104 = imgt.get(23), imgt.get(104)
    cys = [i + 1 for i, c in enumerate(s) if c == "C"]
    canon_ok = bool(c23 and c104 and s[c23 - 1] == "C" and s[c104 - 1] == "C")
    canon = {c23, c104} if canon_ok else set()
    extra = [p for p in cys if p not in canon]
    pairs: list[dict] = []
    free: list[int] = []
    reg = lambda p: regions[p - 1] if p - 1 < len(regions) else ""
    if canon_ok:
        pairs.append({"a": c23, "b": c104, "regions": (reg(c23), reg(c104)), "source": "parity", "canonical": True})
    if sg_xyz:
        have = {p: sg_xyz[p] for p in extra if p in sg_xyz}
        for a, b in itertools.combinations(sorted(have), 2):
            d = sum((x - y) ** 2 for x, y in zip(have[a], have[b])) ** 0.5
            if d <= SG_SG_MAX:
                pairs.append({"a": a, "b": b, "regions": (reg(a), reg(b)), "source": "distance", "canonical": False,
                              "distance": round(d, 2)})
        paired = {p for pr in pairs for p in (pr["a"], pr["b"])}
        if canon_ok and c23 in sg_xyz and c104 in sg_xyz:
            d = sum((x - y) ** 2 for x, y in zip(sg_xyz[c23], sg_xyz[c104])) ** 0.5
            pairs[0]["source"], pairs[0]["distance"] = "distance", round(d, 2)
        rest = [p for p in extra if p not in paired]
        # Without coordinates for every Cys, fall back to parity for the unresolved ones.
        unresolved = [p for p in rest if p not in sg_xyz]
        free = [p for p in rest if p in sg_xyz]
        extra_pairs = unresolved
    else:
        extra_pairs = extra
    if len(extra_pairs) % 2 == 0:
        # Pair the nearest in sequence, CDR3 partner first: the usual VHH arrangement.
        left = sorted(extra_pairs)
        cdr3 = [p for p in left if reg(p) == "CDR3"]
        other = [p for p in left if reg(p) != "CDR3"]
        while cdr3 and other:
            a = other.pop(0)
            b = cdr3.pop(0)
            pairs.append({"a": min(a, b), "b": max(a, b), "regions": (reg(min(a, b)), reg(max(a, b))), "source": "parity", "canonical": False})
        rest2 = sorted(cdr3 + other)
        for a, b in zip(rest2[::2], rest2[1::2]):
            pairs.append({"a": a, "b": b, "regions": (reg(a), reg(b)), "source": "parity", "canonical": False})
    else:
        free.extend(extra_pairs)
    free = sorted(set(free))
    if not canon_ok:
        status = "canonical pair broken"
        note = "Cys23 or Cys104 (the conserved disulfide) is missing; the domain is unlikely to fold as an immunoglobulin V domain."
    elif free:
        status = "free cysteine"
        note = (f"{len(free)} cysteine(s) without a partner at position(s) {', '.join(map(str, free))}: a free thiol "
                "(dimerisation, scrambling, conjugation side products).")
    elif extra:
        status = "extra disulfide"
        pr = [p for p in pairs if not p["canonical"]]
        note = ("Extra cysteines are paired ("
                + "; ".join(f"{p['a']} ({p['regions'][0]}) with {p['b']} ({p['regions'][1]}), by {p['source']}" for p in pr)
                + "). Extra CDR disulfides are common in camelid VHH with a long CDR3; both Cys are paired, so there is no free thiol."
                + (" Pairing was assumed from parity, confirm with the structure." if any(p["source"] == "parity" for p in pr) else ""))
    else:
        status = "ok"
        note = "Only the conserved Cys23 to Cys104 disulfide."
    return {"canonical": {23: c23, 104: c104, "present": canon_ok}, "extra": [{"pos": p, "region": reg(p)} for p in extra],
            "pairs": pairs, "free": free, "status": status, "note": note}


# --------------------------------------------------------------------------------------------- interface

def cdr_interface_share(regions: Sequence[str], contact_residues: Iterable[int]) -> dict:
    """CDR share of an interface. ``contact_residues`` are 1-based binder residues in contact with the target.

    Returns ``{"n", "counts" (per region), "cdr", "framework", "other", "share", "by_cdr" (share per CDR)}``;
    ``share`` is the fraction of contacting V-domain residues that lie in CDR1 to CDR3 (the page meter
    "CDR interface", warn below 0.6); None when no residue of a V domain is in contact. Residues outside the
    V domain (tags, linkers) are counted in ``other`` and left out of the share.
    """
    counts: dict = {}
    for p in sorted({int(x) for x in contact_residues}):
        if 1 <= p <= len(regions):
            counts[regions[p - 1]] = counts.get(regions[p - 1], 0) + 1
    cdr = sum(v for k, v in counts.items() if k.startswith("CDR"))
    fr = sum(v for k, v in counts.items() if k.startswith("FR"))
    other = sum(v for k, v in counts.items() if not (k.startswith("CDR") or k.startswith("FR")))
    tot = cdr + fr
    return {"n": tot + other, "counts": counts, "cdr": cdr, "framework": fr, "other": other,
            "share": (cdr / tot) if tot else None,
            "by_cdr": {c: (counts.get(c, 0) / tot if tot else None) for c in CDR_NAMES}}


# --------------------------------------------------------------------------------------------- liabilities

def cdr_liabilities(seq: str, regions: Sequence[str], contact_residues: Optional[Iterable[int]] = None,
                    rsa: Optional[Sequence[float]] = None, ssbonds: Optional[Iterable[Sequence[int]]] = None) -> list[dict]:
    """Sequence liabilities (see :func:`liabilities.find_liabilities`) with ``region`` and ``in_cdr`` added.

    Antibody rules on top: an N-glycosylation sequon, an NG or a DG motif inside a CDR is raised to "bad"
    (it cannot be fixed outside the paratope and may change binding); Met and Trp in a CDR stay as the base
    module reports them. Framework liabilities keep their severity. Pass ``ssbonds`` (pairs from
    :func:`cysteine_report`) so that paired extra Cys are not reported as unpaired.
    """
    out = []
    for e in _liab.find_liabilities(seq, contact_residues=contact_residues, rsa=rsa, ssbonds=ssbonds):
        e = dict(e)
        p = e["pos"]
        reg = regions[p - 1] if 1 <= p <= len(regions) else ""
        e["region"] = reg
        e["in_cdr"] = reg.startswith("CDR")
        if e["in_cdr"] and (e["kind"].startswith("N-glycosylation") or e["kind"] in ("Deamidation (NG)", "Isomerisation/cleavage (DG)")) and e["sev"] != "bad":
            e["sev"] = "bad"
            e["reason"] += " Inside a CDR, so it is a risk to binding as well as to stability."
        out.append(e)
    return out


# --------------------------------------------------------------------------------------------- scFv

_LINKER = re.compile(r"[GS]{10,}")


def split_scfv(seq: str, backend: str = "auto") -> Optional[dict]:
    """Split an scFv (VH-linker-VL or VL-linker-VH) at its linker.

    Tries every Gly/Ser-rich stretch of 10 or more residues, numbers the two halves and accepts the first split
    where both halves number and are one heavy and one light chain. Returns ``{"order", "first": (a, b),
    "second": (a, b), "linker": (a, b), "nH"}`` with 1-based inclusive spans, or None.
    """
    s = _clean(seq)
    for m in _LINKER.finditer(s):
        a, b = m.start(), m.end()
        # The linker regex may swallow the closing Ser of FR4 (VSS + GGGGS...): give up to 3 residues back.
        for k in range(4):
            if k and s[a + k - 1] != "S":
                break
            left, right = s[:a + k], s[b:]
            nl = number_chain(left, backend)
            if nl["ok"] and nl["regions"].count("FR4") >= (11 if nl["chain_type"] == "H" else 10):
                break
        else:
            left = s[:a]
            nl = number_chain(left, backend)
        a += len(left) - a
        if len(left) < 90 or len(right) < 90:
            continue
        nr = number_chain(right, backend)
        if not (nl["ok"] and nr["ok"]):
            continue
        if (nl["chain_type"] == "H") == (nr["chain_type"] == "H"):
            continue
        order = "VH-VL" if nl["chain_type"] == "H" else "VL-VH"
        hn = nl if order == "VH-VL" else nr
        nH = hn["domain"][1] - hn["domain"][0] + 1
        return {"order": order, "first": (1, a), "second": (b + 1, len(s)), "linker": (a + 1, b), "nH": nH,
                "_numbering": (nl, nr)}
    return None


# --------------------------------------------------------------------------------------------- profile

def _merge_regions(parts: list, n: int) -> list:
    regions = [FLANK] * n
    for offset, num in parts:
        for i, r in enumerate(num["regions"]):
            regions[offset + i] = r
    return regions


def profile_sequence(seq: str, contact_residues: Optional[Iterable[int]] = None,
                     sg_xyz: Optional[Mapping[int, Sequence[float]]] = None, rsa: Optional[Sequence[float]] = None,
                     backend: str = "auto") -> dict:
    """Full antibody profile of one design sequence: the function the exporter calls.

    ``seq``              the binder sequence as in the PDB (VHH, VH, VL, an scFv with linker, or a Fv with a
                         tag; leaders and tags are allowed).
    ``contact_residues`` 1-based binder residues in contact with the target (enables the CDR share).
    ``sg_xyz``           {1-based position: (x, y, z)} of Cys SG atoms (enables distance-based pairing).
    ``rsa``              per-residue relative SASA for the liabilities.

    Returns ``{"ok", "modality" ("nanobody", "scFv", "Fv-chain"), "backend", "chain_types", "regions"
    (per residue), "region_counts" (residues per region, whole sequence), "cdrs" (list of 3 sequences for a
    single domain, {"H": [...], "L": [...]} for scFv), "spans", "nH" (scFv), "hallmarks", "cysteines",
    "interface" (cdr_interface_share, when contacts are given), "liabilities", "cdr_liabilities" (the subset in
    CDRs)}``. On failure ``ok`` is False with ``error``; nothing raises on odd input.
    """
    s = _clean(seq)
    n = len(s)
    single = number_chain(s, backend)
    sp = split_scfv(s, backend) if _LINKER.search(s) else None
    if sp is None and not single["ok"]:
        return {"ok": False, "error": single.get("error", "no V domain found"), "n": n, "backend": None}
    if contact_residues is not None:
        contact = sorted({int(x) for x in contact_residues})
    else:
        contact = None
    out: dict = {"ok": True, "n": n}
    if sp is None:
        regions = single["regions"][:]
        out.update(modality="nanobody" if single["chain_type"] == "H" else "Fv-chain", backend=single["backend"],
                   chain_types=[single["chain_type"]], spans=dict(single["spans"]),
                   cdrs=[single["cdrs"][c] for c in CDR_NAMES])
        out["hallmarks"] = hallmarks(s, single)
        out["cysteines"] = cysteine_report(s, single, sg_xyz)
        out["cysteine_status"] = out["cysteines"]["status"]
        pairs = [(p["a"], p["b"]) for p in out["cysteines"]["pairs"]]
    else:
        nl, nr = sp["_numbering"]
        off2 = sp["second"][0] - 1
        regions = [FLANK] * n
        for off, num in ((0, nl), (off2, nr)):
            for i, r in enumerate(num["regions"]):
                regions[off + i] = r
        slices = {"first": (nl, 0, s[:sp["first"][1]]), "second": (nr, off2, s[off2:])}
        hk, lk = ("first", "second") if sp["order"] == "VH-VL" else ("second", "first")
        spans = {}
        for tag, k in (("H", hk), ("L", lk)):
            num, off, _ = slices[k]
            for name, (a, b) in num["spans"].items():
                spans[f"{tag}.{name}"] = (a + off, b + off)
        out.update(modality="scFv", backend=nl["backend"], chain_types=[slices[hk][0]["chain_type"], slices[lk][0]["chain_type"]],
                   spans=spans, nH=sp["nH"], order=sp["order"], linker=sp["linker"],
                   cdrs={"H": [slices[hk][0]["cdrs"][c] for c in CDR_NAMES], "L": [slices[lk][0]["cdrs"][c] for c in CDR_NAMES]})
        out["hallmarks"] = hallmarks(slices[hk][2], slices[hk][0])
        out["cysteines"] = {}
        pairs = []
        for tag, k in (("H", hk), ("L", lk)):
            num, off, sub = slices[k]
            rep = cysteine_report(sub, num)
            out["cysteines"][tag] = rep
            pairs += [(p["a"] + off, p["b"] + off) for p in rep["pairs"]]
        out["cysteine_status"] = _worst_cys([r["status"] for r in out["cysteines"].values()])
    counts: dict = {}
    for r in regions:
        counts[r] = counts.get(r, 0) + 1
    out["regions"], out["region_counts"] = regions, counts
    if contact is not None:
        out["interface"] = cdr_interface_share(regions, contact)
    liab = cdr_liabilities(s, regions, contact_residues=contact, rsa=rsa, ssbonds=pairs)
    out["liabilities"] = liab
    out["cdr_liabilities"] = [e for e in liab if e["in_cdr"]]
    return out


_CYS_RANK = {"ok": 0, "extra disulfide": 1, "free cysteine": 2, "canonical pair broken": 3}


def _worst_cys(statuses: Sequence[str]) -> str:
    return max(statuses, key=lambda x: _CYS_RANK.get(x, 0))
