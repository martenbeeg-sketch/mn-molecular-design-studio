"""Sequence and structure liabilities of a designed binder, and its bulk sequence properties.

Everything here is a screen: a flag says "look at this", not "this will fail". Positions are 1-based binder
residues. Severity is "bad" (fix or drop), "warn" (look) or "info" (note). The thresholds sit in the constants
below so a profile can override them later; they are judgement calls, documented in docs/interface.md.

A liability is a dict with the inspector's long keys (``kind``, ``pos``, ``motif``, ``sev``) plus ``reason``
(one line) and ``scope`` ("residue" or "global"). ``to_campaign`` gives the campaign blob's short form
(``k``, ``p``, ``m``, ``s``). The ``kind`` strings keep the wording of the mockup because the pages match on
substrings ("Unpaired", "N-glycosylation", "NG", "DP", "Hydrophobic").

Non-standard residues (X, U, B and so on, any non-letter) never raise: they take part in no motif and count
in no total.
"""
from __future__ import annotations

import re
from typing import Iterable, Mapping, Optional, Sequence

from . import disulfides as D
from . import structure as S

SEVERITIES = ("bad", "warn", "info")
EXPOSED_RSA = 0.25          # relative SASA at or above which a residue counts as exposed
DIBASIC_RSA = 0.40          # dibasic protease sites are only flagged where clearly exposed
HYDROPHOBIC_STRETCH = 6     # consecutive V I L M F W Y
PATCH_WARN_AREA = 500.0     # A^2, the same limit the pages use for "large exposed hydrophobic patch"
PATCH_INFO_AREA = 300.0
CHARGE_INFO = 8.0           # |net charge at pH 7.4|
CHARGE_WARN = 12.0
PH = 7.4

_STRETCH = re.compile(r"(?=([VILMFWY]{%d,}))" % HYDROPHOBIC_STRETCH)


def _clean(seq: str) -> str:
    """Upper case, non-letters become 'X' so positions survive."""
    return "".join(c if c.isalpha() else "X" for c in str(seq).upper())


def _entry(kind: str, pos: int, motif: str, sev: str, reason: str, scope: str = "residue", basis: Optional[str] = None) -> dict:
    e = {"kind": kind, "pos": int(pos), "motif": motif, "sev": sev, "reason": reason, "scope": scope}
    if basis:
        e["basis"] = basis       # what the call rests on; only the disulfide entries carry it
    return e


def _rsa_at(rsa: Optional[Sequence[float]], pos: int) -> Optional[float]:
    if rsa is None or pos - 1 >= len(rsa):
        return None
    v = rsa[pos - 1]
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def sequence_properties(seq: str) -> dict:
    """Bulk properties of a sequence: pI, net charge at pH 7.4, GRAVY and residue counts.

    Only the 20 standard letters count; ``length`` is the number of residues given.
    """
    s = _clean(seq)
    std = S._clean(s)
    return {
        "length": len(s),
        "pI": S.isoelectric_point(std) if std else 0.0,
        "charge": S.net_charge(std, PH) if std else 0.0,
        "gravy": S.gravy(std),
        "n_cys": s.count("C"), "n_met": s.count("M"), "n_trp": s.count("W"), "n_lys": s.count("K"),
    }


def global_flags(seq: str) -> list[dict]:
    """Whole-sequence flags (scope "global", position 1 by convention): unusual net charge."""
    props = sequence_properties(seq)
    z = props["charge"]
    out = []
    if abs(z) >= CHARGE_INFO:
        sev = "warn" if abs(z) >= CHARGE_WARN else "info"
        side = "acidic" if z < 0 else "basic"
        out.append(_entry("Unusual charge", 1, f"{z:+.1f}", sev,
                          f"Net charge {z:+.1f} at pH {PH} (pI {props['pI']:.1f}): a strongly {side} binder can bind "
                          f"non-specifically or behave badly on ion exchange.", "global"))
    return out


def find_liabilities(seq: str, contact_residues: Optional[Iterable[int]] = None,
                     rsa: Optional[Sequence[float]] = None, ssbonds: Optional[Iterable[Sequence[int]]] = None,
                     patches: Optional[Sequence[Mapping]] = None, include_global: bool = False,
                     disulfides: Optional[Mapping] = None) -> list[dict]:
    """Sequence and structure liabilities of one binder.

    ``contact_residues``  1-based binder residues in contact with the target (interface); enables the
                          Lys and Met interface notes.
    ``rsa``               relative SASA per residue (list of length n, 0 to 1); enables the exposure rules
                          for Met, Trp, hydrophobic stretches and dibasic sites. Without it Met and Trp
                          oxidation are reported as info with "exposure unknown".
    ``ssbonds``           disulfide pairs (1-based). When given, every Cys outside a pair is unpaired; when
                          None the parity rule applies (an odd number of Cys leaves one unpaired).
    ``disulfides``        the result of :func:`disulfides.assess` ({basis, pairs, nearest}). Preferred to
                          ``ssbonds``: it also says how the pairs are known. A Cys that the file's record (or, with
                          no record, the model's SG-SG geometry) leaves unpaired is reported; with basis "none"
                          (no record and no sulfur coordinates) the cysteine check is off and nothing is
                          reported or guessed. A Cys with no SG atom in an inferred structure is not judged.
    ``patches``           exposed hydrophobic patches ``[{area, res}]`` from the interface analysis.
    ``include_global``    also return the whole-sequence flags of :func:`global_flags`.

    Returns a list sorted by position then kind.
    """
    s = _clean(seq)
    n = len(s)
    contact = {int(x) for x in contact_residues} if contact_residues is not None else None
    out: list[dict] = []

    # N-glycosylation sequon N-X-S/T, X is not Pro. Matters for mammalian, yeast and insect expression.
    for m in re.finditer(r"N(?=[^P][ST])", s):
        i = m.start()
        out.append(_entry("N-glycosylation sequon", i + 1, s[i:i + 3], "warn",
                          "N-X-S/T motif: N-glycosylated in mammalian, yeast or insect expression, not in E. coli or cell-free."))
    # Deamidation. NG is the fast one, NS is slower.
    for m in re.finditer(r"N(?=G)", s):
        out.append(_entry("Deamidation (NG)", m.start() + 1, "NG", "warn",
                          "Asn-Gly deamidates readily at neutral to basic pH."))
    for m in re.finditer(r"N(?=S)", s):
        out.append(_entry("Deamidation (NS)", m.start() + 1, "NS", "info",
                          "Asn-Ser deamidates slowly; worth knowing for long storage."))
    # Asp isomerisation and acid cleavage.
    for m in re.finditer(r"D(?=G)", s):
        out.append(_entry("Isomerisation/cleavage (DG)", m.start() + 1, "DG", "warn",
                          "Asp-Gly isomerises to iso-Asp through a succinimide, the most prone Asp motif."))
    for m in re.finditer(r"D(?=S)", s):
        out.append(_entry("Isomerisation/cleavage (DS)", m.start() + 1, "DS", "info",
                          "Asp-Ser isomerises slowly."))
    for m in re.finditer(r"D(?=P)", s):
        out.append(_entry("Isomerisation/cleavage (DP)", m.start() + 1, "DP", "info",
                          "Asp-Pro bonds can cleave at low pH."))

    # Cysteines.
    cys = [i + 1 for i, c in enumerate(s) if c == "C"]
    if cys and disulfides is not None:
        out.extend(_disulfide_entries(cys, disulfides))
    elif cys:
        if ssbonds is not None:
            paired = {int(x) for pair in ssbonds for x in pair}
            free = [p for p in cys if p not in paired]
        else:
            free = [cys[-1]] if len(cys) % 2 == 1 else []
        for p in free:
            out.append(_entry("Unpaired cysteine", p, "C", "bad",
                              "A free Cys can form wrong disulfides or dimers; remove it or pair it."))

    # Oxidation of Met and Trp, judged on exposure.
    for letter, kind, exposed_sev in (("M", "Met oxidation", "warn"), ("W", "Trp oxidation", "info")):
        for i, c in enumerate(s):
            if c != letter:
                continue
            r = _rsa_at(rsa, i + 1)
            if r is None:
                out.append(_entry(kind, i + 1, c, "info", f"{letter} oxidises when exposed; exposure unknown here."))
            elif r >= EXPOSED_RSA:
                out.append(_entry(kind, i + 1, c, exposed_sev,
                                  f"Exposed {letter} (relative SASA {r:.2f}) oxidises in storage or under light."))

    # Hydrophobic stretches: long runs of V I L M F W Y. A core run in a helix is normal, so severity
    # follows exposure when it is known.
    for m in _STRETCH.finditer(s):
        i = m.start()
        run = m.group(1)
        rs = [_rsa_at(rsa, p + 1) for p in range(i, i + len(run))]
        known = [x for x in rs if x is not None]
        if not known:
            sev, why = "warn", "Long hydrophobic run: aggregation and sticky binding risk if exposed."
        elif sum(known) / len(known) >= EXPOSED_RSA:
            sev, why = "warn", "Long hydrophobic run that is exposed: aggregation and sticky binding risk."
        else:
            sev, why = "info", "Long hydrophobic run, mostly buried: usually the core, but check it is packed."
        out.append(_entry("Hydrophobic stretch", i + 1, run, sev, why))

    # Exposed hydrophobic surface patches from the structure.
    for p in patches or []:
        area = float(p.get("area", 0.0))
        res = sorted(int(x) for x in p.get("res", []))
        if not res or area < PATCH_INFO_AREA:
            continue
        sev = "warn" if area >= PATCH_WARN_AREA else "info"
        motif = "".join(s[r - 1] for r in res if 1 <= r <= n)
        out.append(_entry("Exposed hydrophobic patch", res[0], motif[:12], sev,
                          f"Connected exposed hydrophobic patch of {area:.0f} A2 over {len(res)} residues: sticky surface."))

    # Protease sites. Furin consensus R-X-K/R-R; dibasic pairs only where clearly exposed.
    for m in re.finditer(r"R(?=[A-Z][KR]R)", s):
        i = m.start()
        r = _rsa_at(rsa, i + 1)
        sev = "warn" if r is None or r >= EXPOSED_RSA else "info"
        out.append(_entry("Protease site (furin)", i + 1, s[i:i + 4], sev, "R-X-K/R-R is a furin consensus; cleaved during secretion in mammalian cells."))
    for m in re.finditer(r"(?=([KR][KR]))", s):
        i = m.start()
        r1, r2 = _rsa_at(rsa, i + 1), _rsa_at(rsa, i + 2)
        if r1 is not None and r2 is not None and min(r1, r2) >= DIBASIC_RSA:
            out.append(_entry("Protease site (dibasic)", i + 1, s[i:i + 2], "info",
                              "Exposed Lys/Arg pair: trypsin-like and PC-type proteases cut here."))

    # Interface Lys and Met: BindCraft caps both (default 3), Lys is a trypsin and conjugation hotspot.
    if contact is not None:
        for i, c in enumerate(s):
            if (i + 1) in contact and c == "K":
                out.append(_entry("Interface Lys", i + 1, "K", "info",
                                  "Lys in the interface: lost to amine conjugation or proteolysis changes binding."))
            elif (i + 1) in contact and c == "M":
                out.append(_entry("Interface Met", i + 1, "M", "info",
                                  "Met in the interface: oxidation would change the binding surface."))

    if include_global:
        out.extend(global_flags(s))
    out.sort(key=lambda e: (e["pos"], e["kind"]))
    return out


def _disulfide_entries(cys: Sequence[int], ds: Mapping) -> list[dict]:
    """Unpaired-Cys liabilities, and the note on loosely paired ones, from a :func:`disulfides.assess` result."""
    basis, nearest = ds.get("basis"), ds.get("nearest") or {}
    if basis not in (D.BASIS_RECORDED, D.BASIS_INFERRED):
        return []                                    # no record and no coordinates: the check has no data, so it is off
    pairs = ds.get("pairs") or []
    paired = {int(x) for p in pairs for x in (p["a"], p["b"])}
    out = []
    if basis == D.BASIS_RECORDED:
        judged = list(cys)                           # the file speaks for every Cys
        why = "The file records disulfides in this chain and none of them involves this Cys."
    else:
        judged = [p for p in cys if p in nearest]    # a Cys without an SG atom cannot be placed, so it is not judged
        why = "The file records no disulfides, so pairing is read from the model's geometry."
    for p in judged:
        if p in paired:
            continue
        gap = nearest.get(p)
        if basis == D.BASIS_INFERRED:
            near = (" No other Cys sulfur is in the model." if gap is None else
                    f" The nearest other Cys sulfur is {gap:.1f} A away, beyond the {D.LOOSE_MAX:.1f} A pairing limit.")
            reason = "A free Cys can form wrong disulfides or dimers; remove it or pair it." + near + " " + why
        else:
            reason = "A free Cys can form wrong disulfides or dimers; remove it or pair it. " + why
        out.append(_entry("Unpaired cysteine", p, "C", "bad", reason, basis=basis))
    if basis == D.BASIS_INFERRED:
        for pr in pairs:
            if pr["band"] == "loose":
                out.append(_entry("Disulfide geometry loose", pr["a"], "C", "info",
                                  f"Cys{pr['a']} and Cys{pr['b']} are paired by geometry at {pr['d']:.1f} A, looser than a closed "
                                  f"bond (about 2 A) and inside the {D.CLOSE_MAX:.1f} to {D.LOOSE_MAX:.1f} A band where the model is "
                                  "least sure. The file records no disulfides, so this rests on the model alone.", basis=basis))
    return out


def to_campaign(liabs: Iterable[Mapping]) -> list[dict]:
    """Campaign blob form of a liability list: ``{k, p, m, s}``, plus ``b`` (the basis) on disulfide entries."""
    return [{"k": e["kind"], "p": int(e["pos"]), "m": e["motif"], "s": e["sev"], **({"b": e["basis"]} if e.get("basis") else {})}
            for e in liabs]


def to_inspector(liabs: Iterable[Mapping]) -> list[dict]:
    """Inspector blob form: ``{kind, pos, motif, sev}`` plus ``reason``, and ``basis`` on disulfide entries."""
    return [{"kind": e["kind"], "pos": int(e["pos"]), "motif": e["motif"], "sev": e["sev"], "reason": e.get("reason", ""),
             **({"basis": e["basis"]} if e.get("basis") else {})} for e in liabs]


def severity_counts(liabs: Iterable[Mapping]) -> dict:
    """``{bad, warn, info}`` counts."""
    out = {k: 0 for k in SEVERITIES}
    for e in liabs:
        out[e["sev"]] = out.get(e["sev"], 0) + 1
    return out
