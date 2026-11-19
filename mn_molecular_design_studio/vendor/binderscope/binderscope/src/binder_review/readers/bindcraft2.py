"""BindCraft 2 reader.

Verified against one real campaign: de novo (`binder` modality) output of BindCraft 2 v1.0.3, from its
trajectory table through the refolded candidates to the ranked shortlist (docs/bindcraft2_reader.md says what was checked and how).
Anything outside that (other versions, VHH, scFv, Fab, cyclic peptides, several targets) is read on the same
assumptions but is untested, and such a run keeps experimental=True with a warning that names the reason.
The reader never raises on columns or metadata keys it does not know; unknown columns land in
DesignRecord.extras, and what the reader had to assume is written to RunData.notes.

Layout read (a campaign folder)::

    1_Trajectories/!_Trajectories.csv   one row per trajectory, finished or not (hash, terminated)
    1_Trajectories/<design>/            only <design>_losses.csv unless archive_trajectories is true
    2_Refolded/!_Refolded.csv           one row per refolded candidate (outcome passed / rejected)
    2_Refolded/Complexes/<design>_candidate<N>.cif                bound prediction (chain A target, the binder in the chains
                                                                  after it: one, or several for an scFv or a homo-oligomer)
    2_Refolded/BinderMonomer/<design>_candidate<N>_monomer.cif    the unbound binder prediction
    3_Ranked/!_Ranked.csv               the final shortlist (rank, ids ending in _seq0) + <design>_seq0.cif
    campaign_metadata.json, summary.csv, .campaign_state.json

Mapping to the shared model:
  * the backbone (family) key is `hash`, the 16-hex trajectory hash that all three tables carry. Neither
    !_Refolded.csv nor !_Ranked.csv has a `trajectory` column. The backbone id is the trajectory's design name.
  * accepted design = a refolded row whose `outcome` is passed (and that names no failed filter). A ranked
    row is the same molecule as one refolded row (matched on hash and Binder_Sequence, the ids differ: _seq0
    against _candidate<N>); it gives that design its rank and its 3_Ranked structure and creates no second design.
  * metrics keep the BindCraft 1 stem and scale where the quantity is the same (METRIC_COLUMNS); columns whose
    definition differs (Interface_BuriedArea, the clash counts, ...) stay in extras under their own names.
  * stems BindCraft 2 does not report (dG, ShapeComplementarity, PackStat, n_Interface*, ...) are simply
    absent, so RunData.caps() is false for them.
  * Interface_Binder_Residues is "<amino acid><number>" ("P1,A25,H26"), not <chain><number>. The letters are
    checked against Binder_Sequence.
  * binder residue numbers are 1-based in v1.0.3 (the structures say so). Older releases are documented to
    write them zero-based; that is handled from the structures but has not been seen (see one_based_copy).
"""

from __future__ import annotations

import ast
import json
import math
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

import pandas as pd

from . import ReaderError
from ..config import FilterFailure
from ..model import (
    METRIC_STEMS, BackboneRecord, DesignRecord, FunnelCounts, InputFile, RunData, RunSettings,
    TargetConstruct,
)
from .bindcraft1 import (
    _gaps, _missing_ranges, _one_letter, _reference_construct, parse_hotspots, sha256_file,
)

NAME = "bindcraft2"

#: What the reader has been checked against. A campaign outside this is still read, but is flagged experimental.
VERIFIED_VERSIONS: tuple[tuple[int, int, int], ...] = ((1, 0, 3),)
VERIFIED_MODALITIES: tuple[str, ...] = ("de novo",)
VERIFIED_SETTINGS_MODALITY = "binder"

TRAJECTORY_DIR = "1_Trajectories"
REFOLDED_DIR = "2_Refolded"
RANKED_DIR = "3_Ranked"
STAGE_DIRS = (TRAJECTORY_DIR, REFOLDED_DIR, RANKED_DIR)
TABLES = {TRAJECTORY_DIR: "!_Trajectories.csv", REFOLDED_DIR: "!_Refolded.csv", RANKED_DIR: "!_Ranked.csv"}
METADATA_GLOB = "campaign_metadata*.json"
STRUCTURE_SUFFIXES = (".cif", ".mmcif")
MONOMER_SUFFIX = "_monomer"

#: Last release that is documented to write binder residue numbers zero-based (1.0.0 and 1.0.1). Not seen on
#: a real campaign: v1.0.3 is one-based, which the structures of that campaign confirm.
LAST_ZERO_BASED_VERSION = (1, 0, 1)

#: BindCraft 1 stems BindCraft 2 does not report; downstream checks that need them switch off.
ABSENT_STEMS: tuple[str, ...] = (
    "dG", "ShapeComplementarity", "PackStat", "n_InterfaceHbonds", "InterfaceHbondsPercentage",
    "n_InterfaceUnsatHbonds", "InterfaceUnsatHbondsPercentage", "Binder_Energy_Score", "Hotspot_RMSD",
    "i_pLDDT", "dG/dSASA", "Relaxed_Clashes", "MPNN_score", "MPNN_seq_recovery", "pAE", "Binder_pTM", "Binder_pAE",
)

#: BindCraft 1 stems for which BindCraft 2 has a column that is NOT the same quantity, so nothing is mapped.
#: stem -> (BindCraft 2 column(s), why not).
UNMAPPED_STEMS: dict[str, tuple[str, str]] = {
    "dSASA": ("Interface_BuriedArea",
              "binder-side buried area only; BindCraft 1's dSASA is the buried area of both partners"),
    "Interface_SASA_%": ("Interface_BuriedArea_Fraction",
                         "binder-side buried area over the free binder's surface; BindCraft 1 divides the buried area of "
                         "both sides by the binder's area in the complex"),
    "Unrelaxed_Clashes": ("Backbone_Clashes, All_Atom_Clashes",
                          "interchain only, CA atoms only (Backbone_Clashes) or all atoms within 2.5 A (All_Atom_Clashes); "
                          "one stem cannot hold two counts and neither is BindCraft 1's count"),
}

#: Source column -> ((stem under which it is stored in DesignRecord.metrics, factor), ...). The first entry is the
#: main one. The factor turns BindCraft 2's fractions into the percent scale of the BindCraft 1 stems that end in %.
#: Columns that are not here are not mapped and stay in extras (see UNMAPPED_STEMS and docs/bindcraft2_reader.md).
METRIC_COLUMNS: dict[str, tuple[tuple[str, float], ...]] = {
    "i_pTM": (("i_pTM", 1.0),), "pLDDT": (("pLDDT", 1.0),), "pTM": (("pTM", 1.0),), "i_pAE": (("i_pAE", 1.0),),
    # the binder predicted alone: stored under its own name and as BindCraft 1's Binder_pLDDT ("binder predicted alone").
    # The campaign's FILTER keeps its own name (see _failures); only this value is shared.
    "Unbound_Binder_pLDDT": (("Unbound_Binder_pLDDT", 1.0), ("Binder_pLDDT", 1.0)),
    "i_pDAE": (("i_pDAE", 1.0),), "Target_pLDDT": (("Target_pLDDT", 1.0),),
    "SS_pLDDT": (("ss_pLDDT", 1.0),),
    "Interface_Residues": (("n_InterfaceResidues", 1.0),),
    "Surface_Hydrophobicity": (("Surface_Hydrophobicity", 1.0),),
    "Interface_Hydrophobicity": (("Interface_Hydrophobicity", 100.0),),
    "Binder_Helix_Fraction": (("Binder_Helix%", 100.0),),
    "Binder_BetaSheet_Fraction": (("Binder_BetaSheet%", 100.0),),
    "Binder_Loop_Fraction": (("Binder_Loop%", 100.0),),
    "Interface_Helix_Fraction": (("Interface_Helix%", 100.0),),
    "Interface_BetaSheet_Fraction": (("Interface_BetaSheet%", 100.0),),
    "Interface_Loop_Fraction": (("Interface_Loop%", 100.0),),
    "Target_RMSD": (("Target_RMSD", 1.0),), "Binder_RMSD": (("Binder_RMSD", 1.0),),
    "Off_Paratope_Contact_Fraction": (("Off_Paratope_Contact_Fraction", 1.0),),
    "Scaffold_Framework_RMSD": (("Scaffold_Framework_RMSD", 1.0),),
    "Scaffold_Sequence_Retained_Fraction": (("Scaffold_Sequence_Retained_Fraction", 1.0),),
    "Interdomain_Contact_Fraction": (("Interdomain_Contact_Fraction", 1.0),),
}
ANTIBODY_COLUMNS = ("Off_Paratope_Contact_Fraction", "Scaffold_Framework_RMSD",
                    "Scaffold_Sequence_Retained_Fraction", "Interdomain_Contact_Fraction")
#: Columns on the 0-1 pLDDT scale. A table where one exceeds 1.5 is read as 0-100 and divided by 100; real v1.0.3
#: output never does that, the guard is for output that has not been seen.
_PLDDT_SOURCES = ("pLDDT", "Unbound_Binder_pLDDT", "Target_pLDDT", "SS_pLDDT")

_FIELD_COLUMNS = {"rank", "trajectory", "design", "length", "binder_sequence", "interface_binder_residues",
                  "failed_filters"}
_LOWER_METRIC = {k.lower(): k for k in METRIC_COLUMNS}
#: Columns seen on the real campaign that the reader keeps in extras on purpose (so they are not "unknown").
_KNOWN_EXTRA = {c.lower() for c in (
    "outcome", "hash", "Interface_Target_Residues", "Timing", "bindcraft_version", "terminated", "autotuned",
    "settings_core", "settings_modality", "Interface_BuriedArea", "Interface_BuriedArea_Fraction", "Backbone_Clashes",
    "All_Atom_Clashes", "Hotspot_Contact_Fraction", "Off_Epitope_Contact_Fraction", "Epitope_Residues_Contacted",
    "Binder_Disulfides", "Binder_Mass_kDa", "Binder_pI", "Binder_Net_Charge", "Binder_Cysteines",
    "Binder_Free_Cysteines", "Binder_Extinction")}
_KNOWN_COLUMNS = _FIELD_COLUMNS | set(_LOWER_METRIC) | _KNOWN_EXTRA
#: What a refolded table needs for the reader's assumptions to hold (checked once per run).
_REQUIRED_REFOLDED = ("design", "hash", "outcome", "failed_filters", "Binder_Sequence")

#: Direction of the filters BindCraft 2 documents, used only when the metadata gives a bare number or no `higher`.
_HIGHER_BETTER = {"i_ptm": True, "plddt": True, "ptm": True, "i_pae": False, "i_pdae": True,
                  "unbound_binder_plddt": True, "target_plddt": True}

NAME_RE = re.compile(r"^(?P<campaign>.+)_(?P<modality>[a-z][a-z0-9]*)_l(?P<length>\d+)_(?P<hash>[0-9a-fA-F]{16})"
                     r"(?:_(?P<kind>candidate|seq)(?P<n>\d+))?$")
_SUFFIX_RE = re.compile(r"_(?:candidate|seq)\d+$")
#: Campaign modalities, as words. Order matters where one word contains another:
#: the lookup takes the first key found in the metadata's modality, so
#: "cyclicpeptide" has to be tried before "peptide" and "largebinder" before
#: "binder", or a cyclic peptide campaign reads as a linear one.
#:
#: Every modality BindCraft 2 ships has a word here. Without one, a campaign
#: fell through to the column heuristics, where two or more binder chains means
#: scFv -- so a homo-oligomer campaign, whose binder is several copies of one
#: chain, was reported as an scFv.
_MODALITY_WORDS = {
    "denovo": "de novo", "de-novo": "de novo", "vhh": "VHH", "nanobody": "VHH", "scfv": "scFv", "fab": "Fab",
    "cyclicpeptide": "cyclic peptide", "cyclic": "cyclic peptide", "cyc": "cyclic peptide",
    "peptide": "peptide", "arp": "ARP", "ankyrin": "ARP",
    "largebinder": "large binder", "homooligomer": "homo-oligomer", "multidomain": "multidomain",
}
_ACCEPT_WORDS = ("accept", "pass", "success", "ranked", "final")


# ---------------------------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------------------------


def parse_name(stem: str) -> Optional[dict[str, Any]]:
    """Split a BindCraft 2 design name <campaign>_<modality>_l<len>_<hash16>[_candidate<N> | _seq<N>].

    Returns campaign, modality (and the raw token), length, hash, and for a refolded candidate or a ranked design
    `kind` ("candidate" / "seq") with its number `n`; None if the name does not follow the pattern. A trajectory
    name has no suffix. Pass the design id, not a file name with a rank prefix or a _monomer suffix.
    """
    m = NAME_RE.match(stem)
    if not m:
        return None
    mod = m["modality"].lower()
    return {"campaign": m["campaign"], "modality_token": mod, "modality": _MODALITY_WORDS.get(mod),
            "length": int(m["length"]), "hash": m["hash"].lower(), "kind": m["kind"],
            "n": int(m["n"]) if m["n"] is not None else None}


def _child(root: Path, name: str) -> Optional[Path]:
    if not root.is_dir():
        return None
    p = root / name
    if p.exists():
        return p
    for q in root.iterdir():
        if q.name.lower() == name.lower():
            return q
    return None


def _table(root: Path, stage: str) -> Optional[Path]:
    d = _child(root, stage)
    if d is None or not d.is_dir():
        return None
    p = _child(d, TABLES[stage])
    if p is not None:
        return p
    hits = sorted(d.glob("!_*.csv"))
    return hits[0] if hits else None


def _metadata_files(root: Path) -> list[Path]:
    return sorted(p for p in root.glob(METADATA_GLOB) if p.is_file()) if root.is_dir() else []


def detect(path: str | Path) -> bool:
    """True when the folder looks like a BindCraft 2 campaign.

    Needs two of the three stage folders, or one stage folder plus a !_*.csv table, or the metadata file
    plus a stage folder. A folder with the three BindCraft 1 stats CSVs is never BindCraft 2.
    """
    root = Path(path)
    if not root.is_dir():
        return False
    if all(_child(root, n) is not None for n in ("final_design_stats.csv", "mpnn_design_stats.csv", "trajectory_stats.csv")):
        return False
    folders = [s for s in STAGE_DIRS if (d := _child(root, s)) is not None and d.is_dir()]
    tables = [s for s in folders if _table(root, s) is not None]
    meta = bool(_metadata_files(root))
    return len(folders) >= 2 or (len(folders) >= 1 and (bool(tables) or meta))


def _num(x: Any) -> Optional[float]:
    if x is None:
        return None
    try:
        v = float(str(x).strip())
    except ValueError:
        return None
    return None if math.isnan(v) else v


def _txt(x: Any) -> str:
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else str(x).strip()


def _extra_value(x: Any) -> Any:
    t = _txt(x)
    if t == "":
        return None
    v = _num(t)
    return v if v is not None else t


def _aggregate(items: list[str], what: str, limit: int = 4) -> str:
    shown = ", ".join(items[:limit]) + (f", ... ({len(items)} in all)" if len(items) > limit else "")
    return f"{what}: {shown}"


def read_table(path: Path) -> pd.DataFrame:
    """Read one !_*.csv with every cell as text (formats are not trusted)."""
    return pd.read_csv(path, dtype=str, keep_default_na=False, skipinitialspace=True)


def split_failed_filters(value: Any) -> list[str]:
    """Names in a failed_filters cell.

    BindCraft 2 v1.0.3 joins them with a comma ("i_pTM,pLDDT,Unbound_Binder_pLDDT"; empty for a passed design).
    Splits on , ; | and newlines and accepts a Python or JSON list, and drops empty markers (none, nan, pass, -).
    Whitespace is not a delimiter: BindCraft writes "<name> (not measured)" for a mandatory filter whose metric it
    could not compute, and that stays one entry.
    """
    t = _txt(value)
    if not t:
        return []
    if t[0] in "[(":
        for loader in (json.loads, ast.literal_eval):
            try:
                seq = loader(t)
            except (ValueError, SyntaxError):
                continue
            if isinstance(seq, (list, tuple)):
                t = ";".join(str(s) for s in seq)
                break
        else:
            t = t.strip("[]()")
    out = []
    for part in re.split(r"[;,|\n]+", t):
        p = part.strip().strip("'\"").strip()
        if p and p.lower() not in {"none", "nan", "null", "pass", "passed", "-", "[]"}:
            out.append(p)
    return out


_RES_TOKEN = re.compile(r"^([A-Za-z]?)(-?\d+)$")


def parse_residue_list(text: Any) -> tuple[list[list[tuple[str, int]]], list[str]]:
    """'P1,A25,H26' -> ([[('P', 1), ('A', 25), ('H', 26)]], []).

    BindCraft 2 writes each residue as its amino-acid letter and its number in the written structure. A '/'
    separates the chains of a multi-chain binder and gives one list per chain. A token with no letter ('25') gets
    ''. Second item: the tokens that did not parse (they are not in the first).
    """
    t = _txt(text).strip("[]()")
    if not t:
        return [], []
    parts: list[list[tuple[str, int]]] = []
    bad: list[str] = []
    for chunk in t.split("/"):
        items: list[tuple[str, int]] = []
        for tok in re.split(r"[,;\s|]+", chunk):
            tok = tok.strip("'\"")
            if not tok:
                continue
            m = _RES_TOKEN.match(tok)
            if m:
                items.append((m.group(1).upper(), int(m.group(2))))
            else:
                bad.append(tok)
        parts.append(items)
    return parts, bad


def interface_residues_from(text: Any, sequence: str, chains: list[str], shift: int = 0) -> tuple[list[str], int, int, int]:
    """Binder interface residues as ["B25", ...] (chain letter plus 1-based number) from a residue-list cell.

    chains    binder chain names in the order of the '/' parts (the structure's, else B, C, ...).
    shift     added to the written numbers to make them 1-based: 1 for a zero-based release, else 0.
    The letter of every token is compared with `sequence` at the position its number names (written number minus
    1 - shift, since the first residue is number 1 - shift in the file). Returns (residues, number of letters
    checked, number that disagree with the sequence, number of tokens that did not parse).
    """
    parts, bad = parse_residue_list(text)
    seqs = sequence.split("/") if sequence else []
    origin = 1 - shift
    out: list[str] = []
    checked = wrong = 0
    for ci, items in enumerate(parts):
        chain = chains[ci] if ci < len(chains) else chr(ord("B") + ci)
        seq = seqs[ci] if ci < len(seqs) else ""
        for letter, num in items:
            if letter and seq:
                checked += 1
                idx = num - origin
                if not (0 <= idx < len(seq)) or seq[idx] != letter:
                    wrong += 1
            out.append(f"{chain}{num + shift}")
    return out, checked, wrong, len(bad)


def _version_tuple(text: Any) -> Optional[tuple[int, int, int]]:
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", _txt(text))
    return (int(m.group(1)), int(m.group(2)), int(m.group(3) or 0)) if m else None


def _clean_version(text: Any) -> Optional[str]:
    """'BindCraft 2 v1.0.3' -> '1.0.3' (the stamp already carries the tool's name); text without a number is kept."""
    t = _txt(text)
    if not t:
        return None
    m = re.search(r"v?(\d+\.\d+(?:\.\d+)*)", t)
    return m.group(1) if m else t


def _find_key(obj: Any, names: Iterable[str], types: tuple = (str, int, float, list, dict)) -> Any:
    """First value, depth first, whose key (case-insensitive, - and space read as _) is one of names."""
    want = [n.lower() for n in names]
    stack = [obj]
    found: dict[str, Any] = {}
    while stack:
        cur = stack.pop(0)
        if isinstance(cur, dict):
            for k, v in cur.items():
                key = str(k).lower().replace("-", "_").replace(" ", "_")
                if key in want and isinstance(v, types) and key not in found and v not in ("", None):
                    found[key] = v
                if isinstance(v, (dict, list)):
                    stack.append(v)
        elif isinstance(cur, list):
            stack.extend(x for x in cur if isinstance(x, (dict, list)))
    for n in want:
        if n in found:
            return found[n]
    return None


# ---------------------------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------------------------


def _threshold(value: Any) -> Optional[float]:
    """A finite threshold, or None. BindCraft writes "no limit" as the string "Infinity" (JSON has no infinity)."""
    v = _num(value)
    return v if v is not None and math.isfinite(v) else None


def _is_infinite(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower().lstrip("+-") in ("inf", "infinity")
    return isinstance(value, float) and math.isinf(value)


def _parse_filters(obj: Any) -> tuple[list[dict[str, Any]], list[str]]:
    """The campaign's filters from the metadata: (filters with a finite threshold, names of the unlimited ones).

    Real shape (settings.filters of campaign_metadata.json)::

        {name: {"higher": bool, "threshold": number | "Infinity", "mandatory": bool}}

    higher says which side passes (value >= threshold if true, <= if false; the threshold itself passes). An
    infinite threshold means no limit, so it is dropped (and returned by name), not read as a threshold.
    `mandatory` is BindCraft's flag: a mandatory filter fails a design whose metric could not be measured; a
    missing key means true (the default in bindcraft/filters.py, and how this campaign's metadata is written: only
    the non-binding filters carry `"mandatory": false`). Each entry of the first list has name, op, threshold,
    higher, mandatory and source. Bare numbers and {min | max | threshold, op} are still understood, for layouts
    that have not been seen.
    """
    out: list[dict[str, Any]] = []
    unlimited: list[str] = []
    if isinstance(obj, list):
        merged: dict[str, Any] = {}
        for item in obj:
            if isinstance(item, dict) and "name" in item:
                merged[str(item["name"])] = {k: v for k, v in item.items() if k != "name"}
        obj = merged
    if not isinstance(obj, dict):
        return out, unlimited
    for name, spec in obj.items():
        low = str(name).lower()
        higher: Optional[bool] = None
        raw_thr: Any = None
        mandatory = True
        if isinstance(spec, dict):
            s = {str(k).lower(): v for k, v in spec.items()}
            if isinstance(s.get("higher"), bool):
                higher = s["higher"]
            if isinstance(s.get("mandatory"), bool):
                mandatory = s["mandatory"]
            if "threshold" in s or "value" in s:
                raw_thr = s.get("threshold", s.get("value"))
            if higher is None and "min" in s:
                raw_thr, higher = s["min"], True
            elif higher is None and "max" in s:
                raw_thr, higher = s["max"], False
            if higher is None:
                op = str(s.get("op", s.get("operator", ""))).strip()
                higher = True if op in (">=", ">") else False if op in ("<=", "<") else None
        else:
            raw_thr = spec
        if higher is None:
            higher = _HIGHER_BETTER.get(low)
        if higher is None:
            continue
        thr = _threshold(raw_thr)
        if thr is None:
            if _is_infinite(raw_thr):
                unlimited.append(str(name))
            continue
        out.append({"name": str(name), "threshold": thr, "op": ">=" if higher else "<=", "higher": higher,
                    "mandatory": mandatory, "source": "campaign_metadata"})
    return out, unlimited


def _scaffold_spec(raw: Any) -> Optional[dict[str, Any]]:
    """The scaffold edit spec of a campaign built on a fixed scaffold (VHH, scFv, Fab, ankyrin repeat), or None.

    BindCraft 2 records it as the params of the Scaffold_Sequence_Retained_Fraction filter::

        settings.filters.Scaffold_Sequence_Retained_Fraction.params = {"scaffold": "<path to the scaffold file>",
                                                                       "scaffold_edits": "A27-33,A35*,A52-57(5-7),..."}

    The filter's own threshold and mandatory flag come with it because the paratope analysis needs the retained fraction
    to have been enforced at 1.0 (docs/paratope.md). Only what the file says is kept; nothing is parsed here.
    """
    flt = ((raw.get("settings") or {}).get("filters") or {}) if isinstance(raw, dict) else {}
    if not isinstance(flt, dict):
        return None
    for name, spec in flt.items():
        if str(name).lower() != "scaffold_sequence_retained_fraction" or not isinstance(spec, dict):
            continue
        params = spec.get("params")
        edits = params.get("scaffold_edits") if isinstance(params, dict) else None
        if not isinstance(edits, str) or not edits.strip():
            return None
        scaffold = params.get("scaffold")
        out: dict[str, Any] = {"scaffold": scaffold if isinstance(scaffold, str) else None, "scaffold_edits": edits,
                               "threshold": _threshold(spec.get("threshold")),
                               "higher": spec.get("higher") if isinstance(spec.get("higher"), bool) else None,
                               "mandatory": spec.get("mandatory") if isinstance(spec.get("mandatory"), bool) else True,
                               "source": "settings.filters.Scaffold_Sequence_Retained_Fraction.params"}
        return out
    return None


def _read_metadata(root: Path, notes: list[str], warnings: list[str], files: list[InputFile],
                   explicit: Optional[Path] = None) -> dict[str, Any]:
    """Resolved settings and provenance from campaign_metadata*.json, or an explicit file (all keys optional)."""
    out: dict[str, Any] = {"raw": {}}
    paths = [explicit] if explicit is not None else _metadata_files(root)
    if explicit is not None and not explicit.is_file():
        warnings.append(f"Metadata file {explicit} not found; settings are unknown.")
        return out
    if not paths:
        notes.append("No campaign_metadata.json: settings, hotspots, version and target hash are unknown.")
        return out
    path = paths[0]
    files.append(InputFile("settings", str(path.resolve()), sha256_file(path), path.stat().st_size))
    if len(paths) > 1:
        notes.append(f"Several metadata files ({', '.join(p.name for p in paths)}); read {path.name}.")
    try:
        raw = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        warnings.append(f"{path.name} could not be parsed ({exc}); settings are unknown.")
        return out
    if not isinstance(raw, dict):
        warnings.append(f"{path.name} is not a JSON object; settings are unknown.")
        return out
    out["raw"] = raw
    out["path"] = str(path.resolve())
    out["campaign"] = _find_key(raw, ("campaign_name", "campaign", "name"), (str,))
    ver = _find_key(raw, ("bindcraft_version", "version", "tool_version", "bindcraft2_version"), (str, int, float))
    out["version"] = _clean_version(ver)
    rev = _find_key(raw, ("revision", "git_revision", "git_commit", "commit", "git_sha"), (str,))
    out["revision"] = _txt(rev) or None
    out["modality"] = _find_key(raw, ("modality", "design_mode", "mode", "binder_type", "binder_modality"), (str,))
    hs = _find_key(raw, ("target_hotspots", "hotspots", "hotspot_residues", "hotspot"), (str, list, int))
    if isinstance(hs, list):
        hs = ",".join(str(x) for x in hs)
    out["hotspots"] = _txt(hs) or None
    out["target_hash"] = _find_key(raw, ("target_hash", "target_file_hash", "target_sha256", "target_md5",
                                         "target_pdb_hash", "input_hash", "pdb_hash", "target_file_sha256", "file_hash", "file_sha256"), (str,))
    out["target_file"] = _find_key(raw, ("target_file", "target_pdb", "target_path", "input_pdb", "target"), (str,))
    out["chains"] = _find_key(raw, ("target_chains", "chains", "target_chain"), (str, list))
    lens = _find_key(raw, ("lengths", "binder_lengths", "length_range", "binder_length"), (list, str, int))
    if isinstance(lens, list) and len(lens) >= 2 and all(_num(x) is not None for x in lens[:2]):
        out["length_range"] = (int(_num(lens[0])), int(_num(lens[1])))
    n_final = _find_key(raw, ("number_of_final_designs", "n_final_designs", "final_designs", "num_final_designs"), (int, float, str))
    out["final_designs"] = int(_num(n_final)) if _num(n_final) is not None else None
    filters, unlimited = _parse_filters(_find_key(raw, ("filters", "filter_thresholds", "thresholds"), (dict, list)))
    # BindCraft judges pLDDT (binder, in the complex) against the setting min_plddt_final, which is not an entry of
    # `filters` (bindcraft/filters.py, design_stage_filters). It is the filter behind "pLDDT" in failed_filters.
    min_plddt = _threshold(_find_key(raw, ("min_plddt_final",), (int, float, str)))
    if min_plddt is not None and not any(f["name"].lower() == "plddt" for f in filters):
        filters.append({"name": "pLDDT", "threshold": min_plddt, "op": ">=", "higher": True, "mandatory": True,
                        "source": "settings.min_plddt_final"})
    out["filters"] = filters
    out["filters_unlimited"] = unlimited
    out["scaffold"] = _scaffold_spec(raw)
    arch = _find_key(raw, ("archive_trajectories",), (bool, int, str))
    out["archive_trajectories"] = arch if isinstance(arch, bool) else None
    return out


# ---------------------------------------------------------------------------------------------
# Structures
# ---------------------------------------------------------------------------------------------


def _structure_index(directory: Optional[Path]) -> dict[str, Path]:
    """File stem -> path for every mmCIF below a stage folder."""
    out: dict[str, Path] = {}
    if directory is None or not directory.is_dir():
        return out
    for p in sorted(directory.rglob("*")):
        if p.suffix.lower() in STRUCTURE_SUFFIXES and p.is_file():
            out.setdefault(p.stem, p.resolve())
    return out


def _stem_matches(stem: str, design: str) -> bool:
    """The file stem is the design id, or the id behind a rank prefix ('1_<id>'). Never a substring match:
    <design>_candidate1 is a prefix of <design>_candidate10."""
    return stem == design or re.fullmatch(r"\d+_" + re.escape(design), stem) is not None


def _match_structure(index: Mapping[str, Path], design: str) -> Optional[Path]:
    if design in index:
        return index[design]
    for stem, path in index.items():
        if _stem_matches(stem, design):
            return path
    return None


def inspect_structure(path: Path, n_binder_chains: int = 1) -> dict[str, Any]:
    """Light check of one BindCraft 2 mmCIF.

    Chains are lettered by position with the targets first, so the last n_binder_chains chains are the
    binder (chain A target, chain B binder in v1.0.3). Returns chain names, target numbers and sequence (first
    target chain), binder sequences, the first residue number of each binder chain, mean binder B-factor and
    whether hydrogens are present.
    """
    import gemmi

    st = gemmi.read_structure(str(path))
    model = st[0]
    names = [c.name for c in model]
    nb = max(1, min(n_binder_chains, len(names) - 1)) if len(names) > 1 else 0
    target_names, binder_names = names[: len(names) - nb], names[len(names) - nb:]
    out: dict[str, Any] = {
        "chains": names, "target_chains": target_names, "binder_chains": binder_names,
        "hydrogens": any(at.element.name == "H" for c in model for r in c for at in r),
    }
    if target_names:
        a = model[target_names[0]]
        out["target_numbers"] = tuple(r.seqid.num for r in a)
        out["target_seq"] = "".join(_one_letter(r.name) for r in a)
    seqs, firsts, bs = [], [], []
    for n in binder_names:
        ch = model[n]
        seqs.append("".join(_one_letter(r.name) for r in ch))
        firsts.append(ch[0].seqid.num if len(ch) else None)
        bs.extend(at.b_iso for r in ch for at in r)
    out["binder_seqs"] = seqs
    out["binder_first_numbers"] = firsts
    out["binder_mean_b"] = sum(bs) / len(bs) if bs else None
    return out


def _read_bindcraft_block(path: Path) -> dict[str, Any]:
    """Items of the _bindcraft* category of a ranked mmCIF, if any, kept verbatim.

    A ranked file of v1.0.3 carries the design's identity, its campaign settings and every metric, also per
    validation model (<name>_model_1_ptm); the refolded complexes carry identity and outcome only.
    """
    import gemmi

    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except (RuntimeError, ValueError, OSError):
        return {}
    out: dict[str, Any] = {}
    for item in block:
        if item.pair is not None and item.pair[0].lower().startswith("_bindcraft"):
            out[item.pair[0]] = item.pair[1].strip("'\"")
        elif item.loop is not None:
            for i, tag in enumerate(item.loop.tags):
                if tag.lower().startswith("_bindcraft"):
                    vals = list(item.loop.values[i:: len(item.loop.tags)])
                    out[tag] = [v.strip("'\"") for v in vals]
    return out


def one_based_copy(run: RunData, design: DesignRecord, out_path: str | Path) -> Path:
    """Write a copy of the design's structure with binder residue numbers shifted to 1-based.

    Only needed when annotations['binder_numbering']['offset_to_one_based'] is 1 (zero-based source); with
    offset 0 the file is copied unchanged. Target chains keep their numbers. Returns the path written.
    """
    import gemmi

    if not design.structure_path:
        raise ValueError(f"{design.id} has no structure")
    st = gemmi.read_structure(design.structure_path)
    offset = int(run.annotations.get("binder_numbering", {}).get("offset_to_one_based", 0))
    nb = int(design.extras.get("n_binder_chains", 1))
    names = [c.name for c in st[0]]
    for name in names[len(names) - nb:]:
        for res in st[0][name]:
            res.seqid = gemmi.SeqId(res.seqid.num + offset, res.seqid.icode)
    out = Path(out_path)
    st.make_mmcif_document().write_file(str(out))
    return out


# ---------------------------------------------------------------------------------------------
# Modality
# ---------------------------------------------------------------------------------------------


def _modality_of(name_info: Optional[dict], columns: set[str], n_binder_chains: int, meta_modality: Optional[str]) -> tuple[str, str]:
    """(modality, source) for one design: the name token, then the metadata, then the columns."""
    if name_info and name_info["modality"]:
        return name_info["modality"], "design name"
    if meta_modality:
        low = meta_modality.lower().replace(" ", "").replace("_", "").replace("-", "")
        for key, val in _MODALITY_WORDS.items():
            if key.replace("-", "") in low:
                return val, "campaign_metadata"
    if "Interdomain_Contact_Fraction" in columns or n_binder_chains >= 2:
        return "scFv", "columns"
    if any(c in columns for c in ANTIBODY_COLUMNS):
        return "VHH", "columns"
    return "de novo", "default"


# ---------------------------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------------------------


def _row_metrics(row: Mapping[str, str], lower_cols: Mapping[str, str], scale_plddt: float) -> dict[str, Optional[float]]:
    out: dict[str, Optional[float]] = {}
    for low, src in _LOWER_METRIC.items():
        col = lower_cols.get(low)
        if col is None:
            continue
        v = _num(row.get(col))
        if v is not None and src in _PLDDT_SOURCES:
            v = v * scale_plddt
        for stem, factor in METRIC_COLUMNS[src]:
            out[stem] = None if v is None else round(v * factor, 6)
    return out


def _extras(row: Mapping[str, str]) -> tuple[dict[str, Any], list[str]]:
    """Columns with no field of their own: (values, names of columns the reader does not know at all)."""
    extras: dict[str, Any] = {}
    unknown: list[str] = []
    for col, val in row.items():
        low = col.strip().lower()
        if low in _FIELD_COLUMNS or low in _LOWER_METRIC:
            continue
        extras[col] = _extra_value(val)
        if low not in _KNOWN_COLUMNS:
            unknown.append(col)
    return extras, unknown


def _plddt_scale(df: pd.DataFrame, lower_cols: Mapping[str, str]) -> float:
    for src in _PLDDT_SOURCES:
        col = lower_cols.get(src.lower())
        if col is None:
            continue
        vals = [v for v in (_num(x) for x in df[col]) if v is not None]
        if vals and max(vals) > 1.5:
            return 0.01
    return 1.0


def _is_accepted_word(outcome: str) -> bool:
    low = outcome.lower()
    return bool(low) and not low.startswith(("not", "fail", "reject")) and any(w in low for w in _ACCEPT_WORDS)


def _row_hash(rec: Mapping[str, str], lower: Mapping[str, str], design: str) -> str:
    """The trajectory hash of a row: its `hash` column, else the one in the design name."""
    h = _txt(rec.get(lower.get("hash", ""), "")).lower()
    if h:
        return h
    info = parse_name(design)
    return info["hash"] if info else ""


def _backbone_name(design: str, hash_: str, campaign: str) -> str:
    """The backbone id of a design whose hash has no trajectory row: the design name without its _candidate<N> /
    _seq<N> suffix (which is the trajectory's own design name), else <campaign>_<hash>."""
    stripped = _SUFFIX_RE.sub("", design)
    return stripped if hash_ and stripped.lower().endswith(hash_) else f"{campaign}_{hash_}"


def _legacy_backbone_id(campaign: str, trajectory: str, design: str) -> str:
    """For tables without a hash column: one family per `trajectory` value (a bare number becomes
    <campaign>_traj<n>), and a row with no trajectory value is its own family."""
    if not trajectory:
        return design
    if re.fullmatch(r"\d+(\.0+)?", trajectory):
        return f"{campaign}_traj{int(float(trajectory))}"
    return trajectory


@dataclass
class _Cand:
    """One design before it becomes a DesignRecord: its refolded row and, if it was ranked, its ranked row."""

    id: str
    row: dict[str, str]
    lower: dict[str, str]
    scale: float
    stage: str
    hash: str
    ranked_row: Optional[dict[str, str]] = None
    ranked_lower: dict[str, str] = field(default_factory=dict)
    ranked_id: Optional[str] = None


def _read_state(root: Path) -> dict[str, Any]:
    p = _child(root, ".campaign_state.json")
    if p is None or not p.is_file():
        return {}
    try:
        data = json.loads(p.read_text())
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


# ---------------------------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------------------------


def read_run(
    path: str | Path,
    run_id: str,
    label: Optional[str] = None,
    settings_path: Optional[str | Path] = None,
    construct_path: Optional[str | Path] = None,
    check_structures: bool = True,
    hash_structures: bool = True,
    **options: Any,
) -> RunData:
    """Read one BindCraft 2 campaign folder into a RunData (see the module docstring).

    settings_path     a campaign metadata file to use instead of campaign_metadata*.json in the folder.
    construct_path    the target file given to BindCraft: deleted regions are reported against it and its hash is
                      compared with the recorded one. BindCraft 2 does not copy it into the campaign, but the
                      metadata records where it was (settings.targets[].target_path); that file is used when
                      this argument is not given and the path exists on this machine.
    check_structures  open every mmCIF once (chains, target numbers, binder numbering, sequence agreement).
    hash_structures   include the structure files in RunData.files with their sha256.
    Unknown keyword options are ignored.
    """
    root = Path(path).resolve()
    if not detect(root):
        raise ReaderError(f"{root} does not look like a BindCraft 2 campaign folder "
                         "(expected 1_Trajectories, 2_Refolded, 3_Ranked with !_*.csv tables)")
    label = label or root.name
    warnings: list[str] = []
    notes: list[str] = []
    files: list[InputFile] = []
    annotations: dict[str, Any] = {}

    def add_file(role: str, p: Path) -> None:
        files.append(InputFile(role=role, path=str(p.resolve()), sha256=sha256_file(p), size=p.stat().st_size))

    meta = _read_metadata(root, notes, warnings, files, Path(settings_path) if settings_path else None)

    # -- tables --------------------------------------------------------------------------------
    tables: dict[str, pd.DataFrame] = {}
    for stage in STAGE_DIRS:
        p = _table(root, stage)
        if p is None:
            notes.append(f"No table in {stage}/ (expected {TABLES[stage]}).")
            continue
        try:
            tables[stage] = read_table(p)
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            warnings.append(f"{stage}/{p.name} could not be read ({exc}); skipped.")
            continue
        add_file("csv", p)
    for extra_name in ("summary.csv", ".campaign_state.json"):
        p = _child(root, extra_name)
        if p is not None and p.is_file():
            add_file("csv" if extra_name.endswith(".csv") else "state", p)
    if REFOLDED_DIR not in tables and RANKED_DIR not in tables:
        warnings.append("Neither !_Refolded.csv nor !_Ranked.csv was read: the run has no designs.")

    all_cols: set[str] = set()
    unknown_cols: dict[str, list[str]] = {}
    for stage, df in tables.items():
        all_cols.update(df.columns)
        unk = [c for c in df.columns if c.strip().lower() not in _KNOWN_COLUMNS]
        if unk:
            unknown_cols[stage] = unk
    if unknown_cols:
        annotations["unknown_columns"] = unknown_cols
        notes.append("Columns the reader does not know are kept in DesignRecord.extras: "
                     + "; ".join(f"{s}: {', '.join(c[:6])}" for s, c in unknown_cols.items()))

    # the version: the metadata's, else the stamp BindCraft 2 writes into every table row
    row_versions: Counter[str] = Counter()
    for df in tables.values():
        col = next((c for c in df.columns if c.strip().lower() == "bindcraft_version"), None)
        if col is not None:
            row_versions.update(v for v in (_txt(x) for x in df[col]) if v)
    stamped = {_clean_version(v) for v in row_versions}
    if len(stamped) > 1:
        warnings.append(f"The table rows carry several BindCraft versions ({', '.join(sorted(row_versions))}): "
                        "the campaign was resumed with another release.")
    elif meta.get("version") and stamped and stamped != {meta["version"]}:
        warnings.append(f"The metadata says BindCraft {meta['version']} but the table rows say {', '.join(sorted(row_versions))}.")
    version = meta.get("version") or (sorted(v for v in stamped if v)[0] if stamped - {None} else None)

    # -- structure files, by stage -------------------------------------------------------------
    s_index = {s: _structure_index(_child(root, s)) for s in STAGE_DIRS}

    # -- candidates: refolded rows, with the ranked row attached to the same molecule -----------------
    plddt_scale: dict[str, float] = {}
    lowers: dict[str, dict[str, str]] = {}
    for stage in (REFOLDED_DIR, RANKED_DIR):
        df = tables.get(stage)
        if df is None:
            continue
        lowers[stage] = {c.strip().lower(): c for c in df.columns}
        plddt_scale[stage] = _plddt_scale(df, lowers[stage])
        if plddt_scale[stage] != 1.0:
            notes.append(f"{stage}: pLDDT columns look like 0-100; divided by 100 to the BindCraft 1 scale (0-1).")
        if "design" not in lowers[stage]:
            warnings.append(f"{stage} table has no 'design' column; its rows were skipped.")
            tables.pop(stage)
            lowers.pop(stage)
    traj_df, ref_df, rank_df = (tables.get(s) for s in STAGE_DIRS)

    cands: dict[str, _Cand] = {}
    order: list[str] = []
    by_key: dict[tuple[str, str], list[str]] = {}
    if ref_df is not None:
        lower = lowers[REFOLDED_DIR]
        for rec in ref_df.to_dict("records"):
            did = _txt(rec.get(lower["design"]))
            if not did or did in cands:
                continue
            h = _row_hash(rec, lower, did)
            cands[did] = _Cand(did, rec, lower, plddt_scale[REFOLDED_DIR], REFOLDED_DIR, h)
            order.append(did)
            seq = _txt(rec.get(lower.get("binder_sequence", ""), ""))
            if h and seq:
                by_key.setdefault((h, seq), []).append(did)
    ranked_only: list[str] = []
    if rank_df is not None:
        lower = lowers[RANKED_DIR]
        for rec in rank_df.to_dict("records"):
            rid = _txt(rec.get(lower["design"]))
            if not rid:
                continue
            h = _row_hash(rec, lower, rid)
            seq = _txt(rec.get(lower.get("binder_sequence", ""), ""))
            # the same molecule has another id in 3_Ranked (_seq0) than in 2_Refolded (_candidate<N>): match on
            # (hash, sequence); an equal id is also taken, for tables that keep one id
            target = rid if rid in cands and cands[rid].ranked_row is None else None
            if target is None and h and seq:
                target = next((d for d in by_key.get((h, seq), []) if cands[d].ranked_row is None), None)
            if target is None:
                if rid in cands:
                    warnings.append(f"Ranked id {rid} appears twice in {TABLES[RANKED_DIR]}; the first row was kept.")
                    continue
                cands[rid] = _Cand(rid, rec, lower, plddt_scale[RANKED_DIR], RANKED_DIR, h, rec, lower, rid)
                order.append(rid)
                ranked_only.append(rid)
            else:
                c = cands[target]
                c.ranked_row, c.ranked_lower, c.ranked_id = rec, lower, rid
    if ranked_only and ref_df is not None:
        warnings.append(_aggregate(ranked_only, f"{len(ranked_only)} ranked designs match no refolded row on (hash, "
                                   "Binder_Sequence); each is read from its ranked row alone"))

    campaign_guess = None
    for did in order:
        info = parse_name(did)
        if info:
            campaign_guess = info["campaign"]
            break
    campaign = _txt(meta.get("campaign")) or campaign_guess or root.name

    # -- structure file of each design ------------------------------------------------------------------
    chain_counts: dict[str, int] = {}
    struct_of: dict[str, Optional[Path]] = {}
    ranked_of: dict[str, Optional[Path]] = {}
    monomer_of: dict[str, Optional[Path]] = {}
    for did in order:
        c = cands[did]
        seq = _txt(c.row.get(c.lower.get("binder_sequence", ""), ""))
        chain_counts[did] = max(1, len(seq.split("/"))) if seq else 1
        ranked_of[did] = _match_structure(s_index[RANKED_DIR], c.ranked_id) if c.ranked_id else None
        # a design's structure is its refolded complex; the ranked file is the same prediction (identical atoms on
        # a real campaign) and stands in for it only when there is no refolded one. A trajectory structure is
        # another sequence and is never used.
        struct_of[did] = (_match_structure(s_index[REFOLDED_DIR], did) if c.stage == REFOLDED_DIR else None) or ranked_of[did]
        monomer_of[did] = s_index[REFOLDED_DIR].get(did + MONOMER_SUFFIX)

    # A binder of several chains (an scFv's VH and VL, a homo-oligomer's copies)
    # is read whole: n_binder_chains travels with the design and
    # structure.load_complex takes that many chains. What such a binder changes
    # about the analysis is handled where it matters rather than here -- fold
    # distance is not reported for it (families.fold_unreliable) and the
    # terminus readout says it answers for one pair of ends (termini.Format).
    inspected: dict[str, dict[str, Any]] = {}
    if check_structures:
        bad_read: list[str] = []
        for did in order:
            p = struct_of[did]
            if p is None:
                continue
            try:
                inspected[did] = inspect_structure(p, chain_counts[did])
            except Exception as exc:  # corrupt file: reported once, the design stays without a check
                bad_read.append(f"{did} ({exc})")
        if bad_read:
            warnings.append(_aggregate(bad_read, f"{len(bad_read)} structures could not be read"))

    # binder numbering origin: structure evidence first, then the version stamp
    firsts = Counter(f for info in inspected.values() for f in info["binder_first_numbers"] if f is not None)
    ver = _version_tuple(version)
    zero_based: Optional[bool]
    basis: str
    if firsts:
        top, _ = firsts.most_common(1)[0]
        zero_based = top == 0
        basis = "structure"
        if len(firsts) > 1:
            warnings.append(f"Binder chains start at different residue numbers ({dict(firsts)}); used the most common.")
        if ver is not None and (ver <= LAST_ZERO_BASED_VERSION) != zero_based:
            notes.append(f"Version stamp {version} and the structures disagree on binder numbering "
                         f"(structures start at {top}); the structures were followed (stamps can lag).")
    elif ver is not None:
        zero_based, basis = ver <= LAST_ZERO_BASED_VERSION, "version"
    else:
        zero_based, basis = False, "assumed"
        notes.append("Binder numbering origin unknown (no structure checked, no version stamp): assumed 1-based.")
    shift = 1 if zero_based else 0
    annotations["binder_numbering"] = {
        "source_zero_based": zero_based, "basis": basis, "offset_to_one_based": shift,
        "note": "Structure files keep the numbering BindCraft wrote; use one_based_copy() to renumber a file.",
    }
    if shift:
        notes.append("BindCraft 2 binder numbering is zero-based in this campaign; interface residues in "
                     "DesignRecord are shifted to 1-based, the structure files are unchanged "
                     f"(source: {basis}).")

    # -- modality --------------------------------------------------------------------------------
    mods: dict[str, str] = {}
    for did in order:
        mods[did], _src = _modality_of(parse_name(did), all_cols, chain_counts[did], meta.get("modality"))
    run_mod_counts = Counter(mods.values())
    if run_mod_counts:
        modality = run_mod_counts.most_common(1)[0][0] if len(run_mod_counts) == 1 else "mixed"
    else:
        modality, _ = _modality_of(None, all_cols, 1, meta.get("modality"))
    annotations["modality"] = modality
    annotations["modalities"] = dict(run_mod_counts)

    # -- filters from the metadata -----------------------------------------------------------------
    parsed_filters = meta.get("filters") or []
    if parsed_filters:
        annotations["filters"] = parsed_filters
    if meta.get("filters_unlimited"):
        annotations["filters_unlimited"] = list(meta["filters_unlimited"])
    if meta.get("scaffold"):
        annotations["scaffold"] = dict(meta["scaffold"])
    filter_by_name = {f["name"].lower(): f for f in parsed_filters}

    # -- backbones: one per trajectory row, keyed by hash -----------------------------------------------
    backbones: dict[str, BackboneRecord] = {}
    bb_of_hash: dict[str, str] = {}
    traj_hashes: set[str] = set()
    traj_info: dict[str, Any] = {"rows": 0, "completed": None, "terminated": None}
    if traj_df is not None:
        lower = {c.strip().lower(): c for c in traj_df.columns}
        scale = _plddt_scale(traj_df, lower)
        recs = traj_df.to_dict("records")
        # generation order: the `trajectory` number when every row has a distinct one (rows are appended as
        # trajectories finish, so the file order is not the order they were started in), else the row order
        nums = [_num(r.get(lower.get("trajectory", ""))) for r in recs]
        by_number = None not in nums and len(set(nums)) == len(nums)
        seq_order = sorted(range(len(recs)), key=lambda i: nums[i]) if by_number else list(range(len(recs)))
        for pos, i in enumerate(seq_order):
            rec = recs[i]
            g = lambda key, rec=rec, lower=lower: _txt(rec.get(lower.get(key, ""), ""))  # noqa: E731
            did = g("design")
            h = _row_hash(rec, lower, did)
            bid = (did or f"{campaign}_{h}") if h else _legacy_backbone_id(campaign, g("trajectory"), did)
            if bid in backbones:
                warnings.append(f"Trajectory {bid} appears twice in {TABLES[TRAJECTORY_DIR]}; the first row was kept.")
                continue
            seq = g("binder_sequence")
            info = parse_name(did)
            extras, _ = _extras(rec)
            extras["design"] = did
            if g("trajectory"):
                tnum = _num(g("trajectory"))
                extras["trajectory"] = int(tnum) if tnum is not None and tnum.is_integer() else _extra_value(g("trajectory"))
            if h:
                extras["hash"] = h
                bb_of_hash.setdefault(h, bid)
                traj_hashes.add(h)
            n_ch = max(1, len(seq.split("/"))) if seq else 1
            res, _chk, _wrong, _bad = interface_residues_from(
                g("interface_binder_residues"), seq, [chr(ord("B") + k) for k in range(n_ch)], shift)
            backbones[bid] = BackboneRecord(
                id=bid, run_id=run_id, binder_name=campaign,
                length=int(_num(g("length")) or (len(seq.replace("/", "")) if seq else (info["length"] if info else 0))),
                seed=0, order=pos, trajectory_metrics=_row_metrics(rec, lower, scale), start_sequence=seq or None,
                interface_residues=res, extras=extras,
            )
        traj_info["rows"] = len(backbones)
        term_col = lower.get("terminated")
        if term_col is not None:
            ended = [bool(_txt(r.get(term_col))) for r in recs]
            traj_info["terminated"] = sum(ended)
            traj_info["completed"] = len(ended) - sum(ended)

    # -- designs ---------------------------------------------------------------------------------
    designs: list[DesignRecord] = []
    checked_tokens = wrong_tokens = unparsed_tokens = 0
    wrong_designs: list[str] = []
    target_checked = target_wrong = 0
    target_wrong_designs: list[str] = []
    no_family: list[str] = []
    ranked_contradiction: list[str] = []
    outcome_words: Counter[str] = Counter()
    not_in_mean = 0
    letters_cache: dict[tuple, dict[int, str]] = {}
    has_outcome_col = REFOLDED_DIR in lowers and "outcome" in lowers[REFOLDED_DIR]
    for did in order:
        c = cands[did]
        row, lower = c.row, c.lower
        g = lambda key: _txt(row.get(lower.get(key, ""), ""))  # noqa: E731
        info = parse_name(did)
        seq_cell = g("binder_sequence")
        n_chains = chain_counts[did]
        # The sequence wins when the binder has several chains: BindCraft 2 writes
        # the `length` column as the FIRST chain's length there (131 for an scFv
        # whose own id says l238 and whose Binder_Sequence holds 131 + 107), and
        # a binder reported at half its size is wrong everywhere it is used --
        # the peptide threshold, the tag note, the sequence column.
        whole = len(seq_cell.replace("/", "")) if seq_cell else 0
        if n_chains > 1 and whole:
            length = whole
        else:
            length = int(_num(g("length")) or whole or (info["length"] if info else 0))
        metrics = _row_metrics(row, lower, c.scale)
        extras, _unk = _extras(row)
        if c.ranked_row is not None and c.stage == REFOLDED_DIR:
            # columns only the ranked row has (settings_core, ...) are added; the metrics are the refolded row's
            r_extras, _ = _extras(c.ranked_row)
            for k, v in r_extras.items():
                extras.setdefault(k, v)
        # the family: the hash, joined to its trajectory row
        if c.hash:
            backbone = bb_of_hash.get(c.hash) or _backbone_name(did, c.hash, campaign)
        else:
            traj = g("trajectory")
            backbone = _legacy_backbone_id(campaign, traj, did)
            if not traj:
                no_family.append(did)
        names = split_failed_filters(g("failed_filters"))
        extras["failed_filter_names"] = names
        extras["n_binder_chains"] = n_chains
        extras["modality"] = mods[did]
        if c.hash:
            extras["hash"] = c.hash
        if info and info["kind"] == "candidate":
            extras["candidate"] = info["n"]
        stru = inspected.get(did)
        if stru is not None:
            extras["n_target_chains"] = len(stru["target_chains"])
        outcome = _txt(extras.get("outcome"))
        if outcome:
            outcome_words[outcome] += 1
        ranked = c.ranked_row is not None
        if has_outcome_col and c.stage == REFOLDED_DIR:
            accepted = _is_accepted_word(outcome) and not names
        elif rank_df is not None:
            accepted = ranked                     # no outcome to go by: the ranked shortlist is all that is known
        else:
            accepted = _is_accepted_word(outcome) and not names
        if ranked and c.stage == REFOLDED_DIR and not accepted:
            ranked_contradiction.append(did)
        rank_val = _num(c.ranked_row.get(c.ranked_lower.get("rank", ""))) if ranked else None
        extras["ranked"] = ranked
        if ranked:
            extras["ranked_id"] = c.ranked_id
        if monomer_of[did] is not None:
            extras["unbound_binder_path"] = str(monomer_of[did])
        # Interface_Binder_Residues: letters are amino acids, numbers are positions in the written structure
        bchains = list(stru["binder_chains"]) if stru and len(stru["binder_chains"]) == n_chains else [chr(ord("B") + k) for k in range(n_chains)]
        res, chk, wrong, bad = interface_residues_from(g("interface_binder_residues"), seq_cell, bchains, shift)
        checked_tokens += chk
        wrong_tokens += wrong
        unparsed_tokens += bad
        if wrong or bad:
            wrong_designs.append(did)
            extras["Interface_Binder_Residues"] = g("interface_binder_residues")   # kept as written, not read
            res = []
        # Interface_Target_Residues is checked the same way, against the target chain of the structure
        tcell = _txt(extras.get("Interface_Target_Residues"))
        if stru is not None and tcell and "/" not in tcell and "target_numbers" in stru:
            key = (stru["target_numbers"], stru["target_seq"])
            tmap = letters_cache.setdefault(key, dict(zip(*key)))
            tparts, _tbad = parse_residue_list(tcell)
            hit = False
            for letter, num in (t for part in tparts for t in part):
                if letter:
                    target_checked += 1
                    if tmap.get(num) != letter:
                        target_wrong += 1
                        hit = True
            if hit:
                target_wrong_designs.append(did)
        rec = DesignRecord(
            id=did, run_id=run_id, backbone_id=backbone, status="accepted" if accepted else "rejected",
            binder_name=campaign, length=length, seed=0, mpnn_index=None, sequence=seq_cell,
            rank=int(rank_val) if rank_val is not None else None, protocol=mods[did],
            target_hotspot=meta.get("hotspots"), metrics=metrics, interface_residues=res,
            structure_path=str(struct_of[did]) if struct_of[did] else None,
            ranked_path=str(ranked_of[did]) if ranked_of[did] else None, extras=extras,
        )
        rec.failed_filters, unrep = _failures(names, filter_by_name, row, lower, c.scale)
        if unrep:
            extras["failed_filter_not_in_mean"] = unrep
            not_in_mean += 1
        if ranked_of[did] is not None:
            block = _read_bindcraft_block(ranked_of[did])
            if block:
                rec.extras["cif_bindcraft"] = block
        designs.append(rec)
    if no_family:
        warnings.append(_aggregate(no_family, f"{len(no_family)} designs have neither a hash nor a trajectory value, so each is "
                                   "its own family and the family counts are not meaningful"))
    if ranked_contradiction:
        warnings.append(_aggregate(ranked_contradiction, f"{len(ranked_contradiction)} designs are in {TABLES[RANKED_DIR]} but "
                                   "their refolded outcome is not passed; they are not counted as accepted"))
    odd_outcomes = sorted(o for o in outcome_words if o.lower() not in ("passed", "rejected"))
    if odd_outcomes:
        notes.append(f"outcome values other than passed / rejected: {', '.join(odd_outcomes)}; read by their wording.")
    if ref_df is not None and not has_outcome_col:
        warnings.append("!_Refolded.csv has no outcome column: only the designs of !_Ranked.csv are counted as accepted."
                        if rank_df is not None else
                        "!_Refolded.csv has no outcome column and there is no ranked table: no design can be called accepted.")
    if wrong_designs:
        warnings.append(_aggregate(
            wrong_designs, f"{len(wrong_designs)} designs have interface residue letters that do not match Binder_Sequence "
            f"({wrong_tokens} of {checked_tokens} letters checked, {unparsed_tokens} tokens unparsed): the format of "
            "Interface_Binder_Residues may have changed, so their interface residues were not read"))
    if target_wrong:
        warnings.append(_aggregate(
            target_wrong_designs, f"{len(target_wrong_designs)} designs have Interface_Target_Residues letters that do not match "
            f"the target in the structure ({target_wrong} of {target_checked} letters)"))
    differing = _ranked_vs_refolded(cands)
    if differing:
        warnings.append(_aggregate(sorted(differing), f"{len(differing)} metric columns differ between ranked rows and the refolded "
                                   "rows they were matched to (they are expected to be the same prediction)"))

    # -- structure checks ---------------------------------------------------------------------------
    target = TargetConstruct(run_id=run_id)
    construct = Path(construct_path) if construct_path else None
    if construct is None and meta.get("target_file"):
        recorded = Path(str(meta["target_file"])).expanduser()
        if recorded.is_absolute() and recorded.is_file():
            construct = recorded
            notes.append(f"Target construct taken from the metadata: {recorded}.")
        elif recorded.is_absolute():
            notes.append(f"The metadata records the target file {recorded}, which is not on this machine; "
                         "pass construct_path to compare the structures with it.")
    if check_structures:
        _check_structures(designs, inspected, target, warnings, notes, construct, root, meta)
    _check_table_structure_agreement(designs, cands, warnings, s_index)
    _note_trajectory_folders(root, meta, notes)

    # -- backbones: designs join their trajectory --------------------------------------------------------
    orphans: list[str] = []
    for d in designs:
        bb = backbones.get(d.backbone_id)
        if bb is None:
            orphans.append(d.backbone_id)
            bb = backbones[d.backbone_id] = BackboneRecord(
                id=d.backbone_id, run_id=run_id, binder_name=campaign, length=d.length, seed=0, order=-1)
        bb.design_ids.append(d.id)
        if d.status == "accepted":
            bb.accepted_ids.append(d.id)
    if orphans and traj_df is not None:
        warnings.append(_aggregate(sorted(set(orphans)), "backbones with designs but no row in !_Trajectories.csv"))
    backbone_list = sorted(backbones.values(), key=lambda b: (b.order < 0, b.order, b.id))

    # -- settings -------------------------------------------------------------------------------------
    settings = _settings(run_id, label, campaign, modality, meta, designs, parsed_filters, notes, target)
    for n in settings.hotspots:
        letter = target.residue_at(n)
        if letter is not None:
            target.hotspot_residues[n] = letter
    if target.numbers and settings.hotspots:
        absent = [n for n in settings.hotspots if n not in target.numbers]
        if absent:
            warnings.append(f"hotspot residues {absent} are not present in the target of this run")

    # -- capability map, provenance ---------------------------------------------------------------------
    have = sorted({k for d in designs if d.status == "accepted" for k, v in d.metrics.items() if v is not None})
    annotations["capabilities"] = {
        "absent_bindcraft1_stems": [s for s in ABSENT_STEMS if s in METRIC_STEMS],
        "unmapped_bindcraft1_stems": {s: {"bindcraft2": c, "why": w} for s, (c, w) in UNMAPPED_STEMS.items()},
        "present": have,
        "note": "Absent stems are not reported by BindCraft 2; checks that need them are switched off. Unmapped stems have "
                "a BindCraft 2 column of a different definition, which stays in extras under its own name.",
    }
    outcomes = Counter(_txt(d.extras.get("outcome")) for d in designs if d.extras.get("outcome") is not None)
    if outcomes:
        annotations["outcomes"] = dict(outcomes)
    annotations["tiers"] = {
        "trajectories": traj_info["rows"] if traj_df is not None else None,
        "refolded": sum(1 for c in cands.values() if c.stage == REFOLDED_DIR),
        "passed": sum(1 for d in designs if d.status == "accepted"),
        "ranked": sum(1 for d in designs if d.extras.get("ranked")),
    }
    if not_in_mean:
        annotations["designs_with_failure_not_in_mean"] = not_in_mean
    if hash_structures:
        for d in designs:
            if d.structure_path:
                add_file(f"{d.status}_structure", Path(d.structure_path))
    if version or meta.get("revision"):
        tool_version = "2" + (f" v{version}" if version else "") + (f" (revision {meta['revision']})" if meta.get("revision") else "")
    else:
        tool_version = "2 (version not recorded)"
    if meta.get("revision") is None and version:
        notes.append("Only a version stamp was found; stamps can lag the code, the revision is unknown.")

    funnel = _funnel(designs, backbone_list, traj_df, settings, traj_info, _read_state(root), traj_hashes, notes)

    # -- how far this campaign is covered by what was verified ------------------------------------------
    reasons = _untested_reasons(version, modality, meta, tables, lowers)
    if reasons:
        warnings.insert(0, untested_warning(reasons))
    return RunData(
        run_id=run_id, label=label, reader=NAME, source_dir=str(root), settings=settings, target=target,
        designs=designs, backbones=backbone_list, funnel=funnel, failure_counters=_failure_counters(designs),
        files=files, warnings=warnings, notes=notes, annotations=annotations, tool_version=tool_version,
        experimental=bool(reasons),
    )


def _untested_reasons(version: Optional[str], modality: str, meta: Mapping[str, Any], tables: Mapping[str, pd.DataFrame],
                      lowers: Mapping[str, Mapping[str, str]]) -> list[str]:
    """Why this campaign is outside what the reader was verified on (empty: verified layout, version and modality).

    The decision needs the tables and the metadata only, so it is the same with and without check_structures."""
    reasons: list[str] = []
    ver = _version_tuple(version)
    if ver is None:
        reasons.append("its BindCraft 2 version is not recorded")
    elif ver not in VERIFIED_VERSIONS:
        reasons.append(f"it was written by BindCraft 2 v{version}")
    if modality not in VERIFIED_MODALITIES:
        reasons.append(f"its modality is {modality}")
    settings_mod = _find_key(meta.get("raw") or {}, ("settings_modality",), (str,)) or meta.get("modality")
    if modality in VERIFIED_MODALITIES and settings_mod is not None and str(settings_mod).lower() != VERIFIED_SETTINGS_MODALITY:
        reasons.append(f"its designs are named de novo but its settings name the modality '{settings_mod}'")   # else the modality above says it
    if REFOLDED_DIR in tables:
        missing = [c for c in _REQUIRED_REFOLDED if c.lower() not in lowers.get(REFOLDED_DIR, {})]
        if missing:
            reasons.append(f"!_Refolded.csv lacks the column(s) {', '.join(missing)}")
    if TRAJECTORY_DIR in tables and "hash" not in {c.strip().lower() for c in tables[TRAJECTORY_DIR].columns}:
        reasons.append("!_Trajectories.csv has no hash column")
    if not meta.get("raw"):
        reasons.append("there is no readable campaign_metadata.json")
    return reasons


def untested_warning(reasons: list[str]) -> str:
    """Why this campaign is outside what the reader was checked on, and what that does and does not mean.

    It used to list every modality nobody had checked -- "VHH, scFv, Fab,
    cyclic-peptide and multi-target campaigns" -- on every campaign, whichever
    one this was, and then ask the reader to check the numbers against the CSVs
    by hand. The review does that itself now, so the warning says what the
    automatic check covers and what it cannot.
    """
    verified = ", ".join(f"v{'.'.join(map(str, v))}" for v in VERIFIED_VERSIONS)
    return (f"UNTESTED campaign: this reader was checked number by number on de novo (binder modality) campaigns "
            f"written by BindCraft 2 {verified}; here {'; '.join(reasons)}. A campaign outside that is read on the "
            "same assumptions, so a column name, a value format or a metadata key that BindCraft writes differently "
            "for it would be read wrongly. The arithmetic is checked for you: every count in this review is compared "
            "with the campaign's own tables (manifest validation.counts_match_csv, and the build log says so in "
            "words). That catches a table read wrongly; it cannot catch a column that means something else here.")


def _ranked_vs_refolded(cands: Mapping[str, _Cand]) -> set[str]:
    """Metric columns whose value differs between a ranked row and the refolded row it was matched to."""
    differ: set[str] = set()
    for c in cands.values():
        if c.ranked_row is None or c.stage != REFOLDED_DIR:
            continue
        for low, src in _LOWER_METRIC.items():
            a, b = c.lower.get(low), c.ranked_lower.get(low)
            if a is None or b is None:
                continue
            x, y = _num(c.row.get(a)), _num(c.ranked_row.get(b))
            if (x is None) != (y is None) or (x is not None and y is not None and abs(x - y) > 1e-9):
                differ.add(src)
    return differ


def _failures(names: list[str], filters: Mapping[str, dict], row: Mapping[str, str], lower: Mapping[str, str],
              plddt_scale: float) -> tuple[list[FilterFailure], list[str]]:
    """One FilterFailure for every name BindCraft put in failed_filters whose threshold is known and whose value the
    row carries, and the names whose miss the row does not show.

    The verdict is BindCraft's, not the CSV's. Confidence filters (pLDDT, pTM, i_pTM, Unbound_Binder_pLDDT,
    Target_pLDDT) are judged on each validation model in turn, while the CSV holds the mean of the models run,
    rounded to two decimals; a design can fail on one model with a CSV value at or above the threshold (this occurs
    for pTM). Such a failure is kept with scope "Model" and gap 0.0, which is a lower bound; a miss that the
    mean shows has scope "Average" and its gap. The stem is the campaign's filter name, so Unbound_Binder_pLDDT and
    pLDDT stay two filters.
    """
    out: list[FilterFailure] = []
    unreproduced: list[str] = []
    for n in names:
        base = re.sub(r"\s*\(not measured\)\s*$", "", n).partition(".")[0]
        f = filters.get(base.lower())
        col = lower.get(base.lower())
        if not f or col is None:
            continue
        value = _num(row.get(col))
        if value is None:
            continue
        if base in _PLDDT_SOURCES:
            value *= plddt_scale
        higher = f["op"] == ">="
        gap = (f["threshold"] - value) if higher else (value - f["threshold"])
        if gap > 0:
            out.append(FilterFailure(f["name"], "Average", value, f["threshold"], higher, gap))
        else:
            out.append(FilterFailure(f["name"], "Model", value, f["threshold"], higher, 0.0))
            unreproduced.append(n)
    return out, unreproduced


def _failure_counters(designs: list[DesignRecord]) -> dict[str, int]:
    """Not-accepted designs per named filter, from failed_filters (a design counts once per filter it names)."""
    c: Counter[str] = Counter()
    for d in designs:
        if d.status == "rejected":
            c.update(set(d.extras.get("failed_filter_names", [])))
    return dict(c)


def _check_structures(
    designs: list[DesignRecord], inspected: Mapping[str, dict[str, Any]], target: TargetConstruct,
    warnings: list[str], notes: list[str], construct: Optional[Path], root: Path, meta: Mapping[str, Any],
) -> None:
    sets: Counter[tuple[int, ...]] = Counter()
    seqs: dict[tuple[int, ...], str] = {}
    set_of: dict[str, tuple[int, ...]] = {}
    bad_len, bad_seq, h_files, multi = [], [], [], 0
    by_id = {d.id: d for d in designs}
    for did, info in inspected.items():
        d = by_id[did]
        if "target_numbers" not in info or len(info["chains"]) < 2:
            warnings.append(f"{did}: structure has chains {info['chains']}; expected target chains then binder chains.")
            continue
        nums = tuple(info["target_numbers"])
        sets[nums] += 1
        seqs.setdefault(nums, info["target_seq"])
        set_of[did] = nums
        if len(info["target_chains"]) > 1:
            multi += 1
        if info["hydrogens"]:
            h_files.append(did)
        if sum(len(s) for s in info["binder_seqs"]) != d.length:
            bad_len.append(did)
        if d.sequence and "/".join(info["binder_seqs"]) != d.sequence:
            bad_seq.append(did)
    target.n_structures_checked = len(inspected)
    if bad_len:
        warnings.append(_aggregate(bad_len, f"{len(bad_len)} structures whose binder chains do not add up to the CSV length"))
    if bad_seq:
        warnings.append(_aggregate(bad_seq, f"{len(bad_seq)} structures whose binder sequence differs from the CSV Binder_Sequence"))
    if h_files:
        warnings.append(_aggregate(h_files, f"{len(h_files)} structures contain hydrogens (BindCraft 2 documents heavy atoms only)"))
    if multi:
        notes.append(f"{multi} structures have several target chains; the target construct describes the first one (chain A).")
    if not sets:
        return
    majority, _ = sets.most_common(1)[0]
    target.numbers = list(majority)
    target.sequence = seqs[majority]
    target.first, target.last = majority[0], majority[-1]
    target.consistent = len(sets) == 1
    if not target.consistent:
        odd = [i for i, s in set_of.items() if s != majority]
        warnings.append(_aggregate(odd, f"{len(odd)} structures have a different target residue set than the majority ({len(majority)} residues)"))
    ref = _reference_construct(construct) if construct is not None and construct.is_file() else None
    if ref is not None and construct is not None:
        ref_numbers, ref_seq = ref
        target.reference_path = str(construct.resolve())
        target.reference_length = len(ref_numbers)
        by_num = dict(zip(ref_numbers, ref_seq))
        mism = [n for n, aa in zip(target.numbers, target.sequence) if by_num.get(n) != aa]
        target.sequence_matches_reference = not mism
        if mism:
            warnings.append(f"{len(mism)} target residues differ from the target file {construct.name} (first at {mism[0]})")
        target.deleted_ranges = _missing_ranges([x for x in ref_numbers if x not in set(target.numbers)])
        rec = _txt(meta.get("target_hash")).lower()
        if rec:
            import hashlib

            data = construct.read_bytes()
            hit = [a for a in ("sha256", "sha1", "md5") if hashlib.new(a, data).hexdigest().startswith(rec[:32])]
            notes.append(f"Target file {construct.name} " + (f"matches the hash recorded in the metadata ({hit[0]})." if hit
                         else "does not match the hash recorded in the metadata (or the hash uses another algorithm)."))
    else:
        target.deleted_ranges = _gaps(target.numbers)
        if meta.get("target_hash"):
            notes.append("The target file is not copied into the campaign folder; its recorded hash is in "
                         "settings.settings_raw. Pass construct_path to compare.")


def _check_table_structure_agreement(designs, cands, warnings, s_index) -> None:
    acc_missing = [d.id for d in designs if d.status == "accepted" and not d.structure_path]
    if acc_missing:
        warnings.append(_aggregate(acc_missing, f"{len(acc_missing)} accepted designs have a CSV row but no mmCIF"))
    rank_missing = [c.ranked_id for c in cands.values() if c.ranked_id and not _match_structure(s_index[RANKED_DIR], c.ranked_id)]
    if rank_missing:
        warnings.append(_aggregate(rank_missing, f"{len(rank_missing)} ranked designs have a CSV row but no file in {RANKED_DIR}"))
    known = {c.ranked_id for c in cands.values() if c.ranked_id}
    orphan = sorted(s for s in s_index[RANKED_DIR] if not any(_stem_matches(s, k) for k in known | set(cands)))
    if orphan:
        warnings.append(_aggregate(orphan, f"{len(orphan)} mmCIFs in {RANKED_DIR} have no CSV row"))


def _note_trajectory_folders(root: Path, meta: Mapping[str, Any], notes: list[str]) -> None:
    """Trajectory structures are written only when archive_trajectories is true; without it a trajectory folder holds
    <name>_losses.csv alone. A missing trajectory structure is therefore normal and is not a warning."""
    d = _child(root, TRAJECTORY_DIR)
    if d is None or not d.is_dir():
        return
    folders = [p for p in d.iterdir() if p.is_dir()]
    if not folders:
        return
    with_structure = sum(1 for p in folders if any(q.suffix.lower() in STRUCTURE_SUFFIXES for q in p.iterdir()))
    if with_structure == 0:
        why = " (archive_trajectories is false in the settings)" if meta.get("archive_trajectories") is False else ""
        notes.append(f"{len(folders)} trajectory folders hold no structures{why}; trajectory structures are not read.")
    else:
        notes.append(f"{with_structure} of {len(folders)} trajectory folders hold structures; they are not read "
                     "(a design's structure is its refolded complex).")


def _settings(run_id, label, campaign, modality, meta, designs, parsed_filters, notes, target) -> RunSettings:
    st = RunSettings(run_id=run_id, label=label, binder_name=campaign, protocol=modality)
    raw = meta.get("raw") or {}
    st.settings_raw = raw
    st.settings_path = meta.get("path")
    st.notes.append(f"Modality: {modality}.")
    if meta.get("hotspots"):
        st.hotspot_string = meta["hotspots"]
        st.hotspots = parse_hotspots(meta["hotspots"])
        st.hotspot_source = "campaign_metadata"
        if re.search(r"[A-Za-z]", meta["hotspots"]) and len(re.findall(r"[A-Za-z]+", meta["hotspots"])) > 1:
            st.notes.append("Hotspot chain letters were dropped; hotspots are residue numbers of the first target chain.")
    if meta.get("chains"):
        c = meta["chains"]
        st.chains = ",".join(c) if isinstance(c, list) else str(c)
    if meta.get("length_range"):
        st.length_min, st.length_max = meta["length_range"]
        st.length_source = "campaign_metadata"
    elif designs:
        ls = [d.length for d in designs if d.length]
        if ls:
            st.length_min, st.length_max, st.length_source = min(ls), max(ls), "csv"
    if meta.get("final_designs") is not None:
        st.final_designs_target = meta["final_designs"]
        st.final_designs_target_source = "campaign_metadata"
    if meta.get("target_file"):
        st.starting_pdb_recorded = str(meta["target_file"])
        st.starting_pdb_name = Path(str(meta["target_file"])).name
    if meta.get("target_hash"):
        st.notes.append(f"Target file hash recorded in the metadata: {meta['target_hash']} (the file itself is not copied).")
    if parsed_filters:
        st.filters_name = "campaign_metadata"
        n_binding = sum(1 for f in parsed_filters if f["mandatory"])
        st.notes.append(f"{len(parsed_filters)} filters with a finite threshold were read from the metadata ({n_binding} "
                        f"mandatory); {len(meta.get('filters_unlimited') or [])} with an infinite threshold (no limit) were dropped.")
    notes.extend(st.notes)
    return st


def _funnel(designs, backbones, traj_df, settings, traj_info, state, traj_hashes, notes) -> FunnelCounts:
    f = FunnelCounts()
    st, nt = f.status, f.notes
    rows, completed, terminated = traj_info["rows"], traj_info["completed"], traj_info["terminated"]
    if traj_df is None:
        st["trajectory_rows"] = "unknown"
        nt["trajectory_rows"] = "No !_Trajectories.csv."
    elif completed is not None:
        # the table lists every trajectory that was started, finished or not (`terminated` names the stage it stopped in)
        f.trajectories_attempted = rows
        st["trajectories_attempted"] = "observed"
        nt["trajectories_attempted"] = ("Rows in !_Trajectories.csv: every trajectory that was started is listed, with the "
                                        "stage it stopped in (`terminated`).")
        f.trajectories_terminated = terminated
        st["trajectories_terminated"] = "observed"
        nt["trajectories_terminated"] = "Rows whose `terminated` column names a stage; they carry no metrics."
        f.trajectory_rows = completed
        st["trajectory_rows"] = "observed"
        nt["trajectory_rows"] = "Rows with an empty `terminated`: the trajectory passed every design stage and was refolded."
        listed = state.get("attempted")
        if isinstance(listed, list) and listed:
            if {str(x) for x in listed} == traj_hashes:
                notes.append(f".campaign_state.json lists the same {len(listed)} trajectory hashes as the table.")
            else:
                notes.append(f".campaign_state.json lists {len(listed)} attempted trajectories and the table {rows} "
                             "(a campaign still running, or a table edited by hand).")
    else:
        f.trajectory_rows = rows
        st["trajectory_rows"] = "observed"
        nt["trajectory_rows"] = "Rows in !_Trajectories.csv; whether unfinished trajectories are listed there is not recorded."
        f.trajectories_attempted_estimate = rows
        st["trajectories_attempted_estimate"] = "derived"
        nt["trajectories_attempted_estimate"] = "Derived: assumes every attempted trajectory has a row in !_Trajectories.csv."
    if f.trajectories_attempted is None:
        st["trajectories_attempted"] = "unknown"
    f.mpnn_rows = len(designs)
    st["mpnn_rows"] = "observed"
    nt["mpnn_rows"] = "Rows in !_Refolded.csv (refolded candidates; a ranked row is not counted again)."
    f.accepted_designs = sum(1 for d in designs if d.status == "accepted")
    f.mpnn_rows_rejected = len(designs) - f.accepted_designs
    f.backbones_with_mpnn_rows = sum(1 for b in backbones if b.design_ids)
    f.accepted_backbones = sum(1 for b in backbones if b.accepted_ids)
    for k in ("mpnn_rows_rejected", "backbones_with_mpnn_rows", "accepted_designs", "accepted_backbones"):
        st[k] = "observed"
    nt["accepted_designs"] = "Refolded rows whose outcome is passed (the ranked shortlist is a subset of them)."
    nt["accepted_backbones"] = "A family is one trajectory; counts of families, never of sequences."
    f.final_designs_target = settings.final_designs_target
    if f.final_designs_target is not None:
        st["final_designs_target"] = "observed"
    return f
