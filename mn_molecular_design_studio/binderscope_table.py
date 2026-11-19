"""Small, source-bound BinderScope summaries for the saved-design table."""
from __future__ import annotations

import json
import math


PREFIX = "Analysis · BinderScope · "
VIEWER_COLUMN = "BinderScope viewer"
DEFAULT_FIELDS = ("Status", "Result · binderscope_buried_area_angstrom2", "Result · binderscope_hydrogen_bond_count",
    "Result · binderscope_salt_bridge_count", "Result · binderscope_buried_unpaired_polar_atom_count", "Hotspots")
OPTIONAL_FIELDS = ("Assessed model", "Analysis version", "Mean binder pLDDT", "N tag accessibility", "N terminal pLDDT",
    "N terminal confidence", "C tag accessibility", "C terminal pLDDT", "C terminal confidence", "Target-reference RMSD (Å)")
NUMBER_FORMATS = {PREFIX + field: "%.1f" for field in ("Result · binderscope_buried_area_angstrom2",
    "Mean binder pLDDT", "N terminal pLDDT", "C terminal pLDDT")}
NUMBER_FORMATS[PREFIX + "Target-reference RMSD (Å)"] = "%.3f"
SUMMARY_PATHS = {
    "digest": "source.input_structure_digest", "model": "source.source.structure_label",
    "encoding": "source.assessment.b_factor_interpretation", "plddt": "analysis.plddt_res", "residues": "analysis.n",
    "hotspot_status": "hotspot_compliance.status", "requested": "hotspot_compliance.requested_count",
    "contacted": "hotspot_compliance.contacted_count", "unmapped": "hotspot_compliance.unmapped_count",
    "reference_status": "target_reference_comparison.status", "rmsd": "target_reference_comparison.target_ca_rmsd_angstrom",
    "termini_status": "termini.status", "termini_applies": "termini.applies", "n": "termini.n", "c": "termini.c",
}


def _finite(value):
    return type(value) in {float, int} and math.isfinite(value)


def _tags(terminal):
    tags = terminal.get("g") or {}
    names = {"his6": "His6", "strep": "Twin-Strep", "domain": "GFP/SNAP/mCherry", "large": "Fc/MBP"}
    if not isinstance(tags, dict) or not tags or terminal.get("conf") == "off" or not all(isinstance(value, dict) for value in tags.values()):
        return "Unavailable"
    if set(tags) == set(names) and all(tags[key].get("s") == "free" for key in names):
        return "Free · all 4 families"
    return "; ".join(names.get(key, key) + ": " + str(value.get("s") or "unknown") +
        (f" ({value['l']} aa)" if value.get("s") == "linker" and value.get("l") is not None else "")
        for key, value in sorted(tags.items()))


def summary_columns(summary):
    row = {PREFIX + "Assessed model": summary.get("model") or "Saved complex"}
    values = summary.get("plddt") or []
    if (summary.get("encoding") == "plddt_0_100" and isinstance(values, list) and values
            and len(values) == summary.get("residues") and all(_finite(value) and 0 <= value <= 100 for value in values)):
        row[PREFIX + "Mean binder pLDDT"] = sum(values) / len(values)
    requested, contacted, unmapped = (summary.get(key) for key in ("requested", "contacted", "unmapped"))
    if (summary.get("hotspot_status") in {"assessed", "partial"} and type(requested) is int and requested > 0
            and type(contacted) is int and 0 <= contacted <= requested and type(unmapped) is int and 0 <= unmapped <= requested - contacted):
        row[PREFIX + "Hotspots"] = f"{contacted} / {requested}" + (f" · {unmapped} unmapped" if unmapped else "")
    else:
        row[PREFIX + "Hotspots"] = "Not specified" if summary.get("hotspot_status") == "not_specified" else "Not assessed"
    if summary.get("reference_status") == "assessed" and _finite(summary.get("rmsd")):
        row[PREFIX + "Target-reference RMSD (Å)"] = summary["rmsd"]
    for end in ("n", "c"):
        terminal = summary.get(end) or {}
        if not isinstance(terminal, dict):
            terminal = {}
        prefix = PREFIX + end.upper() + " "
        if summary.get("termini_applies") is False or summary.get("termini_applies") == 0:
            row[prefix + "tag accessibility"] = "Not applicable"
        elif summary.get("termini_status") == "unavailable" or summary.get("encoding") != "plddt_0_100":
            row[prefix + "tag accessibility"] = "Unavailable"
        else:
            row[prefix + "tag accessibility"] = _tags(terminal)
            if terminal.get("conf") in {"low", "ok"}:
                row[prefix + "terminal confidence"] = {"low": "Low", "ok": "OK"}[terminal["conf"]]
                if summary.get("encoding") == "plddt_0_100" and _finite(terminal.get("plddt")) and 0 <= terminal["plddt"] <= 100:
                    row[prefix + "terminal pLDDT"] = terminal["plddt"]
    return row


def inspector_choices(row):
    try:
        records = json.loads(row.get(PREFIX + "Records") or "[]")
    except (ValueError, TypeError):
        return []
    if not isinstance(records, list):
        return []
    return [record for record in records if isinstance(record, dict) and record.get("observation_id")
        and record.get("evaluation_id") and record.get("relationship") == "Selected structure"
        and str(record.get("status") or "").casefold() in {"completed", "complete", "success", "succeeded", "ok"}
        and record.get("report_available")]


def inspector_label(row):
    choices = inspector_choices(row)
    if len(choices) == 1:
        return "View BinderScope"
    return f"Choose assessment ({len(choices)})" if choices else None
