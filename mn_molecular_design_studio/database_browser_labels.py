"""Human-facing labels for saved SQL evidence, without changing its identity."""
from __future__ import annotations

import re


def tool_label(value):
    key = str(value or "").casefold().replace("-", "_")
    if key in {"binderscope", "analysis.binderscope"}:
        return "BinderScope"
    aliases = {"udp": "Uni-Dock Pro", "unidock_pro": "Uni-Dock Pro", "vina": "Vina", "gnina": "GNINA",
        "boltz_2": "Boltz-2", "boltz2": "Boltz-2", "optimized_prediction_boltz2": "Boltz-2",
        "alphafold_3": "AlphaFold 3", "optimized_prediction_alphafold3": "AlphaFold 3",
        "nesso_1": "Nesso-1", "nesso_affinity": "Nesso-1", "plip": "PLIP", "md.analysis": "MD Analysis",
        "md.production": "MD production", "md.mmgbsa.gb": "MM/GBSA", "md.mmgbsa.pb": "MM/PBSA"}
    if key in aliases:
        return aliases[key]
    for marker, label in (("biotite", "Biotite"), ("ipsae", "ipSAE"), ("pep_patch", "PEP-Patch"),
                          ("peppatch", "PEP-Patch"), ("rfd4", "RFD4-Proteina"), ("bindcraft", "BindCraft2"),
                          ("boltzgen", "BoltzGen"), ("ab_mpnn", "AbMPNN"), ("rfantibody", "RFAntibody"),
                          ("validation_acceptance", "Filter decision")):
        if marker in key:
            return label
    return str(value or "Unknown tool").replace("_", " ")


def readable_name(value):
    text = str(value or "").strip()
    return bool(text and not re.fullmatch(r"(?:(?:mol|target|structure|evaluation)-)?[a-f0-9-]{24,}", text, re.I)
                and text.casefold() not in {"target", "ligand", "two_chain", "complex", "idr_sequence_carrier", "sequence_input"})


def descriptive_target_name(value):
    text = str(value or "").strip()
    return bool(readable_name(text) and text.casefold() not in {
        "protein", "receptor", "binder", "molecule", "sequence-defined target", "recorded target", "unnamed target",
    } and not re.fullmatch(r"(?:unnamed\s+)?(?:protein|dna|rna)(?:\s+target)?\s*[·(]?\s*\d+\s*(?:aa|nt)\)?", text, re.I))


def run_label(value):
    text = str(value or "")
    if re.fullmatch(r"\d{8}-\d{6}-[A-Za-z0-9_-]+", text):
        return f"{text[:4]}-{text[4:6]}-{text[6:8]} {text[9:11]}:{text[11:13]}:{text[13:15]}"
    return text or "Run not recorded"


def calculation_label(result):
    metrics = result.get("metrics") or {}
    replicate = metrics.get("replicate")
    if replicate is None:
        matched = re.search(r"legacy-replicate-(\d+)$", result.get("candidate") or "")
        replicate = int(matched[1]) if matched else None
    model = result.get("model")
    parts = [tool_label(result.get("tool"))]
    if model is not None and str(model) != "":
        parts.append("Model / pose " + str(model))
    elif replicate is not None:
        parts.append("Replicate " + (format(replicate, "g") if isinstance(replicate, float) else str(replicate)))
    attempt = metrics.get("attempt_index")
    if attempt is not None:
        parts.append("Attempt " + (format(attempt, "g") if isinstance(attempt, float) else str(attempt)))
    run = result.get("run_id") or ""
    if run and not run.startswith("mnligand-"):
        parts.append(run_label(run))
    return " · ".join(parts)


def target_label(index, target_id):
    molecule = index["molecules"].get(target_id)
    if molecule is None:
        return "Target name not recorded" if target_id else "Target not linked"
    from mn_molecular_design_studio.database_browser import _component_identity, _key
    from mn_molecular_design_studio.entity_names import target_names
    if len(molecule.get("components", [])) == 1:
        identity = _component_identity(molecule["components"][0])
        labels = target_names(index, _key(["component", identity]), molecule.get("group_ids", set())) if identity else []
        if labels:
            return " / ".join(sorted({label["name"] for label in labels}))
    names = molecule.get("target_names", set()) | molecule.get("import_names", set())
    names = [name for name in names if readable_name(name)]
    if names:
        return min(names, key=lambda value: (value.casefold(), value))
    if readable_name(molecule.get("name")):
        return molecule["name"]
    components = molecule.get("components") or []
    proteins = [component for component in components if component.get("type") == "protein"]
    if proteins:
        return f"Unnamed protein target ({sum(len(component.get('sequence') or '') for component in proteins)} aa)"
    return "Unnamed ligand" if any(component.get("type") == "small_molecule" for component in components) else "Unnamed target"


def target_construct_label(index, target_id, *, reference_label=""):
    name = reference_label or target_label(index, target_id)
    components = index["molecules"].get(target_id, {}).get("components", [])
    proteins = [component for component in components if component.get("type") == "protein"]
    length = sum(len(component.get("sequence") or "") for component in proteins)
    if length and not re.search(r"\b" + str(length) + r"\s*aa\b", name, re.I):
        name += " · " + str(length) + " aa"
    if len(proteins) > 1:
        name += " · " + str(len(proteins)) + " protein chains"
    return name


def analysis_inventory(results):
    inventory = {}
    for row in results:
        label = tool_label(row["tool"])
        if label == "Filter decision":
            continue
        bucket = inventory.setdefault(label, {"same_structure": 0, "other_structure": 0, "failed": 0, "unknown": 0})
        status = str(row.get("status") or "").casefold()
        if status in {"failed", "error", "unavailable", "not_available", "cancelled", "canceled"}:
            bucket["failed"] += 1
        elif status not in {"completed", "complete", "success", "succeeded", "ok"}:
            bucket["unknown"] += 1
        elif row.get("structure_relation") in {"Other structure of this design", "Structure link unverified"}:
            bucket["other_structure"] += 1
        else:
            bucket["same_structure"] += 1
    return inventory


def analysis_status(inventory, label):
    bucket = inventory.get(label, {})
    if bucket.get("same_structure"):
        return "Recorded" + (f" ({bucket['same_structure']})" if bucket["same_structure"] > 1 else "")
    if bucket.get("other_structure"):
        return "Other model"
    if bucket.get("failed"):
        return "Failed / unavailable"
    if bucket.get("unknown"):
        return "Status unknown"
    return "—"
