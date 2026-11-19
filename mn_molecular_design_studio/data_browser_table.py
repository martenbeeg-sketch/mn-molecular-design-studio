"""Read-only, evaluation-scoped table data and complete folder CSV exports."""
from __future__ import annotations

import csv
import io
import json

from mn_molecular_design_studio.data_browser import admission_measurements, concise_result_label
from mn_molecular_design_studio.database_browser_labels import run_label, target_construct_label, tool_label
from mn_molecular_design_studio.binderscope_table import DEFAULT_FIELDS, OPTIONAL_FIELDS, PREFIX


BASE_COLUMNS = ["Name", "Project", "Target", "Target construct", "Target sequence", "Method", "Activity", "Replicate", "Attempt", "Model / pose",
                "MD engine", "Starting pose method", "Production", "Replicas", "Production segments", "Saved analyses",
                "Structures", "Evaluations", "Admission", "Admission basis", "Post-analysis", "Status", "Production status", "Molecule type", "Sequence", "SMILES", "Target PDB",
                "Created", "Workstation", "Collection", "Run", "Candidate",
                "Engine version", "Result ID", "Molecule ID", "Evaluation ID", "Observation ID", "Structure ID", "Tool ID"]
RESULT_PRIORITY = ("best_score_kcal_mol", "cnn_score", "cnn_affinity", "ranking_score", "confidence_score",
                   "iptm", "i_pTM", "oracle_ppi.iptm", "plddt_mean", "pLDDT", "oracle_ppi.plddt",
                   "affinity_pred_value", "affinity_probability_binary", "Binder_RMSD")
INPUT_PRIORITY = ("seed", "exhaustiveness", "binder_length", "num_recycles", "engine_options.sampling_steps",
                  "samples", "temperature", "production_ns")
ROW_FIELDS = ("replicate", "attempt_index", "seed", "pose_index", "pose_selection_criterion")
COLUMN_PRESETS = ("Automatic", "Design", "Docking", "Structure prediction", "Binding affinity", "Molecular dynamics", "Compact", "All columns")
PROFILES = {
    "Design": {"activities": {"Protein design", "Small-molecule design", "Nucleic-acid design"},
        "results": ("oracle_ppi.iptm_coi", "oracle_ppi.iptm", "i_pTM", "iptm", "oracle_ppi.plddt_coi", "oracle_ppi.plddt", "pLDDT", "plddt_mean", "oracle_ppi.ipae_coi", "i_pAE", "Binder_RMSD", "refolded_backbone_rmsd.rmsd"),
        "inputs": ("binder_length", "folding_engine", "seed", "temperature")},
    "Docking": {"activities": {"Docking"}, "results": ("best_score_kcal_mol", "cnn_score", "cnn_affinity"),
        "inputs": ("seed", "exhaustiveness", "pose_selection_criterion")},
    "Structure prediction": {"activities": {"Structure prediction", "Cofolding"},
        "results": ("iptm", "ranking_score", "confidence_score", "ptm", "plddt_mean", "complex_plddt"),
        "inputs": ("num_recycles", "samples", "seed", "engine_options.sampling_steps")},
    "Binding affinity": {"activities": {"Binding affinity"},
        "results": ("affinity_pred_value", "affinity_probability_binary", "predicted_affinity", "pKd", "kd_nm"),
        "inputs": ("seed", "temperature")},
    "Molecular dynamics": {"activities": {"Molecular dynamics"},
        "results": ("production_actual_duration_ns", "backbone_rmsd_mean_angstrom", "ligand_rmsd_mean_angstrom",
                    "backbone_rmsd_tail_mean_angstrom", "ligand_rmsd_tail_mean_angstrom", "radius_of_gyration_mean_angstrom"),
        "inputs": ("production_ns", "temperature_kelvin", "temperature", "timestep_fs", "seed")},
}
ANALYSIS_PRIORITY = {
    "ipSAE": ("ipsae_min", "ipSAE_min", "ipsae", "ipSAE", "ipsae_binder_to_target"),
    "Biotite": ("interface_residue_contact_pairs", "interface_contact_pairs", "binder_interface_residues", "binder_ca_radius_of_gyration"),
    "PEP-Patch": ("peppatch_total_area", "peppatch_hydrophobic_area"),
    "PLIP": ("hydrogen_bond_count", "hydrophobic_contact_count", "interaction_count"),
    "MD Analysis": ("backbone_rmsd_tail_mean_angstrom", "backbone_rmsd_mean_angstrom", "ligand_rmsd_mean_angstrom"),
}


def _value(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list)) else value


def _first(*values):
    return next((value for value in values if value is not None and value != ""), None)


def _parameters(values, prefix=""):
    for key, value in sorted(values.items()):
        name = prefix + key
        if isinstance(value, dict) and value:
            yield from _parameters(value, name + ".")
        else:
            yield name, _value(value)


def evaluation_entries(entries):
    """Expand only originals already in the selected folder, without adding analyses."""
    expanded, seen = [], set()
    for entry in entries:
        for original in entry["outcome"]["originals"]:
            if not original.get("evaluation_id"):
                continue
            key = (entry["outcome"]["id"], original["id"])
            if key in seen:
                continue
            seen.add(key)
            expanded.append({**entry, "result": original, "structure_id": "", "evaluation_row": True})
    return expanded


def accessible_entries(index, entries):
    return [entry for entry in entries if entry["result"]["group_id"] in index["groups"] and
            entry["result"]["molecule_id"] in index["molecules"] and entry["result"]["group_id"] in
            index["molecules"][entry["result"]["molecule_id"]]["group_ids"]]


def result_table(database, index, entries, *, cache=None):
    permitted = accessible_entries(index, entries)
    requests = sorted({(entry["result"]["observation_id"], entry["result"]["molecule_id"], entry["result"]["group_id"])
                       for entry in permitted if entry["result"].get("observation_id")})
    members = sorted({(entry["result"]["molecule_id"], entry["result"]["group_id"])
                      for entry in permitted if not entry["result"].get("observation_id")})
    evidence, admission = {}, {}
    with database.connection() as connection:
        connection.execute("BEGIN")
        if requests:
            projection = ",".join("'$.row." + field + "'" for field in ROW_FIELDS)
            for saved in connection.execute(
                "SELECT o.id,o.metrics,json_extract(o.context," + projection + ") AS row_fields "
                "FROM json_each(?) requested JOIN molecule_observations o "
                "ON o.id=json_extract(requested.value,'$[0]') AND o.molecule_id=json_extract(requested.value,'$[1]') "
                "AND o.group_id=json_extract(requested.value,'$[2]')", (json.dumps(requests),)
            ):
                evidence[saved["id"]] = (json.loads(saved["metrics"]), dict(zip(ROW_FIELDS, json.loads(saved["row_fields"]))))
        if members:
            for saved in connection.execute(
                "SELECT mm.molecule_id,mm.group_id,json_extract(mm.metadata,'$.validation_acceptance',"
                "'$.source_candidate_id','$.native_candidate_id','$.candidate_id','$.studio_naming') AS fields "
                "FROM json_each(?) requested JOIN molecule_members mm "
                "ON mm.molecule_id=json_extract(requested.value,'$[0]') AND mm.group_id=json_extract(requested.value,'$[1]')",
                (json.dumps(members),)
            ):
                values = json.loads(saved["fields"])
                metadata = dict(zip(("validation_acceptance", "source_candidate_id", "native_candidate_id", "candidate_id", "studio_naming"), values))
                admission[(saved["molecule_id"], saved["group_id"])] = admission_measurements({"metadata": metadata})
    rows = []
    for entry in permitted:
        result, outcome = entry["result"], entry["outcome"]
        evaluation = result.get("evaluation") or {}
        metrics, source = evidence.get(result.get("observation_id"), ({}, {}))
        parameters = dict(_parameters(evaluation.get("native_parameters") or {}))
        seed = _first(source.get("seed"), metrics.get("seed"))
        if seed is not None:
            if parameters.get("seed") is not None and parameters["seed"] != seed:
                parameters["requested_seed"] = parameters["seed"]
            parameters["seed"] = seed
        for field in ("pose_index", "pose_selection_criterion"):
            if source.get(field) is not None:
                parameters[field] = source[field]
        group = index["groups"][result["group_id"]]
        project = group.get("project_id") or ""
        molecule = index["molecules"][result["molecule_id"]]
        components = molecule.get("components") or []
        target_id = outcome.get("target_id") or ""
        target_components = index["molecules"].get(target_id, {}).get("components", [])
        structure_ids = [entry["structure_id"]] if entry.get("structure_id") else result.get("structure_ids", [])
        row = {"Name": concise_result_label(outcome),
            "Project": index.get("project_names", {}).get(project, "Project name not recorded" if project else "Personal"),
            "Target": outcome.get("target_reference_label") or outcome.get("target_name") or None,
            "Target construct": target_construct_label(index, target_id, reference_label=outcome.get("target_reference_label") or "")
                if any(component.get("type") == "protein" for component in target_components) else None,
            "Target sequence": "; ".join(str(component.get("chain") or "?") + ": " + component["sequence"]
                for component in target_components if component.get("sequence")) or None,
            "Method": tool_label(result["tool"]) if result.get("tool") else ", ".join(sorted({tool_label(tool) for tool in outcome["origin_tools"]})),
            "Activity": outcome["activity"], "Replicate": _first(source.get("replicate"), metrics.get("replicate"), (result.get("md") or {}).get("replica_index")),
            "Attempt": _first(source.get("attempt_index"), metrics.get("attempt_index")),
            "Model / pose": result.get("model") or None, "Structures": len(structure_ids),
            "Evaluations": len({original["evaluation_id"] for original in outcome["originals"] if original.get("evaluation_id")}),
            "Admission": outcome["eligibility"],
            "Admission basis": result.get("acceptance_level") or (result.get("receipt") or {}).get("acceptance_level") or None,
            "Status": result.get("status") or None,
            "Production status": (result.get("md") or {}).get("production_status") or (result.get("md") or {}).get("source_production_status") or None,
            "Molecule type": molecule.get("type"),
            "Sequence": "; ".join(str(component.get("chain") or "?") + ": " + component["sequence"]
                for component in components if component.get("sequence")) or None,
            "SMILES": "; ".join(component["smiles"] for component in components if component.get("smiles")) or None,
            "Target PDB": outcome.get("target_pdb_id") or result.get("target_pdb_id") or None,
            "Created": outcome.get("design_created_at") or result.get("created_at") or None,
            "Workstation": evaluation.get("source_workstation") or result.get("workstation") or outcome.get("generating_workstation") or None,
            "Collection": group["name"], "Run": run_label(result.get("run_id")),
            "Candidate": result.get("candidate") or outcome.get("native_candidate_id") or None,
            "Engine version": evaluation.get("engine_version") or None,
            "Result ID": outcome["id"], "Molecule ID": result["molecule_id"],
            "Target ID": outcome.get("target_id") or None,
            "Evaluation ID": result.get("evaluation_id") or None, "Observation ID": result.get("observation_id") or None,
            "Structure ID": "; ".join(structure_ids), "Tool ID": result.get("tool") or None}
        row.update({"Result · " + key: _value(value) for key, value in metrics.items() if key not in ROW_FIELDS})
        row.update({"Input · " + key: _value(value) for key, value in parameters.items()})
        if outcome["outcome_kind"] == "Dynamics":
            from mn_molecular_design_studio.md_browser import study_summary
            row.update(study_summary(outcome["originals"]))
            if not entry.get("evaluation_row"):
                row["Method"] = row["MD engine"]
        for measurement in admission.get((result["molecule_id"], result["group_id"]), []):
            row["Admission evidence · " + measurement["Tool"] + " · " + measurement["Measurement"]] = measurement["Value"]
        rows.append(row)
    from mn_molecular_design_studio.data_browser_analyses import linked_analysis_columns
    for row, analyses in zip(rows, linked_analysis_columns(database, index, permitted, cache=cache)):
        row.update(analyses)
    return rows


def table_columns(rows):
    available = {key for row in rows for key in row}
    return [key for key in BASE_COLUMNS if key in available] + sorted(available - set(BASE_COLUMNS))


def table_profile(rows, preset="Automatic"):
    if preset != "Automatic":
        return preset
    activities = {row.get("Activity") for row in rows}
    return next((name for name, profile in PROFILES.items() if activities and activities <= profile["activities"]), "Overview")


def default_columns(rows, preset="Automatic", *, binderscope_details=False):
    if preset == "All columns":
        return table_columns(rows)
    available = {key for row in rows for key, value in row.items() if value is not None and value != ""}
    if preset == "Compact":
        return [field for field in ("Name", "Method", "Activity", "Admission", "Post-analysis") if field in available]
    name = table_profile(rows, preset)
    profile = PROFILES.get(name, {"results": RESULT_PRIORITY, "inputs": INPUT_PRIORITY})
    results = ["Result · " + field for field in profile["results"] if "Result · " + field in available]
    if not results:
        results = sorted(field for field in available if field.startswith("Result · ") and any(
            isinstance(row.get(field), (int, float)) and not isinstance(row.get(field), bool) for row in rows))
        results = results or sorted(field for field in available if field.startswith("Admission evidence · "))
    inputs = ["Input · " + field for field in profile["inputs"] if "Input · " + field in available]
    analysis_fields = []
    for status in sorted(field for field in available if field.startswith("Analysis · ") and field.endswith(" · Status")):
        prefix = status.removesuffix("Status")
        tool = prefix.split(" · ")[1]
        if prefix == PREFIX:
            fields = (*DEFAULT_FIELDS, *(OPTIONAL_FIELDS if binderscope_details else ()))
            analysis_fields.extend(prefix + field for field in fields if prefix + field in available)
            continue
        metrics = [prefix + "Result · " + field for field in ANALYSIS_PRIORITY.get(tool, ()) if prefix + "Result · " + field in available]
        if not metrics:
            metrics = sorted(field for field in available if field.startswith(prefix + "Result · ") and any(
                isinstance(row.get(field), (int, float)) and not isinstance(row.get(field), bool) for row in rows))
        analysis_fields.extend([status, *metrics[:1]])
    md_columns = [field for field in ("Production", "Replicas", "Production segments", "Saved analyses", "Starting pose method") if field in available]
    return ["Name", "Method", *(md_columns if name == "Molecular dynamics" else []),
            *(["Admission", "Admission basis"] if name == "Design" else []), *results[:4],
            *analysis_fields, *(["Post-analysis"] if not analysis_fields and "Post-analysis" in available and
                any(row.get("Post-analysis") != "Not applicable" for row in rows) else []),
            "Project", *(["Target"] if "Target" in available else []),
            *(["Target construct"] if len({row.get("Target ID") for row in rows if row.get("Target construct")}) > 1 else []),
            *(["Production status"] if name == "Molecular dynamics" and "Production status" in available else []),
            *(field for field in ("Replicate", "Attempt", "Model / pose") if field in available),
            *inputs[:2]]


def column_label(field):
    if field.startswith("Analysis · "):
        parts = field.split(" · ")
        return parts[1] + ": " + (column_label("Result · " + " · ".join(parts[3:])) if parts[2] == "Result" else " · ".join(parts[2:]))
    if field.startswith("Result · "):
        native = field.removeprefix("Result · ")
        labels = {"best_score_kcal_mol": "Docking score (kcal/mol)", "cnn_score": "CNN pose score",
                  "cnn_affinity": "CNN affinity (native)", "iptm": "ipTM (native)", "i_pTM": "i_pTM (native)",
                  "plddt_mean": "Mean pLDDT (native)", "pLDDT": "pLDDT (native)",
                  "production_actual_duration_ns": "Selected production segment (ns)",
                  "backbone_rmsd_mean_angstrom": "Mean backbone RMSD (Å)",
                  "ligand_rmsd_mean_angstrom": "Mean ligand RMSD (Å)",
                  "radius_of_gyration_mean_angstrom": "Mean radius of gyration (Å)",
                  "binderscope_buried_area_angstrom2": "Buried interface area (Å²)",
                  "binderscope_hydrogen_bond_count": "H-bonds",
                  "binderscope_salt_bridge_count": "Salt bridges",
                  "binderscope_buried_unpaired_polar_atom_count": "Buried unpaired polar atoms"}
        return labels.get(native, native.replace("_", " "))
    if field.startswith("Input · "):
        return "Input: " + field.removeprefix("Input · ").replace("_", " ")
    return field


def table_csv(rows):
    """Keep native numbers, blank missing cells and neutralize spreadsheet formulas."""
    output = io.StringIO(newline="")
    columns = table_columns(rows)
    writer = csv.writer(output)

    def safe(value):
        if isinstance(value, str) and (value.lstrip().startswith(("=", "+", "-", "@")) or value.startswith(("\t", "\r", "\n"))):
            return "'" + value
        return value

    writer.writerow([safe(column) for column in columns])
    writer.writerows([safe(row.get(column)) for column in columns] for row in rows)
    return output.getvalue().encode("utf-8-sig")
