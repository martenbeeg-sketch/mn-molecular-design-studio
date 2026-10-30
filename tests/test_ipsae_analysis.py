import csv
import inspect
from types import SimpleNamespace

from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import (
    WORKFLOWS,
    _get_function,
    workflow_schema,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _parse_ipsae_row
from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
from mn_molecular_design_studio.validation_acceptance import METRICS, apply_acceptance, available_metrics


def test_ipdae_cutoff_is_configurable_in_api_and_graph_node():
    function = _get_function(WORKFLOWS["analysis.ipsae"])
    assert inspect.signature(function).parameters["ipdae_contact_cutoff"].default == 8.0
    assert inspect.signature(function).parameters["structure_contact_cutoff"].default == 4.0
    assert inspect.signature(function).parameters["hotspot_contact_cutoff"].default == 5.0
    schema = workflow_schema("analysis.ipsae")
    cutoff = next(row for row in schema["parameters"] if row["name"] == "ipdae_contact_cutoff")
    assert cutoff["default"] == 8.0
    assert {row["name"] for row in schema["parameters"]} >= {"structure_contact_cutoff", "hotspot_contact_cutoff"}
    assert NODE_CATALOG["ipsae"]["label"] == "ipSAE & Structure Analysis"
    node_cutoff = next(row for row in NODE_CATALOG["ipsae"]["fields"] if row["name"] == "ipdae_contact_cutoff")
    assert node_cutoff["default"] == 8.0


def test_ipdae_container_row_fields_are_parsed_without_losing_error_status(tmp_path):
    path = tmp_path / "row.csv"
    headers = ["ID", "ipsae", "ipsae_d0chn", "ipsae_d0dom", "iptm_af", "iptm_d0chn",
               "binder_chain", "target_chain", "type", "engine", "pae_format", "error",
               "ipae_binder_to_target", "ipae_target_to_binder", "ipae_bidirectional",
               "ipae_binder_to_target_norm31", "ipae_target_to_binder_norm31",
               "ipae_bidirectional_norm31", "ipae_contact_pair_count", "ipae_contact_cutoff",
               "bindcraft2_ipdae_recomputed", "bindcraft2_ipdae_binder_anchor_max",
               "bindcraft2_ipdae_target_anchor_max", "bindcraft2_ipdae_contact_pair_count",
               "bindcraft2_ipdae_contact_cutoff", "pae_measurement_error",
               "structure_interface_atom_contact_pair_count", "structure_interface_residue_contact_pair_count",
               "structure_binder_interface_residue_count", "structure_target_interface_residue_count",
               "structure_ca_clash_pair_count_2_5a", "structure_binder_ca_bfactor_mean",
               "structure_target_ca_bfactor_mean", "structure_complex_ca_bfactor_mean", "structure_contact_cutoff",
               "structure_hotspot_count", "structure_hotspot_found_count", "structure_hotspot_contacted_count",
               "structure_hotspot_contact_fraction_found", "structure_hotspot_contact_fraction_requested",
               "structure_hotspot_min_distance", "structure_hotspot_mean_min_distance", "structure_hotspot_contact_cutoff",
               "structure_hotspots_contacted", "structure_hotspots_missing", "structure_metrics_status",
               "structure_metrics_error", "design_vs_cofold_rmsd_status", "design_vs_cofold_rmsd_error",
               "structure_binder_design_to_cofold_fitted_ca_rmsd", "structure_binder_rmsd_ca_count",
               "structure_binder_rmsd_alignment_coverage", "hotspot_analysis_status", "hotspot_analysis_error"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writerow({"ID": "candidate", "ipsae": "0.72", "error": "", "ipae_binder_to_target": "7.5",
                         "ipae_target_to_binder": "8.5", "ipae_bidirectional": "8.0",
                         "bindcraft2_ipdae_recomputed": "0.74", "bindcraft2_ipdae_binder_anchor_max": "0.74",
                         "bindcraft2_ipdae_target_anchor_max": "0.61", "pae_measurement_error": ""})
    result = _parse_ipsae_row(path)
    assert result["ipsae"] == 0.72
    assert result["ipae_binder_to_target"] == 7.5
    assert result["ipae_target_to_binder"] == 8.5
    assert result["ipae_bidirectional"] == 8.0
    assert result["bindcraft2_ipdae_recomputed"] == 0.74

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writerow({"ID": "candidate", "structure_interface_atom_contact_pair_count": "25",
                         "structure_binder_ca_bfactor_mean": "84.976", "structure_metrics_status": "measured",
                         "structure_hotspots_contacted": "A25", "structure_hotspot_contact_fraction_requested": "1.0",
                         "design_vs_cofold_rmsd_status": "unavailable",
                         "design_vs_cofold_rmsd_error": "no separate designed/reference binder structure was supplied"})
    structure = _parse_ipsae_row(path)
    assert structure["structure_interface_atom_contact_pair_count"] == 25
    assert structure["structure_binder_ca_bfactor_mean"] == 84.976
    assert structure["structure_metrics_status"] == "measured"
    assert structure["structure_hotspots_contacted"] == "A25"
    assert structure["design_vs_cofold_rmsd_status"] == "unavailable"

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writerow({"ID": "candidate", "error": "PAE scorer failed", "pae_measurement_error": "PAE axes unavailable"})
    failed = _parse_ipsae_row(path)
    assert failed["ipsae_error"] == "PAE scorer failed"
    assert failed["pae_measurement_error"] == "PAE axes unavailable"


def test_structure_only_candidate_is_reported_as_partial_measurement(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis

    finished = {}
    monkeypatch.setattr(analysis, "write_candidates", lambda *_args, **_kwargs: {"candidates_jsonl": "candidates.jsonl"})
    monkeypatch.setattr(analysis, "finish_job", lambda _run_dir, success, payload: finished.update(success=success, payload=payload))
    candidate = {"candidate_id": "candidate", "binder_chains": ["B"], "target_chains": ["A"],
                 "metrics": {"ipsae_error": "PAE file missing", "structure_metrics_status": "measured",
                             "structure_interface_atom_contact_pair_count": 12}}
    analysis._finish_metric_analysis(
        SimpleNamespace(run_dir=tmp_path), [candidate], tool="ipsae", table_name="ipsae_metrics.csv",
        metric_names=["structure_interface_atom_contact_pair_count"],
    )
    assert finished["success"] is True
    assert finished["payload"]["scientific_status"] == "partial"
    assert finished["payload"]["metrics"]["structure_measured_count"] == 1
    assert candidate["metrics"]["analysis_status"] == "partial"
    assert candidate["raw_metadata"]["measurement_status"]["status"] == "partial"


def test_design_contact_preservation_is_selectable_and_missing_data_never_passes():
    metric_field = METRICS["ipsae_design_contact_preservation_fraction"]["field"]
    available = dict(available_metrics("miniprotein", "fasta_idr", "complex", engine="ipsae"))
    other_engine = dict(available_metrics("miniprotein", "fasta_idr", "complex", engine="boltz2"))
    assert available["ipsae_design_contact_preservation_fraction"]["field"] == metric_field
    assert "ipsae_design_contact_preservation_fraction" not in other_engine

    policy = {
        "binder_format": "miniprotein",
        "target_type": "fasta_idr",
        "mandatory_filters": [{"field": metric_field, "operator": ">=", "value": 0.5}],
    }
    rows, summary = apply_acceptance([
        {"candidate_id": "preserved", "source_tool": "ipsae", "stage": "analysis",
         "metrics": {"structure_design_interaction_preservation_fraction": 0.75}},
        {"candidate_id": "changed", "source_tool": "ipsae", "stage": "analysis",
         "metrics": {"structure_design_interaction_preservation_fraction": 0.25}},
        {"candidate_id": "no-design-complex", "source_tool": "ipsae", "stage": "analysis", "metrics": {}},
    ], policy)
    status = {row["candidate_id"]: row["raw_metadata"]["validation_acceptance"]["status"] for row in rows}
    assert status == {"preserved": "accepted", "changed": "rejected", "no-design-complex": "unresolved"}
    assert summary["accepted"] == summary["rejected"] == summary["unresolved"] == 1
