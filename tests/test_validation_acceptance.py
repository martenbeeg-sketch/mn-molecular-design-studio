"""Format-aware independent validation and acceptance contracts."""
import pytest

from mn_molecular_design_studio import design_validation
from mn_molecular_design_studio.validation_acceptance import (
    METRICS,
    apply_acceptance,
    available_metrics,
    bindcraft2_reference_policy,
    evaluate_candidate,
    infer_binder_format,
    infer_target_type,
)


def test_source_campaign_provenance_selects_binder_format_and_target_kind():
    candidate = {"source_tool": "bindcraft2", "binder_sequence": "A" * 120, "raw_metadata": {}}
    parameters = {"modality": "VHH", "targets": [{"target_path": "inputs/target.fasta"}]}
    assert infer_binder_format(candidate, parameters) == "vhh"
    assert infer_target_type(candidate, parameters) == "fasta_idr"


def test_format_is_not_guessed_from_sequence_length():
    candidate = {"source_tool": "unknown_generator", "binder_sequence": "A" * 15, "raw_metadata": {}}
    assert infer_binder_format(candidate, {}) is None


def test_candidate_format_takes_precedence_over_a_multi_modality_campaign_setting():
    candidate = {"source_tool": "bindcraft2", "raw_metadata": {"modality": "VHH"}}
    campaign = {"modality": ["binder", "VHH", "peptide"]}
    assert infer_binder_format(candidate, campaign) == "vhh"
    assert infer_binder_format({"source_tool": "bindcraft2"}, {"modality": "binder"}) == "miniprotein"
    assert infer_binder_format({"source_tool": "bindcraft2"}, campaign) is None


def test_bindcraft2_miniprotein_and_vhh_presets_use_compatible_filters_only():
    metrics = dict(available_metrics("vhh", "structured", "complex", engine="boltz2"))
    source_metrics = dict(available_metrics("vhh", "structured", "complex", engine="boltz2",
                                           source_tools={"bindcraft2"}, informational=True))
    mixed_source_metrics = dict(available_metrics("vhh", "structured", "complex", engine="boltz2",
                                                  source_tools={"bindcraft2", "pxdesign"}, informational=True))
    rules, unavailable = bindcraft2_reference_policy("vhh", "structured", "complex", metrics)
    by_id = {row["metric_id"]: row for row in rules}
    assert {"ptm", "iptm", "interface_residues_4a", "backbone_clashes_2_5a"} <= set(by_id)
    assert by_id["ptm"]["value"] == 0.55
    assert by_id["iptm"]["value"] == 0.70
    assert by_id["interface_residues_4a"]["value"] == 7
    assert by_id["backbone_clashes_2_5a"]["value"] == 0
    assert any("Scaffold_Sequence_Retained_Fraction" in item for item in unavailable)
    assert any("Unbound_Binder_pLDDT" in item for item in unavailable)
    assert "bindcraft2_scaffold_sequence_retained" in source_metrics
    assert "bindcraft2_scaffold_sequence_retained" not in mixed_source_metrics


def test_idr_preset_does_not_gate_on_unconstrained_target_confidence_or_pose_rmsd():
    metrics = dict(available_metrics("peptide", "fasta_idr", "complex", engine="boltz2"))
    rules, unavailable = bindcraft2_reference_policy("peptide", "fasta_idr", "complex", metrics)
    by_id = {row["metric_id"]: row for row in rules}
    assert "ptm" not in by_id
    assert by_id["iptm"]["value"] == 0.70
    assert "target_plddt" not in by_id
    assert any("sequence-defined IDR targets" in item for item in unavailable)
    assert any("Interface PAE normalized by 31 Å ≤0.33" in item for item in unavailable)
    assert any("Unbound_Binder_pLDDT" in item for item in unavailable)


def test_rfd4_idr_policy_uses_its_native_openfold3_metrics_and_saved_geometry():
    available = dict([*available_metrics(
        "miniprotein", "fasta_idr", "complex", engine="rfd4_proteina",
        source_tools={"rfd4_proteina"},
    ), *available_metrics(
        "miniprotein", "fasta_idr", "complex", engine="rfd4_proteina",
        source_tools={"rfd4_proteina"}, informational=True,
    )])
    rules, unavailable = bindcraft2_reference_policy(
        "miniprotein", "fasta_idr", "complex", available,
    )
    by_id = {row["metric_id"]: row for row in rules}
    assert set(by_id) == {
        "rfd4_native_ptm", "rfd4_native_iptm",
        "interface_pae_norm31",
        "rfd4_interface_residues_4a", "rfd4_backbone_clashes_2_5a",
    }
    assert by_id["rfd4_native_ptm"]["value"] == 0.55
    assert by_id["rfd4_native_iptm"]["value"] == 0.70
    assert by_id["interface_pae_norm31"]["value"] == 0.35
    assert by_id["interface_pae_norm31"]["operator"] == "<="
    assert not any("normalized by 31 Å" in item for item in unavailable)


def test_rfd4_native_oracle_aliases_are_normalized_without_losing_raw_scores():
    from mn_molecular_design_studio.output_vocabulary import PREFIX, candidate_output

    candidate = {
        "candidate_id": "native-rfd4-1", "source_tool": "rfd4_proteina",
        "stage": "generation.backbone_sequence", "binder_sequence": "AAAA",
        "binder_chains": ["B"], "target_chains": ["A"],
        "metrics": {"oracle_ppi.overall_ptm": 73.0, "oracle_ppi.overall_iptm": 78.0,
                    "oracle_ppi.overall_plddt": 81.65, "oracle_ppi.ipae_coi": 17.48},
    }
    output = candidate_output(candidate)
    assert output["metrics"]["confidence.generator.complex.ptm"] == 0.73
    assert output["metrics"]["confidence.generator.complex.iptm"] == 0.78
    assert output["metrics"]["confidence.generator.complex.plddt"] == 81.65
    assert output["metrics"]["confidence.generator.interface.ipae_coi"] == 17.48
    assert output["metrics"]["confidence.generator.interface.pae_norm31"] == pytest.approx(17.48 / 31)
    pae_provenance = output["provenance"]["confidence.generator.interface.pae_norm31"]
    assert pae_provenance["unit"] == "PAE / 31 Å"
    assert pae_provenance["aliases"][0]["raw_value"] == 17.48
    assert pae_provenance["aliases"][0]["conversion"] == "divide_31"
    assert PREFIX == "common."


def test_bindcraft_and_rfd4_interface_pae_share_31a_normalization_but_keep_provenance():
    policy = {
        "binder_format": "miniprotein", "target_type": "fasta_idr",
        "mandatory_filters": [{"field": METRICS["interface_pae_norm31"]["field"],
                               "operator": "<=", "value": 0.35}],
    }
    bindcraft, _ = apply_acceptance([{
        "candidate_id": "bc2", "source_tool": "bindcraft2", "stage": "generation.backbone_sequence",
        "metrics": {"i_pAE": 0.35}, "raw_metadata": {},
    }], policy)
    rfd4, _ = apply_acceptance([{
        "candidate_id": "rfd4", "source_tool": "rfd4_proteina", "stage": "generation.backbone_sequence",
        "metrics": {"oracle_ppi.ipae_coi": 10.85}, "raw_metadata": {},
    }], policy)
    bc2_measurement = bindcraft[0]["raw_metadata"]["validation_acceptance"]["measurements"][METRICS["interface_pae_norm31"]["field"]]
    rfd4_measurement = rfd4[0]["raw_metadata"]["validation_acceptance"]["measurements"][METRICS["interface_pae_norm31"]["field"]]
    assert bc2_measurement["value"] == pytest.approx(0.35)
    assert rfd4_measurement["value"] == pytest.approx(0.35)
    assert bc2_measurement["source_tool"] == "bindcraft2"
    assert rfd4_measurement["source_tool"] == "rfd4_proteina"
    assert "binder-to-target directional mean" in bc2_measurement["method"]
    assert "averaged across both directions" in rfd4_measurement["method"]
    assert bindcraft[0]["raw_metadata"]["validation_acceptance"]["status"] == "accepted"
    assert rfd4[0]["raw_metadata"]["validation_acceptance"]["status"] == "accepted"


def test_saved_evidence_acceptance_preserves_rfd4_ids_and_never_predicts(tmp_path, monkeypatch, design_source_factory):
    import json
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis

    source, candidate = design_source_factory()
    candidates_path = source / "artifacts/normalized_candidates/candidates.jsonl"
    rows = []
    for identifier, iptm in (("rfd4-pass", 78.0), ("rfd4-reject", 33.0)):
        rows.append({
            **candidate,
            "candidate_id": identifier,
            "source_tool": "rfd4_proteina",
            "stage": "generation.backbone_sequence",
            "metrics": {"oracle_ppi.overall_ptm": 73.0,
                        "oracle_ppi.overall_iptm": iptm,
                        "oracle_ppi.overall_plddt": 81.65,
                        "oracle_ppi.ipae_coi": 17.48,
                        "ppi.pass": 0.0},
            "raw_metadata": {"modality": "miniprotein", "target_type": "fasta_idr"},
        })
    candidates_path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    run_root = tmp_path / "analysis-runs"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(run_root))

    from mn_molecular_design_studio.validation_acceptance import METRICS
    policy = {
        "binder_format": "miniprotein", "target_type": "fasta_idr",
        "mandatory_filters": [
            {"field": METRICS["rfd4_native_ptm"]["field"], "operator": ">=", "value": 0.55},
            {"field": METRICS["rfd4_native_iptm"]["field"], "operator": ">=", "value": 0.70},
        ],
        "informational_scores": [METRICS["rfd4_native_plddt"]["field"],
                                 METRICS["rfd4_native_interface_pae"]["field"]],
    }
    output_dir = analysis.run_validation_acceptance_analysis(
        source, candidates_path, policy,
    )
    result = read_json(output_dir / "result.json")
    filtered = read_candidates(output_dir)
    by_id = {row["candidate_id"]: row for row in filtered}
    summary = result["metrics"]["validation_acceptance"]
    assert summary["accepted"] == 1 and summary["rejected"] == 1
    assert summary["prediction_performed"] is False
    assert set(by_id) == {"rfd4-pass", "rfd4-reject"}
    assert all(row["source_tool"] == "rfd4_proteina" for row in filtered)
    assert by_id["rfd4-pass"]["raw_metadata"]["validation_acceptance"]["status"] == "accepted"
    assert by_id["rfd4-reject"]["raw_metadata"]["validation_acceptance"]["status"] == "rejected"
    assert by_id["rfd4-pass"]["raw_metadata"]["validation_acceptance"]["native_filter_status"] == "unavailable · native PPI metrics incomplete"
    assert "confidence.generator.interface.ipae_coi" in by_id["rfd4-pass"]["normalized_output"]["metrics"]
    assert by_id["rfd4-pass"]["metrics"]["ppi.pass"] == 0.0  # native gate remains distinct
    assert (output_dir / "artifacts/validation/acceptance.csv").is_file()
    assert (output_dir / "artifacts/analysis/acceptance_measurements.csv").is_file()


def test_ipae_and_ipsae_are_not_offered_as_interchangeable_filters():
    boltz = dict(available_metrics("miniprotein", "structured", "complex", engine="boltz2"))
    af3 = dict(available_metrics("miniprotein", "structured", "complex", engine="alphafold3"))
    ipsae_analysis = dict(available_metrics("miniprotein", "structured", "complex", engine="ipsae"))
    proteina_complexa = dict(available_metrics("miniprotein", "structured", "complex", engine="proteina_complexa"))
    assert "ipsae_min" in boltz
    assert "ipsae_min" in ipsae_analysis
    assert "ipsae_bindcraft2_ipdae" in ipsae_analysis
    assert "ipsae_min" in proteina_complexa
    assert "ipsae_min" not in af3
    rules, unavailable = bindcraft2_reference_policy("miniprotein", "structured", "complex", boltz)
    assert all(row["metric_id"] != "ipsae_min" for row in rules)
    assert any("Interface PAE normalized by 31 Å" in item for item in unavailable)
    assert "predicted_pae_bidirectional_norm31" in boltz


def test_ipdae_recomputed_by_ipsae_is_filterable_without_being_called_native_bc2_score():
    ipdae = dict(available_metrics("miniprotein", "fasta_idr", "complex", engine="ipsae"))
    native = dict(available_metrics("miniprotein", "fasta_idr", "complex", engine="bindcraft2",
                                    source_tools={"bindcraft2"}, informational=True))
    assert "ipsae_bindcraft2_ipdae" in ipdae
    assert "bindcraft2_ipdae" in native
    assert ipdae["ipsae_bindcraft2_ipdae"]["field"] != native["bindcraft2_ipdae"]["field"]


def test_interface_measurements_use_bindcraft2_contact_and_ca_clash_definitions(design_source_factory):
    run, _ = design_source_factory()
    observed = design_validation.interface_filter_measurements(run / "designed.cif", ["B"], ["A"])
    assert observed == {
        "validation_interface_residue_count_4a": 1,
        "validation_backbone_ca_clashes_2_5a": 0,
    }


def test_configured_acceptance_separates_independent_metrics_from_source_ipdae():
    candidate = {
        "candidate_id": "independent-1",
        "source_tool": "optimized_prediction_boltz2",
        "stage": "complex_refolding",
        "metrics": {
            "ptm": 0.61,
            "iptm": 0.78,
            "validation_interface_residue_count_4a": 9,
            "validation_backbone_ca_clashes_2_5a": 0,
        },
        "raw_metadata": {
            "source_candidate": {
            "candidate_id": "native-1", "source_tool": "bindcraft2", "stage": "complex_refolding",
                "metrics": {"i_pDAE": 0.91, "i_pAE": 0.33,
                            "Scaffold_Sequence_Retained_Fraction": 1.0}, "raw_metadata": {},
            },
        },
        "parents": ["native-1"],
    }
    fields = {key: METRICS[key]["field"] for key in ("ptm", "iptm", "interface_residues_4a", "backbone_clashes_2_5a", "bindcraft2_ipdae", "bindcraft2_ipae", "bindcraft2_scaffold_sequence_retained")}
    candidate["metrics"]["ipae"] = 4.5
    policy = {
        "binder_format": "vhh",
        "target_type": "structured",
        "mandatory_filters": [
            {"field": fields["ptm"], "operator": ">=", "value": 0.55},
            {"field": fields["iptm"], "operator": ">=", "value": 0.7},
            {"field": fields["interface_residues_4a"], "operator": ">=", "value": 7},
            {"field": fields["backbone_clashes_2_5a"], "operator": "<=", "value": 0},
        ],
        "informational_scores": [fields["bindcraft2_ipdae"], fields["bindcraft2_ipae"], fields["bindcraft2_scaffold_sequence_retained"]],
        "ranking_metric": fields["bindcraft2_ipdae"],
    }
    rows, summary = apply_acceptance([candidate], policy)
    decision = rows[0]["raw_metadata"]["validation_acceptance"]
    assert decision["status"] == "accepted"
    assert decision["ranking_value"] == 0.91
    assert decision["measurements"][fields["bindcraft2_ipdae"]]["source_tool"] == "bindcraft2"
    assert decision["measurements"][fields["bindcraft2_ipae"]]["value"] == 0.33
    assert decision["measurements"][fields["bindcraft2_ipae"]]["unit"] == "0–1"
    assert decision["measurements"][fields["bindcraft2_scaffold_sequence_retained"]]["value"] == 1.0
    assert decision["measurements"][fields["iptm"]]["source_tool"] == "optimized_prediction_boltz2"
    assert summary["accepted"] == 1


def test_acceptance_export_keeps_source_identity_metrics_and_decision(tmp_path):
    from mn_molecular_design_studio.validation_acceptance import write_report

    candidate = {
        "candidate_id": "refolded-vhh", "parents": ["native-vhh-7"],
        "source_tool": "optimized_prediction_boltz2", "stage": "complex_refolding",
        "metrics": {"iptm": 0.4}, "raw_metadata": {},
    }
    accepted, summary = apply_acceptance([candidate], {
        "binder_format": "vhh", "target_type": "structured",
        "mandatory_filters": [{"field": METRICS["iptm"]["field"], "operator": ">=", "value": 0.7}],
    })
    report = write_report(tmp_path, accepted, summary)
    text = report.read_text()
    assert "native-vhh-7" in text
    assert "rejected" in text
    assert "vhh" in text
    assert "= 0.7" in text


def test_prediction_api_preflight_persists_selected_format_and_target(monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction

    captured = {}

    def prepare(**parameters):
        captured.update(parameters)

    monkeypatch.setattr(optimized_prediction, "_prepare_prediction", prepare)
    optimized_prediction.validate_prediction_request({
        "engine": "boltz2", "gpu_device": "0", "validation_mode": "complex",
        "validation_sources": [{"source_run_dir": "/source/run", "candidate_ids": ["vhh-1"]}],
        "validation_target_type": "fasta_idr",
        "acceptance_policy": {"binder_format": "vhh", "target_type": "fasta_idr",
                              "mandatory_filters": [], "informational_scores": []},
    })
    assert captured["validation_binder_format"] == "vhh"
    assert captured["validation_target_type"] == "fasta_idr"
    assert captured["validation_require_passing"] is False


def test_prediction_api_preflight_normalizes_remote_candidate_count(monkeypatch, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction

    captured = {}

    def prepare(**parameters):
        captured.update(parameters)

    source_run = tmp_path / "design-run"
    candidate_index = source_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    monkeypatch.setattr(optimized_prediction, "_prepare_prediction", prepare)
    optimized_prediction.validate_prediction_request({
        "engine": "protenix_v2", "gpu_device": "0", "validation_mode": "complex",
        "validation_target_type": "fasta_idr", "validation_candidate_limit": 10,
        "source_run_dir": source_run, "candidates_jsonl": candidate_index,
        "acceptance_policy": {"binder_format": "miniprotein", "target_type": "fasta_idr",
                              "mandatory_filters": [], "informational_scores": []},
    })

    assert captured["validation_sources"] == [{
        "source_run_dir": str(source_run),
        "candidates_jsonl": str(candidate_index),
        "candidate_limit": 10,
    }]
    assert "validation_candidate_limit" not in captured
    assert "source_run_dir" not in captured
    assert "candidates_jsonl" not in captured
    assert captured["validation_require_passing"] is False


def test_prediction_api_preflight_normalizes_remote_candidate_ids(monkeypatch, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction

    captured = {}
    source_run = tmp_path / "design-run"
    candidate_index = source_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    monkeypatch.setattr(optimized_prediction, "_prepare_prediction", lambda **kwargs: captured.update(kwargs))
    optimized_prediction.validate_prediction_request({
        "engine": "protenix_v2", "gpu_device": "0", "validation_mode": "complex",
        "validation_candidate_ids": ["candidate-a", "candidate-b"],
        "source_run_dir": source_run, "candidates_jsonl": candidate_index,
    })
    assert captured["validation_sources"] == [{
        "source_run_dir": str(source_run),
        "candidates_jsonl": str(candidate_index),
        "candidate_ids": ["candidate-a", "candidate-b"],
    }]
    assert "validation_candidate_ids" not in captured


def test_unavailable_selected_measurement_is_unresolved_not_passed():
    field = METRICS["iptm"]["field"]
    candidate = {"candidate_id": "missing-iptm", "source_tool": "optimized_prediction_boltz2",
                 "stage": "complex_refolding", "metrics": {}, "raw_metadata": {}}
    decision = evaluate_candidate(candidate, {
        "binder_format": "miniprotein", "target_type": "structured",
        "mandatory_filters": [{"field": field, "operator": ">=", "value": 0.7}],
    })
    assert decision["status"] == "unresolved"
    assert decision["measurements"][field]["passed"] is None


def test_acceptance_ui_uses_source_format_defaults_and_ipdae_ranking():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.design_validation_ui import _acceptance_controls
policy, target, errors = _acceptance_controls(
    "acceptance-test", "boltz2", "complex", source_tools={"bindcraft2"},
    inferred_formats=["vhh"], inferred_target_types=["fasta_idr"])
st.session_state["acceptance_policy"] = policy
st.session_state["acceptance_target"] = target
st.session_state["acceptance_errors"] = errors
''').run()
    assert not app.exception
    policy = app.session_state["acceptance_policy"]
    assert policy["binder_format"] == "vhh"
    assert policy["target_type"] == "fasta_idr"
    assert {row["metric_id"] for row in policy["mandatory_filters"]} >= {
        "iptm", "interface_residues_4a", "backbone_clashes_2_5a",
    }
    assert "target_plddt" not in {row["metric_id"] for row in policy["mandatory_filters"]}
    assert policy["ranking_metric"] == METRICS["bindcraft2_ipdae"]["field"]
    assert METRICS["bindcraft2_ipdae"]["field"] in policy["informational_scores"]
    assert policy["reference_preset"]["source"] == "BindCraft2 source 4a56313"
    assert any("Unbound_Binder_pLDDT" in item for item in policy["reference_preset"]["inactive_criteria"])
    assert app.session_state["acceptance_errors"] == []
