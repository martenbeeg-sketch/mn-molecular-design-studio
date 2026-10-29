"""Format-aware independent validation and acceptance contracts."""
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


def test_idr_and_peptide_preset_disables_free_fold_filters_but_requires_target_confidence():
    metrics = dict(available_metrics("peptide", "fasta_idr", "complex", engine="boltz2"))
    rules, unavailable = bindcraft2_reference_policy("peptide", "fasta_idr", "complex", metrics)
    by_id = {row["metric_id"]: row for row in rules}
    assert "ptm" not in by_id
    assert "target_plddt" in by_id and by_id["target_plddt"]["value"] == 60.0
    assert by_id["iptm"]["value"] == 0.70
    assert any("i_pAE ≤0.33" in item for item in unavailable)
    assert any("Unbound_Binder_pLDDT" in item for item in unavailable)


def test_ipae_and_ipsae_are_not_offered_as_interchangeable_filters():
    boltz = dict(available_metrics("miniprotein", "structured", "complex", engine="boltz2"))
    af3 = dict(available_metrics("miniprotein", "structured", "complex", engine="alphafold3"))
    assert "ipsae_min" in boltz
    assert "ipsae_min" not in af3
    rules, unavailable = bindcraft2_reference_policy("miniprotein", "structured", "complex", boltz)
    assert all(row["metric_id"] != "ipsae_min" for row in rules)
    assert any("i_pAE" in item for item in unavailable)


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
        "iptm", "interface_residues_4a", "backbone_clashes_2_5a", "target_plddt",
    }
    assert policy["ranking_metric"] == METRICS["bindcraft2_ipdae"]["field"]
    assert METRICS["bindcraft2_ipdae"]["field"] in policy["informational_scores"]
    assert policy["reference_preset"]["source"] == "BindCraft2 source 4a56313"
    assert any("Unbound_Binder_pLDDT" in item for item in policy["reference_preset"]["inactive_criteria"])
    assert app.session_state["acceptance_errors"] == []
