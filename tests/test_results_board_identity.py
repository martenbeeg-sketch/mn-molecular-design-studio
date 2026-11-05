"""Results board identity fields preserve explicit model zero and route choice."""
from mn_molecular_design_studio.complex_analysis_results_ui import (
    _first_present,
    _measurement_label,
    _source_candidate_id,
    board_entries,
)
from mn_molecular_design_studio.complex_analysis_ui import _biotite_item_matches_source, initialize_workspace_view


def test_model_zero_is_not_treated_as_missing_and_manual_route_is_visible():
    assert _first_present(None, 0, "later") == 0
    row = {"candidate_id": "candidate-model-0", "source_tool": "boltz2",
           "model_index": 0, "metrics": {}, "raw_metadata": {}}
    entry = board_entries(row, label="candidate-model-0", run_id="run", digest="sha256",
                          manual_selection=True)[0]
    assert entry["manual_selection"] is True


def test_results_board_displays_native_identity_and_explains_missing_independent_measurement():
    row = {"candidate_id": "cs_studio-id", "source_tool": "bindcraft2", "metrics": {},
           "raw_metadata": {"candidate_set": {"native_candidate_id": "bindcraft2_00003"}}}
    entry = board_entries(row, label="candidate", run_id="run", digest="sha256")[0]
    assert _source_candidate_id(row) == "bindcraft2_00003"
    assert _measurement_label(row) == "No independent measurement recorded"
    assert entry["source_candidate_id"] == "bindcraft2_00003"
    assert entry["status"] == "not_measured"
    assert entry["measurement_label"] == "No independent measurement recorded"


def test_biotite_receipt_joins_candidate_set_output_to_native_results_identity():
    item = {"candidate_id": "cs-variant-model-0", "structure_digest": "sha256:abc",
            "binder_chains": ["B"], "target_chains": ["A"]}
    outputs = {"cs-variant-model-0": {"candidate_id": "cs-variant-model-0", "raw_metadata": {
        "complex_analysis_entry": {"candidate_id": "cs-variant-model-0",
            "native_candidate_id": "bindcraft2_00003", "source_run": "r3-run",
            "structure_sha256": "sha256:abc", "model_id": "0"}}}}
    args = {"candidate_id": "bindcraft2_00003", "structure_digest": "sha256:abc",
            "binder_chains": ["B"], "target_chains": ["A"], "source_run_id": "r3-run"}
    matched, receipt = _biotite_item_matches_source(item, outputs, **args)
    assert matched is True and receipt["model_id"] == "0"
    assert _biotite_item_matches_source(item, outputs, **{**args, "structure_digest": "sha256:other"}) == (False, {})
    assert _biotite_item_matches_source(item, outputs, **{**args, "source_run_id": "other-run"}) == (False, {})


def test_complex_analysis_initialization_preserves_explicit_results_board_deep_link():
    state = {"complex-analysis-workspace-view": "results_board"}
    initialize_workspace_view(state)
    assert state == {"complex-analysis-workspace-view": "results_board",
                     "complex-analysis-dashboard-version": "run-first-v1"}
