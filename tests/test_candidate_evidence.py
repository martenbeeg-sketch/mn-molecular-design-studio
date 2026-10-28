import copy

from mn_molecular_design_studio.output_vocabulary import candidate_evidence, evidence_rows, metric_value


def final_candidate():
    parent = {"candidate_id": "parental-antibody", "source_tool": "import", "stage": "sequence_design",
              "binder_sequence": "ACDE", "metrics": {}, "parents": []}
    prediction = {"candidate_id": "prediction", "source_tool": "boltz2", "stage": "complex_refolding",
                  "binder_sequence": "ACDF", "metrics": {"iptm": .8}, "parents": ["parental-antibody"],
                  "raw_metadata": {"source_candidate": parent}}
    ipsae = {"candidate_id": "interface", "source_tool": "ipsae", "stage": "analysis",
             "binder_sequence": "ACDF", "metrics": {"ipsae": .7, "analysis_status": "measured"},
             "parents": ["prediction"], "raw_metadata": {"source_candidate": prediction}}
    return {"candidate_id": "qc", "source_tool": "sequence_qc", "stage": "analysis", "binder_sequence": "ACDF",
            "metrics": {"qc_sequence_gravy": -.2}, "parents": ["interface"], "raw_metadata": {"source_candidate": ipsae}}


def test_final_qc_can_select_prediction_and_analysis_without_copying_metrics():
    row = final_candidate(); before = copy.deepcopy(row)
    assert metric_value(row, "evidence/complex_refolding/boltz2/iptm") == .8
    assert metric_value(row, "evidence/analysis/ipsae/ipsae") == .7
    assert metric_value(row, "evidence/analysis/sequence_qc/qc_sequence_gravy") == -.2
    assert metric_value(row, "iptm") is None  # legacy current-stage lookup is unchanged
    assert row == before


def test_lineage_includes_parental_comparator_but_cannot_borrow_its_scores():
    row = final_candidate()
    parent = row["raw_metadata"]["source_candidate"]["raw_metadata"]["source_candidate"]["raw_metadata"]["source_candidate"]
    parent["metrics"] = {"iptm": .99}
    records = candidate_evidence(row)
    assert [record["candidate_id"] for record in records] == ["qc", "interface", "prediction", "parental-antibody"]
    assert records[-1]["role"] == "comparator"
    assert records[-1]["candidate"]["binder_sequence"] == "ACDE"
    assert metric_value(row, "evidence/sequence_design/import/iptm") is None
    assert {entry["measurement_candidate_id"] for entry in evidence_rows([row])} == {"qc", "interface", "prediction", "parental-antibody"}


def test_distributed_campaign_evidence_candidates_are_visible_without_relabeling():
    evidence = {"candidate_id": "prediction-1", "source_tool": "optimized_prediction_boltz2",
                "stage": "complex_refolding", "binder_sequence": "KLMNP", "metrics": {"iptm": .61},
                "parents": ["campaign-1"], "raw_metadata": {"prediction_protocol": {"engine": "boltz2"}}}
    candidate = {"candidate_id": "campaign-1", "source_tool": "bindcraft2",
                 "stage": "generation.backbone_sequence", "binder_sequence": "KLMNP",
                 "metrics": {}, "raw_metadata": {"evidence_candidates": [evidence]}}
    records = candidate_evidence(candidate)
    assert any(record["candidate_id"] == "prediction-1"
               and record["stage"] == "complex_refolding"
               and record["source_tool"] == "optimized_prediction_boltz2"
               and record["same_sequence"] is True for record in records)
    assert metric_value(candidate, "evidence/complex_refolding/optimized_prediction_boltz2/iptm") == .61


def test_identical_metric_names_and_ambiguous_stages_are_never_overwritten():
    row = final_candidate()
    row["metrics"]["iptm"] = .1
    prediction = row["raw_metadata"]["source_candidate"]["raw_metadata"]["source_candidate"]
    prediction["raw_metadata"]["source_candidate"] = {
        **copy.deepcopy(prediction), "candidate_id": "earlier-prediction", "raw_metadata": {}}
    assert metric_value(row, "iptm") == .1
    assert metric_value(row, "evidence/complex_refolding/boltz2/iptm") is None
    assert metric_value(row, "evidence/complex_refolding/boltz2/prediction/iptm") == .8


def test_failed_stage_does_not_supply_stale_measurements():
    row = final_candidate()
    row["raw_metadata"]["source_candidate"]["metrics"]["ipsae_error"] = "original scorer cause"
    assert metric_value(row, "evidence/analysis/ipsae/ipsae") is None
    assert candidate_evidence(row)[1]["error"] == "original scorer cause"


def test_optimized_prediction_engine_has_scoped_common_metrics():
    row = {"stage": "complex_refolding", "source_tool": "optimized_prediction_boltz2", "metrics": {"iptm": .8}}
    assert metric_value(row, "common.confidence.validation.complex.iptm") == .8


def test_real_analysis_adapter_does_not_relabel_inherited_scores(tmp_path, monkeypatch, design_source_factory):
    import json
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    source, candidate = design_source_factory()
    candidate["metrics"]["iptm"] = .8
    path = source / "artifacts/normalized_candidates/candidates.jsonl"
    path.write_text(json.dumps(candidate) + "\n")
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    # Missing PAE fails naturally before any Docker execution.
    run = analysis.run_ipsae_analysis_contract(source, path)
    child = read_candidates(run)[0]
    assert "iptm" not in child["metrics"]
    assert child["parents"] == [candidate["candidate_id"]]
    assert child["raw_metadata"]["source_candidate"]["metrics"]["iptm"] == .8
    assert metric_value(child, "evidence/complex_refolding/bindcraft2/iptm") == .8


def test_results_materializes_provenance_table_without_changing_candidate(tmp_path):
    import pandas as pd
    from pathlib import Path
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import _local_source
    row = final_candidate(); before = copy.deepcopy(row)
    output = _local_source(tmp_path, "results", [row], final_evidence=True)
    table = pd.read_csv(Path(output["candidates_jsonl"]).with_name("final_evidence.csv"))
    assert .8 in table.loc[table["field"] == "iptm", "value"].astype(float).tolist()
    assert row == before
