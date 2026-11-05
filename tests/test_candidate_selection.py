import copy
import json

import pytest

from mn_molecular_design_studio.candidate_selection import evaluate_selection, selected_candidates, write_selection_artifacts


def candidate(identifier, score=.8, **kwargs):
    return {"candidate_id": identifier, "stage": "analysis", "source_tool": "assay", "binder_sequence": "ACD:EFG",
            "metrics": {"score": score, **kwargs}, "raw_metadata": {"analysis_parameters": {"protocol": "one"}, "target_structure_digest": "target"}}


def policy(**kwargs):
    return {"rules": [{"field": "score", "operator": ">=", "value": .5}], "ranking_metric": "score", "direction": "higher", **kwargs}


def test_pass_reject_missing_and_failure_are_distinct_and_sources_unchanged():
    rows = [candidate("good"), candidate("reject", .2), candidate("missing", None), candidate("failed", .9, qc_error="original failure")]
    before = copy.deepcopy(rows)
    outcomes = evaluate_selection(rows, policy())
    assert [row["raw_metadata"]["selection"]["status"] for row in outcomes] == ["pass", "reject", "missing", "failed"]
    assert [row["candidate_id"] for row in selected_candidates(outcomes)] == ["good"]
    assert "original failure" in outcomes[-1]["raw_metadata"]["selection"]["reasons"][0]
    assert rows == before


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), float("-inf"), True, "bad"])
def test_invalid_values_never_pass(value):
    row = evaluate_selection([candidate("bad", value)], policy())[0]
    assert row["raw_metadata"]["selection"]["status"] == "missing"
    assert not selected_candidates([row])


def test_finite_but_invalid_metric_domains_never_pass():
    row = candidate("invalid"); row["metrics"] = {"iptm": 4.}
    assert not selected_candidates(evaluate_selection([row], {"required_metrics": ["iptm"]}))


def test_prediction_protocol_differences_prevent_cross_protocol_ranking():
    rows = [candidate("a"), candidate("b")]
    rows[0]["raw_metadata"]["prediction_protocol"] = {"num_recycles": 1}
    rows[1]["raw_metadata"]["prediction_protocol"] = {"num_recycles": 10}
    assert len(selected_candidates(evaluate_selection(rows, policy(limit=1)))) == 2


def test_common_measurement_direction_uses_declared_contract():
    rows = [candidate("a", rosetta_interface_dG=-10.), candidate("b", rosetta_interface_dG=-20.)]
    settings = {"ranking_metric": "common.energy.rosetta.interface_dg", "limit": 1}
    assert [row["candidate_id"] for row in selected_candidates(evaluate_selection(rows, settings))] == ["b"]
    with pytest.raises(ValueError, match="lower values first"):
        evaluate_selection(rows, {**settings, "direction": "higher"})


def test_ties_are_deterministic_and_limits_are_not_rejections():
    rows = [candidate("c"), candidate("b"), candidate("a")]
    first = evaluate_selection(rows, policy(limit=1))
    second = evaluate_selection(list(reversed(rows)), policy(limit=1))
    assert [row["candidate_id"] for row in selected_candidates(first)] == ["a"]
    assert [row["candidate_id"] for row in selected_candidates(second)] == ["a"]
    assert all(row["raw_metadata"]["selection"]["status"] == "pass" for row in first)
    assert len(first) == 3


def test_ranking_and_limits_apply_within_comparable_protocols():
    rows = [candidate("a", .8), candidate("b", .9), candidate("c", .1)]
    rows[-1]["raw_metadata"]["analysis_parameters"] = {"protocol": "two"}
    outcomes = evaluate_selection(rows, policy(rules=[], limit=1, direction="lower"))
    assert {row["candidate_id"] for row in selected_candidates(outcomes)} == {"a", "c"}
    assert outcomes[0]["raw_metadata"]["selection"]["cohort"] != outcomes[-1]["raw_metadata"]["selection"]["cohort"]


def test_unknown_protocols_are_not_assumed_comparable():
    rows = [candidate("a"), candidate("b")]
    for row in rows:
        row["raw_metadata"] = {}
    assert len(selected_candidates(evaluate_selection(rows, policy(limit=1)))) == 2


def test_diversity_quota_uses_existing_cluster_evidence():
    rows = [candidate("a", .9, diversity_cluster=1), candidate("b", .8, diversity_cluster=1), candidate("c", .7, diversity_cluster=2)]
    outcomes = evaluate_selection(rows, policy(diversity_field="diversity_cluster", max_per_cluster=1))
    assert [row["candidate_id"] for row in selected_candidates(outcomes)] == ["a", "c"]
    assert outcomes[1]["raw_metadata"]["selection"]["reasons"] == ["Diversity cluster quota reached"]


def test_no_policy_does_not_establish_scientific_acceptance_and_manual_selection_is_available():
    rows = [candidate("a"), candidate("b")]
    assert all(row["raw_metadata"]["selection"]["status"] == "not_evaluated" for row in evaluate_selection(rows, {}))
    manual = selected_candidates(evaluate_selection(rows, {"manual_ids": ["b"]}))
    assert [row["candidate_id"] for row in manual] == ["b"]
    assert manual[0]["raw_metadata"]["selection"]["selection_method"] == "manual"


def test_between_and_any_logic_separate_failure_from_unresolved():
    rows = [candidate("in-range", .7, quality=.1), candidate("out-of-range", .9, quality=.1),
            candidate("any-pass", .9, quality=.6), candidate("unresolved", None)]
    outcomes = evaluate_selection(rows, {
        "metric": "score", "operator": "between", "threshold": .5, "threshold_max": .8,
        "rules": [{"field": "quality", "operator": ">=", "value": .5}], "rule_logic": "any",
    })
    status = {row["candidate_id"]: row["raw_metadata"]["selection"]["status"] for row in outcomes}
    assert status == {"in-range": "pass", "out-of-range": "reject", "any-pass": "pass", "unresolved": "missing"}
    # The primary range and additional rule are OR alternatives under Any.
    assert {row["candidate_id"] for row in selected_candidates(outcomes)} == {"in-range", "any-pass"}


def test_quantile_shortlist_is_protocol_and_parent_scoped():
    rows = [candidate("a1", .9), candidate("a2", .7), candidate("b1", .8), candidate("b2", .2)]
    rows[0]["parents"] = ["parent-a"]; rows[1]["parents"] = ["parent-a"]
    rows[2]["parents"] = ["parent-b"]; rows[3]["parents"] = ["parent-b"]
    selected = selected_candidates(evaluate_selection(rows, {
        "ranking_metric": "score", "direction": "higher", "quantile": .5,
        "quantile_parent_based": True,
    }))
    assert {row["candidate_id"] for row in selected} == {"a1", "b1"}


def test_stage_specific_rules_combine_prediction_and_qc():
    from test_candidate_evidence import final_candidate
    rules = [{"field": "evidence/complex_refolding/boltz2/iptm", "operator": ">=", "value": .75},
             {"field": "evidence/analysis/sequence_qc/qc_sequence_gravy", "operator": "<=", "value": 0.}]
    row = final_candidate()
    assert selected_candidates(evaluate_selection([row], {"rules": rules}))
    rules[1]["value"] = -.5
    rejected = evaluate_selection([row], {"rules": rules})
    assert rejected[0]["raw_metadata"]["selection"]["status"] == "reject"
    assert rejected[0]["raw_metadata"]["source_candidate"]["parents"] == ["prediction"]


def test_export_preserves_sequences_structures_provenance_and_rejected_rows(tmp_path):
    structure = tmp_path / "source.cif"; structure.write_bytes(b"exact structure bytes\n")
    good, rejected = candidate("good"), candidate("reject", .1)
    good["complex_pdb"] = str(structure)
    outcomes = evaluate_selection([good, rejected], policy())
    destination = tmp_path / "export"
    write_selection_artifacts(destination, outcomes)
    assert len((destination / "all_candidates.jsonl").read_text().splitlines()) == 2
    selected = json.loads((destination / "selected_candidates.jsonl").read_text())
    assert (destination / selected["complex_pdb"]).read_bytes() == structure.read_bytes()
    assert selected["binder_sequence"] == "ACD:EFG"
    assert selected["raw_metadata"]["analysis_parameters"] == {"protocol": "one"}
    assert "\nACD\n" in (destination / "selected_sequences.fasta").read_text()
    assert "\nEFG\n" in (destination / "selected_sequences.fasta").read_text()
    from mn_molecular_design_studio.result_presentation import candidate_state_rows
    visible = candidate_state_rows(outcomes)
    assert visible[0]["Selected"] is True
    assert visible[1]["Selection"] == "reject"
    assert "fails" in visible[1]["Selection reasons"]


@pytest.mark.parametrize("settings", [{"rules": [{"field": "x", "operator": ">=", "value": float("nan")}]},
    {"manual_ids": ["unknown"]}, {"ranking_metric": "analysis_score"}, {"limit": -1}])
def test_invalid_policy_is_rejected_before_selection(settings):
    with pytest.raises(ValueError):
        evaluate_selection([candidate("one")], settings)


def test_pipeline_filter_results_keeps_rejected_review_rows(tmp_path, monkeypatch, design_source_factory):
    import uuid
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, SCHEMA
    from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    source, row = design_source_factory()
    rows = [{**row, "candidate_id": "good", "metrics": {"iptm": .8}},
            {**row, "candidate_id": "rejected", "metrics": {"iptm": .2}}]
    path = source / "artifacts/normalized_candidates/candidates.jsonl"
    path.write_text("".join(json.dumps(item) + "\n" for item in rows))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "load_organization", lambda home=None: {"projects": [{"id": "test-project"}]})
    monkeypatch.setattr(graph_pipeline, "scope_root", lambda *args: tmp_path / "scope")
    nodes = []
    for index, kind in enumerate(("project_design_source", "filter", "results")):
        spec = NODE_CATALOG[kind]
        params = {field["name"]: copy.deepcopy(field["default"]) for field in spec["fields"]}
        if kind == "project_design_source":
            params.update(source_run_dir=str(source), selected_candidate_ids=["good", "rejected"])
        if kind == "filter":
            params.update(metric="iptm", threshold=.5)
        nodes.append({"id": kind, "position": {"x": index * 300, "y": 100},
                      "data": {"kind": kind, "label": spec["label"], "params": params}})
    graph = {"schema_version": SCHEMA, "id": uuid.uuid4().hex, "name": "selection", "nodes": nodes,
        "edges": [{"id": "first", "source": "project_design_source", "target": "filter"},
                  {"id": "second", "source": "filter", "target": "results"}]}
    run = graph_pipeline.run_graph_pipeline(graph, "selection", "test-user", "test-project", reuse_completed=False, gpu_device="none")
    assert read_json(run / "result.json")["success"] is True
    final = read_candidates(run)
    assert len(final) == 2
    assert sorted(item["raw_metadata"]["selection"]["status"] for item in final) == ["pass", "reject"]
    assert sum(item["raw_metadata"]["selection"]["selected"] for item in final) == 1
    selection_dir = run / "artifacts/pipeline_nodes/filter/artifacts/selection"
    assert len((selection_dir / "selected_candidates.jsonl").read_text().splitlines()) == 1
    assert (run / "artifacts/normalized_candidates/final_evidence.csv").is_file()


def test_filter_routes_pass_fail_and_unresolved_candidates_through_canvas(tmp_path, monkeypatch, design_source_factory):
    import uuid
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, SCHEMA
    from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

    source, row = design_source_factory()
    rows = [{**row, "candidate_id": "good", "metrics": {"iptm": .8}},
            {**row, "candidate_id": "rejected", "metrics": {"iptm": .2}},
            {**row, "candidate_id": "missing", "metrics": {}}]
    (source / "artifacts/normalized_candidates/candidates.jsonl").write_text(
        "".join(json.dumps(item) + "\n" for item in rows))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "load_organization", lambda home=None: {"projects": [{"id": "test-project"}]})
    monkeypatch.setattr(graph_pipeline, "scope_root", lambda *args: tmp_path / "scope")

    nodes = []
    for index, (kind, node_id) in enumerate((
        ("project_design_source", "source"), ("filter", "filter"),
        ("results", "pass-results"), ("results", "fail-results"), ("results", "unresolved-results"),
    )):
        spec = NODE_CATALOG[kind]
        params = {field["name"]: copy.deepcopy(field["default"]) for field in spec["fields"]}
        if kind == "project_design_source":
            params.update(source_run_dir=str(source), selected_candidate_ids=["good", "rejected", "missing"])
        if kind == "filter":
            params.update(metric="iptm", threshold=.5)
        nodes.append({"id": node_id, "position": {"x": index * 300, "y": 100},
                      "data": {"kind": kind, "label": spec["label"], "params": params}})
    graph = {"schema_version": SCHEMA, "id": uuid.uuid4().hex, "name": "filter-routing", "nodes": nodes,
        "edges": [
            {"id": "source-filter", "source": "source", "target": "filter"},
            {"id": "pass", "source": "filter", "target": "pass-results", "sourceHandle": "out"},
            {"id": "fail", "source": "filter", "target": "fail-results", "sourceHandle": "fail"},
            {"id": "unresolved", "source": "filter", "target": "unresolved-results", "sourceHandle": "unresolved"},
        ]}
    run = graph_pipeline.run_graph_pipeline(graph, "filter-routing", "test-user", "test-project",
                                             reuse_completed=False, gpu_device="none")
    state = read_json(run / "graph_execution.json")["nodes"]
    routed = {}
    for result_id in ("pass-results", "fail-results", "unresolved-results"):
        output = state[result_id]["output"]
        result_rows = read_candidates(output["candidates_jsonl"])
        routed[result_id] = [item["raw_metadata"].get("pipeline_original_candidate_id") for item in result_rows]
    # A Results node reached by the pass branch retains the complete review
    # table; fail/unresolved branches remain separately inspectable.
    assert routed == {"pass-results": ["good", "rejected", "missing"],
                      "fail-results": ["rejected"], "unresolved-results": ["missing"]}
    filter_output = state["filter"]["output"]
    assert [item["candidate_id"] for item in read_candidates(filter_output["candidates_jsonl"])] == ["good"]
    all_rows = read_candidates(filter_output["review_candidates_jsonl"])
    assert {item["raw_metadata"]["selection"]["routed_to"] for item in all_rows} == {"out", "fail", "unresolved"}
    assert {item["candidate_id"] for item in all_rows} == {"good", "rejected", "missing"}
