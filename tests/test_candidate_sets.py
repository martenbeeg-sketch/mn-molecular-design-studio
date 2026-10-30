"""Phase 7A CPU orchestration: scientific engines are deterministic workers here."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace
from zipfile import ZipFile
from io import BytesIO

import pytest

from mn_molecular_design_studio import candidate_sets as cs, candidate_pipeline as cp, pipeline_graph as pg
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as gp
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json, read_json


def structure(path, binder="GGG", target="AAA"):
    names = {"A": "ALA", "G": "GLY", "S": "SER", "V": "VAL"}
    lines = []
    for chain, seq, y in [("A", target, 0.), ("B", binder, 5.)]:
        for i, aa in enumerate(seq, 1):
            lines.append(f"ATOM  {len(lines)+1:5d}  CA  {names[aa]} {chain}{i:4d}    {i*3.8:8.3f}{y:8.3f}{0.:8.3f}  1.00 80.00           C  \n")
    path.write_text("".join(lines) + "END\n")
    return str(path)


@pytest.fixture
def collection_harness(tmp_path, monkeypatch):
    from mn_molecular_design_studio import workspace, cache_provenance, runtime
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *a, **k: None)
    monkeypatch.setattr(pg, "app_home", lambda: tmp_path / "home")
    monkeypatch.setattr(runtime, "runs_root", lambda: tmp_path / "children")
    monkeypatch.setattr(cache_provenance, "execution_identity", lambda *a, **k: {"cpu_fixture_runtime": "v1"})
    # graph_pipeline only uses time.sleep. Replace its module-local time
    # reference so the fixture does not monkeypatch Python's process-wide
    # time.sleep and interfere with Streamlit AppTest's runner waits.
    monkeypatch.setattr(gp, "time", SimpleNamespace(sleep=lambda *a, **k: None))
    target = structure(tmp_path / "target.pdb")
    from mn_molecular_design_studio.design_validation import read_structure, protein_chains, sequence
    binder_sequence = sequence(protein_chains(read_structure(Path(target)))["B"])
    parents, submissions = [], []
    def create(*args, **kwargs):
        path = tmp_path / f"pipeline-{len(parents)}"
        path.mkdir()
        write_json(path / "metadata.json", {"status": "queued", "tool": "visual_pipeline"})
        write_json(path / "input.json", {"params": {"pipeline_name": "Fixture collection graph"}})
        parents.append(path)
        return SimpleNamespace(run_dir=path)
    monkeypatch.setattr(gp, "create_job", create)
    def source(name, tool, size, *, rejected=(), backbone=()):
        path = tmp_path / name
        path.mkdir()
        rows = []
        for i in range(size):
            rows.append(dict(candidate_id=f"{name}-{i}", stage="generation.backbone" if i in backbone else "generation.backbone_sequence",
                source_tool=tool, binder_sequence=None if i in backbone else binder_sequence, complex_pdb=target, target_pdb=target,
                binder_chains=["B"], target_chains=["A"], parents=[],
                metrics={"native_pass_filters": i not in rejected, "native_score": i}, raw_metadata={"fixture": True}))
        normalized = path / "artifacts/normalized_candidates/candidates.jsonl"
        normalized.parent.mkdir(parents=True)
        normalized.write_text("".join(json.dumps(r) + "\n" for r in rows))
        write_json(path / "metadata.json", dict(status="completed", tool=tool))
        write_json(path / "input.json", dict(params={"seed": 42, "native_batch_size": 3}))
        write_json(path / "result.json", dict(success=True, metrics={"candidate_count": size}))
        write_json(path / "molecular_design_studio.json", dict(project_id="project"))
        return path
    def graph(n=5, m=2, **source_options):
        data = cp.new_pipeline()
        data["nodes"][0]["data"]["params"].update(target_pdb=target, target_chains="A")
        a = source("rfd4-source", "rfd4_proteina", n, **source_options)
        b = source("bindcraft-source", "bindcraft2", m)
        for node, run in zip(data["nodes"][1:3], [a, b]):
            node["data"]["retained_run_dir"] = str(run)
        for node in data["nodes"]:
            if node["data"]["kind"] in {"boltz2", "prediction_protenix_v2"}:
                node["data"].update(enabled=True, execution=dict(mode="sharded", batch_size=5, max_in_flight=2, retries=1, failure_policy="partial"))
                node["data"]["params"].update(seed=42)
        return data
    failures = set()
    def submit(request, parent, node, **kwargs):
        submissions.append(copy.deepcopy((request, node, kwargs)))
        path = tmp_path / "children" / "design" / ("child-" + str(len(submissions)))
        path.mkdir(parents=True)
        params = request["parameters"]
        kind = node["data"]["kind"]
        rows = []
        if kind in {"ligandmpnn", "foundry_mpnn"}:
            src = dict(run_dir=params["source_run_dir"], candidates_jsonl=params["candidates_jsonl"])
            for row in cs.iter_candidates(src):
                if row["stage"] != "generation.backbone" and node["data"]["params"].get("sequence_stage") == "initial":
                    rows.append(row)
                    continue
                for i in range(params.get("num_seq_per_target", 1)):
                    child = copy.deepcopy(row)
                    child.update(candidate_id=row["candidate_id"] + f"_{node['id']}_mpnn_{i}", parents=[row["candidate_id"]],
                                 stage="sequence_design", source_tool="protein_mpnn", binder_sequence=binder_sequence, metrics={"mpnn_nll": 1.0})
                    rows.append(child)
        else:
            for i, task in enumerate(params["tasks"]):
                rows.append(dict(candidate_id=f"{path.name}_sample_{task['source_candidate_id']}",
                    parents=[task["source_candidate_id"]], binder_sequence=task["binder_sequence"], source_tool="optimized_prediction_" + params["engine"],
                    complex_pdb=target, target_pdb=target, binder_chains=["B"], target_chains=["A"], stage="complex_refolding",
                    metrics={"iptm": .8, "plddt_mean": 80.}, raw_metadata=dict(engine=params["engine"], seed=params["seed"],
                        prediction_protocol=dict(engine=params["engine"], seed=params["seed"], checkpoint="fixture-v1", samples=1), fixture=True)))
        output = path / "artifacts/normalized_candidates/candidates.jsonl"
        output.parent.mkdir(parents=True)
        output.write_text("".join(json.dumps(r) + "\n" for r in reversed(rows)))
        fail = any(value in node["id"] for value in failures)
        write_json(path / "metadata.json", dict(status="failed" if fail else "completed", resource_allocation={"gpu_device": "fixture-no-gpu"}))
        write_json(path / "result.json", dict(success=not fail))
        write_json(path / "input.json", dict(params=params))
        return path
    monkeypatch.setattr(gp, "_submit", submit)
    return SimpleNamespace(graph=graph, submissions=submissions, parents=parents, failures=failures,
                           target=target, tmp_path=tmp_path, source=source)


@pytest.mark.parametrize("size", [0, 1, 2, 5, 10, 25, 50, 100])
def test_real_graph_collection_sizes_and_predictor_matrix(collection_harness, size):
    h = collection_harness
    graph = h.graph(size, 0)
    run = gp.run_graph_pipeline(graph, "collection sizes", "owner", "project", reuse_completed=False)
    result = read_json(run / "result.json")
    assert result["success"] is True
    rows = result["outputs"]["candidates"]
    assert len(rows) == size
    assert len({cs.metadata(r)["variant_id"] for r in rows}) == size
    assert all(len(cs.metadata(r)["predictions"]) == 2 for r in rows)
    assert all(cs.metadata(r)["original_sequence"] == r["binder_sequence"] == "GGG" for r in rows)
    execution = read_json(run / "graph_execution.json")
    for branch in ("predictor-0", "predictor-1"):
        batches = execution["nodes"][branch]["batches"]
        assert [b["unit_count"] for b in batches] == [min(5, size - i) for i in range(0, size, 5)]
        ids = [v for b in batches for v in read_json(Path(b["manifest"]))["input_candidate_ids"]]
        assert len(ids) == len(set(ids)) == size
    assert all(s[1]["data"]["kind"] not in {"ligandmpnn", "foundry_mpnn"} for s in h.submissions)


def test_different_generator_sizes_native_rejection_duplicates_and_export(collection_harness):
    h = collection_harness
    graph = h.graph(10, 3, rejected=[0])
    ranking = next(n for n in graph["nodes"] if n["data"]["kind"] == "candidate_rank")
    ranking["data"]["params"]["ranking_rules"] = [{"field": "prediction/predictor-0/iptm", "direction": "higher"}]
    run = gp.run_graph_pipeline(graph, "mixed sources", "owner", "project", reuse_completed=False)
    result = read_json(run / "result.json")
    assert result["success"] is True
    rows = result["outputs"]["candidates"]
    assert len(rows) == 13
    rejected = [r for r in rows if cs.metadata(r)["native_filter"] == "failed"]
    assert len(rejected) == 1 and not cs.metadata(rejected[0])["predictions"]
    assert all(cs.metadata(r)["deduplication"]["members"] == 13 for r in rows)
    assert sum(r["raw_metadata"]["selection"]["selected"] for r in rows) == 1
    execution = read_json(run / "graph_execution.json")
    assert [b["unit_count"] for b in execution["nodes"]["predictor-0"]["batches"]] == [5, 5, 2]
    from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import build_run_zip
    with ZipFile(BytesIO(build_run_zip(run))) as archive:
        names = set(archive.namelist())
        assert {"artifacts/candidate_set/manifest.json", "artifacts/candidate_set/candidate_index.csv", "artifacts/candidate_set/predictions.jsonl", "artifacts/candidate_set/graph.json"} <= names
        assert len(archive.read("artifacts/candidate_set/candidates.jsonl").splitlines()) == 13
        assert any(name.startswith("artifacts/candidate_set/structures/") for name in names)


def test_backbone_one_to_many_without_score_transfer(collection_harness):
    h = collection_harness
    graph = h.graph(2, 1, backbone=[0])
    node = next(n for n in graph["nodes"] if n["id"] == "sequence-preparation")
    node["data"]["params"]["sequences"] = 3
    run = gp.run_graph_pipeline(graph, "variants", "owner", "project", reuse_completed=False)
    rows = read_json(run / "result.json")["outputs"]["candidates"]
    assert len(rows) == 6  # original backbone + two sequence-ready originals + 3 children
    children = [r for r in rows if cs.metadata(r)["parent_variant_ids"]]
    assert len(children) == 3 and len({cs.metadata(r)["variant_id"] for r in children}) == 3
    assert all("native_score" not in r["metrics"] for r in children)
    assert all(len(cs.metadata(r)["predictions"]) == 2 for r in children)


def test_partial_chunks_retry_and_restart_are_idempotent(collection_harness, monkeypatch):
    h = collection_harness
    graph = h.graph(12, 0)
    h.failures.add("predictor-1_b000005")
    run = gp.run_graph_pipeline(graph, "partial", "owner", "project", reuse_completed=False)
    state = read_json(run / "graph_execution.json")
    result = read_json(run / "result.json")
    assert result["success"] is True and result["outputs"]["comparison_status"] == "partial"
    assert len(result["outputs"]["candidates"]) == 12
    # Failed collection workers can already have completed member outputs.
    # Preserve those exact results instead of rerunning every member.
    assert sum(len(cs.metadata(r)["predictions"]) for r in result["outputs"]["candidates"]) == 24
    assert state["nodes"]["predictor-1"]["candidate_progress"]["retry_count"] == 0
    before = len(h.submissions)
    monkeypatch.setattr(gp, "create_job", lambda *a, **k: SimpleNamespace(run_dir=run))
    gp.run_graph_pipeline(graph, "partial", "owner", "project", reuse_completed=False)
    assert len(h.submissions) == before
    again = read_json(run / "result.json")["outputs"]["candidates"]
    assert {r["candidate_id"] for r in again} == {r["candidate_id"] for r in result["outputs"]["candidates"]}


def test_disabled_branch_and_analyzer_change_reuses_science(collection_harness):
    h = collection_harness
    graph = h.graph(2, 3)
    next(n for n in graph["nodes"] if n["id"] == "generator-1")["data"]["enabled"] = False
    first = gp.run_graph_pipeline(graph, "first", "owner", "project")
    before = len(h.submissions)
    next(n for n in graph["nodes"] if n["id"] == "ranking")["data"]["params"]["ranking_rules"] = [{"field": "prediction/predictor-0/iptm", "direction": "higher"}]
    second = gp.run_graph_pipeline(graph, "new ranking", "owner", "project")
    assert len(h.submissions) == before
    assert len(read_json(second / "result.json")["outputs"]["candidates"]) == 2
    assert read_json(first / "result.json")["outputs"]["candidates"][0]["raw_metadata"]["selection"]["status"] == "not_evaluated"


def test_settings_and_old_workflows_are_additive():
    assert "candidate_execution" not in pg.clean_graph(pg.new_graph())
    graph = pg.new_graph("multi_generator_binder")
    assert pg.validate_graph(graph)["order"][-1] == "results"
    assert graph["nodes"][1]["data"]["kind"] == "rfd4"
    assert graph["nodes"][2]["data"]["kind"] == "bindcraft2"
    node = next(n for n in graph["nodes"] if n["id"] == "ranking")
    node["data"]["params"]["ranking_rules"] = [{"field": "iptm", "direction": "bogus"}]
    with pytest.raises(ValueError, match="Ranking"):
        pg.validate_graph(graph)


def test_corrupt_records_are_operational_failures(tmp_path):
    path = tmp_path / "records.jsonl"
    path.write_text('{"candidate_id": "ok"}\ncorrupt\n')
    with pytest.raises(ValueError, match="Invalid candidate record"):
        list(cs.iter_candidates(dict(run_dir=str(tmp_path), candidates_jsonl=str(path))))


def test_measurement_validity_and_scientific_rejection_are_separate(collection_harness):
    h = collection_harness
    graph = h.graph(5, 0)
    source = Path(graph["nodes"][1]["data"]["retained_run_dir"])
    path = source / "artifacts/normalized_candidates/candidates.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    for row, score in zip(rows, [.8, .1, None, float("nan"), .8]):
        row["metrics"]["iptm"] = score
    rows[-1]["complex_pdb"] = str(h.tmp_path / "missing.pdb")
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    # Preserve the operationally missing artifact for shared filtering; generator
    # preflight correctly refuses it, so test the transformation directly.
    target = dict(target_pdb=h.target, target_chains=["A"])
    target["identity"] = cs.target_identity(target)
    merged = cs.merge(h.tmp_path, "merged", [("source", dict(run_dir=str(source), candidates_jsonl=str(path), target=target))], {"deduplication": "mark_duplicates"})
    output = cs.shared_filter(h.tmp_path, "filtered", merged, dict(exclude_native_rejected=True, rules=[{"field": "iptm", "operator": ">=", "value": .5}]))
    outcomes = list(cs.iter_candidates(output, review=True))
    statuses = {cs.metadata(row)["native_candidate_id"]: cs.metadata(row)["shared_filter"]["status"] for row in outcomes}
    assert statuses == {"rfd4-source-0": "passed", "rfd4-source-1": "failed", "rfd4-source-2": "missing_evidence",
                        "rfd4-source-3": "invalid_metadata", "rfd4-source-4": "missing_evidence"}
    assert output["candidate_count"] == 1


def test_multiple_seeds_remain_independent_and_ambiguous_scores_are_missing(collection_harness):
    h = collection_harness
    graph = h.graph(1, 0)
    run = gp.run_graph_pipeline(graph, "seed sources", "owner", "project", reuse_completed=False)
    state = read_json(run / "graph_execution.json")["nodes"]
    variants = state["shared-filter"]["output"]
    source = state["predictor-0"]["output"]
    prediction = next(cs.iter_candidates(source))
    second = copy.deepcopy(prediction)
    second["candidate_id"] += "-seed2"
    second["raw_metadata"]["prediction_protocol"]["seed"] = 43
    cs.metadata(second)["prediction_id"] = cs.identity(cs.metadata(second)["prediction_id"], 43)
    two_seeds = cs.write_collection(h.tmp_path, "two-seeds", [second, prediction], variants["target"])
    output = cs.analyze(h.tmp_path, "join-seeds", [("variants", variants), ("boltz", two_seeds)], {}, {"predictor-0": {"enabled": True}})
    row = next(cs.iter_candidates(output))
    assert len(cs.metadata(row)["predictions"]) == 2
    assert "prediction/predictor-0/iptm" not in row["metrics"]
    assert len([key for key in row["metrics"] if key.startswith("prediction/predictor-0/") and key.endswith("/iptm")]) == 2
    bad = copy.deepcopy(second)
    bad["binder_sequence"] = "SSS"
    wrong = cs.write_collection(h.tmp_path, "wrong-sequence", [bad], variants["target"])
    with pytest.raises(ValueError, match="cross-sequence"):
        cs.analyze(h.tmp_path, "bad-join", [("variants", variants), ("prediction", wrong)], {}, {})
    conflict = copy.deepcopy(prediction)
    conflict["metrics"]["iptm"] = .2
    wrong = cs.write_collection(h.tmp_path, "conflicting-receipts", [prediction, conflict], variants["target"])
    with pytest.raises(ValueError, match="Conflicting prediction"):
        cs.analyze(h.tmp_path, "conflict-join", [("variants", variants), ("prediction", wrong)], {}, {})


def test_predictor_failure_keeps_other_evidence_and_expected_samples(collection_harness):
    h = collection_harness
    graph = h.graph(2, 0)
    predictor = next(n for n in graph["nodes"] if n["id"] == "predictor-0")
    predictor["data"]["params"]["samples"] = 2
    # Worker deliberately returns only one output per input.
    h.failures.add("predictor-1")
    run = gp.run_graph_pipeline(graph, "missing samples", "owner", "project", reuse_completed=False)
    state = read_json(run / "graph_execution.json")["nodes"]
    batch = state["predictor-0"]["batches"][0]
    assert batch["expected_output_count"] == 4 and batch["returned_output_count"] == 2
    assert batch["missing_output_count"] == 2 and batch["status"] == "partial"
    assert batch["completed_units"] == 0 and batch["returned_input_units"] == 2
    assert all(values["status"] == "partial" for values in batch["candidate_states"].values())
    rows = read_json(run / "result.json")["outputs"]["candidates"]
    assert len(rows) == 2
    assert all(len(cs.metadata(row)["predictions"]) == 2 for row in rows)
    assert all(cs.metadata(row)["predictor_matrix"]["predictor-0"] == ["completed", "missing_evidence"] for row in rows)
    assert all(cs.metadata(row)["predictor_matrix"]["predictor-1"] == ["completed"] for row in rows)


def test_preflight_rejects_wrong_target_and_project_before_jobs(collection_harness):
    h = collection_harness
    graph = h.graph(1, 0)
    wrong_target = structure(h.tmp_path / "wrong-target.pdb", target="SSS")
    graph["nodes"][0]["data"]["params"]["target_pdb"] = wrong_target
    with pytest.raises(ValueError, match="target-chain sequence"):
        gp.run_graph_pipeline(graph, "wrong target", "owner", "project", reuse_completed=False)
    assert h.parents == [] and h.submissions == []
    graph["nodes"][0]["data"]["params"]["target_pdb"] = h.target
    with pytest.raises(ValueError, match="selected project"):
        gp.run_graph_pipeline(graph, "wrong project", "owner", "unrelated-project", reuse_completed=False)
    assert h.parents == [] and h.submissions == []


def test_collection_api_catalog_save_preflight_and_plan_without_submission(collection_harness, monkeypatch):
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio import api
    h = collection_harness
    graph = h.graph(2, 1)
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(h.tmp_path / "cluster"))
    monkeypatch.setattr(api, "_request_payload_and_project", lambda payload, user: ({k: v for k, v in payload.items() if k != "project_id"}, payload["project_id"]))
    monkeypatch.setattr(api, "_validate_user_path_values", lambda *a, **k: None)
    guards = []
    monkeypatch.setattr(api, "_guard_source_job", lambda payload, *a, **k: guards.append(payload))
    source_paths = {Path(n['data']['retained_run_dir']).name: n['data']['retained_run_dir'] for n in graph['nodes'] if n['data'].get('retained_run_dir')}
    monkeypatch.setattr(api, "_resolve_api_job", lambda reference, user: {"run_dir": source_paths[reference.split('/')[-1]]})
    api.app.dependency_overrides[api._current_user] = lambda: {"id": "f" * 32, "role": "user"}
    try:
        with TestClient(api.app) as client:
            assert "candidate_analyzer" in client.get("/api/v1/pipelines/catalog").json()
            payload = dict(project_id="project", graph=graph)
            saved = client.post("/api/v1/pipelines/save", json=payload)
            assert saved.status_code == 200, saved.text
            response = client.post("/api/v1/pipelines/validate", json=payload)
            assert response.status_code == 200, response.text
            assert len(response.json()["collection_execution_plan"]) == len(graph["nodes"])
            assert guards and all("source_job" in guard for guard in guards)
            assert h.parents == [] and h.submissions == []
    finally:
        api.app.dependency_overrides.clear()


def test_collection_builder_and_result_projection_render(collection_harness):
    from streamlit.testing.v1 import AppTest
    h = collection_harness
    graph = h.graph(2, 1)
    run = gp.run_graph_pipeline(graph, "UI collection", "owner", "project", reuse_completed=False)
    fixture = h.tmp_path / "ui-graph.json"
    fixture.write_text(json.dumps(graph))
    app = AppTest.from_string(f'''
import json
from pathlib import Path
import streamlit as st
from mn_molecular_design_studio.pipeline_graph_ui import render_pipeline_builder
st.session_state['visual-pipeline-project-graph'] = json.loads(Path({str(fixture)!r}).read_text())
st.session_state['visual-pipeline-project-reset'] = 'fixture'
st.session_state['visual-pipeline-project-name'] = 'UI collection'
render_pipeline_builder(owner_id='owner', project_selector=lambda _: 'project',
    submit_api=lambda *a, **k: (_ for _ in ()).throw(AssertionError('No submission expected')),
    search_jobs=lambda **k: {{'items': [{{'run_dir': {str(run)!r}, 'run_id': 'ui-run', 'job_code': 'UIRUN', 'status': 'completed'}}]}})
from mn_molecular_design_studio.candidate_sets import iter_candidates
from mn_molecular_design_studio.result_presentation import candidate_state_rows
rows = list(iter_candidates({{'run_dir': {str(run)!r}, 'candidates_jsonl': {str(run / 'artifacts/normalized_candidates/candidates.jsonl')!r}}}))
st.dataframe(candidate_state_rows(rows))
''').run(timeout=45)
    assert not app.exception
    assert any(s.label == "Generator operations" for s in app.multiselect)
    assert any(s.label == "Enable branch" for s in app.checkbox)
    assert any("Native filter" in table.value.columns and len(table.value) == 3 for table in app.dataframe)


def test_analysis_evidence_has_its_own_identity_and_does_not_replace_prediction(collection_harness):
    h = collection_harness
    graph = h.graph(1, 0)
    run = gp.run_graph_pipeline(graph, "analysis source", "owner", "project", reuse_completed=False)
    state = read_json(run / "graph_execution.json")["nodes"]
    variants, prediction = state["shared-filter"]["output"], state["predictor-0"]["output"]
    original = next(cs.iter_candidates(prediction))
    measurement = copy.deepcopy(original)
    measurement.update(stage="analysis", source_tool="ipsae", metrics={"ipsae": .65})
    measurement["raw_metadata"]["analysis_parameters"] = {"pae_cutoff": 10, "distance_cutoff": 10}
    output = gp._local_source(h.tmp_path, "analysis-native", [measurement], variants["target"])
    attached = cs.attach_stage(h.tmp_path, "ipsae", output, prediction, category="Analysis", scientific_key="a" * 64, settings={})
    result = cs.analyze(h.tmp_path, "analyzed", [("variants", variants), ("prediction", prediction), ("ipsae", attached)], {}, {"predictor-0": {"enabled": True}})
    row = next(cs.iter_candidates(result))
    assert len(cs.metadata(row)["predictions"]) == len(cs.metadata(row)["analysis_results"]) == 1
    assert row["metrics"]["prediction/predictor-0/iptm"] == .8
    assert row["metrics"]["analysis/ipsae/ipsae"] == .65


def test_rfd4_explicit_native_refolding_option_is_preserved(collection_harness):
    h = collection_harness
    graph = h.graph(1, 0)
    node = graph["nodes"][1]
    node["data"]["params"]["options_json"] = '{"refold": true}'
    request = gp.node_request("rfd4", node["data"]["params"], {"target": dict(target_pdb=h.target, target_chains=["A"])}, "native option")
    assert request["parameters"]["options"]["refold"] is True
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_managed_paths
    assert resolve_managed_paths({"options_json": '{"refold": true}'})["options_json"] == '{"refold": true}'


def test_existing_scientific_selector_retains_collection_and_review_export(collection_harness):
    h = collection_harness
    graph = h.graph(2, 1)
    graph['nodes'].append(dict(id='existing-filter', type='studio', position=dict(x=2900,y=400),
        data=dict(kind='filter', label='Existing scientific selector', params=dict(
            ranking_metric='evidence/complex_refolding/optimized_prediction_boltz2/iptm', limit=1))))
    graph['edges'] = [edge for edge in graph['edges'] if edge['target'] != 'results']
    for source,target in [('ranking','existing-filter'),('existing-filter','results')]:
        graph['edges'].append(dict(id=source+'-'+target,source=source,target=target,sourceHandle='out',targetHandle='in'))
    run = gp.run_graph_pipeline(graph, 'existing selector', 'owner', 'project', reuse_completed=False)
    result = read_json(run/'result.json')
    assert result['success'] is True
    assert len(result['outputs']['candidates']) == 3
    assert sum(row['raw_metadata']['selection']['selected'] for row in result['outputs']['candidates']) == 1
    assert (run/'artifacts/candidate_set/manifest.json').is_file()


def test_analyzer_recovers_target_identity_for_saved_design_source(collection_harness):
    h = collection_harness
    run = gp.run_graph_pipeline(h.graph(1, 0), 'saved source', 'owner', 'project', reuse_completed=False)
    state = read_json(run / 'graph_execution.json')['nodes']
    variants = copy.deepcopy(state['shared-filter']['output'])
    prediction = copy.deepcopy(state['predictor-0']['output'])
    variants['target'].pop('identity')
    prediction['target'].pop('identity')
    sources = gp._candidate_sources_with_target_identity([('variants', variants), ('predictor-0', prediction)])
    output = cs.analyze(h.tmp_path, 'saved-source-join', sources, {}, {})
    row = next(cs.iter_candidates(output))
    assert output['target']['identity']['target_id'] == cs.metadata(row)['target_id']
    assert len(cs.metadata(row)['predictions']) == 1
    assert 'identity' not in variants['target']


def test_campaign_validation_is_projected_into_candidate_set_without_overwriting_native_measurements(collection_harness):
    h = collection_harness
    source = h.source("campaign-source", "rfd4_proteina", 1)
    target_structure = h.target
    native = next(json.loads(line) for line in (source / "artifacts/normalized_candidates/candidates.jsonl").read_text().splitlines())
    native.update(candidate_id="native-candidate-1", metrics={"iptm": .22, "native_score": 3.0})
    native["raw_metadata"] = {
        "distributed_campaign": {
            "campaign_id": "campaign-1", "attempt_id": "campaign-1:a00000001",
            "attempt_seed": 1234, "generator_id": "rfd4_proteina", "objective_type": "accepted",
            "native_acceptance": "unavailable", "validated_acceptance": "accepted",
            "validation_decision_reason": "Boltz-2 iPTM 0.72 met the configured threshold.",
        },
        "evidence_candidates": [{
            "candidate_id": "boltz2-result-1", "stage": "complex_refolding",
            "source_tool": "optimized_prediction_boltz2", "binder_sequence": native["binder_sequence"],
            "complex_pdb": target_structure, "target_pdb": target_structure,
            "binder_chains": ["B"], "target_chains": ["A"], "metrics": {"iptm": .72},
            "raw_metadata": {"distributed_campaign_validation": {
                "campaign_id": "campaign-1", "attempt_id": "campaign-1:a00000001",
                "stage_index": 0, "workflow": "refolding.optimized_prediction", "status": "completed",
            }},
        }, {
            "candidate_id": "ipsae-result-1", "stage": "analysis",
            "source_tool": "ipsae", "binder_sequence": native["binder_sequence"],
            "complex_pdb": target_structure, "target_pdb": target_structure,
            "binder_chains": ["B"], "target_chains": ["A"], "metrics": {"ipsae": .61},
            "raw_metadata": {"distributed_campaign_validation": {
                "campaign_id": "campaign-1", "attempt_id": "campaign-1:a00000001",
                "stage_index": 1, "workflow": "analysis.ipsae", "status": "completed",
            }},
        }],
    }
    candidate_file = source / "artifacts/normalized_candidates/candidates.jsonl"
    candidate_file.write_text(json.dumps(native) + "\n")
    target = {"target_pdb": target_structure, "target_chains": ["A"]}
    target["identity"] = cs.target_identity(target)
    merged = cs.merge(h.tmp_path, "campaign-merged", [
        ("campaign", dict(run_dir=str(source), candidates_jsonl=str(candidate_file), target=target)),
    ], {"deduplication": "mark_duplicates"})
    campaign_row = next(cs.iter_candidates(merged, review=True))
    info = cs.metadata(campaign_row)
    assert info["native_filter"] == "not_evaluated"
    assert info["native_measurements"] == {"iptm": .22, "native_score": 3.0}
    assert [(record["branch"], record["status"], record["measurements"]["iptm"]["value"])
            for record in info["predictions"]] == [("campaign-validation-0", "completed", .72)]
    assert [(record["branch"], record["status"], record["measurements"]["ipsae"]["value"])
            for record in info["analysis_results"]] == [("campaign-validation-1", "completed", .61)]
    assert cs.observation(campaign_row, "prediction/campaign-validation-0/iptm")["value"] == .72
    assert cs.observation(campaign_row, "analysis/campaign-validation-1/ipsae")["value"] == .61

    analyzed = cs.analyze(h.tmp_path, "campaign-analyzed", [("campaign", merged)], {}, {
        "campaign-validation-0": {"enabled": True},
        "campaign-validation-1": {"enabled": True},
    })
    result = next(cs.iter_candidates(analyzed, review=True))
    assert result["metrics"]["iptm"] == .22
    assert result["metrics"]["prediction/campaign-validation-0/iptm"] == .72
    assert result["metrics"]["analysis/campaign-validation-1/ipsae"] == .61
    assert len(cs.metadata(result)["predictions"]) == 1
    assert len(cs.metadata(result)["analysis_results"]) == 1
    from mn_molecular_design_studio.result_presentation import candidate_state_rows
    visible = candidate_state_rows([result])[0]
    assert visible["Campaign acceptance"] == "accepted"
    assert visible["Campaign acceptance level"] == "validated"
    assert visible["Selection"] == "not_evaluated" and visible["Selected"] is False


def test_node_bound_acceptance_reuses_exact_independent_decision_without_reading_native_score(collection_harness):
    from mn_molecular_design_studio.validation_acceptance import METRICS, evaluate_candidate, campaign_policy_from_node
    h = collection_harness
    source = h.source('designs', 'rfd4_proteina', 1)
    path = source / 'artifacts/normalized_candidates/candidates.jsonl'
    candidate = json.loads(path.read_text())
    params = {'binder_format': 'miniprotein', 'target_type': 'structured',
              'rules': [{'field': METRICS['iptm']['field'], 'operator': '>=', 'value': .7}],
              'exclude_native_rejected': False}
    policy = campaign_policy_from_node(params)
    candidate['metrics']['iptm'] = .1  # Native measurement has a different meaning/source.
    prediction = {**candidate, 'source_tool': 'optimized_prediction_boltz2', 'stage': 'complex_refolding',
                  'metrics': {'iptm': .8}}
    candidate['raw_metadata']['validation_acceptance'] = evaluate_candidate(prediction, policy)
    candidate['raw_metadata']['validation_acceptance_provenance'] = {'acceptance_node': 'filter', 'validation_run_dir': 'independent-run'}
    path.write_text(json.dumps(candidate) + '\n')
    target = {'target_pdb': h.target, 'target_chains': ['A']}
    target['identity'] = cs.target_identity(target)
    merged = cs.merge(h.tmp_path, 'merged', [('generator', {'run_dir': str(source), 'candidates_jsonl': str(path), 'target': target})], {'deduplication': 'keep_all'})
    output = cs.shared_filter(h.tmp_path, 'filter', merged, params)
    rows = list(cs.iter_candidates(output))
    assert len(rows) == 1
    assert rows[0]['metrics']['iptm'] == .1
    measured = rows[0]['raw_metadata']['validation_acceptance']['measurements'][METRICS['iptm']['field']]
    assert measured['value'] == .8
    assert cs.metadata(rows[0])['shared_filter']['status'] == 'passed'
