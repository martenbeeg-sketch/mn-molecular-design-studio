import json
from pathlib import Path

import pytest

from test_api_integrations import api_data, headers
from test_api_extended import candidate_run


def sources(data):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    from mn_molecular_design_studio.job_catalog import _upsert
    from mn_molecular_design_studio.runtime import runs_root
    paths = [candidate_run(data, tool="genie3"), candidate_run(data, tool="pxdesign")]
    for path in paths:
        write_json(path / "metadata.json", {**read_json(path / "metadata.json"),
            "campaign_group_id": "r3-design", "campaign_group_name": "R3 design"})
        _upsert(runs_root(), path, refresh_artifacts=False)
    return paths


def test_api_lists_one_run_exports_both_duplicate_native_ids_and_enforces_access(api_data):
    data = api_data
    paths = sources(data)
    client = data["client"]
    project = data["project"]["id"]
    response = client.get("/api/v1/runs", params={"project_id": project})
    assert response.status_code == 200, response.text
    run = next(r for r in response.json()["items"] if r["name"] == "R3 design")
    assert run["job_count"] == 2 and run["outputs"][0]["source_count"] == 2
    endpoint = "/api/v1/runs/" + run["id"] + "/results"
    result = client.get(endpoint, params={"project_id": project})
    assert result.status_code == 200, result.text
    payload = result.json()
    assert payload["complete"] and payload["candidate_count"] == 2
    assert len({c["run_candidate_id"] for c in payload["candidates"]}) == 2
    assert {c["candidate_id"] for c in payload["candidates"]} == {"generated-1"}
    exported = client.get(endpoint, params={"project_id": project, "export": "jsonl"})
    assert len(exported.text.splitlines()) == 2
    assert all(json.loads(line)["metrics"]["rosetta_total_score"] == -9.5 for line in exported.text.splitlines())
    denied = client.get(endpoint, params={"project_id": project}, headers=headers(data, "other"))
    assert denied.status_code in {403, 404}


def test_api_whole_run_filter_submits_one_named_output_and_replays_no_jobs(api_data, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as workflows
    from mn_molecular_design_studio.validation_acceptance import METRICS
    data = api_data
    sources(data)
    client = data["client"]
    project = data["project"]["id"]
    run = next(r for r in client.get("/api/v1/runs", params={"project_id": project}).json()["items"] if r["name"] == "R3 design")
    original = workflows.enqueue_workflow_call
    def queued(*args, **kwargs):
        kwargs["dispatch"] = False
        return original(*args, **kwargs)
    monkeypatch.setattr(workflows, "enqueue_workflow_call", queued)
    request = dict(schema_version=1, project_id=project, workflow="analysis.validation_acceptance",
        source_job=run["id"], parameters={"job_name": "IDR acceptance", "acceptance_policy": {
            "binder_format": "miniprotein", "target_type": "fasta_idr",
            "mandatory_filters": [{"field": METRICS["iptm"]["field"], "operator": ">=", "value": 0.7}]}})
    validation = client.post("/api/v1/workflows/validate", json=request)
    assert validation.status_code == 200 and validation.json()["source_count"] == 2, validation.text
    response = client.post("/api/v1/workflows", json=request, headers={"Idempotency-Key": "whole-run-filter"})
    assert response.status_code == 202, response.text
    result = response.json()
    assert result["submitted_job_count"] == 2 and not result["dispatch_errors"]
    assert client.post("/api/v1/workflows", json=request, headers={"Idempotency-Key": "whole-run-filter"}).json() == result
    # Direct reindex closes the fixture's short catalog debounce interval.
    from mn_molecular_design_studio.job_catalog import _upsert
    from mn_molecular_design_studio.runtime import runs_root
    for job in result["jobs"]:
        _upsert(runs_root(), Path(job["run_dir"]), refresh_artifacts=False)
    updated = next(r for r in client.get("/api/v1/runs", params={"project_id": project}).json()["items"] if r["id"] == run["id"])
    assert updated["job_count"] == 4
    assert len(updated["outputs"]) == 2
    analyzed = next(o for o in updated["outputs"] if o["reference"] == result["run_output"])
    assert analyzed["name"] == "IDR acceptance" and analyzed["source_count"] == 2


def test_run_graph_freezes_one_native_request_per_source_without_copying_data(api_data):
    from mn_molecular_design_studio import logical_runs, run_graph, candidate_execution
    from test_pipeline_graph import node as new_node
    data = api_data
    sources(data)
    run = next(r for r in logical_runs.catalog(data["admin"], project_id=data["project"]["id"]) if r["name"] == "R3 design")
    source = run_graph.load_source(run["id"], owner_id=data["admin"]["id"], project_id=data["project"]["id"])
    node = new_node("ipsae", "analysis")
    batches = run_graph.initialize_batches(data["root"] / "graph", node, source, node["data"]["params"], "abc", candidate_execution.policy())
    assert len(batches) == 2
    manifests = [json.loads(Path(b["manifest"]).read_text()) for b in batches]
    assert all(m["request"]["workflow"] == "analysis.ipsae" for m in manifests)
    assert {m["request"]["parameters"]["source_run_dir"] for m in manifests} == {str(p) for p in sources_paths(data)}


def sources_paths(data):
    from mn_molecular_design_studio.runtime import runs_root
    return [runs_root() / "design" / ("api-" + tool) for tool in ("genie3", "pxdesign")]


def test_graph_controller_consumes_complete_run_and_combines_sources(api_data, monkeypatch):
    from mn_molecular_design_studio import logical_runs, pipeline_graph
    from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from test_pipeline_graph import node, edge
    data = api_data
    sources(data)
    project = data["project"]["id"]
    run = next(r for r in logical_runs.catalog(data["admin"], project_id=project) if r["name"] == "R3 design")
    source = node("run_source", "input")
    source["data"]["params"]["run_reference"] = run["id"]
    analysis, result = node("ipsae", "analysis"), node("results", "results")
    graph = dict(schema_version=pipeline_graph.SCHEMA, id="b" * 32, name="Whole Run analysis",
                 nodes=[source, analysis, result], edges=[edge("input", "analysis"), edge("analysis", "results")])
    calls = []
    def submit(request, parent, n, **kwargs):
        calls.append({**request, "_gpu_device": kwargs.get("gpu_device")})
        job = create_job("analysis", "analysis", "ipsae", request["parameters"], {})
        rows = runner._absolute_candidates({"run_dir": request["parameters"]["source_run_dir"],
                                           "candidates_jsonl": request["parameters"]["candidates_jsonl"]})
        # The native calculations are simulated, but collection identity,
        # controller dispatch, recovery records and fan-in are the real code.
        output = runner._local_source(job.run_dir, "calculated", rows)
        destination = job.run_dir / "artifacts/normalized_candidates"
        destination.mkdir(parents=True)
        (destination / "candidates.jsonl").write_text(Path(output["candidates_jsonl"]).read_text())
        finish_job(job.run_dir, success=True, result={})
        return job.run_dir
    monkeypatch.setattr(runner, "_submit", submit)
    parent = runner.run_graph_pipeline(graph, "Whole Run analysis", data["admin"]["id"], project, reuse_completed=False, gpu_device="0")
    state = json.loads((parent / "graph_execution.json").read_text())
    assert state["nodes"]["results"]["status"] == "completed", state
    assert len(calls) == 2
    assert all(call["_gpu_device"] == "none" for call in calls)
    rows = read_candidates(parent)
    assert len(rows) == 2 and len({r["candidate_id"] for r in rows}) == 2


def test_offline_source_is_explicit_and_prevents_complete_export(api_data, monkeypatch):
    from mn_molecular_design_studio import logical_runs
    data = api_data
    paths = sources(data)
    project = data["project"]["id"]
    catalog = logical_runs.catalog(data["admin"], project_id=project)
    run = next(r for r in catalog if r["name"] == "R3 design")
    remote = run["outputs"][0]["sources"][1]
    remote.update(node_id="00000000-0000-0000-0000-000000000099", cluster_job_id="remote-source", availability="offline")
    monkeypatch.setattr(logical_runs, "catalog", lambda *a, **kw: [run])
    endpoint = "/api/v1/runs/" + run["id"] + "/results"
    response = data["client"].get(endpoint, params={"project_id": project})
    assert response.status_code == 200
    payload = response.json()
    assert not payload["complete"] and len(payload["source_errors"]) == 1
    assert payload["candidate_count"] == 1
    assert data["client"].get(endpoint, params={"project_id": project, "export": "jsonl"}).status_code == 409


def test_run_input_picker_restores_saved_output_and_keeps_identity_on_refresh(api_data, monkeypatch):
    from mn_molecular_design_studio import run_inputs_ui
    from streamlit.testing.v1 import AppTest
    choices = [{"reference": "run:a~results", "label": "Run A", "row": {}},
               {"reference": "run:b~analysis", "label": "Run B analysis", "row": {}}]
    monkeypatch.setattr(run_inputs_ui, "options", lambda *a, **kw: choices)
    app = AppTest.from_string('from mn_molecular_design_studio.run_inputs_ui import select\n'
        'import streamlit as st\n'
        'chosen=select("project", "source", reference="run:b~analysis")\n'
        'st.text(chosen["reference"])').run()
    assert not app.exception and app.selectbox[0].value == "run:b~analysis"

    choices[1]["label"] = "Run B analysis · 120 candidates"
    app.run()
    assert not app.exception and app.selectbox[0].value == "run:b~analysis"


def test_authorized_inspection_can_read_hidden_native_batches(api_data):
    from mn_molecular_design_studio.cluster.agent import Agent
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    data = api_data
    path = candidate_run(data, tool="genie3")
    write_json(path / "metadata.json", {**read_json(path / "metadata.json"), "hidden": True, "parent_role": "native_batch"})
    response, error = Agent().execute({"id": "inspect-hidden", "job_id": "authorized-source", "kind": "inspect",
        "payload": {"local_reference": "design/" + path.name}, "owner_id": data["admin"]["id"], "project_id": data["project"]["id"]})
    assert error is None, error
    assert response["results"]["candidate_count"] == 1


def test_direct_analysis_preserves_source_ownership_and_groups_under_run(api_data):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _analysis_job
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, submission_tags
    data = api_data
    path = candidate_run(data, tool="genie3")
    job = _analysis_job(path, path / "artifacts/normalized_candidates/candidates.jsonl",
                        tool="ipsae", params={}, note="test", mode="internal")
    metadata = read_json(job.run_dir / "metadata.json")
    assert metadata["studio_owner_id"] == data["admin"]["id"]
    assert metadata["project_id"] == data["project"]["id"]
    assert metadata["campaign_id"] == path.name
    with submission_tags(studio_owner_id="explicit", project_id="explicit-project"):
        other = _analysis_job(path, path / "artifacts/normalized_candidates/candidates.jsonl",
                              tool="ipsae", params={}, note="test", mode="internal")
    assert read_json(other.run_dir / "metadata.json")["studio_owner_id"] == "explicit"
    from mn_molecular_design_studio.run_reconciliation import analysis_identity_patch
    original = read_json(job.run_dir / "metadata.json")
    (job.run_dir / "metadata.json").write_text(json.dumps({k:v for k,v in original.items() if k not in {"studio_owner_id", "project_id"}}))
    repair = analysis_identity_patch(job.run_dir)
    assert repair["studio_owner_id"] == data["admin"]["id"] and repair["project_id"] == data["project"]["id"]
    with pytest.raises(ValueError, match="ownership differs"):
        analysis_identity_patch(other.run_dir)
