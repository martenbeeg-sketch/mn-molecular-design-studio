import asyncio
import copy
import hashlib
import json
import os
from pathlib import Path

import httpx
import pytest
from fastapi.testclient import TestClient

from test_binderscope import write_complex


@pytest.fixture
def preserved_api(tmp_path, monkeypatch):
    for name, value in {
        "MN_MOLECULAR_DESIGN_STUDIO_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_RESULTS_PARENT": tmp_path / "results",
        "MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR": tmp_path / "results/runs",
        "MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR": tmp_path / "references",
        "MDS_CLUSTER_HOME": tmp_path / "cluster",
        "MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE": "production",
    }.items():
        monkeypatch.setenv(name, str(value))
    from mn_molecular_design_studio import api, auth, molecular_database_api, workspace
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    workspace.save_settings({"results_directory": str(tmp_path / "results/runs"),
                             "reference_directory": str(tmp_path / "references")})
    workspace.ensure_studio_workspace()
    config.save({"mode": "disabled"})
    monkeypatch.setattr(molecular_database_api, "sync_workspace_async", lambda *args, **kwargs: [])
    actor = auth.create_initial_admin("research_instructor", "Instructor", "temporary-password-42", preserve_development_owner=True)
    auth.change_password(actor["id"], "instructor-password-42", forced=True)
    actor = next(row for row in auth.list_users() if row["id"] == actor["id"])
    assert actor["id"] == "00000000000000000000000000000001"
    project = workspace.create_project("Saved IDR designs", owner_id=actor["id"])["id"]
    other = auth.create_user("otherstudent", "Other", "temporary-password-42", actor=actor)
    auth.change_password(other["id"], "student-password-42", forced=True)
    source = tmp_path / "accepted.pdb"
    write_complex(source)
    database = MolecularDatabase()
    imported = database.import_records(actor, [{"name": "Accepted fixture", "components": [
        {"type": "protein", "chain": "X", "role": "target", "sequence": "A" * 16},
        {"type": "protein", "chain": "Y", "role": "binder", "sequence": "S" * 16}],
        "files": [{"name": "accepted.pdb", "content": source.read_bytes(), "binder_chains": ["Y"], "target_chains": ["X"]}]}],
        name="Accepted designs", project_id=project, import_method="tool_output",
        origin={"collection_kind": "generated", "task_group": "design", "tool": "bindcraft2", "run_id": "fixture-design"})
    record = database.list_molecules(actor, molecule_ids=imported["molecule_ids"], collection="all")["items"][0]
    structure = record["representations"][0]
    payload = {"schema_version": 1, "workflow": "analysis.binderscope", "project_id": project,
        "parameters": {"cpu_cores": 32, "use_recorded_settings": True}, "resources": {"gpu_device": "none"},
        "binderscope_collection": {"schema": "studio.binderscope_selection.v1", "sources": [{
            "kind": "database", "molecule_id": record["id"], "group_id": imported["group_id"],
            "structure_id": structure["id"], "structure_digest": hashlib.sha256(source.read_bytes()).hexdigest(),
            "binder_chains": ["Y"], "target_chains": ["X"], "assessment": {
                "component_format": "linear protein", "b_factor_interpretation": "plddt_0_100"}}]}}
    pending = []
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda run: pending.append(Path(run)))
    with TestClient(api.app) as client:
        login = client.post("/api/v1/auth/login", json={"username": "research_instructor", "password": "instructor-password-42"})
        assert login.status_code == 200
        token = login.json()["access_token"]
        other_login = client.post("/api/v1/auth/login", json={"username": "otherstudent", "password": "student-password-42"})
        assert other_login.status_code == 200
        studio = StudioClient("http://studio.test", token=token, transport=httpx.ASGITransport(app=api.app))
        yield {"client": client, "studio": studio, "actor": actor, "database": database, "payload": payload,
               "pending": pending, "token": token, "other_token": other_login.json()["access_token"], "root": tmp_path}


def test_authenticated_preserved_owner_queues_505_selectors_once(preserved_api, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_sources
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

    case = preserved_api
    payload = copy.deepcopy(case["payload"])
    selector = payload["binderscope_collection"]["sources"][0]
    with case["database"].connection() as connection:
        for position in range(1, 505):
            identifier = "http-cohort-" + str(position)
            connection.execute("INSERT INTO molecule_structures SELECT ?,molecule_id,group_id,path,label,format,context,created_at "
                               "FROM molecule_structures WHERE id=?", (identifier, selector["structure_id"]))
            payload["binderscope_collection"]["sources"].append({**selector, "structure_id": identifier})
    monkeypatch.setattr(complex_analysis_sources, "build_saved_structure_collection",
                        lambda *args, **kwargs: pytest.fail("HTTP must leave coordinate preparation to the worker"))

    async def submit():
        identity = await case["studio"].request("GET", "/api/v1/auth/me")
        assert identity["id"] == case["actor"]["id"] and identity["development_mode"] is False
        validation = await case["studio"].request("POST", "/api/v1/workflows/validate", body=payload)
        assert validation["designs"] == 505 and validation["jobs"] == 1
        return await case["studio"].request("POST", "/api/v1/workflows", body=payload, headers={"Idempotency-Key": "http-cohort"})

    submitted = asyncio.run(submit())
    assert len(case["pending"]) == 1
    run = case["pending"][0]
    assert run.is_relative_to(case["root"])
    request = read_json(run / "worker_request.json")
    assert request["function"] == "run_saved_binderscope" and request["resources"]["cpu_cores"] == 32
    assert request["resources"]["gpu_device"] == "none"
    assert read_json(run / "metadata.json")["studio_owner_id"] == case["actor"]["id"]
    assert read_json(run / "metadata.json")["studio_project_id"] == payload["project_id"]
    headers = {"Authorization": "Bearer " + case["token"], "Idempotency-Key": "http-cohort"}
    retry = case["client"].post("/api/v1/workflows", json=payload, headers=headers)
    assert retry.status_code == 202 and retry.json() == submitted
    changed = copy.deepcopy(payload)
    changed["parameters"]["cpu_cores"] = 4
    assert case["client"].post("/api/v1/workflows", json=changed, headers=headers).status_code == 409
    assert len(case["pending"]) == 1


def test_binderscope_http_rejects_unauthenticated_and_inaccessible_selectors(preserved_api):
    case = preserved_api
    for route in ("/api/v1/workflows/validate", "/api/v1/workflows"):
        assert case["client"].post(route, json=case["payload"]).status_code == 401
        assert case["client"].post(route, json=case["payload"], headers={"Authorization": "Bearer invalid"}).status_code == 401
        assert case["client"].post(route, json=case["payload"], headers={"Authorization": "Bearer " + case["other_token"]}).status_code == 404
    changed = copy.deepcopy(case["payload"])
    changed["binderscope_collection"]["sources"][0]["structure_id"] = "missing"
    headers = {"Authorization": "Bearer " + case["token"]}
    assert case["client"].post("/api/v1/workflows", json=changed, headers=headers).status_code == 403
    assert not case["pending"]


def test_mcp_http_accepts_preserved_production_instructor(preserved_api):
    pytest.importorskip("mcp")
    from mn_molecular_design_studio.mcp_server import http_app

    case = preserved_api
    with TestClient(http_app(client=case["studio"]), base_url="http://localhost") as client:
        response = client.post("/mcp", headers={"Authorization": "Bearer " + case["token"],
            "Accept": "application/json, text/event-stream", "MCP-Protocol-Version": "2025-11-25"},
            json={"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": "studio_identity", "arguments": {}}})
        assert response.status_code == 200, response.text
        messages = [json.loads(line[6:]) for line in response.text.splitlines() if line.startswith("data: ")]
        assert messages
        result = messages[-1]["result"]
        assert not result.get("isError")
        assert result["structuredContent"]["id"] == case["actor"]["id"]
        assert result["structuredContent"]["development_mode"] is False


@pytest.mark.skipif(os.environ.get("MDS_TEST_BINDERSCOPE_IMAGE") != "1", reason="Explicit offline authenticated HTTP integration")
def test_authenticated_http_worker_and_exact_sql_attachment(preserved_api):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import run_worker_job
    from mn_molecular_design_studio.molecular_database_ingest import record_completed_run

    case = preserved_api
    payload = copy.deepcopy(case["payload"])
    payload["parameters"]["cpu_cores"] = 1
    selector = payload["binderscope_collection"]["sources"][0]
    with case["database"].connection() as connection:
        original_members = [tuple(row) for row in connection.execute("SELECT * FROM molecule_members")]
        original_structures = [tuple(row) for row in connection.execute("SELECT * FROM molecule_structures")]

    async def submit():
        await case["studio"].request("POST", "/api/v1/workflows/validate", body=payload)
        return await case["studio"].request("POST", "/api/v1/workflows", body=payload, headers={"Idempotency-Key": "real-http-worker"})

    submitted = asyncio.run(submit())
    assert len(case["pending"]) == 1
    run = case["pending"][0]
    assert run_worker_job(run) == 0
    assert read_json(run / "result.json")["success"] is True
    record_completed_run(run)
    record_completed_run(run)
    with case["database"].connection() as connection:
        rows = connection.execute("SELECT * FROM molecule_observations WHERE namespace='binderscope' AND run_id=?", (run.name,)).fetchall()
        assert len(rows) == 1
        observation = dict(rows[0])
        context = json.loads(observation["context"])
        assert observation["molecule_id"] == selector["molecule_id"]
        assert context["representation_id"] == selector["structure_id"]
        assert context["source_group_id"] == selector["group_id"]
        assert context["source_structure_digest"] == selector["structure_digest"]
        assert [tuple(row) for row in connection.execute("SELECT * FROM molecule_structures")] == original_structures
        members = [tuple(row) for row in connection.execute("SELECT * FROM molecule_members WHERE group_id=?", (selector["group_id"],))]
        assert members == original_members
    retry = case["client"].post("/api/v1/workflows", json=payload, headers={"Authorization": "Bearer " + case["token"], "Idempotency-Key": "real-http-worker"})
    assert retry.status_code == 202 and retry.json() == submitted and len(case["pending"]) == 1
