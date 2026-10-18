from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mn_molecular_design_studio.cluster import capabilities, config, service, transport
from mn_molecular_design_studio.cluster.agent import Agent
from mn_molecular_design_studio.cluster.store import Store


@pytest.fixture
def api_cluster(monkeypatch, tmp_path):
    for name, value in {
        "MDS_CLUSTER_HOME": tmp_path / "cluster", "MN_MOLECULAR_DESIGN_STUDIO_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_RESULTS_PARENT": tmp_path / "results", "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR": tmp_path / "results/runs", "MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR": tmp_path / "references",
    }.items():
        monkeypatch.setenv(name, str(value))
    from mn_molecular_design_studio import workspace, api
    workspace.ensure_studio_workspace()
    project = workspace.create_project("Shared project", owner_id="owner")
    config.save({"mode": "coordinator", "coordinator_url": "http://coordinator.test:8766", "default_placement": "auto"})
    store = Store()
    origin = store.enroll("current", config.identity())
    remote = store.enroll("4032")
    actor = {"id": "owner", "role": "admin", "must_change_password": False}
    snapshot = {"images": {}, "references": {}, "code_revision": "test-code", "accept_jobs": True, "scheduler_online": True,
                "resources": {"free_gpus": [], "cpu_free": 16, "ram_available_gb": 64, "scratch_free_gb": 100}}
    needs = {"images": {}, "references": [], "code_revision": "test-code", "gpu_count": 0, "cpu_threads": 1}
    store.heartbeat(origin["node_id"], snapshot, [])
    store.heartbeat(remote["node_id"], snapshot, [])
    monkeypatch.setattr(capabilities, "requirements", lambda *a, **k: needs)
    api.app.dependency_overrides[api._current_user] = lambda: actor
    try:
        with TestClient(api.app) as client:
            yield client, store, origin, remote, actor, project, snapshot, needs
    finally:
        api.app.dependency_overrides.pop(api._current_user, None)


def test_api_to_remote_agent_to_local_scheduler_preserves_source_and_ownership(api_cluster, monkeypatch, tmp_path):
    client, store, origin, remote, actor, project, snapshot, needs = api_cluster
    source = {"task_group": "design", "run_id": "source", "owner_id": "owner", "project_id": project["id"], "status": "completed"}
    store.heartbeat(remote["node_id"], snapshot, [source])
    source_id = config.job_id(remote["node_id"], "design/source")
    payload = {"schema_version": 1, "project_id": project["id"], "workflow": "analysis.ipsae", "parameters": {"pae_cutoff": 10.0},
               "source_job": "cluster:" + source_id, "execution": {"placement": "auto", "idempotency_key": "click"}}
    # Even though the origin has capacity, analysis must follow the source.
    validated = client.post("/api/v1/workflows/validate", json=payload)
    assert validated.status_code == 200 and validated.json()["file_validation"]
    response = client.post("/api/v1/workflows", json=payload)
    assert response.status_code == 202, response.text
    job = response.json()["job"]
    assert job["node_id"] == remote["node_id"]
    assert client.post("/api/v1/workflows", json=payload).json()["job"]["cluster_job_id"] == job["cluster_job_id"]
    assert not list((tmp_path / "results").rglob("metadata.json"))
    headers = {"Authorization": "Bearer " + remote["node_token"]}
    command = client.post("/api/v1/cluster/node/claim", headers=headers).json()["command"]
    assert command["job_id"] == job["cluster_job_id"]
    assert command["payload"]["payload"]["source_job"].startswith("cluster://" + remote["node_id"])

    # Simulate only the scientific adapter; exercise the real agent, receipt,
    # ownership marker, local run creation, and coordinator acknowledgement.
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    remote_root = tmp_path / "remote-runs"; remote_root.mkdir()
    monkeypatch.setattr(jobs, "runs_root", lambda: remote_root)
    monkeypatch.setattr(config, "identity", lambda: remote["node_id"])
    captured = []
    monkeypatch.setattr(workflow_automation, "validate_workflow_request", lambda path: {"valid": True})
    def submit(path):
        request = json.loads(Path(path).read_text()); captured.append(request)
        run = jobs.create_job("analysis", "analysis", "ipsae", {}, {})
        return {"job": {"run_dir": str(run.run_dir), "task_group": "analysis", "run_id": run.run_id, "status": "queued", "job_code": "LOCAL"}}
    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit)
    agent = Agent(); monkeypatch.setattr(agent, "snapshot", lambda: snapshot)
    result, error = agent.execute(command)
    assert error is None, error
    assert captured[0]["source_job"] == "design/source"
    local_run = Path(result["job"]["run_dir"])
    metadata = json.loads((local_run / "metadata.json").read_text())
    assert metadata["project_id"] == project["id"] and metadata["studio_owner_id"] == "owner"
    assert metadata["execution_node_id"] == remote["node_id"]
    assert metadata["cluster_job_id"] == job["cluster_job_id"]
    acknowledgement = client.post("/api/v1/cluster/node/commands/" + command["id"] + "/finish", headers=headers,
                                  json={"lease_token": command["lease_token"], "response": result})
    assert acknowledgement.status_code == 200, acknowledgement.text
    shared = client.get("/api/v1/cluster/jobs/" + job["cluster_job_id"]).json()
    assert shared["node_id"] == remote["node_id"] and shared["run_dir"].startswith("cluster://")
    assert shared["location"] == str(local_run)
    agent.execute(command)
    assert len(captured) == 1


def test_remote_source_cannot_be_forced_onto_a_different_workstation(api_cluster):
    client, store, origin, remote, actor, project, snapshot, _ = api_cluster
    store.heartbeat(remote["node_id"], snapshot, [{"task_group": "design", "run_id": "remote", "owner_id": "owner", "project_id": project["id"], "status": "completed"}])
    payload = {"schema_version": 1, "workflow": "analysis.ipsae", "project_id": project["id"], "source_job": "cluster:" + config.job_id(remote["node_id"], "design/remote"),
               "execution": {"placement": origin["node_id"]}}
    result = client.post("/api/v1/workflows", json=payload)
    assert result.status_code == 422 and "owns" in result.text


def test_local_job_code_becomes_node_qualified_before_automatic_dispatch(api_cluster, monkeypatch):
    client, store, origin, _, actor, project, snapshot, _ = api_cluster
    from mn_molecular_design_studio import job_catalog
    from mn_molecular_design_studio.engine.workflow_engine.services import local_automation
    source = {"task_group": "design", "run_id": "local-source", "owner_id": actor["id"], "project_id": project["id"], "status": "completed"}
    monkeypatch.setattr(local_automation, "resolve_job", lambda ref: source)
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda *a: source)
    store.heartbeat(origin["node_id"], snapshot, [source])
    response = client.post("/api/v1/workflows", json={"schema_version": 1, "workflow": "analysis.ipsae", "project_id": project["id"], "source_job": "LOCAL-CODE"})
    assert response.status_code == 202, response.text
    assert response.json()["job"]["node_id"] == origin["node_id"]
    command = store.claim(origin["node_id"])
    assert command["payload"]["payload"]["source_job"] == config.artifact_uri(origin["node_id"], "design", "local-source")


def test_uploaded_input_reaches_only_the_assigned_agent(api_cluster, monkeypatch, tmp_path):
    client, store, origin, remote, _, project, snapshot, _ = api_cluster
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    structure = tmp_path / "uploaded.pdb"; structure.write_text("TEST INPUT COORDINATES")
    response = client.post("/api/v1/workflows", json={"schema_version": 1, "workflow": "sequence.esm_if1", "project_id": project["id"],
        "parameters": {"structure_path": str(structure)}, "execution": {"placement": remote["node_id"]}})
    assert response.status_code == 202, response.text
    command = store.claim(remote["node_id"])
    digest = command["payload"]["bundle_id"]
    assert digest and command["payload"]["payload"]["parameters"]["structure_path"].startswith("cluster-input:///")
    assert client.get("/api/v1/cluster/node/blobs/" + digest, headers={"Authorization": "Bearer " + origin["node_token"]}).status_code == 404
    assert client.get("/api/v1/cluster/node/blobs/" + digest, headers={"Authorization": "Bearer " + remote["node_token"]}).status_code == 200
    remote_root = tmp_path / "remote-results"; remote_root.mkdir()
    monkeypatch.setattr(jobs, "runs_root", lambda: remote_root)
    monkeypatch.setattr(config, "identity", lambda: remote["node_id"])
    captured = []
    def validate(path):
        payload = json.loads(Path(path).read_text())
        staged = Path(payload["parameters"]["structure_path"])
        assert staged != structure and staged.read_text() == structure.read_text()
        captured.append(staged)
        return {"valid": True}
    def submit(path):
        job = jobs.create_job("sequence", "test", "esm_if1", {}, {})
        return {"job": {"run_dir": str(job.run_dir), "task_group": "sequence", "run_id": job.run_id, "status": "queued"}}
    monkeypatch.setattr(workflow_automation, "validate_workflow_request", validate)
    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit)
    agent = Agent(); monkeypatch.setattr(agent, "snapshot", lambda: snapshot)
    result, error = agent.execute(command)
    assert error is None, error
    assert len(captured) == 1 and result["job"]["node_id"] == remote["node_id"]


@pytest.mark.parametrize("path", ["/private/data", "../../private/data", "local/input.pdb", "local.pdb"])
def test_direct_cluster_endpoint_rejects_unstaged_machine_paths(api_cluster, path):
    client, _, origin, _, _, project, _, needs = api_cluster
    response = client.post("/api/v1/cluster/submissions", json={"kind": "workflow", "project_id": project["id"], "origin_node": origin["node_id"], "requirements": needs,
        "payload": {"schema_version": 1, "workflow": "analysis.ipsae", "parameters": {"source_run_dir": path, "candidates_jsonl": "/private/data/candidates.jsonl"}}})
    assert response.status_code == 422 and "machine-qualified" in response.text


def test_input_and_artifact_requests_are_bound_to_job_access(api_cluster):
    client, store, origin, remote, actor, project, snapshot, _ = api_cluster
    store.heartbeat(remote["node_id"], snapshot, [{"task_group": "design", "run_id": "private", "owner_id": "owner", "project_id": project["id"], "status": "completed"}])
    jid = config.job_id(remote["node_id"], "design/private")
    operation = client.post(f"/api/v1/cluster/jobs/{jid}/actions", json={"action": "artifact", "parameters": {"path": "result.json"}})
    assert operation.status_code == 202
    from mn_molecular_design_studio import api
    api.app.dependency_overrides[api._current_user] = lambda: {"id": "outsider", "role": "user", "project_ids": [], "_cluster_project_ids": []}
    assert client.get(f"/api/v1/cluster/jobs/{jid}").status_code == 404
    assert client.get("/api/v1/cluster/commands/" + operation.json()["command_id"]).status_code == 404
    assert client.post(f"/api/v1/cluster/jobs/{jid}/actions", json={"action": "inspect"}).status_code == 404


def test_claiming_coordinator_origin_does_not_bypass_input_bundle_ownership(api_cluster):
    client, store, origin, _, _, project, _, needs = api_cluster
    from mn_molecular_design_studio.cluster import bundles
    digest = "a" * 64
    bundles.blob_path(digest).write_bytes(b"private input bundle")
    with store.connect() as db:
        db.execute("CREATE TABLE IF NOT EXISTS blob_grants(digest TEXT,owner_id TEXT,PRIMARY KEY(digest,owner_id))")
        db.execute("INSERT INTO blob_grants VALUES (?,?)", (digest, "another-user"))
    response = client.post("/api/v1/cluster/submissions", json={"kind": "workflow", "project_id": project["id"], "origin_node": origin["node_id"], "requirements": needs,
        "bundle_id": digest, "payload": {"schema_version": 1, "workflow": "sequence.esm_if1", "parameters": {"structure_path": "cluster-input:///inputs/target.pdb"}}})
    assert response.status_code == 403 and "does not belong" in response.text


def test_lost_storage_mount_does_not_create_an_empty_replacement(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    from mn_molecular_design_studio import runtime
    root = tmp_path / "runs"; root.mkdir()
    monkeypatch.setattr(runtime, "runs_root", lambda: root)
    agent = Agent()
    assert agent.storage_online()
    (root / ".studio-workstation").unlink()
    assert not agent.storage_online()
    assert not (root / ".studio-workstation").exists()


def test_agent_never_reads_artifact_outside_the_selected_run(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    run = tmp_path / "run"; run.mkdir()
    outside = tmp_path / "secret"; outside.write_text("private")
    (run / "symlink").symlink_to(outside)
    from mn_molecular_design_studio.engine.workflow_engine.services import local_automation
    monkeypatch.setattr(local_automation, "resolve_job", lambda ref: {"run_dir": str(run)})
    for path in ("../secret", "symlink", str(outside)):
        with pytest.raises(ValueError):
            Agent().artifact({"payload": {"local_reference": "design/run", "path": path}})


def test_worker_session_uses_shared_user_id_and_offline_cache_until_original_expiry(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    from datetime import datetime, timezone
    from mn_molecular_design_studio.cluster import identity_auth
    user = {"id": "global-user", "username": "scientist", "role": "user", "project_ids": ["global-project"],
            "session_expires_at": datetime.fromtimestamp(time.time() + 600, timezone.utc).isoformat()}
    monkeypatch.setattr(transport, "request", lambda *a, **k: user.copy())
    assert identity_auth.resolve("central-token")["_cluster_project_ids"] == ["global-project"]
    def offline(*a, **k): raise transport.ClusterUnavailable("offline")
    monkeypatch.setattr(transport, "request", offline)
    cached = identity_auth.resolve("central-token")
    assert cached["id"] == "global-user" and cached["coordinator_offline"]
    assert identity_auth.resolve("different-token") is None
    def revoked(*a, **k): raise ValueError("invalid session")
    monkeypatch.setattr(transport, "request", revoked)
    assert identity_auth.resolve("central-token") is None
    monkeypatch.setattr(transport, "request", offline)
    assert identity_auth.resolve("central-token") is None
