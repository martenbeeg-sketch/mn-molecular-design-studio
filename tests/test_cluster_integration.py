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


def test_candidate_batch_poll_persists_pending_candidate_artifact_transfer(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import candidate_batches

    remote_root = "/worker/runs/design/batch"
    structure_rel = "artifacts/rfd4_proteina/output/complex.cif"
    collection = "artifacts/distributed_campaign/accepted_candidates.jsonl"
    artifact_names = [
        "metadata.json", "input.json", "result.json",
        collection, structure_rel,
    ]

    class FakeStore:
        def job(self, _job_id, _user):
            return {
                "status": "completed", "location": remote_root,
                "artifacts": [{"relative_path": f"{remote_root}/{name}"} for name in artifact_names],
                "candidate_summary": {
                    "status": "complete", "candidate_count": 50, "eligible_count": 50,
                    "native_filter_counts": {"passed": 50, "failed": 0, "unavailable": 0},
                    "acceptance_counts": {"accepted": 1, "rejected": 49, "unresolved": 0, "not_evaluated": 0},
                },
            }

    monkeypatch.setattr(candidate_batches, "Store", FakeStore)
    monkeypatch.setattr(candidate_batches, "actor", lambda *_args, **_kwargs: {"id": "owner"})
    peer = {
        "cluster_job_id": "cluster-job",
        "manifest": str(tmp_path / "peer.json"),
        "transfers": {},
    }
    artifact_requests = []
    fetched_artifacts = []

    def artifact(state, _job_id, relative, _user, directory, *, persist_transfers=None):
        fetched_artifacts.append(relative)
        path = Path(directory) / relative
        if relative == structure_rel:
            transfer = state.setdefault("transfers", {}).setdefault(relative, {})
            if "command_id" not in transfer:
                transfer["command_id"] = "structure-transfer-1"
                artifact_requests.append(relative)
                if persist_transfers is not None:
                    persist_transfers(state["transfers"])
                return None
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("verified complex\n", encoding="utf-8")
            transfer.update(sha256="fixture", bytes=path.stat().st_size)
            return path

        path.parent.mkdir(parents=True, exist_ok=True)
        if relative == "metadata.json":
            path.write_text(json.dumps({"status": "completed"}), encoding="utf-8")
        elif relative == "input.json":
            path.write_text("{}", encoding="utf-8")
        elif relative == "result.json":
            path.write_text(json.dumps({"success": True, "metrics": {}}), encoding="utf-8")
        elif relative == collection:
            path.write_text(json.dumps({"candidate_id": "candidate-1", "complex_pdb": structure_rel}) + "\n", encoding="utf-8")
        return path

    monkeypatch.setattr(candidate_batches, "_artifact", artifact)
    persisted = []
    status, output, _job = candidate_batches.poll(
        peer, owner_id="owner", project_id="project",
        persist_transfers=lambda transfers: persisted.append(dict(transfers)),
        candidate_collection=collection,
        require_result=False,
    )
    assert (status, output) == ("collecting", None)
    assert artifact_requests == [structure_rel]
    assert persisted[-1][structure_rel]["command_id"] == "structure-transfer-1"

    status, output, _job = candidate_batches.poll(
        peer, owner_id="owner", project_id="project",
        persist_transfers=lambda transfers: persisted.append(dict(transfers)),
        candidate_collection=collection,
        require_result=False,
    )
    assert status == "completed"
    assert output["candidate_count"] == 1
    assert Path(output["candidates_jsonl"]).read_text(encoding="utf-8").find("candidate-1") >= 0
    assert artifact_requests == [structure_rel]
    assert "result.json" not in fetched_artifacts
    assert "native_execution.json" not in fetched_artifacts
    assert output["source_provenance"]["candidate_summary"]["acceptance_counts"]["accepted"] == 1


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
    owner_id = uuid.uuid4().hex
    project = workspace.create_project("Shared project", owner_id=owner_id)
    config.save({"mode": "coordinator", "coordinator_url": "http://coordinator.test:8766", "default_placement": "auto"})
    store = Store()
    origin = store.enroll("current", config.identity())
    remote = store.enroll("4032")
    actor = {"id": owner_id, "role": "admin", "must_change_password": False}
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
    source = {"task_group": "design", "run_id": "source", "owner_id": actor["id"], "project_id": project["id"], "status": "completed"}
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
    assert metadata["project_id"] == project["id"] and metadata["studio_owner_id"] == actor["id"]
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


def test_api_purge_removes_terminal_job_database_and_local_run_files(api_cluster, tmp_path):
    client, store, origin, _remote, actor, project, snapshot, _needs = api_cluster
    run_id = "purge-test-run"
    task_group = "ligand-docking"
    from mn_molecular_design_studio.workspace import ensure_studio_workspace
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import runs_root
    ensure_studio_workspace()
    run_dir = runs_root() / task_group / run_id
    run_dir.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps({"status": "completed", "run_id": run_id}))
    jid = str(uuid.uuid4())
    store.heartbeat(origin["node_id"], snapshot, [{
        "cluster_job_id": jid,
        "task_group": task_group,
        "run_id": run_id,
        "job_type": "ligand_docking",
        "tool": "vina",
        "status": "completed",
        "owner_id": actor["id"],
        "project_id": project["id"],
        "run_dir": str(run_dir),
    }])

    response = client.delete(f"/api/v1/cluster/jobs/{jid}/purge")

    assert response.status_code == 200, response.text
    result = response.json()
    assert result["status"] == "purged"
    assert result["run_files_removed"] is True
    assert result["database"]["preserved_saved_selections"] is True
    assert not run_dir.exists()
    assert store.search({"id": actor["id"], "role": "admin", "project_ids": [project["id"]]}, project_id=project["id"])["total"] == 0
    assert client.get(f"/api/v1/cluster/jobs/{jid}").status_code == 404


def test_api_purge_rejects_active_jobs_and_queues_remote_run_cleanup(api_cluster):
    client, store, origin, remote, actor, project, snapshot, _needs = api_cluster
    local_active_id = str(uuid.uuid4())
    store.heartbeat(origin["node_id"], snapshot, [{
        "cluster_job_id": local_active_id, "task_group": "design", "run_id": "still-running",
        "status": "running", "owner_id": actor["id"], "project_id": project["id"],
    }])
    active = client.delete(f"/api/v1/cluster/jobs/{local_active_id}/purge")
    assert active.status_code == 409
    assert client.get(f"/api/v1/cluster/jobs/{local_active_id}").status_code == 200

    # The coordinator catalog can lag a workstation's local metadata. Never
    # delete a present run directory if it still identifies itself as active.
    from mn_molecular_design_studio.workspace import ensure_studio_workspace
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import runs_root
    ensure_studio_workspace()
    mismatched_id = str(uuid.uuid4())
    mismatched_dir = runs_root() / "design" / "catalog-says-failed"
    mismatched_dir.mkdir(parents=True)
    (mismatched_dir / "metadata.json").write_text(json.dumps({"status": "running"}))
    store.heartbeat(origin["node_id"], snapshot, [{
        "cluster_job_id": mismatched_id, "task_group": "design", "run_id": "catalog-says-failed",
        "status": "failed", "owner_id": actor["id"], "project_id": project["id"],
    }])
    mismatch_response = client.delete(f"/api/v1/cluster/jobs/{mismatched_id}/purge")
    assert mismatch_response.status_code == 409
    assert mismatched_dir.exists()
    import shutil
    shutil.rmtree(mismatched_dir)

    remote_id = str(uuid.uuid4())
    store.heartbeat(remote["node_id"], snapshot, [{
        "cluster_job_id": remote_id, "task_group": "design", "run_id": "remote-failed",
        "status": "failed", "owner_id": actor["id"], "project_id": project["id"],
    }])
    remote_job = store.job(remote_id, {"id": actor["id"], "role": "admin", "project_ids": [project["id"]]})
    assert remote_job["local_reference"] == "design/remote-failed"
    remote_response = client.delete(f"/api/v1/cluster/jobs/{remote_id}/purge")
    assert remote_response.status_code == 202
    operation = remote_response.json()
    assert operation["status"] == "purge_queued"
    command = store.claim(remote["node_id"])
    assert command["id"] == operation["command_id"]
    assert command["kind"] == "purge"
    store.finish(remote["node_id"], command["id"], command["lease_token"], response={"status": "purged"})
    assert client.get(f"/api/v1/cluster/jobs/{remote_id}").status_code == 404
    assert client.get(f"/api/v1/cluster/commands/{command['id']}").json()["status"] == "completed"


def test_api_purge_removes_paused_orphan_when_workstation_confirms_run_missing(api_cluster):
    client, store, origin, remote, actor, project, snapshot, _needs = api_cluster
    from mn_molecular_design_studio.workspace import ensure_studio_workspace
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import runs_root
    ensure_studio_workspace()

    # A stale missing-reference report must never authorize deletion of a
    # directory whose own metadata still says it is paused.
    protected_id = str(uuid.uuid4())
    protected_dir = runs_root() / "analysis" / "paused-still-present"
    protected_dir.mkdir(parents=True)
    (protected_dir / "metadata.json").write_text(json.dumps({"status": "paused"}))
    store.heartbeat(origin["node_id"], snapshot, [{
        "cluster_job_id": protected_id, "task_group": "analysis", "run_id": "paused-still-present",
        "status": "paused", "owner_id": actor["id"], "project_id": project["id"],
    }])
    store.heartbeat(origin["node_id"], snapshot, [], known_references=[])
    protected_response = client.delete(f"/api/v1/cluster/jobs/{protected_id}/purge")
    assert protected_response.status_code == 409
    assert protected_dir.exists()
    import shutil
    shutil.rmtree(protected_dir)
    assert client.delete(f"/api/v1/cluster/jobs/{protected_id}/purge").status_code == 200

    local_id = str(uuid.uuid4())
    store.heartbeat(origin["node_id"], snapshot, [{
        "cluster_job_id": local_id, "task_group": "analysis", "run_id": "paused-orphan-local",
        "status": "paused",
        "owner_id": actor["id"], "project_id": project["id"],
    }])
    store.heartbeat(origin["node_id"], snapshot, [], known_references=[])
    local_response = client.delete(f"/api/v1/cluster/jobs/{local_id}/purge")
    assert local_response.status_code == 200, local_response.text
    assert local_response.json()["run_files_removed"] is False
    assert client.get(f"/api/v1/cluster/jobs/{local_id}").status_code == 404

    remote_id = str(uuid.uuid4())
    store.heartbeat(remote["node_id"], snapshot, [{
        "cluster_job_id": remote_id, "task_group": "analysis", "run_id": "paused-orphan-remote",
        "status": "paused",
        "owner_id": actor["id"], "project_id": project["id"],
    }])
    store.heartbeat(remote["node_id"], snapshot, [], known_references=[])
    cancel = store.action(remote_id, "cancel", {}, actor)
    assert cancel["status"] == "queued"

    remote_response = client.delete(f"/api/v1/cluster/jobs/{remote_id}/purge")
    assert remote_response.status_code == 202, remote_response.text
    with store.connect() as db:
        cancel_row = db.execute("SELECT status FROM commands WHERE id=?", (cancel["command_id"],)).fetchone()
    assert cancel_row["status"] == "cancelled"

    command = store.claim(remote["node_id"])
    assert command["kind"] == "purge"
    agent = Agent()
    result, error = agent.execute(command)
    assert error is None, error
    assert result["status"] == "purged"
    store.finish(remote["node_id"], command["id"], command["lease_token"], response=result)
    assert client.get(f"/api/v1/cluster/jobs/{remote_id}").status_code == 404


def test_placement_options_api_previews_candidates_without_creating_a_job(api_cluster, monkeypatch):
    client, store, _origin, _remote, _actor, project, _snapshot, _needs = api_cluster
    monkeypatch.setattr(service, "schema_check", lambda *_args: None)
    response = client.post("/api/v1/workflows/placement-options", json={
        "schema_version": 1, "workflow": "analysis.ipsae", "parameters": {"pae_cutoff": 10.0},
        "project_id": project["id"], "execution": {"placement": "auto", "idempotency_key": "preview"},
    })

    assert response.status_code == 200, response.text
    assert len(response.json()["candidates"]) == 2
    with store.connect() as db:
        assert db.execute("SELECT count(*) FROM submissions").fetchone()[0] == 0


def test_cluster_preflight_enforces_local_input_access_before_schema_validation(api_cluster, monkeypatch, tmp_path):
    client, _store, _origin, _remote, _actor, project, _snapshot, _needs = api_cluster
    from mn_molecular_design_studio import api

    user = {"id": "1234567890abcdef1234567890abcdef", "role": "user"}
    upload_root = tmp_path / "users" / user["id"] / "inputs"
    upload_root.mkdir(parents=True)
    monkeypatch.setattr(api, "user_input_root", lambda _user_id: upload_root)
    monkeypatch.setattr(api, "accessible_project_ids", lambda _user: [project["id"]])
    monkeypatch.setattr(api, "require_project_access", lambda *_args: None)
    monkeypatch.setattr(api, "path_is_within_accessible_job", lambda *_args, **_kwargs: False)
    api.app.dependency_overrides[api._current_user] = lambda: user
    try:
        response = client.post("/api/v1/workflows/validate", json={
            "schema_version": 1,
            "workflow": "detection.pesto",
            "project_id": project["id"],
            "parameters": {"target_pdb": "/etc/passwd", "chain_ids": ["A"]},
        })
    finally:
        api.app.dependency_overrides[api._current_user] = lambda: _actor

    assert response.status_code == 403
    assert "uploaded" in response.json()["detail"]


def test_remote_source_cannot_be_forced_onto_a_different_workstation(api_cluster):
    client, store, origin, remote, actor, project, snapshot, _ = api_cluster
    store.heartbeat(remote["node_id"], snapshot, [{"task_group": "design", "run_id": "remote", "owner_id": actor["id"], "project_id": project["id"], "status": "completed"}])
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
    store.heartbeat(remote["node_id"], snapshot, [{"task_group": "design", "run_id": "private", "owner_id": actor["id"], "project_id": project["id"], "status": "completed"}])
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
