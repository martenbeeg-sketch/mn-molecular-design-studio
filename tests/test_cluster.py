"""Workstation federation contracts. No models, Docker runs or SSH connections."""
from __future__ import annotations

import copy
import hashlib
import io
import json
import time
import uuid
import zipfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mn_molecular_design_studio.cluster import bundles, config, service
from mn_molecular_design_studio.cluster.agent import Agent
from mn_molecular_design_studio.cluster.store import Store


@pytest.fixture
def cluster(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    clock = [1000.0]
    store = Store(clock=lambda: clock[0])
    nodes = [store.enroll(name) for name in ("current", "4032", "4033")]
    user = {"id": "user-a", "role": "user", "project_ids": ["project-a"]}
    needs = {"images": {"engine:1": "sha256:exact"}, "references": ["weights/model.pt"], "code_revision": "same-code",
             "gpu_count": 1, "min_vram_gb": 8, "cpu_threads": 4, "ram_gb": 10, "scratch_gb": 1}
    snapshot = {"images": {"engine:1": "sha256:exact"}, "references": {"weights/model.pt": 100}, "code_revision": "same-code",
                "accept_jobs": True, "scheduler_online": True,
                "resources": {"free_gpus": [{"gpu_id": 0, "free_vram_gb": 24}], "cpu_free": 16,
                              "ram_available_gb": 64, "scratch_free_gb": 100}}
    for node in nodes:
        store.heartbeat(node["node_id"], snapshot, [])
    request = {"kind": "workflow", "payload": {"schema_version": 1, "workflow": "analysis.ipsae", "parameters": {}},
               "project_id": "project-a", "bundle_id": None}
    return store, nodes, user, needs, snapshot, request, clock


def test_identity_is_stable_and_node_ids_do_not_depend_on_storage_paths(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "a"))
    first = config.identity()
    assert config.identity() == first
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "b"))
    second = config.identity()
    assert first != second
    assert config.job_id(first, "design/same-run") != config.job_id(second, "design/same-run")


def test_automatic_placement_uses_other_ready_workstation_when_origin_is_busy(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    busy = copy.deepcopy(snapshot); busy["resources"]["free_gpus"] = []
    store.heartbeat(nodes[0]["node_id"], busy, [])
    result = store.submit(request, needs, user, origin_node=nodes[0]["node_id"])
    assert result["job"]["node_id"] in {nodes[1]["node_id"], nodes[2]["node_id"]}


@pytest.mark.parametrize("change", [
    {"images": {"engine:1": "sha256:other"}}, {"images": {}}, {"references": {}},
    {"code_revision": "old-code"}, {"scheduler_online": False}, {"accept_jobs": False},
    {"resources": {"free_gpus": [{"gpu_id": 0, "free_vram_gb": 2}], "cpu_free": 16, "ram_available_gb": 64, "scratch_free_gb": 100}},
])
def test_incompatible_workstations_are_not_selected(cluster, change):
    store, nodes, user, needs, snapshot, request, _ = cluster
    for node in nodes:
        store.heartbeat(node["node_id"], snapshot | change, [])
    result = store.submit(request, needs, user, origin_node=nodes[0]["node_id"])
    assert result["job"]["status"] == "waiting_for_workstation"
    assert result["job"]["node_id"] is None


def test_automatic_dispatch_does_not_overbook_unacknowledged_work(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    jobs = [store.submit(request, needs, user, origin_node=nodes[0]["node_id"])["job"] for _ in range(4)]
    assert len({j["node_id"] for j in jobs[:3]}) == 3
    assert jobs[3]["node_id"] is None


def test_work_request_for_all_gpus_requires_all_local_devices_to_be_free(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    partial = copy.deepcopy(snapshot)
    partial["resources"]["gpus"] = [{"gpu_id": 0}, {"gpu_id": 1}]
    store.heartbeat(nodes[0]["node_id"], partial, [])
    job = store.submit(request, needs | {"all_gpus": True}, user, origin_node=nodes[0]["node_id"])["job"]
    assert job["node_id"] != nodes[0]["node_id"]


def test_analysis_affinity_waits_for_owning_machine_and_never_migrates(cluster):
    store, nodes, user, needs, snapshot, request, clock = cluster
    clock[0] += 180
    store.heartbeat(nodes[1]["node_id"], snapshot, [])
    result = store.submit(request, needs, user, origin_node=nodes[1]["node_id"], affinity_node=nodes[0]["node_id"])
    jid = result["job"]["cluster_job_id"]
    assert result["job"]["status"] == "waiting_for_workstation"
    with pytest.raises(ValueError, match="owns"):
        store.submit(request, needs, user, origin_node=nodes[1]["node_id"], affinity_node=nodes[0]["node_id"], placement=nodes[1]["node_id"])
    store.heartbeat(nodes[0]["node_id"], snapshot, [])
    assert store.job(jid, user)["node_id"] == nodes[0]["node_id"]


def test_finished_job_retains_status_location_and_metrics_when_machine_is_off(cluster):
    store, nodes, user, _, snapshot, _, clock = cluster
    row = {"task_group": "design", "run_id": "run-1", "owner_id": user["id"], "project_id": "project-a", "status": "completed", "metrics": {"score": 0.9}, "run_dir": "/same/path/run-1"}
    store.heartbeat(nodes[0]["node_id"], snapshot, [row])
    jid = config.job_id(nodes[0]["node_id"], "design/run-1")
    clock[0] += 180
    job = store.job(jid, user)
    assert job["status"] == "completed" and job["availability"] == "offline"
    assert job["metrics"]["score"] == .9
    assert job["location"] == "/same/path/run-1"
    assert job["run_dir"].startswith("cluster://" + nodes[0]["node_id"])


def test_same_local_run_name_on_two_nodes_stays_distinct_and_private(cluster):
    store, nodes, user, _, snapshot, _, _ = cluster
    row = {"task_group": "design", "run_id": "same", "owner_id": user["id"], "project_id": "project-a", "status": "completed"}
    for node in nodes[:2]:
        store.heartbeat(node["node_id"], snapshot, [row])
    jobs = store.search(user)["items"]
    assert len(jobs) == 2 and len({j["cluster_job_id"] for j in jobs}) == 2
    outsider = {"id": "outsider", "role": "user", "project_ids": []}
    assert store.search(outsider)["items"] == []
    with pytest.raises(LookupError):
        store.job(jobs[0]["cluster_job_id"], outsider)


def test_submission_idempotency_and_changed_input_conflicts(cluster):
    store, nodes, user, needs, _, request, _ = cluster
    first = store.submit(request, needs, user, origin_node=nodes[0]["node_id"], idempotency_key="one-click")
    again = store.submit(request, needs, user, origin_node=nodes[0]["node_id"], idempotency_key="one-click")
    assert first["job"]["cluster_job_id"] == again["job"]["cluster_job_id"]
    with pytest.raises(ValueError, match="different inputs"):
        store.submit(request | {"bundle_id": "changed"}, needs, user, origin_node=nodes[0]["node_id"], idempotency_key="one-click")


def test_cancel_waiting_submission_remembers_idempotency_and_does_not_dispatch(cluster):
    store, nodes, user, needs, snapshot, request, clock = cluster
    clock[0] += 180
    job = store.submit(request, needs, user, origin_node=nodes[0]["node_id"], idempotency_key="cancel-me")["job"]
    assert job["status"] == "waiting_for_workstation"
    assert store.action(job["cluster_job_id"], "cancel", {}, user)["status"] == "cancelled"
    store.heartbeat(nodes[0]["node_id"], snapshot, [])
    assert store.claim(nodes[0]["node_id"]) is None
    again = store.submit(dict(reversed(list(request.items()))), needs, user, origin_node=nodes[0]["node_id"], idempotency_key="cancel-me")["job"]
    assert again["cluster_job_id"] == job["cluster_job_id"] and again["status"] == "cancelled"


def test_cancel_during_acceptance_is_preserved_until_local_cancel_command(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    job = store.submit(request, needs, user, origin_node=nodes[0]["node_id"], placement="local")["job"]
    command = store.claim(nodes[0]["node_id"])
    assert store.action(job["cluster_job_id"], "cancel", {}, user)["cancel_requested"]
    row = {"cluster_job_id": job["cluster_job_id"], "task_group": "design", "run_id": "accepted", "owner_id": user["id"],
           "project_id": "project-a", "status": "queued"}
    store.heartbeat(nodes[0]["node_id"], snapshot, [row])
    assert store.job(job["cluster_job_id"], user)["cancel_requested"]
    store.finish(nodes[0]["node_id"], command["id"], command["lease_token"], {"job": dict(row, run_dir="/node/results/accepted")})
    cancel = store.claim(nodes[0]["node_id"])
    assert cancel["kind"] == "cancel" and cancel["payload"]["local_reference"] == "design/accepted"


def test_missing_results_are_unavailable_without_erasing_completed_status(cluster):
    store, nodes, user, _, snapshot, _, _ = cluster
    row = {"task_group": "design", "run_id": "saved", "owner_id": user["id"], "project_id": "project-a", "status": "completed"}
    store.heartbeat(nodes[0]["node_id"], snapshot, [row], known_references=["design/saved"])
    jid = config.job_id(nodes[0]["node_id"], "design/saved")
    store.heartbeat(nodes[0]["node_id"], snapshot | {"storage_available": False}, [], known_references=[])
    assert store.job(jid, user)["availability"] == "storage_unavailable"
    store.heartbeat(nodes[0]["node_id"], snapshot, [], known_references=[])
    assert store.job(jid, user)["availability"] == "files_unavailable"
    assert store.job(jid, user)["status"] == "completed"
    store.heartbeat(nodes[0]["node_id"], snapshot, [row], known_references=["design/saved"])
    assert store.job(jid, user)["availability"] == "online"


def test_reference_probes_remain_until_agent_has_accepted_submission(cluster):
    store, nodes, user, needs, _, request, _ = cluster
    store.submit(request, needs, user, origin_node=nodes[0]["node_id"], placement="local")
    command = store.claim(nodes[0]["node_id"])
    assert store.reference_probes() == needs["references"]
    store.finish(nodes[0]["node_id"], command["id"], command["lease_token"], error="Missing input")
    assert store.reference_probes() == []


def test_reference_probe_handles_symlinked_reference_directories(monkeypatch, tmp_path):
    from mn_molecular_design_studio import runtime
    from mn_molecular_design_studio.cluster import capabilities
    refs = tmp_path / "references"; refs.mkdir()
    models = tmp_path / "model-disk"; models.mkdir()
    (models / "model.pt").write_bytes(b"checkpoint")
    (refs / "models").symlink_to(models, target_is_directory=True)
    monkeypatch.setattr(runtime, "reference_root", lambda: refs)
    inventory = capabilities.reference_inventory(["models/model.pt", "../secret", "/etc/passwd"])
    assert inventory == {"models": True, "models/model.pt": 10}


def test_transfer_cleanup_keeps_inputs_awaiting_workstation(cluster):
    import os
    store, nodes, user, needs, _, request, clock = cluster
    clock[0] += 180
    incoming = bundles.blob_path("a" * 64); incoming.write_bytes(b"input")
    expired = bundles.blob_path("b" * 64); expired.write_bytes(b"download")
    for path in (incoming, expired):
        os.utime(path, (clock[0] - 90000, clock[0] - 90000))
    job = store.submit(request | {"bundle_id": incoming.name}, needs, user, origin_node=nodes[0]["node_id"])["job"]
    assert store.cleanup_transfers() == 1 and incoming.exists() and not expired.exists()
    store.action(job["cluster_job_id"], "cancel", {}, user)
    assert store.cleanup_transfers() == 1 and not incoming.exists()


@pytest.mark.parametrize(("engine", "reference"), [
    ("af3_torch", "mn-cofolding/af3_torch"), ("openfold3_ob0", "openfold3/of3-ob-2025-06-30-174k.pt"),
    ("rf3_foundry", "foundry/rf3_foundry_01_24_latest_remapped.ckpt"), ("atlasfold", "atlasfold"),
    ("colabfold", "alphafold_models"), ("opendde_native", "opendde"),
])
def test_cofolding_placement_uses_the_engines_actual_reference_location(engine, reference):
    from mn_molecular_design_studio.cluster import capabilities
    requirement = capabilities.requirements({"workflow": "refolding.optimized_prediction", "parameters": {"engine": engine}}, installed_images={})
    assert reference in requirement["references"]


def test_lost_command_reply_is_reclaimed_only_by_original_node_with_fencing(cluster):
    store, nodes, user, needs, _, request, clock = cluster
    result = store.submit(request, needs, user, origin_node=nodes[0]["node_id"], placement="local")
    command = store.claim(nodes[0]["node_id"])
    clock[0] += 200
    assert store.claim(nodes[1]["node_id"]) is None
    retry = store.claim(nodes[0]["node_id"])
    assert retry["id"] == command["id"] and retry["lease_token"] != command["lease_token"]
    with pytest.raises(PermissionError):
        store.finish(nodes[0]["node_id"], command["id"], command["lease_token"], error="old reply")
    store.finish(nodes[0]["node_id"], retry["id"], retry["lease_token"], {"job": {"task_group": "design", "run_id": "actual", "status": "queued", "run_dir": "/remote/results/actual"}})
    job = store.job(result["job"]["cluster_job_id"], user)
    assert job["node_id"] == nodes[0]["node_id"] and job["local_reference"] == "design/actual"


def test_agent_receipt_prevents_duplicate_inference_after_lost_reply(cluster, monkeypatch):
    store, nodes, user, needs, _, request, _ = cluster
    store.submit(request, needs, user, origin_node=nodes[0]["node_id"], placement="local")
    command = store.claim(nodes[0]["node_id"])
    calls = []
    monkeypatch.setattr(Agent, "submit", lambda self, cmd: calls.append(cmd["id"]) or {"job": {"run_id": "once"}})
    result, error = Agent().execute(command)
    again, again_error = Agent().execute(command)
    assert calls == [command["id"]] and result == again and error == again_error is None


def test_agent_catalog_uses_cluster_submission_id_after_scheduler_rewrites_metadata(monkeypatch, tmp_path):
    from mn_molecular_design_studio import job_catalog

    node_id = str(uuid.uuid4())
    submission_id = str(uuid.uuid4())
    run = tmp_path / "design" / "run-1"
    run.mkdir(parents=True)
    (run / "metadata.json").write_text(json.dumps({"status": "completed", "cluster_submission_id": submission_id}))
    (run / "molecular_design_studio.json").write_text(json.dumps({"cluster_submission_id": submission_id}))
    (run / "result.json").write_text(json.dumps({"success": True, "metrics": {"candidate_count": 1}}))
    monkeypatch.setattr(job_catalog, "search_jobs", lambda **kwargs: {
        "items": [{"task_group": "design", "run_id": "run-1", "run_dir": str(run), "status": "completed",
                   "metadata_mtime_ns": 1, "input_mtime_ns": 1, "result_mtime_ns": 1}],
        "total": 1,
    })
    monkeypatch.setattr(job_catalog, "job_artifacts", lambda *_args: [])

    agent = Agent()
    agent.node_id = node_id
    item = agent.catalog()[0]

    assert item["cluster_job_id"] == submission_id
    assert item["status"] == "completed"


def test_submission_context_writes_both_cluster_job_id_fields(monkeypatch, tmp_path):
    from mn_molecular_design_studio.cluster import context

    node_id = str(uuid.uuid4())
    submission_id = str(uuid.uuid4())
    monkeypatch.setattr(config, "identity", lambda: node_id)
    run = tmp_path / "design" / "run-1"
    run.mkdir(parents=True)
    metadata = {"run_id": "run-1"}
    command = {"job_id": submission_id, "project_id": "project-a", "owner_id": "owner-a"}

    with context.submission(command):
        context.tag_created_job(run, metadata)

    marker = json.loads((run / "molecular_design_studio.json").read_text())
    assert metadata["cluster_job_id"] == metadata["cluster_submission_id"] == submission_id
    assert marker["cluster_job_id"] == marker["cluster_submission_id"] == submission_id


def test_heartbeat_promotes_legacy_fallback_job_to_submission_id_and_replicates_tombstone(cluster, monkeypatch):
    from mn_molecular_design_studio.cluster.replication import Replica

    store, nodes, user, needs, snapshot, request, _clock = cluster
    origin, remote = nodes[:2]
    accepted = store.submit(request, needs, user, origin_node=origin["node_id"], placement=remote["node_id"])["job"]
    run_id = "finished-run"
    local_ref = "design/" + run_id
    fallback_id = config.job_id(remote["node_id"], local_ref)
    fallback = {
        "cluster_job_id": fallback_id,
        "node_id": remote["node_id"],
        "owner_id": user["id"],
        "project_id": "project-a",
        "task_group": "design",
        "run_id": run_id,
        "status": "completed",
        "campaign_name": "smoke",
    }
    with store.connect() as db:
        db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?)", (
            fallback_id, remote["node_id"], local_ref, user["id"], "project-a", json.dumps(fallback), "2026-10-18T00:00:00+00:00",
        ))
    monkeypatch.setattr(config, "identity", lambda: remote["node_id"])
    Replica().capture_jobs()

    completed = fallback | {"cluster_job_id": accepted["cluster_job_id"], "metrics": {"candidate_count": 1}}
    store.heartbeat(remote["node_id"], snapshot, [completed])
    Replica().capture_jobs()

    result = store.job(accepted["cluster_job_id"], user)
    assert result["status"] == "completed"
    assert result["local_reference"] == local_ref
    with store.connect() as db:
        assert db.execute("SELECT 1 FROM jobs WHERE id=?", (fallback_id,)).fetchone() is None
    with Replica().db() as db:
        tombstone = db.execute("SELECT value FROM records WHERE key=?", ("job:" + fallback_id,)).fetchone()
    assert tombstone is not None and json.loads(tombstone["value"]) is None


def test_agent_crash_window_is_uncertain_and_never_submits_again(cluster, monkeypatch):
    store, nodes, user, needs, _, request, _ = cluster
    store.submit(request, needs, user, origin_node=nodes[0]["node_id"], placement="local")
    command = store.claim(nodes[0]["node_id"])
    directory = config.home() / "receipts"; directory.mkdir()
    config.atomic_json(directory / (command["id"] + ".json"), {"state": "accepting"})
    monkeypatch.setattr(Agent, "submit", lambda *a: pytest.fail("Cannot repeat a possibly accepted scientific job"))
    _, error = Agent().execute(command)
    assert "uncertain" in error


def test_machine_qualified_path_never_resolves_on_other_workstation(cluster, monkeypatch, tmp_path):
    _, nodes, *_ = cluster
    root = tmp_path / "runs"; run = root / "design/same"; run.mkdir(parents=True)
    uri = config.artifact_uri(nodes[0]["node_id"], "design", "same")
    assert bundles.resolve_artifact(uri, node_id=nodes[0]["node_id"], root=root) == run
    with pytest.raises(ValueError, match="another workstation"):
        bundles.resolve_artifact(uri, node_id=nodes[1]["node_id"], root=root)
    with pytest.raises(ValueError):
        bundles.resolve_artifact(uri + "/../../../outside", node_id=nodes[0]["node_id"], root=root)


def test_input_staging_copies_only_inputs_and_keeps_runs_and_weights_in_place(cluster, monkeypatch, tmp_path):
    _, nodes, *_ = cluster
    from mn_molecular_design_studio import runtime
    run_root = tmp_path / "runs"; run = run_root / "design/run"; run.mkdir(parents=True)
    refs = tmp_path / "refs"; refs.mkdir(); weight = refs / "model.pt"; weight.write_bytes(b"weights")
    structure = tmp_path / "target.pdb"; structure.write_text("TARGET")
    monkeypatch.setattr(runtime, "runs_root", lambda: run_root)
    monkeypatch.setattr(runtime, "reference_root", lambda: refs)
    from mn_molecular_design_studio.engine.workflow_engine.core import portable_paths
    monkeypatch.setattr(portable_paths, "runs_root", lambda: run_root)
    monkeypatch.setattr(portable_paths, "reference_root", lambda: refs)
    prepared, digest, references = bundles.pack({"parameters": {"target_pdb": str(structure), "source_run_dir": str(run), "checkpoint_path": str(weight)}},
                                               base_dir=tmp_path, node_id=nodes[0]["node_id"])
    assert prepared["parameters"]["target_pdb"].startswith("cluster-input:///")
    assert prepared["parameters"]["source_run_dir"].startswith("cluster://")
    assert prepared["parameters"]["checkpoint_path"] == "reference:///model.pt"
    assert references == ["model.pt"]
    with zipfile.ZipFile(bundles.blob_path(digest)) as archive:
        assert len(archive.namelist()) == 1
    destination = tmp_path / "staged"
    bundles.unpack(bundles.blob_path(digest), destination, expected_digest=digest)
    restored = bundles.restore(prepared, staging=destination, node_id=nodes[0]["node_id"])
    assert Path(restored["parameters"]["target_pdb"]).read_text() == "TARGET"
    assert restored["parameters"]["source_run_dir"] == str(run)


@pytest.mark.parametrize("member", ["../outside", "/absolute", "inputs/../../../outside"])
def test_input_archive_rejects_path_traversal(cluster, tmp_path, member):
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(member, "bad")
    with pytest.raises(ValueError, match="Unsafe"):
        bundles.unpack(archive, tmp_path / "stage")
    assert not (tmp_path / "outside").exists()


def test_input_bundle_digest_mismatch_rejected(cluster, tmp_path):
    archive = tmp_path / "bad.zip"; archive.write_bytes(b"damaged")
    with pytest.raises(ValueError, match="SHA-256"):
        bundles.unpack(archive, tmp_path / "stage", expected_digest="0" * 64)


def test_node_token_cannot_publish_another_nodes_jobs(cluster):
    store, nodes, user, _, snapshot, _, _ = cluster
    jid = config.job_id(nodes[0]["node_id"], "design/run")
    row = {"task_group": "design", "run_id": "run", "cluster_job_id": jid, "owner_id": user["id"], "project_id": "project-a", "status": "completed"}
    store.heartbeat(nodes[0]["node_id"], snapshot, [row])
    with pytest.raises(PermissionError):
        store.heartbeat(nodes[1]["node_id"], snapshot, [row])
    assert store.authenticate_node(nodes[0]["node_token"])["id"] == nodes[0]["node_id"]
    with pytest.raises(PermissionError):
        store.authenticate_node("incorrect")


def test_remote_job_marker_does_not_write_to_local_files(cluster, tmp_path):
    from mn_molecular_design_studio.workspace import mark_studio_job
    result = {"job": {"cluster_job_id": str(uuid.uuid4()), "run_dir": str(tmp_path / "must-not-create")}}
    assert mark_studio_job(result) == result
    assert not (tmp_path / "must-not-create").exists()


def test_api_enrollment_catalog_and_agent_commands_are_authenticated(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    from mn_molecular_design_studio import workspace
    workspace.ensure_studio_workspace()
    from mn_molecular_design_studio import api
    config.save({"mode": "coordinator", "coordinator_url": "http://127.0.0.1:8766"})
    user = {"id": "owner", "role": "admin", "must_change_password": False}
    api.app.dependency_overrides[api._current_user] = lambda: user
    try:
        with TestClient(api.app) as client:
            enrollment = client.post("/api/v1/cluster/nodes", json={"name": "4032"})
            assert enrollment.status_code == 201, enrollment.text
            token = enrollment.json()["node_token"]
            assert client.post("/api/v1/cluster/node/heartbeat", json={}).status_code == 401
            assert client.post("/api/v1/cluster/node/heartbeat", json={}, headers={"Authorization": "Bearer wrong"}).status_code == 403
            response = client.post("/api/v1/cluster/node/heartbeat", json={"snapshot": {"accept_jobs": True}, "jobs": []}, headers={"Authorization": "Bearer " + token})
            assert response.status_code == 200
            nodes = client.get("/api/v1/cluster/nodes").json()["nodes"]
            assert nodes[0]["availability"] == "online"
            assert "token_hash" not in nodes[0] and "node_token" not in nodes[0]
            assert client.get("/api/v1/cluster/jobs").json()["items"] == []
    finally:
        api.app.dependency_overrides.pop(api._current_user, None)
