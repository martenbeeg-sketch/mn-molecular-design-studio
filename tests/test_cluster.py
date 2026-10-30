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
from mn_molecular_design_studio.cluster.agent import Agent, _candidate_summary
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


def test_liveness_pulse_refreshes_node_without_running_catalog_dispatch(cluster, monkeypatch):
    store, nodes, _, _, _, _, clock = cluster
    dispatches = []
    monkeypatch.setattr(store, "dispatch", lambda: dispatches.append(True))
    clock[0] += 180

    store.pulse(nodes[0]["node_id"], {"scheduler_online": True})

    current = next(node for node in store.nodes() if node["id"] == nodes[0]["node_id"])
    assert current["availability"] == "online"
    assert current["snapshot"] == {"scheduler_online": True}
    assert dispatches == []


def test_agent_heartbeat_interval_stays_below_unreachable_threshold(monkeypatch):
    agent = Agent()
    monkeypatch.setattr(config, "load", lambda: {"heartbeat_seconds": 300})

    assert agent._heartbeat_interval() == 15


def test_worker_heartbeat_publishes_counts_without_candidate_identities_or_structures(tmp_path):
    run = tmp_path / "design" / "run-1"
    candidates = run / "artifacts/normalized_candidates"
    candidates.mkdir(parents=True)
    (run / "design.pdb").write_text("private coordinates")
    rows = [
        {"candidate_id": "private-design-1", "source_tool": "rfd4", "stage": "complex_refolding",
         "binder_sequence": "ACDE", "metrics": {"native_pass_filters": True}, "raw_metadata": {}},
        {"candidate_id": "private-design-2", "source_tool": "rfd4", "stage": "complex_refolding",
         "metrics": {"native_pass_filters": False}, "raw_metadata": {}},
        {"candidate_id": "private-design-3", "source_tool": "rfd4", "stage": "complex_refolding",
         "metrics": {}, "raw_metadata": {"validation_acceptance": {"status": "accepted"}}},
    ]
    (candidates / "candidates.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))

    summary = _candidate_summary(run, {"tool": "rfd4"})

    assert summary["candidate_count"] == 3
    assert summary["eligible_count"] == 3
    assert summary["native_filter_counts"] == {"passed": 1, "failed": 1, "unavailable": 1}
    assert summary["acceptance_counts"]["accepted"] == 1
    assert not ({"candidates", "artifact_sha256", "artifact_path"} & set(summary))
    assert "private-design-1" not in json.dumps(summary)
    assert "ACDE" not in json.dumps(summary)
    assert "private coordinates" not in json.dumps(summary)


def test_automatic_placement_uses_other_ready_workstation_when_origin_is_busy(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    busy = copy.deepcopy(snapshot); busy["resources"]["free_gpus"] = []
    store.heartbeat(nodes[0]["node_id"], busy, [])
    result = store.submit(request, needs, user, origin_node=nodes[0]["node_id"])
    assert result["job"]["node_id"] in {nodes[1]["node_id"], nodes[2]["node_id"]}


def test_synchronized_running_job_reserves_its_gpu_for_automatic_placement(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    busy = copy.deepcopy(snapshot)
    busy["resources"]["cpu_total"] = 16
    running = {
        "task_group": "design", "run_id": "active-gpu-job", "owner_id": user["id"],
        "project_id": "project-a", "status": "running", "queue_resource": "gpu:0",
        "resource_allocation": "4 CPU cores · 0",
    }
    store.heartbeat(nodes[0]["node_id"], busy, [running])

    candidates = store.placement_options(needs, origin_node=nodes[0]["node_id"])
    candidate_ids = {candidate["node_id"] for candidate in candidates}
    assert nodes[0]["node_id"] not in candidate_ids
    assert candidate_ids == {nodes[1]["node_id"], nodes[2]["node_id"]}

    result = store.submit(request, needs, user, origin_node=nodes[0]["node_id"])
    assert result["job"]["node_id"] in candidate_ids


def test_explicit_workstation_selection_is_rechecked_before_dispatch(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    busy = copy.deepcopy(snapshot)
    running = {
        "task_group": "design", "run_id": "active-gpu-job", "owner_id": user["id"],
        "project_id": "project-a", "status": "running", "queue_resource": "gpu:0",
        "resource_allocation": "4 CPU cores · 0",
    }
    store.heartbeat(nodes[1]["node_id"], busy, [running])

    result = store.submit(request, needs, user, origin_node=nodes[0]["node_id"], placement=nodes[1]["node_id"])
    assert result["job"]["status"] == "waiting_for_workstation"
    assert result["job"]["node_id"] == nodes[1]["node_id"]


def test_explicit_origin_selection_survives_stale_package_revision(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    current = copy.deepcopy(snapshot)
    current["code_revision"] = "updated-code"
    store.heartbeat(nodes[0]["node_id"], current, [])

    result = store.submit(request, needs, user, origin_node=nodes[0]["node_id"],
                         placement=nodes[0]["node_id"])

    assert result["job"]["status"] == "dispatching"
    with store.connect() as db:
        row = db.execute("SELECT s.requirements,c.payload FROM submissions s JOIN commands c ON c.id=s.command_id WHERE s.job_id=?",
                         (result["job"]["cluster_job_id"],)).fetchone()
    original = json.loads(row["requirements"])
    command_payload = json.loads(row["payload"])
    assert original["code_revision"] == needs["code_revision"]
    assert command_payload["allow_code_revision_refresh"] is True
    assert command_payload["origin_node"] == nodes[0]["node_id"]
    claimed = store.claim(nodes[0]["node_id"])
    assert claimed["node_id"] == nodes[0]["node_id"]


def test_retry_recovers_job_rejected_only_for_stale_compatibility_fingerprint(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    submitted = store.submit(request, needs, user, origin_node=nodes[0]["node_id"],
                             placement=nodes[0]["node_id"])["job"]
    error = "The workstation's engine images, code, or checkpoints no longer match this submission."
    with store.connect() as db:
        submission = db.execute("SELECT command_id FROM submissions WHERE job_id=?",
                                (submitted["cluster_job_id"],)).fetchone()
        db.execute("UPDATE commands SET status='failed',error=? WHERE id=?", (error, submission["command_id"]))
        job = json.loads(db.execute("SELECT payload FROM jobs WHERE id=?",
                                    (submitted["cluster_job_id"],)).fetchone()["payload"])
        job.update(status="failed", error=error)
        db.execute("UPDATE jobs SET payload=? WHERE id=?",
                   (json.dumps(job), submitted["cluster_job_id"]))

    retried = store.action(submitted["cluster_job_id"], "retry", {}, user)

    assert retried["status"] == "dispatching"
    assert not retried.get("error")


def test_waiting_local_job_can_be_moved_to_an_available_gpu(cluster):
    store, nodes, user, needs, snapshot, request, _ = cluster
    selected_needs = dict(needs, gpu_count=1, gpu_ids=[1], min_vram_gb=8)
    selected_request = copy.deepcopy(request)
    selected_request["payload"]["parameters"] = {"gpu_device": "1"}
    submitted = store.submit(selected_request, selected_needs, user,
                             origin_node=nodes[0]["node_id"], placement=nodes[0]["node_id"])["job"]
    assert submitted["status"] == "waiting_for_workstation"

    changed = store.action(submitted["cluster_job_id"], "change_gpu", {"gpu_id": 0}, user)

    assert changed["status"] == "dispatching"
    assert changed["request_summary"]["parameters"]["gpu_device"] == "0"
    with store.connect() as db:
        row = db.execute("SELECT s.requirements,c.payload FROM submissions s JOIN commands c ON c.id=s.command_id WHERE s.job_id=?",
                         (submitted["cluster_job_id"],)).fetchone()
    req = json.loads(row["requirements"])
    command_payload = json.loads(row["payload"])
    assert req["gpu_ids"] == [0]
    assert command_payload["payload"]["parameters"]["gpu_device"] == "0"


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


def test_reference_probes_drop_after_a_submission_is_rejected(cluster):
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


def test_reference_probes_remain_available_after_a_workflow_is_accepted(cluster):
    store, nodes, user, needs, _snapshot, request, _clock = cluster
    submitted = store.submit(
        request, needs, user, origin_node=nodes[0]["node_id"], placement="local",
    )["job"]
    command = store.claim(nodes[0]["node_id"])
    assert store.reference_probes() == ["weights/model.pt"]

    store.finish(
        nodes[0]["node_id"], command["id"], command["lease_token"],
        {"job": {"task_group": "design", "run_id": "accepted", "status": "queued",
                 "run_dir": "/remote/results/design/accepted"}},
    )

    assert store.job(submitted["cluster_job_id"], user)["status"] == "queued"
    assert store.reference_probes() == ["weights/model.pt"]

    store.heartbeat(nodes[0]["node_id"], _snapshot, [{
        "cluster_job_id": submitted["cluster_job_id"], "task_group": "design", "run_id": "accepted",
        "owner_id": user["id"], "project_id": "project-a", "status": "completed",
    }])

    assert store.reference_probes() == []


def test_reference_fingerprint_cache_skips_unchanged_weights_and_refreshes_changes(
    monkeypatch, tmp_path,
):
    from mn_molecular_design_studio import cache_provenance, runtime
    from mn_molecular_design_studio.cluster import capabilities

    root = tmp_path / "reference_files"
    weight = root / "weights/model.pt"
    weight.parent.mkdir(parents=True)
    weight.write_bytes(b"checkpoint-v1")
    monkeypatch.setattr(runtime, "reference_root", lambda: root)
    original_digest = cache_provenance.file_digest
    digests = []

    def counted_digest(path):
        digests.append(Path(path))
        return original_digest(path)

    monkeypatch.setattr(cache_provenance, "file_digest", counted_digest)
    cache = {}
    first = capabilities.reference_fingerprints(["weights/model.pt"], cache)
    second = capabilities.reference_fingerprints(["weights/model.pt"], cache)
    assert first == second and len(digests) == 1

    weight.write_bytes(b"checkpoint-v2-with-new-size")
    refreshed = capabilities.reference_fingerprints(["weights/model.pt"], cache)
    assert refreshed != first and len(digests) == 2


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
    from mn_molecular_design_studio.auth import development_user

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
    assert metadata["studio_development_mode"] is False
    assert marker["studio_development_mode"] is False

    dev_run = tmp_path / "design" / "dev-run"
    dev_run.mkdir(parents=True)
    dev_metadata = {"run_id": "dev-run"}
    dev_command = {"job_id": str(uuid.uuid4()), "project_id": "project-a", "owner_id": development_user()["id"]}
    with context.submission(dev_command):
        context.tag_created_job(dev_run, dev_metadata)
    dev_marker = json.loads((dev_run / "molecular_design_studio.json").read_text())
    assert dev_metadata["studio_development_mode"] is True
    assert dev_marker["studio_development_mode"] is True


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
    assert result["wait_reason"] == ""
    with store.connect() as db:
        assert db.execute("SELECT 1 FROM jobs WHERE id=?", (fallback_id,)).fetchone() is None
    with Replica().db() as db:
        tombstone = db.execute("SELECT value FROM records WHERE key=?", ("job:" + fallback_id,)).fetchone()
    assert tombstone is not None and json.loads(tombstone["value"]) is None


def test_heartbeat_backfills_legacy_empty_ownership_for_same_local_run(cluster):
    store, nodes, user, _needs, snapshot, _request, _clock = cluster
    node = nodes[0]
    local_ref = "analysis/legacy-pipeline-child"
    jid = config.job_id(node["node_id"], local_ref)
    legacy = {
        "cluster_job_id": jid,
        "node_id": node["node_id"],
        "local_reference": local_ref,
        "task_group": "analysis",
        "run_id": "legacy-pipeline-child",
        "owner_id": "",
        "project_id": "",
        "status": "completed",
    }
    with store.connect() as db:
        db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?)", (
            jid, node["node_id"], local_ref, "", "", json.dumps(legacy), "2026-10-23T00:00:00+00:00",
        ))

    updated = legacy | {"owner_id": user["id"], "project_id": "project-a"}
    store.heartbeat(node["node_id"], snapshot, [updated])

    with store.connect() as db:
        row = db.execute("SELECT owner_id,project_id FROM jobs WHERE id=?", (jid,)).fetchone()
    assert row["owner_id"] == user["id"]
    assert row["project_id"] == "project-a"


def test_heartbeat_still_rejects_conflicting_ownership_for_same_local_run(cluster):
    store, nodes, user, _needs, snapshot, _request, _clock = cluster
    node = nodes[0]
    local_ref = "analysis/owned-pipeline-child"
    jid = config.job_id(node["node_id"], local_ref)
    existing = {
        "cluster_job_id": jid,
        "node_id": node["node_id"],
        "local_reference": local_ref,
        "task_group": "analysis",
        "run_id": "owned-pipeline-child",
        "owner_id": user["id"],
        "project_id": "project-a",
        "status": "completed",
    }
    with store.connect() as db:
        db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?)", (
            jid, node["node_id"], local_ref, user["id"], "project-a", json.dumps(existing), "2026-10-23T00:00:00+00:00",
        ))

    conflicting = existing | {"owner_id": "other-user", "project_id": "other-project"}
    with pytest.raises(PermissionError, match="cannot change job ownership or project"):
        store.heartbeat(node["node_id"], snapshot, [conflicting])


def test_replica_import_promotes_canonical_job_before_fallback_tombstone(cluster, monkeypatch):
    from mn_molecular_design_studio.cluster.replication import Replica

    store, nodes, user, needs, _snapshot, request, _clock = cluster
    origin, remote = nodes[:2]
    accepted = store.submit(request, needs, user, origin_node=origin["node_id"], placement=remote["node_id"])["job"]
    local_ref = "design/finished-run"
    fallback_id = config.job_id(remote["node_id"], local_ref)
    fallback = {"cluster_job_id": fallback_id, "node_id": remote["node_id"], "owner_id": user["id"],
                "project_id": "project-a", "status": "completed", "local_reference": local_ref}
    with store.connect() as db:
        db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?)", (
            fallback_id, remote["node_id"], local_ref, user["id"], "project-a", json.dumps(fallback), "2026-10-18T00:00:00+00:00",
        ))
    monkeypatch.setattr(config, "identity", lambda: remote["node_id"])
    canonical = fallback | {"cluster_job_id": accepted["cluster_job_id"]}

    Replica()._apply("job", "job:" + accepted["cluster_job_id"], canonical, remote["node_id"])

    with store.connect() as db:
        assert db.execute("SELECT 1 FROM jobs WHERE id=?", (fallback_id,)).fetchone() is None
        row = db.execute("SELECT local_ref,payload FROM jobs WHERE id=?", (accepted["cluster_job_id"],)).fetchone()
    assert row["local_ref"] == local_ref
    assert json.loads(row["payload"])["cluster_job_id"] == accepted["cluster_job_id"]


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
    # A local directory named "pdb" must not make a prediction output format
    # look like a remote filesystem input during request packaging.
    (tmp_path / "pdb").mkdir()
    refs = tmp_path / "refs"; refs.mkdir(); weight = refs / "model.pt"; weight.write_bytes(b"weights")
    structure = tmp_path / "target.pdb"; structure.write_text("TARGET")
    monkeypatch.setattr(runtime, "runs_root", lambda: run_root)
    monkeypatch.setattr(runtime, "reference_root", lambda: refs)
    from mn_molecular_design_studio.engine.workflow_engine.core import portable_paths
    monkeypatch.setattr(portable_paths, "runs_root", lambda: run_root)
    monkeypatch.setattr(portable_paths, "reference_root", lambda: refs)
    prepared, digest, references = bundles.pack({"parameters": {"target_pdb": str(structure), "source_run_dir": str(run), "checkpoint_path": str(weight), "engine_options": {"output_format": "pdb"}}},
                                               base_dir=tmp_path, node_id=nodes[0]["node_id"])
    assert prepared["parameters"]["target_pdb"].startswith("cluster-input:///")
    assert prepared["parameters"]["source_run_dir"].startswith("cluster://")
    assert prepared["parameters"]["checkpoint_path"] == "reference:///model.pt"
    assert prepared["parameters"]["engine_options"]["output_format"] == "pdb"
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


def test_results_can_hide_and_restore_jobs_without_losing_scheduler_record(cluster):
    store, nodes, user, _, snapshot, _, _ = cluster
    jid = config.job_id(nodes[0]["node_id"], "design/run-visibility")
    row = {
        "task_group": "design", "run_id": "run-visibility", "cluster_job_id": jid,
        "owner_id": user["id"], "project_id": "project-a", "status": "completed",
    }
    store.heartbeat(nodes[0]["node_id"], snapshot, [row])
    assert store.search(user)["total"] == 1

    assert store.action(jid, "hide", {}, user)["hidden"] is True
    assert store.search(user)["total"] == 0
    hidden = store.search(user, include_hidden=True)["items"][0]
    assert hidden["hidden"] is True
    assert hidden["run_id"] == "run-visibility"

    # A heartbeat from the owner does not undo the user's Results cleanup.
    store.heartbeat(nodes[0]["node_id"], snapshot, [row])
    assert store.search(user)["total"] == 0
    assert store.action(jid, "restore", {}, user)["hidden"] is False
    assert store.search(user)["total"] == 1


def test_remote_purge_command_removes_coordinator_result_after_workstation_finishes(cluster):
    store, nodes, user, *_ = cluster
    jid = config.job_id(nodes[1]["node_id"], "design/failed-run")
    store.heartbeat(nodes[1]["node_id"], {}, [{
        "cluster_job_id": jid, "task_group": "design", "run_id": "failed-run",
        "owner_id": user["id"], "project_id": "project-a", "status": "failed",
    }])

    operation = store.action(jid, "purge", {"actor": user}, user)
    command = store.claim(nodes[1]["node_id"])
    assert command["id"] == operation["command_id"]
    assert command["kind"] == "purge"

    store.finish(nodes[1]["node_id"], command["id"], command["lease_token"], response={"status": "purged"})

    with pytest.raises(LookupError):
        store.job(jid, user)
    assert store.command(command["id"], user)["status"] == "completed"


def test_agent_purge_deletes_local_run_folder_and_accepted_inputs(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    from mn_molecular_design_studio import job_catalog, molecular_database_results
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    runs = tmp_path / "runs"
    run_dir = runs / "vina" / "failed-run"
    run_dir.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps({"status": "completed"}))
    accepted_id = str(uuid.uuid4())
    accepted_dir = config.home() / "accepted-inputs" / accepted_id
    accepted_dir.mkdir(parents=True)
    (accepted_dir / "input.json").write_text("{}")
    receipt_dir = config.home() / "receipts"
    receipt_dir.mkdir(exist_ok=True)
    (receipt_dir / f"{accepted_id}.json").write_text("{}")
    monkeypatch.setattr(jobs, "runs_root", lambda: runs)
    monkeypatch.setattr(job_catalog, "reconcile", lambda _root: None)
    monkeypatch.setattr(molecular_database_results, "remove_run_records", lambda *args, **kwargs: {
        "removed_sources": 2, "removed_members": 3, "removed_groups": 0,
        "preserved_saved_selections": True,
    })

    result = Agent().purge({
        "payload": {
            "local_reference": "vina/failed-run", "accepted_command_id": accepted_id,
            "actor": {"id": "owner", "role": "admin"},
            "run_owner_id": "owner", "run_project_id": "project-a",
        },
    })

    assert result["run_files_removed"] is True
    assert result["staging_inputs_removed"] is True
    assert result["database"]["preserved_saved_selections"] is True
    assert not run_dir.exists()
    assert not accepted_dir.exists()
    assert not (receipt_dir / f"{accepted_id}.json").exists()


def test_synchronized_results_page_has_status_filter_and_row_selection(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio.cluster import ui
    from mn_molecular_design_studio import workspace

    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    requests = []
    row = {
        "cluster_job_id": "job-1", "node_id": "node-1", "node_name": "pc002975",
        "task_group": "vina", "run_id": "run-1", "owner_id": "owner",
        "project_id": "", "status": "completed", "availability": "online",
    }
    monkeypatch.setattr(ui, "catalog", lambda user, token, **kwargs: ({"items": [row], "total": 1}, False))
    def request(path, token, payload=None, method="GET", **kwargs):
        requests.append((path, method, payload))
        return {"status": "purged", "database": {"removed_sources": 0, "removed_members": 0}}
    monkeypatch.setattr(ui, "_request", request)
    monkeypatch.setattr(ui.config, "identity", lambda: "node-2")
    monkeypatch.setattr(ui.config, "display_name", lambda value: str(value or "All workstations"))
    monkeypatch.setattr(ui.config, "display_text", lambda value: str(value or ""))
    monkeypatch.setattr(ui, "_row_progress", lambda row: None)
    monkeypatch.setattr(workspace, "load_organization", lambda: {"projects": []})
    user = {"id": "owner", "role": "user", "project_ids": []}
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.cluster.ui import _render_runs_content, _update_results_table_selection
if st.session_state.pop("test-table-row-edit", False):
    _update_results_table_selection("test-table-editor", [[{row!r}]])
_render_runs_content({user!r}, "token", lambda **kwargs: {{"items": [], "total": 0}})
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert app.multiselect(key="cluster-results-status-filter").options == ["Completed"]
    assert "Select" in app.dataframe[0].value.columns
    # AppTest cannot directly edit data_editor cells. Exercise the widget's
    # real on_change callback with the same edited_rows event it receives.
    app.session_state["test-table-editor"] = {"edited_rows": {0: {"Select": True}}}
    app.session_state["test-table-row-edit"] = True
    app.run()
    assert not app.exception
    assert app.button(key="cluster-results-remove-and-purge-request")
    assert not any(button.key in {"cluster-results-hide-selected", "cluster-results-database-removal-request"}
                   for button in app.button)

    app.button(key="cluster-results-remove-and-purge-request").click().run()
    assert not app.exception
    assert app.button(key="cluster-results-remove-and-purge-confirm")
    app.button(key="cluster-results-remove-and-purge-confirm").click().run()
    assert not app.exception
    assert ("/jobs/job-1/purge", "DELETE", None) in requests
    assert app.session_state["cluster-results-selected"] == {}
    assert not app.dataframe[0].value["Select"].any()


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
            pulse = client.post("/api/v1/cluster/node/pulse", json={"snapshot": {"accept_jobs": True}}, headers={"Authorization": "Bearer " + token})
            assert pulse.status_code == 200
            assert pulse.json()["node_id"] == enrollment.json()["node_id"]
            nodes = client.get("/api/v1/cluster/nodes").json()["nodes"]
            assert nodes[0]["availability"] == "online"
            assert "token_hash" not in nodes[0] and "node_token" not in nodes[0]
            assert client.get("/api/v1/cluster/jobs").json()["items"] == []
    finally:
        api.app.dependency_overrides.pop(api._current_user, None)
