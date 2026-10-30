from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows import distributed_campaign
from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import CampaignLedger


def ledger(tmp_path, *, max_attempts=12, target=4):
    return CampaignLedger.create(
        tmp_path / "campaign.json",
        {"accepted_target": target, "max_attempts": max_attempts, "batch_size": 2, "seed": 101},
        campaign_id="campaign-test",
    )


def test_campaign_authorization_is_atomic_and_seeds_are_unique(tmp_path):
    campaign = ledger(tmp_path, max_attempts=7)
    with ThreadPoolExecutor(max_workers=8) as pool:
        reservations = list(pool.map(lambda index: campaign.authorize(3, batch_id=f"b{index}"), range(8)))
    attempts = [attempt for group in reservations for attempt in group]
    assert len(attempts) == 7
    assert len({item["attempt_id"] for item in attempts}) == 7
    assert len({item["seed"] for item in attempts}) == 7
    assert campaign.counters()["attempts_authorized"] == 7
    assert campaign.authorize(1, batch_id="last") == []


def test_result_receipts_are_idempotent_and_acceptance_levels_are_independent(tmp_path):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="b1")[0]
    campaign.transition(attempt["attempt_id"], "running")
    row = {"candidate_id": "campaign-test-a00000001-c1", "metrics": {"native_pass_filters": True}}
    native_state = {"trajectories": 1, "accepted": 1}
    assert campaign.record_receipt(attempt["attempt_id"], "receipt-1", candidates=[row], native_state=native_state)
    assert not campaign.record_receipt(attempt["attempt_id"], "receipt-1", candidates=[row], native_state=native_state)
    assert not campaign.record_receipt(attempt["attempt_id"], "receipt-2", candidates=[row], native_state=native_state)
    counters = campaign.counters()
    assert counters["attempts_completed"] == 1
    assert counters["native_rejected_trajectories"] == 0
    assert counters["native_accepted_candidates"] == 1
    assert counters["validated_accepted_candidates"] == 0
    campaign.set_candidate_acceptance(row["candidate_id"], validated="accepted", diversity="pending")
    counters = campaign.counters()
    assert counters["native_accepted_candidates"] == 1
    assert counters["validated_accepted_candidates"] == 1
    assert counters["diversity_qualified_accepted_candidates"] == 0


def test_deferred_worker_receipt_counts_designs_without_recording_candidate_payloads(tmp_path):
    campaign = ledger(tmp_path, max_attempts=6, target=5)
    attempt = campaign.authorize(1, batch_id="remote-batch")[0]
    campaign.transition(attempt["attempt_id"], "running")
    assert campaign.record_deferred_receipt(
        attempt["attempt_id"], "deferred-receipt",
        native_state={"output_storage": "worker_local", "output_candidates": 3},
        native_output_ids=[], candidate_records=[], candidate_count=3,
        evidence={"kind": "worker_published_candidate_count", "source_job": "cluster://worker/design/run"},
    )
    counters = campaign.counters()
    receipt = campaign.snapshot()["deferred_outputs"]["deferred-receipt"]
    assert counters["generated_candidates"] == 3
    assert counters["deferred_candidates"] == 3
    assert counters["attempts_completed"] == 1
    assert counters["native_accepted_candidates"] == 0
    assert receipt["candidate_count"] == 3 and receipt["candidates"] == []


def test_retry_preserves_attempt_identity_and_seed(tmp_path):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="b1")[0]
    campaign.transition(attempt["attempt_id"], "failed", failure="worker")
    retried = campaign.retry(attempt["attempt_id"])
    assert retried["attempt_id"] == attempt["attempt_id"]
    assert retried["seed"] == attempt["seed"]
    assert retried["retry_count"] == 1
    assert campaign.counters()["attempts_authorized"] == 1


def test_rfd4_retry_is_only_safe_before_native_container_launch(tmp_path):
    child = tmp_path / "rfd4-child"
    child.mkdir()
    assert distributed_campaign._rfd4_retry_safe(child)
    (child / "command.json").write_text(json.dumps({"mode": "local_worker"}))
    assert distributed_campaign._rfd4_retry_safe(child)
    (child / "command.json").write_text(json.dumps({"mode": "docker"}))
    assert not distributed_campaign._rfd4_retry_safe(child)
    (child / "command.json").write_text("not-json")
    assert not distributed_campaign._rfd4_retry_safe(child)


def test_gpu_compute_receipts_are_durable_and_idempotent(tmp_path):
    campaign = ledger(tmp_path, max_attempts=1, target=1)
    attempt = campaign.authorize(1, batch_id="b1")[0]
    assert campaign.record_compute(attempt["attempt_id"], "worker-run", 12.5)
    assert not campaign.record_compute(attempt["attempt_id"], "worker-run", 12.5)
    assert campaign.gpu_compute_seconds() == 12.5


def test_budget_termination_is_explicit_and_stable(tmp_path):
    campaign = ledger(tmp_path, max_attempts=1, target=1)
    campaign.authorize(1, batch_id="b1")
    state = campaign.terminate("budget_exhausted", detail="one trajectory rejected")
    assert state["status"] == "completed"
    assert state["termination_reason"] == "budget_exhausted"
    assert campaign.authorize(1, batch_id="b2") == []
    assert campaign.terminate("target_reached")["termination_reason"] == "budget_exhausted"


def test_campaign_finalization_promotes_budget_exhaustion_when_target_was_reached(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    parent = tmp_path / "parent"
    (parent / "artifacts/distributed_campaign").mkdir(parents=True)
    jobs.write_json(parent / "metadata.json", {})
    campaign = CampaignLedger.create(
        parent / "artifacts/distributed_campaign/ledger.json",
        {"accepted_target": 1, "max_attempts": 1, "acceptance_level": "native"},
        campaign_id="budget-target-campaign",
    )
    attempt = campaign.authorize(1, batch_id="batch")[0]
    candidate = {"candidate_id": "accepted-candidate", "metrics": {"native_pass_filters": True}}
    campaign.record_receipt(
        attempt["attempt_id"], "receipt", candidates=[candidate],
        native_state={"trajectories": 1, "accepted": 1},
    )
    campaign.terminate("budget_exhausted")
    captured = {}
    monkeypatch.setattr(jobs, "finish_job", lambda _path, success, result: captured.update(success=success, result=result))

    distributed_campaign._finish_campaign(parent, campaign, [candidate])

    metrics = captured["result"]["metrics"]
    assert metrics["termination_reason"] == "target_reached"
    assert metrics["campaign_success"] is True
    assert metrics["target_met"] is True
    assert metrics["generation_budget_exhausted"] is True
    assert metrics["attempts_remaining"] == 0
    assert campaign.snapshot()["termination_reason"] == "target_reached"


@pytest.mark.parametrize(
    "config",
    [
        {"accepted_target": 0, "max_attempts": 2},
        {"accepted_target": 3, "max_attempts": 2},
        {"accepted_target": 1, "max_attempts": 1, "acceptance_level": "shortlist"},
        {"accepted_target": 1, "max_attempts": 1, "diversity": {"sequence_identity_max": 0}},
        {"accepted_target": 1, "max_attempts": 1, "max_compute_seconds": float("nan")},
    ],
)
def test_invalid_campaign_policy_is_rejected(tmp_path, config):
    with pytest.raises(ValueError):
        CampaignLedger.create(tmp_path / "campaign.json", config)


def test_validated_acceptance_requires_explicit_stage_measurements(tmp_path):
    with pytest.raises(ValueError, match="explicit measurements"):
        CampaignLedger.create(
            tmp_path / "campaign.json",
            {
                "accepted_target": 1, "max_attempts": 1, "acceptance_level": "validated",
                "validation_stages": [{"workflow": "refolding.optimized_prediction"}],
            },
        )


@pytest.mark.parametrize("workstation_selection, busy", [("auto", False), ("eligible", False), ("eligible", True)])
@pytest.mark.parametrize("input_mode", ["guided_binder", "condition_spec"])
def test_rfd4_remote_attempt_persists_peer_identity_without_starting_local_worker(
    tmp_path, monkeypatch, workstation_selection, busy, input_mode,
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs, workflow_queue, local_worker

    parent = tmp_path / "parent"
    target = parent / "artifacts/rfd4_proteina/input/PD-L1.pdb.cif"
    target.parent.mkdir(parents=True)
    target.write_text("data_target\n")
    params = {
        "task": "protein_binder", "input_mode": input_mode, "target_structure": str(target),
        "target_chains": ["A"], "num_designs": 1, "batch_size": 1, "nsteps": 60,
        "length_min": 20, "length_max": 30,
        "options": {"deterministic": True, "filter_conditions": True, "refold": False},
    }
    request = {"target_structure": "/work/artifacts/rfd4_proteina/input/PD-L1.pdb.cif"}
    if input_mode == "condition_spec":
        params.update(
            input_files=[str(target)],
            condition_spec="""idr_binder:
  input: PD-L1.pdb.cif
  contig: \"20-30, /0, A1-20\"
  conditions:
    C_SEQ:
      true:
        - select: \"sel('A')\"
""",
        )
        request["input_files"] = ["/work/artifacts/rfd4_proteina/input/PD-L1.pdb.cif"]
    jobs.write_json(parent / "input.json", {"params": params})
    jobs.write_json(parent / "artifacts/rfd4_proteina/request.json", request)
    jobs.write_json(parent / "metadata.json", {
        "studio_owner_id": "owner", "studio_project_id": "project", "studio_development_mode": True,
    })
    campaign_config = distributed_campaign.validate_campaign_config({
        "generator_id": "rfd4_proteina", "accepted_target": 1, "max_attempts": 1,
        "acceptance_level": "validated",
        "validation_stages": [{"workflow": "refolding.optimized_prediction", "required_metrics": ["iptm"]}],
        "workstation_selection": workstation_selection, "max_retries": 0,
    })
    ledger = CampaignLedger.create(parent / "ledger.json", campaign_config, campaign_id="remote-campaign")
    attempt = ledger.authorize(1, batch_id="b000001")[0]
    monkeypatch.setattr(distributed_campaign, "_peer_campaign_mode", lambda _config: True)

    child = tmp_path / "design" / "proxy-child"

    def enqueue(_function, *, task_group, tool, job_type, kwargs, **_options):
        child.mkdir(parents=True)
        jobs.write_json(child / "input.json", {"task_group": task_group, "tool": tool,
                                                "job_type": job_type, "params": kwargs})
        jobs.write_json(child / "metadata.json", {
            "status": "queued",
            "resource_request": {
                "cpu_cores": 8, "gpu": True, "gpu_device": "auto", "min_vram_gb": 16.0,
                "ram_gb": 0.0, "scratch_gb": 0.0,
            },
        })
        jobs.write_json(child / "worker_request.json", {"kind": "workflow_call"})
        jobs.write_json(child / "command.json", {"mode": "internal", "command": []})
        return child

    monkeypatch.setattr(workflow_queue, "enqueue_workflow_call", enqueue)
    monkeypatch.setattr(
        local_worker, "spawn_worker_for_run",
        lambda *_args, **_kwargs: pytest.fail("Remote attempts must not start a local worker."),
    )
    submissions = []

    def submit(request, manifest, submission_key, settings, *, owner_id, project_id, state, persist, allow_development):
        assert request["workflow"] == "design.rfd4_proteina"
        assert request["parameters"]["gpu_device"] == "auto"
        assert request["parameters"]["input_mode"] == input_mode
        if input_mode == "condition_spec":
            assert request["parameters"]["condition_spec"].startswith("idr_binder:")
            assert request["parameters"]["input_files"] == [str(target)]
        else:
            assert request["parameters"]["input_files"] == []
        assert "task_group" not in request["parameters"]
        assert request["resources"]["gpu_device"] == "auto"
        assert "gpu" not in request["resources"]
        if busy and settings["placement"] == "eligible":
            from mn_molecular_design_studio.cluster.candidate_batches import NoEligiblePeerAvailable
            raise NoEligiblePeerAvailable("peers busy")
        assert settings["placement"] == ("auto" if busy else workstation_selection)
        if workstation_selection == "eligible" and not busy:
            assert settings["placement_ordinal"] == 1
            assert settings["exclude_node_ids"] == []
        assert owner_id == "owner" and project_id == "project"
        assert allow_development is True
        assert submission_key == state["submission_key"]
        state["prepared_envelope"] = {"idempotency_key": submission_key}
        persist()
        submissions.append((manifest, submission_key, ledger.snapshot()["attempts"][attempt["attempt_id"]]))
        state["cluster_job_id"] = "cluster-job"
        persist()
        return {"cluster_job_id": "cluster-job"}

    monkeypatch.setattr("mn_molecular_design_studio.cluster.candidate_batches.submit", submit)
    distributed_campaign._create_rfd4_attempt_child(parent, ledger, attempt, campaign_config)

    latest = ledger.snapshot()["attempts"][attempt["attempt_id"]]
    assert latest["status"] == "queued"
    assert latest["cluster_job_id"] == "cluster-job"
    if busy:
        assert latest["cluster_submission"]["placement_fallback"] == "auto_queue"
    assert latest["cluster_submission"]["prepared_envelope"]["idempotency_key"] == latest["remote_submission_key"]
    assert submissions and submissions[0][2]["cluster_submission"]["prepared_envelope"]
    child_params = jobs.read_json(child / "input.json")["params"]
    assert child_params["input_mode"] == input_mode
    assert child_params["target_structure"] == str(target)
    metadata = jobs.read_json(child / "metadata.json")
    assert metadata["cluster_remote"] is True
    assert metadata["cluster_job_id"] == "cluster-job"


def test_candidate_batch_actor_can_restore_api_development_identity_in_scheduler_worker(monkeypatch):
    from mn_molecular_design_studio import auth
    from mn_molecular_design_studio.cluster import candidate_batches, service

    development_user = auth.development_user()
    monkeypatch.setattr(auth, "list_users", lambda: [])
    monkeypatch.setattr(auth, "is_development_mode", lambda: False)
    monkeypatch.setattr(auth, "require_project_access", lambda user, project: None)
    monkeypatch.setattr(service, "user_record", lambda user: user)

    with pytest.raises(PermissionError, match="account is unavailable"):
        candidate_batches.actor(development_user["id"], "project")
    restored = candidate_batches.actor(
        development_user["id"], "project", allow_development=True,
    )
    assert restored["id"] == development_user["id"]
    assert restored["role"] == "admin"


def test_dispatch_failure_marks_created_child_terminal(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    parent = tmp_path / "parent"
    parent.mkdir()
    jobs.write_json(parent / "input.json", {"params": {}})
    config = distributed_campaign.validate_campaign_config({
        "generator_id": "rfd4_proteina", "accepted_target": 1, "max_attempts": 1,
        "acceptance_level": "validated", "max_retries": 0,
        "validation_stages": [{
            "workflow": "refolding.optimized_prediction", "required_metrics": ["iptm"],
        }],
    })
    campaign = CampaignLedger.create(
        parent / "artifacts/distributed_campaign/ledger.json", config,
        campaign_id="dispatch-failure-campaign",
    )
    child = tmp_path / "design" / "failed-proxy"

    def create_child(_parent, ledger, attempt, _config):
        child.mkdir(parents=True)
        jobs.write_json(child / "metadata.json", {"status": "queued", "cluster_remote": True})
        ledger.transition(
            attempt["attempt_id"], "queued", child_run_dir=str(child), cluster_remote=True,
        )
        raise PermissionError("peer actor unavailable")

    monkeypatch.setattr(distributed_campaign, "_create_attempt_child", create_child)
    distributed_campaign._dispatch_attempts(
        parent, campaign, config, active_count=0, accepted_count=0,
    )

    attempt_id, attempt = next(iter(campaign.snapshot()["attempts"].items()))
    assert attempt["status"] == "failed"
    metadata = jobs.read_json(child / "metadata.json")
    assert metadata["status"] == "failed"
    assert "peer actor unavailable" in metadata["worker_error"]
    children, active, candidates = distributed_campaign._reconcile_attempts(campaign, parent, config)
    assert children[attempt_id] == child
    assert active == []
    assert candidates == []


def test_rfd4_remote_completion_materializes_verified_artifacts_and_provenance(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    from mn_molecular_design_studio.cluster import candidate_batches

    parent = tmp_path / "parent"
    parent.mkdir()
    jobs.write_json(parent / "metadata.json", {
        "studio_owner_id": "owner", "studio_project_id": "project", "studio_development_mode": True,
    })
    campaign = CampaignLedger.create(
        parent / "ledger.json",
        {"accepted_target": 1, "max_attempts": 1},
        campaign_id="remote-collection-campaign",
    )
    attempt = campaign.authorize(1, batch_id="batch")[0]
    campaign.transition(attempt["attempt_id"], "queued", cluster_job_id="cluster-job")
    child = tmp_path / "design" / "child"
    child.mkdir(parents=True)
    jobs.write_json(child / "metadata.json", {"status": "queued", "cluster_remote": True})
    manifest = tmp_path / "peer" / "submission.json"
    collected = manifest.with_suffix(".collected")
    artifacts = collected / "artifacts"
    rfd4 = artifacts / "rfd4_proteina"
    (rfd4 / "input").mkdir(parents=True)
    jobs.write_json(rfd4 / "request.json", {"target_structure": "/work/artifacts/rfd4_proteina/input/benchmark_target.pdb.cif"})
    (rfd4 / "input" / "benchmark_target.pdb.cif").write_text("target fixture\n")
    jobs.write_json(rfd4 / "native_command.json", {"argv": ["rfd4-generate"]})
    (rfd4 / "effective_config.yaml").write_text("seed: 19\n")
    jobs.write_json(rfd4 / "prepared_inputs.json", {"target_chain": "A"})
    jobs.write_json(collected / "metadata.json", {
        "status": "completed", "resource_allocation": {"gpu_device": "0", "gpu_ids": ["GPU-5090"]},
        "scheduler_worker_id": "worker-5090", "worker_started_at": "2026-10-28T10:00:00Z",
        "completed_at": "2026-10-28T10:02:00Z",
    })
    jobs.write_json(collected / "result.json", {"success": True, "metrics": {"checkpoint_revision": "checkpoint-1", "image": "sha256:image"}})
    normalized = artifacts / "normalized_candidates" / "candidates.jsonl"
    normalized.parent.mkdir(parents=True)
    normalized.write_text(json.dumps({"candidate_id": "native-1", "complex_structure_path": str(artifacts / "rfd4_proteina/output/complex.cif")}) + "\n")
    remote_structure = rfd4 / "output" / "complex.cif"
    remote_structure.parent.mkdir(parents=True)
    remote_structure.write_text("complex fixture\n")
    collected_normalized = artifacts / "normalized_candidates" / "collected.jsonl"
    collected_normalized.write_text(normalized.read_text())

    monkeypatch.setattr(candidate_batches, "actor", lambda *_args, **_kwargs: {"id": "owner"})
    monkeypatch.setattr(
        candidate_batches, "_artifact",
        lambda _peer, _job, relative, _actor, directory: Path(directory) / relative,
    )
    peer = {
        "manifest": str(manifest), "cluster_job_id": "cluster-job", "input_transfer_bytes": 42,
        "output_transfer_bytes": 84, "transfers": {"result.json": {"sha256": "abc", "bytes": 10}},
    }
    output = {"candidates_jsonl": str(collected_normalized)}
    remote_job = {"node_id": "node-5090", "node_name": "pc004032", "metadata": {}}

    assert distributed_campaign._materialize_rfd4_peer_output(
        parent, campaign, attempt, child, peer, output, remote_job,
    )
    local_rows = (child / "artifacts/normalized_candidates/candidates.jsonl").read_text()
    assert str(child / "artifacts/rfd4_proteina/output/complex.cif") in local_rows
    assert (child / "artifacts/rfd4_proteina/input/benchmark_target.pdb.cif").is_file()
    provenance = jobs.read_json(child / "artifacts/distributed_campaign/remote_execution.json")
    assert provenance["worker_node_id"] == "node-5090"
    assert provenance["resource_allocation"]["gpu_ids"] == ["GPU-5090"]
    assert provenance["checkpoint_revision"] == "checkpoint-1"
    assert provenance["input_transfer_bytes"] == 42
    assert jobs.read_json(child / "metadata.json")["scheduler_worker_id"] == "worker-5090"
    assert campaign.snapshot()["attempts"][attempt["attempt_id"]]["worker_node_id"] == "node-5090"


def test_rfd4_pending_artifact_command_is_persisted_and_reused(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    from mn_molecular_design_studio.cluster import candidate_batches

    parent = tmp_path / "parent"
    parent.mkdir()
    jobs.write_json(parent / "metadata.json", {
        "studio_owner_id": "owner", "studio_project_id": "project", "studio_development_mode": True,
    })
    campaign = CampaignLedger.create(
        parent / "artifacts/distributed_campaign/ledger.json",
        {"accepted_target": 1, "max_attempts": 1},
        campaign_id="persist-artifact-transfer-campaign",
    )
    attempt = campaign.authorize(1, batch_id="batch")[0]
    child = tmp_path / "design" / "child"
    child.mkdir(parents=True)
    manifest = tmp_path / "peer" / "submission.json"
    peer = {"manifest": str(manifest), "cluster_job_id": "cluster-job", "transfers": {}}
    campaign.transition(
        attempt["attempt_id"], "queued", child_run_dir=str(child),
        cluster_job_id="cluster-job", cluster_submission=peer,
    )
    monkeypatch.setattr(candidate_batches, "actor", lambda *_args, **_kwargs: {"id": "owner"})
    requests = []

    def pending_artifact(state, _job_id, relative, _user, _directory):
        transfer = state.setdefault("transfers", {}).setdefault(relative, {})
        if "command_id" not in transfer:
            transfer["command_id"] = f"transfer-{len(requests) + 1}"
            requests.append(relative)
        return None

    monkeypatch.setattr(candidate_batches, "_artifact", pending_artifact)
    output = {"candidates_jsonl": str(tmp_path / "not-needed.jsonl")}

    for _ in range(2):
        stored = campaign.snapshot()["attempts"][attempt["attempt_id"]]["cluster_submission"]
        assert not distributed_campaign._materialize_rfd4_peer_output(
            parent, campaign, attempt, child, stored, output, {"node_id": "peer-4032"},
        )

    stored = campaign.snapshot()["attempts"][attempt["attempt_id"]]["cluster_submission"]
    relative = "artifacts/rfd4_proteina/request.json"
    assert stored["transfers"][relative]["command_id"] == "transfer-1"
    assert requests == [
        relative,
        "artifacts/rfd4_proteina/native_command.json",
        "artifacts/rfd4_proteina/effective_config.yaml",
        "artifacts/rfd4_proteina/prepared_inputs.json",
    ]


def test_collection_does_not_hold_generation_slot_but_reserves_objective_result(tmp_path, monkeypatch):
    config = distributed_campaign.validate_campaign_config({
        "generator_id": "rfd4_proteina", "objective_type": "generated",
        "accepted_target": 10, "max_attempts": 20, "batch_size": 4,
        "max_concurrent_batches": 1, "max_concurrent_generation_jobs": 4,
    })
    campaign = CampaignLedger.create(tmp_path / "campaign.json", config, campaign_id="refill-campaign")
    launched = []
    monkeypatch.setattr(
        distributed_campaign, "_create_attempt_child",
        lambda _parent, _ledger, attempt, _config: launched.append(attempt),
    )

    active = [
        {"attempt_id": "a1", "worker_status": "running"},
        {"attempt_id": "a2", "worker_status": "running"},
        {"attempt_id": "a3", "worker_status": "queued"},
        {"attempt_id": "a4", "worker_status": "collecting"},
    ]
    assert distributed_campaign._generation_slot_count(active) == 3
    assert distributed_campaign._collecting_attempt_count(active) == 1
    distributed_campaign._dispatch_attempts(
        tmp_path, campaign, config,
        active_count=distributed_campaign._generation_slot_count(active),
        accepted_count=0,
        objective_reserved_count=distributed_campaign._collecting_attempt_count(active),
    )
    assert len(launched) == 1

    # The pending result still reserves the last generated objective unit.
    before = len(launched)
    distributed_campaign._dispatch_attempts(
        tmp_path, campaign, config, active_count=0, accepted_count=9,
        objective_reserved_count=1,
    )
    assert len(launched) == before


def test_runtime_cancellation_drains_completed_remote_result_collection(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    campaign = CampaignLedger.create(
        tmp_path / "campaign.json", {"accepted_target": 2, "max_attempts": 4},
        campaign_id="deadline-collection-campaign",
    )
    attempt = campaign.authorize(1, batch_id="batch")[0]
    child = tmp_path / "child"
    child.mkdir()
    jobs.write_json(child / "metadata.json", {
        "status": "queued", "cluster_remote": True, "result_collection_pending": True,
        "worker_started_at": "2026-10-29T00:00:00Z",
        "resource_allocation": {"gpu_ids": ["GPU-test"]},
    })
    campaign.transition(
        attempt["attempt_id"], "queued", child_run_dir=str(child),
        cluster_submission={"cluster_job_id": "already-completed"},
    )

    distributed_campaign._cancel_queued_attempts(
        campaign, {attempt["attempt_id"]: child}, "runtime budget reached",
    )

    assert campaign.snapshot()["attempts"][attempt["attempt_id"]]["status"] == "queued"
    assert distributed_campaign._campaign_compute_seconds(
        campaign, {attempt["attempt_id"]: child},
        [{**attempt, "worker_status": "collecting"}],
    ) == 0


def test_campaign_cancellation_uses_parent_ledger_and_preserves_completed_peer_job(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import candidate_batches
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    parent = tmp_path / "parent"
    parent.mkdir()
    jobs.write_json(parent / "metadata.json", {
        "studio_owner_id": "owner", "studio_project_id": "project", "studio_development_mode": True,
    })
    campaign = CampaignLedger.create(
        parent / "artifacts/distributed_campaign/ledger.json",
        {"accepted_target": 1, "max_attempts": 1},
        campaign_id="cancel-completed-peer-campaign",
    )
    attempt = campaign.authorize(1, batch_id="batch")[0]
    child = tmp_path / "design" / "child"
    child.mkdir(parents=True)
    jobs.write_json(child / "metadata.json", {"status": "queued"})
    peer = {"cluster_job_id": "completed-job"}
    campaign.transition(
        attempt["attempt_id"], "queued", child_run_dir=str(child), cluster_submission=peer,
    )

    class Store:
        def job(self, job_id, user):
            assert job_id == "completed-job"
            assert user == {"id": "owner"}
            return {"status": "completed"}

    monkeypatch.setattr(candidate_batches, "Store", Store)
    monkeypatch.setattr(candidate_batches, "actor", lambda *_args, **_kwargs: {"id": "owner"})
    monkeypatch.setattr(
        candidate_batches, "cancel",
        lambda *_args, **_kwargs: pytest.fail("A completed inference job must not be cancelled."),
    )

    distributed_campaign._cancel_campaign_attempt(campaign, attempt, child, "controller stopped")

    assert jobs.read_json(child / "metadata.json")["status"] == "cancelled"
    assert campaign.snapshot()["attempts"][attempt["attempt_id"]]["cluster_submission"] == peer


def test_failed_validation_stage_retries_same_run_within_campaign_limit(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    stage = {
        "workflow": "refolding.optimized_prediction",
        "required_metrics": ["iptm"],
        "rules": [{"field": "iptm", "operator": ">=", "value": 0.5}],
    }
    config = distributed_campaign.validate_campaign_config({
        "accepted_target": 1,
        "max_attempts": 1,
        "acceptance_level": "validated",
        "max_retries": 1,
        "validation_stages": [stage],
    })
    campaign = CampaignLedger.create(tmp_path / "campaign.json", config, campaign_id="retry-campaign")
    attempt = campaign.authorize(1, batch_id="batch")[0]
    candidate = {
        "candidate_id": "retry-campaign-a00000001-c1",
        "binder_sequence": "ACDE",
        "metrics": {"native_pass_filters": True},
    }
    campaign.record_receipt(
        attempt["attempt_id"], "native-receipt", candidates=[candidate],
        native_state={"trajectories": 1, "accepted": 1},
    )
    validation_run = tmp_path / "validation-run"
    validation_run.mkdir()
    jobs.write_json(validation_run / "metadata.json", {"status": "failed", "worker_error": "transient worker error"})
    jobs.write_json(validation_run / "worker_request.json", {"kind": "optimized_prediction_boltz2", "kwargs": {}})
    campaign.update_candidate(
        candidate["candidate_id"],
        validation_stages={"0": {"status": "failed", "run_dir": str(validation_run), "retry_count": 0}},
    )
    resumed = []

    def resume(run_dir, *, gpu_device=None):
        resumed.append((run_dir, gpu_device))
        metadata = jobs.read_json(run_dir / "metadata.json")
        metadata["status"] = "queued"
        jobs.write_json(run_dir / "metadata.json", metadata)

    monkeypatch.setattr(jobs, "resume_job", resume)
    distributed_campaign._drive_validation(tmp_path, campaign, [candidate], config)
    state = campaign.snapshot()["candidates"][candidate["candidate_id"]]
    assert resumed == [(validation_run, "auto")]
    assert state["validation_stages"]["0"]["run_dir"] == str(validation_run)
    assert state["validation_stages"]["0"]["retry_count"] == 1
    assert state["validated_acceptance"] == "pending"

    metadata = jobs.read_json(validation_run / "metadata.json")
    metadata["status"] = "failed"
    jobs.write_json(validation_run / "metadata.json", metadata)
    distributed_campaign._drive_validation(tmp_path, campaign, [candidate], config)
    state = campaign.snapshot()["candidates"][candidate["candidate_id"]]
    assert len(resumed) == 1
    assert state["validated_acceptance"] == "failed"


def test_requeued_validation_reopens_failed_acceptance_until_evidence_is_reconciled(tmp_path):
    config = distributed_campaign.validate_campaign_config({
        "accepted_target": 1,
        "max_attempts": 1,
        "acceptance_level": "validated",
        "validation_stages": [{
            "workflow": "refolding.optimized_prediction",
            "required_metrics": ["iptm"],
        }],
    })
    campaign = CampaignLedger.create(tmp_path / "campaign.json", config, campaign_id="requeue-campaign")
    attempt = campaign.authorize(1, batch_id="batch")[0]
    candidate = {
        "candidate_id": "requeue-candidate",
        "binder_sequence": "ACDE",
        "metrics": {"native_pass_filters": True},
    }
    campaign.record_receipt(
        attempt["attempt_id"], "native-receipt", candidates=[candidate],
        native_state={"trajectories": 1, "accepted": 1},
    )
    validation_run = tmp_path / "validation-run"
    validation_run.mkdir()
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json

    write_json(validation_run / "metadata.json", {"status": "running"})
    campaign.update_candidate(
        candidate["candidate_id"],
        validated_acceptance="failed",
        validation_decision_reason="previous submission failed",
        decision_reason="previous submission failed",
        validation_stages={"0": {"status": "queued", "run_dir": str(validation_run)}},
    )

    distributed_campaign._drive_validation(
        tmp_path, campaign, [candidate], config, dispatch=False,
    )

    state = campaign.snapshot()
    record = state["candidates"][candidate["candidate_id"]]
    assert record["validated_acceptance"] == "pending"
    assert record["validation_decision_reason"] == ""
    assert distributed_campaign._validation_pending_candidates(state, config) == 1


def test_validation_rules_are_stage_and_candidate_qualified():
    source = {
        "candidate_id": "dc_00000001_a", "stage": "generation.backbone_sequence",
        "source_tool": "bindcraft2", "binder_sequence": "ACDE", "metrics": {"native_pass_filters": True},
    }
    prediction = {
        "candidate_id": "prediction_boltz2_001", "stage": "complex_refolding",
        "source_tool": "optimized_prediction_boltz2", "binder_sequence": "ACDE",
        "parents": [source["candidate_id"]], "metrics": {"iptm": 0.72},
        "raw_metadata": {
            "source_candidate": source,
            "measurement_status": {"status": "measured"},
            "prediction_protocol": {"engine": "boltz2", "num_recycles": 10},
            "identity_check": {"sequence_chain_identity": "pass"},
        },
    }
    decision = distributed_campaign._stage_selectors(
        source, prediction,
        {"required_metrics": ["iptm"], "rules": [{"field": "iptm", "operator": ">=", "value": 0.5}]},
    )
    assert decision["status"] == "pass"
    assert decision["measurements"][
        "evidence/complex_refolding/optimized_prediction_boltz2/prediction_boltz2_001/iptm"
    ]["value"] == 0.72

    invalid = {**prediction, "metrics": {"iptm": float("nan")}}
    failed = distributed_campaign._stage_selectors(
        source, invalid, {"required_metrics": ["iptm"], "rules": []},
    )
    assert failed["status"] == "missing"


def test_optimized_prediction_validation_request_keeps_source_sequence_and_resource_contract(
    tmp_path, monkeypatch, design_source_factory
):
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates

    root, parent = _parent_run(tmp_path, monkeypatch, {})
    source_run, _source = design_source_factory(root / "design" / "source")
    row = read_candidates(source_run / "artifacts/normalized_candidates/candidates.jsonl")[0]
    row["raw_metadata"]["source_run_dir"] = str(source_run)
    payload, _resources = distributed_campaign._validation_request(
        parent, row,
        {
            "workflow": "refolding.optimized_prediction",
            "parameters": {"engine": "boltz2", "samples": 1, "num_recycles": 6},
            "resources": {"gpu_device": "0", "cpu_cores": 3},
        },
        0, "campaign-test", None,
    )
    assert payload["workflow"] == "refolding.optimized_prediction"
    assert payload["resources"] == {"cpu_cores": 3, "gpu_device": "0"}
    task = payload["parameters"]["tasks"][0]
    assert task["target_protein_sequences"] == "ACDEFG"
    assert task["binder_sequence"] == "KLMNP"
    assert task["source_candidate_id"] == row["candidate_id"]
    assert task["source_run_dir"] == str(parent)
    assert task["molecule_chain_options"] == [
        {"id": "A", "use_msa": True}, {"id": "B", "use_msa": False},
    ]


def test_validation_dispatch_respects_two_job_concurrency_limit(tmp_path, monkeypatch, design_source_factory):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    stage = {
        "workflow": "refolding.optimized_prediction",
        "parameters": {"engine": "boltz2", "samples": 1, "num_recycles": 1},
        "resources": {"gpu_device": "auto", "cpu_cores": 4, "min_vram_gb": 16},
        "required_metrics": ["iptm"],
    }
    config = distributed_campaign.validate_campaign_config({
        "generator_id": "rfd4_proteina",
        "acceptance_level": "validated",
        "accepted_target": 2,
        "max_attempts": 3,
        "max_concurrent_validations": 2,
        "validation_stages": [stage],
    })
    root, parent = _parent_run(tmp_path, monkeypatch, {
        "generator_id": "rfd4_proteina",
        "acceptance_level": "validated",
        "accepted_target": 2,
        "max_attempts": 3,
        "validation_stages": [stage],
    })
    campaign = CampaignLedger.create(parent / "artifacts/distributed_campaign/ledger.json", config)
    candidates = []
    for ordinal in range(1, 4):
        source_run, candidate = design_source_factory(root / "design" / f"source-{ordinal}")
        candidate["candidate_id"] = f"source-candidate-{ordinal}"
        candidate["raw_metadata"]["source_run_dir"] = str(source_run)
        attempt = campaign.authorize(1, batch_id=f"retained-{ordinal}")[0]
        campaign.transition(attempt["attempt_id"], "running")
        campaign.record_receipt(
            attempt["attempt_id"], f"source-receipt-{ordinal}", candidates=[candidate],
            native_state={"engine": "rfd4_proteina", "authoritative": False, "output_candidates": 1},
        )
        candidates.append(candidate)

    submitted = []

    def submit_validation(_parent, _ledger, candidate, submitted_stage, stage_index, _previous_run):
        assert submitted_stage["resources"]["gpu_device"] == "auto"
        run_dir = root / "refolding-validation" / f"validation-{candidate['candidate_id']}"
        run_dir.mkdir(parents=True)
        jobs.write_json(run_dir / "metadata.json", {
            "run_id": run_dir.name, "task_group": "refolding-validation", "status": "queued",
        })
        submitted.append((candidate["candidate_id"], stage_index, run_dir))
        return run_dir

    monkeypatch.setattr(distributed_campaign, "_submit_validation_stage", submit_validation)

    distributed_campaign._drive_validation(parent, campaign, candidates, config)
    first_state = campaign.snapshot()
    assert len(submitted) == 2
    assert distributed_campaign._validation_active_count(first_state) == 2
    assert distributed_campaign._validation_pending_candidates(first_state, config) == 3
    assert {candidate_id for candidate_id, _index, _path in submitted} == {
        "source-candidate-1", "source-candidate-2",
    }

    distributed_campaign._drive_validation(parent, campaign, candidates, config)
    assert len(submitted) == 2
    assert distributed_campaign._validation_active_count(campaign.snapshot()) == 2


def test_ipsae_stage_request_depends_on_its_preceding_prediction(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import candidates_jsonl_path

    parent = tmp_path / "campaign-parent"
    parent.mkdir()
    prediction_run = tmp_path / "prediction-run"
    prediction_run.mkdir()
    row = {"candidate_id": "dc-1", "binder_sequence": "KLMNP"}
    payload, _resources = distributed_campaign._validation_request(
        parent, row,
        {"workflow": "analysis.ipsae", "parameters": {"pae_cutoff": 8.0}},
        1, "campaign-test", prediction_run,
    )
    assert payload["workflow"] == "analysis.ipsae"
    assert payload["resources"]["gpu_device"] == "none"
    assert payload["parameters"]["source_run_dir"] == str(prediction_run)
    assert payload["parameters"]["candidates_jsonl"] == str(candidates_jsonl_path(prediction_run))


def test_global_diversity_deduplicates_sequences_across_attempts(tmp_path):
    campaign = ledger(tmp_path, max_attempts=2, target=2)
    candidates = []
    for ordinal in (1, 2):
        attempt = campaign.authorize(1, batch_id=f"b{ordinal}")[0]
        candidate_id = f"candidate-{ordinal}"
        row = {
            "candidate_id": candidate_id,
            "binder_sequence": "ACDEFG",
            "metrics": {"native_pass_filters": True},
        }
        campaign.record_receipt(
            attempt["attempt_id"], f"receipt-{ordinal}", candidates=[row],
            native_state={"trajectories": 1, "accepted": 1, "rejections": {}},
        )
        candidates.append(row)
    distributed_campaign._apply_diversity(
        tmp_path, campaign, candidates,
        {"acceptance_level": "diversity_qualified", "diversity": {"deduplicate_exact": True}},
    )
    assert campaign.counters()["diversity_qualified_accepted_candidates"] == 1
    assert campaign.snapshot()["candidates"]["candidate-2"]["diversity_acceptance"] == "rejected"


def test_controller_can_stop_on_validated_global_acceptance(
    tmp_path, monkeypatch, design_source_factory
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates

    config = {
        "accepted_target": 1,
        "acceptance_level": "validated",
        "max_attempts": 2,
        "batch_size": 1,
        "max_concurrent_batches": 2,
        "validation_stages": [{
            "workflow": "refolding.optimized_prediction",
            "parameters": {"engine": "boltz2"},
            "required_metrics": ["iptm"],
            "rules": [{"field": "iptm", "operator": ">=", "value": 0.5}],
        }],
        "poll_interval_seconds": 0.25,
        "max_retries": 0,
        "seed": 102,
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    _controlled_workers(monkeypatch, root, design_source_factory, accepted={1})
    _mock_finish(monkeypatch)

    def submit_validation(request_path):
        request = json.loads(request_path.read_text())
        source_id = request["parameters"]["tasks"][0]["source_candidate_id"]
        source = next(row for row in read_candidates(parent / "artifacts/normalized_candidates/candidates.jsonl")
                      if row["candidate_id"] == source_id)
        run_dir = root / "refolding-validation" / "fake-validation"
        run_dir.mkdir(parents=True)
        jobs.write_json(run_dir / "metadata.json", {
            "run_id": run_dir.name, "task_group": "refolding-validation", "job_type": "complex_refolding",
            "tool": "optimized_prediction_boltz2", "status": "completed",
            "resource_allocation": {"gpu_device": "none"},
        })
        prediction = {
            "candidate_id": "fake-prediction-1", "stage": "complex_refolding",
            "source_tool": "optimized_prediction_boltz2", "binder_sequence": source["binder_sequence"],
            "target_chains": source["target_chains"], "binder_chains": source["binder_chains"],
            "complex_pdb": source["complex_pdb"], "parents": [source_id], "metrics": {"iptm": 0.8},
            "raw_metadata": {
                "source_candidate": source, "measurement_status": {"status": "measured"},
                "prediction_protocol": {"engine": "boltz2", "num_recycles": 6},
                "identity_check": {"sequence_chain_identity": "pass"},
            },
        }
        write_candidates(run_dir, "optimized_prediction_boltz2", [prediction])
        jobs.write_json(run_dir / "result.json", {"success": True})
        return {"job": {"run_dir": str(run_dir)}}

    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit_validation)
    distributed_campaign.run_bindcraft2_distributed_campaign(parent)

    result = jobs.read_json(parent / "result.json")
    assert result["metrics"]["termination_reason"] == "target_reached"
    assert result["metrics"]["native_accepted_candidates"] == 1
    assert result["metrics"]["validated_accepted_candidates"] == 1
    final_candidate = read_candidates(parent / "artifacts/normalized_candidates/candidates.jsonl")[0]
    assert final_candidate["raw_metadata"]["evidence_candidates"][0]["candidate_id"] == "fake-prediction-1"


def test_controller_continues_generation_after_validated_candidate_is_rejected(
    tmp_path, monkeypatch, design_source_factory
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation

    config = {
        "accepted_target": 1,
        "acceptance_level": "validated",
        "max_attempts": 2,
        "batch_size": 1,
        "max_concurrent_batches": 1,
        "validation_stages": [{
            "workflow": "refolding.optimized_prediction",
            "parameters": {"engine": "boltz2"},
            "required_metrics": ["iptm"],
            "rules": [{"field": "iptm", "operator": ">=", "value": 0.5}],
        }],
        "poll_interval_seconds": 0.25,
        "max_retries": 0,
        "seed": 104,
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    made = _controlled_workers(monkeypatch, root, design_source_factory, accepted={1, 2})
    _mock_finish(monkeypatch)
    validation_attempts = []

    def submit_validation(request_path):
        request = json.loads(request_path.read_text())
        source_id = request["parameters"]["tasks"][0]["source_candidate_id"]
        source = next(
            row for row in read_candidates(parent / "artifacts/normalized_candidates/candidates.jsonl")
            if row["candidate_id"] == source_id
        )
        attempt_id = source["raw_metadata"]["distributed_campaign"]["attempt_id"]
        validation_attempts.append(attempt_id)
        ordinal = int(attempt_id.rsplit("a", 1)[1])
        run_dir = root / "refolding-validation" / f"validation-{ordinal}"
        run_dir.mkdir(parents=True)
        jobs.write_json(run_dir / "metadata.json", {
            "run_id": run_dir.name, "task_group": "refolding-validation",
            "job_type": "complex_refolding", "tool": "optimized_prediction_boltz2",
            "status": "completed", "resource_allocation": {"gpu_device": "none"},
        })
        prediction = {
            "candidate_id": f"prediction-{ordinal}", "stage": "complex_refolding",
            "source_tool": "optimized_prediction_boltz2", "binder_sequence": source["binder_sequence"],
            "target_chains": source["target_chains"], "binder_chains": source["binder_chains"],
            "complex_pdb": source["complex_pdb"], "parents": [source_id],
            "metrics": {"iptm": 0.2 if ordinal == 1 else 0.8},
            "raw_metadata": {
                "source_candidate": source, "measurement_status": {"status": "measured"},
                "prediction_protocol": {"engine": "boltz2", "num_recycles": 6},
                "identity_check": {"sequence_chain_identity": "pass"},
            },
        }
        write_candidates(run_dir, "optimized_prediction_boltz2", [prediction])
        jobs.write_json(run_dir / "result.json", {"success": True})
        return {"job": {"run_dir": str(run_dir)}}

    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit_validation)
    distributed_campaign.run_bindcraft2_distributed_campaign(parent)

    result = jobs.read_json(parent / "result.json")
    campaign = CampaignLedger(parent / "artifacts/distributed_campaign/ledger.json")
    state = campaign.snapshot()
    rows = read_candidates(parent / "artifacts/normalized_candidates/candidates.jsonl")
    assert result["metrics"]["termination_reason"] == "target_reached"
    assert result["metrics"]["validated_accepted_candidates"] == 1
    assert result["metrics"]["attempts_authorized"] == 2
    assert result["metrics"]["rejected_candidates"] == 1
    assert len(made) == 2
    assert validation_attempts == [made[0][0], made[1][0]]
    candidates = state["candidates"]
    rejected = next(candidate for candidate in candidates.values() if candidate["validated_acceptance"] == "rejected")
    accepted = next(candidate for candidate in candidates.values() if candidate["validated_acceptance"] == "accepted")
    assert rejected["validation_stages"]["0"]["selection"]["status"] == "reject"
    assert rejected["validation_stages"]["0"]["evidence_candidate"]["metrics"]["iptm"] == 0.2
    assert accepted["validation_stages"]["0"]["selection"]["status"] == "pass"
    assert accepted["validation_stages"]["0"]["evidence_candidate"]["metrics"]["iptm"] == 0.8
    by_id = {row["candidate_id"]: row for row in rows}
    assert by_id[rejected["candidate_id"]]["raw_metadata"]["distributed_campaign"]["validated_acceptance"] == "rejected"
    assert by_id[accepted["candidate_id"]]["raw_metadata"]["distributed_campaign"]["validated_acceptance"] == "accepted"
    for row in (by_id[rejected["candidate_id"]], by_id[accepted["candidate_id"]]):
        evidence = row["raw_metadata"]["evidence_candidates"][0]
        validation = evidence["raw_metadata"]["distributed_campaign_validation"]
        assert validation == {
            "campaign_id": state["campaign_id"],
            "attempt_id": row["raw_metadata"]["distributed_campaign"]["attempt_id"],
            "stage_index": 0,
            "workflow": "refolding.optimized_prediction",
            "status": "completed",
        }


def test_late_validation_recovery_promotes_budget_reason_after_target_is_met(
    tmp_path, monkeypatch, design_source_factory
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation

    config = {
        "accepted_target": 1,
        "acceptance_level": "validated",
        "max_attempts": 1,
        "batch_size": 1,
        "max_concurrent_batches": 1,
        "validation_stages": [{
            "workflow": "refolding.optimized_prediction",
            "parameters": {"engine": "boltz2"},
            "required_metrics": ["iptm"],
            "rules": [{"field": "iptm", "operator": ">=", "value": 0.5}],
        }],
        "poll_interval_seconds": 0.25,
        "max_retries": 0,
        "seed": 103,
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    _controlled_workers(monkeypatch, root, design_source_factory, accepted={1})
    _mock_finish(monkeypatch)
    validation_run = root / "refolding-validation" / "late-validation"

    def submit_validation(request_path):
        request = json.loads(request_path.read_text())
        source_id = request["parameters"]["tasks"][0]["source_candidate_id"]
        source = next(
            row for row in read_candidates(parent / "artifacts/normalized_candidates/candidates.jsonl")
            if row["candidate_id"] == source_id
        )
        validation_run.mkdir(parents=True)
        jobs.write_json(validation_run / "metadata.json", {
            "run_id": validation_run.name, "task_group": "refolding-validation",
            "job_type": "complex_refolding", "tool": "optimized_prediction_boltz2",
            "status": "running", "resource_allocation": {"gpu_device": "none"},
        })
        jobs.write_json(validation_run / "result.json", {})
        return {"job": {"run_dir": str(validation_run)}}

    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit_validation)

    class StopAfterPendingValidation(Exception):
        pass

    monkeypatch.setattr(distributed_campaign.time, "sleep", lambda _seconds: (_ for _ in ()).throw(StopAfterPendingValidation()))
    with pytest.raises(StopAfterPendingValidation):
        distributed_campaign._run_distributed_campaign(parent, "bindcraft2")

    campaign = CampaignLedger(parent / "artifacts/distributed_campaign/ledger.json")
    pending_state = campaign.snapshot()
    assert pending_state["termination_reason"] is None
    assert pending_state["generation_budget_exhausted"] is True
    pending_progress = jobs.read_json(parent / "metadata.json")["campaign_progress"]
    assert pending_progress["attempts_remaining"] == 0
    assert pending_progress["generation_budget_exhausted"] is True
    assert pending_progress["validation_pending_candidates"] == 1

    # Simulate recovery from a legacy controller that finalized before validation settled.
    campaign.terminate("budget_exhausted", detail="legacy premature finalization")
    source = read_candidates(parent / "artifacts/normalized_candidates/candidates.jsonl")[0]
    prediction = {
        "candidate_id": "late-validation-prediction", "stage": "complex_refolding",
        "source_tool": "optimized_prediction_boltz2", "binder_sequence": source["binder_sequence"],
        "target_chains": source["target_chains"], "binder_chains": source["binder_chains"],
        "complex_pdb": source["complex_pdb"], "parents": [source["candidate_id"]],
        "metrics": {"iptm": 0.8},
        "raw_metadata": {
            "source_candidate": source, "measurement_status": {"status": "measured"},
            "prediction_protocol": {"engine": "boltz2", "num_recycles": 6},
            "identity_check": {"sequence_chain_identity": "pass"},
        },
    }
    jobs.write_json(validation_run / "metadata.json", {
        **jobs.read_json(validation_run / "metadata.json"), "status": "completed",
    })
    write_candidates(validation_run, "optimized_prediction_boltz2", [prediction])
    jobs.write_json(validation_run / "result.json", {"success": True})
    monkeypatch.setattr(distributed_campaign.time, "sleep", lambda _seconds: None)

    distributed_campaign.run_bindcraft2_distributed_campaign(parent)

    final_state = campaign.snapshot()
    result = jobs.read_json(parent / "result.json")
    assert final_state["termination_reason"] == "target_reached"
    assert final_state["generation_budget_exhausted"] is True
    assert result["metrics"]["termination_reason"] == "target_reached"
    assert result["metrics"]["target_met"] is True
    assert result["metrics"]["generation_budget_exhausted"] is True
    assert result["metrics"]["validated_accepted_candidates"] == 1
    summary = jobs.read_json(parent / "artifacts/distributed_campaign/summary.json")
    assert summary["termination_reason"] == "target_reached"
    assert summary["generation_budget_exhausted"] is True


def _parent_run(tmp_path, monkeypatch, config):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    from mn_molecular_design_studio.engine.workflow_engine import runtime

    root = tmp_path / "runs"
    monkeypatch.setattr(jobs, "runs_root", lambda: root)
    monkeypatch.setattr(runtime, "runs_root", lambda: root)
    parent = root / "design" / "campaign-parent"
    (parent / "artifacts/raw/bindcraft2").mkdir(parents=True)
    jobs.write_json(parent / "metadata.json", {
        "run_id": parent.name, "task_group": "design", "job_type": "bindcraft2_campaign",
        "tool": "bindcraft2_campaign", "status": "running", "campaign_name": "test-campaign",
    })
    jobs.write_json(parent / "input.json", {
        "job_type": "bindcraft2_campaign",
        "tool": "bindcraft2_campaign",
        "inputs": {},
        "params": {
            "campaign_name": "test-campaign",
            "target_type": "structure",
            "target_file": "target.pdb",
            "target_chains": ["A"],
            "cpu_cores": 2,
            "distributed_campaign": {"enabled": True, **config},
        },
    })
    jobs.write_json(parent / "artifacts/raw/bindcraft2/settings.json", {"campaign_seed": 11})
    return root, parent


def test_campaign_summary_reads_remote_child_metadata_without_resolving_cluster_paths(
    tmp_path, monkeypatch,
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    root, parent = _parent_run(tmp_path, monkeypatch, {
        "generator_id": "rfd4_proteina", "accepted_target": 1, "max_attempts": 1,
        "objective_type": "generated", "acceptance_level": "native",
        "remote_artifact_policy": "defer_generated",
    })
    config = distributed_campaign.validate_campaign_config({
        "generator_id": "rfd4_proteina", "accepted_target": 1, "max_attempts": 1,
        "objective_type": "generated", "acceptance_level": "native",
        "remote_artifact_policy": "defer_generated",
    })
    campaign = CampaignLedger.create(
        parent / "artifacts/distributed_campaign/ledger.json", config,
    )
    attempt = campaign.authorize(1, batch_id="remote-batch")[0]
    child = root / "design" / "remote-child"
    child.mkdir(parents=True)
    (child / "metadata.json").write_text(json.dumps({
        "status": "completed", "scheduler_worker_id": "worker:0",
        "resource_allocation": {"gpu_device": "0"},
        "run_dir": "cluster://remote-node/design/remote-run",
    }))
    campaign.transition(
        attempt["attempt_id"], "completed", child_run_dir=str(child),
        worker_node_id="remote-node", resource_allocation={"gpu_device": "0"},
    )

    distributed_campaign._sync_parent_results(parent, campaign, {attempt["attempt_id"]: child}, [])

    summary = jobs.read_json(parent / "artifacts/distributed_campaign/summary.json")
    assert summary["attempts"][0]["worker"] == "worker:0"
    assert summary["attempts"][0]["resource_allocation"] == {"gpu_device": "0"}


def _controlled_workers(monkeypatch, root, design_source_factory, *, accepted):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    made = []

    def create_child(parent, ledger, attempt, config):
        child = root / "design" / f"controlled-{attempt['ordinal']:08d}"
        passed = attempt["ordinal"] in accepted
        design_source_factory(child, passed=passed)
        metadata = jobs.read_json(child / "metadata.json")
        metadata.update(
            distributed_campaign_id=ledger.snapshot()["campaign_id"],
            distributed_attempt_id=attempt["attempt_id"],
            distributed_batch_id=attempt["batch_id"],
            resource_allocation={"gpu_device": distributed_campaign._gpu_for_attempt(config, attempt["ordinal"])},
        )
        jobs.write_json(child / "metadata.json", metadata)
        params = {
            "target_chains": ["A"],
            "gpu_device": distributed_campaign._gpu_for_attempt(config, attempt["ordinal"]),
        }
        jobs.write_json(child / "input.json", {"tool": "bindcraft2", "params": params})
        native = child / "artifacts/raw/bindcraft2/output/.campaign_state.json"
        native.parent.mkdir(parents=True, exist_ok=True)
        native.write_text(json.dumps({
            "trajectories": 1,
            "accepted": int(passed),
            "rejections": {"terminated": {"screen": int(not passed)}},
        }))
        ledger.transition(attempt["attempt_id"], "queued", child_run_dir=str(child), trajectory_index=1)
        made.append((attempt["attempt_id"], attempt["seed"], params["gpu_device"]))
        return child

    monkeypatch.setattr(distributed_campaign, "_create_attempt_child", create_child)
    return made


def _mock_finish(monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    def finish(run_dir, success, result):
        jobs.write_json(run_dir / "result.json", {"success": success, **result})
        jobs.update_status(run_dir, "completed" if success else "failed")

    monkeypatch.setattr(jobs, "finish_job", finish)


def test_attempt_child_rebases_target_and_native_settings_to_its_run(tmp_path, monkeypatch):
    import json

    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    config = distributed_campaign.validate_campaign_config({
        "accepted_target": 1,
        "max_attempts": 1,
        "seed": 44,
        "gpu_selection": "0",
        "resources": {"cpu_cores": 2},
    })
    _root, parent = _parent_run(tmp_path, monkeypatch, config)
    parent_raw = parent / "artifacts/raw/bindcraft2"
    (parent_raw / "target.pdb").write_text("parent target")
    jobs.write_json(parent / "input.json", {
        "tool": "bindcraft2_campaign",
        "inputs": {
            "target_type": "structure",
            "target_structure": "artifacts/raw/bindcraft2/target.pdb",
            "target_pdb": "artifacts/raw/bindcraft2/target.pdb",
            "target_chains": ["A"],
        },
        "params": {
            "campaign_name": "test-campaign",
            "target_type": "structure",
            "target_file": "target.pdb",
            "target_chains": ["A"],
            "cpu_cores": 2,
            "distributed_campaign": {"enabled": True, **config},
        },
    })
    jobs.write_json(parent_raw / "settings.json", {
        "campaign_seed": 44,
        "project_folder": "output",
        "targets": [{"name": "app_target", "target_path": "target.pdb", "chains": "A"}],
    })
    ledger_path = parent / "artifacts/distributed_campaign/ledger.json"
    campaign = CampaignLedger.create(ledger_path, config, campaign_id="localized-child")
    attempt = campaign.authorize(1, batch_id="batch-1")[0]
    monkeypatch.setattr(
        "mn_molecular_design_studio.engine.workflow_engine.core.local_worker.spawn_worker_for_run",
        lambda _run_dir: None,
    )

    child = distributed_campaign._create_attempt_child(parent, campaign, attempt, config)

    assert (child / "artifacts/raw/bindcraft2/target.pdb").read_text() == "parent target"
    child_input = jobs.read_json(child / "input.json")
    assert child_input["inputs"]["target_pdb"] == str(child / "artifacts/raw/bindcraft2/target.pdb")
    assert child_input["inputs"]["target_structure"] == str(child / "artifacts/raw/bindcraft2/target.pdb")
    assert Path(child_input["params"]["target_file"]).name == "target.pdb"
    settings = json.loads((child / "artifacts/raw/bindcraft2/settings.json").read_text())
    assert settings["project_folder"] == "output"
    assert settings["targets"][0]["target_path"] == "target.pdb"
    assert settings["campaign_seed"] == attempt["seed"]


def test_rfd4_attempt_child_queues_one_seeded_gpu_co_design(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    root = tmp_path / "runs"
    monkeypatch.setattr(jobs, "runs_root", lambda: root)
    parent = root / "design" / "rfd4-parent"
    staged_target = parent / "artifacts/rfd4_proteina/input/target.cif"
    staged_target.parent.mkdir(parents=True)
    staged_target.write_text("target")
    jobs.write_json(parent / "metadata.json", {
        "run_id": parent.name, "task_group": "design", "status": "running",
        "tool": "rfd4_proteina", "job_type": "rfd4_proteina_campaign",
    })
    jobs.write_json(parent / "input.json", {
        "tool": "rfd4_proteina", "inputs": {}, "params": {
            "task": "protein_binder", "input_mode": "guided_binder",
            "target_structure": str(staged_target), "target_chains": ["A"],
            "hotspots": "A40", "length_min": 60, "length_max": 120,
            "num_designs": 8, "batch_size": 1, "nsteps": 60, "seed": 55,
            "options": {"refold": False}, "cpu_cores": 3, "gpu_device": "0",
            "job_name": "PD-L1 campaign", "task_group": "design",
            "distributed_campaign": {"enabled": True, "generator_id": "rfd4_proteina"},
        },
    })
    jobs.write_json(parent / "artifacts/rfd4_proteina/request.json", {
        "target_structure": "/work/artifacts/rfd4_proteina/input/target.cif",
    })
    config = distributed_campaign.validate_campaign_config({
        "generator_id": "rfd4_proteina", "acceptance_level": "validated",
        "max_attempts": 1, "accepted_target": 1, "seed": 101,
        "gpu_selection": "0,1", "resources": {"cpu_cores": 3},
        "validation_stages": [{"workflow": "refolding.optimized_prediction",
                                "required_metrics": ["model_1_iptm"]}],
    })
    ledger = CampaignLedger.create(parent / "artifacts/distributed_campaign/ledger.json", config)
    attempt = ledger.authorize(1, batch_id="batch-rfd4")[0]
    dispatched = []
    monkeypatch.setattr(
        "mn_molecular_design_studio.engine.workflow_engine.core.local_worker.spawn_worker_for_run",
        lambda run_dir: dispatched.append(Path(run_dir)),
    )

    child = distributed_campaign._create_attempt_child(parent, ledger, attempt, config)

    child_params = jobs.read_json(child / "input.json")["params"]
    worker = jobs.read_json(child / "worker_request.json")
    metadata = jobs.read_json(child / "metadata.json")
    assert child_params["num_designs"] == child_params["batch_size"] == 1
    assert child_params["seed"] == attempt["seed"]
    assert child_params["target_structure"] == str(staged_target)
    assert child_params["distributed_campaign"] is None
    assert child_params["task"] == "protein_binder"
    assert worker["resources"]["gpu_device"] == "0"
    assert worker["resources"]["min_vram_gb"] == 16.0
    assert worker["resources"]["cpu_cores"] == 3
    assert metadata["distributed_attempt_id"] == attempt["attempt_id"]
    assert metadata["parent_role"] == "distributed_rfd4_attempt"
    assert dispatched == [child]
    assert ledger.snapshot()["attempts"][attempt["attempt_id"]]["requested_outputs"] == 1


def test_rfd4_candidate_reconciliation_preserves_co_designed_sequence_and_target_mapping(
    tmp_path, design_source_factory
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    child, candidate = design_source_factory(tmp_path / "rfd4-child")
    target_path = child / "target.cif"
    target_path.write_bytes((child / "designed.cif").read_bytes())
    jobs.write_json(child / "input.json", {"params": {
        "target_structure": str(target_path), "target_chains": ["A"],
        "length_min": 5, "length_max": 20,
    }})
    jobs.write_json(child / "metadata.json", {
        "resource_allocation": {"gpu_device": "1", "gpu_ids": [1]},
        "scheduler_worker_id": "local-test-worker",
    })
    candidate.update(
        source_tool="rfd4_proteina", target_pdb=str(target_path),
        raw_metadata={"engine_version": "test"},
    )
    (child / "artifacts/normalized_candidates/candidates.jsonl").write_text(
        json.dumps(candidate) + "\n"
    )
    attempt = {"ordinal": 1, "attempt_id": "a00000001", "batch_id": "batch-1", "seed": 4242}
    native_state = {
        "engine": "rfd4_proteina", "authoritative": False,
        "attempt_unit": "one native design request", "requested_designs": 1,
        "output_candidates": 1, "seed": 4242,
    }

    rows = distributed_campaign._rfd4_candidate_rows(child, "campaign-rfd4", attempt, native_state)

    assert len(rows) == 1
    row = rows[0]
    lineage = row["raw_metadata"]["distributed_campaign"]
    mapping = row["raw_metadata"]["source_to_prediction_components"]
    assert row["source_tool"] == "rfd4_proteina"
    assert row["binder_sequence"] == "KLMNP"
    assert row["binder_chains"] == ["B"] and row["target_chains"] == ["A"]
    assert mapping["target"] == [{
        "source_chain": "A", "prediction_chain": "A",
        "source_sequence": "ACDEFG", "prediction_sequence": "ACDEFG",
    }]
    assert lineage["native_acceptance"] == "unavailable"
    assert lineage["attempt_seed"] == 4242
    assert lineage["sequence_structure_mode"] == "co_design"


def test_rfd4_candidate_reconciliation_maps_mmcif_author_chain_to_label_chain(
    tmp_path, design_source_factory
):
    from Bio.PDB import MMCIFIO
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    child, candidate = design_source_factory(tmp_path / "rfd4-numeric-chain-child")
    target_path = child / "target.cif"
    payload = MMCIF2Dict(str(child / "designed.cif"))
    payload["_atom_site.auth_asym_id"] = [
        "1" if str(value) == "A" else str(value)
        for value in payload["_atom_site.auth_asym_id"]
    ]
    writer = MMCIFIO()
    writer.set_dict(payload)
    writer.save(str(target_path))
    jobs.write_json(child / "input.json", {"params": {
        "target_structure": str(target_path), "target_chains": ["1"],
        "length_min": 5, "length_max": 20,
    }})
    jobs.write_json(child / "metadata.json", {"resource_allocation": {"gpu_device": "0"}})
    candidate.update(source_tool="rfd4_proteina", raw_metadata={})
    (child / "artifacts/normalized_candidates/candidates.jsonl").write_text(
        json.dumps(candidate) + "\n"
    )

    rows = distributed_campaign._rfd4_candidate_rows(
        child, "campaign-rfd4-labels",
        {"ordinal": 1, "attempt_id": "a00000001", "batch_id": "b1", "seed": 9},
        {"authoritative": False, "output_candidates": 1},
    )

    assert rows[0]["raw_metadata"]["source_to_prediction_components"]["target"] == [{
        "source_chain": "1", "prediction_chain": "A",
        "source_sequence": "ACDEFG", "prediction_sequence": "ACDEFG",
    }]


def test_rfd4_campaign_stop_routes_through_global_child_cancellation(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler, jobs

    run_dir = tmp_path / "design" / "campaign-parent"
    run_dir.mkdir(parents=True)
    jobs.write_json(run_dir / "metadata.json", {
        "run_id": run_dir.name, "task_group": "design", "job_type": "rfd4_proteina_campaign",
        "tool": "rfd4_proteina", "status": "running",
    })
    jobs.write_json(run_dir / "result.json", {"success": True})
    calls = []
    monkeypatch.setattr(
        distributed_campaign, "cancel_distributed_campaign",
        lambda path, reason: calls.append((Path(path), reason)),
    )
    monkeypatch.setattr(central_scheduler, "cancel_scheduled_run", lambda _path: None)
    monkeypatch.setattr(jobs, "_signal_worker_process", lambda *_args: None)
    monkeypatch.setattr(jobs, "_stop_run_docker_containers", lambda _path: [])
    monkeypatch.setattr(jobs, "_release_resource_lock", lambda _path: None)

    jobs.stop_job(run_dir, "cancel RFD4 campaign")

    assert calls == [(run_dir, "cancel RFD4 campaign")]
    assert jobs.read_json(run_dir / "metadata.json")["status"] == "cancelled"


def test_controller_reaches_one_global_target_and_keeps_unique_seed_and_gpu_lineage(
    tmp_path, monkeypatch, design_source_factory
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    config = {
        "accepted_target": 2,
        "acceptance_level": "native",
        "max_attempts": 5,
        "batch_size": 2,
        "max_concurrent_batches": 2,
        "gpu_selection": "0,1",
        "poll_interval_seconds": 0.25,
        "max_retries": 0,
        "seed": 1028,
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    made = _controlled_workers(monkeypatch, root, design_source_factory, accepted={1, 2})
    _mock_finish(monkeypatch)

    distributed_campaign.run_bindcraft2_distributed_campaign(parent)

    result = jobs.read_json(parent / "result.json")
    assert result["success"] is True
    assert result["metrics"]["termination_reason"] == "target_reached"
    assert result["metrics"]["native_accepted_candidates"] == 2
    assert result["metrics"]["attempts_authorized"] == 2
    assert len({seed for _attempt, seed, _gpu in made}) == 2
    assert {gpu for _attempt, _seed, gpu in made} == {"0", "1"}
    rows = (parent / "artifacts/normalized_candidates/candidates.jsonl").read_text().splitlines()
    assert len(rows) == 2
    campaign_lineage = [json.loads(row)["raw_metadata"]["distributed_campaign"] for row in rows]
    assert {row["attempt_id"] for row in campaign_lineage} == {attempt for attempt, _seed, _gpu in made}
    assert all(row["native_acceptance"] == "accepted" for row in campaign_lineage)


def test_explicit_generation_job_cap_bounds_dispatched_attempts(tmp_path, monkeypatch):
    config = {
        "accepted_target": 6,
        "max_attempts": 6,
        "batch_size": 2,
        "max_concurrent_batches": 3,
        "max_concurrent_generation_jobs": 2,
        "poll_interval_seconds": 0.25,
        "max_retries": 0,
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    normalized = distributed_campaign.validate_campaign_config(config)
    campaign = CampaignLedger.create(
        parent / "artifacts/distributed_campaign/ledger.json", normalized,
    )
    children = []

    def create_child(_parent, ledger, attempt, _config):
        child = root / "design" / f"capped-{attempt['ordinal']:08d}"
        child.mkdir(parents=True)
        ledger.transition(attempt["attempt_id"], "queued", child_run_dir=str(child))
        children.append(child)
        return child

    monkeypatch.setattr(distributed_campaign, "_create_attempt_child", create_child)
    distributed_campaign._dispatch_attempts(
        parent, campaign, normalized, active_count=0, accepted_count=0,
    )

    attempts = campaign.snapshot()["attempts"]
    assert len(attempts) == len(children) == 2
    assert {attempt["batch_id"] for attempt in attempts.values()} == {"b000001"}


def test_controller_continues_after_native_rejections_then_reports_budget_exhausted(
    tmp_path, monkeypatch, design_source_factory
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    config = {
        "accepted_target": 1,
        "acceptance_level": "native",
        "max_attempts": 2,
        "batch_size": 2,
        "max_concurrent_batches": 2,
        "poll_interval_seconds": 0.25,
        "max_retries": 0,
        "seed": 9,
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    made = _controlled_workers(monkeypatch, root, design_source_factory, accepted=set())
    _mock_finish(monkeypatch)

    distributed_campaign.run_bindcraft2_distributed_campaign(parent)

    result = jobs.read_json(parent / "result.json")
    assert result["success"] is True
    assert result["metrics"]["campaign_success"] is False
    assert result["metrics"]["termination_reason"] == "budget_exhausted"
    assert result["metrics"]["native_accepted_candidates"] == 0
    assert result["metrics"]["attempts_authorized"] == 2
    assert result["metrics"]["attempts_completed"] == 2
    assert result["metrics"]["rejected_candidates"] == 2
    assert result["metrics"]["generation_budget_exhausted"] is True
    assert result["metrics"]["attempts_remaining"] == 0
    assert len(made) == 2


def test_controller_generates_more_after_global_diversity_rejection(
    tmp_path, monkeypatch, design_source_factory
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    config = {
        "accepted_target": 2,
        "acceptance_level": "diversity_qualified",
        "max_attempts": 3,
        "batch_size": 2,
        "max_concurrent_batches": 1,
        "poll_interval_seconds": 0.25,
        "max_retries": 0,
        "seed": 77,
        "diversity": {"deduplicate_exact": True, "max_per_cluster": 1},
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    made = _controlled_workers(monkeypatch, root, design_source_factory, accepted={1, 2, 3})
    _mock_finish(monkeypatch)

    distributed_campaign.run_bindcraft2_distributed_campaign(parent)

    result = jobs.read_json(parent / "result.json")
    assert result["metrics"]["termination_reason"] == "budget_exhausted"
    assert result["metrics"]["native_accepted_candidates"] == 3
    assert result["metrics"]["diversity_qualified_accepted_candidates"] == 1
    assert result["metrics"]["attempts_authorized"] == 3
    assert len(made) == 3


def test_controller_waits_for_newly_dispatched_child_before_budget_exhaustion(
    tmp_path, monkeypatch
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    config = {
        "accepted_target": 2,
        "acceptance_level": "native",
        "max_attempts": 2,
        "batch_size": 1,
        "max_concurrent_batches": 2,
        "gpu_selection": "0,1",
        "poll_interval_seconds": 0.25,
        "max_retries": 0,
        "seed": 88,
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    children = {}
    sleeps = []

    def create_child(_parent, campaign, attempt, _config):
        child = root / "design" / f"queued-{attempt['ordinal']:08d}"
        (child / "artifacts/raw/bindcraft2/output").mkdir(parents=True)
        jobs.write_json(child / "metadata.json", {
            "run_id": child.name, "task_group": "design", "job_type": "bindcraft2_design",
            "tool": "bindcraft2", "status": "queued",
            "resource_allocation": {"gpu_device": str(attempt["ordinal"] - 1)},
        })
        jobs.write_json(child / "input.json", {"params": {"gpu_device": str(attempt["ordinal"] - 1)}})
        jobs.write_json(child / "result.json", {})
        (child / "artifacts/normalized_candidates").mkdir(parents=True)
        (child / "artifacts/normalized_candidates/candidates.jsonl").write_text("")
        campaign.transition(attempt["attempt_id"], "queued", child_run_dir=str(child))
        children[attempt["ordinal"]] = child
        return child

    def advance_scheduler(_seconds):
        sleeps.append(None)
        if len(sleeps) == 1:
            metadata = jobs.read_json(children[1] / "metadata.json")
            metadata.update(status="failed", worker_error="simulated worker failure")
            jobs.write_json(children[1] / "metadata.json", metadata)
        elif len(sleeps) == 2:
            metadata = jobs.read_json(children[2] / "metadata.json")
            metadata.update(status="completed", completed_at=jobs.utc_now())
            jobs.write_json(children[2] / "metadata.json", metadata)
            jobs.write_json(children[2] / "result.json", {"success": True})
            (children[2] / "artifacts/raw/bindcraft2/output/.campaign_state.json").write_text(json.dumps({
                "trajectories": 1, "accepted": 0, "rejections": {"terminated": {"screen": 1}},
            }))

    monkeypatch.setattr(distributed_campaign, "_create_attempt_child", create_child)
    monkeypatch.setattr(distributed_campaign.time, "sleep", advance_scheduler)
    _mock_finish(monkeypatch)

    distributed_campaign.run_bindcraft2_distributed_campaign(parent)

    state = CampaignLedger(parent / "artifacts/distributed_campaign/ledger.json").snapshot()
    assert len(children) == 2
    assert len(sleeps) == 2
    assert state["termination_reason"] == "budget_exhausted"
    second_attempt_id = f"{state['campaign_id']}:a00000002"
    assert state["attempts"][second_attempt_id]["status"] == "completed"
    assert state["candidates"] == {}
    assert CampaignLedger(parent / "artifacts/distributed_campaign/ledger.json").counters()[
        "native_rejected_trajectories"
    ] == 1


def test_generated_objective_allows_rfd4_without_claiming_native_acceptance():
    config = distributed_campaign.validate_campaign_config({
        'generator_id': 'rfd4_proteina', 'objective_type': 'generated',
        'accepted_target': 1000, 'max_attempts': 1100,
    })
    assert config['objective_type'] == 'generated'
    assert config['validation_stages'] == []
    assert config['accepted_target'] == 1000


def test_large_objectives_count_generated_outputs_and_unique_accepted_sequences():
    candidates = {
        str(i): {'native_acceptance': 'accepted', 'candidate': {'binder_sequence': f'SEQ{i // 2}'}}
        for i in range(2000)
    }
    candidates['missing'] = {'native_acceptance': 'accepted', 'candidate': {}}
    state = {'config': {'objective_type': 'generated', 'acceptance_level': 'native'}, 'candidates': candidates}
    assert distributed_campaign._state_accepted_count(state) == 2001
    state['config'].update(objective_type='accepted', unique_accepted=True)
    assert distributed_campaign._state_accepted_count(state) == 1000
    candidates['0']['native_acceptance'] = 'rejected'
    assert distributed_campaign._state_accepted_count(state) == 1000
    candidates['1']['native_acceptance'] = 'rejected'
    assert distributed_campaign._state_accepted_count(state) == 999


def test_generated_controller_stops_on_actual_outputs_without_native_acceptance(tmp_path, monkeypatch, design_source_factory):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    config = {
        'objective_type': 'generated', 'accepted_target': 2, 'max_attempts': 4,
        'acceptance_level': 'native', 'batch_size': 1, 'max_concurrent_batches': 2,
        'poll_interval_seconds': 0.25, 'max_retries': 0, 'seed': 21,
    }
    root, parent = _parent_run(tmp_path, monkeypatch, config)
    made = _controlled_workers(monkeypatch, root, design_source_factory, accepted=set())
    _mock_finish(monkeypatch)
    distributed_campaign.run_bindcraft2_distributed_campaign(parent)
    result = jobs.read_json(parent / 'result.json')
    assert result['metrics']['termination_reason'] == 'target_reached'
    assert result['metrics']['objective_type'] == 'generated'
    assert result['metrics']['objective_count'] == 2
    assert result['metrics']['generated_candidates'] == 2
    assert result['metrics']['native_accepted_candidates'] == 0
    assert len(made) == 2


def test_existing_ledger_without_objective_fields_remains_resumable(tmp_path):
    path = tmp_path / 'campaign.json'
    config = {'accepted_target': 1, 'max_attempts': 1}
    CampaignLedger.create(path, config, campaign_id='legacy')
    state = json.loads(path.read_text())
    state['config'].pop('objective_type')
    state['config'].pop('unique_accepted')
    path.write_text(json.dumps(state))
    assert CampaignLedger.create(path, config).snapshot()['campaign_id'] == 'legacy'


def test_resumed_legacy_campaign_normalizes_config_before_dispatch(tmp_path, monkeypatch):
    """Old parent inputs and ledgers acquire current defaults before scheduling."""
    config = {
        'accepted_target': 1,
        'max_attempts': 3,
        'batch_size': 2,
        'max_concurrent_batches': 2,
        'poll_interval_seconds': 0.25,
        'max_retries': 0,
    }
    _root, parent = _parent_run(tmp_path, monkeypatch, config)
    ledger_path = parent / 'artifacts/distributed_campaign/ledger.json'
    legacy_config = {**config, 'enabled': True, 'generator_id': 'bindcraft2'}
    campaign = CampaignLedger.create(ledger_path, legacy_config, campaign_id='legacy-resume')
    legacy_state = campaign.snapshot()
    legacy_state['config'].pop('max_concurrent_generation_jobs')
    legacy_state['config'].pop('objective_type')
    legacy_state['config'].pop('unique_accepted')
    ledger_path.write_text(json.dumps(legacy_state))

    dispatched = []
    monkeypatch.setattr(
        distributed_campaign,
        '_dispatch_attempts',
        lambda _parent, _ledger, normalized, **kwargs: dispatched.append((normalized, kwargs)),
    )

    class StopAfterDispatch(Exception):
        pass

    monkeypatch.setattr(
        distributed_campaign.time,
        'sleep',
        lambda _seconds: (_ for _ in ()).throw(StopAfterDispatch()),
    )
    with pytest.raises(StopAfterDispatch):
        distributed_campaign._run_distributed_campaign(parent, 'bindcraft2')

    assert len(dispatched) == 1
    normalized, _kwargs = dispatched[0]
    assert normalized['max_concurrent_generation_jobs'] == 4
    assert normalized['objective_type'] == 'accepted'
    assert normalized['unique_accepted'] is False
    assert distributed_campaign.validate_campaign_config(
        CampaignLedger(ledger_path).snapshot()['config']
    ) == normalized


def test_validation_boolean_policy_records_each_rule_and_missing_data():
    source = {'candidate_id': 'source', 'stage': 'generation.backbone_sequence',
              'source_tool': 'rfd4_proteina', 'binder_sequence': 'ACDE', 'metrics': {}}
    prediction = {'candidate_id': 'prediction', 'stage': 'complex_refolding',
                  'source_tool': 'optimized_prediction_boltz2', 'binder_sequence': 'ACDE',
                  'parents': ['source'], 'metrics': {'iptm': 0.4, 'ptm': 0.6},
                  'raw_metadata': {'source_candidate': source, 'measurement_status': {'status': 'measured'},
                                   'prediction_protocol': {'engine': 'boltz2'},
                                   'identity_check': {'sequence_chain_identity': 'pass'}}}
    rules = [{'field': 'iptm', 'operator': '>=', 'value': .7},
             {'field': 'ptm', 'operator': '>=', 'value': .5}]
    any_result = distributed_campaign._stage_selectors(source, prediction, {'match': 'any', 'rules': rules})
    all_result = distributed_campaign._stage_selectors(source, prediction, {'match': 'all', 'rules': rules})
    assert any_result['status'] == 'pass'
    assert all_result['status'] == 'reject'
    assert [item['status'] for item in any_result['condition_decisions']] == ['reject', 'pass']
    required = {'match': 'any', 'required_metrics': ['plddt_mean'], 'rules': rules}
    assert distributed_campaign._stage_selectors(source, prediction, required)['status'] == 'reject'
    assert distributed_campaign._stage_selectors(source, prediction, {**required, 'missing_data': 'unresolved'})['status'] == 'missing'
    assert distributed_campaign._stage_selectors(source, prediction, {**required, 'missing_data': 'ignore'})['status'] == 'pass'
    no_data = {'match': 'any', 'missing_data': 'ignore', 'rules': [{'field': 'plddt_mean', 'operator': '>=', 'value': 70}]}
    assert distributed_campaign._stage_selectors(source, prediction, no_data)['status'] == 'missing'


def test_bindcraft_peer_submission_retains_full_native_settings_and_one_attempt(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import candidate_batches
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    parent = tmp_path / 'parent'; child = tmp_path / 'child'
    raw = child / 'artifacts/raw/bindcraft2'; raw.mkdir(parents=True)
    (raw / 'target.cif').write_text('data_target\n')
    settings = {'targets': [{'target_path': 'target.cif', 'chains': 'A'}],
                'max_trajectories': 1, 'number_of_final_designs': 1, 'campaign_seed': 19,
                'trajectory_only': False, 'filters': {'i_pTM': {'threshold': .7, 'higher': True}},
                'screen_steps': 50}
    jobs.write_json(raw / 'settings.json', settings)
    jobs.write_json(parent / 'metadata.json', {'studio_owner_id': 'owner', 'studio_project_id': 'project'})
    jobs.write_json(child / 'input.json', {'params': {'target_file': 'target.cif', 'target_chains': ['A'],
                    'target_type': 'structure', 'binder_length': '60', 'campaign_seed': 19,
                    'max_trajectories': 1, 'number_of_final_designs': 1, 'workers_per_gpu': 1}})
    jobs.write_json(child / 'metadata.json', {'status': 'queued', 'resource_request': {'gpu': True, 'cpu_cores': 4}})
    config = {'generator_id': 'bindcraft2', 'workstation_selection': 'eligible',
              'accepted_target': 1, 'max_attempts': 2}
    ledger = CampaignLedger.create(parent / 'ledger.json', config, campaign_id='peer-bindcraft')
    completed = ledger.authorize(1, batch_id='previous')[0]
    ledger.transition(completed['attempt_id'], 'queued', cluster_submission={'placement_node_id': 'old-node'})
    ledger.record_receipt(completed['attempt_id'], 'previous-result', candidates=[],
                          native_state={'trajectories': 1, 'accepted': 0, 'rejections': {}})
    attempt = ledger.authorize(1, batch_id='batch')[0]
    captured = {}
    def submit(request, manifest, submission_key, placement, **kwargs):
        captured.update(request=request, placement=placement, submission_key=submission_key)
        kwargs['persist']()
        return {'cluster_job_id': 'peer-job', 'node_id': 'peer-node'}
    monkeypatch.setattr(candidate_batches, 'submit', submit)
    peer = distributed_campaign._submit_rfd4_peer_attempt(parent, ledger, attempt, child)
    parameters = captured['request']['parameters']
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
    allowed = workflow_automation._parameters_for(workflow_automation.WORKFLOWS['design.bindcraft2'])
    assert set(parameters) <= set(allowed)
    assert captured['request']['workflow'] == 'design.bindcraft2'
    assert captured['placement']['placement'] == 'eligible'
    assert 'old-node' not in captured['placement']['exclude_node_ids']
    assert parameters['target_pdb'] == str(raw / 'target.cif')
    assert 'trajectory_only' not in parameters  # Public full-workflow API defaults to the complete native pipeline.
    assert parameters['max_trajectories'] == 1
    assert parameters['engine_overrides']['filters'] == settings['filters']
    assert parameters['engine_overrides']['screen_steps'] == 50
    assert 'targets' not in parameters['engine_overrides']
    assert peer['cluster_job_id'] == 'peer-job'
    assert ledger.snapshot()['attempts'][attempt['attempt_id']]['cluster_submission']['cluster_job_id'] == 'peer-job'


def test_bindcraft_peer_collection_requires_authoritative_native_state(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import candidate_batches
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    parent = tmp_path / 'parent'; child = tmp_path / 'child'; child.mkdir()
    jobs.write_json(parent / 'metadata.json', {'studio_owner_id': 'owner', 'studio_project_id': 'project'})
    jobs.write_json(child / 'metadata.json', {'status': 'queued'})
    ledger = CampaignLedger.create(parent / 'ledger.json', {'accepted_target': 1, 'max_attempts': 1}, campaign_id='peer')
    attempt = ledger.authorize(1, batch_id='batch')[0]
    ledger.transition(attempt['attempt_id'], 'queued', child_run_dir=str(child))
    peer = {'manifest': str(tmp_path / 'transfer.json'), 'cluster_job_id': 'remote-job'}
    collected = Path(peer['manifest']).with_suffix('.collected')
    normalized = collected / 'artifacts/normalized_candidates/collected.jsonl'; normalized.parent.mkdir(parents=True)
    normalized.write_text('')
    jobs.write_json(collected / 'metadata.json', {'status': 'completed', 'resource_allocation': {'gpu_id': 0}})
    jobs.write_json(collected / 'result.json', {'success': True})
    native = collected / 'artifacts/raw/bindcraft2/output/.campaign_state.json'
    jobs.write_json(native, {'trajectories': 1, 'accepted': 0, 'rejections': {'screen': 1}})
    monkeypatch.setattr(candidate_batches, 'actor', lambda *_a, **_k: {'id': 'owner'})
    monkeypatch.setattr(candidate_batches, '_artifact', lambda *_a, **_k: None)
    assert not distributed_campaign._materialize_bindcraft2_peer_output(parent, ledger, attempt, child, peer, {'candidates_jsonl': str(normalized)}, {'node_id': 'remote'})
    monkeypatch.setattr(candidate_batches, '_artifact', lambda *_a, **_k: native)
    assert distributed_campaign._materialize_bindcraft2_peer_output(parent, ledger, attempt, child, peer, {'candidates_jsonl': str(normalized)}, {'node_id': 'remote'})
    from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import read_bindcraft2_campaign_state
    assert read_bindcraft2_campaign_state(child)['trajectories'] == 1
    assert read_bindcraft2_campaign_state(child)['accepted'] == 0
    assert jobs.read_json(child / 'metadata.json')['execution_node_id'] == 'remote'


def test_thousand_design_dispatch_stays_bounded_and_never_reuses_attempt_identity(tmp_path, monkeypatch):
    class SimulatedLedger:
        def __init__(self):
            self.state = {'next_ordinal': 1, 'attempts': {}}
        def assert_production_jobs(self):
            pass
        def snapshot(self):
            return self.state
        def authorize(self, requested, *, batch_id, **_kwargs):
            assert requested == 1
            ordinal = self.state['next_ordinal']
            attempt = {'attempt_id': f'a{ordinal:04d}', 'ordinal': ordinal, 'seed': 10000 + ordinal, 'batch_id': batch_id}
            self.state['attempts'][attempt['attempt_id']] = attempt
            self.state['next_ordinal'] += 1
            return [attempt]
    ledger = SimulatedLedger()
    launched = []
    monkeypatch.setattr(distributed_campaign, '_create_attempt_child',
                        lambda _parent, _ledger, attempt, _config: launched.append(attempt))
    config = distributed_campaign.validate_campaign_config({
        'accepted_target': 1000, 'max_attempts': 1100,
        'batch_size': 2, 'max_concurrent_batches': 2,
    })
    for completed in range(0, 1000, 4):
        before = len(launched)
        distributed_campaign._dispatch_attempts(tmp_path, ledger, config, active_count=0, accepted_count=completed)
        assert len(launched) - before <= 4
    assert len(launched) == 1000
    assert len({a['attempt_id'] for a in launched}) == 1000
    assert len({a['seed'] for a in launched}) == 1000
    distributed_campaign._dispatch_attempts(tmp_path, ledger, config, active_count=0, accepted_count=1000)
    assert len(launched) == 1000


def test_validation_job_metadata_tracks_actual_campaign_generator(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    monkeypatch.setattr(jobs, 'runs_root', lambda: tmp_path)
    parent = tmp_path / 'design' / 'campaign-parent'
    for generator, expected_role in [('rfd4_proteina', 'distributed_rfd4_validation'),
                                     ('bindcraft2', 'distributed_bindcraft2_validation')]:
        with jobs.submission_tags(distributed_campaign_id='campaign',
                                  distributed_generator_id=generator,
                                  distributed_validation_stage=0,
                                  parent_campaign_run_dir=str(parent)):
            child = jobs.create_job('refolding-validation', 'complex_refolding',
                                    'optimized_prediction', {}, {})
        metadata = jobs.read_json(child.run_dir / 'metadata.json')
        assert metadata['benchmark_engine'] == generator
        assert metadata['parent_role'] == expected_role
