"""Deterministic accounting/recovery fixtures. No native or validation inference."""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import json

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows import distributed_campaign as dc
from mn_molecular_design_studio.engine.workflow_engine.workflows import execution_job_ledger as job_contract
from mn_molecular_design_studio.engine.workflow_engine.workflows.execution_job_ledger import LedgerConflict, digest
from mn_molecular_design_studio.engine.workflow_engine.workflows.execution_planning import plan_campaign_execution


@pytest.fixture(autouse=True)
def fixed_ledger_clock(monkeypatch):
    monkeypatch.setattr(dc, "utc_now", lambda: "2026-10-28T10:00:00+00:00")
    monkeypatch.setattr(job_contract, "utc_now", lambda: "2026-10-28T10:00:00+00:00")


def ledger(tmp_path, budget=6, **config):
    return dc.CampaignLedger.create(tmp_path / "ledger.json", {"max_attempts": budget, "seed": 17, **config}, campaign_id="fixture")


def reserve(campaign, count=3, key="reservation"):
    return campaign.reserve_execution_job(key, [{"input_id": f"input-{i % 2}", "sha256": digest(i % 2)} for i in range(count)],
        native_configuration={"engine": "fixture", "parameters": {"unchanged": True}},
        resource_requirements={"gpu_count": 1, "min_vram_gb": 16})


def execution(job):
    return job["executions"][-1]


def result(campaign, job, index=0, *, receipt=None, candidates=None, native_state=None, outcome="completed"):
    return campaign.record_attempt_result(job["execution_job_id"], job["current_execution_id"], receipt or f"receipt-{index}",
        job["manifest"][index]["attempt_id"], manifest_digest=job["manifest_digest"], outcome=outcome,
        native_state={"trajectories": 1, "accepted": 1} if native_state is None else native_state,
        native_output_ids=[f"native-{index}"], candidates=candidates or [], evidence={"verified_artifact_sha256": digest(index)})


def finish(campaign, job, **kwargs):
    return campaign.finish_execution(job["execution_job_id"], job["current_execution_id"],
        evidence={"worker_terminal_record": "fixture"}, **kwargs)


def test_atomic_manifest_and_idempotent_authorization_across_controllers(tmp_path):
    campaign = ledger(tmp_path, budget=5)
    with ThreadPoolExecutor(max_workers=8) as pool:
        jobs = list(pool.map(lambda _: reserve(dc.CampaignLedger(campaign.path)), range(8)))
    assert all(j == jobs[0] for j in jobs)
    snapshot = campaign.snapshot()
    assert len(snapshot["attempts"]) == 3
    assert len(snapshot["execution_jobs"]) == 1
    assert snapshot["next_ordinal"] == 4
    before = campaign.path.read_bytes()
    with pytest.raises(ValueError, match="entire manifest"):
        reserve(campaign, 3, "over-budget")
    with pytest.raises(LedgerConflict, match="request key"):
        reserve(campaign, 2)
    assert campaign.path.read_bytes() == before
    assert [m["job_position"] for m in jobs[0]["manifest"]] == [0, 1, 2]
    assert [m["seed"] for m in jobs[0]["manifest"]] == [dc.CampaignLedger._seed("fixture", 17, i) for i in (1, 2, 3)]


def test_authorization_write_failure_has_no_partial_members(tmp_path, monkeypatch):
    campaign = ledger(tmp_path)
    before = campaign.path.read_bytes()
    with monkeypatch.context() as patch:
        patch.setattr(dc.os, "replace", lambda *_: (_ for _ in ()).throw(OSError("crash before atomic replacement")))
        with pytest.raises(OSError, match="crash"):
            reserve(campaign)
    assert campaign.path.read_bytes() == before
    assert not list(tmp_path.glob("*.tmp"))
    recovered = dc.CampaignLedger.create(campaign.path, {"max_attempts": 6, "seed": 17})
    assert reserve(recovered)["manifest"][0]["ordinal"] == 1


def test_crash_before_submission_and_duplicate_submission_keep_one_identity(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    assert execution(job)["submission_state"] == "reserved"
    recovered = dc.CampaignLedger(campaign.path)
    first = recovered.prepare_execution_submission(job["execution_job_id"])
    before = campaign.path.read_bytes()
    assert recovered.prepare_execution_submission(job["execution_job_id"]) == first
    assert campaign.path.read_bytes() == before
    args = dict(submission_key=first["submission_key"], status="queued", scheduler_identifiers={"scheduler_job_id": "scheduler-42"})
    confirmed = recovered.observe_execution(job["execution_job_id"], first["execution_id"], **args)
    assert recovered.observe_execution(job["execution_job_id"], first["execution_id"], **args) == confirmed
    with pytest.raises(LedgerConflict, match="scheduler_identifiers"):
        recovered.observe_execution(job["execution_job_id"], first["execution_id"], **{**args, "scheduler_identifiers": {"scheduler_job_id": "duplicate-99"}})
    assert len(campaign.snapshot()["attempts"]) == 3
    assert len(campaign.snapshot()["execution_jobs"][job["execution_job_id"]]["executions"]) == 1


def test_partial_completion_missing_receipts_and_late_reconnection(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    intent = campaign.prepare_execution_submission(job["execution_job_id"])
    args = dict(submission_key=intent["submission_key"], status="running", scheduler_identifiers={"peer_job_id": "remote-1"},
                worker_identity={"node_id": "peer-A", "gpu_uuid": "GPU-A"}, resource_allocation={"lease_id": "lease-A", "gpu_count": 1})
    campaign.observe_execution(job["execution_job_id"], intent["execution_id"], **args)
    campaign.observe_execution(job["execution_job_id"], intent["execution_id"], **args, contact_state="disconnected")
    recovered = dc.CampaignLedger(campaign.path)
    assert recovered.counters()["native_attempts_started"] == 0  # worker start proves no per-attempt start
    recovered.observe_execution(job["execution_job_id"], intent["execution_id"], **args, contact_state="connected")
    recovered.record_native_start(job["execution_job_id"], job["current_execution_id"], job["manifest"][1]["attempt_id"],
                                 manifest_digest=job["manifest_digest"], evidence={"native_event": "attempt-2-started"})
    rows = [{"candidate_id": f"c{i}", "metrics": {"native_pass_filters": True}, "binder_sequence": "ACDE"} for i in range(2)]
    assert result(recovered, job, candidates=rows)
    assert finish(recovered, job, status="interrupted")
    counters = recovered.counters()
    assert counters["attempts_authorized"] == 3
    assert counters["native_attempts_started"] == 2
    assert counters["native_attempts_completed"] == 1
    assert counters["native_attempts_unresolved"] == 2
    assert counters["native_accepted_candidates"] == 2
    assert result(recovered, job, 1, native_state={"trajectories": 1, "accepted": 0})
    assert recovered.counters()["native_rejected_trajectories"] == 1
    assert recovered.counters()["native_attempts_unresolved"] == 1
    assert len(recovered.snapshot()["candidates"]) == 2
    with pytest.raises(ValueError, match="post-launch"):
        recovered.retry_execution_job(job["execution_job_id"], retry_key="unsafe-replay")


def test_worker_start_time_is_observed_measurement_not_immutable_identity(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign, count=1)
    intent = campaign.prepare_execution_submission(job["execution_job_id"])
    campaign.observe_execution(
        job["execution_job_id"], intent["execution_id"],
        submission_key=intent["submission_key"], status="running",
        worker_identity={"scheduler_worker_id": "worker-A", "execution_node_id": "node-A",
                         "worker_started_at": "2026-10-29T20:10:24+00:00"},
    )
    campaign.observe_execution(
        job["execution_job_id"], intent["execution_id"],
        submission_key=intent["submission_key"], status="running",
        worker_identity={"scheduler_worker_id": "worker-A", "execution_node_id": "node-A",
                         "worker_started_at": "2026-10-29T20:27:35+00:00"},
    )
    campaign.observe_worker_metadata(job["manifest"][0]["attempt_id"], {
        "scheduler_worker_id": "worker-A", "execution_node_id": "node-A",
        "worker_started_at": "2026-10-29T20:10:24+00:00",
    })
    campaign.observe_worker_metadata(job["manifest"][0]["attempt_id"], {
        "scheduler_worker_id": "worker-A", "execution_node_id": "node-A",
        "worker_started_at": "2026-10-29T20:27:35+00:00",
    })
    saved = execution(campaign.snapshot()["execution_jobs"][job["execution_job_id"]])
    assert saved["worker_identity"] == {
        "scheduler_worker_id": "worker-A", "execution_node_id": "node-A",
    }
    assert saved["worker_timing"]["latest_worker_started_at"] == "2026-10-29T20:27:35+00:00"
    assert [row["worker_started_at"] for row in saved["worker_timing"]["worker_started_at_observations"]] == [
        "2026-10-29T20:10:24+00:00", "2026-10-29T20:27:35+00:00",
    ]


def test_remote_worker_node_reassignment_is_recorded_without_failing_campaign(tmp_path):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="dispatch")[0]
    attempt_id = attempt["attempt_id"]
    details = {"child_run_dir": "/runs/child-1", "cluster_job_id": "cluster-1"}
    campaign.transition(attempt_id, "queued", **details, worker_node_id="node-a")
    campaign.transition(attempt_id, "running", **details, worker_node_id="node-b")

    snapshot = campaign.snapshot()
    attempt = snapshot["attempts"][attempt_id]
    saved = execution(snapshot["execution_jobs"][attempt["execution_job_id"]])
    assert attempt["worker_node_id"] == "node-b"
    assert saved["worker_identity"]["worker_node_id"] == "node-b"
    assert [(row["worker_node_id"], row["status"]) for row in saved["worker_node_id_observations"]] == [
        ("node-a", "queued"), ("node-b", "running"),
    ]


def test_repeated_attempt_transition_is_a_noop(tmp_path):
    campaign = ledger(tmp_path, budget=1)
    attempt = campaign.authorize(1, batch_id="dispatch")[0]
    details = {"child_run_dir": "/runs/child-1", "cluster_job_id": "cluster-1", "worker_node_id": "node-a"}

    campaign.transition(attempt["attempt_id"], "queued", **details)
    before = campaign.path.read_bytes()
    event_count = len(campaign.snapshot()["events"])

    campaign.transition(attempt["attempt_id"], "queued", **details)

    assert campaign.path.read_bytes() == before
    assert len(campaign.snapshot()["events"]) == event_count


def test_worker_metadata_snapshot_skips_reloading_unchanged_ledger(tmp_path, monkeypatch):
    campaign = ledger(tmp_path, budget=1)
    attempt = campaign.authorize(1, batch_id="dispatch")[0]
    metadata = {
        "scheduler_worker_id": "worker-A",
        "execution_node_id": "node-A",
        "worker_started_at": "2026-10-29T20:10:24+00:00",
        "scheduler_job_id": "scheduler-A",
        "cluster_job_id": "cluster-A",
        "remote_submission_key": "submission-A",
        "resource_allocation": {"gpu_device": "0", "gpu_ids": [0]},
        "resource_request": {"gpu": True},
    }
    campaign.observe_worker_metadata(attempt["attempt_id"], metadata)
    known_state = campaign.snapshot()

    def unexpected_read(_path):
        raise AssertionError("unchanged worker metadata should not reload the ledger")

    monkeypatch.setattr(dc, "_read", unexpected_read)
    campaign.observe_worker_metadata(attempt["attempt_id"], metadata, known_state=known_state)


def test_exit_zero_does_not_complete_attempts_or_infer_acceptance(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    campaign.prepare_execution_submission(job["execution_job_id"])
    finish(campaign, job, status="completed")
    assert campaign.counters()["native_attempts_completed"] == 0
    assert campaign.counters()["native_attempts_unresolved"] == 3
    assert campaign.counters()["native_rejected_trajectories"] == 0
    assert campaign.counters()["native_accepted_candidates"] == 0


def test_receipt_replays_are_strict_and_conflicts_are_atomic(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    rows = [{"candidate_id": "same-id", "metrics": {"native_pass_filters": True}}]
    assert result(campaign, job, candidates=rows)
    before = campaign.path.read_bytes()
    assert not result(campaign, job, candidates=rows)
    assert not result(campaign, job, candidates=rows, receipt="same-content-alias")
    assert campaign.path.read_bytes() == before
    for call in (lambda: result(campaign, job, candidates=rows, native_state={"accepted": 99}),
                 lambda: result(campaign, job, 1, receipt="receipt-0"),
                 lambda: result(campaign, job, 1, candidates=rows),
                 lambda: campaign.record_attempt_result(job["execution_job_id"], job["current_execution_id"], "bad-manifest",
                    job["manifest"][1]["attempt_id"], manifest_digest="bad", native_state={}, native_output_ids=[], candidates=[], evidence={"proof": True})):
        with pytest.raises(LedgerConflict):
            call()
        assert campaign.path.read_bytes() == before


def test_cross_job_receipt_membership_conflicts(tmp_path):
    campaign = ledger(tmp_path)
    a, b = reserve(campaign, key="a"), reserve(campaign, key="b")
    before = campaign.path.read_bytes()
    with pytest.raises(LedgerConflict, match="not a member"):
        campaign.record_attempt_result(a["execution_job_id"], a["current_execution_id"], "wrong-owner", b["manifest"][0]["attempt_id"],
            manifest_digest=a["manifest_digest"], native_state={}, native_output_ids=[], candidates=[], evidence={"proof": True})
    assert campaign.path.read_bytes() == before


@pytest.mark.parametrize("field", ["seed", "source_input", "job_position", "execution_job_id", "native_rng_mapping"])
def test_immutable_membership_rejects_edits_and_tampering(tmp_path, field):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="dispatch-group")[0]
    before = campaign.path.read_bytes()
    with pytest.raises(ValueError, match="immutable"):
        campaign.transition(attempt["attempt_id"], "queued", **{field: "changed"})
    assert campaign.path.read_bytes() == before
    state = campaign.snapshot()
    state["attempts"][attempt["attempt_id"]][field] = "changed"
    campaign.path.write_text(json.dumps(state))
    with pytest.raises(LedgerConflict, match="Immutable"):
        campaign.snapshot()


def test_configuration_digest_and_order_are_checked_on_restart(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    state = campaign.snapshot()
    state["execution_jobs"][job["execution_job_id"]]["manifest"].reverse()
    campaign.path.write_text(json.dumps(state))
    with pytest.raises(LedgerConflict, match="manifest digest"):
        dc.CampaignLedger(campaign.path).snapshot()


def test_prelaunch_retry_keeps_attempt_identity_budget_and_rng(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    original = deepcopy(job["manifest"])
    campaign.prepare_execution_submission(job["execution_job_id"])
    finish(campaign, job, status="failed", not_started=True)
    retry = campaign.retry_execution_job(job["execution_job_id"], retry_key="retry-1")
    before = campaign.path.read_bytes()
    assert campaign.retry_execution_job(job["execution_job_id"], retry_key="retry-1") == retry
    assert campaign.path.read_bytes() == before
    saved = campaign.snapshot()["execution_jobs"][job["execution_job_id"]]
    assert saved["manifest"] == original
    assert retry["submission_key"] != execution(job)["submission_key"]
    assert retry["execution_id"] != job["current_execution_id"]
    assert campaign.counters()["attempts_authorized"] == 3
    assert campaign.counters()["execution_retries"] == 1
    assert campaign.counters()["native_attempts_authorized"] == 3


@pytest.mark.parametrize("status", ["failed", "interrupted", "cancelled", "completed"])
def test_unresolved_postlaunch_work_never_automatically_retries(tmp_path, status):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    campaign.prepare_execution_submission(job["execution_job_id"])
    finish(campaign, job, status=status)
    with pytest.raises(ValueError, match="post-launch"):
        campaign.retry_execution_job(job["execution_job_id"], retry_key="unsafe")
    assert len(campaign.snapshot()["attempts"]) == 3


def test_queued_cancellation_requires_no_launch_evidence(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    intent = campaign.prepare_execution_submission(job["execution_job_id"])
    campaign.observe_execution(job["execution_job_id"], intent["execution_id"], submission_key=intent["submission_key"], status="queued")
    assert campaign.cancel_execution_job(job["execution_job_id"])["cancel_requested"]
    assert campaign.counters()["native_attempts_cancelled"] == 0
    finish(campaign, job, status="cancelled", not_started=True)
    assert campaign.counters()["native_attempts_cancelled"] == 3
    assert campaign.counters()["native_attempts_started"] == 0
    with pytest.raises(LedgerConflict, match="never started"):
        result(campaign, job)
    assert campaign.counters()["attempts_authorized"] == 3


def test_unsubmitted_cancellation_needs_no_scheduler_action(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    first = campaign.cancel_execution_job(job["execution_job_id"])
    before = campaign.path.read_bytes()
    assert campaign.cancel_execution_job(job["execution_job_id"]) == first
    assert campaign.path.read_bytes() == before
    assert campaign.counters()["native_attempts_cancelled"] == 3
    assert campaign.counters()["native_attempts_started"] == 0


def test_failed_native_attempt_is_distinct_from_unresolved_worker_failure(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    result(campaign, job, outcome="failed", native_state={"native_error": "deterministic-fixture"})
    finish(campaign, job, status="failed")
    assert campaign.counters()["native_attempts_failed"] == 1
    assert campaign.counters()["native_attempts_unresolved"] == 2
    assert campaign.counters()["native_rejected_trajectories"] == 0


def test_job_lease_counted_once_with_many_attempts_and_validation_independent(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    assert campaign.record_job_compute(job["execution_job_id"], job["current_execution_id"], "gpu-lease", 20)
    assert not campaign.record_job_compute(job["execution_job_id"], job["current_execution_id"], "gpu-lease", 20)
    with pytest.raises(LedgerConflict, match="duration"):
        campaign.record_job_compute(job["execution_job_id"], job["current_execution_id"], "gpu-lease", 21)
    member = job["manifest"][0]["attempt_id"]
    with pytest.raises(ValueError, match="recorded once"):
        campaign.record_compute(member, "per-member-duplication", 20)
    campaign.record_compute(member, "validation-1", 7, scope="validation")
    assert campaign.gpu_compute_seconds() == 27
    assert all(not a.get("compute_runs") for a in list(campaign.snapshot()["attempts"].values())[1:])


def test_unavailable_native_acceptance_late_validation_and_surplus_selection(tmp_path):
    stage = {"workflow": "refolding.optimized_prediction", "required_metrics": ["iptm"], "parameters": {}}
    campaign = ledger(tmp_path, budget=3, generator_id="rfd4_proteina", acceptance_level="validated", validation_stages=[stage])
    job = reserve(campaign)
    rows = [{"candidate_id": f"c{i}", "binder_sequence": "ACDE", "source_structure": f"native-{i}.pdb", "metrics": {"native_pass_filters": True}} for i in range(2)]
    result(campaign, job, candidates=rows, native_state={"authoritative": False, "output_candidates": 2})
    result(campaign, job, 1, native_state={"authoritative": False})
    result(campaign, job, 2, native_state={"authoritative": False})
    campaign.mark_generation_budget_exhausted()
    assert campaign.snapshot()["status"] == "running"
    assert campaign.counters()["native_accepted_candidates"] == 0
    assert campaign.counters()["native_acceptance_unavailable_candidates"] == 2
    campaign.terminate("budget_exhausted")
    recovered = dc.CampaignLedger(campaign.path)
    recovered.set_candidate_acceptance("c0", validated="accepted")
    recovered.set_candidate_acceptance("c1", validated="accepted")
    recovered.update_candidate("c0", selection_state="selected")
    recovered.update_candidate("c1", selection_state="returned_unselected", selection_reason="target already filled")
    assert recovered.terminate("target_reached")["termination_reason"] == "target_reached"
    state = recovered.snapshot()
    assert state["candidates"]["c1"]["candidate"] == rows[1]
    assert recovered.counters()["attempts_authorized"] == 3
    assert recovered.counters()["validated_accepted_candidates"] == 2


def test_old_ledgers_remain_schema_one_and_batch_identity_stays_dispatch_grouping(tmp_path):
    campaign = dc.CampaignLedger.create(tmp_path / "old.json", {"max_attempts": 3}, campaign_id="old", schema_version=1)
    attempts = campaign.authorize(2, batch_id="b-old")
    campaign.transition(attempts[0]["attempt_id"], "failed")
    campaign.retry(attempts[0]["attempt_id"])
    original = campaign.path.read_bytes()
    recovered = dc.CampaignLedger.create(campaign.path, {"max_attempts": 3})
    assert recovered.path.read_bytes() == original
    assert recovered.snapshot()["schema_version"] == 1
    assert "execution_jobs" not in recovered.snapshot()
    assert all(a["batch_id"] == "b-old" and "execution_job_id" not in a for a in recovered.snapshot()["attempts"].values())
    assert recovered.execution_job_summary()["execution_jobs"] is None
    with pytest.raises(ValueError, match="schema v2"):
        reserve(recovered)


def test_production_authorization_stays_one_attempt_per_job_and_rejects_accounting_jobs(tmp_path, monkeypatch):
    campaign = ledger(tmp_path)
    attempts = campaign.authorize(3, batch_id="same-dispatch-group")
    assert len({a["execution_job_id"] for a in attempts}) == 3
    assert all(len(j["manifest"]) == 1 for j in campaign.snapshot()["execution_jobs"].values())
    campaign.assert_production_jobs()
    job = reserve(campaign)
    monkeypatch.setattr(dc, "_create_rfd4_attempt_child", lambda *a: pytest.fail("native dispatch"))
    with pytest.raises(ValueError, match="Accounting-only"):
        dc._create_attempt_child(tmp_path, campaign, campaign.snapshot()["attempts"][job["manifest"][0]["attempt_id"]], {"generator_id": "rfd4_proteina"})
    with pytest.raises(ValueError, match="Accounting-only"):
        campaign.transition(job["manifest"][0]["attempt_id"], "queued")


def test_production_postlaunch_retries_preserve_rfd4_guard(tmp_path):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="b1")[0]
    child = tmp_path / "child"
    child.mkdir()
    (child / "command.json").write_text(json.dumps({"mode": "docker"}))
    campaign.transition(attempt["attempt_id"], "running", child_run_dir=str(child))
    campaign.transition(attempt["attempt_id"], "failed")
    assert not campaign.retry_is_safe(campaign.snapshot()["attempts"][attempt["attempt_id"]])
    with pytest.raises(ValueError, match="post-launch"):
        campaign.retry(attempt["attempt_id"])
    assert campaign.counters()["attempts_authorized"] == 1


def test_planner_exact_multi_input_manifest_is_read_only_and_unsupported(monkeypatch, tmp_path):
    monkeypatch.setattr(dc.CampaignLedger, "create", lambda *a, **kw: pytest.fail("ledger created by planner"))
    monkeypatch.setattr(dc, "_dispatch_attempts", lambda *a, **kw: pytest.fail("planner dispatched"))
    config = {"max_attempts": 13, "seed": 17}
    original = deepcopy(config)
    plan = plan_campaign_execution(config, native_attempts_per_job=5, input_count=2)
    preview = plan["execution_job_manifest_preview"]
    assert plan == plan_campaign_execution(config, native_attempts_per_job=5, input_count=2)
    assert config == original
    assert not plan["submission_authorized"]
    assert not plan["campaign_configuration_supported"]
    assert plan["supported_execution_plan"] is None
    assert not preview["authorization"] and not preview["truncated"]
    assert [j["proposed_attempt_count"] for j in preview["jobs"]] == [5, 2, 5, 1]
    members = [m for j in preview["jobs"] for m in j["manifest"]]
    assert [m["ordinal"] for m in members] == list(range(1, 14))
    assert len({m["attempt_id"] for m in members}) == 13
    assert [sum(m["source_input"]["input_index"] == i for m in members) for i in range(2)] == [7, 6]
    assert all(digest(j["manifest"]) == j["manifest_digest"] for j in preview["jobs"])
    assert list(tmp_path.iterdir()) == []


def test_huge_plans_have_bounded_explicit_previews_and_exact_totals():
    preview = plan_campaign_execution({"max_attempts": 10**12}, native_attempts_per_job=10**9)["execution_job_manifest_preview"]
    assert preview["total_jobs"] == 1000
    assert preview["total_attempts"] == 10**12
    assert preview["shown_attempts"] == 100
    assert preview["truncated"]
    assert preview["jobs"][0]["manifest_digest"] is None
    assert preview["jobs"][0]["omitted_members"] == 10**9 - 100


def test_no_launch_claim_conflicts_with_native_start_evidence(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    campaign.record_native_start(job["execution_job_id"], job["current_execution_id"], job["manifest"][0]["attempt_id"],
                                manifest_digest=job["manifest_digest"], evidence={"native_event": "start-1"})
    before = campaign.path.read_bytes()
    with pytest.raises(LedgerConflict, match="No-launch"):
        finish(campaign, job, status="failed", not_started=True)
    assert campaign.path.read_bytes() == before
    result(campaign, job)
    attempt = campaign.snapshot()["attempts"][job["manifest"][0]["attempt_id"]]
    assert attempt["native_start_evidence"] == {"native_event": "start-1"}
    assert attempt["result_evidence"] == {"verified_artifact_sha256": digest(0)}


def test_scheduler_lease_cannot_move_between_jobs(tmp_path):
    campaign = ledger(tmp_path)
    a, b = reserve(campaign, key="a"), reserve(campaign, key="b")
    campaign.record_job_compute(a["execution_job_id"], a["current_execution_id"], "lease-1", 9)
    with pytest.raises(LedgerConflict, match="multiple executions"):
        campaign.record_job_compute(b["execution_job_id"], b["current_execution_id"], "lease-1", 9)
    assert campaign.gpu_compute_seconds() == 9


def test_malformed_manifest_and_changed_configuration_refuse_recovery(tmp_path):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    good = campaign.snapshot()
    for mutate in (lambda s: s["execution_jobs"][job["execution_job_id"]].pop("manifest"),
                   lambda s: s["execution_jobs"][job["execution_job_id"]]["configuration"]["native"].update(engine="other")):
        state = deepcopy(good)
        mutate(state)
        campaign.path.write_text(json.dumps(state))
        with pytest.raises(LedgerConflict):
            dc.CampaignLedger(campaign.path).snapshot()


def test_single_large_job_preview_reports_omitted_members():
    preview = plan_campaign_execution({"max_attempts": 1000}, native_attempts_per_job=1000)["execution_job_manifest_preview"]
    assert preview["truncated"]
    assert preview["omitted_jobs"] == 0
    assert preview["omitted_attempts"] == 900
    assert preview["jobs"][0]["manifest_complete"] is False
    assert preview["jobs"][0]["native_configuration_resolved"] is False


def test_prelaunch_production_recovery_records_execution_retry_before_resuming(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    campaign = ledger(tmp_path, budget=1)
    child = tmp_path / "child"
    child.mkdir()
    (child / "command.json").write_text(json.dumps({"mode": "local_worker"}))
    jobs.write_json(child / "metadata.json", {"status": "failed", "worker_error": "setup failed", "completed_at": "2026-10-28T10:00:00+00:00"})
    jobs.write_json(child / "result.json", {"success": False})
    attempt = campaign.authorize(1, batch_id="b1")[0]
    campaign.transition(attempt["attempt_id"], "running", child_run_dir=str(child))
    monkeypatch.setattr(dc, "_attempt_child_paths", lambda _, state=None: {attempt["attempt_id"]: child})
    resumed = []
    def resume(path, **kwargs):
        state = campaign.snapshot()
        job = state["execution_jobs"][attempt["execution_job_id"]]
        assert len(job["executions"]) == 2
        assert execution(job)["submission_state"] == "queued"
        assert state["attempts"][attempt["attempt_id"]]["retry_count"] == 1
        jobs.write_json(path / "metadata.json", {"status": "queued"})
        resumed.append(path)
    monkeypatch.setattr(jobs, "resume_job", resume)
    _, active, _ = dc._reconcile_attempts(campaign, tmp_path, campaign.snapshot()["config"])
    assert resumed == [child]
    assert len(active) == 1
    assert campaign.counters()["attempts_authorized"] == 1
    assert campaign.counters()["execution_retries"] == 1
    assert campaign.counters()["native_attempts_started"] == 0


def test_preexisting_schema_one_controller_does_not_retry_missing_or_native_launched_work(tmp_path):
    campaign = dc.CampaignLedger.create(tmp_path / "old.json", {"max_attempts": 2}, schema_version=1)
    a, b = campaign.authorize(2, batch_id="legacy")
    child = tmp_path / "child"
    child.mkdir()
    (child / "command.json").write_text(json.dumps({"mode": "docker"}))
    campaign.transition(a["attempt_id"], "running", child_run_dir=str(child))
    campaign.transition(b["attempt_id"], "running", child_run_dir=str(tmp_path / "missing"))
    assert not campaign.retry_is_safe(campaign.snapshot()["attempts"][a["attempt_id"]])
    assert not campaign.retry_is_safe(campaign.snapshot()["attempts"][b["attempt_id"]])
    assert campaign.snapshot()["schema_version"] == 1


@pytest.mark.parametrize("changes", [
    {"request_key": True}, {"request_key": " "}, {"source_inputs": []}, {"source_inputs": [{}]},
    {"native_configuration": []}, {"resource_requirements": []}, {"resource_requirements": "gpu"},
])
def test_invalid_reservations_never_consume_authorization(tmp_path, changes):
    campaign = ledger(tmp_path)
    request = {"request_key": "test", "source_inputs": [{"input_id": "target-A"}]}
    before = campaign.path.read_bytes()
    with pytest.raises(ValueError):
        campaign.reserve_execution_job(**{**request, **changes})
    assert campaign.path.read_bytes() == before
    assert campaign.counters()["attempts_authorized"] == 0


def test_single_adapter_receipts_keep_native_identifiers_and_full_lineage(tmp_path):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="dispatch-group")[0]
    candidate = {"candidate_id": "global-normalized-id", "binder_sequence": "ACDE", "complex_pdb": "native-structure.cif",
                 "raw_metadata": {"native_candidate_id": "native-local-id"}, "metrics": {}}
    campaign.record_receipt(attempt["attempt_id"], "receipt-1", native_state={"authoritative": False}, candidates=[candidate])
    record = campaign.snapshot()["candidates"]["global-normalized-id"]
    assert record["native_output_ids"] == ["native-local-id"]
    assert record["candidate"] == candidate
    assert record["native_acceptance"] == "unavailable"
    before = campaign.path.read_bytes()
    with pytest.raises(ValueError, match="provenance"):
        campaign.update_candidate("global-normalized-id", candidate={"binder_sequence": "DIFFERENT"})
    assert campaign.path.read_bytes() == before


def test_worker_failure_cannot_erase_verified_native_start_or_enable_retry(tmp_path):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="dispatch-group")[0]
    job = campaign.snapshot()["execution_jobs"][attempt["execution_job_id"]]
    child = tmp_path / "child"
    child.mkdir()
    (child / "command.json").write_text(json.dumps({"mode": "local_worker"}))
    campaign.transition(attempt["attempt_id"], "running", child_run_dir=str(child))
    campaign.record_native_start(job["execution_job_id"], job["current_execution_id"], attempt["attempt_id"],
                                manifest_digest=job["manifest_digest"], evidence={"native_event": "verified-start"})
    campaign.transition(attempt["attempt_id"], "failed")
    saved = campaign.snapshot()["attempts"][attempt["attempt_id"]]
    assert saved["native_started"] is True
    assert saved["native_status"] == "unresolved"
    assert not campaign.retry_is_safe(saved)
    with pytest.raises(ValueError, match="Started scientific"):
        campaign.retry(attempt["attempt_id"])
    assert campaign.counters()["native_attempts_started"] == 1


@pytest.mark.parametrize("evidence_kind", ["native_start", "native_result"])
def test_native_evidence_before_scheduler_ack_prevents_unstarted_cancellation(tmp_path, evidence_kind):
    campaign = ledger(tmp_path)
    job = reserve(campaign)
    member = job["manifest"][0]["attempt_id"]
    if evidence_kind == "native_start":
        campaign.record_native_start(job["execution_job_id"], job["current_execution_id"], member,
                                    manifest_digest=job["manifest_digest"], evidence={"native_event": "verified-start"})
    else:
        result(campaign, job)
    assert campaign.cancel_execution_job(job["execution_job_id"])["cancel_requested"]
    assert campaign.snapshot()["attempts"][member]["native_started"] is True
    assert campaign.counters()["native_attempts_cancelled"] == 0
    finish(campaign, job, status="interrupted")
    assert campaign.snapshot()["attempts"][member]["native_started"] is True
