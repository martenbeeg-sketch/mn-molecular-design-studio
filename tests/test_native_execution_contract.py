"""Characterize legacy compatibility and simulate grouping; no native dispatch."""
import json

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import (
    CampaignLedger, validate_campaign_config,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import list_generator_capabilities


def test_clear_dispatch_alias_keeps_legacy_ledger_recovery_compatible(tmp_path):
    legacy = {"batch_size": 5, "max_concurrent_batches": 2, "max_attempts": 10}
    ledger = CampaignLedger.create(tmp_path / "ledger.json", legacy)
    before = ledger.path.read_bytes()
    alias = {"dispatch_batch_size": 5, "native_attempts_per_job": 1,
             "max_concurrent_batches": 2, "max_attempts": 10}
    assert validate_campaign_config(alias) == validate_campaign_config(legacy)
    assert CampaignLedger.create(ledger.path, alias).snapshot() == ledger.snapshot()
    assert ledger.path.read_bytes() == before
    assert alias["native_attempts_per_job"] == 1


@pytest.mark.parametrize("value", [0, 5, 10, 25, 50, 100, "auto", True, 1.0])
def test_unverified_native_campaign_choices_are_rejected(value):
    with pytest.raises(ValueError, match="supports only 1"):
        validate_campaign_config({"native_attempts_per_job": value})


@pytest.mark.parametrize("value", [2, True, 1.0])
def test_conflicting_dispatch_alias_is_rejected(value):
    with pytest.raises(ValueError, match="conflicts"):
        validate_campaign_config({"batch_size": 1, "dispatch_batch_size": value})


def test_capabilities_do_not_promote_source_reuse_into_measured_batch_support():
    rows = {row["engine_id"]: row for row in list_generator_capabilities()}
    json.dumps(rows)
    for engine in ("bindcraft2", "rfd4_proteina", "rfdiffusion", "rfdiffusion3"):
        native = rows[engine]["native_execution"]
        assert native["multi_attempt_invocation"] is True
        assert native["model_reuse_source_verified"] is True
        assert native["model_reuse_runtime_measured"] is False
        assert native["max_verified_campaign_attempts_per_job"] == (1 if rows[engine]["campaign_adapter"] else 0)
    for engine in ("boltzgen", "pxdesign", "genie3"):
        assert rows[engine]["native_execution"]["multi_attempt_invocation"] is None


def test_simulated_shared_job_keeps_partial_attempt_results_and_surplus_candidates(tmp_path):
    # Historical schema-v1 diagnostic metadata is not a native job identity.
    config = {"max_attempts": 3, "accepted_target": 1, "seed": 91}
    ledger = CampaignLedger.create(tmp_path / "ledger.json", config, campaign_id="simulated", schema_version=1)
    attempts = ledger.authorize(5, batch_id="dispatch-group")
    assert len(attempts) == 3
    for attempt in attempts:
        ledger.transition(attempt["attempt_id"], "running", execution_job_id="simulated-job")
    row = lambda name: {"candidate_id": name, "metrics": {"native_pass_filters": True}}
    first, second, third = [attempt["attempt_id"] for attempt in attempts]
    assert ledger.record_receipt(first, "receipt-1", candidates=[row("c1"), row("c2")],
                                 native_state={"trajectories": 1, "accepted": 2})
    ledger.transition(second, "failed", failure="partial job failure")
    ledger.transition(third, "interrupted", failure="container interrupted")
    assert not ledger.authorize(1, batch_id="extra")
    assert not ledger.record_receipt(first, "receipt-1", candidates=[row("c1")], native_state={})
    assert not ledger.record_receipt(first, "replayed-result", candidates=[row("c1")], native_state={})
    recovered = CampaignLedger.create(ledger.path, config)
    retry = recovered.retry(second)
    assert retry["seed"] == attempts[1]["seed"]
    assert len({a["seed"] for a in attempts}) == 3
    assert recovered.counters()["attempts_authorized"] == 3
    assert recovered.counters()["attempts_started"] == 3
    assert recovered.counters()["attempts_completed"] == 1
    assert recovered.counters()["native_accepted_candidates"] == 2
    assert recovered.counters()["failed_interrupted_attempts"] == 2
    state = recovered.snapshot()
    assert len({a["execution_job_id"] for a in state["attempts"].values()}) == 1
    assert state["candidates"]["c1"]["attempt_id"] == first
    assert state["candidates"]["c2"]["attempt_id"] == first


def test_simulated_cancelled_job_does_not_invent_started_attempts(tmp_path):
    ledger = CampaignLedger.create(tmp_path / "ledger.json", {"max_attempts": 2})
    attempts = ledger.authorize(2, batch_id="simulated-job")
    ledger.transition(attempts[0]["attempt_id"], "running")
    ledger.transition(attempts[0]["attempt_id"], "interrupted")
    ledger.transition(attempts[1]["attempt_id"], "not_started")
    ledger.terminate("cancelled")
    assert ledger.counters()["attempts_authorized"] == 2
    assert ledger.counters()["attempts_started"] == 1
    assert ledger.counters()["attempts_completed"] == 0
    assert ledger.authorize(1, batch_id="after-cancel") == []
