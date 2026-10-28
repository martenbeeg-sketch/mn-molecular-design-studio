from concurrent.futures import ThreadPoolExecutor

import pytest

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
    assert campaign.record_receipt(attempt["attempt_id"], "receipt-1", candidates=[row], native_state={"trajectories": 1})
    assert not campaign.record_receipt(attempt["attempt_id"], "receipt-1", candidates=[row], native_state={"trajectories": 1})
    assert not campaign.record_receipt(attempt["attempt_id"], "receipt-2", candidates=[row], native_state={"trajectories": 1})
    counters = campaign.counters()
    assert counters["attempts_completed"] == 1
    assert counters["native_accepted_candidates"] == 1
    assert counters["validated_accepted_candidates"] == 0
    campaign.set_candidate_acceptance(row["candidate_id"], validated="accepted", diversity="pending")
    counters = campaign.counters()
    assert counters["native_accepted_candidates"] == 1
    assert counters["validated_accepted_candidates"] == 1
    assert counters["diversity_qualified_accepted_candidates"] == 0


def test_retry_preserves_attempt_identity_and_seed(tmp_path):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="b1")[0]
    campaign.transition(attempt["attempt_id"], "failed", failure="worker")
    retried = campaign.retry(attempt["attempt_id"])
    assert retried["attempt_id"] == attempt["attempt_id"]
    assert retried["seed"] == attempt["seed"]
    assert retried["retry_count"] == 1
    assert campaign.counters()["attempts_authorized"] == 1


def test_budget_termination_is_explicit_and_stable(tmp_path):
    campaign = ledger(tmp_path, max_attempts=1, target=1)
    campaign.authorize(1, batch_id="b1")
    state = campaign.terminate("budget_exhausted", detail="one trajectory rejected")
    assert state["status"] == "completed"
    assert state["termination_reason"] == "budget_exhausted"
    assert campaign.authorize(1, batch_id="b2") == []
    assert campaign.terminate("target_reached")["termination_reason"] == "budget_exhausted"


@pytest.mark.parametrize(
    "config",
    [
        {"accepted_target": 0, "max_attempts": 2},
        {"accepted_target": 3, "max_attempts": 2},
        {"accepted_target": 1, "max_attempts": 1, "acceptance_level": "shortlist"},
        {"accepted_target": 1, "max_attempts": 1, "diversity": {"sequence_identity_max": 0}},
    ],
)
def test_invalid_campaign_policy_is_rejected(tmp_path, config):
    with pytest.raises(ValueError):
        CampaignLedger.create(tmp_path / "campaign.json", config)
