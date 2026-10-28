import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import (
    generator_capabilities,
    list_generator_capabilities,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import (
    CampaignLedger,
    validate_campaign_config,
)


def test_generator_capability_registry_describes_attempt_and_science_semantics():
    rows = {row["engine_id"]: row for row in list_generator_capabilities()}
    assert rows["rfd4_proteina"]["sequence_structure_mode"] == "co_design"
    assert rows["rfd4_proteina"]["native_acceptance_available"] is False
    assert rows["rfd4_proteina"]["attempt_unit"] == "one native design request"
    assert rows["rfdiffusion"]["sequence_structure_mode"] == "backbone_only"
    assert rows["rfdiffusion"]["campaign_adapter"] is False
    assert "protein_mpnn" in rows["rfdiffusion"]["downstream_validation"]
    assert rows["bindcraft2"]["native_acceptance_available"] is True
    assert generator_capabilities("RFD4_Proteina").engine_id == "rfd4_proteina"


def test_rfd4_campaign_rejects_native_acceptance_and_missing_external_criteria():
    with pytest.raises(ValueError, match="does not support 'native'"):
        validate_campaign_config({"generator_id": "rfd4_proteina"})
    with pytest.raises(ValueError, match="no authoritative native acceptance"):
        validate_campaign_config({
            "generator_id": "rfd4_proteina",
            "acceptance_level": "validated",
            "validation_stages": [{"workflow": "refolding.optimized_prediction"}],
        })


def test_rfdiffusion_campaign_is_not_exposed_without_sequence_handoff_adapter():
    with pytest.raises(ValueError, match="cannot run an acceptance-driven campaign"):
        validate_campaign_config({"generator_id": "rfdiffusion", "acceptance_level": "validated"})


def test_rfd4_validated_campaign_requires_external_measurements_and_gpu_memory():
    config = validate_campaign_config({
        "generator_id": "rfd4_proteina",
        "acceptance_level": "validated",
        "max_attempts": 2,
        "accepted_target": 1,
        "validation_stages": [{
            "workflow": "refolding.optimized_prediction",
            "parameters": {"engine": "boltz2"},
            "rules": [{"field": "model_1_iptm", "operator": ">=", "value": 0.5}],
        }],
    })
    assert config["generator_id"] == "rfd4_proteina"
    assert config["resources"]["min_vram_gb"] == 16.0
    assert config["acceptance_level"] == "validated"


def test_rfd4_receipts_record_acceptance_unavailable_and_still_allow_external_validation(tmp_path):
    config = validate_campaign_config({
        "generator_id": "rfd4_proteina",
        "acceptance_level": "validated",
        "max_attempts": 1,
        "accepted_target": 1,
        "validation_stages": [{
            "workflow": "refolding.optimized_prediction",
            "required_metrics": ["model_1_iptm"],
        }],
    })
    ledger = CampaignLedger.create(tmp_path / "campaign.json", config, campaign_id="rfd4-test")
    attempt = ledger.authorize(1, batch_id="b1")[0]
    row = {
        "candidate_id": "rfd4-one",
        "binder_sequence": "ACDEFG",
        "metrics": {"generated_residue_count": 6},
    }
    assert ledger.record_receipt(
        attempt["attempt_id"], "receipt-1", candidates=[row],
        native_state={"authoritative": False, "requested_designs": 1, "output_candidates": 1},
    )
    saved = ledger.snapshot()["candidates"]["rfd4-one"]
    assert saved["native_acceptance"] == "unavailable"
    assert ledger.counters()["native_accepted_candidates"] == 0
    assert ledger.counters()["native_acceptance_unavailable_candidates"] == 1
    assert ledger.counters()["native_rejected_trajectories"] == 0
    ledger.set_candidate_acceptance("rfd4-one", validated="accepted", reason="External evidence passed.")
    assert ledger.counters()["validated_accepted_candidates"] == 1
