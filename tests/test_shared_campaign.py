"""Durable cross-generator accounting used by normal graph campaigns."""
from __future__ import annotations

import json
import copy
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from mn_molecular_design_studio.engine.workflow_engine.workflows import distributed_campaign
from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import CampaignLedger
from mn_molecular_design_studio.engine.workflow_engine.workflows.shared_campaign import SharedCampaignLedger


def shared_config(**overrides):
    config = {
        "enabled": True,
        "generator_nodes": ["bindcraft", "rfd4"],
        "target_node": "target",
        "objective_type": "accepted",
        "unique_accepted": True,
        "accepted_target": 2,
        "max_attempts": 4,
        "max_concurrent_jobs": 2,
        "max_concurrent_generation_jobs": 2,
        "max_concurrent_validations": 1,
        "generator_limits": {
            "bindcraft": {"max_attempts": 3, "max_concurrent_generation_jobs": 1},
            "rfd4": {"max_attempts": 3, "max_concurrent_generation_jobs": 1},
        },
        "validation_stages": [{"workflow": "refolding.optimized_prediction", "parameters": {"engine": "boltz2"},
                               "required_metrics": ["iptm"],
                               "rules": [{"field": "iptm", "operator": ">=", "value": 0.7}]}],
    }
    return {**config, **overrides}


def candidate(identifier, sequence):
    return {"candidate_id": identifier, "binder_sequence": sequence, "metrics": {"native_pass_filters": True},
            "raw_metadata": {"generator": identifier.split("-")[0]}}


def test_two_generators_share_unique_acceptance_after_rejection_and_stop_at_target(tmp_path):
    shared = SharedCampaignLedger.create(tmp_path / "shared.json", shared_config())
    assert shared.reserve_attempt("bc-campaign", "bc-a1", "bindcraft", "bindcraft2")
    assert shared.reserve_attempt("rfd-campaign", "rf-a1", "rfd4", "rfd4_proteina")
    # Global concurrency is held across both independent campaign ledgers.
    assert not shared.reserve_attempt("bc-campaign", "bc-a2", "bindcraft", "bindcraft2")

    first = candidate("bindcraft-one", "ACDEFGHIKLMNPQRSTVWY")
    claim = shared.claim_candidate("bc-campaign", "bindcraft", first, native_acceptance="accepted")
    assert claim["role"] == "owner"
    shared.finish_candidate("bc-campaign", first["candidate_id"], "rejected", "Boltz-2 iPTM 0.62 < 0.70")

    second = candidate("rfd4-one", "YYYYVVVVTTTTSSSSRRRR")
    assert shared.claim_candidate("rfd-campaign", "rfd4", second, native_acceptance="unavailable")["role"] == "owner"
    shared.finish_candidate("rfd-campaign", second["candidate_id"], "accepted", "Boltz-2 iPTM 0.81 >= 0.70")

    duplicate = candidate("rfd4-duplicate", " acdefghiklmnpqrstvwy ")
    duplicate_claim = shared.claim_candidate("rfd-campaign", "rfd4", duplicate, native_acceptance="unavailable")
    assert duplicate_claim["role"] == "duplicate"
    assert duplicate_claim["duplicate_of"] == first["candidate_id"]
    # A failed first measurement can be recovered by a later branch; a real
    # scientific rejection remains final and is never reevaluated.
    retry_sequence = candidate("bindcraft-two", "ACDEFGHIKLMNPQRSTVWY")
    assert shared.claim_candidate("bc-campaign", "bindcraft", retry_sequence, native_acceptance="accepted")["role"] == "duplicate"

    summary = shared.summary()
    assert summary["unique_accepted"] == 1
    assert summary["unique_rejected"] == 1
    assert summary["duplicate_candidates"] == 2
    assert summary["attempts_active"] == 2

    # Finish native attempts; one new unique candidate can complete the target.
    shared.sync_campaign("bc-campaign", "bindcraft", {"config": {"generator_id": "bindcraft2"}, "attempts": {
        "bc-a1": {"status": "completed"}}}, 4.0)
    shared.sync_campaign("rfd-campaign", "rfd4", {"config": {"generator_id": "rfd4_proteina"}, "attempts": {
        "rf-a1": {"status": "completed"}}}, 5.0)
    assert shared.reserve_attempt("bc-campaign", "bc-a2", "bindcraft", "bindcraft2")
    third = candidate("bindcraft-three", "MNPQRSTVWYACDEFGHIKL")
    assert shared.claim_candidate("bc-campaign", "bindcraft", third, native_acceptance="accepted")["role"] == "owner"
    shared.finish_candidate("bc-campaign", third["candidate_id"], "accepted", "common predictor criteria passed")
    assert shared.summary()["status"] == "target_reached"
    assert not shared.reserve_attempt("rfd-campaign", "rf-a2", "rfd4", "rfd4_proteina")


def test_shared_attempt_budget_and_interrupted_reservation_recovery(tmp_path):
    shared = SharedCampaignLedger.create(tmp_path / "shared.json", shared_config(
        accepted_target=3, max_attempts=2, generator_limits={
            "bindcraft": {"max_attempts": 2, "max_concurrent_generation_jobs": 1},
            "rfd4": {"max_attempts": 2, "max_concurrent_generation_jobs": 1},
        },
    ))
    assert shared.reserve_attempt("bc-campaign", "bc-a1", "bindcraft", "bindcraft2")
    assert shared.reserve_attempt("rfd-campaign", "rf-a1", "rfd4", "rfd4_proteina")
    assert not shared.reserve_attempt("rfd-campaign", "rf-a2", "rfd4", "rfd4_proteina")

    # Recovery sees no local authorization after interruption between the
    # shared reservation and the native campaign ledger write, so it releases
    # that reservation and permits a safe retry with a fresh identity.
    shared.sync_campaign("bc-campaign", "bindcraft", {
        "config": {"generator_id": "bindcraft2"}, "attempts": {}, "candidates": {},
    })
    state = shared.snapshot()
    assert state["attempts"]["bc-campaign:bc-a1"]["status"] == "released"
    assert shared.reserve_attempt("bc-campaign", "bc-a2", "bindcraft", "bindcraft2")

    shared.sync_campaign("bc-campaign", "bindcraft", {
        "config": {"generator_id": "bindcraft2"}, "attempts": {"bc-a2": {"status": "running"}}, "candidates": {},
    })
    reloaded = SharedCampaignLedger(tmp_path / "shared.json")
    assert reloaded.summary()["attempts_active"] == 2
    reloaded.sync_campaign("bc-campaign", "bindcraft", {
        "config": {"generator_id": "bindcraft2"}, "attempts": {"bc-a2": {"status": "completed"}}, "candidates": {},
    })
    reloaded.sync_campaign("rfd-campaign", "rfd4", {
        "config": {"generator_id": "rfd4_proteina"}, "attempts": {"rf-a1": {"status": "completed"}}, "candidates": {},
    })
    assert reloaded.summary()["status"] == "budget_exhausted"
    unevaluated = reloaded.claim_candidate(
        "rf-campaign", "rfd4", candidate("late-candidate", "ACDEFGHIKLMNPQRSTVWY"),
        native_acceptance="accepted",
    )
    assert unevaluated["role"] == "not_evaluated"


def test_concurrent_generator_controllers_cannot_oversubscribe_shared_job_limit(tmp_path):
    shared = SharedCampaignLedger.create(tmp_path / "shared.json", shared_config(
        accepted_target=20, max_attempts=20, max_concurrent_jobs=3,
        max_concurrent_generation_jobs=3, max_concurrent_validations=2,
        generator_limits={
            "bindcraft": {"max_attempts": 20, "max_concurrent_generation_jobs": 3},
            "rfd4": {"max_attempts": 20, "max_concurrent_generation_jobs": 3},
        },
    ))

    def reserve(index):
        node = "bindcraft" if index % 2 == 0 else "rfd4"
        generator = "bindcraft2" if node == "bindcraft" else "rfd4_proteina"
        return shared.reserve_attempt(f"campaign-{node}-{index}", f"attempt-{index}", node, generator)

    with ThreadPoolExecutor(max_workers=12) as pool:
        accepted_reservations = list(pool.map(reserve, range(12)))
    assert sum(accepted_reservations) == 3
    assert shared.summary()["attempts_active"] == 3


def test_failed_validation_can_be_reclaimed_after_interrupted_downstream_evaluation(tmp_path):
    shared = SharedCampaignLedger.create(tmp_path / "shared.json", shared_config())
    original = candidate("bc-design", "ACDEFGHIKLMNPQRSTVWY")
    assert shared.claim_candidate("bc-campaign", "bindcraft", original,
                                  native_acceptance="accepted")["role"] == "owner"
    shared.finish_candidate("bc-campaign", "bc-design", "failed", "prediction worker disappeared")

    recovered = candidate("rf-recovered", " acdefghiklmnpqrstvwy ")
    claim = shared.claim_candidate("rf-campaign", "rfd4", recovered,
                                   native_acceptance="unavailable")
    assert claim["role"] == "owner"
    assert shared.snapshot()["sequences"][claim["sequence_identity"]]["retry_claimed_at"]
    assert shared.reserve_validation("rf-campaign", "rfd4", "rf-recovered", 0)
    shared.sync_campaign("rf-campaign", "rfd4", {
        "config": {"generator_id": "rfd4_proteina"}, "attempts": {}, "candidates": {},
    })
    recovered_ledger = SharedCampaignLedger(tmp_path / "shared.json")
    assert recovered_ledger.snapshot()["validations"]["rf-campaign:rf-recovered:0"]["status"] == "released"
    assert recovered_ledger.summary()["unique_validation_failed"] == 0
    assert recovered_ledger.reserve_validation("rf-campaign", "rfd4", "rf-recovered", 0)
    recovered_ledger.sync_campaign("rf-campaign", "rfd4", {
        "config": {"generator_id": "rfd4_proteina"}, "attempts": {},
        "candidates": {"rf-recovered": {"validation_stages": {"0": {"status": "submitting"}}}},
    })
    assert recovered_ledger.snapshot()["validations"]["rf-campaign:rf-recovered:0"]["status"] == "submitting"
    assert recovered_ledger.reserve_validation("rf-campaign", "rfd4", "rf-recovered", 0)


def test_shared_validation_reservations_enforce_overall_concurrency_and_deduplicate_evaluation(tmp_path):
    shared = SharedCampaignLedger.create(tmp_path / "shared.json", shared_config(
        max_concurrent_jobs=1, max_concurrent_generation_jobs=1,
    ))
    assert shared.reserve_validation("bc-campaign", "bindcraft", "candidate-a", 0)
    assert not shared.reserve_validation("rf-campaign", "rfd4", "candidate-b", 0)
    shared.sync_campaign("bc-campaign", "bindcraft", {
        "config": {"generator_id": "bindcraft2"}, "attempts": {},
        "candidates": {"candidate-a": {"native_acceptance": "accepted", "validated_acceptance": "pending",
                                         "validation_stages": {"0": {"status": "running", "run_dir": "/tmp/prediction"}}}},
    })
    assert not shared.reserve_attempt("rf-campaign", "rf-a1", "rfd4", "rfd4_proteina")

    first = candidate("candidate-a", "ACDEFGHIKLMNPQRSTVWY")
    assert shared.claim_candidate("bc-campaign", "bindcraft", first, native_acceptance="accepted")["role"] == "owner"
    second = candidate("candidate-b", "acdefghiklmnpqrstvwy")
    duplicate = shared.claim_candidate("rf-campaign", "rfd4", second, native_acceptance="unavailable")
    assert duplicate["role"] == "duplicate"
    assert not shared.reserve_validation("rf-campaign", "rfd4", "candidate-b", 0)


def test_dispatch_attempts_obeys_graph_shared_limit_through_campaign_ledgers(tmp_path, monkeypatch):
    graph_parent = tmp_path / "graph-run"
    shared_config_value = shared_config(max_concurrent_jobs=1, max_concurrent_generation_jobs=1)
    shared_path = graph_parent / "artifacts/shared_campaign/ledger.json"
    shared = SharedCampaignLedger.create(shared_path, shared_config_value)
    submitted = []

    def campaign_parent(node_id, campaign_id, generator_id):
        parent = graph_parent / f"branch-{node_id}"
        parent.mkdir(parents=True, exist_ok=True)
        (parent / "input.json").write_text(json.dumps({"params": {}}))
        (parent / "metadata.json").write_text(json.dumps({
            "pipeline_parent_run_dir": str(graph_parent), "graph_node_id": node_id,
        }))
        config = {
            "enabled": True, "generator_id": generator_id,
            "accepted_target": 2, "max_attempts": 4,
            "acceptance_level": "validated", "unique_accepted": True,
            "max_concurrent_generation_jobs": 1,
            "validation_stages": shared_config_value["validation_stages"],
            "_shared_objective": {"ledger_path": str(shared_path), "generator_node": node_id},
        }
        ledger = CampaignLedger.create(parent / "artifacts/distributed_campaign/ledger.json", config, campaign_id=campaign_id)
        submitted.append((parent, ledger, distributed_campaign.validate_campaign_config(config)))

    campaign_parent("bindcraft", "bc-campaign", "bindcraft2")
    campaign_parent("rfd4", "rf-campaign", "rfd4_proteina")
    monkeypatch.setattr(distributed_campaign, "_create_attempt_child",
                        lambda _parent, _ledger, attempt, _config: attempt)
    first_parent, first_ledger, first_config = submitted[0]
    second_parent, second_ledger, second_config = submitted[1]
    distributed_campaign._dispatch_attempts(first_parent, first_ledger, first_config, active_count=0, accepted_count=0)
    distributed_campaign._dispatch_attempts(second_parent, second_ledger, second_config, active_count=0, accepted_count=0)

    assert first_ledger.counters()["attempts_authorized"] == 1
    assert second_ledger.counters()["attempts_authorized"] == 0
    assert shared.summary()["attempts_active"] == 1


def test_normal_campaign_controllers_continue_from_bindcraft_rejection_to_rfd4_acceptance(
    tmp_path, monkeypatch, design_source_factory,
):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    root = tmp_path / "runs"
    monkeypatch.setattr(jobs, "runs_root", lambda: root)
    monkeypatch.setattr("mn_molecular_design_studio.engine.workflow_engine.runtime.runs_root", lambda: root)
    graph_parent = root / "campaign" / "shared-parent"
    graph_parent.mkdir(parents=True)
    shared_options = shared_config(
        accepted_target=1, max_attempts=2,
        max_concurrent_jobs=2, max_concurrent_generation_jobs=2,
        max_concurrent_validations=1,
        generator_limits={
            "bindcraft": {"max_attempts": 1, "max_concurrent_generation_jobs": 1},
            "rfd4": {"max_attempts": 1, "max_concurrent_generation_jobs": 1},
        },
    )
    shared = SharedCampaignLedger.create(graph_parent / "artifacts/shared_campaign/ledger.json", shared_options)
    branch_specs = [
        ("bindcraft", "bindcraft2", "bc-campaign", "ACDEFG", 0.2),
        ("rfd4", "rfd4_proteina", "rf-campaign", "KLMNA", 0.8),
    ]
    parents = {}
    candidates = {}
    for node_id, generator_id, campaign_id, _sequence, _score in branch_specs:
        parent = root / "design" / f"branch-{node_id}"
        parent.mkdir(parents=True)
        config = {
            "enabled": True, "generator_id": generator_id,
            "accepted_target": 1, "max_attempts": 1,
            "acceptance_level": "validated", "unique_accepted": True,
            "max_concurrent_generation_jobs": 1, "max_concurrent_validations": 1,
            "validation_stages": shared_options["validation_stages"],
            "poll_interval_seconds": 0.25, "max_retries": 0,
            "_shared_objective": {"ledger_path": str(shared.path), "generator_node": node_id},
        }
        jobs.write_json(parent / "metadata.json", {
            "run_id": campaign_id, "task_group": "design", "tool": generator_id,
            "status": "running", "pipeline_parent_run_dir": str(graph_parent), "graph_node_id": node_id,
        })
        jobs.write_json(parent / "input.json", {"params": {"distributed_campaign": config}})
        local = CampaignLedger.create(
            distributed_campaign._campaign_ledger_path(parent), config,
            campaign_id=campaign_id, schema_version=1,
        )
        source_dir, candidate_row = design_source_factory(root / "fixtures" / node_id, missing_occupancy=False)
        candidate_row = copy.deepcopy(candidate_row)
        sequence = "KLMNP" if node_id == "bindcraft" else "KLMNA"
        candidate_row.update(candidate_id=f"{node_id}-design", binder_sequence=sequence,
                              source_tool=generator_id, metrics={"native_pass_filters": True})
        if node_id == "rfd4":
            # Keep the second design structurally sequence-consistent so the
            # normal validation input contract can identify the same binder.
            from Bio.PDB import MMCIFIO, MMCIFParser
            structure = MMCIFParser(QUIET=True).get_structure("rfd4", source_dir / "designed.cif")
            structure[0]["B"][5].resname = "ALA"
            output = MMCIFIO(); output.set_structure(structure); output.save(str(source_dir / "designed.cif"))
            candidate_row["complex_pdb"] = str(source_dir / "designed.cif")
        candidates[node_id] = candidate_row
        parents[node_id] = (parent, local, config)

    def create_attempt(parent, ledger, attempt, _config):
        child = parent / "artifacts/native_attempts" / attempt["attempt_id"]
        child.mkdir(parents=True)
        jobs.write_json(child / "metadata.json", {"status": "completed", "resource_allocation": {"gpu_device": "none"}})
        jobs.write_json(child / "result.json", {"success": True})
        ledger.transition(attempt["attempt_id"], "queued", child_run_dir=str(child))
        return child

    def reconcile(ledger, parent, _config):
        state = ledger.snapshot()
        children = {}
        for attempt_id, attempt in state["attempts"].items():
            child = Path(attempt["child_run_dir"]) if attempt.get("child_run_dir") else None
            if child is None:
                continue
            children[attempt_id] = child
            if attempt_id not in state["receipts"] and attempt.get("status") != "completed":
                node_id = jobs.read_json(parent / "metadata.json")["graph_node_id"]
                row = copy.deepcopy(candidates[node_id])
                ledger.record_receipt(attempt_id, f"receipt-{node_id}", candidates=[row],
                                      native_state={"accepted": 1, "trajectories": 1})
        latest = ledger.snapshot()
        rows = [record["candidate"] for record in latest["candidates"].values()]
        return children, [], rows

    validation_runs = []
    def submit_validation(parent, _ledger, candidate_row, _stage, stage_index, _previous):
        run_dir = parent / "artifacts/validation" / f"{candidate_row['candidate_id']}-{stage_index}"
        run_dir.mkdir(parents=True)
        jobs.write_json(run_dir / "metadata.json", {"status": "completed", "resource_allocation": {"gpu_device": "none"}})
        jobs.write_json(run_dir / "result.json", {"success": True})
        validation_runs.append((candidate_row["candidate_id"], candidate_row["binder_sequence"]))
        return run_dir

    def evidence(_run_dir, candidate_id, sequence):
        node_id = candidate_id.split("-", 1)[0]
        return {
            "candidate_id": f"boltz-{node_id}", "stage": "complex_refolding",
            "source_tool": "boltz2", "binder_sequence": sequence,
            "metrics": {"iptm": 0.2 if node_id == "bindcraft" else 0.8},
            "raw_metadata": {"measurement_status": {"status": "measured"},
                             "identity_check": {"sequence_chain_identity": "pass"}},
        }

    def finish(run_dir, success, result):
        jobs.write_json(run_dir / "result.json", {"success": success, **result})
        jobs.update_status(run_dir, "completed" if success else "failed")

    monkeypatch.setattr(distributed_campaign, "_create_attempt_child", create_attempt)
    monkeypatch.setattr(distributed_campaign, "_reconcile_attempts", reconcile)
    monkeypatch.setattr(distributed_campaign, "_submit_validation_stage", submit_validation)
    monkeypatch.setattr(distributed_campaign, "_validation_evidence_row", evidence)
    monkeypatch.setattr(distributed_campaign.time, "sleep", lambda _seconds: None)
    monkeypatch.setattr(jobs, "finish_job", finish)

    bindcraft_parent = parents["bindcraft"][0]
    distributed_campaign.run_bindcraft2_distributed_campaign(bindcraft_parent)
    rejected = CampaignLedger(distributed_campaign._campaign_ledger_path(bindcraft_parent)).snapshot()
    assert next(iter(rejected["candidates"].values()))["validated_acceptance"] == "rejected"
    assert shared.summary()["unique_rejected"] == 1
    assert shared.summary()["unique_accepted"] == 0

    rfd4_parent = parents["rfd4"][0]
    distributed_campaign.run_rfd4_proteina_distributed_campaign(rfd4_parent)
    assert shared.summary()["status"] == "target_reached"
    assert shared.summary()["unique_accepted"] == 1
    assert shared.summary()["attempts_authorized"] == 2
    assert [sequence for _candidate_id, sequence in validation_runs] == ["KLMNP", "KLMNA"]
    rfd4_state = CampaignLedger(distributed_campaign._campaign_ledger_path(rfd4_parent)).snapshot()
    accepted = next(iter(rfd4_state["candidates"].values()))
    assert accepted["validated_acceptance"] == "accepted"
    assert accepted["validation_stages"]["0"]["evidence_candidate"]["metrics"]["iptm"] == 0.8
