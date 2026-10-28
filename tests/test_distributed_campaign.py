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


def test_retry_preserves_attempt_identity_and_seed(tmp_path):
    campaign = ledger(tmp_path)
    attempt = campaign.authorize(1, batch_id="b1")[0]
    campaign.transition(attempt["attempt_id"], "failed", failure="worker")
    retried = campaign.retry(attempt["attempt_id"])
    assert retried["attempt_id"] == attempt["attempt_id"]
    assert retried["seed"] == attempt["seed"]
    assert retried["retry_count"] == 1
    assert campaign.counters()["attempts_authorized"] == 1


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
            "resources": {"gpu_device": "auto", "cpu_cores": 3},
        },
        0, "campaign-test", None,
    )
    assert payload["workflow"] == "refolding.optimized_prediction"
    assert payload["resources"] == {"cpu_cores": 3, "gpu_device": "auto"}
    task = payload["parameters"]["tasks"][0]
    assert task["target_protein_sequences"] == "ACDEFG"
    assert task["binder_sequence"] == "KLMNP"
    assert task["source_candidate_id"] == row["candidate_id"]
    assert task["source_run_dir"] == str(parent)
    assert task["molecule_chain_options"] == [
        {"id": "A", "use_msa": True}, {"id": "B", "use_msa": False},
    ]


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
