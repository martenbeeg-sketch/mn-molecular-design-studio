from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler, jobs
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import (
    AutomationError,
    validate_workflow_request,
    workflow_catalog,
    workflow_schema,
)
from mn_molecular_design_studio.job_catalog import _job_payload


def test_boltz2_workflow_is_registered_with_seed_and_job_name_fields():
    workflow_id = "refolding.optimized_prediction"
    catalog = workflow_catalog()
    assert any(row["workflow"] == workflow_id for row in catalog)

    schema = workflow_schema(workflow_id)
    parameters = {row["name"]: row for row in schema["parameters"]}
    assert parameters["engine"]["required"] is True
    assert "boltz2" in parameters["engine"]["type"]
    assert parameters["job_name"]["type"] == "str"
    assert parameters["seed"]["type"] == "int"


def test_boltz2_validation_preserves_user_job_name_random_seed_and_normal_priority(tmp_path: Path):
    request = {
        "schema_version": 1,
        "workflow": "refolding.optimized_prediction",
        "parameters": {
            "engine": "boltz2",
            "optimization_mode": "fast",
            "job_name": "boltz2-prediction",
            "seed": 0,
            "samples": 5,
            "tasks": [{"protein_sequences": "ACDEFGHIK"}],
            "engine_options": {"use_msa_server": True},
        },
        "priority": 0,
    }
    request_path = tmp_path / "workflow-request.json"
    request_path.write_text(json.dumps(request), encoding="utf-8")

    validated = validate_workflow_request(request_path)

    assert validated["valid"] is True
    assert validated["workflow"] == "refolding.optimized_prediction"
    assert validated["priority"] == 0
    assert validated["parameters"]["job_name"] == "boltz2-prediction"
    assert validated["parameters"]["seed"] == 0


def test_retired_workflow_id_is_rejected(tmp_path: Path):
    request_path = tmp_path / "stale-request.json"
    request_path.write_text(json.dumps({
        "schema_version": 1,
        "workflow": "refolding.optimized_prediction_retired",
        "parameters": {"engine": "boltz2"},
    }), encoding="utf-8")

    with pytest.raises(AutomationError, match="Unknown workflow"):
        validate_workflow_request(request_path)


@pytest.mark.parametrize(
    ("metadata_name", "params_name", "expected"),
    [
        ("boltz2-prediction", "fallback-name", "boltz2-prediction"),
        ("", "boltz2-prediction", "boltz2-prediction"),
    ],
)
def test_job_catalog_uses_user_job_name(
    tmp_path: Path, metadata_name: str, params_name: str, expected: str,
):
    run_dir = tmp_path / "refolding-validation" / "20261016-test-run"
    run_dir.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps({
        "status": "completed",
        "tool": "optimized_prediction_boltz2",
        "campaign_name": metadata_name,
    }), encoding="utf-8")
    (run_dir / "input.json").write_text(json.dumps({
        "params": {"campaign_name": params_name},
    }), encoding="utf-8")
    (run_dir / "result.json").write_text(json.dumps({"success": True}), encoding="utf-8")

    assert _job_payload(run_dir)["campaign_name"] == expected


def test_scheduler_submission_uses_normal_priority_without_campaign_override(monkeypatch, tmp_path: Path):
    runs_root = tmp_path / "runs"
    run_dir = runs_root / "refolding-validation" / "20261016-test-run"
    run_dir.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps({
        "run_id": run_dir.name,
        "task_group": run_dir.parent.name,
        "tool": "optimized_prediction_boltz2",
        "job_type": "complex_refolding",
        "status": "queued",
    }), encoding="utf-8")
    (run_dir / "input.json").write_text(json.dumps({"params": {}}), encoding="utf-8")
    (run_dir / "worker_request.json").write_text(json.dumps({"kind": "workflow_call"}), encoding="utf-8")

    class FakeSchedulerClient:
        submitted: dict[str, object] | None = None
        submission_count = 0

        def submit_job(self, **kwargs):
            self.submission_count += 1
            self.submitted = kwargs
            return {"job_id": "fake-scheduler-job", "status": "queued", "created_at": "2026-10-16T00:00:00Z"}

        def get_job(self, job_id):
            assert job_id == "fake-scheduler-job"
            return {"job_id": job_id, "status": "queued"}

    class FakeSchedulerAPIError(Exception):
        status = 500

    client = FakeSchedulerClient()
    monkeypatch.setattr(jobs, "runs_root", lambda: runs_root)
    monkeypatch.setattr(central_scheduler, "scheduler_mode", lambda: "central")
    monkeypatch.setattr(central_scheduler, "_central_resource_request", lambda _run_dir: {"gpu": True, "gpu_count": 1})
    monkeypatch.setattr(central_scheduler, "_client", lambda: client)
    monkeypatch.setattr(
        central_scheduler,
        "scheduler_client_module",
        lambda: SimpleNamespace(SchedulerAPIError=FakeSchedulerAPIError),
    )
    monkeypatch.setattr(central_scheduler, "ensure_scheduler_worker_service", lambda: None)

    with central_scheduler.submission_priority(0):
        central_scheduler.submit_run(run_dir)
        central_scheduler.submit_run(run_dir)

    assert client.submitted is not None
    assert client.submission_count == 1
    assert client.submitted["name"] == "optimized_prediction_boltz2"
    assert client.submitted["priority"] == 0
    assert client.submitted["workflow"] == central_scheduler.WORKFLOW_ID
    assert client.submitted["campaign_id"] == ""
    metadata = json.loads((run_dir / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["scheduler_backend"] == central_scheduler.BACKEND
    assert metadata["scheduler_priority"] == 0
    assert metadata["scheduler_job_id"] == "fake-scheduler-job"
