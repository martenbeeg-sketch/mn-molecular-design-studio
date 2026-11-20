from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler_worker
import json
import sys
from types import SimpleNamespace

import pytest


def test_scheduler_bridge_keeps_capacity_for_pipeline_child_jobs(monkeypatch):
    monkeypatch.delenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_WORKERS", raising=False)
    assert central_scheduler_worker.configured_worker_count() == 6


def test_scheduler_bridge_worker_count_can_be_configured_and_is_bounded(monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_WORKERS", "4")
    assert central_scheduler_worker.configured_worker_count() == 4
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_WORKERS", "50")
    assert central_scheduler_worker.configured_worker_count() == 32
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_WORKERS", "0")
    assert central_scheduler_worker.configured_worker_count() == 1


def test_local_preview_cannot_start_or_register_a_central_bridge(monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_MODE', 'local')
    def forbidden(*args, **kwargs):
        raise AssertionError('Local preview must not launch or register a central worker')
    monkeypatch.setattr(central_scheduler.subprocess, 'Popen', forbidden)
    monkeypatch.setattr(central_scheduler_worker, '_refresh_workspace_paths', forbidden)
    monkeypatch.setattr(central_scheduler_worker, '_client', forbidden)
    central_scheduler.ensure_scheduler_worker_service()
    assert central_scheduler_worker.run_service() == 0


def test_queued_python_and_workspace_override_a_stale_demo_bridge(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import _submission_environment
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "real"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "real"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "results"))
    submitted = _submission_environment()
    (tmp_path / "command.json").write_text(json.dumps({"mode":"local_worker", "command":[
        sys.executable,"-m","mn_molecular_design_studio.engine.workflow_engine.core.local_worker"],
        "environment":{**submitted,"UNRELATED_SECRET":"must not be propagated"}}))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", "/tmp/stale-demo/home")
    executable, environment = central_scheduler_worker._submitted_execution_context(tmp_path)
    assert executable == sys.executable
    assert environment["MN_MOLECULAR_DESIGN_STUDIO_APP_HOME"] == str(tmp_path / "real")
    assert environment["MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR"] == str(tmp_path / "results")
    assert "UNRELATED_SECRET" not in environment


def test_direct_job_descriptive_python_uses_the_worker_interpreter(tmp_path):
    path = tmp_path / "command.json"
    command = {"mode": "local_worker", "command": ["python", "-m",
        "mn_molecular_design_studio.engine.workflow_engine.core.local_worker", "--run-dir", str(tmp_path)]}
    path.write_text(json.dumps(command))
    assert central_scheduler_worker._submitted_execution_context(tmp_path) == (sys.executable, {})

    command["command"][0] = "other-python"
    path.write_text(json.dumps(command))
    with pytest.raises(ValueError, match="Python executable is unavailable"):
        central_scheduler_worker._submitted_execution_context(tmp_path)


@pytest.fixture
def bridge_run(tmp_path, monkeypatch):
    from mn_molecular_design_studio import job_dependencies
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    run = tmp_path / "analysis" / "test"
    run.mkdir(parents=True)
    jobs.write_json(run / "worker_request.json", {"kind": "workflow_call"})
    jobs.write_json(run / "metadata.json", {"status": "running", "scheduler_job_id": "job"})
    state = {"time": 100.0, "return_code": None, "finished": [], "cleaned": False}
    process = SimpleNamespace(pid=99999999, poll=lambda: state["return_code"])
    monkeypatch.setattr(central_scheduler_worker, "runs_root", lambda: tmp_path)
    monkeypatch.setattr(central_scheduler_worker, "_apply_allocation", lambda *args: {})
    monkeypatch.setattr(central_scheduler_worker, "_submitted_execution_context", lambda *args: (sys.executable, {}))
    monkeypatch.setattr(central_scheduler_worker.subprocess, "Popen", lambda *args, **kwargs: process)
    monkeypatch.setattr(central_scheduler_worker.time, "monotonic", lambda: state["time"])
    monkeypatch.setattr(central_scheduler_worker.time, "sleep", lambda seconds: state.update(time=state["time"] + 10))
    monkeypatch.setattr(central_scheduler_worker, "_terminate", lambda current, **kwargs: state.update(return_code=-15))
    monkeypatch.setattr(job_dependencies, "cancellation_requested", lambda *args: False)
    def cleanup(path, *, strict=False):
        assert strict and path == run and state["return_code"] == -15
        state["cleaned"] = True
        return ["owned-container"]
    monkeypatch.setattr(jobs, "_stop_run_docker_containers", cleanup)
    client = SimpleNamespace(finish_job=lambda *args, **kwargs: state["finished"].append(kwargs))
    job = {"job_id": "job", "lease_token": "private-lease-token", "payload": {"run_dir": str(run)}}
    return run, state, client, job


def test_worker_liveness_failure_cannot_block_healthy_job_lease(bridge_run):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    run, state, client, job = bridge_run
    state["renewals"] = 0
    def renew(*args, **kwargs):
        state["renewals"] += 1
        if state["renewals"] == 3:
            jobs.write_json(run / "metadata.json", {"status": "completed"})
            jobs.write_json(run / "result.json", {"success": True})
            state["return_code"] = 0
        return {}
    def failed_worker(**kwargs):
        raise TimeoutError("worker registry unavailable")
    client.heartbeat_job, client.heartbeat_worker = renew, failed_worker
    central_scheduler_worker._run_claimed_job(client, job, "worker")
    assert state["renewals"] == 3 and not state["cleaned"]
    assert state["finished"][0]["success"] is True


def test_lost_lease_stops_owned_container_and_records_real_terminal_state(bridge_run):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    run, state, client, job = bridge_run
    def failed_lease(*args, **kwargs):
        raise TimeoutError("lease timeout private-lease-token")
    client.heartbeat_job = failed_lease
    central_scheduler_worker._run_claimed_job(client, job, "worker")
    metadata = jobs.read_json(run / "metadata.json")
    assert state["cleaned"] and state["finished"][0]["success"] is False
    assert metadata["status"] == "failed" and metadata["scheduler_status"] == "failed"
    assert metadata["scheduler_stopped_containers"] == ["owned-container"]
    assert "TimeoutError" in metadata["scheduler_heartbeat_error"]
    assert "private-lease-token" not in json.dumps(metadata)
    assert jobs.read_json(run / "result.json")["success"] is False


def test_successful_durable_result_wins_over_late_heartbeat_error(bridge_run):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    run, state, client, job = bridge_run
    def late_error(*args, **kwargs):
        jobs.write_json(run / "metadata.json", {"status": "completed"})
        jobs.write_json(run / "result.json", {"success": True})
        state["return_code"] = 0
        raise TimeoutError("response was lost after completion")
    client.heartbeat_job = late_error
    central_scheduler_worker._run_claimed_job(client, job, "worker")
    assert state["finished"][0]["success"] is True and not state["cleaned"]


def test_near_expiry_lease_stops_before_full_failure_grace(bridge_run):
    from datetime import datetime, timedelta, timezone

    run, state, client, job = bridge_run
    job["lease_expires_at"] = (datetime.now(timezone.utc) + timedelta(seconds=10)).isoformat()
    def failed_lease(*args, **kwargs):
        raise TimeoutError("lease renewal unavailable")
    client.heartbeat_job = failed_lease
    central_scheduler_worker._run_claimed_job(client, job, "worker")
    assert state["cleaned"] and state["time"] == 110.0


def test_unconfirmed_container_cleanup_is_not_a_voluntary_resource_release(bridge_run, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    run, state, client, job = bridge_run
    def unavailable(*args, **kwargs):
        raise TimeoutError("unavailable")
    client.heartbeat_job = unavailable
    monkeypatch.setattr(jobs, "_stop_run_docker_containers", unavailable)
    central_scheduler_worker._run_claimed_job(client, job, "worker")
    assert state["finished"] == []
    metadata = jobs.read_json(run / "metadata.json")
    assert metadata["status"] == "failed" and metadata["scheduler_cleanup_error"]


def test_heartbeat_diagnostics_cannot_overwrite_concurrent_workflow_completion(bridge_run, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    run, state, client, job = bridge_run
    original_write = central_scheduler_worker.write_json
    def complete_during_diagnostic_write(path, value):
        if path.name == "scheduler_heartbeat.json":
            original_write(run / "metadata.json", {"status": "completed"})
            original_write(run / "result.json", {"success": True})
            state["return_code"] = 0
        original_write(path, value)
    monkeypatch.setattr(central_scheduler_worker, "write_json", complete_during_diagnostic_write)
    client.heartbeat_job = lambda *args, **kwargs: {}
    client.heartbeat_worker = lambda **kwargs: None
    central_scheduler_worker._run_claimed_job(client, job, "worker")
    assert state["finished"][0]["success"] is True
    assert jobs.read_json(run / "metadata.json")["status"] == "completed"
    assert jobs.read_json(run / "scheduler_heartbeat.json")["scheduler_heartbeat_failures"] == 0


def test_docker_cleanup_matches_exact_mount_not_path_substring(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    commands = []
    def docker(command, **kwargs):
        commands.append(command)
        if command[1] == "ps":
            return SimpleNamespace(returncode=0, stdout="owned\nneighbor\n")
        if command[1] == "inspect":
            source = str(tmp_path) if command[2] == "owned" else str(tmp_path) + "-other-run"
            return SimpleNamespace(returncode=0, stdout=json.dumps([{"Mounts": [{"Source": source}]}]))
        return SimpleNamespace(returncode=0, stdout="")
    monkeypatch.setattr(jobs.subprocess, "run", docker)
    assert jobs._stop_run_docker_containers(tmp_path, strict=True) == ["owned"]
    assert [command[-1] for command in commands if command[1] == "stop"] == ["owned"]
