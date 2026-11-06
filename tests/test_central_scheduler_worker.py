from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler_worker
import json
import sys


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
