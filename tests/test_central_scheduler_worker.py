from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler_worker


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
