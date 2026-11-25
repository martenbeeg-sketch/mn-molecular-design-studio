import json

from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import _docker_cpu_quota_supported, queued_worker_runs


def _write_run(root, name, metadata):
    run_dir = root / "design" / name
    run_dir.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps(metadata))
    (run_dir / "worker_request.json").write_text("{}")
    return run_dir


def test_local_scheduler_ignores_distributed_proxy_runs(tmp_path):
    local_run = _write_run(tmp_path, "local", {"status": "queued"})
    cluster_proxy = _write_run(tmp_path, "cluster-proxy", {
        "status": "queued", "cluster_remote": True,
    })
    central_scheduler_run = _write_run(tmp_path, "central-scheduler", {
        "status": "queued", "scheduler_backend": "mn-compute-scheduler",
    })

    assert queued_worker_runs(tmp_path) == [local_run]
    assert cluster_proxy not in queued_worker_runs(tmp_path)
    assert central_scheduler_run not in queued_worker_runs(tmp_path)


def test_cgroup_v2_root_can_support_quotas_without_a_cpu_max_file(tmp_path):
    controllers = tmp_path / "cgroup.controllers"
    controllers.write_text("cpuset cpu io memory\n")
    assert _docker_cpu_quota_supported(tmp_path)
    controllers.write_text("cpuset io memory\n")
    assert not _docker_cpu_quota_supported(tmp_path)
