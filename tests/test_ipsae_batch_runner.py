import importlib.util
import json
import sys
import threading
from pathlib import Path


RUNNER_PATH = Path(__file__).parents[1] / "mn_molecular_design_studio/engine/workflow_engine/workflows/ipsae_batch_runner.py"
SPEC = importlib.util.spec_from_file_location("mn_ipsae_batch_runner", RUNNER_PATH)
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def test_ipsae_batch_runner_parallelizes_independent_candidates_with_bounded_workers(tmp_path, monkeypatch):
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    status_path = work_dir / "batch_status.json"
    manifest_path = work_dir / "manifest.json"
    manifest_path.write_text(json.dumps({
        "work_dir": str(work_dir),
        "status_path": str(status_path),
        "jobs": [{"id": f"candidate-{index}", "command": ["fake"]} for index in range(3)],
    }))
    barrier = threading.Barrier(3)
    active = 0
    maximum_active = 0
    lock = threading.Lock()

    def run_one(job, _work_dir):
        nonlocal active, maximum_active
        with lock:
            active += 1
            maximum_active = max(maximum_active, active)
        barrier.wait(timeout=3)
        with lock:
            active -= 1
        return {"id": job["id"], "return_code": 0}

    monkeypatch.setattr(RUNNER, "_run_one", run_one)
    monkeypatch.setattr(sys, "argv", ["ipsae_batch_runner.py", "--manifest", str(manifest_path), "--workers", "3"])

    assert RUNNER.main() == 0
    status = json.loads(status_path.read_text())
    assert status["workers"] == 3
    assert maximum_active == 3
    assert [job["id"] for job in status["jobs"]] == ["candidate-0", "candidate-1", "candidate-2"]
