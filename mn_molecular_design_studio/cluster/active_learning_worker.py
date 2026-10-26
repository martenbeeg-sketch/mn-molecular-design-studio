"""Remote launcher for one scheduled active-learning docking shard.

This small file is copied to a selected workstation for the duration of a
round. The Uni-Dock Pro workflow itself is queued through that workstation's
normal Studio scheduler, so its CPU and GPU allocations remain visible to the
local scheduler and job table.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def run(request_path: Path) -> dict:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import enqueue_workflow_call
    from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import run_ligand_unidock_pro
    from mn_molecular_design_studio.workspace import ensure_studio_workspace, mark_studio_job

    request = json.loads(Path(request_path).read_text(encoding="utf-8"))
    if request.get("schema_version") != "studio.active-learning-docking-shard.v1":
        raise ValueError("Unsupported active-learning docking shard request.")
    # Remote launchers run in a clean SSH shell, outside the long-lived Studio
    # API process that normally sets these paths. Apply the machine's saved
    # workspace settings before create_job() chooses the run directory.
    ensure_studio_workspace()
    params = {**(request.get("parameters") or {}), "import_to_database": False}
    run_dir = enqueue_workflow_call(
        run_ligand_unidock_pro,
        task_group="ligand-docking",
        tool="unidock_pro",
        job_type="ligand_docking",
        kwargs=params,
        cpu_cores=int(params["cpu_cores"]),
        gpu_device=str(params["gpu_device"]),
    )
    mark_studio_job({"job": {"run_dir": str(run_dir)}},
                    owner_id=str(request["owner_id"]), project_id=str(request["project_id"]))
    metadata = read_json(run_dir / "metadata.json")
    metadata.update(
        active_learning_parent_run=str(request["parent_run_dir"]),
        active_learning_round=int(request["round_number"]),
        active_learning_worker_node=str(request["node_id"]),
        active_learning_worker_name=str(request["node_name"]),
    )
    write_json(run_dir / "metadata.json", metadata)

    cancel_marker = Path(request_path).parent / "cancel"
    cancellation_sent = False
    while True:
        if cancel_marker.exists() and not cancellation_sent:
            from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
            stop_job(run_dir, reason="Active-learning campaign was cancelled by the organizer")
            cancellation_sent = True
        metadata = read_json(run_dir / "metadata.json")
        status = str(metadata.get("status") or "queued")
        scheduler_status = str(metadata.get("scheduler_status") or "").casefold()
        scheduler_error = str(metadata.get("scheduler_error") or "").strip()
        if scheduler_status == "failed" and status not in {"completed", "failed", "cancelled"}:
            return {
                "status": "failed",
                "run_dir": str(run_dir),
                "job_code": str(metadata.get("job_code") or ""),
                "error": scheduler_error or "The machine scheduler rejected the docking job.",
                "score_file": "artifacts/docking_scores.csv",
            }
        if status in {"completed", "failed", "cancelled"}:
            result = read_json(run_dir / "result.json")
            return {
                "status": status,
                "run_dir": str(run_dir),
                "job_code": str(metadata.get("job_code") or ""),
                "error": str(result.get("error") or metadata.get("worker_error") or ""),
                "score_file": "artifacts/docking_scores.csv",
            }
        time.sleep(2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    args = parser.parse_args()
    try:
        payload = run(args.request)
    except BaseException as exc:
        payload = {"status": "failed", "error": f"{type(exc).__name__}: {exc}"}
    print("MDS_ACTIVE_DOCKING_RESULT=" + json.dumps(payload, separators=(",", ":")), flush=True)
    return 0 if payload.get("status") == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
