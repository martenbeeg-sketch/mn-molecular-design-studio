from __future__ import annotations

"""Worker bridge from MN Compute Scheduler claims to app-owned run folders."""

import argparse
import os
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Any
from uuid import uuid4

from mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler import (
    APP_ID,
    BACKEND,
    SERVICE_DIR_NAME,
    WORKFLOW_ID,
    scheduler_client_module,
)
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, utc_now, write_json
from mn_molecular_design_studio.engine.workflow_engine.runtime import app_home, runs_root


HEARTBEAT_SECONDS = max(5.0, float(os.getenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_HEARTBEAT_SECONDS", "10")))
API_FAILURE_STOP_SECONDS = max(20.0, float(os.getenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_FAILURE_GRACE_SECONDS", "40")))
CPU_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
)


def _client():
    return scheduler_client_module().SchedulerClient(timeout=4.0)


def _service_dir() -> Path:
    # Match central_scheduler.py: this path is stable even if Results changes.
    return app_home() / SERVICE_DIR_NAME


def _refresh_workspace_paths() -> Path:
    """Reload the saved Studio directories and return the active Results root."""
    from mn_molecular_design_studio.workspace import ensure_studio_workspace

    ensure_studio_workspace()
    return runs_root().resolve()


def _process_alive(pid: object) -> bool:
    try:
        value = int(pid or 0)
    except (TypeError, ValueError):
        return False
    if value <= 0:
        return False
    try:
        os.kill(value, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _acquire_service_lock() -> bool:
    service_dir = _service_dir()
    service_dir.mkdir(parents=True, exist_ok=True)
    lock_dir = service_dir / "service.lock"
    owner_path = lock_dir / "owner.json"
    for _ in range(5):
        try:
            lock_dir.mkdir()
        except FileExistsError:
            owner = read_json(owner_path)
            if _process_alive(owner.get("pid")):
                return False
            time.sleep(0.2)
            owner = read_json(owner_path)
            if _process_alive(owner.get("pid")):
                return False
            shutil.rmtree(lock_dir, ignore_errors=True)
            continue
        write_json(owner_path, {"pid": os.getpid(), "started_at": utc_now()})
        return True
    return False


def _allocated_gpu_device(job: dict[str, Any]) -> str:
    allocation = job.get("allocation") if isinstance(job.get("allocation"), dict) else {}
    raw = allocation.get("gpu_ids") if isinstance(allocation.get("gpu_ids"), list) else []
    gpu_ids = [str(value) for value in raw]
    if not gpu_ids and allocation.get("gpu_id") is not None:
        gpu_ids = [str(allocation["gpu_id"])]
    return ",".join(gpu_ids) if gpu_ids else "cpu"


def _patch_gpu_selectors(value: Any, gpu_device: str) -> Any:
    if isinstance(value, dict):
        result: dict[str, Any] = {}
        for key, child in value.items():
            normalized = str(key).casefold()
            if normalized == "gpu_device" or normalized.endswith("_gpu_device") or normalized in {"capacity_device", "child_device"}:
                previous = str(child or "").strip().casefold()
                if gpu_device != "cpu" and previous not in {"none", "cpu", "off", "false", "0-gpu"}:
                    result[key] = gpu_device
                else:
                    result[key] = child
            else:
                result[key] = _patch_gpu_selectors(child, gpu_device)
        return result
    if isinstance(value, list):
        return [_patch_gpu_selectors(child, gpu_device) for child in value]
    return value


def _effective_resource_allocation(job: dict[str, Any]) -> dict[str, Any]:
    allocation = job.get("allocation") if isinstance(job.get("allocation"), dict) else {}
    gpu_device = _allocated_gpu_device(job)
    return {
        "scheduler": BACKEND,
        "scheduler_job_id": str(job.get("job_id") or ""),
        "cpu_cores": int(allocation.get("cpu_threads") or 0),
        "cpu_slots": allocation.get("cpu_slots") if isinstance(allocation.get("cpu_slots"), list) else [],
        "gpu_device": gpu_device,
        "gpu_devices": [part for part in gpu_device.split(",") if part and part != "cpu"],
        "gpu_ids": allocation.get("gpu_ids") if isinstance(allocation.get("gpu_ids"), list) else [],
        "ram_gb": float(allocation.get("ram_gb") or 0.0),
        "scratch_gb": float(allocation.get("scratch_gb") or 0.0),
        "allocated_at": utc_now(),
    }


def _apply_allocation(run_dir: Path, job: dict[str, Any]) -> dict[str, Any]:
    allocation = _effective_resource_allocation(job)
    gpu_device = str(allocation["gpu_device"])

    input_path = run_dir / "input.json"
    input_payload = read_json(input_path)
    patched_input = _patch_gpu_selectors(input_payload, gpu_device)
    if patched_input != input_payload:
        write_json(input_path, patched_input)

    request_path = run_dir / "worker_request.json"
    request = read_json(request_path)
    patched_request = _patch_gpu_selectors(request, gpu_device)
    if patched_request != request:
        write_json(request_path, patched_request)

    metadata = read_json(run_dir / "metadata.json")
    metadata.update(
        {
            "scheduler_backend": BACKEND,
            "scheduler_job_id": str(job.get("job_id") or ""),
            "scheduler_status": "running",
            "scheduler_attempt": int(job.get("attempt") or 1),
            "resource_allocation": allocation,
            "resource_request": job.get("resource_request") if isinstance(job.get("resource_request"), dict) else {},
            "worker_started_at": utc_now(),
            "scheduler_wait_reason": "",
            "updated_at": utc_now(),
        }
    )
    write_json(run_dir / "metadata.json", metadata)
    return allocation


def _terminate(process: subprocess.Popen, *, grace_seconds: float = 20.0) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=grace_seconds)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        return
    process.wait(timeout=5)


def _run_claimed_job(client, job: dict[str, Any], worker_id: str) -> None:
    job_id = str(job.get("job_id") or "")
    lease_token = str(job.get("lease_token") or "")
    payload = job.get("payload") if isinstance(job.get("payload"), dict) else {}
    run_dir_text = str(payload.get("run_dir") or job.get("run_dir") or "").strip()
    run_dir = Path(run_dir_text).expanduser().resolve()
    base = runs_root().resolve()

    def finish_failure(message: str, *, traceback_text: str = "") -> None:
        if run_dir.is_relative_to(base) and run_dir.is_dir():
            metadata = read_json(run_dir / "metadata.json")
            metadata.update({"scheduler_backend": BACKEND, "scheduler_job_id": job_id, "updated_at": utc_now()})
            if str(metadata.get("status") or "") not in {"cancelled", "paused"}:
                metadata.update({"scheduler_status": "failed", "scheduler_error": message})
            if traceback_text:
                metadata["scheduler_traceback"] = traceback_text[-12000:]
            write_json(run_dir / "metadata.json", metadata)
        client.finish_job(job_id, worker_id=worker_id, lease_token=lease_token, success=False, error=message)

    if not job_id or not lease_token:
        finish_failure("Scheduler claim is missing its job ID or lease token.")
        return
    if not run_dir.is_relative_to(base) or not run_dir.is_dir():
        finish_failure(f"Scheduled run folder is missing or outside Molecular Design Studio Results: {run_dir}")
        return
    if not (run_dir / "worker_request.json").is_file():
        finish_failure(f"Scheduled run has no worker_request.json: {run_dir}")
        return
    saved_metadata = read_json(run_dir / "metadata.json")
    if str(saved_metadata.get("scheduler_job_id") or "") not in {"", job_id}:
        finish_failure("Scheduler job ID does not match the app run metadata.")
        return
    if str(saved_metadata.get("status") or "") in {"cancelled", "paused"}:
        finish_failure(f"Studio run is {saved_metadata['status']}; it will not be started.")
        return

    try:
        allocation = _apply_allocation(run_dir, job)
        resource_request = job.get("resource_request") if isinstance(job.get("resource_request"), dict) else {}
        cpu_cores = max(1, int(resource_request.get("cpu_threads") or allocation.get("cpu_cores") or 1))
        env = os.environ.copy()
        for key in CPU_THREAD_ENV_VARS:
            env[key] = str(cpu_cores)
        env["MN_MOLECULAR_DESIGN_STUDIO_JOB_CPU_CORES"] = str(cpu_cores)
        env["MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_JOB_ID"] = job_id
        stdout = (run_dir / "worker_stdout.log").open("a", encoding="utf-8")
        stderr = (run_dir / "worker_stderr.log").open("a", encoding="utf-8")
        try:
            process = subprocess.Popen(
                [
                    sys.executable,
                    "-m",
                    "mn_molecular_design_studio.engine.workflow_engine.core.local_worker",
                    "--run-dir",
                    str(run_dir),
                ],
                cwd=str(Path(__file__).resolve().parents[2]),
                stdout=stdout,
                stderr=stderr,
                env=env,
                start_new_session=True,
                close_fds=True,
            )
        finally:
            stdout.close()
            stderr.close()
        metadata = read_json(run_dir / "metadata.json")
        metadata.update({"worker_pid": process.pid, "updated_at": utc_now()})
        write_json(run_dir / "metadata.json", metadata)
    except Exception as exc:
        finish_failure(f"{type(exc).__name__}: {exc}", traceback_text=traceback.format_exc())
        return

    last_heartbeat = 0.0
    heartbeat_failure_started: float | None = None
    cancellation_requested = False
    while process.poll() is None:
        time.sleep(1.0)
        now = time.monotonic()
        if now - last_heartbeat < HEARTBEAT_SECONDS:
            continue
        last_heartbeat = now
        try:
            client.heartbeat_worker(app_id=APP_ID, worker_id=worker_id, pid=os.getpid())
            state = client.heartbeat_job(
                job_id,
                worker_id=worker_id,
                lease_token=lease_token,
                pid=process.pid,
            )
            heartbeat_failure_started = None
            if state.get("cancel_requested"):
                cancellation_requested = True
                _terminate(process)
                break
        except Exception as exc:
            if heartbeat_failure_started is None:
                heartbeat_failure_started = now
            if now - heartbeat_failure_started >= API_FAILURE_STOP_SECONDS:
                _terminate(process)
                break

    return_code = process.poll()
    if return_code is None:
        return_code = process.wait(timeout=30)

    metadata = read_json(run_dir / "metadata.json")
    result = read_json(run_dir / "result.json")
    app_status = str(metadata.get("status") or "")
    if cancellation_requested:
        try:
            from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job

            stop_job(run_dir, "Cancelled through MN Compute Scheduler")
        except Exception:
            pass
        metadata = read_json(run_dir / "metadata.json")
        metadata.update({"scheduler_status": "cancelled", "updated_at": utc_now()})
        write_json(run_dir / "metadata.json", metadata)
        client.finish_job(
            job_id,
            worker_id=worker_id,
            lease_token=lease_token,
            success=False,
            error="Cancelled through MN Compute Scheduler.",
        )
        return

    if heartbeat_failure_started is not None and time.monotonic() - heartbeat_failure_started >= API_FAILURE_STOP_SECONDS:
        error = "Scheduler API heartbeats failed; stopped the worker to avoid running without a valid resource lease."
        finish_failure(error)
        return

    # The Studio result files are the durable record of workflow completion.
    # A cleanup/launcher exit after those files were committed must not turn a
    # successful scientific run into a failed scheduler job.
    success = app_status == "completed" and result.get("success") is True
    error = ""
    if not success:
        error = str(result.get("error") or metadata.get("worker_error") or "")
        if not error:
            if app_status == "completed":
                error = (
                    "Studio recorded completed metadata but did not record a successful result "
                    f"(worker exit status {return_code})."
                )
            else:
                error = (
                    f"Studio workflow status was {app_status or 'missing'} "
                    f"(worker exit status {return_code})."
                )
    metadata = read_json(run_dir / "metadata.json")
    metadata.update({"scheduler_status": "completed" if success else "failed", "updated_at": utc_now()})
    metadata["scheduler_worker_exit_code"] = return_code
    metadata.pop("scheduler_error", None)
    metadata.pop("scheduler_warning", None)
    if error:
        metadata["scheduler_error"] = error
    elif return_code != 0:
        metadata["scheduler_warning"] = (
            "Studio wrote a successful result before its worker exited "
            f"with status {return_code}."
        )
    write_json(run_dir / "metadata.json", metadata)
    scheduler_result = {
        "run_dir": str(run_dir),
        "result_path": str(run_dir / "result.json"),
        "metadata_path": str(run_dir / "metadata.json"),
        "app_status": app_status,
        "worker_exit_code": return_code,
        "success": success,
    }
    client.finish_job(
        job_id,
        worker_id=worker_id,
        lease_token=lease_token,
        success=success,
        result=scheduler_result,
        error=error,
    )


def _worker_loop(worker_index: int, stop: threading.Event) -> None:
    results_root = _refresh_workspace_paths()
    client = _client()
    worker_id = f"{socket.gethostname()}-{os.getpid()}-{worker_index}-{uuid4().hex[:8]}"

    def register_worker(root: Path) -> None:
        client.register_worker(
            app_id=APP_ID,
            worker_id=worker_id,
            workflows=[WORKFLOW_ID],
            pid=os.getpid(),
            hostname=socket.gethostname(),
            labels={"app_home": str(app_home()), "results_root": str(root), "worker": "molecular-design-studio"},
        )

    while not stop.is_set():
        try:
            register_worker(results_root)
            break
        except Exception:
            stop.wait(HEARTBEAT_SECONDS)
    if stop.is_set():
        return
    last_worker_heartbeat = 0.0
    while not stop.is_set():
        now = time.monotonic()
        if now - last_worker_heartbeat >= HEARTBEAT_SECONDS:
            try:
                client.heartbeat_worker(app_id=APP_ID, worker_id=worker_id, pid=os.getpid())
                last_worker_heartbeat = now
            except Exception:
                stop.wait(HEARTBEAT_SECONDS)
                continue
        try:
            # Read the persisted Results directory before each claim. This
            # covers Settings changes made while the bridge is already alive.
            refreshed_root = _refresh_workspace_paths()
            if refreshed_root != results_root:
                register_worker(refreshed_root)
                results_root = refreshed_root
            job = client.claim_job(app_id=APP_ID, worker_id=worker_id)
        except Exception:
            stop.wait(HEARTBEAT_SECONDS)
            continue
        if not job:
            stop.wait(max(0.5, float(os.getenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_POLL_SECONDS", "2"))))
            continue
        try:
            _run_claimed_job(client, job, worker_id)
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
            try:
                client.finish_job(
                    str(job.get("job_id") or ""),
                    worker_id=worker_id,
                    lease_token=str(job.get("lease_token") or ""),
                    success=False,
                    error=error,
                )
            except Exception:
                pass


def configured_worker_count() -> int:
    """Return the number of scheduler bridge workers to register."""
    worker_count = int(os.getenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_WORKERS", "6"))
    return max(1, min(worker_count, 32))


def run_service() -> int:
    results_root = _refresh_workspace_paths()
    service_dir = _service_dir()
    service_dir.mkdir(parents=True, exist_ok=True)
    if not _acquire_service_lock():
        return 0
    print(
        f"Molecular Design Studio scheduler bridge started: app_home={app_home()} results_root={results_root}",
        flush=True,
    )
    stop = threading.Event()
    previous_handlers: dict[int, Any] = {}

    def request_stop(signum, frame) -> None:
        stop.set()

    for signum in (signal.SIGTERM, signal.SIGINT):
        previous_handlers[signum] = signal.signal(signum, request_stop)
    try:
        # Pipeline coordinator jobs wait for GPU-backed child workflows. Keep
        # enough bridge workers available for those children while several
        # coordinators are active; the scheduler still enforces CPU/GPU leases.
        worker_count = configured_worker_count()
        threads = [
            threading.Thread(target=_worker_loop, args=(index + 1, stop), name=f"mds-scheduler-{index + 1}")
            for index in range(worker_count)
        ]
        for thread in threads:
            thread.start()
        while any(thread.is_alive() for thread in threads):
            for thread in threads:
                thread.join(timeout=1.0)
    finally:
        stop.set()
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
        shutil.rmtree(service_dir / "service.lock", ignore_errors=True)
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Molecular Design Studio workers for MN Compute Scheduler.")
    parser.add_argument("--service", action="store_true", required=True)
    args = parser.parse_args()
    if args.service:
        raise SystemExit(run_service())


if __name__ == "__main__":
    main()
