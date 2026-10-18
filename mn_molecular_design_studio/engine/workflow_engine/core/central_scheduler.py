from __future__ import annotations

"""Submit Molecular Design Studio runs to the shared MN Compute Scheduler."""

import os
import importlib
import subprocess
import sys
import threading
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from typing import Any, Iterator

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, utc_now, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import request_for_run
from mn_molecular_design_studio.engine.workflow_engine.runtime import app_home, runs_root


APP_ID = "mn-molecular-design-studio"
WORKFLOW_ID = "molecular_design_studio.run"
BACKEND = "mn-compute-scheduler"
SERVICE_DIR_NAME = "_mn_compute_scheduler_worker"
_submission_priority: ContextVar[int] = ContextVar("mnd_studio_scheduler_priority", default=0)
_scheduler_import_lock = threading.Lock()


def scheduler_client_module():
    """Load the sibling scheduler client even if a vendored package shadows it."""
    with _scheduler_import_lock:
        module = sys.modules.get("mn_compute_scheduler")
        if module is not None and hasattr(module, "SchedulerClient"):
            return module

        client_package_roots = []
        for entry in sys.path:
            package_root = Path(entry or os.getcwd()) / "mn_compute_scheduler"
            if (package_root / "client.py").is_file():
                client_package_roots.append(str(package_root.parent.resolve()))
        if client_package_roots:
            preferred = client_package_roots[0]
            sys.path[:] = [preferred, *(path for path in sys.path if path != preferred)]

        # The app also ships resource helpers under this top-level name when
        # its engine directory is on sys.path. Drop that module so the actual
        # scheduler client and its submodules load consistently from one root.
        module = sys.modules.get("mn_compute_scheduler")
        if module is not None and not hasattr(module, "SchedulerClient"):
            for name in tuple(sys.modules):
                if name == "mn_compute_scheduler" or name.startswith("mn_compute_scheduler."):
                    del sys.modules[name]
        return importlib.import_module("mn_compute_scheduler")


@contextmanager
def submission_priority(priority: int) -> Iterator[None]:
    """Apply a priority to the scheduler submission created in this call context."""
    if type(priority) is not int or not 0 <= priority <= 1000:
        raise ValueError("Scheduler priority must be an integer from 0 to 1000.")
    token = _submission_priority.set(priority)
    try:
        yield
    finally:
        _submission_priority.reset(token)


def scheduler_mode() -> str:
    return os.getenv("MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_MODE", "central").strip().lower() or "central"


def _client(*, timeout: float = 3.0):
    try:
        SchedulerClient = scheduler_client_module().SchedulerClient
    except ImportError as exc:
        raise RuntimeError(
            "The MN Compute Scheduler client is not installed. Install the sibling "
            "mn-compute-scheduler project into this app environment."
        ) from exc
    return SchedulerClient(timeout=timeout)


def _alive(pid: object) -> bool:
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


def scheduler_status() -> dict[str, Any]:
    client = _client(timeout=2.0)
    health = client.health()
    status = client.status()
    app_workers = [
        worker
        for worker in status.get("workers", [])
        if worker.get("app_id") == APP_ID and worker.get("online")
    ]
    return {
        "url": client.base_url,
        "health": health,
        "workers": app_workers,
    }


def _central_resource_request(run_dir: Path) -> dict[str, Any]:
    local = request_for_run(run_dir)
    request: dict[str, Any] = {
        "gpu": bool(local.gpu_resource),
        "cpu_threads": max(0, int(local.cpu_cores)),
        "min_vram_gb": float(local.min_vram_gb),
        "ram_gb": float(local.ram_gb),
        "scratch_gb": float(local.scratch_gb),
    }
    if not local.gpu_resource:
        return request

    raw = str(local.gpu_device or "all").strip().lower().removeprefix("device=")
    if raw in {"cpu", "none", "off", "false"}:
        request["gpu"] = False
        request["min_vram_gb"] = 0.0
        return request
    if raw in {"", "auto", "any"}:
        request["gpu_count"] = 1
    elif raw == "all":
        from mn_molecular_design_studio.engine.mn_compute_scheduler.resources import discover_gpu_ids

        ids = list(discover_gpu_ids())
        request["gpu_count"] = len(ids) or 1
        if ids:
            request["gpu_ids"] = ids
    else:
        try:
            ids = [int(part.strip()) for part in raw.split(",") if part.strip()]
        except ValueError as exc:
            raise ValueError(f"Invalid GPU selector in queued run: {local.gpu_device!r}") from exc
        if not ids:
            request["gpu_count"] = 1
        else:
            request["gpu_ids"] = ids
            request["gpu_count"] = len(ids)
    return request


def _service_dir() -> Path:
    # Worker identity must not move when the Results directory changes in
    # Settings. Keep the singleton lock/log in the stable Studio workspace.
    return app_home() / SERVICE_DIR_NAME


def ensure_scheduler_worker_service() -> None:
    # Settings may have changed since this API process started. Refresh the
    # environment before launching a worker so it inherits the current paths.
    from mn_molecular_design_studio.workspace import ensure_studio_workspace

    ensure_studio_workspace()
    service_dir = _service_dir()
    service_dir.mkdir(parents=True, exist_ok=True)
    owner_path = service_dir / "service.lock" / "owner.json"
    owner = read_json(owner_path)
    if _alive(owner.get("pid")):
        return
    log_path = service_dir / "worker.log"
    with log_path.open("a", encoding="utf-8") as log:
        subprocess.Popen(
            [sys.executable, "-m", "mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler_worker", "--service"],
            # The engine directory contains a bundled resource package also
            # named ``mn_compute_scheduler``. Starting there puts it at
            # sys.path[0] and shadows the external SchedulerClient package
            # supplied through PYTHONPATH. Launch from the project root so
            # the installed/shared scheduler client resolves first.
            cwd=str(Path(__file__).resolve().parents[4]),
            stdout=log,
            stderr=log,
            start_new_session=True,
            close_fds=True,
        )


def submit_run(run_dir: Path) -> dict[str, Any]:
    """Idempotently submit a prepared app run and wake this app's scheduler workers."""
    run_dir = run_dir.expanduser().resolve()
    metadata = read_json(run_dir / "metadata.json")
    if scheduler_mode() == "local":
        from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import ensure_worker_service

        ensure_worker_service()
        return {"backend": "local-worker-service", "job_id": ""}

    priority = _submission_priority.get()

    # Prevent the legacy local scanner from racing a job accepted by the central queue.
    metadata.update(
        {
            "scheduler_backend": BACKEND,
            "scheduler_priority": priority,
            "scheduler_dispatch_requested": True,
            "scheduler_dispatch_requested_at": utc_now(),
            "updated_at": utc_now(),
        }
    )
    write_json(run_dir / "metadata.json", metadata)

    client = _client()
    SchedulerAPIError = scheduler_client_module().SchedulerAPIError

    existing_id = str(metadata.get("scheduler_job_id") or "").strip()
    if existing_id:
        try:
            existing = client.get_job(existing_id)
        except SchedulerAPIError as exc:
            if exc.status != 404:
                metadata = read_json(run_dir / "metadata.json")
                metadata["scheduler_error"] = f"Could not confirm existing scheduler job: {exc}"
                metadata["scheduler_wait_reason"] = "Could not safely check the existing scheduler submission."
                metadata["updated_at"] = utc_now()
                write_json(run_dir / "metadata.json", metadata)
                raise RuntimeError(f"Could not confirm existing scheduler job {existing_id}: {exc}") from exc
            existing = {}
        if existing and str(existing.get("status") or "") in {"queued", "running"}:
            ensure_scheduler_worker_service()
            return existing

    generation = int(metadata.get("scheduler_submission_generation") or 0) + 1
    resource_request = _central_resource_request(run_dir)
    payload = {
        "run_dir": str(run_dir),
        "request_kind": str((read_json(run_dir / "worker_request.json")).get("kind") or ""),
    }
    try:
        job = client.submit_job(
            app_id=APP_ID,
            workflow=WORKFLOW_ID,
            run_id=run_dir.name,
            run_dir=str(run_dir),
            name=str(metadata.get("tool") or metadata.get("job_type") or run_dir.name),
            payload=payload,
            resource_request=resource_request,
            priority=priority,
            idempotency_key=f"{APP_ID}:{run_dir.parent.name}:{run_dir.name}:{generation}",
            campaign_id=str(metadata.get("campaign_id") or ""),
            max_attempts=2,
        )
    except Exception as exc:
        metadata = read_json(run_dir / "metadata.json")
        metadata["scheduler_error"] = f"{type(exc).__name__}: {exc}"
        metadata["scheduler_wait_reason"] = "Could not submit this run to MN Compute Scheduler."
        metadata["updated_at"] = utc_now()
        write_json(run_dir / "metadata.json", metadata)
        raise RuntimeError(f"Could not submit this run to MN Compute Scheduler: {exc}") from exc

    metadata = read_json(run_dir / "metadata.json")
    metadata.update(
        {
            "scheduler_backend": BACKEND,
            "scheduler_priority": priority,
            "scheduler_job_id": str(job.get("job_id") or ""),
            "scheduler_submission_generation": generation,
            "scheduler_resource_request": resource_request,
            "scheduler_status": str(job.get("status") or "queued"),
            "scheduler_submitted_at": str(job.get("created_at") or utc_now()),
            "scheduler_wait_reason": str(job.get("wait_reason") or ""),
            "scheduler_error": "",
            "updated_at": utc_now(),
        }
    )
    write_json(run_dir / "metadata.json", metadata)
    ensure_scheduler_worker_service()
    return job


def cancel_scheduled_run(run_dir: Path) -> None:
    metadata = read_json(run_dir / "metadata.json")
    job_id = str(metadata.get("scheduler_job_id") or "").strip()
    if metadata.get("scheduler_backend") != BACKEND or not job_id:
        return
    try:
        _client().cancel_job(job_id)
    except Exception as exc:
        metadata["scheduler_error"] = f"Could not cancel scheduler job: {exc}"
        metadata["updated_at"] = utc_now()
        write_json(run_dir / "metadata.json", metadata)
