from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Sequence
from uuid import uuid4


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text())
        return payload if isinstance(payload, dict) else {}
    except (OSError, ValueError, TypeError):
        return {}


def _atomic_write(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    temporary.replace(path)


def _pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _heartbeat_is_stale(payload: dict[str, object], *, stale_after_seconds: float) -> bool:
    raw = str(payload.get("heartbeat_at") or payload.get("acquired_at") or "")
    try:
        heartbeat = datetime.fromisoformat(raw)
        if heartbeat.tzinfo is None:
            heartbeat = heartbeat.replace(tzinfo=timezone.utc)
    except ValueError:
        return True
    age = (datetime.now(timezone.utc) - heartbeat).total_seconds()
    return age > stale_after_seconds and not _pid_is_alive(int(payload.get("pid") or 0))


def gpu_ids_from_command(command: Sequence[str]) -> tuple[int, ...] | None:
    """Return explicit Docker GPU IDs, or None when `all`/unspecified."""
    try:
        value = str(command[command.index("--gpus") + 1]).strip().lower()
    except (ValueError, IndexError):
        return None
    if value == "all":
        return None
    value = value.removeprefix("device=")
    parts = tuple(part.strip() for part in value.split(",") if part.strip())
    if not parts or any(not part.isdigit() for part in parts):
        raise ValueError(f"Invalid Docker GPU request: {value!r}")
    return tuple(dict.fromkeys(int(part) for part in parts))


def select_gpu_in_command(command: Sequence[str], gpu_id: int) -> list[str]:
    selected = [str(value) for value in command]
    try:
        index = selected.index("--gpus")
    except ValueError:
        return selected
    if index + 1 >= len(selected):
        raise ValueError("Docker --gpus option is missing its value")
    selected[index + 1] = f"device={int(gpu_id)}"
    return selected


def discover_gpu_ids(
    *,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> tuple[int, ...]:
    configured = os.getenv("MN_COMPUTE_SCHEDULER_GPU_IDS", "").strip()
    if not configured:
        configured = os.getenv("MN_LIGAND_GPU_IDS", "").strip()
    if configured:
        parts = tuple(part.strip() for part in configured.split(",") if part.strip())
        if not parts or any(not part.isdigit() for part in parts):
            raise ValueError(
                "MN_COMPUTE_SCHEDULER_GPU_IDS (or legacy MN_LIGAND_GPU_IDS) must be "
                "a comma-separated list of non-negative integers"
            )
        return tuple(dict.fromkeys(int(part) for part in parts))
    try:
        process = runner(
            ["nvidia-smi", "--query-gpu=index", "--format=csv,noheader,nounits"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ()
    if process.returncode != 0:
        return ()
    ids = tuple(int(line.strip()) for line in process.stdout.splitlines() if line.strip().isdigit())
    return tuple(dict.fromkeys(ids))


@dataclass(frozen=True)
class GPUCapacity:
    gpu_id: int
    free_vram_gb: float
    total_vram_gb: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "gpu_id": self.gpu_id,
            "free_vram_gb": round(self.free_vram_gb, 3),
            "total_vram_gb": round(self.total_vram_gb, 3),
        }


@dataclass(frozen=True)
class ResourceSnapshot:
    cpu_threads_total: int
    ram_available_gb: float
    ram_total_gb: float
    scratch_free_gb: float
    scratch_total_gb: float
    gpus: tuple[GPUCapacity, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "cpu_threads_total": self.cpu_threads_total,
            "ram_available_gb": round(self.ram_available_gb, 3),
            "ram_total_gb": round(self.ram_total_gb, 3),
            "scratch_free_gb": round(self.scratch_free_gb, 3),
            "scratch_total_gb": round(self.scratch_total_gb, 3),
            "gpus": [gpu.to_dict() for gpu in self.gpus],
        }


@dataclass(frozen=True)
class AdmissionRequest:
    gpu: bool = False
    min_vram_gb: float = 0.0
    cpu_threads: int = 0
    ram_gb: float = 0.0
    scratch_gb: float = 0.0

    @classmethod
    def from_dict(cls, payload: dict[str, Any], *, gpu_default: bool = False) -> AdmissionRequest:
        try:
            gpu_value = payload.get("gpu", gpu_default)
            if not isinstance(gpu_value, bool) and gpu_value not in {0, 1}:
                raise ValueError("gpu must be a boolean")
            request = cls(
                gpu=bool(gpu_value),
                min_vram_gb=float(payload.get("min_vram_gb") or 0.0),
                cpu_threads=int(payload.get("cpu_threads") or 0),
                ram_gb=float(payload.get("ram_gb") or 0.0),
                scratch_gb=float(payload.get("scratch_gb") or 0.0),
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Resource values must be numeric: {exc}") from exc
        if request.cpu_threads < 0 or min(
            request.min_vram_gb, request.ram_gb, request.scratch_gb
        ) < 0:
            raise ValueError("Resource values cannot be negative")
        if not all(
            math.isfinite(value)
            for value in (request.min_vram_gb, request.ram_gb, request.scratch_gb)
        ):
            raise ValueError("Resource values must be finite")
        if not request.gpu and request.min_vram_gb:
            raise ValueError("CPU jobs cannot request VRAM")
        return request

    def to_dict(self) -> dict[str, Any]:
        return {
            "gpu": self.gpu,
            "min_vram_gb": self.min_vram_gb,
            "cpu_threads": self.cpu_threads,
            "ram_gb": self.ram_gb,
            "scratch_gb": self.scratch_gb,
        }


@dataclass(frozen=True)
class AdmissionDecision:
    allowed: bool
    permanent: bool = False
    reasons: tuple[str, ...] = ()
    eligible_gpu_ids: tuple[int, ...] = ()


def _memory_capacity(meminfo_path: Path = Path("/proc/meminfo")) -> tuple[float, float]:
    values: dict[str, int] = {}
    try:
        for line in meminfo_path.read_text().splitlines():
            key, raw = line.split(":", 1)
            token = raw.strip().split()[0]
            values[key] = int(token) * 1024
    except (OSError, ValueError, IndexError):
        return 0.0, 0.0
    divisor = float(1024**3)
    total = values.get("MemTotal", 0) / divisor
    available = values.get("MemAvailable", values.get("MemFree", 0)) / divisor
    return available, total


def discover_gpu_capacity(
    *,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> tuple[GPUCapacity, ...]:
    try:
        process = runner(
            [
                "nvidia-smi",
                "--query-gpu=index,memory.free,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ()
    if process.returncode != 0:
        return ()
    capacities: list[GPUCapacity] = []
    for line in process.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != 3:
            continue
        try:
            gpu_id, free_mib, total_mib = int(parts[0]), float(parts[1]), float(parts[2])
        except ValueError:
            continue
        capacities.append(
            GPUCapacity(
                gpu_id=gpu_id,
                free_vram_gb=free_mib / 1024.0,
                total_vram_gb=total_mib / 1024.0,
            )
        )
    return tuple(capacities)


def capture_resource_snapshot(
    scratch_path: Path,
    *,
    gpu_runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> ResourceSnapshot:
    try:
        cpu_threads = len(os.sched_getaffinity(0))
    except AttributeError:
        cpu_threads = int(os.cpu_count() or 1)
    ram_available, ram_total = _memory_capacity()
    disk = shutil.disk_usage(scratch_path)
    divisor = float(1024**3)
    return ResourceSnapshot(
        cpu_threads_total=max(1, cpu_threads),
        ram_available_gb=ram_available,
        ram_total_gb=ram_total,
        scratch_free_gb=disk.free / divisor,
        scratch_total_gb=disk.total / divisor,
        gpus=discover_gpu_capacity(runner=gpu_runner),
    )


def assess_resource_admission(
    request: AdmissionRequest,
    snapshot: ResourceSnapshot,
    *,
    candidate_gpu_ids: Sequence[int] = (),
) -> AdmissionDecision:
    permanent: list[str] = []
    waiting: list[str] = []
    if request.cpu_threads > snapshot.cpu_threads_total:
        permanent.append(
            f"requires {request.cpu_threads} CPU threads; host capacity is {snapshot.cpu_threads_total}"
        )
    if snapshot.ram_total_gb and request.ram_gb > snapshot.ram_total_gb:
        permanent.append(
            f"requires {request.ram_gb:g} GiB RAM; host capacity is {snapshot.ram_total_gb:.2f} GiB"
        )
    elif request.ram_gb > snapshot.ram_available_gb:
        waiting.append(
            f"requires {request.ram_gb:g} GiB RAM; {snapshot.ram_available_gb:.2f} GiB is available"
        )
    if request.scratch_gb > snapshot.scratch_total_gb:
        permanent.append(
            f"requires {request.scratch_gb:g} GiB scratch; filesystem capacity is {snapshot.scratch_total_gb:.2f} GiB"
        )
    elif request.scratch_gb > snapshot.scratch_free_gb:
        waiting.append(
            f"requires {request.scratch_gb:g} GiB scratch; {snapshot.scratch_free_gb:.2f} GiB is free"
        )

    eligible_gpu_ids: tuple[int, ...] = ()
    if request.gpu:
        candidates = tuple(dict.fromkeys(int(value) for value in candidate_gpu_ids))
        if not candidates:
            waiting.append("no requested GPU is available in this worker's device scope")
        elif request.min_vram_gb <= 0:
            eligible_gpu_ids = candidates
        else:
            by_id = {gpu.gpu_id: gpu for gpu in snapshot.gpus}
            observed = [by_id[gpu_id] for gpu_id in candidates if gpu_id in by_id]
            eligible_gpu_ids = tuple(
                gpu.gpu_id for gpu in observed if gpu.free_vram_gb >= request.min_vram_gb
            )
            if not eligible_gpu_ids:
                if observed and len(observed) == len(candidates) and all(
                    gpu.total_vram_gb < request.min_vram_gb for gpu in observed
                ):
                    permanent.append(
                        f"requires {request.min_vram_gb:g} GiB VRAM; requested GPU capacity is insufficient"
                    )
                elif not observed:
                    waiting.append("GPU VRAM availability could not be measured")
                else:
                    free = ", ".join(
                        f"GPU {gpu.gpu_id}: {gpu.free_vram_gb:.2f} GiB free" for gpu in observed
                    )
                    waiting.append(f"requires {request.min_vram_gb:g} GiB VRAM; {free}")

    reasons = tuple(permanent + waiting)
    return AdmissionDecision(
        allowed=not reasons,
        permanent=bool(permanent),
        reasons=reasons,
        eligible_gpu_ids=eligible_gpu_ids,
    )


@dataclass
class FileLease:
    path: Path
    run_id: str
    owner_id: str
    payload: dict[str, object]

    @classmethod
    def acquire(
        cls,
        path: Path,
        *,
        run_id: str,
        owner_id: str,
        extra: dict[str, object] | None = None,
        stale_after_seconds: float = 120.0,
    ) -> FileLease | None:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload: dict[str, object] = {
            "run_id": run_id,
            "owner_id": owner_id,
            "pid": os.getpid(),
            "acquired_at": _utc_now_iso(),
            "heartbeat_at": _utc_now_iso(),
            **(extra or {}),
        }
        for attempt in range(2):
            try:
                descriptor = os.open(str(path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError:
                existing = _read_json(path)
                if attempt == 0 and _heartbeat_is_stale(existing, stale_after_seconds=stale_after_seconds):
                    try:
                        path.unlink()
                    except FileNotFoundError:
                        pass
                    continue
                return None
            with os.fdopen(descriptor, "w") as handle:
                handle.write(json.dumps(payload, indent=2) + "\n")
            return cls(path=path, run_id=run_id, owner_id=owner_id, payload=payload)
        return None

    def heartbeat(self) -> None:
        current = _read_json(self.path)
        if str(current.get("run_id") or "") != self.run_id or str(current.get("owner_id") or "") != self.owner_id:
            raise RuntimeError(f"Lease ownership changed: {self.path}")
        current["heartbeat_at"] = _utc_now_iso()
        _atomic_write(self.path, current)
        self.payload = current

    def set_owner_pid(self, pid: int) -> None:
        current = _read_json(self.path)
        if str(current.get("run_id") or "") != self.run_id or str(current.get("owner_id") or "") != self.owner_id:
            raise RuntimeError(f"Lease ownership changed: {self.path}")
        current["pid"] = int(pid)
        current["heartbeat_at"] = _utc_now_iso()
        _atomic_write(self.path, current)
        self.payload = current

    def release(self) -> bool:
        current = _read_json(self.path)
        if str(current.get("run_id") or "") != self.run_id or str(current.get("owner_id") or "") != self.owner_id:
            return False
        try:
            self.path.unlink()
        except FileNotFoundError:
            return False
        return True


@dataclass
class CPULease:
    """A single job's reservation from the shared CPU-slot pool."""

    leases: tuple[FileLease, ...]
    capacity: int

    @property
    def threads(self) -> int:
        return len(self.leases)

    @property
    def slot_ids(self) -> tuple[int, ...]:
        return tuple(int(lease.payload.get("cpu_slot") or 0) for lease in self.leases)

    def heartbeat(self) -> None:
        for lease in self.leases:
            lease.heartbeat()

    def set_owner_pid(self, pid: int) -> None:
        for lease in self.leases:
            lease.set_owner_pid(pid)

    def release(self) -> bool:
        released = [lease.release() for lease in self.leases]
        return bool(released) and all(released)


def shared_state_dir() -> Path:
    configured = os.getenv("MN_COMPUTE_SCHEDULER_STATE_DIR", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    xdg_state = os.getenv("XDG_STATE_HOME", "").strip()
    base = Path(xdg_state).expanduser() if xdg_state else Path.home() / ".local" / "state"
    return (base / "mn-compute-scheduler").resolve()


def worker_state_dir(runs_dir: Path | None = None, *, state_root: Path | None = None) -> Path:
    """Return a worker state directory, optionally shared by all MN apps."""
    if state_root is not None:
        return Path(state_root).expanduser().resolve() / ".worker"
    if runs_dir is not None:
        return Path(runs_dir) / ".worker"
    return shared_state_dir() / ".worker"


def acquire_job_claim(
    run_dir: Path,
    *,
    run_id: str,
    worker_id: str,
    stale_after_seconds: float = 120.0,
) -> FileLease | None:
    return FileLease.acquire(
        run_dir / ".worker-claim.json",
        run_id=run_id,
        owner_id=worker_id,
        stale_after_seconds=stale_after_seconds,
    )


def cpu_pool_capacity(
    host_cpu_threads: int | None = None,
    *,
    fallback_capacity: int | None = None,
) -> int:
    """Return one host-wide CPU pool size, constrained by the process affinity."""
    raw = os.getenv("MN_COMPUTE_SCHEDULER_CPU_SLOTS", "").strip()
    if raw:
        try:
            configured = int(raw)
        except ValueError as exc:
            raise ValueError("MN_COMPUTE_SCHEDULER_CPU_SLOTS must be a positive integer") from exc
    elif fallback_capacity is not None:
        configured = int(fallback_capacity)
    else:
        try:
            visible = len(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            visible = os.cpu_count() or 1
        configured = max(1, visible - 1)
    if configured < 1:
        raise ValueError("Scheduler CPU slots must be at least one")
    if host_cpu_threads is None:
        try:
            host_cpu_threads = len(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            host_cpu_threads = os.cpu_count() or 1
    return max(1, min(configured, int(host_cpu_threads)))


def acquire_cpu_lease(
    threads: int,
    *,
    run_id: str,
    worker_id: str,
    workflow: str = "",
    app_id: str = "",
    runs_dir: Path | None = None,
    state_root: Path | None = None,
    capacity: int | None = None,
    stale_after_seconds: float = 120.0,
) -> CPULease | None:
    """Atomically reserve CPU slots shared by CPU and GPU workers.

    A short allocator lease prevents workers from each retaining a partial
    reservation. Individual slot leases remain held and heartbeated for the
    complete job lifetime.
    """

    requested = max(1, int(threads))
    pool_capacity = cpu_pool_capacity() if capacity is None else max(1, int(capacity))
    if requested > pool_capacity:
        return None
    lock_dir = worker_state_dir(runs_dir, state_root=state_root) / "locks"
    allocator = FileLease.acquire(
        lock_dir / "cpu-pool-allocator.json",
        run_id=run_id,
        owner_id=worker_id,
        extra={"workflow": workflow, "requested_cpu_threads": requested, "app_id": app_id},
        stale_after_seconds=stale_after_seconds,
    )
    if allocator is None:
        return None
    acquired: list[FileLease] = []
    try:
        for slot_id in range(pool_capacity):
            lease = FileLease.acquire(
                lock_dir / f"cpu-{slot_id}.json",
                run_id=run_id,
                owner_id=worker_id,
                extra={
                    "cpu_slot": slot_id,
                    "cpu_pool_capacity": pool_capacity,
                    "requested_cpu_threads": requested,
                    "workflow": workflow,
                    "app_id": app_id,
                },
                stale_after_seconds=stale_after_seconds,
            )
            if lease is not None:
                acquired.append(lease)
                if len(acquired) == requested:
                    return CPULease(tuple(acquired), pool_capacity)
        for lease in acquired:
            lease.release()
        return None
    finally:
        allocator.release()


def acquire_gpu_lease(
    gpu_id: int,
    *,
    run_id: str,
    worker_id: str,
    workflow: str = "",
    app_id: str = "",
    runs_dir: Path | None = None,
    state_root: Path | None = None,
    stale_after_seconds: float = 120.0,
) -> FileLease | None:
    return FileLease.acquire(
        worker_state_dir(runs_dir, state_root=state_root) / "locks" / f"gpu-{int(gpu_id)}.json",
        run_id=run_id,
        owner_id=worker_id,
        extra={"gpu_id": int(gpu_id), "workflow": workflow, "app_id": app_id},
        stale_after_seconds=stale_after_seconds,
    )


def acquire_first_gpu_lease(
    gpu_ids: Sequence[int],
    *,
    run_id: str,
    worker_id: str,
    workflow: str = "",
    app_id: str = "",
    runs_dir: Path | None = None,
    state_root: Path | None = None,
    stale_after_seconds: float = 120.0,
) -> tuple[int, FileLease] | None:
    for gpu_id in gpu_ids:
        lease = acquire_gpu_lease(
            int(gpu_id),
            run_id=run_id,
            worker_id=worker_id,
            workflow=workflow,
            app_id=app_id,
            runs_dir=runs_dir,
            state_root=state_root,
            stale_after_seconds=stale_after_seconds,
        )
        if lease is not None:
            return int(gpu_id), lease
    return None


@dataclass
class MemoryLease:
    """A cross-application RAM and scratch-space reservation for one job."""

    path: Path
    run_id: str
    owner_id: str
    request: AdmissionRequest
    payload: dict[str, object]

    def update_allocation(self, **fields: object) -> None:
        current = _read_json(self.path)
        if str(current.get("run_id") or "") != self.run_id or str(current.get("owner_id") or "") != self.owner_id:
            raise RuntimeError(f"Memory lease ownership changed: {self.path}")
        current.update(fields)
        current["heartbeat_at"] = _utc_now_iso()
        _atomic_write(self.path, current)
        self.payload = current

    def set_owner_pid(self, pid: int) -> None:
        current = _read_json(self.path)
        if str(current.get("run_id") or "") != self.run_id or str(current.get("owner_id") or "") != self.owner_id:
            raise RuntimeError(f"Memory lease ownership changed: {self.path}")
        current["pid"] = int(pid)
        current["heartbeat_at"] = _utc_now_iso()
        _atomic_write(self.path, current)
        self.payload = current

    def heartbeat(self) -> None:
        current = _read_json(self.path)
        if str(current.get("run_id") or "") != self.run_id or str(current.get("owner_id") or "") != self.owner_id:
            raise RuntimeError(f"Memory lease ownership changed: {self.path}")
        current["heartbeat_at"] = _utc_now_iso()
        _atomic_write(self.path, current)
        self.payload = current

    def release(self) -> bool:
        current = _read_json(self.path)
        if str(current.get("run_id") or "") != self.run_id or str(current.get("owner_id") or "") != self.owner_id:
            return False
        try:
            self.path.unlink()
        except FileNotFoundError:
            return False
        return True


def _memory_reservation_dir(state_root: Path) -> Path:
    return Path(state_root) / ".worker" / "reservations"


def active_memory_reservations(
    state_root: Path | None = None,
    *,
    stale_after_seconds: float = 120.0,
) -> list[dict[str, Any]]:
    root = _memory_reservation_dir(state_root or shared_state_dir())
    if not root.is_dir():
        return []
    active: list[dict[str, Any]] = []
    for path in sorted(root.glob("*.json")):
        payload = _read_json(path)
        if _heartbeat_is_stale(payload, stale_after_seconds=stale_after_seconds):
            try:
                path.unlink()
            except OSError:
                pass
            continue
        payload["reservation_path"] = str(path)
        active.append(payload)
    return active


def active_resource_reservations(
    state_root: Path | None = None,
    *,
    stale_after_seconds: float = 120.0,
) -> list[dict[str, Any]]:
    """Return live per-job CPU, GPU, RAM, and scratch reservations.

    CPU and GPU are file leases, while RAM and scratch use reservation records.
    Merge them by app, run, and worker so either MN app can show the shared
    machine allocation, including jobs started by the other app.
    """

    root = Path(state_root or shared_state_dir()).expanduser().resolve()
    records: dict[tuple[str, str, str], dict[str, Any]] = {}

    def record_for(payload: dict[str, Any]) -> dict[str, Any]:
        app_id = str(payload.get("app_id") or "unknown")
        run_id = str(payload.get("run_id") or "")
        owner_id = str(payload.get("owner_id") or "")
        key = (app_id, run_id, owner_id)
        return records.setdefault(
            key,
            {
                "app_id": app_id,
                "run_id": run_id,
                "owner_id": owner_id,
                "run_dir": str(payload.get("run_dir") or ""),
                "scratch_path": str(payload.get("scratch_path") or ""),
                "cpu_threads": 0,
                "cpu_slots": [],
                "gpu": False,
                "gpu_devices": [],
                "min_vram_gb": 0.0,
                "ram_gb": 0.0,
                "scratch_gb": 0.0,
                "heartbeat_at": "",
            },
        )

    for payload in active_memory_reservations(root, stale_after_seconds=stale_after_seconds):
        row = record_for(payload)
        for field in ("run_dir", "scratch_path", "heartbeat_at"):
            if payload.get(field):
                row[field] = payload[field]
        row["cpu_threads"] = int(payload.get("cpu_threads") or 0)
        row["gpu"] = bool(payload.get("gpu"))
        row["min_vram_gb"] = float(payload.get("min_vram_gb") or 0.0)
        row["ram_gb"] = float(payload.get("ram_gb") or 0.0)
        row["scratch_gb"] = float(payload.get("scratch_gb") or 0.0)
        row["gpu_devices"] = [str(value) for value in payload.get("gpu_devices") or []]

    lock_dir = root / ".worker" / "locks"
    if lock_dir.is_dir():
        for pattern in ("cpu-*.json", "gpu-*.json"):
            for path in sorted(lock_dir.glob(pattern)):
                payload = _read_json(path)
                if _heartbeat_is_stale(payload, stale_after_seconds=stale_after_seconds):
                    try:
                        path.unlink()
                    except OSError:
                        pass
                    continue
                row = record_for(payload)
                if payload.get("heartbeat_at") and str(payload["heartbeat_at"]) > str(row["heartbeat_at"]):
                    row["heartbeat_at"] = payload["heartbeat_at"]
                if path.name.startswith("cpu-"):
                    try:
                        slot = int(payload.get("cpu_slot"))
                    except (TypeError, ValueError):
                        continue
                    if slot not in row["cpu_slots"]:
                        row["cpu_slots"].append(slot)
                    row["cpu_threads"] = max(
                        int(row["cpu_threads"]), int(payload.get("requested_cpu_threads") or 0)
                    )
                else:
                    try:
                        device = str(int(payload.get("gpu_id")))
                    except (TypeError, ValueError):
                        continue
                    if device not in row["gpu_devices"]:
                        row["gpu_devices"].append(device)

    result = list(records.values())
    for row in result:
        row["cpu_slots"].sort()
        row["gpu_devices"].sort(key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))
    return sorted(result, key=lambda row: (row["app_id"], row["run_id"], row["owner_id"]))


def _exclusive_file_lock(path: Path):
    """Open and exclusively lock a tiny allocator file on POSIX systems."""
    try:
        import fcntl
    except ImportError as exc:  # pragma: no cover - MN deployments are Linux
        raise RuntimeError("Shared resource admission requires POSIX file locks") from exc
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a+")
    fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
    return handle


def acquire_memory_lease(
    request: AdmissionRequest,
    snapshot: ResourceSnapshot,
    *,
    run_id: str,
    owner_id: str,
    app_id: str = "",
    run_dir: Path | None = None,
    state_root: Path | None = None,
    scratch_path: Path | None = None,
    stale_after_seconds: float = 120.0,
) -> tuple[MemoryLease | None, AdmissionDecision]:
    """Atomically check and reserve declared RAM and scratch requirements.

    The input snapshot is captured immediately before admission. Existing shared
    reservations are deducted while holding a cross-process lock, so two apps
    cannot both reserve the same remaining RAM or scratch capacity.
    """
    root = Path(state_root or shared_state_dir()).expanduser().resolve()
    requested_scratch_path = Path(scratch_path or tempfile.gettempdir()).expanduser().resolve()
    try:
        requested_scratch_device = int(requested_scratch_path.stat().st_dev)
    except OSError:
        requested_scratch_device = None
    reservation_dir = _memory_reservation_dir(root)
    reservation_dir.mkdir(parents=True, exist_ok=True)
    lock_handle = _exclusive_file_lock(reservation_dir / "allocator.lock")
    try:
        reservations = active_memory_reservations(root, stale_after_seconds=stale_after_seconds)
        held_ram = sum(float(row.get("ram_gb") or 0.0) for row in reservations)
        held_scratch = 0.0
        for row in reservations:
            try:
                reserved_device = int(row.get("scratch_device"))
            except (TypeError, ValueError):
                reserved_device = None
            try:
                reserved_path = Path(str(row.get("scratch_path") or "")).expanduser().resolve()
            except (OSError, RuntimeError):
                reserved_path = None
            same_filesystem = (
                requested_scratch_device is not None
                and reserved_device is not None
                and requested_scratch_device == reserved_device
            ) or (
                requested_scratch_device is None
                and reserved_device is None
                and reserved_path == requested_scratch_path
            )
            if same_filesystem:
                held_scratch += float(row.get("scratch_gb") or 0.0)
        adjusted = ResourceSnapshot(
            cpu_threads_total=snapshot.cpu_threads_total,
            ram_available_gb=max(0.0, snapshot.ram_available_gb - held_ram),
            ram_total_gb=snapshot.ram_total_gb,
            scratch_free_gb=max(0.0, snapshot.scratch_free_gb - held_scratch),
            scratch_total_gb=snapshot.scratch_total_gb,
            gpus=snapshot.gpus,
        )
        memory_request = AdmissionRequest(
            gpu=False,
            cpu_threads=0,
            ram_gb=request.ram_gb,
            scratch_gb=request.scratch_gb,
        )
        decision = assess_resource_admission(memory_request, adjusted)
        if not decision.allowed:
            return None, decision
        safe_run_id = "".join(character if character.isalnum() or character in "._-" else "_" for character in run_id)
        path = reservation_dir / f"{safe_run_id}-{uuid4().hex[:12]}.json"
        now = _utc_now_iso()
        payload: dict[str, object] = {
            "schema_version": 1,
            "run_id": run_id,
            "owner_id": owner_id,
            "pid": os.getpid(),
            "app_id": app_id or os.getenv("MN_COMPUTE_SCHEDULER_APP_ID", ""),
            "run_dir": str(Path(run_dir).expanduser().resolve()) if run_dir is not None else "",
            "acquired_at": now,
            "heartbeat_at": now,
            "cpu_threads": request.cpu_threads,
            "gpu": request.gpu,
            "min_vram_gb": request.min_vram_gb,
            "ram_gb": request.ram_gb,
            "scratch_gb": request.scratch_gb,
            "scratch_path": str(requested_scratch_path),
            "scratch_device": requested_scratch_device,
        }
        _atomic_write(path, payload)
        return MemoryLease(path, run_id, owner_id, request, payload), decision
    finally:
        lock_handle.close()
