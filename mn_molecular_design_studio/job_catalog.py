"""Molecular Design Studio's searchable SQLite index for filesystem jobs.

Run directories remain the durable source of job inputs and outputs. This
catalog belongs to the Studio workspace and stores searchable metadata plus
relative artifact links beneath the configured Results directory.
"""

from __future__ import annotations

import hashlib
import json
import logging
import mmap
import os
import re
import sqlite3
import threading
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from mn_molecular_design_studio.engine.workflow_engine.core.job_identity import job_display_tool, job_name, job_tool

_SCHEMA_LOCK = threading.Lock()
_RECONCILE_LOCK = threading.Lock()
_RECONCILE_SCHEDULE_LOCK = threading.Lock()
_RECONCILE_SCHEDULED: set[str] = set()
_SCHEMA_READY_PATH: Path | None = None
# Reconciliation stats every run directory and refreshes artifacts only when
# inputs/results change. A short interval made large Results trees contend with
# every synchronous API request, so keep it bounded and let concurrent callers
# use the last indexed snapshot while one scan is in progress.
_ROOT_SCAN_INTERVAL_SECONDS = 10.0
_RESULT_SUCCESS_UNSET = object()
_RESULT_INLINE_READ_LIMIT = 8 * 1024 * 1024
_RESULT_TAIL_READ_BYTES = 256 * 1024
_LOG = logging.getLogger(__name__)


def catalog_path() -> Path:
    from mn_molecular_design_studio.workspace import _studio_home

    return (_studio_home() / "database" / "job_catalog.sqlite3").resolve()


def _results_root() -> Path:
    from mn_molecular_design_studio.runtime import runs_root

    return runs_root().expanduser().resolve()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _connect() -> sqlite3.Connection:
    global _SCHEMA_READY_PATH
    path = catalog_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path, timeout=20)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA busy_timeout = 20000")
    connection.execute("PRAGMA foreign_keys = ON")
    try:
        if _SCHEMA_READY_PATH != path:
            with _SCHEMA_LOCK:
                if _SCHEMA_READY_PATH != path:
                    _initialize(connection)
                    _SCHEMA_READY_PATH = path
    except Exception:
        connection.close()
        raise
    return connection


@contextmanager
def _connection():
    connection = _connect()
    try:
        with connection:
            yield connection
    finally:
        connection.close()


def _initialize(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA journal_mode = WAL")
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS storage_roots (
            id TEXT PRIMARY KEY,
            path TEXT NOT NULL UNIQUE,
            created_at TEXT NOT NULL DEFAULT '',
            last_scanned_at TEXT
        );
        CREATE TABLE IF NOT EXISTS jobs (
            id INTEGER PRIMARY KEY,
            storage_root_id TEXT NOT NULL REFERENCES storage_roots(id) ON DELETE CASCADE,
            task_group TEXT NOT NULL,
            run_id TEXT NOT NULL,
            run_relative_path TEXT NOT NULL,
            job_code TEXT NOT NULL DEFAULT '',
            job_type TEXT NOT NULL DEFAULT '',
            tool TEXT NOT NULL DEFAULT '',
            status TEXT NOT NULL DEFAULT 'unknown',
            queue_resource TEXT NOT NULL DEFAULT '',
            resource_allocation TEXT NOT NULL DEFAULT '',
            scheduler_wait_reason TEXT NOT NULL DEFAULT '',
            current_phase TEXT NOT NULL DEFAULT '',
            current_engine TEXT NOT NULL DEFAULT '',
            warning TEXT NOT NULL DEFAULT '',
            warning_level TEXT NOT NULL DEFAULT '',
            campaign_name TEXT NOT NULL DEFAULT '',
            campaign_id TEXT NOT NULL DEFAULT '',
            pipeline_name TEXT NOT NULL DEFAULT '',
            pipeline_parent_node_id TEXT NOT NULL DEFAULT '',
            pipeline_iteration INTEGER NOT NULL DEFAULT 0,
            campaign_step TEXT NOT NULL DEFAULT '',
            campaign_step_index TEXT NOT NULL DEFAULT '',
            created_at TEXT NOT NULL DEFAULT '',
            updated_at TEXT NOT NULL DEFAULT '',
            queue_wait_started_at TEXT NOT NULL DEFAULT '',
            worker_started_at TEXT NOT NULL DEFAULT '',
            completed_at TEXT NOT NULL DEFAULT '',
            success INTEGER,
            hidden INTEGER NOT NULL DEFAULT 0,
            capacity_parent INTEGER NOT NULL DEFAULT 0,
            owner_id TEXT,
            project_id TEXT,
            metadata_mtime_ns INTEGER NOT NULL DEFAULT 0,
            input_mtime_ns INTEGER NOT NULL DEFAULT 0,
            result_mtime_ns INTEGER NOT NULL DEFAULT 0,
            indexed_at TEXT NOT NULL,
            UNIQUE(storage_root_id, task_group, run_id)
        );
        CREATE INDEX IF NOT EXISTS studio_jobs_created
            ON jobs(storage_root_id, created_at DESC, run_id DESC);
        CREATE INDEX IF NOT EXISTS studio_jobs_status_created
            ON jobs(storage_root_id, status, created_at DESC);
        CREATE INDEX IF NOT EXISTS studio_jobs_tool_created
            ON jobs(storage_root_id, tool, created_at DESC);
        CREATE INDEX IF NOT EXISTS studio_jobs_project_created
            ON jobs(storage_root_id, project_id, created_at DESC);
        CREATE INDEX IF NOT EXISTS studio_jobs_owner_created
            ON jobs(storage_root_id, owner_id, created_at DESC);
        CREATE TABLE IF NOT EXISTS job_artifacts (
            id INTEGER PRIMARY KEY,
            job_id INTEGER NOT NULL REFERENCES jobs(id) ON DELETE CASCADE,
            relative_path TEXT NOT NULL,
            name TEXT NOT NULL,
            kind TEXT NOT NULL DEFAULT 'file',
            size_bytes INTEGER NOT NULL DEFAULT 0,
            modified_at TEXT NOT NULL DEFAULT '',
            UNIQUE(job_id, relative_path)
        );
        CREATE INDEX IF NOT EXISTS studio_artifacts_by_job ON job_artifacts(job_id, name);
        """
    )
    existing_columns = {row[1] for row in connection.execute("PRAGMA table_info(jobs)")}
    added_columns = False
    if "run_identity_json" not in existing_columns:
        connection.execute("ALTER TABLE jobs ADD COLUMN run_identity_json TEXT NOT NULL DEFAULT '{}'")
        added_columns = True
    if "pipeline_name" not in existing_columns:
        connection.execute("ALTER TABLE jobs ADD COLUMN pipeline_name TEXT NOT NULL DEFAULT ''")
        added_columns = True
    for column in ("queue_wait_started_at", "worker_started_at", "completed_at", "pipeline_parent_node_id"):
        if column not in existing_columns:
            connection.execute(f"ALTER TABLE jobs ADD COLUMN {column} TEXT NOT NULL DEFAULT ''")
            added_columns = True
    if "pipeline_iteration" not in existing_columns:
        connection.execute("ALTER TABLE jobs ADD COLUMN pipeline_iteration INTEGER NOT NULL DEFAULT 0")
        added_columns = True
    if added_columns:
        # Existing records must be re-read once so newly indexed metadata can
        # be derived from durable run metadata without requiring job changes.
        connection.execute("UPDATE jobs SET metadata_mtime_ns=0")
    from mn_molecular_design_studio.logical_runs import IDENTITY_FIELDS
    connection.execute("CREATE TABLE IF NOT EXISTS catalog_projection (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    revision = json.dumps([*IDENTITY_FIELDS, "database_analysis:v1"])
    previous = connection.execute("SELECT value FROM catalog_projection WHERE key='run_identity'").fetchone()
    if previous is None or previous[0] != revision:
        connection.execute("UPDATE jobs SET metadata_mtime_ns=0")
        connection.execute("INSERT OR REPLACE INTO catalog_projection VALUES ('run_identity', ?)", (revision,))
    connection.execute("CREATE INDEX IF NOT EXISTS studio_jobs_analysis_source ON jobs "
                       "(storage_root_id,json_extract(run_identity_json,'$.database_analysis.source_key'),tool)")
    try:
        connection.executescript(
            """
            CREATE VIRTUAL TABLE IF NOT EXISTS jobs_fts USING fts5(
                campaign_name, tool, job_type, job_code, run_id, task_group, status,
                content='jobs', content_rowid='id'
            );
            CREATE TRIGGER IF NOT EXISTS studio_jobs_fts_insert AFTER INSERT ON jobs BEGIN
                INSERT INTO jobs_fts(rowid, campaign_name, tool, job_type, job_code, run_id, task_group, status)
                VALUES (new.id, new.campaign_name, new.tool, new.job_type, new.job_code, new.run_id, new.task_group, new.status);
            END;
            CREATE TRIGGER IF NOT EXISTS studio_jobs_fts_delete AFTER DELETE ON jobs BEGIN
                INSERT INTO jobs_fts(jobs_fts, rowid, campaign_name, tool, job_type, job_code, run_id, task_group, status)
                VALUES ('delete', old.id, old.campaign_name, old.tool, old.job_type, old.job_code, old.run_id, old.task_group, old.status);
            END;
            CREATE TRIGGER IF NOT EXISTS studio_jobs_fts_update AFTER UPDATE ON jobs BEGIN
                INSERT INTO jobs_fts(jobs_fts, rowid, campaign_name, tool, job_type, job_code, run_id, task_group, status)
                VALUES ('delete', old.id, old.campaign_name, old.tool, old.job_type, old.job_code, old.run_id, old.task_group, old.status);
                INSERT INTO jobs_fts(rowid, campaign_name, tool, job_type, job_code, run_id, task_group, status)
                VALUES (new.id, new.campaign_name, new.tool, new.job_type, new.job_code, new.run_id, new.task_group, new.status);
            END;
            """
        )
        fts_count = connection.execute("SELECT count(*) FROM jobs_fts").fetchone()[0]
        jobs_count = connection.execute("SELECT count(*) FROM jobs").fetchone()[0]
        if jobs_count and not fts_count:
            connection.execute("INSERT INTO jobs_fts(jobs_fts) VALUES ('rebuild')")
    except sqlite3.OperationalError:
        # SQLite builds without FTS5 still use indexed filters and LIKE search.
        pass


def _root_id(root: Path) -> str:
    return hashlib.sha256(str(root.expanduser().resolve()).encode("utf-8")).hexdigest()[:24]


def _mtime(path: Path) -> int:
    try:
        return path.stat().st_mtime_ns
    except OSError:
        return 0


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _read_result_index_fields(path: Path) -> tuple[Any, list[Any]]:
    """Read only the compact fields needed by the job catalog.

    Scientific result files can include every candidate structure and exceed
    100 MB. Loading those files just to index ``success`` made the API retain a
    large allocator high-water mark during Results reconciliation. Keep the
    normal full read for small library jobs (which provide the worker detail
    rows), and inspect only the result tail for the root-level success field in
    large, pretty-printed results.
    """
    try:
        size = path.stat().st_size
    except OSError:
        return None, []
    if size <= _RESULT_INLINE_READ_LIMIT:
        result = _read_json(path)
        metrics = result.get("metrics") if isinstance(result.get("metrics"), dict) else {}
        workers = metrics.get("workers") if isinstance(metrics.get("workers"), list) else []
        return result.get("success"), workers

    try:
        with path.open("rb") as handle:
            if size > _RESULT_TAIL_READ_BYTES:
                handle.seek(-_RESULT_TAIL_READ_BYTES, os.SEEK_END)
            tail = handle.read(_RESULT_TAIL_READ_BYTES)
    except OSError:
        return None, []

    # Studio writes indented JSON with two spaces per object depth. A root key
    # is therefore at exactly two spaces, unlike nested success fields inside
    # candidate records. Success is emitted after the potentially large output
    # payload, so this bounded tail read normally finds it without touching the
    # candidate data. mmap provides a low-allocation fallback for older files
    # whose success key is not in the final 256 KiB.
    success_pattern = re.compile(rb'(?m)^  "success"\s*:\s*(true|false|null)\b')
    match = success_pattern.search(tail)
    if match is None:
        try:
            with path.open("rb") as handle:
                with mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
                    match = success_pattern.search(mapped)
        except (OSError, ValueError):
            match = None
    if match is None:
        return None, []
    token = match.group(1)
    success = True if token == b"true" else False if token == b"false" else None
    return success, []


def analysis_source_key(source):
    fields = ("molecule_id", "structure_id", "group_id", "structure_digest", "model_id")
    values = [str(source.get(key)) if source.get(key) is not None else "" for key in fields]
    if not all(values):
        raise ValueError("Exact Database analysis source identity is incomplete.")
    return hashlib.sha256(json.dumps(values).encode()).hexdigest()


def _database_analysis_projection(run_dir, metadata, inputs):
    """Rebuildable compact execution metadata, authored during normal indexing.

    Only the exact Database subject is indexed, never a complete population or
    coordinate payload. Reading a manifest here is independent of UI rendering.
    """
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    params = inputs.get("params") if isinstance(inputs.get("params"), dict) else {}
    source = metadata.get("database_source") if isinstance(metadata.get("database_source"), dict) else {}
    if not source:
        source_path = (inputs.get("inputs") or {}).get("source_run_dir") or params.get("source_run_dir")
        if source_path:
            try:
                resolved = resolve_stored_path(source_path, run_dir=run_dir, must_exist=True)
            except (OSError, TypeError, ValueError):
                resolved = None
            if resolved:
                source = _read_json(Path(resolved) / "complex_analysis_source.json")
    if not source or not source.get("structure_id"):
        return None
    try:
        key = analysis_source_key(source)
    except ValueError:
        return {"source_key": "", "invalid": True}
    batch = _read_json(run_dir / "artifacts/analysis/batch_manifest.json")
    selection = _read_json(run_dir / "artifacts/analysis/complex_analysis_selection.json")
    selected = batch.get("selection") or selection or params
    candidate = str(source.get("candidate_id") or next(iter(params.get("candidate_ids") or []), ""))
    item = next((row for row in batch.get("items") or [] if str(row.get("candidate_id") or "") == candidate), {})
    binder = selected.get("binder_chains") or []
    target = selected.get("target_chains") or []
    valid = not batch.get("items") or bool(item)
    if item:
        valid = valid and str(item.get("source_structure_digest") or item.get("structure_digest") or "") == source["structure_digest"]
        if str(metadata.get("tool") or inputs.get("tool")) == "biotite_interface":
            valid = valid and sorted(item.get("binder_chains", binder)) == sorted(binder)
            valid = valid and sorted(item.get("target_chains", target)) == sorted(target)
        if str(metadata.get("tool") or inputs.get("tool")) in {"plip", "pandamap", "posebusters"}:
            valid = valid and item.get("ligand_selection") == selected.get("ligand_selection")
    outcome = next((row for row in _read_json(run_dir / "artifacts/analysis/batch_status.json").get("items") or []
                    if row.get("id") == item.get("id")), {}) if item else {}
    return {"source_key": key, "source": {k: v for k, v in source.items() if k in {
                "molecule_id", "structure_id", "group_id", "structure_digest", "model_id", "pae_digest", "ligand_selection",
                "reference_context_id", "reference_molecule_id", "reference_structure_id", "reference_structure_digest",
                "reference_model_id", "reference_chain_id", "placement_run_id"}},
            "candidate_id": candidate, "model_id": str(selected.get("model_id") or selected.get("selected_model_id") or source["model_id"]),
            "binder_chains": binder, "target_chains": target, "ligand_selection": selected.get("ligand_selection") or source.get("ligand_selection") or {},
            "parameters": batch.get("settings") or selection.get("parameters") or params,
            "protocol": params.get("assessment_protocol") or {k: v for k, v in (batch.get("settings") or params).items()
                if k in {"method_version", "protocol_revision", "engine_version", "image", "image_digest"}},
            "item_status": outcome.get("status"), "error": outcome.get("error") or metadata.get("error") or "",
            "invalid": not valid}


def database_analysis_jobs(user, sources):
    """Read exact source-qualified indexed jobs; never reconcile or read Results.

    SQL pages are bounded to 500 rows and filtered by source keys before paging.
    ACLs are checked afresh, without a global first-N historical job limit.
    """
    from mn_molecular_design_studio.auth import accessible_project_ids
    keys = sorted({analysis_source_key(source) for source in sources})
    if not keys:
        return []
    root = _results_root()
    parameters = [_root_id(root)]
    access = "j.storage_root_id=? AND j.task_group='analysis'"
    if user.get("role") != "admin":
        access += " AND (j.owner_id=? OR j.project_id IN (SELECT value FROM json_each(?)))"
        parameters += [str(user.get("id") or ""), json.dumps(accessible_project_ids(user))]
    with _connection() as connection:
        if not connection.execute("SELECT 1 FROM storage_roots WHERE id=?", (parameters[0],)).fetchone():
            raise RuntimeError("Analysis job index is not initialized; refresh through the existing Results catalog before submitting.")
        pending = connection.execute("SELECT 1 FROM jobs j WHERE " + access +
            " AND (j.metadata_mtime_ns=0 OR json_type(j.run_identity_json,'$.database_analysis') IS NULL) LIMIT 1", parameters).fetchone()
        if pending:
            raise RuntimeError("Analysis job index upgrade is incomplete; wait for the existing catalog reconciliation.")
        rows, offset = [], 0
        while True:
            page = connection.execute("SELECT j.* FROM jobs j WHERE " + access +
                " AND json_extract(j.run_identity_json,'$.database_analysis.source_key') IN (SELECT value FROM json_each(?)) "
                "ORDER BY j.created_at DESC,j.id DESC LIMIT 500 OFFSET ?", [*parameters, json.dumps(keys), offset]).fetchall()
            rows.extend(page)
            if len(page) < 500:
                break
            offset += 500
    result = []
    for row in rows:
        item = dict(row)
        item["database_analysis"] = json.loads(item.pop("run_identity_json"))["database_analysis"]
        item["run_dir"] = root / item["run_relative_path"]
        result.append(item)
    return result


def _job_payload(
    run_dir: Path,
    *,
    cached_result_success: object = _RESULT_SUCCESS_UNSET,
) -> dict[str, Any]:
    metadata = _read_json(run_dir / "metadata.json")
    inputs = _read_json(run_dir / "input.json")
    worker_metrics: list[Any] = []
    if cached_result_success is _RESULT_SUCCESS_UNSET:
        result_success, worker_metrics = _read_result_index_fields(run_dir / "result.json")
    else:
        result_success = cached_result_success
    params = inputs.get("params") if isinstance(inputs.get("params"), dict) else {}
    status = str(metadata.get("status") or "").casefold()
    if not status:
        status = "completed" if result_success else "unknown"
    if result_success is False and status not in {
        "running", "queued", "preparing", "pending", "starting", "cancelling", "cancel_requested",
        "paused", "holding", "stopped", "cancelled",
    }:
        status = "failed"
    allocation = metadata.get("resource_allocation") if isinstance(metadata.get("resource_allocation"), dict) else {}
    resource_allocation = ""
    if allocation:
        try:
            cpu = int(allocation.get("cpu_cores") or 0)
        except (TypeError, ValueError):
            cpu = 0
        gpu = str(allocation.get("gpu_device") or "cpu")
        resource_allocation = (
            f"{cpu} CPU core{'s' if cpu != 1 else ''} · {gpu}" if cpu else f"coordination · {gpu}"
        )
    marker = _read_json(run_dir / "molecular_design_studio.json")
    pipeline_parent = str(metadata.get("pipeline_parent_run_dir") or metadata.get("active_learning_parent_run") or "").strip()
    is_visual_pipeline = (
        str(metadata.get("tool") or inputs.get("tool") or "") == "visual_pipeline"
        or str(metadata.get("job_type") or inputs.get("job_type") or "") == "visual_pipeline"
    )
    pipeline_name = (
        str(metadata.get("pipeline_name") or params.get("pipeline_name") or metadata.get("campaign_name") or "").strip()
        if pipeline_parent or is_visual_pipeline
        else ""
    )
    # Keep the worker rows used by the distributed-library dialog available
    # in the compact job catalog. Avoid copying arbitrary result metrics into
    # every catalog response; some workflows retain large result summaries.
    from mn_molecular_design_studio.logical_runs import IDENTITY_FIELDS
    identity = {key: metadata[key] for key in IDENTITY_FIELDS if key in metadata}
    for key in ("campaign_group_id", "campaign_group_name"):
        if not identity.get(key) and params.get(key):
            identity[key] = params[key]
    collection = run_dir / "artifacts/normalized_candidates/candidates.jsonl"
    identity["has_candidate_collection"] = collection.is_file() and collection.stat().st_size > 0
    identity["has_docking_results"] = (run_dir / "artifacts/docking_scores.csv").is_file()
    if run_dir.parent.name == "analysis":
        identity["database_analysis"] = _database_analysis_projection(run_dir, metadata, inputs)
    return {
        "run_identity_json": json.dumps(identity, sort_keys=True),
        "run_id": run_dir.name,
        "task_group": str(metadata.get("task_group") or run_dir.parent.name),
        "job_code": str(metadata.get("job_code") or run_dir.name),
        "job_type": str(metadata.get("job_type") or inputs.get("job_type") or ""),
        "tool": job_tool(metadata.get("tool") or inputs.get("tool"), params),
        "status": status,
        "queue_resource": str(metadata.get("queue_resource") or ""),
        "resource_allocation": resource_allocation,
        "scheduler_wait_reason": str(metadata.get("scheduler_wait_reason") or ""),
        "current_phase": str(metadata.get("current_phase") or ""),
        "worker_plan": metadata.get("worker_plan") if isinstance(metadata.get("worker_plan"), list) else [],
        "metrics": {"workers": worker_metrics} if worker_metrics else {},
        "current_engine": str(metadata.get("current_engine") or ""),
        "warning": str(metadata.get("warning") or ""),
        "warning_level": str(metadata.get("warning_level") or ("degraded" if metadata.get("warning") else "")),
        "campaign_name": job_name(metadata, params),
        "campaign_id": str(metadata.get("campaign_id") or (Path(pipeline_parent.rstrip("/")).name if pipeline_parent else "")),
        "pipeline_parent_node_id": str(metadata.get("pipeline_parent_node_id") or ""),
        "pipeline_iteration": int(metadata.get("pipeline_iteration") or metadata.get("active_learning_round") or 0),
        "pipeline_name": pipeline_name,
        "campaign_step": str(metadata.get("campaign_step") or ""),
        "campaign_step_index": str(metadata.get("campaign_step_index") or ""),
        "created_at": str(metadata.get("created_at") or ""),
        "updated_at": str(metadata.get("updated_at") or ""),
        "queue_wait_started_at": str(metadata.get("queue_wait_started_at") or ""),
        "worker_started_at": str(metadata.get("worker_started_at") or ""),
        "completed_at": str(metadata.get("completed_at") or ""),
        "success": result_success,
        "hidden": bool(metadata.get("hidden")),
        "capacity_parent": bool(metadata.get("capacity_parent_run_id") or metadata.get("capacity_parent")),
        # Only Studio identity metadata is used here. Scheduler/worker owner IDs
        # are operational lease IDs, not user identities.
        "owner_id": metadata.get("studio_owner_id") or marker.get("studio_owner_id"),
        "project_id": metadata.get("project_id") or marker.get("project_id"),
        "metadata_mtime_ns": _mtime(run_dir / "metadata.json"),
        "input_mtime_ns": _mtime(run_dir / "input.json"),
        "result_mtime_ns": _mtime(run_dir / "result.json"),
        "run_dir": str(run_dir.resolve()),
    }


def _artifact_rows(root: Path, run_dir: Path) -> Iterator[dict[str, Any]]:
    for base, dirs, files in os.walk(run_dir, followlinks=False):
        dirs[:] = [name for name in dirs if not (Path(base) / name).is_symlink()]
        for name in files:
            path = Path(base) / name
            try:
                resolved = path.resolve(strict=True)
                relative = resolved.relative_to(root).as_posix()
                stat = resolved.stat()
            except (OSError, ValueError):
                continue
            yield {
                "relative_path": relative,
                "name": resolved.name,
                "kind": resolved.suffix.lower().lstrip(".") or "file",
                "size_bytes": stat.st_size,
                "modified_at": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(timespec="seconds"),
            }


def _upsert(root: Path, run_dir: Path, *, refresh_artifacts: bool) -> None:
    root = root.expanduser().resolve()
    run_dir = run_dir.expanduser().resolve()
    relative = run_dir.relative_to(root).as_posix()
    if len(Path(relative).parts) < 2:
        return
    root_id = _root_id(root)
    task_group, run_id = run_dir.parent.name, run_dir.name
    result_mtime_ns = _mtime(run_dir / "result.json")
    cached_result_success: object = _RESULT_SUCCESS_UNSET
    with _connection() as connection:
        previous = connection.execute(
            "SELECT success, result_mtime_ns FROM jobs "
            "WHERE storage_root_id=? AND task_group=? AND run_id=?",
            (root_id, task_group, run_id),
        ).fetchone()
    if previous and int(previous["result_mtime_ns"]) == result_mtime_ns:
        cached_result_success = (
            None if previous["success"] is None else bool(previous["success"])
        )
    # Result files can contain very large scientific outputs. Metadata changes
    # frequently while a run is active; don't parse an unchanged result.json on
    # every such update. Its success value is already indexed in SQLite.
    row = _job_payload(run_dir, cached_result_success=cached_result_success)
    values = (
        root_id, row["task_group"], row["run_id"], relative, row["job_code"], row["job_type"], row["tool"],
        row["status"], row["queue_resource"], row["resource_allocation"], row["scheduler_wait_reason"],
        row["current_phase"], row["current_engine"], row["warning"], row["warning_level"], row["campaign_name"],
        row["campaign_id"], row["pipeline_name"], row["pipeline_parent_node_id"], row["pipeline_iteration"], row["campaign_step"], row["campaign_step_index"], row["created_at"], row["updated_at"],
        row["queue_wait_started_at"], row["worker_started_at"], row["completed_at"],
        int(row["success"]) if isinstance(row["success"], bool) else None, int(row["hidden"]), int(row["capacity_parent"]),
        row["owner_id"], row["project_id"], row["metadata_mtime_ns"], row["input_mtime_ns"], row["result_mtime_ns"], _now(),
    )
    columns = (
        "run_relative_path", "job_code", "job_type", "tool", "status", "queue_resource", "resource_allocation",
        "scheduler_wait_reason", "current_phase", "current_engine", "warning", "warning_level", "campaign_name",
        "campaign_id", "pipeline_name", "pipeline_parent_node_id", "pipeline_iteration", "campaign_step", "campaign_step_index", "created_at", "updated_at",
        "queue_wait_started_at", "worker_started_at", "completed_at", "success", "hidden",
        "capacity_parent", "owner_id", "project_id", "metadata_mtime_ns", "input_mtime_ns", "result_mtime_ns", "indexed_at",
    )
    with _connection() as connection:
        connection.execute(
            "INSERT INTO storage_roots(id, path, created_at) VALUES (?, ?, ?) ON CONFLICT(path) DO NOTHING",
            (root_id, str(root), _now()),
        )
        connection.execute(
            "INSERT INTO jobs (storage_root_id, task_group, run_id, run_relative_path, job_code, job_type, tool, status, "
            "queue_resource, resource_allocation, scheduler_wait_reason, current_phase, current_engine, warning, warning_level, "
            "campaign_name, campaign_id, pipeline_name, pipeline_parent_node_id, pipeline_iteration, campaign_step, campaign_step_index, created_at, updated_at, "
            "queue_wait_started_at, worker_started_at, completed_at, success, hidden, "
            "capacity_parent, owner_id, project_id, metadata_mtime_ns, input_mtime_ns, result_mtime_ns, indexed_at) "
            "VALUES (" + ",".join("?" for _ in values) + ") "
            "ON CONFLICT(storage_root_id, task_group, run_id) DO UPDATE SET "
            + ",".join(f"{column}=excluded.{column}" for column in columns),
            values,
        )
        job_id = connection.execute(
            "SELECT id FROM jobs WHERE storage_root_id=? AND task_group=? AND run_id=?",
            (root_id, row["task_group"], row["run_id"]),
        ).fetchone()[0]
        connection.execute("UPDATE jobs SET run_identity_json=? WHERE id=?", (row["run_identity_json"], job_id))
        if refresh_artifacts:
            connection.execute("DELETE FROM job_artifacts WHERE job_id=?", (job_id,))
            connection.executemany(
                "INSERT OR REPLACE INTO job_artifacts(job_id, relative_path, name, kind, size_bytes, modified_at) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (
                    (job_id, a["relative_path"], a["name"], a["kind"], a["size_bytes"], a["modified_at"])
                    for a in _artifact_rows(root, run_dir)
                ),
            )


def _reconcile_locked(root: Path) -> None:
    root = (root or _results_root()).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    root_id = _root_id(root)
    with _connection() as connection:
        connection.execute(
            "INSERT INTO storage_roots(id, path, created_at) VALUES (?, ?, ?) ON CONFLICT(path) DO NOTHING",
            (root_id, str(root), _now()),
        )
        scanned_row = connection.execute(
            "SELECT last_scanned_at FROM storage_roots WHERE id=?", (root_id,)
        ).fetchone()
        if scanned_row and scanned_row["last_scanned_at"]:
            try:
                scanned_at = datetime.fromisoformat(scanned_row["last_scanned_at"])
                if (datetime.now(timezone.utc) - scanned_at).total_seconds() < _ROOT_SCAN_INTERVAL_SECONDS:
                    return
            except (TypeError, ValueError):
                pass
        previous = {
            (str(row["task_group"]), str(row["run_id"])): (
                int(row["metadata_mtime_ns"]), int(row["input_mtime_ns"]), int(row["result_mtime_ns"])
            )
            for row in connection.execute(
                "SELECT task_group, run_id, metadata_mtime_ns, input_mtime_ns, result_mtime_ns "
                "FROM jobs WHERE storage_root_id=?", (root_id,)
            )
        }
    found: set[tuple[str, str]] = set()
    for group_dir in root.iterdir():
        if not group_dir.is_dir() or group_dir.is_symlink() or group_dir.name.startswith("_"):
            continue
        for run_dir in group_dir.iterdir():
            if not run_dir.is_dir() or run_dir.is_symlink() or not (run_dir / "metadata.json").is_file():
                continue
            key = (group_dir.name, run_dir.name)
            found.add(key)
            signature = (
                _mtime(run_dir / "metadata.json"),
                _mtime(run_dir / "input.json"),
                _mtime(run_dir / "result.json"),
            )
            old = previous.get(key)
            if old != signature:
                _upsert(root, run_dir, refresh_artifacts=(old is None or old[1:] != signature[1:]))
    with _connection() as connection:
        for group, run_id in previous.keys() - found:
            connection.execute(
                "DELETE FROM jobs WHERE storage_root_id=? AND task_group=? AND run_id=?",
                (root_id, group, run_id),
            )
        connection.execute("UPDATE storage_roots SET last_scanned_at=? WHERE id=?", (_now(), root_id))


def _reconcile_in_background(root: Path, root_key: str) -> None:
    try:
        if not _RECONCILE_LOCK.acquire(blocking=False):
            return
        try:
            _reconcile_locked(root)
        finally:
            _RECONCILE_LOCK.release()
    except Exception:
        _LOG.exception("Results catalog reconciliation failed for %s", root)
    finally:
        with _RECONCILE_SCHEDULE_LOCK:
            _RECONCILE_SCHEDULED.discard(root_key)


def reconcile(root: Path | None = None, *, wait: bool = True) -> None:
    selected_root = (root or _results_root()).expanduser().resolve()
    if not wait:
        root_key = str(selected_root)
        with _RECONCILE_SCHEDULE_LOCK:
            if root_key in _RECONCILE_SCHEDULED:
                return
            _RECONCILE_SCHEDULED.add(root_key)
        thread = threading.Thread(
            target=_reconcile_in_background,
            args=(selected_root, root_key),
            name="studio-results-catalog",
            daemon=True,
        )
        try:
            thread.start()
        except RuntimeError:
            with _RECONCILE_SCHEDULE_LOCK:
                _RECONCILE_SCHEDULED.discard(root_key)
            _LOG.exception("Could not start Results catalog reconciliation for %s", selected_root)
        return
    # search_jobs is called by both the UI and cluster heartbeats. Waiting on a
    # long artifact reconciliation here can occupy every FastAPI sync worker.
    # Concurrent callers can serve the current catalog while an explicit
    # maintenance operation completes its requested synchronous scan.
    if not _RECONCILE_LOCK.acquire(blocking=False):
        return
    try:
        _reconcile_locked(selected_root)
    finally:
        _RECONCILE_LOCK.release()


def _reconcile_for_read(root: Path) -> None:
    """Keep read paths responsive when a populated catalog needs refreshing."""
    with _connection() as connection:
        indexed = connection.execute(
            "SELECT 1 FROM jobs WHERE storage_root_id=? LIMIT 1", (_root_id(root),)
        ).fetchone()
    # A first index must be visible to the caller. Once the catalog has a
    # usable snapshot, refresh it off-thread so a large Results tree cannot
    # block the UI/API request that triggered the refresh.
    reconcile(root, wait=indexed is None)


def update_job_project(run_dir: str | Path) -> None:
    """Refresh a job's indexed project assignment after its marker is written."""
    selected_root = _results_root()
    selected_run = Path(run_dir).expanduser().resolve()
    try:
        selected_run.relative_to(selected_root)
    except ValueError:
        return
    if not (selected_run / "metadata.json").is_file():
        return
    _upsert(selected_root, selected_run, refresh_artifacts=False)


def _search_expression(query: str) -> str | None:
    tokens = re.findall(r"[\w]+", query, flags=re.UNICODE)
    return " AND ".join(f'"{token}"*' for token in tokens) or None


def set_job_hidden(task_group: str, run_id: str, hidden: bool = True) -> None:
    """Hide or restore a run in Results without deleting its files."""
    if not isinstance(task_group, str) or not isinstance(run_id, str):
        raise ValueError("Invalid job reference.")
    if not task_group or not run_id or any(part in {".", ".."} or "/" in part or "\\" in part for part in (task_group, run_id)):
        raise ValueError("Invalid job reference.")
    root = _results_root()
    run_dir = (root / task_group / run_id).resolve()
    try:
        run_dir.relative_to(root)
    except ValueError as exc:
        raise ValueError("Invalid job reference.") from exc
    metadata_path = run_dir / "metadata.json"
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FileNotFoundError("Job not found.") from exc
    if not isinstance(metadata, dict):
        raise ValueError("Invalid job metadata.")
    metadata["hidden"] = bool(hidden)
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json
    write_json(metadata_path, metadata)
    with _connection() as connection:
        connection.execute(
            "UPDATE jobs SET hidden=?, metadata_mtime_ns=? WHERE storage_root_id=? AND task_group=? AND run_id=?",
            (int(bool(hidden)), _mtime(metadata_path), _root_id(root), task_group, run_id),
        )


def search_jobs(
    *,
    query: str = "",
    statuses: list[str] | None = None,
    tools: list[str] | None = None,
    task_group: str | None = None,
    owner_id: str | None = None,
    project_id: str | None = None,
    access_owner_id: str | None = None,
    accessible_project_ids: list[str] | None = None,
    access_admin: bool = False,
    limit: int = 100,
    offset: int = 0,
    sort_by: str = "created_desc",
    cursor: str | None = None,
    include_hidden: bool = False,
) -> dict[str, Any]:
    if not isinstance(query, str):
        raise ValueError("query must be a string")
    if type(limit) is not int or not 1 <= limit <= 1000:
        raise ValueError("limit must be between 1 and 1000")
    if type(offset) is not int or offset < 0:
        raise ValueError("offset must be a non-negative integer")
    order_by = {
        "created_desc": "j.created_at DESC, j.run_id DESC",
        "active_desc": (
            "CASE WHEN lower(j.status) IN ('completed','succeeded','failed','cancelled','canceled') "
            "THEN 1 ELSE 0 END ASC, "
            "COALESCE(NULLIF(j.updated_at,''), j.created_at) DESC, "
            "j.created_at DESC, j.run_id DESC"
        ),
        "finished_desc": (
            "CASE WHEN lower(j.status) IN ('completed','succeeded','failed','cancelled','canceled') "
            "THEN 0 ELSE 1 END ASC, "
            "COALESCE(NULLIF(j.completed_at,''), NULLIF(j.updated_at,''), j.created_at) DESC, "
            "j.created_at DESC, j.run_id DESC"
        ),
        "name_asc": (
            "lower(COALESCE(NULLIF(j.campaign_name,''), NULLIF(j.tool,''), NULLIF(j.job_type,''), j.run_id)) ASC, "
            "j.created_at DESC, j.run_id DESC"
        ),
        "status_asc": "replace(lower(j.status), '_', ' ') ASC, j.created_at DESC, j.run_id DESC",
    }.get(sort_by)
    if order_by is None:
        raise ValueError("sort_by must be created_desc, active_desc, finished_desc, name_asc, or status_asc")
    if cursor is not None:
        try:
            cursor_offset = int(cursor)
        except (TypeError, ValueError) as exc:
            raise ValueError("cursor must be a non-negative offset returned by this endpoint") from exc
        if cursor_offset < 0:
            raise ValueError("cursor must be a non-negative offset returned by this endpoint")
        offset = cursor_offset
    root = _results_root()
    _reconcile_for_read(root)
    clauses = ["j.storage_root_id=?"]
    params: list[Any] = [_root_id(root)]
    if not include_hidden:
        clauses.append("j.hidden=0 AND j.capacity_parent=0")
    if task_group:
        clauses.append("j.task_group=?")
        params.append(task_group)
    if not access_admin:
        access_clauses: list[str] = []
        if access_owner_id:
            access_clauses.append("j.owner_id=?")
            params.append(access_owner_id)
        selected_access_projects = sorted({value for value in (accessible_project_ids or []) if value})
        if selected_access_projects:
            access_clauses.append("j.project_id IN (" + ",".join("?" for _ in selected_access_projects) + ")")
            params.extend(selected_access_projects)
        clauses.append("(" + " OR ".join(access_clauses) + ")" if access_clauses else "0=1")
    if owner_id:
        clauses.append("j.owner_id=?")
        params.append(owner_id)
    if project_id == "__none__":
        clauses.append("(j.project_id IS NULL OR j.project_id='')")
    elif project_id:
        clauses.append("j.project_id=?")
        params.append(project_id)
    if statuses:
        clauses.append("lower(j.status) IN (" + ",".join("?" for _ in statuses) + ")")
        params.extend(value.casefold() for value in statuses)
    if tools:
        placeholders = ",".join("?" for _ in tools)
        clauses.append(f"(j.tool IN ({placeholders}) OR (j.tool='' AND j.job_type IN ({placeholders})))")
        params.extend(tools)
        params.extend(tools)
    expression = _search_expression(query)
    if expression:
        clauses.append("j.id IN (SELECT rowid FROM jobs_fts WHERE jobs_fts MATCH ?)")
        params.append(expression)
    where = " AND ".join(clauses)
    try:
        with _connection() as connection:
            total = int(connection.execute("SELECT count(*) FROM jobs j WHERE " + where, params).fetchone()[0])
            rows = connection.execute(
                "SELECT j.* FROM jobs j WHERE " + where
                + " ORDER BY " + order_by + " LIMIT ? OFFSET ?",
                [*params, limit, offset],
            ).fetchall()
    except sqlite3.OperationalError as exc:
        if "jobs_fts" not in str(exc):
            raise
        like_clauses = clauses[:-1] + [
            "lower(j.campaign_name || ' ' || j.tool || ' ' || j.job_type || ' ' || j.job_code || ' ' "
            "|| j.run_id || ' ' || j.task_group || ' ' || j.status) LIKE ?"
        ]
        like_params = params[:-1] + [f"%{query.casefold()}%"]
        where = " AND ".join(like_clauses)
        with _connection() as connection:
            total = int(connection.execute("SELECT count(*) FROM jobs j WHERE " + where, like_params).fetchone()[0])
            rows = connection.execute(
                "SELECT j.* FROM jobs j WHERE " + where
                + " ORDER BY " + order_by + " LIMIT ? OFFSET ?",
                [*like_params, limit, offset],
            ).fetchall()
    items = []
    for row in rows:
        item = dict(row)
        item.update(json.loads(item.pop("run_identity_json", "{}")))
        item["run_dir"] = str((root / item.pop("run_relative_path")).resolve())
        item["success"] = None if item["success"] is None else bool(item["success"])
        item["hidden"] = bool(item["hidden"])
        item["capacity_parent"] = bool(item["capacity_parent"])
        if item.get("job_type") in {"library_preparation", "library_centroids"}:
            # Keep the distributed-library worker details in the compact
            # catalog response without adding their operational metrics to
            # the SQLite index schema.
            run_dir = Path(item["run_dir"])
            metadata = _read_json(run_dir / "metadata.json")
            input_payload = _read_json(run_dir / "input.json")
            input_params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
            result_payload = _read_json(run_dir / "result.json")
            result_metrics = result_payload.get("metrics") if isinstance(result_payload.get("metrics"), dict) else {}
            worker_metrics = result_metrics.get("workers") if isinstance(result_metrics.get("workers"), list) else []
            item["campaign_name"] = item.get("campaign_name") or str(input_params.get("library_name") or "")
            item["worker_plan"] = metadata.get("worker_plan") if isinstance(metadata.get("worker_plan"), list) else []
            item["metrics"] = {"workers": worker_metrics} if worker_metrics else {}
        display_tool_params = {}
        if item["tool"] == "ligandmpnn":
            input_payload = _read_json(Path(item["run_dir"]) / "input.json")
            display_tool_params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
        item["display_tool"] = job_display_tool(item["tool"], display_tool_params)
        items.append(item)
    return {
        "items": items,
        "total": total,
        "next_cursor": str(offset + limit) if offset + limit < total else None,
    }


def job_list_facets(
    *,
    access_owner_id: str | None = None,
    accessible_project_ids: list[str] | None = None,
    access_admin: bool = False,
) -> dict[str, list[str]]:
    root = _results_root()
    _reconcile_for_read(root)
    params: list[Any] = [_root_id(root)]
    clauses = ["storage_root_id=?", "hidden=0", "capacity_parent=0"]
    if not access_admin:
        access_clauses: list[str] = []
        if access_owner_id:
            access_clauses.append("owner_id=?")
            params.append(access_owner_id)
        selected_access_projects = sorted({value for value in (accessible_project_ids or []) if value})
        if selected_access_projects:
            access_clauses.append("project_id IN (" + ",".join("?" for _ in selected_access_projects) + ")")
            params.extend(selected_access_projects)
        clauses.append("(" + " OR ".join(access_clauses) + ")" if access_clauses else "0=1")
    with _connection() as connection:
        rows = connection.execute(
            "SELECT DISTINCT status, tool, job_type FROM jobs WHERE " + " AND ".join(clauses),
            params,
        ).fetchall()
    return {
        "statuses": sorted({str(row["status"] or "unknown").title() for row in rows}),
        "tools": sorted({str(row["tool"] or row["job_type"] or "Unknown") for row in rows}),
    }


def get_job_access_record(task_group: str, run_id: str) -> dict[str, Any] | None:
    """Read the catalog identity fields used to authorize a particular run."""
    root = _results_root()
    _reconcile_for_read(root)
    with _connection() as connection:
        row = connection.execute(
            "SELECT owner_id, project_id, task_group, run_id, run_relative_path FROM jobs "
            "WHERE storage_root_id=? AND task_group=? AND run_id=?",
            (_root_id(root), task_group, run_id),
        ).fetchone()
    return dict(row) if row else None


def path_is_within_accessible_job(
    path: str | Path,
    *,
    access_owner_id: str,
    accessible_project_ids: list[str] | None = None,
    access_admin: bool = False,
) -> bool:
    """Check whether a local path belongs to a job visible to this account."""
    selected_path = Path(path).expanduser().resolve()
    if access_admin:
        return True
    root = _results_root()
    _reconcile_for_read(root)
    project_ids = sorted({value for value in (accessible_project_ids or []) if value})
    access_clauses = ["j.owner_id=?"]
    clauses = ["j.storage_root_id=?"]
    params: list[Any] = [_root_id(root), access_owner_id]
    if project_ids:
        access_clauses.append("j.project_id IN (" + ",".join("?" for _ in project_ids) + ")")
        params.extend(project_ids)
    clauses.append("(" + " OR ".join(access_clauses) + ")")
    with _connection() as connection:
        rows = connection.execute(
            "SELECT s.path AS storage_path, j.run_relative_path FROM jobs j "
            "JOIN storage_roots s ON s.id=j.storage_root_id WHERE " + " ".join(clauses),
            params,
        ).fetchall()
    for row in rows:
        try:
            run_dir = (Path(row["storage_path"]) / str(row["run_relative_path"])).resolve()
            selected_path.relative_to(run_dir)
            return True
        except (OSError, ValueError):
            continue
    return False


def job_artifacts(task_group: str, run_id: str) -> list[dict[str, Any]]:
    root = _results_root()
    _reconcile_for_read(root)
    with _connection() as connection:
        rows = connection.execute(
            "SELECT a.* FROM job_artifacts a JOIN jobs j ON j.id=a.job_id "
            "WHERE j.storage_root_id=? AND j.task_group=? AND j.run_id=? ORDER BY a.relative_path",
            (_root_id(root), task_group, run_id),
        ).fetchall()
        job = connection.execute(
            "SELECT run_relative_path FROM jobs WHERE storage_root_id=? AND task_group=? AND run_id=?",
            (_root_id(root), task_group, run_id),
        ).fetchone()
    # Runs imported from archives or copied into the results directory after
    # their first catalog scan can have a job row but no artifact rows. Rebuild
    # that run's artifact index lazily so its files remain discoverable.
    if not rows and job:
        run_dir = (root / str(job["run_relative_path"])).resolve()
        try:
            run_dir.relative_to(root)
        except ValueError:
            return []
        if run_dir.is_dir():
            _upsert(root, run_dir, refresh_artifacts=True)
            with _connection() as connection:
                rows = connection.execute(
                    "SELECT a.* FROM job_artifacts a JOIN jobs j ON j.id=a.job_id "
                    "WHERE j.storage_root_id=? AND j.task_group=? AND j.run_id=? ORDER BY a.relative_path",
                    (_root_id(root), task_group, run_id),
                ).fetchall()
    return [dict(row) for row in rows]
