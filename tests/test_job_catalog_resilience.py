from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path

from mn_molecular_design_studio import job_catalog


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_reconcile_does_not_block_when_another_scan_is_running(tmp_path, monkeypatch):
    lock = threading.Lock()
    lock.acquire()
    monkeypatch.setattr(job_catalog, "_RECONCILE_LOCK", lock)
    monkeypatch.setattr(
        job_catalog,
        "_reconcile_locked",
        lambda _root: (_ for _ in ()).throw(AssertionError("redundant scan ran")),
    )

    job_catalog.reconcile(tmp_path / "runs")

    lock.release()


def test_metadata_refresh_reuses_success_when_result_file_is_unchanged(tmp_path, monkeypatch):
    root = tmp_path / "runs"
    run = root / "design" / "run-1"
    run.mkdir(parents=True)
    _write_json(run / "metadata.json", {"task_group": "design", "status": "completed"})
    _write_json(run / "input.json", {"tool": "rfd4_proteina", "params": {}})
    _write_json(run / "result.json", {"success": True, "metrics": {"large": [1, 2, 3]}})
    monkeypatch.setattr(job_catalog, "catalog_path", lambda: tmp_path / "catalog.sqlite3")
    monkeypatch.setattr(job_catalog, "_SCHEMA_READY_PATH", None)

    job_catalog._upsert(root, run, refresh_artifacts=False)
    result_mtime = (run / "result.json").stat().st_mtime_ns
    metadata_path = run / "metadata.json"
    _write_json(metadata_path, {"task_group": "design", "status": "completed", "current_phase": "done"})
    os.utime(metadata_path, ns=(metadata_path.stat().st_atime_ns, metadata_path.stat().st_mtime_ns + 1_000_000_000))

    read_json = job_catalog._read_json

    def reject_unchanged_result_read(path):
        if Path(path).name == "result.json":
            raise AssertionError("unchanged result.json was parsed again")
        return read_json(path)

    monkeypatch.setattr(job_catalog, "_read_json", reject_unchanged_result_read)
    job_catalog._upsert(root, run, refresh_artifacts=False)

    with job_catalog._connection() as connection:
        row = connection.execute(
            "SELECT success, result_mtime_ns, current_phase FROM jobs WHERE task_group='design' AND run_id='run-1'"
        ).fetchone()
    assert row["success"] == 1
    assert row["result_mtime_ns"] == result_mtime
    assert row["current_phase"] == "done"


def test_large_result_indexes_root_success_without_loading_candidate_payload(tmp_path, monkeypatch):
    result_path = tmp_path / "result.json"
    monkeypatch.setattr(job_catalog, "_RESULT_INLINE_READ_LIMIT", 128)
    monkeypatch.setattr(job_catalog, "_RESULT_TAIL_READ_BYTES", 256)
    result_path.write_text(
        json.dumps(
            {
                "outputs": {"candidate": {"success": False, "structure": "x" * 2_000_000}},
                "success": True,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    success, workers = job_catalog._read_result_index_fields(result_path)

    assert success is True
    assert workers == []


def test_large_result_without_root_success_does_not_use_nested_success(tmp_path, monkeypatch):
    result_path = tmp_path / "result.json"
    monkeypatch.setattr(job_catalog, "_RESULT_INLINE_READ_LIMIT", 1)
    monkeypatch.setattr(job_catalog, "_RESULT_TAIL_READ_BYTES", 128)
    result_path.write_text(
        json.dumps({"outputs": {"candidate": {"success": True, "structure": "x" * 1_000}}, "status": "done"}, indent=2),
        encoding="utf-8",
    )

    success, workers = job_catalog._read_result_index_fields(result_path)

    assert success is None
    assert workers == []


def test_indexed_search_returns_while_catalog_refresh_runs_in_background(tmp_path, monkeypatch):
    root = tmp_path / "runs"
    run = root / "design" / "run-1"
    run.mkdir(parents=True)
    _write_json(run / "metadata.json", {"task_group": "design", "status": "completed"})
    _write_json(run / "input.json", {"tool": "rfd4_proteina", "params": {}})
    _write_json(run / "result.json", {"success": True})
    monkeypatch.setattr(job_catalog, "catalog_path", lambda: tmp_path / "catalog.sqlite3")
    monkeypatch.setattr(job_catalog, "_results_root", lambda: root)
    monkeypatch.setattr(job_catalog, "_SCHEMA_READY_PATH", None)
    job_catalog._upsert(root, run, refresh_artifacts=False)

    started = threading.Event()
    finish = threading.Event()

    def slow_scan(_root):
        started.set()
        finish.wait(timeout=3)

    monkeypatch.setattr(job_catalog, "_reconcile_locked", slow_scan)
    root_key = str(root.resolve())
    with job_catalog._RECONCILE_SCHEDULE_LOCK:
        job_catalog._RECONCILE_SCHEDULED.discard(root_key)

    before = time.monotonic()
    page = job_catalog.search_jobs(access_admin=True, limit=10)
    elapsed = time.monotonic() - before
    try:
        assert page["total"] == 1
        assert elapsed < 0.5
        assert started.wait(timeout=1)
    finally:
        finish.set()
        deadline = time.monotonic() + 1
        while time.monotonic() < deadline:
            with job_catalog._RECONCILE_SCHEDULE_LOCK:
                if root_key not in job_catalog._RECONCILE_SCHEDULED:
                    break
            time.sleep(0.01)
