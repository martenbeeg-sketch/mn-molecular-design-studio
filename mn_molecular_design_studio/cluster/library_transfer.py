"""Resumable, integrity-checked transfer of completed external-library outputs."""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
from pathlib import Path

from . import config
from .store import TERMINAL, Store, encoded, stamp

CHUNK_BYTES = 8 * 1024 * 1024
_LIBRARY_ID = re.compile(r"lib-[0-9a-f]{24}\Z")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_SHARD = re.compile(r"features-[0-9]{6}\.parquet\Z")
_CENTROID_FILES = {"centroids.parquet", "assignments.parquet"}


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _row(db, job_id, node_id):
    transfer = db.execute("SELECT * FROM library_transfers WHERE job_id=?", (job_id,)).fetchone()
    if not transfer:
        raise LookupError("No organizer output transfer is registered for this job.")
    if transfer["node_id"] != node_id:
        raise PermissionError("Only the compute workstation assigned to this transfer may upload it.")
    job = db.execute("SELECT node_id,payload FROM jobs WHERE id=?", (job_id,)).fetchone()
    if not job or job["node_id"] != node_id:
        raise PermissionError("This workstation does not own the library preparation job.")
    payload = json.loads(job["payload"])
    if payload.get("tool") not in {"library_preparation", "library_centroids"} or payload.get("status") not in TERMINAL:
        raise ValueError("Fingerprint outputs can be transferred after the library preparation job finishes.")
    return dict(transfer), payload


def _safe_files(raw):
    if not isinstance(raw, list) or not raw:
        raise ValueError("The transfer needs a list of generated output files.")
    files = []
    seen = set()
    for item in raw:
        if not isinstance(item, dict) or set(item) != {"path", "size_bytes", "sha256"}:
            raise ValueError("Invalid generated output file metadata.")
        name, size, digest = item["path"], item["size_bytes"], item["sha256"]
        if (not isinstance(name, str) or name in seen
                or not (name in {"manifest.json", "rejected.parquet", *_CENTROID_FILES} or _SHARD.fullmatch(name))):
            raise ValueError("The output transfer contains an unsupported file name.")
        if type(size) is not int or size < 0 or not isinstance(digest, str) or not _SHA256.fullmatch(digest):
            raise ValueError("Invalid output file size or checksum.")
        seen.add(name)
        files.append({"path": name, "size_bytes": size, "sha256": digest})
    if "manifest.json" not in seen:
        raise ValueError("Transfer the generated library manifest.")
    if not ("rejected.parquet" in seen or "centroids.parquet" in seen):
        raise ValueError("Transfer the fingerprint outputs or the centroid table.")
    return sorted(files, key=lambda item: item["path"])


def _stage_path(transfer):
    destination = Path(transfer["destination"])
    return destination.parent / f".{destination.name}.studio-transfer-{transfer['job_id']}.partial"


def _sizes(files, directory):
    offsets = []
    transferred = 0
    for item in files:
        path = directory / item["path"]
        offset = path.stat().st_size if path.is_file() else 0
        if offset > item["size_bytes"]:
            raise ValueError(f"Partial file {item['path']} is larger than its declared size.")
        transferred += offset
        offsets.append({"path": item["path"], "offset": offset, "size_bytes": item["size_bytes"]})
    return offsets, transferred


def _response(transfer, files, directory):
    offsets, transferred = _sizes(files, directory)
    return {"status": transfer["status"], "files": offsets,
            "total_bytes": sum(item["size_bytes"] for item in files),
            "transferred_bytes": transferred}


def start(job_id, node_id, payload):
    library_id = payload.get("library_id") if isinstance(payload, dict) else None
    if not isinstance(library_id, str) or not _LIBRARY_ID.fullmatch(library_id):
        raise ValueError("Invalid external library ID in transfer request.")
    files = _safe_files(payload.get("files"))
    total_bytes = sum(item["size_bytes"] for item in files)
    store = Store()
    with store.connect() as db:
        db.execute("BEGIN IMMEDIATE")
        transfer, job = _row(db, job_id, node_id)
        if transfer["status"] == "completed":
            if transfer["library_id"] != library_id or json.loads(transfer["files"]) != files:
                raise ValueError("The completed transfer cannot be changed.")
            return _response(transfer, files, Path(transfer["destination"]))
        previous_files = json.loads(transfer["files"])
        if previous_files and (transfer["library_id"] != library_id or previous_files != files):
            raise ValueError("A resumed transfer must use the original library and file manifest.")
        if not previous_files:
            db.execute("UPDATE library_transfers SET library_id=?,files=?,total_bytes=?,status='receiving',error=NULL,updated_at=? WHERE job_id=?",
                       (library_id, encoded(files), total_bytes, stamp(), job_id))
            transfer.update(library_id=library_id, files=encoded(files), total_bytes=total_bytes, status="receiving")
        destination = Path(transfer["destination"])
        if destination.exists() and transfer["status"] not in {"verifying", "failed"}:
            raise ValueError("The organizer destination already exists; Studio will not overwrite it.")
        destination.parent.mkdir(parents=True, exist_ok=True)
        stage = _stage_path(transfer)
        data_directory = destination if destination.exists() else stage
        if not destination.exists():
            stage.mkdir(parents=True, exist_ok=True)
        offsets, transferred = _sizes(files, data_directory)
        free = shutil.disk_usage(destination.parent).free
        if not destination.exists() and free < total_bytes - transferred:
            db.execute("UPDATE library_transfers SET status='waiting_for_space',transferred_bytes=?,error=?,updated_at=? WHERE job_id=?",
                       (transferred, "Not enough free disk space on the organizer.", stamp(), job_id))
            return {"status": "waiting_for_space", "files": offsets, "total_bytes": total_bytes,
                    "transferred_bytes": transferred, "error": "Not enough free disk space on the organizer."}
        status = "verifying" if destination.exists() else "receiving"
        db.execute("UPDATE library_transfers SET status=?,total_bytes=?,transferred_bytes=?,error=NULL,updated_at=? WHERE job_id=?",
                   (status, total_bytes, transferred, stamp(), job_id))
        transfer.update(status=status, total_bytes=total_bytes, transferred_bytes=transferred)
        return _response(transfer, files, destination if status == "verifying" else stage)


def status(job_id, node_id):
    store = Store()
    with store.connect() as db:
        transfer, _ = _row(db, job_id, node_id)
    files = json.loads(transfer["files"])
    destination = Path(transfer["destination"])
    directory = destination if destination.exists() else _stage_path(transfer)
    return _response(transfer, files, directory)


async def put_chunk(job_id, node_id, name, offset, request):
    store = Store()
    with store.connect() as db:
        transfer, _ = _row(db, job_id, node_id)
    if transfer["status"] == "completed":
        raise ValueError("The output transfer is already complete.")
    files = json.loads(transfer["files"])
    item = next((entry for entry in files if entry["path"] == name), None)
    if not item:
        raise ValueError("This file is not part of the registered output transfer.")
    if type(offset) is not int or offset < 0:
        raise ValueError("Invalid upload offset.")
    stage = _stage_path(transfer)
    stage.mkdir(parents=True, exist_ok=True)
    target = stage / name
    current = target.stat().st_size if target.is_file() else 0
    if current != offset:
        raise ValueError(f"Upload offset mismatch; organizer currently has {current} bytes of {name}.")
    written = 0
    try:
        mode = "r+b" if target.exists() else "wb"
        with target.open(mode) as handle:
            handle.seek(offset)
            async for block in request.stream():
                written += len(block)
                if written > CHUNK_BYTES or offset + written > item["size_bytes"]:
                    raise ValueError("An upload request must contain at most one 8 MiB chunk and stay within the declared file size.")
                handle.write(block)
            if not written:
                raise ValueError("Empty output chunks are not accepted.")
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        if target.exists():
            with target.open("r+b") as handle:
                handle.truncate(offset)
        raise
    _offsets, transferred = _sizes(files, stage)
    with store.connect() as db:
        db.execute("UPDATE library_transfers SET status='receiving',transferred_bytes=?,error=NULL,updated_at=? WHERE job_id=? AND node_id=?",
                   (transferred, stamp(), job_id, node_id))
    return {"received_bytes": written, "offset": offset + written, "transferred_bytes": transferred}


def _verify(directory, transfer, files):
    by_name = {item["path"]: item for item in files}
    for item in files:
        path = directory / item["path"]
        if not path.is_file() or path.stat().st_size != item["size_bytes"] or _digest(path) != item["sha256"]:
            raise ValueError(f"Generated output {item['path']} failed its size or SHA-256 integrity check.")
    try:
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("The transferred fingerprint manifest is invalid.") from exc
    if (manifest.get("schema_version") not in {"studio.library-fingerprints.v2", "studio.library-centroids.v2"}
            or manifest.get("status") != "completed" or manifest.get("library_id") != transfer["library_id"]):
        raise ValueError("The transferred manifest does not describe the completed library preparation.")
    expected_names = {"manifest.json"}
    if manifest["schema_version"] == "studio.library-fingerprints.v2":
        shard_rows = manifest.get("shards")
        rejected = manifest.get("rejections")
        if not isinstance(shard_rows, list) or not isinstance(rejected, dict):
            raise ValueError("The fingerprint manifest is missing its shard or rejected-row records.")
        expected_names.add("rejected.parquet")
        for shard in shard_rows:
            if not isinstance(shard, dict):
                raise ValueError("Invalid shard metadata in the fingerprint manifest.")
            item = by_name.get(shard.get("path"))
            if (not item or item["sha256"] != shard.get("sha256")
                    or item["size_bytes"] != shard.get("size_bytes")):
                raise ValueError("The transferred shard list disagrees with the fingerprint manifest.")
            expected_names.add(item["path"])
        rejected_item = by_name.get("rejected.parquet")
        if not rejected_item or rejected_item["sha256"] != rejected.get("sha256") or rejected.get("path") != "rejected.parquet":
            raise ValueError("The rejected-row table disagrees with the fingerprint manifest.")
    else:
        declared = manifest.get("files")
        if not isinstance(declared, list):
            raise ValueError("The centroid manifest is missing its output files.")
        for entry in declared:
            if not isinstance(entry, dict):
                raise ValueError("Invalid output file metadata in the centroid manifest.")
            item = by_name.get(entry.get("path"))
            if not item or item["sha256"] != entry.get("sha256"):
                raise ValueError("The transferred centroid files disagree with their manifest.")
            expected_names.add(item["path"])
    if expected_names != set(by_name):
        raise ValueError("The output transfer contains files not listed in the fingerprint manifest.")


def complete(job_id, node_id):
    store = Store()
    with store.connect() as db:
        transfer, _ = _row(db, job_id, node_id)
    files = json.loads(transfer["files"])
    destination = Path(transfer["destination"])
    stage = _stage_path(transfer)
    directory = destination if destination.exists() else stage
    try:
        _verify(directory, transfer, files)
        if directory == stage:
            if destination.exists():
                raise ValueError("The organizer destination appeared during transfer; no files were overwritten.")
            with store.connect() as db:
                db.execute("UPDATE library_transfers SET status='verifying',error=NULL,updated_at=? WHERE job_id=? AND node_id=?",
                           (stamp(), job_id, node_id))
            os.replace(stage, destination)
        from mn_molecular_design_studio.external_libraries import LibraryCatalog
        LibraryCatalog().attach_transferred_preparation(transfer["library_id"], destination / "manifest.json",
                                                         source_node_id=node_id)
    except Exception as exc:
        with store.connect() as db:
            db.execute("UPDATE library_transfers SET status='failed',error=?,updated_at=? WHERE job_id=? AND node_id=?",
                       (str(exc)[:4000], stamp(), job_id, node_id))
        raise
    with store.connect() as db:
        db.execute("UPDATE library_transfers SET status='completed',transferred_bytes=total_bytes,error=NULL,updated_at=? WHERE job_id=? AND node_id=?",
                   (stamp(), job_id, node_id))
    return {"status": "completed", "directory": str(destination),
            "total_bytes": sum(item["size_bytes"] for item in files)}
