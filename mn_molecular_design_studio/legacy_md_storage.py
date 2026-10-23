"""Safe, run-scoped copies of legacy molecular-dynamics artifacts.

This module only stages bytes under a caller-selected Results run. It does not
register jobs, change ownership, or interact with the Molecular Database.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import socket
import stat
import tempfile
from pathlib import Path, PurePosixPath
from typing import Any


_MANIFEST_NAME = "legacy_copy_manifest.json"
_CHUNK_SIZE = 1024 * 1024


class LegacyMDStorageError(ValueError):
    """Raised when a legacy run cannot be copied safely and consistently."""


def _safe_run_id(value: str) -> str:
    run_id = str(value or "")
    if (not run_id or run_id in {".", ".."} or "/" in run_id or "\\" in run_id
            or Path(run_id).name != run_id or "\x00" in run_id):
        raise LegacyMDStorageError("Source run ID must be a single safe path component.")
    return run_id


def _contained(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def _stat_signature(value: os.stat_result) -> tuple[int, ...]:
    return (value.st_dev, value.st_ino, value.st_mode, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)


def _entry_inventory(source_root: Path) -> list[dict[str, Any]]:
    """Inventory every regular file, following only safe in-run file links."""
    entries: list[dict[str, Any]] = []
    for current, dirnames, filenames in os.walk(source_root, followlinks=False):
        current_path = Path(current)
        for dirname in list(dirnames):
            path = current_path / dirname
            if path.is_symlink():
                try:
                    target = path.resolve(strict=True)
                except (OSError, RuntimeError) as exc:
                    raise LegacyMDStorageError(f"Broken or invalid symlink: {path.relative_to(source_root)}") from exc
                if not _contained(target, source_root):
                    raise LegacyMDStorageError(f"Symlink escapes source run: {path.relative_to(source_root)}")
                if not target.is_dir():
                    raise LegacyMDStorageError(f"Symlinked directory is not a directory: {path.relative_to(source_root)}")
                raise LegacyMDStorageError(f"Symlinked directories are not supported: {path.relative_to(source_root)}")
            mode = path.lstat().st_mode
            if not stat.S_ISDIR(mode):
                raise LegacyMDStorageError(f"Unsupported filesystem entry: {path.relative_to(source_root)}")
        for filename in filenames:
            path = current_path / filename
            relative = path.relative_to(source_root)
            posix = PurePosixPath(relative.as_posix())
            if posix.is_absolute() or any(part in {"", ".", ".."} for part in posix.parts):
                raise LegacyMDStorageError(f"Unsafe relative path in source run: {relative}")
            if relative.as_posix() == _MANIFEST_NAME:
                raise LegacyMDStorageError(f"Source run uses reserved manifest path: {_MANIFEST_NAME}")
            link_stat = path.lstat()
            is_link = stat.S_ISLNK(link_stat.st_mode)
            try:
                actual_path = path.resolve(strict=True) if is_link else path
            except (OSError, RuntimeError) as exc:
                raise LegacyMDStorageError(f"Broken or invalid symlink: {relative}") from exc
            if is_link and not _contained(actual_path, source_root):
                raise LegacyMDStorageError(f"Symlink escapes source run: {relative}")
            actual_stat = actual_path.stat()
            if not stat.S_ISREG(actual_stat.st_mode):
                raise LegacyMDStorageError(f"Unsupported non-regular file: {relative}")
            entries.append({
                "relative_path": relative.as_posix(),
                "size_bytes": actual_stat.st_size,
                "_source_path": actual_path,
                "_link_path": path,
                "_link_signature": _stat_signature(link_stat),
                "_source_signature": _stat_signature(actual_stat),
                "_is_link": is_link,
            })
    entries.sort(key=lambda row: row["relative_path"])
    return entries


def legacy_md_run_inventory(source_run_dir: str | Path) -> dict[str, Any]:
    """Return a dry-run inventory and byte total without creating files."""
    source_root = Path(source_run_dir).expanduser().resolve(strict=True)
    if not source_root.is_dir():
        raise LegacyMDStorageError("Source run path must be a directory.")
    entries = _entry_inventory(source_root)
    return {
        "source_run_path": str(source_root),
        "files": [{"relative_path": row["relative_path"], "size_bytes": row["size_bytes"]} for row in entries],
        "total_size_bytes": sum(row["size_bytes"] for row in entries),
    }


def _sha256_file(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK_SIZE):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def _manifest_payload(source_root: Path, source_run_id: str, source_workstation: str,
                      files: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "source_run_id": source_run_id,
        "source_run_path": str(source_root),
        "source_workstation": source_workstation,
        "files": files,
    }


def _validate_existing(destination: Path, expected: dict[str, Any]) -> None:
    manifest_path = destination / _MANIFEST_NAME
    try:
        actual = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise LegacyMDStorageError("Existing legacy copy has no valid manifest.") from exc
    if actual != expected:
        raise LegacyMDStorageError("Existing legacy copy conflicts with the requested source manifest.")

    expected_paths = {row["relative_path"] for row in expected["files"]}
    actual_paths: set[str] = set()
    for current, dirnames, filenames in os.walk(destination, followlinks=False):
        current_path = Path(current)
        for dirname in dirnames:
            entry = current_path / dirname
            if entry.is_symlink() or not entry.is_dir():
                raise LegacyMDStorageError("Existing legacy copy contains an unsafe directory entry.")
        for filename in filenames:
            entry = current_path / filename
            if entry.is_symlink() or not entry.is_file():
                raise LegacyMDStorageError("Existing legacy copy contains an unsafe file entry.")
            relative = entry.relative_to(destination).as_posix()
            if relative != _MANIFEST_NAME:
                actual_paths.add(relative)
    if actual_paths != expected_paths:
        raise LegacyMDStorageError("Existing legacy copy has missing or unexpected files.")
    for record in expected["files"]:
        path = destination / record["relative_path"]
        digest, size = _sha256_file(path)
        if size != record["size_bytes"] or digest != record["sha256"]:
            raise LegacyMDStorageError(f"Existing legacy copy failed byte verification: {record['relative_path']}")


def stage_legacy_md_run(
    source_run_dir: str | Path,
    destination_run_dir: str | Path,
    *,
    source_run_id: str,
    source_workstation: str | None = None,
) -> dict[str, Any]:
    """Copy a full legacy run into ``artifacts/legacy/<source_run_id>``.

    The source is never modified. Existing identical copies are fully verified
    and returned with ``skipped=True``; conflicts and unsafe sources fail.
    """
    run_id = _safe_run_id(source_run_id)
    source_root = Path(source_run_dir).expanduser().resolve(strict=True)
    if not source_root.is_dir():
        raise LegacyMDStorageError("Source run path must be a directory.")

    destination_run = Path(destination_run_dir).expanduser().resolve(strict=False)
    destination = destination_run / "artifacts" / "legacy" / run_id
    if _contained(destination_run, source_root) or _contained(source_root, destination_run):
        raise LegacyMDStorageError("Source and destination run directories must not overlap.")
    if not _contained(destination.parent.resolve(strict=False), destination_run):
        raise LegacyMDStorageError("Destination artifact path escapes the destination run.")
    # Inventory and expected hashes are computed before accepting an existing
    # destination; never trust its manifest alone.
    entries = _entry_inventory(source_root)

    source_workstation = str(source_workstation or socket.gethostname())
    files: list[dict[str, Any]] = []
    for row in entries:
        source_path = row["_source_path"]
        link_path = row["_link_path"]
        if _stat_signature(link_path.lstat()) != row["_link_signature"]:
            raise LegacyMDStorageError(f"Source changed during inventory: {row['relative_path']}")
        if row["_is_link"]:
            resolved = link_path.resolve(strict=True)
            if not _contained(resolved, source_root):
                raise LegacyMDStorageError(f"Symlink escapes source run: {row['relative_path']}")
            source_path = resolved
        if _stat_signature(source_path.stat()) != row["_source_signature"]:
            raise LegacyMDStorageError(f"Source changed during inventory: {row['relative_path']}")
        digest, size = _sha256_file(source_path)
        if _stat_signature(source_path.stat()) != row["_source_signature"] or size != row["size_bytes"]:
            raise LegacyMDStorageError(f"Source changed while hashing: {row['relative_path']}")
        if _stat_signature(link_path.lstat()) != row["_link_signature"]:
            raise LegacyMDStorageError(f"Source link changed while hashing: {row['relative_path']}")
        files.append({"relative_path": row["relative_path"], "size_bytes": size, "sha256": digest})

    expected = _manifest_payload(source_root, run_id, source_workstation, files)
    if destination.exists() or destination.is_symlink():
        if destination.is_symlink() or not destination.is_dir():
            raise LegacyMDStorageError("Destination already exists and is not a regular directory.")
        _validate_existing(destination, expected)
        return {"destination": destination, "files": files, "total_size_bytes": sum(x["size_bytes"] for x in files), "skipped": True}

    parent = destination.parent
    parent.mkdir(parents=True, exist_ok=True)
    if not _contained(parent.resolve(strict=True), destination_run):
        raise LegacyMDStorageError("Destination artifact path escapes the destination run.")
    staging = Path(tempfile.mkdtemp(prefix=f".{run_id}.staging-", dir=parent))
    try:
        for row, record in zip(entries, files, strict=True):
            source_path = row["_source_path"]
            link_path = row["_link_path"]
            if row["_is_link"]:
                source_path = link_path.resolve(strict=True)
                if not _contained(source_path, source_root):
                    raise LegacyMDStorageError(f"Symlink escapes source run: {row['relative_path']}")
            relative = Path(record["relative_path"])
            staged_file = staging / relative
            staged_file.parent.mkdir(parents=True, exist_ok=True)
            before_link = _stat_signature(link_path.lstat())
            before_source = _stat_signature(source_path.stat())
            if before_link != row["_link_signature"] or before_source != row["_source_signature"]:
                raise LegacyMDStorageError(f"Source changed before copy: {record['relative_path']}")
            digest = hashlib.sha256()
            size = 0
            with source_path.open("rb") as source, staged_file.open("xb") as target:
                while chunk := source.read(_CHUNK_SIZE):
                    target.write(chunk)
                    digest.update(chunk)
                    size += len(chunk)
                target.flush()
                os.fsync(target.fileno())
            if (_stat_signature(source_path.stat()) != row["_source_signature"]
                    or _stat_signature(link_path.lstat()) != row["_link_signature"]):
                raise LegacyMDStorageError(f"Source changed during copy: {record['relative_path']}")
            if size != record["size_bytes"] or digest.hexdigest() != record["sha256"]:
                raise LegacyMDStorageError(f"Copied bytes failed SHA-256 verification: {record['relative_path']}")
            persisted_digest, persisted_size = _sha256_file(staged_file)
            if persisted_size != record["size_bytes"] or persisted_digest != record["sha256"]:
                raise LegacyMDStorageError(f"Staged bytes failed readback verification: {record['relative_path']}")

        # Recheck every source after the full copy, since an earlier file could
        # have changed while a later file was being staged.
        for row in entries:
            source_path = row["_source_path"]
            if row["_is_link"]:
                source_path = row["_link_path"].resolve(strict=True)
                if not _contained(source_path, source_root):
                    raise LegacyMDStorageError(f"Symlink escapes source run: {row['relative_path']}")
            if (_stat_signature(row["_link_path"].lstat()) != row["_link_signature"]
                    or _stat_signature(source_path.stat()) != row["_source_signature"]):
                raise LegacyMDStorageError(f"Source changed before commit: {row['relative_path']}")

        manifest_path = staging / _MANIFEST_NAME
        with manifest_path.open("x", encoding="utf-8") as handle:
            json.dump(expected, handle, ensure_ascii=False, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        _validate_existing(staging, expected)
        if destination.exists() or destination.is_symlink():
            raise LegacyMDStorageError("Destination appeared while staging; refusing to overwrite it.")
        os.rename(staging, destination)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    return {"destination": destination, "files": files, "total_size_bytes": sum(x["size_bytes"] for x in files), "skipped": False}
