"""Verified, immutable transport snapshots for retained scientific artifacts."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile


def authorize_result_artifact(user, structure_path, *, source_run_dir, runs_root):
    """Require authorization and catalogue membership for linked run artifacts."""
    from mn_molecular_design_studio.auth import can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record, job_artifacts

    try:
        relative = structure_path.relative_to(runs_root)
    except ValueError as exc:
        raise ValueError("The selected structure is outside Studio's managed Results storage.") from exc
    if structure_path.is_relative_to(source_run_dir):
        return source_run_dir.parent.name, source_run_dir.name, structure_path.relative_to(source_run_dir).as_posix()
    if len(relative.parts) < 3:
        raise ValueError("The selected structure is not inside a registered Results artifact.")
    task_group, run_id = relative.parts[:2]
    access = get_job_access_record(task_group, run_id)
    if not access or not can_access_job(user, access):
        raise PermissionError("The structure-producing Results job is unavailable or you do not have access to it.")
    for artifact in job_artifacts(task_group, run_id):
        artifact_relative = Path(str(artifact.get("relative_path") or ""))
        if artifact_relative.is_absolute() or ".." in artifact_relative.parts:
            continue
        try:
            path = (runs_root / artifact_relative).resolve(strict=True)
            path.relative_to(runs_root / task_group / run_id)
        except (OSError, ValueError):
            continue
        if path == structure_path:
            return task_group, run_id, f"runs:///{relative.as_posix()}"
    raise ValueError("The selected structure is not a registered artifact of its producing job.")


def verified_snapshot(path, expected):
    from mn_molecular_design_studio.runtime import app_home
    if not isinstance(expected, str) or len(expected) != 64 or any(char not in "0123456789abcdef" for char in expected):
        raise ValueError("The retained artifact lacks a valid SHA-256 identity.")
    path = Path(path)
    root = Path(app_home()) / "verified-artifacts"
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    destination = root / (expected + path.suffix.lower())
    temporary_path = None
    try:
        digest = hashlib.sha256()
        with path.open("rb") as source, tempfile.NamedTemporaryFile(dir=root, delete=False) as temporary:
            temporary_path = Path(temporary.name)
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(chunk)
                temporary.write(chunk)
            temporary.flush()
            os.fchmod(temporary.fileno(), 0o400)
        if digest.hexdigest() != expected:
            raise ValueError("Retained scientific artifact digest mismatch; refusing changed coordinates.")
        os.replace(temporary_path, destination)
        return destination
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def retained_artifact_file(run_dir, path):
    """Enforce recorded identities, retaining raw access for unversioned history."""
    root = Path(run_dir).resolve()
    selected = Path(path).resolve()
    manifest_path = root / "artifacts.json"
    if not manifest_path.is_file():
        return selected
    if not manifest_path.resolve().is_relative_to(root):
        raise ValueError("Artifact manifest escapes its Results run.")
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest, dict) or manifest.get("schema_version") != 1 or manifest.get("run_id") != root.name:
        raise ValueError("Artifact manifest schema or run identity is invalid.")
    if not isinstance(manifest.get("artifacts"), list) or any(not isinstance(artifact, dict) for artifact in manifest["artifacts"]):
        raise ValueError("Artifact manifest entries are invalid.")
    matches = []
    for artifact in manifest.get("artifacts") or []:
        relative = Path(str(artifact.get("path") or ""))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Artifact manifest contains an unsafe path.")
        if (root / relative).resolve() == selected:
            matches.append(artifact)
    digests = {str(artifact.get("sha256") or "") for artifact in matches if artifact.get("sha256")}
    if len(digests) > 1:
        raise ValueError("Artifact manifest contains conflicting scientific identities.")
    if any(artifact.get("run_id") and artifact["run_id"] != root.name for artifact in matches):
        raise ValueError("Artifact run identity does not match its manifest.")
    if any(artifact.get("size_bytes") is not None and int(artifact["size_bytes"]) != selected.stat().st_size
           for artifact in matches):
        raise ValueError("Retained scientific artifact size mismatch.")
    return verified_snapshot(selected, next(iter(digests))) if digests else selected
