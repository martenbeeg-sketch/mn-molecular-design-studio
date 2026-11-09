"""Materialize and validate the scientific artifacts retained in a run.

Container paths remain untouched while a tool executes. At Results collection,
candidate structure paths under the conventional ``/work`` mount are mapped
back to the owning run; other container mounts are resolved from the recorded
Docker command and copied into that run when they are candidate outputs.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path, PurePosixPath
from typing import Any, Iterable
from urllib.parse import urlsplit

from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import (
    is_portable_path,
    portable_path,
    resolve_stored_path,
)


_CONTAINER_PREFIXES = ("/work", "/workspace", "/boltz-context")
_SCIENTIFIC_PATH_FIELDS = ("target_pdb", "complex_pdb", "binder_pdb")


def _is_container_path(value: str) -> bool:
    return any(value == prefix or value.startswith(prefix.rstrip("/") + "/") for prefix in _CONTAINER_PREFIXES)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError, TypeError):
        return {}
    return value if isinstance(value, dict) else {}


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _command_lists(value: Any) -> Iterable[list[str]]:
    if isinstance(value, dict):
        for item in value.values():
            yield from _command_lists(item)
    elif isinstance(value, list):
        if value and all(isinstance(item, (str, int, float)) for item in value):
            yield [str(item) for item in value]
        else:
            for item in value:
                yield from _command_lists(item)


def _container_mounts(run_dir: Path) -> list[tuple[str, Path]]:
    """Read bind mounts from saved commands without changing execution args."""
    metadata = _read_json(run_dir / "metadata.json")
    command = _read_json(run_dir / "command.json")
    commands = [*_command_lists(metadata.get("executed_command")),
                *_command_lists(metadata.get("executed_commands")),
                *_command_lists(command.get("command"))]
    mounts: list[tuple[str, Path]] = []
    for argv in commands:
        for index, item in enumerate(argv):
            spec = ""
            if item in {"-v", "--volume"} and index + 1 < len(argv):
                spec = argv[index + 1]
            elif item.startswith("--volume="):
                spec = item.split("=", 1)[1]
            if spec:
                pieces = spec.rsplit(":", 2)
                if len(pieces) >= 2 and pieces[0].startswith("/") and pieces[1].startswith("/"):
                    mounts.append((pieces[1], Path(pieces[0]).expanduser()))
            if item.startswith("--mount="):
                values = {}
                for component in item.split("=", 1)[1].split(","):
                    key, sep, value = component.partition("=")
                    if sep:
                        values[key] = value
                if values.get("type") == "bind" and values.get("source", "").startswith("/") and values.get("target", "").startswith("/"):
                    mounts.append((values["target"], Path(values["source"]).expanduser()))
            if item == "--mount" and index + 1 < len(argv):
                values = {}
                for component in argv[index + 1].split(","):
                    key, sep, value = component.partition("=")
                    if sep:
                        values[key] = value
                if values.get("type") == "bind" and values.get("source", "").startswith("/") and values.get("target", "").startswith("/"):
                    mounts.append((values["target"], Path(values["source"]).expanduser()))
    # /work is the standard writable run bind mount. Prefer an explicit mount
    # if the command recorded one.
    if not any(target == "/work" for target, _ in mounts):
        mounts.append(("/work", run_dir))
    return sorted(mounts, key=lambda item: len(item[0]), reverse=True)


def _map_container_path(value: str, mounts: list[tuple[str, Path]]) -> Path | None:
    raw = str(value or "").strip()
    if not raw.startswith("/"):
        return None
    for target, source in mounts:
        if raw == target or raw.startswith(target.rstrip("/") + "/"):
            suffix = raw[len(target):].lstrip("/")
            return source / suffix if suffix else source
    return None


def _resolve_host_or_container_path(
    value: str,
    *,
    run_dir: Path,
    mounts: list[tuple[str, Path]],
    allow_directory: bool = False,
    container_symlink_targets: list[str] | None = None,
) -> Path | None:
    """Resolve one file, translating absolute container symlink targets too."""
    raw = str(value or "").strip()
    if not raw:
        return None
    if is_portable_path(raw):
        resolved = resolve_stored_path(raw, run_dir=run_dir, must_exist=True)
        return Path(resolved) if resolved is not None else None
    lexical = Path(raw).expanduser()
    if not lexical.is_absolute():
        lexical = run_dir / lexical
    mapped = _map_container_path(str(lexical), mounts)
    if mapped is not None:
        lexical = mapped

    seen: set[str] = set()
    for _ in range(40):
        token = os.path.abspath(lexical)
        if token in seen:
            return None
        seen.add(token)
        try:
            if lexical.is_symlink():
                target = os.readlink(lexical)
                if Path(target).is_absolute():
                    if container_symlink_targets is not None and _is_container_path(target):
                        container_symlink_targets.append(target)
                    translated = _map_container_path(target, mounts)
                    lexical = translated if translated is not None else Path(target)
                else:
                    lexical = lexical.parent / target
                continue
            if lexical.is_file() or (allow_directory and lexical.is_dir()):
                return lexical.resolve(strict=True)
            return None
        except OSError:
            return None
    return None


def materialize_candidate_paths(run_dir: Path, candidate: dict[str, Any]) -> dict[str, Any]:
    """Persist candidate structures and retain the original path as provenance."""
    run_dir = Path(run_dir).expanduser().resolve()
    row = dict(candidate)
    raw_metadata = dict(row.get("raw_metadata") or {})
    original_paths = dict(raw_metadata.get("container_artifact_paths") or {})
    materialization_errors = list(raw_metadata.get("artifact_materialization_errors") or [])
    mounts = _container_mounts(run_dir)

    for field in _SCIENTIFIC_PATH_FIELDS:
        value = row.get(field)
        if not isinstance(value, str) or not value.strip():
            continue
        original = value.strip()
        parsed = urlsplit(original)
        is_container_path = _is_container_path(original)
        if is_container_path:
            original_paths.setdefault(field, original)
        elif is_portable_path(original):
            # Existing runs:/// and cluster:// references remain portable
            # references owned by the source job; do not duplicate them.
            continue
        elif parsed.scheme:
            message = f"{field} uses an unsupported artifact URI: {original}"
            if message not in materialization_errors:
                materialization_errors.append(message)
            continue

        symlink_container_targets: list[str] = []
        source = _resolve_host_or_container_path(
            original, run_dir=run_dir, mounts=mounts,
            container_symlink_targets=symlink_container_targets,
        )
        if symlink_container_targets:
            raw_metadata.setdefault("container_artifact_provenance", {})[field] = {
                "source_reference": original,
                "absolute_container_symlink_targets": symlink_container_targets,
            }
        if source is None:
            if is_container_path or Path(original).is_absolute():
                message = f"{field} output is inaccessible after execution: {original}"
                if message not in materialization_errors:
                    materialization_errors.append(message)
            continue

        try:
            source.relative_to(run_dir)
            source_is_local = True
        except ValueError:
            source_is_local = False
        if source_is_local:
            row[field] = portable_path(source, run_dir=run_dir)
            continue

        managed_reference = portable_path(source)
        if is_portable_path(managed_reference):
            row[field] = managed_reference
            raw_metadata.setdefault("materialized_source_paths", {})[field] = original
            continue

        digest = _file_digest(source)
        extension = source.suffix.lower() or ".dat"
        destination = run_dir / "artifacts" / "materialized_candidates" / field / f"{digest}{extension}"
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.is_file():
            shutil.copy2(source, destination)
        row[field] = portable_path(destination, run_dir=run_dir)
        raw_metadata.setdefault("materialized_source_paths", {})[field] = original

    if original_paths:
        raw_metadata["container_artifact_paths"] = original_paths
    if materialization_errors:
        raw_metadata["artifact_materialization_errors"] = materialization_errors
    row["raw_metadata"] = raw_metadata
    return row


def _candidate_references(run_dir: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    path = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    required: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    candidate_file = path
    if path.is_symlink():
        try:
            candidate_file = path.resolve(strict=True)
            candidate_file.relative_to(run_dir)
        except (OSError, ValueError):
            errors.append({"path": str(path.relative_to(run_dir)), "reason": "Normalized candidate collection symlink is missing or escapes its Results run"})
            return required, errors
    if not candidate_file.is_file():
        if path.is_symlink():
            errors.append({"path": str(path.relative_to(run_dir)), "reason": "Normalized candidate collection is a dangling symlink"})
        return required, errors
    mounts = _container_mounts(run_dir)
    for line_number, line in enumerate(candidate_file.read_text(errors="replace").splitlines(), 1):
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except ValueError:
            errors.append({"path": str(path.relative_to(run_dir)), "reason": f"Invalid candidate JSON at line {line_number}"})
            continue
        if not isinstance(row, dict):
            continue
        raw = row.get("raw_metadata") if isinstance(row.get("raw_metadata"), dict) else {}
        for message in raw.get("artifact_materialization_errors") or []:
            errors.append({"path": str(row.get("candidate_id") or "candidate"), "reason": str(message)})
        for field in _SCIENTIFIC_PATH_FIELDS:
            value = row.get(field)
            if not isinstance(value, str) or not value.strip():
                continue
            resolved = _resolve_host_or_container_path(value, run_dir=run_dir, mounts=mounts)
            if resolved is None:
                errors.append({"path": value, "reason": f"Required candidate {field} is inaccessible"})
            else:
                try:
                    relative = resolved.relative_to(run_dir).as_posix()
                except ValueError:
                    relative = value
                try:
                    size = resolved.stat().st_size
                    digest = _file_digest(resolved)
                except OSError:
                    errors.append({"path": value, "reason": f"Required candidate {field} became unreadable during validation"})
                    continue
                required.append({"path": relative, "candidate_id": str(row.get("candidate_id") or ""),
                                 "field": field, "size_bytes": str(size), "sha256": digest})
    return required, errors


def _dangling_links(run_dir: Path) -> list[dict[str, str]]:
    links: list[dict[str, str]] = []
    for root, dirs, files in os.walk(run_dir, followlinks=False):
        base = Path(root)
        for name in [*dirs, *files]:
            path = base / name
            try:
                if path.is_symlink() and not path.exists():
                    links.append({
                        "path": path.relative_to(run_dir).as_posix(),
                        "target": os.readlink(path),
                    })
            except OSError:
                continue
        # Do not recurse into directory symlinks.
        dirs[:] = [name for name in dirs if not (base / name).is_symlink()]
    return links


def validate_run_artifacts(
    run_dir: Path,
    result: dict[str, Any] | None = None,
    *,
    strict_dangling_links: bool = True,
) -> dict[str, Any]:
    """Check canonical structure references and declared artifact refs.

    New-job finalization is strict about links it cannot classify. Historical
    backup validation can set ``strict_dangling_links=False`` to inventory
    those unchanged links while still requiring canonical artifacts to resolve.
    This check never rewrites or dereferences a link's contents.
    """
    run_dir = Path(run_dir).expanduser().resolve()
    required, errors = _candidate_references(run_dir)
    result = result if isinstance(result, dict) else _read_json(run_dir / "result.json")
    mounts = _container_mounts(run_dir)
    downstream = result.get("downstream_artifacts") if isinstance(result.get("downstream_artifacts"), dict) else {}
    downstream_refs: list[dict[str, str]] = []

    def collect_downstream(value: Any, key: str) -> None:
        if isinstance(value, dict):
            for child_key, child in value.items():
                collect_downstream(child, str(child_key))
        elif isinstance(value, list):
            for child in value:
                collect_downstream(child, key)
        elif isinstance(value, str) and value.strip():
            resolved = _resolve_host_or_container_path(value, run_dir=run_dir, mounts=mounts, allow_directory=True)
            if resolved is None:
                errors.append({"path": value, "reason": f"Declared downstream artifact {key} is inaccessible"})
                return
            try:
                stored_path = resolved.relative_to(run_dir).as_posix()
            except ValueError:
                stored_path = value
            downstream_refs.append({"path": stored_path, "artifact": key})

    collect_downstream(downstream, "downstream artifact")
    declared: list[dict[str, str]] = []
    optional_paths: set[str] = set()
    manifest_path = run_dir / "artifacts.json"
    if manifest_path.is_symlink():
        try:
            manifest_target = manifest_path.resolve(strict=True)
            manifest_target.relative_to(run_dir)
        except (OSError, ValueError):
            errors.append({"path": "artifacts.json", "reason": "Artifact manifest symlink is missing or escapes its Results run"})
    if manifest_path.is_file() and not (manifest_path.is_symlink() and any(item["path"] == "artifacts.json" for item in errors)):
        try:
            payload = json.loads(manifest_path.read_text())
            if int(payload.get("schema_version") or 0) != 1 or str(payload.get("run_id") or "") != run_dir.name:
                raise ValueError("unsupported schema or run identity mismatch")
            for artifact in payload.get("artifacts") or []:
                if artifact.get("run_id") and str(artifact.get("run_id")) != run_dir.name:
                    errors.append({"path": str(artifact.get("path") or ""), "reason": "Artifact run identity does not match its manifest"})
                    continue
                relative = PurePosixPath(str(artifact.get("path") or ""))
                if not relative.parts or relative.is_absolute() or ".." in relative.parts:
                    errors.append({"path": str(artifact.get("path") or ""), "reason": "Artifact path is not a safe run-relative reference"})
                    continue
                artifact_metadata = artifact.get("metadata") if isinstance(artifact.get("metadata"), dict) else {}
                is_required = artifact.get("required", not bool(artifact_metadata.get("optional")))
                if not is_required:
                    optional_paths.add(relative.as_posix())
                    continue
                target = run_dir.joinpath(*relative.parts)
                try:
                    resolved = target.resolve(strict=True)
                    resolved.relative_to(run_dir)
                except (OSError, ValueError):
                    resolved = None
                if resolved is None or not resolved.is_file():
                    errors.append({"path": relative.as_posix(), "reason": f"Required {artifact.get('artifact_type', 'scientific')} artifact is inaccessible"})
                else:
                    size = artifact.get("size_bytes")
                    digest = str(artifact.get("sha256") or "")
                    if size is not None and resolved.stat().st_size != int(size):
                        errors.append({"path": relative.as_posix(), "reason": "Required artifact size differs from its manifest"})
                        continue
                    if digest and _file_digest(resolved) != digest:
                        errors.append({"path": relative.as_posix(), "reason": "Required artifact digest differs from its manifest"})
                        continue
                    declared.append({"path": relative.as_posix(), "artifact_type": str(artifact.get("artifact_type") or "")})
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
            errors.append({"path": "artifacts.json", "reason": f"Artifact manifest is invalid: {exc}"})

    required_paths = {item["path"] for item in [*required, *declared, *downstream_refs]}
    dangling = _dangling_links(run_dir)
    for item in dangling:
        if item["path"] in required_paths:
            item["classification"] = "required_scientific_artifact"
            continue
        if item["path"] in optional_paths:
            item["classification"] = "optional_execution_artifact"
            continue
        # Known caches, inputs and container-provided runtimes are execution
        # references, not retained scientific outputs. Preserve and report
        # them, but do not follow them during validation or backup.
        normalized = (item["path"] + " " + item["target"]).casefold()
        optional_markers = ("/input/", "/inputs/", "/cache/", "/.cache/", "/tmp/", "/temp/",
                            "/processed/", "/weights/", "/checkpoint/", "/model_params/",
                            "/.venv/", "/community_models/", "/assets/")
        if any(marker in normalized for marker in optional_markers):
            item["classification"] = "optional_execution_artifact"
            continue
        suffix = Path(item["path"]).suffix.casefold()
        in_output = any(part.casefold() in {"output", "outputs", "poses", "predictions"}
                        for part in Path(item["path"]).parts)
        if suffix in {".pdb", ".cif", ".mmcif", ".sdf", ".pdbqt", ".mol2"} and in_output:
            item["classification"] = "unresolved_scientific_output"
            if strict_dangling_links:
                errors.append({"path": item["path"], "reason": "Dangling scientific output is not declared or materialized as a canonical Results artifact"})
            continue
        item["classification"] = "unclassified_dangling_link"
        if strict_dangling_links:
            errors.append({"path": item["path"], "reason": "Dangling link cannot be classified as an optional execution artifact or a declared scientific output"})
    return {
        "status": "incomplete" if errors else "complete",
        "required_artifact_count": len(required) + len(declared) + len(downstream_refs),
        "required_missing": errors,
        "dangling_symlinks": dangling,
        "symlinks_followed": False,
    }
