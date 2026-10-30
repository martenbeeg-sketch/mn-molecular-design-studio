"""Explicit transfer and Database promotion for scientifically accepted designs."""
from __future__ import annotations

import fcntl
import hashlib
import io
import json
import zipfile
from copy import deepcopy
from pathlib import Path, PurePosixPath

MAX_ACCEPTED_BUNDLE_BYTES = 512 * 1024 * 1024
_STRUCTURES = {"complex_pdb", "binder_pdb", "validation_complex_pdb"}


def _candidate_rows(run_dir: Path) -> list[dict]:
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    return read_candidates(Path(run_dir))


def _resolve_structure(value: object, run_dir: Path, runs_root: Path) -> Path | None:
    if not isinstance(value, str) or not value.strip():
        return None
    path = Path(value)
    if not path.is_absolute():
        try:
            from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
            path = resolve_stored_path(value, run_dir=run_dir) or (run_dir / path)
        except (ValueError, OSError):
            path = run_dir / path
    try:
        path = path.resolve(strict=True)
        if not path.is_file() or not path.is_relative_to(runs_root.resolve()):
            return None
    except (OSError, RuntimeError):
        return None
    return path


def _worker_uri(value: object, node_id: str, runs_root: Path) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        relative = Path(value).expanduser().resolve().relative_to(runs_root.resolve())
    except (OSError, ValueError):
        return None
    if len(relative.parts) >= 2:
        # Plain provenance token, deliberately not a filesystem URI. Candidate
        # path resolvers must never try to open it on the importing workstation.
        return f"{node_id}/" + relative.as_posix()
    return None


def build_accepted_design_bundle(run_dir: Path, *, node_id: str, validation_job: str) -> tuple[bytes, dict]:
    """Package accepted rows and their structures only; rejected rows never leave the worker."""
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root

    run_dir = Path(run_dir).resolve()
    runs = Path(runs_root()).resolve()
    source_rows = _candidate_rows(run_dir)
    identities = [str(row.get("candidate_id") or "") for row in source_rows]
    if len(set(identities)) != len(source_rows):
        raise ValueError("The validation run contains duplicate candidate identities.")
    accepted = []
    for row in source_rows:
        raw = row.get("raw_metadata") if isinstance(row.get("raw_metadata"), dict) else {}
        decision = raw.get("validation_acceptance") if isinstance(raw.get("validation_acceptance"), dict) else {}
        if decision.get("status") != "accepted":
            continue
        source = raw.get("source_candidate")
        promoted = deepcopy(source) if isinstance(source, dict) else deepcopy(row)
        promoted_raw = promoted.get("raw_metadata") if isinstance(promoted.get("raw_metadata"), dict) else {}
        promoted_raw = dict(promoted_raw)
        promoted_raw.pop("source_candidate", None)
        promoted_raw["validation_acceptance"] = deepcopy(decision)
        promoted_raw["validation_measurements"] = deepcopy(row.get("metrics") or {})
        source_run = raw.get("source_run_dir") or (raw.get("validation_reference") or {}).get("source_run_dir")
        promoted_raw["validation_provenance"] = {
            "validation_job_id": validation_job.removeprefix("cluster:"),
            "validation_candidate_id": str(row.get("candidate_id") or ""),
            "validation_engine": str((row.get("metrics") or {}).get("engine") or ""),
            "worker_reference": _worker_uri(source_run, node_id, runs) or str(source_run or ""),
        }
        promoted["raw_metadata"] = promoted_raw
        promoted_metrics = dict(promoted.get("metrics") or {})
        for metric, value in (row.get("metrics") or {}).items():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                promoted_metrics[f"validation_{metric}"] = value
        promoted["metrics"] = promoted_metrics
        validation_path = next((row.get(key) for key in ("complex_pdb", "binder_pdb") if row.get(key)), None)
        if validation_path:
            promoted_raw["validation_complex_pdb"] = validation_path
        if not promoted.get("candidate_id"):
            promoted["candidate_id"] = str(row.get("candidate_id") or "")
        accepted.append(promoted)
    if not accepted:
        raise ValueError("This validation run has no explicitly accepted designs to transfer.")

    archive = io.BytesIO()
    files: dict[Path, str] = {}
    artifact_manifest = []

    def rewrite_paths(value, key=""):
        if isinstance(value, dict):
            return {name: rewrite_paths(child, str(name)) for name, child in value.items()}
        if isinstance(value, list):
            return [rewrite_paths(child, key) for child in value]
        if key in _STRUCTURES:
            path = _resolve_structure(value, run_dir, runs)
            if path is None:
                return None
            relative = files.get(path)
            if relative is None:
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                suffix = path.suffix.lower()
                if suffix not in {".pdb", ".cif", ".mmcif", ".gz"}:
                    raise ValueError(f"Unsupported accepted structure format: {path.name}")
                relative = f"artifacts/accepted/{digest[:20]}{suffix}"
                files[path] = relative
                artifact_manifest.append({"path": relative, "sha256": digest, "size_bytes": path.stat().st_size})
            return relative
        if key == "target_pdb":
            return None
        if key == "source_run_dir":
            # Kept as a cluster-qualified provenance value above, but it must
            # not remain an operational path that path resolvers try to open.
            return None
        if isinstance(value, str):
            from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import is_operational_path_key
            if is_operational_path_key(key) and value.startswith(("/", "runs://", "cluster://", "reference://", "app://")):
                return None
        return value

    packaged = [rewrite_paths(candidate) for candidate in accepted]
    if any(not (candidate.get("complex_pdb") or candidate.get("binder_pdb")
                or (candidate.get("raw_metadata") or {}).get("validation_complex_pdb"))
           for candidate in packaged):
        raise ValueError("An accepted design has no transferable structure file.")
    manifest = {
        "schema_version": 1,
        "validation_job": validation_job,
        "validation_run_id": run_dir.name,
        "candidate_count": len(packaged),
        "accepted_candidate_ids": [str(row.get("candidate_id") or "") for row in packaged],
        "artifacts": artifact_manifest,
    }
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as output:
        output.writestr("manifest.json", json.dumps(manifest, sort_keys=True))
        output.writestr("candidates.jsonl", "".join(json.dumps(row, sort_keys=True) + "\n" for row in packaged))
        for path, relative in files.items():
            if path.stat().st_size > MAX_ACCEPTED_BUNDLE_BYTES:
                raise ValueError("An accepted structure exceeds the 512 MB transfer limit.")
            output.write(path, relative)
    payload = archive.getvalue()
    if len(payload) > MAX_ACCEPTED_BUNDLE_BYTES:
        raise ValueError("Accepted structures exceed the 512 MB transfer limit.")
    return payload, manifest


def _read_bundle(payload: bytes) -> tuple[dict, list[dict], dict[str, bytes]]:
    if not isinstance(payload, bytes) or len(payload) > MAX_ACCEPTED_BUNDLE_BYTES:
        raise ValueError("Accepted design bundle is invalid or exceeds the transfer limit.")
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        members = archive.infolist()
        total = sum(item.file_size for item in members)
        if total > MAX_ACCEPTED_BUNDLE_BYTES or len(members) > 4096:
            raise ValueError("Accepted design bundle exceeds the extraction limit.")
        contents = {}
        for item in members:
            name = PurePosixPath(item.filename)
            if name.is_absolute() or ".." in name.parts or item.is_dir():
                raise ValueError("Accepted design bundle contains an unsafe path.")
            if item.filename not in {"manifest.json", "candidates.jsonl"} and not item.filename.startswith("artifacts/accepted/"):
                raise ValueError("Accepted design bundle contains an unexpected file.")
            contents[item.filename] = archive.read(item)
    manifest = json.loads(contents.pop("manifest.json"))
    if not isinstance(manifest, dict) or manifest.get("schema_version") != 1:
        raise ValueError("Unsupported accepted design bundle schema.")
    rows = [json.loads(line) for line in contents.pop("candidates.jsonl").decode("utf-8").splitlines() if line.strip()]
    row_ids = [str(row.get("candidate_id") or "") for row in rows if isinstance(row, dict)]
    if (not rows or len(row_ids) != len(rows) or len(rows) != manifest.get("candidate_count")
            or len(set(row_ids)) != len(rows) or set(row_ids) != set(manifest.get("accepted_candidate_ids") or [])):
        raise ValueError("Accepted design bundle candidate identities do not match its manifest.")
    if any(((row.get("raw_metadata") or {}).get("validation_acceptance") or {}).get("status") != "accepted" for row in rows):
        raise ValueError("The transfer includes a design without an explicit accepted decision.")
    declared_paths = set()
    for item in manifest.get("artifacts", []):
        data = contents.get(item.get("path"))
        if (data is None or len(data) != item.get("size_bytes")
                or hashlib.sha256(data).hexdigest() != item.get("sha256")):
            raise ValueError("An accepted structure failed its SHA-256 integrity check.")
        declared_paths.add(str(item["path"]))
    if set(contents) != declared_paths:
        raise ValueError("Accepted design bundle contains unmanifested structure data.")
    for row in rows:
        raw = row.get("raw_metadata") if isinstance(row.get("raw_metadata"), dict) else {}
        structure_paths = {row.get("complex_pdb"), row.get("binder_pdb"), raw.get("validation_complex_pdb")}
        structure_paths.discard(None)
        if not structure_paths or not structure_paths <= declared_paths:
            raise ValueError("An accepted design is missing its declared structure data.")
    return manifest, rows, contents


def promote_accepted_design_bundle(payload: bytes, *, user: dict, project_id: str) -> dict:
    """Create one idempotent local Results run and explicitly ingest accepted designs."""
    from mn_molecular_design_studio.collection_policy import set_run_policy
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, write_json
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root
    from mn_molecular_design_studio.workspace import mark_studio_job

    manifest, rows, contents = _read_bundle(payload)
    source_job = str(manifest.get("validation_job") or "")
    key_text = json.dumps({"source_job": source_job, "candidate_ids": sorted(manifest["accepted_candidate_ids"])}, sort_keys=True)
    import_key = hashlib.sha256(key_text.encode()).hexdigest()
    group_root = Path(runs_root()) / "candidate-import"
    registry_root = group_root / ".accepted-promotions"
    registry_root.mkdir(parents=True, exist_ok=True)
    lock_path = registry_root / "lock"
    with lock_path.open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        registry_path = registry_root / f"{import_key}.json"
        previous = None
        try:
            previous = json.loads(registry_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            previous = None
        if isinstance(previous, dict) and previous.get("run_dir"):
            prior_run = Path(previous["run_dir"])
            if read_json(prior_run / "metadata.json").get("status") == "completed":
                return {"run_dir": str(prior_run), "candidate_count": len(rows), "already_imported": True}
        # Recover if the Database ingestion completed before the registry write.
        if group_root.is_dir():
            for candidate_run in group_root.iterdir():
                if not candidate_run.is_dir() or candidate_run.name.startswith("."):
                    continue
                metadata = read_json(candidate_run / "metadata.json")
                if metadata.get("accepted_import_key") == import_key and metadata.get("status") == "completed":
                    write_json(registry_path, {"run_dir": str(candidate_run), "candidate_count": len(rows)})
                    return {"run_dir": str(candidate_run), "candidate_count": len(rows), "already_imported": True}

        job = create_job(
            "candidate-import", "accepted_design_import", "accepted_design_import",
            inputs={"validation_job": source_job, "accepted_candidate_ids": manifest["accepted_candidate_ids"]},
            params={"project_id": project_id, "accepted_import_key": import_key},
        )
        metadata = read_json(job.run_dir / "metadata.json")
        metadata.update(accepted_import_key=import_key, project_id=project_id, owner_id=str(user.get("id") or ""))
        write_json(job.run_dir / "metadata.json", metadata)
        mark_studio_job({"job": {"run_dir": str(job.run_dir)}}, project_id=project_id, owner_id=str(user.get("id") or ""))
        for relative, data in contents.items():
            target = (job.run_dir / relative).resolve()
            if not target.is_relative_to(job.run_dir.resolve()):
                raise ValueError("Accepted artifact path escapes the import run.")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        candidates = write_candidates(job.run_dir, "accepted_design_import", rows)
        set_run_policy(job.run_dir, "automatic")
        finish_job(job.run_dir, True, {
            "outputs": {"candidate_count": len(candidates), "source_validation_job": source_job,
                        "accepted_candidate_ids": manifest["accepted_candidate_ids"]},
            "metrics": {"candidate_count": len(candidates), "accepted_candidate_count": len(candidates),
                        "explicit_researcher_promotion": True},
            "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl"},
        })
        write_json(registry_path, {"run_dir": str(job.run_dir), "candidate_count": len(candidates)})
        return {"run_dir": str(job.run_dir), "candidate_count": len(candidates), "already_imported": False}
