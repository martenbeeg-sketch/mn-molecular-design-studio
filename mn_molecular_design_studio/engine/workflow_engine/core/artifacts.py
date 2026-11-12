from __future__ import annotations

from io import BytesIO
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import portable_path


@dataclass(frozen=True)
class Artifact:
    name: str
    path: Path
    type: str
    description: str = ""

    def to_json(self, run_dir: Path) -> dict:
        return {
            "name": self.name,
            "path": portable_path(self.path, run_dir=run_dir),
            "type": self.type,
            "description": self.description,
        }


def artifact_path(run_dir: Path, *parts: str) -> Path:
    path = run_dir / "artifacts" / Path(*parts)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def build_run_zip(run_dir: Path, *, user: dict | None = None) -> bytes:
    """Create a fresh in-memory zip for a run without writing it to disk."""
    from mn_molecular_design_studio.artifact_integrity import retained_artifact_file

    run_dir = Path(run_dir).resolve(strict=True)
    if user is not None:
        from mn_molecular_design_studio.auth import can_access_job
        from mn_molecular_design_studio.job_catalog import get_job_access_record
        from mn_molecular_design_studio.runtime import runs_root
        relative = run_dir.relative_to(Path(runs_root()).resolve())
        access = get_job_access_record(*relative.parts) if len(relative.parts) == 2 else None
        if not access or not can_access_job(user, access):
            raise PermissionError("The source Results job is unavailable or you do not have access to it.")
    contract_names = ["input.json", "metadata.json", "command.json", "stdout.log", "stderr.log", "result.json",
                      "artifacts.json", "graph_execution.json", "graph_snapshot.json", "native_execution.json"]
    skipped_links = []
    paths = [run_dir / name for name in contract_names]
    for directory, folders, files in os.walk(run_dir / "artifacts", followlinks=False):
        for folder in list(folders):
            path = Path(directory) / folder
            if path.is_symlink():
                skipped_links.append(path.relative_to(run_dir).as_posix())
                folders.remove(folder)
        paths.extend(Path(directory) / name for name in sorted(files))
    buffer = BytesIO()
    with ZipFile(buffer, "w", compression=ZIP_DEFLATED) as archive:
        for path in paths:
            if path.is_symlink() and (not path.exists() or not path.resolve().is_relative_to(run_dir)):
                skipped_links.append(path.relative_to(run_dir).as_posix())
                continue
            if path.is_file():
                verified = retained_artifact_file(run_dir, path)
                archive.writestr(path.relative_to(run_dir).as_posix(), verified.read_bytes())
        if user is not None:
            _write_scientific_bundle(archive, run_dir, user)
        archive.writestr("export_manifest.json", json.dumps({"schema": "studio.run_export.v1",
            "skipped_symlinks": sorted(skipped_links), "scientific_bundle": user is not None}, indent=2))
    return buffer.getvalue()


def _write_scientific_bundle(archive: ZipFile, run_dir: Path, user: dict) -> None:
    from mn_molecular_design_studio.artifact_integrity import authorize_result_artifact, retained_artifact_file
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    from mn_molecular_design_studio.runtime import runs_root

    root = Path(runs_root()).resolve()
    candidate_file = run_dir / "artifacts/normalized_candidates/candidates.jsonl"
    if not candidate_file.exists():
        return
    if not candidate_file.resolve().is_relative_to(run_dir):
        raise ValueError("Candidate collection escapes its Results run.")
    rows, structures = [], {}
    for line in retained_artifact_file(run_dir, candidate_file).read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        references = {}
        for field in ("complex_pdb", "binder_pdb", "target_pdb"):
            original = row.get(field)
            if not original:
                continue
            path = resolve_stored_path(original, run_dir=run_dir, must_exist=True)
            if path is None or not path.is_file():
                raise ValueError(f"Required candidate {field} is unavailable for export.")
            path = path.resolve(strict=True)
            group, run_id, _ = authorize_result_artifact(user, path, source_run_dir=run_dir, runs_root=root)
            contents = retained_artifact_file(root / group / run_id, path).read_bytes()
            digest = hashlib.sha256(contents).hexdigest()
            name = f"scientific_bundle/structures/{digest}{path.suffix.lower()}"
            if name not in structures:
                archive.writestr(name, contents)
                structures[name] = {"path": name, "sha256": digest, "size_bytes": len(contents)}
            references[field] = original
            row[field] = name
        row["export_source_references"] = references
        rows.append(row)
    archive.writestr("scientific_bundle/candidates.jsonl", "".join(json.dumps(row) + "\n" for row in rows))
    archive.writestr("scientific_bundle/manifest.json", json.dumps({"schema": "studio.scientific_bundle.v1",
        "paths_relative_to": "archive_root", "candidate_count": len(rows), "structures": list(structures.values()),
        "scope": "final candidate coordinates and measurements; not all intermediate child-run files"}, indent=2))
