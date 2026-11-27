"""Recover preparation steps from exact MN-Ligand archive artifacts."""
from __future__ import annotations

import hashlib
import json
import os
from collections import defaultdict
from pathlib import Path

from mn_molecular_design_studio.engine.workflow_engine.vendor.ligand_runtime.mn_ligand.core.jobs import JobRecord
from mn_molecular_design_studio.engine.workflow_engine.vendor.ligand_runtime.mn_ligand.core.provenance import (
    LINEAGE_METADATA_KEYS, modification_history, ordered_job_lineage,
)
from mn_molecular_design_studio.preparation_history import COORDINATE_SUFFIXES, validate_history

_GROUPS = ("protein-cleaning", "terminal-repair", "target-trimming", "target-orientation",
           "structure-jobs", "prepared-structures", "selected-complexes", "energy-minimization")


def archive_roots():
    root = os.environ.get("MDS_MN_LIGAND_RUNS", "").strip()
    if root:
        return (Path(root).expanduser(),)
    return tuple(path for base in ("/mnt/data/RESULTS", "/data/RESULTS") if
                 (path := Path(base) / "mn-ligand-workdir/workdir/runs").is_dir())


def _load_lineage(root, run_dir, jobs):
    run_id = run_dir.name
    if run_id in jobs:
        return jobs[run_id]
    job = JobRecord.load(run_dir, load_artifacts=False)
    if not job.metadata:
        return None
    jobs[run_id] = job
    identifiers = [job.parent_run_id, job.workflow_parent_run_id, *(job.metadata.get(key) for key in LINEAGE_METADATA_KEYS)]
    for value in identifiers:
        if isinstance(value, str) and Path(value).name == value and value and value not in jobs:
            candidates = sorted(root.glob(f"*/{value}/metadata.json")) if all(c.isalnum() or c in "-_" for c in value) else []
            if len(candidates) == 1 and candidates[0].resolve().is_relative_to(root):
                _load_lineage(root, candidates[0].parent, jobs)
    return job


def _steps(job, jobs):
    steps = []
    changes = defaultdict(list)
    for item in modification_history(job, jobs):
        changes[item["run_id"]].append(item)
    for ancestor in ordered_job_lineage(job, jobs):
        if ancestor.status != "completed":
            continue
        if ancestor.metadata.get("pdb_id") and not ancestor.parent_run_id:
            steps.append({"kind": "pdb_import", "label": "PDB import", "run_id": ancestor.run_id,
                          "summary": str(ancestor.metadata["pdb_id"])})
        steps.extend(changes.get(ancestor.run_id, []))
        if ancestor.job_type == "target_orientation":
            steps.append({"kind": "target_orientation", "label": "Target orientation",
                          "run_id": ancestor.run_id, "summary": ancestor.tool or "Rigid coordinate transform"})
    return steps


def recover_histories(digests_by_size, roots):
    """Only associate completed production artifacts with identical file bytes."""
    wanted = {size: set(digests) for size, digests in digests_by_size.items()}
    found = defaultdict(list)
    if not any(wanted.values()):
        return {}
    for root in roots:
        root = Path(root).resolve()
        jobs = {}
        for group in _GROUPS:
            for path in sorted((root / group).rglob("*")):
                if path.suffix.lower() not in COORDINATE_SUFFIXES or not path.is_file():
                    continue
                relative = path.relative_to(root)
                if len(relative.parts) < 3 or set(relative.parts[2:-1]) & {"input", "inputs", "tmp"}:
                    continue
                if any(tag in path.stem.casefold() for tag in ("_raw", "source_", "template", "ligand_ref")):
                    continue
                try:
                    if not path.resolve().is_relative_to(root) or (size := path.stat().st_size) not in wanted:
                        continue
                    with path.open("rb") as handle:
                        digest = hashlib.file_digest(handle, "sha256").hexdigest()
                except OSError:
                    continue
                if digest not in wanted[size]:
                    continue
                job = _load_lineage(root, root.joinpath(*relative.parts[:2]), jobs)
                if job is None or job.status != "completed":
                    continue
                steps = _steps(job, jobs)
                if any(step["kind"] != "pdb_import" for step in steps):
                    found[digest].append((steps, {"run_id": job.run_id, "path": relative.as_posix(), "sha256": digest}))
    histories = {}
    for digest, matches in found.items():
        signatures = {
            json.dumps([{key: step[key] for key in ("kind", "label", "summary", "details") if key in step}
                        for step in steps], sort_keys=True)
            for steps, _ in matches
        }
        if len(signatures) == 1:
            histories[digest] = validate_history({
                "version": 1, "structure_sha256": digest, "application": "MN-Ligand",
                "verification": "archive_sha256", "steps": matches[0][0],
                "artifacts": [artifact for _, artifact in matches],
            }, digest)
    return histories
