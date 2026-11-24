"""Stage molecular structures and schedule TMol using the existing CUDA image."""
from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path
from typing import Literal

from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import run_docker_job
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path


WORKFLOW = "structure.tmol"
MODE_LABELS = {
    "score": "Score input geometry",
    "minimize": "Cartesian minimization",
    "fast_relax": "FastRelax",
}
FORMATS = {".pdb": "pdb", ".ent": "pdb", ".cif": "cif", ".mmcif": "cif", ".sdf": "sdf", ".mol2": "mol2"}
_RUNTIME_SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "tmol_structure.py"


def validate_inputs(parameters, *, connected=False):
    uploads = bool(parameters.get("structure_paths"))
    source = bool(parameters.get("source_run_dir"))
    collection = bool(parameters.get("candidates_jsonl"))
    if connected:
        if uploads or source or collection:
            raise ValueError("Choose either the connected collection or explicit structure inputs.")
    elif source != collection or uploads == source:
        raise ValueError("Choose uploaded molecular structures or one saved structure collection.")
    if parameters.get("mode", "fast_relax") not in MODE_LABELS:
        raise ValueError("Choose a supported TMol mode.")
    repeats = parameters.get("repeats", 1)
    if type(repeats) is not int or repeats < 1:
        raise ValueError("Choose at least one FastRelax repeat.")
    if str(parameters.get("gpu_device", "")).casefold() in {"none", "cpu", "off", "false"}:
        raise ValueError("TMol requires a CUDA GPU.")


def _structure_format(path):
    kind = FORMATS.get(path.suffix.lower())
    if kind in {"pdb", "cif"}:
        import gemmi

        if len(gemmi.read_structure(str(path))) != 1:
            raise ValueError("Choose a structure file containing one coordinate model.")
    elif kind == "sdf":
        if sum(line.strip() == "$$$$" for line in path.read_text().splitlines()) > 1:
            raise ValueError("Choose one molecule per SDF file.")
    elif kind == "mol2":
        if path.read_text().count("@<TRIPOS>MOLECULE") != 1:
            raise ValueError("Choose one molecule per MOL2 file.")
    else:
        raise ValueError("Use a PDB, mmCIF, SDF or MOL2 structure.")
    return kind


def run_tmol(
    structure_paths: list[Path] | None = None,
    source_run_dir: Path | None = None,
    candidates_jsonl: Path | None = None,
    mode: Literal["score", "minimize", "fast_relax"] = "fast_relax",
    prepare_ligands: bool = False,
    ligand_params_files: list[Path] | None = None,
    repeats: int = 1,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    job_name: str = "TMol structure analysis",
) -> Path:
    """Score, minimize, or FastRelax an explicit collection of fixed structures."""
    ligand_params_files = [Path(path) for path in ligand_params_files or []]
    if any(path.suffix.lower() != ".tmol" or not path.is_file() for path in ligand_params_files):
        raise ValueError("Choose existing .tmol ligand parameter files.")
    inputs = {
        "structure_paths": [str(path) for path in structure_paths or []],
        "ligand_params_files": [str(path) for path in ligand_params_files],
        "source_run_dir": str(source_run_dir) if source_run_dir else None,
        "candidates_jsonl": str(candidates_jsonl) if candidates_jsonl else None,
    }
    parameters = {"mode": mode, "prepare_ligands": prepare_ligands, "repeats": repeats,
                  "cpu_cores": cpu_cores, "gpu_device": gpu_device, "job_name": job_name}
    validate_inputs({**inputs, **parameters})
    if type(cpu_cores) is not int or cpu_cores < 1:
        raise ValueError("Choose at least one CPU core.")
    if candidates_jsonl:
        candidates = read_candidates(candidates_jsonl)
        if not candidates:
            raise ValueError("The selected collection contains no structures.")
        sources = [(candidate, resolve_stored_path(candidate.get("complex_pdb") or candidate.get("binder_pdb"),
                                                    run_dir=source_run_dir, must_exist=True))
                   for candidate in candidates]
    else:
        sources = [({}, Path(path)) for path in structure_paths or []]

    job = create_job("analysis", "tmol_structure_analysis", "tmol", inputs, parameters)
    output = job.artifacts_dir / "tmol"
    staged = job.artifacts_dir / "input"
    output.mkdir(parents=True, exist_ok=True)
    staged.mkdir(parents=True, exist_ok=True)
    items = []
    for index, (candidate, path) in enumerate(sources, 1):
        identifier = f"structure_{index:05d}"
        name = str(candidate.get("candidate_id") or (path.name if path else identifier))
        item = {"id": identifier, "name": name, "input_sha256": None}
        try:
            if not path or not path.is_file():
                raise ValueError("The selected structure is unavailable.")
            kind = _structure_format(path)
            destination = staged / f"{identifier}.{kind}"
            shutil.copyfile(path, destination)
            with destination.open("rb") as handle:
                item["input_sha256"] = hashlib.file_digest(handle, "sha256").hexdigest()
            item.update(format=kind, source=f"/work/artifacts/input/{destination.name}",
                        output=f"/work/artifacts/tmol/{identifier}.cif")
        except (OSError, RuntimeError, ValueError) as exc:
            item["input_error"] = str(exc)
        items.append(item)

    ligand_parameters = []
    for index, source in enumerate(ligand_params_files, start=1):
        destination = staged / f"ligand_parameters_{index:03d}.tmol"
        shutil.copyfile(source, destination)
        with destination.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        ligand_parameters.append({"name": source.name, "sha256": digest,
                                  "path": f"/work/artifacts/input/{destination.name}"})
    shutil.copy2(_RUNTIME_SCRIPT, staged / _RUNTIME_SCRIPT.name)
    write_json(staged / "tmol_request.json", {**parameters, "items": items,
               "ligand_parameters": ligand_parameters, "summary": "/work/artifacts/tmol/summary.json"})
    return_code = run_docker_job(
        job.run_dir, load_manifest("tmol"),
        [f"/work/artifacts/input/{_RUNTIME_SCRIPT.name}", "/work/artifacts/input/tmol_request.json"],
        gpu_device=gpu_device, use_host_user=True,
        environment={"OMP_NUM_THREADS": str(cpu_cores), "PYTHONUNBUFFERED": "1"},
    )
    summary = read_json(output / "summary.json")
    by_id = {row["id"]: row for row in summary.get("records", [])}
    rows, records = [], []
    for item, (source_candidate, _path) in zip(items, sources):
        record = by_id.get(item["id"]) or {**item, "status": "failed",
            "error": item.get("input_error") or f"TMol worker exited with code {return_code} before saving this item."}
        filename = f'{item["id"]}.cif'
        relative = f"artifacts/tmol/{filename}"
        if record.get("status") == "completed" and (output / filename).is_file():
            rows.append({
                "candidate_id": item["id"], "source_tool": "tmol", "stage": "structure_refinement",
                "complex_pdb": relative, "parents": [source_candidate["candidate_id"]] if source_candidate else [],
                "binder_sequence": source_candidate.get("binder_sequence"),
                "binder_chains": source_candidate.get("binder_chains", []),
                "target_chains": source_candidate.get("target_chains", []),
                "metrics": {key: record[key] for key in
                            ("tmol_energy_before", "tmol_energy_after", "tmol_energy_change")},
                "raw_metadata": {"name": item["name"], "input_sha256": item["input_sha256"],
                                 "tmol_mode": mode, "tmol_version": summary.get("tmol_version"),
                                 "source_candidate_id": source_candidate.get("candidate_id")},
            })
            record["structure"] = relative
        elif record.get("status") == "completed":
            record.update(status="failed", error="The worker did not publish the coordinate file.")
        records.append(record)
    summary.update(mode=mode, records=records, ligand_parameters=ligand_parameters)
    write_json(output / "summary.json", summary)
    columns = ["id", "name", "status", "tmol_energy_before", "tmol_energy_after", "tmol_energy_change",
               "structure", "ligand_parameterization", "input_waters", "output_waters", "input_sha256", "error"]
    with (output / "measurements.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    if rows:
        write_candidates(job.run_dir, "tmol", rows)
    finish_job(job.run_dir, bool(rows), {
        "outputs": {"summary": "artifacts/tmol/summary.json", "measurements": "artifacts/tmol/measurements.csv"},
        "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl"} if rows else {},
        "metrics": {"mode": mode, "completed_count": len(rows), "failed_count": len(records) - len(rows),
                    "tmol_version": summary.get("tmol_version")},
        **({"error": records[0]["error"]} if not rows and records else {}),
    })
    return job.run_dir
