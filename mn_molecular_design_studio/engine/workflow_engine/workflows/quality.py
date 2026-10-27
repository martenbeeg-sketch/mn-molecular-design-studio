"""One scheduled batch per collection for QC and diversity analysis."""
from __future__ import annotations

import copy
import hashlib
import os
from pathlib import Path
import shutil
import subprocess

import pandas as pd

from mn_molecular_design_studio.quality_tools import TOOLS
from mn_molecular_design_studio.measurement_status import measurement_outcome, measurement_summary
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import finish_job, read_json, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import apply_docker_cpu_limits_to_steps
from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _analysis_job, _allocated_cpu_cores, _resolve_path


def _digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def docker_user_args():
    """Rootless Docker's container root already maps to the invoking host user."""
    info = subprocess.run(["docker", "info", "--format", "{{json .SecurityOptions}}"], capture_output=True, text=True, check=False)
    return [] if "rootless" in info.stdout else ["--user", f"{os.getuid()}:{os.getgid()}"]


def run_quality(source_run_dir, candidates_jsonl, *, tool, settings, cpu_cores=4, gpu_device=None):
    if tool not in TOOLS:
        raise ValueError("Unknown QC tool.")
    for field in ("ph", "contact_distance", "similarity_threshold", "coverage"):
        if field in settings:
            value = float(settings[field])
            low, high = (0., 14.) if field == "ph" else (.1, 20.) if field == "contact_distance" else (.01, 1.)
            if not low <= value <= high:
                raise ValueError(f"{field} must be between {low} and {high}.")
    source_run_dir = Path(source_run_dir)
    profile = TOOLS[tool]
    if candidates_jsonl is None and profile["kind"] == "ligand":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_candidates
        candidates = docking_candidates(source_run_dir)
        candidates_jsonl = source_run_dir / "artifacts/normalized_candidates/candidates.jsonl"
    else:
        candidates = read_candidates(Path(candidates_jsonl))
    if not candidates:
        raise ValueError("Select a nonempty collection of candidates.")
    if len(candidates) > 100000:
        raise ValueError("QC batches support at most 100,000 candidates.")
    workers = _allocated_cpu_cores(cpu_cores)
    settings = {**settings, "workers": workers}
    references = []
    if profile.get("references"):
        from mn_molecular_design_studio.runtime import reference_root
        from mn_molecular_design_studio.evaluations import file_digest
        for relative in profile["references"]:
            path = reference_root() / relative
            if not path.is_file():
                raise ValueError("Required analysis checkpoint is missing: " + relative)
            references.append(path)
        settings["checkpoint_sha256"] = {path.name: file_digest(path) for path in references}
    module = Path(__file__).parent
    runtime_files = [module / "quality.py", module / "quality_runner.py", module / "quality_ligands.py", *sorted((module / "qc_ovo").glob("*.py")), *sorted((module / "qc_ovo").glob("*.json"))]
    protocol = hashlib.sha256()
    for path in runtime_files:
        protocol.update(path.name.encode())
        protocol.update(path.read_bytes())
    settings["protocol_revision"] = protocol.hexdigest()
    job = _analysis_job(source_run_dir, Path(candidates_jsonl), tool=tool, params={**settings, "cpu_cores": workers},
                        note="One collection, one batch job; per-item progress is in artifacts/analysis/batch_status.json.", mode="docker-batch")
    run = job.run_dir
    output = run / "artifacts/analysis"
    output.mkdir(parents=True, exist_ok=True)
    items, enhanced = [], copy.deepcopy(candidates)
    for index, candidate in enumerate(enhanced):
        identifier = f"item_{index:06d}"
        raw = candidate.get("raw_metadata") or {}
        item = {"id": identifier, "candidate_id": str(candidate.get("candidate_id") or identifier),
                "sequence": candidate.get("binder_sequence") or "",
                "binder_chains": candidate.get("binder_chains") or [], "target_chains": candidate.get("target_chains") or []}
        for key in ("complex_pdb", "binder_pdb", "target_pdb"):
            path = _resolve_path(source_run_dir, str(candidate.get(key) or ""))
            if path and path.is_file():
                staged = run / "artifacts/input" / identifier / (key + path.suffix)
                staged.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, staged)
                candidate[key] = str(staged)
            else:
                candidate[key] = None
        path = candidate.get("complex_pdb") or candidate.get("binder_pdb")
        if path and Path(path).is_file():
            item["structure"] = "/work/" + str(Path(path).relative_to(run))
            item["structure_digest"] = _digest(path)
        elif profile["kind"] != "sequence":
            item["input_error"] = "The selected candidate has no accessible structure."
        # Do not assume equal residue numbers across unrelated structures.
        target_path = candidate.get("target_pdb")
        item["target_key"] = _digest(target_path) if target_path and Path(target_path).is_file() else raw.get("target_structure_digest", "")
        if profile["kind"] == "interface" and (not item["binder_chains"] or not item["target_chains"]):
            item["input_error"] = "Declare binder and target chain groups before interface analysis."
        if profile["kind"] == "ligand" and not item.get("input_error"):
            from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import prepare
            try:
                prepare(item, candidate, run)
                item["structure_digest"] = _digest(candidate["complex_pdb"])
                item["target_key"] = _digest(candidate["target_pdb"])
            except (OSError, ValueError, RuntimeError, StopIteration) as exc:
                item["input_error"] = str(exc)
        candidate["raw_metadata"] = {**(candidate.get("raw_metadata") or raw), "source_candidate": copy.deepcopy(candidates[index]),
                                     "analysis_tool": tool, "analysis_subject": "ligand" if profile["kind"] == "ligand" else "protein", "result_kind": "app_reanalysis"}
        items.append(item)
    manifest = {"schema": "studio.qc.v1", "tool": tool, "items": items, "settings": settings, "output": "/work/artifacts/analysis"}
    write_json(output / "batch_manifest.json", manifest)
    runtime = run / "artifacts/runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(module / "quality_runner.py", runtime / "quality_runner.py")
    shutil.copyfile(module / "quality_ligands.py", runtime / "quality_ligands.py")
    shutil.copyfile(Path(__file__), runtime / "quality.py")
    shutil.copytree(module / "qc_ovo", runtime / "qc_ovo", dirs_exist_ok=True, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    gpu_args = []
    reference_args = []
    if profile.get("gpu"):
        from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args
        gpu_args = docker_gpu_args(gpu_device or "0")
        from mn_molecular_design_studio.runtime import reference_root
        reference_args = ["-v", f"{reference_root()}/esm_models:/reference/esm_models:ro", "-e", "TORCH_HOME=/reference"]
    command = ["docker", "run", "--rm", "--network", "none", "--entrypoint", "python3", *docker_user_args(), *gpu_args, *reference_args,
               "-e", "OMP_NUM_THREADS=1", "-e", "OPENBLAS_NUM_THREADS=1", "-v", f"{run}:/work", "-v", f"{runtime}:/opt/studio-qc:ro",
               "-w", "/work/artifacts/analysis", profile["image"], "/opt/studio-qc/quality_runner.py", "--manifest", "/work/artifacts/analysis/batch_manifest.json"]
    command = apply_docker_cpu_limits_to_steps(run, [{"name": tool + "-batch", "command": command}])[0]["command"]
    inspected = subprocess.run(["docker", "image", "inspect", profile["image"], "--format", "{{.Id}}"], capture_output=True, text=True, check=False)
    image_digest = inspected.stdout.strip() if inspected.returncode == 0 else ""
    if image_digest:
        command[command.index(profile["image"])] = image_digest
    write_json(run / "command.json", {"mode": "docker-batch", "command": command, "schema": manifest["schema"],
               "cpu_workers": workers, "candidate_count": len(items), "image": profile["image"], "image_digest": image_digest})
    with (run / "stdout.log").open("a") as stdout, (run / "stderr.log").open("a") as stderr:
        try:
            return_code = subprocess.run(command, stdout=stdout, stderr=stderr, check=False).returncode
        except OSError as exc:
            stderr.write(str(exc))
            return_code = 127
    status = read_json(output / "batch_status.json")
    results = {item["id"]: item for item in status.get("items", [])}
    table = []
    outcomes = []
    for source, item in zip(enhanced, items):
        result = results.get(item["id"], {"status": "failed", "metrics": {}, "error": f"Container exited {return_code} without an item result; see Logs."})
        measured = dict(result.get("metrics") or {})
        outcome = measurement_outcome(measured, profile["metrics"],
            result.get("error") or ("QC item did not complete: " + str(result["status"]) if result["status"] != "completed" else ""),
            boolean_fields=("posebusters_passed_all",))
        outcomes.append(outcome)
        measured.update(qc_status=result["status"], qc_error=outcome["error"], analysis_status=outcome["status"])
        source.update(candidate_id=item["candidate_id"] + "_" + tool, source_tool=tool, tool=tool, stage="analysis", metrics=measured,
                      parents=[item["candidate_id"]])
        source["raw_metadata"].update(analysis_parameters=settings, image_digest=image_digest, measurement_status=outcome,
                                      clustering_run_id=str(run) if profile.get("cohort") else "",
                                      structure_digest=item.get("structure_digest", ""), target_structure_digest=item["target_key"])
        table.append({"candidate_id": source["candidate_id"], "name": source["raw_metadata"].get("compound_name") or source["raw_metadata"].get("name") or item["candidate_id"],
                      "compound_id": source["raw_metadata"].get("compound_id") or "", "smiles": item.get("smiles") or "",
                      "source_candidate_id": item["candidate_id"], "status": result["status"], **measured})
    pd.DataFrame(table).to_csv(output / "qc_metrics.csv", index=False)
    normalized = write_candidates(run, tool, enhanced)
    counts = measurement_summary(outcomes)
    write_json(output / "qc_summary.json", {"schema": "studio.qc.v1", "tool": tool, **counts,
               "completed_count": counts["scored_count"], "failed_count": len(items) - counts["scored_count"], "settings": settings,
               "image_digest": image_digest, "representative_policy": "first source item in each cluster; native ranks are unchanged"})
    artifacts = [{"name": path.name, "path": str(path.relative_to(run)), "type": "analysis"} for path in output.rglob("*") if path.is_file()]
    finish_job(run, counts["scored_count"] > 0, {"outputs": {"artifacts": artifacts, "candidates": normalized},
        "execution_status": "completed" if return_code == 0 else "failed", "scientific_status": counts["scientific_status"],
        "metrics": counts,
        "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl"}})
    return run


def run_sequence_qc(source_run_dir: Path, candidates_jsonl: Path, ph: float = 7.4, cpu_cores: int = 4) -> Path:
    """Sequence composition, charge at the declared pH, pI, entropy, GRAVY and aromaticity."""
    return run_quality(source_run_dir, candidates_jsonl, tool="sequence_qc", settings={"ph": ph}, cpu_cores=cpu_cores)


def run_proteinsol(source_run_dir: Path, candidates_jsonl: Path, cpu_cores: int = 4) -> Path:
    """Protein-Sol sequence solubility prediction; minimum sequence length 21."""
    return run_quality(source_run_dir, candidates_jsonl, tool="proteinsol", settings={}, cpu_cores=cpu_cores)


def run_biotite_interface(source_run_dir: Path, candidates_jsonl: Path, contact_distance: float = 5.0, cpu_cores: int = 4) -> Path:
    """Heavy-atom contacts, residue maps, geometric salt bridges and explicit-H hydrogen bonds."""
    return run_quality(source_run_dir, candidates_jsonl, tool="biotite_interface", settings={"contact_distance": contact_distance}, cpu_cores=cpu_cores)


def run_sequence_clustering(source_run_dir: Path, candidates_jsonl: Path, similarity_threshold: float = .8, cpu_cores: int = 4) -> Path:
    """Single-linkage clusters from global aligned identity divided by the longer sequence length."""
    return run_quality(source_run_dir, candidates_jsonl, tool="sequence_clustering", settings={"similarity_threshold": similarity_threshold}, cpu_cores=cpu_cores)


def run_interface_clustering(source_run_dir: Path, candidates_jsonl: Path, similarity_threshold: float = .5, contact_distance: float = 5., cpu_cores: int = 4) -> Path:
    """Cluster identical-target binding epitopes using Jaccard similarity of contacted target residues."""
    return run_quality(source_run_dir, candidates_jsonl, tool="interface_clustering", settings={"similarity_threshold": similarity_threshold, "contact_distance": contact_distance}, cpu_cores=cpu_cores)


def run_foldseek(source_run_dir: Path, candidates_jsonl: Path, similarity_threshold: float = .5, coverage: float = .8, cpu_cores: int = 4) -> Path:
    """Foldseek 10-941cd33 structure clustering on the declared binder chains."""
    return run_quality(source_run_dir, candidates_jsonl, tool="foldseek", settings={"similarity_threshold": similarity_threshold, "coverage": coverage}, cpu_cores=cpu_cores)


def run_pep_patch(source_run_dir: Path, candidates_jsonl: Path, ph: float = 7.4, cpu_cores: int = 4) -> Path:
    """Hydrophobic and electrostatic surface patches on the declared binder chains."""
    return run_quality(source_run_dir, candidates_jsonl, tool="pep_patch", settings={"ph": ph}, cpu_cores=cpu_cores)


def run_esm1v_score(source_run_dir: Path, candidates_jsonl: Path, cpu_cores: int = 4, gpu_device: str = "0") -> Path:
    """Five-model ensemble mean native amino-acid probability (unmasked wild-type marginals; not pseudo-likelihood)."""
    return run_quality(source_run_dir, candidates_jsonl, tool="esm1v_score", settings={}, cpu_cores=cpu_cores, gpu_device=gpu_device)


def run_esmif_score(source_run_dir: Path, candidates_jsonl: Path, cpu_cores: int = 4, gpu_device: str = "0") -> Path:
    """Native sequence conditional log likelihood using the selected complex backbone as context."""
    return run_quality(source_run_dir, candidates_jsonl, tool="esmif_score", settings={}, cpu_cores=cpu_cores, gpu_device=gpu_device)


def run_plip(source_run_dir: Path, candidates_jsonl: Path | None = None, cpu_cores: int = 4) -> Path:
    """PLIP interactions on the best docking pose per compound or declared cofolded complexes."""
    return run_quality(source_run_dir, candidates_jsonl, tool="plip", settings={}, cpu_cores=cpu_cores)


def run_pandamap(source_run_dir: Path, candidates_jsonl: Path | None = None, cpu_cores: int = 4) -> Path:
    """PandaMap interactions on selected ligand poses, preserving source receptor residue numbering."""
    return run_quality(source_run_dir, candidates_jsonl, tool="pandamap", settings={}, cpu_cores=cpu_cores)


def run_posebusters(source_run_dir: Path, candidates_jsonl: Path | None = None, cpu_cores: int = 4) -> Path:
    """PoseBusters geometry/chemistry checks using a SMILES-matched native ligand pose and receptor."""
    return run_quality(source_run_dir, candidates_jsonl, tool="posebusters", settings={}, cpu_cores=cpu_cores)
