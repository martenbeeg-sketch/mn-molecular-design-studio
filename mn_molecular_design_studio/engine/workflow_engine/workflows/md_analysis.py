"""Analysis workflow for trajectories imported from legacy MN-Ligand studies."""
from __future__ import annotations

import copy
import json
import math
import shutil
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import run_docker_job
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job,
    finish_job,
    read_json,
    write_json,
)
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root


_RUNTIME_SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "md_analysis.py"


def _source_study(
    source: Path,
    source_segments: list[dict[str, Any]],
    source_run_relative_path: str,
) -> dict[str, Any] | None:
    study_ids = {str(segment.get("study_id") or "") for segment in source_segments}
    study_ids.discard("")
    if len(study_ids) != 1:
        return None
    manifest_path = source / "artifacts" / "md" / "study.json"
    try:
        resolved = manifest_path.resolve(strict=True)
        resolved.relative_to(source)
        manifest = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError, RuntimeError):
        return None
    if not isinstance(manifest, dict) or manifest.get("schema_version") != 1:
        return None
    studies = manifest.get("studies")
    if not isinstance(studies, list):
        return None
    matches = [item for item in studies if isinstance(item, dict) and str(item.get("study_id") or "") in study_ids]
    if len(matches) != 1:
        return None
    study = copy.deepcopy(matches[0])
    history = copy.deepcopy(study.get("source_history")) if isinstance(study.get("source_history"), dict) else {}
    history.update({
        "source_job": {
            "task_group": str(source_run_relative_path).split("/", 1)[0],
            "run_id": str(source_run_relative_path).split("/", 1)[1] if "/" in str(source_run_relative_path) else "",
            "run_relative_path": source_run_relative_path,
        },
        "source_segments": copy.deepcopy(source_segments),
    })
    for key in ("production_runs", "mmgbsa_runs", "preparation_runs", "analysis_runs"):
        value = study.get(key)
        if isinstance(value, list):
            history[key] = copy.deepcopy(value)
        study[key] = []
    study["source_history"] = history
    return study


def _contained(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def _resolve_source_paths(
    source_run_dir: str | Path,
    topology_path: str | Path,
    trajectory_path: str | Path,
) -> tuple[Path, Path, Path]:
    runs = runs_root().expanduser().resolve(strict=False)
    source = resolve_stored_path(str(source_run_dir), must_exist=True)
    if source is None:
        raise ValueError("The selected source MD run is unavailable on this workstation.")
    source = Path(source).expanduser().resolve(strict=True)
    if not source.is_dir() or not _contained(source, runs):
        raise ValueError("The source run must be a managed local Results run owned by this workstation.")
    if len(source.relative_to(runs).parts) != 2:
        raise ValueError("The source run path must identify one Results job directory.")
    metadata = read_json(source / "metadata.json")
    if not metadata or str(metadata.get("status") or "").lower() != "completed":
        raise ValueError("MD analysis requires a completed source MD run.")

    paths: list[Path] = []
    for label, raw in (("Topology", topology_path), ("Trajectory", trajectory_path)):
        candidate = resolve_stored_path(str(raw), run_dir=source, must_exist=True)
        if candidate is None:
            raise ValueError(f"{label} file is unavailable on this workstation.")
        candidate = Path(candidate).expanduser().resolve(strict=True)
        if not candidate.is_file() or not _contained(candidate, source):
            raise ValueError(f"{label} must be a file inside the selected source Results run.")
        paths.append(candidate)
    if paths[0] == paths[1]:
        raise ValueError("Topology and trajectory must be different files.")
    return source, paths[0], paths[1]


def run_md_analysis(
    *,
    source_run_dir: str | Path,
    topology_path: str | Path,
    trajectory_path: str | Path,
    ligand_selection: str,
    protein_selection: str = "protein",
    stride: int = 1,
    start_ns: float = 0.0,
    end_ns: float | None = None,
    contact_cutoff_angstrom: float = 4.5,
    pbc_mode: str = "auto",
    dcd_frame_interval_ps: float | None = None,
    dcd_first_time_ps: float | None = None,
    chunk_size: int = 100,
    cpu_cores: int = 8,
    source_segments: list[dict[str, Any]] | None = None,
) -> Path:
    """Run streamed, protein-aligned RMSD and geometric contact analysis.

    The source topology and trajectory remain in the owning Results job and are
    mounted read-only for the duration of this analysis.  ``source_segments``
    stores caller-provided simulation/checkpoint lineage; duration fields are
    preserved exactly as supplied and are never added or reinterpreted here.
    """
    ligand_selection = str(ligand_selection or "").strip()
    protein_selection = str(protein_selection or "protein").strip()
    if not ligand_selection:
        raise ValueError("Enter an MDTraj atom selection for the ligand.")
    if not protein_selection:
        raise ValueError("Enter an MDTraj atom selection for the protein.")
    if type(stride) is not int or stride < 1:
        raise ValueError("Stride must be a positive integer.")
    if type(chunk_size) is not int or not 1 <= chunk_size <= 10000:
        raise ValueError("Chunk size must be between 1 and 10000 frames.")
    if type(cpu_cores) is not int or not 1 <= cpu_cores <= 256:
        raise ValueError("CPU cores must be between 1 and 256.")
    if pbc_mode not in {"auto", "none", "require"}:
        raise ValueError("PBC mode must be 'auto', 'none', or 'require'.")
    if (dcd_frame_interval_ps is None) != (dcd_first_time_ps is None):
        raise ValueError("DCD timing override requires both interval and first-frame time.")
    if dcd_frame_interval_ps is not None:
        try:
            dcd_frame_interval_ps = float(dcd_frame_interval_ps)
            dcd_first_time_ps = float(dcd_first_time_ps)
        except (TypeError, ValueError) as exc:
            raise ValueError("DCD timing overrides must be numeric picosecond values.") from exc
        if not math.isfinite(dcd_frame_interval_ps) or dcd_frame_interval_ps <= 0:
            raise ValueError("DCD frame interval override must be finite and greater than zero.")
        if not math.isfinite(dcd_first_time_ps) or dcd_first_time_ps < 0:
            raise ValueError("DCD first-frame time override must be finite and non-negative.")
    try:
        start_ns = float(start_ns)
        end_ns = None if end_ns is None else float(end_ns)
        contact_cutoff_angstrom = float(contact_cutoff_angstrom)
    except (TypeError, ValueError) as exc:
        raise ValueError("Time window and contact cutoff must be numeric.") from exc
    if not math.isfinite(start_ns) or start_ns < 0:
        raise ValueError("Start time must be finite and non-negative.")
    if end_ns is not None and (not math.isfinite(end_ns) or end_ns <= start_ns):
        raise ValueError("End time must be finite and greater than start time.")
    if not math.isfinite(contact_cutoff_angstrom) or contact_cutoff_angstrom <= 0:
        raise ValueError("Contact cutoff must be finite and greater than zero.")
    if source_segments is not None and (
        not isinstance(source_segments, list)
        or any(not isinstance(segment, dict) for segment in source_segments)
    ):
        raise ValueError("Source segment lineage must be a list of metadata objects.")
    source, topology, trajectory = _resolve_source_paths(source_run_dir, topology_path, trajectory_path)
    source_ref = source.relative_to(runs_root().expanduser().resolve()).as_posix()
    trajectory_rel = trajectory.relative_to(source).as_posix()
    topology_rel = topology.relative_to(source).as_posix()
    params = {
        "ligand_selection": ligand_selection,
        "protein_selection": protein_selection,
        "stride": stride,
        "start_ns": start_ns,
        "end_ns": end_ns,
        "contact_cutoff_angstrom": contact_cutoff_angstrom,
        "pbc_mode": pbc_mode,
        "dcd_frame_interval_ps": dcd_frame_interval_ps,
        "dcd_first_time_ps": dcd_first_time_ps,
        "chunk_size": chunk_size,
        "cpu_cores": cpu_cores,
        "source_segments": source_segments or [],
    }
    inputs = {
        "source_run_dir": str(source),
        "topology_path": str(topology),
        "trajectory_path": str(trajectory),
    }
    job = create_job("md-analysis", "md_analysis", "md_analysis", inputs, params)
    artifact_root = job.artifacts_dir / "md_analysis"
    artifact_root.mkdir(parents=True, exist_ok=True)
    input_root = job.artifacts_dir / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    staged_runtime = input_root / "md_analysis.py"
    shutil.copy2(_RUNTIME_SCRIPT, staged_runtime)
    config = {
        **params,
        "topology_path": f"/source/{topology_rel}",
        "trajectory_path": f"/source/{trajectory_rel}",
        "topology_source_path": str(topology),
        "trajectory_source_path": str(trajectory),
        "source_run_relative_path": source_ref,
        "output_dir": "/work/artifacts/md_analysis",
    }
    config_path = input_root / "md_analysis_config.json"
    write_json(config_path, config)
    metadata_path = job.run_dir / "metadata.json"
    metadata = read_json(metadata_path)
    metadata.update({
        "workflow": "ligand.md.analysis",
        "engine": "MDTraj (mn-md:cu128)",
        "result_type": "md_trajectory_analysis_result/v1",
        "artifact_contract_version": 1,
        "source_run_dir": str(source),
    })
    write_json(metadata_path, metadata)
    write_json(job.run_dir / "artifact_contract.json", {
        "schema_version": 1,
        "workflow": "ligand.md.analysis",
        "input_type": "md_trajectory_source/v1",
        "output_type": "md_trajectory_analysis_result/v1",
        "inputs": [
            {"name": "source_run", "type": "md_source_run/v1", "path": str(source)},
            {"name": "topology", "type": "md_topology/v1", "path": str(topology)},
            {"name": "trajectory", "type": "md_trajectory/v1", "path": str(trajectory)},
        ],
        "outputs": [
            {"name": "analysis_summary", "type": "md_analysis_summary/v1", "path": "artifacts/md_analysis/summary.json"},
            {"name": "timeseries", "type": "md_analysis_timeseries_csv/v1", "path": "artifacts/md_analysis/timeseries.csv"},
            {"name": "residue_contacts", "type": "md_residue_contact_map_csv/v1", "path": "artifacts/md_analysis/residue_contacts.csv"},
            {"name": "representative_frame", "type": "md_representative_frame_pdb/v1", "path": "artifacts/md_analysis/representative_frame.pdb"},
        ],
        "structure_output": True,
    })
    run_dir_for_mount = source
    return_code = run_docker_job(
        job.run_dir,
        load_manifest("md_analysis"),
        ["python", "/work/artifacts/input/md_analysis.py", "--config", "/work/artifacts/input/md_analysis_config.json"],
        gpu_device="none",
        extra_mounts=[(run_dir_for_mount, "/source", True)],
        environment={"PYTHONUNBUFFERED": "1"},
    )
    summary_path = artifact_root / "summary.json"
    summary = read_json(summary_path) if summary_path.is_file() else {}
    success = bool(return_code == 0 and summary.get("schema_version") == 1 and summary.get("success") is True)
    summary["success"] = success
    source_study = _source_study(source, source_segments or [], source_ref) if success else None
    if source_study is not None:
        summary["source_study_id"] = source_study["study_id"]
        analysis_manifest = {
            "schema_version": 1,
            "studies": [source_study],
        }
        analysis_record = {
            "run_id": job.run_id,
            "status": "completed",
            "result_path": "artifacts/md_analysis/summary.json",
            "source_production_run_id": (source_segments or [{}])[0].get("run_id"),
        }
        source_study["analysis_runs"] = [analysis_record]
        write_json(job.artifacts_dir / "md" / "study.json", analysis_manifest)
    write_json(summary_path, summary)
    error = str(summary.get("error") or "")
    if not success and not error:
        stderr_path = job.run_dir / "stderr.log"
        error = stderr_path.read_text(encoding="utf-8", errors="replace")[-4000:] if stderr_path.is_file() else "MD analysis did not produce a valid summary."
    expected_outputs = {
        "analysis_summary": "artifacts/md_analysis/summary.json",
        "timeseries": "artifacts/md_analysis/timeseries.csv",
        "residue_contacts": "artifacts/md_analysis/residue_contacts.csv",
    }
    if (artifact_root / "representative_frame.pdb").is_file():
        expected_outputs["representative_frame"] = "artifacts/md_analysis/representative_frame.pdb"
    finish_job(job.run_dir, success, {
        "outputs": expected_outputs if success else {},
        "metrics": {
            "engine": "MDTraj",
            "source_run_relative_path": source_ref,
            "source_trajectory_relative_path": trajectory_rel,
            "trajectory_atom_count": summary.get("trajectory_atom_count"),
            "selected_frame_count": summary.get("selected_frame_count"),
            "selected_time_range_ns": summary.get("selected_time_range_ns"),
            "return_code": return_code,
            "error": error,
        },
        "downstream_artifacts": {"md_trajectory_analysis": "artifacts/md_analysis/summary.json"} if success else {},
        "summary": summary,
    })
    return job.run_dir
