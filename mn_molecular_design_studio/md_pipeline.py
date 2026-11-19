"""Typed, read-only MD sources and analysis outputs for visual pipelines."""
from __future__ import annotations

import hashlib
import math
from pathlib import Path


def validate_analysis_parameters(params):
    params = dict(params)
    for field in ("end_ns", "dcd_frame_interval_ps", "dcd_first_time_ps"):
        try:
            if isinstance(params.get(field), bool):
                raise ValueError("Boolean timing value")
            params[field] = float(params[field]) if str(params.get(field) or "").strip() else None
        except (ValueError, TypeError):
            raise ValueError("MD Analysis " + field + " must be a number or blank.") from None
    for field in ("ligand_selection", "protein_selection"):
        if not str(params.get(field) or "").strip():
            raise ValueError("MD Analysis needs explicit protein and ligand selections.")
    for field, maximum in (("stride", 100000), ("chunk_size", 10000), ("cpu_cores", 256)):
        value = params.get(field)
        if type(value) is not int or not 1 <= value <= maximum:
            raise ValueError("MD Analysis " + field + " must be a positive integer within its supported range.")
    for field in ("start_ns", "contact_cutoff_angstrom", "end_ns", "dcd_frame_interval_ps", "dcd_first_time_ps"):
        value = params.get(field)
        if value is None and field in {"end_ns", "dcd_frame_interval_ps", "dcd_first_time_ps"}:
            continue
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value < 0:
            raise ValueError("MD Analysis " + field + " must be finite and non-negative.")
    if params["contact_cutoff_angstrom"] <= 0 or params.get("end_ns") is not None and params["end_ns"] <= params["start_ns"]:
        raise ValueError("MD Analysis needs a positive contact cutoff and an end time greater than its start time.")
    if (params.get("dcd_frame_interval_ps") is None) != (params.get("dcd_first_time_ps") is None):
        raise ValueError("MD Analysis DCD timing overrides require both interval and first-frame time.")
    if params.get("dcd_frame_interval_ps") is not None and params["dcd_frame_interval_ps"] <= 0:
        raise ValueError("MD Analysis DCD interval must be positive.")
    if params.get("pbc_mode") not in {"auto", "require", "none"}:
        raise ValueError("MD Analysis periodic-boundary handling must be auto, require or none.")
    return params


def _digest(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_source(params, *, owner_id, project_id):
    from mn_molecular_design_studio.auth import development_user, is_development_mode, list_users, require_job_access
    from mn_molecular_design_studio.engine.workflow_engine.workflows.md_analysis import _resolve_source_paths
    from mn_molecular_design_studio.job_catalog import get_job_access_record
    from mn_molecular_design_studio.md_analysis_ui import _source_options_for_job
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

    resolved = resolve_stored_path(params["source_run_dir"], must_exist=True)
    if resolved is None:
        raise ValueError("The MD source job is unavailable on this workstation.")
    root = Path(resolved).resolve(strict=True)
    relative = root.relative_to(runs_root().resolve())
    if len(relative.parts) != 2:
        raise ValueError("Choose one managed local MD Results job.")
    access = get_job_access_record(*relative.parts)
    actor = next((user for user in list_users() if user["id"] == owner_id and not user.get("disabled")), None)
    if actor is None and is_development_mode() and development_user()["id"] == owner_id:
        actor = development_user()
    if actor is None:
        raise PermissionError("The pipeline owner is no longer an active account.")
    require_job_access(actor, access)
    if not access or str(access.get("project_id") or "") != str(project_id):
        raise ValueError("The source MD job must belong to this pipeline project.")
    options = _source_options_for_job({"run_dir": str(root)})
    option = next((row for row in options if row["option_id"] == params["option_id"]), None)
    if option is None:
        raise ValueError("The recorded MD production segment is no longer available.")
    topology = root / params["topology_relative_path"]
    if topology.resolve() not in option["topology_paths"]:
        raise ValueError("The topology does not belong to the selected MD trajectory.")
    source, topology, trajectory = _resolve_source_paths(root, topology, option["trajectory_path"])
    return source, topology, trajectory, relative, option


def source_output(params, *, owner_id, project_id):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

    source, topology, trajectory, relative, option = resolve_source(params, owner_id=owner_id, project_id=project_id)
    segment = option.get("segment") or {}
    lineage = {**segment, "study_id": (option.get("study") or {}).get("study_id"),
        "source_run_relative_path": relative.as_posix()}
    paths = [source / "metadata.json", topology, trajectory]
    manifest = source / "artifacts/md/study.json"
    if manifest.is_file():
        paths.append(manifest)
    digests = {str(path): _digest(path) for path in paths}
    return {"md_source": {"source_run_dir": str(source), "topology_path": str(topology),
        "trajectory_path": str(trajectory), "source_segments": [lineage],
        "expected_source_sha256": {"topology": digests[str(topology)], "trajectory": digests[str(trajectory)]}},
        "source_run_dir": str(source), "target": {},
        "content_manifest": digests,
        "source_metadata": {"run_id": source.name, "task_group": source.parent.name,
            "status": read_json(source / "metadata.json").get("status"), "label": option["label"]}}


def analysis_output(child):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, read_json_raw
    from mn_molecular_design_studio.md_results_ui import _safe_artifact_path

    child = Path(child).resolve()
    result = read_json(child / "result.json")
    summary = read_json_raw(child / "artifacts/md_analysis/summary.json")
    if result.get("success") is not True or summary.get("success") is not True:
        raise ValueError("MD analysis did not produce a successful numerical result.")
    paths = [child / "result.json", child / "artifacts/md_analysis/summary.json"]
    for relative in summary.get("artifacts", {}).values():
        if relative:
            paths.append(_safe_artifact_path(child, relative))
    return {"run_dir": str(child), "target": {}, "analysis_run": {
        "kind": "md_analysis", "task_group": child.parent.name, "run_id": child.name,
        "run_dir": str(child), "success": True, "summary_path": "artifacts/md_analysis/summary.json",
        "source_study_id": summary.get("source_study_id"), "source": summary.get("source"),
        "source_segments": summary.get("source_segments"), "metrics": summary.get("metrics") or {}},
        "content_manifest": {str(path): _digest(path) for path in paths}}
