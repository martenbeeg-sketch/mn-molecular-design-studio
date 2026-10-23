"""Reusable setup tile for streamed analysis of an existing MD trajectory."""
from __future__ import annotations

import json
import os
from pathlib import Path, PurePosixPath
from typing import Any, Callable

import streamlit as st


_TRAJECTORY_SUFFIXES = {".dcd", ".xtc", ".trr", ".nc", ".netcdf", ".h5", ".hdf5"}
_TOPOLOGY_SUFFIXES = {".pdb", ".ent", ".gro", ".cif", ".mmcif", ".prmtop", ".parm7", ".psf", ".top"}
_TRAJECTORY_KEYS = ("production_trajectory", "native_production_trajectory", "trajectory", "trajectory_path")
_TOPOLOGY_KEYS = (
    "production_pdb", "production_topology", "production_gro", "system_pdb", "npt_pdb",
    "protein_prepared", "prepared_complex", "topology", "topology_path", "system_topology",
)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _finite_number(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        import math
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def _safe_run_file(root: Path, relative: str) -> Path | None:
    """Resolve a relative artifact path without following it outside the job."""
    value = str(relative or "")
    if not value or "\\" in value or "\x00" in value:
        return None
    logical = PurePosixPath(value)
    if logical.is_absolute() or any(part in {"", ".", ".."} for part in logical.parts):
        return None
    candidate = root.joinpath(*logical.parts)
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root.resolve(strict=True))
    except (OSError, RuntimeError, ValueError):
        return None
    return resolved if resolved.is_file() else None


def _path_candidates_from_output(
    source_root: Path,
    value: Any,
    *,
    segment_run_id: str,
    known_legacy_run_ids: set[str],
) -> list[Path]:
    """Map an old container output path to an exact in-run artifact path.

    Legacy output references include paths such as ``/output/<run-id>/x.dcd``
    and ``/prepared-system/<prep-run-id>/system.pdb``. The ID and full suffix
    are honored; a same-named file elsewhere is never substituted.
    """
    raw = str(value or "").strip()
    if not raw or "\\" in raw or "\x00" in raw:
        return []
    path = PurePosixPath(raw)
    choices: list[str] = []
    if not path.is_absolute():
        relative = path.as_posix()
        if relative.startswith("artifacts/"):
            choices.append(relative)
        if segment_run_id:
            choices.append(f"artifacts/legacy/{segment_run_id}/{relative}")
            choices.append(f"artifacts/md/{relative}")
            choices.append(f"artifacts/{relative}")
    else:
        parts = path.parts[1:]
        if parts[:1] == ("work",):
            choices.append(PurePosixPath(*parts[1:]).as_posix())
        if len(parts) >= 3 and parts[0] in {"output", "prepared-system", "source", "prepared-system-output"}:
            run_id = parts[1]
            if run_id in known_legacy_run_ids and parts[2:]:
                choices.append(f"artifacts/legacy/{run_id}/" + PurePosixPath(*parts[2:]).as_posix())
            # Native output directories can live at the job's artifact root.
            if run_id == segment_run_id and parts[2:]:
                choices.append("artifacts/" + PurePosixPath(*parts[2:]).as_posix())
        # Accept an original host path only if it points inside this owner job.
        try:
            resolved = Path(raw).expanduser().resolve(strict=True)
            resolved.relative_to(source_root.resolve(strict=True))
            if resolved.is_file():
                return [resolved]
        except (OSError, RuntimeError, ValueError):
            pass
    found: list[Path] = []
    seen: set[Path] = set()
    for relative in choices:
        candidate = _safe_run_file(source_root, relative)
        if candidate is not None and candidate not in seen:
            found.append(candidate)
            seen.add(candidate)
    return found


def _walk_files(root: Path, suffixes: set[str], *, limit: int = 20_000) -> list[Path]:
    """Enumerate candidate paths only; never open coordinate files in the UI."""
    try:
        safe_root = root.resolve(strict=True)
        if not safe_root.is_dir():
            return []
    except (OSError, RuntimeError):
        return []
    found: list[Path] = []
    for current, dirnames, filenames in os.walk(safe_root, followlinks=False):
        current_path = Path(current)
        dirnames[:] = [name for name in dirnames if not (current_path / name).is_symlink()]
        for filename in filenames:
            path = current_path / filename
            if path.is_symlink() or path.suffix.lower() not in suffixes:
                continue
            try:
                resolved = path.resolve(strict=True)
                resolved.relative_to(safe_root)
            except (OSError, RuntimeError, ValueError):
                continue
            if resolved.is_file():
                found.append(resolved)
                if len(found) >= limit:
                    return sorted(set(found))
    return sorted(set(found))


def _result_outputs(result: dict[str, Any]) -> dict[str, Any]:
    md_result = result.get("md_result") if isinstance(result.get("md_result"), dict) else {}
    outputs = md_result.get("output_files") if isinstance(md_result.get("output_files"), dict) else {}
    if not outputs and isinstance(result.get("output_files"), dict):
        outputs = result["output_files"]
    return outputs


def _unique_paths(paths: list[Path]) -> list[Path]:
    result: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            result.append(resolved)
            seen.add(resolved)
    return result


def _artifact_choices(
    source_root: Path,
    *,
    result: dict[str, Any],
    segment_run_id: str,
    legacy_run_ids: set[str],
    native: bool,
) -> tuple[list[Path], list[Path]]:
    outputs = _result_outputs(result)
    trajectories: list[Path] = []
    topologies: list[Path] = []
    trajectory_keys = _TRAJECTORY_KEYS if native else ("production_trajectory",)
    topology_keys = _TOPOLOGY_KEYS if native else (
        "production_pdb", "production_topology", "production_gro", "system_pdb",
        "npt_pdb", "protein_prepared", "prepared_complex",
    )
    for key in trajectory_keys:
        raw = outputs.get(key)
        if raw:
            trajectories.extend(_path_candidates_from_output(
                source_root, raw, segment_run_id=segment_run_id,
                known_legacy_run_ids=legacy_run_ids,
            ))
    for key in topology_keys:
        raw = outputs.get(key)
        if raw:
            topologies.extend(_path_candidates_from_output(
                source_root, raw, segment_run_id=segment_run_id,
                known_legacy_run_ids=legacy_run_ids,
            ))
    if native:
        search_roots = [source_root / "artifacts" / "md", source_root / "artifacts" / "input"]
    else:
        # Imported segments must be resolved only through the exact paths in
        # this production result. Never substitute a file found in this or
        # another segment directory, even when its basename happens to match.
        search_roots = []
    # Native jobs may lack a normalized output path, so a scan within their
    # own artifacts remains a bounded fallback. Imported jobs may not.
    if native and not trajectories:
        for root in search_roots:
            trajectories.extend(_walk_files(root, _TRAJECTORY_SUFFIXES))
    if native and not topologies:
        for root in search_roots:
            topologies.extend(_walk_files(root, _TOPOLOGY_SUFFIXES))
    trajectories = _unique_paths(trajectories)
    topologies = _unique_paths(topologies)
    return trajectories, topologies


def _output_duration(result: dict[str, Any]) -> float | None:
    md_result = result.get("md_result") if isinstance(result.get("md_result"), dict) else {}
    production = ((md_result.get("equilibration_stats") or {}).get("production")
                  if isinstance(md_result.get("equilibration_stats"), dict) else None)
    production = production if isinstance(production, dict) else {}
    for value in (production.get("duration_ns"), md_result.get("production_duration_ns"), result.get("production_ns")):
        number = _finite_number(value)
        if number is not None and number > 0:
            return number
    return None


def _segment_record(run_id: str, production: dict[str, Any], result: dict[str, Any]) -> dict[str, Any]:
    duration = _finite_number(production.get("segment_duration_ns"))
    if duration is None:
        duration = _output_duration(result)
    if duration is not None and duration <= 0:
        duration = None
    return {
        "run_id": str(run_id),
        "continuation_of_run_id": str(production.get("continuation_of_run_id") or "") or None,
        "previous_duration_ns": _finite_number(production.get("previous_duration_ns")),
        "target_duration_ns": _finite_number(production.get("target_duration_ns")),
        "segment_duration_ns": duration,
        "trajectory_appended": production.get("trajectory_appended"),
        "duration_status": "recorded" if duration is not None else "unknown",
    }


def _identity_label(study: dict[str, Any], job: dict[str, Any]) -> str:
    compound = str(study.get("compound_name") or study.get("compound_id") or job.get("campaign_name") or "compound unknown")
    target = str(study.get("target_pdb_id") or job.get("pdb_id") or "target unknown")
    pose = str(study.get("source_pose_id") or "pose unknown")
    tool = str(study.get("source_tool") or job.get("tool") or "source tool unknown")
    rank = study.get("source_pose_rank")
    rank_label = f" · rank {rank}" if rank is not None else ""
    return f"{compound} · {target} · {pose}{rank_label} · {tool}"


def _study_source_options(source_root: Path, job: dict[str, Any]) -> list[dict[str, Any]]:
    """Build one independently selectable option per recorded production segment."""
    try:
        source_root = source_root.resolve(strict=True)
        manifest_path = _safe_run_file(source_root, "artifacts/md/study.json")
    except (OSError, RuntimeError):
        return []
    if manifest_path is None:
        return []
    manifest = _read_json(manifest_path)
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("studies"), list):
        return []
    output: list[dict[str, Any]] = []
    for study_index, study in enumerate(manifest["studies"]):
        if not isinstance(study, dict):
            continue
        production_runs = study.get("production_runs")
        if not isinstance(production_runs, list):
            continue
        known_ids = {
            str(item.get("run_id")) for item in production_runs if isinstance(item, dict) and item.get("run_id")
        }
        preparations = study.get("preparation_runs") if isinstance(study.get("preparation_runs"), list) else []
        known_ids.update(str(item.get("run_id")) for item in preparations if isinstance(item, dict) and item.get("run_id"))
        for production in production_runs:
            if not isinstance(production, dict) or str(production.get("status") or "").lower() not in {"completed", "success", "succeeded"}:
                continue
            segment_id = str(production.get("run_id") or "")
            result_path = _safe_run_file(source_root, str(production.get("result_path") or ""))
            # New normalized manifests store exact run-relative trajectory and
            # topology paths. Use them directly and avoid opening potentially
            # large legacy result files on every UI rerun.
            normalized_trajectory = production.get("trajectory_path")
            normalized_topology = production.get("topology_path")
            production_result = {}
            if normalized_trajectory:
                trajectory = _safe_run_file(source_root, str(normalized_trajectory))
                trajectories = [trajectory] if trajectory else []
            else:
                production_result = _read_json(result_path) if result_path else {}
                trajectories, _unused = _artifact_choices(
                    source_root, result=production_result, segment_run_id=segment_id,
                    legacy_run_ids=known_ids, native=False,
                )
            if normalized_topology:
                topology = _safe_run_file(source_root, str(normalized_topology))
                topologies = [topology] if topology else []
            else:
                if not production_result:
                    production_result = _read_json(result_path) if result_path else {}
                _unused, topologies = _artifact_choices(
                    source_root, result=production_result, segment_run_id=segment_id,
                    legacy_run_ids=known_ids, native=False,
                )
            if not segment_id or not trajectories or not topologies:
                continue
            segment = _segment_record(segment_id, production, production_result)
            for trajectory in trajectories:
                output.append({
                    "option_id": f"study:{study_index}:{segment_id}:{trajectory.relative_to(source_root).as_posix()}",
                    "label": f"{_identity_label(study, job)} · segment {segment_id} · {trajectory.name}",
                    "study": study,
                    "segment": segment,
                    "trajectory_path": trajectory,
                    "topology_paths": topologies,
                    "source_result_path": result_path.relative_to(source_root).as_posix() if result_path else "",
                    "lineage_known": True,
                })
    return output


def _native_source_options(source_root: Path, job: dict[str, Any]) -> list[dict[str, Any]]:
    result = _read_json(source_root / "result.json")
    metadata = _read_json(source_root / "metadata.json")
    md_result = result.get("md_result") if isinstance(result.get("md_result"), dict) else {}
    production_stats = ((md_result.get("equilibration_stats") or {}).get("production")
                        if isinstance(md_result.get("equilibration_stats"), dict) else None)
    production_stats = production_stats if isinstance(production_stats, dict) else {}
    production = {
        "run_id": str(metadata.get("run_id") or job.get("run_id") or source_root.name),
        "continuation_of_run_id": metadata.get("continuation_of_run_id"),
        "previous_duration_ns": metadata.get("previous_duration_ns"),
        "target_duration_ns": metadata.get("target_duration_ns"),
        "segment_duration_ns": production_stats.get("duration_ns"),
        "trajectory_appended": production_stats.get("trajectory_appended"),
    }
    segment_id = production["run_id"]
    trajectories, topologies = _artifact_choices(
        source_root, result=result, segment_run_id=segment_id,
        legacy_run_ids=set(), native=True,
    )
    segment = _segment_record(segment_id, production, result)
    if not segment.get("segment_duration_ns") and not segment.get("target_duration_ns"):
        segment["duration_status"] = "unknown"
    identity = {
        "compound_name": metadata.get("campaign_name") or job.get("campaign_name"),
        "target_pdb_id": metadata.get("pdb_id") or result.get("pdb_id"),
        "source_pose_id": metadata.get("source_pose_id") or metadata.get("selected_ligand_key"),
        "source_tool": metadata.get("tool") or job.get("tool"),
    }
    output = []
    for trajectory in trajectories:
        output.append({
            "option_id": f"native:{trajectory.relative_to(source_root).as_posix()}",
            "label": f"{_identity_label(identity, job)} · segment {segment_id} · {trajectory.name}",
            "study": identity,
            "segment": segment,
            "trajectory_path": trajectory,
            "topology_paths": topologies,
            "source_result_path": "result.json",
            "lineage_known": False,
        })
    return output


def _source_options_for_job(job: dict[str, Any]) -> list[dict[str, Any]]:
    source_root = Path(str(job.get("run_dir") or "")).expanduser()
    try:
        source_root = source_root.resolve(strict=True)
    except (OSError, RuntimeError):
        return []
    study_manifest = source_root / "artifacts" / "md" / "study.json"
    options = _study_source_options(source_root, job)
    if options or study_manifest.exists():
        # A present aggregate manifest is authoritative. Do not reinterpret a
        # broken or incomplete aggregate as one untracked native trajectory.
        return options
    # Native Studio MD jobs have no cross-run study manifest. Their one
    # trajectory is handled as one segment; an unrecorded duration stays unknown.
    return _native_source_options(source_root, job)


@st.cache_data(ttl=45, max_entries=64, show_spinner=False)
def _cached_source_option_rows(
    source_root_text: str,
    manifest_signature: tuple[int, int] | None,
    result_signature: tuple[int, int] | None,
    metadata_signature: tuple[int, int] | None,
    campaign_name: str,
    tool: str,
) -> list[dict[str, Any]]:
    """Cache path-only discovery; coordinate files are never opened here."""
    source_root = Path(source_root_text)
    row = {"run_dir": source_root_text, "campaign_name": campaign_name, "tool": tool}
    options = _source_options_for_job(row)
    serialized: list[dict[str, Any]] = []
    for option in options:
        entry = {key: value for key, value in option.items() if key not in {"trajectory_path", "topology_paths"}}
        entry["trajectory_relative_path"] = Path(option["trajectory_path"]).relative_to(source_root).as_posix()
        entry["topology_relative_paths"] = [
            Path(path).relative_to(source_root).as_posix() for path in option.get("topology_paths", [])
        ]
        serialized.append(entry)
    return serialized


def _cached_source_options_for_job(job: dict[str, Any]) -> list[dict[str, Any]]:
    source_root = Path(str(job.get("run_dir") or "")).expanduser().resolve(strict=True)

    def signature(path: Path) -> tuple[int, int] | None:
        try:
            info = path.stat()
            return info.st_mtime_ns, info.st_size
        except OSError:
            return None

    study_path = source_root / "artifacts" / "md" / "study.json"
    result_path = source_root / "result.json"
    metadata_path = source_root / "metadata.json"
    records = _cached_source_option_rows(
        str(source_root), signature(study_path), signature(result_path), signature(metadata_path),
        str(job.get("campaign_name") or ""), str(job.get("tool") or ""),
    )
    output = []
    for record in records:
        option = dict(record)
        option["trajectory_path"] = _safe_run_file(source_root, record["trajectory_relative_path"])
        option["topology_paths"] = [
            path for relative in record.get("topology_relative_paths", [])
            if (path := _safe_run_file(source_root, relative)) is not None
        ]
        if option["trajectory_path"] is not None and option["topology_paths"]:
            output.append(option)
    return output


def md_source_options_for_job(job: dict[str, Any], user: dict[str, Any]) -> dict[str, Any]:
    """Return safe, serializable options for an authorized source MD job.

    Paths are run-relative rather than filesystem paths so an authenticated
    peer API can present these choices without transferring coordinate data.
    Callers still submit the owning job and exact relative paths back to the
    owner workstation for execution.
    """
    from mn_molecular_design_studio.auth import require_job_access
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root
    from mn_molecular_design_studio.job_catalog import get_job_access_record

    group = str(job.get("task_group") or "")
    run_id = str(job.get("run_id") or "")
    if group != "ligand-md" or not run_id:
        raise ValueError("The source must be a cataloged ligand-md job.")
    record = get_job_access_record(group, run_id)
    if not record:
        raise PermissionError("MD source job is not available in the job catalog.")
    require_job_access(user, record)
    if str(record.get("project_id") or "") != str(job.get("project_id") or ""):
        raise PermissionError("MD source job project does not match its access record.")
    source_root = Path(str(job.get("run_dir") or "")).expanduser().resolve(strict=True)
    managed_root = runs_root().expanduser().resolve(strict=True)
    try:
        relative_source = source_root.relative_to(managed_root)
    except ValueError as exc:
        raise ValueError("MD source job must be inside managed Results.") from exc
    if len(relative_source.parts) != 2:
        raise ValueError("MD source must identify one managed Results job directory.")
    if not source_root.is_dir():
        raise ValueError("MD source job Results directory is unavailable.")
    metadata = _read_json(source_root / "metadata.json")
    if str(metadata.get("status") or job.get("status") or "").lower() not in {"completed", "succeeded"}:
        raise ValueError("MD source options are available only for completed jobs.")

    safe_options = []
    for option in _cached_source_options_for_job({**job, "run_dir": str(source_root)}):
        topology_paths = option.get("topology_paths") or []
        trajectory = Path(option["trajectory_path"]).resolve(strict=True)
        try:
            trajectory_rel = trajectory.relative_to(source_root).as_posix()
        except ValueError:
            continue
        study = option.get("study") if isinstance(option.get("study"), dict) else {}
        identity_keys = (
            "study_id", "compound_id", "compound_name", "target_pdb_id", "target_receptor",
            "source_pose_id", "source_pose_rank", "source_tool", "ligand_resname",
            "identity_status", "selection_status",
        )
        for topology_path in topology_paths:
            topology = Path(topology_path).resolve(strict=True)
            try:
                topology_rel = topology.relative_to(source_root).as_posix()
            except ValueError:
                continue
            safe_options.append({
                "option_id": f"{option['option_id']}|topology:{topology_rel}",
                "label": f"{option['label']} · topology {topology_rel}",
                "study": {key: study[key] for key in identity_keys if key in study and study[key] is not None},
                "study_id": str(study.get("study_id") or "") or None,
                "segment": dict(option.get("segment") or {}),
                "trajectory_relative_path": trajectory_rel,
                "topology_relative_path": topology_rel,
                "source_result_path": str(option.get("source_result_path") or ""),
                "lineage_known": bool(option.get("lineage_known")),
            })
    return {
        "task_group": group,
        "run_id": run_id,
        "project_id": str(record.get("project_id") or ""),
        "owner_id": str(record.get("owner_id") or ""),
        "source_label": str(job.get("campaign_name") or job.get("display_tool") or job.get("tool") or run_id),
        "source_options": safe_options,
    }


def _list_accessible_md_jobs(user: dict[str, Any], project_id: str) -> list[dict[str, Any]]:
    from mn_molecular_design_studio.auth import accessible_project_ids, require_job_access
    from mn_molecular_design_studio.job_catalog import get_job_access_record, search_jobs
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root

    managed_runs_root = runs_root().expanduser().resolve(strict=True)

    rows = search_jobs(
        task_group="ligand-md",
        statuses=["completed", "succeeded"],
        project_id=project_id,
        access_owner_id=str(user.get("id") or ""),
        accessible_project_ids=accessible_project_ids(user),
        access_admin=user.get("role") == "admin",
        limit=1000,
    ).get("items", [])
    available: list[dict[str, Any]] = []
    for row in rows:
        group, run_id = str(row.get("task_group") or ""), str(row.get("run_id") or "")
        record = get_job_access_record(group, run_id)
        try:
            require_job_access(user, record)
        except PermissionError:
            continue
        if not record or str(record.get("project_id") or "") != str(project_id or ""):
            continue
        source_root = Path(str(row.get("run_dir") or ""))
        try:
            source_root = source_root.resolve(strict=True)
            source_root.relative_to(managed_runs_root)
        except (OSError, RuntimeError, ValueError):
            # Catalog summaries for another workstation intentionally carry no
            # local run_dir. Also reject a catalog path outside managed Results.
            continue
        if not source_root.is_dir():
            # Catalog summaries for another workstation intentionally carry no
            # local run_dir. The UI never pulls a large trajectory across peers.
            continue
        available.append(row)
    return available


def _list_remote_md_jobs(user: dict[str, Any], project_id: str, token: str) -> list[dict[str, Any]]:
    from mn_molecular_design_studio.cluster import config as cluster_config
    from mn_molecular_design_studio.cluster.service import shared_jobs

    if not cluster_config.enabled():
        return []
    payload = shared_jobs(
        user, token=token, project_id=project_id,
        statuses=["completed", "succeeded"], limit=1000,
    )
    local_node = cluster_config.identity()
    output = []
    for row in payload.get("items", []):
        if str(row.get("task_group") or "") != "ligand-md":
            continue
        node_id = str(row.get("node_id") or "")
        if not node_id or node_id == local_node or str(row.get("project_id") or "") != project_id:
            continue
        local_reference = str(row.get("local_reference") or f"{row['task_group']}/{row['run_id']}")
        row["_md_remote"] = True
        row["_md_cluster_job_id"] = str(row.get("cluster_job_id") or cluster_config.job_id(node_id, local_reference))
        row["_md_source_node"] = node_id
        row["_md_source_reference"] = local_reference
        output.append(row)
    return output


def _remote_source_options_fragment(row: dict[str, Any], token: str, operation_key: str, options_key: str) -> None:
    from mn_molecular_design_studio.cluster import ui as cluster_ui

    command_id = st.session_state.get(operation_key)
    if not command_id:
        return

    @st.fragment(run_every=2)
    def poll_source_metadata() -> None:
        try:
            operation = cluster_ui._request("/commands/" + str(command_id), token)
        except Exception as exc:
            st.warning(f"Could not reach the source workstation: {exc}")
            return
        if operation.get("status") == "failed":
            st.error(operation.get("error") or "The source workstation could not list this MD source.")
            return
        if operation.get("status") != "completed":
            st.info("Waiting for the owning workstation to inspect the recorded MD source paths…")
            return
        response = operation.get("response") or {}
        options = response.get("source_options") if isinstance(response, dict) else None
        if not isinstance(options, list):
            st.error("The owner returned an unsupported MD source-options response.")
            return
        st.session_state[options_key] = options
        st.session_state.pop(operation_key, None)
        st.rerun(scope="app")

    poll_source_metadata()


def _portable_source_path(path: Path, source_run_dir: Path, row: dict[str, Any]) -> str:
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root

    resolved = path.resolve(strict=True)
    root = runs_root().expanduser().resolve(strict=True)
    try:
        relative = resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError("The selected MD file is outside this workstation's managed Results.") from exc
    source_root = source_run_dir.expanduser().resolve(strict=True)
    if len(relative.parts) < 3 and resolved != source_root:
        raise ValueError("Select a file inside a managed Results job.")
    try:
        from mn_molecular_design_studio.cluster import config as cluster_config
        if cluster_config.enabled():
            return cluster_config.artifact_uri(cluster_config.identity(), relative.parts[0], relative.parts[1], "/".join(relative.parts[2:]))
    except (ImportError, OSError, ValueError):
        pass
    return str(resolved)


def render_md_analysis_setup(
    schema: dict[str, Any],
    *,
    user: dict[str, Any],
    owner_id: str | None = None,
    project_selector: Callable[..., str | None],
    submit: Callable[..., None],
    title: str = "MD Analysis",
) -> None:
    """Render the setup form and submit through the Studio workflow callback."""
    from mn_molecular_design_studio.auth import accessible_project_ids, is_development_mode
    owner_id = str(owner_id or user.get("id") or "")
    workflow_id = "ligand.md.analysis"
    project_id = project_selector("md-analysis")
    st.write("Analyze one selected MD production segment at a time. The source trajectory stays on its owning workstation and is mounted read-only; it is never downloaded into the analysis job.")
    st.caption("Each input file is treated as its own recorded segment. Previous and target durations are shown as lineage metadata and are not added together or presented as one continuous trajectory.")
    if not project_id:
        st.info("Select a project before choosing a source MD job.")
        return
    if project_id not in set(accessible_project_ids(user)) and user.get("role") != "admin":
        st.error("You do not have access to the selected project.")
        return
    try:
        jobs = _list_accessible_md_jobs(user, str(project_id))
    except Exception as exc:
        st.error(f"Could not load accessible MD Results jobs: {exc}")
        return
    token = str(st.session_state.get("studio-session-token") or ("development-mode" if is_development_mode() else ""))
    try:
        remote_jobs = _list_remote_md_jobs(user, str(project_id), token)
    except Exception as exc:
        st.warning(f"Could not list MD jobs on other workstations: {exc}")
        remote_jobs = []
    if not jobs and not remote_jobs:
        st.info("No completed MD source jobs are available in this project.")
        return
    modes = (["This workstation"] if jobs else []) + (["Other workstation"] if remote_jobs else [])
    mode = st.radio("MD source location", modes, horizontal=True, key="md-analysis-source-location") if len(modes) > 1 else modes[0]
    is_remote = mode == "Other workstation"
    if not is_remote:
        job_labels = [
            f"{row.get('campaign_name') or row.get('display_tool') or row.get('tool') or 'MD run'} · {row.get('job_code') or row.get('run_id')}"
            for row in jobs
        ]
        selected_job_label = st.selectbox("Source MD job", job_labels, key="md-analysis-source-job")
        row = jobs[job_labels.index(selected_job_label)]
        try:
            job_options = _cached_source_options_for_job(row)
        except (OSError, ValueError, RuntimeError) as exc:
            st.error(f"Could not inspect source Results metadata: {exc}")
            return
        if not job_options:
            st.info("This completed MD job has no exact, usable trajectory and topology pair in its recorded outputs.")
            return
        option_labels = [option["label"] for option in job_options]
        selected_option_label = st.selectbox("Production segment and trajectory", option_labels, key=f"md-analysis-source-segment-{row.get('run_id')}")
        selected = job_options[option_labels.index(selected_option_label)]
        if len(selected["topology_paths"]) > 1:
            topology_labels = [path.relative_to(Path(row["run_dir"]).resolve()).as_posix() for path in selected["topology_paths"]]
            default_index = next((i for i, label in enumerate(topology_labels) if Path(label).name.startswith("production")), 0)
            topology_label = st.selectbox("Topology for this trajectory", topology_labels, index=default_index,
                                          key=f"md-analysis-topology-{selected['option_id']}")
            topology_path = selected["topology_paths"][topology_labels.index(topology_label)]
        else:
            topology_path = selected["topology_paths"][0]
        trajectory_path = selected["trajectory_path"]
        relative_trajectory = trajectory_path.relative_to(Path(row["run_dir"]).resolve()).as_posix()
        relative_topology = topology_path.relative_to(Path(row["run_dir"]).resolve()).as_posix()
        source_run_dir = Path(str(row["run_dir"])).resolve(strict=True)
    else:
        from mn_molecular_design_studio.cluster import config as cluster_config
        from mn_molecular_design_studio.cluster import ui as cluster_ui

        remote_labels = [
            f"{row.get('campaign_name') or row.get('tool') or 'MD run'} · {row.get('job_code') or row.get('run_id')} · {row.get('node_name') or row.get('_md_source_node')}"
            for row in remote_jobs
        ]
        selected_remote_label = st.selectbox("Source MD job", remote_labels, key="md-analysis-remote-source-job")
        row = remote_jobs[remote_labels.index(selected_remote_label)]
        operation_key = f"md-analysis-source-options-op-{row['_md_cluster_job_id']}"
        options_key = f"md-analysis-source-options-{row['_md_cluster_job_id']}"
        remote_options = st.session_state.get(options_key)
        if not remote_options and not st.session_state.get(operation_key):
            if st.button("Load recorded segments from owning workstation", key=f"md-analysis-request-source-{row['_md_cluster_job_id']}"):
                try:
                    operation = cluster_ui._request(
                        "/jobs/" + row["_md_cluster_job_id"] + "/actions",
                        token,
                        {"action": "md_source_options"},
                        method="POST",
                    )
                    st.session_state[operation_key] = operation["command_id"]
                except Exception as exc:
                    st.error(f"Could not request source metadata: {exc}")
        if st.session_state.get(operation_key):
            _remote_source_options_fragment(row, token, operation_key, options_key)
        remote_options = st.session_state.get(options_key)
        if not remote_options:
            st.caption("Only small source metadata is requested. Trajectory and topology files stay on the owning workstation.")
            return
        option_labels = [str(option.get("label") or option.get("option_id")) for option in remote_options]
        selected_option_label = st.selectbox("Production segment and trajectory", option_labels, key=f"md-analysis-remote-segment-{row['_md_cluster_job_id']}")
        selected = dict(remote_options[option_labels.index(selected_option_label)])
        study = selected.get("study") if isinstance(selected.get("study"), dict) else {}
        selected["study"] = study
        selected["topology_paths"] = [selected["topology_relative_path"]]
        selected["trajectory_path"] = selected["trajectory_relative_path"]
        segment = selected.get("segment") if isinstance(selected.get("segment"), dict) else {}
        selected["segment"] = segment
        relative_trajectory = str(selected["trajectory_relative_path"])
        relative_topology = str(selected["topology_relative_path"])
        owner_node = str(row["_md_source_node"])
        source_run_dir = cluster_config.artifact_uri(owner_node, str(row["task_group"]), str(row["run_id"]))
        trajectory_path = cluster_config.artifact_uri(owner_node, str(row["task_group"]), str(row["run_id"]), relative_trajectory)
        topology_path = cluster_config.artifact_uri(owner_node, str(row["task_group"]), str(row["run_id"]), relative_topology)
        st.caption("Source metadata was resolved by the owner. No trajectory bytes were transferred.")
    study = selected["study"]
    segment = selected["segment"]
    st.markdown("**Source identity and ownership**")
    st.write({
        "Compound": study.get("compound_name") or study.get("compound_id") or row.get("campaign_name") or "not recorded",
        "Target": study.get("target_pdb_id") or row.get("pdb_id") or "not recorded",
        "Pose": study.get("source_pose_id") or study.get("selected_ligand_key") or "not recorded",
        "Data owner": str(row.get("owner_id") or owner_id or "unknown"),
        "Project": str(row.get("project_id") or project_id),
        "Source workstation": row.get("node_name") or ("This workstation" if not is_remote else row.get("_md_source_node")),
    })
    duration = segment.get("segment_duration_ns")
    if duration is None:
        st.warning("The segment duration is unknown in the source records. The analysis window will use trajectory frame times and will not infer a duration from target duration or checkpoint lineage.")
    else:
        st.caption(f"Selected file segment: {duration:g} ns (recorded).")
    lineage = {
        "Original production run": segment.get("run_id"),
        "Continues from": segment.get("continuation_of_run_id") or "not recorded",
        "Previous total duration (ns)": segment.get("previous_duration_ns") if segment.get("previous_duration_ns") is not None else "not recorded",
        "Target total duration (ns)": segment.get("target_duration_ns") if segment.get("target_duration_ns") is not None else "not recorded",
        "Current file segment duration (ns)": duration if duration is not None else "unknown",
        "Trajectory appended in source": segment.get("trajectory_appended") if segment.get("trajectory_appended") is not None else "not recorded",
    }
    st.dataframe([{"Lineage field": key, "Recorded value": value} for key, value in lineage.items()], hide_index=True, use_container_width=True)
    if not selected.get("lineage_known"):
        st.info("This native MD job has no study manifest. The selected output is kept as one segment; unrecorded duration and continuation history remain unknown.")
    st.caption(f"Trajectory: {relative_trajectory}")

    selected_ligand = str(study.get("ligand_resname") or study.get("selected_ligand_resname") or "LIG")
    ligand_selection = st.text_input(
        "Ligand atom selection (MDTraj syntax)", value=f"resname {selected_ligand}",
        key=f"md-analysis-ligand-selection-{selected['option_id']}",
        help="Use the residue name or an MDTraj selection such as `resname LIG and not element H`.",
    ).strip()
    protein_selection = st.text_input(
        "Protein atom selection (MDTraj syntax)", value="protein",
        key=f"md-analysis-protein-selection-{selected['option_id']}",
    ).strip()
    control_cols = st.columns(4)
    with control_cols[0]:
        stride = int(st.number_input("Frame stride", min_value=1, max_value=100_000, value=1, step=1, key=f"md-analysis-stride-{selected['option_id']}"))
    with control_cols[1]:
        start_ns = float(st.number_input("Window start (ns)", min_value=0.0, value=0.0, step=0.1, key=f"md-analysis-start-{selected['option_id']}"))
    with control_cols[2]:
        cutoff = float(st.number_input("Contact cutoff (Å)", min_value=0.1, max_value=20.0, value=4.5, step=0.1, key=f"md-analysis-cutoff-{selected['option_id']}"))
    with control_cols[3]:
        chunk_size = int(st.number_input("Frames per chunk", min_value=1, max_value=10000, value=100, step=10, key=f"md-analysis-chunk-{selected['option_id']}"))
    full_window = st.checkbox("Analyze through the last recorded frame", value=True, key=f"md-analysis-full-window-{selected['option_id']}")
    if full_window:
        end_ns = None
    else:
        # Segment duration is elapsed length, while the analysis time axis can
        # be absolute (for example, OpenMM DCD firstStep includes equilibration).
        # Never use elapsed duration as an absolute end-time default.
        default_end = start_ns + 1.0
        end_ns = float(st.number_input("Window end (ns)", min_value=start_ns + 0.001, value=max(start_ns + 0.001, default_end), step=0.1, key=f"md-analysis-end-{selected['option_id']}"))
    st.caption("The analysis window uses the trajectory’s absolute time axis. For OpenMM DCD files this comes from the CHARMM header first step and interval, so it can include equilibration time before the production segment. The workflow does not reset timestamps to zero or join separate segment files.")
    pbc_mode = st.selectbox(
        "Periodic-boundary handling",
        ["auto", "require", "none"],
        format_func=lambda value: {
            "auto": "Use unit cell and bonded topology when available; otherwise keep stored coordinates",
            "require": "Require unit-cell vectors and bonded connectivity",
            "none": "Use stored coordinates without periodic correction",
        }[value],
        key=f"md-analysis-pbc-{selected['option_id']}",
    )
    cpu_cores = int(st.number_input("CPU threads", min_value=1, max_value=256, value=8, step=1, key="md-analysis-cpu-cores"))
    if pbc_mode == "auto":
        st.caption("If the source has a periodic cell but its topology cannot safely unwrap protein/ligand molecules, the workflow fails instead of silently analyzing broken coordinates.")

    placement_ok = True
    settings: dict[str, Any] = {}
    try:
        from mn_molecular_design_studio.cluster import config as cluster_config
        if cluster_config.enabled():
            owner_node = str(row.get("_md_source_node") or cluster_config.identity())
            selected_placement = str(st.session_state.get("cluster-placement", cluster_config.load().get("default_placement", "auto")))
            allowed_placements = {"auto", owner_node} if is_remote else {"auto", "local", owner_node}
            if selected_placement not in allowed_placements:
                where = row.get("node_name") or "the source workstation" if is_remote else "this workstation"
                st.warning(f"The selected source data is on {where}. Change the workstation selector to Automatic or the source owner; the trajectory will not be copied to another machine.")
                placement_ok = False
            else:
                settings["_execution"] = {"placement": selected_placement}
    except (ImportError, OSError, ValueError):
        pass

    disabled = not ligand_selection or not protein_selection or not project_id or not placement_ok
    if st.button("Run MD analysis", type="primary", disabled=disabled, key=f"md-analysis-run-{selected['option_id']}"):
        try:
            if is_remote:
                source_run_parameter = str(source_run_dir)
                topology_parameter = str(topology_path)
                trajectory_parameter = str(trajectory_path)
            else:
                source_run_dir = Path(source_run_dir).resolve(strict=True)
                source_run_parameter = _portable_source_path(source_run_dir, source_run_dir, row)
                topology_parameter = _portable_source_path(topology_path, source_run_dir, row)
                trajectory_parameter = _portable_source_path(trajectory_path, source_run_dir, row)
            parameters = {
                "source_run_dir": source_run_parameter,
                "topology_path": topology_parameter,
                "trajectory_path": trajectory_parameter,
                "ligand_selection": ligand_selection,
                "protein_selection": protein_selection,
                "stride": stride,
                "start_ns": start_ns,
                "end_ns": end_ns,
                "contact_cutoff_angstrom": cutoff,
                "pbc_mode": pbc_mode,
                "chunk_size": chunk_size,
                "cpu_cores": cpu_cores,
                "source_segments": [{
                    **segment,
                    "study_id": study.get("study_id"),
                    # These are provenance references only. The operational
                    # source paths above carry owner affinity and are validated
                    # on the owning workstation.
                    "trajectory_ref": relative_trajectory,
                    "topology_ref": relative_topology,
                    "result_ref": selected.get("source_result_path") or "",
                }],
            }
            submit(schema, workflow_id, parameters, str(project_id), settings, title)
        except (OSError, ValueError, RuntimeError) as exc:
            st.error(f"Could not start MD analysis: {exc}")
