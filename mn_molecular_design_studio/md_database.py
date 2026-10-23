"""Idempotent Molecular Database observations from MD study summaries.

Only observations are written. Molecule identity must already exist in an
authorized collection; MD files and structures remain in the Results run.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import socket
from pathlib import Path, PurePosixPath
from typing import Any

from mn_molecular_design_studio.molecular_database import MolecularDatabase


_ENERGY_TERMS = {
    "g_bind_total": "delta_g_bind_kcal_mol",
    "g_bind": "delta_g_bind_kcal_mol",
    "total": "delta_g_bind_kcal_mol",
    "mm": "delta_mm_kcal_mol",
    "mechanical": "delta_mm_kcal_mol",
    "gbsa": "delta_solvation_kcal_mol",
    "pbsa": "delta_solvation_kcal_mol",
    "solvation": "delta_solvation_kcal_mol",
    "polar_solvation": "delta_solvation_kcal_mol",
    "nonpolar": "delta_nonpolar_kcal_mol",
    "nonpolar_solvation": "delta_nonpolar_kcal_mol",
}
_ANALYSIS_SCALARS = {
    "backbone_rmsd_tail_mean_angstrom",
    "backbone_rmsd_mean_angstrom",
    "ligand_rmsd_tail_mean_angstrom",
    "ligand_rmsd_mean_angstrom",
    "backbone_rmsd_max_angstrom",
    "ligand_rmsd_max_angstrom",
    "reference_site_retained_fraction",
    "ligand_centroid_displacement_max_angstrom",
    "performance_ns_per_day",
    "selected_frame_count",
    "trajectory_atom_count",
    "trajectory_frames_seen_after_stride",
    "selected_duration_ns",
    "contact_count_mean",
}
_ANALYSIS_UNITS = {
    "backbone_rmsd_tail_mean_angstrom": "angstrom",
    "backbone_rmsd_mean_angstrom": "angstrom",
    "ligand_rmsd_tail_mean_angstrom": "angstrom",
    "ligand_rmsd_mean_angstrom": "angstrom",
    "backbone_rmsd_max_angstrom": "angstrom",
    "ligand_rmsd_max_angstrom": "angstrom",
    "reference_site_retained_fraction": "fraction",
    "ligand_centroid_displacement_max_angstrom": "angstrom",
    "performance_ns_per_day": "ns/day",
    "selected_frame_count": "count",
    "trajectory_atom_count": "count",
    "trajectory_frames_seen_after_stride": "count",
    "selected_duration_ns": "ns",
    "contact_count_mean": "count",
}
_OCCUPANCY_FIELDS = (
    "contact_occupancy",
    "contact_backbone_occupancy",
    "contact_sidechain_occupancy",
    "hydrogen_bond_occupancy",
    "hydrogen_bond_backbone_occupancy",
    "hydrogen_bond_sidechain_occupancy",
    "hydrophobic_occupancy",
    "water_bridge_occupancy",
    "salt_bridge_occupancy",
)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"Could not read JSON file: {path.name}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON file is not an object: {path.name}")
    return value


def _safe_path(root: Path, value: Any) -> Path:
    raw = str(value or "")
    logical = PurePosixPath(raw)
    if not raw or "\\" in raw or "\x00" in raw or logical.is_absolute() or any(part in {"", ".", ".."} for part in logical.parts):
        raise ValueError("MD result path must be a safe relative path.")
    try:
        path = (root / Path(*logical.parts)).resolve(strict=True)
    except RuntimeError as exc:
        raise ValueError("MD result path contains a symlink loop.") from exc
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValueError("MD result path resolves outside the source Results run.") from exc
    if not path.is_file():
        raise ValueError("MD result path does not identify a regular file.")
    return path


def _load_run_identity(root: Path, metadata: dict[str, Any]) -> tuple[str, str]:
    marker = {}
    marker_path = root / "molecular_design_studio.json"
    if marker_path.is_file():
        try:
            marker = _read_json(marker_path)
        except ValueError:
            marker = {}
    embedded = metadata.get("molecular_design_studio")
    embedded = embedded if isinstance(embedded, dict) else {}
    owner = str(marker.get("studio_owner_id") or metadata.get("owner_id") or embedded.get("studio_owner_id") or "")
    project = str(marker.get("project_id") or metadata.get("project_id") or embedded.get("project_id") or "")
    return owner, project


def _number(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _summary_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _member_identity(database: MolecularDatabase, group_id: str, molecule_id: str) -> tuple[dict[str, Any] | None, str | None]:
    with database.connection() as connection:
        row = connection.execute(
            "SELECT mm.components,mm.metadata FROM molecule_members AS mm WHERE mm.group_id=? AND mm.molecule_id=?",
            (group_id, molecule_id),
        ).fetchone()
    if not row:
        return None, "molecule is not a member of the stated group"
    try:
        components = json.loads(row["components"])
        member_metadata = json.loads(row["metadata"])
    except (TypeError, ValueError):
        return None, "stored molecule identity is unreadable"
    return {"components": components, "metadata": member_metadata}, None


def _db_identity(database: MolecularDatabase, group_id: str, source_group_id: str, molecule_id: str, user: dict[str, Any], identity_smiles: Any, compound_id: Any) -> tuple[dict[str, Any] | None, str | None]:
    try:
        group = database.group(group_id, user, write=True)
    except (KeyError, PermissionError, ValueError) as exc:
        return None, f"molecule group is missing or inaccessible: {exc}"
    member, error = _member_identity(database, group_id, molecule_id)
    if error:
        return None, error
    if source_group_id and source_group_id != group_id:
        try:
            database.group(source_group_id, user)
        except (KeyError, PermissionError, ValueError) as exc:
            return None, f"source dataset group is missing or inaccessible: {exc}"
        source_member, error = _member_identity(database, source_group_id, molecule_id)
        if error:
            return None, f"source dataset membership could not be verified: {error}"
    else:
        source_member = member
    components = member["components"]
    member_metadata = member["metadata"]
    known_compound_id = member_metadata.get("compound_id") if isinstance(member_metadata, dict) else None
    if known_compound_id not in (None, "") and str(known_compound_id) != str(compound_id or ""):
        return None, "compound ID conflicts with the explicitly referenced Database molecule"
    if identity_smiles:
        try:
            from mn_molecular_design_studio.molecule_imports import canonical_component

            source = canonical_component({"type": "small_molecule", "chain": "L", "smiles": identity_smiles}).get("smiles")
            known = {
                canonical_component({"type": "small_molecule", "chain": item.get("chain") or "L", "smiles": item.get("smiles")}).get("smiles")
                for item in components if isinstance(item, dict) and item.get("type") == "small_molecule" and item.get("smiles")
            }
            if not source or source not in known:
                return None, "ligand SMILES conflicts with the explicitly referenced Database molecule"
        except (ImportError, ValueError, TypeError):
            return None, "ligand SMILES could not be verified against the Database molecule"
    # If the source group differs, it must resolve to the exact same existing molecule.
    source_known_compound = (source_member["metadata"] or {}).get("compound_id")
    if source_known_compound not in (None, "") and str(source_known_compound) != str(compound_id or ""):
        return None, "compound ID conflicts with the source dataset membership"
    return group, None


def _identity(study: dict[str, Any], database: MolecularDatabase, user: dict[str, Any]) -> tuple[dict[str, Any] | None, str | None]:
    if not study.get("study_id") or not study.get("compound_id"):
        return None, "explicit study_id and compound_id are required"
    group_id = str(study.get("group_id") or "")
    molecule_id = str(study.get("molecule_id") or "")
    if not group_id or not molecule_id:
        return None, "explicit group_id and molecule_id are required"
    if study.get("identity_status") != "verified_source_dataset":
        return None, "identity_status is not verified_source_dataset"
    identity_smiles = str(study.get("identity_smiles") or "").strip()
    if not identity_smiles:
        return None, "identity_smiles from the source dataset is required to verify molecule identity"
    source_group_id = str(study.get("source_dataset_group_id") or group_id)
    group, error = _db_identity(database, group_id, source_group_id, molecule_id, user, identity_smiles, study.get("compound_id"))
    if error:
        return None, error
    return {"group_id": group_id, "source_dataset_group_id": source_group_id, "molecule_id": molecule_id, "group": group}, None


def _unit(value: Any) -> str | None:
    normalized = re.sub(r"\s+", "", str(value or "").casefold()).replace("·", "")
    normalized = normalized.replace("per", "/").replace("mol^-1", "/mol").replace("mol-1", "/mol")
    if normalized in {"kcal/mol", "kcalmol", "kcal"}:
        return "kcal"
    if normalized in {"kj/mol", "kjmol", "kj"}:
        return "kj"
    return None


def _energy_term(key: str) -> tuple[str | None, str | None]:
    name = str(key).casefold()
    unit = None
    for suffix, candidate in (("_kcal_mol", "kcal"), ("_kj_mol", "kj"), ("_kcalmol", "kcal"), ("_kjmol", "kj")):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            unit = candidate
            break
    name = name.removeprefix("delta_").replace("-", "_")
    return _ENERGY_TERMS.get(name), unit


def _normalize_endpoint_branch(mmgbsa: dict[str, Any], branch_name: str, branch: dict[str, Any]) -> tuple[dict[str, float] | None, str | None]:
    units = [value for value in (mmgbsa.get("units"), mmgbsa.get("units_secondary")) if value not in (None, "")]
    declared = [_unit(value) for value in units]
    if any(value is None for value in declared):
        return None, "endpoint declares unknown energy units"
    default_unit = _unit(mmgbsa.get("units"))
    candidates: dict[str, list[tuple[str, float]]] = {}
    for key, raw in branch.items():
        metric, explicit_unit = _energy_term(str(key))
        if not metric:
            continue
        value = _number(raw)
        if value is None:
            continue
        effective_unit = explicit_unit or default_unit
        if effective_unit is None:
            return None, f"ambiguous energy units for {branch_name}.{key}"
        normalized = value if effective_unit == "kcal" else value / 4.184
        candidates.setdefault(metric, []).append((effective_unit, normalized))
    result = {}
    for metric, values in candidates.items():
        first = values[0][1]
        if any(not math.isclose(first, value, rel_tol=1e-4, abs_tol=0.02) for _, value in values[1:]):
            return None, f"inconsistent kcal/kJ duplicate values for {branch_name}.{metric}"
        result[metric] = first
    if not result:
        return None, "endpoint contains no supported finite energy terms"
    return result, None


def _endpoint_observations(endpoint: dict[str, Any], method_hint: str = "") -> tuple[list[tuple[str, dict[str, float]]], str | None]:
    mmgbsa = endpoint.get("mmgbsa")
    if not isinstance(mmgbsa, dict) or endpoint.get("success") is not True or mmgbsa.get("status") != "success":
        return [], "endpoint result is not successful"
    branches = []
    for branch_name in ("gb", "pb"):
        branch = mmgbsa.get(branch_name)
        if isinstance(branch, dict) and isinstance(branch.get("delta"), dict):
            normalized, error = _normalize_endpoint_branch(mmgbsa, branch_name, branch["delta"])
            if error:
                return [], error
            branches.append((f"md.mmgbsa.{branch_name}", normalized))
    if not branches and isinstance(mmgbsa.get("delta"), dict):
        method = str(method_hint or mmgbsa.get("method") or "").casefold()
        if "pb" in method or not ("mmgbsa" in method or "gbsa" in method):
            return [], "single endpoint delta has ambiguous or PB method"
        normalized, error = _normalize_endpoint_branch(mmgbsa, "delta", mmgbsa["delta"])
        if error:
            return [], error
        branches.append(("md.mmgbsa.gb", normalized))
    if not branches:
        return [], "endpoint lacks distinct GB/PB branches"
    return branches, None


def _analysis_metrics(replica: dict[str, Any]) -> dict[str, float]:
    return {key: value for key in _ANALYSIS_SCALARS if (value := _number(replica.get(key))) is not None}


def _direct_summary_metrics(summary: dict[str, Any]) -> dict[str, float]:
    metrics: dict[str, float] = {}
    source_metrics = summary.get("metrics") if isinstance(summary.get("metrics"), dict) else {}
    for key in _ANALYSIS_SCALARS:
        raw = source_metrics.get(key, summary.get(key))
        if isinstance(raw, dict):
            unit = str(raw.get("unit") or raw.get("units") or "").casefold().replace("å", "angstrom")
            expected = _ANALYSIS_UNITS.get(key, "")
            if unit and unit != expected:
                continue
            raw = raw.get("value")
        value = _number(raw)
        if value is not None:
            metrics[key] = value
    time_range = summary.get("selected_time_range_ns")
    if isinstance(time_range, dict):
        start, end = _number(time_range.get("start")), _number(time_range.get("end"))
        if start is not None and end is not None and end >= start:
            metrics["selected_duration_ns"] = end - start
    return metrics


def _target_identifier(study: dict[str, Any], provenance: dict[str, Any]) -> str:
    subtype = str(study.get("target_subtype") or study.get("target_receptor_subtype") or study.get("target_receptor") or "target")
    pdb_id = str(study.get("target_pdb_id") or "unknown-pdb")
    prepared_id = str(provenance.get("prepared_target_run_id") or provenance.get("prepared_target_id") or "")
    digest = str(provenance.get("prepared_target_sha256") or provenance.get("sha256") or "")
    prepared = prepared_id or (digest[:12] if digest else "unversioned")
    return f"{subtype}:{pdb_id}:{prepared}"


def _base_context(root: Path, study: dict[str, Any], identity: dict[str, Any], *, source_hash: str, relative_path: str, current_run_id: str, owner_id: str, run_metadata: dict[str, Any]) -> dict[str, Any]:
    target = study.get("target_provenance") if isinstance(study.get("target_provenance"), dict) else {}
    source_metadata = study.get("source_metadata") if isinstance(study.get("source_metadata"), dict) else {}
    return {
        "type": "job",
        "run_id": current_run_id,
        "task_group": root.parent.name,
        "run_dir": str(root),
        "tool": str(run_metadata.get("tool") or source_metadata.get("tool") or study.get("tool") or "md_study"),
        "workstation": socket.gethostname(),
        "target_id": _target_identifier(study, target),
        "target_label": _target_identifier(study, target),
        "source_owner_id": owner_id,
        "source_run_id": current_run_id,
        "source_task_group": root.parent.name,
        "source_workstation": str(study.get("source_workstation") or source_metadata.get("source_workstation") or socket.gethostname()),
        "source_metadata": source_metadata,
        "study_id": str(study.get("study_id") or ""),
        "compound_id": str(study.get("compound_id") or ""),
        "dataset_name": str(study.get("dataset_name") or ""),
        "group_id": identity["group_id"],
        "source_dataset_group_id": identity["source_dataset_group_id"],
        "molecule_id": identity["molecule_id"],
        "identity_status": "verified_source_dataset",
        "identity_smiles": str(study.get("identity_smiles") or ""),
        "ligand_smiles": str(study.get("ligand_smiles") or ""),
        "identity_match_basis": "source_dataset_identity_smiles_and_existing_membership",
        "target_pdb_id": str(study.get("target_pdb_id") or ""),
        "target_receptor": study.get("target_receptor"),
        "target_provenance": target,
        "source_pose_id": str(study.get("source_pose_id") or ""),
        "source_tool": str(study.get("source_tool") or ""),
        "source_pose_rank": study.get("source_pose_rank"),
        "selection_status": str(study.get("selection_status") or ""),
        "selected_complex_run_id": str(study.get("selected_complex_run_id") or ""),
        "source_summary_sha256": source_hash,
        "source_summary_relative_path": relative_path,
        "is_representative": False,
    }


def _write_observation(database: MolecularDatabase, identity: dict[str, Any], namespace: str, run_id: str, source_key: str, metrics: dict[str, Any], context: dict[str, Any]) -> str:
    from mn_molecular_design_studio.molecular_database import _now

    created_at = context.get("source_completed_at") or context.get("source_created_at") or context.get("owning_job_completed_at") or context.get("owning_job_created_at") or _now()
    return database.add_observation(
        identity["group_id"], identity["molecule_id"], namespace, run_id, metrics,
        context=context, source_key=source_key, created_at=str(created_at),
    )


def _source_dates(entry: dict[str, Any], run_metadata: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_created_at": entry.get("source_created_at"),
        "source_completed_at": entry.get("source_completed_at"),
        "owning_job_created_at": run_metadata.get("created_at"),
        "owning_job_completed_at": run_metadata.get("completed_at"),
    }


def _production_details(root: Path, entry: dict[str, Any]) -> dict[str, Any]:
    details = {
        "status": str(entry.get("status") or "unknown"),
        "actual_duration_ns": _number(entry.get("actual_duration_ns") or entry.get("duration_ns") or entry.get("segment_duration_ns")),
        "duration_basis": str(entry.get("segment_duration_basis") or "recorded production duration"),
        "trajectory_frames": _number(entry.get("trajectory_frames")),
    }
    relative = entry.get("result_path")
    if not relative:
        return details
    try:
        payload = _read_json(_safe_path(root, relative))
    except (OSError, ValueError):
        return details
    md_result = payload.get("md_result") if isinstance(payload.get("md_result"), dict) else {}
    stats_root = md_result.get("equilibration_stats") if isinstance(md_result.get("equilibration_stats"), dict) else {}
    production = stats_root.get("production") if isinstance(stats_root.get("production"), dict) else {}
    if payload.get("success") is True and str(md_result.get("status") or "").casefold() not in {"failed", "cancelled", "canceled", "error"}:
        details["actual_duration_ns"] = _number(production.get("duration_ns")) or details["actual_duration_ns"]
        details["trajectory_frames"] = _number(production.get("trajectory_frames")) or details["trajectory_frames"]
    return details


def _ingest_study(root: Path, study: dict[str, Any], database: MolecularDatabase, user: dict[str, Any], current_run_id: str, owner_id: str, run_metadata: dict[str, Any]) -> dict[str, Any]:
    study_id = str(study.get("study_id") or "")
    identity, error = _identity(study, database, user)
    if error:
        return {"study_id": study_id, "status": "withheld", "reason": error, "observations": 0}
    identity_ok = identity
    report = {"study_id": study_id, "status": "imported", "observations": 0, "withheld": []}

    for endpoint_entry in study.get("mmgbsa_runs") or []:
        if not isinstance(endpoint_entry, dict):
            continue
        if not endpoint_entry.get("result_path"):
            continue
        endpoint_run_id = str(endpoint_entry.get("run_id") or "")
        relative = str(endpoint_entry.get("result_path") or "")
        if endpoint_entry.get("status") not in {"completed", "success"}:
            report["withheld"].append({"run_id": endpoint_run_id, "kind": "mmgbsa", "reason": "manifest endpoint status is not completed"})
            continue
        try:
            path = _safe_path(root, relative)
            endpoint = _read_json(path)
            method_hint = str(endpoint_entry.get("method") or "")
            observations, error = _endpoint_observations(endpoint, method_hint)
            if error:
                raise ValueError(error)
            digest = _summary_hash(path)
        except (OSError, ValueError) as exc:
            report["withheld"].append({"run_id": endpoint_run_id, "kind": "mmgbsa", "reason": str(exc)})
            continue
        parameters = endpoint_entry.get("parameters") if isinstance(endpoint_entry.get("parameters"), dict) else {}
        context = _base_context(root, study, identity_ok, source_hash=digest, relative_path=relative, current_run_id=current_run_id, owner_id=owner_id, run_metadata=run_metadata)
        production_id = str(endpoint_entry.get("source_production_run_id") or "")
        production = next((item for item in study.get("production_runs") or [] if isinstance(item, dict) and str(item.get("run_id") or "") == production_id), {})
        production_details = _production_details(root, production) if production else {}
        previous_id = str(production.get("continuation_of_run_id") or "")
        previous = next((item for item in study.get("production_runs") or [] if isinstance(item, dict) and str(item.get("run_id") or "") == previous_id), {})
        previous_details = _production_details(root, previous) if previous else {}
        method = str(endpoint_entry.get("method") or endpoint.get("mmgbsa", {}).get("method") or "")
        target_duration = _number(production.get("target_duration_ns"))
        context.update({
            "endpoint_run_id": endpoint_run_id,
            "source_production_run_id": production_id,
            "md_engine": str(study.get("md_engine") or ""),
            "backend": str(endpoint_entry.get("backend") or parameters.get("backend") or ""),
            "method": method,
            "source_units": endpoint.get("mmgbsa", {}).get("units"),
            "source_units_secondary": endpoint.get("mmgbsa", {}).get("units_secondary"),
            "sampling_parameters": parameters,
            "recorded_trajectory_sha256": str(endpoint_entry.get("recorded_trajectory_sha256") or ""),
            "copied_trajectory_sha256": str(endpoint_entry.get("copied_trajectory_sha256") or ""),
            "trajectory_snapshot_matches": endpoint_entry.get("trajectory_snapshot_matches"),
            "source_md_job": endpoint_entry.get("source_md_job") if isinstance(endpoint_entry.get("source_md_job"), dict) else {},
            "legacy_source_md_job": endpoint_entry.get("legacy_source_md_job") if isinstance(endpoint_entry.get("legacy_source_md_job"), dict) else {},
            "source_trajectory_ref": str(endpoint_entry.get("source_trajectory_ref") or ""),
            "replica_index": production.get("repeat_index"),
            "continuation_of_run_id": str(production.get("continuation_of_run_id") or ""),
            "source_production_status": production_details.get("status", "unknown"),
            "source_production_segment_duration_ns": production_details.get("actual_duration_ns"),
            "source_production_segment_duration_basis": production_details.get("duration_basis"),
            "source_production_previous_duration_ns": previous_details.get("actual_duration_ns"),
            "source_production_requested_target_duration_ns": target_duration,
            "endpoint_superseded_by_run_id": str(endpoint_entry.get("superseded_by_run_id") or ""),
            "is_representative": not bool(endpoint_entry.get("superseded_by_run_id")),
            **_source_dates(endpoint_entry, run_metadata),
        })
        if endpoint_entry.get("trajectory_snapshot_matches") is False:
            context["trajectory_lineage_warning"] = "The copied MD trajectory bytes differ from the bytes recorded when this endpoint result was produced."
        context["settings"] = {
            "method": method,
            "sampling_parameters": parameters,
            "md_engine": str(study.get("md_engine") or ""),
            "segment_duration_ns": production_details.get("actual_duration_ns"),
            "previous_duration_ns": previous_details.get("actual_duration_ns"),
            "target_duration_ns": target_duration,
        }
        method = context["method"] or context["backend"] or "unspecified"
        for namespace, metrics in observations:
            source_key = ":".join(("md-study", owner_id, current_run_id, study_id, identity_ok["group_id"], identity_ok["molecule_id"], endpoint_run_id, context["backend"], method, namespace))
            _write_observation(database, identity_ok, namespace, current_run_id, source_key, metrics, context)
            report["observations"] += 1

    for analysis_entry in study.get("analysis_runs") or []:
        if not isinstance(analysis_entry, dict):
            continue
        if not analysis_entry.get("result_path"):
            continue
        relative = str(analysis_entry.get("result_path") or "")
        analysis_id = str(analysis_entry.get("run_id") or "")
        if analysis_entry.get("status") not in {"completed", "success"}:
            continue
        try:
            path = _safe_path(root, relative)
            analysis = _read_json(path)
            digest = _summary_hash(path)
        except (OSError, ValueError) as exc:
            report["withheld"].append({"run_id": analysis_id, "kind": "analysis", "reason": str(exc)})
            continue
        replica_rows = analysis.get("replicas") if isinstance(analysis.get("replicas"), list) else []
        is_direct_summary = bool(
            analysis.get("source_study_id")
            or analysis_entry.get("result_type") == "md_trajectory_analysis_result/v1"
            or analysis_entry.get("type") == "md_trajectory_analysis_result/v1"
            or analysis.get("result_type") == "md_trajectory_analysis_result/v1"
        )
        if analysis.get("success") is not True and not is_direct_summary:
            continue
        if is_direct_summary:
            if analysis.get("success") is not True:
                report["withheld"].append({"run_id": analysis_id, "kind": "analysis", "reason": "MD trajectory analysis summary did not report success"})
                continue
            if str(analysis.get("source_study_id") or "") != study_id:
                report["withheld"].append({"run_id": analysis_id, "kind": "analysis", "reason": "MD trajectory analysis summary source_study_id does not match this study"})
                continue
            metrics = _direct_summary_metrics(analysis)
            if not metrics:
                report["withheld"].append({"run_id": analysis_id, "kind": "analysis", "reason": "MD trajectory analysis summary has no supported finite scalar metrics"})
                continue
            context = _base_context(root, study, identity_ok, source_hash=digest, relative_path=relative, current_run_id=current_run_id, owner_id=owner_id, run_metadata=run_metadata)
            source = analysis.get("source") if isinstance(analysis.get("source"), dict) else {}
            topology = source.get("topology") if isinstance(source.get("topology"), dict) else {}
            trajectory = source.get("trajectory") if isinstance(source.get("trajectory"), dict) else {}
            parameters = analysis.get("parameters") if isinstance(analysis.get("parameters"), dict) else {}
            context.update({
                "analysis_run_id": analysis_id,
                "analysis_result_type": "md_trajectory_analysis_result/v1",
                "analysis_engine": str(analysis_entry.get("engine") or analysis.get("engine") or "MDTraj"),
                "source_study_id": str(analysis.get("source_study_id") or ""),
                "source_segments": analysis.get("source_segments") if isinstance(analysis.get("source_segments"), list) else [],
                "source_topology": topology,
                "source_trajectory": trajectory,
                "topology_sha256": str(topology.get("sha256") or ""),
                "trajectory_sha256": str(trajectory.get("sha256") or ""),
                "selected_time_range_ns": analysis.get("selected_time_range_ns"),
                "analysis_parameters": parameters,
                "analysis_selection": analysis.get("selection") if isinstance(analysis.get("selection"), dict) else {},
                "analysis_pbc": {key: value for key, value in parameters.items() if str(key).startswith("pbc_")},
                "settings": {
                    "analysis_engine": str(analysis_entry.get("engine") or analysis.get("engine") or "MDTraj"),
                    "analysis_parameters": parameters,
                    "selection": analysis.get("selection") if isinstance(analysis.get("selection"), dict) else {},
                    "time_window_ns": analysis.get("selected_time_range_ns"),
                },
                **_source_dates(analysis_entry, run_metadata),
            })
            artifacts = analysis.get("artifacts") if isinstance(analysis.get("artifacts"), dict) else {}
            contact_relative = str(artifacts.get("residue_contacts_csv") or "")
            if contact_relative:
                try:
                    contact_path = _safe_path(root, contact_relative)
                    context["contact_map_artifact_relative_path"] = contact_relative
                    context["contact_map_artifact_sha256"] = _summary_hash(contact_path)
                except (OSError, ValueError):
                    context["contact_map_artifact_unavailable"] = True
            source_key = ":".join(("md-study", owner_id, current_run_id, study_id, identity_ok["group_id"], identity_ok["molecule_id"], analysis_id, "analysis-summary"))
            _write_observation(database, identity_ok, "md.analysis", current_run_id, source_key, metrics, context)
            report["observations"] += 1
            continue
        for replica in replica_rows:
            if not isinstance(replica, dict):
                continue
            if replica.get("status") not in {"completed", "success"}:
                continue
            metrics = _analysis_metrics(replica)
            if not metrics:
                continue
            replica_run_id = str(replica.get("run_id") or "")
            production = next((item for item in study.get("production_runs") or [] if isinstance(item, dict) and str(item.get("run_id") or "") == replica_run_id), {})
            context = _base_context(root, study, identity_ok, source_hash=digest, relative_path=relative, current_run_id=current_run_id, owner_id=owner_id, run_metadata=run_metadata)
            context.update({
                "analysis_run_id": analysis_id,
                "source_production_run_id": replica_run_id,
                "replica_index": replica.get("replica", production.get("repeat_index")),
                "production_status": str(production.get("status") or replica.get("status") or ""),
                "continuation_of_run_id": str(production.get("continuation_of_run_id") or ""),
                "requested_target_duration_ns": _number(production.get("target_duration_ns")),
                **_source_dates(analysis_entry, run_metadata),
                "settings": {
                    "analysis_parameters": analysis.get("parameters") if isinstance(analysis.get("parameters"), dict) else {},
                    "time_window_ns": analysis.get("selected_time_range_ns"),
                    "replica_index": replica.get("replica", production.get("repeat_index")),
                    "md_engine": str(study.get("md_engine") or ""),
                },
            })
            if isinstance(analysis.get("contact_matrix"), dict):
                context["contact_map_summary_relative_path"] = relative
                context["contact_map_summary_sha256"] = digest
                context["contact_map_source"] = "inline contact_matrix in analysis summary; full residue data remains in Results artifact"
            source_key = ":".join(("md-study", owner_id, current_run_id, study_id, identity_ok["group_id"], identity_ok["molecule_id"], analysis_id, replica_run_id, "analysis"))
            _write_observation(database, identity_ok, "md.analysis", current_run_id, source_key, metrics, context)
            report["observations"] += 1

    for production_entry in study.get("production_runs") or []:
        if not isinstance(production_entry, dict):
            continue
        relative = str(production_entry.get("result_path") or "")
        if not relative:
            continue
        production_id = str(production_entry.get("run_id") or "")
        try:
            path = _safe_path(root, relative)
            production_result = _read_json(path)
            if production_entry.get("status") not in {"completed", "success"} or production_result.get("success") is not True:
                continue
            md_result = production_result.get("md_result") if isinstance(production_result.get("md_result"), dict) else {}
            if str(md_result.get("status") or "").casefold() in {"failed", "cancelled", "canceled", "error"}:
                continue
            equilibration = md_result.get("equilibration_stats") if isinstance(md_result.get("equilibration_stats"), dict) else {}
            production_stats = equilibration.get("production") if isinstance(equilibration.get("production"), dict) else {}
            digest = _summary_hash(path)
        except (OSError, ValueError) as exc:
            report["withheld"].append({"run_id": production_id, "kind": "production", "reason": str(exc)})
            continue
        metrics = {}
        duration = _number(production_stats.get("duration_ns"))
        frames = _number(production_stats.get("trajectory_frames"))
        if duration is not None:
            metrics["production_actual_duration_ns"] = duration
        elif _number(production_entry.get("segment_duration_ns")) is not None:
            metrics["production_recorded_segment_duration_ns"] = _number(production_entry["segment_duration_ns"])
        if frames is not None:
            metrics["production_trajectory_frames"] = frames
        if not metrics:
            continue
        context = _base_context(root, study, identity_ok, source_hash=digest, relative_path=relative, current_run_id=current_run_id, owner_id=owner_id, run_metadata=run_metadata)
        context.update({
            "source_production_run_id": production_id,
            "replica_index": production_entry.get("repeat_index"),
            "continuation_of_run_id": str(production_entry.get("continuation_of_run_id") or ""),
            "requested_target_duration_ns": _number(production_entry.get("target_duration_ns")),
            "production_status": str(production_entry.get("status") or ""),
            "segment_duration_basis": str(production_entry.get("segment_duration_basis") or "recorded production duration"),
            **_source_dates(production_entry, run_metadata),
            "settings": {
                "md_engine": str(study.get("md_engine") or ""),
                "segment_duration_ns": duration if duration is not None else _number(production_entry.get("segment_duration_ns")),
                "requested_target_duration_ns": _number(production_entry.get("target_duration_ns")),
                "replica_index": production_entry.get("repeat_index"),
                "continuation_of_run_id": str(production_entry.get("continuation_of_run_id") or ""),
            },
        })
        source_key = ":".join(("md-study", owner_id, current_run_id, study_id, identity_ok["group_id"], identity_ok["molecule_id"], production_id, "production_stats"))
        _write_observation(database, identity_ok, "md.production", current_run_id, source_key, metrics, context)
        report["observations"] += 1

    if report["withheld"] and not report["observations"]:
        report["status"] = "withheld"
    elif report["withheld"]:
        report["status"] = "partial"
    return report


def ingest_md_studies(run_dir: str | Path, *, database: MolecularDatabase | None = None, user: dict[str, Any] | None = None) -> dict[str, Any]:
    """Ingest validated MD study metrics; never create or infer molecule IDs."""
    root = Path(run_dir).expanduser().resolve(strict=True)
    if not root.is_dir():
        raise ValueError("MD Results run must be a directory.")
    database = database or MolecularDatabase()
    manifest_path = root / "artifacts" / "md" / "study.json"
    try:
        resolved_manifest = manifest_path.resolve(strict=True)
        resolved_manifest.relative_to(root)
        manifest_bytes = resolved_manifest.read_bytes()
        manifest = json.loads(manifest_bytes)
    except FileNotFoundError:
        return {"status": "no_manifest", "studies": [], "observations": 0}
    except (OSError, RuntimeError, ValueError) as exc:
        return {"status": "error", "reason": f"Invalid MD study manifest: {exc}", "studies": [], "observations": 0}
    if not isinstance(manifest, dict) or manifest.get("schema_version") != 1 or not isinstance(manifest.get("studies"), list):
        return {"status": "error", "reason": "Unsupported MD study manifest schema.", "studies": [], "observations": 0}
    metadata_path = root / "metadata.json"
    try:
        metadata = _read_json(metadata_path)
    except ValueError:
        metadata = {}
    owner, project = _load_run_identity(root, metadata)
    if not owner:
        return {"status": "error", "reason": "Trusted Studio owner identity is missing.", "studies": [], "observations": 0}
    if user is None:
        user = {"id": owner, "role": "admin", "project_id": project}
    elif not user.get("id"):
        return {"status": "error", "reason": "Database user identity is missing.", "studies": [], "observations": 0}
    records = []
    for study in manifest["studies"]:
        if not isinstance(study, dict):
            records.append({"study_id": "", "status": "withheld", "reason": "study record is not an object", "observations": 0})
            continue
        records.append(_ingest_study(root, study, database, user, root.name, owner, metadata))
    return {"status": "ok", "studies": records, "observations": sum(row.get("observations", 0) for row in records)}


def ingest_completed_md_run(run_dir, *, database=None, force=False):
    """Apply the normal ingestion cache to curated MD summaries, never bulk data."""
    from mn_molecular_design_studio.molecular_database import _json, _now
    root = Path(run_dir).resolve(strict=True)
    database = database or MolecularDatabase()
    metadata = _read_json(root / 'metadata.json')
    result = _read_json(root / 'result.json')
    if metadata.get('status') not in {'completed', 'succeeded'} or result.get('success') is not True:
        return 0
    manifest = _read_json(_safe_path(root, 'artifacts/md/study.json'))
    paths = [root / name for name in ('metadata.json', 'result.json', 'molecular_design_studio.json', 'artifacts/md/study.json')]
    for study in manifest.get('studies', []):
        for field in ('production_runs', 'mmgbsa_runs', 'analysis_runs'):
            for entry in study.get(field) or []:
                if entry.get('result_path'):
                    try:
                        paths.append(_safe_path(root, entry['result_path']))
                    except (OSError, ValueError):
                        pass  # The observation adapter reports withheld unsafe/missing artifacts.
    signature = _json({'adapter_version': 2, 'files': [(str(p.relative_to(root)), p.stat().st_size, p.stat().st_mtime_ns)
                         for p in paths if p.is_file()]})
    key = f'{socket.gethostname()}/{root.parent.name}/{root.name}'
    with database.connection() as con:
        existing = con.execute('SELECT signature FROM molecular_ingestions WHERE source_key=?', (key,)).fetchone()
    if not force and existing and existing[0] == signature:
        return 0
    report = ingest_md_studies(root, database=database)
    if report.get('status') == 'error':
        raise ValueError(report.get('reason') or 'MD summary ingestion failed')
    temporary = root / '.md_ingestion_report.json.tmp'
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(root / 'md_ingestion_report.json')
    with database.connection() as con:
        con.execute('INSERT INTO molecular_ingestions VALUES (?,?,?) ON CONFLICT(source_key) DO UPDATE SET signature=excluded.signature,updated_at=excluded.updated_at',
                    (key, signature, _now()))
    return int(report.get('observations') or 0)
