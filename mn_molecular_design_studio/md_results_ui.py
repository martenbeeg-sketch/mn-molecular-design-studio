"""Results view for aggregate and legacy molecular-dynamics analyses."""
from __future__ import annotations

import json
import math
import csv
import io
from pathlib import Path, PurePosixPath
from typing import Any, Callable


def _run_root(run_dir: str | Path) -> Path:
    root = Path(run_dir).expanduser().resolve(strict=True)
    if not root.is_dir():
        raise ValueError("Results run path is not a directory.")
    return root


def _safe_artifact_path(run_root: Path, relative_path: Any) -> Path:
    """Resolve a manifest path only when its target remains inside this run."""
    logical = _safe_relative_parts(relative_path)
    path = (run_root / Path(*logical)).resolve(strict=True)
    try:
        path.relative_to(run_root)
    except ValueError as exc:
        raise ValueError("Artifact path resolves outside this Results run.") from exc
    if not path.is_file():
        raise ValueError("Artifact path does not name a file.")
    return path


def _safe_relative_parts(relative_path: Any) -> tuple[str, ...]:
    value = str(relative_path or "")
    if not value or "\\" in value or "\x00" in value:
        raise ValueError("Artifact path is empty or malformed.")
    logical = PurePosixPath(value)
    if logical.is_absolute() or any(part in {"", ".", ".."} for part in logical.parts):
        raise ValueError("Artifact path must be a safe relative path.")
    return logical.parts


def _read_json_artifact(run_root: Path, relative_path: Any, *, json_loader: Callable[[str, str], dict[str, Any] | None] | None = None, description: str = "JSON result") -> dict[str, Any]:
    _safe_relative_parts(relative_path)
    if json_loader is not None:
        payload = json_loader(str(relative_path), description)
        if payload is None:
            raise ValueError(f"{description} is still loading or unavailable.")
        if not isinstance(payload, dict):
            raise ValueError(f"{description} JSON is not an object.")
        return payload
    path = _safe_artifact_path(run_root, relative_path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"Could not read JSON result artifact: {relative_path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"JSON result artifact is not an object: {relative_path}")
    return payload


def load_md_study_manifest(run_dir: str | Path | None = None, *, preloaded: dict[str, Any] | None = None, json_loader: Callable[[str, str], dict[str, Any] | None] | None = None) -> dict[str, Any] | None:
    """Read the normalized study manifest; missing manifests mean legacy mode."""
    if preloaded is not None:
        payload = preloaded
        if not isinstance(payload, dict) or payload.get("schema_version") != 1 or not isinstance(payload.get("studies"), list) or any(not isinstance(item, dict) for item in payload["studies"]):
            raise ValueError("Unsupported MD study manifest schema.")
        return payload
    if json_loader is not None:
        payload = json_loader("artifacts/md/study.json", "MD study manifest")
        if payload is None:
            return None
        if not isinstance(payload, dict) or payload.get("schema_version") != 1:
            raise ValueError("Unsupported MD study manifest schema.")
        if not isinstance(payload.get("studies"), list) or any(not isinstance(item, dict) for item in payload["studies"]):
            raise ValueError("MD study manifest must contain a list of study records.")
        return payload
    if run_dir is None:
        raise ValueError("A local Results path is required without preloaded manifest data.")
    root = _run_root(run_dir)
    path = root / "artifacts" / "md" / "study.json"
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError:
        return None
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError("MD study manifest resolves outside this Results run.") from exc
    if not resolved.is_file():
        raise ValueError("MD study manifest is not a file.")
    try:
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError("Could not read MD study manifest JSON.") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("Unsupported MD study manifest schema.")
    studies = payload.get("studies")
    if not isinstance(studies, list) or any(not isinstance(item, dict) for item in studies):
        raise ValueError("MD study manifest must contain a list of study records.")
    return payload


def _finite(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        converted = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return converted if math.isfinite(converted) else None


def _finite_series(values: Any) -> list[float | None]:
    if not isinstance(values, list):
        return []
    return [_finite(value) for value in values]


def _finite_pairs(x_values: Any, y_values: Any) -> list[tuple[float, float]]:
    xs, ys = _finite_series(x_values), _finite_series(y_values)
    return [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]


def _study_label(study: dict[str, Any], index: int) -> str:
    study_id = str(study.get("study_id") or f"study-{index + 1}")
    compound = str(study.get("compound_name") or study.get("compound_id") or "compound unknown")
    target = str(study.get("target_pdb_id") or "target unknown")
    engine = str(study.get("md_engine") or "engine unknown")
    return f"{study_id} · {compound} · {target} · {engine}"


def _field_rows(study: dict[str, Any]) -> list[dict[str, str]]:
    return [
        {"Field": "Compound ID", "Value": str(study.get("compound_id") or "—")},
        {"Field": "Compound", "Value": str(study.get("compound_name") or "—")},
        {"Field": "Target PDB", "Value": str(study.get("target_pdb_id") or "—")},
        {"Field": "Target receptor", "Value": str(study.get("target_receptor") or "—")},
        {"Field": "MD engine", "Value": str(study.get("md_engine") or "—")},
        {"Field": "Source pose ID", "Value": str(study.get("source_pose_id") or "—")},
        {"Field": "Source tool", "Value": str(study.get("source_tool") or "—")},
        {"Field": "Source pose rank", "Value": str(study.get("source_pose_rank") if study.get("source_pose_rank") is not None else "—")},
        {"Field": "Selection", "Value": str(study.get("selection_status") or "—")},
        {"Field": "Selected complex run", "Value": str(study.get("selected_complex_run_id") or "—")},
    ]


def _rows(value: Any) -> list[dict[str, Any]]:
    return [row for row in value if isinstance(row, dict)] if isinstance(value, list) else []


def _safe_optional_json(st, run_root: Path, relative_path: Any, description: str, *, json_loader=None) -> dict[str, Any] | None:
    if not relative_path:
        st.warning(f"{description} result path is missing.")
        return None
    try:
        return _read_json_artifact(run_root, relative_path, json_loader=json_loader, description=description)
    except (OSError, ValueError) as exc:
        st.warning(f"{description} result is unavailable: {exc}")
        return None


def _line_chart(st, series: list[dict[str, Any]], *, label: str = "Replica") -> None:
    for index, row in enumerate(series):
        replica = row.get("replica", index + 1)
        time = _finite_series(row.get("time_ns"))
        backbone = _finite_series(row.get("backbone_rmsd_angstrom"))
        ligand = _finite_series(row.get("ligand_rmsd_angstrom"))
        count = min(len(time), len(backbone), len(ligand))
        points = []
        for frame in range(count):
            t, bb, lig = time[frame], backbone[frame], ligand[frame]
            if t is not None and (bb is not None or lig is not None):
                points.append({"time_ns": t, "backbone_rmsd_angstrom": bb, "ligand_rmsd_angstrom": lig})
        if points:
            st.caption(f"{label} {replica}: RMSD (Å) over time")
            st.line_chart(points, x="time_ns", y=["backbone_rmsd_angstrom", "ligand_rmsd_angstrom"])


def _render_structural_dynamics(st, structural: dict[str, Any], rmsd: dict[str, Any], *, label: str) -> None:
    """Show only stored OpenMM structural-dynamics outputs; perform no analysis."""
    time_ps = _finite_series(structural.get("time_ps"))
    time_ns = [value / 1000.0 if value is not None else None for value in time_ps]
    rmsd_time = _finite_series(rmsd.get("time_ps"))
    if not time_ns and rmsd_time:
        time_ns = [value / 1000.0 if value is not None else None for value in rmsd_time]
    series = [{
        "replica": label,
        "time_ns": time_ns,
        "backbone_rmsd_angstrom": rmsd.get("backbone_rmsd_angstrom"),
        "ligand_rmsd_angstrom": rmsd.get("ligand_rmsd_angstrom"),
    }]
    _line_chart(st, series, label="Production")

    rmsf = structural.get("rmsf") if isinstance(structural.get("rmsf"), dict) else {}
    residues = rmsf.get("residues")
    ca_values = rmsf.get("ca_rmsf_angstrom")
    if isinstance(residues, list) and isinstance(ca_values, list) and residues:
        st.markdown("**Protein Cα RMSF (Å)**")
        st.dataframe([
            {"Residue": str(residue), "Cα RMSF (Å)": _finite(value)}
            for residue, value in zip(residues, ca_values)
        ], use_container_width=True, hide_index=True)

    ligand_rmsf = structural.get("ligand_rmsf") if isinstance(structural.get("ligand_rmsf"), dict) else {}
    atoms = ligand_rmsf.get("atoms")
    atom_values = ligand_rmsf.get("rmsf_angstrom")
    if isinstance(atoms, list) and isinstance(atom_values, list) and atoms:
        st.markdown("**Ligand atom RMSF (Å)**")
        st.dataframe([
            {"Atom": str(atom), "RMSF (Å)": _finite(value)}
            for atom, value in zip(atoms, atom_values)
        ], use_container_width=True, hide_index=True)

    gyration = structural.get("radius_of_gyration")
    if isinstance(gyration, dict):
        columns = {key: _finite_series(gyration.get(key)) for key in (
            "protein_angstrom", "ligand_angstrom", "complex_angstrom",
        )}
        points = []
        for index, time in enumerate(time_ns):
            if time is None:
                continue
            row = {"time_ns": time}
            for key, values in columns.items():
                if index < len(values) and values[index] is not None:
                    row[key.removesuffix("_angstrom").replace("_", " ") + " Rg (Å)"] = values[index]
            if len(row) > 1:
                points.append(row)
        if points:
            st.markdown("**Radius of gyration (Å)**")
            st.line_chart(points, x="time_ns")

    secondary = structural.get("secondary_structure")
    if isinstance(secondary, dict):
        points = []
        for index, time in enumerate(time_ns):
            if time is None:
                continue
            row = {"time_ns": time}
            for key, values in secondary.items():
                if not isinstance(values, list) or index >= len(values):
                    continue
                value = _finite(values[index])
                if value is not None:
                    row[key.replace("_", " ")] = value
            if len(row) > 1:
                points.append(row)
        if points:
            st.markdown("**Secondary-structure fractions (recorded)**")
            st.line_chart(points, x="time_ns")

    pocket = structural.get("pocket") if isinstance(structural.get("pocket"), dict) else {}
    pocket_details = {key: pocket.get(key) for key in (
        "retained_fraction", "site_displacement_cutoff_angstrom", "contact_cutoff_angstrom",
    ) if key in pocket}
    if pocket_details:
        st.markdown("**Reference-site retention**")
        st.json(pocket_details)
    displacement = _finite_series(pocket.get("ligand_centroid_displacement_angstrom"))
    retained = pocket.get("reference_site_retained")
    if time_ns and (any(value is not None for value in displacement) or isinstance(retained, list)):
        points = []
        for index, time in enumerate(time_ns):
            if time is None:
                continue
            row = {"time_ns": time}
            if index < len(displacement) and displacement[index] is not None:
                row["centroid displacement (Å)"] = displacement[index]
            if isinstance(retained, list) and index < len(retained) and isinstance(retained[index], bool):
                row["reference site retained"] = float(retained[index])
            if len(row) > 1:
                points.append(row)
        if points:
            st.line_chart(points, x="time_ns")

    contacts = structural.get("contacts") if isinstance(structural.get("contacts"), dict) else {}
    contact_rows = _rows(contacts.get("residues"))
    hotspot = contacts.get("interaction_hotspot_score")
    if contact_rows:
        st.markdown("**Recorded protein–ligand contacts and interaction occupancy**")
        st.caption("Occupancies and scores are the values stored by the source analytics.")
        st.dataframe(contact_rows, use_container_width=True, hide_index=True)
    if isinstance(hotspot, dict) and hotspot:
        st.caption(str(hotspot.get("description") or "Recorded interaction score metadata"))
        st.json(hotspot)
    distance_series = contacts.get("distance_series")
    if isinstance(distance_series, dict) and distance_series:
        labels = list(distance_series)
        residue = st.selectbox(f"Recorded contact distance series · {label}", labels,
                               key=f"md-contact-distance-{label}")
        distances = _finite_series(distance_series.get(residue))
        points = [{"time_ns": time_ns[index], "minimum distance (Å)": distances[index]}
                  for index in range(min(len(time_ns), len(distances)))
                  if time_ns[index] is not None and distances[index] is not None]
        if points:
            st.line_chart(points, x="time_ns")

    for section, title in (("water_bridges", "Water bridges"), ("salt_bridges", "Salt bridges")):
        payload = structural.get(section)
        if not isinstance(payload, dict):
            continue
        st.markdown(f"**{title}**")
        meta = {key: value for key, value in payload.items() if key != "residues" and value not in (None, [], {})}
        if meta:
            st.json(meta)
        rows = _rows(payload.get("residues"))
        if rows:
            st.dataframe(rows, use_container_width=True, hide_index=True)

    network = structural.get("interaction_network")
    if isinstance(network, dict):
        nodes, edges = _rows(network.get("nodes")), _rows(network.get("edges"))
        if nodes or edges:
            st.markdown("**Recorded interaction network**")
            if nodes:
                st.caption("Network nodes")
                st.dataframe(nodes, use_container_width=True, hide_index=True)
            if edges:
                st.caption("Network edges")
                st.dataframe(edges, use_container_width=True, hide_index=True)

    interface = structural.get("interface_rin")
    if isinstance(interface, dict) and (interface.get("applicable") or interface.get("interfaces") or interface.get("residue_pairs")):
        st.markdown("**Recorded interface residue interaction network**")
        st.json({key: value for key, value in interface.items() if key not in {"interfaces", "residue_pairs", "nodes"}})
        for key, title in (("interfaces", "Interfaces"), ("residue_pairs", "Residue pairs"), ("nodes", "Nodes")):
            rows = _rows(interface.get(key))
            if rows:
                st.caption(title)
                st.dataframe(rows, use_container_width=True, hide_index=True)

    warnings = structural.get("warnings")
    if warnings:
        for warning in warnings if isinstance(warnings, list) else [warnings]:
            st.warning(str(warning))


def _render_aggregate_analysis(st, run_root: Path, study: dict[str, Any], *, json_loader=None, text_loader=None, key_prefix="md-results") -> None:
    analyses = [item for item in _rows(study.get("analysis_runs")) if item.get("result_path")]
    if not analyses:
        st.info("No aggregate analysis runs are listed for this study.")
        return
    analyses.sort(key=lambda item: (
        item.get("status") in {"completed", "success"},
        str(item.get("source_completed_at") or item.get("source_created_at") or ""),
    ), reverse=True)
    labels = [
        f"{item.get('run_id') or 'analysis'} · {item.get('status') or 'status unknown'}"
        for item in analyses
    ]
    selected_label = st.selectbox("Analysis result", labels, key=f"{key_prefix}-analysis-{study.get('study_id', 'study')}")
    analysis_record = analyses[labels.index(selected_label)]
    st.caption(f"Analysis run status: {analysis_record.get('status') or 'unknown'}")
    report = _safe_optional_json(st, run_root, analysis_record.get("result_path"), "Analysis", json_loader=json_loader)
    if report is None:
        return
    _render_analysis_report(st, report, run_root=run_root, text_loader=text_loader, key_prefix=key_prefix)


def _render_analysis_report(st, report: dict[str, Any], *, run_root: Path | None = None, text_loader=None, key_prefix="md-results") -> None:
    if report.get("schema_version") == 1 and "source" in report and "units" in report:
        _render_md_analysis_summary(st, report, run_root=run_root, text_loader=text_loader, key_prefix=key_prefix)
        return
    replicas = _rows(report.get("replicas"))
    if replicas:
        st.markdown("**Per-replica summary**")
        st.dataframe(replicas, use_container_width=True, hide_index=True)
        failed = [row for row in replicas if row.get("status") not in {"completed", "success"}]
        if failed:
            st.warning(f"{len(failed)} production replica(s) did not report a completed status.")
    series = _rows(report.get("replica_series"))
    if series:
        st.markdown("**RMSD time series**")
        _line_chart(st, series)
    consensus = _rows(report.get("contact_consensus"))
    if consensus:
        st.markdown("**Contact occupancy consensus**")
        st.dataframe(consensus, use_container_width=True, hide_index=True)
    matrix = report.get("contact_matrix")
    if isinstance(matrix, dict):
        residues = matrix.get("residues")
        replicas_axis = matrix.get("replicas")
        occupancy = matrix.get("contact_occupancy")
        if isinstance(residues, list) and isinstance(replicas_axis, list) and isinstance(occupancy, list):
            rows = []
            for residue, values in zip(residues, occupancy):
                if not isinstance(values, list):
                    continue
                row = {"Residue": str(residue)}
                for index, replica in enumerate(replicas_axis):
                    value = _finite(values[index]) if index < len(values) else None
                    row[f"Replica {replica}"] = value
                rows.append(row)
            if rows:
                st.markdown("**Contact occupancy by replica**")
                st.dataframe(rows, use_container_width=True, hide_index=True)
    aggregate = report.get("aggregate")
    if isinstance(aggregate, dict) and aggregate:
        st.markdown("**Available scalar summaries**")
        st.json(aggregate)
    if not any((replicas, series, consensus, isinstance(matrix, dict), isinstance(aggregate, dict) and aggregate)):
        st.info("The analysis result contains no supported RMSD, contact, or scalar summary fields.")


def _artifact_text(run_root: Path | None, relative_path: str, *, text_loader=None, description="artifact") -> str | None:
    parts = _safe_relative_parts(relative_path)
    safe_relative = PurePosixPath(*parts).as_posix()
    if text_loader is not None:
        return text_loader(safe_relative, description)
    if run_root is None:
        return None
    path = _safe_artifact_path(run_root, safe_relative)
    return path.read_text(encoding="utf-8")


def _csv_rows(text: str | None, *, limit: int = 20000) -> list[dict[str, str]]:
    if not text:
        return []
    try:
        return [row for _, row in zip(range(limit), csv.DictReader(io.StringIO(text)))]
    except (csv.Error, UnicodeError):
        return []


def _render_md_analysis_summary(st, summary: dict[str, Any], *, run_root: Path | None = None, text_loader=None, key_prefix="md-results") -> None:
    """Render the standalone MDTraj summary and its small CSV outputs."""
    if summary.get("success") is not True:
        st.warning(f"MD trajectory analysis did not report success: {summary.get('error') or 'status unavailable'}")
        return
    source = summary.get("source") if isinstance(summary.get("source"), dict) else {}
    units = summary.get("units") if isinstance(summary.get("units"), dict) else {}
    parameters = summary.get("parameters") if isinstance(summary.get("parameters"), dict) else {}
    st.markdown("**Trajectory analysis summary**")
    st.json({
        "Analysis": summary.get("analysis", "not recorded"),
        "Selected frames": summary.get("selected_frame_count", "not recorded"),
        "Selected time range (ns)": summary.get("selected_time_range_ns", "not recorded"),
        "Trajectory atoms": summary.get("trajectory_atom_count", "not recorded"),
        "Units": units,
        "Sampling and selections": {key: parameters.get(key) for key in (
            "stride", "chunk_size", "start_ns", "end_ns", "contact_cutoff_angstrom", "pbc_mode_requested",
        ) if key in parameters} | (summary.get("selection") if isinstance(summary.get("selection"), dict) else {}),
        "Source topology": source.get("topology", {}),
        "Source trajectory": source.get("trajectory", {}),
        "Time origin and source": summary.get("source_time", {}),
        "Source segments": summary.get("source_segments", []),
    })
    warnings = summary.get("warnings")
    if warnings:
        for warning in warnings if isinstance(warnings, list) else [warnings]:
            st.warning(str(warning))
    time_path = "artifacts/md_analysis/timeseries.csv"
    timeseries = _csv_rows(_artifact_text(run_root, time_path, text_loader=text_loader, description="MD analysis timeseries CSV"))
    if timeseries:
        points = []
        for row in timeseries:
            time = _finite(row.get("time_ns"))
            backbone = _finite(row.get("protein_aligned_backbone_rmsd_angstrom", row.get("backbone_rmsd_angstrom")))
            ligand = _finite(row.get("ligand_rmsd_angstrom"))
            if time is not None and (backbone is not None or ligand is not None):
                points.append({"time_ns": time, "backbone_rmsd_angstrom": backbone, "ligand_rmsd_angstrom": ligand})
        if points:
            st.markdown("**RMSD over time (ns, Å)**")
            st.line_chart(points, x="time_ns", y=["backbone_rmsd_angstrom", "ligand_rmsd_angstrom"])
    contacts_path = "artifacts/md_analysis/residue_contacts.csv"
    contacts = _csv_rows(_artifact_text(run_root, contacts_path, text_loader=text_loader, description="MD residue contact CSV"))
    if contacts:
        st.markdown("**Geometric ligand contact occupancy**")
        st.caption(str(summary.get("contact_map_semantics") or "Geometric contact occupancy; not hydrogen-bond analysis."))
        st.dataframe(contacts, use_container_width=True, hide_index=True)
    if not timeseries:
        st.info("MD RMSD timeseries CSV is unavailable.")
    if not contacts:
        st.info("MD residue contact occupancy CSV is unavailable.")


def _mmgbsa_payload(result: dict[str, Any]) -> dict[str, Any] | None:
    payload = result.get("mmgbsa")
    return payload if isinstance(payload, dict) else None


def _endpoint_rows(endpoint: dict[str, Any], manifest_row: dict[str, Any]) -> list[dict[str, Any]]:
    mmgbsa = _mmgbsa_payload(endpoint) or endpoint
    gb = mmgbsa.get("gb") if isinstance(mmgbsa.get("gb"), dict) else {}
    pb = mmgbsa.get("pb") if isinstance(mmgbsa.get("pb"), dict) else {}
    gb_delta = gb.get("delta") if isinstance(gb.get("delta"), dict) else mmgbsa.get("delta")
    pb_delta = pb.get("delta") if isinstance(pb.get("delta"), dict) else {}
    fields = (
        ("GB ΔG bind", gb_delta, "delta_g_bind_total_kcal_mol"),
        ("GB ΔMM", gb_delta, "delta_mm_kcal_mol"),
        ("GB ΔGBSA", gb_delta, "delta_gbsa_kcal_mol"),
        ("GB Δ nonpolar", gb_delta, "delta_nonpolar_kcal_mol"),
        ("PB ΔG bind", pb_delta, "delta_g_bind_total_kcal_mol"),
        ("PB ΔMM", pb_delta, "delta_mm_kcal_mol"),
        ("PB ΔPBSA", pb_delta, "delta_gbsa_kcal_mol"),
        ("PB Δ nonpolar", pb_delta, "delta_nonpolar_kcal_mol"),
    )
    row = {
        "Run ID": str(manifest_row.get("run_id") or "—"),
        "Manifest status": str(manifest_row.get("status") or "unknown"),
        "Result status": str(mmgbsa.get("status") or endpoint.get("success") or "unknown"),
        "Superseded by": str(manifest_row.get("superseded_by_run_id") or "—"),
    }
    for label, payload, key in fields:
        row[f"{label} (kcal/mol)"] = _finite(payload.get(key)) if isinstance(payload, dict) else None
    return [row]


def _render_endpoint_details(st, endpoint: dict[str, Any], record: dict[str, Any], production: dict[str, Any] | None = None) -> None:
    mmgbsa = _mmgbsa_payload(endpoint) or endpoint
    metadata = mmgbsa.get("metadata") if isinstance(mmgbsa.get("metadata"), dict) else {}
    params = record.get("parameters") if isinstance(record.get("parameters"), dict) else {}
    details = {
        "Method": mmgbsa.get("method") or record.get("method") or metadata.get("method") or params.get("backend") or "not recorded",
        "Backend": mmgbsa.get("backend") or record.get("backend") or metadata.get("backend") or params.get("backend") or "not recorded",
        "Units": mmgbsa.get("units") or record.get("units") or "kcal/mol (reported terms)",
        "Frames analyzed": metadata.get("n_frames_analyzed", "not recorded"),
        "Frames total": metadata.get("n_frames_total", "not recorded"),
        "Sampling stride": params.get("stride", metadata.get("stride", "not recorded")),
        "Sampling window start (%)": params.get("start_pct", "not recorded"),
        "Sampling window end (%)": params.get("end_pct", "not recorded"),
        "Source production run": record.get("source_production_run_id") or "not recorded",
        "Recorded trajectory SHA-256": record.get("recorded_trajectory_sha256") or "not recorded",
        "Copied trajectory SHA-256": record.get("copied_trajectory_sha256") or "not recorded",
        "Trajectory snapshot matches": record.get("trajectory_snapshot_matches", "not recorded"),
        "Source MD job": record.get("source_md_job") or "not recorded",
        "Legacy source MD job": record.get("legacy_source_md_job") or "not recorded",
        "Source trajectory reference": record.get("source_trajectory_ref") or "not recorded",
    }
    if production:
        details.update({
            "Source production status": production.get("status") or "not recorded",
            "Source segment duration (ns)": _finite(production.get("segment_duration_ns")),
            "Source previous duration (ns)": _finite(production.get("previous_duration_ns")),
            "Source requested total (ns)": _finite(production.get("target_duration_ns")),
        })
    st.json(details)
    if record.get("trajectory_snapshot_matches") is False:
        st.warning("The endpoint input records a different trajectory snapshot from the copied production trajectory.")
    lineage_warnings = record.get("trajectory_lineage_warnings")
    if lineage_warnings:
        for warning in lineage_warnings if isinstance(lineage_warnings, list) else [lineage_warnings]:
            st.warning(str(warning))
    warnings = mmgbsa.get("warnings")
    if warnings:
        st.warning(str(warnings))


def _render_study(st, run_root: Path, study: dict[str, Any], *, json_loader=None, text_loader=None, key_prefix="md-results", remote=False) -> None:
    st.markdown("**Study identity and selection**")
    st.dataframe(_field_rows(study), use_container_width=True, hide_index=True)

    production = _rows(study.get("production_runs"))
    if production:
        st.markdown("**Production lineage**")
        lineage = []
        depended_on = {str(row.get("continuation_of_run_id") or "") for row in production}
        for item in production:
            status = str(item.get("status") or "unknown")
            lineage.append({
                "Run ID": item.get("run_id"),
                "Status": status,
                "Replica": item.get("repeat_index"),
                "This segment (ns)": _finite(item.get("segment_duration_ns")),
                "Previous (ns)": _finite(item.get("previous_duration_ns")),
                "Requested total (ns)": _finite(item.get("target_duration_ns")),
                "Trajectory appended": item.get("trajectory_appended", "not recorded"),
                "Continues from": item.get("continuation_of_run_id") or "—",
            })
        st.dataframe(lineage, use_container_width=True, hide_index=True)
        failed = [row for row in production if row.get("status") not in {"completed", "success"}]
        if failed:
            dependency_only = [row for row in failed if str(row.get("run_id") or "") in depended_on]
            if dependency_only:
                st.warning(
                    f"{len(dependency_only)} failed production run(s) are retained as continuation prerequisites; "
                    "their files support lineage only and are not presented as successful production metrics."
                )
            other_failed = len(failed) - len(dependency_only)
            if other_failed:
                st.warning(f"{other_failed} production run(s) have failed or non-completed status.")
    else:
        st.warning("No production runs are recorded for this study.")

    endpoints = [item for item in _rows(study.get("mmgbsa_runs")) if item.get("result_path")]
    if endpoints:
        st.markdown("**Endpoint energy runs**")
        endpoint_tables = []
        loaded_endpoints = []
        selected_endpoints = endpoints
        if remote and len(endpoints) > 1:
            labels = [f"{item.get('run_id') or 'MM/GBSA'} · {item.get('status') or 'status unknown'}" for item in endpoints]
            selected = st.selectbox("Endpoint result", labels, key=f"{key_prefix}-endpoint-{study.get('study_id', 'study')}")
            selected_endpoints = [endpoints[labels.index(selected)]]
        for item in selected_endpoints:
            endpoint_result = _safe_optional_json(st, run_root, item.get("result_path"), f"MM/GBSA {item.get('run_id') or ''}", json_loader=json_loader)
            if endpoint_result is None:
                continue
            loaded_endpoints.append((item, endpoint_result))
            endpoint_tables.extend(_endpoint_rows(endpoint_result, item))
        if endpoint_tables:
            st.dataframe(endpoint_tables, use_container_width=True, hide_index=True)
        for item, endpoint_result in loaded_endpoints:
            label = item.get("run_id") or "MM/GBSA"
            st.caption(f"{label} · endpoint method, units and sampling")
            source_id = str(item.get("source_production_run_id") or "")
            source_production = next((row for row in production if str(row.get("run_id") or "") == source_id), None)
            _render_endpoint_details(st, endpoint_result, item, source_production)
        superseded = [item for item in endpoints if item.get("superseded_by_run_id")]
        if superseded:
            st.info(f"{len(superseded)} endpoint result(s) are marked superseded; their original values remain visible above.")
    else:
        st.info("No MM/GBSA endpoint runs are listed for this study.")

    _render_aggregate_analysis(st, run_root, study, json_loader=json_loader, text_loader=text_loader, key_prefix=key_prefix)

    # Show one stored production report per selection. This remains available
    # when an aggregate report exists because aggregates omit rich per-run maps.
    production_with_results = [item for item in production if item.get("result_path")]
    history = study.get("source_history")
    if isinstance(history, dict) and history:
        st.markdown("**Source analysis history**")
        st.caption("History metadata references other source jobs; these paths are shown but are not opened here.")
        st.json(history)
    if production_with_results:
        production_with_results.sort(key=lambda item: (
            item.get("status") in {"completed", "success"},
            str(item.get("source_completed_at") or item.get("source_created_at") or ""),
        ), reverse=True)
        labels = [
            f"{item.get('run_id') or 'production'} · {item.get('status') or 'status unknown'}"
            + (f" · {item.get('source_completed_at')}" if item.get("source_completed_at") else "")
            for item in production_with_results
        ]
        st.markdown("**Production analysis details**")
        selected = st.selectbox("Production result", labels, key=f"{key_prefix}-production-{study.get('study_id', 'study')}")
        item = production_with_results[labels.index(selected)]
        st.caption(f"Production {item.get('run_id') or ''} · status {item.get('status') or 'unknown'}")
        if item.get("status") not in {"completed", "success"}:
            st.warning("This production run is not marked completed; its stored analytics are not shown as successful results.")
        else:
            production_result = _safe_optional_json(
                st, run_root, item.get("result_path"),
                f"Production {item.get('run_id') or ''}", json_loader=json_loader,
            )
            if production_result is not None:
                md_result = production_result.get("md_result") if isinstance(production_result.get("md_result"), dict) else {}
                if production_result.get("success") is not True or md_result.get("status") in {"failed", "cancelled", "error"}:
                    st.warning("The source result does not report successful production; stored values are not rendered as successful analysis.")
                else:
                    analytics = md_result.get("analytics") if isinstance(md_result.get("analytics"), dict) else {}
                    rmsd = analytics.get("rmsd") if isinstance(analytics.get("rmsd"), dict) else {}
                    structural = analytics.get("structural_dynamics") if isinstance(analytics.get("structural_dynamics"), dict) else {}
                    if structural or rmsd:
                        _render_structural_dynamics(
                            st, structural, rmsd,
                            label=str(item.get("repeat_index") or item.get("run_id") or "production"),
                        )
                    if analytics.get("error"):
                        st.warning(f"Production analytics reported: {analytics['error']}")
                    if not structural and not rmsd and not analytics.get("error"):
                        st.info("This production result has no recorded structural-dynamics or RMSD fields.")


def _render_legacy_result(st, run_root: Path, result: dict[str, Any]) -> None:
    """Render direct result JSONs used by non-aggregate MD jobs."""
    candidates = [result]
    try:
        md_result_path = run_root / "artifacts" / "md" / "result.json"
        if md_result_path.resolve(strict=True).is_relative_to(run_root):
            md_payload = json.loads(md_result_path.read_text(encoding="utf-8"))
            if isinstance(md_payload, dict):
                candidates.insert(0, md_payload)
    except (OSError, ValueError, TypeError):
        pass
    rendered = False
    for candidate in candidates:
        md_result = candidate.get("md_result") if isinstance(candidate.get("md_result"), dict) else {}
        analytics = md_result.get("analytics") if isinstance(md_result.get("analytics"), dict) else {}
        if analytics:
            st.markdown("**Production analytics**")
            rmsd = analytics.get("rmsd") if isinstance(analytics.get("rmsd"), dict) else {}
            structural = analytics.get("structural_dynamics") if isinstance(analytics.get("structural_dynamics"), dict) else {}
            times_ps = structural.get("time_ps") or rmsd.get("time_ps") or []
            series = [{
                "replica": candidate.get("run_id") or "production",
                "time_ns": [(_finite(value) / 1000.0) if _finite(value) is not None else None for value in times_ps] if isinstance(times_ps, list) else [],
                "backbone_rmsd_angstrom": rmsd.get("backbone_rmsd_angstrom"),
                "ligand_rmsd_angstrom": rmsd.get("ligand_rmsd_angstrom"),
            }]
            _line_chart(st, series, label="Production")
            if analytics.get("error"):
                st.warning(f"Production analytics reported: {analytics['error']}")
            rendered = True
        mmgbsa = _mmgbsa_payload(candidate) or _mmgbsa_payload(md_result)
        if mmgbsa:
            st.markdown("**MM/GBSA endpoint**")
            st.caption(f"Status: {mmgbsa.get('status') or 'unknown'}")
            _render_endpoint_details(st, candidate, {})
            st.dataframe(_endpoint_rows(candidate, {"run_id": candidate.get("run_id"), "status": mmgbsa.get("status")}), use_container_width=True, hide_index=True)
            rendered = True
    if not rendered:
        st.info("No MD analysis or endpoint values are available in this result.")


def render_md_results(
    st,
    run_dir: str | Path | None = None,
    result: dict[str, Any] | None = None,
    *,
    manifest: dict[str, Any] | None = None,
    json_loader: Callable[[str, str], dict[str, Any] | None] | None = None,
    text_loader: Callable[[str, str], str | None] | None = None,
    key_prefix: str | None = None,
    remote: bool = False,
) -> None:
    """Render one selected MD study or the direct legacy MD result payload."""
    try:
        root = _run_root(run_dir) if run_dir is not None else Path(".").resolve()
        manifest = load_md_study_manifest(run_dir, preloaded=manifest, json_loader=json_loader)
    except (OSError, ValueError) as exc:
        st.warning(f"MD analysis metadata is unavailable: {exc}")
        if run_dir is not None and json_loader is None:
            _render_legacy_result(st, Path(run_dir).expanduser().resolve(strict=False), result if isinstance(result, dict) else {})
        return

    if manifest is None:
        if run_dir is not None and json_loader is None:
            summary_path = root / "artifacts" / "md_analysis" / "summary.json"
            try:
                _safe_artifact_path(root, "artifacts/md_analysis/summary.json")
                summary = _read_json_artifact(root, "artifacts/md_analysis/summary.json")
            except (OSError, ValueError):
                summary = None
            if summary is not None:
                _render_md_analysis_summary(st, summary, run_root=root, text_loader=text_loader, key_prefix=key_prefix or f"md-results-{root.name}")
                return
        _render_legacy_result(st, root, result if isinstance(result, dict) else {}) if run_dir is not None and json_loader is None else st.info("No MD study manifest is available.")
        return
    studies = _rows(manifest.get("studies"))
    if not studies:
        st.info("The MD study manifest contains no studies.")
        return
    labels = [_study_label(study, index) for index, study in enumerate(studies)]
    prefix = key_prefix or f"md-results-{root.name}"
    selected = st.selectbox("MD study", labels, key=f"{prefix}-study")
    _render_study(st, root, studies[labels.index(selected)], json_loader=json_loader, text_loader=text_loader, key_prefix=prefix, remote=remote)
