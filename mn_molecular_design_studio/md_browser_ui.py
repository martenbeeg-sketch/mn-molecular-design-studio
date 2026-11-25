"""MD study exploration without transferring trajectories into the SQL browser."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.result_details_ui import render_result_details
from mn_molecular_design_studio.md_browser import (
    MDArtifactReader, ligand_atom_block, matching_study, refresh_requests, report_choices, report_label, snapshot_references, study_observations,
)
from mn_molecular_design_studio.md_results_ui import (
    _endpoint_rows, _finite, _render_analysis_report, _render_endpoint_details, _render_structural_dynamics,
)


def _reader(database, user, index, outcome, observation):
    original = next(row for row in outcome["originals"] if row["observation_id"] == observation["id"])
    run = index["runs"][original["run_key"]]

    def fetch(source, relative, description):
        from mn_molecular_design_studio.cluster import config
        from mn_molecular_design_studio.cluster.ui import _fetch_remote_artifact_text

        token = st.session_state.get("studio-session-token") or ""
        if not token:
            return None
        identifier = config.job_id(source[0], source[1] + "/" + source[2])
        key = hashlib.sha256(json.dumps([user["id"], *source, relative,
            observation["context"].get("source_summary_sha256")]).encode()).hexdigest()
        return _fetch_remote_artifact_text(identifier, relative, token, "md-browser-fetch-" + outcome["id"] + "-" + key, description)

    return MDArtifactReader(database, user, observation, run, fetch_text=fetch)


def _summary_table(observations):
    return [{"Analysis": row["namespace"], "Replica": str(row["context"].get("replica_index") or "—"),
        "Measurement": key.replace("_", " "), "Value": value}
        for row in observations for key, value in row["metrics"].items()]


def _section(payload, field):
    value = payload.get(field)
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError("Stored MD " + field + " must be a JSON object.")
    return value


def _thermodynamics(analytics):
    thermo = _section(analytics, "thermodynamics")
    for field, label in (("temperature_k", "Temperature (K)"), ("density_gcm3", "Density (g/cm³)"),
                         ("potential_energy_kjmol", "Potential energy (kJ/mol)"), ("volume_nm3", "Volume (nm³)")):
        times, values = thermo.get("time_ps"), thermo.get(field)
        if not isinstance(times, list) or not isinstance(values, list):
            continue
        stride = max(1, (min(len(times), len(values)) + 1999) // 2000)
        points = [{"Time (ns)": _finite(time) / 1000, label: _finite(value)}
            for time, value in list(zip(times, values))[::stride] if _finite(time) is not None and _finite(value) is not None]
        if points:
            st.markdown("**" + label + "**")
            st.line_chart(points, x="Time (ns)", y=label)
    st.caption("Long thermodynamic series are sampled for display only; the saved report is unchanged.")


def _contact_plot(structural, prefix):
    rows = _section(structural, "contacts").get("residues") or []
    if not isinstance(rows, list):
        raise ValueError("Stored contact residues must be a list.")
    fields = sorted({field for row in rows if isinstance(row, dict) for field, value in row.items()
        if "occupancy" in field and _finite(value) is not None})
    if fields:
        field = st.selectbox("Contact occupancy plot", fields,
            index=fields.index("contact_occupancy") if "contact_occupancy" in fields else 0,
            format_func=lambda value: value.replace("_", " "), key=prefix + "-contact-plot")
        points = [{"Residue": str(row["residue"]), "Occupancy": _finite(row[field])}
            for row in rows if isinstance(row, dict) and row.get("residue") and _finite(row.get(field)) is not None]
        if points:
            st.bar_chart(points, x="Residue", y="Occupancy")


def _render_study_overview(database, user, index, outcome, observations, prefix):
    from mn_molecular_design_studio.md_plot_matrix import validate_report_sources

    choices = report_choices(observations, category="analysis")
    if not choices:
        st.info("No saved multi-replica analysis is attached to this study. Individual production reports remain available under Dynamics and Contacts.")
        return
    lookup = {row["id"]: row for row in choices}
    chosen = st.selectbox("Stored overview report", list(lookup),
        format_func=lambda value: str(lookup[value]["context"].get("source_created_at") or lookup[value]["created_at"])[:16].replace("T", " ") + " · saved multi-replica analysis",
        key=prefix + "-overview-report")
    observation = lookup[chosen]
    reader = _reader(database, user, index, outcome, observation)
    report = reader.report()
    if report is None:
        st.info("The overview report is loading or unavailable. Individual production reports and SQL summaries remain accessible.")
        return
    productions = {str(row["context"]["source_production_run_id"]): str(row["context"].get("replica_index"))
        for row in observations if row["namespace"] == "md.production" and row["context"].get("source_production_run_id")}
    validate_report_sources(report, productions)
    st.caption("All replicas below belong to this one saved analysis report. Earlier sampling windows remain separately selectable; replicas and continuation segments are not silently stitched together.")
    st.caption("The replica / production selectors above control the structure viewer. This overview compares all selected report replicas; Expand overview gives the plots more space.")
    _render_analysis_report(st, report, text_loader=reader.text, key_prefix=prefix + "-overview-" + observation["id"])
    with st.expander("Overview source and production IDs"):
        render_result_details(st, {"Source report": observation["context"], "Replicas": report.get("replicas"), "Analysis window": report.get("analysis_window")})


def render_md_study(database, user, index, entry, selection):
    outcome = entry["outcome"]
    from mn_molecular_design_studio.data_browser import concise_result_label

    st.markdown("### Molecular dynamics")
    st.markdown("**" + concise_result_label(outcome) + "**")
    st.caption("Stored SQL measurements + original analysis artifacts. No trajectory is loaded or analysis rerun.")
    observations = study_observations(database, user, outcome)
    if not observations:
        st.info("No MD observations remain accessible. Refresh the SQL view.")
        return None
    prefix = "md-browser-" + outcome["id"]
    replicas = sorted({str(row["context"]["replica_index"]) for row in observations
        if row["context"].get("replica_index") is not None}, key=lambda value: (len(value), value))
    chosen_replica = selection.get("md_replica")
    replica_column, segment_column = st.columns([.24, .76])
    with replica_column:
        replica = st.selectbox("Replica", replicas, index=replicas.index(chosen_replica) if chosen_replica in replicas else 0,
            key=prefix + "-replica-" + str(chosen_replica)) if replicas else None
    if not replicas:
        st.caption("Replica identity was not recorded for this run.")
    productions = report_choices(observations, replica=replica)
    production = None
    if productions:
        lookup = {row["id"]: row for row in productions}
        with segment_column:
            chosen = st.selectbox("Production segment", list(lookup), format_func=lambda value: report_label(lookup[value]),
                key=prefix + "-production-" + str(replica))
        production = lookup[chosen]
    analysis_column, refresh_column = st.columns([.72, .28], vertical_alignment="bottom")
    with analysis_column:
        categories = ["Dynamics", "Contacts", "Thermodynamics", "Endpoint energy", "SQL summary", "Provenance"]
        if report_choices(observations, category="analysis"):
            categories.insert(0, "Overview · all replicas")
        category = st.selectbox("Analysis", categories,
            key=prefix + "-view")
    with refresh_column:
        refresh = st.button("Refresh", key=prefix + "-refresh", help="Check pending artifact transfers and retry failed transfers.")
    current = [row for row in observations if replica is None or str(row["context"].get("replica_index")) == replica]
    production_id = production["context"].get("source_production_run_id") if production else None
    if refresh:
        from mn_molecular_design_studio.cluster.ui import _request
        from mn_molecular_design_studio.cluster.transport import ClusterUnavailable

        try:
            refresh_requests(st.session_state, "md-browser-fetch-" + outcome["id"] + "-",
                lambda command: _request("/commands/" + command, st.session_state.get("studio-session-token") or ""))
        except (ValueError, ClusterUnavailable) as exc:
            st.caption("The owning workstation could not be checked: " + str(exc))
    structure_context = {"outcome": outcome, "observation": production, "study": None, "reader": None}
    try:
        if production:
            reader = _reader(database, user, index, outcome, production)
            structure_context["reader"] = reader
            try:
                structure_context["study"] = matching_study(reader, outcome)
            except (OSError, ValueError, PermissionError, UnicodeError) as exc:
                st.caption("Structure-link manifest unavailable: " + str(exc))
        if category == "Overview · all replicas":
            _render_study_overview(database, user, index, outcome, observations, prefix)
        elif category == "SQL summary":
            st.dataframe(_summary_table(current), hide_index=True, use_container_width=True)
        elif category == "Provenance":
            st.caption("Segments and replicas remain separate. Requested duration is not proof of completed sampling.")
            for row in current:
                with st.expander(row["namespace"] + " · " + report_label(row)):
                    render_result_details(st, row["context"])
        elif category == "Endpoint energy":
            choices = report_choices(observations, replica=replica, category="endpoint", production_id=production_id)
            if not choices:
                st.info("No endpoint-energy report is attached to this production segment / replica.")
            else:
                lookup = {row["id"]: row for row in choices}
                chosen = st.selectbox("Endpoint calculation", list(lookup), format_func=lambda value: report_label(lookup[value]), key=prefix + "-endpoint-" + str(replica))
                endpoint = lookup[chosen]
                report = _reader(database, user, index, outcome, endpoint).report()
                if report is not None:
                    context = endpoint["context"]
                    st.dataframe(_endpoint_rows(report, context), hide_index=True, use_container_width=True)
                    st.caption("Endpoint calculation linked to the selected production segment. Historical calculations remain separately selectable.")
                    with st.expander("Method, sampling and source provenance"):
                        _render_endpoint_details(st, report, context)
                    if context.get("endpoint_superseded_by_run_id"):
                        st.warning("This historical endpoint calculation is marked superseded.")
                else:
                    st.info("Endpoint report is loading or unavailable. SQL values remain available under SQL summary.")
        elif production:
            report = reader.report()
            if report is None:
                st.info("The report is loading or unavailable on its workstation. Use Refresh MD artifacts; SQL summary remains available.")
            elif report.get("success") is not True or production["context"].get("production_status") in {"failed", "cancelled", "error"}:
                st.warning("The production report does not record success. Its values are not presented as successful analysis.")
            else:
                payload = _section(report, "md_result") or _section(report, "production")
                if payload.get("status") in {"failed", "cancelled", "error"}:
                    st.warning("The production payload records failure; successful analysis is not asserted.")
                else:
                    analytics = _section(payload, "analytics")
                    structural = _section(analytics, "structural_dynamics")
                    rmsd = _section(analytics, "rmsd")
                    local = _section(rmsd, "per_phase_local").get("production")
                    if isinstance(local, dict) and category == "Dynamics":
                        reference = st.selectbox("RMSD timeline", ["Production-local", "Original recorded timeline"], key=prefix + "-rmsd-reference")
                        if reference == "Production-local":
                            rmsd = local
                    st.caption("Each curve uses its own recorded time axis. Continuation segments are not concatenated. The displayed structure is not automatically the RMSD reference.")
                    if category == "Thermodynamics":
                        _thermodynamics(analytics)
                    elif structural or rmsd:
                        contacts = {"contacts", "water_bridges", "salt_bridges", "interaction_network", "interface_rin", "pocket"}
                        selected = {key: value for key, value in structural.items() if key in ({"time_ps", *contacts} if category == "Contacts" else set(structural) - contacts)}
                        if category == "Contacts":
                            from mn_molecular_design_studio.md_atom_interactions import render_interaction_diagram

                            context = production["context"]
                            try:
                                atom_block, atom_note = ligand_atom_block(reader, structure_context.get("study"), report, production_id)
                            except (OSError, ValueError, PermissionError, UnicodeError) as exc:
                                atom_block, atom_note = "", "Ligand atom source could not be verified: " + str(exc)
                            render_interaction_diagram(st, structural,
                                ligand_smiles=context.get("ligand_smiles") or context.get("identity_smiles") or "",
                                ligand_pdb=atom_block, atom_source_note=atom_note,
                                ligand_label=context.get("compound_id") or "Ligand", key_prefix=prefix,
                                scope="Replica " + str(replica or "not recorded") + " · production " + str(production_id or "not recorded"))
                            _contact_plot(structural, prefix)
                        _render_structural_dynamics(st, selected, {} if category == "Contacts" else rmsd,
                            label="replica " + str(replica or "not recorded"), key_prefix=prefix)
                    else:
                        st.info("No rich production curves are recorded in this report. Check the attached trajectory analysis or SQL summary.")
                    with st.expander("Recorded analysis window and reference metadata"):
                        render_result_details(st, {"Analysis window": structural.get("analysis_window"),
                            "Residue numbering": structural.get("residue_numbering"),
                            "Reference": {key: value for key, value in rmsd.items() if not isinstance(value, (list, dict))},
                            "Source": production["context"]})
                    if analytics.get("error"):
                        st.warning(str(analytics["error"]))
        else:
            st.info("No production report is attached at this level. SQL measurements remain available.")
        if category in {"Dynamics", "Contacts"}:
            analyses = report_choices(observations, replica=replica, category="analysis", production_id=production_id)
            if analyses:
                with st.expander("Attached trajectory-analysis reports"):
                    lookup = {row["id"]: row for row in analyses}
                    chosen = st.selectbox("Trajectory analysis", list(lookup), format_func=lambda value: report_label(lookup[value]), key=prefix + "-analysis-" + str(replica))
                    if st.button("Open stored trajectory analysis", key=prefix + "-open-analysis-" + str(replica)):
                        st.session_state[prefix + "-active-analysis"] = chosen
                    if st.session_state.get(prefix + "-active-analysis") == chosen:
                        analysis_reader = _reader(database, user, index, outcome, lookup[chosen])
                        report = analysis_reader.report()
                        if report is not None:
                            _render_analysis_report(st, report, text_loader=analysis_reader.text, key_prefix=prefix + "-aggregate")
    except (OSError, ValueError, PermissionError, UnicodeError) as exc:
        st.warning("The stored MD artifact could not be verified: " + str(exc))
        st.caption("Missing or changed artifacts do not erase the saved SQL measurements.")
    return structure_context


def render_md_structure(database, user, entry, context):
    st.markdown("#### MD structure context")
    snapshots = snapshot_references(context.get("study"), (context.get("observation") or {}).get("context", {}).get("source_production_run_id", ""))
    if snapshots and context.get("reader"):
        lookup = {row["path"]: row for row in snapshots}
        remembered = st.session_state.get("md-browser-snapshot-choice") or {}
        default = remembered.get("chosen") if remembered.get("outcome") == entry["outcome"]["id"] else None
        chosen = st.selectbox("Structure snapshot", list(lookup), index=list(lookup).index(default) if default in lookup else 0,
            format_func=lambda value: lookup[value]["label"], key="md-browser-snapshot-" + entry["outcome"]["id"])
        st.session_state["md-browser-snapshot-choice"] = {"outcome": entry["outcome"]["id"], "chosen": chosen}
        snapshot = lookup[chosen]
        st.caption(snapshot["note"])
        try:
            text = context["reader"].text(chosen, snapshot["label"], expected=snapshot["sha256"])
            if text is not None:
                from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
                from mn_molecular_design_studio.docking_results import pose_visualizations

                molstar_custom_component(pose_visualizations(text, Path(chosen).name, include_pocket_sidechains=False),
                    key="md-browser-snapshot-viewer-" + snapshot["sha256"],
                    height=720 if st.session_state.get("data-browser-layout") == "Structure" else 430,
                    show_controls=False, show_structure_tools=False)
                st.download_button("Download verified snapshot", text, file_name=snapshot["label"].replace(" ", "_") + Path(chosen).suffix)
            else:
                st.info("The selected snapshot is loading or unavailable. No substitute pose is shown.")
        except (OSError, ValueError, PermissionError, UnicodeError) as exc:
            st.warning("Structure snapshot unavailable: " + str(exc))
        return
    source = entry["result"].get("source_structure_id")
    if source:
        from mn_molecular_design_studio.data_browser import load_saved_record
        from mn_molecular_design_studio.data_browser_ui import _render_structure

        record = load_saved_record(database, user, {**entry["result"], "structure_ids": [source]})
        if record and record["structures"]:
            try:
                path = database.structure_file(user, source)
                expected = entry["result"].get("source_structure_digest")
                if expected and hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                    raise ValueError("The SQL source pose no longer matches this MD input.")
                st.caption("Original SQL source pose. Not asserted to be the production-start or RMSD-reference frame.")
                _render_structure(database, user, record, source)
            except (OSError, ValueError, PermissionError) as exc:
                st.warning("Original source pose unavailable: " + str(exc))
            return
    st.info("No verified structure link is available for this MD study. An unrelated pose is not substituted.")
