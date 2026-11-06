"""Run-output selection used by native tiles and graph source nodes."""
import streamlit as st


@st.cache_data(ttl=15, max_entries=32, show_spinner="Loading Run outputs…")
def options_for_user(project_id, user, token="", require_complex=True, stages=None, allow_docking=False):
    from mn_molecular_design_studio.logical_runs import input_options
    return input_options(project_id, user, token=token, require_complex=require_complex,
                         stages=stages, allow_docking=allow_docking)


def options(project_id, *, require_complex=True, stages=None, allow_docking=False):
    from mn_molecular_design_studio.auth import development_user, is_development_mode
    from mn_molecular_design_studio.logical_runs import input_options
    user = st.session_state.get("_molecular_database_user")
    if not user and is_development_mode():
        user = development_user()
    if not user:
        raise ValueError("The signed-in user context is unavailable.")
    return options_for_user(project_id, user, token=str(st.session_state.get("studio-session-token") or ""),
                            require_complex=require_complex, stages=stages, allow_docking=allow_docking)


def _docking_pose_options(row):
    from mn_molecular_design_studio.logical_runs import reference as source_reference
    from mn_molecular_design_studio.cluster import config
    source_ref = source_reference(row)
    if row.get("node_id") == config.identity():
        from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_pose_options
        return docking_pose_options(row.get("run_dir"), source_job_reference=source_ref,
                                    source_node_id=str(row.get("node_id") or ""))
    # Reuse the authorized Results artifact command already used by Scientific
    # Inspector; fetch only the small score and compound manifests for the UI.
    # Indexed Results rows do not carry a complete artifact listing. The score
    # rows supply candidate pose paths, which are validated as safe relative
    # paths here and resolved inside the authorized source job during execution.
    artifacts = [item for item in row.get("artifacts", []) if isinstance(item, dict)]
    from mn_molecular_design_studio.cluster.ui import _fetch_remote_artifact_text
    from mn_molecular_design_studio.auth import development_user, is_development_mode
    user = st.session_state.get("_molecular_database_user")
    if not user and is_development_mode():
        user = development_user()
    if not user:
        return []
    jid = str(row.get("cluster_job_id") or "")
    if not jid:
        local_ref = str(row.get("local_reference") or f"{row.get('task_group', '')}/{row.get('run_id', '')}")
        jid = config.job_id(str(row.get("node_id") or ""), local_ref)
    token = str(st.session_state.get("studio-session-token") or "")
    scores = _fetch_remote_artifact_text(jid, "artifacts/docking_scores.csv", token,
        f"pipeline-docking-scores-{jid}", label="docking score table")
    compound_set = _fetch_remote_artifact_text(jid, "artifacts/input/compound_set.json", token,
        f"pipeline-docking-compounds-{jid}", label="docking compound set")
    if scores is None or compound_set is None:
        return []
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_pose_options
    return docking_pose_options(source_job_reference=source_ref, source_node_id=str(row.get("node_id") or ""),
        scores_text=scores, compound_set_payload=compound_set, artifact_paths=artifacts,
        artifact_manifest_complete=row.get("artifact_manifest_complete"))


def select(project_id, key, *, require_complex=True, stages=None, reference=None, allow_docking=False,
           selected_docking_poses=None):
    if not project_id:
        st.info("Select a project to choose a Run output.")
        return None
    choices = options(project_id, require_complex=require_complex, stages=stages, allow_docking=allow_docking)
    if not choices:
        st.info("No compatible completed Run outputs are available in this project.")
        return None
    by_reference = {choice["reference"]: choice for choice in choices}
    labels = {reference: str(choice["label"]) for reference, choice in by_reference.items()}
    references = list(by_reference)
    selected_ref = st.selectbox("Run results", references, format_func=lambda value: labels[value],
        index=references.index(reference) if reference in references else 0, key=key + "-" + project_id)
    selected = by_reference[selected_ref]
    st.caption("Uses the saved Run output across its source workstations. Select individual docking poses below when routing raw ligand results to analysis.")
    if not allow_docking:
        return selected
    source_rows = selected["row"].get("_run_sources", [])
    docking_rows = [row for row in source_rows if row.get("has_docking_results")]
    if not docking_rows:
        return selected
    pose_choices = {}
    unavailable_scores = []
    source_errors = []
    for row in docking_rows:
        source_ref = __import__("mn_molecular_design_studio.logical_runs", fromlist=["reference"]).reference(row)
        try:
            poses = _docking_pose_options(row)
        except (OSError, ValueError) as exc:
            source_errors.append(f"{row.get('node_name') or row.get('run_id')}: {exc}")
            poses = []
        for pose in poses:
            if not pose.get("available"):
                unavailable_scores.append({
                    "Source": str(row.get("node_name") or row.get("run_id") or "Docking run"),
                    "Compound": str(pose.get("compound_name") or pose.get("compound_id") or ""),
                    "Pose ID": str(pose.get("pose_id") or "not assigned"),
                    "Score": str(pose.get("label") or ""),
                    "Recorded pose path": str(pose.get("pose_relative_path") or "not recorded"),
                    "Status": str(pose.get("unavailable_reason") or "Pose artifact unavailable"),
                })
                continue
            token = source_ref + "\n" + str(pose["pose_id"])
            pose_choices[token] = {"source_job_reference": source_ref, "pose_id": str(pose["pose_id"]),
                                   "label": str(row.get("node_name") or "Source workstation") + " · " +
                                           str(row.get("job_code") or row.get("run_id") or "Docking") + " · " + pose["label"]}
    if source_errors:
        for error in source_errors[:10]:
            st.warning(f"Could not verify docking artifacts: {error}")
    if unavailable_scores:
        with st.expander(f"Docking score records without an executable pose ({len(unavailable_scores):,})"):
            st.caption("These source score records remain in Results, but are not offered to analysis because their pose artifact or required chemistry is unavailable.")
            st.dataframe(unavailable_scores[:100], hide_index=True, use_container_width=True)
            if len(unavailable_scores) > 100:
                st.caption("Showing the first 100 unavailable records; the full score table remains in the source Run's Files tab.")
    if not pose_choices:
        st.info("No verified docking poses are selectable. A pose is offered only when its authorized receptor, pose coordinates and recorded ligand chemistry are available.")
        return {**selected, "selected_docking_poses": None}
    labels = {value: entry["label"] for value, entry in pose_choices.items()}
    saved = {(str(item.get("source_job_reference")), str(item.get("pose_id")))
             for item in (selected_docking_poses or []) if isinstance(item, dict)}
    default = [value for value, entry in pose_choices.items()
               if (entry["source_job_reference"], entry["pose_id"]) in saved]
    chosen = st.multiselect("Docking poses to analyze", list(pose_choices), default=default,
        format_func=lambda value: labels[value], key=key + "-docking-poses-" + selected_ref)
    if chosen:
        st.caption("The selected poses will be sent to a compatible protein–ligand analysis node (for example, PLIP).")
        from mn_molecular_design_studio.cluster import config
        from mn_molecular_design_studio.logical_runs import reference as source_reference
        remote_refs = {source_reference(row) for row in docking_rows if row.get("node_id") != config.identity()}
        if any(pose_choices[value]["source_job_reference"] in remote_refs for value in chosen):
            st.caption("Remote pose coordinates are checked on the owning workstation before analysis. If a score row refers to a pose that was not retained, the workflow stops before PLIP and reports that source job.")
    return {**selected, "selected_docking_poses": [
        {"source_job_reference": pose_choices[value]["source_job_reference"], "pose_id": pose_choices[value]["pose_id"]}
        for value in chosen] or None}
