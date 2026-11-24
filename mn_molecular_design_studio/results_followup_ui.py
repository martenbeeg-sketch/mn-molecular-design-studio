"""Review follow-up analyses for exact selections made in My Results."""
import hashlib
import json
from pathlib import Path
from uuid import uuid4

import streamlit as st


def open_followup(kind, job, selection, *, remote=False):
    fields = ("task_group", "run_id", "project_id", "status", "node_id", "cluster_job_id", "node_name")
    source_job = {field: job.get(field) for field in fields}
    if not remote:
        source_job.pop("cluster_job_id", None)
    st.session_state["results-followup"] = {"kind": kind, "job": source_job, "selection": selection}
    st.session_state["page"] = "Results follow-up"
    st.rerun()


def selected_docking_option(options, selection):
    """A preview without a current score-row identity cannot launch analysis."""
    matches = [option for option in options if option["row_index"] == selection.get("row_index")
               and option["compound_id"] == selection.get("compound_id")
               and option["pose_relative_path"] == selection.get("pose_path")]
    if len(matches) != 1:
        raise ValueError("The selected pose could not be matched to its exact saved docking score row.")
    option = matches[0]
    if not option["available"]:
        raise ValueError(option["unavailable_reason"] or "The selected pose inputs are unavailable.")
    return option


def _render_docking(job, selection, user, workflow_schema, submit):
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_pose_options
    from mn_molecular_design_studio.quality_tools import WORKFLOWS

    st.title("Analyze this pose")
    project = str(job.get("project_id") or "")
    if not project or job.get("status") not in {"completed", "succeeded"}:
        st.info("Choose a completed docking job assigned to a project.")
        return
    try:
        if jid := job.get("cluster_job_id"):
            from mn_molecular_design_studio.cluster.ui import _fetch_remote_artifact_text
            source_job = "cluster:" + jid
            source_ref = selection["source_identity"]
            node = str(job["node_id"])
            token = str(st.session_state.get("studio-session-token") or "development-mode")
            scores = _fetch_remote_artifact_text(jid, "artifacts/docking_scores.csv", token,
                                                f"cluster-docking-scores-{jid}", label="docking score table")
            compounds = _fetch_remote_artifact_text(jid, "artifacts/input/compound_set.json", token,
                                                   f"followup-compounds-{jid}", label="consumed compound set")
            if scores is None or compounds is None:
                st.button("Refresh source files", key="followup-dock-refresh")
                return
            options = docking_pose_options(source_job_reference=source_ref, source_node_id=node,
                scores_text=scores, compound_set_payload=compounds, artifact_paths=selection.get("artifacts", []))
        else:
            from mn_molecular_design_studio.auth import require_job_access
            from mn_molecular_design_studio.job_catalog import get_job_access_record
            from mn_molecular_design_studio.runtime import runs_root
            group, run = job["task_group"], job["run_id"]
            record = get_job_access_record(group, run)
            require_job_access(user, record)
            if not record or str(record.get("project_id") or "") != project:
                raise ValueError("The source job's project has changed. Reopen it from My Results.")
            source_job = f"{group}/{run}"
            root = Path(runs_root()).resolve(strict=True)
            run_dir = (root / record["run_relative_path"]).resolve(strict=True)
            if not run_dir.is_relative_to(root):
                raise ValueError("The source Results directory is outside managed storage.")
            options = docking_pose_options(run_dir)
        pose = selected_docking_option(options, selection)
    except (OSError, PermissionError, ValueError) as exc:
        st.warning(f"Pose analysis unavailable. {exc}")
        return
    st.markdown(f"**{pose['label']}**")
    st.caption(f"Pose file: {pose['pose_relative_path']} · receptor: {Path(pose['receptor_path']).name}")
    st.caption(f"Source run: {job.get('run_id')} · {job.get('node_name') or 'This workstation'}")
    methods = {workflow: WORKFLOWS[workflow]["name"]
               for workflow in ("analysis.native_ligand_interactions", "analysis.plip", "analysis.pandamap", "analysis.posebusters")}
    workflow = st.selectbox("Analysis method", list(methods), format_func=methods.get, key="followup-docking-method")
    if st.button("Queue pose analysis", type="primary", key="followup-docking-submit"):
        try:
            settings = {"_execution": {"placement": "auto"}} if config.enabled() else {}
            submit(workflow_schema(workflow), workflow, {
                "source_job": source_job, "candidate_ids": [pose["pose_id"]],
            }, project, settings, methods[workflow])
        except (OSError, ValueError, RuntimeError) as exc:
            st.error(f"Could not queue pose analysis: {exc}")


def render_followup(user, *, workflow_schema, project_selector, submit):
    source = st.session_state.get("results-followup")
    if st.button("← Back to Results", key="followup-back"):
        st.session_state.pop("results-followup", None)
        st.session_state["page"] = "Job details"
        st.rerun()
    if not source:
        st.info("Choose a result to analyze from My Results.")
        return

    def submit_followup(schema, workflow, parameters, project, settings, label):
        identity = json.dumps([user["id"], project, workflow, parameters, settings], sort_keys=True, default=str)
        attempt = hashlib.sha256(identity.encode()).hexdigest()
        attempts = st.session_state.setdefault("results-followup-attempts", {})
        record = attempts.setdefault(attempt, {"idempotency_key": uuid4().hex})
        if record.get("reference"):
            source["receipt"] = {**record, "attempt": attempt}
            st.info("This analysis was already submitted with these settings.")
            return
        execution = settings.get("_execution") or {"placement": "local"}
        settings = {**settings, "_execution": {**execution, "idempotency_key": record["idempotency_key"]}}
        try:
            result = submit(schema, workflow, parameters, project, settings, label)
        except (OSError, ValueError, RuntimeError):
            st.caption("If the submission reply timed out, check My Results for the job. Retrying here uses the same request key.")
            raise
        job = (result or {}).get("job") or {}
        if job.get("cluster_job_id"):
            reference = "cluster:" + str(job["cluster_job_id"])
        elif job.get("task_group") and job.get("run_id"):
            reference = f"{job['task_group']}/{job['run_id']}"
        else:
            return
        record.update(reference=reference, label=label, job_code=str(job.get("job_code") or job.get("run_id") or ""))
        source["receipt"] = {**record, "attempt": attempt}

    if source["kind"] == "dock":
        _render_docking(source["job"], source["selection"], user, workflow_schema, submit_followup)
    else:
        from mn_molecular_design_studio.md_analysis_ui import render_md_analysis_setup
        st.title("Analyze trajectory")
        render_md_analysis_setup(workflow_schema("ligand.md.analysis"), user=user,
            project_selector=project_selector, submit=submit_followup, source=source)
    if receipt := source.get("receipt"):
        if st.button("Open analysis Results", key="followup-open-analysis-job",
                     help=f"{receipt['label']} · {receipt['job_code']}"):
            st.session_state["selected-job-reference"] = receipt["reference"]
            st.session_state["page"] = "Job details"
            st.rerun()
        if st.button("Prepare another run", key="followup-repeat-analysis"):
            st.session_state["results-followup-attempts"].pop(receipt["attempt"], None)
            source.pop("receipt")
            st.rerun()
