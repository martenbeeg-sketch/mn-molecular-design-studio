"""Explicit project lifecycle actions; no mutation occurs on ordinary rendering."""

import json
import sqlite3

import streamlit as st

from mn_molecular_design_studio.projects import TERMINAL_STATUSES


def render_pending(user, coordinator, *, read_only=False):
    if coordinator is None or read_only or user.get("role") != "admin":
        return
    pending = coordinator.pending(user)
    if not pending:
        return
    with st.expander("Cleanup operations needing attention", expanded=True):
        for operation in pending:
            st.write(operation["project_name"] + " · " + operation["status"])
            if operation.get("error"):
                st.caption(operation["error"])
            if st.button("Continue confirmed cleanup · " + operation["project_name"],
                         key="projects-resume-" + operation["id"], disabled=not operation.get("confirmed_at")):
                try:
                    with st.spinner("Resuming the recorded scope; completed workstation steps are not repeated…"):
                        coordinator.resume(user, operation["id"], confirmation=operation["project_name"], backup_acknowledged=True)
                    st.session_state.pop("projects-inventory", None)
                    st.rerun()
                except (OSError, ValueError, RuntimeError, sqlite3.Error) as exc:
                    st.error(str(exc))


def render_actions(user, item, *, coordinator=None, on_archive=None, read_only=False):
    project = item["project"]
    project_id = str(project["id"])
    busy = (project.get("lifecycle") or {}).get("status") in {"preparing", "running", "attention"}
    archived = bool(project.get("archived_at"))
    permitted = user.get("role") == "admin" or project.get("created_by") == user.get("id")
    st.markdown("#### Project visibility")
    st.caption("Archive hides a project from the default Projects list and new-job selectors. Its Results and Database data remain available; Restore reverses the archive.")
    if st.button("Restore archived project" if archived else "Archive project", key="projects-archive-" + project_id,
                 disabled=read_only or busy or not permitted or on_archive is None):
        try:
            on_archive(project_id, not archived)
            st.rerun()
        except (OSError, ValueError, RuntimeError, sqlite3.Error) as exc:
            st.error(str(exc))
    st.markdown("#### Remove project data")
    st.caption("Finished run folders move out of workdir into a verified recovery store. Selected SQL memberships, scores and representations are removed; shared molecular identities and shared structure files are retained.")
    action = st.radio("Cleanup action", ["Clean selected runs", "Delete whole project"],
                      horizontal=True, key="projects-cleanup-action-" + project_id, disabled=busy)
    mode = "delete" if action == "Delete whole project" else "clean"
    runs = {row["id"]: row for row in item["logical_runs"]}
    selection_key = "projects-cleanup-runs-" + project_id
    if any(value not in runs for value in st.session_state.get(selection_key, [])):
        st.session_state[selection_key] = [value for value in st.session_state[selection_key] if value in runs]
    chosen = st.multiselect("Runs to remove", list(runs),
        format_func=lambda value: runs[value]["name"] + " · " + str(runs[value]["updated_at"] or "date not recorded") + f" · {runs[value]['job_count']} jobs",
        key=selection_key, disabled=mode == "delete" or busy)
    references = sorted({str(row.get("node_id") or "") + ":" + str(row.get("task_group") or "") + "/" + str(row.get("run_id") or "")
                         for value in chosen for row in runs[value]["jobs"]})
    include = st.checkbox("Also remove this project's imported and explicitly saved collections",
                          value=False, key="projects-cleanup-collections-" + project_id, disabled=busy)
    if mode == "delete":
        st.warning("Whole-project deletion includes every execution and downstream analysis assigned to this project. Shared dependencies, active jobs and excluded saved/imported collections block deletion. The project is hidden only after every workstation confirms completion.")
    elif chosen and any(str(row.get("status") or "").casefold() not in TERMINAL_STATUSES for value in chosen for row in runs[value]["jobs"]):
        st.warning("The selection includes active or unresolved jobs. They will block cleanup.")
    if read_only:
        st.info("Research snapshot: destructive actions are disabled. Use the disposable cleanup sandbox to test the complete workflow.")
    elif user.get("role") != "admin":
        st.info("An administrator must prepare and confirm data cleanup.")
    operation_key = "projects-cleanup-operation-" + project_id
    request_key = (mode, tuple(references), include)
    operation = st.session_state.get(operation_key)
    if operation and operation.get("status") == "review" and operation.get("request_key") != request_key:
        st.session_state.pop(operation_key, None)
        operation = None
    disabled = read_only or coordinator is None or user.get("role") != "admin" or busy or (mode == "clean" and not references and not include)
    if st.button("Prepare cleanup review", key="projects-prepare-cleanup-" + project_id, disabled=disabled):
        try:
            with st.spinner("Checking live SQL, run folders, dependencies and every trusted workstation…"):
                operation = coordinator.prepare(user, project_id, mode=mode, selected_runs=references, include_collections=include)
            operation["request_key"] = request_key
            st.session_state[operation_key] = operation
        except (OSError, ValueError, RuntimeError, sqlite3.Error) as exc:
            st.error("Cleanup review could not establish a safe scope: " + str(exc))
    if not operation and coordinator is not None and not read_only and user.get("role") == "admin":
        pending = coordinator.pending(user, project_id)
        if pending:
            operation = pending[-1]
            st.session_state[operation_key] = operation
    if not operation:
        return
    plans = list(operation["plans"].values())
    st.markdown("#### Review before confirmation")
    st.write(f"{sum(len(plan['runs']) for plan in plans)} executions · "
             f"{sum(len(group['members']) for plan in plans for group in plan['groups'])} SQL memberships · "
             f"{sum(len(plan['structures']) for plan in plans)} structures · {len(plans)} workstations")
    st.caption("These are memberships and representations, not counts of unique passing designs. No smoke/test naming heuristic selects or deletes data.")
    for node_id, plan in operation["plans"].items():
        with st.expander("Workstation " + node_id):
            st.write(f"{len(plan['runs'])} run folders; {len(plan['retained'])} saved/imported collections kept.")
            if plan["runs"]:
                st.dataframe([{"Run reference": row["reference"], "Bytes in workdir": sum(value[2] for value in row["tree"] if value[1] == "file")}
                              for row in plan["runs"]], hide_index=True, use_container_width=True)
    st.download_button("Download reviewed deletion plan", json.dumps(operation, indent=2, default=list),
                       file_name="project-cleanup-plan.json", mime="application/json", key="projects-cleanup-plan-" + operation["id"])
    if operation["blocked"]:
        for problem in operation["blocked"]:
            st.error(problem)
        return
    if operation.get("error"):
        st.error("Cleanup needs attention: " + operation["error"])
    if operation["status"] == "completed":
        st.success("Cleanup completed on every workstation. Recovery files remain available; no shared-file garbage collection was performed.")
        return
    if operation["status"] == "aborted":
        st.info("No data was removed. Prepare a fresh review; the previous scope changed.")
        return
    confirmation = st.text_input("Type the exact project name to confirm", key="projects-cleanup-name-" + operation["id"])
    acknowledged = st.checkbox("Create and verify a recovery backup before removal. I understand recovery files still use disk space.",
                               key="projects-cleanup-backup-" + operation["id"])
    label = "Continue interrupted cleanup" if operation["status"] != "review" else "Confirm project deletion" if mode == "delete" else "Confirm data cleanup"
    if st.button(label, type="primary", key="projects-confirm-cleanup-" + operation["id"],
                 disabled=read_only or coordinator is None or confirmation != operation["project_name"] or not acknowledged):
        try:
            with st.spinner("Verifying recovery files and coordinating cleanup. Do not close the app until status returns…"):
                result = coordinator.resume(user, operation["id"], confirmation=confirmation, backup_acknowledged=acknowledged)
            st.session_state[operation_key] = result
            st.session_state.pop("projects-inventory", None)
            st.rerun()
        except (OSError, ValueError, RuntimeError, sqlite3.Error) as exc:
            st.error(str(exc))
