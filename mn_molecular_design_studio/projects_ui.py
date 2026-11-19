"""Project-first research browsing and explicit cleanup review."""

import json
import math
import sqlite3
import time

import streamlit as st

from mn_molecular_design_studio.cluster.ui import display_tool_label
from mn_molecular_design_studio.projects import cleanup_plan, project_inventory


def render_projects(user, projects, *, load_inventory, on_create, on_results, on_database,
                    read_only=False, physical_purge=True, on_archive=None, coordinator=None, on_cancel=None, on_rename=None,
                    on_datasets=None):
    title, create, refresh = st.columns([0.50, 0.28, 0.22], vertical_alignment="center")
    title.title("Projects")
    from mn_molecular_design_studio.entity_properties_ui import render_notice
    render_notice()
    if create.button("New project", icon=":material/add:", key="projects-new", disabled=read_only):
        on_create()
    reload_requested = refresh.button("Refresh", icon=":material/refresh:", key="projects-refresh")
    st.caption("Start with a research project. Browse its runs and saved designs, or review its test data before cleanup.")
    from mn_molecular_design_studio.project_lifecycle_ui import render_pending
    render_pending(user, coordinator, read_only=read_only)
    if name := st.session_state.pop("projects-last-created", None):
        st.success("Project created: " + name)
    if notice := st.session_state.pop("projects-cancel-notice", None):
        st.success(notice)
    show_archived = st.toggle("Show archived projects", key="projects-show-archived")
    projects = [project for project in projects if not project.get("deleted_at")]
    if not projects:
        st.info("No accessible projects. Create a project to group your research.")
        return
    scope = (str(user.get("id") or ""), user.get("role"), tuple(sorted(str(row["id"]) for row in projects)))
    snapshot = st.session_state.get("projects-inventory")
    if reload_requested or not snapshot or snapshot["scope"] != scope or time.monotonic() - snapshot["loaded_at"] > 60:
        try:
            data = load_inventory()
        except (OSError, ValueError, RuntimeError, sqlite3.Error) as exc:
            st.error(f"Could not load project data: {exc}. No cleanup is available until the inventory can be checked.")
            return
        snapshot = {"scope": scope, "loaded_at": time.monotonic(), "data": data}
        st.session_state["projects-inventory"] = snapshot
    data = snapshot["data"]
    if data.get("stale"):
        st.warning("Workstation catalog is offline or incomplete. Counts are a last-known snapshot, not a complete deletion scope.")
    try:
        items = project_inventory(projects, data.get("jobs", []), data.get("collections", []))
    except ValueError as exc:
        st.error(f"Project run relationships need review: {exc}. No cleanup scope has been inferred.")
        return
    search, order_column = st.columns([0.68, 0.32])
    with search:
        query = st.text_input("Find a project", placeholder="Name or research question…", key="projects-search").strip().casefold()
    order = order_column.selectbox("Order", ["Recent activity", "Name", "Most runs"], key="projects-order")
    visible = [item for item in items.values() if (show_archived or not item["project"].get("archived_at")) and (not query or query in (
        str(item["project"].get("name") or "") + " " + str(item["project"].get("description") or "")
        + " " + " ".join(item["project"].get("name_aliases") or [])
    ).casefold())]
    visible.sort(key=lambda item: (str(item["project"]["name"]).casefold(), str(item["project"]["id"])))
    if order == "Recent activity":
        visible.sort(key=lambda item: item["updated_at"], reverse=True)
    elif order == "Most runs":
        visible.sort(key=lambda item: item["runs"], reverse=True)
    if not visible:
        st.info("No matching projects.")
        return
    options = [str(item["project"]["id"]) for item in visible]
    if st.session_state.get("projects-selected") not in options:
        st.session_state["projects-selected"] = options[0]
    selected = st.selectbox("Choose a project", options, key="projects-selected", format_func=lambda value: (
        f"{items[value]['project']['name']} · {items[value]['runs']} runs"
    ))
    st.caption(f"{len(visible)} projects · only metadata loads here; no structure files or workdir scans.")
    item = items[selected]
    project = item["project"]
    with st.container(border=True):
        st.subheader(project["name"])
        from mn_molecular_design_studio.entity_names_ui import render_project_rename
        render_project_rename(user, project, on_rename, read_only=read_only)
        st.write(project.get("description") or "No research description yet.")
        if (project.get("library_admission_review") or {}).get("status") == "pending":
            st.info("Legacy results stay in My results pending admission review. Completion and computed analyses do not mean a design passed. Input targets and reference molecules are separate from admitted results.")
        metrics = st.columns(4)
        for column, label, value in zip(metrics, ("Runs", "Finished jobs", "Active / unclear", "SQL collections"),
                                        (item["runs"], len(item["finished"]), len(item["blocked"]), len(item["collections"]))):
            column.metric(label, value)
        st.caption("A run can contain several tool executions. SQL collections are saved memberships, not a count of unique passing designs.")
        results, database, datasets = st.columns(3)
        if results.button("Open project results", icon=":material/query_stats:", key="projects-results", use_container_width=True):
            on_results(selected)
        if database.button("Open saved designs & molecules", icon=":material/storage:", key="projects-database", use_container_width=True):
            on_database(selected)
        if datasets.button("Datasets", icon=":material/library_books:", key="projects-datasets",
                disabled=on_datasets is None, use_container_width=True):
            on_datasets(selected)
        overview, cleanup = st.tabs(["Runs & data", "Cleanup review"])
        with overview:
            if item["jobs"]:
                tools = sorted({display_tool_label(row) for row in item["jobs"]})
                st.caption("Tools: " + ", ".join(tools))
            if item["logical_runs"]:
                runs = sorted(item["logical_runs"], key=lambda row: row["updated_at"], reverse=True)
                page_count = math.ceil(len(runs) / 25)
                page_key = "projects-runs-page-" + selected
                if st.session_state.get(page_key, 1) > page_count:
                    st.session_state[page_key] = page_count
                run_page = int(st.number_input("Runs page", min_value=1, max_value=page_count, step=1,
                                              key=page_key)) if page_count > 1 else 1
                st.caption(f"{len(runs):,} runs · page {run_page} of {page_count}. Open project results for files and individual tool executions.")
                st.dataframe([{"Run": row["name"], "Status": row["status"], "Executions": row["job_count"],
                               "Workstations": ", ".join(sorted({str(job.get("node_name") or job.get("node_id") or "Not recorded") for job in row["jobs"]})),
                               "Updated": row["updated_at"]} for row in runs[(run_page - 1) * 25:run_page * 25]],
                             hide_index=True, height=330, use_container_width=True)
                with st.expander("Linked executions & cancellation"):
                    render_execution_dependencies(runs, selected, on_cancel=on_cancel, read_only=read_only)
            else:
                st.info("No indexed runs are assigned to this project.")
            with st.expander(f"Saved / imported collections · {len(item['protected'])}"):
                if item["protected"]:
                    st.dataframe([{"Collection": row["name"], "Records": row.get("molecule_count", 0)}
                                  for row in item["protected"]], hide_index=True, use_container_width=True)
                else:
                    st.caption("No saved or imported collections are assigned here. Unassigned legacy data is not inferred to belong to this project.")
            st.caption("Datasets are reusable inputs; Data browser shows admitted designs and their attached analyses. Automatic run collections stay out of this overview.")
        with cleanup:
            st.markdown("#### Review first — nothing is deleted without confirmation")
            st.write(f"{len(item['finished'])} finished executions to review; {len(item['blocked'])} active or unresolved executions are blocked.")
            with st.expander("Technical inventory & existing Results controls"):
                st.caption("This inventory may be cached. Confirmed project actions always obtain a fresh, authoritative review. Counts include hidden executions.")
                st.download_button("Download cleanup review", json.dumps(cleanup_plan(item, catalog_stale=bool(data.get('stale'))), indent=2),
                                   file_name="project-cleanup-review.json", mime="application/json", key="projects-cleanup-download")
                if st.button("Review and select runs in My results", key="projects-cleanup-results"):
                    on_results(selected)
            if read_only:
                st.info("Development review copy: cleanup execution, project deletion and job submission are disabled. Your real data is unchanged.")
            from mn_molecular_design_studio.project_lifecycle_ui import render_actions
            render_actions(user, item, coordinator=coordinator, on_archive=on_archive, read_only=read_only)


def render_execution_dependencies(runs, project_id, *, on_cancel=None, read_only=False):
    from mn_molecular_design_studio.job_dependencies import TERMINAL, identity, parents
    by_id = {run["id"]: run for run in runs}
    key = "projects-inspect-run-" + project_id
    if st.session_state.get(key) not in by_id:
        st.session_state[key] = runs[0]["id"]
    selected = st.selectbox("Inspect a connected run", list(by_id), key=key,
                            format_func=lambda value: by_id[value]["name"] + f" · {by_id[value]['job_count']} jobs")
    run = by_id[selected]
    edges = parents(run["jobs"])
    labels = {identity(row): str(row.get("campaign_name") or row.get("run_id")) for row in run["jobs"]}
    st.caption("Recorded source/parent links connect these executions. Cancelling a job stops its downstream branch, not its upstream source or unrelated branches. Completed results are kept.")
    st.dataframe([{"Job": labels[identity(row)], "Tool": display_tool_label(row), "Status": row.get("status"),
                   "Depends on": ", ".join(labels.get(source, source) for source in sorted(edges[identity(row)])) or "—",
                   "Workstation": row.get("node_name") or row.get("node_id") or "Not recorded"} for row in run["jobs"]],
                 hide_index=True, use_container_width=True)
    active = {identity(row): row for row in run["jobs"] if str(row.get("status") or "").lower() not in TERMINAL}
    if not active:
        st.caption("All recorded executions are finished. Cancellation does not remove saved results or SQL data.")
        return
    target_key = "projects-cancel-execution-" + selected
    if st.session_state.get(target_key) not in active:
        st.session_state[target_key] = next(iter(active))
    target = st.selectbox("Job to cancel", list(active), key=target_key,
                          format_func=lambda value: labels[value] + " · " + str(active[value].get("status")))
    st.caption("Remote cancellations remain requested until the owning workstation confirms. Offline machines cannot be stopped immediately.")
    if st.button("Cancel selected job + downstream", key="projects-cancel-" + selected, disabled=read_only or on_cancel is None):
        try:
            result = on_cancel(active[target])
            st.session_state["projects-cancel-notice"] = f"Cancellation requested, including {result.get('downstream_count', 0)} recorded downstream jobs. Completed results are preserved."
            st.session_state.pop("projects-inventory", None)
            st.rerun()
        except (OSError, ValueError, RuntimeError, sqlite3.Error, LookupError) as exc:
            st.error(str(exc))
