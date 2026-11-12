"""Molecule/sequence and run/campaign perspectives over saved Database records."""
from __future__ import annotations

import hashlib
import json
import math

import streamlit as st

from mn_molecular_design_studio.database_browser import (
    browse_subjects,
    build_browser_index,
    molecular_subjects,
    run_subjects,
)
from mn_molecular_design_studio.molecular_database import IMPORT_METHODS


def render_browser(database, user, projects, view):
    from mn_molecular_design_studio.database_ui import _projects, render_import_dialog, sync_workspace

    with st.container(key="molecular-database"):
        title, new, refresh, access = st.columns([5, 1.5, 1, 2])
        title.markdown("### " + ("Molecules & sequences" if view == "Molecules" else "Runs & campaigns"))
        if new.button("Import / new", key="molecular-database-new", use_container_width=True):
            st.session_state["molecular-import-open"] = True
            st.session_state["molecular-import-step"] = 0
        if refresh.button("↻", key="molecular-database-refresh", help="Explicitly refresh saved job indexing"):
            st.session_state["molecular-indexing-errors"] = sync_workspace(user, database=database, force=True)
            st.rerun()
        with access.popover("Access & source filters", use_container_width=True):
            collection_label = st.segmented_control("Collection", ["Library", "All indexed"], default="Library", key="molecular-collection-view") or "Library"
            scope_options = ["Mine", "Shared", "All"]
            if st.session_state.get("molecular-scope") not in scope_options:
                st.session_state["molecular-scope"] = "All"
            scope_label = st.segmented_control("Collection access", scope_options, key="molecular-scope") or "All"
            project_id = _projects(projects, "molecular-project-filter", label="Project filter", empty_label="All projects")
            hide_smoke = st.checkbox("Hide confirmed development runs", value=True, key="molecular-hide-smoke-test-groups",
                                     help="Only the previously confirmed owner/run pairs, never names or single-entry groups.")
            methods = [value for value, label in IMPORT_METHODS.items() if st.checkbox(label, value=True, key="molecular-method-" + value)] or ["__none__"]
        scope = {"Mine": "mine", "Shared": "shared", "All": "all"}[scope_label]
        collection = "library" if collection_label == "Library" else "all"
        project_name = next((row["name"] for row in projects if row["id"] == project_id), "All projects")
        st.caption(f"{collection_label} · {scope_label} access · {project_name}. Related history is included only where you have access.")
        search, grouping, ordering = st.columns([4, 2, 2])
        query = search.text_input("Search grouped records", placeholder="Name, sequence, library, job or workstation…", key="database-browser-query")
        if view == "Molecules":
            grouping_label = grouping.selectbox("Group by", ["Whole molecule", "Individual sequence"], key="database-browser-molecule-grouping")
        else:
            grouping_label = grouping.selectbox("Group by", ["Campaigns + runs", "Individual runs"], key="database-browser-run-grouping")
        sort = ordering.selectbox("Order", ["Most jobs", "Recently added", "Name A–Z", "Most molecules"], key="database-browser-sort")
        index = build_browser_index(database, user, project_id=project_id, scope=scope, methods=methods,
                                    collection=collection, hide_smoke_test=hide_smoke)
        if view == "Molecules":
            subjects = molecular_subjects(index, individual_sequences=grouping_label == "Individual sequence")
            st.caption("One row per exact molecular composition or individual sequence across jobs. Different complexes, identity variants, models, targets and calculation settings remain separate records—not averaged or merged.")
        else:
            subjects = run_subjects(index, combine_campaigns=grouping_label == "Campaigns + runs")
            st.caption("Collections from the same recorded job appear together. Campaigns combine jobs only through explicit campaign IDs or parent-run links; names, shared targets and registration batches never imply a common campaign.")
        filtered = browse_subjects(subjects, index, query=query, sort=sort)
        overview = st.columns(3)
        overview[0].metric("Grouped items", f"{len(filtered):,}")
        overview[1].metric("Recorded jobs", f"{len(set().union(*(row['run_keys'] for row in filtered))):,}")
        overview[2].metric("Source collections", f"{len(set().union(*(row['group_ids'] for row in filtered))):,}")
        signature = hashlib.sha256(json.dumps([view, grouping_label, query, sort, project_id, scope, collection, methods, hide_smoke,
                                               [row["id"] for row in filtered]]).encode()).hexdigest()[:16]
        if st.session_state.get("database-browser-signature") != signature:
            st.session_state["database-browser-page"] = 1
            st.session_state["database-browser-signature"] = signature
        page_count = max(1, math.ceil(len(filtered) / 50))
        page_number = 1
        if page_count > 1:
            page_control, page_summary = st.columns([1, 6])
            page_number = int(page_control.number_input("Grouped page", min_value=1, max_value=page_count, step=1, key="database-browser-page"))
            page_summary.caption(f"Page {page_number} of {page_count} · 50 grouped items per page")
        page = filtered[(page_number - 1) * 50:page_number * 50]
        if page:
            rows = [{"Name": row["name"], "Jobs": len(row["run_keys"]),
                     "Molecule records": len(row["molecule_ids"]), "Collections": len(row["group_ids"]),
                     **({"Type": row["type"], "Sequence / chemistry": row["sequence"][:100] + ("…" if len(row["sequence"]) > 100 else "")} if view == "Molecules"
                        else {"Tools": ", ".join(sorted(row["tools"])), "Workstations": ", ".join(sorted(row["workstations"]))})}
                    for row in page]
            event = st.dataframe(rows, hide_index=True, use_container_width=True, on_select="rerun", selection_mode="single-row",
                                 key="database-browser-table-" + signature + str(page_number))
            st.caption(f"{len(filtered):,} grouped items · page {page_number} of {page_count}. Select a row to see its jobs and original records.")
            subject_lookup = {row["id"]: row for row in page}
            selected = st.selectbox("Inspect grouped item", ["", *subject_lookup],
                                    format_func=lambda value: subject_lookup[value]["name"] + " · " + value[:8] if value else "Select a table row or choose an item…",
                                    key="database-browser-inspect-" + signature + str(page_number))
            if not selected and event.selection.rows:
                row_number = event.selection.rows[0]
                if 0 <= row_number < len(page):
                    selected = page[row_number]["id"]
            if selected:
                _render_subject(database, user, index, subject_lookup[selected], project_id, scope)
        else:
            st.info("No matching grouped records. Try another search or adjust the access/source filters. Original collections remain available under Collections.")
    if st.session_state.get("molecular-import-open"):
        render_import_dialog(user, projects)


def _render_subject(database, user, index, subject, project_id, scope):
    from mn_molecular_design_studio.database_ui import render_details

    st.subheader(subject["name"])
    st.caption(f"{len(subject['run_keys']):,} jobs · {len(subject['molecule_ids']):,} original molecule records · {len(subject['group_ids']):,} collections. This is a navigation group, not a merged scientific result.")
    if sequence := subject.get("sequence"):
        with st.expander("Exact sequence / molecular contents"):
            st.code(sequence, language=None)
    run_rows = [index["runs"][key] for key in sorted(subject["run_keys"])]
    if run_rows:
        st.markdown("**Jobs contributing to this item**")
        st.dataframe([{"Name": row["name"], "Run": row["run_id"], "Task": row["task_group"] or "Unknown",
                       "Tools": ", ".join(sorted(row["tools"])), "Workstations": ", ".join(sorted(row["workstations"])),
                       "Node": row["node"], "Owner": row["owner_id"], "Project": row["project_id"] or "Personal"}
                      for row in run_rows], hide_index=True, use_container_width=True)
    molecule_ids = sorted(subject["molecule_ids"])
    molecule_id = st.selectbox("Original molecule record", molecule_ids,
                              format_func=lambda value: index["molecules"][value]["name"] + " · " + value,
                              key="database-browser-record-" + subject["id"])
    source_groups = sorted(index["molecules"][molecule_id]["group_ids"])
    group_id = st.selectbox("Original source collection", source_groups,
                           format_func=lambda value: index["groups"][value]["name"] + " · " + value[:12],
                           key="database-browser-source-" + subject["id"] + molecule_id)
    if st.button("Open original collection", key="database-browser-open-" + subject["id"]):
        database.group(group_id, user)
        st.session_state["molecular-preselect-group"] = group_id
        st.session_state["molecular-preselect-molecule"] = molecule_id
        st.session_state["molecular-preselect-collection"] = "all"
        st.rerun()
    if st.toggle("Load molecule inspector", key="database-browser-load-inspector-" + subject["id"] + molecule_id,
                 help="Load the selected original record, models and evaluation history. The grouped overview itself reads metadata only."):
        records = database.list_molecules(user, molecule_ids=[molecule_id], project_id=project_id,
                                         scope=scope, collection="all", limit=1)["items"]
        if records:
            render_details(database, user, records[0])
