"""Grouped Database browsing in the established filter/table/inspector layout."""
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
    from mn_molecular_design_studio.database_ui import (
        _database_view_selector, _projects, render_details, render_import_dialog, sync_workspace,
    )
    from mn_molecular_design_studio.molecule_table import render_table

    with st.container(key="molecular-database"):
        sidebar, main, details = st.columns([.25, .50, .25], gap="medium")
        navigation = sidebar.container(key="database-groups")
        with navigation:
            st.markdown("### Browse")
            _database_view_selector(view)
            new, refresh = st.columns([4, 1])
            if new.button("Import / new", icon=":material/add:", key="molecular-database-new", use_container_width=True):
                st.session_state["molecular-import-open"] = True
                st.session_state["molecular-import-step"] = 0
            if refresh.button("↻", key="molecular-database-refresh", help="Explicitly refresh saved job indexing"):
                st.session_state["molecular-indexing-errors"] = sync_workspace(user, database=database, force=True)
                st.rerun()
            query = st.text_input("Search grouped records", placeholder="Name, sequence, library, job…", key="database-browser-query")
            if view == "Molecules":
                grouping_label = st.selectbox("Group by", ["Whole molecule", "Individual sequence"], key="database-browser-molecule-grouping")
            else:
                grouping_label = st.selectbox("Group by", ["Campaigns + runs", "Individual runs"], key="database-browser-run-grouping")
            sort = st.selectbox("Order", ["Most jobs", "Recently added", "Name A–Z", "Most molecules"], key="database-browser-sort")
            with st.popover("Access & source filters", use_container_width=True):
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
            st.caption(f"{collection_label} · {scope_label} access · {project_name}")
            with st.expander("How grouping works"):
                st.caption("Grouping changes navigation only. Exact sequences and recorded molecular contents link jobs without merging identities, structures or scores. Campaigns require explicit links. Related history is included only where you have access.")
            index = build_browser_index(database, user, project_id=project_id, scope=scope, methods=methods,
                                        collection=collection, hide_smoke_test=hide_smoke)
            subjects = (molecular_subjects(index, individual_sequences=grouping_label == "Individual sequence")
                        if view == "Molecules" else run_subjects(index, combine_campaigns=grouping_label == "Campaigns + runs"))
            filtered = browse_subjects(subjects, index, query=query, sort=sort)
            signature = hashlib.sha256(json.dumps([view, grouping_label, query, sort, project_id, scope, collection, methods, hide_smoke,
                                                   [row["id"] for row in filtered]]).encode()).hexdigest()[:16]
            if st.session_state.get("database-browser-signature") != signature:
                st.session_state["database-browser-page"] = 1
                st.session_state["database-browser-signature"] = signature
            page_count = max(1, math.ceil(len(filtered) / 50))
            page_number = 1
            if page_count > 1:
                page_number = int(st.number_input("Grouped page", min_value=1, max_value=page_count, step=1, key="database-browser-page"))
            st.caption(f"{len(filtered):,} grouped items · page {page_number} of {page_count}")
            page = filtered[(page_number - 1) * 50:page_number * 50]
            subject_lookup = {row["id"]: row for row in page}
            picker_key = "database-browser-inspect-" + signature + str(page_number)
            pending = st.session_state.pop("database-browser-pending-subject", None)
            if pending and pending.get("picker_key") == picker_key and pending.get("subject_id") in subject_lookup:
                st.session_state[picker_key] = pending["subject_id"]
            selected = ""
            if page:
                selected = st.selectbox("Inspect grouped item", ["", *subject_lookup],
                                        format_func=lambda value: subject_lookup[value]["name"] + " · " + value[:8] if value else "First item on this page",
                                        key=picker_key) or page[0]["id"]
        subject = subject_lookup.get(selected)
        with main:
            st.markdown("### " + ("Molecules" if view == "Molecules" else "Runs & campaigns"))
            st.caption(f"{len(filtered):,} grouped items · {len(index['runs']):,} recorded jobs. Select a row to update the molecule viewer.")
            if page:
                rows = [{"id": row["id"], "values": {
                    "Name": row["name"], "Jobs": len(row["run_keys"]), "Records": len(row["molecule_ids"]),
                    "Collections": len(row["group_ids"]),
                    **({"Type": row["type"], "Sequence / chemistry": row["sequence"][:100] + ("…" if len(row["sequence"]) > 100 else "")}
                       if view == "Molecules" else {"Tools": ", ".join(sorted(row["tools"])), "Workstations": ", ".join(sorted(row["workstations"]))})}}
                        for row in page]
                event = render_table(rows, [{"key": field, "label": field, "kind": "value"} for field in rows[0]["values"]],
                                     selected_ids=[], focus_id=selected, key="database-browser-table-" + signature + str(page_number), allow_selection=False)
                event_key = "database-browser-event-" + signature + str(page_number)
                if event and event.get("token") != st.session_state.get(event_key):
                    st.session_state[event_key] = event.get("token")
                    if event.get("focusId") in subject_lookup and event["focusId"] != selected:
                        st.session_state["database-browser-pending-subject"] = {"picker_key": picker_key, "subject_id": event["focusId"]}
                        st.rerun()
                _render_subject_summary(index, subject)
            else:
                st.info("No matching grouped records. Try another search or adjust the access/source filters. Original collections remain available under Collections.")
        molecule_id = ""
        group_id = ""
        if subject:
            with navigation:
                molecule_ids = sorted(subject["molecule_ids"])
                molecule_id = st.selectbox("Original molecule record", molecule_ids,
                                          format_func=lambda value: index["molecules"][value]["name"] + " · " + value[-8:],
                                          key="database-browser-record-" + subject["id"])
                source_groups = sorted(index["molecules"][molecule_id]["group_ids"] & subject["group_ids"])
                group_id = st.selectbox("Original source collection", source_groups,
                                       format_func=lambda value: index["groups"][value]["name"] + " · " + value[:12],
                                       key="database-browser-source-" + subject["id"] + molecule_id)
                if st.button("Open original collection", key="database-browser-open-" + subject["id"]):
                    database.group(group_id, user)
                    st.session_state["molecular-preselect-group"] = group_id
                    st.session_state["molecular-preselect-molecule"] = molecule_id
                    st.session_state["molecular-preselect-collection"] = "all"
                    st.rerun()
        with details, st.container(key="database-details"):
            if molecule_id:
                records = database.list_molecules(user, group_id=group_id, molecule_ids=[molecule_id], project_id=project_id,
                                                 scope=scope, collection="all", hide_smoke_test=hide_smoke, limit=1)["items"]
                if records:
                    st.session_state["database-focus"] = molecule_id
                    render_details(database, user, records[0])
                else:
                    st.caption("This record is no longer accessible. Choose another item or refresh the view.")
            else:
                st.caption("Select a molecule to inspect its structure, chains, and scores.")
    if st.session_state.get("molecular-import-open"):
        render_import_dialog(user, projects)


def _render_subject_summary(index, subject):
    st.caption(f"{len(subject['run_keys']):,} jobs · {len(subject['molecule_ids']):,} original molecule records · {len(subject['group_ids']):,} collections. Scores and structures remain attached to their original evaluations.")
    with st.expander("Jobs and provenance", expanded=False):
        if sequence := subject.get("sequence"):
            st.code(sequence, language=None)
        run_rows = [index["runs"][key] for key in sorted(subject["run_keys"])]
        if run_rows:
            st.dataframe([{"Name": row["name"], "Run": row["run_id"], "Task": row["task_group"] or "Unknown",
                           "Tools": ", ".join(sorted(row["tools"])), "Workstations": ", ".join(sorted(row["workstations"])),
                           "Node": row["node"], "Owner": row["owner_id"], "Project": row["project_id"] or "Personal"}
                          for row in run_rows], hide_index=True, use_container_width=True)
