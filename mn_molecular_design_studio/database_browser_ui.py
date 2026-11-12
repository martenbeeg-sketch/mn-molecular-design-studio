"""Choose a molecule/run on the left, its evidence centrally, and inspect on the right."""
from __future__ import annotations

import hashlib
import json
import math

import streamlit as st

from mn_molecular_design_studio.database_browser import (
    browse_subjects,
    build_browser_index,
    component_subjects,
    load_result_record,
    molecular_subjects,
    run_subjects,
    subject_results,
)
from mn_molecular_design_studio.molecular_database import IMPORT_METHODS


def _signature(values):
    return hashlib.sha256(json.dumps(values).encode()).hexdigest()[:16]


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
            st.caption("Choose on the left → results in the middle → viewer on the right.")
            new, refresh = st.columns([4, 1])
            if new.button("Import / new", icon=":material/add:", key="molecular-database-new", use_container_width=True):
                st.session_state["molecular-import-open"] = True
                st.session_state["molecular-import-step"] = 0
            if refresh.button("↻", key="molecular-database-refresh", help="Explicitly refresh saved job indexing"):
                st.session_state["molecular-indexing-errors"] = sync_workspace(user, database=database, force=True)
                st.rerun()
            query = st.text_input("Find a molecule" if view == "Molecules" else "Find a run / campaign",
                                  placeholder="Name, sequence, library, job…", key="database-browser-query")
            chooser = st.container()
            with st.popover("Browse & access filters", use_container_width=True):
                if view == "Molecules":
                    grouping_label = st.selectbox("Group by", ["Ligands & targets", "Whole molecule", "Individual sequence"],
                                                  key="database-browser-molecule-grouping-v2")
                else:
                    grouping_label = st.selectbox("Group by", ["Campaigns + runs", "Individual runs"], key="database-browser-run-grouping")
                sort = st.selectbox("Order", ["Most jobs", "Recently added", "Name A–Z", "Most molecules"], key="database-browser-sort")
                collection_label = st.segmented_control("Starting collection", ["Library", "All indexed"], default="Library",
                                                        key="molecular-collection-view") or "Library"
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
            with st.expander("Molecules, runs or collections?"):
                st.caption("Molecules: one exact ligand or polymer, including complexes that contain it and designs with a recorded target link. Whole molecule keeps complete complexes separate.")
                st.caption("Runs: results of one recorded execution. Campaigns combine only explicitly linked executions, not similarly named jobs.")
                st.caption("Collections: imported libraries or folders of selected molecules. A molecule can belong to several folders; a folder is not a new calculation.")
                st.caption("Only indexed, accessible evidence is shown. Missing historical links are not guessed from names. No identities, structures or scores are merged.")
            index = build_browser_index(database, user, project_id=project_id, scope=scope, methods=methods,
                                        collection=collection, hide_smoke_test=hide_smoke, include_related=True)
            if view == "Runs":
                subjects = run_subjects(index, combine_campaigns=grouping_label == "Campaigns + runs")
            elif grouping_label == "Ligands & targets":
                subjects = component_subjects(index)
            else:
                subjects = molecular_subjects(index, individual_sequences=grouping_label == "Individual sequence")
            if view == "Molecules":
                types = sorted({row["type"] for row in subjects})
                type_filter = st.selectbox("Molecule type", ["All types", *types], key="database-browser-entity-type-" + grouping_label)
                if type_filter != "All types":
                    subjects = [row for row in subjects if row["type"] == type_filter]
            filtered = browse_subjects(subjects, index, query=query, sort=sort)
            signature = _signature([view, grouping_label, query, sort, project_id, scope, collection, methods, hide_smoke,
                                    [row["id"] for row in filtered]])
            if st.session_state.get("database-browser-signature") != signature:
                st.session_state["database-browser-page"] = 1
                st.session_state["database-browser-signature"] = signature
            page_count = max(1, math.ceil(len(filtered) / 50))
            page_number = 1
            with chooser:
                if page_count > 1:
                    page_number = int(st.number_input("Browse page", min_value=1, max_value=page_count, step=1, key="database-browser-page"))
                st.caption(f"{len(filtered):,} {'molecules' if view == 'Molecules' else 'runs / campaigns'} · page {page_number} of {page_count}")
                page = filtered[(page_number - 1) * 50:page_number * 50]
                subject_lookup = {row["id"]: row for row in page}
                selected = ""
                if page:
                    with st.container(height=200, border=False):
                        selected = st.radio("Choose a molecule" if view == "Molecules" else "Choose a run / campaign", list(subject_lookup),
                            format_func=lambda value: subject_lookup[value]["name"] + f" · {len(subject_lookup[value]['run_keys']):,} jobs",
                            key="database-browser-subject-" + signature + str(page_number))
        subject = subject_lookup.get(selected)
        selected_result = None
        result_page = []
        result_signature = ""
        if subject:
            results = subject_results(database, index, subject, run_view=view == "Runs")
            with navigation:
                st.markdown("#### Filter connected results")
                result_query = st.text_input("Find a result", placeholder="Name, tool, job, model…",
                                             key="database-results-query-" + subject["id"])
                kinds = sorted({row["kind_label"] for row in results})
                kind = st.selectbox("Result type", ["All result types", *kinds], key="database-results-kind-" + subject["id"])
                tools = sorted({row["tool"] for row in results if row["tool"]})
                tool = st.selectbox("Result tool", ["All tools", *tools], key="database-results-tool-" + subject["id"])
                terms = result_query.casefold().split()
                results = [row for row in results if (kind == "All result types" or row["kind_label"] == kind)
                           and (tool == "All tools" or row["tool"] == tool)
                           and all(term in " ".join(str(value) for value in row.values()).casefold() for term in terms)]
                result_signature = _signature([signature, subject["id"], result_query, kind, tool, [row["id"] for row in results]])
                page_key = "database-results-page-" + result_signature
                result_page_count = max(1, math.ceil(len(results) / 50))
                result_page_number = 1
                if result_page_count > 1:
                    result_page_number = int(st.number_input("Results page", min_value=1, max_value=result_page_count, step=1, key=page_key))
                st.caption(f"{len(results):,} connected results · page {result_page_number} of {result_page_count}")
                result_page = results[(result_page_number - 1) * 50:result_page_number * 50]
                result_lookup = {row["id"]: row for row in result_page}
                result_picker_key = "database-results-inspect-" + result_signature + str(result_page_number)
                pending = st.session_state.pop("database-results-pending", None)
                if pending and pending.get("picker_key") == result_picker_key and pending.get("result_id") in result_lookup:
                    st.session_state[result_picker_key] = pending["result_id"]
                if result_page:
                    focused_result = st.selectbox("Inspect result", list(result_lookup),
                        index=next((offset for offset, row in enumerate(result_page) if row["structure_ids"]), 0),
                        format_func=lambda value: result_lookup[value]["name"] + " · " + result_lookup[value]["kind_label"]
                        + " · " + (result_lookup[value]["run_id"] or result_lookup[value]["collection"]) + " · " + value[-6:],
                        key=result_picker_key)
                    selected_result = result_lookup[focused_result]
                    if st.button("Open original collection", key="database-browser-open-" + subject["id"]):
                        database.group(selected_result["group_id"], user)
                        st.session_state["molecular-preselect-group"] = selected_result["group_id"]
                        st.session_state["molecular-preselect-molecule"] = selected_result["molecule_id"]
                        st.session_state["molecular-preselect-collection"] = "all"
                        st.rerun()
        with main:
            st.markdown("### Connected results")
            if subject:
                st.markdown("**" + subject["name"] + "**")
                st.caption(f"{len(results):,} results · {len(subject['run_keys']):,} recorded jobs. Select a result below to update the viewer.")
                if result_page:
                    rows = [{"id": row["id"], "values": {
                        "Name": row["name"], "Result": row["kind_label"], "Tool": row["tool"], "Job": row["run_id"],
                        "Candidate": row.get("candidate", ""), "Model / pose": row["model"], "Structures": len(row["structure_ids"]),
                        "Status": row["status"], "Collection": row["collection"],
                        "Workstation": row["workstation"], "Connection": row["relation"],
                    }} for row in result_page]
                    table_key = "database-results-table-" + result_signature + str(result_page_number)
                    event = render_table(rows, [{"key": field, "label": field, "kind": "value"} for field in rows[0]["values"]],
                        selected_ids=[], focus_id=selected_result["id"], key=table_key, allow_selection=False)
                    event_key = table_key + "-event"
                    if event and event.get("token") != st.session_state.get(event_key):
                        st.session_state[event_key] = event.get("token")
                        if event.get("focusId") in result_lookup and event["focusId"] != selected_result["id"]:
                            st.session_state["database-results-pending"] = {"picker_key": result_picker_key, "result_id": event["focusId"]}
                            st.rerun()
                else:
                    st.info("No matching connected results. Adjust the result filters on the left.")
                _render_subject_summary(index, subject)
            else:
                st.info("No matching molecules or runs. Try another search or adjust the access/source filters.")
        with details, st.container(key="database-details"):
            if selected_result:
                record = load_result_record(database, user, selected_result, project_id=project_id, scope=scope, hide_smoke_test=hide_smoke)
                if record:
                    st.session_state["database-focus"] = record["id"]
                    if selected_result["evaluation_id"]:
                        st.session_state["database-evaluation-" + record["id"]] = selected_result["evaluation_id"]
                    st.caption(selected_result["kind_label"] + " · " + (selected_result["tool"] or selected_result["collection"]))
                    if not record["structures"]:
                        st.caption("No structure is recorded for this result. Choose a structure-bearing result in the middle for 3D viewing.")
                    render_details(database, user, record)
                else:
                    st.caption("This result is no longer accessible. Choose another result or refresh the view.")
            else:
                st.caption("Select a result to inspect its structure, chains, and scores.")
    if st.session_state.get("molecular-import-open"):
        render_import_dialog(user, projects)


def _render_subject_summary(index, subject):
    with st.expander("Recorded connections & provenance", expanded=False):
        st.caption(f"{len(subject['run_keys']):,} jobs · {len(subject['molecule_ids']):,} original molecule records · {len(subject['group_ids']):,} collections. Each result keeps its own evaluation and scores.")
        if sequence := subject.get("sequence"):
            st.code(sequence, language=None)
        run_rows = [index["runs"][key] for key in sorted(subject["run_keys"])]
        if run_rows:
            st.dataframe([{"Name": row["name"], "Run": row["run_id"], "Task": row["task_group"] or "Unknown",
                           "Tools": ", ".join(sorted(row["tools"])), "Workstations": ", ".join(sorted(row["workstations"])),
                           "Node": row["node"], "Owner": row["owner_id"], "Project": row["project_id"] or "Personal"}
                          for row in run_rows], hide_index=True, use_container_width=True)
