"""Reusable input datasets: browse members, import, annotate and reference projects."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.result_details_ui import render_result_details
from mn_molecular_design_studio import datasets
from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
from mn_molecular_design_studio.database_actions_ui import ERRORS, read_only_workspace, render_library_actions, render_selection_actions
from mn_molecular_design_studio.folder_tree import render_tree
from mn_molecular_design_studio.navigation import open_dataset_results, open_project_database, request_dataset_import


def _dismiss_references():
    st.session_state.pop("datasets-project-editor", None)


@st.dialog("Dataset project references", width="large", on_dismiss=_dismiss_references)
def _reference_dialog(database, user, projects, request):
    group = request["group"]
    known = {project["id"]: project["name"] for project in projects if not project.get("deleted_at")}
    existing = datasets.project_references(group)
    unavailable = set(existing) - set(known)
    st.caption("One dataset can be referenced by several projects. Molecules are not copied. "
        "References do not grant access or transfer ownership; existing private/project-sharing rules remain in force.")
    if unavailable:
        st.warning("Some existing project references are unavailable to this account. An authorized administrator must review those links.")
    discard_unavailable = st.checkbox("Remove the unavailable project references", key="datasets-remove-unavailable-references",
        help="Explicitly detach missing or inaccessible project references without removing the dataset or its members.") if unavailable else False
    selected = st.multiselect("Used by projects", list(known), default=[value for value in existing if value in known],
        format_func=known.get, key="datasets-project-links-" + request["revision"])
    if group["project_id"]:
        st.info("This older dataset still has its original project-based access. "
            "Additional references do not change that access. Cleanup of its source project is blocked while another project references it.")
    if st.button("Save project references", type="primary", disabled=(bool(unavailable) and not discard_unavailable) or read_only_workspace(), key="datasets-save-project-links"):
        try:
            datasets.save_project_references(database, user, group["id"], selected, expected_revision=request["revision"])
        except ERRORS as exc:
            st.error(str(exc))
        else:
            _dismiss_references()
            st.session_state["properties-saved-notice"] = "Saved project references. Ownership, access and scientific data are unchanged."
            st.session_state.pop("projects-inventory", None)
            st.rerun(scope="app")


def dataset_tree(groups, *, by_project=False, projects=()):
    names = {project["id"]: project["name"] for project in projects}
    buckets, choices = {}, {}
    for group in groups:
        categories = ([(project_id, names.get(project_id, "Unavailable project")) for project_id in group["project_references"]]
            or [("independent", "Independent libraries")]) if by_project else [(group["category"], group["category"])]
        for category_id, category in categories:
            root_id = "dataset-category:" + hashlib.sha256(category_id.encode()).hexdigest()[:16]
            identifier = root_id + ":" + group["id"]
            root = buckets.setdefault(root_id, {"id": root_id, "label": category, "children": []})
            root["children"].append({"id": identifier, "label": group.get("display_label", group["name"]), "count": group["molecule_count"]})
            choices[identifier] = group["id"]
    return sorted(buckets.values(), key=lambda node: node["label"].casefold()), choices


def render_datasets(user, projects, *, database=None):
    if database is None:
        configured = os.environ.get("MDS_DATA_BROWSER_DATABASE_ROOT", "").strip()
        database = ReadOnlyLibrary(Path(configured) if configured else None)
    read_only = read_only_workspace()
    st.title("Datasets")
    st.caption("Import reusable molecular inputs or save a selection from Data browser for another study.")
    from mn_molecular_design_studio.entity_properties_ui import render_folder_properties, render_notice, render_record_properties
    render_notice()
    if return_to := st.session_state.get("dataset-import-return"):
        if st.button("Return to setup", key="datasets-return-to-setup"):
            st.session_state["page"] = return_to
            st.session_state.pop("dataset-import-return", None)
            st.session_state.pop("molecular-import-open", None)
            st.rerun()
    columns = st.columns(4)
    for column, label, kind, icon in zip(columns, ("Import proteins", "Import compounds", "Import sequences", "Other formats"),
            ("protein", "compound", "sequence", ""), (":material/biotech:", ":material/science:", ":material/notes:", ":material/upload:")):
        if column.button(label, icon=icon, disabled=read_only, width="stretch", key=f"datasets-import{('-' + kind) if kind else ''}"):
            request_dataset_import(st.session_state, kind, project_id=st.session_state.get("datasets-project", ""))
    new_column, results_column, _ = st.columns([1, 2, 2])
    if new_column.button("New dataset", icon=":material/add:", disabled=read_only, key="datasets-new"):
        st.session_state["data-browser-management"] = {"kind": "group", "ids": [],
            "project_id": st.session_state.get("datasets-project", "")}
    if results_column.button("Start from saved results", icon=":material/folder_copy:", key="datasets-open-browser"):
        open_project_database(st.session_state, st.session_state.get("datasets-project", ""))
        st.session_state["dataset-reuse-request"] = True
        st.rerun()
    render_library_actions(database, user, projects, dataset_mode=True, reference_project=st.session_state.get("datasets-project", ""))
    if not database.path.is_file():
        st.info("No datasets yet. Import a molecule library, sequence collection or structure, without creating a project first.")
        return
    try:
        groups = datasets.catalog(database, user)
    except ERRORS as exc:
        st.error("Datasets could not be loaded: " + str(exc))
        return
    requested = st.session_state.pop("datasets-open-group", "")
    if requested:
        if requested in {group["id"] for group in groups}:
            st.session_state["datasets-selected"] = requested
            st.session_state["datasets-query"] = ""
            st.session_state["datasets-member-query"] = ""
        else:
            st.warning("The requested dataset is no longer accessible. No other dataset has been substituted.")
            st.session_state["datasets-selected"] = ""
    if not groups:
        st.info("Import inputs above, or start from a selection in Data browser.")
        return
    labels = datasets.display_labels(groups, projects)
    groups = [{**group, "display_label": labels[group["id"]]} for group in groups]
    left, middle, right = st.columns([.24, .48, .28], gap="medium")
    project_names = {project["id"]: project["name"] for project in projects if not project.get("deleted_at")}
    with left:
        st.markdown("### Browse datasets")
        query = st.text_input("Find a dataset", placeholder="Library name or alias…", key="datasets-query").casefold().strip()
        options = ["", *project_names]
        if st.session_state.get("datasets-project") not in options:
            st.session_state["datasets-project"] = ""
        project = st.selectbox("Project use", options, format_func=lambda value: project_names.get(value, "All datasets — no project required"),
            key="datasets-project")
        by_project = st.selectbox("Organise folders", ["Input type", "Project references"], key="datasets-tree-order") == "Project references"
        visible = [group for group in groups if (not project or project in group["project_references"]) and all(term in
            " ".join([group["name"], *group.get("name_aliases", [])]).casefold() for term in query.split())]
        lookup = {group["id"]: group for group in visible}
        selected = st.session_state.get("datasets-selected", "")
        if selected not in lookup:
            selected = next(iter(lookup), "") if not requested else ""
        st.session_state["datasets-selected"] = selected
        nodes, choices = dataset_tree(visible, by_project=by_project, projects=projects)
        selected_node = next((identifier for identifier, group_id in choices.items() if group_id == selected), "")
        event = render_tree(nodes, selected_id=selected_node, expanded_ids=[node["id"] for node in nodes], key="datasets-tree")
        if event and event.get("token") != st.session_state.get("datasets-tree-event"):
            st.session_state["datasets-tree-event"] = event.get("token")
            if event.get("selectedId") in choices:
                st.session_state["datasets-selected"] = choices[event["selectedId"]]
                st.session_state["datasets-member-query"] = ""
                st.session_state.pop("datasets-page-selection", None)
                st.rerun()
        st.caption(f"{len(visible):,} datasets")
    if not selected:
        with middle:
            st.info("No dataset matches this selection. Clear the dataset search or choose All datasets.")
        return
    group = lookup[selected]
    view = datasets.DatasetView(database, selected)
    with middle:
        st.markdown("### " + group["name"])
        access = "Private to its owner" if not group["project_id"] else "Shared through its original project: " + project_names.get(group["project_id"], "Project")
        st.caption(f"{group['molecule_count']:,} input members · {access}")
        scope = {"id": "dataset:" + selected, "label": group["name"], "entity": {"kind": "dataset", "group_id": selected}}
        render_folder_properties(database, user, {}, scope, projects, read_only=read_only, key_prefix="datasets-folder")
        labels = [project_names.get(project_id, "Unavailable project") for project_id in group["project_references"]]
        st.caption("Referenced by projects: " + (", ".join(labels) or "None — reusable independent library"))
        if st.button("Edit project references", key="datasets-project-references", disabled=read_only or group["read_only"] or (
                user.get("role") != "admin" and group["owner_id"] != user["id"])):
            snapshot = database.group(selected, user)
            st.session_state["datasets-project-editor"] = {"group": snapshot, "revision": datasets.source_revision(snapshot)}
        request = st.session_state.get("datasets-project-editor")
        if request and request["group"]["id"] == selected:
            _reference_dialog(database, user, projects, request)
        member_query = st.text_input("Find a member", placeholder="Compound ID, name, sequence or SMILES…", key="datasets-member-query")
        signature = hashlib.sha256(json.dumps([selected, member_query]).encode()).hexdigest()[:20]
        page_key = "datasets-member-page-" + signature
        page = int(st.number_input("Member page", min_value=1, step=1, key=page_key))
        listing = view.list_molecules(user, query=member_query, limit=100, offset=(page - 1) * 100)
        records = listing["items"]
        context = (selected, member_query, page)
        retained = st.session_state.get("datasets-page-selection") or {}
        selected_ids = retained.get("ids", []) if retained.get("context") == context else []
        table = st.dataframe([{"Name": record["name"], "Type": record["type"].replace("small_molecule", "Small molecule"),
            "Sequence / SMILES": " : ".join(component.get("sequence") or component.get("smiles") or component.get("ccd", "") for component in record["components"]),
            "Structures": len(record["structures"])} for record in records], hide_index=True, width="stretch", height=530,
            on_select="rerun", selection_mode="multi-row", selection_default={"selection": {"rows": [position for position, record in enumerate(records) if record["id"] in selected_ids]}},
            key="datasets-members-" + signature + str(page))
        selected_ids = [records[position]["id"] for position in table.selection.rows if 0 <= position < len(records)]
        st.session_state["datasets-page-selection"] = {"context": context, "ids": selected_ids}
        st.caption(f"Page {page} · {len(records):,} input members shown. " + ("More members are available on the next page." if listing["next_offset"] is not None else "End of this member search."))
        focus = st.session_state.pop("datasets-focus-molecule", "") or (selected_ids[-1] if selected_ids else records[0]["id"] if records else "")
        chosen = [{"id": identifier, "result": {"group_id": selected, "molecule_id": identifier}} for identifier in selected_ids]
        render_selection_actions(view, user, projects, chosen)
        if st.button("Open connected results in Data browser", key="datasets-connected-results"):
            open_dataset_results(st.session_state, selected)
            st.rerun()
    with right:
        if not focus:
            st.caption("Select an input member to inspect its sequence, chemistry and saved coordinates.")
        else:
            record = next(iter(view.list_molecules(user, molecule_ids=[focus], limit=1)["items"]), None)
            if record:
                st.markdown("### " + record["name"])
                render_record_properties(database, user, record, read_only=read_only, key_prefix="datasets-member")
                for component in record["components"]:
                    st.caption(component["chain"] + " · " + component["type"].replace("small_molecule", "Small molecule"))
                    if component.get("smiles"):
                        from mn_molecular_design_studio.cofolding_input import _smiles_preview
                        preview = _smiles_preview(component["smiles"])
                        if preview:
                            st.image(preview)
                    st.code(component.get("sequence") or component.get("smiles") or component.get("ccd", ""), wrap_lines=True)
                if record["structures"]:
                    from mn_molecular_design_studio.data_browser_ui import _render_structure
                    _render_structure(view, user, record, "")
                else:
                    st.caption("No coordinates are stored in this input dataset. This does not mean that no calculated structures exist elsewhere.")
                if st.button("Open this member’s results in Data browser", key="datasets-member-results"):
                    open_dataset_results(st.session_state, selected, record["id"])
                    st.rerun()
                with st.expander("Input annotations & identifiers"):
                    render_result_details(st, record["metadata"])
