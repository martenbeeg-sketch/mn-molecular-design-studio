"""Explicit library management alongside the scientific Database browser."""
from __future__ import annotations

import hashlib
import os
import sqlite3

import streamlit as st

from mn_molecular_design_studio.molecular_database import MolecularDatabase


ERRORS = (ValueError, PermissionError, OSError, RuntimeError, sqlite3.Error)


def read_only_workspace():
    return (os.environ.get("MDS_READ_ONLY_PREVIEW") == "1"
            or os.environ.get("MDS_API_BASE_URL", "").rstrip("/") == "http://127.0.0.1:1")


def _dismiss_actions():
    st.session_state.pop("data-browser-management", None)


def _open_collection(group_id, *, historical=False):
    st.session_state.update({"molecular-preselect-group": group_id,
        "molecular-preselect-collection": "all" if historical else "library"})
    _dismiss_actions()
    st.rerun(scope="app")


@st.dialog("Manage collections", width="large", on_dismiss=_dismiss_actions)
def _collections_dialog(database, user, projects):
    historical = st.checkbox("Include historical indexed collections", key="data-browser-manage-history")
    groups = database.groups(user, collection="all" if historical else "library")
    groups.sort(key=lambda group: (group["name"].casefold(), group["id"]))
    if not groups:
        st.info("No accessible collections. Import molecules or create a collection from Library actions.")
        return
    lookup = {group["id"]: group for group in groups}
    if st.session_state.get("data-browser-manage-group") not in lookup:
        st.session_state["data-browser-manage-group"] = next(iter(lookup))
    project_names = {project["id"]: project["name"] for project in projects}
    group_id = st.selectbox("Collection", list(lookup), key="data-browser-manage-group",
        format_func=lambda value: lookup[value]["name"] + " · " + project_names.get(lookup[value]["project_id"], "Personal")
            + " · " + value[:8])
    group = lookup[group_id]
    st.caption(f"{group['molecule_count']:,} molecules · membership and measurements remain attached to this source collection.")
    if st.button("Open collection", type="primary", key="data-browser-manage-open"):
        _open_collection(group_id, historical=historical)
    can_edit = user.get("role") == "admin" or group["owner_id"] == user["id"]
    if st.button("Add metadata column", disabled=not can_edit or read_only_workspace() or bool(group["origin"].get("_replica_node")),
            key="data-browser-manage-column", help="Metadata columns are edited on the collection’s source workstation."):
        st.session_state["data-browser-management"] = {"kind": "column", "group_id": group_id}
        st.rerun(scope="app")
    page = int(st.number_input("Collection page", min_value=1, max_value=max(1, (group["molecule_count"] + 99) // 100),
        value=1, key="data-browser-manage-page-" + group_id))
    records = database.list_molecules(user, group_id=group_id, collection="all", limit=100, offset=(page - 1) * 100)["items"]
    rows_key = "data-browser-manage-selection-" + group_id + "-" + str(page)
    selected_ids = st.session_state.get(rows_key, [])
    retained = [position for position, record in enumerate(records) if record["id"] in selected_ids]
    snapshot = hashlib.sha256(repr([record["id"] for record in records]).encode()).hexdigest()[:12]
    event = st.dataframe([{"Name": record["name"], "Type": record["type"], "ID": record["id"]} for record in records],
        hide_index=True, width="stretch", on_select="rerun", selection_mode="multi-row",
        selection_default={"selection": {"rows": retained}}, key=f"data-browser-manage-rows-{group_id}-{page}-{snapshot}")
    selected = [records[position]["id"] for position in event.selection.rows if 0 <= position < len(records)]
    st.session_state[rows_key] = selected
    st.caption("Removal affects only this collection and its attached structures, scores and provenance. "
        "Other collections and original run folders are not removed. Use Projects for whole-run or project cleanup.")
    selection_key = hashlib.sha256(repr((group_id, sorted(selected))).encode()).hexdigest()[:20]
    confirmed = st.checkbox("Remove the selected entries and their data from this collection", key="data-browser-manage-confirm-" + selection_key)
    if st.button("Remove selected entries", disabled=not selected or not confirmed or not can_edit or read_only_workspace(),
            key="data-browser-manage-remove"):
        try:
            count = MolecularDatabase(database.root).remove_from_group(user, group_id, selected)
        except ERRORS as exc:
            st.error(str(exc))
        else:
            st.session_state["molecular-database-flash"] = f"Removed {count:,} entries from “{group['name']}”."
            _dismiss_actions()
            st.rerun(scope="app")


def render_library_actions(database, user, projects):
    from mn_molecular_design_studio.database_ui import render_column_dialog, render_group_dialog, render_import_dialog

    read_only = read_only_workspace()
    with st.popover("Library actions", icon=":material/settings:"):
        if read_only:
            st.caption("Read-only preview: importing and collection changes are disabled.")
        if st.button("Import molecules", key="data-browser-import", disabled=read_only):
            st.session_state.update({"molecular-import-open": True, "molecular-import-step": 0})
        if st.button("New collection", key="data-browser-new-collection", disabled=read_only):
            st.session_state["data-browser-management"] = {"kind": "group", "ids": []}
        if st.button("Manage collections", key="data-browser-manage-collections", disabled=not database.path.is_file()):
            st.session_state["data-browser-management"] = {"kind": "collections"}
        st.checkbox("Include historical indexed data", key="data-browser-indexed-history",
            help="Advanced inspection of saved SQL history. This does not admit unreviewed designs or import run folders.")
    if message := st.session_state.pop("molecular-database-flash", None):
        st.success(message)
    if st.session_state.get("molecular-import-open"):
        if not read_only:
            render_import_dialog(user, projects, database=MolecularDatabase(database.root))
        return
    request = st.session_state.get("data-browser-management") or {}
    if request.get("kind") == "collections":
        _collections_dialog(database, user, projects)
    elif not read_only and request.get("kind") in {"group", "column"}:
        if request["kind"] == "group":
            render_group_dialog(user, projects, request["ids"], request.get("project_id", ""),
                database=MolecularDatabase(database.root))
        else:
            render_column_dialog(user, request["group_id"], database=MolecularDatabase(database.root))


def _dismiss_selection():
    st.session_state.pop("data-browser-selection-dialog", None)


@st.dialog("Selected data", width="large", on_dismiss=_dismiss_selection)
def _selection_dialog(database, user, projects, entries, project_id):
    from mn_molecular_design_studio.database_ui import render_use_selection
    from mn_molecular_design_studio.molecule_exports import export_molecules

    ids = list(dict.fromkeys(entry["result"]["molecule_id"] for entry in entries))
    st.caption(f"{len(entries):,} saved results · {len(ids):,} distinct molecules. No admission decisions are changed.")
    read_only = read_only_workspace()
    if st.button("Create collection from these molecules", disabled=read_only, key="data-browser-collect-selection"):
        st.session_state["data-browser-management"] = {"kind": "group", "ids": ids, "project_id": project_id}
        _dismiss_selection()
        st.rerun(scope="app")
    if not read_only:
        render_use_selection(database, user, projects, ids, project_id)
    export_format = st.selectbox("Molecule export format", ["CSV", "FASTA", "SMILES", "SDF", "Structures ZIP"],
        key="data-browser-molecule-export-format")
    if st.button("Prepare molecule download", key="data-browser-molecule-export-prepare"):
        try:
            content, suffix, mime = export_molecules(database, user, ids, export_format)
        except ERRORS as exc:
            st.error(str(exc))
        else:
            st.download_button("Download " + export_format, content, file_name="selected-molecules." + suffix,
                mime=mime, key="data-browser-molecule-export-download")
    st.caption("Molecule exports include accessible saved molecular data. For exact displayed results and their linked analyses, use Table columns & CSV export.")
    if st.toggle("Compare these saved results", key="data-browser-compare-selection"):
        from mn_molecular_design_studio.comparison_ui import render_comparison
        from mn_molecular_design_studio.data_browser import load_saved_record

        records = [record for entry in entries if (record := load_saved_record(database, user, entry["result"])) is not None]
        render_comparison(records, key="data-browser-comparison", user=user, project_id=project_id,
            database=None if read_only else MolecularDatabase(database.root))


def render_selection_actions(database, user, projects, selected_entries, *, project_id=""):
    signature = (str(database.path.resolve()), user["id"], user.get("role"), tuple(
        (entry.get("id"), entry.get("structure_id"), *(entry["result"].get(field) for field in
            ("group_id", "molecule_id", "evaluation_id", "observation_id"))) for entry in selected_entries))
    if st.session_state.get("data-browser-selection-dialog") != signature:
        _dismiss_selection()
    if st.button("Data actions", icon=":material/more_horiz:", disabled=not selected_entries,
            help="Select table rows to collect, export, compare or reuse their molecules.", key="data-browser-selection-actions"):
        st.session_state["data-browser-selection-dialog"] = signature
    if selected_entries and st.session_state.get("data-browser-selection-dialog") == signature:
        _selection_dialog(database, user, projects, selected_entries, project_id)
