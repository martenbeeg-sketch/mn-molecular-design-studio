import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import data_browser_ui, navigation
from mn_molecular_design_studio.database_browser import build_browser_index
from mn_molecular_design_studio.data_browser import source_collection_scope
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from test_database_browser import OTHER, USER, database


def _app(database, monkeypatch, *, state=None, projects=(), inputs=False):
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    monkeypatch.setattr(data_browser_ui, "_render_details", lambda *args, **kwargs: None)
    renderer = "datasets_ui import render_datasets" if inputs else "data_browser_ui import render_data_browser"
    function = "render_datasets" if inputs else "render_data_browser"
    return AppTest.from_string("import streamlit as st\n"
        f"for key,value in {(state or {})!r}.items(): st.session_state.setdefault(key,value)\n"
        f"from mn_molecular_design_studio.{renderer}\n"
        f"{function}({USER!r}, {list(projects)!r})", default_timeout=30).run()


def _library(database, *, user=USER, name="Research compounds", **kwargs):
    return database.import_records(user, [
        {"name": "Ethanol", "components": [{"chain": "L", "type": "small_molecule", "smiles": "CCO"}]},
        {"name": "Methanol", "components": [{"chain": "L", "type": "small_molecule", "smiles": "CO"}]},
    ], name=name, **kwargs)


def _select_rows(app, marker, rows):
    table = next(frame for frame in app.dataframe if marker in frame.proto.id)
    states = app._tree.get_widget_states()
    state = next((state for state in states.widgets if state.id == table.proto.id), None)
    if state is None:
        state = states.widgets.add(id=table.proto.id)
    state.string_value = json.dumps({"selection": {"rows": rows, "columns": [], "cells": []}})
    app._run(states)
    assert not app.exception


def test_single_database_navigation_accepts_old_browser_sessions():
    assert [page for page in navigation.WORKSPACE_PAGES if page[2] in {"Database", "Data Browser"}] == [
        ("Data browser", ":material/folder_copy:", "Data Browser")]
    assert any(page[2] == "Datasets" for page in navigation.WORKSPACE_PAGES)
    app = AppTest.from_string("import streamlit as st\n"
        "st.session_state.setdefault('page','Database')\n"
        "from mn_molecular_design_studio.navigation import render_workspace_navigation\n"
        f"render_workspace_navigation({USER!r})").run()
    assert not app.exception
    assert app.session_state["page"] == "Data Browser"
    assert app.button(key="nav-data-browser").proto.type == "primary"
    assert not any(button.label == "Database" for button in app.button)


def test_sidebar_omits_separate_complex_analysis_without_removing_legacy_handoffs():
    app = AppTest.from_string("import streamlit as st\n"
        "st.session_state.setdefault('page','Complex Analysis')\n"
        "from mn_molecular_design_studio.navigation import render_workspace_navigation\n"
        f"render_workspace_navigation({USER!r})").run()
    assert not app.exception
    assert not any(button.label == "Complex Analysis" for button in app.button)
    assert app.button(key="nav-all-tools") and app.button(key="nav-data-browser")
    assert app.session_state["page"] == "Complex Analysis"
    source = (Path(data_browser_ui.__file__).parent / "run_app.py").read_text()
    assert 'elif page == "Complex Analysis":' in source


def test_both_database_routes_dispatch_only_the_new_browser():
    source = (Path(data_browser_ui.__file__).parent / "run_app.py").read_text()
    assert 'elif page in {"Database", "Data Browser"}:\n    _render_database()' in source
    assert "database_ui.render_database(" not in source


def test_project_handoff_resets_new_browser_and_keeps_cached_metadata():
    cache = object()
    state = {"data-browser-cache": cache, "data-browser-source-group": "wrong",
        "data-browser-scope-path": ["collections", "wrong"], "data-browser-page-selection": {"selected": ["wrong"]},
        "data-browser-type": "protein", "data-browser-result-query": "old", "data-browser-indexed-history": True}
    navigation.open_project_database(state, "research")
    assert navigation.database_destination(state) == {"project_id": "research"}
    assert state["data-browser-cache"] is cache
    assert state["data-browser-result-query"] == "" and state["data-browser-type"] == ""
    assert not state["data-browser-indexed-history"]
    assert "data-browser-source-group" not in state and "data-browser-page-selection" not in state


@pytest.mark.parametrize("kind", ["curated", "generated", "inputs"])
def test_explicit_source_links_open_exact_authorized_collection(database, monkeypatch, kind):
    saved = _library(database, import_method="tool_inputs" if kind == "inputs" else "user_uploads",
        origin={"collection_kind": kind})
    app = _app(database, monkeypatch, state={"molecular-preselect-group": saved["group_id"],
        "molecular-preselect-molecule": saved["molecule_ids"][-1], "molecular-preselect-collection": "all"})
    assert not app.exception
    assert app.session_state["data-browser-entry-v2"] == "collections"
    assert app.session_state["data-browser-scope"] == "collection:" + saved["group_id"]
    assert app.session_state["data-browser-indexed-history"]
    assert "molecular-preselect-group" not in app.session_state
    assert not app.session_state["data-browser-page-selection"]["selected"]
    assert not any("existing Database page" in item.value for item in app.info)


def test_inaccessible_source_link_never_opens_a_different_collection(database, monkeypatch):
    _library(database)
    private = _library(database, user=OTHER, name="Secret collection")
    before = database.list_molecules(USER)["total"]
    app = _app(database, monkeypatch, state={"molecular-preselect-group": private["group_id"]})
    assert not app.exception
    assert app.session_state["data-browser-scope"] == "collections"
    assert any("No other collection has been substituted" in item.value for item in app.warning)
    assert private["group_id"] not in app.session_state["data-browser-scopes"]
    assert database.list_molecules(USER)["total"] == before
    assert source_collection_scope(build_browser_index(database, USER), private["group_id"]) is None


def test_project_link_opens_project_in_new_tree(database, monkeypatch):
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["research"])
    _library(database, project_id="research")
    app = _app(database, monkeypatch, state={"data-browser-open": {"project_id": "research"}},
        projects=[{"id": "research", "name": "Research"}])
    assert not app.exception
    assert app.session_state["data-browser-scope"] == "project:research"
    assert app.session_state["data-browser-entry-v2"] == "projects"


def test_empty_database_keeps_import_and_creation_accessible(database, monkeypatch):
    calls = []
    monkeypatch.setattr("mn_molecular_design_studio.database_ui.render_import_dialog",
        lambda user, projects, **kwargs: calls.append(kwargs["database"].root))
    app = _app(database, monkeypatch, inputs=True)
    assert not app.exception
    assert not app.button(key="data-browser-import").disabled
    app.button(key="data-browser-import").click().run()
    assert not app.exception and calls == [database.root]
    assert not database.path.exists()


def test_import_into_configured_database_opens_the_new_collection(database, monkeypatch):
    app = _app(database, monkeypatch, inputs=True)
    app.button(key="data-browser-import").click().run()
    app.text_input(key="molecular-import-group-name").set_value("New research library").run()
    app.button(key="molecular-import-next").click().run()
    app.radio(key="molecular-import-method").set_value("SMILES").run()
    assert not app.exception
    assert app.text_area
    app.text_area[0].set_value("CCO ethanol").run()
    next(item for item in app.checkbox if item.label.startswith("I reviewed the report")).set_value(True).run()
    app.button(key="molecular-import-commit").click().run()
    assert not app.exception
    assert database.groups(USER)[0]["name"] == "New research library"
    assert app.session_state["datasets-selected"] == database.groups(USER)[0]["id"]


def test_new_empty_collection_is_not_lost_and_metadata_columns_are_accessible(database, monkeypatch):
    app = _app(database, monkeypatch, inputs=True)
    app.button(key="data-browser-new-collection").click().run()
    app.text_input(key="molecular-group-name").set_value("New empty dataset").run()
    app.button(key="molecular-group-save").click().run()
    assert not app.exception
    group = database.groups(USER)[0]
    assert group["molecule_count"] == 0
    assert app.session_state["datasets-selected"] == group["id"]
    app.button(key="data-browser-manage-collections").click().run()
    app.button(key="data-browser-manage-column").click().run()
    app.text_input(key="molecular-column-name").set_value("Research notes").run()
    app.button(key="molecular-column-save").click().run()
    assert not app.exception
    assert database.group(group["id"], USER)["schema"]["fields"][0]["name"] == "Research notes"


def test_read_only_preview_does_not_open_write_dialogs(database, monkeypatch):
    _library(database)
    monkeypatch.setenv("MDS_READ_ONLY_PREVIEW", "1")
    app = _app(database, monkeypatch, inputs=True, state={"molecular-import-open": True,
        "data-browser-management": {"kind": "group", "ids": []}})
    assert not app.exception
    assert app.button(key="data-browser-import").disabled
    assert app.button(key="data-browser-new-collection").disabled
    assert not app.get("dialog")


def test_focus_handoff_uses_correct_page_without_selecting_for_analysis():
    source = ("import streamlit as st\nfrom mn_molecular_design_studio.data_browser_ui import _listing_window,_page_selection\n"
        "st.session_state['data-browser-rows-shown']='50'\n"
        "entries=[{'id':str(number),'result':{'molecule_id':str(number)}} for number in range(506)]\n"
        "visible,page,signature=_listing_window(entries,'folder',focus_molecule='450')\n"
        "st.session_state['test-page']=page\n_page_selection(signature,page,visible)")
    app = AppTest.from_string(source).run()
    assert not app.exception
    assert app.session_state["test-page"] == 10
    assert app.session_state["data-browser-page-selection"]["selected"] == []


def test_collection_manager_only_removes_confirmed_group_members(database, monkeypatch):
    saved = _library(database)
    duplicate = _library(database, name="Retained elsewhere")
    app = _app(database, monkeypatch)
    app.button(key="data-browser-manage-collections").click().run()
    app.selectbox(key="data-browser-manage-group").set_value(saved["group_id"]).run()
    _select_rows(app, "data-browser-manage-rows-", [0])
    confirm = next(item for item in app.checkbox if "Remove the selected entries" in item.label)
    confirm.set_value(True).run()
    assert not app.button(key="data-browser-manage-remove").disabled
    app.button(key="data-browser-manage-remove").click().run()
    assert not app.exception
    counts = {group["id"]: group["molecule_count"] for group in database.groups(USER)}
    assert counts[saved["group_id"]] == 1
    assert counts[duplicate["group_id"]] == 2


def test_collection_removal_confirmation_is_bound_to_the_exact_selection(database, monkeypatch):
    _library(database)
    app = _app(database, monkeypatch)
    app.button(key="data-browser-manage-collections").click().run()
    _select_rows(app, "data-browser-manage-rows-", [0])
    next(item for item in app.checkbox if "Remove the selected entries" in item.label).set_value(True).run()
    assert not app.button(key="data-browser-manage-remove").disabled
    _select_rows(app, "data-browser-manage-rows-", [1])
    assert app.button(key="data-browser-manage-remove").disabled
    assert database.groups(USER)[0]["molecule_count"] == 2


def test_selected_molecules_keep_export_and_tool_reuse(database, monkeypatch):
    saved = _library(database)
    entries = [{"result": {"molecule_id": molecule_id}} for molecule_id in saved["molecule_ids"]]
    app = AppTest.from_string("from pathlib import Path\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.database_actions_ui import render_selection_actions\n"
        f"render_selection_actions(ReadOnlyLibrary(Path({str(database.root)!r})), {USER!r}, [], {entries!r})").run()
    app.button(key="data-browser-selection-actions").click().run()
    assert not app.exception
    assert app.selectbox(key="database-use-tool").options[:3] == ["AutoDock Vina", "GNINA", "Uni-Dock Pro"]
    app.button(key="data-browser-molecule-export-prepare").click().run()
    assert not app.exception
    assert any(item.label == "Download CSV" for item in app.get("download_button"))
    assert database.groups(USER)[0]["molecule_count"] == 2


def test_selection_action_never_implicitly_selects_the_current_folder(database):
    app = AppTest.from_string("from pathlib import Path\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.database_actions_ui import render_selection_actions\n"
        f"render_selection_actions(ReadOnlyLibrary(Path({str(database.root)!r})), {USER!r}, [], [])").run()
    assert not app.exception
    assert app.button(key="data-browser-selection-actions").disabled
