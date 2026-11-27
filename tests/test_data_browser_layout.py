from copy import deepcopy
import json
from types import SimpleNamespace

import streamlit as st
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import data_browser_ui
from test_data_browser import inventory
from test_database_browser import USER, database
from test_database_design_evaluations import design
from test_database_outcomes import _pass_design


def _layout(app, value):
    next(element for element in app.get("button_group") if element.key == "data-browser-layout").set_value(value).run()
    assert not app.exception


def test_workspace_can_focus_either_pane_or_unpin_only_folders():
    app = AppTest.from_string("from mn_molecular_design_studio.data_browser_ui import _workspace\n_workspace()").run()
    assert not app.exception
    assert app.session_state["data-browser-layout"] == "Split"
    assert app.session_state["data-browser-pin-folders"]
    assert len(app.columns) == 6
    app.slider(key="data-browser-split-details-width").set_value(46).run()
    assert [round(column.proto.weight * 100) for column in app.columns[-3:]] == [27, 27, 46]
    app.toggle[0].set_value(False).run()
    assert len(app.columns) == 5
    assert len(app.get("popover")) == 2
    for mode in ("Table", "Structure"):
        _layout(app, mode)
        assert not app.session_state["data-browser-pin-folders"]
        assert len(app.columns) == 3
        assert not app.slider
        app.toggle[0].set_value(True).run()
        assert len(app.columns) == 5
        assert app.slider(key=f"data-browser-{mode.lower()}-pinned-folders-width").value == 27
    _layout(app, "Split")
    assert app.session_state["data-browser-pin-folders"]
    assert len(app.columns) == 6
    assert app.slider(key="data-browser-split-details-width").value == 46
    app.button(key="data-browser-split-reset-widths").click().run()
    assert app.slider(key="data-browser-split-details-width").value == 30
    _layout(app, None)
    assert app.session_state["data-browser-layout"] == "Split"
    assert app.session_state["data-browser-pin-folders"]


def test_layout_roundtrip_retains_folder_filters_focus_and_batch_without_sql_changes(database, design, monkeypatch):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    outcome = next(row for row in outcomes if row["outcome_kind"] == "Design")
    second = {**deepcopy(outcome), "id": "second-design", "name": "Second design"}
    monkeypatch.setattr("mn_molecular_design_studio.data_browser.saved_results", lambda *args, **kwargs: [outcome, second])
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    details = []
    monkeypatch.setattr(data_browser_ui, "_render_details", lambda database, user, cache, index, entry: details.append(entry["id"]))
    dataframe = st.dataframe
    defaults = []

    def selected_dataframe(*args, **kwargs):
        rendered = dataframe(*args, **kwargs)
        if str(kwargs.get("key") or "").startswith("data-browser-table-"):
            defaults.append(kwargs["selection_default"]["selection"]["rows"])
            injected = st.session_state.get("test-selected-rows")
            if injected is not None:
                return SimpleNamespace(selection=SimpleNamespace(rows=injected))
        return rendered

    monkeypatch.setattr(st, "dataframe", selected_dataframe)
    with database.connection() as connection:
        original = [tuple(row) for row in connection.execute("SELECT * FROM molecule_observations ORDER BY id")]
    app = AppTest.from_string("from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    next(element for element in app.text_input if element.label == "Search this folder").set_value("run-1").run()
    next(element for element in app.selectbox if element.key == "data-browser-type").set_value("protein").run()
    app.session_state["test-selected-rows"] = [1]
    app.run()
    assert not app.exception
    scope = app.session_state["data-browser-scope"]
    retained = deepcopy(app.session_state["data-browser-page-selection"])
    assert retained["selected"] == ["second-design"]
    assert retained["focused"] == "second-design"
    del app.session_state["test-selected-rows"]
    _layout(app, "Structure")
    assert app.session_state["data-browser-page-selection"] == retained
    assert next(element for element in app.selectbox if element.label == "Result on this page").value == "second-design"
    assert details[-1] == "second-design"
    count = len(details)
    _layout(app, "Table")
    assert len(details) == count
    assert defaults[-1] == [1]
    assert app.session_state["data-browser-page-selection"] == retained
    _layout(app, "Split")
    assert app.session_state["data-browser-page-selection"] == retained
    assert app.session_state["data-browser-scope"] == scope
    assert app.session_state["data-browser-result-query"] == "run-1"
    assert app.session_state["data-browser-type"] == "protein"
    assert details[-1] == "second-design"
    with database.connection() as connection:
        assert original == [tuple(row) for row in connection.execute("SELECT * FROM molecule_observations ORDER BY id")]


def test_stale_selected_ids_are_pruned_without_selecting_another_result():
    app = AppTest.from_string("import streamlit as st\n"
        "from mn_molecular_design_studio.data_browser_ui import _page_selection\n"
        "st.session_state.setdefault('data-browser-page-selection', {'context':('folder',1),'selected':['removed'],'focused':'removed'})\n"
        "_page_selection('folder',1,[{'id':'available'}])").run()
    assert not app.exception
    assert app.session_state["data-browser-page-selection"]["selected"] == []
    assert app.session_state["data-browser-page-selection"]["focused"] == "available"


def test_native_table_checked_subset_survives_analysis_button_rerun(database, design, monkeypatch):
    from mn_molecular_design_studio import binderscope_ui

    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    outcome = next(row for row in outcomes if row["outcome_kind"] == "Design")
    second = {**deepcopy(outcome), "id": "second-design", "name": "Second design"}
    monkeypatch.setattr("mn_molecular_design_studio.data_browser.saved_results", lambda *args, **kwargs: [outcome, second])
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    monkeypatch.setattr(data_browser_ui, "_render_details", lambda *args: None)
    opened = []
    monkeypatch.setattr(binderscope_ui, "_browser_dialog",
        lambda database, user, index, chosen, token, **kwargs: opened.append([entry["id"] for entry in chosen]))
    app = AppTest.from_string("import streamlit as st\n"
        "st.session_state.setdefault('data-browser-layout','Table')\n"
        "from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    table = next(frame for frame in app.dataframe if frame.proto.id)
    states = app._tree.get_widget_states()
    state = next((state for state in states.widgets if state.id == table.proto.id), None)
    if state is None:
        state = states.widgets.add(id=table.proto.id)
    state.string_value = json.dumps({"selection": {"rows": [1], "columns": [], "cells": []}})
    app._run(states)
    assert not app.exception
    retained = app.session_state["data-browser-page-selection"]["selected"]
    assert retained == ["second-design"]
    next(button for button in app.button if button.label == "Analyze · BinderScope").click().run()
    assert not app.exception and opened == [retained]
    assert app.session_state["data-browser-page-selection"]["selected"] == retained


def test_empty_folder_has_no_stale_focus():
    app = AppTest.from_string("from mn_molecular_design_studio.data_browser_ui import _page_selection\n"
        "_page_selection('empty',1,[])").run()
    assert not app.exception
    assert app.session_state["data-browser-page-selection"]["focused"] == ""


def test_whole_view_displays_506_rows_and_paging_can_be_restored():
    source = ("import streamlit as st\nfrom mn_molecular_design_studio.data_browser_ui import _listing_window\n"
        "visible,page,signature=_listing_window([{'id':str(number)} for number in range(506)],'folder')\n"
        "st.session_state['test-visible']=[entry['id'] for entry in visible]\n"
        "st.session_state['test-view-signature']=signature")
    app = AppTest.from_string(source).run()
    assert not app.exception and len(app.session_state["test-visible"]) == 506
    assert not app.number_input
    whole_signature = app.session_state["test-view-signature"]
    app.selectbox(key="data-browser-rows-shown").set_value("50").run()
    assert not app.exception and len(app.session_state["test-visible"]) == 50
    assert app.number_input[0].max == 11
    assert app.session_state["test-view-signature"] != whole_signature
    app.number_input[0].set_value(3).run()
    assert app.session_state["test-visible"] == [str(number) for number in range(100,150)]
    app.selectbox(key="data-browser-rows-shown").set_value("All").run()
    assert not app.exception and len(app.session_state["test-visible"]) == 506
    assert not app.number_input
    assert app.session_state["test-view-signature"] == whole_signature


def test_whole_view_is_limited_to_1000_rows_and_large_folders_restore_paging():
    app = AppTest.from_string("import streamlit as st\n"
        "from mn_molecular_design_studio.data_browser_ui import _listing_window\n"
        "count=st.number_input('Test total',value=1000,step=1)\n"
        "visible,page,signature=_listing_window([{'id':str(number)} for number in range(count)],'folder')\n"
        "st.session_state['test-visible']=len(visible)").run()
    assert not app.exception and app.session_state["test-visible"] == 1000
    app.number_input[0].set_value(1001).run()
    assert not app.exception and app.session_state["test-visible"] == 50
    assert app.selectbox(key="data-browser-rows-shown").value == "50"
    assert "All" not in app.selectbox(key="data-browser-rows-shown").options


def test_structure_choice_survives_hiding_details_and_focus_increases_viewer_height(tmp_path, monkeypatch):
    from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer

    path = tmp_path / "complex.pdb"
    path.write_text("HEADER    test structure\nEND\n")
    monkeypatch.setattr(ReadOnlyLibrary, "structure_file", lambda *args: path)
    heights = []
    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda structures, **kwargs: heights.append(kwargs["height"]))
    record = {"id": "saved-design", "name": "Saved design", "type": "protein", "selected_evaluation": {"id": "evaluation"},
        "structures": [{"id": "refolded", "label": "Refolded validation complex"}, {"id": "generated", "label": "Generated design complex"}]}
    app = AppTest.from_string("import streamlit as st\n"
        "from pathlib import Path\nfrom mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.data_browser_ui import _render_structure\n"
        "mode=st.radio('Test workspace',['Split','Table','Structure'])\n"
        "st.session_state['data-browser-layout']=mode\n"
        f"if mode!='Table': _render_structure(ReadOnlyLibrary(Path({str(tmp_path)!r})),{USER!r},{record!r},'')").run()
    assert not app.exception and heights[-1] == 320
    app.selectbox[0].set_value("generated").run()
    app.radio[0].set_value("Table").run()
    assert not app.selectbox
    app.radio[0].set_value("Structure").run()
    assert not app.exception
    assert app.selectbox[0].value == "generated"
    assert heights[-1] == 720
