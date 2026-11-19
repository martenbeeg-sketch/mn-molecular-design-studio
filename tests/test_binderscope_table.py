from copy import deepcopy
from contextlib import contextmanager
import json

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.binderscope_table import (
    DEFAULT_FIELDS, OPTIONAL_FIELDS, PREFIX, VIEWER_COLUMN, inspector_choices, inspector_label, summary_columns,
)
from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
from mn_molecular_design_studio.data_browser_table import column_label, default_columns, result_table, table_csv
from test_data_browser_analyses import listing
from test_database_browser import OTHER, USER, database
from test_database_design_evaluations import _analysis, design


def report(design):
    return {"source": {"input_structure_digest": design["digest"],
        "source": {"structure_label": "Refolded validation complex"},
        "assessment": {"b_factor_interpretation": "plddt_0_100"}},
        "analysis": {"n": 3, "plddt_res": [60., 80., 100.]},
        "hotspot_compliance": {"status": "assessed", "requested_count": 3, "contacted_count": 0, "unmapped_count": 0},
        "target_reference_comparison": {"status": "assessed", "target_ca_rmsd_angstrom": 0.},
        "termini": {"applies": True, "n": {"conf": "low", "plddt": 68.5,
            "g": {key: {"s": "free", "l": 0} for key in ("his6", "strep", "domain", "large")}},
            "c": {"conf": "ok", "plddt": 77.6, "g": {"his6": {"s": "free"}, "large": {"s": "linker", "l": 5}}}}}


def add_saved(database, design, *, run="binder-review", saved=None, context=None):
    return _analysis(database, design, "binderscope", run=run,
        metrics={"binderscope_buried_area_angstrom2": 1465.22693, "binderscope_hydrogen_bond_count": 7,
            "binderscope_salt_bridge_count": 3, "binderscope_buried_unpaired_polar_atom_count": 0},
        context={"candidate": {"source_tool": "binderscope", "stage": "analysis",
            "raw_metadata": {"engine_version": "review-v2", "binderscope": saved if saved is not None else report(design)}},
            **(context or {})})


def test_saved_summary_defaults_optional_details_and_csv_are_source_bound_and_read_only(database, design):
    identifier = add_saved(database, design)
    index, entries = listing(database, design)
    before = database.path.read_bytes()
    row = result_table(ReadOnlyLibrary(database.root), index, entries)[0]
    assert row[PREFIX + "Status"] == "Recorded"
    assert row[PREFIX + "Analysis version"] == "review-v2"
    assert row[PREFIX + "Hotspots"] == "0 / 3"
    assert row[PREFIX + "Mean binder pLDDT"] == 80.
    assert row[PREFIX + "N terminal pLDDT"] == 68.5
    assert row[PREFIX + "N terminal confidence"] == "Low"
    assert row[PREFIX + "C terminal confidence"] == "OK"
    assert row[PREFIX + "N tag accessibility"] == "Free · all 4 families"
    assert row[PREFIX + "C tag accessibility"] == "His6: free; Fc/MBP: linker (5 aa)"
    assert row[PREFIX + "Target-reference RMSD (Å)"] == 0.
    assert row[PREFIX + "Result · binderscope_buried_unpaired_polar_atom_count"] == 0
    assert row["Result · confidence"] == .9
    defaults = default_columns([row])
    assert all(PREFIX + field in defaults for field in DEFAULT_FIELDS)
    assert not any(PREFIX + field in defaults for field in OPTIONAL_FIELDS)
    detailed = default_columns([row], binderscope_details=True)
    assert all(PREFIX + field in detailed for field in OPTIONAL_FIELDS)
    assert column_label(PREFIX + "Result · binderscope_buried_area_angstrom2") == "BinderScope: Buried interface area (Å²)"
    assert column_label(PREFIX + "Result · binderscope_hydrogen_bond_count") == "BinderScope: H-bonds"
    assert inspector_choices(row)[0]["observation_id"] == identifier
    exported = table_csv([row]).decode("utf-8-sig")
    assert "1465.22693" in exported and "N terminal confidence" in exported
    assert database.path.read_bytes() == before


@pytest.mark.parametrize("context,status", [
    ({"status": "failed"}, "Failed / unavailable"),
    ({"status": "pending"}, "Status unknown / incomplete"),
    ({"source_structure_id": "other", "representation_id": "other", "source_structure_digest": "b" * 64,
        "structure_digest": "b" * 64}, "Other model"),
    ({"target_molecule_id": "other-target"}, "Not assessed"),
])
def test_unqualified_reports_never_supply_inspector_or_summary_cells(database, design, context, status):
    add_saved(database, design, context=context)
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row[PREFIX + "Status"] == status
    assert not inspector_choices(row)
    assert PREFIX + "Hotspots" not in row and PREFIX + "Mean binder pLDDT" not in row
    assert PREFIX + "Result · binderscope_buried_area_angstrom2" not in row


def test_conflicting_report_checksum_is_not_shown_as_a_scored_assessment(database, design):
    saved = report(design)
    saved["source"]["input_structure_digest"] = "f" * 64
    add_saved(database, design, saved=saved)
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row[PREFIX + "Status"] == "Report checksum mismatch"
    assert not inspector_choices(row)
    assert not any(field.startswith(PREFIX + "Result · ") for field in row)


@pytest.mark.parametrize("context,expected", [
    ({}, 1), ({"status": "failed"}, 0), ({"status": "pending"}, 0),
    ({"representation_id": "other-model"}, 0),
    ({"source_structure_digest": "f" * 64}, 0),
    ({"calculation_structure_digest": "f" * 64}, 0),
    ({"binder_chains": ["A"], "target_chains": ["B"]}, 0),
    ({"target_molecule_id": "other-target"}, 0),
])
def test_preflight_skips_only_verified_completed_same_complex_assessments(database, design, context, expected):
    from mn_molecular_design_studio.binderscope import existing_assessments

    add_saved(database, design, context={"binder_chains": ["B"], "target_chains": ["A"], **context})
    source = {"molecule_id": design["molecule_id"], "structure_id": design["structure"],
        "group_id": design["source_group_id"], "structure_digest": design["digest"],
        "binder_chains": ["B"], "target_chains": ["A"],
        "chain_role_evidence": {"target_molecule_id": design["target_id"]}}
    key = (source["molecule_id"], source["structure_id"], source["group_id"])
    before = database.path.read_bytes()
    assert existing_assessments(database, USER, [source]) == {key: expected}
    assert existing_assessments(database, OTHER, [source]) == {}
    assert database.path.read_bytes() == before


def test_preflight_does_not_skip_unverified_reports_or_another_project(database, design):
    from mn_molecular_design_studio.binderscope import existing_assessments

    source = {"molecule_id": design["molecule_id"], "structure_id": design["structure"],
        "group_id": design["source_group_id"], "structure_digest": design["digest"],
        "binder_chains": ["B"], "target_chains": ["A"]}
    key = (source["molecule_id"], source["structure_id"], source["group_id"])
    roles = {"binder_chains": ["B"], "target_chains": ["A"]}
    _analysis(database, design, "binderscope", run="legacy", context=roles)
    saved = report(design)
    saved["source"]["input_structure_digest"] = "f" * 64
    add_saved(database, design, run="corrupt", saved=saved, context=roles)
    assert existing_assessments(database, USER, [source]) == {key: 0}
    add_saved(database, design, context=roles)
    assert existing_assessments(database, USER, [source]) == {key: 1}
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET project_id='another-project' WHERE id=?", (design["group_id"],))
    member = {**USER, "_cluster_project_ids": ["another-project"]}
    assert existing_assessments(database, member, [source]) == {key: 0}


def test_legacy_metrics_without_report_remain_recorded_but_cannot_open_a_guessed_report(database, design):
    _analysis(database, design, "binderscope", metrics={"binderscope_hydrogen_bond_count": 0})
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row[PREFIX + "Status"] == "Recorded"
    assert row[PREFIX + "Result · binderscope_hydrogen_bond_count"] == 0
    assert not inspector_choices(row)


def test_recorded_analysis_method_version_is_available_without_an_engine_version(database, design):
    identifier = add_saved(database, design)
    with database.connection() as connection:
        connection.execute("UPDATE molecule_observations SET context=json_set(json_remove(context,"
            "'$.candidate.raw_metadata.engine_version'),'$.candidate.raw_metadata.analysis_parameters.method_version',?) WHERE id=?",
            ("studio.binderscope.saved-cohort.v4", identifier))
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row[PREFIX + "Analysis version"] == "studio.binderscope.saved-cohort.v4"
    assert inspector_choices(row)[0]["analysis_version"] == row[PREFIX + "Analysis version"]


def test_multiple_assessments_stay_explicit_and_never_fill_numeric_cells(database, design):
    first = add_saved(database, design, run="first")
    second = add_saved(database, design, run="second")
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row[PREFIX + "Status"] == "Multiple evaluations (2)"
    assert {choice["observation_id"] for choice in inspector_choices(row)} == {first, second}
    assert not any(field.startswith(PREFIX + "Result · ") for field in row)
    assert PREFIX + "Mean binder pLDDT" not in row


def test_summary_does_not_guess_confidence_missing_hotspots_or_free_cyclic_termini():
    summary = {"encoding": "uninterpreted", "plddt": [90, 95], "residues": 2,
        "hotspot_status": "unavailable", "n": {"conf": "off", "plddt": 95},
        "c": {}, "termini_status": "unavailable", "reference_status": "unavailable", "rmsd": 0}
    row = summary_columns(summary)
    assert row[PREFIX + "Hotspots"] == "Not assessed"
    assert row[PREFIX + "N tag accessibility"] == "Unavailable"
    assert PREFIX + "Mean binder pLDDT" not in row and PREFIX + "N terminal pLDDT" not in row
    assert PREFIX + "Target-reference RMSD (Å)" not in row
    row = summary_columns({**summary, "termini_applies": False, "hotspot_status": "not_specified"})
    assert row[PREFIX + "N tag accessibility"] == "Not applicable"
    assert row[PREFIX + "C tag accessibility"] == "Not applicable"
    assert row[PREFIX + "Hotspots"] == "Not specified"
    row = summary_columns({**summary, "termini_status": None,
        "n": {"conf": "ok", "plddt": 95, "g": {"his6": {"s": "free"}}}})
    assert row[PREFIX + "N tag accessibility"] == "Unavailable"
    assert PREFIX + "N terminal confidence" not in row


@pytest.mark.parametrize("values", [[True, 90], [float("nan"), 90], [-1, 101], [90], "90,95"])
def test_invalid_or_incomplete_plddt_is_not_averaged(values):
    assert PREFIX + "Mean binder pLDDT" not in summary_columns({"encoding": "plddt_0_100", "plddt": values, "residues": 2})


@pytest.mark.parametrize("records", ["bad json", "{}", '[null,42,"bad"]'])
def test_bad_report_records_offer_no_inspector(records):
    assert inspector_choices({PREFIX + "Records": records}) == []


def browser_action(database, design, monkeypatch, *, user=USER):
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    opened = []
    monkeypatch.setattr("mn_molecular_design_studio.binderscope_inspector.inspection_dialog",
        lambda result, database, user, key: opened.append(result))
    source = ("from pathlib import Path\nfrom mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.binderscope_ui import render_saved_browser_action\n"
        f"render_saved_browser_action(ReadOnlyLibrary(Path({str(database.root)!r})),{user!r},{index!r},{entries[0]!r},{row!r},signature='folder')")
    app = AppTest.from_string(source).run()
    assert not app.exception
    return app, opened


def test_table_shortcut_opens_saved_report_without_launching_an_analysis(database, design, monkeypatch):
    identifier = add_saved(database, design)
    app, opened = browser_action(database, design, monkeypatch)
    before = database.path.read_bytes()
    assert not opened and not app.button[0].disabled
    assert not any("Queue" in button.label or "Analyze" in button.label for button in app.button)
    app.button[0].click().run()
    assert not app.exception and len(opened) == 1
    assert opened[0]["_studio_sql_binding"]["observation_id"] == identifier
    assert database.path.read_bytes() == before


def test_table_shortcut_requires_an_explicit_choice_for_multiple_assessments(database, design, monkeypatch):
    add_saved(database, design, run="first")
    second = add_saved(database, design, run="second")
    app, opened = browser_action(database, design, monkeypatch)
    assert app.selectbox[0].value is None and app.button[0].disabled and not opened
    app.selectbox[0].set_value(second).run()
    app.button[0].click().run()
    assert not app.exception
    assert opened[0]["_studio_sql_binding"]["observation_id"] == second


@pytest.mark.parametrize("change", ["delete", "failed", "private"])
def test_open_rechecks_current_status_deletion_and_permissions(database, design, monkeypatch, change):
    identifier = add_saved(database, design)
    app, opened = browser_action(database, design, monkeypatch, user=OTHER if change == "private" else USER)
    with database.connection() as connection:
        if change == "delete":
            connection.execute("DELETE FROM molecule_observations WHERE id=?", (identifier,))
        elif change == "failed":
            connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.status','failed') WHERE id=?", (identifier,))
    app.button[0].click().run()
    assert not app.exception and not opened
    assert "no longer accessible" in app.warning[0].value


def test_optional_details_toggle_changes_columns_without_changing_data(database, design):
    add_saved(database, design)
    index, entries = listing(database, design)
    source = ("from pathlib import Path\nfrom mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.data_browser_ui import _table_controls\nimport streamlit as st\n"
        f"rows,columns=_table_controls(ReadOnlyLibrary(Path({str(database.root)!r})),{index!r},{entries!r},{entries!r},'snapshot','folder',1)\n"
        "st.session_state['test-columns']=columns\nst.session_state['test-rows']=rows")
    app = AppTest.from_string(source).run()
    assert not app.exception
    before = deepcopy(app.session_state["test-rows"])
    assert PREFIX + "N terminal pLDDT" not in app.session_state["test-columns"]
    app.toggle(key="data-browser-binderscope-details").set_value(True).run()
    assert not app.exception
    assert PREFIX + "N terminal pLDDT" in app.session_state["test-columns"]
    assert app.session_state["test-rows"] == before
    app.toggle(key="data-browser-binderscope-details").set_value(False).run()
    assert PREFIX + "N terminal pLDDT" not in app.session_state["test-columns"]


def test_page_summary_is_one_batched_sql_read_and_never_loads_coordinates(database, design, monkeypatch):
    add_saved(database, design)
    index, entries = listing(database, design)
    library = ReadOnlyLibrary(database.root)
    queries = []
    original = library.connection

    @contextmanager
    def traced():
        with original() as connection:
            connection.set_trace_callback(queries.append)
            yield connection

    def unexpected_coordinates(*args, **kwargs):
        pytest.fail("Table summaries must not read coordinates")

    monkeypatch.setattr(library, "connection", traced)
    monkeypatch.setattr(library, "structure_file", unexpected_coordinates)
    rows = result_table(library, index, entries * 50)
    assert len(rows) == 50 and all(row[PREFIX + "Hotspots"] == "0 / 3" for row in rows)
    assert sum(" AS summary," in query for query in queries) == 1


def click_viewer_row(app, position=0, label="View BinderScope"):
    table = next(frame for frame in app.dataframe if VIEWER_COLUMN in frame.value.columns)
    states = app._tree.get_widget_states()
    trigger = states.widgets.add(id=table.proto.button_click_widgets[VIEWER_COLUMN])
    trigger.string_trigger_value.data = json.dumps({"row": position, "label": label})
    return app._run(states)


def test_saved_row_button_works_in_table_workspace_without_rendering_details_or_changing_selection(database, design, monkeypatch):
    from mn_molecular_design_studio import data_browser_ui

    add_saved(database, design)
    listing(database, design)
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    details, opened = [], []
    monkeypatch.setattr(data_browser_ui, "_render_details", lambda *args: details.append(True))
    monkeypatch.setattr("mn_molecular_design_studio.binderscope_inspector.render_inspector",
        lambda result, **kwargs: opened.append(result))
    source = ("import streamlit as st\nfrom mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"st.session_state.setdefault('data-browser-layout','Table')\nrender_data_browser({USER!r}, [])")
    app = AppTest.from_string(source, default_timeout=30).run()
    assert not app.exception and not details
    retained = deepcopy(app.session_state["data-browser-page-selection"])
    assert not any(button.label == "Open saved BinderScope · 3D & interface" for button in app.button)
    table = next(frame for frame in app.dataframe if VIEWER_COLUMN in frame.value.columns)
    assert table.value[VIEWER_COLUMN].tolist() == ["View BinderScope"]
    columns = list(table.proto.column_order)
    assert columns.index(VIEWER_COLUMN) == columns.index("Name") + 1
    assert PREFIX + "Records" not in table.value.columns
    assert PREFIX + "Records" in app.session_state["data-browser-page-table"][1][0]
    launch = next(button for button in app.button if button.label.startswith("Analyze · BinderScope"))
    assert not launch.disabled and launch.label == "Analyze · BinderScope"
    assert all(not any(button.key == launch.key for button in popover.get("button")) for popover in app.get("popover"))
    click_viewer_row(app)
    assert not app.exception and len(opened) == 1 and not details
    assert app.session_state["data-browser-layout"] == "Table"
    assert app.session_state["data-browser-page-selection"] == retained
    assert "binderscope-row-request" not in app.session_state
    app.run()
    assert len(opened) == 1


def test_row_button_labels_distinguish_multiple_and_missing_assessments(database, design):
    assert inspector_label({}) is None
    add_saved(database, design, run="one")
    index, entries = listing(database, design)
    assert inspector_label(result_table(database, index, entries)[0]) == "View BinderScope"
    add_saved(database, design, run="two")
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert inspector_label(row) == "Choose assessment (2)"
    assert VIEWER_COLUMN not in table_csv([row]).decode("utf-8-sig")


@pytest.mark.parametrize("position,label,context,expected", [
    (0, "View BinderScope", "folder:1", True),
    (1, "View BinderScope", "folder:1", False),
    (True, "View BinderScope", "folder:1", False),
    (-1, "View BinderScope", "folder:1", False),
    (0, "not the row action", "folder:1", False),
    (0, "View BinderScope", "old-folder:1", False),
])
def test_row_request_checks_position_label_and_current_folder(position, label, context, expected, monkeypatch):
    opened = []
    monkeypatch.setattr("mn_molecular_design_studio.binderscope_ui._saved_row_dialog", lambda *args: opened.append(args))
    request = {"entry_id": "design", "label": "View BinderScope", "assessments": [("evaluation", "observation")]}
    app = AppTest.from_string("import streamlit as st\n"
        "from mn_molecular_design_studio.binderscope_ui import request_saved_row,render_saved_row_request\n"
        f"st.session_state['test-click']={{'row':{position!r},'label':{label!r}}}\n"
        f"request_saved_row('test-click',{context!r},[{request!r}])\n"
        "render_saved_row_request(None,None,None,[{'id':'design','outcome':{'outcome_kind':'Design'}}],context='folder:1')").run()
    assert not app.exception and bool(opened) is expected


@pytest.mark.parametrize("change", ["delete", "failed", "private", "other-model"])
def test_row_dialog_rechecks_bound_report_and_permissions(database, design, monkeypatch, change):
    identifier = add_saved(database, design)
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    choices = [(choice["evaluation_id"], choice["observation_id"]) for choice in inspector_choices(row)]
    opened = []
    monkeypatch.setattr("mn_molecular_design_studio.binderscope_inspector.render_inspector", lambda *args, **kwargs: opened.append(True))
    with database.connection() as connection:
        if change == "delete":
            connection.execute("DELETE FROM molecule_observations WHERE id=?", (identifier,))
        elif change == "failed":
            connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.status','failed') WHERE id=?", (identifier,))
        elif change == "other-model":
            connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.representation_id','other',"
                "'$.source_structure_id','other','$.source_structure_digest',?,'$.structure_digest',?) WHERE id=?", ("f" * 64, "f" * 64, identifier))
    user = OTHER if change == "private" else USER
    source = ("from pathlib import Path\nfrom mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.binderscope_ui import _saved_row_dialog\n"
        f"_saved_row_dialog(ReadOnlyLibrary(Path({str(database.root)!r})),{user!r},{index!r},{entries[0]!r},{choices!r},'row-test')")
    app = AppTest.from_string(source).run()
    assert not app.exception and not opened and "no longer accessible" in app.warning[0].value


def test_row_dialog_requires_a_choice_even_if_only_one_of_two_assessments_survives(database, design, monkeypatch):
    first = add_saved(database, design, run="first")
    second = add_saved(database, design, run="second")
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    choices = [(choice["evaluation_id"], choice["observation_id"]) for choice in inspector_choices(row)]
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_observations WHERE id=?", (first,))
    opened = []
    monkeypatch.setattr("mn_molecular_design_studio.binderscope_inspector.render_inspector", lambda result, **kwargs: opened.append(result))
    source = ("from pathlib import Path\nfrom mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.binderscope_ui import _saved_row_dialog\n"
        f"_saved_row_dialog(ReadOnlyLibrary(Path({str(database.root)!r})),{USER!r},{index!r},{entries[0]!r},{choices!r},'row-test')")
    app = AppTest.from_string(source).run()
    assert not app.exception and not opened and app.selectbox[0].value is None
    app.selectbox[0].set_value(second).run()
    assert not app.exception and opened[0]["_studio_sql_binding"]["observation_id"] == second
