from __future__ import annotations

import json

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.database_browser import (
    browse_subjects,
    build_browser_index,
    component_subjects,
    load_result_record,
    molecular_subjects,
    run_subjects,
    subject_results,
)
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def database(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "workspace"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "workspace"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    return MolecularDatabase(tmp_path / "database")


def _record(sequence="ACDE", **values):
    return {"name": "Binder", "components": [{"chain": "A", "type": "protein", "sequence": sequence}], **values}


def _app_source(view="Molecules", mode="All molecular entities"):
    return (f"import streamlit as st\nst.session_state.setdefault('database-browser-view', {view!r})\n"
            f"st.session_state.setdefault('database-browser-molecule-grouping-v3', {mode!r})\n"
            f"from mn_molecular_design_studio.database_ui import render_database\nrender_database({USER!r}, [])")


def _import(database, *, name="Design", run_id="run-1", task="design", user=USER, origin=None, records=None, method="tool_output", kind="generated", project_id=""):
    return database.import_records(user, records or [_record()], name=name, import_method=method, project_id=project_id,
                                   origin={"collection_kind": kind, "run_id": run_id, "task_group": task,
                                           "tool": "generator", **(origin or {})})


def test_exact_molecule_collects_different_jobs_without_mutating_scientific_evidence(database, monkeypatch):
    first = _import(database, name="First design", run_id="first")
    second = _import(database, name="Second design", run_id="second")
    molecule_id = first["molecule_ids"][0]
    assert second["molecule_ids"] == [molecule_id]
    database.add_observation(first["group_id"], molecule_id, "tool", "first", {"score": 1},
                             context={"task_group": "design", "run_id": "first", "settings": {"seed": 1}})
    database.add_observation(second["group_id"], molecule_id, "tool", "second", {"score": 8},
                             context={"task_group": "design", "run_id": "second", "settings": {"seed": 2}})
    with database.connection() as connection:
        before = [tuple(row) for row in connection.execute("SELECT * FROM molecule_observations ORDER BY id")]
        revision = connection.execute("SELECT local_revision FROM molecular_database_sync_state").fetchone()[0]
    monkeypatch.setattr(database, "list_molecules", lambda *args, **kwargs: pytest.fail("Overview must not hydrate evaluations or coordinates"))
    index = build_browser_index(database, USER)
    subjects = molecular_subjects(index)
    assert len(subjects) == 1
    assert len(subjects[0]["run_keys"]) == 2
    assert subjects[0]["group_ids"] == {first["group_id"], second["group_id"]}
    assert "scores" not in subjects[0]
    with database.connection() as connection:
        assert before == [tuple(row) for row in connection.execute("SELECT * FROM molecule_observations ORDER BY id")]
        assert revision == connection.execute("SELECT local_revision FROM molecular_database_sync_state").fetchone()[0]


def test_sequence_view_links_complexes_but_whole_molecule_view_preserves_partner_and_stoichiometry(database):
    ligand = {"chain": "L", "type": "small_molecule", "smiles": "CCO"}
    _import(database, run_id="alone")
    _import(database, run_id="complex", records=[_record(components=[{"chain": "B", "type": "protein", "sequence": "ACDE"}, ligand])])
    _import(database, run_id="dimer", records=[_record(components=[{"chain": chain, "type": "protein", "sequence": "ACDE"} for chain in ("A", "B")])])
    index = build_browser_index(database, USER)
    assert len(molecular_subjects(index)) == 3
    sequences = molecular_subjects(index, individual_sequences=True)
    assert len(sequences) == 1
    assert sequences[0]["sequence"] == "ACDE"
    assert len(sequences[0]["molecule_ids"]) == 3
    assert len(sequences[0]["run_keys"]) == 3


def test_composition_overview_links_identity_variants_without_replacing_record_ids(database):
    first = _import(database, records=[_record(identity_context={"structure": "first"})])
    second = _import(database, records=[_record(identity_context={"structure": "second"})], run_id="other")
    assert first["molecule_ids"] != second["molecule_ids"]
    subject = molecular_subjects(build_browser_index(database, USER))[0]
    assert subject["molecule_ids"] == set(first["molecule_ids"] + second["molecule_ids"])


def test_incomplete_legacy_contents_are_never_combined_by_missing_identity():
    index = {"molecules": {key: {"name": key, "names": {key}, "components": [], "group_ids": set(),
                                "run_keys": set(), "created_at": ""} for key in ("first", "second")}}
    assert len(molecular_subjects(index)) == 2


def test_sequence_types_and_ligand_stereochemistry_remain_distinct(database):
    for kind in ("protein", "dna", "rna"):
        _import(database, records=[_record(components=[{"chain": "A", "type": kind, "sequence": "AAA"}])])
    for smiles in ("N[C@H](C)C(=O)O", "N[C@@H](C)C(=O)O"):
        _import(database, records=[_record(components=[{"chain": "L", "type": "small_molecule", "smiles": smiles}])])
    index = build_browser_index(database, USER)
    assert len(molecular_subjects(index)) == 5
    assert len(molecular_subjects(index, individual_sequences=True)) == 3


def test_shared_evaluation_group_uses_row_provenance_not_stale_group_origin(database):
    source = _import(database)
    first_id = source["molecule_ids"][0]
    second = _import(database, run_id="second", records=[_record("FGHI")])
    second_id = second["molecule_ids"][0]
    shared = database.create_group(USER, "Evaluated structures", import_method="tool_output",
                                    origin={"run_id": "stale", "task_group": "old-task", "collection_kind": "evaluation"})
    database.add_to_group(USER, shared, [first_id, second_id])
    for molecule_id, run_id in ((first_id, "eval-first"), (second_id, "eval-second")):
        database.add_observation(shared, molecule_id, "tool", run_id, {"score": 1}, context={"task_group": "validation", "run_id": run_id})
    index = build_browser_index(database, USER)
    runs = {row["run_id"]: row for row in index["runs"].values()}
    assert "stale" not in runs
    assert runs["eval-first"]["molecule_ids"] == {first_id}
    assert runs["eval-second"]["molecule_ids"] == {second_id}


def test_run_view_combines_collections_not_homonymous_jobs_or_private_data(database):
    first = _import(database, name="Generated")
    second = _import(database, name="Selected", method="user_uploads", kind="curated")
    peer = _import(database, name="Peer", origin={"workstation": "Peer"})
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET origin=json_set(origin,'$._replica_node','peer-node') WHERE id=?", (peer["group_id"],))
    _import(database, name="Other task", task="validation")
    foreign = _import(database, name="Private foreign job", user=OTHER)
    index = build_browser_index(database, USER)
    assert foreign["group_id"] not in index["groups"]
    runs = run_subjects(index, combine_campaigns=False)
    assert len(runs) == 3
    local = next(row for row in runs if row["group_ids"] == {first["group_id"], second["group_id"]})
    assert len(local["run_keys"]) == 1
    assert len(local["molecule_ids"]) == 1


def test_only_explicit_parent_lineage_combines_campaign_jobs(database):
    _import(database, name="Campaign", run_id="parent")
    _import(database, name="Child", run_id="child", origin={"parent_campaign_run_dir": "/not/read/runs/design/parent"})
    _import(database, name="Campaign", run_id="unrelated", origin={"campaign_name": "Campaign", "registration_batch": "same"})
    index = build_browser_index(database, USER)
    subjects = run_subjects(index)
    assert len(subjects) == 2
    combined = next(row for row in subjects if len(row["run_keys"]) == 2)
    assert combined["name"] == "Campaign"
    assert len(run_subjects(index, combine_campaigns=False)) == 3


def test_unknown_node_ids_do_not_combine_different_workstations(database):
    _import(database, origin={"workstation": "older-peer-one"})
    _import(database, origin={"workstation": "older-peer-two"})
    assert len(run_subjects(build_browser_index(database, USER))) == 2


def test_conflicting_or_cyclic_campaign_links_are_not_guessed(database):
    _import(database, run_id="first", origin={"parent_campaign_run_dir": "/runs/design/second"})
    _import(database, run_id="second", origin={"parent_campaign_run_dir": "/runs/design/first"})
    _import(database, run_id="third", origin={"campaign_id": "campaign"})
    _import(database, run_id="third", origin={"campaign_id": "other-campaign"})
    index = build_browser_index(database, USER)
    assert len(run_subjects(index)) == 3
    assert all(len(row["run_keys"]) == 1 for row in run_subjects(index))


def test_search_happens_after_grouping_and_includes_all_aliases_and_sequences(database):
    _import(database, name="Original library", records=[_record(name="First alias")])
    _import(database, name="Renamed result", run_id="late-run", records=[_record(name="Other alias")])
    index = build_browser_index(database, USER)
    subjects = molecular_subjects(index)
    assert len(browse_subjects(subjects, index, query="OTHER ACDE late-run")) == 1
    assert len(browse_subjects(subjects, index, query="Original library")) == 1
    assert not browse_subjects(subjects, index, query="absent")


def test_hidden_collections_and_revoked_projects_do_not_leak_into_subject_history(database):
    member = {**USER, "_cluster_project_ids": ["project"]}
    _import(database, user=member, project_id="project", name="Restricted project")
    public = _import(database, name="Personal")
    index = build_browser_index(database, USER)
    assert set(index["groups"]) == {public["group_id"]}
    assert len(index["runs"]) == 1
    assert not browse_subjects(molecular_subjects(index), index, query="Restricted")
    assert not build_browser_index(database, USER, methods=["__none__"])["molecules"]


def test_grouped_browser_default_and_original_collection_drilldown(database, monkeypatch):
    from mn_molecular_design_studio import database_ui

    first = _import(database, name="First campaign")
    _import(database, name="Second campaign", run_id="second")
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    monkeypatch.setattr(database_ui, "sync_workspace", lambda *args, **kwargs: pytest.fail("Reads must not start ingestion"))
    monkeypatch.setattr(database_ui, "render_details", lambda *args, **kwargs: None)
    hydration_calls = []
    original_list = database.list_molecules

    def list_molecules(*args, **kwargs):
        hydration_calls.append(kwargs)
        return original_list(*args, **kwargs)

    monkeypatch.setattr(database, "list_molecules", list_molecules)
    app = AppTest.from_string(_app_source(), default_timeout=15).run()
    assert not app.exception
    assert app.session_state["database-browser-view"] == "Molecules"
    assert any("1 molecules" in item.value for item in app.caption)
    picker = next(item for item in app.radio if item.label == "Choose a molecule")
    subject = component_subjects(build_browser_index(database, USER))[0]
    picker.set_value(subject["id"]).run()
    assert not app.exception
    assert any("2 results" in item.value for item in app.caption)
    assert hydration_calls
    assert all(call["molecule_ids"] == first["molecule_ids"] and call["limit"] == 1 for call in hydration_calls)
    app.button(key="database-browser-open-" + subject["id"]).click().run()
    assert not app.exception
    assert app.session_state["database-browser-view"] == "Collections"
    assert app.session_state["database-selected-group"] in subject["group_ids"]
    assert app.session_state["database-focus"] == first["molecule_ids"][0]


def test_grouped_browser_run_switch_search_and_empty_state(database, monkeypatch):
    from mn_molecular_design_studio import database_ui

    _import(database, run_id="parent")
    _import(database, run_id="child", origin={"parent_campaign_run_dir": "/runs/design/parent"})
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    app = AppTest.from_string(_app_source("Runs"), default_timeout=15).run()
    assert not app.exception
    assert any("1 runs / campaigns" in item.value for item in app.caption)
    app.selectbox(key="database-browser-run-grouping").set_value("Individual runs").run()
    assert not app.exception
    assert any("2 runs / campaigns" in item.value for item in app.caption)
    app.text_input(key="database-browser-query").set_value("no-match").run()
    assert not app.exception
    assert any("No matching molecules or runs" in item.value for item in app.info)


def test_grouped_pagination_and_search_include_late_page_records(database, monkeypatch):
    from mn_molecular_design_studio import database_ui

    _import(database, records=[_record("ACDE" * (index + 1), name=f"Binder {index:03}") for index in range(60)])
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    app = AppTest.from_string(_app_source(), default_timeout=15).run()
    assert not app.exception
    assert any("60 molecules" in item.value for item in app.caption)
    app.number_input(key="database-browser-page").set_value(2).run()
    assert not app.exception
    assert len(next(item for item in app.radio if item.label == "Choose a molecule").options) == 10
    app.text_input(key="database-browser-query").set_value("Binder 059").run()
    assert not app.exception
    assert app.session_state["database-browser-page"] == 1
    assert any("1 molecules" in item.value for item in app.caption)


@pytest.mark.parametrize("view", ["Molecules", "Runs", "Collections"])
def test_grouped_views_keep_filters_left_and_automatic_viewer_right(database, monkeypatch, view):
    import streamlit as st
    from mn_molecular_design_studio import database_ui

    imported = _import(database)
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    monkeypatch.setattr(database_ui, "render_details", lambda database, user, record: st.markdown("VIEWER_ACTIVE: " + record["id"]))
    app = AppTest.from_string(_app_source(view), default_timeout=15).run()
    assert not app.exception
    columns = app.get("column")
    left, center, right = columns[0], next(column for column in columns if column.weight == .50), columns[-1]
    assert left.weight == .25 and right.weight == .25
    if view == "Collections":
        assert left.text_input(key="molecular-query").label == "Search molecules"
        assert left.selectbox(key="molecular-sort").label == "Sort"
        assert left.selectbox(key="molecular-type-filter").label == "Molecule type"
        assert {"Target", "Tool", "Run"}.issubset({item.label for item in left.selectbox})
    else:
        assert left.text_input(key="database-browser-query").label in {"Find a molecule", "Find a run / campaign"}
        assert any(item.label in {"Group by", "Browse molecules"} for item in left.selectbox)
        assert any(item.label == "Inspect result" for item in left.selectbox)
        assert any(item.label == "Find a result" for item in left.text_input)
        assert any(item.label in {"Choose a molecule", "Choose a run / campaign"} for item in left.radio)
        assert not center.radio
    assert not center.text_input
    assert {item.label for item in center.selectbox}.issubset({"Format"})
    assert any("VIEWER_ACTIVE: " + imported["molecule_ids"][0] == item.value for item in right.markdown)
    assert not any(item.label == "Load molecule inspector" for item in app.toggle)
    if view != "Collections":
        app.text_input(key="database-browser-query").set_value("no-match").run()
        assert not app.exception
        assert not any("VIEWER_ACTIVE" in item.value for item in app.markdown)
        assert any("Select a result to inspect" in item.value for item in app.caption)


def test_grouped_sheet_focus_changes_the_inspector_without_stale_selection(database, monkeypatch):
    import streamlit as st
    from mn_molecular_design_studio import database_ui, molecule_table

    imported = _import(database, run_id="first", records=[_record("ACDE", name="First")])
    _import(database, run_id="second", records=[_record("ACDE", name="Second")])
    for run_id, score in (("first", 1), ("second", 2)):
        group = next(group for group in database.groups(USER) if group["origin"].get("run_id") == run_id)
        database.add_observation(group["id"], imported["molecule_ids"][0], "tool", run_id, {"score": score}, context={"task_group": "design"})
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    monkeypatch.setattr(database_ui, "render_details", lambda database, user, record: st.markdown("VIEWER_ACTIVE: " + str(record["scores"].get("tool.score"))))
    index = build_browser_index(database, USER)
    subject = component_subjects(index)[0]
    second = next(row for row in subject_results(database, index, subject) if row["run_id"] == "second")

    def render_table(*args, **kwargs):
        assert kwargs["allow_selection"] is False
        return {"token": "clicked-second", "focusId": second["id"], "selectedIds": []}

    monkeypatch.setattr(molecule_table, "render_table", render_table)
    app = AppTest.from_string(_app_source(), default_timeout=15).run()
    assert not app.exception
    assert app.session_state["database-focus"] == imported["molecule_ids"][0]
    assert any("VIEWER_ACTIVE: 2" == item.value for item in app.markdown)
    next(item for item in app.text_input if item.label == "Find a result").set_value("first").run()
    assert not app.exception
    assert app.session_state["database-focus"] == imported["molecule_ids"][0]
    assert any("VIEWER_ACTIVE: 1" == item.value for item in app.markdown)
    assert not any("VIEWER_ACTIVE: 2" == item.value for item in app.markdown)


@pytest.mark.parametrize("allow_selection", [True, False])
def test_molecular_sheet_preserves_selection_default_and_supports_focus_only(monkeypatch, allow_selection):
    from mn_molecular_design_studio import molecule_table

    arguments = []
    monkeypatch.setattr(molecule_table, "_component", lambda **kwargs: arguments.append(kwargs))
    options = {} if allow_selection else {"allow_selection": False}
    molecule_table.render_table([], [], selected_ids=[], focus_id="", key="sheet-test", **options)
    assert arguments[0]["allowSelection"] is allow_selection


def test_component_navigation_connects_ligands_complexes_and_explicit_target_designs(database):
    target = _import(database, records=[_record("ACDE", name="Target")], method="user_uploads", kind="curated", run_id="")
    target_id = target["molecule_ids"][0]
    ligand = {"chain": "L", "type": "small_molecule", "smiles": "CCO"}
    compound = _import(database, records=[_record(name="Ligand", components=[ligand])], run_id="ligand")
    complex_record = _import(database, records=[_record(name="Complex", components=[
        {"chain": "T", "type": "protein", "sequence": "ACDE"}, ligand])], run_id="complex")
    design = _import(database, records=[_record("FGHI", name="Designed binder", metadata={"design_target_molecule_id": target_id})], run_id="design")
    unrelated = _import(database, records=[_record("KLMN", name="Target")], run_id="unrelated")
    index = build_browser_index(database, USER, include_related=True)
    subjects = component_subjects(index)
    protein = next(row for row in subjects if row["sequence"] == "ACDE")
    assert protein["name"] == "Protein target · 4 aa · ACDE"
    assert protein["molecule_ids"] == set(target["molecule_ids"] + complex_record["molecule_ids"] + design["molecule_ids"])
    assert not protein["molecule_ids"] & set(unrelated["molecule_ids"])
    compound_subject = next(row for row in subjects if row["type"] == "small_molecule")
    assert compound_subject["molecule_ids"] == set(compound["molecule_ids"] + complex_record["molecule_ids"])
    results = subject_results(database, index, protein)
    assert set(row["molecule_id"] for row in results) == protein["molecule_ids"] - set(target["molecule_ids"])
    assert next(row for row in results if row["molecule_id"] == design["molecule_ids"][0])["relation"] == "Recorded target link"


def test_target_history_includes_accessible_nonlibrary_results_without_other_target_evaluations(database):
    target = _import(database, records=[_record("ACDE", name="Target")], method="user_uploads", kind="curated", run_id="")
    other_target = _import(database, records=[_record("KLMN", name="Other target")], method="user_uploads", kind="curated", run_id="")
    design = _import(database, records=[_record("FGHI", name="Binder")], run_id="score", kind="evaluation")
    for target_id, model in ((target["molecule_ids"][0], "0"), (other_target["molecule_ids"][0], "1")):
        database.add_observation(design["group_id"], design["molecule_ids"][0], "validator", "score", {"score": int(model)},
            context={"task_group": "design", "target_molecule_id": target_id, "model_id": model})
    private = _import(database, user=OTHER, records=[_record("PQRS", metadata={"design_target_molecule_id": target["molecule_ids"][0]})])
    index = build_browser_index(database, USER, methods=["user_uploads"], include_related=True)
    subject = next(row for row in component_subjects(index) if row["sequence"] == "ACDE")
    assert design["molecule_ids"][0] not in index["visible_ids"]
    assert design["molecule_ids"][0] in subject["molecule_ids"]
    assert not set(private["molecule_ids"]) & subject["molecule_ids"]
    rows = subject_results(database, index, subject)
    evaluations = [row for row in rows if row["evaluation_id"]]
    assert len(evaluations) == 1 and evaluations[0]["model"] == "0"
    record = load_result_record(database, USER, evaluations[0])
    assert record["selected_evaluation"]["target_id"] == target["molecule_ids"][0]
    assert record["selected_evaluation"]["model_id"] == "0"
    assert len(record["observations"]) == 1
    assert load_result_record(database, OTHER, evaluations[0]) is None


def test_run_results_do_not_include_other_jobs_of_the_same_molecule(database):
    original = _import(database, name="Shared evaluations", run_id="", kind="evaluation")
    molecule_id = original["molecule_ids"][0]
    for run_id, workstation, model in (("first", "peer-one", "0"), ("second", "peer-one", "1"), ("first", "peer-two", "2")):
        database.add_observation(original["group_id"], molecule_id, "predictor", run_id, {"score": int(model)},
            context={"task_group": "validation", "workstation": workstation, "model_id": model})
    index = build_browser_index(database, USER, collection="all", include_related=True)
    subjects = run_subjects(index, combine_campaigns=False)
    assert len(subjects) == 3
    selected = next(row for row in subjects if any(index["runs"][key]["run_id"] == "first"
                    and "peer-one" in index["runs"][key]["workstations"] for key in row["run_keys"]))
    rows = subject_results(database, index, selected, run_view=True)
    assert len(rows) == 1
    assert rows[0]["run_id"] == "first" and rows[0]["model"] == "0"
    assert rows[0]["workstation"] == "peer-one"


def test_selected_evaluations_and_structures_keep_their_own_scores(database):
    imported = _import(database, records=[_record(files=[{"name": "model.pdb", "content": b"END\n", "model_id": "0"}])])
    molecule_id = imported["molecule_ids"][0]
    database.add_observation(imported["group_id"], molecule_id, "predictor", "run-1", {"score": 10},
                             context={"task_group": "design", "model_id": "0"})
    database.add_observation(imported["group_id"], molecule_id, "predictor", "run-1", {"score": 20},
                             context={"task_group": "design", "model_id": "1"})
    unscored = _import(database, run_id="unscored", records=[_record(files=[{"name": "unscored.pdb", "content": b"END\n\n"}])])
    index = build_browser_index(database, USER)
    subject = component_subjects(index)[0]
    with database.connection() as connection:
        before = [tuple(row) for row in connection.execute("SELECT * FROM molecule_observations ORDER BY id")]
        revision = connection.execute("SELECT local_revision FROM molecular_database_sync_state").fetchone()[0]
    rows = subject_results(database, index, subject)
    evaluation_rows = [row for row in rows if row["evaluation_id"]]
    assert len(evaluation_rows) == 2
    for row in evaluation_rows:
        record = load_result_record(database, USER, row)
        assert record["selected_evaluation"]["id"] == row["evaluation_id"]
        assert record["scores"]["predictor.score"] == (10 if row["model"] == "0" else 20)
        assert len(record["structures"]) == (1 if row["model"] == "0" else 0)
    structure_row = next(row for row in rows if row["kind"] == "Structure" and row["group_id"] == unscored["group_id"])
    record = load_result_record(database, USER, structure_row)
    assert len(record["structures"]) == 1
    assert record["selected_evaluation"] is None and not record["scores"] and not record["observations"]
    with database.connection() as connection:
        assert before == [tuple(row) for row in connection.execute("SELECT * FROM molecule_observations ORDER BY id")]
        assert revision == connection.execute("SELECT local_revision FROM molecular_database_sync_state").fetchone()[0]


def test_source_collections_do_not_duplicate_one_run_output(database):
    first = _import(database, name="Generated")
    second = _import(database, name="Selected", method="user_uploads", kind="curated")
    index = build_browser_index(database, USER, include_related=True)
    subject = component_subjects(index)[0]
    rows = subject_results(database, index, subject)
    assert len(rows) == 1
    assert rows[0]["group_ids"] == {first["group_id"], second["group_id"]}
    assert rows[0]["collection"] == "Generated, Selected"
    database.add_observation(first["group_id"], first["molecule_ids"][0], "generator", "run-1", {"score": 4},
                             context={"task_group": "design"})
    index = build_browser_index(database, USER, include_related=True)
    rows = subject_results(database, index, component_subjects(index)[0])
    assert len(rows) == 1 and rows[0]["kind"] == "Evaluation"


def test_result_pages_search_before_pagination_and_hydrate_only_one_record(database, monkeypatch):
    from mn_molecular_design_studio import database_ui

    imported = _import(database)
    for model in range(60):
        database.add_observation(imported["group_id"], imported["molecule_ids"][0], "predictor", "run-1", {"score": model},
                                 context={"task_group": "design", "model_id": str(model), "candidate_id": f"candidate-{model:03}"})
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    calls = []
    original = database.list_molecules

    def hydrate(*args, **kwargs):
        calls.append(kwargs)
        return original(*args, **kwargs)

    monkeypatch.setattr(database, "list_molecules", hydrate)
    app = AppTest.from_string(_app_source(), default_timeout=15).run()
    assert not app.exception
    assert any("60 connected results" in item.value for item in app.caption)
    next(item for item in app.number_input if item.label == "Results page").set_value(2).run()
    assert not app.exception
    assert len(next(item for item in app.selectbox if item.label == "Inspect result").options) == 10
    next(item for item in app.text_input if item.label == "Find a result").set_value("candidate-059").run()
    assert not app.exception
    assert any("1 connected results" in item.value for item in app.caption)
    assert all(call["limit"] == 1 and call["molecule_ids"] == imported["molecule_ids"] for call in calls)


def test_native_evaluation_selector_follows_center_result_for_same_molecule(database, monkeypatch):
    from mn_molecular_design_studio import database_ui

    imported = _import(database)
    for model in ("0", "1"):
        database.add_observation(imported["group_id"], imported["molecule_ids"][0], "predictor", "run-1", {"score": int(model)},
                                 context={"task_group": "design", "model_id": model})
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    index = build_browser_index(database, USER)
    rows = subject_results(database, index, component_subjects(index)[0])
    app = AppTest.from_string(_app_source(), default_timeout=15).run()
    assert not app.exception
    for row in rows:
        next(item for item in app.selectbox if item.label == "Inspect result").set_value(row["id"]).run()
        assert not app.exception
        selector = next(item for item in app.selectbox if item.label == "Evaluation")
        assert selector.value == row["evaluation_id"] and len(selector.options) == 1


def test_targets_use_recorded_names_and_hide_repeated_sequence_carrier_inputs(database):
    target = _import(database, name="Protein Library", run_id="run-1", method="tool_inputs", kind="reference",
        origin={"input_role": "receptor", "collection_kind": "inputs"},
        records=[_record(name="idr_sequence_carrier", metadata={"target_identity_type": "sequence_only"})])
    target_id = target["molecule_ids"][0]
    design = _import(database, run_id="run-1", records=[_record("FGHI", name="bindcraft2_00001",
        metadata={"design_target_molecule_id": target_id, "design_target_name": "R3 (test target)"})])
    database.add_observation(design["group_id"], design["molecule_ids"][0], "designer", "run-1", {"score": 1},
        context={"task_group": "design", "target_molecule_id": target_id, "candidate_id": "candidate"})
    index = build_browser_index(database, USER, include_related=True)
    targets = component_subjects(index, role="targets")
    assert len(targets) == 1 and targets[0]["name"] == "R3 (test target)" and targets[0]["role"] == "Target"
    assert target_id in targets[0]["target_record_ids"]
    assert design["molecule_ids"][0] not in targets[0]["target_record_ids"]
    outputs = subject_results(database, index, targets[0])
    assert len(outputs) == 1 and outputs[0]["evaluation_id"]
    assert outputs[0]["name"] == "bindcraft2_00001"
    assert all(row["kind_label"] != "Input / reference" for row in outputs)
    history = subject_results(database, index, targets[0], include_inputs=True)
    assert any(row["kind_label"] == "Input / reference" for row in history)
    assert len(component_subjects(index, role="designs")) == 1


def test_target_component_never_borrows_a_binder_alias_and_preserves_member_roles(database):
    target = _import(database, name="References", run_id="", kind="reference", method="user_uploads",
        records=[_record(name="Receptor", metadata={"target_identity_type": "sequence_only"})])
    complex_record = _import(database, records=[_record(name="bindcraft2_00001_ipsae_001", components=[
        {"chain": "A", "type": "protein", "sequence": "ACDE"},
        {"chain": "B", "type": "protein", "sequence": "FGHI"}], metadata={"source_to_prediction_components": {
            "target": [{"prediction_chain": "A", "prediction_sequence": "ACDE"}],
            "designed_binder": [{"prediction_chain": "B", "prediction_sequence": "FGHI"}]}})])
    index = build_browser_index(database, USER, include_related=True)
    targets = component_subjects(index, role="targets")
    assert len(targets) == 1 and targets[0]["name"] == "Receptor"
    assert targets[0]["sequence"] == "ACDE"
    assert complex_record["molecule_ids"][0] in targets[0]["molecule_ids"]
    assert component_subjects(index, role="designs")[0]["sequence"] == "FGHI"
    assert target["molecule_ids"][0] in targets[0]["target_record_ids"]


@pytest.mark.parametrize("ambiguous", [False, True])
def test_legacy_design_links_use_only_unambiguous_origin_qualified_target_inputs(database, ambiguous):
    target = _import(database, name="Target reference", method="tool_inputs", kind="reference",
        origin={"input_role": "receptor"}, records=[_record(name="Target reference")])
    design = _import(database, records=[_record("FGHI", name="Generated binder")])
    if ambiguous:
        _import(database, name="Other target", method="tool_inputs", kind="reference",
            origin={"input_role": "receptor"}, records=[_record("KLMN", name="Other target")])
    foreign = _import(database, name="Peer output", records=[_record("PQRS")], origin={"workstation": "other-peer"})
    index = build_browser_index(database, USER, include_related=True)
    selected = next(row for row in component_subjects(index, role="targets") if row["sequence"] == "ACDE")
    outputs = subject_results(database, index, selected)
    assert any(row["molecule_id"] == design["molecule_ids"][0] for row in outputs) is not ambiguous
    assert not any(row["molecule_id"] == foreign["molecule_ids"][0] for row in outputs)
    assert target["molecule_ids"][0] in selected["target_record_ids"]


def test_protein_targets_and_ligands_keep_distinct_poses_and_designed_binders(database):
    receptor = _import(database, name="Protein receptor", run_id="", method="user_uploads", kind="reference",
        records=[_record(name="Protein receptor", metadata={"target_identity_type": "sequence_only"})])
    ligand = _import(database, name="Ligand pose records", run_id="", task="docking", records=[
        _record(name="Compound", components=[{"chain": "L", "type": "small_molecule", "smiles": "CCO"}])])
    for pose in ("1", "2"):
        database.add_observation(ligand["group_id"], ligand["molecule_ids"][0], "dock", "docking", {"score": int(pose)},
            context={"task_group": "docking", "target_molecule_id": receptor["molecule_ids"][0], "pose_id": pose})
    designed = _import(database, name="Ligand binder", records=[_record("FGHI", metadata={
        "design_target_molecule_id": ligand["molecule_ids"][0], "design_target_name": "Ligand target"})])
    database.add_observation(designed["group_id"], designed["molecule_ids"][0], "designer", "run-1", {"score": 8},
        context={"task_group": "design", "target_molecule_id": ligand["molecule_ids"][0]})
    index = build_browser_index(database, USER, include_related=True)
    subjects = component_subjects(index, role="targets")
    assert len(subjects) == 2
    protein = next(row for row in subjects if row["type"] == "protein")
    compound = next(row for row in subjects if row["type"] == "small_molecule")
    protein_results = [row for row in subject_results(database, index, protein) if row["evaluation_id"]]
    assert len(protein_results) == 2 and {row["model"] for row in protein_results} == {"1", "2"}
    compound_results = [row for row in subject_results(database, index, compound) if row["evaluation_id"]]
    assert len(compound_results) == 3
    assert not any(row["sequence"] == "FGHI" for row in subjects)


def test_default_ui_selects_targets_not_designs_and_real_outputs_not_input_history(database, monkeypatch):
    import streamlit as st
    from mn_molecular_design_studio import database_ui

    target = _import(database, name="Reference", method="tool_inputs", kind="reference", origin={"input_role": "receptor"},
        records=[_record(name="idr_sequence_carrier", metadata={"target_identity_type": "sequence_only"})])
    designed = _import(database, records=[_record("FGHI", name="Binder", metadata={
        "design_target_molecule_id": target["molecule_ids"][0], "design_target_name": "R3 reference"})])
    database.add_observation(designed["group_id"], designed["molecule_ids"][0], "designer", "run-1", {"score": 1},
        context={"task_group": "design", "target_molecule_id": target["molecule_ids"][0]})
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    monkeypatch.setattr(database_ui, "render_details", lambda database, user, record: st.markdown("VIEWER: " + record["name"]))
    app = AppTest.from_string(_app_source(mode="Targets & ligands"), default_timeout=15).run()
    assert not app.exception
    assert any("R3 reference" in item.value for item in app.markdown)
    assert any("VIEWER: Binder" == item.value for item in app.markdown)
    assert not any("idr_sequence_carrier" in option for item in app.radio for option in item.options)
    next(item for item in app.checkbox if item.label == "Include input / reference history").check().run()
    assert not app.exception
    assert any("Input / reference" in option for item in app.selectbox if item.label == "Result type" for option in item.options)


def test_reused_binder_without_evaluation_target_uses_its_own_source_collection_link(database):
    first = _import(database, name="First target", method="user_uploads", kind="reference", run_id="", records=[_record("ACDE")])
    second = _import(database, name="Second target", method="user_uploads", kind="reference", run_id="", records=[_record("KLMN")])
    for target, run_id in ((first, "first"), (second, "second")):
        design = _import(database, run_id=run_id, records=[_record("FGHI", metadata={"design_target_molecule_id": target["molecule_ids"][0]})])
        database.add_observation(design["group_id"], design["molecule_ids"][0], "predictor", run_id, {"score": 1},
                                 context={"task_group": "design"})
    index = build_browser_index(database, USER, include_related=True)
    selected = next(row for row in component_subjects(index, role="targets") if row["sequence"] == "ACDE")
    rows = [row for row in subject_results(database, index, selected) if row["evaluation_id"]]
    assert len(rows) == 1 and rows[0]["run_id"] == "first"
