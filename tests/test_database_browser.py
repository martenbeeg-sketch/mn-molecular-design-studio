from __future__ import annotations

import json

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.database_browser import (
    browse_subjects,
    build_browser_index,
    molecular_subjects,
    run_subjects,
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
    app = AppTest.from_string(f"from mn_molecular_design_studio.database_ui import render_database\nrender_database({USER!r}, [])", default_timeout=15).run()
    assert not app.exception
    assert app.session_state["database-browser-view"] == "Molecules"
    assert next(item for item in app.metric if item.label == "Grouped items").value == "1"
    picker = next(item for item in app.selectbox if item.label == "Inspect grouped item")
    subject = molecular_subjects(build_browser_index(database, USER))[0]
    picker.set_value(subject["id"]).run()
    assert not app.exception
    assert any("2 jobs" in item.value for item in app.caption)
    assert not hydration_calls
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
    app = AppTest.from_string(f"import streamlit as st\nst.session_state.setdefault('database-browser-view', 'Runs')\nfrom mn_molecular_design_studio.database_ui import render_database\nrender_database({USER!r}, [])", default_timeout=15).run()
    assert not app.exception
    assert next(item for item in app.metric if item.label == "Grouped items").value == "1"
    app.selectbox(key="database-browser-run-grouping").set_value("Individual runs").run()
    assert not app.exception
    assert next(item for item in app.metric if item.label == "Grouped items").value == "2"
    app.text_input(key="database-browser-query").set_value("no-match").run()
    assert not app.exception
    assert any("No matching grouped records" in item.value for item in app.info)


def test_grouped_pagination_and_search_include_late_page_records(database, monkeypatch):
    from mn_molecular_design_studio import database_ui

    _import(database, records=[_record("ACDE" * (index + 1), name=f"Binder {index:03}") for index in range(60)])
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    app = AppTest.from_string(f"from mn_molecular_design_studio.database_ui import render_database\nrender_database({USER!r}, [])", default_timeout=15).run()
    assert not app.exception
    assert next(item for item in app.metric if item.label == "Grouped items").value == "60"
    app.number_input(key="database-browser-page").set_value(2).run()
    assert not app.exception
    assert len(next(item for item in app.selectbox if item.label == "Inspect grouped item").options) == 11
    app.text_input(key="database-browser-query").set_value("Binder 059").run()
    assert not app.exception
    assert app.session_state["database-browser-page"] == 1
    assert next(item for item in app.metric if item.label == "Grouped items").value == "1"
