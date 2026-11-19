import json
import sqlite3

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import entity_names as names, workspace
from mn_molecular_design_studio.data_browser import ReadOnlyLibrary, folder_scopes
from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects, molecular_subjects
from mn_molecular_design_studio.database_browser_cache import BrowserCache
from mn_molecular_design_studio.database_browser_labels import target_label
from test_database_browser import OTHER, USER, _import, _record, database


def target(database, *, user=USER, sequence="ACDE", name="Old receptor", project_id=""):
    imported = _import(database, user=user, method="user_uploads", kind="reference", project_id=project_id,
        records=[_record(sequence, name=name, metadata={"target_identity_type": "sequence"})])
    index = build_browser_index(database, user, hide_smoke_test=False, include_related=True)
    subject = next(row for row in component_subjects(index, role="targets") if row["sequence"] == sequence)
    return imported, subject


def test_target_rename_preserves_scientific_rows_and_resolves_all_views(database):
    imported, subject = target(database)
    molecule_id = imported["molecule_ids"][0]
    database.add_observation(imported["group_id"], molecule_id, "example", "run", {"score": 1})
    tables = ["molecules", "molecule_members", "molecule_observations", "molecule_sources", "molecule_structures"]
    with database.connection() as connection:
        before = {table: [tuple(row) for row in connection.execute("SELECT * FROM " + table)] for table in tables}
    renamed = names.rename_target(database, USER, subject["id"], "owner:" + USER["id"], "Tau R3", expected_revision="")
    assert renamed["aliases"] == ["Old receptor"]
    with database.connection() as connection:
        assert before == {table: [tuple(row) for row in connection.execute("SELECT * FROM " + table)] for table in tables}
    index = build_browser_index(database, USER, hide_smoke_test=False, include_related=True)
    assert component_subjects(index, role="targets")[0]["name"] == "Tau R3"
    assert molecular_subjects(index)[0]["name"] == "Tau R3"
    assert target_label(index, molecule_id) == "Tau R3"
    folder = next(row for row in folder_scopes(index, []).values() if row["parent"] == "targets")
    assert folder["label"] == "Tau R3" and "Old receptor" in folder["search"]
    for query in ["Tau R3", "Old receptor"]:
        result = database.list_molecules(USER, query=query)
        assert result["total"] == 1
        assert result["items"][0]["name"] == "Tau R3"
        assert result["items"][0]["original_name"] == "Old receptor"


def test_target_cas_restart_and_ingestion_do_not_overwrite_manual_name(database):
    imported, subject = target(database)
    scope = "owner:" + USER["id"]
    first = names.rename_target(database, USER, subject["id"], scope, "R3", expected_revision="")
    with pytest.raises(names.NameConflictError):
        names.rename_target(database, USER, subject["id"], scope, "Stale", expected_revision="")
    target(database, name="New import alias")
    second = names.rename_target(database, USER, subject["id"], scope, "Tau R3", expected_revision=first["revision"])
    assert {"R3", "Old receptor", "New import alias"} <= set(second["aliases"])
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    restarted = MolecularDatabase(database.root)
    assert component_subjects(build_browser_index(restarted, USER), role="targets")[0]["name"] == "Tau R3"
    with restarted.connection() as connection:
        assert connection.execute("SELECT count(*) FROM entity_name_history").fetchone()[0] == 2


def test_target_names_are_private_exact_and_cannot_rename_design_only_entities(database):
    imported, subject = target(database, name="Public receptor")
    _, cropped = target(database, sequence="ACD", name="Different crop")
    _, private = target(database, user=OTHER, name="Private receptor")
    names.rename_target(database, OTHER, private["id"], "owner:" + OTHER["id"], "Secret target", expected_revision="")
    with pytest.raises(PermissionError):
        names.rename_target(database, USER, subject["id"], "owner:" + OTHER["id"], "Stolen", expected_revision="")
    names.rename_target(database, USER, subject["id"], "owner:" + USER["id"], "Public target", expected_revision="")
    index = build_browser_index(database, USER)
    assert {row["name"] for row in component_subjects(index, role="targets")} == {"Public target", "Different crop"}
    assert "Secret" not in json.dumps(index["entity_names"])
    assert database.list_molecules(USER, query="Secret target")["total"] == 0
    _import(database, records=[_record("FGHI", metadata={"output_kind": "design"})])
    with pytest.raises(PermissionError):
        names.rename_target(database, USER, names.component_id({"type": "protein", "sequence": "FGHI"}),
                            "owner:" + USER["id"], "Not a target", expected_revision="")


def test_scope_conflicts_are_not_silently_resolved_alphabetically(database, monkeypatch):
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["p", "q"])
    _, subject = target(database, project_id="p")
    target(database, project_id="q")
    for scope, label in [("project:p", "Tau R3"), ("project:q", "Tau repeat 3")]:
        names.rename_target(database, USER, subject["id"], scope, label, expected_revision="")
    assert component_subjects(build_browser_index(database, USER), role="targets")[0]["name"] == "Tau R3 / Tau repeat 3"
    assert component_subjects(build_browser_index(database, USER, project_id="p"), role="targets")[0]["name"] == "Tau R3"


def test_cache_invalidates_on_sql_rename_and_legacy_reader_does_not_migrate(database):
    _, subject = target(database)
    cache = BrowserCache(database.path)
    reader = ReadOnlyLibrary(database.root)
    before = cache.get_index(reader, USER, scope="all", collection="library", hide_smoke_test=False, include_related=True)
    names.rename_target(database, USER, subject["id"], "owner:" + USER["id"], "Renamed", expected_revision="")
    after = cache.get_index(reader, USER, scope="all", collection="library", hide_smoke_test=False, include_related=True)
    assert after is not before
    assert component_subjects(after, role="targets")[0]["name"] == "Renamed"
    cache.close()
    with database.connection() as connection:
        connection.execute("DROP TABLE entity_names")
        connection.execute("DROP TABLE entity_name_history")
    legacy = build_browser_index(reader, USER)
    assert component_subjects(legacy, role="targets")[0]["name"] == "Old receptor"
    with reader.connection() as connection:
        assert not connection.execute("SELECT 1 FROM sqlite_master WHERE name='entity_names'").fetchone()


def test_project_migration_rename_aliases_and_stale_json_cannot_undo(database, tmp_path):
    home = tmp_path / "organization-home"
    project = workspace.create_project("Original project", home=home, owner_id=USER["id"])
    project_id = project["id"]
    folder = home / "organization" / "projects" / project_id
    (folder / "keep.txt").write_text("research")
    before = workspace.load_organization(home)
    assert before["projects"][0]["name_revision"] == ""
    renamed = names.rename_project(USER, project_id, "  Tau   binders  ", expected_revision="", home=home)
    assert renamed["name"] == "Tau binders" and renamed["name_aliases"] == ["Original project"]
    assert renamed["id"] == project_id and (folder / "keep.txt").read_text() == "research"
    workspace.save_organization(before, home)
    current = workspace.load_organization(home)["projects"][0]
    assert current["name"] == "Tau binders" and current["created_by"] == USER["id"]
    with pytest.raises(names.NameConflictError):
        names.rename_project(USER, project_id, "Stale form", expected_revision="", home=home)
    with sqlite3.connect(home / "database/molecular_library.sqlite3") as connection:
        assert connection.execute("SELECT name FROM entity_names WHERE kind='project'").fetchone()[0] == "Tau binders"


def test_project_owner_duplicate_and_cleanup_guards(database, tmp_path):
    home = tmp_path / "home"
    project = workspace.create_project("First", home=home, owner_id=USER["id"])
    workspace.create_project("Second", home=home, owner_id=USER["id"])
    with pytest.raises(PermissionError):
        names.rename_project(OTHER, project["id"], "Other", expected_revision="", home=home)
    with pytest.raises(ValueError, match="already exists"):
        names.rename_project(USER, project["id"], "second", expected_revision="", home=home)
    organization = workspace.load_organization(home)
    organization["projects"][0]["lifecycle"] = {"status": "running"}
    workspace.save_organization(organization, home)
    with pytest.raises(ValueError, match="cleanup"):
        names.rename_project(USER, project["id"], "Other", expected_revision="", home=home)


@pytest.mark.parametrize("value", ["", "\n", "a\tb", "x" * 201, None, 12])
def test_invalid_names_do_not_write(database, value):
    _, subject = target(database)
    with pytest.raises(ValueError):
        names.rename_target(database, USER, subject["id"], "owner:" + USER["id"], value, expected_revision="")
    with database.connection() as connection:
        assert not names.records(connection)


def test_peer_name_authority_aliases_and_replay(database):
    _, subject = target(database)
    value = names.rename_target(database, USER, subject["id"], "owner:" + USER["id"], "First", expected_revision="")
    key = "entity_name:" + value["id"]
    with database.connection() as connection:
        with pytest.raises(PermissionError, match="authority"):
            names.apply_record(connection, key, {**value, "authority_node": "another-node"}, "another-node")
        next_value = {**value, "name": "Second", "revision": "a" * 32, "aliases": []}
        names.apply_record(connection, key, next_value, value["authority_node"])
        names.apply_record(connection, key, next_value, value["authority_node"])
        assert {"First", "Old receptor"} <= set(names.records(connection)[0]["aliases"])
        assert connection.execute("SELECT count(*) FROM entity_name_history").fetchone()[0] == 2
        with pytest.raises(ValueError):
            names.apply_record(connection, key, None, value["authority_node"])


def test_project_rename_ui_uses_saved_revision_and_searches_alias(database, monkeypatch):
    from mn_molecular_design_studio import projects_ui
    project = {"id": "p", "name": "Current", "name_aliases": ["Former"], "created_by": USER["id"], "name_revision": "base"}
    calls = []
    monkeypatch.setattr(projects_ui, "_names_test_project", project, raising=False)
    monkeypatch.setattr(projects_ui, "_names_test_rename", lambda *args, **kwargs: calls.append((args, kwargs)), raising=False)
    app = AppTest.from_string(f'''
from mn_molecular_design_studio import projects_ui
projects_ui.render_projects({USER!r}, [projects_ui._names_test_project], load_inventory=lambda: {{"jobs":[],"collections":[]}},
    on_create=lambda: None, on_results=lambda value: None, on_database=lambda value: None,
    on_rename=projects_ui._names_test_rename)
''').run()
    assert not app.exception
    app.text_input(key="projects-search").set_value("Former").run()
    assert app.selectbox(key="projects-selected").options == ["Current · 0 runs"]
    next(button for button in app.button if button.label == "Project properties").click().run()
    next(button for button in app.button if button.label == "Edit metadata").click().run()
    assert not app.exception
    project["name_revision"] = "peer-update"
    next(field for field in app.text_input if field.label == "Name").set_value("New name")
    next(button for button in app.button if button.label == "Save changes").click().run()
    assert not app.exception
    assert calls[0][0] == ("p", "New name")
    assert calls[0][1]["expected_revision"] == "base"


def test_target_rename_ui_writes_only_on_save_and_rejects_stale_form(database, monkeypatch):
    from mn_molecular_design_studio import entity_names_ui
    _, subject = target(database)
    monkeypatch.setattr(entity_names_ui, "_names_test_database", database, raising=False)
    app = AppTest.from_string(f'''
from mn_molecular_design_studio import entity_names_ui as ui
from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects
from mn_molecular_design_studio.entity_properties_ui import render_notice
render_notice()
index = build_browser_index(ui._names_test_database, {USER!r})
subject = component_subjects(index, role="targets")[0]
ui.render_target_rename(ui._names_test_database, {USER!r}, index, subject)
''').run()
    assert not app.exception
    with database.connection() as connection:
        assert not names.records(connection)
    next(button for button in app.button if button.label == "Folder properties").click().run()
    next(button for button in app.button if button.label == "Edit metadata").click().run()
    assert not app.exception
    names.rename_target(database, USER, subject["id"], "owner:" + USER["id"], "Peer edit", expected_revision="")
    next(field for field in app.text_input if field.label == "Name").set_value("Stale edit")
    next(button for button in app.button if button.label == "Save changes").click().run()
    assert not app.exception
    assert any("another session" in element.value for element in app.error)
    next(button for button in app.button if button.label == "Reload current values").click().run()
    assert next(field for field in app.text_input if field.label == "Name").value == "Peer edit"
    next(field for field in app.text_input if field.label == "Name").set_value("New reviewed name")
    next(button for button in app.button if button.label == "Save changes").click().run()
    assert not app.exception
    assert component_subjects(build_browser_index(database, USER), role="targets")[0]["name"] == "New reviewed name"
