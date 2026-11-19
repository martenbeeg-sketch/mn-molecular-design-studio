import json

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import entity_names as names, entity_properties as properties, entity_properties_ui as ui, workspace
from mn_molecular_design_studio.data_browser import ReadOnlyLibrary, expand_scope, folder_scopes, listed_results, saved_results
from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects, subject_results
from test_database_browser import OTHER, USER, _import, _record, database
from test_database_design_evaluations import design
from test_database_outcomes import _pass_design
from test_entity_names import target


def scientific_rows(database):
    with database.connection() as connection:
        return {table: [tuple(row) for row in connection.execute("SELECT * FROM " + table + " ORDER BY 1,2")]
                for table in ["molecules", "molecule_groups", "molecule_members", "molecule_observations", "molecule_structures", "molecule_sources"]}


def save(database, snapshot, *, name=None, values=None, actor=USER, outcome=None):
    return properties.save_group_properties(database, actor, snapshot, name or snapshot["name"], values or {}, [],
        expected_revision=snapshot["revision"], outcome=outcome)


def inventory(database):
    index = build_browser_index(database, USER, hide_smoke_test=False, include_related=True)
    subject = {"id": "all", "group_ids": set(index["groups"]), "molecule_ids": set(index["molecules"])}
    return index, saved_results(index, subject_results(database, index, subject))


def test_dataset_properties_project_names_without_rewriting_source_data(database):
    imported, subject = target(database)
    before = scientific_rows(database)
    original = properties.dataset_snapshot(database, USER, imported["group_id"])
    save(database, original, name="HY-L126 study", values={"description": "Imported screening collection", "citation": "Study A", "notes": "Keep original source"})
    assert scientific_rows(database) == before
    assert database.group(imported["group_id"], USER)["name"] == "HY-L126 study"
    index, outcomes = inventory(database)
    dataset = next(scope for scope in folder_scopes(index, []).values() if (scope.get("entity") or {}).get("kind") == "dataset")
    assert dataset["label"] == "HY-L126 study" and "Study A" in dataset["search"]
    with pytest.raises(names.NameConflictError):
        save(database, original, values={"notes": "Stale"})
    current = properties.dataset_snapshot(ReadOnlyLibrary(database.root), USER, imported["group_id"])
    assert current["properties"]["citation"] == "Study A"
    assert current["aliases"] == [original["name"]]


def test_entry_annotations_are_local_to_one_explicit_membership(database):
    first, subject = target(database)
    second, _ = target(database, name="Another entry")
    molecule_id = first["molecule_ids"][0]
    before = scientific_rows(database)
    entry = properties.entry_snapshot(database, USER, first["group_id"], molecule_id)
    save(database, entry, name="Sample A", values={"batch": "Lot 12", "dataset_identifier": "Catalog-1", "notes": "Dataset-local note"})
    assert scientific_rows(database) == before
    assert properties.entry_snapshot(database, USER, first["group_id"], molecule_id)["properties"]["batch"] == "Lot 12"
    assert properties.entry_snapshot(database, USER, second["group_id"], molecule_id)["properties"].get("batch") is None
    assert component_subjects(build_browser_index(database, USER), role="targets")[0]["name"] != "Sample A"


def test_original_notes_changed_by_an_old_writer_reject_stale_entry_form(database):
    imported, subject = target(database)
    entry = properties.entry_snapshot(database, USER, imported["group_id"], imported["molecule_ids"][0])
    database.update_metadata(USER, imported["group_id"], imported["molecule_ids"][0], {"notes": "Another version updated this"})
    with pytest.raises(names.NameConflictError, match="source annotations"):
        save(database, entry, values={"notes": "Stale copy"})


def test_target_and_compound_annotations_preserve_chemical_identity_and_scope(database):
    component = {"type": "small_molecule", "chain": "L", "smiles": "N[C@H](C)C(=O)O"}
    imported = _import(database, method="user_uploads", kind="curated", records=[_record(name="Compound A", components=[component])])
    stereo = _import(database, method="user_uploads", kind="curated", records=[_record(name="Other stereo", components=[{**component, "smiles": "N[C@@H](C)C(=O)O"}])])
    index = build_browser_index(database, USER)
    subject = next(row for row in component_subjects(index, role="targets") if imported["molecule_ids"][0] in row["direct_molecule_ids"])
    before = scientific_rows(database)
    names.rename_target(database, USER, subject["id"], "owner:" + USER["id"], "Alanine reference", expected_revision="",
        properties={"identifiers": "Reference annotation", "notes": "Supplier note"})
    assert scientific_rows(database) == before
    snapshot = properties.target_snapshot(database, USER, subject["id"], "owner:" + USER["id"])
    assert snapshot["type_label"] == "Compound" and snapshot["properties"]["notes"] == "Supplier note"
    index, outcomes = inventory(database)
    assert listed_results(outcomes, query="Supplier note")
    assert all(row["outcome"].get("display_name") != "Alanine reference" for row in listed_results(outcomes)
               if row["outcome"]["molecule_id"] == stereo["molecule_ids"][0])


@pytest.mark.parametrize("kind,forbidden", [("dataset", "project_id"), ("entry", "design_target_molecule_id"),
    ("result", "acceptance_status"), ("target", "sequence"), ("project", "created_by"), ("target", "smiles")])
def test_only_typed_descriptive_fields_are_accepted(kind, forbidden):
    with pytest.raises(ValueError, match="protected"):
        properties.validate_properties(kind, {forbidden: "changed"})


def test_permissions_source_ownership_deletion_and_wrong_member_block_edits(database):
    imported, subject = target(database)
    snapshot = properties.dataset_snapshot(database, USER, imported["group_id"])
    with pytest.raises(PermissionError):
        save(database, snapshot, actor=OTHER, values={"notes": "Not mine"})
    with pytest.raises(PermissionError):
        properties.entry_snapshot(database, USER, imported["group_id"], "missing")
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET origin=? WHERE id=?", (json.dumps({"_replica_node": "peer"}), imported["group_id"]))
    assert not properties.dataset_snapshot(database, USER, imported["group_id"])["writable"]
    with pytest.raises(PermissionError, match="source workstation"):
        save(database, snapshot, values={"notes": "Wrong machine"})


def test_saved_design_annotations_are_separate_from_its_scores_and_admission(database, design):
    _pass_design(database, design)
    index, outcomes = inventory(database)
    outcome = next(row for row in outcomes if row["outcome_kind"] == "Design")
    before = scientific_rows(database)
    snapshot = properties.result_snapshot(database, USER, outcome)
    save(database, snapshot, name="Promising binder", values={"notes": "Review annotation only"}, outcome=outcome)
    assert scientific_rows(database) == before
    index, after = inventory(database)
    result = next(row for row in after if row["id"] == outcome["id"])
    assert result["display_name"] == "Promising binder" and result["eligibility"] == outcome["eligibility"]
    assert listed_results(after, query="Review annotation only")
    with pytest.raises(ValueError):
        save(database, properties.result_snapshot(database, USER, result), values={"score": "100"}, outcome=result)


def test_project_properties_are_atomic_and_old_json_descriptions_cannot_undo_them(database, tmp_path):
    home = tmp_path / "project-home"
    project = workspace.create_project("Original", "Old description", home=home, owner_id=USER["id"])
    original = workspace.load_organization(home)
    snapshot = properties.project_snapshot(USER, original["projects"][0])
    renamed = names.rename_project(USER, project["id"], "New name", home=home, expected_revision=snapshot["revision"],
        expected_source_revision=snapshot["source_revision"], properties={"description": "New description", "notes": "Project notes"})
    workspace.save_organization(original, home)
    current = workspace.load_organization(home)["projects"][0]
    assert current["name"] == "New name" and current["description"] == "New description"
    assert current["entity_properties"]["notes"] == "Project notes"
    with pytest.raises(ValueError):
        names.rename_project(USER, project["id"], "Rejected", home=home, expected_revision=renamed["name_revision"], properties={"owner_id": "someone"})
    assert workspace.load_organization(home)["projects"][0]["name"] == "New name"


def test_virtual_folders_are_not_editable_and_nested_target_keeps_real_identity(database, monkeypatch):
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["p"])
    imported, subject = target(database, project_id="p")
    index = build_browser_index(database, USER)
    scopes = folder_scopes(index, [{"id": "p", "name": "Research"}])
    expand_scope(scopes, index, "project:p")
    assert scopes["project:p:targets"]["entity"] is None
    expand_scope(scopes, index, "project:p:targets")
    nested = next(row for row in scopes.values() if row["parent"] == "project:p:targets")
    assert nested["entity"]["subject_id"] == subject["id"]
    assert nested["subject"]["id"] != subject["id"]


def test_properties_dialog_is_explicit_cancelable_and_does_not_write_on_open(database, monkeypatch):
    imported, subject = target(database)
    monkeypatch.setattr(ui, "_properties_test_database", database, raising=False)
    monkeypatch.setattr(ui, "_properties_test_group", imported["group_id"], raising=False)
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio import entity_properties_ui as ui
ui.render_notice()
st.session_state.setdefault("data-browser-scope", "keep-folder")
st.session_state.setdefault("data-browser-folder", "keep-selection")
scope = {{"id":"dataset", "label":"Dataset", "entity":{{"kind":"dataset", "group_id":ui._properties_test_group}}}}
ui.render_folder_properties(ui._properties_test_database, {USER!r}, {{}}, scope, [])
''').run()
    assert not app.exception and not app.text_area and not app.text_input
    before = scientific_rows(database)
    next(button for button in app.button if button.label == "Folder properties").click().run()
    assert not app.exception and not app.text_input
    next(button for button in app.button if button.label == "Edit metadata").click().run()
    assert not app.exception
    next(field for field in app.text_input if field.label == "Name").set_value("Unsaved name")
    next(button for button in app.button if button.label == "Cancel editing").click().run()
    assert not app.text_input
    next(button for button in app.button if button.label == "Close").click().run()
    assert app.session_state["data-browser-scope"] == "keep-folder"
    assert app.session_state["data-browser-folder"] == "keep-selection"
    assert scientific_rows(database) == before
    with database.connection() as connection:
        assert not names.records(connection)


def test_dialog_revalidates_session_before_reading_or_saving(database, monkeypatch):
    imported, subject = target(database)
    session_user = [USER]
    monkeypatch.setattr("mn_molecular_design_studio.auth.resolve_session", lambda token: session_user[0])
    monkeypatch.setattr(ui, "_properties_test_database", database, raising=False)
    monkeypatch.setattr(ui, "_properties_test_group", imported["group_id"], raising=False)
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio import entity_properties_ui as ui
st.session_state['studio-session-token'] = 'test-session'
ui.render_notice()
scope = {{"id":"dataset", "label":"Dataset", "entity":{{"kind":"dataset", "group_id":ui._properties_test_group}}}}
ui.render_folder_properties(ui._properties_test_database, {USER!r}, {{}}, scope, [])
''').run()
    next(button for button in app.button if button.label == "Folder properties").click().run()
    next(button for button in app.button if button.label == "Edit metadata").click().run()
    session_user[0] = None
    next(field for field in app.text_input if field.label == "Name").set_value("Unauthorized edit")
    next(button for button in app.button if button.label == "Save changes").click().run()
    assert not app.exception
    assert any("session expired" in entry.value for entry in app.error)
    with database.connection() as connection:
        assert not names.records(connection)


def test_project_source_revision_catches_old_writer_description_change(database, tmp_path):
    home = tmp_path / "home"
    project = workspace.create_project("Research", "Original description", home=home, owner_id=USER["id"])
    organization = workspace.load_organization(home)
    snapshot = properties.project_snapshot(USER, organization["projects"][0])
    organization["projects"][0]["description"] = "Updated by an old application"
    workspace.save_organization(organization, home)
    with pytest.raises(names.NameConflictError, match="Project details"):
        names.rename_project(USER, project["id"], "Research", home=home, expected_revision=snapshot["revision"],
            expected_source_revision=snapshot["source_revision"], properties={"description": "Stale copy"})
    assert workspace.load_organization(home)["projects"][0]["description"] == "Updated by an old application"


def test_structure_properties_offer_explicitly_linked_target_for_binder_only_records(database, design, monkeypatch):
    from mn_molecular_design_studio.data_browser import load_saved_record
    _pass_design(database, design)
    index, outcomes = inventory(database)
    outcome = next(row for row in outcomes if row["outcome_kind"] == "Design")
    record = load_saved_record(database, USER, outcome["primary"])
    assert len(record["components"]) == 1 and record["components"][0]["sequence"] == "ACD"
    record["metadata"]["target_molecule_id"] = {"unverified": "Not an identity link"}
    captured = []
    monkeypatch.setattr(ui, "_launch", lambda user, choices, **kwargs: captured.extend(choices))
    monkeypatch.setattr(ui.st, "button", lambda *args, **kwargs: True)
    ui.render_record_properties(database, USER, record, index=index, outcome=outcome)
    target_choices = [choice for choice in captured if choice["label"].startswith("Target ·")]
    assert target_choices
    assert target_choices[0]["load"](USER)["source"]["Identity"] == "EFG"
