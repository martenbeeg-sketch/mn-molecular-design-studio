from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import datasets, workspace
from mn_molecular_design_studio.datasets_ui import dataset_tree
from test_database_browser import USER, database
from test_database_workspace import _app, _library, _select_rows


def test_visible_import_requires_no_project_and_preserves_chemical_identity(database, monkeypatch):
    app = _app(database, monkeypatch, inputs=True)
    assert not app.exception
    assert app.title[0].value == "Datasets"
    app.button(key="datasets-import").click().run()
    app.text_input(key="molecular-import-group-name").set_value("General screening compounds").run()
    assert app.selectbox(key="molecular-import-project").value == ""
    app.button(key="molecular-import-next").click().run()
    app.radio(key="molecular-import-method").set_value("SMILES").run()
    app.text_area(key="molecular-import-smiles").set_value("CCO ethanol\nCCO duplicate").run()
    next(item for item in app.checkbox if item.label.startswith("I reviewed the report")).set_value(True).run()
    app.button(key="molecular-import-commit").click().run()
    assert not app.exception
    group = database.groups(USER)[0]
    assert group["project_id"] == "" and group["molecule_count"] == 1
    assert datasets.is_input_dataset(group)
    assert app.session_state["datasets-selected"] == group["id"]
    assert database.list_molecules(USER)["items"][0]["metadata"]["dataset_input"] is True


def test_project_shortcut_prefills_reference_without_owning_new_dataset(database, monkeypatch):
    project = workspace.create_project("Docking project", owner_id=USER["id"])
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: [project["id"]])
    app = _app(database, monkeypatch, inputs=True, state={"datasets-project": project["id"]}, projects=[project])
    app.button(key="datasets-import").click().run()
    app.text_input(key="molecular-import-group-name").set_value("Reusable inputs").run()
    assert app.selectbox(key="molecular-import-project").value == project["id"]
    app.button(key="molecular-import-next").click().run()
    app.radio(key="molecular-import-method").set_value("Paste sequences").run()
    app.text_area(key="molecular-import-sequences").set_value(">Target\nACDEFGHIKLMNPQRSTVWY").run()
    app.button(key="molecular-import-commit").click().run()
    assert not app.exception
    group = database.groups(USER)[0]
    assert group["project_id"] == ""
    assert datasets.project_references(group) == [project["id"]]


def test_dataset_member_selection_reuses_original_ids_and_links_to_results(database, monkeypatch):
    saved = _library(database)
    app = _app(database, monkeypatch, inputs=True)
    _select_rows(app, "datasets-members-", [1])
    assert app.session_state["datasets-page-selection"]["ids"]
    assert not app.button(key="data-browser-selection-actions").disabled
    app.button(key="datasets-member-results").click().run()
    assert not app.exception
    assert app.session_state["page"] == "Data Browser"
    assert app.session_state["molecular-preselect-group"] == saved["group_id"]
    assert app.session_state["molecular-preselect-molecule"] in saved["molecule_ids"]
    assert database.groups(USER)[0]["molecule_count"] == 2


def test_reference_editor_and_project_filter_keep_one_dataset(database, monkeypatch):
    first, second = [workspace.create_project(name, owner_id=USER["id"]) for name in ("Docking", "Cofolding")]
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: [first["id"], second["id"]])
    saved = _library(database)
    app = _app(database, monkeypatch, inputs=True, projects=[first, second])
    app.button(key="datasets-project-references").click().run()
    key = next(item.key for item in app.multiselect if item.key.startswith("datasets-project-links-"))
    app.multiselect(key=key).set_value([first["id"], second["id"]]).run()
    app.button(key="datasets-save-project-links").click().run()
    assert not app.exception
    app.selectbox(key="datasets-project").set_value(second["id"]).run()
    assert not app.exception
    assert app.session_state["datasets-selected"] == saved["group_id"]
    assert database.group(saved["group_id"], USER)["project_id"] == ""
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecules").fetchone()[0] == 2


def test_project_tree_uses_project_identity_not_identical_display_names():
    group = {"id": "dataset", "name": "Compounds", "molecule_count": 1000,
        "category": "Small-molecule libraries", "project_references": ["first", "second"]}
    nodes, choices = dataset_tree([group], by_project=True, projects=[
        {"id": "first", "name": "Research"}, {"id": "second", "name": "Research"}])
    assert len(nodes) == 2 and nodes[0]["id"] != nodes[1]["id"]
    assert len(choices) == 2 and set(choices.values()) == {"dataset"}


def test_data_browser_has_a_visible_handoff_instead_of_a_competing_input_importer(database, monkeypatch):
    app = _app(database, monkeypatch)
    assert not app.exception
    assert not any(button.key == "data-browser-import" for button in app.button)
    app.button(key="data-browser-open-datasets").click().run()
    assert not app.exception
    assert app.session_state["page"] == "Datasets"


def test_selected_input_export_does_not_offer_scientific_ranking(database):
    saved = _library(database)
    entries = [{"id": molecule_id, "result": {"molecule_id": molecule_id}} for molecule_id in saved["molecule_ids"]]
    app = AppTest.from_string("from pathlib import Path\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.datasets import DatasetView\n"
        "from mn_molecular_design_studio.database_actions_ui import render_selection_actions\n"
        f"view=DatasetView(ReadOnlyLibrary(Path({str(database.root)!r})), {saved['group_id']!r})\n"
        f"render_selection_actions(view, {USER!r}, [], {entries!r})").run()
    app.button(key="data-browser-selection-actions").click().run()
    assert not app.exception
    assert not any(toggle.key == "data-browser-compare-selection" for toggle in app.toggle)
    app.button(key="data-browser-molecule-export-prepare").click().run()
    assert not app.exception and app.get("download_button")


def test_unstructured_input_members_are_paginated_not_hidden(database, monkeypatch):
    database.import_records(USER, [{"name": "Sequence " + str(position), "components": [
        {"chain": "A", "type": "protein", "sequence": "A" * (position + 4)}]} for position in range(121)],
        name="Reusable sequences", origin=datasets.input_origin(USER))
    app = _app(database, monkeypatch, inputs=True)
    assert not app.exception
    table = next(frame for frame in app.dataframe if "datasets-members-" in frame.proto.id)
    assert len(table.value) == 100 and not table.value["Structures"].any()
    key = next(item.key for item in app.number_input if item.key.startswith("datasets-member-page-"))
    app.number_input(key=key).set_value(2).run()
    assert not app.exception
    assert len(next(frame for frame in app.dataframe if "datasets-members-" in frame.proto.id).value) == 21
