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


def test_input_chemistry_and_sequence_previews_are_visible_without_calculations(database, monkeypatch):
    _library(database)
    app = _app(database, monkeypatch, inputs=True)
    assert not app.exception
    assert app.image
    assert any(code.value in {"CCO", "CO"} for code in app.code)
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    imported = database.import_records(USER, [{"name": "Reference sequence", "components": [
        {"chain": "A", "type": "protein", "sequence": sequence}]}], name="Sequence inputs", origin=datasets.input_origin(USER))
    app = _app(database, monkeypatch, inputs=True, state={"datasets-open-group": imported["group_id"]})
    assert not app.exception
    assert any(code.value == sequence for code in app.code)


def test_saved_reference_structure_is_rendered_as_cartoon_without_result_scores(database, monkeypatch):
    from test_molecular_database_policy import _pdb_acd
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
    imported = database.import_records(USER, [{"name": "Reference structure", "components": [
        {"chain": "A", "type": "protein", "sequence": "ACD"}]}], name="Input structures", origin=datasets.input_origin(USER))
    structure = database.attach_structure(imported["group_id"], imported["molecule_ids"][0], {
        "name": "reference.pdb", "content": _pdb_acd().encode()})
    rendered = []
    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda visualizations, **kwargs: rendered.extend(visualizations))
    app = _app(database, monkeypatch, inputs=True, state={"datasets-open-group": imported["group_id"]})
    assert not app.exception
    assert any(item.label == "Saved structure" and item.value == structure for item in app.selectbox)
    assert len(rendered) == 1 and rendered[0].representation_type == "cartoon"
    assert app.get("download_button")
    assert not any("Original result scores" in item.value for item in app.markdown)


def test_input_ligand_preview_does_not_borrow_a_calculated_receptor(database, monkeypatch):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from mn_molecular_design_studio import data_browser_ui
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    imported = _library(database)
    database.attach_structure(imported["group_id"], imported["molecule_ids"][0], {
        "name": "reference.pdb", "content": Chem.MolToPDBBlock(molecule).encode()},
        context={"database_structure_policy": "best_docking_pose"})
    rendered = []
    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda visualizations, **kwargs: rendered.extend(visualizations))
    monkeypatch.setattr(data_browser_ui, "saved_receptor_file", lambda *args: (_ for _ in ()).throw(
        AssertionError("Input previews must not hydrate calculation contexts")))
    app = _app(database, monkeypatch, inputs=True, state={"datasets-focus-molecule": imported["molecule_ids"][0]})
    assert not app.exception
    assert rendered and all(item.representation_type != "cartoon" for item in rendered)


def test_protein_import_is_registered_and_automatic_reference_groups_are_hidden(database, monkeypatch):
    from test_molecular_database_policy import _pdb_acd

    automatic = database.create_group(USER, "Automatic Protein Library", import_method="tool_inputs",
        origin={"type": "job", "collection_kind": "reference"})
    monkeypatch.setattr("mn_molecular_design_studio.molecule_imports.fetch_identifier",
        lambda identifier, source, *, structure_format: (identifier.lower() + ".pdb", _pdb_acd().encode()))
    app = _app(database, monkeypatch, inputs=True)
    app.button(key="datasets-import-protein").click().run()
    app.text_input(key="molecular-import-group-name").set_value("Imported protein").run()
    app.button(key="molecular-import-next").click().run()
    assert app.radio(key="molecular-import-method").value == "PDB accession"
    app.text_input(key="molecular-import-identifier").set_value("1LYZ").run()
    app.button(key="molecular-import-fetch").click().run()
    assert not app.button(key="molecular-import-commit").disabled
    app.button(key="molecular-import-commit").click().run()
    assert not app.exception
    groups = datasets.catalog(database, USER)
    assert len(groups) == 1 and groups[0]["name"] == "Imported protein"
    assert database.group(automatic, USER)["name"] == "Automatic Protein Library"
    from mn_molecular_design_studio.target_editor import saved_protein_structure_options
    assert len(saved_protein_structure_options(database, USER)) == 1
