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


def test_structure_picker_uses_registered_inputs_and_keeps_the_original_bytes(database, monkeypatch, tmp_path):
    from mn_molecular_design_studio import dataset_structure_input, molecule_imports
    from test_molecular_database_policy import _pdb_acd

    content = _pdb_acd().encode()
    database.import_records(USER, molecule_imports.parse_structure(content, "intermediate.pdb"),
                            name="Unregistered output", import_method="tool_output", origin={"type": "job"})
    example = tmp_path / "public-example.pdb"
    example.write_bytes(content)
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(dataset_structure_input, "MolecularDatabase", lambda: database)
    app = AppTest.from_string(
        "import streamlit as st\n"
        "from mn_molecular_design_studio.dataset_structure_input import register_example, select_dataset_files\n"
        f"st.session_state['_molecular_database_user'] = {USER!r}\n"
        f"if st.button('Load example'): register_example({str(example)!r}, {USER['id']!r})\n"
        "selection = select_dataset_files('Coordinates', key='coordinates')\n"
        "if selection: st.session_state['input_bytes'] = selection.getvalue()\n"
    ).run()
    assert not app.exception and not app.get("file_uploader")
    for _ in range(2):
        app.button[0].click().run()
        assert not app.exception
    options = dataset_structure_input.options(database, USER)
    assert len(options) == 1
    assert options[0]["name"] == "public-example.pdb"
    app.selectbox[0].set_value(options[0]["id"]).run()
    assert not app.exception and app.session_state["input_bytes"] == content

    app.button(key="coordinates-datasets").click().run()
    assert app.session_state["page"] == "Datasets"
    assert app.session_state["molecular-import-kind"] == "structure"


def test_preparation_saves_a_separate_trim_with_original_residue_numbers(database, monkeypatch):
    from mn_molecular_design_studio import dataset_preparation_ui, molecule_imports
    from test_molecular_database_policy import _pdb_acd

    expander = dataset_preparation_ui.st.expander
    monkeypatch.setattr(dataset_preparation_ui.st, "expander", lambda label, **kwargs: expander(label, expanded=True, **kwargs))
    content = _pdb_acd().encode()
    imported = database.import_records(USER, datasets.input_records(molecule_imports.parse_structure(content, "target.pdb")),
                                       name="Target inputs", origin=datasets.input_origin(USER))
    group_id = imported["group_id"]
    app = AppTest.from_string(
        "from pathlib import Path\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.datasets import DatasetView\n"
        "from mn_molecular_design_studio.dataset_preparation_ui import render_preparation\n"
        f"user = {USER!r}\n"
        f"view = DatasetView(ReadOnlyLibrary(Path({str(database.root)!r})), {group_id!r})\n"
        f"record = view.list_molecules(user, molecule_ids={imported['molecule_ids']!r}, limit=1)['items'][0]\n"
        "render_preparation(view, user, record)\n"
    ).run()
    assert not app.exception
    app.slider[0].set_value((2, 3)).run()
    next(button for button in app.button if button.label == "Save prepared variant").click().run()
    assert not app.exception
    view = datasets.DatasetView(database, group_id)
    original, cropped = sorted(view.list_molecules(USER)["items"], key=lambda row: row["name"])
    source = original["structures"][0]
    prepared = cropped["structures"][0]
    assert view.structure_file(USER, source["id"]).read_bytes() == content
    assert cropped["components"][0]["sequence"] == "CD"
    assert prepared["context"]["derived_from_structure_id"] == source["id"]
    assert "CYS A   2" in view.structure_file(USER, prepared["id"]).read_text()


def test_import_records_explicit_smiles_replacement_and_original(database, monkeypatch):
    app = _app(database, monkeypatch, inputs=True)
    app.button(key="datasets-import-compound").click().run()
    app.text_input(key="molecular-import-group-name").set_value("Corrected screening compounds").run()
    app.button(key="molecular-import-next").click().run()
    app.radio(key="molecular-import-method").set_value("SMILES").run()
    app.text_area(key="molecular-import-smiles").set_value("C1CC incomplete\nCCO ethanol").run()
    correction = next(item for item in app.text_area if item.label == "Replacement SMILES")
    correction.set_value("1 C1CC1").run()
    next(button for button in app.button if button.label == "Validate corrections").click().run()
    assert not app.exception and not app.error
    next(item for item in app.checkbox if item.label.startswith("I reviewed the report")).set_value(True).run()
    app.button(key="molecular-import-commit").click().run()
    assert not app.exception and not app.error
    with database.connection() as connection:
        import json

        rows = [json.loads(row["provenance"]) for row in connection.execute("SELECT provenance FROM molecule_sources")]
    assert any(row.get("raw_smiles") == "C1CC" and row.get("corrected_smiles") == "C1CC1" for row in rows)


def test_completed_preparation_can_be_saved_in_datasets_without_reupload(database, monkeypatch, tmp_path):
    import json

    from mn_molecular_design_studio.preparation_history import sha256
    from test_molecular_database_policy import _pdb_acd

    prepared = tmp_path / "prepared_target.pdb"
    prepared.write_text(_pdb_acd(), encoding="utf-8")
    digest = sha256(prepared.read_bytes())
    report = tmp_path / "repair_report.json"
    report.write_text(json.dumps({
        "output_sha256": {prepared.name: digest},
        "loop_remodelling": {"enabled": True, "chain": "A", "start": 2, "end": 2},
    }), encoding="utf-8")
    app = AppTest.from_string(
        "import streamlit as st\n"
        "from pathlib import Path\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.dataset_preparation_ui import render_prepared_outputs\n"
        "from mn_molecular_design_studio.datasets_ui import render_datasets\n"
        f"user = {USER!r}\n"
        "if st.session_state.get('page') == 'Datasets':\n"
        f"    render_datasets(user, [], database=ReadOnlyLibrary(Path({str(database.root)!r})))\n"
        "else:\n"
        f"    render_prepared_outputs([({{}}, Path({str(prepared)!r})), ({{}}, Path({str(report)!r}))], 'repair-1', '')\n"
    ).run()
    app.button(key="prepared-output-save-repair-1").click().run()
    assert not app.exception
    app.text_input(key="molecular-import-group-name").set_value("Repaired targets").run()
    app.button(key="molecular-import-next").click().run()
    app.button(key="molecular-import-commit").click().run()
    assert not app.exception and not app.error
    group_id = database.groups(USER)[0]["id"]
    saved = datasets.DatasetView(database, group_id).list_molecules(USER)["items"][0]
    context = saved["structures"][0]["context"]
    assert context["preparation_run_id"] == "repair-1"
    history = context["preparation_history"]
    assert history["structure_sha256"] == digest
    assert history["steps"] == [{
        "label": "MODELLER existing-loop remodelling", "details": {"chain": "A", "start": 2, "end": 2},
    }]
    assert history["artifacts"][0]["report_sha256"] == sha256(report.read_bytes())
