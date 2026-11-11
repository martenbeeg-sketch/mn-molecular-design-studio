"""Target/system → campaign navigation uses saved identities, not names or files."""
import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_catalog import (
    collection_system_facets, retained_system_contexts, saved_relationship_page,
)

OWNER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


def test_exact_target_sequence_groups_across_engines_without_using_designed_binder():
    def group(engine, binder, label):
        return {"origin": {"tool": engine}, "scientific_context": {"targets": [{
            "label": label, "sequences_by_chain": {"A": "ACDE"}, "components": [
                {"type": "protein", "chain": "A", "sequence": "ACDE"},
                {"type": "protein", "chain": "B", "sequence": binder}]}]}}
    a, = collection_system_facets(group("rfd4", "GGG", "RFD4 input"))
    b, = collection_system_facets(group("bindcraft2", "WWWWW", "BindCraft input"))
    assert a["key"] == b["key"]
    assert a["input_kind"] == "Sequence"
    assert a["target"]["sequences_by_chain"] == {"A": "ACDE"}
    assert "components" not in a["target"]


def test_construct_preparation_structure_and_unknown_inputs_stay_distinct():
    def facet(**target):
        return collection_system_facets({"scientific_context": {"targets": [target]}})
    assert not facet(label="Tau")
    assert not facet(label="Tau", model_id="model_0")
    variants = [facet(label="Tau", target_id="tau", sequences_by_chain={"A": "ACDE"}),
                facet(label="Tau", target_id="tau", sequences_by_chain={"A": "ACDEF"}),
                facet(label="Tau", target_id="tau", preparation="prep-1"),
                facet(label="Tau", target_id="tau", preparation="prep-2"),
                facet(label="Tau", target_structure_digest="a" * 64),
                facet(label="Tau", target_structure_digest="b" * 64)]
    assert len({rows[0]["key"] for rows in variants}) == 6
    both, = facet(sequences_by_chain={"A": "ACDE"}, target_structure_digest="a" * 64)
    assert both["input_kind"] == "Sequence + structure"
    described, = facet(target_id="tau", preparation={"type": "engine_specific"})
    assert described["label"].endswith("preparation recorded")
    named, = facet(target_id="tau", preparation={"run_id": "prep-1"})
    assert named["label"].endswith("prep prep-1")


def test_equal_protein_sequence_cannot_erase_recorded_ligand_or_cofactor():
    def key(smiles):
        facets = collection_system_facets({"scientific_context": {"targets": [{
            "sequences_by_chain": {"A": "ACDE"}, "components": [
                {"type": "protein", "chain": "A", "sequence": "ACDE"},
                {"type": "ligand", "chain": "L", "smiles": smiles}]}]}})
        return facets[0]
    a, b = key("CCO"), key("CCN")
    assert a["key"] != b["key"]
    assert a["target"]["components"][0]["smiles"] == "CCO"


@pytest.fixture
def systems(tmp_path):
    db = MolecularDatabase(tmp_path / "db")
    imported = db.import_records(OWNER, molecule_imports.parse_smiles("CCO compound"), name="MD systems")
    group, molecule = imported["group_id"], imported["molecule_ids"][0]
    for i, (target, prep) in enumerate((("alpha", "prep-1"), ("alpha", "prep-2"), ("beta", "prep-3"))):
        db.add_observation(group, molecule, "md.mmgbsa.gb", "run-" + str(i), {"energy": -1},
            context={"target_id": target, "prepared_target_run_id": prep,
                     "study_id": "study-" + str(i), "target_receptor": "Same name"})
    return db, group


def test_system_facet_navigation_preserves_each_study_and_rechecks_acl(systems):
    db, group = systems
    contexts = retained_system_contexts(db, OWNER, [group])[group]
    assert len(contexts) == 3
    facets = collection_system_facets({"retained_system_contexts": contexts})
    assert len({facet["key"] for facet in facets}) == 3
    scope = [contexts[1]["scope"]]
    page = saved_relationship_page(db, OWNER, [group], system_scope=scope)
    assert len(page["items"]) == 1
    assert page["items"][0]["preparation_id"] == scope[0]["preparation_id"]
    assert page["totals"]["studies"] == 1
    assert not saved_relationship_page(db, OWNER, [group], system_scope=scope, query="nonmatching")["items"]
    assert not saved_relationship_page(db, OWNER, [group], system_scope=[])["items"]
    for call in (retained_system_contexts, saved_relationship_page):
        with pytest.raises(PermissionError):
            call(db, OTHER, [group])
    with pytest.raises(ValueError):
        saved_relationship_page(db, OWNER, [group], system_scope=[{**scope[0], "group_id": "unrelated"}])


def test_overview_and_scoped_evidence_do_not_read_coordinates(systems, monkeypatch):
    from streamlit.testing.v1 import AppTest
    db, group = systems
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: db)
    monkeypatch.setattr(db, "structure_file", lambda *args: pytest.fail("No coordinate loading during navigation"))
    scope = [retained_system_contexts(db, OWNER, [group])[group][2]["scope"]]
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_collection_ui import render_collection_workspace
from mn_molecular_design_studio.complex_analysis_evidence_ui import open_saved_work
if 'complex-analysis-evidence-selection' not in st.session_state:
    open_saved_work(st.session_state, user_id={OWNER['id']!r}, group_ids=[{group!r}], label='Saved MD',
                    system_label='Beta preparation', system_scope={scope!r})
render_collection_workspace({OWNER!r})
''').run()
    assert not app.exception
    assert any("Campaign / study overview" in item.value for item in app.markdown)
    assert any("Beta preparation" in item.value for item in app.caption)
    app.button(key="overview-open-models").click().run()
    assert not app.exception
    assert app.session_state["complex-analysis-evidence-selection"]["relationship"]["target_ref"] == "beta"
    app.button(key="complex-analysis-step-analyze").click().run()
    assert not app.exception
    selected = app.session_state["complex-analysis-evidence-selection"]
    assert selected["relationship"]["target_ref"] == "beta"
    assert selected["relationship"]["preparation_id"] == "prep-3"


def test_scope_resolution_is_bounded_to_visible_model_ids_and_preserves_relationship_contract(systems, tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_catalog as catalog
    db, group = systems
    molecule = db.list_molecules(OWNER)["items"][0]["id"]
    model = db.attach_structure(group, molecule, {"name": "model.cif", "content": b"data_fixture\n"})
    db.add_observation(group, molecule, "biotite_interface", "assessed-model", {"contacts": 0},
                       context={"representation_id": model, "target_id": "beta", "prepared_target_run_id": "prep-3"})
    scope = [{"group_id": group, "target_ref": "beta", "target_revision": "", "preparation_id": "prep-3"}]
    old = set()
    for row in catalog.saved_relationship_page(db, OWNER, [group], system_scope=scope)["items"]:
        old.update(catalog.relationship_structure_ids(db, OWNER, [group], row["relationship"]))
    monkeypatch.setattr(catalog, "saved_relationship_page", lambda *a, **kw: pytest.fail("all-relationship paging"))
    monkeypatch.setattr(db, "structure_file", lambda *a, **kw: pytest.fail("coordinate hydration"))
    assert catalog.scoped_structure_ids(db, OWNER, [group], scope, [model]) == old == {model}
    assert catalog.scoped_structure_ids(db, OWNER, [group], scope, ["another-visible-model"]) == set()
    with pytest.raises(PermissionError):
        catalog.scoped_structure_ids(db, OTHER, [group], scope, [model])


def test_campaign_pages_and_exact_system_selection_are_metadata_only(monkeypatch):
    from streamlit.testing.v1 import AppTest
    groups = [{"id": f"g-{i}", "name": f"Campaign {i:02}", "owner_id": OWNER["id"], "project_id": "p",
               "origin": {"type": "job", "tool": "rfd4_proteina" if i % 2 else "bindcraft2",
                          "run_id": f"run-{i}", "target_molecule_id": "target-b" if i == 26 else "target-a"}}
              for i in range(27)]
    class Database:
        def groups(self, *args, **kwargs):
            return groups
        def collection_statistics(self, user, ids):
            return {key: {"molecule_count": 1, "structure_count": 1} for key in ids}
        def cohort_statistics(self, user, cohorts):
            return {key: {"molecule_count": len(ids), "structure_count": len(ids)} for key, ids in cohorts.items()}
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", Database)
    monkeypatch.setattr("mn_molecular_design_studio.complex_analysis_catalog.saved_evidence_summary", lambda *args: {})
    monkeypatch.setattr("mn_molecular_design_studio.complex_analysis_collection_ui._results_source_jobs",
                        lambda *args: pytest.fail("Campaign browsing must not query Results"))
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.complex_analysis_collection_ui import _database_collection_samples
_database_collection_samples({OWNER!r})
''').run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 25
    next(widget for widget in app.number_input if widget.label == "Campaign page").set_value(2).run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 2
    picker = app.selectbox(key="complex-analysis-work-system")
    assert len(picker.options) == 3  # All systems + two exact identities, across both engines.
    b = next(option for option in picker.options if "target-b" in option)
    app.session_state["complex-analysis-collections-open"] = True
    picker.select(b).run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 1
    assert not app.session_state["complex-analysis-collections-open"]


def test_project_picker_uses_names_and_distinct_owner_scope(monkeypatch):
    from streamlit.testing.v1 import AppTest
    project_id = '1a4f7ce287d348d9af6a79bc80cbfeff'
    groups = [{'id': f'g-{i}', 'name': name, 'owner_id': OWNER['id'], 'project_id': project,
               'origin': {'type': 'job', 'tool': 'bindcraft2', 'run_id': f'run-{i}'}}
              for i, (name, project) in enumerate([('Named campaign', project_id), ('Personal campaign', '')])]
    class Database:
        def groups(self, *a, **kw):
            return groups
        def collection_statistics(self, user, ids):
            return {key: {'molecule_count': 1, 'structure_count': 1} for key in ids}
        def cohort_statistics(self, user, cohorts):
            return {key: {'molecule_count': len(ids), 'structure_count': len(ids)} for key, ids in cohorts.items()}
    monkeypatch.setattr('mn_molecular_design_studio.molecular_database.MolecularDatabase', Database)
    monkeypatch.setattr('mn_molecular_design_studio.complex_analysis_catalog.saved_evidence_summary', lambda *a: {})
    monkeypatch.setattr('mn_molecular_design_studio.workspace.load_organization', lambda: {
        'projects': [{'id': project_id, 'name': 'Michele_tau'}]})
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.complex_analysis_collection_ui import _database_collection_samples
_database_collection_samples({OWNER!r})
''').run()
    assert not app.exception
    picker = app.selectbox(key='complex-analysis-work-project')
    assert picker.options == ['All accessible projects', 'Michele_tau', 'Personal / owner-scoped']
    picker.select(project_id).run()
    assert not app.exception
    assert app.session_state['complex-analysis-work-project-value'] == project_id
    assert app.dataframe[0].value['Work group'].tolist() == ['Named campaign']
    app.selectbox(key='complex-analysis-work-project').select('__owner_scoped__').run()
    assert not app.exception
    assert app.dataframe[0].value['Work group'].tolist() == ['Personal campaign']
