"""Evidence-first selection must not depend on coordinate availability."""
import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_catalog import saved_evidence_page, saved_evidence_summary
from mn_molecular_design_studio.complex_analysis_collection_ui import _is_design_collection
from mn_molecular_design_studio.complex_analysis_evidence_ui import evidence_row, open_saved_work
from mn_molecular_design_studio.complex_analysis_batch import offered_assessment_methods
from mn_molecular_design_studio.evaluations import evaluation_descriptor

OWNER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def records(tmp_path):
    db = MolecularDatabase(tmp_path / "db")
    imported = db.import_records(OWNER, molecule_imports.parse_smiles("CCO compound"), name="Study")
    group, molecule = imported["group_id"], imported["molecule_ids"][0]
    context = {"study_id": "study-1", "replica_index": 0, "target_id": "target-1",
               "model_id": "model-1", "source_structure_digest": "a" * 64,
               "_replica_node": "pc004033", "settings": {"method": "gb", "frames": 100}}
    for i in range(3):
        db.add_observation(group, molecule, "md.mmgbsa.gb", "run-" + str(i),
                           {"delta_g_bind_kcal_mol": 0 if i == 0 else -20.0}, context=context)
    db.add_observation(group, molecule, "md.mmgbsa.pb", "run-pb", {"delta_g_bind_kcal_mol": 2.0},
                       context={**context, "settings": {"method": "pb"}})
    return db, group, molecule


def test_structureless_md_evidence_is_discoverable_and_paged(records):
    db, group, molecule = records
    summary = saved_evidence_summary(db, OWNER, [group])[group]
    assert summary == {"observation_count": 4, "methods": {"md.mmgbsa.gb": 3, "md.mmgbsa.pb": 1}}
    assert _is_design_collection({"origin": {"type": "job", "collection_kind": "intermediate"}, "stats": summary})
    first = saved_evidence_page(db, OWNER, [group], namespace="md.mmgbsa.gb", limit=2)
    second = saved_evidence_page(db, OWNER, [group], namespace="md.mmgbsa.gb", limit=2, offset=2)
    assert first["total"] == 3 and first["next_offset"] == 2
    assert len(second["items"]) == 1 and second["next_offset"] is None
    assert not {row["id"] for row in first["items"]} & {row["id"] for row in second["items"]}
    assert all(row["molecule_id"] == molecule for row in [*first["items"], *second["items"]])


def test_evidence_keeps_canonical_identity_provenance_and_measured_zero(records):
    db, group, _ = records
    rows = saved_evidence_page(db, OWNER, [group])["items"]
    for row in rows:
        assert row["evaluation"] == evaluation_descriptor(row)
        assert row["evaluation"]["model_id"] == "model-1"
        assert row["context"]["_replica_node"] == "pc004033"
    zero = next(row for row in rows if row["metrics"]["delta_g_bind_kcal_mol"] == 0)
    assert evidence_row(zero)["delta_g_bind_kcal_mol"] == 0
    assert evidence_row(zero)["Method"] == "MD · MM/GBSA"
    pb = next(row for row in rows if row["namespace"] == "md.mmgbsa.pb")
    assert evidence_row(pb)["Method"] == "MD · MM/PBSA"


def test_evidence_rechecks_authorization_on_every_page(records):
    db, group, _ = records
    for call in (saved_evidence_summary, saved_evidence_page):
        with pytest.raises(PermissionError):
            call(db, OTHER, [group])
        with pytest.raises(PermissionError):
            call(db, OWNER, [group, "unrelated-group"])
    with pytest.raises(ValueError):
        saved_evidence_page(db, OWNER, [group], limit=10000)


def test_open_work_is_metadata_only_and_does_not_invent_acceptance():
    state = {"complex-analysis-workflow-samples": ["old-model"], "complex-analysis-workflow-selection": {}}
    open_saved_work(state, user_id="owner", group_ids=["md-study", "md-study"], label="MD study")
    assert state["complex-analysis-evidence-selection"]["group_ids"] == ["md-study"]
    assert "complex-analysis-workflow-samples" not in state
    assert state["complex-analysis-flow-step"] == "analyze"
    assert "acceptance" not in str(state)


def test_method_choices_follow_molecular_inputs_not_generator_name():
    assert offered_assessment_methods([{"kind": "protein-ligand", "engine": "Boltz-2"}]) == []
    assert offered_assessment_methods([{"kind": "protein-ligand", "engine": "UniDock"}]) == []
    protein = offered_assessment_methods([{"kind": "protein"}])
    assert "analysis.sequence_qc" in protein
    assert "analysis.biotite_interface" not in protein
    ppi = offered_assessment_methods([{"kind": "protein-protein", "engine": "new-engine"}])
    assert "analysis.biotite_interface" in ppi


def test_saved_md_dashboard_opens_and_filters_without_coordinate_access(records, monkeypatch):
    from streamlit.testing.v1 import AppTest
    db, group, _ = records
    monkeypatch.setattr('mn_molecular_design_studio.molecular_database.MolecularDatabase', lambda: db)
    monkeypatch.setattr(db, 'structure_file', lambda *args: pytest.fail('Evidence must not retrieve coordinates'))
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_evidence_ui import render_saved_evidence
render_saved_evidence(MolecularDatabase(), {OWNER!r}, {{"group_ids": [{group!r}], "label": "Saved MD study"}})
''').run()
    assert not app.exception
    assert [metric.value for metric in app.metric] == ['3', '1']
    app.selectbox[0].select('md.mmgbsa.gb').run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 3
    app.radio[0].set_value('Sources & provenance').run()
    assert not app.exception
    assert group in str(app.json[0].value)


def test_project_choice_with_all_projects_survives_navigation():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_collection_ui import _stateful_selectbox
if st.checkbox('Show selection', value=True):
    _stateful_selectbox('Project', [None, 'project-1'], key='project-widget', state_key='project-value', default=None)
''').run()
    assert not app.exception
    app.selectbox[0].set_value('project-1').run()
    app.checkbox[0].uncheck().run()
    app.checkbox[0].check().run()
    assert not app.exception
    assert app.selectbox[0].value == 'project-1'


def test_existing_evaluation_identity_and_target_references_are_reused(records):
    from mn_molecular_design_studio.complex_analysis_catalog import retained_target_references
    db, group, _ = records
    # Normal Database ingestion/query installs its existing Evaluation projection.
    db.list_molecules(OWNER, group_id=group)
    with db.connection() as connection:
        import json
        descriptors = {row['observation_id']: json.loads(row['descriptor']) for row in connection.execute(
            'SELECT observation_id,descriptor FROM molecule_evaluations WHERE group_id=?', (group,))}
    assert descriptors
    assert retained_target_references(db, OWNER, [group]) == {group: ['target-1']}
    for row in saved_evidence_page(db, OWNER, [group])['items']:
        assert row['evaluation'] == descriptors[row['id']]
    with pytest.raises(PermissionError):
        retained_target_references(db, OTHER, [group])


def test_open_saved_work_button_reaches_evidence_dashboard(records, monkeypatch):
    import streamlit as st
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import complex_analysis_collection_ui as ui
    db, group, _ = records
    monkeypatch.setattr('mn_molecular_design_studio.molecular_database.MolecularDatabase', lambda: db)
    monkeypatch.setattr('mn_molecular_design_studio.complex_analysis_catalog.collection_source_contexts', lambda *args: {})
    original = ui._database_analysis_cohorts
    def selected_cohort(groups):
        cohorts = original(groups)
        st.session_state['complex-analysis-collections-open'] = True
        st.session_state['complex-analysis-work-table-scope'] = (None, 'All types', '')
        st.session_state['complex-analysis-selected-cohorts'] = [cohorts[0]['id']]
        st.session_state['complex-analysis-selected-collections'] = [group]
        return cohorts
    monkeypatch.setattr(ui, '_database_analysis_cohorts', selected_cohort)
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.complex_analysis_collection_ui import render_collection_workspace
render_collection_workspace({OWNER!r})
''').run()
    assert not app.exception
    app.button(key='complex-analysis-load-collection').click().run()
    assert not app.exception
    assert app.session_state['complex-analysis-flow-step'] == 'analyze'
    assert [metric.value for metric in app.metric] == ['3', '1']


def test_native_metric_names_cannot_replace_evidence_identity(records):
    db, group, molecule = records
    item = saved_evidence_page(db, OWNER, [group])['items'][0]
    item['metrics'].update({'Method': 'other engine', 'Target': 'unrelated', 'Molecule / design': 'unrelated'})
    row = evidence_row(item)
    assert row['Method'] in {'MD · MM/GBSA', 'MD · MM/PBSA'}
    assert row['Target'] == 'target-1'
    assert row['Molecule / design'] == 'compound'
    assert item['metrics']['Target'] == 'unrelated'  # original evidence is unchanged


def test_ligand_design_work_keeps_its_declared_activity():
    from mn_molecular_design_studio.complex_analysis_collection_ui import _analysis_activity
    assert _analysis_activity({'origin': {'scientific_activity': 'ligand_design'}}) == 'Ligand design'
    assert _analysis_activity({'origin': {'scientific_activity': 'prediction'}}) == 'Folding and prediction'


def test_summary_describes_work_instead_of_treating_all_records_as_designs():
    from mn_molecular_design_studio.complex_analysis_evidence_ui import evidence_summary_cards
    cards = evidence_summary_cards({'molecule_count': 1, 'structure_count': 0},
        {'md.mmgbsa.gb': 23, 'md.mmgbsa.pb': 23, 'md.analysis': 18, 'md.production': 15}, 1)
    assert dict(cards) == {'MM/GBSA records': 23, 'MM/PBSA records': 23,
                           'trajectory analysis records': 18, 'production information records': 15}
    assert all(value != 79 for _, value in cards)
    generic = evidence_summary_cards({'molecule_count': 3, 'structure_count': 16},
                                    {'vina': 10, 'boltz2': 6}, 1)
    assert dict(generic) == {'Molecules / designs': 3, 'Saved structures': 16,
                             'Evidence types': 2, 'Selected collections': 1}
