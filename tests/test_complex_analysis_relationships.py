"""Read-only scientific relationship navigation, independent of producing engine."""
import json
import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_catalog import (
    saved_relationship_page, saved_evidence_page, saved_evidence_summary,
    saved_execution_groups, relationship_structure_ids,
)

OWNER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def studies(tmp_path):
    db = MolecularDatabase(tmp_path / 'db')
    imported = db.import_records(OWNER, molecule_imports.parse_smiles('CCO shared-ligand'), name='Studies')
    group, molecule = imported['group_id'], imported['molecule_ids'][0]
    for study, target, prep in [('study-a', 'alpha', 'prep-a'), ('study-b', 'alpha', 'prep-a'),
                                ('study-c', 'beta', 'prep-b')]:
        for replica in range(3):
            context = {'study_id': study, 'target_id': target, 'replica_index': replica,
                       'prepared_target_run_id': prep, 'source_production_run_id': f'{study}-production-{replica}',
                       'target_receptor': 'Same display name', '_replica_node': 'pc004033'}
            for method in ['md.production', 'md.analysis', 'md.mmgbsa.gb', 'md.mmgbsa.pb']:
                db.add_observation(group, molecule, method, 'reused-normalized-run', {'value': replica}, context=context)
        if study == 'study-a':
            db.add_observation(group, molecule, 'md.production', 'reused-normalized-run', {'value': 0}, context={
                **context, 'replica_index': 0, 'source_production_run_id': 'continued-production',
                'continuation_of_run_id': 'study-a-production-0'})
            db.add_observation(group, molecule, 'md.mmgbsa.gb', 'reused-normalized-run', {'value': -999}, context={
                **context, 'replica_index': 0, 'source_production_run_id': 'study-a-production-0',
                'endpoint_superseded_by_run_id': 'replacement-endpoint'})
    return db, group, molecule


def test_one_ligand_two_preparations_three_studies_are_distinct(studies):
    db, group, molecule = studies
    page = saved_relationship_page(db, OWNER, [group])
    assert len(page['items']) == 3
    assert page['totals']['molecules'] == 1
    assert page['totals']['target_contexts'] == 2
    assert page['totals']['studies'] == 3
    assert {row['study_id'] for row in page['items']} == {'study-a', 'study-b', 'study-c'}
    assert all(row['molecule_id'] == molecule and row['replica_count'] == 3 for row in page['items'])
    assert len({row['target_name'] for row in page['items']}) == 1  # labels do not merge targets
    assert saved_relationship_page(db, OWNER, [group], query='beta')['totals']['studies'] == 1
    assert saved_relationship_page(db, OWNER, [group], limit=2)['next_offset'] == 2
    assert len(saved_relationship_page(db, OWNER, [group], offset=2, limit=2)['items']) == 1


def test_scope_replica_zero_history_and_continuation(studies):
    db, group, _ = studies
    scope = saved_relationship_page(db, OWNER, [group])['items'][0]['relationship']
    current = saved_evidence_page(db, OWNER, [group], relationship=scope, replica=0,
                                  namespace='md.mmgbsa.gb', include_history=False)
    assert len(current['items']) == 1
    assert current['items'][0]['metrics']['value'] == 0
    assert current['items'][0]['context']['study_id'] == 'study-a'
    history = saved_evidence_page(db, OWNER, [group], relationship=scope, replica=0, namespace='md.mmgbsa.gb')
    assert len(history['items']) == 2
    summary = saved_evidence_summary(db, OWNER, [group], relationship=scope, include_history=False)
    assert summary[group]['methods']['md.mmgbsa.gb'] == 3
    execution = saved_execution_groups(db, OWNER, [group], scope)
    continuation = next(item for item in execution if item['production_run'] == 'continued-production')
    assert continuation['replica'] == 0 and continuation['recorded_parents'] == 'study-a-production-0'
    assert sum(item['historical_count'] for item in execution) == 1
    # Nothing was deleted or force-finalized by filtering the view.
    assert saved_evidence_page(db, OWNER, [group])['total'] == 38


def test_acl_and_relationship_revalidation(studies):
    db, group, _ = studies
    scope = saved_relationship_page(db, OWNER, [group])['items'][0]['relationship']
    for call in [saved_relationship_page]:
        with pytest.raises(PermissionError): call(db, OTHER, [group])
    for call in [saved_execution_groups, relationship_structure_ids]:
        with pytest.raises(PermissionError): call(db, OTHER, [group], scope)
    with pytest.raises(ValueError):
        saved_evidence_page(db, OWNER, [group], relationship={**scope, 'group_id': 'unrelated'})
    with pytest.raises(ValueError):
        saved_relationship_page(db, OWNER, [group], limit=10000)


def test_models_require_explicit_study_link_not_shared_molecule_or_run(studies):
    db, group, molecule = studies
    page = saved_relationship_page(db, OWNER, [group])
    a, b = page['items'][0]['relationship'], page['items'][1]['relationship']
    structure = db.attach_structure(group, molecule, {'name': 'linked.cif', 'content': b'data_test\n'},
        context={'run_id': 'reused-normalized-run', 'target_id': 'alpha', 'prepared_target_run_id': 'prep-a'})
    assert relationship_structure_ids(db, OWNER, [group], a) == set()
    db.add_observation(group, molecule, 'test-model-evidence', 'reused-normalized-run', {'value': 1}, context={
        'study_id': 'study-a', 'target_id': 'alpha', 'prepared_target_run_id': 'prep-a', 'representation_id': structure})
    assert relationship_structure_ids(db, OWNER, [group], a) == {structure}
    assert relationship_structure_ids(db, OWNER, [group], b) == set()


def test_shared_docking_cofolding_context_keeps_models_scores_and_preparations(studies):
    db, _, _ = studies
    imported = db.import_records(OWNER, molecule_imports.parse_smiles('CCN ligand'), name='Docking and cofolding')
    group, molecule = imported['group_id'], imported['molecule_ids'][0]
    for method, model, metrics in [('vina', 'pose-1', {'energy': -7}), ('Boltz-2', 'model-1', {'confidence': 0.9})]:
        db.add_observation(group, molecule, method, 'run-1', metrics,
            context={'target_id': 'receptor', 'prepared_target_run_id': 'prep-1', 'model_id': model})
    db.add_observation(group, molecule, 'vina', 'run-1', {'energy': -8},
        context={'target_id': 'receptor', 'prepared_target_run_id': 'prep-2', 'model_id': 'pose-2'})
    page = saved_relationship_page(db, OWNER, [group])
    assert page['totals']['relationships'] == 2 and page['totals']['studies'] == 0
    scope = next(row['relationship'] for row in page['items'] if row['preparation_id'] == 'prep-1')
    records = saved_evidence_page(db, OWNER, [group], relationship=scope)['items']
    assert {row['namespace'] for row in records} == {'vina', 'Boltz-2'}
    assert {row['evaluation']['model_id'] for row in records} == {'pose-1', 'model-1'}
    assert len({row['evaluation']['id'] for row in records}) == 2


def test_ui_switches_study_without_loading_coordinates(studies, monkeypatch):
    from streamlit.testing.v1 import AppTest
    db, group, _ = studies
    monkeypatch.setattr('mn_molecular_design_studio.molecular_database.MolecularDatabase', lambda: db)
    monkeypatch.setattr(db, 'structure_file', lambda *args: pytest.fail('No coordinates for relationship navigation'))
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_evidence_ui import render_saved_evidence
import streamlit as st
selection = st.session_state.setdefault('complex-analysis-evidence-selection', {{'group_ids': [{group!r}], 'label': 'Study cohort'}})
render_saved_evidence(MolecularDatabase(), {OWNER!r}, selection)
''').run()
    assert not app.exception
    assert [item.value for item in app.metric] == ['1', '2', '3']
    app.session_state['complex-analysis-workflow-samples'] = ['model-from-another-target']
    app.session_state['complex-analysis-assessment-open'] = True
    next(item for item in app.selectbox if item.label == 'Open molecular system / study').set_value(2).run()
    assert not app.exception
    assert all(value.startswith('study-c') for value in app.dataframe[-1].value['Study / replica'])
    assert app.session_state['complex-analysis-evidence-selection']['relationship']['target_ref'] == 'beta'
    assert 'complex-analysis-workflow-samples' not in app.session_state
    assert app.session_state['complex-analysis-assessment-open'] is False


def test_missing_references_do_not_merge_independent_source_memberships(studies):
    db, _, _ = studies
    ids = []
    for _ in range(2):
        imported = db.import_records(OWNER, molecule_imports.parse_sequences('>same-protein\nACDEFGHIK'), name='Same source label')
        group, molecule = imported['group_id'], imported['molecule_ids'][0]
        ids.append(group)
        db.add_observation(group, molecule, 'prediction', 'same-run', {'confidence': 0.9}, context={'model_id': 'model-1'})
    page = saved_relationship_page(db, OWNER, ids)
    assert page['totals']['molecules'] == 1
    assert page['totals']['relationships'] == 2
    assert page['totals']['target_contexts'] == 0
    assert all(row['study_id'] == '' and row['target_ref'] == '' for row in page['items'])
    assert len({row['group_id'] for row in page['items']}) == 2


def test_explicit_target_context_can_link_model_without_a_study(studies):
    db, group, molecule = studies
    model = db.attach_structure(group, molecule, {'name': 'other.cif', 'content': b'data_other\n'},
        context={'target_id': 'receptor', 'prepared_target_run_id': 'prep-1'})
    db.add_observation(group, molecule, 'prediction', 'run', {'confidence': 0.8},
        context={'target_id': 'receptor', 'prepared_target_run_id': 'prep-1'})
    row = next(row for row in saved_relationship_page(db, OWNER, [group])['items'] if row['study_id'] == '')
    assert relationship_structure_ids(db, OWNER, [group], row['relationship']) == {model}


def test_false_supersession_flag_is_not_a_recorded_replacement(studies):
    db, group, molecule = studies
    scope = saved_relationship_page(db, OWNER, [group])['items'][0]['relationship']
    db.add_observation(group, molecule, 'new-method', 'run', {'value': 1}, context={
        'study_id': scope['study_id'], 'target_id': scope['target_ref'], 'prepared_target_run_id': scope['preparation_id'],
        'endpoint_superseded_by_run_id': False})
    page = saved_evidence_page(db, OWNER, [group], relationship=scope, namespace='new-method', include_history=False)
    assert page['total'] == 1
