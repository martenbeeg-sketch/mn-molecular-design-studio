"""Results-only Evaluation federation uses saved records, never Database promotion."""
import copy
import hashlib
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import run_results
from mn_molecular_design_studio.complex_analysis_results_ui import _linked_measurements, board_entries, _remote_board_options
from mn_molecular_design_studio.evaluations import results_evaluation_descriptor

FIXTURE = Path(__file__).parent / 'fixtures' / 'remote_biotite_results.json'
USER = {'id': 'test-user', 'role': 'user', '_cluster_project_ids': ['test-project']}


def evidence():
    return json.loads(FIXTURE.read_text())


def configure_transport(tmp_path, monkeypatch, saved):
    from mn_molecular_design_studio.cluster import config, transport
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(tmp_path / 'studio'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_APP_HOME', str(tmp_path / 'studio'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR', str(tmp_path / 'studio/runs'))
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    config.save({'mode': 'peer', 'coordinator_url': 'http://127.0.0.1:1'})
    state = {'allowed': True, 'actions': 0, 'reads': 0, 'pending': False}
    jobs = {kind: saved[kind]['catalog'] for kind in ('source', 'analysis')}

    def request(path, payload=None, **kwargs):
        assert kwargs['token'] == 'test-session'
        for kind, job in jobs.items():
            base = '/api/v1/cluster/jobs/' + job['cluster_job_id']
            if path == base:
                if not state['allowed']:
                    raise PermissionError('This source project is no longer accessible.')
                return job
            if path == base + '/actions':
                assert payload == {'action': 'inspect'}
                state['actions'] += 1
                return {'command_id': kind}
            if path == '/api/v1/cluster/commands/' + kind:
                state['reads'] += 1
                return {'status': 'pending'} if state['pending'] else {'status': 'completed', 'response': {
                    'job': saved[kind]['job'], 'metadata': saved[kind]['metadata'],
                    'results': {'candidates': [saved[kind]['candidate']]}}}
        raise AssertionError(path)
    monkeypatch.setattr(transport, 'request', request)
    return jobs, state


def test_saved_pc004033_result_is_joined_through_the_federated_results_reader(tmp_path, monkeypatch):
    saved = evidence()
    jobs, state = configure_transport(tmp_path, monkeypatch, saved)
    source = run_results.read_source(jobs['source'], USER, token='test-session')
    selected = {'candidate': source['candidates'][0], 'digest': saved['digest'], 'label': 'Original design'}
    def load(job):
        return run_results.read_source(job, USER, token='test-session')
    rows = _linked_measurements(USER, jobs['source'], [selected], jobs=list(jobs.values()), load=load)
    again = _linked_measurements(USER, jobs['source'], [selected], jobs=[jobs['analysis'], jobs['analysis']], load=load)
    assert len(rows) == len(again) == 1
    assert rows[0]['evaluation_id'] == again[0]['evaluation_id']
    assert state['actions'] == 2  # One metadata request for each source, reused thereafter.
    result = rows[0]
    assert result['metrics']['qc_interface_contact_pairs'] == 430
    assert result['metrics']['qc_interface_residue_pairs'] == 63
    assert result['metrics']['qc_interface_salt_bridges'] == 1
    assert result['source_candidate_id'] == 'rfd4_000050_protein_binder_rep9'
    assert result['source_structure_digest'] == saved['digest']
    assert result['source_model_id'] == ''  # An artifact filename is not a scientific model ID.
    assert result['parameters']['contact_distance'] == 5.0
    assert result['binder_chains'] == ['B'] and result['target_chains'] == ['A']
    assert result['source_workstation'] == 'PC004033'
    assert result['provenance']['evaluation']['scope'] == 'results'
    assert not list(tmp_path.rglob('molecular_library.sqlite3'))


def test_cached_inspect_is_not_an_authorization_grant(tmp_path, monkeypatch):
    saved = evidence()
    jobs, state = configure_transport(tmp_path, monkeypatch, saved)
    assert run_results.read_source(jobs['analysis'], USER, token='test-session')['candidate_count'] == 1
    reads = state['reads']
    state['allowed'] = False
    denied = run_results.read_source(jobs['analysis'], USER, token='test-session')
    assert denied['candidates'] == [] and not denied['complete']
    assert 'no longer accessible' in denied['source_errors'][0]['error']
    assert state['reads'] == reads


def test_coordinator_rejects_a_results_source_from_an_unrelated_project(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.cluster.store import Store
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(tmp_path / 'studio'))
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    config.save({'mode': 'coordinator', 'coordinator_url': 'http://127.0.0.1:1'})
    store = Store()
    node = store.enroll('source-workstation')
    source = {**evidence()['source']['catalog'], 'node_id': node['node_id']}
    source['cluster_job_id'] = config.job_id(node['node_id'], source['local_reference'])
    store.heartbeat(node['node_id'], {'accept_jobs': False}, [source])
    outsider = {'id': 'another-user', 'role': 'user', '_cluster_project_ids': []}
    result = run_results.read_source(source, outsider)
    assert result['candidate_count'] == 0 and result['source_errors']
    with store.connect() as connection:
        assert connection.execute('select count(*) from commands').fetchone()[0] == 0


def test_pending_results_do_not_become_zero_measurements(tmp_path, monkeypatch):
    saved = evidence()
    jobs, state = configure_transport(tmp_path, monkeypatch, saved)
    state['pending'] = True
    result = run_results.read_source(jobs['analysis'], USER, token='test-session')
    assert result['candidates'] == [] and result['pending_sources'] and not result['complete']


@pytest.mark.parametrize('field,value', [('node_id', 'different-workstation'), ('upstream_run_id', 'unrelated-run')])
def test_other_source_identity_is_not_joined(tmp_path, monkeypatch, field, value):
    saved = evidence()
    jobs, _ = configure_transport(tmp_path, monkeypatch, saved)
    other = {**jobs['analysis'], field: value}
    assert _linked_measurements(USER, jobs['source'], [{'candidate': saved['source']['candidate'],
        'digest': saved['digest'], 'label': 'Design'}], jobs=[other], load=lambda _: pytest.fail('Unrelated job must not be read')) == []


def test_grouped_graph_handoff_preserves_original_remote_source(tmp_path, monkeypatch):
    saved = evidence()
    for kind in ('source', 'analysis'):
        saved[kind]['catalog']['logical_run_id'] = 'existing-graph-run'
    saved['analysis']['catalog']['upstream_run_id'] = 'intermediate-handoff'
    jobs, _ = configure_transport(tmp_path, monkeypatch, saved)
    selected = {'candidate': saved['source']['candidate'], 'digest': saved['digest'], 'label': 'Design'}
    rows = _linked_measurements(USER, jobs['source'], [selected], jobs=[jobs['analysis']],
        load=lambda row: run_results.read_source(row, USER, token='test-session'))
    assert len(rows) == 1 and rows[0]['metrics']['qc_interface_contact_pairs'] == 430
    assert rows[0]['source_job'].endswith(':design/20261102-081701-97f7d42e')


def test_descriptor_is_stable_across_sync_paths_and_preserves_model_zero():
    saved = evidence(); candidate = saved['analysis']['candidate']; job = saved['analysis']['catalog']
    descriptor = results_evaluation_descriptor(candidate, job)
    moved = results_evaluation_descriptor(candidate, {**job, 'run_dir': '/other/workstation/cache', 'node_name': '4033'})
    assert descriptor['id'] == moved['id']
    runtime_only = copy.deepcopy(candidate)
    runtime_only['raw_metadata']['analysis_parameters']['workers'] = 12
    runtime_only['binder_chains'].reverse()
    assert results_evaluation_descriptor(runtime_only, job)['id'] == descriptor['id']
    assert descriptor['id'].startswith('results-evaluation-')
    assert descriptor['observation_id'] == '' and descriptor['model_id'] == ''
    model = copy.deepcopy(candidate); model['model_index'] = 0; model['metrics']['zero_measurement'] = 0
    entry = board_entries(model, label='Design', run_id=job['run_id'], digest=saved['digest'], job=job)[0]
    assert entry['source_model_id'] == '0' and entry['metrics']['zero_measurement'] == 0
    assert entry['metrics'].get('unavailable_measurement') is None
    assert entry['evaluation_id'] != descriptor['id']
    changed = copy.deepcopy(candidate); changed['raw_metadata']['analysis_parameters']['contact_distance'] = 4.0
    assert results_evaluation_descriptor(changed, job)['id'] != descriptor['id']


def test_remote_options_reject_a_structure_outside_its_source_job():
    saved = evidence(); source = saved['source']; row = copy.deepcopy(source['candidate'])
    row['complex_pdb'] = '/external/private.cif'
    collection = {'candidates': [row], 'source_details': {'cluster:source-job': source}}
    assert _remote_board_options(source['catalog'], collection) == []


def test_remote_board_renders_saved_metrics_and_preserves_inspector_route(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_dashboard, complex_analysis_ui
    from mn_molecular_design_studio.cluster import ui
    saved = evidence()
    content = b'HEADER    isolated structure artifact\nEND\n'
    digest = hashlib.sha256(content).hexdigest()
    saved['digest'] = digest
    saved['analysis']['candidate']['raw_metadata']['structure_digest'] = digest
    jobs, _ = configure_transport(tmp_path, monkeypatch, saved)
    run = {'id': 'test-run', 'jobs': list(jobs.values())}
    monkeypatch.setattr(complex_analysis_dashboard, '_run_catalog', lambda *a, **k: [run])
    monkeypatch.setattr(ui, '_fetch_remote_artifact_bytes', lambda *a, **k: content)
    opened = []
    monkeypatch.setattr(complex_analysis_ui, 'open_remote_result_complex', lambda *a, **k: opened.append((a,k)))
    script = f'''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_results_ui import render_results_board
st.session_state['studio-session-token'] = 'test-session'
st.session_state.setdefault('results-analysis-source', {{'reference': 'cluster:source-job', 'candidate_ids': ['rfd4_000050_protein_binder_rep9']}})
render_results_board({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=20).run()
    assert not app.exception
    boards = [table.value for table in app.dataframe if 'Evaluation ID' in table.value.columns]
    assert len(boards) == 1
    board = boards[0]; analysis = board[board['Method'] == 'biotite_interface']
    assert len(analysis) == 1
    assert analysis.iloc[0]['qc_interface_contact_pairs'] == 430
    assert analysis.iloc[0]['Source candidate'] == 'rfd4_000050_protein_binder_rep9'
    button = next(button for button in app.button if button.label == 'Open in Scientific Inspector')
    button.click().run()
    assert not app.exception and len(opened) == 1
    assert opened[0][0][0]['cluster_job_id'] == 'source-job'
    assert opened[0][1]['structure_sha256'] == digest
    navigated = []
    monkeypatch.setattr(ui, '_navigate_to_job_details', lambda *a, **k: navigated.append((a,k)))
    next(button for button in app.button if button.label == 'Open source Results').click().run()
    assert navigated[0][0] == ('cluster:source-job',)
    assert app.session_state['cluster-candidate-default-source-job'] == 'rfd4_000050_protein_binder_rep9'


def test_production_api_requires_real_session_and_rejects_development_bearer(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio.cluster import config
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(tmp_path / 'studio'))
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE', 'production')
    config.save({'mode': 'disabled'})
    from mn_molecular_design_studio import api
    account = {**USER, 'username': 'test-researcher', 'display_name': 'Researcher', 'must_change_password': False}
    monkeypatch.setattr(api, 'resolve_session', lambda token: account if token == 'valid-session' else None)
    client = TestClient(api.app)
    assert client.get('/api/v1/auth/me').status_code == 401
    assert client.get('/api/v1/auth/me', headers={'Authorization': 'Bearer development-mode'}).status_code == 401
    assert client.get('/api/v1/cluster/jobs').status_code == 401
    signed_in = client.get('/api/v1/auth/me', headers={'Authorization': 'Bearer valid-session'})
    assert signed_in.status_code == 200
    assert signed_in.json()['id'] == USER['id']
    assert signed_in.json()['project_ids'] == ['test-project']


@pytest.mark.parametrize('model_id', [None, 0])
def test_remote_inspector_recognizes_saved_measurement_and_does_not_offer_duplicate_biotite(tmp_path, monkeypatch, model_id):
    from mn_molecular_design_studio import complex_analysis_dashboard, complex_analysis_ui
    saved = evidence()
    if model_id is not None:
        for kind in ('source', 'analysis'):
            saved[kind]['candidate']['model_index'] = model_id
    content = b'HEADER    saved complex artifact\nEND\n'
    saved['analysis']['candidate']['raw_metadata']['structure_digest'] = hashlib.sha256(content).hexdigest()
    jobs, _ = configure_transport(tmp_path, monkeypatch, saved)
    monkeypatch.setattr(complex_analysis_dashboard, '_run_catalog', lambda *a, **k: [{'id': 'run', 'jobs': list(jobs.values())}])
    monkeypatch.setattr(complex_analysis_ui, '_authorized_remote_source_project', lambda _: {'id': 'test-project', 'name': 'Project'})
    stage = tmp_path / 'stage'
    candidate_file = stage / 'artifacts/normalized_candidates/candidates.jsonl'
    candidate_file.parent.mkdir(parents=True)
    candidate_file.write_text(json.dumps(saved['source']['candidate']) + '\n')
    structure = stage / 'complex.cif'
    structure.write_bytes(content)
    selected = {'path': str(structure), 'relative_path': saved['source']['candidate']['complex_pdb'],
                'candidate': saved['source']['candidate']}
    job = {**jobs['source'], 'run_dir': str(stage)}
    script = f'''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import _render_biotite_workspace
st.session_state['studio-session-token'] = 'test-session'
_render_biotite_workspace({USER!r}, {job!r}, {selected!r},
    [{{'type': 'protein', 'chain': 'A'}}, {{'type': 'protein', 'chain': 'B'}}],
    workflow_schema=lambda _: {{'parameters': []}}, project_selector=None, submit=None, viewer_key='saved')
'''
    app = AppTest.from_string(script, default_timeout=20).run()
    assert not app.exception
    assert any('Saved Biotite Interface measurements' in message.value for message in app.success)
    assert not any(str(button.key).startswith('complex-analysis-run-biotite-') for button in app.button)
