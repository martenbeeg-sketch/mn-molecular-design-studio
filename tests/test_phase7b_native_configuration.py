"""PD-L1 native graph configuration contracts; no model or GPU execution."""
import copy
import json

import pytest
from test_api_integrations import api_data

from mn_molecular_design_studio import candidate_pipeline as cp, pipeline_graph as pg
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request


def request(graph):
    graph = pg.validate_graph(graph)['graph']
    node = graph['nodes'][1]
    target = dict(target_pdb='/retained/PD-L1.cif', target_chains=['A'], hotspots='', coldspots='')
    return node_request(node['data']['kind'], node['data']['params'], {'target': target}, 'PD-L1 native contract')


@pytest.mark.parametrize('mode', ['native_complete', 'native_comparative'])
def test_native_configuration_survives_saved_graph_and_whole_invocation(mode, tmp_path, monkeypatch):
    graph = cp.new_pipeline(['bindcraft2'], mode=mode)
    overrides = {'screen_steps': 1, 'refine_steps': 1, 'anneal_steps': 1, 'harden_steps': 1,
                 'mutate_steps': 1, 'sequence_candidates': 1, 'validation_models': ['model_1_ptm']}
    params = graph['nodes'][1]['data']['params']
    params.update(num_designs=1, max_trajectories=2, workers_per_gpu='1',
                  engine_overrides_json=json.dumps(overrides), campaign_seed=20261027)
    monkeypatch.setattr(pg, 'app_home', lambda: tmp_path)
    pg.save_definition(graph, 'owner', 'project', version=True)
    reopened = pg.list_definitions('owner', 'project')[0]['graph']
    p = request(reopened)['parameters']
    assert p['number_of_final_designs'] == 1 and p['max_trajectories'] == 2
    assert p['engine_overrides'] == overrides and p['workers_per_gpu'] == '1'
    assert p['campaign_seed'] == 20261027
    assert 'trajectory_only' not in p  # Complete native adapter keeps its False default.
    assert reopened['nodes'][1]['data']['execution']['mode'] == 'native'
    assert pg.list_versions('owner', 'project', graph['id'])


def test_old_graph_keeps_existing_native_limits_and_defaults():
    graph = cp.new_pipeline(['bindcraft2'], mode='native_complete')
    graph['nodes'][1]['data']['params']['num_designs'] = 7
    p = request(graph)['parameters']
    assert p['max_trajectories'] == p['number_of_final_designs'] == 7
    assert 'workers_per_gpu' not in p and 'engine_overrides' not in p


@pytest.mark.parametrize('settings', ['[]', 'null', '{', '{"x": NaN}', '{"x": 1e999}',
                                    '{"max_trajectories": 100}', '{"trajectory_only": true}',
                                    '{"number_of_final_designs": 50}', '{"workers_per_gpu": 4}'])
def test_invalid_or_budget_changing_native_settings_fail_before_submission(settings):
    graph = cp.new_pipeline(['bindcraft2'], mode='native_complete')
    graph['nodes'][1]['data']['params']['engine_overrides_json'] = settings
    with pytest.raises(ValueError):
        pg.validate_graph(graph)


def test_native_settings_change_scientific_cache_identity():
    original = cp.new_pipeline(['bindcraft2'], mode='native_complete')['nodes'][1]['data']['params']
    changed = copy.deepcopy(original)
    changed['engine_overrides_json'] = '{"screen_steps": 1}'
    assert pg.fingerprint('bindcraft2', original, []) != pg.fingerprint('bindcraft2', changed, [])


def test_protenix_compatibility_probe_matches_native_launcher(monkeypatch):
    from types import SimpleNamespace
    from mn_molecular_design_studio.cluster import gpu_runtime as gr
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(stdout='{"cuda": "12.8", "targets": ["sm_89", "sm_120"]}')
    monkeypatch.setattr(gr.subprocess, 'run', run)
    gr._probe.cache_clear()
    try:
        profile = gr._probe('sha256:protenix-fixture', True)
        assert commands[0][:7] == ['docker', 'run', '--rm', '--network', 'none', '--entrypoint', '/usr/bin/env']
        assert commands[0][8:11] == ['-u', 'PROTENIX_OPT_MODE', 'python']
        assert '--gpus' not in commands[0]
        assert gr.compatible({'compute_capability': '12.0', 'driver_cuda_api': 13020}, {'protenix': profile})
        assert not gr.compatible({'compute_capability': '12.0', 'driver_cuda_api': 12000}, {'protenix': profile})
        gr._probe('sha256:other-fixture')
        assert commands[1][6] == 'python' and '-u' not in commands[1]
    finally:
        gr._probe.cache_clear()


def test_native_tile_exposes_bounded_controls_in_existing_form():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio import candidate_pipeline as cp
from mn_molecular_design_studio.pipeline_node_ui import render_node_form
graph = cp.new_pipeline(['bindcraft2'], mode='native_complete')
values, label, action, replacement = render_node_form(graph['nodes'][1], graph,
    key='native-settings', render_input=lambda: {}, owner_id='owner')
if action == 'apply':
    st.session_state['applied'] = values
''').run(timeout=30)
    assert not app.exception
    app.number_input(key='native-settings-max_trajectories').set_value(1)
    app.selectbox(key='native-settings-workers_per_gpu').select('1')
    app.text_area(key='native-settings-engine_overrides_json').set_value('{"screen_steps": 1}')
    app.button(key='native-settings-apply').click().run(timeout=30)
    assert not app.exception
    assert app.session_state['applied']['max_trajectories'] == 1
    assert app.session_state['applied']['workers_per_gpu'] == '1'
    assert json.loads(app.session_state['applied']['engine_overrides_json']) == {'screen_steps': 1}


@pytest.mark.parametrize('development', [False, True])
def test_native_graph_api_submission_preserves_long_settings_and_retry_identity(api_data, monkeypatch, development):
    """The real Phase 7B failure occurred after validation, during submission."""
    import shutil
    from pathlib import Path
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue, jobs
    from mn_molecular_design_studio import auth
    if development:
        monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE', 'development')
    target = api_data['root'] / 'PD-L1.pdb'
    shutil.copyfile(Path(__file__).parent / 'fixtures/pdl1_p4_ca.pdb', target)
    graph = cp.new_pipeline(['bindcraft2'], mode='native_complete')
    graph['nodes'][0]['data']['params'].update(target_pdb=str(target), target_chains='A')
    overrides = {f'{stage}_steps': 1 for stage in ('screen', 'refine', 'anneal', 'harden', 'mutate')}
    overrides.update(design_models=['model_1_multimer_v3'], validation_models=['model_1_ptm'],
                     sequence_candidates=1, kept_sequences=1, enough_passing_sequences=1,
                     save_failed_trajectories=True, save_design_animations=False, autotune=False)
    content = json.dumps(overrides)
    assert len(content) > 255
    graph['nodes'][1]['data']['params'].update(num_designs=1, max_trajectories=1,
        workers_per_gpu='1', engine_overrides_json=content)
    launches = []
    monkeypatch.setattr(workflow_queue, 'spawn_worker_for_run', lambda run: launches.append(run))
    payload = dict(project_id=api_data['project']['id'], graph=graph, reuse_completed=False,
                   execution={'placement': 'local'}, gpu_device='all')
    headers = {'Idempotency-Key': 'phase7b-inline-native-fixture'}
    response = api_data['client'].post('/api/v1/pipelines/run', json=payload, headers=headers)
    assert response.status_code == 202, response.text
    assert len(launches) == 1
    metadata = jobs.read_json(launches[0] / 'metadata.json')
    assert metadata['studio_owner_id'] == (auth.development_user()['id'] if development else api_data['admin']['id'])
    assert metadata['studio_project_id'] == api_data['project']['id']
    assert metadata['studio_development_mode'] is development
    restored = jobs.read_json(launches[0] / 'worker_request.json')['kwargs']['graph']
    assert restored['nodes'][1]['data']['params']['engine_overrides_json'] == content
    replay = api_data['client'].post('/api/v1/pipelines/run', json=payload, headers=headers)
    assert replay.status_code == 202 and replay.json() == response.json()
    assert len(launches) == 1  # One controller queue entry, no scientific worker execution.
