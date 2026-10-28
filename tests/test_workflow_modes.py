"""Phase 7A.2 normal graph/API/UI integration; PD-L1 routing fixtures only."""
import copy
import json
import shutil
from pathlib import Path

import pytest
from mn_molecular_design_studio import candidate_pipeline as cp, candidate_sets as cs, pipeline_graph as pg
from mn_molecular_design_studio import workflow_modes as wm
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as gp
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
import test_candidate_sets as base


@pytest.fixture
def pdl1_harness(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import config
    monkeypatch.setattr(config, "enabled", lambda: False)
    fixture = Path(__file__).parent / 'fixtures/pdl1_p4_ca.pdb'
    monkeypatch.setattr(base, 'structure', lambda path, **kwargs: str(shutil.copyfile(fixture, path)))
    return base.collection_harness.__wrapped__(tmp_path, monkeypatch)


def configured(h, mode, count, backbone=(), rejected=()):
    graph = cp.new_pipeline(generators=['rfd4'] if mode == 'modular_design' else ['bindcraft2'], mode=mode)
    graph['nodes'][0]['data']['params'].update(target_pdb=h.target, target_chains='A')
    source = h.source('source', 'rfd4_proteina' if mode == 'modular_design' else 'bindcraft2', count,
                      backbone=backbone, rejected=rejected)
    graph['nodes'][1]['data']['retained_run_dir'] = str(source)
    write_json(source / 'result.json', dict(success=True, metrics={'native_rejected_attempts': 2}))
    for node in graph['nodes']:
        if node['id'].startswith('predictor-'):
            node['data']['enabled'] = True
        if node['id'] == 'ranking':
            node['data']['params']['ranking_rules'] = [{'field': 'prediction/predictor-0/iptm', 'direction': 'higher'}]
    return graph


@pytest.mark.parametrize('mode', wm.MODES)
def test_modes_are_editable_serializable_graphs(mode, tmp_path, monkeypatch):
    monkeypatch.setattr(pg, 'app_home', lambda: tmp_path)
    original = pg.new_graph(mode)
    graph = pg.validate_graph(json.loads(json.dumps(original)))['graph']
    assert graph['workflow_mode'] == wm.describe(graph)['effective_mode'] == mode
    pg.save_definition(graph, 'owner', 'project', version=True)
    assert pg.list_definitions('owner', 'project')[0]['graph'] == graph
    assert pg.list_versions('owner', 'project', graph['id'])
    # Preset metadata never constrains graph editing or rewrites native settings.
    next(n for n in graph['nodes'] if n['id'] == 'generator-0')['data']['params']['num_designs'] = 7
    assert pg.validate_graph(graph)['graph']['nodes'][1]['data']['params']['num_designs'] == 7
    if mode == 'native_complete':
        assert all(pg.NODE_CATALOG[n['data']['kind']]['category'] not in {'Structure prediction', 'Sequence design'} for n in graph['nodes'])
    if mode == 'native_comparative':
        assert not any(pg.NODE_CATALOG[n['data']['kind']]['category'] == 'Sequence design' for n in graph['nodes'])


@pytest.mark.parametrize('mode', wm.MODES)
@pytest.mark.parametrize('size', [0, 1, 12])
def test_normal_execution_preserves_native_outcomes_and_runs_downstream(pdl1_harness, mode, size):
    h = pdl1_harness
    graph = configured(h, mode, size, rejected=[0])
    run = gp.run_graph_pipeline(graph, 'PD-L1 fixture modes', 'owner', 'project', reuse_completed=False)
    result = read_json(run / 'result.json')
    assert result['success'] is True
    rows = result['outputs']['candidates']
    assert len(rows) == size
    summary = result['outputs']['workflow_summary']
    assert summary['effective_mode'] == mode and summary['candidate_count'] == size
    assert summary['stage_counts']['collection']['native_rejected_attempts'] == 2
    assert len(cs.target_identity(dict(target_pdb=h.target, target_chains='A'))['chains']['A']['sequence']) == 115
    if mode == 'native_complete':
        assert h.submissions == []
        assert all('selection' not in row['raw_metadata'] for row in rows)
        assert all(wm.native_ranking(row)['status'] == 'unavailable' for row in rows)
    else:
        assert all(pg.NODE_CATALOG[s[1]['data']['kind']]['category'] == 'Structure prediction' for s in h.submissions)
        assert len(h.submissions) == 2 * ((max(0, size - 1) + 4) // 5)
        assert all(len(cs.metadata(r)['predictions']) == (0 if cs.metadata(r)['native_filter'] == 'failed' else 2) for r in rows)
    from mn_molecular_design_studio.result_presentation import candidate_state_rows
    if rows:
        assert candidate_state_rows(rows)[0]['Workflow mode'] == wm.MODES[mode]
        assert all(cs.metadata(r)['native_measurements']['native_score'] == r['metrics']['native_score'] for r in rows)
    from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import build_run_zip
    from zipfile import ZipFile
    from io import BytesIO
    with ZipFile(BytesIO(build_run_zip(run))) as archive:
        summary = json.loads(archive.read('artifacts/candidate_set/workflow_summary.json'))
        assert summary['effective_mode'] == mode
        assert len(json.loads(archive.read('artifacts/candidate_set/native_sources.json'))) == 1


def test_modular_backbones_initial_design_then_opt_in_child_redesign(pdl1_harness):
    h = pdl1_harness
    graph = configured(h, 'modular_design', 3, backbone=[0])
    redesign = next(n for n in graph['nodes'] if n['id'] == 'optional-redesign')
    redesign['data']['params'].update(optimization_enabled=True, model_type='soluble_mpnn', sequences=1)
    run = gp.run_graph_pipeline(graph, 'PD-L1 deterministic sequence children', 'owner', 'project', reuse_completed=False)
    result = read_json(run / 'result.json')
    assert result['success'] is True
    rows = result['outputs']['candidates']
    originals = [r for r in rows if not cs.metadata(r)['parent_variant_ids']]
    children = [r for r in rows if cs.metadata(r)['parent_variant_ids']]
    assert len(originals) == 3 and len(children) == 4  # initial + explicit redesign children
    assert any(cs.sequence_state(r) == 'backbone_only' for r in originals)
    assert all(r['binder_sequence'] == cs.metadata(r)['original_sequence'] for r in originals if r['binder_sequence'])
    assert all('native_score' not in r['metrics'] for r in children)
    assert all(cs.metadata(r)['native_filter'] == 'not_evaluated' for r in children)
    requests = [s[0]['parameters'] for s in h.submissions if s[1]['id'] == 'optional-redesign']
    assert requests and all(p['fix_interface_residues'] is True for p in requests)
    assert len({cs.metadata(r)['variant_id'] for r in rows}) == len(rows)
    assert any(s[1]['id'] == 'sequence-preparation' for s in h.submissions)
    assert any(s[1]['id'] == 'optional-redesign' for s in h.submissions)


@pytest.mark.parametrize('mode', ['native_comparative', 'modular_design'])
def test_failed_prediction_branch_keeps_missing_evidence(pdl1_harness, mode):
    h = pdl1_harness
    graph = configured(h, mode, 2)
    h.failures.add('predictor-1')
    submit = gp._submit
    def failed_without_rows(*args, **kwargs):
        child = submit(*args, **kwargs)
        if 'predictor-1' in args[2]['id']:
            (child / 'artifacts/normalized_candidates/candidates.jsonl').write_text('')
        return child
    from unittest.mock import patch
    patcher = patch.object(gp, '_submit', failed_without_rows)
    patcher.start()
    h.patcher = patcher
    run = gp.run_graph_pipeline(graph, 'PD-L1 partial predictor fixture', 'owner', 'project', reuse_completed=False)
    h.patcher.stop()
    result = read_json(run / 'result.json')
    assert result['success'] is True
    assert all(cs.metadata(row)['predictor_matrix']['predictor-1'] != ['completed'] for row in result['outputs']['candidates'])
    assert all('prediction/predictor-1/iptm' not in row['metrics'] for row in result['outputs']['candidates'])


def test_capability_selection_and_legacy_graphs():
    for mode in ('native_complete', 'native_comparative'):
        assert {'bindcraft2', 'boltzgen', 'pxdesign_pipeline', 'genie3_pipeline'} <= set(wm.methods(mode))
        for invalid in ('rfd4', 'rfdiffusion', 'esm_binder', 'bindcraft2_generator'):
            with pytest.raises(ValueError, match='unsupported'):
                cp.new_pipeline([invalid], mode=mode)
    with pytest.raises(ValueError, match='unsupported'):
        cp.new_pipeline(['bindcraft2'], mode='modular_design')
    assert wm.capability('rfd4')['complete_native_workflow'] is False
    assert 'workflow_mode' not in pg.clean_graph(cp.new_pipeline())
    assert wm.native_ranking({'metrics': {'native_rank': float('nan')}})['rank'] is None
    assert wm.native_ranking({'metrics': {'native_rank': 2}})['rank'] == 2


@pytest.mark.parametrize('mode', wm.MODES)
def test_modes_through_setup_ui_and_results(pdl1_harness, mode):
    from streamlit.testing.v1 import AppTest
    h = pdl1_harness
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.scientific_tasks_ui import render_task_setup
graph = render_task_setup('mode-test')
if graph is not None:
    st.session_state['created-mode-graph'] = graph
''').run(timeout=45)
    app.selectbox(key='mode-test-workflow-mode').select(mode).run()
    assert not app.exception
    generator = next(m for m in app.multiselect if m.label == 'Generator operations')
    assert all(k in wm.methods(mode) for k in generator.value)
    app.button(key='mode-test-create-mode').click().run()
    assert not app.exception
    graph = app.session_state['created-mode-graph']
    assert pg.validate_graph(graph)['graph']['workflow_mode'] == mode
    # Render the same function used by normal Results, including empty native results.
    graph = configured(h, mode, 0)
    run = gp.run_graph_pipeline(graph, 'PD-L1 empty Results', 'owner', 'project', reuse_completed=False)
    panel = AppTest.from_string(f'''
from mn_molecular_design_studio.result_presentation import render_workflow_results_summary
render_workflow_results_summary({str(run)!r})
''').run(timeout=45)
    assert not panel.exception
    assert any('Candidates: 0' in caption.value for caption in panel.caption)


@pytest.mark.parametrize('mode', wm.MODES)
def test_modes_through_normal_api_save_validate(pdl1_harness, monkeypatch, mode):
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio import api
    h = pdl1_harness
    graph = configured(h, mode, 2)
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(h.tmp_path / 'cluster'))
    monkeypatch.setattr(api, '_request_payload_and_project', lambda payload, user: ({k:v for k,v in payload.items() if k != 'project_id'}, payload['project_id']))
    monkeypatch.setattr(api, '_validate_user_path_values', lambda *a, **k: None)
    monkeypatch.setattr(api, '_guard_source_job', lambda *a, **k: None)
    monkeypatch.setattr(api, '_resolve_api_job', lambda *a, **k: {'run_dir': graph['nodes'][1]['data']['retained_run_dir']})
    api.app.dependency_overrides[api._current_user] = lambda: {'id': 'f'*32, 'role': 'user'}
    try:
        with TestClient(api.app) as client:
            payload = dict(project_id='project', graph=graph, version=True)
            saved = client.post('/api/v1/pipelines/save', json=payload)
            assert saved.status_code == 200, saved.text
            loaded = client.get('/api/v1/pipelines', params={'project_id':'project'})
            assert loaded.json()[0]['graph']['workflow_mode'] == mode
            response = client.post('/api/v1/pipelines/validate', json=dict(project_id='project', graph=graph))
            assert response.status_code == 200, response.text
            assert response.json()['valid'] is True
            catalog = client.get('/api/v1/pipelines/catalog').json()
            assert catalog['rfd4']['scientific_workflow']['complete_native_workflow'] is False
            assert h.parents == [] and h.submissions == []
    finally:
        api.app.dependency_overrides.clear()


def test_four_gpu_inventory_and_stale_workstations():
    nodes = [dict(id=str(i), name='fixture-worker-'+str(i), availability='online', snapshot={'resources':dict(
        gpus=[{'gpu_id':str(j), 'name':'fixture-gpu'} for j in range(count)],
        free_gpus=[{'gpu_id':str(j)} for j in range(count) if not (i == 0 and j == 0)])}) for i,count in enumerate([2,1,1])]
    rows = wm.gpu_monitor_rows(nodes)
    assert len(rows) == 4
    assert [r['State'] for r in rows].count('available') == 3
    nodes[1]['availability'] = 'offline'
    assert wm.gpu_monitor_rows(nodes)[2]['State'] == 'unknown / workstation offline'


@pytest.mark.parametrize('mode', ['native_comparative', 'modular_design'])
def test_independent_generators_and_prediction_overlap(pdl1_harness, monkeypatch, mode):
    h = pdl1_harness
    generators = ['bindcraft2', 'boltzgen'] if mode == 'native_comparative' else ['rfd4', 'rfdiffusion']
    graph = cp.new_pipeline(generators, mode=mode)
    graph['nodes'][0]['data']['params'].update(target_pdb=h.target, target_chains='A')
    for node in graph['nodes']:
        if pg.NODE_CATALOG[node['data']['kind']]['category'] == 'Structure prediction':
            node['data']['enabled'] = True
    # Ordinary request validation is independently covered; these workers are
    # deterministic scheduler stand-ins and run no containers/models.
    monkeypatch.setattr(cp, 'preflight', lambda *a: None)
    submit = gp._submit
    active, events = {}, []
    second = None
    polls = [0]
    def worker(request, parent, node, **kwargs):
        nonlocal second
        category = pg.NODE_CATALOG[node['data']['kind']]['category']
        if category in {'Design engine', 'Design pipeline'}:
            index = len(active)
            tool = 'bindcraft2' if index == 0 else 'boltzgen'
            child = h.source('native-' + str(index), tool, 12)
            write_json(child / 'metadata.json', dict(status='running', tool=tool))
            active[str(child)] = child
            if index == 1:
                second = child
            events.append('generator-'+str(index))
            return child
        events.append('prediction-while-generating' if second and read_json(second/'metadata.json')['status'] == 'running' else 'downstream')
        return submit(request, parent, node, **kwargs)
    def advance(*args):
        polls[0] += 1
        if len(active) == 2:
            first = next(iter(active.values()))
            write_json(first/'metadata.json', dict(status='completed', tool='bindcraft2'))
            if any(e == 'prediction-while-generating' for e in events):
                write_json(second/'metadata.json', dict(status='completed', tool='boltzgen'))
        assert polls[0] < 20, 'Graph stalled behind a global generation barrier'
    monkeypatch.setattr(gp, '_submit', worker)
    monkeypatch.setattr(gp.time, 'sleep', advance)
    run = gp.run_graph_pipeline(graph, 'PD-L1 overlapping fixture workers', 'owner', 'project', reuse_completed=False)
    result = read_json(run/'result.json')
    assert result['success'] is True
    assert len(result['outputs']['candidates']) == 24
    assert events[:2] == ['generator-0','generator-1']
    assert 'prediction-while-generating' in events
    assert all(cs.metadata(row)['predictor_matrix']['predictor-0'] == ['completed'] for row in result['outputs']['candidates'])
    assert all(cs.metadata(row)['predictor_matrix']['predictor-1'] == ['completed'] for row in result['outputs']['candidates'])
    state = read_json(run/'graph_execution.json')['nodes']
    assert len(state['generator-0']['batches']) == len(state['generator-0-g1']['batches']) == 1
    assert all(n['data']['execution']['placement'] == 'auto' for n in graph['nodes'] if n['data']['kind'] in generators)


def test_four_independent_chunks_can_be_in_flight(pdl1_harness, monkeypatch):
    h = pdl1_harness
    graph = configured(h, 'modular_design', 12)
    submit = gp._submit
    active, maximum = [], [0]
    def worker(*args, **kwargs):
        child = submit(*args, **kwargs)
        write_json(child/'metadata.json', dict(status='running'))
        active.append(child)
        maximum[0] = max(maximum[0],len(active))
        assert len(active) <= 4  # Two branches, each max_in_flight=2.
        return child
    def advance(*args):
        assert active
        # Four independent chunks enter existing scheduler queues before polling
        # releases them. The four device leases are simulated, not GPU evidence.
        for i, child in enumerate(active):
            write_json(child/'metadata.json', dict(status='completed',resource_allocation={'gpu_device':str(i), 'worker':'fixture'}))
        active.clear()
    monkeypatch.setattr(gp, '_submit', worker)
    monkeypatch.setattr(gp.time, 'sleep', advance)
    run = gp.run_graph_pipeline(graph, 'PD-L1 four independent fixture chunks', 'owner', 'project', reuse_completed=False)
    assert read_json(run/'result.json')['success'] is True
    assert maximum[0] == 4


@pytest.mark.parametrize('mode', wm.MODES)
def test_failed_generator_keeps_other_branch_results(pdl1_harness, monkeypatch, mode):
    h = pdl1_harness
    graph = cp.new_pipeline(['rfd4','rfd4'] if mode=='modular_design' else ['bindcraft2','bindcraft2'], mode=mode)
    graph['nodes'][0]['data']['params'].update(target_pdb=h.target,target_chains='A')
    source=h.source('success', 'rfd4_proteina' if mode=='modular_design' else 'bindcraft2',2)
    graph['nodes'][1]['data']['retained_run_dir']=str(source)
    monkeypatch.setattr(cp,'preflight',lambda *a:None)
    submit=gp._submit
    def fail_branch(*args,**kwargs):
        if pg.NODE_CATALOG[args[2]['data']['kind']]['category'] in {'Design engine','Design pipeline'}:
            raise ValueError('deterministic unavailable generator branch')
        return submit(*args,**kwargs)
    monkeypatch.setattr(gp,'_submit',fail_branch)
    run=gp.run_graph_pipeline(graph,'PD-L1 failed generator fixture','owner','project',reuse_completed=False)
    result=read_json(run/'result.json')
    assert result['success'] is True
    assert len(result['outputs']['candidates'])==2
    assert result['outputs']['comparison_status']=='partial'


@pytest.mark.parametrize('kind', ['bindcraft2','boltzgen','rfd4','rfd3'])
def test_native_distributed_manifest_preserves_one_whole_invocation(pdl1_harness, kind):
    from mn_molecular_design_studio import candidate_execution as ce
    h=pdl1_harness
    graph=cp.new_pipeline([kind],mode='native_complete' if kind in wm.COMPLETE else 'modular_design')
    node=graph['nodes'][1]
    node['native_invocation']=True
    node['data']['params']['num_designs']=7
    if 'seed' in node['data']['params']:
        node['data']['params']['seed']=71
    original=copy.deepcopy(node['data']['params'])
    captured=[]
    def request(kind,params,source,name):
        captured.append(copy.deepcopy(params))
        return gp.node_request(kind,params,source,name)
    batches=ce.initialize(h.tmp_path,node,{'target':dict(target_pdb=h.target,target_chains=['A'],hotspots='',coldspots='')},
        original,'a'*64,ce.policy(node['data']['execution']),make_source=gp._local_source,make_request=request,read_rows=gp._absolute_candidates)
    assert len(batches)==1 and batches[0]['unit_count']==1
    assert ce.read_manifest(batches[0])['unit_kind']=='native_invocation'
    assert captured==[original]  # Native design counts, seed and GPU microbatch controls untouched.


def test_optional_interface_fixed_redesign_setup():
    from streamlit.testing.v1 import AppTest
    app=AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.scientific_tasks_ui import render_task_setup
graph=render_task_setup('redesign-test')
if graph is not None:
    st.session_state['created']=graph
""").run(timeout=45)
    app.selectbox(key='redesign-test-workflow-mode').select('modular_design').run()
    app.checkbox(key='redesign-test-mode-interface-redesign').check().run()
    app.number_input(key='redesign-test-mode-redesign-cutoff').set_value(5.0).run()
    app.button(key='redesign-test-create-mode').click().run()
    assert not app.exception
    graph=pg.validate_graph(app.session_state['created'])['graph']
    params=[n['data']['params'] for n in graph['nodes'] if n['data']['params'].get('sequence_stage')=='soluble_optimization']
    assert params and all(p['optimization_enabled'] and p['fix_interface_residues'] and p['interface_cutoff']==5.0 for p in params)
