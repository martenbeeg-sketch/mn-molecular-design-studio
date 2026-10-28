import copy
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import candidate_execution as ce
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as gp
from mn_molecular_design_studio.pipeline_graph import clean_graph, new_graph


def batches(tmp_path, rows=None):
    rows = rows or [dict(candidate_id=f'c{i}', binder_sequence='AAA', parents=['backbone']) for i in range(5)]
    source = gp._local_source(tmp_path, 'source', rows)
    node = dict(id='predict', data=dict(kind='boltz2', label='Prediction'))
    def request(kind, params, source, name):
        return dict(workflow='mock', parameters=dict(source_run_dir=source['run_dir'], candidates_jsonl=source['candidates_jsonl'], **params))
    return ce.initialize(tmp_path, node, source, {'seed': 42}, 'a' * 64, ce.policy({'mode': 'sharded', 'batch_size': 2}),
                         make_source=gp._local_source, make_request=request, read_rows=gp._absolute_candidates)


def test_shards_preserve_order_lineage_and_immutable_inputs(tmp_path):
    shards = batches(tmp_path)
    assert [b['unit_count'] for b in shards] == [2, 2, 1]
    assert [ce.read_manifest(b)['input_candidate_ids'] for b in shards] == [['c0', 'c1'], ['c2', 'c3'], ['c4']]
    assert all(ce.read_manifest(b)['request']['parameters']['seed'] == 42 for b in shards)
    assert batches(tmp_path) == shards
    path = Path(ce.read_manifest(shards[0])['request']['parameters']['candidates_jsonl'])
    path.write_text('changed')
    with pytest.raises(ValueError, match='changed'):
        ce.read_manifest(shards[0])


def test_collect_and_aggregate_ignore_completion_order_and_retain_provenance(tmp_path):
    shards = batches(tmp_path)
    for batch in reversed(shards):
        manifest = ce.read_manifest(batch)
        structure = tmp_path / (batch['batch_id'] + '.pdb'); structure.write_text('structure')
        rows = [dict(candidate_id='native-' + c, complex_pdb=str(structure), parents=[c]) for c in manifest['input_candidate_ids']]
        batch['attempts'] = [dict(worker='worker-1', allocation={'gpu_device': '1'}, software={'image': 'digest'})]
        output = gp._local_source(tmp_path, 'out-' + batch['batch_id'], rows)
        result = ce.collect(batch, output, read_rows=gp._absolute_candidates)
        path = tmp_path / (batch['batch_id'] + '.json'); path.write_text(json.dumps(result))
        batch['rows_path'] = str(path)
        batch['output_hashes'][str(path)] = ce.file_digest(path)
    result = ce.aggregate(list(reversed(shards)))
    assert [r['parents'][0] for r in result] == ['c0', 'c1', 'c2', 'c3', 'c4']
    assert len({r['candidate_id'] for r in result}) == 5
    assert result[0]['raw_metadata']['candidate_execution']['worker'] == 'worker-1'
    assert ce.progress(shards)['completed_units'] == 5
    Path(result[0]['complex_pdb']).write_text('corrupt')
    with pytest.raises(ValueError, match='checksum'):
        ce.aggregate(shards)


def test_duplicate_and_unrelated_results_are_rejected(tmp_path):
    shard = batches(tmp_path)[0]
    shard['attempts'] = [{}]
    output = gp._local_source(tmp_path, 'out', [dict(candidate_id='bad', parents=['other'])])
    with pytest.raises(ValueError, match='lineage'):
        ce.collect(shard, output, read_rows=gp._absolute_candidates)
    output = gp._local_source(tmp_path, 'out', [dict(candidate_id='same')] * 2)
    with pytest.raises(ValueError, match='duplicated'):
        ce.collect(shard, output, read_rows=gp._absolute_candidates)


def test_policy_is_opt_in_and_validated():
    assert ce.policy()['mode'] == 'native'
    for change in ({'batch_size': 0}, {'retries': True}, {'min_vram_gb': float('nan')}, {'placement': 'unknown'}, {'bogus': 1}):
        with pytest.raises((ValueError, TypeError)):
            ce.policy(change)


def test_tick_runs_bounded_batches_retries_once_and_replays_without_duplicates(tmp_path, monkeypatch):
    from mn_molecular_design_studio import runtime
    monkeypatch.setattr(runtime, 'runs_root', lambda: tmp_path / 'runs')
    source = gp._local_source(tmp_path, 'source', [dict(candidate_id=f'c{i}') for i in range(3)])
    node = dict(id='predict', data=dict(kind='boltz2', label='Prediction'))
    monkeypatch.setattr(gp, 'node_request', lambda kind, params, source, name:
        dict(workflow='mock', parameters=dict(source_run_dir=source['run_dir'], candidates_jsonl=source['candidates_jsonl'], seed=42)))
    submissions = []
    def submit(request, parent, node, **kwargs):
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json
        path = tmp_path / 'runs' / 'prediction' / node['id']; path.mkdir(parents=True)
        failed = not submissions
        rows = gp._absolute_candidates({'run_dir': request['parameters']['source_run_dir'], 'candidates_jsonl': request['parameters']['candidates_jsonl']})
        submissions.append(copy.deepcopy(request))
        normalized = path / 'artifacts/normalized_candidates'; normalized.mkdir(parents=True)
        (normalized / 'candidates.jsonl').write_text(''.join(json.dumps(dict(candidate_id='result-' + r['candidate_id'], parents=[r['candidate_id']])) + '\n' for r in rows))
        write_json(path / 'metadata.json', {'status': 'failed' if failed else 'completed'})
        write_json(path / 'result.json', {'success': not failed})
        write_json(path / 'input.json', {'params': request['parameters']})
        return path
    monkeypatch.setattr(gp, '_submit', submit)
    current = {}
    settings = ce.policy(dict(mode='sharded', batch_size=1, max_in_flight=1, retries=1))
    arguments = dict(graph=gp, owner_id='user', project_id='project', gpu_device='all', persist=lambda: None)
    result = None
    for _ in range(10):
        result = ce.tick(tmp_path, node, current, source, {'seed': 42}, 'a' * 64, settings, **arguments)
        assert sum(b['status'] not in ce.TERMINAL | {'pending', 'retrying'} for b in current['batches']) <= 1
        if result:
            break
    assert result[1] == 'completed'
    assert result[0]['candidate_count'] == 3
    assert len(submissions) == 4
    assert submissions[0] == submissions[1]
    assert current['candidate_progress']['retry_count'] == 1
    original = gp._absolute_candidates(result[0])
    replay, status = ce.tick(tmp_path, node, current, source, {'seed': 42}, 'a' * 64, settings, **arguments)
    assert len(submissions) == 4
    assert gp._absolute_candidates(replay) == original


def test_scheduler_preserves_auto_and_admits_single_gpu(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core import scheduler, central_scheduler
    path = tmp_path / 'job'; path.mkdir()
    (path / 'metadata.json').write_text(json.dumps({'tool': 'boltz2', 'queue_resource': 'gpu'}))
    (path / 'worker_request.json').write_text(json.dumps({'resources': {'gpu_device': 'auto', 'min_vram_gb': 12}}))
    request = scheduler.request_for_run(path)
    assert request.gpu_device == 'auto'
    assert scheduler._gpu_candidates(request, ['0', '1']) == (['0', '1'], False)
    assert central_scheduler._central_resource_request(path)['gpu_count'] == 1
    assert central_scheduler._central_resource_request(path)['min_vram_gb'] == 12


def test_missing_candidate_outputs_are_partial_and_progress_is_real(tmp_path):
    shard = batches(tmp_path)[0]
    shard['attempts'] = [{}]
    output = gp._local_source(tmp_path, 'out', [dict(candidate_id='prediction', parents=['c0'])])
    result = ce.collect(shard, output, read_rows=gp._absolute_candidates)
    assert shard['status'] == 'partial'
    assert shard['completed_units'] == 1
    assert shard['missing_candidate_ids'] == ['c1']
    assert ce.progress([shard])['completed_units'] == 1


def test_cancellation_stops_queued_and_running_batches(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    stopped = []
    monkeypatch.setattr(jobs, 'stop_job', lambda path, **kwargs: stopped.append(path))
    state = {'nodes': {'p': {'batches': [
        dict(status='queued', job_run_dir=str(tmp_path / 'q')),
        dict(status='running', job_run_dir=str(tmp_path / 'r')),
        dict(status='completed', job_run_dir=str(tmp_path / 'c'))]}}}
    ce.cancel(tmp_path, state)
    assert stopped == [tmp_path / 'q', tmp_path / 'r']
    assert [b['status'] for b in state['nodes']['p']['batches']] == ['cancelled', 'cancelled', 'completed']


def test_disconnected_worker_is_quarantined_without_duplicate_execution(tmp_path, monkeypatch):
    import time
    from mn_molecular_design_studio.cluster import candidate_batches as peers
    batch = batches(tmp_path)[0]
    batch.update(cluster_job_id='remote', status='running', attempts=[{}], disconnected_since=time.time() - 120)
    current = {'batches': [batch], 'fingerprint': 'a' * 64}
    cancelled = []
    monkeypatch.setattr(peers, 'poll', lambda *a, **k: ('disconnected', None, {'node_id': 'offline-worker'}))
    monkeypatch.setattr(peers, 'cancel', lambda *a, **k: cancelled.append('remote'))
    monkeypatch.setattr(peers, 'submit', lambda *a, **k: pytest.fail('Must not duplicate a job whose worker may still be alive'))
    node = dict(id='predict', data=dict(kind='boltz2', label='Prediction'))
    result = ce.tick(tmp_path, node, current, {}, {}, 'a' * 64,
        ce.policy(dict(mode='sharded', placement='auto', stall_timeout_seconds=60)), graph=gp,
        owner_id='user', project_id='project', gpu_device='all', persist=lambda: None)
    assert result[1] == 'failed'
    assert batch['status'] == 'lost' and batch['quarantined'] is True
    assert len(batch['attempts']) == 1 and cancelled == ['remote']


def test_incremental_cpu_analysis_advances_completed_batches_before_prediction_finishes(tmp_path, monkeypatch):
    from mn_molecular_design_studio import cache_provenance
    monkeypatch.setattr(cache_provenance, 'execution_identity', lambda *a, **k: {'version': 'tested'})
    assert ce.incremental_cpu_kind('sequence_qc')
    assert ce.incremental_cpu_kind('ipsae')
    assert not ce.incremental_cpu_kind('boltz2')
    calls = []
    def tick(parent, node, state, source, *args, **kwargs):
        calls.append(node['id'])
        assert kwargs['gpu_device'] == 'none'
        return source, 'completed'
    monkeypatch.setattr(ce, 'tick', tick)
    upstream = {'status': 'running', 'batches': []}
    for i in range(2):
        path = tmp_path / f'rows{i}.json'
        path.write_text(json.dumps([dict(candidate_id=f'prediction{i}', parents=[f'input{i}'])]))
        upstream['batches'].append(dict(batch_id=f'b{i}', status='completed' if i == 0 else 'running',
                                       rows_path=str(path), manifest_sha256=str(i)))
    source = gp._local_source(tmp_path, 'probe', [dict(candidate_id='prediction0')])
    node = dict(id='qc', data=dict(kind='sequence_qc', label='QC', params={}))
    current = {}
    args = dict(graph=gp, owner_id='user', project_id='project', persist=lambda: None)
    settings = ce.policy(dict(mode='sharded', incremental_cpu=True))
    assert ce.tick_incremental(tmp_path, node, current, upstream, source, settings, **args) is None
    assert calls == ['qc_b0']  # CPU work advanced while the second GPU batch is still running.
    upstream['batches'][1]['status'] = 'completed'; upstream['status'] = 'completed'
    output, key, status = ce.tick_incremental(tmp_path, node, current, upstream, source, settings, **args)
    assert status == 'completed'
    assert calls == ['qc_b0', 'qc_b1']
    assert [r['candidate_id'] for r in gp._absolute_candidates(output)] == ['prediction0', 'prediction1']


def test_saved_graphs_retain_native_semantics_and_explicit_execution_policy():
    graph = new_graph('blank')
    assert 'candidate_execution' not in clean_graph(graph)
    settings = dict(mode='sharded', failure_policy='partial', placement='local', batch_size=2, incremental_cpu=True)
    graph['candidate_execution'] = settings
    assert clean_graph(graph)['candidate_execution'] == ce.policy(settings)


def test_controller_retry_waits_for_scheduler_retry_and_live_process_termination(monkeypatch):
    from types import SimpleNamespace
    from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler, scheduler
    state = {'status': 'queued'}
    monkeypatch.setattr(central_scheduler, '_client', lambda: SimpleNamespace(get_job=lambda job: state))
    monkeypatch.setattr(scheduler, '_pid_alive', lambda pid: pid == 123)
    metadata = dict(status='failed', scheduler_job_id='original', worker_pid=999)
    assert ce.execution_status(metadata) == 'scheduler_retrying'
    state['status'] = 'failed'
    assert ce.execution_status(metadata) == 'failed'
    metadata['worker_pid'] = 123
    assert ce.execution_status(metadata) == 'terminating'
    def unavailable():
        raise RuntimeError('scheduler disconnected')
    monkeypatch.setattr(central_scheduler, '_client', unavailable)
    assert ce.execution_status(metadata) == 'reconciling'


def test_partial_pool_context_survives_successful_intermediate_analysis():
    prediction = {'candidate_execution': {'missing_batches': ['b000002'], 'progress': {'completed_units': 2}}}
    analysis = {'partial_pool': ce.inherit_partial_pool({}, [prediction])}
    selection = {'partial_pool': ce.inherit_partial_pool({}, [analysis])}
    result = ce.inherit_partial_pool({}, [selection, analysis], missing_branches=['generator-b'])
    assert result['missing_branches'] == ['generator-b']
    assert len(result['candidate_batches']) == 1
    assert result['candidate_batches'][0]['missing_batches'] == ['b000002']
    assert ce.inherit_partial_pool({}, [{}]) is None


@pytest.mark.parametrize('development', [False, True])
def test_peer_batch_lifecycle_uses_durable_controller_identity(tmp_path, monkeypatch, development):
    import time
    from mn_molecular_design_studio import runtime
    from mn_molecular_design_studio.cluster import candidate_batches as peers, config
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    jobs.write_json(tmp_path / 'metadata.json', dict(studio_owner_id='owner',
        studio_project_id='project', studio_development_mode=development))
    monkeypatch.setattr(runtime, 'runs_root', lambda: tmp_path / 'runs')
    monkeypatch.setattr(config, 'enabled', lambda: True)
    batch = batches(tmp_path)[0]
    current = {'batches': [batch]}
    calls = []
    def submit(*args, **kwargs):
        calls.append(('submit', kwargs['allow_development']))
        return {'cluster_job_id': 'remote'}
    def poll(*args, **kwargs):
        calls.append(('poll', kwargs['allow_development']))
        return 'disconnected', None, {'node_id': 'offline'}
    def cancel(*args, **kwargs):
        calls.append(('cancel', kwargs['allow_development']))
        assert kwargs['owner_id'] == 'owner' and kwargs['project_id'] == 'project'
    monkeypatch.setattr(peers, 'submit', submit)
    monkeypatch.setattr(peers, 'poll', poll)
    monkeypatch.setattr(peers, 'cancel', cancel)
    node = dict(id='predict', data=dict(kind='boltz2', label='Prediction'))
    arguments = dict(graph=gp, owner_id='owner', project_id='project', gpu_device='all', persist=lambda: None)
    settings = ce.policy(dict(mode='sharded', placement='auto', stall_timeout_seconds=60))
    assert ce.tick(tmp_path, node, current, {}, {}, 'a' * 64, settings, **arguments) is None
    batch['disconnected_since'] = time.time() - 120
    assert ce.tick(tmp_path, node, current, {}, {}, 'a' * 64, settings, **arguments)[1] == 'failed'
    ce.cancel(tmp_path, {'nodes': {'predict': {'batches': [dict(status='queued', cluster_job_id='other')]}}})
    assert calls == [('submit', development), ('poll', development), ('cancel', development), ('cancel', development)]
