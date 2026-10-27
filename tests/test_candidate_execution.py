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
