"""Small real scheduler/container benchmarks; no downloads or remote changes."""
import copy
import json
import os
import subprocess
import time
from pathlib import Path

import pytest

from mn_molecular_design_studio import candidate_execution as ce, mpnn_input, workspace
from mn_molecular_design_studio.engine.workflow_engine.core import jobs, scheduler, workflow_queue
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as gp, optimized_prediction
from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG


@pytest.mark.gpu_e2e
@pytest.mark.parametrize('kind', ['ligandmpnn', 'boltz2'])
def test_real_candidate_scheduler_benchmark(gpu_runtime, monkeypatch, kind):
    runtime = gpu_runtime
    root = runtime['root']
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(root / 'cluster'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_APP_HOME', str(root / 'app'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(root / 'app'))
    monkeypatch.setattr(workspace, 'validate_project_id', lambda value: value)
    # This test drives the production local scheduler explicitly so it neither
    # starts nor reconfigures a user's long-lived scheduler service.
    monkeypatch.setattr(workflow_queue, 'spawn_worker_for_run', lambda path: None)
    source_path = os.environ.get('MDS_CANDIDATE_BENCHMARK_PREDICTION_SOURCE') if kind == 'boltz2' else os.environ.get('MDS_CANDIDATE_BENCHMARK_SOURCE')
    if not source_path:
        pytest.skip('Set MDS_CANDIDATE_BENCHMARK_SOURCE to a retained small generation run')
    retained = Path(source_path).resolve()
    assert jobs.read_json(retained / 'result.json')['success'] is True
    source_rows = gp._absolute_candidates({'run_dir': str(retained), 'candidates_jsonl': str(retained / 'artifacts/normalized_candidates/candidates.jsonl')})
    count = int(os.environ.get('MDS_CANDIDATE_BENCHMARK_COUNT', '4'))
    rows = []
    for i in range(count):
        row = copy.deepcopy(source_rows[i % len(source_rows)])
        row['candidate_id'] = f'benchmark-input-{i:03d}'
        rows.append(row)
    source = gp._local_source(root, 'benchmark-input', rows)
    runtime['require_image'](mpnn_input.IMAGE if kind == 'ligandmpnn' else optimized_prediction.ENGINE_PROFILES['boltz2']['image'])
    params = {f['name']: f['default'] for f in NODE_CATALOG[kind]['fields']}
    if kind == 'ligandmpnn':
        assert mpnn_input.checkpoint('protein_mpnn').is_file()
        params.update(sequences=1, seed=20261027, backend='shared_ligandmpnn', sampling_temp=.4)
    else:
        for name in ('boltz2_conf.ckpt', 'boltz2_aff.ckpt', 'ccd.pkl'):
            assert (runtime['reference'] / 'boltz_models' / name).is_file()
        params.update(samples=1, num_recycles=1, sampling_steps=20, seed=20261027,
                      target_use_msa=False, use_target_template=False, generate_ipsae=False)
    node = dict(id=kind, data=dict(kind=kind, label=kind))
    root.joinpath('benchmark.json').write_text('{}')
    evidence = dict(kind=kind, retained_source=str(retained), input_ids=[r['candidate_id'] for r in rows],
                    effective_parameters=params, gpu_inventory=json.loads((root / 'runtime.json').read_text()), trials=[])
    for label, batch_size, device in [('single_gpu_native_batch', count, '0'), ('single_gpu_shards', 2, '0'),
                                       ('two_gpus_shards', 2, 'all')]:
        parent = root / label; parent.mkdir()
        settings = ce.policy(dict(mode='sharded', batch_size=batch_size, max_in_flight=2, retries=0, min_vram_gb=8))
        current = {}
        service = scheduler.WorkerScheduler(root=root / 'runs', cpu_slots=8, gpu_devices=['0', '1'])
        start = time.monotonic(); telemetry = []
        controller_seconds = scheduler_seconds = 0.
        def persist():
            jobs.write_json(parent / 'graph_execution.json', {'nodes': {kind: current}})
        try:
            result = None
            while result is None and time.monotonic() - start < 900:
                tick_started = time.monotonic()
                result = ce.tick(parent, node, current, source, params, ce.digest([kind, label, params]), settings,
                    graph=gp, owner_id='benchmark-user', project_id='benchmark-project', gpu_device=device, persist=persist)
                controller_seconds += time.monotonic() - tick_started
                scheduler_started = time.monotonic()
                service.schedule_once()
                scheduler_seconds += time.monotonic() - scheduler_started
                telemetry.append({'elapsed': time.monotonic()-start, 'gpus': subprocess.run(
                    ['nvidia-smi', '--query-gpu=index,utilization.gpu,memory.used', '--format=csv,noheader,nounits'],
                    capture_output=True, text=True, check=True).stdout})
                if result is None:
                    time.sleep(.5)
            assert result is not None, current
            assert result[1] == 'completed', current
            output = gp._absolute_candidates(result[0])
            assert len(output) == count
            assert {p for row in output for p in row.get('parents', [])} == set(evidence['input_ids'])
            assert all(row['binder_sequence'] for row in output)
            elapsed = time.monotonic()-start
            allocations = [b['attempts'][-1]['allocation'] for b in current['batches']]
            if label == 'two_gpus_shards':
                assert {str(a['gpu_device']) for a in allocations} == {'0', '1'}, allocations
            trial = dict(label=label, wall_seconds=elapsed, candidates_per_hour=count*3600/elapsed,
                         candidate_count=len(output), allocations=allocations, state=current, telemetry=telemetry,
                         controller_seconds=controller_seconds, scheduler_seconds=scheduler_seconds,
                         transfer_bytes=0, transfer_seconds=0, retries=current['candidate_progress']['retry_count'],
                         model_initialization_seconds=None, native_cpu_preparation_seconds=None,
                         instrumentation_note='Native model initialization and native CPU preparation are included in wall time, without separate engine timers.')
            evidence['trials'].append(trial)
            root.joinpath('benchmark.json').write_text(json.dumps(evidence, indent=2))
        finally:
            service._refresh_active()
            for worker in list(service.active.values()):
                if worker.process.poll() is None:
                    if jobs.read_json(worker.run_dir / 'metadata.json').get('status') not in {'completed', 'failed', 'cancelled'}:
                        jobs.stop_job(worker.run_dir, reason='Candidate benchmark teardown')
                    worker.process.wait(timeout=60)
            service._refresh_active()

    if kind == 'boltz2':
        _run_cpu_analysis(runtime, result[0], count, evidence)


def _run_cpu_analysis(runtime, cpu_source, count, evidence):
    from mn_molecular_design_studio.quality_tools import WORKFLOWS
    runtime['require_image'](WORKFLOWS['analysis.sequence_qc']['image'])
    cpu_node = dict(id='sequence_qc', data=dict(kind='sequence_qc', label='Sequence QC',
        params={f['name']: f['default'] for f in NODE_CATALOG['sequence_qc']['fields']}))
    root = runtime['root']
    parent = root / 'cpu-analysis'; parent.mkdir()
    current = {}
    service = scheduler.WorkerScheduler(root=root / 'runs', cpu_slots=8, gpu_devices=['0', '1'])
    start = time.monotonic()
    result = None
    try:
        while result is None and time.monotonic() - start < 180:
            result = ce.tick(parent, cpu_node, current, cpu_source, cpu_node['data']['params'], ce.digest(['cpu-qc']),
                ce.policy(dict(mode='sharded', batch_size=2, retries=0)), graph=gp,
                owner_id='benchmark-user', project_id='benchmark-project', gpu_device='all', persist=lambda: None)
            service.schedule_once()
            if result is None:
                time.sleep(.5)
        assert result is not None and result[1] == 'completed', current
        assert result[0]['candidate_count'] == count
        assert all(b['attempts'][-1]['allocation']['gpu_device'] in {'none', 'cpu'} for b in current['batches'])
        assert all(not b['attempts'][-1]['allocation'].get('gpu_devices') for b in current['batches'])
        evidence['cpu_analysis'] = dict(wall_seconds=time.monotonic()-start, state=current, output=result[0])
        root.joinpath('benchmark.json').write_text(json.dumps(evidence, indent=2))
    finally:
        service._refresh_active()
        for worker in list(service.active.values()):
            if worker.process.poll() is None:
                if jobs.read_json(worker.run_dir / 'metadata.json').get('status') not in {'completed', 'failed', 'cancelled'}:
                    jobs.stop_job(worker.run_dir, reason='CPU benchmark teardown')
                worker.process.wait(timeout=60)
        service._refresh_active()


@pytest.mark.gpu_e2e
def test_real_cpu_analysis_of_retained_gpu_candidates(gpu_runtime, monkeypatch):
    runtime = gpu_runtime
    root = runtime['root']
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(root / 'cluster'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_APP_HOME', str(root / 'app'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(root / 'app'))
    monkeypatch.setattr(workspace, 'validate_project_id', lambda value: value)
    monkeypatch.setattr(workflow_queue, 'spawn_worker_for_run', lambda path: None)
    retained = os.environ.get('MDS_CANDIDATE_BENCHMARK_CPU_SOURCE')
    if not retained:
        pytest.skip('Set MDS_CANDIDATE_BENCHMARK_CPU_SOURCE to a retained real prediction collection')
    source = Path(retained).resolve()
    assert jobs.read_json(source / 'result.json')['success'] is True
    collection = dict(run_dir=str(source), candidates_jsonl=str(source / 'artifacts/normalized_candidates/candidates.jsonl'))
    count = len(gp._absolute_candidates(collection))
    assert count > 0
    _run_cpu_analysis(runtime, collection, count, dict(retained_gpu_source=str(source)))
