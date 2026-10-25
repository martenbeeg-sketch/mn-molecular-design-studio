"""Production batch allocation, cumulative retention, and bounded preparation."""
from pathlib import Path
from types import SimpleNamespace
import json
import subprocess

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows import active_learning_controller as controller


def test_multi_gpu_expands_one_host_budget_without_duplicate_cpu_reservations(monkeypatch):
    from mn_molecular_design_studio.cluster import config, distributed_library
    monkeypatch.setattr(config, 'identity', lambda: 'local')
    monkeypatch.setattr(distributed_library, '_node_specs', lambda *_args, **_kwargs: [{
        'node_id': 'local', 'name': 'Organizer', 'local': True, 'cpu_cores': 32,
        'available_cpu': 62, 'free_gpus': [{'gpu_id': '0'}, {'gpu_id': '1'}],
    }])
    workers = controller._selected_worker_specs({'worker_nodes': [{
        'node_id': 'local', 'cpu_cores': 32, 'gpu_devices': ['0', '1'], 'batch_share': 1,
    }]})
    assert [w['cpu_cores'] for w in workers] == [16, 16]
    assert [w['gpu_device'] for w in workers] == ['0', '1']
    assert len({w['lane_id'] for w in workers}) == 2


def test_weighted_allocation_preserves_whole_10000_launches():
    rows = list(range(200_000))
    workers = [{'batch_share': share} for share in [1, 1, 3, 3]]
    assigned = controller._weighted_batch_assignments(rows, workers, range(len(rows)), 10_000)
    assert [len(group) for group in assigned] == [30_000, 30_000, 70_000, 70_000]
    assert sorted(value for group in assigned for value in group) == rows
    assert all(len(group) % 10_000 == 0 for group in assigned)


def test_retained_winners_survive_worker_cleanup_and_are_pruned_when_overtaken(tmp_path):
    source = tmp_path / 'worker/artifacts/poses'
    source.mkdir(parents=True)
    first, second = source / 'one.pdbqt', source / 'two.pdbqt'
    first.write_text('POSE ONE')
    second.write_text('POSE TWO')
    def row(key, name, score):
        return dict(library_id='lib', compound_key=key, compound_id=name,
            best_score=score, pose_path=f'artifacts/poses/{name}.pdbqt',
            worker_run_dir=str(tmp_path / 'worker'), worker_job_code='JOB01')
    one, two = row(1, 'one', -9), row(2, 'two', -8)
    campaign = tmp_path / 'campaign'
    kept = controller._collect_best_poses(campaign, [one], 1, ssh_user='compute')
    first.unlink()
    assert controller._collect_best_poses(campaign, [one, two], 1, ssh_user='compute') == kept
    two['best_score'] = -10
    updated = controller._collect_best_poses(campaign, [two, one], 1, ssh_user='compute')
    assert (campaign / updated[0]).read_text() == 'POSE TWO'
    assert not (campaign / kept[0]).exists()


def test_batch_runner_accepts_partial_output_without_any_rescue(tmp_path, monkeypatch):
    from rdkit import Chem
    from mn_molecular_design_studio.docking_batch_preparation import PREPARATION_RUNNER
    from mn_molecular_design_studio.docking_preparation_protocols import scrub_invocation
    for name in ('input', 'prepared', 'logs', 'ligands'):
        (tmp_path / 'artifacts' / name).mkdir(parents=True)
    inputs = tmp_path / 'artifacts/input'
    (inputs / 'library.smi').write_text('CCO good\nCCN bad\n')
    (inputs / 'preparation_config.json').write_text(json.dumps({
        'route': scrub_invocation('molscrub_0_1_0', ph=7.4, skip_tautomers=True, cpu=2), 'workers': 2}))
    calls = []
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(subprocess, 'check_output', lambda *_a, **_k: '0.1.0\n')
    def run(command, **kwargs):
        calls.append(command)
        if command[0].endswith('/scrub.py'):
            molecule = Chem.MolFromSmiles('CCO')
            molecule.SetProp('_Name', 'good')
            writer = Chem.SDWriter(command[command.index('-o') + 1])
            writer.write(molecule)
            writer.close()
            return SimpleNamespace(returncode=1)
        assert command[0] == 'mk_prepare_ligand.py'
        output = Path(command[command.index('--multimol_outdir') + 1])
        (output / 'good.pdbqt').write_text('ROOT\nATOM\nTORSDOF 1\n')
        return SimpleNamespace(returncode=1)
    monkeypatch.setattr(subprocess, 'run', run)
    exec(compile(PREPARATION_RUNNER, '<preparation>', 'exec'), {})
    summary = json.loads((tmp_path / 'artifacts/preparation_summary.json').read_text())
    assert summary['retry_count'] == 0
    assert summary['failed_input_count'] == 1
    assert summary['output_pdbqt_count'] == 1
    assert len(calls) == 2
    assert calls[0][calls[0].index('--cpu') + 1] == '2'


def test_docking_round_runs_both_local_gpu_workers_and_steals_pending_batches(tmp_path, monkeypatch):
    import threading
    import time
    rows = [dict(library_id='lib', compound_key=i, compound_id=f'c{i}', smiles='CCO',
                 fingerprint=bytes(128)) for i in range(16)]
    workers = [dict(node_id='local', node_name='Organizer', local=True, cpu_cores=16,
                    gpu_device=str(i), lane_id=f'local:{i}', batch_share=1) for i in range(2)]
    fast_finished = threading.Event()
    calls, simultaneous, maximum = [], set(), [0]
    lock = threading.Lock()
    monkeypatch.setattr(controller, '_docking_params', lambda *args, worker, **kw: {
        'cpu_cores': worker['cpu_cores'], 'gpu_device': worker['gpu_device']})
    def dispatch(*, worker, request, **kwargs):
        with lock:
            calls.append(worker['gpu_device'])
            simultaneous.add(worker['gpu_device'])
            maximum[0] = max(maximum[0], len(simultaneous))
        if worker['gpu_device'] == '0':
            assert fast_finished.wait(5)
        else:
            time.sleep(0.01)
        with lock:
            simultaneous.discard(worker['gpu_device'])
        return dict(status='completed', worker_run_dir='/run', job_code='TEST', gpu=worker['gpu_device'])
    def scores(job, mapping):
        return list(mapping.values()), []
    def progress(items):
        if any(item['gpu_device'] == '1' and item['completed_batches'] >= 7 for item in items):
            fast_finished.set()
    monkeypatch.setattr(controller, '_dispatch_docking_shard', dispatch)
    monkeypatch.setattr(controller, '_read_worker_scores', scores)
    scored, failed = controller._dock_round(campaign=tmp_path, round_number=1, rows=rows, workers=workers,
        docking={'worker_batch_size': 2}, target={}, owner_id='owner', project_id='project',
        ssh_user='compute', work_directory=tmp_path / 'round', cancel_check=lambda: False, worker_progress=progress)
    assert len(scored) == 16 and not failed
    assert len({row['compound_key'] for row in scored}) == 16
    assert maximum[0] == 2
    assert calls.count('1') == 7 and calls.count('0') == 1
