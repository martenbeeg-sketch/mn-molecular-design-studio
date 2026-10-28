"""Real Studio adapters with deterministic launches; no Docker/GPU inference.

Image-side tests use inspected third-party adapter sources, never copied into this
repository. Set MDS_BATCH_API_SOURCE to the pinned image source snapshot.
"""
import copy
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio import candidate_pipeline, native_execution
from mn_molecular_design_studio.engine.workflow_engine.core import jobs
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design as sequence, mpnn
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as prediction
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as graph

SIZES = [1, 5, 10, 25, 50, 100]


def evidence(name, **counts):
    directory = os.environ.get('MDS_NATIVE_EXECUTION_EVIDENCE_DIR')
    if directory:
        path = Path(directory); path.mkdir(parents=True, exist_ok=True)
        (path / (name + '.json')).write_text(json.dumps(dict(
            scientific_execution='deterministic process fixture; no inference',
            inputs='retained PD-L1 coordinates' if os.environ.get('MDS_RETAINED_PDL1') == '1' else 'deterministic collection',
            **counts), indent=2) + '\n')


@pytest.fixture
def collection(tmp_path, design_source_factory):
    source, original = design_source_factory()
    originals = [original]
    if os.environ.get('MDS_RETAINED_PDL1') == '1':
        from mn_molecular_design_studio.candidate_sets import iter_candidates
        evidence = json.loads((Path(__file__).resolve().parents[1] / 'docs/testing/evidence/phase7a-retained-graph.json').read_text())
        retained = Path(evidence['run_dir'])
        originals = list(iter_candidates(dict(run_dir=str(retained), candidates_jsonl=str(retained / 'artifacts/normalized_candidates/candidates.jsonl'))))
        assert len(originals) == 3
        for i, row in enumerate(originals):
            path = source / f'retained-{i}.pdb'
            native_path = Path(row['complex_pdb'])
            if not native_path.is_file():
                native_path = Path(str(native_path).replace('/home/user/.local/share/mn-molecular-design-studio/workdir/runs/', '/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs/'))
            mapping = mpnn._pdb_copy(native_path, path)
            row.update(complex_pdb=str(path), target_chains=[mapping.get(c, c) for c in row['target_chains']],
                       binder_chains=[mapping.get(c, c) for c in row['binder_chains']], raw_metadata={})
    structure = source / 'input.pdb'
    mpnn._pdb_copy(source / 'designed.cif', structure)
    structure.write_text('\n'.join(sequence._ligandmpnn_complete_backbone_lines(structure.read_text().splitlines())) + '\n')
    original.update(complex_pdb=str(structure), stage='generation.backbone')
    def build(n):
        rows = [copy.deepcopy(originals[i % len(originals)]) | {'candidate_id': f'input_{i:04d}', 'stage': 'generation.backbone'} for i in range(n)]
        if len(originals) == 1:
            for row in rows: row['complex_pdb'] = str(structure)
        path = source / 'collection.jsonl'
        path.write_text(''.join(json.dumps(r) + '\n' for r in rows))
        return source, path, rows
    return build


def host_path(run, value):
    return run / value.removeprefix('/work/')


@pytest.mark.parametrize('n', SIZES)
@pytest.mark.parametrize('backend', ['shared_ligandmpnn', 'foundry'])
def test_sequence_native_process_grouping(n, backend, collection, monkeypatch, tmp_path):
    source, path, rows = collection(n)
    weights = tmp_path / 'weights.pt'; weights.write_bytes(b'fixture')
    from mn_molecular_design_studio import mpnn_input, prediction_provenance
    monkeypatch.setattr(mpnn_input, 'checkpoint', lambda *a, **k: weights)
    monkeypatch.setattr(prediction_provenance, 'image_identity', lambda *a: {})
    launches, outputs = [], []
    for offset in range(0, n, 5):
        run = tmp_path / f'run-{offset}'; run.mkdir()
        jobs.write_json(run / 'metadata.json', {'resource_allocation': {'gpu_devices': ['0']}})
        chunk = source / f'chunk-{offset}.jsonl'
        chunk.write_text(''.join(json.dumps(r) + '\n' for r in rows[offset:offset+5]))
        def launch(command, **kwargs):
            assert command[:2] == ['docker', 'run']
            launches.append(command)
            def argument(key): return command[command.index(key) + 1]
            if backend == 'foundry':
                output = host_path(run, argument('--out_directory')); output.mkdir(parents=True, exist_ok=True)
                shutil.copy2(host_path(run, argument('--structure_path')), output / 'design.cif')
            else:
                paths = json.loads(host_path(run, argument('--pdb_path_multi')).read_text()) if '--pdb_path_multi' in command else [argument('--pdb_path')]
                if '--pdb_path_multi' in command:
                    fixed = json.loads(host_path(run, argument('--fixed_residues_multi')).read_text())
                    assert set(fixed) == set(paths)
                output = host_path(run, argument('--out_folder'))
                for value in paths:
                    pdb = host_path(run, value)
                    for folder in ('seqs', 'backbones'): (output / folder).mkdir(parents=True, exist_ok=True)
                    shutil.copy2(pdb, output / 'backbones' / f'{pdb.stem}_1.pdb')
                    sequences = sequence.refolding_workflow._sequences_by_chain(pdb)
                    (output / 'seqs' / f'{pdb.stem}.fa').write_text('>fixture, id=1\n' + ':'.join(sequences[c] for c in sorted(sequences)) + '\n')
            return SimpleNamespace(returncode=0)
        monkeypatch.setattr(sequence.subprocess, 'run', launch)
        mpnn.run_mpnn_sequence_design(source, chunk, backend=backend, sequence_stage='initial', num_seq_per_target=1,
                                      seed=42, existing_job=SimpleNamespace(run_dir=run))
        result = jobs.read_json(run / 'result.json')
        assert result['success'], result
        outputs.extend(read_candidates(run))
        receipt = jobs.read_json(run / 'native_execution.json')
        expected = min(5, n-offset) if backend == 'foundry' else 1
        assert receipt['docker_processes'] == receipt['native_invocations'] == expected
        assert receipt['gpu_allocation_events'] == 1
    assert len(launches) == (n if backend == 'foundry' else math.ceil(n/5))
    assert sorted(r['parents'][0] for r in outputs) == [r['candidate_id'] for r in rows]
    evidence(f'{backend}-{n}', candidates=n, chunks=math.ceil(n/5), docker_processes=len(launches),
             native_invocations=len(launches), gpu_allocation_events=math.ceil(n/5),
             gpu_allocation_basis='deterministic scheduler metadata fixture', outputs=len(outputs),
             exact_parent_ids=[r['parents'][0] for r in outputs])


@pytest.fixture
def image_adapter(monkeypatch, request):
    source = Path(os.environ.get('MDS_BATCH_API_SOURCE', Path(__file__).resolve().parents[2] / 'mn-tool-containers/shared/batch_api'))
    if (source / 'boltz/__init__.py').is_file():
        source = source / ('protenix' if request.node.callspec.params['engine'] == 'protenix_v2' else 'boltz')
    if not (source / '__init__.py').is_file():
        pytest.skip('image-side source unavailable; set MDS_BATCH_API_SOURCE')
    spec = importlib.util.spec_from_file_location('mn_batch_api', source / '__init__.py', submodule_search_locations=[str(source)])
    package = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, 'mn_batch_api', package); spec.loader.exec_module(package)
    for name in ('contract', 'adapters', 'legacy_adapters'):
        monkeypatch.delitem(sys.modules, 'mn_batch_api.' + name, raising=False)
    import mn_batch_api.adapters as adapter
    import mn_batch_api.contract as contract
    return adapter, contract


@pytest.mark.parametrize('n', SIZES)
@pytest.mark.parametrize('engine,seeds', [('boltz2', 1), ('boltz2', 3), ('protenix_v2', 1)])
def test_prediction_reaches_image_native_directory_launcher(n, engine, seeds, collection, image_adapter, monkeypatch, tmp_path):
    adapter, contract = image_adapter
    source, path, rows = collection(n)
    from mn_molecular_design_studio import prediction_provenance
    monkeypatch.setattr(prediction_provenance, 'image_identity', lambda *a: {})
    monkeypatch.setattr(prediction, 'reference_root', lambda: tmp_path / 'reference')
    cache = tmp_path / 'reference/boltz_models'
    cache.mkdir(parents=True)
    for name in ('ccd.pkl', 'boltz2_conf.ckpt', 'boltz2_aff.ckpt', 'mols.tar'):
        (cache / name).write_bytes(b'fixture-only; never loaded')
    monkeypatch.setattr(prediction, 'update_status', lambda *a, **k: None)
    containers, native, results = [], [], []
    for offset in range(0, n, 5):
        run = tmp_path / f'prediction-{offset}'; run.mkdir()
        jobs.write_json(run / 'metadata.json', {'resource_allocation': {'gpu_devices': ['0']}})
        monkeypatch.setattr(prediction, 'create_job', lambda *a, **k: SimpleNamespace(run_id=run.name, run_dir=run))
        chunk = source / f'prediction-{offset}.jsonl'
        chunk.write_text(''.join(json.dumps(r) + '\n' for r in rows[offset:offset+5]))
        request = graph.node_request('boltz2' if engine == 'boltz2' else 'prediction_protenix_v2',
                  dict(samples=2, seed_count=seeds, seed=42, generate_ipsae=False, use_target_template=False, target_use_msa=False),
                  dict(run_dir=str(source), candidates_jsonl=str(chunk)), 'fixture')
        active = {}
        def launch(command, **kwargs):
            if command[0] == 'docker':
                containers.append(command)
                manifest = host_path(run, command[command.index('--manifest') + 1])
                payload = json.loads(manifest.read_text())
                for item in payload['items']:
                    item['input'] = str(host_path(run, item['input']))
                fixture_manifest = run / 'fixture-manifest.json'; fixture_manifest.write_text(json.dumps(payload))
                native_request = contract.load_request(fixture_manifest, expected_engine=payload['engine'])
                output = host_path(run, command[command.index('--output-dir') + 1])
                active.update(request=native_request, output=output)
                item_rows, execution = adapter.run_batch(native_request, output)
                for row in item_rows:
                    row['outputs'] = ['/work/' + Path(value).relative_to(run).as_posix() for value in row.get('outputs', [])]
                (output / 'mn-batch-result.json').write_text(json.dumps({'items': item_rows, 'execution': execution}))
                return SimpleNamespace(returncode=0)
            native.append(command)
            r, output = active['request'], active['output']
            if engine == 'boltz2':
                out = Path(next(v.split('=', 1)[1] for v in command if v.startswith('--out_dir='))) if any(v.startswith('--out_dir=') for v in command) else Path(command[command.index('--out_dir') + 1])
                staged = output / '.mn-predict-boltz2-inputs'
                assert sorted(p.stem for p in staged.glob('*.yaml')) == sorted(i.item_id for i in r.items)
            else:
                staged = Path(command[command.index('-i')+1])
                assert sorted(p.stem for p in staged.glob('*.json')) == sorted(i.item_id for i in r.items)
                out = output
            for item in r.items:
                target = out / f'boltz_results_{staged.stem}' / 'predictions' / item.item_id if engine == 'boltz2' else out / item.item_id
                target.mkdir(parents=True, exist_ok=True)
                original = rows[offset + list(r.items).index(item)]
                mapping = {c: chr(65+i) for i, c in enumerate(original['target_chains'] + original['binder_chains'])}
                lines = Path(original['complex_pdb']).read_text().splitlines()
                text = '\n'.join(line[:21] + mapping.get(line[21], line[21]) + line[22:] if line.startswith(('ATOM', 'HETATM', 'TER   ')) and len(line)>21 else line for line in lines) + '\n'
                for sample in range(2): (target / f'{item.item_id}_{sample}.pdb').write_text(text)
            return SimpleNamespace(returncode=0)
        monkeypatch.setattr(prediction.subprocess, 'run', launch)
        prediction.run_optimized_prediction(**request['parameters'])
        result = jobs.read_json(run / 'result.json')
        assert result['success'], result
        results.extend(read_candidates(run))
        receipt = jobs.read_json(run / 'native_execution.json')
        assert receipt['docker_processes'] == 1
        assert receipt['native_invocations'] == seeds
        assert receipt['prediction_outputs'] == min(5, n-offset) * seeds * 2
    assert len(containers) == math.ceil(n/5)
    assert len(native) == math.ceil(n/5) * seeds
    assert len(results) == n * seeds * 2
    assert {p for r in results for p in r['parents']} == {r['candidate_id'] for r in rows}
    evidence(f'{engine}-seeds{seeds}-{n}', candidates=n, chunks=len(containers), docker_processes=len(containers),
             native_invocations=len(native), gpu_allocation_events=len(containers),
             gpu_allocation_basis='deterministic scheduler metadata fixture', prediction_outputs=len(results),
             exact_parent_ids=sorted({p for r in results for p in r['parents']}))


def test_read_only_plan_separates_seeds_chunks_and_concurrency():
    graph = candidate_pipeline.new_pipeline()
    predictor = next(n for n in graph['nodes'] if n['id'] == 'predictor-0')
    predictor['data']['params'].update(seed_count=3, samples=2)
    before = copy.deepcopy(graph)
    plan = {p['node']: p for p in candidate_pipeline.execution_plan(graph, {'predictor-0': 25})}
    assert graph == before
    row = plan['predictor-0']
    assert (row['candidates'], row['collection_chunks'], row['native_invocations'], row['docker_processes']) == (25, 5, 15, 5)
    assert (row['gpu_allocation_events'], row['gpu_slots'], row['prediction_outputs'], row['model_initializations']) == (5, 5, 150, 15)
    assert row['max_concurrent_jobs'] == 2
    assert plan['predictor-1']['native_invocations'] == 0
    assert plan['optional-redesign']['expected_jobs'] == plan['optional-redesign']['native_invocations'] == 0


@pytest.mark.parametrize('n', SIZES + [101])
def test_rfd4_runtime_launches_replicates_in_one_native_process_until_wrapper_limit(n, tmp_path):
    import ast
    import subprocess
    runtime = Path(sequence.__file__).with_name('rfd4_runtime.py')
    module = ast.parse(runtime.read_text())
    function = next(node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == 'run')
    launches = []
    context = dict(Path=Path, json=json, os=os,
        prepare=lambda request: [{}], compose_settings=lambda command, p, root: command,
        sampling_overrides=lambda p: [], metric_overrides=lambda p, **k: [], search_overrides=lambda p: [],
        subprocess=SimpleNamespace(run=lambda command, **kwargs: launches.append(command)))
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(runtime), 'exec'), context)
    with pytest.raises(RuntimeError, match='metrics table'):
        context['run']({'root': str(tmp_path), 'parameters': {'task': 'protein_binder', 'seed': 42, 'num_designs': n, 'batch_size': 2}})
    assert len(launches) == math.ceil(n/100)
    assert sum(int(next(v.split('=')[1] for v in cmd if v.startswith('num_replicates='))) for cmd in launches) == n
    assert all('batch_size=2' in cmd for cmd in launches)
    assert all(row['status'] == 'returned' for row in json.loads((tmp_path / 'native_command.json').read_text()))
    evidence(f'rfd4-{n}', requested_attempts=n, gpu_microbatch=2, native_invocations=len(launches),
             tested_boundary='actual RFD4 runtime function; deterministic native launch and metrics-table stop')


@pytest.mark.parametrize('n', SIZES)
def test_bindcraft_ordinary_trajectory_budget_uses_one_container(n, collection, monkeypatch, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import bindcraft2
    source, _, _ = collection(1)
    monkeypatch.setattr(jobs, 'runs_root', lambda: tmp_path / 'runs')
    monkeypatch.setattr(bindcraft2, 'missing_bindcraft2_parameters', lambda *a: [])
    monkeypatch.setattr(bindcraft2, 'bindcraft2_image_available', lambda *a: True)
    monkeypatch.setattr(bindcraft2, 'normalize_bindcraft2_candidates', lambda *a: [])
    launches = []
    def launch(command, **kwargs):
        launches.append(command)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(bindcraft2.subprocess, 'run', launch)
    run = bindcraft2.run_bindcraft2_campaign(parent_run_dir=tmp_path, target_pdb=source / 'input.pdb', target_chains=['A'],
        binder_length='60', number_of_final_designs=n, max_trajectories=n, workers_per_gpu=1)
    assert len(launches) == 1
    params = jobs.read_json(run / 'input.json')['params']
    assert params['max_trajectories'] == params['number_of_final_designs'] == n
    assert jobs.read_json(run / 'native_execution.json')['native_invocations'] == 1
    evidence(f'bindcraft2-{n}', requested_trajectory_limit=n, accepted_design_target=n, docker_processes=len(launches),
             native_invocations=1, model_initializations=None, worker_count=1,
             tested_boundary='ordinary BindCraft2 staging and Docker launch; native trajectory loop verified in source')


def test_native_launch_window_and_retry_accounting_stay_unknown(tmp_path):
    def crash(command, **kwargs): raise KeyboardInterrupt()
    with pytest.raises(KeyboardInterrupt):
        native_execution.record_process(tmp_path, ['docker', 'run', 'fixture'], native_invocations=1, runner=crash)
    receipt = jobs.read_json(tmp_path / 'native_execution.json')
    assert receipt['docker_processes'] is None and receipt['native_invocations'] is None
    summary = native_execution.execution_totals([dict(unit_count=5, candidate_count=3, attempts=[{'native_execution': receipt}, {}])])
    assert summary['candidates'] == 5 and summary['collection_chunks'] == 1 and summary['returned_records'] == 3
    assert summary['native_invocations'] is None and summary['docker_processes'] is None


@pytest.mark.parametrize('terminal', ['failed', 'lost', 'cancelled'])
def test_interrupted_mpnn_chunk_recovers_exact_completed_pairs_without_resubmission(collection, monkeypatch, tmp_path, terminal):
    from mn_molecular_design_studio import candidate_execution as ce, runtime
    from mn_molecular_design_studio.ligandmpnn_batching import recover
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR', str(tmp_path / 'runs'))
    source, path, rows = collection(5)
    run = tmp_path / 'runs/design/partial'; run.mkdir(parents=True)
    launches = []
    def launch(command, **kwargs):
        launches.append(command)
        paths = json.loads(host_path(run, command[command.index('--pdb_path_multi')+1]).read_text())
        output = host_path(run, command[command.index('--out_folder')+1])
        for folder in ('seqs', 'backbones'): (output / folder).mkdir(parents=True)
        for i, value in enumerate(paths[:3]):
            pdb = host_path(run, value)
            if i < 2: shutil.copy2(pdb, output / 'backbones' / f'{pdb.stem}_1.pdb')
            sequences = sequence.refolding_workflow._sequences_by_chain(pdb)
            (output / 'seqs' / f'{pdb.stem}.fa').write_text('>fixture, id=1\n' + ':'.join(sequences[c] for c in sorted(sequences)) + '\n')
        return SimpleNamespace(returncode=1)
    monkeypatch.setattr(sequence.subprocess, 'run', launch)
    sequence.run_ligandmpnn_sequence_design(source, path, existing_job=SimpleNamespace(run_dir=run))
    returned = read_candidates(run)
    assert [r['parents'][0] for r in returned] == [r['candidate_id'] for r in rows[:2]]
    expected_ids = [r['candidate_id'] for r in returned]
    (run / 'artifacts/normalized_candidates/candidates.jsonl').unlink()
    assert [r['candidate_id'] for r in recover(run)] == expected_ids
    jobs.write_json(run / 'metadata.json', {'status': terminal})
    monkeypatch.setattr(runtime, 'runs_root', lambda: tmp_path / 'runs')
    input_source = graph._local_source(tmp_path, 'source', rows)
    node = dict(id='sequence', data=dict(kind='ligandmpnn', label='sequence'))
    settings = ce.policy(dict(mode='sharded', batch_size=5, failure_policy='partial'))
    batches = ce.initialize(tmp_path, node, input_source, {'seed': 42, 'sequences': 1}, 'a'*64, settings,
        make_source=graph._local_source, make_request=graph.node_request, read_rows=graph._absolute_candidates)
    batches[0].update(job_run_dir=str(run), status='failed', attempts=[dict(index=0)])
    current = dict(batches=batches)
    monkeypatch.setattr(graph, '_submit', lambda *a, **k: pytest.fail('Completed members must not be resampled'))
    result, status = ce.tick(tmp_path, node, current, input_source, {}, 'a'*64, settings,
        graph=graph, owner_id='fixture', project_id=None, gpu_device='0', persist=lambda: None)
    assert status == 'partial'
    assert result['candidate_count'] == 2
    assert batches[0]['missing_candidate_ids'] == [r['candidate_id'] for r in rows[2:]]
    assert len(launches) == 1
    assert current['native_execution']['native_invocations'] == 1
    before = copy.deepcopy(current)
    ce.tick(tmp_path, node, current, input_source, {}, 'a'*64, settings,
        graph=graph, owner_id='fixture', project_id=None, gpu_device='0', persist=lambda: None)
    assert current == before


@pytest.mark.parametrize('engine', ['esm_if1', 'caliby', 'boltzgen_inverse'])
@pytest.mark.parametrize('n', SIZES)
def test_other_sequence_adapters_keep_multi_input_native_process(engine, n, collection, monkeypatch, tmp_path):
    from mn_molecular_design_studio import inverse_folding_input as inputs
    from mn_molecular_design_studio.engine.workflow_engine.workflows import inverse_folding as workflow
    source, path, rows = collection(n)
    run = tmp_path / 'inverse'; run.mkdir()
    monkeypatch.setattr(inputs, 'readiness', lambda *a: {'ready': True})
    monkeypatch.setattr(inputs, 'weights', lambda *a: [workflow.reference_root() / 'fixture.pt'])
    launches = []
    def launch(command, **kwargs):
        launches.append(command)
        payload = json.loads(host_path(run, command[-1]).read_text())
        assert len(payload['items']) == n
        outputs = []
        for item in payload['items']:
            sequences = {}
            for residue in item['mapping']:
                c = residue['engine_chain']
                sequences[c] = sequences.get(c, '') + residue['aa']
            for _ in range(payload['parameters']['num_sequences']):
                outputs.append(dict(key=item['key'], sequences=sequences, metrics={}))
        host_path(run, payload['output']).write_text(json.dumps(outputs))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(sequence.subprocess, 'run', launch)
    function = getattr(workflow, 'run_' + engine)
    settings = {'ensemble_mode': 'single'} if engine == 'caliby' else {}
    function(source_run_dir=source, candidates_jsonl=path, num_sequences=2, seed=42,
             existing_job=SimpleNamespace(run_dir=run), **settings)
    result = jobs.read_json(run / 'result.json')
    assert result['success'], result
    outputs = read_candidates(run)
    assert len(launches) == 1 and len(outputs) == 2*n
    assert {p for r in outputs for p in r['parents']} == {r['candidate_id'] for r in rows}
    receipt = jobs.read_json(run / 'native_execution.json')
    assert receipt['docker_processes'] == receipt['native_invocations'] == 1
    evidence(f'{engine}-{n}', candidates=n, chunks=1, docker_processes=len(launches), native_invocations=1,
             sequence_outputs=len(outputs), model_initializations_expected=n if engine == 'boltzgen_inverse' else 1,
             exact_parent_ids=sorted({p for r in outputs for p in r['parents']}))


def test_inventory_covers_every_category_and_rfd4_operations():
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, validate_graph
    inventory = native_execution.tool_inventory()
    keys = {r['tile'] for r in inventory}
    assert {k for k, s in NODE_CATALOG.items() if s['category'] in {'Design engine', 'Design pipeline', 'Sequence design', 'Structure prediction'}} <= keys
    assert 'rfd4/operation/inverse_folding' in keys
    rfd4_inverse = next(r for r in inventory if r['tile'] == 'rfd4/operation/inverse_folding')
    assert rfd4_inverse['category'] == 'Sequence design'
    assert rfd4_inverse['collection_socket'] is False  # honest integration limit
    for engine in ('ligandmpnn', 'abmpnn', 'esm_if1', 'caliby', 'boltzgen_inverse'):
        graph = candidate_pipeline.new_pipeline(sequence_designer=engine)
        validate_graph(graph)
        node = next(n for n in graph['nodes'] if n['id'] == 'sequence-preparation')
        assert node['data']['params']['sequence_stage'] == 'initial'
        assert any(r['tile'] == engine for r in inventory)


def test_mpnn_native_maps_preserve_distinct_constraints_and_global_groups(tmp_path):
    from mn_molecular_design_studio.ligandmpnn_batching import group_steps, materialize
    steps = []
    for i, source_id in enumerate(['c1', 'c1_2', 'other']):
        bias = tmp_path / f'bias-{i}.json'; bias.write_text(json.dumps({'Z2': {'G': i + 1}}))
        command = ['docker', 'run', 'fixture', 'python', '/opt/LigandMPNN/run.py', '--pdb_path', f'/work/{source_id}.pdb',
                   '--out_folder', f'/work/output/{source_id}', '--fixed_residues', f'Z{i+1}',
                   '--bias_AA_per_residue', '/work/' + bias.name, '--temperature', '.1' if i < 2 else '.5']
        steps.append(dict(command=command, name=source_id, candidate_ids=[source_id]))
    grouped, mappings = group_steps(tmp_path, steps)
    assert len(grouped) == 2
    assert grouped[0]['candidate_ids'] == ['c1', 'c1_2']
    command = grouped[0]['command']
    masks = json.loads(host_path(tmp_path, command[command.index('--fixed_residues_multi')+1]).read_text())
    biases = json.loads(host_path(tmp_path, command[command.index('--bias_AA_per_residue_multi')+1]).read_text())
    assert masks == {'/work/c1.pdb': 'Z1', '/work/c1_2.pdb': 'Z2'}
    assert biases == {'/work/c1.pdb': {'Z2': {'G': 1}}, '/work/c1_2.pdb': {'Z2': {'G': 2}}}
    output = Path(mappings[0]['output']) / 'backbones'; output.mkdir(parents=True)
    (output / 'c1_1.pdb').write_text('first')
    (output / 'c1_2_1.pdb').write_text('second')
    materialize(tmp_path, mappings)
    root = tmp_path / 'artifacts/raw/ligandmpnn/output'
    assert (root / 'c1/backbones/c1_1.pdb').read_text() == 'first'
    assert not (root / 'c1/backbones/c1_2_1.pdb').exists()
    assert (root / 'c1_2/backbones/c1_2_1.pdb').read_text() == 'second'


def test_gpu_allocations_are_not_counted_per_container(tmp_path):
    launcher = lambda command, **kwargs: SimpleNamespace(returncode=0)
    jobs.write_json(tmp_path / 'metadata.json', {'resource_allocation': {'gpu_devices': ['0'], 'allocated_at': 'first'}})
    for _ in range(2): native_execution.record_process(tmp_path, ['docker', 'run', 'fixture'], native_invocations=1, runner=launcher)
    value = jobs.read_json(tmp_path / 'native_execution.json')
    assert value['docker_processes'] == 2 and value['gpu_allocation_events'] == 1 and value['gpu_slots'] == 1
    jobs.write_json(tmp_path / 'metadata.json', {'resource_allocation': {'gpu_devices': ['0'], 'allocated_at': 'retry'}})
    native_execution.record_process(tmp_path, ['docker', 'run', 'fixture'], native_invocations=1, runner=launcher)
    value = jobs.read_json(tmp_path / 'native_execution.json')
    assert value['docker_processes'] == 3 and value['gpu_allocation_events'] == 2 and value['gpu_slots'] == 2


def test_legacy_unconfirmed_docker_job_is_not_blindly_resampled(collection, monkeypatch, tmp_path):
    from mn_molecular_design_studio import candidate_execution as ce, runtime
    source, path, rows = collection(1)
    input_source = graph._local_source(tmp_path, 'source', rows)
    node = dict(id='sequence', data=dict(kind='ligandmpnn', label='sequence'))
    settings = ce.policy(dict(mode='sharded', batch_size=5, retries=1))
    batches = ce.initialize(tmp_path, node, input_source, {'sequences': 1, 'seed': 42}, 'a'*64, settings,
        make_source=graph._local_source, make_request=graph.node_request, read_rows=graph._absolute_candidates)
    run = tmp_path / 'runs/design/legacy'; run.mkdir(parents=True)
    jobs.write_json(run / 'metadata.json', {'status': 'failed'})
    jobs.write_json(run / 'command.json', {'mode': 'docker', 'command': ['docker', 'run', 'legacy']})
    batches[0].update(job_run_dir=str(run), status='failed', attempts=[{}])
    monkeypatch.setattr(runtime, 'runs_root', lambda: tmp_path / 'runs')
    monkeypatch.setattr(graph, '_submit', lambda *a, **k: pytest.fail('Unconfirmed native work cannot be duplicated'))
    current = dict(batches=batches)
    _, status = ce.tick(tmp_path, node, current, input_source, {}, 'a'*64, settings,
        graph=graph, owner_id='fixture', project_id=None, gpu_device='0', persist=lambda: None)
    assert status == 'failed' and batches[0]['quarantined'] is True
    assert current['native_execution']['native_invocations'] is None
