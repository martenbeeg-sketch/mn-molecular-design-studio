"""Independent Complexa generation and safe downstream adapters; no inference."""
import inspect
import json
import shutil
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio import complexa_input as inputs, complexa_ui as ui
from mn_molecular_design_studio import generation_io, pipeline_graph as graph
from mn_molecular_design_studio.complexa_container import ordered_motifs
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation
from mn_molecular_design_studio.engine.workflow_engine.workflows import design, graph_pipeline as runner, refolding, sequence_design


@pytest.mark.parametrize('task', list(inputs.TASKS))
def test_generator_needs_only_own_weights_and_disables_downstream_models(task, tmp_path):
    parameters = inputs.prepare_request(ui.example_parameters(task) | {'generator_only': True})
    parameters['run_name'] = 'generator-test'
    config = inputs.native_overrides(parameters, 'A1-115')
    assert config['generation']['reward_model'] is None
    assert config['generation']['search']['algorithm'] == 'single-pass'
    assert config['generation']['search']['best_of_n']['replicas'] == 1
    for name in inputs.CHECKPOINTS[task]:
        path = inputs.resource_paths(tmp_path)['checkpoints'] / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'weights')
    assert inputs.readiness(task, tmp_path, generator_only=True)['ready']
    assert not inputs.readiness(task, tmp_path)['ready']


def test_motif_order_survives_sorted_json_and_matches_atom_selections():
    parameters = inputs.prepare_request(ui.example_parameters('motif_ligand'))
    parameters['run_name'] = 'motif-test'
    roundtrip = json.loads(json.dumps(parameters, sort_keys=True))
    assert ordered_motifs(roundtrip) == list(parameters['motif_atoms'])
    custom = inputs.native_overrides(roundtrip, '')['generation']['motif_target_dict_cfg']['MN_APP_TARGET']
    for i, original in enumerate(ordered_motifs(roundtrip), 1):
        assert f"B{i}: [{', '.join(parameters['motif_atoms'][original])}]" in custom['contig_atoms']
    with pytest.raises(ValueError, match='motif residue order'):
        ordered_motifs(roundtrip | {'motif_residue_order': ['B64', 'B64']})


def test_generator_schema_and_runner_do_not_mount_folding_dependencies(monkeypatch, tmp_path):
    refs = tmp_path / 'refs'
    for name in inputs.CHECKPOINTS['protein_binder']:
        path = inputs.resource_paths(refs)['checkpoints'] / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'weights')
    run = tmp_path / 'run'
    run.mkdir()
    jobs, commands = [], []
    def create(*args):
        jobs.append(args)
        return SimpleNamespace(run_dir=run)
    monkeypatch.setattr(design, 'reference_root', lambda: refs)
    monkeypatch.setattr(design, 'create_job', create)
    monkeypatch.setattr(design, 'update_status', lambda *a, **k: None)
    monkeypatch.setattr(design, '_finish_design_job', lambda *a, **k: None)
    monkeypatch.setattr(design, '_run_shell_steps', lambda run, steps: commands.extend(steps) or 1)
    spec = automation.WORKFLOWS['design.proteina_complexa_generate']
    parameters = ui.example_parameters('protein_binder')
    allowed = inspect.signature(design.run_proteina_complexa_generation).parameters
    design.run_proteina_complexa_generation(**{k: v for k, v in parameters.items() if k in allowed})
    assert jobs[0][1:3] == ('design_generation', 'proteina_complexa_generator')
    command = commands[0]['command']
    mounts = [command[i+1] for i, arg in enumerate(command[:-1]) if arg == '-v']
    assert mounts and not any(any(x in mount for x in ['/rf3.ckpt', '/model_params', '/mpnn_weights', '/studio-tools']) for mount in mounts)
    config = json.loads((run / 'artifacts/raw/proteina_complexa/input/target_config.json').read_text())
    assert config['generator_only'] is True
    assert config['overrides']['generation']['reward_model'] is None
    fields = {p['name'] for p in automation.workflow_schema(spec.workflow_id)['parameters']}
    assert {'pipeline_task', 'ligand', 'motif_atom_spec'} <= fields
    assert not fields & {'generator_only', 'search_strategy', 'replicas'}
    automation._validate_workflow_inputs(spec, parameters)


def test_generated_outputs_are_separate_from_successful_design_validation(tmp_path):
    from mn_molecular_design_studio.design_validation import design_candidates
    folder = tmp_path / 'artifacts/raw/proteina_complexa/inference/test/job_0'
    folder.mkdir(parents=True)
    source = inputs.ASSETS / 'M0024_1nzy_v3.pdb'
    for name in ['sample.pdb', 'duplicate.pdb']:
        shutil.copy(source, folder / name)
        (folder / name).with_suffix('.studio.json').write_text(json.dumps({
            'fixed_residues': ['B64'], 'fixed_motif_required': True, 'motif_mapping_verified': True}))
    parameters = inputs.prepare_request(ui.example_parameters('motif_ligand') | {'generator_only': True})
    rows = design._normalize_proteina_complexa_candidates(tmp_path, parameters, source)
    assert len(rows) == 1
    assert rows[0]['metrics']['pass_filters'] is None
    assert rows[0]['raw_metadata']['fixed_residues'] == ['B64']
    assert rows[0]['raw_metadata']['target_entity_kind'] == 'ligand'
    for name, value in {
        'metadata.json': {'status': 'completed', 'tool': 'proteina_complexa_generator'},
        'input.json': {'params': parameters}, 'result.json': {'success': True},
        'molecular_design_studio.json': {'project_id': 'project'},
    }.items():
        (tmp_path / name).write_text(json.dumps(value))
    assert len(generation_io.generated_candidates(tmp_path)) == 1
    assert design_candidates(tmp_path) == []
    assert runner.source_generations(str(tmp_path), 'project')[1][0]['candidate_id'] == rows[0]['candidate_id']
    for project, selection in [('other', None), ('project', []), ('project', ['unknown'])]:
        with pytest.raises(ValueError):
            runner.source_generations(str(tmp_path), project, selection)
    (tmp_path / 'result.json').write_text('{"success": false}')
    assert generation_io.generated_candidates(tmp_path) == []


def test_generator_template_is_independent_and_accepts_direct_prediction():
    g = graph.new_graph('complexa')
    assert [n['data']['kind'] for n in g['nodes']] == ['complexa_generator', 'ligandmpnn', 'rf3', 'results']
    graph.validate_graph(g)
    with pytest.raises(ValueError, match='Load a structure'):
        graph.validate_graph(g, require_inputs=True)
    allowed = {f['name'] for f in graph.NODE_CATALOG['complexa_generator']['fields']}
    g['nodes'][0]['data']['params'].update({k: v for k, v in ui.example_parameters('protein_binder').items() if k in allowed})
    g['nodes'][0]['data']['params']['target_chains'] = 'A'
    graph.validate_graph(g, require_inputs=True)
    source = {'target': {'target_pdb': 'inherited.pdb', 'target_chains': ['A'], 'hotspots': 'A37'}}
    p = g['nodes'][0]['data']['params'] | {'hotspots': ''}
    request = runner.node_request('complexa_generator', p, source, 'test')
    assert request['parameters']['target_pdb'] == p['target_pdb']
    assert request['parameters']['hotspots'] == ''
    inspect.signature(design.run_proteina_complexa_generation).bind(**request['parameters'])
    # Remove sequence redesign: the native generator already supplies a sequence.
    g['nodes'].pop(1)
    g['edges'] = [dict(g['edges'][1], source=g['nodes'][0]['id']), g['edges'][2]]
    graph.validate_graph(g)
    g['nodes'][1]['data']['kind'] = 'af2'
    g['nodes'][1]['data']['params'] = {f['name']: f['default'] for f in graph.NODE_CATALOG['af2']['fields']}
    g['nodes'][0]['data']['params']['pipeline_task'] = 'ligand_binder'
    with pytest.raises(ValueError, match='does not predict ligand'):
        graph.validate_graph(g)


def atom(serial, name, residue, chain, index, *, hetero=False):
    return f"{'HETATM' if hetero else 'ATOM  '}{serial:5d} {name:>4s} {residue:>3s} {chain}{index:4d}    {float(serial):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{90.:6.2f}          {name[0]:>2s}"


@pytest.fixture
def ligand_candidate(tmp_path):
    source = tmp_path / 'complex.pdb'
    source.write_text('\n'.join([
        atom(1, 'N', 'ALA', 'B', 1), atom(2, 'CA', 'ALA', 'B', 1),
        atom(3, 'C', 'ALA', 'B', 1), atom(4, 'O', 'ALA', 'B', 1),
        atom(10, 'C1', 'L:0', 'A', 1, hetero=True), atom(11, 'O1', 'L:0', 'A', 1, hetero=True),
        'CONECT   10   11   11', 'CONECT   11   10   10', 'END', '']))
    return {'candidate_id': 'generated', 'complex_pdb': str(source), 'binder_sequence': 'G',
            'binder_chains': ['B'], 'target_chains': ['A'], 'stage': 'generation.backbone_sequence',
            'raw_metadata': {'target_entity_kind': 'ligand', 'ligand_original_code': 'IAI',
                             'fixed_motif_required': True, 'motif_mapping_verified': True, 'fixed_residues': ['B1']}}


def test_mpnn_preserves_ligand_graph_and_remaps_fixed_motif(ligand_candidate, tmp_path):
    source = tmp_path / 'complex.pdb'
    staged = tmp_path / 'staged.pdb'
    original = source.read_bytes()
    roles = sequence_design._stage_ligandmpnn_role_normalized_input(
        source_pdb=source, staged_pdb=staged, candidate=ligand_candidate, requested_design_chains='B')
    assert source.read_bytes() == original
    lines = staged.read_text().splitlines()
    ligand = [line for line in lines if line.startswith('HETATM')]
    assert len(ligand) == 2
    serials = [int(line[6:11]) for line in ligand]
    bonds = [generation_io._connections(line) for line in lines if line.startswith('CONECT')]
    assert [serials[0], serials[1], serials[1]] in bonds  # Double bond preserved.
    fixed = sequence_design._fixed_residues_for_candidate(ligand_candidate, [], roles)
    assert fixed == [roles['binder_chains'][0] + '1']
    ligand_candidate['raw_metadata']['motif_mapping_verified'] = False
    with pytest.raises(ValueError, match='verified mapping'):
        sequence_design._fixed_residues_for_candidate(ligand_candidate, [], roles)


def test_rf3_receives_redesigned_sequence_and_exact_ligand_graph(ligand_candidate, tmp_path):
    input_dir = tmp_path / 'prediction/artifacts/raw/rf3/inputs'
    staged = refolding._stage_complex_inputs(tmp_path, [ligand_candidate], input_dir)
    paths, metrics = refolding._write_rf3_json_inputs(
        source_run_dir=tmp_path, source_candidates=[ligand_candidate], staged=staged,
        input_dir=input_dir, benchmark_run_csv=None, use_target_msa=False, use_target_template=True)
    payload = json.loads(paths['generated'].read_text())
    assert payload['components'][0] == {'seq': 'G', 'chain_id': 'B'}
    assert metrics['rf3_ligand_component_count'] == 1
    assert 'template_selection' not in payload  # Ligand coordinates are not a protein template.
    ligand_path = tmp_path / 'prediction' / payload['components'][1]['path'].removeprefix('/work/')
    text = ligand_path.read_text()
    assert text.count('HETATM') == 2 and 'CONECT   10   11   11' in text


@pytest.mark.parametrize("with_run_dir", [False, True])
def test_ligand_prediction_adapters_accept_only_supported_chemistry(ligand_candidate, monkeypatch, tmp_path, with_run_dir):
    source = {'target': {'pipeline_task': 'ligand_binder'}}
    if with_run_dir:
        source['run_dir'] = str(tmp_path)
    monkeypatch.setattr(runner, '_absolute_candidates', lambda _: [ligand_candidate])
    params = {f['name']: f['default'] for f in graph.NODE_CATALOG['boltz2']['fields']}
    request = runner.node_request('boltz2', params, source, 'ligand-test')
    assert request['parameters']['tasks'][0]['ligand_ccd_codes'] == 'IAI'
    assert request['parameters']['tasks'][0]['binder_sequence'] == 'G'
    assert request['parameters']['target_use_msa'] is False
    assert request['parameters']['engine_options']['use_msa'] is False
    automation._validate_workflow_inputs(automation.WORKFLOWS[request['workflow']], request['parameters'])
    ligand_candidate['raw_metadata']['ligand_original_code'] = 'L:1'
    with pytest.raises(ValueError, match='custom ligand fragment'):
        runner.node_request('boltz2', params, source, 'fragment-test')
    with pytest.raises(ValueError, match='does not predict ligand'):
        runner.node_request('af2', {}, source, 'af2-test')


def test_generator_shared_form_supports_all_tasks_without_submission():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
from mn_molecular_design_studio.complexa_ui import render_setup
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
render_setup(workflow_schema("design.proteina_complexa_generate"), owner_id="test", project_selector=lambda _:None, submit=lambda *a:(_ for _ in ()).throw(AssertionError("No job expected")))
''').run()
    assert not app.exception
    for task in inputs.TASKS:
        app.button_group(key='complexa-generator-task').set_value(task).run()
        assert not app.exception
        assert any(b.label == 'Submit 1 job' and b.disabled for b in app.button)
        assert not any('search' in s.label.lower() for s in app.selectbox)


def test_generator_node_uses_shared_task_form():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
from mn_molecular_design_studio.pipeline_graph import new_graph
from mn_molecular_design_studio.pipeline_node_ui import render_node_form
g=new_graph('complexa')
render_node_form(g['nodes'][0], g, key='generator-node', owner_id='test', render_input=lambda:{})
''').run()
    assert not app.exception
    assert app.button_group(key='generator-node-complexa-task')
    assert not any(n.label == 'Number of jobs' for n in app.number_input)


def test_api_validates_generator_without_separate_target_node(monkeypatch, tmp_path):
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio import api
    actor = {'id': 'test-user', 'role': 'user'}
    monkeypatch.setattr(api, '_request_payload_and_project', lambda payload, user: (
        {k: v for k, v in payload.items() if k != 'project_id'}, payload['project_id']))
    uploaded = tmp_path / 'inputs'
    uploaded.mkdir()
    shutil.copy2(inputs.ASSETS / 'PD-L1.pdb', uploaded / 'target.pdb')
    monkeypatch.setattr(api, 'user_input_root', lambda _: uploaded)
    monkeypatch.setattr(api, 'reference_root', lambda: tmp_path / 'refs')
    monkeypatch.setattr(api, 'path_is_within_accessible_job', lambda *a, **k: False)
    monkeypatch.setattr(api, 'accessible_project_ids', lambda _: [])
    monkeypatch.setattr(runner, 'enqueue_graph_pipeline', lambda *a, **k: pytest.fail('Validation must not enqueue'))
    g = graph.new_graph('complexa')
    g['nodes'][0]['data']['params'].update(target_pdb=str(uploaded / 'target.pdb'), target_chains='A')
    api.app.dependency_overrides[api._current_user] = lambda: actor
    try:
        response = TestClient(api.app).post('/api/v1/pipelines/validate', json={'project_id': 'project', 'graph': g})
        assert response.status_code == 200, response.text
        assert response.json()['valid'] is True
    finally:
        api.app.dependency_overrides.clear()
