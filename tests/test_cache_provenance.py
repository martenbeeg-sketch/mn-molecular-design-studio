"""F9 cache regressions through the production graph scheduler."""
import copy
import json

import pytest
from mn_molecular_design_studio import pipeline_graph as graph

from test_pipeline_graph import execution_harness
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner


def test_cached_output_mutation_is_detected(execution_harness):
    g, calls, _ = execution_harness
    first = runner.run_graph_pipeline(g, 'first', 'owner', 'project')
    state = json.loads((first / 'graph_execution.json').read_text())['nodes']
    output = state['foundry_mpnn-3']['output']
    assert runner._cache_valid(output)
    path = runner._path(output['candidates_jsonl'])
    path.write_text(path.read_text().replace('KLMNP', 'KLMNA'))
    assert not runner._cache_valid(output)
    runner.run_graph_pipeline(g, 'changed', 'owner', 'project')
    assert calls == ['rfd3', 'foundry_mpnn', 'rf3', 'foundry_mpnn']


def test_fresh_sampling_and_reuse_have_recorded_resolved_seeds(execution_harness):
    g, calls, _ = execution_harness
    first = runner.run_graph_pipeline(g, 'first', 'owner', 'project')
    second = runner.run_graph_pipeline(g, 'reuse', 'owner', 'project')
    a = json.loads((first / 'graph_execution.json').read_text())['nodes']
    b = json.loads((second / 'graph_execution.json').read_text())['nodes']
    assert a['rfd3-2']['resolved_seed'] == b['rfd3-2']['resolved_seed'] > 0
    assert b['rfd3-2']['status'] == 'cached'
    fresh = runner.run_graph_pipeline(g, 'fresh', 'owner', 'project', reuse_completed=False)
    c = json.loads((fresh / 'graph_execution.json').read_text())['nodes']
    assert c['rfd3-2']['resolved_seed'] != a['rfd3-2']['resolved_seed']
    assert len(calls) == 6


def test_interface_setting_invalidates_only_sequence_design_and_descendants(execution_harness):
    g, calls, _ = execution_harness
    runner.run_graph_pipeline(g, 'first', 'owner', 'project')
    changed = copy.deepcopy(g)
    mpnn = next(n for n in changed['nodes'] if n['data']['kind'] == 'foundry_mpnn')
    mpnn['data']['params'].update(fixed_residues='B1', interface_cutoff=6.0)
    runner.run_graph_pipeline(changed, 'changed', 'owner', 'project')
    assert calls == ['rfd3', 'foundry_mpnn', 'rf3', 'foundry_mpnn', 'rf3']


def test_changed_runtime_checkpoint_invalidates_mpnn_but_reuses_generation(execution_harness, monkeypatch, tmp_path):
    from mn_molecular_design_studio import mpnn_input
    weights = tmp_path / 'checkpoint.pt'
    weights.write_bytes(b'checkpoint one')
    monkeypatch.setattr(mpnn_input, 'checkpoint', lambda *a, **k: weights)
    g, calls, _ = execution_harness
    runner.run_graph_pipeline(g, 'first', 'owner', 'project')
    weights.write_bytes(b'checkpoint two')
    runner.run_graph_pipeline(g, 'changed weights', 'owner', 'project')
    assert calls == ['rfd3', 'foundry_mpnn', 'rf3', 'foundry_mpnn', 'rf3']


def test_changed_structure_contents_invalidates_output_manifest(execution_harness):
    g, calls, _ = execution_harness
    first = runner.run_graph_pipeline(g, 'first', 'owner', 'project')
    state = json.loads((first / 'graph_execution.json').read_text())['nodes']
    output = state['foundry_mpnn-3']['output']
    structure = runner._path(runner._absolute_candidates(output)[0]['complex_pdb'])
    assert runner._cache_valid(output)
    structure.write_text(structure.read_text() + '\n# changed bytes\n')
    assert not runner._cache_valid(output)


@pytest.mark.parametrize('kind,key', [('cofold_input', 'sequence_file'), ('cofold_input', 'compound_set_file'), ('rfd4', 'input_json')])
def test_fingerprint_tracks_referenced_file_contents(tmp_path, kind, key):
    path = tmp_path / 'input.txt'
    path.write_text('first')
    before = graph.fingerprint(kind, {key: str(path)}, [])
    path.write_text('other')
    assert graph.fingerprint(kind, {key: str(path)}, []) != before


def test_nested_condition_file_tracks_relative_structure_contents(tmp_path):
    structure = tmp_path / 'target.pdb'
    structure.write_text('first')
    conditions = tmp_path / 'conditions.yaml'
    conditions.write_text('global_args:\n  input: target.pdb\n')
    before = graph.fingerprint('rfd4', {'input_files': [str(conditions)]}, [])
    structure.write_text('other')
    assert graph.fingerprint('rfd4', {'input_files': [str(conditions)]}, []) != before


def test_resumed_graph_does_not_publish_old_children_under_new_fingerprint(execution_harness, monkeypatch):
    from types import SimpleNamespace
    g, calls, _ = execution_harness
    first = runner.run_graph_pipeline(g, 'first', 'owner', 'project')
    before = json.loads((first / 'graph_execution.json').read_text())['nodes']
    old_target = runner._path(before['target-0']['output']['target']['target_pdb'])
    original_bytes = old_target.read_bytes()
    target = runner._path(g['nodes'][0]['data']['params']['target_pdb'])
    target.write_bytes(target.read_bytes() + b'\n# new source content\n')
    monkeypatch.setattr(runner, 'create_job', lambda *a, **k: SimpleNamespace(run_dir=first))
    runner.run_graph_pipeline(g, 'first', 'owner', 'project')
    after = json.loads((first / 'graph_execution.json').read_text())['nodes']
    assert calls == ['rfd3', 'foundry_mpnn', 'rf3', 'rfd3', 'foundry_mpnn', 'rf3']
    assert len(after['rfd3-2']['superseded_attempts']) == 1
    assert old_target.read_bytes() == original_bytes
