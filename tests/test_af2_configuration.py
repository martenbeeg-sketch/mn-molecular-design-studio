"""F4/F6 regressions at graph, request, staging and native invocation boundaries."""
import copy
from pathlib import Path

import pytest

from mn_molecular_design_studio import pipeline_graph as catalog
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as graph
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as prediction


def test_af2_modes_are_forwarded_without_forcing_initial_guess(tmp_path):
    source = {'run_dir': str(tmp_path), 'candidates_jsonl': 'candidates.jsonl'}
    requests = [graph.node_request('af2', {'af2_validation_mode': mode}, source, 'validation')['parameters']
                for mode in ('conditioned_self_consistency', 'target_only')]
    assert requests[0]['use_initial_guess'] is True
    assert requests[1]['use_initial_guess'] is False


def test_saved_af2_without_mode_keeps_conditioned_behavior(tmp_path):
    request = graph.node_request('af2', {}, {'run_dir': str(tmp_path), 'candidates_jsonl': 'candidates.jsonl'}, 'old')
    assert request['parameters']['use_initial_guess'] is True
    assert request['parameters']['af2_validation_mode'] == 'conditioned_self_consistency'



def test_target_only_staging_excludes_binder_coordinates_and_keeps_reference(design_source_factory, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import refolding
    source, candidate = design_source_factory()
    original = copy.deepcopy(candidate)
    conditioned = refolding._stage_complex_inputs(source, [candidate], tmp_path / 'conditioned')
    target = refolding._stage_target_template_inputs(source, [candidate], tmp_path / 'target')
    a, b = conditioned[candidate['candidate_id']], target[candidate['candidate_id']]
    assert set(refolding._structure_chains(a)) == {'A', 'B'}
    assert set(refolding._structure_chains(b)) == {'A'}
    assert refolding._pdb_sequences_by_chain(b) == {'A': 'ACDEFG'}
    assert b.with_suffix('.binder_sequence.txt').read_text().strip() == 'KLMNP'
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    chain_map = read_json(b.with_suffix('.chain_map.json'))
    assert chain_map['binder'] == [{'original_chain': 'B', 'engine_chain': 'B', 'role': 'binder'}]
    assert chain_map['targets'][0]['original_chain'] == 'A'
    assert candidate == original


@pytest.mark.parametrize('mode,initial', [('target_only', True), ('conditioned_self_consistency', False)])
def test_af2_conflicting_conditioning_is_rejected(mode, initial):
    from mn_molecular_design_studio.prediction_configuration import af2_mode
    with pytest.raises(ValueError, match='Conflicting'):
        af2_mode(mode, initial)


def test_target_plus_sequence_needs_no_designed_binder_coordinates(design_source_factory, tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import refolding
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    source, candidate = design_source_factory()
    pdb = tmp_path / 'target.pdb'
    refolding._cif_to_pdb(source / candidate['complex_pdb'], pdb)
    direct = {'candidate_id': 'target-sequence', 'stage': 'complex_refolding',
              'target_pdb': str(pdb), 'target_chains': ['A'], 'binder_chains': ['Z'], 'binder_sequence': 'KLMNP'}
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR', str(tmp_path / 'runs'))
    monkeypatch.setattr(refolding, '_run_shell_steps', lambda *args, **kwargs: 1)
    run = refolding.run_af2_initial_guess_complex_refolding(source, None, require_monomer_success=False,
        direct_candidates=[direct], af2_validation_mode='target_only')
    staged = run / 'artifacts/raw/af2_initial_guess/inputs/target-sequence.pdb'
    assert refolding._pdb_sequences_by_chain(staged) == {'A': 'ACDEFG'}
    assert not list(staged.parent.joinpath('references').glob('*.pdb'))
    assert read_json(run / 'result.json')['success'] is False


def test_direct_af2_request_forwards_visible_settings(tmp_path):
    target = tmp_path / 'target.pdb'; target.write_text('END\n')
    source = {'molecule_task': {'target_pdb': str(target), 'partner_protein_sequences': 'ACDEFG'}}
    p = graph.node_request('af2', {'af2_validation_mode': 'target_only', 'seed': 42,
        'msa_source': 'reference_then_colabfold', 'msa_max_depth': 128}, source, 'target')['parameters']
    assert (p['seed'], p['msa_source'], p['msa_max_depth']) == (42, 'reference_then_colabfold', 128)


@pytest.mark.parametrize('template', ['rfdiffusion', 'genie3', 'pxdesign'])
def test_af2_dependent_starters_preserve_conditioned_default(template):
    graph_data = catalog.new_graph(template)
    checked = catalog.clean_graph(graph_data)
    af2 = next(n for n in checked['nodes'] if n['data']['kind'] == 'af2')
    assert af2['data']['params']['af2_validation_mode'] == 'conditioned_self_consistency'


def test_af2_editor_has_unique_controls():
    fields = [f['name'] for f in catalog.NODE_CATALOG['af2']['fields']]
    assert len(fields) == len(set(fields))
