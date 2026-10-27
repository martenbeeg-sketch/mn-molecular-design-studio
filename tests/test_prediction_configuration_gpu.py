"""Small real AF2 and Boltz runs with native feature/manifest evidence retained."""
import json
import math
from pathlib import Path

import pytest

from mn_molecular_design_studio import design_validation as comparison
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

REPO = Path(__file__).resolve().parents[1]
REFERENCE = REPO / 'mn_molecular_design_studio/assets/generators/protpardelle/motifs/nanobody/7eow_AB_atom.pdb'


def small_complex(root):
    from Bio.PDB import PDBIO
    root.mkdir(parents=True, exist_ok=True)
    structure = comparison.read_structure(REFERENCE)
    model = next(structure.get_models())
    for chain in model:
        for residue in list(chain)[24:]:
            chain.detach_child(residue.id)
    model['A'].id = 'T'
    model['B'].id = 'V'
    path = root / 'reference.pdb'
    writer = PDBIO(); writer.set_structure(structure); writer.save(str(path))
    seqs = {chain: comparison.sequence(value) for chain, value in comparison.protein_chains(structure).items()}
    candidate = {'candidate_id': '7eow-fragments', 'stage': 'generation.backbone_sequence',
                 'source_tool': 'reference_fixture', 'complex_pdb': str(path), 'target_pdb': str(path),
                 'binder_sequence': seqs['V'], 'target_chains': ['T'], 'binder_chains': ['V'],
                 'raw_metadata': {'reference_origin': '24-residue fragments from 7EOW; execution test, not binding-quality benchmark'}}
    write_candidates(root, 'reference_fixture', [candidate])
    return candidate


@pytest.mark.gpu
@pytest.mark.parametrize('mode', ['conditioned_self_consistency', 'target_only'])
def test_real_af2_conditioning_modes(gpu_runtime, mode):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import refolding
    gpu_runtime['require_image']('mn-colabdesign:latest')
    if not (refolding.ALPHAFOLD_MODELS_DIR / 'params_model_1_multimer_v3.npz').is_file():
        pytest.skip('Installed AF2 multimer parameters unavailable')
    source = gpu_runtime['root'] / 'source'
    original = small_complex(source)
    run = refolding.run_af2_initial_guess_complex_refolding(
        source, source / 'artifacts/normalized_candidates/candidates.jsonl', require_monomer_success=False,
        num_recycles=1, model_count=1, multimer=True, af2_validation_mode=mode,
        msa_source='single_sequence', gpu_device=gpu_runtime['device'], seed=20261027)
    (gpu_runtime['root'] / 'prediction_run.txt').write_text(str(run))
    assert read_json(run / 'result.json')['success'], run
    rows = read_candidates(run)
    assert len(rows) == 1
    row = rows[0]; raw = row['raw_metadata']; native = raw['effective_settings']
    assert native['use_initial_guess'] == (mode == 'conditioned_self_consistency')
    assert native['redesign'] == (mode == 'conditioned_self_consistency')
    assert native['binder_coordinate_atom_count'] > 0 if mode == 'conditioned_self_consistency' else native['binder_coordinate_atom_count'] == 0
    assert native['binder_template_residue_count'] == 0
    assert native['binder_sequence'] == original['binder_sequence']
    assert native['num_recycles'] == 1 and native['model_seeds'] == [20261027]
    assert raw['conditioning_mode'] == mode
    staged = run / raw['prediction_input_pdb']
    assert set(refolding._structure_chains(staged)) == ({'T', 'V'} if mode == 'conditioned_self_consistency' else {'T'})
    assert set(refolding._structure_chains(run / row['complex_pdb'])) == {'T', 'V'}
    assert row['target_chains'] == ['T'] and row['binder_chains'] == ['V']
    assert refolding._pdb_sequences_by_chain(run / row['complex_pdb'])['V'] == original['binder_sequence']
    assert math.isfinite(row['metrics']['binder_rmsd_target_aligned'])
    assert math.isfinite(row['metrics']['iptm'])
    assert (run / raw['reference_structure_path']).is_file()
    assert raw['validation_chain_map'] == {'T': 'T', 'V': 'V'}


@pytest.mark.gpu
@pytest.mark.parametrize('samples,recycles,pae', [(1, 1, True), (2, 2, True)])
def test_real_boltz_sampling_and_recycling(gpu_runtime, samples, recycles, pae):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as predictor, graph_pipeline
    gpu_runtime['require_image'](predictor.ENGINE_PROFILES['boltz2']['image'])
    for name in ('boltz2_conf.ckpt', 'boltz2_aff.ckpt', 'ccd.pkl'):
        if not (gpu_runtime['reference'] / 'boltz_models' / name).is_file():
            pytest.skip(f'Installed Boltz cache missing: {name}')
    source = gpu_runtime['root'] / 'source'; original = small_complex(source)
    tasks = graph_pipeline._optimized_prediction_tasks('boltz2', {
        'run_dir': str(source), 'candidates_jsonl': 'artifacts/normalized_candidates/candidates.jsonl',
    }, {'target_use_msa': False, 'use_target_template': False})
    request = predictor.boltz2_prediction_parameters(tasks, settings={
        'diffusion_samples': samples, 'recycling_steps': recycles, 'seed': 20261027,
        'use_msa': False, 'use_msa_server': False, 'sampling_steps': 20,
        'generate_ipsae': False, 'write_full_pae': pae, 'gpu_device': gpu_runtime['device'],
    }, job_name='boltz-effective-configuration')
    run = predictor.run_optimized_prediction(**request, use_target_template=False, target_use_msa=False)
    (gpu_runtime['root'] / 'prediction_run.txt').write_text(str(run))
    assert read_json(run / 'result.json')['success'], run
    native = read_json(run / 'artifacts/raw/optimized_prediction/output/mn-batch-result.json')
    assert native['samples_per_item'] == samples
    assert native['execution']['return_code'] == 0
    for command in native['execution']['commands']:
        assert command[command.index('--diffusion_samples') + 1] == str(samples)
        assert command[command.index('--recycling_steps') + 1] == str(recycles)
        assert ('--write_full_pae' in command) == pae
    rows = read_candidates(run)
    assert len(rows) == samples
    for row in rows:
        raw = row['raw_metadata']
        assert raw['effective_settings']['samples_per_item'] == samples
        assert raw['effective_settings']['recycling_steps'] == recycles
        assert raw['effective_settings']['write_full_pae'] == pae
        assert bool(raw['pae_path']) == pae
        assert row['parents'] == [original['candidate_id']]
        assert row['binder_sequence'] == original['binder_sequence']
        assert math.isfinite(row['metrics']['iptm'])
        assert raw['identity_check']['sequence_chain_identity'] == 'pass'
        assert raw['native_execution']['commands'] == native['execution']['commands']
