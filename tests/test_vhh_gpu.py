"""Opt-in antibody-specific redesign and native VHH generation, using installed models."""
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import design_validation as dv, mpnn_input, sabr_runtime
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json


@pytest.mark.gpu
def test_real_abmpnn_preserves_framework_cysteines_and_unselected_cdrs(gpu_runtime):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.mpnn import run_mpnn_sequence_design
    gpu_runtime['require_image'](mpnn_input.IMAGE)
    if not mpnn_input.checkpoint('ab_mpnn').is_file():
        pytest.skip('installed AbMPNN checkpoint unavailable')
    source = mpnn_input.ASSETS / '9F6O.pdb'
    run = run_mpnn_sequence_design(
        structure_path=source, model_type='ab_mpnn', design_chains='B',
        auto_detect_cdrs=True, cdr_regions=['CDR3'], num_seq_per_target=1,
        omit_aas='', seed=20261027, gpu_device=gpu_runtime['device'],
    )
    assert read_json(run / 'result.json').get('success') is True, run
    rows = read_candidates(run)
    assert len(rows) == 1
    reference = dv.protein_chains(dv.read_structure(source))['B']
    annotation = sabr_runtime.annotate_file(source, ['B'])['chains']['B']['cdr_residues']
    fixed, mutable = sabr_runtime.fixed_residues_from_annotation(source, ['B'], ['CDR3'], annotation)
    labels = [sabr_runtime.residue_label('B', r) for r in reference]
    sequence = rows[0]['binder_sequence']
    original = dv.sequence(reference)
    assert len(sequence) == len(original)
    for label, before, after in zip(labels, original, sequence):
        if label in fixed or before == 'C':
            assert before == after, label
    assert rows[0]['raw_metadata']['fixed_sequence_validation']['status'] == 'verified'
    actual = dv.protein_chains(dv.read_structure(Path(rows[0]['complex_pdb'])))
    assert sequence == ''.join(dv.sequence(actual[chain]) for chain in rows[0]['binder_chains'])


@pytest.mark.gpu_e2e
def test_real_esmfold2_vhh_generation_boltz2(gpu_runtime):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import esm_binder
    from mn_molecular_design_studio.antibody_constraints import verify_masked_sequence
    from test_scientific_gpu import predict_source
    gpu_runtime['require_image'](esm_binder.ESMFOLD2_BINDER_IMAGE)
    model = esm_binder.ESMFOLD2_BINDER_MODEL_ROOT / esm_binder.DEFAULT_BINDER_MODEL
    if not (model / 'config.json').is_file() or not (esm_binder.BIOHUB_ESM_ROOT / 'ESMC-6B').is_dir():
        pytest.skip('installed native ESMFold2 VHH/ESMC checkpoints unavailable')
    # A short real PD-L1 sequence fragment limits GPU cost; no binding-quality claim.
    target = Path(__file__).resolve().parents[1] / 'mn_molecular_design_studio/assets/proteina_complexa/PD-L1.pdb'
    target_sequence = dv.sequence(dv.protein_chains(dv.read_structure(target))['A'])[:24]
    run = esm_binder.run_esmfold2_native_binder_design(
        target_sequence=target_sequence, binder_type='vhh', binder_framework='ozoralizumab_tnf_vhh',
        num_designs=1, optimization_steps=1, seed=20261027, gpu_device=gpu_runtime['device'],
        compile_model=False, checkpoint_lm=True,
    )
    result = read_json(run / 'result.json')
    assert result.get('success') is True, (run, result)
    rows = read_candidates(run)
    assert len(rows) == 1
    row = rows[0]
    check = verify_masked_sequence(row['binder_sequence'], row['raw_metadata']['sampled_binder_prompt'])
    assert check['status'] == 'verified'
    # Real antibody numbering must recognize the emitted VHH.
    sabr_runtime.annotate_file(Path(row['complex_pdb']), row['binder_chains'])
    prediction, predicted = predict_source(gpu_runtime, run)
    assert predicted[0]['binder_sequence'] == row['binder_sequence']
    assert predicted[0]['parents'] == [row['candidate_id']]
