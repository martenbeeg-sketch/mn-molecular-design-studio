"""Small real-container integrations. No mocks and no numerical reproducibility claims."""
import json
import math
from pathlib import Path

import pytest

from mn_molecular_design_studio import design_validation as dv
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as predictor
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as graph
from mn_molecular_design_studio.molecule_imports import parse_structure

REPO = Path(__file__).resolve().parents[1]
NANOBODY = REPO / 'mn_molecular_design_studio/assets/generators/protpardelle/motifs/nanobody/7eow_AB_atom.pdb'


def reference_candidate(root):
    """An existing experimental reference encoded with the RFD4 candidate schema."""
    root.mkdir(parents=True, exist_ok=True)
    from Bio.PDB import MMCIFIO
    structure = dv.read_structure(NANOBODY)
    # Exercise source IDs that differ from the predictor's canonical A/B IDs.
    model = next(structure.get_models()); model['A'].id = 'T'; model['B'].id = 'V'
    path = root / 'reference.cif'
    writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(path))
    sequences = {key: dv.sequence(value) for key, value in dv.protein_chains(structure).items()}
    candidate = {'candidate_id': '7eow-reference', 'source_tool': 'rfd4_proteina',
        'stage': 'generation.backbone_sequence', 'complex_pdb': str(path),
        'target_chains': ['T'], 'binder_chains': ['V'], 'binder_sequence': sequences['V'],
        'raw_metadata': {'studio_molecule_id': '7eow-reference',
                         'reference_origin': 'experimental structure fixture; not an RFD4-generated sample',
                         'sequences_by_chain': sequences,
                         'components': parse_structure(path.read_bytes(), path.name)[0]['components']}}
    write_candidates(root, 'rfd4_proteina', [candidate])
    return candidate


def predict_source(runtime, source_run):
    runtime['require_image'](predictor.ENGINE_PROFILES['boltz2']['image'])
    for name in ('boltz2_conf.ckpt', 'boltz2_aff.ckpt', 'ccd.pkl'):
        if not (runtime['reference'] / 'boltz_models' / name).is_file():
            pytest.skip(f'installed Boltz checkpoint/chemistry cache unavailable: {name}')
    tasks = graph._optimized_prediction_tasks('boltz2', {
        'run_dir': str(source_run), 'candidates_jsonl': 'artifacts/normalized_candidates/candidates.jsonl',
    }, {'target_use_msa': False, 'use_target_template': False})
    run = predictor.run_optimized_prediction(
        'boltz2', tasks=tasks, samples=1, seed=20261027, num_recycles=1,
        gpu_device=runtime['device'], cpu_cores=4, use_target_template=False, target_use_msa=False,
        engine_options={'use_msa': False, 'use_msa_server': False, 'sampling_steps': 20, 'generate_ipsae': False},
        job_name='scientific-contract-reference',
    )
    result = read_json(run / 'result.json')
    assert result.get('success') is True, (run, result)
    candidates = read_candidates(run)
    assert len(candidates) == len(tasks)
    by_source = {task['source_candidate_id']: task for task in tasks}
    for candidate in candidates:
        task = by_source[candidate['parents'][0]]
        normalized = predictor._prepare_prediction('boltz2', tasks=[task], engine_options={'use_msa': False}, resolve_random_seed=False).normalized_tasks[0]
        from mn_molecular_design_studio.scientific_contracts import verify_prediction_structure
        verify_prediction_structure(candidate['complex_pdb'], normalized)
        assert candidate['target_chains'] == normalized['target_protein_ids']
        assert candidate['binder_chains'] == normalized['binder_protein_ids']
        assert candidate['binder_sequence'] == normalized['binder_sequence']
        assert candidate['raw_metadata']['identity_check']['sequence_chain_identity'] == 'pass'
        assert math.isfinite(candidate['metrics']['iptm'])
        assert math.isfinite(candidate['metrics']['plddt_mean'])
    return run, candidates


@pytest.mark.gpu_e2e
def test_rfd4_candidate_boltz2_interface_analysis(gpu_runtime):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis
    gpu_runtime['require_image'](analysis.IPSAE_IMAGE)
    gpu_runtime['require_image'](analysis.PYROSETTA_IMAGE)
    source_run = gpu_runtime['root'] / 'source'
    original = reference_candidate(source_run)
    run, candidates = predict_source(gpu_runtime, source_run)
    assert candidates[0]['parents'] == [original['candidate_id']]
    analysis_run = analysis.run_ipsae_analysis_contract(
        run, run / 'artifacts/normalized_candidates/candidates.jsonl', pae_cutoff=30., dist_cutoff=20.,
    )
    rows = read_candidates(analysis_run)
    assert rows and read_json(analysis_run / 'result.json')['success'], (analysis_run, rows)
    assert not rows[0]['metrics'].get('ipsae_error'), rows
    assert math.isfinite(rows[0]['metrics']['ipsae'])
    provenance = rows[0]['raw_metadata']['measurement_provenance']
    assert provenance['source_candidate_id'] == candidates[0]['candidate_id']
    assert provenance['binder_chains'] == ['B'] and provenance['target_chains'] == ['A']
    rosetta = analysis.run_rosetta_metrics_contract(run, run / 'artifacts/normalized_candidates/candidates.jsonl', nprocs=2)
    scores = read_candidates(rosetta)
    assert scores and read_json(rosetta / 'result.json')['success'], (rosetta, scores)
    assert not scores[0]['metrics'].get('rosetta_error'), scores


@pytest.mark.gpu_e2e
def test_real_rfd4_generation_boltz2(gpu_runtime):
    from mn_molecular_design_studio import rfd4_proteina_input as inputs
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import run_rfd4_proteina_design
    gpu_runtime['require_image'](inputs.IMAGE)
    snapshot = gpu_runtime['reference'] / 'rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots' / inputs.HF_REVISION
    if not (snapshot / 'model.safetensors').is_file():
        pytest.skip('installed RFD4 checkpoint unavailable')
    target = REPO / 'mn_molecular_design_studio/assets/proteina_complexa/PD-L1.pdb'
    condition_spec = '''global_args:
  input: PD-L1.pdb
binder:
  contig: "30-30,/0,A1-20"
  conditions:
    C_SEQ:
      true:
        - select: "sel('A')"
    C_CRDGRP:
      1:
        - select: "sel('A/*/*/[N,CA,C,O]')"
'''
    run = run_rfd4_proteina_design(
        task='custom_conditions', input_mode='condition_spec', input_files=[target],
        condition_spec=condition_spec, num_designs=1, batch_size=1, nsteps=60, seed=20261027,
        options={'refold': False, 'filter_conditions': False}, cpu_cores=4, gpu_device=gpu_runtime['device'],
    )
    assert read_json(run / 'result.json').get('success') is True, run
    generated = read_candidates(run)
    assert len(generated) == 1
    assert generated[0]['target_chains'] and generated[0]['binder_chains']
    assert len(generated[0]['binder_sequence']) == 30
    target_sequence = dv.sequence(dv.protein_chains(dv.read_structure(target))['A'])[:20]
    sequences = generated[0]['raw_metadata']['sequences_by_chain']
    assert [sequences[chain] for chain in generated[0]['target_chains']] == [target_sequence]
    prediction, rows = predict_source(gpu_runtime, run)
    assert rows[0]['parents'] == [generated[0]['candidate_id']]


@pytest.mark.gpu
def test_rfd4_multichain_ccd_components_boltz2(gpu_runtime):
    """Reference-sequence fragments plus installed CCD chemistry exercise heterogeneous mapping."""
    import pickle
    import numpy as np
    from Bio.PDB import Atom, Chain, Residue, MMCIFIO
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import sequence_structure
    reference = dv.protein_chains(dv.read_structure(NANOBODY))
    segments = [dv.sequence(reference['B'])[:30], dv.sequence(reference['A'])[:12], dv.sequence(reference['A'])[12:24]]
    source = gpu_runtime['root'] / 'heterogeneous-source'; source.mkdir()
    path = sequence_structure(':'.join(segments), source / 'source.cif')
    structure = dv.read_structure(path); model = next(structure.get_models())
    for chain, name in zip(list(model), ['V', 'T', 'U']):
        chain.id = name
    chemistry = Chain.Chain('L'); model.add(chemistry)
    for index, code in enumerate(['ATP', 'ZN'], 1):
        ccd = gpu_runtime['reference'] / 'boltz_models/mols' / (code + '.pkl')
        if not ccd.is_file():
            pytest.skip(f'installed CCD unavailable: {code}')
        molecule = pickle.loads(ccd.read_bytes())
        residue = Residue.Residue(('H_' + code, index, ' '), code, ' '); chemistry.add(residue)
        for i, atom in enumerate(molecule.GetAtoms()):
            if atom.GetAtomicNum() == 1 or (atom.HasProp('leaving_atom') and atom.GetBoolProp('leaving_atom')):
                continue
            coordinate = list(molecule.GetConformer().GetAtomPosition(i)) if molecule.GetNumConformers() else [i * 1.5, 0., 0.]
            name = atom.GetProp('name')
            residue.add(Atom.Atom(name, np.array(coordinate), 0., 1., ' ', name, i + 1, element=atom.GetSymbol().upper()))
    writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(path))
    candidate = {'candidate_id': 'heterogeneous-reference', 'source_tool': 'rfd4_proteina',
        'stage': 'generation.backbone_sequence', 'complex_pdb': str(path),
        'target_chains': ['U', 'T'], 'binder_chains': ['V'], 'binder_sequence': segments[0],
        'raw_metadata': {'studio_molecule_id': 'heterogeneous-reference',
                         'components': parse_structure(path.read_bytes(), path.name)[0]['components']}}
    write_candidates(source, 'rfd4_proteina', [candidate])
    run, rows = predict_source(gpu_runtime, source)
    row = rows[0]
    assert row['target_chains'] == ['A', 'B'] and row['binder_chains'] == ['C']
    contract = row['raw_metadata']['scientific_contract']
    assert contract['chain_map'] == {'U': 'A', 'T': 'B', 'V': 'C', 'L_1': 'D', 'L_2': 'E'}
    assert {item.get('ccd') for item in contract['components'] if item['type'] == 'small_molecule'} == {'ATP', 'ZN'}
    assert all(item.get('heavy_atoms') for item in contract['components'] if item['type'] == 'small_molecule')
