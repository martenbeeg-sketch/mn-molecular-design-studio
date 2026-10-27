"""Real representative generators and both unified MPNN stages; installed models only."""
import json
import os
from pathlib import Path

import pytest

from mn_molecular_design_studio import design_validation as dv, mpnn_input
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.workflows import mpnn, sequence_design as native

REPO = Path(__file__).resolve().parents[1]


@pytest.mark.gpu_e2e
@pytest.mark.parametrize('engine', ['rfdiffusion', 'rfd4'])
def test_real_two_stage_sequence_design(gpu_runtime, engine):
    runtime = gpu_runtime
    runtime['require_image'](mpnn_input.IMAGE)
    for model in ('protein_mpnn', 'soluble_mpnn'):
        if not mpnn_input.checkpoint(model).is_file():
            pytest.skip(f'installed {model} checkpoint unavailable')
    target = REPO / 'mn_molecular_design_studio/assets/proteina_complexa/PD-L1.pdb'
    if engine == 'rfdiffusion':
        from mn_molecular_design_studio.engine.workflow_engine.workflows.generators import run_rfdiffusion_generation
        runtime['require_image']('mn-rfdiffusion:latest')
        if not (runtime['reference'] / 'rfdiffusion_models/Complex_base_ckpt.pt').is_file():
            pytest.skip('installed RFdiffusion complex checkpoint unavailable')
        from Bio.PDB import PDBIO, Select
        class Fragment(Select):
            def accept_residue(self, residue):
                return residue.parent.id == 'A' and 1 <= residue.id[1] <= 40
        cropped = runtime['root'] / 'target-fragment.pdb'
        writer = PDBIO(); writer.set_structure(dv.read_structure(target)); writer.save(str(cropped), Fragment())
        generation = Path(os.environ['MDS_PHASE4_RFDIFFUSION_GENERATION']) if os.environ.get('MDS_PHASE4_RFDIFFUSION_GENERATION') else run_rfdiffusion_generation(target_pdb=cropped, target_chains=['A'], binder_length='60',
            task='binder', num_designs=1, timesteps=15,
            seed=20261027, cpu_cores=4, gpu_device=runtime['device'])
    else:
        from mn_molecular_design_studio import rfd4_proteina_input as inputs
        from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import run_rfd4_proteina_design
        runtime['require_image'](inputs.IMAGE)
        snapshot = runtime['reference'] / 'rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots' / inputs.HF_REVISION
        if not (snapshot / 'model.safetensors').is_file():
            pytest.skip('installed RFD4 checkpoint unavailable')
        conditions = '''global_args:
  input: PD-L1.pdb
binder:
  contig: "60-60,/0,A1-40"
  conditions:
    C_SEQ:
      true:
        - select: "sel('A')"
    C_CRDGRP:
      1:
        - select: "sel('A/*/*/[N,CA,C,O]')"
'''
        generation = Path(os.environ['MDS_PHASE4_RFD4_GENERATION']) if os.environ.get('MDS_PHASE4_RFD4_GENERATION') else run_rfd4_proteina_design(task='custom_conditions', input_mode='condition_spec',
            input_files=[target], condition_spec=conditions, num_designs=1, batch_size=1, nsteps=60,
            seed=20261027, options={'refold': False, 'filter_conditions': False}, cpu_cores=4, gpu_device=runtime['device'])
    assert read_json(generation / 'result.json')['success'] is True, generation
    generated = read_candidates(generation)
    assert len(generated) == 1
    initial = mpnn.run_mpnn_sequence_design(generation, generation / 'artifacts/normalized_candidates/candidates.jsonl',
        sequence_stage='initial', num_seq_per_target=1, backend='shared_ligandmpnn', seed=20261027, gpu_device=runtime['device'])
    assert read_json(initial / 'result.json')['success'] is True, initial
    parent = read_candidates(initial)[0]
    if engine == 'rfd4':
        assert parent['binder_sequence'] == generated[0]['binder_sequence']
        assert not (initial / 'artifacts/raw/ligandmpnn').exists()
    else:
        assert parent['source_tool'] == 'protein_mpnn'
        assert parent['parents'] == [generated[0]['candidate_id']]
    path = Path(parent['complex_pdb'])
    chains = dv.protein_chains(dv.read_structure(path))
    binder = parent['binder_chains']
    labels = [f'{c}{r.id[1]}{r.id[2].strip()}' for c in binder for r in chains[c]]
    contacts = native._interface_contacts_for_structure(path, binder, parent['target_chains'], cutoff=4.0)
    # Also exercise an independently fixed motif away from the geometric interface.
    free = [label for label, aa in zip(labels, parent['binder_sequence']) if label not in contacts and aa != 'C']
    assert len(free) > 1
    motif = free[0]
    parent.setdefault('raw_metadata', {})['fixed_residues'] = [motif]
    source = runtime['root'] / 'optimization-input'
    write_candidates(source, 'sequence_stage_control', [parent])
    variants_run = mpnn.run_mpnn_sequence_design(source, source / 'artifacts/normalized_candidates/candidates.jsonl',
        sequence_stage='soluble_optimization', optimization_enabled=True, num_seq_per_target=2, sampling_temp=.4,
        seed=20261027, gpu_device=runtime['device'])
    assert read_json(variants_run / 'result.json')['success'] is True, variants_run
    rows = read_candidates(variants_run)
    assert len(rows) == 3
    assert rows[0]['candidate_id'] == parent['candidate_id']
    assert rows[0]['binder_sequence'] == parent['binder_sequence']
    fixed = set(contacts) | {motif} | {label for label, aa in zip(labels, parent['binder_sequence']) if aa == 'C'}
    changes = []
    for child in rows[1:]:
        assert child['parents'] == [parent['candidate_id']]
        assert child['source_tool'] == 'soluble_mpnn'
        assert child['raw_metadata']['fixed_sequence_validation']['status'] == 'verified'
        changed = {label for label, before, after in zip(labels, parent['binder_sequence'], child['binder_sequence']) if before != after}
        assert not changed & fixed
        changes.append(len(changed))
        actual = dv.protein_chains(dv.read_structure(Path(child['complex_pdb'])))
        assert child['binder_sequence'] == ''.join(dv.sequence(actual[c]) for c in child['binder_chains'])
        assert [dv.sequence(chains[c]) for c in parent['target_chains']] == [dv.sequence(actual[c]) for c in child['target_chains']]
        if engine == 'rfd4':
            from mn_molecular_design_studio.scientific_contracts import component_prediction_task
            task = component_prediction_task(child, variants_run, target_use_msa=False)
            assert task['scientific_contract']['binder_sequence'] == child['binder_sequence']
    if engine == 'rfd4':
        from test_scientific_gpu import predict_source
        from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline
        prediction, predicted = predict_source(runtime, variants_run)
        assert {r['binder_sequence'] for r in predicted} == {r['binder_sequence'] for r in rows}
        output = graph_pipeline._attach_refold_input_references(
            {'run_dir': str(prediction), 'candidates_jsonl': str(prediction / 'artifacts/normalized_candidates/candidates.jsonl')},
            {'run_dir': str(variants_run), 'candidates_jsonl': str(variants_run / 'artifacts/normalized_candidates/candidates.jsonl')},
            runtime['root'], 'prediction')
        results = graph_pipeline._local_source(runtime['root'], 'results', graph_pipeline._absolute_candidates(output), final_evidence=True)
        assert len(read_candidates(Path(results['run_dir']))) == 3
    assert any(changes), 'Allowed non-interface positions should change under real sampling'
    assert contacts, 'Reference should exercise actual interface contacts'
    evidence = {'engine': engine, 'generation': str(generation), 'initial': str(initial), 'optimization': str(variants_run),
        'interface_definition': {'atoms': 'heavy', 'cutoff_angstrom': 4.0}, 'interface_residues': contacts,
        'fixed_motif': motif, 'changed_non_interface_counts': changes, 'parent_sequence': parent['binder_sequence'],
        'variant_sequences': [r['binder_sequence'] for r in rows[1:]], 'seed': 20261027}
    (runtime['root'] / 'execution-evidence.json').write_text(json.dumps(evidence, indent=2))
