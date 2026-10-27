"""Scientific identity contracts, independent of predictor numerical results."""
import copy
import json

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as graph
from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import _prepare_prediction
from mn_molecular_design_studio.pipeline_graph import new_graph, validate_graph


def component_source(design_source_factory):
    run, candidate = design_source_factory(missing_occupancy=False)
    candidate.update(source_tool="rfd4_proteina", stage="generation.backbone_sequence")
    candidate["raw_metadata"] = {
        "studio_molecule_id": "rfd4-reference", "sequences_by_chain": {"A": "ACDEFG", "B": "KLMNP"},
        "components": [{"type": "protein", "chain": "B", "sequence": "KLMNP"},
                       {"type": "protein", "chain": "A", "sequence": "ACDEFG"}],
    }
    path = run / "artifacts/normalized_candidates/candidates.jsonl"
    path.write_text(json.dumps(candidate) + "\n")
    return {"run_dir": str(run), "candidates_jsonl": str(path)}


def test_rfd4_component_roles_and_sequences_survive_normalization(design_source_factory):
    source = component_source(design_source_factory)
    tasks = graph._optimized_prediction_tasks("boltz2", source, {"target_use_msa": True})
    task = _prepare_prediction("boltz2", tasks=tasks, resolve_random_seed=False).normalized_tasks[0]
    assert task["target_protein_ids"] == ["A"]
    assert task["binder_protein_ids"] == ["B"]
    assert task["binder_sequence"] == "KLMNP"
    assert task["proteins"] == [{"id": "A", "sequence": "ACDEFG"}, {"id": "B", "sequence": "KLMNP"}]
    assert task["molecule_chain_options"]["A"]["use_msa"] is True
    assert task["molecule_chain_options"]["B"]["use_msa"] is False
    assert task["source_candidate_id"] == "passing-design"
    assert task["scientific_contract"]["chain_map"] == {"A": "A", "B": "B"}


def test_rfd4_declared_sequence_cannot_disagree_with_structure(design_source_factory):
    source = component_source(design_source_factory)
    from pathlib import Path
    path = Path(source["candidates_jsonl"])
    candidate = json.loads(path.read_text())
    candidate["raw_metadata"]["components"][0]["sequence"] = "AAAAA"
    path.write_text(json.dumps(candidate) + "\n")
    with pytest.raises(ValueError, match="sequence"):
        graph._optimized_prediction_tasks("boltz2", source, {})


def test_analysis_missing_roles_are_not_assumed():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _chain_role_groups
    assert _chain_role_groups({}) == ([], [])



def test_multichain_target_nontrivial_order_and_msa_mapping(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import sequence_structure
    from mn_molecular_design_studio.scientific_contracts import component_prediction_task, verify_prediction_structure
    path = sequence_structure('KLMNP:ACDEFG:QRST', tmp_path / 'source.cif')
    from mn_molecular_design_studio import design_validation as dv
    from Bio.PDB import MMCIFIO
    structure = dv.read_structure(path)
    for chain, name in zip(next(structure.get_models()), ['Z', 'X', 'Y']):
        chain.id = name
    writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(path))
    candidate = {'candidate_id': 'multi', 'complex_pdb': str(path), 'binder_sequence': 'KLMNP',
                 'target_chains': ['Y', 'X'], 'binder_chains': ['Z'],
                 'raw_metadata': {'components': [
                     {'type': 'protein', 'chain': key, 'sequence': seq}
                     for key, seq in [('Z', 'KLMNP'), ('X', 'ACDEFG'), ('Y', 'QRST')]],
                     'fixed_residue_identities': {'Z1': 'K'},
                     'molecule_chain_options': {'X': {'use_msa': True, 'msa_a3m_text': '>query\nACDEFG\n'}}}}
    raw = component_prediction_task(candidate, tmp_path)
    task = _prepare_prediction('boltz2', tasks=[raw], resolve_random_seed=False).normalized_tasks[0]
    assert task['target_protein_ids'] == ['A', 'B']
    assert task['binder_protein_ids'] == ['C']
    assert task['scientific_contract']['chain_map'] == {'Y': 'A', 'X': 'B', 'Z': 'C'}
    assert task['molecule_chain_options']['B']['msa_a3m_text'].strip() == '>query\nACDEFG'
    assert task['molecule_chain_options']['C']['use_msa'] is False
    assert task['scientific_contract']['reference_metadata']['fixed_residue_identities'] == {'Z1': 'K'}
    predicted = sequence_structure('QRST:ACDEFG:KLMNP', tmp_path / 'prediction.cif')
    assert verify_prediction_structure(predicted, task)['sequence_chain_identity'] == 'pass'
    wrong = sequence_structure('QRST:ACDEFG:AAAAA', tmp_path / 'wrong.cif')
    with pytest.raises(ValueError, match='mismatch'):
        verify_prediction_structure(wrong, task)


def test_single_atom_ion_and_ccd_ligand_mapping(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import sequence_structure
    from mn_molecular_design_studio import design_validation as dv
    from mn_molecular_design_studio.molecule_imports import parse_structure
    from mn_molecular_design_studio.scientific_contracts import component_prediction_task
    from Bio.PDB import Atom, Chain, Residue, MMCIFIO
    import numpy as np
    path = sequence_structure('ACDEFG', tmp_path / 'complex.cif')
    structure = dv.read_structure(path)
    chain = Chain.Chain('Z'); next(structure.get_models()).add(chain)
    for i, code in enumerate(['ZN', 'ATP'], 1):
        residue = Residue.Residue(('H_' + code, i, ' '), code, ' '); chain.add(residue)
        residue.add(Atom.Atom('ZN' if i == 1 else 'P', np.zeros(3), 0., 1., ' ', 'ZN' if i == 1 else ' P  ', i, element='ZN' if i == 1 else 'P'))
    writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(path))
    components = parse_structure(path.read_bytes(), path.name)[0]['components']
    assert {row.get('ccd') for row in components if row['type'] == 'small_molecule'} == {'ZN', 'ATP'}
    raw = component_prediction_task({'candidate_id': 'chemistry', 'complex_pdb': str(path),
        'binder_chains': ['A'], 'binder_sequence': 'ACDEFG', 'raw_metadata': {'components': components}}, tmp_path)
    task = _prepare_prediction('boltz2', tasks=[raw], resolve_random_seed=False).normalized_tasks[0]
    assert [row['code'] for row in task['ligand_ccd_codes']] == ['ZN', 'ATP']
    assert set(task['scientific_contract']['chain_map']) == {'A', 'Z_1', 'Z_2'}


def test_unsupported_bonds_fail_explicitly(design_source_factory):
    from pathlib import Path
    source = component_source(design_source_factory)
    path = Path(source['candidates_jsonl']); row = json.loads(path.read_text())
    row['raw_metadata']['bonds'] = [{'atom1': 'A1:SG', 'atom2': 'B1:C'}]
    path.write_text(json.dumps(row) + '\n')
    with pytest.raises(ValueError, match='Unsupported chemistry.*bonds'):
        graph._optimized_prediction_tasks('boltz2', source, {})


@pytest.mark.parametrize('value', [None, float('nan'), float('inf')])
def test_missing_or_nonfinite_analysis_never_passes(value):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _passes_thresholds
    assert not _passes_thresholds({'ipsae': value}, {'min_ipsae': 0.0})[0]


def test_analysis_error_overrides_inherited_native_acceptance():
    from mn_molecular_design_studio.design_validation import candidate_passes_filters
    assert not candidate_passes_filters({'metrics': {'native_pass_filters': True, 'ipsae_error': 'PAE missing'}})


def test_analysis_without_measurement_finishes_unsuccessfully(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis
    captured = []
    monkeypatch.setattr(analysis, 'finish_job', lambda run, success, result: captured.append((success, result)))
    analysis._finish_metric_analysis(SimpleNamespace(run_dir=tmp_path),
        [{'candidate_id': 'no-pae', 'metrics': {'ipsae_error': 'PAE missing'}}],
        tool='ipsae', table_name='metrics.csv', metric_names=['ipsae'])
    assert captured[0][0] is False
    assert captured[0][1]['summary']['scored_count'] == 0


def test_rosetta_reports_multichain_targets_as_unsupported(design_source_factory, tmp_path, monkeypatch):
    from Bio.PDB import MMCIFIO
    from mn_molecular_design_studio import design_validation as dv
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates, read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    run, row = design_source_factory(missing_occupancy=False)
    path = run / 'designed.cif'; structure = dv.read_structure(path)
    model = next(structure.get_models()); chain = model['A'].copy(); chain.id = 'C'; model.add(chain)
    writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(path))
    row['target_chains'] = ['C', 'A']
    write_candidates(run, 'reference', [row])
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR', str(tmp_path / 'runs'))
    scored = analysis.run_rosetta_metrics_contract(run, run / 'artifacts/normalized_candidates/candidates.jsonl')
    assert read_json(scored / 'result.json')['success'] is False
    output = read_candidates(scored)[0]
    assert 'one binder chain and one target' in output['metrics']['rosetta_error']
    assert output['raw_metadata']['measurement_provenance']['status'] == 'error'


def test_generic_ligand_labels_restore_only_with_verified_ccd_atoms(tmp_path):
    from mn_molecular_design_studio.scientific_contracts import restore_prediction_ccd_names
    path = tmp_path / 'native.pdb'
    path.write_text('HETATM    1 ZN   LIG A   1       0.000   0.000   0.000  1.00  0.00          ZN  \nEND\n')
    task = {'scientific_contract': {'chain_map': {'metal': 'A'}, 'components': [
        {'type': 'small_molecule', 'chain': 'metal', 'ccd': 'ZN', 'heavy_atoms': {'ZN': 'ZN'}}]}}
    restored = restore_prediction_ccd_names(path, task)
    assert restored != path
    assert 'LIG' in path.read_text() and ' ZN ' in restored.read_text()
    task['scientific_contract']['components'][0]['heavy_atoms'] = {'FE': 'FE'}
    with pytest.raises(ValueError, match='atom identities'):
        restore_prediction_ccd_names(path, task)
