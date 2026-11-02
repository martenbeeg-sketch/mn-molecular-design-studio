import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import cofolding_results, localfold
from mn_molecular_design_studio import cofolding_input
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as prediction


def prepare(model='mol1', **kwargs):
    options = {'model': model, 'use_msa': False, **kwargs.pop('engine_options', {})}
    return prediction._prepare_prediction(
        engine='localfold', protein_sequences='ACDEFGHIKLMNPQRSTVWY', samples=1,
        seed=42, accept_alphafold3_terms=True, engine_options=options, **kwargs)


@pytest.mark.parametrize('model', localfold.MODELS)
def test_models_validate_through_the_actual_prediction_boundary(model):
    value = prepare(model, engine_options={'accept_helixfold3_terms': True})
    assert value.selected_options['model'] == model
    assert value.profile['image'] == localfold.IMAGE
    assert value.requested_gpu_count == 1
    assert value.recycle_count == localfold.MODELS[model]['recycles']


@pytest.mark.parametrize('model', localfold.MODELS)
def test_model_parameter_defaults_match_the_pinned_localfold_registry(model):
    value = prediction._prepare_prediction(
        engine='localfold', protein_sequences='ACDEFGHIKLMNPQRSTVWY', samples=1,
        seed=42, engine_options={'model': model, 'accept_helixfold3_terms': True},
        accept_alphafold3_terms=True,
        resolve_random_seed=False)
    defaults = localfold.MODELS[model]
    assert value.recycle_count == defaults['recycles']
    assert value.selected_options['num_sampling_steps'] == (defaults['diffusion_steps'] or 0)
    assert value.selected_options['max_msa_seqs'] == defaults['msa_depth']
    assert value.selected_options['msa_mode'] == (
        'single_sequence' if model in {'esmfold2_lm600m', 'esmfold2_lm300m'}
        else 'automatic')
    options = prediction._batch_options(
        value.profile, seed=42, num_recycles=value.recycle_count, engine_options=value.selected_options)
    assert options['num_msa'] == (defaults['msa_depth'] if value.selected_options['use_msa'] else 1)
    assert 'msa_mode' not in options and 'max_msa_seqs' not in options


def test_official_weight_terms_are_conditional():
    prediction._prepare_prediction(engine='localfold', protein_sequences='ACDEFG', samples=1,
                                  engine_options={'model': 'mol1', 'use_msa': False})
    with pytest.raises(ValueError, match='accept AlphaFold3'):
        prediction._prepare_prediction(engine='localfold', protein_sequences='ACDEFG', samples=1)


def test_omitted_sample_count_uses_model_default_for_af2():
    value = prediction._prepare_prediction(
        engine='localfold', protein_sequences='ACDEFG', seed=42,
        engine_options={'model': 'af2_ptm'}, resolve_random_seed=False)
    assert value.profile['default_samples'] == 1


@pytest.mark.parametrize('precision', ['int8', 'fp16', 'fp32'])
def test_precision_is_passed_to_the_image_without_implicit_fallback(precision):
    value = prepare(engine_options={'weights_precision': precision})
    options = prediction._batch_options(value.profile, seed=42, num_recycles=value.recycle_count,
                                        engine_options=value.selected_options)
    assert options['weights_precision'] == precision
    assert options['num_recycles'] == 3
    assert 'accept_helixfold3_terms' not in options


def test_localfold_uses_its_own_shared_reference_tree_without_copying_old_weights(tmp_path):
    old = tmp_path / 'mn-cofolding/cache/alphafold3/weights'
    for name in ['mol1-int8', 'esmfold2_lm600m-int8', 'esmc_600m']:
        folder = old / name
        folder.mkdir(parents=True)
        (folder / f'{name}.bin.zst').write_bytes(b'existing reference')
    root = localfold.prepare_references(tmp_path)
    assert root == tmp_path / 'localfold'
    linked = root / 'weights/mol1-int8/mol1-int8.bin.zst'
    assert linked.is_symlink()
    assert linked.read_bytes() == b'existing reference'
    assert (root / 'weights/esmc_600m/esmc_600m.bin.zst').is_symlink()
    assert not (root / 'weights/esmfold2_lm600m-int8').exists()
    assert localfold.prepare_references(tmp_path) == root


def test_localfold_native_confidence_json_is_the_pae_artifact(tmp_path):
    structure = tmp_path / 'seed-42_sample-0' / 'probe_seed-42_sample-0_model.cif'
    structure.parent.mkdir(parents=True)
    sidecar = structure.with_name('probe_seed-42_sample-0_confidences.json')
    sidecar.write_text(json.dumps({'pae': [[0.0, 4.0], [4.0, 0.0]]}))
    assert prediction._localfold_pae_path(structure, 'mol1') == sidecar
    assert prediction._localfold_pae_path(structure, 'bioemu2') is None
    assert sidecar in cofolding_results.confidence_files(structure)
    matrix = cofolding_results._pae_matrix(cofolding_results.read_payload(sidecar))
    assert matrix.tolist() == [[0.0, 4.0], [4.0, 0.0]]


def test_localfold_chain_preferences_are_normalized_for_alignment_and_templates():
    value = prediction._prepare_prediction(
        engine='localfold', tasks=[{
            'protein_sequences': 'ACDEFG',
            'molecule_chain_options': [
                {'id': 'A', 'use_msa': False, 'use_template': False},
            ],
        }], samples=1, seed=42,
        engine_options={'model': 'alphafold3', 'msa_mode': 'mmseqs2_server'},
        accept_alphafold3_terms=True, resolve_random_seed=False)
    chain = value.normalized_tasks[0]['molecule_chain_options']['A']
    assert chain['use_msa'] is False
    assert chain['use_template'] is False


def test_editor_emits_localfold_per_chain_alignment_and_template_preferences():
    model = cofolding_input.default_model()
    model['cards'][0].update({
        'value': 'ACDEFG', 'localfold_alignment': 'single_sequence',
        'localfold_use_template': 'false',
    })
    model['settings']['msa_mode'] = 'automatic'
    payload = cofolding_input.model_to_payload(model)
    chain = payload['molecule_chain_options'][0]
    assert chain['use_msa'] is False
    assert chain['use_template'] is False
    assert chain['msa_a3m_text'] == ''


@pytest.mark.parametrize('model', ['af2_ptm', 'af2_multimer'])
def test_af2_rejects_nonprotein_inputs_instead_of_dropping_them(model):
    with pytest.raises(ValueError, match='does not support these molecule inputs'):
        prepare(model, dna_sequences='ACGT')


def test_bioemu2_rejects_multistep_fk_propagation_before_starting_a_job():
    with pytest.raises(ValueError, match='only one propagation step'):
        prepare('bioemu2', engine_options={'bioemu2_sampler': 'edm_fk', 'bioemu2_propagation_steps': 2})


def test_invalid_contacts_and_empty_potentials_are_preflight_errors():
    with pytest.raises(ValueError, match='JSON list'):
        prepare('bioemu2', engine_options={'bioemu2_contacts': '{}'})
    with pytest.raises(ValueError, match='must not be empty'):
        prepare('bioemu2', engine_options={'bioemu2_potentials': '# empty'})


def test_propagation_accepts_mmcif_comments_before_the_data_block():
    cif = '# LocalFold output notice\ndata_sample\nloop_\n_atom_site.Cartn_x\n0.0\n'
    value = prepare('bioemu2', engine_options={'bioemu2_structure_text': cif})
    assert value.selected_options['bioemu2_structure_text'] == cif


def test_esm_models_without_an_msa_encoder_require_the_language_model():
    with pytest.raises(ValueError, match='language model or an MSA'):
        prepare('esmfold2_lm600m', engine_options={'use_esm_embeddings': False})


def test_modified_residues_reach_the_native_af3_json(tmp_path):
    task = prepare(tasks=[{
        'protein_sequences': 'ACDEFG',
        'molecule_chain_options': [{'id': 'A', 'modifications_text': '3, MSE'}],
    }]).normalized_tasks[0]
    path, _ = prediction._stage_native_input_without_msa(
        tmp_path, profile=prediction.ENGINE_PROFILES['localfold'], item_id='probe', job_name='probe',
        proteins=task['proteins'], dna=[], rna=[], ligand_ccd_codes=[], ligand_smiles=[],
        molecule_chain_options=task['molecule_chain_options'], seed=42, samples=1,
        target_pdb_text='', target_chain_id='A', target_msa_a3m='', binder_sequence='', num_recycles=20)
    assert json.loads(path.read_text())['sequences'][0]['protein']['modifications'] == [
        {'ptmPosition': 3, 'ptmType': 'MSE'}]


def test_bioemu2_cannot_gain_a_confidence_score_from_a_zero_bfactor(tmp_path):
    row = cofolding_results.enrich_candidate(tmp_path, {
        'metrics': {'plddt_mean': 0., 'ptm': 0., 'physicality_passed': True},
        'raw_metadata': {'confidence_available': False, 'confidence_outputs': {'ranking_score': 0.}},
    })
    assert row['metrics'] == {'physicality_passed': True}


def test_localfold_is_available_in_the_api_and_pipeline_catalog():
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, OPTIMIZED_PREDICTION_NODE_ENGINES
    schema = workflow_schema('refolding.optimized_prediction')
    assert 'localfold' in str(schema)
    assert NODE_CATALOG['prediction_localfold']['label'] == 'LocalFold'
    assert OPTIMIZED_PREDICTION_NODE_ENGINES['prediction_localfold'] == 'localfold'


@pytest.mark.parametrize('options,samples,expected', [
    ({'model': 'mol1', 'seed_count': 2}, 3, 12),
    ({'model': 'af2_ptm', 'af2_num_models': 5, 'seed_count': 2}, 1, 20),
    ({'model': 'bioemu2', 'bioemu2_structure_text': 'cif', 'bioemu2_propagation_steps': 4}, 3, 24),
])
def test_native_execution_counts_seeds_parameter_sets_and_propagation_frames(options, samples, expected):
    from mn_molecular_design_studio.native_execution import request_contract
    contract = request_contract({'workflow': 'refolding.optimized_prediction', 'parameters': {
        'engine': 'localfold', 'engine_options': options, 'samples': samples}}, 2)
    assert contract['native_invocations'] == 1
    assert contract['prediction_outputs'] == expected


def test_entity_rows_reach_native_json_without_consuming_chain_ids(tmp_path):
    model = cofolding_input.default_model()
    model['settings'].update(model='mol1', msa_mode='automatic')
    model['cards'] = [
        {'kind': 'crosslink', 'name': 'DSSO', 'residue1': ['A', 2], 'residue2': ['B', 2]},
        {'kind': 'protein', 'value': 'AKSC', 'copies': 2, 'modifications': []},
        {'kind': 'modification', 'chain_id': 'B', 'position': 3, 'code': 'SEP'},
        {'kind': 'dna', 'value': 'ACGT', 'copies': 1, 'modifications': []},
        {'kind': 'modification', 'chain_id': 'C', 'position': 2, 'code': '5CM'},
        {'kind': 'bond', 'atom1': ['A', 4, 'SG'], 'atom2': ['B', 4, 'SG']},
    ]
    payload = cofolding_input.model_to_payload(model)
    assert payload['input_error'] == ''
    assert [row['id'] for row in payload['chain_catalog']] == ['A', 'B', 'C']
    value = prepare(tasks=[{key: val for key, val in payload.items() if key not in {'input_error', 'chain_catalog'}}])
    task = value.normalized_tasks[0]
    assert task['molecule_chain_options']['A']['modifications'] == [{'position': 3, 'ccd': 'SEP'}]
    assert task['molecule_chain_options']['B']['modifications'] == [{'position': 3, 'ccd': 'SEP'}]
    path, _ = prediction._stage_native_input_without_msa(
        tmp_path, profile=value.profile, item_id='probe', job_name='probe',
        proteins=task['proteins'], dna=task['dna'], rna=[], ligand_ccd_codes=[], ligand_smiles=[],
        molecule_chain_options=task['molecule_chain_options'], restraints=task['restraints'],
        seed=42, samples=1, target_pdb_text='', target_chain_id='A', target_msa_a3m='',
        binder_sequence='', num_recycles=3)
    native = json.loads(path.read_text())
    assert native['bondedAtomPairs'] == [[['A', 4, 'SG'], ['B', 4, 'SG']]]
    assert native['crosslinks'] == [{'name': 'DSSO', 'residue_pairs': [[['A', 2], ['B', 2]]]}]
    assert native['sequences'][2]['dna']['modifications'] == [{'basePosition': 2, 'modificationType': '5CM'}]


def test_incomplete_modification_is_a_form_error_not_an_exception():
    model = cofolding_input.default_model()
    model['settings']['msa_mode'] = 'automatic'
    model['cards'][0]['value'] = 'ACDE'
    model['cards'].append({'kind': 'modification', 'chain_id': 'A', 'position': '', 'code': 'MSE'})
    assert 'position' in cofolding_input.model_to_payload(model)['input_error']


@pytest.mark.parametrize('restraint,match', [
    ({'bond': {'atom1': ['A', 2, ''], 'atom2': ['A', 3, 'SG']}}, 'atom name'),
    ({'crosslink': {'name': 'DSSO', 'residue1': ['A', 999], 'residue2': ['A', 2]}}, 'outside'),
    ({'crosslink': {'name': 'unknown', 'residue1': ['A', 1], 'residue2': ['A', 2]}}, 'supported'),
    ({'crosslink': {'name': 'DSSO', 'residue1': ['A', 1], 'residue2': ['A', 1]}}, 'two crosslinks'),
])
def test_invalid_entity_references_are_rejected_before_execution(restraint, match):
    with pytest.raises(ValueError, match=match):
        prepare(tasks=[{'protein_sequences': 'AKSC', 'restraints': [restraint]}])


@pytest.mark.parametrize('model', ['af2_ptm', 'esmfold2_lm600m', 'esmfold2_lm300m'])
def test_models_without_covalent_inputs_refuse_retained_bonds(model):
    with pytest.raises(ValueError, match='does not support bond'):
        prepare(model, tasks=[{'protein_sequences': 'ACSC', 'restraints': [
            {'bond': {'atom1': ['A', 2, 'SG'], 'atom2': ['A', 4, 'SG']}}
        ]}])


def test_localfold_controls_are_unique_and_switchable_without_server():
    profiles = cofolding_input._localfold_editor_profiles()
    for name, profile in profiles.items():
        keys = [field['key'] for field in profile['prediction']['fields']]
        assert len(keys) == len(set(keys)), name
        assert keys.count('seed') == 1
    assert profiles['mol1']['bonds'] and profiles['mol1']['crosslinks']
    assert profiles['bioemu2']['bonds'] and profiles['bioemu2']['contacts']
    assert not profiles['bioemu2']['crosslinks'] and not profiles['bioemu2']['cyclic']
    assert not profiles['esmfold2_lm600m']['bonds']
    assert not profiles['esmfold2_lm600m']['crosslinks']
    assert profiles['esmfold2_fast']['bonds']
    assert profiles['bioemu2']['msa']
    assert profiles['esmfold2_lm600m']['prediction']['samples'] == 5


def test_bioemu2_all_atom_inputs_and_contact_rows_reach_native_json(tmp_path):
    model = cofolding_input.default_model()
    model['settings'].update(model='bioemu2', msa_mode='automatic')
    model['compound_name'] = 'Glycerol complex'
    model['cards'] = [
        {'kind': 'protein', 'value': 'AKSC', 'copies': 1},
        {'kind': 'dna', 'value': 'ACGT', 'copies': 1},
        {'kind': 'rna', 'value': 'ACGU', 'copies': 1},
        {'kind': 'ligand', 'mode': 'ccd', 'value': 'GOL', 'copies': 1},
        {'kind': 'modification', 'chain_id': 'A', 'position': 3, 'code': 'SEP'},
        {'kind': 'bond', 'atom1': ['A', 4, 'SG'], 'atom2': ['D', 1, 'C1']},
        {'kind': 'contact', 'token1': ['B', 2, ''], 'token2': ['D', 1, 'O1'], 'max_distance': 8},
    ]
    payload = cofolding_input.model_to_payload(model)
    assert payload['input_error'] == ''
    value = prepare('bioemu2', tasks=[{k: v for k, v in payload.items() if k not in {'input_error', 'chain_catalog'}}])
    task = value.normalized_tasks[0]
    path, _ = prediction._stage_native_input_without_msa(
        tmp_path, profile=value.profile, item_id='probe', job_name='probe',
        proteins=task['proteins'], dna=task['dna'], rna=task['rna'],
        ligand_ccd_codes=task['ligand_ccd_codes'], ligand_smiles=task['ligand_smiles'],
        molecule_chain_options=task['molecule_chain_options'], restraints=task['restraints'],
        seed=42, samples=1, target_pdb_text='', target_chain_id='A', target_msa_a3m='',
        binder_sequence='', num_recycles=3)
    native = json.loads(path.read_text())
    assert [next(iter(entity)) for entity in native['sequences']] == ['protein', 'dna', 'rna', 'ligand']
    assert native['sequences'][0]['protein']['modifications'] == [{'ptmPosition': 3, 'ptmType': 'SEP'}]
    assert native['bondedAtomPairs'] == [[['A', 4, 'SG'], ['D', 1, 'C1']]]
    assert native['constraints'] == [{'contact': {'token1': ['B', 2], 'token2': ['D', 1, 'O1'], 'maxDistance': 8}}]


def test_localfold_pocket_endpoints_keep_native_ligand_atom_format():
    value = prepare('boltz2', tasks=[{'protein_sequences': 'AKSC', 'ligand_ccd_codes': 'GOL', 'compound_name': 'Glycerol', 'restraints': [
        {'pocket': {'binder': 'A', 'contacts': [['B', 1, 'O1']], 'max_distance': 6}},
    ]}])
    pocket = value.normalized_tasks[0]['restraints'][0]['pocket']
    assert pocket['contacts'] == [['B', 1, 'O1']]


@pytest.mark.parametrize('model,restraint,match', [
    ('esmfold2_lm600m', {'contact': {'token1': ['A', 1], 'token2': ['A', 2], 'max_distance': 8}}, 'contact conditioning'),
    ('bioemu2', {'pocket': {'binder': 'A', 'contacts': [['A', 1]], 'max_distance': 8}}, 'pocket conditioning'),
    ('chai1', {'pocket': {'binder': 'B', 'contacts': [['A', 1]], 'max_distance': 8}}, 'SMILES ligand'),
])
def test_conditioning_is_rejected_for_unsupported_model_inputs(model, restraint, match):
    with pytest.raises(ValueError, match=match):
        prepare(model, tasks=[{'protein_sequences': 'AKSC', 'ligand_ccd_codes': 'GOL', 'compound_name': 'Glycerol', 'restraints': [restraint]}])


def test_schema_discovery_imports_only_the_requested_workflow(monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation
    original = automation._get_function
    seen = []
    def get(spec):
        seen.append(spec.workflow_id)
        return original(spec)
    monkeypatch.setattr(automation, '_get_function', get)
    rows = automation.workflow_catalog(workflow_id='refolding.optimized_prediction')
    assert [row['workflow'] for row in rows] == ['refolding.optimized_prediction']
    assert set(seen) == {'refolding.optimized_prediction'}


def test_modified_residue_identity_is_retained_in_results(tmp_path):
    import gemmi
    from mn_molecular_design_studio.scientific_contracts import verify_prediction_structure
    structure = gemmi.Structure()
    model = gemmi.Model('1')
    chain = gemmi.Chain('A')
    for number, name in enumerate(['ALA', 'LYS', 'SEP', 'CYS'], 1):
        residue = gemmi.Residue()
        residue.name = name
        residue.seqid = gemmi.SeqId(number, ' ')
        atom = gemmi.Atom()
        atom.name = 'CA'
        atom.element = gemmi.Element('C')
        atom.pos = gemmi.Position(number * 3.8, 0, 0)
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    path = tmp_path / 'modified.cif'
    structure.make_mmcif_document().write_file(str(path))
    task = {'proteins': [{'id': 'A', 'sequence': 'AKSC'}], 'molecule_chain_options': {
        'A': {'modifications': [{'position': 3, 'ccd': 'SEP'}]}}}
    assert verify_prediction_structure(path, task)['sequence_chain_identity'] == 'pass'
    task['molecule_chain_options']['A']['modifications'][0]['ccd'] = 'TPO'
    with pytest.raises(ValueError, match='modified residue A3'):
        verify_prediction_structure(path, task)
