from pathlib import Path
import json

import numpy as np
import pytest

from mn_molecular_design_studio import cofolding_results as results
from mn_molecular_design_studio.cofolding_input import prediction_settings_for
from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES, _confidence_sidecar


@pytest.mark.parametrize('engine,structure_name,data_name,summary_name,plddt_key,pae_key', [
    ('af3_torch', 'job_seed-1_sample-0_model.cif', 'job_seed-1_sample-0_confidences.json', 'job_seed-1_sample-0_summary_confidences.json', 'atom_plddts', 'pae'),
    ('openfold3_ob0', 'job_seed_1_sample_0_model.cif', 'job_seed_1_sample_0_confidences.json', 'job_seed_1_sample_0_confidences_aggregated.json', 'plddt', 'pae'),
    ('protenix_v2', 'job_sample_0.cif', 'job_full_data_sample_0.json', 'job_summary_confidence_sample_0.json', 'atom_plddt', 'token_pair_pae'),
    ('opendde_native', 'job_sample_0.cif', 'job_full_data_sample_0.json', 'job_summary_confidence_sample_0.json', 'atom_plddt', 'token_pair_pae'),
    ('rf3_foundry', 'job_model.cif', 'job_confidences.json', 'job_summary_confidences.json', 'atom_plddts', 'pae'),
    ('colabfold', 'job_unrelaxed_rank_001_model_1_seed_001.pdb', 'job_scores_rank_001_model_1_seed_001.json', 'job_scores_rank_001_model_1_seed_001.json', 'plddt', 'pae'),
])
def test_native_per_model_confidence_pairs_do_not_borrow_neighbour_scores(tmp_path, engine, structure_name, data_name, summary_name, plddt_key, pae_key):
    structure = tmp_path / structure_name
    structure.touch()
    payload = {plddt_key: [91., 83.], pae_key: [[0., 3.], [5., 0.]]}
    (tmp_path / data_name).write_text(json.dumps(payload))
    summary_path = tmp_path / summary_name
    summary = payload if data_name == summary_name else {}
    summary_path.write_text(json.dumps({**summary, 'ptm': .82, 'iptm': .76, 'ranking_score': .79}))
    # An unrelated model has a different score and must never be picked by a broad glob.
    (tmp_path / 'other_summary_confidences.json').write_text(json.dumps({'iptm': .01}))
    data = results.read_confidence(structure, tuple(results.confidence_files(structure)))
    np.testing.assert_allclose(data['plddt'], [91, 83])
    np.testing.assert_allclose(data['pae'], [[0, 3], [5, 0]])
    assert data['metrics']['iptm'] == .76
    assert _confidence_sidecar(structure)['ranking_score'] == .79


@pytest.mark.parametrize('structure_name,npz_name', [
    ('job_model_0.pdb', 'pae_job_model_0.npz'),
    ('job__full__s1_x0.cif', 'job__full__s1_x0_pae.npz'),
    ('pred.model_idx_0.cif', 'scores.model_idx_0.npz'),
    ('job_seed-1_sample-0_model.cif', 'job_seed-1_sample-0_confidence.npz'),
    ('job.pdb', 'job.npz'),
])
def test_native_npz_confidence_arrays_and_probability_scale(tmp_path, structure_name, npz_name):
    structure = tmp_path / structure_name
    structure.touch()
    np.savez(tmp_path / npz_name, plddt=[.91, .83], pae=[[0, 3], [5, 0]], ptm=.82)
    data = results.read_confidence(structure, tuple(results.confidence_files(structure)))
    np.testing.assert_allclose(data['plddt'], [91, 83])
    assert data['pae'].shape == (2, 2)
    assert data['metrics']['ptm'] == .82


def test_rankings_are_per_complex_and_missing_primary_scores_sort_last():
    def candidate(item, score=None, plddt=90):
        return {'raw_metadata': {'batch_item_id': item}, 'metrics': {'ranking_score': score, 'plddt_mean': plddt}}
    rows = [candidate('A', .7), candidate('B', .1), candidate('A', .9), candidate('A', None, 99), candidate('B', .5)]
    ranked = results.ranked_candidates(rows)
    assert [row['metrics']['ranking_score'] for row in ranked] == [.9, .7, None, .5, .1]


def test_batch_summary_keeps_compounds_separate_and_picks_each_best_model():
    candidates = [
        {'candidate_id': 'bepristat_a_low', 'raw_metadata': {'batch_item_id': 'task-a', 'source_compound_id': 'Bepristat_1a'},
         'metrics': {'ranking_score': .61}},
        {'candidate_id': 'rutin_best', 'raw_metadata': {'batch_item_id': 'task-r', 'source_compound_id': 'Rutin'},
         'metrics': {'ranking_score': .92}},
        {'candidate_id': 'bepristat_a_best', 'raw_metadata': {'batch_item_id': 'task-a', 'source_compound_id': 'Bepristat_1a'},
         'metrics': {'ranking_score': .87}},
        {'candidate_id': 'rutin_low', 'raw_metadata': {'batch_item_id': 'task-r', 'source_compound_id': 'Rutin'},
         'metrics': {'ranking_score': .45}},
    ]

    best = results.best_candidates_by_complex(candidates)

    assert [results.complex_label(item) for item in best] == ['Bepristat_1a', 'Rutin']
    assert [item['candidate_id'] for item in best] == ['bepristat_a_best', 'rutin_best']


def test_complex_label_uses_task_name_when_no_compound_id_exists():
    assert results.complex_label({'raw_metadata': {'input_name': 'sequence_7'}}, 'Complex 7') == 'sequence_7'


def test_compound_identity_keeps_id_name_and_smiles_separate():
    candidate = {"raw_metadata": {
        "compound_id": "CID-5291", "compound_name": "Imatinib", "ligand_smiles": ["CCN", "CO"],
    }}
    assert results.compound_identity(candidate) == {
        "compound_id": "CID-5291", "compound_name": "Imatinib", "smiles": "CCN; CO",
    }
    assert results.complex_label(candidate) == "CID-5291 · Imatinib"


def test_summary_groups_old_candidates_by_input_name_without_batch_ids():
    candidates = [
        {'candidate_id': 'seq-a-model-2', 'raw_metadata': {'input_name': 'sequence_A'}, 'metrics': {'ranking_score': .5}},
        {'candidate_id': 'seq-b-model-1', 'raw_metadata': {'input_name': 'sequence_B'}, 'metrics': {'ranking_score': .7}},
        {'candidate_id': 'seq-a-model-1', 'raw_metadata': {'input_name': 'sequence_A'}, 'metrics': {'ranking_score': .8}},
    ]

    best = results.best_candidates_by_complex(candidates)

    assert [item['candidate_id'] for item in best] == ['seq-a-model-1', 'seq-b-model-1']


def test_prediction_input_labels_fill_names_missing_from_older_candidate_metadata(tmp_path):
    path = tmp_path / 'artifacts/raw/optimized_prediction/input_mapping.json'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({'items': [
        {'id': 'task-a', 'name': 'Bepristat_1a'},
        {'id': 'task-b', 'name': 'Rutin'},
    ]}))
    labels = results.prediction_input_labels(tmp_path)

    candidate = {'raw_metadata': {'batch_item_id': 'task-b'}}
    assert labels == {'task-a': 'Bepristat_1a', 'task-b': 'Rutin'}
    assert results.complex_label(candidate, 'Complex 2', labels) == 'Rutin'


def test_model_sidecars_are_exact_and_job_paths_cannot_escape(tmp_path):
    model = tmp_path / 'job_sample_1.cif'
    model.touch()
    (tmp_path / 'job_summary_confidence_sample_0.json').write_text('{"iptm":0.99}')
    assert results.confidence_files(model) == []
    assert results.local_artifact(tmp_path, '../outside.json') is None
    assert results.local_artifact(tmp_path, '/work/job_sample_1.cif') == model


def test_msa_preview_is_one_per_chain_and_prefers_staged_reference_copy(tmp_path):
    staged = tmp_path / 'artifacts/raw/af3_jax_cofolding/msas/chain_A.a3m'
    generated = tmp_path / 'artifacts/raw/af3_jax_cofolding/output/uniref.a3m'
    for path in (generated, staged):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('>query\nACD\n>hit\nA-D\n')
    candidate = {'raw_metadata': {'protein_chains': {'A': 'ACD'}}}
    assert results.alignment_files(tmp_path, candidate) == [('A · A3M alignment', staged)]


@pytest.mark.parametrize('engine', ENGINE_PROFILES)
def test_shared_input_controls_keep_each_engine_supported_options(engine):
    profile = ENGINE_PROFILES[engine]
    config = prediction_settings_for(profile)
    fields = {field['key']: field for field in config['fields']}
    assert ('seed' in fields) == profile.get('supports_seed', True)
    if 'seed' in fields:
        assert fields['seed']['default'] == 0
    assert fields['num_recycles']['default'] == profile['default_recycles']
    assert ('optimization_mode' in fields) == (len(profile['modes']) > 1)
    if 'optimization_mode' in fields:
        assert [mode for mode, _ in fields['optimization_mode']['choices']] == list(profile['modes'])
    assert set(field['key'] for field in profile.get('settings', ())) <= fields.keys()
    assert ('accept_terms' in fields) == bool(profile.get('terms'))


def test_boltzgen_shared_archive_uses_the_selected_sample_for_arrays_and_metrics(tmp_path):
    root = tmp_path / 'artifacts/raw/optimized_prediction'
    structure = root / 'output/item/item_sample2.cif'
    structure.parent.mkdir(parents=True)
    structure.touch()
    archive = root / 'boltzgen_data/fold_out_npz/item.npz'
    archive.parent.mkdir(parents=True)
    np.savez(archive, plddt=[[.91, .83], [.61, .53]], pae=[[[0, 3], [5, 0]], [[0, 13], [15, 0]]], iptm=[.86, .36])
    assert archive in results.confidence_files(structure)
    data = results.read_confidence(structure, tuple(results.confidence_files(structure)))
    np.testing.assert_allclose(data['plddt'], [61, 53])
    np.testing.assert_allclose(data['pae'], [[0, 13], [15, 0]])
    assert data['metrics']['iptm'] == .36
    assert _confidence_sidecar(structure)['iptm'] == .36


def test_esmfold2_rows_match_exact_structure_filename(tmp_path):
    structure = tmp_path / 'cif_all/job__full__s1_x1.cif'
    structure.parent.mkdir()
    structure.touch()
    (tmp_path / 'pred_rows.jsonl').write_text('\n'.join(json.dumps(row) for row in [
        {'pred_file': 'job__full__s1_x0.cif', 'iptm': .9},
        {'pred_file': structure.name, 'iptm': .7, 'ptm': .8},
    ]))
    assert _confidence_sidecar(structure)['iptm'] == .7


@pytest.mark.parametrize('structure_name', ['rank1_item.cif', 'item.cif'])
def test_boltzgen_ranked_complex_uses_native_best_refolding_sample(tmp_path, structure_name):
    structure = tmp_path / 'final_ranked_designs/final_1_designs' / structure_name
    structure.parent.mkdir(parents=True)
    structure.touch()
    archive = tmp_path / 'intermediate_designs_inverse_folded/fold_out_npz/item.npz'
    archive.parent.mkdir(parents=True)
    # Highest iPTM alone is sample 0; the native weighted score selects sample 1.
    np.savez(archive, plddt=[[.91, .83], [.61, .53]],
             pae=[[[0, 3], [5, 0]], [[0, 13], [15, 0]]],
             design_to_target_iptm=[.81, .8], design_ptm=[.1, .9], iptm=[.86, .36])
    data = results.read_confidence(structure, (archive,))
    np.testing.assert_allclose(data['pae'], [[0, 13], [15, 0]])
    np.testing.assert_allclose(data['plddt'], [61, 53])
    assert data['metrics']['iptm'] == .36
    # An unassociated structure must not borrow this archive's best sample.
    assert results.sample_payload(structure.with_name('other.cif'), archive,
                                  results.read_payload(archive))['pae'].shape == (2, 2, 2)


def test_boltzgen_native_numbered_refold_sample_is_zero_based(tmp_path):
    structure = tmp_path / 'refold_cif/item_sample0.cif'
    archive = tmp_path / 'fold_out_npz/item.npz'
    payload = {'pae': np.asarray([[[0, 3], [5, 0]], [[0, 13], [15, 0]]])}
    np.testing.assert_allclose(results.sample_payload(structure, archive, payload)['pae'],
                               [[0, 3], [5, 0]])
