from pathlib import Path
import json

import numpy as np
import pytest

from mn_molecular_design_studio import cofolding_results as results
from mn_molecular_design_studio.cofolding_input import prediction_settings_for
from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES, _confidence_sidecar


def test_prediction_viewer_loads_coordinates_without_native_occupancy(tmp_path):
    import gemmi

    structure = gemmi.read_pdb_string(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 82.00           C\nEND\n"
    )
    document = structure.make_mmcif_document()
    block = document.sole_block()
    atoms = block.get_mmcif_category('_atom_site.', raw=True)
    atoms.pop('occupancy')
    block.set_mmcif_category('_atom_site.', atoms, raw=True)
    source = tmp_path / 'native.cif'
    document.write_file(str(source))
    original = source.read_bytes()
    assert not len(gemmi.read_structure(str(source)))
    scratch = tmp_path / 'preview'
    scratch.mkdir()

    preview = results.viewer_structure(source, scratch)
    atom = gemmi.read_structure(str(preview))[0][0][0][0]
    assert (atom.pos.x, atom.pos.y, atom.pos.z) == (1.0, 2.0, 3.0)
    assert atom.b_iso == 82.0
    assert source.read_bytes() == original
    assert results.viewer_structure(preview, scratch) == preview


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


@pytest.mark.parametrize("job_type", ["complex_refolding", "monomer_refolding", "rfd4_structure_prediction"])
def test_legacy_standalone_prediction_routes_without_replacing_design_results(job_type):
    candidate_job = {"tool": "rfd4_proteina", "job_type": job_type}
    assert results.is_prediction_run(candidate_job)
    assert not results.is_prediction_run(candidate_job, {"workflow": "design.rfd4_proteina"})
    assert not results.is_prediction_run({"tool": "rfd4_proteina", "job_type": "rfd4_proteina_design"})


def test_model_picker_separates_inputs_and_disambiguates_repeated_model_indices():
    def candidate(candidate_id, reference, seed, rank):
        return {"candidate_id": candidate_id, "complex_pdb": "prediction_model_0.cif",
                "raw_metadata": {"validation_reference": {"candidate_id": reference}, "seed": seed},
                "metrics": {"native_rank": rank}}

    rows = results.ranked_candidates([
        candidate("a-second", "input-a", 0, 2), candidate("b", "input-b", 0, 1),
        candidate("a-first", "input-a", 0, 1),
    ])
    assert [row["candidate_id"] for row in rows] == ["a-first", "a-second", "b"]
    assert results.model_labels(rows[:2]) == ["Model 0 · seed 0 · entry 1 (best)", "Model 0 · seed 0 · entry 2"]
    assert results.model_scores({"raw_metadata": {"confidence_available": False}, "metrics": {"ptm": .9}}) == {}


def test_peer_confidence_manifest_pairs_only_files_for_the_selected_model(tmp_path):
    for index in range(2):
        (tmp_path / f"job_sample_{index}.cif").touch()
        (tmp_path / f"job_full_data_sample_{index}.json").write_text("{}")
    rows = [{"candidate_id": f"model-{index}", "complex_pdb": f"job_sample_{index}.cif"} for index in range(2)]
    assert results.confidence_artifacts(tmp_path, rows) == {
        "model-0": [{"path": "job_full_data_sample_0.json", "size_bytes": 2}],
        "model-1": [{"path": "job_full_data_sample_1.json", "size_bytes": 2}],
    }


def test_external_refold_shows_its_own_model_confidence_and_never_inherits_input_scores(tmp_path):
    native = tmp_path / "job_seed_6_sample_1_model.cif"
    native.touch()
    (tmp_path / "role_normalized.pdb").touch()
    wrong_sidecar = tmp_path / "job_seed_6_sample_0_confidences.json"
    wrong_sidecar.write_text('{"pae":[[0,30],[30,0]]}')
    correct_sidecar = tmp_path / "job_seed_6_sample_1_confidences.json"
    correct_sidecar.write_text('{"pae":[[0,2],[3,0]]}')
    (tmp_path / "job_seed_6_sample_1_confidences_aggregated.json").write_text('{"iptm":0.41,"avg_plddt":89}')
    row = {"candidate_id": "refold", "complex_pdb": "role_normalized.pdb",
           "metrics": {"complex_refolding_backend": "openfold3", "iptm": .95, "ptm": .99, "openfold3_iptm": .1},
           "raw_metadata": {"source_candidate": {"candidate_id": "parent", "metrics": {"iptm": .95, "ptm": .99}},
                            "seed": 1, "confidence_outputs": {"iptm": .95},
                            "original_prediction": native.name, "role_normalized_prediction": "role_normalized.pdb",
                            "pae_source_path": wrong_sidecar.name, "pae_path": wrong_sidecar.name}}

    presented = results.presentation_candidate(tmp_path, row)
    assert results.model_seed(presented) == 6
    assert results.model_index(presented, 0) == 1
    assert results.complex_group_key(presented) == "parent"
    assert results.model_scores(presented) == {"Mean pLDDT (%)": 89, "Ranking score": None, "pTM": None, "ipTM": .41}
    structure, sidecars = results.prediction_files(tmp_path, presented)
    assert wrong_sidecar not in sidecars
    assert "pae_path" not in presented["raw_metadata"]
    assert results.read_confidence(structure, tuple(sidecars))["pae"].tolist() == [[0, 2], [3, 0]]
    assert row["metrics"]["iptm"] == .95


def test_connected_prediction_structure_resolves_from_managed_source_run(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import portable_paths

    runs = tmp_path / 'runs'
    campaign = runs / 'campaign' / 'campaign-1'
    prediction = runs / 'refolding-validation' / 'prediction-1' / 'artifacts' / 'model.pdb'
    campaign.mkdir(parents=True)
    prediction.parent.mkdir(parents=True)
    prediction.write_text('ATOM\n')
    monkeypatch.setattr(portable_paths, 'runs_root', lambda: runs)

    assert results.local_artifact(campaign, 'runs:///refolding-validation/prediction-1/artifacts/model.pdb') == prediction
    assert results.local_artifact(campaign, str(prediction)) == prediction
    assert results.local_artifact(campaign, str(tmp_path / 'outside.pdb')) is None


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
