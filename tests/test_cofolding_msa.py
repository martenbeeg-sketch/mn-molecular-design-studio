"""MSA policy contracts, with no network searches or prediction processes."""
from pathlib import Path
import hashlib
import json

import pytest
import yaml

from mn_molecular_design_studio import cofolding_msa as msa, cofolding_results as results
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as workflow

SEQUENCE = 'ACDEFGHIK'
ALIGNMENT = f'>query\n{SEQUENCE}\n>homolog\nACDE-GHIK\n>insertion\nACDEaaFGHIK\n'


@pytest.fixture
def reference(tmp_path, monkeypatch):
    root = tmp_path / 'reference_files'
    monkeypatch.setattr(msa, 'reference_root', lambda: root)
    monkeypatch.setattr(workflow, 'reference_root', lambda: root)
    return root


def cache(root, sequence=SEQUENCE, text=ALIGNMENT):
    path = root / 'msa_cache' / f'{hashlib.sha256(sequence.encode()).hexdigest()}.a3m'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def forbid_fetch(*args, **kwargs):
    pytest.fail('A cache hit, disabled MSA, or unsupported engine must not contact a server')


def test_cache_hit_stages_capped_input_and_keeps_full_reusable_alignment(tmp_path, reference, monkeypatch):
    saved = cache(reference)
    monkeypatch.setattr(msa, 'fetch_colabfold', forbid_fetch)
    staged = msa.prepare_chain(tmp_path / 'job', SEQUENCE, 'A', depth=2)
    assert len(msa.records(staged.read_text())) == 2
    assert saved.read_text() == ALIGNMENT
    manifest = json.loads((tmp_path / 'job/artifacts/raw/cofolding_msa/manifest.json').read_text())
    assert manifest['alignments'][0]['source'] == 'reference_cache'
    assert manifest['alignments'][0]['cached_rows'] == 3
    assert results.manifest_alignment_files(tmp_path / 'job', {})[0][1] == staged


def test_cache_miss_fetches_once_saves_both_and_reuses_across_engines(tmp_path, reference, monkeypatch):
    calls = []
    monkeypatch.setattr(msa, 'fetch_colabfold', lambda *args: calls.append(args) or ALIGNMENT)
    first = msa.prepare_chain(tmp_path / 'first-engine', SEQUENCE, 'A', depth=2)
    second = msa.prepare_chain(tmp_path / 'second-engine', SEQUENCE, 'B', depth=3)
    assert len(calls) == 1
    assert len(msa.records(first.read_text())) == 2
    assert second.read_text() == ALIGNMENT
    assert cache_path(reference).read_text() == ALIGNMENT


def cache_path(root):
    return root / 'msa_cache' / f'{hashlib.sha256(SEQUENCE.encode()).hexdigest()}.a3m'


@pytest.mark.parametrize('text', ['>wrong\nVVVVVVVVV\n', f'>query\n{SEQUENCE}\n>bad\nACD\n', ALIGNMENT + '\x00'])
def test_invalid_server_alignment_is_not_published(tmp_path, reference, monkeypatch, text):
    monkeypatch.setattr(msa, 'fetch_colabfold', lambda *args: text)
    with pytest.raises(ValueError):
        msa.prepare_chain(tmp_path / 'job', SEQUENCE, 'A')
    assert not cache_path(reference).exists()


def test_custom_alignment_cached_without_overwriting_deeper_reference(tmp_path, reference):
    cache(reference)
    custom = f'>query\n{SEQUENCE}\n'
    staged = msa.prepare_chain(tmp_path / 'job', SEQUENCE, 'A', supplied=custom)
    assert staged.read_text() == custom
    assert cache_path(reference).read_text() == ALIGNMENT


def test_legacy_reference_repository_is_reused_and_copied_to_canonical_cache(tmp_path, reference, monkeypatch):
    legacy = reference / 'boltz_models/msa_repository/legacy.a3m'
    legacy.parent.mkdir(parents=True)
    legacy.write_text(ALIGNMENT)
    monkeypatch.setattr(msa, 'fetch_colabfold', forbid_fetch)
    assert msa.prepare_chain(tmp_path / 'job', SEQUENCE, 'A').read_text() == ALIGNMENT
    assert cache_path(reference).read_text() == ALIGNMENT


def native_args(engine):
    return dict(profile=workflow.ENGINE_PROFILES[engine], item_id='complex1', job_name='example',
                proteins=[{'id': 'A', 'sequence': SEQUENCE}], dna=[], rna=[], ligand_ccd_codes=[], ligand_smiles=[],
                seed=1, samples=5, target_pdb_text='', target_chain_id='A', target_msa_a3m='',
                binder_sequence='', num_recycles=3)


@pytest.mark.parametrize('engine', [key for key in workflow.ENGINE_PROFILES if key != 'boltzgen_fold'])
def test_every_native_adapter_consumes_job_local_alignment_when_capable(tmp_path, reference, monkeypatch, engine):
    cache(reference)
    monkeypatch.setattr(msa, 'fetch_colabfold', forbid_fetch)
    run = tmp_path / 'job'
    profile = workflow.ENGINE_PROFILES[engine]
    path, details = workflow._stage_native_input(run, **native_args(engine), prepare_msas=True, msa_depth=2)
    if not profile['supports_msa']:
        assert details['msa_source'] == 'single_sequence'
        assert not (run / 'artifacts/raw/cofolding_msa/manifest.json').exists()
        return
    staged = run / details['msa_files_by_chain']['A']
    assert len(msa.records(staged.read_text())) == 2
    assert results.manifest_alignment_files(run, {'raw_metadata': {'batch_item_id': 'complex1'}})
    assert not results.manifest_alignment_files(run, {'raw_metadata': {'batch_item_id': 'other-complex'}})
    native_path = '/work/' + staged.relative_to(run).as_posix()
    kind = profile['input_kind']
    payload = yaml.safe_load(path.read_text()) if kind == 'boltz_yaml' else json.loads(path.read_text()) if path.suffix == '.json' else None
    if kind == 'af3_json':
        assert payload['sequences'][0]['protein']['unpairedMsa'] == staged.read_text()
    elif kind == 'esmfold2_json':
        assert payload['sequences'][0]['msa'] == native_path
        assert workflow._batch_options(profile, seed=1, num_recycles=3)['variant'] == 'full_msa'
    elif kind == 'chai1_json':
        assert payload['msa_dir'] == str(Path(native_path).parent)
        assert workflow._batch_options(profile, seed=1, num_recycles=3)['msa_depth'] == 1024
    elif kind in {'protenix_json', 'opendde_json'}:
        chain = payload[0]['sequences'][0]['proteinChain']
        assert chain['unpairedMsaPath'] == native_path
        assert (run / chain['pairedMsaPath'].removeprefix('/work/')).is_file()
        assert workflow._batch_options(profile, seed=1, num_recycles=3)['use_msa'] is True
    elif kind == 'mol1_json':
        chain = payload['sequences'][0]['proteinChain']
        assert chain['unpairedMsaPath'] == native_path
        assert (run / chain['pairedMsaPath'].removeprefix('/work/')).is_file()
        assert workflow._batch_options(profile, seed=1, num_recycles=3)['use_msa'] is True
    elif kind == 'rf3_json':
        assert payload['msa_paths']['A'] == native_path
    elif kind == 'openfold3_ob0_json':
        query = payload['queries']['example']
        assert query['use_msas'] and query['use_main_msas'] and not query['use_paired_msas']
        assert query['chains'][0]['main_msa_file_paths'] == native_path
        assert workflow._batch_options(profile, seed=1, num_recycles=3)['use_msa_server'] is False
    elif kind == 'boltz_yaml':
        assert payload['sequences'][0]['protein']['msa'] == native_path
    elif kind == 'multimer_fasta':
        assert path.suffix == '.a3m'
        assert path.read_text().startswith('#9\t1\n')
        assert workflow._batch_options(profile, seed=1, num_recycles=3)['msa_mode'] != 'single_sequence'


@pytest.mark.parametrize('engine', [key for key in workflow.ENGINE_PROFILES if key != 'boltzgen_fold'])
def test_disabled_msa_never_fetches(tmp_path, reference, monkeypatch, engine):
    monkeypatch.setattr(msa, 'fetch_colabfold', forbid_fetch)
    path, details = workflow._stage_native_input(tmp_path / 'job', **native_args(engine), prepare_msas=True, use_msa=False)
    assert details['msa_source'] == 'single_sequence'
    assert not (tmp_path / 'job/artifacts/raw/cofolding_msa/manifest.json').exists()


def test_colabfold_homomer_and_heteromer_preserve_unique_lengths_and_cardinality(tmp_path):
    one, two = tmp_path / 'one.a3m', tmp_path / 'two.a3m'
    one.write_text(ALIGNMENT)
    two.write_text('>query\nVVV\n>homolog\nV-V\n')
    chains = [{'id': 'A', 'sequence': SEQUENCE}, {'id': 'B', 'sequence': SEQUENCE}, {'id': 'C', 'sequence': 'VVV'}]
    text = msa.complex_a3m(chains, {'A': one, 'B': one, 'C': two})
    assert text.startswith(f'#9,3\t2,1\n>query\n{SEQUENCE}VVV\n')
    assert all(len(__import__('re').sub('[a-z.]', '', row)) == 12 for _, row in msa.records(text))


def test_mol1_stages_native_input_controls_and_capability_graph_node(tmp_path):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, OPTIMIZED_PREDICTION_NODE_ENGINES

    profile = workflow.ENGINE_PROFILES['mol1']
    path, _ = workflow._stage_native_input(
        tmp_path / 'job', **native_args('mol1'), prepare_msas=False,
        molecule_chain_options={'A': {'modifications': [{'ptmType': 'CCD_HY3', 'ptmPosition': 2}]}},
    )
    payload = json.loads(path.read_text())
    protein = payload['sequences'][0]['proteinChain']
    assert protein['id'] == ['A']
    assert protein['modifications'] == [{'ptmType': 'CCD_HY3', 'ptmPosition': 2}]
    options = workflow._batch_options(profile, seed=7, num_recycles=12,
                                      engine_options={'use_msa': True, 'step': 64, 'use_template': True,
                                                      'need_atom_confidence': True})
    assert options == {
        'model_name': 'mol-1-base', 'checkpoint': '/reference/mol1/checkpoint/mol-1-base.pt',
        'cycle': 12, 'step': 64, 'dtype': 'bf16', 'use_msa': True, 'use_template': True,
        'need_atom_confidence': True,
    }
    kind = next(kind for kind, engine in OPTIMIZED_PREDICTION_NODE_ENGINES.items() if engine == 'mol1')
    assert kind == 'prediction_mol1'
    assert NODE_CATALOG[kind]['workflow'] == 'refolding.optimized_prediction'
    assert NODE_CATALOG[kind]['formats'] == ['pdb', 'cif']
    assert profile['supports_templates'] is True


def test_mol1_native_confidence_and_pae_are_preserved_as_measured(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _pae_matrix_from_json

    structure = tmp_path / 'complex_sample_0.pdb'
    structure.write_text(
        'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 81.00           C\n'
        'ATOM      2  CA  CYS B   1       5.000   0.000   0.000  1.00 72.00           C\n'
        'END\n'
    )
    full_data = structure.with_name('complex_full_data_sample_0.json')
    full_data.write_text(json.dumps({
        'token_pair_pae': [[0.1, 4.0], [6.0, 0.2]], 'token_asym_id': [0, 1],
    }))
    task = {'proteins': [{'id': 'A'}, {'id': 'B'}], 'molecule_chain_options': {}}
    pae_path = workflow._mol1_pae_path(structure, task)
    assert pae_path == tmp_path / 'complex_sample_0_pae.json'
    payload = json.loads(pae_path.read_text())
    assert payload['source'] == 'Mol-1 token_pair_pae'
    assert _pae_matrix_from_json(pae_path) == [[0.1, 4.0], [6.0, 0.2]]
    assert workflow._mean_plddt(structure) == pytest.approx(76.5)
    assert workflow._mean_plddt(structure, ['B']) == pytest.approx(72.0)
    assert workflow._mol1_pae_path(structure, {'proteins': [{'id': 'B'}, {'id': 'A'}]}) is None


def test_historical_single_sequence_job_explains_why_there_is_no_alignment(tmp_path):
    path = tmp_path / 'artifacts/raw/af3_jax_cofolding/input_mapping.json'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({'msa_source': 'single_sequence'}))
    assert 'without an MSA' in results.alignment_status(tmp_path, {})
    assert results.manifest_alignment_files(tmp_path, {}) == []


def test_prefix_reuse_crops_alignment_columns_without_losing_insertions(tmp_path, reference, monkeypatch):
    cache(reference)
    monkeypatch.setattr(msa, 'fetch_colabfold', forbid_fetch)
    path = msa.prepare_chain(tmp_path / 'job', 'ACDEFG', 'A')
    rows = msa.records(path.read_text())
    assert rows[0][1] == 'ACDEFG'
    assert rows[2][1] == 'ACDEaaFG'
    assert all(len(__import__('re').sub('[a-z.]', '', row)) == 6 for _, row in rows)


def test_boltzgen_target_and_binder_use_shared_cache_without_required_upload(tmp_path, reference, monkeypatch):
    from Bio.SeqUtils import seq3
    cache(reference)
    binder_sequence = 'VVVVVVVVV'
    cache(reference, binder_sequence, f'>query\n{binder_sequence}\n>homolog\nVVVV-VVVV\n')
    monkeypatch.setattr(msa, 'fetch_colabfold', forbid_fetch)
    pdb = '\n'.join(f'ATOM  {index:5d}  CA  {seq3(letter).upper()} A{index:4d}    {index*3.8:8.3f}{0.:8.3f}{0.:8.3f}  1.00 90.00           C' for index, letter in enumerate(SEQUENCE, 1)) + '\nEND\n'
    args = native_args('boltzgen_fold')
    args.update(target_pdb_text=pdb, binder_sequence=binder_sequence)
    path, details = workflow._stage_native_input(tmp_path / 'job', **args, prepare_msas=True, msa_depth=2)
    target = Path(details['target_msa'])
    assert len(msa.records(target.read_text())) == 2
    assert msa.records(target.with_name('Z.a3m').read_text())[0][1] == binder_sequence
    previews = results.manifest_alignment_files(tmp_path / 'job', {'raw_metadata': {'batch_item_id': 'complex1'}})
    assert len(previews) == 2


@pytest.mark.parametrize('family', ['af3_jax_cofolding', 'esmfold2_cofolding'])
def test_family_workflows_default_to_reference_then_colabfold_before_prediction(tmp_path, reference, monkeypatch, family):
    from importlib import import_module
    from types import SimpleNamespace
    module = import_module('mn_molecular_design_studio.engine.workflow_engine.workflows.' + family)
    cache(reference)
    monkeypatch.setattr(msa, 'fetch_colabfold', forbid_fetch)
    run = tmp_path / 'job'
    run.mkdir()
    captured = []
    monkeypatch.setattr(module, 'create_job', lambda *args, **kwargs: captured.append(kwargs) or SimpleNamespace(run_dir=run))
    monkeypatch.setattr(module, 'update_status', lambda *args, **kwargs: None)
    monkeypatch.setattr(module, 'finish_job', lambda *args, **kwargs: None)
    def stop_before_prediction(**kwargs):
        raise RuntimeError('test: stop after staging; never start prediction')
    monkeypatch.setattr(module, '_prediction_command', stop_before_prediction)
    kwargs = dict(protein_sequences=SEQUENCE, msa_max_depth=2)
    if family == 'af3_jax_cofolding':
        kwargs['accept_alphafold3_terms'] = True
    getattr(module, 'run_' + family)(**kwargs)
    assert captured[0]['params']['msa_source'] == 'reference_then_colabfold'
    manifest = run / 'artifacts/raw/cofolding_msa/manifest.json'
    assert manifest.is_file(), (run / 'stderr.log').read_text()
    paths = results.manifest_alignment_files(run, {})
    assert len(paths) == 1
    assert len(msa.records(paths[0][1].read_text())) == 2
    assert 'test: stop after staging' in (run / 'stderr.log').read_text()


def test_colabfold_fetch_command_is_msa_only_and_validates_exported_af3_json(tmp_path, monkeypatch):
    from types import SimpleNamespace
    def run(command, **kwargs):
        assert '--msa-only' in command and '--af3-json' in command
        assert '--gpus' not in command
        assert kwargs['timeout'] == 3600
        output = next((tmp_path / 'artifacts/raw/cofolding_msa/colabfold_server').glob('*/output'))
        (output / 'query.json').write_text(json.dumps({'sequences': [{'protein': {'sequence': SEQUENCE, 'unpairedMsa': ALIGNMENT}}]}))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(msa.subprocess, 'run', run)
    assert msa.fetch_colabfold(tmp_path, SEQUENCE, 'A') == ALIGNMENT
