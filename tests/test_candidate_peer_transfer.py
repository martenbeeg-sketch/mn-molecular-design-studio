import json
import zipfile
from pathlib import Path

import pytest

from mn_molecular_design_studio.cluster import bundles, candidate_batches as cb


def test_portable_candidate_snapshot_moves_only_selected_inputs_and_verifies_hashes(tmp_path, monkeypatch):
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR', str(tmp_path / 'reference'))
    source = tmp_path / 'source'; source.mkdir()
    structure = source / 'binder.pdb'; structure.write_text('ATOM')
    (source / 'large-unrelated-database').write_bytes(b'not submitted')
    collection = source / 'selected.jsonl'
    collection.write_text(json.dumps({'candidate_id': 'stable', 'parents': ['original'], 'complex_pdb': str(structure),
                                     'raw_metadata': {'fixed_residues': ['B1'], 'seed': 42}}) + '\n')
    request = dict(workflow='sequence.mpnn', parameters={'source_run_dir': str(source), 'candidates_jsonl': str(collection), 'seed': 42})
    archive = tmp_path / 'snapshot.zip'
    prepared, checksum, refs, size = cb.pack_request(request, archive)
    assert not refs
    with zipfile.ZipFile(archive) as z:
        assert not any('database' in name for name in z.namelist())
    destination = tmp_path / 'other-workstation'
    bundles.unpack(archive, destination, expected_digest=checksum)
    restored = bundles.restore(prepared, staging=destination, node_id='other')
    p = restored['parameters']
    assert Path(p['source_run_dir']).is_dir()
    moved = json.loads(Path(p['candidates_jsonl']).read_text())
    assert Path(moved['complex_pdb']).read_text() == 'ATOM'
    assert moved['candidate_id'] == 'stable' and moved['parents'] == ['original']
    assert moved['raw_metadata']['fixed_residues'] == ['B1']
    assert not str(source) in json.dumps(prepared)
    # Blob integrity and per-file manifest integrity are independently enforced.
    archive.write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='integrity'):
        bundles.unpack(archive, tmp_path / 'bad', expected_digest=checksum)


def test_candidate_bundle_is_content_stable_and_does_not_copy_checkpoints(tmp_path, monkeypatch):
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    reference = tmp_path / 'reference'; reference.mkdir()
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR', str(reference))
    weights = reference / 'weights.ckpt'; weights.write_bytes(b'checkpoint')
    request = dict(workflow='mock', parameters={'model_weights': str(weights), 'seed': 17})
    first = cb.pack_request(request, tmp_path / 'a.zip')
    second = cb.pack_request(request, tmp_path / 'b.zip')
    assert first == second
    assert first[0]['parameters']['model_weights'] == 'reference:///weights.ckpt'
    assert first[2] == ['weights.ckpt']
    assert first[3] == 0


def test_output_transfer_resumes_and_rejects_corruption(tmp_path, monkeypatch):
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    data = b'scientific artifact'
    import hashlib
    checksum = hashlib.sha256(data).hexdigest()
    bundles.blob_path(checksum).write_bytes(data)
    calls = []
    class Store:
        def action(self, *args):
            calls.append(args); return dict(command_id='command')
        def command(self, *args):
            return dict(status='completed', response={'blob_id': checksum})
    monkeypatch.setattr(cb, 'Store', Store)
    state = {}
    args = (state, 'job', 'artifacts/design.pdb', {}, tmp_path / 'result')
    assert cb._artifact(*args) is None
    assert cb._artifact(*args).read_bytes() == data
    assert cb._artifact(*args).read_bytes() == data
    assert len(calls) == 1
    bundles.blob_path(checksum).write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='integrity'):
        cb._artifact(*args)


def test_cuda_targets_distinguish_4090_and_5090_without_assuming_equal_support():
    from mn_molecular_design_studio.cluster.gpu_runtime import compatible
    old = {'torch-image': {'cuda': '12.1', 'targets': ['sm_80', 'sm_86', 'sm_90']}}
    assert compatible({'compute_capability': '8.9', 'driver_cuda_api': 13000}, old)
    assert not compatible({'compute_capability': '12.0'}, old)
    assert not compatible({}, old)
    assert compatible({'compute_capability': '12.0'}, {'new': {'targets': ['sm_120']}})
    assert compatible({'compute_capability': '12.0'}, {'ptx': {'targets': ['compute_90']}})


def test_peer_prediction_does_not_require_or_transfer_upstream_models(tmp_path, monkeypatch):
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    reference = tmp_path / 'reference'; reference.mkdir()
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR', str(reference))
    weights = reference / 'upstream-generator.ckpt'; weights.write_bytes(b'prior model')
    structure = tmp_path / 'candidate.pdb'; structure.write_text('ATOM')
    old_collection = tmp_path / 'prior-full-collection.jsonl'; old_collection.write_text('{}\n')
    provenance = {'checkpoint_provenance': {'path': str(weights), 'sha256': bundles.digest_file(weights)},
                  'candidate_execution': {'effective_request': {'candidates_jsonl': str(old_collection)}}}
    collection = tmp_path / 'selected.jsonl'
    collection.write_text(json.dumps({'candidate_id': 'selected', 'complex_pdb': str(structure), 'raw_metadata': provenance})+'\n')
    request = {'workflow': 'mock', 'parameters': {'source_run_dir': str(tmp_path), 'candidates_jsonl': str(collection)}}
    archive = tmp_path / 'inputs.zip'
    prepared, checksum, refs, size = cb.pack_request(request, archive)
    assert refs == []
    with zipfile.ZipFile(archive) as z:
        assert not any('upstream-generator' in n or 'prior-full-collection' in n for n in z.namelist())
    destination = tmp_path / 'worker'
    bundles.unpack(archive, destination, expected_digest=checksum)
    restored = bundles.restore(prepared, staging=destination, node_id='worker')
    row = json.loads(Path(restored['parameters']['candidates_jsonl']).read_text())
    for key, original in provenance.items():
        assert json.loads(row['raw_metadata'][key]['record']) == original
