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


def test_eligible_campaign_placement_pins_distinct_compatible_peers(tmp_path, monkeypatch):
    class Connection:
        def __enter__(self):
            return self
        def __exit__(self, *_args):
            return False
        def execute(self, *_args):
            return None

    options = [
        {'node_id': 'origin'},
        {'node_id': 'peer-4032'},
        {'node_id': 'peer-4033'},
    ]

    class Store:
        def connect(self):
            return Connection()
        def placement_options(self, needs, *, origin_node):
            assert needs['gpu_count'] == 1
            assert needs['gpu_runtime_profiles'] == {}
            assert origin_node == 'origin'
            return options

    monkeypatch.setattr(cb.config, 'load', lambda: {'mode': 'coordinator'})
    monkeypatch.setattr(cb.config, 'identity', lambda: 'origin')
    monkeypatch.setattr(cb, 'Store', Store)
    monkeypatch.setattr(cb, 'pack_request', lambda *_args: ({'workflow': 'mock'}, 'bundle', [], 7))
    monkeypatch.setattr(cb.capabilities, 'requirements', lambda *_args, **_kwargs: {'gpu_count': 1})
    from mn_molecular_design_studio.cluster import gpu_runtime, service
    monkeypatch.setattr(gpu_runtime, 'profiles', lambda _images: {})
    monkeypatch.setattr(cb, 'actor', lambda *_args, **_kwargs: {'id': 'owner'})
    selected = []

    def accept(envelope, _user):
        selected.append(envelope['placement'])
        return {'job': {'cluster_job_id': f'job-{len(selected)}', 'node_id': envelope['placement']}}

    monkeypatch.setattr(service, 'accept_submission', accept)
    first = {}
    second = {}
    common = {
        'workflow': 'mock', 'parameters': {},
    }
    for ordinal, state, exclusions in (
        (1, first, []),
        (2, second, [first.get('placement_node_id')]),
    ):
        cb.submit(
            common, str(tmp_path / f'{ordinal}.json'), f'key-{ordinal}',
            {'placement': 'eligible', 'placement_ordinal': ordinal, 'exclude_node_ids': exclusions},
            owner_id='owner', project_id='project', state=state, persist=lambda: None,
        )

    assert selected == ['peer-4032', 'peer-4033']
    assert first['placement_policy'] == second['placement_policy'] == 'eligible'
    assert first['placement_reference_probe_pending'] is False
    assert second['placement_reference_probe_pending'] is False
    assert first['placement_node_id'] == first['prepared_envelope']['placement'] == 'peer-4032'
    assert second['placement_node_id'] == second['prepared_envelope']['placement'] == 'peer-4033'


def test_eligible_placement_waits_for_exact_reference_probe_before_dispatch(monkeypatch):
    requirements_seen = []

    class Store:
        def placement_options(self, needs, *, origin_node):
            assert origin_node == 'origin'
            requirements_seen.append(dict(needs))
            if needs.get('references') or needs.get('reference_hashes'):
                return []
            return [{'node_id': 'peer-4032'}]

    monkeypatch.setattr(cb, 'Store', Store)
    requirements = {
        'code_revision': 'reviewed-source',
        'references': ['rfd4/model.safetensors'],
        'reference_hashes': {'rfd4/model.safetensors': 'expected-sha256'},
    }

    selection = cb.select_eligible_peer(requirements, origin_node='origin')

    assert selection == {'node_id': 'peer-4032', 'reference_probe_pending': True}
    assert requirements_seen[0] == requirements
    assert requirements_seen[1] == {
        **requirements, 'references': [], 'reference_hashes': {},
    }


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
                  'candidate_execution': {'effective_request': {'candidates_jsonl': str(old_collection)}},
                  'candidate_set': {'native_measurements': {'ckpt_path': '/worker/old/checkpoint'},
                                    'source_provenance': {'uri': 'cluster://old-node/design/run/result.json'}}}
    collection = tmp_path / 'selected.jsonl'
    collection.write_text(json.dumps({'candidate_id': 'selected', 'complex_pdb': str(structure), 'raw_metadata': provenance, 'metrics': {'ckpt_path': 'runs:///design/remote/nvidia/RFD4-Proteina', 'iptm': 0.7}})+'\n')
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
    for key in ('checkpoint_provenance', 'candidate_execution'):
        assert json.loads(row['raw_metadata'][key]['record']) == provenance[key]
    assert row['raw_metadata']['candidate_set'] == provenance['candidate_set']
    assert row['metrics'] == {'ckpt_path': 'runs:///design/remote/nvidia/RFD4-Proteina', 'iptm': 0.7}


@pytest.mark.parametrize('status', ['completed', 'failed', 'lost'])
def test_poll_preserves_native_receipts_and_completed_partial_members(tmp_path, monkeypatch, status):
    directory = tmp_path / 'batch.collected'
    directory.mkdir()
    collection = 'artifacts/normalized_candidates/candidates.jsonl'
    content = {'metadata.json': {'status': status}, 'input.json': {},
               collection: {'candidate_id': 'exact-member', 'parents': ['original'], 'complex_pdb': 'artifacts/output.pdb'},
               'native_execution.json': {'native_invocations': 1, 'docker_processes': 1}}
    if status != 'lost': content['result.json'] = {'success': status == 'completed'}
    for name, value in content.items():
        path = directory / name; path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value) + '\n')
    (directory / 'artifacts/output.pdb').write_text('ATOM')
    remote = dict(status=status, node_id='worker', artifacts=[{'relative_path': p} for p in content])
    class Store:
        def job(self, *_): return remote
    monkeypatch.setattr(cb, 'Store', Store)
    monkeypatch.setattr(cb, 'actor', lambda *a, **k: {})
    monkeypatch.setattr(cb, '_artifact', lambda batch, job, path, user, destination: directory / path)
    outcome, output, job = cb.poll({'cluster_job_id': 'job', 'manifest': str(tmp_path / 'batch.json')}, owner_id='owner', project_id=None)
    assert outcome == ('completed' if status == 'completed' else 'partial')
    row = json.loads(Path(output['candidates_jsonl']).read_text())
    assert row['candidate_id'] == 'exact-member' and row['parents'] == ['original']
    assert job['native_execution']['native_invocations'] == 1
    assert Path(row['complex_pdb']).read_text() == 'ATOM'
