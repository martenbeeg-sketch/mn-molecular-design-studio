import json
import zipfile
import hashlib
from pathlib import Path
from types import SimpleNamespace

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


def test_candidate_bundle_leaves_saved_source_on_its_owning_peer(tmp_path, monkeypatch):
    from mn_molecular_design_studio import runtime

    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    monkeypatch.setattr(runtime, 'runs_root', lambda: tmp_path / 'peer-runs')
    source = 'cluster://peer-workstation/refolding-validation/saved'
    candidate = source + '/artifacts/normalized_candidates/candidates.jsonl'
    request = {'workflow': 'analysis.validation_acceptance',
               'parameters': {'source_run_dir': source, 'candidates_jsonl': candidate}}
    archive = tmp_path / 'peer.zip'
    prepared, checksum, references, size = cb.pack_request(request, archive)
    assert prepared == request and not references and size == 0
    with zipfile.ZipFile(archive) as bundle:
        assert bundle.namelist() == ['input-manifest.json']

    restored = bundles.restore(prepared, staging=tmp_path, node_id='peer-workstation')['parameters']
    assert Path(restored['source_run_dir']) == tmp_path / 'peer-runs/refolding-validation/saved'
    assert Path(restored['candidates_jsonl']) == Path(restored['source_run_dir']) / 'artifacts/normalized_candidates/candidates.jsonl'

    request['parameters']['candidates_jsonl'] = candidate.replace('/saved/', '/different-run/')
    with pytest.raises(ValueError, match='owning source run'):
        cb.pack_request(request, archive)


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


def test_poll_summary_reads_counts_without_requesting_candidate_artifacts(monkeypatch):
    summary = {
        "schema_version": 1, "status": "complete", "candidate_count": 4, "eligible_count": 3,
        "structure_count": 2, "complex_count": 1,
        "target_systems": [{"identity": "sequence:" + hashlib.sha256(b"ACDE").hexdigest(),
                            "identity_kind": "sequence", "sequence_length": 4,
                            "label": "Sequence-defined target · 4 aa · " + hashlib.sha256(b"ACDE").hexdigest()[:8]}],
        "native_filter_counts": {"passed": 2, "failed": 1, "unavailable": 0},
        "acceptance_counts": {"accepted": 0, "rejected": 0, "unresolved": 0, "not_evaluated": 4},
    }

    class FakeStore:
        def job(self, job_id, user):
            assert job_id == "cluster-job"
            return {"status": "completed", "node_id": "worker", "candidate_summary": summary}

        def nodes(self):
            raise AssertionError("A completed job summary does not need an artifact or node request.")

    monkeypatch.setattr(cb, "Store", FakeStore)
    monkeypatch.setattr(cb, "actor", lambda *args, **kwargs: {"id": "researcher"})
    batch = {"cluster_job_id": "cluster-job"}

    status, job = cb.poll_summary(batch, owner_id="researcher", project_id="project")

    assert status == "completed"
    assert job["candidate_summary"] == summary
    assert batch["candidate_summary"] == summary
    assert "candidates" not in batch["candidate_summary"]


def test_poll_summary_accepts_count_only_filter_results_without_design_metadata(monkeypatch):
    summary = {
        "schema_version": 1, "status": "complete", "candidate_count": 4, "eligible_count": 4,
        "native_filter_counts": {"passed": 0, "failed": 0, "unavailable": 4},
        "acceptance_counts": {"accepted": 2, "rejected": 2, "unresolved": 0, "not_evaluated": 0},
        "native_filter_vs_acceptance": {
            "passed": {"accepted": 0, "rejected": 0, "unresolved": 0, "not_evaluated": 0},
            "failed": {"accepted": 0, "rejected": 0, "unresolved": 0, "not_evaluated": 0},
            "unavailable": {"accepted": 2, "rejected": 2, "unresolved": 0, "not_evaluated": 0},
        },
        "acceptance_evaluation": "worker_local_saved_evidence",
    }

    class FakeStore:
        def job(self, _job_id, _user):
            return {"status": "completed", "node_id": "worker", "candidate_summary": summary}

        def nodes(self):
            raise AssertionError("Aggregate acceptance counts must not trigger artifact access.")

    monkeypatch.setattr(cb, "Store", FakeStore)
    monkeypatch.setattr(cb, "actor", lambda *args, **kwargs: {"id": "researcher"})
    batch = {"cluster_job_id": "cluster-job"}

    status, job = cb.poll_summary(batch, owner_id="researcher", project_id="project")

    assert status == "completed"
    assert job["candidate_summary"] == summary
    assert set(summary).isdisjoint({"candidates", "accepted_candidate_ids", "candidate_sequences"})


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
    persisted = []
    persist = lambda transfers: persisted.append({key: dict(value) for key, value in transfers.items()})
    assert cb._artifact(*args, persist_transfers=persist) is None
    assert persisted[-1]['artifacts/design.pdb']['command_id'] == 'command'
    assert cb._artifact(*args, persist_transfers=persist).read_bytes() == data
    assert persisted[-1]['artifacts/design.pdb']['sha256'] == checksum
    assert cb._artifact(*args, persist_transfers=persist).read_bytes() == data
    assert len(calls) == 1
    bundles.blob_path(checksum).write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='integrity'):
        cb._artifact(*args)


def test_eligible_campaign_placement_uses_origin_and_live_capacity_options(tmp_path, monkeypatch):
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
    for ordinal, state in ((1, first), (2, second)):
        cb.submit(
            common, str(tmp_path / f'{ordinal}.json'), f'key-{ordinal}',
            {'placement': 'eligible', 'placement_ordinal': ordinal, 'include_origin_node': True},
            owner_id='owner', project_id='project', state=state, persist=lambda: None,
        )

    assert selected == ['origin', 'peer-4032']
    assert first['placement_policy'] == second['placement_policy'] == 'eligible'
    assert first['placement_reference_probe_pending'] is False
    assert second['placement_reference_probe_pending'] is False
    assert first['placement_node_id'] == first['prepared_envelope']['placement'] == 'origin'
    assert second['placement_node_id'] == second['prepared_envelope']['placement'] == 'peer-4032'


def test_eligible_peer_placement_excludes_origin_by_default(monkeypatch):
    class Store:
        def placement_options(self, _needs, *, origin_node):
            assert origin_node == 'origin'
            return [{'node_id': 'origin'}, {'node_id': 'peer'}]

    monkeypatch.setattr(cb, 'Store', Store)
    assert cb.select_eligible_peer({}, origin_node='origin')['node_id'] == 'peer'


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
                  'chaperone_handoffs': [{'handler_metrics': {'studio': {'ckpt_path': '/worker/old/checkpoint'}}}],
                  'effective_settings': {'checkpoint': '/reference/mol1/checkpoint/mol-1-base.pt',
                                         'need_atom_confidence': True},
                  'scientific_contract': {'source_structure': '/worker/old/staged.cif', 'schema_version': 1},
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
    assert row['raw_metadata']['effective_settings'] == provenance['effective_settings']
    assert row['raw_metadata']['candidate_set'] == provenance['candidate_set']
    assert row['raw_metadata']['chaperone_handoffs'] == provenance['chaperone_handoffs']
    assert row['raw_metadata']['scientific_contract'] == provenance['scientific_contract']
    assert row['metrics'] == {'ckpt_path': 'runs:///design/remote/nvidia/RFD4-Proteina', 'iptm': 0.7}
    new_prediction = {'workflow': 'refolding.optimized_prediction',
                      'parameters': {'tasks': [{'scientific_contract': provenance['scientific_contract']}]}}
    with pytest.raises(ValueError, match='Worker-specific path cannot be transferred: source_structure'):
        cb.pack_request(new_prediction, tmp_path / 'invalid.zip')


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
    monkeypatch.setattr(cb, '_artifact', lambda batch, job, path, user, destination, **_kwargs: directory / path)
    outcome, output, job = cb.poll({'cluster_job_id': 'job', 'manifest': str(tmp_path / 'batch.json')}, owner_id='owner', project_id=None)
    assert outcome == ('completed' if status == 'completed' else 'partial')
    row = json.loads(Path(output['candidates_jsonl']).read_text())
    assert row['candidate_id'] == 'exact-member' and row['parents'] == ['original']
    assert job['native_execution']['native_invocations'] == 1
    assert Path(row['complex_pdb']).read_text() == 'ATOM'


def test_remote_ensemble_collection_fetches_receipts_without_fetching_trajectories(tmp_path, monkeypatch):
    batch = {"cluster_job_id": "ensemble", "manifest": str(tmp_path / "batch.json"), "unit_count": 1}
    xtc = "artifacts/bioemu2/system/samples.xtc"
    result = {"success": True, "metrics": {"systems": 1}}
    files = {name: tmp_path / name for name in ("metadata.json", "result.json")}
    files["metadata.json"].write_text(json.dumps({"status": "completed"}))
    files["result.json"].write_text(json.dumps(result))
    catalog = [*files, xtc, "artifacts/bioemu2/system/topology.cif"]
    job = {"status": "completed", "node_id": "worker",
           "artifacts": [{"relative_path": name} for name in catalog]}
    monkeypatch.setattr(cb, "Store", lambda: SimpleNamespace(job=lambda *_: job))
    monkeypatch.setattr(cb, "actor", lambda *args, **kwargs: {})
    monkeypatch.setattr(cb, "_artifact", lambda batch, job, path, *args, **kwargs: files[path])
    status, output, _ = cb.poll(batch, owner_id="owner", project_id="project", candidate_collection=None)
    assert status == "completed"
    from mn_molecular_design_studio.candidate_execution import collect_ensemble
    collect_ensemble(batch, output["result"], output["result_path"], output["artifact_paths"], "cluster:ensemble")
    assert batch["tile_run"]["trajectory_paths"] == [xtc]
    assert batch["tile_run"]["job_reference"] == "cluster:ensemble"
    assert batch["completed_units"] == 1
    with pytest.raises(ValueError, match="paired trajectory and topology"):
        collect_ensemble(batch, result, files["result.json"], [xtc], "cluster:ensemble")


def test_poll_collects_the_exact_structure_from_the_authorized_source_job(tmp_path, monkeypatch):
    from contextlib import nullcontext

    directory = tmp_path / 'batch.collected'
    collection = 'artifacts/normalized_candidates/candidates.jsonl'
    source_ref = 'refolding-validation/source'
    structure = 'artifacts/output.cif'
    pae = 'artifacts/prediction_evidence/model/1/model_pae.json'
    row = {'candidate_id': 'native-id', 'complex_pdb': f'runs:///{source_ref}/{structure}',
           'raw_metadata': {'source_run_dir': f'runs:///{source_ref}',
                            'prediction_pae_measurement': {'source': pae},
                            'input_provenance': f'runs:///{source_ref}/execution_inputs.json'}}
    for name, value in {'metadata.json': {'status': 'completed'}, 'input.json': {},
                        'result.json': {'success': True}, collection: row}.items():
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value) + '\n')

    class Store:
        def connect(self):
            return nullcontext(self)

        def execute(self, query, bindings):
            assert bindings == ('peer', source_ref)
            return self

        def fetchone(self):
            return {'id': 'source-job'}

        def job(self, job_id, _user):
            if job_id == 'source-job':
                return {'project_id': 'project'}
            return {'status': 'completed', 'node_id': 'peer', 'local_reference': 'analysis/child'}

    transfers = []

    def artifact(state, job_id, relative, _user, destination, **_kwargs):
        if job_id == 'source-job':
            assert relative in {structure, pae} and destination == directory / 'source'
            transfers.append(relative)
            if transfers.count(relative) == 1:
                return None
            path = destination / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('data_exact_source\n')
            state.setdefault('transfers', {})[relative] = {'bytes': path.stat().st_size}
            return path
        return directory / relative

    monkeypatch.setattr(cb, 'Store', Store)
    monkeypatch.setattr(cb, 'actor', lambda *a, **k: {})
    monkeypatch.setattr(cb, '_artifact', artifact)
    batch = {'cluster_job_id': 'child-job', 'manifest': str(tmp_path / 'batch.json'),
             'prepared_envelope': {'payload': {'parameters': {'source_run_dir': f'cluster://peer/{source_ref}'}}}}
    assert cb.poll(batch, owner_id='owner', project_id='project')[:2] == ('collecting', None)
    status, output, _ = cb.poll(batch, owner_id='owner', project_id='project')
    assert status == 'completed'
    fetched = json.loads(Path(output['candidates_jsonl']).read_text())
    assert fetched['candidate_id'] == 'native-id'
    assert Path(fetched['complex_pdb']).read_text() == 'data_exact_source\n'
    assert Path(fetched['raw_metadata']['prediction_pae_measurement']['source']).read_text() == 'data_exact_source\n'
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates

    normalized = read_candidates(Path(output['candidates_jsonl']))[0]
    assert normalized['raw_metadata']['source_run_dir'] == f'worker-path://peer/{source_ref}'
    assert normalized['raw_metadata']['input_provenance'] == f'worker-path://peer/{source_ref}/execution_inputs.json'
    normalized['raw_metadata']['candidate_execution'] = {
        'effective_request': batch['prepared_envelope']['payload']}
    replay = tmp_path / 'replay.jsonl'
    replay.write_text(json.dumps(normalized) + '\n')
    assert read_candidates(replay)[0]['raw_metadata']['candidate_execution'] == normalized['raw_metadata']['candidate_execution']
    assert batch['output_transfer_bytes'] == 2 * len(b'data_exact_source\n')
