"""Portable candidate inputs and verified result collection over trusted peer commands.

Uses the existing submission database, leases, blobs, authorization and artifact
commands. Never reads a remote absolute pathname on the originating machine.
"""
from __future__ import annotations

import copy
import hashlib
import json
import time
import zipfile
from pathlib import Path

from . import bundles, capabilities, config
from .store import Store


def actor(owner_id, project_id):
    from mn_molecular_design_studio import auth
    user = next((u for u in auth.list_users() if u['id'] == owner_id and not u.get('disabled')), None)
    if user is None and auth.is_development_mode() and owner_id == auth.development_user()['id']:
        user = auth.development_user()
    if user is None:
        raise PermissionError('The candidate pipeline account is unavailable.')
    auth.require_project_access(user, project_id)
    from .service import user_record
    return user_record(user)


def pack_request(request, destination):
    """Stage selected input files, recursively rebasing candidate JSONL contents.

    References/checkpoints stay on workers. A source directory is represented
    by just the candidate snapshot, rather than copying a run or its database.
    """
    from mn_molecular_design_studio.runtime import reference_root
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import is_operational_path_key, resolve_stored_path
    reference = reference_root().resolve()
    files, source_paths, checksums, refs = {}, {}, {}, set()
    total = 0
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    def stage(path, base):
        nonlocal total
        path = path.resolve(strict=True)
        if path in source_paths:
            return source_paths[path]
        if not path.is_file():
            raise ValueError('Candidate transfer requires files, not arbitrary directories.')
        relative = f'inputs/{len(source_paths):05d}/{path.name}'
        source_paths[path] = relative
        if path.suffix == '.jsonl':
            content = ''.join(json.dumps(transform(json.loads(line), path.parent)) + '\n'
                              for line in path.read_text().splitlines() if line.strip()).encode()
        else:
            if path.stat().st_size > bundles.MAX_BYTES:
                raise ValueError('Candidate input file exceeds the 512 MB transfer limit.')
            content = path.read_bytes()
        total += len(content)
        if total > bundles.MAX_BYTES:
            raise ValueError('Candidate input bundle exceeds the 512 MB transfer limit; reduce batch size.')
        files[relative] = content
        checksums[relative] = hashlib.sha256(content).hexdigest()
        return relative

    def transform(value, base, key=''):
        if isinstance(value, dict):
            result = {k: transform(v, base, k) for k, v in value.items() if k != 'source_run_dir'}
            if value.get('source_run_dir') and value.get('candidates_jsonl'):
                uri = result['candidates_jsonl']
                if not uri.startswith('cluster-input:///'):
                    raise ValueError('A portable candidate snapshot is required.')
                result['source_run_dir'] = uri.rsplit('/', 1)[0]
            return result
        if isinstance(value, list):
            return [transform(v, base, key) for v in value]
        if not isinstance(value, str) or not value:
            return value
        try:
            path = resolve_stored_path(value, run_dir=base, must_exist=True)
        except (OSError, ValueError):
            path = None
        if path is not None and path.is_file():
            path = path.resolve()
            if path.is_relative_to(reference) and path.suffix not in {'.pdb', '.cif', '.mmcif', '.a3m', '.fasta', '.fa'}:
                relative = path.relative_to(reference).as_posix(); refs.add(relative)
                return 'reference:///' + relative
            return 'cluster-input:///' + stage(path, base)
        if is_operational_path_key(key) and value.startswith('/'):
            raise ValueError(f'Worker-specific path cannot be transferred: {key}')
        return value

    prepared = transform(copy.deepcopy(request), destination.parent)
    with zipfile.ZipFile(destination, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for relative, content in sorted(files.items()):
            archive.writestr(zipfile.ZipInfo(relative, date_time=(1980, 1, 1, 0, 0, 0)), content)
        archive.writestr(zipfile.ZipInfo('input-manifest.json', date_time=(1980, 1, 1, 0, 0, 0)), json.dumps({'schema_version': 1, 'sha256': checksums}, sort_keys=True))
    checksum = bundles.digest_file(destination)
    blob = bundles.blob_path(checksum)
    if not blob.exists():
        blob.write_bytes(destination.read_bytes()); blob.chmod(0o600)
    return prepared, checksum, sorted(refs), total


def submit(request, manifest, submission_key, settings, *, owner_id, project_id, state):
    if config.load()['mode'] not in {'peer', 'coordinator'}:
        raise ValueError('Candidate peer dispatch requires an enrolled peer or coordinator origin.')
    user = actor(owner_id, project_id)
    if 'prepared_envelope' not in state:
        prepared, checksum, references, size = pack_request(request, Path(manifest).with_suffix('.inputs.zip'))
        needs = capabilities.requirements(request, references=references)
        from .gpu_runtime import profiles, request_images
        needs['gpu_runtime_profiles'] = profiles(request_images(request))
        state['prepared_envelope'] = dict(kind='workflow', payload=prepared, bundle_id=checksum,
            project_id=project_id, requirements=needs, origin_node=config.identity(),
            placement=settings['placement'], idempotency_key=submission_key)
        state['input_transfer_bytes'] = size
    envelope = state['prepared_envelope']
    with Store().connect() as db:
        db.execute('CREATE TABLE IF NOT EXISTS blob_grants(digest TEXT,owner_id TEXT,PRIMARY KEY(digest,owner_id))')
        db.execute('INSERT OR IGNORE INTO blob_grants VALUES (?,?)', (envelope['bundle_id'], owner_id))
    from .service import accept_submission
    return accept_submission(envelope, user)['job']


def _artifact(state, job_id, relative, user, directory):
    """Resume one authenticated artifact command and verify its returned blob."""
    store = Store()
    transfers = state.setdefault('transfers', {})
    transfer = transfers.setdefault(relative, {})
    if 'command_id' not in transfer:
        transfer.update(store.action(job_id, 'artifact', {'path': relative}, user))
        return None
    command = store.command(transfer['command_id'], user)
    if command['status'] == 'failed':
        raise ValueError('Remote artifact transfer failed: ' + str(command.get('error')))
    if command['status'] != 'completed':
        return None
    checksum = command['response'].get('blob_id')
    blob = bundles.blob_path(checksum)
    if not blob.is_file():
        # Existing peer artifact commands upload blobs to their originating
        # store. An expired copy can be requested again without rerunning inference.
        transfer.pop('command_id', None)
        return None
    if bundles.digest_file(blob) != checksum:
        raise ValueError('Remote output failed its SHA-256 integrity check.')
    relative_path = Path(relative)
    directory = Path(directory).resolve()
    path = (directory / relative_path).resolve()
    if relative_path.is_absolute() or not path.is_relative_to(directory):
        raise ValueError('Remote artifact path escapes the candidate collection.')
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.is_file() or bundles.digest_file(path) != checksum:
        temporary = path.with_suffix(path.suffix + '.transfer')
        temporary.write_bytes(blob.read_bytes()); temporary.replace(path)
    transfer.update(sha256=checksum, bytes=path.stat().st_size)
    return path


def poll(batch, *, owner_id, project_id):
    """Return a local verified snapshot only after every candidate file arrives."""
    user = actor(owner_id, project_id)
    store = Store()
    job = store.job(batch['cluster_job_id'], user)
    batch['worker'] = job.get('node_id')
    status = job.get('status', 'queued')
    if status in {'failed', 'cancelled', 'lost'}:
        return status, None, job
    if status != 'completed':
        node = next((n for n in store.nodes() if n['id'] == job.get('node_id')), None)
        if node and node['availability'] != 'online':
            batch['wait_reason'] = 'Execution workstation disconnected; retaining its job identity until reconciliation.'
            return 'disconnected', None, job
        return status, None, job
    directory = Path(batch['manifest']).with_suffix('.collected')
    collection = 'artifacts/normalized_candidates/candidates.jsonl'
    paths = {name: _artifact(batch, batch['cluster_job_id'], name, user, directory)
             for name in ['metadata.json', 'input.json', 'result.json', collection]}
    if not all(paths.values()):
        return 'collecting', None, job
    metadata = json.loads(paths['metadata.json'].read_text())
    result = json.loads(paths['result.json'].read_text())
    if metadata.get('status') != 'completed' or result.get('success') is not True:
        raise ValueError('Remote job completion disagrees with its scientific artifacts.')
    rows = [json.loads(line) for line in paths[collection].read_text().splitlines() if line.strip()]
    missing = set()
    remote_root = str(job.get('location') or '').rstrip('/')
    # Portable paths belong to the execution worker; resolve only run-relative
    # files from its artifact catalog. Never interpret foreign absolute paths here.
    artifacts = {str(a.get('relative_path') or a.get('path') or '').removeprefix(remote_root + '/') for a in job.get('artifacts', [])}
    def rewrite(value, key=''):
        if isinstance(value, dict):
            return {k: rewrite(v, k) for k, v in value.items()}
        if isinstance(value, list):
            return [rewrite(v, key) for v in value]
        if not isinstance(value, str):
            return value
        relative = None
        if value.startswith('runs:///'):
            reference = value.removeprefix('runs:///')
            local_ref = str(job.get('local_reference') or '') + '/'
            if reference.startswith(local_ref):
                relative = reference[len(local_ref):]
            else:
                return 'cluster://' + str(job['node_id']) + '/' + reference
        elif value.startswith('run:///'):
            relative = value.removeprefix('run:///')
        elif remote_root and value.startswith(remote_root + '/'):
            relative = value[len(remote_root) + 1:]
        elif value in artifacts:
            relative = value
        elif not value.startswith('/') and value.startswith(('artifacts/', 'outputs/')) and Path(value).suffix:
            relative = value
        if relative:
            if relative in {collection, 'metadata.json', 'input.json', 'result.json'}:
                return str(directory / relative)
            path = _artifact(batch, batch['cluster_job_id'], relative, user, directory)
            if path is None:
                missing.add(relative)
                return value
            return str(path)
        if value.startswith('/'):
            # Retain machine-qualified provenance for non-output paths; never
            # treat an execution worker's accepted-input/cache paths as local.
            return 'worker-path://' + str(job['node_id']) + value
        return value
    rewritten = rewrite(rows)
    if missing:
        return 'collecting', None, job
    # Keep the verified native collection intact and publish a separately rebased one.
    rebased = directory / 'artifacts/normalized_candidates/collected.jsonl'
    rebased.write_text(''.join(json.dumps(r) + '\n' for r in rewritten))
    batch['output_transfer_bytes'] = sum(t.get('bytes', 0) for t in batch.get('transfers', {}).values())
    return 'completed', {'run_dir': str(directory), 'candidates_jsonl': str(rebased),
                         'candidate_count': len(rows)}, job | {'metadata': metadata}


def cancel(batch, *, owner_id, project_id):
    return Store().action(batch['cluster_job_id'], 'cancel', {}, actor(owner_id, project_id))
