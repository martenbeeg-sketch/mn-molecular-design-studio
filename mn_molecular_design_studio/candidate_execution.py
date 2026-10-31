"""Durable candidate batches for the existing graph controller and job scheduler.

No execution service lives here: each attempt is an ordinary Studio workflow job.
Saved graphs without an execution policy retain their native collection behavior.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import re
from itertools import islice
from pathlib import Path

from mn_molecular_design_studio.cache_provenance import file_digest, referenced_files
from mn_molecular_design_studio.cluster.config import atomic_json
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, utc_now

# Only adapters whose unit of work is an independent candidate are enabled.
COLLECTION_KINDS = {
    'ligandmpnn', 'foundry_mpnn', 'abmpnn', 'af2', 'af3', 'rf3', 'boltz2',
    'protenix', 'ipsae', 'rosetta_metrics',
}
TERMINAL = {'completed', 'partial', 'failed', 'cancelled', 'lost'}


def policy(value=None):
    value = value or {}
    if not isinstance(value, dict):
        raise ValueError('Candidate execution must be an object.')
    defaults = dict(mode='native', placement='local', batch_size=4, max_in_flight=2,
                    retries=1, failure_policy='all', min_vram_gb=8., ram_gb=0.,
                    cpu_cores=4, stall_timeout_seconds=1800, incremental_cpu=False)
    optional = {'gpu_device', 'max_runtime_seconds', 'deadline_at'}
    if set(value) - defaults.keys() - optional:
        raise ValueError('Unsupported candidate execution settings: ' + ', '.join(sorted(set(value) - defaults.keys() - optional)))
    result = defaults | value
    if any(value.get(key) is not None for key in ('max_runtime_seconds', 'deadline_at')):
        from mn_molecular_design_studio.campaign_objectives import normalize_objective
        if result['mode'] != 'sharded':
            raise ValueError('Dispatch time limits require partitioned batches; a native invocation cannot be split by this controller.')
        limits = normalize_objective({'objective_type': 'generated', **{
            key: value[key] for key in ('max_runtime_seconds', 'deadline_at') if value.get(key) is not None}})
        result.update({key: limits[key] for key in ('max_runtime_seconds', 'deadline_at') if key in limits})
    if 'gpu_device' in value:
        from mn_molecular_design_studio.engine.workflow_engine.core.gpu import normalize_gpu_device
        selector = value['gpu_device']
        if not isinstance(selector, str) or selector not in {'auto', 'all'} and not re.fullmatch(r'\d+(?:,\d+)*', selector):
            raise ValueError('Choose auto, all, or explicit GPU indices for this processing step.')
        normalize_gpu_device(selector)
    if result['mode'] not in {'native', 'sharded'} or result['failure_policy'] not in {'all', 'partial'}:
        raise ValueError('Choose native/sharded execution and all/partial failure policy.')
    if type(result['incremental_cpu']) is not bool:
        raise ValueError('incremental_cpu must be true or false.')
    import uuid
    if result['placement'] not in {'auto', 'local'}:
        uuid.UUID(str(result['placement']))
    for key, low, high in [('batch_size', 1, 10000), ('max_in_flight', 1, 32), ('retries', 0, 5),
                           ('cpu_cores', 1, 256), ('stall_timeout_seconds', 60, 86400)]:
        if type(result[key]) is not int or not low <= result[key] <= high:
            raise ValueError(f'{key} must be an integer in [{low}, {high}].')
    for key in ('min_vram_gb', 'ram_gb'):
        if isinstance(result[key], bool) or not isinstance(result[key], (int, float)) or not math.isfinite(result[key]) or result[key] < 0:
            raise ValueError(f'{key} must be finite and nonnegative.')
    return result


def supports(kind, source):
    if source.get("run_sources"):
        return True  # The workflow contract is validated by source_request.
    from mn_molecular_design_studio.pipeline_graph import OPTIMIZED_PREDICTION_NODE_ENGINES
    from mn_molecular_design_studio.quality_tools import WORKFLOWS
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.inverse_folding_input import ENGINES as INVERSE_ENGINES
    workflow = NODE_CATALOG[kind].get('workflow')
    from mn_molecular_design_studio.workflow_modes import capability
    generator = capability(kind)
    return ((bool(source.get('candidates_jsonl')) and
             (kind in COLLECTION_KINDS or kind in OPTIMIZED_PREDICTION_NODE_ENGINES or kind in INVERSE_ENGINES or workflow in WORKFLOWS))
            or generator.get('generator_only', False))


def dispatch_limits(kind, settings):
    """New dispatch limits apply to independent work; preserve existing cohort operations."""
    from mn_molecular_design_studio.quality_tools import TOOLS
    if TOOLS.get(kind, {}).get('cohort'):
        return {}
    return {key: settings[key] for key in ('max_runtime_seconds', 'deadline_at') if settings.get(key) is not None}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def runs_with_submission_key(root, submission_key):
    """Find local run receipts without resolving unrelated managed artifacts.

    Run metadata may include portable paths owned by another workstation. This
    lookup only needs the durable submission key, so parsing raw metadata keeps
    cluster URIs from unrelated runs out of the local recovery path.
    """
    matches = []
    for path in root.glob('*/*/metadata.json'):
        try:
            metadata = json.loads(path.read_text())
        except (OSError, UnicodeError, json.JSONDecodeError):
            continue
        if isinstance(metadata, dict) and metadata.get('candidate_submission_key') == submission_key:
            matches.append(path.parent)
    return matches


def data_hashes(value):
    """Hash data artifacts without repeatedly streaming model/reference stores.

    Installed model identity belongs to the resolved execution specification and
    native checkpoint provenance. It is not an output artifact of each candidate.
    Scientific reference structures/MSAs still receive ordinary content checks.
    """
    from mn_molecular_design_studio.runtime import reference_root
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    root = reference_root().resolve()
    seen, result = set(), {}
    scientific = {'.pdb', '.cif', '.mmcif', '.a3m', '.fasta', '.fa', '.sdf', '.mol2', '.pdbqt'}
    def visit(item, base=None):
        if isinstance(item, dict):
            for key, child in item.items():
                if key not in {'candidate_execution', 'checkpoint_provenance'}:
                    visit(child, base)
        elif isinstance(item, list):
            for child in item:
                visit(child, base)
        elif isinstance(item, (str, Path)) and item:
            try:
                path = resolve_stored_path(str(item), run_dir=base, must_exist=True)
                if path is None or not path.is_file():
                    return
                path = path.resolve()
                if path in seen or path.is_relative_to(root) and path.suffix.lower() not in scientific:
                    return
            except (OSError, ValueError):
                return
            seen.add(path)
            result[str(path)] = file_digest(path)
            try:
                if path.suffix == '.json':
                    visit(json.loads(path.read_text()), path.parent)
                elif path.suffix == '.jsonl':
                    with path.open() as handle:
                        for line in handle:
                            if line.strip():
                                visit(json.loads(line), path.parent)
            except (ValueError, UnicodeError):
                pass
    visit(value)
    return result


def initialize(parent, node, source, params, key, settings, *, make_source, make_request, read_rows):
    """Freeze ordered units and complete effective requests before submitting any work."""
    if source.get("run_sources"):
        from mn_molecular_design_studio.run_graph import initialize_batches
        return initialize_batches(parent, node, source, params, key, settings)
    collection = bool(source.get('candidates_jsonl'))
    native_invocation = node.get('native_invocation', False)
    from mn_molecular_design_studio.workflow_modes import capability
    route = capability(node['data']['kind'])
    generation = (route.get('generator_only', False) or native_invocation) and not collection
    count_field = route.get('design_count_field') or 'num_designs'
    total = 1 if native_invocation else int(params.get(count_field, 1)) if generation else source.get('candidate_count')
    if source.get('candidate_set_manifest'):
        from mn_molecular_design_studio.candidate_sets import iter_candidates
        iterator = iter_candidates(source)
    else:
        rows = read_rows(source) if collection else []
        iterator = iter(rows)
        if not generation:
            total = len(rows)
    if total is None:
        raise ValueError('Candidate collection must declare its count before sharding.')
    if total <= 0:
        return []
    directory = Path(parent) / 'candidate_batches' / node['id'] / key
    directory.mkdir(parents=True, exist_ok=True)
    batches = []
    seen = set()
    seed = int(params.get('seed') or 0)
    if not seed:
        # Persisted scientific key -> stable positive seed on replay/retry.
        seed = int(key[:8], 16) % 2147483646 + 1
    for offset in range(0, total, settings['batch_size']):
        count = min(settings['batch_size'], total - offset)
        batch_id = f'b{offset:06d}'
        selected = list(islice(iterator, count)) if collection else []
        ids = [str(r.get('candidate_id') or '') for r in selected]
        if collection and (len(selected) != count or not all(ids) or len(set(ids)) != len(ids) or seen.intersection(ids)):
            raise ValueError('Candidate sharding requires nonempty unique input candidate IDs and an exact collection count.')
        seen.update(ids)
        batch_source = make_source(parent, f'{node["id"]}_{batch_id}_inputs', selected, source.get('target')) if collection else source
        effective = copy.deepcopy(params)
        if generation and not native_invocation:
            effective[count_field] = count
            # Never repeat the same trajectory seed on separate shards.
            effective['seed'] = (seed + offset) % 2147483646 or 1
            if 'batch_size' in effective:
                effective['batch_size'] = min(int(effective['batch_size']), count)
        elif not native_invocation:
            effective['seed'] = seed if 'seed' in params else params.get('seed')
            if effective['seed'] is None:
                effective.pop('seed')
        request = make_request(node['data']['kind'], effective, batch_source, f'{node["data"]["label"]} · {batch_id}')
        # Bind the node's CPU choice before freezing the request. Prediction
        # adapters may supply a default, which must agree with scheduler resources.
        if 'cpu_cores' in request.get('parameters', {}):
            request['parameters']['cpu_cores'] = settings['cpu_cores']
        manifest = dict(schema_version=1, pipeline=str(parent), node=node['id'], batch_id=batch_id,
                        scientific_key=key, offset=offset, input_candidate_ids=ids,
                        unit_count=count, request=request, input_hashes=data_hashes(request), policy=settings)
        if native_invocation:
            manifest['unit_kind'] = 'native_invocation'
        if request['workflow'] == 'refolding.optimized_prediction':
            prediction = request['parameters']
            manifest['expected_outputs_per_input'] = int(prediction.get('samples', 1)) * int(prediction.get('engine_options', {}).get('seed_count', 1))
            manifest['expected_output_count'] = count * manifest['expected_outputs_per_input']
        if source.get('candidate_set_manifest'):
            manifest['candidate_set_schema'] = 'studio.candidate-set.v1'
        manifest['sha256'] = digest(manifest)
        path = directory / f'{batch_id}.json'
        if path.exists():
            if json.loads(path.read_text()) != manifest:
                raise ValueError('Immutable candidate batch inputs changed; submit a new pipeline.')
        else:
            atomic_json(path, manifest)
        batches.append(dict(batch_id=batch_id, offset=offset, unit_count=count, manifest=str(path),
                            manifest_sha256=manifest['sha256'], status='pending', attempts=[]))
    if collection and next(iterator, None) is not None:
        raise ValueError('Candidate collection contains more records than its declared count.')
    return batches


def read_manifest(batch):
    value = json.loads(Path(batch['manifest']).read_text())
    checksum = value.pop('sha256')
    if checksum != batch['manifest_sha256'] or digest(value) != checksum:
        raise ValueError('Candidate input manifest failed its SHA-256 integrity check.')
    for path, expected in value['input_hashes'].items():
        if not Path(path).is_file() or file_digest(path) != expected:
            raise ValueError(f'Candidate input changed after submission: {path}')
    value['sha256'] = checksum
    return value


def output_valid(batch):
    manifest = batch.get('output_hashes') or {}
    return bool(manifest) and all(Path(p).is_file() and file_digest(p) == h for p, h in manifest.items())


def collect(batch, output, *, read_rows):
    """Assign output identities once, independently of attempt/worker/native run IDs."""
    manifest = read_manifest(batch)
    rows = read_rows({**output, 'candidates_jsonl': output['review_candidates_jsonl']}
                     if manifest.get('unit_kind') == 'native_invocation' and output.get('review_candidates_jsonl') else output)
    for row in rows:
        structure = row.get('complex_pdb') or row.get('binder_pdb')
        if structure and not Path(structure).is_file():
            raise ValueError('Worker returned an unavailable candidate structure.')
    native_ids = [str(r.get('candidate_id') or '') for r in rows]
    if not all(native_ids) or len(set(native_ids)) != len(native_ids):
        raise ValueError('Worker returned missing or duplicated candidate identities.')
    inputs = set(manifest['input_candidate_ids'])
    # Native IDs may include timestamps. The immutable batch and native output
    # order define identity; retain every original identity in provenance.
    mapping = {v: (v if v in inputs else f'{manifest["node"]}:{batch["batch_id"]}:{i:06d}')
               for i, v in enumerate(native_ids)}
    if manifest.get('candidate_set_schema'):
        # Native result ordering is not scientific identity. Run-local prefixes
        # vary on retries; tool-specific output/sample suffixes remain explicit.
        mapping = {native_id: native_id if native_id in inputs else 'cs_output_' + digest([
            manifest['node'], manifest['scientific_key'], sorted(row.get('parents') or []),
            native_id.removeprefix(Path(output['run_dir']).name)]) for row, native_id in zip(rows, native_ids)}
    attempt = batch['attempts'][-1]
    for row, native_id in zip(rows, native_ids):
        native_parents = list(row.get('parents') or [])
        row['candidate_id'] = mapping[native_id]
        row['parents'] = [mapping.get(v, v) for v in row.get('parents') or []]
        if inputs and native_id not in inputs and not (set(row['parents']) & inputs):
            raise ValueError('Worker result lost its input candidate lineage.')
        row.setdefault('raw_metadata', {})['candidate_execution'] = dict(
            batch_id=batch['batch_id'], native_candidate_id=native_id, native_parents=native_parents,
            unit_kind=manifest.get('unit_kind', 'candidate_collection'),
            input_manifest_sha256=manifest['sha256'], scientific_key=manifest['scientific_key'],
            pipeline=manifest['pipeline'], node=manifest['node'],
            effective_request=manifest['request'], attempts=copy.deepcopy(batch['attempts']),
            worker=attempt.get('worker'), allocation=attempt.get('allocation'),
            software=attempt.get('software'), hardware=attempt.get('hardware'), output_hashes=data_hashes(row))
    batch['source_provenance'] = output.get('source_provenance', {})
    batch['source_run_id'] = output.get('source_run_id')
    batch['output_hashes'] = data_hashes(rows)
    accounted = {p for row in rows for p in row.get('parents', [])} | {row['candidate_id'] for row in rows}
    batch['missing_candidate_ids'] = sorted(inputs - accounted)
    expected = manifest.get('expected_outputs_per_input')
    if expected is not None:
        membership = {key: {'expected_outputs': expected, 'returned_outputs': 0, 'valid_outputs': 0} for key in inputs}
        for row in rows:
            represented = ({row['candidate_id']} | set(row.get('parents') or [])) & inputs
            outcome = row.get('raw_metadata', {}).get('measurement_status', {})
            outcome = outcome.get('status') if isinstance(outcome, dict) else outcome
            valid = bool(row.get('complex_pdb') or row.get('binder_pdb')) and outcome not in {'error', 'failed', 'invalid', 'missing'}
            for member in represented:
                membership[member]['returned_outputs'] += 1
                membership[member]['valid_outputs'] += valid
        for values in membership.values():
            values['status'] = 'completed' if values['valid_outputs'] >= expected else 'partial' if values['valid_outputs'] else 'missing'
        batch['candidate_states'] = membership
        batch['expected_output_count'] = manifest['expected_output_count']
        batch['returned_output_count'] = len(rows)
        batch['missing_output_count'] = sum(max(0, v['expected_outputs'] - v['valid_outputs']) for v in membership.values())
    batch['returned_input_units'] = len(inputs & accounted) if inputs else batch['unit_count'] if rows else 0
    batch['completed_units'] = (sum(value['status'] == 'completed' for value in batch['candidate_states'].values())
                                if expected is not None else batch['returned_input_units'])
    if manifest.get('unit_kind') == 'native_invocation':
        batch['completed_units'] = 1  # Completion of one job, never inferred scientific acceptance.
    batch['candidate_count'] = len(rows)
    batch['status'] = 'partial' if batch['missing_candidate_ids'] or batch.get('missing_output_count') else 'completed'
    batch['completed_at'] = utc_now()
    return rows


def aggregate(batches, *, stream=False):
    """Stable input-batch ordering, refusing rather than dropping duplicate results."""
    def rows():
        from mn_molecular_design_studio.candidate_sets import index
        with index() as db:
            db.execute('CREATE TABLE seen (id TEXT PRIMARY KEY)')
            for batch in sorted(batches, key=lambda b: b['offset']):
                if batch['status'] not in {'completed', 'partial'}:
                    continue
                if not output_valid(batch):
                    raise ValueError('Completed candidate artifacts failed their checksum verification.')
                values = json.loads(Path(batch['rows_path']).read_text())
                for row in values:
                    if db.execute('SELECT 1 FROM seen WHERE id=?', (row['candidate_id'],)).fetchone():
                        raise ValueError('Duplicate candidate identity during batch aggregation.')
                    db.execute('INSERT INTO seen VALUES (?)', (row['candidate_id'],))
                    yield row
    return rows() if stream else list(rows())


def progress(batches):
    return dict(total_units=sum(b['unit_count'] for b in batches),
                total_batches=len(batches),
                stopped_units=sum(b['unit_count'] for b in batches if b.get('objective_stop')),
                completed_units=sum(b.get('completed_units', 0) for b in batches if b['status'] in {'completed', 'partial'}),
                completed_batches=sum(b['status'] == 'completed' for b in batches),
                failed_batches=sum(b['status'] in {'failed', 'lost'} for b in batches),
                retry_count=sum(max(0, len(b['attempts']) - 1) for b in batches))


def execution_status(metadata):
    """Fence controller retries against an existing scheduler's own attempts."""
    from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import _pid_alive
    status = metadata.get('status') or 'lost'
    pid = metadata.get('worker_pid')
    if status == 'running' and pid and not _pid_alive(pid):
        status = 'lost'
    if metadata.get('scheduler_status') in {'failed', 'lost'} and status not in TERMINAL:
        status = 'lost'
    if metadata.get('scheduler_job_id') and status in {'failed', 'lost'}:
        from mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler import _client
        try:
            scheduled = _client().get_job(metadata['scheduler_job_id'])
        except Exception:
            return 'reconciling'
        if scheduled.get('status') not in {'failed', 'lost', 'cancelled', 'completed'}:
            return 'scheduler_retrying'
        if pid and _pid_alive(pid):
            return 'terminating'
    return status


def tick(parent, node, current, source, params, key, settings, *, graph, owner_id, project_id, gpu_device, persist,
         shared_campaign_ledger=None, generator_node=None):
    """Advance bounded ordinary scheduler jobs; called by the graph's existing loop."""
    import socket
    import time
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import submission_tags
    allow_development = bool(read_json(Path(parent) / 'metadata.json').get('studio_development_mode'))
    if 'batches' not in current:
        current['batches'] = initialize(parent, node, source, params, key, settings,
            make_source=graph._local_source, make_request=graph.node_request, read_rows=graph._absolute_candidates)
        current.update(status='running', fingerprint=key)
        current.setdefault('started_at', utc_now())
        persist()
    batches = current['batches']
    limits = dispatch_limits(node['data']['kind'], settings)
    for batch in batches:
        if (shared_campaign_ledger and batch['status'] in {'pending', 'retrying'}
                and not batch.get('job_run_dir') and not batch.get('cluster_job_id')
                and shared_campaign_ledger.snapshot().get('status') != 'running'):
            batch.update(status='cancelled', objective_stop=True,
                         stop_reason=shared_campaign_ledger.snapshot().get('termination_reason'))
            continue
        # The controller can stop after persisting the reservation but before
        # recording the child path. Re-enter the existing stable-submission-key
        # reconciliation path; it will attach a created child or submit once.
        if batch.get('status') == 'submitting' and not batch.get('job_run_dir') and not batch.get('cluster_job_id'):
            batch['status'] = 'retrying'
            batch['recovery_reason'] = 'Reconciled an interrupted scheduler submission reservation.'
        if batch['status'] in {'completed', 'partial'}:
            if shared_campaign_ledger:
                attempt = (batch.get('attempts') or [{}])[-1]
                shared_id = batch.pop('shared_job_id', None) or attempt.get('shared_job_id')
                if shared_id:
                    shared_campaign_ledger.finish_job(shared_id, batch['status'])
            continue
        if batch.get('quarantined'):
            continue
        status = None
        if batch.get('cluster_job_id'):
            from mn_molecular_design_studio.cluster import candidate_batches
            status, output, remote_job = candidate_batches.poll(batch, owner_id=owner_id, project_id=project_id,
                                                              allow_development=allow_development)
            batch['status'] = status
            if status == 'disconnected':
                batch.setdefault('disconnected_since', time.time())
                if time.time() - batch['disconnected_since'] >= settings['stall_timeout_seconds']:
                    # A disconnected job may still be computing. Queue cancellation
                    # and quarantine its identity; never sample it again elsewhere
                    # until the original execution can be reconciled.
                    candidate_batches.cancel(batch, owner_id=owner_id, project_id=project_id,
                                             allow_development=allow_development)
                    batch.update(status='lost', quarantined=True,
                        error='Worker disconnected beyond the recovery window; cancellation queued. No duplicate attempt was launched.')
            else:
                batch.pop('disconnected_since', None)
            attempt = batch['attempts'][-1]
            attempt.update(status=status, worker=remote_job.get('node_id'),
                           allocation=remote_job.get('metadata', {}).get('resource_allocation'))
            if remote_job.get('native_execution'):
                attempt['native_execution'] = remote_job['native_execution']
            if output is not None:
                attempt['software'] = remote_job.get('metadata', {}).get('software_environment', {})
                rows = collect(batch, output, read_rows=graph._absolute_candidates)
                if status == 'partial':
                    batch['status'] = 'partial'
                path = Path(batch['manifest']).with_suffix('.results.json')
                atomic_json(path, rows)
                batch['rows_path'] = str(path)
                batch['output_hashes'][str(path)] = file_digest(path)
            elif status in {'failed', 'lost'}:
                attempt.update(error=remote_job.get('wait_reason') or 'Remote candidate job failed.', finished_at=utc_now())
                if not batch.get('quarantined') and len(batch['attempts']) <= settings['retries']:
                    batch.pop('cluster_job_id', None)
                    batch.pop('prepared_envelope', None)
                    batch.pop('transfers', None)
                    batch['status'] = 'retrying'
            if status in {'completed', 'partial', 'failed', 'lost', 'cancelled'} and shared_campaign_ledger:
                attempt = (batch.get('attempts') or [{}])[-1]
                shared_id = batch.pop('shared_job_id', None) or attempt.get('shared_job_id')
                if shared_id:
                    shared_campaign_ledger.finish_job(shared_id, status)
        if batch.get('job_run_dir'):
            child = Path(batch['job_run_dir'])
            metadata = read_json(child / 'metadata.json')
            batch['wait_reason'] = metadata.get('scheduler_wait_reason') or ''
            status = execution_status(metadata)
            if status in {'reconciling', 'scheduler_retrying', 'terminating'}:
                batch['wait_reason'] = 'Reconciling the existing scheduler attempt before any controller retry.'
            attempt = batch['attempts'][-1]
            attempt.update(status=status, worker=metadata.get('execution_node_id') or socket.gethostname(),
                           allocation=metadata.get('resource_allocation'),
                           scheduler_job_id=metadata.get('scheduler_job_id'))
            receipt = read_json(child / 'native_execution.json')
            if receipt:
                attempt['native_execution'] = receipt
            batch['status'] = status
            if status == 'completed':
                result = read_json(child / 'result.json')
                if result.get('success') is not True:
                    status = batch['status'] = 'failed'
                else:
                    from mn_molecular_design_studio.cluster.capabilities import gpu_details
                    attempt['hardware'] = gpu_details()
                    attempt['worker_id'] = metadata.get('scheduler_worker_id') or f'{socket.gethostname()}:{metadata.get("worker_pid", "")}'
                    attempt['software'] = current.get('execution_specification', {}).get('runtime_identity') or read_json(child / 'input.json').get('params', {})
                    output = graph._child_output(child, source.get('target', {}))
                    rows = collect(batch, output, read_rows=graph._absolute_candidates)
                    path = Path(batch['manifest']).with_suffix('.results.json')
                    atomic_json(path, rows)
                    batch['rows_path'] = str(path)
                    batch['output_hashes'][str(path)] = file_digest(path)
            if status in {'failed', 'lost', 'cancelled'}:
                attempt.update(status=status, finished_at=utc_now(),
                               error=metadata.get('error') or metadata.get('scheduler_error') or 'Worker did not complete successfully.')
                # Retain completed native outputs from an interrupted/failed
                # collection job. Never sample completed members again merely
                # because another member failed. Missing members remain explicit.
                result_rows = child / 'artifacts/normalized_candidates/candidates.jsonl'
                if not result_rows.is_file() and receipt:
                    from mn_molecular_design_studio.ligandmpnn_batching import recover
                    recover(child)
                    if not result_rows.is_file():
                        from mn_molecular_design_studio.engine.workflow_engine.workflows.inverse_folding import recover as recover_inverse
                        recover_inverse(child)
                recoverable = result_rows.is_file() and result_rows.stat().st_size and any(
                    (r.get('complex_pdb') or r.get('binder_pdb')) and Path(r.get('complex_pdb') or r.get('binder_pdb')).is_file()
                    for r in graph._absolute_candidates(dict(run_dir=str(child), candidates_jsonl=str(result_rows))))
                if recoverable:
                    output = graph._child_output(child, source.get('target', {}))
                    rows = collect(batch, output, read_rows=graph._absolute_candidates)
                    path = Path(batch['manifest']).with_suffix('.results.json')
                    atomic_json(path, rows)
                    batch['rows_path'] = str(path)
                    batch['output_hashes'][str(path)] = file_digest(path)
                    batch['status'] = 'partial'
                    batch['error'] = attempt['error']
                elif status == 'cancelled':
                    batch['error'] = 'Candidate batch was cancelled; its native artifacts are retained.'
                elif (receipt and receipt.get('launch_count')) or read_json(child / 'command.json').get('mode') == 'docker':
                    batch['quarantined'] = True
                    batch['error'] = 'Interrupted native execution retained for reconciliation; no completed member was resampled.'
                elif len(batch['attempts']) <= settings['retries']:
                    batch.pop('job_run_dir', None)
                    batch['status'] = 'retrying'
                else:
                    batch['error'] = attempt['error']
            elif status == 'cancelled':
                batch['error'] = 'Candidate batch was cancelled.'
        terminal_batch_status = batch['status'] if batch['status'] in {'completed', 'partial', 'failed', 'lost', 'cancelled'} else (
            status if status in {'completed', 'partial', 'failed', 'lost', 'cancelled'} else None)
        if terminal_batch_status and shared_campaign_ledger:
            attempt = (batch.get('attempts') or [{}])[-1]
            shared_id = batch.pop('shared_job_id', None) or attempt.get('shared_job_id')
            if shared_id:
                shared_campaign_ledger.finish_job(shared_id, terminal_batch_status)
        if batch['status'] in {'pending', 'retrying'}:
            from mn_molecular_design_studio.campaign_objectives import time_progress
            if limits and time_progress(limits, current['started_at'])['time_expired']:
                # Recover an existing submission read-only. Repeating a submit
                # here could start another calculation after the time limit.
                pending_key = batch.get('submission_key')
                recovered = None
                matches = []
                if pending_key:
                    from mn_molecular_design_studio.runtime import runs_root
                    matches = runs_with_submission_key(runs_root(), pending_key)
                    if len(matches) > 1:
                        raise ValueError('Duplicate candidate submission receipt.')
                    if settings['placement'] != 'local':
                        from mn_molecular_design_studio.cluster import candidate_batches, config
                        if config.enabled():
                            recovered = candidate_batches.recover_submission(
                                pending_key, owner_id=owner_id, project_id=project_id,
                                allow_development=allow_development)
                if recovered or matches:
                    attempt = dict(index=len(batch['attempts']), status='queued',
                                   submission_key=pending_key, recovered_after_deadline=True)
                    if recovered:
                        attempt['cluster_job_id'] = recovered['cluster_job_id']
                        batch['cluster_job_id'] = recovered['cluster_job_id']
                    else:
                        attempt['run_dir'] = str(matches[0])
                        batch['job_run_dir'] = str(matches[0])
                    batch['attempts'].append(attempt)
                    batch['status'] = 'queued'
                else:
                    batch.update(status='cancelled', objective_stop=True, stop_reason='runtime_exhausted')
                    shared_id = batch.pop('shared_job_id', None)
                    if shared_campaign_ledger and shared_id:
                        shared_campaign_ledger.finish_job(shared_id, 'cancelled')
                persist()
                continue
            active = sum(b['status'] not in TERMINAL | {'pending', 'retrying'} for b in batches)
            if active >= settings['max_in_flight']:
                continue
            manifest = read_manifest(batch)
            attempt_index = len(batch['attempts'])
            submission_key = digest([str(parent), node['id'], key, batch['batch_id'], attempt_index])
            if shared_campaign_ledger and not shared_campaign_ledger.reserve_job(
                    submission_key, str(generator_node or node['id'])):
                batch['wait_reason'] = 'Shared campaign generation concurrency is full.'
                continue
            if shared_campaign_ledger:
                batch['shared_job_id'] = submission_key
            batch['submission_key'] = submission_key
            # A durable creation tag closes the submission-response crash window.
            from mn_molecular_design_studio.runtime import runs_root
            matches = runs_with_submission_key(runs_root(), submission_key)
            if len(matches) > 1:
                raise ValueError('Duplicate candidate submission receipt.')
            request = copy.deepcopy(manifest['request'])
            request['resources'] = dict(cpu_cores=settings['cpu_cores'], min_vram_gb=settings['min_vram_gb'], ram_gb=settings['ram_gb'])
            selected_gpu = settings.get('gpu_device', gpu_device)
            from mn_molecular_design_studio.quality_tools import WORKFLOWS as quality_workflows
            if request['workflow'].startswith('analysis.') and not quality_workflows.get(request['workflow'], {}).get('gpu', False):
                selected_gpu = 'none'
                request['resources']['min_vram_gb'] = 0
            if str(selected_gpu) == 'all':
                selected_gpu = 'auto'
                if request['workflow'] == 'refolding.optimized_prediction':
                    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
                    p = request['parameters']
                    if int(ENGINE_PROFILES[p['engine']].get('gpu_count', {}).get(p.get('optimization_mode', 'off'), 1)) > 1:
                        selected_gpu = 'all'
            batch['status'] = 'submitting'
            persist()
            started = time.monotonic()
            if settings['placement'] != 'local':
                from mn_molecular_design_studio.cluster import candidate_batches, config
                if config.enabled():
                    request['resources']['gpu_device'] = selected_gpu
                    dispatch_settings = {**settings, 'placement': manifest.get('source_node_id') or settings['placement']}
                    result = candidate_batches.submit(request, batch['manifest'], submission_key, dispatch_settings,
                        owner_id=owner_id, project_id=project_id, state=batch,
                        allow_development=allow_development)
                    batch['attempts'].append(dict(index=attempt_index, cluster_job_id=result['cluster_job_id'],
                        status='queued', started_at=utc_now(), submission_key=submission_key,
                        **({'shared_job_id': submission_key} if shared_campaign_ledger else {})))
                    batch.update(status='queued', cluster_job_id=result['cluster_job_id'])
                    persist()
                    continue
                if settings['placement'] != 'auto':
                    raise ValueError('Configure trusted peers before selecting a remote workstation.')
            if matches:
                child = matches[0]
                if not (child / 'worker_request.json').is_file():
                    raise ValueError('Interrupted job preparation; execution was not resubmitted to avoid duplicate sampling.')
            else:
                with submission_tags(candidate_submission_key=submission_key,
                                     pipeline_parent_run_dir=str(parent), graph_node_id=node['id'],
                                     candidate_batch_id=batch['batch_id'], input_manifest_sha256=manifest['sha256']):
                    child = graph._submit(request, parent, {**node, 'id': f'{node["id"]}_{batch["batch_id"]}_a{attempt_index}', 'candidate_batch': True, 'batch_offset': batch['offset']},
                                          owner_id=owner_id, project_id=project_id, gpu_device=selected_gpu)
            batch['attempts'].append(dict(index=attempt_index, run_dir=str(child), status='queued',
                                          started_at=utc_now(), submission_key=submission_key,
                                          scheduler_overhead_seconds=time.monotonic() - started,
                                          **({'shared_job_id': submission_key} if shared_campaign_ledger else {})))
            batch.update(status='queued', job_run_dir=str(child))
            persist()
    current['candidate_progress'] = progress(batches)
    from mn_molecular_design_studio.native_execution import execution_totals
    current['native_execution'] = execution_totals(batches)
    current['candidate_count'] = sum(b.get('candidate_count', 0) for b in batches if b['status'] in {'completed', 'partial'})
    persist()
    if any(b['status'] not in TERMINAL for b in batches):
        return None
    failures = [b['batch_id'] for b in batches if b['status'] != 'completed' and not b.get('objective_stop')]
    if source.get('candidate_set_manifest'):
        from mn_molecular_design_studio.candidate_sets import write_collection
        output = write_collection(parent, node['id'] + '_batches', aggregate(batches, stream=True), source['target'])
        output['collection_stage_pending'] = True
    else:
        output = graph._local_source(parent, node['id'] + '_batches', aggregate(batches), source.get('target'))
    if not source.get('candidate_set_manifest') and graph.NODE_CATALOG[node['data']['kind']]['category'] == 'Structure prediction':
        output = graph._attach_refold_input_references(output, source, parent, node['id'])
    output['candidate_execution'] = dict(policy=settings, progress=progress(batches), missing_batches=failures,
                                       native_execution=execution_totals(batches))
    if limits:
        from mn_molecular_design_studio.campaign_objectives import time_progress
        output['candidate_execution']['time_limit'] = time_progress(limits, current['started_at'])
        output['candidate_execution']['started_at'] = current['started_at']
    if node.get('native_invocation') and len(batches) == 1:
        output.update(source_provenance=batches[0].get('source_provenance', {}),
                      source_run_id=batches[0].get('source_run_id'))
    current['partial'] = bool(failures)
    if failures:
        current['error'] = 'Unsuccessful candidate batches: ' + ', '.join(failures)
    return output, ('partial' if output['candidate_count'] and settings['failure_policy'] == 'partial' else 'failed') if failures else 'completed'


def cancel(parent, state):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
    metadata = read_json(Path(parent) / 'metadata.json')
    for node in state.get('nodes', {}).values():
        for group in node.get('incremental_groups', {}).values():
            cancel(parent, {'nodes': {'group': group['state']}})
        for batch in node.get('batches', []):
            if batch['status'] in {'pending', 'retrying'}:
                batch['status'] = 'cancelled'
            if batch.get('cluster_job_id') and batch['status'] not in TERMINAL:
                from mn_molecular_design_studio.cluster import candidate_batches
                marker = read_json(Path(parent) / 'molecular_design_studio.json')
                candidate_batches.cancel(batch, owner_id=metadata.get('studio_owner_id') or marker.get('studio_owner_id'),
                    project_id=metadata.get('studio_project_id') or marker.get('project_id'),
                    allow_development=bool(metadata.get('studio_development_mode')))
                batch['status'] = 'cancelled'
            if batch.get('job_run_dir') and batch['status'] not in TERMINAL:
                stop_job(Path(batch['job_run_dir']), reason='Parent pipeline was cancelled')
                batch['status'] = 'cancelled'


def incremental_cpu_kind(kind):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.quality_tools import WORKFLOWS
    workflow = NODE_CATALOG[kind].get('workflow') or ''
    return (workflow.startswith('analysis.') and not WORKFLOWS.get(workflow, {}).get('gpu', False)
            and NODE_CATALOG[kind]['outputs'][0]['type'] != 'analysis')


def tick_incremental(parent, node, current, upstream_state, upstream_output, settings, *, graph, owner_id, project_id, persist):
    """Overlap allowlisted CPU analyses with completed upstream GPU batches.

    Complete-pool selection/Results and multi-parent nodes retain their barriers.
    Groups share the same controller loop and scheduler; no executor is added.
    """
    groups = current.setdefault('incremental_groups', {})
    current.setdefault('started_at', utc_now())
    kind, params = node['data']['kind'], node['data']['params']
    from mn_molecular_design_studio.cache_provenance import execution_identity
    if 'incremental_runtime' not in current:
        request = graph.node_request(kind, params, upstream_output, 'incremental-analysis')
        current['incremental_runtime'] = execution_identity(kind, request['parameters'], request['workflow'])
    for batch in upstream_state.get('batches', []):
        if batch['status'] not in {'completed', 'partial'}:
            continue
        if batch['batch_id'] not in groups:
            rows = json.loads(Path(batch['rows_path']).read_text())
            source = graph._local_source(parent, f'{node["id"]}_{batch["batch_id"]}_stream', rows, upstream_output.get('target'))
            groups[batch['batch_id']] = dict(source=source, state={'started_at': current['started_at']}, upstream_manifest=batch['manifest_sha256'])
            persist()
    for batch_id, group in groups.items():
        if group.get('output'):
            if not all(output_valid(b) for b in group['state'].get('batches', []) if b['status'] in {'completed', 'partial'}):
                raise ValueError('Incremental analysis artifacts failed checksum verification.')
            continue
        active = sum(b['status'] not in TERMINAL | {'pending', 'retrying'}
                     for g in groups.values() for b in g['state'].get('batches', []))
        state = group['state']
        has_active = any(b['status'] not in TERMINAL | {'pending', 'retrying'} for b in state.get('batches', []))
        if active >= settings['max_in_flight'] and not has_active:
            continue
        local_node = copy.deepcopy(node); local_node['id'] = f'{node["id"]}_{batch_id}'
        key = digest([kind, params, group['upstream_manifest'], current['incremental_runtime']])
        result = tick(parent, local_node, state, group['source'], params, key,
            settings | {'max_in_flight': 1, 'placement': 'local'}, graph=graph,
            owner_id=owner_id, project_id=project_id, gpu_device='none', persist=persist)
        if result:
            group['output'], group['status'] = result
            persist()
    current['status'] = 'running'
    current['candidate_progress'] = progress([b for g in groups.values() for b in g['state'].get('batches', [])])
    upstream_terminal = upstream_state.get('status') in {'completed', 'cached', 'partial', 'failed'}
    if not upstream_terminal or any(not g.get('output') for g in groups.values()):
        return None
    if upstream_state.get('status') == 'failed' and settings['failure_policy'] != 'partial':
        current.update(status='blocked', error='The complete upstream candidate collection is required.')
        persist()
        return None
    rows = [row for batch_id, g in sorted(groups.items()) for row in graph._absolute_candidates(g['output'])]
    output = graph._local_source(parent, node['id'] + '_incremental', rows, upstream_output.get('target'))
    failed = [k for k, g in groups.items() if g['status'] != 'completed']
    if upstream_state.get('partial') or upstream_state.get('status') == 'failed':
        current['partial'] = True
        failed.extend('upstream:' + b['batch_id'] for b in upstream_state.get('batches', []) if b['status'] != 'completed')
    output['candidate_execution'] = dict(policy=settings, missing_batches=failed, incremental_cpu=True,
                                        progress=current['candidate_progress'])
    key = digest([kind, params, current['incremental_runtime'], [g['upstream_manifest'] for k, g in sorted(groups.items())]])
    return output, key, ('partial' if rows and settings['failure_policy'] == 'partial' else 'failed') if failed else 'completed'


def inherit_partial_pool(output, upstream=(), *, missing_branches=()):
    """Carry incomplete-pool context through ordinary successful middle stages."""
    branches = set(missing_branches)
    batches = {}
    for source in [output, *upstream]:
        pool = source.get('partial_pool', {})
        branches.update(pool.get('missing_branches', []))
        records = list(pool.get('candidate_batches', []))
        execution = source.get('candidate_execution', {})
        if execution.get('missing_batches'):
            records.append(execution)
        for record in records:
            batches[digest(record)] = record
    return dict(missing_branches=sorted(branches), candidate_batches=[batches[k] for k in sorted(batches)]) if branches or batches else None
