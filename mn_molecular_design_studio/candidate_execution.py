"""Durable candidate batches for the existing graph controller and job scheduler.

No execution service lives here: each attempt is an ordinary Studio workflow job.
Saved graphs without an execution policy retain their native collection behavior.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path

from mn_molecular_design_studio.cache_provenance import file_digest, referenced_files
from mn_molecular_design_studio.cluster.config import atomic_json
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, utc_now

# Only adapters whose unit of work is an independent candidate are enabled.
COLLECTION_KINDS = {
    'ligandmpnn', 'foundry_mpnn', 'abmpnn', 'af2', 'af3', 'rf3', 'boltz2',
    'protenix', 'ipsae', 'rosetta_metrics',
}
GENERATION_KINDS = {'rfdiffusion', 'rfd3', 'genie3', 'pxdesign', 'rfd4'}
TERMINAL = {'completed', 'failed', 'cancelled', 'lost'}


def policy(value=None):
    value = value or {}
    if not isinstance(value, dict):
        raise ValueError('Candidate execution must be an object.')
    defaults = dict(mode='native', placement='local', batch_size=4, max_in_flight=2,
                    retries=1, failure_policy='all', min_vram_gb=8., ram_gb=0.,
                    cpu_cores=4, stall_timeout_seconds=1800)
    if set(value) - defaults.keys():
        raise ValueError('Unsupported candidate execution settings: ' + ', '.join(sorted(set(value) - defaults.keys())))
    result = defaults | value
    if result['mode'] not in {'native', 'sharded'} or result['failure_policy'] not in {'all', 'partial'}:
        raise ValueError('Choose native/sharded execution and all/partial failure policy.')
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
    from mn_molecular_design_studio.pipeline_graph import OPTIMIZED_PREDICTION_NODE_ENGINES
    from mn_molecular_design_studio.quality_tools import WORKFLOWS
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    workflow = NODE_CATALOG[kind].get('workflow')
    return ((bool(source.get('candidates_jsonl')) and
             (kind in COLLECTION_KINDS or kind in OPTIMIZED_PREDICTION_NODE_ENGINES or workflow in WORKFLOWS))
            or kind in GENERATION_KINDS)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def initialize(parent, node, source, params, key, settings, *, make_source, make_request, read_rows):
    """Freeze ordered units and complete effective requests before submitting any work."""
    rows = read_rows(source) if source.get('candidates_jsonl') else []
    ids = [str(r.get('candidate_id') or '') for r in rows]
    if rows and (not all(ids) or len(ids) != len(set(ids))):
        raise ValueError('Candidate sharding requires nonempty unique input candidate IDs.')
    generation = node['data']['kind'] in GENERATION_KINDS and not rows
    total = int(params.get('num_designs', 1)) if generation else len(rows)
    if total <= 0:
        return []
    directory = Path(parent) / 'candidate_batches' / node['id'] / key
    directory.mkdir(parents=True, exist_ok=True)
    batches = []
    seed = int(params.get('seed') or 0)
    if not seed:
        # Persisted scientific key -> stable positive seed on replay/retry.
        seed = int(key[:8], 16) % 2147483646 + 1
    for offset in range(0, total, settings['batch_size']):
        count = min(settings['batch_size'], total - offset)
        batch_id = f'b{offset:06d}'
        selected = rows[offset:offset + count]
        batch_source = make_source(parent, f'{node["id"]}_{batch_id}_inputs', selected, source.get('target')) if rows else source
        effective = copy.deepcopy(params)
        if generation:
            effective['num_designs'] = count
            # Never repeat the same trajectory seed on separate shards.
            effective['seed'] = (seed + offset) % 2147483646 or 1
            if 'batch_size' in effective:
                effective['batch_size'] = min(int(effective['batch_size']), count)
        else:
            effective['seed'] = seed if 'seed' in params else params.get('seed')
            if effective['seed'] is None:
                effective.pop('seed')
        request = make_request(node['data']['kind'], effective, batch_source, f'{node["data"]["label"]} · {batch_id}')
        manifest = dict(schema_version=1, pipeline=str(parent), node=node['id'], batch_id=batch_id,
                        scientific_key=key, offset=offset, input_candidate_ids=ids[offset:offset + count],
                        unit_count=count, request=request, input_hashes=referenced_files(request), policy=settings)
        manifest['sha256'] = digest(manifest)
        path = directory / f'{batch_id}.json'
        if path.exists():
            if json.loads(path.read_text()) != manifest:
                raise ValueError('Immutable candidate batch inputs changed; submit a new pipeline.')
        else:
            atomic_json(path, manifest)
        batches.append(dict(batch_id=batch_id, offset=offset, unit_count=count, manifest=str(path),
                            manifest_sha256=manifest['sha256'], status='pending', attempts=[]))
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
    rows = read_rows(output)
    native_ids = [str(r.get('candidate_id') or '') for r in rows]
    if not all(native_ids) or len(set(native_ids)) != len(native_ids):
        raise ValueError('Worker returned missing or duplicated candidate identities.')
    inputs = set(manifest['input_candidate_ids'])
    # Native IDs may include timestamps. The immutable batch and native output
    # order define identity; retain every original identity in provenance.
    mapping = {v: (v if v in inputs else f'{manifest["node"]}:{batch["batch_id"]}:{i:06d}')
               for i, v in enumerate(native_ids)}
    attempt = batch['attempts'][-1]
    for row, native_id in zip(rows, native_ids):
        row['candidate_id'] = mapping[native_id]
        row['parents'] = [mapping.get(v, v) for v in row.get('parents') or []]
        if inputs and native_id not in inputs and not (set(row['parents']) & inputs):
            raise ValueError('Worker result lost its input candidate lineage.')
        row.setdefault('raw_metadata', {})['candidate_execution'] = dict(
            batch_id=batch['batch_id'], native_candidate_id=native_id,
            input_manifest_sha256=manifest['sha256'], scientific_key=manifest['scientific_key'],
            pipeline=manifest['pipeline'], node=manifest['node'],
            effective_request=manifest['request'], attempts=copy.deepcopy(batch['attempts']),
            worker=attempt.get('worker'), allocation=attempt.get('allocation'),
            software=attempt.get('software'), output_hashes=referenced_files(row))
    batch['output_hashes'] = referenced_files(rows)
    batch['completed_units'] = batch['unit_count']
    batch['candidate_count'] = len(rows)
    batch['status'] = 'completed'
    batch['completed_at'] = utc_now()
    return rows


def aggregate(batches):
    """Stable input-batch ordering, refusing rather than dropping duplicate results."""
    rows, seen = [], set()
    for batch in sorted(batches, key=lambda b: b['offset']):
        if batch['status'] != 'completed':
            continue
        if not output_valid(batch):
            raise ValueError('Completed candidate artifacts failed their checksum verification.')
        values = json.loads(Path(batch['rows_path']).read_text())
        for row in values:
            if row['candidate_id'] in seen:
                raise ValueError('Duplicate candidate identity during batch aggregation.')
            seen.add(row['candidate_id']); rows.append(row)
    return rows


def progress(batches):
    return dict(total_units=sum(b['unit_count'] for b in batches),
                completed_units=sum(b.get('completed_units', 0) for b in batches if b['status'] == 'completed'),
                completed_batches=sum(b['status'] == 'completed' for b in batches),
                failed_batches=sum(b['status'] in {'failed', 'lost'} for b in batches),
                retry_count=sum(max(0, len(b['attempts']) - 1) for b in batches))


def tick(parent, node, current, source, params, key, settings, *, graph, owner_id, project_id, gpu_device, persist):
    """Advance bounded ordinary scheduler jobs; called by the graph's existing loop."""
    import socket
    import time
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import submission_tags
    if 'batches' not in current:
        current['batches'] = initialize(parent, node, source, params, key, settings,
            make_source=graph._local_source, make_request=graph.node_request, read_rows=graph._absolute_candidates)
        current.update(status='running', fingerprint=key, started_at=utc_now())
        persist()
    batches = current['batches']
    for batch in batches:
        if batch['status'] == 'completed':
            continue
        if batch.get('job_run_dir'):
            child = Path(batch['job_run_dir'])
            metadata = read_json(child / 'metadata.json')
            status = metadata.get('status') or 'lost'
            if metadata.get('scheduler_status') in {'failed', 'lost'} and status not in TERMINAL:
                status = 'lost'
            attempt = batch['attempts'][-1]
            attempt.update(status=status, worker=metadata.get('execution_node_id') or socket.gethostname(),
                           allocation=metadata.get('resource_allocation'),
                           scheduler_job_id=metadata.get('scheduler_job_id'))
            batch['status'] = status
            if status == 'completed':
                result = read_json(child / 'result.json')
                if result.get('success') is not True:
                    status = batch['status'] = 'failed'
                else:
                    attempt['software'] = read_json(child / 'input.json').get('params', {})
                    output = graph._child_output(child, source.get('target', {}))
                    rows = collect(batch, output, read_rows=graph._absolute_candidates)
                    path = Path(batch['manifest']).with_suffix('.results.json')
                    atomic_json(path, rows)
                    batch['rows_path'] = str(path)
                    batch['output_hashes'][str(path)] = file_digest(path)
            if status in {'failed', 'lost'}:
                attempt.update(status=status, finished_at=utc_now(),
                               error=metadata.get('error') or metadata.get('scheduler_error') or 'Worker did not complete successfully.')
                if len(batch['attempts']) <= settings['retries']:
                    batch.pop('job_run_dir', None)
                    batch['status'] = 'retrying'
                else:
                    batch['error'] = attempt['error']
            elif status == 'cancelled':
                batch['error'] = 'Candidate batch was cancelled.'
        if batch['status'] in {'pending', 'retrying'}:
            active = sum(b['status'] not in TERMINAL | {'pending', 'retrying'} for b in batches)
            if active >= settings['max_in_flight']:
                continue
            manifest = read_manifest(batch)
            attempt_index = len(batch['attempts'])
            submission_key = digest([str(parent), node['id'], key, batch['batch_id'], attempt_index])
            # A durable creation tag closes the submission-response crash window.
            from mn_molecular_design_studio.runtime import runs_root
            matches = [p.parent for p in runs_root().glob('*/*/metadata.json')
                       if read_json(p).get('candidate_submission_key') == submission_key]
            if len(matches) > 1:
                raise ValueError('Duplicate candidate submission receipt.')
            request = copy.deepcopy(manifest['request'])
            request['resources'] = dict(cpu_cores=settings['cpu_cores'], min_vram_gb=settings['min_vram_gb'], ram_gb=settings['ram_gb'])
            selected_gpu = gpu_device
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
            if matches:
                child = matches[0]
                if not (child / 'worker_request.json').is_file():
                    raise ValueError('Interrupted job preparation; execution was not resubmitted to avoid duplicate sampling.')
            else:
                with submission_tags(candidate_submission_key=submission_key,
                                     pipeline_parent_run_dir=str(parent), graph_node_id=node['id'],
                                     candidate_batch_id=batch['batch_id'], input_manifest_sha256=manifest['sha256']):
                    child = graph._submit(request, parent, {'id': f'{node["id"]}_{batch["batch_id"]}_a{attempt_index}'},
                                          owner_id=owner_id, project_id=project_id, gpu_device=selected_gpu)
            batch['attempts'].append(dict(index=attempt_index, run_dir=str(child), status='queued',
                                          started_at=utc_now(), submission_key=submission_key,
                                          scheduler_overhead_seconds=time.monotonic() - started))
            batch.update(status='queued', job_run_dir=str(child))
            persist()
    current['candidate_progress'] = progress(batches)
    current['candidate_count'] = sum(b.get('candidate_count', 0) for b in batches if b['status'] == 'completed')
    persist()
    if any(b['status'] not in TERMINAL for b in batches):
        return None
    failures = [b['batch_id'] for b in batches if b['status'] != 'completed']
    rows = aggregate(batches)
    output = graph._local_source(parent, node['id'] + '_batches', rows, source.get('target'))
    if graph.NODE_CATALOG[node['data']['kind']]['category'] == 'Structure prediction':
        output = graph._attach_refold_input_references(output, source, parent, node['id'])
    output['candidate_execution'] = dict(policy=settings, progress=progress(batches), missing_batches=failures)
    current['partial'] = bool(failures)
    if failures:
        current['error'] = 'Unsuccessful candidate batches: ' + ', '.join(failures)
    return output, ('partial' if rows and settings['failure_policy'] == 'partial' else 'failed') if failures else 'completed'


def cancel(parent, state):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
    for node in state.get('nodes', {}).values():
        for batch in node.get('batches', []):
            if batch.get('job_run_dir') and batch['status'] not in TERMINAL:
                stop_job(Path(batch['job_run_dir']), reason='Parent pipeline was cancelled')
                batch['status'] = 'cancelled'
