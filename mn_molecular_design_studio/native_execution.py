"""Read-only ordinary-workflow execution contracts, never campaign authorization.

Counts describe model entrypoint calls, not GPU forward batches. Unknown values
stay unknown; a container-level batch API is not evidence of model reuse.
"""
from __future__ import annotations

from pathlib import Path
import json


def prediction_engine(kind):
    from mn_molecular_design_studio.pipeline_graph import OPTIMIZED_PREDICTION_NODE_ENGINES
    return OPTIMIZED_PREDICTION_NODE_ENGINES.get(kind) or {'boltz2': 'boltz2', 'af3': 'alphafold3'}.get(kind)


def record_process(run_dir, command, *, native_invocations, candidate_ids=(), gpu_slots=1, runner, **kwargs):
    """Persist process receipts around the existing synchronous launcher.

    A crash in the launch window is explicitly indeterminate. Allocation counts
    come from the scheduler receipt, never from a Docker command's --gpus flag.
    """
    from mn_molecular_design_studio.cluster.config import atomic_json
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    path = Path(run_dir) / "native_execution.json"
    value = read_json(path) or {"schema": "studio.native-execution.v1", "processes": [], "prediction_outputs": 0}
    metadata = read_json(Path(run_dir) / 'metadata.json')
    row = dict(command=command, native_invocations=native_invocations,
               candidate_ids=list(candidate_ids), status="launch_window", requested_gpu_slots=gpu_slots,
               allocation=metadata.get('resource_allocation'), worker_started_at=metadata.get('worker_started_at'))
    value["processes"].append(row)
    def save():
        confirmed = all(p["status"] == "returned" for p in value["processes"])
        allocations = {json.dumps([p.get('worker_started_at'), p['allocation']], sort_keys=True): p['allocation']
                       for p in value['processes'] if p.get('allocation')}
        slot_counts = []
        for allocation in allocations.values():
            devices = allocation.get('gpu_devices')
            selector = str(allocation.get('gpu_device', '')).removeprefix('device=')
            slot_counts.append(len(devices) if isinstance(devices, list) else len(selector.split(',')) if selector and all(d.isdigit() for d in selector.split(',')) else 0 if selector in {'cpu', 'none'} else None)
        known = all(n is not None for n in slot_counts)
        value.update(docker_processes=len(value["processes"]) if confirmed else None,
                     native_invocations=sum(p["native_invocations"] for p in value["processes"]) if confirmed and all(p["native_invocations"] is not None for p in value["processes"]) else None,
                     gpu_allocation_events=sum(n > 0 for n in slot_counts) if known else None,
                     gpu_slots=sum(slot_counts) if known else None,
                     allocation_receipts=list(allocations.values()), launch_count=len(value["processes"]),
                     model_initialization_basis="Expected per native invocation; not a runtime measurement.")
        atomic_json(path, value)
    save()
    try:
        result = runner(command, **kwargs)
    except BaseException:
        save()
        raise
    row.update(status="returned", return_code=result.returncode)
    save()
    return result


def prediction_contract(engine, count, *, seeds=1, samples=1, mode="off"):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
    profile = ENGINE_PROFILES[engine]
    # Inspected pinned image-side adapters. Other selectable tiles remain explicit
    # unknowns until their native launcher/model lifecycle has been verified.
    known = engine in {"boltz2", "protenix_v2"}
    per_job = seeds if engine == "boltz2" else 1 if known else None
    gpus = int(profile.get("gpu_count", {}).get(mode, 1))
    return dict(native_multi_input=True if known else None,
                native_invocations=per_job if count else 0,
                docker_processes=1 if count else 0,
                gpu_allocation_events=1 if count else 0,
                gpu_slots=gpus if count else 0,
                prediction_outputs=count * samples * (seeds if engine == "boltz2" else 1),
                model_initializations=(per_job if gpus == 1 else None) if count else 0,
                count_basis="source-verified expectation" if known else "native lifecycle unverified",
                model_reuse=("One model entrypoint per seed for the entire input directory; initialization repeats per seed/chunk."
                             if engine == "boltz2" else
                             "One model process for the input collection; initialization repeats per chunk."
                             if known else
                             "One container receives the collection; internal model process count is unknown."))


def request_contract(request, count):
    """Expected execution of one scheduler child; conditional groups are bounds."""
    workflow = request["workflow"]
    p = request.get("parameters", {})
    if workflow == 'refolding.rfd4_proteina':
        replicates = int(p.get('num_designs', 1))
        invocations = (replicates + 99) // 100 if count else 0
        return dict(native_multi_input=True, native_invocations=invocations,
                    docker_processes=1 if count else 0, gpu_allocation_events=1 if count else 0,
                    gpu_slots=1 if count else 0, prediction_outputs=count * replicates,
                    model_initializations=invocations, count_basis='RFD4 runtime and model lifecycle source verified',
                    model_reuse='One RFD4 model per native replicate command for all input files; input names and lineage must pass existing staging checks.')
    if workflow in {'sequence.esm_if1', 'sequence.caliby', 'sequence.boltzgen_inverse'}:
        ensemble = workflow == 'sequence.caliby' and p.get('ensemble_mode', 'generate') == 'generate'
        repeats = count if workflow == 'sequence.boltzgen_inverse' else 1
        return dict(native_multi_input=True, native_invocations=1 if count else 0,
                    docker_processes=1 if count else 0, gpu_allocation_events=1 if count else 0,
                    gpu_slots=1 if count else 0, prediction_outputs=0,
                    model_initializations=None if ensemble else repeats if count else 0,
                    count_basis='Studio native runner source verified',
                    model_reuse=('One Python process, but BoltzGen constructs a new native task/model per input; initialization repeats within the process.'
                                 if workflow == 'sequence.boltzgen_inverse' else
                                 'Caliby main model loads once; generated ensembles can initialize additional models per input.' if ensemble else
                                 'Model loads before the input loop; one initialization per chunk.'))
    if workflow == "refolding.optimized_prediction":
        return prediction_contract(p["engine"], count, seeds=p.get("engine_options", {}).get("seed_count", 1),
                                   samples=p.get("samples", 1), mode=p.get("optimization_mode", "off"))
    if workflow in {"sequence.mpnn", "sequence.ligandmpnn", "sequence.foundry_mpnn", "sequence.abmpnn"}:
        foundry = workflow == "sequence.foundry_mpnn" or p.get("backend") in {"auto", "foundry", None}
        return dict(native_multi_input=None if foundry else True,
                    native_invocations=None, docker_processes=None,
                    native_invocations_min=1 if count else 0, native_invocations_max=count,
                    gpu_allocation_events=1 if count else 0, gpu_slots=1 if count else 0,
                    prediction_outputs=0, model_initializations=None,
                    count_basis="backend and compatible constraint groups resolved during staging",
                    model_reuse="LigandMPNN reuses one model per compatible multi-PDB group. Foundry uses one process per structure. Initial co-design and disabled redesign use zero processes.")
    return dict(native_multi_input=None, native_invocations=None, docker_processes=None,
                gpu_allocation_events=None, gpu_slots=None, prediction_outputs=None,
                model_initializations=None, count_basis="native lifecycle unverified",
                model_reuse="Current tool contract retained; internal process count is unknown.")


def tool_inventory():
    """Inventory every registered scientific tile, plus operation/model variants.

    Registration/collection compatibility does not certify scientific eligibility
    or native batching. Uninspected adapters stay explicitly unverified.
    """
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, OPTIMIZED_PREDICTION_NODE_ENGINES
    from mn_molecular_design_studio import mpnn_input, inverse_folding_input, rfd4_proteina_input
    result = []
    for kind, spec in NODE_CATALOG.items():
        if spec['category'] not in {'Design engine', 'Design pipeline', 'Sequence design', 'Structure prediction'}:
            continue
        params = {f['name']: f['default'] for f in spec['fields']}
        engine = prediction_engine(kind)
        contract = prediction_contract(engine, 5) if engine else request_contract(dict(workflow=spec['workflow'], parameters=params), 5)
        row = dict(tile=kind, category=spec['category'], workflow=spec['workflow'], label=spec['label'], hidden=bool(spec.get('hidden')),
                   collection_socket=any('candidates' in s['accepts'] for s in spec['inputs']),
                   native_multi_input=contract['native_multi_input'], count_basis=contract['count_basis'], model_reuse=contract['model_reuse'])
        if kind in {'rfd4', 'rfd4_fold', 'bindcraft2', 'rfdiffusion', 'rfd3'}:
            row.update(native_multi_input=kind in {'rfd4', 'rfd4_fold', 'rfd3'}, count_basis='native lifecycle source verified; generator fixtures and prior Phase 6C evidence',
                model_reuse='Ordinary attempts/samples reuse models within an invocation; BindCraft2 reuse is per worker, validation can rebuild. Campaign limit remains one verified attempt per child.')
        if kind in {'ligandmpnn', 'foundry_mpnn', 'abmpnn'}:
            row['variants'] = ', '.join(mpnn_input.MODELS) if kind != 'abmpnn' else 'ab_mpnn'
        if kind == 'caliby':
            row['variants'] = ', '.join(inverse_folding_input.CALIBY_MODELS) + '; single/generated/provided ensembles'
        result.append(row)
        if kind == 'rfd4':
            for operation, label in rfd4_proteina_input.OPERATIONS.items():
                result.append(dict(row, tile='rfd4/operation/' + operation, label='RFD4 ' + label,
                    category='Sequence design' if operation == 'inverse_folding' else 'Structure prediction' if operation == 'structure_prediction' else 'Design engine',
                    collection_socket=False,
                    model_reuse='Existing RFD4 existing-structure operation accepts several input files and replicates per process. Dedicated CandidateSet sequence-stage handoff is not yet verified; do not substitute structure prediction for inverse folding.'))
    return result


def execution_totals(batches):
    """Receipts across every attempt, including failures/retries; never count plans."""
    receipts = [a.get("native_execution") for b in batches for a in b.get("attempts", [])]
    result = dict(candidates=sum(b["unit_count"] for b in batches), collection_chunks=len(batches),
                  scheduler_attempts=len(receipts), receipt_count=sum(r is not None for r in receipts),
                  returned_records=sum(b.get("candidate_count", 0) for b in batches))
    for key in ("native_invocations", "docker_processes", "gpu_allocation_events", "gpu_slots", "prediction_outputs"):
        result[key] = sum(r[key] for r in receipts) if receipts and all(r is not None and r.get(key) is not None for r in receipts) else (0 if not receipts else None)
    return result
