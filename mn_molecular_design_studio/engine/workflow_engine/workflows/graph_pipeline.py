"""Execute a saved visual graph through the normal scheduler job adapters."""
from __future__ import annotations

import copy
import hashlib
import inspect
import json
import operator
import re
import shutil
import signal
import tempfile
import time
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, fingerprint, scope_root, validate_graph
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, update_status, utc_now, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _path(value, run_dir=None) -> Path:
    path = resolve_stored_path(value, run_dir=run_dir, must_exist=True)
    if path is None:
        raise ValueError(f"Pipeline input is unavailable: {value}")
    return path


def _absolute_candidates(source: dict) -> list[dict]:
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    run_dir = _path(source["run_dir"])
    rows = copy.deepcopy(read_candidates(_path(source["candidates_jsonl"], run_dir)))
    for row in rows:
        for key in ("target_pdb", "complex_pdb", "binder_pdb"):
            if row.get(key):
                resolved = resolve_stored_path(row[key], run_dir=run_dir)
                if resolved:
                    row[key] = str(resolved)
    return rows


def _local_source(parent: Path, node_id: str, candidates: list[dict], target: dict | None = None) -> dict:
    """Write a new candidate collection, keeping upstream artifacts read-only."""
    destination = parent / "artifacts" / "pipeline_nodes" / node_id
    normalized = destination / "artifacts" / "normalized_candidates"
    normalized.mkdir(parents=True, exist_ok=True)
    candidate_path = normalized / "candidates.jsonl"
    candidate_path.write_text("".join(json.dumps(c) + "\n" for c in candidates))
    write_json(destination / "metadata.json", {"status": "completed", "tool": "pipeline_collection", "run_id": node_id})
    write_json(destination / "input.json", {"inputs": {}, "params": {}})
    write_json(destination / "result.json", {"success": True, "outputs": {"candidates": candidates}})
    return {"run_dir": str(destination), "candidates_jsonl": str(candidate_path), "target": target or {}, "candidate_count": len(candidates)}


def source_designs(source_run_dir: str, project_id: str, selected_ids=None) -> tuple[Path, list[dict]]:
    from mn_molecular_design_studio.design_validation import design_candidates
    run_dir = _path(source_run_dir)
    marker = read_json(run_dir / "molecular_design_studio.json")
    if str(marker.get("project_id") or "") != project_id:
        raise ValueError("The selected design job must belong to this project.")
    candidates = design_candidates(run_dir)
    if selected_ids is not None:
        wanted = set(selected_ids)
        if not wanted.issubset({str(c.get("candidate_id")) for c in candidates}):
            raise ValueError("Only filter-passing designs from the selected successful job can be used.")
        candidates = [c for c in candidates if str(c.get("candidate_id")) in wanted]
    if not candidates:
        raise ValueError("This job has no successful, filter-passing designs.")
    return run_dir, candidates


def source_generations(source_run_dir: str, project_id: str, selected_ids=None, *, stage=None) -> tuple[Path, list[dict]]:
    from mn_molecular_design_studio.generation_io import generated_candidates
    run_dir=_path(source_run_dir)
    if str(read_json(run_dir/'molecular_design_studio.json').get('project_id') or '')!=project_id:
        raise ValueError('The selected generator job must belong to this project.')
    candidates=generated_candidates(run_dir,stage=stage)
    if selected_ids is not None:
        wanted=set(selected_ids)
        if not wanted.issubset({c['candidate_id'] for c in candidates}):
            raise ValueError('Select generated designs from the completed generator job.')
        candidates=[c for c in candidates if c['candidate_id'] in wanted]
    if not candidates:raise ValueError('This job has no completed generation outputs.')
    return run_dir,candidates


def _generation_target(params, source=None):
    target=copy.deepcopy((source or {}).get('target') or {})
    own_input=bool(params.get('target_pdb'))
    target.update({k:v for k,v in params.items() if k in {'target_pdb','target_chains','hotspots','pipeline_task','ligand','motif_atom_spec'} and (own_input or v not in ('',None,[]))})
    if params.get('ligand_residue'):
        target['ligand_original_code']=params['ligand_residue']['name']
    chains=target.get('target_chains',[])
    if isinstance(chains,str):chains=[c for c in re.split(r'[,;\s]+',chains) if c]
    target['target_chains']=chains
    return target


def filter_candidates(rows: list[dict], settings: dict) -> list[dict]:
    from mn_molecular_design_studio.design_validation import candidate_passes_filters
    metric = str(settings.get("metric") or "").strip()
    compare = {">=": operator.ge, "<=": operator.le, ">": operator.gt, "<": operator.lt}[settings["operator"]]
    def score(row):
        value = row.get("metrics", {}).get(metric)
        try:
            return float(value) if value is not None and not isinstance(value, bool) else None
        except (TypeError, ValueError):
            return None
    result = [row for row in rows if not settings["passing_only"] or candidate_passes_filters(row)]
    if metric:
        import math
        result = [row for row in result if score(row) is not None and math.isfinite(score(row)) and compare(score(row), settings["threshold"])]
        result.sort(key=score, reverse=settings["higher_is_better"])
    limit = int(settings["limit"])
    return result[:limit] if limit else result


def node_request(kind: str, params: dict, source: dict, name: str) -> dict:
    """Only allowlisted adapters can turn graph data into a workflow call."""
    workflow = NODE_CATALOG[kind]["workflow"]
    p = copy.deepcopy(params)
    from mn_molecular_design_studio import inverse_folding_input
    if kind in inverse_folding_input.ENGINES:
        if source.get("candidates_jsonl"):
            p.update(source_run_dir=source["run_dir"], candidates_jsonl=source["candidates_jsonl"])
        else:
            target = source.get("target") or {}
            p["structure_path"] = target.get("target_pdb")
            if not p.get("design_chains"):
                chains = target.get("target_chains") or []
                p["design_chains"] = chains if isinstance(chains, str) else ",".join(chains)
        p["job_name"] = name
        return {"schema_version": 1, "workflow": workflow, "parameters": p}
    from mn_molecular_design_studio import generator_input
    if kind in generator_input.ENGINES:
        inherited=source.get('target') or {}
        own=bool(p.get('target_pdb'))
        p['target_pdb']=p.get('target_pdb') or inherited.get('target_pdb') or None
        chains=p.get('target_chains') if own else inherited.get('target_chains',p.get('target_chains'))
        p['target_chains']=[v for v in re.split(r'[,;\s]+',chains) if v] if isinstance(chains,str) else list(chains or [])
        if not own and not p.get('hotspots'):p['hotspots']=inherited.get('hotspots','')
        if p.get('task') in {'unconditional','symmetry'}:p.update(target_pdb=None,target_chains=[],hotspots='')
        for key in ('scaffold_pdb','scaffold_ss_file','scaffold_adj_file','target_ss_file','target_adj_file'):
            if key in p:p[key]=p[key] or None
        p['campaign_name']=name
        return {'schema_version':1,'workflow':workflow,'parameters':p}
    if kind=='complexa_generator':
        target=_generation_target(p,source)
        p.update(target_pdb=target.get('target_pdb',''),target_chains=target['target_chains'],hotspots=target.get('hotspots',''),campaign_name=name)
        if p['pipeline_task']!='protein_binder':p.update(target_chains=[],hotspots='')
        return {"schema_version":1,"workflow":workflow,"parameters":p}
    if NODE_CATALOG[kind]["category"] in {"Design engine", "Design pipeline"}:
        target = source.get("target") or {}
        if not target.get("target_pdb"):
            raise ValueError("This design step needs a target structure.")
        p.update(target_pdb=target["target_pdb"], target_chains=target["target_chains"], hotspots=target.get("hotspots", ""), campaign_name=name)
        if kind == "rfd3":
            length = str(p.pop("binder_length"))
            match = re.fullmatch(r"(\d+)(?:-(\d+))?", length)
            if not match or int(match[1]) < 1 or int(match[2] or match[1]) < int(match[1]):
                raise ValueError("Binder length must be a positive length or range, such as 74 or 60-100.")
            p.update(task="protein_binder", target_structure=p.pop("target_pdb"), length_min=int(match[1]), length_max=int(match[2] or match[1]), job_name=p.pop("campaign_name"))
        elif kind == "rfdiffusion":
            from mn_molecular_design_studio.engine.workflow_engine.workflows.design import default_target_contig
            p["contig"] = default_target_contig(Path(p["target_pdb"]), p["target_chains"], p["binder_length"])
            if not p["contig"]:
                raise ValueError("RFdiffusion requires a PDB target with readable selected chains.")
            p["run_vanilla_pipeline"] = False
        elif kind in {"genie3", "genie3_pipeline"}:
            p["run_mode"] = "generation_only" if kind == "genie3" else "full_vanilla_pipeline"
        elif kind in {"pxdesign", "pxdesign_pipeline"}:
            p["run_mode"] = "generation_only" if kind == "pxdesign" else "pipeline"
        elif kind == "bindcraft2":
            attempts = int(p.pop("num_designs")); p.update(number_of_final_designs=attempts, max_trajectories=attempts)
        elif kind == "bindcraft2_generator":
            trajectories = int(p.pop("num_designs")); p.update(number_of_final_designs=1, max_trajectories=trajectories)
        elif kind == "boltzgen_generator":
            # The generator workflow pins the generation-only contract in the
            # workflow API adapter; this node exposes only generation controls.
            p.pop("budget", None)
        elif kind == "protpardelle_pipeline":
            p["num_mpnn_seqs"] = int(p.pop("sequences"))
        elif kind == "esm_binder":
            chains = p.pop("target_chains")
            if len(chains) != 1:
                raise ValueError("ESMFold2 Binder Design currently accepts one target chain.")
            p["target_chain"] = chains[0]
    else:
        ligand_context=source.get('target',{}).get('pipeline_task') in {'ligand_binder','motif_ligand'}
        if ligand_context and kind=='af2':
            raise ValueError('AlphaFold2 does not predict ligand complexes. Use RosettaFold3, or ESMFold for the binder alone.')
        if ligand_context and kind=='boltz2':
            tasks=[]
            for candidate in _absolute_candidates(source):
                code=(candidate.get('raw_metadata') or {}).get('ligand_original_code','')
                if not re.fullmatch(r'[A-Z0-9]{1,5}',code):
                    raise ValueError('Boltz-2 needs a CCD or SMILES ligand definition. This motif uses a custom ligand fragment; connect RosettaFold3 to preserve its input atom graph.')
                if not candidate.get('binder_sequence'):raise ValueError('The generated design has no binder sequence.')
                tasks.append({'name':candidate['candidate_id'],'binder_sequence':candidate['binder_sequence'],
                              'ligand_ccd_codes':code,'molecule_chain_options':[{'id':'A','use_msa':False}]})
            return {'schema_version':1,'workflow':'refolding.optimized_prediction','parameters':{
                'engine':'boltz2','job_name':name,'tasks':tasks,'samples':int(p['diffusion_samples']),
                'num_recycles':int(p['recycling_steps']),
                'engine_options':{'sampling_steps':int(p['sampling_steps']),'use_msa_server':False}}}
        p.update(source_run_dir=source["run_dir"], candidates_jsonl=source["candidates_jsonl"])
        if kind in {"foundry_mpnn", "ligandmpnn"}:
            from mn_molecular_design_studio import mpnn_input
            p["num_seq_per_target"] = int(p.pop("sequences")) * int(p.pop("number_of_batches", 1))
            p["noise_level"] = float(p["noise_level"]) if p.get("noise_level") not in (None, "") else None
            p = mpnn_input.prepare(p)
            from mn_molecular_design_studio.engine.workflow_engine.workflows.mpnn import run_mpnn_sequence_design
            allowed = inspect.signature(run_mpnn_sequence_design).parameters
            p = {k: v for k, v in p.items() if k in allowed}
        elif kind in {"rf3", "boltz2", "af2", "protenix"}:
            p["require_monomer_success"] = False
            if kind == "af2":
                p["use_initial_guess"] = True
            if kind == "boltz2":
                p.pop("use_target_msa", None)  # This adapter always uses shared target MSA preparation.
    return {"schema_version": 1, "workflow": workflow, "parameters": p}


def _submit(request: dict, parent: Path, node: dict, *, owner_id: str, project_id: str, gpu_device: str) -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import submit_workflow_request
    from mn_molecular_design_studio.workspace import mark_studio_job
    request = copy.deepcopy(request)
    request["resources"] = {"gpu_device": gpu_device, "cpu_cores": 4}
    path = parent / "requests" / f"{node['id']}.json"
    # Workflow submission resolves relative paths against the request file's
    # parent directory. `write_json` stores managed paths relative to the
    # enclosing run directory instead, which is a different base here because
    # requests live in a subdirectory. Keep these operational requests as JSON
    # with their absolute/portable input paths intact.
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(request, indent=2, sort_keys=True) + "\n")
    result = submit_workflow_request(path)
    mark_studio_job(result, owner_id=owner_id, project_id=project_id)
    child = _path(result["job"]["run_dir"])
    metadata = read_json(child / "metadata.json")
    metadata.update(pipeline_parent_run_dir=str(parent), campaign_id=parent.name, campaign_name=read_json(parent / "input.json").get("params", {}).get("pipeline_name", "Pipeline"), graph_node_id=node["id"])
    write_json(child / "metadata.json", metadata)
    return child


def _child_output(child: Path, target: dict) -> dict:
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    candidate_path = child / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    if not candidate_path.is_file():
        raise ValueError(f"Job {child.name} completed without a candidate collection.")
    candidates = read_candidates(candidate_path)
    params=read_json(child/'input.json').get('params',{})
    if params.get('generator_only'):
        target=_generation_target(params)
    # Bundled tools retain intermediate stage records; the next node receives
    # only final outputs, rather than redesigning/refolding all three pools.
    final = [c for c in candidates if c.get("metrics", {}).get("candidate_pool_level") == "native_pipeline"]
    if final and len(final) != len(candidates):
        output = _local_source(child, "final_outputs", _absolute_candidates({"run_dir": str(child), "candidates_jsonl": str(candidate_path)}), target)
        rows = [c for c in _absolute_candidates(output) if c.get("metrics", {}).get("candidate_pool_level") == "native_pipeline"]
        return _local_source(child, "final_outputs", rows, target)
    return {"run_dir": str(child), "candidates_jsonl": str(candidate_path), "target": target, "candidate_count": len(candidates)}


def _attach_refold_input_references(output: dict, source: dict, parent: Path, node_id: str) -> dict:
    """Keep the exact input backbone alongside each refolded result for pose RMSD."""
    source_rows = _absolute_candidates(source) if source.get("candidates_jsonl") else []
    by_id = {str(row.get("candidate_id") or ""): row for row in source_rows}
    rows = _absolute_candidates(output)
    for row in rows:
        raw = dict(row.get("raw_metadata") or {})
        direct = raw.get("source_candidate") if isinstance(raw.get("source_candidate"), dict) else {}
        possible_ids = [str(value) for value in row.get("parents") or [] if str(value)]
        if direct.get("candidate_id"):
            possible_ids.insert(0, str(direct["candidate_id"]))
        reference = next((by_id[value] for value in possible_ids if value in by_id), None)
        if reference is None and direct:
            reference = direct
        if not reference:
            continue
        structure = reference.get("complex_pdb") or reference.get("binder_pdb")
        if not structure:
            continue
        raw["pipeline_refold_input"] = {
            "complex_pdb": str(structure),
            "binder_chains": list(reference.get("binder_chains") or []),
            "target_chains": list(reference.get("target_chains") or []),
            "source_candidate_id": str(reference.get("candidate_id") or ""),
        }
        row["raw_metadata"] = raw
    return _local_source(parent, f"{node_id}_pose_inputs", rows, output.get("target"))


def _add_target_aligned_pose_rmsd(candidate: dict, parent: Path) -> dict:
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    reference = raw.get("pipeline_refold_input") if isinstance(raw.get("pipeline_refold_input"), dict) else None
    if not reference:
        return candidate
    metrics = dict(candidate.get("metrics") or {})
    reference_path = Path(str(reference.get("complex_pdb") or ""))
    prediction_path = Path(str(candidate.get("complex_pdb") or ""))
    if not reference_path.is_file() or not prediction_path.is_file():
        metrics["target_aligned_binder_rmsd"] = None
        metrics["target_aligned_binder_rmsd_status"] = "missing reference or Boltz-2 complex"
        candidate["metrics"] = metrics
        return candidate
    try:
        from mn_molecular_design_studio.engine.workflow_engine.workflows import mpnn as mpnn_workflow
        from mn_molecular_design_studio.engine.workflow_engine.workflows.benchmark import _target_aligned_binder_pose_metrics

        with tempfile.TemporaryDirectory(prefix="pipeline-pose-rmsd-") as temporary:
            temp = Path(temporary)
            reference_pdb, prediction_pdb = temp / "reference.pdb", temp / "prediction.pdb"
            reference_map = mpnn_workflow._pdb_copy(reference_path, reference_pdb)
            prediction_map = mpnn_workflow._pdb_copy(prediction_path, prediction_pdb)
            ref_binder = [reference_map[chain] for chain in reference.get("binder_chains") or [] if chain in reference_map]
            ref_target = [reference_map[chain] for chain in reference.get("target_chains") or [] if chain in reference_map]
            pred_binder = [prediction_map[chain] for chain in candidate.get("binder_chains") or [] if chain in prediction_map]
            pred_target = [prediction_map[chain] for chain in candidate.get("target_chains") or [] if chain in prediction_map]
            result = _target_aligned_binder_pose_metrics(
                reference_pdb, prediction_pdb,
                reference_binder_chains=ref_binder,
                reference_target_chains=ref_target,
                prediction_binder_chains=pred_binder or ref_binder,
                prediction_target_chains=pred_target or ref_target,
            )
        metrics.update(result)
        metrics["target_aligned_binder_rmsd_status"] = result.get("pose_rmsd_status", "")
    except Exception as exc:
        metrics["target_aligned_binder_rmsd"] = None
        metrics["target_aligned_binder_rmsd_status"] = f"calculation failed: {type(exc).__name__}: {exc}"
    candidate["metrics"] = metrics
    return candidate


def _cache_valid(output: dict) -> bool:
    try:
        path = _path(output["candidates_jsonl"])
        root = _path(output["run_dir"])
        if read_json(root / "metadata.json").get("status") != "completed" or read_json(root / "result.json").get("success") is not True:
            return False
        for row in _absolute_candidates(output):
            structure = row.get("complex_pdb") or row.get("binder_pdb")
            if structure and not Path(structure).is_file():
                return False
        if output.get("target", {}).get("target_pdb") and not Path(output["target"]["target_pdb"]).is_file():
            return False
        return path.is_file()
    except (OSError, ValueError, KeyError, TypeError):
        return False


def cancel_graph_children(run_dir: Path) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
    for node in read_json(run_dir / "graph_execution.json").get("nodes", {}).values():
        raw = node.get("job_run_dir")
        if not raw:
            continue
        child = Path(raw)
        metadata = read_json(child / "metadata.json")
        if metadata.get("pipeline_parent_run_dir") == str(run_dir) and metadata.get("status") in {"queued", "running", "pending"}:
            stop_job(child, reason="Parent pipeline was cancelled")


def enqueue_graph_pipeline(graph: dict, *, pipeline_name: str, owner_id: str, project_id: str, reuse_completed=True, gpu_device="all") -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import enqueue_workflow_call
    graph = validate_graph(graph, require_inputs=True)["graph"]
    return enqueue_workflow_call(run_graph_pipeline, task_group="campaign", tool="visual_pipeline", job_type="visual_pipeline",
        kwargs={"graph": graph, "pipeline_name": pipeline_name, "owner_id": owner_id, "project_id": project_id, "reuse_completed": reuse_completed, "gpu_device": gpu_device}, cpu_cores=1, gpu_device="none")


def run_graph_pipeline(graph: dict, pipeline_name: str, owner_id: str, project_id: str, reuse_completed: bool = True, gpu_device: str = "all") -> Path:
    plan = validate_graph(graph, require_inputs=True)
    graph = plan["graph"]
    job = create_job("campaign", "visual_pipeline", "visual_pipeline", {}, {"pipeline_name": pipeline_name, "graph": graph, "reuse_completed": reuse_completed, "gpu_device": "none"})
    parent = job.run_dir
    from mn_molecular_design_studio.workspace import mark_studio_job
    mark_studio_job({"job": {"run_dir": str(parent)}}, owner_id=owner_id, project_id=project_id)
    snapshot = parent / "graph_snapshot.json"
    if snapshot.exists() and read_json(snapshot) != graph:
        raise ValueError("A submitted pipeline snapshot cannot be edited. Submit the changed draft as a new run.")
    if not snapshot.exists():
        with snapshot.open("x") as handle:
            json.dump(graph, handle, indent=2)
    state = read_json(parent / "graph_execution.json") or {"nodes": {}, "graph_id": graph["id"]}
    node_states, outputs = state["nodes"], {}
    nodes = {n["id"]: n for n in graph["nodes"]}
    cache_root = scope_root(owner_id, project_id) / "cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    cancelled = [False]
    old_handler = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGTERM, lambda *_: cancelled.__setitem__(0, True))
    def persist():
        state["updated_at"] = utc_now(); write_json(parent / "graph_execution.json", state)
    def complete(node_id, output, key, status="completed"):
        output["fingerprint"] = fingerprint(key, {"candidate_digest": file_digest(Path(output["candidates_jsonl"]))}, []) if output.get("candidates_jsonl") else key
        outputs[node_id] = output
        node_states.setdefault(node_id, {}).update(status=status, output=output, fingerprint=key, candidate_count=output.get("candidate_count", 0), completed_at=utc_now())
        if output.get("candidates_jsonl") and status != "cached":
            write_json(cache_root / f"{key}.json", {"output": output})
        persist()
    try:
        update_status(parent, "running", campaign_name=pipeline_name, current_phase="Scheduling pipeline steps")
        while len(outputs) < len(nodes):
            if cancelled[0] or read_json(parent / "metadata.json").get("status") == "cancelled":
                cancel_graph_children(parent)
                update_status(parent, "cancelled", progress_label="Pipeline cancelled; completed results preserved")
                return parent
            progressed = False
            for node_id in plan["order"]:
                if node_id in outputs:
                    continue
                node, parents = nodes[node_id], plan["incoming"][node_id]
                current = node_states.setdefault(node_id, {"status": "pending"})
                if current.get("status") in {"failed", "blocked"}:
                    continue
                if any(node_states.get(p, {}).get("status") in {"failed", "blocked"} for p in parents):
                    current.update(status="blocked", error="An upstream step failed."); persist(); progressed=True; continue
                if not all(p in outputs for p in parents):
                    continue
                kind, params = node["data"]["kind"], node["data"]["params"]
                upstream = [outputs[p] for p in parents]
                source = upstream[0] if upstream else {}
                if kind in {"ligandmpnn", "foundry_mpnn"} and len(upstream) > 1:
                    merged_rows = []
                    merged_target = next((item.get("target") for item in upstream if item.get("target")), {})
                    for source_id, output in zip(parents, upstream):
                        for row in _absolute_candidates(output):
                            candidate = copy.deepcopy(row)
                            original_id = str(candidate.get("candidate_id") or "candidate")
                            candidate["candidate_id"] = f"{source_id}:{original_id}"
                            raw = dict(candidate.get("raw_metadata") or {})
                            raw["pipeline_generator_node"] = source_id
                            raw["pipeline_original_candidate_id"] = original_id
                            candidate["raw_metadata"] = raw
                            merged_rows.append(candidate)
                    source = _local_source(parent, f"{node_id}_fan_in", merged_rows, merged_target)
                key_params = copy.deepcopy(params)
                if kind == "target" or kind=='complexa_generator' and params.get('target_pdb'):
                    target_path = _path(params["target_pdb"])
                    key_params["structure_digest"] = file_digest(target_path)
                if kind in {"design_source","generation_source","backbone_source"}:
                    loader=source_generations if kind in {'generation_source','backbone_source'} else source_designs
                    extra={'stage':'generation.backbone' if kind=='backbone_source' else 'generation.backbone_sequence'} if kind!='design_source' else {}
                    source_dir, _ = loader(params["source_run_dir"], project_id, params.get("selected_candidate_ids"),**extra)
                    key_params["source_digest"] = file_digest(source_dir / "artifacts/normalized_candidates/candidates.jsonl")
                if kind=='complexa_generator':
                    from mn_molecular_design_studio import complexa_input, complexa_container
                    key_params['input_adapter_digest']=file_digest(Path(complexa_input.__file__))
                    key_params['container_adapter_digest']=file_digest(Path(complexa_container.__file__))
                from mn_molecular_design_studio import generator_input
                from mn_molecular_design_studio import inverse_folding_input
                if kind in inverse_folding_input.ENGINES:
                    from mn_molecular_design_studio.engine.workflow_engine.workflows import inverse_folding
                    key_params["input_adapter_digest"] = file_digest(Path(inverse_folding_input.__file__))
                    key_params["container_adapter_digest"] = file_digest(inverse_folding.NATIVE)
                if kind in generator_input.ENGINES:
                    key_params['input_adapter_digest']=file_digest(Path(generator_input.__file__))
                    for key in ('target_pdb','scaffold_pdb','scaffold_ss_file','scaffold_adj_file','target_ss_file','target_adj_file'):
                        if params.get(key):key_params[key+'_digest']=file_digest(_path(params[key]))
                if NODE_CATALOG[kind]["workflow"]:
                    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function
                    implementation = _get_function(WORKFLOWS[NODE_CATALOG[kind]["workflow"]])
                    key_params["implementation_digest"] = file_digest(Path(inspect.getfile(implementation)))
                key = fingerprint(kind, key_params, [o["fingerprint"] for o in upstream])
                try:
                    if current.get("job_run_dir"):
                        child = Path(current["job_run_dir"])
                        metadata = read_json(child / "metadata.json")
                        current["status"] = str(metadata.get("status") or "queued")
                        current["job_code"] = metadata.get("job_code")
                        if current["status"] == "completed" and read_json(child / "result.json").get("success") is True:
                            output = _child_output(child, source.get("target", {}))
                            if NODE_CATALOG[kind]["category"] == "Structure prediction":
                                output = _attach_refold_input_references(output, source, parent, node_id)
                            complete(node_id, output, key); progressed=True
                        elif current["status"] in {"failed", "cancelled"}:
                            current.update(status="failed", error=f"Child job {current.get('job_code') or child.name} did not complete successfully.");persist();progressed=True
                        continue
                    cached = read_json(cache_root / f"{key}.json").get("output") if reuse_completed else None
                    if cached and _cache_valid(cached):
                        if NODE_CATALOG[kind]["category"] == "Structure prediction":
                            cached = _attach_refold_input_references(cached, source, parent, node_id)
                        complete(node_id, cached, key, "cached");progressed=True;continue
                    if kind == "target":
                        chains = [x for x in re.split(r"[\s,;]+", params["target_chains"]) if x]
                        if not chains:
                            raise ValueError("Select at least one target chain.")
                        staged_target = parent / "artifacts" / "pipeline_inputs" / node_id / ("target" + target_path.suffix)
                        staged_target.parent.mkdir(parents=True, exist_ok=True)
                        if not staged_target.exists():
                            shutil.copy2(target_path, staged_target)
                        output = {"target": {"target_pdb": str(staged_target), "target_chains": chains, "hotspots": params["hotspots"]}}
                        complete(node_id, output, key);progressed=True;continue
                    if kind in {"design_source","generation_source","backbone_source"}:
                        loader=source_generations if kind in {'generation_source','backbone_source'} else source_designs
                        extra={'stage':'generation.backbone' if kind=='backbone_source' else 'generation.backbone_sequence'} if kind!='design_source' else {}
                        source_dir, selected = loader(params["source_run_dir"], project_id, params.get("selected_candidate_ids"),**extra)
                        all_rows = _absolute_candidates({"run_dir": str(source_dir), "candidates_jsonl": str(source_dir / "artifacts/normalized_candidates/candidates.jsonl")})
                        selected_ids = {c["candidate_id"] for c in selected}
                        if kind in {'generation_source','backbone_source'}:
                            source_input = read_json(source_dir/'input.json')
                            source_params = source_input.get('params', {}) if isinstance(source_input.get('params'), dict) else {}
                            source_inputs = source_input.get('inputs', {}) if isinstance(source_input.get('inputs'), dict) else {}
                            target_inputs = dict(source_inputs)
                            if not target_inputs.get('target_pdb'):
                                target_inputs['target_pdb'] = target_inputs.get('target_structure')
                            target = _generation_target({**target_inputs, **source_params})
                            # Candidate paths have already been resolved relative to the
                            # source run, so prefer them over external/original upload paths.
                            candidate_target = next((row.get('target_pdb') for row in all_rows if row.get('target_pdb')), None)
                            if candidate_target:
                                target['target_pdb'] = candidate_target
                            elif target.get('target_pdb'):
                                resolved_target = resolve_stored_path(target['target_pdb'], run_dir=source_dir)
                                if resolved_target:
                                    target['target_pdb'] = str(resolved_target)
                            if not target.get('target_chains'):
                                target['target_chains'] = next((row.get('target_chains') for row in all_rows if row.get('target_chains')), [])
                            if not target.get('hotspots'):
                                target['hotspots'] = next((row.get('hotspots') for row in all_rows if row.get('hotspots')), [])
                        else:
                            target = {}
                        complete(node_id, _local_source(parent, node_id, [c for c in all_rows if c["candidate_id"] in selected_ids],target), key);progressed=True;continue
                    if kind in {"filter", "results"}:
                        rows = []
                        for source_id, output in zip(parents, upstream):
                            for row in _absolute_candidates(output):
                                if kind == "results":
                                    row = _add_target_aligned_pose_rmsd(row, parent)
                                    row.setdefault("raw_metadata", {}).update(pipeline_source_node=source_id, pipeline_original_candidate_id=row["candidate_id"])
                                    row["candidate_id"] = f"{source_id}:{row['candidate_id']}"
                                rows.append(row)
                        if kind == "filter":
                            rows = filter_candidates(rows, params)
                        complete(node_id, _local_source(parent, node_id, rows, source.get("target")), key);progressed=True;continue
                    if source.get("candidates_jsonl") and not source.get("candidate_count"):
                        complete(node_id, _local_source(parent, node_id, [], source.get("target")), key, "skipped");progressed=True;continue
                    request = node_request(kind, params, source, f"{pipeline_name} · {node['data']['label']}")
                    child = _submit(request, parent, node, owner_id=owner_id, project_id=project_id, gpu_device=gpu_device)
                    metadata = read_json(child / "metadata.json")
                    current.update(status="queued", job_run_dir=str(child), job_code=metadata.get("job_code"), fingerprint=key, started_at=utc_now())
                    persist();progressed=True
                except Exception as exc:
                    current.update(status="failed", error=str(exc));persist();progressed=True
            unresolved = [n for n in nodes if n not in outputs and node_states.get(n, {}).get("status") not in {"failed", "blocked"}]
            if not unresolved:
                break
            if not progressed:
                persist();time.sleep(2)
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
        candidates=[]
        for node_id, output in outputs.items():
            if nodes[node_id]["data"]["kind"] != "results":
                continue
            for candidate in _absolute_candidates(output):
                candidate["candidate_id"] = f"{node_id}:{candidate['candidate_id']}"
                candidate.setdefault("raw_metadata", {})["pipeline_result_node"] = node_id
                candidates.append(candidate)
        write_candidates(parent, "visual_pipeline", candidates)
        failures = [v for v in node_states.values() if v.get("status") in {"failed", "blocked"}]
        state["status"] = "failed" if failures else "completed";persist()
        finish_job(parent, not failures, {"outputs": {"candidates": candidates, "pipeline_nodes": node_states}, "metrics": {"candidate_count": len(candidates), "node_count": len(nodes), "cached_steps": sum(v.get("status")=="cached" for v in node_states.values()), "failed_steps": len(failures)}, "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl", "graph_snapshot": "graph_snapshot.json", "graph_execution": "graph_execution.json"}})
        return parent
    except Exception:
        cancel_graph_children(parent)
        raise
    finally:
        signal.signal(signal.SIGTERM, old_handler)
