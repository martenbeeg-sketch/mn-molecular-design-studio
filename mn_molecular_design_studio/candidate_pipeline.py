"""Composition and retained-evidence handoff for the ordinary Studio graph."""
from __future__ import annotations

import copy
import json
import re
import tempfile
import uuid
from pathlib import Path

from mn_molecular_design_studio import candidate_sets as cs
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path


def new_pipeline(generators=None, predictors=None, sequence_designer="ligandmpnn", *, mode=None, analyses=None, redesign_engine="ligandmpnn"):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, SCHEMA
    from mn_molecular_design_studio.workflow_modes import MODES, capability, cofolder_methods, methods, sequence_count_field, sequence_methods
    if mode is not None and mode not in MODES:
        raise ValueError("Unknown scientific workflow mode.")
    defaults = ["rfd4", "rfdiffusion"] if mode == "modular_design" else ["bindcraft2"] if mode else ["rfd4", "bindcraft2"]
    generators = defaults if generators is None else generators
    if mode and any(kind not in methods(mode) for kind in generators):
        raise ValueError("Selected operation is unsupported for this workflow mode.")
    predictors = ["boltz2", "prediction_protenix_v2"] if predictors is None else predictors
    if mode == "native_complete":
        predictors = []
    if not generators or (not predictors and mode != "native_complete"):
        raise ValueError("Choose at least one generator and predictor.")
    for kind in generators:
        if kind not in NODE_CATALOG or NODE_CATALOG[kind]["category"] not in {"Design engine", "Design pipeline"}:
            raise ValueError("Choose existing generator tiles.")
        operation = capability(kind)
        if not operation['target_connected'] or not operation['modular_source']:
            raise ValueError(f"{NODE_CATALOG[kind]['label']} requires its own input workflow; it cannot receive a shared structure target.")
    if mode == 'native_complete' and analyses:
        raise ValueError('Complete Native preserves the selected workflows unchanged; add independent analyses in Native + Comparative or Modular mode.')
    for kind in predictors:
        if kind not in NODE_CATALOG or NODE_CATALOG[kind]["category"] != "Structure prediction":
            raise ValueError("Choose existing predictor tiles.")
        if mode and kind not in cofolder_methods('protein_protein'):
            raise ValueError(f"{NODE_CATALOG[kind]['label']} does not declare compatible protein–protein co-folding inputs.")
    analyses = list(analyses or [])
    for kind in analyses:
        if kind not in NODE_CATALOG or NODE_CATALOG[kind]['category'] != 'Analysis' or not NODE_CATALOG[kind].get('workflow'):
            raise ValueError('Choose registered Analysis operations.')
    if NODE_CATALOG[sequence_designer]["category"] != "Sequence design":
        raise ValueError("Choose an existing sequence-design tile.")
    if mode == 'modular_design' and sequence_designer not in sequence_methods():
        raise ValueError('Choose a registered sequence designer compatible with normalized candidate inputs.')
    if NODE_CATALOG[redesign_engine]['category'] != 'Sequence design':
        raise ValueError('Choose an existing sequence-design tile for optional redesign.')
    if mode == 'modular_design' and redesign_engine not in sequence_methods():
        raise ValueError('Choose a registered sequence designer compatible with optional candidate redesign.')
    if mode in {"native_comparative", "modular_design"} and len(generators) > 1:
        # Each completed generator feeds its own stages immediately. Only final
        # comparative selection waits for fan-in; no generator-wide barrier.
        branches = [new_pipeline([kind], predictors, sequence_designer, mode=mode, analyses=analyses,
                                 redesign_engine=redesign_engine) for kind in generators]
        branch_outputs = [next(edge['source'] for edge in branch['edges'] if edge['target'] == 'ranking')
                          for branch in branches]
        graph = branches[0]
        graph['nodes'] = [n for n in graph['nodes'] if n['id'] not in {'ranking', 'results'}]
        graph['edges'] = [e for e in graph['edges'] if not {e['source'], e['target']} & {'ranking', 'results'}]
        analyzers = [branch_outputs[0]]
        for i, branch in enumerate(branches[1:], 1):
            rename = {n['id']: n['id'] + '-g' + str(i) for n in branch['nodes']}
            rename['target'] = 'target'
            for node in branch['nodes']:
                old = node['id']
                if old in {'target', 'ranking', 'results'}:
                    continue
                if NODE_CATALOG[node['data']['kind']]['category'] == 'Structure prediction':
                    node['data']['evidence_namespace'] = old
                node['id'] = rename[old]
                node['position']['y'] += i * 560
                graph['nodes'].append(node)
            for edge in branch['edges']:
                if {edge['source'], edge['target']} & {'ranking', 'results'}:
                    continue
                edge.update(id=edge['id'] + '-g' + str(i), source=rename[edge['source']], target=rename[edge['target']])
                graph['edges'].append(edge)
            analyzers.append(rename[branch_outputs[i]])
        def final_node(kind, identifier, x):
            spec = NODE_CATALOG[kind]
            graph['nodes'].append(dict(id=identifier, type='studio', position=dict(x=x, y=170),
                data=dict(kind=kind, label=spec['label'], params={f['name']: copy.deepcopy(f['default']) for f in spec['fields']},
                          execution=dict(failure_policy='partial'))))
        final_node('candidate_merge', 'evaluated-collection', 2740)
        final_node('candidate_rank', 'ranking', 3080)
        final_node('results', 'results', 3420)
        for source, target in [*((node, 'evaluated-collection') for node in analyzers), ('evaluated-collection', 'ranking'), ('ranking', 'results')]:
            graph['edges'].append(dict(id=source+'--'+target, source=source, target=target, sourceHandle='out', targetHandle='in'))
        return graph
    nodes, edges = [], []
    def add(kind, node_id, x, y, **params):
        spec = NODE_CATALOG[kind]
        defaults = {f["name"]: copy.deepcopy(f["default"]) for f in spec["fields"]}
        data = dict(kind=kind, label=spec["label"], params=defaults | params)
        nodes.append(dict(id=node_id, type="studio", position=dict(x=x, y=y), data=data))
        return data
    def connect(source, target):
        edges.append(dict(id=source + "--" + target, source=source, target=target, sourceHandle="out", targetHandle="in"))
    add("target", "target", 0, 170)
    for i, kind in enumerate(generators):
        add(kind, "generator-" + str(i), 340, i * 240)["execution"] = dict(mode="native", placement="auto", failure_policy="partial")
        if mode and kind == 'rfd4':
            nodes[-1]['data']['execution'].update(cpu_cores=8, min_vram_gb=16)
        connect("target", "generator-" + str(i))
        connect("generator-" + str(i), "collection")
    add("candidate_merge", "collection", 700, 170)["execution"] = dict(failure_policy="partial")
    if mode == "native_complete":
        nodes[-1]["data"]["label"] = "Native outcomes"
        add("results", "results", 1040, 170)["label"] = "Native Results"
        connect("collection", "results")
        return dict(schema_version=SCHEMA, id=uuid.uuid4().hex, name=MODES[mode],
                    template_id=mode, workflow_mode=mode, nodes=nodes, edges=edges,
                    viewport=dict(x=0, y=0, zoom=.65))
    initial_count = {sequence_count_field(sequence_designer): 1}
    redesign_count = {sequence_count_field(redesign_engine): 1}
    add(sequence_designer, "sequence-preparation", 1040, 170, sequence_stage="initial", **initial_count)
    add(redesign_engine, "optional-redesign", 1380, 170, sequence_stage="soluble_optimization",
        optimization_enabled=False, **redesign_count)
    if mode == 'modular_design' and redesign_engine == 'ligandmpnn':
        nodes[-1]['data']['params'].update(model_type='soluble_mpnn', fix_interface_residues=True)
    add("shared_filter", "shared-filter", 1720, 170)
    connect("collection", "sequence-preparation")
    connect("sequence-preparation", "optional-redesign")
    connect("optional-redesign", "shared-filter")
    for i, kind in enumerate(predictors):
        data = add(kind, "predictor-" + str(i), 2060, i * 240)
        data["execution"] = dict(mode="sharded", batch_size=5, max_in_flight=2, failure_policy="partial", placement="auto")
        # A second compatible adapter is available for explicit opt-in; no
        # extra GPU work is implied by opening the template.
        data["enabled"] = i == 0
        if mode and "seed" in data["params"]:
            data["params"]["seed"] = 42  # Same explicit default protocol across generator branches.
        connect("shared-filter", "predictor-" + str(i))
        connect("predictor-" + str(i), "analyzer")
    add("candidate_analyzer", "analyzer", 2400, 170)["execution"] = dict(failure_policy="partial")
    connect("collection", "analyzer")
    connect("shared-filter", "analyzer")
    analysis_nodes = []
    for index, kind in enumerate(analyses):
        node_id = f'analysis-{index}'
        add(kind, node_id, 2580 + index * 320, 170)
        connect('analyzer', node_id)
        analysis_nodes.append(node_id)
    final_source = 'analyzer'
    rank_x = 2740
    if analysis_nodes:
        add('candidate_analyzer', 'analysis-analyzer', 2580 + len(analysis_nodes) * 320, 170)['execution'] = dict(failure_policy='partial')
        for source in ['analyzer', *analysis_nodes]:
            connect(source, 'analysis-analyzer')
        final_source = 'analysis-analyzer'
        rank_x = 2940 + len(analysis_nodes) * 320
    add("candidate_rank", "ranking", rank_x, 170)
    add("results", "results", rank_x + 340, 170)
    connect(final_source, "ranking")
    connect("ranking", "results")
    if mode is None:
        for node in nodes:
            if 'execution' in node['data']:
                node['data']['execution'].pop('placement', None)
                if NODE_CATALOG[node['data']['kind']]['category'] in {'Design engine', 'Design pipeline'}:
                    node['data']['execution'] = dict(mode='native')
        next(n for n in nodes if n['id'] == 'collection')['data'].pop('execution', None)
    if mode == "native_comparative":
        nodes = [n for n in nodes if n["id"] not in {"sequence-preparation", "optional-redesign"}]
        edges = [e for e in edges if not {e["source"], e["target"]} & {"sequence-preparation", "optional-redesign"}]
        connect("collection", "shared-filter")
    if mode:
        return dict(schema_version=SCHEMA, id=uuid.uuid4().hex, name=MODES[mode],
                    template_id=mode, workflow_mode=mode, nodes=nodes, edges=edges,
                    viewport=dict(x=0, y=0, zoom=.45))
    return dict(schema_version=SCHEMA, id=uuid.uuid4().hex,
                name="Multi-Generator Binder Design & Evaluation", template_id="multi_generator_binder",
                nodes=nodes, edges=edges, viewport=dict(x=0, y=0, zoom=.45))


def retained_source(node, project_id):
    """Read-only explicit reuse, separate from native controls and graph caching."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, OPTIMIZED_PREDICTION_NODE_ENGINES
    path = resolve_stored_path(node["data"]["retained_run_dir"], must_exist=True)
    if path is None or not path.is_dir():
        raise ValueError("The retained source run is unavailable.")
    marker = read_json(path / "molecular_design_studio.json")
    if marker.get("project_id") != project_id:
        raise ValueError("Retained evidence must belong to the selected project.")
    meta, request, result = (read_json(path / f) for f in ("metadata.json", "input.json", "result.json"))
    if meta.get("status") not in {"completed", "succeeded"} or result.get("success") is not True:
        raise ValueError("Choose a completed retained run with successful execution; native rejections remain valid evidence.")
    kind = node["data"]["kind"]
    workflow = NODE_CATALOG[kind]["workflow"]
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
    workflow_spec = WORKFLOWS.get(workflow)
    if workflow_spec is None:
        raise ValueError(f"The {NODE_CATALOG[kind]['label']} tile has no registered native workflow route.")
    expected_tool = workflow_spec.tool
    params = request.get("params", {})
    predictor = NODE_CATALOG[kind]["category"] == "Structure prediction"
    candidate_path = path / "artifacts/normalized_candidates/candidates.jsonl"
    if not candidate_path.is_file():
        raise ValueError("Retained source has no normalized candidate collection (including an explicit empty collection).")
    tool = meta.get("tool") or request.get("tool")
    embedded = predictor and tool in {"rfd4_proteina", "bindcraft2", "distributed_campaign"}
    engine = OPTIMIZED_PREDICTION_NODE_ENGINES.get(kind, kind)
    if predictor and not embedded:
        from mn_molecular_design_studio.output_vocabulary import engine_name
        if params.get("engine") != engine and engine_name(tool) != engine:
            raise ValueError("Retained evidence does not match the predictor tile.")
    elif not predictor and tool != expected_tool:
        raise ValueError("Retained evidence does not match the generator tile.")
    provenance = dict(run_id=path.name, workflow=workflow, retained=True,
        input=request, metadata=meta, result_metrics=result.get("metrics", {}),
        source_files={str(path / name): cs.file_digest(path / name) for name in ("input.json", "metadata.json", "result.json") if (path / name).is_file()})
    return dict(run_dir=str(path), candidates_jsonl=str(candidate_path), source_run_id=path.name,
                source_provenance=provenance, embedded_predictions=embedded, expected_predictor_engine=engine)


def preflight(graph, plan, project_id):
    """Check the common target and generator handoffs before creating any jobs."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
    if not any(n["data"]["kind"] in cs.KINDS for n in graph["nodes"]):
        return
    nodes = {n["id"]: n for n in graph["nodes"]}
    contexts = {}
    for node_id in plan["order"]:
        node = nodes[node_id]
        kind, params = node["data"]["kind"], node["data"]["params"]
        parents = plan["incoming"][node_id]
        target = copy.deepcopy(contexts.get(parents[0], {})) if parents else {}
        if kind == "target":
            target = copy.deepcopy(params)
            target["target_chains"] = [s for s in re.split(r"[,;\s]+", target["target_chains"]) if s]
            target["identity"] = cs.target_identity(target)
        if kind == "binder_context":
            for key in ("hotspots", "coldspots"):
                target[key] = params.get(key) if params.get("clear_unselected") else params.get(key) or target.get(key, "")
        contexts[node_id] = target
        if node["data"].get("enabled") is False:
            continue
        if node["data"].get("retained_run_dir"):
            retained = retained_source(node, project_id)
            if NODE_CATALOG[kind]["category"] in {"Design engine", "Design pipeline"}:
                if not target.get("identity"):
                    raise ValueError("A retained generator branch requires the common target input.")
                for row in cs.iter_candidates(retained):
                    verify_target(row, target)
            continue
        if NODE_CATALOG[kind]["category"] not in {"Design engine", "Design pipeline"}:
            continue
        if not target.get("identity"):
            raise ValueError("Connect the common target to every generator branch.")
        for own_key in ("target_structure", "target_pdb"):
            if params.get(own_key) and cs.file_digest(params[own_key]) != target["identity"]["source_sha256"]:
                raise ValueError("A generator overrides the common target with a different structure.")
        for constraint, support in (("hotspots", "hotspot_support"), ("coldspots", "coldspot_support")):
            supported = NODE_CATALOG[kind].get(support)
            if kind == "rfd4" and constraint == "hotspots":
                supported = "native"
            if target.get(constraint) and supported not in {"native", "mapped", "conditional"}:
                raise ValueError(f"{NODE_CATALOG[kind]['label']} does not accept the selected {constraint}; configure its branch context explicitly.")
        request = node_request(kind, params, {"target": target}, "pipeline-preflight")
        # Native request validation performs the existing engine-specific input
        # checks, without containers, scheduling or scientific output generation.
        from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import validate_workflow_request
        with tempfile.TemporaryDirectory(prefix="studio-pipeline-preflight-") as directory:
            path = Path(directory) / "request.json"
            path.write_text(json.dumps(request))
            validate_workflow_request(path)


def verify_target(row, target):
    """Compare exact selected target sequences with the existing verified role reader."""
    from mn_molecular_design_studio.design_validation import read_structure, protein_chains, sequence
    if not (row.get("complex_pdb") or row.get("binder_pdb")):
        return
    chains = protein_chains(read_structure(Path(row.get("complex_pdb") or row["binder_pdb"])))
    targets, binders = row.get("target_chains", []), row.get("binder_chains", [])
    if not targets or not binders or set(targets) & set(binders) or not set(targets + binders) <= chains.keys():
        raise ValueError("Retained candidate has unresolved target/binder chain roles.")
    sequences = [sequence(chains[chain]) for chain in targets]
    expected = [v["sequence"] for v in target["identity"]["chains"].values()]
    if sequences != expected:
        raise ValueError("Retained candidate target-chain sequence differs from the shared target.")
    if cs.sequence_state(row) == "sequence_ready" and ":".join(sequence(chains[chain]) for chain in binders).replace(":", "") != row["binder_sequence"].replace(":", ""):
        raise ValueError("Retained candidate binder sequence differs from its structure.")


def replay(parent, node, source, retained, scientific_key):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    kind, node_id = node["data"]["kind"], node["id"]
    predictor = NODE_CATALOG[kind]["category"] == "Structure prediction"
    target = source["target"]
    if not predictor:
        def rows():
            for row in cs.iter_candidates(retained):
                verify_target(row, target)
                yield row
        # Deliberately no candidate-set annotation here: the merge is the
        # generic output adapter for ordinary and retained generator sources.
        result = cs.write_collection(parent, node_id + "_retained", rows(), target)
        result.pop("candidate_set_manifest")
        result.update(source_run_id=retained["source_run_id"], source_provenance=retained["source_provenance"])
        result["stage_counts"].update(cs.native_counts(retained["source_provenance"]))
        return result
    with cs.index() as db:
        db.execute("CREATE TABLE inputs (native_id TEXT, seq TEXT, payload TEXT)")
        for row in cs.iter_candidates(source):
            info = cs.metadata(row)
            db.execute("INSERT INTO inputs VALUES (?,?,?)", (info["native_candidate_id"], row.get("binder_sequence"), json.dumps(row)))
        def predictions():
            for original in cs.iter_candidates(retained):
                rows = original.get("raw_metadata", {}).get("evidence_candidates", []) if retained["embedded_predictions"] else [original]
                for row in rows:
                    raw = row.get("raw_metadata", {})
                    engine = retained["expected_predictor_engine"]
                    if engine not in {raw.get("engine"), str(row.get("source_tool", "")).removeprefix("optimized_prediction_")}:
                        continue
                    row = copy.deepcopy(row)
                    # Resolve nested evidence paths through the ordinary portable-path reader.
                    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_managed_paths
                    row = resolve_managed_paths(row, run_dir=Path(retained["run_dir"]))
                    for key in ("complex_pdb", "binder_pdb", "target_pdb"):
                        if row.get(key):
                            path = resolve_stored_path(row[key], run_dir=Path(retained["run_dir"]), must_exist=True)
                            if path:
                                row[key] = str(path)
                    ids = [*row.get("parents", []), raw.get("source_candidate", {}).get("candidate_id")]
                    if retained["embedded_predictions"]:
                        ids.append(original["candidate_id"])
                    matches = db.execute("SELECT payload FROM inputs WHERE native_id IN (" + ",".join("?" for _ in ids) + ") AND seq=?", (*ids, row.get("binder_sequence"))).fetchall() if ids else []
                    if not matches:
                        continue  # unrelated retained predictions are never borrowed by sequence alone
                    if len(matches) != 1:
                        raise ValueError("Retained prediction has ambiguous source scientific ancestry.")
                    reference = json.loads(matches[0][0])
                    verify_target(row, target)
                    info = copy.deepcopy(cs.metadata(reference))
                    info.update(record_kind="prediction", predictor_node=node["data"].get("evidence_namespace", node_id),
                                predictor_execution_node=node_id,
                                prediction_id=cs.identity(node_id, retained["source_run_id"], row["candidate_id"], info["variant_id"]),
                                predictor_settings=retained["source_provenance"]["input"].get("params", {}),
                                retained_provenance=retained["source_provenance"])
                    row.setdefault("raw_metadata", {})["candidate_set"] = info
                    row["raw_metadata"]["pipeline_refold_input"] = dict(
                        complex_pdb=reference.get("complex_pdb") or reference.get("binder_pdb"),
                        binder_chains=reference.get("binder_chains", []), target_chains=reference.get("target_chains", []),
                        source_candidate_id=reference["candidate_id"])
                    yield row
        return cs.write_collection(parent, node_id + "_retained", predictions(), target)


def execution_plan(graph, candidate_counts=None):
    """Read-only native execution forecast. Counts are eligible inputs per tile.

    Unknown generator yields, sequence routing and compatibility groups stay
    unknown. Preview counts never change saved scientific or scheduler settings.
    """
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, OPTIMIZED_PREDICTION_NODE_ENGINES
    from mn_molecular_design_studio.candidate_execution import policy, supports, dispatch_limits
    from mn_molecular_design_studio.native_execution import prediction_contract, request_contract, prediction_engine
    from mn_molecular_design_studio.engine.workflow_engine.workflows.execution_planning import describe_native_invocation
    from mn_molecular_design_studio.batching import batch_plan
    counts, plan = candidate_counts or {}, []
    for node in graph["nodes"]:
        data = node["data"]
        kind, params = data["kind"], data.get("params", {})
        spec = NODE_CATALOG[kind]
        settings = policy(data.get("execution", graph.get("candidate_execution")))
        sharded = settings["mode"] == "sharded" and supports(kind, {"candidates_jsonl": "plan"})
        count = counts.get(node["id"])
        if count is not None and (type(count) is not int or count < 0):
            raise ValueError("Preview candidate counts must be nonnegative integers.")
        preserve_sequences = params.get('sequence_stage') == 'soluble_optimization' and not params.get('optimization_enabled')
        inactive = not data.get("enabled", True) or bool(data.get("retained_run_dir")) or not spec["workflow"] or preserve_sequences
        chunks = 0 if inactive else None if count is None else batch_plan(count, settings['batch_size'])['batch_count'] if sharded else int(count > 0)
        row = dict(node=node["id"], tool=spec["label"], enabled=data.get("enabled", True),
            execution="disabled" if not data.get("enabled", True) else "preserve sequences" if preserve_sequences else "retained evidence" if data.get("retained_run_dir") else "collection chunks" if sharded else "ordinary native workflow" if spec["workflow"] else "collection transformation",
            candidates=count, collection_chunks=chunks, expected_jobs=chunks,
            candidates_per_chunk=settings["batch_size"] if sharded else None,
            items_per_batch=settings['batch_size'] if sharded else None,
            batch_count=chunks,
            dispatch_limits=dispatch_limits(kind, settings) if sharded else {},
            max_concurrent_jobs=settings["max_in_flight"] if sharded else 1,
            placement=settings["placement"] if sharded or spec["workflow"] else "local",
            native_invocations=None, docker_processes=None, gpu_allocation_events=None,
            gpu_slots=None, prediction_outputs=None, model_initializations=None,
            note="Input counts unknown until upstream outputs and sequence routing resolve.")
        if inactive:
            row.update(native_invocations=0, docker_processes=0, gpu_allocation_events=0, gpu_slots=0,
                       prediction_outputs=0, model_initializations=0, note="No new scientific execution.")
        elif kind in {"rfd4", "bindcraft2", "rfdiffusion", "rfd3"}:
            engine = {"rfd4": "rfd4_proteina", "rfd3": "rfdiffusion3"}.get(kind, kind)
            attempts = int(params.get("num_designs", 1))
            input_files = params.get('input_files') or []
            input_count = max(1, sum(bool(str(value).strip()) for value in (input_files.splitlines() if isinstance(input_files, str) else input_files))) if kind == 'rfd4' else 1
            attempt_chunks = (attempts + settings['batch_size'] - 1) // settings['batch_size'] if sharded else 1
            native = describe_native_invocation(engine, attempts_per_input=attempts, input_count=input_count,
                         native_batch_size=int(params.get("batch_size", params.get("diffusion_batch_size", 1))))
            row.update(requested_generator_attempts=attempts * input_count, attempts_per_input=attempts, generator_input_count=input_count,
                       gpu_microbatch=params.get("batch_size", params.get("diffusion_batch_size")),
                       expected_jobs=attempt_chunks, batch_count=attempt_chunks, generator_attempt_chunks=attempt_chunks, collection_chunks=0, docker_processes=attempt_chunks,
                       native_invocations=(sum((min(settings['batch_size'], attempts-offset)+99)//100 for offset in range(0, attempts, settings['batch_size'])) if sharded and kind == 'rfd4' else attempt_chunks if sharded else native["native_python_invocations"]),
                       gpu_allocation_events=attempt_chunks, gpu_slots=attempt_chunks, prediction_outputs=0,
                       model_initializations=None,
                       note="Ordinary native attempts retained. BindCraft2 initialization is per worker; auto worker count and validation rebuilds are unknown. RFD4 search/oracle multiplicity is excluded from simple counts.")
            if kind == "rfd4" and (any(json.loads(params.get("options_json") or "{}").get(k) for k in ("refold", "apo", "reward_search", "cluster")) or params.get('avoid_structures') or params.get('also_bind_structures')):
                row.update(native_invocations=None, docker_processes=None, model_initializations=None,
                           note="RFD4 oracle/search enabled: extra process count depends on native execution.")
            if params.get("distributed_campaign"):
                row.update(native_invocations=None, docker_processes=None, expected_jobs=None,
                           note="Acceptance-driven campaign uses its existing planner and one verified attempt per child; this preview does not authorize native campaign batching.")
            elif kind != 'bindcraft2':
                row['model_initializations'] = row['native_invocations']
        else:
            engine = prediction_engine(kind)
            native = (prediction_contract(engine, count or 1,
                        seeds=params.get("seed_count", 1), samples=params.get("samples", 1),
                        mode=params.get("optimization_mode", "off")) if engine else
                      request_contract(dict(workflow=spec["workflow"], parameters=params), count or 1))
            row.update(native_multi_input=native["native_multi_input"], note=native["model_reuse"], count_basis=native["count_basis"])
            if count is not None:
                for key in ("native_invocations", "docker_processes", "gpu_allocation_events", "gpu_slots", "model_initializations"):
                    value = native.get(key)
                    row[key] = 0 if not count else value * chunks if value is not None else None
                row["prediction_outputs"] = (count * int(params.get("samples", 1)) * (int(params.get("seed_count", 1)) if engine == "boltz2" else 1) if engine else native.get('prediction_outputs'))
                if "native_invocations_max" in native:
                    row.update(native_invocations_min=chunks, native_invocations_max=count)
        plan.append(row)
    return plan
