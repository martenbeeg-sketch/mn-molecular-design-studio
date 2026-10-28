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


def new_pipeline(generators=None, predictors=None):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, SCHEMA
    generators = ["rfd4", "bindcraft2"] if generators is None else generators
    predictors = ["boltz2", "prediction_protenix_v2"] if predictors is None else predictors
    if not generators or not predictors:
        raise ValueError("Choose at least one generator and predictor.")
    for kind in generators:
        if NODE_CATALOG[kind]["category"] not in {"Design engine", "Design pipeline"}:
            raise ValueError("Choose existing generator tiles.")
        if not any("target" in socket["accepts"] for socket in NODE_CATALOG[kind]["inputs"]):
            raise ValueError(f"{NODE_CATALOG[kind]['label']} requires its own input workflow; it cannot receive a shared structure target.")
    for kind in predictors:
        if NODE_CATALOG[kind]["category"] != "Structure prediction":
            raise ValueError("Choose existing predictor tiles.")
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
        add(kind, "generator-" + str(i), 340, i * 240)["execution"] = dict(mode="native")
        connect("target", "generator-" + str(i))
        connect("generator-" + str(i), "collection")
    add("candidate_merge", "collection", 700, 170)
    add("ligandmpnn", "sequence-preparation", 1040, 170, sequence_stage="initial", sequences=1)
    add("ligandmpnn", "optional-redesign", 1380, 170, sequence_stage="soluble_optimization", optimization_enabled=False)
    add("shared_filter", "shared-filter", 1720, 170)
    connect("collection", "sequence-preparation")
    connect("sequence-preparation", "optional-redesign")
    connect("optional-redesign", "shared-filter")
    for i, kind in enumerate(predictors):
        data = add(kind, "predictor-" + str(i), 2060, i * 240)
        data["execution"] = dict(mode="sharded", batch_size=5, max_in_flight=2, failure_policy="partial")
        # A second compatible adapter is available for explicit opt-in; no
        # extra GPU work is implied by opening the template.
        data["enabled"] = i == 0
        connect("shared-filter", "predictor-" + str(i))
        connect("predictor-" + str(i), "analyzer")
    add("candidate_analyzer", "analyzer", 2400, 170)["execution"] = dict(failure_policy="partial")
    connect("collection", "analyzer")
    connect("shared-filter", "analyzer")
    add("candidate_rank", "ranking", 2740, 170)
    add("results", "results", 3080, 170)
    connect("analyzer", "ranking")
    connect("ranking", "results")
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
    expected = {"rfd4": "rfd4_proteina", "rfd3": "rfdiffusion3", "complexa": "proteina_complexa",
                "complexa_generator": "proteina_complexa_generate", "protpardelle": "protpardelle_1c"}.get(kind, kind)
    params = request.get("params", {})
    predictor = NODE_CATALOG[kind]["category"] == "Structure prediction"
    candidate_path = path / "artifacts/normalized_candidates/candidates.jsonl"
    if not candidate_path.is_file():
        raise ValueError("Retained source has no normalized candidate collection (including an explicit empty collection).")
    tool = meta.get("tool") or request.get("tool")
    embedded = predictor and tool in {"rfd4_proteina", "bindcraft2", "distributed_campaign"}
    engine = OPTIMIZED_PREDICTION_NODE_ENGINES.get(kind, kind)
    if predictor and not embedded:
        if params.get("engine") != engine and tool not in {expected, "optimized_prediction_" + engine}:
            raise ValueError("Retained evidence does not match the predictor tile.")
    elif not predictor and tool != expected:
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
                    info.update(record_kind="prediction", predictor_node=node_id,
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
    """Read-only collection plan; native attempt batching remains a separate planner."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.candidate_execution import policy, supports
    counts = candidate_counts or {}
    plan = []
    for node in graph["nodes"]:
        kind = node["data"]["kind"]
        spec = NODE_CATALOG[kind]
        settings = policy(node["data"].get("execution", graph.get("candidate_execution")))
        sharded = settings["mode"] == "sharded" and supports(kind, {"candidates_jsonl": "plan"})
        count = counts.get(node["id"])
        plan.append(dict(node=node["id"], tool=spec["label"], enabled=node["data"].get("enabled", True),
            execution="retained evidence" if node["data"].get("retained_run_dir") else "collection chunks" if sharded else "ordinary native workflow" if spec["workflow"] else "collection transformation",
            candidates_per_chunk=settings["batch_size"] if sharded else None,
            expected_jobs=(count + settings["batch_size"] - 1) // settings["batch_size"] if sharded and count is not None else None,
            max_concurrent_jobs=settings["max_in_flight"] if sharded else None,
            placement=settings["placement"] if sharded else None,
            note="Counts remain unknown until generators finish; native sampling and GPU microbatch settings remain on each tile."))
    return plan
