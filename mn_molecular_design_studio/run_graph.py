"""Run-output inputs for the existing graph controller and candidate ledger."""
import copy
from pathlib import Path


def compatible(kind):
    import inspect
    from mn_molecular_design_studio.quality_tools import TOOLS
    if TOOLS.get(kind, {}).get("cohort"):
        return False
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function
    if kind == "shared_filter":
        return True
    spec = WORKFLOWS.get(NODE_CATALOG[kind].get("workflow"))
    return bool(spec and {"source_run_dir", "candidates_jsonl"}.issubset(inspect.signature(_get_function(spec)).parameters))


def deferred_source(row, source):
    from mn_molecular_design_studio.cluster import config
    native_dir = str(row.get("run_dir") or "")
    if row.get("node_id") != config.identity():
        native_dir = config.artifact_uri(row["node_id"], row["task_group"], row["run_id"])
    return dict(run_dir=native_dir, candidates_jsonl=native_dir.rstrip("/") + "/artifacts/normalized_candidates/candidates.jsonl",
                deferred_run_source=True, target={}, require_native_passing=source.get("require_native_passing", False))


def source_request(kind, params, source, name):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
    if kind == "shared_filter":
        from mn_molecular_design_studio.validation_acceptance import campaign_policy_from_node
        return dict(schema_version=1, workflow="analysis.validation_acceptance", parameters={
            "source_run_dir": source["run_dir"], "candidates_jsonl": source["candidates_jsonl"],
            "acceptance_policy": campaign_policy_from_node(params), "job_name": name})
    if not NODE_CATALOG[kind].get("workflow"):
        raise ValueError("Connect Run results to a compatible analysis, sequence-design or prediction tool before cohort operations.")
    return node_request(kind, params, source, name)


def initialize_batches(parent, node, source, params, key, settings):
    from mn_molecular_design_studio.candidate_execution import digest, data_hashes
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.cluster.config import atomic_json
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function
    import inspect
    destination = Path(parent) / "candidate_batches" / node["id"] / key
    destination.mkdir(parents=True, exist_ok=True)
    batches = []
    offset = 0
    for index, row in enumerate(source["run_sources"]):
        local = row.get("node_id") == config.identity()
        deferred = deferred_source(row, source)
        request = source_request(node["data"]["kind"], copy.deepcopy(params), deferred, node["data"]["label"])
        spec = WORKFLOWS[request["workflow"]]
        if not {"source_run_dir", "candidates_jsonl"}.issubset(inspect.signature(_get_function(spec)).parameters):
            raise ValueError("This tool has no compatible saved-collection input contract.")
        batch_id = f"source_{index:06d}"
        count = row.get("source_candidate_count") or (row.get("candidate_summary") or {}).get("candidate_count") or 1
        manifest = dict(schema_version=1, pipeline=str(parent), node=node["id"], batch_id=batch_id,
            scientific_key=key, offset=offset, input_candidate_ids=[], unit_count=count, request=request,
            input_hashes=data_hashes(request) if local else {}, policy=settings,
            source_node_id=row.get("node_id"), source_run_id=row["run_id"], unit_kind="source_collection")
        manifest["sha256"] = digest(manifest)
        path = destination / (batch_id + ".json")
        if path.is_file():
            import json
            if json.loads(path.read_text()) != manifest:
                raise ValueError("The frozen Run source collection changed; submit a new graph Run.")
        else:
            atomic_json(path, manifest)
        batches.append(dict(batch_id=batch_id, offset=offset, unit_count=count, manifest=str(path),
                            manifest_sha256=manifest["sha256"], status="pending", attempts=[]))
        offset += count
    return batches


def load_source(reference, *, owner_id, project_id, allow_development=False):
    from mn_molecular_design_studio import logical_runs
    from mn_molecular_design_studio.cluster.candidate_batches import actor
    user = actor(owner_id, project_id, allow_development=allow_development)
    run, output = logical_runs.resolve(reference, logical_runs.catalog(user, project_id=project_id))
    sources = [r for r in output["sources"] if r.get("status") in logical_runs.TERMINAL
               and (r.get("has_candidate_collection") or int((r.get("candidate_summary") or {}).get("candidate_count") or 0) > 0)]
    if not sources:
        raise ValueError("This Run output has no completed saved candidate collections.")
    from mn_molecular_design_studio.cluster import config
    for row in sources:
        if row.get("node_id") == config.identity():
            path = Path(row["run_dir"]) / "artifacts/normalized_candidates/candidates.jsonl"
            stat = path.stat()
            counts, _ = logical_runs._collection_properties(str(path), stat.st_mtime_ns, stat.st_size)
            row["source_candidate_count"] = counts["total"]
    return dict(run_sources=sources, run_reference=output["reference"], target={})
