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


def compatible_docking_pose_analysis(kind):
    """Raw ligand poses are inputs only to registered ligand-analysis tools."""
    from mn_molecular_design_studio.quality_tools import TOOLS
    return bool(TOOLS.get(kind, {}).get("kind") == "ligand" and compatible(kind))


def deferred_source(row, source):
    from mn_molecular_design_studio.cluster import config
    native_dir = str(row.get("run_dir") or "")
    if row.get("node_id") != config.identity():
        native_dir = config.artifact_uri(row["node_id"], row["task_group"], row["run_id"])
    selected_pose_ids = list(row.get("selected_pose_ids") or [])
    return dict(run_dir=native_dir,
                candidates_jsonl=None if selected_pose_ids else native_dir.rstrip("/") + "/artifacts/normalized_candidates/candidates.jsonl",
                deferred_run_source=True, target={}, require_native_passing=source.get("require_native_passing", False),
                selected_pose_ids=selected_pose_ids,
                source_job_reference=str(row.get("source_job_reference") or ""),
                source_workstation_id=str(row.get("node_id") or ""),
                source_workstation=str(row.get("node_name") or ""))


def source_request(kind, params, source, name):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
    if source.get("selected_pose_ids") and not compatible_docking_pose_analysis(kind):
        raise ValueError("Selected raw docking poses can connect only to a registered protein–ligand analysis such as PLIP.")
    if kind == "shared_filter":
        from mn_molecular_design_studio.validation_acceptance import campaign_policy_from_node
        return dict(schema_version=1, workflow="analysis.validation_acceptance", parameters={
            "source_run_dir": source["run_dir"], "candidates_jsonl": source["candidates_jsonl"],
            "acceptance_policy": campaign_policy_from_node(params), "job_name": name})
    if source.get("selected_pose_ids"):
        params = copy.deepcopy(params)
        params["candidate_ids"] = list(source["selected_pose_ids"])
        params["source_job_reference"] = source.get("source_job_reference") or None
        params["source_workstation_id"] = source.get("source_workstation_id") or None
        params["source_workstation"] = source.get("source_workstation") or None
        source = {**source, "candidates_jsonl": None}
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


def load_source(reference, *, owner_id, project_id, selected_docking_poses=None, allow_development=False):
    from mn_molecular_design_studio import logical_runs
    from mn_molecular_design_studio.cluster.candidate_batches import actor
    from mn_molecular_design_studio.cluster import config
    user = actor(owner_id, project_id, allow_development=allow_development)
    run, output = logical_runs.resolve(reference, logical_runs.catalog(user, project_id=project_id))
    if selected_docking_poses is not None:
        if (not isinstance(selected_docking_poses, list) or not selected_docking_poses
                or len(selected_docking_poses) > 1000):
            raise ValueError("Select one or more completed docking poses from the saved Run output.")
        selected_by_source = {}
        for item in selected_docking_poses:
            if not isinstance(item, dict) or set(item) != {"source_job_reference", "pose_id"}:
                raise ValueError("Docking selections must identify an authorized source job and exact pose.")
            source_reference, pose_id = str(item["source_job_reference"] or ""), str(item["pose_id"] or "")
            if not source_reference or not pose_id or pose_id in selected_by_source.setdefault(source_reference, []):
                raise ValueError("Docking selections contain an empty or duplicate source pose identity.")
            selected_by_source[source_reference].append(pose_id)
        sources = []
        for row in output["sources"]:
            source_reference = logical_runs.reference(row)
            pose_ids = selected_by_source.get(source_reference)
            if not pose_ids:
                continue
            if str(row.get("status") or "").lower() not in {"completed", "succeeded"} or not row.get("has_docking_results"):
                raise ValueError("A selected source job is not a completed docking result in this Run output.")
            if str(row.get("project_id") or "") != str(project_id):
                raise PermissionError("The selected docking source belongs to a different project.")
            local = row.get("node_id") == config.identity()
            if local:
                from mn_molecular_design_studio.runtime import runs_root
                source_run = Path(str(row.get("run_dir") or "")).expanduser().resolve(strict=True)
                root = Path(runs_root()).resolve()
                if not source_run.is_relative_to(root):
                    raise ValueError("The selected docking source path is outside the authorized Results store.")
                if source_run.name != str(row.get("run_id")) or source_run.parent.name != str(row.get("task_group")):
                    raise ValueError("The selected docking source path does not match its authorized job identity.")
                from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_pose_options
                options = docking_pose_options(source_run, source_job_reference=source_reference,
                                               source_node_id=str(row.get("node_id") or ""))
                available = {option["pose_id"] for option in options if option.get("available")}
                missing = [pose_id for pose_id in pose_ids if pose_id not in available]
                if missing:
                    raise ValueError("A selected docking pose is unavailable or lacks recorded ligand chemistry.")
            row["source_job_reference"] = source_reference
            row["selected_pose_ids"] = list(pose_ids)
            row["source_candidate_count"] = len(pose_ids)
            sources.append(row)
        if set(selected_by_source) != {logical_runs.reference(row) for row in sources}:
            raise ValueError("One or more selected docking source jobs are not part of this authorized Run output.")
    else:
        sources = [r for r in output["sources"] if r.get("status") in logical_runs.TERMINAL
                   and (r.get("has_candidate_collection") or int((r.get("candidate_summary") or {}).get("candidate_count") or 0) > 0)]
    if not sources:
        if selected_docking_poses is None and any(r.get("has_docking_results") for r in output["sources"]):
            raise ValueError("This Run output contains raw docking results; select one or more completed poses before connecting it to PLIP.")
        raise ValueError("This Run output has no completed saved candidate collections or selected docking poses.")
    for row in sources:
        if row.get("selected_pose_ids"):
            continue
        if row.get("node_id") == config.identity():
            path = Path(row["run_dir"]) / "artifacts/normalized_candidates/candidates.jsonl"
            if path.is_file():
                stat = path.stat()
                counts, _ = logical_runs._collection_properties(str(path), stat.st_mtime_ns, stat.st_size)
                row["source_candidate_count"] = counts["total"]
    return dict(run_sources=sources, run_reference=output["reference"], target={})
