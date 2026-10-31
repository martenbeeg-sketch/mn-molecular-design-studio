"""Finite whole-Run operations submitted through the existing workflow API.

No scheduler/controller runs here. The existing queue owns each native job;
this manifest freezes input membership and records idempotent dispatch receipts.
"""
from __future__ import annotations

import copy
import hashlib
import inspect
import json
import uuid

from mn_molecular_design_studio import logical_runs as runs
from mn_molecular_design_studio.cluster.config import atomic_json


def prepare(payload, catalog):
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function
    run, output = runs.resolve(payload["source_job"], catalog)
    spec = WORKFLOWS.get(payload.get("workflow"))
    if spec is None:
        raise ValueError("Unknown workflow.")
    signature = inspect.signature(_get_function(spec)).parameters
    if not {"source_run_dir", "candidates_jsonl"}.issubset(signature) and spec.workflow_id != "campaign.validation_sequence":
        raise ValueError("This tool does not accept existing result collections. Use a compatible workflow input.")
    parameters = payload.get("parameters", {})
    if not isinstance(parameters, dict):
        raise ValueError("A Run operation needs one shared parameters object.")
    if {"source_run_dir", "candidates_jsonl", "initial_source"}.intersection(parameters):
        raise ValueError("Choose the Run output or explicit source paths, not both.")
    selected = []
    excluded = []
    for source in output["sources"]:
        has_results = source.get("has_candidate_collection")
        if has_results is None:
            count = (source.get("candidate_summary") or {}).get("candidate_count")
            has_results = count is None or int(count) > 0
        if str(source.get("status") or "").lower() not in runs.TERMINAL or not has_results:
            excluded.append({"source_job": runs.source_reference(source), "status": source.get("status"),
                             "reason": "No finished saved candidate collection"})
        else:
            selected.append(source)
    if not selected:
        raise ValueError("This Run output has no completed candidate collections to process.")
    from mn_molecular_design_studio.quality_tools import WORKFLOWS as quality_workflows
    if len(selected) > 1 and quality_workflows.get(spec.workflow_id, {}).get("cohort"):
        raise ValueError("This analysis compares the complete cohort and cannot run independently per source batch. Connect a compatible analysis step and merge its candidates before this cohort operation.")
    # Candidate IDs are local to a native collection. A plain candidate-ID
    # selection spanning multiple batches would select unrelated duplicates.
    if len(selected) > 1 and any(parameters.get(key) for key in ("candidate_ids", "selected_candidate_ids", "validation_candidate_ids", "candidate_limit", "validation_candidate_limit")):
        raise ValueError("Select the complete Run output, or a saved subset with source-qualified candidate identities.")
    return run, output, selected, excluded


def requests(payload, sources, context=None, operation_id=None):
    result = []
    for index, row in enumerate(sources):
        item = copy.deepcopy(payload)
        item["source_job"] = runs.source_reference(row)
        if context:
            item["run_context"] = context
        execution = item.setdefault("execution", {})
        if execution.get("placement") not in {None, "auto", "local", row.get("node_id")}:
            raise ValueError("Whole-Run processing uses each source workstation. Choose automatic placement.")
        # source_job already supplies affinity. Native batching is preserved.
        from mn_molecular_design_studio.cluster import config
        execution["placement"] = "auto" if config.load()["mode"] != "disabled" else "local"
        if operation_id:
            execution["idempotency_key"] = operation_id + "-" + str(index)
        result.append(item)
    return result


def submit(payload, catalog, *, owner_id, validate, dispatch, idempotency_key=None):
    import fcntl
    from mn_molecular_design_studio.workspace import _studio_home
    key = idempotency_key or str(uuid.uuid4())
    lock_id = hashlib.sha256(f"{owner_id}:{payload.get('source_job')}:{key}".encode()).hexdigest()
    directory = _studio_home() / "run_operations"
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / (lock_id + ".lock")).open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        return _submit(payload, catalog, owner_id=owner_id, validate=validate, dispatch=dispatch, idempotency_key=key)


def _submit(payload, catalog, *, owner_id, validate, dispatch, idempotency_key):
    from mn_molecular_design_studio.workspace import _studio_home
    run, output, sources, excluded = prepare(payload, catalog)
    request_digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    operation = "op-" + hashlib.sha256(f"{owner_id}:{run['id']}:{idempotency_key or uuid.uuid4()}".encode()).hexdigest()[:24]
    path = _studio_home() / "run_operations" / (operation + ".json")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        state = json.loads(path.read_text())
        if state["request_digest"] != request_digest:
            raise ValueError("This Run operation key was already used with different settings.")
        # Recovery uses exactly the original source snapshot.
        children = state["requests"]
    else:
        context = dict(studio_run_id=run["id"], studio_run_name=run["name"][:200], studio_run_output_id=operation,
                       studio_run_output_name=str(payload.get("parameters", {}).get("job_name") or payload["workflow"])[:200])
        children = requests(payload, sources, context, operation)
        # Validate all sources before submitting any scientific work.
        for child in children:
            validate(child)
        state = dict(schema_version=1, run_id=run["id"], source_output=output["reference"], output_id=operation,
                     owner_id=owner_id, project_id=run["project_id"], request_digest=request_digest,
                     requests=children, payload=payload, idempotency_key=idempotency_key,
                     excluded=excluded, receipts={}, errors={})
        atomic_json(path, state)
    for index, child in enumerate(children):
        key = str(index)
        if key in state["receipts"]:
            continue
        try:
            state["receipts"][key] = dispatch(child)
            state["errors"].pop(key, None)
        except Exception as exc:
            # Preserve successful receipts; replaying with the same key safely
            # retries dispatch failures through existing API idempotency.
            state["errors"][key] = str(exc)
        atomic_json(path, state)
    receipts = list(state["receipts"].values())
    if not receipts:
        raise ValueError("No Run jobs were submitted: " + "; ".join(state["errors"].values()))
    return {"job": receipts[0]["job"], "jobs": [r["job"] for r in receipts], "run_id": run["id"],
            "operation_id": operation,
            "run_output": run["id"] + "~" + operation, "source_output": output["reference"],
            "submitted_job_count": len(receipts), "excluded_sources": state["excluded"],
            "dispatch_errors": state["errors"], "status": "partial" if state["errors"] else "submitted"}
