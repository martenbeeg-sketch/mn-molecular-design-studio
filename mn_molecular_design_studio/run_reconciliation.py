"""Recover historical Run membership from exact execution-ledger receipts.

This is an explicit metadata repair, never a name-based ownership migration.
The saved projection is small; Results need not re-read campaign ledgers.
"""
from pathlib import Path

from mn_molecular_design_studio import logical_runs as runs


def ledger_memberships(parent, ledger):
    if str(parent.get("distributed_campaign_id")) != str(ledger.get("campaign_id")):
        raise ValueError("The ledger does not belong to this campaign.")
    scope = dict(owner_id=str(parent.get("owner_id") or ""), project_id=str(parent.get("project_id") or ""))
    if not all(scope.values()) or not parent.get("node_id"):
        raise ValueError("A ledger repair requires an owned project and organizer identity.")
    tags = dict(parent_run_id=parent["run_id"], parent_node_id=parent["node_id"])
    result = {}
    for attempt in ledger.get("attempts", {}).values():
        peer = attempt.get("cluster_submission") or {}
        job_id = attempt.get("cluster_job_id") or peer.get("cluster_job_id")
        if job_id:
            result["cluster:" + str(job_id)] = dict(**scope, tags=tags, evidence=ledger["campaign_id"])
        child = Path(attempt.get("child_run_dir") or "")
        if child.name and str(child) != ".":
            key = f"{parent['node_id']}:{child.parent.name}/{child.name}"
            result[key] = dict(**scope, tags={**tags, "execution_proxy": bool(attempt.get("cluster_remote") and job_id)},
                               evidence=ledger["campaign_id"])
    return result


def merge_memberships(existing, incoming):
    result = dict(existing)
    for key, item in incoming.items():
        previous = result.get(key)
        if previous and previous != item:
            raise ValueError("Conflicting Run membership for " + key)
        result[key] = item
    return result


def analysis_identity_patch(run_dir):
    """Recover absent analysis identity only from its recorded local input."""
    import json
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    from mn_molecular_design_studio.runtime import runs_root
    run_dir = Path(run_dir)
    metadata = json.loads((run_dir / "metadata.json").read_text())
    if metadata.get("task_group", run_dir.parent.name) != "analysis":
        return {}
    inputs = json.loads((run_dir / "input.json").read_text())
    source = resolve_stored_path((inputs.get("inputs") or {}).get("source_run_dir"), run_dir=run_dir, must_exist=True)
    if source is None or not source.resolve().is_relative_to(runs_root().resolve()):
        return {}
    source_metadata = json.loads((source / "metadata.json").read_text())
    marker = source / "molecular_design_studio.json"
    identity = {**(json.loads(marker.read_text()) if marker.is_file() else {}), **source_metadata}
    if not identity.get("studio_owner_id") or not (identity.get("project_id") or identity.get("studio_project_id")):
        return {}
    patch = {}
    for field, value in {"studio_owner_id": identity["studio_owner_id"],
                         "project_id": identity.get("project_id") or identity.get("studio_project_id")}.items():
        if metadata.get(field) and metadata[field] != value:
            raise ValueError("Analysis ownership differs from its recorded source: " + str(run_dir))
        if not metadata.get(field):
            patch[field] = value
    parent_id = source_metadata.get("campaign_id") or source.name
    if not metadata.get("campaign_id"):
        patch["campaign_id"] = parent_id
    if not metadata.get("parent_node_id"):
        patch["parent_node_id"] = config.identity()
    name = (inputs.get("params") or {}).get("job_name")
    if name and metadata.get("campaign_name") != name:
        patch["campaign_name"] = name
    return patch
