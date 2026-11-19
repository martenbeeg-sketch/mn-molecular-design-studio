"""Metadata-only project inventory and non-destructive cleanup plans."""

from datetime import datetime, timezone

from mn_molecular_design_studio.logical_runs import apply_operation_links, group_runs


TERMINAL_STATUSES = frozenset({"completed", "succeeded", "failed", "cancelled", "canceled"})


def load_project_jobs(local_search, *, remote_catalog=None, local_node_id="", local_node_name=""):
    page = local_search(limit=1000, include_hidden=True)
    local = list(page["items"])
    stale = False
    while len(local) < page["total"]:
        more = local_search(limit=1000, offset=len(local), include_hidden=True)["items"]
        if not more:
            stale = True
            break
        local.extend(more)
    remote, remote_stale = remote_catalog() if remote_catalog else ({"items": [], "total": 0}, False)
    stale = stale or remote_stale or len(remote["items"]) < remote["total"]
    merged = {execution_identity(row): dict(row) for row in remote["items"]}
    for row in local:
        row = dict(row, node_id=local_node_id, node_name=local_node_name)
        identity = execution_identity(row)
        merged[identity] = merged.get(identity, {}) | row
    return {"jobs": apply_operation_links(merged.values()), "stale": stale}


def execution_identity(row):
    return (str(row.get("node_id") or ""), str(row.get("task_group") or ""),
            str(row.get("run_id") or row.get("cluster_job_id") or ""))


def project_inventory(projects, jobs, collections):
    inventory = {}
    for project in projects:
        project_id = str(project["id"])
        executions = {}
        for row in jobs:
            if str(row.get("project_id") or "") == project_id:
                executions[execution_identity(row)] = dict(row)
        rows = list(executions.values())
        groups = [dict(row) for row in collections if str(row.get("project_id") or "") == project_id]
        finished = [row for row in rows if str(row.get("status") or "").casefold() in TERMINAL_STATUSES]
        protected = [row for row in groups if (row.get("origin") or {}).get("type") != "job"]
        runs = group_runs([dict(row, hidden=False) for row in rows])
        inventory[project_id] = {
            "project": dict(project), "jobs": rows, "collections": groups,
            "runs": len(runs), "logical_runs": runs, "finished": finished,
            "blocked": [row for row in rows if str(row.get("status") or "").casefold() not in TERMINAL_STATUSES], "protected": protected,
            "updated_at": max([str(project.get("created_at") or ""), str(project.get("updated_at") or ""),
                               *[str(row.get("updated_at") or row.get("completed_at") or row.get("created_at") or "") for row in rows]]),
        }
    return inventory


def cleanup_plan(item, *, catalog_stale=True):
    def reference(row):
        return {key: row.get(key) for key in (
            "cluster_job_id", "node_id", "node_name", "task_group", "run_id",
            "project_id", "owner_id", "status", "campaign_name", "hidden",
        )}

    return {
        "schema": "studio-project-cleanup-review-v1", "destructive": False,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "catalog_complete": not catalog_stale,
        "project": {key: item["project"].get(key) for key in ("id", "name")},
        "finished_executions_to_review": [reference(row) for row in item["finished"]],
        "blocked_executions": [reference(row) for row in item["blocked"]],
        "protected_collections": [{key: row.get(key) for key in ("id", "name", "molecule_count")}
                                  for row in item["protected"]],
        "policy": {
            "automatic_selection_or_deletion": False,
            "preserve_imports_and_explicit_saved_selections": True,
            "preserve_project_metadata": True,
            "revalidate_permissions_and_finished_state_at_purge": True,
            "backup_required_before_purge": True,
            "folder_deletion_alone_does_not_remove_sql_or_peer_records": True,
        },
    }
