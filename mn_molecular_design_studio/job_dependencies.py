"""Qualified execution dependencies shared by browsing and cancellation."""

import json
from pathlib import Path
from urllib.parse import urlsplit


TERMINAL = frozenset({"completed", "succeeded", "failed", "lost", "cancelled", "canceled", "stopped", "partial"})
CANCELLED = frozenset({"cancelled", "canceled", "stopped"})
FIELDS = ("dependency_refs", "upstream_node_id", "parent_task_group", "batch_parent_run_id")
SOURCE_KEYS = frozenset({"source_job", "source_run_dir", "candidates_jsonl", "candidates_csv"})


def identity(row):
    if row.get("task_group") and row.get("run_id"):
        return str(row.get("node_id") or "") + ":" + row["task_group"] + "/" + row["run_id"]
    return "cluster:" + str(row.get("cluster_job_id") or "")


def scope(row):
    return str(row.get("project_id") or ""), str(row.get("owner_id") or "")


def cancellation_requested(run_dir, metadata=None):
    return bool((metadata or {}).get("cancel_requested_at") or (run_dir / "cancellation_request.json").is_file())


def projected_references(metadata, inputs, *, root=None):
    references = list(metadata.get("dependency_refs") or [])

    def add(value, node=""):
        if not isinstance(value, str) or not value or "\\" in value:
            return
        if value.startswith("cluster:") and not value.startswith("cluster://"):
            references.append({"cluster_job_id": value.removeprefix("cluster:")})
            return
        parsed = urlsplit(value)
        if parsed.scheme == "cluster":
            node, parts = parsed.netloc, parsed.path.strip("/").split("/")
        elif parsed.scheme == "runs":
            parts = (parsed.netloc + parsed.path).strip("/").split("/")
        elif parsed.scheme:
            return
        elif Path(value).is_absolute():
            if root is None:
                return
            try:
                parts = Path(value).relative_to(root).parts
            except ValueError:
                return
        else:
            parts = value.split("/")
        if len(parts) >= 2 and all(part not in {"", ".", ".."} for part in parts):
            references.append({"node_id": str(node), "task_group": parts[0], "run_id": parts[1]})

    for field in ("pipeline_parent_run_dir", "active_learning_parent_run", "parent_run_dir"):
        raw = metadata.get(field)
        node = metadata.get("pipeline_parent_node_id") or metadata.get("parent_node_id") or ""
        if raw and Path(str(raw)).is_absolute() and node:
            parent = Path(str(raw))
            add(parent.parent.name + "/" + parent.name, node)
        else:
            add(raw, node)

    def walk(value):
        if isinstance(value, dict):
            for key, child in value.items():
                if key in SOURCE_KEYS:
                    add(child)
                elif isinstance(child, (dict, list)):
                    walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)

    walk(inputs)
    return [json.loads(value) for value in sorted({json.dumps(value, sort_keys=True) for value in references if isinstance(value, dict)})]


def parents(rows):
    rows = list(rows)
    indexed = {}
    cluster_ids = {}
    for row in rows:
        for alias in {str(row.get("run_id") or ""), str(row.get("distributed_campaign_id") or "")} - {""}:
            indexed.setdefault((*scope(row), alias), []).append(row)
        if row.get("cluster_job_id"):
            cluster_ids[(*scope(row), row["cluster_job_id"])] = row
    result = {identity(row): set() for row in rows}
    for row in rows:
        refs = list(row.get("dependency_refs") or [])
        for prefix, run_id, node in (
            ("parent", row.get("parent_run_id") or row.get("campaign_id") or row.get("batch_parent_run_id"),
             row.get("pipeline_parent_node_id") or row.get("parent_node_id") or row.get("organizer_node_id")),
            ("upstream", row.get("upstream_run_id"), row.get("upstream_node_id") or row.get("node_id")),
        ):
            if run_id:
                refs.append({"run_id": str(run_id), "task_group": row.get(prefix + "_task_group"), "node_id": node,
                             "parent": prefix == "parent"})
        summary = row.get("request_summary")
        if isinstance(summary, dict):
            refs.extend(projected_references({}, summary))
        for ref in refs:
            if not isinstance(ref, dict):
                continue
            if ref.get("cluster_job_id"):
                candidates = [cluster_ids.get((*scope(row), ref["cluster_job_id"]))]
                candidates = [candidate for candidate in candidates if candidate]
            else:
                candidates = indexed.get((*scope(row), str(ref.get("run_id") or "")), [])
                node = ref.get("node_id") or (None if ref.get("parent") else row.get("dependency_origin_node") or row.get("node_id"))
                if node:
                    candidates = [candidate for candidate in candidates if candidate.get("node_id") == node]
                if ref.get("task_group"):
                    candidates = [candidate for candidate in candidates if candidate.get("task_group") == ref["task_group"]]
                if not node and ref.get("parent") and len(candidates) > 1:
                    candidates = [candidate for candidate in candidates if candidate.get("node_id") == row.get("node_id")]
            if len(candidates) == 1 and identity(candidates[0]) != identity(row):
                result[identity(row)].add(identity(candidates[0]))
    return result


def descendants(rows, roots):
    rows = list(rows)
    edges = parents(rows)
    selected = {identity(row) for row in roots}
    children = {}
    for child, sources in edges.items():
        for source in sources:
            children.setdefault(source, set()).add(child)
    queue = list(selected)
    while queue:
        source = queue.pop()
        for child in children.get(source, set()) - selected:
            selected.add(child)
            queue.append(child)
    return [row for row in rows if identity(row) in selected]


def local_row(run_dir, *, node_id, root):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json_raw
    metadata = read_json_raw(run_dir / "metadata.json")
    marker = read_json_raw(run_dir / "molecular_design_studio.json")
    row = {**marker, **metadata, "node_id": node_id, "run_id": run_dir.name,
           "task_group": run_dir.parent.name, "run_dir": str(run_dir)}
    row["owner_id"] = row.get("studio_owner_id") or row.get("owner_id")
    row["project_id"] = row.get("studio_project_id") or row.get("project_id")
    row["dependency_refs"] = projected_references(row, read_json_raw(run_dir / "input.json"), root=root)
    return row


def local_rows(root):
    from mn_molecular_design_studio.cluster import config
    node = config.identity()
    result = []
    for path in root.glob("*/*/metadata.json"):
        if path.parent.is_symlink() or path.parent.parent.is_symlink():
            continue
        result.append(local_row(path.parent, node_id=node, root=root))
    return result


def cancelled_dependency(run_dir):
    """Check durable local ancestors at dispatch and immediately before execution."""
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import runs_root
    root = runs_root().resolve()
    if not run_dir.is_relative_to(root):
        return None
    node = config.identity()
    initial = local_row(run_dir, node_id=node, root=root)
    queue, seen = [initial], {identity(initial)}
    while queue:
        current = queue.pop()
        refs = list(current.get("dependency_refs") or [])
        for prefix in ("parent", "upstream"):
            if current.get(prefix + "_run_id") and current.get(prefix + "_task_group"):
                refs.append({"run_id": current[prefix + "_run_id"], "task_group": current[prefix + "_task_group"],
                             "node_id": current.get(prefix + "_node_id") or current.get("pipeline_parent_node_id") or node})
        for ref in refs:
            if ref.get("node_id") not in {None, "", node}:
                continue
            group, run = str(ref.get("task_group") or ""), str(ref.get("run_id") or "")
            if not group or not run or any(value in {".", ".."} or "/" in value or "\\" in value for value in (group, run)):
                continue
            path = root / group / run
            if path.is_symlink() or path.parent.is_symlink() or not (path / "metadata.json").is_file():
                continue
            parent = local_row(path, node_id=node, root=root)
            if scope(parent) != scope(initial) or identity(parent) in seen:
                continue
            if str(parent.get("status") or "").lower() in CANCELLED or cancellation_requested(path, parent):
                return identity(parent)
            seen.add(identity(parent))
            queue.append(parent)
    return None


def guard_execution(run_dir):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, stop_job
    metadata = read_json(run_dir / "metadata.json")
    if str(metadata.get("status") or "").lower() in CANCELLED:
        return False
    if cancellation_requested(run_dir, metadata):
        stop_job(run_dir, reason=metadata.get("cancel_reason") or "Recorded cancellation request")
        return False
    upstream = cancelled_dependency(run_dir)
    if upstream:
        stop_job(run_dir, reason="Upstream job cancelled: " + upstream)
        return False
    return True


def request_cancellation(row, user, *, token=""):
    """Route through current access checks rather than trusting a UI snapshot."""
    from mn_molecular_design_studio.cluster import config, service, transport
    from mn_molecular_design_studio.cluster.store import Store
    settings = config.load()
    if settings["mode"] != "disabled" and row.get("cluster_job_id"):
        if settings["mode"] in {"coordinator", "peer"}:
            return Store().action(row["cluster_job_id"], "cancel", {}, service.user_record(user))
        return transport.request("/api/v1/cluster/jobs/" + row["cluster_job_id"] + "/actions", token=token,
                                 payload={"action": "cancel"}, method="POST")
    if row.get("node_id") and row["node_id"] != config.identity():
        raise ValueError("The owning workstation must be connected before cancellation can be requested.")
    from mn_molecular_design_studio.auth import require_job_access
    from mn_molecular_design_studio.job_catalog import get_job_access_record
    from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import resolve_job
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
    record = get_job_access_record(row["task_group"], row["run_id"])
    require_job_access(user, record)
    if scope(record or {}) != scope(row):
        raise ValueError("Job ownership changed; refresh the project before cancellation.")
    current = resolve_job(row["task_group"] + "/" + row["run_id"], include_hidden=True)
    summary = stop_job(Path(current["run_dir"]), reason="Cancelled from Projects")
    return {**summary, "downstream_count": max(0, len(summary["cancelled_runs"]) - 1)}
