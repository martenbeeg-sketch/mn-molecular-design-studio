"""Scientific Runs and their outputs, independent of scheduler job boundaries.

A Run is a projection over durable job identities. Batches and analysis jobs keep
their existing directories, ownership, scheduling and recovery records.
"""
from __future__ import annotations

import hashlib
import re
import json
from functools import lru_cache
from pathlib import Path

TERMINAL = {"completed", "succeeded", "failed", "lost", "cancelled", "canceled", "stopped", "partial"}
SUCCESS = {"completed", "succeeded"}
IDENTITY_FIELDS = ("studio_run_id", "studio_run_name", "studio_run_output_id", "studio_run_output_name",
                   "parent_run_id", "parent_node_id", "parent_role", "distributed_campaign_id",
                   "campaign_group_id", "campaign_group_name", "has_candidate_collection", "has_docking_results",
                   "cluster_job_id", "organizer_node_id", "graph_node_id", "execution_proxy")


def reference(row):
    if row.get("task_group") and row.get("run_id"):
        return f"{row.get('node_id', '')}:{row['task_group']}/{row['run_id']}"
    return str(row.get("cluster_job_id") or "")


def source_reference(row):
    from mn_molecular_design_studio.cluster import config
    local = str(row.get("local_reference") or f"{row.get('task_group', '')}/{row.get('run_id', '')}")
    if row.get("task_group") == "pending" and row.get("cluster_job_id"):
        return "cluster:" + str(row["cluster_job_id"])
    if not row.get("node_id") or row["node_id"] == config.identity():
        return local
    return "cluster:" + str(row.get("cluster_job_id") or config.job_id(row["node_id"], local))


def stable_id(*values):
    return "run:" + hashlib.sha256(repr(values).encode()).hexdigest()[:24]


def legacy_batch(row):
    """Compatibility with the batch names already recognized by Results."""
    name = str(row.get("campaign_name") or "")
    match = re.search(r"(?i)\bbatch\s*#?(\d+)\b", name)
    if not match or not name[:match.start()].strip(" ·|-"):
        return None
    label = name[:match.start()].strip(" ·|-")
    date = re.search(r"\b20\d{2}-\d{2}-\d{2}\b", name[match.end():])
    return {"label": label + (" · " + date[0] if date else ""), "batch_number": int(match[1])}


def status(rows):
    states = [str(r.get("status") or "unknown").lower() for r in rows]
    if any(s not in TERMINAL for s in states):
        return "running"
    if all(s in SUCCESS for s in states):
        return "completed"
    if any(s in SUCCESS or s == "partial" for s in states):
        return "partial"
    return "failed" if any(s in {"failed", "lost"} for s in states) else "cancelled"


def group_runs(rows):
    """Resolve parents before filtering or pagination; never cross access scopes.

    Explicit parent/workstation identities take precedence. Historical name
    grouping is restricted to the same project, owner and producer tool.
    """
    normalized = []
    for original in rows:
        row = dict(original)
        if not row.get("run_id") and row.get("cluster_job_id"):
            row.update(run_id=row["cluster_job_id"], task_group=row.get("task_group") or "pending")
        if row.get("task_group") == "pending":
            row["has_candidate_collection"] = False
        if row.get("run_id") and not row.get("capacity_parent"):
            normalized.append(row)
    unique = {reference(row): row for row in normalized}
    rows = list(unique.values())
    scope = lambda r: (str(r.get("project_id") or ""), str(r.get("owner_id") or ""))
    aliases = {}
    for row in rows:
        for identifier in {row["run_id"], row.get("distributed_campaign_id")} - {None, ""}:
            aliases.setdefault((*scope(row), str(identifier)), []).append(row)
    parents = {}
    for row in rows:
        parent_id = row.get("parent_run_id") or row.get("campaign_id")
        if not parent_id or str(parent_id) == str(row["run_id"]):
            continue
        candidates = [r for r in aliases.get((*scope(row), str(parent_id)), []) if reference(r) != reference(row)]
        node = row.get("pipeline_parent_node_id") or row.get("parent_node_id") or row.get("organizer_node_id")
        if node:
            candidates = [r for r in candidates if r.get("node_id") == node]
        elif len(candidates) > 1:
            candidates = [r for r in candidates if r.get("node_id") == row.get("node_id")]
        if len(candidates) == 1:
            parents[reference(row)] = candidates[0]

    def root(row):
        seen = set()
        while reference(row) in parents:
            if reference(row) in seen:
                raise ValueError("Cyclic Run parent identities.")
            seen.add(reference(row))
            row = parents[reference(row)]
        return row

    groups = {}
    for row in rows:
        parent = root(row)
        row["_logical_root_ref"] = reference(parent)
        legacy = legacy_batch(parent)
        explicit = parent.get("campaign_group_id")
        rid = row.get("studio_run_id") or parent.get("studio_run_id")
        if not rid:
            rid = stable_id(*scope(parent), "group", explicit) if explicit else (
                stable_id(*scope(parent), "legacy", parent.get("tool"), legacy["label"]) if legacy else
                stable_id(*scope(parent), "job", reference(parent)))
        # Persisted links cannot escape the project and owner of a Run.
        key = (*scope(row), rid)
        item = groups.setdefault(key, {"id": rid, "root": parent, "jobs": [], "parented": False, "project_id": row.get("project_id"),
            "name": str(parent.get("studio_run_name") or parent.get("campaign_group_name") or
                        (legacy or {}).get("label") or parent.get("campaign_name") or parent.get("tool") or "Run")})
        item["jobs"].append(row)
        if reference(parent) != reference(row):
            item["parented"] = True
            item["root"] = parent
        elif not item["parented"] and not parent.get("studio_run_output_id"):
            item["root"] = parent
    result = []
    for item in groups.values():
        if all(r.get("hidden") for r in item["jobs"]):
            continue
        item["jobs"].sort(key=lambda r: (str(r.get("created_at") or ""), reference(r)))
        item["status"] = status(item["jobs"])
        item["job_count"] = len(item["jobs"])
        item["updated_at"] = max((str(r.get("updated_at") or r.get("created_at") or "") for r in item["jobs"]), default="")
        item["outputs"] = run_outputs(item)
        for row in item["jobs"]:
            row.pop("_logical_root_ref", None)
        result.append(item)
    return sorted(result, key=lambda r: r["updated_at"], reverse=True)


def run_outputs(run):
    groups = {}
    root = run["root"]
    for row in run["jobs"]:
        output_id = row.get("studio_run_output_id")
        output_name = row.get("studio_run_output_name")
        stage_id = row.get("graph_node_id") or (row.get("campaign_step_index") if row.get("pipeline_parent_node_id") else None)
        if not output_id and stage_id and not row.get("pipeline_iteration") and root.get("tool") == "visual_pipeline":
            stage_id = re.sub(r"_b\d+_a\d+$", "", str(stage_id))
            output_id = "step-" + hashlib.sha256(stage_id.encode()).hexdigest()[:16]
            output_name = row.get("campaign_step") or row.get("tool") or "Workflow output"
        # Historical analyses already record their upstream job. Keep each
        # launch separate; identical tool names do not imply identical policy.
        if not output_id and row.get("task_group") == "analysis" and reference(row) != reference(root):
            output_id = "analysis-" + hashlib.sha256(reference(row).encode()).hexdigest()[:16]
        output_id = output_id or "results"
        output = groups.setdefault(output_id, {"id": output_id,
            "name": output_name or ("Results" if output_id == "results" else row.get("campaign_name") or row.get("tool")),
            "jobs": []})
        output["jobs"].append(row)
    for output in groups.values():
        jobs = output["jobs"]
        # A parent's normalized collection is authoritative. Do not include
        # its batches again and thereby count each candidate twice.
        authoritative = [r for r in jobs if reference(r) == reference(root) and r.get("has_candidate_collection")
                         and len({r.get("_logical_root_ref", reference(r)) for r in jobs}) == 1]
        output["sources"] = authoritative or [r for r in jobs if not r.get("execution_proxy") and not any(
            child.get("parent_run_id") == r.get("run_id") and child.get("node_id") == r.get("node_id") for child in jobs)
            and not (r.get("parent_role") or "").endswith("_validation")]
        output["status"] = status(jobs)
        output["dispatch_error_count"] = max((int(r.get("run_dispatch_error_count") or 0) for r in jobs), default=0)
        output["expected_source_count"] = max((int(r.get("run_expected_source_count") or 0) for r in jobs), default=len(output["sources"]))
        if output["dispatch_error_count"]:
            output["status"] = "partial"
        output["reference"] = run["id"] + "~" + output["id"]
    return sorted(groups.values(), key=lambda o: (o["id"] != "results", o["id"]))


def catalog(user, *, project_id=None, token=""):
    """Use the existing indexed local and federated catalogs, with full paging."""
    from mn_molecular_design_studio.auth import accessible_project_ids
    from mn_molecular_design_studio.job_catalog import search_jobs
    from mn_molecular_design_studio.cluster import config, service
    all_rows = {}
    def collect(search):
        offset = 0
        while True:
            page = search(limit=1000, offset=offset)
            for row in page["items"]:
                row = dict(row)
                row["task_group"] = row.get("task_group") or "pending"
                row["run_id"] = row.get("run_id") or row.get("cluster_job_id")
                if not row["run_id"]:
                    continue
                row.setdefault("node_id", config.identity())
                row.setdefault("node_name", config.load().get("name"))
                all_rows[(row["node_id"], row["task_group"], row["run_id"])] = row
            offset += len(page["items"])
            if offset >= int(page.get("total", offset + len(page["items"]))) or not page["items"]:
                break
    if config.load()["mode"] != "disabled":
        collect(lambda **kw: service.shared_jobs(user, token=token, project_id=project_id, include_hidden=True, **kw))
    collect(lambda **kw: search_jobs(project_id=project_id, include_hidden=True, access_owner_id=user["id"],
        accessible_project_ids=accessible_project_ids(user), access_admin=user.get("role") == "admin", **kw))
    return group_runs(apply_operation_links(all_rows.values()))


def apply_operation_links(rows):
    """Organizer receipts also link jobs executed by older compatible workers."""
    import json
    from mn_molecular_design_studio.workspace import _studio_home
    links = {}
    for path in (_studio_home() / "run_operations").glob("op-*.json"):
        state = json.loads(path.read_text())
        for index, receipt in state.get("receipts", {}).items():
            job = receipt.get("job") or {}
            request = state["requests"][int(index)]
            context = request.get("run_context")
            if context:
                links[str(job.get("cluster_job_id") or job.get("run_id"))] = (context, state)
    membership_path = _studio_home() / "run_membership.json"
    memberships = {}
    if membership_path.is_file():
        stat = membership_path.stat()
        memberships = _saved_memberships(str(membership_path), stat.st_mtime_ns, stat.st_size)
    result = []
    for row in rows:
        row = dict(row)
        for key in ("cluster:" + str(row.get("cluster_job_id") or ""), reference(row)):
            saved = memberships.get(key)
            if saved and (str(row.get("owner_id") or ""), str(row.get("project_id") or "")) == (saved["owner_id"], saved["project_id"]):
                row.update(saved["tags"])
        link = links.get(str(row.get("cluster_job_id") or "")) or links.get(str(row.get("run_id") or ""))
        if link and str(row.get("project_id")) == str(link[1]["project_id"]) and str(row.get("owner_id")) == str(link[1]["owner_id"]):
            row.update(link[0])
            row.update(run_dispatch_error_count=len(link[1].get("errors", {})), run_expected_source_count=len(link[1]["requests"]))
        result.append(row)
    return result


@lru_cache(maxsize=8)
def _saved_memberships(path, mtime_ns, size):
    return json.loads(Path(path).read_text())["memberships"]


def resolve(reference_text, runs):
    run_id, _, output_id = str(reference_text).partition("~")
    matches = [run for run in runs if run["id"] == run_id]
    if len(matches) != 1:
        raise LookupError("Run is unavailable or ambiguous in this project.")
    run = matches[0]
    outputs = [o for o in run["outputs"] if o["id"] == (output_id or "results")]
    if len(outputs) != 1:
        raise LookupError("The selected Run output is unavailable.")
    return run, outputs[0]


def public(run):
    """Compact metadata only. Scientific files are fetched on explicit request."""
    return {**{k: v for k, v in run.items() if k not in {"jobs", "root", "outputs", "parented"}},
            "outputs": [{k: v for k, v in o.items() if k not in {"jobs", "sources"}} | {"source_count": len(o["sources"])} for o in run["outputs"]]}


def context_tags(value):
    if value is None:
        return {}
    keys = {"studio_run_id", "studio_run_name", "studio_run_output_id", "studio_run_output_name"}
    if not isinstance(value, dict) or set(value) != keys or any(not isinstance(v, str) or not v or len(v) > 200 for v in value.values()):
        raise ValueError("Run context requires a Run identity, name and output identity/name.")
    if not re.fullmatch(r"run:[a-f0-9]{24}", value["studio_run_id"]) or not re.fullmatch(r"[a-zA-Z0-9_-]+", value["studio_run_output_id"]):
        raise ValueError("Invalid Run or output identity.")
    return dict(value)


def native_context(run_dir):
    """Carry a native controller's logical identity to remote execution jobs."""
    from mn_molecular_design_studio.cluster import config
    metadata = json.loads((Path(run_dir) / "metadata.json").read_text())
    metadata.setdefault("run_id", Path(run_dir).name)
    metadata.setdefault("task_group", Path(run_dir).parent.name)
    row = dict(metadata, owner_id=metadata.get("studio_owner_id"),
               project_id=metadata.get("studio_project_id") or metadata.get("project_id"),
               node_id=config.identity(), hidden=False, capacity_parent=False)
    run = group_runs([row])[0]
    return context_tags(dict(studio_run_id=run["id"], studio_run_name=run["name"][:200],
        studio_run_output_id=metadata.get("studio_run_output_id") or "results",
        studio_run_output_name=metadata.get("studio_run_output_name") or "Results"))


@lru_cache(maxsize=1024)
def _collection_properties(path, mtime, size):
    counts = {"complex": 0, "sequence": 0, "structure": 0, "total": 0}
    stages = set()
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            counts["total"] += 1
            counts["complex"] += bool(row.get("complex_pdb"))
            counts["structure"] += bool(row.get("complex_pdb") or row.get("binder_pdb"))
            counts["sequence"] += bool(row.get("binder_sequence"))
            stages.add(row.get("stage"))
    return counts, stages


def input_options(project_id, user, *, token="", require_complex=True, allow_docking=False, stages=None):
    """Whole-output choices shared by analysis and sequence-design tiles."""
    from mn_molecular_design_studio.cluster import config
    options = []
    for run in catalog(user, project_id=project_id, token=token):
        for output in run["outputs"]:
            count = 0
            eligible = []
            for row in output["sources"]:
                if row.get("status") not in TERMINAL:
                    continue
                local = not row.get("node_id") or row["node_id"] == config.identity()
                if local:
                    path = Path(str(row.get("run_dir") or "")) / "artifacts/normalized_candidates/candidates.jsonl"
                    if path.is_file():
                        stat = path.stat()
                        counts, available_stages = _collection_properties(str(path), stat.st_mtime_ns, stat.st_size)
                        n = counts["complex"] if require_complex else max(counts["sequence"], counts["structure"])
                        if n and (not stages or set(stages).intersection(available_stages)):
                            count += n
                            eligible.append(row)
                    elif allow_docking and row.get("has_docking_results"):
                        eligible.append(row)
                elif row.get("has_candidate_collection") or int((row.get("candidate_summary") or {}).get("candidate_count") or 0) > 0:
                    eligible.append(row)
                    count += int((row.get("candidate_summary") or {}).get("candidate_count") or 0)
            if not eligible:
                continue
            label = run["name"] + " · " + str(output["name"])
            label += f" · {len(eligible)} source collection(s)"
            if count:
                label += f" · {count:,} candidates"
            options.append({"reference": output["reference"], "label": label,
                            "row": {**run["root"], "_run_sources": eligible, "_logical_run_id": run["id"],
                                    "_run_output_id": output["id"]}})
    return options
