"""Combined, source-qualified candidate metadata; structures stay on their owner."""
import copy
import fcntl
import hashlib
import json
import math
import time
from contextlib import contextmanager, nullcontext
from pathlib import Path
from functools import lru_cache

from mn_molecular_design_studio import logical_runs


@lru_cache(maxsize=128)
def _local_rows(run_dir, mtime_ns, size):
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    rows = read_candidates(Path(run_dir))
    for row in rows:
        for key in ("complex_pdb", "target_pdb", "binder_pdb"):
            if row.get(key):
                path = resolve_stored_path(row[key], run_dir=Path(run_dir))
                if path is not None:
                    row[key] = str(path)
    return rows


@contextmanager
def _inspection_lock(receipt):
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with receipt.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def _remote_inspection(source, ref, user, token, *, refresh, retry_failed):
    from mn_molecular_design_studio.cluster import config, service, transport
    from mn_molecular_design_studio.cluster.store import Store
    from mn_molecular_design_studio.workspace import _studio_home

    authorized = service.lookup_job(ref, user, token=token)
    if (str(authorized.get("node_id") or "") != str(source.get("node_id") or "")
            or str(authorized.get("local_reference") or "") != str(
                source.get("local_reference") or f"{source['task_group']}/{source['run_id']}")):
        raise PermissionError("The saved Results source identity changed.")
    key = hashlib.sha256(f"{user['id']}:{ref}:{source.get('updated_at')}".encode()).hexdigest()
    cache = _studio_home() / "run_result_receipts" / (key + ".json")
    coordinator = config.load()["mode"] == "coordinator"

    def poll(command_id):
        if coordinator:
            return Store().command(command_id, service.user_record(user))
        return transport.request("/api/v1/cluster/commands/" + command_id, token=token)

    with (_inspection_lock(cache) if refresh else nullcontext()):
        try:
            receipt = json.loads(cache.read_text())
        except FileNotFoundError:
            receipt = None
        if receipt is not None:
            operation = {"command_id": receipt["command_id"], **poll(receipt["command_id"])}
            if operation["status"] not in {"failed", "cancelled"}:
                return operation
            remaining = max(0, math.ceil(2 - (time.time() - receipt.get("requested_at", 0))))
            if not refresh or not retry_failed or remaining or source.get("availability") != "online":
                return {**operation, "retry_after_seconds": remaining}
        if not refresh:
            return None
        if source.get("availability") != "online":
            raise ValueError("Source workstation is offline; its saved results are still part of this Run.")
        job_id = ref.removeprefix("cluster:")
        if coordinator:
            operation = Store().action(job_id, "inspect", {}, service.user_record(user))
        else:
            operation = transport.request("/api/v1/cluster/jobs/" + job_id + "/actions",
                                          {"action": "inspect"}, method="POST", token=token)
        command_id = operation["command_id"]
        config.atomic_json(cache, {"command_id": command_id, "requested_at": time.time()})
        return {"command_id": command_id, **poll(command_id)}


def collect(run, output, user, *, token="", refresh=False, retry_failed=False):
    from mn_molecular_design_studio.cluster import transport
    candidates, pending, errors, excluded, inspection_required, source_details = [], [], [], [], [], {}
    if output.get("dispatch_error_count"):
        errors.append({"source_job": output["reference"], "error": f"{output['dispatch_error_count']} source collection(s) have not been submitted; retry their dispatch before exporting a complete output."})
    for source in output["sources"]:
        ref = logical_runs.source_reference(source)
        if source.get("status") not in logical_runs.TERMINAL:
            excluded.append({"source_job": ref, "status": source.get("status")})
            continue
        if source.get("has_candidate_collection") is False:
            excluded.append({"source_job": ref, "status": source.get("status"), "reason": "No saved candidate output"})
            continue
        try:
            if ref.startswith("cluster:"):
                operation = _remote_inspection(source, ref, user, token, refresh=refresh, retry_failed=retry_failed)
                if operation is None:
                    inspection_required.append({"source_job": ref, "availability": source.get("availability")})
                    continue
                if operation["status"] in {"failed", "cancelled"}:
                    errors.append({"source_job": ref, "error": operation.get("error") or "Could not read the saved source results.",
                                   "retryable": source.get("availability") == "online",
                                   "retry_after_seconds": operation.get("retry_after_seconds", 0)})
                    continue
                if operation["status"] != "completed":
                    pending.append({"source_job": ref, "command_id": operation["command_id"]})
                    continue
                rows = operation["response"].get("results", {}).get("candidates")
                if not isinstance(rows, list):
                    raise ValueError("The source inspection did not return a candidate collection.")
                details = operation["response"]
            else:
                root = Path(source["run_dir"])
                stat = (root / "artifacts/normalized_candidates/candidates.jsonl").stat()
                rows = _local_rows(str(root), stat.st_mtime_ns, stat.st_size)
                details = {"job": source, "metadata": json.loads((root / "metadata.json").read_text())
                           if (root / "metadata.json").is_file() else {}}
            source_details[ref] = {"job": details.get("job") or source,
                                   "metadata": details.get("metadata") or {}}
            ids = [row.get("candidate_id") for row in rows]
            if not all(ids) or len(ids) != len(set(ids)):
                raise ValueError("Source collection contains missing or duplicate candidate identities.")
            for original in rows:
                row = copy.deepcopy(original)
                # Preserve native identity and every measured value. This
                # separate ID makes repeated design_0 names unambiguous.
                row["run_candidate_id"] = hashlib.sha256(f"{ref}:{row['candidate_id']}".encode()).hexdigest()[:24]
                row["run_id"] = run["id"]
                row["run_output_id"] = output["id"]
                row["source_job"] = ref
                row["source_node_id"] = source.get("node_id")
                from mn_molecular_design_studio.candidate_sets import json_safe
                candidates.append(json_safe(row))
        except (OSError, ValueError, LookupError, transport.ClusterUnavailable) as exc:
            errors.append({"source_job": ref, "error": str(exc)})
    return {"run_id": run["id"], "output": output["reference"], "candidates": candidates,
            "candidate_count": len(candidates), "complete": not inspection_required and not pending and not errors and not any(
                item.get("status") not in logical_runs.TERMINAL for item in excluded),
            "pending_sources": pending, "source_errors": errors, "excluded_sources": excluded,
            "inspection_required": inspection_required, "source_details": source_details}


def read_source(source, user, *, token=""):
    """Read one saved collection through the same authorized Results transport."""
    reference = logical_runs.source_reference(source)
    return collect({"id": reference}, {"id": "results", "reference": reference,
                                      "sources": [source]}, user, token=token, refresh=True)


def completed_sources(user, *, token="", runs=None):
    """Flatten the existing federated catalog without losing producer identity."""
    sources = {}
    for run in logical_runs.catalog(user, token=token) if runs is None else runs:
        for source in run["jobs"]:
            if (str(source.get("status") or "").casefold() not in logical_runs.SUCCESS
                    or source.get("success") is False or source.get("hidden")
                    or source.get("has_candidate_collection") is False):
                continue
            sources[logical_runs.reference(source)] = {**source, 'logical_run_id': run['id']}
    return sorted(sources.values(), key=lambda row: str(row.get("updated_at") or row.get("created_at") or ""), reverse=True)
