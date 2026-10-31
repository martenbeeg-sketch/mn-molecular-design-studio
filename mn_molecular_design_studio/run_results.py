"""Combined, source-qualified candidate metadata; structures stay on their owner."""
import copy
import hashlib
import json
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


def collect(run, output, user, *, token=""):
    from mn_molecular_design_studio.cluster import config, service, transport
    from mn_molecular_design_studio.cluster.store import Store
    from mn_molecular_design_studio.workspace import _studio_home
    candidates, pending, errors, excluded = [], [], [], []
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
                key = hashlib.sha256(f"{user['id']}:{ref}:{source.get('updated_at')}".encode()).hexdigest()
                cache = _studio_home() / "run_result_receipts" / (key + ".json")
                cache.parent.mkdir(parents=True, exist_ok=True)
                if cache.is_file():
                    operation_id = json.loads(cache.read_text())["command_id"]
                else:
                    if source.get("availability") != "online":
                        raise ValueError("Source workstation is offline; its saved results are still part of this Run.")
                    if config.load()["mode"] in {"coordinator", "peer"}:
                        operation = Store().action(ref.removeprefix("cluster:"), "inspect", {}, service.user_record(user))
                    else:
                        operation = transport.request("/api/v1/cluster/jobs/" + ref.removeprefix("cluster:") + "/actions",
                            {"action": "inspect"}, method="POST", token=token)
                    config.atomic_json(cache, operation)
                    operation_id = operation["command_id"]
                operation = (Store().command(operation_id, service.user_record(user)) if config.load()["mode"] in {"coordinator", "peer"}
                             else transport.request("/api/v1/cluster/commands/" + operation_id, token=token))
                if operation["status"] == "failed":
                    raise ValueError(operation.get("error") or "Could not read the saved source results.")
                if operation["status"] != "completed":
                    pending.append({"source_job": ref, "command_id": operation_id})
                    continue
                rows = operation["response"].get("results", {}).get("candidates", [])
            else:
                root = Path(source["run_dir"])
                stat = (root / "artifacts/normalized_candidates/candidates.jsonl").stat()
                rows = _local_rows(str(root), stat.st_mtime_ns, stat.st_size)
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
        except (OSError, ValueError, LookupError) as exc:
            errors.append({"source_job": ref, "error": str(exc)})
    return {"run_id": run["id"], "output": output["reference"], "candidates": candidates,
            "candidate_count": len(candidates), "complete": not pending and not errors and not any(
                item.get("status") not in logical_runs.TERMINAL for item in excluded),
            "pending_sources": pending, "source_errors": errors, "excluded_sources": excluded}
