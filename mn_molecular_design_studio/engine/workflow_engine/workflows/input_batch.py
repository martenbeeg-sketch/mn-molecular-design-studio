"""Run repeated selections under one scheduler allocation and Results entry."""
from __future__ import annotations

import copy
import inspect
from pathlib import Path

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    JobPaths, adopt_next_job_creation, create_job, finish_job, read_json, update_status, write_json,
)


def _compound_batch(items, run_dir):
    """Send compatible compound selections to the engine in one container."""
    from mn_molecular_design_studio.ligand_compounds import load_compound_set, write_compound_set
    ignored = {"compound_set_file", "compound_limit", "job_name"}
    common = {key: value for key, value in items[0].items() if key not in ignored}
    if not all(item.get("compound_set_file") and common == {key: value for key, value in item.items() if key not in ignored}
               for item in items):
        return None
    rows, seen = [], set()
    for index, item in enumerate(items, 1):
        selected = load_compound_set(Path(item["compound_set_file"]))["compounds"]
        limit = item.get("compound_limit", 0)
        for row in selected[:limit] if limit else selected:
            row = copy.deepcopy(row)
            identifier = row["compound_id"]
            if identifier.casefold() in seen:
                row.setdefault("properties", {})["source_compound_id"] = identifier
                identifier = f"selection_{index}_{len(rows)}_{identifier[:65]}"
            while identifier.casefold() in seen:
                identifier = "copy_" + identifier
            seen.add(identifier.casefold())
            row["compound_id"] = identifier
            rows.append(row)
    path = run_dir / "artifacts/input/selected_compounds.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import utc_now
    write_compound_set(path, name="Selected compounds", compounds=rows,
                       owner_id="", project_id="", created_at=utc_now())
    return {**items[0], "compound_set_file": path, "compound_limit": 0}


def run_input_batch(workflow_id: str, items: list[dict], cpu_cores: int = 4,
                    gpu_device: str = "none", job_name: str = "Selected inputs") -> Path:
    """Prefer native collections; contain scalar-only adapters in one batch run.

    Scalar adapters share the parent's allocation and run in order. They never
    submit independent workers or reserve each item's resources simultaneously.
    """
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function, _jsonable
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
    spec = WORKFLOWS[workflow_id]
    if spec.mode != "queued":
        raise ValueError(f"{workflow_id} must use its native collection input.")
    function = _get_function(spec)
    job = create_job(spec.task_group, spec.job_type, spec.tool, {"items": _jsonable(items)},
                     {"job_name": job_name, "cpu_cores": cpu_cores, "gpu_device": gpu_device})
    if items and items[0].get("compound_set_file"):
        merged = _compound_batch(items, job.run_dir)
        if merged is not None:
            if "gpu_device" in inspect.signature(function).parameters:
                merged["gpu_device"] = gpu_device
            with adopt_next_job_creation(job):
                return Path(function(**merged))
    parent = read_json(job.run_dir / "metadata.json")
    rows = [{"id": f"input_{index:06d}", "name": item.get("job_name") or f"input_{index:06d}",
             "status": "pending", "run_dir": str(job.run_dir / "artifacts/batch_items" / f"{job.run_id}_input_{index:06d}")}
            for index, item in enumerate(items, 1)]
    candidates = []
    manifest = {"schema_version": 1, "automatic": True, "execution": "managed",
                "workflow": workflow_id, "item_count": len(items), "items": rows}
    signature = inspect.signature(function).parameters
    for index, parameters in enumerate(items, 1):
        if read_json(job.run_dir / "metadata.json").get("status") in {"cancelled", "paused"}:
            return job.run_dir
        item_id = f"input_{index:06d}"
        child = job.run_dir / "artifacts/batch_items" / f"{job.run_id}_{item_id}"
        row = rows[index - 1]
        row["status"] = "running"
        write_json(job.run_dir / "batch_manifest.json", manifest)
        update_status(job.run_dir, "running", batch={"automatic": True, "item_count": len(items),
                      "completed": sum(r["status"] == "completed" for r in rows)},
                      current_phase=f"Input {index} of {len(items)}: {row['name']}")
        previous = read_json(child / "result.json")
        if previous.get("success") is not True:
            child.mkdir(parents=True, exist_ok=True)
            (child / "artifacts").mkdir(exist_ok=True)
            for name in ("stdout.log", "stderr.log"):
                (child / name).touch()
            inherited = {key: parent[key] for key in ("scheduler_backend", "scheduler_job_id", "resource_allocation", "resource_request") if key in parent}
            write_json(child / "metadata.json", {**inherited, "run_id": child.name, "task_group": spec.task_group,
                       "status": "preparing", "hidden": True, "batch_parent_run_id": job.run_id})
            call = copy.deepcopy(parameters)
            if "gpu_device" in signature:
                call["gpu_device"] = gpu_device
            for key in ("cpu_cores", "cpu_threads"):
                if key in signature:
                    call[key] = min(int(call.get(key) or cpu_cores), cpu_cores)
            try:
                with adopt_next_job_creation(JobPaths(spec.task_group, child.name, child)):
                    function(**call)
            except Exception as exc:
                finish_job(child, False, {"error": str(exc)})
        stopped = read_json(job.run_dir / "metadata.json").get("status")
        if stopped in {"cancelled", "paused"}:
            for unfinished in rows[index - 1:]:
                unfinished["status"] = stopped
            write_json(job.run_dir / "batch_manifest.json", manifest)
            return job.run_dir
        result = read_json(child / "result.json")
        row.update(status="completed" if result.get("success") is True else "failed",
                   error=str(result.get("error") or result.get("metrics", {}).get("error") or ""),
                   metrics=result.get("metrics", {}))
        for candidate in read_candidates(child):
            candidate["candidate_id"] = f"{item_id}_{candidate['candidate_id']}"
            candidate.setdefault("raw_metadata", {})["batch_item_id"] = item_id
            candidates.append(candidate)
        write_json(job.run_dir / "batch_manifest.json", manifest)
    write_candidates(job.run_dir, spec.tool, candidates)
    failed = sum(row["status"] == "failed" for row in rows)
    metadata = read_json(job.run_dir / "metadata.json")
    metadata["batch"] = {"automatic": True, "execution": "managed", "item_count": len(rows),
                         "completed": len(rows) - failed, "failed": failed}
    write_json(job.run_dir / "metadata.json", metadata)
    finish_job(job.run_dir, bool(rows) and failed < len(rows), {
        "outputs": {"batch_items": rows},
        **({"warning": f"{failed} of {len(rows)} selected input calculations failed; successful results are retained."} if failed else {}),
        "metrics": {"item_count": len(rows), "completed_items": len(rows) - failed,
                    "failed_items": failed, "candidate_count": len(candidates)},
        "downstream_artifacts": {"batch_manifest": "batch_manifest.json",
                                 "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl"},
    })
    return job.run_dir
