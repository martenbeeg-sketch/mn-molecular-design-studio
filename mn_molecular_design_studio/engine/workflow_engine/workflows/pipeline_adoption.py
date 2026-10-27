"""Attach completed jobs to a display-only pipeline without rerunning them."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job,
    finish_job,
    read_json,
    runs_root,
    utc_now,
    write_json,
)


_ADOPTION_KEYS = (
    "campaign_name",
    "campaign_id",
    "campaign_step",
    "campaign_step_index",
    "pipeline_parent_run_dir",
    "upstream_task_group",
    "upstream_run_id",
    "upstream_job_code",
)


def _signature(project_id: str, owner_id: str, groups: list[dict[str, Any]]) -> str:
    identity = {
        "project_id": project_id,
        "owner_id": owner_id,
        "pipelines": [
            {"name": p["name"], "runs": [f"{s['task_group']}/{s['run_id']}" for s in p["steps"]]}
            for p in groups
        ],
    }
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()


def _graph_for(name: str, steps: list[dict[str, Any]]) -> dict[str, Any]:
    """Build an inspectable graph snapshot that describes, but cannot launch, adoption."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, SCHEMA, clean_graph

    kinds = ("generation_source", "abmpnn", "boltz2")
    nodes = []
    for index, (kind, step) in enumerate(zip(kinds, steps)):
        fields = NODE_CATALOG[kind]["fields"]
        params = {field["name"]: field["default"] for field in fields}
        if kind == "generation_source":
            params.update(source_run_dir=step["run_dir"], selected_candidate_ids=step["selected_candidate_ids"])
            label = step["stage_label"]
            description = "Existing completed generator run; adopted without regeneration."
        elif kind == "abmpnn":
            p = step["input"].get("params", {})
            params.update(
                sequences=int(p.get("num_seq_per_target") or 1),
                sampling_temp=float(p.get("sampling_temp") or 0.1),
                omit_aas=str(p.get("omit_aas") or "C"),
                seed=int(p.get("seed") or 0),
            )
            for key in ("redesign_cdr1", "redesign_cdr2", "redesign_cdr3"):
                params[key] = True
            label = "SAbR + AbMPNN"
            description = "SAbR annotation and AbMPNN sequence redesign are both recorded inside this existing run."
        else:
            p = step["input"].get("params", {})
            engine_options = p.get("engine_options") if isinstance(p.get("engine_options"), dict) else {}
            params.update(
                samples=int(p.get("samples") or 1),
                num_recycles=int(p.get("num_recycles") or 3),
                sampling_steps=int(engine_options.get("sampling_steps") or 200),
                seed_count=int(engine_options.get("seed_count") or 1),
                use_target_template=False,
                write_full_pae=bool(engine_options.get("write_full_pae", True)),
                generate_ipsae=bool(engine_options.get("generate_ipsae", True)),
                prediction_config_version=2,
            )
            label = "Boltz-2 target-MSA refold"
            description = "Existing completed Boltz-2 run; target MSA enabled, structural templates off."
        nodes.append({
            "id": f"adopted-{kind}-{index}",
            "type": "studio",
            "position": {"x": index * 330, "y": 180},
            "data": {"kind": kind, "label": label, "description": description, "params": params},
        })
    edges = [
        {"id": "adopted-edge-1", "source": nodes[0]["id"], "sourceHandle": "out", "target": nodes[1]["id"], "targetHandle": "in"},
        {"id": "adopted-edge-2", "source": nodes[1]["id"], "sourceHandle": "out", "target": nodes[2]["id"], "targetHandle": "in"},
    ]
    graph = {
        "schema_version": SCHEMA,
        "id": hashlib.md5(name.encode()).hexdigest(),
        "name": name,
        "nodes": nodes,
        "edges": edges,
        "viewport": {"x": 0, "y": 0, "zoom": 0.72},
    }
    return clean_graph(graph)


def _run_dir(step: dict[str, Any]) -> Path:
    task_group, run_id = step["task_group"], step["run_id"]
    path = (runs_root() / task_group / run_id).resolve()
    path.relative_to(runs_root().resolve())
    return path


def _child_metadata(run_dir: Path, *, name: str, pipeline_id: str, pipeline_dir: Path,
                    step_label: str, index: int, previous: dict[str, Any]) -> dict[str, Any]:
    metadata = read_json(run_dir / "metadata.json")
    for key in _ADOPTION_KEYS:
        if key in metadata:
            previous[key] = metadata[key]
    metadata.update({
        "campaign_name": name,
        "campaign_id": pipeline_id,
        "campaign_step": step_label,
        "campaign_step_index": index,
        "pipeline_parent_run_dir": str(pipeline_dir),
    })
    write_json(run_dir / "metadata.json", metadata)
    return metadata


def adopt_completed_pipelines(groups: list[dict[str, Any]], *, owner_id: str, project_id: str,
                              resolved_rows: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    """Create completed summary roots and tag existing runs. No worker is queued."""
    from mn_molecular_design_studio.workspace import mark_studio_job

    created: list[Path] = []
    changed: list[tuple[Path, dict[str, Any]]] = []
    output: list[dict[str, Any]] = []
    try:
        for group in groups:
            name = group["name"]
            steps = group["steps"]
            signature = _signature(project_id, owner_id, [group])
            # Retry-safe: reuse an existing parent for this exact group.
            existing = None
            for candidate in (runs_root() / "campaign").glob("*"):
                adopted = read_json(candidate / "metadata.json").get("pipeline_adoption")
                if isinstance(adopted, dict) and adopted.get("signature") == signature:
                    existing = candidate
                    break
            if existing is not None:
                output.append({
                    "name": read_json(existing / "input.json").get("params", {}).get("pipeline_name"),
                    "pipeline_run_dir": str(existing),
                    "pipeline_run_id": existing.name,
                    "pipeline_job_code": read_json(existing / "metadata.json").get("job_code"),
                    "step_run_refs": read_json(existing / "adoption_manifest.json").get("step_run_refs", []),
                })
                continue
            refs = [f"{step['task_group']}/{step['run_id']}" for step in steps]
            step_paths = [_run_dir(step) for step in steps]
            # A job may be adopted once. A same-signature retry was returned
            # above; changing its parent would silently steal it from a group.
            for ref, path in zip(refs, step_paths):
                child_metadata = read_json(path / "metadata.json")
                if child_metadata.get("campaign_id") or child_metadata.get("pipeline_parent_run_dir"):
                    raise ValueError(f"{ref} already belongs to another pipeline.")
            step_data = []
            for spec, path in zip(steps, step_paths):
                input_payload = read_json(path / "input.json")
                step_data.append({**spec, "run_dir": str(path), "input": input_payload})
            graph = _graph_for(name, step_data)

            parent = create_job(
                "campaign", "visual_pipeline", "visual_pipeline",
                {"adoption_mode": "existing_completed_runs", "step_run_refs": refs},
                {"pipeline_name": name, "adopted_existing_jobs": True, "adoption_signature": signature},
            )
            created.append(parent.run_dir)
            mark_studio_job({"job": {"run_dir": str(parent.run_dir)}}, project_id=project_id, owner_id=owner_id)
            parent_meta = read_json(parent.run_dir / "metadata.json")
            parent_meta["campaign_name"] = name
            parent_meta["pipeline_adoption"] = {
                "mode": "existing_completed_runs",
                "signature": signature,
                "step_run_refs": refs,
                "created_at": utc_now(),
            }
            write_json(parent.run_dir / "metadata.json", parent_meta)
            write_json(parent.run_dir / "graph_snapshot.json", graph)
            node_states = {}
            for index, (spec, data, path) in enumerate(zip(steps, step_data, step_paths)):
                ref = refs[index]
                step_id = graph["nodes"][index]["id"]
                metadata = read_json(path / "metadata.json")
                before = {key: metadata[key] for key in _ADOPTION_KEYS if key in metadata}
                changed.append((path, before))
                _child_metadata(
                    path, name=name, pipeline_id=parent.run_id, pipeline_dir=parent.run_dir,
                    step_label=spec["stage_label"], index=index + 1, previous=before,
                )
                node_states[step_id] = {
                    "status": "completed",
                    "job_run_dir": str(path),
                    "job_code": str(resolved_rows[ref].get("job_code") or metadata.get("job_code") or ""),
                    "candidate_count": int(read_json(path / "result.json").get("metrics", {}).get("candidate_count") or 1),
                    "completed_at": str(metadata.get("completed_at") or metadata.get("updated_at") or ""),
                    "adopted_existing_run": ref,
                    "output": {
                        "run_dir": str(path),
                        "run_id": path.name,
                        "task_group": spec["task_group"],
                        "job_code": str(metadata.get("job_code") or ""),
                        "candidate_count": int(read_json(path / "result.json").get("metrics", {}).get("candidate_count") or 1),
                        "adopted_existing_run": ref,
                    },
                }
            write_json(parent.run_dir / "graph_execution.json", {
                "graph_id": graph["id"],
                "status": "completed",
                "adopted_existing_runs": True,
                "nodes": node_states,
                "updated_at": utc_now(),
            })
            write_json(parent.run_dir / "adoption_manifest.json", {
                "mode": "existing_completed_runs",
                "signature": signature,
                "project_id": project_id,
                "owner_id": owner_id,
                "pipeline_name": name,
                "step_run_refs": refs,
                "step_original_metadata": [
                    {"task_group": spec["task_group"], "run_id": spec["run_id"], "fields": before}
                    for spec, (_, before) in zip(steps, changed[-len(steps):])
                ],
            })
            finish_job(parent.run_dir, True, {
                "outputs": {"pipeline_nodes": node_states, "adopted_existing_runs": refs},
                "metrics": {"node_count": len(steps), "completed_steps": len(steps), "adopted_existing_runs": True},
                "downstream_artifacts": {"graph_snapshot": "graph_snapshot.json", "graph_execution": "graph_execution.json", "adoption_manifest": "adoption_manifest.json"},
            })
            # The summary parent represents the historical work interval, not
            # the time at which these completed runs were grouped.
            child_metadata_rows = [read_json(path / "metadata.json") for path in step_paths]
            historical_starts = [
                str(item.get("worker_started_at") or item.get("created_at") or "")
                for item in child_metadata_rows
                if item.get("worker_started_at") or item.get("created_at")
            ]
            historical_finishes = [
                str(item.get("completed_at") or item.get("updated_at") or "")
                for item in child_metadata_rows
                if item.get("completed_at") or item.get("updated_at")
            ]
            parent_metadata = read_json(parent.run_dir / "metadata.json")
            if historical_starts:
                parent_metadata["worker_started_at"] = min(historical_starts)
            if historical_finishes:
                parent_metadata["completed_at"] = max(historical_finishes)
            parent_metadata["adopted_at"] = utc_now()
            write_json(parent.run_dir / "metadata.json", parent_metadata)
            output.append({
                "name": name,
                "pipeline_run_dir": str(parent.run_dir),
                "pipeline_run_id": parent.run_id,
                "pipeline_job_code": read_json(parent.run_dir / "metadata.json").get("job_code"),
                "step_run_refs": refs,
            })
    except Exception:
        for path, previous in reversed(changed):
            metadata = read_json(path / "metadata.json")
            for key in _ADOPTION_KEYS:
                metadata.pop(key, None)
            metadata.update(previous)
            metadata["updated_at"] = utc_now()
            write_json(path / "metadata.json", metadata)
        for path in reversed(created):
            shutil.rmtree(path, ignore_errors=True)
        from mn_molecular_design_studio.job_catalog import reconcile
        reconcile()
        raise

    from mn_molecular_design_studio.job_catalog import reconcile
    reconcile()
    return output


def undo_pipeline_adoption(run_dir: Path, *, owner_id: str, project_id: str) -> dict[str, Any]:
    """Restore saved metadata fields and remove only the display-only parent run."""
    manifest = read_json(run_dir / "adoption_manifest.json")
    metadata = read_json(run_dir / "metadata.json")
    if manifest.get("mode") != "existing_completed_runs" or metadata.get("tool") != "visual_pipeline":
        raise ValueError("This is not an adopted existing-jobs pipeline.")
    if manifest.get("owner_id") != owner_id or manifest.get("project_id") != project_id:
        raise PermissionError("The adopted pipeline belongs to a different owner or project.")
    for original in manifest.get("step_original_metadata", []):
        path = _run_dir(original)
        child = read_json(path / "metadata.json")
        if str(child.get("campaign_id") or "") != run_dir.name:
            raise ValueError(f"Pipeline membership changed for {original['task_group']}/{original['run_id']}.")
        for key in _ADOPTION_KEYS:
            child.pop(key, None)
        child.update(original.get("fields") or {})
        write_json(path / "metadata.json", child)
    shutil.rmtree(run_dir)
    from mn_molecular_design_studio.job_catalog import reconcile
    reconcile()
    return {"status": "undone", "removed_pipeline_run_id": run_dir.name, "restored_step_count": len(manifest.get("step_run_refs") or [])}
