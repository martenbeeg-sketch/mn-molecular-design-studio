"""External active-learning results and bounded, explicit final promotion.

The Parquet files are the complete scientific record. The molecular Database
receives at most 1,000 final ranked rows, never the loop's training dataset.
"""
from __future__ import annotations

import hashlib
import json
import math
import socket
from pathlib import Path

from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id, _json, _now

MAX_FINAL_IMPORT = 1000
RESULT_DIRECTORY = "artifacts/active_learning"


def _limit(value):
    if type(value) is not int or not 0 <= value <= MAX_FINAL_IMPORT:
        raise ValueError("Final Database import must be between 0 and 1,000 compounds.")
    return value


def _head(path, limit):
    """Read a bounded prefix of the controller's score-ordered final table."""
    import pyarrow.parquet as pq
    table = pq.ParquetFile(path)
    rows = []
    if limit:
        for batch in table.iter_batches(batch_size=min(limit, MAX_FINAL_IMPORT)):
            rows.extend(batch.to_pylist()[:limit - len(rows)])
            if len(rows) >= limit:
                break
    return table.metadata.num_rows, rows


def external_results(run_dir, user, *, preview_limit=100):
    from mn_molecular_design_studio.molecular_database_results import _authorized_run
    selected = _authorized_run(run_dir, user)[0]
    return screening_results_payload(selected, preview_limit=preview_limit)


def screening_results_payload(run_dir, *, preview_limit=100):
    """Bounded trusted job-inspection payload; public callers authorize first."""
    selected = Path(run_dir).resolve()
    path = selected / RESULT_DIRECTORY / "docking_scores.parquet"
    if not path.is_file():
        return None
    count, rows = _head(path, _limit(preview_limit))
    manifest_path = selected / RESULT_DIRECTORY / "external_results.json"
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    manifest = read_json(manifest_path) if manifest_path.is_file() else {}
    history = read_json(selected / RESULT_DIRECTORY / "round_history.json")
    convergence = read_json(selected / RESULT_DIRECTORY / "score_convergence.json")
    summaries = {row["round"]: row["score_summaries"] for row in convergence.get("rounds", [])}
    rounds = [{**row, "score_summaries": row.get("score_summaries") or summaries.get(row.get("round"), [])}
              for row in history.get("rounds", [])]
    graph = read_json(selected / "graph_snapshot.json")
    decision = next((node.get("data", {}).get("params", {}) for node in graph.get("nodes", [])
                     if node.get("data", {}).get("kind") == "active_decision"), {})
    return {"compound_count": count, "rows": rows, "scores_path": str(path),
            "rounds": rounds, "stop_reason": history.get("stop_reason"),
            "score_field": history.get("score_field"), "score_units": history.get("score_units"),
            "score_direction": history.get("score_direction"),
            "stopping_criteria": {key: decision[key] for key in
                                  ("top_n", "score_tolerance", "stable_rounds", "max_iterations") if key in decision},
            "manifest_path": str(manifest_path), "selection": manifest.get("database_selection", {}),
            "files": manifest.get("files", {})}


def publish_final_selection(run_dir, user, limit, *, database=None):
    """Trusted controller adapter; call after final Parquet export, before finish.

    Retrying this adapter is idempotent. Reopening an immutable final result may
    expand its curated shortlist up to 1,000 rows; it never removes user edits.
    Public UI calls must use import_final_results for job-access/status checks.
    """
    limit = _limit(limit)
    run_dir = Path(run_dir).resolve()
    data = run_dir / RESULT_DIRECTORY
    count, rows = _head(data / "docking_scores.parquet", limit)
    marker = json.loads((run_dir / "molecular_design_studio.json").read_text())
    owner, project = str(marker["studio_owner_id"]), str(marker.get("project_id") or "")
    if user["id"] != owner and user.get("role") != "admin":
        raise PermissionError("Only the campaign owner can import its final shortlist.")
    history_path = data / "round_history.json"
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    history = read_json(history_path) if history_path.is_file() else {}
    graph_path = run_dir / "graph_snapshot.json"
    graph = read_json(graph_path) if graph_path.is_file() else {}
    params = {n["data"]["kind"]: n["data"].get("params", {}) for n in graph.get("nodes", [])}
    docking = params.get("active_docking", {})
    target = params.get("active_target", {})
    field = str(history.get("score_field") or docking.get("score_field") or "score")
    units = str(history.get("score_units") or docking.get("score_units") or "")
    direction = str(history.get("score_direction") or docking.get("score_direction") or "Lower is better")
    stamp = _now()
    group_id, imported = None, 0
    if rows:
        # Validate the bounded selection before writing any Database metadata.
        records = []
        origin = {"type": "job", "collection_kind": "curated", "saved_selection": True,
                  "final_screening_selection": True, "task_group": "campaign", "run_id": run_dir.name,
                  "run_dir": str(run_dir), "external_scores": str(data / "docking_scores.parquet"),
                  "screening_library_id": history.get("screening_library_id"), "source_node": socket.gethostname()}
        for rank, row in enumerate(rows, 1):
            score = float(row.get("score", row.get("best_score")))
            if not math.isfinite(score) or not row.get("smiles") or not row.get("compound_id"):
                raise ValueError("A selected final result lacks a valid score, ID or SMILES.")
            metadata = {key: row.get(key) for key in ("compound_id", "compound_key", "library_id", "compound_name",
                        "first_round", "best_round", "worker_job_code", "worker_run_dir", "worker_node_id")}
            metadata.update(final_rank=rank, source_smiles=row["smiles"], external_scores=origin["external_scores"])
            record = {"name": str(row.get("compound_name") or row.get("name") or row["compound_id"]),
                      "components": [{"type": "small_molecule", "chain": "L", "smiles": row["smiles"]}],
                      "metadata": metadata, "files": []}
            if row.get("pose_path"):
                pose = (run_dir / row["pose_path"]).resolve()
                # Only the verified, retained organizer poses can enter this snapshot.
                if pose.is_relative_to((data / "retained_poses").resolve()) and pose.is_file():
                    content = pose.read_bytes()
                    record["files"] = [{"name": pose.name, "content": content, "label": "Final retained docking pose",
                                        "portable": True, "portable_sha256": hashlib.sha256(content).hexdigest(),
                                        "external_path": str(pose), "dimensions": 3}]
            records.append((row, rank, score, record))
        database = database or MolecularDatabase()
        from mn_molecular_design_studio.bound_docking import screening_receptor, register_receptor
        receptor = screening_receptor(run_dir, target) if any(record["files"] for _, _, _, record in records) else None
        target_context = register_receptor(database, user, receptor, project_id=project,
            name=str(target.get("target_name") or "Docking receptor")) if receptor else {}
        group_id = database.create_group(user, f"{graph.get('name') or run_dir.name} · Final shortlist",
            project_id=project, import_method="tool_output", origin=origin,
            group_id=_id("screening-final-", str(run_dir)))
        with database.connection() as connection:
            for row, rank, score, record in records:
                source_key = f"{socket.gethostname()}/campaign/{run_dir.name}:final:{row.get('library_id')}:{row.get('compound_key', rank)}"
                files = record["files"]
                molecule_id = database._put(connection, group_id, {**record, "files": []}, source_key=source_key,
                    origin={**origin, "compound_id": row["compound_id"], "final_rank": rank}, stamp=stamp,
                    preserve_existing=True)
                representation = None
                for item in files:
                    path = database._store_file(item["name"], item["content"])
                    structure_context = {**origin, "compound_id": row["compound_id"], "final_rank": rank,
                        **{key: value for key, value in item.items() if key not in {"name", "content", "label"}},
                        **target_context, "database_structure_policy": "best_docking_pose"}
                    existing = connection.execute("SELECT id FROM molecule_structures WHERE molecule_id=? AND group_id=? AND path=? "
                        "AND json_extract(context,'$.run_id')=? AND json_extract(context,'$.compound_id')=?",
                        (molecule_id, group_id, str(path), run_dir.name, row["compound_id"])).fetchone()
                    representation = existing[0] if existing else _id("structure-", _json([molecule_id, group_id, str(path), structure_context]))
                    connection.execute("INSERT INTO molecule_structures VALUES (?,?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET context=excluded.context", (
                        representation, molecule_id, group_id, str(path), item["label"], path.suffix.lstrip("."), _json(structure_context), stamp))
                context = {**origin, "compound_id": row["compound_id"], "units": units,
                           "primary_ranking_metric": field, "ranking_direction": direction,
                           **target_context, "representation_id": representation, "settings": docking,
                           "best_round": row.get("best_round"), "first_round": row.get("first_round")}
                observation_id = _id("score-", source_key)
                connection.execute("INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?) "
                    "ON CONFLICT(id) DO UPDATE SET metrics=excluded.metrics,context=excluded.context", (
                    observation_id, molecule_id, group_id, str(docking.get("engine") or "docking"),
                    run_dir.name, _json({field: score, "final_rank": rank}), _json(context), stamp))
            imported = connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=?", (group_id,)).fetchone()[0]
    selection = {"requested_limit": limit, "imported_count": imported, "group_id": group_id,
                 "maximum_import": MAX_FINAL_IMPORT, "ranking": "final cumulative docking scores"}
    files = {name: str(path) for name in ("docking_scores.parquet", "failed_compounds.parquet", "round_history.json",
             "score_convergence.json", "receptor.pdb", "retained_poses.json", "additional_pose_archives.json", "model.pt")
             if (path := data / name).is_file()}
    manifest = {"schema_version": "studio.external-screening-results.v1", "compound_count": count,
                "database_collection_policy": "manual", "database_selection": selection, "files": files,
                "screening_library_id": history.get("screening_library_id"),
                "screening_feature_manifest": history.get("screening_feature_manifest"),
                "centroid_manifest": history.get("centroid_manifest"), "updated_at": stamp}
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json
    write_json(data / "external_results.json", manifest)
    return selection


def import_final_results(run_dir, user, limit, *, database=None):
    from mn_molecular_design_studio.molecular_database_results import _authorized_run
    selected = _authorized_run(run_dir, user)[0]
    metadata = json.loads((selected / "metadata.json").read_text())
    result = json.loads((selected / "result.json").read_text())
    if metadata.get("status") != "completed" or not result.get("success"):
        raise ValueError("Import a final shortlist only after the campaign has completed successfully.")
    return publish_final_selection(selected, user, limit, database=database)
