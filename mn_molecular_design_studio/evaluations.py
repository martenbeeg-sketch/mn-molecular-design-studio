"""Evaluation identity and rebuildable typed projections of scientific observations.

Raw observations are replicated. These indexes are local and may always be rebuilt.
An evaluation is one model/pose/calculation, not a bag of independently latest scores.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from functools import lru_cache

from mn_molecular_design_studio.output_vocabulary import PREFIX, VERSION, normalize_metrics, normalize_parameters

PROJECTION_VERSION = "evaluation.v6:" + VERSION


def _projectable_metrics(metrics, namespace):
    if namespace != "rosetta_metrics":
        return metrics
    # The installed standalone scorer mislabeled total interface dSASA /
    # bound binder SASA as a binder percentage. Keep historical raw evidence,
    # but exclude it from displayed scores and query/filter indexes.
    return {key: value for key, value in metrics.items() if key not in {
        "rosetta_interface_percent_of_binder_in_intf", "interface_percent_of_binder_in_intf",
    }}


def encoded(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, default=str, allow_nan=False)


def identifier(prefix, value):
    return prefix + hashlib.sha256(encoded(value).encode()).hexdigest()[:32]


@lru_cache(maxsize=4096)
def _file_digest(path, size, modified):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def file_digest(path):
    try:
        source = Path(path)
        stat = source.stat()
        return _file_digest(str(source.resolve()), stat.st_size, stat.st_mtime_ns) if source.is_file() else ""
    except (OSError, TypeError):
        return ""


def evaluation_descriptor(observation):
    context = observation.get("context") or {}
    if isinstance(context, str):
        context = json.loads(context)
    candidate = context.get("candidate") or {}
    row = context.get("row") or {}
    run = str(observation.get("run_id") or context.get("run_id") or "")
    candidate_id = str(context.get("candidate_id") or candidate.get("candidate_id") or row.get("compound_id") or "")
    target = str(context.get("target_molecule_id") or context.get("target_id") or context.get("target") or "")
    target_revision = context.get("target_structure_digest") or ""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
    pose = context.get("pose_identity") or candidate.get("pose_identity") or row.get("pose_identity") or {}
    pose_id = str(pose.get("id") or context.get("pose_id") or "") if isinstance(pose, dict) else str(pose or "")
    model_id = str(context.get("model_id") or candidate.get("model_id") or candidate.get("model_index")
                   or row.get("model_id") or selection.get("model_id") or selection.get("model_index")
                   or raw.get("model_id") or raw.get("model_index") or "")
    binder_chains = (candidate.get("binder_chains") or selection.get("binder_chains")
                     or raw.get("binder_chains") or context.get("binder_chains") or [])
    target_chains = (candidate.get("target_chains") or selection.get("target_chains")
                     or raw.get("target_chains") or context.get("target_chains") or [])
    ligand_selection = (selection.get("ligand_selection") or raw.get("ligand_selection")
                        or context.get("ligand_selection") or {})
    selected_components = selection.get("components") or selection.get("selected_components") or []
    scientific_selection = {"binder_chains": sorted(map(str, binder_chains)),
                            "target_chains": sorted(map(str, target_chains)),
                            "ligand_selection": ligand_selection,
                            "components": selected_components}
    selection_digest = identifier("selection-", scientific_selection)
    source_structure_digest = str(context.get("source_structure_digest") or raw.get("source_structure_digest") or "")
    calculation_structure_digest = str(context.get("calculation_structure_digest") or context.get("structure_digest")
                                       or context.get("pose_digest") or row.get("pose_digest") or "")
    structure_digest = calculation_structure_digest or source_structure_digest
    representation = (structure_digest or context.get("_replica_representation") or context.get("representation_id")
                      or pose_id or row.get("pose_path") or candidate.get("complex_pdb") or candidate.get("binder_pdb") or "")
    settings = context.get("settings") or raw.get("analysis_parameters") or {}
    runtime_fields = {"source_run_dir", "candidates_jsonl", "job_name", "campaign_name", "cpu_cores", "cpu_threads", "nprocs", "workers", "gpu_device", "backend", "output_dir"}
    scientific_settings = {k: v for k, v in settings.items() if k not in runtime_fields}
    parameter_projection = normalize_parameters(scientific_settings, candidate.get("source_tool") or observation.get("namespace"))
    parameters_digest = identifier("params-", scientific_settings)
    instance = identifier("instance-", [context.get("node_id"), context.get("task_group"), run, candidate_id,
                                        candidate.get("stage"), model_id, pose_id, structure_digest,
                                        source_structure_digest, selection_digest, parameters_digest])
    subject = identifier("subject-", [observation.get("molecule_id"), target, target_revision, representation,
                                     selection_digest, context.get("prepared_variant_id"), model_id, pose_id])
    canonical_id = identifier("evaluation-", [context.get("workstation"), context.get("_replica_observation") or observation["id"]])
    return {"id": canonical_id, "observation_id": observation["id"], "instance_id": instance, "subject_id": subject,
            "target_id": target, "target_revision": str(target_revision or ""),
            "representation_id": str(context.get("representation_id") or ""),
            "representation": str(representation), "structure_digest": structure_digest,
            "source_structure_digest": source_structure_digest,
            "calculation_structure_digest": calculation_structure_digest,
            "model_id": model_id, "pose_id": pose_id, "selection_id": selection_digest,
            "binder_chains": sorted(map(str, binder_chains)), "target_chains": sorted(map(str, target_chains)),
            "ligand_selection": ligand_selection, "candidate_id": candidate_id,
            "evaluator": str(candidate.get("source_tool") or observation.get("namespace") or ""),
            "run_id": run, "stage": str(candidate.get("stage") or context.get("stage") or ""),
            "parameters_digest": parameters_digest, "native_parameters": scientific_settings,
            "normalized_parameters": parameter_projection["normalized"],
            "parameter_provenance": parameter_projection["provenance"],
            "parameter_conflicts": parameter_projection["conflicts"],
            "engine_version": str(raw.get("engine_version") or settings.get("engine_version") or ""),
            "source_workstation": str(context.get("workstation") or context.get("source_workstation") or ""),
            "image_digest": str(context.get("image_digest") or settings.get("image_digest") or ""),
            "status": str(context.get("status") or "completed"),
            "preferred": bool(context.get("is_representative") or context.get("representation_id")) and context.get("status", "completed") == "completed",
            "context_complete": bool(target and representation), "created_at": observation.get("created_at", "")}


def install_schema(connection):
    connection.executescript('''
        CREATE TABLE IF NOT EXISTS molecule_evaluations (
            observation_id TEXT PRIMARY KEY REFERENCES molecule_observations(id) ON DELETE CASCADE,
            molecule_id TEXT NOT NULL, group_id TEXT NOT NULL, subject_id TEXT NOT NULL,
            target_id TEXT NOT NULL, evaluator TEXT NOT NULL, run_id TEXT NOT NULL,
            preferred INTEGER NOT NULL, created_at TEXT NOT NULL, descriptor TEXT NOT NULL,
            projection_version TEXT NOT NULL);
        CREATE INDEX IF NOT EXISTS molecular_evaluation_selection ON molecule_evaluations(molecule_id,target_id,evaluator,run_id,created_at);
        CREATE TABLE IF NOT EXISTS molecule_metric_values (
            observation_id TEXT NOT NULL REFERENCES molecule_evaluations(observation_id) ON DELETE CASCADE,
            field TEXT NOT NULL, value_type TEXT NOT NULL, value_number REAL, value_json TEXT NOT NULL,
            definition TEXT NOT NULL, PRIMARY KEY(observation_id,field));
        CREATE INDEX IF NOT EXISTS molecular_metric_filter ON molecule_metric_values(field,value_number,observation_id);
        CREATE TRIGGER IF NOT EXISTS molecular_evaluation_invalidate AFTER UPDATE ON molecule_observations BEGIN
            DELETE FROM molecule_metric_values WHERE observation_id=OLD.id;
            DELETE FROM molecule_evaluations WHERE observation_id=OLD.id;
        END;
        CREATE TRIGGER IF NOT EXISTS molecular_evaluation_delete AFTER DELETE ON molecule_observations BEGIN
            DELETE FROM molecule_metric_values WHERE observation_id=OLD.id;
            DELETE FROM molecule_evaluations WHERE observation_id=OLD.id;
        END;
    ''')


def rebuild_missing(connection, *, batch_size=500):
    """Recompute missing or stale projections in bounded, versioned batches."""
    if type(batch_size) is not int or not 1 <= batch_size <= 5000:
        raise ValueError("Projection rebuild batches must be between 1 and 5,000 rows.")
    rows = connection.execute("SELECT o.* FROM molecule_observations o LEFT JOIN molecule_evaluations e ON e.observation_id=o.id "
                              "WHERE e.observation_id IS NULL OR e.projection_version!=? ORDER BY o.id LIMIT ?",
                              (PROJECTION_VERSION, batch_size)).fetchall()
    for raw in rows:
        observation = dict(raw)
        descriptor = evaluation_descriptor(observation)
        metrics, context = json.loads(observation["metrics"]), json.loads(observation["context"])
        metrics = _projectable_metrics(metrics, observation["namespace"])
        normalized = normalize_metrics(metrics, observation["namespace"], context)
        connection.execute("DELETE FROM molecule_metric_values WHERE observation_id=?", (observation["id"],))
        connection.execute("INSERT OR REPLACE INTO molecule_evaluations VALUES (?,?,?,?,?,?,?,?,?,?,?)", (
            observation["id"], observation["molecule_id"], observation["group_id"], descriptor["subject_id"], descriptor["target_id"],
            descriptor["evaluator"], descriptor["run_id"], int(descriptor["preferred"]), observation["created_at"], encoded(descriptor), PROJECTION_VERSION))
        fields = [(observation["namespace"] + "." + key, value, {}) for key, value in metrics.items()]
        fields.extend((PREFIX + key, value, normalized["provenance"].get(key, {})) for key, value in normalized["metrics"].items())
        for key, value, definition in fields:
            number = float(value) if type(value) in {int, float} and math.isfinite(value) else None
            kind = "boolean" if isinstance(value, bool) else "number" if number is not None else "categorical" if isinstance(value, str) else "array" if isinstance(value, list) else "object"
            connection.execute("INSERT OR REPLACE INTO molecule_metric_values VALUES (?,?,?,?,?,?)", (
                observation["id"], key, kind, number, encoded(value), encoded(definition)))
    return len(rows)


def selection_clause(filters, parameters, *, alias="e"):
    filters = filters or {}
    allowed = {"target_id", "evaluator", "run_id", "evaluation_id", "subject_id", "evaluation_ids"}
    if set(filters) - allowed:
        raise ValueError("Unsupported evaluation selection fields.")
    clauses = []
    for key in sorted(filters):
        value = filters[key]
        if value in (None, "", []):
            continue
        column = "json_extract(" + alias + ".descriptor,'$.id')" if key in {"evaluation_id", "evaluation_ids"} else f"{alias}.{key}"
        if key == "evaluation_ids":
            if not isinstance(value, list) or len(value) > 100000:
                raise ValueError("Evaluation IDs must be a list of at most 100,000 entries.")
            clauses.append(f"{column} IN (SELECT value FROM json_each(?))")
            parameters.append(encoded(value))
        else:
            if not isinstance(value, str):
                raise ValueError("Evaluation selection values must be text.")
            clauses.append(f"{column}=?")
            parameters.append(value)
    return " AND ".join(clauses) or "1"


def project_observation(observation):
    metrics = _projectable_metrics(observation["metrics"], observation["namespace"])
    normalized = normalize_metrics(metrics, observation["namespace"], observation["context"])
    scores = {observation["namespace"] + "." + key: value for key, value in metrics.items()}
    scores.update({PREFIX + key: value for key, value in normalized["metrics"].items()})
    sources = {field: {"observation_id": observation["id"], "run_id": observation["run_id"],
                       "namespace": observation["namespace"], "created_at": observation["created_at"],
                       **normalized["provenance"].get(field.removeprefix(PREFIX), {})} for field in scores}
    return scores, sources, normalized
