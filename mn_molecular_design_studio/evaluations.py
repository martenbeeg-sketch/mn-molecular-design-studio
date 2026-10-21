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

from mn_molecular_design_studio.output_vocabulary import PREFIX, VERSION, normalize_metrics

PROJECTION_VERSION = "evaluation.v4:" + VERSION


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
    representation = context.get("structure_digest") or context.get("_replica_representation") or context.get("representation_id") or row.get("pose_path") or candidate.get("complex_pdb") or candidate.get("binder_pdb") or ""
    settings = context.get("settings") or {}
    runtime_fields = {"source_run_dir", "candidates_jsonl", "job_name", "campaign_name", "cpu_cores", "cpu_threads", "nprocs", "workers", "gpu_device", "backend", "output_dir"}
    scientific_settings = {k: v for k, v in settings.items() if k not in runtime_fields}
    instance = identifier("instance-", [context.get("node_id"), context.get("task_group"), run, candidate_id, candidate.get("stage")])
    subject = identifier("subject-", [observation.get("molecule_id"), target, target_revision, representation,
                                     candidate.get("binder_chains"), candidate.get("target_chains"), context.get("prepared_variant_id")])
    canonical_id = identifier("evaluation-", [context.get("workstation"), context.get("_replica_observation") or observation["id"]])
    return {"id": canonical_id, "observation_id": observation["id"], "instance_id": instance, "subject_id": subject,
            "target_id": target, "target_revision": str(target_revision or ""),
            "representation_id": str(context.get("representation_id") or ""),
            "representation": str(representation), "candidate_id": candidate_id,
            "evaluator": str(candidate.get("source_tool") or observation.get("namespace") or ""),
            "run_id": run, "stage": str(candidate.get("stage") or context.get("stage") or ""),
            "parameters_digest": identifier("params-", scientific_settings),
            "engine_version": str((candidate.get("raw_metadata") or {}).get("engine_version") or settings.get("engine_version") or ""),
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


def rebuild_missing(connection):
    rows = connection.execute("SELECT o.* FROM molecule_observations o LEFT JOIN molecule_evaluations e ON e.observation_id=o.id "
                              "WHERE e.observation_id IS NULL OR e.projection_version!=?", (PROJECTION_VERSION,)).fetchall()
    for raw in rows:
        observation = dict(raw)
        descriptor = evaluation_descriptor(observation)
        metrics, context = json.loads(observation["metrics"]), json.loads(observation["context"])
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
    normalized = normalize_metrics(observation["metrics"], observation["namespace"], observation["context"])
    scores = {observation["namespace"] + "." + key: value for key, value in observation["metrics"].items()}
    scores.update({PREFIX + key: value for key, value in normalized["metrics"].items()})
    sources = {field: {"observation_id": observation["id"], "run_id": observation["run_id"],
                       "namespace": observation["namespace"], "created_at": observation["created_at"],
                       **normalized["provenance"].get(field.removeprefix(PREFIX), {})} for field in scores}
    return scores, sources, normalized
