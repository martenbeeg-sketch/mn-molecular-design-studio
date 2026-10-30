"""Versioned scientific collections on the existing normalized-candidate contract.

Collections are JSONL artifacts, not scheduler resources. SQLite scratch indexes
bound the memory used by fan-in, identity joins and deduplication. No engine or
campaign is executed here.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import operator
import math
import sqlite3
import tempfile
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

from mn_molecular_design_studio.cache_provenance import file_digest
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import normalize_candidate
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_managed_paths, resolve_stored_path
from mn_molecular_design_studio.measurement_status import finite_number, RANGES

SCHEMA = "studio.candidate-set.v1"
KINDS = {"candidate_merge", "shared_filter", "candidate_analyzer", "candidate_rank"}
COMPARISONS = {">=": operator.ge, "<=": operator.le, ">": operator.gt, "<": operator.lt}


def identity(*values):
    return hashlib.sha256(json.dumps(values, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def json_safe(value):
    """Retain nonfinite evidence explicitly without writing non-standard JSON."""
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def metadata(row):
    return row.get("raw_metadata", {}).get("candidate_set", {})


def iter_candidates(source, *, review=False):
    """Strict streaming reader: corruption is an operational error, never rejection."""
    path = Path(source.get("review_candidates_jsonl") if review and source.get("review_candidates_jsonl") else source["candidates_jsonl"])
    base = Path(source["run_dir"])
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                payload = json.loads(line)
                if not isinstance(payload, dict) or not payload.get("candidate_id"):
                    raise ValueError("candidate_id is required")
                row = resolve_managed_paths(normalize_candidate(payload), run_dir=base)
                info = metadata(row)
                if info and (info.get("schema_version") != SCHEMA or any(not info.get(key) for key in ("target_id", "backbone_id", "variant_id"))):
                    raise ValueError("Unsupported or incomplete CandidateSet identity metadata")
            except (ValueError, TypeError) as exc:
                raise ValueError(f"Invalid candidate record {path}:{number}: {exc}") from exc
            for key in ("complex_pdb", "binder_pdb", "target_pdb"):
                if row.get(key):
                    resolved = resolve_stored_path(row[key], run_dir=base)
                    if resolved:
                        row[key] = str(resolved)
            yield row


@contextmanager
def index():
    with tempfile.TemporaryDirectory(prefix="studio-collection-") as directory:
        db = sqlite3.connect(str(Path(directory) / "index.sqlite"))
        try:
            db.execute("PRAGMA cache_size=-4096")
            yield db
        finally:
            db.close()


def target_identity(target):
    """Original file, exact author residue IDs and selected-chain sequences."""
    from mn_molecular_design_studio.design_validation import read_structure, protein_chains, sequence
    path = Path(target["target_pdb"])
    chains = target.get("target_chains", [])
    if isinstance(chains, str):
        import re
        chains = [s for s in re.split(r"[,;\s]+", chains) if s]
    available = protein_chains(read_structure(path))
    if not chains or len(chains) != len(set(chains)) or not set(chains) <= available.keys():
        raise ValueError("Select unambiguous existing protein target chains.")
    mapping = {chain: {"sequence": sequence(available[chain]),
                       "residues": [list(res.id) for res in available[chain]]} for chain in chains}
    digest = file_digest(path)
    return dict(target_id=identity(digest, chains), source_sha256=digest,
                chains=mapping, selected_chains=chains,
                constraints={k: copy.deepcopy(target.get(k, "")) for k in ("hotspots", "coldspots")})


def sequence_state(row):
    seq = row.get("binder_sequence")
    if row.get("stage") == "generation.backbone":
        return "backbone_only"
    if not isinstance(seq, str) or not seq or any(c not in "ACDEFGHIKLMNPQRSTVWY:" for c in seq) or any(not p for p in seq.split(":")):
        return "unresolved"
    return "sequence_ready"


def native_filter(row):
    metrics = row.get("metrics", {})
    value = metrics.get("native_pass_filters")
    if type(value) is bool:
        return "passed" if value else "failed"
    if value is not None:
        return "invalid_metadata"
    campaign = row.get("raw_metadata", {}).get("distributed_campaign", {})
    value = campaign.get("native_acceptance")
    if value in {"accepted", "rejected"}:
        return "passed" if value == "accepted" else "failed"
    return "not_evaluated" if value in {None, "unavailable"} else "unresolved"


def artifact_evidence(row):
    result = {}
    for key in ("complex_pdb", "binder_pdb", "target_pdb"):
        if row.get(key):
            path = Path(row[key])
            result[key] = {"path": str(path), "status": "present" if path.is_file() else "missing_evidence",
                           "sha256": file_digest(path) if path.is_file() else None}
    return result


def _campaign_acceptance(row):
    raw = row.get("raw_metadata") or {}
    campaign = raw.get("distributed_campaign") or {}
    if not campaign:
        return None
    validated = campaign.get("validated_acceptance")
    native = campaign.get("native_acceptance")
    level = campaign.get("acceptance_level")
    if level not in {"native", "validated"}:
        level = "validated" if validated not in {None, "pending", "unavailable"} else "native"
    status = validated if level == "validated" else native
    if status == "unavailable":
        status = "not_evaluated"
    return dict(
        campaign_id=campaign.get("campaign_id"),
        attempt_id=campaign.get("attempt_id"),
        attempt_seed=campaign.get("attempt_seed"),
        generator_id=campaign.get("generator_id"),
        objective_type=campaign.get("objective_type"),
        acceptance_level=level,
        status=status or "pending",
        reason=campaign.get("validation_decision_reason") or campaign.get("decision_reason", ""),
    )


def _campaign_predictions(row):
    """Project completed controller validations into the existing CandidateSet prediction schema."""
    raw = row.get("raw_metadata") or {}
    campaign = raw.get("distributed_campaign") or {}
    records = []
    evidence_rows = raw.get("evidence_candidates")
    if not isinstance(evidence_rows, list):
        return records
    for ordinal, evidence in enumerate(evidence_rows):
        if not isinstance(evidence, dict) or evidence.get("stage") not in {"complex_refolding", "analysis"}:
            continue
        if not row.get("binder_sequence") or evidence.get("binder_sequence") != row.get("binder_sequence"):
            continue
        evidence_raw = evidence.get("raw_metadata") or {}
        validation = evidence_raw.get("distributed_campaign_validation") or {}
        stage_index = validation.get("stage_index", ordinal)
        branch = f"campaign-validation-{stage_index}"
        key = "analysis_results" if evidence.get("stage") == "analysis" else "predictions"
        identity_field = "analysis_result_id" if key == "analysis_results" else "prediction_id"
        record_id = identity(
            "distributed-campaign-evidence", campaign.get("campaign_id"),
            campaign.get("attempt_id"), stage_index, evidence.get("candidate_id"),
        )
        metrics = evidence.get("metrics") if isinstance(evidence.get("metrics"), dict) else {}
        status = validation.get("status") or ("completed" if metrics else "missing_evidence")
        if status not in {"completed", "failed", "missing_evidence", "unresolved"}:
            status = "unresolved"
        candidate = copy.deepcopy(evidence)
        if key == "analysis_results":
            analysis_info = candidate.setdefault("raw_metadata", {}).setdefault("candidate_set", {})
            analysis_info.update(
                schema_version=SCHEMA,
                target_id=metadata(row).get("target_id"),
                backbone_id=metadata(row).get("backbone_id"),
                variant_id=metadata(row).get("variant_id"),
                record_kind="measurement",
                analysis_result_id=record_id,
                measurement_node=branch,
            )
        records.append({
            identity_field: record_id,
            "branch": branch,
            "predictor": evidence.get("source_tool") or evidence.get("tool"),
            "status": status,
            "candidate": candidate,
            "measurements": {
                metric: {"status": "valid", "value": value}
                for metric, value in metrics.items()
                if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
            },
        })
    return records


def _attach_campaign_evidence(row):
    info = metadata(row)
    if not info:
        return
    acceptance = _campaign_acceptance(row)
    if acceptance:
        info["campaign_acceptance"] = acceptance
    for record in _campaign_predictions(row):
        key = "analysis_results" if "analysis_result_id" in record else "predictions"
        entries = info.setdefault(key, [])
        record_id = record.get("prediction_id") or record.get("analysis_result_id")
        if not any((item.get("prediction_id") or item.get("analysis_result_id")) == record_id for item in entries):
            entries.append(record)
    if info.get("predictions"):
        matrix = info.setdefault("predictor_matrix", {})
        for record in info["predictions"]:
            matrix.setdefault(record["branch"], []).append(record.get("status", "unresolved"))
        for branch, statuses in matrix.items():
            matrix[branch] = list(dict.fromkeys(statuses))


def source_record(row, source, node_id, target):
    if metadata(row):
        if metadata(row)["target_id"] != target["target_id"]:
            raise ValueError("Candidate collections refer to different targets.")
        _attach_campaign_evidence(row)
        return row
    raw = row.setdefault("raw_metadata", {})
    invocation = raw.get("candidate_execution", {})
    native_id = invocation.get("native_candidate_id") if invocation.get("unit_kind") == "native_invocation" else row["candidate_id"]
    native_id = native_id or row["candidate_id"]
    campaign = raw.get("distributed_campaign", {})
    origin = raw.get("retained_origin", {})
    run_id = origin.get("run_id") or source.get("source_run_id") or Path(source["run_dir"]).name
    backbone_native = raw.get("backbone_id") or row.get("metrics", {}).get("native_design_index")
    backbone_native = backbone_native or raw.get("source_candidate", {}).get("candidate_id") or native_id
    backbone = identity(target["target_id"], row["source_tool"], run_id, backbone_native)
    variant = identity(backbone, native_id, row.get("binder_sequence"))
    row["candidate_id"] = "cs_" + variant
    raw["candidate_set"] = dict(schema_version=SCHEMA, target_id=target["target_id"],
        generator=row["source_tool"], source_run=run_id, source_run_dir=source["run_dir"],
        source_node=node_id, native_candidate_id=native_id,
        native_attempt_id=campaign.get("attempt_id") or raw.get("trajectory_id"),
        backbone_id=backbone, variant_id=variant, parent_variant_ids=[], record_kind="variant",
        native_backbone_identity=backbone_native, native_parents=invocation.get("native_parents", list(row.get("parents", []))),
        original_sequence=row.get("binder_sequence"), original_structure=row.get("complex_pdb") or row.get("binder_pdb"),
        native_filter=native_filter(row), sequence_state=sequence_state(row),
        native_measurements=copy.deepcopy(row.get("metrics", {})),
        artifacts=artifact_evidence(row), source_provenance=origin or source.get("source_provenance", {}))
    _attach_campaign_evidence(row)
    if (source.get("review_candidates_jsonl") or row.get("metrics", {}).get("candidate_pool_level")) and row.get("metrics", {}).get("candidate_pool_level") != "native_pipeline":
        raw["candidate_set"]["native_output_role"] = "intermediate"
    return row


def write_collection(parent, node_id, rows, target, *, details=None):
    destination = Path(parent) / "artifacts/pipeline_nodes" / node_id
    normalized = destination / "artifacts/normalized_candidates"
    normalized.mkdir(parents=True, exist_ok=True)
    path = normalized / "candidates.jsonl"
    temporary = path.with_suffix(".tmp")
    counts = Counter()
    with temporary.open("w") as handle:
        for row in rows:
            handle.write(json.dumps(json_safe(normalize_candidate(row)), sort_keys=True, allow_nan=False) + "\n")
            counts["normalized_candidates"] += 1
            info = metadata(row)
            counts[info.get("sequence_state", sequence_state(row))] += 1
            counts["native_rejected_candidates"] += info.get("native_filter") == "failed"
            counts["sequence_variants"] += bool(info.get("parent_variant_ids"))
            counts["completed_predictions"] += sum(p.get("status") == "completed" for p in info.get("predictions", [])) + (info.get("record_kind") == "prediction")
            counts["analyzed_candidates"] += "analysis" in info
            counts["selected_candidates"] += row.get("raw_metadata", {}).get("selection", {}).get("selected") is True
    temporary.replace(path)
    manifest = destination / "artifacts/candidate_set/manifest.json"
    write_json(manifest, dict(schema_version=SCHEMA, candidates_jsonl=str(path),
                             sha256=file_digest(path), counts=dict(counts), target=target,
                             **(details or {})))
    write_json(destination / "metadata.json", dict(status="completed", tool="pipeline_collection", run_id=node_id))
    write_json(destination / "input.json", {"inputs": {}, "params": {}})
    write_json(destination / "result.json", {"success": True, "outputs": {"candidates_jsonl": str(path)}})
    return dict(run_dir=str(destination), candidates_jsonl=str(path), candidate_count=counts["normalized_candidates"],
                target=target, candidate_set_manifest=str(manifest), stage_counts=dict(counts))


def native_counts(provenance):
    metrics = provenance.get("result_metrics", {})
    progress = provenance.get("metadata", {}).get("campaign_progress", {})
    counts = {}
    for name, choices in {"native_attempts": [progress.get("attempts_started"), metrics.get("attempts_started"), metrics.get("trajectory_attempt_count"), metrics.get("native_attempts")],
                          "native_rejected_attempts": [progress.get("native_rejected_trajectories"), metrics.get("native_rejected_trajectories"), metrics.get("native_rejected_attempts")],
                          "raw_designs": [metrics.get("candidate_count"), metrics.get("raw_designs")]}.items():
        value = next((v for v in choices if type(v) is int), None)
        if value is not None:
            counts[name] = value
    return counts


def merge(parent, node_id, sources, params):
    target = next((s["target"] for _, s in sources if s.get("target")), {})
    target_info = target.get("identity") or target_identity(target)
    policy = params.get("deduplication", "mark_duplicates")
    with index() as db:
        db.execute("CREATE TABLE records (id TEXT PRIMARY KEY, seq TEXT, payload TEXT)")
        for branch, source in sources:
            other = source.get("target", {}).get("identity")
            if other and other["target_id"] != target_info["target_id"]:
                raise ValueError("Cannot merge candidate sets from different targets.")
            for row in iter_candidates(source, review=True):
                row = source_record(row, source, branch, target_info)
                seq = row.get("binder_sequence") if sequence_state(row) == "sequence_ready" else None
                payload = json.dumps(row, sort_keys=True)
                previous = db.execute("SELECT payload FROM records WHERE id=?", (row["candidate_id"],)).fetchone()
                if previous and previous[0] != payload:
                    raise ValueError("Conflicting records for the same scientific candidate identity.")
                db.execute("INSERT OR IGNORE INTO records VALUES (?,?,?)", (row["candidate_id"], seq, payload))
        db.execute("CREATE INDEX sequence_index ON records(seq)")
        def merged():
            for key, seq, payload in db.execute("SELECT * FROM records ORDER BY id"):
                row = json.loads(payload)
                if seq and policy == "mark_duplicates":
                    count = db.execute("SELECT count(*) FROM records WHERE seq=?", (seq,)).fetchone()[0]
                    metadata(row)["deduplication"] = dict(policy=policy, group=identity(target_info["target_id"], seq), members=count)
                yield row
        result = write_collection(parent, node_id, merged(), target, details=dict(deduplication=policy,
            sources=[dict(node=branch, manifest=source.get("candidate_set_manifest"),
                          run_dir=source["run_dir"], provenance=source.get("source_provenance", {})) for branch, source in sources]))
        totals = Counter()
        seen_runs = set()
        for _, source in sources:
            if source.get("source_run_id") not in seen_runs:
                totals.update(native_counts(source.get("source_provenance", {})))
                seen_runs.add(source.get("source_run_id"))
        result["stage_counts"].update(totals)
    return result


def validate_settings(kind, params):
    if kind == "candidate_merge":
        if params.get("deduplication") not in {"keep_all", "mark_duplicates"}:
            raise ValueError("Choose keep_all or mark_duplicates; histories are never deleted.")
    if kind in {"shared_filter", "candidate_rank"}:
        for rule in params.get("rules", []):
            if not isinstance(rule, dict) or set(rule) != {"field", "operator", "value"} or not isinstance(rule["field"], str) or not rule["field"] or rule["operator"] not in COMPARISONS or finite_number(rule["value"]) is None:
                raise ValueError("Each scientific rule needs a field, comparison operator and finite value.")
    if kind == "shared_filter":
        from mn_molecular_design_studio.validation_acceptance import normalize_policy
        normalize_policy({"binder_format": params.get("binder_format", "miniprotein"),
                          "target_type": params.get("target_type", "structured"),
                          "mandatory_filters": params.get("rules", []),
                          "informational_scores": params.get("informational_scores", []),
                          "ranking_metric": params.get("ranking_metric", "")})
    if kind == "candidate_rank":
        if len(params.get("ranking_rules", [])) > 16:
            raise ValueError("At most sixteen ordered ranking measurements are supported.")
        for rule in params.get("ranking_rules", []):
            if not isinstance(rule, dict) or set(rule) != {"field", "direction"} or not isinstance(rule["field"], str) or not rule["field"] or rule["direction"] not in {"higher", "lower"}:
                raise ValueError("Ranking rules need an exact field and higher/lower direction.")


def observation(row, field, output=None):
    from mn_molecular_design_studio.output_vocabulary import metric_value, candidate_output
    value = metric_value(row, field)
    if value is None and str(field).startswith(("prediction/", "analysis/")):
        parts = str(field).split("/")
        if len(parts) not in {3, 4}:
            raise ValueError("Use prediction/<branch>/<metric> or analysis/<branch>/<metric>, optionally naming a result ID.")
        family, branch, *rest = parts
        evidence_key = "predictions" if family == "prediction" else "analysis_results"
        records = [item for item in metadata(row).get(evidence_key, [])
                   if item.get("branch") == branch and item.get("status") == "completed"
                   and (len(rest) == 1 or (item.get("prediction_id") or item.get("analysis_result_id")) == rest[0])]
        if len(records) == 1:
            value = (records[0].get("candidate") or {}).get("metrics", {}).get(rest[-1])
    if value is None:
        return dict(status="missing_evidence", value=None)
    output = candidate_output(row) if output is None else output
    if any(issue["field"] == field for issue in output.get("issues", [])) or any(
        any(alias["source_field"] == field for alias in aliases) for aliases in output.get("conflicts", {}).values()
    ):
        return dict(status="invalid_metadata", value=None, raw_value=value)
    number = finite_number(value)
    low, high = RANGES.get(field.split("/")[-1], (None, None))
    if number is None or low is not None and number < low or high is not None and number > high:
        return dict(status="invalid_metadata", value=None, raw_value=value)
    return dict(status="valid", value=number)


def shared_filter(parent, node_id, source, params):
    from mn_molecular_design_studio.validation_acceptance import evaluate_candidate, normalize_policy
    policy = normalize_policy({"binder_format": params.get("binder_format", "miniprotein"),
                               "target_type": params.get("target_type", "structured"),
                               "mandatory_filters": params.get("rules", []),
                               "informational_scores": params.get("informational_scores", []),
                               "ranking_metric": params.get("ranking_metric", "")})
    acceptance_enabled = bool(policy["mandatory_filters"] or policy["informational_scores"] or policy["ranking_metric"])
    review_path = Path(parent) / "artifacts/pipeline_nodes" / node_id / "artifacts/candidate_set/all_candidates.jsonl"
    review_path.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    def passed():
        with review_path.open("w") as review:
            for row in iter_candidates(source):
                info = metadata(row)
                if not info:
                    raise ValueError("Shared filters require a CandidateSet merge first.")
                decisions = []
                retained = (row.get("raw_metadata") or {}).get("validation_acceptance") or {}
                provenance = (row.get("raw_metadata") or {}).get("validation_acceptance_provenance") or {}
                reuse_decision = (provenance.get("acceptance_node") == node_id and retained.get("policy") == policy)
                acceptance = (copy.deepcopy(retained) if reuse_decision else evaluate_candidate(row, policy)) if acceptance_enabled else None
                if acceptance:
                    row.setdefault("raw_metadata", {})["validation_acceptance"] = acceptance
                    for rule in policy["mandatory_filters"]:
                        measured = acceptance["measurements"][rule["field"]]
                        status = ("passed" if measured["passed"] else "failed") if measured["status"] == "measured" else measured["status"]
                        decisions.append(dict(rule=rule, observation=measured, status=status))
                if not (row.get("complex_pdb") or row.get("binder_pdb")) or any(a.get("status") != "present" for a in artifact_evidence(row).values()):
                    status = "missing_evidence"
                elif info.get("native_output_role") == "intermediate":
                    status = "not_evaluated"
                elif info.get("native_filter") == "failed" and params["exclude_native_rejected"]:
                    status = "failed"
                elif info.get("native_filter") == "invalid_metadata":
                    status = "invalid_metadata"
                elif sequence_state(row) != "sequence_ready":
                    status = "unresolved"
                elif acceptance and acceptance["status"] == "rejected":
                    status = "failed"
                elif acceptance and acceptance["status"] == "unresolved":
                    status = next((v for v in ("invalid_metadata", "missing_evidence")
                                   if any(d["status"] == v for d in decisions)), "unresolved")
                elif acceptance and acceptance["status"] == "accepted":
                    status = "passed"
                else:
                    status = next((v for v in ("invalid_metadata", "missing_evidence", "failed") if any(d["status"] == v for d in decisions)), "passed" if decisions else "not_evaluated")
                info["shared_filter"] = dict(status=status, decisions=decisions, policy=copy.deepcopy(params),
                                              acceptance_status=acceptance["status"] if acceptance else "not_evaluated")
                info["eligible_for_prediction"] = status in {"passed", "not_evaluated"} and info.get("native_output_role") != "intermediate"
                counts[status] += 1
                review.write(json.dumps(row, sort_keys=True) + "\n")
                if info["eligible_for_prediction"]:
                    yield row
    result = write_collection(parent, node_id, passed(), source["target"])
    result["review_candidates_jsonl"] = str(review_path)
    result["stage_counts"].update(candidates_filtered=counts["failed"], **{"filter_" + k: v for k, v in counts.items()})
    return result


def sequence_inputs(parent, node_id, source, params):
    """Only compatible records enter native MPNN; every original stays in fan-in."""
    initial = params.get("sequence_stage") == "initial"
    disabled = params.get("sequence_stage") == "soluble_optimization" and not params.get("optimization_enabled")
    def eligible():
        if disabled:
            return
        for row in iter_candidates(source):
            info = metadata(row)
            if info.get("native_filter") in {"failed", "invalid_metadata"} or info.get("native_output_role") == "intermediate":
                continue
            state = sequence_state(row)
            allowed_states = {"backbone_only"} if initial else {"sequence_ready"} if params.get("sequence_stage") == "soluble_optimization" else {"backbone_only", "sequence_ready"}
            if state not in allowed_states:
                continue
            if not row.get("binder_chains") or not (row.get("complex_pdb") or row.get("binder_pdb")):
                continue
            if any(item["status"] != "present" for item in artifact_evidence(row).values()):
                continue
            yield row
    return write_collection(parent, node_id + "_sequence_inputs", eligible(), source["target"])


def rank(parent, node_id, source, params):
    """Ordered measurements within recorded predictor protocols, then exact diversity.

    No implicit score averaging or cross-engine confidence equivalence. Unknown
    protocols are singleton cohorts, as in the existing scientific selector.
    """
    from mn_molecular_design_studio.candidate_selection import _protocol
    rules = params["ranking_rules"]
    with index() as db:
        columns = ",".join(f"r{i} REAL" for i in range(len(rules)))
        db.execute("CREATE TABLE ranked (id TEXT PRIMARY KEY, cohort TEXT, status TEXT, seq TEXT, payload TEXT" + ("," + columns if columns else "") + ")")
        for row in iter_candidates(source, review=True):
            info = metadata(row)
            observed = {rule["field"]: observation(row, rule["field"]) for rule in [*params["rules"], *rules]}
            reasons, protocols = [], []
            status = "passed" if observed else "not_evaluated"
            if not info.get("eligible_for_prediction", False):
                status = info.get("shared_filter", {}).get("status", "unresolved")
                reasons.append("Candidate was not eligible after shared filtering.")
            for field, value in observed.items():
                if value["status"] != "valid":
                    status = value["status"]
                    reasons.append(field + ": " + value["status"])
                if field.startswith(("prediction/", "analysis/")):
                    parts = field.split("/")
                    evidence_key = "predictions" if parts[0] == "prediction" else "analysis_results"
                    candidates = [p["candidate"] for p in info.get(evidence_key, []) if p["branch"] == parts[1]
                                  and (len(parts) == 3 or p["prediction_id"] == parts[2]) and p["status"] == "completed"]
                    protocols.append([_protocol(p) for p in candidates])
                else:
                    protocols.append(_protocol(row))
            for rule in params["rules"]:
                value = observed[rule["field"]]["value"]
                if value is not None and not COMPARISONS[rule["operator"]](value, rule["value"]):
                    if status in {"passed", "not_evaluated"}:
                        status = "failed"
                    reasons.append("Configured scientific threshold failed: " + rule["field"])
            cohort = identity(info["target_id"], protocols)
            decision = dict(status=status, selected=False, reasons=reasons, rank=None, cohort=cohort,
                            ranking_metric=", ".join(r["field"] for r in rules), direction="ordered_rules",
                            measurements=observed, policy=params)
            row["raw_metadata"]["selection"] = decision
            values = [observed[r["field"]]["value"] for r in rules]
            values = [v * (-1 if r["direction"] == "higher" else 1) if v is not None else None for r, v in zip(rules, values)]
            payload = [row["candidate_id"], cohort, status, row.get("binder_sequence"), json.dumps(row), *values]
            db.execute("INSERT INTO ranked VALUES (" + ",".join("?" for _ in payload) + ")", payload)
        db.execute("CREATE TABLE diversity (cohort TEXT, seq TEXT, PRIMARY KEY(cohort,seq))")
        ordering = "cohort," + ",".join([*(f"r{i}" for i in range(len(rules))), "id"])
        def ordered():
            previous, position, kept = None, 0, 0
            for _, cohort, status, seq, payload in db.execute("SELECT id,cohort,status,seq,payload FROM ranked ORDER BY " + ordering):
                if cohort != previous:
                    previous, position, kept = cohort, 0, 0
                row = json.loads(payload)
                decision = row["raw_metadata"]["selection"]
                if status == "passed":
                    position += 1
                    decision["rank"] = position if rules else None
                    duplicate = params["unique_sequences"] and db.execute("SELECT 1 FROM diversity WHERE cohort=? AND seq=?", (cohort, seq)).fetchone()
                    if duplicate:
                        decision["reasons"].append("Exact-sequence diversity representative already selected; history retained.")
                    elif params["limit"] and kept >= params["limit"]:
                        decision["reasons"].append("Configured shortlist limit within comparable protocol.")
                    else:
                        decision["selected"] = True
                        kept += 1
                        db.execute("INSERT OR IGNORE INTO diversity VALUES (?,?)", (cohort, seq))
                metadata(row)["selection_id"] = identity(node_id, metadata(row)["variant_id"], params)
                yield row
        return write_collection(parent, node_id, ordered(), source["target"], details=dict(selection_policy=params))


def export_collection(parent, source, graph, node_states):
    """Self-contained scientific export inside the existing run ZIP artifact tree."""
    import shutil
    destination = Path(parent) / "artifacts/candidate_set"
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source["candidates_jsonl"], destination / "candidates.jsonl")
    write_json(destination / "graph.json", graph)
    write_json(destination / "execution.json", node_states)
    write_json(destination / "manifest.json", read_json(Path(source["candidate_set_manifest"])))
    from mn_molecular_design_studio.workflow_modes import describe, observed_execution_plan
    summary = describe(graph)
    write_json(destination / "workflow_summary.json", dict(**summary,
        stage_counts={key: value.get('stage_counts', {}) for key, value in node_states.items()},
        execution_plan=observed_execution_plan(graph, node_states)))
    # Keep empty native outcomes and rejection provenance in the export too.
    write_json(destination / "native_sources.json", [dict(operation=operation,
        output=node_states.get(operation['node'], {}).get('output', {})) for operation in summary['generator_operations']])
    selected_files = {}
    with (destination / "candidate_index.csv").open("w") as handle, (destination / "predictions.jsonl").open("w") as predictions, (destination / "selected_sequences.fasta").open("w") as fasta, (destination / "native_sequences.fasta").open("w") as native_fasta:
        writer = csv.DictWriter(handle, fieldnames=["candidate_id", "generator", "source_run", "native_attempt_id", "backbone_id", "variant_id", "parent_variant_ids", "sequence_state", "native_filter", "shared_filter", "validation_acceptance", "acceptance_reasons", "predictor_matrix", "rank", "selected"])
        writer.writeheader()
        for row in iter_candidates(source, review=True):
            info = metadata(row)
            selection = row.get("raw_metadata", {}).get("selection", {})
            index_row = {key: info.get(key) for key in writer.fieldnames
                         if key not in {"candidate_id", "shared_filter", "validation_acceptance",
                                        "acceptance_reasons", "rank", "selected"}}
            index_row.update(
                candidate_id=row["candidate_id"],
                shared_filter=info.get("shared_filter", {}).get("status"),
                validation_acceptance=(row.get("raw_metadata", {}).get("validation_acceptance") or {}).get("status"),
                acceptance_reasons="; ".join((row.get("raw_metadata", {}).get("validation_acceptance") or {}).get("reasons", [])),
                rank=selection.get("rank"), selected=selection.get("selected", False),
            )
            writer.writerow(index_row)
            for prediction in info.get("predictions", []):
                predictions.write(json.dumps(dict(variant_id=info["variant_id"], **prediction), sort_keys=True) + "\n")
            if selection.get("selected") or summary['effective_mode'] == 'native_complete':
                if row.get('binder_sequence'):
                    (native_fasta if summary["effective_mode"] == "native_complete" else fasta).write(">" + row["candidate_id"] + "\n" + row["binder_sequence"] + "\n")
                structures = [row, *(p["candidate"] for p in info.get("predictions", []) if p["status"] == "completed")]
                for structure in structures:
                    for key in ("complex_pdb", "binder_pdb", "target_pdb"):
                        if not structure.get(key):
                            continue
                        path = Path(structure[key])
                        if not path.is_file() and summary["effective_mode"] == "native_complete":
                            continue  # Missing native artifacts remain explicit evidence, never a new failure.
                        if not path.is_file():
                            raise ValueError("Selected artifact is unavailable: " + str(path))
                        digest = file_digest(path)
                        local = destination / "structures" / (digest + path.suffix)
                        local.parent.mkdir(exist_ok=True)
                        if not local.exists():
                            shutil.copyfile(path, local)
                        selected_files[str(path)] = dict(path=str(local.relative_to(destination)), sha256=digest)
    write_json(destination / ("native_artifacts.json" if summary["effective_mode"] == "native_complete" else "selected_artifacts.json"), selected_files)
    with (destination / "analysis_results.jsonl").open("w") as analyses:
        for row in iter_candidates(source, review=True):
            analyses.write(json.dumps(dict(candidate_id=row["candidate_id"], variant_id=metadata(row)["variant_id"],
                analysis=metadata(row).get("analysis"), analysis_results=metadata(row).get("analysis_results", []),
                selection=row.get("raw_metadata", {}).get("selection")), sort_keys=True) + "\n")
    write_json(destination / "checksums.json", {str(p.relative_to(destination)): file_digest(p) for p in destination.rglob("*") if p.is_file() and p.name != "checksums.json"})
    return destination


def attach_stage(parent, node_id, output, source, *, category, scientific_key, settings, evidence_namespace=None):
    """Carry stable variant identities through existing MPNN/prediction adapters."""
    if not source.get("candidate_set_manifest") or output.get("candidate_set_manifest") and not output.get("collection_stage_pending"):
        return output
    with index() as db:
        db.execute("CREATE TABLE parents (id TEXT PRIMARY KEY, payload TEXT)")
        for row in iter_candidates(source):
            db.execute("INSERT INTO parents VALUES (?,?)", (row["candidate_id"], json.dumps(row)))
        def transformed():
            if category == "Sequence design":
                yield from iter_candidates(source)
            for row in iter_candidates(output):
                ids = [row["candidate_id"], *row.get("parents", [])]
                direct = row.get("raw_metadata", {}).get("source_candidate", {})
                if direct.get("candidate_id"):
                    ids.append(direct["candidate_id"])
                matches = {key: json.loads(found[0]) for key in ids
                           if (found := db.execute("SELECT payload FROM parents WHERE id=?", (key,)).fetchone())}
                if len(matches) != 1:
                    raise ValueError("Collection stage lost or ambiguously mapped its input candidate.")
                original = next(iter(matches.values()))
                info = copy.deepcopy(metadata(original))
                native_id = row["candidate_id"]
                if category == "Sequence design" and row["candidate_id"] == original["candidate_id"]:
                    continue
                if category == "Sequence design":
                    info.update(parent_variant_ids=[info["variant_id"]],
                        variant_id=identity(info["variant_id"], scientific_key, native_id, row.get("binder_sequence")),
                        sequence_state=sequence_state(row), record_kind="variant",
                        sequence_design=dict(node=node_id, settings=settings, native_candidate_id=native_id))
                    # Never carry validation or original native filtering onto a child sequence.
                    for key in ("predictions", "analysis", "shared_filter", "eligible_for_prediction"):
                        info.pop(key, None)
                    info["native_filter"] = "not_evaluated"
                elif category == "Structure prediction":
                    if row.get("binder_sequence") != original.get("binder_sequence"):
                        raise ValueError("Prediction changed the source sequence variant.")
                    info.update(record_kind="prediction", prediction_id=identity(node_id, scientific_key, info["variant_id"], native_id),
                                predictor_execution_node=node_id,
                                predictor_node=evidence_namespace or node_id, predictor_settings=settings,
                                prediction_execution=output.get("source_provenance", {}))
                    row.setdefault("raw_metadata", {})["pipeline_refold_input"] = dict(
                        complex_pdb=original.get("complex_pdb") or original.get("binder_pdb"),
                        binder_chains=original.get("binder_chains", []), target_chains=original.get("target_chains", []),
                        source_candidate_id=original["candidate_id"])
                else:
                    info["record_kind"] = "measurement"
                    info["analysis_result_id"] = identity(node_id, scientific_key, info["variant_id"], native_id)
                    info["measurement_node"] = node_id
                row.setdefault("raw_metadata", {})["candidate_set"] = info
                yield row
        result = write_collection(parent, node_id + "_collection", transformed(), source["target"])
    for key in ("candidate_execution", "partial_pool"):
        if key in output:
            result[key] = output[key]
    return result


def analyze(parent, node_id, sources, params, expected_predictors):
    """Join by variant identity, retaining an incomplete predictor matrix and raw evidence."""
    target = next(s["target"] for _, s in sources if s.get("target"))
    with index() as db:
        db.execute("CREATE TABLE variants (id TEXT PRIMARY KEY, payload TEXT)")
        db.execute("CREATE TABLE predictions (id TEXT PRIMARY KEY, variant TEXT, branch TEXT, payload TEXT)")
        for branch, source in sources:
            for row in iter_candidates(source, review=True):
                info = metadata(row)
                if not info:
                    raise ValueError("Analyzer requires CandidateSet lineage on every record.")
                if info["target_id"] != target["identity"]["target_id"]:
                    raise ValueError("Prediction refers to a different target.")
                if info.get("record_kind") in {"prediction", "measurement"}:
                    key = info.get("analysis_result_id") if info.get("record_kind") == "measurement" else info.get("prediction_id")
                    key = key or identity(branch, row["candidate_id"])
                    payload = json.dumps(row, sort_keys=True)
                    old = db.execute("SELECT payload FROM predictions WHERE id=?", (key,)).fetchone()
                    if old and old[0] != payload:
                        raise ValueError("Conflicting prediction receipt for one prediction identity.")
                    db.execute("INSERT OR IGNORE INTO predictions VALUES (?,?,?,?)", (key, info["variant_id"], info.get("measurement_node") or info.get("predictor_node", branch), payload))
                else:
                    # A shared-filter decision supersedes the earlier unfiltered view.
                    old = db.execute("SELECT payload FROM variants WHERE id=?", (info["variant_id"],)).fetchone()
                    if old and json.loads(old[0]).get("binder_sequence") != row.get("binder_sequence"):
                        raise ValueError("Conflicting sequences for one sequence variant identity.")
                    if not old or "shared_filter" in info:
                        db.execute("INSERT OR REPLACE INTO variants VALUES (?,?)", (info["variant_id"], json.dumps(row)))
                    # Controller-run validation evidence is already a completed
                    # prediction attached to the candidate. Make it available
                    # to the same analyzer fan-in as graph prediction nodes.
                    for evidence_key in ("predictions", "analysis_results"):
                        for record in info.get(evidence_key, []):
                            evidence_candidate = record.get("candidate")
                            if not isinstance(evidence_candidate, dict):
                                continue
                            evidence_id = record.get("prediction_id") or record.get("analysis_result_id")
                            if not evidence_id:
                                continue
                            prediction_payload = json.dumps(evidence_candidate, sort_keys=True)
                            previous = db.execute("SELECT payload FROM predictions WHERE id=?", (evidence_id,)).fetchone()
                            if previous and previous[0] != prediction_payload:
                                raise ValueError("Conflicting attached controller evidence for one prediction identity.")
                            db.execute("INSERT OR IGNORE INTO predictions VALUES (?,?,?,?)", (
                                evidence_id, info["variant_id"], record.get("branch") or branch, prediction_payload,
                            ))
        if db.execute("SELECT count(*) FROM predictions p LEFT JOIN variants v ON p.variant=v.id WHERE v.id IS NULL").fetchone()[0]:
            raise ValueError("Prediction has no matching sequence variant in the analyzer inputs.")
        db.execute("CREATE INDEX prediction_variants ON predictions(variant)")
        def joined():
            for variant, payload in db.execute("SELECT * FROM variants ORDER BY id"):
                row = json.loads(payload)
                info = metadata(row)
                predictions, analyses, measurements, matrix = [], [], {}, {}
                for prediction_id, branch, pred_json in db.execute("SELECT id,branch,payload FROM predictions WHERE variant=? ORDER BY id", (variant,)):
                    pred = json.loads(pred_json)
                    is_analysis = metadata(pred).get("record_kind") == "measurement"
                    if params.get("calculate_self_consistency") and not is_analysis and pred.get("raw_metadata", {}).get("pipeline_refold_input"):
                        from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import _add_target_aligned_pose_rmsd
                        pred = _add_target_aligned_pose_rmsd(pred, Path(parent))
                    if pred.get("binder_sequence") != row.get("binder_sequence"):
                        raise ValueError("Analyzer refuses cross-sequence prediction evidence.")
                    from mn_molecular_design_studio.output_vocabulary import candidate_evidence
                    status = candidate_evidence(pred)[0]["status"]
                    if status in {"error", "failed", "invalid", "missing"}:
                        status = "failed" if status in {"error", "failed"} else "unresolved"
                    else:
                        status = "completed" if Path(pred.get("complex_pdb") or "").is_file() else "missing_evidence"
                    from mn_molecular_design_studio.output_vocabulary import candidate_output
                    output_contract = candidate_output(pred)
                    record = dict(prediction_id=prediction_id, branch=branch, predictor=pred["source_tool"],
                                  status=status, candidate=pred, measurements={k: observation(pred, k, output_contract) for k in pred.get("metrics", {})
                                                                            if isinstance(pred["metrics"][k], (int, float))})
                    (analyses if is_analysis else predictions).append(record)
                    if not is_analysis:
                        matrix.setdefault(branch, []).append(status)
                    # Sample-qualified keys avoid averaging different seeds or model outputs.
                    for metric, observed in record["measurements"].items():
                        if status == "completed" and observed["status"] == "valid":
                            prefix = "analysis" if is_analysis else "prediction"
                            measurements[f"{prefix}/{branch}/{prediction_id}/{metric}"] = observed["value"]
                for branch in expected_predictors:
                    branch_info = expected_predictors[branch]
                    if branch not in matrix:
                        matrix[branch] = ["disabled" if not branch_info.get("enabled", True) else
                            "failed" if branch_info.get("status") in {"failed", "blocked"} else
                            "missing_evidence" if info.get("eligible_for_prediction", False) else "not_requested"]
                    elif not branch_info.get("retained") and info.get("eligible_for_prediction", False):
                        settings = branch_info.get("params", {})
                        expected = int(settings.get("samples", settings.get("diffusion_batch_size", settings.get("model_count", 1)))) * int(settings.get("seed_count", 1))
                        missing = max(0, expected - len(matrix[branch]))
                        matrix[branch].extend(["missing_evidence"] * missing)
                # An unambiguous branch/metric alias is convenient for configured ranking.
                for prefix, evidence in (("prediction", predictions), ("analysis", analyses)):
                    for branch in {record["branch"] for record in evidence}:
                        records = [p for p in evidence if p["branch"] == branch and p["status"] == "completed"]
                        if len(records) == 1:
                            for metric, observed in records[0]["measurements"].items():
                                if observed["status"] == "valid":
                                    measurements[f"{prefix}/{branch}/{metric}"] = observed["value"]
                row["metrics"].update(measurements)
                info.update(predictions=predictions, analysis_results=analyses, predictor_matrix=matrix,
                    analysis=dict(analysis_id=identity(node_id, variant, params, [p["prediction_id"] for p in [*predictions, *analyses]]),
                                  measurements=measurements, settings=params,
                                  unavailable_measurements="Only supplied valid measurements are collected; no synthetic scores."))
                # Existing selection/export and evidence tables understand this additive list.
                row["raw_metadata"]["evidence_candidates"] = [p["candidate"] for p in [*predictions, *analyses]]
                yield row
        result = write_collection(parent, node_id, joined(), target, details=dict(expected_predictors=expected_predictors))
    return result
