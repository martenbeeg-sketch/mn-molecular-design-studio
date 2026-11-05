"""Opt-in, auditable selection over existing candidate collections.

Thresholds are supplied by the task. Ranking is a single declared metric within
an identical recorded protocol; no composite scientific score is constructed.
"""
from __future__ import annotations

import copy
import hashlib
import json
import operator
from pathlib import Path
import shutil

from mn_molecular_design_studio.measurement_status import finite_number, RANGES
from mn_molecular_design_studio.output_vocabulary import DEFINITIONS, candidate_evidence, candidate_output, metric_value

OPERATORS = {">=": operator.ge, "<=": operator.le, ">": operator.gt, "<": operator.lt}
RULE_OPERATORS = {*OPERATORS, "between"}


def selection_policy(settings):
    policy = copy.deepcopy(settings)
    rules = policy.setdefault("rules", [])
    if not isinstance(rules, list):
        raise ValueError("Selection rules must be a list.")
    metric = str(policy.get("metric") or "").strip()
    if metric:
        primary = {"field": metric, "operator": policy.get("operator", ">="), "value": policy.get("threshold")}
        if primary["operator"] == "between":
            primary["upper_value"] = policy.get("threshold_max")
        rules.insert(0, primary)
    for rule in rules:
        if not isinstance(rule, dict) or not rule.get("field") or rule.get("operator") not in RULE_OPERATORS or finite_number(rule.get("value")) is None:
            raise ValueError("Each selection rule needs a field, comparison operator and finite numeric value.")
        if rule["operator"] == "between":
            upper = finite_number(rule.get("upper_value", rule.get("max_value")))
            if upper is None or upper < float(rule["value"]):
                raise ValueError("A between rule needs an upper_value greater than or equal to its lower value.")
    for key in ("required_metrics", "manual_ids"):
        if not isinstance(policy.setdefault(key, []), list) or not all(isinstance(value, str) and value for value in policy[key]):
            raise ValueError(f"{key} must be a list of nonempty strings.")
    policy["ranking_metric"] = policy.get("ranking_metric") or metric
    ranking_field = policy["ranking_metric"].split("/")[-1]
    definition = DEFINITIONS.get(ranking_field.removeprefix("common.")) if ranking_field.startswith("common.") else None
    known_direction = definition.direction if definition else ""
    configured_direction = policy.get("direction")
    if not configured_direction and "higher_is_better" in policy:
        configured_direction = "higher" if policy["higher_is_better"] else "lower"
    policy["direction"] = configured_direction or known_direction or "higher"
    if policy["direction"] not in {"higher", "lower"}:
        raise ValueError("Ranking direction must be higher or lower.")
    if known_direction and policy["direction"] != known_direction:
        raise ValueError(f"The declared output contract ranks {ranking_field} with {known_direction} values first.")
    if policy["ranking_metric"] == "analysis_score":
        raise ValueError("Choose an interpretable measurement instead of the legacy composite analysis_score.")
    logic = str(policy.get("rule_logic", "all")).casefold()
    if logic not in {"all", "any"}:
        raise ValueError("Rule logic must be all or any.")
    policy["rule_logic"] = logic
    quantile = finite_number(policy.get("quantile", 0.0))
    if quantile is None or quantile < 0 or quantile > 1:
        raise ValueError("Quantile must be between 0 (disabled) and 1.")
    policy["quantile"] = quantile
    if quantile and not policy["ranking_metric"]:
        raise ValueError("Choose a ranking measurement before enabling quantile selection.")
    for key, default in (("limit", 0), ("max_per_cluster", 1)):
        value = policy.setdefault(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or value < (1 if key == "max_per_cluster" else 0):
            raise ValueError(f"{key} must be a nonnegative integer (max_per_cluster at least 1).")
    return policy


def _observation(candidate, field):
    """Resolve the exact source for validity, ranking and protocol checks."""
    if field.startswith("evidence/"):
        parts = field.split("/")
        if len(parts) not in {4, 5}:
            raise ValueError("Invalid evidence selector.")
        records = [record for record in candidate_evidence(candidate) if record["same_sequence"]
                   and record["stage"] == parts[1] and record["source_tool"] == parts[2]
                   and (len(parts) == 4 or record["candidate_id"] == parts[3])]
        if len(records) != 1:
            return None, "missing", "Measurement source is missing or ambiguous", None
        record = records[0]
        item, native_field = record["candidate"], parts[-1]
    else:
        record = candidate_evidence(candidate)[0]
        item, native_field = candidate, field
    if record["status"] in {"error", "failed", "missing", "invalid"}:
        return None, "failed", record["error"] or record["status"], item
    normalized = candidate_output(item)
    if any(issue["field"] == native_field for issue in normalized.get("issues", [])):
        return None, "missing", "Measurement violates its declared output contract", item
    if any(any(entry["source_field"] == native_field for entry in entries) for entries in normalized.get("conflicts", {}).values()):
        return None, "missing", "Conflicting measurement aliases", item
    value = finite_number(metric_value(item, native_field))
    if value is None:
        return None, "missing", "Measurement is missing, invalid or nonfinite", item
    low, high = RANGES.get(native_field, (None, None))
    if low is not None and value < low or high is not None and value > high:
        return None, "missing", "Measurement is outside its physical output domain", item
    return value, "valid", "", item


def _protocol(item):
    raw = item.get("raw_metadata") or {}
    parameters = raw.get("analysis_parameters") or raw.get("prediction_protocol")
    if not parameters:
        # Unknown protocols are not evidence of comparability.
        return {"unrecorded_protocol_candidate": item.get("candidate_id")}
    # PAE filenames identify individual observations, not different protocols.
    parameters = {key: value for key, value in parameters.items() if key not in {"pae_file", "workers", "nprocs"}}
    return {"stage": item.get("stage"), "tool": item.get("source_tool"), "parameters": parameters,
            "prediction_protocol": raw.get("prediction_protocol"), "prediction_engine": raw.get("engine"),
            "image": raw.get("image_digest") or raw.get("image"),
            "target": raw.get("target_structure_digest") or raw.get("target_sequences") or
                      item.get("target_pdb"),
            "target_chains": item.get("target_chains")}


def _digest(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def evaluate_selection(candidates, settings):
    policy = selection_policy(settings)
    required = list(dict.fromkeys([*policy["required_metrics"], *[rule["field"] for rule in policy["rules"]],
                                  *([policy["ranking_metric"]] if policy["ranking_metric"] else []),
                                  *([policy["diversity_field"]] if policy.get("diversity_field") else [])]))
    explicit = bool(required or policy.get("passing_only") or policy["manual_ids"] or policy.get("limit") or policy.get("quantile"))
    identifiers = [str(row.get("candidate_id") or "") for row in candidates]
    if len(set(identifiers)) != len(identifiers) or not all(identifiers):
        raise ValueError("Selection requires unique nonempty candidate IDs.")
    if set(policy["manual_ids"]) - set(identifiers):
        raise ValueError("Manual selection contains unknown candidate IDs.")
    outputs, groups = [], {}
    for candidate in candidates:
        row = copy.deepcopy(candidate)
        observations = {field: _observation(candidate, field) for field in required}
        reasons = []
        notes = []
        hard_measurement_fields = set(policy["required_metrics"])
        if policy["ranking_metric"]:
            hard_measurement_fields.add(policy["ranking_metric"])
        if policy.get("diversity_field"):
            hard_measurement_fields.add(policy["diversity_field"])
        hard_measurements = [observations[field] for field in hard_measurement_fields]
        status = ("failed" if any(value[1] == "failed" for value in hard_measurements)
                  else "missing" if any(value[1] != "valid" for value in hard_measurements)
                  else "pass" if explicit else "not_evaluated")
        for field in hard_measurement_fields:
            value = observations[field]
            if value[1] != "valid":
                reasons.append(f"{field}: {value[2]}")
        if policy.get("passing_only"):
            from mn_molecular_design_studio.design_validation import candidate_passes_filters
            if not candidate_passes_filters(candidate):
                reasons.append("No affirmative native filter pass")
                if status == "pass":
                    status = "reject"
        rule_results = []
        for rule in policy["rules"]:
            value, observation_status, error, _ = observations[rule["field"]]
            if observation_status != "valid":
                rule_results.append(None)
                notes.append(f"{rule['field']}: {error}")
                continue
            if rule["operator"] == "between":
                passed = float(rule["value"]) <= value <= float(rule.get("upper_value", rule.get("max_value")))
                comparison = f"between {float(rule['value']):g} and {float(rule.get('upper_value', rule.get('max_value'))):g}"
            else:
                passed = OPERATORS[rule["operator"]](value, float(rule["value"]))
                comparison = f"{rule['operator']} {float(rule['value']):g}"
            rule_results.append(passed)
            if not passed:
                notes.append(f"{rule['field']} = {value:g} fails {comparison}")
        if policy["rules"] and status == "pass":
            valid_results = [result for result in rule_results if result is not None]
            missing_rule = len(valid_results) != len(rule_results)
            if policy["rule_logic"] == "all":
                if any(result is False for result in valid_results):
                    status = "reject"
                    reasons.extend(note for note in notes if "fails" in note)
                elif missing_rule:
                    status = "failed" if any(observations[rule["field"]][1] == "failed" for rule in policy["rules"]) else "missing"
                    reasons.extend(note for note in notes if "fails" not in note)
            else:
                if any(result is True for result in valid_results):
                    if missing_rule:
                        notes.append("One or more criteria were unavailable; another criterion passed under Any logic.")
                elif missing_rule:
                    status = "failed" if any(observations[rule["field"]][1] == "failed" for rule in policy["rules"]) else "missing"
                    reasons.extend(note for note in notes if "fails" not in note)
                else:
                    status = "reject"
                    reasons.extend(note for note in notes if "fails" in note)
        rank_field = policy["ranking_metric"]
        cohort = _digest(_protocol(observations[rank_field][3])) if rank_field and observations[rank_field][3] else "unranked"
        decision = {"status": status, "selected": False, "reasons": reasons, "rank": None, "cohort": cohort,
                    "ranking_metric": rank_field, "direction": policy["direction"], "policy": policy,
                    "logic": policy["rule_logic"], "criteria_notes": notes,
                    "measurements": {field: {"value": value[0], "status": value[1], "error": value[2],
                        "candidate_id": value[3].get("candidate_id") if value[3] else None} for field, value in observations.items()}}
        row.setdefault("raw_metadata", {})["selection"] = decision
        outputs.append(row)
        if status == "pass":
            parent_id = (candidate.get("parents") or (candidate.get("raw_metadata") or {}).get("source_candidate", {}).get("parents") or [candidate.get("candidate_id")])[0]
            group_key = (cohort, str(parent_id)) if policy.get("quantile_parent_based") else (cohort, "")
            groups.setdefault(group_key, []).append(row)
    for group in groups.values():
        field = policy["ranking_metric"]
        group.sort(key=lambda row: ((-1 if policy["direction"] == "higher" else 1) * row["raw_metadata"]["selection"]["measurements"][field]["value"] if field else 0,
                                    str(row["candidate_id"])))
        clusters, kept = {}, 0
        limit = policy["limit"]
        if policy["quantile"]:
            import math
            quantile_limit = max(1, math.ceil(len(group) * policy["quantile"]))
            limit = min(limit, quantile_limit) if limit else quantile_limit
        for rank, row in enumerate(group, 1):
            decision = row["raw_metadata"]["selection"]
            decision["rank"] = rank if field else None
            if policy["manual_ids"] and row["candidate_id"] not in policy["manual_ids"]:
                decision["reasons"].append("Not manually selected")
                continue
            cluster_key = None
            if policy.get("diversity_field"):
                value, _, _, item = _observation(row, policy["diversity_field"])
                raw = item.get("raw_metadata") or {}
                cluster_key = (_digest(_protocol(item)), raw.get("clustering_run_id", ""), value)
                if clusters.get(cluster_key, 0) >= policy["max_per_cluster"]:
                    decision["reasons"].append("Diversity cluster quota reached")
                    continue
            if limit and kept >= limit:
                decision["reasons"].append("Shortlist or quantile limit reached within comparable protocol")
                continue
            decision["selected"] = True
            decision["selection_method"] = "manual" if policy["manual_ids"] else "configured_rules"
            kept += 1
            if cluster_key is not None:
                clusters[cluster_key] = clusters.get(cluster_key, 0) + 1
    return outputs


def selected_candidates(outcomes):
    return sorted([row for row in outcomes if row["raw_metadata"]["selection"]["selected"]],
                  key=lambda row: (row["raw_metadata"]["selection"]["cohort"], row["raw_metadata"]["selection"]["rank"] or 0, row["candidate_id"]))


def write_selection_artifacts(destination, outcomes):
    """Retain every decision plus a self-contained selected sequence/structure bundle."""
    import pandas as pd
    from mn_molecular_design_studio.output_vocabulary import evidence_rows
    destination = Path(destination); destination.mkdir(parents=True, exist_ok=True)
    (destination / "all_candidates.jsonl").write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in outcomes))
    pd.DataFrame([{"candidate_id": row["candidate_id"], **{key: row["raw_metadata"]["selection"][key] for key in
                  ("status", "selected", "reasons", "rank", "cohort", "ranking_metric", "direction")}} for row in outcomes]).to_csv(destination / "selection_decisions.csv", index=False)
    pd.DataFrame(evidence_rows(outcomes)).to_csv(destination / "final_evidence.csv", index=False)
    selected, fasta = [], []
    for index, candidate in enumerate(selected_candidates(outcomes), 1):
        row = copy.deepcopy(candidate)
        files = []
        for field in ("complex_pdb", "binder_pdb", "target_pdb"):
            value = candidate.get(field)
            if not value:
                continue
            source = Path(value)
            if not source.is_file():
                raise ValueError(f"Selected structure is unavailable: {value}")
            target = destination / "structures" / f"{index:06d}_{field}{source.suffix}"
            target.parent.mkdir(exist_ok=True)
            shutil.copyfile(source, target)
            row[field] = str(target.relative_to(destination))
            files.append({"field": field, "path": row[field], "source": str(source), "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
        row["raw_metadata"]["selection_export"] = {"structures": files, "sequence": row.get("binder_sequence"), "sequence_chain_separator": ":"}
        sequence = row.get("binder_sequence")
        if sequence:
            for chain, part in enumerate(str(sequence).split(":"), 1):
                fasta.append(f">selected_{index:06d}_chain_{chain}\n{part}\n")
        selected.append(row)
    (destination / "selected_sequences.fasta").write_text("".join(fasta))
    (destination / "selected_candidates.jsonl").write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in selected))
    (destination / "selection_manifest.json").write_text(json.dumps({"schema": "studio.selection.v1", "candidate_count": len(outcomes),
        "selected_count": len(selected), "selection_purpose": "downstream experimental evaluation; not experimental validation",
        "policy": outcomes[0]["raw_metadata"]["selection"]["policy"] if outcomes else {},
        "selected_candidates": "selected_candidates.jsonl", "all_candidates": "all_candidates.jsonl"}, indent=2))
