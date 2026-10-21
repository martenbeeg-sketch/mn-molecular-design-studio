"""Reproducible cohort selection and ranking, independent of presentation."""
from __future__ import annotations

import math
import json
import pandas as pd


def prepare_rows(rows, *, combine_same_subject=False, evaluator_runs=None):
    """Explicitly join distinct evaluators only for the identical molecular subject.

    Conflicting repetitions are omitted, rather than resolved by latest-value
    guessing. Users can restrict the contributing runs to resolve a conflict.
    """
    from mn_molecular_design_studio.evaluations import identifier
    if type(combine_same_subject) is not bool or evaluator_runs is not None and (not isinstance(evaluator_runs, list) or not all(isinstance(run, str) for run in evaluator_runs)):
        raise ValueError("Invalid evaluation composition settings.")
    selected = []
    for row in rows:
        evaluation = row.get("evaluation") or row.get("selected_evaluation") or {}
        if evaluator_runs is not None and evaluation.get("evaluator", "") + " / " + evaluation.get("run_id", "") not in evaluator_runs:
            continue
        selected.append(row)
    if not combine_same_subject:
        return selected
    groups = {}
    for row in selected:
        evaluation = row.get("evaluation") or row.get("selected_evaluation") or {}
        subject = evaluation.get("subject_id") if evaluation.get("subject_id") and evaluation.get("representation") and evaluation.get("status", "completed") == "completed" else evaluation.get("id")
        groups.setdefault((row.get("molecule_id") or row.get("id"), subject), []).append(row)
    result = []
    for group in groups.values():
        if len(group) == 1:
            result.append(group[0])
            continue
        evaluations = [r.get("evaluation") or r.get("selected_evaluation") or {} for r in group]
        ids = sorted(e["id"] for e in evaluations)
        scores, sources, conflicts = {}, {}, []
        fields = {field for row in group for field in row.get("scores", {})}
        for field in fields:
            values = [(row["scores"][field], row.get("score_sources", {}).get(field, {}), evaluation)
                      for row, evaluation in zip(group, evaluations) if field in row.get("scores", {})]
            signatures = {(json.dumps(value, sort_keys=True, default=str), source.get("method"), source.get("unit"), source.get("scope"),
                           source.get("parameters_digest") or evaluation.get("parameters_digest"),
                           source.get("image_digest") or evaluation.get("image_digest")) for value, source, evaluation in values}
            if len(signatures) != 1:
                conflicts.append(field)
                continue
            value, source, evaluation = values[0]
            scores[field] = value
            sources[field] = {**source, "parameters_digest": evaluation.get("parameters_digest"),
                              "image_digest": evaluation.get("image_digest"), "evaluation_ids": [e["id"] for _, _, e in values]}
        evaluation = {**evaluations[0], "id": identifier("comparison-subject-", ids), "source_evaluation_ids": ids,
                      "evaluator": " + ".join(sorted({e.get("evaluator", "") for e in evaluations})),
                      "run_id": " + ".join(sorted({e.get("run_id", "") for e in evaluations})), "parameters_digest": ""}
        result.append({**group[0], "evaluation": evaluation, "scores": scores, "score_sources": sources, "conflicting_metrics": conflicts})
    return result


def compare(rows, *, fields=None, directions=None, weights=None, thresholds=None, method="mean", combine_same_subject=False, evaluator_runs=None):
    rows = prepare_rows(rows, combine_same_subject=combine_same_subject, evaluator_runs=evaluator_runs)
    if fields is not None and (not isinstance(fields, list) or not all(isinstance(field, str) for field in fields)):
        raise ValueError("Comparison fields must be a list of metric IDs.")
    if not isinstance(directions or {}, dict) or not isinstance(weights or {}, dict):
        raise ValueError("Directions and weights must be mappings by metric ID.")
    if thresholds is not None and (not isinstance(thresholds, list) or not all(isinstance(item, dict) and "field" in item and "value" in item for item in thresholds)):
        raise ValueError("Thresholds must be a list of field/operator/value rules.")
    fields, directions, weights = list(dict.fromkeys(fields or [])), directions or {}, weights or {}
    if method not in {"mean", "worst", "product", "pareto"}:
        raise ValueError("Choose mean, worst, product, or pareto ranking.")
    if len(fields) > 32 or len(rows) > 100000 or (method == "pareto" and len(rows) > 5000):
        raise ValueError("Comparison exceeds its size limit (Pareto supports 5,000 evaluations).")
    if any(directions.get(field) not in {"higher", "lower"} for field in fields):
        raise ValueError("Every ranking field needs an explicit direction.")
    if any(not math.isfinite(float(weights.get(field, 1))) or float(weights.get(field, 1)) <= 0 for field in fields):
        raise ValueError("Ranking weights must be positive finite numbers.")
    tests = []
    for item in thresholds or []:
        operator = item.get("operator", ">=")
        if operator not in {">=", "<=", ">", "<", "="} or not math.isfinite(float(item["value"])):
            raise ValueError("Invalid comparison threshold.")
        tests.append((str(item["field"]), operator, float(item["value"])))
    records = []
    for source in rows:
        scores = source.get("scores") or {}
        values = {key: float(value) for key, value in scores.items() if type(value) in {int, float} and math.isfinite(value)}
        failures = []
        for key, operator, threshold in tests:
            value = values.get(key)
            if value is None:
                failures.append(key + " missing")
            elif not {">=": value >= threshold, "<=": value <= threshold, ">": value > threshold, "<": value < threshold, "=": value == threshold}[operator]:
                failures.append(f"{key} {operator} {threshold:g}")
        missing = [key for key in fields if key not in values]
        evaluation = source.get("evaluation") or source.get("selected_evaluation") or {}
        if evaluation.get("status", "completed") != "completed":
            failures.append("evaluation " + str(evaluation.get("status")))
        metadata = source.get("metadata") or {}
        row = {"Molecule ID": source.get("molecule_id") or source.get("id"), "Name": source.get("name"),
               "Compound ID": metadata.get("compound_id", ""), "Evaluation ID": evaluation.get("id", ""),
               "Evaluator": evaluation.get("evaluator", ""), "Target": evaluation.get("target_id", ""),
               "Run": evaluation.get("run_id", ""), "Native rank": next((v for k, v in values.items() if k.endswith(".rank")), None),
               "Source evaluations": evaluation.get("source_evaluation_ids") or [evaluation.get("id", "")],
               "Conflicting metrics": ", ".join(source.get("conflicting_metrics", [])),
               "Pass": not failures and evaluation.get("status", "completed") == "completed", "Filter failures": "; ".join(failures),
               "Metric coverage": f"{len(fields)-len(missing)}/{len(fields)}", "Missing metrics": ", ".join(missing),
               "Comparison rank": None, **values}
        row["_eligible"] = not missing and not failures and bool(fields) and evaluation.get("status", "completed") == "completed"
        records.append(row)
    frame = pd.DataFrame(records)
    if frame.empty:
        return frame
    eligible = frame.index[frame["_eligible"]]
    if len(eligible):
        ranks = pd.DataFrame({field: frame.loc[eligible, field].rank(ascending=directions[field] == "lower", method="average")
                              * float(weights.get(field, 1)) for field in fields})
        if method == "mean":
            score = ranks.sum(axis=1) / sum(float(weights.get(field, 1)) for field in fields)
        elif method == "worst":
            score = ranks.max(axis=1)
        elif method == "product":
            score = sum(frame.loc[eligible, field].rank(ascending=directions[field] == "lower", method="average").map(math.log)
                        * float(weights.get(field, 1)) for field in fields)
        else:
            # Pareto fronts use the actual measures; plotting embeddings are never scores.
            remaining, fronts, front = list(eligible), {}, 0
            oriented = {i: [frame.loc[i, f] * (-1 if directions[f] == "higher" else 1) for f in fields] for i in eligible}
            while remaining:
                current = [i for i in remaining if not any(j != i and all(a <= b for a, b in zip(oriented[j], oriented[i]))
                           and any(a < b for a, b in zip(oriented[j], oriented[i])) for j in remaining)]
                front += 1
                fronts.update({i: front for i in current})
                remaining = [i for i in remaining if i not in fronts]
            score = pd.Series(fronts)
        frame.loc[eligible, "Comparison rank"] = score.rank(method="min") if method != "pareto" else score
    return frame.drop(columns="_eligible").sort_values(["Comparison rank", "Name", "Evaluation ID"], na_position="last", kind="stable").reset_index(drop=True)


def compatibility_warnings(rows, fields):
    warnings = []
    targets = {(row.get("evaluation") or row.get("selected_evaluation") or {}).get("target_id") for row in rows}
    target_fields = [f for f in fields if not f.startswith(("common.sequence.", "common.geometry.binder.", "common.secondary_structure.binder.", "sequence_qc.", "proteinsol.", "common.design.sequence_length"))]
    if target_fields and len(targets - {None, ""}) > 1:
        warnings.append("This comparison includes different targets. Select one target for binding-score ranking.")
    for field in fields:
        definitions = {(row.get("score_sources", {}).get(field, {}).get("method"),
                        row.get("score_sources", {}).get(field, {}).get("unit"),
                        row.get("score_sources", {}).get(field, {}).get("scope")) for row in rows if field in row.get("scores", {})}
        if len(definitions) > 1:
            warnings.append(field + " includes different calculation methods or scopes. Select a single evaluator for ranking.")
        protocols = {str(row.get("score_sources", {}).get(field, {}).get("parameters_digest") or (row.get("evaluation") or row.get("selected_evaluation") or {}).get("parameters_digest") or "")
                     for row in rows if field in row.get("scores", {})}
        if len(protocols - {""}) > 1:
            warnings.append(field + " was calculated with different parameters. Select a single protocol/run for ranking.")
        images = {str(row.get("score_sources", {}).get(field, {}).get("image_digest") or (row.get("evaluation") or row.get("selected_evaluation") or {}).get("image_digest") or "")
                  for row in rows if field in row.get("scores", {})}
        if len(images - {""}) > 1:
            warnings.append(field + " was calculated with different container versions. Select a single version for ranking.")
    return list(dict.fromkeys(warnings))
