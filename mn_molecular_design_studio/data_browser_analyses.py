"""Conservative, batched post-analysis columns for saved scientific results."""
from __future__ import annotations

from collections import Counter, defaultdict
import json

from mn_molecular_design_studio.database_browser import _key, subject_results
from mn_molecular_design_studio.database_browser_labels import tool_label
from mn_molecular_design_studio.database_design_evaluations import (
    _design_identity, _target_identity, analysis_results, connected_evaluations, design_subject,
)
from mn_molecular_design_studio.binderscope_table import PREFIX, SUMMARY_PATHS, summary_columns


COMPLETED = {"completed", "complete", "success", "succeeded", "ok"}
FAILED = {"failed", "error", "unavailable", "not_available", "cancelled", "canceled"}


def _anchor(entry):
    original = entry["result"]
    evaluation = {**original.get("evaluation", {})}
    evaluation["target_id"] = evaluation.get("target_id") or entry["outcome"].get("target_id") or ""
    selected = {**original, "evaluation": evaluation}
    structure = entry.get("structure_id")
    if structure:
        selected["structure_ids"] = [structure]
        if structure != evaluation.get("representation_id") and original.get("structure_ids") != [structure]:
            for field in ("structure_digest", "source_structure_digest", "calculation_structure_digest"):
                selected[field] = ""
                evaluation[field] = ""
            selected["source_structure_id"] = ""
    return selected


def _candidates(database, index, entries, cache):
    molecules, groups = set(), set()
    for entry in entries:
        if entry["outcome"]["outcome_kind"] in {"Reference molecule", "Imported molecule"}:
            continue
        cohort = design_subject(index, entry["result"])
        if cohort:
            molecules.update(cohort["molecule_ids"])
            groups.update(cohort["group_ids"])
    if not molecules:
        return {}
    subject = {"id": _key(["table-analysis", sorted(molecules), sorted(groups)]),
               "molecule_ids": molecules, "direct_molecule_ids": molecules, "group_ids": groups}
    results = (cache.get_results(database, index, subject, include_metrics=True) if cache else
               subject_results(database, index, subject, include_metrics=True))
    candidates = defaultdict(list)
    for result in results:
        receipt = result.get("receipt") or {}
        if not result.get("evaluation_id") or receipt.get("role") == "design_validation" or receipt.get("design_refold"):
            continue
        if tool_label(result.get("tool")) == "Filter decision" or result.get("tool") == "md.production":
            continue
        if receipt.get("role") != "post_analysis" and not analysis_results([result]) and not result["tool"].startswith(("md.analysis", "md.mmgbsa")):
            continue
        identity = _design_identity(index["molecules"][result["molecule_id"]], result["group_id"])
        candidates[identity].append(result)
    return candidates


def _relation(index, entry, selected, analysis):
    if entry["outcome"]["outcome_kind"] == "Dynamics":
        anchor_context = selected.get("md") or {}
        context = analysis.get("md") or {}
        production = anchor_context.get("source_production_run_id") or selected.get("run_id")
        if not production or context.get("source_production_run_id") != production:
            return "Other MD segment"
        replica = anchor_context.get("replica_index")
        if replica is None or context.get("replica_index") is None or str(replica) != str(context["replica_index"]):
            return "Other / unverified replica"
        return "Selected MD segment"
    evaluation = analysis.get("evaluation") or {}
    anchor = selected["evaluation"]
    target, candidate_target = anchor.get("target_id"), evaluation.get("target_id")
    if bool(target) != bool(candidate_target) or target and _target_identity(index, target) != _target_identity(index, candidate_target):
        return "Target link unverified"
    for field in ("binder_chains", "target_chains", "ligand_selection"):
        first, second = anchor.get(field), evaluation.get(field)
        if first and second and first != second:
            return "Different component selection"
    return analysis["structure_relation"]


def _binderscope_summaries(database, candidates):
    requested = sorted({(analysis["observation_id"], analysis["molecule_id"], analysis["group_id"])
        for analyses in candidates.values() for analysis in analyses
        if analysis["tool"] == "binderscope" and analysis.get("observation_id")})
    if not requested:
        return {}
    paths = ",".join("'$.candidate.raw_metadata.binderscope." + path + "'" for path in SUMMARY_PATHS.values())
    with database.connection() as connection:
        return {row["id"]: {**dict(zip(SUMMARY_PATHS, json.loads(row["summary"]))), "bound_digest": row["bound_digest"],
                "method_version": row["method_version"]}
            for row in connection.execute("SELECT o.id,json_extract(o.context," + paths + ") AS summary,"
                "json_extract(o.context,'$.source_structure_digest') AS bound_digest,"
                "json_extract(o.context,'$.candidate.raw_metadata.analysis_parameters.method_version') AS method_version "
                "FROM json_each(?) requested JOIN molecule_observations o "
                "ON o.id=json_extract(requested.value,'$[0]') AND o.molecule_id=json_extract(requested.value,'$[1]') "
                "AND o.group_id=json_extract(requested.value,'$[2]')", (json.dumps(requested),))}


def linked_analysis_columns(database, index, entries, *, cache=None):
    candidates = _candidates(database, index, entries, cache)
    summaries_by_id = _binderscope_summaries(database, candidates)
    columns = []
    for entry in entries:
        if entry["outcome"]["outcome_kind"] in {"Reference molecule", "Imported molecule"}:
            columns.append({"Post-analysis": "Not applicable"})
            continue
        selected = _anchor(entry)
        identity = _design_identity(index["molecules"][selected["molecule_id"]], selected["group_id"])
        related = candidates.get(identity, [])
        if entry["outcome"]["outcome_kind"] == "Dynamics":
            linked, unresolved = [], 0
            source_run = index["runs"].get(selected.get("run_key"), {})
            for analysis in related:
                candidate_run = index["runs"].get(analysis.get("run_key"), {})
                if any(source_run.get(field) != candidate_run.get(field) for field in ("owner_id", "project_id", "node")):
                    continue
                study = (selected.get("md") or {}).get("study_id")
                if not study or study != (analysis.get("md") or {}).get("study_id"):
                    continue
                if selected.get("target_pdb_id") != analysis.get("target_pdb_id"):
                    continue
                target = selected["evaluation"].get("target_id")
                candidate_target = analysis["evaluation"].get("target_id")
                if target and candidate_target and _target_identity(index, target) != _target_identity(index, candidate_target):
                    continue
                if _relation(index, entry, selected, analysis) == "Selected MD segment":
                    linked.append({**analysis, "structure_relation": "Selected MD segment"})
        else:
            linked, unresolved = connected_evaluations(index, selected, related)
        tools = defaultdict(dict)
        for analysis in linked:
            if analysis["evaluation_id"] == selected.get("evaluation_id"):
                continue
            relation = _relation(index, entry, selected, analysis)
            tools[tool_label(analysis["tool"])][(analysis["evaluation_id"], analysis.get("observation_id"))] = {**analysis, "structure_relation": relation}
        row = {"Post-analysis": "Not recorded"}
        if entry["outcome"]["outcome_kind"] == "Design":
            row[PREFIX + "Status"] = "Not assessed"
        summaries = []
        for tool, unique in sorted(tools.items()):
            analyses = list(unique.values())
            exact = [analysis for analysis in analyses if analysis["structure_relation"] in {"Selected structure", "Selected MD segment"}]
            if tool == "BinderScope":
                for analysis in analyses:
                    summary = summaries_by_id.get(analysis.get("observation_id"), {})
                    analysis["report_available"] = bool(summary.get("digest") and summary["digest"] == summary.get("bound_digest"))
                    analysis["report_conflict"] = bool(summary.get("digest") and summary["digest"] != summary.get("bound_digest"))
                    analysis["analysis_version"] = summary.get("method_version") or analysis["evaluation"].get("engine_version")
            measured = [analysis for analysis in exact if str(analysis.get("status") or "").casefold() in COMPLETED
                and not analysis.get("report_conflict")]
            if len(measured) > 1:
                status = f"Multiple evaluations ({len(measured)})"
            elif measured:
                status = "Recorded" if measured[0]["structure_relation"] == "Selected structure" else "Recorded MD segment"
                extra = Counter("Report checksum mismatch" if analysis.get("report_conflict") else
                    "Failed / unavailable" if str(analysis.get("status") or "").casefold() in FAILED else
                    "Other model" if analysis["structure_relation"] == "Other structure of this design" else analysis["structure_relation"]
                    for analysis in analyses if analysis not in measured)
                if extra:
                    status += "; " + "; ".join(f"{label} ({count})" for label, count in sorted(extra.items()))
            elif any(analysis.get("report_conflict") for analysis in exact):
                status = "Report checksum mismatch"
            elif exact and all(str(analysis.get("status") or "").casefold() in FAILED for analysis in exact):
                status = "Failed / unavailable"
            elif exact:
                status = "Status unknown / incomplete"
            else:
                relations = {analysis["structure_relation"] for analysis in analyses}
                status = "; ".join(sorted("Other model" if relation == "Other structure of this design" else relation for relation in relations))
            prefix = "Analysis · " + tool + " · "
            row[prefix + "Status"] = status
            row[prefix + "Evaluations"] = len(analyses)
            row[prefix + "Evaluation IDs"] = "; ".join(sorted({analysis["evaluation_id"] for analysis in analyses}))
            summaries.append(tool + ": " + status)
            row[prefix + "Records"] = json.dumps([{
                "evaluation_id": analysis["evaluation_id"], "observation_id": analysis.get("observation_id"),
                "tool": analysis["tool"], "run_id": analysis.get("run_id"), "status": analysis.get("status"),
                "relationship": analysis["structure_relation"], "evaluation": analysis["evaluation"],
                "molecule_id": analysis["molecule_id"], "group_id": analysis["group_id"],
                "source_molecule_id": analysis.get("source_molecule_id"), "source_structure_digest": analysis.get("source_structure_digest"),
                "source_structure_id": analysis.get("source_structure_id"), "metrics": analysis.get("metrics") or {},
                **({"report_available": analysis.get("report_available", False), "analysis_version": analysis.get("analysis_version")}
                    if tool == "BinderScope" else {}),
            } for analysis in analyses], ensure_ascii=False, sort_keys=True)
            if len(measured) == 1:
                analysis = measured[0]
                evaluation = analysis["evaluation"]
                row[prefix + "Scored evaluation ID"] = analysis["evaluation_id"]
                row[prefix + "Engine version"] = evaluation.get("engine_version") or None
                for name, value in (analysis.get("metrics") or {}).items():
                    row[prefix + "Result · " + name] = json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list)) else value
                if tool == "BinderScope" and analysis.get("report_available"):
                    row.update(summary_columns(summaries_by_id[analysis["observation_id"]]))
                    row[PREFIX + "Analysis version"] = analysis.get("analysis_version")
        if unresolved:
            summaries.append(f"Unverified links ({unresolved})")
        if summaries:
            row["Post-analysis"] = "; ".join(summaries)
        columns.append(row)
    return columns
