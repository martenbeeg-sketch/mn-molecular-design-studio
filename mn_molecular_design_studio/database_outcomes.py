"""Read-only scientific outcomes above immutable molecules and evaluations."""
from __future__ import annotations

from collections import defaultdict

from mn_molecular_design_studio.database_browser import _key
from mn_molecular_design_studio.database_browser_labels import analysis_inventory, target_label
from mn_molecular_design_studio.database_design_evaluations import _design_identity, _target_identity, analysis_results, connected_evaluations


def _object(value):
    return value if isinstance(value, dict) else {}


def _status(value):
    status = str(value or "").strip().casefold()
    if status in {"accepted", "passed", "pass", "accept"}:
        return "Passed"
    if status in {"rejected", "reject", "failed", "fail", "filtered", "discarded", "not_accepted"}:
        return "Rejected"
    return "Needs review"


def outcome_receipt(context, *, source=None, member=None, origin=None, native_flags=None):
    source, member, origin = _object(source), _object(member), _object(origin)
    candidate = _object(context.get("candidate"))
    raw = _object(candidate.get("raw_metadata"))
    direct = {**origin, **member, **source, **context}
    validation = _object(raw.get("validation_acceptance") or direct.get("validation_acceptance"))
    campaign = _object(raw.get("distributed_campaign") or direct.get("distributed_campaign"))
    decision = _object(direct.get("acceptance_decision"))
    evidence = _object(decision.get("evidence"))
    validation = validation or _object(evidence.get("validation_acceptance"))
    campaign = campaign or _object(evidence.get("distributed_campaign"))
    flags = {**_object(evidence.get("native_filter_metrics")), **(native_flags or {})}
    if campaign.get("validation_required") or campaign.get("acceptance_level") == "validated":
        state, level = _status(campaign.get("validated_acceptance")), "Campaign acceptance"
    elif validation:
        state, level = _status(validation.get("status")), "Final acceptance"
    elif decision:
        state, level = _status(decision.get("status") or decision.get("cohort_status")), str(direct.get("acceptance_stage") or "Recorded acceptance")
    elif direct.get("acceptance_status"):
        state, level = _status(direct["acceptance_status"]), "Saved acceptance"
    else:
        known = [value for value in flags.values() if isinstance(value, (bool, int, float)) and value in (0, 1)]
        state = "Rejected" if any(value == 0 for value in known) else "Passed" if known else "Needs review"
        level = "Native acceptance" if known else "Not recorded"
    refold = _object(raw.get("pipeline_refold_input") or direct.get("pipeline_refold_input"))
    return {"state": state, "acceptance_level": level,
            "saved_selection": direct.get("saved_selection") is True,
            "role": str(direct.get("workflow_role") or raw.get("workflow_role") or ""),
            "workflow": str(direct.get("workflow_id") or direct.get("workflow") or ""),
            "generator_engine": str(direct.get("generator_engine") or ""),
            "prediction_engine": str(direct.get("prediction_engine") or ""),
            "source_stage": str(direct.get("source_stage") or member.get("stage") or ""),
            "calculation_status": str(direct.get("status") or _object(direct.get("row")).get("status") or "").casefold(),
            "calculation_recorded": direct.get("calculation_recorded") is True,
            "target_id": str(direct.get("design_target_molecule_id") or direct.get("target_molecule_id") or direct.get("target_id") or ""),
            "design_refold": bool(refold.get("source_candidate_id")),
            "representation_role": str(direct.get("design_representation_role") or _object(raw.get("design_representation")).get("role") or ""),
            "candidate_id": str(direct.get("candidate_id") or candidate.get("candidate_id") or "")}


def workflow_role(metadata, inputs, candidate, *, mode, analysis_only=False, prediction=False):
    raw = _object(candidate.get("raw_metadata"))
    explicit = metadata.get("workflow_role") or inputs.get("workflow_role") or raw.get("workflow_role")
    if explicit:
        return str(explicit)
    if analysis_only:
        return "post_analysis"
    if mode == "generated":
        return "sequence_redesign" if candidate.get("stage") == "sequence_design" else "design_generation"
    workflow = str(inputs.get("workflow") or metadata.get("workflow") or "")
    if workflow.startswith("ligand.docking.") or metadata.get("task_group") == "ligand-docking":
        return "docking"
    if workflow in {"nesso_affinity", "ligand.affinity.nesso"}:
        return "binding_affinity"
    if prediction:
        if _object(raw.get("pipeline_refold_input")).get("source_candidate_id"):
            return "design_validation"
        return "independent_prediction"
    return ""


def _prediction_kind(index, result):
    components = index["molecules"][result["molecule_id"]]["components"]
    return "Cofolding" if len(components) > 1 or (_target(index, result) and any(item["type"] == "small_molecule" for item in components)) else "Structure prediction"


def outcome_kind(index, result):
    receipt = result.get("receipt") or {}
    stage = result.get("evaluation", {}).get("stage") or receipt.get("source_stage") or ""
    role = receipt.get("role") or ""
    run = index["runs"].get(result.get("run_key"), {})
    group = index["groups"][result["group_id"]]
    from mn_molecular_design_studio.datasets import is_input_dataset
    if is_input_dataset(group) or result["group_id"] in index["molecules"][result["molecule_id"]].get("input_dataset_groups", set()):
        return "Reference"
    if role == "post_analysis" or analysis_results([result]):
        return "Analysis"
    if role == "design_validation" or receipt.get("design_refold") and role != "independent_prediction":
        return "Design validation"
    if role == "docking" or run.get("task_group") == "ligand-docking":
        return "Binding/docking"
    if role == "binding_affinity":
        return "Binding affinity"
    if role == "independent_prediction":
        return _prediction_kind(index, result)
    if role in {"design_generation", "sequence_redesign"} or receipt.get("generator_engine"):
        return "Design"
    if (stage.startswith(("generation.", "sequence_design")) or
            group.get("collection_kind") == "generated" and run.get("task_group") == "design"):
        return "Design"
    if receipt.get("saved_selection") and receipt.get("prediction_engine"):
        return _prediction_kind(index, result)
    if result.get("kind_label") == "Input / reference":
        return "Reference"
    return "Unclassified"


def _target(index, result):
    target = result.get("evaluation", {}).get("target_id") or result.get("receipt", {}).get("target_id")
    if not target:
        molecule = index["molecules"][result["molecule_id"]]
        links = molecule.get("metadata_links_by_group", {}).get(result["group_id"], set())
        if len(links) == 1:
            target = next(iter(links))
    return target or ""


def _matches_structure(original, analysis):
    first = set(original["structure_ids"])
    second = set(analysis["structure_ids"])
    if analysis.get("source_structure_id"):
        second.add(analysis["source_structure_id"])
    first_digests = {value for value in (original.get("structure_digest"), original.get("source_structure_digest")) if value}
    second_digests = {value for value in (analysis.get("source_structure_digest"), analysis.get("structure_digest")) if value}
    if first_digests and second_digests and not first_digests & second_digests:
        return False
    return bool(first & second or first_digests & second_digests)


def scientific_outcomes(index, results):
    """Project retained outcomes; analysis additions cannot create outcome rows."""
    grouped = defaultdict(list)
    analyses = defaultdict(list)
    analyses_by_identity = defaultdict(list)
    for result in results:
        identity = _design_identity(index["molecules"][result["molecule_id"]], result["group_id"]) or ("unknown", result["molecule_id"])
        target_id = _target(index, result)
        target = _target_identity(index, target_id) if target_id else None
        context_key = _key([identity, target, result["target_pdb_id"]]) if result.get("target_pdb_id") else _key([identity, target])
        kind = outcome_kind(index, result)
        if kind == "Design" and index["groups"][result["group_id"]].get("collection_kind") not in {"generated", "curated"} and not result.get("receipt", {}).get("saved_selection"):
            analyses[context_key].append(result)
            continue
        if kind in {"Analysis", "Design validation"}:
            analyses[context_key].append(result)
            analyses_by_identity[identity].append(result)
            continue
        if kind == "Reference":
            continue
        if kind == "Unclassified" and result["kind"] == "Molecule" and index["groups"][result["group_id"]].get("collection_kind") in {"evaluation", "inputs"}:
            continue
        if kind == "Design":
            key = (kind, context_key if target else _key([context_key, result["molecule_id"], result["run_key"]]))
        elif kind in {"Binding/docking", "Structure prediction", "Cofolding", "Binding affinity"}:
            evaluation = result.get("evaluation") or {}
            revision = evaluation.get("target_revision") or _key(["unknown-target-revision", result["run_key"]])
            key = (kind, _key([context_key, revision, result.get("receipt", {}).get("role"),
                               result["run_key"] if kind == "Structure prediction" else ""]))
        else:
            key = (kind, _key([context_key, result["run_key"], result["molecule_id"]]))
        grouped[(key, context_key)].append(result)
    outcomes = []
    for ((kind, identity), context_key), originals in grouped.items():
        retained = []
        for original in originals:
            receipt = original.get("receipt") or {}
            state = receipt.get("state", "Needs review")
            level = receipt.get("acceptance_level", "Not recorded")
            if kind == "Design" and level not in {"Final acceptance", "Campaign acceptance", "Saved acceptance"}:
                final = [row["receipt"]["state"] for row in analyses[context_key]
                         if row.get("receipt", {}).get("acceptance_level") in {"Final acceptance", "Campaign acceptance"}
                         and _matches_structure(original, row)
                         and (not row.get("source_molecule_id") or row["source_molecule_id"] in {item["molecule_id"] for item in originals})]
                if final:
                    state = next(iter(set(final))) if len(set(final)) == 1 else "Needs review"
                    level = "Final acceptance" if len(set(final)) == 1 else "Conflicting acceptance decisions"
            if state == "Rejected" and not receipt.get("saved_selection"):
                continue
            if state != "Passed" and receipt.get("saved_selection"):
                state = "Manually selected"
            if kind in {"Binding/docking", "Structure prediction", "Cofolding"} and state == "Needs review" and any(row["structure_ids"] for row in originals):
                state = "Recorded"
            if kind == "Binding affinity" and state == "Needs review" and receipt.get("calculation_recorded") and receipt.get("calculation_status") in {
                "completed", "complete", "success", "succeeded"}:
                state, level = "Recorded", "Completed affinity calculation; no filter pass asserted"
            if kind == "Unclassified":
                state = "Needs review"
            retained.append({**original, "eligibility": state, "acceptance_level": level})
        if not retained:
            continue
        eligibility = next((state for state in ("Passed", "Recorded", "Manually selected", "Needs review")
                            if any(row["eligibility"] == state for row in retained)), "Needs review")
        selected = [row for row in retained if row["eligibility"] == eligibility]
        selected.sort(key=lambda row: (row.get("receipt", {}).get("representation_role") != "refolded_validation",
                                      not bool(row["structure_ids"]), not bool(row["evaluation_id"]),
                                      row.get("run_key") or "", row.get("candidate") or "", row["model"], row["id"]))
        primary = selected[0]
        admitted = [row for row in retained if row["eligibility"] in {"Passed", "Recorded", "Manually selected"}] or selected
        evaluated_runs = {row["run_key"] for row in admitted if row["evaluation_id"]}
        sources = [row for row in admitted if row["evaluation_id"] or row["run_key"] not in evaluated_runs]
        target_id = _target(index, primary)
        target_name = target_label(index, target_id)
        chemical_identity = _design_identity(index["molecules"][primary["molecule_id"]], primary["group_id"])
        linked, unresolved = connected_evaluations(index, primary, analyses_by_identity.get(chemical_identity, []))
        outcomes.append({**primary, "id": _key(["scientific-outcome", kind, identity]), "primary": primary,
            "outcome_kind": kind, "kind_label": kind, "eligibility": eligibility, "originals": sources,
            "target_id": target_id, "target_name": target_name,
            "origin_tools": sorted({row.get("receipt", {}).get("generator_engine") or row["tool"] for row in sources if row["tool"]}),
            "source_runs": len({row["run_key"] for row in sources if row["run_key"]}),
            "analysis_inventory": analysis_inventory(analysis_results(linked)),
            "retained_in_library": any(index["groups"][row["group_id"]].get("collection_kind") in {"curated", "generated"} for row in sources),
            "attached": analyses[context_key], "structure_count": len({structure for row in selected for structure in row["structure_ids"]})})
    outcomes.sort(key=lambda row: (row["outcome_kind"], row["name"].casefold(), row["id"]))
    from mn_molecular_design_studio.entity_properties import project_results
    return project_results(index, outcomes)


def filter_outcomes(outcomes, *, kind="All outcomes", eligibility="Admitted results", query="", tool="All tools"):
    terms = query.casefold().split()
    states = ({"Passed", "Recorded", "Manually selected"} if eligibility == "Admitted results" else
              {"Passed", "Recorded"} if eligibility == "Passed / recorded" else {eligibility}
              if eligibility not in {"All retained", "Saved in library"} else None)
    return [row for row in outcomes if (kind == "All outcomes" or row["outcome_kind"] == kind)
            and (eligibility != "Saved in library" or row.get("retained_in_library") or row["eligibility"] == "Recorded")
            and (states is None or row["eligibility"] in states)
            and (tool == "All tools" or tool in row["origin_tools"])
            and all(term in " ".join([row["name"], row.get("display_name", ""), row["target_name"], *row.get("aliases", []), *row["origin_tools"],
                                      *(str(source.get(field) or "") for source in row["originals"]
                                        for field in ("run_id", "candidate", "model", "collection"))]).casefold() for term in terms)]
