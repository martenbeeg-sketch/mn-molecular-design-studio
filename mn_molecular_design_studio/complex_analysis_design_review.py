"""Metadata-only review of exact Database models and their saved evidence.

This is a presentation projection, not a new candidate or Evaluation identity.
Coordinates and residue artifacts are retrieved only by existing explicit actions.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
import math
import re
from mn_molecular_design_studio.scientific_settings import model_identity


def short_design_label(candidate_id: str, fallback: str = "Design") -> str:
    """Short labels for display only; original identities remain untouched."""
    match = re.fullmatch(r"rfd4_0*(\d+)_.+_rep(\d+)", candidate_id, re.I)
    if match:
        return f"Design {int(match[1])} · sample {match[2]}"
    match = re.fullmatch(r"rfd4_0*(\d+)_.+", candidate_id, re.I)
    if match:
        return f"Design {int(match[1])}"
    label = candidate_id or fallback
    return label if len(label) <= 44 else label[:30] + "…" + label[-9:]


def _digest(context):
    return str(context.get("portable_sha256") or context.get("content_sha256") or
               context.get("sha256") or context.get("structure_digest") or "")


def design_review_page(database, user, group_ids, *, offset=0, limit=250, structure_ids=None):
    """Read a bounded model page with linked evidence from authorized memberships.

    An Evaluation links by exact representation, or by a recorded digest with
    compatible candidate/model identity. Shared names, runs and molecule IDs
    alone are insufficient. Multiple saved parameter sets remain separate.
    """
    from mn_molecular_design_studio.evaluations import evaluation_descriptor, PROJECTION_VERSION
    from mn_molecular_design_studio.complex_analysis_collection import acceptance_decisions, evaluation_table_rows, component_kind
    from mn_molecular_design_studio.complex_analysis_collection_ui import _database_candidate_view

    if not 1 <= limit <= 250 or offset < 0:
        raise ValueError("Invalid design review page.")
    if structure_ids is not None and (not isinstance(structure_ids, list) or len(structure_ids) > 250):
        raise ValueError("Select at most 250 exact models per metadata request.")
    ids = list(dict.fromkeys(group_ids))
    database._authorize_group_ids_for_query(user, ids)
    visible = database._visible_ids(user, collection="all")
    with database.connection() as connection:
        membership = ("s.group_id IN (SELECT value FROM json_each(?)) AND EXISTS (SELECT 1 FROM molecule_members selected "
                      "WHERE selected.molecule_id=s.molecule_id AND selected.group_id IN (SELECT value FROM json_each(?)) "
                      "AND (selected.group_id=s.group_id OR NOT EXISTS (SELECT 1 FROM molecule_structures owned "
                      "WHERE owned.group_id=selected.group_id AND owned.molecule_id=selected.molecule_id)))")
        scope = (json.dumps(visible), json.dumps(ids))
        if structure_ids is not None:
            membership += " AND s.id IN (SELECT value FROM json_each(?))"
            scope += (json.dumps(structure_ids),)
        total = connection.execute("SELECT COUNT(*) FROM molecule_structures s WHERE " + membership, scope).fetchone()[0]
        rows = connection.execute("SELECT s.*,m.name,m.components,m.metadata,g.project_id,g.origin "
            "FROM molecule_structures s JOIN molecule_members m ON m.group_id=s.group_id AND m.molecule_id=s.molecule_id "
            "JOIN molecule_groups g ON g.id=s.group_id WHERE " + membership +
            " ORDER BY m.name,s.created_at,s.id LIMIT ? OFFSET ?", (*scope, limit, offset)).fetchall()
        molecule_ids = list({row["molecule_id"] for row in rows})
        observations = connection.execute("SELECT * FROM molecule_observations WHERE molecule_id IN "
            "(SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?))",
            (json.dumps(molecule_ids), json.dumps(visible))).fetchall() if rows else []
        has_evaluations = connection.execute("SELECT 1 FROM sqlite_master WHERE name='molecule_evaluations' AND type='table'").fetchone()
        saved_evaluations = ({row["observation_id"]: json.loads(row["descriptor"]) for row in connection.execute(
            "SELECT observation_id,descriptor FROM molecule_evaluations WHERE molecule_id IN "
            "(SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?)) AND projection_version=?",
            (json.dumps(molecule_ids), json.dumps(visible), PROJECTION_VERSION))} if has_evaluations and rows else {})
        aliases = connection.execute("SELECT id,molecule_id,context FROM molecule_structures WHERE molecule_id IN "
            "(SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?))",
            (json.dumps(molecule_ids), json.dumps(visible))).fetchall() if rows else []
        sources = connection.execute("SELECT molecule_id,source_key,provenance FROM molecule_sources WHERE molecule_id IN "
            "(SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?))",
            (json.dumps(molecule_ids), json.dumps(visible))).fetchall() if rows else []
        reviews = connection.execute("SELECT subject_id,structure_ref,status,note FROM complex_review_annotations "
            "WHERE source_kind='database' AND user_id=? AND evaluation_id='' AND structure_ref IN "
            "(SELECT value FROM json_each(?))", (str(user["id"]), json.dumps([row["id"] for row in rows]))).fetchall() if rows else []
    annotation = {row["structure_ref"]: dict(row) for row in reviews}
    aliases = {row["id"]: {**dict(row), "context": json.loads(row["context"])} for row in aliases}
    by_molecule = defaultdict(list)
    source_by_molecule = defaultdict(list)
    for row in sources:
        source_by_molecule[row["molecule_id"]].append({"source_key": row["source_key"], **json.loads(row["provenance"])})
    for row in observations:
        observation = {**dict(row), "metrics": json.loads(row["metrics"]), "context": json.loads(row["context"])}
        by_molecule[row["molecule_id"]].append((observation, saved_evaluations.get(row["id"]) or evaluation_descriptor(observation)))
    samples = []
    for row in rows:
        context = json.loads(row["context"])
        digest = _digest(context)
        model = model_identity(context.get("model_id"), context.get("candidate_model_id"), context.get("source_model_id"), context.get("model_index"))
        native_id = str(context.get("candidate_id") or context.get("native_candidate_id") or "")
        associated = []
        for observation, descriptor in by_molecule[row["molecule_id"]]:
            representation = descriptor.get("representation_id")
            # Database assessments have their own execution candidate/model
            # namespace. Their explicit source binding, not those execution
            # IDs, identifies the design whose coordinates were assessed.
            source = observation["context"]
            source_structure = str(source.get("source_structure_id") or "")
            source_digest = str(source.get("source_structure_digest") or "")
            source_model = model_identity(source.get("model_id"))
            if (representation == row["id"] and source_structure in {
                    row["id"], str(context.get("_replica_structure") or "")}
                    and source_structure and digest and source_digest == digest
                    and source.get("source_molecule_id") == row["molecule_id"]
                    and str(descriptor.get("model_id") or "") in {
                        "", row["id"], model, str(context.get("_replica_structure") or "")}
                    and source_model in {"", row["id"], model, str(context.get("_replica_structure") or "")}):
                associated.append((observation, descriptor))
                continue
            if descriptor.get("model_id") and model and str(descriptor["model_id"]) != model:
                continue
            if descriptor.get("candidate_id") and native_id and descriptor["candidate_id"] != native_id:
                continue
            if representation == row["id"]:
                associated.append((observation, descriptor))
                continue
            alias = aliases.get(representation)
            if representation and (not alias or alias["molecule_id"] != row["molecule_id"]):
                continue
            if alias:
                alias_model = model_identity(alias["context"].get("model_id"), alias["context"].get("candidate_model_id"), alias["context"].get("source_model_id"), alias["context"].get("model_index"))
                alias_candidate = alias["context"].get("candidate_id") or alias["context"].get("native_candidate_id")
                if (alias_model and model and alias_model != model) or (alias_candidate and native_id and alias_candidate != native_id):
                    continue
            evidence_digest = (_digest(alias["context"]) if alias else
                               descriptor.get("source_structure_digest") or descriptor.get("structure_digest"))
            if digest and evidence_digest == digest:
                associated.append((observation, descriptor))
        components = json.loads(row["components"])
        structure = {key: row[key] for key in ("id", "molecule_id", "group_id", "label", "format", "created_at", "path")}
        structure["context"] = context
        record = {"id": row["molecule_id"], "name": row["name"], "components": components,
                  "metadata": json.loads(row["metadata"]), "structures": [structure], "representations": [structure],
                  "observations": [observation for observation, _ in associated], "sources": source_by_molecule[row["molecule_id"]]}
        candidate = _database_candidate_view(record, structure, components, [], associated)
        # Use recorded sequence fields; no chain-role or sequence inference.
        for observation, _ in associated:
            saved = observation["context"].get("candidate") or {}
            if not candidate.get("binder_sequence") and saved.get("candidate_id") == candidate["candidate_id"]:
                candidate["binder_sequence"] = saved.get("binder_sequence") or ""
        sample = {"identity": row["id"], "structure_id": row["id"], "molecule_id": row["molecule_id"],
                  "structure": structure, "record": record, "candidate": candidate,
                  "label": short_design_label(native_id, row["name"]), "model_id": model,
                  "digest": digest, "context": context, "components": components, "kind": component_kind(components),
                  "project_id": row["project_id"], "source_kind": "Database", "source_reference": "database:" + row["molecule_id"],
                  "associated": associated, "decisions": acceptance_decisions(candidate),
                  "review_annotation": annotation.get(row["id"], {"status": "Unreviewed", "note": ""})}
        sample["evaluation_entries"] = evaluation_table_rows([
            {"structure_id": row["id"], "structure_label": sample["label"], "molecule_id": row["molecule_id"],
             "observation": observation, "evaluation": descriptor, "context": observation["context"]}
            for observation, descriptor in associated])
        samples.append(sample)
    # Short names may recur across workers. Add a display suffix, never replace IDs.
    counts = Counter(sample["label"] for sample in samples)
    sources = sorted({(s["context"].get("_replica_node") or s["context"].get("workstation") or "",
                       s["context"].get("run_id") or "") for s in samples})
    used = Counter()
    for sample in samples:
        if counts[sample["label"]] > 1:
            source = (sample["context"].get("_replica_node") or sample["context"].get("workstation") or "",
                      sample["context"].get("run_id") or "")
            sample["label"] += " · source " + str(sources.index(source) + 1)
            used[sample["label"]] += 1
            if used[sample["label"]] > 1:
                sample["label"] += " · representation " + str(used[sample["label"]])
    return {"items": samples, "total": total, "next_offset": offset + limit if offset + limit < total else None}


def review_measurements(samples):
    """Method/version/parameter-specific numeric series; conflicting values stay blank."""
    series, values = {}, defaultdict(dict)
    for sample in samples:
        collected = defaultdict(list)
        for observation, descriptor in sample["associated"]:
            if str(descriptor.get("status") or "completed") not in {"completed", "succeeded", "success"}:
                continue
            for key, value in observation["metrics"].items():
                if type(value) not in {int, float} or not math.isfinite(value):
                    continue
                identity = (observation["namespace"], descriptor.get("engine_version") or "",
                            descriptor.get("parameters_digest") or "", key)
                token = hashlib.sha256(json.dumps(identity).encode()).hexdigest()[:20]
                series[token] = {"method": identity[0], "version": identity[1], "parameters_digest": identity[2], "field": key,
                                 "parameters": descriptor.get("native_parameters") or {}}
                collected[token].append(value)
        for token, numbers in collected.items():
            if len(set(numbers)) == 1:
                values[sample["structure_id"]][token] = numbers[0]
    return series, dict(values)


def selected_review_models(database, user, group_ids, structure_ids, *, system_scope=None):
    ids = list(dict.fromkeys(structure_ids))
    if len(ids) > 1000:
        raise ValueError("Select at most 1,000 models for metadata review.")
    selected = []
    for offset in range(0, len(ids), 250):
        selected.extend(design_review_page(database, user, group_ids, structure_ids=ids[offset:offset + 250])["items"])
    if system_scope is not None:
        from mn_molecular_design_studio.complex_analysis_catalog import scoped_structure_ids
        allowed = scoped_structure_ids(database, user, group_ids, system_scope, [sample["structure_id"] for sample in selected])
        selected = [sample for sample in selected if sample["structure_id"] in allowed]
    return selected


def remember_review_selection(state, token, ids, allowed):
    """One ID-based selection shared by plot/table/actions, scoped to current access."""
    key = "design-review-selected-" + token
    previous = state.get(key, [])
    selected = list(dict.fromkeys([value for value in previous if value not in allowed] +
                                 [value for value in ids if value in allowed])) if allowed else []
    if selected != state.get(key, []):
        # Previously hydrated coordinates and downstream row filters belong to
        # the old exact model selection, even within the same campaign.
        state.pop("complex-analysis-workflow-samples", None)
        state.pop("complex-analysis-workflow-selection", None)
        state["complex-analysis-assessment-open"] = False
        for old_key in list(state):
            if old_key.startswith("complex-analysis-cohort-") and not old_key.startswith("complex-analysis-cohort-table-"):
                del state[old_key]
    state[key] = selected
    return selected
