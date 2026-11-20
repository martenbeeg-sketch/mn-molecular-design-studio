"""Read-only binder Results evidence joined to the authorized scientific library."""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import re

from mn_molecular_design_studio.structure_labels import ranked_bindcraft2_labels

_BINDER_TOOLS = {"bindcraft", "bindcraft2", "kejia_idr", "peptide_trap", "rfd4_binder_acceptance_campaign"}
_BINDER_TASKS = {"targeted_protein_binder", "template_peptide_idr", "protein_binder"}
_SHA = re.compile(r"[a-f0-9]{64}")
_COMPLETED = {"completed", "complete", "success", "succeeded", "ok"}


def _record(value):
    return value if isinstance(value, dict) else {}


def _status(value):
    if value is True:
        return "passed"
    if value is False:
        return "failed"
    return str(value or "not_evaluated").strip().casefold()


def binder_candidate(candidate):
    raw = _record(candidate.get("raw_metadata"))
    return bool(raw.get("pipeline_task") in _BINDER_TASKS or candidate.get("source_tool") in _BINDER_TOOLS
                or candidate.get("binder_sequence") and candidate.get("binder_chains") and candidate.get("target_chains")
                and (str(candidate.get("stage") or "").startswith(("generation", "sequence_design"))
                     or any(raw.get(key) for key in ("candidate_set", "distributed_campaign",
                                                     "validation_acceptance", "design_representation"))))


def binder_results(job, candidates, inputs=None):
    if job.get("tool") in _BINDER_TOOLS or any(binder_candidate(row) for row in candidates):
        return True
    params = _record(_record(inputs).get("params"))
    graph = _record(params.get("graph"))
    if params.get("task") in _BINDER_TASKS or _record(graph.get("task_configuration")).get("task_id") in _BINDER_TASKS:
        return True
    return any(node.get("data", {}).get("kind") in _BINDER_TOOLS
               or (node.get("data", {}).get("kind") == "rfd4"
                   and node["data"].get("params", {}).get("task") in _BINDER_TASKS)
               for node in graph.get("nodes") or [] if isinstance(node, dict))


def binder_decisions(candidate):
    from mn_molecular_design_studio.native_filter_metrics import native_filter_evidence

    raw = _record(candidate.get("raw_metadata"))
    collection = _record(raw.get("candidate_set"))
    native = native_filter_evidence(candidate)
    if native["source"] is None and ("native_filter" in collection or "native_filter" in raw):
        recorded = collection.get("native_filter", raw.get("native_filter"))
        native = {"status": _status(recorded), "source": "Preserved native decision", "reason": ""}
    validation = _record(raw.get("validation_acceptance"))
    campaign = _record(raw.get("distributed_campaign"))
    shared = _record(collection.get("shared_filter"))
    if campaign.get("validation_required") or campaign.get("acceptance_level") == "validated":
        final = {"status": _status(campaign.get("validated_acceptance")), "basis": "Campaign validation",
                 "reasons": [campaign["validation_decision_reason"]] if campaign.get("validation_decision_reason") else []}
    elif validation.get("status"):
        final = {"status": _status(validation["status"]), "basis": "Validation and acceptance",
                 "reasons": validation.get("reasons") or []}
    elif _status(shared.get("acceptance_status")) not in {"not_evaluated", "unavailable"}:
        final = {"status": _status(shared["acceptance_status"]), "basis": "Shared acceptance policy", "reasons": []}
    elif (campaign.get("acceptance_level") == "native" and campaign.get("objective_type") != "generated"
          and campaign.get("native_acceptance") in {"accepted", "rejected"}):
        final = {"status": campaign["native_acceptance"], "basis": "Native campaign policy",
                 "reasons": [campaign["decision_reason"]] if campaign.get("decision_reason") else []}
    else:
        final = {"status": "not_evaluated", "basis": "", "reasons": []}
    return native, final


def structure_digests(job, candidates, user):
    """Read only the candidate's primary complex on this Results workstation."""
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    from mn_molecular_design_studio.evaluations import file_digest
    from mn_molecular_design_studio.job_catalog import get_job_access_record
    from mn_molecular_design_studio.runtime import runs_root
    from mn_molecular_design_studio.auth import can_access_job

    root = Path(job["run_dir"]).resolve()
    managed_root = Path(runs_root()).resolve()
    allowed_runs = {root: root.is_relative_to(managed_root)}
    digests = {}
    for row in candidates:
        candidate_id = row.get("candidate_id")
        if not candidate_id or not binder_candidate(row) or not row.get("complex_pdb"):
            continue
        path = resolve_stored_path(row["complex_pdb"], run_dir=root, must_exist=True)
        if path is None:
            continue
        path = path.resolve()
        if not path.is_relative_to(managed_root):
            continue
        relative = path.relative_to(managed_root)
        if len(relative.parts) < 3 or path.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            continue
        producing_run = managed_root.joinpath(*relative.parts[:2])
        if producing_run not in allowed_runs:
            record = get_job_access_record(*relative.parts[:2])
            allowed_runs[producing_run] = bool(record and can_access_job(user, record))
        if allowed_runs[producing_run]:
            digests[candidate_id] = file_digest(path)
    return digests


def _anchors(job, candidate):
    raw = _record(candidate.get("raw_metadata"))
    collection = _record(raw.get("candidate_set"))
    original = _record(raw.get("source_candidate"))
    origin = _record(raw.get("retained_origin"))
    anchors = {
        (str(job.get("run_id") or ""), str(candidate.get("candidate_id") or "")),
        (str(collection.get("source_run") or ""), str(collection.get("native_candidate_id") or "")),
        (str(origin.get("run_id") or ""), str(raw.get("native_candidate_id") or candidate.get("candidate_id") or "")),
    }
    # Validation wrappers may preserve a distinct underlying native candidate.
    for ref, candidate_id in (
        (raw.get("source_run_dir") or _record(raw.get("validation_reference")).get("source_run_dir"),
         original.get("candidate_id")),
        (candidate.get("source_job"), raw.get("pipeline_original_candidate_id")),
    ):
        if isinstance(ref, str) and ref and candidate_id:
            anchors.add((Path(ref).name, str(candidate_id)))
    return {(run_id, candidate_id) for run_id, candidate_id in anchors if run_id and candidate_id}


def _saved_digest(database, row, context):
    digest = context.get("portable_sha256") or context.get("sha256") or context.get("content_sha256")
    if digest:
        return digest if isinstance(digest, str) and _SHA.fullmatch(digest) else ""
    path = Path(row["path"])
    # Legacy content-addressed copies use the digest as the filename.
    if path.parent.parent == database.root / "molecular_files" and path.parent.name == path.stem[:2] and _SHA.fullmatch(path.stem):
        return path.stem
    return ""


def _saved_models(database, user, job, anchors):
    """Use scoped candidate lineage, keeping exact representations separate."""
    with database.connection() as connection:
        initialized = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='molecule_groups'"
        ).fetchone()
    if not initialized:
        return {}
    project = str(job.get("project_id") or "")
    groups = database.groups(user, project_id=project, scope="all", collection="library")
    groups = {group["id"]: group for group in groups if database._collection_kind(group) in {"generated", "curated"}
              and str(group.get("project_id") or "") == project
              and (project or group.get("owner_id") == job.get("owner_id") == user["id"])}
    if not groups or not anchors:
        return {}
    candidate_ids = json.dumps(sorted({candidate_id for _, candidate_id in anchors}))
    args = (json.dumps(list(groups)), candidate_ids)
    matches = defaultdict(list)
    with database.connection() as connection:
        # Provenance can also be attached as a run alias on a curated model.
        sources = defaultdict(set)
        for row in connection.execute(
            "SELECT group_id,molecule_id,provenance FROM molecule_sources "
            "WHERE group_id IN (SELECT value FROM json_each(?)) "
            "AND json_extract(provenance,'$.candidate_id') IN (SELECT value FROM json_each(?))", args,
        ):
            provenance = json.loads(row["provenance"])
            pair = (str(provenance.get("run_id") or ""), str(provenance.get("candidate_id") or ""))
            if pair in anchors:
                sources[(row["group_id"], row["molecule_id"], pair[1])].add(pair)
        structures = connection.execute(
            "SELECT s.* FROM molecule_structures s JOIN molecule_members m "
            "ON s.group_id=m.group_id AND s.molecule_id=m.molecule_id "
            "WHERE s.group_id IN (SELECT value FROM json_each(?)) "
            "AND json_extract(s.context,'$.candidate_id') IN (SELECT value FROM json_each(?))", args,
        ).fetchall()
        labels = ranked_bindcraft2_labels(connection, [row["id"] for row in structures if row["label"] == "Complex Pdb"])
        for row in structures:
            context = json.loads(row["context"])
            candidate_id = str(context.get("candidate_id") or "")
            origin = groups[row["group_id"]].get("origin") or {}
            run_id = str(context.get("run_id") or origin.get("run_id") or "")
            pairs = {(run_id, candidate_id)} | sources[(row["group_id"], row["molecule_id"], candidate_id)]
            for pair in pairs & anchors:
                matches[pair].append({
                    "structure_id": row["id"], "molecule_id": row["molecule_id"], "group_id": row["group_id"],
                    "label": labels.get(row["id"], row["label"]), "sha256": _saved_digest(database, row, context),
                    "source_sha256": context.get("source_structure_digest") or "",
                    "representation_role": str(context.get("design_representation_role") or ""),
                    "binder_chains": sorted(context.get("binder_chains") or []),
                    "target_chains": sorted(context.get("target_chains") or []),
                    "target_molecule_id": context.get("design_target_molecule_id") or context.get("target_molecule_id") or "",
                })
    return matches


def _assessments(database, user, models):
    from mn_molecular_design_studio.evaluations import evaluation_descriptor

    wanted = {(row["molecule_id"], row["structure_id"]): row for row in models}
    if not wanted:
        return {}
    visible = database.groups(user, scope="all", collection="all")
    projects = {group["id"]: group["project_id"] for group in visible}
    assessments = defaultdict(list)
    with database.connection() as connection:
        for row in connection.execute(
            "SELECT id,group_id,molecule_id,run_id,namespace,created_at,"
            "json_remove(context,'$.candidate.raw_metadata.binderscope') context,"
            "json_extract(context,'$.candidate.raw_metadata.binderscope.source.input_structure_digest') AS report_digest "
            "FROM molecule_observations WHERE namespace='binderscope' "
            "AND group_id IN (SELECT value FROM json_each(?)) "
            "AND molecule_id IN (SELECT value FROM json_each(?)) ORDER BY created_at DESC,id",
            (json.dumps(list(projects)), json.dumps(sorted({key[0] for key in wanted}))),
        ):
            evaluation = evaluation_descriptor(dict(row))
            model = wanted.get((row["molecule_id"], evaluation["representation_id"]))
            if not model or evaluation["status"].casefold() not in _COMPLETED:
                continue
            digest = model["sha256"]
            if (not digest or row["report_digest"] != digest or evaluation["source_structure_digest"] != digest
                    or evaluation["calculation_structure_digest"] not in {"", digest}
                    or model["group_id"] not in projects or row["group_id"] not in projects
                    or projects[model["group_id"]] != projects[row["group_id"]]
                    or (model["target_molecule_id"] and evaluation["target_id"] != model["target_molecule_id"])):
                continue
            assessments[model["structure_id"]].append({
                "observation_id": row["id"], "molecule_id": row["molecule_id"], "group_id": row["group_id"],
                "run_id": row["run_id"], "structure_id": model["structure_id"], "sha256": digest,
                "model_label": model["label"], "representation_role": model["representation_role"],
                "binder_chains": evaluation["binder_chains"], "target_chains": evaluation["target_chains"],
            })
    return assessments


def binder_outcomes(database, user, job, candidates, *, digests=None):
    """Link decisions to saved models through scoped lineage and digests."""
    if not candidates:
        return []
    if digests is None:
        digests = structure_digests(job, candidates, user)
    anchors = [_anchors(job, row) if digests.get(row.get("candidate_id")) else set() for row in candidates]
    wanted = set().union(*anchors)
    saved = _saved_models(database, user, job, wanted) if wanted and database.path.is_file() else {}
    matches = []
    for row, pairs in zip(candidates, anchors):
        digest = digests.get(row.get("candidate_id"), "")
        models = {model["structure_id"]: model for pair in pairs for model in saved.get(pair, [])}
        exact = [model for model in models.values() if digest and model["sha256"] == digest]
        # Historical refolds retain their own bytes and a digest of the run's
        # source model. They are separate viewer choices, with separate results.
        derived = [model for model in models.values() if digest and model["sha256"]
                   and model["sha256"] != digest and model["source_sha256"] == digest]
        matches.append((digest, models, sorted(exact, key=lambda model: model["structure_id"]),
                        sorted(derived, key=lambda model: model["structure_id"])))
    assessments = _assessments(database, user, [model for _, _, exact, derived in matches for model in exact + derived])
    results = []
    for row, (digest, models, exact, derived) in zip(candidates, matches):
        native, final = binder_decisions(row)
        linked = [item for model in exact for item in assessments.get(model["structure_id"], [])
                  if all(row.get(key) and sorted(row[key]) == item[key] for key in ("binder_chains", "target_chains"))]
        derived_assessments = [
            item for model in derived for item in assessments.get(model["structure_id"], [])
            if all((model.get(key) or sorted(row.get(key) or [])) == item[key]
                   for key in ("binder_chains", "target_chains"))
        ]
        preferred = min(linked + derived_assessments,
                        key=lambda item: item["representation_role"] != "refolded_validation", default=None)
        registered = bool(exact) or any(model["representation_role"] == "refolded_validation" for model in derived)
        database_status = ("registered" if registered else "structure_unverified" if not digest else
                           "different_saved_representation" if models else "not_registered")
        results.append({"candidate_id": row["candidate_id"], "native_filter": native,
                        "final_acceptance": final, "database": {"status": database_status,
                            "structures": exact, "derived_structures": derived},
                        "binderscope": linked, "binderscope_derived": derived_assessments,
                        "preferred_binderscope_observation_id": preferred["observation_id"] if preferred else ""})
    return results
