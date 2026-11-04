"""Run-scoped result access and exact-subset promotion into the molecular library."""
from __future__ import annotations

import json
import hashlib
import socket
from pathlib import Path

from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id, _json, _now


def _authorized_run(run_dir, user):
    from mn_molecular_design_studio.auth import can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record
    from mn_molecular_design_studio.runtime import runs_root
    from mn_molecular_design_studio.workspace import is_studio_job

    root = Path(runs_root()).expanduser().resolve()
    selected = Path(run_dir).expanduser().resolve()
    try:
        relative = selected.relative_to(root)
    except ValueError as exc:
        raise FileNotFoundError("Studio run not found.") from exc
    if len(relative.parts) != 2 or any(part in {"", ".", ".."} for part in relative.parts):
        raise FileNotFoundError("Studio run not found.")
    task_group, run_id = relative.parts
    if not selected.is_dir() or not is_studio_job(selected):
        raise FileNotFoundError("Studio run not found.")
    job = get_job_access_record(task_group, run_id)
    if not can_access_job(user, job):
        raise PermissionError("This Studio run is not accessible.")
    if job:
        catalog_path = Path(job.get("run_relative_path") or "")
        if catalog_path.parts != (task_group, run_id):
            raise FileNotFoundError("Studio run not found.")
    marker_path = selected / "molecular_design_studio.json"
    try:
        marker = json.loads(marker_path.read_text())
    except (OSError, ValueError):
        marker = {}
    owner_id = str(marker.get("studio_owner_id") or "")
    project_id = str(marker.get("project_id") or "")
    if job:
        if owner_id and job.get("owner_id") and owner_id != str(job["owner_id"]):
            raise PermissionError("This Studio run is not accessible.")
        if project_id and job.get("project_id") and project_id != str(job["project_id"]):
            raise PermissionError("This Studio run is not accessible.")
        owner_id = str(job.get("owner_id") or owner_id)
        project_id = str(job.get("project_id") or project_id or "")
    if not owner_id:
        metadata_path = selected / "metadata.json"
        try:
            metadata = json.loads(metadata_path.read_text())
        except (OSError, ValueError):
            metadata = {}
        owner_id = str(metadata.get("owner_id") or (metadata.get("molecular_design_studio") or {}).get("studio_owner_id") or "")
        project_id = project_id or str(metadata.get("project_id") or (metadata.get("molecular_design_studio") or {}).get("project_id") or "")
    if not owner_id:
        raise FileNotFoundError("Studio run not found.")
    return selected, task_group, run_id, owner_id, project_id


def _run_sources(database, task_group, run_id):
    prefix = f"{socket.gethostname()}/{task_group}/{run_id}:"
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT source_key,molecule_id,group_id,provenance FROM molecule_sources "
            "WHERE substr(source_key,1,?)=? AND json_extract(provenance,'$.saved_selection') IS NOT 1 ORDER BY source_key",
            (len(prefix), prefix),
        ).fetchall()
        groups = connection.execute(
            "SELECT id,name,owner_id,project_id,import_method,origin FROM molecule_groups "
            "WHERE id IN (SELECT DISTINCT group_id FROM molecule_sources WHERE substr(source_key,1,?)=? AND json_extract(provenance,'$.saved_selection') IS NOT 1)",
            (len(prefix), prefix),
        ).fetchall()
    return [dict(row) for row in rows], [dict(row) for row in groups]


def _result_candidate_token(task_group, run_id, candidate_id):
    return _id("result-", _json([socket.gethostname(), str(task_group), str(run_id), str(candidate_id)]))


def _candidate_for_explicit_selection(candidate, run_dir):
    """Resolve one authorized Results candidate into its selected library identity."""
    from mn_molecular_design_studio.molecular_database_ingest import _candidate_acceptance_state, _record

    record = _record(candidate, run_dir)
    if record is None:
        return None
    identity = dict(record)
    identity["files"] = list(record.get("files") or [])
    sequence = candidate.get("binder_sequence") or candidate.get("sequence")
    binder_chains = {str(value) for value in candidate.get("binder_chains") or []}
    if sequence:
        binders = [component for component in record.get("components") or []
                   if component.get("role") == "binder" or str(component.get("chain") or "") in binder_chains]
        if binders:
            identity["components"] = binders
    state = _candidate_acceptance_state(candidate)
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    identity["metadata"] = {
        **(identity.get("metadata") or {}),
        "result_selection": True,
        "original_acceptance_status": state,
        "acceptance_evidence": raw.get("validation_acceptance") or raw.get("distributed_campaign") or {},
    }
    return identity


def _unindexed_result_rows(run_dir, task_group, run_id, owner_id, project_id, indexed_sources, items):
    """Expose valid candidates in Results without promoting or ingesting them."""
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.molecular_database_ingest import _candidate_acceptance_state, _record
    from mn_molecular_design_studio.molecule_imports import canonical_components

    indexed_candidates = set()
    for source in indexed_sources:
        try:
            provenance = json.loads(source["provenance"])
        except (TypeError, ValueError):
            continue
        if provenance.get("candidate_id"):
            indexed_candidates.add(str(provenance["candidate_id"]))
    for candidate in read_candidates(run_dir):
        candidate_id = str(candidate.get("candidate_id") or "")
        if not candidate_id or candidate_id in indexed_candidates:
            continue
        record = _record(candidate, run_dir)
        if record is None:
            continue
        identity = _candidate_for_explicit_selection(candidate, run_dir)
        if identity is None:
            continue
        components, fingerprint, molecule_type = canonical_components(identity["components"])
        if identity.get("identity_context"):
            fingerprint = hashlib.sha256((fingerprint + _json(identity["identity_context"])).encode()).hexdigest()
        molecule_id = "mol-" + fingerprint[:32]
        files = identity.get("files") or []
        representations = []
        for index, item in enumerate(files):
            path = Path(item["path"])
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            representation_id = _id("result-structure-", _json([candidate_id, digest, index]))
            representations.append({
                "id": representation_id,
                "label": str(item.get("label") or path.name),
                "format": path.suffix.lower().lstrip("."),
                "digest": digest,
                "available_in_results": True,
            })
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        items.append({
            "id": _result_candidate_token(task_group, run_id, candidate_id),
            "molecule_id": molecule_id,
            "name": str(identity.get("name") or candidate_id),
            "type": molecule_type,
            "components": components,
            "metadata": identity.get("metadata") or {},
            "groups": [],
            "representations": representations,
            "structures": representations,
            "observations": [],
            "scores": dict(candidate.get("metrics") or {}),
            "run_group_ids": [],
            "run_groups": [],
            "source_group_id": None,
            "result_only": True,
            "candidate_id": candidate_id,
            "acceptance_status": _candidate_acceptance_state(candidate),
            "source_run": {"task_group": task_group, "run_id": run_id,
                           "owner_id": owner_id, "project_id": project_id},
            "has_structure": bool(representations),
            "native_metadata": {key: raw[key] for key in (
                "validation_acceptance", "distributed_campaign", "engine_version", "structure_kind") if key in raw},
        })


def _list_authorized_run(run_dir, user, database):
    selected, task_group, run_id, owner_id, project_id = _authorized_run(run_dir, user)
    # This is a read path. Candidate previews remain in Results until a user
    # explicitly promotes them or a separate completion/index operation runs.
    sources, source_groups = _run_sources(database, task_group, run_id)
    molecule_ids = list(dict.fromkeys(row["molecule_id"] for row in sources))
    result = (database.list_molecules(user, molecule_ids=molecule_ids, collection="all",
                                      limit=min(100_000, max(1, len(molecule_ids))))
              if molecule_ids else {"items": []})
    accessible_groups = {}
    for group in source_groups:
        if database._can_read(group, user):
            try:
                origin = json.loads(group["origin"])
            except (TypeError, ValueError):
                origin = {}
            accessible_groups[group["id"]] = {
                "id": group["id"], "name": group["name"], "project_id": group["project_id"],
                "collection_kind": database._collection_kind({"origin": origin}),
                "origin": origin,
            }
    by_molecule = {}
    for source in sources:
        if source["group_id"] in accessible_groups:
            by_molecule.setdefault(source["molecule_id"], set()).add(source["group_id"])
    items = []
    for record in result["items"]:
        from mn_molecular_design_studio.evaluations import project_observation, evaluation_descriptor
        run_observations = [observation for observation in record["observations"] if observation["run_id"] == run_id]
        if run_observations:
            run_observations.sort(key=lambda obs: (not evaluation_descriptor(obs)["preferred"], obs["id"]))
            chosen = run_observations[0]
            record["scores"], record["score_sources"], _ = project_observation(chosen)
            record["selected_evaluation"] = evaluation_descriptor(chosen)
            record["observations"] = run_observations
            record["structures"] = [s for s in record.get("representations", []) if s["context"].get("run_id") == run_id]
        run_group_ids = sorted(by_molecule.get(record["id"], set()))
        run_group_ids.sort(key=lambda group_id: (
            {"generated": 0, "curated": 1, "inputs": 2, "evaluation": 3, "intermediate": 4}.get(
                accessible_groups[group_id]["collection_kind"], 5), group_id))
        record["run_group_ids"] = run_group_ids
        record["run_groups"] = [accessible_groups[group_id] for group_id in run_group_ids]
        record["source_group_id"] = run_group_ids[0] if run_group_ids else None
        items.append(record)
    _unindexed_result_rows(selected, task_group, run_id, owner_id, project_id, sources, items)
    return selected, task_group, run_id, owner_id, project_id, items, accessible_groups


def list_run_molecules(run_dir, user, *, database=None) -> list[dict]:
    """Return normal molecular records associated with this accessible Studio run."""
    database = database or MolecularDatabase()
    return _list_authorized_run(run_dir, user, database)[5]


def list_saved_run_collections(run_dir, user, *, database=None) -> list[dict]:
    """Find accessible explicit selections linked to an authorized source run.

    Saved selections remain separate from automatic run outputs: a Results
    cleanup must preserve them, while its navigation may still show the link.
    """
    database = database or MolecularDatabase()
    _, task_group, run_id, _, _ = _authorized_run(run_dir, user)
    from mn_molecular_design_studio.cluster import config

    with database.connection() as connection:
        rows = connection.execute(
            "SELECT g.*,s.molecule_id,s.provenance FROM molecule_groups g "
            "JOIN molecule_sources s ON s.group_id=g.id "
            "WHERE json_extract(s.provenance,'$.task_group')=? "
            "AND json_extract(s.provenance,'$.run_id')=? "
            "AND json_extract(s.provenance,'$.saved_selection')=1 "
            "ORDER BY g.name,g.id,s.molecule_id", (task_group, run_id),
        ).fetchall()
    groups = {}
    for row in rows:
        if not database._can_read(row, user):
            continue
        provenance = json.loads(row["provenance"])
        source_node = provenance.get("source_node")
        if source_node and source_node != config.identity():
            continue
        if not source_node and provenance.get("workstation") not in (None, "", socket.gethostname()):
            continue
        item = groups.setdefault(row["id"], {
            "id": row["id"], "name": row["name"], "project_id": row["project_id"],
            "origin": json.loads(row["origin"]), "molecule_ids": [],
        })
        if row["molecule_id"] not in item["molecule_ids"]:
            item["molecule_ids"].append(row["molecule_id"])
    return list(groups.values())


def remove_run_records(task_group, run_id, user, *, owner_id=None, project_id=None, database=None) -> dict:
    """Remove only automatically indexed Database records originating in a run.

    Run files and explicitly promoted/saved molecule selections are preserved.
    The job identity is authorized before any Database rows are changed.
    """
    from mn_molecular_design_studio.auth import can_access_job

    task_group, run_id = str(task_group or ""), str(run_id or "")
    if not task_group or not run_id or any(value in {".", ".."} or "/" in value or "\\" in value for value in (task_group, run_id)):
        raise ValueError("Invalid job reference.")
    if not can_access_job(user, {"owner_id": owner_id, "project_id": project_id}):
        raise PermissionError("This job is not accessible.")
    database = database or MolecularDatabase()
    with database.connection() as connection:
        source_rows = connection.execute(
            "SELECT source_key,molecule_id,group_id FROM molecule_sources "
            "WHERE json_extract(provenance,'$.task_group')=? AND json_extract(provenance,'$.run_id')=? "
            "AND json_extract(provenance,'$.saved_selection') IS NOT 1",
            (task_group, run_id),
        ).fetchall()
        group_rows = {
            str(row["id"]): dict(row)
            for row in connection.execute(
                "SELECT * FROM molecule_groups WHERE id IN (SELECT DISTINCT group_id FROM molecule_sources "
                "WHERE json_extract(provenance,'$.task_group')=? AND json_extract(provenance,'$.run_id')=? "
                "AND json_extract(provenance,'$.saved_selection') IS NOT 1)",
                (task_group, run_id),
            ).fetchall()
        }
    # Never remove a collaborator's Database rows merely because the job is
    # visible through a shared project.
    for group_id in group_rows:
        database.group(group_id, user, write=True)

    pairs = {(str(row["group_id"]), str(row["molecule_id"])) for row in source_rows}
    removed_members = 0
    with database.connection() as connection:
        for row in source_rows:
            connection.execute("DELETE FROM molecule_sources WHERE source_key=?", (row["source_key"],))
        for group_id, molecule_id in pairs:
            for table in ("molecule_structures", "molecule_observations"):
                connection.execute(
                    f"DELETE FROM {table} WHERE group_id=? AND molecule_id=? "
                    "AND json_extract(context,'$.task_group')=? AND json_extract(context,'$.run_id')=?",
                    (group_id, molecule_id, task_group, run_id),
                )
            origin_row = group_rows.get(group_id)
            if origin_row:
                try:
                    origin = json.loads(origin_row["origin"])
                except (TypeError, ValueError):
                    origin = {}
                if origin.get("type") == "job" and origin.get("task_group") == task_group and origin.get("run_id") == run_id:
                    relations_remain = any(
                        connection.execute(f"SELECT 1 FROM {table} WHERE group_id=? AND molecule_id=? LIMIT 1", (group_id, molecule_id)).fetchone()
                        for table in ("molecule_sources", "molecule_structures", "molecule_observations")
                    )
                    if not relations_remain:
                        cursor = connection.execute("DELETE FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id))
                        removed_members += max(0, cursor.rowcount)
        removable_groups = []
        for group_id, group in group_rows.items():
            try:
                origin = json.loads(group["origin"])
            except (TypeError, ValueError):
                origin = {}
            if (origin.get("type") == "job" and origin.get("task_group") == task_group
                    and origin.get("run_id") == run_id):
                count = connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=?", (group_id,)).fetchone()[0]
                if not count:
                    connection.execute("DELETE FROM molecule_groups WHERE id=?", (group_id,))
                    removable_groups.append(group_id)
    return {
        "removed_sources": len(source_rows),
        "removed_members": removed_members,
        "removed_groups": len(removable_groups),
        "preserved_saved_selections": True,
    }


def _copy_run_associations(database, task_group, run_id, source_keys, source_group_ids, destination_group_id, molecule_ids):
    if not molecule_ids:
        return
    id_json = _json(molecule_ids)
    group_json = _json(source_group_ids)
    run_context = {"run_id": run_id, "task_group": task_group}
    with database.connection() as connection:
        structures = connection.execute(
            "SELECT * FROM molecule_structures WHERE molecule_id IN (SELECT value FROM json_each(?)) "
            "AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at,id", (id_json, group_json)
        ).fetchall()
        observations = connection.execute(
            "SELECT * FROM molecule_observations WHERE molecule_id IN (SELECT value FROM json_each(?)) "
            "AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at,id", (id_json, group_json)
        ).fetchall()
        sources = connection.execute(
            "SELECT * FROM molecule_sources WHERE source_key IN (SELECT value FROM json_each(?)) "
            "AND molecule_id IN (SELECT value FROM json_each(?))", (_json(source_keys), id_json)
        ).fetchall()
        for row in structures:
            context = json.loads(row["context"])
            if context.get("run_id") != run_id or context.get("task_group") != task_group:
                continue
            identifier = _id("structure-", _json([row["id"], destination_group_id, run_id, task_group]))
            connection.execute("INSERT OR IGNORE INTO molecule_structures VALUES (?,?,?,?,?,?,?,?)", (
                identifier, row["molecule_id"], destination_group_id, row["path"], row["label"], row["format"], row["context"], row["created_at"]))
        for row in observations:
            context = json.loads(row["context"])
            if context.get("run_id") != run_id or context.get("task_group") != task_group:
                continue
            identifier = _id("score-", f"saved:{destination_group_id}:{row['id']}")
            connection.execute("INSERT OR IGNORE INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)", (
                identifier, row["molecule_id"], destination_group_id, row["namespace"], row["run_id"],
                row["metrics"], row["context"], row["created_at"]))
        for row in sources:
            identifier = row["source_key"] + ":saved:" + destination_group_id
            provenance = json.loads(row["provenance"])
            provenance["saved_selection"] = True
            provenance["library_source_key"] = row["source_key"]
            connection.execute("INSERT OR IGNORE INTO molecule_sources VALUES (?,?,?,?)", (
                identifier, row["molecule_id"], destination_group_id, _json(provenance)))


def save_run_selection(run_dir, user, molecule_ids, *, database=None, group_id=None,
                       name="Saved results", project_id=None) -> dict:
    """Promote exactly the selected run records into a curated user collection."""
    database = database or MolecularDatabase()
    selected, task_group, run_id, owner_id, run_project_id, items, source_groups = _list_authorized_run(run_dir, user, database)
    requested = list(dict.fromkeys(str(value) for value in molecule_ids or []))
    if not requested:
        raise ValueError("Select at least one run molecule to save.")
    records_by_id = {row["id"]: row for row in items}
    if set(requested) - set(records_by_id):
        raise PermissionError("Some selected molecules are not associated with this run or are inaccessible.")
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    candidates_by_id = {
        str(candidate.get("candidate_id") or ""): candidate
        for candidate in read_candidates(selected) if candidate.get("candidate_id")
    }
    result_selections = {}
    for selection_id in requested:
        preview = records_by_id[selection_id]
        if not preview.get("result_only"):
            continue
        candidate = candidates_by_id.get(str(preview.get("candidate_id") or ""))
        if candidate is None:
            raise FileNotFoundError("The selected Results candidate is no longer available.")
        identity_record = _candidate_for_explicit_selection(candidate, selected)
        if identity_record is None:
            raise ValueError("The selected Results candidate has no valid molecular identity.")
        result_selections[selection_id] = (candidate, identity_record, preview)
    from mn_molecular_design_studio.molecule_schema import validate_record

    requested_project = None if project_id is None else str(project_id or "")
    destination_project = requested_project
    from mn_molecular_design_studio.auth import require_project_access
    if destination_project:
        require_project_access(user, destination_project)
    if group_id:
        group = database.group(group_id, user, write=True)
        if requested_project is not None and str(group["project_id"] or "") != requested_project:
            raise ValueError("The destination group belongs to a different project.")
        require_project_access(user, group["project_id"] or None)
        for identifier in requested:
            record = records_by_id[identifier]
            validate_record({"components": record["components"], "metadata": record.get("metadata") or {}}, group["schema"])
        origin = group["origin"]
        editable_curated = (origin.get("collection_kind") == "curated" or origin.get("type") in {"user", "legacy"}
                            or group["import_method"] in {"user_uploads", "agent_inputs"})
        if not editable_curated:
            raise ValueError("Choose a curated or imported destination group.")
        destination_project = str(group["project_id"] or "")
        if origin.get("collection_kind") != "curated":
            updated_origin = {**origin, "collection_kind": "curated"}
            with database.connection() as connection:
                connection.execute("UPDATE molecule_groups SET origin=?,updated_at=? WHERE id=?",
                                   (_json(updated_origin), _now(), group_id))
    else:
        if destination_project is None:
            destination_project = run_project_id
        require_project_access(user, destination_project or None)
        for identifier in requested:
            record = records_by_id[identifier]
            validate_record({"components": record["components"], "metadata": record.get("metadata") or {}}, {})
        group_id = database.create_group(
            user, name, project_id=destination_project, import_method="user_uploads",
            origin={"type": "user", "user_id": user["id"], "collection_kind": "curated"},
        )

    promoted_ids = []
    result_source_keys = []
    for selection_id, (candidate, identity_record, preview) in result_selections.items():
        candidate_id = str(candidate.get("candidate_id") or "")
        provenance = {
            "type": "explicit_results_selection",
            "saved_selection": True,
            "collection_kind": "curated",
            "owner_id": owner_id,
            "project_id": run_project_id,
            "task_group": task_group,
            "run_id": run_id,
            "candidate_id": candidate_id,
            "workstation": socket.gethostname(),
            "source_workstation": socket.gethostname(),
            "acceptance_status": preview.get("acceptance_status", "unknown"),
            "acceptance_evidence": (candidate.get("raw_metadata") or {}).get("validation_acceptance")
                                   or (candidate.get("raw_metadata") or {}).get("distributed_campaign") or {},
            "native_metrics_preserved": True,
        }
        source_key = _id("results-selection-", _json([group_id, owner_id, task_group, run_id, candidate_id]))
        with database.connection() as connection:
            molecule_id = database._put(connection, group_id, identity_record, source_key=source_key,
                                         origin=provenance, stamp=_now(), preserve_existing=True)
        promoted_ids.append(molecule_id)
        result_source_keys.append(source_key)
        with database.connection() as connection:
            structure_ids = [str(row[0]) for row in connection.execute(
                "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=? "
                "AND json_extract(context,'$.candidate_id')=?",
                (group_id, molecule_id, candidate_id),
            ).fetchall()]
        for structure_id in structure_ids:
            database.make_structure_portable(user, structure_id)
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        if metrics:
            database.add_observation(
                group_id, molecule_id, str(candidate.get("source_tool") or "results_candidate"), run_id, metrics,
                context={**provenance, "candidate": candidate, "settings": {},
                         "native_metrics": metrics, "is_representative": True},
                source_key=source_key + ":metrics",
            )

    existing_ids = [selection_id for selection_id in requested if not records_by_id[selection_id].get("result_only")]
    if existing_ids:
        database.add_to_group(user, group_id, existing_ids)
    source_keys = []
    source_group_ids = set()
    sources, _ = _run_sources(database, task_group, run_id)
    for source in sources:
        if source["molecule_id"] in existing_ids:
            source_keys.append(source["source_key"])
            if source["group_id"] in source_groups:
                source_group_ids.add(source["group_id"])
    # The helper's accessible group map is intentionally limited to this run.
    _copy_run_associations(database, task_group, run_id, source_keys, sorted(source_group_ids), group_id, existing_ids)
    promoted_ids.extend(existing_ids)
    promoted_ids = list(dict.fromkeys(promoted_ids))
    return {"group_id": group_id, "molecule_ids": promoted_ids, "saved_count": len(promoted_ids)}
