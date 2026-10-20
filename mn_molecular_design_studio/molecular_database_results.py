"""Run-scoped result access and exact-subset promotion into the molecular library."""
from __future__ import annotations

import json
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


def _list_authorized_run(run_dir, user, database):
    selected, task_group, run_id, owner_id, project_id = _authorized_run(run_dir, user)
    # Authorization is complete before trusted ingestion sees the filesystem path.
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    ingest_run(selected, database=database)
    sources, source_groups = _run_sources(database, task_group, run_id)
    molecule_ids = list(dict.fromkeys(row["molecule_id"] for row in sources))
    if not molecule_ids:
        return selected, task_group, run_id, owner_id, project_id, [], {}
    result = database.list_molecules(user, molecule_ids=molecule_ids, collection="all", limit=min(100_000, max(1, len(molecule_ids))))
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
        run_group_ids = sorted(by_molecule.get(record["id"], set()))
        run_group_ids.sort(key=lambda group_id: (
            {"generated": 0, "curated": 1, "inputs": 2, "evaluation": 3, "intermediate": 4}.get(
                accessible_groups[group_id]["collection_kind"], 5), group_id))
        record["run_group_ids"] = run_group_ids
        record["run_groups"] = [accessible_groups[group_id] for group_id in run_group_ids]
        record["source_group_id"] = run_group_ids[0] if run_group_ids else None
        items.append(record)
    return selected, task_group, run_id, owner_id, project_id, items, accessible_groups


def list_run_molecules(run_dir, user, *, database=None) -> list[dict]:
    """Return normal molecular records associated with this accessible Studio run."""
    database = database or MolecularDatabase()
    return _list_authorized_run(run_dir, user, database)[5]


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

    database.add_to_group(user, group_id, requested)
    source_keys = []
    source_group_ids = set()
    sources, _ = _run_sources(database, task_group, run_id)
    for source in sources:
        if source["molecule_id"] in requested:
            source_keys.append(source["source_key"])
            if source["group_id"] in source_groups:
                source_group_ids.add(source["group_id"])
    # The helper's accessible group map is intentionally limited to this run.
    _copy_run_associations(database, task_group, run_id, source_keys, sorted(source_group_ids), group_id, requested)
    return {"group_id": group_id, "molecule_ids": requested, "saved_count": len(requested)}
