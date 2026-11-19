"""Reusable SQL datasets and project references, separate from result admission."""
from __future__ import annotations

import hashlib
import json
import uuid


class DatasetConflictError(ValueError):
    pass


def registration(group):
    origin = group.get("origin") or {}
    origin = json.loads(origin) if isinstance(origin, str) else origin
    value = origin.get("dataset_registration") or {}
    return value if isinstance(value, dict) and value.get("version") == 1 else {}


def is_input_dataset(group):
    return registration(group).get("role") == "inputs"


def is_dataset(group):
    from mn_molecular_design_studio.scientific_browser_tree import is_dataset as legacy_dataset
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    kind = group.get("collection_kind") or MolecularDatabase._collection_kind(group)
    return is_input_dataset(group) or kind == "reference" or legacy_dataset({**group, "collection_kind": kind})


def project_references(group):
    values = registration(group).get("project_ids") or []
    if not isinstance(values, list) or any(not isinstance(value, str) for value in values):
        return []
    return sorted(set(values) | ({group["project_id"]} if group.get("project_id") else set()))


def install_schema(connection):
    missing = not connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='dataset_project_links'").fetchone()
    connection.execute("CREATE TABLE IF NOT EXISTS dataset_project_links ("
        "group_id TEXT NOT NULL REFERENCES molecule_groups(id) ON DELETE CASCADE, "
        "project_id TEXT NOT NULL, PRIMARY KEY(group_id,project_id))")
    connection.execute("CREATE INDEX IF NOT EXISTS dataset_project_lookup ON dataset_project_links(project_id,group_id)")
    projection = ("DELETE FROM dataset_project_links WHERE group_id=NEW.id; "
        "INSERT OR IGNORE INTO dataset_project_links SELECT NEW.id,value FROM json_each(NEW.origin,'$.dataset_registration.project_ids') "
        "WHERE json_extract(NEW.origin,'$.dataset_registration.version')=1 AND type='text' AND length(value)>0;")
    for action in ("INSERT", "UPDATE"):
        connection.execute("CREATE TRIGGER IF NOT EXISTS dataset_projects_" + action.lower() + " AFTER " + action
            + " ON molecule_groups BEGIN " + projection + " END")
    connection.execute("CREATE TRIGGER IF NOT EXISTS dataset_projects_delete AFTER DELETE ON molecule_groups BEGIN "
        "DELETE FROM dataset_project_links WHERE group_id=OLD.id; END")
    if missing:
        connection.execute("INSERT OR IGNORE INTO dataset_project_links SELECT g.id,reference.value FROM molecule_groups g, "
            "json_each(g.origin,'$.dataset_registration.project_ids') reference "
            "WHERE json_extract(g.origin,'$.dataset_registration.version')=1 AND reference.type='text' AND length(reference.value)>0")


def validate_projects(user, project_ids):
    from mn_molecular_design_studio.auth import require_project_access
    from mn_molecular_design_studio.project_lifecycle import require_mutable_project
    from mn_molecular_design_studio.workspace import load_organization
    if not isinstance(project_ids, (list, tuple)) or any(not isinstance(value, str) for value in project_ids):
        raise ValueError("Project references must be a list of project IDs.")
    selected = sorted(set(project_ids))
    if len(selected) > 100:
        raise ValueError("Link at most 100 projects at a time.")
    projects = {project["id"] for project in load_organization().get("projects", []) if not project.get("deleted_at")}
    if set(selected) - projects:
        raise ValueError("A referenced project no longer exists.")
    for project_id in selected:
        require_project_access(user, project_id)
        require_mutable_project(project_id)
    return selected


def input_origin(user, project_ids=()):
    if not user.get("id") or user.get("disabled"):
        raise PermissionError("An active signed-in account is required.")
    return {"type": "user", "user_id": user["id"], "collection_kind": "curated", "dataset_registration": {
        "version": 1, "role": "inputs", "project_ids": validate_projects(user, project_ids) if project_ids else [],
        "revision": uuid.uuid4().hex}}


def input_records(records):
    return [{**record, "metadata": {**(record.get("metadata") or {}), "dataset_input": True}} for record in records]


def source_revision(group):
    return hashlib.sha256(json.dumps([group["project_id"], group["owner_id"], group.get("origin")], sort_keys=True).encode()).hexdigest()


def create_dataset(database, user, name, *, molecule_ids=(), project_ids=()):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    writable = MolecularDatabase(database.root)
    origin = input_origin(user, project_ids)
    if molecule_ids:
        writable.molecules_by_ids(user, list(molecule_ids))
    group_id = writable.create_group(user, name, origin=origin)
    if molecule_ids:
        writable.add_to_group(user, group_id, list(molecule_ids))
    return writable.group(group_id, user)


def save_project_references(database, user, group_id, project_ids, *, expected_revision):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase, _json, _now
    writable = MolecularDatabase(database.root)
    group = writable.group(group_id, user, write=True)
    if not is_dataset(group):
        raise ValueError("Project references can only be edited for reusable datasets, not automatic result collections.")
    selected = validate_projects(user, project_ids)
    with writable.connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        current = connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group_id,)).fetchone()
        if current is None or not writable._can_read(current, user):
            raise PermissionError("The dataset is no longer accessible.")
        current = {**dict(current), "origin": json.loads(current["origin"])}
        if source_revision(current) != expected_revision:
            raise DatasetConflictError("The dataset changed in another session. Reopen its project references before saving.")
        if current["origin"].get("_replica_node"):
            raise PermissionError("Edit project references on the dataset’s source workstation.")
        if connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='project_lifecycle_locks'").fetchone():
            if connection.execute("SELECT 1 FROM project_lifecycle_locks WHERE project_id IN (SELECT value FROM json_each(?)) "
                    "AND state IN ('locked','deleted')", (json.dumps(selected),)).fetchone():
                raise PermissionError("A referenced project is locked for cleanup or deleted.")
        previous = registration(current)
        updated = {**previous, "version": 1, "project_ids": selected, "revision": uuid.uuid4().hex}
        origin = {**current["origin"], "dataset_registration": updated}
        connection.execute("UPDATE molecule_groups SET origin=?,updated_at=? WHERE id=?", (_json(origin), _now(), group_id))
    return writable.group(group_id, user)


def catalog(database, user, *, project_id="", query="", kind=""):
    groups = [group for group in database.groups(user, collection="library", scope="all") if is_dataset(group)]
    if project_id:
        with database.connection() as connection:
            if connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='dataset_project_links'").fetchone():
                linked = {row[0] for row in connection.execute("SELECT group_id FROM dataset_project_links WHERE project_id=?", (project_id,))}
                groups = [group for group in groups if group["id"] in linked or group["project_id"] == project_id]
    ids = [group["id"] for group in groups]
    with database.connection() as connection:
        counts = {row["group_id"]: dict(row) for row in connection.execute(
            "SELECT mm.group_id,COUNT(*) AS members,COUNT(DISTINCT m.type) AS types,MIN(m.type) AS molecule_type,"
            "SUM(m.type='small_molecule') AS compounds FROM molecule_members mm JOIN molecules m ON m.id=mm.molecule_id "
            "WHERE mm.group_id IN (SELECT value FROM json_each(?)) GROUP BY mm.group_id", (json.dumps(ids),))}
        with_coordinates = {row[0] for row in connection.execute("SELECT DISTINCT group_id FROM molecule_structures "
            "WHERE group_id IN (SELECT value FROM json_each(?))", (json.dumps(ids),))}
    result = []
    for group in groups:
        summary = counts.get(group["id"], {})
        category = ("Empty datasets" if not summary else "Small-molecule libraries" if summary["types"] == 1 and summary["molecule_type"] == "small_molecule"
            else "Reference structures" if group["id"] in with_coordinates
            else "Protein & sequence libraries" if summary["types"] == 1 and summary["molecule_type"] in {"protein", "dna", "rna"}
            else "Mixed molecular datasets")
        searchable = " ".join([group["name"], *group.get("name_aliases", []), *group.get("entity_properties", {}).values()]).casefold()
        if (not project_id or project_id in project_references(group)) and (not kind or category == kind) and all(
                term in searchable for term in query.casefold().split()):
            result.append({**group, "category": category, "project_references": project_references(group),
                "reference_revision": source_revision(group)})
    return sorted(result, key=lambda group: (group["name"].casefold(), group["id"]))


def cleanup_references(group, project_id):
    return sorted(set(project_references(group)) - {project_id}) if is_dataset(group) else []


def display_labels(groups, projects):
    from collections import Counter
    names = {project["id"]: project["name"] for project in projects}
    duplicates = Counter(group["name"].casefold() for group in groups)
    labels = {}
    for group in groups:
        label = group["name"]
        if duplicates[label.casefold()] > 1:
            context = names.get(group.get("project_id"), "Independent library")
            label += " · " + context + " · " + str(group.get("molecule_count", 0)) + " members"
        labels[group["id"]] = label
    return labels


class DatasetView:
    input_only = True

    def __init__(self, database, group_id):
        self.database = database
        self.group_id = group_id
        self.root = database.root
        self.path = database.path

    def list_molecules(self, user, **kwargs):
        from mn_molecular_design_studio.entity_names import component_id, group_scope, resolve_record, visible_records
        from mn_molecular_design_studio.entity_properties import entry_subject, label_for
        group = self.database.group(self.group_id, user)
        limit, offset = int(kwargs.get("limit", 100)), int(kwargs.get("offset", 0))
        if not 1 <= limit <= 100_000 or offset < 0:
            raise ValueError("Invalid dataset page size or offset.")
        clause, parameters = "mm.group_id=?", [self.group_id]
        if kwargs.get("molecule_ids") is not None:
            clause += " AND m.id IN (SELECT value FROM json_each(?))"
            parameters.append(json.dumps(kwargs["molecule_ids"]))
        with self.database.connection() as connection:
            connection.execute("BEGIN")
            labels = visible_records(connection, {group["id"]: group})
            if kwargs.get("query"):
                needle = "%" + kwargs["query"] + "%"
                match = "mm.name LIKE ? OR mm.components LIKE ? OR mm.metadata LIKE ?"
                parameters.extend([needle] * 3)
                if labels:
                    connection.create_function("studio_component_id", 1, component_id, deterministic=True)
                    connection.create_function("studio_dataset_entry", 1, lambda identifier: entry_subject(group, identifier), deterministic=True)
                    match += (" OR EXISTS (SELECT 1 FROM entity_names label WHERE label.scope=? AND (label.name LIKE ? OR label.aliases LIKE ?) "
                        "AND ((label.kind='entry' AND label.subject_id=studio_dataset_entry(m.id)) OR (label.kind='target' AND EXISTS "
                        "(SELECT 1 FROM json_each(mm.components) component WHERE label.subject_id=studio_component_id(component.value)))))")
                    parameters.extend([group_scope(group), needle, needle])
                clause += " AND (" + match + ")"
            total = connection.execute("SELECT count(*) FROM molecules m JOIN molecule_members mm ON mm.molecule_id=m.id WHERE "
                + clause, parameters).fetchone()[0]
            rows = connection.execute("SELECT m.*,mm.name,mm.metadata,mm.components AS member_components "
                "FROM molecules m JOIN molecule_members mm ON mm.molecule_id=m.id WHERE " + clause
                + " ORDER BY mm.added_at,m.id LIMIT ? OFFSET ?", [*parameters, limit, offset]).fetchall()
            records = {row["id"]: {**dict(row), "components": json.loads(row["member_components"]), "metadata": json.loads(row["metadata"]),
                "groups": [{"id": group["id"], "name": group["name"]}], "origin": group["origin"], "schema": group["schema"],
                "structures": [], "sources": [], "scores": {}, "score_sources": {}, "observations": [], "selected_evaluation": None}
                for row in rows}
            selected = json.dumps(list(records))
            for structure in connection.execute("SELECT * FROM molecule_structures WHERE group_id=? "
                    "AND molecule_id IN (SELECT value FROM json_each(?)) ORDER BY created_at DESC,id", (self.group_id, selected)):
                records[structure["molecule_id"]]["structures"].append({**dict(structure), "context": json.loads(structure["context"])})
            for source in connection.execute("SELECT * FROM molecule_sources WHERE group_id=? "
                    "AND molecule_id IN (SELECT value FROM json_each(?))", (self.group_id, selected)):
                records[source["molecule_id"]]["sources"].append({"source_key": source["source_key"], **json.loads(source["provenance"])})
        for record in records.values():
            resolve_record(record, labels, {group["id"]: group})
            entry = label_for(labels, "entry", entry_subject(group, record["id"]), group_scope(group))
            if entry:
                record["original_name"] = record.get("original_name", record["name"])
                record["name"] = entry["name"]
                record["name_aliases"] = sorted(set(record.get("name_aliases", [])) | set(entry["aliases"]))
                record["metadata"].update(entry.get("properties", {}))
            record["representations"] = list(record["structures"])
        return {"items": list(records.values()), "total": total,
            "next_offset": offset + len(rows) if offset + len(rows) < total else None}

    def molecules_by_ids(self, user, ids):
        records = self.list_molecules(user, molecule_ids=ids, limit=100_000)["items"]
        if {record["id"] for record in records} != set(ids):
            raise PermissionError("Some selected molecules are not in this dataset or are inaccessible.")
        return records

    def structure_file(self, user, structure_id):
        self.database.group(self.group_id, user)
        with self.database.connection() as connection:
            if not connection.execute("SELECT 1 FROM molecule_structures WHERE id=? AND group_id=?", (structure_id, self.group_id)).fetchone():
                raise PermissionError("This structure does not belong to the selected input dataset.")
        return self.database.structure_file(user, structure_id)
