"""Canonical, scoped SQL labels independent of scientific identities and files."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import uuid

PROTOCOL = "entity-names.v1"
PROPERTIES_PROTOCOL = "entity-properties.v1"


class NameConflictError(ValueError):
    pass


def clean_name(value):
    if not isinstance(value, str) or any(ord(character) < 32 for character in value):
        raise ValueError("Enter a plain-text name without control characters.")
    name = " ".join(value.split())
    if not name or len(name) > 200:
        raise ValueError("Names must contain between 1 and 200 characters.")
    return name


def identifier(kind, subject_id, scope):
    return hashlib.sha256(json.dumps([kind, subject_id, scope], separators=(",", ":")).encode()).hexdigest()


def component_id(component):
    if isinstance(component, str):
        component = json.loads(component)
    field = next((field for field in ("sequence", "smiles", "ccd") if component.get(field)), "")
    if not field:
        return ""
    identity = [component["type"], field, component[field]]
    return hashlib.sha256(json.dumps(["component", identity], sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def group_scope(group):
    return "project:" + str(group["project_id"]) if group.get("project_id") else "owner:" + str(group["owner_id"])


def install_schema(connection):
    connection.execute("""CREATE TABLE IF NOT EXISTS entity_names (
        id TEXT PRIMARY KEY, kind TEXT NOT NULL, subject_id TEXT NOT NULL, scope TEXT NOT NULL,
        name TEXT NOT NULL, aliases TEXT NOT NULL, revision TEXT NOT NULL,
        authority_node TEXT NOT NULL, updated_by TEXT NOT NULL, updated_at TEXT NOT NULL,
        UNIQUE(kind,subject_id,scope))""")
    connection.execute("CREATE INDEX IF NOT EXISTS entity_names_subject ON entity_names(kind,subject_id)")
    connection.execute("""CREATE TABLE IF NOT EXISTS entity_name_history (
        id TEXT NOT NULL, revision TEXT NOT NULL, record_json TEXT NOT NULL,
        PRIMARY KEY(id,revision))""")
    connection.execute("CREATE TABLE IF NOT EXISTS entity_name_details (id TEXT PRIMARY KEY, properties TEXT NOT NULL)")


def records(connection):
    if not connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='entity_names'").fetchone():
        return []
    details = {}
    if connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='entity_name_details'").fetchone():
        details = {row["id"]: json.loads(row["properties"]) for row in connection.execute("SELECT * FROM entity_name_details")}
    return [{**dict(row), "aliases": json.loads(row["aliases"]),
             **({"properties": details[row["id"]]} if row["id"] in details else {})}
            for row in connection.execute("SELECT * FROM entity_names")]


@contextmanager
def _project_connection(home):
    root = Path(home) / "database"
    root.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(root / "molecular_library.sqlite3", timeout=30)
    connection.row_factory = sqlite3.Row
    try:
        install_schema(connection)
        with connection:
            yield connection
    finally:
        connection.close()


def project_projection(organization, home):
    projects = [dict(project) for project in organization.get("projects", []) if isinstance(project, dict)]
    if not projects:
        return organization
    with _project_connection(home) as connection:
        existing = {row["subject_id"]: row for row in records(connection) if row["kind"] == "project"}
        for project in projects:
            subject = str(project.get("id") or "")
            if not subject or not project.get("name"):
                continue
            if subject not in existing:
                name = str(project["name"])
                scope = "project:" + subject
                connection.execute("INSERT OR IGNORE INTO entity_names VALUES (?,?,?,?,?,?,?,?,?,?)", (
                    identifier("project", subject, scope), "project", subject, scope, name, "[]", "", "", "",
                    str(project.get("created_at") or "")))
            elif not existing[subject]["revision"] and existing[subject]["name"] != project["name"]:
                aliases = sorted(set(existing[subject]["aliases"]) | {existing[subject]["name"]})
                connection.execute("UPDATE entity_names SET name=?,aliases=? WHERE id=? AND revision=''",
                    (str(project["name"]), json.dumps(aliases), existing[subject]["id"]))
        existing = {row["subject_id"]: row for row in records(connection) if row["kind"] == "project"}
    for project in projects:
        label = existing.get(str(project.get("id") or ""))
        if label:
            project.update(name=label["name"], name_aliases=label["aliases"], name_revision=label["revision"],
                           name_authority=label["authority_node"])
            project["entity_properties"] = label.get("properties", {})
            if "description" in project["entity_properties"]:
                project["description"] = project["entity_properties"]["description"]
    return {**organization, "projects": projects}


def _write(connection, kind, subject_id, scope, name, user, expected_revision, aliases=(), *, properties=None):
    from mn_molecular_design_studio.cluster import config
    if config.is_worker():
        raise ValueError("Rename shared entities on their coordinator or source workstation.")
    authority = config.identity()
    if not isinstance(aliases, (list, tuple, set)) or any(not isinstance(alias, str) or len(alias) > 1000 for alias in aliases):
        raise ValueError("Aliases must be a list of short text values.")
    key = identifier(kind, subject_id, scope)
    previous = next((row for row in records(connection) if row["id"] == key), None)
    if (previous["revision"] if previous else "") != expected_revision:
        raise NameConflictError("The name changed in another session. Refresh and review the current name before saving.")
    if previous and previous["authority_node"] not in {"", authority}:
        raise PermissionError("Rename this shared name on its source workstation: " + previous["authority_node"])
    name = clean_name(name)
    if properties is not None:
        from mn_molecular_design_studio.entity_properties import validate_properties
        properties = validate_properties(kind, properties)
        properties = {**(previous.get("properties", {}) if previous else {}), **properties}
    elif previous and "properties" in previous:
        properties = previous["properties"]
    if (previous and previous["name"] == name and set(aliases) <= set(previous["aliases"])
            and properties == previous.get("properties")):
        return previous
    old_names = set(aliases) | set(previous["aliases"] if previous else [])
    if previous:
        old_names.add(previous["name"])
    old_names.discard(name)
    value = {"id": key, "kind": kind, "subject_id": subject_id, "scope": scope, "name": name,
             "aliases": sorted({str(alias) for alias in old_names if alias}, key=lambda alias: (alias.casefold(), alias)),
             "revision": uuid.uuid4().hex, "authority_node": authority, "updated_by": str(user["id"]),
             "updated_at": datetime.now(timezone.utc).isoformat(timespec="microseconds")}
    if properties is not None:
        value["properties"] = properties
    _store(connection, value)
    return value


def _store(connection, value):
    connection.execute("INSERT INTO entity_names VALUES (?,?,?,?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET "
        "name=excluded.name,aliases=excluded.aliases,revision=excluded.revision,authority_node=excluded.authority_node,"
        "updated_by=excluded.updated_by,updated_at=excluded.updated_at", tuple(
            json.dumps(value[field], ensure_ascii=False) if field == "aliases" else value[field]
            for field in ("id", "kind", "subject_id", "scope", "name", "aliases", "revision", "authority_node", "updated_by", "updated_at")))
    if "properties" in value:
        connection.execute("INSERT INTO entity_name_details VALUES (?,?) ON CONFLICT(id) DO UPDATE SET properties=excluded.properties",
                           (value["id"], json.dumps(value["properties"], ensure_ascii=False, sort_keys=True)))
    connection.execute("INSERT OR IGNORE INTO entity_name_history VALUES (?,?,?)",
        (value["id"], value["revision"], json.dumps(value, ensure_ascii=False, sort_keys=True)))


def rename_project(user, project_id, name, *, expected_revision, home=None, properties=None, aliases=(), expected_source_revision=None):
    from mn_molecular_design_studio.workspace import _studio_home, load_organization, save_organization
    from mn_molecular_design_studio.cluster.replication import locked
    from mn_molecular_design_studio.project_lifecycle import BUSY
    if not user.get("id") or user.get("disabled"):
        raise PermissionError("Sign in with an active account to rename a project.")
    from mn_molecular_design_studio.cluster import config
    if config.is_worker():
        from mn_molecular_design_studio.cluster import identity_auth, transport
        payload = {"name": name, "expected_revision": expected_revision}
        if properties is not None:
            payload.update(properties=properties, aliases=list(aliases), expected_source_revision=expected_source_revision)
        project = transport.request("/api/v1/cluster/projects/" + project_id + "/name",
            payload, token=identity_auth.current_token(), method="POST")
        organization = transport.request("/api/v1/cluster/organization", token=identity_auth.current_token())
        config.atomic_json(config.home() / "organization-cache.json", organization)
        return project
    home = Path(home or _studio_home())
    name = clean_name(name)
    with locked():
        organization = load_organization(home)
        project = next((row for row in organization.get("projects", []) if row.get("id") == project_id), None)
        if not project or (user.get("role") != "admin" and project.get("created_by") != user["id"]):
            raise PermissionError("Only the project owner or an administrator can rename this project.")
        if project.get("deleted_at") or (project.get("lifecycle") or {}).get("status") in BUSY:
            raise ValueError("Finish project cleanup before renaming; deleted projects cannot be renamed.")
        if expected_source_revision is not None:
            from mn_molecular_design_studio.entity_properties import project_snapshot
            if project_snapshot(user, project)["source_revision"] != expected_source_revision:
                raise NameConflictError("Project details changed in another session. Reload current values before saving.")
        if any(row.get("id") != project_id and not row.get("deleted_at") and str(row.get("name") or "").casefold() == name.casefold()
               for row in organization.get("projects", [])):
            raise ValueError("An active project with this name already exists.")
        with _project_connection(home) as connection:
            connection.execute("BEGIN IMMEDIATE")
            _write(connection, "project", project_id, "project:" + project_id, name, user, expected_revision, aliases, properties=properties)
        organization = project_projection(organization, home)
        save_organization(organization, home)
        publish(home=home)
        return next(row for row in organization["projects"] if row["id"] == project_id)


def target_options(index, subject, user):
    options = {}
    direct_groups = {group_id for molecule_id in subject["direct_molecule_ids"] - subject["designed_record_ids"]
                     for group_id in index["molecules"][molecule_id]["group_ids"]}
    for group_id in sorted(direct_groups):
        group = index["groups"][group_id]
        if user.get("role") != "admin" and group.get("owner_id") != user.get("id"):
            continue
        if (group.get("origin") or {}).get("project_deleted"):
            continue
        scope = group_scope(group)
        options.setdefault(scope, {"scope": scope, "group_id": group_id})
    for scope, option in options.items():
        option["record"] = next((row for row in index.get("entity_names", [])
            if row["kind"] == "target" and row["subject_id"] == subject["id"] and row["scope"] == scope), None)
    return options


def rename_target(database, user, subject_id, scope, name, *, expected_revision, properties=None, aliases=()):
    from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects
    from mn_molecular_design_studio.cluster.replication import locked
    if not user.get("id") or user.get("disabled"):
        raise PermissionError("Sign in with an active account to rename a target.")
    with locked():
        index = build_browser_index(database, user, collection="all", hide_smoke_test=False, include_related=True)
        subject = next((row for row in component_subjects(index, role="targets") if row["id"] == subject_id), None)
        option = target_options(index, subject, user).get(scope) if subject else None
        if not option:
            raise PermissionError("Only the owner of an accessible target collection or an administrator can rename it.")
        if scope.startswith("project:"):
            from mn_molecular_design_studio.project_lifecycle import require_mutable_project
            require_mutable_project(scope.split(":", 1)[1])
        aliases = set(aliases)
        for molecule_id in subject["direct_molecule_ids"]:
            molecule = index["molecules"][molecule_id]
            if len(molecule["components"]) == 1 and molecule_id not in subject["designed_record_ids"]:
                aliases.update(value for group_id, value in molecule.get("member_names", {}).items()
                               if group_scope(index["groups"][group_id]) == scope)
        with database.connection() as connection:
            connection.execute("BEGIN IMMEDIATE")
            value = _write(connection, "target", subject_id, scope, name, user, expected_revision, aliases, properties=properties)
        publish(database=database)
        return value


def visible_records(connection, groups):
    scopes = {group_scope(group) for group in groups.values()}
    return [row for row in records(connection) if row["scope"] in scopes]


def resolve_record(record, labels, groups):
    if len(record.get("components", [])) != 1 or (record.get("metadata") or {}).get("output_kind"):
        return
    group_ids = {group["id"] for group in record.get("groups", [])}
    allowed = [groups[group_id] for group_id in group_ids if group_id in groups]
    if not (record.get("metadata", {}).get("target_identity_type") or any(
            group.get("import_method") in {"user_uploads", "tool_inputs", "agent_inputs"} for group in allowed)):
        return
    subject_id = component_id(record["components"][0])
    scopes = {group_scope(group) for group in allowed}
    names = [row for row in labels if row["kind"] == "target" and row["subject_id"] == subject_id and row["scope"] in scopes]
    if names:
        record["original_name"] = record["name"]
        record["name"] = " / ".join(sorted({row["name"] for row in names}))
        record["name_aliases"] = sorted({alias for row in names for alias in row["aliases"]})


def target_names(index, subject_id, group_ids):
    if not index.get("entity_names"):
        return []
    scopes = {group_scope(index["groups"][group_id]) for group_id in group_ids}
    if "_target_name_index" not in index:
        lookup = {}
        for row in index["entity_names"]:
            if row["kind"] == "target":
                lookup.setdefault(row["subject_id"], []).append(row)
        index["_target_name_index"] = lookup
    return sorted((row for row in index["_target_name_index"].get(subject_id, []) if row["scope"] in scopes), key=lambda row: row["id"])


def publish(*, database=None, home=None):
    from mn_molecular_design_studio.cluster import config
    if not config.is_peer():
        return
    from mn_molecular_design_studio.cluster.replication import Replica
    if database is None:
        from mn_molecular_design_studio.molecular_database import MolecularDatabase
        database = MolecularDatabase(Path(home) / "database" if home else None)
    with database.connection() as connection:
        values = {"entity_name:" + row["id"]: row for row in records(connection)
                  if row["authority_node"] == config.identity() and row["revision"]}
    replica = Replica()
    with replica.db() as connection:
        for key, value in values.items():
            previous = connection.execute("SELECT value FROM records WHERE key=?", (key,)).fetchone()
            if previous and json.loads(previous["value"]) == value:
                continue
            if connection.execute("SELECT 1 FROM conflicts WHERE key=?", (key,)).fetchone():
                continue
            replica._write(connection, key, "entity_name", value)


def apply_record(connection, key, value, author):
    required = {"id", "kind", "subject_id", "scope", "name", "aliases", "revision", "authority_node", "updated_by", "updated_at"}
    if not isinstance(value, dict) or set(value) not in (required, required | {"properties"}):
        raise ValueError("Invalid canonical name record; name deletion is not supported.")
    if any(not isinstance(value[field], str) for field in value if field not in {"aliases", "properties"}):
        raise ValueError("Invalid canonical name fields.")
    if (value["kind"] not in {"target", "project", "dataset", "entry", "result"} or value["authority_node"] != author
            or key != "entity_name:" + identifier(value["kind"], value["subject_id"], value["scope"])
            or key != "entity_name:" + value["id"] or not value["revision"]
            or not isinstance(value["aliases"], list) or any(not isinstance(alias, str) for alias in value["aliases"])
            or not str(value["scope"]).startswith(("project:", "owner:"))):
        raise ValueError("Invalid canonical name identity or authority.")
    if value["kind"] == "project" and value["scope"] != "project:" + value["subject_id"]:
        raise ValueError("Project name scope must match its stable project ID.")
    if clean_name(value["name"]) != value["name"] or not value["scope"].split(":", 1)[1]:
        raise ValueError("Invalid canonical name or empty visibility scope.")
    if "properties" in value:
        from mn_molecular_design_studio.entity_properties import validate_properties
        validate_properties(value["kind"], value["properties"])
    previous = next((row for row in records(connection) if row["id"] == value["id"]), None)
    if previous and previous["authority_node"] not in {"", author}:
        raise PermissionError("Conflicting name authority: review on the source workstation; no name was overwritten.")
    if previous:
        if "properties" not in value and "properties" in previous:
            value = {**value, "properties": previous["properties"]}
        value = {**value, "aliases": sorted(set(value["aliases"]) | set(previous["aliases"]) |
            ({previous["name"]} if previous["name"] != value["name"] else set()))}
    value = {**value, "aliases": [alias for alias in value["aliases"] if alias != value["name"]]}
    _store(connection, value)
