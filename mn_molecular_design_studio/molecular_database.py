"""Durable molecular identities, collections, structures and score observations.

Group membership grants visibility. Scientific observations and representations
retain their own group/project context, even when molecular identities coincide.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sqlite3
import tempfile
import time
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from mn_molecular_design_studio import molecule_imports as imports
from mn_molecular_design_studio.runtime import app_home
from mn_molecular_design_studio.output_vocabulary import PREFIX, normalize_metrics, sql_common_metrics

IMPORT_METHODS = {"user_uploads": "User uploads", "tool_inputs": "Tool inputs", "tool_output": "Tool output", "agent_inputs": "Agent inputs"}
COLLECTION_KINDS = {"curated", "generated", "inputs", "evaluation", "intermediate"}
_SYNC_TIMES: dict[str, float] = {}


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False, default=str)


def _now():
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def _id(prefix, value):
    return prefix + hashlib.sha256(str(value).encode()).hexdigest()[:32]


class MolecularDatabase:
    def __init__(self, root: Path | None = None):
        self.root = Path(root) if root is not None else app_home() / "database"
        self.root.mkdir(parents=True, exist_ok=True)
        self.path = self.root / "molecular_library.sqlite3"

    @contextmanager
    def connection(self):
        connection = sqlite3.connect(self.path, timeout=30)
        connection.row_factory = sqlite3.Row
        connection.create_function("studio_common_metrics", 3, sql_common_metrics, deterministic=True)
        connection.execute("PRAGMA foreign_keys=ON")
        connection.execute("PRAGMA busy_timeout=30000")
        connection.execute("PRAGMA journal_mode=WAL")
        connection.executescript('''
            CREATE TABLE IF NOT EXISTS molecules (
                id TEXT PRIMARY KEY, fingerprint TEXT NOT NULL UNIQUE,
                type TEXT NOT NULL, components TEXT NOT NULL, created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS molecule_groups (
                id TEXT PRIMARY KEY, name TEXT NOT NULL, owner_id TEXT NOT NULL,
                project_id TEXT NOT NULL DEFAULT '', import_method TEXT NOT NULL,
                origin TEXT NOT NULL DEFAULT '{}', schema_json TEXT NOT NULL DEFAULT '{}',
                created_at TEXT NOT NULL, updated_at TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS molecular_groups_access ON molecule_groups(project_id,owner_id);
            CREATE TABLE IF NOT EXISTS molecule_members (
                group_id TEXT NOT NULL REFERENCES molecule_groups(id),
                molecule_id TEXT NOT NULL REFERENCES molecules(id),
                name TEXT NOT NULL, components TEXT NOT NULL, metadata TEXT NOT NULL DEFAULT '{}',
                added_at TEXT NOT NULL, PRIMARY KEY(group_id,molecule_id)
            );
            CREATE INDEX IF NOT EXISTS molecular_members_by_molecule ON molecule_members(molecule_id);
            CREATE TABLE IF NOT EXISTS molecule_structures (
                id TEXT PRIMARY KEY, molecule_id TEXT NOT NULL REFERENCES molecules(id),
                group_id TEXT NOT NULL REFERENCES molecule_groups(id), path TEXT NOT NULL,
                label TEXT NOT NULL, format TEXT NOT NULL, context TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS molecular_structure_lookup ON molecule_structures(molecule_id,group_id);
            CREATE TABLE IF NOT EXISTS molecule_observations (
                id TEXT PRIMARY KEY, molecule_id TEXT NOT NULL REFERENCES molecules(id),
                group_id TEXT NOT NULL REFERENCES molecule_groups(id), namespace TEXT NOT NULL,
                run_id TEXT NOT NULL, metrics TEXT NOT NULL, context TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS molecular_score_lookup ON molecule_observations(molecule_id,namespace,created_at);
            CREATE TABLE IF NOT EXISTS molecule_sources (
                source_key TEXT PRIMARY KEY, molecule_id TEXT NOT NULL REFERENCES molecules(id),
                group_id TEXT NOT NULL REFERENCES molecule_groups(id), provenance TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS molecular_ingestions (
                source_key TEXT PRIMARY KEY, signature TEXT NOT NULL, updated_at TEXT NOT NULL
            );
        ''')
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    @staticmethod
    def _can_read(group, user):
        if user.get("role") == "admin":
            return True
        project = str(group["project_id"] or "")
        if project:
            from mn_molecular_design_studio.auth import accessible_project_ids
            return project in set(accessible_project_ids(user))
        return str(group["owner_id"]) == str(user["id"])

    def group(self, group_id: str, user: dict, *, write=False) -> dict:
        with self.connection() as connection:
            row = connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group_id,)).fetchone()
        if not row or not self._can_read(row, user):
            raise PermissionError("This molecule group is not accessible.")
        if write and user.get("role") != "admin" and row["owner_id"] != user["id"]:
            raise PermissionError("Only the group owner or an administrator can change this group.")
        if write and json.loads(row["origin"]).get("_replica_node"):
            raise PermissionError("Edit this group on its source workstation, or create a local group from the selection.")
        result = dict(row)
        result["origin"] = json.loads(result["origin"])
        result["read_only"] = bool(result["origin"].get("_replica_node"))
        result["schema"] = json.loads(result.pop("schema_json"))
        return result

    @staticmethod
    def _collection_kind(group) -> str:
        origin = group.get("origin") if isinstance(group, dict) else None
        if isinstance(origin, str):
            try:
                origin = json.loads(origin)
            except (TypeError, ValueError):
                origin = {}
        origin = origin if isinstance(origin, dict) else {}
        kind = origin.get("collection_kind")
        if kind in COLLECTION_KINDS:
            return kind
        # Old automatically indexed job groups are kept in the all-indexed
        # view. Groups created or imported by users remain in the library.
        return "intermediate" if origin.get("type") == "job" else "curated"

    def groups(self, user: dict, *, project_id="", scope="organization", methods=None, collection="all") -> list[dict]:
        if collection not in {"all", "library"}:
            raise ValueError("Collection must be 'all' or 'library'.")
        with self.connection() as connection:
            rows = connection.execute("SELECT g.*, COUNT(mm.molecule_id) AS molecule_count FROM molecule_groups g "
                                      "LEFT JOIN molecule_members mm ON mm.group_id=g.id GROUP BY g.id ORDER BY g.updated_at DESC").fetchall()
        result = []
        for row in rows:
            item = {**dict(row), "origin": json.loads(row["origin"]), "schema": json.loads(row["schema_json"])}
            item["read_only"] = bool(item["origin"].get("_replica_node"))
            item["collection_kind"] = self._collection_kind(item)
            if not self._can_read(row, user) or (project_id and row["project_id"] != project_id):
                continue
            if scope == "mine" and row["owner_id"] != user["id"]:
                continue
            if not row["molecule_count"] and row["import_method"] not in {"user_uploads", "agent_inputs"}:
                continue
            if collection == "library" and item["collection_kind"] not in {"curated", "generated"}:
                continue
            if methods and row["import_method"] not in methods:
                continue
            result.append(item)
        return result

    def create_group(self, user: dict, name: str, *, project_id="", import_method="user_uploads", origin=None, group_id=None, schema=None) -> str:
        from mn_molecular_design_studio.molecule_schema import validate_schema
        validate_schema(schema or {})
        from mn_molecular_design_studio.auth import require_project_access
        require_project_access(user, project_id or None)
        if not str(name).strip():
            raise ValueError("Name the group.")
        if import_method not in IMPORT_METHODS:
            raise ValueError("Unknown import method.")
        group_id = group_id or str(uuid.uuid4())
        if len(group_id) > 160:
            raise ValueError("Invalid group ID.")
        stamp = _now()
        with self.connection() as connection:
            existing = connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group_id,)).fetchone()
            if existing:
                if not self._can_read(existing, user) or (existing["owner_id"] != user["id"] and user.get("role") != "admin"):
                    raise PermissionError("The group ID belongs to another owner.")
                return group_id
            connection.execute("INSERT INTO molecule_groups VALUES (?,?,?,?,?,?,?,?,?)", (
                group_id, str(name).strip()[:250], str(user["id"]), str(project_id or ""), import_method,
                _json(origin or {"type": "user", "user_id": user["id"]}), _json(schema or {}), stamp, stamp))
        return group_id

    def _put(self, connection, group_id: str, record: dict, *, source_key="", origin=None, stamp=None):
        components, fingerprint, kind = imports.canonical_components(record["components"])
        if record.get("identity_context"):
            fingerprint = hashlib.sha256((fingerprint + _json(record["identity_context"])).encode()).hexdigest()
        molecule_id = "mol-" + fingerprint[:32]
        stamp = stamp or _now()
        connection.execute("INSERT OR IGNORE INTO molecules VALUES (?,?,?,?,?)",
                           (molecule_id, fingerprint, kind, _json(components), stamp))
        existing = connection.execute("SELECT metadata FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone()
        metadata = {**(json.loads(existing[0]) if existing else {}), **(record.get("metadata") or {})}
        connection.execute("INSERT INTO molecule_members VALUES (?,?,?,?,?,?) ON CONFLICT(group_id,molecule_id) "
                           "DO UPDATE SET metadata=excluded.metadata", (
                               group_id, molecule_id, str(record.get("name") or molecule_id)[:250], _json(components), _json(metadata), stamp))
        if source_key:
            connection.execute("INSERT INTO molecule_sources VALUES (?,?,?,?) ON CONFLICT(source_key) DO UPDATE SET "
                               "molecule_id=excluded.molecule_id,group_id=excluded.group_id,provenance=excluded.provenance",
                               (source_key, molecule_id, group_id, _json(origin or {})))
        for item in record.get("files") or []:
            path = Path(item["path"]) if item.get("path") else self._store_file(item["name"], item["content"])
            context = {**(origin or {}), **{key: value for key, value in item.items() if key not in {"path", "content", "name", "label"}}}
            structure_id = _id("structure-", _json([molecule_id, group_id, str(path), context]))
            connection.execute("INSERT OR IGNORE INTO molecule_structures VALUES (?,?,?,?,?,?,?,?)", (
                structure_id, molecule_id, group_id, str(path), str(item.get("label") or path.name),
                path.suffix.lower().lstrip("."), _json(context), stamp))
        connection.execute("UPDATE molecule_groups SET updated_at=? WHERE id=?", (stamp, group_id))
        return molecule_id

    def _store_file(self, name: str, content: bytes) -> Path:
        suffix = Path(name).suffix.lower()
        if suffix not in {".pdb", ".ent", ".cif", ".mmcif", ".sdf", ".mol", ".fasta", ".fa"}:
            raise ValueError("Unsupported molecular structure file.")
        digest = hashlib.sha256(content).hexdigest()
        directory = self.root / "molecular_files" / digest[:2]
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / (digest + suffix)
        if not path.exists():
            with tempfile.NamedTemporaryFile(dir=directory, delete=False) as temporary:
                temporary.write(content)
                temporary_path = Path(temporary.name)
            os.replace(temporary_path, path)
        return path

    def import_records(self, user: dict, records: list[dict], *, name="Imported molecules", project_id="", group_id=None, import_method="user_uploads", origin=None) -> dict:
        from mn_molecular_design_studio.molecule_schema import validate_record
        if not records or len(records) > imports.MAX_ROWS:
            raise ValueError("Import between 1 and 100,000 molecules.")
        for record in records:
            imports.canonical_components(record["components"])
        if group_id:
            group = self.group(group_id, user, write=True)
            for record in records:
                validate_record(record, group["schema"])
        else:
            group_id = self.create_group(user, name, project_id=project_id, import_method=import_method, origin=origin)
        with self.connection() as connection:
            origin = origin or {"type": "user", "user_id": user["id"]}
            ids = [self._put(connection, group_id, record, source_key=record.get("source_key") or _id("import-", _json([group_id, index, record["name"]])), origin=origin) for index, record in enumerate(records)]
        return {"group_id": group_id, "molecule_ids": list(dict.fromkeys(ids)), "imported_count": len(records), "unique_count": len(set(ids))}

    def add_observation(self, group_id, molecule_id, namespace, run_id, metrics, *, context=None, source_key=None, created_at=None):
        clean = {str(key): value for key, value in metrics.items() if value is not None and
                 not (isinstance(value, float) and not math.isfinite(value))}
        observation_id = _id("score-", source_key or _json([group_id, molecule_id, namespace, run_id, context or {}]))
        with self.connection() as connection:
            if not connection.execute("SELECT 1 FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone():
                raise ValueError("The observation's molecule must belong to its source group.")
            connection.execute("INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET "
                               "molecule_id=excluded.molecule_id,group_id=excluded.group_id,namespace=excluded.namespace,"
                               "run_id=excluded.run_id,metrics=excluded.metrics,context=excluded.context,created_at=excluded.created_at", (
                                   observation_id, molecule_id, group_id, namespace, run_id, _json(clean), _json(context or {}), created_at or _now()))
        return observation_id

    def attach_structure(self, group_id, molecule_id, item, *, context=None, created_at=None):
        """Attach a run representation to an existing identity without editing its annotations."""
        with self.connection() as connection:
            if not connection.execute("SELECT 1 FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone():
                raise ValueError("The structure's molecule must belong to its source group.")
            path = Path(item["path"]) if item.get("path") else self._store_file(item["name"], item["content"])
            structure_context = {**(context or {}), **{key: value for key, value in item.items() if key not in {"path", "content", "name", "label"}}}
            structure_id = _id("structure-", _json([molecule_id, group_id, str(path), structure_context]))
            stamp = created_at or _now()
            connection.execute("INSERT OR IGNORE INTO molecule_structures VALUES (?,?,?,?,?,?,?,?)", (
                structure_id, molecule_id, group_id, str(path), str(item.get("label") or path.name),
                path.suffix.lower().lstrip("."), _json(structure_context), stamp))
            connection.execute("UPDATE molecule_groups SET updated_at=? WHERE id=?", (stamp, group_id))
        return structure_id

    def _visible_ids(self, user, group_id="", project_id="", scope="organization", methods=None, collection="all"):
        groups = self.groups(user, project_id=project_id, scope=scope, methods=methods, collection=collection)
        if group_id:
            groups = [row for row in groups if row["id"] == group_id]
        return [row["id"] for row in groups]

    def metric_fields(self, user, *, group_id="", project_id="", scope="organization", methods=None, collection="all"):
        ids = self._visible_ids(user, group_id, project_id, scope, methods, collection)
        # Collection selection controls which molecules are listed, while
        # scientific history may accumulate through any group the user can read.
        accessible = self._visible_ids(user, project_id=project_id, collection="all")
        with self.connection() as connection:
            rows = connection.execute("SELECT DISTINCT o.namespace,j.key FROM molecule_observations o, json_each(o.metrics) j "
                                      "WHERE o.group_id IN (SELECT value FROM json_each(?)) AND o.molecule_id IN "
                                      "(SELECT molecule_id FROM molecule_members WHERE group_id IN (SELECT value FROM json_each(?))) "
                                      "AND j.type IN ('real','integer') ORDER BY 1,2", (_json(accessible), _json(ids))).fetchall()
            common = connection.execute("SELECT DISTINCT j.key FROM molecule_observations o, "
                                        "json_each(studio_common_metrics(o.metrics,o.namespace,o.context)) j "
                                        "WHERE o.group_id IN (SELECT value FROM json_each(?)) AND o.molecule_id IN "
                                        "(SELECT molecule_id FROM molecule_members WHERE group_id IN (SELECT value FROM json_each(?))) "
                                        "ORDER BY j.key", (_json(accessible), _json(ids))).fetchall()
        return [*[PREFIX + row[0] for row in common], *[f"{row[0]}.{row[1]}" for row in rows]]

    def list_molecules(self, user, *, group_id="", project_id="", scope="organization", methods=None, collection="all",
                       query="", molecule_type="", limit=100, offset=0, sort="added", descending=True, metric_filters=None, molecule_ids=None):
        if not 1 <= int(limit) <= 100_000 or int(offset) < 0:
            raise ValueError("Invalid molecular page size or offset.")
        if metric_filters is not None and (not isinstance(metric_filters, list) or any(not isinstance(value, dict) for value in metric_filters)):
            raise ValueError("Metric filters must be a list of comparisons.")
        visible = self._visible_ids(user, group_id, project_id, scope, methods, collection)
        visibility = _json(visible)
        # Keep observations and representations from every accessible group,
        # even when the displayed collection is filtered to curated library.
        accessible = _json(self._visible_ids(user, project_id=project_id, collection="all"))
        parameters: list[Any] = [visibility]
        clause = "EXISTS (SELECT 1 FROM molecule_members mm WHERE mm.molecule_id=m.id AND mm.group_id IN (SELECT value FROM json_each(?)))"
        if query:
            clause += " AND (m.components LIKE ? OR EXISTS (SELECT 1 FROM molecule_members mm WHERE mm.molecule_id=m.id AND mm.group_id IN (SELECT value FROM json_each(?)) AND (mm.name LIKE ? OR mm.metadata LIKE ?)))"
            needle = "%" + query + "%"
            parameters.extend([needle, visibility, needle, needle])
        if molecule_type:
            clause += " AND m.type=?"
            parameters.append(molecule_type)
        if molecule_ids is not None:
            clause += " AND m.id IN (SELECT value FROM json_each(?))"
            parameters.append(_json(molecule_ids))
        def metric_expression(field, params):
            if field.startswith(PREFIX):
                key = '$.' + json.dumps(field[len(PREFIX):])
                params.extend([key, accessible, key])
                return ("(SELECT json_extract(studio_common_metrics(o.metrics,o.namespace,o.context),?) "
                        "FROM molecule_observations o WHERE o.molecule_id=m.id "
                        "AND o.group_id IN (SELECT value FROM json_each(?)) "
                        "AND json_type(studio_common_metrics(o.metrics,o.namespace,o.context),?) IN ('integer','real') "
                        "ORDER BY o.created_at DESC,o.id DESC LIMIT 1)")
            namespace, separator, key = field.partition(".")
            if not separator:
                raise ValueError("Invalid metric field.")
            params.extend(['$.' + json.dumps(key), accessible, namespace])
            return ("(SELECT json_extract(o.metrics,?) FROM molecule_observations o WHERE o.molecule_id=m.id "
                    "AND o.group_id IN (SELECT value FROM json_each(?)) AND o.namespace=? ORDER BY o.created_at DESC,o.id DESC LIMIT 1)")
        for item in metric_filters or []:
            operator = item.get("operator", ">=")
            if operator not in {">", ">=", "<", "<=", "="}:
                raise ValueError("Unsupported score comparison.")
            expression = metric_expression(str(item["field"]), parameters)
            clause += f" AND CAST({expression} AS REAL) {operator} ?"
            parameters.append(float(item["value"]))
        with self.connection() as connection:
            total = connection.execute("SELECT COUNT(*) FROM molecules m WHERE " + clause, parameters).fetchone()[0]
            ordering_params = []
            if sort == "name":
                ordering = "(SELECT mm.name FROM molecule_members mm WHERE mm.molecule_id=m.id AND mm.group_id IN (SELECT value FROM json_each(?)) ORDER BY mm.added_at DESC LIMIT 1) COLLATE NOCASE"
                ordering_params.append(visibility)
            elif sort == "added":
                ordering = "m.created_at"
            else:
                ordering = metric_expression(sort, ordering_params)
            rows = connection.execute("SELECT m.* FROM molecules m WHERE " + clause + " ORDER BY " + ordering +
                                      (" DESC" if descending else " ASC") + ",m.id LIMIT ? OFFSET ?",
                                      [*parameters, *ordering_params, int(limit), int(offset)]).fetchall()
            ids = _json([row["id"] for row in rows])
            members = connection.execute("SELECT mm.*,g.name AS group_name,g.origin,g.schema_json FROM molecule_members mm JOIN molecule_groups g ON g.id=mm.group_id "
                                         "WHERE mm.molecule_id IN (SELECT value FROM json_each(?)) AND mm.group_id IN (SELECT value FROM json_each(?)) ORDER BY mm.added_at DESC", (ids, visibility)).fetchall()
            # Accumulate only observations and representations accessible to this user,
            # including other groups when viewing a particular collection.
            structures = connection.execute("SELECT * FROM molecule_structures WHERE molecule_id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at DESC", (ids, accessible)).fetchall()
            observations = connection.execute("SELECT * FROM molecule_observations WHERE molecule_id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at DESC,id DESC", (ids, accessible)).fetchall()
            sources = connection.execute("SELECT * FROM molecule_sources WHERE molecule_id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?))", (ids, accessible)).fetchall()
        result = {row["id"]: {**dict(row), "components": json.loads(row["components"]), "groups": [], "structures": [], "observations": [], "sources": [], "scores": {}} for row in rows}
        for member in members:
            record = result[member["molecule_id"]]
            record["groups"].append({"id": member["group_id"], "name": member["group_name"]})
            if "name" not in record:
                record.update(name=member["name"], metadata=json.loads(member["metadata"]), components=json.loads(member["components"]), origin=json.loads(member["origin"]), schema=json.loads(member["schema_json"]))
        for row in structures:
            result[row["molecule_id"]]["structures"].append({**dict(row), "context": json.loads(row["context"])})
        for row in observations:
            record = result[row["molecule_id"]]
            observation = {**dict(row), "metrics": json.loads(row["metrics"]), "context": json.loads(row["context"])}
            observation["normalized_output"] = normalize_metrics(observation["metrics"], row["namespace"], observation["context"])
            record["observations"].append(observation)
            for key, value in observation["metrics"].items():
                record["scores"].setdefault(row["namespace"] + "." + key, value)
            for key, value in observation["normalized_output"]["metrics"].items():
                field = PREFIX + key
                if field not in record["scores"]:
                    record["scores"][field] = value
                    record.setdefault("score_sources", {})[field] = {"observation_id": row["id"], "run_id": row["run_id"],
                        "namespace": row["namespace"], "created_at": row["created_at"],
                        **observation["normalized_output"]["provenance"][key]}
        for row in sources:
            result[row["molecule_id"]]["sources"].append({"source_key": row["source_key"], **json.loads(row["provenance"])})
        from mn_molecular_design_studio.molecule_schema import evaluate
        for record in result.values():
            for field in record.get("schema", {}).get("fields", []):
                if field.get("formula"):
                    record["metadata"][field["name"]] = evaluate(field["formula"], {**record["metadata"], **record["scores"]})
        return {"items": list(result.values()), "total": total, "next_offset": offset + len(rows) if offset + len(rows) < total else None}

    def molecules_by_ids(self, user, ids: list[str]) -> list[dict]:
        if not ids:
            raise ValueError("Select at least one molecule.")
        visible = self._visible_ids(user)
        with self.connection() as connection:
            rows = connection.execute("SELECT m.id,mm.name,mm.components,mm.metadata,mm.group_id FROM molecules m JOIN molecule_members mm ON mm.molecule_id=m.id "
                                      "WHERE m.id IN (SELECT value FROM json_each(?)) AND mm.group_id IN (SELECT value FROM json_each(?)) ORDER BY mm.added_at DESC", (_json(ids), _json(visible))).fetchall()
        result = {}
        for row in rows:
            result.setdefault(row["id"], {**dict(row), "components": json.loads(row["components"]), "metadata": json.loads(row["metadata"])})
        if set(ids) != set(result):
            raise PermissionError("Some selected molecules are missing or inaccessible.")
        return [result[value] for value in dict.fromkeys(ids)]

    def add_to_group(self, user, group_id, ids):
        group = self.group(group_id, user, write=True)
        rows = self.molecules_by_ids(user, ids)
        from mn_molecular_design_studio.molecule_schema import validate_record
        for row in rows:
            validate_record({"components": row["components"], "metadata": row["metadata"]}, group["schema"])
        with self.connection() as connection:
            for row in rows:
                # Reuse identity exactly, including unsequenced backbone identity.
                connection.execute("INSERT OR IGNORE INTO molecule_members VALUES (?,?,?,?,?,?)", (
                    group_id, row["id"], row["name"], _json(row["components"]), _json(row["metadata"]), _now()))
            connection.execute("UPDATE molecule_groups SET updated_at=? WHERE id=?", (_now(), group_id))
        return len(rows)

    def remove_from_group(self, user, group_id, ids):
        """Remove selected molecule rows and their group-scoped data only."""
        self.group(group_id, user, write=True)
        selected = list(dict.fromkeys(str(value) for value in ids or []))
        if not selected:
            raise ValueError("Select at least one molecule.")
        with self.connection() as connection:
            present = {row[0] for row in connection.execute(
                "SELECT molecule_id FROM molecule_members WHERE group_id=? AND molecule_id IN (SELECT value FROM json_each(?))",
                (group_id, _json(selected)),
            ).fetchall()}
            if present != set(selected):
                raise PermissionError("Some selected molecules are not in this group or are inaccessible.")
            for table in ("molecule_sources", "molecule_structures", "molecule_observations"):
                connection.execute(
                    f"DELETE FROM {table} WHERE group_id=? AND molecule_id IN (SELECT value FROM json_each(?))",
                    (group_id, _json(selected)),
                )
            connection.execute(
                "DELETE FROM molecule_members WHERE group_id=? AND molecule_id IN (SELECT value FROM json_each(?))",
                (group_id, _json(selected)),
            )
            connection.execute("UPDATE molecule_groups SET updated_at=? WHERE id=?", (_now(), group_id))
        return len(selected)

    def update_metadata(self, user, group_id, molecule_id, metadata, *, name=None):
        group = self.group(group_id, user, write=True)
        with self.connection() as connection:
            row = connection.execute("SELECT metadata FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone()
            if not row:
                raise ValueError("The molecule is not in this group.")
            values = {**json.loads(row[0]), **metadata}
            from mn_molecular_design_studio.molecule_schema import validate_record
            components = json.loads(connection.execute("SELECT components FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone()[0])
            validate_record({"metadata": values, "components": components}, group["schema"])
            connection.execute("UPDATE molecule_members SET metadata=? WHERE group_id=? AND molecule_id=?", (_json(values), group_id, molecule_id))
            if name is not None:
                if not name.strip():
                    raise ValueError("A molecule name cannot be empty.")
                connection.execute("UPDATE molecule_members SET name=? WHERE group_id=? AND molecule_id=?", (name.strip(), group_id, molecule_id))

    def set_schema(self, user, group_id, schema):
        self.group(group_id, user, write=True)
        from mn_molecular_design_studio.molecule_schema import validate_schema, validate_record
        validate_schema(schema)
        with self.connection() as connection:
            members = connection.execute("SELECT components,metadata FROM molecule_members WHERE group_id=?", (group_id,)).fetchall()
            for member in members:
                validate_record({"components": json.loads(member[0]), "metadata": json.loads(member[1])}, schema)
            connection.execute("UPDATE molecule_groups SET schema_json=?,updated_at=? WHERE id=?", (_json(schema), _now(), group_id))

    def set_chain_role(self, user, group_id, molecule_id, chain, role):
        self.group(group_id, user, write=True)
        if role not in {"", "lead", "target", "cofactor", "heavy", "light"}:
            raise ValueError("Unknown chain role.")
        with self.connection() as connection:
            row = connection.execute("SELECT components FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone()
            if not row:
                raise ValueError("The molecule is not in this group.")
            components = json.loads(row[0])
            if chain not in {component["chain"] for component in components}:
                raise ValueError("This chain is not part of the molecule.")
            for component in components:
                if role == "lead" and component.get("role") == "lead":
                    component["role"] = ""
                if component["chain"] == chain:
                    component["role"] = role
            connection.execute("UPDATE molecule_members SET components=? WHERE group_id=? AND molecule_id=?", (_json(components), group_id, molecule_id))

    def compound_snapshot(self, user, ids, *, name="Selected compounds", project_id="") -> dict:
        from mn_molecular_design_studio import ligand_compounds
        from mn_molecular_design_studio.auth import require_project_access
        require_project_access(user, project_id or None)
        rows = self.molecules_by_ids(user, ids)
        compounds, seen = [], set()
        for row in rows:
            components = row["components"]
            if len(components) != 1 or components[0]["type"] != "small_molecule" or not components[0].get("smiles"):
                raise ValueError("This selection must contain only compounds with known SMILES.")
            compound_id = re_safe(row["name"])[:80] or row["id"]
            if compound_id.casefold() in seen:
                compound_id += "-" + row["id"][-8:]
            seen.add(compound_id.casefold())
            compounds.append({"compound_id": compound_id, "smiles": components[0]["smiles"],
                              "properties": {**row["metadata"], "studio_molecule_id": row["id"]}})
        digest = hashlib.sha256(_json([user["id"], project_id, compounds]).encode()).hexdigest()
        from mn_molecular_design_studio.auth import user_input_root
        path = user_input_root(user["id"]) / "database_selections" / (digest + ".json")
        if not path.is_file():
            ligand_compounds.write_compound_set(path, name=name, compounds=compounds, owner_id=user["id"],
                project_id=project_id, created_at=_now(), compound_set_id="selection-" + digest[:24])
        return {"path": str(path), "name": name, "compound_count": len(compounds), "molecule_ids": [row["id"] for row in rows]}

    def compound_options(self, user, *, collection="library"):
        result = []
        for group in self.groups(user, collection=collection):
            if not group["molecule_count"]:
                continue
            with self.connection() as connection:
                rows = connection.execute("SELECT molecule_id,components FROM molecule_members WHERE group_id=? ORDER BY added_at,molecule_id", (group["id"],)).fetchall()
            if any(len(c := json.loads(row["components"])) != 1 or not c[0].get("smiles") for row in rows):
                continue
            snapshot = self.compound_snapshot(user, [row["molecule_id"] for row in rows], name=group["name"], project_id=group["project_id"])
            result.append({**snapshot, "group_id": "library:" + group["id"]})
        return result

    def structure_file(self, user, structure_id, *, materialize=False):
        with self.connection() as connection:
            row = connection.execute("SELECT * FROM molecule_structures WHERE id=?", (structure_id,)).fetchone()
        if not row:
            raise FileNotFoundError("Structure not found.")
        self.group(row["group_id"], user)
        if str(row["path"]).startswith("cluster-library://"):
            from mn_molecular_design_studio.molecular_database_peer import fetch_structure
            path = fetch_structure(row, user)
        else:
            path = Path(row["path"])
        if not path.is_file():
            raise FileNotFoundError("The source structure is unavailable on this workstation.")
        if not materialize:
            return path
        from mn_molecular_design_studio.auth import user_input_root
        content = path.read_bytes()
        destination = user_input_root(user["id"]) / "database_structures" / (hashlib.sha256(content).hexdigest() + path.suffix)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists():
            destination.write_bytes(content)
        return destination


def save_legacy_groups(path: Path, groups: list[dict], user: dict):
    """Compatibility writes preserve other owners and concurrent additions."""
    import fcntl
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_suffix(".lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            existing = json.loads(path.read_text())
        except (OSError, ValueError):
            existing = []
        by_id = {row["id"]: row for row in existing if isinstance(row, dict) and row.get("id")}
        for group in groups:
            previous = by_id.get(group["id"], {})
            if user.get("role") != "admin" and (group.get("owner_id") != user["id"] or (previous and previous.get("owner_id") != user["id"])):
                raise PermissionError("Cannot update another user's molecule group.")
            members = {_json(row): row for row in [*(previous.get("members") or []), *(group.get("members") or [])]}
            by_id[group["id"]] = {**previous, **group, "members": list(members.values())}
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as temporary:
            json.dump(list(by_id.values()), temporary, indent=2)
            temporary_path = Path(temporary.name)
        os.replace(temporary_path, path)


def re_safe(value):
    import re
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(value)).strip("-._")
