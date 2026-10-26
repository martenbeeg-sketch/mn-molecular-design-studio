"""Durable molecular identities, collections, structures and score observations.

Group membership grants visibility. Scientific observations and representations
retain their own group/project context, even when molecular identities coincide.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
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
_SMOKE_TEST_GROUP = re.compile(r"(?<![a-z0-9])(smoke|test|testing)(?![a-z0-9])", re.IGNORECASE)
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
            CREATE INDEX IF NOT EXISTS molecular_source_lookup ON molecule_sources(molecule_id,group_id);
            CREATE INDEX IF NOT EXISTS molecular_observation_lookup ON molecule_observations(molecule_id,group_id);
            CREATE TABLE IF NOT EXISTS molecular_ingestions (
                source_key TEXT PRIMARY KEY, signature TEXT NOT NULL, updated_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS molecule_import_batches (
                id TEXT PRIMARY KEY,
                group_id TEXT NOT NULL REFERENCES molecule_groups(id),
                created_by TEXT NOT NULL,
                source_name TEXT NOT NULL,
                source_sha256 TEXT NOT NULL,
                source_path TEXT NOT NULL DEFAULT '',
                mapping_json TEXT NOT NULL,
                summary_json TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS molecular_import_batches_group ON molecule_import_batches(group_id,created_at DESC);
            CREATE TABLE IF NOT EXISTS molecule_import_rows (
                batch_id TEXT NOT NULL REFERENCES molecule_import_batches(id) ON DELETE CASCADE,
                row_number INTEGER NOT NULL,
                source_id TEXT NOT NULL DEFAULT '',
                name TEXT NOT NULL DEFAULT '',
                raw_smiles TEXT NOT NULL DEFAULT '',
                canonical_smiles TEXT NOT NULL DEFAULT '',
                connectivity_smiles TEXT NOT NULL DEFAULT '',
                status TEXT NOT NULL,
                duplicate_of_row INTEGER,
                stereo_variant_of_row INTEGER,
                warnings_json TEXT NOT NULL DEFAULT '[]',
                errors_json TEXT NOT NULL DEFAULT '[]',
                descriptors_json TEXT NOT NULL DEFAULT '{}',
                source_data_json TEXT NOT NULL DEFAULT '{}',
                molecule_id TEXT REFERENCES molecules(id),
                matched_existing INTEGER NOT NULL DEFAULT 0,
                PRIMARY KEY(batch_id,row_number)
            );
            CREATE INDEX IF NOT EXISTS molecular_import_rows_status ON molecule_import_rows(batch_id,status,row_number);
            CREATE INDEX IF NOT EXISTS molecular_import_rows_molecule ON molecule_import_rows(molecule_id);
            CREATE TABLE IF NOT EXISTS external_library_records (
                id TEXT PRIMARY KEY, record_json TEXT NOT NULL,
                created_at TEXT NOT NULL, updated_at TEXT NOT NULL
            );
        ''')
        from mn_molecular_design_studio.evaluations import install_schema
        install_schema(connection)
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

    @staticmethod
    def _is_smoke_test_group(group) -> bool:
        origin = group.get("origin") if isinstance(group, dict) else None
        if isinstance(origin, str):
            try:
                origin = json.loads(origin)
            except (TypeError, ValueError):
                origin = {}
        origin = origin if isinstance(origin, dict) else {}
        if any(origin.get(key) is True for key in ("is_test", "test", "smoke", "is_smoke_test")):
            return True
        labels = [str(group.get("name") or "")]
        labels.extend(str(origin.get(key) or "") for key in ("name", "label", "run_name", "campaign", "workflow", "description"))
        return any(_SMOKE_TEST_GROUP.search(label) for label in labels)

    def groups(self, user: dict, *, project_id="", scope="organization", methods=None, collection="all",
               hide_smoke_test=False) -> list[dict]:
        if collection not in {"all", "library"}:
            raise ValueError("Collection must be 'all' or 'library'.")
        if scope not in {"organization", "all", "mine", "shared"}:
            raise ValueError("Scope must be 'mine', 'shared', or 'all'.")
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
            if scope == "mine" and (row["owner_id"] != user["id"] or bool(row["project_id"])):
                continue
            if scope == "shared" and not row["project_id"]:
                continue
            if hide_smoke_test and self._is_smoke_test_group(item):
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

    def set_group_project(self, user: dict, group_id: str, project_id: str) -> dict:
        """Move an existing molecule collection into a project without changing identities."""
        from mn_molecular_design_studio.auth import require_project_access
        from mn_molecular_design_studio.workspace import validate_project_id

        selected_project = validate_project_id(project_id)
        if not selected_project:
            raise ValueError("Choose a valid project for this molecule collection.")
        require_project_access(user, selected_project)
        group = self.group(group_id, user, write=True)
        stamp = _now()
        with self.connection() as connection:
            connection.execute(
                "UPDATE molecule_groups SET project_id=?,updated_at=? WHERE id=?",
                (selected_project, stamp, group_id),
            )
        return self.group(group_id, user)

    def _put(self, connection, group_id: str, record: dict, *, source_key="", origin=None, stamp=None, preserve_existing=False):
        components, fingerprint, kind = imports.canonical_components(record["components"])
        if record.get("identity_context"):
            fingerprint = hashlib.sha256((fingerprint + _json(record["identity_context"])).encode()).hexdigest()
        molecule_id = "mol-" + fingerprint[:32]
        stamp = stamp or _now()
        connection.execute("INSERT OR IGNORE INTO molecules VALUES (?,?,?,?,?)",
                           (molecule_id, fingerprint, kind, _json(components), stamp))
        existing = connection.execute("SELECT name,metadata FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone()
        supplied_metadata = record.get("metadata") or {}
        if preserve_existing:
            metadata = dict(supplied_metadata)
            if existing:
                for key, value in json.loads(existing["metadata"]).items():
                    if value not in (None, ""):
                        metadata[key] = value
            saved_name = existing["name"] if existing else str(record.get("name") or molecule_id)[:250]
        else:
            metadata = {**(json.loads(existing["metadata"]) if existing else {}), **supplied_metadata}
            saved_name = str(record.get("name") or molecule_id)[:250]
        connection.execute("INSERT INTO molecule_members VALUES (?,?,?,?,?,?) ON CONFLICT(group_id,molecule_id) "
                           "DO UPDATE SET metadata=excluded.metadata", (
                               group_id, molecule_id, saved_name, _json(components), _json(metadata), stamp))
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
        if suffix not in {".pdb", ".ent", ".cif", ".mmcif", ".pdbqt", ".sdf", ".mol", ".fasta", ".fa"}:
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

    def _store_import_source(self, name: str, content: bytes) -> Path:
        """Store the original tabular/SMILES source by content hash."""
        suffix = Path(name).suffix.lower()
        if suffix not in {".csv", ".tsv", ".txt", ".xlsx", ".xlsm", ".xls", ".smi", ".smiles"}:
            suffix = ".bin"
        digest = hashlib.sha256(content).hexdigest()
        directory = self.root / "molecular_import_sources" / digest[:2]
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / (digest + suffix)
        if not path.exists():
            with tempfile.NamedTemporaryFile(dir=directory, delete=False) as temporary:
                temporary.write(content)
                temporary_path = Path(temporary.name)
            os.replace(temporary_path, path)
        return path

    def group_identity_map(self, user, group_id):
        """Return canonical structure fingerprints already present in a group."""
        self.group(group_id, user)
        with self.connection() as connection:
            rows = connection.execute(
                "SELECT m.fingerprint,m.id,mm.name FROM molecules m JOIN molecule_members mm ON mm.molecule_id=m.id WHERE mm.group_id=?",
                (group_id,),
            ).fetchall()
        return {row["fingerprint"]: {"id": row["id"], "name": row["name"]} for row in rows}

    def external_library_record(self, library_id):
        """Read metadata for an external compound library; compounds stay in their source files."""
        with self.connection() as connection:
            row = connection.execute("SELECT record_json FROM external_library_records WHERE id=?", (str(library_id),)).fetchone()
        if not row:
            raise KeyError("External library metadata record not found.")
        return json.loads(row["record_json"])

    def external_library_records(self):
        with self.connection() as connection:
            rows = connection.execute("SELECT record_json FROM external_library_records ORDER BY updated_at DESC,id").fetchall()
        return [json.loads(row["record_json"]) for row in rows]

    def save_external_library_record(self, record):
        """Store one metadata-only library record in the app's SQLite database."""
        library_id = str(record.get("id") or "")
        if not library_id.startswith("lib-"):
            raise ValueError("External library records require a valid library ID.")
        stamp = _now()
        with self.connection() as connection:
            existing = connection.execute("SELECT created_at FROM external_library_records WHERE id=?", (library_id,)).fetchone()
            created_at = str(record.get("created_at") or (existing["created_at"] if existing else stamp))
            record = {**record, "created_at": created_at, "updated_at": stamp}
            connection.execute(
                "INSERT INTO external_library_records VALUES (?,?,?,?) "
                "ON CONFLICT(id) DO UPDATE SET record_json=excluded.record_json,updated_at=excluded.updated_at",
                (library_id, _json(record), created_at, stamp),
            )
        return record

    def _save_import_audit(self, connection, user, group_id, audit, records, molecule_ids, source_content=None, batch_id=None):
        if not audit:
            return
        batch_id = str(batch_id or uuid.uuid4().hex)
        if not re.fullmatch(r"[A-Za-z0-9_-]{1,100}", batch_id):
            raise ValueError("Invalid import batch ID.")
        source_name = Path(str(audit.get("source_name") or "source file")).name[:250]
        source_sha256 = str(audit.get("source_sha256") or "")
        if not re.fullmatch(r"[0-9a-f]{64}", source_sha256):
            raise ValueError("The import source checksum is invalid.")
        if source_content is not None and hashlib.sha256(source_content).hexdigest() != source_sha256:
            raise ValueError("The uploaded source no longer matches its validation report.")
        if source_content is not None and len(source_content) > imports.MAX_BYTES:
            raise ValueError("The original import file exceeds 100 MB.")
        mapping = audit.get("mapping") or {}
        summary = audit.get("summary") or {}
        rows = audit.get("rows") or []
        if len(rows) > imports.MAX_ROWS:
            raise ValueError("An import report can contain at most 100,000 rows.")
        source_path = ""
        if source_content is not None:
            source_path = str(self._store_import_source(source_name, source_content))
        molecule_by_smiles = {}
        for record, molecule_id in zip(records, molecule_ids):
            smiles = next((component.get("smiles") for component in record.get("components", [])
                           if component.get("type") == "small_molecule"), None)
            if smiles:
                molecule_by_smiles.setdefault(smiles, molecule_id)
        batch_key = "import-batch-" + batch_id
        created_at = _now()
        connection.execute(
            "INSERT INTO molecule_import_batches VALUES (?,?,?,?,?,?,?,?,?)",
            (batch_key, group_id, str(user["id"]), source_name, source_sha256, source_path,
             _json(mapping), _json(summary), created_at),
        )
        for row in rows:
            row_number = int(row.get("row_number") or 0)
            if row_number < 1:
                raise ValueError("Invalid source row number in the import report.")
            canonical = str(row.get("canonical_smiles") or "")
            molecule_id = molecule_by_smiles.get(canonical)
            if molecule_id and row.get("status") in {"valid", "duplicate"}:
                row["molecule_id"] = molecule_id
            provenance = {
                "type": "table_import_row",
                "import_batch_id": batch_key,
                "source_name": source_name,
                "source_sha256": source_sha256,
                "source_row": row_number,
                "source_id": str(row.get("source_id") or ""),
                "validation_status": str(row.get("status") or ""),
                "raw_smiles": str(row.get("raw_smiles") or ""),
                "canonical_smiles": canonical,
                "warnings": row.get("warnings") or [],
                "errors": row.get("errors") or [],
            }
            if molecule_id:
                source_key = f"import-row:{batch_id}:{row_number}"
                connection.execute(
                    "INSERT INTO molecule_sources VALUES (?,?,?,?) ON CONFLICT(source_key) DO UPDATE SET "
                    "molecule_id=excluded.molecule_id,group_id=excluded.group_id,provenance=excluded.provenance",
                    (source_key, molecule_id, group_id, _json(provenance)),
                )
            source_data = row.get("source_data") or {}
            descriptors = dict(row.get("descriptors") or {})
            descriptors["normalized_compound_id"] = row.get("compound_id") or ""
            descriptors["annotation_stripped_candidate"] = row.get("annotation_stripped_candidate") or ""
            descriptors["candidate_status"] = row.get("candidate_status") or ""
            connection.execute(
                "INSERT INTO molecule_import_rows VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (batch_key, row_number, str(row.get("source_id") or "")[:250], str(row.get("name") or "")[:250],
                 str(row.get("raw_smiles") or ""), canonical, str(row.get("connectivity_smiles") or ""),
                 str(row.get("status") or "unknown"), row.get("duplicate_of_row"), row.get("stereo_variant_of_row"),
                 _json(row.get("warnings") or []), _json(row.get("errors") or []), _json(descriptors),
                 _json(source_data), molecule_id, int(bool(row.get("matched_existing")))),
            )

    def import_batches(self, user, group_id):
        self.group(group_id, user)
        with self.connection() as connection:
            rows = connection.execute(
                "SELECT * FROM molecule_import_batches WHERE group_id=? ORDER BY created_at DESC,id DESC", (group_id,)
            ).fetchall()
        return [{**dict(row), "mapping": json.loads(row["mapping_json"]),
                 "summary": json.loads(row["summary_json"])} for row in rows]

    def import_batch_rows(self, user, batch_id):
        with self.connection() as connection:
            batch = connection.execute("SELECT group_id FROM molecule_import_batches WHERE id=?", (batch_id,)).fetchone()
            if not batch:
                raise ValueError("Import report not found.")
            self.group(batch["group_id"], user)
            rows = connection.execute(
                "SELECT * FROM molecule_import_rows WHERE batch_id=? ORDER BY row_number", (batch_id,)
            ).fetchall()
        result = []
        for row in rows:
            descriptors = json.loads(row["descriptors_json"])
            result.append({
                **dict(row),
                "warnings": json.loads(row["warnings_json"]),
                "errors": json.loads(row["errors_json"]),
                "compound_id": descriptors.get("normalized_compound_id", ""),
                "annotation_stripped_candidate": descriptors.get("annotation_stripped_candidate", ""),
                "candidate_status": descriptors.get("candidate_status", ""),
                "descriptors": descriptors,
                "source_data": json.loads(row["source_data_json"]),
            })
        return result

    def import_source_file(self, user, batch_id):
        with self.connection() as connection:
            batch = connection.execute("SELECT group_id,source_path FROM molecule_import_batches WHERE id=?", (batch_id,)).fetchone()
        if not batch:
            raise ValueError("Import source not found.")
        self.group(batch["group_id"], user)
        path = Path(batch["source_path"])
        if not path.is_file():
            raise FileNotFoundError("The original import file is unavailable on this workstation.")
        return path

    def import_records(self, user: dict, records: list[dict], *, name="Imported molecules", project_id="", group_id=None, import_method="user_uploads", origin=None, import_audit=None, source_content=None) -> dict:
        from mn_molecular_design_studio.molecule_schema import validate_record
        if (not records and not import_audit) or len(records) > imports.MAX_ROWS:
            raise ValueError("Import between 1 and 100,000 molecules, or save a validation report.")
        for record in records:
            imports.canonical_components(record["components"])
        if group_id:
            group = self.group(group_id, user, write=True)
            for record in records:
                validate_record(record, group["schema"])
        else:
            group_id = self.create_group(user, name, project_id=project_id, import_method=import_method, origin=origin)
        audit_batch_id = uuid.uuid4().hex if import_audit else None
        with self.connection() as connection:
            origin = origin or {"type": "user", "user_id": user["id"]}
            ids = []
            for index, record in enumerate(records):
                source_key = record.get("source_key") or _id("import-", _json([group_id, index, record["name"]]))
                if audit_batch_id and record.get("source_row_number"):
                    source_key = f"import-row:{audit_batch_id}:{int(record['source_row_number'])}"
                ids.append(self._put(connection, group_id, record, source_key=source_key, origin=origin,
                                     preserve_existing=bool(import_audit)))
            self._save_import_audit(connection, user, group_id, import_audit, records, ids,
                                    source_content=source_content, batch_id=audit_batch_id)
        return {"group_id": group_id, "molecule_ids": list(dict.fromkeys(ids)), "imported_count": len(records), "unique_count": len(set(ids)),
                "import_batch_id": "import-batch-" + audit_batch_id if audit_batch_id else None}

    def add_observation(self, group_id, molecule_id, namespace, run_id, metrics, *, context=None, source_key=None, created_at=None):
        from mn_molecular_design_studio.evaluations import file_digest
        context = dict(context or {})
        if context.get("target_path"):
            context.setdefault("target_structure_digest", file_digest(context["target_path"]))
        clean = {str(key): value for key, value in metrics.items() if value is not None and
                 not (isinstance(value, float) and not math.isfinite(value))}
        observation_id = _id("score-", source_key or _json([group_id, molecule_id, namespace, run_id, context or {}]))
        with self.connection() as connection:
            if not connection.execute("SELECT 1 FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone():
                raise ValueError("The observation's molecule must belong to its source group.")
            if not context.get("representation_id") and context.get("is_representative") is not False:
                structures = connection.execute("SELECT * FROM molecule_structures WHERE molecule_id=? AND group_id=? ORDER BY created_at DESC,id DESC",
                                                (molecule_id, group_id)).fetchall()
                for structure in structures:
                    provenance = json.loads(structure["context"])
                    if provenance.get("run_id") != run_id or (context.get("candidate_id") and provenance.get("candidate_id") != context["candidate_id"]):
                        continue
                    context["representation_id"] = structure["id"]
                    context["structure_digest"] = file_digest(structure["path"])
                    break
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

    def make_structure_portable(self, user, structure_id):
        """Keep a structure in the Database store so peer workstations can replicate it."""
        with self.connection() as connection:
            row = connection.execute("SELECT * FROM molecule_structures WHERE id=?", (structure_id,)).fetchone()
        if not row:
            raise FileNotFoundError("Structure not found.")
        self.group(row["group_id"], user, write=True)
        source = Path(row["path"])
        if not source.is_file():
            raise FileNotFoundError("The source structure is unavailable on this workstation.")
        content = source.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        destination = self._store_file(source.name, content)
        context = json.loads(row["context"])
        context.update(portable=True, portable_sha256=digest)
        with self.connection() as connection:
            connection.execute("UPDATE molecule_structures SET path=?,context=? WHERE id=?",
                               (str(destination), _json(context), structure_id))
        return {"structure_id": structure_id, "path": str(destination), "sha256": digest}

    def _visible_ids(self, user, group_id="", project_id="", scope="organization", methods=None, collection="all",
                     hide_smoke_test=False):
        groups = self.groups(user, project_id=project_id, scope=scope, methods=methods, collection=collection,
                             hide_smoke_test=hide_smoke_test)
        if group_id:
            groups = [row for row in groups if row["id"] == group_id]
        return [row["id"] for row in groups]

    def metric_fields(self, user, *, group_id="", project_id="", scope="organization", methods=None, collection="all",
                      hide_smoke_test=False):
        ids = self._visible_ids(user, group_id, project_id, scope, methods, collection, hide_smoke_test)
        # Collection selection controls which molecules are listed, while
        # scientific history may accumulate through any group the user can read.
        accessible = self._visible_ids(user, project_id=project_id, collection="all", hide_smoke_test=hide_smoke_test)
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
                       query="", molecule_type="", limit=100, offset=0, sort="added", descending=True, metric_filters=None, molecule_ids=None,
                       evaluation_filters=None, hide_smoke_test=False):
        from mn_molecular_design_studio.evaluations import rebuild_missing, selection_clause, project_observation
        if group_id and evaluation_filters is None:
            group = next((g for g in self.groups(user, collection="all") if g["id"] == group_id), {})
            saved = (group.get("origin") or {}).get("comparison_spec") or {}
            evaluation_filters = saved.get("evaluation_filters")
        if not 1 <= int(limit) <= 100_000 or int(offset) < 0:
            raise ValueError("Invalid molecular page size or offset.")
        if metric_filters is not None and (not isinstance(metric_filters, list) or any(not isinstance(value, dict) for value in metric_filters)):
            raise ValueError("Metric filters must be a list of comparisons.")
        visible = self._visible_ids(user, group_id, project_id, scope, methods, collection, hide_smoke_test)
        visibility = _json(visible)
        # Keep observations and representations from every accessible group,
        # even when the displayed collection is filtered to curated library.
        accessible = _json(self._visible_ids(user, project_id=project_id, collection="all", hide_smoke_test=hide_smoke_test))
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
        def selected_evaluation(params):
            params.append(accessible)
            selected = selection_clause(evaluation_filters, params)
            return ("(SELECT e.observation_id FROM molecule_evaluations e WHERE e.molecule_id=m.id "
                    "AND e.group_id IN (SELECT value FROM json_each(?)) AND " + selected +
                    " ORDER BY e.preferred DESC,e.created_at DESC,e.observation_id DESC LIMIT 1)")

        def metric_expression(field, params):
            if "." not in field:
                raise ValueError("Invalid metric field.")
            params.append(field)
            return "(SELECT v.value_number FROM molecule_metric_values v WHERE v.field=? AND v.observation_id=" + selected_evaluation(params) + ")"
        if evaluation_filters:
            clause += " AND " + selected_evaluation(parameters) + " IS NOT NULL"
        for item in metric_filters or []:
            operator = item.get("operator", ">=")
            if operator not in {">", ">=", "<", "<=", "="}:
                raise ValueError("Unsupported score comparison.")
            expression = metric_expression(str(item["field"]), parameters)
            clause += f" AND CAST({expression} AS REAL) {operator} ?"
            parameters.append(float(item["value"]))
        with self.connection() as connection:
            rebuild_missing(connection)
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
                                      (" DESC" if descending else " ASC") + " NULLS LAST,m.id LIMIT ? OFFSET ?",
                                      [*parameters, *ordering_params, int(limit), int(offset)]).fetchall()
            ids = _json([row["id"] for row in rows])
            members = connection.execute("SELECT mm.*,g.name AS group_name,g.origin,g.schema_json FROM molecule_members mm JOIN molecule_groups g ON g.id=mm.group_id "
                                         "WHERE mm.molecule_id IN (SELECT value FROM json_each(?)) AND mm.group_id IN (SELECT value FROM json_each(?)) ORDER BY mm.added_at DESC", (ids, visibility)).fetchall()
            # Accumulate only observations and representations accessible to this user,
            # including other groups when viewing a particular collection.
            structures = connection.execute("SELECT * FROM molecule_structures WHERE molecule_id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at DESC", (ids, accessible)).fetchall()
            observations = connection.execute("SELECT * FROM molecule_observations WHERE molecule_id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at DESC,id DESC", (ids, accessible)).fetchall()
            sources = connection.execute("SELECT * FROM molecule_sources WHERE molecule_id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?))", (ids, accessible)).fetchall()
            selection_parameters = []
            selection_sql = selected_evaluation(selection_parameters)
            selected = dict(connection.execute("SELECT m.id," + selection_sql + " FROM molecules m WHERE m.id IN (SELECT value FROM json_each(?))",
                                               [*selection_parameters, ids]).fetchall())
            evaluations = {row["observation_id"]: json.loads(row["descriptor"]) for row in connection.execute(
                "SELECT * FROM molecule_evaluations WHERE observation_id IN (SELECT value FROM json_each(?))", (_json(list(selected.values())),)).fetchall()}
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
            if selected.get(row["molecule_id"]) == row["id"]:
                record["scores"], record["score_sources"], _ = project_observation(observation)
                record["selected_evaluation"] = evaluations.get(row["id"])
        for row in sources:
            result[row["molecule_id"]]["sources"].append({"source_key": row["source_key"], **json.loads(row["provenance"])})
        from mn_molecular_design_studio.molecule_schema import evaluate
        for record in result.values():
            record["evaluation_count"] = len(record["observations"])
            record["representations"] = list(record["structures"])
            evaluation = record.get("selected_evaluation")
            if evaluation:
                if evaluation["representation_id"]:
                    record["structures"] = [s for s in record["structures"] if s["id"] == evaluation["representation_id"]]
                else:
                    matched = [s for s in record["structures"] if s["context"].get("run_id") == evaluation["run_id"]
                               and (not evaluation["candidate_id"] or s["context"].get("candidate_id") == evaluation["candidate_id"])]
                    record["structures"] = matched[:1]
            for field in record.get("schema", {}).get("fields", []):
                if field.get("formula"):
                    record["metadata"][field["name"]] = evaluate(field["formula"], {**record["metadata"], **record["scores"]})
        return {"items": list(result.values()), "total": total, "next_offset": offset + len(rows) if offset + len(rows) < total else None}

    def evaluation_options(self, user, *, group_id="", project_id="", collection="all", scope="organization",
                           hide_smoke_test=False):
        from mn_molecular_design_studio.evaluations import rebuild_missing
        visible = self._visible_ids(user, group_id=group_id, project_id=project_id, scope=scope, collection=collection,
                                    hide_smoke_test=hide_smoke_test)
        accessible = self._visible_ids(user, project_id=project_id, hide_smoke_test=hide_smoke_test)
        with self.connection() as connection:
            rebuild_missing(connection)
            rows = connection.execute("SELECT DISTINCT e.target_id,e.evaluator,e.run_id FROM molecule_evaluations e "
                "WHERE e.group_id IN (SELECT value FROM json_each(?)) AND e.molecule_id IN "
                "(SELECT molecule_id FROM molecule_members WHERE group_id IN (SELECT value FROM json_each(?))) ORDER BY 1,2,3",
                (_json(accessible), _json(visible))).fetchall()
        return {key: sorted({row[key] for row in rows if row[key]}) for key in ("target_id", "evaluator", "run_id")}

    def evaluation_rows(self, user, *, evaluation_filters=None, **filters):
        from mn_molecular_design_studio.evaluations import evaluation_descriptor, project_observation
        if filters.get("group_id") and evaluation_filters is None:
            saved = self.group(filters["group_id"], user)["origin"].get("comparison_spec") or {}
            evaluation_filters = saved.get("evaluation_filters")
        records = self.list_molecules(user, limit=100000, evaluation_filters=evaluation_filters, **filters)["items"]
        result = []
        for record in records:
            for observation in record["observations"]:
                evaluation = evaluation_descriptor(observation)
                if any(value and (evaluation["id"] not in value if key == "evaluation_ids" else evaluation.get("id" if key == "evaluation_id" else key) != value)
                       for key, value in (evaluation_filters or {}).items()):
                    continue
                scores, sources, _ = project_observation(observation)
                result.append({"molecule_id": record["id"], "name": record["name"], "metadata": record.get("metadata", {}),
                               "components": record["components"], "evaluation": evaluation, "scores": scores, "score_sources": sources})
        return result

    def save_comparison(self, user, name, molecule_ids, *, project_id="", evaluation_filters=None, settings=None):
        from mn_molecular_design_studio.comparisons import compare, prepare_rows, compatibility_warnings
        if not isinstance(molecule_ids, list) or not molecule_ids or len(molecule_ids) > 100000:
            raise ValueError("Select between 1 and 100,000 molecules.")
        if not isinstance(settings or {}, dict) or set(settings or {}) - {"fields", "directions", "weights", "thresholds", "method", "combine_same_subject", "evaluator_runs"}:
            raise ValueError("Unsupported comparison settings.")
        rows = self.list_molecules(user, molecule_ids=molecule_ids, project_id=project_id,
                                  evaluation_filters=evaluation_filters, limit=100000)["items"]
        if set(molecule_ids) != {row["id"] for row in rows}:
            raise ValueError("Some selected molecules have no accessible matching evaluation.")
        evaluations = self.evaluation_rows(user, molecule_ids=molecule_ids, project_id=project_id, evaluation_filters=evaluation_filters)
        compare(evaluations, **(settings or {}))
        prepared = prepare_rows(evaluations, combine_same_subject=(settings or {}).get("combine_same_subject", False), evaluator_runs=(settings or {}).get("evaluator_runs"))
        frozen = sorted({evaluation_id for row in prepared for evaluation_id in (row["evaluation"].get("source_evaluation_ids") or [row["evaluation"]["id"]])})
        if not frozen:
            raise ValueError("Select at least one completed evaluation to save a comparison.")
        if set(molecule_ids) != {row["molecule_id"] for row in prepared}:
            raise ValueError("Some selected molecules have no evaluation in the contributing runs.")
        warnings = compatibility_warnings(prepared, (settings or {}).get("fields") or [])
        if warnings:
            raise ValueError("; ".join(warnings))
        specification = {"version": "studio.comparison.v1", "evaluation_filters": {"evaluation_ids": frozen},
                         "settings": settings or {}, "selection_policy": "frozen evaluations"}
        group = self.create_group(user, name, project_id=project_id, origin={"type": "comparison", "collection_kind": "curated", "comparison_spec": specification})
        self.add_to_group(user, group, molecule_ids)
        return {"group_id": group, "molecule_count": len(rows), "comparison_spec": specification}

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
            metadata = row.get("metadata") or {}
            compound_id = re_safe(metadata.get("compound_id") or metadata.get("pubchem_cid")
                                  or metadata.get("catalog_id") or row["name"])[:80] or row["id"]
            if compound_id.casefold() in seen:
                compound_id += "-" + row["id"][-8:]
            seen.add(compound_id.casefold())
            compounds.append({"compound_id": compound_id, "compound_name": row["name"], "smiles": components[0]["smiles"],
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
