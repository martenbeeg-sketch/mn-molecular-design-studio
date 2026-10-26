"""Molecular metadata replication and authorized retrieval of source files.

Only metadata enters peer synchronization. Files remain on the source node and
are fetched by an authorized user when previewed or reused.
"""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path
from fastapi import Request

from mn_molecular_design_studio.molecular_database import IMPORT_METHODS, MolecularDatabase, _id, _json
from mn_molecular_design_studio import molecule_imports

KINDS = {"molecule_group", "molecule_entry", "external_library"}
_DIGESTS = {}
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _object_json(value, label):
    """Decode a database JSON object and reject malformed peer metadata."""
    try:
        decoded = json.loads(value) if isinstance(value, str) else value
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid {label} metadata.") from exc
    if not isinstance(decoded, dict):
        raise ValueError(f"Invalid {label} metadata.")
    return decoded


def _canonical_components(value, label):
    try:
        return molecule_imports.canonical_components(value)
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {label} components.") from exc


def _owned_group(connection, group_id, author, original_group):
    row = connection.execute("SELECT origin FROM molecule_groups WHERE id=?", (group_id,)).fetchone()
    if not row:
        return False
    origin = _object_json(row["origin"], "group origin")
    if origin.get("_replica_node") != author or origin.get("_replica_group") != original_group:
        raise PermissionError("A peer cannot modify a molecular group owned by another source.")
    return True


def _remove_entry_relations(connection, group_id, molecule_id, author=None, *, remove_member=True):
    """Remove a member's replicated relations, keeping unrelated peer data.

    Entry tombstones address a whole (source-group, molecule) snapshot, so their
    relations can be removed by that exact pair. Reconciliation uses the author
    marker to remove only records copied from this source when a local record
    happens to share the same pair.
    """
    for table in ("molecule_sources", "molecule_structures", "molecule_observations"):
        rows = connection.execute(f"SELECT rowid,* FROM {table} WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchall()
        for row in rows:
            if table == "molecule_structures":
                metadata = _object_json(row["context"], "structure context")
            elif table == "molecule_observations":
                metadata = _object_json(row["context"], "observation context")
            else:
                metadata = _object_json(row["provenance"], "source provenance")
            marked_author = metadata.get("_replica_node")
            # Accept the source_node marker written by the earlier adapter for
            # stale-row cleanup during the first full-snapshot reconciliation.
            if author is None or marked_author == author or (marked_author is None and metadata.get("source_node") == author):
                connection.execute(f"DELETE FROM {table} WHERE rowid=?", (row["rowid"],))
    if remove_member:
        connection.execute("DELETE FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id))


def _validate_entry(value, author, original_group, molecule_id):
    if not isinstance(value, dict) or value.get("node_id") != author or value.get("group_id") != original_group:
        raise ValueError("Invalid molecular entry identity.")
    molecule = value.get("molecule")
    member = value.get("member")
    if not isinstance(molecule, dict) or not isinstance(member, dict) or molecule.get("id") != molecule_id:
        raise ValueError("Invalid molecular entry identity.")
    required_molecule = {"id", "fingerprint", "type", "components", "created_at"}
    required_member = {"molecule_id", "name", "components", "metadata", "added_at"}
    if not required_molecule.issubset(molecule) or not required_member.issubset(member):
        raise ValueError("Incomplete molecular entry snapshot.")
    if member.get("molecule_id") != molecule_id:
        raise ValueError("Invalid molecular entry membership.")
    if not isinstance(molecule.get("fingerprint"), str) or not _SHA256.fullmatch(molecule["fingerprint"]) or molecule_id != "mol-" + molecule["fingerprint"][:32]:
        raise ValueError("Molecular identity does not match its fingerprint.")
    if any(not isinstance(molecule.get(field), str) or not molecule[field] for field in ("type", "created_at")):
        raise ValueError("Invalid molecular identity metadata.")
    if any(not isinstance(member.get(field), str) or not member[field] for field in ("name", "added_at")):
        raise ValueError("Invalid molecular membership metadata.")
    try:
        components = json.loads(molecule["components"])
        member_components = json.loads(member["components"])
        member_metadata = json.loads(member["metadata"])
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValueError("Invalid molecular entry JSON.") from exc
    if not isinstance(components, list) or not isinstance(member_components, list) or not isinstance(member_metadata, dict):
        raise ValueError("Invalid molecular entry JSON.")
    canonical, molecule_composition, molecule_type = _canonical_components(components, "molecule")
    member_canonical, member_composition, member_type = _canonical_components(member_components, "member")
    if (_json(canonical) != molecule["components"] or _json(member_canonical) != member["components"]
            or molecule_type != member_type or molecule_composition != member_composition):
        raise ValueError("Molecular components are not in canonical form.")
    try:
        _json(member_metadata)
    except (TypeError, ValueError) as exc:
        raise ValueError("Invalid molecular member properties.") from exc
    structures, observations, sources = (value.get(name) for name in ("structures", "observations", "sources"))
    if not all(isinstance(rows, list) for rows in (structures, observations, sources)):
        raise ValueError("A molecular entry must include complete structure, score, and source snapshots.")

    clean_structures = []
    seen_structure_ids = set()
    for row in structures:
        if not isinstance(row, dict) or not {"id", "label", "format", "context", "created_at"}.issubset(row):
            raise ValueError("Invalid molecular structure metadata.")
        if row.get("molecule_id") != molecule_id or row.get("group_id") != original_group:
            raise ValueError("Molecular structure belongs to a different source entry.")
        source_id = row.get("id")
        fmt = row.get("format")
        if not isinstance(source_id, str) or not source_id or source_id in seen_structure_ids:
            raise ValueError("Invalid molecular structure identity.")
        if not isinstance(fmt, str) or not re.fullmatch(r"[a-zA-Z0-9]{1,16}", fmt):
            raise ValueError("Invalid molecular structure format.")
        if not isinstance(row.get("label"), str) or not isinstance(row.get("created_at"), str):
            raise ValueError("Invalid molecular structure metadata.")
        seen_structure_ids.add(source_id)
        context = _object_json(row.get("context"), "structure context")
        digest = row.get("sha256") or ""
        if not isinstance(digest, str) or (digest and not _SHA256.fullmatch(digest)):
            raise ValueError("Invalid molecular structure checksum.")
        clean_structures.append((row, context, digest))

    clean_observations = []
    seen_observation_ids = set()
    for row in observations:
        if not isinstance(row, dict) or not {"id", "namespace", "run_id", "metrics", "context", "created_at"}.issubset(row):
            raise ValueError("Invalid molecular observation metadata.")
        if row.get("molecule_id") != molecule_id or row.get("group_id") != original_group:
            raise ValueError("Molecular observation belongs to a different source entry.")
        source_id = row.get("id")
        if not isinstance(source_id, str) or not source_id or source_id in seen_observation_ids:
            raise ValueError("Invalid molecular observation identity.")
        seen_observation_ids.add(source_id)
        try:
            metrics = json.loads(row.get("metrics"))
        except (TypeError, json.JSONDecodeError) as exc:
            raise ValueError("Invalid molecular observation metrics.") from exc
        if not isinstance(metrics, dict):
            raise ValueError("Invalid molecular observation metrics.")
        try:
            _json(metrics)
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid molecular observation metrics.") from exc
        if any(not isinstance(row.get(field), str) for field in ("namespace", "run_id", "created_at")):
            raise ValueError("Invalid molecular observation metadata.")
        context = _object_json(row.get("context"), "observation context")
        clean_observations.append((row, context))

    clean_sources = []
    seen_source_keys = set()
    for row in sources:
        if not isinstance(row, dict) or not {"source_key", "molecule_id", "group_id", "provenance"}.issubset(row):
            raise ValueError("Invalid molecular source metadata.")
        if row.get("molecule_id") != molecule_id or row.get("group_id") != original_group:
            raise ValueError("Molecular source belongs to a different source entry.")
        source_key = row.get("source_key")
        if not isinstance(source_key, str) or not source_key or source_key in seen_source_keys:
            raise ValueError("Invalid molecular source identity.")
        seen_source_keys.add(source_key)
        provenance = _object_json(row.get("provenance"), "source provenance")
        clean_sources.append((row, provenance))
    return molecule, member, clean_structures, clean_observations, clean_sources


def _digest(path):
    path = Path(path)
    if not path.is_file():
        return ""
    stat = path.stat()
    key = (str(path), stat.st_size, stat.st_mtime_ns)
    if key not in _DIGESTS:
        with path.open("rb") as handle:
            _DIGESTS[key] = hashlib.file_digest(handle, "sha256").hexdigest()
    return _DIGESTS[key]


def capture(replica, *, database=None, node_id=None):
    from mn_molecular_design_studio.cluster import config
    database = database or MolecularDatabase()
    node_id = node_id or config.identity()
    values = {}
    with database.connection() as connection:
        for row in connection.execute("SELECT id,record_json FROM external_library_records"):
            record = _object_json(row["record_json"], "external library")
            author = record.get("_replica_node")
            if author not in (None, node_id):
                continue
            values[f"external_library:{node_id}:{row['id']}"] = (
                "external_library", {"node_id": node_id, "record": record},
            )
        groups = connection.execute("SELECT * FROM molecule_groups").fetchall()
        local_groups = [row for row in groups if not json.loads(row["origin"]).get("_replica_node")]
        for group in local_groups:
            values["molecule_group:" + node_id + ":" + group["id"]] = ("molecule_group", {"node_id": node_id, "group": dict(group)})
            members = connection.execute("SELECT * FROM molecule_members WHERE group_id=?", (group["id"],)).fetchall()
            for member in members:
                molecule = connection.execute("SELECT * FROM molecules WHERE id=?", (member["molecule_id"],)).fetchone()
                structures = []
                for row in connection.execute("SELECT * FROM molecule_structures WHERE group_id=? AND molecule_id=?", (group["id"], member["molecule_id"])):
                    context = json.loads(row["context"])
                    target_path = context.pop("target_path", None)
                    if target_path and not context.get("target_structure_id"):
                        target = connection.execute("SELECT id FROM molecule_structures WHERE path=? LIMIT 1", (target_path,)).fetchone()
                        if target:
                            context["target_structure_id"] = target[0]
                    structures.append({**dict(row), "path": "", "context": _json(context), "sha256": _digest(row["path"])})
                observations = [dict(row) for row in connection.execute("SELECT * FROM molecule_observations WHERE group_id=? AND molecule_id=?", (group["id"], member["molecule_id"]))]
                sources = [dict(row) for row in connection.execute("SELECT * FROM molecule_sources WHERE group_id=? AND molecule_id=?", (group["id"], member["molecule_id"]))]
                value = {"node_id": node_id, "group_id": group["id"], "molecule": dict(molecule), "member": dict(member), "structures": structures, "observations": observations, "sources": sources}
                values[f"molecule_entry:{node_id}:{group['id']}:{member['molecule_id']}"] = ("molecule_entry", value)
    with replica.db() as connection:
        previous = {row["key"]: row for row in connection.execute("SELECT * FROM records WHERE kind IN ('molecule_group','molecule_entry','external_library') AND author=?", (node_id,))}
        for key, (kind, value) in values.items():
            if key not in previous or json.loads(previous[key]["value"]) != value:
                replica._write(connection, key, kind, value, author=node_id)
        for key, row in previous.items():
            if key not in values and json.loads(row["value"]) is not None:
                replica._write(connection, key, row["kind"], None, author=node_id)


def apply_record(kind, key, value, author, *, database=None):
    database = database or MolecularDatabase()
    if not isinstance(kind, str) or kind not in KINDS or not isinstance(author, str) or not author:
        raise ValueError("Invalid molecular peer record authority.")
    prefix = kind + ":" + author + ":"
    if not isinstance(key, str) or not key.startswith(prefix):
        raise PermissionError("Only the source workstation may publish a molecular record.")
    identifier = key[len(prefix):]
    if kind == "external_library":
        if not identifier.startswith("lib-") or len(identifier) != 28 or any(c not in "0123456789abcdef" for c in identifier[4:]):
            raise ValueError("Invalid external library identity.")
        if value is None:
            with database.connection() as connection:
                row = connection.execute("SELECT record_json FROM external_library_records WHERE id=?", (identifier,)).fetchone()
                if row:
                    record = _object_json(row["record_json"], "external library")
                    if record.get("_replica_node") != author:
                        raise PermissionError("A peer cannot remove an external library owned by another workstation.")
                    connection.execute("DELETE FROM external_library_records WHERE id=?", (identifier,))
            return
        if (not isinstance(value, dict) or value.get("node_id") != author
                or not isinstance(value.get("record"), dict) or value["record"].get("id") != identifier):
            raise ValueError("Invalid external library snapshot.")
        record = dict(value["record"])
        if record.get("schema_version") != "studio.external-library.v1" or not isinstance(record.get("sources"), list):
            raise ValueError("Incomplete external library metadata.")
        if len(json.dumps(record, ensure_ascii=False, allow_nan=False, default=str)) > 1_000_000:
            raise ValueError("External library metadata exceeds the synchronization limit.")
        with database.connection() as connection:
            existing = connection.execute("SELECT record_json FROM external_library_records WHERE id=?", (identifier,)).fetchone()
            if existing:
                previous = _object_json(existing["record_json"], "external library")
                if previous.get("_replica_node") not in (None, author):
                    raise PermissionError("A peer cannot replace an external library owned by another workstation.")
                from mn_molecular_design_studio.cluster import config
                local_preparations = [item for item in previous.get("preparations", [])
                                      if item.get("owner_node_id") == config.identity()]
                known = {item.get("sha256") for item in record.get("preparations", [])}
                record.setdefault("preparations", []).extend(
                    item for item in local_preparations if item.get("sha256") not in known
                )
                if local_preparations:
                    latest = local_preparations[-1]
                    if latest.get("schema_version") == "studio.library-fingerprints.v2":
                        record["prepared_compound_count"] = latest.get("rows")
                        record["rejected_row_count"] = latest.get("rejected_rows")
                        record["fingerprint_spec"] = latest.get("fingerprint_spec")
                    elif latest.get("schema_version") == "studio.library-centroids.v2":
                        record["centroid_count"] = latest.get("rows")
        record["_replica_node"] = author
        database.save_external_library_record(record)
        return
    if not identifier or (kind == "molecule_entry" and ":" not in identifier):
        raise ValueError("Invalid molecular peer record key.")
    original_group = identifier if kind == "molecule_group" else identifier.rsplit(":", 1)[0]
    molecule_id = "" if kind == "molecule_group" else identifier.rsplit(":", 1)[-1]
    if not original_group or (kind == "molecule_entry" and not molecule_id):
        raise ValueError("Invalid molecular peer record key.")
    group_id = _id("peer-", author + ":" + original_group)
    with database.connection() as connection:
        if kind == "molecule_group":
            if value is None:
                if not _owned_group(connection, group_id, author, original_group):
                    return
                for table in ("molecule_sources", "molecule_structures", "molecule_observations", "molecule_members"):
                    connection.execute(f"DELETE FROM {table} WHERE group_id=?", (group_id,))
                connection.execute("DELETE FROM molecule_groups WHERE id=?", (group_id,))
                return
            if not isinstance(value, dict) or value.get("node_id") != author or not isinstance(value.get("group"), dict):
                raise ValueError("Invalid molecular group metadata.")
            group = value["group"]
            if (group.get("id") != original_group or not isinstance(group.get("owner_id"), str) or not group["owner_id"]
                    or not isinstance(group.get("name"), str) or not group["name"].strip()):
                raise ValueError("Invalid molecular group identity.")
            project_id = group.get("project_id") or ""
            import_method = group.get("import_method")
            if not isinstance(project_id, str) or not isinstance(import_method, str) or import_method not in IMPORT_METHODS:
                raise ValueError("Invalid molecular group access or collection type.")
            raw_origin = _object_json(group.get("origin"), "group origin")
            schema = _object_json(group.get("schema_json"), "group schema")
            try:
                from mn_molecular_design_studio.molecule_schema import validate_schema
                validate_schema(schema)
            except (AttributeError, TypeError, ValueError) as exc:
                raise ValueError("Invalid molecular group schema.") from exc
            origin = {**raw_origin, "_replica_node": author, "_replica_group": original_group}
            exists = _owned_group(connection, group_id, author, original_group)
            if exists:
                existing_owner = connection.execute("SELECT owner_id FROM molecule_groups WHERE id=?", (group_id,)).fetchone()[0]
                if str(existing_owner) != str(group["owner_id"]):
                    raise PermissionError("Peer synchronization cannot transfer a molecular group to another owner.")
            connection.execute("INSERT INTO molecule_groups VALUES (?,?,?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET name=excluded.name,project_id=excluded.project_id,import_method=excluded.import_method,origin=excluded.origin,schema_json=excluded.schema_json,updated_at=excluded.updated_at",
                (group_id, str(group["name"]).strip()[:250], str(group["owner_id"]), project_id, import_method, _json(origin), _json(schema), group["created_at"], group["updated_at"]))
            return
        if value is None:
            if _owned_group(connection, group_id, author, original_group):
                _remove_entry_relations(connection, group_id, molecule_id, author)
            return
        molecule, member, structures, observations, sources = _validate_entry(value, author, original_group, molecule_id)
        if not _owned_group(connection, group_id, author, original_group):
            raise ValueError("The source molecular group has not synchronized yet.")
        _, composition_fingerprint, molecule_type = _canonical_components(json.loads(molecule["components"]), "molecule")
        if molecule.get("type") != molecule_type:
            raise ValueError("Molecular type does not match its components.")
        existing_molecule = connection.execute("SELECT * FROM molecules WHERE id=?", (molecule_id,)).fetchone()
        if existing_molecule:
            try:
                _, existing_composition, existing_type = _canonical_components(json.loads(existing_molecule["components"]), "existing molecule")
            except (TypeError, json.JSONDecodeError, ValueError) as exc:
                raise ValueError("An existing molecular identity has invalid components.") from exc
            if (existing_molecule["fingerprint"] != molecule["fingerprint"] or
                    existing_molecule["type"] != molecule["type"] or
                    existing_type != molecule_type or existing_composition != composition_fingerprint):
                raise ValueError("A synchronized molecule identity conflicts with an existing record.")
        for table, id_column, relation_rows, metadata_column in (
            ("molecule_structures", "id", [(row, context) for row, context, _ in structures], "context"),
            ("molecule_observations", "id", observations, "context"),
            ("molecule_sources", "source_key", sources, "provenance"),
        ):
            for row, metadata in relation_rows:
                if table == "molecule_structures":
                    relation_id = _id("peer-structure-", author + ":" + row["id"])
                elif table == "molecule_observations":
                    relation_id = _id("peer-score-", author + ":" + row["id"])
                else:
                    relation_id = _id("peer-source-", author + ":" + row["source_key"])
                existing = connection.execute(
                    f"SELECT group_id,molecule_id,{metadata_column} FROM {table} WHERE {id_column}=?", (relation_id,)
                ).fetchone()
                if existing:
                    existing_metadata = _object_json(existing[metadata_column], "replicated relation metadata")
                    source_owned = existing_metadata.get("_replica_node") == author or (
                        existing_metadata.get("_replica_node") is None and existing_metadata.get("source_node") == author
                    )
                    if existing["group_id"] != group_id or existing["molecule_id"] != molecule_id or not source_owned:
                        raise PermissionError("A synchronized molecular relation conflicts with a local or other-peer record.")
        connection.execute("INSERT OR IGNORE INTO molecules VALUES (?,?,?,?,?)", (molecule_id, molecule["fingerprint"], molecule["type"], molecule["components"], molecule["created_at"]))
        connection.execute("INSERT OR REPLACE INTO molecule_members VALUES (?,?,?,?,?,?)", (group_id, molecule_id, member["name"], member["components"], member["metadata"], member["added_at"]))
        current_structures = set()
        for row, context, digest in structures:
            current_structures.add(_id("peer-structure-", author + ":" + row["id"]))
            context.update(_replica_node=author, _replica_structure=row["id"], sha256=digest)
            if context.get("target_structure_id"):
                context["target_structure_id"] = _id("peer-structure-", author + ":" + context["target_structure_id"])
            structure_id = _id("peer-structure-", author + ":" + row["id"])
            structure_path = "cluster-library://" + author + "/" + row["id"]
            existing_structure = connection.execute("SELECT path,context FROM molecule_structures WHERE id=?", (structure_id,)).fetchone()
            if context.get("portable") is True and existing_structure and not str(existing_structure["path"]).startswith("cluster-library://"):
                existing_context = _object_json(existing_structure["context"], "existing structure context")
                existing_path = Path(existing_structure["path"])
                if (existing_context.get("portable") is True and existing_path.is_file()
                        and _digest(existing_path) == digest):
                    structure_path = str(existing_path)
            connection.execute("INSERT OR REPLACE INTO molecule_structures VALUES (?,?,?,?,?,?,?,?)", (structure_id, molecule_id, group_id, structure_path, row["label"], row["format"], _json(context), row["created_at"]))
        current_observations = set()
        for row, context in observations:
            current_observations.add(_id("peer-score-", author + ":" + row["id"]))
            context.pop("target_path", None)
            if context.get("representation_id"):
                context["_replica_representation"] = context["representation_id"]
                context["representation_id"] = _id("peer-structure-", author + ":" + context["representation_id"])
            context.update(source_node=author, _replica_node=author, _replica_observation=row["id"])
            connection.execute("INSERT OR REPLACE INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)", (_id("peer-score-", author + ":" + row["id"]), molecule_id, group_id, row["namespace"], row["run_id"], row["metrics"], _json(context), row["created_at"]))
        current_sources = set()
        for row, provenance in sources:
            current_sources.add(_id("peer-source-", author + ":" + row["source_key"]))
            provenance.update(source_node=author, _replica_node=author, _replica_source_key=row["source_key"])
            connection.execute("INSERT OR REPLACE INTO molecule_sources VALUES (?,?,?,?)", (_id("peer-source-", author + ":" + row["source_key"]), molecule_id, group_id, _json(provenance)))

        # A record is a full snapshot for this source's one group/molecule. Drop
        # stale copied relations while leaving local rows and other-source rows
        # untouched. Existing records from the first adapter used source_node
        # as their ownership marker, so also reconcile those during migration.
        for table, id_column, expected, metadata_column in (
            ("molecule_structures", "id", current_structures, "context"),
            ("molecule_observations", "id", current_observations, "context"),
            ("molecule_sources", "source_key", current_sources, "provenance"),
        ):
            rows = connection.execute(f"SELECT rowid,{id_column},{metadata_column} FROM {table} WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchall()
            for row in rows:
                metadata = _object_json(row[metadata_column], "replicated relation metadata")
                from_author = metadata.get("_replica_node") == author or (metadata.get("_replica_node") is None and metadata.get("source_node") == author)
                if from_author and row[id_column] not in expected:
                    connection.execute(f"DELETE FROM {table} WHERE rowid=?", (row["rowid"],))


def fetch_structure(row, user):
    from mn_molecular_design_studio.cluster import config, transport
    from mn_molecular_design_studio.auth import user_input_root
    context = json.loads(row["context"])
    node_id, structure_id = context["_replica_node"], context["_replica_structure"]
    peer = next((peer for peer in config.load().get("peers", []) if peer["id"] == node_id), None)
    if not peer:
        raise FileNotFoundError("Connect the source workstation to retrieve this structure.")
    digest = context.get("sha256")
    path = user_input_root(user["id"]) / "database_remote" / ((digest or _id("file-", row["id"])) + "." + row["format"])
    if path.is_file() and (not digest or _digest(path) == digest):
        return path
    try:
        content = transport.node_request("/molecules/structures/" + structure_id, {"user_id": user["id"]}, method="POST", base_url=peer["url"], raw=True, timeout=60)
    except transport.ClusterUnavailable as exc:
        raise FileNotFoundError("The source workstation is unavailable. Its molecular metadata remains visible.") from exc
    if digest and hashlib.sha256(content).hexdigest() != digest:
        raise ValueError("The source structure changed; refresh the molecular records before using it.")
    path.parent.mkdir(parents=True, exist_ok=True)
    import tempfile
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        handle.write(content)
        temporary = Path(handle.name)
    temporary.replace(path)
    return path


def materialize_portable_structures(*, database=None, peer=None):
    """Cache portable peer structures in this workstation's shared Database store."""
    from mn_molecular_design_studio.cluster import config, transport

    database = database or MolecularDatabase()
    peers = [peer] if peer else [item for item in config.load().get("peers", []) if item.get("id") and item.get("url")]
    peers = [item for item in peers if item and item.get("id") and item.get("url")]
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT s.*,g.owner_id FROM molecule_structures s "
            "JOIN molecule_groups g ON g.id=s.group_id WHERE s.path LIKE 'cluster-library://%'"
        ).fetchall()
    materialized = 0
    for row in rows:
        context = _object_json(row["context"], "structure context")
        if context.get("portable") is not True:
            continue
        source_node = str(context.get("_replica_node") or "")
        source_structure = str(context.get("_replica_structure") or "")
        digest = str(context.get("portable_sha256") or context.get("sha256") or "")
        if not source_node or not source_structure or not _SHA256.fullmatch(digest):
            continue
        local_id = _id("peer-structure-", source_node + ":" + source_structure)
        content = None
        for remote_peer in peers:
            remote_structure_id = source_structure if remote_peer["id"] == source_node else local_id
            try:
                content = transport.node_request(
                    "/molecules/structures/" + remote_structure_id,
                    {"user_id": row["owner_id"]}, method="POST", base_url=remote_peer["url"], raw=True, timeout=60,
                )
            except (transport.ClusterUnavailable, OSError, ValueError):
                continue
            if hashlib.sha256(content).hexdigest() == digest:
                break
            content = None
        if content is None:
            continue
        destination = database._store_file("portable." + str(row["format"] or "pdb"), content)
        with database.connection() as connection:
            current = connection.execute("SELECT context,path FROM molecule_structures WHERE id=?", (row["id"],)).fetchone()
            if not current:
                continue
            current_context = _object_json(current["context"], "structure context")
            if (current_context.get("portable") is True
                    and current_context.get("portable_sha256", current_context.get("sha256")) == digest):
                connection.execute("UPDATE molecule_structures SET path=? WHERE id=?", (str(destination), row["id"]))
                materialized += 1
    return materialized


def install_routes(app):
    from fastapi import APIRouter, Body, HTTPException, Request
    from fastapi.responses import FileResponse
    from mn_molecular_design_studio.cluster.store import Store
    from mn_molecular_design_studio import auth
    router = APIRouter(prefix="/api/v1/cluster/node/molecules", tags=["molecular database peers"])

    @router.post("/structures/{structure_id}")
    def file(structure_id: str, request: Request, payload: dict = Body(...)):
        authorization = request.headers.get("authorization", "")
        try:
            if not authorization.startswith("Bearer "):
                raise PermissionError("A workstation credential is required.")
            Store().authenticate_node(authorization[7:])
            user = next((row for row in auth.list_users() if row["id"] == payload.get("user_id") and not row.get("disabled")), None)
            if not user and auth.is_development_mode() and payload.get("user_id") == auth.development_user()["id"]:
                user = auth.development_user()
            if not user:
                raise PermissionError("Unknown or disabled user.")
            database = MolecularDatabase()
            path = database.structure_file(user, structure_id)
            return FileResponse(path, filename=path.name)
        except (PermissionError, FileNotFoundError) as exc:
            raise HTTPException(404, "Structure not found or inaccessible.") from exc

    @router.post("/structures/{structure_id}/portable")
    def mark_portable(structure_id: str, request: Request, payload: dict = Body(...)):
        authorization = request.headers.get("authorization", "")
        try:
            if not authorization.startswith("Bearer "):
                raise PermissionError("A workstation credential is required.")
            Store().authenticate_node(authorization[7:])
            user = next((row for row in auth.list_users() if row["id"] == payload.get("user_id") and not row.get("disabled")), None)
            if not user and auth.is_development_mode() and payload.get("user_id") == auth.development_user()["id"]:
                user = auth.development_user()
            if not user:
                raise PermissionError("Unknown or disabled user.")
            result = MolecularDatabase().make_structure_portable(user, structure_id)
            from mn_molecular_design_studio.cluster.replication import Replica
            capture(Replica())
            return result
        except (PermissionError, FileNotFoundError, ValueError) as exc:
            raise HTTPException(404, "Structure not found or inaccessible.") from exc
    app.include_router(router)
