"""Molecular metadata replication and authorized retrieval of source files.

Only metadata enters peer synchronization. Files remain on the source node and
are fetched by an authorized user when previewed or reused.
"""
from __future__ import annotations
import hashlib
import json
import re
from contextlib import contextmanager
from pathlib import Path
from fastapi import Request

from mn_molecular_design_studio.molecular_database import IMPORT_METHODS, MolecularDatabase, _id, _json
from mn_molecular_design_studio import molecule_imports

KINDS = {"molecule_group", "molecule_entry", "external_library", "complex_review_annotation"}
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


def _review_structure_identity(connection, annotation, node_id):
    """Return a workstation-independent reference for the reviewed structure."""
    if annotation["source_kind"] == "database":
        row = connection.execute(
            "SELECT id,path,context FROM molecule_structures WHERE id=? AND molecule_id=?",
            (annotation["structure_ref"], annotation["subject_id"]),
        ).fetchone()
        if not row:
            raise ValueError("The reviewed Database structure is no longer present.")
        context = _object_json(row["context"], "reviewed structure context")
        source_node = str(context.get("_replica_node") or node_id)
        source_structure = str(context.get("_replica_structure") or row["id"])
        digest = str(context.get("sha256") or context.get("portable_sha256")
                     or context.get("content_sha256") or _digest(row["path"]) or "")
        model = context.get("model_id") or context.get("model") or context.get("model_index")
        return {
            "structure_identity": f"database:{source_node}:{source_structure}",
            "structure_digest": digest,
            "model_identity": str(model or ""),
            "structure_context": {
                key: context.get(key) for key in (
                    "candidate_id", "source_candidate_id", "model_id", "model",
                    "model_index", "target_chains", "binder_chains", "component_roles",
                ) if context.get(key) is not None
            },
        }
    # Results structures use job/run identity plus a safe relative artifact
    # path, not a local filesystem path or SQLite identifier.
    task_group, separator, run_id = str(annotation["subject_id"]).partition("/")
    relative = Path(str(annotation["structure_ref"]))
    if not separator or not task_group or not run_id or relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Invalid Results review identity.")
    digest = ""
    try:
        from mn_molecular_design_studio.job_catalog import get_job_access_record
        from mn_molecular_design_studio.runtime import runs_root
        job = get_job_access_record(task_group, run_id)
        if job:
            selected = (Path(runs_root()).expanduser().resolve() / str(job.get("run_relative_path") or f"{task_group}/{run_id}") / relative).resolve()
            if selected.is_file():
                digest = _digest(str(selected))
    except (OSError, RuntimeError, ValueError):
        pass
    return {
        "structure_identity": f"results:{annotation['subject_id']}:{relative.as_posix()}",
        "structure_digest": digest,
        "model_identity": "",
        "structure_context": {},
    }


def _review_record_key(annotation):
    identity = [
        annotation["source_kind"], annotation["subject_id"],
        annotation["structure_identity"], annotation.get("evaluation_id", ""),
        annotation["user_id"],
    ]
    return "complex_review_annotation:" + hashlib.sha256(_json(identity).encode()).hexdigest()


def _review_snapshot(database, row, node_id, connection=None):
    annotation = dict(row)
    if connection is None:
        with database.connection() as opened:
            structure = _review_structure_identity(opened, annotation, node_id)
    else:
        structure = _review_structure_identity(connection, annotation, node_id)
    return {
        **{key: annotation[key] for key in (
            "source_kind", "subject_id", "evaluation_id", "user_id", "status",
            "note", "created_at", "updated_at",
        )},
        **structure,
    }


def _canonical_components(value, label):
    try:
        return molecule_imports.canonical_components(value)
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {label} components.") from exc


def _source_components(value, label):
    """Validate a peer's stored components without rewriting small-molecule SMILES.

    Canonical SMILES are produced by RDKit, whose traversal and stereochemical
    serialization can change between releases. A molecule's persisted
    fingerprint is a hash of the source serialization, so receiving peers must
    validate that serialization and preserve it instead of canonicalizing it
    again with their local RDKit version.
    """
    if not isinstance(value, list) or not value:
        raise ValueError(f"Invalid {label} components.")

    normalized = []
    identities = []
    types = set()
    chains = set()
    for source_component in value:
        if not isinstance(source_component, dict):
            raise ValueError(f"Invalid {label} components.")
        kind = str(source_component.get("type") or "protein").lower()
        if kind not in molecule_imports.KINDS:
            raise ValueError(f"Invalid {label} components.")

        if kind == "small_molecule":
            component = {
                "type": kind,
                "chain": str(source_component.get("chain") or "A"),
                "role": str(source_component.get("role") or ""),
            }
            smiles = str(source_component.get("smiles") or "").strip()
            if smiles:
                try:
                    from rdkit import Chem
                    parsed = Chem.MolFromSmiles(smiles)
                except (ImportError, TypeError, ValueError) as exc:
                    raise ValueError(f"Invalid SMILES for chain {component['chain']}.") from exc
                if parsed is None:
                    raise ValueError(f"Invalid SMILES for chain {component['chain']}.")
                # Do not call MolToSmiles here: the exact source string is part
                # of this identity and may encode a different stereochemical
                # traversal under a newer RDKit release.
                component["smiles"] = smiles
            elif (source_component.get("ccd")
                  and re.fullmatch(r"[A-Za-z0-9_-]{1,20}", str(source_component["ccd"]))):
                component["ccd"] = str(source_component["ccd"]).upper()
            else:
                raise ValueError("A small molecule needs valid SMILES or a CCD component code.")
            if source_component.get("residue_key"):
                component["residue_key"] = str(source_component["residue_key"])
        else:
            try:
                component = molecule_imports.canonical_component(source_component)
            except (AttributeError, KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"Invalid {label} components.") from exc

        try:
            if _json(component) != _json(source_component):
                raise ValueError(f"Invalid {label} component serialization.")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid {label} component serialization.") from exc

        chain = component["chain"]
        if chain in chains:
            raise ValueError("Chain IDs must be unique within a molecule.")
        chains.add(chain)
        types.add(component["type"])
        identities.append((component["type"], component.get("sequence") or component.get("smiles")
                           or "ccd:" + component["ccd"]))
        normalized.append(component)

    fingerprint = hashlib.sha256(json.dumps(sorted(identities), separators=(",", ":")).encode()).hexdigest()
    molecule_type = next(iter(types)) if len(types) == 1 else "complex"
    return normalized, fingerprint, molecule_type


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


def _review_user(user_id):
    from mn_molecular_design_studio import auth
    user = next((item for item in auth.list_users()
                 if str(item.get("id") or "") == str(user_id) and not item.get("disabled")), None)
    if not user and auth.is_development_mode() and str(user_id) == str(auth.development_user()["id"]):
        user = auth.development_user()
    if not user:
        raise PermissionError("The shortlist reviewer is unknown or disabled on this workstation.")
    return user


def _local_review_structure(connection, annotation, node_id):
    for row in connection.execute(
        "SELECT id,path,context FROM molecule_structures WHERE molecule_id=?",
        (annotation["subject_id"],),
    ):
        context = _object_json(row["context"], "reviewed structure context")
        owner = str(context.get("_replica_node") or node_id)
        source_id = str(context.get("_replica_structure") or row["id"])
        if f"database:{owner}:{source_id}" != annotation["structure_identity"]:
            continue
        digest = str(context.get("sha256") or context.get("portable_sha256")
                     or context.get("content_sha256") or _digest(row["path"]) or "")
        if annotation.get("structure_digest") and digest != annotation["structure_digest"]:
            raise ValueError("The synchronized review refers to a different structure digest.")
        return str(row["id"])
    return None


def _apply_review_annotation(database, key, value, author, shared_connection=None):
    """Apply one replicated human decision after resolving local ACL-bound refs."""
    from mn_molecular_design_studio.cluster import config
    if value is None:
        # A null replication value is an explicit clear tombstone. Resolve its
        # stable key against current rows; never interpret a missing source as
        # permission to delete a different user's annotation.
        node_id = config.identity()
        with _record_transaction(database, shared_connection) as connection:
            rows = connection.execute(
                "SELECT source_kind,subject_id,structure_ref,evaluation_id,user_id,status,note,created_at,updated_at "
                "FROM complex_review_annotations"
            ).fetchall()
            for row in rows:
                try:
                    snapshot = _review_snapshot(database, row, node_id, connection)
                except (KeyError, TypeError, ValueError):
                    continue
                if _review_record_key(snapshot) == key:
                    connection.execute(
                        "DELETE FROM complex_review_annotations WHERE source_kind=? AND subject_id=? "
                        "AND structure_ref=? AND evaluation_id=? AND user_id=?",
                        (row["source_kind"], row["subject_id"], row["structure_ref"],
                         row["evaluation_id"], row["user_id"]),
                    )
        return
    required = {
        "source_kind", "subject_id", "structure_identity", "structure_digest", "model_identity",
        "structure_context", "evaluation_id", "user_id", "status", "note", "created_at", "updated_at",
    }
    if not isinstance(value, dict) or not required.issubset(value):
        raise ValueError("Incomplete synchronized shortlist annotation.")
    if (value.get("source_kind") not in {"database", "results"}
            or not isinstance(value.get("subject_id"), str) or not value["subject_id"]
            or not isinstance(value.get("structure_identity"), str) or not value["structure_identity"]
            or not isinstance(value.get("structure_digest"), str)
            or not isinstance(value.get("model_identity"), str)
            or not isinstance(value.get("structure_context"), dict)
            or not isinstance(value.get("evaluation_id"), str)
            or not isinstance(value.get("user_id"), str) or not value["user_id"]
            or value.get("status") not in {
                "Unreviewed", "Needs review", "Reviewed", "Shortlisted for experimental testing",
                "Not shortlisted for experimental testing", "Not applicable",
            }
            or not isinstance(value.get("note"), str) or len(value["note"]) > 10_000
            or not isinstance(value.get("created_at"), str) or not isinstance(value.get("updated_at"), str)):
        raise ValueError("Invalid synchronized shortlist annotation.")
    if _review_record_key(value) != key:
        raise PermissionError("The synchronized shortlist key does not match its scientific identity.")
    reviewer = _review_user(value["user_id"])
    node_id = config.identity()
    local_structure_ref = ""
    if value["source_kind"] == "database":
        with _record_transaction(database, shared_connection) as connection:
            local_structure_ref = _local_review_structure(connection, value, node_id) or ""
        if not local_structure_ref:
            raise ValueError("The reviewed structure has not synchronized to this workstation.")
        database._validate_complex_review_subject(
            reviewer, "database", value["subject_id"], local_structure_ref, value["evaluation_id"],
        )
    else:
        task_group, separator, run_id = value["subject_id"].partition("/")
        artifact_ref = str(value["structure_identity"]).removeprefix(f"results:{value['subject_id']}:")
        relative = Path(artifact_ref)
        if not separator or not task_group or not run_id or not artifact_ref or relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Invalid synchronized Results structure reference.")
        from mn_molecular_design_studio.auth import can_access_job
        from mn_molecular_design_studio.job_catalog import get_job_access_record
        job = get_job_access_record(task_group, run_id)
        if not job or not can_access_job(reviewer, job):
            raise PermissionError("The Results source is not authorized for this reviewer.")
        local_structure_ref = relative.as_posix()
    row_id = _id("review-", _json([
        value["source_kind"], value["subject_id"], local_structure_ref,
        value["evaluation_id"], value["user_id"],
    ]))
    with _record_transaction(database, shared_connection) as connection:
        connection.execute(
            "INSERT INTO complex_review_annotations "
            "(id,source_kind,subject_id,structure_ref,evaluation_id,user_id,status,note,created_at,updated_at) "
            "VALUES (?,?,?,?,?,?,?,?,?,?) ON CONFLICT(source_kind,subject_id,structure_ref,evaluation_id,user_id) "
            "DO UPDATE SET status=excluded.status,note=excluded.note,updated_at=excluded.updated_at",
            (row_id, value["source_kind"], value["subject_id"], local_structure_ref,
             value["evaluation_id"], value["user_id"], value["status"], value["note"],
             value["created_at"], value["updated_at"]),
        )


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
    source_components, molecule_composition, molecule_type = _source_components(components, "molecule")
    member_source_components, member_composition, member_type = _source_components(member_components, "member")
    if (_json(source_components) != molecule["components"]
            or _json(member_source_components) != member["components"]
            or molecule_composition != molecule["fingerprint"]
            or member_composition != molecule["fingerprint"]
            or molecule_type != member_type or molecule_composition != member_composition):
        raise ValueError("Molecular identity fingerprint does not match its source components.")
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


def _database_signature(database, replica):
    """Return a revision that changes only when this workstation edits its library."""
    with database.connection() as connection:
        row = connection.execute(
            "SELECT local_revision FROM molecular_database_sync_state WHERE singleton=1"
        ).fetchone()
    return {
        "path": str(database.path),
        "local_revision": int(row["local_revision"] if row else 0),
        "replica_epoch": replica.epoch,
    }


def capture(replica, *, database=None, node_id=None):
    from mn_molecular_design_studio.cluster import config
    database = database or MolecularDatabase()
    node_id = node_id or config.identity()
    cache_path = config.home() / "molecular-database-peer-capture.json"
    cache_key = hashlib.sha256(f"{node_id}\0{database.path}".encode()).hexdigest()
    try:
        cache = json.loads(cache_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        cache = {}
    signature = _database_signature(database, replica)
    previous = cache.get(cache_key) if isinstance(cache, dict) else None
    if isinstance(previous, dict) and previous.get("signature") == signature:
        return False

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
        for row in connection.execute(
            "SELECT source_kind,subject_id,structure_ref,evaluation_id,user_id,status,note,created_at,updated_at "
            "FROM complex_review_annotations"
        ):
            snapshot = _review_snapshot(database, row, node_id, connection)
            values[_review_record_key(snapshot)] = ("complex_review_annotation", snapshot)
    with replica.db() as connection:
        for key, (kind, value) in values.items():
            previous = connection.execute("SELECT value,author,counter FROM records WHERE key=?", (key,)).fetchone()
            previous_value = json.loads(previous["value"]) if previous else None
            if kind == "complex_review_annotation" and previous_value is not None:
                if previous_value.get("_conflicts"):
                    value = {**value, "_conflicts": previous_value["_conflicts"]}
                previous_content = {k: v for k, v in previous_value.items() if k != "_base_revision"}
                if previous_content == value:
                    continue
                value = {
                    **value,
                    "_base_revision": {"author": previous["author"], "counter": previous["counter"]},
                }
            elif previous_value == value:
                continue
            replica._write(connection, key, kind, value, author=node_id)

        # Stream the prior source snapshot in bounded pages. A stale replica may
        # have hundreds of thousands of rows, so holding its full JSON payloads
        # in a Python dict during a routine capture caused multi-gigabyte peaks.
        last_key = ""
        while True:
            previous_rows = connection.execute(
                "SELECT key,kind FROM records WHERE kind IN ('molecule_group','molecule_entry','external_library','complex_review_annotation') "
                "AND author=? AND key>? AND value!='null' ORDER BY key LIMIT 500",
                (node_id, last_key),
            ).fetchall()
            if not previous_rows:
                break
            for row in previous_rows:
                last_key = row["key"]
                if row["key"] not in values:
                    replica._write(connection, row["key"], row["kind"], None, author=node_id)

    # The source database can change while it is being scanned. Save a marker
    # only when the before/after file state matches; otherwise the next peer
    # heartbeat will reconcile it again.
    after = _database_signature(database, replica)
    if after == signature:
        cache[cache_key] = {"signature": signature}
        config.atomic_json(cache_path, cache)
    return True


@contextmanager
def _record_transaction(database, shared_connection):
    """Use one page connection while keeping each replicated row atomic."""
    if shared_connection is None:
        with database.connection() as connection:
            yield connection
        return
    with shared_connection:
        yield shared_connection


def _resolve_owner_qualified_structure(connection, context):
    """Translate a source-owned structure reference into this node's row ID."""
    owner_node = str(context.get("source_structure_owner_node") or context.get("source_node") or "")
    owner_group = str(context.get("source_structure_owner_group_id") or "")
    owner_structure = str(context.get("source_structure_owner_structure_id") or "")
    if not (owner_node and owner_group and owner_structure):
        return False
    group = connection.execute("SELECT id,origin FROM molecule_groups WHERE id=?", (owner_group,)).fetchone()
    local_group = None
    if group:
        origin = _object_json(group["origin"], "source structure group origin")
        if not origin.get("_replica_node"):
            local_group = group
    if local_group is None:
        local_group = connection.execute(
            "SELECT id,origin FROM molecule_groups WHERE json_extract(origin,'$._replica_node')=? "
            "AND json_extract(origin,'$._replica_group')=? LIMIT 1",
            (owner_node, owner_group),
        ).fetchone()
    if local_group is None:
        context["representation_id"] = ""
        context["source_structure_id"] = ""
        context["source_group_id"] = ""
        context["source_structure_reference_status"] = "owner_record_unavailable"
        return True
    structure = connection.execute(
        "SELECT id FROM molecule_structures WHERE id=? AND group_id=? AND molecule_id=?",
        (owner_structure, local_group["id"], context.get("source_molecule_id") or ""),
    ).fetchone()
    if structure is None:
        structure = connection.execute(
            "SELECT s.id FROM molecule_structures s WHERE s.group_id=? AND s.molecule_id=? "
            "AND json_extract(s.context,'$._replica_node')=? "
            "AND json_extract(s.context,'$._replica_structure')=? LIMIT 1",
            (local_group["id"], context.get("source_molecule_id") or "", owner_node, owner_structure),
        ).fetchone()
    if structure is None:
        context["representation_id"] = ""
        context["source_structure_id"] = ""
        context["source_group_id"] = str(local_group["id"])
        context["source_structure_reference_status"] = "owner_structure_unavailable"
        return True
    context["representation_id"] = str(structure["id"])
    context["source_structure_id"] = str(structure["id"])
    context["source_group_id"] = str(local_group["id"])
    context["source_structure_reference_status"] = "resolved"
    return True


def _reconcile_owner_qualified_observations(connection):
    """Resolve pending Evaluation links after a source group or entry arrives."""
    rows = connection.execute(
        "SELECT id,context FROM molecule_observations "
        "WHERE COALESCE(json_extract(context,'$.source_structure_owner_node'),'')!='' "
        "AND COALESCE(json_extract(context,'$.source_structure_owner_group_id'),'')!='' "
        "AND COALESCE(json_extract(context,'$.source_structure_owner_structure_id'),'')!='' "
        "AND COALESCE(json_extract(context,'$.source_structure_reference_status'),'')!='resolved'"
    ).fetchall()
    for row in rows:
        context = _object_json(row["context"], "observation context")
        previous = _json(context)
        if _resolve_owner_qualified_structure(connection, context) and _json(context) != previous:
            connection.execute("UPDATE molecule_observations SET context=? WHERE id=?", (_json(context), row["id"]))


def apply_record(kind, key, value, author, *, database=None, connection=None):
    database = database or MolecularDatabase()
    if not isinstance(kind, str) or kind not in KINDS or not isinstance(author, str) or not author:
        raise ValueError("Invalid molecular peer record authority.")
    if kind == "complex_review_annotation":
        # Reviews are mutable from any trusted peer. Their key is stable across
        # workstations and the shared Replica applies the established Lamport
        # ordering and reports concurrent edits as synchronization conflicts.
        if not isinstance(key, str) or not key.startswith(kind + ":") or len(key) > 300:
            raise ValueError("Invalid shortlist annotation key.")
        _apply_review_annotation(database, key, value, author, connection)
        return
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
    with _record_transaction(database, connection) as connection:
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
            _reconcile_owner_qualified_observations(connection)
            return
        if value is None:
            if _owned_group(connection, group_id, author, original_group):
                _remove_entry_relations(connection, group_id, molecule_id, author)
            return
        molecule, member, structures, observations, sources = _validate_entry(value, author, original_group, molecule_id)
        if not _owned_group(connection, group_id, author, original_group):
            raise ValueError("The source molecular group has not synchronized yet.")
        _, composition_fingerprint, molecule_type = _source_components(json.loads(molecule["components"]), "molecule")
        if molecule.get("type") != molecule_type:
            raise ValueError("Molecular type does not match its components.")
        existing_molecule = connection.execute("SELECT * FROM molecules WHERE id=?", (molecule_id,)).fetchone()
        if existing_molecule:
            try:
                _, existing_composition, existing_type = _source_components(json.loads(existing_molecule["components"]), "existing molecule")
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
            has_owner_qualified_reference = _resolve_owner_qualified_structure(connection, context)
            if context.get("representation_id") and not has_owner_qualified_reference:
                context["_replica_representation"] = context["representation_id"]
                context["representation_id"] = _id("peer-structure-", author + ":" + context["representation_id"])
            context.setdefault("source_node", author)
            context.update(_replica_node=author, _replica_observation=row["id"])
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
        _reconcile_owner_qualified_observations(connection)


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


def _removal_records(database, group_id, selected):
    """Publish only the authoritative group policy and selected tombstones."""
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.cluster.replication import Replica
    node_id = config.identity()
    with database.connection() as connection:
        group = dict(connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group_id,)).fetchone())
    if json.loads(group["origin"]).get("_replica_node"):
        raise PermissionError("A replica cannot publish source membership changes.")
    values = [("molecule_group:" + node_id + ":" + group_id, "molecule_group", {"node_id": node_id, "group": group})]
    values.extend((f"molecule_entry:{node_id}:{group_id}:{mid}", "molecule_entry", None) for mid in selected)
    replica, records = Replica(), []
    with replica.db() as connection:
        for key, kind, value in values:
            replica._write(connection, key, kind, value, author=node_id)
            counter = connection.execute("SELECT counter FROM records WHERE key=?", (key,)).fetchone()[0]
            records.append(dict(key=key, kind=kind, author=node_id, counter=counter, value=value))
    return records


def publish_member_removal(database, user, group_id, selected):
    """Persist a bounded mutation outbox even while bulk library sharing is paused."""
    from mn_molecular_design_studio.cluster import config
    # Explicit user removals are much smaller than a full library catch-up.
    # A source-authored tombstone still uses the normal Lamport/authority checks.
    directory = config.home() / "molecular-removals"
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    for start in range(0, len(selected), 499):
        batch = selected[start:start + 499]
        records.extend(_removal_records(database, group_id, batch))
        payload = {"user_id": user["id"], "source_group_id": group_id,
                   "molecule_ids": batch, "records": records[-(len(batch) + 1):]}
        key = hashlib.sha256(_json(payload).encode()).hexdigest()
        config.atomic_json(directory / (key + ".json"), {"payload": payload, "delivered": []})
    for peer in config.load().get("peers", []):
        flush_member_removals(peer)
    return records


def flush_member_removals(peer):
    """Retry explicit removal events without changing the bulk sharing policy."""
    from mn_molecular_design_studio.cluster import config, transport
    for path in sorted((config.home() / "molecular-removals").glob("*.json")):
        event = json.loads(path.read_text())
        if peer["id"] in event.get("delivered", []):
            continue
        try:
            response = transport.node_request("/molecules/removals", event["payload"], method="POST",
                                              base_url=peer["url"], timeout=15)
            if response.get("node_id") != peer["id"]:
                raise ValueError("The peer API belongs to a different workstation.")
        except (transport.ClusterUnavailable, ValueError, OSError):
            continue
        event["delivered"] = sorted(set(event.get("delivered", [])) | {peer["id"]})
        config.atomic_json(path, event)


def remove_remote_members(database, user, group, selected):
    """Forward an authorized request; never edit a replica optimistically."""
    from mn_molecular_design_studio.cluster import config, transport
    node_id = group["origin"]["_replica_node"]
    original_group = group["origin"]["_replica_group"]
    peer = next((p for p in config.load().get("peers", []) if p["id"] == node_id), None)
    if not peer:
        raise transport.ClusterUnavailable("Connect the source workstation to remove this Database membership.")
    if len(selected) > 499:
        return sum(remove_remote_members(database, user, group, selected[start:start + 499])
                   for start in range(0, len(selected), 499))
    response = transport.node_request("/molecules/groups/" + original_group + "/members/remove",
        {"user_id": user["id"], "molecule_ids": selected}, method="POST", base_url=peer["url"], timeout=60)
    if response.get("node_id") != node_id:
        raise ValueError("The mutation response belongs to a different workstation.")
    receive_member_removal(database, user, node_id, response["mutation"])
    return int(response["removed_count"])


def receive_member_removal(database, user, author, payload):
    """Validate the narrow source mutation before normal replication applies it."""
    from mn_molecular_design_studio.cluster.replication import Replica
    group_id = str(payload.get("source_group_id") or "")
    ids = payload.get("molecule_ids")
    records = payload.get("records")
    if (not group_id or not isinstance(ids, list) or not ids or len(ids) > 499
            or any(not isinstance(mid, str) or not mid.startswith("mol-") for mid in ids)
            or not isinstance(records, list) or len(records) != len(set(ids)) + 1):
        raise ValueError("Invalid source membership-removal event.")
    replica_group = _id("peer-", author + ":" + group_id)
    with database.connection() as connection:
        existing = connection.execute("SELECT id FROM molecule_groups WHERE id=?", (replica_group,)).fetchone()
    if existing:
        group = database.group(replica_group, user)
        if user.get("role") != "admin" and group["owner_id"] != user["id"]:
            raise PermissionError("Only the group owner or an administrator can remove membership.")
    expected = {f"molecule_entry:{author}:{group_id}:{mid}" for mid in ids}
    group_key = f"molecule_group:{author}:{group_id}"
    if {record.get("key") for record in records} != expected | {group_key}:
        raise ValueError("Removal records do not match the selected source group and molecules.")
    for record in records:
        if record.get("author") != author:
            raise PermissionError("Only the authoritative source may publish a removal.")
        if record["key"] in expected and (record.get("kind") != "molecule_entry" or record.get("value") is not None):
            raise ValueError("Membership removals require source entry tombstones.")
        if record["key"] == group_key:
            source_group = (record.get("value") or {}).get("group") or {}
            if user.get("role") != "admin" and source_group.get("owner_id") != user["id"]:
                raise PermissionError("Only the group owner or an administrator can remove membership.")
    return Replica().receive(records)


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

    def mutation_user(request, payload):
        authorization = request.headers.get("authorization", "")
        if not authorization.startswith("Bearer "):
            raise PermissionError("A workstation credential is required.")
        actor = Store().authenticate_node(authorization[7:])
        user = next((row for row in auth.list_users()
                     if row["id"] == payload.get("user_id") and not row.get("disabled")), None)
        if not user and auth.is_development_mode() and payload.get("user_id") == auth.development_user()["id"]:
            user = auth.development_user()
        if not user:
            raise PermissionError("Unknown or disabled user.")
        return actor, user

    @router.post("/groups/{group_id}/members/remove")
    def remove_members(group_id: str, request: Request, payload: dict = Body(...)):
        try:
            _, user = mutation_user(request, payload)
            selected = payload.get("molecule_ids")
            if (not isinstance(selected, list) or not 1 <= len(selected) <= 499
                    or any(not isinstance(mid, str) or not mid.startswith("mol-") for mid in selected)):
                raise ValueError("Select between 1 and 499 molecular identities.")
            selected = list(dict.fromkeys(selected))
            database = MolecularDatabase()
            # The write check rejects replicas and rechecks the actual source ACL.
            group = database.group(group_id, user, write=True)
            with database.connection() as connection:
                present = {row[0] for row in connection.execute(
                    "SELECT molecule_id FROM molecule_members WHERE group_id=? AND molecule_id IN (SELECT value FROM json_each(?))",
                    (group_id, _json(selected)))}
            removed = set(group["origin"].get("removed_molecule_ids") or [])
            if set(selected) - present - removed:
                raise PermissionError("Some selected molecules are not in the source group.")
            if present:
                database.remove_from_group(user, group_id, sorted(present))
            from mn_molecular_design_studio.cluster import config
            mutation = {"user_id": user["id"], "source_group_id": group_id,
                        "molecule_ids": selected, "records": _removal_records(database, group_id, selected)}
            return {"node_id": config.identity(), "removed_count": len(selected), "mutation": mutation}
        except PermissionError as exc:
            raise HTTPException(403, str(exc)) from exc
        except (KeyError, TypeError, ValueError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/removals")
    def apply_removal(request: Request, payload: dict = Body(...)):
        try:
            actor, user = mutation_user(request, payload)
            result = receive_member_removal(MolecularDatabase(), user, actor["id"], payload)
            from mn_molecular_design_studio.cluster import config
            return {"node_id": config.identity(), "accepted": result}
        except PermissionError as exc:
            raise HTTPException(403, str(exc)) from exc
        except (KeyError, TypeError, ValueError) as exc:
            raise HTTPException(422, str(exc)) from exc

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
        except (ValueError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

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
