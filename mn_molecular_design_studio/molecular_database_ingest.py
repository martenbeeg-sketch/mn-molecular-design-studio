"""Idempotent adapters for Studio's existing molecular artifact contracts."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import logging
import socket
import time
from pathlib import Path

from mn_molecular_design_studio import molecule_imports as imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id, _json, _now

log = logging.getLogger(__name__)
_LAST_SYNC = {}
_LAST_SYNC_ERRORS = {}


def _read(path):
    try:
        value = json.loads(Path(path).read_text())
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _path(value, run_dir):
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    if not value:
        return None
    resolved = resolve_stored_path(value, run_dir=run_dir)
    path = Path(resolved) if resolved else Path(str(value))
    if not path.is_absolute():
        path = run_dir / path
    return path if path.is_file() else None


def _record(candidate, run_dir):
    files = []
    components = None
    # A candidate has one molecular identity. Do not attach target-only or
    # binder-only files as if they represented a different complex's identity.
    for key in ("complex_pdb", "binder_pdb", "target_pdb"):
        path = _path(candidate.get(key), run_dir)
        if not path:
            continue
        try:
            components = imports.parse_structure(path.read_bytes(), path.name)[0]["components"]
        except ValueError:
            continue
        files = [{"path": str(path), "label": key.replace("_", " ").title(), "dimensions": 3}]
        break
    sequence = candidate.get("binder_sequence") or candidate.get("sequence")
    raw_components = (candidate.get("raw_metadata") or {}).get("components")
    if components is None and isinstance(raw_components, list):
        try:
            components = imports.canonical_components(raw_components)[0]
        except (KeyError, TypeError, ValueError):
            components = None
    if components is None and sequence:
        components = imports.parse_sequences(str(sequence))[0]["components"]
    if not components:
        return None
    # An unsequenced poly-glycine backbone is a structure proposal, not an
    # assertion that every generated backbone is the same designed molecule.
    unresolved = str(candidate.get("stage") or "") == "generation.backbone" and not sequence
    identity_context = None
    if unresolved and files:
        identity_context = {"unsequenced_backbone": hashlib.sha256(Path(files[0]["path"]).read_bytes()).hexdigest()}
    return {"name": candidate.get("candidate_id") or "Molecule", "components": components, "files": files,
            "identity_context": identity_context,
            "metadata": {"stage": candidate.get("stage"), "status": (candidate.get("raw_metadata") or {}).get("status") or ("Unsequenced backbone" if unresolved else ""),
                         "source_tool": candidate.get("source_tool"), **(candidate.get("raw_metadata") or {})}}


def _candidate_rejected(candidate):
    raw = candidate.get("raw_metadata") or {}
    metrics = candidate.get("metrics") or {}
    status = str(candidate.get("status") or raw.get("status") or raw.get("filter_status") or metrics.get("status") or "").strip().casefold()
    if status in {"rejected", "failed", "invalid", "filtered", "discarded", "not accepted", "not_accepted"}:
        return True
    return candidate.get("accepted") is False or raw.get("accepted") is False or metrics.get("pass_filters") is False


def _run_mode(run_dir, metadata, inputs):
    if run_dir.parent.name == "candidate-import":
        return "curated"
    values = [run_dir.parent.name]
    for mapping in (metadata, inputs, inputs.get("params") or {}):
        values.extend(str(mapping.get(key) or "") for key in ("tool", "workflow", "workflow_id", "job_type", "engine", "task"))
    text = " ".join(values).casefold().replace("-", "_")
    if any(token in text for token in ("optimized_prediction", "refolding.optimized", "vina", "docking", "affinity", "rescoring", "evaluation", "boltz2", "boltz_2", "molecule_selection")):
        return "evaluation"
    if any(token in text for token in (
        "ligand.generation.", "design_generation", "bindcraft2", "rfantibody", "proteina_complexa_generator",
        "complexa_generator", "proteina_complexa", "boltzgen", "esmfold2_binder_design", "proteinmpnn",
        "ligandmpnn", "caliby", "esm_if1", "esmif1")):
        return "generated"
    return "auto"


def _canonical_fingerprint(record):
    _, fingerprint, _ = imports.canonical_components(record["components"])
    if record.get("identity_context"):
        fingerprint = hashlib.sha256((fingerprint + _json(record["identity_context"])).encode()).hexdigest()
    return fingerprint


def _explicit_molecule_ids(candidate):
    values = []
    if str(candidate.get("candidate_id") or "").startswith("mol-"):
        values.append(str(candidate["candidate_id"]))
    raw = candidate.get("raw_metadata") or {}
    for mapping in (candidate, raw):
        for key in ("studio_molecule_id", "molecule_id", "database_molecule_id", "source_molecule_id"):
            if mapping.get(key):
                values.append(str(mapping[key]))
    for parent in candidate.get("parents") or []:
        if isinstance(parent, dict):
            for key in ("studio_molecule_id", "molecule_id", "database_molecule_id"):
                if parent.get(key):
                    values.append(str(parent[key]))
    return list(dict.fromkeys(values))


def _find_original_group(database, owner, project, record, candidate=None):
    """Match only exact identities in a local group owned in this project."""
    fingerprint = _canonical_fingerprint(record)
    explicit_ids = _explicit_molecule_ids(candidate or {})
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT m.id,mm.group_id,g.origin,g.import_method,g.updated_at FROM molecules m "
            "JOIN molecule_members mm ON mm.molecule_id=m.id JOIN molecule_groups g ON g.id=mm.group_id "
            "WHERE m.fingerprint=? AND g.owner_id=? AND g.project_id=?",
            (fingerprint, str(owner), str(project or "")),
        ).fetchall()
    local = [dict(row) for row in rows if not json.loads(row["origin"]).get("_replica_node")]
    local.sort(key=lambda row: (
        row["id"] not in explicit_ids,
        MolecularDatabase._collection_kind({"origin": row["origin"]}) not in {"curated", "generated"},
        row["import_method"] not in {"user_uploads", "agent_inputs"},
        row["updated_at"],
        row["group_id"],
    ))
    return (local[0]["id"], local[0]["group_id"], True) if local else (None, None, False)


def _ensure_group(database, user, group_id, name, project, import_method, origin, collection_kind):
    from mn_molecular_design_studio.molecular_database import COLLECTION_KINDS
    if collection_kind not in COLLECTION_KINDS:
        raise ValueError("Unknown automatic collection kind.")
    with database.connection() as connection:
        row = connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group_id,)).fetchone()
        if row:
            if row["owner_id"] != user["id"] or row["project_id"] != str(project or ""):
                raise PermissionError("The automatic molecule group belongs to a different owner or project.")
            stored = json.loads(row["origin"])
            if stored.get("_replica_node") or stored.get("type") != "job":
                raise PermissionError("The automatic molecule group cannot replace a user collection.")
            stored["collection_kind"] = collection_kind
            connection.execute("UPDATE molecule_groups SET origin=? WHERE id=?", (_json(stored), group_id))
            return group_id
    return database.create_group(user, name, project_id=project, import_method=import_method,
                                 origin={**origin, "collection_kind": collection_kind}, group_id=group_id)


def _put_source(database, group_id, molecule_id, source_key, provenance):
    with database.connection() as connection:
        connection.execute("INSERT INTO molecule_sources VALUES (?,?,?,?) ON CONFLICT(source_key) DO UPDATE SET "
                           "molecule_id=excluded.molecule_id,group_id=excluded.group_id,provenance=excluded.provenance",
                           (source_key, molecule_id, group_id, _json(provenance)))


def _find_or_create_input(database, user, project, run_name, source_key, record, *, provenance, source_key_suffix, candidate=None):
    molecule_id, group_id, matched = _find_original_group(database, user["id"], project, record, candidate)
    if group_id:
        _put_source(database, group_id, molecule_id, source_key + source_key_suffix, provenance)
        return molecule_id, group_id, True
    group_id = _ensure_group(database, user, _id("inputs-", f"{user['id']}:{project or ''}"),
        "Tool inputs", project, "tool_inputs",
        {"type": "job", "owner_id": user["id"], "project_id": project}, "inputs")
    with database.connection() as connection:
        molecule_id = database._put(connection, group_id, record, source_key=source_key + source_key_suffix,
                                    origin=provenance)
    return molecule_id, group_id, False


def _is_completed_success(metadata, result):
    status = str(metadata.get("status") or "").strip().casefold()
    return status in {"completed", "complete", "success", "succeeded"} and result.get("success") is True


def _numeric(values):
    result = {}
    for key, value in values.items():
        if isinstance(value, (int, float, bool)):
            result[key] = value
        elif isinstance(value, str) and value.strip():
            try:
                result[key] = float(value)
            except ValueError:
                pass
    return result


def ingest_run(run_dir: Path, *, database=None, force=False):
    from mn_molecular_design_studio.workspace import is_studio_job
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    run_dir = Path(run_dir)
    if not is_studio_job(run_dir):
        return 0
    database = database or MolecularDatabase()
    marker = _read(run_dir / "molecular_design_studio.json")
    metadata, inputs, result = (_read(run_dir / name) for name in ("metadata.json", "input.json", "result.json"))
    embedded = metadata.get("molecular_design_studio") or {}
    owner = str(marker.get("studio_owner_id") or metadata.get("owner_id") or embedded.get("studio_owner_id") or "legacy-unassigned")
    project = str(marker.get("project_id") or metadata.get("project_id") or embedded.get("project_id") or "")
    user = {"id": owner, "role": "admin"}  # Trusted job adapter; all reads still enforce group access.
    source_key = f"{socket.gethostname()}/{run_dir.parent.name}/{run_dir.name}"
    observed_paths = [run_dir / name for name in (
        "input.json", "result.json", "metadata.json", "molecular_design_studio.json", "selection_snapshot.json",
        "artifacts/normalized_candidates/candidates.jsonl", "artifacts/input/compound_set.json",
        "artifacts/generated_compound_set.json", "artifacts/generated_compounds.sdf",
        "artifacts/docking_scores.csv", "artifacts/nesso_affinity.csv", "artifacts/input/receptor.pdb",
        "artifacts/input/target.pdb", "artifacts/input/complex.pdb")]
    signature = _json({"ingestion_version": 2, "files": [
        (str(path.relative_to(run_dir)), path.stat().st_mtime_ns, path.stat().st_size)
        for path in observed_paths if path.is_file()]})
    with database.connection() as connection:
        previous = connection.execute("SELECT signature FROM molecular_ingestions WHERE source_key=?", (source_key,)).fetchone()
    if not force and previous and previous[0] == signature:
        return 0
    if not _is_completed_success(metadata, result):
        with database.connection() as connection:
            connection.execute("INSERT INTO molecular_ingestions VALUES (?,?,?) ON CONFLICT(source_key) DO UPDATE SET signature=excluded.signature,updated_at=excluded.updated_at",
                               (source_key, signature, _now()))
        return 0

    tool = str(metadata.get("tool") or inputs.get("tool") or run_dir.parent.name)
    name = str(metadata.get("campaign_name") or (inputs.get("params") or {}).get("job_name") or metadata.get("job_name") or run_dir.name)
    stamp = str(metadata.get("completed_at") or metadata.get("created_at") or _now())
    origin = {"type": "job", "run_id": run_dir.name, "task_group": run_dir.parent.name, "run_dir": str(run_dir.resolve()),
              "tool": tool, "workstation": socket.gethostname(), "status": metadata.get("status")}
    mode = _run_mode(run_dir, metadata, inputs)
    output_method = "user_uploads" if mode == "curated" else "tool_output"
    seen_molecules = set()
    candidate_groups = {}

    def output_group(kind):
        if kind in {"generated", "curated"}:
            group_id = _id("job-", source_key)
        elif kind == "intermediate":
            group_id = _id("intermediate-", source_key)
        elif kind == "evaluation":
            group_id = _id("evaluation-", f"{user['id']}:{project or ''}")
        else:
            raise ValueError("Unsupported automatic result collection.")
        return _ensure_group(database, user, group_id, name if kind != "evaluation" else "Evaluated structures",
                             project, output_method, origin, kind)

    candidates = read_candidates(run_dir)
    candidate_records = []
    for candidate in candidates:
        record = _record(candidate, run_dir)
        if record is not None:
            candidate_records.append((candidate, record))

    input_set_path = run_dir / "artifacts/input/compound_set.json"
    generated_set_path = run_dir / "artifacts/generated_compound_set.json"
    score_paths = [run_dir / "artifacts" / filename for filename in ("docking_scores.csv", "nesso_affinity.csv")]
    md_metrics_supported = "md" in run_dir.parent.name and bool(result.get("metrics"))
    if not candidate_records and not generated_set_path.is_file() and not (input_set_path.is_file() and any(path.is_file() for path in score_paths)) and not md_metrics_supported:
        with database.connection() as connection:
            connection.execute("INSERT INTO molecular_ingestions VALUES (?,?,?) ON CONFLICT(source_key) DO UPDATE SET signature=excluded.signature,updated_at=excluded.updated_at",
                               (source_key, signature, _now()))
        return 0

    # Register generated and explicitly imported candidates in the library.
    # Backbone-only, rejected, and unmatched intermediate proposals stay in
    # the all-indexed view; evaluation outputs reuse an exact input identity.
    for candidate, record in candidate_records:
        stage = str(candidate.get("stage") or "")
        sequence = candidate.get("binder_sequence") or candidate.get("sequence")
        is_backbone = stage == "generation.backbone" and not sequence
        rejected = _candidate_rejected(candidate)
        provenance = {**origin, "candidate_id": candidate.get("candidate_id"), "parents": candidate.get("parents") or []}
        candidate_id = str(candidate.get("candidate_id") or _canonical_fingerprint(record))
        source_for_molecule = f"{source_key}:{candidate_id}"

        if mode == "curated":
            destination = output_group("curated")
            with database.connection() as connection:
                molecule_id = database._put(connection, destination, record, source_key=source_for_molecule, origin=provenance, stamp=stamp)
        elif not rejected and not is_backbone and (mode == "generated" or (mode == "auto" and stage in {"sequence_design", "generation.backbone_sequence"})):
            destination = output_group("generated")
            with database.connection() as connection:
                molecule_id = database._put(connection, destination, record, source_key=source_for_molecule, origin=provenance, stamp=stamp)
        elif not rejected and not is_backbone:
            molecule_id, destination, matched = _find_original_group(database, user["id"], project, record, candidate)
            if not matched:
                if mode == "evaluation" or stage in {"complex_refolding", "analysis", "benchmark", "monomer_refolding"}:
                    destination = output_group("evaluation")
                    with database.connection() as connection:
                        molecule_id = database._put(connection, destination, record, source_key=source_for_molecule,
                                                    origin=provenance, stamp=stamp)
                else:
                    destination = output_group("intermediate")
                    with database.connection() as connection:
                        molecule_id = database._put(connection, destination, record, source_key=source_for_molecule,
                                                    origin=provenance, stamp=stamp)
            else:
                _put_source(database, destination, molecule_id, source_for_molecule, provenance)
                for item in record.get("files") or []:
                    database.attach_structure(destination, molecule_id, item, context=provenance, created_at=stamp)
        else:
            destination = output_group("intermediate")
            with database.connection() as connection:
                molecule_id = database._put(connection, destination, record, source_key=source_for_molecule, origin=provenance, stamp=stamp)

        candidate_groups[candidate_id] = (molecule_id, destination)
        seen_molecules.add(molecule_id)
        metrics = _numeric(candidate.get("metrics") or {})
        if metrics:
            database.add_observation(destination, molecule_id, tool, run_dir.name, metrics,
                context=provenance, source_key=f"{source_for_molecule}:metrics", created_at=stamp)

    # The ligand contracts preserve exactly the consumed and generated sets.
    compound_rows = {}
    set_paths = [(input_set_path, False), (generated_set_path, True)]
    for path, generated in set_paths:
        if not path.is_file():
            continue
        from mn_molecular_design_studio.ligand_compounds import load_compound_set
        payload = load_compound_set(path)
        destination = output_group("generated") if generated else None
        for index, compound in enumerate(payload["compounds"]):
            record = {"name": compound["compound_id"], "components": [{"type": "small_molecule", "chain": "L", "smiles": compound["smiles"]}],
                      "metadata": {**(compound.get("properties") or {}), "source_tool": tool, "stage": "generated" if generated else "input"}}
            coordinates = run_dir / "artifacts/generated_compounds.sdf"
            if generated and coordinates.is_file():
                record["files"] = [{"path": str(coordinates), "label": "Generated coordinates", "molecule_index": index}]
            row_source_key = f"{source_key}:{'generated-' if generated else 'input-'}compound:{compound['compound_id']}"
            provenance = {**origin, "compound_id": compound["compound_id"], "collection_kind": "generated" if generated else "inputs"}
            if generated:
                with database.connection() as connection:
                    molecule_id = database._put(connection, destination, record, source_key=row_source_key, origin=provenance, stamp=stamp)
            else:
                molecule_id, destination, matched = _find_original_group(database, user["id"], project, record)
                if not matched:
                    destination = _ensure_group(database, user, _id("inputs-", f"{user['id']}:{project or ''}"),
                        "Tool inputs", project, "tool_inputs",
                        {"type": "job", "owner_id": user["id"], "project_id": project}, "inputs")
                    with database.connection() as connection:
                        molecule_id = database._put(connection, destination, record, source_key=row_source_key, origin=provenance, stamp=stamp)
                else:
                    _put_source(database, destination, molecule_id, row_source_key, provenance)
            compound_rows[compound["compound_id"]] = (molecule_id, destination, record)
            seen_molecules.add(molecule_id)
    # Preserve target context for docking/affinity and other structure-based tools.
    target_path = next((path for name in ("artifacts/input/receptor.pdb", "artifacts/input/target.pdb", "artifacts/input/complex.pdb") if (path := run_dir / name).is_file()), None)
    if target_path is None:
        input_values = inputs.get("inputs") or {}
        target_path = next((path for key in ("complex_pdb", "target_pdb", "protein_pdb", "target_structure") if (path := _path(input_values.get(key), run_dir)) is not None), None)
    target_id = None
    target_group_id = None
    target_record = None
    if target_path:
        try:
            target_record = imports.parse_structure(target_path.read_bytes(), target_path.name)[0]
            target_record["files"] = [{"path": str(target_path), "label": "Input receptor", "dimensions": 3}]
            target_id, target_group_id, matched = _find_original_group(database, user["id"], project, target_record)
            target_provenance = {**origin, "input_role": "receptor", "collection_kind": "inputs"}
            if not matched:
                target_group_id = _ensure_group(database, user, _id("inputs-", f"{user['id']}:{project or ''}"),
                    "Tool inputs", project, "tool_inputs",
                    {"type": "job", "owner_id": user["id"], "project_id": project}, "inputs")
                with database.connection() as connection:
                    target_id = database._put(connection, target_group_id, target_record, source_key=source_key + ":target", origin=target_provenance, stamp=stamp)
            else:
                _put_source(database, target_group_id, target_id, source_key + ":target", target_provenance)
                for item in target_record.get("files") or []:
                    database.attach_structure(target_group_id, target_id, item, context=target_provenance, created_at=stamp)
            seen_molecules.add(target_id)
        except ValueError as exc:
            log.info("Could not index receptor for %s: %s", run_dir, exc)
    params = inputs.get("params") or {}
    context = {**origin, "target_molecule_id": target_id, "target_path": str(target_path) if target_path else None,
               "settings": params, "units": (result.get("metrics") or {}).get("score_units"),
               "primary_ranking_metric": (result.get("metrics") or {}).get("primary_ranking_metric"),
               "ranking_direction": (result.get("metrics") or {}).get("ranking_direction")}
    for filename in ("docking_scores.csv", "nesso_affinity.csv"):
        path = run_dir / "artifacts" / filename
        if not path.is_file():
            continue
        text = path.read_text(errors="replace")
        reader = csv.DictReader(io.StringIO(text), delimiter="\t" if "\t" in text.partition("\n")[0] else ",")
        for index, row in enumerate(reader):
            match = compound_rows.get(str(row.get("compound_id") or ""))
            if not match:
                continue
            molecule_id, destination, record = match
            pose = _path(row.get("pose_path"), run_dir)
            if pose:
                database.attach_structure(destination, molecule_id, {"path": str(pose), "label": f"{tool} pose {row.get('pose_rank') or index+1}"},
                                          context={**context, "pose_rank": row.get("pose_rank")}, created_at=stamp)
            database.add_observation(destination, molecule_id, tool, run_dir.name, _numeric({key: value for key, value in row.items() if key != "compound_id"}),
                context={**context, "row": row}, source_key=f"{source_key}:{filename}:{index}", created_at=stamp)
    # MD produces observations and large artifacts for one already prepared complex.
    if "md" in run_dir.parent.name and target_id and target_group_id:
        database.add_observation(target_group_id, target_id, tool, run_dir.name, _numeric(result.get("metrics") or {}),
            context=context, source_key=f"{source_key}:target:metrics", created_at=stamp)
    count = len(seen_molecules)
    with database.connection() as connection:
        connection.execute("INSERT INTO molecular_ingestions VALUES (?,?,?) ON CONFLICT(source_key) DO UPDATE SET signature=excluded.signature,updated_at=excluded.updated_at", (source_key, signature, _now()))
    return count


def migrate_legacy_groups(database, user):
    from mn_molecular_design_studio.runtime import runs_root
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    path = database.root / "groups.json"
    try:
        groups = json.loads(path.read_text())
    except (OSError, ValueError):
        groups = []
    referenced_sets = {
        str(member.get("compound_set_id") or "")
        for group in groups if isinstance(group, dict)
        for member in (group.get("members") or []) if isinstance(member, dict)
        if member.get("record_type") == "compound"
    }
    for group in groups:
        if not isinstance(group, dict) or not group.get("id"):
            continue
        access = {"owner_id": group.get("owner_id") or "legacy-unassigned", "project_id": group.get("project_id") or ""}
        if not database._can_read(access, user):
            continue
        source_key = "legacy-group:" + str(group["id"])
        signature = hashlib.sha256(_json({"migration_version": 2, "group": group}).encode()).hexdigest()
        with database.connection() as connection:
            prior = connection.execute("SELECT signature FROM molecular_ingestions WHERE source_key=?", (source_key,)).fetchone()
        if prior and prior[0] == signature:
            continue
        actor = {"id": access["owner_id"], "role": "admin"}
        group_id = _id("legacy-", group["id"])
        records = []
        for member in group.get("members") or []:
            if member.get("record_type") == "compound":
                from mn_molecular_design_studio.ligand_compounds import load_compound_set
                set_id = str(member.get("compound_set_id") or "")
                if not set_id or any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for character in set_id):
                    continue
                set_path = database.root / "compound_sets" / (set_id + ".json")
                if not set_path.is_file():
                    continue
                payload = load_compound_set(set_path)
                # Access to a legacy group does not grant access to a referenced set.
                set_access = {"owner_id": payload.get("owner_id") or "legacy-unassigned", "project_id": payload.get("project_id") or ""}
                if not database._can_read(set_access, user):
                    continue
                compound = next((row for row in payload["compounds"] if row["compound_id"] == member.get("compound_id")), None)
                if compound:
                    records.append({"name": compound["compound_id"], "components": [{"type": "small_molecule", "chain": "L", "smiles": compound["smiles"]}], "metadata": compound.get("properties") or {}})
            else:
                task, run = str(member.get("task_group") or ""), str(member.get("run_id") or "")
                if not task or not run or Path(task).name != task or Path(run).name != run:
                    continue
                from mn_molecular_design_studio.job_catalog import get_job_access_record
                from mn_molecular_design_studio.auth import can_access_job
                if not can_access_job(user, get_job_access_record(task, run)):
                    continue
                run_dir = runs_root() / task / run
                candidate = next((row for row in read_candidates(run_dir) if row["candidate_id"] == member.get("candidate_id")), None)
                if candidate:
                    record = _record(candidate, run_dir)
                    if record:
                        records.append(record)
        if records:
            database.create_group(actor, group["name"], group_id=group_id, project_id=access["project_id"],
                                  import_method="user_uploads", origin={"type": "user", "legacy_group_id": group["id"]})
            database.import_records(actor, records, group_id=group_id)
        with database.connection() as connection:
            connection.execute("INSERT OR REPLACE INTO molecular_ingestions VALUES (?,?,?)", (source_key, signature, _now()))

    # Older compound imports lived in root-level immutable compound-set files.
    # Tool snapshots and pipeline selections are stored elsewhere and are not
    # promoted by this migration.
    sets_dir = database.root / "compound_sets"
    if not sets_dir.is_dir():
        return
    from mn_molecular_design_studio.ligand_compounds import load_compound_set
    for set_path in sorted(sets_dir.glob("*.json")):
        try:
            payload = load_compound_set(set_path)
        except (OSError, ValueError):
            continue
        set_id = str(payload.get("compound_set_id") or set_path.stem)
        if set_id in referenced_sets:
            continue
        access = {"owner_id": str(payload.get("owner_id") or "legacy-unassigned"),
                  "project_id": str(payload.get("project_id") or "")}
        if not database._can_read(access, user):
            continue
        source_key = "legacy-compound-set:" + set_id
        signature = hashlib.sha256(_json({"migration_version": 2, "payload": payload}).encode()).hexdigest()
        with database.connection() as connection:
            prior = connection.execute("SELECT signature FROM molecular_ingestions WHERE source_key=?", (source_key,)).fetchone()
        if prior and prior[0] == signature:
            continue
        actor = {"id": access["owner_id"], "role": "admin"}
        group_id = _id("legacy-cs-", set_id)
        database.create_group(actor, str(payload.get("name") or "Imported compounds"), group_id=group_id,
                              project_id=access["project_id"], import_method="user_uploads",
                              origin={"type": "legacy", "legacy_compound_set_id": set_id})
        records = [{"name": compound["compound_id"],
                    "components": [{"type": "small_molecule", "chain": "L", "smiles": compound["smiles"]}],
                    "metadata": compound.get("properties") or {}}
                   for compound in payload["compounds"]]
        database.import_records(actor, records, group_id=group_id)
        with database.connection() as connection:
            connection.execute("INSERT OR REPLACE INTO molecular_ingestions VALUES (?,?,?)", (source_key, signature, _now()))


def sync_workspace(user, *, force=False, database=None):
    from mn_molecular_design_studio.auth import accessible_project_ids
    from mn_molecular_design_studio.job_catalog import search_jobs
    database = database or MolecularDatabase()
    key = (str(database.path), str(user["id"]))
    if not force and time.monotonic() - _LAST_SYNC.get(key, 0) < 30:
        # A throttled refresh must not erase an indexing warning from the
        # previous scan while the underlying problem is still unresolved.
        return list(_LAST_SYNC_ERRORS.get(key, []))
    errors = []
    try:
        migrate_legacy_groups(database, user)
    except Exception as exc:
        log.exception("Molecular legacy migration failed for %s", user["id"])
        errors.append({"run_id": "legacy-migration", "error": str(exc)})
    offset = 0
    while True:
        page = search_jobs(access_owner_id=user["id"], accessible_project_ids=accessible_project_ids(user), access_admin=user.get("role") == "admin", limit=1000, offset=offset)
        for job in page["items"]:
            try:
                ingest_run(Path(job["run_dir"]), database=database, force=force)
            except Exception as exc:
                log.exception("Molecular indexing failed for %s", job["run_dir"])
                errors.append({"run_id": job["run_id"], "error": str(exc)})
        if not page["next_cursor"]:
            break
        offset = int(page["next_cursor"])
    _LAST_SYNC[key] = time.monotonic()
    _LAST_SYNC_ERRORS[key] = list(errors)
    return errors


def record_completed_run(run_dir):
    """Indexing cannot change the scientific job's completion status."""
    try:
        from mn_molecular_design_studio.runtime import runs_root
        if not Path(run_dir).resolve().is_relative_to(runs_root().resolve()):
            return
        ingest_run(Path(run_dir))
    except Exception:
        log.exception("Molecular Database could not index completed run %s; refresh will retry", run_dir)
