"""Reusable source proteins and structures, independent of prediction outcomes.

Targets remain roles in run provenance. Sequence-derived PDB references remain
separate from supplied templates and from generated/predicted representations.
"""
from __future__ import annotations

import hashlib
import json
import logging
import re
import socket
import sqlite3
import tempfile
import time
from collections import Counter
from pathlib import Path

from mn_molecular_design_studio import molecule_imports as imports
from mn_molecular_design_studio.target_identity import is_sequence_carrier_reference, protein_sequence_reference

log = logging.getLogger(__name__)
_FETCH_FAILURES: dict[str, float] = {}
_FILE_FIELDS = ("target_structure", "target_pdb", "protein_pdb", "complex_pdb", "source_pdb",
                "structure_path", "input_pdb", "pdb_path")


class _NotProteinInput(ValueError):
    pass


def _json_file(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def input_structure_paths(input_data: dict, run_dir: Path) -> list[tuple[str, Path]]:
    """Resolve declared structure inputs, without scanning output artifacts."""
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

    found = {}
    for mapping in (input_data, input_data.get("inputs") or {}, input_data.get("params") or {}):
        if not isinstance(mapping, dict):
            continue
        for key in _FILE_FIELDS:
            value = mapping.get(key)
            if not isinstance(value, str) or not value or is_sequence_carrier_reference(value, run_dir=run_dir):
                continue
            path = resolve_stored_path(value, run_dir=run_dir, must_exist=True)
            if path and Path(path).suffix.lower() in {".pdb", ".ent", ".cif", ".mmcif"}:
                found.setdefault(str(Path(path).resolve()), (key, Path(path)))
    # These are explicit archived input locations used by older adapters.
    for filename, role in (("receptor.pdb", "protein_pdb"), ("target.pdb", "target_pdb"),
                           ("complex.pdb", "complex_pdb")):
        path = Path(run_dir) / "artifacts" / "input" / filename
        if path.is_file() and not is_sequence_carrier_reference(path, run_dir=run_dir):
            found.setdefault(str(path.resolve()), (role, path))
    return list(found.values())


def _reference_file(database, pdb_id: str) -> Path:
    """Fetch a public PDB reference once, at registration rather than page render."""
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb_id):
        raise ValueError("Invalid source PDB accession.")
    path = database.root / "reference_sources" / (pdb_id + ".cif")
    if path.is_file():
        return path
    if time.monotonic() - _FETCH_FAILURES.get(pdb_id, -1000) < 300:
        raise OSError("The reference download is pending after a recent network failure.")
    try:
        filename, content = imports.fetch_identifier(pdb_id, "PDB")
        parsed = imports.parse_structure(content, filename)
        if not any(component.get("type") == "protein" for component in parsed[0]["components"]):
            raise ValueError("The source reference contains no protein coordinates.")
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            stream.write(content)
            temporary = Path(stream.name)
        temporary.replace(path)
    except (OSError, ValueError):
        _FETCH_FAILURES[pdb_id] = time.monotonic()
        raise
    return path


def collect_run_inputs(run_dir: Path, *, database=None, fetch_references=True) -> dict:
    """Idempotently retain authorized inputs; never import or evaluate outputs."""
    from mn_molecular_design_studio.collection_policy import automatic_collection_allowed
    from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id, _now
    from mn_molecular_design_studio.molecular_database_ingest import _ensure_group, _target_sequence_record
    from mn_molecular_design_studio.target_identity import target_reference_type
    from mn_molecular_design_studio.workspace import is_studio_job

    run_dir = Path(run_dir).resolve()
    report = {"molecule_ids": [], "structures": 0, "errors": [], "skipped": []}
    if not is_studio_job(run_dir) or not automatic_collection_allowed(run_dir):
        return report
    if run_dir.parent.name == "analysis":
        # Analysis consumes existing scientific representations. Its inputs
        # must not promote an intermediate design into the reusable library.
        return report
    marker = _json_file(run_dir / "molecular_design_studio.json")
    metadata = _json_file(run_dir / "metadata.json")
    embedded = metadata.get("molecular_design_studio") or {}
    owner = str(marker.get("studio_owner_id") or metadata.get("owner_id") or embedded.get("studio_owner_id") or "")
    if not owner:
        return report
    project = str(marker.get("project_id") or metadata.get("project_id") or embedded.get("project_id") or "")
    user = {"id": owner, "role": "admin"}
    database = database or MolecularDatabase()
    input_data = _json_file(run_dir / "input.json")
    input_values = input_data.get("inputs") if isinstance(input_data.get("inputs"), dict) else {}
    parameters = input_data.get("params") if isinstance(input_data.get("params"), dict) else {}
    derived_inputs = bool(parameters.get("validation_sources") or parameters.get("validation_mode")
                          or parameters.get("acceptance_policy") or input_data.get("source_job"))
    origin = {"type": "job", "task_group": run_dir.parent.name, "run_id": run_dir.name,
              "run_dir": str(run_dir), "workstation": socket.gethostname(), "collection_kind": "reference"}
    source_prefix = f"{socket.gethostname()}/{run_dir.parent.name}/{run_dir.name}:protein-input"
    group_id = None

    def register(record: dict, suffix: str, context: dict) -> str:
        nonlocal group_id
        if group_id is None:
            group_id = _ensure_group(database, user, _id("protein-library-", f"{owner}:{project}"),
                                     "Protein Library", project, "tool_inputs",
                                     {"type": "job", "owner_id": owner, "project_id": project}, "reference")
        _, fingerprint, _ = imports.canonical_components(record["components"])
        molecule_id = "mol-" + fingerprint[:32]
        # Prefer an existing human-supplied name in this same ownership scope.
        with database.connection() as connection:
            named = connection.execute(
                "SELECT mm.name FROM molecule_members mm JOIN molecule_groups g ON g.id=mm.group_id "
                "WHERE mm.molecule_id=? AND g.owner_id=? AND g.project_id=? "
                "AND g.import_method IN ('user_uploads','agent_inputs') "
                "ORDER BY (g.import_method='user_uploads') DESC,mm.added_at LIMIT 1",
                (molecule_id, owner, project)).fetchone()
            if named:
                record["name"] = named[0]
            latest = connection.execute(
                "SELECT mm.metadata FROM molecule_members mm JOIN molecule_groups g ON g.id=mm.group_id "
                "WHERE mm.molecule_id=? AND g.owner_id=? AND g.project_id=? ORDER BY mm.added_at DESC LIMIT 1",
                (molecule_id, owner, project)).fetchone()
            if latest:
                record["metadata"] = {**json.loads(latest[0]), **(record.get("metadata") or {})}
            molecule_id = database._put(connection, group_id, record, source_key=source_prefix + suffix,
                                         origin={**origin, **context}, stamp=_now())
        if molecule_id not in report["molecule_ids"]:
            report["molecule_ids"].append(molecule_id)
        report["structures"] += len(record.get("files") or [])
        return molecule_id

    def structure_record(path: Path, name: str) -> dict:
        content = path.read_bytes()
        record = imports.parse_structure(content, path.name)[0]
        if not any(component.get("type") == "protein" for component in record["components"]):
            raise _NotProteinInput("The input has no protein component.")
        record.update(name=name, files=[{"name": path.name, "content": content,
            "label": "Source reference structure", "structure_role": "source_reference", "portable": True,
            "portable_sha256": hashlib.sha256(content).hexdigest(), "dimensions": 3}])
        record["metadata"] = {"source_structure_available": True, "library_category": "protein_reference"}
        return record

    if target_reference_type(input_data, metadata, run_dir=run_dir) == "sequence_only":
        target_record = _target_sequence_record(input_data, metadata, run_dir)
        if target_record:
            target_record["metadata"]["library_category"] = "protein_reference"
            register(target_record, ":sequence-target", {"input_role": "target", "reference_type": "sequence_only"})

    for field, path in input_structure_paths(input_data, run_dir):
        if derived_inputs and field in {"complex_pdb", "source_pdb", "structure_path", "input_pdb", "pdb_path"}:
            continue
        try:
            record = structure_record(path, path.stem)
            register(record, ":file:" + hashlib.sha256(str(path).encode()).hexdigest()[:16],
                     {"input_field": field, "input_role": "target" if field.startswith("target_") else "source_structure",
                      "source_filename": path.name, "source_structure_digest": hashlib.sha256(path.read_bytes()).hexdigest()})
        except _NotProteinInput as exc:
            report["skipped"].append(f"{path.name}: {exc}")
        except (OSError, ValueError) as exc:
            report["errors"].append(f"{path.name}: {exc}")
    for index, task in enumerate(input_values.get("tasks") or []):
        if not isinstance(task, dict):
            continue
        if derived_inputs:
            continue
        target_text = str(task.get("source_target_pdb_text") or "")
        if target_text.strip():
            try:
                content = target_text.encode("utf-8")
                target_record = imports.parse_structure(content, "target.pdb")[0]
                target_record["name"] = str(task.get("name") or "Supplied protein structure")
                target_record["files"] = [{"name": "target.pdb", "content": content,
                    "label": "Source reference structure", "structure_role": "source_reference", "portable": True,
                    "portable_sha256": hashlib.sha256(content).hexdigest(), "dimensions": 3}]
                target_record["metadata"] = {"source_structure_available": True, "library_category": "protein_reference"}
                register(target_record, f":task-target:{index}", {"input_role": "target", "input_task": index})
            except ValueError as exc:
                report["errors"].append(f"Task {index + 1} target: {exc}")
        if any(task.get(key) for key in ("dna_chains", "rna_chains", "ligand_ccd_codes", "ligand_smiles")):
            continue
        chains = task.get("protein_chains") or []
        if not chains or not all(isinstance(chain, dict) and chain.get("sequence") for chain in chains):
            continue
        reference = protein_sequence_reference([chain["sequence"] for chain in chains])
        record = {"name": (f"{reference['pdb_id']} system · {len(chains)} protein chains" if reference else
                           str(task.get("name") or "Protein input system")),
                  "components": [{"type": "protein", "chain": str(chain.get("id") or chr(65 + i)),
                                  "sequence": chain["sequence"]} for i, chain in enumerate(chains)],
                  "metadata": {"library_category": "protein_reference"}, "files": []}
        if reference:
            record["metadata"]["sequence_source_reference"] = reference
            if fetch_references:
                try:
                    reference_path = _reference_file(database, reference["pdb_id"])
                    pdb_record = structure_record(reference_path, reference["pdb_id"] + " reference")
                    pdb_record["metadata"].update(source_pdb_id=reference["pdb_id"],
                                                  source_url=f"https://www.rcsb.org/structure/{reference['pdb_id']}")
                    reference_id = register(pdb_record, ":pdb:" + reference["pdb_id"],
                                            {"input_role": "sequence_source_reference", "pdb_id": reference["pdb_id"]})
                    record["metadata"]["sequence_source_reference_molecule_id"] = reference_id
                    record["metadata"]["reference_structure_status"] = "available"
                except (OSError, ValueError) as exc:
                    record["metadata"]["reference_structure_status"] = "download_pending"
                    report["errors"].append(f"{reference['pdb_id']}: {exc}")
        register(record, f":system:{index}", {"input_role": "prediction_system", "input_task": index,
                 "sequence_source_reference": reference, "chain_roles": {str(chain.get("id") or i):
                     str(chain.get("role") or "unspecified") for i, chain in enumerate(chains)}})
    return report


def record_run_inputs(run_dir: Path) -> None:
    """Library collection failure cannot change a scientific job's status."""
    try:
        report = collect_run_inputs(run_dir)
        if report["errors"]:
            log.warning("Protein Library input collection needs retry for %s: %s", run_dir, report["errors"])
    except Exception:
        log.exception("Protein Library input collection failed for %s; source inputs remain in the run", run_dir)


def audit_database(path: Path) -> dict:
    """Read a consistent inventory without migrations, reindexing, or writes."""
    connection = sqlite3.connect(f"file:{Path(path).resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("BEGIN")
    types = dict(connection.execute("SELECT type,COUNT(*) FROM molecules GROUP BY type"))
    rows = connection.execute("SELECT id,type,components FROM molecules WHERE type IN ('protein','complex')").fetchall()
    members = connection.execute("SELECT mm.molecule_id,mm.name,mm.metadata,g.origin,g.import_method FROM molecule_members mm JOIN molecule_groups g ON g.id=mm.group_id").fetchall()
    structures = connection.execute("SELECT id,molecule_id,path,context FROM molecule_structures").fetchall()
    connection.close()
    by_member, by_structure = {}, {}
    for member in members:
        by_member.setdefault(member["molecule_id"], []).append(dict(member))
    for structure in structures:
        by_structure.setdefault(structure["molecule_id"], []).append(dict(structure))
    counts = Counter()
    records = []
    for row in rows:
        components = json.loads(row["components"])
        if not any(component.get("type") == "protein" for component in components):
            continue
        counts["protein_containing_identities"] += 1
        memberships = by_member.get(row["id"], [])
        reps = by_structure.get(row["id"], [])
        counts["with_membership" if memberships else "without_membership"] += 1
        counts["with_structures" if reps else "without_structures"] += 1
        reference = (protein_sequence_reference([component["sequence"] for component in components])
                     if all(component.get("type") == "protein" and component.get("sequence") for component in components) else {})
        if reference:
            counts["known_pdb_sequence_reference_matches"] += 1
            counts["pdb_" + reference["pdb_id"]] += 1
        recorded = any((json.loads(member["metadata"]).get("sequence_source_reference") or {}).get("pdb_id")
                       for member in memberships)
        analysis_only = bool(memberships) and all(json.loads(member["origin"]).get("collection_kind") == "evaluation"
                                                  for member in memberships)
        if reference and not recorded and not analysis_only:
            counts["known_pdb_reference_missing_metadata"] += 1
        elif reference and not recorded:
            counts["known_pdb_analysis_only_matches"] += 1
        library = any(json.loads(member["origin"]).get("collection_kind") == "reference" for member in memberships)
        counts["in_protein_library"] += bool(library)
        local = remote = missing = 0
        for rep in reps:
            if str(rep["path"]).startswith("cluster-library://"):
                remote += 1
            elif Path(rep["path"]).is_file():
                local += 1
            else:
                missing += 1
        counts["local_structure_representations"] += local
        counts["remote_structure_references"] += remote
        counts["unavailable_local_structure_paths"] += missing
        records.append({"molecule_id": row["id"], "name": memberships[0]["name"] if memberships else "Unassociated identity",
                        "memberships": len(memberships), "local_structures": local, "remote_references": remote,
                        "unavailable_paths": missing, "pdb_sequence_reference": reference,
                        "reference_recorded": bool(recorded), "analysis_only": analysis_only,
                        "in_protein_library": bool(library)})
    return {"molecule_types": types, "counts": dict(counts), "records": records}
