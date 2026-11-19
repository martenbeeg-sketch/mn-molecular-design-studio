"""Exact molecular evidence for target registration; names never establish identity."""
from __future__ import annotations

import hashlib
import io
import json
import re
from pathlib import Path

from mn_molecular_design_studio.molecule_imports import canonical_components


def target_components(components, metadata=None, candidate=None):
    metadata, candidate = metadata or {}, candidate or {}
    chains = set(candidate.get("target_chains") or metadata.get("target_chains") or [])
    binder_chains = set(candidate.get("binder_chains") or metadata.get("binder_chains") or [])
    if chains & binder_chains:
        return []
    declared = [row for row in components if row.get("role") in {"target", "receptor"} or row.get("chain") in chains]
    if any(row.get("role") in {"binder", "designed_binder", "design"} or row.get("chain") in binder_chains for row in declared):
        return []
    mapping = (metadata.get("source_to_prediction_components") or {}).get("target") or []
    sequences = metadata.get("sequences_by_chain") or {}
    if not declared and chains and all(sequences.get(chain) for chain in chains):
        declared = [{"type": "protein", "chain": chain, "sequence": sequences[chain]} for chain in sorted(chains)]
    if not declared and mapping and all(row.get("prediction_sequence") and row.get("prediction_chain") for row in mapping):
        declared = [{"type": "protein", "chain": row["prediction_chain"], "sequence": row["prediction_sequence"]} for row in mapping]
    if not declared:
        return []
    try:
        normalized = canonical_components([{**row, "role": "target"} for row in declared])[0]
    except (KeyError, TypeError, ValueError):
        return []
    if chains and {row["chain"] for row in normalized} != chains:
        return []
    return normalized


def contents_identity(components):
    return tuple(sorted((row["type"], "sequence" if row.get("sequence") else "smiles" if row.get("smiles") else "ccd",
                        row.get("sequence") or row.get("smiles") or row.get("ccd") or "") for row in components))


def registration_evidence(components, original_target):
    return {"policy": "studio.exact-target-registration.v1", "original_target_id": str(original_target or ""),
        "basis": "explicit target chains and exact molecular contents",
        "contents_sha256": hashlib.sha256(json.dumps(contents_identity(components)).encode()).hexdigest()}


def reference_context(*sources, path=None):
    """Keep explicit reference accessions separate from exact molecular identity."""
    accessions = set()
    receptors = []
    for source in sources:
        if not isinstance(source, dict):
            continue
        reference = source.get("target_reference")
        reference = reference if isinstance(reference, dict) else {}
        values = [source.get(key) for key in ("target_pdb_id", "source_pdb_id", "pdb_id")]
        reference_ids = reference.get("pdb_ids")
        values.extend([reference.get("pdb_id"), *(reference_ids if isinstance(reference_ids, list) else [])])
        for value in values:
            if isinstance(value, str) and re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value.strip()):
                accessions.add(value.strip().upper())
        receptor = source.get("target_receptor") or source.get("receptor") or reference.get("receptor")
        if isinstance(receptor, str) and receptor.strip() and receptor.strip() not in receptors:
            receptors.append(receptor.strip())
    if path is not None:
        source_path = Path(path)
        if source_path.stat().st_size <= 100_000_000:
            if source_path.suffix.lower() in {".cif", ".mmcif"}:
                from Bio.PDB.MMCIF2Dict import MMCIF2Dict
                try:
                    document = MMCIF2Dict(io.StringIO(source_path.read_text()))
                    values = document.get("_entry.id", [])
                except (OSError, UnicodeError, ValueError, KeyError):
                    values = []
            else:
                with source_path.open() as stream:
                    header = stream.read(262144)
                values = [line[62:66].strip() for line in header.splitlines() if line.startswith("HEADER")]
            accessions.update(value.upper() for value in values
                if isinstance(value, str) and re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value))
    status = "recorded" if len(accessions) == 1 else "conflicting" if accessions else "unspecified"
    pdb_id = next(iter(accessions)) if len(accessions) == 1 else ""
    return {"target_pdb_id": pdb_id, "target_receptor": receptors[0] if receptors else "",
        "target_reference": {"schema": "studio.target-reference.v1", "status": status,
            "pdb_id": pdb_id, "pdb_ids": sorted(accessions), "receptor": receptors[0] if receptors else ""}}


def saved_reference_context(database, user, project_id, *, structure_id="", molecule_id="", path=None):
    """Read reference metadata only from an accessible, selected or byte-identical representation."""
    digest = ""
    if path is not None:
        with Path(path).open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if not structure_id and not digest:
        return reference_context()
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT s.*,mm.metadata FROM molecule_structures s "
            "JOIN molecule_members mm ON mm.group_id=s.group_id AND mm.molecule_id=s.molecule_id "
            "JOIN molecule_groups g ON g.id=s.group_id WHERE g.project_id=? AND "
            "((?<>'' AND s.id=? AND (?='' OR s.molecule_id=?)) OR (?='' AND ("
            "json_extract(s.context,'$.content_sha256')=? OR json_extract(s.context,'$.portable_sha256')=?)))",
            (str(project_id or ""), structure_id, structure_id, molecule_id, molecule_id,
             structure_id, digest, digest),
        ).fetchall()
    references = []
    for row in rows:
        try:
            database.group(row["group_id"], user)
            verified = database.structure_file(user, row["id"], materialize=True)
            if not structure_id:
                with verified.open("rb") as stream:
                    if hashlib.file_digest(stream, "sha256").hexdigest() != digest:
                        continue
            context = json.loads(row["context"])
            references.append(context)
            if reference_context(context)["target_reference"]["status"] == "unspecified":
                references.append(json.loads(row["metadata"]))
        except (OSError, ValueError, PermissionError, FileNotFoundError):
            continue
    return reference_context(*references)


def register_exact_target(database, user, project_id, components, *, name, source_key, origin, reference=None):
    """Register a target input once; never admit it as a generated design."""
    from mn_molecular_design_studio.molecular_database_ingest import _find_or_create_input

    components = canonical_components(components)[0]
    reference = reference_context(reference or {})
    evidence = registration_evidence(components, (reference or {}).get("target_pdb_id"))
    record = {"name": name, "components": components, "files": [],
        "metadata": {"target_identity_type": "sequence_only" if all(row["type"] == "protein" for row in components) else "molecular_identity", "structural_reference": False,
            **reference, "target_registration": evidence}}
    molecule_id, _, _ = _find_or_create_input(database, user, project_id, name, source_key, record,
        provenance={**origin, **reference, "input_role": "modeled_target", "target_registration": evidence},
        source_key_suffix=":target")
    return {**reference, "target_molecule_id": molecule_id, "target_registration": evidence}
