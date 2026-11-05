"""Idempotent adapters for Studio's existing molecular artifact contracts."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import logging
import math
import re
import socket
import threading
import time
from pathlib import Path

from mn_molecular_design_studio import molecule_imports as imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id, _json, _now

log = logging.getLogger(__name__)
_LAST_SYNC = {}
_LAST_SYNC_ERRORS = {}
_SYNC_WORKERS = set()
_SYNC_WORKERS_LOCK = threading.Lock()


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


def _target_structure_path(inputs, run_dir):
    """Resolve a supplied target structure, excluding sequence carrier CIFs."""
    from mn_molecular_design_studio.target_identity import is_sequence_carrier_reference

    for relative in ("artifacts/input/receptor.pdb", "artifacts/input/target.pdb", "artifacts/input/complex.pdb"):
        path = Path(run_dir) / relative
        if path.is_file() and not is_sequence_carrier_reference(path, run_dir=run_dir):
            return path
    mappings = [inputs.get("inputs") or {}, inputs.get("params") or {}, inputs]
    for mapping in mappings:
        if not isinstance(mapping, dict):
            continue
        for key in ("complex_pdb", "target_pdb", "protein_pdb", "target_structure"):
            path = _path(mapping.get(key), run_dir)
            if path is not None and not is_sequence_carrier_reference(path, run_dir=run_dir):
                return path
    return None


def _remove_sequence_carrier_structures(database, user, run_dir):
    """Avoid misclassified IDR carrier structures without touching curated rows."""
    from mn_molecular_design_studio.target_identity import is_sequence_carrier_reference, sequence_carrier_path

    if not sequence_carrier_path(_read(Path(run_dir) / "input.json"), run_dir=run_dir):
        return 0
    removed = 0
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT s.id,s.group_id,s.path,s.context,g.origin FROM molecule_structures s "
            "JOIN molecule_groups g ON g.id=s.group_id WHERE g.owner_id=?",
            (str(user.get("id") or ""),),
        ).fetchall()
    stale = []
    for row in rows:
        try:
            context = json.loads(row["context"] or "{}")
        except (TypeError, ValueError):
            continue
        if (str(context.get("run_id") or "") == Path(run_dir).name
                and is_sequence_carrier_reference(row["path"], run_dir=run_dir)):
            try:
                group_origin = json.loads(row["origin"] or "{}")
            except (TypeError, ValueError):
                group_origin = {}
            if group_origin.get("collection_kind") == "curated":
                continue
            stale.append((row["id"], row["group_id"]))
    if stale:
        with database.connection() as connection:
            for structure_id, group_id in stale:
                connection.execute("DELETE FROM molecule_structures WHERE id=? AND group_id=?",
                                   (structure_id, group_id))
                connection.execute("UPDATE molecule_groups SET updated_at=? WHERE id=?", (_now(), group_id))
                removed += 1
    return removed


def _legacy_attempt_target_context(run_dir, attempt_id, run_metadata):
    """Resolve the copied receptor and provenance for one bundled MN-Ligand attempt."""
    attempt_id = str(attempt_id or "").strip()
    if not attempt_id or Path(attempt_id).name != attempt_id:
        return {}
    attempt_dir = Path(run_dir) / "artifacts" / "replicates" / attempt_id
    if not attempt_dir.is_dir():
        return {}

    target_path = next(
        (path for relative in (
            "input/receptor.pdb",
            "artifacts/input/receptor.pdb",
            "input/target.pdb",
            "artifacts/input/target.pdb",
            "artifacts/prepared/receptor.pdb",
        ) if (path := attempt_dir / relative).is_file()),
        None,
    )
    if target_path is None:
        return {}

    attempt_inputs = _read(attempt_dir / "input.json")
    attempt_metadata = _read(attempt_dir / "metadata.json")
    target_ref = attempt_inputs.get("target_artifact") or attempt_inputs.get("prepared_target")
    if not isinstance(target_ref, dict):
        nested_inputs = attempt_inputs.get("inputs")
        nested_inputs = nested_inputs if isinstance(nested_inputs, dict) else {}
        target_ref = nested_inputs.get("target_artifact") or nested_inputs.get("prepared_target")
    target_ref = target_ref if isinstance(target_ref, dict) else {}
    artifact_metadata = target_ref.get("metadata")
    artifact_metadata = artifact_metadata if isinstance(artifact_metadata, dict) else {}
    migration_metadata = run_metadata.get("legacy_migration")
    migration_metadata = migration_metadata if isinstance(migration_metadata, dict) else {}
    attempt_migration = attempt_metadata.get("legacy_migration")
    attempt_migration = attempt_migration if isinstance(attempt_migration, dict) else {}
    pdb_id = str(
        run_metadata.get("pdb_id")
        or migration_metadata.get("pdb_id")
        or attempt_metadata.get("pdb_id")
        or attempt_migration.get("pdb_id")
        or ""
    ).upper()
    campaign_parts = [part.strip() for part in str(run_metadata.get("campaign_name") or "").split(" · ")]
    campaign_receptor = campaign_parts[0] if len(campaign_parts) > 1 and campaign_parts[1].upper() == pdb_id else ""
    receptor = str(
        run_metadata.get("receptor")
        or migration_metadata.get("receptor")
        or attempt_metadata.get("receptor")
        or attempt_migration.get("receptor")
        or campaign_receptor
    )
    role = str(target_ref.get("role") or target_ref.get("artifact_type") or "")
    digest = hashlib.sha256(target_path.read_bytes()).hexdigest()
    display_parts = [part for part in (receptor, pdb_id, role.replace("_", " ")) if part]
    if not display_parts:
        display_parts.append("MN-Ligand target")
    target_label = " · ".join(display_parts) + f" · {digest[:8]}"
    referenced_digest = str(target_ref.get("sha256") or "")

    return {
        "target_id": target_label,
        "target_label": target_label,
        "target_path": str(target_path.resolve()),
        "target_structure_digest": digest,
        "target_variant_sha256": digest,
        "target_pdb_id": pdb_id,
        "target_receptor": receptor,
        "target_source_application": "MN-Ligand",
        "target_source_attempt_run_id": attempt_id,
        "prepared_target_run_id": str(
            target_ref.get("run_id")
            or attempt_metadata.get("prepared_target_run_id")
            or attempt_metadata.get("parent_run_id")
            or ""
        ),
        "target_artifact_id": str(target_ref.get("artifact_id") or ""),
        "target_artifact_role": role,
        "target_artifact_sha256": referenced_digest,
        "target_artifact_sha256_matches_input": referenced_digest == digest if referenced_digest else None,
        "target_coordinate_transform": artifact_metadata.get("coordinate_transform"),
        "target_artifact": target_ref,
    }


def _record(candidate, run_dir):
    files = []
    components = None
    sequence = candidate.get("binder_sequence") or candidate.get("sequence")
    raw_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    rejected_files = []
    representation_type = str(raw_metadata.get("structure_kind") or raw_metadata.get("representation_type") or "unknown")
    # A candidate has one molecular identity. Do not attach target-only or
    # binder-only files as if they represented a different complex's identity.
    # Keep every available candidate representation on that identity, while
    # choosing the first valid structure (complex first) for component parsing.
    for key in ("complex_pdb", "binder_pdb", "target_pdb"):
        if key == "target_pdb" and sequence:
            # A target-only file cannot identify the newly designed sequence.
            continue
        path = _path(candidate.get(key), run_dir)
        if not path:
            continue
        try:
            structure_bytes = path.read_bytes()
            selection = raw_metadata.get("complex_analysis_selection") if isinstance(
                raw_metadata.get("complex_analysis_selection"), dict) else {}
            selected_model = _coordinate_model_selection(candidate, raw_metadata, selection)
            parsed = imports.parse_structure(structure_bytes, path.name, model_id=selected_model,
                                             require_model_selection=True)[0]
            parsed_components = parsed["components"]
        except (OSError, ValueError) as exc:
            rejected_files.append({"field": key, "path": str(path), "status": "unresolved_model_or_coordinates",
                                   "reason": str(exc)})
            continue
        valid, reason = _candidate_structure_correspondence(key, parsed_components, candidate, components)
        if not valid:
            rejected_files.append({"field": key, "path": str(path), "status": "molecular_mismatch", "reason": reason})
            continue
        if components is None:
            components = parsed_components
        if all(existing.get("path") != str(path) for existing in files):
            files.append({"path": str(path), "label": key.replace("_", " ").title(), "dimensions": 3,
                          "source_file_field": key, "representation_type": representation_type,
                          "model_id": parsed.get("model_id"),
                          "validated_content_sha256": hashlib.sha256(structure_bytes).hexdigest()})
    raw_components = (candidate.get("raw_metadata") or {}).get("components")
    validation_path = _path((candidate.get("raw_metadata") or {}).get("validation_complex_pdb"), run_dir)
    if validation_path and validation_path.is_file() and all(row.get("path") != str(validation_path) for row in files):
        try:
            validation_bytes = validation_path.read_bytes()
            selection = raw_metadata.get("complex_analysis_selection") if isinstance(
                raw_metadata.get("complex_analysis_selection"), dict) else {}
            selected_model = _coordinate_model_selection(candidate, raw_metadata, selection)
            validation_components = imports.parse_structure(
                validation_bytes, validation_path.name, model_id=selected_model,
                require_model_selection=True,
            )[0]["components"]
        except (OSError, ValueError) as exc:
            validation_components = None
            rejected_files.append({"field": "validation_complex_pdb", "path": str(validation_path),
                                   "status": "unresolved_model_or_coordinates", "reason": str(exc)})
        if validation_components:
            valid, reason = _candidate_structure_correspondence(
                "validation_complex_pdb", validation_components, candidate, components,
            )
            if not valid:
                rejected_files.append({"field": "validation_complex_pdb", "path": str(validation_path),
                                       "status": "molecular_mismatch", "reason": reason})
            elif components is None:
                components = validation_components
            if valid:
                files.append({"path": str(validation_path), "label": "Independent validation complex", "dimensions": 3,
                              "source_file_field": "validation_complex_pdb",
                              "representation_type": str(raw_metadata.get("validation_structure_kind") or "unknown"),
                              "model_id": str(selected_model or "0"),
                              "validated_content_sha256": hashlib.sha256(validation_bytes).hexdigest()})
    if components is None and isinstance(raw_components, list):
        try:
            components = imports.canonical_components(raw_components)[0]
        except (KeyError, TypeError, ValueError):
            components = None
    small_molecules = [component for component in components or [] if component.get("type") == "small_molecule"]
    raw_smiles = raw_metadata.get("ligand_smiles") or []
    if isinstance(raw_smiles, (str, dict)):
        raw_smiles = [raw_smiles]
    smiles_values = [
        str(value.get("smiles") or value.get("value") or "").strip() if isinstance(value, dict)
        else str(value or "").strip()
        for value in raw_smiles if isinstance(value, (str, dict))
    ]
    smiles_values = [value for value in smiles_values if value]
    if small_molecules and len(smiles_values) == len(small_molecules):
        for component, smiles in zip(small_molecules, smiles_values):
            component["smiles"] = smiles
    raw_ccd = raw_metadata.get("ligand_ccd_codes") or []
    if isinstance(raw_ccd, str):
        raw_ccd = [value for value in re.split(r"[\s,;]+", raw_ccd) if value]
    if small_molecules and not smiles_values and isinstance(raw_ccd, list) and len(raw_ccd) == len(small_molecules):
        for component, ccd in zip(small_molecules, raw_ccd):
            if isinstance(ccd, str) and ccd.strip():
                component["ccd"] = ccd.strip()
    if components is None and sequence:
        components = imports.parse_sequences(str(sequence))[0]["components"]
        designed_chains = list(candidate.get("binder_chains") or [])
        if len(components) == len(designed_chains):
            for component, chain in zip(components, designed_chains):
                component.update(chain=chain, role="binder")
    if not components:
        return None
    binder_chains = {str(chain) for chain in (candidate.get("binder_chains") or [])}
    target_chains = {str(chain) for chain in (candidate.get("target_chains") or [])}
    for component in components:
        # Structure parsers give polymer chains their literal chain ID, but
        # split non-polymers into synthetic IDs such as ``C_108`` to keep
        # residues unique. Their residue_key retains the source chain ID
        # (``C:108``), which is what candidate annotations reference.
        aliases = {str(component.get("chain") or "")}
        residue_key = str(component.get("residue_key") or "")
        if ":" in residue_key:
            aliases.add(residue_key.split(":", 1)[0])
        binder_match = bool(aliases & binder_chains)
        target_match = bool(aliases & target_chains)
        if binder_match != target_match:
            component["role"] = "binder" if binder_match else "target"
    # An unsequenced poly-glycine backbone is a structure proposal, not an
    # assertion that every generated backbone is the same designed molecule.
    unresolved = str(candidate.get("stage") or "") == "generation.backbone" and not sequence
    identity_context = None
    if unresolved and files:
        identity_context = {"unsequenced_backbone": hashlib.sha256(Path(files[0]["path"]).read_bytes()).hexdigest()}
    return {"name": raw_metadata.get("compound_name") or raw_metadata.get("source_compound_name")
            or raw_metadata.get("source_compound_id") or raw_metadata.get("input_name")
            or candidate.get("candidate_id") or "Molecule", "components": components, "files": files,
            "identity_context": identity_context,
            "metadata": {"status": raw_metadata.get("status") or ("Unsequenced backbone" if unresolved else ""),
                         "structure_validation": rejected_files,
                         **raw_metadata, "stage": candidate.get("stage"), "source_tool": candidate.get("source_tool"),
                         "output_kind": (candidate.get("normalized_output") or {}).get("kind")}}


def _candidate_structure_correspondence(field, components, candidate, primary_components=None):
    """Check a coordinate file against declared binder/target chains and sequences."""
    protein_by_chain = {str(item.get("chain") or ""): str(item.get("sequence") or "")
                        for item in components if item.get("type") == "protein"}
    if not protein_by_chain:
        return False, "No protein chain was parsed from the candidate representation."
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    ligand_analysis = str(raw.get("analysis_subject") or "").casefold() == "ligand"
    declared_components = raw.get("components") if isinstance(raw.get("components"), list) else []
    declared_proteins = {str(item.get("chain") or ""): str(item.get("sequence") or "").replace(" ", "").upper()
                         for item in declared_components
                         if isinstance(item, dict) and item.get("type", "protein") == "protein"
                         and item.get("sequence")}
    if field in {"complex_pdb", "validation_complex_pdb"} and declared_proteins \
            and protein_by_chain != declared_proteins:
        return False, "Protein chain IDs or sequences do not match the candidate's declared molecular composition."
    sequences_by_chain = raw.get("sequences_by_chain") if isinstance(raw.get("sequences_by_chain"), dict) else {}
    binder_chains = [str(value) for value in candidate.get("binder_chains") or []]
    target_chains = [str(value) for value in candidate.get("target_chains") or []]
    binder_sequence = str(candidate.get("binder_sequence") or candidate.get("sequence") or "")
    expected_binder = [part.upper() for part in re.split(r"[:/]", binder_sequence) if part]
    if binder_chains and field != "target_pdb" and not ligand_analysis:
        missing = sorted(set(binder_chains) - set(protein_by_chain))
        if missing:
            return False, "Declared binder chains are absent: " + ", ".join(missing)
        selected_binder = [protein_by_chain[chain] for chain in binder_chains]
        declared = [str(sequences_by_chain[chain]).replace(" ", "").upper()
                    for chain in binder_chains if chain in sequences_by_chain]
        expected = declared or expected_binder
        if expected and sorted(selected_binder) != sorted(expected):
            return False, "Binder chain sequence does not match the declared design sequence."
    elif expected_binder and field != "target_pdb" and not ligand_analysis \
            and sorted(protein_by_chain.values())[:len(expected_binder)] != sorted(expected_binder):
        # Without chain metadata, exact sequence matching is the only safe
        # way to identify which parsed component is the designed binder.
        if not all(sequence in protein_by_chain.values() for sequence in expected_binder):
            return False, "No parsed chain has the declared binder sequence."
    if target_chains and field != "binder_pdb":
        missing = sorted(set(target_chains) - set(protein_by_chain))
        if missing:
            return False, "Declared target chains are absent: " + ", ".join(missing)
        declared_targets = [str(sequences_by_chain[chain]).replace(" ", "").upper()
                            for chain in target_chains if chain in sequences_by_chain]
        if declared_targets and sorted(protein_by_chain[chain] for chain in target_chains) != sorted(declared_targets):
            return False, "Target chain sequence does not match the declared modeled target."
    if field == "binder_pdb":
        expected = expected_binder or [str(sequences_by_chain[chain]).replace(" ", "").upper()
                                       for chain in binder_chains if chain in sequences_by_chain]
        if expected and sorted(protein_by_chain.values()) != sorted(expected):
            return False, "Binder-only representation contains a different protein composition."
        if binder_chains and set(protein_by_chain) != set(binder_chains):
            return False, "Binder-only representation has different chain identifiers."
    if field in {"complex_pdb", "validation_complex_pdb"} and primary_components:
        primary = {str(item.get("chain") or ""): str(item.get("sequence") or "")
                   for item in primary_components if item.get("type") == "protein"}
        if primary and primary != protein_by_chain:
            return False, "Complex representation does not preserve the primary structure's protein composition."
    return True, ""


def _candidate_calculation_digest(record):
    """Return the validated digest for the exact coordinates analyzed by a candidate."""
    files = (record or {}).get("files") or []
    for field in ("complex_pdb", "binder_pdb", "validation_complex_pdb", "target_pdb"):
        match = next((item for item in files if item.get("source_file_field") == field
                      and item.get("validated_content_sha256")), None)
        if match:
            return str(match["validated_content_sha256"])
    return ""


def _coordinate_model_selection(candidate, raw_metadata, selection):
    """Return an actual coordinate-model selector, excluding representation IDs."""
    for mapping, keys in (
        (candidate, ("model_index",)),
        (raw_metadata, ("model_index",)),
        (selection, ("model_index",)),
        (candidate, ("model_id",)),
        (raw_metadata, ("model_id",)),
        (selection, ("model_id",)),
    ):
        if not isinstance(mapping, dict):
            continue
        for key in keys:
            value = mapping.get(key)
            if value is not None and re.fullmatch(r"\d+", str(value).strip()):
                return value
    return None


def _prediction_tasks_by_id(run_dir, inputs):
    """Associate older prediction candidates with their saved task inputs."""
    nested = inputs.get("inputs") if isinstance(inputs.get("inputs"), dict) else {}
    tasks = nested.get("tasks") if isinstance(nested.get("tasks"), list) else inputs.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        return {}
    tasks_by_name = {}
    for task in tasks:
        if not isinstance(task, dict):
            continue
        for key in ("name", "source_input_name", "source_compound_id"):
            if task.get(key):
                tasks_by_name[str(task[key])] = task
    result = {}
    for relative in (
        "artifacts/raw/optimized_prediction/input_mapping.json",
        "artifacts/raw/af3_jax_cofolding/input_mapping.json",
    ):
        mapping = _read(Path(run_dir) / relative)
        rows = mapping.get("items") if isinstance(mapping.get("items"), list) else []
        for row in rows:
            if not isinstance(row, dict) or not row.get("id"):
                continue
            task = next(
                (tasks_by_name[str(row[key])] for key in ("source_input_name", "name")
                 if row.get(key) and str(row[key]) in tasks_by_name),
                None,
            )
            if task is not None:
                result[str(row["id"])] = task
    return result


def _attach_prediction_input(candidate, tasks_by_id):
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    task = tasks_by_id.get(str(raw.get("batch_item_id") or ""))
    if not task:
        return candidate
    enriched_raw = dict(raw)
    task_name = str(task.get("source_input_name") or task.get("name") or "").strip()
    if task_name:
        enriched_raw.setdefault("input_name", task_name)
    compound_name = str(task.get("compound_name") or task_name).strip()
    if compound_name:
        enriched_raw.setdefault("compound_name", compound_name)
        enriched_raw.setdefault("source_compound_name", compound_name)
    compound_id = str(task.get("compound_id") or task.get("source_compound_id") or "").strip()
    if compound_id:
        enriched_raw.setdefault("compound_id", compound_id)
        enriched_raw.setdefault("source_compound_id", compound_id)
    ligand_smiles = task.get("ligand_smiles") or []
    if isinstance(ligand_smiles, (str, dict)):
        ligand_smiles = [ligand_smiles]
    if isinstance(ligand_smiles, list) and ligand_smiles:
        enriched_raw.setdefault("ligand_smiles", ligand_smiles)
        if len(ligand_smiles) == 1 and task_name:
            enriched_raw.setdefault("source_compound_id", task_name)
    ligand_ccd_codes = task.get("ligand_ccd_codes") or []
    if ligand_ccd_codes:
        enriched_raw.setdefault("ligand_ccd_codes", ligand_ccd_codes)
    return {**candidate, "raw_metadata": enriched_raw}


def _candidate_rejected(candidate):
    raw = candidate.get("raw_metadata") or {}
    metrics = candidate.get("metrics") or {}
    validation = raw.get("validation_acceptance")
    if isinstance(validation, dict) and validation.get("status") == "rejected":
        return True
    status = str(candidate.get("status") or raw.get("status") or raw.get("filter_status") or metrics.get("status") or "").strip().casefold()
    if status in {"rejected", "failed", "invalid", "filtered", "discarded", "not accepted", "not_accepted"}:
        return True
    campaign = raw.get("distributed_campaign") if isinstance(raw.get("distributed_campaign"), dict) else {}
    if campaign.get("validation_required") and campaign.get("validated_acceptance") != "accepted":
        # Campaign candidates remain visible in Results, but only designs that
        # passed the configured campaign filter belong in the reusable library.
        return True
    return (candidate.get("accepted") is False or raw.get("accepted") is False
            or metrics.get("pass_filters") is False or metrics.get("native_pass_filters") is False)


def _candidate_acceptance_state(candidate):
    """Return only explicit candidate decisions; missing evidence stays unknown."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
    campaign = raw.get("distributed_campaign") if isinstance(raw.get("distributed_campaign"), dict) else {}
    if campaign.get("validation_required"):
        campaign_status = str(campaign.get("validated_acceptance") or "").strip().casefold()
        if campaign_status in {"accepted", "accept", "passed", "pass"}:
            return "accepted"
        if campaign_status in {"rejected", "reject", "failed", "filtered", "discarded", "not_accepted"}:
            return "rejected"
        if campaign_status:
            return "unresolved"
        return "unknown"
    decision = raw.get("validation_acceptance") if isinstance(raw.get("validation_acceptance"), dict) else {}
    status = str(decision.get("status") or candidate.get("status") or raw.get("status") or "").strip().casefold()
    if status in {"accepted", "accept", "passed", "pass"}:
        return "accepted"
    if status in {"rejected", "reject", "failed", "filtered", "discarded", "not_accepted"}:
        return "rejected"
    if status in {"unresolved", "unknown", "not_evaluated", "pending"}:
        return "unresolved"
    if candidate.get("accepted") is True or raw.get("accepted") is True:
        return "accepted"
    if candidate.get("accepted") is False or raw.get("accepted") is False:
        return "rejected"
    if metrics.get("native_pass_filters") is True or metrics.get("pass_filters") is True:
        return "accepted"
    if metrics.get("native_pass_filters") is False or metrics.get("pass_filters") is False:
        return "rejected"
    return "unknown"


_ANALYSIS_WORKFLOW_MARKERS = (
    "analysis.", "analysis_", "biotite_interface", "pep_patch", "ipsae",
    "sequence_qc", "sequence_liability", "plip", "pandamap", "posebusters", "foldseek",
    "tag_fusion", "panda_map",
)


def _analysis_only_workflow(run_dir, metadata, inputs):
    values = [Path(run_dir).parent.name]
    for mapping in (metadata, inputs, inputs.get("params") or {}):
        if isinstance(mapping, dict):
            values.extend(str(mapping.get(key) or "") for key in
                          ("tool", "workflow", "workflow_id", "job_type", "engine", "task"))
    text = " ".join(values).casefold().replace("-", "_")
    return Path(run_dir).parent.name.casefold() == "analysis" or any(
        marker in text for marker in _ANALYSIS_WORKFLOW_MARKERS
    )


def _run_mode(run_dir, metadata, inputs):
    if run_dir.parent.name == "candidate-import":
        return "curated"
    values = [run_dir.parent.name]
    for mapping in (metadata, inputs, inputs.get("params") or {}):
        values.extend(str(mapping.get(key) or "") for key in ("tool", "workflow", "workflow_id", "job_type", "engine", "task"))
    text = " ".join(values).casefold().replace("-", "_")
    # Analysis workflows are Evaluation producers even when their registered
    # workflow ID does not contain the word "evaluation". Keep this explicit
    # so their candidates cannot fall through into identity-creating modes.
    if _analysis_only_workflow(run_dir, metadata, inputs):
        return "evaluation"
    if any(token in text for token in ("optimized_prediction", "refolding.optimized", "validation_acceptance", "vina", "docking", "affinity", "rescoring", "evaluation", "boltz2", "boltz_2", "molecule_selection")):
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


def _find_accessible_group(database, user, project, record):
    """Find the exact identity in any ACL-readable group, including replicas."""
    fingerprint = _canonical_fingerprint(record)
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT m.id,mm.group_id,g.project_id,g.updated_at FROM molecules m "
            "JOIN molecule_members mm ON mm.molecule_id=m.id JOIN molecule_groups g ON g.id=mm.group_id "
            "WHERE m.fingerprint=? AND g.project_id=? ORDER BY g.updated_at,mm.group_id",
            (fingerprint, str(project or "")),
        ).fetchall()
    for row in rows:
        try:
            database.group(str(row["group_id"]), user)
        except (FileNotFoundError, PermissionError, ValueError):
            continue
        return str(row["id"]), str(row["group_id"]), True
    return None, None, False


def _database_analysis_subject(database, user, owner, project, source):
    """Resolve a completed Database analysis to its exact original structure."""
    if not isinstance(source, dict) or source.get("schema") != "studio.complex_analysis_database_source.v1":
        return None
    molecule_id = str(source.get("molecule_id") or "")
    structure_id = str(source.get("structure_id") or "")
    group_id = str(source.get("group_id") or "")
    digest = str(source.get("structure_digest") or "")
    if not all((molecule_id, structure_id, group_id, digest)):
        return None
    try:
        group = database.group(group_id, user)
        if str(group.get("project_id") or "") != str(project or ""):
            return None
        if not str(project or "") and str(group.get("owner_id") or "") != str(owner):
            return None
        page = database.list_molecules(user, molecule_ids=[molecule_id], collection="all", limit=1)
        rows = page.get("items") if isinstance(page, dict) else []
        record = rows[0] if rows else None
        if record is None or not any(str(item.get("id") or "") == group_id for item in record.get("groups", [])):
            return None
        representation = next((item for item in record.get("representations", [])
                               if str(item.get("id") or "") == structure_id
                               and str(item.get("group_id") or "") == group_id), None)
        if representation is None:
            return None
        path = database.structure_file(user, structure_id, materialize=True)
        if _file_sha256(path) != digest:
            return None
        context = representation.get("context") or {}
        represented_context = dict(context)
        if (isinstance(source.get("ligand_selection"), dict)
                and not represented_context.get("ligand_selection")):
            represented_context["ligand_selection"] = source["ligand_selection"]
        representation_components, selected_model_id = _selected_representation_components(
            path.read_bytes(), path.name, represented_context,
        )
        _apply_representation_identity_annotations(representation_components, record.get("components") or [])
        target_id = context.get("design_target_molecule_id") if (
            source.get("target_chains") and
            set(source["target_chains"]) == set(context.get("target_chains") or [])
        ) else None
        origin = group.get("origin") if isinstance(group.get("origin"), dict) else {}
        if origin.get("_replica_node"):
            source_owner_node = str(origin.get("_replica_node") or "")
            source_owner_group_id = str(origin.get("_replica_group") or group_id)
            source_owner_structure_id = str(context.get("_replica_structure") or structure_id)
        else:
            try:
                from mn_molecular_design_studio.cluster.config import identity as node_identity
                source_owner_node = str(origin.get("workstation") or origin.get("source_workstation") or node_identity())
            except (OSError, RuntimeError, ValueError):
                source_owner_node = str(origin.get("workstation") or origin.get("source_workstation") or "")
            source_owner_group_id = group_id
            source_owner_structure_id = structure_id
        return {"molecule_id": molecule_id, "structure_id": structure_id, "group_id": group_id,
                "structure_digest": digest, "record": record,
                "representation": {**representation, "context": represented_context},
                "representation_path": str(path), "representation_components": representation_components,
                "representation_model_id": selected_model_id,
                "source_reference_kind": "database",
                "target_molecule_id": target_id, "source_owner_id": str(group.get("owner_id") or ""),
                "source_node": source_owner_node,
                "source_structure_owner_node": source_owner_node,
                "source_structure_owner_group_id": source_owner_group_id,
                "source_structure_owner_structure_id": source_owner_structure_id}
    except (FileNotFoundError, OSError, PermissionError, ValueError, TypeError):
        return None


def _selected_representation_components(content, filename, context):
    """Read composition from the exact saved coordinate model, then apply its provenance."""
    context = context if isinstance(context, dict) else {}
    model_id = next((context.get(key) for key in ("model_id", "model_index", "model")
                     if context.get(key) is not None and str(context.get(key)).strip()), None)
    parsed = imports.parse_structure(
        content, filename, include_ligands=True, model_id=model_id,
        require_model_selection=True,
    )[0]
    components = [dict(item) for item in parsed.get("components") or []]
    binder_chains = {str(value) for value in context.get("binder_chains") or []}
    target_chains = {str(value) for value in context.get("target_chains") or []}
    if binder_chains & target_chains:
        raise ValueError("Saved representation has overlapping binder and target chain assignments.")
    ligand_selection = context.get("ligand_selection") if isinstance(context.get("ligand_selection"), dict) else {}
    for component in components:
        chain = str(component.get("chain") or "")
        if component.get("type") == "protein":
            if chain in binder_chains:
                component["role"] = "binder"
            elif chain in target_chains:
                component["role"] = "target"
        if component.get("type") == "small_molecule":
            if (ligand_selection and
                    str(component.get("residue_key") or "") == str(ligand_selection.get("residue_key") or "") and
                    (not ligand_selection.get("ccd") or
                     str(component.get("ccd") or "").casefold() == str(ligand_selection["ccd"]).casefold())):
                component["smiles"] = str(ligand_selection.get("smiles") or context.get("ligand_isomeric_smiles") or "")
    return components, str(parsed.get("model_id") or "")


def _apply_representation_identity_annotations(components, identity_components):
    """Carry stored chemistry labels onto coordinates only by exact residue correspondence."""
    for component in components:
        if component.get("type") != "small_molecule" or component.get("smiles"):
            continue
        matches = [identity for identity in identity_components
                   if identity.get("type") == "small_molecule"
                   and str(identity.get("chain") or "") == str(component.get("chain") or "")
                   and str(identity.get("residue_key") or "") == str(component.get("residue_key") or "")
                   and (not identity.get("ccd") or not component.get("ccd")
                        or str(identity.get("ccd") or "").casefold() == str(component.get("ccd") or "").casefold())
                   and identity.get("smiles")]
        if len(matches) == 1:
            component["smiles"] = str(matches[0]["smiles"])


def _canonical_smiles(value):
    if not value:
        return ""
    try:
        from rdkit import Chem
        molecule = Chem.MolFromSmiles(str(value))
        return Chem.MolToSmiles(molecule, isomericSmiles=True) if molecule else ""
    except (ImportError, TypeError, ValueError):
        return ""


def _representation_contains_identity(identity_components, representation_components):
    """Confirm canonical identity components occur in the selected representation."""
    remaining = [dict(item) for item in representation_components]
    for identity in identity_components:
        kind = str(identity.get("type") or "")
        if kind not in {"protein", "small_molecule", "dna", "rna"}:
            return False
        candidates = [item for item in remaining if str(item.get("type") or "") == kind]
        if kind in {"protein", "dna", "rna"}:
            sequence = str(identity.get("sequence") or "").replace(" ", "").upper()
            candidates = [item for item in candidates
                          if str(item.get("sequence") or "").replace(" ", "").upper() == sequence]
            role = str(identity.get("role") or "").casefold()
            if role and any(str(item.get("role") or "").casefold() == role for item in candidates):
                candidates = [item for item in candidates if str(item.get("role") or "").casefold() == role]
        else:
            smiles = _canonical_smiles(identity.get("smiles"))
            if smiles:
                candidates = [item for item in candidates if _canonical_smiles(item.get("smiles")) == smiles]
            elif identity.get("ccd"):
                ccd = str(identity.get("ccd") or "").casefold()
                candidates = [item for item in candidates if str(item.get("ccd") or "").casefold() == ccd]
            else:
                return False
        identity_chain = str(identity.get("chain") or "")
        exact_chain = [item for item in candidates if identity_chain and
                       str(item.get("chain") or "") == identity_chain]
        if exact_chain:
            candidates = exact_chain
        if not candidates:
            return False
        remaining.remove(candidates[0])
    return True


def _protein_chain_composition(components):
    result = {}
    for item in components:
        if item.get("type") != "protein":
            continue
        chain = str(item.get("chain") or "")
        sequence = str(item.get("sequence") or "").replace(" ", "").upper()
        if not chain or not sequence or chain in result:
            return {}
        result[chain] = sequence
    return result


def _mapped_protein_compositions_match(source_map, calculation_map, chain_map):
    if not source_map or not calculation_map:
        return False
    if chain_map is None:
        return source_map == calculation_map
    if not isinstance(chain_map, dict) or set(map(str, chain_map)) != set(source_map):
        return False
    mapped = {str(chain_map.get(chain) or ""): sequence for chain, sequence in source_map.items()}
    return bool(all(mapped) and len(mapped) == len(set(mapped)) and mapped == calculation_map)


def _ligand_geometry(path, residue_key):
    """Return heavy-atom composition and intramolecular distances for one mapped residue."""
    from Bio.PDB import PDBParser, MMCIFParser
    try:
        path = Path(path)
        filename = path.name[:-3] if path.name.casefold().endswith(".gz") else path.name
        content = path.read_bytes()
        if path.name.casefold().endswith(".gz"):
            content = gzip.decompress(content)
        parser = MMCIFParser(QUIET=True) if Path(filename).suffix.casefold() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
        structure = parser.get_structure("ligand-correspondence", io.StringIO(content.decode("utf-8")))
        models = list(structure.get_models())
        if len(models) != 1:
            return None
        chain_id, _, position = str(residue_key).partition(":")
        match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", position)
        if not match:
            return None
        number, insertion = int(match.group(1)), match.group(2)
        chain = models[0].child_dict.get(chain_id)
        if chain is None:
            return None
        residues = [residue for residue in chain
                    if int(residue.id[1]) == number and str(residue.id[2]).strip() == insertion]
        if len(residues) != 1:
            return None
        atoms = []
        for atom in residues[0].get_atoms():
            element = str(atom.element or atom.get_name()[:1]).strip().upper()
            if element in {"H", "D"}:
                continue
            xyz = tuple(float(value) for value in atom.coord)
            if not all(math.isfinite(value) for value in xyz):
                return None
            atoms.append((element, xyz))
        if not atoms:
            return None
        distances = []
        for index, (_, left) in enumerate(atoms):
            for _, right in atoms[index + 1:]:
                distances.append(round(math.sqrt(sum((left[axis] - right[axis]) ** 2 for axis in range(3))), 1))
        return tuple(sorted(element for element, _ in atoms)), tuple(sorted(distances))
    except (OSError, ValueError, TypeError, KeyError):
        return None


def _analysis_calculation_matches_source(source, candidate, record):
    """Compare calculation coordinates with the exact selected representation, not its owner identity alone."""
    if not isinstance(source, dict) or not isinstance(candidate, dict):
        return False
    selected_calculation_field = next((field for field in ("complex_pdb", "binder_pdb")
                                       if candidate.get(field)), "")
    rejected_structures = ((record or {}).get("metadata") or {}).get("structure_validation") or []
    if selected_calculation_field and any(
            str(item.get("field") or "") == selected_calculation_field
            for item in rejected_structures if isinstance(item, dict)):
        return False
    source_components = source.get("representation_components") or []
    identity_components = (source.get("record") or {}).get("components") or []
    calculation_components = (record or {}).get("components") or []
    if (not source_components or not calculation_components
            or not _representation_contains_identity(identity_components, source_components)):
        return False
    source_proteins = [item for item in source_components if item.get("type") == "protein"]
    calculation_proteins = [item for item in calculation_components if item.get("type") == "protein"]
    if not source_proteins or not calculation_proteins:
        return False
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
    source_representation = source.get("representation") if isinstance(source.get("representation"), dict) else {}
    source_context = source_representation.get("context") if isinstance(source_representation.get("context"), dict) else {}
    ligand_selection = selection.get("ligand_selection") if isinstance(selection.get("ligand_selection"), dict) else {}
    if not ligand_selection and isinstance(source_context.get("ligand_selection"), dict):
        ligand_selection = source_context["ligand_selection"]
    ligand_analysis = str(raw.get("analysis_subject") or "").casefold() == "ligand" or bool(ligand_selection)
    chain_map = raw.get("source_to_calculation_chain_map")
    source_protein_map = _protein_chain_composition(source_proteins)
    calculation_protein_map = _protein_chain_composition(calculation_proteins)
    if not _mapped_protein_compositions_match(source_protein_map, calculation_protein_map, chain_map):
        return False
    if not ligand_analysis:
        return True

    residue_key = str(ligand_selection.get("residue_key") or "")
    if not residue_key:
        return False
    source_ligand_chain, _, residue_position = residue_key.partition(":")
    calculation_ligand_chain = str(chain_map.get(source_ligand_chain) or "") if isinstance(chain_map, dict) else ""
    calculation_residue_key = (
        calculation_ligand_chain + ":" + residue_position
        if calculation_ligand_chain else residue_key
    )
    source_ligands = [item for item in source_components if item.get("type") == "small_molecule"
                      and str(item.get("residue_key") or "") == residue_key]
    calculation_ligands = [item for item in calculation_components if item.get("type") == "small_molecule"
                            and str(item.get("residue_key") or "") == calculation_residue_key]
    requested_ccd = str(ligand_selection.get("ccd") or "").casefold()
    if requested_ccd:
        source_ligands = [item for item in source_ligands
                          if str(item.get("ccd") or "").casefold() == requested_ccd]
        calculation_ligands = [item for item in calculation_ligands
                                if str(item.get("ccd") or "").casefold() == requested_ccd]
    if len(source_ligands) != 1 or len(calculation_ligands) != 1:
        return False
    ligand_smiles = str(ligand_selection.get("smiles") or source_ligands[0].get("smiles") or "")
    if not ligand_smiles:
        return False
    calculation_smiles = str(calculation_ligands[0].get("smiles") or "")
    if calculation_smiles and _canonical_smiles(calculation_smiles) != _canonical_smiles(ligand_smiles):
        return False
    try:
        from rdkit import Chem
        molecule = Chem.MolFromSmiles(ligand_smiles)
        if molecule is None:
            return False
        expected_elements = tuple(sorted(atom.GetSymbol().upper() for atom in molecule.GetAtoms()))
    except (ImportError, TypeError, ValueError):
        return False
    source_path = source.get("representation_path")
    calculation_path = next((item.get("path") for item in (record.get("files") or [])
                             if item.get("source_file_field") == selected_calculation_field and item.get("path")), None)
    if not calculation_path:
        calculation_path = next((item.get("path") for item in (record.get("files") or []) if item.get("path")), None)
    source_geometry = _ligand_geometry(source_path, residue_key) if source_path else None
    calculation_geometry = _ligand_geometry(calculation_path, calculation_residue_key) if calculation_path else None
    return bool(source_geometry and calculation_geometry
                and source_geometry == calculation_geometry
                and source_geometry[0] == expected_elements)


def _results_analysis_subject(database, user, owner, project, analysis_run, inputs, candidate, metadata=None):
    """Resolve Results analysis to one exact source structure or sequence identity."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    source_candidate = raw.get("source_candidate") if isinstance(raw.get("source_candidate"), dict) else {}
    source_candidate = source_candidate or candidate
    params = inputs.get("params") if isinstance(inputs.get("params"), dict) else {}
    nested = inputs.get("inputs") if isinstance(inputs.get("inputs"), dict) else {}
    workflow_values = []
    for mapping in (metadata or {}, inputs, params):
        if isinstance(mapping, dict):
            workflow_values.extend(str(mapping.get(key) or "") for key in
                                   ("tool", "workflow", "workflow_id", "job_type", "engine", "task"))
    workflow_text = " ".join(workflow_values).casefold().replace("-", "_")
    sequence_only = any(name in workflow_text for name in ("sequence_qc", "sequence_liability"))
    source_ref = inputs.get("source_run_dir") or params.get("source_run_dir") or nested.get("source_run_dir")
    if not source_ref:
        # Sequence analyses may be launched directly against an already
        # selected Database molecule without a Results run directory. Accept
        # that explicit source only when its saved protein sequence exactly
        # matches the analyzed binder; never infer an identity from the
        # analysis output's own coordinates.
        molecule_id = str(raw.get("studio_molecule_id") or source_candidate.get("studio_molecule_id") or "")
        expected_sequence = str(source_candidate.get("binder_sequence") or source_candidate.get("sequence") or "")
        expected_sequences = [part.replace(" ", "").upper() for part in re.split(r"[:/]", expected_sequence) if part]
        if sequence_only and molecule_id and expected_sequences:
            selected_chains = {str(value) for value in source_candidate.get("binder_chains") or []}
            with database.connection() as connection:
                source_rows = connection.execute(
                    "SELECT m.components,g.id AS group_id,g.owner_id,g.project_id,g.origin "
                    "FROM molecules m JOIN molecule_members mm ON mm.molecule_id=m.id "
                    "JOIN molecule_groups g ON g.id=mm.group_id WHERE m.id=? "
                    "AND g.project_id=? ORDER BY g.updated_at, g.id",
                    (molecule_id, str(project or "")),
                ).fetchall()
            for source in source_rows:
                if not database._can_read({"owner_id": source["owner_id"], "project_id": source["project_id"]}, user):
                    continue
                try:
                    source_components = json.loads(source["components"])
                except (TypeError, ValueError):
                    continue
                proteins = [item for item in source_components if item.get("type") == "protein"]
                selected = ([item for item in proteins if str(item.get("chain") or "") in selected_chains]
                            if selected_chains else [item for item in proteins
                                                     if str(item.get("role") or "").casefold() == "binder"])
                if not selected and not selected_chains:
                    selected = [item for item in proteins if str(item.get("sequence") or "").upper()
                                in expected_sequences]
                actual_sequences = [str(item.get("sequence") or "").replace(" ", "").upper()
                                    for item in selected]
                if sorted(actual_sequences) != sorted(expected_sequences):
                    continue
                try:
                    source_origin = json.loads(source["origin"] or "{}")
                except (TypeError, ValueError):
                    source_origin = {}
                return {
                    "molecule_id": molecule_id, "group_id": str(source["group_id"]),
                    "structure_id": "", "structure_digest": "", "source_reference_kind": "sequence",
                    "source_sequence_sha256": hashlib.sha256(":".join(sorted(expected_sequences)).encode()).hexdigest(),
                    "target_molecule_id": str(raw.get("target_molecule_id")
                                              or source_candidate.get("target_molecule_id") or ""),
                    "source_owner_id": str(source["owner_id"] or ""),
                    "source_node": str(source_origin.get("_replica_node") or ""),
                }
        return None
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    source_dir = resolve_stored_path(source_ref, must_exist=True)
    if source_dir is None:
        source_path = Path(str(source_ref)).expanduser()
        if not source_path.is_absolute():
            source_path = Path(analysis_run) / source_path
        try:
            source_dir = source_path.resolve(strict=True)
        except OSError:
            return None
    source_dir = Path(source_dir)
    if not source_dir.is_dir():
        return None
    source_metadata = _read(source_dir / "metadata.json")
    source_run_id = str(source_metadata.get("run_id") or source_dir.name)
    source_task_group = str(source_metadata.get("task_group") or source_dir.parent.name)
    source_candidate_id = str(source_candidate.get("candidate_id") or raw.get("analysis_source_candidate_id") or "")
    if not source_candidate_id:
        return None

    selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
    digest = str(selection.get("source_structure_digest") or selection.get("structure_digest")
                 or raw.get("source_structure_digest") or "")
    if not sequence_only and not digest:
        for key in ("complex_pdb", "binder_pdb"):
            path = _path(source_candidate.get(key), source_dir)
            if path and path.is_file():
                digest = _file_sha256(path)
                if digest:
                    break
    if not sequence_only and len(digest) != 64:
        return None
    expected_sequence = str(source_candidate.get("binder_sequence") or source_candidate.get("sequence") or "")
    expected_sequences = [part.replace(" ", "").upper() for part in re.split(r"[:/]", expected_sequence) if part]
    if sequence_only and not expected_sequences:
        return None

    source_marker = _read(source_dir / "molecular_design_studio.json")
    workstation = str(source_marker.get("workstation") or source_metadata.get("workstation") or "")
    with database.connection() as connection:
        sources = connection.execute(
            "SELECT s.molecule_id,s.group_id,s.provenance,g.owner_id,g.project_id,g.origin "
            "FROM molecule_sources s JOIN molecule_groups g ON g.id=s.group_id "
            "WHERE g.owner_id=? AND g.project_id=? AND json_extract(s.provenance,'$.run_id')=? "
            "AND json_extract(s.provenance,'$.task_group')=? AND json_extract(s.provenance,'$.candidate_id')=?",
            (str(owner), str(project or ""), source_run_id, source_task_group, source_candidate_id),
        ).fetchall()
        matches = []
        for source in sources:
            provenance = json.loads(source["provenance"])
            if workstation and str(provenance.get("workstation") or "") != workstation:
                continue
            group_access = {"owner_id": source["owner_id"], "project_id": source["project_id"]}
            if not database._can_read(group_access, user):
                continue
            if sequence_only:
                molecule = connection.execute(
                    "SELECT components FROM molecules WHERE id=?", (source["molecule_id"],)
                ).fetchone()
                if molecule is None:
                    continue
                try:
                    components = json.loads(molecule["components"])
                except (TypeError, ValueError):
                    continue
                proteins = [item for item in components if item.get("type") == "protein"]
                selected_chains = [str(value) for value in source_candidate.get("binder_chains") or []]
                if selected_chains:
                    selected = [item for item in proteins if str(item.get("chain") or "") in selected_chains]
                else:
                    role_binders = [item for item in proteins if str(item.get("role") or "").casefold() == "binder"]
                    selected = role_binders or [
                        item for item in proteins if str(item.get("sequence") or "").upper() in expected_sequences
                    ]
                selected_sequences = [str(item.get("sequence") or "").replace(" ", "").upper()
                                      for item in selected]
                if sorted(selected_sequences) != sorted(expected_sequences):
                    continue
                matches.append({
                    "molecule_id": str(source["molecule_id"]), "group_id": str(source["group_id"]),
                    "structure_id": "", "structure_digest": "", "source_reference_kind": "sequence",
                    "source_sequence_sha256": hashlib.sha256(":".join(sorted(expected_sequences)).encode()).hexdigest(),
                    "target_molecule_id": str(raw.get("target_molecule_id")
                                              or source_candidate.get("target_molecule_id") or ""),
                    "source_owner_id": str(source["owner_id"] or ""),
                    "source_node": str((json.loads(source["origin"] or "{}")).get("_replica_node") or ""),
                })
                continue
            molecule = connection.execute(
                "SELECT components FROM molecules WHERE id=?", (source["molecule_id"],)
            ).fetchone()
            if molecule is None:
                continue
            try:
                source_components = json.loads(molecule["components"])
            except (TypeError, ValueError):
                continue
            structures = connection.execute(
                "SELECT id,path,context FROM molecule_structures WHERE molecule_id=? AND group_id=?",
                (source["molecule_id"], source["group_id"]),
            ).fetchall()
            for structure in structures:
                context = json.loads(structure["context"])
                actual_digest = str(context.get("portable_sha256") or "") or _file_sha256(structure["path"])
                if actual_digest == digest:
                    represented_context = dict(context)
                    selected_ligand = raw.get("complex_analysis_selection")
                    selected_ligand = selected_ligand if isinstance(selected_ligand, dict) else {}
                    if selected_ligand.get("ligand_selection") and not represented_context.get("ligand_selection"):
                        represented_context["ligand_selection"] = selected_ligand["ligand_selection"]
                    try:
                        stored_path = str(structure["path"])
                        representation_path = (
                            database.structure_file(user, str(structure["id"]))
                            if stored_path.startswith("cluster-library://") else Path(stored_path)
                        )
                        representation_components, selected_model_id = _selected_representation_components(
                            representation_path.read_bytes(), representation_path.name, represented_context,
                        )
                        _apply_representation_identity_annotations(representation_components, source_components)
                    except (OSError, ValueError, TypeError):
                        representation_path = None
                        representation_components, selected_model_id = [], ""
                    source_origin = json.loads(source["origin"] or "{}")
                    if source_origin.get("_replica_node"):
                        source_owner_node = str(source_origin["_replica_node"])
                        source_owner_group = str(source_origin.get("_replica_group") or source["group_id"])
                        source_owner_structure = str(context.get("_replica_structure") or structure["id"])
                    else:
                        source_owner_node = str(provenance.get("workstation") or workstation)
                        if not source_owner_node:
                            try:
                                from mn_molecular_design_studio.cluster.config import identity as node_identity
                                source_owner_node = str(node_identity())
                            except (OSError, RuntimeError, ValueError):
                                source_owner_node = ""
                        source_owner_group = str(source["group_id"])
                        source_owner_structure = str(structure["id"])
                    matches.append({
                        "molecule_id": str(source["molecule_id"]), "group_id": str(source["group_id"]),
                        "structure_id": str(structure["id"]), "structure_digest": digest,
                        "source_reference_kind": "structure", "source_sequence_sha256": "",
                        "target_molecule_id": str(raw.get("target_molecule_id")
                                                  or source_candidate.get("target_molecule_id") or ""),
                        "source_owner_id": str(source["owner_id"] or ""),
                        "source_node": source_owner_node,
                        "source_structure_owner_node": source_owner_node,
                        "source_structure_owner_group_id": source_owner_group,
                        "source_structure_owner_structure_id": source_owner_structure,
                        "record": {"components": source_components},
                        "representation": {"id": str(structure["id"]), "context": represented_context},
                        "representation_path": str(representation_path) if representation_path else "",
                        "representation_components": representation_components,
                        "representation_model_id": selected_model_id,
                    })
    return matches[0] if len(matches) == 1 else None


def _file_sha256(path):
    try:
        with Path(path).open("rb") as stream:
            return hashlib.file_digest(stream, "sha256").hexdigest()
    except OSError:
        return ""


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


def _merge_member_roles(database, group_id, molecule_id, components):
    """Fill missing chain roles when an equivalent candidate joins a group."""
    if not components:
        return
    roles = {str(row.get("chain") or ""): str(row.get("role") or "")
             for row in components if row.get("role")}
    if not roles:
        return
    with database.connection() as connection:
        row = connection.execute(
            "SELECT components FROM molecule_members WHERE group_id=? AND molecule_id=?",
            (group_id, molecule_id),
        ).fetchone()
        if not row:
            return
        current = json.loads(row["components"])
        changed = False
        for component in current:
            role = roles.get(str(component.get("chain") or ""))
            if role and not component.get("role"):
                component["role"] = role
                changed = True
        if changed:
            connection.execute(
                "UPDATE molecule_members SET components=? WHERE group_id=? AND molecule_id=?",
                (_json(current), group_id, molecule_id),
            )


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


def _target_sequence_components(candidate=None, inputs=None, metadata=None, run_dir=None):
    """Return target protein chains only when their role is explicit."""
    raw = candidate.get("raw_metadata") if isinstance(candidate, dict) else {}
    raw = raw if isinstance(raw, dict) else {}
    target_chains = candidate.get("target_chains") if isinstance(candidate, dict) else None
    target_chains = target_chains or raw.get("target_chains") or []
    if not isinstance(target_chains, (list, tuple)):
        target_chains = []
    sequences_by_chain = raw.get("sequences_by_chain")
    if isinstance(sequences_by_chain, dict) and target_chains:
        rows = [{"type": "protein", "chain": str(chain), "role": "target",
                 "sequence": sequences_by_chain.get(str(chain))}
                for chain in target_chains if sequences_by_chain.get(str(chain))]
        if rows:
            try:
                return imports.canonical_components(rows)[0]
            except (KeyError, TypeError, ValueError):
                pass

    raw_components = raw.get("components")
    if isinstance(raw_components, list):
        rows = [dict(component, role="target") for component in raw_components
                if isinstance(component, dict) and str(component.get("role") or "").casefold() == "target"
                and str(component.get("type") or "protein").casefold() == "protein"]
        if rows:
            try:
                return imports.canonical_components(rows)[0]
            except (KeyError, TypeError, ValueError):
                pass

    # Only standalone run inputs are considered here. Do not interpret a
    # generic sequence field or infer target identity from a binder sequence.
    mappings = [metadata or {}, inputs or {}, (inputs or {}).get("inputs") or {},
                (inputs or {}).get("params") or {}]
    for mapping in mappings:
        if not isinstance(mapping, dict):
            continue
        value = next((mapping.get(key) for key in
                      ("target_sequence", "target_protein_sequences", "target_fasta")
                      if mapping.get(key) is not None), None)
        if value is None:
            continue
        rows = []
        if (isinstance(value, dict) and any(key in value for key in ("sequence", "fasta", "path"))):
            text = str(value.get("sequence") or value.get("fasta") or value.get("path") or "")
            source = _path(text, Path(run_dir)) if run_dir is not None else None
            if source and source.is_file():
                try:
                    text = source.read_text(encoding="utf-8-sig")
                except (OSError, UnicodeError):
                    continue
            try:
                parsed = imports.parse_sequences(text)
            except (TypeError, ValueError):
                continue
            rows = [{**component, "role": "target"} for parsed_row in parsed
                    for component in parsed_row.get("components") or []]
        elif isinstance(value, dict) and value and all(isinstance(sequence, str) for sequence in value.values()):
            rows = [{"type": "protein", "chain": str(chain), "role": "target", "sequence": sequence}
                    for chain, sequence in value.items()]
        elif isinstance(value, list):
            rows = [{"type": "protein", "chain": str(item.get("chain") or chr(65 + index)),
                     "role": "target", "sequence": item.get("sequence") or item.get("fasta") or ""}
                    if isinstance(item, dict) else {"type": "protein", "chain": chr(65 + index),
                                                    "role": "target", "sequence": str(item)}
                    for index, item in enumerate(value)]
        else:
            text = str(value.get("sequence") or value.get("fasta") or value.get("path") or "") if isinstance(value, dict) else str(value)
            source = _path(text, Path(run_dir)) if run_dir is not None else None
            if source and source.is_file():
                try:
                    text = source.read_text(encoding="utf-8-sig")
                except (OSError, UnicodeError):
                    continue
            try:
                parsed = imports.parse_sequences(text)
            except (TypeError, ValueError):
                continue
            rows = [{**component, "role": "target"} for parsed_row in parsed
                    for component in parsed_row.get("components") or []]
        if rows:
            try:
                return imports.canonical_components(rows)[0]
            except (KeyError, TypeError, ValueError):
                continue

    # RFD4's IDR preset stages an ideal backbone CIF only so the native
    # condition file can address target residues by sequence. Preserve that
    # sequence as the target identity, but never import its carrier coordinates
    # as a structural reference.
    if run_dir is not None:
        from mn_molecular_design_studio.target_identity import sequence_carrier_path

        carrier = sequence_carrier_path(inputs or {}, run_dir=run_dir)
        if carrier:
            try:
                parsed = imports.parse_structure(carrier.read_bytes(), carrier.name)[0]
                rows = [dict(component, role="target") for component in parsed.get("components") or []
                        if component.get("type") == "protein" and component.get("sequence")]
                if rows:
                    return imports.canonical_components(rows)[0]
            except (OSError, ValueError, KeyError, TypeError):
                pass
    return []


def _validated_target_structure(path, inputs, metadata, run_dir, expected_record=None):
    """Parse one explicit target model and verify its declared chains/sequences."""
    mappings = [metadata or {}, inputs or {}, (inputs or {}).get("inputs") or {},
                (inputs or {}).get("params") or {}]
    requested_chains = []
    selected_model = None
    for mapping in mappings:
        if not isinstance(mapping, dict):
            continue
        if not requested_chains and mapping.get("target_chains"):
            requested_chains = [str(chain) for chain in mapping["target_chains"]]
        if selected_model is None:
            for key in ("target_model_index", "target_model_id", "model_index", "model_id"):
                value = mapping.get(key)
                if value is not None and re.fullmatch(r"\d+", str(value).strip()):
                    selected_model = value
                    break
    structure_bytes = path.read_bytes()
    parsed = imports.parse_structure(structure_bytes, path.name, model_id=selected_model,
                                     require_model_selection=True)[0]
    components = parsed.get("components") or []
    proteins = {str(item.get("chain") or ""): str(item.get("sequence") or "")
                for item in components if item.get("type") == "protein" and item.get("sequence")}
    if not proteins:
        raise ValueError("The target structure has no protein coordinates.")
    if requested_chains and set(requested_chains) - set(proteins):
        raise ValueError("Declared target chains are absent from the selected coordinate model.")
    selected_sequences = ([proteins[chain] for chain in requested_chains] if requested_chains
                          else list(proteins.values()))
    expected = sorted(str(item.get("sequence") or "").replace(" ", "").upper()
                      for item in (expected_record or {}).get("components", [])
                      if item.get("type") == "protein")
    if expected:
        if requested_chains:
            if sorted(selected_sequences) != expected:
                raise ValueError("Target coordinate chain sequences do not match the declared target sequence.")
        else:
            unmatched = list(proteins.values())
            mapped = []
            for sequence in expected:
                matches = [value for value in unmatched if value == sequence]
                if len(matches) != 1:
                    raise ValueError("Target sequence does not map uniquely to the coordinate model.")
                unmatched.remove(sequence)
                mapped.append(sequence)
            selected_sequences = mapped
    elif requested_chains:
        selected_sequences = [proteins[chain] for chain in requested_chains]
    else:
        requested_chains = list(proteins)
    return {**parsed, "files": [{"path": str(path), "label": "Input receptor", "dimensions": 3,
                                 "validated_content_sha256": hashlib.sha256(structure_bytes).hexdigest()}],
            "metadata": {"target_identity_type": "structure_backed", "structural_reference": True,
                         "validated_model_id": str(parsed.get("model_id") or ""),
                         "validated_target_chains": requested_chains,
                         "validated_target_sequence_sha256": hashlib.sha256(
                             ":".join(sorted(selected_sequences)).encode()).hexdigest()}}


def _target_sequence_record(inputs, metadata, run_dir):
    components = _target_sequence_components(inputs=inputs, metadata=metadata, run_dir=run_dir)
    if not components:
        return None
    label = next((str(mapping.get(key) or "").strip() for mapping in
                  (metadata or {}, inputs or {}, (inputs or {}).get("inputs") or {}, (inputs or {}).get("params") or {})
                  if isinstance(mapping, dict) for key in ("target_name", "target_label", "target_id")
                  if mapping.get(key)), "Sequence-defined target")
    sequence_hash = hashlib.sha256(":".join(sorted(component["sequence"] for component in components)).encode()).hexdigest()
    return {"name": label, "components": components, "files": [],
            "metadata": {"target_identity_type": "sequence_only", "target_sequence_sha256": sequence_hash,
                         "structural_reference": False}}


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


def _is_cofolding_prediction(tool, metadata, inputs, mode):
    if mode != "evaluation":
        return False
    text = " ".join(str(mapping.get(key) or "") for mapping in (metadata, inputs, inputs.get("params") or {})
                    for key in ("tool", "workflow", "workflow_id", "engine", "task")).casefold().replace("-", "_")
    text = f"{text} {tool}".casefold().replace("-", "_")
    return not any(token in text for token in ("docking", "vina", "gnina", "rosetta_ligand", "unidock", "autodock")) and any(
        token in text for token in ("optimized_prediction", "boltz", "af3", "alphafold3", "chai", "cofold", "esmfold")
    )


def _prediction_engine_label(tool):
    value = str(tool or "").casefold().replace("_", "-")
    if "boltz" in value:
        return "Boltz-2"
    if "alphafold3" in value or "af3" in value:
        return "AlphaFold 3"
    return str(tool or "Prediction").replace("_", " ").title()


def _enrich_prediction_candidates(run_dir, candidates):
    """Read the same confidence sidecars used by the cofolding results page."""
    from mn_molecular_design_studio.cofolding_results import enrich_candidate

    return [enrich_candidate(run_dir, candidate) for candidate in candidates]


def _best_docking_pose_indices(rows, pose_paths):
    """Choose one pose per compound, preferring explicit pose rank then score."""
    def number(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def sort_key(index):
        row = rows[index]
        rank = next((number(row.get(key)) for key in ("pose_rank", "rank") if number(row.get(key)) is not None), None)
        score = None
        score_direction = 1
        for key in ("estimated_dg_reu", "best_score_kcal_mol", "score", "total_score_reu", "cnn_score", "cnn_affinity"):
            value = number(row.get(key))
            if value is not None:
                score = value
                score_direction = -1 if key in {"cnn_score", "cnn_affinity"} else 1
                break
        # Engine-provided pose rank is the primary ordering. A score resolves
        # replicate ties; without a rank, use the tool's score direction.
        if rank is not None:
            return (0, rank, score * score_direction if score is not None else 0.0, index)
        return (1, score * score_direction if score is not None else 0.0, index)

    candidates = {}
    for index, path in enumerate(pose_paths):
        compound_id = str(rows[index].get("compound_id") or "").strip()
        if path is not None and compound_id:
            row = rows[index]
            source_run_id = str(row.get("legacy_run_id") or "").strip()
            replicate = str(row.get("replicate") or "").strip()
            # Grouped MN-Ligand jobs retain source attempt IDs in their
            # merged score table. Keep one best pose per compound per attempt
            # (and per explicit intra-attempt replicate when provided).
            key = (compound_id, source_run_id, replicate) if source_run_id else (compound_id,)
            candidates.setdefault(key, []).append(index)
    return {min(indices, key=sort_key) for indices in candidates.values()}


def _prune_stale_run_structures(database, task_group, run_id, keep_ids, *, candidate_ids=None, docking=False):
    """Drop superseded indexed models/poses while preserving saved curated copies."""
    candidate_ids = {str(value) for value in (candidate_ids or []) if value}
    keep_ids = set(keep_ids)
    if not candidate_ids and not docking:
        return
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT s.id,s.context,g.origin FROM molecule_structures s "
            "JOIN molecule_groups g ON g.id=s.group_id"
        ).fetchall()
        for row in rows:
            context = json.loads(row["context"])
            if (str(context.get("task_group") or "") != str(task_group)
                    or str(context.get("run_id") or "") != str(run_id)):
                continue
            try:
                group_origin = json.loads(row["origin"])
            except (TypeError, ValueError):
                group_origin = {}
            if group_origin.get("collection_kind") == "curated":
                continue
            is_candidate = str(context.get("candidate_id") or "") in candidate_ids
            is_pose = ("pose_rank" in context or "compound_id" in (context.get("row") or {})
                       or context.get("database_structure_policy") == "best_docking_pose")
            if (is_candidate or (docking and is_pose)) and row["id"] not in keep_ids:
                connection.execute("DELETE FROM molecule_structures WHERE id=?", (row["id"],))


def _prune_empty_evaluation_members(database, memberships):
    """Remove old automatic evaluation entries after their data is rehomed."""
    if not memberships:
        return
    with database.connection() as connection:
        for group_id, molecule_id in memberships:
            group = connection.execute(
                "SELECT origin FROM molecule_groups WHERE id=?", (group_id,)
            ).fetchone()
            if not group:
                continue
            try:
                origin = json.loads(group["origin"])
            except (TypeError, ValueError):
                origin = {}
            if origin.get("type") != "job" or origin.get("collection_kind") != "evaluation":
                continue
            remaining = any(connection.execute(
                f"SELECT 1 FROM {table} WHERE group_id=? AND molecule_id=? LIMIT 1",
                (group_id, molecule_id),
            ).fetchone() for table in ("molecule_sources", "molecule_structures", "molecule_observations"))
            if not remaining:
                connection.execute(
                    "DELETE FROM molecule_members WHERE group_id=? AND molecule_id=?",
                    (group_id, molecule_id),
                )


def _run_input_compounds(run_dir):
    """Read a run's staged compound set, including TSV-only Rosetta inputs."""
    json_path = Path(run_dir) / "artifacts/input/compound_set.json"
    if json_path.is_file():
        try:
            from mn_molecular_design_studio.ligand_compounds import load_compound_set
            compounds = list(load_compound_set(json_path).get("compounds") or [])
            if compounds:
                return compounds
        except (OSError, ValueError, KeyError, TypeError):
            log.info("Could not read staged compound set %s", json_path)
    for path in (Path(run_dir) / "artifacts/input/compounds.tsv", Path(run_dir) / "input/compounds.tsv"):
        if not path.is_file():
            continue
        try:
            with path.open(newline="", encoding="utf-8-sig") as handle:
                rows = csv.DictReader(handle, delimiter="\t")
                compounds = []
                for row in rows:
                    compound_id = str(row.get("compound_id") or row.get("id") or "").strip()
                    smiles = str(row.get("smiles") or row.get("SMILES") or "").strip()
                    if not compound_id or not smiles:
                        continue
                    compounds.append({
                        "compound_id": compound_id,
                        "compound_name": str(row.get("compound_name") or row.get("name") or compound_id).strip(),
                        "smiles": smiles,
                        "properties": {},
                    })
                if compounds:
                    return compounds
        except (OSError, UnicodeError, csv.Error):
            log.info("Could not read run compound table %s", path)
    return []


def _candidate_structure_signature(run_dir):
    """Track referenced structure files so late-written outputs are re-indexed."""
    manifest = Path(run_dir) / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    try:
        lines = manifest.read_text(errors="ignore").splitlines()
    except OSError:
        return []
    signature = []
    for line in lines:
        try:
            candidate = json.loads(line)
        except (TypeError, ValueError):
            continue
        if not isinstance(candidate, dict):
            continue
        references = [(field, candidate.get(field)) for field in ("complex_pdb", "binder_pdb", "target_pdb")]
        raw_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        references.append(("validation_complex_pdb", raw_metadata.get("validation_complex_pdb")))
        for field, reference in references:
            if not reference:
                continue
            path = _path(reference, Path(run_dir))
            try:
                stat = path.stat() if path else None
            except OSError:
                stat = None
            signature.append({
                "candidate_id": str(candidate.get("candidate_id") or ""),
                "field": field,
                "reference": str(reference),
                "path": str(path) if path else None,
                "size": stat.st_size if stat else None,
                "mtime_ns": stat.st_mtime_ns if stat else None,
            })
    return signature


def ingest_run(run_dir: Path, *, database=None, force=False):
    from mn_molecular_design_studio.workspace import is_studio_job
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    run_dir = Path(run_dir)
    if not is_studio_job(run_dir):
        return 0
    from mn_molecular_design_studio.collection_policy import automatic_collection_allowed
    if not automatic_collection_allowed(run_dir):
        return 0
    if (run_dir / "artifacts/md/study.json").is_file():
        from mn_molecular_design_studio.md_database import ingest_completed_md_run
        return ingest_completed_md_run(run_dir, database=database, force=force)
    database = database or MolecularDatabase()
    marker = _read(run_dir / "molecular_design_studio.json")
    metadata, inputs, result = (_read(run_dir / name) for name in ("metadata.json", "input.json", "result.json"))
    embedded = metadata.get("molecular_design_studio") or {}
    owner = str(marker.get("studio_owner_id") or metadata.get("owner_id") or embedded.get("studio_owner_id") or "legacy-unassigned")
    project = str(marker.get("project_id") or metadata.get("project_id") or embedded.get("project_id") or "")
    user = {"id": owner, "role": "admin"}  # Trusted job adapter; all reads still enforce group access.
    tool = str(metadata.get("tool") or inputs.get("tool") or run_dir.parent.name)
    mode = _run_mode(run_dir, metadata, inputs)
    analysis_only = _analysis_only_workflow(run_dir, metadata, inputs)
    needs_compound_identity_refresh = run_dir.parent.name == "ligand-docking" or _is_cofolding_prediction(tool, metadata, inputs, mode)
    is_legacy_aggregate = (run_dir / "mn_ligand_aggregate.json").is_file()
    source_key = f"{socket.gethostname()}/{run_dir.parent.name}/{run_dir.name}"
    observed_paths = [run_dir / name for name in (
        "input.json", "result.json", "metadata.json", "molecular_design_studio.json", "selection_snapshot.json",
        "artifacts/normalized_candidates/candidates.jsonl", "artifacts/input/compound_set.json",
        "artifacts/input/compounds.tsv", "input/compounds.tsv",
        "artifacts/generated_compound_set.json", "artifacts/generated_compounds.sdf",
        "artifacts/docking_scores.csv", "artifacts/nesso_affinity.csv", "artifacts/input/receptor.pdb",
        "artifacts/input/target.pdb", "artifacts/input/complex.pdb",
        "artifacts/raw/optimized_prediction/input_mapping.json",
        "artifacts/raw/af3_jax_cofolding/input_mapping.json")]
    if is_legacy_aggregate:
        replicate_root = run_dir / "artifacts" / "replicates"
        if replicate_root.is_dir():
            for replicate_dir in sorted(path for path in replicate_root.iterdir() if path.is_dir()):
                observed_paths.extend(replicate_dir / relative for relative in (
                    "input.json",
                    "metadata.json",
                    "input/receptor.pdb",
                    "artifacts/input/receptor.pdb",
                    "input/target.pdb",
                    "artifacts/input/target.pdb",
                    "artifacts/prepared/receptor.pdb",
                ))
    database_analysis_run = str(metadata.get("upstream_source_kind") or "") == "database"
    ingestion_version = 11 if database_analysis_run else (10 if is_legacy_aggregate else (9 if needs_compound_identity_refresh else 8))
    signature = _json({"ingestion_version": ingestion_version, "files": [
        (str(path.relative_to(run_dir)), path.stat().st_mtime_ns, path.stat().st_size)
        for path in observed_paths if path.is_file()],
        "candidate_structures": _candidate_structure_signature(run_dir)})
    with database.connection() as connection:
        previous = connection.execute("SELECT signature FROM molecular_ingestions WHERE source_key=?", (source_key,)).fetchone()
    if not force and previous and previous[0] == signature:
        return 0
    if not _is_completed_success(metadata, result):
        with database.connection() as connection:
            connection.execute("INSERT INTO molecular_ingestions VALUES (?,?,?) ON CONFLICT(source_key) DO UPDATE SET signature=excluded.signature,updated_at=excluded.updated_at",
                               (source_key, signature, _now()))
        return 0

    name = str(metadata.get("campaign_name") or (inputs.get("params") or {}).get("job_name") or metadata.get("job_name") or run_dir.name)
    stamp = str(metadata.get("completed_at") or metadata.get("created_at") or _now())
    origin = {"type": "job", "run_id": run_dir.name, "task_group": run_dir.parent.name, "run_dir": str(run_dir.resolve()),
              "tool": tool, "workstation": socket.gethostname(), "status": metadata.get("status")}
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

    task_inputs_by_id = _prediction_tasks_by_id(run_dir, inputs)
    candidates = [
        _attach_prediction_input(candidate, task_inputs_by_id)
        for candidate in read_candidates(run_dir)
    ]
    prediction_run = _is_cofolding_prediction(tool, metadata, inputs, mode)
    best_prediction_ids = set()
    if prediction_run and candidates:
        candidates = _enrich_prediction_candidates(run_dir, candidates)
    candidate_records = []
    for candidate in candidates:
        record = _record(candidate, run_dir)
        if record is not None:
            candidate_records.append((candidate, record))
    input_target_record = _target_sequence_record(inputs, metadata, run_dir)
    _remove_sequence_carrier_structures(database, user, run_dir)
    target_path = _target_structure_path(inputs, run_dir)
    source_target_sequences = []
    target_structure_record = None
    if target_path:
        try:
            target_structure_record = _validated_target_structure(
                target_path, inputs, metadata, run_dir, expected_record=input_target_record,
            )
            source_target_sequences = [str(component.get("sequence") or "")
                for component in target_structure_record.get("components") or []
                if component.get("type") == "protein" and component.get("sequence")
                and (not target_structure_record["metadata"]["validated_target_chains"]
                     or str(component.get("chain") or "") in target_structure_record["metadata"]["validated_target_chains"])]
        except (OSError, ValueError, KeyError, TypeError) as exc:
            log.warning("Skipped target structure for %s because its model or molecular correspondence is invalid: %s",
                        run_dir, exc)
            target_path = None
            target_structure_record = None
            source_target_sequences = []
    if prediction_run and candidate_records:
        from mn_molecular_design_studio.cofolding_results import best_candidates_by_complex
        with_structures = [candidate for candidate, record in candidate_records if record.get("files")]
        best_prediction_ids = {id(candidate) for candidate in best_candidates_by_complex(with_structures)}

    input_compounds = _run_input_compounds(run_dir)
    score_paths = [run_dir / "artifacts" / filename for filename in ("docking_scores.csv", "nesso_affinity.csv")]
    md_metrics_supported = "md" in run_dir.parent.name and bool(result.get("metrics"))
    if (not candidate_records and not input_target_record
            and not (input_compounds and any(path.is_file() for path in score_paths)) and not md_metrics_supported):
        with database.connection() as connection:
            connection.execute("INSERT INTO molecular_ingestions VALUES (?,?,?) ON CONFLICT(source_key) DO UPDATE SET signature=excluded.signature,updated_at=excluded.updated_at",
                               (source_key, signature, _now()))
        return 0

    # Register generated and explicitly imported candidates in the library.
    # Backbone-only, rejected, and unmatched intermediate proposals stay in
    # the all-indexed view; evaluation outputs reuse an exact input identity.
    preferred_prediction_structure_ids = set()
    prediction_candidate_ids = set()
    generated_molecule_ids = set()
    stale_evaluation_memberships = set()
    with database.connection() as connection:
        removed_group = connection.execute(
            "SELECT origin FROM molecule_groups WHERE id=?", (_id("job-", source_key),)
        ).fetchone()
    removed_ids = set(json.loads(removed_group[0]).get("removed_molecule_ids") or []) if removed_group else set()
    for candidate, record in candidate_records:
        stage = str(candidate.get("stage") or "")
        sequence = candidate.get("binder_sequence") or candidate.get("sequence")
        is_backbone = stage == "generation.backbone" and not sequence
        raw_candidate = candidate.get("raw_metadata") or {}
        acceptance = raw_candidate.get("validation_acceptance")
        acceptance = acceptance if isinstance(acceptance, dict) else {}
        acceptance_state = _candidate_acceptance_state(candidate)
        # A saved Studio acceptance decision is a separate policy from native
        # generator diagnostics (for example, the native PPI diagnostic may
        # fail while a saved campaign policy accepts the candidate).
        rejected = acceptance_state == "rejected"
        is_library_design = (
            mode in {"generated", "curated"}
            or (mode == "auto" and stage in {"sequence_design", "generation.backbone_sequence"})
        )
        if ((mode == "generated" or (mode == "auto" and is_library_design))
                and acceptance_state != "accepted"):
            # Full generator history remains in Results. Unknown, rejected,
            # and unresolved outputs do not become reusable Database rows.
            continue
        if (acceptance_state == "accepted" and candidate.get("binder_sequence")
                and is_library_design and not (record.get("files") or [])):
            # An accepted binder must retain its design coordinates in a
            # reusable collection. Keep the run's evidence intact, but defer
            # Database promotion until its referenced structure is available.
            # _candidate_structure_signature includes unresolved references,
            # so the next workspace sync retries once the source file lands.
            log.warning(
                "Deferred accepted design %s from %s: its source structure is not available yet.",
                candidate.get("candidate_id") or "<unknown>", run_dir,
            )
            continue
        database_source = raw_candidate.get("database_source") if isinstance(raw_candidate, dict) else None
        raw_candidate_metadata = raw_candidate.get("raw_metadata") if isinstance(raw_candidate, dict) else None
        raw_candidate_metadata = raw_candidate_metadata if isinstance(raw_candidate_metadata, dict) else {}
        if database_source is None:
            database_source = raw_candidate_metadata.get("database_source")
        if database_source is None and isinstance(metadata.get("database_source"), dict):
            database_source = metadata["database_source"]
        database_subject = None
        if database_source is not None:
            database_subject = _database_analysis_subject(database, user, owner, project, database_source)
            selection = raw_candidate.get("complex_analysis_selection") or raw_candidate_metadata.get("complex_analysis_selection")
            selection = selection if isinstance(selection, dict) else {}
            source_selection_digest = str(
                selection.get("source_structure_digest") or selection.get("structure_digest") or ""
            )
            # Evaluators may protonate, normalize, or convert the selected
            # structure before emitting their analysis complex. That output
            # is expected to differ byte-for-byte from the saved Database
            # representation. Validate lineage against the original source
            # digest instead; _database_analysis_subject has already checked
            # the exact accessible representation and its stored digest.
            if (database_subject is None or not (record.get("files") or [])
                    or (source_selection_digest and source_selection_digest != database_subject["structure_digest"])):
                log.warning("Skipped Database evaluation projection for %s: its source structure is no longer an exact accessible representation.", run_dir)
                continue
        elif (mode == "evaluation" and not prediction_run
              and (analysis_only or run_dir.parent.name == "analysis" or stage in {"analysis", "benchmark"})):
            database_subject = _results_analysis_subject(database, user, owner, project, run_dir, inputs, candidate,
                                                         metadata=metadata)
            if database_subject is None:
                log.warning("Skipped Database evaluation projection for %s: its Results source did not resolve to one exact saved molecule, structure, or sequence.", run_dir)
                continue
            if database_subject.get("source_reference_kind") != "sequence":
                if not (record.get("files") or []):
                    log.warning("Skipped Database evaluation projection for %s: its source structure has no validated analysis representation.", run_dir)
                    continue
                selection = raw_candidate.get("complex_analysis_selection") or raw_candidate_metadata.get("complex_analysis_selection")
                selection = selection if isinstance(selection, dict) else {}
                source_selection_digest = str(
                    selection.get("source_structure_digest") or selection.get("structure_digest") or ""
                )
                if source_selection_digest and source_selection_digest != database_subject["structure_digest"]:
                    log.warning("Skipped Database evaluation projection for %s: its selected structure digest does not match the saved source representation.", run_dir)
                    continue
        if database_subject is not None and database_subject.get("source_reference_kind") != "sequence":
            if not _analysis_calculation_matches_source(database_subject, candidate, record):
                log.warning("Skipped Database evaluation projection for %s: analyzed molecular chains do not correspond to the exact selected source molecule.", run_dir)
                continue
        if analysis_only and database_subject is None:
            # Analysis output belongs in Results unless an exact existing
            # Database or Results subject was resolved above. In particular,
            # do not fall through to an Evaluation/intermediate identity.
            log.warning("Skipped Database association for %s: analysis source is unresolved.", run_dir)
            continue
        campaign_filter = raw_candidate.get("distributed_campaign") if isinstance(raw_candidate.get("distributed_campaign"), dict) else {}
        if mode == "generated" and campaign_filter.get("validation_required") and rejected:
            # Campaign validation is an eligibility gate, not an intermediate
            # library observation. Keep failures in campaign Results only.
            continue
        acceptance_evidence = {
            "validation_acceptance": raw_candidate.get("validation_acceptance"),
            "distributed_campaign": raw_candidate.get("distributed_campaign"),
            "native_filter_metrics": {key: (candidate.get("metrics") or {}).get(key)
                                      for key in ("native_pass_filters", "pass_filters")
                                      if key in (candidate.get("metrics") or {})},
        }
        provenance = {**origin, "candidate_id": candidate.get("candidate_id"),
                      "parents": candidate.get("parents") or [],
                      "acceptance_decision": {"status": acceptance_state, "evidence": acceptance_evidence}}
        candidate_id = str(candidate.get("candidate_id") or _canonical_fingerprint(record))
        source_for_molecule = f"{source_key}:{candidate_id}"
        with database.connection() as connection:
            previous_source = connection.execute(
                "SELECT s.molecule_id,s.group_id,s.provenance,mm.metadata FROM molecule_sources s "
                "LEFT JOIN molecule_members mm ON mm.group_id=s.group_id AND mm.molecule_id=s.molecule_id "
                "WHERE s.source_key=?",
                (source_for_molecule,),
            ).fetchone()
        if previous_source:
            previous_provenance = json.loads(previous_source["provenance"] or "{}")
            previous_decision = previous_provenance.get("acceptance_decision")
            if not previous_decision and previous_source["metadata"]:
                previous_metadata = json.loads(previous_source["metadata"] or "{}")
                prior_validation = previous_metadata.get("validation_acceptance")
                prior_campaign = previous_metadata.get("distributed_campaign")
                if prior_validation or prior_campaign:
                    previous_decision = {"status": prior_validation.get("status") if prior_validation else None,
                                         "evidence": {"validation_acceptance": prior_validation,
                                                      "distributed_campaign": prior_campaign}}
            history = list(previous_provenance.get("acceptance_history") or [])
            if previous_decision and previous_decision != provenance["acceptance_decision"]:
                if previous_decision not in history:
                    history.append(previous_decision)
                history.append(provenance["acceptance_decision"])
            if history:
                provenance["acceptance_history"] = history
        prediction_candidate_ids.add(candidate_id)
        identity_record = {**record, "files": []}
        designed_target_components = []
        ligand_analysis = (candidate.get("raw_metadata") or {}).get("analysis_subject") == "ligand"
        protein_analysis = (candidate.get("raw_metadata") or {}).get("analysis_subject") == "protein"
        prediction_ligand = prediction_run and sum(
            component.get("type") == "small_molecule" for component in record.get("components") or []
        ) == 1
        if prediction_ligand:
            # Match the prediction to its imported ligand identity, while the
            # full complex remains attached below as the prediction structure.
            identity_record["components"] = [
                component for component in record["components"]
                if component.get("type") == "small_molecule"
            ]
        if ligand_analysis:
            identity_record["components"] = [component for component in record["components"] if component.get("type") == "small_molecule"]
        if protein_analysis:
            designed = [component for component in record["components"] if component.get("chain") in (candidate.get("binder_chains") or []) and component.get("type") == "protein"]
            selected_components = (candidate.get("raw_metadata") or {}).get("components")
            if _explicit_molecule_ids(candidate) and isinstance(selected_components, list) and selected_components:
                # Frozen Database selections preserve whether the original
                # identity was a complete complex or just its designed chains.
                identity_record["components"] = selected_components
            elif designed:
                identity_record["components"] = designed
        if mode == "generated" and sequence and not rejected:
            # Generated designs are identified by their designed chain(s).
            # A later structure parse may add the target chain to the same
            # candidate; that complex should be attached as a representation
            # without creating a second molecule identity for the binder.
            binder_chains = {str(chain) for chain in candidate.get("binder_chains") or []}
            binder_components = [
                component for component in record.get("components") or []
                if component.get("role") == "binder" or str(component.get("chain")) in binder_chains
            ]
            if binder_components:
                identity_record["components"] = binder_components
                designed_target_components = [
                    component for component in record.get("components") or []
                    if component.get("role") == "target" and component.get("type") == "protein"
                    and component.get("sequence")
                ]
                if not designed_target_components:
                    designed_target_components = _target_sequence_components(
                        candidate=candidate, inputs=inputs, metadata=metadata, run_dir=run_dir,
                    )
                provenance.update(
                    binder_chains=[str(component["chain"]) for component in binder_components],
                    target_chains=[str(component["chain"]) for component in designed_target_components],
                )

        # An explicit owner/admin removal applies only to this run's group.
        # Preserve its Results candidate and memberships in all other runs.
        if ("mol-" + _canonical_fingerprint(identity_record)[:32]) in removed_ids:
            continue

        if designed_target_components:
            # Preserve a molecule relationship as well as the complex coordinates.
            # This is the exact modeled target sequence, which may be a crop of
            # the campaign input; it is not a native structural reference.
            target_record = {"name": f"{name} · modeled target", "components": designed_target_components,
                             "files": [], "metadata": {"target_identity_type": "sequence_only",
                                                         "structural_reference": False}}
            target_id, _, _ = _find_or_create_input(
                database, user, project, name, source_for_molecule, target_record,
                provenance={**provenance, "input_role": "modeled_target"}, source_key_suffix=":target",
            )
            relationship = {
                "design_target_molecule_id": target_id,
                "design_target_relationship": "exact modeled target sequence",
                "design_target_sequence_sha256": hashlib.sha256(
                ":".join(sorted(component["sequence"] for component in designed_target_components)).encode()
                ).hexdigest(),
                "design_target_is_structural_reference": bool(
                    source_target_sequences
                    and sorted(sequence for sequence in source_target_sequences)
                    == sorted(component["sequence"] for component in designed_target_components)
                ),
            }
            provenance.update(relationship)
            identity_record["metadata"] = {**identity_record.get("metadata", {}), **relationship}

        if database_subject is not None:
            molecule_id, destination = database_subject["molecule_id"], database_subject["group_id"]
            # Evaluations authored on this workstation belong to its local
            # Evaluation collection, even when the selected molecule and
            # representation are owned by a remote peer. Reuse the exact
            # canonical identity; never write an observation/source row into
            # the source owner's replicated group.
            destination = output_group("evaluation")
            database.add_to_group(user, destination, [molecule_id])
            provenance.update(
                source_kind="database",
                source_molecule_id=molecule_id,
                source_structure_id=database_subject["structure_id"],
                source_group_id=database_subject["group_id"],
                source_structure_digest=database_subject["structure_digest"],
                source_owner_id=database_subject.get("source_owner_id") or "",
                source_node=database_subject.get("source_node") or "",
                source_structure_owner_node=database_subject.get("source_structure_owner_node") or "",
                source_structure_owner_group_id=database_subject.get("source_structure_owner_group_id") or "",
                source_structure_owner_structure_id=database_subject.get("source_structure_owner_structure_id") or "",
                source_reference_kind=database_subject.get("source_reference_kind") or "structure",
                source_sequence_sha256=database_subject.get("source_sequence_sha256") or "",
            )
            _put_source(database, destination, molecule_id, source_for_molecule, provenance)
        elif mode == "curated":
            destination = output_group("curated")
            with database.connection() as connection:
                molecule_id = database._put(connection, destination, identity_record, source_key=source_for_molecule, origin=provenance, stamp=stamp)
        elif not rejected and not is_backbone and (mode == "generated" or (mode == "auto" and stage in {"sequence_design", "generation.backbone_sequence"})):
            destination = output_group("generated")
            with database.connection() as connection:
                molecule_id = database._put(connection, destination, identity_record, source_key=source_for_molecule, origin=provenance, stamp=stamp)
            if mode == "generated":
                generated_molecule_ids.add(molecule_id)
        elif not rejected and not is_backbone:
            lookup_record = identity_record if ligand_analysis or protein_analysis or prediction_ligand else record
            molecule_id, destination, matched = _find_original_group(database, user["id"], project, lookup_record, candidate)
            if not matched:
                if mode == "evaluation" or stage in {"complex_refolding", "analysis", "benchmark", "monomer_refolding"}:
                    destination = output_group("evaluation")
                    with database.connection() as connection:
                        molecule_id = database._put(connection, destination, identity_record, source_key=source_for_molecule,
                                                    origin=provenance, stamp=stamp)
                else:
                    destination = output_group("intermediate")
                    with database.connection() as connection:
                        molecule_id = database._put(connection, destination, identity_record, source_key=source_for_molecule,
                                                    origin=provenance, stamp=stamp)
            else:
                _put_source(database, destination, molecule_id, source_for_molecule, provenance)
        else:
            destination = output_group("intermediate")
            with database.connection() as connection:
                molecule_id = database._put(connection, destination, identity_record, source_key=source_for_molecule, origin=provenance, stamp=stamp)

        _merge_member_roles(database, destination, molecule_id, record.get("components"))

        if previous_source and (previous_source["molecule_id"], previous_source["group_id"]) != (molecule_id, destination):
            stale_evaluation_memberships.add((previous_source["group_id"], previous_source["molecule_id"]))

        attach_candidate_structure = (database_subject is None) and (not prediction_run or id(candidate) in best_prediction_ids)
        candidate_files = (record.get("files") or []) if attach_candidate_structure else []
        candidate_structure_ids = {}
        for item in candidate_files:
            selected_item = dict(item)
            if prediction_run:
                selected_item["database_structure_policy"] = "best_prediction"
                label = f"{_prediction_engine_label(tool)} best complex"
                selected_item["label"] = label
                # Include the label in the structure identity so old generic
                # "Complex Pdb" records are replaced during re-indexing.
                selected_item["database_structure_label"] = label
            structure_id = database.attach_structure(destination, molecule_id, selected_item,
                                                    context=provenance, created_at=stamp)
            candidate_structure_ids[str(item.get("source_file_field") or "")] = structure_id
            # These are the representative structures users expect to open
            # from any workstation. Keep a content-addressed copy in the
            # shared Database store before peer replication captures the row.
            database.make_structure_portable(user, structure_id)
            if prediction_run:
                preferred_prediction_structure_ids.add(structure_id)

        candidate_groups[candidate_id] = (molecule_id, destination)
        seen_molecules.add(molecule_id)
        metrics = _numeric(candidate.get("metrics") or {})
        if metrics or (candidate.get("normalized_output") or {}).get("metrics") or (candidate.get("metrics") or {}).get("qc_status"):
            # Include the actual producing adapter, rather than only the run's
            # parent engine (a generator may contain several sequence methods).
            raw = candidate.get("raw_metadata") or {}
            candidate_context = {key: candidate.get(key) for key in ("candidate_id", "source_tool", "stage", "binder_sequence", "binder_chains", "target_chains", "complex_pdb", "binder_pdb", "target_pdb", "model_id", "model_index")}
            # String QC markers are deliberately excluded from native numeric
            # observation metrics, but shared output measures may depend on
            # them for provenance (for example, heavy-atom geometry status).
            if (candidate.get("metrics") or {}).get("binder_geometry_status") is not None:
                candidate_context["binder_geometry_status"] = candidate["metrics"]["binder_geometry_status"]
            candidate_context["raw_metadata"] = {key: raw[key] for key in (
                "fasta_header", "engine_version", "operation", "structure_kind", "image_digest",
                "validation_acceptance", "validation_measurements", "validation_provenance",
                "database_source", "complex_analysis_selection", "reference_placement",
                "reference_environment", "tag_fusion_suitability", "sequence_liability",
                "measurement_provenance", "analysis_parameters",
            ) if key in raw}
            settings = inputs.get("params") or {}
            from mn_molecular_design_studio.evaluations import file_digest
            target_source = _path(candidate.get("target_pdb"), run_dir)
            target_components = [component for component in record.get("components", []) if component.get("role") == "target"]
            evaluation_target = {"target_id": _id("target-", _json(target_components))} if target_components else {}
            if provenance.get("design_target_molecule_id"):
                evaluation_target["target_id"] = provenance["design_target_molecule_id"]
            elif database_subject is not None and database_subject.get("target_molecule_id"):
                evaluation_target["target_id"] = database_subject["target_molecule_id"]
            analysis_source_context = ({
                "source_kind": "database" if database_source is not None else "results",
                "source_molecule_id": database_subject["molecule_id"],
                "source_structure_id": database_subject["structure_id"],
                "source_group_id": database_subject["group_id"],
                "source_structure_digest": database_subject["structure_digest"],
                "representation_id": database_subject["structure_id"],
                "structure_digest": _candidate_calculation_digest(record) or database_subject["structure_digest"],
                "calculation_structure_digest": _candidate_calculation_digest(record),
                "source_reference_kind": database_subject.get("source_reference_kind") or "structure",
                "source_sequence_sha256": database_subject.get("source_sequence_sha256") or "",
                "is_representative": bool(database_subject["structure_id"]),
            } if database_subject is not None else {})
            if database_subject is None and candidate_structure_ids:
                selected_structure_id = next((candidate_structure_ids[key] for key in (
                    "complex_pdb", "binder_pdb", "validation_complex_pdb", "target_pdb")
                    if candidate_structure_ids.get(key)), "")
                if selected_structure_id:
                    with database.connection() as connection:
                        selected_row = connection.execute(
                            "SELECT context FROM molecule_structures WHERE id=?", (selected_structure_id,)
                        ).fetchone()
                    selected_context = json.loads(selected_row["context"] or "{}") if selected_row else {}
                    selected_digest = str(selected_context.get("content_sha256") or selected_context.get("portable_sha256") or "")
                    analysis_source_context = {
                        "representation_id": selected_structure_id,
                        "structure_digest": selected_digest,
                        "source_structure_digest": selected_digest,
                        "calculation_structure_digest": selected_digest,
                        "is_representative": True,
                    }
            database.add_observation(destination, molecule_id, tool, run_dir.name, metrics,
                context={**provenance, **evaluation_target, **analysis_source_context,
                         "candidate": candidate_context, "settings": settings,
                         "status": (candidate.get("metrics") or {}).get("qc_status") or "completed",
                         "error": (candidate.get("metrics") or {}).get("qc_error") or "",
                         "image_digest": raw.get("image_digest") or "",
                         "target_structure_digest": raw.get("target_structure_digest") or (file_digest(target_source) if target_source else ""),
                         "structure_digest": analysis_source_context.get("structure_digest") or raw.get("structure_digest") or "",
                         "model_id": str(candidate.get("model_id") or candidate.get("model_index") or
                                         ((raw.get("complex_analysis_selection") or {}).get("model_id") if isinstance(raw.get("complex_analysis_selection"), dict) else "") or ""),
                         "is_representative": not prediction_run or id(candidate) in best_prediction_ids},
                source_key=f"{source_for_molecule}:metrics", created_at=stamp)

    # Re-indexing is additive. A changed or incomplete historical decision
    # must not erase a previously retained candidate; explicit user removal
    # remains represented by the run-group tombstone checked above.

    if prediction_run:
        _prune_stale_run_structures(database, run_dir.parent.name, run_dir.name,
                                    preferred_prediction_structure_ids, candidate_ids=prediction_candidate_ids)
    _prune_empty_evaluation_members(database, stale_evaluation_memberships)

    # Screening inputs stay in Results. Docking rows may still attach to an
    # already imported canonical compound, but never create one just because
    # it appeared in a large campaign input set.
    compound_rows = {}
    for compound in (() if analysis_only else input_compounds):
        record = {"name": compound.get("compound_name") or compound["compound_id"],
                  "components": [{"type": "small_molecule", "chain": "L", "smiles": compound["smiles"]}],
                  "metadata": {**(compound.get("properties") or {}), "compound_id": compound["compound_id"],
                               "source_tool": tool, "stage": "input"}}
        molecule_id, destination, matched = _find_original_group(database, user["id"], project, record)
        remote_match = False
        if not matched:
            molecule_id, source_group_id, remote_match = _find_accessible_group(database, user, project, record)
            if remote_match:
                destination = output_group("evaluation")
                database.add_to_group(user, destination, [molecule_id])
        if not matched:
            if not remote_match:
                continue
        row_source_key = f"{source_key}:input-compound:{compound['compound_id']}"
        _put_source(database, destination, molecule_id, row_source_key,
                    {**origin, "compound_id": compound["compound_id"],
                     "source_group_id": source_group_id if remote_match else destination,
                     "collection_kind": "evaluation" if remote_match else "inputs"})
        compound_rows[compound["compound_id"]] = (molecule_id, destination, record)
        seen_molecules.add(molecule_id)
    # Preserve the explicit target structure parsed before candidate ingestion.
    target_id = None
    target_group_id = None
    target_record = None
    if target_path and not analysis_only:
        try:
            target_record = target_structure_record
            target_id, target_group_id, matched = _find_original_group(database, user["id"], project, target_record)
            target_provenance = {**origin, "input_role": "receptor", "collection_kind": "inputs"}
            if not matched:
                target_group_id = _ensure_group(database, user, _id("inputs-", f"{user['id']}:{project or ''}"),
                    "Tool inputs", project, "tool_inputs",
                    {"type": "job", "owner_id": user["id"], "project_id": project}, "inputs")
                with database.connection() as connection:
                    target_id = database._put(connection, target_group_id, target_record, source_key=source_key + ":target", origin=target_provenance, stamp=stamp)
                with database.connection() as connection:
                    target_structures = connection.execute(
                        "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=?",
                        (target_group_id, target_id),
                    ).fetchall()
                for structure in target_structures:
                    database.make_structure_portable(user, str(structure["id"]))
            else:
                _put_source(database, target_group_id, target_id, source_key + ":target", target_provenance)
                for item in target_record.get("files") or []:
                    structure_id = database.attach_structure(
                        target_group_id, target_id, item, context=target_provenance, created_at=stamp,
                    )
                    database.make_structure_portable(user, structure_id)
                database.update_metadata(user, target_group_id, target_id, {
                    "target_identity_type": "structure_backed", "structural_reference": True,
                })
            seen_molecules.add(target_id)
        except ValueError as exc:
            log.info("Could not index receptor for %s: %s", run_dir, exc)
    elif input_target_record and not analysis_only:
        target_provenance = {**origin, "input_role": "sequence_defined_target",
                             "collection_kind": "inputs", "structural_reference": False}
        target_id, target_group_id, _ = _find_or_create_input(
            database, user, project, name, source_key, input_target_record,
            provenance=target_provenance, source_key_suffix=":target-sequence",
        )
        database.update_metadata(user, target_group_id, target_id, {
            "target_identity_type": "sequence_only", "structural_reference": False,
        })
        seen_molecules.add(target_id)
    params = inputs.get("params") or {}
    context = {**origin, "target_molecule_id": target_id, "target_path": str(target_path) if target_path else None,
               "settings": params, "units": (result.get("metrics") or {}).get("score_units"),
               "primary_ranking_metric": (result.get("metrics") or {}).get("primary_ranking_metric"),
               "ranking_direction": (result.get("metrics") or {}).get("ranking_direction")}
    if target_path and (run_dir / "artifacts/docking_scores.csv").is_file():
        from mn_molecular_design_studio.bound_docking import register_receptor
        try:
            context.update(register_receptor(database, user, target_path, project_id=project,
                name=str((target_record or {}).get("name") or "Docking receptor")))
        except (ValueError, OSError) as exc:
            log.info("Could not retain docking receptor snapshot for %s: %s", run_dir, exc)
    target_context_by_attempt = {}
    if is_legacy_aggregate:
        replicate_root = run_dir / "artifacts" / "replicates"
        if replicate_root.is_dir():
            for replicate_dir in sorted(path for path in replicate_root.iterdir() if path.is_dir()):
                attempt_context = _legacy_attempt_target_context(run_dir, replicate_dir.name, metadata)
                if attempt_context:
                    target_context_by_attempt[replicate_dir.name] = attempt_context
    preferred_docking_structure_ids = set()
    has_docking_scores = False
    for filename in ("docking_scores.csv", "nesso_affinity.csv"):
        path = run_dir / "artifacts" / filename
        if not path.is_file():
            continue
        text = path.read_text(errors="replace")
        reader = csv.DictReader(io.StringIO(text), delimiter="\t" if "\t" in text.partition("\n")[0] else ",")
        rows = list(reader)
        pose_paths = [_path(row.get("pose_path"), run_dir) for row in rows]
        best_pose_indices = _best_docking_pose_indices(rows, pose_paths) if filename == "docking_scores.csv" else set()
        if filename == "docking_scores.csv":
            has_docking_scores = True
        for index, row in enumerate(rows):
            match = compound_rows.get(str(row.get("compound_id") or ""))
            if not match:
                continue
            molecule_id, destination, record = match
            legacy_run_id = str(row.get("legacy_run_id") or "").strip()
            row_target_context = target_context_by_attempt.get(legacy_run_id, {})
            pose = pose_paths[index]
            pose_digest = _file_sha256(pose) if pose and pose.is_file() else ""
            pose_identity = {
                "id": _id("pose-", _json([tool, str(row.get("compound_id") or ""),
                                          legacy_run_id, row.get("attempt_index"), row.get("pose_rank"), pose_digest])),
                "compound_id": str(row.get("compound_id") or ""),
                "pose_rank": str(row.get("pose_rank") or ""),
                "legacy_run_id": legacy_run_id,
                "attempt_index": row.get("attempt_index"),
                "structure_digest": pose_digest,
            }
            if filename == "docking_scores.csv" and index in best_pose_indices and pose:
                attempt_label = f" · attempt {legacy_run_id[:8]}" if legacy_run_id else ""
                structure_id = database.attach_structure(
                    destination, molecule_id,
                    {"path": str(pose), "label": f"{tool} best pose{attempt_label}", "database_structure_policy": "best_docking_pose"},
                    context={**context, **row_target_context, "pose_rank": row.get("pose_rank"), "legacy_run_id": legacy_run_id,
                             "attempt_index": row.get("attempt_index"), "row": row}, created_at=stamp)
                database.make_structure_portable(user, structure_id)
                preferred_docking_structure_ids.add(structure_id)
            database.add_observation(destination, molecule_id, tool, run_dir.name, _numeric({key: value for key, value in row.items() if key != "compound_id"}),
                context={**context, **row_target_context, "row": row, "legacy_run_id": legacy_run_id,
                         "attempt_index": row.get("attempt_index"), "pose_identity": pose_identity,
                         "pose_id": pose_identity["id"], "pose_digest": pose_digest,
                         "structure_digest": pose_digest, "model_id": row.get("model_id") or "",
                         "representation_id": structure_id if filename == "docking_scores.csv" and index in best_pose_indices and pose else "",
                         "is_representative": index in best_pose_indices},
                source_key=f"{source_key}:{filename}:{index}", created_at=stamp)
    if has_docking_scores:
        _prune_stale_run_structures(database, run_dir.parent.name, run_dir.name,
                                    preferred_docking_structure_ids, docking=True)
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
        signature = hashlib.sha256(_json({"migration_version": 3, "payload": payload}).encode()).hexdigest()
        with database.connection() as connection:
            prior = connection.execute("SELECT signature FROM molecular_ingestions WHERE source_key=?", (source_key,)).fetchone()
        if prior and prior[0] == signature:
            continue
        actor = {"id": access["owner_id"], "role": "admin"}
        group_id = _id("legacy-cs-", set_id)
        database.create_group(actor, str(payload.get("name") or "Imported compounds"), group_id=group_id,
                              project_id=access["project_id"], import_method="user_uploads",
                              origin={"type": "legacy", "legacy_compound_set_id": set_id})
        records = [{"name": compound.get("compound_name") or compound["compound_id"],
                    "components": [{"type": "small_molecule", "chain": "L", "smiles": compound["smiles"]}],
                    "metadata": {**(compound.get("properties") or {}), "compound_id": compound["compound_id"]}}
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


def sync_workspace_async(user, *, force=False, database=None):
    """Refresh rebuildable job indexes without blocking the saved Library view."""
    database = database or MolecularDatabase()
    key = (str(database.path), str(user["id"]))
    with _SYNC_WORKERS_LOCK:
        if key in _SYNC_WORKERS or (not force and time.monotonic() - _LAST_SYNC.get(key, 0) < 30):
            return list(_LAST_SYNC_ERRORS.get(key, []))
        _SYNC_WORKERS.add(key)
    actor = dict(user)

    def refresh():
        try:
            sync_workspace(actor, force=force, database=database)
        except Exception as exc:
            log.exception("Molecular workspace refresh failed for %s", actor["id"])
            _LAST_SYNC_ERRORS[key] = [{"run_id": "workspace-index", "error": str(exc)}]
        finally:
            with _SYNC_WORKERS_LOCK:
                _LAST_SYNC[key] = time.monotonic()
                _SYNC_WORKERS.discard(key)

    try:
        threading.Thread(target=refresh, name="studio-molecular-index", daemon=True).start()
    except RuntimeError:
        with _SYNC_WORKERS_LOCK:
            _SYNC_WORKERS.discard(key)
        log.exception("Could not start molecular workspace refresh")
    return list(_LAST_SYNC_ERRORS.get(key, []))


def workspace_sync_running(user, *, database):
    with _SYNC_WORKERS_LOCK:
        return (str(database.path), str(user["id"])) in _SYNC_WORKERS


def record_completed_run(run_dir):
    """Indexing cannot change the scientific job's completion status."""
    try:
        from mn_molecular_design_studio.workspace import is_studio_job
        run_dir = Path(run_dir).resolve()
        if not is_studio_job(run_dir):
            return
        from mn_molecular_design_studio.collection_policy import automatic_collection_allowed, policy_for_run
        from mn_molecular_design_studio.protein_library import record_run_inputs

        if policy_for_run(run_dir) != "none":
            # Source capture is independent of output promotion. A user's
            # Protein Library input remains useful even when the prediction
            # fails or the run's output policy is Results-only.
            record_run_inputs(run_dir)
        if not automatic_collection_allowed(run_dir):
            return
        metadata = _read(run_dir / "metadata.json")
        result = _read(run_dir / "result.json")
        if not _is_completed_success(metadata, result):
            return
        # Completion is the authoritative indexing event. Bypass the change
        # signature so a prior queued-run scan cannot suppress final outputs.
        for attempt in range(3):
            try:
                ingest_run(run_dir, force=True)
                return
            except Exception:
                if attempt == 2:
                    log.exception("Molecular Database could not index completed run %s; refresh will retry", run_dir)
                    return
                log.warning("Molecular Database indexing attempt %d failed for %s; retrying", attempt + 1, run_dir, exc_info=True)
                time.sleep(0.25 * (2 ** attempt))
    except Exception:
        log.exception("Molecular Database could not index completed run %s; refresh will retry", run_dir)
