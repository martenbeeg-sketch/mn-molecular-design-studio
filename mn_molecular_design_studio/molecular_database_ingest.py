"""Idempotent adapters for Studio's existing molecular artifact contracts."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import logging
import re
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
    # A candidate has one molecular identity. Do not attach target-only or
    # binder-only files as if they represented a different complex's identity.
    for key in ("complex_pdb", "binder_pdb", "target_pdb"):
        if key == "target_pdb" and sequence:
            # A target-only file cannot identify the newly designed sequence.
            continue
        path = _path(candidate.get(key), run_dir)
        if not path:
            continue
        try:
            components = imports.parse_structure(path.read_bytes(), path.name)[0]["components"]
        except ValueError:
            continue
        files = [{"path": str(path), "label": key.replace("_", " ").title(), "dimensions": 3}]
        break
    raw_components = (candidate.get("raw_metadata") or {}).get("components")
    if components is None and isinstance(raw_components, list):
        try:
            components = imports.canonical_components(raw_components)[0]
        except (KeyError, TypeError, ValueError):
            components = None
    raw_metadata = candidate.get("raw_metadata") or {}
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
            component.pop("ccd", None)
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
                         **raw_metadata, "stage": candidate.get("stage"), "source_tool": candidate.get("source_tool"),
                         "output_kind": (candidate.get("normalized_output") or {}).get("kind")}}


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
    status = str(candidate.get("status") or raw.get("status") or raw.get("filter_status") or metrics.get("status") or "").strip().casefold()
    if status in {"rejected", "failed", "invalid", "filtered", "discarded", "not accepted", "not_accepted"}:
        return True
    return (candidate.get("accepted") is False or raw.get("accepted") is False
            or metrics.get("pass_filters") is False or metrics.get("native_pass_filters") is False)


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


def _prune_generated_group(database, group_id, keep_ids):
    """Remove superseded or rejected entries from one run's generated group."""
    keep_ids = set(keep_ids)
    with database.connection() as connection:
        members = connection.execute(
            "SELECT molecule_id FROM molecule_members WHERE group_id=?", (group_id,)
        ).fetchall()
        for member in members:
            molecule_id = member["molecule_id"]
            if molecule_id in keep_ids:
                continue
            for table in ("molecule_sources", "molecule_structures", "molecule_observations"):
                connection.execute(f"DELETE FROM {table} WHERE group_id=? AND molecule_id=?", (group_id, molecule_id))
            connection.execute("DELETE FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id))


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
        for field in ("complex_pdb", "binder_pdb", "target_pdb"):
            reference = candidate.get(field)
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
    ingestion_version = 8 if is_legacy_aggregate else (7 if needs_compound_identity_refresh else 6)
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
    if prediction_run and candidate_records:
        from mn_molecular_design_studio.cofolding_results import best_candidates_by_complex
        with_structures = [candidate for candidate, record in candidate_records if record.get("files")]
        best_prediction_ids = {id(candidate) for candidate in best_candidates_by_complex(with_structures)}

    input_set_path = run_dir / "artifacts/input/compound_set.json"
    input_compounds = _run_input_compounds(run_dir)
    generated_set_path = run_dir / "artifacts/generated_compound_set.json"
    score_paths = [run_dir / "artifacts" / filename for filename in ("docking_scores.csv", "nesso_affinity.csv")]
    md_metrics_supported = "md" in run_dir.parent.name and bool(result.get("metrics"))
    if not candidate_records and not generated_set_path.is_file() and not (input_compounds and any(path.is_file() for path in score_paths)) and not md_metrics_supported:
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
    for candidate, record in candidate_records:
        stage = str(candidate.get("stage") or "")
        sequence = candidate.get("binder_sequence") or candidate.get("sequence")
        is_backbone = stage == "generation.backbone" and not sequence
        rejected = _candidate_rejected(candidate)
        raw_candidate = candidate.get("raw_metadata") or {}
        preserve_failed_native_refold = (
            mode == "generated" and rejected and stage == "complex_refolding"
            and raw_candidate.get("result_kind") == "native_pipeline"
            and raw_candidate.get("candidate_pool_level") == "native_pipeline"
        )
        # A failed generated design remains visible in its job results and
        # native tables, but it is not a reusable library candidate. Native
        # pipeline refolds with failed tool filters remain as intermediate
        # observations so their structures and stage metrics are auditable.
        if mode == "generated" and rejected and not preserve_failed_native_refold:
            continue
        provenance = {**origin, "candidate_id": candidate.get("candidate_id"), "parents": candidate.get("parents") or []}
        candidate_id = str(candidate.get("candidate_id") or _canonical_fingerprint(record))
        source_for_molecule = f"{source_key}:{candidate_id}"
        with database.connection() as connection:
            previous_source = connection.execute(
                "SELECT molecule_id,group_id FROM molecule_sources WHERE source_key=?",
                (source_for_molecule,),
            ).fetchone()
        prediction_candidate_ids.add(candidate_id)
        identity_record = {**record, "files": []}
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

        if mode == "curated":
            destination = output_group("curated")
            with database.connection() as connection:
                molecule_id = database._put(connection, destination, identity_record, source_key=source_for_molecule, origin=provenance, stamp=stamp)
        elif preserve_failed_native_refold:
            destination = output_group("intermediate")
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

        attach_candidate_structure = not prediction_run or id(candidate) in best_prediction_ids
        candidate_files = (record.get("files") or []) if attach_candidate_structure else []
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
            candidate_context = {key: candidate.get(key) for key in ("candidate_id", "source_tool", "stage", "binder_sequence", "binder_chains", "target_chains", "complex_pdb", "binder_pdb", "target_pdb")}
            # String QC markers are deliberately excluded from native numeric
            # observation metrics, but shared output measures may depend on
            # them for provenance (for example, heavy-atom geometry status).
            if (candidate.get("metrics") or {}).get("binder_geometry_status") is not None:
                candidate_context["binder_geometry_status"] = candidate["metrics"]["binder_geometry_status"]
            candidate_context["raw_metadata"] = {key: raw[key] for key in ("fasta_header", "engine_version", "operation", "structure_kind", "image_digest") if key in raw}
            settings = inputs.get("params") or {}
            from mn_molecular_design_studio.evaluations import file_digest
            target_source = _path(candidate.get("target_pdb"), run_dir)
            target_components = [component for component in record.get("components", []) if component.get("role") == "target"]
            evaluation_target = {"target_id": _id("target-", _json(target_components))} if target_components else {}
            database.add_observation(destination, molecule_id, tool, run_dir.name, metrics,
                context={**provenance, **evaluation_target, "candidate": candidate_context, "settings": settings,
                         "status": (candidate.get("metrics") or {}).get("qc_status") or "completed",
                         "error": (candidate.get("metrics") or {}).get("qc_error") or "",
                         "image_digest": raw.get("image_digest") or "",
                         "target_structure_digest": raw.get("target_structure_digest") or (file_digest(target_source) if target_source else ""),
                         "structure_digest": raw.get("structure_digest") or "",
                         "is_representative": not prediction_run or id(candidate) in best_prediction_ids},
                source_key=f"{source_for_molecule}:metrics", created_at=stamp)

    if mode == "generated":
        _prune_generated_group(database, _id("job-", source_key), generated_molecule_ids)

    if prediction_run:
        _prune_stale_run_structures(database, run_dir.parent.name, run_dir.name,
                                    preferred_prediction_structure_ids, candidate_ids=prediction_candidate_ids)
    _prune_empty_evaluation_members(database, stale_evaluation_memberships)

    # The ligand contracts preserve exactly the consumed and generated sets.
    compound_rows = {}
    set_payloads = [(input_compounds, False)]
    if generated_set_path.is_file():
        from mn_molecular_design_studio.ligand_compounds import load_compound_set
        try:
            set_payloads.append((list(load_compound_set(generated_set_path).get("compounds") or []), True))
        except (OSError, ValueError, KeyError, TypeError):
            log.info("Could not read generated compound set %s", generated_set_path)
    for compounds, generated in set_payloads:
        if not compounds:
            continue
        destination = output_group("generated") if generated else None
        for index, compound in enumerate(compounds):
            record = {"name": compound.get("compound_name") or compound["compound_id"], "components": [{"type": "small_molecule", "chain": "L", "smiles": compound["smiles"]}],
                      "metadata": {**(compound.get("properties") or {}), "compound_id": compound["compound_id"],
                                   "source_tool": tool, "stage": "generated" if generated else "input"}}
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
                context={**context, **row_target_context, "row": row, "legacy_run_id": row.get("legacy_run_id"),
                         "attempt_index": row.get("attempt_index"), "is_representative": index in best_pose_indices},
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


def record_completed_run(run_dir):
    """Indexing cannot change the scientific job's completion status."""
    try:
        from mn_molecular_design_studio.workspace import is_studio_job
        run_dir = Path(run_dir).resolve()
        if not is_studio_job(run_dir):
            return
        from mn_molecular_design_studio.collection_policy import automatic_collection_allowed
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
