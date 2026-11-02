"""Authorized adapters from saved structures into existing analysis workflows.

Database structures are represented as analysis inputs, not prediction jobs. The
adapter makes a small, user-owned source bundle consumed by Studio's registered
workflow functions and scheduler; it does not create a Results job or scientific
measurements of its own.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.evaluations import evaluation_descriptor, file_digest


_WORKFLOWS = {
    "analysis.biotite_interface",
    "analysis.ipsae",
    "analysis.rosetta",
    "analysis.pep_patch",
    "analysis.interface_clustering",
    "analysis.foldseek",
    "analysis.plip",
    "analysis.pandamap",
    "analysis.posebusters",
}
_STRUCTURE_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif"}
_SELECTION_FIELDS = {
    "molecule_id", "structure_id", "group_id", "evaluation_id", "model_id",
    "binder_chains", "target_chains", "ligand_selection",
}
_LIGAND_WORKFLOWS = {"analysis.plip", "analysis.pandamap", "analysis.posebusters"}


def _string(value: Any, field: str) -> str:
    text = str(value or "").strip()
    if not text or len(text) > 512:
        raise ValueError(f"Database analysis source requires a valid {field}.")
    return text


def _chain_ids(components: list[dict]) -> tuple[list[str], list[str]]:
    proteins, all_chains = set(), set()
    for component in components:
        if not isinstance(component, dict):
            continue
        residue_key = str(component.get("residue_key") or "")
        chain = str(component.get("auth_asym_id") or component.get("chain") or "").strip()
        if residue_key and ":" in residue_key:
            chain = residue_key.split(":", 1)[0].strip()
        if not chain:
            continue
        all_chains.add(chain)
        if str(component.get("type") or "").casefold() == "protein":
            proteins.add(chain)
    return sorted(proteins), sorted(all_chains)


def authorize_database_structure(
    user: dict,
    selection: dict,
    *,
    project_id: str | None,
    database=None,
    workflow_id: str | None = None,
) -> dict[str, Any]:
    """Resolve stable Database IDs and validate read access before staging paths."""
    if workflow_id and workflow_id not in _WORKFLOWS:
        raise ValueError("Database structures are not enabled for this analysis workflow.")
    if not isinstance(selection, dict) or set(selection) - _SELECTION_FIELDS:
        raise ValueError("Database analysis source contains unsupported fields.")
    if not str(user.get("id") or ""):
        raise PermissionError("A signed-in user is required to analyze a Database structure.")

    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    database = database or MolecularDatabase()
    molecule_id = _string(selection.get("molecule_id"), "record ID")
    structure_id = _string(selection.get("structure_id"), "structure ID")
    group_id = _string(selection.get("group_id"), "source group ID")
    rows = database.list_molecules(user, molecule_ids=[molecule_id], collection="all", limit=1)["items"]
    if not rows:
        raise PermissionError("This Database record is unavailable or you do not have access to it.")
    record = rows[0]
    structure = next((item for item in record.get("representations", [])
                      if str(item.get("id") or "") == structure_id), None)
    if structure is None or str(structure.get("group_id") or "") != group_id:
        raise PermissionError("The selected Database structure is not attached to this accessible record.")
    group = database.group(group_id, user)
    source_project_id = str(group.get("project_id") or "")
    requested_project = str(project_id or "")
    if source_project_id != requested_project:
        raise PermissionError(
            "Analysis must remain in the source Database group's project; private Database records require a private analysis job."
        )

    evaluation_id = str(selection.get("evaluation_id") or "")
    selected_evaluation = None
    if evaluation_id:
        known = {}
        for observation in record.get("observations", []):
            try:
                known[evaluation_descriptor(observation)["id"]] = observation
            except (KeyError, TypeError, ValueError):
                continue
        if evaluation_id not in known:
            raise PermissionError("The selected Database evaluation is unavailable or inaccessible.")
        selected_evaluation = known[evaluation_id]

    context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
    saved_model = str(context.get("model_id") or context.get("model_index") or context.get("model") or "")
    model_id = str(selection.get("model_id") or saved_model or structure_id)
    if saved_model and model_id != saved_model:
        raise ValueError("The selected model does not match this Database structure representation.")
    path = database.structure_file(user, structure_id, materialize=True).resolve()
    if path.suffix.casefold() not in _STRUCTURE_SUFFIXES:
        raise ValueError("This analysis workflow requires a PDB or mmCIF Database structure.")
    digest = file_digest(path)
    if not digest:
        raise FileNotFoundError("The selected Database structure could not be read.")
    if selected_evaluation:
        context = selected_evaluation.get("context") if isinstance(selected_evaluation.get("context"), dict) else {}
        selected_representation = str(context.get("representation_id") or "")
        # Evaluation outputs may describe a normalized/protonated analysis
        # representation in structure_digest. When the evaluator records its
        # input lineage, compare that original digest to the selected saved
        # Database representation instead. The representation ID check above
        # remains mandatory, so this does not allow an evaluation from another
        # model/pose to authorize an analysis of the selected structure.
        selected_digest = str(context.get("source_structure_digest") or context.get("structure_digest") or "")
        if selected_representation and selected_representation != structure_id:
            raise ValueError("Choose a structure associated with the selected Database evaluation.")
        if selected_digest and selected_digest != digest:
            raise ValueError("The selected evaluation was calculated for a different structure representation.")

    from mn_molecular_design_studio.molecule_imports import parse_structure

    parsed = parse_structure(path.read_bytes(), path.name, include_ligands=True)
    components = list(parsed[0].get("components") or []) if parsed else []
    protein_chains, _ = _chain_ids(components)
    ligand_selection = selection.get("ligand_selection")
    if workflow_id in _LIGAND_WORKFLOWS:
        if not isinstance(ligand_selection, dict) or set(ligand_selection) - {"chain", "residue_key", "ccd", "smiles"}:
            raise ValueError("Select one protein-bound ligand with a supported identity.")
        selected_residue = str(ligand_selection.get("residue_key") or "").strip()
        selected_chain = str(ligand_selection.get("chain") or "").strip()
        selected_ccd = str(ligand_selection.get("ccd") or "").strip().upper()
        if not selected_residue or not selected_chain:
            raise ValueError("A ligand selection must include its source chain and residue key.")
        parsed_ligands = [item for item in components if str(item.get("type") or "").casefold() == "small_molecule"]
        matches = [item for item in parsed_ligands
                   if str(item.get("residue_key") or "") == selected_residue
                   and str(item.get("residue_key") or "").split(":", 1)[0] == selected_chain
                   and (not selected_ccd or str(item.get("ccd") or "").upper() == selected_ccd)]
        if len(matches) != 1:
            raise ValueError("The selected ligand residue is not unique in this saved structure.")
        matching_record = next((item for item in (record.get("components") or [])
            if str(item.get("type") or "").casefold() == "small_molecule"
            and str(item.get("residue_key") or "") == selected_residue), {})
        stored_smiles = str(matching_record.get("smiles") or "").strip()
        requested_smiles = str(ligand_selection.get("smiles") or "").strip()
        if not requested_smiles and stored_smiles:
            requested_smiles = stored_smiles
        if not str(ligand_selection.get("smiles") or "").strip():
            if not requested_smiles:
                raise ValueError("The saved ligand has no isomeric SMILES identity. Add its verified source identity before running PLIP, PandaMap, or PoseBusters.")
        from rdkit import Chem
        molecule = Chem.MolFromSmiles(requested_smiles)
        if molecule is None:
            raise ValueError("The selected ligand SMILES is invalid.")
        if stored_smiles:
            stored_molecule = Chem.MolFromSmiles(stored_smiles)
            if (stored_molecule is None or Chem.MolToSmiles(molecule, isomericSmiles=True)
                    != Chem.MolToSmiles(stored_molecule, isomericSmiles=True)):
                raise ValueError("The selected ligand identity does not match the saved Database component.")
        ligand_selection = {"chain": selected_chain, "residue_key": selected_residue,
                            "ccd": selected_ccd or str(matches[0].get("ccd") or "").upper(),
                            "smiles": Chem.MolToSmiles(molecule, isomericSmiles=True)}
    elif ligand_selection is not None:
        raise ValueError("Ligand selections are accepted only by PLIP, PandaMap, and PoseBusters.")
    binder_chains = [str(value).strip() for value in (selection.get("binder_chains") or []) if str(value).strip()]
    target_chains = [str(value).strip() for value in (selection.get("target_chains") or []) if str(value).strip()]
    for field, values in (("binder_chains", binder_chains), ("target_chains", target_chains)):
        if len(values) != len(set(values)) or any(value not in protein_chains for value in values):
            raise ValueError(f"{field} must contain unique protein chains in the selected structure.")
    if set(binder_chains) & set(target_chains):
        raise ValueError("The selected protein chain groups must not overlap.")
    if workflow_id == "analysis.pep_patch":
        if not binder_chains:
            raise ValueError("PEP-Patch requires at least one selected protein surface chain.")
    elif workflow_id == "analysis.foldseek":
        if len(binder_chains) != 1:
            raise ValueError("Foldseek requires exactly one selected protein chain per structure.")
    elif workflow_id in _LIGAND_WORKFLOWS:
        if binder_chains or target_chains:
            raise ValueError("Protein-ligand workflows use a selected ligand residue and do not accept protein interface chain groups.")
    elif not binder_chains or not target_chains:
        raise ValueError("Select nonempty, disjoint protein chain groups for this analysis.")
    if workflow_id in {"analysis.rosetta", "analysis.ipsae"} and (
        len(binder_chains) != 1 or len(target_chains) != 1
    ):
        raise ValueError("The registered workflow currently requires one selected protein chain in each group.")

    return {
        "kind": "database",
        "molecule_id": molecule_id,
        "structure_id": structure_id,
        "group_id": group_id,
        "project_id": source_project_id,
        "evaluation_id": evaluation_id,
        "model_id": model_id,
        "structure_digest": digest,
        "structure_path": str(path),
        "structure_label": str(structure.get("label") or path.name),
        "structure_format": str(structure.get("format") or path.suffix.lstrip(".")),
        "binder_chains": binder_chains,
        "target_chains": target_chains,
        "ligand_selection": ligand_selection or {},
        "record_name": str(record.get("name") or molecule_id),
        "record_components": record.get("components") or components,
        "parsed_components": components,
        "structure_context": context,
        "group_origin": group.get("origin") or {},
    }


def build_database_analysis_source(
    user: dict,
    selection: dict,
    *,
    workflow_id: str,
    project_id: str | None,
    database=None,
) -> dict[str, Any]:
    """Create a deterministic structure-input bundle for a registered evaluator."""
    source = authorize_database_structure(
        user, selection, project_id=project_id, database=database, workflow_id=workflow_id,
    )
    identity = {key: source[key] for key in (
        "kind", "molecule_id", "structure_id", "group_id", "project_id", "evaluation_id",
        "model_id", "structure_digest", "binder_chains", "target_chains", "ligand_selection",
    )}
    bundle_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    bundle_dir = user_input_root(str(user["id"])) / "complex_analysis_sources" / bundle_id
    bundle_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        bundle_dir.chmod(0o700)
    except OSError:
        pass
    candidate_id = f"database-{hashlib.sha256((source['molecule_id'] + ':' + source['structure_id']).encode()).hexdigest()[:24]}"
    source_path = Path(source["structure_path"])
    staged_structure = bundle_dir / f"selected_structure{source_path.suffix.casefold()}"
    if not staged_structure.is_file() or file_digest(staged_structure) != source["structure_digest"]:
        shutil.copyfile(source_path, staged_structure)
    staged_structure = staged_structure.resolve()
    selection_context = {
        "candidate_id": candidate_id,
        "structure_path": str(staged_structure),
        "structure_digest": source["structure_digest"],
        "model_id": source["model_id"],
        "binder_chains": source["binder_chains"],
        "target_chains": source["target_chains"],
        "ligand_selection": source.get("ligand_selection") or {},
    }
    source_manifest = {
        "schema": "studio.complex_analysis_database_source.v1",
        **identity,
        "candidate_id": candidate_id,
        "record_name": source["record_name"],
        "structure_path": str(staged_structure),
        "structure_label": source["structure_label"],
        "structure_format": source["structure_format"],
        "ligand_selection": source.get("ligand_selection") or {},
    }
    (bundle_dir / "complex_analysis_source.json").write_text(
        json.dumps(source_manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8",
    )
    (bundle_dir / "metadata.json").write_text(json.dumps({
        "schema": "studio.complex_analysis_source_bundle.v1",
        "source_kind": "database",
        "source_bundle_id": bundle_id,
        "owner_id": str(user["id"]),
        "studio_owner_id": str(user["id"]),
        "project_id": source["project_id"],
        "source_molecule_id": source["molecule_id"],
        "source_structure_id": source["structure_id"],
        "source_group_id": source["group_id"],
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    (bundle_dir / "molecular_design_studio.json").write_text(json.dumps({
        "application": "mn-molecular-design-studio",
        "studio_owner_id": str(user["id"]),
        "project_id": source["project_id"],
        "source_kind": "database",
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    raw_metadata = {
        "analysis_subject": "complex_structure",
        "components": source["record_components"],
        "database_source": source_manifest,
        "complex_analysis_selection": selection_context,
        "model_index": source["model_id"],
    }
    candidate = {
        "candidate_id": candidate_id,
        "source_molecule_id": source["molecule_id"],
        "stage": "database_structure",
        "source_tool": "molecular_database",
        "tool": "molecular_database",
        "complex_pdb": str(staged_structure),
        "binder_chains": source["binder_chains"],
        "target_chains": source["target_chains"],
        "model_index": source["model_id"],
        "raw_metadata": raw_metadata,
        "metrics": {},
    }
    candidate_path = bundle_dir / "analysis_subjects.jsonl"
    candidate_path.write_text(json.dumps(candidate, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "source": source,
        "source_manifest": source_manifest,
        "candidate_id": candidate_id,
        "bundle_id": bundle_id,
        "source_run_dir": str(bundle_dir.resolve()),
        "candidates_jsonl": str(candidate_path.resolve()),
        "selected_structure_path": str(staged_structure),
        "selected_model_id": source["model_id"],
    }


def _collection_structure_suffix(path: Path) -> str:
    name = path.name.casefold()
    if name.endswith(".pdb.gz"):
        return ".pdb.gz"
    if name.endswith(".cif.gz"):
        return ".cif.gz"
    if name.endswith(".mmcif.gz"):
        return ".mmcif.gz"
    if name.endswith(".ent.gz"):
        return ".pdb.gz"
    suffix = path.suffix.casefold()
    return ".pdb" if suffix == ".ent" else suffix


def _collection_structure_bytes(path: Path) -> bytes:
    return path.read_bytes()


def saved_structure_collection_selection_key(selection: dict, project_id: str | None) -> str:
    """Stable lookup key shared by Compare and saved-collection adapters."""
    sources = selection.get("sources") if isinstance(selection, dict) else None
    canonical = []
    for source in sources if isinstance(sources, list) else []:
        if not isinstance(source, dict):
            continue
        item = {key: value for key, value in source.items() if key != "label"}
        if isinstance(item.get("evaluation_ids"), list):
            item["evaluation_ids"] = sorted({str(value) for value in item["evaluation_ids"] if str(value)})
        canonical.append(item)
    payload = {
        "schema": str(selection.get("schema") or ""),
        "project_id": str(project_id or ""),
        "sources": canonical,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def interface_clustering_selection_key(selection: dict, project_id: str | None) -> str:
    """Backward-compatible name for Interface Clustering selection keys."""
    return saved_structure_collection_selection_key(selection, project_id)


def _collection_evaluation_ids(user: dict, source: dict, *, database=None) -> list[str]:
    """Recheck Database evaluation IDs and their association with the structure."""
    requested = sorted({str(value).strip() for value in source.get("evaluation_ids", []) if str(value).strip()})
    if not requested:
        return []
    from mn_molecular_design_studio.complex_analysis_ui import _evaluation_matches_structure, _evaluation_options

    rows = database.list_molecules(user, molecule_ids=[source["molecule_id"]], collection="all", limit=1)["items"]
    if not rows:
        raise PermissionError("The selected Database structure or its evaluations are no longer accessible.")
    record = rows[0]
    structure = next((item for item in record.get("representations", [])
                      if str(item.get("id") or "") == str(source["structure_id"])), None)
    if structure is None:
        raise PermissionError("The selected Database structure is no longer attached to this record.")
    _, by_id = _evaluation_options(record)
    for evaluation_id in requested:
        pair = by_id.get(evaluation_id)
        if pair is None or not _evaluation_matches_structure(pair[1], structure):
            raise PermissionError("A selected Database evaluation is unavailable or belongs to another structure.")
    return requested


def build_saved_structure_collection(
    user: dict,
    selection: dict,
    *,
    workflow_id: str,
    project_id: str | None,
    database=None,
) -> dict[str, Any]:
    """Authorize and stage saved structures for a registered cohort workflow.

    This adapter assembles existing structures and provenance only. Scientific
    calculations remain in the registered workflow. Interface Clustering has
    a shared-target requirement; Foldseek's cohort adapter only requires one
    selected protein chain per input and does not align target groups.
    """
    contracts = {
        "analysis.interface_clustering": {
            "selection_schema": "studio.interface_clustering_selection.v1",
            "collection_schema": "studio.interface_clustering_collection.v1",
            "item_schema": "studio.interface_clustering_collection_item.v1",
            "maximum": 2000,
            "label": "Interface Clustering",
        },
        "analysis.foldseek": {
            "selection_schema": "studio.foldseek_selection.v1",
            "collection_schema": "studio.foldseek_collection.v1",
            "item_schema": "studio.foldseek_collection_item.v1",
            "maximum": 100000,
            "label": "Foldseek",
        },
    }
    contract = contracts.get(workflow_id)
    if contract is None:
        raise ValueError("Saved structure collections are not enabled for this workflow.")
    label = contract["label"]
    if not isinstance(selection, dict) or selection.get("schema") != contract["selection_schema"]:
        raise ValueError(f"{label} requires a versioned saved-structure collection.")
    sources = selection.get("sources")
    if not isinstance(sources, list) or not 2 <= len(sources) <= contract["maximum"]:
        raise ValueError(f"{label} requires between 2 and {contract['maximum']:,} saved structures.")
    if not str(user.get("id") or ""):
        raise PermissionError("A signed-in user is required to analyze saved structures.")

    from mn_molecular_design_studio.complex_analysis_comparison import (
        compatibility_report,
        protein_chain_numbering,
        target_structure_fingerprint,
    )
    from mn_molecular_design_studio.evaluations import file_digest
    from mn_molecular_design_studio import complex_analysis_ui as workspace

    requested_project = str(project_id or "")
    jobs = workspace._accessible_result_jobs(user)
    jobs_by_identity = {
        (str(job.get("task_group") or ""), str(job.get("run_id") or "")): job
        for job in jobs
    }
    if database is None:
        from mn_molecular_design_studio.molecular_database import MolecularDatabase
        database = MolecularDatabase()

    resolved: list[dict[str, Any]] = []
    seen_identities = set()
    for raw_selector in sources:
        if not isinstance(raw_selector, dict):
            raise ValueError(f"Every {label} input must be a saved structure selector.")
        kind = str(raw_selector.get("kind") or "").casefold()
        target_chains = [str(value).strip() for value in raw_selector.get("target_chains", []) if str(value).strip()]
        binder_chains = [str(value).strip() for value in raw_selector.get("binder_chains", []) if str(value).strip()]
        target_required = workflow_id == "analysis.interface_clustering"
        if ((target_required and not target_chains) or len(target_chains) != len(set(target_chains))
                or len(binder_chains) != len(set(binder_chains)) or set(target_chains) & set(binder_chains)):
            raise ValueError(f"Each {label} structure needs unique protein chain selections with no overlap.")
        if workflow_id == "analysis.foldseek" and len(binder_chains) != 1:
            raise ValueError("Foldseek requires exactly one selected protein chain per structure.")
        if not binder_chains:
            raise ValueError(f"Each {label} structure needs at least one selected protein chain.")

        if kind == "results":
            task_group = str(raw_selector.get("task_group") or "")
            run_id = str(raw_selector.get("run_id") or "")
            job = jobs_by_identity.get((task_group, run_id))
            if job is None:
                raise PermissionError("A selected Results structure is unavailable or access was denied.")
            if str(job.get("project_id") or "") != requested_project:
                raise PermissionError("All selected Results structures must belong to the selected analysis project.")
            structure_ref = str(raw_selector.get("structure_ref") or "")
            if not structure_ref or Path(structure_ref).is_absolute() or ".." in Path(structure_ref).parts:
                raise ValueError("The selected Results structure reference is invalid.")
            selected = next((item for item in workspace._job_structure_options(user, job)
                             if str(item.get("relative_path") or "") == structure_ref), None)
            if selected is None:
                raise PermissionError("A selected Results structure is unavailable or access was denied.")
            candidate = selected.get("candidate") if isinstance(selected.get("candidate"), dict) else {}
            candidate_id = str(candidate.get("candidate_id") or "")
            requested_candidate_id = str(raw_selector.get("candidate_id") or "")
            if requested_candidate_id and candidate_id != requested_candidate_id:
                raise ValueError("The selected Results structure no longer matches its candidate identity.")
            raw_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
            model_id = str(candidate.get("model_index") or raw_metadata.get("model_index") or selected["path"].name)
            identity = f"{task_group}/{run_id}:{structure_ref}"
            source_path = Path(selected["path"]).resolve(strict=True)
            source_project_id = str(job.get("project_id") or "")
            source_record = {
                "kind": "results", "task_group": task_group, "run_id": run_id,
                "structure_ref": structure_ref, "candidate_id": candidate_id,
                "model_id": model_id, "source_tool": str(candidate.get("source_tool") or job.get("tool") or ""),
                "evaluation_ids": [candidate_id] if candidate.get("metrics") and candidate_id else [],
                "job_code": str(job.get("job_code") or ""),
            }
        elif kind == "database":
            from mn_molecular_design_studio.complex_analysis_sources import authorize_database_structure

            database_selector = {
                "molecule_id": raw_selector.get("molecule_id"),
                "structure_id": raw_selector.get("structure_id"),
                "group_id": raw_selector.get("group_id"),
                "evaluation_id": "",
                "model_id": raw_selector.get("model_id"),
                "binder_chains": binder_chains,
                "target_chains": target_chains,
            }
            authorized = authorize_database_structure(
                user, database_selector, project_id=project_id, database=database,
                workflow_id=workflow_id,
            )
            source_path = Path(authorized["structure_path"]).resolve(strict=True)
            if (str(raw_selector.get("structure_digest") or "")
                    and str(raw_selector.get("structure_digest")) != authorized["structure_digest"]):
                raise ValueError("The selected Database representation changed after comparison was opened.")
            if str(raw_selector.get("model_id") or authorized["model_id"]) != authorized["model_id"]:
                raise ValueError("The selected Database representation no longer matches its model identity.")
            evaluation_ids = _collection_evaluation_ids(user, raw_selector, database=database)
            model_id = authorized["model_id"]
            identity = f"database:{authorized['molecule_id']}:{authorized['structure_id']}"
            source_project_id = str(authorized["project_id"] or "")
            candidate_id = ""
            candidate = {}
            raw_metadata = {}
            source_record = {
                "kind": "database", "molecule_id": authorized["molecule_id"],
                "structure_id": authorized["structure_id"], "group_id": authorized["group_id"],
                "model_id": model_id, "record_name": authorized["record_name"],
                "structure_label": authorized["structure_label"],
                "evaluation_ids": evaluation_ids,
            }
        else:
            raise ValueError(f"{label} accepts saved Results or Database structures.")

        if source_project_id != requested_project:
            raise PermissionError("All selected structures must share the analysis project's access boundary.")
        if identity in seen_identities:
            raise ValueError("Each saved structure representation may appear only once in the collection.")
        seen_identities.add(identity)
        suffix = _collection_structure_suffix(source_path)
        if suffix not in {".pdb", ".cif", ".mmcif", ".pdb.gz", ".cif.gz", ".mmcif.gz"}:
            raise ValueError(f"{label} requires saved PDB or mmCIF structures.")
        source_digest = file_digest(source_path)
        requested_digest = str(raw_selector.get("structure_digest") or "")
        if requested_digest and requested_digest != source_digest:
            raise ValueError("A selected structure changed after comparison was opened.")
        if raw_selector.get("model_id") and str(raw_selector["model_id"]) != model_id:
            raise ValueError("A selected structure no longer matches its model identity.")

        try:
            target_numbering = protein_chain_numbering(source_path, target_chains) if target_chains else []
            binder_numbering = protein_chain_numbering(source_path, binder_chains)
        except ValueError as exc:
            raise ValueError(f"Could not verify the selected protein chains: {exc}") from exc
        if (target_required and not target_numbering) or not binder_numbering:
            raise ValueError(f"Selected protein chain groups must contain readable residues for {label}.")
        target_digest = target_structure_fingerprint(source_path, target_chains) if target_chains else ""
        resolved.append({
            "identity": identity, "label": str(raw_selector.get("label") or identity),
            "kind": kind, "path": source_path, "source_digest": source_digest,
            "target_digest": target_digest, "model_id": model_id,
            "target_chains": target_chains, "binder_chains": binder_chains,
            "target_numbering": target_numbering,
            "target_chains_sequence": [item["sequence"] for item in target_numbering],
            "binder_numbering": binder_numbering,
            "binder_chains_sequence": [item["sequence"] for item in binder_numbering],
            "complex_type": "protein-protein", "candidate_id": candidate_id,
            "candidate": candidate, "candidate_metadata": raw_metadata,
            "source_record": source_record, "source_project_id": source_project_id,
            "evaluation_ids": list(source_record.get("evaluation_ids") or []),
            "database_source": None,
        })
        if kind == "database":
            resolved[-1]["database_source"] = {
                "schema": "studio.complex_analysis_database_source.v1",
                "kind": "database", "molecule_id": authorized["molecule_id"],
                "structure_id": authorized["structure_id"], "group_id": authorized["group_id"],
                "evaluation_id": resolved[-1]["evaluation_ids"][0] if resolved[-1]["evaluation_ids"] else "",
                "evaluation_ids": resolved[-1]["evaluation_ids"], "model_id": model_id,
                "structure_digest": source_digest, "record_name": authorized["record_name"],
                "structure_label": authorized["structure_label"],
                "binder_chains": binder_chains, "target_chains": target_chains,
            }
            resolved[-1]["candidate_metadata"] = {"components": authorized.get("record_components") or []}

    reference = resolved[0]
    if workflow_id == "analysis.interface_clustering":
        for item in resolved[1:]:
            report = compatibility_report(reference, item)
            if not report["compatible"]:
                raise ValueError("Selected structures are not compatible for Interface Clustering: " + "; ".join(report["errors"]))
            if item["target_digest"] != reference["target_digest"]:
                raise ValueError("Interface Clustering requires the identical target structure, author chain IDs, and numbering across the collection.")
            if item["target_chains"] != reference["target_chains"]:
                raise ValueError("Interface Clustering requires the same target author chain IDs in every structure.")

    manifest_sources = [{
        "identity": item["identity"], "kind": item["kind"], "source": item["source_record"],
        "source_structure_digest": item["source_digest"], "target_structure_digest": item["target_digest"],
        "model_id": item["model_id"], "binder_chains": item["binder_chains"],
        "target_chains": item["target_chains"], "evaluation_ids": item["evaluation_ids"],
    } for item in resolved]
    collection_identity = {"schema": contract["collection_schema"],
                           "project_id": requested_project, "sources": manifest_sources}
    collection_id = hashlib.sha256(json.dumps(collection_identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    bundle_dir = user_input_root(str(user["id"])) / "complex_analysis_collections" / collection_id
    bundle_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        bundle_dir.chmod(0o700)
    except OSError:
        pass
    structures_dir = bundle_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    candidates = []
    collection_items = []
    for index, item in enumerate(resolved):
        candidate_id = "comparison-" + hashlib.sha256(item["identity"].encode()).hexdigest()[:24]
        staged_path = structures_dir / (candidate_id + _collection_structure_suffix(item["path"]))
        content = _collection_structure_bytes(item["path"])
        staged_digest = hashlib.sha256(content).hexdigest()
        if not staged_path.is_file() or file_digest(staged_path) != staged_digest:
            staged_path.write_bytes(content)
        selection_item = {
            "schema": contract["item_schema"],
            "collection_id": collection_id, "candidate_id": candidate_id,
            "source_identity": item["identity"], "source_kind": item["kind"],
            "source": item["source_record"], "source_structure_digest": item["source_digest"],
            "input_structure_digest": staged_digest, "target_structure_digest": item["target_digest"],
            "model_id": item["model_id"], "binder_chains": item["binder_chains"],
            "target_chains": item["target_chains"], "evaluation_ids": item["evaluation_ids"],
        }
        raw_metadata = {
            **item["candidate_metadata"],
            "analysis_subject": "protein", "model_index": item["model_id"],
            "target_structure_digest": item["target_digest"],
            "complex_analysis_selection": selection_item,
        }
        if item.get("database_source"):
            raw_metadata["database_source"] = item["database_source"]
        candidate = {
            "candidate_id": candidate_id, "source_tool": "complex_analysis_saved_structure",
            "tool": "complex_analysis_saved_structure", "stage": "analysis",
            "complex_pdb": str(staged_path.resolve()), "binder_chains": item["binder_chains"],
            "target_chains": item["target_chains"], "binder_sequence": ":".join(item["binder_chains_sequence"]),
            "metrics": dict(item["candidate"].get("metrics") or {}),
            "parents": [item["candidate_id"]] if item["candidate_id"] else [item["identity"]],
            "raw_metadata": raw_metadata,
        }
        if item.get("database_source"):
            candidate["source_molecule_id"] = item["database_source"]["molecule_id"]
        candidates.append(candidate)
        collection_items.append(selection_item)

    candidates_path = bundle_dir / "analysis_subjects.jsonl"
    candidates_path.write_text("".join(json.dumps(candidate, sort_keys=True) + "\n" for candidate in candidates), encoding="utf-8")
    collection_manifest = {
        **collection_identity, "collection_id": collection_id,
        "selection_key": saved_structure_collection_selection_key(selection, project_id),
        "owner_id": str(user["id"]),
        "items": collection_items,
    }
    (bundle_dir / "complex_analysis_collection.json").write_text(
        json.dumps(collection_manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8",
    )
    (bundle_dir / "metadata.json").write_text(json.dumps({
        "schema": "studio.complex_analysis_collection_source.v1",
        "source_kind": "complex_analysis_collection", "source_collection_id": collection_id,
        "owner_id": str(user["id"]), "studio_owner_id": str(user["id"]),
        "project_id": requested_project,
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    (bundle_dir / "molecular_design_studio.json").write_text(json.dumps({
        "application": "mn-molecular-design-studio", "studio_owner_id": str(user["id"]),
        "project_id": requested_project, "source_kind": "complex_analysis_collection",
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return {
        "collection_id": collection_id, "source_run_dir": str(bundle_dir.resolve()),
        "candidates_jsonl": str(candidates_path.resolve()),
        "collection_manifest": collection_manifest,
    }


def build_interface_clustering_collection(user: dict, selection: dict, *, project_id: str | None,
                                          database=None) -> dict[str, Any]:
    """Backward-compatible Interface Clustering collection adapter."""
    return build_saved_structure_collection(
        user, selection, workflow_id="analysis.interface_clustering", project_id=project_id,
        database=database,
    )


def build_foldseek_collection(user: dict, selection: dict, *, project_id: str | None,
                              database=None) -> dict[str, Any]:
    """Authorize and stage a Foldseek cohort without imposing target compatibility."""
    return build_saved_structure_collection(
        user, selection, workflow_id="analysis.foldseek", project_id=project_id,
        database=database,
    )


def apply_database_analysis_source(payload: dict, user: dict, project_id: str | None) -> dict:
    """Turn a stable source selector into parameters accepted by the workflow registry."""
    selection = payload.get("database_source")
    if selection is None:
        return payload
    workflow_id = str(payload.get("workflow") or "")
    if workflow_id not in _WORKFLOWS:
        raise ValueError("Database structure submission is supported only for registered contextual analysis workflows.")
    if payload.get("source_job"):
        raise ValueError("Use either a Results source job or a Database structure source, not both.")
    raw_parameters = payload.get("parameters")
    if not isinstance(raw_parameters, dict):
        raise ValueError("A Database structure analysis requires one parameter object.")
    adapter = build_database_analysis_source(
        user, selection, workflow_id=workflow_id, project_id=project_id,
    )
    parameters = dict(raw_parameters)
    parameters.update({
        "source_run_dir": adapter["source_run_dir"],
        "candidates_jsonl": adapter["candidates_jsonl"],
        "candidate_ids": [adapter["candidate_id"]],
        "selected_structure_path": adapter["selected_structure_path"],
        "selected_model_id": adapter["selected_model_id"],
    })
    if workflow_id in _LIGAND_WORKFLOWS:
        parameters.pop("binder_chains", None)
        parameters.pop("target_chains", None)
        parameters["ligand_selection"] = dict(adapter["source"]["ligand_selection"])
    else:
        parameters["binder_chains"] = list(adapter["source"]["binder_chains"])
        if workflow_id != "analysis.pep_patch":
            parameters["target_chains"] = list(adapter["source"]["target_chains"])
    result = {key: value for key, value in payload.items() if key != "database_source"}
    result["parameters"] = parameters
    return result
