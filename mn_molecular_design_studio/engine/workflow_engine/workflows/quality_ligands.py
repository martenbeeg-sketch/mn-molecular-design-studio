"""Prepare an explicit ligand pose and receptor for the existing analysis images."""
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from pathlib import PurePosixPath


def docking_pose_options(source_run_dir=None, *, source_job_reference="", source_node_id="",
                         scores_text=None, compound_set_payload=None, artifact_paths=None):
    """Describe individually selectable poses from a completed docking result.

    The result's score table and consumed compound set are the only sources of
    pose/chemistry identity.  Remote callers may supply those two small saved
    artifacts and the source job's artifact manifest; coordinate files are
    still resolved and hashed by the owning worker at execution time.
    """
    from mn_molecular_design_studio.ligand_compounds import (
        compound_content_sha256, load_compound_set, normalize_compounds,
    )

    source = Path(source_run_dir).expanduser().resolve(strict=True) if source_run_dir else None
    if source and (not source_job_reference or not source_node_id):
        metadata_path = source / "metadata.json"
        try:
            source_metadata = json.loads(metadata_path.read_text())
        except (OSError, ValueError):
            source_metadata = {}
        source_job_reference = source_job_reference or str(source_metadata.get("source_job_reference") or
            f"{source_metadata.get('task_group') or source.parent.name}:{source_metadata.get('run_id') or source.name}")
        if not source_node_id:
            from mn_molecular_design_studio.cluster import config
            source_node_id = config.identity()
    table = source / "artifacts/docking_scores.csv" if source else None
    compound_path = source / "artifacts/input/compound_set.json" if source else None
    if scores_text is None and table and table.is_file():
        scores_text = table.read_text(encoding="utf-8-sig")
    if compound_set_payload is None and compound_path and compound_path.is_file():
        compound_set_payload = load_compound_set(compound_path)
    if not isinstance(scores_text, str) or not scores_text.strip() or compound_set_payload is None:
        raise ValueError("Docking pose selection requires the saved docking score table and consumed compound set.")
    if isinstance(compound_set_payload, (str, bytes)):
        try:
            compound_set_payload = json.loads(compound_set_payload)
        except (ValueError, UnicodeDecodeError) as exc:
            raise ValueError("The saved compound set is invalid JSON.") from exc
    if not isinstance(compound_set_payload, dict) or compound_set_payload.get("record_type") != "compound_set":
        raise ValueError("The selected docking result does not contain a typed compound set.")
    compounds = normalize_compounds(compound_set_payload.get("compounds") or [])
    recorded_digest = str(compound_set_payload.get("content_sha256") or "")
    if recorded_digest and recorded_digest != compound_content_sha256(compounds):
        if compound_set_payload.get("schema_version") == 1:
            legacy = [{key: row[key] for key in ("compound_id", "smiles", "properties") if key in row}
                      for row in compounds]
            legacy_digest = hashlib.sha256(json.dumps(legacy, ensure_ascii=False, sort_keys=True,
                separators=(",", ":")).encode("utf-8")).hexdigest()
            if recorded_digest != legacy_digest:
                raise ValueError("The consumed compound set failed its content-integrity check.")
        else:
            raise ValueError("The consumed compound set failed its content-integrity check.")
    compound_by_id = {str(row["compound_id"]): row for row in compounds}

    header = scores_text.splitlines()[0]
    rows = list(csv.DictReader(io.StringIO(scores_text), delimiter="\t" if "\t" in header else ","))
    declared = {str(value).replace("\\", "/").lstrip("/") for value in (artifact_paths or []) if value}
    def is_declared(relative):
        return not declared or relative in declared or any(value.endswith("/" + relative) for value in declared)
    receptor_relative = next((name for name in ("artifacts/input/receptor.pdb", "artifacts/input/target.pdb")
                              if (source and (source / name).is_file()) or is_declared(name)), "")
    if not receptor_relative:
        raise ValueError("This docking result has no authorized saved receptor structure.")
    receptor = _resolve_source_artifact(source, receptor_relative) if source else None
    if source and (not receptor or not receptor.is_file()):
        raise ValueError("This docking result's receptor artifact is missing or outside its source job.")

    options = []
    for index, row in enumerate(rows):
        compound_id = str(row.get("compound_id") or row.get("compound") or "").strip()
        pose_text = str(row.get("pose_path") or "").strip()
        status = str(row.get("status") or "completed").strip().casefold()
        if not compound_id or not pose_text or status not in {"completed", "succeeded", "success", ""}:
            continue
        relative = _relative_source_artifact(source, pose_text)
        if not relative:
            raise ValueError("A docking score row refers to a pose outside its source job.")
        if not is_declared(relative):
            raise ValueError("A selected pose is not present in the authorized source artifact manifest.")
        pose_path = _resolve_source_artifact(source, relative) if source else None
        if source and (pose_path is None or not pose_path.is_file()):
            continue
        compound = compound_by_id.get(compound_id)
        native_pose_id = str(row.get("pose_id") or "").strip()
        identity_parts = [source_job_reference, source_node_id, compound_id, native_pose_id,
                          row.get("pose_rank") or row.get("rank") or "", row.get("replicate") or "",
                          row.get("attempt_index") or "", relative, index]
        pose_id = "dockpose-" + hashlib.sha256(json.dumps(identity_parts, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()[:32]
        label_parts = [compound.get("compound_name", compound_id) if compound else compound_id]
        rank = str(row.get("pose_rank") or row.get("rank") or "").strip()
        if rank:
            label_parts.append("pose " + rank)
        score = next((f"{key} {float(row[key]):g}" for key in
                      ("estimated_dg_reu", "best_score_kcal_mol", "score", "total_score_reu", "cnn_score", "cnn_affinity")
                      if row.get(key) not in (None, "") and _finite_float(row[key]) is not None), "")
        if score:
            label_parts.append(score)
        options.append({"pose_id": pose_id, "candidate_id": pose_id, "compound_id": compound_id,
                        "compound_name": compound.get("compound_name", compound_id) if compound else compound_id,
                        "smiles": str(compound.get("smiles") or "") if compound else "",
                        "pose_path": str(pose_path) if pose_path else relative,
                        "pose_relative_path": relative, "receptor_path": str(receptor) if receptor else receptor_relative,
                        "source_job_reference": source_job_reference, "source_node_id": source_node_id,
                        "score_row": dict(row), "compound_record": dict(compound) if compound else None,
                        "available": bool(compound and compound.get("smiles")),
                        "row_index": index,
                        "unavailable_reason": "The consumed compound set has no chemistry for this compound." if not compound else
                                              "The compound has no recorded SMILES." if not compound.get("smiles") else "",
                        "label": " · ".join(label_parts)})
    return options


def docking_candidates(source_run_dir, candidate_ids=None, *, source_job_reference="", source_node_id="",
                       source_workstation=""):
    """Adapt selected raw docking poses to the shared normalized-candidate row.

    No protein sequence is invented: ligand chemistry and pose provenance live
    in typed component metadata, while receptor and pose remain structure paths.
    The no-selection legacy path retains its historical best-pose-per-compound
    behavior for standalone PLIP.
    """
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import normalize_candidate
    from mn_molecular_design_studio.evaluations import file_digest
    from mn_molecular_design_studio.molecular_database_ingest import _best_docking_pose_indices

    source = Path(source_run_dir).expanduser().resolve(strict=True)
    options = docking_pose_options(source, source_job_reference=source_job_reference, source_node_id=source_node_id)
    if candidate_ids is not None:
        selected_ids = [str(value).strip() for value in candidate_ids if str(value).strip()]
        if not selected_ids or len(selected_ids) != len(set(selected_ids)):
            raise ValueError("Select one or more unique docking pose identities.")
        by_id = {row["pose_id"]: row for row in options if row["available"]}
        missing = [value for value in selected_ids if value not in by_id]
        if missing:
            unavailable = {row["pose_id"]: row.get("unavailable_reason") for row in options}
            detail = next((unavailable[value] for value in missing if value in unavailable), None)
            raise ValueError(detail or "One or more selected docking poses are unavailable in the authorized source job.")
        chosen = [by_id[value] for value in selected_ids]
    else:
        from mn_molecular_design_studio.ligand_compounds import load_compound_set
        compound_path = source / "artifacts/input/compound_set.json"
        compound_set = load_compound_set(compound_path)
        with (source / "artifacts/docking_scores.csv").open(encoding="utf-8-sig") as handle:
            header = handle.readline()
            handle.seek(0)
            score_rows = list(csv.DictReader(handle, delimiter="\t" if "\t" in header else ","))
        pose_paths = [_resolve_source_artifact(source, _relative_source_artifact(source, str(row.get("pose_path") or "")) or "")
                      for row in score_rows]
        best = set(_best_docking_pose_indices(score_rows, pose_paths))
        chosen = [row for row in options if row["available"] and row["row_index"] in best]

    result_path = source / "result.json"
    if result_path.is_file() and json.loads(result_path.read_text()).get("success") is not True:
        raise ValueError("Only a successfully completed docking result can be analyzed.")
    metadata = json.loads((source / "metadata.json").read_text()) if (source / "metadata.json").is_file() else {}
    input_data = json.loads((source / "input.json").read_text()) if (source / "input.json").is_file() else {}
    parameters = input_data.get("params") if isinstance(input_data.get("params"), dict) else {}
    artifact_contract = {}
    contract_path = source / "artifact_contract.json"
    if contract_path.is_file():
        try:
            artifact_contract = json.loads(contract_path.read_text())
        except (OSError, ValueError):
            artifact_contract = {}
    source_node_id = source_node_id or str(metadata.get("source_node_id") or "")
    if not source_node_id:
        from mn_molecular_design_studio.cluster import config
        source_node_id = config.identity()
    source_job_reference = source_job_reference or str(metadata.get("source_job_reference") or
        f"{metadata.get('task_group') or source.parent.name}:{source.name}")
    receptor = Path(chosen[0]["receptor_path"]).resolve(strict=True)
    project_id = str(metadata.get("studio_project_id") or metadata.get("project_id") or "")
    target_chains = parameters.get("receptor_chains") or parameters.get("target_chains") or parameters.get("chain_ids")
    if isinstance(target_chains, str):
        target_chains = [value for value in target_chains.replace(",", " ").split() if value]
    target_chains = [str(value) for value in target_chains] if isinstance(target_chains, list) else []
    prepared_records = _load_preparation_records(source)
    candidates = []
    for option in chosen:
        pose_path = Path(option["pose_path"]).resolve(strict=True)
        score_row = option["score_row"]
        compound_id = option["compound_id"]
        compound_path = source / "artifacts/input/compound_set.json"
        compound = option.get("compound_record")
        if compound is None or not compound.get("smiles"):
            raise ValueError(f"Docking pose {option['pose_id']} has no recorded ligand chemistry.")
        score_metadata = {key: value for key, value in score_row.items() if value not in (None, "")}
        metrics = {key: number for key, value in score_metadata.items() if (number := _finite_float(value)) is not None}
        pose_digest, receptor_digest = file_digest(pose_path), file_digest(receptor)
        from mn_molecular_design_studio.cluster import config
        task_group = str(metadata.get("task_group") or source.parent.name)
        source_run_id = str(metadata.get("run_id") or source.name)
        source_prefix = config.artifact_uri(source_node_id, task_group, source_run_id)
        artifact_references = {"docking_scores": source_prefix + "/artifacts/docking_scores.csv",
            "compound_set": source_prefix + "/artifacts/input/compound_set.json",
            "receptor": source_prefix + "/artifacts/input/" + receptor.name,
            "pose": source_prefix + "/" + option["pose_relative_path"]}
        preparation = prepared_records.get(compound_id)
        ligand_properties = dict(compound.get("properties") or {})
        if preparation:
            ligand_properties["docking_preparation"] = preparation
        candidate = normalize_candidate({
            "candidate_id": option["pose_id"], "source_tool": str(score_row.get("engine") or metadata.get("tool") or "ligand_docking"),
            "stage": "analysis", "target_pdb": str(receptor), "complex_pdb": str(pose_path),
            "target_chains": target_chains, "binder_chains": [], "parents": [compound_id],
            "metrics": metrics,
            "raw_metadata": {
                "artifact_kind": "protein_ligand_pose", "artifact_schema": "ligand_docking_result/v1",
                "source_job_reference": source_job_reference, "source_job_id": str(metadata.get("run_id") or source.name),
                "source_workstation_id": source_node_id, "source_workstation": source_workstation or str(metadata.get("source_workstation") or ""),
                "source_project_id": project_id, "source_candidate_id": compound_id,
                "source_artifact_references": artifact_references,
                "source_artifact_sha256": {"docking_scores": file_digest(source / "artifacts/docking_scores.csv"),
                    "compound_set": file_digest(compound_path), "receptor": receptor_digest, "pose": pose_digest},
                "source_structure_digest": pose_digest, "pose_id": option["pose_id"],
                "native_pose_id": score_row.get("pose_id") or None, "pose_rank": score_row.get("pose_rank") or score_row.get("rank") or None,
                "replicate": score_row.get("replicate") or None, "attempt_index": score_row.get("attempt_index") or None,
                "model_id": score_row.get("receptor_model_id") or score_row.get("model_id") or parameters.get("receptor_model_id") or None,
                "ligand_chain": score_row.get("ligand_chain") or (compound.get("properties") or {}).get("ligand_chain") or None,
                "ligand_residue_key": score_row.get("ligand_residue_key") or (compound.get("properties") or {}).get("ligand_residue_key") or None,
                "ligand_ccd": score_row.get("ligand_ccd") or (compound.get("properties") or {}).get("ligand_ccd") or None,
                "ligand_atom_mapping": score_row.get("atom_mapping") or (compound.get("properties") or {}).get("atom_mapping") or None,
                "receptor_identity": {"role": "receptor", "structure_sha256": receptor_digest,
                    "selected_chains": target_chains or None, "chain_selection_status": "recorded" if target_chains else "not_recorded",
                    "model_id": score_row.get("receptor_model_id") or score_row.get("model_id") or parameters.get("receptor_model_id") or None},
                "ligand_identity": {"role": "small_molecule", "compound_id": compound_id,
                    "name": compound.get("compound_name") or compound_id, "smiles": compound["smiles"],
                    "properties": ligand_properties, "atom_mapping": score_row.get("atom_mapping") or None,
                    "stereochemistry": "encoded_in_recorded_smiles" if any(mark in compound["smiles"] for mark in ("@", "/", "\\")) else "not_recorded"},
                "docking_engine": score_row.get("engine") or metadata.get("engine") or metadata.get("tool") or "not_recorded",
                "docking_protocol": score_row.get("protocol") or parameters.get("protocol") or metadata.get("workflow") or "not_recorded",
                "docking_protocol_parameters": parameters,
                "docking_scores": score_metadata,
                "docking_score_units": next((entry.get("score_units") for entry in artifact_contract.get("outputs", [])
                    if isinstance(entry, dict) and entry.get("name") in {"docking_scores", "scores"}), None),
                "source_artifact_contract": artifact_contract,
                "source_compound_id": compound_id, "compound_id": compound_id,
                "compound_name": compound.get("compound_name") or compound_id, "name": compound.get("compound_name") or compound_id,
                "ligand_smiles": compound["smiles"], "ligand_preparation": preparation,
                "components": [{"type": "protein", "role": "receptor", "structure_sha256": receptor_digest,
                                "selected_chains": target_chains or None},
                               {"type": "small_molecule", "role": "ligand", "compound_id": compound_id,
                                "smiles": compound["smiles"], "properties": ligand_properties,
                                "pose_sha256": pose_digest}],
            },
        })
        candidates.append(candidate)
    return candidates


def _finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _relative_source_artifact(source, value):
    text = str(value or "").strip().replace("\\", "/")
    if not text or text.startswith(("cluster://", "runs://", "http://", "https://")):
        return ""
    path = Path(text).expanduser()
    if path.is_absolute():
        if source is None:
            return ""
        try:
            return path.resolve(strict=False).relative_to(source.resolve()).as_posix()
        except ValueError:
            return ""
    relative = PurePosixPath(text.lstrip("/"))
    if not relative.parts or ".." in relative.parts or "." in relative.parts:
        return ""
    if source is not None:
        try:
            resolved = (source / Path(*relative.parts)).resolve(strict=False)
            if not resolved.is_relative_to(source.resolve()):
                return ""
        except OSError:
            return ""
    return relative.as_posix()


def _resolve_source_artifact(source, relative):
    if source is None:
        return None
    safe = _relative_source_artifact(source, relative)
    if not safe:
        return None
    path = (source / Path(*PurePosixPath(safe).parts)).resolve(strict=False)
    return path if path.is_relative_to(source.resolve()) else None


def _load_preparation_records(source):
    path = source / "artifacts/preparation_summary.json"
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    if isinstance(payload, dict):
        records = payload.get("compounds") or payload.get("results") or []
        if isinstance(records, list):
            return {str(row.get("compound_id")): row for row in records if isinstance(row, dict) and row.get("compound_id")}
    return {}


def prepare(item, candidate, run_dir, *, ligand_selection=None):
    """Select exactly one ligand, then match its graph to the declared SMILES.

    Ambiguous complexes and failed topology assignment are explicit failures;
    atom order is never assumed to imply matching chemical connectivity.
    """
    import gemmi
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from mn_molecular_design_studio.docking_results import pdbqt_to_pdb
    raw = candidate.get("raw_metadata") or {}
    components = raw.get("components") or []
    ligands = [c for c in components if c.get("type") == "small_molecule"]
    saved_selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
    ligand_selection = ligand_selection or saved_selection.get("ligand_selection") or {}
    if not isinstance(ligand_selection, dict):
        raise ValueError("The explicit ligand selection is invalid.")
    selected_component = None
    if ligand_selection:
        selected_component = next((component for component in ligands
            if (not ligand_selection.get("residue_key") or str(component.get("residue_key") or "") == str(ligand_selection["residue_key"]))
            and (not ligand_selection.get("ccd") or str(component.get("ccd") or "").upper() == str(ligand_selection["ccd"]).upper())
            and (not ligand_selection.get("chain") or str(component.get("residue_key") or "").split(":", 1)[0] == str(ligand_selection["chain"])
                 or str(component.get("chain") or "") == str(ligand_selection["chain"]))), None)
        if ligands and selected_component is None and not ligand_selection.get("smiles"):
            raise ValueError("The selected ligand identity is not present in the source candidate metadata.")
    ligand_smiles = ligand_selection.get("smiles") or (selected_component or {}).get("smiles")
    smiles = ligand_smiles or raw.get("ligand_smiles") or raw.get("smiles") or (ligands[0].get("smiles") if len(ligands) == 1 else "")
    if isinstance(smiles, list):
        smiles = smiles[0] if len(smiles) == 1 else ""
    if not isinstance(smiles, str) or not smiles:
        raise ValueError("Declare one ligand SMILES in the source prediction; multiple ligands require an explicit ligand selection.")
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError("The source ligand SMILES is invalid.")
    directory = Path(run_dir) / "artifacts/input" / item["id"]
    directory.mkdir(parents=True, exist_ok=True)
    source = Path(candidate["complex_pdb"])
    source_digest = item.get("source_structure_digest") or _digest(source)
    receptor = Path(candidate["target_pdb"]) if candidate.get("target_pdb") else None
    requested_chain = ligand_selection.get("chain") or raw.get("ligand_chain") or (
        (selected_component.get("residue_key") or "").split(":", 1)[0] if selected_component else
        ((ligands[0].get("residue_key") or "").split(":", 1)[0] if len(ligands) == 1 else ""))
    requested_residue = ligand_selection.get("residue_key") or (selected_component or {}).get("residue_key") or ""
    requested_ccd = ligand_selection.get("ccd") or (selected_component or {}).get("ccd") or ""
    selected_model_id = str(item.get("model_id") or saved_selection.get("model_id") or "").strip()
    if source.suffix.lower() in {".sdf", ".mol"}:
        poses = [m for m in Chem.SDMolSupplier(str(source), removeHs=False) if m is not None]
        if len(poses) != 1:
            raise ValueError("Select one coordinate-bearing ligand pose per analysis input.")
        molecule = poses[0]
        if molecule is None:
            raise ValueError("The docked SDF pose could not be read.")
        protein_text = receptor.read_text() if receptor else ""
    else:
        if source.suffix.lower() == ".pdbqt":
            text = pdbqt_to_pdb(source.read_text())
            protein_text = receptor.read_text() if receptor else ""
        elif source.suffix.lower() in {".cif", ".mmcif", ".pdb", ".ent"}:
            text = _selected_model_pdb(source, selected_model_id)
            protein_text = ""
        else:
            text = source.read_text()
            protein_text = ""
        protein_lines = [line for line in text.splitlines() if line.startswith("ATOM  ") and gemmi.find_tabulated_residue(line[17:20].strip()).is_amino_acid()]
        if protein_lines:
            protein_text = "\n".join(protein_lines) + "\nEND\n"
        elif not protein_text and receptor:
            protein_text = receptor.read_text()
        ligand_rows = [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM")) and not gemmi.find_tabulated_residue(line[17:20].strip()).is_amino_acid() and line[17:20].strip() not in {"HOH", "WAT", "DOD"}]
        residues = {(line[21:22], line[22:27], line[17:20]) for line in ligand_rows}
        if requested_chain:
            residues = {r for r in residues if r[0].strip() == str(requested_chain)}
        if requested_residue:
            residues = {r for r in residues if f"{r[0].strip()}:{r[1].strip()}" == str(requested_residue)}
        if requested_ccd and len(str(requested_ccd).strip()) <= 3:
            residues = {r for r in residues if r[2].strip().upper() == str(requested_ccd).strip().upper()}
        # PDB's residue-name field is limited to three characters. Gemmi
        # truncates longer mmCIF component IDs when exporting this selected
        # model (for example, LIG_L becomes LIG). In that case, the source's
        # validated chain/residue selector identifies the pose; the subsequent
        # graph-to-template match validates the exact chemical identity.
        # Ignore isolated ions, but never choose the largest of several organic ligands.
        residues = {r for r in residues if sum((line[21:22], line[22:27], line[17:20]) == r for line in ligand_rows) > 1}
        if len(residues) != 1:
            raise ValueError("This complex has multiple or no unambiguous organic ligand residues. Select a single ligand in the source input.")
        selected = next(iter(residues))
        block = "\n".join(line for line in ligand_rows if (line[21:22], line[22:27], line[17:20]) == selected) + "\nEND\n"
        molecule = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False)
    if molecule is None or not protein_text:
        raise ValueError("A coordinate-bearing ligand and receptor are both required.")
    molecule = AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(template), Chem.RemoveHs(molecule, sanitize=False))
    if molecule.GetNumConformers() != 1:
        raise ValueError("Ligand coordinates are missing.")
    ligand_file = directory / "ligand.sdf"
    with Chem.SDWriter(str(ligand_file)) as writer:
        writer.write(molecule)
    receptor_file = directory / "receptor.pdb"
    receptor_file.write_text(protein_text)
    from mn_molecular_design_studio.evaluations import file_digest
    ligand_pose_digest = file_digest(ligand_file)
    receptor_digest = file_digest(receptor_file)
    canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
    ligand_provenance = {
        "identity": {"ccd": str(requested_ccd or (selected_component or {}).get("ccd") or ""),
                     "smiles": canonical_smiles, "residue_key": str(requested_residue or (selected_component or {}).get("residue_key") or ""),
                     "chain": str(requested_chain or "")},
        "source_structure_digest": source_digest,
        "pose_coordinates_sha256": ligand_pose_digest,
        "receptor_coordinates_sha256": receptor_digest,
        "model_id": str(item.get("model_id") or saved_selection.get("model_id") or ""),
    }
    item.update(source_structure_digest=source_digest, ligand_pose_digest=ligand_pose_digest,
                receptor_digest=receptor_digest, ligand_identity=ligand_provenance["identity"],
                ligand_provenance=ligand_provenance)
    item.update(smiles=smiles, ligand_file="/work/" + str(ligand_file.relative_to(run_dir)), receptor_file="/work/" + str(receptor_file.relative_to(run_dir)))
    protein_lines = [line for line in protein_text.splitlines() if line.startswith("ATOM  ")]
    chains = {line[21:22] for line in protein_lines}
    ligand_chain = next(c for c in "LZBCDEFGHIJKLMNOPQRSTUVWXYAabcdefghijklmnopqrstuvwxyz" if c not in chains)
    ligand_lines = []
    for serial, line in enumerate((line for line in Chem.MolToPDBBlock(molecule).splitlines() if line.startswith(("ATOM  ", "HETATM"))), len(protein_lines) + 1):
        ligand_lines.append(f"HETATM{serial:5d}" + line[11:17] + "LIG" + line[20:21] + ligand_chain + "   1" + line[26:])
    complex_file = directory / "analyzed_complex.pdb"
    complex_file.write_text("\n".join([*protein_lines, "TER", *ligand_lines, "END", ""]))
    candidate.update(complex_pdb=str(complex_file), target_pdb=str(receptor_file), binder_chains=[ligand_chain], target_chains=sorted(chains))
    candidate["raw_metadata"] = {**raw, "analysis_subject": "ligand", "ligand_smiles": smiles,
                                 "ligand_provenance": ligand_provenance,
                                 "source_compound_id": raw.get("source_compound_id") or item["candidate_id"],
                                 "compound_id": raw.get("source_compound_id") or raw.get("compound_id") or item["candidate_id"],
                                 "compound_name": raw.get("compound_name") or raw.get("name") or item["candidate_id"]}
    item["structure"] = "/work/" + str(complex_file.relative_to(run_dir))
    return ligand_file, receptor_file


def _digest(path):
    from mn_molecular_design_studio.evaluations import file_digest
    return file_digest(path)


def _selected_model_pdb(source, selected_model_id):
    """Export exactly one source model; never concatenate alternate models."""
    import gemmi
    structure = gemmi.read_structure(str(source))
    if not len(structure):
        raise ValueError("The selected structure contains no coordinate model.")
    if len(structure) > 1:
        models = list(structure)
        selected = next((model for model in models if selected_model_id and str(model.name) == selected_model_id), None)
        if selected is None and str(selected_model_id).isdigit():
            index = int(selected_model_id) - 1
            selected = models[index] if 0 <= index < len(models) else None
        if selected is None:
            raise ValueError("The saved structure contains multiple models and the selected model identity cannot be resolved unambiguously.")
        one_model = gemmi.Structure()
        one_model.cell = structure.cell
        one_model.spacegroup_hm = structure.spacegroup_hm
        one_model.add_model(selected.clone())
        structure = one_model
    return structure.make_pdb_string()
