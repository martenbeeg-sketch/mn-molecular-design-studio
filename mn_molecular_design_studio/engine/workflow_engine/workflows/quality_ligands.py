"""Prepare an explicit ligand pose and receptor for the existing analysis images."""
import csv
import json
from pathlib import Path


def docking_candidates(source_run_dir):
    from mn_molecular_design_studio.ligand_compounds import load_compound_set
    from mn_molecular_design_studio.molecular_database_ingest import _best_docking_pose_indices, _path
    source = Path(source_run_dir)
    table = source / "artifacts/docking_scores.csv"
    compound_path = source / "artifacts/input/compound_set.json"
    if not table.is_file() or not compound_path.is_file():
        raise ValueError("Select docking results with saved poses and the consumed compound set, or normalized cofolding candidates.")
    compounds = {c["compound_id"]: c for c in load_compound_set(compound_path)["compounds"]}
    with table.open() as handle:
        header = handle.readline()
        handle.seek(0)
        rows = list(csv.DictReader(handle, delimiter="\t" if "\t" in header else ","))
    poses = [_path(row.get("pose_path"), source) for row in rows]
    receptor = next((source / name for name in ("artifacts/input/receptor.pdb", "artifacts/input/target.pdb") if (source / name).is_file()), None)
    if receptor is None:
        raise ValueError("This docking result has no saved receptor.")
    candidates = []
    for i in sorted(_best_docking_pose_indices(rows, poses)):
        row, compound = rows[i], compounds.get(rows[i]["compound_id"])
        if not compound or str(row.get("status", "completed")).lower() not in {"completed", "succeeded", ""}:
            continue
        candidates.append({"candidate_id": compound["compound_id"], "source_tool": "docking_selection", "stage": "analysis",
            "complex_pdb": str(poses[i]), "target_pdb": str(receptor), "binder_chains": [], "target_chains": [],
            "raw_metadata": {"name": compound.get("compound_name") or compound["compound_id"], "compound_name": compound.get("compound_name") or compound["compound_id"],
                             "source_compound_id": compound["compound_id"], "ligand_smiles": compound["smiles"]}, "metrics": {}})
    if not candidates:
        raise ValueError("No successful docking poses are available for analysis.")
    return candidates


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
