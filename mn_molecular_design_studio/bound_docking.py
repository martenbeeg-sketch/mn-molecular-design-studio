"""Reusable docking poses with a shared receptor and exact-coordinate MD inputs."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def register_receptor(database, user, path, *, project_id="", name="Docking receptor"):
    """Store one portable receptor snapshot, shared by all poses in its frame."""
    from mn_molecular_design_studio import molecule_imports
    from mn_molecular_design_studio.molecular_database import _id
    from mn_molecular_design_studio.molecular_database_ingest import _find_original_group

    path = Path(path)
    content = path.read_bytes()
    record = molecule_imports.parse_structure(content, path.name, include_ligands=False)[0]
    record.update(name=name, files=[])
    molecule, group, matched = _find_original_group(database, user["id"], project_id, record)
    if not matched:
        group = database.create_group(user, "Docking targets", project_id=project_id, import_method="tool_output",
            group_id=_id("docking-targets-", f"{user['id']}:{project_id}"),
            origin={"type": "user", "collection_kind": "curated", "purpose": "Reusable docking receptors"})
        with database.connection() as connection:
            molecule = database._put(connection, group, record, preserve_existing=True)
    digest = hashlib.sha256(content).hexdigest()
    stored = database._store_file(path.name, content)
    structure = database.attach_structure(group, molecule,
        {"path": str(stored), "label": "Docking receptor · pose coordinate frame", "dimensions": 3,
         "portable": True, "portable_sha256": digest, "database_structure_policy": "docking_receptor"})
    return {"target_molecule_id": molecule, "target_structure_id": structure, "target_path": str(stored),
            "target_structure_digest": digest, "target_variant_sha256": digest}


def screening_receptor(run_dir, target):
    """Keep the actual docked receptor, including any rigid axis alignment."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.runtime import runs_root

    run_dir = Path(run_dir)
    destination = run_dir / "artifacts/active_learning/receptor.pdb"
    if destination.is_file():
        return destination
    history = read_json(run_dir / "artifacts/active_learning/round_history.json")
    for row in history.get("rounds", []):
        for worker in row.get("workers", []):
            if worker.get("host"):
                continue
            for reference in worker.get("run_dirs") or [worker.get("run_dir", "")]:
                folder = runs_root() / "ligand-docking" / Path(reference).name
                snapshot = folder / "artifacts/input/receptor.pdb"
                if snapshot.is_file():
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    destination.write_bytes(snapshot.read_bytes())
                    return destination
    source = target.get("target_pdb")
    if not source or not Path(source).is_file():
        return None
    text = Path(source).read_text()
    if target.get("align_longest_axis_to_x"):
        from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import (
            ligand_longest_axis_transform, transform_pdb_data)
        text = transform_pdb_data(text, ligand_longest_axis_transform(Path(target["alignment_ligand_sdf"])))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text)
    return destination


def receptor_file(database, user, structure):
    context = structure.get("context") or {}
    if context.get("target_structure_id"):
        # Resolve through the Database ACL and replica mapping before host paths.
        return database.structure_file(user, context["target_structure_id"])
    if context.get("target_path") and Path(context["target_path"]).is_file():
        return Path(context["target_path"])
    if context.get("target_molecule_id"):
        rows = database.list_molecules(user, molecule_ids=[context["target_molecule_id"]], limit=1)["items"]
        for target in rows:
            for item in target["structures"]:
                if item["format"] in {"pdb", "ent", "cif", "mmcif"}:
                    return database.structure_file(user, item["id"])
    return None


def is_bound_pose(structure):
    context = structure.get("context") or {}
    return bool(context.get("target_structure_id") or context.get("target_path") or context.get("target_molecule_id")) and (
        context.get("database_structure_policy") in {"best_docking_pose", "best_prediction"})


def _canonical_isomeric_smiles(smiles):
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(str(smiles or ""))
    if molecule is None:
        return ""
    molecule = Chem.RemoveHs(molecule)
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(molecule)
    return Chem.MolToSmiles(molecule, isomericSmiles=True, canonical=True)


def pose_molecule(path, context, smiles, *, return_provenance=False):
    """Reconstruct chemistry from the preparation map, preserving docked atoms."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    path = Path(path)
    text = path.read_text()
    source_identity = _canonical_isomeric_smiles(smiles)
    if not source_identity:
        raise ValueError("The saved Database record lacks a valid isomeric SMILES identity for this docking pose.")
    provenance = {"atom_mapping": [], "pose_isomeric_smiles": "", "source_isomeric_smiles": str(smiles)}
    if path.suffix.lower() == ".pdbqt":
        # UniDock/Meeko records explicit atom indices; do not infer atom order
        # from the original library SMILES (protonation/tautomers can differ).
        lines = text.split("ENDMDL", 1)[0].splitlines()
        mapped_smiles = next((line.removeprefix("REMARK SMILES ") for line in lines
                              if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX ")), "")
        molecule = Chem.MolFromSmiles(mapped_smiles) if mapped_smiles else None
        numbers = [int(value) for line in lines if line.startswith("REMARK SMILES IDX ") for value in line.split()[3:]]
        if molecule is None or not numbers or len(numbers) % 2:
            raise ValueError("The docking pose lacks a valid Meeko SMILES atom map; use a chemically annotated SDF pose.")
        mapping = dict(zip((index - 1 for index in numbers[::2]), numbers[1::2]))
        if len(mapping) != molecule.GetNumAtoms() or len(set(mapping.values())) != len(mapping):
            raise ValueError("The docking pose has an incomplete or ambiguous atom map.")
        coordinates = {int(line[6:11]): tuple(float(line[start:start+8]) for start in (30, 38, 46))
                       for line in lines if line.startswith(("ATOM  ", "HETATM"))}
        if not all(math.isfinite(value) for xyz in coordinates.values() for value in xyz):
            raise ValueError("The docking pose contains non-finite coordinates.")
        elements = {int(line[6:11]): line.split()[-1].upper() for line in lines if line.startswith(("ATOM  ", "HETATM"))}
        from mn_molecular_design_studio.docking_results import _PDBQT_ELEMENTS
        conformer = Chem.Conformer(molecule.GetNumAtoms())
        conformer.Set3D(True)
        for atom in molecule.GetAtoms():
            serial = mapping.get(atom.GetIdx())
            if serial not in coordinates:
                raise ValueError("The docking pose lacks a complete chemistry-to-coordinate atom map.")
            element = _PDBQT_ELEMENTS.get(elements[serial], elements[serial].title())
            if atom.GetSymbol().casefold() != element.casefold():
                raise ValueError("The docking pose atom map disagrees with its atom elements.")
            conformer.SetAtomPosition(atom.GetIdx(), coordinates[serial])
            provenance["atom_mapping"].append({
                "ligand_atom_index": atom.GetIdx(), "pose_atom_serial": serial, "element": element,
            })
        molecule.AddConformer(conformer)
    elif path.suffix.lower() in {".sdf", ".mol"}:
        index = int(context.get("molecule_index") or 0)
        molecules = list(Chem.ForwardSDMolSupplier(str(path), removeHs=False))
        molecule = molecules[index] if index < len(molecules) else None
    else:
        molecule = Chem.MolFromPDBBlock(text, removeHs=False)
        template = Chem.MolFromSmiles(smiles) if smiles else None
        if molecule is not None and template is not None:
            molecule = AllChem.AssignBondOrdersFromTemplate(template, molecule)
    if molecule is None or molecule.GetNumConformers() != 1:
        raise ValueError("This pose does not provide one usable ligand conformation for MD.")
    if path.suffix.lower() != ".pdbqt":
        # For chemically annotated coordinate formats, the source molecule's
        # atom order is the coordinate map. Record that mapping explicitly so
        # the persisted complex can be traced back to the source pose atoms.
        heavy_index = 0
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue
            provenance["atom_mapping"].append({
                "ligand_atom_index": heavy_index,
                "pose_atom_index": atom.GetIdx(),
                "element": atom.GetSymbol(),
            })
            heavy_index += 1
    molecule = Chem.RemoveHs(molecule)
    try:
        Chem.SanitizeMol(molecule)
    except Exception as exc:
        raise ValueError(f"The selected pose does not define valid ligand chemistry: {exc}") from exc
    pose_identity = _canonical_isomeric_smiles(Chem.MolToSmiles(molecule, isomericSmiles=True))
    if pose_identity != source_identity:
        raise ValueError("The selected docking pose chemistry does not match the saved Database ligand identity.")
    provenance["pose_isomeric_smiles"] = pose_identity
    provenance["source_canonical_isomeric_smiles"] = source_identity
    provenance["reconstruction_method"] = (
        "Meeko SMILES/IDX atom-to-coordinate map" if path.suffix.lower() == ".pdbqt"
        else "chemically annotated coordinate structure"
    )
    return (molecule, provenance) if return_provenance else molecule


def _resolve_pose_receptor(database, user, structure):
    """Resolve the exact receptor representation associated with a saved pose."""
    context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
    expected = str(context.get("target_structure_digest") or context.get("target_variant_sha256") or "")
    structure_id = str(context.get("target_structure_id") or "")
    target_molecule_id = str(context.get("target_molecule_id") or "")
    if structure_id:
        if target_molecule_id:
            targets = database.list_molecules(user, molecule_ids=[target_molecule_id], collection="all", limit=1)["items"]
            if not targets or not any(str(item.get("id") or "") == structure_id
                                      for item in targets[0].get("representations") or []):
                raise ValueError("The linked receptor structure is not attached to the saved target molecule.")
        path = database.structure_file(user, structure_id)
        target = next((item for item in (targets[0].get("representations") or [])
                       if str(item.get("id") or "") == structure_id), {}) if target_molecule_id else {}
        target_context = target.get("context") if isinstance(target.get("context"), dict) else {}
        recorded = str(target_context.get("portable_sha256") or target_context.get("target_structure_digest") or "")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if not recorded and not expected:
            raise ValueError("The saved docking pose and linked receptor have no digest to verify their coordinate frame.")
        if recorded and recorded != digest:
            raise ValueError("The saved receptor file no longer matches its Database structure digest.")
        if expected and expected != digest:
            raise ValueError("The saved docking pose and linked receptor have different coordinate-frame digests.")
        return path, digest, structure_id

    target_path = str(context.get("target_path") or "")
    if target_path:
        path = Path(target_path)
        if not path.is_file():
            raise FileNotFoundError("The linked receptor path is unavailable on this workstation.")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if not expected:
            raise ValueError("The saved pose has no receptor digest to verify its coordinate frame.")
        if expected != digest:
            raise ValueError("The saved docking pose and receptor path have different coordinate-frame digests.")
        return path, digest, ""

    if target_molecule_id:
        targets = database.list_molecules(user, molecule_ids=[target_molecule_id], collection="all", limit=1)["items"]
        available = [item for item in (targets[0].get("representations") or [])
                     if str(item.get("format") or "").casefold() in {"pdb", "ent", "cif", "mmcif"}] if targets else []
        if len(available) != 1:
            raise ValueError("The saved pose does not identify one unambiguous receptor structure and coordinate frame.")
        path = database.structure_file(user, available[0]["id"])
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if not expected or expected != digest:
            raise ValueError("The saved pose does not verify the available receptor's coordinate-frame digest.")
        return path, digest, str(available[0]["id"])
    raise ValueError("The saved docking pose has no linked receptor structure or coordinate frame.")


def prepare_bound_inputs(database, user, molecule_id, *, structure_id=None):
    """Stage PDB + matching SDF for the selected evaluation, without starting MD."""
    from rdkit import Chem
    from mn_molecular_design_studio.auth import user_input_root
    from mn_molecular_design_studio.docking_preparations import pdb_text

    records = database.list_molecules(user, molecule_ids=[molecule_id], limit=1, collection="all")["items"]
    if not records:
        raise PermissionError("This molecule is not accessible.")
    record = records[0]
    representations = record.get("representations") or record["structures"]
    structure = next((s for s in (representations if structure_id else record["structures"])
                      if (not structure_id or s["id"] == structure_id) and is_bound_pose(s)), None)
    if structure is None:
        raise ValueError("Select a docking pose with its saved target receptor.")
    pose = database.structure_file(user, structure["id"])
    receptor, receptor_digest, receptor_structure_id = _resolve_pose_receptor(database, user, structure)
    ligand_components = [item for item in record.get("components", [])
                         if str(item.get("type") or "").casefold() == "small_molecule"]
    if len(ligand_components) != 1:
        raise ValueError("The selected Database record must identify exactly one small-molecule ligand for pose reconstruction.")
    smiles = str(ligand_components[0].get("smiles") or "")
    ligand, ligand_provenance = pose_molecule(pose, structure.get("context") or {}, smiles, return_provenance=True)
    receptor_text = pdb_text(receptor.read_text()) if receptor.suffix.lower() in {".cif", ".mmcif"} else receptor.read_text()
    atoms = [line for line in receptor_text.splitlines() if line.startswith(("ATOM  ", "HETATM", "TER   "))]
    if not any(line.startswith("ATOM  ") for line in atoms):
        raise ValueError("The saved receptor has no protein coordinates.")
    chains = {line[21:22] for line in atoms if line.startswith(("ATOM  ", "HETATM"))}
    chain = next((c for c in "LZYXWVUTSRQPONMKJIHGFEDCBA0123456789" if c not in chains), None)
    if not chain:
        raise ValueError("No free PDB chain ID is available for the bound ligand.")
    serial = max((int(line[6:11]) for line in atoms if line.startswith(("ATOM  ", "HETATM"))), default=0)
    if serial + ligand.GetNumAtoms() > 99999:
        raise ValueError("This complex exceeds the PDB atom-number format supported by MD.")
    for atom in ligand.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{atom.GetIdx()+1}".ljust(4)[:4])
        info.SetResidueName("LIG")
        info.SetResidueNumber(1)
        info.SetChainId(chain)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    ligand.SetProp("_Name", record["name"])
    block = Chem.MolToPDBBlock(ligand)
    ligand_lines = []
    for line in block.splitlines():
        if line.startswith("HETATM"):
            ligand_lines.append(f"{line[:6]}{int(line[6:11])+serial:5d}{line[11:]}")
        elif line.startswith("CONECT"):
            numbers = [int(line[start:start+5]) + serial for start in range(6, len(line), 5) if line[start:start+5].strip()]
            ligand_lines.append("CONECT" + "".join(f"{number:5d}" for number in numbers))
    complex_text = "\n".join([*atoms, "TER", *ligand_lines, "END", ""])
    complex_bytes = complex_text.encode("utf-8")
    complex_content_digest = hashlib.sha256(complex_bytes).hexdigest()
    digest = hashlib.sha256((complex_text + Chem.MolToMolBlock(ligand)).encode()).hexdigest()
    folder = user_input_root(user["id"]) / "bound_complexes" / digest
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "complex.pdb").write_text(complex_text)
    with Chem.SDWriter(str(folder / "ligand.sdf")) as writer:
        writer.write(ligand)
    context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
    ligand_selection = {"chain": chain, "residue_key": f"{chain}:1", "ccd": "LIG",
                        "smiles": str(smiles)}
    payload = {"molecule_id": molecule_id, "structure_id": structure["id"], "name": record["name"],
               "complex_pdb": str(folder / "complex.pdb"), "ligand_sdf": str(folder / "ligand.sdf"),
               "ligand_key": f"LIG|{chain}|1|_", "source_pose": structure["id"],
               "target_structure_id": context.get("target_structure_id"),
               "receptor_structure_id": receptor_structure_id,
               "target_structure_digest": receptor_digest,
               "source_pose_digest": hashlib.sha256(pose.read_bytes()).hexdigest(),
               "complex_structure_digest": complex_content_digest,
               "coordinate_frame_id": f"receptor-sha256:{receptor_digest}",
               "ligand_selection": ligand_selection,
               "ligand_isomeric_smiles": str(smiles),
               "ligand_canonical_isomeric_smiles": ligand_provenance["source_canonical_isomeric_smiles"],
               "pose_atom_mapping": ligand_provenance["atom_mapping"],
               "reconstruction_method": ligand_provenance["reconstruction_method"],
               "source_pose_context": context}
    (folder / "provenance.json").write_text(json.dumps(payload, indent=2))
    return payload


def register_bound_complex(database, user, molecule_id, *, structure_id=None):
    """Attach a validated evaluator-compatible complex to its existing molecule."""
    prepared = prepare_bound_inputs(database, user, molecule_id, structure_id=structure_id)
    records = database.list_molecules(user, molecule_ids=[molecule_id], collection="all", limit=1)["items"]
    if not records:
        raise PermissionError("This molecule is not accessible.")
    record = records[0]
    source_pose = next((item for item in (record.get("representations") or record.get("structures") or [])
                        if str(item.get("id") or "") == prepared["source_pose"]), None)
    if source_pose is None:
        raise PermissionError("The source docking pose is no longer attached to this Database record.")
    group_id = str(source_pose.get("group_id") or "")
    database.group(group_id, user, write=True)
    complex_bytes = Path(prepared["complex_pdb"]).read_bytes()
    stored = database._store_file("bound-docking-complex.pdb", complex_bytes)
    context = {
        "database_structure_policy": "bound_docking_complex",
        "target_molecule_id": (prepared["source_pose_context"] or {}).get("target_molecule_id") or "",
        "source_pose_structure_id": prepared["source_pose"],
        "source_pose_digest": prepared["source_pose_digest"],
        "source_receptor_structure_id": prepared["receptor_structure_id"],
        "target_structure_id": prepared["receptor_structure_id"] or prepared.get("target_structure_id") or "",
        "target_structure_digest": prepared["target_structure_digest"],
        "coordinate_frame_id": prepared["coordinate_frame_id"],
        "complex_structure_digest": prepared["complex_structure_digest"],
        "ligand_selection": prepared["ligand_selection"],
        "ligand_isomeric_smiles": prepared["ligand_isomeric_smiles"],
        "ligand_canonical_isomeric_smiles": prepared["ligand_canonical_isomeric_smiles"],
        "pose_atom_mapping": prepared["pose_atom_mapping"],
        "reconstruction_method": prepared["reconstruction_method"],
        "source_docking_provenance": prepared["source_pose_context"],
    }
    structure_id = database.attach_structure(group_id, molecule_id, {
        "path": str(stored), "label": "Docking pose · receptor–ligand complex",
        "dimensions": 3, "portable": True,
        "portable_sha256": prepared["complex_structure_digest"],
        "database_structure_policy": "bound_docking_complex",
    }, context=context)
    return {**prepared, "group_id": group_id, "complex_structure_id": structure_id}
