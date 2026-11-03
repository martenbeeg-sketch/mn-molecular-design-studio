"""Reference-frame steric and structural-environment evidence.

This is a static geometry evaluator. It reports atom distances and van-der-Waals
envelope overlaps as evidence; it does not infer binding, accessibility or
biological incompatibility. A placement Evaluation and exact Database identities
are required by the registered workflow before this calculation is called.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from mn_molecular_design_studio.complex_analysis_placement import PLACEMENT_VERSION

METHOD = "Static atom-pair proximity and configured van der Waals envelope overlap"
METHOD_VERSION = "studio.reference-environment.v1"
PROXIMITY_CUTOFF_ANGSTROM = 8.0
OVERLAP_EPSILON_ANGSTROM = 0.4
VDW_RADII_ANGSTROM = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
    "F": 1.47, "P": 1.80, "S": 1.80, "CL": 1.75,
    "BR": 1.85, "I": 1.98, "SE": 1.90, "MG": 1.73, "CA": 2.31,
    "FE": 2.00, "ZN": 1.39, "NA": 2.27, "K": 2.75,
}
_SUGAR_RESIDUES = {
    "NAG", "NDG", "BMA", "MAN", "GLC", "GAL", "FUC", "FUL", "SIA",
    "NAN", "NEU", "XYS", "XYP", "RIB", "ARA", "FRU", "G6D", "GCU",
}


def _read_models(path: Path):
    from mn_molecular_design_studio.complex_analysis_comparison import _read_biopython_structure
    structure = _read_biopython_structure(path)
    models = list(structure.get_models())
    if len(models) != 1:
        raise ValueError("Reference-environment analysis requires the exact single-model representations used by placement.")
    return models[0]


def _atoms(residue, *, chain_id: str, transformed: np.ndarray | None = None,
           unsupported_elements: dict[str, int] | None = None) -> list[dict[str, Any]]:
    result = []
    for atom in residue.get_atoms():
        element = str(getattr(atom, "element", "") or "").strip().upper()
        if not element:
            element = "".join(char for char in str(atom.get_name()) if char.isalpha()).upper()[:1]
        if element in {"H", "D"} or not element:
            continue
        radius = VDW_RADII_ANGSTROM.get(element)
        if radius is None:
            if unsupported_elements is not None:
                unsupported_elements[element] = unsupported_elements.get(element, 0) + 1
            continue
        xyz = np.asarray(atom.get_coord(), dtype=float)
        if transformed is not None:
            xyz = transformed[:3, :3] @ xyz + transformed[:3, 3]
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise ValueError("A selected heavy atom has non-finite coordinates.")
        number = int(residue.id[1])
        insertion = str(residue.id[2] or "").strip()
        result.append({
            "xyz": xyz, "radius": radius, "chain_id": str(chain_id),
            "is_protein_residue": str(residue.resname).strip().upper() in _AA_NAMES,
            "author_number": number, "insertion_code": insertion,
            "author_residue": f"{number}{insertion}",
            "residue_name": str(residue.resname).strip().upper(),
            "atom_name": str(atom.get_name()).strip(), "element": element,
            "occupancy": float(atom.get_occupancy()) if atom.get_occupancy() is not None else None,
            "bfactor": float(atom.get_bfactor()) if atom.get_bfactor() is not None else None,
        })
    return result


def _cif_assembly_membership(path: Path, assembly_id: str, auth_chain: str) -> dict[str, Any]:
    """Confirm only identity-operator assembly membership; never expand operators."""
    if not path.name.casefold().endswith((".cif", ".mmcif")):
        return _pdb_assembly_membership(path, assembly_id, auth_chain)
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    data = MMCIF2Dict(str(path))
    def col(key):
        value = data.get(key, [])
        return list(value) if isinstance(value, (list, tuple)) else ([value] if value else [])
    label_chains = set()
    auths, labels = col("_pdbx_poly_seq_scheme.pdb_strand_id"), col("_pdbx_poly_seq_scheme.asym_id")
    if len(auths) == len(labels):
        for auth, label in zip(auths, labels):
            if str(auth_chain).strip() in {part.strip() for part in str(auth).split(",")}:
                label_chains.add(str(label).strip())
    if not label_chains:
        atom_auths, atom_labels = col("_atom_site.auth_asym_id"), col("_atom_site.label_asym_id")
        if len(atom_auths) == len(atom_labels):
            label_chains.update(str(label).strip() for auth, label in zip(atom_auths, atom_labels)
                                if str(auth).strip() == str(auth_chain).strip())
        if not label_chains:
            label_chains.add(str(auth_chain).strip())
    gen_assembly = col("_pdbx_struct_assembly_gen.assembly_id")
    expressions = col("_pdbx_struct_assembly_gen.oper_expression")
    asym_lists = col("_pdbx_struct_assembly_gen.asym_id_list")
    operators = col("_pdbx_struct_oper_list.id")
    matrices = [[col(f"_pdbx_struct_oper_list.matrix[{r}][{c}]") for c in range(1, 4)] for r in range(1, 4)]
    vectors = [col(f"_pdbx_struct_oper_list.vector[{r}]") for r in range(1, 4)]
    found = False
    for aid, expr, asym in zip(gen_assembly, expressions, asym_lists):
        chains = {part.strip() for part in str(asym).split(",")}
        if str(aid) != str(assembly_id) or not chains.intersection(label_chains):
            continue
        found = True
        expression = str(expr).strip().strip("()")
        op_ids = [item.strip() for item in expression.split(",")]
        if len(op_ids) != 1 or not op_ids[0].isdigit():
            return {"status": "unsupported_operator_expression", "assembly_id": assembly_id}
        try:
            index = operators.index(op_ids[0])
            matrix = np.array([[float(matrices[r][c][index]) for c in range(3)] for r in range(3)])
            vector = np.array([float(vectors[r][index]) for r in range(3)])
        except (ValueError, IndexError, TypeError):
            return {"status": "assembly_operator_unavailable", "assembly_id": assembly_id}
        if not np.allclose(matrix, np.eye(3), atol=1e-6) or not np.allclose(vector, 0, atol=1e-6):
            return {"status": "non_identity_assembly_operator_not_supported", "assembly_id": assembly_id}
    if found:
        return {"status": "confirmed_identity_operator_member", "assembly_id": str(assembly_id)}
    return {"status": "not_listed_in_selected_assembly", "assembly_id": str(assembly_id)}


def _pdb_assembly_membership(path: Path, assembly_id: str, auth_chain: str) -> dict[str, Any]:
    """Read PDB REMARK 350 membership and BIOMT, accepting identity only."""
    import re
    text = path.read_text(encoding="utf-8", errors="replace")
    current_assembly = ""
    chains: set[str] = set()
    operators: dict[str, dict[int, tuple[list[float], float]]] = {}
    assembly_chains: dict[str, set[str]] = {}
    assembly_operators: dict[str, dict[str, dict[int, tuple[list[float], float]]]] = {}
    for line in text.splitlines():
        if not line.startswith("REMARK 350"):
            continue
        if "BIOMOLECULE:" in line:
            current_assembly = line.split("BIOMOLECULE:", 1)[1].strip().split(",", 1)[0]
            chains = set()
            operators = {}
            assembly_chains.setdefault(current_assembly, set())
            assembly_operators.setdefault(current_assembly, {})
        elif "APPLY THE FOLLOWING TO CHAINS:" in line or "AND CHAINS:" in line:
            value = line.split(":", 1)[1]
            chains.update(part.strip().strip(".") for part in value.split(",") if part.strip().strip(".") not in {"NULL", ""})
            assembly_chains.setdefault(current_assembly, set()).update(chains)
        else:
            match = re.search(r"BIOMT([123])\s+(\d+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)", line)
            if match:
                row, op = int(match.group(1)), match.group(2)
                operators.setdefault(op, {})[row] = ([float(match.group(i)) for i in (3, 4, 5)], float(match.group(6)))
                assembly_operators.setdefault(current_assembly, {})[op] = operators[op]
    if str(auth_chain).strip() not in assembly_chains.get(str(assembly_id), set()):
        return {"status": "not_listed_in_selected_assembly", "assembly_id": str(assembly_id)}
    operators = assembly_operators.get(str(assembly_id), {})
    if not operators:
        return {"status": "assembly_operator_unavailable", "assembly_id": str(assembly_id)}
    for rows in operators.values():
        if set(rows) != {1, 2, 3}:
            return {"status": "assembly_operator_unavailable", "assembly_id": str(assembly_id)}
        matrix = np.asarray([rows[index][0] for index in (1, 2, 3)], dtype=float)
        vector = np.asarray([rows[index][1] for index in (1, 2, 3)], dtype=float)
        if not np.allclose(matrix, np.eye(3), atol=1e-6) or not np.allclose(vector, 0, atol=1e-6):
            return {"status": "non_identity_assembly_operator_not_supported", "assembly_id": str(assembly_id)}
    return {"status": "confirmed_identity_operator_member", "assembly_id": str(assembly_id)}


def _reference_residue_groups(model, *, omitted_chain: str, mapped_reference_residues: set[tuple[int, str]],
                              selected_chains: list[str], reference_path: Path,
                              assembly_id: str | None,
                              unsupported_elements: dict[str, int] | None = None) -> list[dict[str, Any]]:
    if omitted_chain not in model:
        raise ValueError("The mapped reference chain is absent from the selected reference model.")
    selected = set(str(x) for x in selected_chains)
    if omitted_chain in selected:
        selected.remove(omitted_chain)
    groups: list[dict[str, Any]] = []
    assembly_status = {}
    if assembly_id and selected:
        assembly_status[omitted_chain] = _cif_assembly_membership(reference_path, assembly_id, omitted_chain)
        if assembly_status[omitted_chain]["status"] == "non_identity_assembly_operator_not_supported":
            raise ValueError("The associated reference chain requires a non-identity operator; assembly expansion is not implemented, so no placement-context geometry was calculated.")
        if assembly_status[omitted_chain]["status"] != "confirmed_identity_operator_member":
            raise ValueError(f"The associated reference chain is not confirmed as an identity-operator member of assembly {assembly_id!r}.")
        for chain_id in selected:
            if chain_id not in model:
                raise ValueError(f"Selected reference component chain {chain_id!r} is absent from the associated reference model.")
            assembly_status[chain_id] = _cif_assembly_membership(reference_path, assembly_id, chain_id)
            if assembly_status[chain_id]["status"] == "non_identity_assembly_operator_not_supported":
                raise ValueError("The selected biological-assembly component requires a non-identity operator; assembly expansion is not implemented, so no placement-context geometry was calculated.")
            if assembly_status[chain_id]["status"] != "confirmed_identity_operator_member":
                raise ValueError(f"Selected chain {chain_id!r} is not confirmed as an identity-operator member of assembly {assembly_id!r}.")
    branched_components_by_chain: dict[str, set[str]] = {}
    if reference_path.name.casefold().endswith((".cif", ".mmcif")):
        try:
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict
            data = MMCIF2Dict(str(reference_path))
            def col(name):
                value = data.get(name, [])
                return list(value) if isinstance(value, (list, tuple)) else ([value] if value else [])
            branch_entities = col("_pdbx_entity_branch_list.entity_id")
            branch_components = col("_pdbx_entity_branch_list.comp_id")
            asym_ids, entity_ids = col("_struct_asym.id"), col("_struct_asym.entity_id")
            auth_asym, label_asym = col("_atom_site.auth_asym_id"), col("_atom_site.label_asym_id")
            author_by_label = {}
            for author, label in zip(auth_asym, label_asym):
                author_by_label.setdefault(str(label), str(author))
            branch_names_by_entity: dict[str, set[str]] = {}
            for entity, component in zip(branch_entities, branch_components):
                branch_names_by_entity.setdefault(str(entity), set()).add(str(component).strip().upper())
            for asym, entity in zip(asym_ids, entity_ids):
                if str(entity) in branch_names_by_entity:
                    author_chain = author_by_label.get(str(asym), str(asym))
                    branched_components_by_chain.setdefault(author_chain, set()).update(branch_names_by_entity[str(entity)])
        except (OSError, ValueError, TypeError):
            branched_components_by_chain = {}
    for chain in model:
        chain_id = str(chain.id)
        if chain_id != omitted_chain and chain_id not in selected:
            continue
        for residue in chain:
            hetero, number, insertion = residue.id
            key = (int(number), str(insertion or "").strip())
            is_protein = str(residue.resname).strip().upper() in _AA_NAMES
            is_sugar = (str(residue.resname).strip().upper() in _SUGAR_RESIDUES
                        or str(residue.resname).strip().upper() in branched_components_by_chain.get(chain_id, set()))
            if is_sugar:
                kind = "explicitly_modeled_sugar_residue"
            elif chain_id == omitted_chain:
                if not is_protein or key in mapped_reference_residues:
                    continue
                kind = "omitted_reference_region"
            elif chain_id not in selected:
                continue
            else:
                kind = ("confirmed_assembly_component" if assembly_status.get(chain_id, {}).get("status")
                        == "confirmed_identity_operator_member" else "selected_reference_component")
            atoms = _atoms(residue, chain_id=chain_id, unsupported_elements=unsupported_elements)
            if atoms:
                groups.append({"kind": kind, "chain_id": chain_id, "residue": residue, "atoms": atoms,
                               "assembly_membership": assembly_status.get(chain_id, {"status": "not_assessed"})})
    for chain_id in selected:
        if chain_id not in model:
            raise ValueError(f"Selected reference component chain {chain_id!r} is absent from the associated reference model.")
    return groups


_AA_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL",
}


def calculate_reference_environment(source_path: Path | str, reference_path: Path | str, *,
                                    placement: dict[str, Any], source_identity: dict[str, Any],
                                    reference_identity: dict[str, Any], reference_context_id: str,
                                    selected_reference_chains: list[str] | None = None,
                                    assembly_id: str | None = None) -> dict[str, Any]:
    """Measure static atom proximity from the placed complex to selected context."""
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import validate_rigid_transform_matrix
    from mn_molecular_design_studio.evaluations import file_digest
    source_path, reference_path = Path(source_path), Path(reference_path)
    if placement.get("method_version") != PLACEMENT_VERSION:
        raise ValueError("The saved placement uses an unsupported method version.")
    if not str(placement.get("mapping_method") or "").strip():
        raise ValueError("The saved placement does not record its residue-mapping method.")
    if placement.get("transformation_direction") != "selected complex/source coordinates -> associated reference coordinates":
        raise ValueError("The saved placement transformation direction is unsupported.")
    if str(placement.get("reference_context_id") or "") != str(reference_context_id):
        raise ValueError("The saved placement is not bound to this exact reference association.")
    if assembly_id and not selected_reference_chains:
        raise ValueError("An assembly ID is meaningful only when additional reference components are explicitly selected.")
    if placement.get("source") != source_identity or placement.get("reference") != reference_identity:
        raise ValueError("The saved placement identities do not exactly match the selected Database structures.")
    if file_digest(source_path) != str(source_identity.get("structure_digest") or ""):
        raise ValueError("The selected source structure digest is stale.")
    if file_digest(reference_path) != str(reference_identity.get("structure_digest") or ""):
        raise ValueError("The associated reference structure digest is stale.")
    matrix = np.asarray(validate_rigid_transform_matrix(placement.get("matrix_4x4")), dtype=float)
    if matrix.shape != (4, 4):
        raise ValueError("The saved placement matrix is invalid.")
    source_model, reference_model = _read_models(source_path), _read_models(reference_path)
    source_atoms = []
    unsupported_elements: dict[str, int] = {}
    for chain in source_model:
        for residue in chain:
            source_atoms.extend(_atoms(residue, chain_id=str(chain.id), transformed=matrix,
                                       unsupported_elements=unsupported_elements))
    if not source_atoms:
        raise ValueError("The selected complex has no supported heavy-atom coordinates.")
    mapping = placement.get("mapped_residue_pairs")
    if not isinstance(mapping, list) or len(mapping) < 3:
        raise ValueError("The saved placement has insufficient validated residue correspondence.")
    source_keys, reference_keys = [], []
    for pair in mapping:
        if not isinstance(pair, dict) or not isinstance(pair.get("source"), dict) or not isinstance(pair.get("reference"), dict):
            raise ValueError("The saved placement residue correspondence is malformed.")
        try:
            source_keys.append((int(pair["source"]["number"]), str(pair["source"].get("insertion") or "").strip()))
            reference_keys.append((int(pair["reference"]["number"]), str(pair["reference"].get("insertion") or "").strip()))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("The saved placement residue correspondence is malformed.") from exc
    if len(set(source_keys)) != len(source_keys) or len(set(reference_keys)) != len(reference_keys):
        raise ValueError("The saved placement residue correspondence is ambiguous because residue keys are repeated.")
    if int(placement.get("aligned_pair_count") or 0) != len(mapping):
        raise ValueError("The saved placement pair count does not match its residue correspondence.")
    mapped_reference = {
        key for key in reference_keys
    }
    reference_chain = str(reference_identity.get("chain_id") or "")
    selected_reference_chains = list(dict.fromkeys(str(chain) for chain in (selected_reference_chains or [])))
    if reference_chain in selected_reference_chains:
        raise ValueError("The associated reference chain is already assessed for omitted residues and cannot also be selected as an additional component.")
    if assembly_id and not selected_reference_chains:
        raise ValueError("An assembly ID requires at least one explicitly selected additional reference component.")
    context_groups = _reference_residue_groups(
        reference_model, omitted_chain=reference_chain, mapped_reference_residues=mapped_reference,
        selected_chains=selected_reference_chains, reference_path=reference_path,
        assembly_id=assembly_id, unsupported_elements=unsupported_elements,
    )
    if not context_groups:
        raise ValueError("No resolved omitted-reference or selected context heavy atoms are available for this evaluation.")
    source_xyz = np.asarray([atom["xyz"] for atom in source_atoms], dtype=float)
    tree = cKDTree(source_xyz)
    source_component_chain = str(source_identity.get("chain_id") or "")
    role_indexes = {
        "associated_source_component": [index for index, atom in enumerate(source_atoms)
                                        if atom["chain_id"] == source_component_chain and atom["is_protein_residue"]],
        "other_source_complex_components": [index for index, atom in enumerate(source_atoms)
                                             if atom["chain_id"] != source_component_chain or not atom["is_protein_residue"]],
    }
    if not role_indexes["associated_source_component"]:
        raise ValueError("The associated source component has no supported heavy atoms in the selected structure.")
    role_atoms = {role: [source_atoms[index] for index in indexes] for role, indexes in role_indexes.items()}
    role_trees = {role: cKDTree(np.asarray([atom["xyz"] for atom in atoms], dtype=float))
                  for role, atoms in role_atoms.items() if atoms}
    rows, closest = [], []
    source_by_index = source_atoms
    for group in context_groups:
        context_atoms = group["atoms"]
        context_xyz = np.asarray([atom["xyz"] for atom in context_atoms], dtype=float)
        pairs = tree.query_ball_point(context_xyz, PROXIMITY_CUTOFF_ANGSTROM)
        pair_rows = []
        min_pair = None
        for context_index, source_indices in enumerate(pairs):
            context_atom = context_atoms[context_index]
            for source_index in source_indices:
                source_atom = source_by_index[source_index]
                distance = float(np.linalg.norm(source_atom["xyz"] - context_atom["xyz"]))
                overlap = float(source_atom["radius"] + context_atom["radius"] - distance)
                pair_rows.append((distance, overlap, source_atom, context_atom))
                if min_pair is None or distance < min_pair[0]:
                    min_pair = (distance, overlap, source_atom, context_atom)
        if min_pair is None:
            # Distinguish no pair within the reporting cutoff from no geometry.
            distances, indices = tree.query(context_xyz, k=1)
            index = int(np.argmin(distances))
            source_atom = source_by_index[int(indices[index])]
            context_atom = context_atoms[index]
            min_pair = (float(distances[index]), float(source_atom["radius"] + context_atom["radius"] - distances[index]), source_atom, context_atom)
        def evidence(pair):
            distance, overlap, source_atom, context_atom = pair
            return {
                "distance_angstrom": round(float(distance), 4),
                "vdw_envelope_overlap_angstrom": round(float(overlap), 4),
                "source_role": ("associated_source_component" if source_atom["chain_id"] == source_component_chain and source_atom["is_protein_residue"]
                                else "other_source_complex_component"),
                "source": {key: source_atom[key] for key in ("chain_id", "author_number", "insertion_code", "author_residue", "residue_name", "atom_name", "element")},
                "context": {key: context_atom[key] for key in ("chain_id", "author_number", "insertion_code", "author_residue", "residue_name", "atom_name", "element")},
            }
        close = [item for item in pair_rows if item[0] <= PROXIMITY_CUTOFF_ANGSTROM]
        overlaps = [item for item in close if item[1] > OVERLAP_EPSILON_ANGSTROM]
        role_evidence = {}
        for role, indexes in role_indexes.items():
            role_pairs = [item for item in close if item[2]["chain_id"] == source_component_chain and item[2]["is_protein_residue"]] if role == "associated_source_component" else [item for item in close if item[2]["chain_id"] != source_component_chain or not item[2]["is_protein_residue"]]
            role_overlaps = [item for item in role_pairs if item[1] > OVERLAP_EPSILON_ANGSTROM]
            role_nearest = None
            if role_pairs:
                role_nearest = min(role_pairs, key=lambda item: item[0])
            elif role in role_trees:
                distances, indices = role_trees[role].query(context_xyz, k=1)
                context_index = int(np.argmin(distances))
                local_index = int(indices[context_index])
                context_atom = context_atoms[context_index]
                source_atom = role_atoms[role][local_index]
                distance = float(distances[context_index])
                role_nearest = (distance, float(source_atom["radius"] + context_atom["radius"] - distance), source_atom, context_atom)
            role_evidence[role] = {
                "closest_distance_angstrom": round(float(role_nearest[0]), 4) if role_nearest else None,
                "atom_pairs_within_proximity_cutoff": len(role_pairs),
                "atom_pairs_with_vdw_overlap_above_epsilon": len(role_overlaps),
                "closest_atom_pair": evidence(role_nearest) if role_nearest else None,
                "closest_overlap_pair": evidence(min(role_overlaps, key=lambda item: item[0])) if role_overlaps else None,
            }
        rows.append({
            "context_kind": group["kind"], "chain_id": group["chain_id"],
            "author_number": int(group["residue"].id[1]),
            "insertion_code": str(group["residue"].id[2] or "").strip(),
            "residue_name": str(group["residue"].resname).strip().upper(),
            "assembly_membership": group.get("assembly_membership", {"status": "not_assessed"}),
            "closest_distance_angstrom": round(min_pair[0], 4),
            "atom_pairs_within_proximity_cutoff": len(close),
            "atom_pairs_with_vdw_overlap_above_epsilon": len(overlaps),
            "closest_atom_pair": evidence(min_pair),
            "closest_overlap_pair": evidence(min(overlaps, key=lambda item: item[0])) if overlaps else None,
            "source_component_geometry": role_evidence["associated_source_component"],
            "other_source_complex_geometry": role_evidence["other_source_complex_components"],
        })
        closest.append((min_pair[0], evidence(min_pair)))
    rows.sort(key=lambda row: (row["closest_distance_angstrom"], row["context_kind"], row["chain_id"], row["author_number"]))
    closest.sort(key=lambda item: item[0])
    reference_origin = _coordinate_origin(reference_path)
    sequence_sequons = _sequence_only_sequons(reference_path, reference_chain, mapping)
    limitations = [
        "Static coordinate geometry only; positive van-der-Waals envelope overlap is a potential overlap, not an experimental clash or binding-quality score.",
        "Reference conformational uncertainty, source-model uncertainty, alternate conformers, solvent, flexibility and missing assembly operators can change these measurements.",
        "Selected additional chains are not asserted to be biological partners unless explicitly confirmed in the selected assembly metadata.",
        "Sequence-only N-X-S/T motifs are not measured glycans; explicit modeled sugar coordinates are reported separately by residue identity.",
        "Sugar residue identity is reported from recognized component codes; covalent glycan attachment/linkage is not validated by this evaluator.",
        "No glycan occupancy, glycan conformational ensemble, solvent accessibility or biological accessibility is inferred.",
    ]
    if unsupported_elements:
        limitations.append("Heavy atoms with unsupported configured radii were excluded from atom-pair measurements; counts by element are recorded in parameters.")
    return {
        "schema": "studio.complex_analysis_reference_environment.v1",
        "method": METHOD, "method_version": METHOD_VERSION,
        "reference_context_id": str(reference_context_id),
        "source": dict(source_identity), "reference": dict(reference_identity),
        "reference_coordinate_origin": reference_origin,
        "placement": {
            "method_version": placement.get("method_version"),
            "mapping_method": str(placement.get("mapping_method") or ""),
            "evaluation_id": placement.get("evaluation_id"),
            "run_id": placement.get("run_id"),
            "aligned_pair_count": int(placement.get("aligned_pair_count") or 0),
            "source_coverage": float(placement.get("source_coverage") or 0),
            "reference_coverage": float(placement.get("reference_coverage") or 0),
            "rmsd_angstrom": float(placement.get("rmsd_angstrom") or 0),
            "transformation_direction": placement.get("transformation_direction"),
        },
        "parameters": {
            "heavy_atoms_only": True, "proximity_cutoff_angstrom": PROXIMITY_CUTOFF_ANGSTROM,
            "vdw_radii": "Configured element-specific radii; Bondi values are used for H, C, N, O, F, P, S, Cl, Br, I, and Se; remaining listed elements use explicit configured values",
            "vdw_radii_angstrom": dict(VDW_RADII_ANGSTROM),
            "vdw_radius_sources": {
                "bondi_1964_elements": ["H", "C", "N", "O", "F", "P", "S", "CL", "BR", "I", "SE"],
                "studio_configured_elements": [element for element in VDW_RADII_ANGSTROM
                                                if element not in {"H", "C", "N", "O", "F", "P", "S", "CL", "BR", "I", "SE"}],
            },
            "overlap_epsilon_angstrom": OVERLAP_EPSILON_ANGSTROM,
            "selected_reference_chains": list(selected_reference_chains or []),
            "assembly_id": assembly_id, "assembly_operator_policy": "identity operators only; non-identity expansion is unsupported",
            "residue_correspondence_method": str(placement.get("mapping_method") or ""),
            "sugar_residue_identification": {
                "recognized_residue_names": sorted(_SUGAR_RESIDUES),
                "mmcif_branch_components": "also recognized by exact branch entity component IDs when available",
                "covalent_attachment_validated": False,
            },
            "alternate_conformers": "Bio.PDB selected conformer; unresolved alternatives are not separately enumerated",
            "unsupported_heavy_atom_elements": dict(sorted(unsupported_elements.items())),
        },
        "measurements": {
            "context_residue_count": len(rows),
            "residue_count_with_atom_pairs_within_proximity_cutoff": sum(bool(row["atom_pairs_within_proximity_cutoff"]) for row in rows),
            "residue_count_with_vdw_overlap_above_epsilon": sum(bool(row["atom_pairs_with_vdw_overlap_above_epsilon"]) for row in rows),
            "closest_distance_angstrom": rows[0]["closest_distance_angstrom"] if rows else None,
            "atom_pairs_within_proximity_cutoff": sum(row["atom_pairs_within_proximity_cutoff"] for row in rows),
            "atom_pairs_with_vdw_overlap_above_epsilon": sum(row["atom_pairs_with_vdw_overlap_above_epsilon"] for row in rows),
            "associated_source_component_residues_with_pairs_within_proximity_cutoff": sum(bool(row["source_component_geometry"]["atom_pairs_within_proximity_cutoff"]) for row in rows),
            "associated_source_component_closest_distance_angstrom": min((row["source_component_geometry"]["closest_distance_angstrom"] for row in rows if row["source_component_geometry"]["closest_distance_angstrom"] is not None), default=None),
            "associated_source_component_atom_pairs_within_proximity_cutoff": sum(row["source_component_geometry"]["atom_pairs_within_proximity_cutoff"] for row in rows),
            "associated_source_component_residues_with_vdw_overlap_above_epsilon": sum(bool(row["source_component_geometry"]["atom_pairs_with_vdw_overlap_above_epsilon"]) for row in rows),
            "associated_source_component_atom_pairs_with_vdw_overlap_above_epsilon": sum(row["source_component_geometry"]["atom_pairs_with_vdw_overlap_above_epsilon"] for row in rows),
            "other_source_complex_residues_with_pairs_within_proximity_cutoff": sum(bool(row["other_source_complex_geometry"]["atom_pairs_within_proximity_cutoff"]) for row in rows),
            "other_source_complex_closest_distance_angstrom": min((row["other_source_complex_geometry"]["closest_distance_angstrom"] for row in rows if row["other_source_complex_geometry"]["closest_distance_angstrom"] is not None), default=None),
            "other_source_complex_atom_pairs_within_proximity_cutoff": sum(row["other_source_complex_geometry"]["atom_pairs_within_proximity_cutoff"] for row in rows),
            "other_source_complex_residues_with_vdw_overlap_above_epsilon": sum(bool(row["other_source_complex_geometry"]["atom_pairs_with_vdw_overlap_above_epsilon"]) for row in rows),
            "other_source_complex_atom_pairs_with_vdw_overlap_above_epsilon": sum(row["other_source_complex_geometry"]["atom_pairs_with_vdw_overlap_above_epsilon"] for row in rows),
            "unsupported_heavy_atom_count": sum(unsupported_elements.values()),
        },
        "residue_evidence": rows,
        "closest_atom_pairs": [item[1] for item in closest[:100]],
        "sequence_only_potential_n_linked_sequons": sequence_sequons,
        "limitations": limitations,
    }


def _coordinate_origin(path: Path) -> dict[str, Any]:
    """Report only experimental/predicted labels explicit in the coordinate file."""
    if path.name.casefold().endswith((".cif", ".mmcif")):
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        data = MMCIF2Dict(str(path))
        method = data.get("_exptl.method", [])
        methods = list(method) if isinstance(method, (list, tuple)) else ([method] if method else [])
        if methods:
            classification = "predicted" if any(
                token in str(value).casefold() for value in methods for token in ("theoretical", "predict", "model")
            ) else "experimental"
            return {"classification": classification, "file_reported_method": [str(value) for value in methods]}
        model_type = data.get("_ma_model_list.model_type", [])
        models = list(model_type) if isinstance(model_type, (list, tuple)) else ([model_type] if model_type else [])
        if any("pred" in str(value).casefold() for value in models):
            return {"classification": "predicted", "file_reported_method": [str(value) for value in models]}
    else:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("EXPDTA"):
                value = line[10:].strip()
                if "THEORETICAL MODEL" in value.upper() or "PREDICT" in value.upper():
                    return {"classification": "predicted", "file_reported_method": [value]}
                return {"classification": "experimental", "file_reported_method": [value]}
    return {"classification": "unknown", "file_reported_method": []}


def _sequence_only_sequons(path: Path, chain_id: str, residue_mapping: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """List N-X-S/T sequence motifs without inferring glycan occupancy or geometry."""
    from mn_molecular_design_studio.complex_analysis_reference_context import reference_component_track
    track = reference_component_track(path, chain_id, residue_mapping)
    sequence = str(track.get("sequence") or "")
    rows = list(track.get("residues") or [])
    result = []
    for index in range(max(0, len(sequence) - 2)):
        motif = sequence[index:index + 3]
        if len(motif) == 3 and motif[0] == "N" and motif[1] != "P" and motif[2] in {"S", "T"}:
            residue = rows[index] if index < len(rows) else {}
            result.append({
                "kind": "sequence_only_potential_n_linked_sequon",
                "sequence_motif": motif,
                "sequence_start_one_based": index + 1,
                "asn_author_residue": residue.get("author_residue"),
                "asn_author_number": residue.get("author_number"),
                "insertion_code": residue.get("insertion_code") or "",
                "sequence_only": True,
                "occupancy_or_glycan_geometry_measured": False,
            })
    return result
