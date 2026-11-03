"""Validated rigid placement of a saved complex in an associated reference frame.

This module fits coordinates only after the caller has an explicit, digest-bound
Database association and a unique exact residue correspondence. It does not
infer biological relationships or calculate accessibility/context scores.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np


PLACEMENT_METHOD = "Bio.PDB.Superimposer Kabsch fit on mapped C-alpha atoms"
PLACEMENT_VERSION = "studio-reference-placement.v1"
COORDINATE_RANK_RELATIVE_TOLERANCE = 1e-6
COORDINATE_RANK_ABSOLUTE_TOLERANCE_ANGSTROM = 1e-7


def fit_ca_atom_superimposer(reference_atoms: list[Any], source_atoms: list[Any], *,
                              require_non_collinear: bool = False):
    """Fit source atoms onto reference atoms with Studio's shared Kabsch fitter.

    The returned Bio.PDB.Superimposer retains its native ``apply`` behavior for
    existing benchmark callers. Placement callers additionally request explicit
    coordinate and degeneracy validation before trusting the transform.
    """
    from Bio.PDB import Superimposer

    if len(reference_atoms) != len(source_atoms):
        raise ValueError("Reference and source C-alpha atom counts do not match.")
    if len(reference_atoms) < 3:
        raise ValueError("At least three corresponding C-alpha atoms are required.")
    fixed = np.asarray([atom.get_coord() for atom in reference_atoms], dtype=float)
    moving = np.asarray([atom.get_coord() for atom in source_atoms], dtype=float)
    if fixed.shape != moving.shape or fixed.ndim != 2 or fixed.shape[1] != 3:
        raise ValueError("Corresponding C-alpha coordinates must have matching N×3 shapes.")
    if not np.isfinite(fixed).all() or not np.isfinite(moving).all():
        raise ValueError("C-alpha coordinates contain non-finite values.")
    if require_non_collinear:
        for label, coordinates in (("reference", fixed), ("source", moving)):
            centered = coordinates - coordinates.mean(axis=0)
            singular_values = np.linalg.svd(centered, compute_uv=False)
            tolerance = max(float(singular_values[0]) * COORDINATE_RANK_RELATIVE_TOLERANCE,
                            COORDINATE_RANK_ABSOLUTE_TOLERANCE_ANGSTROM) if len(singular_values) else COORDINATE_RANK_ABSOLUTE_TOLERANCE_ANGSTROM
            if int(np.count_nonzero(singular_values > tolerance)) < 2:
                raise ValueError(f"The {label} C-alpha coordinates are collinear or geometrically degenerate.")

    superimposer = Superimposer()
    superimposer.set_atoms(reference_atoms, source_atoms)
    rotation, translation = superimposer.rotran
    rotation = np.asarray(rotation, dtype=float)
    translation = np.asarray(translation, dtype=float)
    if not np.isfinite(rotation).all() or not np.isfinite(translation).all():
        raise ValueError("The fitted rigid transformation contains non-finite values.")
    if (not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-5)
            or not math.isclose(float(np.linalg.det(rotation)), 1.0, rel_tol=0.0, abs_tol=1e-5)):
        raise ValueError("The fitted transformation is not a proper rigid rotation.")
    # Mol* uses column vectors while Bio.PDB applies row coordinates as x @ R + t.
    # Transposing the rotation preserves the same source-to-reference direction.
    matrix = np.eye(4, dtype=float)
    matrix[:3, :3] = rotation.T
    matrix[:3, 3] = translation
    return superimposer, matrix.tolist()


def _ca_atom_for_residue(residue, label: str):
    if "CA" not in residue:
        raise ValueError(f"Mapped {label} residue {residue.id[1]}{str(residue.id[2]).strip()} has no C-alpha atom.")
    atom = residue["CA"]
    if getattr(atom, "is_disordered", lambda: 0)() == 2:
        alternatives = list(atom.child_dict.values())
        occupancies = [float(child.get_occupancy() or 0.0) for child in alternatives]
        if not occupancies or not all(math.isfinite(value) and value >= 0.0 for value in occupancies):
            raise ValueError(f"Mapped {label} residue has invalid alternate C-alpha occupancies.")
        maximum = max(occupancies)
        if occupancies.count(maximum) != 1:
            raise ValueError(f"Mapped {label} residue has tied alternate C-alpha atoms; correspondence is ambiguous.")
        atom = alternatives[occupancies.index(maximum)]
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if element and element != "C":
        raise ValueError(f"Mapped {label} residue has a C-alpha atom with non-carbon element {element!r}.")
    occupancy = atom.get_occupancy()
    if occupancy is not None and (not math.isfinite(float(occupancy)) or float(occupancy) <= 0.0):
        raise ValueError(f"Mapped {label} residue has an invalid or zero-occupancy C-alpha atom.")
    coordinates = np.asarray(atom.get_coord(), dtype=float)
    if coordinates.shape != (3,) or not np.isfinite(coordinates).all():
        raise ValueError(f"Mapped {label} residue has non-finite C-alpha coordinates.")
    return atom


def _residue_for_mapping(chain, residue_record: dict, label: str):
    number = int(residue_record["number"])
    insertion = str(residue_record.get("insertion") or " ")[:1]
    residue_id = (" ", number, insertion)
    if residue_id not in chain:
        raise ValueError(f"Mapped {label} residue {number}{insertion.strip()} is absent from the selected chain.")
    residue = chain[residue_id]
    amino_acid = _amino_acid_letter(str(residue.resname).strip().upper())
    if not amino_acid:
        raise ValueError(f"Mapped {label} residue {number}{insertion.strip()} is not a supported amino acid.")
    if amino_acid != str(residue_record.get("amino_acid") or ""):
        raise ValueError(f"Mapped {label} residue identity does not match the saved sequence correspondence.")
    return residue


def _amino_acid_letter(residue_name: str) -> str:
    try:
        from Bio.Data.PDBData import protein_letters_3to1_extended
        return str(protein_letters_3to1_extended.get(residue_name) or "")
    except ImportError:
        return _THREE_TO_ONE.get(residue_name, "")


_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def calculate_reference_placement(source_path: Path | str, reference_path: Path | str, *,
                                  source_chain: str, reference_chain: str,
                                  source_model_id: str, reference_model_id: str,
                                  source_digest: str, reference_digest: str,
                                  mapping: dict[str, Any], source_identity: dict[str, str],
                                  reference_identity: dict[str, str]) -> dict[str, Any]:
    """Fit the explicitly mapped source component and return provenance + matrix."""
    from mn_molecular_design_studio.complex_analysis_comparison import _read_biopython_structure
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import validate_rigid_transform_matrix

    if str(mapping.get("mapping_status") or "") != "validated":
        raise ValueError("Reference placement requires a uniquely validated exact sequence correspondence.")
    if int(mapping.get("alignment_count") or 0) != 1:
        raise ValueError("Reference placement requires exactly one sequence alignment.")
    source_component = mapping.get("source_component") if isinstance(mapping.get("source_component"), dict) else {}
    reference_component = mapping.get("reference_component") if isinstance(mapping.get("reference_component"), dict) else {}
    if str(source_component.get("chain_id") or "") != str(source_chain):
        raise ValueError("The saved source mapping does not match the selected source chain.")
    if str(reference_component.get("chain_id") or "") != str(reference_chain):
        raise ValueError("The saved reference mapping does not match the selected reference chain.")
    source_rows = source_component.get("residues") if isinstance(source_component.get("residues"), list) else []
    reference_rows = reference_component.get("residues") if isinstance(reference_component.get("residues"), list) else []
    residue_pairs = mapping.get("residue_mapping") if isinstance(mapping.get("residue_mapping"), list) else []
    if not source_rows or not reference_rows or len(source_rows) != int(mapping.get("source_residue_count") or -1):
        raise ValueError("The saved sequence correspondence has incomplete source residue provenance.")
    if len(residue_pairs) != len(source_rows) or len({int((item.get("source") or {}).get("sequence_index", -1)) for item in residue_pairs}) != len(source_rows):
        raise ValueError("The saved mapping does not cover each source residue exactly once.")
    if any(not bool(item.get("sequence_match")) for item in residue_pairs):
        raise ValueError("The saved mapping includes non-identical residues; exact-correspondence placement is not supported.")
    if len({int((item.get("reference") or {}).get("sequence_index", -1)) for item in residue_pairs}) != len(residue_pairs):
        raise ValueError("The saved mapping reuses a reference residue.")
    source_sequence = str(source_component.get("sequence") or "")
    reference_sequence = str(reference_component.get("sequence") or "")
    if not source_sequence or len(source_sequence) != len(source_rows) or not reference_sequence or len(reference_sequence) != len(reference_rows):
        raise ValueError("The saved sequence correspondence is missing complete sequence provenance.")
    if len(residue_pairs) != len(source_sequence):
        raise ValueError("Source mapping coverage must be complete before placement.")
    if {int((item.get("source") or {}).get("sequence_index", -1)) for item in residue_pairs} != set(range(len(source_rows))):
        raise ValueError("The saved mapping does not cover the source sequence in a one-to-one correspondence.")
    reference_indices_in_source_order = [
        int(item["reference"]["sequence_index"])
        for item in sorted(residue_pairs, key=lambda row: int(row["source"]["sequence_index"]))
    ]
    if any(right != left + 1 for left, right in zip(reference_indices_in_source_order, reference_indices_in_source_order[1:])):
        raise ValueError("The saved correspondence is not a contiguous exact subsequence; indel or construct-modification mapping is not supported yet.")
    for item in residue_pairs:
        source_index = int(item["source"]["sequence_index"])
        reference_index = int(item["reference"]["sequence_index"])
        if not 0 <= reference_index < len(reference_rows):
            raise ValueError("The saved mapping refers to a reference residue outside its sequence.")
        if item["source"] != source_rows[source_index] or item["reference"] != reference_rows[reference_index]:
            raise ValueError("The saved residue-pair records disagree with the component numbering provenance.")
        if source_sequence[source_index] != reference_sequence[reference_index]:
            raise ValueError("The saved exact mapping contains a residue identity mismatch.")

    source_structure = _read_biopython_structure(Path(source_path))
    reference_structure = _read_biopython_structure(Path(reference_path))
    source_models, reference_models = list(source_structure.get_models()), list(reference_structure.get_models())
    if len(source_models) != 1 or len(reference_models) != 1:
        raise ValueError("Reference placement requires the exact single-model source and reference representations.")
    if str(source_model_id or "").strip() != str(source_identity.get("model_id") or "").strip():
        raise ValueError("The selected source model does not match the saved Database association.")
    if str(reference_model_id or "").strip() != str(reference_identity.get("model_id") or "").strip():
        raise ValueError("The selected reference model does not match the saved Database association.")
    if str(source_digest or "") != str(source_identity.get("structure_digest") or ""):
        raise ValueError("The source structure digest does not match the saved Database association.")
    if str(reference_digest or "") != str(reference_identity.get("structure_digest") or ""):
        raise ValueError("The reference structure digest does not match the saved Database association.")
    source_model, reference_model = source_models[0], reference_models[0]
    if source_chain not in source_model or reference_chain not in reference_model:
        raise ValueError("The selected source or reference chain is absent from its selected model.")
    source_atoms, reference_atoms, mapped_pairs = [], [], []
    for pair in sorted(residue_pairs, key=lambda item: int((item.get("source") or {}).get("sequence_index", -1))):
        source_residue_info = pair.get("source") if isinstance(pair.get("source"), dict) else {}
        reference_residue_info = pair.get("reference") if isinstance(pair.get("reference"), dict) else {}
        source_residue = _residue_for_mapping(source_model[source_chain], source_residue_info, "source")
        reference_residue = _residue_for_mapping(reference_model[reference_chain], reference_residue_info, "reference")
        source_atoms.append(_ca_atom_for_residue(source_residue, "source"))
        reference_atoms.append(_ca_atom_for_residue(reference_residue, "reference"))
        mapped_pairs.append({
            "source": {key: source_residue_info.get(key) for key in ("sequence_index", "number", "insertion", "amino_acid")},
            "reference": {key: reference_residue_info.get(key) for key in ("sequence_index", "number", "insertion", "amino_acid")},
        })

    superimposer, matrix = fit_ca_atom_superimposer(reference_atoms, source_atoms, require_non_collinear=True)
    matrix = validate_rigid_transform_matrix(matrix)
    pair_count = len(source_atoms)
    source_count, reference_count = len(source_rows), len(reference_rows)
    return {
        "schema": "studio.complex_analysis_reference_placement.v1",
        "method": PLACEMENT_METHOD,
        "method_version": PLACEMENT_VERSION,
        "fitter": "Bio.PDB.Superimposer.set_atoms(reference_Calpha, source_Calpha)",
        "fitter_version": _biopython_version(),
        "transformation_direction": "selected complex/source coordinates -> associated reference coordinates",
        "matrix_4x4": matrix,
        "aligned_pair_count": pair_count,
        "source_residue_count": source_count,
        "reference_residue_count": reference_count,
        "source_coverage": pair_count / source_count if source_count else 0.0,
        "reference_coverage": pair_count / reference_count if reference_count else 0.0,
        "rmsd_angstrom": float(superimposer.rms),
        "mapping_method": str(mapping.get("mapping_method") or ""),
        "mapping_parameters": mapping.get("mapping_parameters") or {},
        "mapped_residue_pairs": mapped_pairs,
        "source": dict(source_identity),
        "reference": dict(reference_identity),
        "parameters": {
            "atom_name": "CA", "atom_element": "C", "minimum_pairs": 3,
            "require_non_collinear_coordinates": True,
            "alternate_calpha_policy": "unique highest occupancy; tied occupancy rejected",
            "coordinate_rank_relative_tolerance": COORDINATE_RANK_RELATIVE_TOLERANCE,
            "coordinate_rank_absolute_tolerance_angstrom": COORDINATE_RANK_ABSOLUTE_TOLERANCE_ANGSTROM,
            "mapping_status_required": "validated", "maximum_rmsd_cutoff": None,
            "require_contiguous_reference_subsequence": True,
        },
        "limitations": [
            "The transform places coordinates in the selected reference frame; it does not establish experimental binding or biological accessibility.",
            "RMSD reports conformational difference over mapped C-alpha atoms and has no universal acceptance cutoff.",
            "Unmapped reference residues remain present in the reference structure and are not represented by source coordinates.",
            "No clash, glycan, steric, or full-length accessibility evaluation is calculated.",
        ],
    }


def _biopython_version() -> str:
    try:
        from Bio import __version__
        return str(__version__)
    except ImportError:
        return "unavailable"
