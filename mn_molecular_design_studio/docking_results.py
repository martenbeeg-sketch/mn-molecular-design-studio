"""Structure helpers for showing protein–ligand docking poses."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import (
    ChainVisualization,
    StructureVisualization,
)


_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASX", "GLX", "SEC", "PYL", "HYP", "MSE",
}
_WATER_RESIDUES = {"HOH", "WAT", "DOD"}
_PDBQT_ELEMENTS = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "O": "O", "OA": "O",
    "S": "S", "SA": "S", "P": "P", "F": "F", "CL": "Cl", "BR": "Br",
    "NS": "N", "OS": "O",
    "I": "I", "HD": "H", "HS": "H", "MG": "Mg", "CA": "Ca", "MN": "Mn",
    "FE": "Fe", "ZN": "Zn", "CU": "Cu", "CO": "Co", "NI": "Ni",
}
_PROTEIN_COLOR = "0x318c70"
_LIGAND_COLOR = "0xe08b3e"


def display_docking_rows(
    rows: list[dict[str, Any]],
    compounds: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Add stable compound identity columns to a docking score table.

    Docking engines intentionally emit their native score schema, which usually
    contains only ``compound_id``.  The staged compound set is the source of
    truth for the display name and SMILES, so this helper enriches the UI/export
    view without rewriting the engine's raw artifact.
    """
    by_id = {
        str(row.get("compound_id") or "").strip(): row
        for row in compounds or []
        if isinstance(row, dict) and str(row.get("compound_id") or "").strip()
    }
    displayed: list[dict[str, Any]] = []
    for raw in rows:
        row = dict(raw)
        compound_id = str(row.get("compound_id") or row.get("Compound ID") or row.get("compound") or "").strip()
        source = by_id.get(compound_id, {})
        source_properties = source.get("properties") if isinstance(source.get("properties"), dict) else {}
        name = str(
            row.get("compound_name") or row.get("Compound name") or source.get("compound_name")
            or source_properties.get("name") or compound_id
        ).strip()
        smiles = str(row.get("smiles") or row.get("SMILES") or source.get("smiles") or "").strip()
        identity = {
            "Compound ID": compound_id or "—",
            "Compound name": name or "—",
            "SMILES": smiles or "—",
        }
        # Keep native score columns, but present the identity fields first and
        # avoid duplicate ID/name/SMILES columns with inconsistent spellings.
        remainder = {
            key: value for key, value in row.items()
            if str(key).casefold() not in {
                "compound_id", "compound", "compound_name", "compound name", "smiles",
            }
        }
        displayed.append({**identity, **remainder})
    return displayed


def compound_identity_map(compounds: list[dict[str, Any]] | None = None) -> dict[str, dict[str, str]]:
    """Build the compound ID-to-name/SMILES map used by pose selectors."""
    return {
        str(row.get("compound_id") or "").strip(): {
            "compound_name": str(row.get("compound_name") or (row.get("properties") or {}).get("name") or row.get("compound_id") or "").strip(),
            "smiles": str(row.get("smiles") or "").strip(),
        }
        for row in compounds or []
        if isinstance(row, dict) and str(row.get("compound_id") or "").strip()
    }


def pdbqt_to_pdb(text: str) -> str:
    """Convert PDBQT atom records to PDB records, retaining only the first pose."""
    lines: list[str] = []
    saw_model = False
    in_first_model = True
    for line in str(text).splitlines():
        if line.startswith("MODEL"):
            if saw_model:
                break
            saw_model = True
            in_first_model = True
            continue
        if line.startswith("ENDMDL"):
            if saw_model:
                break
            continue
        if not in_first_model or not line.startswith(("ATOM  ", "HETATM")):
            continue
        fields = line.split()
        atom_type = fields[-1].upper() if fields else ""
        element = _PDBQT_ELEMENTS.get(atom_type)
        if element is None:
            element = next(
                (symbol for prefix, symbol in (("C", "C"), ("N", "N"), ("O", "O"), ("S", "S"), ("P", "P"))
                 if atom_type.startswith(prefix)),
                atom_type[:1].upper() + atom_type[1:2].lower(),
            )
        lines.append(line[:66].ljust(76) + element.rjust(2))
    return "\n".join([*lines, "END", ""])


def _protein_residue_and_ligand_counts(text: str) -> tuple[int, int]:
    protein_residues: set[tuple[str, str, str]] = set()
    ligand_atoms = 0
    for line in str(text).splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = line[17:20].strip().upper()
        if residue in _AMINO_ACIDS:
            protein_residues.add((line[21:22], line[22:26].strip(), line[26:27]))
        elif residue not in _WATER_RESIDUES:
            ligand_atoms += 1
    return len(protein_residues), ligand_atoms


def _protein_only_pdb(text: str) -> str:
    lines = [
        line for line in str(text).splitlines()
        if line.startswith(("ATOM  ", "TER   "))
        or (line.startswith("HETATM") and line[17:20].strip().upper() in _AMINO_ACIDS)
    ]
    return "\n".join([*lines, "END", ""])


def _ligand_only_pdb(text: str) -> str:
    lines = [
        line for line in str(text).splitlines()
        if line.startswith("HETATM")
        and line[17:20].strip().upper() not in (_AMINO_ACIDS | _WATER_RESIDUES)
    ]
    return "\n".join([*lines, "END", ""])


def _element(line: str) -> str:
    element = line[76:78].strip().upper()
    if element:
        return element
    atom_name = line[12:16].strip().upper().lstrip("0123456789")
    return atom_name[:1]


def _coordinates(line: str) -> tuple[float, float, float] | None:
    try:
        return float(line[30:38]), float(line[38:46]), float(line[46:54])
    except (TypeError, ValueError):
        return None


def _nearby_residue_visualizations(
    protein_text: str,
    ligand_text: str,
    radius_angstrom: float = 5.0,
) -> list[ChainVisualization]:
    """Select whole protein residues near the ligand for a pocket stick overlay."""
    ligand_positions = []
    for line in str(ligand_text).splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = line[17:20].strip().upper()
        if residue in _AMINO_ACIDS or residue in _WATER_RESIDUES or _element(line) in {"H", "D"}:
            continue
        point = _coordinates(line)
        if point is not None:
            ligand_positions.append(point)
    if not ligand_positions:
        return []

    protein_atoms: dict[tuple[str, str, str], list[tuple[float, float, float]]] = {}
    for line in str(protein_text).splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip().upper() not in _AMINO_ACIDS or _element(line) in {"H", "D"}:
            continue
        point = _coordinates(line)
        if point is None:
            continue
        residue_key = (line[21:22], line[22:26].strip(), line[26:27])
        protein_atoms.setdefault(residue_key, []).append(point)

    cutoff_squared = float(radius_angstrom) ** 2
    residues_by_chain: dict[str, set[int]] = {}
    for (chain_id, residue_number, _insertion_code), atom_positions in protein_atoms.items():
        try:
            residue_index = int(residue_number)
        except ValueError:
            continue
        nearby = any(
            (x - ligand_x) ** 2 + (y - ligand_y) ** 2 + (z - ligand_z) ** 2 <= cutoff_squared
            for x, y, z in atom_positions
            for ligand_x, ligand_y, ligand_z in ligand_positions
        )
        if nearby:
            residues_by_chain.setdefault(chain_id, set()).add(residue_index)
    return [
        ChainVisualization(
            chain_id=chain_id,
            residues=sorted(residues),
            color="uniform",
            color_params={"value": _PROTEIN_COLOR},
            representation_type="ball-and-stick",
        )
        for chain_id, residues in sorted(residues_by_chain.items())
    ]


def pose_visualizations(
    pose_text: str,
    pose_name: str,
    receptor_text: str | None = None,
    receptor_name: str = "",
    *,
    include_pocket_sidechains: bool = True,
) -> list[StructureVisualization]:
    """Return Mol* structures for a full complex or a ligand pose plus receptor."""
    suffix = Path(pose_name).suffix.casefold()
    if suffix in {".cif", ".mmcif"}:
        # Use one common PDB path for complex splitting so predicted CIFs render
        # with the same protein/ligand styling as docking poses.
        from mn_molecular_design_studio.docking_preparations import pdb_text
        converted = pdb_text(pose_text)
        if not any(line.startswith(("ATOM  ", "HETATM")) for line in converted.splitlines()):
            raise ValueError("Could not convert this prediction structure for the Database viewer.")
        pose_text = converted
        suffix = ".pdb"

    ligand_text = pdbqt_to_pdb(pose_text) if suffix == ".pdbqt" else pose_text
    protein_residues, ligand_atoms = _protein_residue_and_ligand_counts(ligand_text)
    if suffix != ".pdbqt" and protein_residues >= 5 and ligand_atoms > 0:
        protein_pdb = _protein_only_pdb(pose_text)
        ligand_pdb = _ligand_only_pdb(pose_text)
        nearby_residues = _nearby_residue_visualizations(protein_pdb, ligand_pdb) if include_pocket_sidechains else []
        return [
            StructureVisualization(
                pdb=protein_pdb,
                color="uniform",
                color_params={"value": _PROTEIN_COLOR},
                representation_type="cartoon",
                chains=nearby_residues,
            ),
            StructureVisualization(
                pdb=ligand_pdb,
                color="uniform",
                color_params={"value": _LIGAND_COLOR},
                representation_type="ball-and-stick",
            ),
        ]

    if not receptor_text:
        return [StructureVisualization(pdb=ligand_text, color="uniform", color_params={"value": "0xe08b3e"}, representation_type="ball-and-stick")]

    if Path(receptor_name).suffix.casefold() == ".pdbqt":
        receptor_pdb = pdbqt_to_pdb(receptor_text)
    else:
        receptor_pdb = _protein_only_pdb(receptor_text)
    structures = []
    if any(line.startswith("ATOM  ") for line in receptor_pdb.splitlines()):
        nearby_residues = _nearby_residue_visualizations(receptor_pdb, ligand_text) if include_pocket_sidechains else []
        structures.append(StructureVisualization(
            pdb=receptor_pdb,
            color="uniform",
            color_params={"value": _PROTEIN_COLOR},
            representation_type="cartoon",
            chains=nearby_residues,
        ))
    structures.append(StructureVisualization(
        pdb=ligand_text,
        color="uniform",
        color_params={"value": _LIGAND_COLOR},
        representation_type="ball-and-stick",
    ))
    return structures
