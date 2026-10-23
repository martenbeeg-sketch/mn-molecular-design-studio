"""Small adapters around the vendored SAbR antibody-numbering package."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio.PDB import MMCIFIO, PDBIO

from mn_molecular_design_studio.design_validation import read_structure
from mn_molecular_design_studio.vendor.sabr import cdr_residue_records, number_chain, renumber_structure
from mn_molecular_design_studio.vendor.sabr import constants as sabr_constants


CDR_NAMES = tuple(sabr_constants.IMGT_LOOPS)


def residue_label(chain_id: str, residue) -> str:
    """Return the author-style residue label used by Studio's MPNN masks."""
    number = int(residue.id[1])
    insertion = str(residue.id[2] or "").strip()
    return f"{chain_id}{number}{insertion}"


def annotate_chain(
    structure,
    chain_id: str,
    *,
    chain_type: str = "auto",
    noise_level: float = 0.0,
    mode: str = "sabr",
    scheme: str = "imgt",
) -> dict[str, Any]:
    """Number one chain and map IMGT CDRs back to input-structure residue IDs."""
    model = next(structure.get_models())
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id!r} is not present in the structure.")
    chain = model[chain_id]
    residue_rows = list(chain.get_list())
    cdr_records = cdr_residue_records(
        structure, chain_id, chain_type=chain_type, noise_level=noise_level, mode=mode
    )
    cdr_by_index: dict[int, str] = {}
    imgt_by_index: dict[int, int] = {}
    cdr_labels: dict[str, list[str]] = {name: [] for name in CDR_NAMES}
    for region, records in cdr_records.items():
        for index, imgt_position, insertion, amino_acid in records:
            if index < 0 or index >= len(residue_rows):
                raise ValueError("SAbR returned a residue index outside the input chain.")
            residue = residue_rows[index]
            cdr_by_index[index] = region
            imgt_by_index[index] = int(imgt_position)
            cdr_labels[region].append(residue_label(chain_id, residue))

    numbering = number_chain(
        structure, chain_id, scheme=scheme, chain_type=chain_type,
        noise_level=noise_level, mode=mode,
    )
    numbered_by_index = {
        index: {"number": int(number), "insertion_code": str(insertion or "").strip(), "amino_acid": aa}
        for index, number, insertion, aa in numbering
    }
    residues = []
    for index, residue in enumerate(residue_rows):
        label = residue_label(chain_id, residue)
        number = numbered_by_index.get(index)
        residues.append({
            "input_index": index,
            "input_residue": label,
            "residue_name": str(residue.resname),
            "cdr": cdr_by_index.get(index),
            "imgt_position": imgt_by_index.get(index),
            "numbered_position": number["number"] if number else None,
            "numbered_insertion_code": number["insertion_code"] if number else "",
            "amino_acid": number["amino_acid"] if number else "X",
        })
    return {
        "chain_id": chain_id,
        "chain_type": chain_type,
        "scheme": scheme,
        "mode": mode,
        "noise_level": float(noise_level),
        "cdr_residues": cdr_labels,
        "residues": residues,
        "numbered_residue_count": len(numbering),
    }


def annotate_file(
    path: Path,
    chain_ids: list[str] | None = None,
    *,
    chain_type: str = "auto",
    noise_level: float = 0.0,
    mode: str = "sabr",
    scheme: str = "imgt",
) -> dict[str, Any]:
    structure = read_structure(Path(path))
    model = next(structure.get_models())
    selected = list(chain_ids or [chain.id for chain in model])
    return {
        "structure": str(Path(path).resolve()),
        "chains": {
            chain_id: annotate_chain(
                structure, chain_id, chain_type=chain_type, noise_level=noise_level,
                mode=mode, scheme=scheme,
            )
            for chain_id in selected
        },
    }


def cdr_fixed_residues(
    path: Path,
    chain_ids: list[str],
    cdr_regions: list[str],
    *,
    chain_type: str = "auto",
    noise_level: float = 0.0,
    mode: str = "sabr",
) -> tuple[list[str], dict[str, list[str]]]:
    """Return all non-selected-CDR residues as fixed plus selected CDR labels."""
    unknown = set(cdr_regions) - set(CDR_NAMES)
    if unknown:
        raise ValueError(f"Unknown CDR region(s): {', '.join(sorted(unknown))}.")
    if not cdr_regions:
        raise ValueError("Choose at least one CDR to redesign.")
    structure = read_structure(Path(path))
    model = next(structure.get_models())
    fixed: list[str] = []
    selected_by_region = {name: [] for name in CDR_NAMES}
    for chain_id in chain_ids:
        if chain_id not in model:
            raise ValueError(f"Binder chain {chain_id!r} is not present in {Path(path).name}.")
        annotation = annotate_chain(
            structure, chain_id, chain_type=chain_type, noise_level=noise_level, mode=mode
        )
        chosen = {
            label
            for region in cdr_regions
            for label in annotation["cdr_residues"][region]
        }
        for region in CDR_NAMES:
            selected_by_region[region].extend(
                label for label in annotation["cdr_residues"][region] if region in cdr_regions
            )
        chain = model[chain_id]
        for residue in chain.get_list():
            if residue.resname in sabr_constants.AA_3TO1:
                label = residue_label(chain_id, residue)
                if label not in chosen:
                    fixed.append(label)
        if not chosen:
            raise ValueError(f"SAbR found no selected CDR residues in chain {chain_id!r}.")
    return list(dict.fromkeys(fixed)), selected_by_region


def fixed_residues_from_annotation(
    path: Path,
    chain_ids: list[str],
    cdr_regions: list[str],
    residues_by_region: dict[str, list[str]],
) -> tuple[list[str], dict[str, list[str]]]:
    """Use upstream SAbR CDR mapping to fix every non-selected antibody residue."""
    unknown = set(cdr_regions) - set(CDR_NAMES)
    if unknown:
        raise ValueError(f"Unknown CDR region(s): {', '.join(sorted(unknown))}.")
    if not cdr_regions:
        raise ValueError("Choose at least one CDR to redesign.")
    structure = read_structure(Path(path))
    model = next(structure.get_models())
    selected = {name: [] for name in CDR_NAMES}
    selected_labels: set[str] = set()
    for region in cdr_regions:
        for label in residues_by_region.get(region, []):
            chain_id = next((chain for chain in chain_ids if str(label).startswith(chain)), None)
            if chain_id:
                selected[region].append(str(label))
                selected_labels.add(str(label))
    fixed: list[str] = []
    for chain_id in chain_ids:
        if chain_id not in model:
            raise ValueError(f"Binder chain {chain_id!r} is not present in {Path(path).name}.")
        for residue in model[chain_id].get_list():
            if residue.resname in sabr_constants.AA_3TO1:
                label = residue_label(chain_id, residue)
                if label not in selected_labels:
                    fixed.append(label)
    if not selected_labels:
        raise ValueError("The SAbR annotation contains no residues in the selected CDRs.")
    return list(dict.fromkeys(fixed)), selected


def write_renumbered(path: Path, destination: Path, chain_ids: list[str], *, scheme: str,
                     chain_type: str, noise_level: float, mode: str) -> None:
    structure = read_structure(Path(path))
    for chain_id in chain_ids:
        structure = renumber_structure(
            structure, chain_id, scheme=scheme, chain_type=chain_type,
            noise_level=noise_level, mode=mode,
        )
    writer = MMCIFIO() if Path(destination).suffix.lower() in {".cif", ".mmcif"} else PDBIO()
    writer.set_structure(structure)
    writer.save(str(destination))
