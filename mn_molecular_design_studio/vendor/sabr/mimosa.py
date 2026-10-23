"""Prepare antibody CDR motifs for MIMOSA-guided design runs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.Align import PairwiseAligner
from scipy.spatial.distance import cdist

from . import constants
from .api import cdr_residue_records
from .structure import extract_chain

AA3_TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def _load_structure(path: str | Path):
    path = Path(path)
    parser = (
        MMCIFParser(QUIET=True)
        if path.suffix.lower() in {".cif", ".mmcif"}
        else PDBParser(QUIET=True)
    )
    return parser.get_structure(path.stem, str(path))


def _chain_atoms(chain) -> tuple[list, np.ndarray]:
    atoms = [
        atom
        for residue in chain
        for atom in residue.get_atoms()
        if residue.id[0] == " "
        and residue.resname.upper() in AA3_TO1
        and atom.element.upper() != "H"
    ]
    if not atoms:
        raise ValueError(f"Chain {chain.id!r} has no heavy atoms.")
    return atoms, np.asarray([atom.coord for atom in atoms], dtype=np.float64)


def _sequence_alignment_pairs(
    source_sequence: str, target_sequence: str
) -> list[tuple[int, int]]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(source_sequence, target_sequence)[0]
    pairs = []
    source_blocks, target_blocks = alignment.aligned
    for (source_start, source_end), (target_start, target_end) in zip(
        source_blocks, target_blocks
    ):
        for offset in range(min(source_end - source_start, target_end - target_start)):
            pairs.append((int(source_start + offset), int(target_start + offset)))
    return pairs


def prepare_mimosa_motif(
    complex_path: str | Path,
    target_path: str | Path,
    *,
    binder_chain: str,
    complex_target_chain: str,
    target_chain: str,
    cdr: str | None = None,
    contact_cutoff: float = 8.0,
    chain_type: str = "auto",
    selected_imgt_positions: list[str] | None = None,
    selected_source_residues: list[str] | None = None,
) -> dict[str, Any]:
    """Extract an interface CDR motif and align it to the design target.

    The output is JSON-serializable. Each motif residue carries its source
    chain and author residue ID, source offset, amino-acid identity, and
    C-alpha coordinate transformed into the supplied design target's frame.
    Source partner chains may be comma-separated when a motif spans chains.

    ``contact_cutoff`` is an all-heavy-atom minimum-distance cutoff in
    angstroms. Target-chain residues are aligned by sequence and their C-alpha
    coordinates are superposed before transferring binder coordinates.
    """
    normalized_cdr = cdr.strip().upper() if cdr else None
    if normalized_cdr is not None and normalized_cdr not in constants.IMGT_LOOPS:
        raise ValueError(f"cdr must be one of {', '.join(constants.IMGT_LOOPS)}.")
    if not np.isfinite(contact_cutoff) or contact_cutoff <= 0:
        raise ValueError("contact_cutoff must be a positive finite number.")

    complex_structure = _load_structure(complex_path)
    target_structure = _load_structure(target_path)
    source_model = next(iter(complex_structure))
    target_model = next(iter(target_structure))
    binder_chain_ids = [value.strip() for value in str(binder_chain).split(",") if value.strip()]
    if not binder_chain_ids or len(binder_chain_ids) != len(set(binder_chain_ids)):
        raise ValueError("binder_chain must contain one or more unique comma-separated chain IDs.")
    try:
        source_binder_chains = [source_model[chain_id] for chain_id in binder_chain_ids]
        source_target_chain = source_model[complex_target_chain]
        destination_target_chain = target_model[target_chain]
    except KeyError as error:
        raise ValueError(f"Requested chain is absent: {error}.") from error

    if complex_target_chain in binder_chain_ids:
        raise ValueError("The source partner and source target must use different chains.")
    chain_residues = []
    protein_residues = []
    for chain_id, chain in zip(binder_chain_ids, source_binder_chains):
        for residue in chain:
            global_index = len(chain_residues)
            chain_residues.append(residue)
            if residue.id[0] == " " and residue.resname.upper() in AA3_TO1:
                protein_residues.append((global_index, chain_id, residue))
    cdr_records = []
    if selected_source_residues is None and normalized_cdr is not None:
        if len(binder_chain_ids) != 1:
            raise ValueError("Source CDR extraction requires exactly one antibody partner chain.")
        cdr_records = cdr_residue_records(
            complex_structure,
            binder_chain_ids[0],
            chain_type=chain_type,
        )[normalized_cdr]
        if not cdr_records:
            raise ValueError(
                f"SAbR found no residues in {normalized_cdr} of chain {binder_chain_ids[0]!r}."
            )

    source_target = extract_chain(complex_structure, complex_target_chain, None)
    destination_target = extract_chain(target_structure, target_chain, None)
    residue_pairs = _sequence_alignment_pairs(
        source_target.sequence, destination_target.sequence
    )
    source_ca_atoms = []
    destination_ca_atoms = []
    source_target_chain_residues = list(source_target_chain)
    destination_target_chain_residues = list(destination_target_chain)
    for source_idx, destination_idx in residue_pairs:
        source_chain_idx = source_target.residue_indices[source_idx]
        destination_chain_idx = destination_target.residue_indices[destination_idx]
        source_residue = source_target_chain_residues[source_chain_idx]
        destination_residue = destination_target_chain_residues[destination_chain_idx]
        if "CA" not in source_residue or "CA" not in destination_residue:
            continue
        if source_target.sequence[source_idx] != destination_target.sequence[destination_idx]:
            continue
        source_ca_atoms.append(source_residue["CA"])
        destination_ca_atoms.append(destination_residue["CA"])
    if len(source_ca_atoms) < 3:
            raise ValueError(
                "Target chains need at least three sequence-matched C-alpha atoms to align."
            )

    superimposer = Superimposer()
    superimposer.set_atoms(destination_ca_atoms, source_ca_atoms)
    rotation, translation = superimposer.rotran
    target_alignment = []
    for source_idx, destination_idx in residue_pairs:
        if source_target.sequence[source_idx] != destination_target.sequence[destination_idx]:
            continue
        destination_chain_idx = destination_target.residue_indices[destination_idx]
        destination_residue = destination_target_chain_residues[destination_chain_idx]
        if "CA" not in destination_residue:
            continue
        target_alignment.append({
            "target_offset": int(destination_idx),
            "target_chain_id": target_chain,
            "target_residue_id": f"{destination_residue.id[1]}{str(destination_residue.id[2] or '').strip()}",
            "insertion_code": str(destination_residue.id[2] or "").strip(),
            "amino_acid": AA3_TO1.get(destination_residue.resname.upper(), "X"),
            "ca": np.asarray(destination_residue["CA"].coord, dtype=np.float64).tolist(),
        })
    if len(target_alignment) < 3:
        raise ValueError("MIMOSA target fit needs at least three sequence-matched target C-alpha atoms.")

    _, target_tree_coords = _chain_atoms(source_target_chain)
    if cdr_records:
        source_records = [
            (
                int(record[0]), binder_chain_ids[0], int(record[1]),
                str(record[2] or ""), str(record[3]), offset,
            )
            for offset, record in enumerate(cdr_records)
        ]
    else:
        residue_offsets = {
            chain_index: offset
            for offset, (chain_index, _chain_id, _residue) in enumerate(protein_residues)
        }
        source_records = [
            (
                chain_index,
                chain_id,
                int(residue.id[1]),
                str(residue.id[2] or "").strip(),
                AA3_TO1[residue.resname.upper()],
                residue_offsets[chain_index],
            )
            for chain_index, chain_id, residue in protein_residues
        ]
    source_chain_length = len(source_records)
    source_residue_ids = []
    for _index, chain_id, number, insertion_code, _aa, _offset in source_records:
        residue_id = f"{number}{insertion_code}".strip()
        source_residue_ids.append(f"{chain_id}:{residue_id}" if len(binder_chain_ids) > 1 else residue_id)
    all_source_records = list(source_records)
    requested_source: list[str] | None = None
    if selected_source_residues is not None:
        requested_source = [str(value).strip().upper() for value in selected_source_residues]
        if not requested_source:
            raise ValueError("Select at least one contacting amino acid from the source partner chain.")
        if len(set(requested_source)) != len(requested_source):
            raise ValueError("Selected source partner residue positions must be unique.")
        by_selector = {}
        for record in all_source_records:
            _index, chain_id, number, insertion, _aa, _offset = record
            residue_id = f"{number}{insertion}".upper()
            by_selector[f"{chain_id}:{residue_id}".upper()] = record
            if len(binder_chain_ids) == 1:
                by_selector[residue_id] = record
        resolved_source = []
        missing = []
        for value in requested_source:
            selector = value if ":" in value or len(binder_chain_ids) == 1 else ""
            record = by_selector.get(selector)
            if record is None:
                missing.append(value)
            else:
                resolved_source.append(record)
        if missing:
            raise ValueError("Selected source partner residues are absent from the selected chain(s): " + ", ".join(missing))
        requested_source = [
            f"{record[1]}:{record[2]}{record[3]}" if len(binder_chain_ids) > 1 else f"{record[2]}{record[3]}"
            for record in resolved_source
        ]

    motif = []
    for chain_index, source_chain_id, residue_number, insertion_code, amino_acid, source_offset in all_source_records:
        try:
            residue = chain_residues[chain_index]
        except IndexError as error:
            raise ValueError("SAbR returned a source-chain index outside the parsed chain.") from error
        if "CA" not in residue:
            continue
        binder_heavy = [
            atom.coord for atom in residue.get_atoms() if atom.element.upper() != "H"
        ]
        if not binder_heavy:
            continue
        distances = cdist(np.asarray(binder_heavy), target_tree_coords)
        if float(distances.min()) > contact_cutoff:
            continue
        transformed_ca = (
            np.asarray(residue["CA"].coord, dtype=np.float64) @ rotation
            + translation
        )
        motif.append(
            {
                "source_chain_index": int(chain_index),
                "source_chain_id": source_chain_id,
                "source_residue_id": f"{residue_number}{insertion_code}",
                "source_selector": (
                    f"{source_chain_id}:{residue_number}{insertion_code}"
                    if len(binder_chain_ids) > 1 else f"{residue_number}{insertion_code}"
                ),
                "imgt_position": int(residue_number) if cdr_records else None,
                "insertion_code": insertion_code,
                "cdr_offset": int(source_offset),
                "source_offset": int(source_offset),
                "amino_acid": amino_acid,
                "ca": transformed_ca.tolist(),
            }
        )
    contact_motif = list(motif)
    if requested_source is not None:
        by_selector = {str(row["source_selector"]).upper(): row for row in contact_motif}
        if len(binder_chain_ids) == 1:
            by_selector.update({str(row["source_residue_id"]).upper(): row for row in contact_motif})
        missing = [value for value in requested_source if value not in by_selector]
        if missing:
            raise ValueError(
                "Selected source partner residues do not contact the source target within "
                f"{contact_cutoff:g} angstroms: " + ", ".join(missing)
            )
        motif = sorted(
            [by_selector[value] for value in requested_source],
            key=lambda row: int(row["source_offset"]),
        )
    if selected_imgt_positions:
        def position_key(number: int, insertion_code: str) -> str:
            return f"{int(number)}{str(insertion_code or '').strip().upper()}"

        requested = [str(value).strip().upper() for value in selected_imgt_positions]
        if len(set(requested)) != len(requested):
            raise ValueError("Selected IMGT motif positions must be unique.")
        available = {position_key(row["imgt_position"], row["insertion_code"]): row for row in contact_motif}
        missing = [value for value in requested if value not in available]
        if missing:
            raise ValueError(
                "Selected IMGT positions are not contacting residues in the chosen source CDR: "
                + ", ".join(missing)
            )
        motif = [available[value] for value in requested]
    if not motif:
        raise ValueError(
            f"No selected source-partner residues contact the target within {contact_cutoff:g} angstroms."
        )

    return {
        "schema": "sabr-mimosa-motif-v2",
        "algorithm_revision": "mimosa-si-a1-a5-v1",
        "guidance_progress_convention": "linear_reverse_step_index_0_to_1",
        "selection": {
            "method": "selected_source_residues" if selected_source_residues is not None else ("selected_imgt_positions" if selected_imgt_positions else "all_contacting_source_cdr_residues" if cdr_records else "all_contacting_source_partner_residues"),
            "requested_source_residues": list(selected_source_residues or []),
            "available_contact_source_residues": [row["source_selector"] for row in contact_motif],
            "requested_imgt_positions": list(selected_imgt_positions or []),
            "available_contact_imgt_positions": [
                f"{row['imgt_position']}{row['insertion_code']}".strip()
                for row in contact_motif
            ],
        },
        "source": {
            "complex_path": str(Path(complex_path)),
            "target_path": str(Path(target_path)),
            "binder_chain": binder_chain,
            "binder_chains": binder_chain_ids,
            "complex_target_chain": complex_target_chain,
            "target_chain": target_chain,
            "target_sequence": destination_target.sequence,
            "cdr": normalized_cdr,
            "source_chain_length": source_chain_length,
            "source_residue_ids": source_residue_ids,
            "cdr_sequence": "".join(record[3] for record in cdr_records),
            "cdr_imgt_positions": [int(record[1]) for record in cdr_records],
            "cdr_residues": [
                {
                    "cdr_offset": int(offset),
                    "imgt_position": int(record[1]),
                    "insertion_code": record[2],
                    "amino_acid": record[3],
                }
                for offset, record in enumerate(cdr_records)
            ],
            "target_sequence_identity_fraction": float(
                np.mean(
                    [
                        source_target.sequence[i]
                        == destination_target.sequence[j]
                        for i, j in residue_pairs
                    ]
                )
            ),
            "target_alignment_rmsd_angstrom": float(superimposer.rms),
            "contact_cutoff_angstrom": float(contact_cutoff),
        },
        "target_alignment": target_alignment,
        "motif": motif,
    }


def write_mimosa_motif(packet: dict[str, Any], output_path: str | Path) -> Path:
    """Write a prepared MIMOSA motif packet as readable JSON."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(packet, indent=2) + "\n", encoding="utf-8")
    return path
