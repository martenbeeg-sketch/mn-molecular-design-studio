"""Local structure helpers for choosing source-interface motif residues."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1
from scipy.spatial.distance import cdist

import streamlit as st

from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
    StructureVisualization,
    molstar_custom_component,
)


PAPER_INTERACTION_TEMPLATES = {
    "3PB1": "uPA · PAI-1 SARMAP loop",
    "2FLU": "KEAP1 · Nrf2 ETGE motif",
    "6WZM": "IL-8 · LY3041658 antibody interaction",
    "8IC0": "IL-8 · CXCR1 receptor interaction",
}

PAPER_MIMOSA_EXAMPLES = {
    "keap1_nrf2": {
        "label": "Fig. 2c · KEAP1 / Nrf2 ETGE motif (2FLU)",
        "target_name": "KEAP1 / Nrf2 MIMOSA paper example",
        "target_id": "2FLU",
        "source_template": "2FLU",
        "target_chain": "X",
        "source_partner_chain": "P",
        "source_target_chain": "X",
        "motif_positions": ["79", "80", "81", "82"],
        "condition": "sequence_fixed",
        "design_loop": "CDR3",
        "rfantibody_regions": ["H1", "H2", "H3"],
    },
    "upa_pai1": {
        "label": "Fig. 2b · uPA / PAI-1 SARMAP loop (3PB1)",
        "target_name": "uPA / PAI-1 MIMOSA paper example",
        "target_id": "3PB1",
        "source_template": "3PB1",
        "target_chain": "E",
        "source_partner_chain": "I",
        "source_target_chain": "E",
        "motif_positions": ["344", "345", "346", "347", "348", "349"],
        "condition": "sequence_fixed",
        "design_loop": "CDR3",
    },
    "il8_cxcr1": {
        "label": "Fig. 3b · IL-8 / CXCR1 receptor motif (8IC0)",
        "target_name": "IL-8 / CXCR1 MIMOSA paper example",
        "target_id": "8IC0",
        "source_template": "8IC0",
        "target_chain": "F",
        "source_partner_chain": "A",
        "source_target_chain": "F",
        "motif_positions": ["21", "22", "23", "26", "27", "28", "29"],
        "condition": "sequence_fixed",
        "design_loop": "all_cdrs",
    },
    "il8_ly3041658": {
        "label": "Fig. 3c · IL-8 / LY3041658 antibody motif (6WZM)",
        "target_name": "IL-8 / LY3041658 MIMOSA paper example",
        "target_id": "6WZM",
        "source_template": "6WZM",
        "target_chain": "E",
        "source_partner_chain": "A,B",
        "source_target_chain": "E",
        "motif_positions": ["A:33", "A:99", "A:102", "A:104", "A:105", "A:106", "A:110", "B:94"],
        "condition": "sequence_fixed",
        "design_loop": "all_cdrs",
    },
}


def _fetch_rcsb_structure(accession: str, extension: str, destination: str | Path) -> Path:
    destination = Path(destination)
    if destination.is_file() and destination.stat().st_size:
        return destination
    request = Request(
        f"https://files.rcsb.org/download/{accession}.{extension}",
        headers={"User-Agent": "MN-Molecular-Design-Studio/0.1"},
    )
    try:
        with urlopen(request, timeout=30) as response:
            content = response.read(50_000_001)
    except (HTTPError, URLError, TimeoutError, OSError) as error:
        raise ValueError(f"Could not fetch RCSB structure {accession}.{extension}: {error}") from error
    if len(content) > 50_000_000:
        raise ValueError(f"RCSB structure {accession}.{extension} exceeds the 50 MB structure limit.")
    if not content.strip():
        raise ValueError(f"RCSB returned an empty structure for {accession}.{extension}.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(content)
    return destination


def fetch_paper_interaction_template(template_id: str, destination: str | Path) -> Path:
    """Download one source interaction complex cited in the MIMOSA paper."""
    normalized = str(template_id).strip().upper()
    if normalized not in PAPER_INTERACTION_TEMPLATES:
        raise ValueError(f"Unsupported MIMOSA interaction template {template_id!r}.")
    return _fetch_rcsb_structure(normalized, "cif", destination)


def fetch_paper_mimosa_example(example_id: str, destination_dir: str | Path) -> dict[str, object]:
    """Load both target and interaction structures for a paper MIMOSA example."""
    normalized = str(example_id).strip().lower()
    try:
        example = PAPER_MIMOSA_EXAMPLES[normalized]
    except KeyError as error:
        raise ValueError(f"Unsupported MIMOSA paper example {example_id!r}.") from error
    destination_dir = Path(destination_dir)
    target_path = _fetch_rcsb_structure(
        str(example["target_id"]), "pdb", destination_dir / f"{example['target_id']}-target.pdb"
    )
    source_path = fetch_paper_interaction_template(
        str(example["source_template"]), destination_dir / f"{example['source_template']}-interaction.cif"
    )
    return {**example, "target_path": target_path, "source_path": source_path}


def _load_model(path: str | Path):
    structure_path = Path(path)
    parser = (
        MMCIFParser(QUIET=True)
        if structure_path.suffix.casefold() in {".cif", ".mmcif"}
        else PDBParser(QUIET=True)
    )
    structure = parser.get_structure(structure_path.stem, str(structure_path))
    return next(iter(structure))


def _protein_residue(residue) -> str | None:
    if residue.id[0] != " ":
        return None
    amino_acid = seq1(residue.resname, custom_map={"MSE": "M"}).upper()
    return amino_acid if len(amino_acid) == 1 and amino_acid in "ACDEFGHIKLMNPQRSTVWY" else None


def source_complex_chains(path: str | Path) -> list[str]:
    """Return chain IDs with standard amino acids in a source complex."""
    model = _load_model(path)
    result = []
    for chain in model:
        if any(_protein_residue(residue) for residue in chain):
            result.append(str(chain.id))
    return result


def contacting_partner_residues(
    path: str | Path,
    partner_chain_id: str | list[str],
    target_chain_id: str | list[str],
    contact_cutoff: float = 8.0,
) -> list[dict[str, str]]:
    """List source partner residues contacting the source target chain.

    The returned ``position`` is the source structure's author residue number
    and insertion code, suitable for passing to the MIMOSA packet preparer.
    """
    model = _load_model(path)
    partner_chain_ids = _normalize_chain_ids(partner_chain_id)
    target_chain_ids = _normalize_chain_ids(target_chain_id)
    if set(partner_chain_ids) & set(target_chain_ids):
        raise ValueError("Choose different source-complex partner and target chains.")
    try:
        partner_chains = [model[chain_id] for chain_id in partner_chain_ids]
        target_chains = [model[chain_id] for chain_id in target_chain_ids]
    except KeyError as error:
        raise ValueError(f"Source complex does not contain chain {error.args[0]!r}.") from error
    target_atoms = [
        atom
        for target_chain in target_chains
        for residue in target_chain
        for atom in residue.get_atoms()
        if atom.element.upper() != "H" and _protein_residue(residue)
    ]
    if not target_atoms:
        raise ValueError(f"Source target chain(s) {', '.join(target_chain_ids)!r} have no standard protein atoms.")
    target_coords = np.asarray([atom.coord for atom in target_atoms], dtype=np.float64)
    results = []
    qualified = len(partner_chain_ids) > 1
    for chain_id, partner_chain in zip(partner_chain_ids, partner_chains):
        for residue in partner_chain:
            amino_acid = _protein_residue(residue)
            if amino_acid is None or "CA" not in residue:
                continue
            atoms = [atom for atom in residue.get_atoms() if atom.element.upper() != "H"]
            if not atoms:
                continue
            coordinates = np.asarray([atom.coord for atom in atoms], dtype=np.float64)
            minimum_distance = float(cdist(coordinates, target_coords).min())
            if minimum_distance > float(contact_cutoff):
                continue
            number = int(residue.id[1])
            insertion = str(residue.id[2] or "").strip().upper()
            residue_id = f"{number}{insertion}"
            position = f"{chain_id}:{residue_id}" if qualified else residue_id
            results.append(
                {
                    "position": position,
                    "label": f"{chain_id}:{amino_acid}{residue_id}" if qualified else f"{amino_acid}{residue_id}",
                    "amino_acid": amino_acid,
                    "chain_id": chain_id,
                    "distance": f"{minimum_distance:.1f} Å",
                }
            )
    return results


def _normalize_chain_ids(value: str | list[str]) -> list[str]:
    if isinstance(value, str):
        chain_ids = [item.strip() for item in value.split(",") if item.strip()]
    else:
        chain_ids = [str(item).strip() for item in value if str(item).strip()]
    if not chain_ids or len(chain_ids) != len(set(chain_ids)):
        raise ValueError("Choose one or more unique source-complex chain IDs.")
    return chain_ids


def render_source_motif_viewer(
    path: str | Path,
    partner_chain_id: str | list[str],
    selected_positions: list[str],
    *,
    key_prefix: str,
    height: int = 440,
) -> None:
    """Show the source interaction complex and highlight selected motif residues."""
    structure_path = Path(path)
    if not structure_path.is_file():
        return

    partner_chain_ids = _normalize_chain_ids(partner_chain_id)
    selections = []
    unhighlightable = False
    for position in selected_positions:
        token = str(position).strip()
        qualified = re.fullmatch(r"([A-Za-z0-9]+):(-?\d+[A-Za-z]?)", token)
        if qualified:
            selections.append(f"{qualified.group(1)}:{qualified.group(2)}")
        elif re.fullmatch(r"-?\d+[A-Za-z]?", token) and len(partner_chain_ids) == 1:
            selections.append(f"{partner_chain_ids[0]}:{token}")
        else:
            unhighlightable = True

    structure_text = structure_path.read_text(errors="replace")
    digest_input = "\0".join(
        (str(structure_path.resolve()), ",".join(partner_chain_ids), ",".join(selected_positions))
    )
    digest = hashlib.sha1(digest_input.encode()).hexdigest()[:12]
    st.markdown("#### Source motif in MolView")
    st.caption(
        f"The source complex is shown below. Selected amino acids on partner chain(s) {', '.join(partner_chain_ids)} "
        "are highlighted; these are the motif residues MIMOSA transfers to the generated loop."
    )
    if unhighlightable:
        st.caption("Residues with insertion codes remain selected in the list but cannot be highlighted in MolView.")
    molstar_custom_component(
        [
            StructureVisualization(
                pdb=structure_text,
                color="uniform",
                color_params={"value": "0xe7b8aa"},
                representation_type="cartoon",
                highlighted_selections=selections,
            )
        ],
        key=f"{key_prefix}-source-motif-viewer-{digest}",
        height=height,
        show_controls=True,
        show_structure_tools=False,
        selection_mode=False,
    )
