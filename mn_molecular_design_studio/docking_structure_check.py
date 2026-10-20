"""Quick receptor integrity checks shared by structure-based ligand tools."""

from __future__ import annotations

import math
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.auth import user_input_root


_PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASX", "GLX", "MSE", "SEC", "PYL", "HID", "HIE", "HIP", "CYX", "CYM",
    "HYP", "SEP", "TPO", "PTR",
}


def inspect_receptor_pdb(text: str) -> dict:
    """Report obvious coordinate, backbone, and chain-continuity problems.

    This is a fast input-page screen, not a full chemical validation or a
    replacement for the Complex Preparation workflow.
    """
    malformed_coordinates = 0
    model_count = sum(line.startswith("MODEL ") for line in text.splitlines())
    residues: dict[tuple[str, int, str, str], dict] = {}
    chain_residues: dict[str, list[tuple[tuple[str, int, str, str], dict]]] = {}

    for line in text.splitlines():
        record = line[:6].strip()
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 54:
            if record == "ATOM":
                malformed_coordinates += 1
            continue
        try:
            coordinates = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
            if not all(math.isfinite(value) for value in coordinates):
                raise ValueError("non-finite coordinate")
        except ValueError:
            if record == "ATOM":
                malformed_coordinates += 1
            continue
        residue_name = line[17:20].strip().upper()
        if residue_name not in _PROTEIN_RESIDUES:
            continue
        try:
            number = int(line[22:26])
        except ValueError:
            malformed_coordinates += 1
            continue
        chain = line[21].strip() or "_"
        insertion = line[26].strip()
        key = (chain, number, insertion, residue_name)
        altloc = line[16].strip()
        # Alternate location B and later are alternatives, not additional
        # atoms for this quick completeness check.
        if altloc not in {"", "A"}:
            continue
        atom_name = line[12:16].strip().upper()
        residue = residues.get(key)
        if residue is None:
            residue = {"atoms": set(), "coordinates": {}}
            residues[key] = residue
            chain_residues.setdefault(chain, []).append((key, residue))
        residue["atoms"].add(atom_name)
        if atom_name in {"N", "CA", "C", "O"}:
            residue["coordinates"].setdefault(atom_name, coordinates)

    findings: list[str] = []
    if malformed_coordinates:
        findings.append(
            f"{malformed_coordinates} protein atom record(s) have incomplete or invalid coordinates."
        )
    if not residues:
        findings.append("No recognizable protein residues were found in the PDB ATOM records.")

    missing_backbone = []
    broken_bonds = []
    numbering_gaps = []
    for chain, rows in sorted(chain_residues.items()):
        for key, residue in rows:
            missing = {"N", "CA", "C", "O"} - residue["atoms"]
            if missing:
                missing_backbone.append(
                    f"Chain {chain}, {key[3]} {key[1]}{key[2]} lacks {', '.join(sorted(missing))}"
                )
        for (left_key, left), (right_key, right) in zip(rows, rows[1:]):
            left_number, right_number = left_key[1], right_key[1]
            if right_number > left_number + 1:
                numbering_gaps.append(
                    f"Chain {chain} jumps from residue {left_number} to {right_number}"
                )
            carbon = left["coordinates"].get("C")
            nitrogen = right["coordinates"].get("N")
            if carbon is not None and nitrogen is not None:
                distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(carbon, nitrogen)))
                if distance > 2.2:
                    broken_bonds.append(
                        f"Chain {chain} has a possible backbone break between residues "
                        f"{left_number} and {right_number} ({distance:.1f} Å C–N distance)"
                    )

    if missing_backbone:
        findings.extend(missing_backbone[:12])
        if len(missing_backbone) > 12:
            findings.append(f"…and {len(missing_backbone) - 12} more residues with missing backbone atoms.")
    if broken_bonds:
        findings.extend(broken_bonds[:8])
        if len(broken_bonds) > 8:
            findings.append(f"…and {len(broken_bonds) - 8} more possible backbone breaks.")
    if numbering_gaps:
        findings.extend(numbering_gaps[:8])
        if len(numbering_gaps) > 8:
            findings.append(f"…and {len(numbering_gaps) - 8} more residue-number gaps.")
    if model_count > 1:
        findings.append(f"The file contains {model_count} models; review which model should be docked.")

    return {
        "chains": len(chain_residues),
        "residues": len(residues),
        "findings": findings,
    }


def render_receptor_check(text: str, source_path: Path, *, owner_id: str, key: str) -> dict:
    """Show the quick receptor check and a direct path to the repair tile."""
    report = inspect_receptor_pdb(text)
    if report["findings"]:
        st.warning(
            f"Structure check found {len(report['findings'])} possible issue(s) "
            f"in {report['chains']} chain(s) / {report['residues']} protein residues. "
            "Review the findings before docking."
        )
        with st.expander("Structure check findings", expanded=True):
            for finding in report["findings"]:
                st.markdown(f"- {finding}")
            st.caption(
                "This quick check covers coordinate records, backbone atoms, chain gaps, "
                "and peptide-bond continuity. It does not validate protonation or all side-chain chemistry."
            )
            st.button(
                "Repair in Complex Preparation",
                key=key,
                help="Open the repair tile with this PDB preloaded. Repair is a separate job; the source remains unchanged.",
                on_click=open_complex_preparation,
                args=(str(source_path),),
                kwargs={"owner_id": owner_id},
                type="primary",
            )
    else:
        st.caption(
            f"Structure check: {report['chains']} chain(s), {report['residues']} protein residues; "
            "no obvious coordinate or backbone issues found."
        )
    return report


def open_complex_preparation(source_path: str, *, owner_id: str) -> None:
    """Route to Complex Preparation with the current private PDB preselected."""
    source = Path(source_path).expanduser().resolve(strict=True)
    source.relative_to(user_input_root(owner_id))
    if source.suffix.lower() not in {".pdb", ".ent"}:
        raise ValueError("Complex Preparation requires a saved PDB structure.")
    st.session_state["ligand-complex-prepare-pending-source"] = str(source)
    st.session_state.pop("selected_workflow_group", None)
    st.session_state["selected_workflow"] = "ligand.complex.prepare"
    st.session_state["page"] = "Workflow setup"
