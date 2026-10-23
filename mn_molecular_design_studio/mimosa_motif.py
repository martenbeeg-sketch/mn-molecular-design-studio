"""Recover and carry MIMOSA sequence-fixed motif identities between workflows."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any


_BOLTZGEN_CANONICAL_AA = "ARNDCQEGHILKMFPSTWYV"
_BOLTZGEN_TOKEN_TO_AA = {index + 2: aa for index, aa in enumerate(_BOLTZGEN_CANONICAL_AA)}


def _structure_residues(structure_path: Path) -> dict[str, list[dict[str, Any]]]:
    from Bio.SeqUtils import seq1

    from mn_molecular_design_studio.design_validation import read_structure

    model = next(read_structure(structure_path).get_models(), None)
    if model is None:
        raise ValueError("The MIMOSA candidate has no coordinate model.")
    chains: dict[str, list[dict[str, Any]]] = {}
    for chain in model:
        rows = []
        for residue in chain:
            if "CA" not in residue:
                continue
            amino_acid = seq1(str(residue.resname).strip(), undef_code="X").upper()
            if amino_acid == "X":
                continue
            number = int(residue.id[1])
            insertion = str(residue.id[2] or "").strip()
            label = f"{chain.id}{number}{insertion}"
            rows.append({"label": label, "number": number, "insertion_code": insertion, "amino_acid": amino_acid})
        if rows:
            chains[str(chain.id)] = rows
    return chains


def _boltzgen_fixed_identities(
    metadata_path: Path,
    structure_path: Path,
    binder_chains: list[str],
    target_chains: list[str],
) -> tuple[dict[str, str], int]:
    import numpy as np

    with np.load(metadata_path, allow_pickle=False) as archive:
        if "mimosa_condition" not in archive.files or "mimosa_fixed_aa" not in archive.files:
            return {}, 0
        condition = str(archive["mimosa_condition"].item()).strip().lower()
        if condition != "sequence_fixed":
            return {}, 0
        fixed = np.asarray(archive["mimosa_fixed_aa"], dtype=np.int64)
        if fixed.ndim != 1:
            raise ValueError("BoltzGen MIMOSA fixed-residue metadata is not a residue vector.")
        positions = [int(index) for index in np.flatnonzero(fixed >= 0)]
        if not positions:
            return {}, 0
        if "source_context" not in archive.files:
            raise ValueError("BoltzGen MIMOSA metadata has no chain-to-residue mapping.")
        context = json.loads(str(archive["source_context"].item()))

    chains = _structure_residues(structure_path)
    binder_ids = [str(chain) for chain in binder_chains]
    target_ids = {str(chain) for chain in target_chains}
    identities: dict[str, str] = {}
    offset = 0
    token_count = 0
    for context_chain in context.get("chains") or []:
        indices = context_chain.get("indices") or []
        chain_tokens = len(indices)
        start, stop = offset, offset + chain_tokens
        selected_positions = [position for position in positions if start <= position < stop]
        offset = stop
        token_count = stop
        if not selected_positions:
            continue

        source_chain = str(context_chain.get("source_chain") or "")
        output_chain = str(context_chain.get("output_chain") or "")
        direct_matches = [chain for chain in (source_chain, output_chain) if chain in binder_ids]
        if len(set(direct_matches)) == 1:
            structure_chain = direct_matches[0]
        elif len(set(direct_matches)) > 1:
            raise ValueError("BoltzGen MIMOSA chain mapping is ambiguous.")
        elif len(binder_ids) == 1 and source_chain not in target_ids and output_chain not in target_ids:
            structure_chain = binder_ids[0]
        else:
            raise ValueError("Could not map a MIMOSA fixed motif chain to the selected antibody chain.")
        residue_rows = chains.get(structure_chain) or []
        if not residue_rows:
            raise ValueError(f"The antibody chain {structure_chain} has no standard amino-acid residues.")
        if len(residue_rows) == chain_tokens:
            token_to_residue = {start + index: row for index, row in enumerate(residue_rows)}
        else:
            # Cropped BoltzGen contexts store source residue indices. Use these
            # only when they identify a unique author residue in the output.
            token_to_residue = {}
            by_number: dict[int, list[dict[str, Any]]] = {}
            for row in residue_rows:
                by_number.setdefault(int(row["number"]), []).append(row)
            for local_index, source_index in enumerate(indices):
                matches = by_number.get(int(source_index), [])
                if len(matches) == 1:
                    token_to_residue[start + local_index] = matches[0]

        for position in selected_positions:
            residue = token_to_residue.get(position)
            amino_acid = _BOLTZGEN_TOKEN_TO_AA.get(int(fixed[position]))
            if residue is None or amino_acid is None:
                raise ValueError(f"Could not map fixed BoltzGen MIMOSA token {position} to a canonical antibody residue.")
            label = str(residue["label"])
            if label in identities:
                raise ValueError(f"BoltzGen MIMOSA mapped two fixed tokens to residue {label}.")
            identities[label] = amino_acid

    if token_count != len(fixed) or len(identities) != len(positions):
        raise ValueError("BoltzGen MIMOSA fixed residues do not all map to the generated antibody structure.")
    return identities, len(positions)


def _rfantibody_pose_residue_map(metadata_path: Path, structure_path: Path) -> dict[int, tuple[str, dict[str, Any]]]:
    """Map RFdiffusion's global pose indices to saved chain-local residues.

    RFantibody's sampler sidecar uses the input framework's author numbering,
    while the QV candidate writer renumbers each saved chain from one. The QV
    CA-atom order is the stable bridge between those two representations.
    """
    stages_dir = metadata_path.parent
    qv_path = stages_dir / "1_rfdiffusion.qv"
    if not qv_path.is_file():
        return {}
    name = metadata_path.name
    prefix = "1_rfdiffusion_"
    suffix = ".mimosa.json"
    if not (name.startswith(prefix) and name.endswith(suffix)):
        return {}
    tag = name[len(prefix):-len(suffix)]
    qv_records: dict[str, list[str]] = {}
    current_tag: str | None = None
    for line in qv_path.read_text(errors="replace").splitlines():
        if line.startswith("QV_TAG "):
            current_tag = line.split(maxsplit=1)[1].strip()
            qv_records[current_tag] = []
        elif current_tag == tag and line.startswith("ATOM") and line[12:16].strip() == "CA":
            qv_records[current_tag].append(line)
    lines = qv_records.get(tag) or []
    if not lines:
        return {}

    chains = _structure_residues(structure_path)
    chain_indices: dict[str, int] = {}
    mapped: dict[int, tuple[str, dict[str, Any]]] = {}
    for pose_index, line in enumerate(lines):
        chain = line[21].strip()
        local_index = chain_indices.get(chain, 0)
        chain_indices[chain] = local_index + 1
        residues = chains.get(chain) or []
        if local_index >= len(residues):
            raise ValueError(f"RFantibody QV pose index {pose_index} exceeds saved chain {chain} length.")
        mapped[pose_index] = (chain, residues[local_index])
    for chain, count in chain_indices.items():
        if count != len(chains.get(chain) or []):
            raise ValueError(f"RFantibody QV and saved structure disagree on chain {chain} residue count.")
    return mapped


def _rfantibody_fixed_identities(metadata_path: Path, structure_path: Path, binder_chains: list[str]) -> tuple[dict[str, str], int]:
    data = json.loads(metadata_path.read_text())
    if str(data.get("condition") or "").lower() != "sequence_fixed":
        return {}, 0
    entries = data.get("fixed_residues") if isinstance(data.get("fixed_residues"), list) else []
    if not entries:
        return {}, 0
    chains = _structure_residues(structure_path)
    binder_ids = {str(chain) for chain in binder_chains}
    pose_map = _rfantibody_pose_residue_map(metadata_path, structure_path)
    identities: dict[str, str] = {}
    for entry in entries:
        position = entry.get("generated_position") if isinstance(entry, dict) else None
        if not isinstance(position, dict):
            raise ValueError("RFantibody MIMOSA fixed-residue metadata has no generated residue position.")
        chain = str(position.get("chain_id") or "")
        number = str(position.get("author_residue_id") or "").strip()
        insertion = str(position.get("insertion_code") or "").strip()
        amino_acid = str(entry.get("amino_acid") or "").strip().upper()
        if chain not in binder_ids or len(amino_acid) != 1:
            raise ValueError("RFantibody MIMOSA fixed-residue metadata does not identify a valid antibody residue.")
        pose_value = entry.get("zero_based_pose_index", entry.get("generated_pose_index"))
        mapped = None
        if pose_value is not None:
            try:
                mapped = pose_map.get(int(pose_value))
            except (TypeError, ValueError):
                mapped = None
            if mapped is None:
                raise ValueError(f"RFantibody MIMOSA pose index {pose_value!r} is absent from the saved QV structure.")
            mapped_chain, residue = mapped
            if mapped_chain != chain:
                raise ValueError(
                    f"RFantibody MIMOSA pose index {pose_value} maps to chain {mapped_chain}, not declared chain {chain}."
                )
            label = str(residue["label"])
        else:
            label = f"{chain}{number}{insertion}"
            observed = {str(row["label"]): row for row in chains.get(chain, [])}
            if label not in observed:
                raise ValueError(f"RFantibody MIMOSA fixed residue {label} is absent from the generated structure.")
        if label in identities:
            raise ValueError(f"RFantibody MIMOSA mapped two fixed identities to residue {label}.")
        identities[label] = amino_acid
    return identities, len(entries)


def _mimosa_metadata_path(candidate: dict[str, Any], source_run_dir: Path) -> tuple[str, Path | None]:
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    tool = str(candidate.get("source_tool") or "").lower()
    if "boltzgen" in tool:
        relative = str(raw.get("metadata_path") or "").strip()
        path = Path(relative)
        if relative and not path.is_absolute():
            path = source_run_dir / path
        return "boltzgen", path if relative and path.is_file() else None
    if "rfantibody" in tool:
        tag = str(raw.get("generation_tag") or "").strip()
        if tag and Path(tag).name == tag:
            stages = source_run_dir / "artifacts/raw/rfantibody/stages"
            for path in (stages / f"1_rfdiffusion_{tag}.mimosa.json", stages / f"{tag}.mimosa.json"):
                if path.is_file():
                    return "rfantibody", path
    return "", None


def hydrate_candidate_fixed_motif(candidate: dict[str, Any], source_run_dir: str | Path) -> dict[str, Any]:
    """Attach verified source labels and intended amino acids for sequence-fixed MIMOSA designs.

    The structure may carry placeholder residue names. Downstream sequence
    design must apply the intended identities before asking MPNN to freeze them.
    """
    row = copy.deepcopy(candidate)
    raw = dict(row.get("raw_metadata") or {})
    if "mimosa_fixed_residue_identities" in raw:
        return row
    run_dir = Path(source_run_dir)
    engine, metadata_path = _mimosa_metadata_path(row, run_dir)
    if not engine or metadata_path is None:
        return row
    structure_path = Path(str(row.get("complex_pdb") or row.get("binder_pdb") or ""))
    if not structure_path.is_absolute():
        structure_path = run_dir / structure_path
    if not structure_path.is_file():
        return row

    try:
        if engine == "boltzgen":
            identities, declared_count = _boltzgen_fixed_identities(
                metadata_path,
                structure_path,
                [str(value) for value in row.get("binder_chains") or []],
                [str(value) for value in row.get("target_chains") or []],
            )
        else:
            identities, declared_count = _rfantibody_fixed_identities(
                metadata_path, structure_path, [str(value) for value in row.get("binder_chains") or []]
            )
        required = declared_count > 0
        raw.update({
            "fixed_motif_required": required,
            "fixed_residues": list(identities),
            "motif_mapping_verified": len(identities) == declared_count,
            "mimosa_fixed_residue_identities": identities,
            "mimosa_fixed_residue_count": declared_count,
            "mimosa_fixed_residue_source": engine,
        })
        if len(identities) != declared_count:
            raw["motif_mapping_warning"] = "Not every sequence-fixed MIMOSA identity mapped to the generated antibody structure. AbMPNN is blocked until that mapping is repaired."
        else:
            observed = {label: residue["amino_acid"] for values in _structure_residues(structure_path).values() for residue in values for label in [residue["label"]]}
            raw["mimosa_fixed_residue_input_mismatches"] = [
                {"residue": label, "expected": amino_acid, "observed": observed.get(label)}
                for label, amino_acid in identities.items()
                if observed.get(label) != amino_acid
            ]
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
        raw.update({
            "fixed_motif_required": True,
            "fixed_residues": [],
            "motif_mapping_verified": False,
            "mimosa_fixed_residue_identities": {},
            "mimosa_fixed_residue_source": engine,
            "motif_mapping_warning": f"Could not map sequence-fixed MIMOSA identities: {type(exc).__name__}: {exc}",
        })
    row["raw_metadata"] = raw
    return row


def apply_fixed_motif_identities(
    candidate: dict[str, Any], staged_pdb: Path, source_to_engine_residues: dict[str, str]
) -> dict[str, str]:
    """Patch the MPNN staging copy to the MIMOSA identities before freezing them."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    if not raw.get("fixed_motif_required"):
        return {}
    if raw.get("motif_mapping_verified") is not True:
        raise ValueError(raw.get("motif_mapping_warning") or "Sequence redesign requires a verified MIMOSA fixed-motif mapping.")
    identities = raw.get("mimosa_fixed_residue_identities")
    if not isinstance(identities, dict) or not identities:
        raise ValueError("Sequence-fixed MIMOSA metadata has no mapped amino-acid identities.")

    from mn_molecular_design_studio.mpnn_input import inventory

    info = inventory(staged_pdb)
    patches: dict[str, str] = {}
    for source_label, amino_acid in identities.items():
        label = source_to_engine_residues.get(str(source_label), str(source_label))
        residue = info["residues"].get(label)
        if not residue or residue.get("kind") != "protein":
            raise ValueError(f"Could not apply fixed MIMOSA residue {source_label} in the AbMPNN staging structure.")
        aa = str(amino_acid).upper()
        if len(aa) != 1 or aa not in "ACDEFGHIKLMNPQRSTVWY":
            raise ValueError(f"Invalid fixed MIMOSA amino acid {amino_acid!r} at {source_label}.")
        patches[label] = aa

    aa3 = {
        "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
        "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
        "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
        "W": "TRP", "Y": "TYR",
    }
    records = {}
    labels_by_record = {}
    for label, aa in patches.items():
        residue = info["residues"][label]
        number = int(residue["number"])
        insertion = str(residue.get("insertion_code") or "")
        chain = label[: -(len(str(number)) + len(insertion))]
        key = (chain, number, insertion)
        records[key] = aa3[aa]
        labels_by_record[key] = label

    lines = staged_pdb.read_text(errors="ignore").splitlines()
    patched = set()
    output = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27:
            chain = line[21].strip() or "_"
            try:
                number = int(line[22:26])
            except ValueError:
                number = None
            insertion = line[26].strip()
            key = (chain, number, insertion)
            if key in records:
                line = line[:17] + records[key] + line[20:]
                patched.add(key)
        output.append(line)
    if len(patched) != len(records):
        missing = [labels_by_record[key] for key in set(records) - patched]
        raise ValueError("Could not write all fixed MIMOSA amino-acid identities into the AbMPNN staging structure: " + ", ".join(missing or patches))
    staged_pdb.write_text("\n".join(output) + "\n")
    return patches
