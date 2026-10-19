"""Structure, residue, and compound inputs for the AutoDock Vina form."""

from __future__ import annotations

import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

from mn_molecular_design_studio import ligand_compounds
from mn_molecular_design_studio.auth import user_input_root

EXAMPLE_ROOT = Path(__file__).parent / "assets" / "vina"
_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL", "HID", "HIE", "HIP", "CYX", "CYM",
}
_MODIFIED = {"MSE": "M", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C"}


def protein_sequence_entries(pdb_text: str) -> dict[str, list[tuple[int, str]]]:
    chains: dict[str, list[tuple[int, str]]] = {}
    seen: set[tuple[str, int]] = set()
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")) or line[17:20].strip() not in _AMINO_ACIDS:
            continue
        try:
            residue = int(line[22:26])
        except ValueError:
            continue
        chain = line[21].strip() or "_"
        if (chain, residue) not in seen:
            seen.add((chain, residue))
            letter = seq1(line[17:20].strip(), custom_map=_MODIFIED)
            chains.setdefault(chain, []).append((residue, letter))
    return chains


def parse_residue_ranges(text: str, valid: set[int]) -> set[int]:
    selected: set[int] = set()
    text = re.sub(r"(?<=\d)\s*([-–:])\s*(?=-?\d)", r"\1", text)
    for part in re.split(r"[,;\s]+", text.strip()):
        if not part:
            continue
        match = re.fullmatch(r"(-?\d+)(?:\s*[-–:]\s*(-?\d+))?", part)
        if not match:
            raise ValueError(f"Invalid residue selection {part!r}. Enter numbers or ranges, for example 225, 230-235.")
        start, end = int(match[1]), int(match[2] or match[1])
        if end < start or end - start > 100_000:
            raise ValueError(f"Invalid residue range {part!r}.")
        requested = set(range(start, end + 1))
        if not requested <= valid:
            missing = sorted(requested - valid)
            raise ValueError(f"Residues {', '.join(map(str, missing[:8]))} are not present in this chain.")
        selected.update(requested)
    return selected


def residue_range_text(residues: set[int]) -> str:
    ranges = []
    ordered = sorted(residues)
    for value in ordered:
        if ranges and value == ranges[-1][1] + 1:
            ranges[-1][1] = value
        else:
            ranges.append([value, value])
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in ranges)


def stage_compound_input(owner_id: str, name: str, content: bytes, *, project_id: str = "") -> dict:
    compounds = ligand_compounds.parse_compound_table(content, filename=name)
    return _stage_compounds(owner_id, name, content, compounds, project_id=project_id)


def _stage_compounds(owner_id: str, name: str, content: bytes, compounds: list[dict], *, project_id: str = "") -> dict:
    digest = hashlib.sha256(content).hexdigest()[:20]
    folder = user_input_root(owner_id) / "ligand_tools" / "compound_inputs" / digest
    folder.mkdir(parents=True, exist_ok=True, mode=0o700)
    source = folder / Path(name).name
    source.write_bytes(content)
    source.chmod(0o600)
    path = folder / "compound_set.json"
    payload = ligand_compounds.write_compound_set(
        path, name=Path(name).stem, compounds=compounds, owner_id=owner_id,
        project_id=project_id, created_at=datetime.now(timezone.utc).isoformat(),
    )
    path.chmod(0o600)
    return {"path": str(path), "source_path": str(source), "name": name,
            "compound_count": payload["compound_count"], "source_sha256": digest}


def stage_example(owner_id: str, *, project_id: str = "") -> tuple[Path, dict]:
    folder = user_input_root(owner_id) / "ligand_tools" / "examples"
    folder.mkdir(parents=True, exist_ok=True, mode=0o700)
    target = folder / "1iep.pdb"
    target.write_bytes((EXAMPLE_ROOT / "1iep.pdb").read_bytes())
    target.chmod(0o600)
    ligand_name = "Conformer3D_COMPOUND_CID_5291.sdf"
    compound = _stage_compounds(
        owner_id, ligand_name, (EXAMPLE_ROOT / ligand_name).read_bytes(),
        json.loads((EXAMPLE_ROOT / "imatinib.json").read_text()), project_id=project_id,
    )
    return target, compound


def predicted_contact_residues(target_text: str, predicted_text: str, contact_distance: float) -> set[tuple[str, int]]:
    """Map predicted probe contacts back to author residues in the uploaded receptor.

    Boltz assigns protein chains A, B, ... in input order and renumbers their
    residues. Sequence alignment maps those contacts into the original frame;
    predicted Cartesian coordinates are never used as the Vina box center.
    """
    if not math.isfinite(contact_distance) or not 1 <= contact_distance <= 20:
        raise ValueError("Contact distance must be between 1 and 20 Å.")
    originals = list(protein_sequence_entries(target_text).items())
    predictions = list(protein_sequence_entries(predicted_text).items())
    if len(predictions) != len(originals):
        raise ValueError("The predicted complex does not contain the same protein chains as the Vina target.")
    atom_coordinates: dict[tuple[str, int], list[tuple[float, float, float]]] = {}
    ligand_points = []
    for line in predicted_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else line[12:16].strip()[:1].upper()
        if element in {"H", "D"}:
            continue
        try:
            number = int(line[22:26])
            xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
        except ValueError:
            continue
        resname = line[17:20].strip()
        if resname in _AMINO_ACIDS:
            atom_coordinates.setdefault((line[21].strip() or "_", number), []).append(xyz)
        elif resname not in {"HOH", "WAT", "DOD", "NA", "CL", "MG", "CA", "ZN", "K"}:
            ligand_points.append(xyz)
    if not ligand_points:
        raise ValueError("No probe-ligand heavy atoms were found in the predicted complex.")
    ligand = np.asarray(ligand_points)
    result: set[tuple[str, int]] = set()
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-2, open_gap_score=-6, extend_gap_score=-1)
    for (original_chain, original), (predicted_chain, predicted) in zip(originals, predictions):
        original_seq = "".join(letter for _, letter in original)
        predicted_seq = "".join(letter for _, letter in predicted)
        alignment = aligner.align(original_seq, predicted_seq)[0]
        mapping = {}
        matches = 0
        for (oa, ob), (pa, pb) in zip(*alignment.aligned):
            for oi, pi in zip(range(oa, ob), range(pa, pb)):
                mapping[predicted[pi][0]] = original[oi][0]
                matches += original_seq[oi] == predicted_seq[pi]
        if matches / max(len(original_seq), len(predicted_seq), 1) < .9:
            raise ValueError("A predicted protein chain does not match the Vina target sequence.")
        for predicted_number, original_number in mapping.items():
            points = atom_coordinates.get((predicted_chain, predicted_number), [])
            if points and np.any(np.sum((np.asarray(points)[:, None, :] - ligand[None, :, :]) ** 2, axis=-1) <= contact_distance ** 2):
                result.add((original_chain, original_number))
    return result
