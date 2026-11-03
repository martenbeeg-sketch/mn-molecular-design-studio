"""Data helpers for the Database Complex Analysis scientific Inspector.

The Inspector only formats saved structure and Evaluation evidence. It does
not calculate new scientific measurements.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any


_CHAIN_ID = re.compile(r"^[^:\s^]+$")


def author_residue_selector(chain_id: str, number: int, insertion_code: str = "") -> str:
    """Encode one exact author residue for Mol*; ``^_`` denotes blank insertion."""
    chain_id = str(chain_id or "_").strip() or "_"
    if not _CHAIN_ID.fullmatch(chain_id):
        return ""
    insertion_code = str(insertion_code or "").strip()
    if len(insertion_code) > 1 or (insertion_code and not re.fullmatch(r"[A-Za-z0-9]", insertion_code)):
        return ""
    try:
        number = int(number)
    except (TypeError, ValueError, OverflowError):
        return ""
    return f"{chain_id}:{number}^{insertion_code or '_'}"


def parse_author_residue_selector(value: str) -> dict[str, Any] | None:
    """Parse the Inspector's exact selector and legacy Mol* numeric selectors."""
    value = str(value or "").strip()
    match = re.fullmatch(r"([^:\s^]+):(-?\d+)(?:\^([A-Za-z0-9_])|([A-Za-z]))?", value)
    if not match:
        match = re.fullmatch(r"([^:\s^]+)(-?\d+)", value)
    if not match:
        return None
    chain = "" if match.group(1) == "_" else match.group(1)
    insertion = (match.group(3) or match.group(4)) if ":" in value and match.lastindex and match.lastindex >= 3 else None
    if insertion == "_":
        insertion = ""
    return {
        "chain_id": chain,
        "number": int(match.group(2)),
        "insertion_code": insertion,
        "exact_insertion": insertion is not None,
    }


def protein_sequence_tracks(path: Path | str, chain_ids: list[str]) -> list[dict[str, Any]]:
    """Return per-chain protein tracks with author residue IDs and selectors."""
    try:
        from Bio.Data.PDBData import protein_letters_3to1_extended
    except ImportError:
        from Bio.Data.IUPACData import protein_letters_3to1_extended as _letters
        protein_letters_3to1_extended = {key.upper(): value for key, value in _letters.items()}
        protein_letters_3to1_extended.setdefault("MSE", "M")

    from mn_molecular_design_studio.complex_analysis_comparison import _read_biopython_structure

    structure = _read_biopython_structure(path)
    models = list(structure.get_models())
    if len(models) != 1:
        return []
    model = models[0]
    tracks = []
    for chain_id in dict.fromkeys(map(str, chain_ids)):
        if chain_id not in model:
            continue
        residues = []
        for residue in model[chain_id]:
            amino_acid = protein_letters_3to1_extended.get(str(residue.resname).strip().upper())
            if not amino_acid:
                continue
            number = int(residue.id[1])
            insertion = str(residue.id[2] or "").strip()
            residues.append({
                "sequence_index": len(residues),
                "amino_acid": str(amino_acid),
                "author_number": number,
                "insertion_code": insertion,
                "author_residue": f"{number}{insertion}",
                "selector": author_residue_selector(chain_id, number, insertion),
            })
        if residues:
            tracks.append({"chain_id": chain_id, "sequence": "".join(row["amino_acid"] for row in residues),
                           "residues": residues})
    return tracks


def resolve_author_residue_selector(
    tracks: list[dict[str, Any]], chain_id: str, number: int | str,
    insertion_code: str | None = None,
) -> str:
    """Return a selector only when artifact residue identity resolves uniquely."""
    try:
        number = int(number)
    except (TypeError, ValueError, OverflowError):
        return ""
    candidates = []
    for track in tracks:
        if str(track.get("chain_id") or "") != str(chain_id or ""):
            continue
        for residue in track.get("residues") or []:
            if int(residue.get("author_number")) != number:
                continue
            if insertion_code is not None and str(residue.get("insertion_code") or "") != str(insertion_code or ""):
                continue
            candidates.append(str(residue.get("selector") or ""))
    return candidates[0] if len(candidates) == 1 else ""


def structure_digest(path: Path | str) -> str:
    """Hash the selected saved representation's exact bytes."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def observation_context(observation: dict) -> dict[str, Any]:
    value = observation.get("context") or {}
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except (TypeError, ValueError):
            return {}
    return value if isinstance(value, dict) else {}


def _selector_for_reference(reference: str, tracks: list[dict[str, Any]]) -> str:
    parsed = parse_author_residue_selector(reference)
    if not parsed:
        return ""
    return resolve_author_residue_selector(
        tracks, parsed["chain_id"], parsed["number"],
        parsed["insertion_code"] if parsed["exact_insertion"] else None,
    )


def inspector_evidence_rows(
    evaluator: str,
    artifact_name: str,
    artifact_rows: list[dict[str, Any]],
    tracks: list[dict[str, Any]],
    *,
    ligand_selector: str = "",
) -> list[dict[str, Any]]:
    """Build clickable, provenance-preserving rows from registered artifacts."""
    evaluator = str(evaluator or "unknown")
    artifact_name = str(artifact_name or "")
    output = []
    for index, raw in enumerate(artifact_rows):
        row = dict(raw)
        selectors: list[str] = []
        residue_refs: list[str] = []
        interaction_type = str(row.get("interaction_type") or "").strip()
        protein_chain = str(row.get("protein_chain") or row.get("chain_id") or row.get("chain") or "").strip()
        protein_number = row.get("protein_residue_number", row.get("residue_number", row.get("residue_no")))
        protein_insertion = row.get("protein_insertion_code")

        if artifact_name in {"contacts.json", "salt_bridges.json"}:
            for field in ("binder_residue", "target_residue"):
                reference = str(row.get(field) or "").strip()
                if not reference:
                    continue
                residue_refs.append(reference)
                selector = _selector_for_reference(reference, tracks)
                if selector:
                    selectors.append(selector)
            kind = "Salt bridge" if artifact_name == "salt_bridges.json" else "Interface contact"
            label = kind + " · " + " ↔ ".join(residue_refs or ["saved pair"])
        elif evaluator == "posebusters":
            if ligand_selector:
                selectors.append(ligand_selector)
            label = "PoseBusters ligand pose checks"
        elif protein_number not in (None, "") and protein_chain:
            protein_ref = f"{protein_chain}:{protein_number}"
            if protein_insertion not in (None, ""):
                protein_ref += f"^{str(protein_insertion).strip()}"
            residue_refs.append(protein_ref)
            selector = resolve_author_residue_selector(
                tracks, protein_chain, protein_number,
                str(protein_insertion).strip() if protein_insertion not in (None, "") else None,
            )
            if selector:
                selectors.append(selector)
            if ligand_selector and evaluator in {"plip", "pandamap"}:
                selectors.append(ligand_selector)
            label = interaction_type or str(row.get("patch_type") or row.get("finding") or "Saved residue observation")
            if residue_refs:
                label += " · " + ", ".join(residue_refs)
        else:
            # Recognize an explicit evaluator residue reference without
            # guessing a chain or renumbering it.
            reference = str(row.get("residue_ref") or row.get("residue_key") or "").strip()
            if reference:
                residue_refs.append(reference)
                selector = _selector_for_reference(reference, tracks)
                if selector:
                    selectors.append(selector)
            if ligand_selector and evaluator == "posebusters":
                selectors.append(ligand_selector)
            label = interaction_type or str(row.get("finding") or row.get("metric") or artifact_name)

        output.append({
            "id": f"{evaluator}:{artifact_name}:{index}",
            "evaluator": evaluator,
            "artifact": artifact_name,
            "label": label,
            "residue_refs": residue_refs,
            "selectors": list(dict.fromkeys(selector for selector in selectors if selector)),
            "raw": row,
        })
    return output
