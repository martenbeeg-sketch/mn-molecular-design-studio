"""Helpers for distinguishing target sequence identity from target structures."""
from __future__ import annotations

from pathlib import Path
from collections import Counter
from functools import lru_cache
import hashlib
import json
import re


_STRUCTURE_KEYS = ("target_structure", "target_pdb", "protein_pdb", "complex_pdb")
_SEQUENCE_KEYS = ("target_sequence", "target_protein_sequences", "target_fasta")


@lru_cache(maxsize=1)
def _pdb_example_sequences() -> dict[str, Counter]:
    path = Path(__file__).parent / "components" / "cofolding_editor" / "examples.json"
    examples = json.loads(path.read_text(encoding="utf-8"))
    return {code: Counter(str(card.get("value") or "").upper() for card in value.get("cards") or []
                          if card.get("kind") == "protein")
            for code, value in examples.items() if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", code)}


def protein_sequence_reference(sequences: list[str]) -> dict:
    """Identify an exact, unique chain-sequence match to a bundled PDB example.

    This is reference evidence, not a claim that coordinates were provided to
    the prediction engine. Chain order and biological roles are never inferred.
    Copy counts matter, and changing even one residue invalidates the match.
    """
    normalized = ["".join(str(value).split()).upper() for value in sequences if value]
    if not normalized:
        return {}
    observed = Counter(normalized)
    matches = [(code, known) for code, known in _pdb_example_sequences().items() if observed <= known]
    if len(matches) != 1:
        return {}
    code, known = matches[0]
    return {"pdb_id": code.upper(), "basis": "exact_bundled_example_chain_sequences",
            "matched_chain_count": len(normalized), "reference_chain_count": sum(known.values()),
            "sequence_sha256": hashlib.sha256(":".join(sorted(normalized)).encode()).hexdigest()}


def prediction_input_systems(input_data: dict) -> list[dict]:
    """Describe saved protein-only prediction systems without inventing a target.

    An explicitly assigned target is handled by the ordinary target resolver.
    For inputs with unspecified roles, retain the whole molecular composition.
    """
    inputs = input_data.get("inputs") if isinstance(input_data.get("inputs"), dict) else input_data
    systems = {}
    for task in inputs.get("tasks") or []:
        if not isinstance(task, dict) or any(task.get(key) for key in
                ("dna_chains", "rna_chains", "ligand_ccd_codes", "ligand_smiles")):
            continue
        chains = task.get("protein_chains") or []
        if not chains or any(not isinstance(chain, dict) for chain in chains):
            continue
        if any(chain.get("role") in {"target", "binder"} for chain in chains):
            continue
        sequences = ["".join(str(chain.get("sequence") or "").split()).upper() for chain in chains]
        if not all(sequences):
            continue
        reference = protein_sequence_reference(sequences)
        digest = hashlib.sha256(":".join(sorted(sequences)).encode()).hexdigest()
        identity = "sequence:" + digest
        systems[identity] = {
            "identity": identity, "identity_kind": "system", "sequence_length": sum(map(len, sequences)),
            "chain_count": len(chains), "label": (reference["pdb_id"] + " system" if reference else
                                                    f"Protein system · {len(chains)} chain(s) · {digest[:8]}"),
            "reference_type": "structure_backed" if reference else "sequence_only",
            "reference_role": "sequence_source" if reference else "prediction_sequence_input",
            "pdb_id": reference.get("pdb_id", ""), "reference_basis": reference.get("basis", ""),
        }
    return list(systems.values())


def _mappings(input_data: dict | None, metadata: dict | None = None) -> list[dict]:
    input_data = input_data if isinstance(input_data, dict) else {}
    metadata = metadata if isinstance(metadata, dict) else {}
    return [metadata, input_data, input_data.get("inputs") or {}, input_data.get("params") or {}]


def is_sequence_carrier_reference(value, *, run_dir: Path | str | None = None) -> bool:
    """Recognize Studio-generated CIFs that carry sequences, not structures.

    RFD4 has emitted both ``idr_sequence_carrier.cif`` and
    ``sequence_input.cif``. Their placeholder coordinates are computational
    conditioning artifacts and are never source/native structure evidence.
    """
    text = str(value or "").strip()
    if not text:
        return False
    filename = text.replace("\\", "/").split("?", 1)[0].rsplit("/", 1)[-1].casefold()
    carrier_markers = {
        "idr_sequence_carrier.cif": ("data_idr_sequence_carrier",),
        "sequence_input.cif": ("data_sequence", "data_studio_sequence_carrier"),
    }
    marker_set = carrier_markers.get(filename, ())
    named_carrier = bool(marker_set)
    path = None
    try:
        from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

        resolved = resolve_stored_path(text, run_dir=run_dir, must_exist=True) if run_dir else None
        if resolved:
            path = Path(resolved)
    except (OSError, TypeError, ValueError):
        path = None
    if path is None and run_dir and named_carrier:
        root = Path(run_dir)
        for relative in (
            "artifacts/rfd4_proteina/input/idr_sequence_carrier.cif",
            "artifacts/input/idr_sequence_carrier.cif",
            "artifacts/sequence_input.cif",
            "sequence_input.cif",
        ):
            candidate = root / relative
            if candidate.is_file():
                path = candidate
                break
    if path is None:
        # Historical/remote job metadata can retain only the portable URI.
        return named_carrier
    try:
        content = path.read_text(encoding="utf-8", errors="ignore")[:4096].casefold()
        return path.is_file() and any(marker in content for marker in marker_set)
    except OSError:
        return False


def sequence_carrier_path(input_data: dict, *, run_dir: Path | str) -> Path | None:
    """Resolve the sequence carrier staged with a completed RFD4 run."""
    for mapping in _mappings(input_data):
        if not isinstance(mapping, dict):
            continue
        for key in _STRUCTURE_KEYS:
            value = mapping.get(key)
            if not is_sequence_carrier_reference(value, run_dir=run_dir):
                continue
            try:
                from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

                resolved = resolve_stored_path(value, run_dir=run_dir, must_exist=True)
            except (OSError, TypeError, ValueError):
                resolved = None
            if resolved and Path(resolved).is_file():
                return Path(resolved)
            for relative in (
                "artifacts/rfd4_proteina/input/idr_sequence_carrier.cif",
                "artifacts/input/idr_sequence_carrier.cif",
                "artifacts/sequence_input.cif",
                "sequence_input.cif",
            ):
                candidate = Path(run_dir) / relative
                if candidate.is_file():
                    return candidate
    return None


def target_reference_type(input_data: dict, metadata: dict | None = None, *, run_dir=None) -> str:
    """Return sequence_only, structure_backed, or unknown from saved inputs."""
    mappings = _mappings(input_data, metadata)
    saw_structure = False
    for mapping in mappings:
        if not isinstance(mapping, dict):
            continue
        for key in _STRUCTURE_KEYS:
            value = mapping.get(key)
            if not value:
                continue
            saw_structure = True
            if is_sequence_carrier_reference(value, run_dir=run_dir):
                return "sequence_only"
    if saw_structure:
        return "structure_backed"
    for mapping in mappings:
        if not isinstance(mapping, dict):
            continue
        if any(mapping.get(key) for key in _SEQUENCE_KEYS):
            return "sequence_only"
    return "unknown"
