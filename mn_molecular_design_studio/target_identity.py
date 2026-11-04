"""Helpers for distinguishing target sequence identity from target structures."""
from __future__ import annotations

from pathlib import Path


_STRUCTURE_KEYS = ("target_structure", "target_pdb", "protein_pdb", "complex_pdb")
_SEQUENCE_KEYS = ("target_sequence", "target_protein_sequences", "target_fasta")


def _mappings(input_data: dict | None, metadata: dict | None = None) -> list[dict]:
    input_data = input_data if isinstance(input_data, dict) else {}
    metadata = metadata if isinstance(metadata, dict) else {}
    return [metadata, input_data, input_data.get("inputs") or {}, input_data.get("params") or {}]


def is_sequence_carrier_reference(value, *, run_dir: Path | str | None = None) -> bool:
    """Recognize the RFD4 CIF generated only to carry an IDR sequence.

    The filename is an explicit marker emitted by Studio; when the referenced
    file is available, also verify its CIF data block before classifying it.
    """
    text = str(value or "").strip()
    if not text:
        return False
    filename = text.replace("\\", "/").split("?", 1)[0].rsplit("/", 1)[-1].casefold()
    named_carrier = filename == "idr_sequence_carrier.cif"
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
        ):
            candidate = root / relative
            if candidate.is_file():
                path = candidate
                break
    if path is None:
        # Historical/remote job metadata can retain only the portable URI.
        return named_carrier
    try:
        return path.is_file() and "data_idr_sequence_carrier" in path.read_text(
            encoding="utf-8", errors="ignore")[:4096].casefold()
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
