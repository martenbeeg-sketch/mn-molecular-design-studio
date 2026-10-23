"""Small input normalizers shared by MIMOSA design pages and workflows."""

from __future__ import annotations

import re


def normalize_mimosa_source_positions(value: object) -> list[str]:
    """Normalize source-complex residue labels such as ``21,22,26A`` or ``H:33``.

    Chain-qualified IDs support motifs that span multiple source partner chains.
    The positions are checked against contacting residues when the motif packet
    is prepared. The UI selects these labels from the uploaded structure.
    """
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        tokens = [str(item).strip().upper() for item in value if str(item).strip()]
    else:
        tokens = [item.upper() for item in re.split(r"[,;\s]+", str(value).strip()) if item]
    normalized: list[str] = []
    for token in tokens:
        match = re.fullmatch(r"(?:([A-Z0-9]+):)?(-?\d{1,5})([A-Z]?)", token)
        if not match:
            raise ValueError(
                f"MIMOSA source positions must be author residue IDs, optionally chain-qualified and with an insertion letter (for example 21, 26A, H:33); got {token!r}."
            )
        chain, number, insertion = match.groups()
        label = f"{chain + ':' if chain else ''}{int(number)}{insertion}"
        if label in normalized:
            raise ValueError(f"MIMOSA source position {label} was selected more than once.")
        normalized.append(label)
    return normalized
