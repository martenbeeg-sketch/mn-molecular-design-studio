"""Optional component-to-structure reference context for Complex Analysis.

Reference associations are user annotations over saved Database structures.
They are deliberately separate from molecular identity, representations, and
scientific Evaluation records.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any


ALIGNMENT_PARAMETERS = {
    "mode": "global",
    "match_score": 2.0,
    "mismatch_score": -1.0,
    "open_gap_score": -8.0,
    "extend_gap_score": -0.5,
    "reference_terminal_gap_score": 0.0,
    "sequence_roles": {"target": "selected reference chain", "query": "selected source component"},
}


def _single_model(path: Path) -> None:
    from mn_molecular_design_studio.complex_analysis_comparison import _read_biopython_structure

    structure = _read_biopython_structure(path)
    model_count = sum(1 for _ in structure.get_models())
    if model_count != 1:
        raise ValueError(
            f"Reference mapping requires one model per selected Database representation; found {model_count} in {path.name}."
        )


def _inventory(path: Path, chain_id: str) -> dict[str, Any]:
    from mn_molecular_design_studio.complex_analysis_comparison import protein_chain_numbering

    rows = protein_chain_numbering(path, [chain_id])
    if len(rows) != 1:
        raise ValueError(f"Protein chain {chain_id!r} is not present in {path.name}.")
    row = rows[0]
    sequence = str(row.get("sequence") or "")
    numbering = list(row.get("numbering") or [])
    if not sequence or len(sequence) != len(numbering):
        raise ValueError(f"Could not preserve residue numbering for protein chain {chain_id!r} in {path.name}.")
    return {
        "chain_id": str(chain_id),
        "sequence": sequence,
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "residue_numbering": "PDB/mmCIF author chain and author residue number with insertion code",
        "residues": [
            {
                "sequence_index": index,
                "number": int(item["number"]),
                "insertion": str(item.get("insertion") or ""),
                "amino_acid": sequence[index],
            }
            for index, item in enumerate(numbering)
        ],
    }


def component_reference_mapping(
    source_path: Path | str,
    source_chain: str,
    reference_path: Path | str,
    reference_chain: str,
) -> dict[str, Any]:
    """Map a selected protein component to a selected reference chain.

    The mapping is sequence-derived and retains both source and reference
    author residue identifiers. It is marked sequence-validated only when the
    unique alignment maps every source residue to the same amino acid. This
    function does not fit coordinates or calculate a context score.
    """
    source_path, reference_path = Path(source_path), Path(reference_path)
    supported = (".pdb", ".ent", ".cif", ".mmcif", ".pdb.gz", ".ent.gz", ".cif.gz", ".mmcif.gz")
    for label, path in (("source", source_path), ("reference", reference_path)):
        if not path.name.casefold().endswith(supported):
            raise ValueError(f"The selected {label} representation must be a PDB or mmCIF structure.")
    _single_model(source_path)
    _single_model(reference_path)
    source = _inventory(source_path, str(source_chain))
    reference = _inventory(reference_path, str(reference_chain))

    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                              open_gap_score=-8, extend_gap_score=-0.5)
    # Free terminal overhangs on the selected reference chain let a shorter
    # construct map into a full-length chain. Internal deletions remain
    # explicit alignment gaps. The source component is always aligned in full.
    aligner.target_left_gap_score = 0.0
    aligner.target_right_gap_score = 0.0
    alignments = aligner.align(reference["sequence"], source["sequence"])
    alignment_count = len(alignments)
    mapping: list[dict[str, Any]] = []
    if alignment_count == 1:
        alignment = alignments[0]
        reference_blocks, source_blocks = alignment.aligned
        for (reference_start, reference_end), (source_start, source_end) in zip(reference_blocks, source_blocks):
            for reference_index, source_index in zip(
                range(int(reference_start), int(reference_end)),
                range(int(source_start), int(source_end)),
            ):
                source_residue = source["residues"][source_index]
                reference_residue = reference["residues"][reference_index]
                mapping.append({
                    "source": dict(source_residue),
                    "reference": dict(reference_residue),
                    "sequence_match": source_residue["amino_acid"] == reference_residue["amino_acid"],
                })

    mapped_source = {int(item["source"]["sequence_index"]) for item in mapping}
    mapped_reference = sorted(int(item["reference"]["sequence_index"]) for item in mapping)
    exact_matches = sum(bool(item["sequence_match"]) for item in mapping)
    source_count = len(source["residues"])
    contiguous_reference = bool(mapped_reference) and mapped_reference == list(
        range(mapped_reference[0], mapped_reference[0] + len(mapped_reference))
    )
    if alignment_count != 1:
        status = "ambiguous"
        reason = f"The selected chains have {alignment_count} equally optimal sequence mappings; residue correspondence is not inferred."
        mapping = []
    elif not mapping:
        status = "unmatched"
        reason = "The selected protein chains have no sequence correspondence."
    elif len(mapped_source) == source_count and exact_matches == source_count and contiguous_reference:
        status = "validated"
        reason = f"A unique exact contiguous subsequence mapping covers all {source_count} source residues."
    elif len(mapped_source) == source_count and exact_matches == source_count:
        status = "partial"
        reason = (
            f"All {source_count} source residues match exactly, but the mapped reference residues are not contiguous; "
            "this is not a truncated exact subsequence mapping."
        )
    else:
        status = "partial"
        reason = (
            f"The unique alignment maps {len(mapped_source)} of {source_count} source residues, "
            f"with {exact_matches} identical amino-acid positions; reference-dependent evaluations are disabled."
        )

    return {
        "mapping_method": "Biopython PairwiseAligner affine semiglobal sequence mapping",
        "mapping_parameters": dict(ALIGNMENT_PARAMETERS),
        "mapping_status": status,
        "mapping_reason": reason,
        "alignment_count": alignment_count,
        "source_component": source,
        "reference_component": reference,
        "residue_mapping": mapping,
        "mapped_source_residue_count": len(mapped_source),
        "exact_sequence_match_count": exact_matches,
        "source_residue_count": source_count,
        "reference_residue_count": len(reference["residues"]),
        "source_coverage": len(mapped_source) / source_count if source_count else 0.0,
        "mapped_sequence_identity": exact_matches / len(mapping) if mapping else 0.0,
    }
