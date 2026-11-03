"""Reference-sequence context for a validated Complex Analysis placement.

This module formats polymer sequence and residue-presence provenance from the
selected Database reference. It does not calculate placement, clashes,
accessibility, glycan occupancy, or other biological context measurements.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


def _values(value: Any) -> list[str]:
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [str(value)]


def _aa_map() -> dict[str, str]:
    try:
        from Bio.Data.PDBData import protein_letters_3to1_extended
        return {str(key).upper(): str(value) for key, value in protein_letters_3to1_extended.items()}
    except ImportError:
        from Bio.Data.IUPACData import protein_letters_3to1_extended
        return {str(key).upper(): str(value) for key, value in protein_letters_3to1_extended.items()}


def _clean(value: Any) -> str:
    value = str(value or "").strip()
    return "" if value in {".", "?"} else value


def _integer(value: Any) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError, OverflowError):
        return None


def _cif_polymer_rows(path: Path, chain_id: str) -> list[dict[str, Any]]:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    data = MMCIF2Dict(str(path))
    keys = (
        "_pdbx_poly_seq_scheme.pdb_strand_id",
        "_pdbx_poly_seq_scheme.asym_id",
        "_pdbx_poly_seq_scheme.seq_id",
        "_pdbx_poly_seq_scheme.auth_seq_num",
        "_pdbx_poly_seq_scheme.pdb_ins_code",
        "_pdbx_poly_seq_scheme.mon_id",
    )
    columns = [_values(data.get(key)) for key in keys]
    if not columns or not columns[0]:
        return []
    if len({len(column) for column in columns}) != 1:
        return []
    authors, label_chains, sequence_ids, author_numbers, insertions, monomers = columns
    letters = _aa_map()
    rows = []
    for auth_chain, label_chain, sequence_id, number, insertion, monomer in zip(
        authors, label_chains, sequence_ids, author_numbers, insertions, monomers,
    ):
        author_chain_ids = {part.strip() for part in str(auth_chain).split(",") if part.strip()}
        if chain_id not in author_chain_ids:
            continue
        aa = letters.get(str(monomer).strip().upper(), "")
        seq_id = _integer(sequence_id)
        if not aa or seq_id is None:
            continue
        author_number = _integer(_clean(number))
        insertion_code = _clean(insertion)
        rows.append({
            "sequence_index": seq_id - 1,
            "amino_acid": aa,
            "author_number": author_number,
            "insertion_code": insertion_code,
            "author_residue": (f"{author_number}{insertion_code}" if author_number is not None else f"seq#{seq_id}"),
            "label_chain_id": str(label_chain),
        })
    return sorted(rows, key=lambda row: row["sequence_index"])


def _pdb_seqres_rows(path: Path, chain_id: str) -> list[dict[str, Any]]:
    """Read a PDB SEQRES track and map observed author IDs only when unique.

    SEQRES has no author residue numbers for unresolved residues. Those tiles
    keep their sequence index and are explicitly labelled as unresolved.
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    letters = _aa_map()
    names: list[str] = []
    for line in text.splitlines():
        if not line.startswith("SEQRES") or (line[11:12].strip() or "_") != (chain_id.strip() or "_"):
            continue
        names.extend(line[19:70].split())
    sequence = "".join(letters.get(name.upper(), "") for name in names)
    if not sequence:
        return []

    from mn_molecular_design_studio.complex_analysis_inspector import protein_sequence_tracks
    observed_tracks = protein_sequence_tracks(path, [chain_id])
    observed = observed_tracks[0]["residues"] if observed_tracks else []
    observed_sequence = "".join(row["amino_acid"] for row in observed)
    positions: dict[int, dict[str, Any]] = {}
    if observed_sequence:
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                                  open_gap_score=-8, extend_gap_score=-0.5)
        aligner.target_left_gap_score = 0.0
        aligner.target_right_gap_score = 0.0
        alignments = aligner.align(sequence, observed_sequence)
        if len(alignments) == 1:
            alignment = alignments[0]
            ref_blocks, observed_blocks = alignment.aligned
            for (ref_start, ref_end), (obs_start, obs_end) in zip(ref_blocks, observed_blocks):
                for ref_index, obs_index in zip(range(int(ref_start), int(ref_end)), range(int(obs_start), int(obs_end))):
                    residue = observed[int(obs_index)]
                    if sequence[int(ref_index)] == observed_sequence[int(obs_index)]:
                        positions[int(ref_index)] = residue

    return [{
        "sequence_index": index,
        "amino_acid": aa,
        "author_number": positions.get(index, {}).get("author_number"),
        "insertion_code": positions.get(index, {}).get("insertion_code", ""),
        "author_residue": positions.get(index, {}).get("author_residue") or f"seq#{index + 1}",
        "label_chain_id": chain_id,
    } for index, aa in enumerate(sequence)]


def reference_component_track(
    path: Path | str,
    chain_id: str,
    residue_mapping: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build the full polymer sequence track for a selected reference chain.

    ``residue_mapping`` must come from the saved, validated placement result.
    Reference-only resolved residues are marked omitted from the source
    construct; sequence positions without atoms are marked unresolved.
    """
    path = Path(path)
    suffix = path.name.casefold()
    if suffix.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        if suffix.endswith(".gz"):
            import gzip
            unpacked = path.with_suffix("")
            temporary = path.with_name(path.name + ".parsed")
            temporary.write_bytes(gzip.decompress(path.read_bytes()))
            try:
                rows = _cif_polymer_rows(temporary, str(chain_id))
            finally:
                temporary.unlink(missing_ok=True)
        else:
            rows = _cif_polymer_rows(path, str(chain_id))
    else:
        rows = _pdb_seqres_rows(path, str(chain_id))
    if not rows:
        return {"chain_id": str(chain_id), "sequence": "", "residues": [],
                "status": "unavailable", "reason": "The selected reference file has no supported full polymer sequence and numbering map for this chain."}

    from mn_molecular_design_studio.complex_analysis_inspector import author_residue_selector, protein_sequence_tracks

    observed_tracks = protein_sequence_tracks(path, [str(chain_id)])
    observed = observed_tracks[0]["residues"] if observed_tracks else []
    resolved: dict[tuple[int, str], dict[str, Any]] = {}
    for residue in observed:
        key = (int(residue["author_number"]), str(residue.get("insertion_code") or ""))
        if key not in resolved:
            resolved[key] = residue

    mapped_keys = set()
    for item in residue_mapping:
        reference_residue = item.get("reference") if isinstance(item.get("reference"), dict) else {}
        number = _integer(reference_residue.get("number"))
        insertion = str(reference_residue.get("insertion") or "").strip()
        if number is not None:
            mapped_keys.add((number, insertion))

    counts = {"mapped": 0, "present": 0, "omitted": 0, "unresolved": 0}
    track_rows = []
    for row in rows:
        number = row.get("author_number")
        insertion = str(row.get("insertion_code") or "")
        match = resolved.get((int(number), insertion)) if number is not None else None
        if match and match["amino_acid"] != row["amino_acid"]:
            match = None
        is_mapped = bool(match and (int(number), insertion) in mapped_keys)
        if is_mapped:
            status = "mapped"
            counts["mapped"] += 1
            counts["present"] += 1
        elif match:
            status = "omitted"
            counts["omitted"] += 1
            counts["present"] += 1
        else:
            status = "unresolved"
            counts["unresolved"] += 1
        selector = author_residue_selector(str(chain_id), int(number), insertion) if match and number is not None else ""
        track_rows.append({
            **row,
            "status": status,
            "reference_coordinates_present": bool(match),
            "source_construct_present": is_mapped,
            "selector": selector,
        })
    return {
        "chain_id": str(chain_id),
        "sequence": "".join(row["amino_acid"] for row in track_rows),
        "residues": track_rows,
        "status": "available",
        "counts": counts,
        "numbering_source": "PDB/mmCIF author residue numbers and insertion codes where available; unresolved PDB SEQRES positions have no author number",
    }


def explicit_reference_components(path: Path | str) -> list[dict[str, Any]]:
    """Summarize declared polymer/branched entities from the exact reference.

    These are file-presence annotations only. A declared carbohydrate entity
    is not evidence of biological glycan occupancy or accessibility.
    """
    path = Path(path)
    if not path.name.casefold().endswith((".cif", ".mmcif")):
        return []
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    data = MMCIF2Dict(str(path))
    entity_ids = _values(data.get("_pdbx_entity_branch.entity_id"))
    branch_counts = _values(data.get("_pdbx_entity_branch_list.entity_id"))
    branch_count_by_entity: dict[str, int] = {}
    for entity in branch_counts:
        branch_count_by_entity[entity] = branch_count_by_entity.get(entity, 0) + 1
    asym_ids = _values(data.get("_struct_asym.id"))
    asym_entities = _values(data.get("_struct_asym.entity_id"))
    entity_description_ids = _values(data.get("_entity.id"))
    descriptions = _values(data.get("_entity.pdbx_description"))
    description_by_entity = dict(zip(entity_description_ids, descriptions))
    label_to_author: dict[str, str] = {}
    labels = _values(data.get("_atom_site.label_asym_id"))
    authors = _values(data.get("_atom_site.auth_asym_id"))
    for label, author in zip(labels, authors):
        label_to_author.setdefault(label, author)
    selected = set(entity_ids)
    return [{
        "entity_id": entity_id,
        "label_chain_id": asym,
        "author_chain_id": label_to_author.get(asym, asym),
        "description": description_by_entity.get(entity_id, "Branched component"),
        "residue_count": branch_count_by_entity.get(entity_id, 0),
    } for asym, entity_id in zip(asym_ids, asym_entities) if entity_id in selected]
