"""Coordinate coverage for structure import review and the structure preview API."""
from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from itertools import pairwise
from pathlib import Path

from Bio.PDB.Polypeptide import is_aa, is_nucleic


MAX_PREVIEW_BYTES = 50_000_000


def _cif_annotations(metadata, model_number):
    prefix = "_pdbx_unobs_or_zero_occ_residues."
    count = len(metadata.get(prefix + "id", ()))
    columns = ("auth_asym_id", "auth_seq_id", "PDB_ins_code", "PDB_model_num", "polymer_flag", "occupancy_flag")
    for chain, number, insertion, model, polymer, occupancy in zip(
        *(metadata.get(prefix + column, ["?"] * count) for column in columns), strict=True
    ):
        if polymer != "Y" or model not in {str(model_number), ".", "?"} or occupancy not in {"0", "1"}:
            continue
        try:
            position = (int(number), "" if insertion in {".", "?"} else insertion)
        except ValueError:
            continue
        yield "_" if chain in {".", "?"} else chain, position, "unobserved" if occupancy == "1" else "zero_occupancy"


def _pdb_annotations(metadata, model_number):
    for row in metadata.get("missing_residues", ()):
        if row["model"] in {None, model_number}:
            yield row["chain"], (row["ssseq"], row["insertion"] or ""), "unobserved"


def _position(first, last, bounds):
    if bounds is None:
        return "undetermined"
    if last < bounds[0]:
        return "before_observed"
    if first > bounds[1]:
        return "after_observed"
    return "internal"


def _segment(chain, first, last, count, kind, position, source):
    return {"chain": chain, "first": f"{first[0]}{first[1]}", "last": f"{last[0]}{last[1]}",
            "count": count, "kind": kind, "position": position, "source": source}


def _annotated_segments(chain, records, observed):
    bounds = (min(observed), max(observed)) if observed else None
    groups = []
    for site, (kind, source) in sorted(records.items()):
        if groups:
            first, last, count, previous_kind, previous_source = groups[-1]
            if (kind, source) == (previous_kind, previous_source) and not (site[1] or last[1]) and site[0] == last[0] + 1:
                groups[-1] = (first, site, count + 1, kind, source)
                continue
        groups.append((site, site, 1, kind, source))
    return [_segment(chain, first, last, count, kind, _position(first, last, bounds), source)
            for first, last, count, kind, source in groups]


def inspect_coordinate_coverage(data: bytes, filename: str):
    """Report first-model coverage by author residue number using the original file."""
    from mn_molecular_design_studio.molecule_imports import read_coordinate_file

    structure, metadata = read_coordinate_file(data, filename)
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("No coordinate models were found.")
    number = model.serial_num or 1
    cif = Path(filename).suffix.lower() in {".cif", ".mmcif"}
    observed, records = defaultdict(set), defaultdict(dict)
    for chain in model:
        chain_id = chain.id.strip() or "_"
        for residue in chain:
            if not (is_aa(residue) or is_nucleic(residue)):
                continue
            site = (residue.id[1], residue.id[2].strip())
            if all(atom.occupancy == 0 for atom in residue.get_unpacked_list()):
                records[chain_id][site] = ("zero_occupancy", "Atom occupancies")
            else:
                observed[chain_id].add(site)

    annotations = _cif_annotations(metadata, number) if cif else _pdb_annotations(metadata, number)
    source = "mmCIF deposition" if cif else "PDB REMARK 465"
    for chain, site, kind in annotations:
        if site not in observed.get(chain, ()):
            records[chain][site] = (kind, source)

    result = {"model_number": number, "chains": [], "segments": []}
    for chain in sorted(observed.keys() | records.keys()):
        observed_sites, annotated_sites = observed[chain], records[chain]
        result["chains"].append({"chain": chain, "observed_count": len(observed_sites)})
        result["segments"].extend(_annotated_segments(chain, annotated_sites, observed_sites))
        annotated_numbers = sorted({site[0] for site in annotated_sites})
        for left, right in pairwise(sorted({site[0] for site in observed_sites})):
            if right <= left + 1:
                continue
            index = bisect_right(annotated_numbers, left)
            if index < len(annotated_numbers) and annotated_numbers[index] < right:
                continue
            result["segments"].append(_segment(chain, (left + 1, ""), (right - 1, ""), None,
                                               "numbering_jump", "internal", "Coordinate numbering"))
    return result


def select_coverage_chains(report, chains):
    if chains is None:
        return report
    keep = set(chains)
    return {**report, "chains": [row for row in report["chains"] if row["chain"] in keep],
            "segments": [row for row in report["segments"] if row["chain"] in keep]}
