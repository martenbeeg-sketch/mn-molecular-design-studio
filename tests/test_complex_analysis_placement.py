from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pytest

from mn_molecular_design_studio.complex_analysis_placement import calculate_reference_placement
from mn_molecular_design_studio.complex_analysis_reference import component_reference_mapping


_ROTATION = np.asarray([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
_TRANSLATION = np.asarray([10.0, -4.0, 2.0])
_BASE = np.asarray([
    [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 2.0, 0.0],
    [3.0, 2.0, 3.0], [5.0, 3.0, 2.0],
])


def _pdb(path: Path, chain: str, sequence: str, coordinates: list[list[float]], *,
         residue_numbers: list[tuple[int, str]] | None = None,
         missing_ca: set[int] | None = None) -> Path:
    three = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
             "G": "GLY", "H": "HIS", "S": "SER", "T": "THR"}
    residue_numbers = residue_numbers or [(index + 1, "") for index in range(len(sequence))]
    missing_ca = missing_ca or set()
    lines, serial = [], 1
    for index, (amino_acid, xyz, (number, insertion)) in enumerate(zip(sequence, coordinates, residue_numbers)):
        x, y, z = xyz
        atom_name, element = ("N", "N") if index in missing_ca else ("CA", "C")
        lines.append(
            f"ATOM  {serial:5d} {atom_name:^4s} {three[amino_acid]:>3s} {chain}{number:4d}{insertion[:1] or ' '}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}  "
        )
        serial += 1
    path.write_text("\n".join([*lines, "TER", "END", ""]))
    return path


def _make_pair(tmp_path: Path, *, source_sequence: str = "CDEFG", reference_sequence: str = "ACDEFGH",
               source_coordinates: np.ndarray | None = None, reference_coordinates: np.ndarray | None = None,
               source_missing_ca: set[int] | None = None, source_numbers=None):
    source_coordinates = np.asarray(source_coordinates if source_coordinates is not None else _BASE[:len(source_sequence)])
    if reference_coordinates is None:
        reference_coordinates = np.asarray([[20.0, 1.0, -2.0], *(_BASE @ _ROTATION + _TRANSLATION), [30.0, 2.0, 4.0]])
    source = _pdb(tmp_path / "source.pdb", "A", source_sequence, source_coordinates.tolist(),
                  residue_numbers=source_numbers, missing_ca=source_missing_ca)
    reference_numbers = [(100, ""), (101, ""), (102, "A"), (103, ""), (110, ""), (111, ""), (112, "")]
    reference = _pdb(tmp_path / "reference.pdb", "R", reference_sequence,
                     np.asarray(reference_coordinates).tolist(), residue_numbers=reference_numbers)
    mapping = component_reference_mapping(source, "A", reference, "R")
    return source, reference, mapping


def _identity(path: Path, *, molecule: str, structure: str, group: str, model: str, chain: str):
    return {
        "molecule_id": molecule, "structure_id": structure, "group_id": group,
        "structure_digest": hashlib.sha256(path.read_bytes()).hexdigest(),
        "model_id": model, "chain_id": chain,
    }


def _fit(source: Path, reference: Path, mapping: dict, *, source_chain="A", reference_chain="R",
         source_model="source-model", reference_model="reference-model", source_identity=None,
         reference_identity=None):
    source_identity = source_identity or _identity(source, molecule="source-record", structure="source-structure",
                                                   group="source-group", model=source_model, chain=source_chain)
    reference_identity = reference_identity or _identity(reference, molecule="reference-record", structure="reference-structure",
                                                         group="reference-group", model=reference_model, chain=reference_chain)
    return calculate_reference_placement(
        source, reference, source_chain=source_chain, reference_chain=reference_chain,
        source_model_id=source_model, reference_model_id=reference_model,
        source_digest=source_identity["structure_digest"], reference_digest=reference_identity["structure_digest"],
        mapping=mapping, source_identity=source_identity, reference_identity=reference_identity,
    )


def test_unique_truncated_construct_fits_known_transform_with_author_insertion_codes(tmp_path):
    source_numbers = [(20, ""), (21, "A"), (22, ""), (25, ""), (26, "")]
    source, reference, mapping = _make_pair(tmp_path, source_numbers=source_numbers)
    original = source.read_bytes()

    result = _fit(source, reference, mapping)

    assert mapping["mapping_status"] == "validated"
    assert result["aligned_pair_count"] == 5
    assert result["source_coverage"] == 1.0
    assert result["reference_coverage"] == pytest.approx(5 / 7)
    assert result["rmsd_angstrom"] == pytest.approx(0.0, abs=1e-5)
    matrix = np.asarray(result["matrix_4x4"])
    assert matrix[:3, :3] == pytest.approx(_ROTATION.T, abs=1e-5)
    assert matrix[:3, 3] == pytest.approx(_TRANSLATION, abs=1e-5)
    assert result["transformation_direction"].startswith("selected complex/source coordinates ->")
    assert any(item["source"]["insertion"] == "A" for item in result["mapped_residue_pairs"])
    assert source.read_bytes() == original


def test_ambiguous_or_nonidentical_mapping_is_not_fitted(tmp_path):
    coordinates = [[0, 0, 0], [2, 0, 0], [0, 2, 0]]
    source = _pdb(tmp_path / "ambiguous-source.pdb", "A", "ACD", coordinates)
    reference = _pdb(tmp_path / "ambiguous-reference.pdb", "R", "ACDACD", coordinates * 2)
    mapping = component_reference_mapping(source, "A", reference, "R")

    assert mapping["mapping_status"] == "ambiguous"
    with pytest.raises(ValueError, match="uniquely validated"):
        _fit(source, reference, mapping)


def test_internal_reference_insertion_is_not_treated_as_an_exact_truncated_subsequence(tmp_path):
    source_coords = [[0, 0, 0], [3, 0, 0], [3, 2, 0], [3, 2, 3]]
    reference_coords = [[0, 0, 0], [3, 0, 0], [1, 1, 1], [3, 2, 0], [3, 2, 3]]
    source = _pdb(tmp_path / "source-without-insertion.pdb", "A", "ACDE", source_coords)
    reference = _pdb(tmp_path / "reference-with-insertion.pdb", "R", "ACGDE", reference_coords)
    mapping = component_reference_mapping(source, "A", reference, "R")

    assert mapping["mapping_status"] == "partial"
    assert "not contiguous" in mapping["mapping_reason"]
    with pytest.raises(ValueError, match="uniquely validated"):
        _fit(source, reference, mapping)


def test_conformational_difference_is_reported_without_arbitrary_rmsd_cutoff(tmp_path):
    reference_coordinates = np.asarray([[20.0, 1.0, -2.0], *(_BASE @ _ROTATION + _TRANSLATION), [30.0, 2.0, 4.0]])
    reference_coordinates[3] += np.asarray([0.0, 0.0, 2.0])
    source, reference, mapping = _make_pair(tmp_path, reference_coordinates=reference_coordinates)

    result = _fit(source, reference, mapping)

    assert result["rmsd_angstrom"] > 0.1
    assert result["parameters"]["maximum_rmsd_cutoff"] is None
    assert any("no universal acceptance cutoff" in item for item in result["limitations"])


def test_insufficient_or_collinear_ca_coordinates_are_rejected(tmp_path):
    source, reference, mapping = _make_pair(tmp_path, source_sequence="CD", reference_sequence="CD",
                                            source_coordinates=np.asarray([[0, 0, 0], [2, 0, 0]]),
                                            reference_coordinates=np.asarray([[0, 0, 0], [2, 0, 0]]))
    with pytest.raises(ValueError, match="At least three"):
        _fit(source, reference, mapping)

    line = np.asarray([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]])
    source, reference, mapping = _make_pair(tmp_path, source_coordinates=line,
                                            reference_coordinates=np.asarray([[20, 1, -2], *(line @ _ROTATION + _TRANSLATION), [30, 2, 4]]))
    with pytest.raises(ValueError, match="collinear or geometrically degenerate"):
        _fit(source, reference, mapping)


def test_missing_ca_wrong_model_or_chain_and_stale_digest_are_rejected(tmp_path):
    source, reference, mapping = _make_pair(tmp_path, source_missing_ca={2})
    with pytest.raises(ValueError, match="no C-alpha atom"):
        _fit(source, reference, mapping)

    source, reference, mapping = _make_pair(tmp_path)
    with pytest.raises(ValueError, match="selected source model"):
        _fit(source, reference, mapping, source_identity={**_identity(source, molecule="m", structure="s", group="g", model="other", chain="A")})
    with pytest.raises(ValueError, match="selected source chain"):
        _fit(source, reference, mapping, source_chain="B")
    with pytest.raises(ValueError, match="selected reference chain"):
        _fit(source, reference, mapping, reference_chain="X")
    damaged_mapping = {**mapping, "residue_mapping": [dict(item) for item in mapping["residue_mapping"]]}
    damaged_mapping["residue_mapping"][0] = {
        **damaged_mapping["residue_mapping"][0],
        "source": {**damaged_mapping["residue_mapping"][0]["source"], "sequence_index": 3},
    }
    with pytest.raises(ValueError, match="cover each source residue exactly once"):
        _fit(source, reference, damaged_mapping)
    stale_identity = _identity(source, molecule="source-record", structure="source-structure",
                                group="source-group", model="source-model", chain="A")
    stale_identity["structure_digest"] = "0" * 64
    with pytest.raises(ValueError, match="source structure digest"):
        calculate_reference_placement(
            source, reference, source_chain="A", reference_chain="R",
            source_model_id="source-model", reference_model_id="reference-model",
            source_digest=hashlib.sha256(source.read_bytes()).hexdigest(),
            reference_digest=hashlib.sha256(reference.read_bytes()).hexdigest(),
            mapping=mapping, source_identity=stale_identity,
            reference_identity=_identity(reference, molecule="reference-record", structure="reference-structure",
                                         group="reference-group", model="reference-model", chain="R"),
        )
