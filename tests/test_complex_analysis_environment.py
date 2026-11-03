from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def _write_atoms(path: Path, *, source: bool, transform: tuple[np.ndarray, np.ndarray] | None = None,
                 include_context: bool = True) -> Path:
    source_positions = [[0, 0, 0], [3, 0, 0], [0, 4, 0]]
    matrix, translation = transform if transform is not None else (np.eye(3), np.zeros(3))
    inverse = matrix.T

    def coord(xyz):
        xyz = np.asarray(xyz, dtype=float)
        return inverse @ (xyz - translation) if source and transform is not None else xyz

    lines = ["HEADER    CONTROLLED REFERENCE ENVIRONMENT FIXTURE", "EXPDTA    X-RAY DIFFRACTION"]
    if not source:
        lines.append("SEQRES   1 A    7  ALA CYS ASP GLU ASN GLY SER")
    serial = 1

    def atom(record, name, residue, chain, number, xyz, element="C"):
        nonlocal serial
        x, y, z = coord(xyz)
        lines.append(
            f"{record:<6s}{serial:5d} {name:>4s} {residue:>3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{30.0:6.2f}          {element:>2s}  "
        )
        serial += 1

    for index, (residue, xyz) in enumerate(zip(("ALA", "CYS", "ASP"), source_positions), 1):
        atom("ATOM", "CA", residue, "A", index if source else index + 9, xyz)
    if source:
        atom("ATOM", "CA", "GLY", "B", 1, [30, 30, 30])
    elif include_context:
        atom("ATOM", "CA", "GLU", "A", 13, [0, 0, 2])
        atom("ATOM", "CA", "GLY", "B", 1, [5, 0, 0])
        atom("HETATM", "C1", "NAG", "A", 501, [0, 1, 0])
    else:
        # Reference contains only residues covered by the placement map.
        pass
    lines.extend(["TER", "END", ""])
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def _identity(path: Path, *, prefix: str, chain: str) -> dict:
    from mn_molecular_design_studio.evaluations import file_digest
    return {"molecule_id": f"{prefix}-molecule", "structure_id": f"{prefix}-structure",
            "group_id": f"{prefix}-group", "structure_digest": file_digest(path),
            "model_id": f"{prefix}-model", "chain_id": chain}


def _placement(source_identity, reference_identity, matrix=None):
    from mn_molecular_design_studio.complex_analysis_placement import PLACEMENT_VERSION
    return {
        "method_version": PLACEMENT_VERSION,
        "mapping_method": "unique_exact_sequence_subsequence",
        "transformation_direction": "selected complex/source coordinates -> associated reference coordinates",
        "matrix_4x4": (np.eye(4) if matrix is None else matrix).tolist(),
        "source": source_identity, "reference": reference_identity,
        "reference_context_id": "association-1", "evaluation_id": "placement-eval-1", "run_id": "placement-run-1",
        "aligned_pair_count": 3, "source_coverage": 1.0, "reference_coverage": 0.75,
        "rmsd_angstrom": 0.2,
        "mapped_residue_pairs": [
            {"source": {"number": index, "insertion": ""},
             "reference": {"number": index + 9, "insertion": ""}}
            for index in range(1, 4)
        ],
    }


def _calculate(source_path, reference_path, source_identity, reference_identity, *, selected=(), assembly_id=None, matrix=None):
    from mn_molecular_design_studio.complex_analysis_environment import calculate_reference_environment
    return calculate_reference_environment(
        source_path, reference_path, placement=_placement(source_identity, reference_identity, matrix),
        source_identity=source_identity, reference_identity=reference_identity,
        reference_context_id="association-1", selected_reference_chains=list(selected), assembly_id=assembly_id,
    )


def test_environment_reports_omitted_residues_glycan_coordinates_and_sequons_separately(tmp_path):
    source = _write_atoms(tmp_path / "source.pdb", source=True)
    text = source.read_text()
    unsupported_atom = (
        f"HETATM{9000:5d} {'XE':>4s} {'XEN':>3s} Z{999:4d}    "
        f"{100.0:8.3f}{100.0:8.3f}{100.0:8.3f}{1.0:6.2f}{30.0:6.2f}          {'XE':>2s}  "
    )
    ligand_atom = (
        f"HETATM{9001:5d} {'C1':>4s} {'LIG':>3s} A{901:4d}    "
        f"{5.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{30.0:6.2f}          {'C':>2s}  "
    )
    source.write_text(text.replace("END\n", unsupported_atom + "\n" + ligand_atom + "\nEND\n"))
    reference = _write_atoms(tmp_path / "reference.pdb", source=False)
    source_identity, reference_identity = _identity(source, prefix="source", chain="A"), _identity(reference, prefix="reference", chain="A")
    result = _calculate(source, reference, source_identity, reference_identity, selected=["B"])

    evidence = result["residue_evidence"]
    assert {row["context_kind"] for row in evidence} == {
        "omitted_reference_region", "selected_reference_component", "explicitly_modeled_sugar_residue",
    }
    omitted = next(row for row in evidence if row["context_kind"] == "omitted_reference_region")
    assert omitted["author_number"] == 13
    assert omitted["closest_distance_angstrom"] == pytest.approx(2.0)
    assert omitted["atom_pairs_with_vdw_overlap_above_epsilon"] == 1
    assert omitted["closest_atom_pair"]["source"]["chain_id"] == "A"
    assert omitted["closest_atom_pair"]["context"]["author_residue"] == "13"
    assert omitted["source_component_geometry"]["atom_pairs_within_proximity_cutoff"] > 0
    assert omitted["other_source_complex_geometry"]["atom_pairs_within_proximity_cutoff"] > 0
    assert omitted["other_source_complex_geometry"]["closest_atom_pair"]["source_role"] == "other_source_complex_component"
    sugar = next(row for row in evidence if row["context_kind"] == "explicitly_modeled_sugar_residue")
    assert sugar["residue_name"] == "NAG"
    partner = next(row for row in evidence if row["context_kind"] == "selected_reference_component")
    assert partner["other_source_complex_geometry"]["atom_pairs_with_vdw_overlap_above_epsilon"] > 0
    assert partner["other_source_complex_geometry"]["closest_atom_pair"]["source_role"] == "other_source_complex_component"
    assert result["reference_coordinate_origin"]["classification"] == "experimental"
    sequon = result["sequence_only_potential_n_linked_sequons"][0]
    assert sequon["sequence_motif"] == "NGS"
    assert sequon["sequence_only"] is True
    assert sequon["occupancy_or_glycan_geometry_measured"] is False
    assert result["measurements"]["unsupported_heavy_atom_count"] == 1
    assert result["parameters"]["unsupported_heavy_atom_elements"] == {"XE": 1}
    assert any("unsupported configured radii" in item for item in result["limitations"])
    assert any("Sequence-only N-X-S/T motifs are not measured glycans" in item for item in result["limitations"])


def test_environment_is_rigid_transform_invariant_and_keeps_identity_provenance(tmp_path):
    reference = _write_atoms(tmp_path / "reference.pdb", source=False)
    identity_source = _write_atoms(tmp_path / "source-identity.pdb", source=True)
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.asarray([12.0, -4.0, 7.0])
    transformed_source = _write_atoms(tmp_path / "source-transformed.pdb", source=True, transform=(rotation, translation))
    reference_identity = _identity(reference, prefix="reference", chain="A")
    identity = _identity(identity_source, prefix="source", chain="A")
    transformed_identity = _identity(transformed_source, prefix="source", chain="A")
    identity_matrix = np.eye(4)
    transformed_matrix = np.eye(4)
    transformed_matrix[:3, :3] = rotation
    transformed_matrix[:3, 3] = translation

    first = _calculate(identity_source, reference, identity, reference_identity, selected=["B"])
    second = _calculate(transformed_source, reference, transformed_identity, reference_identity,
                        selected=["B"], matrix=transformed_matrix)
    assert first["measurements"] == second["measurements"]
    assert first["source"] == identity
    assert second["source"] == transformed_identity
    assert first["placement"]["mapping_method"] == "unique_exact_sequence_subsequence"
    assert first["placement"]["transformation_direction"].startswith("selected complex/source")


def test_environment_reports_confirmed_assembly_membership_and_rejects_nonidentity_operator(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_environment as environment
    source = _write_atoms(tmp_path / "source.pdb", source=True)
    reference = _write_atoms(tmp_path / "reference.pdb", source=False)
    source_identity, reference_identity = _identity(source, prefix="source", chain="A"), _identity(reference, prefix="reference", chain="A")
    monkeypatch.setattr(environment, "_cif_assembly_membership", lambda *_args: {
        "status": "confirmed_identity_operator_member", "assembly_id": "1"})
    result = _calculate(source, reference, source_identity, reference_identity, selected=["B"], assembly_id="1")
    partner = next(row for row in result["residue_evidence"] if row["chain_id"] == "B")
    assert partner["context_kind"] == "confirmed_assembly_component"
    assert partner["assembly_membership"]["status"] == "confirmed_identity_operator_member"

    monkeypatch.setattr(environment, "_cif_assembly_membership", lambda *_args: {
        "status": "non_identity_assembly_operator_not_supported", "assembly_id": "1"})
    with pytest.raises(ValueError, match="non-identity operator"):
        _calculate(source, reference, source_identity, reference_identity, selected=["B"], assembly_id="1")


def test_assembly_membership_parsers_require_identity_operator_evidence(tmp_path):
    from mn_molecular_design_studio.complex_analysis_environment import (
        _cif_assembly_membership, _pdb_assembly_membership,
    )
    mmcif = Path(__file__).parents[1] / "mn_molecular_design_studio" / "assets" / "generators" / "pxdesign" / "5o45.cif"
    assert _cif_assembly_membership(mmcif, "1", "A")["status"] == "confirmed_identity_operator_member"

    pdb = tmp_path / "assembly.pdb"
    pdb.write_text("\n".join([
        "REMARK 350 BIOMOLECULE: 1",
        "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A",
        "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000",
        "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000",
        "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000",
    ]) + "\n")
    assert _pdb_assembly_membership(pdb, "1", "A")["status"] == "confirmed_identity_operator_member"
    pdb.write_text(pdb.read_text().replace("1.000000  0.000000  0.000000", "2.000000  0.000000  0.000000", 1))
    assert _pdb_assembly_membership(pdb, "1", "A")["status"] == "non_identity_assembly_operator_not_supported"


def test_environment_fails_closed_for_missing_context_wrong_identity_and_stale_digest(tmp_path):
    source = _write_atoms(tmp_path / "source.pdb", source=True)
    reference = _write_atoms(tmp_path / "mapped-only-reference.pdb", source=False, include_context=False)
    source_identity, reference_identity = _identity(source, prefix="source", chain="A"), _identity(reference, prefix="reference", chain="A")
    with pytest.raises(ValueError, match="No resolved omitted-reference"):
        _calculate(source, reference, source_identity, reference_identity)
    with pytest.raises(ValueError, match="already assessed for omitted residues"):
        _calculate(source, reference, source_identity, reference_identity, selected=["A"])

    with pytest.raises(ValueError, match="identities do not exactly match"):
        bad_placement = _placement(source_identity, reference_identity)
        bad_placement["reference"] = {**reference_identity, "structure_id": "another-reference"}
        from mn_molecular_design_studio.complex_analysis_environment import calculate_reference_environment
        calculate_reference_environment(source, reference, placement=bad_placement,
                                        source_identity=source_identity, reference_identity=reference_identity,
                                        reference_context_id="association-1")

    missing_method = _placement(source_identity, reference_identity)
    missing_method.pop("mapping_method")
    from mn_molecular_design_studio.complex_analysis_environment import calculate_reference_environment
    with pytest.raises(ValueError, match="does not record its residue-mapping method"):
        calculate_reference_environment(source, reference, placement=missing_method,
                                        source_identity=source_identity, reference_identity=reference_identity,
                                        reference_context_id="association-1")

    duplicate_map = _placement(source_identity, reference_identity)
    duplicate_map["mapped_residue_pairs"][1]["reference"] = duplicate_map["mapped_residue_pairs"][0]["reference"]
    with pytest.raises(ValueError, match="ambiguous because residue keys are repeated"):
        calculate_reference_environment(source, reference, placement=duplicate_map,
                                        source_identity=source_identity, reference_identity=reference_identity,
                                        reference_context_id="association-1")

    source.write_text(source.read_text() + "REMARK stale digest\n")
    with pytest.raises(ValueError, match="source structure digest is stale"):
        _calculate(source, reference, source_identity, reference_identity)


def test_predicted_coordinate_origin_is_reported_from_explicit_file_metadata(tmp_path):
    source = _write_atoms(tmp_path / "source.pdb", source=True)
    reference = _write_atoms(tmp_path / "reference-predicted.pdb", source=False)
    reference.write_text(reference.read_text().replace("EXPDTA    X-RAY DIFFRACTION", "EXPDTA    THEORETICAL MODEL"))
    source_identity, reference_identity = _identity(source, prefix="source", chain="A"), _identity(reference, prefix="reference", chain="A")
    result = _calculate(source, reference, source_identity, reference_identity)
    assert result["reference_coordinate_origin"] == {
        "classification": "predicted", "file_reported_method": ["THEORETICAL MODEL"],
    }
