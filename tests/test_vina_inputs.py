from pathlib import Path

import pytest

from mn_molecular_design_studio import ligand_compounds, vina_inputs


def _atom(serial, residue, chain, number, x, *, record="ATOM", name="CA"):
    return f"{record:<6}{serial:5d} {name:^4} {residue:>3} {chain}{number:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C\n"


def _atom_xyz(serial, residue, chain, number, xyz, *, record="ATOM", name="CA"):
    x, y, z = xyz
    return f"{record:<6}{serial:5d} {name:^4} {residue:>3} {chain}{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"


def test_probe_contacts_map_to_original_author_residues_and_frame():
    target = _atom(1, "ALA", "Q", 225, 100) + _atom(2, "GLY", "Q", 230, 120)
    predicted = (_atom(1, "ALA", "A", 1, 0) + _atom(2, "GLY", "A", 2, 20)
                 + _atom(3, "LIG", "B", 1, 21, record="HETATM", name="C1"))
    assert vina_inputs.predicted_contact_residues(target, predicted, 3) == {("Q", 230)}
    assert vina_inputs.predicted_contact_residues(target, predicted, 1) == {("Q", 230)}


def test_probe_contacts_map_selected_target_chain_to_boltz_input_order():
    target = (
        _atom_xyz(1, "ALA", "A", 10, (50, 0, 0))
        + _atom_xyz(2, "ALA", "B", 30, (100, 0, 0))
        + _atom_xyz(3, "GLY", "B", 31, (100, 10, 0))
    )
    # Boltz assigns the one selected PDB chain to A, regardless of its original
    # chain ID. An extra prediction-context chain follows it as B.
    predicted = (
        _atom_xyz(1, "ALA", "A", 1, (0, 0, 0))
        + _atom_xyz(2, "GLY", "A", 2, (0, 10, 0))
        + _atom_xyz(3, "ALA", "B", 1, (50, 0, 0))
        + _atom_xyz(4, "LIG", "C", 1, (1, 10, 0), record="HETATM", name="C1")
    )

    assert vina_inputs.predicted_contact_residues(
        target, predicted, 3, selected_chains=["B"],
    ) == {("B", 31)}


def test_predicted_probe_overlay_aligns_to_selected_original_chain():
    target = (
        _atom_xyz(1, "ALA", "A", 10, (50, 0, 0))
        + _atom_xyz(2, "ALA", "B", 30, (100, 0, 0))
        + _atom_xyz(3, "GLY", "B", 31, (100, 10, 0))
        + _atom_xyz(4, "SER", "B", 32, (100, 0, 10))
    )
    predicted = (
        _atom_xyz(1, "ALA", "A", 1, (0, 0, 0))
        + _atom_xyz(2, "GLY", "A", 2, (0, 10, 0))
        + _atom_xyz(3, "SER", "A", 3, (0, 0, 10))
        + _atom_xyz(4, "LIG", "B", 1, (1, 1, 1), record="HETATM", name="C1")
    )

    overlay = vina_inputs.predicted_ligand_overlay_pdb(
        target, predicted, selected_chains=["B"],
    )
    ligand_line = next(line for line in overlay.splitlines() if line.startswith("HETATM"))
    coordinates = tuple(float(ligand_line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
    assert coordinates == pytest.approx((101, 1, 1))
    assert all(line.startswith(("HETATM", "END")) for line in overlay.splitlines())


def test_residue_ranges_validate_missing_author_residues():
    assert vina_inputs.parse_residue_ranges("225, 230 - 232", {225, 230, 231, 232}) == {225, 230, 231, 232}
    assert vina_inputs.residue_range_text({225, 230, 231, 232}) == "225, 230-232"
    with pytest.raises(ValueError, match="not present"):
        vina_inputs.parse_residue_ranges("225-230", {225, 230})


def test_example_is_offline_portable_and_owned_by_current_user(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    owner = "a" * 32
    target, compound = vina_inputs.stage_example(owner, project_id="project-1")
    payload = ligand_compounds.load_compound_set(Path(compound["path"]))
    assert "STI" in target.read_text()
    assert compound["compound_count"] == 1
    assert payload["owner_id"] == owner
    assert payload["project_id"] == "project-1"
    assert payload["compounds"][0]["compound_id"] == "imatinib-CID5291"
