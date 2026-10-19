from pathlib import Path

import pytest

from mn_molecular_design_studio import ligand_compounds, vina_inputs


def _atom(serial, residue, chain, number, x, *, record="ATOM", name="CA"):
    return f"{record:<6}{serial:5d} {name:^4} {residue:>3} {chain}{number:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C\n"


def test_probe_contacts_map_to_original_author_residues_and_frame():
    target = _atom(1, "ALA", "Q", 225, 100) + _atom(2, "GLY", "Q", 230, 120)
    predicted = (_atom(1, "ALA", "A", 1, 0) + _atom(2, "GLY", "A", 2, 20)
                 + _atom(3, "LIG", "B", 1, 21, record="HETATM", name="C1"))
    assert vina_inputs.predicted_contact_residues(target, predicted, 3) == {("Q", 230)}
    assert vina_inputs.predicted_contact_residues(target, predicted, 1) == {("Q", 230)}


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
