import copy
import xml.etree.ElementTree as ET

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.md_atom_interactions import (
    build_interaction_diagram, interaction_residues, ligand_drawing, render_interaction_diagram,
)
from test_md_results_ui import FakeStreamlit
from test_database_browser import USER, database
from test_md_browser import _add_md, _cached, _reader


def _structural():
    return {"analysis_window": {"analyzed_frames": 100, "trajectory_stride": 3}, "contacts": {"residues": [
        {"residue": "ARG262 · chain A", "contact_occupancy": 1, "hydrogen_bond_occupancy": 0.75,
            "salt_bridge_occupancy": 0.25},
        {"residue": "LEU276 · chain A", "contact_occupancy": 0.8, "hydrophobic_occupancy": 0.5},
        {"residue": "ARG262 · chain B", "contact_occupancy": 0.3},
        {"residue": "GLY278 · chain A", "contact_occupancy": 0},
    ]}}


def _ligand_pdb(smiles="CCO", *, reorder=False, duplicate_names=False):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(molecule, randomSeed=71) == 0
    for index, atom in enumerate(molecule.GetAtoms()):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(("C11" if duplicate_names else f"{atom.GetSymbol()}{index + 11}").rjust(4))
        info.SetResidueName("LIG")
        info.SetResidueNumber(1)
        info.SetChainId("L")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    if reorder:
        molecule = Chem.RenumberAtoms(molecule, list(reversed(range(molecule.GetNumAtoms()))))
    return Chem.MolToPDBBlock(molecule)


@pytest.mark.parametrize("reorder", [False, True])
def test_lines_end_at_named_source_atoms_not_smiles_indices_or_panel_boundary(reorder):
    pdb_text = _ligand_pdb(reorder=reorder)
    structural = {"contacts": {"residues": [
        {"residue": "SER42 · chain A", "contact_occupancy": 0.8, "hydrogen_bond_occupancy": 0.7, "top_ligand_atom": "O13"},
        {"residue": "PHE42 · chain B", "contact_occupancy": 0.9, "hydrophobic_occupancy": 0.8, "top_ligand_atom": "C11"}]}}
    diagram = build_interaction_diagram(structural, ligand_smiles="CCO", ligand_pdb=pdb_text)
    assert diagram.mapped_residues == 2
    drawing = ligand_drawing("CCO", ligand_pdb=pdb_text, label_atoms={"C11", "O13"})
    paths = [node for node in ET.fromstring(diagram.svg).iter() if node.get("data-ligand-atom")]
    assert {path.get("data-ligand-atom") for path in paths} == {"C11", "O13"}
    for path in paths:
        atom_x, atom_y = drawing.atoms[path.get("data-ligand-atom")]
        assert path.get("d").startswith(f"M {atom_x:.2f},{atom_y:.2f} Q ")
    assert "recorded representative" in diagram.svg
    assert diagram.atom_rows[0]["Recorded representative ligand atom"] == "C11"


@pytest.mark.parametrize("pdb_text,smiles", [("", "CCO"), (_ligand_pdb(duplicate_names=True), "CCO"),
    (_ligand_pdb(), "CCN"), (_ligand_pdb(), "CC"), (_ligand_pdb("C[C@H](O)F"), "C[C@@H](O)F")])
def test_absent_ambiguous_or_mismatched_atom_identity_never_gets_an_arbitrary_endpoint(pdb_text, smiles):
    structural = {"contacts": {"residues": [{"residue": "SER42", "contact_occupancy": 0.8, "top_ligand_atom": "C11"}]}}
    diagram = build_interaction_diagram(structural, ligand_smiles=smiles, ligand_pdb=pdb_text)
    assert diagram.chemical_depiction and diagram.mapped_residues == 0
    assert not [node for node in ET.fromstring(diagram.svg).iter() if node.get("data-ligand-atom")]
    assert "not guessed" in diagram.atom_mapping_note


def test_missing_atom_names_are_not_inferred_from_network_edges_or_other_residues():
    structural = {"contacts": {"residues": [
        {"residue": "SER42 · chain A", "contact_occupancy": 0.8, "top_ligand_atom": "O13"},
        {"residue": "SER42 · chain B", "contact_occupancy": 0.8},
        {"residue": "SER43 · chain A", "contact_occupancy": 0.8, "top_ligand_atom": "X99"}]},
        "interaction_network": {"edges": [{"source": "C11", "target": "SER42 · chain B", "weight": 0.8}]}}
    diagram = build_interaction_diagram(structural, ligand_smiles="CCO", ligand_pdb=_ligand_pdb())
    assert diagram.shown == 3 and diagram.mapped_residues == 1
    assert diagram.atom_rows[1]["Atom endpoint"] == "unavailable / ambiguous"
    paths = [node for node in ET.fromstring(diagram.svg).iter() if node.get("data-ligand-atom")]
    assert len(paths) == 1 and paths[0].get("data-residue") == "SER42 · chain A"


def test_diagram_renders_linked_chemical_identity_and_recorded_types_without_mutating_data():
    structural = _structural()
    before = copy.deepcopy(structural)
    diagram = build_interaction_diagram(structural, ligand_smiles="CC(=O)O", ligand_label="Compound A",
        scope="Replica 2 · production run-two")
    root = ET.fromstring(diagram.svg)
    titles = [node.text for node in root.iter() if node.tag.endswith("title")]
    assert diagram.chemical_depiction
    assert diagram.shown == diagram.eligible == 3
    assert diagram.atom_rows[0]["Hydrogen bond residue occupancy (%)"] == 75
    assert diagram.atom_rows[0]["Salt bridge residue occupancy (%)"] == 25
    assert "ARG262 · chain B: Geometric contact 30.0%" in titles
    assert not any(title == "ARG262 · chain A: Geometric contact 100.0%" for title in titles)
    assert "Replica 2 · production run-two" in diagram.svg
    assert "trajectory_stride" in diagram.svg
    assert "not atom-pair occupancies or atom-specific bonds" in diagram.svg
    assert "Compound A" in diagram.svg
    assert structural == before


def test_threshold_types_and_limit_use_exact_occupancies_and_deterministic_order():
    structural = _structural()
    hydrogen_bonds = build_interaction_diagram(structural, interaction="hydrogen_bond", minimum_occupancy=0.75)
    assert hydrogen_bonds.shown == 1
    assert build_interaction_diagram(structural, interaction="hydrogen_bond", minimum_occupancy=0.751).shown == 0
    assert build_interaction_diagram(structural, interaction="water_bridge").shown == 0
    first = build_interaction_diagram(structural, residue_limit=1)
    structural["contacts"]["residues"].reverse()
    second = build_interaction_diagram(structural, residue_limit=1)
    assert first.svg == second.svg
    assert first.shown == 1 and first.eligible == 3


@pytest.mark.parametrize("value", [True, False, "not a number", None, float("nan"), float("inf"), -0.1, 1.01, 75])
def test_invalid_values_are_not_clamped_or_interpreted_as_percentages(value):
    structural = {"contacts": {"residues": [{"residue": "ALA1", "contact_occupancy": value}]}}
    assert build_interaction_diagram(structural, minimum_occupancy=0).shown == 0


def test_missing_contact_rows_use_network_nodes_but_never_invent_atom_or_hydrogen_bond_links():
    structural = {"interaction_network": {"nodes": [{"id": "ALA1 · chain A", "contact_occupancy": 0.4}],
        "edges": [{"source": "C9", "target": "ALA1 · chain A", "interaction": "hydrogen_bond", "weight": 0.9}]}}
    diagram = build_interaction_diagram(structural)
    assert diagram.shown == 1 and not diagram.chemical_depiction
    assert "Chemical depiction unavailable" in diagram.svg
    assert "ALA1 · chain A: Geometric contact 40.0%" in diagram.svg
    assert "ALA1 · chain A: Hydrogen bond" not in diagram.svg
    assert "C9" not in diagram.svg


def test_conflicting_duplicate_residue_rows_are_excluded_instead_of_arbitrarily_merged():
    structural = _structural()
    structural["contacts"]["residues"].append({"residue": "ARG262 · chain A", "contact_occupancy": 0.9})
    labels = [row["residue"] for row in interaction_residues(structural)]
    assert "ARG262 · chain A" not in labels
    assert "ARG262 · chain B" in labels


def test_svg_escapes_external_labels_and_handles_missing_or_invalid_chemical_identity():
    structural = {"contacts": {"residues": [{"residue": '<script>alert("bad")</script>', "contact_occupancy": 0.5}]}}
    for smiles in ("", "not-a-smiles", "C" * 4097):
        diagram = build_interaction_diagram(structural, ligand_smiles=smiles, ligand_label='<img onerror="bad">')
        assert not diagram.chemical_depiction
        assert "<script>" not in diagram.svg and "<img" not in diagram.svg
        root = ET.fromstring(diagram.svg)
        assert not any(node.tag.endswith("script") for node in root.iter())


@pytest.mark.parametrize("settings", [{"interaction": "pi_stacking"}, {"minimum_occupancy": 2}, {"residue_limit": 100}, {"residue_limit": True}])
def test_bad_display_settings_fail_explicitly(settings):
    with pytest.raises(ValueError):
        build_interaction_diagram(_structural(), **settings)


def test_empty_report_is_explicit_not_substituted_with_another_structure():
    streamlit = FakeStreamlit()
    render_interaction_diagram(streamlit, {})
    assert any(name == "info" and "not inferred" in args[0] for name, args, _ in streamlit.calls)
    assert not any(name == "image" for name, _, _ in streamlit.calls)


def test_native_geometric_contact_diagram_never_claims_classified_interactions():
    from mn_molecular_design_studio.md_results_ui import _render_md_analysis_summary

    artifacts = {"artifacts/md_analysis/residue_contacts.csv": "residue,contact_occupancy_fraction\nA:1:ALA,0.6\n"}
    streamlit = FakeStreamlit()
    _render_md_analysis_summary(streamlit, {"success": True, "source": {}, "parameters": {}, "selected_frame_count": 10},
        text_loader=lambda relative, description: artifacts.get(relative))
    image = next(args[0] for name, args, _ in streamlit.calls if name == "image")
    assert "A:1:ALA: Geometric contact 60.0%" in image
    assert "A:1:ALA: Hydrogen bond" not in image
    table = next(args[0] for name, args, _ in streamlit.calls if name == "dataframe" and "Geometric contact (%)" in str(args[0]))
    assert table[0]["Geometric contact (%)"] == 60
    assert "Hydrogen bond (%)" not in table[0]


def test_streamlit_controls_and_enlarged_dialog_render_without_recalculation():
    source = f'''import streamlit as st
from mn_molecular_design_studio.md_atom_interactions import render_interaction_diagram
render_interaction_diagram(st, {_structural()!r}, ligand_smiles="CC(=O)O", ligand_label="Compound A", key_prefix="test")
'''
    app = AppTest.from_string(source).run()
    assert not app.exception
    assert len(app.get("image")) == 1
    app.selectbox[0].set_value("water_bridge").run()
    assert not app.exception
    assert any("No recorded interactions match" in row.value for row in app.info)
    app.selectbox[0].set_value("hydrogen_bond").run()
    next(button for button in app.button if button.label == "Enlarge interaction diagram").click().run()
    assert not app.exception
    assert len(app.get("image")) == 2


def test_database_browser_diagram_follows_selected_replica_and_linked_smiles(database, monkeypatch):
    import mn_molecular_design_studio.md_atom_interactions as diagrams

    reader, index, outcome = _reader(database)
    _add_md(database, replica=2)
    for replica in (1, 2):
        _cached(reader, f"artifacts/legacy/study-one-{replica}/result.json", {"success": True, "ligand_pdb": _ligand_pdb(),
            "md_result": {"status": "completed", "analytics": {"structural_dynamics": {
                "contacts": {"residues": [{"residue": f"ARG{replica} · chain A", "contact_occupancy": 0.8,
                    "water_bridge_occupancy": replica * 0.03,
                    "top_ligand_atom": "O13" if replica == 1 else "C11"}]}}}}})
    with database.connection() as connection:
        connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.ligand_smiles','CCO') WHERE namespace='md.production'")
    drawn = []
    original = diagrams.build_interaction_diagram

    def capture(structural, **settings):
        diagram = original(structural, **settings)
        drawn.append(diagram)
        return diagram

    monkeypatch.setattr(diagrams, "build_interaction_diagram", capture)
    source = f'''from mn_molecular_design_studio.data_browser import ReadOnlyLibrary, saved_results, listed_results
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from mn_molecular_design_studio.md_browser_ui import render_md_study
database=ReadOnlyLibrary({str(database.root)!r})
user={USER!r}
index=build_browser_index(database,user,collection="library",include_related=True)
raw=subject_results(database,index,{{"id":"all","molecule_ids":set(index["molecules"]),"group_ids":set(index["groups"])}})
outcome=next(row for row in saved_results(index,raw) if row["outcome_kind"]=="Dynamics")
render_md_study(database,user,index,listed_results([outcome])[0],{{}})
'''
    app = AppTest.from_string(source).run()
    next(element for element in app.selectbox if element.label == "Analysis").set_value("Contacts").run()
    assert not app.exception
    assert drawn[-1].chemical_depiction
    assert drawn[-1].mapped_residues == 1
    assert drawn[-1].atom_rows[0]["Recorded representative ligand atom"] == "O13"
    assert "ARG1 · chain A" in drawn[-1].svg and "ARG2 · chain A" not in drawn[-1].svg
    occupancy_table = next(table.value for table in app.dataframe if "Water bridge (%)" in table.value.columns)
    assert list(occupancy_table["Residue"]) == ["ARG1 · chain A"]
    assert occupancy_table["Water bridge (%)"].iloc[0] == pytest.approx(3)
    next(element for element in app.selectbox if element.label == "Replica").set_value("2").run()
    assert not app.exception
    assert "ARG2 · chain A" in drawn[-1].svg and "ARG1 · chain A" not in drawn[-1].svg
    assert "Replica 2 · production study-one-production-2" in drawn[-1].svg
    assert drawn[-1].atom_rows[0]["Recorded representative ligand atom"] == "C11"
    occupancy_table = next(table.value for table in app.dataframe if "Water bridge (%)" in table.value.columns)
    assert list(occupancy_table["Residue"]) == ["ARG2 · chain A"]
    assert occupancy_table["Water bridge (%)"].iloc[0] == pytest.approx(6)
