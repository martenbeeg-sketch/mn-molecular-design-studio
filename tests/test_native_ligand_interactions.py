import numpy as np
import pytest

from mn_molecular_design_studio.engine.workflow_engine.tools.ligand_interactions import (
    Atom, Classifier, Occupancy, _matched_ligand_charges,
)


def test_all_interaction_classes_track_exact_pairs_and_unique_residues_per_frame():
    atoms = [
        Atom("CB", "C", "A:1:ALA", "ALA"),
        Atom("NZ", "N", "A:2:LYS", "LYS"),
        Atom("HZ1", "H", "A:2:LYS", "LYS"),
        Atom("OD1", "O", "A:3:ASP", "ASP"),
        Atom("O1", "O", "L:1:LIG", "LIG"),
        Atom("C1", "C", "L:1:LIG", "LIG"),
        Atom("O", "O", "W:1:HOH", "HOH"),
        Atom("H1", "H", "W:1:HOH", "HOH"),
        Atom("H2", "H", "W:1:HOH", "HOH"),
    ]
    # The water sits between ligand O1 and Asp OD1, while the Lys H
    # points directly at O1. The alternate frame moves the ligand away.
    xyz = np.array([
        [0.0, -3.0, 0.0], [-3.0, 0.0, 0.0], [-2.0, 0.0, 0.0],
        [5.6, 0.0, 0.0], [0.0, 0.0, 0.0], [1.25, -0.25, 0.0],
        [2.8, 0.0, 0.0], [1.8, 0.0, 0.0], [3.75, 0.25, 0.0],
    ])
    classifier = Classifier(atoms, [0, 1, 3], [4, 5], [6],
                            {1: [2], 6: [7, 8]}, {3, 4, 6}, {4: -1})
    occupancy = Occupancy(classifier)
    counts = occupancy.update(classifier.frame(xyz))
    assert counts == {"contact_residues": 3, "hydrophobic_residues": 1,
                      "hydrogen_bond_residues": 1, "water_bridge_residues": 1,
                      "salt_bridge_residues": 1}
    xyz[[4, 5]] += 100
    occupancy.update(classifier.frame(xyz))
    table = {row["residue"]: row for row in occupancy.report()["contacts"]["residues"]}
    assert table["A:2:LYS"]["contact_occupancy"] == 0.5
    assert table["A:3:ASP"]["water_bridge_occupancy"] == 0.5
    assert table["A:3:ASP"]["atom_endpoints"]["water_bridge"]["ligand_atom"] == "O1"
    hbond, = [row for row in occupancy.pair_rows() if row["interaction"] == "hydrogen_bond"]
    assert (hbond["protein_atom"], hbond["ligand_atom"], hbond["frames"], hbond["occupancy_percent"]) == ("NZ", "O1", 1, 50)


def test_pose_distance_evidence_is_distinct_and_missing_trajectory_chemistry_is_unavailable():
    atoms = [Atom("NZ", "N", "A:1:LYS", "LYS"), Atom("O1", "O", "L:1:LIG", "LIG")]
    coordinates = np.array([[0.0, 0.0, 0.0], [3.2, 0.0, 0.0]])
    pose = Classifier(atoms, [0], [1], [], {0: []}, {1}, None, pose=True)
    assert ("hydrogen_bond", 0, 1, "distance_only") in pose.frame(coordinates)
    trajectory = Classifier(atoms, [0], [1], [], {}, {1}, None)
    occupancy = Occupancy(trajectory)
    occupancy.update(trajectory.frame(coordinates))
    residue, = occupancy.report()["contacts"]["residues"]
    assert residue["hydrogen_bond_occupancy"] is None
    assert residue["salt_bridge_occupancy"] is None


def test_saved_ligand_charge_mapping_requires_an_exact_protonation_match():
    md = pytest.importorskip("mdtraj")
    topology = md.Topology()
    residue = topology.add_residue("LIG", topology.add_chain("L"))
    oxygen = topology.add_atom("O1", md.element.oxygen, residue)
    carbon = topology.add_atom("C1", md.element.carbon, residue)
    topology.add_bond(oxygen, carbon)
    for index in range(3):
        hydrogen = topology.add_atom(f"H{index + 1}", md.element.hydrogen, residue)
        topology.add_bond(carbon, hydrogen)
    assert _matched_ligand_charges(topology, [0, 1], "C[O-]") == {0: -1}
    assert _matched_ligand_charges(topology, [0, 1], "CO") is None
