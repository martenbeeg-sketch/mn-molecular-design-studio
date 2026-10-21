from mn_molecular_design_studio.docking_results import display_docking_rows, pose_visualizations


def test_docking_score_display_resolves_id_name_and_smiles_from_compound_set():
    rows = display_docking_rows(
        [{"compound_id": "70701262", "best_score_kcal_mol": -8.4, "status": "completed"}],
        [{"compound_id": "70701262", "compound_name": "Bepristat_1a", "smiles": "CCN"}],
    )

    assert rows == [{
        "Compound ID": "70701262",
        "Compound name": "Bepristat_1a",
        "SMILES": "CCN",
        "best_score_kcal_mol": -8.4,
        "status": "completed",
    }]


def test_docking_score_identity_defaults_name_to_id_for_legacy_runs():
    rows = display_docking_rows([{"compound_id": "ML359", "score": -7.1}])

    assert rows[0]["Compound ID"] == "ML359"
    assert rows[0]["Compound name"] == "ML359"
    assert rows[0]["SMILES"] == "—"


def _atom(record, serial, atom, residue, chain, residue_number, x, y, z):
    return (
        f"{record:<6}{serial:5d} {atom:>4} {residue:>3} {chain}{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {atom[0]:>2}"
    )


def test_rosettaligand_pose_preview_splits_protein_and_bound_ligand():
    protein_records = [
        record
        for serial in range(1, 7)
        for record in (
            _atom("ATOM", serial * 2 - 1, "CA", "ALA", "A", serial, serial * 10.0, 0.0, 0.0),
            _atom("ATOM", serial * 2, "CB", "ALA", "A", serial, serial * 10.0 + 1.0, 0.0, 0.0),
        )
    ]
    ligand_records = [
        _atom("HETATM", 101, "C1", "LIG", "Z", 1, 11.0, 1.0, 0.0),
        _atom("HETATM", 102, "O1", "LIG", "Z", 1, 12.0, 1.0, 0.0),
    ]
    pose = "\n".join([*protein_records, *ligand_records, "END", ""])

    structures = pose_visualizations(pose, "docked_pose.pdb")

    assert len(structures) == 2
    protein, ligand = structures
    assert protein.representation_type == "cartoon"
    assert protein.color_params == {"value": "0x318c70"}
    assert protein.chains[0].chain_id == "A"
    assert protein.chains[0].residues == [1]
    assert protein.chains[0].representation_type == "ball-and-stick"
    assert "ATOM  " in protein.pdb
    assert "HETATM" not in protein.pdb
    assert ligand.representation_type == "ball-and-stick"
    assert ligand.color_params == {"value": "0xe08b3e"}
    assert "HETATM" in ligand.pdb
    assert "ATOM  " not in ligand.pdb


def test_pdbqt_pose_preview_pairs_ligand_with_prepared_receptor():
    receptor = "\n".join(
        [
            _atom("ATOM", serial, "CA", "ALA", "A", serial, float(serial), 0.0, 0.0)
            for serial in range(1, 7)
        ]
        + ["END", ""]
    )
    ligand = "ATOM      1  C1  LIG Z   1       2.000   1.000   0.000  0.000 C\nEND\n"

    structures = pose_visualizations(ligand, "vina_pose.pdbqt", receptor, "receptor.pdb")

    assert len(structures) == 2
    assert "ALA" in structures[0].pdb
    assert "LIG" not in structures[0].pdb
    assert "LIG" in structures[1].pdb
    assert structures[1].pdb.rstrip().endswith("END")
