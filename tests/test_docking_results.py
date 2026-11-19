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

    structures = pose_visualizations(pose, "docked_pose.pdb", include_pocket_sidechains=True)

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


def _cif_complex():
    import gemmi

    protein = [_atom("ATOM", serial, "CA", "ALA", "A", serial, float(serial), 0.0, 0.0) for serial in range(1, 7)]
    ligand = [_atom("HETATM", 101, "C1", "LIG", "L", 1, 2.0, 1.0, 0.0),
              _atom("HETATM", 102, "O1", "LIG", "L", 1, 3.4, 1.0, 0.0)]
    structure = gemmi.read_pdb_string("\n".join([*protein, *ligand, "END", ""]))
    structure[0]["L"][0].name = "LIG1"
    structure[0]["L"].name = "ligand_long_chain"
    document = structure.make_mmcif_document()
    bonds = document.sole_block().init_loop("_chem_comp_bond.", ["comp_id", "atom_id_1", "atom_id_2", "value_order"])
    bonds.add_row(["LIG1", "C1", "O1", "DOUB"])
    return document.as_string()


def test_cofolded_cif_retains_native_ligand_names_coordinates_and_bond_tables():
    text = _cif_complex()
    structures = pose_visualizations(text, "boltz_complex.cif")
    assert len(structures) == 1
    assert structures[0].pdb == text
    assert "LIG1" in structures[0].pdb and "ligand_long_chain" in structures[0].pdb
    assert "_chem_comp_bond.value_order" in structures[0].pdb
    assert structures[0].representation_type == "cartoon" and structures[0].chains is None


def test_full_pdb_complex_defaults_to_cartoon_with_all_ligands_and_connectivity():
    protein = [_atom("ATOM", serial, "CA", "ALA", "A", serial, float(serial), 0.0, 0.0) for serial in range(1, 7)]
    ligand = [_atom("HETATM", 101, "C1", "LIG", "L", 1, 2.0, 1.0, 0.0),
              _atom("HETATM", 102, "O1", "LIG", "L", 1, 3.4, 1.0, 0.0),
              _atom("HETATM", 103, "C1", "DRG", "X", 2, 4.0, 1.0, 0.0)]
    text = "\n".join([*protein, *ligand, "CONECT  101  102", "END", ""])
    structures = pose_visualizations(text, "complex.pdb")
    assert len(structures) == 1 and structures[0].pdb == text
    assert structures[0].representation_type == "cartoon" and structures[0].chains is None


def test_native_cif_receptor_can_be_shown_with_its_saved_ligand_pose():
    ligand = "ATOM      1  C1  LIG Z   1       2.000   1.000   0.000  0.000 C\nEND\n"
    receptor = _cif_complex()
    structures = pose_visualizations(ligand, "pose.pdbqt", receptor, "receptor.cif")
    assert len(structures) == 2
    assert structures[0].pdb == receptor and structures[0].representation_type == "cartoon"
    assert not structures[0].chains
    assert structures[1].representation_type == "ball-and-stick"


def test_ligand_only_cif_stays_native_when_combined_with_a_verified_receptor():
    import gemmi

    ligand = gemmi.read_pdb_string(_atom("HETATM", 101, "C1", "LIG", "L", 1, 2.0, 1.0, 0.0) + "\nEND\n")
    text = ligand.make_mmcif_document().as_string()
    receptor = "\n".join([*[_atom("ATOM", serial, "CA", "ALA", "A", serial, float(serial), 0.0, 0.0)
                            for serial in range(1, 7)], "END", ""])
    structures = pose_visualizations(text, "ligand.cif", receptor, "receptor.pdb")
    assert len(structures) == 2
    assert structures[0].representation_type == "cartoon"
    assert structures[1].pdb == text and structures[1].representation_type == "ball-and-stick"
