from __future__ import annotations

import pytest

from mn_molecular_design_studio.molecule_imports import parse_structure


def test_parse_generated_mmcif_without_optional_occupancy():
    cif = """data_model
#
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.B_iso_or_equiv
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.id
ATOM N N . ALA A 1 . 1 ALA A N 20.0 1 0.0 0.0 0.0 1 1
ATOM C CA . ALA A 1 . 1 ALA A CA 20.0 1 1.4 0.0 0.0 1 2
ATOM C C . ALA A 1 . 1 ALA A C 20.0 1 2.8 0.0 0.0 1 3
ATOM O O . ALA A 1 . 1 ALA A O 20.0 1 3.3 1.2 0.0 1 4
#
"""

    record = parse_structure(cif.encode(), "generated.cif")[0]

    assert record["components"] == [{"type": "protein", "chain": "A", "role": "", "sequence": "A"}]
    assert record["files"][0]["name"] == "generated.cif"


def test_parse_generated_mmcif_without_optional_bfactor_and_occupancy_preserves_bytes():
    cif = b"""data_model
#
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.id
ATOM N N . ALA A 1 . 1 ALA A N 1 0.0 0.0 0.0 1 1
ATOM C CA . ALA A 1 . 1 ALA A CA 1 1.4 0.0 0.0 1 2
ATOM C C . ALA A 1 . 1 ALA A C 1 2.8 0.0 0.0 1 3
ATOM O O . ALA A 1 . 1 ALA A O 1 3.3 1.2 0.0 1 4
#
"""

    record = parse_structure(cif, "generated.cif")[0]

    assert record["components"] == [{"type": "protein", "chain": "A", "role": "", "sequence": "A"}]
    assert record["files"][0]["content"] == cif


def test_parse_structure_still_rejects_malformed_mmcif_coordinates():
    cif = b"""data_model
#
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.B_iso_or_equiv
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.id
ATOM N N . ALA A 1 . 1 ALA A N 20.0 1 not-a-coordinate 0.0 0.0 1 1
#
"""

    with pytest.raises(ValueError, match="Could not read generated.cif"):
        parse_structure(cif, "generated.cif")
