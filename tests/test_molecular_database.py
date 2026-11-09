from __future__ import annotations

import io
import json
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from mn_molecular_design_studio import molecule_imports as imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase, save_legacy_groups

A = {"id": "a" * 32, "role": "user", "_cluster_project_ids": ["project-1"]}
B = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}
ADMIN = {"id": "c" * 32, "role": "admin"}


@pytest.fixture
def database(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    return MolecularDatabase()


def test_reimport_and_group_membership_preserve_identity_structures_and_scores(database):
    first = database.import_records(A, imports.parse_smiles("CCO ethanol"), name="Inputs")
    molecule_id = first["molecule_ids"][0]
    database.add_observation(first["group_id"], molecule_id, "vina", "run-1", {"energy": -7.4}, context={"target": "target-1"})
    second = database.import_records(A, imports.parse_smiles("OCC renamed"), name="Shortlist")
    assert second["molecule_ids"] == [molecule_id]
    result = database.list_molecules(A, group_id=second["group_id"])
    assert result["total"] == 1
    assert result["items"][0]["scores"] == {"vina.energy": -7.4}
    assert len(database.groups(A)) == 2
    assert database.metric_fields(A, group_id=second["group_id"]) == ["vina.energy"]
    assert database.list_molecules(A, group_id=second["group_id"], metric_filters=[{"field": "vina.energy", "operator": "<", "value": -7}])["total"] == 1


def test_collection_statistics_aggregate_selected_group_counts_and_enforce_access(database):
    imported = database.import_records(A, imports.parse_sequences(">binder\nACDEFG"), name="Design cohort")
    molecule_id = imported["molecule_ids"][0]
    group_id = imported["group_id"]
    database.attach_structure(
        group_id, molecule_id, {"name": "model.cif", "content": b"data_model\n"},
        context={"model_id": "model-1"},
    )
    database.add_observation(group_id, molecule_id, "biotite-interface-v1", "analysis-1", {"contacts": 12})
    from mn_molecular_design_studio.evaluations import install_schema, rebuild_missing

    with database.connection() as connection:
        install_schema(connection)
        rebuild_missing(connection)

    assert database.collection_statistics(A, [group_id]) == {
        group_id: {"molecule_count": 1, "structure_count": 1, "evaluated_design_count": 1}
    }

    private = database.import_records(B, imports.parse_sequences(">private\nGGGGGG"), name="Private cohort")
    with pytest.raises(PermissionError, match="not accessible"):
        database.collection_statistics(A, [private["group_id"]])


def test_selected_compound_snapshot_contains_exact_members_and_is_immutable(database):
    result = database.import_records(A, imports.parse_smiles("CCO selected\nCC excluded"), name="Original")
    snapshot = database.compound_snapshot(A, [result["molecule_ids"][0]], project_id="project-1")
    payload = json.loads(Path(snapshot["path"]).read_text())
    assert [row["compound_id"] for row in payload["compounds"]] == ["selected"]
    before = Path(snapshot["path"]).read_bytes()
    database.import_records(A, imports.parse_smiles("CCN added"), group_id=result["group_id"])
    assert Path(snapshot["path"]).read_bytes() == before
    assert "/users/" in snapshot["path"]


def test_adding_selection_respects_destination_group_schema(database):
    source = database.import_records(A, imports.parse_smiles("CCO ethanol"), name="Source")
    destination = database.create_group(A, "Constrained group")
    database.set_schema(A, destination, {"fields": [{"name": "affinity", "type": "number", "required": True}]})

    with pytest.raises(ValueError, match="Required property is missing: affinity"):
        database.add_to_group(A, destination, source["molecule_ids"])

    assert database.list_molecules(A, group_id=destination)["total"] == 0


def test_removing_selected_molecules_is_scoped_to_the_current_group(database):
    first = database.import_records(A, imports.parse_smiles("CCO ethanol"), name="First")
    molecule_id = first["molecule_ids"][0]
    second = database.import_records(A, imports.parse_smiles("OCC ethanol"), name="Second")
    database.add_observation(first["group_id"], molecule_id, "vina", "run-1", {"energy": -7.0})

    assert database.remove_from_group(A, first["group_id"], [molecule_id]) == 1
    assert database.list_molecules(A, group_id=first["group_id"])["total"] == 0
    second_record = database.list_molecules(A, group_id=second["group_id"])["items"][0]
    assert second_record["id"] == molecule_id
    assert not second_record["observations"]
    with pytest.raises(PermissionError):
        database.remove_from_group(B, second["group_id"], [molecule_id])


def test_groups_do_not_grant_visibility_to_other_private_scores_or_files(database):
    private = database.import_records(B, imports.parse_smiles("CCO secret-name"), name="Private")
    molecule_id = private["molecule_ids"][0]
    database.add_observation(private["group_id"], molecule_id, "private", "secret-run", {"score": 42})
    public = database.import_records(A, imports.parse_smiles("CCO known-name"), project_id="project-1")
    assert public["molecule_ids"] == [molecule_id]
    row = database.list_molecules(A)["items"][0]
    assert row["name"] == "known-name"
    assert not row["scores"]
    assert not row["observations"]
    assert database.metric_fields(A) == []
    with pytest.raises(PermissionError):
        database.group(private["group_id"], A)
    with pytest.raises(PermissionError):
        database.add_to_group(A, private["group_id"], [molecule_id])


def test_preserves_stereochemistry_stoichiometry_and_group_chain_ids(database):
    first = database.import_records(A, [{"name": "complex", "components": [{"chain": "A", "type": "protein", "sequence": "ACD"}, {"chain": "B", "type": "protein", "sequence": "ACD"}]}])
    renamed = database.import_records(A, [{"name": "renamed chains", "components": [{"chain": "H", "type": "protein", "sequence": "ACD"}, {"chain": "L", "type": "protein", "sequence": "ACD"}]}])
    monomer = database.import_records(A, imports.parse_sequences("ACD"))
    assert first["molecule_ids"] == renamed["molecule_ids"]
    assert first["molecule_ids"] != monomer["molecule_ids"]
    row = database.list_molecules(A, group_id=renamed["group_id"])["items"][0]
    assert [component["chain"] for component in row["components"]] == ["H", "L"]
    stereo = database.import_records(A, imports.parse_smiles("N[C@@H](C)C(=O)O L\nN[C@H](C)C(=O)O D"))
    assert len(stereo["molecule_ids"]) == 2


def test_structure_catalog_exposes_explicit_target_provenance_without_hydrating_complex_identity(database):
    imported = database.import_records(A, imports.parse_sequences(">binder\nACDEFG"), name="Target-conditioned binder")
    molecule_id = imported["molecule_ids"][0]
    database.update_metadata(A, imported["group_id"], molecule_id, {
        "target_fasta": "target.fasta",
        "target_chain_inference": "fasta_sequence_match",
        "target_type": "sequence",
    })
    structure_id = database.attach_structure(
        imported["group_id"], molecule_id,
        {"name": "complex.cif", "content": b"data_test\n", "label": "Saved complex"},
        context={"candidate_id": "candidate-1"},
    )

    record = next(item for item in database.list_structure_records(A)["items"] if item["id"] == molecule_id)
    assert record["source_target_provenance"] is True
    assert len(record["components"]) == 1
    assert record["representations"][0]["id"] == structure_id
    assert database.list_structure_records(B)["items"] == []


def test_structure_catalog_does_not_mark_generic_target_files_as_an_explicit_chain_link(database):
    imported = database.import_records(A, imports.parse_sequences(">binder\nACDEFG"), name="Sequence-only binder")
    molecule_id = imported["molecule_ids"][0]
    database.update_metadata(A, imported["group_id"], molecule_id, {"target_fasta": "target.fasta"})
    database.attach_structure(imported["group_id"], molecule_id,
                              {"name": "binder.cif", "content": b"data_test\n"})

    record = next(item for item in database.list_structure_records(A)["items"] if item["id"] == molecule_id)
    assert record["source_target_provenance"] is False


def test_structure_catalog_can_be_restricted_to_selected_collection_members(database):
    first = database.import_records(A, imports.parse_sequences(">same binder\nACDE"), name="RFD4 run")
    molecule_id = first["molecule_ids"][0]
    first_structure = database.attach_structure(
        first["group_id"], molecule_id, {"name": "rfd4.cif", "content": b"data_rfd4\n"},
        context={"candidate_id": "rfd4-candidate", "model_id": "model-1"},
    )
    second = database.create_group(A, "BindCraft2 run", import_method="user_uploads",
                                   origin={"type": "results_design_collection", "collection_kind": "curated"})
    database.add_to_group(A, second, [molecule_id])
    second_structure = database.attach_structure(
        second, molecule_id, {"name": "bindcraft.cif", "content": b"data_bindcraft\n"},
        context={"candidate_id": "bindcraft-candidate", "model_id": "model-2"},
    )

    selected = database.list_structure_records(
        A, collection="all", group_ids=[first["group_id"]], restrict_structure_groups=True,
    )
    assert [row["id"] for row in selected["items"]] == [molecule_id]
    assert {item["id"] for item in selected["items"][0]["representations"]} == {first_structure}

    # Ordinary Database discovery continues to expose both source models.
    all_models = database.list_structure_records(A, collection="all", group_ids=[first["group_id"]])
    assert {item["id"] for item in all_models["items"][0]["representations"]} == {first_structure, second_structure}


def test_concurrent_legacy_writes_preserve_other_owners_and_additions(tmp_path):
    path = tmp_path / "groups.json"
    path.write_text(json.dumps([{"id": "bob", "name": "Bob", "owner_id": B["id"], "members": []}]))
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda index: save_legacy_groups(path, [{"id": f"a-{index}", "name": str(index), "owner_id": A["id"], "members": []}], A), range(8)))
    assert {row["id"] for row in json.loads(path.read_text())} == {"bob", *[f"a-{index}" for index in range(8)]}
    with pytest.raises(PermissionError):
        save_legacy_groups(path, [{"id": "bob", "name": "changed", "owner_id": B["id"]}], A)


def test_import_table_supports_mixed_chains_and_preserves_numeric_properties():
    data = b"name,sequence,smiles,affinity_nM\ncomplex,ACDE,CCO,12.5\n"
    rows = imports.parse_table(data, "molecules.csv", name_column="name", chain_columns={"sequence": {"chain": "A", "type": "protein", "role": "target"}, "smiles": {"chain": "L", "type": "small_molecule", "role": "lead"}})
    assert rows[0]["metadata"] == {"affinity_nM": 12.5}
    assert imports.canonical_components(rows[0]["components"])[2] == "complex"
    with pytest.raises(ValueError, match="Row 1"):
        imports.parse_table(b"sequence\nINVALID?\n", "bad.csv", chain_columns={"sequence": {"type": "protein", "chain": "A"}})


def test_import_table_keeps_separate_compound_id_name_and_smiles():
    rows = imports.parse_table(
        b"compound_id,compound_name,smiles\n70701262,Bepristat_1a,CCN\n",
        "compounds.csv",
        id_column="compound_id",
        name_column="compound_name",
        chain_columns={"smiles": {"chain": "L", "type": "small_molecule"}},
    )

    assert rows[0]["name"] == "Bepristat_1a"
    assert rows[0]["metadata"]["compound_id"] == "70701262"
    assert rows[0]["components"][0]["smiles"] == "CCN"


def test_import_table_uses_id_as_display_name_when_name_column_is_missing():
    rows = imports.parse_table(
        b"compound_id,smiles\nML359,CCN\n",
        "compounds.csv",
        id_column="compound_id",
        chain_columns={"smiles": {"chain": "L", "type": "small_molecule"}},
    )

    assert rows[0]["name"] == "ML359"
    assert rows[0]["metadata"]["compound_id"] == "ML359"


def test_formula_operations_bound_string_and_numeric_results_before_expansion():
    from mn_molecular_design_studio.molecule_schema import MAX_FORMULA_STRING, evaluate

    assert evaluate('"a" + "b"', {}) == "ab"
    assert evaluate('"x" * 1000000000', {}) == "#VALUE!"
    assert evaluate('[payload] * 1000000000', {"payload": "x"}) == "#VALUE!"
    assert evaluate("[payload]", {"payload": "x" * (MAX_FORMULA_STRING + 1)}) == "#VALUE!"
    assert evaluate("2 ** (2 ** 100)", {}) == "#VALUE!"


def test_fasta_multichain_nucleic_acids_sdf_and_archive(database):
    assert len(imports.parse_sequences(">complex\nACD:EFG\n>single\nKKK")[0]["components"]) == 2
    assert imports.parse_sequences(">dna\nACGTN", kind="dna")[0]["components"][0]["sequence"] == "ACGTN"
    from rdkit import Chem
    molecule = Chem.MolFromSmiles("C[C@H](O)Cl")
    block = Chem.MolToMolBlock(molecule) + "\n$$$$\n"
    records = imports.parse_sdf(block.encode())
    imported = database.import_records(A, records)
    row = database.list_molecules(A)["items"][0]
    assert "@" in row["components"][0]["smiles"]
    assert database.structure_file(A, row["structures"][0]["id"]).is_file()
    with pytest.raises(PermissionError):
        database.structure_file(B, row["structures"][0]["id"])
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("../molecule.sdf", block)
    assert len(imports.parse_archive(buffer.getvalue())) == 1
    assert len(imported["molecule_ids"]) == 1


def test_pagination_metric_sort_and_latest_observations(database):
    result = database.import_records(A, imports.parse_smiles("CCO a\nCC b\nCCN c"))
    for index, molecule in enumerate(result["molecule_ids"]):
        database.add_observation(result["group_id"], molecule, "test", "old", {"score": 100}, created_at="2026-01-01")
        database.add_observation(result["group_id"], molecule, "test", "new", {"score": index}, created_at="2026-02-01")
    page = database.list_molecules(A, sort="test.score", descending=True, limit=2)
    assert [row["name"] for row in page["items"]] == ["c", "b"]
    assert page["next_offset"] == 2
    assert len(page["items"][0]["observations"]) == 2
    assert database.list_molecules(A, metric_filters=[{"field": "test.score", "operator": ">", "value": 10}])["total"] == 0
