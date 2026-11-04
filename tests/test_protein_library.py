from __future__ import annotations

import json
from pathlib import Path

import pytest

from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.protein_library import collect_run_inputs, audit_database

USER = {"id": "a" * 32, "role": "admin"}
PDB = ("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 50.00           C  \n"
       "ATOM      2  CA  CYS A   2       3.800   0.000   0.000  1.00 50.00           C  \n"
       "ATOM      3  CA  ASP A   3       7.600   0.000   0.000  1.00 50.00           C  \nEND\n")


@pytest.fixture
def source_run(tmp_path, monkeypatch):
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID
    from mn_molecular_design_studio import protein_library

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setattr(protein_library, "_FETCH_FAILURES", {})
    run = tmp_path / "runs/design/input-probe"
    run.mkdir(parents=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({
        "application": STUDIO_SOURCE_ID, "studio_owner_id": USER["id"], "project_id": ""}))
    (run / "metadata.json").write_text(json.dumps({"status": "queued"}))
    return run, MolecularDatabase(tmp_path / "database")


def test_source_structure_is_portable_reusable_and_deduplicated_before_job_completes(source_run):
    run, database = source_run
    (run / "source.pdb").write_text(PDB)
    (run / "input.json").write_text(json.dumps({"inputs": {"target_structure": "source.pdb"}}))
    first = collect_run_inputs(run, database=database, fetch_references=False)
    second = collect_run_inputs(run, database=database, fetch_references=False)
    assert first["molecule_ids"] == second["molecule_ids"]
    rows = database.list_molecules(USER, collection="library")["items"]
    assert len(rows) == 1
    assert rows[0]["groups"][0]["name"] == "Protein Library"
    assert len(rows[0]["representations"]) == 1
    representation = rows[0]["representations"][0]
    assert representation["context"]["portable"] is True
    assert representation["context"]["input_role"] == "target"
    assert not rows[0]["observations"]
    (run / "source.pdb").unlink()
    assert database.structure_file(USER, representation["id"]).read_text() == PDB


def test_restricted_collection_and_idr_carriers_are_not_promoted(source_run):
    run, database = source_run
    (run / "source.pdb").write_text(PDB)
    (run / "input.json").write_text(json.dumps({"inputs": {"target_structure": "source.pdb"}}))
    (run / "metadata.json").write_text(json.dumps({"database_collection_policy": "manual"}))
    assert collect_run_inputs(run, database=database)["molecule_ids"] == []
    (run / "metadata.json").write_text("{}")
    (run / "idr_sequence_carrier.cif").write_text("data_idr_sequence_carrier\n")
    (run / "input.json").write_text(json.dumps({"inputs": {"target_structure": "idr_sequence_carrier.cif"}}))
    assert collect_run_inputs(run, database=database)["molecule_ids"] == []


def test_validation_input_cannot_promote_an_unaccepted_design(source_run):
    run, database = source_run
    (run / "input.json").write_text(json.dumps({"inputs": {"tasks": [{"protein_chains": [
        {"id": "A", "sequence": "ACDEFG", "role": "binder"}]}]},
        "params": {"validation_sources": [{"candidate_id": "not-yet-accepted"}]}}))
    assert collect_run_inputs(run, database=database)["molecule_ids"] == []


def test_idr_target_sequence_is_reusable_without_a_structural_template(source_run):
    run, database = source_run
    (run / "input.json").write_text(json.dumps({"inputs": {"target_sequence": "ACDEFGHIKLMNPQ"}}))
    collect_run_inputs(run, database=database, fetch_references=False)
    rows = database.list_molecules(USER, collection="library")["items"]
    assert len(rows) == 1
    assert rows[0]["components"][0]["sequence"] == "ACDEFGHIKLMNPQ"
    assert not rows[0]["representations"]
    assert rows[0]["metadata"]["target_identity_type"] == "sequence_only"


def test_submission_registration_keeps_source_even_when_job_is_queued(source_run, monkeypatch):
    from mn_molecular_design_studio import molecular_database, molecular_database_ingest

    run, database = source_run
    (run / "source.pdb").write_text(PDB)
    (run / "input.json").write_text(json.dumps({"inputs": {"target_structure": "source.pdb"}}))
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(molecular_database_ingest, "ingest_run", lambda *a, **k: pytest.fail("Queued output must not be ingested"))
    molecular_database_ingest.record_completed_run(run)
    assert len(database.list_molecules(USER, collection="library")["items"]) == 1
    assert json.loads((run / "metadata.json").read_text())["status"] == "queued"


def test_source_library_remains_visible_after_peer_replication(source_run, monkeypatch, tmp_path):
    from mn_molecular_design_studio.cluster.replication import Replica
    from mn_molecular_design_studio.molecular_database_peer import capture, apply_record, materialize_portable_structures
    from mn_molecular_design_studio.cluster import transport

    run, database = source_run
    (run / "source.pdb").write_text(PDB)
    (run / "input.json").write_text(json.dumps({"inputs": {"target_structure": "source.pdb"}}))
    collected = collect_run_inputs(run, database=database, fetch_references=False)
    replica = Replica()
    capture(replica, database=database, node_id="source-node")
    with replica.db() as connection:
        records = connection.execute("SELECT key,kind,value FROM records").fetchall()
    target = MolecularDatabase(tmp_path / "peer-database")
    for row in sorted(records, key=lambda r: r["kind"] != "molecule_group"):
        if row["key"].startswith(row["kind"] + ":source-node:"):
            apply_record(row["kind"], row["key"], json.loads(row["value"]), "source-node", database=target)
    rows = target.list_molecules(USER, collection="library")["items"]
    assert [row["id"] for row in rows] == collected["molecule_ids"]
    assert rows[0]["groups"][0]["name"] == "Protein Library"
    monkeypatch.setattr(transport, "node_request", lambda *a, **k: PDB.encode())
    assert materialize_portable_structures(database=target, peer={"id": "source-node", "url": "https://source"}) == 1
    saved = target.list_molecules(USER, collection="library")["items"][0]
    assert target.structure_file(USER, saved["representations"][0]["id"]).read_text() == PDB


def test_pdb_sequence_source_is_separate_from_prediction_and_preserves_metadata(source_run, monkeypatch):
    from mn_molecular_design_studio import molecule_imports
    from mn_molecular_design_studio.target_identity import _pdb_example_sequences

    run, database = source_run
    sequences = list(_pdb_example_sequences()["8JYR"].elements())
    task = {"protein_chains": [{"id": chr(65 + i), "sequence": sequence, "role": "unspecified"}
                               for i, sequence in enumerate(sequences)]}
    (run / "input.json").write_text(json.dumps({"inputs": {"tasks": [task]}}))
    fetched = []
    import gemmi
    cif = gemmi.read_pdb_string(PDB).make_mmcif_document().as_string().encode()
    monkeypatch.setattr(molecule_imports, "fetch_identifier", lambda code, source: fetched.append(code) or ("8jyr.cif", cif))
    first = collect_run_inputs(run, database=database)
    assert not first["errors"]
    rows = database.list_molecules(USER, collection="library")["items"]
    system = next(row for row in rows if row["metadata"].get("sequence_source_reference"))
    reference_id = system["metadata"]["sequence_source_reference_molecule_id"]
    assert reference_id != system["id"]
    assert not system["representations"]
    assert not any(component.get("role") == "target" for component in system["components"])
    database.update_metadata(USER, system["groups"][0]["id"], system["id"], {"design_target_molecule_id": "existing-target", "acceptance_status": "accepted"})
    second = collect_run_inputs(run, database=database)
    assert not second["errors"]
    assert fetched == ["8JYR"]
    saved = database.list_molecules(USER, molecule_ids=[system["id"]])["items"][0]
    assert saved["metadata"]["design_target_molecule_id"] == "existing-target"
    assert saved["metadata"]["acceptance_status"] == "accepted"
    audit = audit_database(database.path)
    assert audit["counts"]["in_protein_library"] == 2
