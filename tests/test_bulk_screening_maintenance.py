import json
from pathlib import Path

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.bulk_screening_maintenance import remove_automatic_bulk_records, capture_bulk_deletions
from mn_molecular_design_studio.molecular_database import MolecularDatabase

USER = {"id": "a" * 32, "role": "admin"}


def test_bulk_cleanup_removes_orphan_tool_inputs_preserves_curated_and_other_runs(tmp_path):
    db = MolecularDatabase(tmp_path / "database")
    group = db.create_group(USER, "Tool inputs", import_method="tool_inputs",
                            origin={"type": "job", "collection_kind": "inputs"})
    saved = db.create_group(USER, "Curated shortlist", origin={"type": "user", "collection_kind": "curated"})
    other = db.create_group(USER, "Other user's inputs", origin={"type": "job", "collection_kind": "inputs"})
    run = {"type": "job", "task_group": "ligand-docking", "run_id": "batch-1"}
    ids = []
    for i, smiles in enumerate(("CCO", "CCN", "CCC")):
        record = molecule_imports.parse_smiles(smiles + f" test-{i}")[0]
        record["source_key"] = f"bulk-{i}"
        result = db.import_records(USER, [record], group_id=group, origin=run)
        ids.append(result["molecule_ids"][0])
        db.add_observation(group, ids[-1], "vina", "batch-1", {"score": -7}, context=run)
    record = molecule_imports.parse_smiles("CCO test-0")[0]
    record["source_key"] = "saved"
    db.import_records(USER, [record], group_id=saved, origin={**run, "saved_selection": True})
    db.add_observation(saved, ids[0], "vina", "batch-1", {"score": -7}, context={**run, "saved_selection": True})
    record = molecule_imports.parse_smiles("CCN test-1")[0]
    record["source_key"] = "other-job"
    db.import_records(USER, [record], group_id=group, origin={**run, "run_id": "batch-2"})
    record = molecule_imports.parse_smiles("CCCC unrelated")[0]
    record["source_key"] = "unrelated"
    unrelated_id = db.import_records(USER, [record], group_id=other, origin={**run, "run_id": "other"})["molecule_ids"][0]
    external = tmp_path / "scores.parquet"
    external.write_bytes(b"external result bytes")
    refs = [("ligand-docking", "batch-1")]
    dry = remove_automatic_bulk_records(refs, USER, owner_id=USER["id"], database=db)
    assert dry["automatic_sources"] == 3 and dry["dry_run"]
    assert db.list_molecules(USER, group_id=group)["total"] == 3
    result = remove_automatic_bulk_records(refs, USER, owner_id=USER["id"], database=db, dry_run=False)
    assert result["removed_members"] == 2 and result["removed_identities"] == 1
    assert result["removed_observations"] == 3
    assert db.list_molecules(USER, group_id=group)["total"] == 1
    assert db.list_molecules(USER, group_id=saved)["items"][0]["id"] == ids[0]
    assert db.list_molecules(USER, group_id=other)["items"][0]["id"] == unrelated_id
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == 1
    assert external.read_bytes() == b"external result bytes"


def test_cleanup_publishes_deletions_without_recapturing_or_removing_other_authors(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.cluster.replication import Replica
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    config.home().mkdir(exist_ok=True, parents=True)
    db = MolecularDatabase(tmp_path / "db")
    group = db.create_group(USER, "Kept")
    mid = db.import_records(USER, molecule_imports.parse_smiles("CCO keep"), group_id=group)["molecule_ids"][0]
    replica = Replica()
    with replica.db() as c:
        for i, (key, author, molecule) in enumerate((("gone", "local", "gone"), ("keep", "local", mid), ("other", "remote", "gone")), 1):
            c.execute("INSERT INTO records VALUES (?, 'molecule_entry', ?, ?, ?, ?)",
                      (key, author, i, json.dumps({"group_id": group, "molecule": {"id": molecule}}), i))
    assert capture_bulk_deletions(replica, database=db, node_id="local") == 1
    with replica.db() as c:
        assert c.execute("SELECT value FROM records WHERE key='gone'").fetchone()[0] == "null"
        assert c.execute("SELECT value FROM records WHERE key='other'").fetchone()[0] != "null"
        assert c.execute("SELECT counter FROM records WHERE key='gone'").fetchone()[0] == 4
    assert capture_bulk_deletions(replica, database=db, node_id="local") == 0
