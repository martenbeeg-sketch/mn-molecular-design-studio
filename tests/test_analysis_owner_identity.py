import hashlib
import json

import pytest

from mn_molecular_design_studio.cluster import config
from mn_molecular_design_studio.cluster.replication import Replica
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_ingest import ingest_run, repair_local_analysis_owner_nodes
from mn_molecular_design_studio.molecular_database_peer import capture
from test_molecular_database_peer import OWNER, _author_peer_entry, _evaluation_run, _sync_author_records, roots


def local_analysis(roots, monkeypatch):
    node = "11111111-2222-4333-8444-555555555555"
    monkeypatch.setattr(config, "identity", lambda: node)
    database = MolecularDatabase(roots / "owner")
    group, molecule, structure, content = _author_peer_entry(database, node)
    with database.connection() as connection:
        origin = json.loads(connection.execute("SELECT origin FROM molecule_groups WHERE id=?", (group,)).fetchone()[0])
        connection.execute("UPDATE molecule_groups SET origin=? WHERE id=?", (json.dumps({**origin, "workstation": "pc-example"}), group))
    digest = hashlib.sha256(content).hexdigest()
    source = {"schema": "studio.complex_analysis_database_source.v1", "group_id": group,
              "molecule_id": molecule, "structure_id": structure, "structure_digest": digest}
    run = _evaluation_run(roots / "evaluation", node, source, digest)
    assert ingest_run(run, database=database) == 1
    with database.connection() as connection:
        observation = dict(connection.execute("SELECT * FROM molecule_observations WHERE run_id=?", (run.name,)).fetchone())
    return node, database, run, observation, content


def test_local_analysis_uses_node_identity_not_workstation_display_name(roots, monkeypatch):
    node, database, run, observation, content = local_analysis(roots, monkeypatch)
    context = json.loads(observation["context"])
    assert context["source_node"] == context["source_structure_owner_node"] == node
    replica = Replica()
    receiver = MolecularDatabase(roots / "receiver")
    capture(replica, database=database, node_id=node)
    _sync_author_records(replica, node, [receiver])
    with receiver.connection() as connection:
        context = json.loads(connection.execute("SELECT context FROM molecule_observations WHERE run_id=?", (run.name,)).fetchone()[0])
        assert context["source_structure_reference_status"] == "resolved"
        assert connection.execute("SELECT 1 FROM molecule_structures WHERE id=? AND molecule_id=?",
            (context["representation_id"], observation["molecule_id"])).fetchone()


@pytest.mark.parametrize("invalid_digest", [False, True])
def test_legacy_owner_repair_is_verified_scoped_idempotent_and_replayable(roots, monkeypatch, invalid_digest):
    node, database, run, observation, content = local_analysis(roots, monkeypatch)
    context = json.loads(observation["context"])
    context["source_node"] = context["source_structure_owner_node"] = "pc-example"
    if invalid_digest:
        context["source_structure_digest"] = "0" * 64
    with database.connection() as connection:
        connection.execute("UPDATE molecule_observations SET context=? WHERE id=?", (json.dumps(context), observation["id"]))
    if invalid_digest:
        with pytest.raises(ValueError, match="exact source structure"):
            repair_local_analysis_owner_nodes(run.name, database=database, dry_run=False)
        return
    replica = Replica()
    receiver = MolecularDatabase(roots / "receiver")
    capture(replica, database=database, node_id=node)
    _sync_author_records(replica, node, [receiver])
    with receiver.connection() as connection:
        before = json.loads(connection.execute("SELECT context FROM molecule_observations WHERE run_id=?", (run.name,)).fetchone()[0])
        assert before["source_structure_reference_status"] == "owner_record_unavailable"
    assert repair_local_analysis_owner_nodes("another-run", database=database)["matched"] == 0
    assert repair_local_analysis_owner_nodes(run.name, database=database)["matched"] == 1
    assert repair_local_analysis_owner_nodes(run.name, database=database, dry_run=False)["updated"] == 1
    assert repair_local_analysis_owner_nodes(run.name, database=database, dry_run=False)["updated"] == 0
    assert ingest_run(run, database=database, force=True) == 1
    with database.connection() as connection:
        rows = connection.execute("SELECT id,metrics FROM molecule_observations WHERE run_id=?", (run.name,)).fetchall()
        assert [(row["id"], row["metrics"]) for row in rows] == [(observation["id"], observation["metrics"])]
    capture(replica, database=database, node_id=node)
    for attempt in range(2):
        _sync_author_records(replica, node, [receiver])
    with receiver.connection() as connection:
        rows = connection.execute("SELECT context FROM molecule_observations WHERE run_id=?", (run.name,)).fetchall()
        assert len(rows) == 1
        assert json.loads(rows[0]["context"])["source_structure_reference_status"] == "resolved"
