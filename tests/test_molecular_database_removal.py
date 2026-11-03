"""Owner-authorized cross-workstation removal, independent of bulk sync policy."""
import json

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id
from mn_molecular_design_studio import molecular_database_peer as peer
from mn_molecular_design_studio.cluster import config, transport

OWNER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
ADMIN = {"id": "b" * 32, "role": "admin"}
OTHER = {"id": "c" * 32, "role": "user", "_cluster_project_ids": []}
SOURCE = "11111111-1111-4111-8111-111111111111"
TARGET = "22222222-2222-4222-8222-222222222222"


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "target"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "target"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setattr(config, "is_peer", lambda: False)
    monkeypatch.setattr(config, "identity", lambda: SOURCE)
    source = MolecularDatabase(tmp_path / "source")
    target = MolecularDatabase()
    record = {"name": "saved", "components": molecule_imports.parse_sequences("ACDE")[0]["components"],
              "files": [{"name": "model.pdb", "content": b"HEADER    saved Results\nEND\n"}]}
    imported = source.import_records(OWNER, [record])
    group, molecule = imported["group_id"], imported["molecule_ids"][0]
    other = source.create_group(OWNER, "Separately pending")
    source.add_to_group(OWNER, other, [molecule])
    from test_molecular_database_peer import _group_record, _entry_record
    group_record = _group_record(source, group)
    entry = _entry_record(source, group, molecule)
    entry["node_id"] = SOURCE
    peer.apply_record("molecule_group", f"molecule_group:{SOURCE}:{group}", {"node_id": SOURCE, "group": group_record}, SOURCE, database=target)
    peer.apply_record("molecule_entry", f"molecule_entry:{SOURCE}:{group}:{molecule}", entry, SOURCE, database=target)
    return source, target, group, molecule, other


def _source_api(monkeypatch, source):
    from mn_molecular_design_studio import auth
    from mn_molecular_design_studio.cluster.store import Store
    monkeypatch.setattr(auth, "list_users", lambda: [OWNER, ADMIN, OTHER])
    monkeypatch.setattr(Store, "authenticate_node", lambda self, token: {"id": TARGET} if token == "trusted" else (_ for _ in ()).throw(PermissionError("Invalid node")))
    monkeypatch.setattr(peer, "MolecularDatabase", lambda: source)
    app = FastAPI()
    peer.install_routes(app)
    return TestClient(app)


@pytest.mark.parametrize("user", [OWNER, ADMIN])
def test_source_api_rechecks_owner_and_admin_and_preserves_other_membership(workspace, monkeypatch, user):
    source, target, group, molecule, other = workspace
    client = _source_api(monkeypatch, source)
    path = source.list_molecules(OWNER, group_id=group)["items"][0]["representations"][0]["path"]
    response = client.post(f"/api/v1/cluster/node/molecules/groups/{group}/members/remove",
        headers={"authorization": "Bearer trusted"}, json={"user_id": user["id"], "molecule_ids": [molecule]})
    assert response.status_code == 200, response.text
    assert source.list_molecules(OWNER, group_id=group)["total"] == 0
    assert source.list_molecules(OWNER, group_id=other)["total"] == 1
    from pathlib import Path
    assert Path(path).is_file()
    # A retry after a lost response is an authorized idempotent removal.
    assert client.post(f"/api/v1/cluster/node/molecules/groups/{group}/members/remove",
        headers={"authorization": "Bearer trusted"}, json={"user_id": user["id"], "molecule_ids": [molecule]}).status_code == 200
    records = response.json()["mutation"]["records"]
    for row in records:
        peer.apply_record(row["kind"], row["key"], row["value"], SOURCE, database=target)
    assert target.list_molecules(OWNER, group_id=_id("peer-", SOURCE + ":" + group))["total"] == 0
    assert molecule in source.group(group, OWNER)["origin"]["removed_molecule_ids"]


def test_nonowner_and_untrusted_workstation_denied_without_changes(workspace, monkeypatch):
    source, _, group, molecule, _ = workspace
    client = _source_api(monkeypatch, source)
    url = f"/api/v1/cluster/node/molecules/groups/{group}/members/remove"
    assert client.post(url, json={"user_id": OWNER["id"], "molecule_ids": [molecule]}).status_code == 403
    assert client.post(url, headers={"authorization": "Bearer trusted"},
                       json={"user_id": OTHER["id"], "molecule_ids": [molecule]}).status_code == 403
    assert source.list_molecules(OWNER, group_id=group)["total"] == 1


def test_disconnected_source_does_not_mutate_replica(workspace, monkeypatch):
    _, target, group, molecule, _ = workspace
    replica_group = _id("peer-", SOURCE + ":" + group)
    monkeypatch.setattr(config, "load", lambda: {"peers": [{"id": SOURCE, "url": "http://offline"}]})
    monkeypatch.setattr(transport, "node_request", lambda *a, **k: (_ for _ in ()).throw(transport.ClusterUnavailable("offline")))
    with pytest.raises(transport.ClusterUnavailable):
        target.remove_from_group(OWNER, replica_group, [molecule])
    assert target.list_molecules(OWNER, group_id=replica_group)["total"] == 1
    with pytest.raises(PermissionError):
        target.remove_from_group(OTHER, replica_group, [molecule])


def test_explicit_removal_outbox_retries_with_bulk_sharing_disabled(workspace, monkeypatch):
    source, _, group, molecule, _ = workspace
    remote = {"id": TARGET, "url": "http://peer"}
    monkeypatch.setattr(config, "load", lambda: {"share_molecular_database": False, "peers": [remote]})
    monkeypatch.setattr(config, "is_peer", lambda: True)
    monkeypatch.setattr(transport, "node_request", lambda *a, **k: (_ for _ in ()).throw(transport.ClusterUnavailable("offline")))
    assert source.remove_from_group(OWNER, group, [molecule]) == 1
    event_path = next((config.home() / "molecular-removals").glob("*.json"))
    assert json.loads(event_path.read_text())["delivered"] == []
    sent = []
    def send(path, payload, **kwargs):
        sent.append(payload)
        return {"node_id": TARGET, "accepted": len(payload["records"])}
    monkeypatch.setattr(transport, "node_request", send)
    peer.flush_member_removals(remote)
    assert json.loads(event_path.read_text())["delivered"] == [TARGET]
    assert {r["kind"] for r in sent[0]["records"]} == {"molecule_group", "molecule_entry"}
