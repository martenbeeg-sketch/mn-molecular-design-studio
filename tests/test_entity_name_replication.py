from fastapi.testclient import TestClient
import pytest

from mn_molecular_design_studio import auth, entity_names as names, workspace
from mn_molecular_design_studio.cluster import config
from mn_molecular_design_studio.cluster.replication import Replica
from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_peer import capture
from test_peer_workstations import PASSWORD, account, peers, transfer


def test_project_sql_names_sync_and_source_ownership_survives_relay(peers):
    nodes, on = peers
    user, project = account(on)
    transfer(on, 0, 1)
    with on(0):
        renamed = names.rename_project(user, project["id"], "Tau design study", expected_revision="")
    pages = transfer(on, 0, 1)
    assert any(row["kind"] == "entity_name" for row in pages)
    with on(1):
        current = workspace.load_organization()["projects"][0]
        assert current["name"] == renamed["name"]
        assert current["name_aliases"] == ["Experiment"]
        assert current["name_revision"] == renamed["name_revision"]
        with pytest.raises(PermissionError, match="source workstation"):
            names.rename_project(user, project["id"], "Remote overwrite", expected_revision=current["name_revision"])
        names.publish()
        assert all(row["value"] is not None for row in Replica().export()["records"] if row["kind"] == "entity_name")
    transfer(on, 1, 2)
    with on(2):
        assert workspace.load_organization()["projects"][0]["name"] == "Tau design study"
        with MolecularDatabase().connection() as connection:
            label = next(row for row in names.records(connection) if row["kind"] == "project")
            assert label["authority_node"] == nodes[0]["id"]


def test_target_sql_names_follow_shared_exact_component_without_changing_members(peers):
    nodes, on = peers
    user, project = account(on)
    with on(0):
        database = MolecularDatabase()
        imported = database.import_records(user, [{"name": "Original receptor", "metadata": {"target_identity_type": "sequence"},
            "components": [{"type": "protein", "chain": "A", "sequence": "ACDE"}]}], name="Reference",
            import_method="user_uploads", project_id=project["id"], origin={"collection_kind": "reference"})
        subject = component_subjects(build_browser_index(database, user), role="targets")[0]
        names.rename_target(database, user, subject["id"], "project:" + project["id"], "Tau R3", expected_revision="")
        capture(Replica())
    transfer(on, 0, 1)
    with on(1):
        database = MolecularDatabase()
        index = build_browser_index(database, user)
        assert component_subjects(index, role="targets")[0]["name"] == "Tau R3"
        assert set(index["molecules"][imported["molecule_ids"][0]]["member_names"].values()) == {"Original receptor"}
        names.publish()
    transfer(on, 1, 2)
    with on(2):
        assert component_subjects(build_browser_index(MolecularDatabase(), user), role="targets")[0]["name"] == "Tau R3"


def test_name_export_negotiation_and_sharing_disabled_receive(peers):
    nodes, on = peers
    user, project = account(on)
    from mn_molecular_design_studio import api
    with on(0):
        names.rename_project(user, project["id"], "Renamed", expected_revision="")
        replica = Replica()
        compatible = replica.export(include_entity_names=False)
        assert all(row["kind"] != "entity_name" for row in compatible["records"])
        names_page = [row for row in replica.export()["records"] if row["kind"] == "entity_name"]
        assert names_page
        with TestClient(api.app) as client:
            headers = {"Authorization": "Bearer " + nodes[1]["token"]}
            assert all(row["kind"] != "entity_name" for row in client.get(
                "/api/v1/cluster/node/peer/records", headers=headers).json()["records"])
            assert any(row["kind"] == "entity_name" for row in client.get(
                "/api/v1/cluster/node/peer/records?include_entity_names=true", headers=headers).json()["records"])
    with on(1):
        config.save({"share_molecular_database": False})
        with TestClient(api.app) as client:
            response = client.post("/api/v1/cluster/node/peer/records", json={"records": names_page},
                headers={"Authorization": "Bearer " + nodes[0]["token"]})
            assert response.status_code == 200 and response.json()["accepted"] == 0


def test_project_rename_api_validates_revision_and_returns_conflict(peers):
    nodes, on = peers
    user, project = account(on)
    from mn_molecular_design_studio import api
    with on(0), TestClient(api.app) as client:
        _, token, _ = auth.authenticate("scientist", PASSWORD)
        headers = {"Authorization": "Bearer " + token}
        path = "/api/v1/cluster/projects/" + project["id"] + "/name"
        assert client.post(path, headers=headers, json={"name": "New"}).status_code == 422
        response = client.post(path, headers=headers, json={"name": "New", "expected_revision": ""})
        assert response.status_code == 200
        stale = client.post(path, headers=headers, json={"name": "Stale", "expected_revision": ""})
        assert stale.status_code == 409


def test_stale_name_replay_after_restart_does_not_undo_newer_label(peers):
    _, on = peers
    user, project = account(on)
    with on(0):
        first = names.rename_project(user, project["id"], "First edit", expected_revision="")
        old = [row for row in Replica().export()["records"] if row["kind"] == "entity_name"]
        names.rename_project(user, project["id"], "Second edit", expected_revision=first["name_revision"])
    transfer(on, 0, 1)
    with on(1):
        assert Replica().receive(old) == 0
        project = workspace.load_organization()["projects"][0]
        assert project["name"] == "Second edit"
        assert {"Experiment", "First edit"} <= set(project["name_aliases"])


def test_dataset_and_membership_properties_replicate_by_source_id(peers):
    from mn_molecular_design_studio import entity_properties as properties
    nodes, on = peers
    user, project = account(on)
    with on(0):
        database = MolecularDatabase()
        imported = database.import_records(user, [{"name": "Compound A", "components": [
            {"type": "small_molecule", "chain": "L", "smiles": "CCO"}]}], name="Library", project_id=project["id"])
        group_id, molecule_id = imported["group_id"], imported["molecule_ids"][0]
        for snapshot, label, values in [
            (properties.dataset_snapshot(database, user, group_id), "HY-L126 study", {"description": "Imported study", "citation": "Paper A"}),
            (properties.entry_snapshot(database, user, group_id, molecule_id), "Catalog A", {"batch": "Lot 12", "notes": "Only this dataset"}),
        ]:
            properties.save_group_properties(database, user, snapshot, label, values, [], expected_revision=snapshot["revision"])
        from mn_molecular_design_studio.data_browser import saved_results
        from mn_molecular_design_studio.database_browser import subject_results
        index = build_browser_index(database, user, hide_smoke_test=False, include_related=True)
        scope = {"id": "all", "group_ids": set(index["groups"]), "molecule_ids": {molecule_id}}
        outcome = saved_results(index, subject_results(database, index, scope))[0]
        snapshot = properties.result_snapshot(database, user, outcome)
        properties.save_group_properties(database, user, snapshot, "Reference outcome", {"notes": "Result-level annotation"}, [],
            expected_revision=snapshot["revision"], outcome=outcome)
        capture(Replica())
        legacy = [row for row in Replica().export(include_entity_properties=False)["records"] if row["kind"] == "entity_name"]
        assert all(row["value"]["kind"] in {"project", "target"} and "properties" not in row["value"] for row in legacy)
    transfer(on, 0, 1)
    with on(1):
        database = MolecularDatabase()
        group = next(row for row in database.groups(user) if row["origin"].get("_replica_group") == group_id)
        assert group["name"] == "HY-L126 study" and group["entity_properties"]["citation"] == "Paper A"
        snapshot = properties.entry_snapshot(database, user, group["id"], molecule_id)
        assert snapshot["name"] == "Catalog A" and snapshot["properties"]["batch"] == "Lot 12"
        assert not snapshot["writable"]
        index = build_browser_index(database, user, hide_smoke_test=False, include_related=True)
        scope = {"id": "all", "group_ids": set(index["groups"]), "molecule_ids": {molecule_id}}
        outcome = saved_results(index, subject_results(database, index, scope))[0]
        assert outcome["display_name"] == "Reference outcome"
        assert properties.result_snapshot(database, user, outcome)["properties"]["notes"] == "Result-level annotation"
        names.publish()
    transfer(on, 1, 2)
    with on(2):
        database = MolecularDatabase()
        group = next(row for row in database.groups(user) if row["origin"].get("_replica_group") == group_id)
        assert properties.entry_snapshot(database, user, group["id"], molecule_id)["properties"]["notes"] == "Only this dataset"


def test_legacy_name_only_replay_does_not_erase_properties(peers):
    _, on = peers
    user, project = account(on)
    with on(0):
        first = names.rename_project(user, project["id"], "Named", expected_revision="", properties={"description": "Research", "notes": "Keep me"})
        modern = [row for row in Replica().export()["records"] if row["kind"] == "entity_name"]
        legacy = [row for row in Replica().export(include_entity_properties=False)["records"] if row["kind"] == "entity_name"]
        assert "properties" in modern[0]["value"] and "properties" not in legacy[0]["value"]
    transfer(on, 0, 1)
    with on(1):
        legacy[0]["counter"] += 100
        legacy[0]["value"]["name"] = "Legacy rename"
        legacy[0]["value"]["revision"] = "b" * 32
        Replica().receive(legacy)
        current = workspace.load_organization()["projects"][0]
        assert current["name"] == "Legacy rename" and current["description"] == "Research"
        assert current["entity_properties"]["notes"] == "Keep me"


def test_equal_local_dataset_ids_on_two_nodes_do_not_share_properties(peers):
    from mn_molecular_design_studio import entity_properties as properties
    nodes, on = peers
    user, project = account(on)
    transfer(on, 0, 1)
    for node in (0, 1):
        with on(node):
            database = MolecularDatabase()
            group_id = database.create_group(user, "Library", group_id="same-local-id", project_id=project["id"])
            snapshot = properties.dataset_snapshot(database, user, group_id)
            properties.save_group_properties(database, user, snapshot, "Library on " + str(node), {"notes": "Node " + str(node)}, [],
                                              expected_revision=snapshot["revision"])
            capture(Replica())
    transfer(on, 0, 1)
    transfer(on, 1, 0)
    for node in (0, 1):
        with on(node):
            database = MolecularDatabase()
            assert database.group("same-local-id", user)["name"] == "Library on " + str(node)
            with database.connection() as connection:
                rows = [row for row in names.records(connection) if row["kind"] == "dataset"]
            assert {row["name"] for row in rows} == {"Library on 0", "Library on 1"}
