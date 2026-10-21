from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id, _json
from mn_molecular_design_studio.molecular_database_peer import apply_record, capture, materialize_portable_structures
from mn_molecular_design_studio.cluster.replication import Replica
from mn_molecular_design_studio import molecule_imports


SOURCE = "source-node-a"
OTHER_SOURCE = "source-node-b"
OWNER = {"id": "owner-1", "role": "user", "_cluster_project_ids": ["project-current", "project-next"]}
OUTSIDER = {"id": "owner-2", "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def roots(tmp_path, monkeypatch):
    # Both workspaces and the replica metadata store must live under pytest's
    # temporary root; avoid reading or writing a developer's live Studio DB.
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "runtime"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    return tmp_path


def _capture_rows(replica, kind=None):
    with replica.db() as connection:
        query = "SELECT key,kind,value FROM records"
        rows = connection.execute(query).fetchall()
    return {row["key"]: (row["kind"], json.loads(row["value"])) for row in rows
            if kind is None or row["kind"] == kind}


def _sync_all(replica, source_database, target_database, author=SOURCE):
    capture(replica, database=source_database, node_id=author)
    records = _capture_rows(replica)
    for key, (kind, value) in sorted(records.items(), key=lambda item: (item[1][0] != "molecule_group", item[0])):
        if key.startswith(kind + ":" + author + ":"):
            apply_record(kind, key, value, author, database=target_database)


def test_saved_comparison_keeps_canonical_evaluations_across_sync(roots):
    source, target = MolecularDatabase(roots / "source-evals"), MolecularDatabase(roots / "target-evals")
    imported = source.import_records(OWNER, molecule_imports.parse_sequences("ACDE"))
    molecule, group = imported["molecule_ids"][0], imported["group_id"]
    source.add_observation(group, molecule, "test", "first", {"score": 1}, context={"target_id": "t"}, created_at="2026-01-01")
    source.add_observation(group, molecule, "test", "second", {"score": 2}, context={"target_id": "t"}, created_at="2026-02-01")
    saved = source.save_comparison(OWNER, "Frozen", [molecule], evaluation_filters={"run_id": "first"})
    _sync_all(Replica(), source, target)
    replica_group = _id("peer-", SOURCE + ":" + saved["group_id"])
    restored = target.list_molecules(OWNER, group_id=replica_group)["items"]
    assert len(restored) == 1
    assert restored[0]["scores"]["test.score"] == 1
    assert restored[0]["selected_evaluation"]["id"] == source.list_molecules(OWNER, group_id=saved["group_id"])["items"][0]["selected_evaluation"]["id"]


def test_synchronized_generator_observations_project_same_common_values(roots):
    source = MolecularDatabase(roots / "generator-source")
    target = MolecularDatabase(roots / "generator-target")
    imported = source.import_records(OWNER, molecule_imports.parse_sequences("ACDE"), name="Generated designs")
    source.add_observation(imported["group_id"], imported["molecule_ids"][0], "bindcraft2_generator", "generation-1",
        {"pTM": .8, "Target_pLDDT": .9}, context={"candidate": {"source_tool": "bindcraft2_generator", "stage": "generation.backbone_sequence", "binder_sequence": "ACDE"}})
    replica = Replica()
    _sync_all(replica, source, target)
    row = target.list_molecules(OWNER)["items"][0]
    assert row["scores"]["common.confidence.generator.complex.ptm"] == .8
    assert row["scores"]["common.confidence.generator.target.plddt"] == 90
    assert row["observations"][0]["metrics"]["Target_pLDDT"] == .9
    assert "common.design.sequence_length" in target.metric_fields(OWNER)


def test_composed_comparison_preserves_subject_and_sources_across_peer_sync(roots):
    from mn_molecular_design_studio.comparisons import compare
    source, target = MolecularDatabase(roots / "source-composition"), MolecularDatabase(roots / "target-composition")
    imported = source.import_records(OWNER, molecule_imports.parse_sequences("ACDE"))
    molecule, group = imported["molecule_ids"][0], imported["group_id"]
    context = {"target_id": "t", "structure_digest": "shared-digest", "target_structure_digest": "target-digest"}
    source.add_observation(group, molecule, "prediction", "p", {"score": .9}, context=context)
    source.add_observation(group, molecule, "quality", "q", {"solubility": .5}, context=context)
    settings = {"combine_same_subject": True, "fields": ["prediction.score", "quality.solubility"],
                "directions": {"prediction.score": "higher", "quality.solubility": "higher"}}
    saved = source.save_comparison(OWNER, "Combined", [molecule], settings=settings)
    expected = compare(source.evaluation_rows(OWNER, group_id=saved["group_id"]), **settings)
    _sync_all(Replica(), source, target)
    replica_group = _id("peer-", SOURCE + ":" + saved["group_id"])
    actual = compare(target.evaluation_rows(OWNER, group_id=replica_group), **settings)
    assert len(actual) == 1
    assert actual.iloc[0]["Source evaluations"] == expected.iloc[0]["Source evaluations"]
    assert actual.iloc[0]["Evaluation ID"] == expected.iloc[0]["Evaluation ID"]
    assert actual.iloc[0]["Comparison rank"] == 1


def test_portable_structure_is_copied_to_peer_database_and_survives_source_offline(roots, monkeypatch):
    source = MolecularDatabase(roots / "portable-source")
    target = MolecularDatabase(roots / "portable-target")
    content = b"HEADER    portable protein structure\nEND\n"
    imported = source.import_records(
        OWNER,
        [{"name": "portable", "components": [{"chain": "A", "type": "protein", "sequence": "ACDE"}],
          "files": [{"name": "model.pdb", "content": content}]}],
        name="Portable structures", project_id="project-current",
    )
    molecule_id = imported["molecule_ids"][0]
    with source.connection() as connection:
        source_structure_id = connection.execute(
            "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=?",
            (imported["group_id"], molecule_id),
        ).fetchone()[0]
    portable = source.make_structure_portable(OWNER, source_structure_id)
    assert Path(portable["path"]).is_file()
    assert portable["sha256"] == hashlib.sha256(content).hexdigest()

    replica = Replica()
    _sync_all(replica, source, target)
    peer_group_id = _id("peer-", SOURCE + ":" + imported["group_id"])
    with target.connection() as connection:
        target_structure = dict(connection.execute(
            "SELECT * FROM molecule_structures WHERE group_id=? AND molecule_id=?",
            (peer_group_id, molecule_id),
        ).fetchone())
    assert target_structure["path"].startswith("cluster-library://")

    from mn_molecular_design_studio.cluster import transport
    monkeypatch.setattr(transport, "node_request", lambda *args, **kwargs: content)
    assert materialize_portable_structures(database=target, peer={"id": SOURCE, "url": "https://source"}) == 1
    with target.connection() as connection:
        target_structure = dict(connection.execute(
            "SELECT * FROM molecule_structures WHERE id=?", (target_structure["id"],)
        ).fetchone())
    assert not target_structure["path"].startswith("cluster-library://")
    assert target.structure_file(OWNER, target_structure["id"]).read_bytes() == content

    # A later metadata refresh must not replace the local portable copy with a
    # remote URI; that copy remains usable if the source workstation is offline.
    with source.connection() as connection:
        connection.execute("UPDATE molecule_groups SET name=? WHERE id=?", ("Renamed portable structures", imported["group_id"]))
    _sync_all(replica, source, target)
    with target.connection() as connection:
        refreshed = dict(connection.execute(
            "SELECT * FROM molecule_structures WHERE id=?", (target_structure["id"],)
        ).fetchone())
    assert not refreshed["path"].startswith("cluster-library://")
    assert target.structure_file(OWNER, refreshed["id"]).read_bytes() == content


def _make_source_with_entry(root):
    database = MolecularDatabase(root)
    imported = database.import_records(
        OWNER,
        [{"name": "ethanol", "components": molecule_imports.parse_smiles("CCO")[0]["components"],
          "files": [{"name": "model.pdb", "content": b"HEADER    source model\nEND\n"},],
          "metadata": {"tag": "original"}}],
        name="Reusable molecules", project_id="project-current", import_method="tool_output",
        origin={"type": "curated_collection", "curated": True, "collection_kind": "curated"},
    )
    group_id, molecule_id = imported["group_id"], imported["molecule_ids"][0]
    observation_id = database.add_observation(
        group_id, molecule_id, "vina", "run-1", {"score": -7.2},
        context={"settings": {"exhaustiveness": 8}}, source_key="score-source-1",
    )
    return database, group_id, molecule_id, observation_id


def _group_record(database, group_id):
    with database.connection() as connection:
        return dict(connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group_id,)).fetchone())


def _entry_record(database, group_id, molecule_id, *, digest=True):
    with database.connection() as connection:
        molecule = dict(connection.execute("SELECT * FROM molecules WHERE id=?", (molecule_id,)).fetchone())
        member = dict(connection.execute("SELECT * FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone())
        structures = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_structures WHERE group_id=? AND molecule_id=? ORDER BY id", (group_id, molecule_id))]
        observations = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_observations WHERE group_id=? AND molecule_id=? ORDER BY id", (group_id, molecule_id))]
        sources = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_sources WHERE group_id=? AND molecule_id=? ORDER BY source_key", (group_id, molecule_id))]
    for structure in structures:
        path = structure["path"]
        data = Path(path).read_bytes() if path else b""
        structure["path"] = ""
        structure["sha256"] = hashlib.sha256(data).hexdigest() if digest and data else ""
        context = json.loads(structure["context"])
        target_path = context.pop("target_path", None)
        if target_path:
            with database.connection() as connection:
                target = connection.execute("SELECT id FROM molecule_structures WHERE path=? LIMIT 1", (target_path,)).fetchone()
            if target:
                context["target_structure_id"] = target[0]
        structure["context"] = _json(context)
    return {"node_id": SOURCE, "group_id": group_id, "molecule": molecule, "member": member,
            "structures": structures, "observations": observations, "sources": sources}


def test_group_snapshot_round_trips_schema_curated_flags_and_project_acl(roots):
    source = MolecularDatabase(roots / "source")
    replica = MolecularDatabase(roots / "replica")
    group_id = source.create_group(
        OWNER, "Curated collection", project_id="project-current", import_method="tool_output",
        origin={"type": "curated_collection", "curated": True, "collection_kind": "curated", "published": True},
        schema={"fields": [{"name": "confidence", "type": "number", "required": True}]},
    )
    transport = Replica()
    capture(transport, database=source, node_id=SOURCE)
    key = f"molecule_group:{SOURCE}:{group_id}"
    kind, value = _capture_rows(transport)[key]
    apply_record(kind, key, value, SOURCE, database=replica)

    peer_group_id = _id("peer-", SOURCE + ":" + group_id)
    group = replica.group(peer_group_id, OWNER)
    assert group["project_id"] == "project-current"
    assert group["import_method"] == "tool_output"
    assert group["schema"] == {"fields": [{"name": "confidence", "type": "number", "required": True}]}
    assert group["origin"]["curated"] is True
    assert group["origin"]["collection_kind"] == "curated"
    assert replica.groups(OUTSIDER, scope="organization") == []

    # Moving the source collection updates the replica's ACL fields and kind.
    with source.connection() as connection:
        connection.execute(
            "UPDATE molecule_groups SET project_id=?,import_method=?,updated_at=? WHERE id=?",
            ("project-next", "user_uploads", "2026-10-20T12:00:00Z", group_id),
        )
    capture(transport, database=source, node_id=SOURCE)
    kind, value = _capture_rows(transport)[key]
    apply_record(kind, key, value, SOURCE, database=replica)
    updated = replica.group(peer_group_id, OWNER)
    assert updated["project_id"] == "project-next"
    assert updated["import_method"] == "user_uploads"
    assert updated["origin"]["curated"] is True
    assert replica.groups(OUTSIDER, scope="organization") == []

    with pytest.raises(PermissionError):
        apply_record(kind, key, {**value, "node_id": OTHER_SOURCE}, OTHER_SOURCE, database=replica)
    assert replica.group(peer_group_id, OWNER)["project_id"] == "project-next"


def test_entry_update_reconciles_deleted_relations_and_keeps_other_records(roots):
    source, group_id, molecule_id, observation_id = _make_source_with_entry(roots / "source")
    replica = MolecularDatabase(roots / "replica")
    transport = Replica()
    _sync_all(transport, source, replica)
    peer_group_id = _id("peer-", SOURCE + ":" + group_id)

    with replica.connection() as connection:
        # Preserve a local relation in the same pair and a second peer group's
        # membership/score for the globally shared molecule identity.
        connection.execute("INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)",
                           ("local-score", molecule_id, peer_group_id, "local", "local-run", "{}", "{}", "now"))
        other_group_id = "other-source-group"
        other_peer_group_id = _id("peer-", OTHER_SOURCE + ":" + other_group_id)
        connection.execute("INSERT INTO molecule_groups VALUES (?,?,?,?,?,?,?,?,?)",
                           (other_peer_group_id, "Other source", OUTSIDER["id"], "", "tool_output",
                            _json({"_replica_node": OTHER_SOURCE, "_replica_group": other_group_id}), "{}", "now", "now"))
        connection.execute("INSERT INTO molecule_members VALUES (?,?,?,?,?,?)",
                           (other_peer_group_id, molecule_id, "Other peer's member", "[]", "{}", "now"))
        connection.execute("INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)",
                           ("other-peer-score", molecule_id, other_peer_group_id, "other", "other-run", "{}",
                            _json({"_replica_node": OTHER_SOURCE, "source_node": OTHER_SOURCE}), "now"))

    # An update to existing relation IDs changes their values without leaving
    # stale duplicates. This also updates the member's source-authoritative name.
    with source.connection() as connection:
        connection.execute("UPDATE molecule_members SET name=?,metadata=? WHERE group_id=? AND molecule_id=?",
                           ("ethanol revised", _json({"tag": "revised"}), group_id, molecule_id))
        connection.execute("UPDATE molecule_structures SET label=? WHERE group_id=? AND molecule_id=?",
                           ("revised model", group_id, molecule_id))
        connection.execute("UPDATE molecule_observations SET metrics=?,context=? WHERE id=?",
                           (_json({"score": -8.1}), _json({"settings": {"exhaustiveness": 16}}), observation_id))
        connection.execute("UPDATE molecule_sources SET provenance=? WHERE group_id=? AND molecule_id=?",
                           (_json({"curated": True}), group_id, molecule_id))
    capture(transport, database=source, node_id=SOURCE)
    entry_key = f"molecule_entry:{SOURCE}:{group_id}:{molecule_id}"
    kind, value = _capture_rows(transport)[entry_key]
    apply_record(kind, entry_key, value, SOURCE, database=replica)
    with replica.connection() as connection:
        member = connection.execute("SELECT name,metadata FROM molecule_members WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()
        structure = connection.execute("SELECT label FROM molecule_structures WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()
        score = connection.execute("SELECT metrics FROM molecule_observations WHERE group_id=? AND molecule_id=? AND namespace='vina'", (peer_group_id, molecule_id)).fetchone()
        provenance = connection.execute("SELECT provenance FROM molecule_sources WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()
    assert member["name"] == "ethanol revised"
    assert json.loads(member["metadata"]) == {"tag": "revised"}
    assert structure["label"] == "revised model"
    assert json.loads(score["metrics"]) == {"score": -8.1}
    assert json.loads(provenance["provenance"])["curated"] is True

    # Source now removes every copied relation. A local relation in the same
    # group/molecule and the other peer group's rows must survive reconciliation.
    with source.connection() as connection:
        for table in ("molecule_sources", "molecule_structures", "molecule_observations"):
            connection.execute(f"DELETE FROM {table} WHERE group_id=? AND molecule_id=?", (group_id, molecule_id))
    capture(transport, database=source, node_id=SOURCE)
    kind, value = _capture_rows(transport)[entry_key]
    apply_record(kind, entry_key, value, SOURCE, database=replica)
    with replica.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_structures WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_sources WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE id='local-score'").fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE id='other-peer-score'").fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id=?", (other_peer_group_id, molecule_id)).fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()[0] == 1


def test_entry_tombstone_and_malformed_snapshot_are_atomic_and_author_scoped(roots):
    source, group_id, molecule_id, _ = _make_source_with_entry(roots / "source")
    replica = MolecularDatabase(roots / "replica")
    transport = Replica()
    _sync_all(transport, source, replica)
    peer_group_id = _id("peer-", SOURCE + ":" + group_id)
    entry_key = f"molecule_entry:{SOURCE}:{group_id}:{molecule_id}"
    snapshot = _entry_record(source, group_id, molecule_id)
    with replica.connection() as connection:
        # Local records in a replicated pair survive a source tombstone.
        connection.execute(
            "INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)",
            ("local-score", molecule_id, peer_group_id, "local", "local-run", "{}", "{}", "now"),
        )

    malformed = json.loads(json.dumps(snapshot))
    malformed["member"]["name"] = "should roll back"
    malformed["structures"][0]["context"] = {"not_finite": float("nan")}
    before = replica.list_molecules(OWNER, group_id=peer_group_id)["items"][0]
    with pytest.raises(ValueError):
        apply_record("molecule_entry", entry_key, malformed, SOURCE, database=replica)
    after = replica.list_molecules(OWNER, group_id=peer_group_id)["items"][0]
    assert after["name"] == before["name"]

    apply_record("molecule_entry", entry_key, None, SOURCE, database=replica)
    assert replica.list_molecules(OWNER, group_id=peer_group_id)["total"] == 0
    with replica.connection() as connection:
        for table in ("molecule_members", "molecule_structures", "molecule_sources"):
            assert connection.execute(f"SELECT count(*) FROM {table} WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE id='local-score'").fetchone()[0] == 1
        assert connection.execute(
            "SELECT count(*) FROM molecule_observations WHERE group_id=? AND molecule_id=? AND id!='local-score'",
            (peer_group_id, molecule_id),
        ).fetchone()[0] == 0

    # A source's tombstone cannot delete a group's row with a different author.
    other_group = "other-source-group"
    other_peer_group = _id("peer-", OTHER_SOURCE + ":" + other_group)
    group_value = {"node_id": OTHER_SOURCE, "group": {
        "id": other_group, "name": "Other", "owner_id": OUTSIDER["id"], "project_id": "", "import_method": "user_uploads",
        "origin": "{}", "schema_json": "{}", "created_at": "now", "updated_at": "now"}}
    apply_record("molecule_group", f"molecule_group:{OTHER_SOURCE}:{other_group}", group_value, OTHER_SOURCE, database=replica)
    with pytest.raises(PermissionError):
        apply_record("molecule_entry", entry_key, None, OTHER_SOURCE, database=replica)
    assert replica.group(other_peer_group, OUTSIDER)["name"] == "Other"


def test_same_identity_with_group_specific_chains_roles_and_order_replicates(roots):
    source = MolecularDatabase(roots / "source")
    target = MolecularDatabase(roots / "target")
    first = source.import_records(OWNER, [{
        "name": "first representation",
        "components": [
            {"chain": "A", "type": "protein", "sequence": "ACDE", "role": "binder"},
            {"chain": "B", "type": "protein", "sequence": "FGHI", "role": "target"},
        ],
    }], name="Source one")
    second = source.import_records(OWNER, [{
        "name": "second representation",
        "components": [
            {"chain": "H", "type": "protein", "sequence": "FGHI", "role": "target"},
            {"chain": "L", "type": "protein", "sequence": "ACDE", "role": "binder"},
        ],
    }], name="Source two")
    molecule_id = first["molecule_ids"][0]
    assert second["molecule_ids"] == [molecule_id]

    # A destination's existing copy may use its own chain names and order.
    target.import_records(OUTSIDER, [{
        "name": "local representation",
        "components": [
            {"chain": "X", "type": "protein", "sequence": "ACDE", "role": "binder"},
            {"chain": "Y", "type": "protein", "sequence": "FGHI", "role": "target"},
        ],
    }], name="Local copy")
    transport = Replica()
    _sync_all(transport, source, target)

    first_peer_group = _id("peer-", SOURCE + ":" + first["group_id"])
    second_peer_group = _id("peer-", SOURCE + ":" + second["group_id"])
    first_row = target.list_molecules(OWNER, group_id=first_peer_group)["items"][0]
    second_row = target.list_molecules(OWNER, group_id=second_peer_group)["items"][0]
    assert first_row["id"] == second_row["id"] == molecule_id
    assert first_row["components"][0]["chain"] == "A"
    assert second_row["components"][0]["chain"] == "H"
    assert [component["sequence"] for component in second_row["components"]] == ["FGHI", "ACDE"]
    assert target.list_molecules(OUTSIDER, group_id=target.groups(OUTSIDER)[0]["id"])["total"] == 1


def test_group_tombstone_removes_only_that_peer_group(roots):
    source, group_id, molecule_id, _ = _make_source_with_entry(roots / "source")
    replica = MolecularDatabase(roots / "replica")
    transport = Replica()
    _sync_all(transport, source, replica)
    peer_group_id = _id("peer-", SOURCE + ":" + group_id)
    with replica.connection() as connection:
        local_group = "local-group"
        connection.execute("INSERT INTO molecule_groups VALUES (?,?,?,?,?,?,?,?,?)",
                           (local_group, "Local", OWNER["id"], "", "user_uploads", "{}", "{}", "now", "now"))
        connection.execute("INSERT INTO molecule_members VALUES (?,?,?,?,?,?)",
                           (local_group, molecule_id, "Local copy", "[]", "{}", "now"))
    apply_record("molecule_group", f"molecule_group:{SOURCE}:{group_id}", None, SOURCE, database=replica)
    with replica.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_groups WHERE id=?", (peer_group_id,)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id=?", (local_group, molecule_id)).fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecules WHERE id=?", (molecule_id,)).fetchone()[0] == 1
