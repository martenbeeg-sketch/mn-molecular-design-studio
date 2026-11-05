"""Isolated regressions for the P0 promotion and curated-visibility fixes."""
from __future__ import annotations

import pytest

from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_visibility import (
    CONFIRMED_DEVELOPMENT_RUNS,
    UNREVIEWED_SCREENING_SHORTLISTS,
    is_unreviewed_screening_shortlist,
)


# Independent literals from the production origin rows. These must not be
# imported from the visibility predicate under test.
PRODUCTION_OWNER_ID = "00000000000000000000000000000001"
PRODUCTION_SOURCE_NODE = "pc002975"
PRODUCTION_SOURCE_NODE_ID = "50e484e4-0ad5-47c9-a8c4-b46d4000a59b"
ADMIN = {"id": PRODUCTION_OWNER_ID, "role": "admin", "_cluster_project_ids": []}
PRODUCTION_SHORTLISTS = (
    ("screening-final-a5ba730c962dc449fb40bbab7702a25b", "20261025-220608-e91da528"),
    ("screening-final-7483c3d8846ffad1263b3b8d78078d3d", "20261026-170220-29a4a53e"),
)


def _add_group(database, user, group_id, name, origin, smiles, *, file_name=None, file_bytes=None):
    database.create_group(
        user, name, import_method="tool_output", origin=origin, group_id=group_id,
    )
    record = {
        "name": name,
        "components": [{"type": "small_molecule", "chain": "L", "smiles": smiles}],
    }
    if file_name:
        record["files"] = [{"name": file_name, "content": file_bytes or b"fixture", "label": "Saved source"}]
    result = database.import_records(user, [record], group_id=group_id, import_method="tool_output")
    return result["molecule_ids"][0]


def test_historical_spiro_shortlists_are_library_ineligible_but_remain_addressable(tmp_path):
    database = MolecularDatabase(tmp_path / "database")
    shared_group_id = "screening-final-a5ba730c962dc449fb40bbab7702a25b"
    hidden_group_id = "screening-final-7483c3d8846ffad1263b3b8d78078d3d"
    shared_run_id = "20261025-220608-e91da528"
    hidden_run_id = "20261026-170220-29a4a53e"

    shared_molecule_id = _add_group(
        database, ADMIN, shared_group_id, "SpiroSpace shortlist 1",
        {"type": "job", "collection_kind": "curated", "run_id": shared_run_id,
         "source_node": PRODUCTION_SOURCE_NODE,
         "saved_selection": True, "final_screening_selection": True},
        "CCO", file_name="shortlist-1.sdf", file_bytes=b"original structure fixture",
    )
    hidden_molecule_id = _add_group(
        database, ADMIN, hidden_group_id, "SpiroSpace shortlist 2",
        {"type": "job", "collection_kind": "curated", "run_id": hidden_run_id,
         "source_node": PRODUCTION_SOURCE_NODE,
         "saved_selection": True, "final_screening_selection": True},
        "CCN", file_name="shortlist-2.sdf", file_bytes=b"original structure fixture 2",
    )
    observation_id = database.add_observation(
        shared_group_id, shared_molecule_id, "vina", shared_run_id,
        {"energy": -7.5}, context={"run_id": shared_run_id, "compound_id": "source-ethanol"},
    )
    independent_import = database.import_records(
        ADMIN,
        [{"name": "Intentional ethanol import", "components": [
            {"type": "small_molecule", "chain": "L", "smiles": "CCO"},
        ]}],
        name="Intentional reusable compounds",
        import_method="user_uploads",
        origin={"type": "user", "collection_kind": "curated", "source_id": "HY-L126-fixture"},
    )
    assert independent_import["molecule_ids"] == [shared_molecule_id]

    curated_groups = database.groups(ADMIN, collection="library")
    historical_groups = database.groups(ADMIN, collection="all")
    assert not ({shared_group_id, hidden_group_id} & {group["id"] for group in curated_groups})
    assert {shared_group_id, hidden_group_id} <= {group["id"] for group in historical_groups}
    # The ethanol identity remains in the curated view because it has an
    # independent intentional import; the other shortlist-only identity does not.
    library_ids = {row["id"] for row in database.list_molecules(ADMIN, collection="library")["items"]}
    assert shared_molecule_id in library_ids
    assert hidden_molecule_id not in library_ids
    visible_smiles = {
        part["smiles"] for row in database.list_molecules(ADMIN, collection="library")["items"]
        for part in row["components"] if part.get("smiles")
    }
    assert "CCO" in visible_smiles
    assert "CCN" not in visible_smiles

    # All indexed and direct-ID history preserve the original Evaluation,
    # source structure and both original memberships for review.
    hidden_history = database.list_molecules(
        ADMIN, molecule_ids=[hidden_molecule_id], collection="all",
    )["items"][0]
    shared_history = database.list_molecules(
        ADMIN, group_id=shared_group_id, collection="all",
    )["items"][0]
    assert hidden_history["id"] == hidden_molecule_id
    assert shared_history["id"] == shared_molecule_id
    assert observation_id in {row["id"] for row in shared_history["observations"]}
    assert any(row["group_id"] == shared_group_id for row in shared_history["representations"])
    assert database.group(hidden_group_id, ADMIN)["origin"]["run_id"] == hidden_run_id
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id IN (?,?)",
                                  (shared_group_id, hidden_group_id)).fetchone()[0] == 2
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE id=?",
                                  (observation_id,)).fetchone()[0] == 1


def test_hide_smoke_control_matches_only_confirmed_owner_run_pairs(tmp_path):
    database = MolecularDatabase(tmp_path / "database")

    confirmed_smiles = ["CCO", "CCN", "CCC", "CCCC", "CCCl", "CCBr", "CCF", "COC"]
    for index, (owner, run_id) in enumerate(sorted(CONFIRMED_DEVELOPMENT_RUNS)):
        group_id = f"confirmed-development-{index}"
        user = {"id": owner, "role": "admin"}
        _add_group(
            database, user, group_id, "Documented development validation",
            {"type": "job", "collection_kind": "generated", "run_id": run_id,
             "studio_development_mode": True},
            confirmed_smiles[index],
        )

    science = [
        ("PC004032", "20261019-181522-21f9a96c", "R3 smoke test campaign"),
        ("PC004033", "20261021-190947-7ea80083", "Tau test campaign"),
        ("PC004033", "20261029-094502-6ad69fc5", "PD-L1 smoke campaign"),
    ]
    for index, (owner, run_id, name) in enumerate(science):
        user = {"id": owner, "role": "admin"}
        _add_group(
            database, user, f"genuine-science-{index}", name,
            {"type": "job", "collection_kind": "generated", "run_id": run_id,
             "studio_development_mode": True},
            ["CC(=O)O", "c1ccccc1", "CC(=O)N"][index],
        )

    ambiguous_id = _add_group(
        database, ADMIN, "ambiguous-smoke-label", "Unreviewed test-name evidence",
        {"type": "job", "collection_kind": "generated", "run_id": "ambiguous-run",
         "is_test": True, "studio_development_mode": True},
        "C=O",
    )

    visible_with_control = {
        group["id"] for group in database.groups(ADMIN, collection="all", hide_smoke_test=True)
    }
    visible_historical = {
        group["id"] for group in database.groups(ADMIN, collection="all", hide_smoke_test=False)
    }
    confirmed_ids = {f"confirmed-development-{index}" for index in range(len(CONFIRMED_DEVELOPMENT_RUNS))}
    science_ids = {f"genuine-science-{index}" for index in range(len(science))}

    assert not (confirmed_ids & visible_with_control)
    assert confirmed_ids <= visible_historical
    assert science_ids <= visible_with_control
    assert "ambiguous-smoke-label" in visible_with_control
    assert ambiguous_id in {row["id"] for row in database.list_molecules(ADMIN, collection="all")["items"]}


def test_screening_collection_filter_requires_exact_audited_collection_owner_and_run(tmp_path):
    database = MolecularDatabase(tmp_path / "database")
    audited_id, identity = next(iter(UNREVIEWED_SCREENING_SHORTLISTS.items()))
    owner, run_id = identity["owner_id"], identity["run_id"]
    _add_group(
        database, ADMIN, audited_id, "Exact audited group",
        {"type": "job", "collection_kind": "curated", "run_id": run_id,
         "source_node": PRODUCTION_SOURCE_NODE}, "CO",
    )
    # Similar names and run IDs do not broaden the exact-ID exclusion.
    _add_group(
        database, ADMIN, "same-label-different-id", "SpiroSpace final shortlist",
        {"type": "job", "collection_kind": "curated", "run_id": run_id}, "CN",
    )
    shown = {group["id"] for group in database.groups(ADMIN, collection="library")}
    assert audited_id not in shown
    assert "same-label-different-id" in shown
    assert owner == PRODUCTION_OWNER_ID
    assert identity["source_node"] == PRODUCTION_SOURCE_NODE
    assert identity["source_node_id"] == PRODUCTION_SOURCE_NODE_ID
    assert not is_unreviewed_screening_shortlist({
        "id": audited_id, "owner_id": owner,
        "origin": {"run_id": "different-run", "source_node": PRODUCTION_SOURCE_NODE},
    })
    assert not is_unreviewed_screening_shortlist({
        "id": audited_id, "owner_id": "different-owner",
        "origin": {"run_id": run_id, "source_node": PRODUCTION_SOURCE_NODE},
    })
    assert not is_unreviewed_screening_shortlist({
        "id": audited_id, "owner_id": owner,
        "origin": {"run_id": run_id, "source_node": "PC004032"},
    })


@pytest.mark.parametrize("group_id,run_id", PRODUCTION_SHORTLISTS)
def test_production_shaped_source_and_peer_rows_match_owner_and_provenance(group_id, run_id):
    source = {
        "id": group_id,
        "owner_id": "00000000000000000000000000000001",
        "origin": {
            "type": "job", "task_group": "campaign", "source_node": "pc002975",
            "run_id": run_id, "saved_selection": True, "final_screening_selection": True,
        },
    }
    peer = {
        "id": "peer-production-shape",
        "owner_id": "00000000000000000000000000000001",
        "origin": {
            "type": "job", "task_group": "campaign", "source_node": "spoofed-source",
            "_replica_group": group_id,
            "_replica_node": "50e484e4-0ad5-47c9-a8c4-b46d4000a59b",
            "run_id": run_id, "saved_selection": True, "final_screening_selection": True,
        },
    }
    assert is_unreviewed_screening_shortlist(source)
    assert is_unreviewed_screening_shortlist(peer)
    assert not is_unreviewed_screening_shortlist({**source, "owner_id": "private-owner"})
    assert not is_unreviewed_screening_shortlist({**source, "origin": {**source["origin"], "run_id": "other-run"}})
    assert not is_unreviewed_screening_shortlist({
        **peer,
        "origin": {**peer["origin"], "_replica_node": "foreign-node", "source_node": "pc002975"},
    })


@pytest.mark.parametrize("source_group_id", sorted(UNREVIEWED_SCREENING_SHORTLISTS))
def test_audited_shortlist_visibility_survives_real_peer_replication(tmp_path, monkeypatch, source_group_id):
    import socket

    from mn_molecular_design_studio.molecular_database import _id
    from mn_molecular_design_studio.molecular_database_peer import apply_record
    from mn_molecular_design_studio.cluster.replication import Replica

    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))

    identity = UNREVIEWED_SCREENING_SHORTLISTS[source_group_id]
    owner, run_id = identity["owner_id"], identity["run_id"]
    assert owner == PRODUCTION_OWNER_ID
    assert identity["source_node"] == PRODUCTION_SOURCE_NODE
    assert identity["source_node_id"] == PRODUCTION_SOURCE_NODE_ID
    user = {"id": owner, "role": "admin", "_cluster_project_ids": []}
    source = MolecularDatabase(tmp_path / "source")
    replica = MolecularDatabase(tmp_path / "replica")
    replication_state = Replica()
    with replication_state.db() as connection:
        replication_state._write(connection, "sentinel:replicated-group", "molecule_group", {"owner": "fixture"},
                                 author="fixture-source-workstation")
        replication_state._write(connection, "sentinel:removed-group", "molecule_group", None,
                                 author="fixture-source-workstation")

    def replication_snapshot():
        with replication_state.db() as connection:
            return {
                table: [tuple(row) for row in connection.execute(f"SELECT * FROM {table} ORDER BY rowid")]
                for table in ("records", "state", "conflicts")
            }

    replication_before = replication_snapshot()
    molecule_id = _add_group(
        source, user, source_group_id, "Audited Spiro shortlist",
        {"type": "job", "collection_kind": "curated", "run_id": run_id,
         "source_node": PRODUCTION_SOURCE_NODE,
         "saved_selection": True, "final_screening_selection": True},
        "CCN",
    )
    # Keep a Results-page source association alongside the audited collection.
    # The Results reader is run-scoped and does not use the Library predicate.
    from mn_molecular_design_studio.molecular_database import _json
    with source.connection() as connection:
        connection.execute(
            "INSERT INTO molecule_sources VALUES (?,?,?,?)",
            (f"{socket.gethostname()}/campaign/{run_id}:historical-results", molecule_id,
             source_group_id, _json({"task_group": "campaign", "run_id": run_id,
                                     "compound_id": "historical-shortlist-record"})),
        )

    def snapshot(database):
        with database.connection() as connection:
            return {table: [tuple(row) for row in connection.execute(f"SELECT * FROM {table} ORDER BY rowid")]
                    for table in ("molecule_groups", "molecule_members", "molecule_sources",
                                  "molecule_structures", "molecule_observations")}

    source_before = snapshot(source)
    with source.connection() as connection:
        source_group = dict(connection.execute(
            "SELECT * FROM molecule_groups WHERE id=?", (source_group_id,),
        ).fetchone())
        molecule = dict(connection.execute(
            "SELECT * FROM molecules WHERE id=?", (molecule_id,),
        ).fetchone())
        member = dict(connection.execute(
            "SELECT * FROM molecule_members WHERE group_id=? AND molecule_id=?",
            (source_group_id, molecule_id),
        ).fetchone())
        observations = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_observations WHERE group_id=?", (source_group_id,),
        )]
        sources = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_sources WHERE group_id=?", (source_group_id,),
        )]

    author = PRODUCTION_SOURCE_NODE_ID
    apply_record("molecule_group", f"molecule_group:{author}:{source_group_id}",
                 {"node_id": author, "group": source_group}, author, database=replica)
    apply_record("molecule_entry", f"molecule_entry:{author}:{source_group_id}:{molecule_id}",
                 {"node_id": author, "group_id": source_group_id, "molecule": molecule,
                  "member": member, "structures": [], "observations": observations,
                  "sources": sources}, author, database=replica)
    replica_group_id = _id("peer-", author + ":" + source_group_id)
    assert snapshot(source) == source_before
    assert replication_snapshot() == replication_before

    # The source and replica use the same audited source identity and owner/run
    # restrictions; the peer-local primary key is not mistaken for another run.
    assert is_unreviewed_screening_shortlist(source.groups(user, collection="all")[0])
    assert is_unreviewed_screening_shortlist(replica.group(replica_group_id, user))
    assert not source.groups(user, collection="library")
    assert not replica.groups(user, collection="library")
    source_history = source.list_molecules(user, molecule_ids=[molecule_id], collection="all")["items"][0]
    assert source_history["id"] == molecule_id
    assert replica.group(replica_group_id, user)["origin"]["_replica_group"] == source_group_id
    assert replica.group(replica_group_id, user)["origin"]["_replica_node"] == author
    assert replica.group(replica_group_id, user)["owner_id"] == owner
    assert replica.group(replica_group_id, user)["origin"]["run_id"] == run_id
    assert replica.list_molecules(user, collection="library")["total"] == 0

    # Exact group and run collisions from a foreign workstation stay eligible.
    foreign_peer = {
        "id": "peer-foreign-collision", "owner_id": owner,
        "origin": {"_replica_group": source_group_id, "_replica_node": "foreign-node",
                   "source_node": PRODUCTION_SOURCE_NODE, "run_id": run_id},
    }
    assert not is_unreviewed_screening_shortlist(foreign_peer)

    # Independent intentional membership restores this molecular identity to
    # Library without making the unreviewed peer collection eligible.
    independent = replica.import_records(
        user, [{"name": "Intentional local selection", "components": [
            {"type": "small_molecule", "chain": "L", "smiles": "CCN"},
        ]}], name="Intentional local selection", import_method="user_uploads",
    )
    assert independent["molecule_ids"] == [molecule_id]
    library_ids = {row["id"] for row in replica.list_molecules(user, collection="library")["items"]}
    assert molecule_id in library_ids
    assert replica_group_id not in {row["id"] for row in replica.groups(user, collection="library")}

    historical = replica.list_molecules(user, molecule_ids=[molecule_id], collection="all")["items"][0]
    assert historical["id"] == molecule_id
    assert replica_group_id in {group["id"] for group in historical["groups"]}
    assert replica.list_molecules(user, group_id=replica_group_id, collection="all")["items"][0]["id"] == molecule_id

    from mn_molecular_design_studio import molecular_database_results
    run_dir = tmp_path / "campaign" / run_id
    monkeypatch.setattr(molecular_database_results, "_authorized_run",
                        lambda _selected, _user: (run_dir, "campaign", run_id, owner, ""))
    result_rows = molecular_database_results.list_run_molecules(run_dir, user, database=source)
    assert molecule_id in {row["id"] for row in result_rows}

    stranger = {"id": "unrelated-user", "role": "user", "_cluster_project_ids": []}
    assert replica_group_id not in {row["id"] for row in replica.groups(stranger, collection="all")}
    with pytest.raises(PermissionError):
        replica.group(replica_group_id, stranger)

    replica_before = snapshot(replica)
    source.groups(user, collection="library")
    replica.groups(user, collection="library")
    replica.list_molecules(user, collection="all")
    assert snapshot(source) == source_before
    assert snapshot(replica) == replica_before
    assert replication_snapshot() == replication_before
