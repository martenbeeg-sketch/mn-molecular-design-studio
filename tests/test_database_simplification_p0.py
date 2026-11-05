"""Isolated regressions for the P0 promotion and curated-visibility fixes."""
from __future__ import annotations

from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_visibility import (
    CONFIRMED_DEVELOPMENT_RUNS,
    UNREVIEWED_SCREENING_SHORTLISTS,
    is_unreviewed_screening_shortlist,
)


ADMIN = {"id": "pc002975", "role": "admin", "_cluster_project_ids": []}


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
         "saved_selection": True, "final_screening_selection": True},
        "CCO", file_name="shortlist-1.sdf", file_bytes=b"original structure fixture",
    )
    hidden_molecule_id = _add_group(
        database, ADMIN, hidden_group_id, "SpiroSpace shortlist 2",
        {"type": "job", "collection_kind": "curated", "run_id": hidden_run_id,
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
    audited_id, (owner, run_id) = next(iter(UNREVIEWED_SCREENING_SHORTLISTS.items()))
    _add_group(
        database, ADMIN, audited_id, "Exact audited group",
        {"type": "job", "collection_kind": "curated", "run_id": run_id}, "CO",
    )
    # Similar names and run IDs do not broaden the exact-ID exclusion.
    _add_group(
        database, ADMIN, "same-label-different-id", "SpiroSpace final shortlist",
        {"type": "job", "collection_kind": "curated", "run_id": run_id}, "CN",
    )
    shown = {group["id"] for group in database.groups(ADMIN, collection="library")}
    assert audited_id not in shown
    assert "same-label-different-id" in shown
    assert owner == "pc002975"
    assert not is_unreviewed_screening_shortlist({
        "id": audited_id, "owner_id": owner,
        "origin": {"run_id": "different-run"},
    })
    assert not is_unreviewed_screening_shortlist({
        "id": audited_id, "owner_id": "different-owner",
        "origin": {"run_id": run_id},
    })
