from __future__ import annotations

import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_results import remove_run_records


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}


def test_run_cleanup_removes_automatic_records_but_keeps_saved_molecules(tmp_path):
    database = MolecularDatabase(tmp_path / "database")
    task_group, run_id = "vina", "run-123"
    run_origin = {
        "type": "job", "task_group": task_group, "run_id": run_id,
        "collection_kind": "generated", "tool": "vina",
    }
    automatic_group = database.create_group(
        USER, "Run output", import_method="tool_output", origin=run_origin,
        group_id="job-output-run-123",
    )
    automatic = database.import_records(
        USER, molecule_imports.parse_smiles("CCO ethanol"),
        group_id=automatic_group, import_method="tool_output", origin=run_origin,
    )
    molecule_id = automatic["molecule_ids"][0]
    database.add_observation(
        automatic_group, molecule_id, "vina", run_id, {"affinity": -7.2},
        context={"task_group": task_group, "run_id": run_id},
    )
    saved_group = database.create_group(USER, "Shortlist", origin={"type": "user"})
    saved_origin = {"type": "job", "task_group": task_group, "run_id": run_id, "saved_selection": True}
    database.import_records(
        USER, molecule_imports.parse_smiles("OCC ethanol"), group_id=saved_group,
        origin=saved_origin,
    )

    result = remove_run_records(
        task_group, run_id, USER, owner_id=USER["id"], project_id="", database=database,
    )

    assert result["removed_sources"] == 1
    assert result["removed_groups"] == 1
    assert result["preserved_saved_selections"]
    assert database.list_molecules(USER, group_id=automatic_group)["total"] == 0
    saved = database.list_molecules(USER, group_id=saved_group)["items"]
    assert [row["id"] for row in saved] == [molecule_id]
    assert not saved[0]["observations"]
    with database.connection() as connection:
        remaining_sources = connection.execute(
            "SELECT count(*) FROM molecule_sources WHERE group_id=?", (saved_group,)
        ).fetchone()[0]
    assert remaining_sources == 1


def test_run_cleanup_requires_job_access(tmp_path):
    database = MolecularDatabase(tmp_path / "database")
    other_user = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}
    with pytest.raises(PermissionError):
        remove_run_records("vina", "run-123", other_user, owner_id=USER["id"], database=database)
