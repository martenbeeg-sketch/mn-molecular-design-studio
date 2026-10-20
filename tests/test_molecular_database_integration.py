from __future__ import annotations
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import ligand_compounds, molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_ingest import ingest_run, migrate_legacy_groups
from mn_molecular_design_studio.molecular_database_peer import apply_record

A = {"id": "a" * 32, "role": "user", "_cluster_project_ids": ["project-1"]}
B = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    # The API import initializes the full Studio workspace, whose canonical
    # root is MN_MOLECULAR_DESIGN_STUDIO_HOME. Isolate that as well as the
    # database runtime path so API middleware cannot reuse a developer library.
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    return MolecularDatabase()


def test_legacy_subset_migrates_exact_members_and_remains_a_tool_subset(isolated):
    path = isolated.root / "compound_sets/set1.json"
    ligand_compounds.write_compound_set(path, name="Original", compounds=[{"compound_id": "chosen", "smiles": "CCO"}, {"compound_id": "excluded", "smiles": "CC"}], owner_id=A["id"], project_id="project-1", created_at="now", compound_set_id="set1")
    (isolated.root / "groups.json").write_text(json.dumps([{"id": "subset", "name": "Only chosen", "owner_id": A["id"], "project_id": "project-1", "members": [{"record_type": "compound", "compound_set_id": "set1", "compound_id": "chosen"}]}]))
    migrate_legacy_groups(isolated, A)
    option = isolated.compound_options(A)[0]
    assert option["compound_count"] == 1
    assert [row["compound_id"] for row in json.loads(Path(option["path"]).read_text())["compounds"]] == ["chosen"]
    migrate_legacy_groups(isolated, A)
    assert len(isolated.groups(A)) == 1


def test_docking_result_index_is_idempotent_and_links_existing_compound(isolated, tmp_path):
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID
    original = isolated.import_records(A, molecule_imports.parse_smiles("CCO ethanol"), name="Library", project_id="project-1")
    run = tmp_path / "runs/ligand-docking/run-1"
    (run / "artifacts/input").mkdir(parents=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({"application": STUDIO_SOURCE_ID, "project_id": "project-1", "studio_owner_id": A["id"]}))
    (run / "metadata.json").write_text(json.dumps({"tool": "vina", "status": "completed", "campaign_name": "Docking", "completed_at": "2026-10-20T12:00:00Z"}))
    (run / "input.json").write_text(json.dumps({"tool": "vina", "params": {"exhaustiveness": 32}}))
    (run / "result.json").write_text(json.dumps({"success": True, "metrics": {"primary_ranking_metric": "best_score_kcal_mol"}}))
    ligand_compounds.write_compound_set(run / "artifacts/input/compound_set.json", name="Consumed", compounds=[{"compound_id": "ethanol", "smiles": "OCC"}], owner_id=A["id"], project_id="project-1", created_at="now")
    (run / "artifacts/docking_scores.csv").write_text("compound_id\tstatus\tbest_score_kcal_mol\tpose_path\nethanol\tcompleted\t-7.5\t\n")
    assert ingest_run(run, database=isolated) == 1
    assert ingest_run(run, database=isolated) == 0
    row = isolated.list_molecules(A, group_id=original["group_id"])["items"][0]
    assert row["id"] == original["molecule_ids"][0]
    assert row["scores"]["vina.best_score_kcal_mol"] == -7.5
    assert row["observations"][0]["context"]["settings"]["exhaustiveness"] == 32
    assert len(row["observations"]) == 1


def test_database_ui_import_wizard_and_score_columns(isolated, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import database_ui
    monkeypatch.setattr(database_ui, "sync_workspace", lambda *args, **kwargs: [])
    script = f'''\
from mn_molecular_design_studio.database_ui import render_database
render_database({A!r}, [])
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    app.button(key="molecular-database-new").click().run()
    app.text_input(key="molecular-import-group-name").set_value("Imported proteins").run()
    app.button(key="molecular-import-next").click().run()
    assert not app.exception
    app.radio(key="molecular-import-method").set_value("Paste sequences").run()
    app.text_area(key="molecular-import-sequences").set_value(">first\nACDEFG\n>second\nACDEFA").run()
    assert not app.button(key="molecular-import-commit").disabled
    app.button(key="molecular-import-commit").click().run()
    assert not app.exception
    assert isolated.list_molecules(A)["total"] == 2
    assert any("Imported 2 rows" in element.value for element in app.success)


def test_api_import_search_export_and_access(isolated, monkeypatch):
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio.api import app, _current_user
    from mn_molecular_design_studio import molecular_database_api
    monkeypatch.setattr(molecular_database_api, "sync_workspace", lambda *args, **kwargs: [])
    app.dependency_overrides[_current_user] = lambda: A
    try:
        with TestClient(app) as client:
            imported = client.post("/api/v1/molecules/imports", json={"method": "smiles", "text": "CCO selected\nCC excluded", "project_id": "project-1"})
            assert imported.status_code == 201, imported.text
            ids = imported.json()["molecule_ids"]
            assert client.get("/api/v1/molecules").json()["total"] == 2
            selection = client.post("/api/v1/molecules/selections", json={"molecule_ids": [ids[0]], "project_id": "project-1"})
            assert selection.status_code == 200, selection.text
            assert len(json.loads(Path(selection.json()["path"]).read_text())["compounds"]) == 1
            private = isolated.import_records(B, molecule_imports.parse_smiles("CCN private"))
            assert client.get("/api/v1/molecules/" + private["molecule_ids"][0]).status_code == 404
            rejected = client.post("/api/v1/molecules/imports", json={"method": "records", "records": [{"name": "read local", "components": [{"type": "protein", "chain": "A", "sequence": "ACD"}], "files": [{"path": "/etc/passwd"}]}]})
            assert rejected.status_code == 201
            molecule = client.get("/api/v1/molecules/" + rejected.json()["molecule_ids"][0]).json()
            assert molecule["structures"] == []
    finally:
        app.dependency_overrides.pop(_current_user, None)


def test_peer_metadata_is_read_only_and_structure_access_follows_group_acl(isolated, tmp_path):
    import hashlib
    from mn_molecular_design_studio.molecular_database import _id

    source = MolecularDatabase(tmp_path / "source")
    content = b"HEADER    protein\nEND\n"
    imported = source.import_records(A, [{"name": "shared", "components": [{"chain": "A", "type": "protein", "sequence": "ACDE"}],
                                         "files": [{"name": "model.pdb", "content": content}]}], project_id="project-1")
    with source.connection() as connection:
        group = dict(connection.execute("SELECT * FROM molecule_groups WHERE id=?", (imported["group_id"],)).fetchone())
        molecule_id = imported["molecule_ids"][0]
        molecule = dict(connection.execute("SELECT * FROM molecules WHERE id=?", (molecule_id,)).fetchone())
        member = dict(connection.execute("SELECT * FROM molecule_members WHERE group_id=? AND molecule_id=?", (imported["group_id"], molecule_id)).fetchone())
        structure = dict(connection.execute("SELECT * FROM molecule_structures WHERE group_id=? AND molecule_id=?", (imported["group_id"], molecule_id)).fetchone())
    node_id = "node-1"
    apply_record("molecule_group", f"molecule_group:{node_id}:{group['id']}", {"node_id": node_id, "group": group}, node_id, database=isolated)
    entry = {"node_id": node_id, "group_id": group["id"], "molecule": molecule, "member": member,
             "structures": [{**structure, "path": "", "sha256": hashlib.sha256(content).hexdigest()}],
             "observations": [], "sources": []}
    apply_record("molecule_entry", f"molecule_entry:{node_id}:{group['id']}:{molecule_id}", entry, node_id, database=isolated)

    replica_group_id = _id("peer-", node_id + ":" + group["id"])
    replicated = isolated.group(replica_group_id, A)
    assert replicated["read_only"] is True
    assert isolated.list_molecules(A, group_id=replica_group_id)["total"] == 1
    assert isolated.list_molecules(B, group_id=replica_group_id)["total"] == 0
    with pytest.raises(PermissionError):
        isolated.group(replica_group_id, A, write=True)
    with isolated.connection() as connection:
        remote_structure_id = connection.execute("SELECT id FROM molecule_structures WHERE group_id=?", (replica_group_id,)).fetchone()[0]
    with pytest.raises(PermissionError):
        isolated.structure_file(B, remote_structure_id)


def test_protein_pipeline_source_freezes_exact_database_selection(isolated):
    from mn_molecular_design_studio.database_selections import freeze_protein_selection
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.workspace import create_project

    project_id = create_project("Molecular selection test", home=isolated.root.parent, owner_id=A["id"])["id"]
    project_user = {**A, "_cluster_project_ids": [project_id]}
    imported = isolated.import_records(project_user, [
        {"name": "lead-target", "components": [{"chain": "A", "type": "protein", "sequence": "ACDE", "role": "lead"},
                                                   {"chain": "B", "type": "protein", "sequence": "FGHI", "role": "target"}]},
        {"name": "excluded", "components": [{"chain": "A", "type": "protein", "sequence": "KLMN"}]},
    ], name="Candidates", project_id=project_id)
    selected = imported["molecule_ids"][0]
    run_dir = freeze_protein_selection(isolated, project_user, [selected], project_id, name="Chosen candidate")
    snapshot = json.loads((run_dir / "selection_snapshot.json").read_text())
    candidates = read_candidates(run_dir)
    assert snapshot["molecule_ids"] == [selected]
    assert [row["id"] for row in snapshot["molecules"]] == [selected]
    assert candidates[0]["candidate_id"] == selected
    assert candidates[0]["binder_sequence"] == "ACDE"
    assert candidates[0]["target_chains"] == ["B"]
