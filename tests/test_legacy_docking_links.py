from __future__ import annotations

import csv
import hashlib
import json

import pytest

from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
from mn_molecular_design_studio.legacy_docking_links import apply_docking_links, plan_docking_links
from mn_molecular_design_studio.molecule_imports import parse_sequences, parse_smiles
from test_molecular_database import ADMIN, database


def _original(database, tmp_path):
    target = database.import_records(ADMIN, parse_sequences("ACDE"), name="Targets")
    compound = database.import_records(ADMIN, parse_smiles("CCO ethanol"), name="Imported compounds")
    directory = tmp_path / "originals/docking/original-attempt"
    (directory / "input").mkdir(parents=True)
    receptor = directory / "input/receptor.pdb"
    residues = ("ALA", "CYS", "ASP", "GLU")
    receptor.write_text("\n".join(f"ATOM  {number:5d}  CA  {residue} A{number:4d}    "
        f"{number * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        for number, residue in enumerate(residues, 1)) + "\nEND\n")
    digest = hashlib.sha256(receptor.read_bytes()).hexdigest()
    pose = directory / "results/replicate_001/CMP1_out.pdbqt"
    pose.parent.mkdir(parents=True)
    pose.write_text("MODEL 1\nREMARK VINA RESULT: -7.100 0.000 0.000\n"
        "ATOM      1  C1  UNL L   1       1.000   2.000   3.000  0.00  0.00     0.000 C\nENDMDL\n")
    (directory / "metadata.json").write_text(json.dumps({"run_id": directory.name, "engine": "vina", "status": "completed"}))
    (directory / "input.json").write_text(json.dumps({"target_artifact": {"role": "receptor", "sha256": digest}}))
    source_row = {"compound_id": "CMP1", "engine": "vina", "replicate": "1", "seed": "1001",
        "best_score_kcal_mol": "-7.1", "pose_file": "results/replicate_001/CMP1_out.pdbqt"}
    with (directory / "scores.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(source_row))
        writer.writeheader()
        writer.writerow(source_row)
    context = {"task_group": "docking", "run_id": "mnligand-imported", "tool": "vina", "is_representative": True,
        "legacy_run_id": directory.name, "target_molecule_id": target["molecule_ids"][0],
        "target_structure_digest": digest, "target_variant_sha256": digest,
        "representation_id": "missing-old-pose", "structure_digest": "f" * 64,
        "row": {**source_row, "legacy_run_id": directory.name,
                "pose_path": "artifacts/replicates/original-attempt/" + source_row["pose_file"]}}
    observation = database.add_observation(compound["group_id"], compound["molecule_ids"][0], "vina", "mnligand-imported",
        {"replicate": 1, "seed": 1001, "best_score_kcal_mol": -7.1}, context=context)
    return directory, observation, pose


def test_original_score_and_pose_are_reconnected_with_exact_receptor_without_new_scores(database, tmp_path):
    directory, observation, pose = _original(database, tmp_path)
    before = database.path.stat()
    readonly = ReadOnlyLibrary(database.root)
    plan = plan_docking_links(readonly, ADMIN, directory.parent.parent)
    assert (before.st_size, before.st_mtime_ns) == (database.path.stat().st_size, database.path.stat().st_mtime_ns)
    assert plan["unresolved"] == []
    assert len(plan["observations"]) == 1
    assert len(plan["structures"]) == 2
    structure_context = next(json.loads(row["context"]) for row in plan["structures"]
                             if json.loads(row["context"]).get("database_structure_policy") == "best_docking_pose")
    assert "snapshots" not in structure_context["docking_link_reconciliation"]
    assert structure_context["docking_link_reconciliation"]["evidence_sha256"]
    backup = tmp_path / "recovery.sqlite3"
    result = apply_docking_links(database, ADMIN, plan, backup_path=backup)
    assert result["observations"] == 1 and backup.is_file()
    with database.connection() as connection:
        rows = connection.execute("SELECT * FROM molecule_observations").fetchall()
        assert len(rows) == 1
        saved = dict(rows[0])
        assert saved["id"] == observation
        assert saved["metrics"] == plan["observations"][0]["before"]["metrics"]
        context = json.loads(saved["context"])
        assert context["structure_digest"] == hashlib.sha256(pose.read_bytes()).hexdigest()
        assert context["docking_link_reconciliation"]["previous"]["representation_id"] == "missing-old-pose"
        for key in ("representation_id", "target_structure_id"):
            structure = connection.execute("SELECT * FROM molecule_structures WHERE id=?", (context[key],)).fetchone()
            assert structure and database.structure_file(ADMIN, structure["id"]).is_file()
    repeated = plan_docking_links(readonly, ADMIN, directory.parent.parent)
    assert repeated["observations"] == [] and repeated["unresolved"] == []


@pytest.mark.parametrize("field,value", [("replicate", "2"), ("engine", "gnina"),
    ("best_score_kcal_mol", "-9.5"), ("pose_path", "artifacts/replicates/other/results/replicate_001/CMP1_out.pdbqt"),
    ("compound_id", "another-compound")])
def test_conflicting_source_evidence_is_not_repaired(database, tmp_path, field, value):
    directory, observation, _ = _original(database, tmp_path)
    with database.connection() as connection:
        context = json.loads(connection.execute("SELECT context FROM molecule_observations WHERE id=?", (observation,)).fetchone()[0])
        context["row"][field] = value
        connection.execute("UPDATE molecule_observations SET context=? WHERE id=?", (json.dumps(context), observation))
    plan = plan_docking_links(ReadOnlyLibrary(database.root), ADMIN, directory.parent.parent)
    assert plan["observations"] == [] and len(plan["unresolved"]) == 1


@pytest.mark.parametrize("change", ["source", "score", "membership", "peer", "plan"])
def test_apply_rejects_changed_sources_concurrent_edits_and_peer_ownership(database, tmp_path, change):
    directory, observation, pose = _original(database, tmp_path)
    plan = plan_docking_links(ReadOnlyLibrary(database.root), ADMIN, directory.parent.parent)
    if change == "source":
        pose.write_text(pose.read_text() + "REMARK changed\n")
    elif change == "plan":
        plan["observations"][0]["after_context"] = "{}"
    else:
        with database.connection() as connection:
            if change == "score":
                connection.execute("UPDATE molecule_observations SET metrics='{}' WHERE id=?", (observation,))
            elif change == "membership":
                connection.execute("UPDATE molecule_members SET name='changed'")
            elif change == "peer":
                connection.execute("UPDATE molecule_groups SET origin=?", (json.dumps({"_replica_node": "other-node"}),))
    with pytest.raises(ValueError):
        apply_docking_links(database, ADMIN, plan, backup_path=tmp_path / "recovery.sqlite3")
    with database.connection() as connection:
        context = json.loads(connection.execute("SELECT context FROM molecule_observations WHERE id=?", (observation,)).fetchone()[0])
        assert context["representation_id"] == "missing-old-pose"
        assert connection.execute("SELECT count(*) FROM molecule_structures").fetchone()[0] == 0


def test_receptor_mismatch_and_escaping_pose_are_rejected(database, tmp_path):
    directory, _, _ = _original(database, tmp_path)
    (directory / "input/receptor.pdb").write_text("HEADER wrong frame\n")
    plan = plan_docking_links(ReadOnlyLibrary(database.root), ADMIN, directory.parent.parent)
    assert plan["observations"] == [] and "receptor frame" in plan["unresolved"][0]["reason"]


def test_same_coordinate_bytes_do_not_merge_separate_attempts_or_replicates(database, tmp_path):
    directory, observation, _ = _original(database, tmp_path)
    imported = database.import_records(ADMIN, parse_sequences("ACDE"), name="Model evidence")
    receptor = directory / "input/receptor.pdb"
    item = {"path": str(receptor)}
    base = {"run_id": "same-import", "tool": "vina", "legacy_run_id": "attempt-1", "replicate": "1"}
    first = database.attach_structure(imported["group_id"], imported["molecule_ids"][0], item, context=base)
    second = database.attach_structure(imported["group_id"], imported["molecule_ids"][0], item,
        context={**base, "legacy_run_id": "attempt-2"})
    third = database.attach_structure(imported["group_id"], imported["molecule_ids"][0], item,
        context={**base, "replicate": "2"})
    assert len({first, second, third}) == 3
    unrelated = database.add_observation(imported["group_id"], imported["molecule_ids"][0], "vina", "same-import", {"score": -7},
        context={"legacy_run_id": "attempt-3", "replicate": "1"})
    with database.connection() as connection:
        context = json.loads(connection.execute("SELECT context FROM molecule_observations WHERE id=?", (unrelated,)).fetchone()[0])
    assert not context.get("representation_id")


def test_identical_replicate_bytes_do_not_infer_a_different_replicate(database, tmp_path):
    directory, _, _ = _original(database, tmp_path)
    imported = database.import_records(ADMIN, parse_sequences("ACDE"), name="Representations")
    database.attach_structure(imported["group_id"], imported["molecule_ids"][0], {"path": str(directory / "input/receptor.pdb")},
        context={"run_id": "same-import", "legacy_run_id": "attempt-1", "row": {"replicate": "1", "compound_id": "CMP1"}})
    observation = database.add_observation(imported["group_id"], imported["molecule_ids"][0], "vina", "same-import", {"score": -7},
        context={"legacy_run_id": "attempt-1", "row": {"replicate": "2", "compound_id": "CMP1"}})
    with database.connection() as connection:
        saved = json.loads(connection.execute("SELECT context FROM molecule_observations WHERE id=?", (observation,)).fetchone()[0])
    assert not saved.get("representation_id")


def test_plan_skips_foreign_owned_groups(database, tmp_path):
    directory, _, _ = _original(database, tmp_path)
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET origin=?", (json.dumps({"_replica_node": "owner-node"}),))
    plan = plan_docking_links(ReadOnlyLibrary(database.root), ADMIN, directory.parent.parent)
    assert plan["examined"] == 0 and plan["observations"] == []


def test_source_pose_symlink_cannot_escape_selected_original_run(database, tmp_path):
    directory, _, pose = _original(database, tmp_path)
    outside = tmp_path / "other-pose.pdbqt"
    outside.write_bytes(pose.read_bytes())
    pose.unlink()
    pose.symlink_to(outside)
    plan = plan_docking_links(ReadOnlyLibrary(database.root), ADMIN, directory.parent.parent)
    assert not plan["observations"] and "escapes" in plan["unresolved"][0]["reason"]


def test_verified_pose_and_receptor_links_survive_peer_replay(database, tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster.replication import Replica
    from mn_molecular_design_studio.evaluations import evaluation_descriptor
    from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id
    from test_molecular_database_peer import SOURCE, _sync_all

    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    directory, observation, _ = _original(database, tmp_path)
    peer = MolecularDatabase(tmp_path / "peer-database")
    replica = Replica()
    _sync_all(replica, database, peer)
    plan = plan_docking_links(ReadOnlyLibrary(database.root), ADMIN, directory.parent.parent)
    apply_docking_links(database, ADMIN, plan, backup_path=tmp_path / "before.sqlite3")
    _sync_all(replica, database, peer)
    _sync_all(replica, database, peer)
    with database.connection() as source, peer.connection() as target:
        original = dict(source.execute("SELECT * FROM molecule_observations WHERE id=?", (observation,)).fetchone())
        replayed = dict(target.execute("SELECT * FROM molecule_observations WHERE id=?", (_id("peer-score-", SOURCE + ":" + observation),)).fetchone())
        assert original["metrics"] == replayed["metrics"]
        assert evaluation_descriptor(original)["id"] == evaluation_descriptor(replayed)["id"]
        context = json.loads(replayed["context"])
        assert target.execute("SELECT 1 FROM molecule_structures WHERE id=?", (context["representation_id"],)).fetchone()
        assert target.execute("SELECT 1 FROM molecule_structures WHERE id=?", (context["target_structure_id"],)).fetchone()
        pose = target.execute("SELECT context FROM molecule_structures WHERE id=?", (context["representation_id"],)).fetchone()
        assert json.loads(pose[0])["target_structure_id"] == context["target_structure_id"]
