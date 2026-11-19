from __future__ import annotations

import hashlib
import json
import sqlite3

import pytest

from mn_molecular_design_studio.library_curation import apply_curation, encoded, plan_curation
from mn_molecular_design_studio.scientific_naming import design_name, meaningful_name
from mn_molecular_design_studio.target_registration import contents_identity, target_components
from test_database_browser import USER, _import, _record, database


ADMIN = {**USER, "role": "admin"}


def creation(sequence="ACDE"):
    from mn_molecular_design_studio.design_creation import sequence_digest
    return {"schema": "studio.design-creation.v1", "at": "2026-11-01T12:34:56+00:00", "precision": "second",
        "basis": "native_cif_creation", "workstation": "pc004033", "source_sha256": "a" * 64,
        "binder_sequence_sha256": sequence_digest(sequence)}


def legacy_design(database):
    target = [{"type": "protein", "chain": "A", "sequence": "EFG", "role": "target"}]
    reference = _import(database, name="References", run_id="", method="user_uploads", kind="reference",
        records=[_record(name="Test target", components=target)])
    generated = _import(database, name="Design run", run_id="20261101-120000-aa11", origin={"tool": "rfd4_proteina"},
        records=[_record(name="rfd4_000001_binder_rep3", metadata={"design_creation": creation("ACD")}, components=[*target,
            {"type": "protein", "chain": "B", "sequence": "ACD", "role": "binder"}])])
    observation = database.add_observation(generated["group_id"], generated["molecule_ids"][0], "rfd4_proteina", "20261101-120000-aa11", {"score": 0.8},
        context={"target_id": "legacy-target", "candidate": {"target_chains": ["A"], "binder_chains": ["B"]}})
    return reference, generated, observation


def test_naming_uses_recorded_creation_time_and_workstation_not_run_start():
    label = design_name("cs_" + "a" * 64, target_name="PD-L1 · 115 aa", tool="rfd4_proteina", run_id="20261101-120000-aa11",
        metadata={"native_candidate_id": "rfd4_000024_protein_binder_rep3", "distributed_campaign": {"batch_id": "b000004"},
            "design_creation": creation()})
    assert label == "Design_pc004033_20261101T123456Z"
    assert not meaningful_name("target-5e949950d8255ec3")
    assert not design_name("unknown", target_name="Target", tool="generator", run_id="run")


def test_future_generated_registration_names_automatically_without_changing_identity(database):
    target = _import(database, run_id="", name="Targets", kind="reference", method="user_uploads",
        records=[_record("EFG", name="Test target", components=[{"type": "protein", "chain": "A", "sequence": "EFG", "role": "target"}])])
    record = _record(name="rfd4_000001_binder_rep3", metadata={"design_target_molecule_id": target["molecule_ids"][0],
        "design_creation": creation()})
    imported = _import(database, records=[record], run_id="20261101-120000-aa11", origin={"tool": "rfd4_proteina"})
    repeated = _import(database, records=[record], run_id="20261101-120001-aa22", origin={"tool": "rfd4_proteina"})
    assert imported["molecule_ids"] == repeated["molecule_ids"]
    with database.connection() as connection:
        rows = connection.execute("SELECT name,metadata FROM molecule_members WHERE molecule_id=?", (imported["molecule_ids"][0],)).fetchall()
    assert len({row["name"] for row in rows}) == 1
    assert all(row["name"] == "Design_pc004033_20261101T123456Z" for row in rows)
    assert all(json.loads(row["metadata"])["studio_naming"]["native_names"] == [record["name"]] for row in rows)


def test_curation_is_backed_up_idempotent_and_keeps_native_evidence(database, tmp_path):
    reference, generated, observation = legacy_design(database)
    plan = plan_curation(database, ADMIN)
    assert not plan["unresolved"]
    assert len(plan["members"]) == 1 and len(plan["observations"]) == 1
    backup = tmp_path / "recovery/library.sqlite3"
    report = apply_curation(database, ADMIN, plan, backup_path=backup)
    assert report["members"] == 1 and backup.is_file()
    with database.connection() as connection:
        row = connection.execute("SELECT * FROM molecule_observations WHERE id=?", (observation,)).fetchone()
        assert json.loads(row["metrics"]) == {"score": 0.8}
        context = json.loads(row["context"])
        assert context["target_id"] == "legacy-target"
        assert context["target_molecule_id"] == reference["molecule_ids"][0]
        assert json.loads(connection.execute("SELECT metadata FROM molecule_members WHERE group_id=? AND molecule_id=?",
            (generated["group_id"], generated["molecule_ids"][0])).fetchone()[0])["studio_naming"]["native_names"] == ["rfd4_000001_binder_rep3"]
    assert not plan_curation(database, ADMIN)["members"]
    assert not plan_curation(database, ADMIN)["observations"]


def test_curation_rejects_stale_plan_without_partial_changes(database, tmp_path):
    reference, generated, observation = legacy_design(database)
    plan = plan_curation(database, ADMIN)
    database.update_metadata(ADMIN, generated["group_id"], generated["molecule_ids"][0], {}, name="Instructor's chosen name")
    with pytest.raises(ValueError, match="changed"):
        apply_curation(database, ADMIN, plan, backup_path=tmp_path / "backup.sqlite3")
    with database.connection() as connection:
        assert not json.loads(connection.execute("SELECT context FROM molecule_observations WHERE id=?", (observation,)).fetchone()[0]).get("target_molecule_id")


def test_names_never_establish_a_target_and_conflicting_roles_are_rejected():
    assert not target_components([{"type": "protein", "chain": "A", "sequence": "ACD", "name": "Target"}])
    assert not target_components([{"type": "protein", "chain": "A", "sequence": "ACD", "role": "binder"}], candidate={"target_chains": ["A"]})
    assert not target_components([], candidate={"target_chains": ["A"], "binder_chains": ["A"]})
    ligand = {"type": "small_molecule", "chain": "L", "smiles": "N[C@H](C)C(=O)O", "role": "target"}
    assert contents_identity(target_components([ligand])) != contents_identity(target_components([{**ligand, "smiles": "N[C@@H](C)C(=O)O"}]))


def test_missing_target_name_or_conflicting_existing_target_is_flagged(database):
    reference, generated, observation = legacy_design(database)
    with database.connection() as connection:
        connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.target_molecule_id',?) WHERE id=?",
            (generated["molecule_ids"][0], observation))
    plan = plan_curation(database, ADMIN)
    assert not plan["members"] and not plan["observations"]
    assert any("conflicts" in row["reason"] for row in plan["unresolved"])


def test_pipeline_name_never_replaces_the_recorded_design_engine(database):
    target = _import(database, run_id="", kind="reference", method="user_uploads",
        records=[_record("EFG", name="Test target")])
    result = _import(database, records=[_record(name="rfd4_000001", metadata={"source_tool": "rfd4_proteina",
        "design_target_molecule_id": target["molecule_ids"][0], "design_creation": creation()})], origin={"tool": "visual_pipeline"}, run_id="20261101-120000-aa11")
    with database.connection() as connection:
        row = connection.execute("SELECT name,metadata FROM molecule_members WHERE group_id=?", (result["group_id"],)).fetchone()
    assert row["name"] == "Design_pc004033_20261101T123456Z"
    assert json.loads(row["metadata"])["studio_naming"]["tool"] == "rfd4_proteina"
