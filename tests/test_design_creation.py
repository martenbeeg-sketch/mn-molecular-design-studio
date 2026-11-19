from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import os

import pytest

from mn_molecular_design_studio import design_creation
from mn_molecular_design_studio.design_creation import creation_from_structure, generating_workstation, preserve_creation_times, valid_creation
from mn_molecular_design_studio.library_design_creation import CreationRecovery, plan_creation_names
from mn_molecular_design_studio.library_curation import apply_curation
from mn_molecular_design_studio.scientific_naming import design_name
from test_database_browser import USER, _import, _record, database
from test_bound_docking import pdb_atom


def evidence(at="2026-11-01T12:34:56+00:00", sequence="ACDE", **fields):
    return {"schema": design_creation.SCHEMA, "at": at, "basis": "native_cif_creation", "precision": "second",
        "workstation": "pc004033", "node_id": "worker", "source_path": "/original/design.cif",
        "source_sha256": "a" * 64, "binder_sequence_sha256": design_creation.sequence_digest(sequence), **fields}


def native_run(tmp_path, *, remote=False):
    run = tmp_path / "runs/design/20261101-120000-job"
    run.mkdir(parents=True)
    metadata = {"created_at": "2026-11-01T11:00:00+00:00", "completed_at": "2026-11-01T13:00:00+00:00",
        "execution_node_id": "worker", "cluster_remote": remote}
    (run / "metadata.json").write_text(json.dumps(metadata))
    (run / "command.json").write_text(json.dumps({"command": ["docker", "run", "-e", "TZ=UTC", "mn-rfd4-proteina:0.0.3-cu130"]}))
    path = run / "artifacts/rfd4_proteina/output/inference_samples/design_generated.cif"
    path.parent.mkdir(parents=True)
    path.write_text("data_design\n_entry.date 2026-11-01\n_entry.time 12:34:56\n")
    candidate = {"candidate_id": "rfd4_000001", "source_tool": "rfd4_proteina", "stage": "generation.backbone_sequence",
        "binder_sequence": "ACDE", "binder_chains": ["B"], "complex_pdb": str(path), "raw_metadata": {"native_structure": str(path)}}
    return run, path, candidate


def test_names_use_creation_and_actual_generator_not_run_or_import_times():
    metadata = {"design_creation": evidence()}
    label = design_name("rfd4_000001", target_name="Tau", tool="rfd4_proteina", run_id="20261101-110000-launch", metadata=metadata)
    assert label == "Design_pc004033_20261101T123456Z"
    assert not design_name("rfd4_000001", target_name="Tau", tool="rfd4_proteina", run_id="20261101-110000-launch")
    assert not design_name("rfd4_000001", target_name="Tau", tool="rfd4_proteina", run_id="ignored",
        metadata={"design_creation": evidence(at="2026-11-01T12:34:56")})


def test_native_stamp_survives_remote_copy_and_matches_exact_binder(tmp_path, monkeypatch):
    run, path, candidate = native_run(tmp_path, remote=True)
    os.utime(path, (2_000_000_000, 2_000_000_000))
    monkeypatch.setattr(design_creation, "_structure_sequences", lambda path: {"A": "EFG", "B": "ACDE"})
    actual = creation_from_structure(path, candidate, run, names={"worker": "pc004033"}, local_node="coordinator", local_name="pc002975")
    assert actual["at"] == "2026-11-01T12:34:56+00:00"
    assert actual["workstation"] == "pc004033"
    assert actual["basis"] == "native_cif_creation"
    assert actual["source_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert creation_from_structure(path, {**candidate, "binder_sequence": "AAAA"}, run,
        names={"worker": "pc004033"}, local_node="coordinator", local_name="pc002975") is None


def test_naive_cif_time_is_not_assumed_utc_after_a_timezone_override(tmp_path, monkeypatch):
    run, path, candidate = native_run(tmp_path, remote=True)
    (run / "command.json").write_text(json.dumps({"command": ["docker", "run", "-e", "TZ=Europe/Berlin"]}))
    monkeypatch.setattr(design_creation, "_structure_sequences", lambda path: {"B": "ACDE"})
    assert creation_from_structure(path, candidate, run, names={"worker": "pc004033"}, local_node="coordinator",
        local_name="pc002975", allow_file_time=True) is None


def test_original_native_file_time_is_not_taken_from_a_coordinator_copy(tmp_path, monkeypatch):
    run, path, candidate = native_run(tmp_path, remote=True)
    path.write_text("data_no_creation_stamp\n")
    monkeypatch.setattr(design_creation, "_structure_sequences", lambda path: {"B": "ACDE"})
    assert creation_from_structure(path, {**candidate, "source_tool": "bindcraft2"}, run,
        names={"worker": "pc004033"}, local_node="coordinator", local_name="pc002975", allow_file_time=True) is None


def test_file_write_timestamp_is_precise_and_its_limitation_is_preserved(tmp_path, monkeypatch):
    run, path, candidate = native_run(tmp_path)
    epoch = int(datetime(2026, 11, 1, 12, 34, 56, tzinfo=timezone.utc).timestamp())
    nanoseconds = epoch * 1_000_000_000 + 123456789
    os.utime(path, ns=(nanoseconds, nanoseconds))
    monkeypatch.setattr(design_creation, "_structure_sequences", lambda path: {"B": "ACDE"})
    actual = creation_from_structure(path, {**candidate, "source_tool": "bindcraft2"}, run,
        names={"worker": "pc004033"}, local_node="worker", local_name="pc004033", allow_file_time=True)
    assert actual["at"] == "2026-11-01T12:34:56.123456+00:00"
    assert actual["source_mtime_ns"] == nanoseconds
    assert actual["basis"] == "native_structure_written" and actual["limitation"]


def test_existing_creation_survives_renormalization_and_post_analysis(tmp_path, monkeypatch):
    run, path, candidate = native_run(tmp_path)
    original = evidence()
    previous = {**candidate, "raw_metadata": {"design_creation": original}}
    destination = run / "artifacts/normalized_candidates/candidates.jsonl"
    destination.parent.mkdir(parents=True)
    destination.write_text(json.dumps(previous) + "\n")
    monkeypatch.setattr(design_creation, "workstation_names", lambda: {"worker": "pc004033"})
    from mn_molecular_design_studio.cluster import config
    monkeypatch.setattr(config, "identity", lambda: "worker")
    monkeypatch.setattr(design_creation, "candidate_creation", lambda *args, **kwargs: pytest.fail("Replaced original time"))
    assert preserve_creation_times(run, [candidate])[0]["raw_metadata"]["design_creation"] == original
    analysis = {**candidate, "source_tool": "biotite_interface", "stage": "analysis", "raw_metadata": {"source_candidate": previous}}
    assert design_creation.candidate_creation is not None
    monkeypatch.undo()
    assert design_creation.candidate_creation(run, analysis, names={}, local_node="worker", local_name="pc004033") == original
    changed = {**analysis, "binder_sequence": "AAAA"}
    assert design_creation.candidate_creation(run, changed, names={}, local_node="worker", local_name="pc004033") is None


def test_deferred_worker_identity_takes_precedence_over_importing_computer():
    node, name = generating_workstation({"execution_node_id": "coordinator", "deferred_source_job": "cluster://worker/design/source"},
        {"worker": "pc004033", "coordinator": "pc002975"}, local_node="coordinator", local_name="pc002975")
    assert (node, name) == ("worker", "pc004033")


def test_future_sql_registration_disambiguates_real_creation_time_collisions(database):
    first = _record("ACDE", name="rfd4_000001", metadata={"design_creation": evidence()})
    second = _record("AAAA", name="rfd4_000002", metadata={"design_creation": evidence(sequence="AAAA", source_sha256="b" * 64)})
    imported = _import(database, records=[first, second])
    with database.connection() as connection:
        rows = connection.execute("SELECT name,metadata FROM molecule_members WHERE group_id=? ORDER BY name", (imported["group_id"],)).fetchall()
    assert len({row["name"] for row in rows}) == 2
    assert all(row["name"].startswith("Design_pc004033_20261101T123456Z") for row in rows)
    assert all(json.loads(row["metadata"])["studio_naming"]["created_at"] == evidence()["at"] for row in rows)


def test_historical_plan_is_backed_up_idempotent_and_does_not_modify_science(database, tmp_path):
    original_file = tmp_path / "original.cif"
    original_file.write_text("native design evidence")
    original = evidence(source_path=str(original_file), source_sha256=hashlib.sha256(original_file.read_bytes()).hexdigest())
    result = _import(database, records=[_record(name="rfd4_000001", metadata={"design_creation": original})])
    with database.connection() as connection:
        connection.execute("UPDATE molecule_members SET name='legacy run-timestamp label' WHERE group_id=?", (result["group_id"],))
        metadata = json.loads(connection.execute("SELECT metadata FROM molecule_members WHERE group_id=?", (result["group_id"],)).fetchone()[0])
        metadata["studio_naming"]["display_name"] = "legacy run-timestamp label"
        connection.execute("UPDATE molecule_members SET metadata=? WHERE group_id=?", (json.dumps(metadata), result["group_id"]))
    recovery = CreationRecovery(names={"worker": "pc004033"}, local_node="worker", local_name="pc004033", runs_root=tmp_path)
    actor = {**USER, "role": "admin"}
    plan = plan_creation_names(database, actor, recovery)
    assert len(plan["members"]) == 1 and not plan["unresolved"]
    apply_curation(database, actor, plan, backup_path=tmp_path / "backup.sqlite3")
    assert not plan_creation_names(database, actor, recovery)["members"]


def test_reimport_never_replaces_creation_with_a_later_run_or_copy(database):
    original = evidence()
    first = _import(database, records=[_record(name="rfd4_000001", metadata={"design_creation": original})])
    changed = evidence(at="2026-11-02T00:00:00+00:00")
    with database.connection() as connection:
        database._put(connection, first["group_id"], _record(name="renormalized", metadata={"design_creation": changed}))
        row = connection.execute("SELECT name,metadata FROM molecule_members WHERE group_id=?", (first["group_id"],)).fetchone()
    assert row["name"] == "Design_pc004033_20261101T123456Z"
    assert json.loads(row["metadata"])["design_creation"] == original


def test_wrong_sequence_creation_metadata_is_never_used_as_a_name(database):
    result = _import(database, records=[_record("AAAA", name="rfd4_000001", metadata={"design_creation": evidence()})])
    with database.connection() as connection:
        name = connection.execute("SELECT name FROM molecule_members WHERE group_id=?", (result["group_id"],)).fetchone()[0]
    assert name == "rfd4_000001"


def test_creation_metadata_does_not_change_scientific_evaluation_identity():
    from mn_molecular_design_studio.evaluations import evaluation_descriptor
    observation = {"id": "one-evaluation", "molecule_id": "one-design", "namespace": "rfd4_proteina", "run_id": "run",
        "context": {"candidate_id": "native-design", "candidate": {"stage": "generation.backbone_sequence", "source_tool": "rfd4_proteina"},
            "target_molecule_id": "one-target", "structure_digest": "a" * 64, "settings": {"seed": 42}}}
    before = evaluation_descriptor(observation)
    observation["context"]["candidate"]["raw_metadata"] = {"design_creation": evidence()}
    assert evaluation_descriptor(observation) == before


def test_native_bindcraft_does_not_replace_missing_creation_with_rank_promotion_time(tmp_path):
    run, path, candidate = native_run(tmp_path)
    promoted = {**candidate, "source_tool": "bindcraft2", "stage": "complex_refolding",
        "raw_metadata": {"design_creation_structure": None}}
    assert design_creation.candidate_creation(run, promoted, names={"worker": "pc004033"}, local_node="worker", local_name="pc004033") is None
