"""Audit regressions at scientific identity and authorized file boundaries."""
import hashlib
import json

import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.evaluations import evaluation_descriptor
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user"}


@pytest.fixture
def subject(tmp_path):
    database = MolecularDatabase(tmp_path / "db")
    imported = database.import_records(USER, molecule_imports.parse_sequences(">binder\nACDEFG"))
    return database, imported["group_id"], imported["molecule_ids"][0]


def attach(subject, tmp_path, name, **identity):
    database, group, molecule = subject
    path = tmp_path / f"{name}.cif"
    path.write_text(f"data_{name}\n#\n")
    structure = database.attach_structure(group, molecule, {"path": str(path)}, context={
        "run_id": "run-1", "candidate_id": "candidate-1", **identity,
    })
    return structure, hashlib.sha256(path.read_bytes()).hexdigest()


def observe(subject, context, *, source_key="same-measurement"):
    database, group, molecule = subject
    oid = database.add_observation(group, molecule, "prediction", "run-1", {"confidence": .9},
                                   context=context, source_key=source_key)
    with database.connection() as connection:
        row = dict(connection.execute("SELECT * FROM molecule_observations WHERE id=?", (oid,)).fetchone())
    row["context"] = json.loads(row["context"])
    return row


@pytest.mark.parametrize("reverse", [False, True])
def test_auto_attachment_requires_the_requested_model_regardless_of_order(subject, tmp_path, reverse):
    ids = {}
    for index in ([2, 1] if reverse else [1, 2]):
        ids[index], _ = attach(subject, tmp_path, f"model{index}", model_id=f"model-{index}")
    context = {"candidate_id": "candidate-1", "model_id": "model-1"}
    observation = observe(subject, context)
    assert observation["context"]["representation_id"] == ids[1]
    descriptor = evaluation_descriptor(observation)
    assert (descriptor["model_id"], descriptor["representation_id"]) == ("model-1", ids[1])
    assert observe(subject, context)["id"] == observation["id"]


@pytest.mark.parametrize("selector", [
    {}, {"candidate_id": "candidate-1"}, {"model_id": "absent"},
    {"candidate_id": "wrong", "model_id": "model-1"},
    {"candidate_id": "candidate-1", "structure_digest": "f" * 64},
])
def test_ambiguous_or_conflicting_identity_remains_unattached(subject, tmp_path, selector):
    attach(subject, tmp_path, "one", model_id="model-1")
    attach(subject, tmp_path, "two", model_id="model-2")
    assert not observe(subject, selector)["context"].get("representation_id")


@pytest.mark.parametrize("selector", [
    {"model_index": 0}, {"candidate": {"model_index": 0}},
    {"candidate": {"raw_metadata": {"complex_analysis_selection": {"model_index": 0}}}},
])
def test_zero_model_index_is_preserved_in_new_observations(subject, tmp_path, selector):
    first, _ = attach(subject, tmp_path, "zero", model_index=0)
    attach(subject, tmp_path, "one", model_index=1)
    observation = observe(subject, selector)
    assert observation["context"]["representation_id"] == first
    assert evaluation_descriptor(observation)["model_id"] == "0"


def test_pose_and_digest_constraints_are_not_replaced_by_recency(subject, tmp_path):
    first, digest = attach(subject, tmp_path, "pose1", pose_id="pose-1", pose_rank=0)
    attach(subject, tmp_path, "pose2", pose_id="pose-2", pose_rank=1)
    for selector in ({"pose_identity": {"id": "pose-1"}}, {"pose_rank": 0},
                     {"structure_digest": digest}, {"source_structure_digest": digest}):
        observation = observe(subject, selector, source_key=json.dumps(selector))
        assert observation["context"]["representation_id"] == first
        assert observation["context"].get("structure_digest") == digest


def test_source_digest_does_not_overwrite_calculation_digest(subject, tmp_path):
    first, digest = attach(subject, tmp_path, "original", model_id="1")
    observation = observe(subject, {"source_structure_digest": digest, "structure_digest": "a" * 64,
                                    "calculation_structure_digest": "a" * 64})
    assert observation["context"]["representation_id"] == first
    assert observation["context"]["structure_digest"] == "a" * 64


def test_numeric_zero_model_id_remains_visible_in_evaluation_identity(subject, tmp_path):
    first, _ = attach(subject, tmp_path, "zero-id", model_id=0)
    attach(subject, tmp_path, "one-id", model_id=1)
    observation = observe(subject, {"model_id": 0})
    assert observation["context"]["representation_id"] == first
    assert evaluation_descriptor(observation)["model_id"] == "0"


def test_identical_models_with_different_source_candidates_are_ambiguous(subject, tmp_path):
    attach(subject, tmp_path, "one", candidate_id="candidate-1", model_id="1")
    attach(subject, tmp_path, "two", candidate_id="candidate-2", model_id="1")
    assert not observe(subject, {"model_id": "1"})["context"].get("representation_id")


def test_peer_qualified_representation_stays_with_its_source_group(subject, tmp_path):
    database, group, molecule = subject
    first, _ = attach(subject, tmp_path, "owner", model_id="1")
    peer = database.create_group(USER, "Peer collection")
    database.add_to_group(USER, peer, [molecule])
    database.attach_structure(peer, molecule, {"name": "peer.cif", "content": b"data_peer\n"},
                              context={"run_id": "run-1", "candidate_id": "candidate-1", "model_id": "1",
                                       "_replica_node": "peer-b", "_replica_structure": "remote-model"})
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET origin=? WHERE id=?", (json.dumps({"_replica_node": "peer-b"}), peer))
    assert observe(subject, {"candidate_id": "candidate-1", "model_id": "1"})["context"]["representation_id"] == first


@pytest.mark.parametrize("operation", ["download", "materialize", "portable"])
def test_changed_source_is_rejected_before_any_scientific_reuse(subject, tmp_path, operation):
    database, _, _ = subject
    sid, _ = attach(subject, tmp_path, "changed", model_id="1")
    (tmp_path / "changed.cif").write_text("data_different_model\n")
    with pytest.raises(ValueError, match="digest"):
        if operation == "portable":
            database.make_structure_portable(USER, sid)
        else:
            database.structure_file(USER, sid, materialize=operation == "materialize")


def test_verified_snapshot_survives_source_replacement(subject, tmp_path):
    database, _, _ = subject
    sid, digest = attach(subject, tmp_path, "snapshot", model_id="1")
    source = tmp_path / "snapshot.cif"
    original = source.read_bytes()
    delivered = database.structure_file(USER, sid)
    replacement = tmp_path / "replacement.cif"
    replacement.write_text("data_replacement\n")
    replacement.replace(source)
    assert delivered != source
    assert delivered.read_bytes() == original
    assert hashlib.sha256(delivered.read_bytes()).hexdigest() == digest


@pytest.mark.parametrize("digest_field", ["sha256", "portable_sha256"])
def test_corrupted_local_peer_copy_is_rejected(subject, tmp_path, digest_field):
    database, _, _ = subject
    sid, digest = attach(subject, tmp_path, "peer-copy", model_id="1")
    with database.connection() as connection:
        row = connection.execute("SELECT context FROM molecule_structures WHERE id=?", (sid,)).fetchone()
        context = {**json.loads(row["context"]), digest_field: digest, "_replica_node": "peer-b"}
        connection.execute("UPDATE molecule_structures SET context=? WHERE id=?", (json.dumps(context), sid))
    (tmp_path / "peer-copy.cif").write_text("data_corrupt\n")
    with pytest.raises(ValueError, match="digest"):
        database.structure_file(USER, sid)


def test_legacy_digestless_external_source_is_explicitly_unverified(subject, tmp_path):
    database, _, _ = subject
    sid, _ = attach(subject, tmp_path, "legacy", model_id="1")
    with database.connection() as connection:
        connection.execute("UPDATE molecule_structures SET context='{}' WHERE id=?", (sid,))
    with pytest.raises(ValueError, match="[Ll]egacy.*digest"):
        database.structure_file(USER, sid)
    with database.connection() as connection:
        assert connection.execute("SELECT context FROM molecule_structures WHERE id=?", (sid,)).fetchone()[0] == "{}"


def test_legacy_content_addressed_file_can_be_verified_without_rewriting_metadata(subject):
    database, group, molecule = subject
    content = b"data_legacy_stored\n"
    sid = database.attach_structure(group, molecule, {"name": "stored.cif", "content": content})
    with database.connection() as connection:
        connection.execute("UPDATE molecule_structures SET context='{}' WHERE id=?", (sid,))
    assert database.structure_file(USER, sid).read_bytes() == content


def test_new_file_import_records_digest_and_rejects_corrupt_cas_content(subject, tmp_path):
    database, _, _ = subject
    source = tmp_path / "external.cif"
    source.write_text("data_external\n")
    imported = database.import_records(USER, [{"name": "external", "components": [
        {"chain": "A", "type": "protein", "sequence": "ACDEFG"}], "files": [{"path": str(source)}]}])
    with database.connection() as connection:
        structure = connection.execute("SELECT * FROM molecule_structures WHERE group_id=?", (imported["group_id"],)).fetchone()
    assert json.loads(structure["context"])["content_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    content = b"data_cas\n"
    stored = database._store_file("cas.cif", content)
    stored.write_text("data_corrupt\n")
    with pytest.raises(ValueError, match="digest"):
        database._store_file("cas.cif", content)


def test_corrupt_materialized_input_is_replaced_only_from_verified_source(subject, tmp_path):
    database, _, _ = subject
    sid, digest = attach(subject, tmp_path, "input-cache", model_id="1")
    snapshot = database.structure_file(USER, sid, materialize=True)
    snapshot.chmod(0o600)
    snapshot.write_text("data_wrong_cached_input\n")
    fresh = database.structure_file(USER, sid, materialize=True)
    assert hashlib.sha256(fresh.read_bytes()).hexdigest() == digest


def test_missing_and_unauthorized_structures_never_reuse_a_prior_snapshot(subject, tmp_path):
    database, _, _ = subject
    sid, _ = attach(subject, tmp_path, "unavailable", model_id="1")
    assert database.structure_file(USER, sid).is_file()
    (tmp_path / "unavailable.cif").unlink()
    with pytest.raises(PermissionError):
        database.structure_file({"id": "b" * 32, "role": "user"}, sid)
    with pytest.raises(FileNotFoundError):
        database.structure_file(USER, sid)


def test_download_api_rejects_changed_bytes_and_serves_a_verified_snapshot(subject, tmp_path, monkeypatch):
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio import molecular_database_api

    database, _, _ = subject
    sid, _ = attach(subject, tmp_path, "api", model_id="1")
    source = tmp_path / "api.cif"
    original = source.read_bytes()
    monkeypatch.setattr(molecular_database_api, "MolecularDatabase", lambda: database)
    app = FastAPI()
    molecular_database_api.install(app, lambda: USER)
    client = TestClient(app)
    source.write_text("data_changed\n")
    assert client.get(f"/api/v1/molecules/structures/{sid}/file").status_code == 422
    assert client.post(f"/api/v1/molecules/structures/{sid}/materialize").status_code == 422
    source.write_bytes(original)
    structure_file = database.structure_file

    def swap_after_read(*args, **kwargs):
        path = structure_file(*args, **kwargs)
        source.write_text("data_changed_after_validation\n")
        return path

    monkeypatch.setattr(database, "structure_file", swap_after_read)
    response = client.get(f"/api/v1/molecules/structures/{sid}/file")
    assert response.status_code == 200
    assert response.content == original


def test_database_analysis_adapter_rejects_changed_source_without_an_evaluation(subject, tmp_path):
    from mn_molecular_design_studio.complex_analysis_sources import authorize_database_structure

    database, group, molecule = subject
    sid, _ = attach(subject, tmp_path, "analysis", model_id="1")
    (tmp_path / "analysis.cif").write_text("data_different\n")
    with pytest.raises(ValueError, match="digest"):
        authorize_database_structure(USER, {"molecule_id": molecule, "group_id": group,
                                              "structure_id": sid, "model_id": "1"},
                                         project_id="", database=database, workflow_id="analysis.biotite_interface")
