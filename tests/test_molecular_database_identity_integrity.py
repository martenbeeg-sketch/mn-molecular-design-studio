"""Audit regressions at scientific identity and authorized file boundaries."""
import hashlib
import json

import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.evaluations import evaluation_descriptor
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "audit-owner", "role": "user"}


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
