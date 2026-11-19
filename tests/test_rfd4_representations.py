import json
from pathlib import Path

import gemmi
import pytest

from mn_molecular_design_studio.rfd4_representations import recorded_oracle_metrics, refolded_representation, recover_saved_refold
from mn_molecular_design_studio.engine.workflow_engine.workflows import rfd4_proteina as workflow
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from mn_molecular_design_studio.data_browser import listed_results, load_saved_record, saved_results
from mn_molecular_design_studio.binderscope import prepare_browser_selection
from test_database_browser import USER, _import, _record, database
from test_molecular_database_policy import _pdb_efg_acd_complex


@pytest.fixture
def models(tmp_path):
    directory = tmp_path / "artifacts/rfd4_proteina/output/inference_samples/example"
    directory.mkdir(parents=True)
    generated = directory / "example_generated.cif"
    predicted = directory / "example_of3_ppi_refolded.cif"
    structure = gemmi.read_pdb_string(_pdb_efg_acd_complex())
    structure.make_mmcif_document().write_file(str(generated))
    for atom in structure[0].all():
        atom.atom.pos.x += 1
        atom.atom.b_iso = 85
    structure.make_mmcif_document().write_file(str(predicted))
    confidence = directory / "example_of3_ppi_confidences.json"
    confidence.write_text(json.dumps({"overall_plddt": 85, "plddt_coi": 85}))
    return generated, predicted, confidence


def test_exact_oracle_receipt_and_composition_are_required(models):
    generated, predicted, confidence = models
    metrics = {"oracle_ppi.overall_plddt": 85}
    result = refolded_representation(generated, metrics)
    assert result["refolded_path"] == str(predicted)
    assert result["generated_sha256"] != result["refolded_sha256"]
    with pytest.raises(ValueError, match="disagrees"):
        refolded_representation(generated, {"oracle_ppi.overall_plddt": 30})
    structure = gemmi.read_structure(str(predicted))
    structure[0]["B"][0].name = "GLY"
    structure.make_mmcif_document().write_file(str(predicted))
    with pytest.raises(ValueError, match="identities"):
        refolded_representation(generated, metrics)


def test_curated_receipt_recovers_only_candidate_specific_native_measurements():
    metadata = {"source_candidate_id": "accepted-design", "validation_acceptance": {"measurements": {
        "common.plddt": {"status": "measured", "measurement_candidate_id": "accepted-design", "value": 0.85,
            "output_provenance": {"plddt": {"aliases": [{"source_field": "oracle_ppi.overall_plddt", "raw_value": 85}]}}},
        "common.iptm": {"status": "measured", "measurement_candidate_id": "another-design", "value": 0.9,
            "output_provenance": {"iptm": {"aliases": [{"source_field": "oracle_ppi.overall_iptm", "raw_value": 0.9}]}}}}}}
    assert recorded_oracle_metrics(metadata) == {"oracle_ppi.overall_plddt": 85}
    metadata["source_candidate_id"] = "missing-design"
    assert recorded_oracle_metrics(metadata) == {}


def test_native_adapter_prefers_refold_and_retains_generated_artifact(models, tmp_path):
    generated, predicted, confidence = models
    root = tmp_path / "artifacts/rfd4_proteina"
    (root / "metrics.csv").write_text("designed_structure_path,oracle_ppi.overall_plddt\n" + str(generated) + ",85\n")
    (root / "result_index.json").write_text(json.dumps({"structures": [{"structure": str(generated),
        "designed_chains": ["B"], "target_chains": ["A"], "sequences": {"A": "EFG", "B": "ACD"}, "residue_count": 6}],
        "metrics_files": [str(root / "metrics.csv")]}))
    prepared = {"root": root, "parameters": {"length_min": 1, "length_max": 10, "options": {"refold": True}}}
    candidate = workflow.normalize_outputs(tmp_path, prepared, "protein_binder", 1)[0]
    assert candidate["complex_pdb"] == str(predicted)
    assert candidate["raw_metadata"]["generated_complex_pdb"] == str(generated)
    from mn_molecular_design_studio.molecular_database_ingest import _record
    record = _record(candidate, tmp_path)
    assert {row["label"] for row in record["files"]} == {"Refolded validation complex", "Generated design complex"}
    assert record["components"][1]["sequence"] == "ACD"


def legacy_design(database, models):
    generated, predicted, confidence = models
    target = _import(database, kind="reference", task="reference", run_id="target", records=[_record("EFG")])["molecule_ids"][0]
    metadata = {"source_tool": "rfd4_proteina", "task": "protein_binder", "native_structure": str(generated),
        "reference_constraints": {"options": {"refold": True}}, "design_target_molecule_id": target,
        "validation_acceptance": {"status": "accepted", "binder_format": "miniprotein"}}
    imported = _import(database, records=[{"name": "Saved design", "components": [
        {"type": "protein", "chain": "A", "sequence": "EFG", "role": "target"},
        {"type": "protein", "chain": "B", "sequence": "ACD", "role": "binder"}], "metadata": metadata,
        "files": [{"path": str(generated), "label": "Generated complex", "binder_chains": ["B"], "target_chains": ["A"]}]}])
    structure = database.list_molecules(USER, molecule_ids=imported["molecule_ids"], collection="all")["items"][0]["structures"][0]
    database.add_observation(imported["group_id"], imported["molecule_ids"][0], "rfd4_proteina", "run-1", {"oracle_ppi.overall_plddt": 85},
        context={"candidate_id": "native-example", "representation_id": structure["id"], "target_id": target,
            "candidate": {"stage": "generation.backbone_sequence", "binder_chains": ["B"], "target_chains": ["A"],
                "raw_metadata": {"validation_acceptance": metadata["validation_acceptance"]}}})
    return imported


def inventory(database):
    index = build_browser_index(database, USER, collection="library", include_related=True)
    rows = subject_results(database, index, {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])})
    return index, saved_results(index, rows)


def test_recovery_is_database_owned_idempotent_and_keeps_one_design_with_both_models(database, models):
    imported = legacy_design(database, models)
    index, before = inventory(database)
    original = before[0]["primary"]
    database.make_structure_portable(USER, original["structure_ids"][0])
    plan = recover_saved_refold(database, USER, original)
    assert plan["status"] == "verified"
    stored = recover_saved_refold(database, USER, original, apply=True)
    assert stored["status"] == "stored"
    repeated = recover_saved_refold(database, USER, original, apply=True)
    assert repeated["structure_id"] == stored["structure_id"]
    assert repeated["observation_id"] == stored["observation_id"]
    with database.connection() as connection:
        path = connection.execute("SELECT path FROM molecule_structures WHERE id=?", (stored["structure_id"],)).fetchone()[0]
    assert Path(path).is_relative_to(database.root)
    index, after = inventory(database)
    assert len(after) == len(before) == 1 and after[0]["id"] == before[0]["id"]
    assert after[0]["primary"]["receipt"]["representation_role"] == "refolded_validation"
    from mn_molecular_design_studio.database_design_evaluations import connected_evaluations
    connected, unresolved = connected_evaluations(index, after[0]["primary"], after[0]["originals"])
    assert not unresolved
    assert next(row for row in connected if row["evaluation_id"] == original["evaluation_id"])["structure_relation"] == "Source model (not selected coordinates)"
    listing = listed_results(after)
    selection, excluded = prepare_browser_selection(database, USER, listing)
    assert not excluded
    assert selection["sources"][0]["structure_id"] == stored["structure_id"]
    assert selection["sources"][0]["assessment"]["b_factor_interpretation"] == "plddt_0_100"
    record = load_saved_record(database, USER, listing[0]["result"])
    assert len(record["structures"]) == 2
    assert {row["label"] for row in record["structures"]} == {"Refolded validation complex", "Generated design complex"}
    models[1].unlink()
    assert database.structure_file(USER, stored["structure_id"]).is_file()
    all_models = listed_results(after, structures=True)
    assert len(all_models) == 2
    generated = next(row for row in all_models if row["structure_id"] == original["structure_ids"][0])
    assert generated["result"]["evaluation_id"] != after[0]["primary"]["evaluation_id"]


def test_missing_refold_is_not_silently_replaced_for_default_binderscope(database, models):
    legacy_design(database, models)
    index, outcomes = inventory(database)
    selection, excluded = prepare_browser_selection(database, USER, listed_results(outcomes))
    assert not selection["sources"]
    assert "refolded validation complex must be recovered" in excluded[0]["Reason"]
    explicit = listed_results(outcomes, structures=True)
    selection, excluded = prepare_browser_selection(database, USER, explicit)
    assert not excluded and len(selection["sources"]) == 1


def test_peer_replay_preserves_both_models_and_prefers_verified_refold(database, models, tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import transport
    from mn_molecular_design_studio.cluster.replication import Replica
    from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id
    from mn_molecular_design_studio.molecular_database_peer import materialize_portable_structures
    from test_molecular_database_peer import SOURCE, _sync_all

    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    legacy_design(database, models)
    index, outcomes = inventory(database)
    original = outcomes[0]["primary"]
    database.make_structure_portable(USER, original["structure_ids"][0])
    stored = recover_saved_refold(database, USER, original, apply=True)
    receiver = MolecularDatabase(tmp_path / "receiver")
    replica = Replica()
    _sync_all(replica, database, receiver)
    monkeypatch.setattr(transport, "node_request", lambda endpoint, *args, **kwargs:
        database.structure_file(USER, endpoint.rsplit("/", 1)[-1]).read_bytes())
    assert materialize_portable_structures(database=receiver, peer={"id": SOURCE, "url": "https://source.invalid"}) == 2
    _sync_all(replica, database, receiver)
    index, restored = inventory(receiver)
    assert len(restored) == 1
    primary = restored[0]["primary"]
    assert primary["receipt"]["representation_role"] == "refolded_validation"
    assert primary["source_structure_id"] == _id("peer-structure-", SOURCE + ":" + original["structure_ids"][0])
    assert primary["structure_ids"] == [_id("peer-structure-", SOURCE + ":" + stored["structure_id"])]
    assert len(load_saved_record(receiver, USER, primary)["structures"]) == 2
    selection, excluded = prepare_browser_selection(receiver, USER, listed_results(restored))
    assert not excluded and selection["sources"][0]["assessment"]["b_factor_interpretation"] == "plddt_0_100"


def test_new_run_ingestion_keeps_both_portable_models_and_binds_scores_to_refold(database, models):
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from test_molecular_database_policy import _make_run

    generated, predicted, confidence = models
    representation = refolded_representation(generated, {"oracle_ppi.overall_plddt": 85})
    candidate = {"candidate_id": "new-native-design", "source_tool": "rfd4_proteina", "stage": "generation.backbone_sequence",
        "complex_pdb": str(predicted), "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
        "metrics": {"oracle_ppi.overall_plddt": 85}, "raw_metadata": {"generated_complex_pdb": str(generated),
            "design_representation": representation, "task": "protein_binder", "reference_constraints": {"options": {"refold": True}},
            "validation_acceptance": {"status": "accepted", "binder_format": "miniprotein"}}}
    run = _make_run(database, "design", "native-refold", project="", tool="rfd4_proteina", workflow="design.rfd4_proteina", candidates=[candidate])
    assert ingest_run(run, database=database) == 1
    index, outcomes = inventory(database)
    assert len(outcomes) == 1
    primary = outcomes[0]["primary"]
    assert primary["receipt"]["representation_role"] == "refolded_validation"
    assert primary["source_structure_id"] not in primary["structure_ids"]
    saved = load_saved_record(database, USER, primary)
    assert len(saved["structures"]) == 2
    assert {row["label"] for row in saved["structures"]} == {"Generated design complex", "Refolded validation complex"}
    assert saved["selected_evaluation"]["structure_digest"] == representation["refolded_sha256"]
    selection, excluded = prepare_browser_selection(database, USER, listed_results(outcomes))
    assert not excluded and selection["sources"][0]["assessment"]["b_factor_interpretation"] == "plddt_0_100"
