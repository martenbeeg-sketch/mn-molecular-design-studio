import copy
import json

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from mn_molecular_design_studio import molecule_imports, molecular_database_api
from mn_molecular_design_studio.comparisons import compare, compatibility_warnings, prepare_rows
from mn_molecular_design_studio.evaluations import evaluation_descriptor
from mn_molecular_design_studio.molecular_database import MolecularDatabase

USER = {"id": "owner", "role": "user"}


def database_with_two_targets(tmp_path):
    database = MolecularDatabase(tmp_path)
    imported = database.import_records(USER, molecule_imports.parse_smiles("CCO ethanol"))
    group, molecule = imported["group_id"], imported["molecule_ids"][0]
    first = database.add_observation(group, molecule, "vina", "run-a", {"energy": -7., "rank": 1}, context={"target_id": "target-a"}, created_at="2026-01-01")
    second = database.add_observation(group, molecule, "vina", "run-b", {"energy": -9.}, context={"target_id": "target-b"}, created_at="2026-02-01")
    return database, group, molecule, first, second


def test_target_selection_cannot_mix_metrics_or_fill_missing_rank(tmp_path):
    database, group, molecule, first, second = database_with_two_targets(tmp_path)
    latest = database.list_molecules(USER)["items"][0]
    assert latest["scores"] == {"vina.energy": -9.}
    assert latest["selected_evaluation"]["target_id"] == "target-b"
    old = database.list_molecules(USER, evaluation_filters={"target_id": "target-a"})["items"][0]
    assert old["scores"] == {"vina.energy": -7., "vina.rank": 1}
    assert database.list_molecules(USER, metric_filters=[{"field": "vina.rank", "operator": ">=", "value": 1}])["total"] == 0
    assert len(database.evaluation_rows(USER)) == 2


def test_frozen_comparison_does_not_follow_a_later_run(tmp_path):
    database, group, molecule, first, second = database_with_two_targets(tmp_path)
    selected = database.evaluation_rows(USER, evaluation_filters={"target_id": "target-a"})[0]
    saved = database.save_comparison(USER, "shortlist", [molecule], evaluation_filters={"evaluation_ids": [selected["evaluation"]["id"]]},
        settings={"fields": ["vina.energy"], "directions": {"vina.energy": "lower"}, "method": "mean"})
    database.add_observation(group, molecule, "vina", "run-c", {"energy": -99}, context={"target_id": "target-a"}, created_at="2026-03-01")
    restored = database.list_molecules(USER, group_id=saved["group_id"])["items"][0]
    assert restored["scores"]["vina.energy"] == -7
    assert len(database.evaluation_rows(USER, group_id=saved["group_id"])) == 1


def test_typed_index_invalidates_on_update_and_delete(tmp_path):
    database, group, molecule, first, second = database_with_two_targets(tmp_path)
    database.list_molecules(USER)
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_evaluations").fetchone()[0] == 2
        connection.execute("UPDATE molecule_observations SET metrics=? WHERE id=?", (json.dumps({"energy": -10}), second))
    assert database.list_molecules(USER)["items"][0]["scores"]["vina.energy"] == -10
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_observations WHERE id=?", (second,))
        assert connection.execute("SELECT count(*) FROM molecule_metric_values WHERE observation_id=?", (second,)).fetchone()[0] == 0


def test_selected_structure_is_the_scored_pose(tmp_path):
    database, group, molecule, first, second = database_with_two_targets(tmp_path / "db")
    pose = tmp_path / "pose.pdb"
    pose.write_text("HEADER posed complex\nEND\n")
    structure = database.attach_structure(group, molecule, {"path": str(pose), "label": "Best pose"}, context={"run_id": "run-pose", "candidate_id": "pose-1"})
    database.add_observation(group, molecule, "vina", "run-pose", {"energy": -11}, context={"candidate_id": "pose-1", "target_id": "target-a", "is_representative": True}, created_at="2026-03-01")
    record = database.list_molecules(USER)["items"][0]
    assert record["selected_evaluation"]["representation_id"] == structure
    assert [s["id"] for s in record["structures"]] == [structure]
    old = database.list_molecules(USER, evaluation_filters={"run_id": "run-a"})["items"][0]
    assert old["structures"] == []  # never show an unrelated later pose


def comparison_rows():
    return [{"molecule_id": str(i), "name": str(i), "scores": scores, "evaluation": {"id": "eval-" + str(i), "status": status, "target_id": "t"}}
            for i, scores, status in [(1, {"x": 1, "y": 1, "native.rank": 9}, "completed"),
                                      (2, {"x": 2}, "completed"), (3, {"x": 3, "y": 3}, "failed")]]


@pytest.mark.parametrize("method", ["mean", "worst", "product", "pareto"])
def test_missing_and_failed_rows_never_receive_a_rank(method):
    rows = comparison_rows()
    original = copy.deepcopy(rows)
    frame = compare(rows, fields=["x", "y"], directions={"x": "higher", "y": "higher"}, method=method).set_index("Molecule ID")
    assert frame.loc["1", "Comparison rank"] == 1
    assert frame.loc["1", "Native rank"] == 9
    assert frame.loc["2", "Missing metrics"] == "y"
    assert frame.loc[["2", "3"], "Comparison rank"].isna().all()
    assert rows == original


def test_thresholds_explain_missing_and_actual_failures():
    frame = compare(comparison_rows(), thresholds=[{"field": "y", "operator": ">=", "value": 2}]).set_index("Molecule ID")
    assert frame.loc["1", "Filter failures"] == "y >= 2"
    assert frame.loc["2", "Filter failures"] == "y missing"


def test_conflicting_targets_and_methods_block_binding_ranking_but_not_sequence_properties():
    rows = comparison_rows()[:2]
    rows[1]["evaluation"]["target_id"] = "another"
    assert compatibility_warnings(rows, ["x"])
    assert compatibility_warnings(rows, ["common.sequence.gravy"]) == []


def test_preview_and_save_use_the_same_authenticated_api(tmp_path, monkeypatch):
    database, group, molecule, first, second = database_with_two_targets(tmp_path)
    monkeypatch.setattr(molecular_database_api, "MolecularDatabase", lambda: database)
    app = FastAPI()
    molecular_database_api.install(app, lambda: USER)
    client = TestClient(app)
    payload = {"molecule_ids": [molecule], "evaluation_filters": {"target_id": "target-a"},
               "settings": {"fields": ["vina.energy"], "directions": {"vina.energy": "lower"}}}
    result = client.post("/api/v1/molecules/comparisons/preview", json={"selection": {k: v for k, v in payload.items() if k != "settings"}, "settings": payload["settings"]})
    assert result.status_code == 200, result.text
    assert result.json()["items"][0]["vina.energy"] == -7
    saved = client.post("/api/v1/molecules/comparisons", json={"name": "API shortlist", **payload})
    assert saved.status_code == 201, saved.text
    assert saved.json()["comparison_spec"]["selection_policy"] == "frozen evaluations"


def test_protocol_digest_ignores_execution_placement_and_paths():
    a = {"id": "obs", "context": {"settings": {"ph": 7.4, "cpu_cores": 4, "source_run_dir": "/one"}}}
    b = {"id": "obs", "context": {"settings": {"ph": 7.4, "cpu_cores": 8, "source_run_dir": "/two"}}}
    assert evaluation_descriptor(a)["parameters_digest"] == evaluation_descriptor(b)["parameters_digest"]
    b["context"]["settings"]["ph"] = 5.5
    assert evaluation_descriptor(a)["parameters_digest"] != evaluation_descriptor(b)["parameters_digest"]


def test_evaluation_identity_keeps_exact_structure_pose_and_native_parameters():
    first = {"id": "obs-1", "molecule_id": "mol-a", "group_id": "group-a", "run_id": "run-1",
             "namespace": "ipsae", "context": {
                 "candidate_id": "candidate-1", "source_structure_digest": "a" * 64,
                 "representation_id": "structure-a", "pose_identity": {"id": "pose-1"},
                 "settings": {"pae_cutoff": 10.0, "dist_cutoff": 8.0},
                 "candidate": {"source_tool": "ipsae", "stage": "analysis", "binder_chains": ["B"], "target_chains": ["A"]},
             }}
    second = {**first, "id": "obs-2", "context": {**first["context"],
              "source_structure_digest": "b" * 64, "representation_id": "structure-b",
              "pose_identity": {"id": "pose-2"}}}
    changed_settings = {**first, "id": "obs-3", "context": {**first["context"],
                       "settings": {"pae_cutoff": 8.0, "dist_cutoff": 8.0}}}

    first_descriptor = evaluation_descriptor(first)
    second_descriptor = evaluation_descriptor(second)
    changed_descriptor = evaluation_descriptor(changed_settings)
    assert first_descriptor["structure_digest"] == "a" * 64
    assert first_descriptor["pose_id"] == "pose-1"
    assert first_descriptor["native_parameters"] == {"pae_cutoff": 10.0, "dist_cutoff": 8.0}
    assert first_descriptor["normalized_parameters"]["analysis.ipsae.pae_cutoff"] == 10.0
    assert first_descriptor["subject_id"] != second_descriptor["subject_id"]
    assert first_descriptor["instance_id"] != second_descriptor["instance_id"]
    assert first_descriptor["parameters_digest"] != changed_descriptor["parameters_digest"]


def test_nested_model_and_ligand_residue_selections_are_part_of_evaluation_identity():
    base = {"id": "obs-1", "molecule_id": "mol-a", "group_id": "group-a", "run_id": "run-1",
            "namespace": "plip", "context": {"source_structure_digest": "a" * 64,
                "candidate": {"candidate_id": "pose-1", "raw_metadata": {
                    "complex_analysis_selection": {"model_id": "2", "ligand_selection": {
                        "chain": "L", "residue_key": "L:501", "ccd": "LIG"}}}}}}
    first = evaluation_descriptor(base)
    changed_model = copy.deepcopy(base)
    changed_model["id"] = "obs-2"
    changed_model["context"]["candidate"]["raw_metadata"]["complex_analysis_selection"]["model_id"] = "3"
    changed_residue = copy.deepcopy(base)
    changed_residue["id"] = "obs-3"
    changed_residue["context"]["candidate"]["raw_metadata"]["complex_analysis_selection"]["ligand_selection"]["residue_key"] = "L:502"

    model = evaluation_descriptor(changed_model)
    residue = evaluation_descriptor(changed_residue)
    assert first["model_id"] == "2"
    assert first["ligand_selection"]["residue_key"] == "L:501"
    assert first["subject_id"] != model["subject_id"]
    assert first["instance_id"] != model["instance_id"]
    assert first["selection_id"] != residue["selection_id"]
    assert first["subject_id"] != residue["subject_id"]
    assert first["instance_id"] != residue["instance_id"]


def test_legacy_evaluation_projection_upgrade_is_bounded_and_preserves_raw_observations(tmp_path):
    from mn_molecular_design_studio.evaluations import PROJECTION_VERSION, rebuild_missing

    database = MolecularDatabase(tmp_path / "projection-upgrade")
    imported = database.import_records(USER, molecule_imports.parse_smiles("CCO"))
    group, molecule = imported["group_id"], imported["molecule_ids"][0]
    first = database.add_observation(group, molecule, "rosetta_metrics", "legacy-rosetta", {
        "interface_percent_of_binder_in_intf": 17.0, "score": -4.0,
    }, context={"structure_digest": "pose-1"})
    second = database.add_observation(group, molecule, "vina", "current", {"energy": -7.0})
    database.evaluation_rows(USER)  # Build current projections before isolating legacy rows.

    with database.connection() as connection:
        raw_before = connection.execute(
            "SELECT metrics,context FROM molecule_observations WHERE id=?", (first,)
        ).fetchone()
        connection.execute("UPDATE molecule_evaluations SET projection_version='evaluation.v5:studio.output.v1'")
        ordered = [row[0] for row in connection.execute(
            "SELECT observation_id FROM molecule_evaluations ORDER BY observation_id")]
        assert set(ordered) == {first, second}
        connection.execute("DELETE FROM molecule_metric_values WHERE observation_id=?", (first,))

        assert rebuild_missing(connection, batch_size=1) == 1
        versions = dict(connection.execute(
            "SELECT observation_id,projection_version FROM molecule_evaluations").fetchall())
        assert versions[ordered[0]] == PROJECTION_VERSION
        assert versions[ordered[1]] == "evaluation.v5:studio.output.v1"
        assert connection.execute(
            "SELECT count(*) FROM molecule_metric_values WHERE observation_id=? AND field='vina.energy'", (second,)
        ).fetchone()[0] == 1
        assert connection.execute(
            "SELECT count(*) FROM molecule_metric_values WHERE observation_id=? "
            "AND field='common.geometry.rosetta_metrics.interface_percent_of_binder_in_intf'", (first,)
        ).fetchone()[0] == 0
        assert tuple(connection.execute(
            "SELECT metrics,context FROM molecule_observations WHERE id=?", (first,)
        ).fetchone()) == tuple(raw_before)

        assert rebuild_missing(connection, batch_size=1) == 1
        assert connection.execute(
            "SELECT count(*) FROM molecule_metric_values WHERE observation_id=? AND field='vina.energy'", (second,)
        ).fetchone()[0] == 1
        with pytest.raises(ValueError, match="between 1 and 5,000"):
            rebuild_missing(connection, batch_size=5001)


def test_comparison_view_restores_exact_saved_rules():
    from streamlit.testing.v1 import AppTest
    settings = {"fields": ["x", "y"], "directions": {"x": "lower", "y": "higher"}, "weights": {"x": 2., "y": 1.},
                "thresholds": [{"field": "x", "operator": "<", "value": 3.}], "method": "worst"}
    app = AppTest.from_string(f'''from mn_molecular_design_studio.comparison_ui import render_comparison
render_comparison({comparison_rows()!r}, key="test", settings={settings!r})''').run()
    assert not app.exception
    assert app.multiselect(key="test-fields").value == ["x", "y"]
    assert app.selectbox(key="test-method").value == "worst"
    assert any(widget.value == "<" for widget in app.selectbox)


def test_analysis_does_not_substitute_complex_confidence_or_guess_scale():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import build_analysis_table
    candidates = [{"candidate_id": "a", "source_tool": "unknown", "metrics": {"plddt": .9, "boltz2_complex_plddt": .8, "boltz2_confidence_score": .7}}]
    row = build_analysis_table(candidates, thresholds={"min_binder_plddt": 70.}).iloc[0]
    assert row["binder_plddt"] is None
    assert "binder_plddt missing" in row["filter_failures"]


def test_analysis_uses_identical_metrics_for_every_derived_score():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import build_analysis_table
    candidates = [{"candidate_id": "a", "source_tool": "x", "metrics": {"binder_plddt": 90., "iptm": .9}},
                  {"candidate_id": "b", "source_tool": "x", "metrics": {"binder_plddt": 80.}}]
    frame = build_analysis_table(candidates)
    assert set(frame["comparison_metrics"]) == {"binder_plddt"}


def same_structure_rows():
    rows = []
    for molecule, confidence, solubility in [("a", .9, .6), ("b", .8, .5)]:
        for tool, scores in [("prediction", {"confidence": confidence}), ("proteinsol", {"solubility": solubility})]:
            rows.append({"molecule_id": molecule, "name": molecule, "scores": scores,
                         "evaluation": {"id": molecule + tool, "subject_id": "subject-" + molecule,
                                        "representation": "digest-" + molecule, "status": "completed", "target_id": "t",
                                        "evaluator": tool, "run_id": tool + "-run", "parameters_digest": tool + "-params"}})
    return rows


def test_explicit_composition_joins_distinct_evaluators_of_the_identical_structure():
    rows = same_structure_rows()
    original = copy.deepcopy(rows)
    frame = compare(rows, fields=["confidence", "solubility"], directions={"confidence": "higher", "solubility": "higher"}, combine_same_subject=True)
    assert len(frame) == 2
    assert frame["Comparison rank"].tolist() == [1, 2]
    assert frame.iloc[0]["Source evaluations"] == ["aprediction", "aproteinsol"]
    assert rows == original
    assert compatibility_warnings(prepare_rows(rows, combine_same_subject=True), ["confidence", "solubility"]) == []


def test_composition_never_joins_different_structures_or_incomplete_context():
    rows = same_structure_rows()[:2]
    rows[1]["evaluation"]["subject_id"] = "other-pose"
    assert len(prepare_rows(rows, combine_same_subject=True)) == 2
    rows[1]["evaluation"]["subject_id"] = rows[0]["evaluation"]["subject_id"] = ""
    assert len(prepare_rows(rows, combine_same_subject=True)) == 2


def test_conflicting_repetitions_are_missing_until_a_contributing_run_is_selected():
    rows = same_structure_rows()[:2]
    rows.append(copy.deepcopy(rows[1]))
    rows[2]["evaluation"].update(id="repeat", run_id="repeat", parameters_digest="different-protocol")
    # Even equal values must not conceal differing protocols.
    combined = prepare_rows(rows, combine_same_subject=True)[0]
    assert "solubility" not in combined["scores"]
    assert combined["conflicting_metrics"] == ["solubility"]
    restored = prepare_rows(rows, combine_same_subject=True, evaluator_runs=["prediction / prediction-run", "proteinsol / proteinsol-run"])
    assert restored[0]["scores"]["solubility"] == .6
    assert prepare_rows(rows, evaluator_runs=[]) == []


def test_saved_composed_comparison_freezes_all_contributing_evaluations(tmp_path):
    database = MolecularDatabase(tmp_path)
    imported = database.import_records(USER, molecule_imports.parse_sequences("ACDE"))
    molecule, group = imported["molecule_ids"][0], imported["group_id"]
    context = {"target_id": "t", "structure_digest": "identical-structure"}
    database.add_observation(group, molecule, "prediction", "p", {"confidence": .9}, context=context)
    database.add_observation(group, molecule, "proteinsol", "q", {"solubility": .5}, context=context)
    database.add_observation(group, molecule, "proteinsol", "excluded", {"solubility": .6}, context=context)
    settings = {"combine_same_subject": True, "evaluator_runs": ["prediction / p", "proteinsol / q"],
                "fields": ["prediction.confidence", "proteinsol.solubility"],
                "directions": {"prediction.confidence": "higher", "proteinsol.solubility": "higher"}}
    saved = database.save_comparison(USER, "Composed", [molecule], settings=settings)
    assert len(saved["comparison_spec"]["evaluation_filters"]["evaluation_ids"]) == 2
    database.add_observation(group, molecule, "proteinsol", "later", {"solubility": 999}, context=context)
    restored = compare(database.evaluation_rows(USER, group_id=saved["group_id"]), **settings)
    assert len(restored) == 1
    assert restored.iloc[0]["proteinsol.solubility"] == .5
    assert restored.iloc[0]["Comparison rank"] == 1


@pytest.mark.parametrize("settings", [{"fields": "bad"}, {"combine_same_subject": "yes"}, {"evaluator_runs": [42]}, {"thresholds": ["bad"]}])
def test_preview_rejects_malformed_settings_as_validation_errors(settings, tmp_path, monkeypatch):
    database, group, molecule, _, _ = database_with_two_targets(tmp_path)
    monkeypatch.setattr(molecular_database_api, "MolecularDatabase", lambda: database)
    app = FastAPI()
    molecular_database_api.install(app, lambda: USER)
    result = TestClient(app).post("/api/v1/molecules/comparisons/preview", json={"selection": {"molecule_ids": [molecule]}, "settings": settings})
    assert result.status_code == 422, result.text


def test_comparison_view_restores_composed_runs_and_combines_rows():
    from streamlit.testing.v1 import AppTest
    settings = {"combine_same_subject": True, "evaluator_runs": ["prediction / prediction-run", "proteinsol / proteinsol-run"],
                "fields": ["confidence", "solubility"], "directions": {"confidence": "higher", "solubility": "higher"}}
    app = AppTest.from_string(f'''from mn_molecular_design_studio.comparison_ui import render_comparison
render_comparison({same_structure_rows()!r}, key="composed", settings={settings!r})''').run()
    assert not app.exception
    assert app.checkbox(key="composed-combine").value is True
    assert app.multiselect(key="composed-runs").value == settings["evaluator_runs"]
    assert len(app.dataframe[0].value) == 2
