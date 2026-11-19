from __future__ import annotations

from collections import Counter
from copy import deepcopy
import json

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects, subject_results
from mn_molecular_design_studio.database_browser_cache import BrowserCache
from mn_molecular_design_studio.database_outcomes import (
    filter_outcomes, outcome_receipt, scientific_outcomes, workflow_role,
)
from test_database_browser import USER, _app_source, _import, _record, database
from test_database_design_evaluations import _analysis, design


def _pass_design(database, design, status="passed"):
    with database.connection() as connection:
        connection.execute("UPDATE molecule_observations SET context=json_set(context,"
            "'$.workflow_role','design_generation','$.candidate.raw_metadata.validation_acceptance',json(?)) WHERE id=?",
            (json.dumps({"status": status, "policy": {"version": 1}}), design["observation"]))


def _inventory(database, sequence="ACD"):
    index = build_browser_index(database, USER, include_related=True)
    subject = next(row for row in component_subjects(index, role="designs") if row["sequence"] == sequence)
    rows = subject_results(database, index, subject, include_metrics=True)
    return index, rows, scientific_outcomes(index, rows)


def test_analyses_never_add_outcomes_or_replace_producing_tool_scores(database, design):
    _pass_design(database, design)
    index, before_rows, before = _inventory(database)
    assert len(filter_outcomes(before)) == 1
    assert before[0]["originals"][0]["metrics"] == {"confidence": 0.9}
    for tool in ("ipsae", "pep_patch", "biotite_interface"):
        _analysis(database, design, tool, run=tool, metrics={"confidence": 0.1, tool: 0.7})
    index, after_rows, after = _inventory(database)
    assert len(after_rows) == len(before_rows) + 3
    assert [row["id"] for row in after] == [row["id"] for row in before]
    assert after[0]["origin_tools"] == ["predictor"]
    assert after[0]["originals"][0]["metrics"] == {"confidence": 0.9}
    assert {row["tool"] for row in after[0]["attached"]} == {"ipsae", "pep_patch", "biotite_interface"}


@pytest.mark.parametrize("status,expected", [("passed", "Passed"), ("pending", "Needs review"), ("rejected", None)])
def test_only_explicitly_passing_designs_are_default_retained(database, design, status, expected):
    _pass_design(database, design, status)
    index, rows, outcomes = _inventory(database)
    if expected is None:
        assert not outcomes
        assert rows
    else:
        assert outcomes[0]["eligibility"] == expected
        assert bool(filter_outcomes(outcomes)) == (expected == "Passed")


def test_manual_unknown_and_rejected_are_not_silently_promoted():
    pending = outcome_receipt({"saved_selection": True, "acceptance_status": "unknown"})
    assert pending["saved_selection"] and pending["state"] == "Needs review"
    assert outcome_receipt({"status": "completed", "score": 1})["state"] == "Needs review"
    assert outcome_receipt({"saved_selection": True, "acceptance_status": "failed"})["state"] == "Rejected"
    final = {"candidate": {"raw_metadata": {"validation_acceptance": {"status": "passed"}}}}
    assert outcome_receipt(final, native_flags={"pass_filters": False})["state"] == "Passed"
    assert outcome_receipt({}, native_flags={"pass_filters": "true"})["state"] == "Needs review"
    final["candidate"]["raw_metadata"]["distributed_campaign"] = {"validation_required": True, "validated_acceptance": "rejected"}
    assert outcome_receipt(final)["state"] == "Rejected"


def test_exact_final_acceptance_can_qualify_design_but_digest_conflicts_cannot(database, design):
    _analysis(database, design, "validation_acceptance", context={"candidate": {
        "stage": "analysis", "raw_metadata": {"validation_acceptance": {"status": "passed"}}}})
    index, rows, outcomes = _inventory(database)
    assert len(filter_outcomes(outcomes)) == 1
    with database.connection() as connection:
        connection.execute("UPDATE molecule_observations SET context=json_set(context,"
            "'$.source_structure_digest','invalid','$.structure_digest','invalid') WHERE namespace='validation_acceptance'")
    index, rows, outcomes = _inventory(database)
    assert not filter_outcomes(outcomes)


def test_conflicting_final_decisions_are_review_not_arbitrary_latest(database, design):
    for status in ("passed", "rejected"):
        _analysis(database, design, "validation_acceptance", run=status, context={"candidate": {
            "stage": "analysis", "raw_metadata": {"validation_acceptance": {"status": status}}}})
    index, rows, outcomes = _inventory(database)
    assert outcomes[0]["eligibility"] == "Needs review"
    assert not filter_outcomes(outcomes)


def test_same_chemical_design_across_runs_has_one_row_but_separate_original_scores(database, design):
    _pass_design(database, design)
    second = _import(database, records=[design["record"]], run_id="other-run")
    database.add_observation(second["group_id"], second["molecule_ids"][0], "other-generator", "other-run", {"confidence": 0.4},
        context={"task_group": "design", "target_molecule_id": design["target_id"], "workflow_role": "design_generation",
                 "candidate": {"raw_metadata": {"validation_acceptance": {"status": "passed"}}}})
    index, rows, outcomes = _inventory(database)
    assert len(outcomes) == 1 and outcomes[0]["source_runs"] == 2
    assert sorted(row["metrics"]["confidence"] for row in outcomes[0]["originals"]) == [0.4, 0.9]
    changed = _import(database, records=[_record("ACF", metadata={"design_target_molecule_id": design["target_id"]})])
    index = build_browser_index(database, USER, include_related=True)
    target = next(row for row in component_subjects(index, role="targets") if row["sequence"] == "EFG")
    outcomes = scientific_outcomes(index, subject_results(database, index, target))
    assert len(outcomes) == 2
    assert changed["molecule_ids"][0] in {row["molecule_id"] for row in outcomes}


def test_workflow_purpose_not_predictor_name_distinguishes_design_refold_and_standalone():
    source = {"stage": "complex_refolding", "source_tool": "boltz2", "raw_metadata": {
        "source_candidate": {"stage": "generation.backbone_sequence"}}}
    assert workflow_role({}, {}, source, mode="evaluation", prediction=True) == "independent_prediction"
    source["raw_metadata"]["pipeline_refold_input"] = {"source_candidate_id": "design-1"}
    assert workflow_role({}, {}, source, mode="evaluation", prediction=True) == "design_validation"
    assert workflow_role({"workflow_role": "independent_prediction"}, {}, source, mode="evaluation", prediction=True) == "independent_prediction"
    assert workflow_role({}, {}, source, mode="generated", prediction=True) == "design_generation"
    assert workflow_role({}, {}, source, mode="evaluation", analysis_only=True) == "post_analysis"


def test_refolding_within_design_is_child_but_independent_prediction_is_an_outcome(database, design):
    _pass_design(database, design)
    index, rows, outcomes = _inventory(database)
    original = next(row for row in rows if row["observation_id"] == design["observation"])
    refold = {**deepcopy(original), "id": "refold", "evaluation_id": "refold", "run_key": "refold-run",
              "receipt": {"role": "design_validation", "state": "Needs review"}}
    prediction = {**deepcopy(original), "id": "prediction", "evaluation_id": "prediction", "run_key": "prediction-run",
                  "receipt": {"role": "independent_prediction", "state": "Needs review"}}
    index["runs"]["prediction-run"] = {"task_group": "refolding-validation"}
    index["runs"]["refold-run"] = {"task_group": "refolding-validation"}
    outcomes = scientific_outcomes(index, [original, refold, prediction])
    assert Counter(row["outcome_kind"] for row in outcomes) == {"Design": 1, "Structure prediction": 1}
    assert len(next(row for row in outcomes if row["outcome_kind"] == "Design")["attached"]) == 1


def test_binding_cases_preserve_target_revision_but_collect_tools_and_poses(database, design):
    _pass_design(database, design)
    index, rows, outcomes = _inventory(database)
    original = next(row for row in rows if row["observation_id"] == design["observation"])
    poses = []
    for offset, (target, revision, tool) in enumerate(((design["target_id"], "first", "vina"),
            (design["target_id"], "first", "diffdock"), (design["target_id"], "second", "vina"), ("other-target", "first", "vina"))):
        poses.append({**deepcopy(original), "id": str(offset), "tool": tool, "model": str(offset),
                      "receipt": {"role": "docking", "state": "Needs review"},
                      "evaluation": {**original["evaluation"], "target_id": target, "target_revision": revision}})
    outcomes = scientific_outcomes(index, poses)
    assert len(outcomes) == 3
    combined = next(row for row in outcomes if len(row["originals"]) == 2)
    assert combined["origin_tools"] == ["diffdock", "vina"]
    assert len(filter_outcomes(outcomes)) == 3


def test_missing_identity_or_target_never_combines_unrelated_records(database, design):
    _pass_design(database, design)
    index, rows, outcomes = _inventory(database)
    original = next(row for row in rows if row["observation_id"] == design["observation"])
    unknown = deepcopy(original)
    unknown["molecule_id"] = "missing-identity"
    index["molecules"]["missing-identity"] = {**deepcopy(index["molecules"][original["molecule_id"]]), "components": []}
    assert len(scientific_outcomes(index, [original, unknown])) == 2
    unknown = deepcopy(original)
    unknown["run_key"] = "different-run"
    unknown["evaluation"]["target_id"] = ""
    unknown["receipt"]["target_id"] = ""
    index["molecules"][original["molecule_id"]]["metadata_links_by_group"] = {}
    original["evaluation"]["target_id"] = ""
    original["receipt"]["target_id"] = ""
    assert len(scientific_outcomes(index, [original, unknown])) == 2


def test_outcome_cache_refreshes_after_analysis_deletion_without_changing_identity(database, design):
    _pass_design(database, design)
    observation = _analysis(database, design)
    cache = BrowserCache(database.path)
    index = cache.get_index(database, USER, scope="all", include_related=True)
    subject = next(row for row in component_subjects(index, role="designs") if row["sequence"] == "ACD")
    before = cache.get_outcomes(database, index, subject)
    assert cache.get_outcomes(database, index, subject) is before
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_observations WHERE id=?", (observation,))
    index = cache.get_index(database, USER, scope="all", include_related=True)
    after = cache.get_outcomes(database, index, subject)
    assert after is not before and after[0]["id"] == before[0]["id"]
    assert not after[0]["attached"]
    cache.close()


def test_default_outcome_ui_keeps_original_and_three_analyses_on_this_molecule_tab(database, design, monkeypatch):
    from mn_molecular_design_studio import database_ui, molecule_table

    table_rows = []

    def capture_table(rows, *args, **kwargs):
        table_rows.extend(rows)

    monkeypatch.setattr(molecule_table, "render_table", capture_table)

    _pass_design(database, design)
    for tool in ("ipsae", "pep_patch", "biotite_interface"):
        _analysis(database, design, tool, run=tool, metrics={tool + "_measurement": 0.7})
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    app = AppTest.from_string(_app_source(mode="Targets & ligands", outcomes=None), default_timeout=20).run()
    assert not app.exception
    assert app.selectbox(key="database-browser-result-mode").value == "Scientific outcomes"
    assert next(item for item in app.selectbox if item.label == "Database entries").value == "Admitted results"
    assert any("1 outcomes" in item.value for item in app.caption)
    assert len(table_rows) == 1
    assert table_rows[0]["values"]["Admission"] == "Run filter passed"
    assert not {"ipSAE", "Biotite", "PEP-Patch"} & table_rows[0]["values"].keys()
    molecule_tab = next(tab for tab in app.tabs if tab.label == "This molecule")
    inventory = next(frame.value for frame in molecule_tab.dataframe if "Analysis" in frame.value.columns)
    assert inventory.set_index("Analysis")["Status"].to_dict() == {
        "ipSAE": "Recorded", "Biotite": "Recorded", "PEP-Patch": "Recorded"}
    assert any("confidence" in str(frame.value) and "0.9" in str(frame.value) for frame in molecule_tab.dataframe)
    for tool in ("ipsae", "pep_patch", "biotite_interface"):
        assert any(tool + "_measurement" in str(frame.value) for frame in molecule_tab.dataframe)
    selector = next(item for item in app.selectbox if item.label == "Evaluation")
    index, rows, outcomes = _inventory(database)
    selected = next(row for row in rows if row["tool"] == "ipsae")
    selector.set_value(selected["evaluation_id"]).run()
    assert not app.exception
    assert any("confidence" in str(frame.value) and "0.9" in str(frame.value)
               for tab in app.tabs if tab.label == "This molecule" for frame in tab.dataframe)
    app.selectbox(key="database-browser-result-mode").set_value("Calculation history").run()
    assert not app.exception
    assert any(f"{len(rows)} calculations / records" in item.value for item in app.caption)


def test_protein_sequence_reference_filter_does_not_borrow_a_design_complex(database, design, monkeypatch):
    from mn_molecular_design_studio import database_ui

    _pass_design(database, design)
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    monkeypatch.setattr(database_ui, "render_details", lambda *args, **kwargs: None)
    app = AppTest.from_string(_app_source(mode="Targets & ligands", outcomes=True), default_timeout=20).run()
    assert not app.exception
    selector = app.selectbox(key="database-browser-representation")
    selector.set_value("Recorded structure").run()
    assert not app.exception
    assert any("No matching molecules" in item.value for item in app.info)
    app.selectbox(key="database-browser-representation").set_value("Sequence only").run()
    assert not app.exception
    assert any("1 molecules" in item.value for item in app.caption)


def test_manual_selection_is_separate_from_passes_and_legacy_review(database, design):
    index, rows, outcomes = _inventory(database)
    original = next(row for row in rows if row["observation_id"] == design["observation"])
    original["receipt"].update(saved_selection=True)
    outcomes = scientific_outcomes(index, [original])
    assert outcomes[0]["eligibility"] == "Manually selected"
    assert len(filter_outcomes(outcomes)) == 1
    assert not filter_outcomes(outcomes, eligibility="Passed / recorded")
    assert len(filter_outcomes(outcomes, eligibility="Manually selected")) == 1
    assert len(filter_outcomes(outcomes, eligibility="All retained")) == 1
    original["receipt"]["state"] = "Rejected"
    retained = scientific_outcomes(index, [original])
    assert retained[0]["eligibility"] == "Manually selected"
    assert retained[0]["receipt"]["state"] == "Rejected"
    assert filter_outcomes(retained) == retained
    assert not filter_outcomes(retained, eligibility="Passed / recorded")


def test_default_admission_excludes_unreviewed_design_despite_completed_post_analysis(database, design, monkeypatch):
    from mn_molecular_design_studio import database_ui

    _analysis(database, design, "biotite_interface", context={"status": "completed"})
    index, rows, outcomes = _inventory(database)
    assert outcomes and not filter_outcomes(outcomes)
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    app = AppTest.from_string(_app_source(mode="Targets & ligands", outcomes=None), default_timeout=20).run()
    assert not app.exception
    assert any("No admitted results" in item.value for item in app.info)
    assert not any("#### Original result scores" == item.value for item in app.markdown)


def test_identical_native_names_select_the_correct_design_and_its_analyses(database, design, monkeypatch):
    from mn_molecular_design_studio import database_ui

    _pass_design(database, design)
    _analysis(database, design)
    second = _import(database, run_id="other-run", records=[_record("ACF", name=design["record"]["name"],
        components=[{"type": "protein", "chain": "B", "sequence": "ACF", "role": "binder"}],
        metadata={"design_target_molecule_id": design["target_id"], "validation_acceptance": {"status": "passed"}})])
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    app = AppTest.from_string(_app_source(mode="Targets & ligands", outcomes=None), default_timeout=20).run()
    assert not app.exception
    picker = next(item for item in app.selectbox if item.label == "Inspect result")
    assert len(picker.options) == len(set(picker.options)) == 2
    cache = app.session_state["database-browser-cache"]
    outcomes = next(iter(cache.outcomes.values()))
    selected = next(row for row in outcomes if row["molecule_id"] == second["molecule_ids"][0])
    picker.set_value(selected["id"]).run()
    assert not app.exception
    assert app.session_state["database-focus"] == second["molecule_ids"][0]
    selected = next(row for row in outcomes if row["molecule_id"] == design["molecule_id"])
    next(item for item in app.selectbox if item.label == "Inspect result").set_value(selected["id"]).run()
    assert not app.exception
    assert app.session_state["database-focus"] == design["molecule_id"]
    assert any("Recorded" in str(frame.value) for tab in app.tabs if tab.label == "This molecule" for frame in tab.dataframe)


def test_explicit_later_selection_shares_one_entry_with_the_passing_design(database, design):
    _pass_design(database, design)
    selected = _import(database, kind="curated", run_id="later-selection", records=[{
        **design["record"], "metadata": {**design["record"]["metadata"], "workflow_role": "design_generation"}}])
    index, rows, outcomes = _inventory(database)
    admitted = filter_outcomes(outcomes)
    assert len(admitted) == 1
    assert admitted[0]["source_runs"] == 2
    assert {row["eligibility"] for row in admitted[0]["originals"]} == {"Passed", "Manually selected"}
    assert selected["group_id"] in {row["group_id"] for row in admitted[0]["originals"]}


def test_curated_candidate_import_records_later_selection_without_inventing_a_pass(database):
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from test_molecular_database_policy import _make_run, _pdb_acd

    run = _make_run(database, "candidate-import", "selected-later", owner=USER["id"], project="",
        tool="candidate_import", workflow="candidate_import", candidates=[{
            "candidate_id": "chosen", "stage": "generation.backbone_sequence", "binder_sequence": "ACD",
            "binder_chains": ["A"], "binder_pdb": "binder.pdb", "metrics": {"confidence": 0.4},
            "raw_metadata": {"validation_acceptance": {"status": "rejected"}}}])
    (run / "binder.pdb").write_text(_pdb_acd())
    assert ingest_run(run, database=database) == 1
    index, rows, outcomes = _inventory(database)
    admitted = filter_outcomes(outcomes)
    assert len(admitted) == 1
    assert admitted[0]["eligibility"] == "Manually selected"
    assert admitted[0]["receipt"]["state"] == "Rejected"
    assert not filter_outcomes(outcomes, eligibility="Passed / recorded")


def test_saved_library_default_exposes_legacy_design_without_claiming_a_pass(database, design):
    index, rows, outcomes = _inventory(database)
    assert len(filter_outcomes(outcomes, eligibility="Saved in library")) == 1
    assert outcomes[0]["eligibility"] == "Needs review"
    assert not filter_outcomes(outcomes, eligibility="Passed / recorded")
    _pass_design(database, design, "rejected")
    index, rows, rejected = _inventory(database)
    assert not filter_outcomes(rejected, eligibility="Saved in library")
    assert rows


def test_outcome_analysis_inventory_excludes_other_targets_and_preserves_failure_status(database, design):
    _analysis(database, design, "ipsae", context={"status": "completed"})
    _analysis(database, design, "biotite_interface", run="failed", context={"status": "failed"})
    _analysis(database, design, "pep_patch", run="wrong-target", context={"target_molecule_id": "wrong-target", "status": "completed"})
    index, rows, outcomes = _inventory(database)
    inventory = outcomes[0]["analysis_inventory"]
    assert inventory["ipSAE"]["same_structure"] == 1
    assert inventory["Biotite"]["failed"] == 1
    assert "PEP-Patch" not in inventory


@pytest.mark.parametrize("refold,expected", [(False, "independent_prediction"), (True, "design_validation")])
def test_new_ingestion_persists_workflow_purpose_without_changing_collection_policy(database, refold, expected):
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from test_molecular_database_policy import _make_run, _pdb_acd

    raw = {"pipeline_refold_input": {"source_candidate_id": "source-design"}} if refold else {}
    run = _make_run(database, "optimized-prediction", expected, owner=USER["id"], project="", tool="optimized_prediction_boltz2",
        candidates=[{"candidate_id": "candidate", "stage": "complex_refolding", "binder_sequence": "ACD",
                     "binder_chains": ["A"], "complex_pdb": "complex.pdb", "raw_metadata": raw, "metrics": {"confidence": 0.9}}])
    (run / "complex.pdb").write_text(_pdb_acd())
    assert ingest_run(run, database=database) == 1
    with database.connection() as connection:
        purposes = [row[0] for row in connection.execute("SELECT json_extract(context,'$.workflow_role') FROM molecule_observations")]
        assert purposes == [expected]
    assert ingest_run(run, database=database) == 0
    metadata = json.loads((run / "metadata.json").read_text())
    metadata.pop("database_collection_policy")
    (run / "metadata.json").write_text(json.dumps(metadata))
    assert ingest_run(run, database=database, force=True) == 0


def test_source_acceptance_is_never_borrowed_from_another_model_or_member(database, design):
    with database.connection() as connection:
        connection.execute("UPDATE molecule_members SET metadata=json_set(metadata,'$.validation_acceptance.status','passed') WHERE molecule_id=?",
                           (design["molecule_id"],))
        connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.candidate_id','candidate','$.model_id',0) WHERE id=?",
                           (design["observation"],))
        connection.execute("INSERT INTO molecule_sources (source_key,molecule_id,group_id,provenance) VALUES (?,?,?,?)",
            ("other-model", design["molecule_id"], design["source_group_id"], json.dumps({"candidate_id": "candidate", "model_id": 1,
             "task_group": "design", "run_id": "run-1", "acceptance_decision": {"status": "accepted"}})))
    index, rows, outcomes = _inventory(database)
    selected = next(row for row in rows if row["observation_id"] == design["observation"])
    assert selected["receipt"]["state"] == "Needs review"
    assert not filter_outcomes(outcomes)


def test_complex_roles_are_scoped_to_the_original_collection(database, design):
    from mn_molecular_design_studio.database_design_evaluations import connected_evaluations, design_subject

    components = [{"type": "protein", "chain": "A", "sequence": "EFG", "role": "target"},
                  {"type": "protein", "chain": "B", "sequence": "ACD", "role": "binder"}]
    first = _import(database, run_id="complex", records=[_record(components=components, metadata={"design_target_molecule_id": design["target_id"]})])
    inverse = [{**component, "role": "binder" if component["role"] == "target" else "target"} for component in components]
    second = _import(database, run_id="inverse", records=[_record(components=inverse)])
    assert first["molecule_ids"] == second["molecule_ids"]
    complex_design = {**design, "molecule_id": first["molecule_ids"][0], "group_id": first["group_id"]}
    included = _analysis(database, complex_design)
    excluded = _analysis(database, {**complex_design, "group_id": second["group_id"]}, run="inverse")
    index, rows, outcomes = _inventory(database)
    selected = next(row for row in rows if row["observation_id"] == design["observation"])
    cohort = design_subject(index, selected)
    evaluations, unresolved = connected_evaluations(index, selected, subject_results(database, index, cohort, include_metrics=True))
    assert included in {row["observation_id"] for row in evaluations}
    assert excluded not in {row["observation_id"] for row in evaluations}


def test_intermediate_passing_sequence_is_not_a_retained_design(database, design):
    intermediate = _import(database, run_id="intermediate", kind="intermediate", records=[
        _record("ACF", components=[{"type": "protein", "chain": "B", "sequence": "ACF", "role": "binder"}],
                metadata={"design_target_molecule_id": design["target_id"]})])
    database.add_observation(intermediate["group_id"], intermediate["molecule_ids"][0], "mpnn", "intermediate",
        {"pass_filters": True}, context={"task_group": "design", "workflow_role": "sequence_redesign",
                                        "candidate": {"stage": "sequence_design"}})
    index = build_browser_index(database, USER, collection="all", include_related=True)
    subject = next(row for row in component_subjects(index, role="designs") if row["sequence"] == "ACF")
    rows = subject_results(database, index, subject)
    outcomes = scientific_outcomes(index, rows)
    assert rows and not outcomes


def test_prediction_without_recorded_coordinates_is_review_not_a_successful_outcome(database, design):
    index, rows, outcomes = _inventory(database)
    original = next(row for row in rows if row["observation_id"] == design["observation"])
    original["receipt"]["role"] = "independent_prediction"
    original["structure_ids"] = []
    outcomes = scientific_outcomes(index, [original])
    assert outcomes[0]["eligibility"] == "Needs review"
    assert not filter_outcomes(outcomes)


def test_independent_complex_prediction_is_cofolding_not_docking_or_another_design(database, design):
    index, rows, outcomes = _inventory(database)
    original = next(row for row in rows if row["observation_id"] == design["observation"])
    original["receipt"].update(role="independent_prediction", design_refold=True)
    index["molecules"][original["molecule_id"]]["components"].append({"type": "protein", "chain": "A", "sequence": "EFG"})
    outcomes = scientific_outcomes(index, [original])
    assert len(outcomes) == 1 and outcomes[0]["outcome_kind"] == "Cofolding"
    assert filter_outcomes(outcomes, kind="Cofolding", query="0 Final design") == outcomes


def test_unannotated_legacy_analysis_membership_reuses_unambiguous_design_roles(database, design):
    from mn_molecular_design_studio.database_design_evaluations import connected_evaluations, design_subject

    unannotated = _import(database, kind="evaluation", task="analysis", run_id="legacy", records=[
        _record("ACD", components=[{"type": "protein", "chain": "B", "sequence": "ACD"}])])
    assert unannotated["molecule_ids"] == [design["molecule_id"]]
    observation = _analysis(database, {**design, "group_id": unannotated["group_id"]})
    index, rows, outcomes = _inventory(database)
    selected = next(row for row in rows if row["observation_id"] == design["observation"])
    cohort = design_subject(index, selected)
    linked, unresolved = connected_evaluations(index, selected, subject_results(database, index, cohort))
    assert observation in {row["observation_id"] for row in linked}


def test_analysis_does_not_hide_original_structure_when_native_scores_are_missing(database, design):
    from mn_molecular_design_studio.database_browser import load_result_record

    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_observations WHERE id=?", (design["observation"],))
        connection.execute("UPDATE molecule_members SET metadata=json_set(metadata,'$.validation_acceptance.status','passed') WHERE molecule_id=? AND group_id=?",
                           (design["molecule_id"], design["source_group_id"]))
    _analysis(database, design)
    index, rows, outcomes = _inventory(database)
    assert len(filter_outcomes(outcomes)) == 1
    selected = filter_outcomes(outcomes)[0]
    assert selected["structure_ids"] == [design["structure"]]
    record = load_result_record(database, USER, selected)
    assert [structure["id"] for structure in record["structures"]] == [design["structure"]]
    assert not record["observations"] and not record["scores"]
