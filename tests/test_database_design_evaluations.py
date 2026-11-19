from __future__ import annotations

import hashlib
import json
import sqlite3
from contextlib import contextmanager

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects, load_result_record, subject_results
from mn_molecular_design_studio.database_browser_cache import BrowserCache
from mn_molecular_design_studio.database_browser_labels import tool_label
from mn_molecular_design_studio.database_design_evaluations import (
    analysis_results, connected_evaluations, design_subject, evaluation_label,
)
from test_database_browser import OTHER, USER, _app_source, _import, _record, database
from test_molecular_database_policy import _pdb_efg_acd_complex


@pytest.fixture
def design(database):
    target = _import(database, name="Target", run_id="", method="user_uploads", kind="reference",
                     records=[_record("EFG", name="Target", components=[
                         {"type": "protein", "chain": "A", "sequence": "EFG", "role": "target"}])])
    target_id = target["molecule_ids"][0]
    record = _record("ACD", metadata={"design_target_molecule_id": target_id}, components=[
        {"type": "protein", "chain": "B", "sequence": "ACD", "role": "binder"}])
    generated = _import(database, name="Final design", records=[record])
    molecule_id = generated["molecule_ids"][0]
    content = _pdb_efg_acd_complex().encode()
    structure = database.attach_structure(generated["group_id"], molecule_id,
                                          {"name": "complex.pdb", "content": content}, context={"run_id": "run-1"})
    digest = hashlib.sha256(content).hexdigest()
    observation = database.add_observation(generated["group_id"], molecule_id, "predictor", "run-1", {"confidence": 0.9},
        context={"task_group": "design", "target_molecule_id": target_id, "representation_id": structure,
                 "structure_digest": digest, "model_id": "0", "candidate": {"stage": "complex_refolding", "source_tool": "predictor"}})
    evaluated = _import(database, name="Evaluated structures", task="analysis", run_id="post", kind="evaluation", records=[record])
    return {"molecule_id": molecule_id, "group_id": evaluated["group_id"], "source_group_id": generated["group_id"],
            "structure": structure, "digest": digest, "target_id": target_id, "observation": observation, "record": record}


def _analysis(database, design, tool="biotite_interface", *, run="post", context=None, metrics=None):
    return database.add_observation(design["group_id"], design["molecule_id"], tool, run, metrics or {"score": 0.7},
        context={"task_group": "analysis", "target_molecule_id": design["target_id"],
                 "source_molecule_id": design["molecule_id"], "source_structure_id": design["structure"],
                 "representation_id": design["structure"], "source_structure_digest": design["digest"],
                 "structure_digest": design["digest"], "model_id": "0",
                 "candidate": {"source_tool": tool, "stage": "analysis", "raw_metadata": {"engine_version": "v2"}},
                 **(context or {})})


def _connected(database, design):
    index = build_browser_index(database, USER, include_related=True)
    subject = next(row for row in component_subjects(index, role="designs") if row["sequence"] == "ACD")
    selected = next(row for row in subject_results(database, index, subject) if row["observation_id"] == design["observation"])
    cohort = design_subject(index, selected)
    rows = subject_results(database, index, cohort, include_metrics=True)
    connected, unresolved = connected_evaluations(index, selected, rows)
    return index, selected, connected, unresolved


def test_post_analyses_are_listed_together_without_merging_scores(database, design):
    for tool, metrics in (("ipsae", {"ipsae": 0.8}), ("pep_patch", {"patch_score": 0.4}), ("biotite_interface", {"contacts": 12})):
        _analysis(database, design, tool, run=tool, metrics=metrics)
    index, selected, connected, unresolved = _connected(database, design)
    assert len(connected) == 4 and unresolved == 0
    assert len(analysis_results(connected)) == 3
    assert {row["structure_relation"] for row in connected} == {"Selected structure"}
    assert "v2" in evaluation_label(next(row for row in connected if row["tool"] == "ipsae"))
    for row in connected:
        record = load_result_record(database, USER, row)
        assert len(record["observations"]) == 1
        assert {entry["id"] for entry in record["structures"]} == {design["structure"]}
        assert all(source["observation_id"] == row["observation_id"] for source in record["score_sources"].values())


def test_binder_and_complex_share_analysis_but_redesigned_sequence_does_not(database, design):
    complex_record = _record(name="Complex", components=[
        {"type": "protein", "chain": "A", "sequence": "EFG", "role": "target"},
        {"type": "protein", "chain": "B", "sequence": "ACD", "role": "binder"}],
        metadata={"design_target_molecule_id": design["target_id"]})
    imported = _import(database, records=[complex_record], run_id="complex", kind="evaluation")
    complex_design = {**design, "molecule_id": imported["molecule_ids"][0], "group_id": imported["group_id"]}
    complex_analysis = _analysis(database, complex_design)
    changed = _import(database, run_id="changed", records=[_record("ACF", metadata={"design_target_molecule_id": design["target_id"]})])
    changed_analysis = _analysis(database, {**design, "molecule_id": changed["molecule_ids"][0], "group_id": changed["group_id"]})
    index, selected, connected, unresolved = _connected(database, design)
    assert complex_design["molecule_id"] != design["molecule_id"]
    assert complex_analysis in {row["observation_id"] for row in connected}
    assert changed_analysis not in {row["observation_id"] for row in connected}


def test_another_target_is_excluded_even_for_the_same_sequence_and_structure(database, design):
    other = _analysis(database, design, context={"target_molecule_id": "other-target"})
    index, selected, connected, unresolved = _connected(database, design)
    assert other not in {row["observation_id"] for row in connected}


def test_other_models_and_prepared_structures_are_labelled_not_merged(database, design):
    other = _analysis(database, design, run="other-model", context={
        "source_structure_id": "other-structure", "representation_id": "other-structure", "model_id": "2",
        "source_structure_digest": "b" * 64, "structure_digest": "b" * 64})
    derived = _analysis(database, design, run="prepared", context={"calculation_structure_digest": "c" * 64})
    index, selected, connected, unresolved = _connected(database, design)
    by_observation = {row["observation_id"]: row for row in connected}
    assert by_observation[other]["structure_relation"] == "Other structure of this design"
    assert by_observation[derived]["structure_relation"] == "Derived from selected structure"


@pytest.mark.parametrize("context", [
    {"source_structure_digest": "bad-digest"},
    {"source_molecule_id": "unrelated-molecule"},
])
def test_conflicting_source_provenance_is_flagged_not_attached(database, design, context):
    observation = _analysis(database, design, context=context)
    index, selected, connected, unresolved = _connected(database, design)
    assert observation not in {row["observation_id"] for row in connected}
    assert unresolved == 1


def test_unknown_target_without_exact_structure_stays_unresolved_but_sequence_qc_is_available(database, design):
    context = {"target_molecule_id": "", "representation_id": "", "source_structure_id": "",
               "source_structure_digest": "", "structure_digest": ""}
    unresolved_id = _analysis(database, design, context=context)
    sequence_id = _analysis(database, design, "sequence_qc", run="sequence", context={**context, "source_reference_kind": "sequence"})
    index, selected, connected, unresolved = _connected(database, design)
    assert unresolved == 1
    assert unresolved_id not in {row["observation_id"] for row in connected}
    sequence = next(row for row in connected if row["observation_id"] == sequence_id)
    assert sequence["structure_relation"] == "Sequence only"
    assert not sequence["structure_ids"]


def test_versions_and_settings_remain_distinct(database, design):
    first = _analysis(database, design, "ipsae", context={"settings": {"pae_cutoff": 10}})
    second = _analysis(database, design, "ipsae", run="new-version", context={"settings": {"pae_cutoff": 15},
        "candidate": {"source_tool": "ipsae", "stage": "analysis", "raw_metadata": {"engine_version": "v3"}}})
    index, selected, connected, unresolved = _connected(database, design)
    rows = [row for row in connected if row["observation_id"] in {first, second}]
    assert len(rows) == 2
    assert len({row["evaluation_id"] for row in rows}) == 2
    assert len({row["evaluation"]["parameters_digest"] for row in rows}) == 2
    assert {row["evaluation"]["engine_version"] for row in rows} == {"v2", "v3"}


@pytest.mark.parametrize("saved", [True, False])
def test_compact_browser_context_matches_raw_scientific_receipt_without_writing(database, design, monkeypatch, saved):
    from mn_molecular_design_studio import evaluations

    with database.connection() as connection:
        context = json.loads(connection.execute("SELECT context FROM molecule_observations WHERE id=?", (design["observation"],)).fetchone()[0])
        context.update(saved_selection=saved, workflow_id="workflow-1", workflow_role="design_generation",
            pipeline_refold_input={"source_candidate_id": "design-1"},
            acceptance_decision={"status": "accepted", "evidence": {"native_filter_metrics": {"pass_filters": True}}},
            acceptance_stage="Final selection")
        connection.execute("UPDATE molecule_observations SET context=? WHERE id=?", (json.dumps(context), design["observation"]))
        while evaluations.rebuild_missing(connection):
            pass
    index = build_browser_index(database, USER, include_related=True)
    subject = next(row for row in component_subjects(index, role="designs") if row["sequence"] == "ACD")

    @contextmanager
    def read_only():
        connection = sqlite3.connect(database.path.as_uri() + "?mode=ro", uri=True)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only=ON")
        try:
            yield connection
            assert connection.total_changes == 0
        finally:
            connection.close()

    monkeypatch.setattr(database, "connection", read_only)
    projected = subject_results(database, index, subject, include_metrics=True)
    monkeypatch.setattr(evaluations, "PROJECTION_VERSION", "future-version")
    raw = subject_results(database, index, subject, include_metrics=True)
    assert projected == raw
    selected = next(row for row in projected if row["observation_id"] == design["observation"])
    assert selected["receipt"]["saved_selection"] is saved
    assert selected["receipt"]["workflow"] == "workflow-1"
    assert selected["receipt"]["design_refold"] is True


def test_ligand_stereochemistry_and_binder_stoichiometry_do_not_cross_link(database):
    first = _import(database, run_id="first", records=[_record(name="Ligand", components=[
        {"type": "small_molecule", "chain": "L", "smiles": "F[C@H](Cl)Br"}])])
    second = _import(database, run_id="second", records=[_record(name="Other isomer", components=[
        {"type": "small_molecule", "chain": "L", "smiles": "F[C@@H](Cl)Br"}])])
    binder = _import(database, run_id="binder", records=[_record("ACD", components=[
        {"type": "protein", "chain": "B", "sequence": "ACD", "role": "binder"}])])
    dimer = _import(database, run_id="dimer", records=[_record(components=[
        {"type": "protein", "chain": chain, "sequence": "ACD", "role": "binder"} for chain in ("B", "C")])])
    index = build_browser_index(database, USER, include_related=True)
    ligand_cohort = design_subject(index, {"molecule_id": first["molecule_ids"][0]})
    assert not ligand_cohort["molecule_ids"] & set(second["molecule_ids"])
    binder_cohort = design_subject(index, {"molecule_id": binder["molecule_ids"][0]})
    assert not binder_cohort["molecule_ids"] & set(dimer["molecule_ids"])


def test_private_analyses_never_enter_design_inventory(database, design):
    private = _import(database, user=OTHER, records=[design["record"]], kind="evaluation")
    secret = _analysis(database, {**design, "group_id": private["group_id"]}, metrics={"secret": 987})
    index, selected, connected, unresolved = _connected(database, design)
    assert secret not in {row["observation_id"] for row in connected}
    assert "secret" not in json.dumps(connected, default=sorted)


def test_cache_distinguishes_metric_inventory_and_reloads_after_analysis_deletion(database, design):
    observation = _analysis(database, design)
    cache = BrowserCache(database.path)
    index = cache.get_index(database, USER, scope="all", include_related=True)
    subject = next(row for row in component_subjects(index, role="designs") if row["sequence"] == "ACD")
    selected = next(row for row in cache.get_results(database, index, subject) if row["observation_id"] == design["observation"])
    cohort = design_subject(index, selected)
    assert not any(row["metrics"] for row in cache.get_results(database, index, cohort))
    with_metrics = cache.get_results(database, index, cohort, include_metrics=True)
    assert any(row["metrics"] for row in with_metrics)
    assert cache.get_results(database, index, cohort, include_metrics=True) is with_metrics
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_observations WHERE id=?", (observation,))
    refreshed = cache.get_index(database, USER, scope="all", include_related=True)
    assert refreshed is not index
    rows = cache.get_results(database, refreshed, design_subject(refreshed, selected), include_metrics=True)
    assert observation not in {row["observation_id"] for row in rows}
    cache.close()


def test_revoked_project_analyses_disappear_from_cached_inventory(database, design):
    member = {**USER, "_cluster_project_ids": ["shared-project"]}
    shared = _import(database, user=member, project_id="shared-project", records=[design["record"]], kind="evaluation")
    observation = _analysis(database, {**design, "group_id": shared["group_id"]})
    cache = BrowserCache(database.path)
    index = cache.get_index(database, member, scope="all", include_related=True)
    subject = next(row for row in component_subjects(index, role="designs") if row["sequence"] == "ACD")
    selected = next(row for row in cache.get_results(database, index, subject) if row["observation_id"] == design["observation"])
    before = cache.get_results(database, index, design_subject(index, selected), include_metrics=True)
    assert observation in {row["observation_id"] for row in before}
    revoked = cache.get_index(database, USER, scope="all", include_related=True)
    after = cache.get_results(database, revoked, design_subject(revoked, selected), include_metrics=True)
    assert observation not in {row["observation_id"] for row in after}
    cache.close()


def test_viewer_shows_all_three_post_analyses_and_switches_exact_evaluation(database, design, monkeypatch):
    from mn_molecular_design_studio import database_browser, database_ui

    for tool in ("ipsae", "pep_patch", "biotite_interface"):
        _analysis(database, design, tool, run=tool, metrics={tool + "_measurement": 0.7})
    index, selected, connected, unresolved = _connected(database, design)
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    loaded = []
    original = database_browser.load_result_record

    def capture(*args, **kwargs):
        record = original(*args, **kwargs)
        if record:
            loaded.append(record)
        return record

    monkeypatch.setattr(database_browser, "load_result_record", capture)
    app = AppTest.from_string(_app_source(), default_timeout=20).run()
    assert not app.exception
    picker = next(item for item in app.selectbox if item.label == "Inspect result")
    picker.set_value(selected["id"]).run()
    assert not app.exception
    assert "Design analyses (3)" in [tab.label for tab in app.tabs]
    for tool in ("ipsae", "pep_patch", "biotite_interface"):
        assert any(item.label.startswith(tool_label(tool) + " · ") for item in app.expander)
        assert any(tool + "_measurement" in str(frame.value) for frame in app.dataframe)
    selector = next(item for item in app.selectbox if item.label == "Evaluation")
    assert len(selector.options) == 4
    analysis = next(row for row in connected if row["tool"] == "ipsae")
    selector.set_value(analysis["evaluation_id"]).run()
    assert not app.exception
    assert loaded[-1]["selected_evaluation"]["id"] == analysis["evaluation_id"]
    assert len(loaded[-1]["observations"]) == 1
    assert all(source["namespace"] == "ipsae" for source in loaded[-1]["score_sources"].values())


def test_viewer_handles_deleted_analysis_and_original_structure_entries(database, design, monkeypatch):
    from mn_molecular_design_studio import database_ui

    observation = _analysis(database, design, "ipsae")
    index, selected, connected, unresolved = _connected(database, design)
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    app = AppTest.from_string(_app_source(), default_timeout=20).run()
    picker = next(item for item in app.selectbox if item.label == "Inspect result")
    picker.set_value(selected["id"]).run()
    analysis = next(row for row in connected if row["observation_id"] == observation)
    next(item for item in app.selectbox if item.label == "Evaluation").set_value(analysis["evaluation_id"]).run()
    assert not app.exception
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_observations WHERE id=?", (observation,))
    app.run()
    assert not app.exception
    assert "Design analyses (0)" in [tab.label for tab in app.tabs]
    selector = next(item for item in app.selectbox if item.label == "Evaluation")
    assert selector.value == selected["evaluation_id"]


def test_structure_entry_retains_original_view_and_exposes_connected_analyses(database, design, monkeypatch):
    from mn_molecular_design_studio import database_ui

    _analysis(database, design, "ipsae")
    structure_id = database.attach_structure(design["source_group_id"], design["molecule_id"],
        {"name": "other.pdb", "label": "Alternate view", "database_structure_label": "alternate",
         "content": _pdb_efg_acd_complex().encode()},
        context={"run_id": "run-1"})
    index = build_browser_index(database, USER, include_related=True)
    subject = next(row for row in component_subjects(index, role="designs") if row["sequence"] == "ACD")
    rows = subject_results(database, index, subject)
    structure_row = next(row for row in rows if row["kind"] == "Structure" and structure_id in row["structure_ids"])
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    app = AppTest.from_string(_app_source(), default_timeout=20).run()
    next(item for item in app.selectbox if item.label == "Inspect result").set_value(structure_row["id"]).run()
    assert not app.exception
    assert "Design analyses (1)" in [tab.label for tab in app.tabs]
    selector = next(item for item in app.selectbox if item.label == "Evaluation")
    assert selector.value == ""
    assert len(selector.options) == 3
