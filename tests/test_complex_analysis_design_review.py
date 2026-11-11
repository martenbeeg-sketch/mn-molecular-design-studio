"""The design workbench joins exact evidence without loading structures or Results."""
import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_design_review import (
    design_review_page, review_measurements, short_design_label, remember_review_selection,
)

OWNER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def designs(tmp_path):
    db = MolecularDatabase(tmp_path / "db")
    created = db.import_records(OWNER, molecule_imports.parse_smiles("CCO example"), name="Saved designs",
        origin={"tool": "rfd4_proteina", "type": "job", "run_id": "design-source", "campaign_id": "mixed-campaign"})
    group, molecule = created["group_id"], created["molecule_ids"][0]
    models = []
    for index, engine in enumerate(("rfd4_proteina", "bindcraft2")):
        candidate = f"rfd4_00000{index+1}_binder_rep0" if index == 0 else "bindcraft_native_2"
        context = {"candidate_id": candidate, "model_id": f"model-{index}", "tool": engine,
                   "run_id": "design-source", "workstation": "PC004033", "target_id": "target",
                   "binder_chains": ["B"], "target_chains": ["A"]}
        model = db.attach_structure(group, molecule, {"name": f"model{index}.cif", "content": f"data_{index}\n".encode()}, context=context)
        models.append(model)
        raw = {"validation_acceptance": {"status": "accepted", "policy": {"schema": "v1"}}} if index == 0 else {}
        db.add_observation(group, molecule, engine, "design-source", {"native_pass_filters": True, "iptm": 0.7 + index/10},
            context={**context, "representation_id": model, "settings": {"cutoff": 4}, "candidate": {
                "candidate_id": candidate, "model_id": f"model-{index}", "source_tool": engine,
                "binder_sequence": "ACDE", "raw_metadata": raw}})
        db.add_observation(group, molecule, "biotite_interface", "analysis", {"qc_interface_contact_pairs": index},
            context={"representation_id": model, "candidate_id": candidate, "model_id": f"model-{index}", "settings": {"cutoff": 4}})
        if index == 0:
            db.add_observation(group, molecule, "validation_acceptance", "validation", {"accepted": True},
                context={"representation_id": model, "candidate_id": candidate, "model_id": "model-0", "candidate": {
                    "candidate_id": candidate, "model_id": "model-0", "source_tool": "validation_acceptance", "raw_metadata": raw}})
    return db, group, molecule, models


def test_review_combines_methods_per_exact_model_with_no_coordinates(designs, monkeypatch):
    db, group, molecule, models = designs
    monkeypatch.setattr(db, "structure_file", lambda *a: pytest.fail("Metadata review must not fetch coordinates"))
    monkeypatch.setattr(db, "list_molecules", lambda *a, **k: pytest.fail("Metadata review must not hydrate full records"))
    page = design_review_page(db, OWNER, [group])
    assert page["total"] == 2
    assert {sample["structure_id"] for sample in page["items"]} == set(models)
    assert all(sample["molecule_id"] == molecule for sample in page["items"])
    assert sorted(len(sample["associated"]) for sample in page["items"]) == [2, 3]
    assert sum(any(d["stage_id"] == "final" for d in sample["decisions"]) for sample in page["items"]) == 1
    assert all(sample["digest"] and sample["candidate"]["binder_sequence"] == "ACDE" for sample in page["items"])
    series, values = review_measurements(page["items"])
    contacts = next(key for key in series if series[key]["field"] == "qc_interface_contact_pairs")
    assert sorted(v[contacts] for v in values.values()) == [0, 1]
    assert len({descriptor["id"] for sample in page["items"] for _, descriptor in sample["associated"]}) == 5


def test_no_run_or_name_only_join_and_model_decisions_cannot_leak(designs):
    db, group, molecule, models = designs
    db.add_observation(group, molecule, "unrelated", "design-source", {"iptm": 1.0},
                       context={"candidate_id": "unrelated", "model_id": "model-0"})
    db.add_observation(group, molecule, "wrong-model", "analysis", {"iptm": 1.0},
                       context={"representation_id": models[0], "model_id": "model-1"})
    samples = design_review_page(db, OWNER, [group])["items"]
    assert sorted(len(sample["associated"]) for sample in samples) == [2, 3]
    assert samples[0]["label"] != samples[1]["label"]


def test_protocol_conflicts_and_failed_evidence_are_not_silent_numbers(designs):
    db, group, molecule, models = designs
    for run, value, cutoff, status in [("repeat", 3, 4, "completed"), ("different", 9, 5, "completed"), ("failed", 99, 4, "failed")]:
        db.add_observation(group, molecule, "biotite_interface", run, {"qc_interface_contact_pairs": value},
            context={"representation_id": models[0], "settings": {"cutoff": cutoff}, "status": status})
    page = design_review_page(db, OWNER, [group])
    series, values = review_measurements(page["items"])
    assert len([key for key in series if series[key]["field"] == "qc_interface_contact_pairs"]) == 2
    assert 99 not in values[models[0]].values()
    assert 3 not in values[models[0]].values()  # conflicting 0/3 is absent, not arbitrary newest.
    assert 9 in values[models[0]].values()


def test_review_rechecks_acl_and_pagination(designs):
    db, group, _, _ = designs
    with pytest.raises(PermissionError):
        design_review_page(db, OTHER, [group])
    with pytest.raises(PermissionError):
        design_review_page(db, OWNER, [group, "unrelated"])
    with pytest.raises(ValueError):
        design_review_page(db, OWNER, [group], limit=1000)
    first = design_review_page(db, OWNER, [group], limit=1)
    second = design_review_page(db, OWNER, [group], limit=1, offset=first["next_offset"])
    assert first["items"][0]["structure_id"] != second["items"][0]["structure_id"]


def test_digest_aliases_link_evidence_but_not_other_models_or_private_groups(designs):
    db, group, molecule, models = designs
    alias_group = db.create_group(OWNER, "Saved exact evaluation")
    db.add_to_group(OWNER, alias_group, [molecule])
    alias = db.attach_structure(alias_group, molecule, {"name": "same.cif", "content": b"data_0\n"},
        context={"model_id": "model-0", "candidate_id": "rfd4_000001_binder_rep0"})
    db.add_observation(alias_group, molecule, "chemistry", "analysis", {"buried_surface": 42}, context={"representation_id": alias})
    private_record = db.import_records(OTHER, molecule_imports.parse_smiles("CCO independently imported"), name="Private evidence")
    private = private_record["group_id"]
    assert private_record["molecule_ids"] == [molecule]
    private_alias = db.attach_structure(private, molecule, {"name": "private.cif", "content": b"data_0\n"})
    db.add_observation(private, molecule, "private_chemistry", "private", {"buried_surface": 99},
                       context={"representation_id": private_alias})
    page = design_review_page(db, OWNER, [group])
    a = next(s for s in page["items"] if s["structure_id"] == models[0])
    b = next(s for s in page["items"] if s["structure_id"] == models[1])
    assert "chemistry" in {o["namespace"] for o, _ in a["associated"]}
    assert "chemistry" not in {o["namespace"] for o, _ in b["associated"]}
    assert "private_chemistry" not in {o["namespace"] for s in page["items"] for o, _ in s["associated"]}


def test_display_labels_and_selection_do_not_replace_scientific_identity():
    assert short_design_label("rfd4_000010_idr_peptide_binder_rep9") == "Design 10 · sample 9"
    assert short_design_label("rfd4_000001_idr_peptide_binder") == "Design 1"
    state = {}
    assert remember_review_selection(state, "cohort", ["model-1", "model-1", "unrelated"], {"model-1"}) == ["model-1"]
    assert remember_review_selection(state, "cohort", state["design-review-selected-cohort"], set()) == []


def test_workbench_browser_contract_is_metadata_first_and_final_default(designs, monkeypatch):
    from streamlit.testing.v1 import AppTest
    db, group, _, _ = designs
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: db)
    monkeypatch.setattr(db, "structure_file", lambda *a: pytest.fail("No coordinate access during review"))
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_design_review_ui import render_design_review
render_design_review(MolecularDatabase(), {OWNER!r}, {{"group_ids": [{group!r}], "label": "Mixed engines"}})
''').run()
    assert not app.exception
    assert app.selectbox[0].value == "Final accepted designs"
    assert len(app.dataframe[0].value) == 1
    app.selectbox[0].select("All designs").run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 2
    assert set(app.dataframe[0].value["Engine"]) == {"RFD4", "BindCraft2"}
    assert "Target / preparation" not in app.dataframe[0].value.columns


def test_selected_design_reuses_existing_shortlist_and_navigation_contract(designs, monkeypatch):
    import hashlib
    import json
    from streamlit.testing.v1 import AppTest
    db, group, _, models = designs
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: db)
    monkeypatch.setattr(db, "structure_file", lambda *a: pytest.fail("Shortlist must not fetch coordinates"))
    token = hashlib.sha256(json.dumps([OWNER["id"], [group], None], sort_keys=True).encode()).hexdigest()[:12]
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_design_review_ui import render_design_review
st.session_state.setdefault('design-review-selected-{token}', [{models[0]!r}])
render_design_review(MolecularDatabase(), {OWNER!r}, {{"group_ids": [{group!r}], "label": "Mixed engines"}})
''').run()
    assert not app.exception
    assert any("1 selected for review" in item.value for item in app.markdown)
    assert app.button(key="review-inspector-" + token).label == "Open focused Scientific Inspector"
    app.button(key="review-shortlist-" + token).click().run()
    assert not app.exception
    assert any(button.label == "Save selected-design review decisions" for button in app.button)
    assert len(app.get("download_button")) == 3


def _select_rows(app, rows):
    """Deliver a real Streamlit dataframe selection event, including its callback."""
    import json
    states = app._tree.get_widget_states()
    widget = states.widgets.add()
    widget.id = app.dataframe[0].proto.id
    widget.string_value = json.dumps({"selection": {"rows": rows, "columns": [], "cells": []}})
    return app._run(states)


def test_campaign_switch_replaces_opened_review_and_assessment_models(designs, monkeypatch):
    from streamlit.testing.v1 import AppTest
    db, group, molecule, models = designs
    second = db.create_group(OWNER, "Other campaign", origin={"tool": "rfd4_proteina", "type": "job",
        "campaign_id": "other-campaign", "run_id": "other-source", "task_group": "design"})
    db.add_to_group(OWNER, second, [molecule])
    db.attach_structure(second, molecule, {"name": "other.cif", "content": b"data_other\n"},
        context={"model_id": "other-model", "candidate_id": "other-candidate", "tool": "rfd4_proteina"})
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: db)
    monkeypatch.setattr(db, "structure_file", lambda *a: pytest.fail("Campaign navigation must not retrieve coordinates"))
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.complex_analysis_collection_ui import render_collection_workspace
render_collection_workspace({OWNER!r})
''').run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 2
    _select_rows(app, [0])
    assert not app.exception
    first_selection = dict(app.session_state["complex-analysis-evidence-selection"])
    app.button(key="complex-analysis-step-overview").click().run()
    assert not app.exception
    assert first_selection["label"] in [item.value for item in app.markdown]
    app.button(key="complex-analysis-step-select").click().run()
    app.session_state["complex-analysis-workflow-samples"] = [{"structure_id": "previous-model"}]
    app.session_state["complex-analysis-workflow-selection"] = {"structure_ids": ["previous-model"]}
    app.session_state["complex-analysis-assessment-open"] = True
    _select_rows(app, [1])
    assert not app.exception
    current = dict(app.session_state["complex-analysis-evidence-selection"])
    assert set(current["group_ids"]) != set(first_selection["group_ids"])
    assert not app.session_state["complex-analysis-assessment-open"]
    assert "complex-analysis-workflow-samples" not in app.session_state
    assert "complex-analysis-workflow-selection" not in app.session_state
    app.button(key="complex-analysis-step-overview").click().run()
    assert not app.exception
    assert current["label"] in [item.value for item in app.markdown]
    expected = 1 if current["group_ids"] == [second] else 2
    assert app.metric[0].value == str(expected)
    app.button(key="complex-analysis-step-select").click().run()
    assert app.session_state["complex-analysis-evidence-selection"] == current  # widget recreation retains source
    _select_rows(app, [])
    assert not app.exception
    assert "complex-analysis-evidence-selection" not in app.session_state
    assert app.button(key="complex-analysis-step-overview").disabled


@pytest.mark.parametrize("action", ["analyze", "compare"])
def test_assessment_action_and_step_use_only_the_review_selected_models(designs, monkeypatch, action):
    import hashlib
    import json
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import complex_analysis_collection_ui as collection
    db, group, _, models = designs
    token = hashlib.sha256(json.dumps([OWNER["id"], [group], None], sort_keys=True).encode()).hexdigest()[:12]
    opened = []
    def load(database, user, records, group_ids, **kwargs):
        opened.append((group_ids, kwargs))
        return []  # no scientific calculations or coordinate retrieval in this test
    monkeypatch.setattr(collection, "_load_collection_structures", load)
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: db)
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_design_review_ui import render_design_review
from mn_molecular_design_studio.molecular_database import MolecularDatabase
st.session_state.setdefault('design-review-selected-{token}', [{models[0]!r}])
render_design_review(MolecularDatabase(), {OWNER!r}, {{"group_ids": [{group!r}], "label": "Mixed engines"}},
                     requested_action=st.session_state.get('request', ''))
''').run()
    assert not app.exception
    app.button(key="review-assess-" + token).click().run()
    assert not app.exception
    assert opened == [([group], {"selected_structure_ids": [models[0]], "passing_cohort": False, "next_step": "analyze", "display_labels": {models[0]: "Design 1 · sample 0"}})]
    app.session_state['request'] = action
    app.run()
    assert not app.exception
    assert opened[-1] == ([group], {"selected_structure_ids": [models[0]], "passing_cohort": False, "next_step": action, "display_labels": {models[0]: "Design 1 · sample 0"}})
    # The second model has only native acceptance; it is never implicitly added.
    assert models[1] not in opened[-1][1]["selected_structure_ids"]


def test_assessment_screen_keeps_source_and_reuses_existing_submission_ui(designs, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import complex_analysis_collection_ui as collection
    db, group, _, _ = designs
    samples = design_review_page(db, OWNER, [group])["items"][:1]
    samples[0].update(kind="protein-protein", components=[
        {"type": "protein", "chain": "A", "sequence": "ACDEFGHIKLMNPQRSTVWYACDE"},
        {"type": "protein", "chain": "B", "sequence": "ACDEFGHIKLMNPQRSTVWYACDE"}])
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: db)
    monkeypatch.setattr(collection, "_load_workspace_samples", lambda user: samples)
    monkeypatch.setattr("mn_molecular_design_studio.complex_analysis_ui._matching_database_analysis_jobs", lambda *a, **kw: [])
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_evidence_ui import open_saved_work
from mn_molecular_design_studio.complex_analysis_collection_ui import render_collection_workspace
if 'complex-analysis-evidence-selection' not in st.session_state:
    open_saved_work(st.session_state, user_id={OWNER['id']!r}, group_ids=[{group!r}],
                    label='Selected mixed-engine campaign', work_type='Protein design', next_step='analyze')
    st.session_state['complex-analysis-assessment-open'] = True
render_collection_workspace({OWNER!r}, workflow_schema=lambda method: {{'parameters': []}},
                            submit=lambda *args: st.error('Unexpected calculation'))
''').run()
    assert not app.exception
    assert app.selectbox(key="complex-analysis-batch-profile").options == [
        "Core protein assessment", "Comprehensive protein assessment", "Custom"]
    assert any(item.value == 'Selected mixed-engine campaign' for item in app.markdown)
    assert any('Analysis coverage' in item.value for item in app.markdown)
    assert app.number_input(key="complex-analysis-batch-biotite-distance").value == 5.0
    assert app.multiselect[1].value == ["B"]
    assert app.multiselect[2].value == ["A"]
    assert not app.error
    assert not any(item.label == 'Acceptance cohort' for item in app.selectbox)
    assert not any(item.label == 'Design/model rows to assess and compare' for item in app.multiselect)
    app.selectbox(key="complex-analysis-batch-profile").select("Comprehensive protein assessment").run()
    assert not app.exception
    assert app.selectbox(key="complex-analysis-batch-tag-format")
    app.selectbox(key="complex-analysis-batch-profile").select("Custom").run()
    assert not app.exception
    app.multiselect(key="complex-analysis-batch-methods").select("analysis.biotite_interface").run()
    assert not app.exception
    assert not any(item.key == 'complex-analysis-batch-tag-format' for item in app.selectbox)
    assert not any(item.key == 'complex-analysis-batch-liability-format' for item in app.selectbox)
    assert not app.error
    app.button(key="complex-analysis-back-to-design-review").click().run()
    assert not app.exception
    assert any(item.value == '### Review designs' for item in app.markdown)
    assert not app.session_state['complex-analysis-assessment-open']


def test_changing_review_models_invalidates_only_the_old_loaded_selection():
    state = {"design-review-selected-cohort": ["model-a"], "complex-analysis-workflow-samples": ["old"],
             "complex-analysis-workflow-selection": {"structure_ids": ["model-a"]},
             "complex-analysis-evidence-selection": {"group_ids": ["campaign"]},
             "complex-analysis-assessment-open": True, "complex-analysis-cohort-selected-value": ["model-a"]}
    remember_review_selection(state, "cohort", ["model-b"], {"model-a", "model-b"})
    assert "complex-analysis-workflow-samples" not in state
    assert "complex-analysis-workflow-selection" not in state
    assert "complex-analysis-cohort-selected-value" not in state
    assert not state["complex-analysis-assessment-open"]
    assert state["complex-analysis-evidence-selection"] == {"group_ids": ["campaign"]}
    state["complex-analysis-workflow-samples"] = ["model-b"]
    remember_review_selection(state, "cohort", ["model-b"], {"model-a", "model-b"})
    assert state["complex-analysis-workflow-samples"] == ["model-b"]


def test_shortlist_step_uses_review_models_without_reselection_or_coordinates(designs, monkeypatch):
    import hashlib
    import json
    from streamlit.testing.v1 import AppTest
    db, group, _, models = designs
    token = hashlib.sha256(json.dumps([OWNER["id"], [group], None], sort_keys=True).encode()).hexdigest()[:12]
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", lambda: db)
    monkeypatch.setattr(db, "structure_file", lambda *a: pytest.fail("Shortlist needs no coordinates"))
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.complex_analysis_design_review_ui import render_design_review
st.session_state.setdefault('design-review-selected-{token}', [{models[0]!r}])
render_design_review(MolecularDatabase(), {OWNER!r}, {{"group_ids": [{group!r}], "label": "My campaign"}},
                     requested_action='shortlist')
''').run()
    assert not app.exception
    assert any(button.label == "Save selected-design review decisions" for button in app.button)
    assert len(app.get("download_button")) == 3
    assert not app.get("plotly_chart")


def test_assessment_navigation_preserves_existing_ligand_adapter(monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import complex_analysis_collection_ui as ui
    calls = []
    monkeypatch.setattr(ui, '_render_batch_assessment', lambda user, samples, **kw: calls.append(('protein', samples)))
    monkeypatch.setattr('mn_molecular_design_studio.complex_analysis_ui._render_database_analysis',
        lambda user, record, structure, path, components, *a, **kw: calls.append(('ligand', structure, path)))
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.complex_analysis_collection_ui import _render_selected_assessments
_render_selected_assessments({OWNER!r}, [{{'kind': 'protein-protein', 'structure_id': 'binder-model'}},
    {{'kind': 'protein-ligand', 'label': 'Pose 1', 'record': {{'id': 'compound'}},
      'structure': {{'id': 'exact-pose'}}, 'path': '/unused/preserved-pose', 'components': []}}])
''').run()
    assert not app.exception
    assert calls == [('protein', [{'kind': 'protein-protein', 'structure_id': 'binder-model'}]),
                     ('ligand', {'id': 'exact-pose'}, '/unused/preserved-pose')]
