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
    created = db.import_records(OWNER, molecule_imports.parse_smiles("CCO example"), name="Saved designs")
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
