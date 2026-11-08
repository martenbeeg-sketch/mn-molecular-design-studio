from __future__ import annotations

from pathlib import Path

from mn_molecular_design_studio.complex_analysis_collection import acceptance_decisions
from mn_molecular_design_studio.complex_analysis_collection_ui import (
    _load_selected_structures,
    _merge_database_record_options,
)


def test_database_record_picker_keeps_explicitly_selected_rows_across_pages_and_searches():
    rfd4 = {"id": "mol-rfd4", "name": "rfd4", "representations": [{"id": "structure-rfd4"}]}
    bindcraft = {"id": "mol-bindcraft", "name": "bindcraft2", "representations": [{"id": "structure-bindcraft"}]}
    unrelated = {"id": "mol-unrelated", "name": "unrelated", "representations": []}

    # The visible page is BindCraft-only, but an explicitly selected RFD4 row
    # remains available from the session cache. Unselected prior-page rows do
    # not leak into the picker.
    options = _merge_database_record_options(
        [bindcraft], {"mol-rfd4": rfd4, "mol-unrelated": unrelated}, ["mol-rfd4"],
    )

    assert set(options) == {"mol-rfd4", "mol-bindcraft"}
    assert options["mol-rfd4"]["representations"][0]["id"] == "structure-rfd4"


def test_selected_database_structure_hydrates_its_saved_evaluations_and_provenance(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_inspector, complex_analysis_ui

    user = {"id": "staging-reviewer", "role": "admin"}
    structure_id = "peer-structure-1"
    candidate_id = "bindcraft2_00001"
    molecule_id = "mol-bindcraft-1"
    structure = {"id": structure_id, "context": {"candidate_id": candidate_id, "task_group": "analysis"}}
    source_candidate = {"candidate_id": candidate_id, "source_tool": "bindcraft2", "metrics": {
        "native_pass_filters": True,
    }}
    observation = {
        "id": "obs-bindcraft-1",
        "namespace": "bindcraft2",
        "run_id": "source-run-1",
        "metrics": {"native_pass_filters": True},
        "context": {"representation_id": structure_id, "candidate": source_candidate},
    }
    descriptor = {
        "id": "evaluation-bindcraft-1",
        "evaluator": "bindcraft2",
        "representation_id": structure_id,
        "candidate_id": candidate_id,
        "run_id": "source-run-1",
    }
    light_record = {"id": molecule_id, "representations": [structure], "structures": [structure]}
    full_record = {
        **light_record,
        "components": [],
        "sources": [{"candidate_id": candidate_id, "source_tool": "bindcraft2"}],
        "observations": [observation],
    }

    class Database:
        def __init__(self):
            self.hydration_calls = []

        def list_molecules(self, selected_user, **kwargs):
            self.hydration_calls.append((selected_user, kwargs))
            return {"items": [full_record]}

        def structure_file(self, selected_user, selected_structure_id):
            assert selected_user is user
            assert selected_structure_id == structure_id
            path = tmp_path / "structure.cif"
            path.write_text("fixture coordinates")
            return path

    database = Database()
    monkeypatch.setattr(complex_analysis_ui, "component_inventory", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(complex_analysis_ui, "_evaluation_options",
                        lambda _record: ([], {descriptor["id"]: (observation, descriptor)}))
    monkeypatch.setattr(complex_analysis_ui, "_database_inspector_data",
                        lambda *_args, **_kwargs: {"evaluations": [], "evidence": [], "tracks": []})
    monkeypatch.setattr(complex_analysis_inspector, "protein_sequence_tracks", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(complex_analysis_inspector, "structure_digest", lambda _path: "sha256:fixture")

    samples = _load_selected_structures(database, user, [light_record], [structure_id])

    assert len(samples) == 1
    assert database.hydration_calls == [(user, {
        "collection": "all", "molecule_ids": [molecule_id], "limit": 1,
    })]
    assert samples[0]["candidate"]["candidate_id"] == candidate_id
    assert samples[0]["candidate"]["source_tool"] == "bindcraft2"
    assert samples[0]["associated"] == [(observation, descriptor)]
    native = [item for item in acceptance_decisions(samples[0]["candidate"]) if item["stage_id"] == "native"]
    assert len(native) == 1
    assert native[0]["cohort_status"] == "passed"
    assert native[0]["source"] == "native_pass_filters"
