from __future__ import annotations

import csv
import json

import pytest

from mn_molecular_design_studio.complex_analysis_shortlist_export import (
    NOT_SHORTLISTED,
    SHORTLISTED,
    build_shortlist_manifest,
    build_shortlist_record,
    manifest_csv,
    manifest_json,
)


def _sample(*, source_kind="Results", final_status="passed", native_pass=True, model="model-2"):
    candidate = {
        "candidate_id": "studio-candidate-7",
        "native_candidate_id": "engine-candidate-7",
        "model_id": model,
        "source_tool": "rfd4_proteina",
        "binder_sequence": "G" * 12,
        "binder_chains": ["B"],
        "target_chains": ["A"],
        "metrics": {"native_pass_filters": native_pass, "sampled_confidence": 0.83},
        "raw_metadata": {
            "generator_engine": "rfd4_proteina",
            "prediction_engine": "boltz2",
            "target_identity": "target-construct-A",
            "target_molecule_id": "target-molecule-1",
            "source_candidate": {"candidate_id": "engine-candidate-7"},
            **({"validation_acceptance": {"status": final_status, "policy": {"version": "accept-v2"}}}
               if final_status else {}),
        },
    }
    return {
        "source_kind": source_kind,
        "candidate": candidate,
        "candidate_id": candidate["candidate_id"],
        "molecule_id": "mol-7" if source_kind == "Database" else "",
        "structure_id": "structure-7",
        "record": {"id": "mol-7", "sources": [{
            "source_key": "pc004033/design/run-7:native-id",
            "candidate_id": "engine-candidate-7",
            "task_group": "design", "run_id": "run-7", "tool": "rfd4_proteina",
            "workstation": "PC004033", "parents": [], "run_dir": "/private/workdir/run-7",
        }]},
        "model_id": model,
        "digest": "a" * 64,
        "source_reference": "cluster:workstation/run-7" if source_kind == "Results" else "database:mol-7",
        "source_job": "prediction/run-7" if source_kind == "Results" else "design/run-source-7",
        "source_workstation": "PC004033",
        "project_id": "project-1",
        "relative_path": "artifacts/prediction/model-2.cif",
        "context": {"task_group": "design", "run_id": "run-source-7", "workstation": "PC004033",
                    "job_reference": "cluster:source-7"},
        "tracks": [
            {"chain_id": "A", "sequence": "A" * 15},
            {"chain_id": "B", "sequence": "G" * 12},
        ],
        "evaluation_entries": [{
            "evaluation_id": "eval-1", "evaluator": "binder_interface_chemistry_v1",
            "method_version": "studio.ppi-interface-chemistry.v1", "status": "completed",
            "run_id": "analysis-1", "metrics": {"interface_total_delta_sasa_angstrom2": 1000.0},
            "provenance": {"source_model_id": model, "source_structure_digest": "a" * 64,
                           "target_chains": ["A"], "binder_chains": ["B"], "candidate": {
                               "raw_metadata": {"binder_interface_chemistry": {
                                   "contacts": [{"target": {"author_number": 8}}],
                                   "hydrogen_bonds": [],
                                   "residue_annotations": [{"chain_id": "A", "author_number": 8}],
                                   "methodological_assumptions": ["geometry-based"],
                               }}},
                           },
        }],
        "assessment_coverage": {"binder_interface_chemistry": "compatible saved Evaluation"},
        "review_annotation": {"status": SHORTLISTED, "note": "Shortlist for pilot validation only."},
    }


def test_export_preserves_candidate_structure_evaluation_and_explicit_acceptance():
    sample = _sample(source_kind="Results")
    row = build_shortlist_record(sample)

    assert row["candidate_id"] == "studio-candidate-7"
    assert row["native_candidate_id"] == "engine-candidate-7"
    assert row["structure_sha256"] == "a" * 64
    assert row["model_id"] == "model-2"
    assert row["source"]["job"] == "prediction/run-7"
    assert row["original_artifact_reference"] == "cluster:source-7"
    assert row["source_lineage"][0]["candidate_id"] == "engine-candidate-7"
    assert "run_dir" not in row["source_lineage"][0]
    assert row["target"]["sequences_by_chain"] == {"A": "A" * 15}
    assert row["binder"]["sequences_by_chain"] == {"B": "G" * 12}
    assert row["acceptance"]["native"] == "passed"
    assert row["acceptance"]["final_computational"] == "passed"
    assert row["evaluations"][0]["method_version"] == "studio.ppi-interface-chemistry.v1"
    assert row["evaluations"][0]["measurements"]["interface_total_delta_sasa_angstrom2"] == 1000.0
    assert row["evaluations"][0]["evidence_summary"] == {
        "contact_pairs": 1, "hydrogen_bonds": 0, "residue_annotations": 1,
    }
    assert row["evaluations"][0]["residue_level_evidence"]["residue_annotations"][0]["author_number"] == 8
    assert row["human_review"]["status"] == SHORTLISTED


def test_native_pass_is_not_promoted_to_final_acceptance_and_unknown_native_id_stays_unknown():
    sample = _sample(source_kind="Database", final_status=None, native_pass=True)
    sample["candidate"].pop("native_candidate_id")
    sample["candidate"]["raw_metadata"].pop("source_candidate")
    row = build_shortlist_record(sample)

    assert row["acceptance"]["native"] == "passed"
    assert row["acceptance"]["final_computational"] == "not recorded"
    assert row["native_candidate_id"] == "Not recorded"
    assert row["source"]["kind"] == "Database"
    assert row["original_artifact_reference"] == "cluster:source-7"


def test_export_does_not_guess_source_kind_from_molecular_complex_type_or_export_absolute_paths():
    sample = _sample(source_kind="Results")
    sample.pop("source_kind")
    sample["kind"] = "protein-protein"
    sample["context"]["job_reference"] = "/home/another-workstation/private/run/artifact.cif"
    row = build_shortlist_record(sample)

    assert row["source"]["kind"] == "Unknown"
    assert row["original_artifact_reference"] == ""


def test_manifest_json_and_csv_are_machine_readable_and_shortlisted_export_is_explicit():
    selected = _sample(source_kind="Results")
    other = _sample(source_kind="Database", final_status=None)
    other["candidate"]["candidate_id"] = "studio-candidate-8"
    other["candidate"]["native_candidate_id"] = "engine-candidate-8"
    other["candidate_id"] = "studio-candidate-8"
    other["molecule_id"] = "mol-8"
    other["structure_id"] = "structure-8"
    other["digest"] = "b" * 64
    other["review_annotation"] = {"status": NOT_SHORTLISTED, "note": "Excluded for this pilot."}
    manifest = build_shortlist_manifest([selected, other], exported_at="2026-11-08T12:00:00+00:00")
    restored = json.loads(manifest_json(manifest))
    rows = list(csv.DictReader(manifest_csv(manifest).decode().splitlines()))
    shortlist_rows = list(csv.DictReader(manifest_csv(manifest, shortlisted_only=True).decode().splitlines()))

    assert restored["record_count"] == 2
    assert restored["shortlisted_count"] == 1
    assert len(rows) == 2 and len(shortlist_rows) == 1
    assert shortlist_rows[0]["candidate_id"] == "studio-candidate-7"
    assert json.loads(rows[1]["recorded_acceptance_json"])[0]["stage_id"] == "native"
    assert rows[1]["final_computational_acceptance"] == "not recorded"
    assert "affinity" in restored["interpretation"]
    assert "universal quality score" in restored["interpretation"]


def test_manifest_rejects_duplicate_source_model_structure_identity():
    sample = _sample()
    with pytest.raises(ValueError, match="duplicate source/model/structure"):
        build_shortlist_manifest([sample, sample])


def test_shortlist_panels_can_render_twice_for_same_cohort_in_distinct_views():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        """
from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export
sample = {
    'source_kind': 'Results', 'candidate_id': 'c-1', 'candidate': {'candidate_id': 'c-1', 'raw_metadata': {}},
    'structure_id': 'structure-1', 'model_id': 'model-1', 'digest': 'a' * 64,
    'source_reference': 'runs:///source', 'source_job': 'design/run-1', 'source_workstation': 'PC004033',
}
render_shortlist_review_export({'id': 'test-user'}, [sample], key_prefix='explorer-shortlist')
render_shortlist_review_export({'id': 'test-user'}, [sample], key_prefix='comparison-shortlist')
""",
        default_timeout=10,
    ).run()
    assert not app.exception


def test_namespaced_shortlist_save_uses_the_rendered_panel_values(monkeypatch):
    from streamlit.testing.v1 import AppTest
    import hashlib
    from mn_molecular_design_studio import molecular_database

    # AppTest shares imported modules with pytest. Register the original class
    # so the script's fake cannot escape into subsequent Database/API tests.
    monkeypatch.setattr(molecular_database, "MolecularDatabase", molecular_database.MolecularDatabase)

    token = hashlib.sha256(b"mol-1:structure-1").hexdigest()[:12]
    app = AppTest.from_string(
        f"""
import streamlit as st
import mn_molecular_design_studio.molecular_database as db_module
from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export

class FakeDatabase:
    def group(self, group_id, user):
        return {{"project_id": "staging"}}
    def get_complex_review(self, user, **kwargs):
        return {{"status": "Unreviewed", "note": ""}}
    def save_complex_review(self, user, **kwargs):
        st.session_state.setdefault("saved_reviews", []).append(kwargs)
        return {{"status": kwargs["status"], "note": kwargs["note"]}}

db_module.MolecularDatabase = FakeDatabase
sample = {{
    "source_kind": "Database", "candidate_id": "c-1", "molecule_id": "mol-1",
    "structure_id": "structure-1", "digest": "a" * 64,
    "record": {{"id": "mol-1"}}, "structure": {{"id": "structure-1", "group_id": "group-1"}},
    "candidate": {{"candidate_id": "c-1", "raw_metadata": {{}}}},
}}
st.session_state["explorer-shortlist-decision-{token}"] = "Shortlisted for experimental testing"
st.session_state["explorer-shortlist-rationale-{token}"] = "Advance for pilot review."
render_shortlist_review_export({{"id": "test-user"}}, [sample], key_prefix="explorer-shortlist")
render_shortlist_review_export({{"id": "test-user"}}, [sample], key_prefix="comparison-shortlist")
""",
        default_timeout=10,
    ).run()
    assert not app.exception
    app.button(key="explorer-shortlist-save-review").click().run()
    assert not app.exception
    saved = app.session_state["saved_reviews"]
    assert len(saved) == 1
    assert saved[0]["status"] == "Shortlisted for experimental testing"
    assert saved[0]["note"] == "Advance for pilot review."
