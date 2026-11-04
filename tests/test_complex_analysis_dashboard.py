from __future__ import annotations

import json

from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.complex_analysis_dashboard import (
    _component_profile,
    _discover_structures,
    _evaluation_matches_candidate,
    _method_status,
    _result_sources,
    _run_action_label,
    _status_frame,
    _table_rows,
    _target_rows,
)


def test_component_profile_requires_recorded_roles_and_chemical_ligand_identity():
    item = {
        "candidate": {"candidate_id": "design-1"},
        "components": [
            {"type": "protein", "source_chain": "A"},
            {"type": "protein", "source_chain": "B"},
        ],
    }
    profile = _component_profile(item)
    assert profile["kind"] == "protein-protein"
    assert profile["roles_valid"] is False
    assert profile["binder_chains"] == []

    item["candidate"].update(binder_chains=["B"], target_chains=["A"])
    assert _component_profile(item)["roles_valid"] is True

    ligand = {
        "candidate": {"candidate_id": "pose-1"},
        "components": [
            {"type": "protein", "source_chain": "A"},
            {"type": "small_molecule", "source_chain": "L", "residue_key": "L:501", "ccd": "LIG"},
        ],
    }
    assert _component_profile(ligand)["ligand_selection"] is None
    ligand["candidate"]["raw_metadata"] = {"ligand_smiles": "CCO"}
    assert _component_profile(ligand)["ligand_selection"] == {
        "chain": "L", "residue_key": "L:501", "ccd": "LIG", "smiles": "CCO",
    }


def test_evaluation_matching_rejects_a_different_model_even_for_same_candidate():
    digest = "a" * 64
    candidate = {"candidate_id": "design-1", "model_index": 2}
    descriptor = {"candidate_id": "design-1_biotite_interface", "representation": digest}
    observation = {
        "context": {
            "source_structure_digest": digest,
            "candidate": {
                "candidate_id": "design-1_biotite_interface",
                "raw_metadata": {"complex_analysis_selection": {"model_id": "2"}},
            },
        },
    }
    assert _evaluation_matches_candidate(observation, descriptor, candidate, digest) == (True, True)
    observation["context"]["candidate"]["raw_metadata"]["complex_analysis_selection"]["model_id"] = "1"
    assert _evaluation_matches_candidate(observation, descriptor, candidate, digest) == (False, True)


def test_complex_analysis_defaults_to_run_first_dashboard():
    user = {"id": "user-1", "role": "admin"}
    run = {
        "id": "run:0123456789abcdef01234567",
        "name": "R3 BindCraft campaign",
        "status": "completed",
        "project_id": "project-1",
        "outputs": [{
            "id": "results",
            "name": "Results",
            "sources": [{
                "task_group": "bindcraft2",
                "run_id": "20261019-181522-21f9a96c",
                "tool": "bindcraft2",
                "status": "completed",
                "project_id": "project-1",
                "candidate_summary": {"candidate_count": 65, "structure_count": 65},
            }],
        }],
        "jobs": [],
    }
    script = f'''import streamlit as st
from mn_molecular_design_studio import logical_runs
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
logical_runs.catalog = lambda user, token="": [{run!r}]
render_complex_analysis_workspace({user!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    assert any(item.value == "Complex Analysis" for item in app.title)
    assert any("Run-first workflow" in item.value for item in app.info)
    assert any(item.label == "Select a run" for item in app.button)
    assert any("R3 BindCraft campaign" in str(item.value) for item in app.dataframe)


def test_inventory_distinguishes_not_applicable_unavailable_and_stale_methods():
    item = {
        "candidate_id": "design-1", "digest": "a" * 64, "model_id": "1",
        "relative_path": "poses/design-1.pdb", "candidate": {
            "candidate_id": "design-1", "binder_chains": [], "target_chains": [],
        },
        "components": [{"type": "protein", "source_chain": "A"},
                       {"type": "protein", "source_chain": "B"}],
    }
    frame = _status_frame([item], {}, {})
    statuses = dict(zip(frame["Method"], frame["Status"]))
    assert statuses["Biotite interface geometry"] == "Unavailable"
    assert statuses["PEP-Patch surface chemistry"] == "Unavailable"
    assert statuses["PLIP ligand interactions"] == "Not applicable"
    assert statuses["PoseBusters pose checks"] == "Not applicable"

    profile = {**_component_profile(item), "binder_chains": ["B"], "target_chains": ["A"]}
    database_index = {("design-1", "biotite_interface"): [{
        "match": "Exact", "descriptor": {"status": "completed"},
        "observation": {"context": {"candidate": {
            "candidate_id": "design-1", "binder_chains": ["A"], "target_chains": ["B"],
            "raw_metadata": {"complex_analysis_selection": {"model_id": "1"}},
        }}},
    }]}
    assert _method_status(database_index, {}, item, profile, "biotite_interface")[0] == "Stale or nonmatching"


def test_result_run_rows_show_only_exact_linked_analysis_activity(monkeypatch):
    from mn_molecular_design_studio import cluster

    monkeypatch.setattr(cluster.config, "identity", lambda: "node-a")
    source = {"task_group": "bindcraft2", "run_id": "run-1", "node_id": "node-a",
              "tool": "bindcraft2", "updated_at": "2026-11-01",
              "candidate_summary": {"candidate_count": 2, "structure_count": 2}}
    run = {"id": "logical-1", "name": "Run one", "jobs": [
        {"task_group": "analysis", "run_id": "analysis-1", "node_id": "node-a",
         "upstream_task_group": "bindcraft2", "upstream_run_id": "run-1",
         "updated_at": "2026-11-03"},
        {"task_group": "analysis", "run_id": "analysis-2", "node_id": "node-b",
         "upstream_task_group": "bindcraft2", "upstream_run_id": "run-1",
         "updated_at": "2026-11-04"},
    ], "outputs": [{"id": "results", "sources": [source]}]}
    row = _result_sources([run])[0]
    assert row["analysis_job_count"] == 1
    assert row["analysis_job_statuses"] == ["unknown"]
    assert row["last_analysis_activity"] == "2026-11-03"
    assert _table_rows([row])[0]["Last analysis activity"] == "2026-11-03"


def test_result_sources_aggregates_child_jobs_and_excludes_empty_runs(monkeypatch):
    from mn_molecular_design_studio import cluster

    monkeypatch.setattr(cluster.config, "identity", lambda: "local")
    sources = [
        {"task_group": "design", "run_id": "batch-1", "node_id": "local", "tool": "rfd4",
         "candidate_summary": {"candidate_count": 50, "structure_count": 50,
                                "target_systems": [{"identity": "sequence:abc", "label": "R3 · 31 aa"}]}},
        {"task_group": "design", "run_id": "batch-2", "node_id": "local", "tool": "rfd4",
         "candidate_summary": {"candidate_count": 50, "structure_count": 1,
                                "target_systems": [{"identity": "sequence:abc", "label": "R3 · 31 aa"}]}},
        {"task_group": "design", "run_id": "batch-empty", "node_id": "local", "tool": "rfd4",
         "candidate_summary": {"candidate_count": 0, "structure_count": 0}},
    ]
    rows = _result_sources([{"id": "logical-run", "name": "R3 · accepted 200",
                             "job_count": 3, "jobs": [],
                             "outputs": [{"id": "results", "sources": sources}]}])
    assert len(rows) == 1
    assert rows[0]["logical_run_job_count"] == 3
    assert rows[0]["source_job_count"] == 3
    assert rows[0]["structure_count"] == 51
    assert rows[0]["candidate_count"] == 100
    assert rows[0]["target_summary"]["identity"] == "sequence:abc"
    assert _table_rows(rows)[0]["Structures"] == 51


def test_run_without_verified_structures_or_sequence_only_target_is_not_an_analysis_row(monkeypatch):
    from mn_molecular_design_studio import cluster

    monkeypatch.setattr(cluster.config, "identity", lambda: "local")
    run = {"id": "empty-run", "name": "Sequence only", "jobs": [], "outputs": [{"id": "results", "sources": [
        {"task_group": "design", "run_id": "sequence-only", "node_id": "local", "tool": "rfd4",
         "candidate_summary": {"candidate_count": 3, "structure_count": 0,
                                "target_systems": [{"identity": "sequence:target", "label": "Target · 42 aa"}]}}
    ]}]}
    assert _result_sources([run]) == []


def test_target_view_groups_runs_by_exact_sequence_identity():
    target = {"identity": "sequence:abc", "label": "Target · 70 aa", "identity_kind": "sequence",
              "reference_type": "sequence_only"}
    structured_target = {**target, "reference_type": "structure_backed"}
    rows = [
        {"logical_run_id": "run-a", "structure_count": 50, "engine_names": ["RFD4"],
         "run_sources": [{"source_summary": {"target_systems": [target]}}]},
        {"logical_run_id": "run-b", "structure_count": 10, "engine_names": ["OpenFold3"],
         "run_sources": [{"source_summary": {"target_systems": [structured_target]}}]},
        {"logical_run_id": "run-c", "structure_count": 1, "engine_names": ["RFD4"], "run_sources": []},
    ]

    result, unrecorded = _target_rows(rows)

    assert len(result) == 1
    assert result[0]["Identity"] == "sequence:abc"
    assert result[0]["Runs"] == 2
    assert result[0]["Result structures"] == 60
    assert result[0]["Engines"] == "OpenFold3, RFD4"
    assert result[0]["Input target"] == "Sequence-only and structure-backed runs"
    assert unrecorded == 1


def test_run_action_label_distinguishes_new_resumable_and_reviewable_runs():
    assert _run_action_label({}) == "Analyze"
    assert _run_action_label({"analysis_job_count": 1, "analysis_job_statuses": ["completed"]}) == "Open dashboard"
    assert _run_action_label({"analysis_job_count": 1, "analysis_job_statuses": ["running"]}) == "Continue analysis"


def test_run_inventory_ignores_unassociated_input_and_reference_structures(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui

    candidate_path = tmp_path / "candidate.pdb"
    input_path = tmp_path / "target.pdb"
    candidate_path.write_text("ATOM\n")
    input_path.write_text("ATOM\n")
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda user, job: [
        {"path": input_path, "relative_path": "artifacts/pipeline_inputs/target.pdb", "candidate": {}},
        {"path": candidate_path, "relative_path": "artifacts/normalized_candidates/candidate.pdb",
         "candidate_id": "candidate-1", "candidate": {"candidate_id": "candidate-1"}},
    ])
    monkeypatch.setattr(complex_analysis_ui, "component_inventory", lambda *args, **kwargs: [])
    found = _discover_structures({"id": "user-1"}, {"task_group": "design", "run_id": "run-1"})
    assert [item["candidate_id"] for item in found] == ["candidate-1"]
