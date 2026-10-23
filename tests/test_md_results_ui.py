from __future__ import annotations

import json
from pathlib import Path

import pytest

from mn_molecular_design_studio.md_results_ui import (
    _endpoint_rows,
    _finite_pairs,
    _read_json_artifact,
    _safe_artifact_path,
    load_md_study_manifest,
    render_md_results,
)


class FakeStreamlit:
    def __init__(self):
        self.calls = []

    def __getattr__(self, name):
        def record(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            if name == "selectbox":
                options = args[1]
                return options[0]
            return None

        return record


def _json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return path


def _study_manifest():
    return {
        "schema_version": 1,
        "studies": [{
            "study_id": "study-1",
            "compound_id": "cmp-123",
            "compound_name": "compound A",
            "target_pdb_id": "1ABC",
            "target_receptor": "PDB chain A",
            "md_engine": "GROMACS",
            "source_pose_id": "pose-4",
            "source_tool": "GNINA",
            "selection_status": "manual_nonbest",
            "selected_complex_run_id": "complex-7",
            "source_pose_rank": 4,
            "production_runs": [
                {"run_id": "prod-1", "status": "completed", "repeat_index": 1, "target_duration_ns": 20, "continuation_of_run_id": "", "source_completed_at": "2025-02-01T00:00:00Z", "result_path": "artifacts/legacy/prod-1/result.json"},
                {"run_id": "prod-2", "status": "failed", "repeat_index": 2, "target_duration_ns": 20, "continuation_of_run_id": "prod-1"},
            ],
            "mmgbsa_runs": [{"run_id": "endpoint-1", "status": "completed", "source_production_run_id": "prod-1", "result_path": "artifacts/legacy/source/endpoint/result.json", "superseded_by_run_id": "endpoint-2", "recorded_trajectory_sha256": "a" * 64, "copied_trajectory_sha256": "b" * 64, "trajectory_snapshot_matches": False, "source_md_job": {"task_group": "ligand-md", "run_id": "studio-md", "workstation": "pc004033"}, "legacy_source_md_job": {"task_group": "bound-ligand-md", "run_id": "prod-1", "workstation": "pc004033"}, "source_trajectory_ref": "cluster://node/ligand-md/studio-md/artifacts/legacy/prod-1/traj.dcd", "parameters": {"backend": "ambertools_mmpbsa", "stride": 10, "start_pct": 20, "end_pct": 80}}],
            "analysis_runs": [{"run_id": "analysis-1", "status": "completed", "result_path": "artifacts/md/analysis.json"}],
        }],
    }


def test_manifest_renderer_shows_study_lineage_replica_analysis_and_separate_endpoint_terms(tmp_path):
    run = tmp_path / "run"
    _json(run / "artifacts/md/study.json", _study_manifest())
    _json(run / "artifacts/md/analysis.json", {
        "replicas": [{"replica": 1, "status": "completed", "delta_g_bind_kcal_mol": -4.5}],
        "replica_series": [{"replica": 1, "time_ns": [0, 1, 2], "backbone_rmsd_angstrom": [0, 1, 2], "ligand_rmsd_angstrom": [0.5, 1.5, 2.5]}],
        "contact_consensus": [{"residue": "A:42", "mean_contact_occupancy": 0.75}],
        "contact_matrix": {"residues": ["A:42"], "replicas": [1], "contact_occupancy": [[0.75]]},
        "aggregate": {"completed": 1, "duration_ns": 20},
    })
    _json(run / "artifacts/legacy/prod-1/result.json", {
        "success": True,
        "md_result": {"status": "completed", "analytics": {
            "rmsd": {"time_ps": [0, 1000], "backbone_rmsd_angstrom": [0.2, 0.4], "ligand_rmsd_angstrom": [0.5, 0.7]},
            "structural_dynamics": {
                "time_ps": [0, 1000],
                "rmsf": {"residues": ["A:1", "A:2"], "ca_rmsf_angstrom": [0.8, 1.2]},
                "ligand_rmsf": {"atoms": ["C1", "O1"], "rmsf_angstrom": [0.6, 1.4]},
                "radius_of_gyration": {"protein_angstrom": [10, 11], "ligand_angstrom": [2, 2.1], "complex_angstrom": [12, 12.5]},
                "secondary_structure": {"helix_fraction": [0.4, 0.5], "sheet_fraction": [0.2, 0.1], "coil_fraction": [0.4, 0.4]},
                "pocket": {"retained_fraction": 0.75, "site_displacement_cutoff_angstrom": 5, "contact_cutoff_angstrom": 4.5,
                           "ligand_centroid_displacement_angstrom": [0.1, 1.2], "reference_site_retained": [True, False]},
                "contacts": {"residues": [{"residue": "A:42", "contact_occupancy": 0.8, "hydrogen_bond_occupancy": 0.4,
                                             "hydrophobic_occupancy": 0.5, "mean_minimum_distance_angstrom": 3.2}],
                             "distance_series": {"A:42": [3.0, 3.4]},
                             "interaction_hotspot_score": {"description": "stored weighted interaction occupancy", "weights": {"hydrogen_bond": 3}}},
                "water_bridges": {"method": "stored geometry method", "residues": [{"residue": "A:42", "occupancy": 0.2}]},
                "salt_bridges": {"applicable": True, "distance_cutoff_angstrom": 4.0, "residues": [{"residue": "A:99", "occupancy": 0.1}]},
                "interaction_network": {"nodes": [{"id": "A:42", "contact_occupancy": 0.8}], "edges": [{"source": "C1", "target": "A:42", "weight": 0.4, "interaction": "hydrogen_bond"}]},
                "interface_rin": {"applicable": True, "interfaces": [{"chain_a": "A", "chain_b": "B", "contact_fraction": 0.4}],
                                  "residue_pairs": [{"residue_a": "A:1", "residue_b": "B:2", "contact_occupancy": 0.5}], "nodes": []},
            },
        }},
    })
    _json(run / "artifacts/legacy/source/endpoint/result.json", {
        "success": True,
        "mmgbsa": {
            "status": "success",
            "method": "MMPBSA.py",
            "units": "kcal/mol",
            "metadata": {"n_frames_analyzed": 60, "n_frames_total": 300},
            "gb": {"delta": {"delta_g_bind_total_kcal_mol": -5.2, "delta_gbsa_kcal_mol": -2.0}},
            "pb": {"delta": {"delta_g_bind_total_kcal_mol": -3.1, "delta_gbsa_kcal_mol": -1.4}},
        },
    })
    st = FakeStreamlit()

    render_md_results(st, run, {"success": True})

    strings = json.dumps(st.calls, default=str)
    assert "cmp-123" in strings
    assert "GNINA" in strings and "pose-4" in strings
    assert "manual_nonbest" in strings
    assert "prod-1" in strings and "prod-2" in strings
    assert "continues" not in strings  # lineage is shown using the source run ID field
    assert "endpoint-2" in strings
    endpoint_details = next(call[1][0] for call in st.calls if call[0] == "json" and "Trajectory snapshot matches" in str(call[1][0]))
    assert endpoint_details["Trajectory snapshot matches"] is False
    assert endpoint_details["Source production status"] == "completed"
    assert endpoint_details["Source MD job"]["run_id"] == "studio-md"
    assert endpoint_details["Legacy source MD job"]["run_id"] == "prod-1"
    assert "different trajectory snapshot" in strings
    lineage = next(call[1][0] for call in st.calls if call[0] == "dataframe" and "This segment (ns)" in str(call[1][0]))
    assert lineage[0]["Requested total (ns)"] == 20
    assert lineage[0]["This segment (ns)"] is None
    assert -5.2 in [value for call in st.calls if call[0] == "dataframe" for row in call[1][0] for value in row.values() if isinstance(value, float)]
    assert -3.1 in [value for call in st.calls if call[0] == "dataframe" for row in call[1][0] for value in row.values() if isinstance(value, float)]
    assert any(call[0] == "line_chart" for call in st.calls)
    assert sum(call[0] == "dataframe" for call in st.calls) >= 4
    rendered = json.dumps(st.calls, default=str, ensure_ascii=False)
    assert "Protein Cα RMSF" in rendered and "Ligand atom RMSF" in rendered
    assert "Reference-site retention" in rendered and "Recorded interaction network" in rendered
    assert "Water bridges" in rendered and "Salt bridges" in rendered
    assert "Interface residue pairs" in rendered or "Residue pairs" in rendered
    production_choices = [call for call in st.calls if call[0] == "selectbox" and call[1][0] == "Production result"]
    assert production_choices and production_choices[0][1][1][0].startswith("prod-1")


def test_safe_artifact_resolution_rejects_parent_traversal_and_external_symlinks(tmp_path):
    root = tmp_path / "run"
    root.mkdir()
    external = tmp_path / "outside.json"
    external.write_text("{}")
    (root / "link.json").symlink_to(external)

    with pytest.raises(ValueError, match="safe relative"):
        _safe_artifact_path(root.resolve(), "../outside.json")
    with pytest.raises(ValueError, match="outside"):
        _safe_artifact_path(root.resolve(), "link.json")
    with pytest.raises(ValueError, match="safe relative"):
        _safe_artifact_path(root.resolve(), "/tmp/outside.json")


def test_manifest_rejects_external_symlink_and_bad_schema(tmp_path):
    root = tmp_path / "run"
    root.mkdir()
    outside = _json(tmp_path / "outside.json", {"schema_version": 1, "studies": []})
    (root / "artifacts/md").mkdir(parents=True)
    (root / "artifacts/md/study.json").symlink_to(outside)
    with pytest.raises(ValueError, match="outside"):
        load_md_study_manifest(root)

    (root / "artifacts/md/study.json").unlink()
    _json(root / "artifacts/md/study.json", {"schema_version": 99, "studies": []})
    with pytest.raises(ValueError, match="Unsupported"):
        load_md_study_manifest(root)


def test_nonfinite_and_mismatched_plot_series_are_filtered():
    assert _finite_pairs([0, 1, 2], [1, float("nan"), 3, 4]) == [(0.0, 1.0), (2.0, 3.0)]


def test_mmgbsa_pb_and_gb_remain_separate_and_failed_status_is_not_inferred_from_file(tmp_path):
    endpoint = {
        "success": False,
        "mmgbsa": {
            "status": "failed",
            "units": "kcal/mol",
            "gb": {"delta": {"delta_g_bind_total_kcal_mol": -4.0}},
            "pb": {"delta": {"delta_g_bind_total_kcal_mol": -2.0}},
        },
    }
    row = _endpoint_rows(endpoint, {"run_id": "ep-1", "status": "failed"})[0]
    assert row["GB ΔG bind (kcal/mol)"] == -4
    assert row["PB ΔG bind (kcal/mol)"] == -2
    assert row["Result status"] == "failed"


def test_legacy_result_fallback_shows_direct_analysis_and_endpoint(tmp_path):
    run = tmp_path / "legacy"
    _json(run / "artifacts/md/result.json", {
        "success": False,
        "md_result": {
            "analytics": {
                "rmsd": {"time_ps": [0, 1000], "backbone_rmsd_angstrom": [0, float("inf")], "ligand_rmsd_angstrom": [1, 2]},
            },
        },
        "mmgbsa": {"status": "not_requested"},
    })
    st = FakeStreamlit()

    render_md_results(st, run, {"success": True})

    assert any(call[0] == "line_chart" for call in st.calls)
    assert any(call[0] == "dataframe" for call in st.calls)
    assert "not_requested" in json.dumps(st.calls, default=str)


def test_preloaded_manifest_uses_validated_loader_and_job_scoped_widget_keys():
    requested = []
    manifest = _study_manifest()

    def load_json(relative, label):
        requested.append(relative)
        if relative == "../outside.json":
            raise AssertionError("unsafe path reached callback")
        return None

    manifest["studies"][0]["mmgbsa_runs"] = [{"run_id": "bad", "status": "completed", "result_path": "../outside.json"}]
    manifest["studies"][0]["analysis_runs"] = [{"run_id": "stage-only", "status": "completed"}]
    st = FakeStreamlit()

    render_md_results(st, manifest=manifest, result={}, json_loader=load_json, key_prefix="remote-job", remote=True)

    assert requested == ["artifacts/legacy/prod-1/result.json"]
    keys = [call[2].get("key") for call in st.calls if call[0] == "selectbox"]
    assert keys == ["remote-job-study", "remote-job-production-study-1"]
    assert any(call[0] == "info" for call in st.calls)


def test_analysis_choices_put_latest_successful_result_first(tmp_path):
    run = tmp_path / "run"
    study = _study_manifest()["studies"][0]
    study["production_runs"] = []
    study["mmgbsa_runs"] = []
    study["analysis_runs"] = [
        {"run_id": "older-success", "status": "completed", "source_completed_at": "2025-01-01T00:00:00Z", "result_path": "artifacts/md/old.json"},
        {"run_id": "newer-success", "status": "success", "source_completed_at": "2025-02-01T00:00:00Z", "result_path": "artifacts/md/new.json"},
        {"run_id": "latest-failed", "status": "failed", "source_completed_at": "2025-03-01T00:00:00Z", "result_path": "artifacts/md/failed.json"},
    ]
    _json(run / "artifacts/md/study.json", {"schema_version": 1, "studies": [study]})
    _json(run / "artifacts/md/new.json", {"aggregate": {"selected": "newer-success"}})
    st = FakeStreamlit()

    render_md_results(st, run, {})

    analysis_choices = next(call[1][1] for call in st.calls if call[0] == "selectbox" and call[1][0] == "Analysis result")
    assert analysis_choices[0].startswith("newer-success")


def test_remote_production_view_loads_only_selected_owned_result_and_never_source_history_paths():
    study = _study_manifest()["studies"][0]
    study["production_runs"] = [
        {"run_id": "old", "status": "completed", "source_completed_at": "2025-01-01T00:00:00Z", "result_path": "artifacts/legacy/old/result.json"},
        {"run_id": "latest", "status": "completed", "source_completed_at": "2025-02-01T00:00:00Z", "result_path": "artifacts/legacy/latest/result.json"},
    ]
    study["analysis_runs"] = []
    study["mmgbsa_runs"] = []
    study["source_history"] = {
        "source_job": {"task_group": "old-group", "run_id": "old-job", "run_relative_path": "../outside/old-job"},
        "production_runs": [{"run_id": "old", "result_path": "../outside/result.json"}],
    }
    manifest = {"schema_version": 1, "studies": [study]}
    requested = []

    def load_json(relative, label):
        requested.append(relative)
        if relative != "artifacts/legacy/latest/result.json":
            raise AssertionError(f"unexpected remote artifact request: {relative}")
        return {"success": True, "md_result": {"status": "completed", "analytics": {
            "rmsd": {"time_ps": [0, 1000], "backbone_rmsd_angstrom": [0.2, 0.4]},
            "structural_dynamics": {"time_ps": [0, 1000], "pocket": {"retained_fraction": 0.8}},
        }}}

    st = FakeStreamlit()
    render_md_results(st, manifest=manifest, result={}, json_loader=load_json, key_prefix="remote", remote=True)

    assert requested == ["artifacts/legacy/latest/result.json"]
    calls = json.dumps(st.calls, default=str)
    assert "Source analysis history" in calls
    assert "../outside/result.json" in calls  # displayed as provenance only
    assert "Reference-site retention" in calls
    choices = [call for call in st.calls if call[0] == "selectbox" and call[1][0] == "Production result"]
    assert choices and choices[0][1][1][0].startswith("latest")


def test_production_duration_labels_and_failed_continuation_prerequisite(tmp_path):
    run = tmp_path / "run"
    study = _study_manifest()["studies"][0]
    study["production_runs"] = [
        {"run_id": "failed-prereq", "status": "failed", "segment_duration_ns": 1, "previous_duration_ns": 0, "target_duration_ns": 100, "trajectory_appended": False},
        {"run_id": "continuation", "status": "completed", "segment_duration_ns": 99, "previous_duration_ns": 1, "target_duration_ns": 100, "continuation_of_run_id": "failed-prereq", "trajectory_appended": True},
    ]
    study["mmgbsa_runs"] = []
    study["analysis_runs"] = []
    _json(run / "artifacts/md/study.json", {"schema_version": 1, "studies": [study]})
    st = FakeStreamlit()

    render_md_results(st, run, {})

    calls = json.dumps(st.calls, default=str)
    assert "This segment (ns)" in calls and "Previous (ns)" in calls and "Requested total (ns)" in calls
    assert "files support lineage only" in calls
