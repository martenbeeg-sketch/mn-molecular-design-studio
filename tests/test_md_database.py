from __future__ import annotations


def test_production_details_retains_explicit_step_derived_segment_duration(tmp_path):
    from mn_molecular_design_studio.md_database import _production_details

    details = _production_details(tmp_path, {
        "status": "completed", "segment_duration_ns": 50.0,
        "segment_duration_basis": "completed production_steps × production_timestep_fs",
    })
    assert details["actual_duration_ns"] == 50.0
    assert details["duration_basis"] == "completed production_steps × production_timestep_fs"

import hashlib
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.md_database import ingest_completed_md_run, ingest_md_studies
from mn_molecular_design_studio.molecular_database import MolecularDatabase


OWNER = {"id": "a" * 32, "role": "user"}


@pytest.fixture
def database(tmp_path):
    return MolecularDatabase(tmp_path / "database")


def _identity(database, smiles="CCO"):
    imported = database.import_records(OWNER, molecule_imports.parse_smiles(f"{smiles} ligand"), name="known ligands")
    return imported["group_id"], imported["molecule_ids"][0]


def _write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True))
    return path


def _endpoint(*, success=True, units="kJ/mol", value_kj=4.184, value_kcal=1.0, status="success"):
    return {
        "success": success,
        "mmgbsa": {
            "status": status,
            "method": "ambertools_mmpbsa_serial",
            "units": units,
            "units_secondary": "kcal/mol",
            "metadata": {"n_frames_analyzed": 20},
            "gb": {"delta": {
                "delta_g_bind_total_kj_mol": value_kj,
                "delta_g_bind_total_kcal_mol": value_kcal,
                "delta_mm_kcal_mol": -2.0,
            }},
            "pb": {"delta": {"delta_g_bind_total_kcal_mol": -3.0}},
        },
    }


def _study(group_id, molecule_id, *, result_path="artifacts/md/endpoint.json", identity_status="verified_source_dataset", smiles="CCO"):
    return {
        "study_id": "study-1",
        "compound_id": "dataset-row-8",
        "compound_name": "ligand",
        "group_id": group_id,
        "molecule_id": molecule_id,
        "identity_status": identity_status,
        "ligand_smiles": smiles,
        "identity_smiles": "CCO",
        "target_pdb_id": "1ABC",
        "target_receptor": "PDB chain A",
        "target_provenance": {"source": "dataset-v2", "sha256": "target-hash"},
        "source_pose_id": "pose-8",
        "source_tool": "GNINA",
        "source_pose_rank": 8,
        "selection_status": "manual_nonbest",
        "selected_complex_run_id": "complex-8",
        "md_engine": "GROMACS",
        "production_runs": [{
            "run_id": "prod-1", "status": "completed", "repeat_index": 1,
            "continuation_of_run_id": "prod-0", "target_duration_ns": 100,
            "result_path": "artifacts/md/prod.json",
        }],
        "mmgbsa_runs": [{
            "run_id": "endpoint-1", "status": "completed", "source_production_run_id": "prod-1",
            "result_path": result_path, "backend": "ambertools_mmpbsa",
            "source_md_job": {"task_group": "ligand-md", "run_id": "studio-md", "workstation": "pc004033"},
            "legacy_source_md_job": {"task_group": "bound-ligand-md", "run_id": "prod-1", "workstation": "pc004033"},
            "source_trajectory_ref": "cluster://node/ligand-md/studio-md/artifacts/legacy/prod-1/traj.dcd",
            "recorded_trajectory_sha256": "a" * 64, "copied_trajectory_sha256": "a" * 64,
            "trajectory_snapshot_matches": True,
            "parameters": {"backend": "ambertools_mmpbsa", "stride": 10, "start_pct": 25, "end_pct": 90},
        }],
        "analysis_runs": [{"run_id": "analysis-1", "status": "completed", "result_path": "artifacts/md/analysis.json"}],
    }


def _run(tmp_path, study, *, owner_id=OWNER["id"]):
    run = tmp_path / "runs" / "md-study" / "aggregate-1"
    run.mkdir(parents=True, exist_ok=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": owner_id, "project_id": ""}))
    (run / "metadata.json").write_text(json.dumps({"owner_id": owner_id, "tool": "md_study"}))
    _write_json(run / "artifacts/md/study.json", {"schema_version": 1, "studies": [study]})
    return run


def _observations(database):
    with database.connection() as connection:
        return [dict(row) for row in connection.execute("SELECT * FROM molecule_observations ORDER BY namespace")]


def test_ingest_is_idempotent_and_preserves_exact_verified_identity_and_provenance(database, tmp_path):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id)
    endpoint_path = _write_json(tmp_path / "runs/md-study/aggregate-1/artifacts/md/endpoint.json", _endpoint())
    analysis_path = _write_json(tmp_path / "runs/md-study/aggregate-1/artifacts/md/analysis.json", {
        "success": True,
        "replicas": [{"replica": 1, "run_id": "prod-1", "status": "completed",
                      "backbone_rmsd_tail_mean_angstrom": 1.3, "ligand_rmsd_max_angstrom": 2.4,
                      "delta_g_bind_kcal_mol": -999}],
        "contact_matrix": {"residues": ["A:42", "B:17"], "replicas": [1], "contact_occupancy": [[0.75], [float("nan")]]},
    })
    _write_json(tmp_path / "runs/md-study/aggregate-1/artifacts/md/prod.json", {
        "success": True,
        "md_result": {"equilibration_stats": {"production": {"duration_ns": 99.0, "trajectory_frames": 120}}},
    })
    run = _run(tmp_path, study)

    first = ingest_md_studies(run, database=database)
    observations_first = _observations(database)
    second = ingest_md_studies(run, database=database)
    observations_second = _observations(database)

    assert first["status"] == second["status"] == "ok"
    assert len(observations_first) == len(observations_second) == 4
    by_namespace = {row["namespace"]: json.loads(row["metrics"]) for row in observations_second}
    assert by_namespace["md.mmgbsa.gb"]["delta_g_bind_kcal_mol"] == pytest.approx(1.0)
    assert by_namespace["md.mmgbsa.gb"]["delta_mm_kcal_mol"] == -2
    assert by_namespace["md.mmgbsa.pb"]["delta_g_bind_kcal_mol"] == -3
    assert by_namespace["md.analysis"]["backbone_rmsd_tail_mean_angstrom"] == pytest.approx(1.3)
    assert "contact_occupancy.residue_1_A_42" not in by_namespace["md.analysis"]
    assert len(by_namespace["md.analysis"]) == 2
    assert "delta_g_bind_kcal_mol" not in by_namespace["md.analysis"]
    assert by_namespace["md.production"]["production_actual_duration_ns"] == 99
    assert by_namespace["md.production"]["production_trajectory_frames"] == 120
    production_context = json.loads(next(row["context"] for row in observations_second if row["namespace"] == "md.production"))
    assert production_context["requested_target_duration_ns"] == 100
    assert production_context["continuation_of_run_id"] == "prod-0"
    endpoint_context = json.loads(next(row["context"] for row in observations_second if row["namespace"] == "md.mmgbsa.gb"))
    assert endpoint_context["source_summary_sha256"] == hashlib.sha256(endpoint_path.read_bytes()).hexdigest()
    assert endpoint_context["source_production_run_id"] == "prod-1"
    assert endpoint_context["source_pose_id"] == "pose-8"
    assert endpoint_context["selection_status"] == "manual_nonbest"
    assert endpoint_context["source_md_job"]["run_id"] == "studio-md"
    assert endpoint_context["legacy_source_md_job"]["run_id"] == "prod-1"
    assert endpoint_context["source_trajectory_ref"].startswith("cluster://node/ligand-md/studio-md/artifacts/legacy/")
    assert endpoint_context["target_provenance"]["source"] == "dataset-v2"
    assert endpoint_context["is_representative"] is True
    assert endpoint_context["type"] == "job"
    assert endpoint_context["run_id"] == "aggregate-1"
    assert endpoint_context["task_group"] == "md-study"
    assert endpoint_context["tool"] == "md_study"
    assert endpoint_context["workstation"]
    assert endpoint_context["target_id"] == "PDB chain A:1ABC:target-hash"
    assert endpoint_context["identity_match_basis"] == "source_dataset_identity_smiles_and_existing_membership"
    assert endpoint_context["ligand_smiles"] == "CCO"
    analysis_context = json.loads(next(row["context"] for row in observations_second if row["namespace"] == "md.analysis"))
    assert analysis_context["contact_map_summary_relative_path"] == "artifacts/md/analysis.json"
    assert analysis_context["contact_map_summary_sha256"] == hashlib.sha256(analysis_path.read_bytes()).hexdigest()
    assert not any("structures" in row for row in observations_second)


def test_identity_must_be_explicit_verified_and_a_group_member(database, tmp_path):
    group_id, molecule_id = _identity(database)
    unverified = _run(tmp_path / "unverified", _study(group_id, molecule_id, identity_status="inferred"))
    missing = _run(tmp_path / "missing", {**_study(group_id, molecule_id), "molecule_id": ""})
    conflicting_smiles = _run(tmp_path / "conflict", {**_study(group_id, molecule_id), "identity_smiles": "CCN"})
    outside_member = _run(tmp_path / "nonmember", _study(group_id, "mol-does-not-exist"))

    reports = [ingest_md_studies(path, database=database) for path in (unverified, missing, conflicting_smiles, outside_member)]
    assert all(report["studies"][0]["status"] == "withheld" for report in reports)
    assert all(not _observations(database) for _ in [0])


def test_inaccessible_group_is_withheld_for_caller_and_trusted_owner_is_derived(database, tmp_path):
    group_id, molecule_id = _identity(database)
    run = _run(tmp_path, _study(group_id, molecule_id))
    _write_json(run / "artifacts/md/endpoint.json", _endpoint())
    stranger = {"id": "b" * 32, "role": "user"}

    denied = ingest_md_studies(run, database=database, user=stranger)
    assert denied["studies"][0]["status"] == "withheld"
    assert not _observations(database)
    accepted = ingest_md_studies(run, database=database)
    assert accepted["observations"] == 2  # GB and PB endpoint records only.


def test_failed_superseded_and_external_results_are_never_ingested(database, tmp_path):
    group_id, molecule_id = _identity(database)
    run = tmp_path / "runs/md-study/aggregate-1"
    run.mkdir(parents=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": OWNER["id"]}))
    (run / "metadata.json").write_text(json.dumps({"owner_id": OWNER["id"]}))
    failed_path = _write_json(run / "artifacts/md/failed.json", _endpoint(success=False))
    old_endpoint_path = _write_json(run / "artifacts/md/old.json", _endpoint())
    outside = _write_json(tmp_path / "outside.json", _endpoint())
    study = _study(group_id, molecule_id)
    study["mmgbsa_runs"] = [
        {"run_id": "bad", "status": "failed", "source_production_run_id": "prod-1", "result_path": "artifacts/md/failed.json", "parameters": {"backend": "ambertools_mmpbsa"}},
        {"run_id": "external", "status": "completed", "result_path": str(outside), "parameters": {"backend": "ambertools_mmpbsa"}},
        {"run_id": "old", "status": "completed", "source_production_run_id": "prod-1", "result_path": "artifacts/md/old.json", "superseded_by_run_id": "new", "parameters": {"backend": "ambertools_mmpbsa"}},
    ]
    study["production_runs"] = []
    study["analysis_runs"] = []
    _write_json(run / "artifacts/md/study.json", {"schema_version": 1, "studies": [study]})

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 2
    assert len(report["studies"][0]["withheld"]) == 2
    rows = _observations(database)
    assert {row["namespace"] for row in rows} == {"md.mmgbsa.gb", "md.mmgbsa.pb"}
    context = json.loads(rows[0]["context"])
    assert context["endpoint_superseded_by_run_id"] == "new"
    assert context["is_representative"] is False
    assert failed_path.is_file() and old_endpoint_path.is_file() and outside.is_file()


@pytest.mark.parametrize(
    ("endpoint", "reason"),
    [
        (_endpoint(units="mystery units"), "unknown energy units"),
        (_endpoint(value_kj=8.368, value_kcal=1.0), "inconsistent kcal/kJ duplicate"),
    ],
)
def test_unknown_units_and_inconsistent_unit_duplicates_are_withheld(database, tmp_path, endpoint, reason):
    group_id, molecule_id = _identity(database)
    run = _run(tmp_path, _study(group_id, molecule_id))
    _write_json(run / "artifacts/md/endpoint.json", endpoint)

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 0
    assert reason in report["studies"][0]["withheld"][0]["reason"]
    assert not _observations(database)


def test_manifest_and_result_paths_cannot_escape_source_run(database, tmp_path):
    group_id, molecule_id = _identity(database)
    run = _run(tmp_path, _study(group_id, molecule_id, result_path="../outside.json"))

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 0
    assert "safe relative" in report["studies"][0]["withheld"][0]["reason"]
    assert not _observations(database)


def test_modeling_smiles_is_not_used_as_source_identity(database, tmp_path):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id, smiles="CCN")
    run = _run(tmp_path, study)
    _write_json(run / "artifacts/md/endpoint.json", _endpoint())

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 2
    context = json.loads(_observations(database)[0]["context"])
    assert context["identity_smiles"] == "CCO"
    assert context["ligand_smiles"] == "CCN"


def test_analysis_status_guards_and_missing_stage_paths_are_silent(database, tmp_path):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id)
    study["mmgbsa_runs"].append({"run_id": "stage-only", "status": "completed"})
    study["analysis_runs"] = [
        {"run_id": "analysis-stage-only", "status": "completed"},
        {"run_id": "analysis-1", "status": "completed", "result_path": "artifacts/md/analysis.json"},
    ]
    study["production_runs"] = []
    run = _run(tmp_path, study)
    _write_json(run / "artifacts/md/endpoint.json", _endpoint())
    _write_json(run / "artifacts/md/analysis.json", {
        "success": True,
        "replicas": [
            {"replica": 1, "status": "failed", "backbone_rmsd_tail_mean_angstrom": 1.0},
            {"replica": 2, "status": "completed", "backbone_rmsd_tail_mean_angstrom": 2.0},
        ],
    })

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 3
    rows = _observations(database)
    assert len([row for row in rows if row["namespace"] == "md.analysis"]) == 1
    analysis_context = json.loads(next(row["context"] for row in rows if row["namespace"] == "md.analysis"))
    assert analysis_context["replica_index"] == 2
    assert not report["studies"][0]["withheld"]


def test_legacy_single_delta_mmgbsa_imports_gb_only_when_method_is_explicit(database, tmp_path):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id)
    study["analysis_runs"] = []
    study["production_runs"] = []
    study["mmgbsa_runs"] = [{
        "run_id": "legacy-endpoint", "status": "completed",
        "method": "single-trajectory_openmm_mmgbsa_script",
        "result_path": "artifacts/md/legacy.json",
    }]
    run = _run(tmp_path, study)
    _write_json(run / "artifacts/md/legacy.json", {
        "success": True,
        "mmgbsa": {"status": "success", "method": "single-trajectory_openmm_mmgbsa_script",
                   "units": "kcal/mol", "delta": {"g_bind_total": -12.5}},
    })

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 1
    row = _observations(database)[0]
    assert row["namespace"] == "md.mmgbsa.gb"
    assert json.loads(row["metrics"])["delta_g_bind_kcal_mol"] == -12.5
    assert json.loads(row["context"])["method"] == "single-trajectory_openmm_mmgbsa_script"


def test_source_dataset_membership_and_local_evaluation_group_are_both_required(database, tmp_path):
    source_group, molecule_id = _identity(database)
    local_group = database.create_group(OWNER, "MD evaluations")
    with database.connection() as connection:
        row = connection.execute("SELECT name,components,metadata FROM molecule_members WHERE group_id=? AND molecule_id=?", (source_group, molecule_id)).fetchone()
        connection.execute("INSERT INTO molecule_members VALUES (?,?,?,?,?,?)", (local_group, molecule_id, row["name"], row["components"], row["metadata"], "test"))
    study = _study(local_group, molecule_id)
    study["source_dataset_group_id"] = source_group
    run = _run(tmp_path, study)
    _write_json(run / "artifacts/md/endpoint.json", _endpoint())

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 2
    context = json.loads(_observations(database)[0]["context"])
    assert context["group_id"] == local_group
    assert context["source_dataset_group_id"] == source_group


def _direct_summary(study_id):
    return {
        "schema_version": 1,
        "success": True,
        "result_type": "md_trajectory_analysis_result/v1",
        "source_study_id": study_id,
        "metrics": {
            "backbone_rmsd_mean_angstrom": {"value": 1.25, "unit": "angstrom"},
            "backbone_rmsd_max_angstrom": {"value": 2.4, "unit": "angstrom"},
            "ligand_rmsd_mean_angstrom": {"value": 2.1, "unit": "angstrom"},
            "selected_frame_count": {"value": 80, "unit": "count"},
            "contact_count_mean": {"value": 3.5, "unit": "count"},
            "delta_g_bind_kcal_mol": -1000,
            "arbitrary_residue_contact_occupancy": 0.99,
        },
        "selected_time_range_ns": {"start": 5.0, "end": 25.0},
        "source": {
            "topology": {"path": "artifacts/input/topol.tpr", "sha256": "a" * 64, "size_bytes": 123},
            "trajectory": {"path": "artifacts/input/traj.xtc", "sha256": "b" * 64, "size_bytes": 456},
        },
        "source_segments": [{"run_id": "prod-1", "continuation_of_run_id": "prod-0"}],
        "parameters": {"stride": 10, "contact_cutoff_angstrom": 4.5, "pbc_mode_requested": "auto", "pbc_imaging_applied": True},
        "selection": {"ligand": "resname LIG", "protein": "protein"},
        "artifacts": {"residue_contacts_csv": "artifacts/md_analysis/residue_contacts.csv"},
    }


def test_direct_mdtraj_summary_imports_curated_metrics_and_retains_hash_lineage(database, tmp_path):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id)
    study["production_runs"] = []
    study["mmgbsa_runs"] = []
    study["analysis_runs"] = [{
        "run_id": "analysis-job", "status": "completed",
        "result_type": "md_trajectory_analysis_result/v1",
        "result_path": "artifacts/md_analysis/summary.json",
        "source_created_at": "2025-03-01T00:00:00Z", "source_completed_at": "2025-03-02T00:00:00Z",
    }]
    run = _run(tmp_path, study)
    _write_json(run / "artifacts/md_analysis/summary.json", _direct_summary("study-1"))
    contacts = run / "artifacts/md_analysis/residue_contacts.csv"
    contacts.parent.mkdir(parents=True, exist_ok=True)
    contacts.write_text("residue,contact_occupancy_fraction\nA:42,0.75\n")

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 1
    row = _observations(database)[0]
    metrics = json.loads(row["metrics"])
    assert metrics == {
        "backbone_rmsd_mean_angstrom": 1.25,
        "backbone_rmsd_max_angstrom": 2.4,
        "ligand_rmsd_mean_angstrom": 2.1,
        "selected_frame_count": 80,
        "contact_count_mean": 3.5,
        "selected_duration_ns": 20,
    }
    context = json.loads(row["context"])
    assert context["analysis_result_type"] == "md_trajectory_analysis_result/v1"
    assert context["source_study_id"] == "study-1"
    assert context["topology_sha256"] == "a" * 64
    assert context["trajectory_sha256"] == "b" * 64
    assert context["selected_time_range_ns"] == {"start": 5.0, "end": 25.0}
    assert context["analysis_parameters"]["pbc_mode_requested"] == "auto"
    assert context["settings"]["time_window_ns"] == {"start": 5.0, "end": 25.0}
    assert context["source_created_at"] == "2025-03-01T00:00:00Z"
    assert context["source_completed_at"] == "2025-03-02T00:00:00Z"
    assert row["created_at"] == "2025-03-02T00:00:00Z"
    assert context["source_segments"][0]["continuation_of_run_id"] == "prod-0"
    assert context["contact_map_artifact_relative_path"] == "artifacts/md_analysis/residue_contacts.csv"
    assert context["contact_map_artifact_sha256"] == hashlib.sha256(contacts.read_bytes()).hexdigest()


@pytest.mark.parametrize(
    ("success", "source_study_id", "reason"),
    [(False, "study-1", "did not report success"), (True, "another-study", "does not match this study")],
)
def test_direct_summary_requires_success_and_matching_source_study(database, tmp_path, success, source_study_id, reason):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id)
    study["production_runs"] = []
    study["mmgbsa_runs"] = []
    study["analysis_runs"] = [{
        "run_id": "analysis-job", "status": "completed",
        "result_type": "md_trajectory_analysis_result/v1",
        "result_path": "artifacts/md_analysis/summary.json",
    }]
    run = _run(tmp_path, study)
    summary = _direct_summary(source_study_id)
    summary["success"] = success
    _write_json(run / "artifacts/md_analysis/summary.json", summary)

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 0
    assert len(report["studies"][0]["withheld"]) == 1
    assert reason in report["studies"][0]["withheld"][0]["reason"]
    assert not _observations(database)


def test_completed_md_run_cache_skips_unchanged_summary(database, tmp_path):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id)
    study["production_runs"] = []
    study["mmgbsa_runs"] = []
    study["analysis_runs"] = [{
        "run_id": "analysis-job", "status": "completed",
        "result_type": "md_trajectory_analysis_result/v1",
        "result_path": "artifacts/md_analysis/summary.json",
    }]
    run = _run(tmp_path, study)
    (run / "metadata.json").write_text(json.dumps({"owner_id": OWNER["id"], "status": "completed", "tool": "md_analysis"}))
    _write_json(run / "result.json", {"success": True})
    _write_json(run / "artifacts/md_analysis/summary.json", _direct_summary("study-1"))

    first = ingest_completed_md_run(run, database=database)
    second = ingest_completed_md_run(run, database=database)

    assert first == 1
    assert second == 0
    assert len(_observations(database)) == 1


def test_original_stage_dates_and_protocol_settings_drive_observation_order_and_digest(database, tmp_path):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id)
    study["production_runs"] = [
        {"run_id": "prod-0", "status": "completed", "repeat_index": 1,
         "target_duration_ns": 50, "result_path": "artifacts/md/prod-0.json",
         "source_created_at": "2025-01-01T00:00:00Z", "source_completed_at": "2025-01-02T00:00:00Z"},
        {"run_id": "prod-1", "status": "completed", "repeat_index": 1,
         "continuation_of_run_id": "prod-0", "target_duration_ns": 100,
         "result_path": "artifacts/md/prod.json", "source_created_at": "2025-01-03T00:00:00Z",
         "source_completed_at": "2025-01-05T00:00:00Z"},
    ]
    study["mmgbsa_runs"] = [
        {"run_id": "endpoint-old", "status": "completed", "source_production_run_id": "prod-1",
         "result_path": "artifacts/md/endpoint.json", "method": "ambertools_mmpbsa_serial",
         "parameters": {"stride": 10, "start_pct": 20}, "source_created_at": "2025-01-06T00:00:00Z",
         "source_completed_at": "2025-01-07T00:00:00Z"},
        {"run_id": "endpoint-new", "status": "completed", "source_production_run_id": "prod-1",
         "result_path": "artifacts/md/endpoint.json", "method": "ambertools_mmpbsa_serial",
         "parameters": {"stride": 5, "start_pct": 20}, "source_created_at": "2025-01-08T00:00:00Z",
         "source_completed_at": "2025-01-09T00:00:00Z"},
    ]
    study["analysis_runs"] = []
    run = _run(tmp_path, study)
    _write_json(run / "artifacts/md/endpoint.json", _endpoint())
    _write_json(run / "artifacts/md/prod-0.json", {
        "success": True, "md_result": {"status": "completed", "equilibration_stats": {"production": {"duration_ns": 50, "trajectory_frames": 25}}},
    })
    _write_json(run / "artifacts/md/prod.json", {
        "success": True, "md_result": {"status": "completed", "equilibration_stats": {"production": {"duration_ns": 99, "trajectory_frames": 120}}},
    })
    (run / "metadata.json").write_text(json.dumps({
        "owner_id": OWNER["id"], "tool": "md_study", "created_at": "2025-02-01T00:00:00Z",
        "completed_at": "2025-02-02T00:00:00Z",
    }))

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 6  # Two endpoint branches per endpoint, plus two production segments.
    observations = _observations(database)
    endpoints = [row for row in observations if row["namespace"] == "md.mmgbsa.gb"]
    assert {row["created_at"] for row in endpoints} == {"2025-01-07T00:00:00Z", "2025-01-09T00:00:00Z"}
    first_context = json.loads(next(row["context"] for row in endpoints if json.loads(row["context"])["endpoint_run_id"] == "endpoint-old"))
    assert first_context["source_created_at"] == "2025-01-06T00:00:00Z"
    assert first_context["source_completed_at"] == "2025-01-07T00:00:00Z"
    assert first_context["source_production_status"] == "completed"
    assert first_context["source_production_segment_duration_ns"] == 99
    assert first_context["source_production_previous_duration_ns"] == 50
    assert first_context["source_production_requested_target_duration_ns"] == 100
    assert first_context["settings"]["sampling_parameters"]["stride"] == 10
    assert first_context["settings"]["segment_duration_ns"] == 99

    from mn_molecular_design_studio.evaluations import rebuild_missing
    with database.connection() as connection:
        rebuild_missing(connection)
        descriptors = [json.loads(row[0]) for row in connection.execute(
            "SELECT e.descriptor FROM molecule_evaluations e JOIN molecule_observations o ON o.id=e.observation_id WHERE o.namespace='md.mmgbsa.gb'"
        )]
    assert len({item["parameters_digest"] for item in descriptors}) == 2


def test_failed_or_cancelled_production_result_does_not_publish_duration(database, tmp_path):
    group_id, molecule_id = _identity(database)
    study = _study(group_id, molecule_id)
    study["mmgbsa_runs"] = []
    study["analysis_runs"] = []
    study["production_runs"] = [{
        "run_id": "failed-prod", "status": "failed", "result_path": "artifacts/md/failed-prod.json",
        "target_duration_ns": 100,
    }]
    run = _run(tmp_path, study)
    _write_json(run / "artifacts/md/failed-prod.json", {
        "success": False, "md_result": {"status": "failed", "equilibration_stats": {"production": {"duration_ns": 99, "trajectory_frames": 100}}},
    })

    report = ingest_md_studies(run, database=database)

    assert report["observations"] == 0
    assert not _observations(database)
