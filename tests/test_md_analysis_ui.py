import json
from pathlib import Path

from mn_molecular_design_studio import md_analysis_ui


def _write(path: Path, payload: bytes | dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(payload, dict):
        path.write_text(json.dumps(payload), encoding="utf-8")
    else:
        path.write_bytes(payload)
    return path


def test_study_scanner_maps_exact_legacy_paths_and_preserves_segment_lineage(tmp_path):
    source = tmp_path / "aggregate"
    trajectory = _write(source / "artifacts/legacy/prod-99/production.xtc", b"trajectory")
    topology = _write(source / "artifacts/legacy/prep-1/prepared.pdb", b"topology")
    # Same basename in a different source folder must never be substituted for
    # an explicit output path from the production result.
    _write(source / "artifacts/legacy/prod-1/production.xtc", b"wrong trajectory")
    _write(source / "artifacts/legacy/prep-2/prepared.pdb", b"wrong topology")
    _write(source / "artifacts/legacy/prod-99/result.json", {
        "md_result": {"output_files": {
            "production_trajectory": "/output/prod-99/production.xtc",
            "production_pdb": "/prepared-system/prep-1/prepared.pdb",
        }}
    })
    _write(source / "artifacts/md/study.json", {
        "schema_version": 1,
        "studies": [{
            "study_id": "study-a",
            "compound_name": "Ligand A",
            "target_pdb_id": "1ABC",
            "source_pose_id": "pose-7",
            "source_tool": "MN-Ligand",
            "preparation_runs": [{"run_id": "prep-1"}, {"run_id": "prep-2"}],
            "production_runs": [{
                "run_id": "prod-99",
                "status": "completed",
                "result_path": "artifacts/legacy/prod-99/result.json",
                "continuation_of_run_id": "prod-1",
                "previous_duration_ns": 1.0,
                "target_duration_ns": 100.0,
                "segment_duration_ns": 99.0,
                "trajectory_appended": True,
            }],
        }],
    })

    options = md_analysis_ui._study_source_options(source, {"campaign_name": "campaign"})
    assert len(options) == 1
    option = options[0]
    assert option["trajectory_path"] == trajectory
    assert option["topology_paths"] == [topology]
    assert option["lineage_known"] is True
    assert option["segment"] == {
        "run_id": "prod-99",
        "continuation_of_run_id": "prod-1",
        "previous_duration_ns": 1.0,
        "target_duration_ns": 100.0,
        "segment_duration_ns": 99.0,
        "trajectory_appended": True,
        "duration_status": "recorded",
    }
    assert "Ligand A" in option["label"] and "pose-7" in option["label"]


def test_unrecorded_duration_stays_unknown_instead_of_using_target_duration():
    segment = md_analysis_ui._segment_record(
        "run-2", {"target_duration_ns": 100.0, "previous_duration_ns": 1.0}, {}
    )
    assert segment["segment_duration_ns"] is None
    assert segment["duration_status"] == "unknown"
    assert segment["target_duration_ns"] == 100.0


def test_absolute_legacy_output_requires_exact_known_run_id_and_suffix(tmp_path):
    source = tmp_path / "aggregate"
    actual = _write(source / "artifacts/legacy/run-a/nested/traj.xtc", b"exact")
    _write(source / "artifacts/legacy/run-b/nested/traj.xtc", b"same basename")

    resolved = md_analysis_ui._path_candidates_from_output(
        source,
        "/output/run-a/nested/traj.xtc",
        segment_run_id="run-a",
        known_legacy_run_ids={"run-a", "run-b"},
    )
    assert resolved == [actual]
    assert md_analysis_ui._path_candidates_from_output(
        source,
        "/output/not-listed/nested/traj.xtc",
        segment_run_id="run-a",
        known_legacy_run_ids={"run-a", "run-b"},
    ) == []


def test_manifest_with_no_resolvable_segment_does_not_fall_back_to_native(tmp_path):
    source = tmp_path / "aggregate"
    _write(source / "artifacts/md/study.json", {"schema_version": 1, "studies": []})
    _write(source / "artifacts/md/trajectory.xtc", b"not represented by study manifest")
    assert md_analysis_ui._source_options_for_job({"run_dir": str(source)}) == []


def test_missing_recorded_trajectory_does_not_fall_back_to_replica_with_same_basename(tmp_path):
    source = tmp_path / "aggregate"
    other_replica = _write(source / "artifacts/legacy/replica-1/production.xtc", b"other replica")
    topology = _write(source / "artifacts/legacy/prep-1/system.pdb", b"topology")
    _write(source / "artifacts/legacy/replica-2/result.json", {
        "md_result": {"output_files": {
            "production_trajectory": "/output/replica-2/production.xtc",
            "system_pdb": "/prepared-system/prep-1/system.pdb",
        }}
    })
    _write(source / "artifacts/md/study.json", {
        "schema_version": 1,
        "studies": [{
            "study_id": "study-a",
            "preparation_runs": [{"run_id": "prep-1"}],
            "production_runs": [{
                "run_id": "replica-2", "status": "completed",
                "result_path": "artifacts/legacy/replica-2/result.json",
            }],
        }],
    })

    options = md_analysis_ui._study_source_options(source, {})
    assert options == []
    assert other_replica.is_file() and topology.is_file()


def test_normalized_study_paths_avoid_reading_per_segment_result_json(tmp_path, monkeypatch):
    source = tmp_path / "aggregate"
    trajectory = _write(source / "artifacts/legacy/prod-1/production.xtc", b"trajectory")
    topology = _write(source / "artifacts/legacy/prep-1/system.pdb", b"topology")
    _write(source / "artifacts/md/study.json", {
        "schema_version": 1,
        "studies": [{
            "study_id": "study-a", "preparation_runs": [{"run_id": "prep-1"}],
            "production_runs": [{
                "run_id": "prod-1", "status": "completed",
                "result_path": "artifacts/legacy/prod-1/result.json",
                "trajectory_path": "artifacts/legacy/prod-1/production.xtc",
                "topology_path": "artifacts/legacy/prep-1/system.pdb",
                "segment_duration_ns": 1.0,
            }],
        }],
    })
    original_read_json = md_analysis_ui._read_json

    def read_json_without_segment_results(path):
        if Path(path).name == "result.json":
            raise AssertionError("normalized source paths should bypass raw result JSON")
        return original_read_json(Path(path))

    monkeypatch.setattr(md_analysis_ui, "_read_json", read_json_without_segment_results)
    options = md_analysis_ui._study_source_options(source, {})
    assert len(options) == 1
    assert options[0]["trajectory_path"] == trajectory
    assert options[0]["topology_paths"] == [topology]


def test_portable_source_path_accepts_the_managed_owner_job_directory(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine import runtime
    from mn_molecular_design_studio.cluster import config as cluster_config

    managed_root = tmp_path / "results"
    source = managed_root / "ligand-md" / "source-job"
    source.mkdir(parents=True)
    monkeypatch.setattr(runtime, "runs_root", lambda: managed_root)
    monkeypatch.setattr(cluster_config, "enabled", lambda: False)

    assert md_analysis_ui._portable_source_path(source, source, {}) == str(source.resolve())


def test_public_source_options_are_serializable_and_run_relative(tmp_path, monkeypatch):
    import json as json_module
    from mn_molecular_design_studio import auth, job_catalog
    from mn_molecular_design_studio.engine.workflow_engine import runtime

    managed_root = tmp_path / "results"
    source = managed_root / "ligand-md" / "source-job"
    trajectory = _write(source / "artifacts/md/production.xtc", b"trajectory")
    topology = _write(source / "artifacts/input/system.pdb", b"topology")
    _write(source / "result.json", {"md_result": {"output_files": {}}})
    _write(source / "metadata.json", {"run_id": "source-job", "status": "completed"})
    record = {"project_id": "project-a", "owner_id": "user-a"}
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda _group, _run_id: record)
    monkeypatch.setattr(auth, "require_job_access", lambda _user, _record: None)
    monkeypatch.setattr(runtime, "runs_root", lambda: managed_root)
    job = {
        "task_group": "ligand-md", "run_id": "source-job", "project_id": "project-a",
        "run_dir": str(source), "campaign_name": "compound-a", "tool": "OpenMM",
    }

    payload = md_analysis_ui.md_source_options_for_job(job, {"id": "user-a"})
    json_module.dumps(payload)
    assert payload["source_options"] == [{
        "option_id": f"native:{trajectory.relative_to(source).as_posix()}|topology:{topology.relative_to(source).as_posix()}",
        "label": "compound-a · target unknown · pose unknown · OpenMM · segment source-job · production.xtc · topology artifacts/input/system.pdb",
        "study": {"compound_name": "compound-a", "source_tool": "OpenMM"},
        "study_id": None,
        "segment": {
            "run_id": "source-job", "continuation_of_run_id": None,
            "previous_duration_ns": None, "target_duration_ns": None,
            "segment_duration_ns": None, "trajectory_appended": None,
            "duration_status": "unknown",
        },
        "trajectory_relative_path": "artifacts/md/production.xtc",
        "topology_relative_path": "artifacts/input/system.pdb",
        "source_result_path": "result.json",
        "lineage_known": False,
    }]


def test_source_study_identity_requires_exact_manifest_study_id(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.md_analysis import _source_study

    source = tmp_path / "source"
    original = {
        "study_id": "study-a", "compound_id": "CMP-1", "compound_name": "Ligand A",
        "target_pdb_id": "1ABC", "production_runs": [{"run_id": "prod-1"}],
    }
    _write(source / "artifacts/md/study.json", {"schema_version": 1, "studies": [original]})

    resolved = _source_study(source, [{"study_id": "study-a", "run_id": "prod-1"}], "ligand-md/source-job")
    assert resolved["study_id"] == "study-a"
    assert resolved["compound_id"] == "CMP-1"
    assert resolved["production_runs"] == []
    assert resolved["mmgbsa_runs"] == []
    assert resolved["preparation_runs"] == []
    assert resolved["analysis_runs"] == []
    assert resolved["source_history"]["production_runs"] == original["production_runs"]
    assert resolved["source_history"]["source_job"] == {
        "task_group": "ligand-md", "run_id": "source-job", "run_relative_path": "ligand-md/source-job",
    }
    assert _source_study(source, [{"study_id": "unknown-study", "run_id": "prod-1"}], "ligand-md/source-job") is None
    assert _source_study(source, [{"run_id": "prod-1"}], "ligand-md/source-job") is None


def test_analysis_provenance_refs_are_not_operational_paths():
    from mn_molecular_design_studio.api import _is_path_input_key
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import is_operational_path_key

    for key in ("trajectory_ref", "topology_ref", "result_ref"):
        assert not is_operational_path_key(key)
        assert not _is_path_input_key(key)
