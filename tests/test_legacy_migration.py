import csv
import json
from pathlib import Path

from mn_molecular_design_studio.legacy_migration import (
    DATASET_PROJECTS,
    build_migration_plan,
    apply_migration_plan,
    normalize_mn_ligand_run,
)
from mn_molecular_design_studio.molecular_database_ingest import (
    _best_docking_pose_indices,
    _legacy_attempt_target_context,
)
from mn_molecular_design_studio.cofolding_results import best_candidates_by_complex


def _write_run(root: Path, task_group: str, run_id: str, *, metadata: dict, inputs=None, result=None):
    run_dir = root / task_group / run_id
    run_dir.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps({"run_id": run_id, **metadata}), encoding="utf-8")
    if inputs is not None:
        (run_dir / "input.json").write_text(json.dumps(inputs), encoding="utf-8")
    if result is not None:
        (run_dir / "result.json").write_text(json.dumps(result), encoding="utf-8")
    return run_dir


def test_build_plan_follows_dataset_and_target_lineage_without_source_changes(tmp_path):
    source = tmp_path / "old-runs"
    dataset = DATASET_PROJECTS[0]
    _write_run(source, "compound-import", "import-1", metadata={
        "job_type": "compound_import",
        "dataset_name": dataset.source_name,
        "compound_count": 12,
        "status": "completed",
    })
    _write_run(source, "compound-selection", "selection-1", metadata={
        "job_type": "docking_parent_selection",
        "dataset_name": dataset.source_name,
        "parent_run_id": "import-1",
        "status": "completed",
    })
    protein = _write_run(source, "protein-import", "protein-1", metadata={
        "job_type": "protein_import",
        "pdb_id": "4LNW",
        "status": "completed",
    }, result={"imported_target": "artifacts/imported/4lnw.pdb"})
    protein_file = protein / "artifacts" / "imported" / "4lnw.pdb"
    protein_file.parent.mkdir(parents=True)
    protein_file.write_text("ATOM\n", encoding="utf-8")
    _write_run(source, "target-orientation", "prepared-1", metadata={
        "job_type": "target_orientation",
        "pdb_id": "4LNW",
        "parent_run_id": "protein-1",
        "status": "completed",
    })
    docking = _write_run(source, "docking", "dock-1", metadata={
        "job_type": "docking_campaign",
        "tool": "vina",
        "prepared_target_run_id": "prepared-1",
        "compound_run_ids": ["selection-1"],
        "compound_count": 12,
        "replicates": 3,
        "status": "completed",
    })
    _write_run(source, "pose-validation", "validation-1", metadata={
        "job_type": "pose_validation",
        "tool": "PoseBusters",
        "parent_run_id": "dock-1",
        "status": "completed",
    })
    before = {path: path.read_bytes() for path in source.rglob("*.json")}

    plan = build_migration_plan(source, dataset_projects=(dataset,))

    jobs = plan.projects[dataset.project_name]
    assert {job["run_id"] for job in jobs} == {"dock-1", "validation-1"}
    docking_job = next(job for job in jobs if job["run_id"] == "dock-1")
    assert docking_job["receptor"] == "TRα"
    assert docking_job["replicates"] == 3
    assert "import-1" in docking_job["dependency_run_ids"]
    assert plan.targets[0]["source_file"].endswith("artifacts/imported/4lnw.pdb")
    assert all(path.read_bytes() == content for path, content in before.items())
    assert (docking / "metadata.json").is_file()


def test_non_thyroid_target_is_not_selected(tmp_path):
    source = tmp_path / "old-runs"
    dataset = DATASET_PROJECTS[0]
    _write_run(source, "compound-import", "import-1", metadata={
        "job_type": "compound_import", "dataset_name": dataset.source_name,
    })
    _write_run(source, "compound-selection", "selection-1", metadata={
        "job_type": "docking_parent_selection", "dataset_name": dataset.source_name,
        "parent_run_id": "import-1",
    })
    _write_run(source, "structure", "protein-1", metadata={"pdb_id": "4WBK"})
    _write_run(source, "docking", "dock-1", metadata={
        "job_type": "docking_campaign", "tool": "vina", "pdb_id": "4WBK",
        "parent_run_id": "protein-1", "compound_run_ids": ["selection-1"],
        "status": "completed",
    })

    plan = build_migration_plan(source, dataset_projects=(dataset,))

    assert plan.projects[dataset.project_name] == []


def test_apply_migration_plan_copies_runs_and_reuses_existing_groups(tmp_path, monkeypatch):
    from mn_molecular_design_studio.auth import development_user
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    home = tmp_path / "studio-home"
    results = tmp_path / "studio-runs"
    source = tmp_path / "old-runs"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(home))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(home))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(results))
    user = development_user()
    database = MolecularDatabase(root=home / "database")
    for dataset in DATASET_PROJECTS:
        database.create_group(user, dataset.group_names[0], import_method="user_uploads")

    for index, dataset in enumerate(DATASET_PROJECTS, start=1):
        import_id = f"import-{index}"
        selection_id = f"selection-{index}"
        protein_id = f"protein-{index}"
        prepared_id = f"prepared-{index}"
        docking_ids = [f"dock-{index}-a", f"dock-{index}-b"]
        _write_run(source, "compound-import", import_id, metadata={
            "job_type": "compound_import", "dataset_name": dataset.source_name,
            "compound_count": 1, "status": "completed",
        })
        _write_run(source, "compound-selection", selection_id, metadata={
            "job_type": "docking_parent_selection", "dataset_name": dataset.source_name,
            "parent_run_id": import_id, "status": "completed",
        })
        protein = _write_run(source, "protein-import", protein_id, metadata={
            "job_type": "protein_import", "pdb_id": "4LNW", "status": "completed",
        }, result={"imported_target": "artifacts/imported/4lnw.pdb"})
        pdb = protein / "artifacts" / "imported" / "4lnw.pdb"
        pdb.parent.mkdir(parents=True)
        pdb.write_text(
            "ATOM      1  N   ALA A   1      11.104  13.207  12.679  1.00 20.00           N  \nTER\nEND\n",
            encoding="utf-8",
        )
        (protein / "artifacts" / "imported" / "docker-only-path").symlink_to("/boltz-context/processed")
        _write_run(source, "target-orientation", prepared_id, metadata={
            "job_type": "target_orientation", "pdb_id": "4LNW",
            "parent_run_id": protein_id, "status": "completed",
        })
        for index, docking_id in enumerate(docking_ids):
            docking_run = _write_run(source, "docking", docking_id, metadata={
                "job_type": "docking_campaign", "tool": "vina", "pdb_id": "4LNW",
                "prepared_target_run_id": prepared_id, "compound_run_ids": [selection_id],
                "compound_count": 1, "replicates": 3, "status": "completed",
            }, inputs={"tool": "vina"}, result={"success": True})
            compounds_input = docking_run / "input" / "compounds.tsv"
            compounds_input.parent.mkdir(parents=True)
            compounds_input.write_text(
                "compound_id\tcompound_name\tsmiles\n"
                f"compound-{index}\tCompound {index}\tC{index}\n",
                encoding="utf-8",
            )
            (docking_run / "poses").mkdir()
            (docking_run / "poses" / "best.pdbqt").write_text("pose", encoding="utf-8")
            (docking_run / "scores.csv").write_text(
                "compound_id,engine,best_score_kcal_mol,pose_file\n"
                f"compound-{index},vina,-7.5,poses/best.pdbqt\n",
                encoding="utf-8",
            )

    plan = build_migration_plan(source)
    summary = apply_migration_plan(
        plan,
        user=user,
        results_root=results,
        workspace_home=home,
        database=database,
    )

    assert summary["copied_run_count"] == len(plan.runs)
    assert summary["source_result_run_count"] == 2 * len(DATASET_PROJECTS)
    assert summary["result_job_count"] == len(DATASET_PROJECTS)
    assert len(summary["projects"]) == len(DATASET_PROJECTS)
    assert len(summary["targets"]) == 1
    aggregate_index = json.loads((results / "_mn_ligand_migration" / "aggregate_index.json").read_text())
    assert aggregate_index["aggregate_job_count"] == len(DATASET_PROJECTS)
    for dataset in DATASET_PROJECTS:
        group = next(group for group in database.groups(user, scope="all") if group["name"] == dataset.group_names[0])
        assert group["project_id"] == summary["projects"][dataset.project_name]
    for job in plan.projects.values():
        run_id = job[0]["run_id"]
        aggregate_relative = aggregate_index["source_run_to_job"][run_id]["path"]
        migrated = results / aggregate_relative
        aggregate_marker = json.loads((migrated / "mn_ligand_aggregate.json").read_text())
        assert len(aggregate_marker["source_run_ids"]) == 2
        assert all((migrated / "artifacts" / "replicates" / source_run_id / "metadata.json").is_file()
                   for source_run_id in aggregate_marker["source_run_ids"])
        with (migrated / "artifacts" / "docking_scores.csv").open(newline="") as handle:
            score_rows = list(csv.DictReader(handle))
        assert len(score_rows) == 2
        assert {row["legacy_run_id"] for row in score_rows} == set(aggregate_marker["source_run_ids"])
        assert all((migrated / row["pose_path"]).is_file() for row in score_rows)
        with (migrated / "artifacts" / "input" / "compounds.tsv").open(newline="") as handle:
            compound_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(compound_rows) == 2
        assert all(row["smiles"] for row in compound_rows)
        marker = json.loads((migrated / "molecular_design_studio.json").read_text())
        assert marker["project_id"] in summary["projects"].values()
    migrated_link = results / "_mn_ligand_migration" / "dependencies" / "protein-import" / "protein-1" / "artifacts" / "imported" / "docker-only-path"
    assert migrated_link.is_symlink()
    assert migrated_link.readlink() == Path("/boltz-context/processed")


def test_best_structure_selection_preserves_each_legacy_attempt(tmp_path):
    rows = [
        {"compound_id": "A", "legacy_run_id": "attempt-a", "best_score_kcal_mol": "-7.0"},
        {"compound_id": "A", "legacy_run_id": "attempt-a", "best_score_kcal_mol": "-8.0"},
        {"compound_id": "A", "legacy_run_id": "attempt-b", "best_score_kcal_mol": "-6.0"},
        {"compound_id": "A", "legacy_run_id": "attempt-b", "best_score_kcal_mol": "-9.0"},
    ]
    pose_paths = []
    for index in range(len(rows)):
        pose = tmp_path / f"pose-{index}.pdbqt"
        pose.write_text("pose", encoding="utf-8")
        pose_paths.append(pose)

    assert _best_docking_pose_indices(rows, pose_paths) == {1, 3}

    predictions = [
        {"candidate_id": "a1", "metrics": {"confidence_score": 0.7},
         "raw_metadata": {"legacy_run_id": "attempt-a", "source_compound_id": "A"}},
        {"candidate_id": "a2", "metrics": {"confidence_score": 0.9},
         "raw_metadata": {"legacy_run_id": "attempt-a", "source_compound_id": "A"}},
        {"candidate_id": "b1", "metrics": {"confidence_score": 0.8},
         "raw_metadata": {"legacy_run_id": "attempt-b", "source_compound_id": "A"}},
    ]
    best = best_candidates_by_complex(predictions)
    assert {candidate["candidate_id"] for candidate in best} == {"a2", "b1"}


def test_legacy_attempt_target_context_uses_results_copy_and_retains_target_provenance(tmp_path):
    aggregate = tmp_path / "runs" / "ligand-docking" / "aggregate-1"
    attempt_id = "attempt-1"
    attempt = aggregate / "artifacts" / "replicates" / attempt_id
    receptor = attempt / "input" / "receptor.pdb"
    receptor.parent.mkdir(parents=True)
    receptor.write_text("ATOM      1  N   ALA A   1      11.104  13.207  12.679  1.00 20.00           N\nEND\n")
    import hashlib
    digest = hashlib.sha256(receptor.read_bytes()).hexdigest()
    (attempt / "input.json").write_text(json.dumps({
        "target_artifact": {
            "run_id": "prepared-target-1",
            "artifact_id": "artifact-receptor-1",
            "artifact_type": "prepared_target",
            "role": "axis_aligned_receptor",
            "sha256": digest,
            "metadata": {"coordinate_transform": {"rotation": "identity"}},
        },
    }))

    context = _legacy_attempt_target_context(
        aggregate,
        attempt_id,
        {"pdb_id": "3GWS", "campaign_name": "TRβ · 3GWS · HY-L126 · vina"},
    )

    assert context["target_path"] == str(receptor.resolve())
    assert context["target_source_attempt_run_id"] == attempt_id
    assert context["prepared_target_run_id"] == "prepared-target-1"
    assert context["target_receptor"] == "TRβ"
    assert context["target_pdb_id"] == "3GWS"
    assert context["target_artifact_sha256_matches_input"] is True
    assert context["target_coordinate_transform"] == {"rotation": "identity"}


def test_normalize_mn_ligand_run_keeps_replicates_and_original_affinity_table(tmp_path):
    run_dir = tmp_path / "docking" / "legacy-run"
    run_dir.mkdir(parents=True)
    (run_dir / "mn_ligand_migration.json").write_text(json.dumps({
        "source_application": "MN-Ligand", "role": "result",
    }), encoding="utf-8")
    (run_dir / "results" / "replicate_001").mkdir(parents=True)
    (run_dir / "results" / "replicate_001" / "A_out.pdbqt").write_text("pose", encoding="utf-8")
    source_scores = (
        "compound_id,engine,replicate,best_score_kcal_mol,pose_file\n"
        "A,vina,1,-7.5,results/replicate_001/A_out.pdbqt\n"
    ).encode()
    (run_dir / "scores.csv").write_bytes(source_scores)
    artifacts = run_dir / "artifacts"
    artifacts.mkdir()
    source_affinity = (
        "candidate_id,replicate,pIC50,affinity_log10_ic50_uM\n"
        "B,1,6.2,-0.2\nB,2,6.4,-0.4\n"
    ).encode()
    (artifacts / "nesso_affinity.csv").write_bytes(source_affinity)
    (artifacts / "boltz2_metrics.csv").write_text(
        "candidate_id,compound_id,model_id,replicate,seed,confidence_score,iptm,structure_file\n"
        "C,C,C_model_0,1,1001,0.8,0.7,models/r1.cif\n"
        "C,C,C_model_1,1,1001,0.9,0.9,models/r1b.cif\n"
        "C,C,C_model_0,2,1002,0.7,0.6,models/r2.cif\n"
        "C,C,C_model_1,2,1002,0.8,0.8,models/r2b.cif\n",
        encoding="utf-8",
    )
    (run_dir / "models").mkdir()
    for name in ("r1.cif", "r1b.cif", "r2.cif", "r2b.cif"):
        (run_dir / "models" / name).write_text("structure", encoding="utf-8")
    compounds = {
        key: {"compound_id": key, "compound_name": key, "smiles": "CC"}
        for key in ("A", "B", "C")
    }

    summary = normalize_mn_ligand_run(
        run_dir, project_name=DATASET_PROJECTS[0].project_name, compounds=compounds,
    )

    assert source_scores == (run_dir / "scores.csv").read_bytes()
    assert source_affinity == (artifacts / "mn_ligand_source" / "nesso_affinity.csv").read_bytes()
    with (artifacts / "docking_scores.csv").open(newline="") as handle:
        docking_rows = list(csv.DictReader(handle))
    assert docking_rows[0]["pose_path"] == "results/replicate_001/A_out.pdbqt"
    with (artifacts / "nesso_affinity.csv").open(newline="") as handle:
        affinity_rows = list(csv.DictReader(handle))
    assert [row["compound_id"] for row in affinity_rows] == ["B", "B"]
    candidates = [json.loads(line) for line in (artifacts / "normalized_candidates" / "candidates.jsonl").read_text().splitlines()]
    assert len(candidates) == 2
    assert [row["metrics"]["replicate"] for row in candidates] == [1.0, 2.0]
    assert round(candidates[0]["metrics"]["confidence_score"], 2) == 0.85
    assert candidates[0]["complex_pdb"] == "models/r1b.cif"
    with (artifacts / "input" / "compounds.tsv").open(newline="") as handle:
        input_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["compound_id"] for row in input_rows} == {"A", "B", "C"}
    assert summary["cofold_candidate_count"] == 2
    assert summary["unresolved_compound_ids"] == 0
