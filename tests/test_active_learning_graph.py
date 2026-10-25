"""The study drawing is an iteration, with convergence before retraining."""
import copy
import json

import pytest

from mn_molecular_design_studio import pipeline_graph as graph
from mn_molecular_design_studio.active_learning_settings import ActiveLearningSettings


def _pairs(study):
    kinds = {n["id"]: n["data"]["kind"] for n in study["nodes"]}
    return {(kinds[e["source"]], kinds[e["target"]]) for e in study["edges"]}


def test_stopping_decision_precedes_retraining_and_selected_batch_returns_to_docking():
    study = graph.new_graph("unidock_active_learning")
    plan = graph.validate_graph(study)
    pairs = _pairs(study)
    assert ("active_accumulate", "active_decision") in pairs
    assert ("active_decision", "active_binder_rule") in pairs
    assert ("active_decision", "active_results") in pairs
    assert ("active_library", "active_predict") in pairs
    assert ("active_library", "active_select") in pairs
    assert ("active_library", "active_seed") not in pairs
    assert ("active_select", "active_docking") in pairs
    assert ("active_select", "active_decision") not in pairs
    assert ("active_decision", "active_docking") not in pairs
    by_kind = {n["data"]["kind"]: n["id"] for n in study["nodes"]}
    assert plan["order"].index(by_kind["active_decision"]) < plan["order"].index(by_kind["active_train"])
    feedback = [e for e in study["edges"] if e.get("loop_back")]
    assert len(feedback) == 1
    assert feedback[0]["source"] == by_kind["active_select"]


def test_legacy_study_migration_keeps_selected_inputs_settings_and_ids():
    study = graph.new_graph("unidock_active_learning")
    study.pop("template_version")
    study["name"] = "My study"
    for node in study["nodes"]:
        node["position"] = {"x": 7, "y": 9}
        if node["data"]["kind"] == "active_library":
            node["data"]["params"]["library_id"] = "registered-million"
        if node["data"]["kind"] == "active_binder_rule":
            node["data"]["params"]["binder_percent_target"] = 0.02
            node["data"]["label"] = "Custom binder schedule"
        if node["data"]["kind"] == "active_docking":
            node["data"]["params"]["batch_size"] = 25_000
    before = copy.deepcopy(study)
    migrated = graph.migrate_template_graph(study)
    assert study == before
    assert migrated["id"] == study["id"]
    assert migrated["name"] == "My study"
    assert [n["id"] for n in migrated["nodes"]] == [n["id"] for n in study["nodes"]]
    for old, new in zip(study["nodes"], migrated["nodes"]):
        assert old["data"] == new["data"]
        position = graph.ACTIVE_LEARNING_LAYOUT[new["data"]["kind"]]
        assert new["position"] == {axis: position[axis] for axis in ("x", "y")}
    assert graph.migrate_template_graph(migrated) is migrated
    graph.validate_graph(migrated)


def test_canvas_round_trip_preserves_template_version_and_recovers_feedback_flag():
    study = graph.new_graph("unidock_active_learning")
    next(e for e in study["edges"] if e.get("loop_back")).pop("loop_back")
    cleaned = graph.clean_graph(study)
    assert cleaned["template_version"] == graph.ACTIVE_LEARNING_TEMPLATE_VERSION
    assert sum(e.get("loop_back", False) for e in cleaned["edges"]) == 1
    graph.validate_graph(cleaned)


def test_old_fixed_round_size_default_migrates_to_later_round_cap():
    study = graph.new_graph("unidock_active_learning")
    study["template_version"] = 4
    docking = next(node for node in study["nodes"] if node["data"]["kind"] == "active_docking")
    docking["data"]["params"]["batch_size"] = 10_000

    migrated = graph.migrate_template_graph(study)

    assert migrated["template_version"] == 7
    migrated_docking = next(node["data"]["params"] for node in migrated["nodes"]
                            if node["data"]["kind"] == "active_docking")
    assert migrated_docking["batch_size"] == 0


def test_legacy_docking_cpu_setting_migrates_to_per_worker_resources():
    study = graph.new_graph("unidock_active_learning")
    study["template_version"] = 4
    docking = next(node for node in study["nodes"] if node["data"]["kind"] == "active_docking")
    docking["data"]["params"].update(cpu_cores=8, worker_nodes=[{"node_id": "worker", "cpu_cores": 6}])

    migrated = graph.migrate_template_graph(study)
    migrated_docking = next(node["data"]["params"] for node in migrated["nodes"]
                            if node["data"]["kind"] == "active_docking")

    assert "cpu_cores" not in migrated_docking
    assert migrated_docking["worker_nodes"] == [{"node_id": "worker", "cpu_cores": 6}]
    assert migrated_docking["worker_batch_size"] == 2_000
    graph.validate_graph(migrated)


def test_v5_default_batch_cap_migrates_to_no_cap():
    study = graph.new_graph("unidock_active_learning")
    study["template_version"] = 5
    docking = next(node for node in study["nodes"] if node["data"]["kind"] == "active_docking")
    docking["data"]["params"]["batch_size"] = 1_000_000

    migrated = graph.migrate_template_graph(study)
    migrated_docking = next(node["data"]["params"] for node in migrated["nodes"]
                            if node["data"]["kind"] == "active_docking")

    assert migrated["template_version"] == 7
    assert migrated_docking["batch_size"] == 0


def test_decision_cannot_bypass_training_prediction_and_selection():
    study = graph.new_graph("unidock_active_learning")
    decision = next(n["id"] for n in study["nodes"] if n["data"]["kind"] == "active_decision")
    next(e for e in study["edges"] if e.get("loop_back"))["source"] = decision
    with pytest.raises(ValueError, match="Only next-batch selection"):
        graph.validate_graph(study)


def test_binder_percentage_changes_logarithmically_each_round_then_holds_final_value():
    settings = ActiveLearningSettings()
    percentages = [settings.percentile(i) for i in range(1, 11)]
    assert percentages[0] == 10
    assert percentages[-1] == 0.01
    assert all(a > b for a, b in zip(percentages, percentages[1:]))
    assert settings.percentile(11) == 0.01
    assert not settings.iteration_limit_reached(1000000)


def test_next_docking_round_uses_half_of_the_actual_centroid_seed_sample():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.active_learning_controller import _next_round_selection_counts

    study = graph.new_graph("unidock_active_learning")
    selection = next(node["data"]["params"] for node in study["nodes"]
                     if node["data"]["kind"] == "active_select")
    docking = next(node["data"]["params"] for node in study["nodes"]
                   if node["data"]["kind"] == "active_docking")
    assert selection["next_round_seed_fraction_percent"] == 50
    assert selection["predicted_fraction_percent"] == 80
    assert docking["batch_size"] == 0
    assert docking["worker_batch_size"] == 2_000
    assert _next_round_selection_counts(
        initial_seed_count=10_000, remaining_compounds=100_000, max_batch_size=10_000,
    ) == {"total": 5_000, "predicted": 4_000, "random": 1_000, "target_total": 5_000}
    assert _next_round_selection_counts(
        initial_seed_count=1_000_000, remaining_compounds=1_000_000, max_batch_size=1_000_000,
    ) == {"total": 500_000, "predicted": 400_000, "random": 100_000, "target_total": 500_000}
    assert _next_round_selection_counts(
        initial_seed_count=10_000_000, remaining_compounds=10_000_000, max_batch_size=0,
    ) == {"total": 5_000_000, "predicted": 4_000_000, "random": 1_000_000, "target_total": 5_000_000}
    assert _next_round_selection_counts(
        initial_seed_count=10_000, remaining_compounds=100_000, max_batch_size=2_500,
    ) == {"total": 2_500, "predicted": 2_000, "random": 500, "target_total": 5_000}
    assert _next_round_selection_counts(
        initial_seed_count=10_000, remaining_compounds=2_000, max_batch_size=10_000,
    ) == {"total": 2_000, "predicted": 1_600, "random": 400, "target_total": 5_000}


def test_docking_assignments_split_into_supported_launch_sizes():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.active_learning_controller import _compound_batches
    from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import MAX_COMPOUNDS_PER_LAUNCH

    rows = [{"compound_id": str(index)} for index in range(12_001)]
    batches = _compound_batches(rows, 5_000)

    assert MAX_COMPOUNDS_PER_LAUNCH == 150_000
    assert [len(batch) for batch in batches] == [5_000, 5_000, 2_001]
    assert [row["compound_id"] for batch in batches for row in batch] == [str(index) for index in range(12_001)]


def test_remote_docking_request_rewrite_uses_worker_conda_python():
    import shlex
    from mn_molecular_design_studio.engine.workflow_engine.workflows.active_learning_controller import _remote_request_rewrite_command

    command = _remote_request_rewrite_command(
        "/home/compute/envs/studio/bin/python3.11", "/tmp/job/request.json", "/tmp/job",
        "/tmp/job/alignment_ligand.sdf",
    )

    assert shlex.split(command)[0] == "/home/compute/envs/studio/bin/python3.11"
    assert "target_pdb" in command
    assert "/tmp/job/alignment_ligand.sdf" in command
    assert "reference_ligand_sdf" in command


def test_active_learning_worker_batch_shares_allocate_exact_randomized_counts():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.active_learning_controller import _weighted_worker_assignments

    rows = [{"compound_id": str(index)} for index in range(30_000)]
    workers = [{"batch_share": 1}, {"batch_share": 1}, {"batch_share": 1}]
    assignments = _weighted_worker_assignments(rows, workers, range(len(rows)))
    assert [len(item) for item in assignments] == [10_000, 10_000, 10_000]
    assert [row["compound_id"] for item in assignments for row in item] == [str(index) for index in range(30_000)]

    weighted = _weighted_worker_assignments(
        rows[:20_000], [{"batch_share": 1}, {"batch_share": 3}, {"batch_share": 3}], range(20_000),
    )
    assert [len(item) for item in weighted] == [2_857, 8_572, 8_571]


def test_remote_docking_launcher_initializes_workspace_and_surfaces_scheduler_failure(monkeypatch, tmp_path):
    from mn_molecular_design_studio.cluster import active_learning_worker
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs, workflow_queue
    from mn_molecular_design_studio import workspace

    request_path = tmp_path / "request.json"
    request_path.write_text(json.dumps({
        "schema_version": "studio.active-learning-docking-shard.v1",
        "parameters": {"cpu_cores": 1, "gpu_device": "0"},
        "round_number": 1,
        "parent_run_dir": str(tmp_path),
        "owner_id": "owner",
        "project_id": "project",
        "node_id": "worker-id",
        "node_name": "worker-name",
    }))
    child_run = tmp_path / "child-run"
    child_run.mkdir()
    metadata = {"status": "queued", "scheduler_status": "failed", "scheduler_error": "run folder is outside Results"}
    calls = []
    monkeypatch.setattr(workspace, "ensure_studio_workspace", lambda: calls.append("workspace"))
    monkeypatch.setattr(workflow_queue, "enqueue_workflow_call", lambda *args, **kwargs: calls.append("enqueue") or child_run)
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(jobs, "read_json", lambda path: metadata)
    monkeypatch.setattr(jobs, "write_json", lambda *args, **kwargs: None)

    result = active_learning_worker.run(request_path)

    assert calls == ["workspace", "enqueue"]
    assert result["status"] == "failed"
    assert result["error"] == "run folder is outside Results"


def test_initial_seed_picker_limit_comes_from_centroid_manifest(tmp_path):
    from mn_molecular_design_studio.pipeline_node_ui import _centroid_manifest_row_count

    manifest = tmp_path / "centroids.json"
    manifest.write_text(json.dumps({"row_count": 1_500_000}))

    assert _centroid_manifest_row_count(manifest) == 1_500_000


def _source_manifests(tmp_path):
    spec = {"algorithm": "morgan", "bits": 1024, "radius": 2, "use_chirality": True, "use_features": True}
    screening = {"library_id": "screening-source", "schema_version": "studio.library-fingerprints.v2", "status": "completed", "fingerprint_spec": spec}
    centroids = {"library_id": "different-centroid-source", "schema_version": "studio.library-centroids.v2", "status": "completed", "fingerprint_spec": dict(spec)}
    feature_path, centroid_path = tmp_path / "features.json", tmp_path / "centroids.json"
    feature_path.write_text(json.dumps(screening))
    centroid_path.write_text(json.dumps(centroids))
    return ({"library_id": screening["library_id"], "feature_manifest": str(feature_path)},
            {"library_id": centroids["library_id"], "centroid_manifest": str(centroid_path)})


def test_initial_centroids_and_screening_library_can_have_different_sources(tmp_path):
    screening, seed = _source_manifests(tmp_path)
    features, centroids = graph.validate_active_library_sources(screening, seed)
    assert features["library_id"] != centroids["library_id"]


@pytest.mark.parametrize("setting,value", [("radius", 3), ("use_chirality", False), ("use_features", False)])
def test_independent_sources_must_still_have_matching_model_inputs(tmp_path, setting, value):
    screening, seed = _source_manifests(tmp_path)
    path = tmp_path / "centroids.json"
    manifest = json.loads(path.read_text())
    manifest["fingerprint_spec"][setting] = value
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="compatible settings"):
        graph.validate_active_library_sources(screening, seed)


def test_legacy_centroid_selection_migrates_to_independent_input():
    study = graph.new_graph("unidock_active_learning")
    study["template_version"] = 2
    library = next(n for n in study["nodes"] if n["data"]["kind"] == "active_library")
    seed = next(n for n in study["nodes"] if n["data"]["kind"] == "active_seed")
    library["data"]["params"].update(library_id="original-source", centroid_manifest="/saved/centroids.json")
    seed["data"]["params"].pop("library_id")
    seed["data"]["params"].pop("centroid_manifest")
    migrated = graph.migrate_template_graph(study)
    seed = next(n for n in migrated["nodes"] if n["data"]["kind"] == "active_seed")
    assert seed["data"]["params"]["library_id"] == "original-source"
    assert seed["data"]["params"]["centroid_manifest"] == "/saved/centroids.json"


def test_template_titles_are_generic_and_score_conventions_are_configurable():
    study = graph.new_graph("unidock_active_learning")
    assert all("THR" not in n["data"]["label"] and "annotated" not in n["data"]["label"].lower() for n in study["nodes"])
    docking = next(n for n in study["nodes"] if n["data"]["kind"] == "active_docking")
    docking["data"]["params"].update(score_field="custom_score", score_units="engine units", score_direction="Higher is better")
    plan = graph.validate_graph(study)
    params = next(n["data"]["params"] for n in plan["graph"]["nodes"] if n["data"]["kind"] == "active_docking")
    assert params["engine"] == "UniDock Pro"
    assert params["score_field"] == "custom_score" and params["score_units"] == "engine units"
    assert params["score_direction"] == "Higher is better"
    assert all("kcal" not in f["label"] for f in graph.NODE_CATALOG["active_decision"]["fields"])


def test_active_learning_graph_queues_the_round_controller(monkeypatch, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    from mn_molecular_design_studio.engine.workflow_engine.workflows import active_learning_controller
    study = graph.new_graph("unidock_active_learning")
    monkeypatch.setattr(runner, "validate_graph", lambda *_args, **_kwargs: {"graph": study})
    calls = []
    monkeypatch.setattr(active_learning_controller, "run_active_learning_pipeline", lambda **kwargs: None)
    monkeypatch.setattr(workflow_queue, "enqueue_workflow_call", lambda function, **kwargs: calls.append((function, kwargs)) or tmp_path / "queued")
    queued = runner.enqueue_graph_pipeline(study, pipeline_name="Smoke", owner_id="test", project_id="test")
    assert queued == tmp_path / "queued"
    assert calls[0][0] is active_learning_controller.run_active_learning_pipeline
    assert calls[0][1]["task_group"] == "campaign"


def test_docking_picker_keeps_stale_registered_computers_for_live_check(monkeypatch):
    from mn_molecular_design_studio import pipeline_node_ui
    from mn_molecular_design_studio.cluster import config as cluster_config
    from mn_molecular_design_studio.cluster.store import Store

    rows = [
        {"id": "organizer", "name": "pc002975", "disabled": 0, "availability": "offline",
         "snapshot": {"accept_jobs": True, "scheduler_online": True, "resources": {
             "cpu_free": 63, "free_gpus": [{"gpu_id": 0}, {"gpu_id": 1}]}}},
        {"id": "remote", "name": "4032", "disabled": 0, "availability": "unreachable",
         "snapshot": {"accept_jobs": True, "scheduler_online": True, "resources": {
             "cpu_free": 23, "free_gpus": [{"gpu_id": 0}]}}},
        {"id": "disabled", "name": "disabled", "disabled": 1, "availability": "online",
         "snapshot": {"accept_jobs": True, "scheduler_online": True, "resources": {
             "cpu_free": 23, "free_gpus": [{"gpu_id": 0}]}}},
    ]
    monkeypatch.setattr(Store, "nodes", lambda _self: rows)
    monkeypatch.setattr(cluster_config, "identity", lambda: "organizer")

    machines = pipeline_node_ui._available_docking_machines(gpu_engine=True)

    assert [machine["id"] for machine in machines] == ["organizer", "remote"]
    assert all(machine["heartbeat_stale"] for machine in machines)


def test_active_docking_uses_live_worker_gpu_resources_after_stale_heartbeat(monkeypatch):
    from mn_molecular_design_studio.cluster import config as cluster_config, distributed_library
    from mn_molecular_design_studio.engine.workflow_engine.workflows import active_learning_controller

    monkeypatch.setattr(cluster_config, "identity", lambda: "organizer")
    monkeypatch.setattr(distributed_library, "_node_specs", lambda _nodes, **_kwargs: [
        {"node_id": "organizer", "name": "pc002975", "local": True, "cpu_cores": 32,
         "available_cpu": 63, "free_gpus": [{"gpu_id": 0}, {"gpu_id": 1}]},
        {"node_id": "remote", "name": "pc004032", "local": False, "host": "pc004032",
         "cpu_cores": 23, "available_cpu": 23, "free_gpus": [{"gpu_id": 0}]},
    ])
    monkeypatch.setattr(distributed_library, "_remote_environment",
                        lambda *_args, **_kwargs: ("/env/bin/python", "/repo", "/state", 23))

    workers = active_learning_controller._selected_worker_specs({
        "ssh_user": "compute",
        "worker_nodes": [
            {"node_id": "organizer", "gpu_device": "1", "cpu_cores": 32, "batch_share": 1},
            {"node_id": "remote", "gpu_device": "0", "cpu_cores": 23, "batch_share": 1},
        ],
    })

    assert [worker["batch_share"] for worker in workers] == [1, 1]
    assert [worker["gpu_device"] for worker in workers] == ["1", "0"]
    assert workers[1]["remote_python"] == "/env/bin/python"


def test_source_picker_accepts_libraries_with_only_the_needed_preparation(tmp_path, monkeypatch):
    from mn_molecular_design_studio import molecular_database, pipeline_graph_ui
    feature, centroid = tmp_path / "fp.json", tmp_path / "centroids.json"
    feature.touch()
    centroid.touch()
    records = [
        {"id": "screen", "preparations": [{"schema_version": "studio.library-fingerprints.v2", "path": str(feature)}]},
        {"id": "seed", "preparations": [{"schema_version": "studio.library-centroids.v2", "path": str(centroid)}]},
    ]
    class FakeDatabase:
        def external_library_records(self):
            return records
    monkeypatch.setattr(molecular_database, "MolecularDatabase", FakeDatabase)
    choices = pipeline_graph_ui._active_library_options()
    assert [row[0]["id"] for row in choices] == ["screen", "seed"]
    assert choices[0][1] and not choices[0][2]
    assert choices[1][2] and not choices[1][1]


def test_centroid_picker_exposes_all_preparations_under_the_source_database(tmp_path, monkeypatch):
    from mn_molecular_design_studio import molecular_database, pipeline_graph_ui

    source = tmp_path / "SpiroSpace_2025_client.csv"
    feature_path = tmp_path / "fingerprints" / "manifest.json"
    feature_path.parent.mkdir()
    feature_path.touch()
    centroid_specs = [
        (10_000, 781, "exact_greedy_molecular_weight_order", "2026-10-23T16:40:22+00:00"),
        (10_000, 509, "random_balanced_sampled_best_first", "2026-10-23T17:49:22+00:00"),
        (100_000, 6_357, "random_balanced_sampled_best_first", "2026-10-23T17:59:23+00:00"),
        (1_000_000, 177_430, "random_balanced_sampled_best_first", "2026-10-23T18:22:55+00:00"),
    ]
    centroid_items = []
    for index, (inputs, centers, algorithm, created_at) in enumerate(centroid_specs):
        path = tmp_path / f"centroids-{index}" / "manifest.json"
        path.parent.mkdir()
        path.write_text(json.dumps({
            "schema_version": "studio.library-centroids.v2",
            "input_rows": inputs,
            "row_count": centers,
            "algorithm": algorithm,
            "similarity_threshold": 0.6,
            "created_at": created_at,
        }))
        centroid_items.append({
            "schema_version": "studio.library-centroids.v2",
            "path": str(path),
            "rows": centers,
        })
    record = {
        "id": "spirospace",
        "name": "SpiroSpace 2025 client · first 10,000 rows",
        "sources": [{"path": str(source)}],
        "preparations": [
            {"schema_version": "studio.library-fingerprints.v2", "path": str(feature_path), "rows": 1_000_000},
            *centroid_items,
        ],
    }

    class FakeDatabase:
        def external_library_records(self):
            return [record]

    monkeypatch.setattr(molecular_database, "MolecularDatabase", FakeDatabase)
    libraries = pipeline_graph_ui._active_library_options()
    assert len(libraries) == 1
    selected_record, features, centroids = libraries[0]
    assert selected_record["id"] == "spirospace"
    assert len(centroids) == 4
    assert pipeline_graph_ui._external_library_label(selected_record, features) == (
        "SpiroSpace_2025_client.csv · 1,000,000 compounds prepared"
    )
    labels = [pipeline_graph_ui._centroid_preparation_label(item) for item in centroids]
    assert "10,000 input → 781 centroids · exact greedy · Tanimoto 0.60 · 2026-10-23" in labels
    assert "10,000 input → 509 centroids · chunked OpenVS-style, sampled · Tanimoto 0.60 · 2026-10-23" in labels
    assert "100,000 input → 6,357 centroids · chunked OpenVS-style, sampled · Tanimoto 0.60 · 2026-10-23" in labels
    assert "1,000,000 input → 177,430 centroids · chunked OpenVS-style, sampled · Tanimoto 0.60 · 2026-10-23" in labels


def test_first_database_receptor_selection_prefills_saved_ligand_centered_box(monkeypatch):
    from mn_molecular_design_studio import pipeline_graph_ui

    box = {
        "center": {"x": 27.121873292170918, "y": 37.707897689855635, "z": 32.25059854627194},
        "size": {"x": 30, "y": 20, "z": 20},
        "units": "angstrom",
        "align_longest_axis_to_x": True,
        "alignment_ligand_sdf": "/database/reference-ligand.sdf",
    }
    option = ("Database · THRα 4LNW", "/database/receptor.pdb", "THRα 4LNW", box)
    monkeypatch.setattr(pipeline_graph_ui, "_active_target_options", lambda: [option])
    monkeypatch.setattr(pipeline_graph_ui.st, "selectbox", lambda *_args, **_kwargs: option[1])
    monkeypatch.setattr(pipeline_graph_ui.st, "caption", lambda *_args, **_kwargs: None)

    selected = pipeline_graph_ui._render_active_target({}, "test-target")

    assert selected["target_pdb"] == option[1]
    assert [selected[f"center_{axis}"] for axis in "xyz"] == [box["center"][axis] for axis in "xyz"]
    assert [selected[f"size_{axis}"] for axis in "xyz"] == [30, 20, 20]
    assert selected["align_longest_axis_to_x"] is True
    assert selected["alignment_ligand_sdf"] == box["alignment_ligand_sdf"]


def test_active_learning_unidock_parameters_include_reference_ligand_alignment(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.active_learning_controller import _docking_params

    ligand = tmp_path / "reference.sdf"
    ligand.write_text("reference ligand")
    target = {
        "target_pdb": "/db/receptor.pdb",
        "center_x": 26.6, "center_y": 38.2, "center_z": 32.7,
        "size_x": 30, "size_y": 20, "size_z": 20,
        "align_longest_axis_to_x": True,
        "alignment_ligand_sdf": str(ligand),
    }
    params = {"engine": "UniDock Pro", "score_field": "vina_score", "score_direction": "Lower is better"}
    worker = {"node_name": "RTX 5090", "cpu_cores": 16, "gpu_device": "0"}

    request = _docking_params(params, target, tmp_path / "batch.json", round_number=1, worker=worker, limit=2)

    assert request["align_longest_axis_to_x"] is True
    assert request["alignment_ligand_sdf"] == str(ligand)
    assert request["preserve_box_size_after_alignment"] is True
    assert [request[f"size_{axis}"] for axis in "xyz"] == [30, 20, 20]
    assert request["mode"] == "hybrid"
    assert request["search_mode"] == "balance"
    assert request["reference_ligand_sdf"] == str(ligand)
