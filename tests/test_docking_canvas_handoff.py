import json
from pathlib import Path

import pytest


def _docking_source(root, *, external_pose=False):
    from mn_molecular_design_studio import ligand_compounds

    source = root / "ligand-docking" / "dock-001"
    inputs = source / "artifacts/input"
    poses = source / "artifacts/poses"
    inputs.mkdir(parents=True)
    poses.mkdir(parents=True)
    receptor = inputs / "receptor.pdb"
    receptor.write_text("HEADER receptor fixture\nEND\n")
    compound_set = inputs / "compound_set.json"
    ligand_compounds.write_compound_set(compound_set, name="Fixture ligands", compounds=[
        {"compound_id": "CMP-A", "compound_name": "Ligand A", "smiles": "C[C@H](O)F",
         "properties": {"protonation_state": "as_supplied", "atom_mapping": {"C1": 1}}},
    ], owner_id="owner", project_id="project", created_at="2026-11-06T00:00:00Z")
    first, second = poses / "cmp-a-pose-1.pdbqt", poses / "cmp-a-pose-2.pdbqt"
    first.write_text("REMARK pose one\nEND\n")
    second.write_text("REMARK pose two\nEND\n")
    external = root / "unrelated.pdbqt"
    external.write_text("REMARK unrelated external pose\nEND\n")
    table = source / "artifacts/docking_scores.csv"
    second_path = str(external) if external_pose else "artifacts/poses/cmp-a-pose-2.pdbqt"
    table.write_text(
        "compound_id,engine,status,pose_rank,best_score_kcal_mol,pose_path,receptor_model_id\n"
        "CMP-A,vina,completed,1,-7.1,artifacts/poses/cmp-a-pose-1.pdbqt,\n"
        f"CMP-A,vina,completed,2,-8.0,{second_path},\n"
    )
    (source / "metadata.json").write_text(json.dumps({
        "run_id": "dock-001", "task_group": "ligand-docking", "tool": "vina",
        "engine": "AutoDock Vina", "workflow": "ligand.docking.vina",
        "studio_project_id": "project", "status": "completed",
    }))
    (source / "input.json").write_text(json.dumps({"params": {"center_x": 1.0, "use_scrub": True}}))
    (source / "result.json").write_text(json.dumps({"success": True}))
    (source / "artifact_contract.json").write_text(json.dumps({
        "schema_version": 1, "output_type": "ligand_docking_result/v1",
        "outputs": [{"name": "docking_scores", "score_units": "kcal/mol"}],
    }))
    return source


def test_selected_docking_pose_becomes_exact_candidate_and_plip_request(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.engine.workflow_engine.workflows import quality
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import (
        docking_candidates, docking_pose_options,
    )
    from mn_molecular_design_studio import run_graph

    source = _docking_source(tmp_path)
    monkeypatch.setattr(config, "identity", lambda: "node-test")
    source_ref = "node-test:ligand-docking/dock-001"
    options = docking_pose_options(source, source_job_reference=source_ref, source_node_id="node-test")
    assert len(options) == 2 and options[0]["pose_id"] != options[1]["pose_id"]
    chosen_id = options[1]["pose_id"]
    candidate = docking_candidates(source, [chosen_id], source_job_reference=source_ref,
                                   source_node_id="node-test", source_workstation="pc002975")[0]
    assert candidate["candidate_id"] == chosen_id
    assert candidate["parents"] == ["CMP-A"]
    assert candidate["binder_sequence"] is None
    assert candidate["complex_pdb"].endswith("cmp-a-pose-2.pdbqt")
    assert candidate["target_pdb"].endswith("receptor.pdb")
    assert candidate["metrics"]["best_score_kcal_mol"] == -8.0
    raw = candidate["raw_metadata"]
    assert raw["artifact_schema"] == "ligand_docking_result/v1"
    assert raw["source_job_reference"] == source_ref
    assert raw["source_workstation"] == "pc002975"
    assert raw["source_artifact_sha256"]["pose"]
    assert raw["source_artifact_sha256"]["receptor"]
    assert raw["source_artifact_references"]["pose"].endswith("artifacts/poses/cmp-a-pose-2.pdbqt")
    assert raw["ligand_identity"]["smiles"] == "C[C@H](O)F"
    assert raw["ligand_identity"]["properties"]["protonation_state"] == "as_supplied"
    assert raw["docking_engine"] == "vina"
    assert raw["docking_score_units"] == "kcal/mol"

    row = {"node_id": "node-test", "node_name": "pc002975", "task_group": "ligand-docking",
           "run_id": "dock-001", "run_dir": str(source), "status": "completed",
           "has_docking_results": True, "selected_pose_ids": [chosen_id], "source_job_reference": source_ref}
    deferred = run_graph.deferred_source(row, {})
    request = run_graph.source_request("plip", {}, deferred, "PLIP selected poses")
    assert request["workflow"] == "analysis.plip"
    assert request["parameters"]["candidate_ids"] == [chosen_id]
    assert request["parameters"]["candidates_jsonl"] is None
    assert request["parameters"]["source_job_reference"] == source_ref
    assert request["parameters"]["source_workstation_id"] == "node-test"

    captured = {}
    monkeypatch.setattr(quality, "run_quality", lambda *args, **kwargs: captured.update(kwargs) or source)
    quality.run_plip(**request["parameters"])
    assert captured["candidate_ids"] == [chosen_id]
    assert captured["source_job_reference"] == source_ref


def test_canvas_source_requires_selected_authorized_pose_and_rejects_unrelated_pose(tmp_path, monkeypatch):
    from mn_molecular_design_studio import logical_runs
    from mn_molecular_design_studio.cluster import candidate_batches, config
    from mn_molecular_design_studio import runtime
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_pose_options
    from mn_molecular_design_studio import run_graph

    source = _docking_source(tmp_path)
    monkeypatch.setattr(config, "identity", lambda: "node-test")
    monkeypatch.setattr(runtime, "runs_root", lambda: tmp_path)
    monkeypatch.setattr(candidate_batches, "actor", lambda *args, **kwargs: {"id": "owner", "role": "user"})
    row = {"node_id": "node-test", "node_name": "pc002975", "task_group": "ligand-docking",
           "run_id": "dock-001", "run_dir": str(source), "status": "completed", "has_docking_results": True,
           "project_id": "project"}
    source_ref = logical_runs.reference(row)
    options = docking_pose_options(source, source_job_reference=source_ref, source_node_id="node-test")
    logical_run = {"id": "run:123", "outputs": [{"id": "results", "reference": "run:123~results",
        "sources": [row]}]}
    monkeypatch.setattr(logical_runs, "catalog", lambda *args, **kwargs: [logical_run])
    selected = [{"source_job_reference": source_ref, "pose_id": options[0]["pose_id"]}]
    loaded = run_graph.load_source("run:123~results", owner_id="owner", project_id="project",
                                  selected_docking_poses=selected)
    assert len(loaded["run_sources"]) == 1
    assert loaded["run_sources"][0]["selected_pose_ids"] == [options[0]["pose_id"]]

    with pytest.raises(ValueError, match="not part of this authorized Run output"):
        run_graph.load_source("run:123~results", owner_id="owner", project_id="project",
            selected_docking_poses=[{"source_job_reference": "other-node:ligand-docking/elsewhere",
                                     "pose_id": options[0]["pose_id"]}])

    external_source = _docking_source(tmp_path / "external", external_pose=True)
    with pytest.raises(ValueError, match="outside its source job"):
        docking_pose_options(external_source, source_job_reference=source_ref, source_node_id="node-test")


def test_docking_pose_without_recorded_chemistry_is_not_selectable(tmp_path):
    from mn_molecular_design_studio import ligand_compounds
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_candidates, docking_pose_options

    source = _docking_source(tmp_path)
    (source / "artifacts/input/compound_set.json").unlink()
    ligand_compounds.write_compound_set(source / "artifacts/input/compound_set.json", name="Incomplete", compounds=[
        {"compound_id": "OTHER", "smiles": "CCO"},
    ], owner_id="owner", project_id="project", created_at="now")
    options = docking_pose_options(source)
    assert options and all(not item["available"] for item in options)
    with pytest.raises(ValueError, match="no chemistry"):
        docking_candidates(source, [options[0]["pose_id"]])


def test_missing_or_oversized_pose_is_retained_as_score_but_not_executable(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import quality_ligands
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import (
        docking_candidates, docking_pose_options,
    )

    source = _docking_source(tmp_path)
    (source / "artifacts/poses/cmp-a-pose-1.pdbqt").unlink()
    missing = docking_pose_options(source)
    assert missing[0]["available"] is False
    assert "Missing source artifact" in missing[0]["unavailable_reason"]
    assert missing[0]["score_row"]["best_score_kcal_mol"] == "-7.1"
    with pytest.raises(ValueError, match="Missing source artifact"):
        docking_candidates(source, [missing[0]["pose_id"]])

    pose = source / "artifacts/poses/cmp-a-pose-2.pdbqt"
    pose.write_text("REMARK " + "large-pose " * 20 + "\nEND\n")
    monkeypatch.setattr(quality_ligands, "MAX_DOCKING_ARTIFACT_BYTES", 64)
    oversized = docking_pose_options(source)
    assert oversized[1]["available"] is False
    assert "Oversized source artifact" in oversized[1]["unavailable_reason"]


def test_remote_pose_absent_from_complete_manifest_is_not_selectable():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_pose_options
    from mn_molecular_design_studio.ligand_compounds import normalize_compounds

    payload = {"record_type": "compound_set", "compounds": normalize_compounds([
        {"compound_id": "CMP-A", "compound_name": "Ligand A", "smiles": "C[C@H](O)F"},
    ])}
    options = docking_pose_options(
        scores_text=("compound_id,status,score,pose_path\n"
                     "CMP-A,completed,-9.1,artifacts/poses/pose-1.pdbqt\n"),
        compound_set_payload=payload,
        artifact_paths=[{"relative_path": "ligand-docking/remote-1/artifacts/input/receptor.pdb", "size_bytes": 128}],
        artifact_manifest_complete=True,
    )

    assert len(options) == 1
    assert options[0]["available"] is False
    assert "Missing source artifact" in options[0]["unavailable_reason"]
    assert options[0]["score_row"]["score"] == "-9.1"


def test_remote_results_selector_uses_authorized_job_artifact_fetch_without_manifest(monkeypatch):
    import streamlit as st
    from mn_molecular_design_studio import run_inputs_ui
    from mn_molecular_design_studio.cluster import config
    compound_path = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/vina/imatinib.json"
    ligand = json.loads(compound_path.read_text())[0]
    from mn_molecular_design_studio.ligand_compounds import normalize_compounds
    payload = {"record_type": "compound_set", "compounds": normalize_compounds([
        {"compound_id": "STI-571", "compound_name": "Imatinib", "smiles": ligand["smiles"]}
    ])}
    score_text = ("compound_id,engine,status,pose_rank,best_score_kcal_mol,pose_path\n"
                  "STI-571,unidock_pro,completed,1,-8.2,artifacts/poses/STI-571-1.pdbqt\n")
    fetched = []
    def fetch(job_id, path, token, state_key, label="saved structure"):
        fetched.append((job_id, path))
        return score_text if path.endswith("docking_scores.csv") else json.dumps(payload)

    monkeypatch.setattr(config, "identity", lambda: "local-node")
    monkeypatch.setattr(st, "session_state", {"_molecular_database_user": {"id": "authorized-user"},
                                               "studio-session-token": "session-token"})
    monkeypatch.setattr("mn_molecular_design_studio.auth.is_development_mode", lambda: False)
    monkeypatch.setattr("mn_molecular_design_studio.cluster.ui._fetch_remote_artifact_text", fetch)
    row = {"node_id": "remote-node", "node_name": "PC004033", "task_group": "ligand-docking",
           "run_id": "remote-run", "cluster_job_id": "authorized-job", "artifact_manifest_complete": True,
           "artifacts": [
               {"relative_path": "ligand-docking/remote-run/artifacts/input/receptor.pdb", "size_bytes": 128},
               {"relative_path": "ligand-docking/remote-run/artifacts/poses/STI-571-1.pdbqt", "size_bytes": 1024},
           ]}

    options = run_inputs_ui._docking_pose_options(row)

    assert len(options) == 1 and options[0]["available"]
    assert options[0]["source_node_id"] == "remote-node"
    assert options[0]["pose_relative_path"] == "artifacts/poses/STI-571-1.pdbqt"
    assert fetched == [("authorized-job", "artifacts/docking_scores.csv"),
                       ("authorized-job", "artifacts/input/compound_set.json")]


def test_run_results_canvas_editor_persists_selected_run_and_poses(monkeypatch):
    from mn_molecular_design_studio import pipeline_graph_ui, run_inputs_ui

    selected = {"id": "source", "data": {"kind": "run_source", "params": {}}}
    choice = {
        "reference": "run:0123456789abcdef01234567~results",
        "selected_docking_poses": [{"source_job_reference": "node:ligand-docking/run",
                                    "pose_id": "pose-001"}],
    }
    received = {}

    def select(project_id, key, **kwargs):
        received.update(project_id=project_id, key=key, **kwargs)
        return choice

    monkeypatch.setattr(run_inputs_ui, "select", select)
    contexts = pipeline_graph_ui._render_node_input(
        selected, owner_id="owner", project_id="project", search_jobs=None, prefix="editor",
    )

    assert received["project_id"] == "project"
    assert received["allow_docking"] is True
    assert contexts == {"source": {
        "run_reference": choice["reference"],
        "selected_docking_poses": choice["selected_docking_poses"],
    }}
    assert selected["data"]["params"]["run_reference"] == choice["reference"]
    assert selected["data"]["params"]["selected_docking_poses"] == choice["selected_docking_poses"]


def test_raw_docking_pose_connections_are_limited_to_ligand_analysis():
    from mn_molecular_design_studio import pipeline_graph

    def node(kind, identifier):
        spec = pipeline_graph.NODE_CATALOG[kind]
        return {"id": identifier, "type": "studio", "position": {"x": 0, "y": 0},
                "data": {"kind": kind, "label": spec["label"],
                         "params": {field["name"]: field["default"] for field in spec["fields"]}}}

    source = node("run_source", "source")
    source["data"]["params"].update(run_reference="run:0123456789abcdef01234567~results",
                                    selected_docking_poses=[{"source_job_reference": "node:job/run",
                                                             "pose_id": "dockpose-example"}])
    plip, protein_analysis, result = node("plip", "plip"), node("biotite_interface", "biotite"), node("results", "results")
    graph = {"schema_version": pipeline_graph.SCHEMA, "id": "a" * 32,
             "name": "Docking to PLIP", "nodes": [source, plip, result],
             "edges": [{"id": "source-plip", "source": "source", "target": "plip",
                        "sourceHandle": "out", "targetHandle": "in"},
                       {"id": "plip-results", "source": "plip", "target": "results",
                        "sourceHandle": "out", "targetHandle": "in"}]}

    pipeline_graph.validate_graph(graph)
    graph["nodes"] = [source, protein_analysis, result]
    graph["edges"][0].update(target="biotite")
    graph["edges"][1].update(source="biotite")
    with pytest.raises(ValueError, match="protein–ligand analysis"):
        pipeline_graph.validate_graph(graph)


def test_remote_artifact_failure_categories_are_distinct():
    from mn_molecular_design_studio.cluster import ui as cluster_ui

    assert cluster_ui._artifact_failure_detail("The requested source artifact is missing") .startswith("Missing source artifact:")
    assert cluster_ui._artifact_failure_detail("The requested source artifact exceeds the 512 MB download limit.").startswith("Oversized source artifact:")
    assert cluster_ui._artifact_failure_detail("Job not found").startswith("Source job unavailable or you are not authorized")
    assert cluster_ui._artifact_failure_detail("invalid path outside its source job").startswith("Rejected artifact path:")
    assert cluster_ui._artifact_failure_detail("socket timed out", connection_failure=True).startswith("Retrieval/connection failure:")


def test_worker_artifact_errors_distinguish_missing_unsafe_and_oversized(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import agent as agent_module, bundles

    run = tmp_path / "ligand-docking" / "dock-1"
    (run / "artifacts").mkdir(parents=True)
    outside = tmp_path / "outside.pdbqt"
    outside.write_text("x")
    monkeypatch.setattr(
        "mn_molecular_design_studio.engine.workflow_engine.services.local_automation.resolve_job",
        lambda _reference: {"run_dir": str(run)},
    )
    worker = object.__new__(agent_module.Agent)

    with pytest.raises(FileNotFoundError, match="Missing source artifact"):
        worker.artifact({"payload": {"local_reference": "ligand-docking/dock-1", "path": "artifacts/missing.pdbqt"}})
    with pytest.raises(ValueError, match="outside its source job"):
        worker.artifact({"payload": {"local_reference": "ligand-docking/dock-1", "path": "../../outside.pdbqt"}})

    pose = run / "artifacts/pose.pdbqt"
    pose.write_text("too large")
    monkeypatch.setattr(bundles, "MAX_BYTES", 1)
    with pytest.raises(ValueError, match="exceeds the 512 MB download limit"):
        worker.artifact({"payload": {"local_reference": "ligand-docking/dock-1", "path": "artifacts/pose.pdbqt"}})
