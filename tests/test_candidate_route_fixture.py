"""Deterministic, CPU-only exercise of heterogeneous graph candidate handoffs."""
from __future__ import annotations

import copy
from pathlib import Path
from types import SimpleNamespace

from mn_molecular_design_studio import candidate_sets, pipeline_graph
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner


def _node(kind, index, **updates):
    spec = pipeline_graph.NODE_CATALOG[kind]
    params = {field["name"]: copy.deepcopy(field["default"]) for field in spec["fields"]}
    params.update(updates)
    return {"id": kind, "type": "studio", "position": {"x": 260 * index, "y": 0},
            "data": {"kind": kind, "label": spec["label"], "params": params}}


def _edge(index, source, target):
    return {"id": f"route-{index}", "source": source, "sourceHandle": "out",
            "target": target, "targetHandle": "in"}


def test_backbone_to_redesign_prediction_filter_analysis_biotite_results(
        tmp_path, monkeypatch, design_source_factory):
    source, _ = design_source_factory(tmp_path / "runs/design/backbone-fixture")
    structure = source / "designed.cif"
    row = {
        "candidate_id": "backbone-001", "source_tool": "fixture-backbone", "stage": "generation.backbone",
        "complex_pdb": str(structure), "target_pdb": str(structure), "binder_pdb": str(structure),
        "target_chains": ["A"], "binder_chains": ["B"], "binder_sequence": "ACDEFG",
        "parents": ["target-input"], "metrics": {}, "raw_metadata": {"retry": 1, "design_method": "fixture"},
    }
    write_candidates(source, "fixture-backbone", [row])
    write_json(source / "metadata.json", {"status": "completed", "task_group": "generation"})
    write_json(source / "result.json", {"success": True})
    write_json(source / "molecular_design_studio.json", {"project_id": "fixture-project"})

    nodes = [
        _node("backbone_source", 0, source_run_dir=str(source)),
        _node("foundry_mpnn", 1),
        _node("boltz2", 2),
        _node("candidate_merge", 3),
        _node("shared_filter", 4, binder_format="miniprotein", target_type="structured",
              policy_profile="custom", rules=[{"field": "iptm", "operator": ">=", "value": 0.7}]),
        _node("complex_analysis", 5, entry_mode="accepted_only"),
        _node("biotite_interface", 6),
        _node("results", 7),
    ]
    graph = pipeline_graph.clean_graph({
        "schema_version": pipeline_graph.SCHEMA, "id": "a" * 32, "name": "Fixture route",
        "nodes": nodes, "edges": [_edge(i, nodes[i]["id"], nodes[i+1]["id"]) for i in range(len(nodes)-1)],
        "viewport": {"x": 0, "y": 0, "zoom": 0.7},
    })
    parent = tmp_path / "runs/campaign/fixture-route"
    parent.mkdir(parents=True)
    write_json(parent / "metadata.json", {"status": "queued", "job_type": "visual_pipeline"})
    write_json(parent / "input.json", {"params": {"pipeline_name": "Fixture route"}})
    monkeypatch.setattr(runner, "validate_graph", pipeline_graph.validate_graph)
    monkeypatch.setattr(runner, "create_job", lambda *a, **k: SimpleNamespace(run_dir=parent))
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *a, **k: None)
    # Native GPU/container preflight is intentionally skipped: this test checks
    # the graph's handoff contract using deterministic job receipts.
    monkeypatch.setattr("mn_molecular_design_studio.candidate_pipeline.preflight", lambda *a: None)
    monkeypatch.setattr(runner, "source_generations", lambda *a, **k: (source, [copy.deepcopy(row)]))

    def submit(request, parent_dir, node, **kwargs):
        kind = node["data"]["kind"]
        params = request["parameters"]
        candidate_path = params.get("candidates_jsonl")
        if candidate_path is None and params.get("tasks"):
            candidate_path = (Path(params["tasks"][0]["source_run_dir"])
                              / "artifacts/normalized_candidates/candidates.jsonl")
        inputs = read_candidates(Path(candidate_path))
        if kind == "foundry_mpnn":
            assert len(inputs) == 1 and inputs[0]["candidate_id"] == "backbone-001"
            assert inputs[0]["stage"] == "generation.backbone" and inputs[0]["binder_sequence"] == "ACDEFG"
            assert inputs[0]["complex_pdb"] == str(structure)
        elif kind == "boltz2":
            assert len(inputs) == 1 and inputs[0]["candidate_id"] == "sequence-001"
            assert inputs[0]["stage"] == "sequence_design" and inputs[0]["binder_sequence"] == "KLMNPQ"
            assert inputs[0]["parents"] == ["backbone-001"]
        elif kind == "biotite_interface":
            assert len(inputs) == 1 and inputs[0]["stage"] == "complex_refolding"
            assert inputs[0]["binder_sequence"] == "KLMNPQ" and inputs[0]["metrics"]["iptm"] == 0.82
            assert inputs[0]["parents"] == ["sequence-001"]
            assert Path(inputs[0]["complex_pdb"]).is_file()
        child = parent_dir / "mock-children" / kind
        child.mkdir(parents=True)
        outputs = []
        for source_row in inputs:
            if kind == "foundry_mpnn":
                item = copy.deepcopy(source_row)
                item.update(candidate_id="sequence-001", source_tool="foundry_mpnn", stage="sequence_design",
                            parents=[source_row["candidate_id"]], binder_sequence="KLMNPQ")
            elif kind == "boltz2":
                item = copy.deepcopy(source_row)
                item.update(candidate_id="prediction-model-0", source_tool="boltz2", stage="complex_refolding",
                            parents=[source_row["candidate_id"]], metrics={"iptm": 0.82})
                item["raw_metadata"]["model_index"] = 0
            elif kind == "biotite_interface":
                item = copy.deepcopy(source_row)
                item.update(candidate_id=source_row["candidate_id"] + "-biotite", source_tool="biotite_interface",
                            stage="analysis", parents=[source_row["candidate_id"]],
                            metrics={**(source_row.get("metrics") or {}), "qc_interface_contact_pairs": 12})
            else:
                raise AssertionError(f"Unexpected execution node: {kind}")
            outputs.append(item)
        write_candidates(child, kind, outputs)
        write_json(child / "metadata.json", {"status": "completed", "tool": kind})
        write_json(child / "result.json", {"success": True, "outputs": {"candidates": outputs}})
        return child

    monkeypatch.setattr(runner, "_submit", submit)
    output_dir = runner.run_graph_pipeline(graph, "Fixture route", "fixture-owner", "fixture-project",
                                           reuse_completed=False, gpu_device="none")
    result = read_json(output_dir / "result.json")
    assert result["success"] is True
    assert len(result["outputs"]["candidates"]) == 1
    final = result["outputs"]["candidates"][0]
    assert final["candidate_id"].endswith("-biotite")
    assert final["binder_sequence"] == "KLMNPQ"
    assert final["parents"]
    assert final["raw_metadata"]["complex_analysis_entry"]["native_candidate_id"] == "prediction-model-0"
    assert final["metrics"]["iptm"] == 0.82
    state = read_json(output_dir / "graph_execution.json")["nodes"]
    assert [state[k]["status"] for k in ("backbone_source", "foundry_mpnn", "boltz2", "candidate_merge",
                                            "shared_filter", "complex_analysis", "biotite_interface", "results")] == ["completed"] * 8
    assert read_candidates(source / "artifacts/normalized_candidates/candidates.jsonl")[0]["candidate_id"] == "backbone-001"
