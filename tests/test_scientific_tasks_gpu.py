"""Replay retained real generators through newly configured sequence/prediction nodes."""
import json
import os
from pathlib import Path

import pytest

from mn_molecular_design_studio import design_validation as dv, mpnn_input
from mn_molecular_design_studio.scientific_tasks import new_task_graph
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as graph, mpnn, sequence_design, optimized_prediction


@pytest.mark.gpu_e2e
@pytest.mark.parametrize("method", ["rfdiffusion", "rfd4"])
def test_real_task_recipe_sequence_variants_prediction_and_results(gpu_runtime, method):
    runtime = gpu_runtime
    retained = os.environ.get("MDS_TASK_GENERATION_" + method.upper())
    if not retained:
        pytest.skip("Set MDS_TASK_GENERATION_" + method.upper() + " to a retained successful real generator run")
    source = Path(retained).resolve()
    assert read_json(source / "result.json")["success"] is True
    runtime["require_image"](mpnn_input.IMAGE)
    runtime["require_image"](optimized_prediction.ENGINE_PROFILES["boltz2"]["image"])
    for model in ("protein_mpnn", "soluble_mpnn"):
        if not mpnn_input.checkpoint(model).is_file():
            pytest.skip("Missing installed " + model + " weights")
    for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "ccd.pkl"):
        if not (runtime["reference"] / "boltz_models" / name).is_file():
            pytest.skip("Missing installed Boltz prerequisite: " + name)
    generated = read_candidates(source)
    recipe = new_task_graph("targeted_protein_binder", {"methods": [method], "soluble_mpnn": True,
        "soluble_settings": {"sequences": 1, "seed": 20261027, "sampling_temp": .4},
        "prediction_settings": {"boltz2": {"samples": 1, "num_recycles": 1, "sampling_steps": 20,
            "seed": 20261027, "target_use_msa": False, "use_target_template": False, "generate_ipsae": False}},
        "selection": {"required_metrics": ["iptm"], "ranking_metric": "iptm", "limit": 2}})
    output = {"run_dir": str(source), "candidates_jsonl": str(source / "artifacts/normalized_candidates/candidates.jsonl")}
    stages = []
    for node in recipe["nodes"]:
        kind = node["data"]["kind"]
        if kind not in {"ligandmpnn", "boltz2"}:
            continue
        request = graph.node_request(kind, node["data"]["params"], output, "phase5-task-" + method)
        parameters = request["parameters"]
        parameters["gpu_device"] = runtime["device"]
        if kind == "ligandmpnn":
            if parameters["sequence_stage"] == "initial":
                parameters.update(seed=20261027, backend="shared_ligandmpnn")
            run = mpnn.run_mpnn_sequence_design(**parameters)
        else:
            assert parameters["target_use_msa"] is False
            assert parameters["engine_options"]["use_msa"] is False
            run = optimized_prediction.run_optimized_prediction(**parameters)
        assert read_json(run / "result.json")["success"] is True, run
        output = {"run_dir": str(run), "candidates_jsonl": str(run / "artifacts/normalized_candidates/candidates.jsonl")}
        stages.append({"node": node["id"], "kind": kind, "run_dir": str(run), "parameters": parameters})
        if kind == "ligandmpnn" and parameters["sequence_stage"] == "soluble_optimization":
            variants = read_candidates(run)
            assert len(variants) == 2
            parent, child = variants
            if method == "rfd4":
                assert parent["binder_sequence"] == generated[0]["binder_sequence"]
            assert child["parents"] == [parent["candidate_id"]]
            path = Path(parent["complex_pdb"])
            structure = dv.protein_chains(dv.read_structure(path))
            contacts = sequence_design._interface_contacts_for_structure(path, parent["binder_chains"], parent["target_chains"], cutoff=4.)
            labels = [f"{chain}{res.id[1]}{res.id[2].strip()}" for chain in parent["binder_chains"] for res in structure[chain]]
            assert contacts
            assert all(a == b for label, a, b in zip(labels, parent["binder_sequence"], child["binder_sequence"]) if label in contacts or a == "C")
    predicted = read_candidates(Path(output["run_dir"]))
    assert {r["binder_sequence"] for r in predicted} == {r["binder_sequence"] for r in variants}
    selection_node = next(n for n in recipe["nodes"] if n["data"]["kind"] == "filter")
    selected = graph.filter_candidates(predicted, selection_node["data"]["params"])
    assert len(selected) == 2
    results = graph._local_source(runtime["root"], "task-results", selected, final_evidence=True)
    assert len(read_candidates(Path(results["run_dir"]))) == 2
    (runtime["root"] / "task-recipe-evidence.json").write_text(json.dumps({"method": method,
        "generation_reused": str(source), "recipe": recipe, "executed_stages": stages, "results": results}, indent=2))
