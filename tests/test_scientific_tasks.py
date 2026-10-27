"""Task recipes must preserve scientific routing and all historical graph entry points."""
import copy
import json

import pytest

from mn_molecular_design_studio import pipeline_graph as graph
from mn_molecular_design_studio.scientific_tasks import (
    TASKS, BACKBONE, NATIVE, PREDICTION_MODES, new_task_graph, task_methods, prediction_methods,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner
from test_pipeline_graph import execution_harness
from test_api_integrations import api_data, headers


def kinds(g):
    return [n["data"]["kind"] for n in g["nodes"]]


@pytest.mark.parametrize("task", TASKS)
def test_each_task_has_a_valid_editable_default(task):
    result = graph.validate_graph(new_task_graph(task))["graph"]
    assert result["task_configuration"]["task_id"] == task
    assert len(result["nodes"]) < 12
    assert sum(kind in {"results", "active_results"} for kind in kinds(result)) == 1
    assert graph.clean_graph(result) == result


def test_catalogue_has_six_standard_and_two_experimental_tasks_and_blank():
    assert len(TASKS) == 8
    assert sum(bool(t.get("experimental")) for t in TASKS.values()) == 2
    assert len(graph.STARTER_TEMPLATES) == 9
    assert "blank" in graph.STARTER_TEMPLATES.values()


@pytest.mark.parametrize("method", task_methods("targeted_protein_binder"))
def test_binder_methods_route_backbones_codesigns_and_native_workflows(method):
    result = new_task_graph("targeted_protein_binder", {"methods": [method]})
    sequence = [n for n in result["nodes"] if n["data"]["kind"] == "ligandmpnn"]
    assert bool(sequence) == (method in BACKBONE)
    if sequence:
        assert sequence[0]["data"]["params"]["sequence_stage"] == "initial"
        assert sequence[0]["data"]["params"]["model_type"] == "protein_mpnn"
    predictions = [n for n in result["nodes"] if graph.NODE_CATALOG[n["data"]["kind"]]["category"] == "Structure prediction"]
    assert len(predictions) == (0 if method in NATIVE else 1)
    assert "rosetta_metrics" not in kinds(result)
    assert "ipsae" not in kinds(result)  # Boltz already supplies this measurement.


@pytest.mark.parametrize("method", ["rfdiffusion", "rfd4", "boltzgen", "esm_binder"])
def test_soluble_stage_is_explicit_and_cannot_disable_interface_or_parent_retention(method):
    result = new_task_graph("targeted_protein_binder", {"methods": [method], "soluble_mpnn": True,
        "soluble_settings": {"optimization_enabled": False, "fix_interface_residues": False, "interface_cutoff": 6.5}})
    stages = [n["data"]["params"] for n in result["nodes"] if n["data"]["kind"] == "ligandmpnn"]
    params = stages[-1]
    assert params["sequence_stage"] == "soluble_optimization"
    assert params["model_type"] == "soluble_mpnn"
    assert params["optimization_enabled"] is params["fix_interface_residues"] is True
    assert params["interface_cutoff"] == 6.5
    assert "boltz2" in kinds(result)  # Native parents and requested variants are both evaluated.


@pytest.mark.parametrize("method", task_methods("targeted_protein_binder", binder_format="vhh"))
def test_vhh_methods_retain_protected_format_and_annotation_before_redesign(method):
    result = new_task_graph("targeted_protein_binder", {"methods": [method], "binder_format": "vhh", "soluble_mpnn": True})
    assert kinds(result).index(method) < kinds(result).index("sabr") < kinds(result).index("ligandmpnn")
    params = next(n["data"]["params"] for n in result["nodes"] if n["data"]["kind"] == method)
    assert params.get("modality") == "VHH" or params.get("task") == "nanobody" or params.get("binder_type") == "vhh"


@pytest.mark.parametrize("mode", PREDICTION_MODES)
def test_prediction_modes_default_to_one_compatible_engine(mode):
    result = new_task_graph("molecular_structure_prediction", {"mode": mode})
    assert kinds(result) == ["cofold_input", "boltz2", "results"]
    params = result["nodes"][0]["data"]["params"]
    assert params["task_type"] == mode
    assert params["cyclic_partner"] == (mode == "protein_macrocyclic")


@pytest.mark.parametrize("mode", PREDICTION_MODES)
def test_every_prediction_method_with_template_support_keeps_engine_budget_constraints(mode):
    inputs = {"target_pdb": "target.pdb"}
    for method in prediction_methods(mode, inputs):
        result = new_task_graph("molecular_structure_prediction", {"mode": mode, "inputs": inputs, "methods": [method]})
        graph.validate_graph(result)
        if method == "af2":
            assert result["nodes"][1]["data"]["params"]["af2_validation_mode"] == "target_only"


@pytest.mark.parametrize("task", TASKS)
def test_all_advertised_methods_generate_valid_graphs(task):
    for method in task_methods(task):
        graph.validate_graph(new_task_graph(task, {"methods": [method]}))


@pytest.mark.parametrize("task,config", [
    ("targeted_protein_binder", {"methods": ["rfd3", "rfd4"]}),
    ("targeted_protein_binder", {"methods": ["protpardelle"]}),
    ("targeted_protein_binder", {"binder_format": "vhh", "methods": ["rfdiffusion"]}),
    ("targeted_protein_binder", {"methods": ["rfd4"], "method_settings": {"rfd4": {"task": "unconditional"}}}),
    ("new_protein_fold", {"soluble_mpnn": True}),
    ("ligand_conditioned_pocket", {"prediction": "af2"}),
    ("molecular_structure_prediction", {"mode": "protein_ligand", "methods": ["prediction_esmfold"]}),
    ("molecular_structure_prediction", {"mode": "protein_ligand", "inputs": {"ligand_input_type": "ccd"}, "methods": ["prediction_chai1"]}),
    ("molecular_structure_prediction", {"mode": "protein_protein", "methods": ["af2"]}),
    ("molecular_structure_prediction", {"mode": "protein_macrocyclic", "methods": ["rf3"]}),
    ("molecular_structure_prediction", {"mode": "protein", "inputs": {"dna_sequences": "AGCT"}}),
])
def test_invalid_scientific_combinations_are_rejected(task, config):
    with pytest.raises(ValueError):
        new_task_graph(task, config)


def test_prediction_compatibility_accounts_for_template_and_ligand_notation():
    assert "af2" not in prediction_methods("protein_protein")
    assert "af2" in prediction_methods("protein_protein", {"target_pdb": "target.pdb"})
    assert "prediction_chai1" not in prediction_methods("protein_ligand", {"ligand_input_type": "ccd"})
    assert "prediction_chai1" in prediction_methods("protein_ligand", {"ligand_input_type": "smiles"})


def test_specialists_and_optional_shortlist_preserve_specific_stages():
    cdr = new_task_graph("antibody_cdr_redesign")
    assert kinds(cdr) == ["project_design_source", "sabr", "abmpnn", "sequence_qc", "boltz2", "results"]
    assert "ligandmpnn" not in kinds(new_task_graph("template_peptide_idr"))
    assert "ligandmpnn" in kinds(new_task_graph("template_peptide_idr", {"methods": ["kejia_idr"]}))
    fold = new_task_graph("new_protein_fold")
    assert fold["nodes"][0]["data"]["params"]["task"] == "unconditional"
    for task in ("ligand_conditioned_pocket", "catalytic_motif_scaffolding"):
        rfd4 = next(n for n in new_task_graph(task)["nodes"] if n["data"]["kind"] == "rfd4")
        assert rfd4["data"]["params"]["options_json"] == "{}"
    selected = new_task_graph("targeted_protein_binder", {"selection": {"ranking_metric": "common.interface_confidence", "limit": 10},
        "additional_validation": {"engine": "rf3", "settings": {"diffusion_batch_size": 1}}})
    assert kinds(selected)[-3:] == ["filter", "rf3", "results"]
    analyses = new_task_graph("targeted_protein_binder", {"interface_analysis": ["ipsae"]})
    assert next(n for n in analyses["nodes"] if n["data"]["kind"] == "boltz2")["data"]["params"]["generate_ipsae"] is False


@pytest.mark.parametrize("template", sorted(graph.SUPPORTED_TEMPLATE_IDS - set(TASKS) - {"blank"}))
def test_legacy_constructors_and_saved_roundtrips_remain_available(template):
    original = graph.new_graph(template)
    cleaned = graph.validate_graph(json.loads(json.dumps(original)))["graph"]
    assert kinds(cleaned) == kinds(original)
    assert cleaned["nodes"] == graph.clean_graph(original)["nodes"]
    if template.endswith("_generated"):
        assert kinds(cleaned)[0] in {"generation_source", "backbone_source"}
        assert template.removesuffix("_generated") not in kinds(cleaned)


@pytest.mark.parametrize("comparison", [True, False])
def test_failed_method_does_not_block_successful_independent_results(execution_harness, monkeypatch, comparison):
    base, calls, _ = execution_harness
    inputs = base["nodes"][0]["data"]["params"]
    result = new_task_graph("targeted_protein_binder", {"methods": ["rfd3", "rfd4"], "comparison": True, "inputs": inputs,
        "method_settings": {"rfd3": {"seed": 42}, "rfd4": {"seed": 43}}})
    if not comparison:
        result.pop("task_configuration")  # Old/custom graphs retain full failure semantics.
    original_submit = runner._submit
    def submit(request, parent, node, **kwargs):
        if node["data"]["kind"] == "rfd3":
            raise RuntimeError("Representative independent generator failure")
        return original_submit(request, parent, node, **kwargs)
    monkeypatch.setattr(runner, "_submit", submit)
    run = runner.run_graph_pipeline(result, "comparison", "owner", "project", reuse_completed=False)
    payload = json.loads((run / "result.json").read_text())
    assert payload["success"] is comparison
    assert payload["outputs"]["candidates"]
    assert payload["outputs"]["comparison_status"] == ("partial" if comparison else "failed")
    assert calls == ["rfd4", "boltz2"]
    evidence = payload["outputs"]["candidates"][0]["raw_metadata"]
    if comparison:
        assert evidence["pipeline_task"] == "targeted_protein_binder"
        assert any(item["kind"] == "rfd4" for item in evidence["pipeline_method_nodes"])
    states = payload["outputs"]["pipeline_nodes"]
    assert states["method-0-rfd3-stage-2"]["status"] == "failed"
    assert states["method-1-rfd4-results"]["status"] == "completed"
    assert states["method-1-rfd4-stage-2"]["execution_specification"]["effective_parameters"]["seed"] == 43


def test_comparison_branches_are_disjoint_and_settings_remain_method_specific():
    result = new_task_graph("targeted_protein_binder", {"methods": ["rfd3", "rfd4"], "comparison": True,
        "num_designs": 3, "method_settings": {"rfd3": {"seed": 19}, "rfd4": {"seed": 23}}})
    assert all(e["source"].split("-", 2)[:2] == e["target"].split("-", 2)[:2] for e in result["edges"])
    generators = [n["data"]["params"] for n in result["nodes"] if n["data"]["kind"] in {"rfd3", "rfd4"}]
    assert [p["seed"] for p in generators] == [19, 23]
    assert [p["num_designs"] for p in generators] == [3, 3]
    assert kinds(result).count("results") == 2


def test_request_preparation_failure_is_isolated_too(execution_harness, monkeypatch):
    base, calls, _ = execution_harness
    result = new_task_graph("targeted_protein_binder", {"methods": ["rfd3", "rfd4"], "comparison": True,
        "inputs": base["nodes"][0]["data"]["params"]})
    prepare = runner.node_request
    def request(kind, *args, **kwargs):
        if kind == "rfd3":
            raise ValueError("Method-specific request preparation failed")
        return prepare(kind, *args, **kwargs)
    monkeypatch.setattr(runner, "node_request", request)
    run = runner.run_graph_pipeline(result, "preparation comparison", "owner", "project", reuse_completed=False)
    payload = json.loads((run / "result.json").read_text())
    assert payload["success"] is True
    assert payload["outputs"]["comparison_status"] == "partial"
    assert payload["outputs"]["candidates"]
    assert calls == ["rfd4", "boltz2"]


def test_comparison_with_all_methods_failed_remains_failed(execution_harness, monkeypatch):
    base, _, _ = execution_harness
    result = new_task_graph("targeted_protein_binder", {"methods": ["rfd3", "rfd4"], "comparison": True,
        "inputs": base["nodes"][0]["data"]["params"]})
    def submit(*args, **kwargs):
        raise RuntimeError("No method succeeded")
    monkeypatch.setattr(runner, "_submit", submit)
    run = runner.run_graph_pipeline(result, "failed comparison", "owner", "project", reuse_completed=False)
    payload = json.loads((run / "result.json").read_text())
    assert payload["success"] is False
    assert payload["outputs"]["candidates"] == []


def test_failed_measurement_review_rows_are_not_a_successful_comparison_branch(execution_harness, monkeypatch):
    base, _, _ = execution_harness
    result = new_task_graph("targeted_protein_binder", {"methods": ["rfd3", "rfd4"], "comparison": True,
        "inputs": base["nodes"][0]["data"]["params"]})
    submit = runner._submit
    def fail_with_review_rows(request, parent, node, **kwargs):
        if node["data"]["kind"] == "rfd3":
            raise RuntimeError("Generator failure")
        child = submit(request, parent, node, **kwargs)
        if node["data"]["kind"] == "boltz2":
            metadata = json.loads((child / "metadata.json").read_text())
            metadata["status"] = "failed"
            (child / "metadata.json").write_text(json.dumps(metadata))
            (child / "result.json").write_text(json.dumps({"success": False, "scientific_status": "failed"}))
        return child
    monkeypatch.setattr(runner, "_submit", fail_with_review_rows)
    run = runner.run_graph_pipeline(result, "measurement failure comparison", "owner", "project", reuse_completed=False)
    payload = json.loads((run / "result.json").read_text())
    assert payload["success"] is False
    assert payload["outputs"]["candidates"]  # Failed rows remain inspectable.
    assert payload["outputs"]["successful_result_nodes"] == []


@pytest.mark.parametrize("fail_method", [False, True])
def test_multi_engine_prediction_has_independent_results_and_engine_settings(execution_harness, monkeypatch, fail_method):
    result = new_task_graph("molecular_structure_prediction", {"methods": ["boltz2", "prediction_esmfold"],
        "comparison": True, "inputs": {"target_protein_sequences": "GGGG"},
        "method_settings": {"boltz2": {"num_recycles": 2, "seed": 42}, "prediction_esmfold": {"num_recycles": 3, "seed": 43}}})
    submit = runner._submit
    def maybe_fail(request, parent, node, **kwargs):
        if fail_method and node["data"]["kind"] == "prediction_esmfold":
            raise RuntimeError("Independent prediction engine failed")
        return submit(request, parent, node, **kwargs)
    monkeypatch.setattr(runner, "_submit", maybe_fail)
    run = runner.run_graph_pipeline(result, "prediction comparison", "owner", "project", reuse_completed=False)
    payload = json.loads((run / "result.json").read_text())
    assert payload["success"] is True
    candidates = payload["outputs"]["candidates"]
    assert len(candidates) == (1 if fail_method else 2)
    assert len({row["candidate_id"] for row in candidates}) == len(candidates)
    assert payload["outputs"]["comparison_status"] == ("partial" if fail_method else "completed")
    boltz = payload["outputs"]["pipeline_nodes"]["method-0-boltz2-prediction"]["execution_specification"]["effective_parameters"]
    assert boltz["seed"] == 42 and boltz["num_recycles"] == 2


def test_api_supports_task_configuration_and_old_identifiers(api_data):
    client, project = api_data["client"], api_data["project"]["id"]
    listing = client.get("/api/v1/pipelines/templates").json()
    assert len(listing["tasks"]) == 8 and len(listing["legacy_items"]) == 22
    request = {"project_id": project, "template": "molecular_structure_prediction", "configuration": {
        "mode": "protein_ligand", "methods": ["boltz2"], "inputs": {"target_protein_sequences": "GGGG", "ligand_smiles": "CCO"}}}
    response = client.post("/api/v1/pipelines/drafts", json=request)
    assert response.status_code == 200, response.text
    assert kinds(response.json()["graph"]) == ["cofold_input", "boltz2", "results"]
    assert client.post("/api/v1/pipelines/drafts", json={"project_id": project, "template": "foundry"}).status_code == 200
    assert client.post("/api/v1/pipelines/drafts", json=request, headers=headers(api_data, "other")).status_code == 404


def test_task_ui_selection_does_not_replace_saved_draft_until_create(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(graph, "app_home", lambda: tmp_path)
    app = AppTest.from_string('''
from mn_molecular_design_studio.pipeline_graph_ui import render_pipeline_builder
render_pipeline_builder(owner_id="test-user", project_selector=lambda _:None,
    submit_api=lambda *a:{}, search_jobs=lambda **k:{"items":[]})
''').run()
    assert not app.exception
    prefix = "visual-pipeline-unassigned"
    original = copy.deepcopy(app.session_state[prefix + "-graph"])
    app.selectbox(key=prefix + "-task").set_value("molecular_structure_prediction").run()
    assert not app.exception
    assert app.session_state[prefix + "-graph"] == original
    app.button(key=prefix + "-create-task").click().run()
    assert not app.exception
    assert kinds(app.session_state[prefix + "-graph"]) == ["cofold_input", "boltz2", "results"]


def test_task_ui_creates_every_task_and_blank_without_starting_jobs():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.scientific_tasks_ui import render_task_setup
graph = render_task_setup("task-test")
if graph is not None:
    st.session_state["created"] = graph
''').run()
    for task in [*TASKS, "blank"]:
        app.selectbox(key="task-test-task").set_value(task).run()
        assert not app.exception, task
        app.button(key="task-test-create-task").click().run()
        assert not app.exception, task
        assert app.session_state["created"]["nodes"]
        if task != "blank":
            assert app.session_state["created"]["task_configuration"]["task_id"] == task
