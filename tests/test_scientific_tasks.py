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
    assert len(graph.MAIN_STARTERS) == 6
    assert graph.MAIN_STARTERS["Molecular dynamics"] == "md_analysis"
    assert "md_simulation" in graph.STARTER_TEMPLATES.values()
    assert graph.STARTER_TEMPLATES["Multi-Generator Binder Design & Evaluation"] == "multi_generator_binder"
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
    assert len(listing["tasks"]) == 8 and len(listing["legacy_items"]) == 23
    tasks = {item["id"]: item for item in listing["items"]}
    assert {k for k, v in tasks.items() if v["section"] == "main"} == set(graph.MAIN_STARTERS.values())
    assert tasks["native_comparative"]["parent_id"] == "targeted_protein_binder"
    assert tasks["md_simulation"]["parent_id"] == "md_analysis"
    request = {"project_id": project, "template": "molecular_structure_prediction", "configuration": {
        "mode": "protein_ligand", "methods": ["boltz2"], "inputs": {"target_protein_sequences": "GGGG", "ligand_smiles": "CCO"}}}
    response = client.post("/api/v1/pipelines/drafts", json=request)
    assert response.status_code == 200, response.text
    assert kinds(response.json()["graph"]) == ["cofold_input", "boltz2", "results"]
    assert client.post("/api/v1/pipelines/drafts", json={"project_id": project, "template": "foundry"}).status_code == 200
    assert client.post("/api/v1/pipelines/drafts", json=request, headers=headers(api_data, "other")).status_code == 404


def test_pipeline_template_and_trash_are_independent_of_saved_versions(api_data):
    client, project = api_data["client"], api_data["project"]["id"]
    root = "/api/v1/pipelines"
    original = graph.new_graph("unidock_active_learning")
    original["name"] = "Screening setup"
    pipeline_id = original["id"]
    saved = client.post(root + "/save", json={"project_id": project, "graph": original, "version": True})
    assert saved.status_code == 200, saved.text
    revision = saved.json()["revision"]
    publish = {"project_id": project, "expected_revision": revision, "name": "Reusable screen"}
    published = client.post(f"{root}/{pipeline_id}/template", json=publish)
    assert published.status_code == 200, published.text
    assert client.post(f"{root}/{pipeline_id}/template", json=publish).json() == published.json()
    template_id = published.json()["id"]
    assert template_id in {row["id"] for row in client.get(root + "/shared-templates", headers=headers(api_data, "other")).json()["items"]}
    assert client.post(f"{root}/{pipeline_id}/template", json=publish, headers=headers(api_data, "student")).status_code == 403

    edit = client.post(root + "/save", json={"project_id": project, "expected_revision": revision,
                                            "graph": {**original, "name": "Edited screening setup"}})
    assert edit.status_code == 200, edit.text
    trash = {"project_id": project, "expected_revision": edit.json()["revision"], "trashed": True}
    assert client.post(f"{root}/{pipeline_id}/trash", json={**trash, "expected_revision": revision}).status_code == 409
    assert client.post(f"{root}/{pipeline_id}/trash", json=trash).status_code == 200
    assert client.get(root, params={"project_id": project}).json() == []
    listed = client.get(root, params={"project_id": project, "include_trashed": True}).json()
    assert listed[0]["trashed_at"] and listed[0]["versions"] == 1
    assert client.post(root + "/run", json={"project_id": project, "graph": original}).status_code == 409
    assert client.get(f"{root}/{pipeline_id}/versions", params={"project_id": project}).json()["items"][0]["revision"] == revision

    copy = client.post(f"{root}/shared-templates/{template_id}/use", headers=headers(api_data, "other"))
    assert copy.status_code == 200, copy.text
    new = copy.json()["graph"]
    assert new["id"] != pipeline_id
    assert new["name"] == "Reusable screen"
    assert new["nodes"] == original["nodes"] and new["edges"] == original["edges"]
    assert client.post(f"{root}/{pipeline_id}/trash", json={**trash, "trashed": False}).status_code == 200
    assert client.get(root, params={"project_id": project}).json()[0]["graph"]["id"] == pipeline_id


def test_shared_starter_update_is_visible_without_project_access(api_data):
    client = api_data["client"]
    base = "/api/v1/pipelines/shared-templates"
    template_id = "unidock_active_learning"
    reader = headers(api_data, "other")
    starter = client.get(f"{base}/{template_id}", headers=reader).json()
    edited = copy.deepcopy(starter["graph"])
    edited.update(name="Shared screening recipe", viewport={"x": 74, "y": -49, "zoom": .68})
    edited["nodes"][0]["position"] = {"x": 156, "y": 252}
    for node in edited["nodes"]:
        if node["data"]["kind"] == "active_seed":
            node["data"]["params"].update(sample_size=137, library_id="private-library",
                                            centroid_manifest="/private/centroids.json")
    request = {"template_id": template_id, "expected_revision": starter["revision"], "graph": edited}
    assert client.post(f"{base}/save", json=request, headers=reader).status_code == 403
    assert not graph.can_edit_shared_templates({"role": "admin", "username": "another_admin"})
    result = client.post(f"{base}/save", json=request)
    assert result.status_code == 200, result.text
    saved = result.json()
    assert saved["revision"] != starter["revision"]
    assert client.post(f"{base}/save", json=request).json() == saved
    assert client.get(f"{base}/{template_id}", headers=reader).json() == saved
    assert [item["name"] for item in client.get(base, headers=reader).json()["items"]
            if item["id"] == template_id] == ["Shared screening recipe"]
    assert client.post(f"{base}/save", json={**request, "name": "Stale edit"}).status_code == 409
    assert client.post(f"{base}/save", json={**request, "project_id": "any"}).status_code == 422
    fork = client.post(f"{base}/{template_id}/use", headers=reader).json()["graph"]
    assert fork["id"] != saved["graph"]["id"]
    assert fork["viewport"] == edited["viewport"]
    assert fork["nodes"][0]["position"] == edited["nodes"][0]["position"]
    assert fork["name"] == edited["name"]
    seed = next(node["data"]["params"] for node in fork["nodes"] if node["data"]["kind"] == "active_seed")
    assert (seed["sample_size"], seed["library_id"], seed["centroid_manifest"]) == (137, "", "")


def test_shared_protein_design_routes_open_and_save_independently(api_data):
    client = api_data["client"]
    reader = headers(api_data, "other")
    base = "/api/v1/pipelines/shared-templates"
    items = client.get(base, headers=reader).json()["items"]
    protein, = [item for item in items if item["id"].startswith("protein_design")]
    assert protein["name"] == "Protein design"
    assert [point["label"] for point in protein["starting_points"]] == [
        "Design against a target", "Design a new fold", "Redesign an existing structure",
    ]
    previews = {}
    for point in protein["starting_points"]:
        for route in point["routes"]:
            template_id = route["id"]
            preview = client.get(f"{base}/{template_id}", headers=reader).json()
            previews[template_id] = preview
            graph.validate_graph(preview["graph"])
            assert kinds(preview["graph"]).count("results") == 1
            assert client.get(f"{base}/{template_id}", headers=reader).json()["revision"] == preview["revision"]
            used = client.post(f"{base}/{template_id}/use", headers=reader).json()["graph"]
            assert used["id"] != preview["graph"]["id"]
            assert used["nodes"] == preview["graph"]["nodes"]

    assert kinds(previews["protein_design_native"]["graph"]) == ["target", "binder_context", "bindcraft2", "results"]
    comparison = kinds(previews["protein_design_comparison"]["graph"])
    assert comparison.count("target") == 1 and "bindcraft2" in comparison and "rfd4" in comparison
    assert kinds(previews["protein_design_new_fold"]["graph"]) == ["protpardelle", "ligandmpnn", "esmfold", "results"]
    redesign = previews["protein_design_redesign"]
    assert kinds(redesign["graph"]) == ["target", "esm_if1", "esmfold", "results"]
    edited = copy.deepcopy(redesign["graph"])
    edited["nodes"][0]["position"]["y"] += 24
    request = {"template_id": redesign["id"], "expected_revision": redesign["revision"], "graph": edited}
    assert client.post(f"{base}/save", json=request, headers=reader).status_code == 403
    saved = client.post(f"{base}/save", json=request)
    assert saved.status_code == 200, saved.text
    assert client.get(f"{base}/{redesign['id']}", headers=reader).json() == saved.json()
    assert client.get(f"{base}/protein_design", headers=reader).json() == previews["protein_design"]
    assert client.get(base, headers=reader).json()["items"] == items


def test_pipeline_templates_and_guided_setup_open_separate_canvases(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import pipeline_graph_ui

    monkeypatch.setattr(graph, "app_home", lambda: tmp_path)

    def canvas(**args):
        action = pipeline_graph_ui.st.session_state.pop("canvas-test-action", None)
        if action:
            assert args["guided_setup"]
            return {"graph": args["graph"], "action": {"type": action, "id": "guided-setup"}}

    monkeypatch.setattr(pipeline_graph_ui, "_component", lambda: canvas)
    app = AppTest.from_string("""
from mn_molecular_design_studio.pipeline_graph_ui import render_pipeline_builder
render_pipeline_builder(owner_id="test-user", project_selector=lambda _:None,
    submit_api=lambda *a: (_ for _ in ()).throw(AssertionError("No request expected")),
    search_jobs=lambda **k: {"items": []})
""").run()
    prefix = "visual-pipeline-unassigned"
    app.segmented_control(key=prefix + "-library-section").set_value("Templates").run()
    assert not app.exception
    assert not app.selectbox
    assert app.button(key=prefix + "-template-unidock_active_learning").label == "Large-library screen"
    app.button(key=prefix + "-new").click().run()
    blank = app.session_state[prefix + "-graph"]
    assert not blank["nodes"]
    app.session_state["canvas-test-action"] = "guided_setup"
    app.run()
    assert not app.exception
    guided = prefix + "-guided-" + blank["id"]
    assert app.selectbox(key=guided + "-task").value is None
    app.selectbox(key=guided + "-task").set_value("advanced").run()
    app.selectbox(key=guided + "-advanced-advanced-recipe").set_value("legacy").run()
    app.selectbox(key=guided + "-advanced-legacy").set_value("Foundry · RFD3 → MPNN → RF3").run()
    app.button(key=guided + "-advanced-open-legacy").click().run()
    assert not app.exception
    composed = app.session_state[prefix + "-graph"]
    assert composed["id"] == blank["id"]
    assert kinds(composed) == ["target", "binder_context", "rfd3", "foundry_mpnn", "rf3", "results"]


def test_task_ui_creates_every_task_and_blank_without_starting_jobs():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.scientific_tasks_ui import render_task_setup
graph = render_task_setup("task-test")
if graph is not None:
    st.session_state["created"] = graph
''').run()
    from mn_molecular_design_studio.pipeline_graph import MAIN_STARTERS
    for task in [*TASKS, "blank"]:
        main = task in MAIN_STARTERS.values()
        app.selectbox(key="task-test-task").set_value(task if main else "advanced").run()
        if not main:
            app.selectbox(key="task-test-advanced-recipe").set_value(task).run()
        assert not app.exception, task
        app.button(key="task-test-create-task").click().run()
        assert not app.exception, task
        assert app.session_state["created"]["nodes"]
        if task != "blank":
            assert app.session_state["created"]["task_configuration"]["task_id"] == task



def test_md_starter_pairs_prepared_complex_and_trajectory_ports():
    simulation = graph.new_graph("md_simulation")
    assert kinds(simulation) == ["tile_complex_prepare", "tile_md_openmm", "md_analysis", "results"]
    assert simulation["edges"][0]["targetHandle"] == "complex"
    plan = graph.validate_graph(simulation)
    assert plan["order"] == kinds(simulation)
