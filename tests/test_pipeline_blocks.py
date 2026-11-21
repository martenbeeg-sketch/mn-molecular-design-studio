import copy

import pytest

from mn_molecular_design_studio.pipeline_blocks import SCREEN, add_block, blank_program, connect_blocks, next_blocks, program
from mn_molecular_design_studio.pipeline_graph import new_graph, validate_graph


def test_screen_reads_as_native_repeat_and_decision():
    graph = new_graph(SCREEN)
    nodes = {node["data"]["kind"]: node for node in graph["nodes"]}
    nodes["active_decision"]["data"]["params"]["max_iterations"] = 7
    nodes["active_accumulate"]["data"]["params"]["retained_pose_count"] = 2000
    before = copy.deepcopy(graph)

    view = program(graph)
    screen = view["screen"]
    assert screen["retention"] == 2000
    assert screen["stop_conditions"][-1] == "The campaign reaches 7 rounds"
    assert screen["decision"] == nodes["active_decision"]["id"]
    assert screen["finish"] == [nodes["active_results"]["id"]]
    assert screen["next_round"][-1] == nodes["active_select"]["id"]
    assert graph == before
    assert next(edge for edge in graph["edges"] if edge.get("loop_back"))["target"] == screen["rounds"][0]
    validate_graph(graph)


def test_blocks_can_build_a_workflow_with_typed_filter_branches():
    graph, source = add_block(blank_program(), "cofold_input")
    graph["nodes"][0]["data"]["params"].update(
        task_type="protein_ligand", target_protein_sequences="ACDE", ligand_smiles="CCO",
    )
    possible, _ = next_blocks(graph, source, "out")
    assert ("boltz2", "in") in possible
    assert ("af2", "in") not in possible
    with pytest.raises(ValueError, match="downstream"):
        validate_graph(graph, require_inputs=True, allow_incomplete=True)

    graph, predictor = add_block(graph, "boltz2", source=source)
    graph, decision = add_block(graph, "filter", source=predictor)
    graph, results = add_block(graph, "results", source=decision)
    _, existing = next_blocks(graph, decision, "fail")
    assert (results, "in") in existing
    assert (predictor, "in") not in existing
    graph = connect_blocks(graph, decision, "fail", results, "in")
    graph = connect_blocks(graph, decision, "unresolved", results, "in")
    assert program(graph)["blocks"][results]["sources"] == [
        "Filter candidates · Pass", "Filter candidates · Fail", "Filter candidates · Unresolved",
    ]
    validate_graph(graph, require_inputs=True, require_files=False)


def test_pipeline_builder_opens_screen_and_blocks_unconfigured_run(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import pipeline_graph
    from mn_molecular_design_studio.navigation import WORKSPACE_PAGES

    pages = [item[2] for item in WORKSPACE_PAGES]
    assert "Pipelines" in pages and "Pipeline Builder" in pages
    monkeypatch.setattr(pipeline_graph, "app_home", lambda: tmp_path)
    app = AppTest.from_string("""
from mn_molecular_design_studio.pipeline_blocks_ui import render_pipeline_blocks
render_pipeline_blocks(owner_id="test", project_selector=lambda _: "project",
    submit_api=lambda *args: (_ for _ in ()).throw(AssertionError("Nothing should run")),
    search_jobs=lambda **kwargs: {"items": []})
""").run()
    assert not app.exception
    app.button(key="pipeline-blocks-project-screen").click().run()
    assert not app.exception
    assert app.session_state["pipeline-blocks-project-graph"]["template_id"] == SCREEN
    assert app.button(key="pipeline-blocks-project-run").disabled
    assert any("fingerprint dataset" in message.value for message in app.info)


def test_pipeline_builder_submits_the_existing_graph_and_opens_its_results(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import pipeline_graph

    monkeypatch.setattr(pipeline_graph, "app_home", lambda: tmp_path)
    requests = []

    def submit(path, payload):
        requests.append((path, payload))
        return {"job": {"cluster_job_id": "test-job", "job_code": "TEST"}}

    app = AppTest.from_function(_ready_program, args=(submit,)).run()
    assert not app.exception
    assert not requests
    app.button(key="pipeline-blocks-project-run").click().run()
    assert not app.exception
    assert len(requests) == 1
    path, payload = requests[0]
    assert path == "/api/v1/pipelines/run"
    assert payload["project_id"] == "project"
    validate_graph(payload["graph"], require_inputs=True, require_files=False)
    app.button(key="pipeline-blocks-project-results").click().run()
    assert not app.exception
    assert app.session_state["selected-job-reference"] == "cluster:test-job"


def _ready_program(submit):
    import streamlit as st
    from mn_molecular_design_studio.pipeline_blocks_ui import render_pipeline_blocks
    from mn_molecular_design_studio.scientific_tasks import new_task_graph

    key = "pipeline-blocks-project-graph"
    if key not in st.session_state:
        st.session_state[key] = new_task_graph("molecular_structure_prediction", {
            "inputs": {"target_protein_sequences": "ACDE"},
        })
    render_pipeline_blocks(owner_id="test", project_selector=lambda _: "project",
                           submit_api=submit, search_jobs=lambda **kwargs: {"items": []})
