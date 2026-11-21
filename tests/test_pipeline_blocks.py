import copy

import pytest

from mn_molecular_design_studio.pipeline_blocks import (
    SCREEN, add_block, blank_program, canvas_program, connect_blocks, drop_block, move_canvas_block,
    next_blocks, place_canvas_block, placements, program,
)
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
    assert screen["selection_id"] == screen["next_round"][-1]
    assert screen["selection"]["example"] == "Up to 5,000 with 10,000 initial centroids"
    assert screen["selection"]["predicted_count"] == "4,000"
    assert screen["selection"]["random_count"] == "1,000"
    assert screen["finish"] == [nodes["active_results"]["id"]]
    assert screen["next_round"][-1] == nodes["active_select"]["id"]
    assert graph == before
    assert next(edge for edge in graph["edges"] if edge.get("loop_back"))["target"] == screen["rounds"][0]
    validate_graph(graph)


def test_screen_reports_edited_seed_mix_size_cap_and_cleanup():
    graph = new_graph(SCREEN)
    nodes = {node["data"]["kind"]: node for node in graph["nodes"]}
    nodes["active_seed"]["data"]["params"].update(sample_size=2_000, sampling_seed=67)
    nodes["active_docking"]["data"]["params"].update(batch_size=500, cleanup_after_round=False)
    nodes["active_select"]["data"]["params"].update(next_round_seed_fraction_percent=75, predicted_fraction_percent=60)
    nodes["active_accumulate"]["data"]["params"]["retained_pose_count"] = 0
    screen = program(graph)["screen"]

    assert screen["selection"] == {
        "size": "75% of the initial centroid sample, limited by the number of unattempted entries and capped at 500.",
        "example": "Up to 500 with 2,000 initial centroids",
        "predicted": "60% highest model probabilities", "predicted_count": "300",
        "random": "40% uniform random from the rest", "random_count": "200",
    }
    assert screen["step_details"][nodes["active_seed"]["id"]] == "Round 1 · 2,000 randomly sampled centroids · seed 67."
    assert screen["step_details"][nodes["active_accumulate"]["id"]].endswith("No campaign poses retained. Keep worker files.")


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


def test_drop_inserts_between_two_blocks_and_preserves_other_wires():
    from mn_molecular_design_studio.scientific_tasks import new_task_graph

    graph = new_task_graph("molecular_structure_prediction", {
        "inputs": {"target_protein_sequences": "ACDE"},
    })
    source, predictor, results = [node["id"] for node in graph["nodes"]]
    original_first, original_last = graph["edges"]
    slots, cards = placements(graph, predictor)
    socket = f"between/{original_last['id']}"
    assert socket in cards["filter/in"]["targets"]
    with pytest.raises(ValueError, match="does not fit"):
        drop_block(graph, slots, cards, socket, "cofold_input/")

    updated, filter_id = drop_block(graph, slots, cards, socket, "filter/in")
    assert graph["edges"] == [original_first, original_last]
    assert original_first in updated["edges"]
    assert [(edge["source"], edge["target"]) for edge in updated["edges"]] == [
        (source, predictor), (predictor, filter_id), (filter_id, results),
    ]
    validate_graph(updated, require_inputs=True, require_files=False)


def test_screen_recipe_snaps_into_an_empty_program():
    empty = blank_program()
    slots, cards = placements(empty, "")
    graph, selected = drop_block(empty, slots, cards, "start", f"{SCREEN}/")
    assert selected in program(graph)["screen"]["setup"]
    assert placements(graph, selected) == ({}, {})
    after_input, selected = drop_block(empty, slots, cards, "start", "cofold_input/")
    assert f"{SCREEN}/" not in placements(after_input, selected)[1]


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


def test_pipeline_builder_receives_a_palette_drop_without_submitting(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import pipeline_graph, pipeline_blocks_ui

    monkeypatch.setattr(pipeline_graph, "app_home", lambda: tmp_path)
    monkeypatch.setattr(pipeline_blocks_ui, "_COMPONENT", lambda **props: {
        "action": "place", "card": "cofold_input", "position": {"x": 20, "y": 50}, "token": "single-drop",
    })
    app = AppTest.from_string("""
from mn_molecular_design_studio.pipeline_blocks_ui import render_pipeline_blocks
render_pipeline_blocks(owner_id="test", project_selector=lambda _: "project",
    submit_api=lambda *args: (_ for _ in ()).throw(AssertionError("Nothing should run")),
    search_jobs=lambda **kwargs: {"items": []})
""").run()
    assert not app.exception
    nodes = app.session_state["pipeline-blocks-project-graph"]["nodes"]
    assert [item["data"]["kind"] for item in nodes] == ["cofold_input"]
    assert nodes[0]["position"] == {"x": 20, "y": 50}


def test_canvas_snap_and_move_preserve_typed_edges_and_saved_canvas_positions():
    from mn_molecular_design_studio.pipeline_graph import clean_graph

    graph, source = place_canvas_block(blank_program(), "cofold_input", {"x": 20, "y": 50})
    graph, predictor = place_canvas_block(graph, "boltz2", {"x": 260, "y": 50}, source=source, output="out")
    graph, results = place_canvas_block(graph, "results", {"x": 200, "y": 260})
    with pytest.raises(ValueError, match="specific input"):
        place_canvas_block(graph, "cofold_input", {"x": 40, "y": 50}, source=predictor, output="out")
    updated = connect_blocks(graph, predictor, "out", results, "in")
    position = {"x": 530, "y": 46}
    moved = clean_graph(move_canvas_block(updated, results, position))
    view = canvas_program(moved)
    assert view["nodes"][results]["position"] == position
    assert view["edges"] == updated["edges"]
    assert [(edge["source"], edge["target"]) for edge in view["edges"]] == [(source, predictor), (predictor, results)]
    assert ("results", "in") in view["compatible"]["boltz2"]["out"]
    with pytest.raises(ValueError, match="valid position"):
        move_canvas_block(updated, results, {"x": True, "y": 0})


def test_canvas_displays_the_native_screen_feedback_edge():
    graph, _ = place_canvas_block(blank_program(), SCREEN, {"x": 20, "y": 50})
    view = canvas_program(graph)
    assert view["managed"] and not view["compatible"]
    loop, = [edge for edge in view["edges"] if edge.get("loop_back")]
    assert view["nodes"][loop["source"]]["kind"] == "active_select"
    assert view["nodes"][loop["target"]]["kind"] == "active_docking"


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
