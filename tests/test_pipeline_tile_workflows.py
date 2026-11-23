import json

import pytest

from mn_molecular_design_studio import pipeline_tile_workflows as tiles
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import (
    _cache_valid, _tile_output, node_request,
)
from mn_molecular_design_studio.pipeline_graph import connection_compatibility
from mn_molecular_design_studio.run_graph import source_request
from mn_molecular_design_studio.tool_connections import connection_catalog


def test_visible_tiles_have_registered_canvas_operations():
    result = connection_catalog(view="tiles", query="", operation="", direction="both",
                                artifact_type="", offset=0, limit=200)
    assert result["total"] == 95
    assert not [tile for tile in result["items"] if not tile["pipeline_operations"]]


def test_compound_and_receptor_ports_wire_exact_artifacts():
    compounds = {"tile_artifact": {"path": "/test/compounds.json"}}
    target = {"target": {"target_pdb": "/test/receptor.pdb", "target_chains": ["A"]}}
    specification = {"workflow_parameters": {"center_x": 0, "center_y": 0, "center_z": 0}}
    request = node_request("tile_docking_vina", specification, {
        "workflow_inputs": {"target": target, "compounds": compounds},
    }, "Vina")
    assert request["parameters"] == {
        "target_pdb": "/test/receptor.pdb", "compound_set_file": "/test/compounds.json",
        "center_x": 0, "center_y": 0, "center_z": 0, "job_name": "Vina",
    }
    assert connection_compatibility(
        "tile_generation_omtra", "tile_docking_vina", target_handle="compounds",
    )["compatible"]
    assert not connection_compatibility(
        "tile_generation_omtra", "tile_docking_vina", target_handle="target",
    )["compatible"]


def test_rosetta_ligand_reference_is_a_separate_canvas_input():
    from mn_molecular_design_studio.pipeline_graph import capability_contract

    kind = "tile_docking_rosetta_ligand"
    assert [port["id"] for port in capability_contract(kind)["inputs"]] == [
        "target", "compounds", "reference_ligand",
    ]
    assert connection_compatibility("ligand_sdf_input", kind, target_handle="reference_ligand")["compatible"]
    assert not connection_compatibility("ligand_sdf_input", kind, target_handle="compounds")["compatible"]
    settings = {"center_x": 0, "center_y": 0, "center_z": 0}
    sources = {"target": {"target": {"target_pdb": "/test/receptor.pdb"}},
               "compounds": {"tile_artifact": {"path": "/test/compounds.json"}},
               "reference_ligand": {"tile_artifact": {"path": "/test/reference.sdf"}}}
    request = node_request(kind, {"workflow_parameters": settings}, {"workflow_inputs": sources}, "Dock")
    assert request["parameters"] == {
        **settings, "target_pdb": "/test/receptor.pdb", "compound_set_file": "/test/compounds.json",
        "reference_ligand_sdf": "/test/reference.sdf", "job_name": "Dock",
    }
    sources.pop("reference_ligand")
    request = node_request(kind, {"workflow_parameters": settings}, {"workflow_inputs": sources}, "Dock")
    assert "reference_ligand_sdf" not in request["parameters"]


def test_ligand_sdf_wires_use_the_workflows_registered_input_fields():
    from mn_molecular_design_studio.pipeline_graph import capability_contract

    ligand = {"tile_artifact": {"path": "/test/ligand.sdf"}}
    complex_ = {"tile_artifact": {"path": "/test/complex.pdb"}}
    for kind in ("tile_md_openmm", "tile_md_gromacs"):
        assert connection_compatibility("ligand_sdf_input", kind, target_handle="ligand")["compatible"]
        assert [port["id"] for port in capability_contract(kind)["inputs"]] == ["complex", "ligand"]
        config = {"ligand_key": "LIG"}
        tiles.validate_parameters(tiles.TILES[kind], config, {"complex", "ligand"}, runnable=True)
        request = node_request(kind, {"workflow_parameters": config},
                               {"workflow_inputs": {"complex": complex_, "ligand": ligand}}, "MD")
        assert request["parameters"]["complex_pdb"] == "/test/complex.pdb"
        assert request["parameters"]["ligand_sdf"] == "/test/ligand.sdf"
    kind = "tile_generation_lddm"
    assert connection_compatibility("ligand_sdf_input", kind, target_handle="reference_ligand")["compatible"]
    assert tiles.connected_parameter_names(tiles.TILES[kind], "reference_ligand") == {"reference_ligand_sdf"}


def test_canvas_workflow_form_reads_the_registry_and_keeps_wired_fallbacks():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.pipeline_node_ui import _render_tile_workflow_settings
from mn_molecular_design_studio.pipeline_tile_workflows import TILES
target = {"id": "target", "data": {"label": "Receptor"}}
node = {"id": "docking", "data": {"label": "Vina"}}
graph = {"nodes": [target, node], "edges": [{"source": "target", "target": "docking", "targetHandle": "target"}]}
settings = {"target_pdb": "/test/fallback.pdb", "center_x": 0.0}
st.session_state["values"], st.session_state["ready"] = _render_tile_workflow_settings(
    TILES["tile_docking_vina"], settings, node, graph, "tool")
""").run()
    assert not app.exception
    assert app.number_input(key="tool-center_x").value == 0.0
    assert app.number_input(key="tool-center_y").value is None
    assert not any(widget.key == "tool-target_pdb" for widget in app.text_input)
    assert app.session_state["values"] == {"target_pdb": "/test/fallback.pdb", "center_x": 0.0}
    app.number_input(key="tool-center_y").set_value(2.5).run()
    assert not app.exception
    assert app.session_state["values"]["center_y"] == 2.5


def test_manifest_handoff_requires_real_artifact_and_checks_its_digest(tmp_path):
    run = tmp_path / "ligand-generation" / "fixture"
    run.mkdir(parents=True)
    data = run / "compound_set.json"
    data.write_text('{"compounds": []}')
    (run / "metadata.json").write_text(json.dumps({"status": "completed"}))
    (run / "result.json").write_text(json.dumps({
        "success": True, "outputs": {"generated_compound_set": "compound_set.json"},
    }))
    output = _tile_output("tile_generation_omtra", run, {}, owner_id="user", project_id="project")
    assert output["tile_artifact"]["path"] == str(data)
    assert _cache_valid(output)
    data.write_text('{"compounds": ["changed"]}')
    assert not _cache_valid(output)
    data.unlink()
    with pytest.raises(ValueError, match="unavailable"):
        _tile_output("tile_generation_omtra", run, {}, owner_id="user", project_id="project")


def test_library_registration_feeds_the_recorded_manifest_and_id(tmp_path):
    output_file = tmp_path / "features.json"
    output_file.write_text('{"library_id": "fixture"}')
    adapter = tiles.TILES["tile_prepare_library"]
    source = tiles.artifact_output(adapter, tmp_path,
                                   {"outputs": {"feature_manifest": str(output_file), "library_id": "fixture"}})
    request = node_request("tile_library_centroids", {
        "workflow_parameters": {"output_directory": "/test/centroids", "library_id": "different"},
    }, {"workflow_inputs": {"features": source}}, "Centroids")
    assert request["parameters"]["library_id"] == "fixture"
    assert request["parameters"]["feature_manifest"] == str(output_file)


def test_pipeline_tools_reject_undeclared_parameters_and_missing_inputs():
    tile = tiles.TILES["tile_md_openmm"]
    with pytest.raises(ValueError, match="unknown_control"):
        tiles.validate_parameters(tile, {"unknown_control": False}, set(), runnable=False)
    with pytest.raises(ValueError, match="ligand_sdf"):
        tiles.validate_parameters(tile, {}, {"complex"}, runnable=True)


@pytest.mark.parametrize("empty", [None, "", "  "])
def test_md_starter_requires_a_ligand_input_before_queuing(empty):
    from mn_molecular_design_studio.pipeline_graph import PipelineInputError, new_graph, validate_graph

    graph = new_graph("md_simulation")
    steps = {node["id"]: node["data"]["params"] for node in graph["nodes"]}
    steps["tile_complex_prepare"]["workflow_parameters"] = {"source_pdb": "complex.pdb"}
    steps["tile_md_openmm"]["workflow_parameters"] = {"ligand_sdf": empty, "ligand_key": "LIG"}

    validate_graph(graph)
    with pytest.raises(PipelineInputError, match="ligand_sdf") as exc:
        validate_graph(graph, require_inputs=True, require_files=False)
    assert exc.value.node_id == "tile_md_openmm"


def test_tool_input_digests_track_explicit_files_and_preserve_inline_values(tmp_path):
    inputs = tmp_path / "inputs.json"
    unrelated = tmp_path / "unrelated.txt"
    unrelated.write_text("unrelated")
    inputs.write_text(json.dumps({"text": str(unrelated)}))
    settings = {"source_files": [str(inputs)], "sequence": str(unrelated), "source_name": str(unrelated)}
    initial = tiles.input_digests(settings, tmp_path)
    assert set(initial) == {str(inputs)}
    inputs.write_text('{"updated": true}')
    assert tiles.input_digests(settings, tmp_path) != initial


def test_tool_input_digests_leave_peer_paths_for_the_owning_worker(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import config

    monkeypatch.setattr(config, "identity", lambda: "local-workstation")
    local = tmp_path / "explicit.cif"
    local.write_text("data_target\n")
    remote = "cluster://peer-workstation/refolding-validation/saved/artifacts/normalized_candidates/candidates.jsonl"
    assert set(tiles.input_digests({"candidates_jsonl": remote, "target_pdb": str(local)}, tmp_path)) == {str(local)}


def test_saved_run_uses_the_exact_candidate_manifest_for_acceptance():
    kind = "tile_validation_acceptance"
    settings = {"acceptance_policy": {"rules": []}}
    source = {"deferred_run_source": True, "run_dir": "/work/design/cohort",
              "candidates_jsonl": "/work/design/cohort/artifacts/normalized_candidates/candidates.jsonl"}
    request = source_request(kind, {"workflow_parameters": settings}, source, "Saved Run")
    assert request["workflow"] == tiles.TILES[kind].workflow
    assert request["parameters"]["source_run_dir"] == source["run_dir"]
    assert request["parameters"]["candidates_jsonl"] == source["candidates_jsonl"]
    assert all(request["parameters"][key] == value for key, value in settings.items())
    assert connection_compatibility("run_source", kind, target_handle="candidates")["compatible"]


def test_saved_run_rejects_outputs_the_batch_collector_cannot_merge():
    assert not connection_compatibility("run_source", "tile_bioemu2", target_handle="candidates")["compatible"]
    assert connection_compatibility("project_design_source", "tile_bioemu2", target_handle="candidates")["compatible"]
    assert not connection_compatibility("run_source", "candidate_merge")["compatible"]


def test_binderscope_requires_an_explicit_native_campaign_on_canvas():
    tile = tiles.TILES["tile_binderscope"]
    for kind in ("run_source", "project_design_source", "tile_validation_acceptance"):
        assert not connection_compatibility(kind, "tile_binderscope", target_handle="candidates")["compatible"]
    assert tile.output == "tool_result"
    with pytest.raises(ValueError, match="Analyze saved database designs from Data Browser"):
        tiles.validate_parameters(tile, {"source_run_dir": "/work/design", "candidates_jsonl": "/work/candidates.jsonl"},
                                  set(), runnable=True)
    settings = {"native_run_dirs": ["/work/bindcraft"], "target_pdb": "/work/target.pdb"}
    assert tiles.workflow_request(tile, settings, {}, "Native campaign")["parameters"] == settings
