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
