import json

import pytest
import yaml

from mn_molecular_design_studio.binder_pipeline import ROUTES, compatible_routes, new_binder_pipeline, stage_idr_target
from mn_molecular_design_studio.pipeline_graph import validate_graph
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request


def _settings(graph, kind):
    return next(node["data"]["params"] for node in graph["nodes"] if node["data"]["kind"] == kind)


@pytest.mark.parametrize("route", ROUTES, ids=lambda route: route.id)
def test_registered_binder_routes_build_valid_graphs(route):
    graph = new_binder_pipeline(route.id, target_input=route.target_input, binder_format=route.format)
    kinds = [node["data"]["kind"] for node in graph["nodes"]]
    assert kinds[0] == "target" and kinds[-1] == "results"
    assert route.method in kinds
    if route.method == "bindcraft2":
        assert "shared_filter" not in kinds
        assert _settings(graph, route.method)["modality"] == ("binder" if route.format == "miniprotein" else route.format)
    else:
        assert _settings(graph, "shared_filter")["binder_format"] == route.format
    if route.method == "rfd4":
        external = route.id == "rfd4-external"
        assert ("boltz2" in kinds) == external
        assert json.loads(_settings(graph, "rfd4")["options_json"])["refold"] == (not external)


def test_route_compatibility_and_branch_preserve_existing_target():
    assert not compatible_routes("idr", "cyclic_peptide")
    assert not compatible_routes("sequence", "vhh")
    with pytest.raises(ValueError, match="does not support"):
        new_binder_pipeline("bindcraft2-cyclic", target_input="idr", binder_format="cyclic_peptide")
    first = new_binder_pipeline("bindcraft2", target_input="structure", binder_format="miniprotein",
                                inputs={"target_chains": "B", "hotspots": "B12"})
    target_id = first["nodes"][0]["id"]
    second = new_binder_pipeline("rfd4-internal", target_input="structure", binder_format="miniprotein",
                                 graph=first, target_node_id=target_id)
    assert first["nodes"] == second["nodes"][:len(first["nodes"])]
    assert first["edges"] == second["edges"][:len(first["edges"])]
    assert sum(edge["source"] == target_id for edge in second["edges"]) == 2
    assert sum(node["data"]["kind"] == "target" for node in second["nodes"]) == 1
    validate_graph(second)


def test_idr_request_regenerates_conditions_when_target_or_binder_length_changes(tmp_path, monkeypatch):
    from mn_molecular_design_studio import auth

    monkeypatch.setattr(auth, "user_input_root", lambda owner: tmp_path)
    target = stage_idr_target(">example\nACDEFGHIKLMN", "owner", hotspots="A3")
    graph = new_binder_pipeline("rfd4-idr", target_input="idr", binder_format="miniprotein", inputs=target)
    params = _settings(graph, "rfd4")
    params.update(length_min=40, length_max=60)
    request = node_request("rfd4", params, {"target": target}, "IDR test")
    submitted = request["parameters"]
    assert submitted["target_chains"] == ["A"]
    assert submitted["input_mode"] == "condition_spec"
    assert yaml.safe_load(submitted["condition_spec"])
    assert "C_SEQ" in submitted["condition_spec"]
    assert "C_CRD" not in submitted["condition_spec"]
    assert "40-60" in submitted["condition_spec"]
    with pytest.raises(ValueError, match="does not match"):
        node_request("rfd4", params, {"target": {**target, "target_sequence": "ACDEFGHIKLMR"}}, "IDR test")
    validate_graph(graph, require_inputs=True)
