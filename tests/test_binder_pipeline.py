import json

import pytest
import yaml

from mn_molecular_design_studio.binder_pipeline import BC2_FORMATS, ROUTES, compatible_routes, new_binder_pipeline, stage_idr_target, write_idr_fasta
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
        assert _settings(graph, route.method)["modality"] == BC2_FORMATS[route.format]
    elif "shared_filter" in kinds:
        assert _settings(graph, "shared_filter")["binder_format"] == route.format
    else:
        assert (route.method.startswith("tile_design_")
                or route.method in {"boltzgen", "complexa", "peptide_trap", "kejia_idr"}
                or route.format in {"antibody", "tcr"} or route.id == "rfd4-cyclic")
    if route.method == "boltzgen":
        assert _settings(graph, route.method)["task"] == {"vhh": "nanobody", "miniprotein": "protein"}.get(route.format, route.format)
    if route.method == "tile_design_rfantibody":
        parameters = _settings(graph, route.method)["workflow_parameters"]
        assert parameters["task"] == ("nanobody" if route.format == "vhh" else "antibody")
    if route.method == "rfd4":
        external = route.id == "rfd4-external"
        assert ("boltz2" in kinds) == external
        settings = _settings(graph, "rfd4")
        if route.format in {"antibody", "tcr"}:
            assert settings["task"] == route.format
        elif route.id == "rfd4-cyclic":
            assert settings["task"] == "cyclic_peptide_binder" and settings["input_mode"] == "annotated_cif"
            assert json.loads(settings["options_json"])["refold"] is False
        else:
            assert json.loads(settings["options_json"])["refold"] == (not external)


def test_design_capabilities_match_registered_workflows_and_picker_routes():
    from mn_molecular_design_studio.design_capabilities import catalogue
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.tool_catalog import TOOL_CATALOG_WORKFLOW_IDS

    data = catalogue()
    assert len({route.id for route in ROUTES}) == len(ROUTES)
    assert set(data["binder_formats"]) == {row["output"] for row in data["items"]
                                            if row["output"] in data["binder_formats"]}
    assert {row["route_id"] for row in data["items"] if row["access"] == "pipeline"} == {route.id for route in ROUTES}
    assert {row["access"] for row in data["items"] if row["output"] in data["binder_formats"]} == {"pipeline"}
    assert all(NODE_CATALOG[row["canvas_node"]]["workflow"] == row["workflow"]
               for row in data["items"] if row["access"] == "canvas")
    assert {row["access"] for row in data["items"]} == {"pipeline", "canvas"}
    assert {row["workflow"] for row in data["items"]} <= set(WORKFLOWS)
    assert {workflow for workflow in TOOL_CATALOG_WORKFLOW_IDS if workflow.startswith(("design.", "ligand.generation."))} <= {
        row["workflow"] for row in data["items"]}
    cyclic = catalogue(target_input="structure", binder_format="cyclic_peptide")["items"]
    assert {row["engine"] for row in cyclic if row["access"] == "pipeline"} == {
        "BindCraft2", "BoltzGen", "RFdiffusion3", "RFD4-Proteina", "RFpeptides",
    }


def test_cyclic_binder_tile_inherits_connected_target():
    graph = new_binder_pipeline("rfpeptides", target_input="structure", binder_format="cyclic_peptide")
    target = {"target_pdb": "fixture.pdb", "target_chains": ["A"], "hotspots": "A3"}
    kind = "tile_design_rfpeptides"
    request = node_request(kind, _settings(graph, kind), {"workflow_inputs": {"target": {"target": target}}}, "fixture")
    assert request["workflow"] == "design.rfpeptides"
    assert all(request["parameters"][key] == value for key, value in target.items())


def test_multidomain_target_reuses_verified_structural_routes():
    from mn_molecular_design_studio.design_capabilities import catalogue

    routes = compatible_routes("multidomain", "miniprotein")
    expected_engines = {"RFD4-Proteina", "BoltzGen", "RFdiffusion3", "BindCraft2", "PXDesign", "Genie3", "Proteina-Complexa"}
    assert {route.engine for route in routes} == expected_engines
    items = catalogue(target_input="multidomain", binder_format="miniprotein")["items"]
    assert {item["route_id"] for item in items} == {route.id for route in routes}
    inputs = {"target_pdb": "fixture.cif", "target_chains": "A,B", "hotspots": "A7,B9"}
    graph = new_binder_pipeline("rfd4-internal", target_input="multidomain", binder_format="miniprotein", inputs=inputs)
    target = graph["nodes"][0]
    assert target["data"]["label"] == "Multidomain target"
    assert target["data"]["params"]["target_input"] == "multidomain"
    request = node_request("rfd4", _settings(graph, "rfd4"), {"target": target["data"]["params"]}, "Fixture")
    assert request["parameters"]["target_chains"] == ["A", "B"]
    assert request["parameters"]["hotspots"] == "A7,B9"
    assert request["parameters"]["task"] == "protein_binder"
    with pytest.raises(ValueError, match="does not support"):
        new_binder_pipeline("esmfold2-protein", target_input="multidomain", binder_format="miniprotein")
    with pytest.raises(ValueError, match="selected Target step"):
        new_binder_pipeline("esmfold2-protein", target_input="structure", binder_format="miniprotein",
                            graph=graph, target_node_id=target["id"])


def test_route_compatibility_and_branch_preserve_existing_target():
    assert [route.id for route in compatible_routes("idr", "cyclic_peptide")] == ["bindcraft2-idr-cyclic-peptide"]
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


@pytest.mark.parametrize("route_id, kind", [
    ("bindcraft2-idr-miniprotein", "bindcraft2"), ("peptide-trap-idr", "peptide_trap"), ("kejia-idr", "kejia_idr"),
])
def test_shared_idr_input_reaches_design_nodes(route_id, kind, tmp_path, monkeypatch):
    from mn_molecular_design_studio import auth

    monkeypatch.setattr(auth, "user_input_root", lambda owner: tmp_path)
    target = stage_idr_target("ACDEFGHIKLMN", "owner", hotspots="A3")
    graph = new_binder_pipeline(route_id, target_input="idr", binder_format="miniprotein", inputs=target)
    settings = _settings(graph, kind)
    if kind == "kejia_idr":
        settings["template_library"] = "walle"
    validate_graph(graph, require_inputs=True, require_files=False)
    target["target_fasta"] = write_idr_fasta(target["target_sequence"], tmp_path)
    request = node_request(kind, settings, {"target": target}, "Sequence pipeline")
    if kind == "bindcraft2":
        assert request["parameters"]["target_type"] == "sequence"
        assert request["parameters"]["target_fasta"] == target["target_fasta"]
        assert request["parameters"]["target_pdb"] is None
        assert request["parameters"]["target_chains"] == []
        assert request["parameters"]["hotspots"] == ""
    else:
        assert request["parameters"]["target_sequence"] == target["target_sequence"]


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
