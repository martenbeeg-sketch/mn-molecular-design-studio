from mn_molecular_design_studio import pipeline_graph
from mn_molecular_design_studio.tool_connections import connection_catalog, editor_connections


def test_registered_nodes_expose_one_uniform_contract():
    required = {"operation", "category", "inputs", "outputs", "parameters", "constraints"}
    for kind, spec in pipeline_graph.NODE_CATALOG.items():
        names = [field["name"] for field in spec["fields"]]
        assert len(names) == len(set(names)), kind
        contract = pipeline_graph.capability_contract(kind, spec)
        assert required <= set(contract)
        assert len(contract["inputs"]) == len(spec.get("inputs") or [])
        assert len(contract["outputs"]) == len(spec.get("outputs") or [])
        assert all("artifact_types" in port and "required" in port and "multiple" in port
                   for port in contract["inputs"])
        assert all("artifact_types" in port and "metrics" in port for port in contract["outputs"])


def test_new_registered_metric_tool_automatically_connects_and_populates_filter(monkeypatch):
    # A hypothetical registration proves that compatibility and Filter metric
    # discovery consume the registry contract without Filter-specific edits.
    spec = {
        "key": "minimal_test_tool", "label": "Minimal test tool", "category": "Analysis",
        "description": "A temporary declared operation used to test modularity.",
        "inputs": [{"id": "in", "label": "Prediction", "accepts": ["prediction"], "multiple": False}],
        "outputs": [{"id": "out", "label": "Prediction with score", "type": "prediction"}],
        "fields": [{"name": "cutoff", "default": 0.0, "min": -1, "max": 1}],
        "formats": ["pdb"], "workflow": "test.minimal_metric",
        "output_metrics": [{"field": "test_tool.score", "label": "Test score", "unit": "fraction",
                            "method": "minimal deterministic test tool", "direction": "higher"}],
    }
    monkeypatch.setitem(pipeline_graph.NODE_CATALOG, "minimal_test_tool", spec)
    assert pipeline_graph.connection_compatibility("minimal_test_tool", "filter")["compatible"]
    graph = {"nodes": [
        {"id": "predictor", "data": {"kind": "boltz2"}},
        {"id": "test", "data": {"kind": "minimal_test_tool"}},
        {"id": "filter", "data": {"kind": "filter"}},
    ], "edges": [
        {"source": "predictor", "target": "test"},
        {"source": "test", "target": "filter"},
    ]}
    declared = pipeline_graph.filter_metric_options(graph, "filter")
    assert any(item["field"] == "test_tool.score" and item["availability"] == "declared" for item in declared)
    observed = pipeline_graph.filter_metric_options(graph, "filter", [
        {"candidate_id": "observed", "source_tool": "minimal_test_tool", "stage": "analysis",
         "metrics": {"runtime_metric": 0.5}, "raw_metadata": {}}
    ])
    assert any(item["field"] == "runtime_metric" and item["availability"] == "observed" for item in observed)



def test_localfold_canvas_form_preserves_numeric_choices():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
from mn_molecular_design_studio.pipeline_node_ui import render_node_form

spec = NODE_CATALOG["prediction_localfold"]
node = {"id": "step", "data": {"kind": "prediction_localfold", "label": spec["label"], "params": {}}}
st.session_state["settings"] = render_node_form(
    node, {"nodes": [node], "edges": []}, key="test", render_input=lambda: {}, owner_id="qa"
)[0]
""").run()
    assert not app.exception
    assert app.selectbox(key="test-af2_num_models").value == 1
    app.selectbox(key="test-af2_num_models").select(3).run()
    assert not app.exception
    assert app.session_state["settings"]["af2_num_models"] == 3


def test_connection_catalog_matches_canvas_for_every_registered_port_pair():
    page = connection_catalog(view="connections", query="", operation="", direction="both",
                              artifact_type="", offset=0, limit=10000)
    advertised = {(row["source"], row["output"], row["target"], row["input"]): row["artifact_type"]
                  for row in page["items"]}
    assert len(advertised) == page["total"]
    for source_kind, source in pipeline_graph.NODE_CATALOG.items():
        for output in source["outputs"]:
            for target_kind, target in pipeline_graph.NODE_CATALOG.items():
                for input_port in target["inputs"]:
                    key = (source_kind, output["id"], target_kind, input_port["id"])
                    result = pipeline_graph.connection_compatibility(
                        source_kind, target_kind, source_handle=output["id"], target_handle=input_port["id"],
                    )
                    assert advertised.get(key) == result.get("artifact_type"), key
                    assert (key in advertised) == result["compatible"], key


def test_active_learning_catalogue_only_offers_connections_in_the_template():
    graph = pipeline_graph.new_graph("unidock_active_learning")
    pipeline_graph.validate_graph(graph)
    kinds = {node["id"]: node["data"]["kind"] for node in graph["nodes"]}
    expected = {(kinds[edge["source"]], kinds[edge["target"]]) for edge in graph["edges"]}
    page = connection_catalog(view="connections", query="", operation="", direction="both",
                              artifact_type="active_learning_state", offset=0, limit=100)
    assert {(row["source"], row["target"]) for row in page["items"]} == expected


def test_editor_and_mcp_use_the_same_outgoing_port_matches():
    page = connection_catalog(view="connections", query="", operation="filter", direction="outgoing",
                              artifact_type="", offset=0, limit=10000, port="fail")
    index = editor_connections()
    assert page["items"]
    assert all(row["output"] == "fail" for row in page["items"])
    rebuilt = {(row["target"], row["input"]): row["artifact_type"] for row in page["items"]}
    assert rebuilt == {(target, port): kind for target, match in index["wires"]["filter"]["fail"].items()
                       for port, kind in match.items()}
    assert index["wires"]["tile_prepare_library"]["out"]["tile_library_centroids"]["features"] == "fingerprint_manifest"
    assert any(row["target"] == "library.seed_dock" for row in index["handoffs"]["library.centroids"])


def test_canvas_groups_round_trip_without_modifying_the_execution_plan():
    graph = pipeline_graph.new_graph("md_simulation")
    baseline = pipeline_graph.validate_graph(graph)
    groups = [{"id": "group-preparation", "name": "Prepare and simulate", "collapsed": True,
               "members": [node["id"] for node in graph["nodes"][:2]]}]
    grouped = pipeline_graph.validate_graph(dict(graph, layout_groups=groups))
    assert grouped["graph"]["layout_groups"] == groups
    assert grouped["graph"]["nodes"] == baseline["graph"]["nodes"]
    assert grouped["graph"]["edges"] == baseline["graph"]["edges"]
    assert grouped["order"] == baseline["order"]
