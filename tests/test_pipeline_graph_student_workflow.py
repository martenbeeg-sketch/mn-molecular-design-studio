from mn_molecular_design_studio.pipeline_graph import new_graph, validate_graph


def test_rfd4_template_is_an_editable_complete_student_analysis_route():
    graph = new_graph("rfd4_protein_binder")
    validate_graph(graph)
    nodes = {node["data"]["kind"]: node["data"] for node in graph["nodes"]}
    route = [node["data"]["kind"] for node in graph["nodes"]]

    assert route == [
        "target", "binder_context", "rfd4", "ligandmpnn", "boltz2", "ipsae",
        "rosetta_metrics", "candidate_merge", "shared_filter", "complex_analysis",
        "biotite_interface", "results",
    ]
    assert nodes["ligandmpnn"]["label"] == "Optional SolubleMPNN optimization"
    assert nodes["ligandmpnn"]["params"]["sequence_stage"] == "soluble_optimization"
    assert nodes["ligandmpnn"]["params"]["optimization_enabled"] is False
    assert nodes["rfd4"]["params"]["task"] == "protein_binder"
    assert nodes["complex_analysis"]["params"]["entry_mode"] == "accepted_only"


def test_backbone_only_generators_have_a_separate_initial_sequence_design_route():
    graph = new_graph("all_generators_boltz2")
    initial = next(node for node in graph["nodes"] if node["data"]["kind"] == "ligandmpnn")
    soluble = next(node for node in graph["nodes"] if node["id"] == "optional-soluble-mpnn")

    assert all(node["data"]["kind"] != "rfd4" for node in graph["nodes"])
    assert initial["data"]["params"]["sequence_stage"] == "initial"
    assert soluble["data"]["params"]["sequence_stage"] == "soluble_optimization"
