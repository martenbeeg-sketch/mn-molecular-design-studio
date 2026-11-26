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
                or route.format in {"antibody", "tcr"} or route.id in {"rfd4-cyclic", "rfd4-generator"})
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
            assert json.loads(settings["options_json"])["refold"] == (route.id not in {"rfd4-generator", "rfd4-external"})


def test_design_capabilities_match_registered_workflows_and_picker_routes():
    from collections import Counter

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
    by_route = {(row["route_id"], row["target_input"]): row["native_verification"]
                for row in data["items"] if row["route_id"]}
    assert Counter(row["status"] for row in by_route.values()) == {
        "archived_downstream": 19, "archived_results": 7, "unverified": 56,
    }
    assert by_route["rfd4-idr", "idr"]["status"] == "archived_downstream"
    assert by_route["rfpeptides", "structure"]["status"] == "archived_results"
    assert by_route["rfd4-internal", "structure"]["status"] == "archived_downstream"
    assert by_route["rfd4-internal", "multidomain"]["status"] == "unverified"
    assert by_route["rfd4-generator", "structure"]["status"] == "archived_results"
    assert by_route["bindcraft2", "structure"]["status"] == "archived_results"
    assert by_route["bindcraft2", "multidomain"]["status"] == "unverified"
    assert all(row["status"] != "unverified" for row in by_route.values() if row["priority_review"])
    assert by_route["rfd4-generator", "structure"]["priority_review"] is True
    assert all("native_verification" not in row for row in data["items"] if not row["route_id"])


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


@pytest.mark.parametrize("refolder", ["boltz2", "rf3", "prediction_openfold3_ob0", "prediction_mol1"])
def test_rfd4_modular_presets_connect_only_the_requested_stages(refolder):
    arguments = {"target_input": "multidomain", "binder_format": "miniprotein"}
    generation = new_binder_pipeline("rfd4-generator", **arguments)
    assert [n["data"]["kind"] for n in generation["nodes"]] == ["target", "rfd4", "results"]
    modular = new_binder_pipeline("rfd4-external", **arguments, refold_engine=refolder,
                                  soluble_mpnn=True, interface_analysis=True)
    assert [n["data"]["kind"] for n in modular["nodes"]] == [
        "target", "rfd4", "ligandmpnn", refolder, "ipsae", "candidate_merge", "shared_filter", "results",
    ]
    assert _settings(modular, "ligandmpnn")["fix_interface_residues"] is True
    assert json.loads(_settings(modular, "rfd4")["options_json"])["refold"] is False
    native = new_binder_pipeline("rfd4-internal", **arguments, interface_analysis=True)
    assert [n["data"]["kind"] for n in native["nodes"]] == [
        "target", "rfd4", "ipsae", "candidate_merge", "shared_filter", "results",
    ]
    assert json.loads(_settings(native, "rfd4")["options_json"])["refold"] is True
    with pytest.raises(ValueError, match="refolding engine"):
        new_binder_pipeline("rfd4-generator", **arguments, refold_engine=refolder)


@pytest.mark.parametrize(("route_id", "kind", "needs_mpnn"), [
    ("rfd3", "rfd3", True),
    ("bindcraft2-generator", "bindcraft2_generator", False),
    ("boltzgen-generator", "boltzgen_generator", True),
    ("complexa-generator", "complexa_generator", False),
])
def test_generator_presets_connect_mol1_with_target_msa_and_protein_mpnn_when_needed(route_id, kind, needs_mpnn):
    graph = new_binder_pipeline(route_id, target_input="structure", binder_format="miniprotein",
                                refold_engine="prediction_mol1", interface_analysis=True)
    expected = ["target", kind, *(["ligandmpnn"] if needs_mpnn else []), "prediction_mol1", "ipsae",
                "candidate_merge", "shared_filter", "results"]
    assert [node["data"]["kind"] for node in graph["nodes"]] == expected
    assert _settings(graph, "prediction_mol1")["target_use_msa"] is True
    if needs_mpnn:
        settings = _settings(graph, "ligandmpnn")
        assert settings["sequence_stage"] == "initial"
        assert settings["backend"] == "shared_ligandmpnn"


@pytest.mark.parametrize("route", ROUTES, ids=lambda route: route.id)
def test_interface_constraints_are_inserted_only_for_supported_routes(route):
    from mn_molecular_design_studio.binder_pipeline import supported_interface_constraints

    arguments = {"target_input": route.target_input, "binder_format": route.format, "interface_constraints": True}
    if not supported_interface_constraints(route):
        with pytest.raises(ValueError, match="target annotations"):
            new_binder_pipeline(route.id, **arguments)
        return
    draft = new_binder_pipeline(route.id, **arguments)
    target, annotation, engine = draft["nodes"][:3]
    assert [node["data"]["kind"] for node in (target, annotation, engine)] == ["target", "binder_context", route.method]
    assert annotation["data"]["label"] == "Interface Constraints"
    assert [(edge["source"], edge["target"]) for edge in draft["edges"][:2]] == [
        (target["id"], annotation["id"]), (annotation["id"], engine["id"]),
    ]


def test_interface_editor_shows_effective_constraints_in_one_view_and_reports_overlap(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import target_editor

    pdb = tmp_path / "target.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY A   2       3.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    calls = []

    def viewer(structures, **kwargs):
        calls.append((structures[0].highlighted_selections, kwargs["selection_mode"]))

    monkeypatch.setattr(target_editor, "molstar_custom_component", viewer)
    app = AppTest.from_string(f"""
from pathlib import Path
from mn_molecular_design_studio.target_editor import TargetSource, render_target_editor
render_target_editor(TargetSource(path=Path({str(pdb)!r})), key_prefix="constraints",
                     allow_trimming=False, inherited_interface={{"hotspots": "A1", "coldspots": "A2"}})
""").run()
    assert not app.exception
    assert calls == [([], True)]
    app.radio(key="constraints-annotation-view").set_value("both").run()
    assert calls[-1] == (["A1", "A2"], False)
    coldspots = next(field for field in app.text_input if field.label == "Chosen coldspot residues")
    coldspots.set_value("A1").run()
    assert "Residues in both lists: A1" in app.error[0].value
    assert calls[-1] == (["A1"], False)
    app.radio(key="constraints-annotation-view").set_value("coldspots").run()
    assert calls[-1] == (["A1"], True)


@pytest.mark.parametrize("prefix", [None, "project-editor"])
def test_capability_cell_opens_rfd4_preset_controls(monkeypatch, prefix):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_capabilities_ui

    def click_cell(**kwargs):
        assert 'aria-label="Build pipeline · RFD4-Proteina · Miniprotein"' in kwargs["data"]
        return {"choose": ["RFD4-Proteina", "miniprotein"]}

    monkeypatch.setattr(design_capabilities_ui, "_matrix", click_cell)

    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.design_capabilities_ui import render_design_tasks
from mn_molecular_design_studio.pipeline_graph_ui import open_capability_graph
prefix = st.session_state["test-prefix"]
pending = st.session_state.get("pipeline-pending-draft", {})
graph = st.session_state.get(f"{prefix}-graph") or pending.get("graph")
if graph:
    st.json([node["data"]["kind"] for node in graph["nodes"]])
else:
    render_design_tasks("cap", lambda draft: open_capability_graph(draft, prefix))
""")
    app.session_state["cap-target"] = "multidomain"
    app.session_state["test-prefix"] = prefix
    app.run()
    key = "cap-multidomain-RFD4-Proteina-miniprotein"
    selector = app.selectbox(key=key + "-variant")
    assert selector.value.id == "rfd4-generator"
    external = next(route for route in compatible_routes("multidomain", "miniprotein") if route.id == "rfd4-external")
    selector.set_value(external).run()
    app.selectbox(key=key + "-refold").set_value("rf3").run()
    app.checkbox(key=key + "-constraints").check().run()
    assert any("Coldspots are unavailable" in caption.value for caption in app.caption)
    app.checkbox(key=key + "-analysis").check().run()
    app.button(key=key + "-create-pipeline").click().run()
    assert not app.exception
    assert json.loads(app.json[0].value) == ["target", "binder_context", "rfd4", "rf3", "ipsae", "candidate_merge", "shared_filter", "results"]
    if prefix is None:
        assert app.session_state["page"] == "Pipelines"
        assert app.session_state["pipeline-pending-draft"]["project_id"] is None
    else:
        assert app.session_state[f"{prefix}-view"] == "editor"


def test_other_design_tasks_show_only_engines_for_the_selected_task(monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_capabilities_ui

    def render_table(**kwargs):
        import streamlit as st

        st.text(kwargs["data"])
        return {}

    monkeypatch.setattr(design_capabilities_ui, "_matrix", render_table)
    app = AppTest.from_string("""
from mn_molecular_design_studio.design_capabilities_ui import render_design_tasks
render_design_tasks("cap", lambda _: None)
""").run()
    assert app.selectbox(key="cap-task").value == "protein_binders"
    assert "Small-molecule generation" not in app.selectbox(key="cap-task").options
    app.selectbox(key="cap-task").set_value("new_fold").run()
    assert not app.exception
    assert 'data-engine="Genie3" data-output="new_fold"' in app.text[0].value
    assert 'data-engine="AbMPNN"' not in app.text[0].value
    assert app.text[0].value.count('scope="col"') == 2


def test_protein_binder_verification_filter_keeps_only_matching_variants(monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_capabilities_ui

    def click_filtered_cell(**kwargs):
        markup = kwargs["data"]
        if kwargs["key"].endswith("-all"):
            assert "Varies by workflow" in markup
            return {}
        if kwargs["key"].endswith("-priority"):
            assert "RFpeptides" not in markup
            assert "Miniprotein" in markup
            return {}
        assert kwargs["key"].endswith("-archived_results")
        assert "Archived Results verified" in markup
        assert "Archived downstream verified" not in markup
        return {"choose": ["BindCraft2", "miniprotein"]}

    monkeypatch.setattr(design_capabilities_ui, "_matrix", click_filtered_cell)
    app = AppTest.from_string("""
from mn_molecular_design_studio.design_capabilities_ui import render_design_tasks
render_design_tasks("cap", lambda _: None)
""").run()
    app.selectbox(key="cap-verification").set_value("priority").run()
    assert not app.exception
    app.selectbox(key="cap-verification").set_value("archived_results").run()
    assert not app.exception
    assert all(widget.label != "Workflow goal" for widget in app.selectbox)
    assert any("Archived Results verified" in item.value and "Priority review" in item.value for item in app.caption)


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
