"""Exercise every advertised setup without submitting scientific jobs."""

import json
from html import escape

import pytest

from mn_molecular_design_studio import design_capabilities_ui as ui
from mn_molecular_design_studio.binder_pipeline import new_binder_pipeline
from mn_molecular_design_studio.capability_pipeline import (
    MD_ENGINES,
    SAVED_DESIGN_REFOLDERS,
    add_followup,
    build_capability_graph,
    followup_choices,
    source_choices,
)
from mn_molecular_design_studio.design_capabilities import catalogue, small_molecule_catalogue
from mn_molecular_design_studio.engine.workflow_engine.workflows.bioemu2 import BIOEMU2_MODES
from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, list_definitions, save_definition, validate_graph
from mn_molecular_design_studio.scientific_tasks import PREDICTION_MODES, prediction_methods
from mn_molecular_design_studio.tool_connections import editor_tool_catalog


def binder_cases():
    for row in catalogue()["items"]:
        if row["access"] != "pipeline":
            continue
        label = f'{row["engine"]}/{row["target_input"]}/{row["output"]}/{row["route_id"]}'
        choices = [{}]
        if row["route_id"] in ui.SEPARATE_REFOLD_ROUTES:
            choices.extend({"refold_engine": engine} for engine in ui.RFD4_REFOLD_ENGINES)
        if row["route_id"] in ui.SEPARATE_REFOLD_ROUTES | {"rfd4-internal"}:
            choices.append({"interface_analysis": True})
        if row["route_id"] in ui.RFD4_MODULAR_ROUTES:
            choices.append({"soluble_mpnn": True})
        for options in choices:
            suffix = "/".join(f"{key}={value}" for key, value in options.items()) or "default"
            yield pytest.param(row, options, id=f"{label}/{suffix}")


def operation_cases():
    for row in catalogue()["items"]:
        task = row["output"]
        if row["access"] == "canvas" and task != "ligand_design":
            yield pytest.param(row["canvas_node"], {"design_task": task}, id=f'{row["engine"]}/{task}')
    for row in small_molecule_catalogue()["items"]:
        for mode in row["modes"]:
            yield pytest.param(row["operation"], {"ligand_design_mode": mode}, id=f'{row["engine"]}/{mode}')
    for mode in BIOEMU2_MODES:
        yield pytest.param("tile_bioemu2", {"bioemu2_mode": mode}, id=f"bioemu2/{mode}")
    yield pytest.param("tile_bioemu2", {"bioemu2_mode": "propagate", "source": "project_design_source"},
                       id="bioemu2/propagate/saved-design")
    operations = editor_tool_catalog()
    for category, categories in ui._OPERATION_CATEGORIES.items():
        if category == ui._SAMPLING_CATEGORY:
            continue
        for kind, spec in operations.items():
            if spec["category"] not in categories or spec.get("hidden"):
                continue
            if category == ui._PREDICTION_CATEGORY:
                if kind not in ui._PREDICTION_ENGINES:
                    continue
                modes = [mode for mode in PREDICTION_MODES if kind in prediction_methods(mode)]
                for mode in modes:
                    yield pytest.param(kind, {"prediction_mode": mode}, id=f"{kind}/{mode}")
                if modes and kind not in SAVED_DESIGN_REFOLDERS:
                    continue
            sources = source_choices(kind, database_available=True)
            if not sources or kind == "tile_tmol":
                yield pytest.param(kind, {}, id=f"{kind}/default")
            for source in sources:
                yield pytest.param(kind, {"source": source}, id=f"{kind}/{source}")
            if kind in (*MD_ENGINES, "tile_library_centroids"):
                yield pytest.param(kind, {"prepare": True}, id=f"{kind}/prepare-first")
    yield pytest.param("boltz2", {"prediction_mode": "protein_ligand", "predict_affinity": True}, id="boltz2/affinity")


def assert_editable(graph, operation):
    saved = json.loads(json.dumps(graph))
    validated = validate_graph(saved)["graph"]
    assert validated["nodes"] == graph["nodes"]
    assert validated["edges"] == graph["edges"]
    kinds = {node["data"]["kind"] for node in validated["nodes"]}
    assert {operation, "results"} <= kinds


def comparable(graph):
    indices = {node["id"]: index for index, node in enumerate(graph["nodes"])}
    nodes = [node["data"] for node in graph["nodes"]]
    edges = sorted(
        (indices[edge["source"]], edge["sourceHandle"], indices[edge["target"]], edge["targetHandle"])
        for edge in graph["edges"]
    )
    return nodes, edges


@pytest.mark.parametrize(("row", "options"), tuple(binder_cases()))
def test_binder_preset_contract(row, options):
    graph = new_binder_pipeline(row["route_id"], target_input=row["target_input"], binder_format=row["output"], **options)
    assert_editable(graph, row["canvas_node"])
    kinds = {node["data"]["kind"] for node in graph["nodes"]}
    if refolder := options.get("refold_engine"):
        assert refolder in kinds
    if options.get("interface_analysis"):
        assert "ipsae" in kinds
    if options.get("soluble_mpnn"):
        assert any(
            node["data"]["kind"] == "ligandmpnn" and node["data"]["params"]["model_type"] == "soluble_mpnn"
            for node in graph["nodes"]
        )


@pytest.mark.parametrize(("row", "options"), tuple(binder_cases()))
def test_binder_cell_creates_and_reopens_its_own_preset(row, options, monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest

    engine, output, context = (row[key] for key in ("engine", "output", "target_input"))
    expected_button = f'data-engine="{escape(engine)}" data-output="{escape(output)}"'

    def choose_cell(**kwargs):
        assert expected_button in kwargs["data"]
        return {"choose": [engine, output]}

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    monkeypatch.setattr(ui, "_matrix", choose_cell)
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.design_capabilities_ui import render_design_tasks
render_design_tasks("coverage", lambda graph: st.session_state.__setitem__("coverage-graph", graph))
""")
    app.session_state["coverage-target"] = context
    app.run()
    assert not app.exception
    prefix = f"coverage-{context}-{engine}-{output}"
    variants = [widget for widget in app.selectbox if widget.key == prefix + "-variant"]
    if variants:
        from mn_molecular_design_studio.binder_pipeline import compatible_routes

        route = next(route for route in compatible_routes(context, output) if route.id == row["route_id"])
        variants[0].set_value(route).run()
    if refolder := options.get("refold_engine"):
        app.selectbox(key=prefix + "-refold").set_value(refolder).run()
    if options.get("interface_analysis"):
        app.checkbox(key=prefix + "-analysis").check().run()
    if options.get("soluble_mpnn"):
        app.checkbox(key=prefix + "-soluble").check().run()
    app.button(key=prefix + "-create-pipeline").click().run()
    assert not app.exception and not app.error
    graph = app.session_state["coverage-graph"]
    assert_editable(graph, row["canvas_node"])
    refolders = [widget for widget in app.selectbox if widget.key == prefix + "-refold"]
    settings = dict(options)
    if refolders:
        settings["refold_engine"] = refolders[0].value
    expected = new_binder_pipeline(row["route_id"], target_input=context, binder_format=output, **settings)
    assert comparable(graph) == comparable(expected)
    saved = save_definition(graph, "coverage-user", "coverage-project")
    reopened, = list_definitions("coverage-user", "coverage-project")
    assert reopened["revision"] == saved["revision"] and reopened["graph"] == graph


@pytest.mark.parametrize(("operation", "options"), tuple(operation_cases()))
def test_operation_preset_contract(operation, options):
    graph = build_capability_graph(operation, **options)
    assert_editable(graph, operation)
    by_kind = {node["data"]["kind"]: node["data"]["params"] for node in graph["nodes"]}
    if mode := options.get("prediction_mode"):
        assert by_kind["cofold_input"]["task_type"] == mode
    if mode := options.get("ligand_design_mode"):
        assert by_kind[operation]["ligand_design_mode"] == mode
    if mode := options.get("bioemu2_mode"):
        assert by_kind[operation]["workflow_parameters"]["mode"] == mode
    if source := options.get("source"):
        assert source in by_kind


@pytest.mark.parametrize(("operation", "options"), tuple(operation_cases()))
def test_operation_setup_saves_and_reopens_preset(operation, options, monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    setup = {k: v for k, v in options.items() if k in
             {"predict_affinity", "ligand_design_mode", "bioemu2_mode", "design_task"}}
    setup["prediction"] = bool(options.get("prediction_mode") and not options.get("predict_affinity")) or (
        operation in SAVED_DESIGN_REFOLDERS and not options.get("predict_affinity"))
    app = AppTest.from_string(f"""
import streamlit as st
from mn_molecular_design_studio.design_capabilities_ui import _operation_setup
_operation_setup({operation!r}, "setup", lambda graph: st.session_state.__setitem__("graph", graph),
                 database_user={{"id": "coverage-user", "role": "user"}}, project_id="coverage-project", **{setup!r})
""")
    if options.get("source") == "molecule_source":
        # Selection, access and provenance are exercised through the live database-input check.
        def ready_source(*args):
            def configure(graph):
                node = next(node for node in graph["nodes"] if node["data"]["kind"] == "molecule_source")
                node["data"]["params"]["source_run_dir"] = "runs:///database-selection/test-snapshot"
            return configure
        monkeypatch.setattr(ui, "_database_source", ready_source)

    app.run()
    key = "setup-prediction-source"
    if not options.get("prediction_mode") and any(control.key == key for control in app.radio):
        app.radio(key=key).set_value("Saved designs").run()
    if options.get("prediction_mode") and not options.get("predict_affinity"):
        app.selectbox(key="setup-molecular-mode").set_value(options["prediction_mode"]).run()
    if options.get("source") and any(control.key == "setup-source" for control in app.selectbox):
        app.selectbox(key="setup-source").set_value(options["source"]).run()
    if options.get("prepare"):
        app.radio(key="setup-prepare").set_value("unprepared").run()
    if options.get("design_task") == "new_fold" and operation != "rfd4":
        app.selectbox(key="setup-goal").set_value("backbones").run()
    app.button(key="setup-create-pipeline").click().run()
    assert not app.exception and not app.error
    graph = app.session_state["graph"]
    expected = build_capability_graph(operation, **options)
    if options.get("source") == "molecule_source":
        ready_source()(expected)
    assert comparable(graph) == comparable(expected)
    saved = save_definition(graph, "coverage-user", "coverage-project")
    reopened, = list_definitions("coverage-user", "coverage-project")
    assert reopened["revision"] == saved["revision"] and reopened["graph"] == graph


@pytest.mark.parametrize(("start", "options", "next_steps"), [
    ("prediction_esmfold2", {"prediction_mode": "protein"}, ("sequence_qc", "proteinsol")),
    ("boltz2", {"prediction_mode": "protein_ligand", "predict_affinity": True}, ("posebusters", "native_ligand_interactions")),
    ("tile_complex_prepare", {}, ("tile_md_openmm", "md_analysis")),
    ("tile_complex_prepare", {}, ("tile_md_gromacs", "md_analysis")),
    ("tile_prepare_library", {}, ("tile_library_centroids",)),
])
def test_guided_optional_steps_assemble_connected_graphs(start, options, next_steps):
    graph = build_capability_graph(start, **options)
    for step in next_steps:
        assert step in followup_choices(graph, (step,))
        graph = add_followup(graph, step, handoff=step in MD_ENGINES)
        assert_editable(graph, step)
    with pytest.raises(ValueError, match="cannot consume"):
        add_followup(graph, "target")


def test_monomer_label_in_prediction_input_editor(monkeypatch):
    from mn_molecular_design_studio import pipeline_node_ui

    def selectbox(label, choices, *, index, format_func, **kwargs):
        assert label == "Prediction task"
        assert format_func(choices[index]) == "Protein monomer"
        return choices[index]

    monkeypatch.setattr(pipeline_node_ui.st, "selectbox", selectbox)
    field = next(field for field in NODE_CATALOG["cofold_input"]["fields"] if field["name"] == "task_type")
    assert pipeline_node_ui._render_field(field, {"task_type": "protein"}, "test") == "protein"


def test_library_editor_accepts_paths_and_keeps_centroids_in_a_separate_step(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest

    reference = tmp_path / "reference"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "workspace"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR", str(reference))
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.pipeline_node_ui import _render_tile_workflow_settings
from mn_molecular_design_studio.pipeline_tile_workflows import TILES
node = {"id": "prep", "data": {"kind": "tile_prepare_library", "label": "Prepare library"}}
values, ready = _render_tile_workflow_settings(
    TILES["tile_prepare_library"], {}, node, {"nodes": [node], "edges": []}, "library",
)
st.session_state["values"] = values
""").run()
    assert not app.exception
    assert "external reference directory" in app.caption[0].value
    assert not any(control.label == "Create centroids" for control in app.checkbox)
    paths = [str(reference / "compounds.csv"), str(reference / "compounds.sdf")]
    app.text_area(key="library-source_files").set_value("\n\n".join(paths)).run()
    assert app.session_state["values"]["source_files"] == paths
