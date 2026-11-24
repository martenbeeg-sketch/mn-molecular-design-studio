"""The visible analysis edge must also govern native campaign acceptance."""
from test_pipeline_graph import node, edge
import pytest
from mn_molecular_design_studio import pipeline_graph as graph


def test_canvas_preserves_typed_prediction_choices():
    definition = graph.new_graph("cofold_protein_protein")
    localfold = next(n for n in definition["nodes"] if n["data"]["kind"] == "prediction_localfold")
    assert localfold["data"]["params"]["af2_num_models"] == 1
    graph.clean_graph(definition)
    localfold["data"]["params"]["af2_num_models"] = True
    with pytest.raises(ValueError, match="Unsupported AF2 model number"):
        graph.clean_graph(definition)


@pytest.mark.parametrize("kind,source_kind,affinity", [("af2", "project_design_source", False),
    ("boltz2", "cofold_input", False), ("boltz2", "cofold_input", True)])
def test_prediction_capabilities_show_main_engines_and_open_connected_drafts(monkeypatch, kind, source_kind, affinity):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_capabilities_ui

    def click_engine(*, data, **_):
        assert data.count("scope='col'") == 2
        assert data.count('data-pipeline="true"') == data.count('type="button"') > 0
        if affinity:
            assert 'Pipeline Builder · Boltz-2 affinity' in data
            assert 'Pipeline Builder · AutoDock Vina' in data
        else:
            assert 'Pipeline Builder · AlphaFold2 Initial Guess' in data
            assert 'data-engine="cofold_input"' not in data
            assert 'data-engine="prediction_esmfold2"' in data
            assert 'data-engine="prediction_esmfold2_fast"' not in data
            assert 'data-engine="tile_af3_jax"' not in data
        return {"choose": [kind, None]}

    monkeypatch.setattr(design_capabilities_ui, "_matrix", click_engine)
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.design_capabilities_ui import render_capabilities
from mn_molecular_design_studio.pipeline_graph_ui import _open_capability_graph
if 'cap-graph' not in st.session_state:
    render_capabilities('cap', lambda draft: _open_capability_graph(draft, 'cap'))
""")
    app.session_state["cap-capability"] = "Docking & affinity" if affinity else "Structure prediction & cofolding"
    app.run()
    assert "cap-graph" not in app.session_state
    category = "Docking & affinity" if affinity else "Structure prediction & cofolding"
    setup_key = f"cap-operations-{category}-{kind}"
    if kind == "boltz2" and not affinity:
        app.selectbox(key=setup_key + "-molecular-mode").set_value("protein_ligand").run()
    app.button(key=setup_key + "-create-pipeline").click().run()
    assert not app.exception
    draft = app.session_state["cap-graph"]
    source, predictor, results = draft["nodes"]
    assert [node["data"]["kind"] for node in draft["nodes"]] == [source_kind, kind, "results"]
    assert [(edge["source"], edge["target"]) for edge in draft["edges"]] == [
        (source["id"], predictor["id"]), (predictor["id"], results["id"]),
    ]
    if kind == "boltz2":
        assert source["data"]["params"]["task_type"] == "protein_ligand"
        assert predictor["data"]["params"]["predict_affinity"] is affinity
    else:
        assert predictor["data"]["params"]["af2_validation_mode"] == "conditioned_self_consistency"


@pytest.mark.parametrize("engine,mode,kind", [("LDDM", "reference_ligand", "tile_generation_lddm"),
                                             ("OMTRA", "pharmacophore", "tile_generation_omtra")])
def test_small_molecule_modes_open_the_matching_pipeline_draft(monkeypatch, engine, mode, kind):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_capabilities_ui

    def click_mode(*, data, **_):
        assert 'data-engine="LDDM" data-output="reference_ligand"' in data
        assert 'data-engine="LDDM" data-output="pharmacophore"' not in data
        for omtra_mode in ("pocket", "reference_ligand", "pharmacophore"):
            assert f'data-engine="OMTRA" data-output="{omtra_mode}"' in data
        assert 'aria-label="Pipeline Builder · LDDM' in data
        assert data.count('data-pipeline="true"') == data.count('type="button"') > 0
        return {"choose": [engine, mode]}

    monkeypatch.setattr(design_capabilities_ui, "_matrix", click_mode)
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.design_capabilities_ui import render_capabilities
from mn_molecular_design_studio.pipeline_graph_ui import _open_capability_graph
if 'cap-graph' not in st.session_state:
    render_capabilities('cap', lambda draft: _open_capability_graph(draft, 'cap'))
""")
    app.session_state["cap-capability"] = "Small-molecule design"
    app.run()
    assert app.selectbox(key="cap-capability").options[0] == "Protein design"
    app.button(key=f"cap-{kind}-{mode}-create-pipeline").click().run()
    assert not app.exception
    draft = app.session_state["cap-graph"]
    nodes = {node["data"]["kind"]: node for node in draft["nodes"]}
    assert nodes[kind]["data"]["params"]["ligand_design_mode"] == mode
    assert ("target" in nodes) == (mode == "reference_ligand")
    assert "results" in nodes
    assert ("ligand_sdf_input" in nodes) == (mode == "reference_ligand")
    graph.validate_graph(draft)


def test_capabilities_distinguish_simulation_sampling_and_analysis_tasks(monkeypatch):
    import re

    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_capabilities_ui

    def render_table(*, data, **_):
        import streamlit as st

        st.text(data)
        return {}

    monkeypatch.setattr(design_capabilities_ui, "_matrix", render_table)
    app = AppTest.from_string("""
from mn_molecular_design_studio.design_capabilities_ui import render_capabilities
render_capabilities('cap', lambda _: None)
""")
    app.session_state["cap-capability"] = "Molecular dynamics"
    app.run()
    assert re.findall(r'data-engine="([^"]+)"', app.text[0].value) == ["tile_md_gromacs", "tile_md_openmm"]
    app.selectbox(key="cap-capability").set_value("Conformational sampling").run()
    assert re.findall(r'data-engine="([^"]+)"', app.text[0].value) == ["BioEmu-2"] * 3
    app.selectbox(key="cap-capability").set_value("Analysis & validation").run()
    assert app.selectbox(key="cap-analysis-task").value == "Protein design evaluation"
    assert 'data-engine="tile_binderscope"' in app.text[0].value
    assert 'data-engine="plip"' not in app.text[0].value
    app.selectbox(key="cap-analysis-task").set_value("Protein–ligand pose analysis").run()
    assert set(re.findall(r'data-engine="([^"]+)"', app.text[0].value)) == {"native_ligand_interactions", "plip", "pandamap", "posebusters"}
    app.selectbox(key="cap-analysis-task").set_value("MD trajectory analysis").run()
    assert 'data-engine="md_analysis"' in app.text[0].value
    assert "radius of gyration" in app.text[0].value
    app.selectbox(key="cap-analysis-task").set_value("Combine, filter & rank candidates").run()
    assert re.findall(r'data-engine="([^"]+)"', app.text[0].value) == [
        "candidate_merge", "candidate_analyzer", "filter", "shared_filter", "candidate_rank",
    ]
    assert "Pass, Fail, or Unresolved" in app.text[0].value
    assert "binder acceptance profile" in app.text[0].value
    assert not app.exception


@pytest.mark.parametrize("mode,project_designs", [
    ("sample-equilibrium", False), ("predict-structure", False), ("propagate", True), ("propagate", False),
])
def test_bioemu2_mode_cells_open_mode_specific_drafts(monkeypatch, mode, project_designs):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_capabilities_ui

    def click_mode(*, data, **_):
        assert all(f'data-output="{value}"' in data for value in ("sample-equilibrium", "predict-structure", "propagate"))
        return {"choose": ["BioEmu-2", mode]}

    monkeypatch.setattr(design_capabilities_ui, "_matrix", click_mode)
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.design_capabilities_ui import render_capabilities
from mn_molecular_design_studio.pipeline_graph_ui import _open_capability_graph
if 'cap-graph' not in st.session_state:
    render_capabilities('cap', lambda draft: _open_capability_graph(draft, 'cap'))
""")
    app.session_state["cap-capability"] = "Conformational sampling"
    app.run()
    if project_designs:
        app.selectbox(key=f"cap-bioemu2-{mode}-source").set_value("project_design_source").run()
    app.button(key=f"cap-bioemu2-{mode}-create-pipeline").click().run()
    assert not app.exception
    draft = app.session_state["cap-graph"]
    expected = ["project_design_source", "tile_bioemu2", "results"] if project_designs else ["tile_bioemu2", "results"]
    assert [node["data"]["kind"] for node in draft["nodes"]] == expected
    assert draft["nodes"][-2]["data"]["params"]["workflow_parameters"]["mode"] == mode
    graph.validate_graph(draft)


@pytest.mark.parametrize("analysis", [False, True])
def test_md_setup_builds_the_selected_dependencies(monkeypatch, analysis):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_capabilities_ui

    monkeypatch.setattr(design_capabilities_ui, "_matrix", lambda **_: {"choose": ["tile_md_openmm", None]})
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.design_capabilities_ui import render_capabilities
from mn_molecular_design_studio.pipeline_graph_ui import _open_capability_graph
if "cap-graph" not in st.session_state:
    render_capabilities("cap", lambda draft: _open_capability_graph(draft, "cap"))
""")
    app.session_state["cap-capability"] = "Molecular dynamics"
    app.run()
    key = "cap-operations-Molecular dynamics-tile_md_openmm"
    app.radio(key=key + "-prepare").set_value("unprepared").run()
    if analysis:
        app.checkbox(key=key + "-trajectory").check().run()
    assert "cap-graph" not in app.session_state
    app.button(key=key + "-create-pipeline").click().run()
    assert not app.exception
    draft = app.session_state["cap-graph"]
    nodes = {item["id"]: item["data"]["kind"] for item in draft["nodes"]}
    assert set(nodes.values()) == {"tile_complex_prepare", "tile_md_openmm", "results"} | (
        {"md_analysis"} if analysis else set()
    )
    links = {(nodes[edge["source"]], nodes[edge["target"]], edge["targetHandle"]) for edge in draft["edges"]}
    assert ("tile_complex_prepare", "tile_md_openmm", "complex") in links
    assert ("tile_md_openmm", "results", "in") in links
    if analysis:
        assert ("tile_md_openmm", "md_analysis", "in") in links
        assert ("md_analysis", "results", "in") in links
    graph.validate_graph(draft)


def test_prediction_pose_checks_are_parallel_and_keep_the_original_prediction():
    from mn_molecular_design_studio.capability_pipeline import add_followup, build_capability_graph

    draft = build_capability_graph("boltz2", prediction_mode="protein_ligand")
    draft = add_followup(draft, "posebusters")
    draft = add_followup(draft, "plip")
    nodes = {item["id"]: item["data"]["kind"] for item in draft["nodes"]}
    links = {(nodes[edge["source"]], nodes[edge["target"]]) for edge in draft["edges"]}
    assert links == {
        ("cofold_input", "boltz2"), ("boltz2", "results"), ("boltz2", "posebusters"),
        ("boltz2", "plip"), ("posebusters", "results"), ("plip", "results"),
    }
    graph.validate_graph(draft)


def test_saved_database_selection_is_frozen_only_when_creating_the_draft(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import database_selections, design_capabilities_ui, molecular_database

    class Database:
        def list_molecules(self, user, *, project_id, query, limit):
            assert user == {"id": "owner"}
            assert (project_id, query, limit) == ("project", "Example", 100)
            return {"items": [{"id": "saved-entry", "name": "Example entry",
                               "components": [{"type": "protein", "sequence": "GGGG"}]}]}

    snapshots = []

    def freeze(database, user, ids, project_id):
        snapshots.append((user["id"], ids, project_id))
        return tmp_path / "selection"

    monkeypatch.setattr(molecular_database, "MolecularDatabase", Database)
    monkeypatch.setattr(database_selections, "freeze_protein_selection", freeze)
    monkeypatch.setattr(design_capabilities_ui, "_matrix", lambda **_: {"choose": ["sequence_qc", None]})
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio.design_capabilities_ui import render_capabilities
from mn_molecular_design_studio.pipeline_graph_ui import _open_capability_graph
if "cap-graph" not in st.session_state:
    render_capabilities("cap", lambda draft: _open_capability_graph(draft, "cap"),
                        database_user={"id": "owner"}, project_id="project")
""")
    app.session_state["cap-capability"] = "Analysis & validation"
    app.session_state["cap-analysis-task"] = "Sequence & protein properties"
    app.run()
    key = "cap-operations-Analysis & validation-Sequence & protein properties-sequence_qc"
    app.selectbox(key=key + "-source").set_value("molecule_source").run()
    assert app.button(key=key + "-create-pipeline").disabled
    app.text_input(key=key + "-database-query").set_value("Example").run()
    app.multiselect(key=key + "-database-selection").set_value(["saved-entry"]).run()
    assert not snapshots
    app.button(key=key + "-create-pipeline").click().run()
    assert not app.exception
    assert snapshots == [("owner", ["saved-entry"], "project")]
    source = app.session_state["cap-graph"]["nodes"][0]["data"]
    assert source["kind"] == "molecule_source"
    assert source["params"]["source_run_dir"] == str(tmp_path / "selection")


def test_campaign_stop_waits_for_the_canvas_ipsae_step_before_counting_passes():
    target, generator = node("target", "target"), node("rfd4", "generator")
    generator["data"]["params"].update(task="protein_binder", num_designs=1000)
    analysis = node("ipsae", "structure-analysis")
    collection = node("candidate_merge", "collection")
    analysis["data"]["params"]["pae_cutoff"] = 8.0
    acceptance = node("shared_filter", "acceptance")
    acceptance["data"]["params"].update(binder_format="miniprotein", target_type="fasta_idr",
        rules=[{"field": "common.confidence.generator.complex.iptm", "operator": ">=", "value": .7}])
    stop, result = node("campaign_stop", "stop"), node("results", "result")
    definition = dict(schema_version=graph.SCHEMA, id="a" * 32, name="Saved canvas analysis before acceptance",
        workflow_mode="native_complete", nodes=[target,generator,collection,analysis,acceptance,stop,result],
        edges=[edge("target","generator"),edge("generator","collection"),edge("collection","structure-analysis"),
               edge("structure-analysis","acceptance"),edge("acceptance","stop"),edge("stop","result")],
        viewport={"x":0,"y":0,"zoom":1})
    with pytest.raises(ValueError, match="analysis node must not be skipped"):
        graph.validate_graph(definition)
