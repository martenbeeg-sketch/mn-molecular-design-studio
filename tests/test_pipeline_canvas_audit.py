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
from mn_molecular_design_studio.pipeline_graph_ui import _open_capability_step
if 'cap-graph' not in st.session_state:
    render_capabilities('cap', lambda _: None, lambda step, **options: _open_capability_step(step, 'cap', **options))
""")
    app.session_state["cap-capability"] = "Docking & affinity" if affinity else "Structure prediction & cofolding"
    app.run()
    assert any(item.label == "Pipeline inputs and outputs" for item in app.expander)
    if kind == "boltz2" and not affinity:
        app.selectbox(key="cap-prediction-input-boltz2").set_value("protein_ligand").run()
    app.button(key=f"cap-open-operation-{kind}").click().run()
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
from mn_molecular_design_studio.pipeline_graph_ui import _open_capability_step
if 'cap-graph' not in st.session_state:
    render_capabilities('cap', lambda _: None, lambda step, **options: _open_capability_step(step, 'cap', **options))
""")
    app.session_state["cap-capability"] = "Small-molecule design"
    app.run()
    assert app.selectbox(key="cap-capability").options[0] == "Protein design"
    app.button(key="cap-open-ligand-mode").click().run()
    assert not app.exception
    draft = app.session_state["cap-graph"]
    nodes = {node["data"]["kind"]: node for node in draft["nodes"]}
    assert nodes[kind]["data"]["params"]["ligand_design_mode"] == mode
    assert "target" in nodes and "results" in nodes
    assert ("ligand_sdf_input" in nodes) == (mode == "reference_ligand")
    graph.validate_graph(draft)


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
