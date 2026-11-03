"""The visible analysis edge must also govern native campaign acceptance."""
from test_pipeline_graph import node, edge
import pytest
from mn_molecular_design_studio import pipeline_graph as graph


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
