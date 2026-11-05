"""Graph nodes write workflow artifacts but never auto-promote into Database."""
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio.collection_policy import _CURRENT
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner


def test_graph_parent_is_enqueued_under_no_database_collection_policy(monkeypatch, tmp_path):
    from mn_molecular_design_studio import pipeline_graph
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue

    monkeypatch.setattr(runner, "validate_graph", lambda graph, require_inputs: {"graph": graph})
    monkeypatch.setattr("mn_molecular_design_studio.candidate_pipeline.preflight", lambda *args: None)
    observed = []
    monkeypatch.setattr(workflow_queue, "enqueue_workflow_call",
                        lambda *args, **kwargs: observed.append(_CURRENT.get()) or tmp_path)
    assert runner.enqueue_graph_pipeline({}, pipeline_name="fixture", owner_id="owner", project_id="project") == tmp_path
    assert observed == ["none"]


@pytest.mark.parametrize("workflow", ["sequence.mpnn", "refolding.optimized_prediction", "analysis.biotite_interface"])
def test_graph_child_request_explicitly_disables_database_collection(monkeypatch, tmp_path, workflow):
    child = tmp_path / "child"
    child.mkdir()
    write_json(child / "metadata.json", {"status": "queued"})
    parent = tmp_path / "parent"
    parent.mkdir()
    write_json(parent / "metadata.json", {"studio_development_mode": True})
    write_json(parent / "input.json", {"params": {"pipeline_name": "fixture"}})
    observed = {}

    def submit(path):
        observed.update(json.loads(Path(path).read_text()))
        return {"job": {"run_dir": str(child)}}

    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit)
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: None)
    runner._submit({"workflow": workflow, "parameters": {}}, parent,
                   {"id": "child-node"}, owner_id="owner", project_id="project", gpu_device="none")
    assert observed["database_collection_policy"] == "none"
