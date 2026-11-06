from __future__ import annotations

from mn_molecular_design_studio import run_results
from mn_molecular_design_studio.cluster import config


def test_peer_results_collect_inspects_remote_source_through_coordinator_api(tmp_path, monkeypatch):
    from mn_molecular_design_studio.cluster import transport
    from mn_molecular_design_studio import workspace

    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    config.save({"mode": "peer", "coordinator_url": "http://127.0.0.1:18766"})
    monkeypatch.setattr(workspace, "_studio_home", lambda: tmp_path / "studio")
    candidate = {
        "candidate_id": "rfd4-design-17",
        "complex_pdb": "artifacts/models/complex.pdb",
        "model_index": "boltz-model-2",
        "metrics": {"iptm": 0.81},
    }
    calls = []

    def request(path, payload=None, **kwargs):
        calls.append((path, payload, kwargs))
        if path == "/api/v1/cluster/jobs/remote-source-17/actions":
            assert payload == {"action": "inspect"}
            return {"command_id": "inspect-command-17"}
        if path == "/api/v1/cluster/commands/inspect-command-17":
            return {"status": "completed", "response": {"results": {"candidates": [candidate]}}}
        raise AssertionError(path)

    monkeypatch.setattr(transport, "request", request)
    run = {"id": "run-17"}
    output = {
        "id": "output-17", "reference": "output-17", "dispatch_error_count": 0,
        "sources": [{
            "cluster_job_id": "remote-source-17", "task_group": "design", "run_id": "source-17",
            "node_id": "remote-node", "local_reference": "design/source-17",
            "updated_at": "2026-11-06T10:00:00Z", "status": "completed",
            "availability": "online", "has_candidate_collection": True,
        }],
    }

    result = run_results.collect(run, output, {"id": "user-17"}, token="session-token")

    assert [call[0] for call in calls] == [
        "/api/v1/cluster/jobs/remote-source-17/actions",
        "/api/v1/cluster/commands/inspect-command-17",
    ]
    assert all(call[2]["token"] == "session-token" for call in calls)
    assert result["complete"] is True
    assert result["candidate_count"] == 1
    row = result["candidates"][0]
    assert row["candidate_id"] == "rfd4-design-17"
    assert row["model_index"] == "boltz-model-2"
    assert row["complex_pdb"] == candidate["complex_pdb"]
    assert row["metrics"] == {"iptm": 0.81}
    assert row["source_job"] == "cluster:remote-source-17"
    assert row["source_node_id"] == "remote-node"

