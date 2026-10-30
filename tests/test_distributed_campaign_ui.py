from __future__ import annotations

import json
from pathlib import Path

from mn_molecular_design_studio.distributed_campaign_ui import campaign_context_for_job
from mn_molecular_design_studio.cluster import config as cluster_config
from mn_molecular_design_studio.cluster import ui as cluster_ui


def _campaign_fixture(root: Path, *, with_index: bool, validation_status: str = "completed", error: str = "") -> Path:
    campaign_id = "follow-up-campaign"
    source_id = "original-campaign"
    candidate_id = "dc_00000001_example"
    (root / "artifacts/distributed_campaign").mkdir(parents=True)
    (root / "artifacts/normalized_candidates").mkdir(parents=True)
    (root / "input.json").write_text(json.dumps({"params": {"batch_size": 4}}), encoding="utf-8")
    config = {
        "acceptance_level": "validated",
        "accepted_target": 2,
        "max_attempts": 3,
        "batch_size": 1,
        "max_concurrent_batches": 2,
        "max_concurrent_validations": 2,
        "validation_stages": [{"required_metrics": ["iptm", "ptm", "plddt_mean"]}],
    }
    stage = {"status": validation_status, "run_dir": "validation-run"}
    if validation_status == "completed":
        stage["evidence_candidate"] = {"metrics": {"iptm": 0.64, "ptm": 0.78, "plddt_mean": 84.5}}
    if error:
        stage["error"] = error
    ledger = {
        "campaign_id": campaign_id,
        "status": "completed" if validation_status == "completed" else "failed",
        "termination_reason": "target_reached",
        "config": config,
        "attempts": {
            "attempt-1": {
                "attempt_id": "attempt-1",
                "ordinal": 1,
                "batch_id": "retained-candidate-validation",
                "candidate_ids": [candidate_id],
                "status": "completed",
                "native_state": {
                    "attempt_unit": "one retained candidate; no new generation",
                    "source_campaign_id": source_id,
                    "source_worker": "worker-uuid",
                    "source_worker_run_id": "source-run-1",
                },
            }
        },
        "candidates": {
            candidate_id: {
                "attempt_id": "attempt-1",
                "candidate_id": candidate_id,
                "candidate": {
                    "candidate_id": candidate_id,
                    "raw_metadata": {
                        "distributed_campaign": {
                            "campaign_id": source_id,
                            "worker": "worker-uuid",
                            "worker_run_id": "source-run-1",
                        }
                    },
                },
                "validated_acceptance": "accepted" if validation_status == "completed" else "failed",
                "validation_decision_reason": error,
                "validation_stages": {"0": stage},
            }
        },
    }
    (root / "artifacts/distributed_campaign/ledger.json").write_text(json.dumps(ledger), encoding="utf-8")
    if with_index:
        row = {"candidate_id": candidate_id, "stage": "generation.backbone_sequence"}
        (root / "artifacts/normalized_candidates/candidates.jsonl").write_text(json.dumps(row) + "\n", encoding="utf-8")
    return root


def test_follow_up_projection_separates_dispatch_inference_measurements_and_acceptance(tmp_path):
    root = _campaign_fixture(tmp_path / "followup", with_index=True)

    context = campaign_context_for_job(root)

    assert context["campaign_kind"] == "Validation-only follow-up"
    assert context["source_campaign_ids"] == ["original-campaign"]
    assert context["controller_batch_size"] == 1
    assert context["controller_dispatch_capacity"] == 2
    assert context["max_concurrent_validations"] == 2
    assert context["accepted_target"] == 2
    assert context["attempt_budget"] == 3
    assert context["attempts_authorized"] == 1
    assert context["attempts_completed"] == 1
    assert context["active_attempts"] == 0
    assert context["native_inference_batch_size"] == 4
    candidate = context["candidates"][0]
    assert candidate["candidate_index_status"] == "Row 1"
    assert candidate["acceptance_label"] == "Accepted"
    assert candidate["validation_stages"][0]["measurement_status"] == "Valid measurements"
    assert context["accepted_count"] == 1
    assert context["valid_measurement_count"] == 1
    assert candidate["worker"] == "worker-uuid"
    assert candidate["worker_run_id"] == "source-run-1"
    assert candidate["attempt_id"] == "attempt-1"


def test_campaign_context_projects_deferred_worker_batches_as_counts_and_references(tmp_path):
    root = _campaign_fixture(tmp_path / "deferred", with_index=False)
    ledger_path = root / "artifacts/distributed_campaign/ledger.json"
    ledger = json.loads(ledger_path.read_text())
    ledger["config"]["remote_artifact_policy"] = "defer_generated"
    ledger["attempts"] = {
        "attempt-remote": {
            "attempt_id": "attempt-remote", "ordinal": 1, "status": "completed",
            "deferred_candidate_count": 4, "native_state": {"output_storage": "worker_local"},
        }
    }
    ledger["candidates"] = {}
    ledger["deferred_outputs"] = {
        "receipt": {
            "receipt_id": "receipt", "attempt_id": "attempt-remote", "candidate_count": 4,
            "candidates": [], "status": "deferred",
            "evidence": {
                "source_job": "cluster://worker/design/attempt", "worker_node_id": "worker",
                "worker_run_id": "design/attempt", "native_filter_counts": {"passed": 2, "failed": 1, "unavailable": 1},
            },
        }
    }
    ledger_path.write_text(json.dumps(ledger))

    context = campaign_context_for_job(root)

    assert context["candidate_count"] == 4
    assert context["deferred_candidate_count"] == 4
    assert len(context["deferred_batches"]) == 1
    assert context["deferred_batches"][0]["candidate_count"] == 4
    assert context["deferred_batches"][0]["source_job"] == "cluster://worker/design/attempt"
    assert context["deferred_candidates"] == []


def test_missing_candidate_index_and_lineage_are_unresolved_not_scientific_rejection(tmp_path):
    root = _campaign_fixture(
        tmp_path / "missing-index",
        with_index=False,
        validation_status="failed",
        error="returned 0 same-sequence descendants; normalized candidate index was not available",
    )

    context = campaign_context_for_job(root)

    candidate = context["candidates"][0]
    assert candidate["candidate_index_status"] == "Missing index file"
    assert candidate["acceptance_label"] == "Unresolved · candidate index/lineage metadata"
    assert candidate["lineage_error"] is True
    assert candidate["metadata_issue"] is True
    assert context["rejected_count"] == 0
    assert context["unresolved_count"] == 1
    assert context["valid_measurement_count"] == 0


def test_gpu_inventory_distinguishes_detected_available_and_allocated_devices():
    node = {
        "name": "pc002975",
        "availability": "online",
        "snapshot": {"resources": {
            "gpus": [
                {"gpu_id": 0, "name": "RTX 4090"},
                {"gpu_id": 1, "name": "RTX 4090"},
            ],
            "free_gpus": [{"gpu_id": 1, "name": "RTX 4090"}],
        }},
    }

    assert cluster_ui._gpu_counts(node) == (2, 1, 1)
    assert [row["Scheduler state"] for row in cluster_ui._gpu_inventory_rows([node])] == [
        "Allocated now", "Available now"
    ]
    node["availability"] = "offline"
    assert [row["Scheduler state"] for row in cluster_ui._gpu_inventory_rows([node])] == [
        "Last reported allocated", "Last reported available"
    ]
    node["availability"] = "online"
    node["snapshot"]["resources"].pop("free_gpus")
    assert cluster_ui._gpu_counts(node) == (2, None, None)
    assert {row["Scheduler state"] for row in cluster_ui._gpu_inventory_rows([node])} == {
        "Allocation status not reported"
    }


def test_workstation_settings_show_the_verified_four_gpu_pool(tmp_path, monkeypatch):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    cluster_config.save({"mode": "peer", "coordinator_url": "http://127.0.0.1:18766"})
    nodes = [
        {"id": "pc002975", "name": "pc002975", "availability": "online", "last_seen_at": "now", "snapshot": {"accept_jobs": True, "resources": {
            "gpus": [{"gpu_id": 0, "name": "RTX 4090"}, {"gpu_id": 1, "name": "RTX 4090"}],
            "free_gpus": [{"gpu_id": 0, "name": "RTX 4090"}], "cpu_free": 16,
        }}},
        {"id": "pc004032", "name": "pc004032", "availability": "online", "last_seen_at": "now", "snapshot": {"accept_jobs": True, "resources": {
            "gpus": [{"gpu_id": 0, "name": "RTX 5090"}], "free_gpus": [], "cpu_free": 16,
        }}},
        {"id": "pc004033", "name": "pc004033", "availability": "online", "last_seen_at": "now", "snapshot": {"accept_jobs": True, "resources": {
            "gpus": [{"gpu_id": 0, "name": "RTX 5090"}], "free_gpus": [{"gpu_id": 0, "name": "RTX 5090"}], "cpu_free": 16,
        }}},
    ]
    monkeypatch.setattr(cluster_ui, "_request", lambda *args, **kwargs: {"nodes": nodes})

    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        "from mn_molecular_design_studio.cluster.ui import settings\n"
        "settings({'id':'admin', 'role':'admin'}, 'token')"
    ).run()

    assert not app.exception
    assert any("4 detected · 2 available · 2 allocated GPUs" in item.value for item in app.caption)
    inventory = next(frame.value for frame in app.dataframe if "Model" in frame.value.columns)
    assert len(inventory) == 4
    assert list(inventory["Model"]).count("RTX 4090") == 2
    assert list(inventory["Model"]).count("RTX 5090") == 2


def test_campaign_results_panel_renders_without_mutating_evidence(tmp_path):
    from streamlit.testing.v1 import AppTest

    root = _campaign_fixture(tmp_path / "rendered", with_index=False, validation_status="failed", error="missing candidate index and lineage")
    context = campaign_context_for_job(root)
    source_before = (root / "artifacts/distributed_campaign/ledger.json").read_bytes()
    app = AppTest.from_string(
        "import streamlit as st\n"
        "from mn_molecular_design_studio.distributed_campaign_ui import render_campaign_context\n"
        "render_campaign_context(st, st.session_state['context'])"
    )
    app.session_state["context"] = context
    rendered = app.run()

    assert not rendered.exception
    assert any("Validation-only follow-up" in item.value for item in rendered.markdown)
    assert any("Global acceptance progress" in item.value for item in rendered.caption)
    assert any("not a scientific candidate rejection" in item.value for item in rendered.warning)
    assert (root / "artifacts/distributed_campaign/ledger.json").read_bytes() == source_before

    from io import BytesIO
    from zipfile import ZipFile
    from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import build_run_zip

    export_root = _campaign_fixture(tmp_path / "exported", with_index=True)
    with ZipFile(BytesIO(build_run_zip(export_root))) as archive:
        names = set(archive.namelist())
    assert "artifacts/distributed_campaign/ledger.json" in names
    assert "artifacts/normalized_candidates/candidates.jsonl" in names


def test_rfd4_tile_campaign_controls_expose_independent_settings_before_submit():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        "import streamlit as st\n"
        "from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG\n"
        "from mn_molecular_design_studio.pipeline_node_ui import render_node_form\n"
        "def make_node(kind, identifier):\n"
        "    spec = NODE_CATALOG[kind]\n"
        "    return {'id':identifier,'type':'studio','position':{'x':0,'y':0},'data':{'kind':kind,'label':spec['label'],'params':{f['name']:f['default'] for f in spec['fields']}}}\n"
        "target = make_node('target','target')\n"
        "tile = make_node('rfd4','rfd4')\n"
        "tile['data']['params']['seed'] = 17\n"
        "graph = {'nodes':[target,tile],'edges':[{'id':'edge','source':'target','target':'rfd4'}]}\n"
        "values, _, _, _ = render_node_form(tile, graph, key='tile', render_input=lambda: {}, owner_id='fixture')\n"
        "st.session_state['effective_config'] = values.get('distributed_campaign')\n"
        "st.session_state['config_ready'] = 'distributed_campaign' in values"
    ).run()
    assert not app.exception
    assert any(widget.label == "Use distributed campaign" for widget in app.checkbox)
    app.checkbox(key="tile-campaign-enabled").check().run()
    assert not app.exception
    assert app.session_state["effective_config"]["objective_type"] == "generated"
    assert app.session_state["effective_config"]["validation_stages"] == []
    app.selectbox(key="tile-campaign-objective-type").set_value("accepted").run()
    assert any(widget.label == "Target unique passing designs" for widget in app.number_input)
    assert any(widget.label == "Maximum campaign attempts" for widget in app.number_input)
    assert any(widget.label == "Scientific acceptance policy" for widget in app.selectbox)
    assert any(widget.label == "Maximum concurrent generation jobs" for widget in app.number_input)
    assert any(widget.label == "Concurrent downstream validations" for widget in app.number_input)
    assert any(widget.label == "Worker workstation placement" for widget in app.selectbox)
    assert any("Effective campaign settings before submission" in item.value for item in app.markdown)

    app.number_input(key="tile-campaign-accepted-target").set_value(3)
    app.number_input(key="tile-campaign-max-attempts").set_value(7)
    app.number_input(key="tile-campaign-generation-job-concurrency").set_value(2)
    app.number_input(key="tile-campaign-validation-concurrency").set_value(4)
    app.selectbox(key="tile-campaign-workstation-placement").set_value("eligible")
    app.run()
    config = app.session_state["effective_config"]
    assert app.session_state["config_ready"] is True
    assert config["accepted_target"] == 3
    assert config["max_attempts"] == 7
    assert config["acceptance_level"] == "validated"
    assert config["max_concurrent_generation_jobs"] == 2
    assert config["batch_size"] == 1
    assert config["max_concurrent_batches"] == 1
    assert config["max_concurrent_validations"] == 4
    assert config["workstation_selection"] == "eligible"
    assert config["seed"] == 17
    assert config["generator_id"] == "rfd4_proteina"
    assert config["validation_stages"][0]["required_metrics"] == ["iptm", "ptm", "plddt_mean"]
    assert config["unique_accepted"] is True

    app.selectbox(key="tile-campaign-objective-type").set_value("generated").run()
    assert not app.exception
    generated = app.session_state["effective_config"]
    assert generated["objective_type"] == "generated"
    assert generated["validation_stages"] == []
    assert generated["unique_accepted"] is False
    assert any(widget.label == "Target designs" for widget in app.number_input)

    app.selectbox(key="tile-campaign-objective-type").set_value("time").run()
    assert not app.exception
    app.number_input(key="tile-campaign-duration-minutes").set_value(30).run()
    timed = app.session_state["effective_config"]
    assert timed["objective_type"] == "time"
    assert timed["accepted_target"] is None
    assert timed["validation_stages"] == []
    assert timed["max_runtime_seconds"] > 0
    assert not any(widget.label == "Target unique passing designs" for widget in app.number_input)

    bindcraft = AppTest.from_string(
        "import streamlit as st\n"
        "from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG\n"
        "from mn_molecular_design_studio.pipeline_node_ui import render_node_form\n"
        "def make_node(kind, identifier):\n"
        "    spec = NODE_CATALOG[kind]\n"
        "    return {'id':identifier,'type':'studio','position':{'x':0,'y':0},'data':{'kind':kind,'label':spec['label'],'params':{f['name']:f['default'] for f in spec['fields']}}}\n"
        "target = make_node('target','target')\n"
        "tile = make_node('bindcraft2','bindcraft2')\n"
        "graph = {'nodes':[target,tile],'edges':[{'id':'edge','source':'target','target':'bindcraft2'}]}\n"
        "values, _, _, _ = render_node_form(tile, graph, key='tile', render_input=lambda: {}, owner_id='fixture')\n"
        "st.session_state['effective_config'] = values.get('distributed_campaign')"
    ).run()
    assert not bindcraft.exception
    bindcraft.checkbox(key="tile-campaign-enabled").check().run()
    assert not bindcraft.exception
    assert bindcraft.session_state["effective_config"]["generator_id"] == "bindcraft2"
    assert bindcraft.session_state["effective_config"]["acceptance_level"] == "native"


def test_selected_shared_branch_uses_graph_objective_instead_of_independent_campaign_form():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        "import streamlit as st\n"
        "from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG\n"
        "from mn_molecular_design_studio.pipeline_node_ui import render_node_form\n"
        "def make_node(kind, identifier):\n"
        "    spec = NODE_CATALOG[kind]\n"
        "    return {'id':identifier,'type':'studio','position':{'x':0,'y':0},'data':{'kind':kind,'label':spec['label'],'params':{f['name']:f['default'] for f in spec['fields']}}}\n"
        "target = make_node('target','target')\n"
        "rfd4 = make_node('rfd4','rfd4')\n"
        "bindcraft = make_node('bindcraft2','bindcraft')\n"
        "graph = {'nodes':[target,rfd4,bindcraft],'edges':[{'id':'e1','source':'target','target':'rfd4'},{'id':'e2','source':'target','target':'bindcraft'}], 'shared_campaign':{'enabled':True,'generator_nodes':['rfd4','bindcraft']}}\n"
        "render_node_form(rfd4, graph, key='tile', render_input=lambda: {}, owner_id='fixture')\n"
    ).run()

    assert not app.exception
    assert any("graph’s shared campaign objective" in item.value for item in app.info)
    assert not any(widget.label == "Use distributed campaign" for widget in app.checkbox)


def test_campaign_lifecycle_controls_offer_pause_resume_and_cancel_only_in_valid_states(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs

    run_dir = tmp_path / "campaign-run"
    run_dir.mkdir()
    (run_dir / "worker_request.json").write_text("{}", encoding="utf-8")
    calls = []
    monkeypatch.setattr(jobs, "pause_job", lambda folder, **kwargs: calls.append(("pause", Path(folder))))
    monkeypatch.setattr(jobs, "resume_job", lambda folder: calls.append(("resume", Path(folder))))

    script = (
        "import streamlit as st\n"
        "from pathlib import Path\n"
        "from mn_molecular_design_studio.distributed_campaign_controls_ui import render_campaign_lifecycle_controls\n"
        "render_campaign_lifecycle_controls(st, reference='local:design:fixture', run_dir=st.session_state['run_dir'], "
        "status=st.session_state['status'], context={'campaign_id':'fixture'}, "
        "cancel=lambda ref: st.session_state.__setitem__('cancelled', ref))"
    )
    active = AppTest.from_string(script)
    active.session_state["run_dir"] = str(run_dir)
    active.session_state["status"] = "running"
    active.run()
    assert not active.exception
    assert any(button.label == "Pause campaign controller" for button in active.button)
    assert any(button.label == "Cancel campaign" for button in active.button)
    active.button(key="campaign-pause-local:design:fixture").click().run()
    assert calls[-1] == ("pause", run_dir)

    paused = AppTest.from_string(script)
    paused.session_state["run_dir"] = str(run_dir)
    paused.session_state["status"] = "paused"
    paused.run()
    assert any(button.label == "Resume campaign" for button in paused.button)
    assert not any(button.label == "Cancel campaign" for button in paused.button)
    paused.button(key="campaign-resume-local:design:fixture").click().run()
    assert calls[-1] == ("resume", run_dir)

    terminal = AppTest.from_string(script)
    terminal.session_state["run_dir"] = str(run_dir)
    terminal.session_state["status"] = "cancelled"
    terminal.run()
    assert not any(button.label in {"Resume campaign", "Cancel campaign", "Pause campaign controller"} for button in terminal.button)
