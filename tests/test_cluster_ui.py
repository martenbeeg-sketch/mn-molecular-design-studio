from __future__ import annotations

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest
from typer.testing import CliRunner

from mn_molecular_design_studio.cluster import cli, config, transport, ui


@pytest.fixture
def isolated(monkeypatch, tmp_path):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    return tmp_path


def test_settings_render_without_enabling_sharing(isolated):
    at = AppTest.from_string('from mn_molecular_design_studio.cluster.ui import settings\nsettings({"id":"admin", "role":"admin"}, "token")').run()
    assert not at.exception
    assert at.selectbox[0].value == "disabled"
    assert config.load()["mode"] == "disabled"


def test_workstation_resources_are_visible_to_an_ordinary_user(isolated, monkeypatch):
    config.save({"mode": "worker", "coordinator_url": "http://127.0.0.1:18766", "default_placement": "auto"})
    node = {"id": "remote", "name": "4032", "availability": "online", "snapshot": {"resources": {
        "free_gpus": [{"gpu_id": 0}, {"gpu_id": 1}], "cpu_free": 16, "ram_available_gb": 64}}}
    monkeypatch.setattr(ui, "_request", lambda *a, **k: {"nodes": [node]})
    at = AppTest.from_string('from mn_molecular_design_studio.cluster.ui import execution_selector\nexecution_selector({"id":"scientist", "role":"user"}, "token")').run()
    assert not at.exception
    assert any("2 free GPUs" in c.value for c in at.caption)
    assert at.selectbox[0].value == "auto"


def test_shared_results_can_page_beyond_first_hundred_jobs(isolated, monkeypatch):
    rows = [{"cluster_job_id": f"job-{i}", "node_id": "remote", "node_name": "4032", "project_id": "p", "campaign_name": f"Design {i}",
             "status": "completed", "availability": "offline", "created_at": f"{i:05d}"} for i in range(121)]
    monkeypatch.setattr(ui, "catalog", lambda *a, **kw: ({"items": rows, "total": len(rows)}, False))
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "load_organization", lambda: {"projects": [{"id": "p", "name": "Project"}]})
    at = AppTest.from_string('from mn_molecular_design_studio.cluster.ui import render_runs\nrender_runs({"id":"scientist"}, "token", lambda **kw: {"items":[], "total":0})').run()
    assert not at.exception
    assert len(at.dataframe) == 1
    table = at.dataframe[0].value
    assert len(table) == 50
    assert "job-120" in table.index
    assert set(table["View"]) == {"View"}
    assert list(table.columns).index("View") == list(table.columns).index("Job") + 1
    at.number_input[0].set_value(3).run()
    assert not at.exception
    assert len(at.dataframe) == 1
    assert len(at.dataframe[0].value) == 21


def test_pipeline_results_use_one_table_and_load_only_requested_steps(isolated, monkeypatch):
    rows = []
    for index in range(6):
        root = {
            "cluster_job_id": f"pipeline-{index}", "run_id": f"run-{index}", "node_id": "remote",
            "node_name": "4032", "project_id": "p", "campaign_name": f"Pipeline {index}",
            "tool": "visual_pipeline", "status": "completed", "created_at": f"{index:05d}",
        }
        rows.append(root)
        rows.extend({
            **root, "cluster_job_id": f"child-{index}-{step}", "run_id": f"step-{index}-{step}",
            "campaign_id": root["run_id"], "tool": "boltz2", "campaign_step": f"Step {step}",
        } for step in range(3))
    rows.append({
        "cluster_job_id": "standalone", "campaign_name": "Standalone", "node_id": "remote",
        "project_id": "p", "tool": "sequence_qc", "status": "completed", "created_at": "99999",
    })
    catalog_calls = []
    elapsed_calls = []

    def catalog(*args, **kwargs):
        catalog_calls.append(True)
        return {"items": rows, "total": len(rows)}, False

    monkeypatch.setattr(ui, "catalog", catalog)
    monkeypatch.setattr(ui, "_test_elapsed", lambda row: elapsed_calls.append(row["cluster_job_id"]) or "1s", raising=False)
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "load_organization", lambda: {"projects": [{"id": "p", "name": "Project"}]})
    at = AppTest.from_string(
        'from mn_molecular_design_studio.cluster import ui\n'
        'ui.render_runs({"id":"scientist"}, "token", lambda **kw: {"items":[], "total":0}, elapsed_label=ui._test_elapsed)'
    ).run()
    assert not at.exception
    assert len(at.dataframe) == 1
    assert not at.expander
    table = at.dataframe[0].value
    assert len(table) == 7
    assert table.loc["pipeline-2", "Steps"] == "3/3 finished"
    assert table["Steps"].isna().loc["standalone"]
    assert not any(reference.startswith("child-") for reference in elapsed_calls)

    elapsed_calls.clear()
    at.session_state["cluster-results-pipeline-details-request"] = "pipeline-2"
    at.run()
    assert not at.exception
    assert len(at.dataframe) == 2
    assert list(at.dataframe[1].value["Step"]) == ["Step 0", "Step 1", "Step 2"]
    assert {reference for reference in elapsed_calls if reference.startswith("child-")} == {
        "child-2-0", "child-2-1", "child-2-2",
    }
    assert len(catalog_calls) == 1


def test_results_table_selection_selects_and_clears_the_whole_pipeline(monkeypatch):
    parent = {"cluster_job_id": "parent"}
    child = {"cluster_job_id": "child"}
    standalone = {"cluster_job_id": "standalone"}
    state = {
        "table": {"edited_rows": {0: {"Select": True}, 1: {"Select": True}}},
        "cluster-results-pending-remove-and-purge": [standalone],
    }
    monkeypatch.setattr(ui.st, "session_state", state)
    ui._update_results_table_selection("table", [[parent, child], [standalone]])
    assert set(state["cluster-results-selected"]) == {"parent", "child", "standalone"}
    assert "cluster-results-pending-remove-and-purge" not in state
    state["table"] = {"edited_rows": {0: {"Select": False}}}
    ui._update_results_table_selection("table", [[parent, child], [standalone]])
    assert set(state["cluster-results-selected"]) == {"standalone"}


def test_pipeline_step_button_uses_the_clicked_row(monkeypatch):
    state = {"steps-click": {"row": 1, "label": "3/3 finished"}}
    monkeypatch.setattr(ui.st, "session_state", state)
    ui._request_results_pipeline_steps("steps-click", [None, "pipeline-2"])
    assert state["cluster-results-pipeline-details-request"] == "pipeline-2"


def test_campaign_batches_group_across_workstations_and_keep_campaign_scopes_separate():
    rows = [
        {"cluster_job_id": "batch-a", "campaign_name": "R3 IDR · RFD4 10-design batch 01 · pc004032 GPU0 · 2026-10-30",
         "node_id": "pc004032", "project_id": "michele", "owner_id": "user", "tool": "rfd4_proteina"},
        {"cluster_job_id": "batch-b", "campaign_name": "R3 IDR · RFD4 10-design batch 01 · pc004033 GPU0 · 2026-10-30",
         "node_id": "pc004033", "project_id": "michele", "owner_id": "user", "tool": "rfd4_proteina"},
        {"cluster_job_id": "batch-c", "campaign_name": "R3 IDR · RFD4 10-design batch 02 · pc002975 GPU1 · 2026-10-30",
         "node_id": "pc002975", "project_id": "michele", "owner_id": "user", "tool": "rfd4_proteina"},
        {"cluster_job_id": "other-date", "campaign_name": "R3 IDR · RFD4 10-design batch 01 · pc004032 GPU0 · 2026-10-31",
         "node_id": "pc004032", "project_id": "michele", "owner_id": "user", "tool": "rfd4_proteina"},
        {"cluster_job_id": "other-project", "campaign_name": "R3 IDR · RFD4 10-design batch 03 · pc004032 GPU0 · 2026-10-30",
         "node_id": "pc004032", "project_id": "other", "owner_id": "user", "tool": "rfd4_proteina"},
    ]
    groups, row_to_group, labels = ui._campaign_batch_groups(rows)
    assert len(groups) == 1
    group_id, members = next(iter(groups.items()))
    assert [row["cluster_job_id"] for row in members] == ["batch-a", "batch-b", "batch-c"]
    assert {row_to_group[row["cluster_job_id"]] for row in members} == {group_id}
    assert labels[group_id] == "R3 IDR · RFD4 10-design · 2026-10-30"
    explicit, _, explicit_labels = ui._campaign_batch_groups([
        {"cluster_job_id": "native-a", "request_summary": {"parameters": {
            "campaign_group_id": "stable-campaign-id", "campaign_group_name": "Native design campaign"}},
         "project_id": "michele", "owner_id": "user", "tool": "boltzgen"},
        {"cluster_job_id": "native-b", "request_summary": {"payload": {"parameters": {
            "campaign_group_id": "stable-campaign-id", "campaign_group_name": "Native design campaign"}}},
         "project_id": "michele", "owner_id": "user", "tool": "boltzgen"},
    ])
    assert len(explicit) == 1
    assert next(iter(explicit_labels.values())) == "Native design campaign"


def test_campaign_results_summarize_shared_and_varying_parameters(monkeypatch):
    members = [
        {"cluster_job_id": "first", "campaign_name": "R3 IDR · batch 01", "campaign_parameters": {
            "target_type": "IDR", "num_designs": 10, "seed": 1}, "status": "completed"},
        {"cluster_job_id": "second", "campaign_name": "R3 IDR · batch 02", "request_summary": {"payload": {"parameters": {
            "target_type": "IDR", "num_designs": 10, "seed": 2}}}, "status": "failed"},
    ]
    table = {row["Parameter"]: row["Value"] for row in ui._campaign_parameter_rows(members)}
    assert table["Target Type"] == "IDR (2 jobs)"
    assert table["Num Designs"] == "10 (2 jobs)"
    assert table["Seed"] == "1 (1 jobs) · 2 (1 jobs)"
    assert ui._campaign_group_status(members) == ("Partial", "1/2 complete · 1 failed")
    state = {"campaign-results-click": {"row": 0}}
    monkeypatch.setattr(ui.st, "session_state", state)
    ui._request_results_pipeline_results("campaign-results-click", [{"kind": "campaign", "members": members}])
    assert state["cluster-results-campaign-details-request"] == members


def test_campaign_candidate_results_show_measured_values_and_keep_native_metrics_distinct():
    records = [
        {"job": {"campaign_name": "Campaign · batch 01", "tool": "protenix_v1", "node_name": "pc004032"},
         "candidate": {"candidate_id": "design-a", "source_tool": "protenix_v1", "stage": "complex_refolding",
                       "binder_sequence": "ACDE", "metrics": {"iptm": 0.72, "plddt_mean": 84.9, "pae": 6.1}}},
        {"job": {"campaign_name": "Campaign · batch 02", "tool": "mol1", "node_name": "pc004033"},
         "candidate": {"candidate_id": "design-b", "source_tool": "mol1", "stage": "complex_refolding",
                       "binder_sequence": "FGHI", "metrics": {"iptm": 0.81, "plddt_mean": 88.314,
                                                                  "complex_plddt": 88.314, "pae": 5.2},
                       "raw_metadata": {"confidence_outputs": {"plddt_mean": 84.976}}}},
        {"job": {"campaign_name": "Campaign · batch 03", "tool": "mol1", "node_name": "pc004032"},
         "candidate": {"candidate_id": "design-c", "source_tool": "mol1", "stage": "complex_refolding",
                       "binder_sequence": "KLMN", "metrics": {"iptm": 0.5}}},
    ]
    rows, definitions = ui._campaign_candidate_result_rows(records)
    by_id = {row["Candidate"]: row for row in rows}
    assert by_id["design-a"]["Validation complex iptm (0–1)"] == 0.72
    assert by_id["design-a"]["Mean structure Cα pLDDT (0–100)"] == 84.9
    assert '"pae": 6.1' in by_id["design-a"]["Native metrics"]
    assert by_id["design-b"]["Mean structure Cα pLDDT (0–100)"] == 88.314
    assert by_id["design-b"]["Mol-1 model summary pLDDT (0–100)"] == 84.976
    assert by_id["design-b"]["Validation complex iptm (0–1)"] == 0.81
    assert by_id["design-c"]["Mean structure Cα pLDDT (0–100)"] == "—"
    assert any("Validation complex iptm" in item["Metric"] and "protenix" in item["Methods"].casefold()
               for item in definitions)
    assert all("pass" not in item["Metric"].casefold() for item in definitions)


def test_catalog_campaign_parameter_summary_keeps_key_settings_without_sequences_or_paths():
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import campaign_result_parameters

    summary = campaign_result_parameters({
        "target_name": "R3",
        "target_type": "IDR",
        "num_designs": 10,
        "seed": 42,
        "hotspots": ["A25", "A30-35"],
        "target_sequence": "ACDEFGHIKLMNPQRSTVWY",
        "target_pdb": "/reference_files/r3.pdb",
        "api_token": "private-token",
    })
    assert summary == {
        "target_name": "R3", "target_type": "IDR", "num_designs": 10,
        "seed": 42, "hotspots": ["A25", "A30-35"],
    }


def test_legacy_campaign_batch_results_collapse_into_one_row(isolated, monkeypatch):
    rows = [
        {"cluster_job_id": f"batch-{index}", "campaign_name": f"PD-L1 · RFD4 10-design batch {index:02d} · pc00403{index} GPU0 · 2026-10-30",
         "node_id": f"pc00403{index}", "node_name": f"pc00403{index}", "project_id": "p", "owner_id": "user",
         "tool": "rfd4_proteina", "status": "completed", "availability": "online",
         "created_at": f"2026-10-30T00:0{index}:00Z", "campaign_parameters": {"target_type": "IDR", "num_designs": 10}}
        for index in range(1, 4)
    ]
    rows.append({"cluster_job_id": "unrelated", "campaign_name": "Standalone design", "node_id": "pc004032",
                 "project_id": "p", "owner_id": "user", "tool": "boltz2", "status": "completed",
                 "availability": "online", "created_at": "2026-10-30T00:00:00Z"})
    monkeypatch.setattr(ui, "catalog", lambda *a, **kw: ({"items": rows, "total": len(rows)}, False))
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "load_organization", lambda: {"projects": [{"id": "p", "name": "Project"}]})
    at = AppTest.from_string(
        'from mn_molecular_design_studio.cluster.ui import render_runs\n'
        'render_runs({"id":"scientist"}, "token", lambda **kw: {"items":[], "total":0})'
    ).run()
    assert not at.exception
    table = at.dataframe[0].value
    assert len(table) == 2
    campaign = table[table["Job"] == "PD-L1 · RFD4 10-design · 2026-10-30"].iloc[0]
    assert campaign["Results"] == "Campaign"
    assert campaign["Steps"] == "3 batch jobs"
    assert campaign["Workstation"] == "3 computers"
    at.session_state["cluster-results-campaign-details-request"] = rows[:3]
    at.run()
    assert not at.exception
    assert any("Parameter" in frame.value.columns for frame in at.dataframe if hasattr(frame.value, "columns"))


def test_campaign_results_load_candidate_metrics_from_the_owning_workstation(isolated, monkeypatch):
    row = {
        "cluster_job_id": "batch-a", "campaign_name": "PD-L1 · batch 01", "node_id": "remote",
        "node_name": "pc004032", "project_id": "p", "owner_id": "user", "tool": "protenix_v1",
        "status": "completed", "availability": "online", "local_reference": "design/run-a",
        "campaign_parameters": {"num_designs": 2},
    }
    candidate = {
        "candidate_id": "design-42", "source_tool": "protenix_v1", "stage": "complex_refolding",
        "binder_sequence": "ACDEFG", "metrics": {"iptm": 0.74, "plddt_mean": 88.2, "pae": 7.4},
        "complex_pdb": "artifacts/design-42.pdb",
    }
    calls = []

    def request(path, token, payload=None, method="GET", **kwargs):
        calls.append((path, payload, method))
        if path == "/jobs/batch-a/actions" and method == "POST":
            return {"command_id": "inspect-a"}
        if path == "/commands/inspect-a":
            return {"status": "completed", "response": {"results": {"candidates": [candidate]}}}
        raise AssertionError(f"Unexpected cluster request: {path} {payload} {method}")

    monkeypatch.setattr(ui, "_request", request)
    at = AppTest.from_string(
        'from mn_molecular_design_studio.cluster.ui import _render_campaign_candidate_measurements\n'
        'member = {"cluster_job_id":"batch-a", "campaign_name":"PD-L1 · batch 01", '
        '"node_id":"remote", "node_name":"pc004032", "tool":"protenix_v1", '
        '"status":"completed", "availability":"online", "local_reference":"design/run-a"}\n'
        '_render_campaign_candidate_measurements([member], "token")'
    ).run()
    assert not at.exception
    next(button for button in at.button if button.label == "Load candidate measurements").click().run()
    assert not at.exception
    assert ("/jobs/batch-a/actions", {"action": "inspect"}, "POST") in calls
    assert any(path == "/commands/inspect-a" for path, _, _ in calls)
    result_tables = [frame.value for frame in at.dataframe if hasattr(frame.value, "columns")]
    candidate_table = next(frame for frame in result_tables if "Candidate" in frame.columns)
    assert candidate_table.iloc[0]["Candidate"] == "design-42"
    assert candidate_table.iloc[0]["Validation complex iptm (0–1)"] == 0.74
    assert candidate_table.iloc[0]["Mean structure Cα pLDDT (0–100)"] == 88.2
    assert '"pae": 7.4' in candidate_table.iloc[0]["Native metrics"]


def test_selected_pipeline_includes_new_and_hidden_steps_and_blocks_active_purge(isolated, monkeypatch):
    parent = {
        "cluster_job_id": "parent", "run_id": "pipeline-run", "tool": "visual_pipeline",
        "node_id": "remote", "project_id": "p", "status": "completed",
    }
    child = {
        **parent, "cluster_job_id": "child", "run_id": "child-run", "tool": "boltz2",
        "campaign_id": parent["run_id"], "status": "running", "hidden": True,
    }
    monkeypatch.setattr(ui, "catalog", lambda *a, **kw: ({"items": [parent, child], "total": 2}, False))
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace, "load_organization", lambda: {"projects": [{"id": "p", "name": "Project"}]})
    at = AppTest.from_string(
        'from mn_molecular_design_studio.cluster.ui import render_runs\n'
        'render_runs({"id":"scientist"}, "token", lambda **kw: {"items":[], "total":0})'
    ).run()
    at.session_state["cluster-results-selected"] = {"parent": parent}
    at.session_state["cluster-results-pending-remove-and-purge"] = [parent]
    at.run()
    assert not at.exception
    assert set(at.session_state["cluster-results-selected"]) == {"parent", "child"}
    assert next(button for button in at.button if button.label == "Remove and purge selected").disabled
    assert not any(button.label == "Confirm removal and purge" for button in at.button)


def test_remote_candidate_results_include_sequence_and_scores():
    rows = ui._candidate_result_rows([{
        "candidate_id": "rf2_model_1",
        "binder_sequence": "QVQLVESGGGLVQPGGSLRLSCAAS",
        "binder_length": 127,
        "metrics": {"pred_lddt": 0.9, "interaction_pae": 9.57, "hotspot_atom_contact_fraction": 1.0},
    }])
    assert rows == [{
        "Rank": 1,
        "Design": "rf2_model_1",
        "Sequence": "QVQLVESGGGLVQPGGSLRLSCAAS",
        "Length": 127,
        "pLDDT": 0.9,
        "Interaction PAE": 9.57,
        "Hotspot atom contacts": 1.0,
    }]


def test_candidate_structure_uses_confidence_colors_only_with_per_atom_scores():
    candidate = {"metrics": {"pred_lddt": 0.9}}
    atom_lines = "\n".join(("ATOM  " + " " * 54 + f"{score:6.2f}") for score in (55.0, 90.0))
    assert ui._candidate_structure_color(candidate) == "chain-id"
    assert ui._candidate_structure_color(candidate, atom_lines) == "plddt"
    assert ui._candidate_structure_color({"metrics": {"binder_length": 127}}, atom_lines) == "chain-id"


def test_remote_candidate_structure_path_is_run_relative_and_confined():
    data = {"job": {"run_dir": "/data/results/runs/job-1"}}
    assert ui._candidate_structure_path({
        "complex_pdb": "/data/results/runs/job-1/artifacts/design/model.pdb",
    }, data) == ("artifacts/design/model.pdb", "complex_pdb")
    assert ui._candidate_structure_path({
        "complex_pdb": "/data/results/runs/other/model.pdb",
        "binder_pdb": "../../outside.pdb",
    }, data) is None


def test_remote_job_displays_sequence_and_loads_a_representative_structure(isolated, monkeypatch):
    import importlib

    candidate = {
        "candidate_id": "rf2_model_1",
        "binder_sequence": "QVQLVESGGGLVQPGGSLRLSCAAS",
        "complex_pdb": "/data/results/runs/job-1/artifacts/design/model.pdb",
        "binder_length": 127,
        "metrics": {"pred_lddt": 0.9, "interaction_pae": 9.57},
    }
    row = {
        "cluster_job_id": "job-1", "campaign_name": "Remote antibody", "tool": "rfantibody",
        "status": "completed", "availability": "online", "node_name": "pc004032",
        "local_reference": "design/run-1", "project_id": "project-1", "artifacts": [],
    }
    structure_requests = []
    viewer_calls = []

    def request(path, token, payload=None, method="GET", **kwargs):
        if path == "/jobs/job-1":
            return row
        if path == "/jobs/job-1/actions" and payload["action"] == "inspect":
            return {"command_id": "inspect-command"}
        if path == "/commands/inspect-command":
            return {"status": "completed", "response": {
                "job": {"run_dir": "/data/results/runs/job-1"},
                "results": {"candidates": [candidate]}, "input": {}, "result": {},
            }}
        if path == "/jobs/job-1/actions" and payload["action"] == "artifact":
            structure_requests.append(payload["parameters"]["path"])
            return {"command_id": "structure-command"}
        if path == "/commands/structure-command":
            return {"status": "completed", "response": {"name": "model.pdb", "size_bytes": 12}}
        if path == "/commands/structure-command/file":
            return b"ATOM sample\n"
        raise AssertionError(f"Unexpected cluster request: {path} {payload}")

    monkeypatch.setattr(ui, "_request", request)
    viewer = importlib.import_module("mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer")
    monkeypatch.setattr(viewer, "molstar_custom_component", lambda *args, **kwargs: viewer_calls.append((args, kwargs)))
    at = AppTest.from_string(
        'from mn_molecular_design_studio.cluster.ui import render_job\n'
        'render_job("cluster:job-1", {"id":"scientist"}, "token", lambda *a, **k: {})'
    ).run()
    assert not at.exception
    next(button for button in at.button if button.label == "Load full results from workstation").click().run()
    assert not at.exception
    assert structure_requests == ["artifacts/design/model.pdb"]
    assert viewer_calls
    visualization = viewer_calls[0][0][0][0]
    assert visualization.color == "chain-id"
    assert visualization.representation_type == "cartoon"
    assert any("QVQLVESGGGLVQPGGSLRLSCAAS" in str(table.value) for table in at.dataframe)
    assert any("How to interpret these scores" in block.value for block in at.markdown)


def test_services_use_the_app_conda_and_private_configuration_home(isolated, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: isolated)
    runner = CliRunner()
    result = runner.invoke(cli.app, ["init", "--name", "workstation"])
    assert result.exit_code == 0, result.output
    result = runner.invoke(cli.app, ["install-services"])
    assert result.exit_code == 0, result.output
    unit = (isolated / ".config/systemd/user/studio-agent.service").read_text()
    import sys
    assert sys.executable in unit and str(config.home()) in unit
    secret = config.load()["node_token"]
    assert secret not in unit and secret not in runner.invoke(cli.app, ["status"]).output


def test_local_placement_refuses_unavailable_results_disk(isolated, monkeypatch):
    from mn_molecular_design_studio.cluster import service
    from mn_molecular_design_studio.cluster.agent import Agent
    config.save({"mode": "coordinator", "coordinator_url": "http://127.0.0.1:8766"})
    monkeypatch.setattr(Agent, "storage_online", lambda self: False)
    with pytest.raises(transport.ClusterUnavailable, match="Results storage"):
        service.submit_from_api({"execution": {"placement": "local"}}, {"id": "user"}, request_dir=isolated, token="", authorize_path=lambda p: None)
