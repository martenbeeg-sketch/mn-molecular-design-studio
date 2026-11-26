import copy
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import logical_runs as runs, run_operations


def row(identifier, **values):
    return dict(dict(run_id=identifier, task_group="design", node_id="00000000-0000-0000-0000-000000000001", owner_id="u", project_id="p",
                tool="generator", status="completed", campaign_name=identifier, has_candidate_collection=True),
                **values)


def test_one_run_contains_native_jobs_and_separate_analysis_outputs():
    parent = row("controller")
    child = row("batch", parent_run_id="controller", parent_node_id="00000000-0000-0000-0000-000000000001")
    analysis = {**row("scored"), "task_group": "analysis", "campaign_id": "controller"}
    result = runs.group_runs([child, analysis, parent])
    assert len(result) == 1 and result[0]["job_count"] == 3
    assert len(result[0]["outputs"]) == 2
    original = result[0]["outputs"][0]
    assert [r["run_id"] for r in original["sources"]] == ["controller"]
    assert result[0]["outputs"][1]["sources"][0]["run_id"] == "scored"


def test_legacy_independent_batches_are_all_sources_even_after_analysis():
    first = row("a", campaign_group_id="r3", campaign_group_name="R3")
    second = {**row("b", campaign_group_id="r3"), "node_id": "00000000-0000-0000-0000-000000000002"}
    score = {**row("scored"), "task_group": "analysis", "campaign_id": "a"}
    run = runs.group_runs([first, score, second])[0]
    assert run["job_count"] == 3
    assert {r["run_id"] for r in run["outputs"][0]["sources"]} == {"a", "b"}


def test_parent_link_does_not_cross_project_owner_or_ambiguous_workstations():
    parent = row("root")
    for child in [row("child", parent_run_id="root", project_id="other"),
                  row("child", parent_run_id="root", owner_id="other"),
                  row("child", parent_run_id="root", parent_node_id="unknown")]:
        assert len(runs.group_runs([parent, child])) == 2
    other = {**parent, "node_id": "00000000-0000-0000-0000-000000000002"}
    child = {**row("child", parent_run_id="root"), "node_id": "00000000-0000-0000-0000-000000000003"}
    assert len(runs.group_runs([parent, other, child])) == 3


def test_singletons_and_distributed_batches_have_stable_ids_and_no_input_mutation():
    rows = [row("a", campaign_group_id="g"), row("b", campaign_group_id="g")]
    original = copy.deepcopy(rows)
    before = runs.group_runs(rows)[0]
    after = runs.group_runs(list(reversed(rows)))[0]
    assert before["id"] == after["id"]
    assert len(before["outputs"][0]["sources"]) == 2
    assert rows == original


def test_graph_steps_are_named_outputs_and_the_combined_result_is_not_duplicated():
    parent = {**row("pipeline"), "tool": "visual_pipeline", "task_group": "campaign"}
    children = [{**row(f"prediction-{i}", campaign_id="pipeline", pipeline_parent_node_id=parent["node_id"],
                        graph_node_id="predictor", campaign_step="Independent prediction"), "task_group": "refolding"} for i in range(2)]
    run = runs.group_runs([*children, parent])[0]
    assert len(run["outputs"]) == 2
    assert run["outputs"][0]["sources"][0]["run_id"] == "pipeline"
    assert run["outputs"][1]["name"] == "Independent prediction"
    assert len(run["outputs"][1]["sources"]) == 2


def test_saved_source_batches_have_one_named_step_output():
    parent = row("pipeline", tool="visual_pipeline", task_group="campaign")
    children = [row(f"ensemble-{i}", tool="bioemu2", task_group="simulation", campaign_id="pipeline",
                    graph_node_id=f"sample_source_{i:06d}_a0") for i in range(2)]
    run, = runs.group_runs([parent, *children])
    results, propagation = run["outputs"]
    assert results["sources"] == [run["root"]]
    assert {job["run_id"] for job in propagation["sources"]} == {"ensemble-0", "ensemble-1"}


def test_pipeline_reusing_a_saved_input_keeps_its_own_run_and_steps():
    original = row("original", task_group="campaign", tool="visual_pipeline")
    saved_trajectory = row("trajectory", task_group="ligand-md", campaign_id="original")
    analysis = row("reused", task_group="campaign", tool="visual_pipeline",
                   dependency_refs=[{"run_id": "trajectory", "task_group": "ligand-md"}])
    worker = row("analysis-worker", task_group="md-analysis", campaign_id="reused",
                 graph_node_id="md-analysis")
    grouped = runs.group_runs([original, saved_trajectory, analysis, worker])
    assert {run["name"]: {job["run_id"] for job in run["jobs"]} for run in grouped} == {
        "original": {"original", "trajectory"},
        "reused": {"reused", "analysis-worker"},
    }


def test_catalog_groups_legacy_pose_worker_by_pipeline_parent_before_source_campaign(isolated, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, submission_tags
    from mn_molecular_design_studio.job_catalog import _job_payload

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(isolated / "runs"))
    node_id = row("source")["node_id"]
    with submission_tags(campaign_id="source-campaign", pipeline_parent_node_id=node_id,
                         pipeline_parent_run_dir="runs:///campaign/analysis-pipeline", graph_node_id="contacts"):
        job = create_job("analysis", "analysis", "native_ligand_interactions", {}, {})
    payload = _job_payload(job.run_dir)
    worker = row(job.run_dir.name, task_group="analysis", campaign_id=payload["campaign_id"],
                 pipeline_parent_node_id=payload["pipeline_parent_node_id"], **json.loads(payload["run_identity_json"]))
    parent = row("analysis-pipeline", tool="visual_pipeline", task_group="campaign")
    source = row("source-campaign", tool="visual_pipeline", task_group="campaign")

    grouped = runs.group_runs([source, worker, parent])
    assert {run["name"]: {item["run_id"] for item in run["jobs"]} for run in grouped} == {
        "analysis-pipeline": {"analysis-pipeline", job.run_dir.name},
        "source-campaign": {"source-campaign"},
    }


def test_completed_pipeline_keeps_worker_retry_history_without_becoming_partial():
    root = row("pipeline", task_group="campaign", tool="visual_pipeline")
    worker = row("replacement", campaign_id="pipeline", graph_node_id="simulation")
    first_attempt = row("first", campaign_id="pipeline", graph_node_id="simulation", status="cancelled")
    grouped, = runs.group_runs([root, worker, first_attempt])
    assert grouped["status"] == "completed"
    assert {job["run_id"] for job in grouped["jobs"]} == {"pipeline", "first", "replacement"}


def test_remote_graph_receipts_are_scoped_to_the_run_owner_and_project(isolated):
    root = dict(row("pipeline"), task_group="campaign", tool="visual_pipeline",
                node_id="organizer", pipeline_batch_jobs={"submitted-id": "refold"})
    child = row("batch", node_id="worker", cluster_job_id="submitted-id")
    run = runs.group_runs([child, root])[0]
    assert run["root"]["run_id"] == "pipeline"
    assert run["job_count"] == 2
    assert [item["run_id"] for item in run["outputs"][0]["sources"]] == ["pipeline"]
    assert [item["run_id"] for item in run["outputs"][1]["sources"]] == ["batch"]
    for other in (dict(child, owner_id="other"), dict(child, project_id="other"),
                  dict(child, cluster_job_id="unrelated")):
        assert len(runs.group_runs([other, root])) == 2
    duplicate = dict(root, run_id="different-launch")
    assert len(runs.group_runs([child, root, duplicate])) == 3


def test_cancelled_unaccepted_batch_is_execution_history_without_a_data_source():
    completed = row("finished", campaign_group_id="run-group")
    unaccepted = dict(cluster_job_id="not-started", campaign_group_id="run-group", tool="generator",
                      campaign_name="Cancelled batch", project_id="p", owner_id="u", status="cancelled")
    run = runs.group_runs([completed, unaccepted])[0]
    assert run["job_count"] == 2
    sources = run["outputs"][0]["sources"]
    cancelled = next(r for r in sources if r.get("cluster_job_id") == "not-started")
    assert cancelled["has_candidate_collection"] is False
    assert runs.source_reference(cancelled) == "cluster:not-started"
    _, _, selected, excluded = run_operations.prepare(dict(workflow="analysis.ipsae", source_job=run["id"], parameters={}), [run])
    assert len(selected) == 1 and len(excluded) == 1


def test_ledger_repair_groups_exact_worker_jobs_and_keeps_proxy_out_of_results(isolated):
    from mn_molecular_design_studio.run_reconciliation import ledger_memberships, merge_memberships
    from mn_molecular_design_studio.cluster.config import atomic_json
    parent = row("controller", distributed_campaign_id="campaign", has_candidate_collection=False)
    ledger = dict(campaign_id="campaign", attempts={"1": dict(cluster_job_id="native-job", cluster_remote=True,
        child_run_dir="/runs/design/proxy")})
    memberships = ledger_memberships(parent, ledger)
    (isolated / "studio").mkdir()
    atomic_json(isolated / "studio/run_membership.json", dict(memberships=memberships))
    worker = row("native", node_id="remote", cluster_job_id="native-job")
    proxy = row("proxy", hidden=True)
    unrelated = row("other", cluster_job_id="native-job", owner_id="other")
    result = runs.group_runs(runs.apply_operation_links([parent, worker, proxy, unrelated]))
    grouped = next(r for r in result if r["job_count"] == 3)
    assert {r["run_id"] for r in grouped["outputs"][0]["sources"]} == {"native"}
    assert len(result) == 2
    assert merge_memberships(memberships, memberships) == memberships
    with pytest.raises(ValueError, match="Conflicting"):
        merge_memberships(memberships, {"cluster:native-job": {}})
    with pytest.raises(ValueError, match="does not belong"):
        ledger_memberships(parent, dict(campaign_id="wrong"))


def test_multi_source_cohort_is_never_silently_split():
    from mn_molecular_design_studio import run_graph
    run = runs.group_runs([row("a", campaign_group_id="g"), row("b", campaign_group_id="g")])[0]
    with pytest.raises(ValueError, match="complete cohort"):
        run_operations.prepare(dict(workflow="analysis.sequence_clustering", source_job=run["id"], parameters={}), [run])
    assert not run_graph.compatible("sequence_clustering")
    assert run_graph.compatible("sequence_qc")


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    return tmp_path


def test_operation_receipts_match_exact_executions_within_scope(isolated):
    from mn_molecular_design_studio.cluster import config

    local = config.identity()
    directory = isolated / "studio" / "run_operations"
    directory.mkdir(parents=True)

    def receipt(output, job, *, owner="u", project="p"):
        path = directory / f"op-{output}.json"
        config.atomic_json(path, dict(owner_id=owner, project_id=project,
            requests=[{"run_context": {"studio_run_output_id": output}}], receipts={"0": {"job": job}}))

    receipt("local", dict(task_group="design", run_id="shared"))
    receipt("worker", dict(cluster_job_id="cluster-job"))
    receipt("incomplete", dict(run_id="bare"))
    receipt("other-scope", dict(task_group="design", run_id="shared"), project="alternate")
    receipt("ambiguous-1", dict(task_group="design", run_id="ambiguous", node_id="worker"))
    receipt("ambiguous-2", dict(task_group="design", run_id="ambiguous", node_id="worker"))
    rows = [row("shared", node_id=local), row("shared", node_id="worker"),
            row("shared", node_id=local, task_group="analysis"), row("shared", node_id=local, owner_id="another"),
            row("shared", node_id=local, project_id="alternate"), row("bare", node_id=local),
            row("worker-run", node_id="worker", cluster_job_id="cluster-job"),
            row("worker-run", node_id="worker", cluster_job_id="cluster-job", owner_id="another"),
            row("ambiguous", node_id="worker")]
    original = copy.deepcopy(rows)
    receipts = {path: path.read_bytes() for path in directory.glob("*.json")}
    linked = runs.apply_operation_links(rows)
    assert [item.get("studio_run_output_id") for item in linked] == [
        "local", None, None, None, "other-scope", None, "worker", None, None,
    ]
    assert rows == original
    assert all(path.read_bytes() == content for path, content in receipts.items())


def test_run_dispatch_freezes_all_sources_recovers_and_keeps_remote_affinity(isolated):
    run = runs.group_runs([row("a", campaign_group_id="g"), row("b", campaign_group_id="g")])[0]
    request = dict(schema_version=1, workflow="analysis.ipsae", source_job=run["id"], parameters={}, project_id="p")
    calls, validations = [], []
    fail = True
    def dispatch(child):
        calls.append(child)
        if len(calls) == 2 and fail:
            raise ValueError("temporary disconnect")
        return {"job": {"run_id": "out" + str(len(calls)), "cluster_job_id": "jid" + str(len(calls))}}
    result = run_operations.submit(request, [run], owner_id="u", validate=validations.append,
                                   dispatch=dispatch, idempotency_key="once")
    assert len(validations) == 2 and result["submitted_job_count"] == 1 and result["dispatch_errors"]
    fail = False
    recovered = run_operations.submit(request, [run], owner_id="u", validate=validations.append,
                                      dispatch=dispatch, idempotency_key="once")
    assert recovered["submitted_job_count"] == 2 and not recovered["dispatch_errors"]
    assert len(calls) == 3 and calls[1]["execution"]["idempotency_key"] == calls[2]["execution"]["idempotency_key"]
    assert calls[0]["source_job"].startswith("cluster:")
    assert calls[0]["run_context"]["studio_run_id"] == run["id"]


def test_validation_failure_submits_no_jobs(isolated):
    run = runs.group_runs([row("a", campaign_group_id="g"), row("b", campaign_group_id="g")])[0]
    calls = []
    checked = []
    def validate(child):
        checked.append(child)
        if len(checked) == 2:
            raise ValueError("incompatible input")
    with pytest.raises(ValueError, match="incompatible"):
        run_operations.submit(dict(workflow="analysis.ipsae", source_job=run["id"], parameters={}), [run],
            owner_id="u", validate=validate, dispatch=calls.append)
    assert not calls


def test_multi_batch_subset_cannot_confuse_duplicate_native_ids():
    run = runs.group_runs([row("a", campaign_group_id="g"), row("b", campaign_group_id="g")])[0]
    with pytest.raises(ValueError, match="source-qualified"):
        run_operations.prepare(dict(workflow="analysis.ipsae", source_job=run["id"],
            parameters={"candidate_ids": ["design_0"]}), [run])


def test_new_jobs_publish_run_identity_and_nested_submission_tags(isolated, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, submission_tags
    from mn_molecular_design_studio.job_catalog import _job_payload, _upsert, search_jobs
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(isolated / "runs"))
    context = dict(studio_run_id=runs.stable_id("r3"), studio_run_name="R3", studio_run_output_id="scored", studio_run_output_name="Analysis")
    with submission_tags(studio_owner_id="u", project_id="p"):
        with submission_tags(**context):
            job = create_job("analysis", "analysis", "ipsae", {}, {})
    payload = _job_payload(job.run_dir)
    assert payload["owner_id"] == "u" and json.loads(payload["run_identity_json"])["studio_run_id"] == context["studio_run_id"]
    _upsert(isolated / "runs", job.run_dir, refresh_artifacts=False)
    result = search_jobs(access_admin=True)["items"][0]
    assert result["studio_run_output_id"] == "scored"
def test_run_context_preserves_count_only_distributed_output_policy():
    from mn_molecular_design_studio.logical_runs import context_tags

    value = {
        "studio_run_id": "run:" + "a" * 24,
        "studio_run_name": "count-only campaign",
        "studio_run_output_id": "results",
        "studio_run_output_name": "Results",
        "distributed_output_policy": "accepted_candidates_only",
    }

    assert context_tags(value)["distributed_output_policy"] == "accepted_candidates_only"
