import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import job_dependencies as dependencies
from mn_molecular_design_studio.cluster import config
from mn_molecular_design_studio.cluster.store import Store
from mn_molecular_design_studio.engine.workflow_engine.core import jobs, local_worker, central_scheduler


def row(run, **extra):
    return {"node_id": "node", "task_group": "design", "run_id": run, "owner_id": "instructor",
            "project_id": "research", "status": "queued", "cluster_job_id": run, **extra}


def test_transitive_cascade_traverses_completed_intermediates_without_siblings():
    root = row("root")
    middle = row("middle", upstream_task_group="design", upstream_run_id="root", status="completed")
    leaf = row("leaf", upstream_task_group="design", upstream_run_id="middle")
    sibling = row("unrelated", campaign_name="root", studio_run_id="same-display-only")
    assert {value["run_id"] for value in dependencies.descendants([root, middle, leaf, sibling], [root])} == {"root", "middle", "leaf"}
    assert {value["run_id"] for value in dependencies.descendants([root, middle, leaf], [middle])} == {"middle", "leaf"}


def test_scope_workstation_and_task_qualified_edges_and_cycle_termination():
    root = row("root")
    child = row("child", upstream_run_id="root", upstream_task_group="design")
    others = [dict(child, run_id="other-owner", owner_id="student"), dict(child, run_id="other-project", project_id="other"),
              dict(child, run_id="wrong-node", upstream_node_id="other"), dict(child, run_id="wrong-group", upstream_task_group="prediction")]
    assert len(dependencies.descendants([root, child, *others], [root])) == 2
    cycle = dict(root, upstream_run_id="child", upstream_task_group="design")
    assert len(dependencies.descendants([cycle, child], [cycle])) == 2


def test_remote_and_multiple_sources_are_explicit_not_name_inferred(tmp_path):
    refs = dependencies.projected_references({}, {"source_job": "cluster:root", "parameters": {
        "source_run_dir": "cluster://second/design/source", "job_name": "unrelated", "candidates_jsonl": "runs:///design/local/artifacts/candidates.jsonl"}}, root=tmp_path)
    child = row("child", dependency_refs=refs)
    sources = [row("root"), row("source", node_id="second"), row("local")]
    assert dependencies.parents([*sources, child])[dependencies.identity(child)] == {dependencies.identity(source) for source in sources}
    assert not dependencies.projected_references({}, {"source_run_dir": "runs:///../outside/run"}, root=tmp_path)
    assert not dependencies.projected_references({}, {"source_run_dir": "/unmanaged/data"}, root=tmp_path)


@pytest.fixture
def local(tmp_path, monkeypatch):
    root = tmp_path / "runs"
    root.mkdir()
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(root))
    monkeypatch.setattr(jobs, "runs_root", lambda: root)
    monkeypatch.setattr(config, "identity", lambda: "node")
    monkeypatch.setattr(central_scheduler, "cancel_scheduled_run", lambda path: None)
    monkeypatch.setattr(jobs, "_signal_worker_process", lambda *args: None)
    monkeypatch.setattr(jobs, "_stop_run_docker_containers", lambda path: [])
    monkeypatch.setattr(jobs, "_release_resource_lock", lambda path: None)

    def create(run, *, group="design", inputs=None, **metadata):
        path = root / group / run
        jobs.write_json(path / "metadata.json", {"status": "queued", "studio_owner_id": "instructor", "project_id": "research", **metadata})
        jobs.write_json(path / "input.json", {"inputs": inputs or {}})
        jobs.write_json(path / "worker_request.json", {"kind": "unknown-test-kind"})
        return path
    return root, create


def test_local_cancellation_preserves_completed_bytes_and_cancels_transitive_children(local):
    root, create = local
    parent = create("parent")
    done = create("done", group="analysis", status="completed", inputs={"source_run_dir": str(parent)})
    jobs.write_json(done / "result.json", {"success": True, "scores": {"ipsae": 0.8}})
    before = (done / "result.json").read_bytes(), (done / "metadata.json").read_bytes()
    leaf = create("leaf", inputs={"source_run_dir": str(done)})
    other = create("other", studio_owner_id="student", inputs={"source_run_dir": str(parent)})
    summary = jobs.stop_job(parent)
    assert jobs.read_json(parent / "metadata.json")["status"] == jobs.read_json(leaf / "metadata.json")["status"] == "cancelled"
    assert before == ((done / "result.json").read_bytes(), (done / "metadata.json").read_bytes())
    assert jobs.read_json(other / "metadata.json")["status"] == "queued"
    assert str(done) in summary["preserved_runs"]
    assert not jobs.stop_job(parent)["cancelled_runs"]


def test_late_child_and_restart_do_not_execute_after_ancestor_cancellation(local):
    root, create = local
    parent = create("parent")
    jobs.stop_job(parent)
    child = create("late", inputs={"source_run_dir": str(parent)})
    assert local_worker.run_worker_job(child) == 0
    assert jobs.read_json(child / "metadata.json")["status"] == "cancelled"
    assert local_worker.spawn_worker_for_run(child) is None
    assert "Upstream job cancelled" in jobs.read_json(child / "metadata.json")["cancel_reason"]


def test_local_queue_ignores_late_child_but_retains_independent_branch(local):
    from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import queued_worker_runs
    root, create = local
    parent = create("parent")
    jobs.stop_job(parent)
    create("late", upstream_task_group="design", upstream_run_id="parent")
    other = create("other")
    assert queued_worker_runs(root) == [other]


def test_explicit_resume_clears_ancestor_cancel_marker(local):
    root, create = local
    parent = create("parent")
    jobs.stop_job(parent)
    jobs.prepare_job_for_resume(parent)
    child = create("late", inputs={"source_run_dir": str(parent)})
    assert dependencies.cancelled_dependency(child) is None


def test_completed_native_worker_replay_remains_available_without_cancellation(local):
    root, create = local
    completed = create('completed', status='completed')
    assert dependencies.guard_execution(completed)


def test_stale_worker_status_cannot_clear_durable_cancellation(local):
    root, create = local
    parent = create('parent', status='running', cancel_requested_at='recorded')
    assert not dependencies.guard_execution(parent)
    assert jobs.read_json(parent / 'metadata.json')['status'] == 'cancelled'


def test_older_metadata_rewrite_cannot_remove_cancellation_receipt(local):
    root, create = local
    parent = create('parent')
    original = jobs.read_json(parent / 'metadata.json')
    jobs.stop_job(parent)
    jobs.write_json(parent / 'metadata.json', original)
    child = create('late', inputs={'source_run_dir': str(parent)})
    assert dependencies.cancelled_dependency(child)
    assert not dependencies.guard_execution(parent)
    jobs.prepare_job_for_resume(parent)
    assert not (parent / 'cancellation_request.json').exists()
    assert dependencies.cancelled_dependency(child) is None


def test_local_scheduler_enforces_cancellation_racing_worker_pid_publication(local, monkeypatch):
    from types import SimpleNamespace
    from mn_molecular_design_studio.engine.workflow_engine.core import scheduler, central_scheduler_worker
    root, create = local
    run = create('racing-worker', status='running')
    jobs.write_json(run / 'cancellation_request.json', {'requested_at': 'recorded'})
    state = {'return_code': None, 'terminated': False}
    process = SimpleNamespace(poll=lambda: state['return_code'])
    def terminate(current):
        assert current is process
        state.update(return_code=-15, terminated=True)
    monkeypatch.setattr(central_scheduler_worker, '_terminate', terminate)
    service = object.__new__(scheduler.WorkerScheduler)
    service.active = {str(run): scheduler.ActiveWorker(process=process, run_dir=run, request=None)}
    service._refresh_active()
    assert state['terminated'] and not service.active
    assert jobs.read_json(run / 'metadata.json')['status'] == 'cancelled'


def test_central_bridge_enforces_local_cancel_without_waiting_for_api_heartbeat(local, monkeypatch):
    import sys
    from types import SimpleNamespace
    from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler_worker
    root, create = local
    run = create('central-race', scheduler_job_id='scheduled')
    state = {'return_code': None, 'finished': []}
    process = SimpleNamespace(pid=99999999, poll=lambda: state['return_code'])
    def launch(*args, **kwargs):
        metadata = jobs.read_json(run / 'metadata.json')
        jobs.write_json(run / 'cancellation_request.json', {'requested_at': 'recorded'})
        jobs.write_json(run / 'metadata.json', metadata)
        return process
    def terminate(current):
        assert current is process
        state['return_code'] = -15
    monkeypatch.setattr(central_scheduler_worker, 'runs_root', lambda: root)
    monkeypatch.setattr(central_scheduler_worker, '_apply_allocation', lambda *args: {})
    monkeypatch.setattr(central_scheduler_worker, '_submitted_execution_context', lambda *args: (sys.executable, {}))
    monkeypatch.setattr(central_scheduler_worker.subprocess, 'Popen', launch)
    monkeypatch.setattr(central_scheduler_worker.time, 'sleep', lambda *args: None)
    monkeypatch.setattr(central_scheduler_worker, '_terminate', terminate)
    client = SimpleNamespace(finish_job=lambda *args, **kwargs: state['finished'].append(kwargs))
    central_scheduler_worker._run_claimed_job(client, {'job_id': 'scheduled', 'lease_token': 'lease',
        'payload': {'run_dir': str(run)}}, 'worker')
    assert state['return_code'] == -15
    assert state['finished'] and state['finished'][0]['success'] is False


def test_central_bridge_does_not_cancel_successful_work_after_checking_dependencies(local, monkeypatch):
    import sys
    from types import SimpleNamespace
    from mn_molecular_design_studio.engine.workflow_engine.core import central_scheduler_worker

    root, create = local
    run = create('central-success', scheduler_job_id='scheduled')
    state = {'return_code': None, 'finished': []}
    process = SimpleNamespace(pid=99999999, poll=lambda: state['return_code'])
    def finish_work(*args):
        state['return_code'] = 0
        jobs.write_json(run / 'metadata.json', {**jobs.read_json(run / 'metadata.json'), 'status': 'completed'})
        jobs.write_json(run / 'result.json', {'success': True})
    monkeypatch.setattr(central_scheduler_worker, 'runs_root', lambda: root)
    monkeypatch.setattr(central_scheduler_worker, '_apply_allocation', lambda *args: {})
    monkeypatch.setattr(central_scheduler_worker, '_submitted_execution_context', lambda *args: (sys.executable, {}))
    monkeypatch.setattr(central_scheduler_worker.subprocess, 'Popen', lambda *args, **kwargs: process)
    monkeypatch.setattr(central_scheduler_worker.time, 'sleep', finish_work)
    client = SimpleNamespace(heartbeat_worker=lambda **kwargs: None, heartbeat_job=lambda *args, **kwargs: {},
        finish_job=lambda *args, **kwargs: state['finished'].append(kwargs))
    central_scheduler_worker._run_claimed_job(client, {'job_id': 'scheduled', 'lease_token': 'lease',
        'payload': {'run_dir': str(run)}}, 'worker')
    assert state['finished'] and state['finished'][0]['success'] is True
    assert jobs.read_json(run / 'metadata.json')['scheduler_status'] == 'completed'
    assert not (run / 'cancellation_request.json').exists()


def test_grouping_keeps_standalone_refolding_output_available():
    from mn_molecular_design_studio.logical_runs import group_runs
    source = row('source', has_candidate_collection=True)
    refolding = row('refolding', task_group='refolding', has_candidate_collection=True,
                    dependency_refs=[{'node_id': 'node', 'task_group': 'design', 'run_id': 'source'}])
    runs = group_runs([source, refolding])
    assert len(runs) == 1
    assert len(runs[0]['outputs']) == 2
    assert {output['sources'][0]['run_id'] for output in runs[0]['outputs']} == {'source', 'refolding'}
    assert not any('_dependency_parented' in value for value in runs[0]['jobs'])


def test_catalog_projects_dependencies_and_projects_group_separate_analyses(local, monkeypatch):
    from mn_molecular_design_studio import job_catalog
    from mn_molecular_design_studio.projects import project_inventory
    root, create = local
    parent = create("parent", campaign_name="R3 campaign")
    create("analysis", group="analysis", inputs={"source_run_dir": str(parent)})
    monkeypatch.setattr(job_catalog, "_results_root", lambda: root)
    rows = [dict(job_catalog._job_payload(path.parent), node_id="node") for path in root.glob("*/*/metadata.json")]
    for value in rows:
        value.update(json.loads(value["run_identity_json"]))
    item = project_inventory([{"id": "research", "name": "Research"}], rows, [])['research']
    assert item['runs'] == 1 and item['logical_runs'][0]['job_count'] == 2


@pytest.fixture
def cluster(tmp_path, monkeypatch):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    store = Store()
    nodes = [store.enroll(name)["node_id"] for name in ("first", "second")]
    user = {"id": "instructor", "role": "user", "project_ids": ["research"]}
    return store, nodes, user


def publish(store, node, run, **extra):
    value = row(run, node_id=node, **extra)
    store.heartbeat(node, {}, [value])
    return value


def test_cluster_cascade_queues_remote_stops_once_and_preserves_results(cluster):
    store, nodes, user = cluster
    root = publish(store, nodes[0], "root", status="running")
    publish(store, nodes[1], "middle", status="completed", upstream_task_group="design", upstream_run_id="root", upstream_node_id=nodes[0])
    publish(store, nodes[1], "leaf", upstream_task_group="design", upstream_run_id="middle")
    publish(store, nodes[1], "unrelated")
    response = store.action("root", "cancel", {}, user)
    assert response['downstream_count'] == 1 and set(response['cancelled_jobs']) == {'root', 'leaf'}
    assert store.job('middle', user)['status'] == 'completed'
    assert not store.job('unrelated', user).get('cancel_requested')
    store.action("root", "cancel", {}, user)
    with store.connect() as connection:
        assert connection.execute("SELECT count(*) FROM commands WHERE kind='cancel'").fetchone()[0] == 2
    command = store.claim(nodes[1])
    assert command['kind'] == 'cancel' and command['payload']['local_reference'] == 'design/leaf'


def test_legacy_remote_parent_matches_grouping_but_ambiguous_identity_does_not():
    parent = row('root', node_id='organizer', distributed_campaign_id='campaign')
    child = row('child', node_id='worker', campaign_id='campaign')
    assert len(dependencies.descendants([parent, child], [parent])) == 2
    duplicate = dict(parent, node_id='other')
    assert dependencies.descendants([parent, duplicate, child], [parent]) == [parent]


def test_cluster_dependency_projection_excludes_artifacts_and_keeps_full_payload(cluster):
    store, nodes, user = cluster
    artifacts = [{'path': 'x' * 10000}]
    publish(store, nodes[0], 'root', status='running', artifacts=artifacts, metrics={'score': 0.8})
    with store.connect() as connection:
        projected = store._dependency_rows(connection, user['id'], 'research')
    assert len(json.dumps(projected)) < 2000
    store.action('root', 'cancel', {}, user)
    current = store.job('root', user)
    assert current['artifacts'] == artifacts and current['metrics'] == {'score': 0.8}


def test_cluster_late_heartbeat_inherits_cancellation_after_restart(cluster):
    store, nodes, user = cluster
    publish(store, nodes[0], "root", status="running")
    store.action("root", "cancel", {}, user)
    restarted = Store(store.root)
    publish(restarted, nodes[1], "late", upstream_run_id="root", upstream_task_group="design", upstream_node_id=nodes[0])
    assert restarted.job('late', user)['cancel_requested']
    assert restarted.claim(nodes[1])['job_id'] == 'late'


def test_pending_cluster_descendant_cancelled_before_placement(cluster):
    store, nodes, user = cluster
    publish(store, nodes[0], "root", status="running")
    request = {"kind": "workflow", "project_id": "research", "payload": {"workflow": "analysis.ipsae", "parameters": {},
               "source_job": "cluster:root"}}
    pending = store.submit(request, {}, user, origin_node=nodes[0])["job"]
    result = store.action("root", "cancel", {}, user)
    assert pending['cluster_job_id'] in result['cancelled_jobs']
    assert store.job(pending['cluster_job_id'], user)['status'] == 'cancelled'
    new = store.submit(request, {}, user, origin_node=nodes[0])["job"]
    assert new['status'] == 'cancelled'


def test_cluster_cascade_never_crosses_owner_or_project(cluster):
    store, nodes, user = cluster
    publish(store, nodes[0], "root", status="running")
    publish(store, nodes[1], "other", upstream_run_id="root", upstream_task_group="design", upstream_node_id=nodes[0], owner_id='student')
    publish(store, nodes[1], "project", upstream_run_id="root", upstream_task_group="design", upstream_node_id=nodes[0], project_id='other')
    response = store.action('root', 'cancel', {}, user)
    assert response['downstream_count'] == 0
    assert not store.job('other', user).get('cancel_requested')


def test_terminal_cancel_race_does_not_overwrite_successful_cluster_result(cluster):
    store, nodes, user = cluster
    publish(store, nodes[0], "root", status="completed", metrics={"score": 0.8})
    before = store.job('root', user)
    assert store.action('root', 'cancel', {}, user)['preserved']
    assert store.job('root', user) == before


def test_current_access_and_ownership_are_rechecked_before_projects_cancel(local, monkeypatch):
    from mn_molecular_design_studio import auth, job_catalog
    root, create = local
    parent = create('parent')
    value = row('parent')
    monkeypatch.setattr(config, 'load', lambda: {'mode': 'disabled'})
    monkeypatch.setattr(job_catalog, 'get_job_access_record', lambda *args: {'owner_id': 'student', 'project_id': 'research'})
    with pytest.raises(ValueError, match='ownership changed'):
        dependencies.request_cancellation(value, {'id': 'instructor', 'role': 'admin'})
    assert jobs.read_json(parent / 'metadata.json')['status'] == 'queued'
    monkeypatch.setattr(auth, 'accessible_project_ids', lambda user: [])
    with pytest.raises(PermissionError):
        dependencies.request_cancellation(value, {'id': 'instructor', 'role': 'user'})


def test_idle_dispatch_does_not_scan_dependency_history(cluster, monkeypatch):
    store, nodes, user = cluster
    publish(store, nodes[0], 'root', status='cancelled')
    monkeypatch.setattr(store, '_reconcile_dependency_cancellations', lambda *args: pytest.fail('Idle dispatch scanned historical dependencies'))
    store.dispatch()


def test_disposable_demo_cancels_only_chain_and_keeps_completed_data(tmp_path, monkeypatch):
    from mn_molecular_design_studio.project_lifecycle_demo import demo_workspace, add_cancellation_demo, cancel_simulated_job
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(tmp_path / 'unrelated'))
    user = {'id': 'instructor', 'role': 'admin'}
    node, coordinator, reader = demo_workspace(tmp_path / 'sandbox', user)
    project_id = add_cancellation_demo(node.home, user)
    rows = [value for value in reader() if value['project_id'] == project_id]
    selected = next(value for value in rows if value['tool'] == 'demo_design')
    result = cancel_simulated_job(node, user, selected)
    assert result['downstream_count'] == 2
    current = {value['tool']: value['status'] for value in reader() if value['project_id'] == project_id}
    assert current == {'demo_design': 'cancelled', 'demo_refolding': 'cancelled', 'demo_analysis': 'cancelled',
                       'demo_completed': 'completed', 'demo_independent': 'queued'}
    assert not list(node.runs.glob('*/*/worker_request.json'))
