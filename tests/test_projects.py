import copy
import json

from mn_molecular_design_studio.navigation import open_project_database, open_project_results
from mn_molecular_design_studio.projects import cleanup_plan, load_project_jobs, project_inventory


def test_project_jobs_page_every_execution_include_hidden_and_deduplicate_local():
    rows = [{"task_group": "design", "run_id": str(index), "project_id": "p", "status": "completed"}
            for index in range(1003)]
    calls = []
    def local_search(**options):
        calls.append(options)
        offset = options.get('offset', 0)
        return {"items": rows[offset:offset + options['limit']], "total": len(rows)}
    remote = {"items": [{**rows[0], "node_id": "local", "cluster_job_id": "linked"},
                        {**rows[0], "node_id": "remote"}], "total": 2}
    snapshot = load_project_jobs(local_search, remote_catalog=lambda: (remote, False), local_node_id='local')
    assert len(snapshot['jobs']) == 1004
    assert all(call['include_hidden'] for call in calls)
    assert next(row for row in snapshot['jobs'] if row['node_id'] == 'local' and row['run_id'] == '0')['cluster_job_id'] == 'linked'
    assert snapshot['stale'] is False


def test_incomplete_catalog_is_not_reported_as_complete():
    snapshot = load_project_jobs(lambda **options: {"items": [], "total": 1})
    assert snapshot['stale'] is True


def test_project_inventory_respects_scope_qualified_saved_run_membership(monkeypatch, tmp_path):
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(tmp_path))
    rows = [{"task_group": "design", "run_id": str(index), "node_id": "local",
             "project_id": "p", "owner_id": "scientist", "status": "completed"} for index in range(2)]
    (tmp_path / 'run_membership.json').write_text(json.dumps({'memberships': {
        'local:design/' + row['run_id']: {'project_id': 'p', 'owner_id': 'scientist',
                                       'tags': {'studio_run_id': 'logical-one', 'studio_run_name': 'One campaign'}}
        for row in rows
    }}))
    snapshot = load_project_jobs(lambda **options: {'items': rows, 'total': 2}, local_node_id='local')
    item = project_inventory([{'id': 'p', 'name': 'Project'}], snapshot['jobs'], [])['p']
    assert item['runs'] == 1
    assert item['logical_runs'][0]['name'] == 'One campaign'
    assert item['logical_runs'][0]['job_count'] == 2
    snapshot = load_project_jobs(lambda **options: {"items": [], "total": 0},
                                 remote_catalog=lambda: ({"items": [], "total": 1}, False))
    assert snapshot['stale'] is True


def test_project_inventory_preserves_machine_identity_and_filters_projects():
    project = {"id": "research", "name": "Smoke is just a name"}
    base = {"project_id": "research", "task_group": "design", "run_id": "same", "status": "completed"}
    jobs = [{**base, "node_id": "alpha"}, {**base, "node_id": "beta"},
            {**base, "node_id": "alpha"}, {**base, "project_id": "other", "node_id": "alpha"}]
    before = copy.deepcopy(jobs)
    item = project_inventory([project], jobs, [])['research']
    assert len(item['jobs']) == 2
    assert len(item['finished']) == 2
    assert jobs == before
    plan = cleanup_plan(item)
    assert plan['destructive'] is False
    assert plan['policy']['automatic_selection_or_deletion'] is False
    assert {row['node_id'] for row in plan['finished_executions_to_review']} == {'alpha', 'beta'}


def test_cleanup_keeps_imports_and_manual_collections_and_blocks_unknown_states():
    jobs = [{"project_id": "p", "run_id": str(index), "status": status} for index, status in enumerate(
        ["completed", "FAILED", "cancelled", "running", "queued", "paused", "unknown", None]
    )]
    collections = [
        {"id": "auto", "name": "run", "project_id": "p", "origin": {"type": "job"}},
        {"id": "import", "name": "Imported research", "project_id": "p", "origin": {"type": "file"}},
        {"id": "saved", "name": "Later selection", "project_id": "p", "origin": {"type": "results_design_collection"}},
        {"id": "legacy", "name": "Unknown origin", "project_id": "p", "origin": {}},
        {"id": "unassigned", "name": "Legacy unassigned", "origin": {}},
    ]
    item = project_inventory([{"id": "p", "name": "Project"}], jobs, collections)['p']
    plan = cleanup_plan(item)
    assert len(plan['finished_executions_to_review']) == 3
    assert len(plan['blocked_executions']) == 5
    assert {row['id'] for row in plan['protected_collections']} == {'import', 'saved', 'legacy'}
    assert len(item['collections']) == 4
    assert plan['policy']['preserve_project_metadata'] is True


def test_hidden_finished_run_still_counts_for_cleanup_without_being_unhidden():
    jobs = [{"project_id": "p", "task_group": "design", "run_id": "old", "hidden": True, "status": "completed"}]
    item = project_inventory([{"id": "p", "name": "Project"}], jobs, [])['p']
    assert item['runs'] == 1
    assert item['jobs'][0]['hidden'] is True
    assert cleanup_plan(item)['finished_executions_to_review'][0]['hidden'] is True


def test_project_results_handoff_clears_old_selection_confirmation_and_filters():
    state = {"cluster-results-pending-remove-and-purge": [{"project_id": "wrong"}],
             "results-pending-remove-and-purge": [{"project_id": "wrong"}],
             "cluster-results-show-removed": True,
             "cluster-results-selected": {"other": {}}, "results-selected-jobs": {"other": {}},
             "cluster-results-machine": "remote", "cluster-results-query": "old search",
             "results-row-select-wrong": True, "cluster-results-purge-operations": [{"command_id": "pending"}]}
    open_project_results(state, 'p')
    assert state['page'] == 'My runs'
    assert state['results-project-filter'] == state['cluster-results-project'] == 'p'
    assert not state['cluster-results-selected'] and not state['results-selected-jobs']
    assert state['cluster-results-machine'] == state['cluster-results-query'] == ''
    assert 'cluster-results-pending-remove-and-purge' not in state
    assert 'results-pending-remove-and-purge' not in state
    assert state['cluster-results-show-removed'] is False
    assert 'results-row-select-wrong' not in state
    assert state['cluster-results-purge-operations'] == [{"command_id": "pending"}]


def test_project_database_handoff_discards_previous_collection_preselection():
    cache = object()
    state = {"molecular-preselect-group": "wrong", "molecular-preselect-molecule": "wrong",
             "molecular-scope": "Mine", "database-browser-query": "wrong search",
             "database-browser-result-mode": "Calculation history", "database-results-query-previous": "wrong result",
             "molecular-method-user_uploads": False, "database-browser-cache": cache}
    open_project_database(state, 'p')
    assert state['molecular-project-filter'] == 'p'
    assert state['database-browser-view'] == 'Molecules'
    assert not any(key.startswith('molecular-preselect-') for key in state)
    assert state['molecular-scope'] == 'All'
    assert state['molecular-collection-view'] == 'Library'
    assert 'database-browser-query' not in state and 'database-browser-result-mode' not in state
    assert 'database-results-query-previous' not in state
    assert state['molecular-method-user_uploads'] is True
    assert state['database-browser-cache'] is cache
