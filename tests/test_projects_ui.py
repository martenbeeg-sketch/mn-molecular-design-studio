from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import projects_ui


PROJECTS = [{"id": "p", "name": "Important research", "description": "Keep imported data"},
            {"id": "q", "name": "Smoke comparison", "description": "Not automatic deletion"}]
SCRIPT = '''
from mn_molecular_design_studio import projects_ui as ui
ui.render_projects({"id": "scientist", "role": "user"}, ui._test_projects,
    load_inventory=ui._test_load, on_create=lambda: None,
    on_results=ui._test_results, on_database=ui._test_database, read_only=True)
'''


def setup_ui(monkeypatch):
    calls = []
    monkeypatch.setattr(projects_ui, '_test_projects', PROJECTS, raising=False)
    monkeypatch.setattr(projects_ui, '_test_load', lambda: calls.append('load') or {
        'jobs': [{"run_id": "one", "task_group": "design", "project_id": "p", "status": "completed"}],
        'collections': [{"id": "g", "name": "Imported research", "project_id": "p", "origin": {"type": "file"}}],
    }, raising=False)
    monkeypatch.setattr(projects_ui, '_test_results', lambda project: calls.append(('results', project)), raising=False)
    monkeypatch.setattr(projects_ui, '_test_database', lambda project: calls.append(('database', project)), raising=False)
    return calls


def test_projects_navigation_cache_and_safe_cleanup_review(monkeypatch):
    calls = setup_ui(monkeypatch)
    app = AppTest.from_string(SCRIPT).run()
    assert not app.exception
    assert next(button for button in app.button if button.label == 'New project').disabled
    assert [metric.value for metric in app.metric] == ['1', '1', '0', '1']
    assert any('nothing is deleted without confirmation' in entry.value for entry in app.markdown)
    assert not any('Confirm' in button.label or 'purge' in button.label for button in app.button)
    app.button(key='projects-results').click().run()
    assert calls == ['load', ('results', 'p')]
    app.button(key='projects-database').click().run()
    assert calls[-1] == ('database', 'p')
    app.button(key='projects-refresh').click().run()
    assert calls.count('load') == 2
    app.selectbox(key='projects-selected').set_value('q').run()
    assert [metric.value for metric in app.metric] == ['0', '0', '0', '0']
    app.button(key='projects-cleanup-results').click().run()
    assert calls[-1] == ('results', 'q')


def test_project_pending_admission_distinguishes_results_from_library(monkeypatch):
    setup_ui(monkeypatch)
    monkeypatch.setattr(projects_ui, '_test_projects', [dict(PROJECTS[0], library_admission_review={'status': 'pending'})])
    app = AppTest.from_string(SCRIPT).run()
    assert not app.exception
    assert any('pending admission review' in item.value and 'do not mean a design passed' in item.value for item in app.info)


def test_project_search_and_access_scope_invalidate_inventory(monkeypatch):
    calls = setup_ui(monkeypatch)
    app = AppTest.from_string(SCRIPT).run()
    app.text_input[0].set_value('Smoke').run()
    assert app.selectbox(key='projects-selected').options == ['Smoke comparison · 0 runs']
    monkeypatch.setattr(projects_ui, '_test_projects', [PROJECTS[0]])
    app.text_input[0].set_value('').run()
    assert not app.exception
    assert app.selectbox(key='projects-selected').options == ['Important research · 1 runs']
    assert calls.count('load') == 2


def test_inventory_failure_never_looks_like_an_empty_project(monkeypatch):
    setup_ui(monkeypatch)
    def unavailable():
        raise OSError('Workstation unavailable')
    monkeypatch.setattr(projects_ui, '_test_load', unavailable)
    app = AppTest.from_string(SCRIPT).run()
    assert not app.exception
    assert app.error
    assert not app.metric
    assert not any(button.label == 'Review and select runs in My results' for button in app.button)


def test_sidebar_order_and_navigation_preserve_account_controls():
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.navigation import render_workspace_navigation
with st.sidebar:
    render_workspace_navigation({"id": "admin", "role": "admin"})
''').run()
    assert [button.label for button in app.sidebar.button] == [
        'All tools', 'Projects', 'Pipelines', 'My results', 'Database', 'Data browser', 'Complex Analysis', 'Users & access', 'Account',
    ]
    assert len(app.sidebar.get('divider')) == 1
    app.button(key='nav-projects').click().run()
    assert app.session_state['page'] == 'Projects'
    app.session_state['selected_workflow'] = 'previous-tool'
    app.session_state['selected_workflow_group'] = 'previous-group'
    app.button(key='nav-all-tools').click().run()
    assert app.session_state['page'] == 'All tools'
    assert 'selected_workflow' not in app.session_state
    assert 'selected_workflow_group' not in app.session_state


def test_project_run_pagination_and_clamping_after_refresh(monkeypatch):
    setup_ui(monkeypatch)
    jobs = [{"task_group": "design", "run_id": str(index), "project_id": "p", "status": "completed"}
            for index in range(53)]
    monkeypatch.setattr(projects_ui, '_test_load', lambda: {'jobs': list(jobs), 'collections': []})
    app = AppTest.from_string(SCRIPT).run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 25
    app.number_input(key='projects-runs-page-p').set_value(3).run()
    assert len(app.dataframe[0].value) == 3
    del jobs[26:]
    app.button(key='projects-refresh').click().run()
    assert not app.exception
    assert app.number_input(key='projects-runs-page-p').value == 2
    assert len(app.dataframe[0].value) == 1


def test_project_handoff_scopes_shared_results_and_drops_stale_purge(monkeypatch, tmp_path):
    from mn_molecular_design_studio import workspace
    from mn_molecular_design_studio.cluster import ui

    setup_ui(monkeypatch)
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(tmp_path / 'studio'))
    monkeypatch.setattr(workspace, 'load_organization', lambda: {'projects': PROJECTS})
    jobs = [{'cluster_job_id': project, 'node_id': 'remote', 'project_id': project,
             'campaign_name': 'Job ' + project, 'status': 'completed'} for project in ('p', 'q')]
    monkeypatch.setattr(ui, 'catalog', lambda *args, **kwargs: ({'items': jobs, 'total': 2}, False))
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio import projects_ui as ui
from mn_molecular_design_studio.navigation import open_project_results
from mn_molecular_design_studio.cluster.ui import render_runs
def results(project_id):
    open_project_results(st.session_state, project_id)
    st.rerun()
if st.session_state.get('page') == 'My runs':
    render_runs({'id': 'scientist'}, 'token', lambda **kwargs: {'items': [], 'total': 0})
else:
    ui.render_projects({'id': 'scientist'}, ui._test_projects, load_inventory=ui._test_load,
        on_create=lambda: None, on_results=results, on_database=lambda project_id: None)
''').run()
    app.session_state['cluster-results-selected'] = {'q': jobs[1]}
    app.session_state['cluster-results-pending-remove-and-purge'] = [jobs[1]]
    app.button(key='projects-cleanup-results').click().run()
    assert not app.exception
    assert app.selectbox(key='cluster-results-project').value == 'p'
    assert list(app.dataframe[0].value.index) == ['p']
    assert not app.session_state['cluster-results-selected']
    assert not any(button.label == 'Confirm removal and purge' for button in app.button)


def test_connected_execution_panel_and_cancellation_preserve_upstream_and_scope(monkeypatch):
    calls = []
    jobs = [{'node_id': 'node', 'task_group': 'design', 'run_id': name, 'project_id': 'p', 'owner_id': 'scientist',
             'status': 'queued', 'campaign_name': name} for name in ('design', 'analysis', 'unrelated')]
    jobs[1].update(upstream_task_group='design', upstream_run_id='design')
    monkeypatch.setattr(projects_ui, '_test_projects', [PROJECTS[0]], raising=False)
    monkeypatch.setattr(projects_ui, '_test_load', lambda: {'jobs': jobs, 'collections': []}, raising=False)
    monkeypatch.setattr(projects_ui, '_test_cancel', lambda row: calls.append(row['run_id']) or {'downstream_count': 1}, raising=False)
    app = AppTest.from_string('''
from mn_molecular_design_studio import projects_ui as ui
ui.render_projects({'id': 'scientist'}, ui._test_projects, load_inventory=ui._test_load,
    on_create=lambda: None, on_results=lambda project: None, on_database=lambda project: None,
    on_cancel=ui._test_cancel)
''').run()
    assert not app.exception
    assert app.metric[0].value == '2'
    from mn_molecular_design_studio.projects import project_inventory
    connected = next(run for run in project_inventory(PROJECTS[:1], jobs, [])['p']['logical_runs'] if run['job_count'] == 2)
    app.selectbox(key='projects-inspect-run-p').set_value(connected['id']).run()
    assert 'Depends on' in app.dataframe[1].value.columns
    next(button for button in app.button if button.label == 'Cancel selected job + downstream').click().run()
    assert not app.exception and len(calls) == 1
    assert any('1 recorded downstream' in entry.value for entry in app.success)
