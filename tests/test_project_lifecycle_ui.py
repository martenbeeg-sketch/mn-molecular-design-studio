import json

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import project_lifecycle_ui
from mn_molecular_design_studio.project_lifecycle import archive_project
from mn_molecular_design_studio.project_lifecycle_demo import demo_workspace
from mn_molecular_design_studio.projects import project_inventory


@pytest.fixture
def demo(tmp_path, monkeypatch):
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(tmp_path / 'unrelated-runtime'))
    owner = {'id': 'instructor', 'role': 'admin'}
    node, coordinator, jobs = demo_workspace(tmp_path / 'disposable', owner)
    def item():
        projects = json.loads((node.home / 'organization/organization.json').read_text())['projects']
        project = next(row for row in projects if row['name'].startswith('Disposable cleanup demo'))
        return project_inventory([project], jobs(), node.database.groups(owner, scope='all'))[project['id']]
    monkeypatch.setattr(project_lifecycle_ui, '_demo_item', item, raising=False)
    monkeypatch.setattr(project_lifecycle_ui, '_demo_coordinator', coordinator, raising=False)
    monkeypatch.setattr(project_lifecycle_ui, '_demo_archive', lambda project_id, archived: archive_project(owner, project_id, archived, home=node.home), raising=False)
    return node, coordinator, jobs, item


SOURCE = '''
from mn_molecular_design_studio import project_lifecycle_ui as ui
ui.render_actions({'id': 'instructor', 'role': 'admin'}, ui._demo_item(),
    coordinator=ui._demo_coordinator, on_archive=ui._demo_archive)
'''


def test_demo_cleanup_requires_name_backup_and_explicit_collection_scope(demo):
    node, coordinator, jobs, item = demo
    app = AppTest.from_string(SOURCE).run()
    assert not app.exception
    app.radio[0].set_value('Delete whole project').run()
    app.button(key='projects-prepare-cleanup-' + item()['project']['id']).click().run()
    assert not app.exception
    assert app.error and not any(button.label == 'Confirm project deletion' for button in app.button)
    app.checkbox(key='projects-cleanup-collections-' + item()['project']['id']).set_value(True).run()
    app.button(key='projects-prepare-cleanup-' + item()['project']['id']).click().run()
    assert not app.exception
    confirm = next(button for button in app.button if button.label == 'Confirm project deletion')
    assert confirm.disabled
    app.text_input[0].set_value(item()['project']['name']).run()
    assert next(button for button in app.button if button.label == 'Confirm project deletion').disabled
    next(widget for widget in app.checkbox if 'recovery backup' in widget.label).set_value(True).run()
    assert not next(button for button in app.button if button.label == 'Confirm project deletion').disabled
    next(button for button in app.button if button.label == 'Confirm project deletion').click().run()
    assert not app.exception
    assert any('completed on every workstation' in entry.value for entry in app.success)
    assert not jobs()
    with node.database.connection() as connection:
        assert connection.execute('SELECT count(*) FROM molecule_members').fetchone()[0] == 1


def test_demo_archive_and_restore_do_not_remove_data(demo):
    node, coordinator, jobs, item = demo
    app = AppTest.from_string(SOURCE).run()
    next(button for button in app.button if button.label == 'Archive project').click().run()
    assert not app.exception and item()['project']['archived_at']
    assert len(jobs()) == 2
    next(button for button in app.button if button.label == 'Restore archived project').click().run()
    assert not item()['project'].get('archived_at')


def test_non_admin_can_read_project_actions_without_operation_access(demo):
    app = AppTest.from_string(SOURCE.replace("'instructor', 'role': 'admin'", "'student', 'role': 'user'")).run()
    assert not app.exception
    assert next(button for button in app.button if button.label == 'Prepare cleanup review').disabled
    assert any('administrator' in entry.value for entry in app.info)
