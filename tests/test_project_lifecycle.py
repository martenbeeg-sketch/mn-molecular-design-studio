import json
import os
from pathlib import Path
import shutil
import sqlite3
import uuid

import pytest

from mn_molecular_design_studio import molecule_imports, workspace
from mn_molecular_design_studio.cluster import config
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.project_lifecycle import (
    ProjectCoordinator, ProjectNode, ProjectUnavailableError, archive_project,
)


ADMIN = {"id": "instructor", "role": "admin"}


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    home = tmp_path / 'studio'
    runs = tmp_path / 'runs'
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(home))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_APP_HOME', str(home))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR', str(runs))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR', str(tmp_path / 'reference'))
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    project = workspace.create_project('Important research', home=home, owner_id=ADMIN['id'])
    database = MolecularDatabase(home / 'database')

    def jobs():
        result = []
        for path in runs.glob('*/*/metadata.json'):
            record = json.loads(path.read_text())
            result.append({**record, 'task_group': path.parent.parent.name, 'run_id': path.parent.name})
        return result

    node = ProjectNode(home=home, runs=runs, database=database, node_id='local', job_reader=jobs)
    coordinator = ProjectCoordinator(node, peers=[])

    def run(name='one', status='completed'):
        directory = runs / 'design' / name
        directory.mkdir(parents=True)
        (directory / 'metadata.json').write_text(json.dumps({'project_id': project['id'], 'owner_id': ADMIN['id'], 'status': status}))
        (directory / 'result.txt').write_text('retained scientific result')
        return directory

    def group(name='Automatic results', automatic=True, reference='one'):
        imported = database.import_records(ADMIN, molecule_imports.parse_sequences('>design\nACDEFG'),
            name=name, project_id=project['id'])
        if automatic:
            with database.connection() as connection:
                connection.execute('UPDATE molecule_groups SET origin=? WHERE id=?',
                    (json.dumps({'type': 'job', 'task_group': 'design', 'run_id': reference}), imported['group_id']))
        return imported

    return home, runs, project, database, node, coordinator, run, group


def confirm(coordinator, operation):
    return coordinator.resume(ADMIN, operation['id'], confirmation=operation['project_name'], backup_acknowledged=True)


def test_archive_is_reversible_and_blocks_new_work(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    directory = run()
    imported = group()
    archive_project(ADMIN, project['id'])
    with pytest.raises(ProjectUnavailableError):
        workspace.validate_project_id(project['id'])
    assert directory.is_dir()
    assert database.list_molecules(ADMIN, group_id=imported['group_id'])['total'] == 1
    archive_project(ADMIN, project['id'], False)
    assert workspace.validate_project_id(project['id']) == project['id']
    with pytest.raises(PermissionError):
        archive_project({'id': 'other', 'role': 'user'}, project['id'])
    from mn_molecular_design_studio.runtime import runs_root
    assert runs_root() == runs


def test_cleanup_recovery_contains_actual_coordinator_database(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    run()
    group()
    path = home / 'cluster/coordinator.sqlite3'
    path.parent.mkdir()
    with sqlite3.connect(path) as connection:
        connection.execute('CREATE TABLE retained_command (id TEXT PRIMARY KEY)')
        connection.execute("INSERT INTO retained_command VALUES ('remote-command')")
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one'])
    result = confirm(coordinator, operation)
    assert result['status'] == 'completed'
    recovery = runs.parent / '.project-recovery' / operation['id']
    with sqlite3.connect(recovery / 'cluster.sqlite3') as connection:
        assert connection.execute('SELECT id FROM retained_command').fetchone()[0] == 'remote-command'


def test_cleanup_moves_files_and_removes_only_automatic_memberships(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    directory = run()
    automatic = group()
    saved = group('Saved research', False)
    database.attach_structure(automatic['group_id'], automatic['molecule_ids'][0], {'name': 'model.pdb', 'content': b'HEADER model\nEND\n'})
    database.add_observation(automatic['group_id'], automatic['molecule_ids'][0], 'ipsae', 'one', {'score': .8}, context={'task_group': 'design', 'run_id': 'one'})
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one'])
    assert not operation['blocked']
    assert directory.is_dir()
    result = confirm(coordinator, operation)
    assert result['status'] == 'completed', result.get('error')
    assert not directory.exists()
    recovery = runs.parent / '.project-recovery' / operation['id']
    assert (recovery / 'runs/design/one/result.txt').read_text() == 'retained scientific result'
    assert database.list_molecules(ADMIN, group_id=automatic['group_id'])['total'] == 0
    assert database.list_molecules(ADMIN, group_id=saved['group_id'])['total'] == 1
    assert workspace.validate_project_id(project['id']) == project['id']
    with sqlite3.connect(recovery / 'molecular_library.sqlite3') as backup:
        assert backup.execute('SELECT count(*) FROM molecule_members WHERE group_id=?', (automatic['group_id'],)).fetchone()[0] == 1
        assert backup.execute('SELECT count(*) FROM molecule_observations').fetchone()[0] == 1
    assert confirm(coordinator, operation)['status'] == 'completed'


def test_deletion_requires_explicit_saved_collection_inclusion_and_keeps_shared_identity(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    run()
    saved = group('Imported important research', False)
    other = database.import_records(ADMIN, molecule_imports.parse_sequences('ACDEFG'), name='Other project')
    blocked = coordinator.prepare(ADMIN, project['id'], mode='delete')
    assert blocked['blocked']
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    result = confirm(coordinator, operation)
    assert result['status'] == 'completed', result.get('error')
    assert next(row for row in workspace.load_organization()['projects'] if row['id'] == project['id'])['deleted_at']
    with pytest.raises(ProjectUnavailableError):
        workspace.validate_project_id(project['id'])
    assert database.list_molecules(ADMIN, group_id=other['group_id'])['total'] == 1
    with database.connection() as connection:
        assert connection.execute('SELECT count(*) FROM molecule_members WHERE group_id=?', (saved['group_id'],)).fetchone()[0] == 0
        assert connection.execute('SELECT count(*) FROM molecules').fetchone()[0] == 1
        with pytest.raises(sqlite3.IntegrityError, match='locked|deleted'):
            connection.execute('INSERT INTO molecule_members VALUES (?,?,?,?,?,?)',
                (saved['group_id'], saved['molecule_ids'][0], 'stale replay', '[]', '{}', 'now'))


def test_active_run_wrong_confirmation_and_stale_review_never_remove_data(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    active = run(status='running')
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    assert operation['blocked']
    with pytest.raises(ValueError, match='blocking'):
        confirm(coordinator, operation)
    metadata = json.loads((active / 'metadata.json').read_text())
    metadata['status'] = 'completed'
    (active / 'metadata.json').write_text(json.dumps(metadata))
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    with pytest.raises(ValueError, match='exact project name'):
        coordinator.resume(ADMIN, operation['id'], confirmation='wrong', backup_acknowledged=True)
    with pytest.raises(ValueError, match='acknowledge'):
        coordinator.resume(ADMIN, operation['id'], confirmation=project['name'], backup_acknowledged=False)
    (active / 'result.txt').write_text('changed after review')
    result = confirm(coordinator, operation)
    assert result['status'] == 'aborted'
    assert active.is_dir()
    assert not (runs.parent / '.project-recovery').exists()


def test_retained_structure_in_another_collection_blocks_run_removal(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    directory = run()
    saved = group('Saved research', False)
    structure = directory / 'complex.pdb'
    structure.write_text('HEADER retained\n')
    database.attach_structure(saved['group_id'], saved['molecule_ids'][0], {'path': str(structure), 'name': 'complex.pdb'})
    with database.connection() as connection:
        connection.execute('UPDATE molecule_structures SET path=? WHERE group_id=?', (str(structure), saved['group_id']))
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one'])
    assert any('retained or shared' in message for message in operation['blocked'])
    assert directory.exists()


def test_sql_guards_block_peer_insert_and_cross_project_run_references(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    directory = run()
    automatic = group()
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one'])
    node.begin(ADMIN, operation['plans']['local']['id'], operation['id'])
    with database.connection() as connection:
        with pytest.raises(sqlite3.IntegrityError, match='locked'):
            connection.execute('UPDATE molecule_members SET name=? WHERE group_id=?', ('stale edit', automatic['group_id']))
    node.abort(ADMIN, operation['id'])
    assert workspace.validate_project_id(project['id']) == project['id']
    assert directory.is_dir()


def test_interruption_after_run_move_can_resume_and_recovery_is_restorable(fixture, monkeypatch, tmp_path):
    home, runs, project, database, node, coordinator, run, group = fixture
    directory = run()
    automatic = group()
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one'])
    original = node.save
    interrupted = []
    def crash(kind, key, value):
        if kind == 'nodes' and value.get('moved') and not interrupted:
            interrupted.append(True)
            raise OSError('simulated process interruption')
        return original(kind, key, value)
    monkeypatch.setattr(node, 'save', crash)
    result = confirm(coordinator, operation)
    assert result['status'] == 'attention'
    assert not directory.exists()
    monkeypatch.setattr(node, 'save', original)
    result = confirm(coordinator, operation)
    assert result['status'] == 'completed', result.get('error')
    recovery = runs.parent / '.project-recovery' / operation['id']
    restore = tmp_path / 'offline-restore-check'
    restore.mkdir()
    shutil.copy2(recovery / 'molecular_library.sqlite3', restore / 'library.sqlite3')
    shutil.copytree(recovery / 'runs', restore / 'runs')
    with sqlite3.connect(restore / 'library.sqlite3') as connection:
        assert connection.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        assert connection.execute('SELECT count(*) FROM molecule_members WHERE group_id=?', (automatic['group_id'],)).fetchone()[0] == 1
    assert (restore / 'runs/design/one/result.txt').read_text() == 'retained scientific result'


def test_unavailable_or_old_peer_blocks_preparation_before_any_lock(fixture, monkeypatch):
    home, runs, project, database, node, coordinator, run, group = fixture
    directory = run()
    peer = {'id': 'remote', 'url': 'http://offline.invalid'}
    original = config.load()
    monkeypatch.setattr(config, 'load', lambda: {**original, 'peers': [peer]})
    coordinator = ProjectCoordinator(node, peers=[peer], remote=lambda *args: (_ for _ in ()).throw(RuntimeError('peer does not implement protocol')))
    with pytest.raises(RuntimeError, match='protocol'):
        coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    assert directory.is_dir()
    assert workspace.validate_project_id(project['id']) == project['id']


def test_symlink_run_cannot_escape_results_or_delete_shared_files(fixture, tmp_path):
    home, runs, project, database, node, coordinator, run, group = fixture
    outside = tmp_path / 'important-files'
    outside.mkdir()
    (runs / 'design').mkdir(parents=True)
    (runs / 'design/unsafe').symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match='escapes'):
        node.run_path('design/unsafe')
    with pytest.raises(ValueError):
        node.run_path('../important-files')
    assert outside.is_dir()


def test_shared_sql_analysis_cohort_cannot_be_partially_removed(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    run('one')
    run('two')
    imported = group()
    molecule = imported['molecule_ids'][0]
    database.add_observation(imported['group_id'], molecule, 'biotite', 'one', {'score': 1}, context={'task_group': 'design', 'run_id': 'one'})
    database.add_observation(imported['group_id'], molecule, 'biotite', 'two', {'score': 2}, context={'task_group': 'design', 'run_id': 'two'})
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one'])
    assert any('share one SQL collection' in problem for problem in operation['blocked'])
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one', 'local:design/two'])
    assert not operation['blocked']
    assert confirm(coordinator, operation)['status'] == 'completed'


def test_unselected_downstream_job_and_legacy_lineage_block_cleanup(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    run('one')
    consumer = run('two')
    (consumer / 'input.json').write_text(json.dumps({'source_job': 'design/one'}))
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one'])
    assert any('downstream' in problem for problem in operation['blocked'])
    imported = group()
    database.add_observation(imported['group_id'], imported['molecule_ids'][0], 'unknown-analysis', 'legacy-run', {'score': 1})
    operation = coordinator.prepare(ADMIN, project['id'], mode='clean', selected_runs=['local:design/one', 'local:design/two'])
    assert any('lineage is incomplete' in problem for problem in operation['blocked'])


def test_concurrent_human_review_and_old_purge_are_blocked_by_cleanup_lock(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    run()
    imported = group()
    molecule = imported['molecule_ids'][0]
    database.attach_structure(imported['group_id'], molecule, {'name': 'model.pdb', 'content': b'HEADER model\n'})
    with database.connection() as connection:
        structure_id = connection.execute('SELECT id FROM molecule_structures').fetchone()[0]
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    node.begin(ADMIN, operation['plans']['local']['id'], operation['id'])
    with database.connection() as connection:
        with pytest.raises(sqlite3.IntegrityError, match='review is locked'):
            connection.execute('INSERT INTO complex_review_annotations VALUES (?,?,?,?,?,?,?,?,?,?)',
                ('review', 'database', molecule, structure_id, '', ADMIN['id'], 'Reviewed', 'concurrent edit', 'now', 'now'))
    with pytest.raises(PermissionError, match='locked'):
        database.group(imported['group_id'], ADMIN, write=True)
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import delete_job_run
    with pytest.raises(PermissionError, match='locked'):
        delete_job_run('design', 'one')


def test_recovery_corruption_stops_before_sql_deletion(fixture, monkeypatch):
    home, runs, project, database, node, coordinator, run, group = fixture
    run()
    imported = group()
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    node.begin(ADMIN, operation['plans']['local']['id'], operation['id'])
    original = node._backup
    def corrupted(operation, plan):
        destination = original(operation, plan)
        Path(operation['backup']['database']).write_bytes(b'corrupted backup')
        raise OSError('interruption before removing original data')
    monkeypatch.setattr(node, '_backup', corrupted)
    result = confirm(coordinator, operation)
    assert result['status'] == 'attention'
    monkeypatch.setattr(node, '_backup', original)
    result = confirm(coordinator, operation)
    assert result['status'] == 'attention' and 'checksum' in result['error']
    assert database.list_molecules(ADMIN, group_id=imported['group_id'])['total'] == 1


def test_two_workstations_remain_pending_until_both_confirm_and_retry_is_idempotent(fixture, tmp_path):
    home, runs, project, database, node, coordinator, run, group = fixture
    run()
    group()
    remote_home = tmp_path / 'remote'
    (remote_home / 'organization').mkdir(parents=True)
    shutil.copy2(home / 'organization/organization.json', remote_home / 'organization/organization.json')
    remote_database = MolecularDatabase(remote_home / 'database')
    remote_import = remote_database.import_records(ADMIN, molecule_imports.parse_sequences('GHIKL'), project_id=project['id'])
    remote_node = ProjectNode(home=remote_home, runs=remote_home / 'runs', database=remote_database,
                              node_id='remote', job_reader=lambda: [], participant_ids=['local', 'remote'])
    node.participant_ids = ['local', 'remote']
    unavailable = []
    counts = {}
    def remote(peer, method, payload):
        counts[method] = counts.get(method, 0) + 1
        if method == 'execute' and not unavailable:
            unavailable.append(True)
            raise OSError('remote workstation went offline')
        return getattr(remote_node, method)(ADMIN, **{key: value for key, value in payload.items() if key != 'actor_id'})
    coordinator = ProjectCoordinator(node, peers=[{'id': 'remote', 'url': 'fixture'}], remote=remote)
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    result = confirm(coordinator, operation)
    assert result['status'] == 'attention'
    assert not node.project(project['id']).get('deleted_at')
    assert remote_database.list_molecules(ADMIN, group_id=remote_import['group_id'])['total'] == 1
    result = confirm(coordinator, operation)
    assert result['status'] == 'completed', result.get('error')
    assert counts['begin'] == 1
    assert node.project(project['id'], allow_deleted=True)['deleted_at']
    assert remote_node.project(project['id'], allow_deleted=True)['deleted_at']
    with remote_database.connection() as connection:
        assert connection.execute('SELECT count(*) FROM molecule_members').fetchone()[0] == 0


def test_completion_receipt_can_recover_after_project_tombstone_write(fixture, monkeypatch):
    home, runs, project, database, node, coordinator, run, group = fixture
    run()
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    original = node.save
    interrupted = []
    def crash(kind, key, value):
        if kind == 'nodes' and value.get('status') == 'completed' and not interrupted:
            interrupted.append(True)
            raise OSError('interruption after project tombstone')
        return original(kind, key, value)
    monkeypatch.setattr(node, 'save', crash)
    assert confirm(coordinator, operation)['status'] == 'attention'
    monkeypatch.setattr(node, 'save', original)
    assert coordinator.pending(ADMIN)
    assert confirm(coordinator, operation)['status'] == 'completed'


def test_old_project_metadata_cannot_resurrect_deleted_identity(fixture):
    home, runs, project, database, node, coordinator, run, group = fixture
    from mn_molecular_design_studio.cluster.replication import Replica
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    assert confirm(coordinator, operation)['status'] == 'completed'
    with pytest.raises(ValueError, match='resurrected'):
        Replica()._apply('project', 'project:' + project['id'], project, 'old-peer')


def test_cleanup_detaches_hardlinks_to_preserve_recovery_bytes(fixture, tmp_path):
    home, runs, project, database, node, coordinator, run, group = fixture
    directory = run()
    external = tmp_path / 'shared-result.txt'
    external.write_text('original shared bytes')
    os.link(external, directory / 'linked.txt')
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    assert confirm(coordinator, operation)['status'] == 'completed'
    external.write_text('later external change')
    recovery = runs.parent / '.project-recovery' / operation['id']
    assert (recovery / 'runs/design/one/linked.txt').read_text() == 'original shared bytes'


def test_unindexed_owned_run_and_conflicting_marker_block_project_deletion(fixture, monkeypatch):
    home, runs, project, database, node, coordinator, run, group = fixture
    directory = run()
    monkeypatch.setattr(node, 'job_reader', lambda: [])
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    assert any('not reconciled' in problem for problem in operation['blocked'])
    monkeypatch.setattr(node, 'job_reader', lambda: [{'project_id': project['id'], 'task_group': 'design', 'run_id': 'one', 'status': 'completed'}])
    (directory / 'molecular_design_studio.json').write_text(json.dumps({'project_id': 'other-project'}))
    operation = coordinator.prepare(ADMIN, project['id'], mode='delete', include_collections=True)
    assert any('disagree' in problem for problem in operation['blocked'])
