from threading import Event
from mn_molecular_design_studio import molecular_database_ingest as ingest
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio import molecule_imports


def test_library_refresh_returns_saved_records_while_one_indexer_runs(tmp_path, monkeypatch):
    database = MolecularDatabase(tmp_path / 'database')
    user = {'id': 'background-test', 'role': 'user'}
    entered, release, finished = Event(), Event(), Event()
    calls = []
    key = (str(database.path), user['id'])

    def blocking_refresh(actor, *, force, database):
        calls.append(actor['id'])
        entered.set()
        assert release.wait(5)
        ingest._LAST_SYNC_ERRORS[key] = [{'run_id': 'unreadable-source', 'error': 'Source unavailable'}]
        finished.set()

    monkeypatch.setattr(ingest, 'sync_workspace', blocking_refresh)
    try:
        assert ingest.sync_workspace_async(user, database=database) == []
        assert entered.wait(2)
        assert ingest.workspace_sync_running(user, database=database)
        assert ingest.sync_workspace_async(user, database=database, force=True) == []
        assert calls == [user['id']]
    finally:
        release.set()
        assert finished.wait(2)
    # Wait on the scheduler lock until the worker's finally block completes.
    import time
    deadline = time.monotonic() + 2
    while ingest.workspace_sync_running(user, database=database) and time.monotonic() < deadline:
        time.sleep(.005)
    assert not ingest.workspace_sync_running(user, database=database)
    assert ingest.sync_workspace_async(user, database=database) == [{'run_id': 'unreadable-source', 'error': 'Source unavailable'}]
    assert calls == [user['id']]


def test_identical_job_refresh_preserves_cached_evaluation_but_changed_scores_invalidate_it(tmp_path):
    database = MolecularDatabase(tmp_path / 'database')
    user = {'id': 'stable-index-test', 'role': 'admin'}
    imported = database.import_records(user, molecule_imports.parse_sequences('ACDEFGHIK'))
    molecule, group = imported['molecule_ids'][0], imported['group_id']
    arguments = dict(source_key='stable-observation', created_at='2026-11-04T12:00:00Z')
    observation = database.add_observation(group, molecule, 'test', 'run', {'score': .9}, **arguments)
    database.list_molecules(user)
    with database.connection() as connection:
        cached = tuple(connection.execute('SELECT * FROM molecule_evaluations WHERE observation_id=?', (observation,)).fetchone())
    database.add_observation(group, molecule, 'test', 'run', {'score': .9}, **arguments)
    with database.connection() as connection:
        assert tuple(connection.execute('SELECT * FROM molecule_evaluations WHERE observation_id=?', (observation,)).fetchone()) == cached
    database.add_observation(group, molecule, 'test', 'run', {'score': .5}, **arguments)
    with database.connection() as connection:
        assert connection.execute('SELECT 1 FROM molecule_evaluations WHERE observation_id=?', (observation,)).fetchone() is None
    assert database.list_molecules(user)['items'][0]['scores']['test.score'] == .5
