"""Durable metadata replication between trusted workstations.

Each record has a Lamport version; clock skew cannot make an old record win.
Execution records belong to the machine holding the run. Sessions and scientific
files are never replicated. Account password hashes are private peer metadata.
"""
from __future__ import annotations

import json
import fcntl
import logging
from functools import wraps
import sqlite3
import threading
import uuid
from contextlib import contextmanager
from contextvars import ContextVar

from . import config
from .store import Store, encoded

_mutex = threading.RLock()
_depth = ContextVar("studio_replication_lock", default=0)
_importing = ContextVar("studio_replication_import", default=False)
IDENTITY_KINDS = {"user", "membership", "project"}
USER_FIELDS = ("id", "username", "display_name", "role", "password_hash", "must_change_password", "disabled", "created_at", "updated_at")


@contextmanager
def locked():
    with _mutex:
        if _depth.get():
            yield
            return
        with (config.home() / "replication.lock").open("a") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            token = _depth.set(1)
            try:
                yield
            finally:
                _depth.reset(token)
                fcntl.flock(handle, fcntl.LOCK_UN)


def identity_edit(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        if not config.is_peer():
            return function(*args, **kwargs)
        with locked():
            result = function(*args, **kwargs)
            capture_identity()
            return result
    return wrapped


class Replica:
    def __init__(self):
        self.path = config.home() / "replication.sqlite3"
        with self.db() as db:
            db.executescript("""
                PRAGMA journal_mode=WAL;
                CREATE TABLE IF NOT EXISTS records (
                    key TEXT PRIMARY KEY, kind TEXT NOT NULL, author TEXT NOT NULL,
                    counter INTEGER NOT NULL, value TEXT NOT NULL, seq INTEGER NOT NULL);
                CREATE INDEX IF NOT EXISTS records_sequence ON records(seq);
                CREATE TABLE IF NOT EXISTS state (key TEXT PRIMARY KEY, value TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS conflicts (key TEXT PRIMARY KEY, record TEXT NOT NULL, error TEXT NOT NULL);
            """)
            db.execute("INSERT OR IGNORE INTO state VALUES ('epoch',?)", (str(uuid.uuid4()),))
        self.path.chmod(0o600)

    @contextmanager
    def db(self):
        db = sqlite3.connect(self.path, timeout=20)
        db.row_factory = sqlite3.Row
        try:
            with db:
                yield db
        finally:
            db.close()

    def state(self, key, default=None):
        with self.db() as db:
            row = db.execute("SELECT value FROM state WHERE key=?", (key,)).fetchone()
        return json.loads(row[0]) if row else default

    def set_state(self, key, value):
        with self.db() as db:
            db.execute("INSERT OR REPLACE INTO state VALUES (?,?)", (key, encoded(value)))

    @property
    def epoch(self):
        with self.db() as db:
            return db.execute("SELECT value FROM state WHERE key='epoch'").fetchone()[0]

    def _write(self, db, key, kind, value, *, author=None, counter=None):
        maximum = db.execute("SELECT coalesce(max(counter),0) FROM records").fetchone()[0]
        clock = db.execute("SELECT value FROM state WHERE key='clock'").fetchone()
        maximum = max(maximum, int(clock[0]) if clock else 0)
        seq = db.execute("SELECT coalesce(max(seq),0)+1 FROM records").fetchone()[0]
        author = author or config.identity()
        counter = counter if counter is not None else maximum + 1
        db.execute("INSERT OR REPLACE INTO state VALUES ('clock',?)", (str(max(maximum, counter)),))
        db.execute("INSERT OR REPLACE INTO records VALUES (?,?,?,?,?,?)", (key, kind, author, counter, encoded(value), seq))

    def capture(self, kind, values):
        with self.db() as db:
            existing = {r['key']: r for r in db.execute("SELECT * FROM records WHERE kind=?", (kind,))}
            for key in sorted(set(values) | set(existing)):
                value = values.get(key)
                if key in existing and json.loads(existing[key]['value']) == value:
                    continue
                # A failed import must remain a conflict, not become a local deletion.
                if key not in values and db.execute("SELECT 1 FROM conflicts WHERE key=?", (key,)).fetchone():
                    continue
                self._write(db, key, kind, value)

    def capture_jobs(self):
        own = config.identity()
        with Store().connect() as db:
            rows = db.execute("SELECT j.*,s.origin_node FROM jobs j LEFT JOIN submissions s ON s.job_id=j.id WHERE (j.node_id=? AND j.local_ref IS NOT NULL) OR (s.origin_node=? AND j.local_ref IS NULL)", (own, own)).fetchall()
        current = {}
        for row in rows:
            value = json.loads(row['payload'])
            if not row['local_ref']:
                value['submission_origin'] = own
            current['job:' + row['id']] = value
        with self.db() as db:
            previous = db.execute("SELECT key,value FROM records WHERE kind='job' AND author=?", (own,)).fetchall()
            for record in previous:
                if record['key'] not in current and json.loads(record['value']) is not None:
                    self._write(db, record['key'], 'job', None)
            for key, value in current.items():
                previous = db.execute("SELECT value FROM records WHERE key=?", (key,)).fetchone()
                if not previous or json.loads(previous[0]) != value:
                    self._write(db, key, 'job', value)

    def export(self, since=0, limit=500):
        if type(since) is not int or since < 0 or type(limit) is not int or not 1 <= limit <= 500:
            raise ValueError("Invalid synchronization cursor.")
        with self.db() as db:
            rows = db.execute("SELECT * FROM records WHERE seq>? ORDER BY seq LIMIT ?", (since, limit)).fetchall()
        records = [dict(key=r['key'], kind=r['kind'], author=r['author'], counter=r['counter'], value=json.loads(r['value'])) for r in rows]
        return {'epoch': self.epoch, 'records': records, 'cursor': rows[-1]['seq'] if rows else since, 'more': len(rows) == limit}

    def receive(self, records):
        if not isinstance(records, list) or len(records) > 500:
            raise ValueError("Synchronize at most 500 metadata records per request.")
        known = {n['id'] for n in Store().nodes()} | {config.identity()}
        accepted = 0
        with locked():
            capture_identity()
            for record in records:
                if not isinstance(record, dict) or set(record) != {'key','kind','author','counter','value'}:
                    raise ValueError("Invalid peer metadata record.")
                key, kind, author, counter, value = (record[k] for k in ('key','kind','author','counter','value'))
                if author not in known or kind not in IDENTITY_KINDS | {'job', 'molecule_group', 'molecule_entry'} or not isinstance(key, str) or not key.startswith(kind + ':') or len(key) > 300 or type(counter) is not int or not 0 < counter < 2**62:
                    raise ValueError("Invalid metadata authority or version.")
                if len(encoded(record)) > 2 * 1024 * 1024:
                    raise ValueError("Metadata record exceeds 2 MB.")
                with self.db() as db:
                    previous = db.execute("SELECT * FROM records WHERE key=?", (key,)).fetchone()
                    old = json.loads(previous['value']) if previous else None
                    # Once accepted, execution-machine ownership outranks origin placeholders.
                    promote = kind == 'job' and value and value.get('local_reference') and not (old or {}).get('local_reference')
                    if kind == 'job' and (old or {}).get('local_reference') and not (value or {}).get('local_reference'):
                        continue
                    if previous and not promote and (counter, author) <= (previous['counter'], previous['author']):
                        continue
                    try:
                        token = _importing.set(True)
                        try:
                            self._apply(kind, key, value, author)
                        finally:
                            _importing.reset(token)
                    except (ValueError, PermissionError, sqlite3.IntegrityError) as exc:
                        db.execute("INSERT OR REPLACE INTO conflicts VALUES (?,?,?)", (key, encoded(record), str(exc)))
                        continue
                    self._write(db, key, kind, value, author=author, counter=counter)
                    db.execute("DELETE FROM conflicts WHERE key=?", (key,))
                    accepted += 1
        return accepted

    def retry_conflicts(self):
        with self.db() as db:
            rows = db.execute("SELECT record FROM conflicts").fetchall()
        for row in rows:
            self.receive([json.loads(row[0])])

    def _apply(self, kind, key, value, author):
        identifier = key.split(':', 1)[1]
        if kind in {'molecule_group', 'molecule_entry'}:
            from mn_molecular_design_studio.molecular_database_peer import apply_record
            apply_record(kind, key, value, author)
            return
        if kind == 'job':
            if value is None:
                with Store().connect() as db:
                    previous = db.execute("SELECT node_id,payload FROM jobs WHERE id=?", (identifier,)).fetchone()
                    if previous:
                        payload = json.loads(previous['payload'])
                        if previous['node_id'] != author and payload.get('submission_origin') != author:
                            raise PermissionError("Only the job's owning workstation may remove its superseded identity.")
                        db.execute("DELETE FROM jobs WHERE id=?", (identifier,))
                return
            if not isinstance(value, dict) or value.get('cluster_job_id') != identifier:
                raise ValueError("Invalid job identity.")
            accepted = bool(value.get('local_reference'))
            if (value.get('node_id') if accepted else value.get('submission_origin')) != author:
                raise PermissionError("Only the job's owning workstation may publish its state.")
            local_ref = value.get('local_reference')
            if accepted:
                parts = str(local_ref).split('/')
                if len(parts) != 2 or any(p in {'', '.', '..'} for p in parts):
                    raise ValueError("Invalid machine-local job reference.")
                fallback_id = config.job_id(author, '/'.join(parts))
                value = dict(value, run_dir=config.artifact_uri(author, *parts))
            # Older runs may have no recorded account/project; preserve the
            # original store's empty ownership fields and admin-only visibility.
            owner_id = str(value.get('owner_id') or '')
            project_id = str(value.get('project_id') or '')
            with Store().connect() as db:
                previous = db.execute("SELECT * FROM jobs WHERE id=?", (identifier,)).fetchone()
                if previous and (previous['owner_id'] != owner_id or previous['project_id'] != project_id or accepted and previous['node_id'] not in {None, author}):
                    raise PermissionError("Synchronization cannot change a job's execution machine or user/project ownership.")
                if previous and previous['local_ref'] and not accepted:
                    return
                if accepted:
                    superseded = db.execute("SELECT * FROM jobs WHERE node_id=? AND local_ref=? AND id<>?", (author, '/'.join(parts), identifier)).fetchone()
                    if superseded:
                        if (superseded['id'] != fallback_id or superseded['owner_id'] != owner_id
                                or superseded['project_id'] != project_id):
                            raise PermissionError("Synchronization found a conflicting identity for this workstation run.")
                        db.execute("DELETE FROM jobs WHERE id=?", (superseded['id'],))
                db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET node_id=excluded.node_id,local_ref=excluded.local_ref,payload=excluded.payload", (identifier, value.get('node_id'), local_ref, owner_id, project_id, encoded(value), value.get('created_at','')))
            return
        from mn_molecular_design_studio import auth, workspace
        if kind == 'project':
            if value is not None and (not isinstance(value, dict) or value.get('id') != identifier or not value.get('name')):
                raise ValueError("Invalid project metadata.")
            org = workspace.load_organization()
            projects = [p for p in org.get('projects',[]) if p.get('id') != identifier]
            if value is not None:
                projects.append(value)
            workspace.save_organization(dict(org, projects=projects))
            return
        with auth._connection() as db:
            if kind == 'user':
                if value is None:
                    db.execute('DELETE FROM users WHERE id=?', (identifier,))
                    return
                if not isinstance(value, dict) or set(value) != set(USER_FIELDS) or value['id'] != identifier or value['role'] not in {'admin','user'} or not str(value['password_hash']).startswith('scrypt$'):
                    raise ValueError("Invalid synchronized account.")
                auth._normalize_username(value['username'])
                previous = db.execute('SELECT password_hash,disabled FROM users WHERE id=?', (identifier,)).fetchone()
                if previous and (previous['password_hash'] != value['password_hash'] or value['disabled']):
                    db.execute('DELETE FROM sessions WHERE user_id=?', (identifier,))
                columns = ','.join(USER_FIELDS)
                updates = ','.join(k+'=excluded.'+k for k in USER_FIELDS if k != 'id')
                db.execute('INSERT INTO users ('+columns+') VALUES ('+','.join('?' for _ in USER_FIELDS)+') ON CONFLICT(id) DO UPDATE SET '+updates, tuple(value[k] for k in USER_FIELDS))
            else:
                parts = identifier.split(':')
                if len(parts) != 2:
                    raise ValueError("Invalid project membership.")
                if value is None:
                    db.execute('DELETE FROM project_memberships WHERE project_id=? AND user_id=?', parts)
                elif value.get('project_id') != parts[0] or value.get('user_id') != parts[1]:
                    raise ValueError("Membership identity cannot change.")
                else:
                    db.execute('INSERT OR REPLACE INTO project_memberships VALUES (?,?,?,?)', tuple(value[k] for k in ('project_id','user_id','granted_by','granted_at')))


def capture_identity():
    if not config.is_peer() or _importing.get():
        return
    from mn_molecular_design_studio import auth, workspace
    with locked():
        db = auth._connect()
        try:
            users = {"user:"+r['id']: dict(r) for r in db.execute('SELECT * FROM users')}
            members = {"membership:"+r['project_id']+':'+r['user_id']: dict(r) for r in db.execute('SELECT * FROM project_memberships')}
        finally:
            db.close()
        projects = {'project:'+p['id']: p for p in workspace.load_organization().get('projects',[])}
        replica = Replica()
        for kind, values in (('user',users), ('project',projects), ('membership',members)):
            replica.capture(kind, values)


def sync_peer(peer, snapshot):
    from . import transport
    replica = Replica()
    from mn_molecular_design_studio.molecular_database_peer import capture
    capture(replica)
    state = replica.state('peer:'+peer['id'], {'push':0, 'pull':0})
    response = transport.node_request('/peer/hello', {'node_id': config.identity(), 'snapshot': snapshot, 'reference_probes': Store().reference_probes()}, method='POST', base_url=peer['url'], timeout=4)
    if response.get('node_id') != peer['id']:
        raise ValueError("The peer API belongs to a different workstation. Check its SSH tunnel.")
    Store().heartbeat(peer['id'], response['snapshot'], [])
    if state.get('epoch') != response['epoch']:
        state = {'push':0, 'pull':0, 'epoch':response['epoch']}
    # Independent durable cursors let a machine catch up after any length of outage.
    for _ in range(20):
        page = replica.export(state['push'])
        if page['records']:
            transport.node_request('/peer/records', {'records':page['records']}, method='POST', base_url=peer['url'], timeout=15)
        state['push'] = page['cursor']
        replica.set_state('peer:'+peer['id'],state)
        if not page['more']:
            break
    for _ in range(20):
        page = transport.node_request('/peer/records?since='+str(state['pull']), base_url=peer['url'], timeout=15)
        replica.receive(page['records'])
        state['pull'] = page['cursor']
        replica.set_state('peer:'+peer['id'],state)
        if not page['more']:
            break
    replica.retry_conflicts()
    try:
        from mn_molecular_design_studio.molecular_database_peer import materialize_portable_structures
        materialize_portable_structures(peer=peer)
    except Exception:
        logging.getLogger(__name__).exception("Could not materialize portable molecular Database structures from %s", peer.get("id"))
    return response.get('reference_probes',[])
