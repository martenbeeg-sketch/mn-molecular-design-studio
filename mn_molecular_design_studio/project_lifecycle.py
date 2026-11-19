"""Confirmed project cleanup with durable workstation receipts and recovery files."""

from datetime import datetime, timezone
from functools import wraps
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import uuid

from mn_molecular_design_studio.cluster import config
from mn_molecular_design_studio.cluster.replication import identity_edit, locked
from mn_molecular_design_studio.projects import TERMINAL_STATUSES


PROTOCOL = "project-lifecycle.v1"
BUSY = {"preparing", "running", "attention"}


class ProjectUnavailableError(ValueError):
    pass


def stamp():
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def project_writable(project):
    return not (project.get("archived_at") or project.get("deleted_at")
                or (project.get("lifecycle") or {}).get("status") in BUSY)


def serialized(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        with locked():
            return function(*args, **kwargs)
    return wrapped


def exclusive_operation(function):
    @wraps(function)
    def wrapped(self, user, operation_id):
        path = self.root / "locks" / (identifier(operation_id) + ".lock")
        path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        with path.open("a") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                return function(self, user, operation_id)
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)
    return wrapped


def require_mutable_project(project_id):
    from mn_molecular_design_studio.workspace import _studio_home, load_organization
    project = next((row for row in load_organization(_studio_home()).get("projects", []) if row.get("id") == project_id), {})
    if project.get("deleted_at") or (project.get("lifecycle") or {}).get("status") in BUSY:
        raise PermissionError("Project is locked for cleanup or has been deleted.")


def identifier(value):
    return str(uuid.UUID(str(value)))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def file_digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from strings(child)


def _admin(user):
    if not user.get("id") or user.get("role") != "admin" or user.get("disabled"):
        raise PermissionError("Only an active administrator can clean up project data.")


def _read(path):
    return json.loads(path.read_text())


@identity_edit
def change_project(user, project_id, changes, *, home=None):
    from mn_molecular_design_studio.workspace import _studio_home, load_organization, save_organization
    home = Path(home or _studio_home())
    with locked():
        organization = load_organization(home)
        project = next((row for row in organization.get("projects", []) if row.get("id") == project_id), None)
        if project is None:
            raise ValueError("Project not found.")
        if user.get("role") != "admin" and project.get("created_by") != user.get("id"):
            raise PermissionError("Only the project owner or an administrator can manage this project.")
        if project.get("deleted_at"):
            raise ProjectUnavailableError("This project has been deleted. Its recovery receipt is retained.")
        project.update(changes, updated_at=stamp())
        save_organization(organization, home)
        return dict(project)


def archive_project(user, project_id, archived=True, *, home=None):
    from mn_molecular_design_studio.workspace import _studio_home, load_organization
    home = Path(home or _studio_home())
    if config.is_worker():
        raise ValueError("Manage project visibility on its coordinator or a peer workstation.")
    with locked():
        project = next((row for row in load_organization(home).get("projects", []) if row.get("id") == project_id), {})
        if (project.get("lifecycle") or {}).get("status") in BUSY:
            raise ValueError("Finish the pending cleanup before changing project visibility.")
        return change_project(user, project_id, {"archived_at": stamp() if archived else None}, home=home)


def install_guards(connection):
    connection.execute("CREATE TABLE IF NOT EXISTS project_lifecycle_locks (project_id TEXT PRIMARY KEY, operation_id TEXT NOT NULL, state TEXT NOT NULL)")
    connection.execute("CREATE TABLE IF NOT EXISTS project_lifecycle_paths (path TEXT PRIMARY KEY, operation_id TEXT NOT NULL)")
    for table in ("molecule_groups", "molecule_members", "molecule_structures", "molecule_observations", "molecule_sources", "molecule_import_batches"):
        exists = connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)).fetchone()
        if not exists:
            continue
        project = "NEW.project_id" if table == "molecule_groups" else "(SELECT project_id FROM molecule_groups WHERE id=NEW.group_id)"
        for operation in ("INSERT", "UPDATE"):
            connection.execute(f"CREATE TRIGGER IF NOT EXISTS project_lifecycle_{table}_{operation.lower()} "
                f"BEFORE {operation} ON {table} WHEN EXISTS (SELECT 1 FROM project_lifecycle_locks "
                f"WHERE project_id={project} AND state IN ('locked','deleted')) "
                "BEGIN SELECT RAISE(ABORT,'Project is locked for cleanup or has been deleted'); END")
    for operation in ("INSERT", "UPDATE"):
        connection.execute(f"CREATE TRIGGER IF NOT EXISTS project_lifecycle_retained_path_{operation.lower()} "
            f"BEFORE {operation} ON molecule_structures WHEN EXISTS (SELECT 1 FROM project_lifecycle_paths "
            "WHERE NEW.path=path OR substr(NEW.path,1,length(path)+1)=path||'/') "
            "BEGIN SELECT RAISE(ABORT,'Source run is locked for cleanup or has been removed'); END")
        connection.execute(f"CREATE TRIGGER IF NOT EXISTS project_lifecycle_review_{operation.lower()} "
            f"BEFORE {operation} ON complex_review_annotations WHEN NEW.source_kind='database' AND EXISTS ("
            "SELECT 1 FROM molecule_structures s JOIN molecule_groups g ON g.id=s.group_id "
            "JOIN project_lifecycle_locks l ON l.project_id=g.project_id "
            "WHERE s.id=NEW.structure_ref AND l.state IN ('locked','deleted')) "
            "BEGIN SELECT RAISE(ABORT,'Project review is locked for cleanup'); END")


class ProjectNode:
    def __init__(self, *, home=None, runs=None, database=None, node_id=None, job_reader=None, participant_ids=None):
        from mn_molecular_design_studio.workspace import _studio_home
        from mn_molecular_design_studio.runtime import runs_root
        from mn_molecular_design_studio.molecular_database import MolecularDatabase
        self.home = Path(home or _studio_home()).resolve()
        self.runs = Path(runs or runs_root()).resolve()
        self.database = database or MolecularDatabase()
        self.node_id = str(node_id or config.identity())
        self.root = self.home / "organization/project-lifecycle"
        self.job_reader = job_reader
        self.participant_ids = participant_ids

    def save(self, kind, key, value):
        path = self.root / kind / (identifier(key) + ".json")
        path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        config.atomic_json(path, value)
        return value

    def read(self, kind, key, user):
        _admin(user)
        value = _read(self.root / kind / (identifier(key) + ".json"))
        if value["actor_id"] != str(user["id"]):
            raise PermissionError("This cleanup belongs to another administrator.")
        return value

    def project(self, project_id, *, allow_deleted=False):
        from mn_molecular_design_studio.workspace import load_organization
        project = next((row for row in load_organization(self.home).get("projects", []) if row.get("id") == project_id), None)
        if project is None or (project.get("deleted_at") and not allow_deleted):
            raise ProjectUnavailableError("Project not found or already deleted.")
        return dict(project)

    def _jobs(self):
        if self.job_reader is not None:
            return self.job_reader()
        from mn_molecular_design_studio.job_catalog import reconcile, search_jobs
        reconcile(self.runs)
        rows, offset = [], 0
        while True:
            page = search_jobs(limit=1000, offset=offset, include_hidden=True, access_admin=True)
            rows.extend(page["items"])
            offset += len(page["items"])
            if offset >= page["total"]:
                return rows
            if not page["items"]:
                raise ValueError("The job catalog is incomplete; cleanup is blocked.")

    def run_path(self, reference):
        parts = str(reference).split("/")
        if len(parts) != 2 or any(part in {"", ".", ".."} or "\\" in part for part in parts):
            raise ValueError("Invalid run reference.")
        path = self.runs.joinpath(*parts)
        if path.is_symlink() or path.parent.is_symlink() or not path.resolve().is_relative_to(self.runs):
            raise ValueError("A run path escapes managed Results storage.")
        return path

    def tree(self, path):
        rows = []
        if not path.exists():
            return rows
        for directory, directories, files in os.walk(path, followlinks=False):
            for name in sorted(directories + files):
                current = Path(directory) / name
                relative = current.relative_to(path).as_posix()
                if current.is_symlink():
                    rows.append([relative, "link", os.readlink(current)])
                elif current.is_file():
                    info = current.stat()
                    rows.append([relative, "file", info.st_size, info.st_mtime_ns])
                elif current.is_dir():
                    rows.append([relative, "directory"])
                else:
                    raise ValueError("A run contains an unsupported special file.")
        return sorted(rows)

    def inspect(self, user, project_id, *, mode, selected_runs=None, include_collections=False, persist=True, participants=None):
        _admin(user)
        if mode not in {"clean", "delete"} or type(include_collections) is not bool:
            raise ValueError("Choose cleanup or project deletion.")
        project = self.project(project_id)
        expected_participants = self.participant_ids or [self.node_id, *[peer["id"] for peer in config.load().get("peers", [])]]
        if participants is not None and sorted(participants) != sorted(expected_participants):
            raise ValueError("The trusted workstation roster differs; cleanup cannot establish complete coverage.")
        jobs = self._jobs()
        project_jobs = [row for row in jobs if row.get("project_id") == project_id]
        requested = set(selected_runs or [])
        selected = [row for row in project_jobs if mode == "delete" or self.node_id + ":" + row["task_group"] + "/" + row["run_id"] in requested]
        references = {row["task_group"] + "/" + row["run_id"] for row in selected}
        blocked, runs, dependencies = [], [], []
        known_references = {row["task_group"] + "/" + row["run_id"] for row in jobs}
        if mode == "delete":
            for metadata_path in self.runs.glob('*/*/metadata.json'):
                reference = metadata_path.parent.relative_to(self.runs).as_posix()
                path = self.run_path(reference)
                metadata = _read(metadata_path)
                marker = _read(path / 'molecular_design_studio.json') if (path / 'molecular_design_studio.json').exists() else {}
                if project_id in {metadata.get('project_id'), marker.get('project_id')} and reference not in known_references:
                    blocked.append(reference + ': project run is not reconciled in the authoritative catalog')
        catalog_jobs = []
        if self.job_reader is None and config.enabled():
            from mn_molecular_design_studio.cluster.store import Store
            offset = 0
            while True:
                page = Store().search(user, project_id=project_id, include_hidden=True, limit=1000, offset=offset)
                for row in page["items"]:
                    if row.get("node_id") not in (None, "", self.node_id):
                        continue
                    reference = str(row.get("local_reference") or "")
                    if mode != "delete" and self.node_id + ":" + reference not in requested:
                        continue
                    if str(row.get("status") or "").casefold() not in TERMINAL_STATUSES:
                        blocked.append("Catalog job " + str(row.get("cluster_job_id") or reference) + ": active, waiting or unresolved")
                    else:
                        catalog_jobs.append(str(row["cluster_job_id"]))
                offset += len(page["items"])
                if offset >= page["total"]:
                    break
                if not page["items"]:
                    raise ValueError("The workstation catalog is incomplete.")
        for row in selected:
            reference = row["task_group"] + "/" + row["run_id"]
            path = self.run_path(reference)
            if str(row.get("status") or "").casefold() not in TERMINAL_STATUSES:
                blocked.append(reference + ": job is active or unresolved")
            if path.exists():
                metadata = _read(path / "metadata.json")
                marker = _read(path / "molecular_design_studio.json") if (path / "molecular_design_studio.json").exists() else {}
                actual_project = metadata.get("project_id") or marker.get("project_id")
                if actual_project != project_id or str(metadata.get("status") or "").casefold() not in TERMINAL_STATUSES:
                    blocked.append(reference + ": live metadata no longer matches the reviewed finished project job")
                if metadata.get('project_id') and marker.get('project_id') and metadata['project_id'] != marker['project_id']:
                    blocked.append(reference + ': metadata and ownership marker disagree about the project')
            runs.append({"reference": reference, "cluster_job_id": row.get("cluster_job_id"), "tree": self.tree(path)})
        for row in jobs:
            reference = row["task_group"] + "/" + row["run_id"]
            upstream = str(row.get("upstream_task_group") or "") + "/" + str(row.get("upstream_run_id") or "")
            if upstream != "/":
                dependencies.append({"source": str(row.get("upstream_node_id") or self.node_id) + ":" + upstream,
                                     "consumer": self.node_id + ":" + reference})
            directory = self.run_path(reference)
            for name in ("input.json", "pipeline.json"):
                if not (directory / name).is_file():
                    continue
                for text in strings(_read(directory / name)):
                    if text.startswith("cluster:"):
                        dependencies.append({"source": text, "consumer": self.node_id + ":" + reference})
                        continue
                    if text.startswith("runs:///"):
                        candidate = self.runs / text.removeprefix("runs:///")
                    elif text.startswith("run:///"):
                        candidate = directory / text.removeprefix("run:///")
                    elif text.startswith("/"):
                        candidate = Path(text)
                    elif text.split("/", 1)[0] in {job["task_group"] for job in jobs}:
                        candidate = self.runs / text
                    else:
                        continue
                    try:
                        parts = candidate.resolve().relative_to(self.runs).parts
                    except (OSError, ValueError):
                        continue
                    if len(parts) >= 2 and "/".join(parts[:2]) != reference:
                        dependencies.append({"source": self.node_id + ":" + "/".join(parts[:2]), "consumer": self.node_id + ":" + reference})
        groups, retained, structures, fingerprint = [], [], [], hashlib.sha256()
        with self.database.connection() as connection:
            local_groups = connection.execute("SELECT * FROM molecule_groups WHERE project_id=? ORDER BY id", (project_id,)).fetchall()
            for raw in local_groups:
                group = dict(raw)
                origin = json.loads(group["origin"])
                if origin.get("_replica_node"):
                    if mode == "delete" and origin["_replica_node"] not in expected_participants:
                        blocked.append("A project collection belongs to an uncovered workstation: " + str(origin["_replica_node"]))
                    continue
                for table in ("molecule_groups", "molecule_members", "molecule_structures", "molecule_observations", "molecule_sources", "molecule_import_batches"):
                    query = "SELECT * FROM " + table + (" WHERE id=?" if table == "molecule_groups" else " WHERE group_id=?") + " ORDER BY 1,2"
                    for record in connection.execute(query, (group["id"],)):
                        fingerprint.update(json.dumps(dict(record), sort_keys=True).encode())
                for record in connection.execute("SELECT a.* FROM complex_review_annotations a JOIN molecule_structures s ON s.id=a.structure_ref WHERE a.source_kind='database' AND s.group_id=? ORDER BY a.id", (group["id"],)):
                    fingerprint.update(json.dumps(dict(record), sort_keys=True).encode())
                members = [row[0] for row in connection.execute("SELECT molecule_id FROM molecule_members WHERE group_id=? ORDER BY molecule_id", (group["id"],))]
                automatic = origin.get("type") == "job"
                recorded_runs = set()
                has_saved = False
                unqualified = False
                for table, column in (("molecule_sources", "provenance"), ("molecule_structures", "context"), ("molecule_observations", "context")):
                    for row in connection.execute(f"SELECT {column} FROM {table} WHERE group_id=?", (group["id"],)):
                        context = json.loads(row[0])
                        has_saved = has_saved or context.get("saved_selection") in (True, 1)
                        if context.get("task_group") and context.get("run_id"):
                            recorded_runs.add(str(context["task_group"]) + "/" + str(context["run_id"]))
                for row in connection.execute("SELECT run_id,context FROM molecule_observations WHERE group_id=?", (group["id"],)):
                    if not json.loads(row["context"]).get("task_group") and row["run_id"] != origin.get("run_id"):
                        unqualified = True
                source_ref = str(origin.get("task_group") or "") + "/" + str(origin.get("run_id") or "")
                relevant = recorded_runs or ({source_ref} if source_ref != "/" else set())
                remove = (mode == "delete" and (automatic or include_collections)) or (
                    mode == "clean" and ((automatic and relevant and relevant.issubset(references)) or (not automatic and include_collections)))
                if automatic and mode == "clean" and relevant & references and not relevant.issubset(references):
                    blocked.append(group["name"] + ": selected and unselected runs share one SQL collection; select the complete recorded cohort or keep it.")
                if automatic and mode == "clean" and relevant & references and unqualified:
                    blocked.append(group["name"] + ": legacy score lineage is incomplete; reconcile it or review the whole project instead.")
                if mode == "delete" and has_saved and not include_collections:
                    remove = False
                if remove:
                    if automatic and mode == "clean":
                        saved = {row[0] for row in connection.execute("SELECT molecule_id FROM molecule_sources WHERE group_id=? AND json_extract(provenance,'$.saved_selection')=1", (group["id"],))}
                        members = [value for value in members if value not in saved]
                    groups.append({"id": group["id"], "name": group["name"], "members": members, "automatic": automatic})
                    for row in connection.execute("SELECT * FROM molecule_structures WHERE group_id=? AND molecule_id IN (SELECT value FROM json_each(?))", (group["id"], json.dumps(members))):
                        structures.append(dict(row))
                elif members:
                    retained.append({"id": group["id"], "name": group["name"], "members": len(members)})
            removing = {(group["id"], member) for group in groups for member in group["members"]}
            for row in connection.execute("SELECT group_id,molecule_id,path FROM molecule_structures"):
                if (row["group_id"], row["molecule_id"]) in removing:
                    continue
                path = Path(row["path"])
                if path.is_absolute():
                    for run in runs:
                        if path.resolve().is_relative_to(self.run_path(run["reference"])):
                            blocked.append(run["reference"] + ": a retained or shared structure still uses its files")
        if mode == "delete" and retained:
            blocked.append("Saved/imported collections remain: explicitly include them or archive the project instead.")
        value = {"id": str(uuid.uuid4()), "protocol": PROTOCOL, "node_id": self.node_id,
                 "actor_id": str(user["id"]), "project_id": project_id, "project_name": project["name"],
                 "mode": mode, "selected_runs": sorted(requested), "include_collections": include_collections,
                 "catalog_jobs": sorted(catalog_jobs),
                 "runs": runs, "groups": groups, "retained": retained, "structures": structures,
                 "dependencies": dependencies, "blocked": sorted(set(blocked)),
                 "project_snapshot": {key: value for key, value in project.items() if key not in {"lifecycle", "updated_at"}},
                 "participants": sorted(expected_participants),
                 "database_revision": fingerprint.hexdigest(), "generated_at": stamp()}
        value["fingerprint"] = digest({key: value[key] for key in ("project_id", "project_name", "project_snapshot", "mode", "runs", "catalog_jobs", "groups", "retained", "structures", "database_revision")})
        return self.save("plans", value["id"], value) if persist else value

    @serialized
    def begin(self, user, plan_id, operation_id):
        plan = self.read("plans", plan_id, user)
        operation_id = identifier(operation_id)
        path = self.root / "nodes" / (operation_id + ".json")
        existing = None
        if path.exists():
            existing = self.read("nodes", operation_id, user)
            if existing["plan_id"] != plan_id:
                raise ValueError("Operation identity was already used for another plan.")
            if existing["status"] == "aborted":
                raise ValueError("This cleanup was aborted; prepare a fresh review.")
            if existing["status"] != "preparing":
                return existing
        project = self.project(plan["project_id"])
        current = project.get("lifecycle") or {}
        if current.get("status") in BUSY and current.get("operation_id") != operation_id:
            raise ValueError("Another cleanup is already in progress.")
        fresh = self.inspect(user, plan["project_id"], mode=plan["mode"], selected_runs=plan["selected_runs"], include_collections=plan["include_collections"], persist=False, participants=plan["participants"])
        if fresh["fingerprint"] != plan["fingerprint"] or fresh["blocked"]:
            raise ValueError("Project data changed or has blocking references. Prepare a fresh review.")
        operation = existing or {"id": operation_id, "plan_id": plan_id, "actor_id": str(user["id"]),
                     "project_id": plan["project_id"], "status": "preparing", "moved": [],
                     "removed_groups": [], "backup": None, "previous_lifecycle": current}
        self.save("nodes", operation_id, operation)
        change_project(user, plan["project_id"], {"lifecycle": {"operation_id": operation_id, "status": "preparing"}}, home=self.home)
        with self.database.connection() as connection:
            install_guards(connection)
            connection.execute("INSERT INTO project_lifecycle_locks VALUES (?,?,?) ON CONFLICT(project_id) DO UPDATE SET operation_id=excluded.operation_id,state=excluded.state", (plan["project_id"], operation_id, "locked"))
            for run in plan["runs"]:
                connection.execute("INSERT OR REPLACE INTO project_lifecycle_paths VALUES (?,?)", (str(self.run_path(run["reference"])), operation_id))
        fresh = self.inspect(user, plan["project_id"], mode=plan["mode"], selected_runs=plan["selected_runs"], include_collections=plan["include_collections"], persist=False, participants=plan["participants"])
        if fresh["fingerprint"] != plan["fingerprint"] or fresh["blocked"]:
            self.abort(user, operation_id)
            raise ValueError("Project changed while the cleanup lock was being established. Review it again.")
        operation["status"] = "prepared"
        return self.save("nodes", operation_id, operation)

    @serialized
    def abort(self, user, operation_id):
        path = self.root / "nodes" / (identifier(operation_id) + ".json")
        if not path.exists():
            _admin(user)
            return {"status": "not_started"}
        operation = self.read("nodes", operation_id, user)
        if operation["moved"] or operation.get("backup") or operation["status"] in {"data_removed", "completed"}:
            raise ValueError("Cleanup has started; continue it instead of clearing its lock.")
        project = self.project(operation["project_id"])
        if (project.get("lifecycle") or {}).get("operation_id") == operation_id:
            change_project(user, operation["project_id"], {"lifecycle": operation["previous_lifecycle"]}, home=self.home)
        with self.database.connection() as connection:
            connection.execute("DELETE FROM project_lifecycle_locks WHERE project_id=? AND operation_id=?", (operation["project_id"], operation_id))
            connection.execute("DELETE FROM project_lifecycle_paths WHERE operation_id=?", (operation_id,))
        operation["status"] = "aborted"
        return self.save("nodes", operation_id, operation)

    def _backup(self, operation, plan):
        recovery_root = self.runs.parent / '.project-recovery'
        recovery_root.mkdir(parents=True, exist_ok=True, mode=0o700)
        destination = recovery_root / operation["id"]
        if destination.is_symlink() or destination.parent.is_symlink():
            raise ValueError("Recovery storage cannot be a symbolic link.")
        destination.mkdir(parents=True, exist_ok=True, mode=0o700)
        snapshot = destination / "molecular_library.sqlite3"
        if not operation.get("backup"):
            linked_bytes = sum((self.run_path(run['reference']) / row[0]).stat().st_size for run in plan['runs'] for row in run['tree']
                               if row[1] == 'file' and (self.run_path(run['reference']) / row[0]).is_file()
                               and (self.run_path(run['reference']) / row[0]).stat().st_nlink > 1)
            required = self.database.path.stat().st_size + linked_bytes + sum(Path(row["path"]).stat().st_size for row in plan["structures"] if Path(row["path"]).is_file())
            if shutil.disk_usage(destination).free < required + 16 * 1024 * 1024:
                raise OSError("Not enough free space for a verified recovery snapshot.")
            temporary = destination / "snapshot.tmp"
            with sqlite3.connect(self.database.path) as source, sqlite3.connect(temporary) as target:
                source.backup(target)
                if target.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
                    raise ValueError("Recovery database verification failed.")
            temporary.replace(snapshot)
            snapshot.chmod(0o600)
            files = []
            for row in plan["structures"]:
                source_path = Path(row["path"])
                if not source_path.is_absolute() or not source_path.is_file():
                    raise ValueError("A retained structure is unavailable for recovery: " + row["id"])
                content_digest = file_digest(source_path)
                context = json.loads(row["context"])
                if any(context.get(key) for key in ("sha256", "portable_sha256", "content_sha256")):
                    self.database._verified_structure_content(row, source_path)
                copied = destination / "structures" / content_digest
                copied.parent.mkdir(exist_ok=True, mode=0o700)
                if not copied.exists():
                    shutil.copyfile(source_path, copied)
                    copied.chmod(0o600)
                if file_digest(copied) != content_digest:
                    raise ValueError("Recovery structure checksum mismatch.")
                files.append({"id": row["id"], "original_path": row["path"], "sha256": content_digest, "backup": str(copied)})
            receipt = {"database": str(snapshot), "database_sha256": file_digest(snapshot), "structures": files,
                       "project_id": plan["project_id"], "project": plan["project_snapshot"],
                       "operation_id": operation["id"], "verified_at": stamp(), "runs": {}}
            receipt["auxiliary_databases"] = []
            cluster_path = (config.home() if self.job_reader is None else self.home / 'cluster') / 'coordinator.sqlite3'
            for name, source_path in (("job_catalog", self.home / "database/job_catalog.sqlite3"), ("cluster", cluster_path)):
                if source_path.is_file():
                    copied = destination / (name + '.sqlite3')
                    with sqlite3.connect(source_path) as source, sqlite3.connect(copied) as target:
                        source.backup(target)
                        if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                            raise ValueError('Auxiliary recovery database verification failed.')
                    copied.chmod(0o600)
                    receipt["auxiliary_databases"].append({"path": str(copied), "sha256": file_digest(copied)})
            config.atomic_json(destination / "receipt.json", receipt)
            operation["backup"] = receipt
            self.save("nodes", operation["id"], operation)
        elif file_digest(snapshot) != operation["backup"]["database_sha256"]:
            raise ValueError("Recovery database checksum changed; cleanup is blocked.")
        for row in operation["backup"]["structures"]:
            if file_digest(Path(row["backup"])) != row["sha256"]:
                raise ValueError("Recovery structure checksum changed; cleanup is blocked.")
        return destination

    @exclusive_operation
    def execute(self, user, operation_id):
        operation = self.read("nodes", operation_id, user)
        if operation["status"] in {"data_removed", "completed"}:
            return operation
        plan = self.read("plans", operation["plan_id"], user)
        project = self.project(plan["project_id"])
        if (project.get("lifecycle") or {}).get("operation_id") != operation_id:
            raise ValueError("Cleanup no longer owns the project lock.")
        with self.database.connection() as connection:
            lock = connection.execute("SELECT operation_id,state FROM project_lifecycle_locks WHERE project_id=?", (plan["project_id"],)).fetchone()
        if not lock or lock["operation_id"] != operation_id or lock["state"] != "locked":
            raise ValueError("The SQL cleanup lock is missing or changed; no data can be removed.")
        if operation["status"] not in {"prepared", "running"}:
            raise ValueError("The workstation cleanup has not been prepared.")
        operation["status"] = "running"
        self.save("nodes", operation_id, operation)
        destination = self._backup(operation, plan)
        for run in plan["runs"]:
            reference = run["reference"]
            source = self.run_path(reference)
            target = destination / "runs" / reference
            if reference not in operation["moved"]:
                if source.exists() and target.exists():
                    raise ValueError("Both live and recovery run directories exist; manual reconciliation is required.")
                if source.exists():
                    if self.tree(source) != run["tree"]:
                        raise ValueError("Run files changed after review; cleanup stopped.")
                    target.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
                    source.rename(target)
                if run["tree"] and (not target.exists() or self.tree(target) != run["tree"]):
                    raise ValueError("Recovery run verification failed.")
                for row in run['tree']:
                    saved = target / row[0]
                    if row[1] == 'file' and saved.stat().st_nlink > 1:
                        temporary = saved.with_name(saved.name + '.recovery-' + uuid.uuid4().hex)
                        shutil.copy2(saved, temporary)
                        temporary.replace(saved)
                operation["backup"]["runs"][reference] = [
                    {"path": row[0], "sha256": file_digest(target / row[0])} for row in run["tree"] if row[1] == "file"]
                operation["moved"].append(reference)
                self.save("nodes", operation_id, operation)
            for row in operation["backup"]["runs"].get(reference, []):
                if file_digest(target / row["path"]) != row["sha256"]:
                    raise ValueError("Recovery run checksum changed; cleanup is blocked.")
        config.atomic_json(destination / "receipt.json", operation["backup"])
        with self.database.connection() as connection:
            connection.execute("BEGIN IMMEDIATE")
            install_guards(connection)
            connection.execute("UPDATE project_lifecycle_locks SET state='mutating' WHERE project_id=? AND operation_id=?", (plan["project_id"], operation_id))
            structure_ids = json.dumps([row["id"] for row in plan["structures"]])
            connection.execute("DELETE FROM complex_review_annotations WHERE source_kind='database' AND structure_ref IN (SELECT value FROM json_each(?))", (structure_ids,))
            for group in plan["groups"]:
                if not group["members"]:
                    continue
                members = json.dumps(group["members"])
                remaining = connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id NOT IN (SELECT value FROM json_each(?))", (group["id"], members)).fetchone()[0]
                if not remaining:
                    connection.execute("DELETE FROM molecule_import_batches WHERE group_id=?", (group["id"],))
                for table in ("molecule_sources", "molecule_structures", "molecule_observations", "molecule_members"):
                    connection.execute(f"DELETE FROM {table} WHERE group_id=? AND molecule_id IN (SELECT value FROM json_each(?))", (group["id"], members))
                current = connection.execute("SELECT origin FROM molecule_groups WHERE id=?", (group["id"],)).fetchone()
                if current:
                    origin = json.loads(current[0])
                    origin["removed_molecule_ids"] = sorted(set(origin.get("removed_molecule_ids") or []) | set(group["members"]))
                    connection.execute("UPDATE molecule_groups SET origin=?,updated_at=? WHERE id=?", (json.dumps(origin), stamp(), group["id"]))
            connection.execute("UPDATE project_lifecycle_locks SET state='locked' WHERE project_id=? AND operation_id=?", (plan["project_id"], operation_id))
        for group in plan["groups"]:
            if config.is_peer() and group["members"]:
                from mn_molecular_design_studio.molecular_database_peer import publish_member_removal
                publish_member_removal(self.database, user, group["id"], group["members"])
        if self.job_reader is None:
            from mn_molecular_design_studio.job_catalog import reconcile
            reconcile(self.runs)
        from mn_molecular_design_studio.cluster.store import Store
        if config.enabled() and self.job_reader is None:
            store = Store()
            job_ids = set(plan["catalog_jobs"]) | {run["cluster_job_id"] for run in plan["runs"] if run.get("cluster_job_id")}
            for job_id in job_ids:
                try:
                    store.purge_metadata(job_id, user)
                except LookupError:
                    pass
        operation.update(status="data_removed", completed_at=stamp(), removed_groups=[group["id"] for group in plan["groups"]])
        return self.save("nodes", operation_id, operation)

    @serialized
    def finish(self, user, operation_id):
        operation = self.read("nodes", operation_id, user)
        if operation["status"] == "completed":
            return operation
        if operation["status"] == "aborted":
            raise ValueError("This cleanup was aborted. Prepare a fresh review.")
        if operation["status"] != "data_removed":
            raise ValueError("Workstation data cleanup has not completed.")
        plan = self.read("plans", operation["plan_id"], user)
        project = self.project(plan["project_id"], allow_deleted=True)
        if (project.get("lifecycle") or {}).get("operation_id") != operation_id:
            raise ValueError("The project cleanup lock changed.")
        changes = {"lifecycle": {"operation_id": operation_id, "status": "completed", "completed_at": stamp()}}
        if plan["mode"] == "delete":
            changes["deleted_at"] = stamp()
        with self.database.connection() as connection:
            connection.execute("BEGIN IMMEDIATE")
            connection.execute("UPDATE project_lifecycle_locks SET state='mutating' WHERE project_id=? AND operation_id=?", (plan["project_id"], operation_id))
            for raw in connection.execute("SELECT id,origin FROM molecule_groups WHERE project_id=?", (plan["project_id"],)).fetchall():
                origin = json.loads(raw["origin"])
                reference = str(origin.get("_replica_node") or "") + ":" + str(origin.get("task_group") or "") + "/" + str(origin.get("run_id") or "")
                remove_replica = origin.get("_replica_node") and (plan["mode"] == "delete" or (
                    origin.get("type") == "job" and reference in plan["selected_runs"]) or (plan["include_collections"] and origin.get("type") != "job"))
                if remove_replica:
                    for table in ("molecule_sources", "molecule_structures", "molecule_observations", "molecule_members"):
                        connection.execute(f"DELETE FROM {table} WHERE group_id=?", (raw["id"],))
                if plan["mode"] == "delete":
                    origin["project_deleted"] = operation_id
                    connection.execute("UPDATE molecule_groups SET origin=? WHERE id=?", (json.dumps(origin), raw["id"]))
            if plan["mode"] == "delete":
                connection.execute("UPDATE project_lifecycle_locks SET state='deleted' WHERE project_id=? AND operation_id=?", (plan["project_id"], operation_id))
            else:
                connection.execute("DELETE FROM project_lifecycle_locks WHERE project_id=? AND operation_id=?", (plan["project_id"], operation_id))
        if not project.get("deleted_at"):
            change_project(user, plan["project_id"], changes, home=self.home)
        operation.update(status="completed")
        return self.save("nodes", operation_id, operation)


class ProjectCoordinator:
    def __init__(self, node=None, peers=None, remote=None):
        self.node = node or ProjectNode()
        self.peers = list(config.load().get("peers", []) if peers is None else peers)
        self.remote = remote or self._remote

    def _remote(self, peer, method, payload):
        from mn_molecular_design_studio.cluster.transport import node_request
        result = node_request("/projects/lifecycle/" + method, payload, method="POST", base_url=peer["url"], timeout=600)
        if result.get("node_id") != peer["id"] or result.get("protocol") != PROTOCOL:
            raise ValueError("Workstation does not confirm the project cleanup protocol.")
        return result["result"]

    def call(self, peer, method, payload, user):
        if peer is None:
            return getattr(self.node, method)(user, **payload)
        return self.remote(peer, method, {**payload, "actor_id": str(user["id"])})

    def pending(self, user, project_id=None):
        _admin(user)
        rows = [_read(path) for path in (self.node.root / "operations").glob("*.json")]
        return sorted([row for row in rows if row["actor_id"] == str(user["id"]) and (project_id is None or row["project_id"] == project_id)
                       and row["status"] not in {"review", "completed", "aborted"}], key=lambda row: row["created_at"])

    def prepare(self, user, project_id, *, mode, selected_runs=None, include_collections=False):
        _admin(user)
        if (self.node.project(project_id).get("lifecycle") or {}).get("status") in BUSY:
            raise ValueError("Continue the existing cleanup before preparing a new operation.")
        if config.load()["mode"] in {"coordinator", "worker"}:
            raise ValueError("Project-wide cleanup requires standalone or upgraded peer mode.")
        plans = {}
        for peer in [None, *self.peers]:
            plan = self.call(peer, "inspect", {"project_id": project_id, "mode": mode, "selected_runs": selected_runs,
                                             "include_collections": include_collections,
                                             "participants": [self.node.node_id, *[peer["id"] for peer in self.peers]]}, user)
            plans[plan["node_id"]] = plan
        all_selected = {node_id + ":" + run["reference"] for node_id, plan in plans.items() for run in plan["runs"]}
        cluster_references = {"cluster:" + str(run["cluster_job_id"]): node_id + ":" + run["reference"]
                              for node_id, plan in plans.items() for run in plan["runs"] if run.get("cluster_job_id")}
        blocked = [problem for plan in plans.values() for problem in plan["blocked"]]
        if mode == "clean":
            blocked.extend("Selected run is missing from the authoritative workstation inventory: " + value for value in sorted(set(selected_runs or []) - all_selected))
            if not all_selected and not any(group["members"] for plan in plans.values() for group in plan["groups"]):
                raise ValueError("No removable data matches this selection.")
        for plan in plans.values():
            for dependency in plan["dependencies"]:
                if cluster_references.get(dependency["source"], dependency["source"]) in all_selected and dependency["consumer"] not in all_selected:
                    blocked.append("An unselected downstream job still depends on " + dependency["source"])
        value = {"id": str(uuid.uuid4()), "actor_id": str(user["id"]), "project_id": project_id,
                 "project_name": self.node.project(project_id)["name"], "mode": mode, "plans": plans,
                 "blocked": sorted(set(blocked)), "status": "review", "nodes": {}, "created_at": stamp(),
                 "peers": [{"id": peer["id"], "url": peer["url"]} for peer in self.peers]}
        return self.node.save("operations", value["id"], value)

    def resume(self, user, operation_id, *, confirmation, backup_acknowledged):
        operation = self.node.read("operations", operation_id, user)
        if confirmation != operation["project_name"] or backup_acknowledged is not True:
            raise ValueError("Type the exact project name and acknowledge the verified recovery backup.")
        if operation["blocked"]:
            raise ValueError("Resolve every blocking reference before cleanup.")
        if operation["status"] == "completed":
            return operation
        if operation["status"] == "aborted":
            raise ValueError("This cleanup was aborted. Prepare a fresh review.")
        if sorted(peer["id"] for peer in self.peers) != sorted(peer["id"] for peer in operation["peers"]):
            raise ValueError("Trusted workstations changed. Reconcile the cleanup plan before continuing.")
        endpoints = [None, *self.peers]
        if not operation.get("confirmed_at"):
            operation["confirmed_at"] = stamp()
            self.node.save("operations", operation_id, operation)
        try:
            for method in ("begin", "execute", "finish"):
                ordered = [*self.peers, None] if method == "finish" else endpoints
                for peer in ordered:
                    node_id = peer["id"] if peer else self.node.node_id
                    stage = operation["nodes"].setdefault(node_id, {})
                    if method in stage:
                        continue
                    arguments = {"operation_id": operation_id}
                    if method == "begin":
                        arguments["plan_id"] = operation["plans"][node_id]["id"]
                    stage[method] = self.call(peer, method, arguments, user)
                    operation.update(status="running", error=None)
                    self.node.save("operations", operation_id, operation)
            operation.update(status="completed", completed_at=stamp())
        except (OSError, ValueError, sqlite3.Error, RuntimeError) as exc:
            operation.update(status="attention", error=str(exc))
            if method == "begin" and isinstance(exc, ValueError):
                released = True
                for peer in endpoints:
                    try:
                        self.call(peer, "abort", {"operation_id": operation_id}, user)
                    except (OSError, ValueError, sqlite3.Error, RuntimeError):
                        released = False
                if released:
                    operation["status"] = "aborted"
        return self.node.save("operations", operation_id, operation)


def runtime_coordinator():
    return ProjectCoordinator()
