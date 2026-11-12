"""Session-local metadata snapshots, invalidated by all SQLite commits and fresh ACLs."""
from __future__ import annotations

from collections import OrderedDict
import json
import sqlite3
import threading

from mn_molecular_design_studio.database_browser import build_browser_index, subject_results


class BrowserCache:
    def __init__(self, path):
        self.path = path.resolve()
        self.connection = None
        self.identity = None
        self.generation = 0
        self.lock = threading.RLock()
        self.index_key = None
        self.index = None
        self.results = OrderedDict()

    def close(self):
        with self.lock:
            if self.connection is not None:
                self.connection.close()
                self.connection = None

    def __del__(self):
        self.close()

    def _token(self):
        stat = self.path.stat()
        identity = (stat.st_dev, stat.st_ino)
        if self.connection is None or identity != self.identity:
            self.close()
            self.connection = sqlite3.connect(self.path.as_uri() + "?mode=ro", uri=True, check_same_thread=False, timeout=15)
            self.connection.execute("PRAGMA query_only=ON")
            self.identity = identity
            self.generation += 1
        return (*identity, self.generation, self.connection.execute("PRAGMA data_version").fetchone()[0])

    def get_index(self, database, user, **filters):
        with self.lock:
            if self.path.exists():
                self._token()
            access = {key: filters[key] for key in ("project_id", "scope", "hide_smoke_test") if key in filters}
            visible = database.groups(user, **{key: value for key, value in filters.items() if key != "include_related"})
            accessible = database.groups(user, collection="all", **access)
            token = self._token()
            key = (token, user["id"], user.get("role"), json.dumps(filters, sort_keys=True),
                   tuple(sorted(row["id"] for row in visible)), tuple(sorted(row["id"] for row in accessible)))
            if key != self.index_key:
                index = build_browser_index(database, user, **filters)
                self.results.clear()
                self.index = index
                self.index_key = key if token == self._token() else None
            return self.index

    def get_results(self, database, index, subject, *, run_view=False, include_inputs=False):
        with self.lock:
            if index is not self.index or not subject["group_ids"] <= index["groups"].keys():
                raise ValueError("This molecular snapshot is no longer current. Refresh the browser.")
            token = self._token()
            key = (token, self.index_key, subject["id"], run_view, include_inputs)
            if key not in self.results:
                results = subject_results(database, index, subject, run_view=run_view, include_inputs=include_inputs)
                if self.index_key is None or token != self._token():
                    return results
                self.results[key] = results
                while len(self.results) > 4:
                    self.results.popitem(last=False)
            self.results.move_to_end(key)
            return self.results[key]
