"""Session-local metadata snapshots, invalidated by scientific edits and fresh ACLs."""
from __future__ import annotations

from collections import OrderedDict
import hashlib
import json
import sqlite3
import threading

from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects, molecular_subjects, run_subjects, subject_results
from mn_molecular_design_studio.database_outcomes import scientific_outcomes


MAX_SAVED_SCOPES = 8
MAX_SAVED_ROWS = 32000
_SHARED_INDEXES = OrderedDict()
_SHARED_INDEX_LOCK = threading.RLock()


def _subject_key(subject):
    return hashlib.sha256(json.dumps(subject, sort_keys=True, default=lambda value: sorted(value)).encode()).hexdigest()


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
        self.outcomes = OrderedDict()
        self.saved = OrderedDict()
        self.subjects = {}

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
        from mn_molecular_design_studio.browser_revision import revision
        current = revision(self.connection)
        token = ("browser", current) if current is not None else ("sqlite", self.connection.execute("PRAGMA data_version").fetchone()[0])
        return (*identity, self.generation, token)

    def get_index(self, database, user, **filters):
        with self.lock:
            if self.path.exists():
                self._token()
            access = {key: filters[key] for key in ("project_id", "scope", "hide_smoke_test") if key in filters}
            visible = database.groups(user, **{key: value for key, value in filters.items()
                                               if key not in {"include_related", "hide_pipeline_batches"}})
            accessible = database.groups(user, collection="all", **access)
            token = self._token()
            key = (token, user["id"], user.get("role"), json.dumps(filters, sort_keys=True),
                   tuple(sorted(row["id"] for row in visible)), tuple(sorted(row["id"] for row in accessible)))
            if key != self.index_key:
                shared_key = (str(self.path), *token[:2], token[-1], *key[1:]) if token[-1][0] == "browser" else None
                with _SHARED_INDEX_LOCK:
                    if shared_key is not None and shared_key in _SHARED_INDEXES:
                        index = dict(_SHARED_INDEXES[shared_key])
                        _SHARED_INDEXES.move_to_end(shared_key)
                    else:
                        index = build_browser_index(database, user, **filters)
                        if shared_key is not None and token == self._token():
                            _SHARED_INDEXES[shared_key] = dict(index)
                            while len(_SHARED_INDEXES) > 2:
                                _SHARED_INDEXES.popitem(last=False)
                self.results.clear()
                self.outcomes.clear()
                self.saved.clear()
                self.subjects.clear()
                self.index = index
                self.index_key = key if token == self._token() else None
            return self.index

    def get_subjects(self, index, *, role=None, individual_sequences=False, combine_campaigns=None):
        with self.lock:
            if index is not self.index:
                raise ValueError("This molecular snapshot is no longer current. Refresh the browser.")
            key = (role, individual_sequences, combine_campaigns)
            if key not in self.subjects:
                self.subjects[key] = (run_subjects(index, combine_campaigns=combine_campaigns)
                    if combine_campaigns is not None else component_subjects(index, role=role)
                    if role is not None else molecular_subjects(index, individual_sequences=individual_sequences))
            return self.subjects[key]

    def get_results(self, database, index, subject, *, run_view=False, include_inputs=False, include_metrics=False):
        with self.lock:
            if index is not self.index or not subject["group_ids"] <= index["groups"].keys():
                raise ValueError("This molecular snapshot is no longer current. Refresh the browser.")
            token = self._token()
            key = (token, self.index_key, _subject_key(subject), run_view, include_inputs, include_metrics)
            if key not in self.results:
                results = subject_results(database, index, subject, run_view=run_view, include_inputs=include_inputs,
                                          include_metrics=include_metrics)
                if self.index_key is None or token != self._token():
                    return results
                self.results[key] = results
                while len(self.results) > 4:
                    self.results.popitem(last=False)
            self.results.move_to_end(key)
            return self.results[key]

    def get_outcomes(self, database, index, subject, *, run_view=False):
        with self.lock:
            token = self._token()
            results = self.get_results(database, index, subject, run_view=run_view)
            key = (token, self.index_key, _subject_key(subject), run_view)
            if key not in self.outcomes:
                outcomes = scientific_outcomes(index, results)
                if self.index_key is None or token != self._token():
                    return outcomes
                self.outcomes[key] = outcomes
                while len(self.outcomes) > 4:
                    self.outcomes.popitem(last=False)
            self.outcomes.move_to_end(key)
            return self.outcomes[key]

    def get_saved_results(self, database, index, subject, *, run_view=False, include_references=False, fallback_references=False):
        from mn_molecular_design_studio.data_browser import saved_results

        with self.lock:
            if index is not self.index or not subject["group_ids"] <= index["groups"].keys():
                raise ValueError("This molecular snapshot is no longer current. Refresh the browser.")
            token = self._token()
            key = (token, self.index_key, _subject_key(subject), run_view, include_references, fallback_references)
            if key in self.saved:
                self.saved.move_to_end(key)
                return self.saved[key][0]
            raw = self.get_results(database, index, subject, run_view=run_view, include_inputs=include_references)
            inventory = saved_results(index, raw, include_references=include_references)
            if not inventory and fallback_references and not include_references:
                raw = self.get_results(database, index, subject, run_view=run_view, include_inputs=True)
                inventory = saved_results(index, raw, include_references=True)
            weight = max(1, len(raw))
            if self.index_key is not None and token == self.index_key[0] and token == self._token() and weight <= MAX_SAVED_ROWS:
                self.saved[key] = (inventory, weight)
                while len(self.saved) > MAX_SAVED_SCOPES or sum(entry[1] for entry in self.saved.values()) > MAX_SAVED_ROWS:
                    self.saved.popitem(last=False)
            return inventory
