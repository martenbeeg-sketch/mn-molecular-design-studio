"""Request-level retry protection, including direct local workflow submission."""
from __future__ import annotations

import hashlib
import json
import sqlite3
from contextlib import contextmanager
from pathlib import Path


class Conflict(ValueError):
    pass


def submission_status(user, kind, key):
    from mn_molecular_design_studio.runtime import app_home
    path = Path(app_home()) / "api-requests.sqlite3"
    if not path.is_file():
        raise FileNotFoundError("Submission not found.")
    db = sqlite3.connect(path)
    try:
        row = db.execute("SELECT status,response FROM requests WHERE owner=? AND kind=? AND key=?",
                         (str(user["id"]), kind, key)).fetchone()
    finally:
        db.close()
    if row is None:
        raise FileNotFoundError("Submission not found.")
    return {"kind": kind, "key": key, "status": row[0], "response": json.loads(row[1]) if row[1] else None}


def submit_once(user, kind, key, payload, operation):
    if key is None:
        return operation()
    if not isinstance(key, str) or not 1 <= len(key.strip()) <= 200:
        raise ValueError("Idempotency keys must contain 1–200 characters.")
    from mn_molecular_design_studio.runtime import app_home
    root = Path(app_home())
    root.mkdir(parents=True, exist_ok=True)
    path = root / "api-requests.sqlite3"
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()

    @contextmanager
    def connect():
        db = sqlite3.connect(path, timeout=30)
        try:
            db.row_factory = sqlite3.Row
            db.execute("CREATE TABLE IF NOT EXISTS requests(owner TEXT,kind TEXT,key TEXT,fingerprint TEXT,status TEXT,response TEXT,PRIMARY KEY(owner,kind,key))")
            path.chmod(0o600)
            with db:
                yield db
        finally:
            db.close()

    identity = (str(user["id"]), kind, key)
    with connect() as db:
        db.execute("BEGIN IMMEDIATE")
        previous = db.execute("SELECT * FROM requests WHERE owner=? AND kind=? AND key=?", identity).fetchone()
        if previous:
            if previous["fingerprint"] != fingerprint:
                raise Conflict("This idempotency key was already used for a different request.")
            if previous["status"] == "completed":
                return json.loads(previous["response"])
            raise Conflict("This submission is pending or its outcome is uncertain. Inspect your jobs before resubmitting; reuse this key when polling the submission.")
        db.execute("INSERT INTO requests VALUES (?,?,?,?,?,NULL)", (*identity, fingerprint, "pending"))
    try:
        response = operation()
    except BaseException:
        # A failure may happen after a scheduler accepted work. Never blindly
        # repeat an uncertain submission after a lost reply or process crash.
        with connect() as db:
            db.execute("UPDATE requests SET status='uncertain' WHERE owner=? AND kind=? AND key=?", identity)
        raise
    with connect() as db:
        db.execute("UPDATE requests SET status='completed',response=? WHERE owner=? AND kind=? AND key=?",
                   (json.dumps(response, allow_nan=False), *identity))
    return response
