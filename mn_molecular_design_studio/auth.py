"""Local account, session, and project-access management for the Studio.

Passwords are scrypt hashes and API/session tokens are random values whose
SHA-256 digests are stored in the workspace database.
"""

from __future__ import annotations

import base64
from contextlib import contextmanager
import hashlib
import hmac
import re
import secrets
import sqlite3
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.runtime import app_home

_SCRYPT_N = 1 << 15
_SCRYPT_R = 8
_SCRYPT_P = 1
_PASSWORD_MIN_LENGTH = 12
_SESSION_HOURS = 12
_LOGIN_FAILURE_LIMIT = 5
_LOGIN_LOCK_MINUTES = 15
_USERNAME_PATTERN = re.compile(r"^[^\s@]+@[^\s@]+\.[^\s@]+$|^[a-zA-Z0-9_.-]{3,64}$")
_DEVELOPMENT_USER_ID = "00000000000000000000000000000001"


def is_development_mode() -> bool:
    """Return whether this process explicitly opted into local dev access."""
    import os

    return os.environ.get("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE", "production").strip().casefold() in {
        "development",
        "dev",
    }


def development_user() -> dict[str, Any]:
    """Synthetic local admin used only when the process explicitly runs in dev mode."""
    now = _stamp()
    return {
        "id": _DEVELOPMENT_USER_ID,
        "username": "developer",
        "display_name": "Developer",
        "role": "admin",
        "must_change_password": False,
        "disabled": False,
        "created_at": now,
        "updated_at": now,
        "development_mode": True,
    }


def _database_path() -> Path:
    path = app_home() / "database" / "user_access.sqlite3"
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        path.parent.chmod(0o700)
    except OSError:
        pass
    return path


def user_input_root(user_id: str) -> Path:
    """Return the private folder for files uploaded by one Studio account."""
    safe_id = re.sub(r"[^a-fA-F0-9_-]", "", str(user_id or ""))
    if not safe_id or safe_id != str(user_id):
        raise ValueError("Invalid user ID.")
    path = app_home() / "users" / safe_id / "inputs"
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        path.parent.chmod(0o700)
        path.chmod(0o700)
    except OSError:
        pass
    return path.resolve()


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _stamp(value: datetime | None = None) -> str:
    return (value or _now()).isoformat(timespec="seconds")


def _connect() -> sqlite3.Connection:
    path = _database_path()
    connection = sqlite3.connect(path, timeout=20)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA busy_timeout=20000")
    connection.execute("PRAGMA foreign_keys=ON")
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS users (
            id TEXT PRIMARY KEY,
            username TEXT NOT NULL UNIQUE COLLATE NOCASE,
            display_name TEXT NOT NULL,
            role TEXT NOT NULL CHECK(role IN ('admin', 'user')),
            password_hash TEXT NOT NULL,
            must_change_password INTEGER NOT NULL DEFAULT 0,
            disabled INTEGER NOT NULL DEFAULT 0,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS sessions (
            token_hash TEXT PRIMARY KEY,
            user_id TEXT NOT NULL REFERENCES users(id) ON DELETE CASCADE,
            created_at TEXT NOT NULL,
            expires_at TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS sessions_by_user ON sessions(user_id);
        CREATE TABLE IF NOT EXISTS project_memberships (
            project_id TEXT NOT NULL,
            user_id TEXT NOT NULL REFERENCES users(id) ON DELETE CASCADE,
            granted_by TEXT NOT NULL DEFAULT '',
            granted_at TEXT NOT NULL,
            PRIMARY KEY(project_id, user_id)
        );
        CREATE INDEX IF NOT EXISTS project_memberships_by_user ON project_memberships(user_id);
        CREATE TABLE IF NOT EXISTS login_attempts (
            username_key TEXT PRIMARY KEY,
            failures INTEGER NOT NULL DEFAULT 0,
            locked_until TEXT NOT NULL DEFAULT ''
        );
        """
    )
    try:
        path.chmod(0o600)
    except OSError:
        pass
    return connection


@contextmanager
def _connection():
    connection = _connect()
    try:
        with connection:
            yield connection
    finally:
        connection.close()


def _normalize_username(username: str) -> str:
    normalized = " ".join(str(username or "").split()).casefold()
    if not normalized or len(normalized) > 254 or not _USERNAME_PATTERN.fullmatch(normalized):
        raise ValueError("Use a valid email address or a username of 3–64 letters, numbers, dots, underscores, or hyphens.")
    return normalized


def _password_hash(password: str) -> str:
    if not isinstance(password, str) or len(password) < _PASSWORD_MIN_LENGTH:
        raise ValueError(f"Passwords must contain at least {_PASSWORD_MIN_LENGTH} characters.")
    if len(password) > 1024:
        raise ValueError("Password is too long.")
    salt = secrets.token_bytes(16)
    derived = hashlib.scrypt(password.encode("utf-8"), salt=salt, n=_SCRYPT_N, r=_SCRYPT_R, p=_SCRYPT_P, dklen=32)
    salt_text = base64.urlsafe_b64encode(salt).decode("ascii")
    digest_text = base64.urlsafe_b64encode(derived).decode("ascii")
    return f"scrypt${_SCRYPT_N}${_SCRYPT_R}${_SCRYPT_P}${salt_text}${digest_text}"


def _verify_password(password: str, encoded: str) -> bool:
    try:
        scheme, n_text, r_text, p_text, salt_text, digest_text = encoded.split("$")
        if scheme != "scrypt":
            return False
        salt = base64.urlsafe_b64decode(salt_text.encode("ascii"))
        expected = base64.urlsafe_b64decode(digest_text.encode("ascii"))
        actual = hashlib.scrypt(
            str(password).encode("utf-8"),
            salt=salt,
            n=int(n_text),
            r=int(r_text),
            p=int(p_text),
            dklen=len(expected),
        )
    except (ValueError, TypeError, UnicodeError):
        return False
    return hmac.compare_digest(actual, expected)


def _public_user(row: sqlite3.Row | dict[str, Any]) -> dict[str, Any]:
    item = dict(row)
    item.pop("password_hash", None)
    item["must_change_password"] = bool(item.get("must_change_password"))
    item["disabled"] = bool(item.get("disabled"))
    return item


def _require_admin_actor(actor: dict[str, Any]) -> None:
    if actor.get("role") != "admin":
        raise PermissionError("Only an admin can manage accounts or project grants.")


def user_count() -> int:
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        return 1  # Accounts are managed by the coordinator.
    with _connection() as connection:
        return int(connection.execute("SELECT count(*) FROM users").fetchone()[0])


def create_initial_admin(username: str, display_name: str, password: str) -> dict[str, Any]:
    """Create the first and only bootstrap admin; intended for the local CLI."""
    normalized = _normalize_username(username)
    clean_name = " ".join(str(display_name or "").split()) or normalized
    if len(clean_name) > 120:
        raise ValueError("Display name must be 120 characters or fewer.")
    encoded_password = _password_hash(password)
    now = _stamp()
    user = {
        "id": uuid.uuid4().hex,
        "username": normalized,
        "display_name": clean_name,
        "role": "admin",
        "password_hash": encoded_password,
        "must_change_password": 0,
        "disabled": 0,
        "created_at": now,
        "updated_at": now,
    }
    with _connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        if connection.execute("SELECT 1 FROM users LIMIT 1").fetchone():
            raise ValueError("An account already exists. Only the existing admin can add users.")
        connection.execute(
            "INSERT INTO users (id, username, display_name, role, password_hash, must_change_password, disabled, created_at, updated_at) "
            "VALUES (:id, :username, :display_name, :role, :password_hash, :must_change_password, :disabled, :created_at, :updated_at)",
            user,
        )
    return {key: value for key, value in user.items() if key != "password_hash"}


def reset_admin_password(username: str, password: str) -> None:
    """Recover the admin account from the machine that hosts the local Studio."""
    normalized = _normalize_username(username)
    encoded_password = _password_hash(password)
    with _connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        row = connection.execute("SELECT id FROM users WHERE username=? AND role='admin'", (normalized,)).fetchone()
        if row is None:
            raise ValueError("No admin account matches that username.")
        connection.execute(
            "UPDATE users SET password_hash=?, must_change_password=0, disabled=0, updated_at=? WHERE id=?",
            (encoded_password, _stamp(), row["id"]),
        )
        connection.execute("DELETE FROM sessions WHERE user_id=?", (row["id"],))


def create_user(
    username: str,
    display_name: str,
    temporary_password: str,
    *,
    actor: dict[str, Any],
) -> dict[str, Any]:
    _require_admin_actor(actor)
    normalized = _normalize_username(username)
    clean_name = " ".join(str(display_name or "").split()) or normalized
    if len(clean_name) > 120:
        raise ValueError("Display name must be 120 characters or fewer.")
    encoded_password = _password_hash(temporary_password)
    now = _stamp()
    user = {
        "id": uuid.uuid4().hex,
        "username": normalized,
        "display_name": clean_name,
        "role": "user",
        "password_hash": encoded_password,
        "must_change_password": 1,
        "disabled": 0,
        "created_at": now,
        "updated_at": now,
    }
    try:
        with _connection() as connection:
            connection.execute(
                "INSERT INTO users (id, username, display_name, role, password_hash, must_change_password, disabled, created_at, updated_at) "
                "VALUES (:id, :username, :display_name, :role, :password_hash, :must_change_password, :disabled, :created_at, :updated_at)",
                user,
            )
    except sqlite3.IntegrityError as exc:
        raise ValueError("That username already exists.") from exc
    return {key: value for key, value in user.items() if key != "password_hash"}


def list_users() -> list[dict[str, Any]]:
    with _connection() as connection:
        rows = connection.execute(
            "SELECT id, username, display_name, role, must_change_password, disabled, created_at, updated_at "
            "FROM users ORDER BY role DESC, username COLLATE NOCASE"
        ).fetchall()
    return [_public_user(row) for row in rows]


def _record_failed_login(username_key: str) -> None:
    now = _now()
    with _connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        row = connection.execute(
            "SELECT failures, locked_until FROM login_attempts WHERE username_key=?", (username_key,)
        ).fetchone()
        if row and row["locked_until"]:
            try:
                if datetime.fromisoformat(row["locked_until"]) > now:
                    return
            except ValueError:
                pass
        failures = (int(row["failures"]) if row else 0) + 1
        locked_until = _stamp(now + timedelta(minutes=_LOGIN_LOCK_MINUTES)) if failures >= _LOGIN_FAILURE_LIMIT else ""
        connection.execute(
            "INSERT INTO login_attempts(username_key, failures, locked_until) VALUES (?, ?, ?) "
            "ON CONFLICT(username_key) DO UPDATE SET failures=excluded.failures, locked_until=excluded.locked_until",
            (username_key, failures, locked_until),
        )


def _new_session(user_id: str) -> tuple[str, str]:
    token = secrets.token_urlsafe(40)
    now = _now()
    expires = now + timedelta(hours=_SESSION_HOURS)
    digest = hashlib.sha256(token.encode("utf-8")).hexdigest()
    with _connection() as connection:
        connection.execute(
            "INSERT INTO sessions(token_hash, user_id, created_at, expires_at) VALUES (?, ?, ?, ?)",
            (digest, user_id, _stamp(now), _stamp(expires)),
        )
    return token, _stamp(expires)


def authenticate(username: str, password: str) -> tuple[dict[str, Any], str, str] | None:
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        from mn_molecular_design_studio.cluster.identity_auth import authenticate as cluster_authenticate
        return cluster_authenticate(username, password)
    try:
        normalized = _normalize_username(username)
    except ValueError:
        normalized = str(username or "").strip().casefold()[:254]
    with _connection() as connection:
        attempt = connection.execute(
            "SELECT failures, locked_until FROM login_attempts WHERE username_key=?", (normalized,)
        ).fetchone()
        user = connection.execute("SELECT * FROM users WHERE username=?", (normalized,)).fetchone()
    encoded = user["password_hash"] if user else _password_hash("dummy-password-that-is-never-valid")
    valid = _verify_password(password, encoded)
    locked = False
    if attempt and attempt["locked_until"]:
        try:
            locked = datetime.fromisoformat(attempt["locked_until"]) > _now()
        except ValueError:
            locked = False
    if not user or not valid or user["disabled"] or locked:
        if not locked:
            _record_failed_login(normalized)
        return None
    with _connection() as connection:
        connection.execute("DELETE FROM login_attempts WHERE username_key=?", (normalized,))
    token, expires_at = _new_session(str(user["id"]))
    return _public_user(user), token, expires_at


def resolve_session(token: object) -> dict[str, Any] | None:
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        from mn_molecular_design_studio.cluster.identity_auth import resolve
        return resolve(token)
    if not isinstance(token, str) or not token:
        return None
    digest = hashlib.sha256(token.encode("utf-8")).hexdigest()
    with _connection() as connection:
        row = connection.execute(
            "SELECT u.*, s.expires_at AS session_expires_at FROM sessions s JOIN users u ON u.id=s.user_id "
            "WHERE s.token_hash=? AND s.expires_at>? AND u.disabled=0",
            (digest, _stamp()),
        ).fetchone()
    return dict(_public_user(row), session_expires_at=row["session_expires_at"]) if row else None


def revoke_session(token: object) -> None:
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        from mn_molecular_design_studio.cluster.identity_auth import revoke
        return revoke(token)
    if not isinstance(token, str) or not token:
        return
    digest = hashlib.sha256(token.encode("utf-8")).hexdigest()
    with _connection() as connection:
        connection.execute("DELETE FROM sessions WHERE token_hash=?", (digest,))


def change_password(user_id: str, new_password: str, *, current_password: str | None = None, forced: bool = False) -> None:
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        from mn_molecular_design_studio.cluster import identity_auth, transport
        transport.request("/api/v1/auth/password", {"new_password": new_password, "current_password": current_password},
                          token=identity_auth.current_token(), method="POST")
        return
    encoded_password = _password_hash(new_password)
    with _connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        row = connection.execute("SELECT password_hash, must_change_password FROM users WHERE id=? AND disabled=0", (user_id,)).fetchone()
        if row is None:
            raise ValueError("This account is unavailable.")
        if not forced and (not current_password or not _verify_password(current_password, row["password_hash"])):
            raise ValueError("Current password is incorrect.")
        if forced and not row["must_change_password"]:
            raise ValueError("A password change is no longer required. Sign in again and use account settings.")
        connection.execute(
            "UPDATE users SET password_hash=?, must_change_password=0, updated_at=? WHERE id=?",
            (encoded_password, _stamp(), user_id),
        )
        connection.execute("DELETE FROM sessions WHERE user_id=?", (user_id,))


def reset_user_password(user_id: str, temporary_password: str, *, actor: dict[str, Any]) -> None:
    _require_admin_actor(actor)
    encoded_password = _password_hash(temporary_password)
    with _connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        user = connection.execute("SELECT role FROM users WHERE id=?", (user_id,)).fetchone()
        if user is None:
            raise ValueError("User account not found.")
        connection.execute(
            "UPDATE users SET password_hash=?, must_change_password=1, updated_at=? WHERE id=?",
            (encoded_password, _stamp(), user_id),
        )
        connection.execute("DELETE FROM sessions WHERE user_id=?", (user_id,))


def set_user_disabled(user_id: str, disabled: bool, *, actor: dict[str, Any]) -> None:
    _require_admin_actor(actor)
    actor_id = str(actor.get("id") or "")
    with _connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        row = connection.execute("SELECT role, disabled FROM users WHERE id=?", (user_id,)).fetchone()
        if row is None:
            raise ValueError("User account not found.")
        if user_id == actor_id:
            raise ValueError("You cannot disable your own admin account.")
        if disabled and row["role"] == "admin":
            active_admins = int(connection.execute("SELECT count(*) FROM users WHERE role='admin' AND disabled=0").fetchone()[0])
            if active_admins <= 1:
                raise ValueError("The last active admin account cannot be disabled.")
        connection.execute("UPDATE users SET disabled=?, updated_at=? WHERE id=?", (int(disabled), _stamp(), user_id))
        if disabled:
            connection.execute("DELETE FROM sessions WHERE user_id=?", (user_id,))


def accessible_project_ids(user: dict[str, Any]) -> list[str]:
    if "_cluster_project_ids" in user:
        return list(user["_cluster_project_ids"])
    if user.get("role") == "admin":
        from mn_molecular_design_studio.workspace import load_organization

        organization = load_organization()
        return sorted(
            str(project["id"])
            for project in organization.get("projects", [])
            if isinstance(project, dict) and project.get("id")
        )
    with _connection() as connection:
        rows = connection.execute(
            "SELECT project_id FROM project_memberships WHERE user_id=? ORDER BY project_id", (user["id"],)
        ).fetchall()
    return [str(row["project_id"]) for row in rows]


def project_member_ids(project_id: str) -> list[str]:
    with _connection() as connection:
        rows = connection.execute(
            "SELECT user_id FROM project_memberships WHERE project_id=? ORDER BY user_id", (project_id,)
        ).fetchall()
    return [str(row["user_id"]) for row in rows]


def grant_project_access(
    project_id: str,
    user_id: str,
    *,
    actor: dict[str, Any],
    allow_self_grant: bool = False,
) -> None:
    if actor.get("role") != "admin" and not (allow_self_grant and str(actor.get("id")) == str(user_id)):
        raise PermissionError("Only an admin can grant project access.")
    from mn_molecular_design_studio.workspace import load_organization

    organization = load_organization()
    project_ids = {
        str(project.get("id")) for project in organization.get("projects", []) if isinstance(project, dict)
    }
    if project_id not in project_ids:
        raise ValueError("Project not found.")
    with _connection() as connection:
        connection.execute(
            "INSERT INTO project_memberships(project_id, user_id, granted_by, granted_at) VALUES (?, ?, ?, ?) "
            "ON CONFLICT(project_id, user_id) DO UPDATE SET granted_by=excluded.granted_by, granted_at=excluded.granted_at",
            (project_id, user_id, str(actor.get("id") or ""), _stamp()),
        )


def revoke_project_access(project_id: str, user_id: str, *, actor: dict[str, Any]) -> None:
    _require_admin_actor(actor)
    with _connection() as connection:
        connection.execute("DELETE FROM project_memberships WHERE project_id=? AND user_id=?", (project_id, user_id))


def require_project_access(user: dict[str, Any], project_id: str | None) -> None:
    if not project_id or user.get("role") == "admin":
        return
    if project_id not in set(accessible_project_ids(user)):
        raise PermissionError("You do not have access to this project.")


def can_access_job(user: dict[str, Any], job: dict[str, Any] | None) -> bool:
    if not job:
        return user.get("role") == "admin"
    if user.get("role") == "admin":
        return True
    if str(job.get("owner_id") or "") == str(user.get("id") or ""):
        return True
    project_id = str(job.get("project_id") or "")
    return bool(project_id and project_id in set(accessible_project_ids(user)))


def require_job_access(user: dict[str, Any], job: dict[str, Any] | None) -> None:
    if not can_access_job(user, job):
        raise PermissionError("You do not have access to this job.")


def project_is_member(project_id: str, user_id: str) -> bool:
    with _connection() as connection:
        row = connection.execute(
            "SELECT 1 FROM project_memberships WHERE project_id=? AND user_id=?", (project_id, user_id)
        ).fetchone()
    return row is not None


def password_min_length() -> int:
    return _PASSWORD_MIN_LENGTH
