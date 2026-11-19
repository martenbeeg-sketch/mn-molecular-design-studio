"""Keep local stdio credentials private and bound to the Studio API origin."""
from __future__ import annotations

import json
import os
import stat
import tempfile
from pathlib import Path
from urllib.parse import urlsplit

import httpx

from mn_molecular_design_studio.api_client import StudioAPIError, StudioClient

DEFAULT_API_URL = "http://127.0.0.1:8766"


def _session_file():
    override = os.environ.get("STUDIO_MCP_SESSION_FILE")
    if override:
        return Path(override).expanduser()
    state = Path(os.environ.get("XDG_STATE_HOME") or Path.home() / ".local/state")
    return state / "mn-molecular-design-studio" / "mcp-session.json"


def saved_token(origin):
    try:
        fd = os.open(_session_file(), os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    except FileNotFoundError:
        return None
    except OSError as exc:
        raise StudioAPIError("Could not read the local Studio MCP session. Run mcp-login again.") from exc
    with os.fdopen(fd, encoding="utf-8") as handle:
        info = os.fstat(handle.fileno())
        if not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid() or info.st_mode & 0o077 or info.st_size > 4096:
            raise StudioAPIError("The Studio MCP session file must be private and owned by this account. Run mcp-login again.")
        try:
            session = json.load(handle)
        except (ValueError, UnicodeError) as exc:
            raise StudioAPIError("The local Studio MCP session is invalid. Run mcp-login again.") from exc
    token = session.get("token") if isinstance(session, dict) and session.get("origin") == origin else None
    return token if isinstance(token, str) and token and token != "development-mode" else None


def sign_in(origin, username, password):
    return _connect(origin, "login", {"username": username, "password": password})


def pair_with_studio(origin, code):
    return _connect(origin, "mcp/redeem", {"code": code.strip()})


def _connect(origin, endpoint, payload):
    pairing = endpoint == "mcp/redeem"
    origin = StudioClient(origin).base_url
    url = urlsplit(origin)
    if url.scheme != "https" and url.hostname not in {"localhost", "127.0.0.1", "::1"}:
        raise StudioAPIError("Use HTTPS or a local SSH tunnel to sign in to Studio MCP.")
    with httpx.Client(base_url=origin, timeout=15, follow_redirects=False, trust_env=False) as client:
        response = client.post(f"/api/v1/auth/{endpoint}", json=payload)
        if response.status_code != 200:
            if pairing:
                if response.status_code == 404:
                    raise StudioAPIError("This Studio release does not support Codex pairing.")
                raise StudioAPIError("Could not pair Codex. Check the API address and generate a new code in Studio Account.")
            raise StudioAPIError("Studio sign-in failed. Check the API address and your Studio account.")
        record = response.json()
        token = record.get("access_token")
        if not isinstance(token, str) or not token or token == "development-mode":
            raise StudioAPIError("Studio MCP sign-in requires a production credential.")
        identity = client.get("/api/v1/auth/me", headers={"Authorization": "Bearer " + token})
        actor = identity.json() if identity.status_code == 200 else {}
        if (actor.get("development_mode") is not False or actor.get("must_change_password")
                or not actor.get("id") or actor["id"] != (record.get("user") or {}).get("id")
                or (actor.get("mcp_session") is not True if pairing else not record.get("expires_at"))):
            raise StudioAPIError("Finish signing in to the Studio app and resolve any required password change first.")
    path = _session_file()
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            json.dump({"origin": origin, "token": token}, handle)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return actor["username"]
