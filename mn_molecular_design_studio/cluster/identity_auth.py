"""Delegate worker sign-in to the coordinator; cache verified sessions locally."""
from __future__ import annotations

import hashlib
import json
import time
from contextvars import ContextVar
from datetime import datetime

from . import config, transport

_token = ContextVar("cluster_user_token", default="")


def current_token():
    return _token.get()


def resolve(token):
    if not isinstance(token, str) or not token:
        return None
    digest = hashlib.sha256(token.encode()).hexdigest()
    directory = config.home() / "verified-sessions"
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    path = directory / (digest + ".json")
    try:
        user = transport.request("/api/v1/auth/me", token=token, timeout=3)
    except ValueError:
        path.unlink(missing_ok=True)
        return None
    except transport.ClusterUnavailable:
        try:
            cached = json.loads(path.read_text())
            if cached["expires_at"] <= time.time():
                return None
            user = dict(cached["user"], coordinator_offline=True)
        except (FileNotFoundError, KeyError, json.JSONDecodeError):
            return None
    else:
        expiry = user.get("session_expires_at")
        expires_at = datetime.fromisoformat(expiry).timestamp() if expiry else time.time() + 3600
        config.atomic_json(path, {"user": user, "expires_at": expires_at})
    _token.set(token)
    user["_cluster_project_ids"] = user.get("project_ids", [])
    return user


def authenticate(username, password):
    try:
        result = transport.request("/api/v1/auth/login", {"username": username, "password": password}, method="POST")
    except ValueError:
        return None
    user = resolve(result["access_token"])
    return (user, result["access_token"], result.get("expires_at")) if user else None


def revoke(token):
    if not isinstance(token, str) or not token:
        return
    path = config.home() / "verified-sessions" / (hashlib.sha256(token.encode()).hexdigest() + ".json")
    path.unlink(missing_ok=True)
    try:
        transport.request("/api/v1/auth/logout", token=token, method="POST")
    except (ValueError, transport.ClusterUnavailable):
        pass
