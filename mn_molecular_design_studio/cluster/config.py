from __future__ import annotations

import json
import os
import socket
import uuid
from pathlib import Path
from urllib.parse import urlsplit


def home() -> Path:
    from mn_molecular_design_studio.workspace import _studio_home
    path = Path(os.environ.get("MDS_CLUSTER_HOME") or _studio_home() / "cluster").expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    return path


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            os.chmod(temporary, 0o600)
            json.dump(payload, handle, indent=2)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def load() -> dict:
    path = home() / "config.json"
    try:
        payload = json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        payload = {}
    return {"mode": "disabled", "name": socket.gethostname(), "coordinator_url": "",
            "node_token": "", "default_placement": "local", "heartbeat_seconds": 15,
            "accept_jobs": True, "peers": [], **payload}


def identity() -> str:
    path = home() / "node-id"
    try:
        return str(uuid.UUID(path.read_text().strip()))
    except FileNotFoundError:
        value = str(uuid.uuid4())
        try:
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        except FileExistsError:
            return identity()
        with os.fdopen(fd, "w") as handle:
            handle.write(value + "\n")
        return value


def save(changes: dict) -> dict:
    config = load() | changes
    if config["mode"] not in {"disabled", "coordinator", "worker", "peer"}:
        raise ValueError("Choose standalone, peer, coordinator, or worker mode.")
    if not str(config["name"]).strip():
        raise ValueError("Enter a workstation name.")
    address = str(config.get("coordinator_url") or "").rstrip("/")
    if config["mode"] != "disabled":
        parsed = urlsplit(address)
        if parsed.scheme not in {"http", "https"} or not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment or parsed.path:
            raise ValueError("Enter the coordinator HTTP(S) address, including its port.")
    if type(config["accept_jobs"]) is not bool:
        raise ValueError("Accept jobs must be enabled or disabled.")
    if config["default_placement"] not in {"local", "auto"}:
        raise ValueError("Default placement must be local or auto.")
    peers = config.get("peers", [])
    if not isinstance(peers, list) or len(peers) > 100:
        raise ValueError("Enter at most 100 trusted peer addresses.")
    seen = set()
    for peer in peers:
        identity = str(uuid.UUID(peer["id"]))
        parsed = urlsplit(str(peer.get("url", "")))
        if identity in seen or parsed.scheme not in {"http", "https"} or not parsed.hostname or parsed.username or parsed.password or parsed.path or parsed.query or parsed.fragment:
            raise ValueError("Each peer needs a unique workstation ID and an HTTP(S) API address.")
        seen.add(identity)
    config["coordinator_url"] = address
    atomic_json(home() / "config.json", config)
    return config


def enabled() -> bool:
    return load()["mode"] != "disabled"


def is_worker() -> bool:
    return load()["mode"] == "worker"


def job_id(node_id: str, local_reference: str) -> str:
    return str(uuid.uuid5(uuid.UUID(node_id), local_reference))


def artifact_uri(node_id: str, group: str, run_id: str, relative: str = "") -> str:
    return f"cluster://{node_id}/{group}/{run_id}" + ("/" + relative.lstrip("/") if relative else "")


def is_peer() -> bool:
    return load()["mode"] == "peer"
