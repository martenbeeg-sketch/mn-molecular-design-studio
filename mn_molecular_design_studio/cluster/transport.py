from __future__ import annotations

import json
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from . import config


class ClusterUnavailable(RuntimeError):
    pass


def request(path: str, payload=None, *, token: str = "", method: str = "GET", binary: bytes | None = None,
            timeout: float = 15, base_url: str | None = None, raw: bool = False):
    address = (base_url or config.load()["coordinator_url"]).rstrip("/")
    if not address:
        raise ClusterUnavailable("No coordinator is configured.")
    headers = {"Authorization": "Bearer " + token}
    data = binary
    if binary is not None:
        headers["Content-Type"] = "application/octet-stream"
    elif payload is not None:
        data = json.dumps(payload, allow_nan=False).encode()
        headers["Content-Type"] = "application/json"
    try:
        with urlopen(Request(address + path, data=data, headers=headers, method=method), timeout=timeout) as response:
            body = response.read()
            return body if raw else json.loads(body)
    except HTTPError as exc:
        try:
            message = json.loads(exc.read()).get("detail", str(exc))
        except (ValueError, UnicodeError):
            message = str(exc)
        if exc.code in {401, 403, 404, 409, 422}:
            raise ValueError(str(message)) from exc
        raise ClusterUnavailable(f"Coordinator returned HTTP {exc.code}: {message}") from exc
    except (URLError, TimeoutError, OSError) as exc:
        raise ClusterUnavailable("Coordinator is unavailable. Existing local jobs continue; remote status is last known.") from exc


def node_request(path: str, payload=None, **kwargs):
    return request("/api/v1/cluster/node" + path, payload, token=config.load()["node_token"], **kwargs)
