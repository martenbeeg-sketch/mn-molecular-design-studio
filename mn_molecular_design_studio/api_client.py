"""Small authenticated REST client, with no imports of scientific/UI services."""
from __future__ import annotations

import json
from urllib.parse import unquote, urlsplit

import httpx


class StudioAPIError(RuntimeError):
    pass


def public_path(path):
    decoded = unquote(str(path))
    if (not decoded.startswith("/api/v1/") or any(c in decoded for c in ("?", "#", "\\", "\0", "%"))
            or any(part in {".", ".."} for part in decoded.split("/"))
            or decoded == "/api/v1/cluster/node" or decoded.startswith("/api/v1/cluster/node/")):
        raise ValueError("Choose a public /api/v1/ operation; node endpoints and external URLs are not client tools.")
    return decoded


class StudioClient:
    def __init__(self, base_url, *, token=None, allow_development=False, transport=None, max_response_bytes=1024 * 1024):
        parsed = urlsplit(base_url)
        if (parsed.scheme not in {"http", "https"} or not parsed.hostname or parsed.username or parsed.password
                or parsed.query or parsed.fragment or parsed.path not in {"", "/"}):
            raise ValueError("Use the Studio API origin, without a path, credentials or query.")
        self.base_url = base_url.rstrip("/")
        self.token = token
        self.allow_development = allow_development
        self.transport = transport
        self.max_response_bytes = max_response_bytes

    async def request(self, method, path, *, body=None, query=None, headers=None,
                      bearer=None, http_client=False, content=False, content_bytes=None):
        path = public_path(path)
        method = str(method).upper()
        if method not in {"GET", "POST", "PUT", "PATCH", "DELETE"}:
            raise ValueError("Unsupported HTTP method.")
        # Incoming HTTP credentials never fall back to a process-wide token.
        token = bearer if http_client else bearer or self.token
        if not token and (http_client or not self.allow_development):
            raise StudioAPIError("A Studio user session token is required.")
        request_headers = dict(headers or {})
        if token:
            request_headers["Authorization"] = "Bearer " + token
        async with httpx.AsyncClient(base_url=self.base_url, transport=self.transport,
                                     timeout=60, follow_redirects=False) as client:
            identity = await client.get("/api/v1/auth/me", headers=request_headers)
            if identity.status_code != 200:
                raise StudioAPIError(f"Studio authentication failed ({identity.status_code}).")
            actor = identity.json()
            development = actor.get("development_mode") or actor.get("id") == "00000000000000000000000000000001"
            if development and (http_client or not self.allow_development):
                raise StudioAPIError("Studio uses development access. Use production accounts for HTTP clients; local stdio may explicitly opt into development access.")
            async with client.stream(method, path, params=query, json=body if content_bytes is None else None,
                                     content=content_bytes, headers=request_headers) as response:
                media_type = response.headers.get("content-type", "application/octet-stream")
                is_json = "application/json" in media_type or "+json" in media_type
                if response.is_success and not is_json and not content:
                    return {"download_url": self.base_url + path, "path": path,
                            "requires_authentication": True, "media_type": media_type,
                            "size_bytes": response.headers.get("content-length"),
                            "filename": response.headers.get("content-disposition"),
                            "resource_hint": "Read the corresponding studio://files/ or studio://structures/ MCP resource for small files."}
                chunks, size = [], 0
                async for chunk in response.aiter_bytes():
                    size += len(chunk)
                    if size > self.max_response_bytes:
                        if content:
                            raise StudioAPIError("This file is too large for an inline MCP resource. Download it through the authenticated API.")
                        return {"status_code": response.status_code, "response_too_large": True,
                                "path": path, "message": "Use a paginated query or authenticated file download. Do not repeat a submission with a new idempotency key."}
                    chunks.append(chunk)
                raw = b"".join(chunks)
                if not response.is_success:
                    try:
                        detail = json.loads(raw).get("detail", "Request rejected.")
                    except (ValueError, AttributeError):
                        detail = "Request rejected."
                    raise StudioAPIError(f"Studio API {response.status_code}: {detail}")
                if content:
                    return raw
                return json.loads(raw) if raw else {"status_code": response.status_code}
