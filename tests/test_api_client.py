import asyncio

import httpx
import pytest

from mn_molecular_design_studio.api_client import StudioAPIError, StudioClient


@pytest.mark.parametrize("identifier", ["00000000000000000000000000000001", "production-student"])
@pytest.mark.parametrize("http_client", [False, True])
def test_production_access_uses_explicit_mode_not_account_id(identifier, http_client):
    calls = []

    def respond(request):
        calls.append(request.url.path)
        assert request.headers["Authorization"] == "Bearer user-session"
        return httpx.Response(200, json={"id": identifier, "development_mode": False})

    client = StudioClient("http://studio.test", token="user-session", transport=httpx.MockTransport(respond))
    result = asyncio.run(client.request("GET", "/api/v1/auth/me", bearer="user-session", http_client=http_client))
    assert result["id"] == identifier
    assert calls == ["/api/v1/auth/me", "/api/v1/auth/me"]


@pytest.mark.parametrize("identifier", ["00000000000000000000000000000001", "another-development-account"])
@pytest.mark.parametrize("allow_development,http_client", [(False, False), (False, True), (True, True)])
def test_development_access_remains_rejected(identifier, allow_development, http_client):
    calls = []

    def respond(request):
        calls.append(request.url.path)
        return httpx.Response(200, json={"id": identifier, "development_mode": True})

    client = StudioClient("http://studio.test", token="user-session", allow_development=allow_development,
                          transport=httpx.MockTransport(respond))
    with pytest.raises(StudioAPIError, match="development access"):
        asyncio.run(client.request("GET", "/api/v1/projects", bearer="user-session", http_client=http_client))
    assert calls == ["/api/v1/auth/me"]


def test_local_development_requires_explicit_opt_in():
    def respond(request):
        return httpx.Response(200, json={"development_mode": True})

    client = StudioClient("http://studio.test", allow_development=True, transport=httpx.MockTransport(respond))
    assert asyncio.run(client.request("GET", "/api/v1/auth/me"))["development_mode"] is True


@pytest.mark.parametrize("actor", [{}, {"development_mode": None}, {"development_mode": 0},
    {"development_mode": 1}, {"development_mode": "false"}, [], None])
@pytest.mark.parametrize("allow_development", [False, True])
def test_missing_or_invalid_authentication_mode_fails_closed(actor, allow_development):
    calls = []

    def respond(request):
        calls.append(request.url.path)
        return httpx.Response(200, json=actor)

    client = StudioClient("http://studio.test", token="user-session", allow_development=allow_development,
                          transport=httpx.MockTransport(respond))
    with pytest.raises(StudioAPIError, match="authentication mode"):
        asyncio.run(client.request("GET", "/api/v1/projects"))
    assert calls == ["/api/v1/auth/me"]
