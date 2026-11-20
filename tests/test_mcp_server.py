"""Public MCP screening import contract remains explicit and bounded."""

import asyncio
import json

import httpx
import pytest


def test_final_shortlist_mcp_tool_defaults_to_no_import_and_accepts_exact_ids():
    from mn_molecular_design_studio.mcp_server import create_server

    server = create_server(client=object())
    tools = asyncio.run(server.list_tools())
    schema = next(tool.input_schema for tool in tools if tool.name == "import_final_shortlist")
    properties = schema["properties"]
    assert properties["limit"]["default"] == 0
    assert properties["selected_compound_ids"]["default"] is None
    assert properties["selected_compound_ids"]["anyOf"][0]["items"]["type"] == "string"
    assert schema["required"] == ["task_group", "run_id"]


def test_mcp_login_saves_only_a_private_origin_bound_session(tmp_path, monkeypatch):
    from mn_molecular_design_studio import mcp_session
    from mn_molecular_design_studio.api_client import StudioAPIError

    path = tmp_path / "private" / "mcp.json"
    monkeypatch.setenv("STUDIO_MCP_SESSION_FILE", str(path))
    client_class = httpx.Client
    def request(req):
        if req.url.path == "/api/v1/auth/login":
            assert json.loads(req.content) == {"username": "fixture", "password": "private-password"}
            return httpx.Response(200, json={"access_token": "fixture-session", "expires_at": "fixture-time", "user": {"id": "u"}})
        assert req.headers["authorization"] == "Bearer fixture-session"
        return httpx.Response(200, json={"id": "u", "username": "fixture", "development_mode": False})

    monkeypatch.setattr(mcp_session.httpx, "Client", lambda **kwargs: client_class(transport=httpx.MockTransport(request), **kwargs))
    with pytest.raises(StudioAPIError, match="HTTPS"):
        mcp_session.sign_in("http://studio.test", "fixture", "private-password")
    assert mcp_session.sign_in("http://127.0.0.1:8766", "fixture", "private-password") == "fixture"
    assert path.stat().st_mode & 0o777 == 0o600
    assert "private-password" not in path.read_text()
    assert mcp_session.saved_token("http://127.0.0.1:8766") == "fixture-session"
    assert mcp_session.saved_token("http://127.0.0.1:8765") is None
    path.chmod(0o644)
    with pytest.raises(StudioAPIError, match="private"):
        mcp_session.saved_token("http://127.0.0.1:8766")


def test_mcp_pairing_uses_the_hidden_cli_prompt_and_saves_no_code(tmp_path, monkeypatch):
    from typer.testing import CliRunner
    from mn_molecular_design_studio import mcp_session
    from mn_molecular_design_studio.cli import app

    path = tmp_path / "mcp.json"
    monkeypatch.setenv("STUDIO_MCP_SESSION_FILE", str(path))
    monkeypatch.setenv("STUDIO_API_URL", "http://127.0.0.1:8766")
    client_class = httpx.Client

    def request(req):
        if req.url.path == "/api/v1/auth/mcp/redeem":
            assert json.loads(req.content) == {"code": "one-time-test-code"}
            return httpx.Response(200, json={"access_token": "paired-session", "user": {"id": "instructor-id"}})
        assert req.url.path == "/api/v1/auth/me"
        assert req.headers["authorization"] == "Bearer paired-session"
        return httpx.Response(200, json={"id": "instructor-id", "username": "research_instructor",
                                         "development_mode": False, "mcp_session": True})

    monkeypatch.setattr(mcp_session.httpx, "Client", lambda **kwargs: client_class(transport=httpx.MockTransport(request), **kwargs))
    result = CliRunner().invoke(app, ["mcp-login", "--pair"], input="one-time-test-code\n")
    assert result.exit_code == 0, result.output
    assert "research_instructor" in result.output
    assert "paired-session" not in result.output and "one-time-test-code" not in result.output
    assert path.stat().st_mode & 0o777 == 0o600
    assert "one-time-test-code" not in path.read_text()
    assert mcp_session.saved_token("http://127.0.0.1:8766") == "paired-session"


def test_results_read_tools_do_not_request_remote_inspection():
    from mcp import Client
    from mn_molecular_design_studio.mcp_server import create_server

    class Recorder:
        async def request(self, method, path, **kwargs):
            return {"method": method, "path": path}

    async def check():
        async with Client(create_server(Recorder())) as client:
            for name, expected in [("job_results", "GET"), ("request_job_inspection", "POST")]:
                result = await client.call_tool(name, {"reference": "cluster:fixture"})
                assert result.structured_content["method"] == expected
            result = await client.call_tool("job_results", {"reference": "ligand-docking/fixture"})
            assert result.structured_content["path"] == "/api/v1/executions/local%3Aligand-docking%3Afixture/results"
            result = await client.call_tool("result_presentation", {"reference": "ligand-docking/fixture"})
            assert result.structured_content == {"method": "GET", "path": "/api/v1/jobs/local%3Aligand-docking%3Afixture/presentation"}
    asyncio.run(check())


def test_result_presentation_retries_large_pages_without_skipping_rows():
    from mcp import Client
    from mn_molecular_design_studio.mcp_server import create_server

    requests = []

    class Recorder:
        async def request(self, method, path, *, query, **kwargs):
            assert method == "GET"
            assert path.endswith("/presentation")
            requests.append(dict(query))
            limit, offset = query["limit"], query["offset"]
            if limit > 25:
                return {"response_too_large": True, "status_code": 200}
            return {"items": list(range(offset, offset + limit)), "offset": offset, "limit": limit,
                    "total": 150, "next_offset": offset + limit}

    async def check():
        async with Client(create_server(Recorder())) as client:
            parameters = {"reference": "design/fixture", "offset": 13}
            result = await client.call_tool("result_presentation", parameters)
            assert result.structured_content["items"] == list(range(13, 38))
            assert result.structured_content["next_offset"] == 38
            assert requests == [{"offset": 13, "limit": size} for size in (100, 50, 25)]

    asyncio.run(check())


def test_stdio_tools_explain_missing_login(monkeypatch, tmp_path):
    from mcp import Client
    from mn_molecular_design_studio.mcp_server import create_server

    monkeypatch.delenv("STUDIO_API_TOKEN", raising=False)
    monkeypatch.delenv("STUDIO_MCP_ALLOW_DEVELOPMENT", raising=False)
    monkeypatch.setenv("STUDIO_MCP_SESSION_FILE", str(tmp_path / "missing.json"))

    async def check():
        async with Client(create_server()) as client:
            result = await client.call_tool("studio_identity", {})
            assert result.is_error and "mcp-login" in result.content[0].text
    asyncio.run(check())
