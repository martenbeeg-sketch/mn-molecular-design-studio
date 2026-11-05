"""Public MCP screening import contract remains explicit and bounded."""

import asyncio


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
