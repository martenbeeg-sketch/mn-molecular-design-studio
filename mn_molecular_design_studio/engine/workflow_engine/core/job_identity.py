"""Resolve the same job identity before dispatch, during execution, and on read."""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any


def job_tool(tool: object, params: Mapping[str, Any]) -> str:
    selected = str(tool or "")
    if selected == "optimized_prediction":
        engine = params.get("engine")
        if isinstance(engine, str) and re.fullmatch(r"[a-z][a-z0-9_]*", engine):
            return f"optimized_prediction_{engine}"
    return selected


def job_display_tool(tool: object, params: Mapping[str, Any]) -> str:
    """Return the user-facing engine identity without changing the runtime tool."""
    selected = str(tool or "")
    model_type = str(params.get("model_type") or "").strip().casefold()
    if selected == "ligandmpnn" and model_type == "ab_mpnn":
        return "abmpnn"
    return job_tool(selected, params)


def job_name(metadata: Mapping[str, Any], params: Mapping[str, Any]) -> str:
    return str(metadata.get("campaign_name") or params.get("campaign_name") or params.get("job_name") or "")
