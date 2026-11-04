"""Interactive, self-contained presentation for the run-first analysis landing page."""
from __future__ import annotations

import hashlib
from pathlib import Path

_components = {}


def render_landing(*, rows: list[dict], view: str, selected_view: str,
                   view_state: dict, key: str) -> dict | None:
    import streamlit.components.v1 as components

    directory = Path(__file__).parent / "components" / "complex_analysis"
    digest = hashlib.sha256((directory / "index.html").read_bytes()).hexdigest()[:12]
    if digest not in _components:
        _components[digest] = components.declare_component(
            "studio_analysis_landing_" + digest, path=str(directory))
    return _components[digest](rows=rows, view=view, selected_view=selected_view,
                               view_state=view_state, key=key, default=None)
