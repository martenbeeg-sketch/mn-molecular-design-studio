"""Compact, accessible molecular sheet with independent focus and selection."""
from pathlib import Path

_component = None


def render_table(rows, columns, *, selected_ids, focus_id, key):
    import streamlit.components.v1 as components
    global _component
    if _component is None:
        _component = components.declare_component("studio_molecular_sheet", path=str(Path(__file__).parent / "components/molecule_table"))
    return _component(rows=rows, columns=columns, selected=selected_ids, focus=focus_id, key=key, default=None)
