"""Locally bundled JSME editor; molecular drawings stay in the Studio session."""
from pathlib import Path

_component = None


def render_sketcher(*, key):
    import streamlit.components.v1 as components
    global _component
    if _component is None:
        _component = components.declare_component("studio_molecule_sketcher", path=str(Path(__file__).parent / "components/molecule_sketcher"))
    return _component(key=key, default=None)
