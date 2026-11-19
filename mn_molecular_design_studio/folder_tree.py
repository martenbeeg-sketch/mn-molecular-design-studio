"""Selectable folder tree for the experimental SQL data browser."""
from pathlib import Path


_component = None


def render_tree(nodes, *, selected_id, expanded_ids, key):
    import streamlit.components.v1 as components
    global _component
    if _component is None:
        _component = components.declare_component("studio_data_folders", path=str(Path(__file__).parent / "components/folder_tree"))
    return _component(nodes=nodes, selectedId=selected_id, expandedIds=expanded_ids, key=key, default=None)
