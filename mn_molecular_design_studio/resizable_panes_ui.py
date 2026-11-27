"""Shared resize controls for Streamlit workspaces with two or three panes."""
from __future__ import annotations

from pathlib import Path

import streamlit as st


_SPLITTER_JS = (Path(__file__).parent / "components/resizable_panes/resize.js").read_text(encoding="utf-8")


def resizable_columns(key, defaults, *, adjustable, toolbar=None):
    """Allocate percentages to adjustable panes and let the remaining pane fill the row.

    Each entry in adjustable maps a pane index to its widget key suffix, label,
    and width limits. The workspace key scopes its preferences to this layout.
    """
    widgets = {pane: (f"{key}-{suffix}-width", label, limits)
               for pane, (suffix, label, limits) in adjustable.items()}
    reset = {widget_key: defaults[pane] for pane, (widget_key, _, _) in widgets.items()}
    preferences_key = f"{key}-column-widths"
    for widget_key, width in st.session_state.get(preferences_key, reset).items():
        st.session_state.setdefault(widget_key, width)

    flexible = next(pane for pane in range(len(defaults)) if pane not in widgets)
    widths = list(defaults)
    controls = st if toolbar is None else toolbar
    with controls.popover("Column widths", icon=":material/view_column:"):
        for pane, (widget_key, label, (low, high)) in widgets.items():
            widths[pane] = st.slider(label, low, high, format="%d%%", key=widget_key)
        widths[flexible] = 100 - sum(widths[pane] for pane in widgets)
        st.caption(f"Remaining pane · {widths[flexible]}%")
        st.button("Reset widths", key=f"{key}-reset-widths", on_click=st.session_state.update, args=(reset,))
    st.session_state[preferences_key] = {widget_key: widths[pane] for pane, (widget_key, _, _) in widgets.items()}

    component_key = f"{key}-width-dividers"
    workspace_key = f"{key}-panes"

    def accept_resize():
        change = st.session_state[component_key].get("resize")
        if not isinstance(change, list) or len(change) != 2 or type(change[0]) is not int or change[0] not in widgets:
            return
        pane, width = change
        widget_key, _, (low, high) = widgets[pane]
        if type(width) is int and low <= width <= high:
            st.session_state[widget_key] = width

    splitter = st.components.v2.component(
        "studio_workspace_splitter",
        html="<span hidden></span>",
        isolate_styles=False,
        css="""
            .studio-pane-split { position: relative; }
            .studio-pane-divider {
                position: absolute; inset-block: 0; width: 16px; z-index: 2;
                transform: translateX(-50%); cursor: col-resize; touch-action: none;
            }
            .studio-pane-divider::before {
                content: ""; position: absolute; left: 7px; top: 0; bottom: 0; width: 2px;
                background: var(--st-border-color, #d7dedb); border-radius: 2px;
            }
            .studio-pane-divider::after {
                content: ""; position: absolute; left: 4px; top: 1.5rem; width: 8px; height: 1.5rem;
                background: var(--st-border-color, #b7c7bf); border-radius: 4px;
            }
            .studio-pane-divider:is(:hover, :focus-visible, :active)::before,
            .studio-pane-divider:is(:hover, :focus-visible, :active)::after {
                background: var(--st-primary-color, #287457);
            }
            .studio-pane-dragging { cursor: col-resize; user-select: none; }
        """,
        js=_SPLITTER_JS,
    )
    with st.container(key=workspace_key, gap=None):
        columns = st.columns(widths, gap="medium")
        splitter(key=component_key, height=1, on_resize_change=accept_resize, data={
            "workspace": workspace_key, "widths": widths, "flexible": flexible,
            "controls": {pane: {"label": label, "limits": limits} for pane, (_, label, limits) in widgets.items()},
        })
    return columns
