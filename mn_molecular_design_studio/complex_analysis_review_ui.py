"""Shared saved-base editor for the existing account-scoped human reviews."""
from __future__ import annotations

import hashlib

REVIEW_STATUSES = ["Unreviewed", "Needs review", "Reviewed", "Shortlisted for experimental testing",
                   "Not shortlisted for experimental testing", "Not applicable"]


def render_review_editor(user, identity, current, *, key_prefix, status_label, note_label, placeholder):
    import streamlit as st

    token = hashlib.sha256("\0".join(str(value or "") for value in (
        user.get("id"), identity["source_kind"], identity["subject_id"], identity["structure_ref"],
        identity.get("evaluation_id", ""))).encode()).hexdigest()[:16]
    status_key, note_key, base_key = (f"{key_prefix}-{part}-{token}" for part in ("status", "note", "base"))
    base = st.session_state.get(base_key)
    if base is None:
        base = dict(current)
        st.session_state[base_key] = base
    dirty = (st.session_state.get(status_key, base["status"]), st.session_state.get(note_key, base["note"])) != (base["status"], base["note"])
    if base.get("revision") != current.get("revision") and not dirty:
        base = dict(current)
        st.session_state[base_key] = base
        st.session_state[status_key], st.session_state[note_key] = base["status"], base["note"]
    st.session_state.setdefault(status_key, base["status"])
    st.session_state.setdefault(note_key, base["note"])
    if base.get("revision") != current.get("revision"):
        st.warning("This review changed elsewhere. Your unsaved draft is preserved; load the latest saved review before reconciling and saving.")
        def reload():
            for key in (status_key, note_key, base_key):
                st.session_state.pop(key, None)
        st.button("Load latest saved review", key=f"{key_prefix}-reload-{token}", on_click=reload)
    st.selectbox(status_label, REVIEW_STATUSES, key=status_key)
    st.text_area(note_label, key=note_key, max_chars=10_000, placeholder=placeholder)
    status, note = st.session_state[status_key], st.session_state[note_key]
    return {**identity, "status": status, "note": note, "expected_revision": base.get("revision")}, base_key, (status, note) != (base["status"], base["note"])


def mark_review_saved(base_key, saved):
    import streamlit as st
    st.session_state[base_key] = dict(saved)
