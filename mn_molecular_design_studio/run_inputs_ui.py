"""Run-output selection used by native tiles and graph source nodes."""
import streamlit as st


@st.cache_data(ttl=15, max_entries=32, show_spinner="Loading Run outputs…")
def options_for_user(project_id, user, token="", require_complex=True, stages=None, allow_docking=False):
    from mn_molecular_design_studio.logical_runs import input_options
    return input_options(project_id, user, token=token, require_complex=require_complex,
                         stages=stages, allow_docking=allow_docking)


def options(project_id, *, require_complex=True, stages=None):
    from mn_molecular_design_studio.auth import development_user, is_development_mode
    from mn_molecular_design_studio.logical_runs import input_options
    user = st.session_state.get("_molecular_database_user")
    if not user and is_development_mode():
        user = development_user()
    if not user:
        raise ValueError("The signed-in user context is unavailable.")
    return options_for_user(project_id, user, token=str(st.session_state.get("studio-session-token") or ""),
                            require_complex=require_complex, stages=stages)


def select(project_id, key, *, require_complex=True, stages=None, reference=None):
    if not project_id:
        st.info("Select a project to choose a Run output.")
        return None
    choices = options(project_id, require_complex=require_complex, stages=stages)
    if not choices:
        st.info("No compatible completed Run outputs are available in this project.")
        return None
    by_reference = {choice["reference"]: choice for choice in choices}
    labels = {reference: str(choice["label"]) for reference, choice in by_reference.items()}
    references = list(by_reference)
    selected_ref = st.selectbox("Run results", references, format_func=lambda value: labels[value],
        index=references.index(reference) if reference in references else 0, key=key + "-" + project_id)
    selected = by_reference[selected_ref]
    st.caption("Uses the entire output across its source workstations. Results are saved as a named output under this Run.")
    return selected
