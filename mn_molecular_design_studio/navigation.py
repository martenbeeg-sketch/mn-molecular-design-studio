"""Shared workspace navigation and project-scoped handoffs."""

import streamlit as st


WORKSPACE_PAGES = (
    ("Projects", ":material/folder_open:", "Projects"),
    ("Datasets", ":material/library_books:", "Datasets"),
    ("Pipelines", ":material/account_tree:", "Pipelines"),
    ("All tools", ":material/apps:", "All tools"),
    ("My results", ":material/query_stats:", "My runs"),
    ("Data browser", ":material/folder_copy:", "Data Browser"),
)


def render_workspace_navigation(user, *, development=False):
    if st.session_state.get("page") == "Database":
        st.session_state["page"] = "Data Browser"
    if st.session_state.get("page") == "Pipeline Builder":
        st.session_state["page"] = "Pipelines"
    st.markdown("""<style>
        [data-testid="stSidebar"] { background: #f7f8f7; }
        [data-testid="stSidebar"] [class*="st-key-nav-"] button {
            min-height: 2.65rem; justify-content: flex-start; border: 0;
            border-radius: 0.65rem; padding: 0.4rem 0.65rem;
            text-align: left; box-shadow: none;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] button p { text-align: left; }
        [data-testid="stSidebar"] [class*="st-key-nav-"] [data-testid="stBaseButton-tertiary"] {
            background: transparent; color: #46514b;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] [data-testid="stBaseButton-tertiary"]:hover {
            background: #edf1ee; color: #1f5d41;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] [data-testid="stBaseButton-primary"] {
            background: #e5ece7; color: #1c6342; font-weight: 600;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] [data-testid="stBaseButton-primary"]:hover {
            background: #dce7df; color: #174f35;
        }
    </style>""", unsafe_allow_html=True)
    pages = list(WORKSPACE_PAGES)
    if user.get("role") == "admin" and not development:
        pages.append(("Users & access", ":material/manage_accounts:", "Users & access"))
    pages.append(("Account", ":material/person:", "Account"))
    for label, icon, page in pages:
        current = st.session_state.get("page", "Projects")
        active = current == page or (page == "My runs" and current == "Job details") or (
            page == "All tools" and current == "Workflow setup"
        )
        if st.button(label, key="nav-" + label.lower().replace(" ", "-"), icon=icon,
                     type="primary" if active else "tertiary", use_container_width=True):
            st.session_state["page"] = page
            if page == "All tools":
                st.session_state.pop("selected_workflow", None)
                st.session_state.pop("selected_workflow_group", None)
            st.rerun()
        if page == "All tools":
            st.divider()


def open_workflow_tool(state, workflow):
    state.pop("selected_workflow_group", None)
    state.pop("selected_workflow_engine", None)
    state.update({"page": "Workflow setup", "selected_workflow": workflow})


def open_project_results(state, project_id):
    state.update({
        "page": "My runs", "results-project-filter": project_id,
        "cluster-results-project": project_id, "results-search": "",
        "results-status-filter": [], "results-tool-filter": [],
        "cluster-results-query": "", "cluster-results-machine": "",
        "cluster-results-status-filter": [], "results-page": 0,
        "cluster-results-show-removed": False,
        "cluster-results-page": 1, "results-selected-jobs": {},
        "cluster-results-selected": {}, "cluster-results-reset-selection": True,
    })
    for key in list(state):
        if str(key).startswith("results-row-select-") or key in {
            "cluster-results-pending-remove-and-purge", "results-pending-remove-and-purge", "cluster-results-filter-state",
            "results-filter-state", "results-export-data", "results-export-count",
        }:
            state.pop(key, None)


def open_project_database(state, project_id):
    for key in list(state):
        if str(key).startswith(("molecular-preselect-", "database-results-", "database-outcomes-")) or (
            str(key).startswith("database-browser-") and key != "database-browser-cache"
        ):
            state.pop(key, None)
        elif str(key).startswith("molecular-method-"):
            state[key] = True
    state.update({"page": "Data Browser", "molecular-project-filter": project_id,
                  "database-browser-view": "Molecules", "molecular-scope": "All",
                  "molecular-collection-view": "Library", "data-browser-open": {"project_id": project_id}})


def open_datasets(state, *, project_id="", group_id="", molecule_id=""):
    state.update({"page": "Datasets", "datasets-project": project_id, "datasets-open-group": group_id,
        "datasets-focus-molecule": molecule_id, "datasets-query": "", "datasets-member-query": ""})
    state.pop("datasets-page-selection", None)


def request_dataset_import(state, kind="", *, project_id="", return_to=""):
    if state.get("page") != "Datasets":
        open_datasets(state, project_id=project_id)
    state.update({"molecular-import-open": True, "molecular-import-step": 0,
                  "molecular-import-kind": kind})
    for key in ("molecular-import-destination", "molecular-import-group", "molecular-import-group-name",
                "molecular-import-method", "molecular-import-project", "molecular-import-fetched", "molecular-import-prepared"):
        state.pop(key, None)
    if return_to:
        state["dataset-import-return"] = return_to


def open_dataset_results(state, group_id, molecule_id=""):
    state.update({"page": "Data Browser", "molecular-preselect-group": group_id,
        "molecular-preselect-molecule": molecule_id, "molecular-preselect-collection": "library",
        "dataset-results-request": {"group_id": group_id, "molecule_id": molecule_id}})


def database_destination(state):
    destination = state.pop("data-browser-open", None)
    group_id = state.pop("molecular-preselect-group", None)
    molecule_id = state.pop("molecular-preselect-molecule", None)
    collection = state.pop("molecular-preselect-collection", None)
    if group_id or molecule_id:
        destination = {"group_id": group_id, "molecule_id": molecule_id, "collection": collection or "library"}
    if destination is not None:
        for key in list(state):
            if str(key).startswith(("data-browser-activity-", "data-browser-construct-")) or key in {
                "data-browser-page-selection", "data-browser-scope-path", "data-browser-expanded",
                "data-browser-source-group", "binderscope-browser-preflight",
                "data-browser-dataset-context", "data-browser-foldseek-selection",
            }:
                state.pop(key, None)
        state.update({"data-browser-result-query": "", "data-browser-folder-query": "",
                      "data-browser-type": "", "data-browser-coordinates": False,
                      "data-browser-indexed-history": destination.get("collection") == "all"})
    return destination
