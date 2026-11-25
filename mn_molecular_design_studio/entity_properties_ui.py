"""Contextual Properties dialogs shared by project and molecular browsers."""
import hashlib
import sqlite3
import uuid

import streamlit as st

from mn_molecular_design_studio.result_details_ui import render_result_details
from mn_molecular_design_studio import entity_names as names, entity_properties as properties

ERRORS = (ValueError, PermissionError, OSError, RuntimeError, sqlite3.Error)


def _key(value):
    return hashlib.sha256(str(value).encode()).hexdigest()[:24]


def _dismiss():
    active = st.session_state.pop("properties-active", None)
    if active:
        st.session_state.pop("properties-edit-" + _key([choice["id"] for choice in active["choices"]]), None)


@st.dialog("Properties", width="large", on_dismiss=_dismiss)
def _properties_dialog(user, choices, *, read_only=False):
    active = st.session_state.get("properties-active", {})
    if active.get("requires_session"):
        from mn_molecular_design_studio.auth import resolve_session
        current = resolve_session(st.session_state.get("studio-session-token"))
        if not current or current["id"] != user["id"]:
            _dismiss()
            st.error("Your session expired or changed. Sign in again before editing metadata.")
            return
        if current.get("role") != user.get("role"):
            _dismiss()
            st.error("Your account permissions changed. Reopen Properties from the current view.")
            return
        user = current
    dialog_key = _key([choice["id"] for choice in choices])
    edit_key = "properties-edit-" + dialog_key
    lookup = {choice["id"]: choice for choice in choices}
    choice_id = st.selectbox("Object to inspect", list(lookup), format_func=lambda value: lookup[value]["label"],
        key="properties-object-" + dialog_key, disabled=bool(st.session_state.get(edit_key)))
    choice = lookup[choice_id]
    try:
        snapshot = choice["load"](user)
    except ERRORS as exc:
        st.error("Properties could not be loaded: " + str(exc))
        return
    st.subheader(snapshot["name"])
    st.caption(snapshot["type_label"])
    st.info(snapshot["scope_label"])
    from mn_molecular_design_studio.cluster import config
    remote = snapshot.get("authority_node") and snapshot["authority_node"] != config.identity() and not config.is_worker()
    writable = snapshot["writable"] and bool(user.get("id")) and not user.get("disabled") and not read_only and not remote
    if not writable:
        st.caption("Read-only here. Editing requires ownership or administration rights and the metadata source workstation; cleanup locks also block editing.")
    state = st.session_state.get(edit_key)
    if state is None:
        for field, label in properties.FIELDS[snapshot["kind"]].items():
            if value := snapshot["properties"].get(field):
                st.markdown("**" + label + "**")
                st.text(value)
        if not any(snapshot["properties"].values()):
            st.caption("No descriptive annotations have been recorded yet.")
        if snapshot["aliases"]:
            st.markdown("**Previous names / aliases**")
            st.write("; ".join(snapshot["aliases"]))
        with st.expander("Identity and original source · read-only"):
            render_result_details(st, snapshot["source"])
        if st.button("Edit metadata", key="properties-edit-button-" + dialog_key, disabled=not writable, type="primary"):
            st.session_state[edit_key] = {"snapshot": snapshot, "instance": uuid.uuid4().hex}
            st.rerun(scope="app")
    else:
        base = state["snapshot"]
        form_key = "properties-form-" + state["instance"]
        st.caption("Only descriptive metadata is editable. Scientific identities, links, scores, parameters and pass/fail decisions are protected.")
        with st.form(form_key):
            name = st.text_input("Dataset entry label" if base["kind"] == "entry" else "Name", value=base["name"],
                                 max_chars=200, key=form_key + "-name", disabled=not writable)
            aliases = st.text_area("Additional aliases (one per line)", key=form_key + "-aliases", height=70, disabled=not writable)
            values = {}
            for field, label in properties.FIELDS[base["kind"]].items():
                if base["type_label"] == "Compound" and field == "description":
                    label = "Description"
                value = base["properties"].get(field, "")
                if field in {"notes", "description", "citation"}:
                    values[field] = st.text_area(label, value=value, key=form_key + "-" + field, height=90,
                        max_chars=12000 if field == "notes" else 4000, disabled=not writable)
                else:
                    values[field] = st.text_input(label, value=value, key=form_key + "-" + field, max_chars=4000, disabled=not writable)
            submitted = st.form_submit_button("Save changes", type="primary", disabled=not writable)
        if submitted:
            try:
                if not writable:
                    raise PermissionError("Editing is no longer authorized in this context.")
                choice["save"](base, name, values, [names.clean_name(alias) for alias in aliases.splitlines() if alias.strip()], actor=user)
            except ERRORS as exc:
                st.error(str(exc))
            else:
                _dismiss()
                st.session_state["properties-saved-notice"] = "Saved descriptive metadata for “" + name.strip() + "”. Scientific data is unchanged."
                st.session_state.pop("projects-inventory", None)
                st.rerun(scope="app")
        cancel, reload_column = st.columns(2)
        if cancel.button("Cancel editing", key=form_key + "-cancel"):
            st.session_state.pop(edit_key, None)
            st.rerun(scope="app")
        if reload_column.button("Reload current values", key=form_key + "-reload"):
            st.session_state[edit_key] = {"snapshot": snapshot, "instance": uuid.uuid4().hex}
            st.rerun(scope="app")
    if st.button("Close", key="properties-close-" + dialog_key):
        _dismiss()
        st.rerun(scope="app")


def _launch(user, choices, *, read_only=False):
    st.session_state.pop("properties-edit-" + _key([choice["id"] for choice in choices]), None)
    st.session_state["properties-active"] = {"user": user, "choices": choices, "read_only": read_only,
                                             "requires_session": bool(st.session_state.get("studio-session-token"))}
    _properties_dialog(user, choices, read_only=read_only)


def _project_choice(user, project, on_save=None):
    def load(actor):
        from mn_molecular_design_studio.workspace import load_organization
        current = next((row for row in load_organization().get("projects", []) if row.get("id") == project["id"]), project)
        return properties.project_snapshot(actor, current)

    def save(base, name, values, aliases, *, actor):
        callback = on_save or (lambda project_id, name, **kwargs: names.rename_project(actor, project_id, name, **kwargs))
        return callback(project["id"], name, expected_revision=base["revision"], properties=values, aliases=aliases,
                        expected_source_revision=base["source_revision"])

    return {"id": "project:" + project["id"], "label": "Project · " + project["name"], "load": load, "save": save}


def _group_choice(database, user, kind, group_id, *, molecule_id=None, outcome=None, label=""):
    def load(actor):
        if kind == "dataset":
            return properties.dataset_snapshot(database, actor, group_id)
        if kind == "entry":
            return properties.entry_snapshot(database, actor, group_id, molecule_id)
        return properties.result_snapshot(database, actor, outcome)

    def save(base, name, values, aliases, *, actor):
        return properties.save_group_properties(database, actor, base, name, values, aliases,
                                               expected_revision=base["revision"], outcome=outcome)

    return {"id": kind + ":" + group_id + ":" + str(molecule_id or (outcome or {}).get("id") or ""),
            "label": label, "load": load, "save": save}


def _target_choices(database, user, index, subject_id, *, group_ids=None):
    from mn_molecular_design_studio.database_browser import component_subjects
    subject = next((row for row in component_subjects(index, role="targets") if row["id"] == subject_id), None)
    if subject is None:
        return []
    groups = subject["group_ids"] & group_ids if group_ids is not None else subject["group_ids"]
    choices = []
    for scope in sorted({names.group_scope(index["groups"][group_id]) for group_id in groups}):
        label = ("Compound" if subject["type"] == "small_molecule" else "Target") + " · " + subject["name"]
        label += " · " + ("Personal library" if scope.startswith("owner:") else index.get("project_names", {}).get(scope.split(":", 1)[1], "Project"))
        def load(actor, scope=scope):
            return properties.target_snapshot(database, actor, subject_id, scope)

        def save(base, name, values, aliases, scope=scope, *, actor):
            from mn_molecular_design_studio.molecular_database import MolecularDatabase
            return names.rename_target(MolecularDatabase(database.root), actor, subject_id, scope, name,
                expected_revision=base["revision"], properties=values, aliases=aliases)

        choices.append({"id": "target:" + subject_id + ":" + scope, "label": label, "load": load, "save": save})
    return choices


def render_project_properties(user, project, on_save=None, *, read_only=False, key_prefix="project"):
    if st.button("Project properties", icon=":material/info:", key=key_prefix + "-properties-" + project["id"]):
        _launch(user, [_project_choice(user, project, on_save)], read_only=read_only)


def render_folder_properties(database, user, index, scope, projects, *, read_only=False, key_prefix="folder"):
    entity = scope.get("entity")
    if not entity:
        return
    choices = []
    if entity["kind"] == "project":
        project = next((row for row in projects if row["id"] == entity["project_id"]), None)
        if project:
            choices = [_project_choice(user, project)]
    elif entity["kind"] == "dataset":
        choices = [_group_choice(database, user, "dataset", entity["group_id"], label="Dataset · " + scope["label"])]
    elif entity["kind"] == "target":
        choices = _target_choices(database, user, index, entity["subject_id"], group_ids=scope["subject"]["group_ids"])
    if choices and st.button("Folder properties", icon=":material/info:", key=key_prefix + "-properties-" + scope["id"]):
        _launch(user, choices, read_only=read_only)


def render_record_properties(database, user, record, *, index=None, outcome=None, read_only=False, key_prefix="record"):
    key = key_prefix + "-properties-" + record["id"] + ((outcome or {}).get("id") or "")
    if not st.button("Result / molecule properties", icon=":material/info:", key=key):
        return
    if index is None:
        from mn_molecular_design_studio.database_browser import build_browser_index
        index = build_browser_index(database, user, collection="all", hide_smoke_test=False, include_related=True)
    if "project_names" not in index:
        from mn_molecular_design_studio.workspace import load_organization
        index["project_names"] = {row["id"]: row["name"] for row in load_organization().get("projects", [])}
    choices = []
    if outcome is not None:
        choices.append(_group_choice(database, user, "result", outcome["group_id"], outcome=outcome, label="Saved result / design"))
    seen = set()
    components = list(record["components"])
    metadata = record.get("metadata") or {}
    target_ids = {value for value in ((outcome or {}).get("target_id"), metadata.get("design_target_molecule_id"),
                  metadata.get("target_molecule_id")) if isinstance(value, str) and value}
    for target_id in target_ids:
        components.extend(index["molecules"].get(target_id, {}).get("components", []))
    for component in components:
        subject_id = names.component_id(component)
        if subject_id not in seen:
            choices.extend(_target_choices(database, user, index, subject_id))
            seen.add(subject_id)
    for group in record["groups"]:
        choices.append(_group_choice(database, user, "entry", group["id"], molecule_id=record["id"],
                                     label="Dataset entry · " + group["name"]))
    if choices:
        _launch(user, choices, read_only=read_only)


def render_notice():
    if notice := st.session_state.pop("properties-saved-notice", None):
        st.success(notice)
    if active := st.session_state.get("properties-active"):
        _properties_dialog(active["user"], active["choices"], read_only=active["read_only"])
