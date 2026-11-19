"""Compatibility entry points for the shared Properties dialogs."""
from mn_molecular_design_studio.entity_properties_ui import render_project_properties


def render_project_rename(user, project, on_rename, *, read_only=False):
    render_project_properties(user, project, on_rename, read_only=read_only)


def render_target_rename(database, user, index, subject, *, key_prefix="target", read_only=False):
    if not subject or subject.get("role") not in {"Target", "Ligand"}:
        return
    from mn_molecular_design_studio.entity_names import component_id
    from mn_molecular_design_studio.entity_properties_ui import render_folder_properties
    identity = subject.get("identity")
    if not identity or len(identity) != 3 or not isinstance(identity[0], str):
        return
    subject_id = component_id({"type": identity[0], identity[1]: identity[2]})
    scope = {"id": subject["id"], "label": subject["name"], "subject": subject,
             "entity": {"kind": "target", "subject_id": subject_id}}
    render_folder_properties(database, user, index, scope, [], read_only=read_only, key_prefix=key_prefix)
