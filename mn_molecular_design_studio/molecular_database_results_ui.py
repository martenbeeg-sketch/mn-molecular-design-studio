"""Molecular Database actions for completed local run result pages."""
from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.molecular_database import MolecularDatabase, IMPORT_METHODS


TYPE_LABELS = {
    "protein": "Protein",
    "small_molecule": "Small molecule",
    "dna": "DNA",
    "rna": "RNA",
    "complex": "Complex",
}
log = logging.getLogger(__name__)


def render_run_comparison(run_dir, user, project_id=""):
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules
    from mn_molecular_design_studio.evaluations import evaluation_descriptor, project_observation
    from mn_molecular_design_studio.comparison_ui import render_comparison
    try:
        database = MolecularDatabase()
        records = list_run_molecules(run_dir, user, database=database)
    except (PermissionError, FileNotFoundError):
        return
    if not records or not any(record.get("observations") for record in records):
        return
    key = "run-comparison-" + hashlib.sha1(str(run_dir).encode()).hexdigest()[:12]
    with st.expander("Compare evaluations", expanded=False):
        all_models = st.checkbox("Include alternative models and poses", key=key + "-all")
        rows = []
        for record in records:
            observations = record.get("observations") or []
            if not all_models:
                chosen = (record.get("selected_evaluation") or {}).get("id")
                observations = [o for o in observations if evaluation_descriptor(o)["id"] == chosen]
            for observation in observations:
                scores, sources, _ = project_observation(observation)
                rows.append({"molecule_id": record["id"], "name": record["name"], "metadata": record.get("metadata", {}),
                             "components": record["components"], "evaluation": evaluation_descriptor(observation), "scores": scores, "score_sources": sources})
        render_comparison(rows, key=key, database=database, user=user, project_id=project_id)


def _record_label(record: dict) -> str:
    name = str(record.get("name") or record.get("id") or "Molecule")
    molecule_type = TYPE_LABELS.get(str(record.get("type") or ""), str(record.get("type") or "Molecule"))
    components = record.get("components") or []
    suffix = f" · {len(components)} chains" if len(components) > 1 else ""
    return f"{name} · {molecule_type}{suffix}"


def _source_group(records: list[dict]) -> dict | None:
    """Return the highest-priority output group explicitly linked to this run."""
    counts: dict[str, int] = {}
    groups: dict[str, dict] = {}
    for record in records:
        candidates = []
        for key in ("source_group", "run_group", "run_groups", "source_group_id", "run_group_id", "source_group_ids", "run_group_ids"):
            value = record.get(key)
            if value:
                candidates.extend(value if isinstance(value, list) else [value])
        visible_groups = record.get("groups") or []
        if isinstance(visible_groups, dict):
            visible_groups = [visible_groups]
        visible_by_id = {str(row.get("id")): row for row in visible_groups if isinstance(row, dict) and row.get("id")}
        seen = set()
        for item in candidates:
            if isinstance(item, str):
                item = {"id": item}
            if not isinstance(item, dict) or not item.get("id"):
                continue
            group_id = str(item["id"])
            if group_id not in seen:
                counts[group_id] = counts.get(group_id, 0) + 1
                seen.add(group_id)
            groups[group_id] = {**visible_by_id.get(group_id, {}), **groups.get(group_id, {}), **item, "id": group_id}
    if not counts:
        return None

    def priority(group: dict) -> int:
        origin = group.get("origin") or {}
        kind = group.get("collection_kind") or (
            origin.get("collection_kind") if isinstance(origin, dict) else None
        )
        # Prefer generated output structures, then evaluation collections.
        # Inputs remain a useful fallback for workflows with no output group.
        return {
            "generated": 0,
            "evaluation": 1,
            "curated": 2,
            "intermediate": 3,
            "inputs": 4,
        }.get(kind, 5)

    best_group_id = min(counts, key=lambda group_id: (priority(groups[group_id]), -counts[group_id]))
    return groups[best_group_id]


def _record_belongs_to_group(record: dict, group_id: str) -> bool:
    """Check the run's group associations, including compact ID-only forms."""
    for key in ("source_group", "run_group", "run_groups", "source_group_id", "run_group_id", "source_group_ids", "run_group_ids", "groups"):
        value = record.get(key)
        if not value:
            continue
        for item in value if isinstance(value, list) else [value]:
            if isinstance(item, dict):
                item = item.get("id")
            if str(item or "") == str(group_id):
                return True
    return False


def _collection_for_group(group: dict | None) -> str:
    if not group:
        return "all"
    origin = group.get("origin") or {}
    kind = group.get("collection_kind") or (origin.get("collection_kind") if isinstance(origin, dict) else None)
    method = group.get("import_method")
    return "library" if kind in {"curated", "generated"} or method in {"user_uploads", "agent_inputs"} else "all"


def _editable_destinations(database, user: dict) -> list[dict] | None:
    try:
        groups = database.groups(user, collection="library")
    except PermissionError:
        return None
    except Exception:
        log.exception("Could not load editable Database groups for a result page")
        return None
    editable = []
    for group in groups:
        origin = group.get("origin") or {}
        collection_kind = origin.get("collection_kind") if isinstance(origin, dict) else None
        is_curated_or_import = (
            collection_kind == "curated"
            or group.get("import_method") in {"user_uploads", "agent_inputs"}
        )
        is_owner = str(group.get("owner_id")) == str(user.get("id"))
        if is_curated_or_import and not group.get("read_only") and (is_owner or user.get("role") == "admin"):
            editable.append(group)
    return editable


def _project_picker(projects: list[dict], run_key: str, default_project: str = "") -> str:
    project_ids = [str(row["id"]) for row in projects if row.get("id")]
    default_project = str(default_project or "")
    options = ["", *project_ids]
    selected_default = default_project if default_project in project_ids else ""
    if default_project and default_project not in project_ids:
        st.caption("The run's project is not available here; choose an accessible project or a personal group.")
    key = f"run-database-project-{run_key}"
    if st.session_state.get(key) not in options:
        st.session_state[key] = selected_default
    if not project_ids:
        st.caption("No project is available; the new group will be personal.")
    return st.selectbox(
        "Project",
        options,
        key=key,
        format_func=lambda value: next((str(row["name"]) for row in projects if str(row.get("id")) == value), "Personal group"),
    )


def render_run_database_actions(run_dir, user: dict, projects: list[dict]) -> None:
    """Render selective promotion and database navigation for a local run."""
    run_path = Path(run_dir).expanduser().resolve()
    try:
        from mn_molecular_design_studio.molecular_database_results import list_run_molecules, save_run_selection

        database = MolecularDatabase()
        records = list_run_molecules(run_path, user, database=database)
    except PermissionError:
        st.caption("Molecular Database actions are unavailable for this run.")
        return
    except Exception:
        log.exception("Molecular Database actions could not load for run %s", run_path)
        st.caption("Molecular Database actions could not load for this run. Refresh the result page to try again.")
        return

    if not isinstance(records, list):
        return
    distinct = {}
    for record in records:
        if isinstance(record, dict) and record.get("id"):
            distinct.setdefault(str(record["id"]), record)
    records = list(distinct.values())
    if not records:
        return

    run_key = hashlib.sha1(str(run_path).encode()).hexdigest()[:12]
    record_by_id = {str(record["id"]): record for record in records}
    ids = list(record_by_id)
    selection_key = f"run-database-selection-{run_key}"
    if selection_key in st.session_state:
        st.session_state[selection_key] = [value for value in st.session_state[selection_key] if value in record_by_id]
    else:
        st.session_state[selection_key] = ids[:1] if len(ids) == 1 else []

    source_group = _source_group(records)
    if source_group:
        try:
            source_group = {**database.group(source_group["id"], user), **source_group}
        except PermissionError:
            source_group = None
        except Exception:
            log.exception("Could not verify source collection for run %s", run_path)
            st.caption("The source collection could not be verified; Database will open in All indexed.")
            source_group = None
    collection = _collection_for_group(source_group)
    st.markdown("#### Molecular Database")
    st.caption(f"{len(records):,} distinct molecular record{'s' if len(records) != 1 else ''} associated with this run.")
    open_col, context_col = st.columns([1.25, 4])
    if source_group:
        group_members = sum(
            any(isinstance(group, dict) and str(group.get("id")) == source_group["id"] for group in (record.get("groups") or []))
            for record in records
        )
        context_col.caption(f"Open source collection · {group_members:,} of {len(records):,} run records are in this group.")
    else:
        context_col.caption("Open All indexed and focus a molecular record from this run.")
    if open_col.button("Open in Database", key=f"run-database-open-{run_key}"):
        st.session_state["molecular-preselect-group"] = source_group.get("id") if source_group else None
        st.session_state["molecular-preselect-collection"] = collection
        focus_id = next(
            (
                str(record["id"])
                for record in records
                if source_group and _record_belongs_to_group(record, source_group["id"])
            ),
            ids[0],
        )
        st.session_state["molecular-preselect-molecule"] = focus_id
        if source_group and database._is_smoke_test_group(source_group):
            # Explicit navigation from this run should reveal its own source
            # collection, even when the browser normally hides smoke groups.
            st.session_state["molecular-hide-smoke-test-groups"] = False
        if not source_group:
            st.session_state["database-selected-group"] = ""
            st.session_state["molecular-project-filter"] = ""
            st.session_state["molecular-scope"] = "Organization"
            st.session_state["molecular-group-search"] = ""
            st.session_state["molecular-query"] = ""
            st.session_state["molecular-type-filter"] = ""
            st.session_state["molecular-score-filters"] = []
            for method in IMPORT_METHODS:
                st.session_state["molecular-method-" + method] = True
        st.session_state["page"] = "Database"
        st.rerun()

    save_message = None
    selected_ids = [value for value in st.session_state.get(selection_key, []) if value in record_by_id]
    with st.expander("Save selected to Database", expanded=False):
        st.caption(
            "Choose distinct molecular records. Structures, poses, scores, and their source run remain linked; "
            "poses are not separate molecules."
        )
        st.multiselect(
            "Molecules to save",
            options=ids,
            key=selection_key,
            format_func=lambda value: _record_label(record_by_id[value]),
            help="Only the selected molecular records are promoted.",
        )
        selected_ids = [value for value in st.session_state.get(selection_key, []) if value in record_by_id]
        st.caption(f"{len(selected_ids):,} selected of {len(records):,} molecular records.")

        destinations = _editable_destinations(database, user)
        groups_unavailable = destinations is None
        destinations = destinations or []
        if groups_unavailable:
            st.caption("Existing groups could not be loaded. You can still save to a new project group.")
        mode_options = (["Existing group"] if destinations else []) + ["New group"]
        destination_mode = st.radio(
            "Destination",
            mode_options,
            index=0,
            horizontal=True,
            key=f"run-database-destination-mode-{run_key}",
        )
        group_id = None
        group_name = "Saved results"
        project_id = None
        if destination_mode == "Existing group":
            project_names = {str(row.get("id")): str(row.get("name")) for row in projects or [] if row.get("id") and row.get("name")}

            def destination_label(value):
                group = next(group for group in destinations if str(group["id"]) == value)
                project_name = project_names.get(str(group.get("project_id") or ""))
                scope = f" · {project_name}" if project_name else ""
                return f"{group['name']}{scope} · {group.get('molecule_count', 0):,} molecules"

            group_id = st.selectbox(
                "Curated or imported group",
                [str(group["id"]) for group in destinations],
                key=f"run-database-destination-group-{run_key}",
                format_func=destination_label,
            )
            group_name = next(group["name"] for group in destinations if str(group["id"]) == group_id)
        else:
            default_project = str(source_group.get("project_id") or "") if source_group else ""
            group_name = st.text_input("New group name", value="Saved results", key=f"run-database-new-name-{run_key}")
            project_id = _project_picker(projects or [], run_key, default_project)

        disabled_reason = ""
        if not selected_ids:
            disabled_reason = "Select at least one molecule before saving."
        elif destination_mode == "New group" and not group_name.strip():
            disabled_reason = "Enter a name for the new group."
        if disabled_reason:
            st.caption(disabled_reason)
        if st.button(
            "Save selected molecules",
            type="primary",
            disabled=bool(disabled_reason),
            key=f"run-database-save-{run_key}",
        ):
            try:
                result = save_run_selection(
                    run_path,
                    user,
                    selected_ids,
                    database=database,
                    group_id=group_id,
                    name=group_name.strip() or "Saved results",
                    project_id=project_id,
                )
                saved_count = int(result.get("saved_count", len(selected_ids)))
                saved_group_id = str(result.get("group_id") or group_id or "")
                save_message = f"Saved {saved_count:,} selected molecule{'s' if saved_count != 1 else ''} to {group_name}."
                if saved_group_id:
                    st.session_state[f"run-database-saved-group-{run_key}"] = saved_group_id
                    st.session_state[f"run-database-saved-count-{run_key}"] = saved_count
            except (PermissionError, ValueError, OSError) as exc:
                st.error(f"Could not save the selected molecules: {exc}")
            except Exception:
                st.error("Could not save the selected molecules to the Database.")
    if save_message:
        st.success(save_message)
    saved_group_id = st.session_state.get(f"run-database-saved-group-{run_key}")
    if saved_group_id and st.button("Open saved group in Database", key=f"run-database-open-saved-{run_key}"):
        st.session_state["molecular-preselect-group"] = saved_group_id
        st.session_state["molecular-preselect-collection"] = "library"
        st.session_state["molecular-preselect-molecule"] = selected_ids[0] if selected_ids else ids[0]
        st.session_state["page"] = "Database"
        st.rerun()
