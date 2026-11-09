"""Molecular Database actions for completed local run result pages."""
from __future__ import annotations

import hashlib
import json
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


def _render_saved_run_collections(run_path, user, database):
    from mn_molecular_design_studio.molecular_database_results import list_saved_run_collections
    try:
        groups = list_saved_run_collections(run_path, user, database=database)
    except (FileNotFoundError, PermissionError):
        return
    for group in groups:
        if st.button("Open collection in Complex Analysis · " + group["name"], key="run-saved-collection-" + group["id"]):
            st.session_state["complex-analysis-selected-collections"] = [group["id"]]
            st.session_state["complex-analysis-collections-open"] = True
            st.session_state["complex-analysis-individual-mode"] = False
            st.session_state["complex-analysis-flow-step"] = "select"
            st.session_state["page"] = "Complex Analysis"
            st.rerun()


def _candidate_ids_on_record(record: dict) -> set[tuple[str, str]]:
    if record.get("result_only"):
        candidate_id = str(record.get("candidate_id") or "")
        model_id = str(record.get("model_id") or "")
        if candidate_id:
            return {(candidate_id, model_id)}
    result = set()
    structures = record.get("structures") or record.get("representations") or []
    for structure in structures:
        context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
        candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
        candidate_id = str(context.get("candidate_id") or context.get("native_candidate_id")
                           or candidate.get("candidate_id") or "")
        model_id = str(context.get("model_id") or context.get("model_index") or context.get("validated_model_id")
                       or candidate.get("model_id") or candidate.get("model_index") or "")
        if candidate_id:
            result.add((candidate_id, model_id))
    return result


def _accepted_records(records: list[dict], candidates: list[dict]) -> tuple[list[str], list[dict]]:
    """Map explicit passing candidate/model decisions to exact run records."""
    selected = []
    unresolved = []
    for candidate in candidates:
        candidate_id = str(candidate.get("candidate_id") or "")
        model_id = str(candidate.get("model_id") or "")
        matches = []
        for record in records:
            pairs = _candidate_ids_on_record(record)
            if (candidate_id, model_id) in pairs or (not model_id and any(item[0] == candidate_id for item in pairs)):
                matches.append(str(record.get("id") or ""))
        matches = [value for value in matches if value]
        if matches:
            selected.extend(matches)
        else:
            unresolved.append(candidate)
    return list(dict.fromkeys(selected)), unresolved


def _collection_target_labels(passing: list[dict], records: list[dict], record_ids: list[str]) -> list[str]:
    """Use explicit target provenance from the pass decision or its selected Results record."""
    labels = [str(row.get("target_identity") or "").strip() for row in passing
              if str(row.get("target_identity") or "").strip()
              and str(row.get("target_identity")) != "Not recorded"]
    by_id = {str(record.get("id") or ""): record for record in records if isinstance(record, dict)}
    for record_id in record_ids:
        record = by_id.get(str(record_id)) or {}
        metadata = record.get("metadata") if isinstance(record.get("metadata"), dict) else {}
        target = str(metadata.get("target_identity") or "").strip()
        if target and target != "Not recorded":
            labels.append(target)
            continue
        sequences = metadata.get("target_sequences_by_chain")
        if isinstance(sequences, dict) and sequences:
            labels.extend(
                f"Target sequence (chain {chain} · {len(str(sequence))} aa)"
                for chain, sequence in sorted(sequences.items()) if str(sequence).strip()
            )
    return list(dict.fromkeys(labels))


def _render_run_passing_collection(run_path: Path, user: dict, database, records: list[dict], run_key: str) -> None:
    """Offer one explicit run-level import of candidates passing a recorded stage."""
    from mn_molecular_design_studio.molecular_database_results import run_acceptance_inventory, save_run_selection

    try:
        inventory = run_acceptance_inventory(run_path, user)
    except (FileNotFoundError, PermissionError, OSError, ValueError, TypeError) as exc:
        st.caption(f"Acceptance-based collection is unavailable for this Results run: {exc}")
        return
    if not inventory.get("candidate_count"):
        return

    st.markdown("#### Collect passing designs into Molecular Database")
    st.caption(f"Run contains {inventory['candidate_count']:,} generated candidate(s) across {inventory['model_count']:,} candidate/model row(s). Passing collection uses only decisions explicitly saved on this run.")
    stages = list(inventory.get("stages") or [])
    if not stages:
        st.info("This run has no recorded acceptance stage. No designs will be treated as passing. Use the advanced custom Results selection only if you deliberately want to collect unfiltered designs.")
        return

    stage_by_id = {str(row["stage_id"]): row for row in stages}
    final_present = "final" in stage_by_id
    placeholder = "Choose a recorded stage explicitly"
    stage_ids = (["final"] if final_present else []) + [str(row["stage_id"]) for row in stages if row["stage_id"] != "final"]
    options = [*stage_ids] if final_present else ["", *stage_ids]
    stage_key = f"run-collect-acceptance-stage-{run_key}"
    if st.session_state.get(stage_key) not in options:
        st.session_state[stage_key] = "final" if final_present else ""
    stage_id = st.selectbox(
        "Acceptance stage to collect",
        options,
        key=stage_key,
        format_func=lambda value: placeholder if not value else stage_by_id[value]["label"],
        help="A native engine pass remains a native-stage decision. It is not relabeled as final computational acceptance.",
    )
    if not stage_id:
        available = []
        for row in stages:
            available.append(
                f"{row['label']}: {int(row.get('passed_design_count') or 0):,} passing design(s) "
                f"({len(row.get('passed') or []):,} model rows), "
                f"{int(row.get('failed_model_count') or 0):,} failed, "
                f"{int(row.get('pending_model_count') or 0):,} unresolved, "
                f"{int(row.get('unrecorded_model_count') or 0):,} without a decision"
            )
        st.caption("Available recorded stage totals — " + "; ".join(available) + ".")
        st.info("No final decision is recorded for this run. Choose the exact earlier stage you want to use; native acceptance is not final computational acceptance.")
        return
    stage = stage_by_id[stage_id]
    passing = list(stage.get("passed") or [])
    passed_designs = int(stage.get("passed_design_count") or 0)
    st.caption(
        f"Stage: {stage['label']} · {passed_designs:,} passing design(s) across {len(passing):,} passing model row(s) · "
        f"{int(stage.get('failed_model_count') or 0):,} failed · {int(stage.get('pending_model_count') or 0):,} unresolved · "
        f"{int(stage.get('unrecorded_model_count') or 0):,} without a decision at this stage."
    )
    if not passing:
        return
    record_ids, unresolved = _accepted_records(records, passing)
    if unresolved:
        st.warning(f"{len(unresolved):,} passing candidate/model row(s) do not have a matching accessible Results structure in this run. They will not be silently substituted or collected.")
    if not record_ids:
        st.info("No passing model has a transferable run record and structure.")
        return

    try:
        from mn_molecular_design_studio.cluster import config as cluster_config
        cluster_settings = cluster_config.load()
        source_workstation = cluster_config.display_name(cluster_settings.get("name")) or "This workstation"
        source_node_id = str(cluster_config.identity())
    except (OSError, RuntimeError, ValueError):
        source_workstation, source_node_id = "This workstation", ""
    targets = sorted(_collection_target_labels(passing, records, record_ids))
    generators = sorted({str(row.get("generator_engine") or "") for row in passing
                         if row.get("generator_engine")})
    source_tools = sorted({str(row.get("source_tool") or "") for row in passing if row.get("source_tool")})
    decision_rows = [dict(item.get("decision") or {}) for item in passing]
    policies = list({json.dumps(row.get("policy") or {}, sort_keys=True): row.get("policy") or {}
                     for row in decision_rows}.values())
    collection_metadata = {
        "acceptance_stage": stage["label"],
        "acceptance_stage_id": stage_id,
        "acceptance_policies": policies,
        "accepted_candidate_models": [
            {"candidate_id": row["candidate_id"], "model_id": row.get("model_id") or "",
             "decision": row.get("decision") or {}, "source_tool": row.get("source_tool") or ""}
            for row in passing
        ],
        "target_identities": targets,
        "target_identity": targets[0] if len(targets) == 1 else (targets if targets else "Not recorded"),
        "generators": generators,
        "generator": generators[0] if len(generators) == 1 else (generators if generators else "Not recorded"),
        "source_tools": source_tools,
        "source_workstation": source_workstation,
        "source_node_id": source_node_id,
        "generated_candidate_count": int(inventory["candidate_count"]),
        "collected_design_count": passed_designs,
        "collected_model_count": len(passing),
    }
    for accepted, inventory_row in zip(collection_metadata["accepted_candidate_models"], passing):
        accepted["target_identity"] = inventory_row.get("target_identity") or ""
        accepted["target_sequences_by_chain"] = inventory_row.get("target_sequences_by_chain") or {}
    target_sequence_maps = {
        json.dumps(row.get("target_sequences_by_chain") or {}, sort_keys=True)
        for row in passing
    }
    if len(target_sequence_maps) == 1:
        collection_metadata["target_sequences_by_chain"] = json.loads(next(iter(target_sequence_maps)))
    target_label = targets[0] if len(targets) == 1 else "mixed targets" if targets else "target not recorded"
    generator_label = ", ".join(generators) if generators else "Design run"
    group_name = f"{generator_label} · {target_label} · {inventory['run_id']} · {stage['label']}"[:250]
    button_key = f"run-collect-passing-{run_key}-{stage_id}"
    if st.button("Collect passing designs into Molecular Database", type="primary",
                 key=button_key, disabled=bool(unresolved)):
        try:
            result = save_run_selection(
                run_path, user, record_ids, database=database, name=group_name,
                project_id=inventory.get("project_id") or None,
                candidate_models=passing, collection_metadata=collection_metadata,
            )
            group_id = str(result.get("group_id") or "")
            st.session_state[f"run-collect-result-{run_key}"] = {
                "group_id": group_id, "count": int(result.get("saved_count") or 0), "name": group_name,
                "design_count": passed_designs, "model_count": len(passing),
            }
        except (PermissionError, ValueError, OSError, FileNotFoundError) as exc:
            st.error(f"Could not collect the passing run designs: {exc}")
        except Exception:
            log.exception("Could not collect passing designs from run %s", run_path)
            st.error("Could not collect the passing run designs into the Molecular Database.")
    completed = st.session_state.get(f"run-collect-result-{run_key}")
    if isinstance(completed, dict) and completed.get("group_id"):
        st.success(f"Saved {completed['design_count']:,} passing design(s) and {completed['model_count']:,} model row(s) in {completed['name']} ({completed['count']:,} unique Database identity/identities). Repeating this collection reuses the same collection and scientific identities.")
        if st.button("Open collection in Complex Analysis", key=f"run-collect-open-{run_key}"):
            st.session_state["complex-analysis-selected-collections"] = [completed["group_id"]]
            st.session_state["complex-analysis-collections-open"] = True
            st.session_state["complex-analysis-individual-mode"] = False
            st.session_state["complex-analysis-flow-step"] = "select"
            st.session_state["page"] = "Complex Analysis"
            st.rerun()


def render_run_database_actions(run_dir, user: dict, projects: list[dict]) -> None:
    """Render selective promotion and database navigation for a local run."""
    run_path = Path(run_dir).expanduser().resolve()
    if (run_path / "artifacts/active_learning/docking_scores.parquet").is_file():
        _render_external_screening_results(run_path, user)
        return
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

    _render_saved_run_collections(run_path, user, database)

    if not isinstance(records, list):
        return
    distinct = {}
    for record in records:
        if isinstance(record, dict) and record.get("id"):
            distinct.setdefault(str(record["id"]), record)
    records = list(distinct.values())
    run_key = hashlib.sha1(str(run_path).encode()).hexdigest()[:12]
    _render_run_passing_collection(run_path, user, database, records, run_key)
    if not records:
        return

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
    from mn_molecular_design_studio.molecular_database_visibility import is_unreviewed_screening_shortlist
    collection = "all" if is_unreviewed_screening_shortlist(source_group) else _collection_for_group(source_group)
    indexed_records = [record for record in records if not record.get("result_only")]
    st.markdown("#### Molecular Database")
    st.caption(f"{len(records):,} distinct molecular record{'s' if len(records) != 1 else ''} associated with this run.")
    open_col, context_col = st.columns([1.25, 4])
    if not indexed_records:
        context_col.caption("These candidates remain in Results until selected for the curated Database.")
    elif source_group and is_unreviewed_screening_shortlist(source_group):
        context_col.caption("Original unreviewed computational shortlist · accessible in All indexed with its saved Results and Evaluations.")
    elif source_group:
        group_members = sum(
            any(isinstance(group, dict) and str(group.get("id")) == source_group["id"] for group in (record.get("groups") or []))
            for record in records
        )
        context_col.caption(f"Open source collection · {group_members:,} of {len(records):,} run records are in this group.")
    else:
        context_col.caption("Open All indexed and focus a molecular record from this run.")
    if indexed_records and open_col.button("Open in Database", key=f"run-database-open-{run_key}"):
        st.session_state["molecular-preselect-group"] = source_group.get("id") if source_group else None
        st.session_state["molecular-preselect-collection"] = collection
        focus_id = next(
            (
                str(record["id"])
                for record in records
                if source_group and _record_belongs_to_group(record, source_group["id"])
            ),
            str(indexed_records[0]["id"]),
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
    with st.expander("Advanced: collect a custom Results subset", expanded=False):
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


def _render_external_screening_results(run_path, user):
    """Preview external bulk tables without ingesting their molecular records."""
    from mn_molecular_design_studio.bulk_screening_results import external_results, import_final_results
    try:
        result = external_results(run_path, user)
    except (PermissionError, FileNotFoundError):
        return
    if result is None:
        return
    key = "external-screening-" + hashlib.sha1(str(run_path).encode()).hexdigest()[:12]
    with st.expander("External screening results and final Database shortlist", expanded=True):
        selection = result["selection"]
        st.caption(f"{result['compound_count']:,} scored compounds remain in external Parquet files. "
                   f"{selection.get('imported_count', 0):,} final compounds imported into Database.")
        st.code(result["scores_path"], language=None)
        for name, path in result["files"].items():
            if name != "docking_scores.parquet":
                st.caption(f"{name}: {path}")
        columns = ("compound_id", "name", "smiles", "score", "best_score", "first_round", "best_round", "pose_path")
        st.dataframe([{k: v for k, v in row.items() if k in columns} for row in result["rows"]],
                     hide_index=True, use_container_width=True)
        st.caption("Preview: first 100 compounds in the final ranking. The complete table retains all scored compounds for later investigation.")
        if selection.get("group_id") and st.button("Open final shortlist in Database", key=key + "-open"):
            st.session_state["molecular-preselect-group"] = selection["group_id"]
            from mn_molecular_design_studio.molecular_database_visibility import is_unreviewed_screening_shortlist_id
            historical = is_unreviewed_screening_shortlist_id(selection["group_id"])
            st.session_state["molecular-preselect-collection"] = "all" if historical else "library"
            st.session_state["molecular-hide-smoke-test-groups"] = False
            st.session_state["page"] = "Database"
            st.rerun()
        from mn_molecular_design_studio.molecular_database_visibility import is_unreviewed_screening_shortlist_id
        if is_unreviewed_screening_shortlist_id(selection.get("group_id")):
            st.caption("This is an original, unreviewed computational shortlist. Open it in All indexed to inspect its saved memberships and Evaluations.")
        if selection.get("error"):
            st.warning("The optional final Database import failed: " + str(selection["error"]))
        by_id = {}
        duplicate_ids = set()
        for row in result["rows"]:
            compound_id = str(row.get("compound_id") or "").strip()
            if not compound_id:
                continue
            if compound_id in by_id:
                duplicate_ids.add(compound_id)
            else:
                by_id[compound_id] = row
        selectable_ids = [value for value in by_id if value not in duplicate_ids]
        if duplicate_ids:
            st.caption("Some preview IDs repeat in the screening library and cannot be selected by ID. Use the reviewed top-ranked option below for those rows.")
        selected_ids = st.multiselect(
            "Select exact compounds to import",
            options=selectable_ids,
            key=key + "-selected-compounds",
            format_func=lambda value: (
                f"{value} · {by_id[value].get('compound_name') or by_id[value].get('name') or 'Compound'} · "
                f"rank {by_id[value].get('_final_rank', selectable_ids.index(value) + 1)}"
            ),
            help="Choose specific IDs from the ranked preview. The selection is saved separately from the automatic shortlist.",
        )
        st.caption("Only the top 100 ranked rows are shown here. Full scores, poses, and files remain available in Results.")
        if st.button("Import selected compounds into Database", key=key + "-import-selected", disabled=not selected_ids):
            try:
                imported = import_final_results(run_path, user, 0, selected_compound_ids=selected_ids)
                st.success(f"Selected {imported['imported_count']:,} compounds into the Database. Full results remain in Results.")
                st.rerun()
            except (ValueError, OSError, PermissionError) as exc:
                st.error(str(exc))

        count = int(st.number_input("Or import a reviewed top-ranked shortlist (0 = none)", min_value=0,
                                   max_value=1000, value=0, step=100, key=key + "-count"))
        st.caption("A positive count is an explicit choice to curate the top-ranked compounds. Repeating it reuses the same selection and preserves annotations.")
        if st.button("Import final shortlist into Database", key=key + "-import", disabled=count == 0):
            try:
                imported = import_final_results(run_path, user, count)
                st.success(f"Final shortlist: {imported['imported_count']:,} compounds in Database. Full results remain external.")
                st.rerun()
            except (ValueError, OSError, PermissionError) as exc:
                st.error(str(exc))
