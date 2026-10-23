"""AutoDock Vina's protein, pocket, and ligand input page."""

from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
import uuid
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import streamlit as st

from mn_molecular_design_studio import docking_preparations, vina_inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.hotspot_sequence import render_hotspot_sequence
from mn_molecular_design_studio.target_editor import TargetSource, _target_storage_path
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
    StructureVisualization, molstar_custom_component,
)


def _save_current_preparation() -> None:
    active = st.session_state.get("docking-active-preparation")
    if active and st.session_state.get("docking-preparation-current-project") == active["project_id"]:
        docking_preparations.save(**active, state=dict(st.session_state))


def _activate_preparation(record: dict, *, owner_id: str, project_id: str, workflow: str) -> None:
    previous = st.session_state.get("docking-active-preparation") or {}
    workflows = {workflow, previous.get("workflow", workflow)}
    for key in list(st.session_state):
        if any(docking_preparations.state_key(key, value) for value in workflows):
            st.session_state.pop(key, None)
        elif key.startswith("vina-sequence-event-") or key in {
            "vina-pending-placement", "vina-pocket-placement-node", "vina-last-target-upload",
            "vina-pending-residue-text", "vina-residue-error", "vina-residue-apply-error",
        }:
            st.session_state.pop(key, None)
    st.session_state.update(docking_preparations.snapshot(record["state"], workflow))
    st.session_state["vina-input-epoch"] = st.session_state.get("vina-input-epoch", 0) + 1
    st.session_state["docking-active-preparation"] = {
        "owner_id": owner_id, "project_id": project_id, "workflow": workflow, "preparation_id": record["id"],
    }


def _new_preparation(owner_id: str, project_id: str, workflow: str, choice_key: str) -> None:
    _save_current_preparation()
    record = docking_preparations.save(owner_id, project_id, workflow, str(uuid.uuid4()), {})
    st.session_state[choice_key] = record["id"]
    _activate_preparation(record, owner_id=owner_id, project_id=project_id, workflow=workflow)


def _render_preparation_picker(owner_id: str, project_id: str | None, workflow: str) -> None:
    st.session_state["docking-preparation-current-project"] = project_id
    if not project_id:
        return
    records = docking_preparations.list_saved(owner_id, project_id, workflow)
    if not records:
        # Preserve the open form when the first project preparation is created.
        values = dict(st.session_state) if not st.session_state.get("docking-active-preparation") else {}
        records = [docking_preparations.save(owner_id, project_id, workflow, str(uuid.uuid4()), values)]
    choice_key = "docking-preparation-choice-" + hashlib.sha256(f"{owner_id}:{project_id}:{workflow}".encode()).hexdigest()[:16]
    ids = [row["id"] for row in records]
    if st.session_state.get(choice_key) not in ids:
        st.session_state[choice_key] = ids[0]
    chooser, new = st.columns([5, 1], vertical_alignment="bottom")
    selected = chooser.selectbox(
        "Saved preparation", ids, key=choice_key,
        format_func=lambda value: next(f"{row['name']} · {row['id'][:6]}" for row in records if row["id"] == value),
        help="Protein, compounds, pocket settings, and the linked prediction are saved for this project on this Studio installation.",
    )
    active = st.session_state.get("docking-active-preparation") or {}
    if active.get("preparation_id") != selected:
        _activate_preparation(next(row for row in records if row["id"] == selected), owner_id=owner_id, project_id=project_id, workflow=workflow)
    else:
        # Streamlit removes widget values when their page is not rendered.
        record = next(row for row in records if row["id"] == selected)
        for key, value in docking_preparations.snapshot(record["state"], workflow).items():
            st.session_state.setdefault(key, value)
    new.button("New preparation", key="docking-new-preparation", on_click=_new_preparation,
               args=(owner_id, project_id, workflow, choice_key), use_container_width=True)
    st.caption("Changes are saved automatically. Reopen this tool, select this project and preparation to continue later.")


def _default_co_crystallized_ligand(path: Path):
    """Return the first bound ligand with a valid default Vina box."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
    from mn_molecular_design_studio.ligand_tools_ui import _vina_box_from_points

    try:
        choices = ligand_tools.bound_ligand_choices(path)
    except (OSError, ValueError):
        return None, None, None
    for ligand in choices:
        try:
            center, size = _vina_box_from_points(
                list(ligand["heavy_atom_coordinates"].values()), padding=5.0,
            )
        except ValueError:
            continue
        return ligand, center, size
    return None, None, None


def _load_example(owner_id: str, project_id: str = "", *, job_name_key: str = "ligand-vina-job-name", example_job_name: str = "autodock-vina-example") -> None:
    target, compounds = vina_inputs.stage_example(owner_id, project_id=project_id)
    example_ligand, center, size = _default_co_crystallized_ligand(target)
    st.session_state.update({
        "vina-target-path": str(target), "vina-compound-input": compounds,
        "vina-ligand-format": "SDF File", "vina-protein-format": "PDB",
        job_name_key: example_job_name,
        "ligand-vina-pocket-mode": "Co-crystallized Ligand" if example_ligand else "Manual",
        "ligand-vina-box-ligand": example_ligand["key"] if example_ligand else "",
        "ligand-vina-selected-residues": [],
        "ligand-vina-applied-residues": [],
        "vina-default-target-signature": hashlib.sha1(target.read_bytes()).hexdigest()[:12],
        "vina-derived-box-signature": None,
    })
    st.session_state.pop("vina-applied-residue-padding", None)
    center = center or dict.fromkeys("xyz", 0.0)
    size = size or dict.fromkeys("xyz", 20.0)
    for axis in "xyz":
        st.session_state[f"vina-box-center-{axis}"] = center[axis]
        st.session_state[f"vina-box-size-{axis}"] = size[axis]
    st.session_state["vina-input-epoch"] = st.session_state.get("vina-input-epoch", 0) + 1
    st.session_state.pop("vina-pending-prediction", None)
    st.session_state.pop("vina-applied-prediction", None)
    st.session_state.pop("vina-residue-error", None)


def _focus_vina_box() -> None:
    st.session_state["vina-viewer-focus-epoch"] = st.session_state.get("vina-viewer-focus-epoch", 0) + 1


def _select_vina_bound_ligand() -> None:
    """Move the box and viewer focus together when a bound ligand is changed."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
    from mn_molecular_design_studio.ligand_tools_ui import _vina_box_from_points

    target_path = Path(st.session_state.get("vina-target-path", ""))
    selected_key = st.session_state.get("ligand-vina-box-ligand", "")
    if not target_path.is_file() or not selected_key:
        return
    try:
        ligand = next(
            (row for row in ligand_tools.bound_ligand_choices(target_path) if row["key"] == selected_key),
            None,
        )
        if ligand is None:
            return
        padding = float(st.session_state.get("ligand-vina-ligand-padding", 5.0))
        center, size = _vina_box_from_points(list(ligand["heavy_atom_coordinates"].values()), padding)
    except (OSError, ValueError):
        return
    for axis in "xyz":
        st.session_state[f"vina-box-center-{axis}"] = center[axis]
        st.session_state[f"vina-box-size-{axis}"] = size[axis]
    st.session_state["vina-derived-box-signature"] = json.dumps([
        st.session_state.get("vina-default-target-signature"),
        "Co-crystallized Ligand", selected_key, padding,
    ])
    _focus_vina_box()


def _apply_vina_selected_residues() -> None:
    """Apply the pending residue selection to the docking box and viewer."""
    from mn_molecular_design_studio.ligand_tools_ui import (
        _vina_box_from_points, _vina_protein_atom_coordinates,
    )

    target_path = Path(st.session_state.get("vina-normalized-target-path", ""))
    selected = sorted(tuple(row) for row in st.session_state.get("ligand-vina-selected-residues", []))
    padding = float(st.session_state.get("ligand-vina-residue-padding", 10.0))
    if not target_path.is_file():
        st.session_state["vina-residue-apply-error"] = "Load a protein structure before applying a residue selection."
        return
    try:
        coordinates = _vina_protein_atom_coordinates(target_path.read_text())
        points = [point for chain, number in selected for point in coordinates.get((chain, number), [])]
        center, size = _vina_box_from_points(points, padding)
    except (OSError, ValueError) as exc:
        st.session_state["vina-residue-apply-error"] = str(exc)
        return
    for axis in "xyz":
        st.session_state[f"vina-box-center-{axis}"] = center[axis]
        st.session_state[f"vina-box-size-{axis}"] = size[axis]
    st.session_state["ligand-vina-applied-residues"] = selected
    st.session_state["vina-applied-residue-padding"] = padding
    st.session_state["vina-derived-box-signature"] = json.dumps([
        st.session_state.get("vina-default-target-signature"),
        "Select Residues", selected, padding,
    ])
    st.session_state.pop("vina-residue-apply-error", None)
    _focus_vina_box()


def _apply_vina_padding_to_confirmed_selection() -> None:
    """Recalculate the box when padding changes for an already-confirmed pocket."""
    selected = sorted(tuple(row) for row in st.session_state.get("ligand-vina-selected-residues", []))
    applied = sorted(tuple(row) for row in st.session_state.get("ligand-vina-applied-residues", []))
    if selected and selected == applied:
        _apply_vina_selected_residues()


def _use_structure(path: str) -> None:
    st.session_state["vina-target-path"] = path
    st.session_state.pop("vina-last-target-upload", None)
    st.session_state["vina-input-epoch"] = st.session_state.get("vina-input-epoch", 0) + 1


def _remove_structure() -> None:
    st.session_state.pop("vina-target-path", None)
    st.session_state.pop("vina-last-target-upload", None)
    st.session_state["vina-input-epoch"] = st.session_state.get("vina-input-epoch", 0) + 1


def _file_row(path: Path) -> None:
    st.markdown(
        f'<div class="vina-file"><span>{html.escape(path.name)}<small>{path.stat().st_size / 1024:.0f} KB</small></span>'
        '<span class="vina-upload-complete">Upload complete<small>Ready to use</small></span></div>',
        unsafe_allow_html=True,
    )


def _render_target_input(owner_id: str) -> Path | None:
    st.markdown('<p class="vina-label">Protein File</p>', unsafe_allow_html=True)
    file_format = st.segmented_control("Protein format", ["PDB", "CIF"], default="PDB",
                                       key="vina-protein-format", label_visibility="collapsed")
    with st.container(horizontal=True, vertical_alignment="center"):
        code = st.text_input("Enter PDB Code", placeholder="Enter PDB Code", key="vina-pdb-code",
                             label_visibility="collapsed", width=190).strip().upper()
        fetch = st.button("Fetch PDB", key="vina-fetch", disabled=not code, width=110)
    if fetch:
        try:
            if len(code) != 4 or not code.isalnum():
                raise ValueError("Enter a four-character PDB code.")
            suffix = ".cif" if file_format == "CIF" else ".pdb"
            try:
                with urlopen(Request(f"https://files.rcsb.org/download/{code}{suffix}", headers={"User-Agent": "MN-Molecular-Design-Studio"}), timeout=30) as response:
                    content = response.read(50_000_001)
            except HTTPError as exc:
                if suffix != ".pdb" or exc.code != 404:
                    raise
                suffix = ".cif"
                with urlopen(f"https://files.rcsb.org/download/{code}.cif", timeout=30) as response:
                    content = response.read(50_000_001)
            if len(content) > 50_000_000:
                raise ValueError("Protein files may be at most 50 MB.")
            path = _target_storage_path(content, code.lower() + suffix, owner_id=owner_id, allow_non_protein=True)
            _use_structure(str(path))
            st.rerun()
        except (OSError, ValueError, URLError) as exc:
            st.error(f"Could not fetch {code}: {exc}")
    with st.container(key="vina-protein-source"):
        upload_col, saved_col = st.columns([12, 1.4], gap="small")
        uploaded = upload_col.file_uploader("Protein file", label_visibility="collapsed",
            type=["cif", "mmcif"] if file_format == "CIF" else ["pdb", "ent"],
            key=f"vina-protein-upload-{st.session_state.get('vina-input-epoch', 0)}")
        with saved_col.popover("My proteins", use_container_width=True):
            saved = sorted(user_input_root(owner_id).rglob("*.pdb"), key=lambda p: p.stat().st_mtime, reverse=True)[:100]
            if saved:
                chosen = st.selectbox("Saved protein inputs", saved, format_func=lambda p: p.name, key="vina-saved-target")
                st.button("Use structure", key="vina-use-target", on_click=_use_structure, args=(str(chosen),))
            else:
                st.caption("Uploaded and fetched proteins will appear here.")
        if uploaded is not None:
            digest = hashlib.sha256(uploaded.getvalue()).hexdigest()
            if digest != st.session_state.get("vina-last-target-upload"):
                try:
                    path = _target_storage_path(uploaded.getvalue(), uploaded.name, owner_id=owner_id, allow_non_protein=True)
                    st.session_state["vina-target-path"] = str(path)
                    st.session_state["vina-last-target-upload"] = digest
                except (ValueError, OSError) as exc:
                    st.error(str(exc))
        raw_path = st.session_state.get("vina-target-path")
        if raw_path:
            path = Path(raw_path).resolve()
            if path.is_file() and path.is_relative_to(user_input_root(owner_id)):
                with upload_col:
                    _file_row(path)
                saved_col.button("Remove", help="Remove protein input", key="vina-remove-target", on_click=_remove_structure, use_container_width=True)
                return path
    return None


def _set_residues(chain: str, residues: set[int]) -> None:
    selected = {tuple(row) for row in st.session_state.get("ligand-vina-selected-residues", [])}
    selected = {row for row in selected if row[0] != chain} | {(chain, number) for number in residues}
    st.session_state["ligand-vina-selected-residues"] = sorted(selected)
    st.session_state[f"vina-residue-text-{chain}"] = vina_inputs.residue_range_text(residues)
    st.session_state.pop("vina-residue-error", None)


def _typed_residues(chain: str, allowed: set[int]) -> None:
    try:
        _set_residues(chain, vina_inputs.parse_residue_ranges(st.session_state.get(f"vina-residue-text-{chain}", ""), allowed))
    except ValueError as exc:
        st.session_state["vina-residue-error"] = str(exc)


def _render_residue_picker(sequence_by_chain: dict) -> float:
    chain_col, padding_col = st.columns([3, 1])
    chain = chain_col.selectbox("Chain", list(sequence_by_chain), key="vina-pocket-chain",
        format_func=lambda key: f"Chain {'(blank)' if key == '_' else key} ({len(sequence_by_chain[key])} residues)")
    padding = float(padding_col.number_input(
        "Box padding (Å)", value=10.0, min_value=0.0, max_value=30.0, step=1.0,
        key="ligand-vina-residue-padding", on_change=_apply_vina_padding_to_confirmed_selection,
    ))
    st.caption("Click and drag or type indices below to select residues to include. Selected residues are highlighted in green.")
    selected = {int(number) for selected_chain, number in st.session_state.get("ligand-vina-selected-residues", []) if selected_chain == chain}
    sequence = sequence_by_chain[chain]
    event = render_hotspot_sequence(chain=chain, sequence=sequence, selected_residues=selected,
        key=f"vina-sequence-{st.session_state.get('vina-default-target-signature')}-{chain}", compact=True)
    if isinstance(event, dict):
        signature = json.dumps(event, sort_keys=True)
        if st.session_state.get(f"vina-sequence-event-{chain}") != signature:
            st.session_state[f"vina-sequence-event-{chain}"] = signature
            requested = {int(n) for n in event.get("residues", [])} & {n for n, _ in sequence}
            if requested != selected:
                _set_residues(chain, requested)
                st.rerun()
    st.session_state.setdefault(f"vina-residue-text-{chain}", vina_inputs.residue_range_text(selected))
    st.text_input("Selected residues", placeholder="e.g. 225, 230-235", key=f"vina-residue-text-{chain}",
                  on_change=_typed_residues, args=(chain, {n for n, _ in sequence}))
    st.caption(f"Valid range: {sequence[0][0]}–{sequence[-1][0]}")
    with st.container(horizontal=True):
        st.button("Select all", key="vina-select-all", on_click=_set_residues, args=(chain, {n for n, _ in sequence}))
        st.button("Clear", key="vina-clear-residues", disabled=not selected, on_click=_set_residues, args=(chain, set()))
    if st.session_state.get("vina-residue-error"):
        st.error(st.session_state["vina-residue-error"])
    pending = sorted(tuple(row) for row in st.session_state.get("ligand-vina-selected-residues", []))
    applied = sorted(tuple(row) for row in st.session_state.get("ligand-vina-applied-residues", []))
    applied_padding = st.session_state.get("vina-applied-residue-padding")
    confirmed = pending == applied and applied_padding is not None and float(applied_padding) == padding
    if confirmed:
        st.caption("The docking box uses the confirmed residues. Changing padding recalculates it automatically.")
    else:
        st.caption("Confirm the selected residues and padding to update the docking box and viewer.")
    if st.session_state.get("vina-residue-apply-error"):
        st.error(st.session_state["vina-residue-apply-error"])
    st.button(
        "Confirm selection and update box", key="vina-apply-selected-residues",
        on_click=_apply_vina_selected_residues,
        disabled=not pending or bool(st.session_state.get("vina-residue-error")),
        use_container_width=True,
    )
    return padding


def _active_compound_input(compound_sets: list[dict]) -> dict | None:
    if st.session_state.get("vina-ligand-format") == "Database compound set":
        selected_id = str(st.session_state.get("ligand-tool-compound-set-ligand.docking.vina") or "")
        return next(
            (row for row in compound_sets if str(row.get("group_id")) == selected_id),
            compound_sets[0] if compound_sets else None,
        )
    selected = st.session_state.get("vina-compound-input")
    return selected if isinstance(selected, dict) else None


def _boltz_probe_defaults(compound_input: dict | None) -> list[dict[str, str]]:
    """Prefill Boltz probes from the selected imported compound input."""
    rows: list[dict[str, str]] = []
    raw_path = (compound_input or {}).get("path")
    if raw_path:
        path = Path(raw_path)
        try:
            stat = path.stat()
            source_id = str((compound_input or {}).get("group_id") or (compound_input or {}).get("source_sha256") or path.resolve())
            for index, (compound_id, compound_name, smiles) in enumerate(_compound_preview_rows(str(path), stat.st_mtime_ns, stat.st_size)[:10]):
                display_name = compound_name or compound_id
                rows.append({
                    "id": f"imported:{source_id}:{index}:{compound_id}", "source": "imported", "chain": "",
                    "label": f"Imported compound {compound_id} · {display_name}",
                    "compound_id": compound_id, "compound_name": display_name, "smiles": smiles,
                })
        except (OSError, ValueError, json.JSONDecodeError):
            pass
    return rows


def _sync_boltz_probes(defaults: list[dict[str, str]]) -> list[dict[str, str]]:
    """Keep editable probes in sync with selected-chain and compound-set defaults."""
    old_rows = st.session_state.get("vina-probe-rows")
    saved_values = dict(st.session_state.get("vina-probe-values") or {})
    dismissed = set(st.session_state.get("vina-probe-dismissed") or ())
    if isinstance(old_rows, list):
        for index, row in enumerate(old_rows):
            if isinstance(row, dict):
                saved_values[str(row.get("id") or "")] = str(
                    st.session_state.get(f"vina-probe-{index}", row.get("smiles") or "")
                )
                row["compound_id"] = str(st.session_state.get(f"vina-probe-id-{index}", row.get("compound_id") or ""))
                row["compound_name"] = str(st.session_state.get(f"vina-probe-name-{index}", row.get("compound_name") or ""))
    else:
        old_rows = []

    previous_by_id = {
        str(row.get("id") or ""): row for row in old_rows if isinstance(row, dict)
    }
    rows: list[dict[str, str]] = []
    for default in defaults[:10]:
        probe_id = str(default["id"])
        if probe_id in dismissed:
            continue
        row = dict(default)
        old = previous_by_id.get(probe_id)
        row["smiles"] = (
            saved_values.get(probe_id, row["smiles"])
            if old is None else saved_values.get(probe_id, str(old.get("smiles") or row["smiles"]))
        )
        if old is not None:
            row["compound_id"] = str(old.get("compound_id") or row.get("compound_id") or "")
            row["compound_name"] = str(old.get("compound_name") or row.get("compound_name") or "")
        rows.append(row)

    for old in old_rows:
        if not isinstance(old, dict) or old.get("source") != "manual":
            continue
        probe_id = str(old.get("id") or "")
        if probe_id and probe_id not in dismissed:
            rows.append({**old, "smiles": saved_values.get(probe_id, str(old.get("smiles") or ""))})

    if "vina-probe-initialized" not in st.session_state and not rows:
        rows.append({"id": "manual:0", "source": "manual", "chain": "", "label": "Probe",
                     "compound_id": "", "compound_name": "", "smiles": ""})
        st.session_state["vina-probe-next-id"] = 1
    if len(rows) > 10:
        rows = rows[:10]

    # Reindex after a chain is deselected, while retaining edited values by stable ID.
    for index in range(10):
        st.session_state.pop(f"vina-probe-{index}", None)
        st.session_state.pop(f"vina-probe-id-{index}", None)
        st.session_state.pop(f"vina-probe-name-{index}", None)
    for index, row in enumerate(rows):
        st.session_state[f"vina-probe-{index}"] = row["smiles"]
        st.session_state[f"vina-probe-id-{index}"] = row.get("compound_id") or ""
        st.session_state[f"vina-probe-name-{index}"] = row.get("compound_name") or ""
        saved_values[row["id"]] = row["smiles"]
    st.session_state["vina-probe-rows"] = rows
    st.session_state["vina-probe-values"] = saved_values
    st.session_state["vina-probe-dismissed"] = sorted(dismissed)
    st.session_state["vina-probe-initialized"] = True
    return rows


def _add_boltz_probe() -> None:
    rows = list(st.session_state.get("vina-probe-rows") or [])
    if len(rows) < 10:
        next_id = int(st.session_state.get("vina-probe-next-id", 0))
        rows.append({"id": f"manual:{next_id}", "source": "manual", "chain": "", "label": "Probe",
                     "compound_id": "", "compound_name": "", "smiles": ""})
        st.session_state["vina-probe-rows"] = rows
        st.session_state["vina-probe-next-id"] = next_id + 1


def _remove_boltz_probe(index: int) -> None:
    rows = list(st.session_state.get("vina-probe-rows") or [])
    if not 0 <= index < len(rows):
        return
    values = [
        str(st.session_state.get(f"vina-probe-{current}", row.get("smiles") or ""))
        for current, row in enumerate(rows)
    ]
    identities = [
        (
            str(st.session_state.get(f"vina-probe-id-{current}", row.get("compound_id") or "")),
            str(st.session_state.get(f"vina-probe-name-{current}", row.get("compound_name") or "")),
        )
        for current, row in enumerate(rows)
    ]
    removed = rows.pop(index)
    values.pop(index)
    identities.pop(index)
    probe_id = str(removed.get("id") or "")
    if removed.get("source") != "manual" and probe_id:
        dismissed = set(st.session_state.get("vina-probe-dismissed") or ())
        dismissed.add(probe_id)
        st.session_state["vina-probe-dismissed"] = sorted(dismissed)
    saved_values = dict(st.session_state.get("vina-probe-values") or {})
    saved_values.pop(probe_id, None)
    for current in range(10):
        st.session_state.pop(f"vina-probe-{current}", None)
        st.session_state.pop(f"vina-probe-id-{current}", None)
        st.session_state.pop(f"vina-probe-name-{current}", None)
    for current, (row, value, (compound_id, compound_name)) in enumerate(zip(rows, values, identities)):
        row["smiles"] = value
        row["compound_id"] = compound_id
        row["compound_name"] = compound_name
        st.session_state[f"vina-probe-{current}"] = value
        st.session_state[f"vina-probe-id-{current}"] = compound_id
        st.session_state[f"vina-probe-name-{current}"] = compound_name
        saved_values[str(row.get("id") or "")] = value
    st.session_state["vina-probe-rows"] = rows
    st.session_state["vina-probe-values"] = saved_values


def _render_boltz_pocket(
    target_text: str,
    project_id: str | None,
    job_name: str,
    submit,
    *,
    docking_engine: str = "AutoDock Vina",
    target_path: Path | None = None,
    compound_input: dict | None = None,
    placement_options=None,
    prediction_state=None,
    active: bool = True,
) -> bool:
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
    from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction

    _parse_proteins = optimized_prediction._parse_proteins

    sequences_by_chain = vina_inputs.protein_sequence_entries(target_text)
    available_chains = list(sequences_by_chain)
    if available_chains:
        chain_selector_options = {
            "default": available_chains,
        } if "vina-boltz-selected-chains" not in st.session_state else {}
        selected_chains = st.multiselect(
            "Protein chains for Boltz",
            available_chains,
            format_func=lambda chain: f"Chain {'(blank)' if chain == '_' else chain} · {len(sequences_by_chain[chain])} residues",
            key="vina-boltz-selected-chains",
            help="Select one chain or a complex of chains. Boltz uses these sequences, and predicted contacts are mapped back onto the same chains in the imported PDB.",
            **chain_selector_options,
        )
        fasta = "\n\n".join(
            f">Chain {'(blank)' if chain == '_' else chain}\n" + "".join(letter for _, letter in sequences_by_chain[chain])
            for chain in selected_chains
        )
        st.markdown("**Selected chain sequence(s)**")
        st.code(fasta)
        additional_sequence_label = "Additional protein sequence(s) (optional)"
        st.caption("The imported PDB remains the receptor for docking. Extra sequences are prediction context; only selected PDB chains can define the transferred pocket.")
    else:
        selected_chains = []
        additional_sequence_label = "Protein sequence(s)"
        st.warning("No protein-chain sequence was recognized in this PDB. Add a valid sequence and ensure it corresponds to a chain in the imported structure so Boltz contacts can map back to the docking receptor.")
    additional_sequences_text = st.text_area(
        additional_sequence_label,
        placeholder="Paste one protein sequence per line. Separate chains with new lines or colons.",
        key="vina-boltz-additional-sequences",
        height=100,
    )
    sequence_error = ""
    try:
        additional_sequences = _parse_proteins(additional_sequences_text, required=False)
    except ValueError as exc:
        additional_sequences = []
        sequence_error = str(exc)
        st.error(sequence_error)
    protein_sequence_rows = [
        {"chain": chain, "sequence": "".join(letter for _, letter in sequences_by_chain[chain])}
        for chain in available_chains if chain in selected_chains
    ]
    selected_proteins = [row["sequence"] for row in protein_sequence_rows] + [row["sequence"] for row in additional_sequences]
    proteins = "\n".join(selected_proteins)

    sync_probes = active or bool(st.session_state.get("vina-probe-initialized"))
    defaults = _boltz_probe_defaults(compound_input) if sync_probes else []
    probe_rows = _sync_boltz_probes(defaults) if sync_probes else []
    if active and len(defaults) > 10:
        st.caption("Showing the first 10 compounds from the selected input set.")
    st.caption(
        "Boltz predicts probe binding poses from the selected protein-chain sequence(s). "
        f"The probe pose is aligned back onto the imported PDB to locate the pocket; "
        f"{docking_engine} then docks against the imported receptor."
    )
    probe_values = []
    probe_identities = []
    for index, row in enumerate(probe_rows):
        columns = st.columns([2.4, 2.8, 5, 1.2])
        compound_id = columns[0].text_input(
            "Compound ID",
            value=str(row.get("compound_id") or ""),
            key=f"vina-probe-id-{index}",
            help="Stable identifier retained in results. Leave blank to use the compound name.",
        )
        compound_name = columns[1].text_input(
            "Compound name *",
            value=str(row.get("compound_name") or ""),
            key=f"vina-probe-name-{index}",
            help="Required; shown with the ID in prediction results and tables.",
        )
        value = columns[2].text_input(
            f"{row['label']} · SMILES",
            key=f"vina-probe-{index}",
            help="Imported compound SMILES are prefilled when available. Edit the SMILES or remove this probe.",
        )
        columns[3].button(
            "Remove", key=f"vina-remove-probe-{index}",
            on_click=_remove_boltz_probe, args=(index,),
        )
        probe_values.append(value)
        probe_identities.append({
            "compound_id": compound_id.strip() or compound_name.strip(),
            "compound_name": compound_name.strip(),
            "smiles": value.strip(),
        })
    st.button("+ Add Probe", disabled=len(probe_rows) >= 10, key="vina-add-probe", on_click=_add_boltz_probe)
    probes = [value.strip() for value in probe_values if value.strip()]
    st.caption(f"{len(probes)}/10 probes")
    distance = float(st.number_input("Contact distance (Å)", value=5.0, min_value=1.0, max_value=20.0, step=1.0, key="vina-probe-distance"))
    missing_requirements = []
    if not project_id:
        missing_requirements.append("select a project")
    if not selected_chains:
        missing_requirements.append("select at least one protein chain")
    if not probes:
        missing_requirements.append("provide at least one probe SMILES")
    if any(probe["smiles"] and not probe["compound_name"] for probe in probe_identities):
        missing_requirements.append("enter a compound name for each probe")
    if sequence_error:
        missing_requirements.append("correct the protein sequence input")
    if missing_requirements:
        st.caption(f"To enable pocket prediction, {', '.join(missing_requirements)}.")
    prediction_context = {
        "target_hash": hashlib.sha256(target_text.encode()).hexdigest(),
        "request_signature": _boltz_request_signature(target_text, proteins, probes, distance, job_name, probe_identities),
        "distance": distance, "project_id": str(project_id or ""),
        "selected_chains": [row["chain"] for row in protein_sequence_rows],
        "protein_sequences": selected_proteins, "probe_smiles": probes,
    }
    # Adopt a prediction submitted before durable preparations were introduced.
    previous_job = st.session_state.get("last_submitted_job") or {}
    if (not st.session_state.get("vina-pending-prediction")
            and previous_job.get("project_id") == project_id
            and previous_job.get("campaign_name") == f"{job_name}-pocket"
            and "boltz2" in str(previous_job.get("tool"))):
        _remember_boltz_prediction(prediction_context)
    predict_clicked = st.button(
        "Predict Binding Site", type="primary", disabled=bool(missing_requirements),
        help=(f"To run Boltz, {', '.join(missing_requirements)}." if missing_requirements else "Run Boltz-2 to predict probe poses and transfer their contacts to the docking box."),
        key="vina-predict-pocket",
    )
    if predict_clicked:
        st.session_state.pop("vina-pending-placement", None)
        st.session_state.pop("vina-pocket-placement-node", None)
        try:
            from rdkit import Chem
            if any(Chem.MolFromSmiles(value) is None for value in probes):
                raise ValueError("Every probe must be a valid SMILES string.")
            if not proteins:
                raise ValueError("Select at least one chain from the imported PDB or add a protein sequence.")
            parameters = optimized_prediction.boltz2_prediction_parameters(
                [
                    {"name": f"probe_{i + 1}", "protein_sequences": proteins,
                     "compound_id": probe_identities[i]["compound_id"],
                     "compound_name": probe_identities[i]["compound_name"],
                     "ligand_smiles": value}
                    for i, value in enumerate(probes)
                ],
                settings={
                    "optimization_mode": "off", "seed": 0, "samples": 1,
                    "num_recycles": optimized_prediction.ENGINE_PROFILES["boltz2"]["default_recycles"],
                    "cpu_cores": 4, "gpu_device": "0",
                    "use_msa": False, "use_msa_server": False, "generate_ipsae": False,
                },
                job_name=f"{job_name}-pocket",
            )
            schema = workflow_schema("refolding.optimized_prediction")
            idempotency_key = str(uuid.uuid4())
            options = placement_options(schema, "refolding.optimized_prediction", parameters, str(project_id), idempotency_key) if placement_options else {}
            candidates = options.get("candidates", []) if isinstance(options, dict) else []
            if len(candidates) > 1:
                st.session_state.pop("vina-pocket-placement-node", None)
                st.session_state["vina-pending-placement"] = {
                    "schema": schema, "workflow": "refolding.optimized_prediction",
                    "parameters": parameters, "project_id": str(project_id),
                    "idempotency_key": idempotency_key, "candidates": candidates,
                    "prediction_context": prediction_context,
                }
            else:
                previous = _submit_boltz_prediction(
                    submit, schema, parameters, str(project_id),
                    options.get("placement") if isinstance(options, dict) else None,
                    idempotency_key,
                )
                _remember_boltz_prediction(prediction_context, previous=previous)
        except (ValueError, OSError) as exc:
            st.error(f"Could not start pocket prediction: {exc}")

    pending_placement = st.session_state.get("vina-pending-placement")
    if isinstance(pending_placement, dict):
        pending_context = pending_placement.get("prediction_context", {})
        current_signature = _boltz_request_signature(target_text, proteins, probes, distance, job_name, probe_identities)
        if pending_context.get("request_signature") != current_signature:
            st.warning("The protein, probe, or pocket settings changed. Start prediction again to refresh workstation availability.")
            if st.button("Dismiss workstation choices", key="vina-dismiss-stale-placement"):
                st.session_state.pop("vina-pending-placement", None)
                st.session_state.pop("vina-pocket-placement-node", None)
                st.rerun()
        else:
            candidates = pending_placement.get("candidates", [])

            def _candidate_label(candidate):
                gpus = candidate.get("free_gpus", [])
                gpu_text = ", ".join(
                    f"GPU {gpu.get('gpu_id')} · {float(gpu.get('free_vram_gb', 0)):.1f} GB VRAM"
                    for gpu in gpus
                ) or "CPU only"
                return f"{candidate.get('name', 'Workstation')} · {gpu_text} · {float(candidate.get('cpu_free', 0)):g} CPU threads"

            st.info("Automatic placement found more than one compatible workstation. Choose where this Boltz-2 prediction should run.")
            selected_node = st.selectbox(
                "Workstation", [candidate["node_id"] for candidate in candidates],
                format_func=lambda node_id: _candidate_label(next(item for item in candidates if item["node_id"] == node_id)),
                key="vina-pocket-placement-node",
            )
            st.caption("Availability is checked again on submission. If the selected workstation becomes busy, the job will wait there.")
            confirm_col, cancel_col = st.columns([1, 1])
            if confirm_col.button("Predict on selected workstation", type="primary", key="vina-confirm-pocket-placement"):
                try:
                    previous = st.session_state.get("last_submitted_job")
                    _submit_boltz_prediction(
                        submit, pending_placement["schema"], pending_placement["parameters"],
                        pending_placement["project_id"], selected_node, pending_placement["idempotency_key"],
                    )
                    if st.session_state.get("last_submitted_job") != previous:
                        _remember_boltz_prediction(pending_context, previous=previous)
                        st.session_state.pop("vina-pending-placement", None)
                except (ValueError, OSError) as exc:
                    st.error(f"Could not start pocket prediction: {exc}")
            if cancel_col.button("Cancel", key="vina-cancel-pocket-placement"):
                st.session_state.pop("vina-pending-placement", None)
                st.rerun()
    pending = st.session_state.get("vina-pending-prediction")
    if (pending and pending.get("target_hash") == prediction_context["target_hash"]
            and pending.get("project_id", str(project_id)) == str(project_id)):
        if (pending.get("protein_sequences") != selected_proteins
                or pending.get("probe_smiles", probes) != probes):
            st.warning("The chain or probe selection changed. Restore the prediction inputs or start a new prediction.")
            return False
        _render_prediction_progress(pending, target_text, distance, prediction_state)
        return st.session_state.get("vina-applied-prediction") == _prediction_application_signature(pending, distance)
    return False


def _prediction_application_signature(pending: dict, distance: float) -> str:
    job = pending.get("job") or {}
    identity = job.get("cluster_job_id") or job.get("run_id") or pending.get("run_dir")
    return json.dumps([identity, distance, pending.get("selected_chains") or []])


@st.fragment(run_every="10s")
def _render_prediction_progress(pending: dict, target_text: str, distance: float, prediction_state) -> None:
    try:
        if prediction_state is None:
            st.info("Open Job Results to follow the saved prediction.")
            return
        result = prediction_state(pending.get("job") or {"run_dir": pending.get("run_dir")})
        job = result["job"]
        pending["job"] = job
        st.session_state["vina-pending-prediction"] = pending
        status = str(job.get("status") or "queued")
        label = job.get("job_code") or job.get("run_id") or str(job.get("cluster_job_id") or "")[:8]
        machine = f" on {result['workstation']}" if result.get("workstation") else ""
        title = f"Boltz-2 · Job {label} · {status.replace('_', ' ').capitalize()}{machine}"
        if status in {"failed", "cancelled", "canceled", "interrupted"}:
            st.error(title)
            st.caption(str(job.get("error") or job.get("wait_reason") or "See Job Results for the prediction log."))
        elif status in {"completed", "succeeded"}:
            st.success(title)
        else:
            st.info(title)
            if job.get("wait_reason"):
                st.caption(job["wait_reason"])
        st.caption("Status refreshes every 10 seconds. You can leave this page; the job and saved preparation will remain available.")
        st.button("Refresh prediction", key="vina-refresh-pocket")
        signature = _prediction_application_signature(pending, distance)
        if status not in {"completed", "succeeded"} or st.session_state.get("vina-applied-prediction") == signature:
            return
        structures = result.get("structures") or []
        if not structures:
            st.info("Loading prediction structures…" if result.get("loading") else "The completed prediction has no usable complex structures. See Job Results for details.")
            return
        residues, overlay = set(), ""
        selected_chains = pending.get("selected_chains") or []
        for structure_text in structures:
            predicted = docking_preparations.pdb_text(structure_text)
            residues |= vina_inputs.predicted_contact_residues(target_text, predicted, distance, selected_chains=selected_chains)
            if not overlay:
                overlay = vina_inputs.predicted_ligand_overlay_pdb(target_text, predicted, selected_chains=selected_chains)
        if not residues:
            raise ValueError("No protein residues contact the predicted probes at this distance.")
        st.session_state["ligand-vina-selected-residues"] = sorted(residues)
        st.session_state["vina-pending-residue-text"] = {
            chain: vina_inputs.residue_range_text({number for selected, number in residues if selected == chain})
            for chain in vina_inputs.protein_sequence_entries(target_text)
        }
        st.session_state["vina-prediction-overlay-pdb"] = overlay
        st.session_state["vina-applied-prediction"] = signature
        _save_current_preparation()
        st.rerun()
    except (ValueError, OSError, LookupError) as exc:
        st.error(f"Could not update the binding-site prediction: {exc}")


def _boltz_request_signature(
    target_text: str, proteins: str, probes: list[str], distance: float, job_name: str,
    probe_identities: list[dict[str, str]] | None = None,
) -> str:
    value = json.dumps([target_text, proteins, probes, probe_identities or [], float(distance), job_name], separators=(",", ":"))
    return hashlib.sha256(value.encode()).hexdigest()


def _submit_boltz_prediction(submit, schema, parameters, project_id, placement, idempotency_key):
    previous = st.session_state.get("last_submitted_job")
    settings = {"_execution": {"placement": placement, "idempotency_key": idempotency_key}} if placement else {}
    submit(
        schema, "refolding.optimized_prediction", parameters, project_id,
        settings,
        "Boltz binding-site prediction",
    )
    return previous


def _remember_boltz_prediction(context: dict, *, previous=None) -> None:
    job = st.session_state.get("last_submitted_job")
    if isinstance(job, dict) and job != previous and (job.get("cluster_job_id") or job.get("run_dir")):
        st.session_state["vina-pending-prediction"] = {"run_dir": job.get("run_dir"), "job": job, **context}
        st.session_state.pop("vina-applied-prediction", None)
        _save_current_preparation()


def _use_database_set(group_id: str) -> None:
    st.session_state["vina-ligand-format"] = "Database compound set"
    st.session_state["ligand-tool-compound-set-ligand.docking.vina"] = group_id


def _remove_ligand_input() -> None:
    st.session_state.pop("vina-compound-input", None)
    st.session_state["vina-input-epoch"] = st.session_state.get("vina-input-epoch", 0) + 1


def _set_ligand_preview_index(state_key: str, index: int) -> None:
    st.session_state[state_key] = index


@st.cache_data(show_spinner=False, max_entries=32)
def _compound_preview_rows(path: str, modified_ns: int, file_size: int) -> tuple[tuple[str, str, str], ...]:
    """Read compound IDs, names, and SMILES for ligand previews and probes."""
    del modified_ns, file_size  # Included in the cache key to invalidate changed files.
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    compounds = payload.get("compounds") or []
    return tuple(
        (
            str(row.get("compound_id") or f"Compound {index}"),
            str(row.get("compound_name") or (row.get("properties") or {}).get("name") or row.get("compound_id") or f"Compound {index}"),
            str(row.get("smiles") or ""),
        )
        for index, row in enumerate(compounds, start=1)
    )


def _render_ligand_2d_preview(compound_set: dict) -> None:
    raw_path = compound_set.get("path")
    if not raw_path:
        return
    path = Path(raw_path)
    if not path.is_file():
        return
    try:
        stat = path.stat()
        compounds = _compound_preview_rows(
            str(path), stat.st_mtime_ns, stat.st_size,
        )
    except (OSError, ValueError, json.JSONDecodeError):
        return
    count = len(compounds)
    if not count:
        return
    preview_key = hashlib.sha1(str(path.resolve()).encode("utf-8")).hexdigest()[:12]
    state_key = f"vina-ligand-preview-index-{preview_key}"
    index = min(max(0, int(st.session_state.get(state_key, 0))), count - 1)
    st.session_state[state_key] = index
    if count > 1:
        previous, position, following = st.columns([1, 2, 1])
        previous.button(
            "← Previous", key=f"vina-ligand-preview-prev-{preview_key}",
            disabled=index == 0,
            on_click=_set_ligand_preview_index, args=(state_key, index - 1),
        )
        position.caption(f"Compound {index + 1:,} of {count:,}")
        following.button(
            "Next →", key=f"vina-ligand-preview-next-{preview_key}",
            disabled=index >= count - 1,
            on_click=_set_ligand_preview_index, args=(state_key, index + 1),
        )
    compound_id, compound_name, smiles = compounds[index]
    if not smiles:
        return
    from mn_molecular_design_studio.cofolding_input import _smiles_preview

    preview = _smiles_preview(smiles)
    label = html.escape(compound_id if compound_id == compound_name else f"{compound_id} · {compound_name}")
    if not preview:
        st.caption(f"2D structure preview is unavailable for {label} in the configured chemistry environment.")
        return
    st.markdown(
        '<div style="margin: 14px 0; padding: 12px; border: 1px solid #dce3df; border-radius: 8px; background: #fff;">'
        f'<div style="font-weight: 600; margin-bottom: 8px;">2D structure preview · {label}</div>'
        f'<img src="{preview}" alt="2D molecular structure preview for {label}" '
        'style="display:block; width:240px; max-width:100%; height:auto; margin:auto;" />'
        '</div>',
        unsafe_allow_html=True,
    )


def _render_ligand_input(owner_id: str, project_id: str | None, compound_sets: list[dict]) -> dict | None:
    fmt = st.selectbox("Ligand Format", ["SDF File", "SMILES (CSV)", "Database compound set"], key="vina-ligand-format", width=300,
        help="Upload SDF coordinates or a CSV containing a SMILES column, or choose a compound set from Studio Database.")
    if fmt == "Database compound set":
        saved_selection = st.session_state.get("vina-compound-input") or {}
        if (str(saved_selection.get("group_id") or "").startswith("selection:")
                and Path(saved_selection.get("path") or "").is_file()
                and saved_selection["group_id"] not in {row["group_id"] for row in compound_sets}):
            compound_sets = [*compound_sets, saved_selection]
        if not compound_sets:
            st.caption("Import compounds in Database, or choose SDF File or SMILES (CSV) above.")
            return None
        group_id = st.selectbox("My compounds", [str(row["group_id"]) for row in compound_sets],
            format_func=lambda value: next(f"{row['name']} · {row['compound_count']:,} compounds" for row in compound_sets if str(row["group_id"]) == value),
            key="ligand-tool-compound-set-ligand.docking.vina")
        selected_set = next(row for row in compound_sets if str(row["group_id"]) == group_id)
        _render_ligand_2d_preview(selected_set)
        return selected_set
    st.markdown('<p class="vina-label">Ligand File</p>', unsafe_allow_html=True)
    st.caption("SDF structure file for your ligand" if fmt == "SDF File" else "CSV file with a SMILES column and optional compound names")
    with st.container(key="vina-ligand-source"):
        upload_col, database_col = st.columns([12, 1.4], gap="small")
        uploaded = upload_col.file_uploader("Ligand File", label_visibility="collapsed", type=["sdf", "mol"] if fmt == "SDF File" else ["csv", "tsv", "smi"],
            key=f"vina-ligand-upload-{fmt}-{st.session_state.get('vina-input-epoch', 0)}", accept_multiple_files=False)
        with database_col.popover("My compounds", use_container_width=True):
            if compound_sets:
                value = st.selectbox("Database compound set", [str(row["group_id"]) for row in compound_sets],
                    format_func=lambda value: next(str(row["name"]) for row in compound_sets if str(row["group_id"]) == value), key="vina-compound-popover")
                st.button("Use compounds", key="vina-use-compounds", on_click=_use_database_set, args=(value,))
            else:
                st.caption("Your Database compound sets will appear here.")
        if uploaded is not None:
            try:
                digest = hashlib.sha256(uploaded.getvalue()).hexdigest()[:20]
                if (st.session_state.get("vina-compound-input") or {}).get("source_sha256") != digest:
                    st.session_state["vina-compound-input"] = vina_inputs.stage_compound_input(owner_id, uploaded.name, uploaded.getvalue(), project_id=str(project_id or ""))
            except (ValueError, OSError) as exc:
                st.error(f"Could not read ligand input: {exc}")
                return None
        selected = st.session_state.get("vina-compound-input")
        if selected:
            source = Path(selected["source_path"])
            if source.is_file() and ((fmt == "SDF File") == (source.suffix.lower() in {".sdf", ".mol"})):
                with upload_col:
                    _file_row(source)
                database_col.button("Remove", key="vina-remove-ligand", on_click=_remove_ligand_input, use_container_width=True)
                _render_ligand_2d_preview(selected)
                return selected
    return None


def render_docking_setup(workflow_id, schema, *, owner_id, compound_sets, project_selector, submit, title, placement_options=None, prediction_state=None):
    from mn_molecular_design_studio.ligand_tools_ui import (
        _vina_target_as_pdb, _vina_protein_atom_coordinates, _selected_vina_protein_residues,
        _vina_box_from_points, _set_vina_derived_box, _save_reference_ligand_upload,
    )
    tool_key = {
        "ligand.docking.vina": "vina",
        "ligand.docking.autodock_gpu": "autodock-gpu",
        "ligand.docking.gnina": "gnina",
        "ligand.docking.unidock_pro": "unidock-pro",
        "ligand.docking.rosetta_ligand": "rosetta",
    }.get(workflow_id)
    if tool_key is None:
        raise ValueError(f"The shared docking page does not support {workflow_id}.")
    is_vina = workflow_id == "ligand.docking.vina"
    is_autodock_gpu = workflow_id == "ligand.docking.autodock_gpu"
    is_gnina = workflow_id == "ligand.docking.gnina"
    is_unidock = workflow_id == "ligand.docking.unidock_pro"
    is_rosetta = workflow_id == "ligand.docking.rosetta_ligand"
    compound_limit_max = 5000
    if is_unidock:
        from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import MAX_COMPOUNDS_PER_LAUNCH

        compound_limit_max = MAX_COMPOUNDS_PER_LAUNCH
    rosetta_readiness = None
    if is_rosetta:
        from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools

        rosetta_readiness = ligand_tools.rosettaligand_runner_readiness()
    defaults = {
        "vina": "autodock-vina", "autodock-gpu": "autodock-gpu", "gnina": "gnina",
        "unidock-pro": "unidock-pro", "rosetta": "rosettaligand",
    }
    example_names = {
        "vina": "autodock-vina-example", "autodock-gpu": "autodock-gpu-example", "gnina": "gnina-example",
        "unidock-pro": "unidock-pro-example", "rosetta": "rosettaligand-example",
    }
    job_name_key = f"ligand-{tool_key}-job-name"
    default_job_name = defaults[tool_key]
    display_names = {
        "vina": "AutoDock Vina", "autodock-gpu": "AutoDock-GPU", "gnina": "GNINA",
        "unidock-pro": "Uni-Dock Pro", "rosetta": "RosettaLigand",
    }
    page_title = title or display_names[tool_key]
    paper_links = {
        "vina": ("https://github.com/ccsb-scripps/AutoDock-Vina", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3041641/"),
        "autodock-gpu": ("https://github.com/ccsb-scripps/AutoDock-GPU", "https://doi.org/10.1021/acs.jctc.0c01006"),
        "gnina": ("https://github.com/gnina/gnina", "https://doi.org/10.1186/s13321-025-00973-x"),
        "unidock-pro": ("https://github.com/dptech-corp/Uni-Dock", "https://doi.org/10.1093/bioinformatics/btac490"),
        "rosetta": ("https://github.com/RosettaCommons/rosetta", "https://docs.rosettacommons.org/docs/latest/rosettaligand-cleaned-dock"),
    }
    info_link_label = "Protocol docs ↗" if is_rosetta else "Paper ↗"
    source_url, paper_url = paper_links[tool_key]
    st.markdown("<style>" + (vina_inputs.EXAMPLE_ROOT / "form.css").read_text() + "</style>", unsafe_allow_html=True)
    with st.container(key="vina-page"):
        st.markdown(
            f'<div class="vina-heading"><h1>{html.escape(page_title)}</h1><p>Protein-ligand docking</p><nav>'
            f'<a href="{html.escape(source_url)}" target="_blank">Source code ↗</a>'
            f'<a href="{html.escape(paper_url)}" target="_blank">{info_link_label}</a>'
            '</nav></div>', unsafe_allow_html=True)
        if is_rosetta:
            st.caption("RosettaLigand uses Studio's bundled OpenVS runner.")
            if not rosetta_readiness["ready"]:
                st.warning(f"The bundled RosettaLigand runner is unavailable at {rosetta_readiness['source']}.")
        name_col, project_col, _ = st.columns([2.4, 1.9, 4])
        with project_col:
            project_id = project_selector(f"ligand-tool-{workflow_id}")
        _render_preparation_picker(owner_id, project_id, workflow_id)
        for chain, text in st.session_state.pop("vina-pending-residue-text", {}).items():
            st.session_state[f"vina-residue-text-{chain}"] = text
        job_name = name_col.text_input("Job Name", value=default_job_name, key=job_name_key).strip() or default_job_name
        source_path = _render_target_input(owner_id)
        target_path = None
        target_text = ""
        if source_path:
            try:
                target_path = _vina_target_as_pdb(owner_id, TargetSource(path=source_path))
                target_text = target_path.read_text() if target_path else ""
            except (ValueError, OSError) as exc:
                st.error(f"Could not read protein input: {exc}")
        if target_path is not None and target_text:
            from mn_molecular_design_studio.ligand_tools_ui import render_receptor_check

            render_receptor_check(
                target_text,
                target_path,
                owner_id=owner_id,
                key=f"ligand-{tool_key}-repair-receptor",
            )
        if target_path is not None:
            st.session_state["vina-normalized-target-path"] = str(target_path)
        else:
            st.session_state.pop("vina-normalized-target-path", None)
        coordinates = _vina_protein_atom_coordinates(target_text)
        sequences = vina_inputs.protein_sequence_entries(target_text)
        signature = hashlib.sha1(target_text.encode()).hexdigest()[:12]
        if target_text and st.session_state.get("vina-default-target-signature") != signature:
            st.session_state["vina-default-target-signature"] = signature
            st.session_state["ligand-vina-selected-residues"] = []
            st.session_state["ligand-vina-applied-residues"] = []
            st.session_state.pop("vina-applied-residue-padding", None)
            st.session_state["vina-derived-box-signature"] = None
            st.session_state.pop("vina-pending-prediction", None)
            st.session_state.pop("vina-applied-prediction", None)
            st.session_state.pop("vina-prediction-overlay-pdb", None)
            for index in range(10):
                st.session_state.pop(f"vina-probe-{index}", None)
            for key in ("vina-probe-rows", "vina-probe-values", "vina-probe-dismissed",
                        "vina-probe-initialized", "vina-probe-next-id"):
                st.session_state.pop(key, None)
            st.session_state["vina-boltz-selected-chains"] = list(sequences)
            st.session_state["vina-boltz-additional-sequences"] = ""
            st.session_state.pop("vina-residue-error", None)
            st.session_state.pop("vina-residue-apply-error", None)
            for chain in sequences:
                st.session_state[f"vina-residue-text-{chain}"] = ""
            default_ligand, ligand_center, ligand_size = _default_co_crystallized_ligand(target_path)
            st.session_state["ligand-vina-pocket-mode"] = "Co-crystallized Ligand" if default_ligand else "Manual"
            st.session_state["ligand-vina-box-ligand"] = default_ligand["key"] if default_ligand else ""
            if default_ligand:
                center, size = ligand_center, ligand_size
            else:
                points = [point for atoms in coordinates.values() for point in atoms]
                center = {
                    axis: round((min(point[i] for point in points) + max(point[i] for point in points)) / 2, 3)
                    if points else 0.0
                    for i, axis in enumerate("xyz")
                }
                size = dict.fromkeys("xyz", 10.0)
            for axis in "xyz":
                st.session_state[f"vina-box-center-{axis}"] = center[axis]
                st.session_state[f"vina-box-size-{axis}"] = size[axis]
        st.markdown('<div class="vina-box-title">Set Docking Box</div>', unsafe_allow_html=True)
        mode = "Manual"
        pocket_ready = True
        padding = 10.0
        selected_ligand = None
        if target_text:
            with st.container(horizontal=True, horizontal_alignment="center"):
                with st.container(key="vina-viewer", width=400):
                    overlay_pdb = st.session_state.get("vina-prediction-overlay-pdb", "")
                    viewer_structures = [StructureVisualization(
                        pdb=target_text, color="uniform", color_params={"value": "318c70"}, representation_type="cartoon",
                        highlighted_selections=[f"{c}:{n}" for c, n in st.session_state.get("ligand-vina-selected-residues", [])],
                        docking_box={"center": [float(st.session_state.get(f"vina-box-center-{a}", 0)) for a in "xyz"],
                                     "size": [float(st.session_state.get(f"vina-box-size-{a}", 10)) for a in "xyz"]},
                    )]
                    if overlay_pdb:
                        viewer_structures.append(StructureVisualization(
                            pdb=overlay_pdb, color="uniform", color_params={"value": "d84a7b"},
                            representation_type="ball-and-stick",
                        ))
                    overlay_signature = hashlib.sha1(overlay_pdb.encode("utf-8")).hexdigest()[:8]
                    value = molstar_custom_component(viewer_structures,
                        key=f"vina-molview-{signature}-{overlay_signature}-{st.session_state.get('vina-viewer-focus-epoch', 0)}",
                        height=400, width="100%", show_controls=False, show_structure_tools=False, selection_mode=True, minimal_viewer=True)
            if value is not None:
                event_signature = json.dumps(value, sort_keys=True)
                if st.session_state.get("vina-viewer-event") != event_signature:
                    st.session_state["vina-viewer-event"] = event_signature
                    picked = _selected_vina_protein_residues(value, coordinates)
                    for chain in sequences:
                        _set_residues(chain, {n for c, n in picked if c == chain})
            with st.container(key="vina-pocket"):
                st.markdown("**Define Binding Pocket**")
                tabs = st.tabs(["Co-crystallized Ligand", "Select Residues", "Predict (Boltz)", "Manual"], key="ligand-vina-pocket-mode", on_change="rerun")
                mode = st.session_state.get("ligand-vina-pocket-mode", "Manual")
                with tabs[0]:
                    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
                    choices = ligand_tools.bound_ligand_choices(target_path)
                    selected_key = st.selectbox(
                        "Bound ligand and chain for box",
                        ["", *[row["key"] for row in choices]],
                        format_func=lambda key: "Select a ligand…" if not key else next(f"{row['resname']} (Chain {row['chain']}, Res {row['resseq']}) — {row['heavy_atom_count']} atoms" for row in choices if row["key"] == key),
                        key="ligand-vina-box-ligand", on_change=_select_vina_bound_ligand,
                    )
                    st.caption("Choose the bound ligand copy that defines the docking box. The protein chains remain in the receptor structure.")
                    selected_ligand = next((row for row in choices if row["key"] == selected_key), None)
                    ligand_padding = st.number_input("Box padding (Å)", value=5.0, min_value=0.0, max_value=30.0, step=1.0, key="ligand-vina-ligand-padding", width=230)
                    if not choices:
                        st.caption("No co-crystallized ligands were found in this protein file.")
                    if selected_ligand:
                        st.caption("Docking box coordinates updated based on the selected ligand position.")
                with tabs[1]:
                    residue_padding = _render_residue_picker(sequences) if sequences else 10.0
                with tabs[2]:
                    prediction_ready = _render_boltz_pocket(
                        target_text, project_id, job_name, submit,
                        docking_engine=page_title,
                        target_path=target_path, compound_input=_active_compound_input(compound_sets),
                        placement_options=placement_options,
                        prediction_state=prediction_state,
                        active=mode == "Predict (Boltz)",
                    )
                with tabs[3]:
                    st.caption("Enter docking box coordinates directly in the fields below.")
            if mode == "Co-crystallized Ligand":
                pocket_ready = selected_ligand is not None
            elif mode == "Select Residues":
                selected = sorted(tuple(row) for row in st.session_state.get("ligand-vina-selected-residues", []))
                applied = sorted(tuple(row) for row in st.session_state.get("ligand-vina-applied-residues", []))
                applied_padding = st.session_state.get("vina-applied-residue-padding")
                pocket_ready = (
                    bool(selected) and selected == applied
                    and applied_padding is not None and float(applied_padding) == float(residue_padding)
                    and not st.session_state.get("vina-residue-error")
                )
            elif mode == "Predict (Boltz)":
                pocket_ready = prediction_ready
            derived = None
            derived_signature = ""
            if mode == "Co-crystallized Ligand" and selected_ligand:
                try:
                    derived = _vina_box_from_points(list(selected_ligand["heavy_atom_coordinates"].values()), float(ligand_padding))
                    derived_signature = json.dumps([signature, mode, selected_ligand["key"], ligand_padding])
                except ValueError as exc:
                    st.error(str(exc))
            elif mode in {"Select Residues", "Predict (Boltz)"}:
                if mode == "Select Residues":
                    padding = float(st.session_state.get("vina-applied-residue-padding", residue_padding))
                    selected = st.session_state.get("ligand-vina-applied-residues", [])
                else:
                    padding = 10.0
                    selected = st.session_state.get("ligand-vina-selected-residues", [])
                points = [point for c, n in selected for point in coordinates.get((c, n), [])]
                if points:
                    try:
                        derived = _vina_box_from_points(points, padding)
                        derived_signature = json.dumps([signature, mode, selected, padding])
                    except ValueError as exc:
                        st.error(str(exc))
            if derived and _set_vina_derived_box(derived_signature, *derived):
                st.rerun()
            if mode == "Manual":
                st.session_state["vina-derived-box-signature"] = None
        with st.container(horizontal=True, horizontal_alignment="center"):
            with st.container(key="vina-box-fields", width=850):
                cols = st.columns(3, gap="large")
                center = {f"center_{axis}": float(cols[i].number_input(f"box{axis.upper()}", value=0.0, step=.1, format="%.3f", key=f"vina-box-center-{axis}", help=f"{axis.upper()} coordinate of bounding box center")) for i, axis in enumerate("xyz")}
                cols = st.columns(3, gap="large")
                size = {f"size_{axis}": float(cols[i].number_input(label, value=10.0, min_value=1.0, max_value=60.0, step=1.0, key=f"vina-box-size-{axis}", help=f"{label.capitalize()} of bounding box in Å")) for i, (axis, label) in enumerate(zip("xyz", ("width", "height", "depth")))}
        if is_rosetta:
            st.caption("RosettaLigand anchors the search at the displayed box center. Its protocol-specific search padding is set below.")
        compounds = _render_ligand_input(owner_id, project_id, compound_sets)
        with st.container(key="vina-optional"):
            with st.expander("View Optional Settings"):
                alignment_upload = None
                if is_rosetta:
                    reference_mode = st.selectbox(
                        "Placement mode", ["pocket_center", "reference_guided"],
                        format_func=lambda value: "Pocket-centered global search" if value == "pocket_center" else "Reference-guided search",
                        key="ligand-rosetta-reference-mode",
                    )
                    protocol = st.selectbox(
                        "Protocol", ["vsh", "vsx", "convergence"],
                        format_func=lambda value: {"vsh": "VSH", "vsx": "VSX", "convergence": "Convergence"}[value],
                        key="ligand-rosetta-protocol",
                    )
                    reference_upload = st.file_uploader(
                        "Reference ligand coordinates (SDF) *" if reference_mode == "reference_guided" else "Reference ligand coordinates (SDF, optional)",
                        type=["sdf"], key="ligand-rosetta-reference-ligand",
                    )
                    cols = st.columns(3)
                    limit = int(cols[0].number_input("Compound limit (0 = all selected)", value=0, min_value=0, max_value=5000, step=1, key="ligand-rosetta-limit"))
                    cores = int(cols[1].number_input("CPU workers", value=8, min_value=1, max_value=128, key="ligand-rosetta-cpu"))
                    replicates = int(cols[2].number_input("Replicates", value=1, min_value=1, max_value=100, key="ligand-rosetta-replicates"))
                    cols = st.columns(3)
                    ph = float(cols[0].number_input("Ligand protonation pH", value=7.4, min_value=0.0, max_value=14.0, step=0.1, key="ligand-rosetta-ph"))
                    conformers = int(cols[1].number_input("Input conformers", value=20, min_value=1, max_value=1000, key="ligand-rosetta-conformers"))
                    minimization_steps = int(cols[2].number_input("Minimization steps", value=2000, min_value=1, max_value=100000, step=100, key="ligand-rosetta-min-steps"))
                    cols = st.columns(2)
                    padding = float(cols[0].number_input("Search padding (Å)", value=4.0, min_value=1.0, max_value=30.0, step=0.5, key="ligand-rosetta-padding"))
                    seed_start = int(cols[1].number_input("Seed start", value=1001, min_value=1, max_value=2_147_483_646, key="ligand-rosetta-seed-start"))
                    st.caption("RosettaLigand scores are reported in REU. Compare scores within the same protocol.")
                    parameters = {
                        **{f"center_{axis}": center[f"center_{axis}"] for axis in "xyz"},
                        "reference_mode": reference_mode, "protocol": protocol,
                        "compound_limit": limit, "cpu_cores": cores,
                        "replicates": replicates, "seed_start": seed_start,
                        "ph": ph, "conformers": conformers,
                        "minimization_steps": minimization_steps, "padding": padding,
                        "job_name": job_name,
                    }
                else:
                    st.markdown("**Ligand preparation**")
                    prep_cols = st.columns(3)
                    use_scrub = prep_cols[0].checkbox(
                        "Use Scrub",
                        value=True,
                        key=f"ligand-{tool_key}-use-scrub",
                        help="Normalize ligand protonation and chemistry with Scrub before Meeko prepares the ligand for docking.",
                    )
                    scrub_ph = float(prep_cols[1].number_input(
                        "Scrub pH",
                        min_value=0.0,
                        max_value=14.0,
                        value=7.4,
                        step=0.1,
                        disabled=not use_scrub,
                        key=f"ligand-{tool_key}-scrub-ph",
                    ))
                    scrub_skip_tautomer = prep_cols[2].checkbox(
                        "Skip tautomer enumeration",
                        value=True,
                        disabled=not use_scrub,
                        key=f"ligand-{tool_key}-skip-tautomer",
                        help="When cleared, Scrub may generate tautomer states before docking.",
                    )
                    st.markdown("**Search settings**")
                    cols = st.columns(4)
                    if is_autodock_gpu:
                        nrun = int(cols[0].number_input(
                            "LGA runs per compound", value=20, min_value=1, max_value=500,
                            key="ligand-autodock-gpu-nrun",
                            help="Number of AutoDock4 genetic-algorithm runs per ligand. AutoDock-GPU returns its best pose.",
                        ))
                        limit = int(cols[1].number_input(
                            "Compound limit (0 = all selected)", value=0, min_value=0, max_value=5000,
                            key="ligand-autodock-gpu-compound-limit",
                        ))
                        cores = int(cols[2].number_input(
                            "CPU cores", value=4, min_value=1, max_value=256,
                            key="ligand-autodock-gpu-cpu-cores",
                        ))
                        gpu_device = cols[3].text_input(
                            "GPU device", value="0", key="ligand-autodock-gpu-gpu",
                            help="Use a device ID such as 0.",
                        )
                    elif is_unidock:
                        search_mode = cols[0].selectbox("Search effort", ["fast", "balance", "detail"], index=1, key="ligand-unidock-pro-search-mode")
                    else:
                        exhaustiveness = int(cols[0].number_input(
                            "Exhaustiveness", value=32 if is_vina else 8, min_value=1, max_value=1024,
                            key=f"ligand-{tool_key}-exhaustiveness",
                            help=(
                                "32 is a thorough starting point for pose-focused docking; expect longer runtimes."
                                if is_vina else "Increasing this to 32 spends more search effort and can give more consistent results."
                            ),
                        ))
                    if not is_autodock_gpu:
                        modes = int(cols[1].number_input(
                            "Poses per compound", value=10 if is_vina else 9, min_value=1, max_value=100,
                            key=f"ligand-{tool_key}-num-modes",
                            help="Maximum binding modes to return. This does not increase search effort, and the engine may return fewer.",
                        ))
                        limit = int(cols[2].number_input("Compound limit (0 = all selected)", value=0, min_value=0,
                                                         max_value=compound_limit_max, key=f"ligand-{tool_key}-compound-limit"))
                        cores = int(cols[3].number_input("CPU cores", value=4, min_value=1, max_value=256, key=f"ligand-{tool_key}-cpu-cores"))
                        gpu_device = st.text_input("GPU device", value="0", key=f"ligand-{tool_key}-gpu", help="Use a device ID such as 0.") if (is_gnina or is_unidock) else "none"
                    seed = int(st.number_input("Random seed", value=42, min_value=1, max_value=2_147_483_646, key=f"ligand-{tool_key}-seed"))
                    cnn_scoring = "rescore"
                    if is_gnina:
                        cnn_scoring = st.selectbox(
                            "CNN scoring", ["rescore", "refinement", "all", "none"],
                            format_func=lambda value: {
                                "rescore": "Rescore poses (recommended)",
                                "refinement": "Refine and rescore",
                                "all": "CNN throughout docking",
                                "none": "Empirical scoring only",
                            }[value],
                            key="ligand-gnina-cnn-scoring",
                            help="Refinement uses more GPU time; CNN scoring throughout docking is the slowest option.",
                        )
                        compound_workers = int(st.number_input(
                            "Concurrent compounds", min_value=1, max_value=max(1, min(4, cores, limit)),
                            value=max(1, min(2, cores, limit)), step=1,
                            key="ligand-gnina-compound-workers",
                            help="Run independent compounds concurrently on the assigned GPU. CPU threads are divided across workers; use 1 if GPU memory is tight.",
                        ))
                    align = st.checkbox("Align a reference ligand's longest axis with global X", key=f"ligand-{tool_key}-align-axis-x")
                    alignment_upload = st.file_uploader("Coordinate-bearing alignment ligand (SDF)", type=["sdf"], key=f"ligand-{tool_key}-alignment-ligand") if align else None
                    if is_gnina:
                        st.caption("CNN scoring changes both runtime and result ranking. GNINA can process independent compounds concurrently in one job; each worker shares the assigned GPU and receives a portion of the CPU allocation.")
                    if is_unidock:
                        st.caption("Uni-Dock Pro requires a CUDA GPU. Search effort controls the speed and sampling depth.")
                    if is_autodock_gpu:
                        st.caption("AutoDock-GPU uses AutoGrid maps and AutoDock4 scoring. The selected GPU runs compounds one at a time; LGA runs control per-compound search effort.")
                    st.caption("Alignment rotates the receptor, box, and returned poses together. Scores and poses remain available to downstream pipelines.")
                    parameters = {
                        **center, **size,
                        **({"nrun": nrun} if is_autodock_gpu else {"search_mode": search_mode} if is_unidock else {"exhaustiveness": exhaustiveness}),
                        **({} if is_autodock_gpu else {"num_modes": modes}),
                        "compound_limit": limit, "cpu_cores": cores,
                        "seed": seed, "gpu_device": gpu_device,
                        "align_longest_axis_to_x": align,
                        "use_scrub": use_scrub,
                        "scrub_ph": scrub_ph,
                        "scrub_skip_tautomer": scrub_skip_tautomer,
                        "job_name": job_name,
                    }
                    if is_gnina:
                        parameters["cnn_scoring"] = cnn_scoring
                        parameters["compound_workers"] = compound_workers
        with st.container(key="vina-example"):
            st.button("Load example inputs", on_click=_load_example,
                kwargs={"owner_id": owner_id, "project_id": str(project_id or ""), "job_name_key": job_name_key, "example_job_name": example_names[tool_key]},
                key="vina-load-example")
        runtime_estimate = (
            "Depends on ligand count, box size, LGA runs per compound, and GPU availability."
            if is_autodock_gpu else
            "Depends on ligand count, concurrent compound workers, CNN mode, and GPU availability."
            if is_gnina else
            "Depends on ligand count, box size, CPU allocation, and exhaustiveness."
        )
        st.markdown(
            f'<div class="vina-runtime"><strong>Estimated Runtime</strong>{html.escape(runtime_estimate)}</div>',
            unsafe_allow_html=True,
        )
        with st.container(key="vina-submit"):
            rosetta_ready = not is_rosetta or rosetta_readiness["ready"]
            reference_required = is_rosetta and reference_mode == "reference_guided"
            run = st.button("Submit", type="primary", disabled=(not project_id or not target_path or not coordinates or not compounds or not pocket_ready or not rosetta_ready or (reference_required and reference_upload is None) or (not is_rosetta and align and alignment_upload is None)), key=f"ligand-{tool_key}-run", use_container_width=True)
            st.caption(f"Submit 1 job · {min(compounds['compound_count'], limit) if compounds else 0} ligand(s)")
        if run:
            try:
                parameters.update({"compound_set_file": str(compounds["path"]), "target_pdb": str(target_path)})
                if is_rosetta and reference_upload is not None:
                    path = _save_reference_ligand_upload(owner_id, reference_upload.name, reference_upload.getvalue())
                    parameters["reference_ligand_sdf"] = str(path)
                if not is_rosetta and alignment_upload is not None:
                    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
                    path = _save_reference_ligand_upload(owner_id, alignment_upload.name, alignment_upload.getvalue())
                    ligand_tools._sdf_heavy_atom_coordinates(path)
                    parameters["alignment_ligand_sdf"] = str(path)
                submit(schema, workflow_id, parameters, str(project_id), {}, title)
            except (ValueError, OSError) as exc:
                st.error(f"Could not start {page_title}: {exc}")
        try:
            _save_current_preparation()
        except (ValueError, OSError) as exc:
            st.error(f"Could not save the docking preparation: {exc}")


# Compatibility for existing tests and callers; docking engines share this page implementation.
render_vina_setup = render_docking_setup
