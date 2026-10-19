"""AutoDock Vina's protein, pocket, and ligand input page."""

from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import streamlit as st

from mn_molecular_design_studio import vina_inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.hotspot_sequence import render_hotspot_sequence
from mn_molecular_design_studio.target_editor import TargetSource, _target_storage_path
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
    StructureVisualization, molstar_custom_component,
)


def _load_example(owner_id: str, project_id: str = "", *, job_name_key: str = "ligand-vina-job-name", example_job_name: str = "autodock-vina-example") -> None:
    target, compounds = vina_inputs.stage_example(owner_id, project_id=project_id)
    st.session_state.update({
        "vina-target-path": str(target), "vina-compound-input": compounds,
        "vina-ligand-format": "SDF File", "vina-protein-format": "PDB",
        job_name_key: example_job_name,
        "ligand-vina-pocket-mode": "Manual", "ligand-vina-selected-residues": [],
        "vina-default-target-signature": hashlib.sha1(target.read_bytes()).hexdigest()[:12],
        "vina-derived-box-signature": None,
    })
    for axis, value in zip("xyz", (15.190, 53.903, 16.917)):
        st.session_state[f"vina-box-center-{axis}"] = value
        st.session_state[f"vina-box-size-{axis}"] = 20.0
    st.session_state["vina-input-epoch"] = st.session_state.get("vina-input-epoch", 0) + 1
    st.session_state.pop("vina-pending-prediction", None)
    st.session_state.pop("vina-applied-prediction", None)
    st.session_state.pop("vina-residue-error", None)


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
    padding = float(padding_col.number_input("Box padding (Å)", value=10.0, min_value=0.0, max_value=30.0, step=1.0, key="ligand-vina-residue-padding"))
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
    return padding


def _render_boltz_pocket(target_text: str, project_id: str | None, job_name: str, submit) -> bool:
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

    st.caption("Enter 1–10 probe SMILES. Boltz will predict where they bind, and the docking box will be set around the predicted binding site.")
    probe_count = st.session_state.get("vina-probe-count", 1)
    probes = [st.text_input(f"Probe SMILES {index + 1}", key=f"vina-probe-{index}") for index in range(probe_count)]
    if st.button("+ Add Probe", disabled=probe_count >= 10, key="vina-add-probe"):
        st.session_state["vina-probe-count"] = probe_count + 1
        st.rerun()
    probes = [value.strip() for value in probes if value.strip()]
    st.caption(f"{len(probes)}/10 probes")
    distance = float(st.number_input("Contact distance (Å)", value=5.0, min_value=1.0, max_value=20.0, step=1.0, key="vina-probe-distance"))
    if st.button("Predict Binding Site", type="primary", disabled=not project_id or not probes, key="vina-predict-pocket"):
        try:
            from rdkit import Chem
            if any(Chem.MolFromSmiles(value) is None for value in probes):
                raise ValueError("Every probe must be a valid SMILES string.")
            proteins = "\n".join("".join(letter for _, letter in sequence) for sequence in vina_inputs.protein_sequence_entries(target_text).values())
            parameters = {"engine": "boltz2", "job_name": f"{job_name}-pocket", "samples": 1,
                "tasks": [{"name": f"probe_{i+1}", "protein_sequences": proteins, "ligand_smiles": value} for i, value in enumerate(probes)],
                "cpu_cores": 4, "gpu_device": "0", "engine_options": {
                    "use_msa": False, "use_msa_server": False, "generate_ipsae": False, "output_format": "pdb"}}
            previous = st.session_state.get("last_submitted_job")
            submit(workflow_schema("refolding.optimized_prediction"), "refolding.optimized_prediction", parameters, str(project_id), {}, "Boltz binding-site prediction")
            job = st.session_state.get("last_submitted_job")
            if isinstance(job, dict) and job != previous and job.get("run_dir"):
                st.session_state["vina-pending-prediction"] = {"run_dir": job["run_dir"], "target_hash": hashlib.sha256(target_text.encode()).hexdigest(), "distance": distance}
        except (ValueError, OSError) as exc:
            st.error(f"Could not start pocket prediction: {exc}")
    pending = st.session_state.get("vina-pending-prediction")
    if pending and pending.get("target_hash") == hashlib.sha256(target_text.encode()).hexdigest():
        run_dir = Path(pending["run_dir"])
        metadata = read_json(run_dir / "metadata.json") if (run_dir / "metadata.json").is_file() else {}
        st.caption(f"Boltz prediction: {metadata.get('status', 'queued')}")
        st.button("Refresh prediction", key="vina-refresh-pocket")
        candidates = read_candidates(run_dir)
        application_signature = json.dumps([str(run_dir), distance])
        if candidates and st.session_state.get("vina-applied-prediction") != application_signature:
            try:
                residues = set()
                for candidate in candidates:
                    structure = Path(candidate.get("complex_pdb") or "").resolve()
                    if structure.is_file() and structure.is_relative_to(run_dir.resolve()):
                        residues |= vina_inputs.predicted_contact_residues(target_text, structure.read_text(), distance)
                if not residues:
                    raise ValueError("No protein residues contact the predicted probes at this distance.")
                for chain in vina_inputs.protein_sequence_entries(target_text):
                    _set_residues(chain, {n for c, n in residues if c == chain})
                st.session_state["vina-applied-prediction"] = application_signature
                st.rerun()
            except (ValueError, OSError) as exc:
                st.error(f"Could not apply the predicted site: {exc}")
        return st.session_state.get("vina-applied-prediction") == application_signature
    return False


def _use_database_set(group_id: str) -> None:
    st.session_state["vina-ligand-format"] = "Database compound set"
    st.session_state["ligand-tool-compound-set-ligand.docking.vina"] = group_id


def _remove_ligand_input() -> None:
    st.session_state.pop("vina-compound-input", None)
    st.session_state["vina-input-epoch"] = st.session_state.get("vina-input-epoch", 0) + 1


def _render_ligand_input(owner_id: str, project_id: str | None, compound_sets: list[dict]) -> dict | None:
    fmt = st.selectbox("Ligand Format", ["SDF File", "SMILES (CSV)", "Database compound set"], key="vina-ligand-format", width=300,
        help="Upload SDF coordinates or a CSV containing a SMILES column, or choose a compound set from Studio Database.")
    if fmt == "Database compound set":
        if not compound_sets:
            st.caption("Import compounds in Database, or choose SDF File or SMILES (CSV) above.")
            return None
        group_id = st.selectbox("My compounds", [str(row["group_id"]) for row in compound_sets],
            format_func=lambda value: next(f"{row['name']} · {row['compound_count']:,} compounds" for row in compound_sets if str(row["group_id"]) == value),
            key="ligand-tool-compound-set-ligand.docking.vina")
        return next(row for row in compound_sets if str(row["group_id"]) == group_id)
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
                return selected
    return None


def render_docking_setup(workflow_id, schema, *, owner_id, compound_sets, project_selector, submit, title):
    from mn_molecular_design_studio.ligand_tools_ui import (
        _vina_target_as_pdb, _vina_protein_atom_coordinates, _selected_vina_protein_residues,
        _vina_box_from_points, _set_vina_derived_box, _save_reference_ligand_upload,
    )
    tool_key = {
        "ligand.docking.vina": "vina",
        "ligand.docking.gnina": "gnina",
        "ligand.docking.unidock_pro": "unidock-pro",
        "ligand.docking.rosetta_ligand": "rosetta",
    }.get(workflow_id)
    if tool_key is None:
        raise ValueError(f"The shared docking page does not support {workflow_id}.")
    is_gnina = workflow_id == "ligand.docking.gnina"
    is_unidock = workflow_id == "ligand.docking.unidock_pro"
    is_rosetta = workflow_id == "ligand.docking.rosetta_ligand"
    rosetta_readiness = None
    if is_rosetta:
        from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools

        rosetta_readiness = ligand_tools.rosettaligand_runner_readiness()
    defaults = {
        "vina": "autodock-vina", "gnina": "gnina",
        "unidock-pro": "unidock-pro", "rosetta": "rosettaligand",
    }
    example_names = {
        "vina": "autodock-vina-example", "gnina": "gnina-example",
        "unidock-pro": "unidock-pro-example", "rosetta": "rosettaligand-example",
    }
    job_name_key = f"ligand-{tool_key}-job-name"
    default_job_name = defaults[tool_key]
    display_names = {
        "vina": "AutoDock Vina", "gnina": "GNINA",
        "unidock-pro": "Uni-Dock Pro", "rosetta": "RosettaLigand",
    }
    page_title = title or display_names[tool_key]
    paper_links = {
        "vina": ("https://github.com/ccsb-scripps/AutoDock-Vina", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3041641/"),
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
        job_name = name_col.text_input("Job Name", value=default_job_name, key=job_name_key).strip() or default_job_name
        with project_col:
            project_id = project_selector(f"ligand-tool-{workflow_id}")
        source_path = _render_target_input(owner_id)
        target_path = None
        target_text = ""
        if source_path:
            try:
                target_path = _vina_target_as_pdb(owner_id, TargetSource(path=source_path))
                target_text = target_path.read_text() if target_path else ""
            except (ValueError, OSError) as exc:
                st.error(f"Could not read protein input: {exc}")
        coordinates = _vina_protein_atom_coordinates(target_text)
        sequences = vina_inputs.protein_sequence_entries(target_text)
        signature = hashlib.sha1(target_text.encode()).hexdigest()[:12]
        if target_text and st.session_state.get("vina-default-target-signature") != signature:
            st.session_state["vina-default-target-signature"] = signature
            st.session_state["ligand-vina-selected-residues"] = []
            st.session_state["vina-derived-box-signature"] = None
            st.session_state.pop("vina-pending-prediction", None)
            st.session_state.pop("vina-applied-prediction", None)
            st.session_state.pop("vina-residue-error", None)
            for chain in sequences:
                st.session_state[f"vina-residue-text-{chain}"] = ""
            points = [point for atoms in coordinates.values() for point in atoms]
            for i, axis in enumerate("xyz"):
                st.session_state[f"vina-box-center-{axis}"] = round((min(p[i] for p in points) + max(p[i] for p in points)) / 2, 3) if points else 0.0
                st.session_state[f"vina-box-size-{axis}"] = 10.0
        st.markdown('<div class="vina-box-title">Set Docking Box</div>', unsafe_allow_html=True)
        mode = "Manual"
        pocket_ready = True
        padding = 10.0
        selected_ligand = None
        if target_text:
            with st.container(horizontal=True, horizontal_alignment="center"):
                with st.container(key="vina-viewer", width=400):
                    value = molstar_custom_component([StructureVisualization(
                        pdb=target_text, color="uniform", color_params={"value": "318c70"}, representation_type="cartoon",
                        highlighted_selections=[f"{c}:{n}" for c, n in st.session_state.get("ligand-vina-selected-residues", [])],
                        docking_box={"center": [float(st.session_state.get(f"vina-box-center-{a}", 0)) for a in "xyz"],
                                     "size": [float(st.session_state.get(f"vina-box-size-{a}", 10)) for a in "xyz"]},
                    )], key=f"vina-molview-{signature}", height=400, width="100%", show_controls=False, show_structure_tools=False, selection_mode=True, minimal_viewer=True)
            if value is not None:
                event_signature = json.dumps(value, sort_keys=True)
                if st.session_state.get("vina-viewer-event") != event_signature:
                    st.session_state["vina-viewer-event"] = event_signature
                    picked = _selected_vina_protein_residues(value, coordinates)
                    for chain in sequences:
                        _set_residues(chain, {n for c, n in picked if c == chain})
            with st.container(key="vina-pocket"):
                st.markdown("**Define Binding Pocket**")
                tabs = st.tabs(["Manual", "Co-crystallized Ligand", "Select Residues", "Predict (Boltz)"], key="ligand-vina-pocket-mode", on_change="rerun")
                mode = st.session_state.get("ligand-vina-pocket-mode", "Manual")
                with tabs[0]:
                    st.caption("Enter docking box coordinates directly in the fields below.")
                with tabs[1]:
                    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
                    choices = ligand_tools.bound_ligand_choices(target_path)
                    selected_key = st.selectbox("Co-crystallized ligand", ["", *[row["key"] for row in choices]], label_visibility="collapsed",
                        format_func=lambda key: "Select a ligand…" if not key else next(f"{row['resname']} (Chain {row['chain']}, Res {row['resseq']}) — {row['heavy_atom_count']} atoms" for row in choices if row["key"] == key), key="ligand-vina-box-ligand")
                    selected_ligand = next((row for row in choices if row["key"] == selected_key), None)
                    ligand_padding = st.number_input("Box padding (Å)", value=5.0, min_value=0.0, max_value=30.0, step=1.0, key="ligand-vina-ligand-padding", width=230)
                    if not choices:
                        st.caption("No co-crystallized ligands were found in this protein file.")
                    if selected_ligand:
                        st.caption("Docking box coordinates updated based on the selected ligand position.")
                with tabs[2]:
                    residue_padding = _render_residue_picker(sequences) if sequences else 10.0
                with tabs[3]:
                    prediction_ready = _render_boltz_pocket(target_text, project_id, job_name, submit)
            if mode == "Co-crystallized Ligand":
                pocket_ready = selected_ligand is not None
            elif mode == "Select Residues":
                pocket_ready = bool(st.session_state.get("ligand-vina-selected-residues")) and not st.session_state.get("vina-residue-error")
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
                padding = float(residue_padding) if mode == "Select Residues" else 10.0
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
                    limit = int(cols[0].number_input("Maximum compounds", value=10, min_value=1, max_value=5000, step=1, key="ligand-rosetta-limit"))
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
                    cols = st.columns(4)
                    if is_unidock:
                        search_mode = cols[0].selectbox("Search effort", ["fast", "balance", "detail"], index=1, key="ligand-unidock-pro-search-mode")
                    else:
                        exhaustiveness = int(cols[0].number_input("Exhaustiveness", value=8, min_value=1, max_value=1024, key=f"ligand-{tool_key}-exhaustiveness", help="Increasing this to 32 spends more search effort and can give more consistent results."))
                    modes = int(cols[1].number_input("Poses per compound", value=9, min_value=1, max_value=100, key=f"ligand-{tool_key}-num-modes"))
                    limit = int(cols[2].number_input("Maximum compounds", value=100, min_value=1, max_value=5000, key=f"ligand-{tool_key}-compound-limit"))
                    cores = int(cols[3].number_input("CPU cores", value=4, min_value=1, max_value=256, key=f"ligand-{tool_key}-cpu-cores"))
                    seed = int(st.number_input("Random seed", value=42, min_value=1, max_value=2_147_483_646, key=f"ligand-{tool_key}-seed"))
                    gpu_device = st.text_input("GPU device", value="0", key=f"ligand-{tool_key}-gpu", help="Use a device ID such as 0.") if (is_gnina or is_unidock) else "none"
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
                    align = st.checkbox("Align a reference ligand's longest axis with global X", key=f"ligand-{tool_key}-align-axis-x")
                    alignment_upload = st.file_uploader("Coordinate-bearing alignment ligand (SDF)", type=["sdf"], key=f"ligand-{tool_key}-alignment-ligand") if align else None
                    if is_gnina:
                        st.caption("CNN scoring changes both runtime and result ranking. With empirical scoring only, rank by docking energy.")
                    if is_unidock:
                        st.caption("Uni-Dock Pro requires a CUDA GPU. Search effort controls the speed and sampling depth.")
                    st.caption("Alignment rotates the receptor, box, and returned poses together. Scores and poses remain available to downstream pipelines.")
                    parameters = {
                        **center, **size,
                        **({"search_mode": search_mode} if is_unidock else {"exhaustiveness": exhaustiveness}),
                        "num_modes": modes, "compound_limit": limit, "cpu_cores": cores,
                        "seed": seed, "gpu_device": gpu_device,
                        "align_longest_axis_to_x": align,
                        "job_name": job_name,
                    }
                    if is_gnina:
                        parameters["cnn_scoring"] = cnn_scoring
        with st.container(key="vina-example"):
            st.button("Load example inputs", on_click=_load_example,
                kwargs={"owner_id": owner_id, "project_id": str(project_id or ""), "job_name_key": job_name_key, "example_job_name": example_names[tool_key]},
                key="vina-load-example")
        st.markdown('<div class="vina-runtime"><strong>Estimated Runtime</strong>Depends on ligand count, box size, CPU allocation, and exhaustiveness.</div>', unsafe_allow_html=True)
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


# Compatibility for existing tests and callers; all four engines share the same page implementation.
render_vina_setup = render_docking_setup
