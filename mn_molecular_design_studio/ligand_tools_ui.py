"""Dedicated setup cards for standalone ligand tool tiles."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Callable

import streamlit as st

from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.docking_structure_check import render_receptor_check
from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
from mn_molecular_design_studio.target_editor import TargetSource


def _save_target_upload(owner_id: str, name: str, content: bytes) -> Path:
    suffix = Path(name).suffix.lower()
    if suffix not in {".pdb", ".ent"}:
        raise ValueError("Upload a PDB or ENT protein target structure.")
    digest = hashlib.sha256(content).hexdigest()[:20]
    target = user_input_root(owner_id) / "ligand_tools" / "targets" / f"target-{digest}.pdb"
    target.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    if not target.exists():
        target.write_bytes(content)
        try:
            target.chmod(0o600)
        except OSError:
            pass
    ligand_tools.validate_target_pdb(target)
    return target


def _vina_target_as_pdb(owner_id: str, source: TargetSource) -> Path | None:
    if source.path is None:
        return None
    path = source.path.resolve()
    content = path.read_bytes()
    if path.suffix.lower() in {".cif", ".mmcif"}:
        from mn_molecular_design_studio.target_editor import _rcsb_mmcif_to_pdb

        content = _rcsb_mmcif_to_pdb(content, allow_non_protein=True)
        name = f"{path.stem}.pdb"
    else:
        name = path.name
    return _save_target_upload(owner_id, name, content)


def _vina_protein_atom_coordinates(pdb_text: str) -> dict[tuple[str, int], list[tuple[float, float, float]]]:
    coordinates: dict[tuple[str, int], list[tuple[float, float, float]]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else line[12:16].strip()[:1].upper()
        if element in {"H", "D"}:
            continue
        try:
            residue = int(line[22:26])
            point = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        chain = line[21].strip() or "_"
        coordinates.setdefault((chain, residue), []).append(point)
    return coordinates


def _vina_box_from_points(
    points: list[tuple[float, float, float]], padding: float = 5.0,
) -> tuple[dict[str, float], dict[str, float]]:
    if not points:
        raise ValueError("Select protein residues or a co-crystallized ligand with coordinates.")
    if not 0 <= float(padding) <= 30:
        raise ValueError("Pocket padding must be between 0 and 30 Å.")
    lower = [min(point[axis] for point in points) for axis in range(3)]
    upper = [max(point[axis] for point in points) for axis in range(3)]
    center = {axis: round((lower[index] + upper[index]) / 2, 3) for index, axis in enumerate("xyz")}
    size = {
        axis: round(upper[index] - lower[index] + 2 * float(padding), 3)
        for index, axis in enumerate("xyz")
    }
    if any(value < 1 for value in size.values()):
        size = {axis: max(1.0, value) for axis, value in size.items()}
    if any(value > 60 for value in size.values()):
        raise ValueError("The calculated pocket box exceeds 60 Å. Reduce the padding or select a smaller pocket.")
    return center, size


def _selected_vina_protein_residues(
    viewer_value: object,
    residue_coordinates: dict[tuple[str, int], list[tuple[float, float, float]]],
) -> set[tuple[str, int]]:
    if isinstance(viewer_value, str):
        try:
            viewer_value = json.loads(viewer_value)
        except json.JSONDecodeError:
            return set()
    if not isinstance(viewer_value, dict):
        return set()
    selected: set[tuple[str, int]] = set()
    for entry in viewer_value.get("sequenceSelections") or []:
        if not isinstance(entry, dict):
            continue
        chain = str(entry.get("chainId") or "").strip() or "_"
        for raw_residue in entry.get("residues") or []:
            try:
                residue = int(raw_residue)
            except (TypeError, ValueError):
                continue
            if (chain, residue) in residue_coordinates:
                selected.add((chain, residue))
    return selected


def _set_vina_derived_box(
    signature: str,
    center: dict[str, float],
    size: dict[str, float],
) -> bool:
    signature_key = "vina-derived-box-signature"
    if st.session_state.get(signature_key) == signature:
        return False
    for axis in "xyz":
        st.session_state[f"vina-box-center-{axis}"] = center[axis]
        st.session_state[f"vina-box-size-{axis}"] = size[axis]
    st.session_state[signature_key] = signature
    return True


def _render_vina_setup(
    workflow_id: str,
    schema: dict[str, Any],
    *,
    owner_id: str,
    compound_sets: list[dict[str, Any]],
    project_selector: Callable[..., str | None],
    submit: Callable[..., None],
    title: str,
    placement_options: Callable[..., dict[str, Any]] | None = None,
    prediction_state: Callable[..., dict[str, Any]] | None = None,
) -> None:
    from mn_molecular_design_studio.vina_ui import render_vina_setup

    render_vina_setup(
        workflow_id, schema, owner_id=owner_id, compound_sets=compound_sets,
        project_selector=project_selector, submit=submit, title=title,
        placement_options=placement_options,
        prediction_state=prediction_state,
    )


def _save_reference_ligand_upload(owner_id: str, name: str, content: bytes) -> Path:
    suffix = Path(name).suffix.lower()
    if suffix != ".sdf":
        raise ValueError("Upload the reference ligand as an SDF file.")
    digest = hashlib.sha256(content).hexdigest()[:20]
    ligand = user_input_root(owner_id) / "ligand_tools" / "reference_ligands" / f"reference-{digest}.sdf"
    ligand.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    if not ligand.exists():
        ligand.write_bytes(content)
        try:
            ligand.chmod(0o600)
        except OSError:
            pass
    ligand_tools.validate_ligand_sdf(ligand)
    return ligand


def _save_pgmg_pharmacophore_upload(owner_id: str, name: str, content: bytes) -> Path:
    if Path(name).suffix.lower() != ".posp":
        raise ValueError("Upload a PGMG .posp pharmacophore file.")
    digest = hashlib.sha256(content).hexdigest()[:20]
    path = user_input_root(owner_id) / "ligand_tools" / "pharmacophores" / f"pharmacophore-{digest}.posp"
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    if not path.exists():
        path.write_bytes(content)
        try:
            path.chmod(0o600)
        except OSError:
            pass
    ligand_tools._validate_pgmg_posp(path)
    return path


def _save_omtra_pharmacophore_upload(owner_id: str, name: str, content: bytes) -> Path:
    if Path(name).suffix.lower() != ".json":
        raise ValueError("Upload an OMTRA Pharmit pharmacophore as JSON.")
    digest = hashlib.sha256(content).hexdigest()[:20]
    path = user_input_root(owner_id) / "ligand_tools" / "pharmacophores" / f"pharmit-{digest}.json"
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    if not path.exists():
        path.write_bytes(content)
        try:
            path.chmod(0o600)
        except OSError:
            pass
    ligand_tools._validate_omtra_pharmacophore_json(path)
    return path


def _render_complex_preparation_setup(
    workflow_id: str,
    schema: dict[str, Any],
    *,
    owner_id: str,
    project_selector: Callable[..., str | None],
    submit: Callable[..., None],
    title: str,
) -> None:
    st.write(
        "Import a protein–ligand complex, run Studio's bundled protein-repair workflow, "
        "and optionally extend or trim the protein before using the result in other tiles."
    )
    st.caption(
        "The original complex is preserved as an input artifact. The prepared outputs "
        "include a protein-only target for docking and a complex that retains the input "
        "ligand and non-water cofactors."
    )
    readiness = ligand_tools.ligand_complex_preparation_readiness()
    if not readiness["ready"]:
        st.warning(
            "Studio's bundled structure-preparation runtime is unavailable at "
            f"{readiness['source']}. Reinstall or repair the Studio package."
        )
    project_id = project_selector(f"ligand-tool-{workflow_id}")
    pending_source = st.session_state.get("ligand-complex-prepare-pending-source")
    pending_path: Path | None = None
    if pending_source:
        try:
            candidate = Path(str(pending_source)).expanduser().resolve(strict=True)
            candidate.relative_to(user_input_root(owner_id))
            if candidate.suffix.lower() not in {".pdb", ".ent"}:
                raise ValueError("The preloaded structure is not a PDB file.")
            pending_path = candidate
            st.info(f"Loaded from the docking tile: {candidate.name}. The source structure will be preserved.")
        except (OSError, ValueError) as exc:
            st.warning(f"The preloaded docking structure is unavailable: {exc}")
            st.session_state.pop("ligand-complex-prepare-pending-source", None)
    uploaded = st.file_uploader(
        "Replace with a different complex (PDB)" if pending_path else "Protein–ligand complex (PDB) *",
        type=["pdb", "ent"],
        key=(
            f"ligand-complex-prepare-input-{workflow_id}-"
            + (hashlib.sha256(str(pending_path).encode()).hexdigest()[:10] if pending_path else "upload")
        ),
        help="Upload a coordinate-bearing complex. Ligand and cofactor records are retained in the complex output.",
    )
    source_name = uploaded.name if uploaded is not None else (pending_path.name if pending_path else "")
    source_content = uploaded.getvalue() if uploaded is not None else None
    if uploaded is None and pending_path is not None:
        try:
            source_content = pending_path.read_bytes()
        except OSError as exc:
            st.error(f"Could not read the preloaded structure: {exc}")
    pdb_data = ""
    chain_ranges: dict[str, dict[str, Any]] = {}
    if source_content is not None:
        try:
            pdb_data = source_content.decode("utf-8")
            chain_ranges = ligand_tools._pdb_polymer_chain_ranges(pdb_data)
        except (UnicodeDecodeError, ValueError) as exc:
            st.error(f"Could not read the complex PDB: {exc}")
    if not chain_ranges and source_content is not None:
        st.warning("No canonical protein chains with readable residue numbers were found.")

    st.subheader("Protein repair")
    repair = st.checkbox(
        "Repair missing atoms and eligible internal gaps",
        value=True,
        key=f"ligand-complex-repair-{workflow_id}",
        help="Uses Studio's bundled MODELLER/PDBFixer preparation workflow. Missing terminal residues are left unchanged unless you add a C-terminal sequence below.",
    )
    repair_cols = st.columns(3)
    map_modified = repair_cols[0].checkbox(
        "Map supported modified residues",
        value=True,
        disabled=not repair,
        key=f"ligand-complex-map-modified-{workflow_id}",
    )
    skip_terminal = repair_cols[1].checkbox(
        "Skip missing terminal residues",
        value=True,
        disabled=not repair,
        key=f"ligand-complex-skip-terminal-{workflow_id}",
    )
    max_gap = int(
        repair_cols[2].number_input(
            "Largest internal gap to rebuild",
            min_value=0,
            max_value=100,
            value=15,
            step=1,
            disabled=not repair,
            key=f"ligand-complex-max-gap-{workflow_id}",
        )
    )

    st.subheader("Optional C-terminal extension")
    add_terminal = st.checkbox(
        "Append a modeled sequence to a protein chain",
        value=False,
        key=f"ligand-complex-add-terminal-{workflow_id}",
    )
    extension_chain = ""
    extension_sequence = ""
    extension_models = 10
    if add_terminal and chain_ranges:
        extension_cols = st.columns(3)
        chains = list(chain_ranges)
        extension_chain = extension_cols[0].selectbox(
            "Chain to extend",
            chains,
            format_func=lambda chain: "Blank chain" if chain == "_" else chain,
            key=f"ligand-complex-extension-chain-{workflow_id}",
        )
        extension_sequence = extension_cols[1].text_input(
            "Residues to append",
            placeholder="e.g. GSG",
            help="Canonical one-letter amino-acid codes, up to 30 residues.",
            key=f"ligand-complex-extension-sequence-{workflow_id}",
        ).strip().upper()
        extension_models = int(
            extension_cols[2].number_input(
                "MODELLER conformations",
                min_value=1,
                max_value=50,
                value=10,
                step=1,
                key=f"ligand-complex-extension-models-{workflow_id}",
            )
        )

    st.subheader("Optional chain trimming")
    trim_protein = st.checkbox(
        "Keep selected residue ranges and trim the rest",
        value=False,
        key=f"ligand-complex-trim-{workflow_id}",
        help="Trimming keeps the selected protein ranges and retains the complex ligand and other non-water heterogens.",
    )
    trim_ranges: dict[str, list[int]] = {}
    if trim_protein and chain_ranges:
        for chain, row in chain_ranges.items():
            st.markdown(f"**Chain {'(blank)' if chain == '_' else chain}**")
            columns = st.columns(3)
            keep_chain = columns[0].checkbox(
                "Keep this chain",
                value=True,
                key=f"ligand-complex-keep-chain-{workflow_id}-{chain}",
            )
            minimum, source_end = int(row["start"]), int(row["end"])
            extra = len(extension_sequence) if add_terminal and chain == extension_chain else 0
            maximum = source_end + extra
            start = int(
                columns[1].number_input(
                    "First residue",
                    min_value=minimum,
                    max_value=maximum,
                    value=minimum,
                    step=1,
                    key=f"ligand-complex-trim-start-{workflow_id}-{chain}",
                )
            )
            end = int(
                columns[2].number_input(
                    "Last residue",
                    min_value=minimum,
                    max_value=maximum,
                    value=maximum,
                    step=1,
                    key=f"ligand-complex-trim-end-{workflow_id}-{chain}",
                )
            )
            if keep_chain:
                trim_ranges[chain] = [start, end]

    project_settings = project_id
    trim_json = json.dumps(trim_ranges, sort_keys=True) if trim_ranges else ""
    parameters = {
        "repair": repair,
        "map_modified_residues": map_modified,
        "skip_terminal_missing_residues": skip_terminal,
        "max_internal_gap": max_gap,
        "extension_chain": extension_chain,
        "extension_sequence": extension_sequence,
        "extension_models": extension_models,
        "trim_ranges_json": trim_json,
    }
    if source_content is not None:
        st.caption(f"Detected protein chains: {', '.join('blank' if key == '_' else key for key in chain_ranges)}")
    disabled = not project_id or source_content is None or not readiness["ready"] or not chain_ranges
    if add_terminal and not extension_sequence:
        disabled = True
    if trim_protein and not trim_ranges:
        disabled = True
    if st.button(
        "Prepare complex",
        type="primary",
        disabled=disabled,
        key=f"ligand-complex-prepare-run-{workflow_id}",
    ) and source_content is not None:
        try:
            previous_job = st.session_state.get("last_submitted_job") or {}
            previous_run_id = previous_job.get("run_id") if isinstance(previous_job, dict) else None
            source_path = _save_target_upload(owner_id, source_name, source_content)
            parameters["source_pdb"] = str(source_path)
            submit(schema, workflow_id, parameters, str(project_settings), {}, title)
            current_job = st.session_state.get("last_submitted_job") or {}
            if isinstance(current_job, dict) and current_job.get("run_id") and current_job.get("run_id") != previous_run_id:
                st.session_state.pop("ligand-complex-prepare-pending-source", None)
        except (OSError, ValueError) as exc:
            st.error(f"Could not start complex preparation: {exc}")


def render_setup(
    workflow_id: str,
    schema: dict[str, Any],
    *,
    owner_id: str,
    compound_sets: list[dict[str, Any]],
    project_selector: Callable[..., str | None],
    submit: Callable[..., None],
    placement_options: Callable[..., dict[str, Any]] | None = None,
    prediction_state: Callable[..., dict[str, Any]] | None = None,
) -> None:
    is_complex_prepare = workflow_id == "ligand.complex.prepare"
    is_vina = workflow_id == "ligand.docking.vina"
    is_autodock_gpu = workflow_id == "ligand.docking.autodock_gpu"
    is_gnina = workflow_id == "ligand.docking.gnina"
    is_unidock = workflow_id == "ligand.docking.unidock_pro"
    is_rosetta = workflow_id == "ligand.docking.rosetta_ligand"
    is_docking = is_vina or is_autodock_gpu or is_gnina or is_unidock or is_rosetta
    is_shared_docking_page = is_vina or is_autodock_gpu or is_gnina or is_unidock or is_rosetta
    is_nesso = workflow_id == "ligand.affinity.nesso"
    is_lddm = workflow_id == "ligand.generation.lddm"
    is_pocketxmol = workflow_id == "ligand.generation.pocketxmol"
    is_generation = workflow_id.startswith("ligand.generation.")
    generator_engine = workflow_id.removeprefix("ligand.generation.") if is_generation else ""
    needs_pocket = generator_engine in {"pocketxmol", "flowr_root", "conditar", "paopt", "drugrpg", "pfm", "pocketflow"}
    needs_reference = generator_engine in {"lddm", "flowr_root", "paopt"}
    is_pgmg = generator_engine == "pgmg"
    is_omtra = generator_engine == "omtra"
    is_gromacs = workflow_id == "ligand.md.gromacs"
    is_openmm = workflow_id == "ligand.md.openmm"
    is_md = is_gromacs or is_openmm
    uses_compound_set = is_docking or is_nesso
    if not (uses_compound_set or is_generation or is_md or is_complex_prepare):
        raise ValueError(f"No standalone tile form is registered for {workflow_id}.")

    generation_titles = {
        "lddm": "LDDM", "pocketxmol": "PocketXMol", "omtra": "OMTRA",
        "flowr_root": "FLOWR.root", "conditar": "conDitar", "paopt": "paOPT",
        "drugrpg": "DrugRPG", "pfm": "PFM", "pocketflow": "PocketFlow", "pgmg": "PGMG",
    }
    title = (
        "Complex Preparation" if is_complex_prepare else
        "AutoDock Vina" if is_vina else
        "AutoDock-GPU" if is_autodock_gpu else
        "GNINA" if is_gnina else
        "Uni-Dock Pro" if is_unidock else
        "RosettaLigand" if is_rosetta else
        "Nesso-1" if is_nesso else
        generation_titles.get(generator_engine, "OpenMM" if is_openmm else "GROMACS")
    )
    st.button("← Back to All tools", key=f"back-to-ligand-tool-{workflow_id}", on_click=lambda: st.session_state.update({"page": "All tools"}))
    if is_shared_docking_page:
        _render_vina_setup(
            workflow_id,
            schema,
            owner_id=owner_id,
            compound_sets=compound_sets,
            project_selector=project_selector,
            submit=submit,
            placement_options=placement_options,
            prediction_state=prediction_state,
            title=title,
        )
        return
    st.title(title)
    if is_complex_prepare:
        _render_complex_preparation_setup(
            workflow_id,
            schema,
            owner_id=owner_id,
            project_selector=project_selector,
            submit=submit,
            title=title,
        )
        return
    if is_docking:
        engine_description = "AutoDock Vina" if is_vina else "AutoDock-GPU with AutoGrid and AutoDock4 scoring" if is_autodock_gpu else "GNINA with CNN pose rescoring" if is_gnina else "GPU-accelerated Uni-Dock Pro" if is_unidock else "RosettaLigand global or reference-guided docking"
        st.write(f"Dock compounds from a Database set using {engine_description}.")
        if is_rosetta:
            readiness = ligand_tools.rosettaligand_runner_readiness()
            st.caption("This tile uses Studio's bundled RosettaLigand/OpenVS runner. Rosetta energy values are reported in REU and should only be ranked within this protocol.")
            if not readiness["ready"]:
                st.warning(f"Studio's bundled RosettaLigand runner is unavailable at {readiness['source']}.")
        else:
            st.caption("Upload a protein PDB, select a compound set, then define the docking box in Ångstroms. The target is prepared with Meeko and ligands are built from SMILES with RDKit.")
        if is_gnina:
            st.caption("GNINA reports an empirical docking score and CNN rescoring metrics; CNN score ranks higher values first.")
        elif is_unidock:
            st.caption("Uni-Dock Pro runs classic receptor docking with fast, balanced, or detailed GPU search.")
        elif is_autodock_gpu:
            st.caption("AutoDock-GPU uses AutoGrid affinity maps and AutoDock4 scoring; rank its kcal/mol scores within this scoring protocol.")
    elif is_nesso:
        st.write("Predict protein–ligand affinity for compounds from a Database set.")
        st.caption("Nesso-1 returns affinity estimates; it does not generate a bound structure or docking pose.")
        readiness = ligand_tools.nesso_readiness()
        if not readiness["ready"]:
            missing = [
                label for key, label in (
                    ("checkpoint_ready", "Nesso checkpoint"),
                    ("ccd_ready", "CCD dictionary"),
                    ("esm_ready", "ESM2 cache"),
                ) if not readiness[key]
            ]
            st.warning("This tile needs reference resources before it can run: " + ", ".join(missing) + ".")
    elif is_generation:
        generation_descriptions = {
            "lddm": "Generate new ligand structures around a known bound ligand with LDDM.",
            "pocketxmol": "Generate de novo 3D ligand structures inside a defined protein pocket with PocketXMol.",
            "omtra": "Generate ligands with OMTRA; an optional pocket, reference ligand, and pharmacophore add structural context.",
            "flowr_root": "Generate or redesign ligands using a pocket and reference ligand to define the binding frame.",
            "conditar": "Generate pocket-conditioned ligands with conDitar.",
            "paopt": "Steer pocket-conditioned conDitar generation toward selected ADMET objectives.",
            "drugrpg": "Generate pocket-conditioned ligands with the physics-guided DrugRPG model.",
            "pfm": "Generate pocket-conditioned ligands with perturbed flow matching.",
            "pocketflow": "Generate chemically constrained molecules autoregressively inside a protein pocket.",
            "pgmg": "Generate compounds from a pharmacophore hypothesis using PGMG.",
        }
        st.write(generation_descriptions[generator_engine])
        if is_lddm:
            st.caption("Uses the MIT-licensed CD+BB checkpoint. Upload a coordinate-bearing reference ligand in SDF format.")
        elif is_pgmg:
            st.caption("Upload a PGMG .posp file with one to eight supported pharmacophore points. Generated structures are also normalized into a Database-ready compound set.")
        elif generator_engine in {"conditar", "paopt"}:
            st.caption("This model uses local reference weights for the permissioned deployment; its container and model sources have redistribution restrictions.")
        elif needs_pocket:
            st.caption("Upload the full target PDB and a pocket-only PDB used by the model to define its generation frame.")
        else:
            st.caption("All generated compound tables and coordinates are saved as typed job artifacts for Database import and downstream pipelines.")
    elif is_md:
        engine_label = "OpenMM" if is_openmm else "GROMACS"
        st.write(f"Prepare a bound protein–ligand complex, run explicit-solvent {engine_label} molecular dynamics, and save its trajectory.")
        st.caption("The ligand SDF must match the selected bound residue's heavy-atom count and element order. The result includes the prepared system, production trajectory, checkpoint, and run summary.")
        readiness = ligand_tools.openmm_md_runner_readiness() if is_openmm else ligand_tools.gromacs_md_runner_readiness()
        if not readiness["ready"]:
            st.warning(f"Studio's bundled {engine_label} runtime is unavailable at {readiness['source']}.")

    project_id = project_selector(f"ligand-tool-{workflow_id}")
    selected_set = None
    if uses_compound_set and compound_sets:
        group_ids = [str(item["group_id"]) for item in compound_sets]
        selected_group = st.selectbox(
            "Database compound set *",
            group_ids,
            format_func=lambda key: next(
                f"{item['name']} · {item['compound_count']:,} compounds"
                for item in compound_sets if str(item["group_id"]) == key
            ),
            key=f"ligand-tool-compound-set-{workflow_id}",
        )
        selected_set = next((item for item in compound_sets if str(item["group_id"]) == selected_group), None)
    elif uses_compound_set:
        st.info("Import a compound set from CSV, TSV, SMI, or SDF in Database before starting this tile.")
    uploaded = st.file_uploader(
        "Protein–ligand complex (PDB) *" if is_md else "Target protein structure (PDB) *",
        type=["pdb", "ent"],
        key=f"ligand-tool-target-{workflow_id}",
    )
    if uploaded is not None and (is_autodock_gpu or is_nesso):
        try:
            checked_target = _save_target_upload(owner_id, uploaded.name, uploaded.getvalue())
            render_receptor_check(
                checked_target.read_text(encoding="utf-8", errors="replace"),
                checked_target,
                owner_id=owner_id,
                key=f"ligand-{workflow_id}-repair-receptor",
            )
        except (OSError, ValueError) as exc:
            st.error(f"Could not inspect the target structure: {exc}")
    bound_ligands: list[dict[str, Any]] = []
    saved_complex: Path | None = None
    selected_ligand_key = ""
    if is_md and uploaded:
        try:
            saved_complex = _save_target_upload(owner_id, uploaded.name, uploaded.getvalue())
            bound_ligands = ligand_tools.bound_ligand_choices(saved_complex)
            if bound_ligands:
                selected_ligand_key = st.selectbox(
                    "Bound ligand residue *",
                    [str(item["key"]) for item in bound_ligands],
                    format_func=lambda key: next(
                        f"{item['resname']} · chain {item['chain']} · residue {item['resseq']} · {item['heavy_atom_count']} heavy atoms"
                        for item in bound_ligands if item["key"] == key
                    ),
                    key="ligand-gromacs-bound-residue",
                )
            else:
                st.warning("No non-water bound ligand residues were found in this complex PDB.")
        except (OSError, ValueError) as exc:
            st.error(f"Could not inspect the complex PDB: {exc}")
    uploaded_md_ligand = st.file_uploader(
        "Matching bound-ligand coordinates (SDF) *",
        type=["sdf"],
        key=f"ligand-tool-md-ligand-{workflow_id}",
    ) if is_md else None
    uploaded_pocket = st.file_uploader(
        "Pocket-only structure (PDB) *" if needs_pocket else "Pocket-only structure (PDB, optional)",
        type=["pdb", "ent"],
        key=f"ligand-tool-pocket-{workflow_id}",
    ) if (needs_pocket or generator_engine == "omtra") else None
    uploaded_reference = st.file_uploader(
        "Reference ligand coordinates (SDF) *" if needs_reference else "Reference ligand coordinates (SDF)" if (is_rosetta or generator_engine == "omtra") else "Reference ligand coordinates (SDF)",
        type=["sdf"],
        key=f"ligand-tool-reference-{workflow_id}",
    ) if (needs_reference or is_rosetta or generator_engine == "omtra") else None
    uploaded_pharmacophore = st.file_uploader(
        "PGMG pharmacophore (.posp) *", type=["posp"], key=f"ligand-tool-pharmacophore-{workflow_id}"
    ) if is_pgmg else None
    uploaded_omtra_pharmacophore = st.file_uploader(
        "OMTRA Pharmit pharmacophore (JSON, optional)", type=["json"], key=f"ligand-tool-omtra-pharmacophore-{workflow_id}"
    ) if is_omtra else None

    parameters: dict[str, Any] = {}
    settings: dict[str, Any] = {}
    alignment_ligand_upload = None
    if is_rosetta:
        st.subheader("RosettaLigand placement")
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
        if reference_mode == "reference_guided":
            st.caption("Upload a coordinate-bearing ligand SDF that defines the placement anchor.")
        center_cols = st.columns(3)
        center = {
            f"center_{axis}": center_cols[index].number_input(f"Pocket center {axis.upper()} (Å)", value=0.0, step=1.0, key=f"ligand-rosetta-center-{axis}")
            for index, axis in enumerate("xyz")
        }
        with st.expander("Preparation and sampling", expanded=False):
            cols = st.columns(3)
            compound_limit = int(cols[0].number_input("Compound limit (0 = all selected)", min_value=0, max_value=5000, value=0, step=1, key="ligand-rosetta-limit"))
            cpu_cores = int(cols[1].number_input("CPU workers", min_value=1, max_value=128, value=8, step=1, key="ligand-rosetta-cpu"))
            replicates = int(cols[2].number_input("Replicates", min_value=1, max_value=100, value=1, step=1, key="ligand-rosetta-replicates"))
            cols = st.columns(3)
            ph = float(cols[0].number_input("Ligand protonation pH", min_value=0.0, max_value=14.0, value=7.4, step=0.1, key="ligand-rosetta-ph"))
            conformers = int(cols[1].number_input("Input conformers", min_value=1, max_value=1000, value=20, step=1, key="ligand-rosetta-conformers"))
            minimization_steps = int(cols[2].number_input("Minimization steps", min_value=1, max_value=100000, value=2000, step=100, key="ligand-rosetta-min-steps"))
            padding = float(st.number_input("Search padding (Å)", min_value=1.0, max_value=30.0, value=4.0, step=0.5, key="ligand-rosetta-padding"))
        parameters = {
            **center, "reference_mode": reference_mode, "protocol": protocol,
            "compound_limit": compound_limit, "cpu_cores": cpu_cores,
            "replicates": replicates, "seed_start": 1001, "ph": ph,
            "conformers": conformers, "minimization_steps": minimization_steps,
            "padding": padding,
        }
    elif is_docking:
        tool_key = "unidock-pro" if is_unidock else "autodock-gpu" if is_autodock_gpu else "gnina" if is_gnina else "vina"
        st.subheader("Docking box")
        center_cols = st.columns(3)
        center = {
            "center_x": center_cols[0].number_input("Center X (Å)", value=0.0, step=1.0, key=f"ligand-{tool_key}-center-x"),
            "center_y": center_cols[1].number_input("Center Y (Å)", value=0.0, step=1.0, key=f"ligand-{tool_key}-center-y"),
            "center_z": center_cols[2].number_input("Center Z (Å)", value=0.0, step=1.0, key=f"ligand-{tool_key}-center-z"),
        }
        box_cols = st.columns(3)
        sizes = {
            "size_x": box_cols[0].number_input("Size X (Å)", min_value=1.0, max_value=60.0, value=20.0, step=1.0, key=f"ligand-{tool_key}-size-x"),
            "size_y": box_cols[1].number_input("Size Y (Å)", min_value=1.0, max_value=60.0, value=20.0, step=1.0, key=f"ligand-{tool_key}-size-y"),
            "size_z": box_cols[2].number_input("Size Z (Å)", min_value=1.0, max_value=60.0, value=20.0, step=1.0, key=f"ligand-{tool_key}-size-z"),
        }
        with st.expander("Search settings", expanded=False):
            search_cols = st.columns(3)
            if is_unidock:
                search_mode = search_cols[0].selectbox("Search effort", ["fast", "balance", "detail"], index=1, key=f"ligand-{tool_key}-search-mode")
            elif is_autodock_gpu:
                nrun = int(search_cols[0].number_input("LGA runs per compound", min_value=1, max_value=500, value=20, step=1, key=f"ligand-{tool_key}-nrun"))
            else:
                exhaustiveness = int(search_cols[0].number_input("Exhaustiveness", min_value=1, max_value=1024, value=8, step=1, key=f"ligand-{tool_key}-exhaustiveness"))
            num_modes = 9 if is_autodock_gpu else int(search_cols[1].number_input("Poses per compound", min_value=1, max_value=100, value=9, step=1, key=f"ligand-{tool_key}-modes"))
            compound_limit = int(search_cols[1 if is_autodock_gpu else 2].number_input(
                "Compound limit (0 = all selected)", min_value=0,
                max_value=ligand_tools.MAX_COMPOUNDS_PER_LAUNCH if is_unidock else 5000,
                value=0, step=100 if is_unidock else 10, key=f"ligand-{tool_key}-limit",
            ))
            cpu_label = "CPU cores / Vina compound workers" if is_vina else "CPU cores"
            cpu_help = (
                "This is the total CPU allocation. Vina docks up to this many compounds concurrently, using one CPU thread per compound."
                if is_vina else None
            )
            cpu_cores = int(st.number_input(
                cpu_label, min_value=1, max_value=256, value=4, step=1,
                key=f"ligand-{tool_key}-cpu", help=cpu_help,
            ))
            gpu_device = st.text_input("GPU device", value="0", key=f"ligand-{tool_key}-gpu", help="Use a device ID such as 0.") if (is_gnina or is_unidock or is_autodock_gpu) else "none"
            compound_workers = 1
            if is_gnina:
                compound_workers = int(st.number_input(
                    "Concurrent compounds", min_value=1,
                    max_value=max(1, min(4, cpu_cores, compound_limit)),
                    value=max(1, min(2, cpu_cores, compound_limit)), step=1,
                    key="ligand-tool-gnina-compound-workers",
                    help="Run independent compounds concurrently on the assigned GPU. CPU threads are divided across workers; use 1 if GPU memory is tight.",
                ))
            use_scrub = True
            scrub_ph = 7.4
            scrub_skip_tautomer = True
            if is_autodock_gpu:
                st.markdown("**Ligand preparation**")
                prep_cols = st.columns(3)
                use_scrub = prep_cols[0].checkbox(
                    "Use Scrub", value=True, key=f"ligand-{tool_key}-use-scrub",
                    help="Normalize ligand protonation and chemistry with Scrub before Meeko prepares the ligand for docking.",
                )
                scrub_ph = float(prep_cols[1].number_input(
                    "Scrub pH", min_value=0.0, max_value=14.0, value=7.4, step=0.1,
                    disabled=not use_scrub, key=f"ligand-{tool_key}-scrub-ph",
                ))
                scrub_skip_tautomer = prep_cols[2].checkbox(
                    "Skip tautomer enumeration", value=True, disabled=not use_scrub,
                    key=f"ligand-{tool_key}-skip-tautomer",
                    help="When cleared, Scrub may generate tautomer states before docking.",
                )
        align_longest_axis_to_x = False
        if not is_rosetta:
            st.subheader("Docking coordinate frame")
            align_longest_axis_to_x = st.checkbox(
                "Align a reference ligand's longest axis with global X",
                value=False,
                key=f"ligand-{tool_key}-align-axis-x",
                help="Rotate the receptor and docking box together around the uploaded ligand's centroid.",
            )
            if align_longest_axis_to_x:
                alignment_ligand_upload = st.file_uploader(
                    "Coordinate-bearing alignment ligand (SDF) *",
                    type=["sdf"],
                    key=f"ligand-{tool_key}-alignment-ligand",
                    help="This ligand defines the axis only; it is not docked unless it is also included in the Database compound set.",
                )
                st.caption("Docking poses and the staged receptor will share the rotated coordinate frame. The transform is saved with the results.")
        parameters = {
            **center, **sizes,
            **({"search_mode": search_mode} if is_unidock else {"nrun": nrun} if is_autodock_gpu else {"exhaustiveness": exhaustiveness}),
            **({} if is_autodock_gpu else {"num_modes": num_modes}),
            "compound_limit": compound_limit, "cpu_cores": cpu_cores,
            "seed": 42, "gpu_device": gpu_device,
        }
        if not is_rosetta:
            parameters["align_longest_axis_to_x"] = align_longest_axis_to_x
        if is_gnina:
            parameters["compound_workers"] = compound_workers
        if is_autodock_gpu:
            parameters.update({
                "use_scrub": use_scrub,
                "scrub_ph": scrub_ph,
                "scrub_skip_tautomer": scrub_skip_tautomer,
            })
    elif is_md:
        with st.expander("Simulation settings", expanded=True):
            md_cols = st.columns(3)
            production_ns = float(md_cols[0].number_input("Production duration (ns)", min_value=0.001, max_value=1000.0, value=1.0, step=0.1, key="ligand-gromacs-duration"))
            timestep_fs = float(md_cols[1].selectbox("Timestep (fs)", [2.0, 4.0], key="ligand-gromacs-timestep"))
            temperature = float(md_cols[2].number_input("Temperature (K)", min_value=200.0, max_value=400.0, value=300.0, step=1.0, key="ligand-gromacs-temperature"))
            md_cols = st.columns(3)
            pressure = float(md_cols[0].number_input("Pressure (bar)", min_value=0.1, max_value=20.0, value=1.0, step=0.1, key="ligand-gromacs-pressure"))
            padding_nm = float(md_cols[1].number_input("Solvent padding (nm)", min_value=0.6, max_value=3.0, value=1.0, step=0.1, key="ligand-gromacs-padding"))
            ionic_strength = float(md_cols[2].number_input("Ionic strength (M)", min_value=0.0, max_value=2.0, value=0.15, step=0.05, key="ligand-gromacs-ionic-strength"))
            md_cols = st.columns(3)
            engine_key = "openmm" if is_openmm else "gromacs"
            ligand_forcefields = ["openff-2.2.0", "gaff2"] if is_openmm else ["gaff2", "gaff"]
            protein_forcefields = ["amber14-all", "amber19-all"] if is_openmm else ["ff14SB", "ff19SB", "ff15ipq", "ff03.r1"]
            forcefield_method = md_cols[0].selectbox("Ligand force field", ligand_forcefields, key=f"ligand-{engine_key}-ligand-ff")
            protein_forcefield_method = md_cols[1].selectbox("Protein force field", protein_forcefields, key=f"ligand-{engine_key}-protein-ff")
            water_model = md_cols[2].selectbox("Water model", ["tip3p", "opc", "spce", "tip4pew"], key=f"ligand-{engine_key}-water")
            hydrogen_mass_repartitioning = st.checkbox("Hydrogen-mass repartitioning for 4 fs", value=True, disabled=timestep_fs != 4.0, key=f"ligand-{engine_key}-hmr")
            use_gpu = st.checkbox("Use GPU acceleration", value=True, key=f"ligand-{engine_key}-use-gpu")
            gpu_device = st.text_input("GPU device", value="0", disabled=not use_gpu, key=f"ligand-{engine_key}-gpu") if use_gpu else "none"
            cpu_cores = int(st.number_input("CPU cores", min_value=1, max_value=256, value=8, step=1, key=f"ligand-{engine_key}-cpu"))
            seed = int(st.number_input("Random seed", min_value=1, max_value=2_147_483_646, value=42, step=1, key=f"ligand-{engine_key}-seed"))
        parameters = {
            "ligand_key": selected_ligand_key,
            "production_ns": production_ns,
            "timestep_fs": timestep_fs,
            "temperature": temperature,
            "pressure": pressure,
            "padding_nm": padding_nm,
            "ionic_strength": ionic_strength,
            "forcefield_method": forcefield_method,
            "protein_forcefield_method": protein_forcefield_method,
            "water_model": water_model,
            "hydrogen_mass_repartitioning": hydrogen_mass_repartitioning,
            "use_gpu": use_gpu,
            "gpu_device": gpu_device.strip() or "0",
            "cpu_cores": cpu_cores,
            "seed": seed,
        }
    elif is_nesso:
        with st.expander("Nesso settings", expanded=False):
            settings_cols = st.columns(3)
            recycling_steps = int(settings_cols[0].number_input("Recycling steps", min_value=1, max_value=64, value=5, step=1, key="ligand-nesso-recycles"))
            seed = int(settings_cols[1].number_input("Seed", min_value=1, max_value=2_147_483_646, value=42, step=1, key="ligand-nesso-seed"))
            compound_limit = int(settings_cols[2].number_input("Compound limit (0 = all selected)", min_value=0, max_value=5000, value=0, step=10, key="ligand-nesso-limit"))
        parameters = {
            "recycling_steps": recycling_steps, "seed": seed,
            "compound_limit": compound_limit, "gpu_device": "0",
        }
    else:
        with st.expander("Generation settings", expanded=False):
            setting_cols = st.columns(3)
            requested_count = int(setting_cols[0].number_input("Molecules to generate", min_value=1, max_value=10000, value=64, step=8, key=f"ligand-{workflow_id}-count"))
            batch_size = int(setting_cols[1].number_input("Batch size", min_value=1, max_value=512, value=8, step=1, key=f"ligand-{workflow_id}-batch"))
            seed = int(setting_cols[2].number_input("Seed", min_value=1, max_value=2_147_483_646, value=42, step=1, key=f"ligand-{workflow_id}-seed"))
            engine_settings: dict[str, Any] = {}
            if is_lddm:
                setting_cols = st.columns(3)
                n_steps = int(setting_cols[0].number_input("Integration steps", min_value=1, max_value=2000, value=100, step=10, key="ligand-lddm-steps"))
                sampler = setting_cols[1].selectbox("Sampler", ["ForwardEuler", "HeunSampler"], key="ligand-lddm-sampler")
                sampling_noise = float(setting_cols[2].number_input("Sampling noise", min_value=0.0, max_value=20.0, value=5.0, step=0.5, key="ligand-lddm-noise"))
            elif is_pocketxmol:
                setting_cols = st.columns(3)
                diffusion_steps = int(setting_cols[0].number_input("Diffusion steps", min_value=1, max_value=2000, value=100, step=10, key="ligand-pocketxmol-steps"))
                ligand_atoms_mean = float(setting_cols[1].number_input("Mean ligand atom count", min_value=5.0, max_value=100.0, value=28.0, step=1.0, key="ligand-pocketxmol-atoms-mean"))
                ligand_atoms_std = float(setting_cols[2].number_input("Atom-count deviation", min_value=0.0, max_value=20.0, value=2.0, step=0.5, key="ligand-pocketxmol-atoms-std"))
            elif generator_engine == "omtra":
                setting_cols = st.columns(3)
                engine_settings["integration_steps"] = int(setting_cols[0].number_input("Integration steps", min_value=25, max_value=2000, value=250, step=25, key="ligand-omtra-steps"))
                engine_settings["noise_scale"] = float(setting_cols[1].number_input("Noise scale", min_value=0.0, max_value=5.0, value=1.0, step=0.1, key="ligand-omtra-noise"))
                engine_settings["epsilon"] = float(setting_cols[2].number_input("Epsilon", min_value=0.0001, max_value=0.2, value=0.01, step=0.001, format="%.4f", key="ligand-omtra-epsilon"))
                engine_settings["stochastic_sampling"] = st.checkbox("Use stochastic sampling", value=False, key="ligand-omtra-stochastic")
            elif generator_engine == "flowr_root":
                setting_cols = st.columns(3)
                engine_settings["integration_steps"] = int(setting_cols[0].number_input("Integration steps", min_value=10, max_value=2000, value=100, step=10, key="ligand-flowr-steps"))
                engine_settings["corrector_steps"] = int(setting_cols[1].number_input("Corrector steps", min_value=0, max_value=100, value=0, step=1, key="ligand-flowr-corrector"))
                engine_settings["solver"] = setting_cols[2].selectbox("Solver", ["euler", "heun"], key="ligand-flowr-solver")
                engine_settings["diversity_threshold"] = float(st.number_input("Diversity threshold", min_value=0.0, max_value=1.0, value=0.9, step=0.05, key="ligand-flowr-diversity"))
                options_cols = st.columns(3)
                engine_settings["use_sde_simulation"] = options_cols[0].checkbox("Use SDE simulation", value=False, key="ligand-flowr-sde")
                engine_settings["sample_molecule_sizes"] = options_cols[1].checkbox("Sample molecule sizes", value=False, key="ligand-flowr-size-sampling")
                engine_settings["filter_diversity"] = options_cols[2].checkbox("Filter for diversity", value=False, key="ligand-flowr-filter-diversity")
            elif generator_engine in {"conditar", "paopt"}:
                setting_cols = st.columns(2)
                engine_settings["diffusion_steps"] = int(setting_cols[0].number_input("Diffusion steps", min_value=100, max_value=5000, value=1000, step=100, key=f"ligand-{generator_engine}-steps"))
                engine_settings["pocket_radius"] = float(setting_cols[1].number_input("Pocket radius (Å)", min_value=4.0, max_value=24.0, value=10.0, step=0.5, key=f"ligand-{generator_engine}-radius"))
                if generator_engine == "paopt":
                    endpoints = ["BBBP", "HIA", "hERG", "Ames", "Carcinogenicity", "CaCo2"]
                    engine_settings["optimize_properties"] = st.multiselect("ADMET endpoints to steer", endpoints, default=["Carcinogenicity"], key="ligand-paopt-properties")
                    engine_settings["minimize_properties"] = st.multiselect("Endpoints to minimize", endpoints, default=["Carcinogenicity"], key="ligand-paopt-minimize")
                    setting_cols = st.columns(3)
                    engine_settings["optimization_steps"] = int(setting_cols[0].number_input("paOPT steering steps", min_value=1, max_value=10, value=1, step=1, key="ligand-paopt-opt-steps"))
                    engine_settings["gradient_estimate_pairs"] = int(setting_cols[1].number_input("Gradient estimate pairs", min_value=1, max_value=16, value=4, step=1, key="ligand-paopt-pairs"))
                    engine_settings["perturbation_size"] = float(setting_cols[2].number_input("Perturbation size", min_value=0.001, max_value=0.2, value=0.03, step=0.005, format="%.3f", key="ligand-paopt-perturbation"))
            elif generator_engine == "drugrpg":
                engine_settings["max_atoms"] = int(st.number_input("Maximum ligand atoms", min_value=5, max_value=100, value=30, step=1, key="ligand-drugrpg-max-atoms"))
            elif generator_engine == "pocketflow":
                setting_cols = st.columns(3)
                engine_settings["atom_temperature"] = float(setting_cols[0].number_input("Atom temperature", min_value=0.1, max_value=5.0, value=1.0, step=0.1, key="ligand-pocketflow-atom-temp"))
                engine_settings["bond_temperature"] = float(setting_cols[1].number_input("Bond temperature", min_value=0.1, max_value=5.0, value=1.0, step=0.1, key="ligand-pocketflow-bond-temp"))
                engine_settings["max_atoms"] = int(setting_cols[2].number_input("Maximum ligand atoms", min_value=5, max_value=100, value=40, step=1, key="ligand-pocketflow-max-atoms"))
                setting_cols = st.columns(3)
                engine_settings["focus_strategy"] = setting_cols[0].selectbox("Pocket focus strategy", ["maximum", "random"], key="ligand-pocketflow-focus")
                engine_settings["focus_threshold"] = float(setting_cols[1].number_input("Focus threshold", min_value=0.0, max_value=1.0, value=0.5, step=0.05, key="ligand-pocketflow-threshold"))
                engine_settings["min_protein_distance"] = float(setting_cols[2].number_input("Minimum protein distance (Å)", min_value=0.0, max_value=10.0, value=3.0, step=0.25, key="ligand-pocketflow-distance"))
            resource_cols = st.columns(2)
            gpu_device = resource_cols[0].text_input("GPU device", value="0", key=f"ligand-{workflow_id}-gpu", help="Use a device ID such as 0.")
            cpu_cores = int(resource_cols[1].number_input("CPU cores", min_value=1, max_value=256, value=8, step=1, key=f"ligand-{workflow_id}-cpu"))
        parameters = {
            "requested_count": requested_count,
            "batch_size": batch_size,
            "seed": seed,
            "gpu_device": gpu_device.strip() or "0",
            "cpu_cores": cpu_cores,
        }
        if is_lddm:
            parameters.update({"n_steps": n_steps, "sampler": sampler, "sampling_noise": sampling_noise})
        elif is_pocketxmol:
            parameters.update({
                "diffusion_steps": diffusion_steps,
                "ligand_atoms_mean": ligand_atoms_mean,
                "ligand_atoms_std": ligand_atoms_std,
            })
        else:
            parameters["engine_settings"] = engine_settings
    disabled = (
        not project_id
        or not uploaded
        or (uses_compound_set and not selected_set)
        or (needs_pocket and not uploaded_pocket)
        or (needs_reference and not uploaded_reference)
        or (is_pgmg and not uploaded_pharmacophore)
        or (is_md and (not uploaded_md_ligand or not selected_ligand_key))
        or (is_gromacs and ligand_tools.gromacs_md_runner_readiness()["ready"] is False)
        or (is_openmm and ligand_tools.openmm_md_runner_readiness()["ready"] is False)
        or (is_rosetta and ligand_tools.rosettaligand_runner_readiness()["ready"] is False)
        or (is_rosetta and parameters.get("reference_mode") == "reference_guided" and not uploaded_reference)
        or (is_docking and not is_rosetta and parameters.get("align_longest_axis_to_x") and not alignment_ligand_upload)
        or (is_nesso and not ligand_tools.nesso_readiness()["ready"])
    )
    if st.button("Run tool", type="primary", disabled=disabled, key=f"ligand-tool-run-{workflow_id}"):
        try:
            if is_md:
                target_path = saved_complex or _save_target_upload(owner_id, uploaded.name, uploaded.getvalue())
                parameters.update({"complex_pdb": str(target_path)})
                ligand_path = _save_reference_ligand_upload(owner_id, uploaded_md_ligand.name, uploaded_md_ligand.getvalue())
                parameters.update({"ligand_sdf": str(ligand_path)})
            else:
                target_path = _save_target_upload(owner_id, uploaded.name, uploaded.getvalue())
                parameters.update({"target_pdb": str(target_path)})
            if is_docking and not is_rosetta and parameters.get("align_longest_axis_to_x") and alignment_ligand_upload:
                alignment_path = _save_reference_ligand_upload(
                    owner_id,
                    alignment_ligand_upload.name,
                    alignment_ligand_upload.getvalue(),
                )
                ligand_tools._sdf_heavy_atom_coordinates(alignment_path)
                parameters["alignment_ligand_sdf"] = str(alignment_path)
            if uses_compound_set:
                parameters["compound_set_file"] = str(selected_set["path"])
            if is_generation and uploaded_pocket:
                pocket_path = _save_target_upload(owner_id, uploaded_pocket.name, uploaded_pocket.getvalue())
                parameters["pocket_pdb"] = str(pocket_path)
            if is_generation and uploaded_reference:
                reference_path = _save_reference_ligand_upload(owner_id, uploaded_reference.name, uploaded_reference.getvalue())
                parameters["reference_ligand_sdf"] = str(reference_path)
            elif is_rosetta and uploaded_reference:
                reference_path = _save_reference_ligand_upload(owner_id, uploaded_reference.name, uploaded_reference.getvalue())
                parameters["reference_ligand_sdf"] = str(reference_path)
            if is_pgmg and uploaded_pharmacophore:
                pharmacophore_path = _save_pgmg_pharmacophore_upload(owner_id, uploaded_pharmacophore.name, uploaded_pharmacophore.getvalue())
                parameters["pharmacophore_posp"] = str(pharmacophore_path)
            if is_omtra and uploaded_omtra_pharmacophore:
                pharmacophore_path = _save_omtra_pharmacophore_upload(owner_id, uploaded_omtra_pharmacophore.name, uploaded_omtra_pharmacophore.getvalue())
                parameters["pharmacophore_json"] = str(pharmacophore_path)
            submit(schema, workflow_id, parameters, str(project_id), settings, title)
        except (OSError, ValueError) as exc:
            st.error(f"Could not start {title}: {exc}")
