"""Studio input page for BioEmu-2 ensemble sampling."""
from __future__ import annotations

import hashlib
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.auth import user_input_root


BIOEMU2_EXAMPLES = {
    "Ubiquitin": (
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    ),
    "Adenylate kinase + ATP + Mg": (
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG",
    ),
    "SHP2 SH2 + peptides (contact)": (
        "GSHMTSRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLNCAD",
        "RVIYFVPLNR",
    ),
}


def _store_structures(uploads, owner_id: str) -> list[str]:
    if not uploads:
        raise ValueError("Upload at least one starting complex structure.")
    if len(uploads) > 32:
        raise ValueError("Upload at most 32 target–design structures in one BioEmu-2 batch.")
    total_bytes = 0
    stored: list[str] = []
    seen_names: set[str] = set()
    for upload in uploads:
        name = Path(upload.name).name
        suffix = Path(name).suffix.lower()
        if suffix not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise ValueError(f"Unsupported structure format for {name}. Upload PDB or mmCIF.")
        if name.casefold() in seen_names:
            raise ValueError("Uploaded structure filenames must be unique.")
        seen_names.add(name.casefold())
        content = upload.getvalue()
        total_bytes += len(content)
        if len(content) > 50_000_000:
            raise ValueError(f"{name} is larger than the 50 MB per-file limit.")
        if total_bytes > 100_000_000:
            raise ValueError("The total upload size cannot exceed 100 MB.")
        digest = hashlib.sha256(content).hexdigest()
        destination = user_input_root(owner_id) / "bioemu2" / digest[:16] / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists():
            destination.write_bytes(content)
        stored.append(str(destination))
    return stored


def _parse_sequences(text: str) -> list[str]:
    lines = [line.strip() for line in str(text or "").splitlines() if line.strip()]
    if not any(line.startswith(">") for line in lines):
        return ["".join(line.split()).upper() for line in lines]
    sequences: list[str] = []
    current: list[str] = []
    for stripped in lines:
        if stripped.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
            continue
        current.append("".join(stripped.split()).upper())
    if current:
        sequences.append("".join(current))
    return sequences


def _new_component_id() -> int:
    component_id = int(st.session_state.get("bioemu2-next-component-id", 1))
    st.session_state["bioemu2-next-component-id"] = component_id + 1
    return component_id


def _component_key(component_id: int, field: str) -> str:
    return f"bioemu2-component-{component_id}-{field}"


def _replace_sequence_components(sequences: list[str]) -> None:
    for component_id in st.session_state.get("bioemu2-component-ids", []):
        st.session_state.pop(_component_key(component_id, "sequence"), None)
        st.session_state.pop(_component_key(component_id, "copies"), None)
    component_ids = []
    for sequence in sequences:
        component_id = _new_component_id()
        component_ids.append(component_id)
        st.session_state[_component_key(component_id, "sequence")] = sequence
        st.session_state[_component_key(component_id, "copies")] = 1
    st.session_state["bioemu2-component-ids"] = component_ids


def _load_sequence_example(name: str) -> None:
    st.session_state["bioemu2-task-mode-v2"] = "Equilibrium Sampling"
    _replace_sequence_components(list(BIOEMU2_EXAMPLES[name]))


def _add_sequence_component() -> None:
    component_id = _new_component_id()
    st.session_state["bioemu2-component-ids"] = [
        *st.session_state.get("bioemu2-component-ids", []), component_id
    ]
    st.session_state[_component_key(component_id, "sequence")] = ""
    st.session_state[_component_key(component_id, "copies")] = 1


def _duplicate_sequence_component(component_id: int) -> None:
    current_ids = list(st.session_state.get("bioemu2-component-ids", []))
    new_id = _new_component_id()
    st.session_state["bioemu2-component-ids"] = [*current_ids, new_id]
    st.session_state[_component_key(new_id, "sequence")] = st.session_state.get(
        _component_key(component_id, "sequence"), ""
    )
    st.session_state[_component_key(new_id, "copies")] = 1


def _remove_sequence_component(component_id: int) -> None:
    st.session_state["bioemu2-component-ids"] = [
        value for value in st.session_state.get("bioemu2-component-ids", [])
        if value != component_id
    ]


def render_setup(schema: dict, *, owner_id: str, user: dict, project_selector, submit, source_job_options) -> None:
    if st.button("← Back to All tools", key="back-to-bioemu2-tools"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.session_state.pop("selected_workflow_group", None)
        st.rerun()

    st.html("""<style>
      .block-container{max-width:800px!important;padding:1rem 1.15rem 2rem!important;margin:0 auto!important}
      [data-testid="stMainBlockContainer"]{background:#fff;border:1px solid #d9ddda;border-radius:11px;box-shadow:0 2px 8px rgba(16,24,40,.08);padding:1.1rem 1.35rem 1.5rem}
      .bioemu2-heading{text-align:center;margin:.15rem auto .8rem}
      .bioemu2-heading h1{font-size:1.9rem;color:#17191b;margin:0 0 .35rem;font-weight:650;letter-spacing:-.025em}
      .bioemu2-heading .badge{display:inline-block;vertical-align:middle;margin-left:.35rem;padding:.16rem .48rem;border-radius:1rem;background:#28674d;color:#fff;font-size:.62rem;font-weight:700;letter-spacing:.02em}
      .bioemu2-heading .hero-link{display:inline-block;vertical-align:middle;margin-left:.4rem;color:#68756e;text-decoration:none;font-size:.84rem}
      .bioemu2-heading .meta{color:#69767f;font-size:.82rem;margin:.55rem 0 .2rem}
      .bioemu2-heading p{color:#52616a;max-width:700px;margin:.55rem auto .65rem;line-height:1.5;font-size:.92rem}
      .bioemu2-heading a{color:#526579;text-decoration:none;font-size:.9rem}
      .bioemu2-note{padding:.75rem .95rem;border-radius:7px;background:#edf5f0;color:#315d47;line-height:1.55}
      .bioemu2-project-reminder{font-size:.78rem;color:#7b858b;margin:-.15rem 0 .65rem}
      .bioemu2-examples-title{font-size:.78rem;font-weight:700;letter-spacing:.05em;color:#33424b;margin:.7rem 0 .45rem}
      .st-key-bioemu2-task-mode-v2 button[aria-pressed="true"]{background:#e8ece9!important;color:#20332a!important;border-color:#cfd8d2!important;box-shadow:none!important}
      .st-key-bioemu2-task-mode-v2 button[aria-pressed="false"]{background:#fff!important;color:#40505a!important;border-color:#d5d9d7!important}
      [data-testid="stVerticalBlockBorderWrapper"]{border-color:#d9ddda!important;border-radius:9px!important}
      [data-testid="stButton"] button{border-radius:7px!important}
      button[kind="primary"]{background:#28674d!important;border-color:#28674d!important;color:white!important;border-radius:7px!important}
      button[kind="primary"]:disabled{background:#91b5a4!important;border-color:#91b5a4!important;color:white!important}
      [data-testid="stTextArea"] textarea{border-radius:7px!important}
    </style>""")
    st.markdown(
        '<div class="bioemu2-heading"><h1>BioEmu 2 <span class="badge">NEW</span>'
        '<a class="hero-link" href="https://github.com/microsoft/bioemu2/" target="_blank" title="Source code" aria-label="Source code">GH</a>'
        '<a class="hero-link" href="https://github.com/microsoft/bioemu2/blob/main/MODEL_CARD.md" target="_blank" title="Paper and model card" aria-label="Paper and model card">▤</a></h1>'
        '<div class="meta">Molecular Dynamics · ~10 min · BioEmu-2 0.1.0rc0</div>'
        '<p>Predict structures, sample ensembles, and propagate complexes in time.</p></div>',
        unsafe_allow_html=True,
    )
    action_cols = st.columns([1, 1, 1, .8])
    with action_cols[0]:
        with st.popover("▤ Show Outputs"):
            st.markdown("**Propagation**")
            st.write("Produces sampled trajectory frames, a topology structure, and a physicality summary for each starting complex.")
            st.markdown("**Sequence modes**")
            st.write("Produce predicted or sampled structures and a run summary. Full files remain available from the job results page.")
    with action_cols[1]:
        with st.popover("⚙ Show Settings"):
            st.write("Choose the sample count and sampler below. Advanced controls include GPU batch size, precision, seed, denoising steps, and physicality profile.")
    with action_cols[2]:
        with st.popover("↗ API Docs"):
            st.markdown("[BioEmu-2 GitHub](https://github.com/microsoft/bioemu2/) · [Model card](https://github.com/microsoft/bioemu2/blob/main/MODEL_CARD.md)")
    with action_cols[3]:
        with st.popover("⌄ Details"):
            st.write(
                "This tile supports equilibrium sampling and structure prediction from protein chain sequences, plus time-lagged propagation from saved or uploaded structures. "
                "The current workflow accepts protein chains only; ligand/cofactor components and contact restraints shown in some upstream examples are not submitted."
            )
    st.divider()
    name_col, project_col = st.columns([1.2, 1.2])
    with name_col:
        job_name = st.text_input("Job Name", value="bioemu2-sampling", key="bioemu2-job-name").strip()
    with project_col:
        project_id = project_selector("bioemu2-project")
    st.markdown(
        '<div class="bioemu2-project-reminder">Recommended: keep related design, prediction, and analysis jobs in the same project so their inputs and results stay together.</div>',
        unsafe_allow_html=True,
    )

    mode_label = st.pills(
        "",
        ["Equilibrium Sampling", "Structure Prediction", "Propagate Structure"],
        selection_mode="single",
        default="Propagate Structure",
        key="bioemu2-task-mode-v2",
        help="Select saved target–design complexes to propagate, or enter protein chain sequences for equilibrium sampling or structure prediction.",
        label_visibility="collapsed",
    ) or "Propagate Structure"
    mode = {
        "Propagate Structure": "propagate",
        "Equilibrium Sampling": "sample-equilibrium",
        "Structure Prediction": "predict-structure",
    }[mode_label]
    st.caption({
        "propagate": "Propagate a saved or uploaded protein complex in time.",
        "sample-equilibrium": "Sample the equilibrium ensemble for the supplied protein sequence(s).",
        "predict-structure": "Predict folded structures for the supplied protein sequence(s).",
    }[mode])
    st.markdown('<div class="bioemu2-examples-title">EXAMPLES</div>', unsafe_allow_html=True)
    example_columns = st.columns(len(BIOEMU2_EXAMPLES))
    for column, example_name in zip(example_columns, BIOEMU2_EXAMPLES):
        with column:
            st.button(
                example_name,
                use_container_width=True,
                key=f"bioemu2-example-{example_name}",
                on_click=_load_sequence_example,
                args=(example_name,),
            )
    st.caption(
        "Examples fill protein sequences. Any ligand, cofactor, or contact restraint named in an example is not included by this protein-only input."
    )
    uploads = []
    sequences: list[str] = []
    selected_paths: list[str] = []
    selected_database_ids: list[str] = []
    database = None
    database_structures: dict[str, dict] = {}
    if mode == "propagate":
        source = st.pills(
            "Starting structures",
            ["Project designs", "Molecular database", "Upload structures"],
            selection_mode="single",
            default="Project designs",
            key="bioemu2-starting-structure-source",
        ) or "Project designs"
        with st.container(border=True):
            st.subheader("Target–design complexes")
            if source == "Upload structures":
                uploads = st.file_uploader(
                    "Upload PDB or mmCIF structures",
                    type=["pdb", "ent", "cif", "mmcif"],
                    accept_multiple_files=True,
                    key="bioemu2-starting-structures",
                    help="Upload one target–design complex per file. Every file is sampled independently within this single queued job. PDB inputs are converted to mmCIF inside the container.",
                ) or []
                has_structures = bool(uploads)
                if uploads:
                    st.caption(f"{len(uploads)} starting complex(es)")
            elif source == "Molecular database":
                try:
                    from mn_molecular_design_studio.molecular_database import MolecularDatabase

                    database = MolecularDatabase()
                    groups = database.groups(user, scope="all", collection="library", hide_smoke_test=True)
                    groups = [group for group in groups if int(group.get("molecule_count") or 0) > 0]
                except (OSError, ValueError, PermissionError) as exc:
                    groups = []
                    st.warning(f"Could not load molecular database structures: {exc}")
                if not groups:
                    st.info("No accessible molecular database groups are available. Structures saved in the Database library will appear here after synchronization.")
                    has_structures = False
                else:
                    selected_group = st.selectbox(
                        "Database group",
                        groups,
                        format_func=lambda group: f"{group['name']} · {group['molecule_count']} molecules",
                        key="bioemu2-database-group",
                    )
                    try:
                        records = database.list_molecules(
                            user,
                            group_id=selected_group["id"],
                            scope="all",
                            collection="library",
                            sort="name",
                            descending=False,
                            limit=100_000,
                        )["items"]
                    except (OSError, ValueError, PermissionError) as exc:
                        records = []
                        st.warning(f"Could not read this database group: {exc}")
                    for record in records:
                        for structure in record.get("structures") or []:
                            if structure.get("group_id") != selected_group["id"]:
                                continue
                            if str(structure.get("format") or "").lower() not in {"pdb", "ent", "cif", "mmcif"}:
                                continue
                            structure_id = str(structure.get("id") or "")
                            if not structure_id:
                                continue
                            database_structures[structure_id] = {
                                "record": record,
                                "structure": structure,
                            }
                    if not database_structures:
                        st.info("This database group has no saved PDB or mmCIF structures.")
                        has_structures = False
                    else:
                        selected_database_ids = st.multiselect(
                            "Structures to propagate",
                            list(database_structures),
                            format_func=lambda structure_id: (
                                f"{database_structures[structure_id]['record']['name']} · "
                                f"{database_structures[structure_id]['structure']['label']}"
                            ),
                            key="bioemu2-database-structures",
                            help="Select saved target–design complexes. Portable Database structures can be used from any connected workstation.",
                        )
                        has_structures = bool(selected_database_ids)
                        if len(selected_database_ids) > 32:
                            st.error("Select at most 32 starting structures per BioEmu-2 batch.")
                            has_structures = False
                        elif has_structures:
                            st.caption(f"{len(selected_database_ids)} database structure(s) selected")
            else:
                try:
                    sources = source_job_options(project_id) if project_id else []
                except (OSError, ValueError) as exc:
                    sources = []
                    st.warning(f"Could not load completed project structures: {exc}")
                if not sources:
                    st.info("No completed project results with saved complex structures are available.")
                    has_structures = False
                else:
                    selected_source = st.selectbox(
                        "Completed design or prediction job",
                        sources,
                        format_func=lambda row: row["label"],
                        key="bioemu2-source-job",
                    )
                    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
                    source_run = Path(str(selected_source["row"].get("run_dir") or ""))
                    candidates = [
                        candidate for candidate in read_candidates(source_run)
                        if candidate.get("complex_pdb")
                        and Path(str(candidate["complex_pdb"])).is_file()
                        and Path(str(candidate["complex_pdb"])).suffix.lower() in {".pdb", ".ent", ".cif", ".mmcif"}
                    ]
                    labels = {
                        str(candidate.get("candidate_id") or f"Candidate {index + 1}"): candidate
                        for index, candidate in enumerate(candidates)
                    }
                    selected_ids = st.multiselect(
                        "Designs to propagate",
                        list(labels),
                        format_func=lambda candidate_id: f"{candidate_id} · {Path(str(labels[candidate_id]['complex_pdb'])).name}",
                        key="bioemu2-source-candidates",
                        help="Select saved complex structures. Each design is propagated independently in one scheduled batch.",
                    )
                    selected_paths = [str(Path(labels[candidate_id]["complex_pdb"]).resolve()) for candidate_id in selected_ids]
                    has_structures = bool(selected_paths)
                    if not candidates:
                        st.info("This result job has no supported complex-structure files.")
            st.markdown(
                '<div class="bioemu2-note">Project designs lists completed jobs in the selected project. Choose Molecular database for saved designs such as Tau R3 100-Pass Binder Design, or upload target–design complexes. Each selected structure is propagated independently in one batch and produces an XTC ensemble, topology mmCIF, and physicality report.</div>',
                unsafe_allow_html=True,
            )
    else:
        with st.container(border=True):
            st.subheader("Protein chain sequences")
            st.caption("Provide one protein chain per entry. Multiple sequences are treated as chains in the same system.")
            if "bioemu2-component-ids" not in st.session_state:
                _replace_sequence_components([""])
            component_ids = list(st.session_state.setdefault("bioemu2-component-ids", []))
            sequences = []
            for component_index, component_id in enumerate(component_ids):
                sequence_key = _component_key(component_id, "sequence")
                copies_key = _component_key(component_id, "copies")
                st.session_state.setdefault(sequence_key, "")
                st.session_state.setdefault(copies_key, 1)
                with st.container(border=True):
                    heading, copy_col, duplicate_col, remove_col = st.columns([1.4, 1, .8, .7])
                    with heading:
                        st.markdown("**Protein**")
                        st.caption(f"Chain {chr(65 + component_index) if component_index < 26 else component_index + 1}")
                    with copy_col:
                        copies = st.number_input(
                            "Copies", min_value=1, max_value=26, step=1,
                            key=copies_key, label_visibility="visible",
                        )
                    with duplicate_col:
                        st.markdown("<div style='height:1.72rem'></div>", unsafe_allow_html=True)
                        st.button(
                            "Duplicate", key=f"bioemu2-duplicate-{component_id}",
                            on_click=_duplicate_sequence_component, args=(component_id,),
                            use_container_width=True,
                        )
                    with remove_col:
                        st.markdown("<div style='height:1.72rem'></div>", unsafe_allow_html=True)
                        st.button(
                            "Remove", key=f"bioemu2-remove-{component_id}",
                            on_click=_remove_sequence_component, args=(component_id,),
                            use_container_width=True,
                        )
                    sequence = st.text_area(
                        "Protein sequence", height=105, key=sequence_key,
                        placeholder="Paste an amino-acid sequence",
                    )
                    clean_sequence = "".join(sequence.split()).upper()
                    if clean_sequence:
                        st.caption(f"{len(clean_sequence)} residues")
                        sequences.extend([clean_sequence] * int(copies))
            action_cols = st.columns([1.2, 1.2, 1.1, 1.1])
            with action_cols[0]:
                st.button("＋ Add component", key="bioemu2-add-component", on_click=_add_sequence_component, use_container_width=True)
            with action_cols[1]:
                with st.popover("Import FASTA", use_container_width=True):
                    fasta_file = st.file_uploader(
                        "FASTA file", type=["fasta", "fa", "faa", "txt"],
                        key="bioemu2-fasta-upload",
                        help="Upload one or more protein sequences. FASTA headers split chains.",
                    )
                    if fasta_file and st.button("Load sequences", key="bioemu2-load-fasta", use_container_width=True):
                        try:
                            fasta_text = fasta_file.getvalue().decode("utf-8")
                            _replace_sequence_components(_parse_sequences(fasta_text))
                            st.rerun()
                        except UnicodeDecodeError:
                            st.error("The FASTA file must be UTF-8 text.")
            with action_cols[2]:
                st.button("Clear all", key="bioemu2-clear-components", on_click=_replace_sequence_components, args=([],), use_container_width=True)
            if sequences:
                st.caption(f"{len(sequences)} chain(s) · {sum(map(len, sequences)):,} total residues")
                if len(sequences) > 26:
                    st.error("BioEmu-2 accepts at most 26 protein chains per system.")
            else:
                st.caption("No protein chains added yet.")
            if mode == "sample-equilibrium":
                st.info("BioEmu-2 samples the equilibrium ensemble for the supplied sequence-defined protein system.")
            else:
                st.info("BioEmu-2 predicts a folded structure for the supplied sequence-defined protein system.")

    with st.container(border=True):
        st.subheader("Sampling")
        left, middle, right = st.columns(3)
        with left:
            num_samples = st.number_input("Samples per complex", min_value=1, max_value=256, value=8, step=1, key="bioemu2-num-samples")
        with middle:
            lag_time_ns = st.number_input(
                "Time lag (ns)", min_value=0.001, max_value=1000.0, value=1.0, step=1.0,
                disabled=mode != "propagate", key="bioemu2-lag-time",
                help="Propagation endpoints at this lag are model-based samples, not atomistic MD frames.",
            )
        with right:
            sampler = st.selectbox("Sampler", ["edm", "edm_gd", "edm_fk", "edm_fkc"], key="bioemu2-sampler", help="EDM is the unsteered default sampler.")
        with st.expander("Advanced settings", expanded=False):
            col1, col2, col3 = st.columns(3)
            with col1:
                batch_size = st.number_input("GPU batch size", min_value=1, max_value=256, value=8, step=1, key="bioemu2-batch-size")
                seed = st.number_input("Random seed", min_value=0, max_value=2_147_483_647, value=0, step=1, key="bioemu2-seed", help="Zero uses a generated deterministic seed recorded with the job.")
            with col2:
                precision = st.selectbox("Trunk precision", ["bf16-mixed", "32-true"], key="bioemu2-precision")
                matmul_precision = st.selectbox("Denoiser matmul precision", ["high", "highest", "medium"], key="bioemu2-matmul-precision")
            with col3:
                num_steps = st.number_input("Denoising steps (optional)", min_value=0, max_value=2000, value=0, step=10, key="bioemu2-num-steps", help="Zero uses the sampler's configured default.")
                physicality_filter = st.selectbox(
                    "Physicality report profile",
                    ["high_quality", "high_quality_per_residue", "no_violations", "no_filter"],
                    format_func=lambda value: {
                        "high_quality": "High quality",
                        "high_quality_per_residue": "High quality per residue",
                        "no_violations": "No violations",
                        "no_filter": "No filter",
                    }[value],
                    key="bioemu2-physicality",
                    help="Controls the per-frame physicality report. All sampled structures are retained.",
                )
                compile_denoiser = st.checkbox("Compile denoiser", value=False, key="bioemu2-compile", help="Can add substantial startup time; generally useful only for large sample batches.")

    prioritize = st.checkbox("Prioritize in scheduler queue", key="bioemu2-prioritize", help="The scheduler places this GPU job on an available workstation/GPU.")
    has_input = has_structures if mode == "propagate" else bool(sequences)
    with st.container(border=True):
        st.subheader("Estimated Runtime")
        st.caption("Runtime depends on the selected workstation and GPU, input length, sample count, and current queue. Follow queued and running status in My Results.")
    run_columns = st.columns([1, 1, 1])
    with run_columns[1]:
        run = st.button("Submit", type="primary", disabled=not project_id or not has_input, key="bioemu2-run")
    if not project_id:
        st.markdown('<div style="text-align:center;color:#748087;font-size:.82rem;font-style:italic">Select a project to submit.</div>', unsafe_allow_html=True)
    elif not has_input:
        missing_text = "Select at least one starting structure." if mode == "propagate" else "Enter at least one protein sequence."
        st.markdown(f'<div style="text-align:center;color:#748087;font-size:.82rem;font-style:italic">{missing_text}</div>', unsafe_allow_html=True)
    if run:
        try:
            structure_paths = (
                _store_structures(uploads, owner_id) if uploads
                else [str(database.structure_file(user, structure_id, materialize=True)) for structure_id in selected_database_ids] if selected_database_ids and database
                else selected_paths if mode == "propagate"
                else []
            )
            parameters = {
                "mode": mode,
                "structure_paths": structure_paths,
                "protein_sequences": sequences,
                "num_samples": int(num_samples),
                "lag_time_ns": float(lag_time_ns),
                "sampler": sampler,
                "physicality_filter": physicality_filter,
                "precision": precision,
                "matmul_precision": matmul_precision,
                "batch_size": int(batch_size),
                "num_steps": int(num_steps) if num_steps else None,
                "seed": int(seed),
                "compile_denoiser": bool(compile_denoiser),
                "job_name": job_name or "bioemu2-sampling",
                "gpu_device": "all",
            }
            submit(
                schema,
                "simulation.bioemu2",
                parameters,
                project_id,
                {"prioritize": prioritize},
                "BioEmu-2",
            )
        except (OSError, ValueError, TypeError) as exc:
            st.error(f"Could not queue BioEmu-2: {exc}")
