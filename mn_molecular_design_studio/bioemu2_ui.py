"""Studio input page for BioEmu-2 ensemble sampling."""
from __future__ import annotations

import hashlib
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.auth import user_input_root


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
    sequences: list[str] = []
    current: list[str] = []
    for line in str(text or "").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
            continue
        current.append("".join(stripped.split()).upper())
    if current:
        sequences.append("".join(current))
    return sequences


def render_setup(schema: dict, *, owner_id: str, project_selector, submit, source_job_options) -> None:
    if st.button("← Back to All tools", key="back-to-bioemu2-tools"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.session_state.pop("selected_workflow_group", None)
        st.rerun()

    st.html("""<style>
      .block-container{max-width:1260px!important;padding-top:1.4rem!important}
      .bioemu2-heading{text-align:center;margin:.3rem auto 1.5rem}
      .bioemu2-heading h1{font-size:2.05rem;color:#173d31;margin:0 0 .4rem}
      .bioemu2-heading p{color:#74857b;max-width:700px;margin:.4rem auto;line-height:1.55}
      .bioemu2-card{padding:1rem 1.2rem;border:1px solid #d8e5dc;border-radius:10px;background:#fff;margin:.5rem 0}
      .bioemu2-note{padding:.8rem 1rem;border-radius:8px;background:#edf5f0;color:#315d47}
    </style>""")
    st.markdown(
        '<div class="bioemu2-heading"><h1>BioEmu-2</h1>'
        '<p>Sample protein conformational ensembles and time-lagged changes with BioEmu-2.</p></div>',
        unsafe_allow_html=True,
    )
    st.caption("BioEmu-2 ensemble sampling · 0.1.0rc0")
    name_col, project_col = st.columns([1.2, 1.2])
    with name_col:
        job_name = st.text_input("Job name", value="bioemu2-sampling", key="bioemu2-job-name").strip()
    with project_col:
        project_id = project_selector("bioemu2-project")

    mode_label = st.radio(
        "Input mode",
        ["Propagate starting complex", "Sample protein sequence(s)", "Predict structure from sequence(s)"],
        horizontal=True,
        key="bioemu2-input-mode",
        help="Propagation starts from one or more predicted or experimental complexes. Sequence modes define one system; multiple sequences are chains of the same complex.",
    )
    mode = {
        "Propagate starting complex": "propagate",
        "Sample protein sequence(s)": "sample-equilibrium",
        "Predict structure from sequence(s)": "predict-structure",
    }[mode_label]
    uploads = []
    sequences: list[str] = []
    selected_paths: list[str] = []
    if mode == "propagate":
        source = st.radio(
            "Starting structures",
            ["Upload structures", "Choose from completed project results"],
            horizontal=True,
            key="bioemu2-structure-source",
        )
        with st.container(border=True):
            st.subheader("Starting complexes")
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
                        "Completed result job",
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
                        "Target–design complexes",
                        list(labels),
                        format_func=lambda candidate_id: f"{candidate_id} · {Path(str(labels[candidate_id]['complex_pdb'])).name}",
                        key="bioemu2-source-candidates",
                        help="Each selected candidate is propagated independently in the same scheduled job.",
                    )
                    selected_paths = [str(Path(labels[candidate_id]["complex_pdb"]).resolve()) for candidate_id in selected_ids]
                    has_structures = bool(selected_paths)
                    if not candidates:
                        st.info("This result job has no supported complex-structure files.")
            st.markdown(
                '<div class="bioemu2-note">For a target with several designed binders, upload each predicted target–binder complex as a separate structure. BioEmu-2 propagates each complex independently and writes an XTC ensemble, topology mmCIF, and physicality report.</div>',
                unsafe_allow_html=True,
            )
    else:
        with st.container(border=True):
            st.subheader("Protein chain sequences")
            sequence_text = st.text_area(
                "Enter one protein chain per line (FASTA headers are optional)",
                height=150,
                key="bioemu2-sequences",
                placeholder="> target\nMST...\n> designed binder\nEVQ...",
            )
            sequences = _parse_sequences(sequence_text)
            if sequences:
                st.caption(f"{len(sequences)} chain(s) · {sum(map(len, sequences)):,} total residues")
            st.info("All entered sequences form one complex. To sample separate target–design predictions, use Propagate starting complex and upload one predicted complex per design.")

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
    if not project_id:
        st.info("Select a project to enable scheduling.")
    run = st.button("Run BioEmu-2", type="primary", disabled=not project_id or not has_input, key="bioemu2-run")
    if run:
        try:
            structure_paths = (
                _store_structures(uploads, owner_id) if uploads
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
