"""Standalone setup and result pages for β-pairing conditioned RFdiffusion."""
from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import rfdiffusion_beta_pairing_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import StructureVisualization, molstar_custom_component
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.target_editor import _target_widget_keys, render_target_editor, render_target_source


STYLE = """<style>
.st-key-beta-pairing-surface{--beta-green:#27674e;--beta-ink:#183e31;color:#193e30}
.st-key-beta-pairing-surface .beta-heading{text-align:center;padding:.5rem 0 1.3rem}
.st-key-beta-pairing-surface .beta-heading h1{font-size:2.1rem;color:#173d31;margin:.2rem 0 .5rem;letter-spacing:-.04em}
.st-key-beta-pairing-surface .beta-heading p{color:#71867a;font-size:.92rem;line-height:1.6;max-width:700px;margin:.4rem auto}
.st-key-beta-pairing-surface .beta-eyebrow,.st-key-beta-pairing-surface .beta-section{font-size:.7rem;font-weight:650;letter-spacing:.11em;color:#658572;text-transform:uppercase}
.st-key-beta-pairing-surface .beta-section{margin:.8rem 0 .6rem}
.st-key-beta-pairing-surface [data-testid=stVerticalBlockBorderWrapper]{border-color:#d7e6dc!important;border-radius:10px!important}
.st-key-beta-pairing-surface [data-testid=stTextInput] input,.st-key-beta-pairing-surface [data-testid=stNumberInput] input,.st-key-beta-pairing-surface textarea{background:#f3f7f4!important;color:#193e30!important}
.st-key-beta-pairing-surface [data-baseweb=input],.st-key-beta-pairing-surface [data-baseweb=textarea],.st-key-beta-pairing-surface [data-baseweb=select]>div{border-color:#d5e2d8!important;background:#f3f7f4!important;border-radius:6px!important}
.st-key-beta-pairing-surface [data-baseweb=tag]{background:#e0eee5!important;color:#27674e!important}
.st-key-beta-pairing-surface [data-testid=stButton] button{border-radius:6px!important;border-color:#b7d0c0!important;color:#27674e;box-shadow:none!important}
.st-key-beta-pairing-surface [data-testid=stButton] button:hover{border-color:#69977c!important;background:#f0f7f2!important;color:#174d32!important}
.st-key-beta-pairing-surface button[kind=primary]{background:#27674e!important;border-color:#27674e!important;color:#fff!important}
.st-key-beta-pairing-surface [data-testid=stExpander]{border-color:#d6e3da!important;border-radius:7px!important;background:#fcfefc}
.st-key-beta-pairing-surface [data-testid=stFileUploaderDropzone]{border:1px dashed #bad1c1!important;background:#f5f9f6!important;border-radius:8px}
</style>"""


def render_setup(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    st.html(STYLE)
    with st.container(key="beta-pairing-surface"):
        if st.button("← Back to All tools", key="beta-pairing-back"):
            st.session_state["page"] = "All tools"
            st.session_state.pop("selected_workflow", None)
            st.rerun()
        st.html('<div class="beta-heading"><div class="beta-eyebrow">BETA-PAIRING BINDER GENERATION</div><h1>RFdiffusion β-pairing</h1><p>Condition RFdiffusion with explicit secondary-structure and adjacency tensors so a designed binder β-strand pairs with an exposed target edge strand.</p></div>')
        with st.expander("How this method works", expanded=False):
            st.write("Studio builds the published SS/adjacency interface tensors for the selected target β-strand(s), then samples RFdiffusion backbones across the possible binder-strand registers. This is a backbone generator: run ProteinMPNN for sequences, then use AlphaFold2 Initial Guess or another structure-prediction step to evaluate them.")
            st.markdown("[Paper: Improved protein binder design using β-pairing targeted RFdiffusion](https://doi.org/10.1038/s41467-025-67866-3) · [RFdiffusion source](https://github.com/RosettaCommons/RFdiffusion) · [Published tensor script](https://files.ipd.uw.edu/pub/strand_pairing_binders_2025/interface_tensors.zip)")
            st.caption("Studio uses its existing RFdiffusion and PyRosetta scoring images. The cited AF2 initial-guess implementation is available as a downstream Studio step; this tile does not duplicate or bundle AF2 model weights.")

        left, right = st.columns([1.5, 1], vertical_alignment="bottom")
        with left:
            job_name = st.text_input("Job name", "rfdiffusion-beta-pairing", key="beta-pairing-job-name").strip()
        with right:
            project = project_selector("beta-pairing")

        st.html('<div class="beta-section">TARGET AND PAIRED EDGE STRAND</div>')
        prefix = "beta-pairing-target"
        example = st.button("Load ALK2 paper example", icon=":material/science:", key="beta-pairing-example")
        if example:
            source = inputs.TENSOR_SCRIPT.parent / "alk2_target.pdb"
            data = source.read_bytes()
            destination = user_input_root(owner_id) / "beta-pairing-examples" / hashlib.sha256(data).hexdigest()[:12] / source.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                shutil.copy2(source, destination)
            st.session_state[prefix + "-fetched-target"] = str(destination)
            st.session_state[prefix + "-fetched-target-name"] = "ALK2 paper example"
            st.session_state[_target_widget_keys(prefix, destination)["hotspots"]] = "B99-103"
            st.rerun()
        source = render_target_source(
            key_prefix=prefix,
            label="Target",
            example_code="2E9W",
            owner_id=owner_id,
            allow_cif=False,
            upload_file_types=["pdb", "ent"],
            upload_help="Choose one protein-chain PDB. The selected chain is renumbered contiguously in the job copy so gaps and insertion codes cannot shift the β-pairing tensors.",
            upload_label="Upload a target PDB",
        )
        selected = None
        error: str | None = None
        if source.path:
            try:
                key = _target_widget_keys(prefix, source.path)["hotspots"]
                if source.path.name == "alk2_target.pdb":
                    st.session_state.setdefault(key, "B99-103")
                selected = render_target_editor(
                    source,
                    key_prefix=prefix,
                    label="Target",
                    residue_selection_title="Exposed target β-strand residues",
                    residue_selection_caption="Select one or more continuous edge-strand segments in the sequence panel or MolView. Each selected segment gets a corresponding designed binder β-strand.",
                    residue_selection_input_label="Target β-strand residue(s)",
                    allow_trimming=True,
                    viewer_height=470,
                    show_structure_tools=False,
                )
                if len(selected.chains) != 1:
                    error = "Choose exactly one protein chain for the β-pairing target."
                elif selected.hotspots:
                    numbers = inputs.parse_residue_selection(selected.hotspots, selected.chains[0])
                    groups = []
                    for number in numbers:
                        if not groups or number != groups[-1][-1] + 1:
                            groups.append([])
                        groups[-1].append(number)
                    if any(len(group) < 3 for group in groups):
                        error = "Each target β-strand segment must contain at least three consecutive residues."
                    else:
                        st.caption(f"Selected {len(groups)} target strand segment(s). Studio will enumerate about {inputs.estimate_register_count(65, 6, len(groups)):,} binder-strand register(s) for the default 65-residue binder.")
            except (ValueError, OSError, KeyError) as exc:
                error = str(exc)

        with st.container(border=True):
            st.html('<div class="beta-section">BACKBONE SAMPLING</div>')
            col1, col2 = st.columns(2)
            with col1:
                binder_length = int(st.number_input("Binder length", 30, 500, 65, key="beta-pairing-binder-length"))
                binder_strand_length = int(st.number_input("Designed β-strand length", 3, 40, 6, key="beta-pairing-strand-length", help="The paper’s ALK2 example conditions a six-residue binder β-strand."))
            with col2:
                num_designs = int(st.number_input("Number of backbones", 1, 10000, 10, key="beta-pairing-num-designs"))
            st.caption("The sampler selects RFdiffusion’s Complex_Fold_base checkpoint, which is trained for secondary-structure and block-adjacency conditioning. The tensor script enumerates all non-overlapping binder-strand placements. Random mode samples across those registers; systematic mode walks them in order. Set the backbone count to at least the register count to cover every register once.")

        with st.expander("Optional settings"):
            col1, col2 = st.columns(2)
            with col1:
                timesteps = int(st.number_input("Diffusion timesteps", 1, 200, 50, key="beta-pairing-timesteps"))
                systematic = st.checkbox("Systematic register sampling", value=False, key="beta-pairing-systematic")
            with col2:
                seed = int(st.number_input("Random seed (0 = random)", 0, 2147483647, 0, key="beta-pairing-seed"))
                cpu_cores = int(st.number_input("CPU cores for tensor generation", 1, 256, 4, key="beta-pairing-cpus"))

        if error:
            st.error(error)
        can_submit = bool(project and selected and len(selected.chains) == 1 and selected.hotspots and not error)
        if st.button("Generate β-pairing backbones", key="beta-pairing-submit", type="primary", use_container_width=True, disabled=not can_submit):
            params = {
                "target_pdb": str(selected.design_path or source.path),
                "target_chain": selected.chains[0],
                "target_beta_strands": selected.hotspots,
                "binder_length": binder_length,
                "binder_strand_length": binder_strand_length,
                "num_designs": num_designs,
                "timesteps": timesteps,
                "systematic_sampling": systematic,
                "seed": seed,
                "cpu_cores": cpu_cores,
                "gpu_device": "0",
                "job_name": job_name or "rfdiffusion-beta-pairing",
            }
            try:
                inputs.prepare_request(params)
                submit(schema, schema["workflow"], params, project, {"prioritize": False}, inputs.TITLE)
            except (TypeError, ValueError, OSError) as exc:
                st.error(f"Could not start β-pairing generation: {exc}")
        if not project:
            st.caption("Select a project to submit.")


def render_results(run_id, run_dir, candidates, parameters, result):
    st.subheader("RFdiffusion β-pairing · generated structures")
    st.metric("Generated backbones", len(candidates))
    st.caption("These are SS/adjacency-conditioned backbones. They do not contain designed binder sequences or prediction scores yet.")
    marker = read_json(Path(run_dir) / "molecular_design_studio.json")
    if candidates and result.get("success") is True and marker.get("project_id"):
        if st.button("Continue on pipeline canvas", key="beta-pairing-continue-" + run_id, type="primary"):
            from mn_molecular_design_studio.pipeline_graph import new_graph
            graph = new_graph("rfdiffusion_beta_pairing_generated")
            graph["nodes"][0]["data"]["params"].update(
                source_run_dir=str(run_dir),
                selected_candidate_ids=[candidate["candidate_id"] for candidate in candidates],
            )
            # The cited method evaluates redesigned sequences with an initial guess.
            for node in graph["nodes"]:
                if node["data"]["kind"] == "af2":
                    node["data"]["params"].update(
                        use_binder_template=False,
                        use_interface_template=False,
                        msa_source="single_sequence",
                        model_count=1,
                        num_recycles=3,
                    )
            project = marker["project_id"]
            st.session_state["pipeline-pending-draft"] = {"project_id": project, "graph": graph}
            st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = project
            st.session_state["page"] = "Pipelines"
            st.rerun()
    if not candidates:
        st.info("Generated structures will appear here when the job finishes.")
        return
    st.dataframe([
        {
            "Design": candidate["candidate_id"],
            "Binder chain": ", ".join(candidate.get("binder_chains") or []),
            "Target chain": ", ".join(candidate.get("target_chains") or []),
            "Conditioned target strand": candidate.get("raw_metadata", {}).get("target_beta_strands", ""),
            "Stage": "Backbone generated",
        }
        for candidate in candidates
    ], hide_index=True, use_container_width=True)
    selected_index = st.selectbox("Design", range(len(candidates)), format_func=lambda index: candidates[index]["candidate_id"], key="beta-pairing-result-" + run_id)
    candidate = candidates[selected_index]
    path = resolve_stored_path(candidate.get("complex_pdb"), run_dir=Path(run_dir), must_exist=True)
    if path:
        molstar_custom_component(
            [StructureVisualization(pdb=path.read_text(), color="chain-id", representation_type="cartoon")],
            key=f"beta-pairing-view-{run_id}-{selected_index}", height=520, show_structure_tools=False,
        )
        st.download_button("Download generated complex", path.read_bytes(), file_name=path.name, key=f"beta-pairing-download-{run_id}-{selected_index}")
    with st.expander("Generation details"):
        st.json(candidate.get("raw_metadata") or {})
