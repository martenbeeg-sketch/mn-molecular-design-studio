"""Standalone SAbR input page for antibody numbering and CDR annotation."""
from __future__ import annotations

from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import mpnn_input
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.target_editor import TargetSource, render_target_source
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
    StructureVisualization,
    molstar_custom_component,
)


def _load_example(owner_id: str) -> None:
    source = mpnn_input.ASSETS / "9F6O.pdb"
    destination = user_input_root(owner_id) / "sabr-examples" / source.name
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(source.read_bytes())
    st.session_state["sabr-example-path"] = str(destination)


def render_setup(schema: dict, *, owner_id: str, project_selector, submit, search_jobs) -> None:
    st.markdown("""<style>
    .st-key-sabr-surface{max-width:1100px;margin:auto;color:#233e31}
    .st-key-sabr-surface h1{color:#183d30;letter-spacing:-.035em}
    .st-key-sabr-surface [data-testid=stButton] button{border:1px solid #b8d0c0;border-radius:7px;color:#29644c;background:white}
    .st-key-sabr-surface [data-testid=stButton] button[kind=primary]{background:#27674e;color:white;border-color:#27674e}
    .st-key-sabr-surface [data-testid=stFileUploaderDropzone]{border:1px dashed #b8cfbf;background:#f6faf7}
    </style>""", unsafe_allow_html=True)
    if st.button("← Back to All tools", key="sabr-back"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()

    with st.container(key="sabr-surface"):
        st.caption("ANTIBODY NUMBERING · CDR ANNOTATION")
        st.title("SAbR")
        st.write("Number antibody chains and map IMGT CDR1, CDR2, and CDR3 positions onto the original structure.")
        st.link_button("SAbR source code", "https://github.com/delalamo/sabr")
        with st.expander("About SAbR"):
            st.write("SAbR uses the antibody backbone to assign residue numbering. CDR selections are mapped back to the input residue IDs, so downstream design tools can keep the framework fixed and act only on chosen loops. The numbering pass does not change the source file.")

        project = project_selector("sabr")
        mode = st.segmented_control("Input source", ["Structure", "Existing results"], default="Structure", key="sabr-input-mode") or "Structure"
        parameters: dict = {
            "chain_type": "auto", "scheme": "imgt", "mode": "sabr",
            "noise_level": 0.0, "write_numbered_structures": True,
        }
        error = None
        if mode == "Structure":
            st.button("Load 9F6O nanobody example", key="sabr-load-example", on_click=_load_example, args=(owner_id,))
            source = render_target_source(
                key_prefix="sabr-source", label="Antibody structure", owner_id=owner_id,
                example_code="9F6O", allow_cif=True, allow_non_protein=True,
                upload_in_expander=False, upload_label="Upload PDB or mmCIF",
                upload_file_types=["pdb", "ent", "cif", "mmcif"],
            )
            if not source.path and st.session_state.get("sabr-example-path"):
                example_path = Path(st.session_state["sabr-example-path"])
                if example_path.is_file():
                    source = TargetSource(path=example_path, name=example_path.stem, owner_id=owner_id)
            if source.path:
                parameters["structure_path"] = str(source.path)
                inventory = mpnn_input.inventory(source.path)
                chains = [chain for chain, data in inventory["chains"].items() if data["kind"] == "protein"]
                default = ["B"] if "B" in chains else chains[:1]
                selected = st.multiselect("Antibody chain(s)", chains, default=default, key="sabr-chain-selection")
                parameters["chain_ids"] = selected
                if selected:
                    st.caption("Only selected protein chains are numbered. Other chains remain untouched.")
                    st.markdown("**Structure preview**")
                    molstar_custom_component(
                        [StructureVisualization(pdb=source.path.read_text(errors="replace"), color="chain-id",
                            representation_type="cartoon")], key="sabr-structure-preview", height=390,
                        show_controls=True, show_structure_tools=False,
                    )
                else:
                    error = "Select at least one antibody chain."
            else:
                error = "Upload or fetch a PDB/mmCIF structure."
        else:
            from mn_molecular_design_studio.mpnn_ui import _existing_results
            source, error = _existing_results(project, search_jobs, "sabr", antibody_backbones_only=True)
            parameters.update(source)
            if not project:
                error = "Select a project."

        st.markdown("### Numbering and loop settings")
        a, b, c = st.columns(3)
        parameters["scheme"] = a.selectbox("Numbering scheme", ["imgt", "chothia", "kabat", "martin", "aho", "wolfguy"], index=0, key="sabr-scheme")
        parameters["chain_type"] = b.selectbox("Chain type", ["auto", "H", "K", "L"], index=0, key="sabr-chain-type")
        parameters["mode"] = c.selectbox("Alignment mode", ["sabr", "softalign"], index=0, key="sabr-mode")
        parameters["noise_level"] = st.select_slider("Backbone noise level", options=[0.0, 0.2, 0.5, 1.0, 2.0], value=0.0, key="sabr-noise",
            help="Use the reference set matching the backbone-noise level.")
        parameters["write_numbered_structures"] = st.checkbox("Save a numbered structure copy", value=True, key="sabr-write-numbered")
        name = st.text_input("Job name", "sabr-antibody-numbering", key="sabr-job-name").strip() or "sabr-antibody-numbering"
        parameters["job_name"] = name
        if error:
            st.caption(error)
        if st.button("Number and identify CDRs", type="primary", key="sabr-submit", disabled=not project or bool(error)):
            allowed = {row["name"] for row in schema["parameters"]}
            submit(schema, schema["workflow"], {key: value for key, value in parameters.items() if key in allowed},
                   project, {"gpu_device": "none"}, name)

