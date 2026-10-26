"""Standalone Kejia IDR inputs, scheduler submission, and pipeline handoff."""
from __future__ import annotations

import hashlib
from pathlib import Path
import shutil

import streamlit as st

from mn_molecular_design_studio import kejia_idr_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.design_validation import protein_chains, read_structure
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path


def _render_template_preview(template: Path, binder: str, peptide: str) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
        ChainVisualization, StructureVisualization, molstar_custom_component,
    )

    chains = protein_chains(read_structure(template))
    content = template.read_text()
    digest = hashlib.sha256(content.encode()).hexdigest()[:16]
    st.subheader("Template complex · MolView")
    st.caption(f"Binder: chain {binder} · {len(chains[binder])} residues · green. "
               f"Template peptide: chain {peptide} · {len(chains[peptide])} residues · orange.")
    molstar_custom_component(
        [StructureVisualization(
            pdb=content,
            representation_type=None,
            chains=[
                ChainVisualization(chain_id=binder, color="uniform", color_params={"value": "0x318c70"},
                                   representation_type="cartoon", label="Binder template"),
                ChainVisualization(chain_id=peptide, color="uniform", color_params={"value": "0xe39c35"},
                                   representation_type="cartoon+ball-and-stick", label="Template peptide (target position)"),
            ],
        )],
        key=f"kejia-template-view-{digest}-{binder}-{peptide}",
        height=520,
        show_controls=True,
        show_structure_tools=False,
    )
    st.caption("Shows the template before threading. Your sequence fragment is applied when the job runs.")


def _open_pipeline(parameters: dict, project: str, *, run_dir: Path | None = None, candidates=None) -> None:
    from mn_molecular_design_studio.pipeline_graph import new_graph
    graph = new_graph("kejia_idr_generated" if run_dir else "kejia_idr")
    if run_dir:
        graph["nodes"][0]["data"]["params"].update(
            source_run_dir=str(run_dir), selected_candidate_ids=[c["candidate_id"] for c in candidates or []])
    else:
        allowed = set(graph["nodes"][0]["data"]["params"])
        graph["nodes"][0]["data"]["params"].update({key: value if value is not None else "" for key, value in parameters.items() if key in allowed})
    st.session_state["pipeline-pending-draft"] = {"project_id": project, "graph": graph}
    st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = project
    st.session_state["page"] = "Pipelines"
    st.rerun()


def render_setup(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    if st.button("← Back to All tools", key="kejia-back"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    st.title(inputs.TITLE)
    st.caption("Template threading and flexible peptide docking for intrinsically disordered targets")
    st.markdown("Enter a sequence fragment and a binder–peptide template. Studio runs native [Kejia threading](https://github.com/drhicks/Kejia_peptide_binders) in a dedicated CPU container. Use the pipeline canvas to add ProteinMPNN and AlphaFold2 Initial Guess.")
    with st.expander("About the method"):
        st.write("The peptide can adopt different conformations across templates and flexible docking. A template is still required; the workflow does not predict an unbound IDR structure. A longer input than the template peptide is evaluated as overlapping fragments, with the matched range recorded for each candidate.")
        st.write("Walle and Walle2 are the authors' recommended starting libraries. Mini and Strand mainly represent strand-pairing geometries. The Studio starter uses existing ProteinMPNN and AF2 adapters; the upstream iterative filtering and modified diffusion branch are not reproduced automatically.")
        st.code(inputs.CODE_VERSION, language=None)

    left, right = st.columns(2)
    with left:
        job_name = st.text_input("Job name", "kejia-idr", key="kejia-job-name")
    with right:
        project = project_selector("kejia-idr")
    if st.button("Load Tamarind example inputs", key="kejia-example"):
        st.session_state["kejia-sequence"] = inputs.EXAMPLE_SEQUENCE
        st.session_state["kejia-source"] = "upload"
        st.session_state["kejia-example-template"] = str(inputs.example_template())
        st.session_state["kejia-upload-revision"] = st.session_state.get("kejia-upload-revision", 0) + 1
        st.session_state["kejia-binder-chain"] = "A"
        st.session_state["kejia-peptide-chain"] = "B"
        st.rerun()
    sequence = st.text_area("Sequence fragment", key="kejia-sequence", height=120,
                            help="One target fragment or one FASTA record, 4–40 residues. Whitespace is ignored.")
    if sequence:
        try:
            st.caption(f"{len(inputs.normalize_sequence(sequence))} residues")
        except ValueError as exc:
            st.error(str(exc))
    counts = inputs.library_counts()
    library = st.selectbox("Template source", ["upload", *inputs.LIBRARIES], key="kejia-source",
                           format_func=lambda value: "Upload a template" if value == "upload" else f"{value.title()} library · {counts[value]} templates")
    template, binder, peptide = None, "A", "B"
    if library == "upload":
        upload = st.file_uploader("Template file", type=["pdb", "ent", "cif", "mmcif"], key=f"kejia-upload-{st.session_state.get('kejia-upload-revision', 0)}",
                                  help="Exactly two protein chains: binder and peptide.")
        if upload:
            data = upload.getvalue()
            template = user_input_root(owner_id) / "kejia-idr" / hashlib.sha256(data).hexdigest()[:16] / Path(upload.name).name
            template.parent.mkdir(parents=True, exist_ok=True)
            if not template.exists():
                template.write_bytes(data)
        elif st.session_state.get("kejia-example-template"):
            source = Path(st.session_state["kejia-example-template"])
            if source.is_file():
                template = user_input_root(owner_id) / "kejia-idr/example" / source.name
                template.parent.mkdir(parents=True, exist_ok=True)
                if not template.exists():
                    shutil.copy2(source, template)
                st.caption("Tamarind example: FNCREWCWN with Walle conf_2.pdb")
            else:
                st.error("The example template is missing. Install the Kejia library first.")
        if template:
            try:
                chains = inputs.template_chains(template)
                a, b = st.columns(2)
                binder = a.selectbox("Binder chain", chains, index=0, key="kejia-binder-chain")
                peptide = b.selectbox("Peptide chain", chains, index=1, key="kejia-peptide-chain")
                if binder != peptide:
                    _render_template_preview(template, binder, peptide)
            except (ValueError, OSError) as exc:
                st.error(str(exc))
    elif library in {"mini", "strand"}:
        st.caption("This library emphasizes strand pairing. The authors recommend starting with Walle or Walle2 unless the target supports that binding geometry.")
    with st.expander("Optional settings"):
        a, b = st.columns(2)
        maximum = int(a.number_input("Maximum library templates", 1, 1000, 10, disabled=library == "upload"))
        repeats = int(b.number_input("Flexible docking repeats per register", 1, 100, 3))
        identity = a.number_input("Minimum sequence identity", 0.0, 1.0, 0.0)
        blosum = b.number_input("Minimum BLOSUM score", -100.0, 100.0, -10.0)
        farep = a.number_input("Maximum repulsive-energy increase", 0.0, 10000.0, 10.0)
        keep = b.checkbox("Retain template bidentate contacts", False)
        seed = int(a.number_input("Random seed (0 = random)", 0, 2147483647, 0))
        cpus = int(b.number_input("CPU cores", 1, 256, 4))
    parameters = dict(target_sequence=sequence, template_library=library, template_structure=str(template) if template else None,
                      binder_chain=binder, peptide_chain=peptide, max_templates=maximum,
                      flex_docking_repeats=repeats, min_seq_identity=float(identity), min_blosum=float(blosum),
                      max_farep=float(farep), keep_bidentates=keep, seed=seed, cpu_cores=cpus, job_name=job_name or "kejia-idr")
    error = None
    try:
        inputs.prepare_request(parameters)
    except (TypeError, ValueError, OSError) as exc:
        error = str(exc)
    a, b = st.columns(2)
    if a.button("Run template threading", type="primary", disabled=bool(error) or not project, key="kejia-submit"):
        submit(schema, schema["workflow"], parameters, project, {}, inputs.TITLE)
    if b.button("Open full pipeline", disabled=bool(error) or not project, key="kejia-open-pipeline"):
        _open_pipeline(parameters, project)
    if error:
        st.caption(error)
    if not project:
        st.caption("Select a project to submit or open a pipeline.")


def render_results(run_id, run_dir, candidates, parameters, result):
    st.subheader("Kejia IDR · threaded structures")
    st.metric("Threaded candidates", len(candidates))
    st.caption("These structures need binder sequence redesign and independent prediction before interpreting binding confidence.")
    project = read_json(Path(run_dir) / "molecular_design_studio.json").get("project_id")
    if candidates and result.get("success") is True and project:
        if st.button("Continue with ProteinMPNN → AF2", type="primary", key="kejia-continue-" + run_id):
            _open_pipeline({}, project, run_dir=Path(run_dir), candidates=candidates)
    if not candidates:
        st.info("No threaded candidates are available. Check the job logs and template filters.")
        return
    st.dataframe([{
        "Candidate": row["candidate_id"],
        "Target range": "–".join(map(str, row.get("raw_metadata", {}).get("target_fragment_range", []))),
        "Repulsive-energy increase": row.get("metrics", {}).get("delta_fa_rep"),
        "Template": Path(row.get("raw_metadata", {}).get("template", {}).get("path", "")).name,
    } for row in candidates], hide_index=True, use_container_width=True)
    selected = st.selectbox("Threaded candidate", range(len(candidates)), format_func=lambda index: candidates[index]["candidate_id"], key="kejia-result-" + run_id)
    row = candidates[selected]
    path = resolve_stored_path(row.get("complex_pdb"), run_dir=Path(run_dir), must_exist=True)
    if path:
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component, StructureVisualization
        molstar_custom_component([StructureVisualization(pdb=path.read_text(), color="chain-id", representation_type="cartoon")],
                                 key=f"kejia-view-{run_id}-{selected}", height=520, show_structure_tools=False)
        st.download_button("Download threaded complex", path.read_bytes(), file_name=path.name, key=f"kejia-download-{run_id}-{selected}")
    with st.expander("Threading details"):
        st.json(row.get("raw_metadata") or {})
