"""Published peptide-trap scaffold inputs and editable pipeline handoff."""
from __future__ import annotations

from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import peptide_trap_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.kejia_idr_ui import _render_template_preview


def _open_pipeline(parameters: dict, project: str, *, run_dir: Path | None = None, candidates=None) -> None:
    from mn_molecular_design_studio.pipeline_graph import new_graph
    graph = new_graph("peptide_trap_generated" if run_dir else "peptide_trap")
    if run_dir:
        graph["nodes"][0]["data"]["params"].update(source_run_dir=str(run_dir),
            selected_candidate_ids=[row["candidate_id"] for row in candidates or []])
    else:
        allowed = set(graph["nodes"][0]["data"]["params"])
        graph["nodes"][0]["data"]["params"].update({key: value for key, value in parameters.items() if key in allowed})
    st.session_state["pipeline-pending-draft"] = {"project_id": project, "graph": graph}
    st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = project
    st.session_state["page"] = "Pipelines"
    st.rerun()


def render_setup(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    if st.button("← Back to All tools", key="peptide-trap-back"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    st.title(inputs.TITLE)
    st.caption("Retarget published β-strand peptide-binding clefts")
    st.markdown("Based on [Design of amyloidogenic peptide traps](https://doi.org/10.1038/s41589-024-01578-5). Matches a target fragment to the paper's buried/exposed residue pattern, threads it into a published scaffold, then redesigns the binder interface with Rosetta.")
    with st.expander("Method and available templates"):
        st.write("This method assumes the target binds as an extended β-strand. The initial release uses C37 and the full C104.1 design model. The other deposited models are retained for inspection; the initial scaffold-building procedure and hairpin retargeting are not included.")
        st.write("The target sequence stays fixed. Binder-interface redesign and relaxation follow an adaptation of the deposited Rosetta protocol. ProteinMPNN is optional; the starter passes the designed sequence directly to AlphaFold2 and Rosetta interface analysis.")
        st.markdown("[Authors' code and model archive](https://doi.org/10.5281/zenodo.10391229)")
    a, b = st.columns(2)
    job_name = a.text_input("Job name", "peptide-trap", key="peptide-trap-job-name")
    with b:
        project = project_selector("peptide-trap")
    if st.button("Load C37 paper example", key="peptide-trap-example"):
        st.session_state["peptide-trap-sequence"] = inputs.EXAMPLE_SEQUENCE
        st.session_state["peptide-trap-scaffolds"] = ["C37"]
        st.rerun()
    target = st.text_area("Target sequence / single-record FASTA", height=120, key="peptide-trap-sequence",
                          help="Matching windows are eight or ten residues long. Longer sequences are scanned for compatible regions.")
    selected = st.multiselect("Scaffold templates", list(inputs.SCAFFOLDS), default=list(inputs.SCAFFOLDS),
                              format_func=lambda key: inputs.SCAFFOLDS[key]["label"], key="peptide-trap-scaffolds")
    error, prepared = None, None
    try:
        prepared = inputs.prepare_request({"target_sequence": target, "scaffold_ids": selected})
    except (ValueError, TypeError, OSError) as exc:
        error = str(exc)
    if selected:
        preview = st.selectbox("Preview template", selected, format_func=lambda key: inputs.SCAFFOLDS[key]["label"], key="peptide-trap-preview")
        path = inputs.template_path(preview)
        if path.is_file():
            _render_template_preview(path, "B", "A")
            st.caption("Burial pattern: " + inputs.SCAFFOLDS[preview]["pattern"] + " · 0 = buried, 1 = exposed")
        else:
            st.info("The template library is not installed on this workstation.")
    if prepared:
        st.dataframe([{"Scaffold": row["scaffold"], "Start": row["start"], "End": row["end"], "Fragment": row["sequence"]}
                      for row in prepared["all_matches"][:100]], hide_index=True, use_container_width=True)
        st.caption(f"{len(prepared['all_matches'])} compatible placement(s); showing up to 100. A matching pattern is an initial geometry screen, not a binding prediction.")
    with st.expander("Optional settings"):
        a, b = st.columns(2)
        maximum = int(a.number_input("Maximum matching placements", 1, 1000, 20))
        designs = int(b.number_input("Designs per placement", 1, 100, 1))
        design_repeats = int(a.number_input("Rosetta design repeats", 1, 5, 2))
        relax_repeats = int(b.number_input("Rosetta relaxation repeats", 1, 5, 2))
        rmsd = float(a.number_input("Maximum template backbone RMSD (Å)", 0.1, 10.0, 2.0))
        seed = int(b.number_input("Random seed (0 = random)", 0, 2147483647, 0))
        cpus = int(a.number_input("CPU cores", 1, 256, 4))
    parameters = dict(target_sequence=target, scaffold_ids=selected, max_registers=maximum, designs_per_register=designs,
                      design_repeats=design_repeats, relax_repeats=relax_repeats,
                      max_backbone_rmsd=rmsd, seed=seed, cpu_cores=cpus, job_name=job_name or "peptide-trap")
    a, b = st.columns(2)
    if a.button("Generate peptide traps", type="primary", disabled=bool(error) or not project, key="peptide-trap-submit"):
        submit(schema, schema["workflow"], parameters, project, {}, inputs.TITLE)
    if b.button("Open full pipeline", disabled=bool(error) or not project, key="peptide-trap-pipeline"):
        _open_pipeline(parameters, project)
    if error:
        st.caption(error)
    if not project:
        st.caption("Select a project to submit or open a pipeline.")


def render_results(run_id, run_dir, candidates, parameters, result):
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    st.subheader("Peptide Trap · redesigned complexes")
    st.metric("Candidates retaining template geometry", len(candidates))
    st.caption("The binder sequences have been redesigned by Rosetta. Independent prediction and binding validation are still needed.")
    project = read_json(Path(run_dir) / "molecular_design_studio.json").get("project_id")
    if candidates and result.get("success") is True and project:
        if st.button("Continue with AF2 → interface analysis", type="primary", key="peptide-trap-continue-" + run_id):
            _open_pipeline({}, project, run_dir=Path(run_dir), candidates=candidates)
    if not candidates:
        st.info("No candidates passed the template-geometry check. Inspect the job logs and retained raw structures.")
        return
    st.dataframe([{"Candidate": row["candidate_id"], "Scaffold": row.get("raw_metadata", {}).get("scaffold"),
                   "Target fragment": row.get("raw_metadata", {}).get("threaded_target_sequence"),
                   "Rosetta score": row.get("metrics", {}).get("rosetta_total_score"),
                   "Complex backbone RMSD": row.get("metrics", {}).get("template_complex_backbone_rmsd")}
                  for row in candidates], hide_index=True, use_container_width=True)
    index = st.selectbox("Candidate", range(len(candidates)), format_func=lambda i: candidates[i]["candidate_id"], key="peptide-trap-result-" + run_id)
    row = candidates[index]
    path = resolve_stored_path(row["complex_pdb"], run_dir=Path(run_dir), must_exist=True)
    if path:
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
            ChainVisualization, StructureVisualization, molstar_custom_component,
        )
        molstar_custom_component([StructureVisualization(pdb=path.read_text(), representation_type=None, chains=[
            ChainVisualization(chain_id=row["binder_chains"][0], color_params={"value": "0x318c70"}, label="Designed binder"),
            ChainVisualization(chain_id=row["target_chains"][0], color_params={"value": "0xe39c35"},
                               representation_type="cartoon+ball-and-stick", label="Target peptide"),
        ])], key="peptide-trap-result-view-" + run_id + "-" + str(index), height=520, show_controls=True, show_structure_tools=False)
        st.download_button("Download designed complex", path.read_bytes(), file_name=path.name, key="peptide-trap-download-" + run_id)
    with st.expander("Design details"):
        st.json(row.get("raw_metadata") or {})
