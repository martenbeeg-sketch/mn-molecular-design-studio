"""Project and design selection inside each cofolding engine's batch controls."""
from pathlib import Path

import numpy as np
import streamlit as st

from mn_molecular_design_studio.design_validation import design_candidates, protein_chains, read_structure
from mn_molecular_design_studio.cofolding_input import prediction_settings_for, render_molecule_editor


def render_batch_controls(engine_id, profile, project_id, job_name, schema, *, search_jobs, submit, allow_new_batch=True):
    key = f"design-validation-{profile.get('form_key', engine_id)}-{project_id or 'no-project'}"
    sources = []
    with st.expander("Batch mode", expanded=False):
        kind = st.radio("Batch input", ["New molecular inputs", "Validate existing designs"],
                        key=f"{key}-kind", horizontal=True)
        if kind == "New molecular inputs":
            count = int(st.number_input("Independent complexes", min_value=1, max_value=32, value=1,
                                        key=f"{key}-count")) if allow_new_batch else 1
            return False, count
        available_modes = ["Complex refolding"] if profile.get("requires_target_template") or profile["input_kind"] == "boltzgen" else ["Complex refolding", "Binder alone"]
        mode_label = st.radio("Prediction", available_modes, key=f"{key}-mode", horizontal=True)
        mode = "complex" if mode_label == "Complex refolding" else "binder"
        if not project_id:
            st.info("Select a project above to see its successful design campaigns with filter-passing designs.")
            return True, 0
        rows, cursor = [], None
        while True:
            page = search_jobs(project_id=project_id, task_group="design", statuses=["completed", "succeeded"],
                               limit=1000, cursor=cursor)
            rows.extend(page["items"])
            cursor = page.get("next_cursor")
            if not cursor:
                break
        eligible = {row["run_id"]: designs for row in rows
                    if (designs := design_candidates(Path(row["run_dir"])))}
        by_id = {row["run_id"]: row for row in rows if row["run_id"] in eligible}
        selected = st.multiselect("Design campaigns", list(by_id), key=f"{key}-runs",
                                 format_func=lambda identifier: f"{by_id[identifier].get('campaign_name') or by_id[identifier]['tool']} · {by_id[identifier].get('job_code') or identifier}",
                                 help="Only successfully completed design jobs in this project with at least one filter-passing design are listed.")
        if not by_id:
            st.info("This project has no successful design campaigns with filter-passing designs.")
        for identifier in selected:
            run_dir = Path(by_id[identifier]["run_dir"])
            candidates = eligible[identifier]
            by_candidate = {c["candidate_id"]: c for c in candidates}
            chosen = st.multiselect(
                f"Designs · {by_id[identifier].get('job_code') or identifier}", list(by_candidate),
                default=list(by_candidate), key=f"{key}-candidates-{identifier}",
                format_func=lambda candidate: f"{candidate} · {len(str(by_candidate[candidate].get('binder_sequence') or '').replace(':', ''))} binder residues",
            )
            if not candidates:
                st.caption("This campaign has no filter-passing designs available for validation.")
            if chosen:
                sources.append({"source_run_dir": str(run_dir), "candidate_ids": chosen})
        target_msa = st.checkbox("Use target MSA", value=True, key=f"{key}-target-msa",
                                 disabled=mode == "binder" or not profile.get("supports_msa"),
                                 help="Reuse the shared reference MSA or fetch it from ColabFold. The designed binder always uses a single sequence.")
        target_template = st.checkbox("Use target structure as a template", value=True, key=f"{key}-target-template",
                                      disabled=bool(mode == "binder" or not profile.get("supports_templates") or profile.get("requires_target_template")),
                                      help="Supply only the target coordinates; the designed binder is excluded from the input template.")
        if mode == "complex" and not profile.get("supports_templates"):
            st.caption("The installed engine accepts sequence inputs without a supplied structure template.")
        elif profile.get("requires_target_template"):
            st.caption("This engine requires the target structure.")
        st.caption(f"{sum(len(s['candidate_ids']) for s in sources)} selected designs. RMSD and overlays appear in the prediction results.")
    config = prediction_settings_for(profile)
    config["fields"] = [field for field in config["fields"] if field["key"] != "use_msa"]
    _, settings, submitted = render_molecule_editor(
        f"{key}-settings", (), settings_only=True, prediction_settings=config,
        can_submit=bool(project_id and sources), modes=tuple(profile["modes"]),
    )
    if submitted:
        try:
            inference_mode = str(settings.get("optimization_mode", "off"))
            parameters = {
                "engine": engine_id, "job_name": job_name or f"{engine_id}-validation",
                "optimization_mode": inference_mode, "seed": int(settings.get("seed") or 0),
                "samples": int(settings["samples"]), "cpu_cores": int(settings["cpu_cores"]),
                "num_recycles": int(settings.get("num_loops", settings.get("num_recycles", profile.get("default_recycles", 10)))),
                "gpu_device": "0,1" if profile.get("gpu_count", {}).get(inference_mode, 1) == 2 else "0",
                "accept_alphafold3_terms": bool(settings.get("accept_terms")),
                "engine_options": {field["key"]: settings[field["key"]] for field in profile.get("settings", ()) if field["key"] in settings},
                "validation_sources": sources, "validation_mode": mode,
                "target_use_msa": bool(target_msa and mode == "complex" and profile.get("supports_msa")),
                "use_target_template": bool(target_template and mode == "complex" and profile.get("supports_templates")),
            }
            submit(schema, "refolding.optimized_prediction", parameters, project_id, settings, profile["label"])
        except (ValueError, TypeError, KeyError) as exc:
            st.error(f"Could not submit design validation: {exc}")
    return True, 0


def comparison_viewer(original, aligned, reference, chain_map):
    """Preserve the fitted frame and accept native CIFs the regular viewer skips."""
    from py2Dmol import view
    viewer = view(size=(760, 510), controls=True, box=True, style="richardson", bg="white", multi=True)
    binder_only = reference["mode"] == "binder"
    for path, name, color, selected in (
        (original, "Designed reference", "#aaaaaa", reference["binder_chains"] if binder_only else None),
        (aligned, "Refolded prediction", "#1599a9", list(chain_map.values()) if binder_only else None),
    ):
        residues = [(chain, residue) for chain, rows in protein_chains(read_structure(path)).items()
                    if selected is None or chain in selected for residue in rows]
        if not residues:
            raise ValueError(f"{name} has no protein coordinates for this overlay.")
        viewer.add(np.array([r["CA"].coord for _, r in residues]),
                   chains=[c for c, _ in residues], position_types=["P"] * len(residues),
                   position_names=[r.resname for _, r in residues], residue_numbers=[r.id[1] for _, r in residues],
                   name=name, color=color, align=False)
    return viewer


def render_comparison(run_dir, candidate, *, key):
    from mn_molecular_design_studio.cofolding_results import local_artifact
    metrics = candidate.get("metrics") or {}
    raw = candidate.get("raw_metadata") or {}
    reference = raw.get("validation_reference")
    if not reference:
        return
    st.subheader("Design comparison")
    st.caption(f"Reference design: {reference['candidate_id']}")
    if metrics.get("validation_status") != "completed":
        st.error(metrics.get("validation_error", "Structure comparison is unavailable."))
        return
    complex_mode = reference["mode"] == "complex"
    columns = st.columns(2)
    value = metrics.get("binder_rmsd_target_aligned" if complex_mode else "binder_rmsd_fitted")
    columns[0].metric("Binder Cα RMSD after target alignment" if complex_mode else "Binder Cα RMSD after binder alignment", f"{value:.3f} Å")
    columns[1].metric("Compared binder residues", str(metrics["binder_rmsd_ca_count"]))
    if complex_mode:
        st.caption(f"Target fit: {metrics['target_alignment_rmsd']:.3f} Å over {metrics['target_alignment_ca_count']} Cα atoms. The binder was measured in that target-aligned frame.")
    else:
        st.caption("The predicted binder was fitted to the designed binder using matched Cα atoms.")
    st.caption(f"Binder coordinate coverage: {metrics['binder_rmsd_coverage']:.1%}")
    aligned = local_artifact(run_dir, raw.get("aligned_structure_path"))
    original = local_artifact(run_dir, raw.get("reference_structure_path"))
    if aligned and original:
        with st.expander("Structure overlay", expanded=True):
            import streamlit.components.v1 as components
            viewer = comparison_viewer(original, aligned, reference, raw.get("validation_chain_map", {}))
            components.html(viewer.to_html(title="Designed and refolded structure overlay"), height=585, scrolling=False)
            st.caption("Grey: designed structure. Teal: aligned prediction.")
            st.download_button("Download aligned prediction", aligned.read_bytes(), file_name=aligned.name,
                               mime="chemical/x-mmcif", key=f"{key}-aligned")
