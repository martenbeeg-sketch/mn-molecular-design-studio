"""Small task setup panel for the existing Streamlit graph editor."""
import json

import streamlit as st

from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, LEGACY_STARTER_TEMPLATES, new_graph
from mn_molecular_design_studio.scientific_tasks import (
    TASKS, NATIVE, BACKBONE, PREDICTION_MODES, method_limitations,
    task_methods, prediction_methods, new_task_graph, workload,
)


def _prediction_sampling(kind, key):
    """Expose controls that the selected adapter actually submits."""
    fields = {field["name"]: field for field in NODE_CATALOG[kind]["fields"]}
    values = {}
    for name in ("samples", "diffusion_batch_size", "model_count", "seed_count", "seed"):
        if name not in fields:
            continue
        field = fields[name]
        values[name] = st.number_input(field["label"], min_value=int(field.get("min", 0)),
            max_value=int(field.get("max", 2147483647)),
            value=0 if name == "seed" else max(1, int(field.get("min", 1))), key=key + "-" + name)
    return values


def render_task_setup(prefix):
    """Return a new graph only on an explicit create action; never overwrite drafts on selection."""
    with st.expander("Create pipeline by scientific objective", expanded=True):
        task = st.selectbox("Scientific objective", [*TASKS, "blank"],
                            format_func=lambda key: "Blank canvas / custom workflow" if key == "blank" else
                            TASKS[key]["label"] + (" · Experimental" if TASKS[key].get("experimental") else " · Standard"),
                            key=prefix + "-task")
        if task == "blank":
            return new_graph("blank") if st.button("Create editable pipeline", key=prefix + "-create-task") else None
        spec = TASKS[task]
        st.caption(spec["description"])
        st.caption("Methods listed here have compatible registered adapters. Installed models, images, terms and chemical inputs are checked by their existing execution preflight.")
        inputs, config = {}, {}
        mode, binder_format = "protein", "protein"
        st.markdown("**Molecular inputs**")
        if task == "molecular_structure_prediction":
            mode = st.selectbox("Molecular mode", list(PREDICTION_MODES), format_func=PREDICTION_MODES.get, key=prefix + "-mode")
            inputs["target_protein_sequences"] = st.text_area("Target protein sequence(s)", key=prefix + "-target-sequences")
            if mode in {"protein_protein", "protein_macrocyclic"}:
                inputs["partner_protein_sequences"] = st.text_area("Partner protein sequence(s)", key=prefix + "-partner-sequences")
            elif mode == "protein_ligand":
                notation = st.selectbox("Ligand notation", ["smiles", "ccd"], key=prefix + "-notation")
                inputs["ligand_input_type"] = notation
                inputs["ligand_smiles" if notation == "smiles" else "ligand_ccd_codes"] = st.text_area("Ligand " + notation.upper(), key=prefix + "-ligand")
            elif mode in {"protein_dna", "protein_rna"}:
                key = "dna_sequences" if mode == "protein_dna" else "rna_sequences"
                inputs[key] = st.text_area("DNA sequence(s)" if mode == "protein_dna" else "RNA sequence(s)", key=prefix + "-" + key)
            inputs["target_pdb"] = st.text_input("Optional target PDB path", key=prefix + "-prediction-target")
            config["mode"] = mode
        elif task == "antibody_cdr_redesign":
            st.caption("Select existing designs here or use the source node's project candidate picker after graph creation.")
            inputs["source_run_dir"] = st.text_input("Parent design run directory", key=prefix + "-source")
            inputs["selected_candidate_ids"] = [v.strip() for v in st.text_area("Parent candidate IDs (one per line)", key=prefix + "-candidate-ids").splitlines() if v.strip()]
        elif task == "template_peptide_idr":
            inputs["target_sequence"] = st.text_area("Target peptide / IDR fragment sequence", key=prefix + "-fragment")
        elif task != "new_protein_fold":
            inputs["target_pdb"] = st.text_input("Target / conditioning structure path", key=prefix + "-design-target")
            if task == "targeted_protein_binder":
                inputs["target_chains"] = st.text_input("Target chains", value="A", key=prefix + "-target-chains")
                inputs["hotspots"] = st.text_input("Hotspots (blank = exploratory site)", key=prefix + "-hotspots")
                inputs["coldspots"] = st.text_input("Coldspots", key=prefix + "-coldspots")
                binder_format = st.selectbox("Binder format", ["protein", "vhh"], key=prefix + "-format")
                config["binder_format"] = binder_format
        else:
            st.caption("Unconditional generation. Set length/topology and protected motifs in the generator controls.")
        allowed = task_methods(task, mode=mode, binder_format=binder_format, inputs=inputs)
        comparison = len(allowed) > 1 and st.checkbox("Compare selected methods independently", key=prefix + "-comparison")
        label = lambda value: NODE_CATALOG.get(value, {}).get("label", value)
        default = spec["default"] if spec["default"] in allowed else allowed[0]
        if comparison:
            methods = st.multiselect("Methods to compare", allowed, default=[default], format_func=label, key=prefix + "-methods-" + task + "-" + mode + "-" + binder_format)
            st.caption("Each method has its own sequence, prediction, selection and Results branch. Equal attempt counts do not imply equal compute or native acceptance budgets.")
        else:
            methods = [st.selectbox("Method", allowed, index=allowed.index(default), format_func=label,
                                    key=prefix + "-method-" + task + "-" + mode + "-" + binder_format)]
        for method in methods:
            st.caption(label(method) + ": " + method_limitations(method))
        if spec.get("experimental"):
            options = st.text_area("RFD4 ligand / catalytic constraints (JSON object)", value="{}",
                help="Use exact input residue and atom names. Catalytic scaffolding requires motif_residues; active_site_atoms specifies partial-atom anchors. Pocket design accepts ligand atom and interaction controls.",
                key=prefix + "-experimental-options-" + task)
            config["method_settings"] = {"rfd4": {"options_json": options}}
        if task == "molecular_structure_prediction":
            config["method_settings"] = {}
            for method in methods:
                with st.expander("Samples and seed · " + label(method)):
                    config["method_settings"][method] = _prediction_sampling(method, prefix + "-sampling-" + method)
        if task != "small_molecule_hit_prioritization":
            config["comparison"] = bool(comparison)
            if task != "molecular_structure_prediction":
                count_label = ("CDR sequence variants per parent" if task == "antibody_cdr_redesign" else
                               "Designs per register (Peptide Trap) / maximum templates (Kejia)" if task == "template_peptide_idr" else "Generation attempts per method")
                config["num_designs"] = st.number_input(count_label, min_value=1, max_value=10000, value=1, key=prefix + "-attempts")
                if set(methods) & BACKBONE:
                    config["sequences"] = st.number_input("Initial ProteinMPNN sequences per backbone", min_value=1, max_value=1000, value=1, key=prefix + "-sequences")
                has_native = bool(set(methods) & NATIVE)
                validate = not has_native or st.checkbox("Optional additional prediction after native workflow", key=prefix + "-native-validation")
                if validate:
                    prediction_mode = "protein" if task == "new_protein_fold" else "protein_ligand" if spec.get("experimental") else "protein_protein"
                    predictors = prediction_methods(prediction_mode, {"target_pdb": inputs.get("target_pdb", "")})
                    if task == "new_protein_fold":
                        predictors = ["esmfold", *predictors]
                    preferred = "esmfold" if task == "new_protein_fold" else "af2" if task == "template_peptide_idr" and "af2" in predictors else "boltz2"
                    config["prediction"] = st.selectbox("Primary prediction engine", predictors,
                        index=predictors.index(preferred) if preferred in predictors else 0,
                        format_func=label, key=prefix + "-predictor-" + task)
                    config["prediction_settings"] = {config["prediction"]: _prediction_sampling(
                        config["prediction"], prefix + "-primary-sampling-" + config["prediction"])}
            if task in {"targeted_protein_binder", "antibody_cdr_redesign", "template_peptide_idr"}:
                config["soluble_mpnn"] = st.checkbox("Optional SolubleMPNN variants · keep original sequences and interface", key=prefix + "-soluble")
                st.caption("SolubleMPNN protects the mapped interface, motifs, frameworks and cysteines. It does not guarantee improved experimental solubility.")
                config["interface_analysis"] = st.multiselect("Optional interface measurements", ["ipsae", "rosetta_metrics"], key=prefix + "-analysis")
            config["sequence_qc"] = st.checkbox("Optional sequence quality analysis", key=prefix + "-qc")
            if st.checkbox("Optional candidate shortlist", key=prefix + "-selection"):
                metric = st.text_input("Ranking measurement (stage-qualified evidence supported)", key=prefix + "-ranking")
                config["selection"] = {"ranking_metric": metric, "required_metrics": [metric] if metric else [],
                    "limit": st.number_input("Shortlist per comparable protocol", min_value=1, max_value=10000, value=10, key=prefix + "-limit"),
                    "higher_is_better": st.checkbox("Higher measurement ranks first", value=True, key=prefix + "-higher")}
                st.caption("Set scientific thresholds, required finite evidence, manual selections and diversity controls in the filter node. Raw rankings remain separated by protocol.")
        with st.expander("Advanced method and stage settings"):
            st.caption("All generated nodes remain editable with their existing engine forms. JSON here can set method_settings, prediction_settings, soluble_settings, selection, or additional_validation before creation.")
            advanced = st.text_area("Advanced task settings (JSON object)", value="{}", key=prefix + "-advanced-" + task)
        config.update(methods=methods, inputs=inputs)
        if st.button("Create editable pipeline", key=prefix + "-create-task"):
            try:
                values = json.loads(advanced)
                allowed_advanced = {"method_settings", "prediction_settings", "soluble_settings", "selection", "additional_validation"}
                if not isinstance(values, dict) or set(values) - allowed_advanced:
                    raise ValueError("Use the named advanced stage settings only.")
                config.update(values)
                graph = new_task_graph(task, config)
                estimate = workload(graph)
                for line in estimate["counts"]:
                    st.caption(line)
                st.caption(estimate["note"])
                return graph
            except (ValueError, TypeError, KeyError) as exc:
                st.error(str(exc))
    with st.expander("Multi-Generator Binder Design & Evaluation"):
        st.caption("Compose complete native generator workflows, conditional sequence preparation, shared filtering and independent predictor branches.")
        generators = [kind for kind, spec in NODE_CATALOG.items() if spec["category"] in {"Design engine", "Design pipeline"} and not spec.get("hidden")
                      and any("target" in socket["accepts"] for socket in spec["inputs"])]
        predictors = [kind for kind, spec in NODE_CATALOG.items() if spec["category"] == "Structure prediction" and not spec.get("hidden")]
        selected_generators = st.multiselect("Generator branches", generators, default=["rfd4", "bindcraft2"],
            format_func=lambda kind: NODE_CATALOG[kind]["label"], key=prefix + "-multi-generators")
        selected_predictors = st.multiselect("Independent predictor branches", predictors, default=["boltz2", "prediction_protenix_v2"],
            format_func=lambda kind: NODE_CATALOG[kind]["label"], key=prefix + "-multi-predictors")
        st.caption("The first predictor starts enabled. Enable additional branches and set each processing chunk/concurrency in Branch and processing settings.")
        if st.button("Create multi-generator pipeline", key=prefix + "-create-multi"):
            from mn_molecular_design_studio.candidate_pipeline import new_pipeline
            try:
                return new_pipeline(selected_generators, selected_predictors)
            except (ValueError, KeyError) as exc:
                st.error(str(exc))
    with st.expander("Legacy recipes and continuation compatibility"):
        st.caption("Saved graphs and existing candidate handoffs retain their original stages. Engines are also available from the canvas node library and Tools.")
        legacy = st.selectbox("Legacy recipe", list(LEGACY_STARTER_TEMPLATES), key=prefix + "-legacy")
        if st.button("Open legacy recipe", key=prefix + "-open-legacy"):
            return new_graph(LEGACY_STARTER_TEMPLATES[legacy])
    return None
