"""Automatic collection inputs shared by all cofolding engine tiles."""
from pathlib import Path

import numpy as np
import streamlit as st

from mn_molecular_design_studio.design_validation import (
    candidate_passes_filters, design_candidates, protein_chains, read_structure,
)
from mn_molecular_design_studio.cofolding_input import prediction_settings_for, render_molecule_editor


def _acceptance_controls(key, engine_id, mode, source_tools=(), inferred_formats=(), inferred_target_types=()):
    from mn_molecular_design_studio.validation_acceptance import (
        BINDER_FORMATS, TARGET_TYPES, BINDCRAFT2_REFERENCE, available_metrics, bindcraft2_reference_policy,
    )

    known_formats = {value for value in inferred_formats if value in BINDER_FORMATS}
    format_errors = []
    if len(known_formats) > 1:
        st.error("The selected designs include different binder formats. Run acceptance separately for each format so each gets the correct criteria.")
        format_errors.append("Selected designs have mixed binder formats.")
    unknown_formats = any(value is None for value in inferred_formats)
    if unknown_formats:
        st.warning("At least one selected design has no unambiguous format in its recorded provenance. Confirm or override the inferred format below.")
    default_format = next(iter(known_formats)) if len(known_formats) == 1 else "other"
    format_signature = "-".join(sorted(known_formats)) or "unknown"
    format_key = f"{key}-binder-format-{format_signature}-unknown{sum(value is None for value in inferred_formats)}"
    binder_format = st.selectbox("Binder format", list(BINDER_FORMATS),
                                 index=list(BINDER_FORMATS).index(default_format),
                                 format_func=lambda value: BINDER_FORMATS[value],
                                 key=format_key)
    if len(known_formats) == 1:
        st.caption("Format inferred from the selected design campaign and candidate metadata. Change it only if that provenance is incorrect.")
    if unknown_formats and inferred_formats:
        if not st.checkbox(f"Confirm all selected designs are {BINDER_FORMATS[binder_format]}",
                           value=False, key=f"{format_key}-confirm"):
            format_errors.append("Confirm the binder format for candidates without recorded format provenance.")
    target_ids = list(TARGET_TYPES)
    if engine_id == "boltzgen_fold":
        target_ids = ["structured"]
    known_targets = {value for value in inferred_target_types if value in TARGET_TYPES}
    if len(known_targets) > 1:
        st.error("The selected designs use different target input types. Run structured and FASTA / IDR targets separately.")
        format_errors.append("Selected designs have mixed target input types.")
    unknown_targets = any(value is None for value in inferred_target_types)
    if unknown_targets:
        st.warning("At least one selected campaign does not record whether its target was supplied as structure or FASTA / IDR. Confirm the target type below.")
    default_target = next(iter(known_targets)) if len(known_targets) == 1 else "structured"
    target_signature = "-".join(sorted(known_targets)) or "unknown"
    target_key = f"{key}-target-type-{target_signature}-unknown{sum(value is None for value in inferred_target_types)}"
    if engine_id == "boltzgen_fold" and default_target == "fasta_idr":
        st.error("BoltzGen Fold requires a structured target and cannot validate this FASTA / IDR design.")
        format_errors.append("BoltzGen Fold cannot validate a FASTA / IDR target.")
    target_type = st.selectbox("Target type", target_ids, format_func=lambda value: TARGET_TYPES[value],
                               index=target_ids.index(default_target) if default_target in target_ids else 0,
                               key=target_key)
    if len(known_targets) == 1:
        st.caption("Target type inferred from the source campaign and candidate metadata.")
    if unknown_targets and inferred_target_types:
        if not st.checkbox(f"Confirm all selected targets are {TARGET_TYPES[target_type]}",
                           value=False, key=f"{target_key}-confirm"):
            format_errors.append("Confirm the target type for candidates without recorded input provenance.")
    mandatory_rows = available_metrics(binder_format, target_type, mode, engine=engine_id)
    information_rows = available_metrics(binder_format, target_type, mode, engine=engine_id,
                                         source_tools=source_tools, informational=True)
    mandatory_options = {spec["field"]: spec for _, spec in mandatory_rows
                         if spec["direction"] in {"higher", "lower"}}
    information_options = {spec["field"]: spec for _, spec in [*mandatory_rows, *information_rows]}
    mandatory_fields = list(mandatory_options)
    criteria_key = f"{key}-{binder_format}-{target_type}-{mode}"
    reference_rules, inactive_reference = bindcraft2_reference_policy(
        binder_format, target_type, mode, dict(mandatory_rows))
    reference_by_field = {rule["field"]: rule for rule in reference_rules}
    default_filters = [field for field in mandatory_fields if field in reference_by_field]
    st.caption("Starting criteria follow the installed BindCraft2 format and target defaults where this independent engine provides a compatible measurement. Values and operators are editable; each result retains the selected engine and exact measurement method.")
    selected_filters = st.multiselect("Mandatory acceptance measurements", mandatory_fields,
                                      default=default_filters,
                                      format_func=lambda field: f"{mandatory_options[field]['label']} · {mandatory_options[field]['unit']}",
                                      key=f"{criteria_key}-acceptance-filters",
                                      help="Only selected mandatory filters can accept or reject a candidate. An unavailable measurement is unresolved.")
    filters, errors = [], list(format_errors)
    for index, field in enumerate(selected_filters):
        spec = mandatory_options[field]
        left, middle, right = st.columns([1.1, .7, 1.8])
        with left:
            st.caption(spec["label"])
        with middle:
            operators = [">=", ">"] if spec["direction"] == "higher" else ["<=", "<"]
            operator = st.selectbox("Rule", operators,
                                    index=operators.index(reference_by_field[field]["operator"]) if field in reference_by_field else 0,
                                    key=f"{criteria_key}-acceptance-operator-{field}", label_visibility="collapsed")
        with right:
            preset = reference_by_field.get(field, {}).get("value")
            raw_threshold = st.text_input(f"Threshold ({spec['unit']})", value="" if preset is None else f"{preset:g}", key=f"{criteria_key}-acceptance-threshold-{field}",
                                           placeholder=f"Enter a {spec['unit']} threshold")
        if not raw_threshold.strip():
            errors.append(f"Enter a threshold for {spec['label']}.")
            continue
        try:
            threshold = float(raw_threshold)
        except ValueError:
            threshold = float("nan")
        if not np.isfinite(threshold):
            errors.append(f"Threshold for {spec['label']} must be a finite number.")
            continue
        filters.append({"field": field, "operator": operator, "value": threshold,
                        "metric_id": reference_by_field[field]["metric_id"] if field in reference_by_field else ""})
    source_only_bindcraft2 = bool(source_tools) and set(source_tools) == {"bindcraft2"}
    ipdae_field = next((spec["field"] for metric, spec in information_rows if metric == "bindcraft2_ipdae"), None)
    default_information = [ipdae_field] if source_only_bindcraft2 and ipdae_field else []
    score_fields = list(information_options)
    informational = st.multiselect("Informational measurements", score_fields,
                                   default=default_information,
                                   format_func=lambda field: f"{information_options[field]['label']} · {information_options[field]['unit']}",
                                   key=f"{criteria_key}-acceptance-information",
                                   help="These values are recorded and never gate acceptance.")
    rank_options = [""] + score_fields
    ranking_metric = st.selectbox("Optional ranking measurement", rank_options,
                                  format_func=lambda field: "No ranking" if not field else f"{information_options[field]['label']} · {information_options[field]['unit']}",
                                  index=rank_options.index(ipdae_field) if source_only_bindcraft2 and ipdae_field else 0,
                                  key=f"{criteria_key}-acceptance-ranking")
    if target_type == "fasta_idr":
        st.caption("FASTA / IDR mode disables target-coordinate templates and target-aligned pose RMSD. Binder self-consistency can still be measured; target-pose RMSD is unavailable.")
    measured_labels = [spec["label"] for _, spec in mandatory_rows]
    st.caption("Independent engine measurements available for this setup: " + (", ".join(measured_labels) or "none"))
    source_labels = [spec["label"] for _, spec in information_rows if spec.get("source_tools")]
    if source_labels:
        st.caption("Source-design measurements available for display or ranking: " + ", ".join(source_labels))
    if inactive_reference:
        with st.expander("BindCraft2 reference filters unavailable in this independent prediction", expanded=False):
            st.caption("These criteria are not replaced with a differently defined metric. Add a compatible prediction or analysis when they are needed.")
            for item in inactive_reference:
                st.write("• " + item)
    st.caption("Results show actual values and source methods. Missing measurements remain unresolved and cannot pass a selected filter.")
    if errors:
        for error in errors:
            st.error(error)
    policy = {"binder_format": binder_format, "target_type": target_type,
              "mandatory_filters": filters, "informational_scores": informational,
              "ranking_metric": ranking_metric,
              "reference_preset": {"tool": "BindCraft2", "source": BINDCRAFT2_REFERENCE["version"],
                                   "format": binder_format, "target_type": target_type,
                                   "compatible_defaults": reference_rules,
                                   "inactive_criteria": inactive_reference}}
    return policy, target_type, errors


def render_batch_controls(engine_id, profile, project_id, job_name, schema, *, search_jobs, submit,
                          allow_new_batch=True, include_acceptance=False, existing_designs_only=False):
    key = f"design-validation-{profile.get('form_key', engine_id)}-{project_id or 'no-project'}"
    sources = []
    with st.container():
        choices = ["Validate existing designs"] if existing_designs_only else [
            "New molecular inputs", "Sequence collection", "Database molecules", "Validate existing designs"]
        if profile["input_kind"] == "boltzgen" and not existing_designs_only:
            choices = ["New molecular inputs", "Validate existing designs"]
        if "ligand_smiles" in profile.get("molecule_inputs", ()):
            choices.insert(1, "Compound collection")
        kind = st.radio("Input source", choices,
                        key=f"{key}-kind", horizontal=True)
        if kind in {"Compound collection", "Sequence collection", "Database molecules"}:
            _render_collection_inputs(engine_id, profile, project_id, job_name, schema, kind, key, submit)
            return True, 0
        if kind == "New molecular inputs":
            count = int(st.number_input("Independent complexes", min_value=1, max_value=32, value=1,
                                        key=f"{key}-count")) if allow_new_batch else 1
            return False, count
        available_modes = ["Complex refolding"] if profile.get("requires_target_template") or profile["input_kind"] == "boltzgen" else ["Complex refolding", "Binder alone"]
        mode_label = st.radio("Prediction", available_modes, key=f"{key}-mode", horizontal=True)
        mode = "complex" if mode_label == "Complex refolding" else "binder"
        if not project_id:
            st.info("Select a project above to see its successfully completed design campaigns.")
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
                    if (designs := design_candidates(Path(row["run_dir"]), require_passing=not include_acceptance))}
        by_id = {row["run_id"]: row for row in rows if row["run_id"] in eligible}
        source_tools = set()
        inferred_formats = []
        inferred_target_types = []
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
        from mn_molecular_design_studio.validation_acceptance import infer_binder_format, infer_target_type
        workflow_parameters_by_run = {}
        selected = st.multiselect("Design campaigns", list(by_id), key=f"{key}-runs",
                                 format_func=lambda identifier: f"{by_id[identifier].get('campaign_name') or by_id[identifier]['tool']} · {by_id[identifier].get('job_code') or identifier}",
                                 help=("Only completed design jobs with filter-passing outputs are listed."
                                       if not include_acceptance else
                                       "Completed design jobs are listed, including candidates that did not pass their source workflow's native filters."))
        if not by_id:
            st.info("This project has no completed design candidates available for the selected input mode.")
        for identifier in selected:
            run_dir = Path(by_id[identifier]["run_dir"])
            candidates = eligible[identifier]
            by_candidate = {c["candidate_id"]: c for c in candidates}
            chosen = st.multiselect(
                f"Designs · {by_id[identifier].get('job_code') or identifier}", list(by_candidate),
                default=list(by_candidate), key=f"{key}-candidates-{identifier}",
                format_func=lambda candidate: (
                    f"{candidate} · {len(str(by_candidate[candidate].get('binder_sequence') or '').replace(':', ''))} binder residues"
                    + (f" · native {'pass' if candidate_passes_filters(by_candidate[candidate]) else 'no pass recorded'}"
                       if include_acceptance else "")
                ),
            )
            if not candidates:
                st.caption("This campaign has no eligible designs available for validation.")
            if chosen:
                sources.append({"source_run_dir": str(run_dir), "candidate_ids": chosen})
                source_tools.update(str(by_candidate[candidate].get("source_tool") or "") for candidate in chosen)
                if identifier not in workflow_parameters_by_run:
                    try:
                        params = read_json(run_dir / "input.json").get("params", {})
                    except (OSError, ValueError, TypeError):
                        params = {}
                    native_settings = run_dir / "artifacts/raw/bindcraft2/settings.json"
                    if native_settings.is_file():
                        try:
                            params = {"input_params": params, "bindcraft2_settings": read_json(native_settings)}
                        except (OSError, ValueError, TypeError):
                            pass
                    workflow_parameters_by_run[identifier] = params if isinstance(params, dict) else {}
                inferred_formats.extend(
                    infer_binder_format(by_candidate[candidate], workflow_parameters_by_run[identifier])
                    for candidate in chosen
                )
                inferred_target_types.extend(
                    infer_target_type(by_candidate[candidate], workflow_parameters_by_run[identifier])
                    for candidate in chosen
                )
        acceptance_policy, acceptance_errors = None, []
        target_type = "structured"
        if include_acceptance:
            with st.expander("Acceptance criteria", expanded=True):
                acceptance_policy, target_type, acceptance_errors = _acceptance_controls(
                    key, engine_id, mode, source_tools=source_tools,
                    inferred_formats=inferred_formats,
                    inferred_target_types=inferred_target_types)
            if acceptance_policy.get("binder_format") == "cyclic_peptide" and not profile.get("supports_cyclic"):
                message = f"{profile['label']} does not preserve cyclic peptide topology. Choose a co-folding engine with cyclic-chain support."
                st.error(message)
                acceptance_errors.append(message)
        target_msa = st.checkbox("Use target MSA", value=True, key=f"{key}-target-msa",
                                 disabled=mode == "binder" or not profile.get("supports_msa"),
                                 help="Reuse the shared reference MSA or fetch it from ColabFold. The designed binder always uses a single sequence.")
        target_template = st.checkbox("Use target structure as a template", value=True, key=f"{key}-target-template",
                                      disabled=bool(mode == "binder" or target_type == "fasta_idr" or not profile.get("supports_templates") or profile.get("requires_target_template")),
                                      help="Supply only the target coordinates; the designed binder is excluded from the input template.")
        target_template = bool(target_template and target_type == "structured")
        if mode == "complex" and not profile.get("supports_templates"):
            st.caption("The installed engine accepts sequence inputs without a supplied structure template.")
        elif profile.get("requires_target_template"):
            st.caption("This engine requires the target structure.")
        st.caption(f"{sum(len(s['candidate_ids']) for s in sources)} selected designs. Validation measurements and acceptance decisions appear in Results.")
    config = prediction_settings_for(profile)
    config["fields"] = [field for field in config["fields"] if field["key"] != "use_msa"]
    _, settings, submitted = render_molecule_editor(
        f"{key}-settings", (), settings_only=True, prediction_settings=config,
        can_submit=bool(project_id and sources and not acceptance_errors), modes=tuple(profile["modes"]),
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
            if include_acceptance:
                parameters["validation_target_type"] = target_type
                parameters["acceptance_policy"] = acceptance_policy
            submit(schema, "refolding.optimized_prediction", parameters, project_id, settings, profile["label"])
        except (ValueError, TypeError, KeyError) as exc:
            st.error(f"Could not submit design validation: {exc}")
    return True, 0


def _render_collection_inputs(engine_id, profile, project_id, job_name, schema, kind, key, submit):
    import hashlib
    from mn_molecular_design_studio import batching
    from mn_molecular_design_studio.auth import user_input_root
    from mn_molecular_design_studio.ligand_compounds import parse_compound_table, write_compound_set
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import utc_now
    user = st.session_state.get("_molecular_database_user")
    if not user or not project_id:
        st.info("Select a project to choose the inputs for this calculation.")
        return
    input_parameters, count, tasks = {}, 0, None
    prefix = f"{key}-collection"
    if kind == "Compound collection":
        from mn_molecular_design_studio.molecular_database import MolecularDatabase
        collections = MolecularDatabase().compound_options(user)
        pending = st.session_state.get(prefix + "-selected")
        if pending:
            collections = [pending, *[row for row in collections if row["path"] != pending["path"]]]
        source = st.radio("Compounds", ["My compounds", "Upload file"], key=prefix + "-source", horizontal=True)
        if source == "My compounds":
            if collections:
                selected = st.selectbox("Compound collection", collections, format_func=lambda row: f"{row['name']} · {row['compound_count']:,} compounds", key=prefix + "-compounds")
                input_parameters["compound_set_file"] = selected["path"]
                count = selected["compound_count"]
            else:
                st.info("Import compounds into Database or upload a compound file here.")
        else:
            uploaded = st.file_uploader("Compound file", type=["sdf", "csv", "tsv", "smi", "smiles"], key=prefix + "-upload")
            if uploaded:
                try:
                    content = uploaded.getvalue()
                    compounds = parse_compound_table(content, filename=uploaded.name)
                    path = user_input_root(user["id"]) / "prediction_collections" / (hashlib.sha256(content).hexdigest() + ".json")
                    write_compound_set(path, name=uploaded.name, compounds=compounds,
                                       owner_id=user["id"], project_id=project_id, created_at=utc_now())
                    input_parameters["compound_set_file"] = str(path)
                    count = len(compounds)
                except ValueError as exc:
                    st.error(str(exc))
        st.caption("Enter the shared target below. Every selected compound is predicted independently in one job.")
    elif kind == "Sequence collection":
        uploaded = st.file_uploader("Sequence file", type=["fasta", "fa", "faa", "txt"], key=prefix + "-sequences")
        text = uploaded.getvalue().decode("utf-8-sig") if uploaded else st.text_area("Sequences", key=prefix + "-sequence-text", placeholder=">sequence_1\nACDEFGHIK\n>sequence_2\nACDEFGHIL")
        if text.strip():
            try:
                count = len(batching.sequence_tasks(text))
                path = user_input_root(user["id"]) / "prediction_collections" / (hashlib.sha256(text.encode()).hexdigest() + ".fasta")
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
                input_parameters["sequence_file"] = str(path)
            except ValueError as exc:
                st.error(str(exc))
        st.caption("Each FASTA record or plain sequence line is an independent prediction. Use ':' between chains belonging to the same complex.")
    else:
        from mn_molecular_design_studio.database_ui import render_molecule_picker
        from mn_molecular_design_studio.database_selections import editor_model
        from mn_molecular_design_studio.cofolding_input import model_to_payload
        allowed = {"protein", *[kind for kind, field in (("dna", "dna_sequences"), ("rna", "rna_sequences"), ("small_molecule", "ligand_smiles")) if field in profile.get("molecule_inputs", ())]}
        selected = render_molecule_picker(user, key=prefix, allowed_types=allowed, multiple=True)
        if selected:
            st.session_state[prefix + "-molecules"] = selected
        records = st.session_state.get(prefix + "-molecules") or []
        tasks = []
        for record in records:
            values = model_to_payload(editor_model(record))
            values.pop("input_error", None)
            values.pop("chain_catalog", None)
            tasks.append({**values, "name": record["name"]})
        count = len(tasks)
    if count:
        st.caption(f"{count:,} selected inputs → one job. Per-input structures, scores and errors are retained.")
        if count > batching.MAX_PREDICTION_ITEMS:
            st.error(f"Select at most {batching.MAX_PREDICTION_ITEMS:,} prediction inputs for one run.")
    config = prediction_settings_for(profile)
    if engine_id == "boltz2" and kind == "Compound collection":
        config["fields"] = [{**field, "default": True} if field["key"] == "use_potentials" else field for field in config["fields"]]
    values, settings, submitted = render_molecule_editor(
        prefix + "-settings", tuple(field for field in profile.get("molecule_inputs", ()) if field not in {"ligand_smiles", "ligand_ccd_codes"}),
        prediction_settings=config, settings_only=kind != "Compound collection",
        supports_templates=bool(profile.get("supports_templates")),
        modes=tuple(profile["modes"]), can_submit=bool(count and count <= batching.MAX_PREDICTION_ITEMS),
    )
    if not submitted:
        return
    if kind == "Compound collection" and (values.get("input_error") or not values.get("protein_sequences", "").strip()):
        st.error("Enter the target protein sequence or select a target from My molecules.")
        return
    mode = str(settings.get("optimization_mode", "off"))
    parameters = {
        "engine": engine_id, "job_name": job_name, "optimization_mode": mode,
        "seed": int(settings.get("seed") or 0), "samples": int(settings.get("samples") or 1),
        "num_recycles": int(settings.get("num_recycles", profile.get("default_recycles", 10))),
        "cpu_cores": int(settings.get("cpu_cores") or 4),
        "gpu_device": "0,1" if profile.get("gpu_count", {}).get(mode, 1) == 2 else "0",
        "accept_alphafold3_terms": bool(settings.get("accept_terms")),
        "engine_options": {field["key"]: settings[field["key"]] for field in profile.get("settings", ()) if field["key"] in settings},
        **input_parameters,
    }
    if kind == "Compound collection":
        parameters["protein_sequences"] = values["protein_sequences"]
        # Preserve target chain templates/MSA choices in every expanded item.
        context = {field: values[field] for field in ("protein_sequences", "dna_sequences", "rna_sequences", "molecule_chain_options", "structure_templates", "boltz_templates") if field in values}
        tasks = batching.prediction_collection_tasks(context, compound_set_file=input_parameters["compound_set_file"])
        parameters.pop("compound_set_file", None)
        parameters.pop("protein_sequences", None)
    if tasks is not None:
        parameters["tasks"] = tasks
    try:
        submit(schema, "refolding.optimized_prediction", parameters, project_id, settings, profile["label"])
    except (ValueError, TypeError, KeyError) as exc:
        st.error(f"Could not submit the selected inputs: {exc}")


def comparison_viewer(original, aligned, reference, chain_map):
    """Preserve the fitted frame and accept native CIFs the regular viewer skips."""
    from py2Dmol import view
    viewer = view(size=(760, 510), controls=True, box=True, style="richardson", bg="white", multi=True)
    binder_only = reference["mode"] == "binder" or reference.get("target_type") == "fasta_idr"
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
    decision = raw.get("validation_acceptance") or {}
    if decision:
        st.subheader("Validation & acceptance")
        st.metric("Configured acceptance", decision.get("status", "unresolved").replace("_", " ").title())
        if decision.get("native_filter_status") is not None:
            st.caption(f"Native workflow filter status: {decision['native_filter_status']} · reported separately from this configured acceptance policy.")
        measurement_rows = []
        for field, item in decision.get("measurements", {}).items():
            threshold = f"{item.get('operator')} {item.get('threshold')}" if item.get("mandatory") else "Informational"
            value = item.get("value")
            measurement_rows.append({"Measurement": field, "Value": "Unavailable" if value is None else value,
                                    "Unit": item.get("unit") or "source-native", "Status": item.get("status"),
                                    "Criterion": threshold, "Passed": item.get("passed"),
                                    "Source": item.get("source_tool") or "unknown", "Method": item.get("method") or "source method retained"})
        if measurement_rows:
            st.dataframe(measurement_rows, hide_index=True, use_container_width=True)
        reasons = decision.get("reasons") or []
        if reasons:
            st.caption("Acceptance reasons")
            for reason in reasons:
                st.write("• " + str(reason))
        elif decision.get("status") == "not_evaluated":
            st.caption("No mandatory acceptance filters were configured; informational measurements do not gate a result.")
    reference = raw.get("validation_reference")
    if not reference:
        return
    st.subheader("Design comparison")
    st.caption(f"Reference design: {reference['candidate_id']}")
    if metrics.get("validation_status") != "completed":
        st.error(metrics.get("validation_error", "Structure comparison is unavailable."))
        return
    complex_mode = reference["mode"] == "complex"
    target_pose_comparable = complex_mode and reference.get("target_type", "structured") != "fasta_idr"
    columns = st.columns(2)
    value = metrics.get("binder_rmsd_target_aligned" if target_pose_comparable else "binder_rmsd_fitted")
    columns[0].metric("Binder Cα RMSD after target alignment" if target_pose_comparable else "Binder Cα self-consistency RMSD", f"{value:.3f} Å")
    columns[1].metric("Compared binder residues", str(metrics["binder_rmsd_ca_count"]))
    if target_pose_comparable:
        st.caption(f"Target fit: {metrics['target_alignment_rmsd']:.3f} Å over {metrics['target_alignment_ca_count']} Cα atoms. The binder was measured in that target-aligned frame.")
    elif complex_mode:
        st.caption("The target is treated as unconstrained. Binder self-consistency is measured after binder-only alignment; relative target-pose RMSD is unavailable.")
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
