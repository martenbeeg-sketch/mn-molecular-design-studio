"""Task-first recipes composed from existing canvas nodes, never saved-graph migrations.

Method IDs are node kinds. A recipe is a starting point, not an engine contract or
a claim of experimental validity. Native workflows remain indivisible nodes.
"""
from __future__ import annotations

import copy
import json
import uuid


TASKS = {
    "targeted_protein_binder": {
        "label": "Targeted Protein Binder Design", "default": "rfd3",
        "methods": ["rfd3", "rfdiffusion", "rfd4", "bindcraft2", "boltzgen", "pxdesign",
                    "genie3", "complexa_generator", "esm_binder", "rfdiffusion_beta_pairing",
                    "bindcraft2_generator", "boltzgen_generator", "complexa", "pxdesign_pipeline", "genie3_pipeline"],
        "description": "Design against a folded protein target. Choose one method or explicitly compare independent methods.",
    },
    "antibody_cdr_redesign": {
        "label": "Antibody CDR Redesign", "default": "abmpnn", "methods": ["abmpnn"],
        "description": "Redesign selected CDRs of existing antibodies; protect frameworks, cysteines and fixed motifs.",
    },
    "template_peptide_idr": {
        "label": "Template-Based Peptide / IDR Recognition", "default": "peptide_trap",
        "methods": ["peptide_trap", "kejia_idr"],
        "description": "Peptide Trap uses compatible beta-strand clefts. Kejia uses flexible template threading; select its template and fragment explicitly.",
    },
    "molecular_structure_prediction": {
        "label": "Molecular Structure Prediction", "default": "boltz2", "methods": [],
        "description": "Choose molecular mode and one compatible predictor. Compare structural conclusions; raw engine ranking scores are not calibrated across methods.",
    },
    "small_molecule_hit_prioritization": {
        "label": "Small-Molecule Hit Prioritization", "default": "unidock_active_learning",
        "methods": ["unidock_active_learning"],
        "description": "Use the existing adaptive docking campaign for large libraries; direct docking remains available in Tools. Surrogate labels represent docking scores.",
    },
    "new_protein_fold": {
        "label": "New Protein Fold Design", "default": "protpardelle", "methods": ["protpardelle"],
        "description": "Unconditional fold generation, sequence design and monomer self-consistency. This task does not establish target binding.",
    },
    "ligand_conditioned_pocket": {
        "label": "Ligand-Conditioned Protein-Pocket Design", "default": "rfd4", "methods": ["rfd4"],
        "experimental": True,
        "description": "Specify exact ligand chemistry and atom constraints. Pocket geometry and prediction do not establish affinity.",
    },
    "catalytic_motif_scaffolding": {
        "label": "Catalytic Motif Scaffolding", "default": "rfd4", "methods": ["rfd4"],
        "experimental": True,
        "description": "Supply reaction-specific motif residues, active-site atoms and cofactors. Scaffold confidence does not establish catalytic activity.",
    },
}

PREDICTION_MODES = {
    "protein": "Protein monomer", "protein_protein": "Protein–protein complex",
    "protein_ligand": "Protein–ligand complex", "protein_dna": "Protein–DNA complex",
    "protein_rna": "Protein–RNA complex", "protein_macrocyclic": "Cyclic protein/peptide partner",
}
NATIVE = {"bindcraft2", "boltzgen", "complexa", "pxdesign_pipeline", "genie3_pipeline"}
BACKBONE = {"rfd3", "rfdiffusion", "pxdesign", "genie3", "protpardelle", "kejia_idr", "rfdiffusion_beta_pairing"}
VHH_METHODS = {"bindcraft2_generator", "boltzgen_generator", "esm_binder"}
TEMPLATES = {"rfd3": "foundry", "rfd4": "rfd4_protein_binder", "esm_binder": "esmfold2_binder",
             "complexa_generator": "complexa"}
CONFIG_KEYS = {"methods", "comparison", "mode", "binder_format", "inputs", "method_settings", "prediction",
               "prediction_settings", "num_designs", "sequences", "initial_model", "soluble_mpnn",
               "soluble_settings", "sequence_qc", "interface_analysis", "selection", "additional_validation"}


def prediction_methods(mode="protein", inputs=None):
    """Use the same capability profiles/preflight as the existing prediction tools."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, _cofold_model_kinds
    if mode not in PREDICTION_MODES:
        raise ValueError("Choose a supported molecular prediction mode.")
    values = inputs or {}
    return [kind for kind in _cofold_model_kinds(mode)
            if not (NODE_CATALOG[kind].get("requires_target_structure") and not values.get("target_pdb"))
            and not (mode == "protein_ligand" and values.get("ligand_input_type", "smiles")
                     not in NODE_CATALOG[kind].get("ligand_input_formats", []))]


def task_methods(task_id, *, mode="protein", binder_format="protein", inputs=None):
    if task_id not in TASKS:
        raise ValueError("Unknown scientific task.")
    if task_id == "molecular_structure_prediction":
        return prediction_methods(mode, inputs)
    if binder_format != "protein" and task_id != "targeted_protein_binder":
        raise ValueError("Binder format selection applies to targeted binder generation.")
    if binder_format not in {"protein", "vhh"}:
        raise ValueError("Choose protein or VHH binder format.")
    methods = TASKS[task_id]["methods"]
    if task_id == "targeted_protein_binder" and binder_format == "vhh":
        return [method for method in methods if method in VHH_METHODS]
    return list(methods)


def catalogue():
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    result = []
    for key, spec in TASKS.items():
        methods = task_methods(key)
        result.append({"id": key, **copy.deepcopy(spec), "experimental": bool(spec.get("experimental")),
                       "methods": [{"id": method, "label": NODE_CATALOG.get(method, {}).get("label", method),
                                    "native": method in NATIVE, "backbone_only": method in BACKBONE,
                                    "limitations": method_limitations(method)} for method in methods],
                       "prediction_modes": PREDICTION_MODES if key == "molecular_structure_prediction" else {},
                       "methods_by_mode": {mode: prediction_methods(mode) for mode in PREDICTION_MODES}
                           if key == "molecular_structure_prediction" else {},
                       "methods_by_binder_format": {fmt: task_methods(key, binder_format=fmt) for fmt in ("protein", "vhh")}
                           if key == "targeted_protein_binder" else {},
                       "availability": "Compatible registered adapters; installed images, checkpoints, terms and input chemistry remain subject to engine preflight."})
    return result


def method_limitations(method):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    specific = {
        "rfdiffusion_beta_pairing": "Requires explicit compatible sheet/strand conditioning and motif mapping.",
        "kejia_idr": "Flexible template strategy; confidence in beta pairing does not establish native IDR accessibility. Exploratory method.",
        "peptide_trap": "Requires a compatible published cleft/register; computational predictions do not establish anti-aggregation activity.",
        "rfd4": "Co-designed sequences are retained; component/atom mapping restrictions are enforced by the existing handoff.",
    }
    return specific.get(method, NODE_CATALOG.get(method, {}).get("description", TASKS["small_molecule_hit_prioritization"]["description"]))


def _node(kind, identifier, params=None, label=None):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    spec = NODE_CATALOG[kind]
    values = {f["name"]: copy.deepcopy(f["default"]) for f in spec["fields"]}
    if params:
        values.update(copy.deepcopy(params))
    return {"id": identifier, "type": "studio", "position": {"x": 0, "y": 0},
            "data": {"kind": kind, "label": label or spec["label"], "params": values}}


def _update(node, values):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    allowed = {f["name"] for f in NODE_CATALOG[node["data"]["kind"]]["fields"]}
    if node["data"]["kind"] in {"generation_source", "backbone_source", "project_design_source"}:
        allowed.add("selected_candidate_ids")
    if not isinstance(values, dict) or set(values) - allowed:
        raise ValueError(f"Unsupported settings for {node['data']['label']}.")
    node["data"]["params"].update(copy.deepcopy(values))


def _predictor(kind, identifier, settings):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    if kind not in NODE_CATALOG or NODE_CATALOG[kind]["category"] != "Structure prediction":
        raise ValueError("Choose a compatible structure-prediction engine.")
    node = _node(kind, identifier)
    definitions = {f["name"]: f for f in NODE_CATALOG[kind]["fields"]}
    # Compact defaults apply only to newly composed recipes. Explicit settings win.
    for key, value in {"samples": 1, "diffusion_samples": 1, "diffusion_batch_size": 1,
                       "num_recycles": 3, "n_recycles": 3, "model_count": 1}.items():
        if key in node["data"]["params"]:
            node["data"]["params"][key] = max(value, definitions[key].get("min", value))
    _update(node, settings)
    return node


def new_task_graph(task_id, configuration=None):
    """Generate independent editable branches; legacy new_graph constructors stay intact."""
    from mn_molecular_design_studio.pipeline_graph import SCHEMA, new_graph, validate_graph
    config = copy.deepcopy(configuration or {})
    if not isinstance(config, dict) or set(config) - CONFIG_KEYS:
        raise ValueError("Unsupported scientific task configuration.")
    inputs = config.get("inputs", {})
    if not isinstance(inputs, dict):
        raise ValueError("Molecular inputs must be an object.")
    mode, binder_format = config.get("mode", "protein"), config.get("binder_format", "protein")
    if task_id == "molecular_structure_prediction":
        molecular_fields = {"protein": set(), "protein_protein": {"partner_protein_sequences"},
                            "protein_macrocyclic": {"partner_protein_sequences"},
                            "protein_ligand": {"ligand_ccd_codes", "ligand_smiles", "compound_set_file"},
                            "protein_dna": {"dna_sequences"}, "protein_rna": {"rna_sequences"}}
        molecule_keys = set().union(*molecular_fields.values())
        if any(inputs.get(key) for key in molecule_keys - molecular_fields.get(mode, set())):
            raise ValueError("Molecular inputs must match the selected prediction mode.")
    allowed = task_methods(task_id, mode=mode, binder_format=binder_format, inputs=inputs)
    default = TASKS[task_id]["default"]
    if default not in allowed:
        default = allowed[0] if allowed else None
    methods = config.get("methods", [default])
    if not isinstance(methods, list) or not methods or any(m not in allowed for m in methods) or len(set(methods)) != len(methods):
        raise ValueError("Choose distinct methods compatible with this scientific task and its molecular inputs.")
    if len(methods) > 1 and config.get("comparison") is not True:
        raise ValueError("Enable explicit comparison mode to run multiple methods.")
    if config.get("comparison") and len(methods) < 2:
        raise ValueError("A comparison requires at least two methods.")
    for flag in ("comparison", "soluble_mpnn", "sequence_qc"):
        if flag in config and type(config[flag]) is not bool:
            raise ValueError(f"{flag} must be true or false.")
    for key in ("num_designs", "sequences"):
        if key in config and (type(config[key]) is not int or not 1 <= config[key] <= 10000):
            raise ValueError(f"{key} must be a positive integer, at most 10000.")
    method_settings = config.get("method_settings", {})
    if not isinstance(method_settings, dict) or set(method_settings) - set(methods):
        raise ValueError("Method settings must name selected methods only.")
    prediction_settings = config.get("prediction_settings", {})
    if not isinstance(prediction_settings, dict):
        raise ValueError("Prediction settings must be an object keyed by prediction engine.")
    soluble_allowed = task_id in {"targeted_protein_binder", "antibody_cdr_redesign", "template_peptide_idr"}
    if config.get("soluble_mpnn") and not soluble_allowed:
        raise ValueError("Interface-preserving SolubleMPNN requires a protein target–binder design task.")
    initial_model = config.get("initial_model", "protein_mpnn")
    if initial_model != "protein_mpnn":
        raise ValueError("Initial backbone preparation currently supports ProteinMPNN. Other deliberate redesign models remain available in the editor.")
    analysis = config.get("interface_analysis", [])
    if not isinstance(analysis, list) or any(a not in {"ipsae", "rosetta_metrics"} for a in analysis):
        raise ValueError("Choose supported optional interface analyses.")
    if analysis and task_id not in {"targeted_protein_binder", "antibody_cdr_redesign", "template_peptide_idr"}:
        raise ValueError("Protein interface analysis requires a protein target–binder task.")
    if task_id == "small_molecule_hit_prioritization":
        if set(config) - {"methods", "inputs", "method_settings"}:
            raise ValueError("Configure optional screening stages in the existing adaptive campaign editor.")
        graph = new_graph("unidock_active_learning")
        _update(next(n for n in graph["nodes"] if n["data"]["kind"] == "active_target"), inputs)
        # Keep its template ID and execution dispatch, including feedback scheduling.
        for kind, values in method_settings.get(methods[0], {}).items():
            node = next((n for n in graph["nodes"] if n["data"]["kind"] == kind), None)
            if node is None:
                raise ValueError("Unknown adaptive campaign node.")
            _update(node, values)
        graph["task_configuration"] = {"task_id": task_id, "version": 1, **config}
        return validate_graph(graph)["graph"]

    nodes, edges = [], []
    def chain(branch, additions):
        for item in additions:
            if branch:
                edges.append({"id": f"edge-{branch[-1]['id']}-{item['id']}", "source": branch[-1]["id"],
                              "target": item["id"], "sourceHandle": "out", "targetHandle": "in"})
            branch.append(item)

    for index, method in enumerate(methods):
        prefix = f"method-{index}-{method}"
        branch = []
        if task_id == "molecular_structure_prediction":
            source = _node("cofold_input", prefix + "-input", {"task_type": mode, "cyclic_partner": mode == "protein_macrocyclic"})
            if "task_type" in inputs or "cyclic_partner" in inputs:
                raise ValueError("Set prediction mode through mode, rather than overriding molecular input types.")
            _update(source, inputs)
            predictor = _predictor(method, prefix + "-prediction", {
                **prediction_settings.get(method, {}), **method_settings.get(method, {})})
            if method == "af2":
                predictor["data"]["params"]["af2_validation_mode"] = "target_only"
            chain(branch, [source, predictor])
            if config.get("sequence_qc"):
                chain(branch, [_node("sequence_qc", prefix + "-qc")])
        else:
            template = ("abmpnn_cdr_redesign" if task_id == "antibody_cdr_redesign" else
                        "rfd4_small_molecule" if task_id == "ligand_conditioned_pocket" else
                        "rfd4_enzyme" if task_id == "catalytic_motif_scaffolding" else
                        "vhh_design" if binder_format == "vhh" else TEMPLATES.get(method, method))
            if method in NATIVE - {"bindcraft2", "boltzgen"} or method in {"bindcraft2_generator", "boltzgen_generator"} and binder_format != "vhh":
                generator = _node(method, prefix + "-generation")
                originals = [_node("target", "target"), _node("binder_context", "context"), generator]
            else:
                originals = new_graph(template)["nodes"]
            generation_kinds = {method}
            retained = [copy.deepcopy(n) for n in originals if n["data"]["kind"] in
                        {"target", "binder_context", "project_design_source", *generation_kinds,
                         *({"sabr"} if task_id == "antibody_cdr_redesign" else set())}]
            if task_id == "targeted_protein_binder" and not any(n["data"]["kind"] == "target" for n in retained):
                retained = [_node("target", "target"), _node("binder_context", "context"), *retained]
            for number, item in enumerate(retained):
                item["id"] = f"{prefix}-stage-{number}"
            generator = next(n for n in retained if n["data"]["kind"] == method)
            # Experimental recipes must not copy tutorial-specific atom/residue selectors.
            if task_id in {"ligand_conditioned_pocket", "catalytic_motif_scaffolding"}:
                generator["data"]["params"]["options_json"] = "{}"
            if "num_designs" in generator["data"]["params"]:
                generator["data"]["params"]["num_designs"] = config.get("num_designs", 1)
            if method == "abmpnn":
                generator["data"]["params"]["sequences"] = config.get("num_designs", 1)
            elif method == "peptide_trap":
                generator["data"]["params"]["designs_per_register"] = config.get("num_designs", 1)
            elif method == "kejia_idr":
                generator["data"]["params"]["max_templates"] = config.get("num_designs", 1)
            if method == "boltzgen":
                generator["data"]["params"]["budget"] = config.get("num_designs", 1)
            _update(generator, method_settings.get(method, {}))
            if method == "rfd4":
                try:
                    options = json.loads(generator["data"]["params"]["options_json"])
                except (TypeError, ValueError) as exc:
                    raise ValueError("RFD4 task constraints must be a JSON object.") from exc
                if not isinstance(options, dict):
                    raise ValueError("RFD4 task constraints must be a JSON object.")
            expected_task = ("small_molecule_binder" if task_id == "ligand_conditioned_pocket" else
                             "enzyme" if task_id == "catalytic_motif_scaffolding" else "protein_binder")
            if method == "rfd4" and generator["data"]["params"]["task"] != expected_task:
                raise ValueError("RFD4 task settings must match the selected scientific objective.")
            if method == "complexa_generator" and generator["data"]["params"]["pipeline_task"] != "protein_binder":
                raise ValueError("Protein binder design requires Complexa's protein_binder task.")
            if method in {"rfdiffusion", "genie3", "pxdesign", "protpardelle"}:
                expected = "unconditional" if task_id == "new_protein_fold" else "binder"
                if generator["data"]["params"].get("task") != expected:
                    raise ValueError("Generator mode must match the scientific objective.")
            if binder_format == "vhh":
                required = {"bindcraft2_generator": ("modality", "VHH"),
                            "boltzgen_generator": ("task", "nanobody"), "esm_binder": ("binder_type", "vhh")}
                key, value = required[method]
                if generator["data"]["params"][key] != value:
                    raise ValueError("VHH method settings must retain the protected antibody mode.")
            # Target/collection inputs belong to the source; standalone generators own their inputs.
            _update(retained[0], inputs)
            chain(branch, retained)
            if method in BACKBONE:
                chain(branch, [_node("ligandmpnn", prefix + "-initial", {
                    "sequence_stage": "initial", "model_type": initial_model,
                    "sequences": config.get("sequences", 1), "fix_interface_residues": False,
                }, "Initial sequences · preserve native or design backbone")])
            if binder_format == "vhh" and method in VHH_METHODS:
                chain(branch, [_node("sabr", prefix + "-antibody-constraints")])
            if config.get("soluble_mpnn"):
                soluble = _node("ligandmpnn", prefix + "-soluble", label="Optional SolubleMPNN · retain parents and protect interface")
                _update(soluble, config.get("soluble_settings", {}))
                soluble["data"]["params"].update(sequence_stage="soluble_optimization", model_type="soluble_mpnn",
                                                  optimization_enabled=True, fix_interface_residues=True)
                chain(branch, [soluble])
            if config.get("sequence_qc") or task_id == "antibody_cdr_redesign":
                chain(branch, [_node("sequence_qc", prefix + "-qc")])
            # Native end-to-end methods need no automatic generic prediction/refinement.
            if method not in NATIVE or config.get("prediction") or config.get("soluble_mpnn"):
                predictor = config.get("prediction", "esmfold" if task_id == "new_protein_fold" else
                                       "af2" if method in {"kejia_idr", "peptide_trap", "rfdiffusion_beta_pairing"} else "boltz2")
                from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
                prediction_mode = "protein" if task_id == "new_protein_fold" else "protein_ligand" if task_id in {"ligand_conditioned_pocket", "catalytic_motif_scaffolding"} else "protein_protein"
                supported = NODE_CATALOG.get(predictor, {}).get("molecule_tasks", [])
                if prediction_mode not in supported and not (predictor == "esmfold" and prediction_mode == "protein"):
                    raise ValueError("Prediction engine does not support this scientific task's molecules.")
                prediction = _predictor(predictor, prefix + "-prediction", prediction_settings.get(predictor, {}))
                if predictor == "boltz2" and "ipsae" in analysis:
                    prediction["data"]["params"]["generate_ipsae"] = False
                chain(branch, [prediction])
            chain(branch, [_node(kind, prefix + "-" + kind) for kind in analysis])
        if config.get("selection") is not None:
            selection = _node("filter", prefix + "-selection", label="Optional candidate selection · declared evidence")
            _update(selection, config["selection"])
            chain(branch, [selection])
        if config.get("additional_validation"):
            if config.get("selection") is None:
                raise ValueError("Additional validation requires an explicit shortlist selection stage.")
            extra = config["additional_validation"]
            if not isinstance(extra, dict) or set(extra) != {"engine", "settings"}:
                raise ValueError("Additional validation needs engine and settings.")
            kind = extra["engine"]
            if task_id == "molecular_structure_prediction" and kind not in allowed:
                raise ValueError("Additional validation engine is incompatible with the molecular inputs.")
            if task_id != "molecular_structure_prediction":
                from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
                required_mode = "protein" if task_id == "new_protein_fold" else "protein_ligand" if TASKS[task_id].get("experimental") else "protein_protein"
                if required_mode not in NODE_CATALOG.get(kind, {}).get("molecule_tasks", []) and not (kind == "esmfold" and required_mode == "protein"):
                    raise ValueError("Additional validation engine is incompatible with this scientific objective.")
            chain(branch, [_predictor(kind, prefix + "-additional-validation", extra["settings"])])
        chain(branch, [_node("results", prefix + "-results", label="Results · " + method)])
        for position, item in enumerate(branch):
            item["position"] = {"x": 50 + position * 330, "y": 80 + index * 260}
        nodes.extend(branch)
    graph = {"schema_version": SCHEMA, "id": uuid.uuid4().hex, "name": TASKS[task_id]["label"],
             "nodes": nodes, "edges": edges, "viewport": {"x": 0, "y": 0, "zoom": .65},
             "template_id": task_id, "task_configuration": {"task_id": task_id, "version": 1, **config}}
    return validate_graph(graph)["graph"]


def workload(graph):
    """Report counts actually encoded by the graph; never fabricate wall-clock cost."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    summaries = []
    for node in graph["nodes"]:
        kind, params = node["data"]["kind"], node["data"]["params"]
        if "num_designs" in params:
            summaries.append(f"{node['data']['label']}: {params['num_designs']} attempt(s)")
        if kind in {"ligandmpnn", "foundry_mpnn"}:
            summaries.append(f"{node['data']['label']}: {params['sequences']} sequence(s) per eligible input")
        if kind == "abmpnn":
            summaries.append(f"AbMPNN: {params['sequences']} sequence variant(s) per parent")
        elif kind == "peptide_trap":
            summaries.append(f"Peptide Trap: up to {params['max_registers']} registers × {params['designs_per_register']} design(s)")
        elif kind == "kejia_idr":
            summaries.append(f"Kejia: up to {params['max_templates']} templates, {params['flex_docking_repeats']} flexible docking repeat(s) per register")
        if NODE_CATALOG[kind]["category"] == "Structure prediction":
            count = params.get("samples", params.get("diffusion_batch_size", params.get("model_count", 1)))
            summaries.append(f"{node['data']['label']}: {count} sample(s)/model(s) per input")
    return {"counts": summaries, "note": "Counts are not equal GPU time or native acceptance budgets. Prediction workload also scales with retained parents, variants and samples. Runtime depends on molecular size and method; no reliable time estimate is available."}
