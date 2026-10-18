"""Versioned visual pipeline definitions; independent of the UI and GPU runners."""
from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.runtime import app_home

SCHEMA = "molecular-design-studio.graph.v1"


def field(name, label, default, *, kind="number", minimum=None, maximum=None, options=None, section="optional", help=None):
    return {"name": name, "label": label, "default": default, "type": kind, "section": section, **({"help": help} if help else {}),
            **({"min": minimum} if minimum is not None else {}),
            **({"max": maximum} if maximum is not None else {}),
            **({"options": options} if options else {})}


def _spec(key, label, category, description, *, output="backbone", accepts=None, workflow=None, fields=(), color="#427f68", formats=("pdb", "cif")):
    return {"key": key, "label": label, "category": category, "description": description,
            "color": color, "workflow": workflow, "formats": list(formats),
            "inputs": [] if accepts is None else [{"id": "in", "label": "Target" if accepts == ["target"] else "Candidates", "accepts": accepts, "multiple": key == "results"}],
            "outputs": [] if output is None else [{"id": "out", "label": {"target": "Target", "backbone": "Backbones", "sequence": "Sequences", "prediction": "Predictions", "candidates": "Candidates"}.get(output, output), "type": output}],
            "fields": list(fields)}


COUNT = field("num_designs", "Design attempts", 10, minimum=1, maximum=10000)
LENGTH = field("binder_length", "Binder length", "60-100", kind="text")
STEPS = field("timesteps", "Diffusion steps", 50, minimum=1, maximum=1000)
SEED = field("seed", "Random seed", 0, minimum=0, maximum=2147483647)
SEQUENCES = field("sequences", "Sequences per backbone", 1, minimum=1, maximum=1000)
MSA = field("use_target_msa", "Use cached target MSA; fetch if missing", True, kind="boolean")
TEMPLATE = field("use_target_template", "Use target structure as template", True, kind="boolean")

NODE_CATALOG = {
    "target": _spec("target", "Target", "Input", "Structure, selected chains, hotspots and shared target context.", output="target", color="#758f84", fields=[
        field("target_pdb", "Target structure", "", kind="text"), field("target_chains", "Target chains", "A", kind="text"), field("hotspots", "Hotspot residues", "", kind="text")]),
    "design_source": _spec("design_source", "Successful project designs", "Input", "Selected filter-passing designs from a successful job in this project.", output="sequence", color="#758f84", fields=[field("source_run_dir", "Source run", "", kind="text")]),
    "generation_source": _spec("generation_source", "Generated project designs", "Input", "Select generated structures and sequences from a completed generator job in this project.", output="sequence", color="#758f84", formats=("pdb", "cif"), fields=[field("source_run_dir", "Generator run", "", kind="text")]),
    "rfd3": _spec("rfd3", "RFdiffusion3", "Design engine", "Generate backbones only; sequence design and refolding are separate nodes.", accepts=["target"], workflow="design.rfdiffusion3", fields=[COUNT, LENGTH, STEPS, SEED], formats=("cif",)),
    "rfdiffusion": _spec("rfdiffusion", "RFdiffusion", "Design engine", "Generate binder backbones without running the bundled downstream pipeline.", accepts=["target"], workflow="design.rfdiffusion_classic", fields=[COUNT, LENGTH, STEPS], formats=("pdb",)),
    "genie3": _spec("genie3", "Genie3", "Design engine", "Generation only; connect sequence design and prediction nodes.", accepts=["target"], workflow="design.genie3", fields=[COUNT, LENGTH, SEED], formats=("pdb",)),
    "pxdesign": _spec("pxdesign", "PXDesign generator", "Design engine", "PXDesign generation only, with downstream steps chosen on the canvas.", accepts=["target"], workflow="design.pxdesign", fields=[COUNT, LENGTH, field("n_steps", "Diffusion steps", 400, minimum=1, maximum=2000)], formats=("cif",)),
    "bindcraft2": _spec("bindcraft2", "BindCraft2", "Design pipeline", "Complete design, MPNN refinement, validation and filtering workflow.", accepts=["target"], output="prediction", workflow="design.bindcraft2", fields=[COUNT, LENGTH], color="#8268a2"),
    "bindcraft2_generator": _spec("bindcraft2_generator", "BindCraft2 Generator", "Design engine", "Run gradient design and emit only trajectories that pass the pre-ProteinMPNN gate. Connect ProteinMPNN and structure prediction as separate steps.", accepts=["target"], output="sequence", workflow="design.bindcraft2_generator", fields=[COUNT, LENGTH], formats=("cif",)),
    "boltzgen": _spec("boltzgen", "BoltzGen", "Design pipeline", "Generation, inverse folding, refolding, analysis and filtering.", accepts=["target"], output="prediction", workflow="design.boltzgen", fields=[COUNT, LENGTH, field("budget", "Final selection budget", 10, minimum=1, maximum=10000)], color="#8268a2", formats=("cif",)),
    "boltzgen_generator": _spec("boltzgen_generator", "BoltzGen Generator", "Design engine", "Run BoltzGen's design stage only. Connect sequence redesign and structure prediction as separate steps.", accepts=["target"], output="sequence", workflow="design.boltzgen_generator", fields=[COUNT, LENGTH, field("task", "Design task", "protein", kind="select", options=["protein", "nanobody", "peptide", "cyclotide", "redesign"]), field("sampling_steps", "Sampling steps", 20, minimum=1, maximum=1000), field("hotspots", "Binding-site residues", "", kind="text"), field("non_binding_residues", "Non-binding residues", "", kind="text")], formats=("cif",)),
    "pxdesign_pipeline": _spec("pxdesign_pipeline", "PXDesign pipeline", "Design pipeline", "Bundled ProteinMPNN and AF2-IG evaluation; Extended also uses Protenix.", accepts=["target"], output="prediction", workflow="design.pxdesign", fields=[COUNT, LENGTH, field("preset", "Evaluation preset", "preview", kind="select", options=["preview", "extended"])], color="#8268a2"),
    "genie3_pipeline": _spec("genie3_pipeline", "Genie3 pipeline", "Design pipeline", "Bundled generation, inverse folding and structure evaluation.", accepts=["target"], output="prediction", workflow="design.genie3", fields=[COUNT, LENGTH, SEED], color="#8268a2"),
    "protpardelle": _spec("protpardelle", "Protpardelle-1C", "Design pipeline", "Bundled sampling, sequence design and ESMFold evaluation.", accepts=["target"], output="prediction", workflow="design.protpardelle_1c", fields=[COUNT, LENGTH, SEQUENCES], color="#8268a2"),
    "complexa": _spec("complexa", "Proteina-Complexa", "Design pipeline", "Native generation, filtering, evaluation and analysis workflow.", accepts=["target"], output="prediction", workflow="design.proteina_complexa", fields=[COUNT, LENGTH, SEED], color="#8268a2"),
    "complexa_generator": _spec("complexa_generator", "Proteina-Complexa Generator", "Design engine", "Generate backbone and sequence. Choose a target here or connect a Target block, then attach sequence design or prediction steps.", accepts=["target"], output="sequence", workflow="design.proteina_complexa_generate", formats=("pdb",), fields=[
        field("pipeline_task", "Design task", "protein_binder", kind="select", options=["protein_binder","ligand_binder","motif_ligand"]),
        field("target_pdb", "Target structure", "", kind="text"), field("target_chains", "Target chains", "", kind="text"),
        field("hotspots", "Hotspots", "", kind="text"), field("ligand", "Ligand", "", kind="text"), field("motif_atom_spec", "Motif atoms", "", kind="text"),
        field("num_designs", "Number of designs", 1, minimum=1, maximum=16),field("binder_length", "Binder length", "70-150", kind="text"),SEED]),
    "esm_binder": _spec("esm_binder", "ESMFold2 Binder Generator", "Design engine", "Generate miniprotein, scFv, or VHH binder sequences against one target chain. ESMFold2 structure prediction is part of the gradient-guided design loop; sequences and predicted complexes remain available for downstream MPNN redesign or independent structure prediction.", accepts=["target"], output="sequence", workflow="design.esmfold2_binder", fields=[field("num_designs", "Number of designs", 1, minimum=1, maximum=128, section="basic"), field("binder_type", "Binder type", "miniprotein", kind="select", options=["miniprotein", "scfv", "vhh"], section="basic"), field("binder_framework", "Framework", "auto", kind="select", options=["auto", "trastuzumab_framework_vhvl", "atezolizumab_framework_vhvl", "ocankitug_framework_vhvl", "caplacizumab_vhh", "ozoralizumab_tnf_vhh", "vobarilizumab_il6r_vhh", "custom"], section="basic"), field("binder_template", "Custom VHH masked template", "", kind="text"), field("binder_length", "Miniprotein length", 80, minimum=40, maximum=200), field("optimization_steps", "Optimization steps", 150, minimum=1, maximum=2000), SEED], color="#427f68"),
    "foundry_mpnn": _spec("foundry_mpnn", "Foundry MPNN", "Sequence design", "Redesign binder chains while preserving the target sequence.", accepts=["backbone"], output="sequence", workflow="sequence.foundry_mpnn", fields=[SEQUENCES, field("model_type", "MPNN model", "protein_mpnn", kind="select", options=["protein_mpnn", "ligand_mpnn"])], color="#a88845", formats=("cif",)),
    "ligandmpnn": _spec("ligandmpnn", "ProteinMPNN / LigandMPNN", "Sequence design", "Sequence design on PDB or CIF candidate backbones; keep target fixed.", accepts=["backbone", "sequence", "prediction"], output="sequence", workflow="sequence.ligandmpnn", fields=[SEQUENCES, field("model_type", "MPNN model", "protein_mpnn", kind="select", options=["protein_mpnn", "ligand_mpnn", "soluble_mpnn"]), field("sampling_temp", "Sampling temperature", 0.1, minimum=0.000001, maximum=2)], color="#a88845"),
    "rf3": _spec("rf3", "RosettaFold3", "Structure prediction", "Refold designed target–binder complexes with Foundry RF3.", accepts=["sequence", "prediction"], output="prediction", workflow="refolding.rf3", fields=[MSA, TEMPLATE], color="#477e9b", formats=("cif",)),
    "boltz2": _spec("boltz2", "Boltz-2", "Structure prediction", "Refold target–binder complexes with an optional target template; target MSA is prepared automatically.", accepts=["sequence", "prediction"], output="prediction", workflow="refolding.boltz2_complex", fields=[TEMPLATE], color="#477e9b", formats=("cif",)),
    "af2": _spec("af2", "AlphaFold2 Initial Guess", "Structure prediction", "Evaluate complexes with a target template and initial guess.", accepts=["sequence", "prediction"], output="prediction", workflow="refolding.af2_initial_guess", fields=[field("num_recycles", "Recycles", 3, minimum=1, maximum=30)], color="#477e9b"),
    "esmfold": _spec("esmfold", "ESMFold · binder alone", "Structure prediction", "Predict each designed binder without the target.", accepts=["sequence", "prediction"], output="prediction", workflow="refolding.esmfold_monomer", fields=[], color="#477e9b"),
    "filter": _spec("filter", "Filter candidates", "Analysis", "Apply a score threshold, require passing filters, or retain the best candidates.", accepts=["backbone", "sequence", "prediction", "candidates"], output="candidates", color="#9a7854", fields=[
        field("passing_only", "Require an explicit filter pass", False, kind="boolean"), field("metric", "Score field", "", kind="text"), field("operator", "Comparison", ">=", kind="select", options=[">=", "<=", ">", "<"]), field("threshold", "Threshold", 0.7), field("limit", "Keep at most (0 = all)", 0, minimum=0, maximum=100000), field("higher_is_better", "Higher scores rank first", True, kind="boolean")]),
    "results": _spec("results", "Results & comparison", "Results", "Collect branch results for the studio's structure overlay and RMSD analysis.", accepts=["backbone", "sequence", "prediction", "candidates"], output=None, color="#5d7b6c"),
    "ipsae": _spec("ipsae", "ipSAE Interface Analysis", "Analysis", "Score predicted interfaces using the saved pairwise PAE matrix and complex structure.", accepts=["prediction"], output="prediction", workflow="analysis.ipsae", color="#9a7854", fields=[field("pae_cutoff", "PAE cutoff (Å)", 10., minimum=0., maximum=100.), field("dist_cutoff", "Interface distance cutoff (Å)", 10., minimum=0., maximum=100.)]),
    "rosetta_metrics": _spec("rosetta_metrics", "Rosetta Interface Analysis", "Analysis", "Calculate PyRosetta interface, packing, SASA, hydrogen-bond, and binder surface metrics.", accepts=["prediction"], output="prediction", workflow="analysis.rosetta", color="#9a7854", fields=[]),
}
NODE_CATALOG["complexa_generator"]["inputs"][0]["optional"] = True


# These are the same public workflow parameters used by the standalone input
# pages. Graph adapters supply structures, candidate collections and execution
# mode from connections; the form exposes the remaining engine controls.
_OPTIONAL_FIELDS = {
    "rfd3": [field("diffusion_batch_size", "Diffusion batch size", 1, minimum=1, maximum=32, section="basic"),
        field("step_scale", "Step scale", 1.5, minimum=1., maximum=5.),
        field("gamma_0", "Gamma 0", .6, minimum=0., maximum=1.),
        field("is_non_loopy", "Condition on non-loopy structures", True, kind="boolean"),
        field("classifier_free_guidance", "Classifier-free guidance", False, kind="boolean"),
        field("cfg_scale", "Guidance scale", 1.5, minimum=0., maximum=10.),
        field("low_memory_mode", "Low memory mode (experimental)", False, kind="boolean", help="The pinned Foundry checkpoint omits some chunked-only weights; this mode may reduce prediction quality.")],
    "rfdiffusion": [field("model_weights", "Model weights", "Complex_base", kind="text"),
        field("rfdiffusion_guidance_preset", "Guidance preset", "none", kind="text"),
        field("deterministic", "Deterministic sampling", False, kind="boolean"),
        field("save_trajectory", "Save diffusion trajectory", False, kind="boolean")],
    "genie3": [field("direction_scale", "Direction scale", 0., minimum=0.),
        field("compile_generation", "Compile generation model", False, kind="boolean")],
    "pxdesign": [field("dtype", "Precision", "bf16", kind="select", options=["bf16", "fp32"]),
        field("use_fast_ln", "Fast layer normalization", True, kind="boolean"),
        field("use_deepspeed_evo_attention", "DeepSpeed Evo attention", False, kind="boolean"),
        field("prepare_target_msa", "Prepare target MSA", True, kind="boolean")],
    "bindcraft2": [field("campaign_seed", "Campaign seed (0 = random)", 0, minimum=0, maximum=2147483647),
        field("redesign_interface", "Redesign interface", False, kind="boolean"),
        field("coldspots", "Coldspot residues", "", kind="text")],
    "boltzgen": [field("sampling_steps", "Sampling steps", 20, minimum=1, maximum=2000),
        field("omit_amino_acids", "Excluded amino acids", "", kind="text"),
        field("non_binding_residues", "Non-binding residues", "", kind="text"),
        field("filter_binding_site_adherence", "Filter binding-site adherence", False, kind="boolean")],
    "protpardelle": [field("step_scale", "Step scale", 1.2, minimum=0.),
        field("schurn", "Sampling churn", 200, minimum=0),
        field("crop_cond_start", "Crop conditioning start", 0., minimum=0., maximum=1.),
        field("batch_size", "Batch size", 1, minimum=1, maximum=128), SEED],
    "complexa": [field("n_steps", "Sampling steps", 400, minimum=50, maximum=1000),
        field("search_strategy", "Search strategy", "auto", kind="select", options=["auto", "best-of-n", "single-pass"]),
        field("guidance_weight", "Guidance weight", 1., minimum=0., maximum=1.),
        field("sampling_noise", "Sampling noise", .1, minimum=0., maximum=1.),
        field("deduplicate_sequences", "Remove duplicate sequences", True, kind="boolean"),
        field("replicas", "Replicas", 2, minimum=1, maximum=128),
        field("batch_size", "Batch size", 1, minimum=1, maximum=128)],
    "esm_binder": [field("learning_rate", "Learning rate", .1, minimum=.00001, maximum=10.),
        field("compile_model", "Compile model", False, kind="boolean"),
        field("checkpoint_lm", "Checkpoint language model", False, kind="boolean")],
    "foundry_mpnn": [field("number_of_batches", "Number of batches", 1, minimum=1, maximum=1000)],
    "ligandmpnn": [field("omit_aas", "Excluded amino acids", "CX", kind="text"),
        field("require_backbone_hotspot_filter_pass", "Require backbone hotspot filter pass", False, kind="boolean")],
    "rf3": [field("n_recycles", "Recycles", 10, minimum=1, maximum=100),
        field("num_steps", "Diffusion steps", 50, minimum=1, maximum=2000),
        field("diffusion_batch_size", "Number of samples", 5, minimum=1, maximum=100, section="basic"), SEED],
    "boltz2": [field("recycling_steps", "Recycling steps", 10, minimum=1, maximum=100),
        field("sampling_steps", "Sampling steps", 200, minimum=1, maximum=2000),
        field("diffusion_samples", "Number of samples", 3, minimum=1, maximum=100, section="basic"),
        field("write_full_pae", "Write full PAE", True, kind="boolean")],
    "af2": [field("model_count", "Number of models", 1, minimum=1, maximum=5, section="basic"),
        field("multimer", "Use multimer model", True, kind="boolean"),
        field("use_binder_template", "Use binder template", False, kind="boolean"),
        field("use_interface_template", "Use whole-complex template", False, kind="boolean"),
        field("msa_source", "MSA source", "single_sequence", kind="select", options=["single_sequence", "reference_then_colabfold"]),
        field("msa_max_depth", "Maximum MSA depth", 256, minimum=1, maximum=10000)],
    "esmfold": [field("min_plddt", "Minimum pLDDT", 70., minimum=0., maximum=100.),
        field("num_recycles", "Recycles", 4, minimum=1, maximum=100)],
}
_OPTIONAL_FIELDS["genie3_pipeline"] = _OPTIONAL_FIELDS["genie3"] + [
    field("inverse_folding_num_seq", "Sequences per backbone", 1, minimum=1, maximum=1000),
    field("folding_num_models", "Folding models", 5, minimum=1, maximum=5),
    field("folding_num_recycles", "Folding recycles", 20, minimum=1, maximum=100)]
_OPTIONAL_FIELDS["pxdesign_pipeline"] = copy.deepcopy(_OPTIONAL_FIELDS["pxdesign"])
_OPTIONAL_FIELDS["complexa_generator"] = [copy.deepcopy(f) for f in _OPTIONAL_FIELDS["complexa"] if f["name"] not in {"search_strategy","replicas"}]
_OPTIONAL_FIELDS["pxdesign"] = [f for f in _OPTIONAL_FIELDS["pxdesign"] if f["name"] != "prepare_target_msa"]
for _kind, _fields in _OPTIONAL_FIELDS.items():
    NODE_CATALOG[_kind]["fields"].extend(copy.deepcopy(_fields))
# Primary controls have the same simple form structure as the tool pages.
for _specification in NODE_CATALOG.values():
    for _field in _specification["fields"]:
        if _field["name"] in {"num_designs", "binder_length", "sequences", "model_type", "optimization_steps"}:
            _field["section"] = "basic"

# One tool contract and form for both saved MPNN kinds. Keep the old kind IDs
# so saved Foundry graphs remain editable without rewriting their history.
from mn_molecular_design_studio import mpnn_input
_mpnn_fields = [field("sequences", "Number of Sequences", 2, minimum=1, maximum=10000, section="basic")]
for _name, _default in mpnn_input.defaults().items():
    if _name == "num_seq_per_target":
        continue
    if _name == "model_type":
        _mpnn_fields.append(field(_name, "Model Type", _default, kind="select", options=list(mpnn_input.MODELS), section="basic"))
    elif _name == "backend":
        _mpnn_fields.append(field(_name, "Execution backend", _default, kind="select", options=["auto", "shared_ligandmpnn", "foundry"]))
    else:
        _default = 0 if _name == "seed" else "" if _default is None or isinstance(_default, dict) else _default
        _mpnn_fields.append(field(
            _name,
            "Fix interface residues from valid generator sequences" if _name == "fix_interface_residues" else _name.replace("_", " ").capitalize(),
            _default,
            kind="boolean" if isinstance(_default, bool) else "number" if isinstance(_default, (float, int)) else "text",
            help="Keep binder residues within 4 Å of the target fixed when a generator candidate has a valid sequence matching its structure; otherwise this setting is ignored." if _name == "fix_interface_residues" else None,
        ))
for _kind in ("ligandmpnn", "foundry_mpnn"):
    _specification = NODE_CATALOG[_kind]
    _specification.update(label="ProteinMPNN", workflow="sequence.mpnn", formats=["pdb", "cif"], fields=copy.deepcopy(_mpnn_fields))
    _specification["inputs"][0]["accepts"] = ["backbone", "sequence", "prediction"]
    _specification["inputs"][0]["multiple"] = True
NODE_CATALOG["foundry_mpnn"]["hidden"] = True
NODE_CATALOG["foundry_mpnn"]["fields"].append(field("number_of_batches", "Number of batches", 1, minimum=1, maximum=1000))

# Existing pipeline kinds retain their original workflow contract for saved
# versions. The four tool tiles and generator kinds have dedicated contracts.
from mn_molecular_design_studio import generator_input
NODE_CATALOG['protpardelle_pipeline']=copy.deepcopy(NODE_CATALOG['protpardelle'])
NODE_CATALOG['protpardelle_pipeline'].update(key='protpardelle_pipeline',label='Protpardelle-1C native pipeline')
for _legacy in ('protpardelle_pipeline','genie3_pipeline','pxdesign_pipeline'):
    NODE_CATALOG[_legacy]['hidden']=True
for _engine,_metadata in generator_input.ENGINES.items():
    NODE_CATALOG[_engine]=_spec(_engine,_metadata['title'],'Design engine',_metadata['description']+' Generation only.',
        accepts=['target'],output='backbone',workflow=generator_input.WORKFLOWS[_engine],
        fields=generator_input.canvas_fields(_engine),formats=('cif',) if _engine=='pxdesign' else ('pdb',))
    NODE_CATALOG[_engine]['inputs'][0]['optional']=True
    NODE_CATALOG[_engine]['inputs'][0]['label']='Target (optional)'
NODE_CATALOG['backbone_source']=_spec('backbone_source','Generated project backbones','Input',
    'Choose completed PXDesign, RFdiffusion, Genie3 or Protpardelle generation outputs in this project.',output='backbone',
    color='#758f84',fields=[field('source_run_dir','Generator run','',kind='text')])
NODE_CATALOG['protenix']=_spec('protenix','Protenix','Structure prediction','Evaluate designed complexes using the shared target MSA preparation.',
    accepts=['sequence','prediction'],output='prediction',workflow='refolding.protenix',color='#477e9b',formats=('cif',),fields=[
        field('use_msa','Use cached target MSA; fetch if missing',True,kind='boolean'),
        field('cycle','Recycles',10,minimum=1,maximum=100),field('diffusion_steps','Diffusion steps',200,minimum=1,maximum=2000),
        field('samples','Number of samples',1,minimum=1,maximum=100,section='basic')])

from mn_molecular_design_studio import inverse_folding_input
for _engine, _metadata in inverse_folding_input.ENGINES.items():
    NODE_CATALOG[_engine] = _spec(_engine, _metadata["title"], "Sequence design", _metadata["description"],
        accepts=["target", "backbone", "sequence", "prediction"], output="sequence", workflow=inverse_folding_input.WORKFLOWS[_engine],
        fields=inverse_folding_input.canvas_fields(_engine), color="#a88845")

TEMPLATE_REFERENCES={
    'pxdesign':('https://github.com/bytedance/PXDesign#-quick-start','Extended evaluation: generator → ProteinMPNN → AF2 Initial Guess and Protenix. The two prediction branches remain independently configurable.'),
    'rfdiffusion':('https://github.com/RosettaCommons/RFdiffusion#sequence-design-for-binders','ProteinMPNN sequence design followed by target-templated AF2 Initial Guess. The authors also permit ProteinMPNN without FastRelax; this template uses that option.'),
    'genie3':('https://github.com/aqlaboratory/genie3#application-3-binder-design','One ProteinMPNN sequence per backbone, then AF2 target-template evaluation with five models and twenty recycles. Studio uses its AF2 Initial Guess adapter; the authors use ColabFold template mode.'),
    'protpardelle':('https://github.com/ProteinDesignLab/protpardelle-1c#sampling','Structure generation, eight ProteinMPNN sequences per backbone, then ESMFold self-consistency evaluation. The starter task is unconditional generation; switch to a conditioned task in the generator form.'),
    'esmfold2_binder':('https://github.com/Biohub/esm/tree/main/cookbook/tutorials','ESMFold2 Binder Design performs gradient-guided sequence optimization with its built-in ESMFold2 and ESMC scoring. Its candidate output can be sent through ProteinMPNN and then any compatible independent structure predictor; edit or replace these downstream starter steps on the canvas.'),
    'all_generators_boltz2':('https://github.com/RosettaCommons/RFdiffusion','Fan every available design generator into one shared ProteinMPNN step, refold the redesigned complexes with Boltz-2, then calculate ipSAE and Rosetta interface metrics before collecting results.'),
}


def new_graph(template="foundry") -> dict[str, Any]:
    if template == "all_generators_boltz2":
        generator_kinds = ["rfd3", "rfdiffusion", "genie3", "pxdesign", "protpardelle", "bindcraft2_generator", "boltzgen_generator", "complexa_generator", "esm_binder"]
        kinds = ["target", *generator_kinds, "ligandmpnn", "boltz2", "ipsae", "rosetta_metrics", "results"]
        nodes = []
        for index, kind in enumerate(kinds):
            spec = NODE_CATALOG[kind]
            params = {f["name"]: copy.deepcopy(f["default"]) for f in spec["fields"]}
            if kind in generator_kinds:
                if "num_designs" in params:
                    params["num_designs"] = 1
                if kind in {"rfdiffusion", "genie3", "protpardelle"}:
                    params["task"] = "binder"
                if kind == "complexa_generator":
                    params["pipeline_task"] = "protein_binder"
            if kind == "ligandmpnn":
                params["sequences"] = 1
            if kind == "boltz2":
                params["write_full_pae"] = True
            nodes.append({
                "id": f"{kind}-{index}", "type": "studio",
                "position": {"x": 0, "y": 0},
                "data": {"kind": kind, "label": spec["label"], "params": params},
            })
        nodes[0]["position"] = {"x": 40, "y": 980}
        for i, kind in enumerate(generator_kinds):
            nodes[i + 1]["position"] = {"x": 390, "y": 50 + i * 235}
        downstream = nodes[len(generator_kinds) + 1:]
        for i, node in enumerate(downstream):
            node["position"] = {"x": 760 + i * 330, "y": 980}
        edges = []
        for i, node in enumerate(nodes[1:len(generator_kinds) + 1]):
            edges.append({"id": f"edge-target-{i}", "source": nodes[0]["id"], "sourceHandle": "out", "target": node["id"], "targetHandle": "in"})
            edges.append({"id": f"edge-mpnn-{i}", "source": node["id"], "sourceHandle": "out", "target": downstream[0]["id"], "targetHandle": "in"})
        for i, (source, dest) in enumerate(zip(downstream, downstream[1:])):
            edges.append({"id": f"edge-analysis-{i}", "source": source["id"], "sourceHandle": "out", "target": dest["id"], "targetHandle": "in"})
        result = {
            "schema_version": SCHEMA, "id": uuid.uuid4().hex,
            "name": "All generators → ProteinMPNN → Boltz-2 + interface analysis",
            "nodes": nodes, "edges": edges,
            "viewport": {"x": 0, "y": 0, "zoom": 0.42},
            "template_id": template,
        }
        return clean_graph(result)
    kinds = ["target", "rfd3", "foundry_mpnn", "rf3", "results"] if template == "foundry" else ["target", "results"]
    if template == "bindcraft2":
        kinds = ["target", "bindcraft2", "results"]
    elif template == "boltzgen":
        kinds = ["target", "boltzgen", "results"]
    elif template in {"complexa","complexa_generated"}:
        kinds = ["complexa_generator" if template=="complexa" else "generation_source","ligandmpnn","rf3","results"]
    elif template == "esmfold2_binder":
        kinds = ["target", "esm_binder", "ligandmpnn", "rf3", "results"]
    engine=template.removesuffix('_generated')
    if engine in generator_input.ENGINES:
        kinds=['backbone_source' if template.endswith('_generated') else engine,'ligandmpnn','esmfold' if engine=='protpardelle' else 'af2','results']
        if engine=='pxdesign':kinds.insert(-1,'protenix')
    nodes = []
    for index, kind in enumerate(kinds):
        spec = NODE_CATALOG[kind]
        positions = [(70, 70), (390, 70), (390, 290), (70, 290), (70, 510)]
        x, y = positions[index] if index<len(positions) else (710,290+220*(index-5))
        nodes.append({"id": f"{kind}-{index}", "type": "studio", "position": {"x": x, "y": y},
                      "data": {"kind": kind, "label": spec["label"], "params": {f["name"]: f["default"] for f in spec["fields"]}}})
    if template == "esmfold2_binder":
        for node, x in zip(nodes, (50, 360, 670, 980, 1290)):
            node["position"] = {"x": x, "y": 180}
    edges = [{"id": f"edge-{i}", "source": nodes[i]["id"], "sourceHandle": "out", "target": nodes[i+1]["id"], "targetHandle": "in"} for i in range(len(nodes)-1)] if template != "blank" else []
    if engine in generator_input.ENGINES:
        sequence=next(n for n in nodes if n['data']['kind']=='ligandmpnn')
        sequence['data']['params']['sequences']=8 if engine=='protpardelle' else 2 if engine=='rfdiffusion' else 1
        af2=next((n for n in nodes if n['data']['kind']=='af2'),None)
        if af2:
            af2['data']['params'].update(model_count=5 if engine=='genie3' else 1,num_recycles=20 if engine=='genie3' else 3,
                                       msa_source='single_sequence',use_binder_template=False,use_interface_template=False)
        if engine=='protpardelle' and not template.endswith('_generated'):
            defaults=generator_input.default_parameters(engine,'unconditional')
            nodes[0]['data']['params'].update({k:v for k,v in defaults.items() if k in nodes[0]['data']['params']})
            nodes[0]['data']['params']['target_pdb']='';nodes[0]['data']['params']['target_chains']=''
        if engine=='pxdesign':
            a,b,c,d,r=nodes
            a['position']={'x':50,'y':170};b['position']={'x':380,'y':170};c['position']={'x':710,'y':50};d['position']={'x':710,'y':300};r['position']={'x':1040,'y':170}
            pairs=[(a,b),(b,c),(b,d),(c,r),(d,r)]
            edges=[{'id':f'edge-{i}','source':a['id'],'sourceHandle':'out','target':b['id'],'targetHandle':'in'} for i,(a,b) in enumerate(pairs)]
    result={"schema_version": SCHEMA, "id": uuid.uuid4().hex, "name": {"foundry": "Foundry binder design", "bindcraft2": "BindCraft2 design", "boltzgen": "BoltzGen design", "complexa":"Proteina-Complexa custom pipeline", "complexa_generated":"Continue generated designs", "esmfold2_binder":"ESMFold2 Binder Design → ProteinMPNN → RosettaFold3"}.get(template, generator_input.ENGINES[engine]['title']+' · author-inspired pipeline' if engine in generator_input.ENGINES else "Untitled pipeline"), "nodes": nodes, "edges": edges, "viewport": {"x": 0, "y": 0, "zoom": 0.75}}
    if template in TEMPLATE_REFERENCES:result['template_id']=template
    return result


def clean_graph(graph: dict) -> dict:
    """Discard render state so saved definitions never contain execution results."""
    if not isinstance(graph, dict):
        raise ValueError("Pipeline graph must be an object.")
    graph = copy.deepcopy(graph)
    graph.setdefault("schema_version", SCHEMA)
    if not re.fullmatch(r"[a-f0-9]{32}", str(graph.get("id", ""))):
        raise ValueError("Pipeline ID must be a UUID hex string.")
    if graph.get("schema_version") != SCHEMA:
        raise ValueError("Unsupported pipeline graph schema.")
    if len(json.dumps(graph)) > 1_000_000:
        raise ValueError("Pipeline definition is too large.")
    nodes, edges = graph.get("nodes"), graph.get("edges")
    if not isinstance(nodes, list) or not isinstance(edges, list) or len(nodes) > 100 or len(edges) > 250:
        raise ValueError("A pipeline supports up to 100 nodes and 250 connections.")
    cleaned = []
    for node in nodes:
        if not isinstance(node, dict) or not re.fullmatch(r"[A-Za-z0-9_-]{1,100}", str(node.get("id", ""))):
            raise ValueError("Each node needs a unique, valid ID.")
        data = node.get("data") or {}
        kind = data.get("kind")
        if kind not in NODE_CATALOG:
            raise ValueError(f"Unknown pipeline node: {kind}.")
        params = data.get("params", {})
        if not isinstance(params, dict):
            raise ValueError("Node settings must be an object.")
        if kind=='protpardelle' and 'sequences' in params:
            kind='protpardelle_pipeline'
        allowed = {f["name"] for f in NODE_CATALOG[kind]["fields"]} | ({"selected_candidate_ids"} if kind in {"design_source","generation_source","backbone_source"} else set())
        # Migrate generation-only kinds from the previous small canvas forms.
        if kind=='rfdiffusion':params.pop('rfdiffusion_guidance_preset',None)
        if kind=='pxdesign':
            params.pop('use_deepspeed_evo_attention',None)
            params.pop('prepare_target_msa',None)
        if set(params) - allowed:
            raise ValueError(f"Unsupported settings for {NODE_CATALOG[kind]['label']}: {', '.join(sorted(set(params)-allowed))}.")
        if "selected_candidate_ids" in params and (not isinstance(params["selected_candidate_ids"], list) or not all(isinstance(v, str) for v in params["selected_candidate_ids"])):
            raise ValueError("Selected designs must be a list of candidate IDs.")
        for f in NODE_CATALOG[kind]["fields"]:
            value = params.setdefault(f["name"], f["default"])
            if f["type"] == "number":
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                    raise ValueError(f"{f['label']} must be a finite number.")
                if "min" in f and value < f["min"] or "max" in f and value > f["max"]:
                    raise ValueError(f"{f['label']} is outside its supported range.")
                if isinstance(f["default"], int) and not float(value).is_integer():
                    raise ValueError(f"{f['label']} must be a whole number.")
            elif f["type"] == "boolean" and type(value) is not bool:
                raise ValueError(f"{f['label']} must be true or false.")
            elif f["type"] in {"text", "select"} and not isinstance(value, str):
                raise ValueError(f"{f['label']} must be text.")
            elif f['type']=='json' and not isinstance(value,list):
                raise ValueError(f"{f['label']} must be a list.")
            if f["type"] == "select" and value not in f["options"]:
                raise ValueError(f"Unsupported {f['label']}.")
        pos = node.get("position") or {}
        position = {axis: float(pos.get(axis, 0)) for axis in ("x", "y")}
        if not all(math.isfinite(v) and abs(v) < 1e7 for v in position.values()):
            raise ValueError("Node positions must be finite.")
        cleaned.append({"id": node["id"], "type": "studio", "position": position, "data": {"kind": kind, "label": str(data.get("label") or NODE_CATALOG[kind]["label"])[:100], "params": params}})
    viewport = graph.get("viewport") or {"x": 0, "y": 0, "zoom": 1}
    if not isinstance(viewport, dict) or any(isinstance(viewport.get(k), bool) or not isinstance(viewport.get(k), (float, int)) or not math.isfinite(viewport[k]) for k in ("x", "y", "zoom")) or not .05 <= viewport["zoom"] <= 5:
        raise ValueError("Invalid canvas viewport.")
    return {"schema_version": SCHEMA, "id": graph["id"], "name": str(graph.get("name") or "Untitled pipeline")[:160], "nodes": cleaned,
            "edges": [{k: edge[k] for k in ("id", "source", "target", "sourceHandle", "targetHandle") if k in edge} for edge in edges if isinstance(edge, dict)],
            "viewport": viewport,**({'template_id':graph['template_id']} if graph.get('template_id') in TEMPLATE_REFERENCES else {})}


def validate_graph(graph: dict, *, require_inputs=False) -> dict:
    graph = clean_graph(graph)
    nodes = {n["id"]: n for n in graph["nodes"]}
    if len(nodes) != len(graph["nodes"]):
        raise ValueError("Node IDs must be unique.")
    if not nodes:
        raise ValueError("Add at least one input and a Results node.")
    incoming = {n: [] for n in nodes}
    outgoing = {n: [] for n in nodes}
    edge_ids = set()
    for edge in graph["edges"]:
        source, target = edge.get("source"), edge.get("target")
        if source not in nodes or target not in nodes or source == target:
            raise ValueError("Connections must join two different existing nodes.")
        if not edge.get("id") or edge["id"] in edge_ids:
            raise ValueError("Connection IDs must be unique.")
        edge_ids.add(edge["id"])
        a, b = (NODE_CATALOG[nodes[i]["data"]["kind"]] for i in (source, target))
        if not a["outputs"] or not b["inputs"] or edge.get("sourceHandle", "out") != "out" or edge.get("targetHandle", "in") != "in":
            raise ValueError("Connection references an unavailable socket.")
        output_type = a["outputs"][0]["type"]
        accepts = b["inputs"][0]["accepts"]
        if output_type == "candidates" and accepts == ["target"] or output_type != "candidates" and output_type not in accepts:
            raise ValueError(f"{a['label']} outputs cannot connect to {b['label']}.")
        incoming[target].append(source)
        outgoing[source].append(target)
        if len(incoming[target]) > 1 and not b["inputs"][0]["multiple"]:
            raise ValueError(f"{b['label']} accepts one upstream connection. Branch after its output instead.")
    queue = sorted(n for n in nodes if not incoming[n])
    order, degrees = [], {n: len(incoming[n]) for n in nodes}
    while queue:
        node_id = queue.pop(0); order.append(node_id)
        for dest in outgoing[node_id]:
            degrees[dest] -= 1
            if degrees[dest] == 0:
                queue.append(dest)
    if len(order) != len(nodes):
        raise ValueError("A pipeline cannot contain a cycle. Connect steps in the order data should flow.")
    output_types = {}
    output_formats = {}
    ligand_context = {}
    for node_id in order:
        kind = nodes[node_id]["data"]["kind"]
        spec = NODE_CATALOG[kind]
        ligand_context[node_id] = (nodes[node_id]['data']['params'].get('pipeline_task') in {'ligand_binder','motif_ligand'}
                                   if kind == 'complexa_generator'
                                   else any(ligand_context[source] for source in incoming[node_id]))
        if kind == 'af2' and ligand_context[node_id]:
            raise ValueError('AlphaFold2 does not predict ligand complexes. Use RosettaFold3, or ESMFold for the binder alone.')
        for source in incoming[node_id]:
            if output_types[source] not in spec["inputs"][0]["accepts"]:
                raise ValueError(f"{nodes[source]['data']['label']} does not supply the data required by {spec['label']}.")
        output_types[node_id] = (output_types[incoming[node_id][0]] if kind == "filter" and incoming[node_id] else spec["outputs"][0]["type"] if spec["outputs"] else "results")
        output_formats[node_id] = output_formats[incoming[node_id][0]] if kind == "filter" and incoming[node_id] else spec["formats"]
    for node_id, node in nodes.items():
        kind = node["data"]["kind"]
        if NODE_CATALOG[kind]["inputs"] and not incoming[node_id] and not NODE_CATALOG[kind]["inputs"][0].get("optional"):
            raise ValueError(f"Connect an input to {node['data']['label']}.")
        if kind != "results" and not outgoing[node_id]:
            raise ValueError(f"Connect {node['data']['label']} to a downstream step or Results node.")
        if require_inputs and kind == "target":
            if not node["data"]["params"].get("target_pdb"):
                raise ValueError("Load a target structure before running the pipeline.")
        if require_inputs and kind == "complexa_generator" and not incoming[node_id] and not node["data"]["params"].get("target_pdb"):
            raise ValueError("Load a structure in Proteina-Complexa Generator or connect a Target block.")
        length = node["data"]["params"].get("binder_length")
        if isinstance(length, str):
            match = re.fullmatch(r"(\d+)(?:-(\d+))?", length)
            if not match or not 1 <= int(match[1]) <= int(match[2] or match[1]) <= 10000:
                raise ValueError("Binder length must be a positive length or range, such as 74 or 60-100.")
        if require_inputs and kind == "design_source" and not node["data"]["params"].get("source_run_dir"):
            raise ValueError("Select a successful design job from this project.")
        if require_inputs and kind == "generation_source" and not node["data"]["params"].get("source_run_dir"):
            raise ValueError("Select a completed generator job from this project.")
        if require_inputs and kind=='backbone_source' and not node['data']['params'].get('source_run_dir'):
            raise ValueError('Select a completed generator job from this project.')
        if require_inputs and kind in generator_input.ENGINES and not incoming[node_id]:
            p=node['data']['params']
            if p.get('task') not in {'unconditional','symmetry','custom','fold'} and not p.get('target_pdb'):
                raise ValueError('Load a structure in '+node['data']['label']+' or connect a Target block.')
    if not any(n["data"]["kind"] == "results" for n in nodes.values()):
        raise ValueError("Add a Results node to collect the pipeline outputs.")
    return {"graph": graph, "order": order, "incoming": incoming, "outgoing": outgoing}


def fingerprint(kind: str, params: dict, upstream: list[str]) -> str:
    return hashlib.sha256(json.dumps({"schema": SCHEMA, "kind": kind, "params": params, "upstream": sorted(upstream)}, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def scope_root(owner_id: str, project_id: str) -> Path:
    scope = hashlib.sha256(f"{owner_id}:{project_id}".encode()).hexdigest()
    return app_home() / "pipelines" / scope


def save_definition(graph: dict, owner_id: str, project_id: str, *, version=False) -> dict:
    graph = clean_graph(graph)
    root = scope_root(owner_id, project_id) / graph["id"]
    root.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(json.dumps(graph, sort_keys=True).encode()).hexdigest()
    payload = {"graph": graph, "revision": digest, "saved_at": datetime.now(timezone.utc).isoformat()}
    draft = root / "draft.json"
    temporary = root / f".draft-{uuid.uuid4().hex}.tmp"
    temporary.write_text(json.dumps(payload, indent=2)); temporary.replace(draft)
    if version:
        versions = root / "versions"; versions.mkdir(exist_ok=True)
        try:
            with (versions / f"{digest}.json").open("x") as handle:
                json.dump(payload, handle, indent=2)
        except FileExistsError:
            pass
    return payload


def list_definitions(owner_id: str, project_id: str) -> list[dict]:
    definitions = []
    for path in scope_root(owner_id, project_id).glob("*/draft.json"):
        try:
            d = json.loads(path.read_text()); d["versions"] = len(list((path.parent / "versions").glob("*.json"))); definitions.append(d)
        except (OSError, ValueError):
            continue
    return sorted(definitions, key=lambda d: d["saved_at"], reverse=True)


def list_versions(owner_id: str, project_id: str, pipeline_id: str) -> list[dict]:
    if not re.fullmatch(r"[a-f0-9]{32}", pipeline_id):
        raise ValueError("Invalid pipeline ID.")
    versions = []
    for path in (scope_root(owner_id, project_id) / pipeline_id / "versions").glob("*.json"):
        try:
            versions.append(json.loads(path.read_text()))
        except (OSError, ValueError):
            continue
    return sorted(versions, key=lambda d: d["saved_at"], reverse=True)


WORKFLOW_ROLES = {
    **{workflow: "Sequence design engine" for workflow in inverse_folding_input.WORKFLOWS.values()},
    **{workflow:'Design engine' for workflow in generator_input.WORKFLOWS.values()},
    "design.rfdiffusion3": "Design engine",
    "design.rfdiffusion_classic": "Design engine · pipeline available",
    "design.bindcraft2": "Design pipeline", "design.boltzgen": "Design pipeline",
    "design.bindcraft2_generator": "Design engine", "design.boltzgen_generator": "Design engine",
    "design.pxdesign": "Design pipeline", "design.proteina_complexa": "Design pipeline",
    "design.proteina_complexa_generate": "Design engine",
    "design.genie3": "Design engine · pipeline available",
    "design.protpardelle_1c": "Design engine · pipeline available",
    "design.esmfold2_binder": "Design engine · downstream pipeline available",
    "sequence.foundry_mpnn": "Sequence design engine", "sequence.ligandmpnn": "Sequence design engine",
    "analysis.ipsae": "Result analysis", "analysis.rosetta": "Result analysis",
}
