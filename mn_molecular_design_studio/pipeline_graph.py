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


def _spec(key, label, category, description, *, output="backbone", accepts=None, workflow=None, fields=(), color="#427f68", formats=("pdb", "cif"), coldspot_support=None):
    return {"key": key, "label": label, "category": category, "description": description,
            "color": color, "workflow": workflow, "formats": list(formats),
            "inputs": [] if accepts is None else [{"id": "in", "label": "Target" if accepts == ["target"] else "Results" if key == "results" else "Candidates", "accepts": accepts, "multiple": key == "results"}],
            "outputs": [] if output is None else [{"id": "out", "label": {"target": "Target", "backbone": "Backbones", "sequence": "Sequences", "prediction": "Predictions", "candidates": "Candidates", "analysis": "Surface results"}.get(output, output), "type": output}],
            "fields": list(fields), **({"coldspot_support": coldspot_support} if coldspot_support else {})}


COUNT = field("num_designs", "Design attempts", 10, minimum=1, maximum=10000)
LENGTH = field("binder_length", "Binder length", "60-100", kind="text")
STEPS = field("timesteps", "Diffusion steps", 50, minimum=1, maximum=1000)
SEED = field("seed", "Random seed", 0, minimum=0, maximum=2147483647)
SEQUENCES = field("sequences", "Sequences per backbone", 1, minimum=1, maximum=1000)
MSA = field("use_target_msa", "Use cached target MSA; fetch if missing", True, kind="boolean")
TEMPLATE = field("use_target_template", "Use target structure as template", True, kind="boolean")

NODE_CATALOG = {
    "cofold_input": _spec("cofold_input", "Cofolding Inputs", "Cofolding input", "Enter a protein target and any partner molecules for the selected prediction task.", output="molecule_task", color="#758f84", fields=[
        field("task_type", "Prediction task", "protein_protein", kind="select", options=["protein", "protein_protein", "protein_ligand", "protein_dna", "protein_rna", "protein_macrocyclic"], section="basic"),
        field("target_protein_sequences", "Target protein sequence(s)", "", kind="text", section="basic", help="Enter one sequence per target chain. Separate chains with a colon or new line."),
        field("partner_protein_sequences", "Protein partner sequence(s)", "", kind="text", section="basic", help="Used for protein–protein and protein–macrocycle tasks. Separate partner chains with a colon or new line."),
        field("dna_sequences", "DNA sequence(s)", "", kind="text", section="basic", help="Enter DNA chains when using the protein–DNA task."),
        field("rna_sequences", "RNA sequence(s)", "", kind="text", help="Enter RNA chains when using the protein–RNA task."),
        field("ligand_input_type", "Ligand notation", "smiles", kind="select", options=["smiles", "ccd"], help="SMILES works across the broadest set of ligand-capable models; CCD-only inputs exclude models that accept SMILES only."),
        field("ligand_ccd_codes", "Ligand CCD code(s)", "", kind="text", help="One or more CCD component IDs, separated by colons or new lines."),
        field("ligand_smiles", "Ligand SMILES", "", kind="text", help="Enter one SMILES string per line."),
        field("cyclic_partner", "Treat the protein partner as cyclic", False, kind="boolean", section="basic", help="Only engines with native cyclic-protein support are included in the macrocyclic template."),
        field("target_pdb", "Optional target PDB path", "", kind="text", help="Required by AlphaFold2 Initial Guess and BoltzGen Fold; optional as a target template for engines that support templates."),
        field("target_chain_id", "Target PDB chain ID", "A", kind="text", help="Single chain selected from the optional target PDB. Used by AlphaFold2 Initial Guess and BoltzGen Fold."),
    ]),
    "target": _spec("target", "Target Structure", "Target inputs", "Shared PDB or mmCIF structure, selected protein chains and optional baseline interface annotations. Add Binder Design Context for branch-specific hotspot or coldspot settings.", output="target", color="#758f84", fields=[
        field("target_pdb", "Target structure", "", kind="text"), field("target_chains", "Target chains", "A", kind="text"), field("hotspots", "Hotspot residues", "", kind="text"), field("coldspots", "Coldspot residues", "", kind="text")]),
    "binder_context": _spec("binder_context", "Binder Design Context", "Target context", "Set or override binder hotspots and coldspots on the connected target. The same structure continues through this branch, while other branches stay unchanged.", accepts=["target"], output="target", color="#668d75", fields=[
        field("hotspots", "Hotspot override", "", kind="text", help="Leave empty to inherit the connected Target Structure hotspots."),
        field("coldspots", "Coldspot override", "", kind="text", help="Leave empty to inherit the connected Target Structure coldspots."),
        field("clear_unselected", "Clear inherited annotations when overrides are empty", False, kind="boolean", help="Enable this to explicitly pass no hotspots or coldspots when the corresponding override is blank.")]),
    "design_source": _spec("design_source", "Successful project designs", "Input", "Selected filter-passing designs from a successful job in this project.", output="sequence", color="#758f84", fields=[field("source_run_dir", "Source run", "", kind="text")]),
    "generation_source": _spec("generation_source", "Generated project designs", "Input", "Select generated structures and sequences from a completed generator job in this project.", output="sequence", color="#758f84", formats=("pdb", "cif"), fields=[field("source_run_dir", "Generator run", "", kind="text")]),
    "rfd3": _spec("rfd3", "RFdiffusion3", "Design engine", "Generate backbones only; sequence design and refolding are separate nodes.", accepts=["target"], workflow="design.rfdiffusion3", fields=[COUNT, LENGTH, STEPS, SEED], formats=("cif",)),
    "rfdiffusion": _spec("rfdiffusion", "RFdiffusion", "Design engine", "Generate binder backbones without running the bundled downstream pipeline.", accepts=["target"], workflow="design.rfdiffusion_classic", fields=[COUNT, LENGTH, STEPS], formats=("pdb",)),
    "genie3": _spec("genie3", "Genie3", "Design engine", "Generation only; connect sequence design and prediction nodes.", accepts=["target"], workflow="design.genie3", fields=[COUNT, LENGTH, SEED], formats=("pdb",)),
    "pxdesign": _spec("pxdesign", "PXDesign generator", "Design engine", "PXDesign generation only, with downstream steps chosen on the canvas.", accepts=["target"], workflow="design.pxdesign", fields=[COUNT, LENGTH, field("n_steps", "Diffusion steps", 400, minimum=1, maximum=2000)], formats=("cif",)),
    "bindcraft2": _spec("bindcraft2", "BindCraft2", "Design pipeline", "Complete design, MPNN refinement, validation and filtering workflow.", accepts=["target"], output="prediction", workflow="design.bindcraft2", fields=[COUNT, LENGTH], color="#8268a2"),
    "bindcraft2_generator": _spec("bindcraft2_generator", "BindCraft2 Generator", "Design engine", "Run gradient design and emit only trajectories that pass the pre-ProteinMPNN gate. Select a binder format such as VHH; connect ProteinMPNN and structure prediction as separate steps.", accepts=["target"], output="sequence", workflow="design.bindcraft2_generator", fields=[COUNT, LENGTH, field("modality", "Binder format", "binder", kind="select", options=["binder", "large_binder", "peptide", "cyclic_peptide", "VHH", "scFv", "Fab", "ARP", "multidomain", "homo_oligomer"], section="basic"), field("scaffold_edits", "VHH scaffold edit regions", "", kind="text", section="basic", help="Only used for VHH, scFv, Fab, or ARP scaffold modes. For VHH, the template supplies the default CDR edit regions.")], formats=("cif",)),
    "boltzgen": _spec("boltzgen", "BoltzGen", "Design pipeline", "Generation, inverse folding, refolding, analysis and filtering.", accepts=["target"], output="prediction", workflow="design.boltzgen", fields=[COUNT, LENGTH, field("budget", "Final selection budget", 10, minimum=1, maximum=10000)], color="#8268a2", formats=("cif",)),
    "boltzgen_generator": _spec("boltzgen_generator", "BoltzGen Generator", "Design engine", "Run BoltzGen's design stage only. Connect sequence redesign and structure prediction as separate steps.", accepts=["target"], output="sequence", workflow="design.boltzgen_generator", fields=[COUNT, LENGTH, field("task", "Design task", "protein", kind="select", options=["protein", "nanobody", "peptide", "cyclotide", "redesign"]), field("sampling_steps", "Sampling steps", 20, minimum=1, maximum=1000), field("hotspots", "Binding-site residues", "", kind="text"), field("non_binding_residues", "Non-binding residues override", "", kind="text", help="Leave empty to inherit coldspots from the connected Target block. BoltzGen maps these residues to its chain-local numbering.")], formats=("cif",)),
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
    "rf3": _spec("rf3", "RosettaFold3", "Structure prediction", "Refold designed target–binder complexes with Foundry RF3.", accepts=["sequence", "prediction", "molecule_task"], output="prediction", workflow="refolding.rf3", fields=[MSA, TEMPLATE], color="#477e9b", formats=("cif",)),
    "boltz2": _spec("boltz2", "Boltz-2", "Structure prediction", "Refold target–binder complexes with an optional target template; target MSA is prepared automatically.", accepts=["sequence", "prediction", "molecule_task"], output="prediction", workflow="refolding.boltz2_complex", fields=[TEMPLATE], color="#477e9b", formats=("cif",)),
    "af3": _spec("af3", "AlphaFold3", "Structure prediction", "Refold target–binder candidates with MN Cofolding's AlphaFold3-family JAX runner. Select the official AF3 weights or another supported AF3-family parameter set.", accepts=["sequence", "prediction", "molecule_task"], output="prediction", workflow="refolding.optimized_prediction", fields=[
        field("model", "AF3-family parameter set", "alphafold3", kind="select", options=["alphafold3", "openfold3", "openbind0", "boltz2", "protenix2", "rosettafold3", "intellifold2", "opendde"], section="basic"),
        field("optimization_mode", "Optimization mode", "off", kind="select", options=["off", "exact", "fast", "big"], section="basic"),
        field("samples", "Samples per candidate", 1, minimum=1, maximum=5, section="basic"),
        field("num_recycles", "Recycles", 10, minimum=1, maximum=100),
        field("target_use_msa", "Use cached target MSA; fetch if missing", True, kind="boolean"),
        field("use_target_template", "Use target structure as template", True, kind="boolean"),
        field("accept_alphafold3_terms", "Accept official AlphaFold3 weight terms", False, kind="boolean",
              help="Required only when the official AlphaFold3 weights are selected. Review those terms before enabling them."),
    ], color="#477e9b", formats=("cif",)),
    "af2": _spec("af2", "AlphaFold2 Initial Guess", "Structure prediction", "Evaluate complexes with a target template and initial guess.", accepts=["sequence", "prediction", "molecule_task"], output="prediction", workflow="refolding.af2_initial_guess", fields=[field("num_recycles", "Recycles", 3, minimum=1, maximum=30), field("model_count", "Models", 1, minimum=1, maximum=5), field("multimer", "Use multimer model", True, kind="boolean")], color="#477e9b"),
    "esmfold": _spec("esmfold", "ESMFold · binder alone", "Structure prediction", "Predict each designed binder without the target.", accepts=["sequence", "prediction"], output="prediction", workflow="refolding.esmfold_monomer", fields=[], color="#477e9b"),
    "filter": _spec("filter", "Filter candidates", "Analysis", "Apply a score threshold, require passing filters, or retain the best candidates.", accepts=["backbone", "sequence", "prediction", "candidates"], output="candidates", color="#9a7854", fields=[
        field("passing_only", "Require an explicit filter pass", False, kind="boolean"), field("metric", "Score field", "", kind="text"), field("operator", "Comparison", ">=", kind="select", options=[">=", "<=", ">", "<"]), field("threshold", "Threshold", 0.7), field("limit", "Keep at most (0 = all)", 0, minimum=0, maximum=100000), field("higher_is_better", "Higher scores rank first", True, kind="boolean")]),
    "surface_analysis": _spec("surface_analysis", "Protein Surface Analysis", "Analysis", "Run PeSTo, ScanNet, Surf2Spot, and MaSIF-seed against one target and compare predicted PPI and hotspot residues.", accepts=["target"], output="analysis", workflow="detection.surface_analysis", color="#397b5b", formats=("pdb",), fields=[
        field("pesto_interface", "PeSTo interface class", "protein", kind="select", options=["protein", "dna_rna", "ion", "ligand", "lipid"]),
        field("scannet_mode", "ScanNet prediction mode", "interface", kind="select", options=["interface", "epitope", "idp"]),
        field("scannet_use_msa", "Use MSA with ScanNet", False, kind="boolean"),
        field("masif_chain_id", "MaSIF-seed chain", "", kind="text")]),
    "results": _spec("results", "Results & comparison", "Results", "Collect candidate structures, metrics, and surface-analysis runs from connected branches.", accepts=["backbone", "sequence", "prediction", "candidates", "analysis"], output=None, color="#5d7b6c"),
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
        field("coldspots", "Coldspot override", "", kind="text", help="Leave empty to inherit coldspots from the connected Target block.")],
    "boltzgen": [field("sampling_steps", "Sampling steps", 20, minimum=1, maximum=2000),
        field("omit_amino_acids", "Excluded amino acids", "", kind="text"),
        field("non_binding_residues", "Non-binding residues override", "", kind="text", help="Leave empty to inherit coldspots from the connected Target block. BoltzGen maps these residues to its chain-local numbering."),
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
# Coldspots are target-wide context, but only BindCraft2 and BoltzGen expose
# matching engine inputs. The canvas shows this support state on every generator.
for _kind in ('rfd3','rfdiffusion','genie3','pxdesign','pxdesign_pipeline','genie3_pipeline',
              'protpardelle','protpardelle_pipeline','complexa','complexa_generator','esm_binder'):
    NODE_CATALOG[_kind]['coldspot_support'] = 'unsupported'
for _kind in ('bindcraft2','bindcraft2_generator'):
    NODE_CATALOG[_kind]['coldspot_support'] = 'native'
for _kind in ('boltzgen','boltzgen_generator'):
    NODE_CATALOG[_kind]['coldspot_support'] = 'mapped'
for _kind in ('rfd3', 'rfdiffusion', 'genie3', 'pxdesign', 'pxdesign_pipeline',
              'genie3_pipeline', 'bindcraft2', 'bindcraft2_generator', 'boltzgen',
              'boltzgen_generator', 'complexa'):
    NODE_CATALOG[_kind]['hotspot_support'] = 'native'
NODE_CATALOG['complexa_generator']['hotspot_support'] = 'conditional'
for _kind in ('protpardelle', 'protpardelle_pipeline', 'esm_binder'):
    NODE_CATALOG[_kind]['hotspot_support'] = 'unsupported'
NODE_CATALOG['surface_analysis']['hotspot_support'] = 'unsupported'
NODE_CATALOG['surface_analysis']['coldspot_support'] = 'unsupported'
NODE_CATALOG['surface_analysis']['target_context_support'] = 'structure_only'


def _uses_binder_context(kind: str) -> bool:
    """Whether a node can use at least one shared binder-site annotation."""
    spec = NODE_CATALOG.get(kind, {})
    return (
        spec.get('hotspot_support') in {'native', 'conditional'}
        or spec.get('coldspot_support') in {'native', 'mapped'}
    )


NODE_CATALOG['backbone_source']=_spec('backbone_source','Generated project backbones','Input',
    'Choose completed PXDesign, RFdiffusion, Genie3 or Protpardelle generation outputs in this project.',output='backbone',
    color='#758f84',fields=[field('source_run_dir','Generator run','',kind='text')])
NODE_CATALOG['protenix']=_spec('protenix','Protenix','Structure prediction','Evaluate designed complexes using the shared target MSA preparation.',
    accepts=['sequence','prediction'],output='prediction',workflow='refolding.protenix',color='#477e9b',formats=('cif',),fields=[
        field('use_msa','Use cached target MSA; fetch if missing',True,kind='boolean'),
        field('cycle','Recycles',10,minimum=1,maximum=100),field('diffusion_steps','Diffusion steps',200,minimum=1,maximum=2000),
        field('samples','Number of samples',1,minimum=1,maximum=100,section='basic')])
NODE_CATALOG['openfold3_legacy']=_spec(
    'openfold3_legacy','OpenFold3 · P2 weights','Structure prediction',
    'Run the standalone OpenFold3 0.5 P2 workflow on connected candidate complexes.',
    accepts=['sequence','prediction'],output='prediction',workflow='refolding.openfold3',color='#477e9b',formats=('cif',),fields=[
        field('use_target_msa','Use cached target MSA; fetch if missing',True,kind='boolean',section='basic'),
        field('num_diffusion_samples','Diffusion samples',5,minimum=1,maximum=100,section='basic'),
        field('num_model_seeds','Model seeds',1,minimum=1,maximum=10),
        field('num_recycles','Recycles',3,minimum=1,maximum=100),
        field('use_msa_server','Allow MSA server fallback',False,kind='boolean')])
NODE_CATALOG['protenix_cli']=_spec(
    'protenix_cli','Protenix 2.0 · CLI','Structure prediction',
    'Run the standalone Protenix 2.0 CLI workflow on connected candidate complexes.',
    accepts=['sequence','prediction'],output='prediction',workflow='refolding.protenix_cli',color='#477e9b',formats=('cif',),fields=[
        field('use_msa','Use cached target MSA; fetch if missing',True,kind='boolean',section='basic'),
        field('use_template','Use templates',False,kind='boolean',section='basic'),
        field('cycle','Recycles',10,minimum=1,maximum=100),
        field('diffusion_steps','Diffusion steps',200,minimum=1,maximum=1000),
        field('samples','Samples per candidate',1,minimum=1,maximum=100,section='basic'),
        field('seeds','Random seeds', '101',kind='text')])
NODE_CATALOG['boltzgen_fold_batch']=_spec(
    'boltzgen_fold_batch','BoltzGen Fold · candidate refolding','Structure prediction',
    'Run the standalone BoltzGen target-template folding workflow on connected candidate complexes.',
    accepts=['sequence','prediction'],output='prediction',workflow='refolding.boltzgen_fold',color='#477e9b',formats=('cif',),fields=[
        field('recycling_steps','Recycling steps',3,minimum=1,maximum=100),
        field('sampling_steps','Sampling steps',200,minimum=1,maximum=1000),
        field('diffusion_samples','Diffusion samples',5,minimum=1,maximum=100,section='basic')])
for _kind in ('openfold3_legacy','protenix_cli','boltzgen_fold_batch'):
    NODE_CATALOG[_kind]['replacement_group']='complex_prediction'

# Expose the cofolding engines already available as standalone tiles in the
# canvas replacement menu. Existing hand-adapted nodes keep their saved kind
# IDs; the remaining optimized profiles use the shared prediction workflow.
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as _optimized_prediction

OPTIMIZED_PREDICTION_NODE_ENGINES: dict[str, str] = {}
_existing_profile_nodes = {
    'alphafold3',  # the AF3 node below carries its terms checkbox explicitly
    'boltz2',      # the existing node keeps its familiar Boltz-2 input mask
    'rf3_foundry', # the existing RF3 node keeps its Foundry options
}


def _prediction_profile_fields(engine: str, profile: dict[str, Any]) -> list[dict[str, Any]]:
    fields = [
        field('optimization_mode', 'Optimization mode', 'off', kind='select',
              options=list(profile['modes']), section='basic'),
        field('samples', 'Samples per candidate', 2 if profile['api_engine'] == 'boltzgen_fold_opt_msa' else 1,
              minimum=1 if profile['api_engine'] != 'boltzgen_fold_opt_msa' else 2,
              maximum=int(profile.get('max_samples', 5)), section='basic'),
        field('num_recycles', 'Recycles', int(profile.get('default_recycles', 10)), minimum=1, maximum=100),
        field('seed', 'Random seed (0 = random)', 0, minimum=0, maximum=2147483647),
    ]
    if profile.get('supports_msa') and any(item['key'] == 'use_msa' for item in profile.get('settings', ())):
        msa_default = next(bool(item.get('default', True)) for item in profile.get('settings', ()) if item['key'] == 'use_msa')
        fields.append(field('target_use_msa', 'Use cached target MSA; fetch if missing', msa_default,
                            kind='boolean', section='basic',
                            help='Uses the target chain MSA. Newly designed binder chains stay single-sequence.'))
    if profile.get('supports_templates'):
        fields.append(field('use_target_template', 'Use target structure as template', True,
                            kind='boolean', section='basic'))
    if profile.get('requires_target_template'):
        fields.extend([
            field('target_pdb', 'Target PDB path', '', kind='text', section='basic',
                  help='BoltzGen Fold needs a PDB target structure; a target sequence alone is not sufficient.'),
            field('target_chain_id', 'Target PDB chain ID', 'A', kind='text', section='basic'),
        ])

    for setting in profile.get('settings', ()):
        key = setting['key']
        if key in {'use_msa', 'use_msa_server'}:
            continue
        kind = setting.get('type')
        if kind == 'bool':
            kind = 'boolean'
        elif kind == 'select':
            kind = 'select'
        else:
            kind = 'number'
        choices = setting.get('choices', ())
        options = [choice[0] if isinstance(choice, (tuple, list)) else choice for choice in choices]
        fields.append(field(
            key, setting.get('label', key.replace('_', ' ').capitalize()), setting.get('default'),
            kind=kind, minimum=setting.get('min'), maximum=setting.get('max'),
            options=options, help=setting.get('help'),
        ))
    return fields


def _prediction_task_support(engine: str, profile: dict[str, Any]) -> set[str]:
    """Return molecule tasks accepted by the shared prediction input adapter."""
    if engine == 'esmfold':
        # The shared ESMFold adapter accepts multichain FASTA; its profile text
        # flags complex use as exploratory because the model was trained on monomers.
        return {'protein', 'protein_protein'}
    if profile.get('input_kind') == 'boltzgen':
        return {'protein_protein'}
    supported = {'protein', 'protein_protein'}
    molecule_inputs = set(profile.get('molecule_inputs', ()))
    if molecule_inputs.intersection({'ligand_ccd_codes', 'ligand_smiles'}):
        supported.add('protein_ligand')
    if 'dna_sequences' in molecule_inputs:
        supported.add('protein_dna')
    if 'rna_sequences' in molecule_inputs:
        supported.add('protein_rna')
    if profile.get('supports_cyclic'):
        supported.add('protein_macrocyclic')
    return supported


for _engine, _profile in _optimized_prediction.ENGINE_PROFILES.items():
    if _engine in _existing_profile_nodes:
        continue
    _kind = f'prediction_{_engine}'
    _description = f"Refold connected target–binder candidates with {_profile['label']}."
    if _profile.get('input_summary'):
        _description += f" {_profile['input_summary']}"
    if _profile.get('fixed_settings'):
        _description += " Fixed settings: " + "; ".join(_profile['fixed_settings']) + "."
    NODE_CATALOG[_kind] = _spec(
        _kind, _profile['label'], 'Structure prediction', _description,
        accepts=['sequence', 'prediction', 'molecule_task'], output='prediction',
        workflow='refolding.optimized_prediction', fields=_prediction_profile_fields(_engine, _profile),
        color='#477e9b', formats=('pdb', 'cif'),
    )
    NODE_CATALOG[_kind]['replacement_group'] = 'complex_prediction'
    NODE_CATALOG[_kind]['molecule_tasks'] = sorted(_prediction_task_support(_engine, _profile))
    _ligand_inputs = set(_profile.get('molecule_inputs', ()))
    NODE_CATALOG[_kind]['ligand_input_formats'] = sorted(
        ({'ccd'} if 'ligand_ccd_codes' in _ligand_inputs else set())
        | ({'smiles'} if 'ligand_smiles' in _ligand_inputs else set())
    )
    NODE_CATALOG[_kind]['requires_target_structure'] = bool(_profile.get('requires_target_template'))
    OPTIMIZED_PREDICTION_NODE_ENGINES[_kind] = _engine

for _kind in ('rf3', 'boltz2', 'af3', 'af2', 'protenix'):
    NODE_CATALOG[_kind]['replacement_group'] = 'complex_prediction'
NODE_CATALOG['esmfold']['replacement_group'] = 'monomer_prediction'
for _kind, _engine in {'rf3': 'rf3_foundry', 'boltz2': 'boltz2', 'af3': 'alphafold3'}.items():
    _profile = _optimized_prediction.ENGINE_PROFILES[_engine]
    NODE_CATALOG[_kind]['molecule_tasks'] = sorted(_prediction_task_support(_engine, _profile))
    _ligand_inputs = set(_profile.get('molecule_inputs', ()))
    NODE_CATALOG[_kind]['ligand_input_formats'] = sorted(
        ({'ccd'} if 'ligand_ccd_codes' in _ligand_inputs else set())
        | ({'smiles'} if 'ligand_smiles' in _ligand_inputs else set())
    )
NODE_CATALOG['af2']['molecule_tasks'] = ['protein_protein']
NODE_CATALOG['af2']['requires_target_structure'] = True

from mn_molecular_design_studio import inverse_folding_input
for _engine, _metadata in inverse_folding_input.ENGINES.items():
    NODE_CATALOG[_engine] = _spec(_engine, _metadata["title"], "Sequence design", _metadata["description"],
        accepts=["target", "backbone", "sequence", "prediction"], output="sequence", workflow=inverse_folding_input.WORKFLOWS[_engine],
        fields=inverse_folding_input.canvas_fields(_engine), color="#a88845")

TEMPLATE_REFERENCES={
    'pxdesign':('https://github.com/bytedance/PXDesign#-quick-start','Binder Design Context feeds supported target-site annotations into PXDesign, followed by ProteinMPNN and independent AF2 Initial Guess and Protenix evaluation branches.'),
    'rfdiffusion':('https://github.com/RosettaCommons/RFdiffusion#sequence-design-for-binders','Binder Design Context supplies supported target-site annotations, then ProteinMPNN designs sequences followed by target-templated AF2 Initial Guess. The authors also permit ProteinMPNN without FastRelax; this template uses that option.'),
    'genie3':('https://github.com/aqlaboratory/genie3#application-3-binder-design','Binder Design Context supplies supported target-site annotations. One ProteinMPNN sequence per backbone is followed by AF2 target-template evaluation with five models and twenty recycles. Studio uses its AF2 Initial Guess adapter; the authors use ColabFold template mode.'),
    'protpardelle':('https://github.com/ProteinDesignLab/protpardelle-1c#sampling','Structure generation, eight ProteinMPNN sequences per backbone, then ESMFold self-consistency evaluation. The starter task is unconditional generation; switch to a conditioned task in the generator form.'),
    'esmfold2_binder':('https://github.com/Biohub/esm/tree/main/cookbook/tutorials','ESMFold2 Binder Design performs gradient-guided sequence optimization with its built-in ESMFold2 and ESMC scoring. Its candidate output can be sent through ProteinMPNN and then any compatible independent structure predictor; edit or replace these downstream starter steps on the canvas.'),
    'all_generators_boltz2':('https://github.com/RosettaCommons/RFdiffusion','Miniprotein binder workflow: branch from one Miniprotein Design Context into each available generator, use compact 60–100 residue designs where that engine exposes a binder-length range, then redesign with ProteinMPNN, refold with Boltz-2, and calculate ipSAE and Rosetta interface metrics. Each generator receives only annotations it supports.'),
    'vhh_design':('https://github.com/Biohub/esm/tree/main/cookbook/tutorials','VHH-focused workflow using only engines with an explicit VHH or nanobody mode: BindCraft2 VHH scaffold with CDR edit regions, BoltzGen nanobody task, and ESMFold2 Binder Design in VHH mode with a selectable framework. Their outputs converge at ProteinMPNN, then Boltz-2 and interface analysis. Hotspot and coldspot annotations are passed only to compatible generator paths.'),
    'cofold_protein':('https://github.com/google-deepmind/alphafold3','Compare the registered cofolding models that accept protein sequence inputs. ESMFold2 and other sequence-only engines are included for exploratory multimers; engines that require a target structure become runnable when one is provided.'),
    'cofold_protein_protein':('https://github.com/google-deepmind/alphafold3','Compare the registered protein–protein cofolding models. Enter target and partner sequences; add a target PDB path for AlphaFold2 Initial Guess and BoltzGen Fold, which require a structure template.'),
    'cofold_protein_ligand':('https://github.com/google-deepmind/alphafold3','Compare protein–ligand models that accept CCD codes or SMILES. Chai-1 accepts SMILES ligands; engines without ligand support are omitted.'),
    'cofold_protein_dna':('https://github.com/google-deepmind/alphafold3','Compare registered models whose molecule input accepts DNA. Protein-only cofolding engines are omitted.'),
    'cofold_protein_rna':('https://github.com/google-deepmind/alphafold3','Compare registered models whose molecule input accepts RNA. Protein-only cofolding engines are omitted.'),
    'cofold_protein_macrocyclic':('https://github.com/MatthiasZep/boltz','Compare models with explicit cyclic-protein support for a protein–macrocyclic-peptide complex. The current shared engine profile advertises native cyclic chains for Boltz-2.'),
}

COFOLDING_TASK_TEMPLATES = {
    'cofold_protein': ('protein', 'Protein monomer', 'Single-chain protein folding'),
    'cofold_protein_protein': ('protein_protein', 'Protein / protein', 'Protein–protein complex prediction'),
    'cofold_protein_ligand': ('protein_ligand', 'Protein / ligand', 'Protein–small-molecule complex prediction'),
    'cofold_protein_dna': ('protein_dna', 'Protein / DNA', 'Protein–DNA complex prediction'),
    'cofold_protein_rna': ('protein_rna', 'Protein / RNA', 'Protein–RNA complex prediction'),
    'cofold_protein_macrocyclic': ('protein_macrocyclic', 'Protein / macrocyclic peptide', 'Cyclic-peptide complex prediction'),
}


def _cofold_model_kinds(task_type: str) -> list[str]:
    return [
        kind for kind, spec in NODE_CATALOG.items()
        if spec.get('category') == 'Structure prediction'
        and 'molecule_task' in (spec.get('inputs') or [{}])[0].get('accepts', [])
        and task_type in spec.get('molecule_tasks', [])
    ]


def new_graph(template="foundry") -> dict[str, Any]:
    if template in COFOLDING_TASK_TEMPLATES:
        task_type, task_label, task_description = COFOLDING_TASK_TEMPLATES[template]
        model_kinds = _cofold_model_kinds(task_type)
        nodes = []
        input_spec = NODE_CATALOG['cofold_input']
        input_params = {f['name']: copy.deepcopy(f['default']) for f in input_spec['fields']}
        input_params['task_type'] = task_type
        if task_type == 'protein_macrocyclic':
            input_params['cyclic_partner'] = True
        nodes.append({
            'id': 'cofold-input', 'type': 'studio', 'position': {'x': 50, 'y': 0},
            'data': {'kind': 'cofold_input', 'label': f'{task_label} Inputs',
                     'description': task_description + '. Set sequences and optional structure or molecule inputs.',
                     'params': input_params},
        })
        edges = []
        branch_gap = 225
        center_y = 70 + max(0, len(model_kinds) - 1) * branch_gap / 2
        nodes[0]['position']['y'] = center_y
        for index, kind in enumerate(model_kinds):
            spec = NODE_CATALOG[kind]
            params = {f['name']: copy.deepcopy(f['default']) for f in spec['fields']}
            y = 70 + index * branch_gap
            model_id = f'{kind}-{index}'
            result_id = f'results-{index}'
            nodes.extend([
                {'id': model_id, 'type': 'studio', 'position': {'x': 430, 'y': y},
                 'data': {'kind': kind, 'label': spec['label'], 'params': params}},
                {'id': result_id, 'type': 'studio', 'position': {'x': 850, 'y': y},
                 'data': {'kind': 'results', 'label': f'Results · {spec["label"]}',
                          'params': {f['name']: copy.deepcopy(f['default']) for f in NODE_CATALOG['results']['fields']}}},
            ])
            edges.extend([
                {'id': f'edge-cofold-input-{index}', 'source': 'cofold-input', 'sourceHandle': 'out', 'target': model_id, 'targetHandle': 'in'},
                {'id': f'edge-cofold-result-{index}', 'source': model_id, 'sourceHandle': 'out', 'target': result_id, 'targetHandle': 'in'},
            ])
        graph = {
            'schema_version': SCHEMA, 'id': uuid.uuid4().hex,
            'name': f'{task_label} · compare compatible cofolding models',
            'nodes': nodes, 'edges': edges,
            'viewport': {'x': 0, 'y': 0, 'zoom': 0.5},
            'template_id': template,
        }
        return clean_graph(graph)
    if template == "all_generators_boltz2":
        generator_kinds = ["rfd3", "rfdiffusion", "genie3", "pxdesign", "protpardelle", "bindcraft2_generator", "boltzgen_generator", "complexa_generator", "esm_binder"]
        kinds = ["target", "binder_context", *generator_kinds, "ligandmpnn", "boltz2", "ipsae", "rosetta_metrics", "results"]
        nodes = []
        for index, kind in enumerate(kinds):
            spec = NODE_CATALOG[kind]
            params = {f["name"]: copy.deepcopy(f["default"]) for f in spec["fields"]}
            if kind in generator_kinds:
                if "num_designs" in params:
                    params["num_designs"] = 1
                if "binder_length" in params and kind != "esm_binder":
                    params["binder_length"] = "60-100"
                if kind in {"rfdiffusion", "genie3", "protpardelle"}:
                    params["task"] = "binder"
                if kind == "complexa_generator":
                    params["pipeline_task"] = "protein_binder"
                if kind == "esm_binder":
                    params["binder_type"] = "miniprotein"
                if kind == "bindcraft2_generator":
                    params["modality"] = "binder"
                    params["scaffold_edits"] = ""
            if kind == "ligandmpnn":
                params["sequences"] = 1
            if kind == "boltz2":
                params["write_full_pae"] = True
            nodes.append({
                "id": f"{kind}-{index}", "type": "studio",
                "position": {"x": 0, "y": 0},
                "data": {"kind": kind, "label": spec["label"], "params": params},
            })
        # Keep shared input and output steps centered on the parallel generator
        # branches. Placing them below the branches made the wires converge in
        # a long overlapping bundle and left most of the canvas unbalanced.
        branch_center_y = 50 + (len(generator_kinds) - 1) * 235 / 2
        nodes[0]["position"] = {"x": 40, "y": branch_center_y}
        nodes[1]["position"] = {"x": 400, "y": branch_center_y}
        nodes[1]["data"].update(
            label="Miniprotein Design Context",
            description="Define the target epitope for de novo miniprotein binders. Set hotspot residues the binder should contact and optional coldspots it should avoid. These annotations are inherited by each connected generator only when that engine supports them.",
        )
        for i, kind in enumerate(generator_kinds):
            nodes[i + 2]["position"] = {"x": 780, "y": 50 + i * 235}
        downstream = nodes[len(generator_kinds) + 2:]
        for i, node in enumerate(downstream):
            node["position"] = {"x": 1210 + i * 330, "y": branch_center_y}
        target, context = nodes[:2]
        edges = [{"id": "edge-target-context", "source": target["id"], "sourceHandle": "out", "target": context["id"], "targetHandle": "in"}]
        for i, node in enumerate(nodes[2:len(generator_kinds) + 2]):
            # Keep all generators in the same parallel stage. The context
            # node carries the shared structure; node_request removes any
            # annotations unsupported by an individual generator.
            edges.append({"id": f"edge-generator-input-{i}", "source": context["id"], "sourceHandle": "out", "target": node["id"], "targetHandle": "in"})
            edges.append({"id": f"edge-mpnn-{i}", "source": node["id"], "sourceHandle": "out", "target": downstream[0]["id"], "targetHandle": "in"})
        for i, (source, dest) in enumerate(zip(downstream, downstream[1:])):
            edges.append({"id": f"edge-analysis-{i}", "source": source["id"], "sourceHandle": "out", "target": dest["id"], "targetHandle": "in"})
        result = {
            "schema_version": SCHEMA, "id": uuid.uuid4().hex,
            "name": "Miniprotein binder design → ProteinMPNN → Boltz-2 + interface analysis",
            "nodes": nodes, "edges": edges,
            "viewport": {"x": 0, "y": 0, "zoom": 0.42},
            "template_id": template, "template_version": 3,
        }
        return clean_graph(result)
    if template == "vhh_design":
        generator_kinds = ["bindcraft2_generator", "boltzgen_generator", "esm_binder"]
        kinds = ["target", "binder_context", *generator_kinds, "ligandmpnn", "boltz2", "ipsae", "rosetta_metrics", "results"]
        nodes = []
        for index, kind in enumerate(kinds):
            spec = NODE_CATALOG[kind]
            params = {f["name"]: copy.deepcopy(f["default"]) for f in spec["fields"]}
            if kind == "bindcraft2_generator":
                from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import BC2_SCAFFOLD_EDIT_DEFAULTS
                params.update(num_designs=1, modality="VHH", binder_length="", scaffold_edits=BC2_SCAFFOLD_EDIT_DEFAULTS["VHH"])
            elif kind == "boltzgen_generator":
                params.update(num_designs=1, task="nanobody", binder_length="100-130")
            elif kind == "esm_binder":
                params.update(num_designs=1, binder_type="vhh", binder_framework="auto")
            elif kind == "ligandmpnn":
                params["sequences"] = 1
            elif kind == "boltz2":
                params["write_full_pae"] = True
            nodes.append({
                "id": f"{kind}-{index}", "type": "studio", "position": {"x": 0, "y": 0},
                "data": {"kind": kind, "label": spec["label"], "params": params},
            })
        center_y = 250
        nodes[0]["position"] = {"x": 40, "y": center_y}
        nodes[1]["position"] = {"x": 360, "y": center_y}
        nodes[1]["data"].update(
            label="VHH Design Context",
            description="Set the target epitope for VHH design. Hotspots and coldspots flow to compatible generator nodes; ESMFold2 VHH mode uses its selected target chain and framework settings instead.",
        )
        nodes[2]["data"]["label"] = "BindCraft2 VHH Generator"
        nodes[3]["data"]["label"] = "BoltzGen Nanobody Generator"
        nodes[4]["data"]["label"] = "ESMFold2 VHH Generator"
        nodes[2]["position"] = {"x": 700, "y": 20}
        nodes[3]["position"] = {"x": 700, "y": 260}
        nodes[4]["position"] = {"x": 700, "y": 500}
        downstream = nodes[len(generator_kinds) + 2:]
        for index, node in enumerate(downstream):
            node["position"] = {"x": 1040 + index * 320, "y": center_y}
        target, context = nodes[:2]
        edges = [{"id": "edge-target-context", "source": target["id"], "sourceHandle": "out", "target": context["id"], "targetHandle": "in"}]
        for index, node in enumerate(nodes[2:len(generator_kinds) + 2]):
            edges.append({"id": f"edge-vhh-input-{index}", "source": context["id"], "sourceHandle": "out", "target": node["id"], "targetHandle": "in"})
            edges.append({"id": f"edge-vhh-mpnn-{index}", "source": node["id"], "sourceHandle": "out", "target": downstream[0]["id"], "targetHandle": "in"})
        for index, (source, dest) in enumerate(zip(downstream, downstream[1:])):
            edges.append({"id": f"edge-vhh-analysis-{index}", "source": source["id"], "sourceHandle": "out", "target": dest["id"], "targetHandle": "in"})
        return clean_graph({
            "schema_version": SCHEMA, "id": uuid.uuid4().hex,
            "name": "VHH design · supported generators → ProteinMPNN → Boltz-2 + analysis",
            "nodes": nodes, "edges": edges,
            "viewport": {"x": 0, "y": 0, "zoom": 0.62},
            "template_id": template,
        })
    kinds = ["target", "binder_context", "rfd3", "foundry_mpnn", "rf3", "results"] if template == "foundry" else ["target", "results"]
    if template == "bindcraft2":
        kinds = ["target", "binder_context", "bindcraft2", "results"]
    elif template == "boltzgen":
        kinds = ["target", "binder_context", "boltzgen", "results"]
    elif template in {"complexa","complexa_generated"}:
        kinds = ["complexa_generator" if template=="complexa" else "generation_source","ligandmpnn","rf3","results"]
    elif template == "esmfold2_binder":
        kinds = ["target", "esm_binder", "ligandmpnn", "rf3", "results"]
    engine=template.removesuffix('_generated')
    if engine in generator_input.ENGINES:
        kinds=['backbone_source' if template.endswith('_generated') else engine,'ligandmpnn','esmfold' if engine=='protpardelle' else 'af2','results']
        if not template.endswith('_generated') and _uses_binder_context(engine):
            kinds = ['target', 'binder_context', *kinds]
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
    elif 'binder_context' in kinds and engine != 'pxdesign':
        for index, node in enumerate(nodes):
            node['position'] = {'x': 50 + index * 330, 'y': 180}
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
            if 'binder_context' in kinds:
                by_kind = {node['data']['kind']: node for node in nodes}
                target = by_kind['target']
                context = by_kind['binder_context']
                generator = by_kind['pxdesign']
                sequence = by_kind['ligandmpnn']
                predictions = [node for node in nodes if node['data']['kind'] in {'af2', 'protenix'}]
                results = by_kind['results']
                target['position']={'x':50,'y':180};context['position']={'x':340,'y':180}
                generator['position']={'x':630,'y':180};sequence['position']={'x':920,'y':180}
                for prediction, y in zip(predictions, (50, 300)):
                    prediction['position']={'x':1210,'y':y}
                results['position']={'x':1510,'y':180}
                pairs=[(target,context),(context,generator),(generator,sequence),
                       *((sequence,prediction) for prediction in predictions),
                       *((prediction,results) for prediction in predictions)]
            else:
                a,b,c,d,r=nodes
                a['position']={'x':50,'y':170};b['position']={'x':380,'y':170};c['position']={'x':710,'y':50};d['position']={'x':710,'y':300};r['position']={'x':1040,'y':170}
                pairs=[(a,b),(b,c),(b,d),(c,r),(d,r)]
            edges=[{'id':f'edge-{i}','source':source['id'],'sourceHandle':'out','target':dest['id'],'targetHandle':'in'} for i,(source,dest) in enumerate(pairs)]
    result={"schema_version": SCHEMA, "id": uuid.uuid4().hex, "name": {"foundry": "Foundry binder design", "bindcraft2": "BindCraft2 design", "boltzgen": "BoltzGen design", "complexa":"Proteina-Complexa custom pipeline", "complexa_generated":"Continue generated designs", "esmfold2_binder":"ESMFold2 Binder Design → ProteinMPNN → RosettaFold3"}.get(template, generator_input.ENGINES[engine]['title']+' · author-inspired pipeline' if engine in generator_input.ENGINES else "Untitled pipeline"), "nodes": nodes, "edges": edges, "viewport": {"x": 0, "y": 0, "zoom": 0.75}}
    if template in TEMPLATE_REFERENCES:result['template_id']=template
    return result


def migrate_template_graph(graph: dict[str, Any]) -> dict[str, Any]:
    """Upgrade old All generators drafts to the miniprotein design template.

    Keep node IDs, settings, and the pipeline name, but restore the template's
    non-overlapping columns and route every generator from Binder Design
    Context. This is a one-time layout migration, recorded on the graph.
    """
    if graph.get("template_id") != "all_generators_boltz2" or graph.get("template_version") == 3:
        return graph

    migrated = copy.deepcopy(graph)
    canonical = new_graph("all_generators_boltz2")
    canonical_positions = {node["data"]["kind"]: node["position"] for node in canonical["nodes"]}
    nodes_by_kind = {}
    for node in migrated.get("nodes", []):
        nodes_by_kind.setdefault(node.get("data", {}).get("kind"), []).append(node)
        position = canonical_positions.get(node.get("data", {}).get("kind"))
        if position is not None:
            node["position"] = copy.deepcopy(position)

    target_nodes = nodes_by_kind.get("target", [])
    context_nodes = nodes_by_kind.get("binder_context", [])
    if not target_nodes or not context_nodes:
        return graph
    target_id, context_id = target_nodes[0]["id"], context_nodes[0]["id"]
    context = context_nodes[0]
    if context.get("data", {}).get("label") == "Binder Design Context":
        context["data"]["label"] = "Miniprotein Design Context"
    context["data"]["description"] = "Define the target epitope for de novo miniprotein binders. Set hotspot residues the binder should contact and optional coldspots it should avoid. These annotations are inherited by each connected generator only when that engine supports them."
    if migrated.get("name") in {
        "All generators → ProteinMPNN → Boltz-2 + analysis",
        "All generators → ProteinMPNN → Boltz-2 + interface analysis",
    }:
        migrated["name"] = "Miniprotein binder design → ProteinMPNN → Boltz-2 + interface analysis"
    generator_kinds = {
        "rfd3", "rfdiffusion", "genie3", "pxdesign", "protpardelle",
        "bindcraft2_generator", "boltzgen_generator", "complexa_generator", "esm_binder",
    }
    generator_ids = {
        node["id"] for kind in generator_kinds for node in nodes_by_kind.get(kind, [])
    }
    # Reuse existing edge IDs so the canvas history remains understandable.
    # If an old graph already has the context edge, remove its direct bypass.
    for generator_id in generator_ids:
        incoming = [edge for edge in migrated.get("edges", []) if edge.get("target") == generator_id]
        has_context_edge = any(edge.get("source") == context_id for edge in incoming)
        if has_context_edge:
            migrated["edges"] = [
                edge for edge in migrated["edges"]
                if not (edge.get("target") == generator_id and edge.get("source") == target_id)
            ]
        else:
            for edge in incoming:
                if edge.get("source") == target_id:
                    edge["source"] = context_id

    migrated["viewport"] = copy.deepcopy(canonical["viewport"])
    # Narrow untouched length defaults to the miniprotein template's compact
    # range while preserving any values the user already changed.
    for node in migrated.get("nodes", []):
        kind = node.get("data", {}).get("kind")
        params = node.get("data", {}).get("params", {})
        if kind not in generator_kinds or kind == "esm_binder" or "binder_length" not in params:
            continue
        default_length = next((item["default"] for item in NODE_CATALOG[kind]["fields"] if item["name"] == "binder_length"), None)
        if params.get("binder_length") == default_length:
            params["binder_length"] = "60-100"
    migrated["template_version"] = 3
    return clean_graph(migrated)


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
        description = data.get("description")
        if description is not None and not isinstance(description, str):
            raise ValueError("Node descriptions must be text.")
        node_data = {"kind": kind, "label": str(data.get("label") or NODE_CATALOG[kind]["label"])[:100], "params": params}
        if description and description.strip():
            node_data["description"] = description.strip()[:1000]
        cleaned.append({"id": node["id"], "type": "studio", "position": position, "data": node_data})
    viewport = graph.get("viewport") or {"x": 0, "y": 0, "zoom": 1}
    if not isinstance(viewport, dict) or any(isinstance(viewport.get(k), bool) or not isinstance(viewport.get(k), (float, int)) or not math.isfinite(viewport[k]) for k in ("x", "y", "zoom")) or not .05 <= viewport["zoom"] <= 5:
        raise ValueError("Invalid canvas viewport.")
    return {"schema_version": SCHEMA, "id": graph["id"], "name": str(graph.get("name") or "Untitled pipeline")[:160], "nodes": cleaned,
            "edges": [{k: edge[k] for k in ("id", "source", "target", "sourceHandle", "targetHandle") if k in edge} for edge in edges if isinstance(edge, dict)],
            "viewport": viewport,**({'template_id':graph['template_id']} if graph.get('template_id') in TEMPLATE_REFERENCES else {}),
            **({'template_version':graph['template_version']} if graph.get('template_id') == 'all_generators_boltz2' and type(graph.get('template_version')) is int and 1 <= graph['template_version'] <= 3 else {})}


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
    cofold_task_context = {}
    cofold_source_context = {}
    for node_id in order:
        kind = nodes[node_id]["data"]["kind"]
        spec = NODE_CATALOG[kind]
        inherited_tasks = {cofold_task_context[source] for source in incoming[node_id] if source in cofold_task_context}
        inherited_sources = {cofold_source_context[source] for source in incoming[node_id] if source in cofold_source_context}
        cofold_task_context[node_id] = (
            nodes[node_id]["data"]["params"].get("task_type") if kind == "cofold_input"
            else next(iter(inherited_tasks)) if len(inherited_tasks) == 1
            else None
        )
        cofold_source_context[node_id] = (
            node_id if kind == 'cofold_input'
            else next(iter(inherited_sources)) if len(inherited_sources) == 1
            else None
        )
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
        if require_inputs and kind == 'cofold_input':
            params = node['data']['params']
            task = params.get('task_type')
            if not params.get('target_protein_sequences', '').strip():
                raise ValueError('Enter at least one target protein sequence in Cofolding Inputs.')
            target_count = len([part for part in re.split(r'[:\n\r]+', params.get('target_protein_sequences', '')) if part.strip()])
            partner_count = len([part for part in re.split(r'[:\n\r]+', params.get('partner_protein_sequences', '')) if part.strip()])
            if task == 'protein' and target_count != 1:
                raise ValueError('The protein monomer task accepts exactly one protein sequence.')
            if task in {'protein_protein', 'protein_macrocyclic'} and not params.get('partner_protein_sequences', '').strip():
                raise ValueError('Enter at least one partner protein sequence for this cofolding task.')
            if task == 'protein_macrocyclic' and partner_count != 1:
                raise ValueError('The protein–macrocyclic peptide task accepts exactly one cyclic partner sequence.')
            if task == 'protein_ligand' and not (params.get('ligand_ccd_codes', '').strip() or params.get('ligand_smiles', '').strip()):
                raise ValueError('Enter a CCD code or SMILES string for the protein–ligand task.')
            ligand_format = params.get('ligand_input_type', 'smiles')
            if task == 'protein_ligand' and ligand_format == 'ccd' and not params.get('ligand_ccd_codes', '').strip():
                raise ValueError('Enter a CCD code for the selected ligand notation.')
            if task == 'protein_ligand' and ligand_format == 'smiles' and not params.get('ligand_smiles', '').strip():
                raise ValueError('Enter a SMILES string for the selected ligand notation.')
            if task == 'protein_dna' and not params.get('dna_sequences', '').strip():
                raise ValueError('Enter a DNA sequence for the protein–DNA task.')
            if task == 'protein_rna' and not params.get('rna_sequences', '').strip():
                raise ValueError('Enter an RNA sequence for the protein–RNA task.')
            if task == 'protein_macrocyclic' and not params.get('cyclic_partner'):
                raise ValueError('Mark the protein partner as cyclic for the macrocyclic task.')
            optional_target_path = str(params.get('target_pdb') or '').strip()
            if optional_target_path and not Path(optional_target_path).expanduser().is_file():
                raise ValueError('The optional target PDB path in Cofolding Inputs does not exist.')
        if require_inputs and NODE_CATALOG[kind].get('requires_target_structure'):
            if cofold_task_context.get(node_id):
                source_id = cofold_source_context.get(node_id)
                source_input = nodes[source_id]['data']['params'] if source_id in nodes else None
                path = str((source_input or {}).get('target_pdb') or '').strip()
                if not path or not Path(path).expanduser().is_file():
                    raise ValueError(f"{node['data']['label']} requires an existing target PDB path in Cofolding Inputs.")
                if kind == 'af2' and Path(path).suffix.lower() != '.pdb':
                    raise ValueError('AlphaFold2 Initial Guess requires a PDB file (mmCIF is not supported by this input mode).')
        if require_inputs and cofold_task_context.get(node_id) and NODE_CATALOG[kind].get('molecule_tasks'):
            requested_task = cofold_task_context[node_id]
            if requested_task not in NODE_CATALOG[kind]['molecule_tasks']:
                raise ValueError(f"{node['data']['label']} does not support the selected {requested_task.replace('_', ' / ')} prediction task.")
            if requested_task == 'protein_ligand' and NODE_CATALOG[kind].get('ligand_input_formats'):
                source_id = cofold_source_context.get(node_id)
                ligand_format = nodes[source_id]['data']['params'].get('ligand_input_type', 'smiles') if source_id in nodes else 'smiles'
                if ligand_format not in NODE_CATALOG[kind]['ligand_input_formats']:
                    raise ValueError(f"{node['data']['label']} does not accept {ligand_format.upper()} ligand inputs.")
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
    "design.rfpeptides": "Design engine · downstream sequence design available",
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
