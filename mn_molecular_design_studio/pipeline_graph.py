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
from mn_molecular_design_studio.docking_preparation_protocols import (
    DEFAULT_SCRUB_PROTOCOL, SCRUB_PROTOCOLS,
)

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
    "peptide_trap": _spec("peptide_trap", "Peptide Trap", "Design engine", "Retarget published β-strand clefts by sequence matching, threading and Rosetta interface redesign.", output="sequence", workflow="design.peptide_trap", fields=[
        field("target_sequence", "Target sequence / FASTA", "", kind="text", section="basic"),
        field("scaffold_ids", "Scaffold templates", ["C37", "C104.1"], kind="json", section="basic"),
        field("max_registers", "Maximum matching placements", 20, minimum=1, maximum=1000),
        field("designs_per_register", "Designs per placement", 1, minimum=1, maximum=100),
        field("design_repeats", "Rosetta design repeats", 2, minimum=1, maximum=5),
        field("relax_repeats", "Rosetta relaxation repeats", 2, minimum=1, maximum=5),
        field("max_backbone_rmsd", "Maximum template backbone RMSD (Å)", 2.0, minimum=0.1, maximum=10.),
        field("cpu_cores", "CPU cores", 4, minimum=1, maximum=256), SEED], formats=("pdb",)),
    "kejia_idr": _spec("kejia_idr", "Kejia IDR threading", "Design engine", "Thread a target sequence onto binder–peptide templates and refine peptide flexibility; connect ProteinMPNN and AlphaFold2 downstream.", workflow="design.kejia_idr", fields=[
        field("target_sequence", "Target fragment / single-record FASTA", "", kind="text", section="basic"),
        field("template_library", "Template source", "upload", kind="select", options=["upload", "walle", "walle2", "mini", "strand"], section="basic"),
        field("template_structure", "Uploaded template path", "", kind="text", section="basic"),
        field("binder_chain", "Binder chain", "A", kind="text"), field("peptide_chain", "Peptide chain", "B", kind="text"),
        field("max_templates", "Maximum library templates", 10, minimum=1, maximum=1000),
        field("flex_docking_repeats", "Flexible docking repeats per register", 3, minimum=1, maximum=100),
        field("min_seq_identity", "Minimum sequence identity", 0.0, minimum=0., maximum=1.),
        field("min_blosum", "Minimum BLOSUM score", -10.0, minimum=-100., maximum=100.),
        field("max_farep", "Maximum repulsive-energy increase", 10.0, minimum=0., maximum=10000.),
        field("keep_bidentates", "Retain template bidentate contacts", False, kind="boolean"),
        field("cpu_cores", "CPU cores", 4, minimum=1, maximum=256), SEED], formats=("pdb",)),
    "cofold_input": _spec("cofold_input", "Cofolding Inputs", "Cofolding input", "Enter a protein target and any partner molecules for the selected prediction task.", output="molecule_task", color="#758f84", fields=[
        field("task_type", "Prediction task", "protein_protein", kind="select", options=["protein", "protein_protein", "protein_ligand", "protein_dna", "protein_rna", "protein_macrocyclic"], section="basic"),
        field("target_protein_sequences", "Target protein sequence(s)", "", kind="text", section="basic", help="Enter one sequence per target chain. Separate chains with a colon or new line."),
        field("partner_protein_sequences", "Protein partner sequence(s)", "", kind="text", section="basic", help="Used for protein–protein and protein–macrocycle tasks. Separate partner chains with a colon or new line."),
        field("dna_sequences", "DNA sequence(s)", "", kind="text", section="basic", help="Enter DNA chains when using the protein–DNA task."),
        field("rna_sequences", "RNA sequence(s)", "", kind="text", help="Enter RNA chains when using the protein–RNA task."),
        field("ligand_input_type", "Ligand notation", "smiles", kind="select", options=["smiles", "ccd"], help="SMILES works across the broadest set of ligand-capable models; CCD-only inputs exclude models that accept SMILES only."),
        field("ligand_ccd_codes", "Ligand CCD code(s)", "", kind="text", help="One or more CCD component IDs, separated by colons or new lines."),
        field("ligand_smiles", "Ligand SMILES", "", kind="text", help="Enter one SMILES string per line."),
        field("compound_set_file", "Compound collection", "", kind="text", help="A Database compound snapshot. Each compound is predicted independently in one batch against the shared target."),
        field("sequence_file", "Sequence collection", "", kind="text", help="FASTA records or one sequence per line. Every record is one independent prediction; colon-separated chains within a record form one complex."),
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
    "molecule_source": _spec("molecule_source", "Database selection", "Input", "An immutable selection of Database molecules, prepared through Use in pipeline.", output="sequence", color="#758f84", formats=("pdb", "cif"), fields=[field("source_run_dir", "Selection snapshot", "", kind="text")]),
    "design_source": _spec("design_source", "Successful project designs", "Input", "Selected filter-passing designs from a successful job in this project.", output="sequence", color="#758f84", fields=[field("source_run_dir", "Source run", "", kind="text")]),
    "generation_source": _spec("generation_source", "Generated project designs", "Input", "Select generated structures and sequences from a completed generator job in this project.", output="sequence", color="#758f84", formats=("pdb", "cif"), fields=[field("source_run_dir", "Generator run", "", kind="text")]),
    "project_design_source": _spec("project_design_source", "Existing project designs", "Input", "Select structure-bearing candidates from any completed design job in this project, including generator and refold outputs.", output="sequence", color="#758f84", formats=("pdb", "cif"), fields=[field("source_run_dir", "Design job", "", kind="text")]),
    "rfd4": _spec("rfd4", "RFD4-Proteina", "Design engine", "Co-design sequence and structure; attach any compatible prediction and analysis nodes.", accepts=["target"], output="sequence", workflow="design.rfd4_proteina", fields=[
        field("task", "Design task", "protein_binder", kind="select", options=["protein_binder","nanobody","antibody","tcr","motif_scaffolding","enzyme","small_molecule_binder","dna_binder","rna_design","existing_structure","custom_conditions","score_structures","unconditional"]),
        field("target_structure", "Input structure", "", kind="text"),
        field("input_mode", "Input method", "guided", kind="select", options=["guided","annotated_cif","condition_spec"]),
        field("condition_spec", "Native condition specification", "", kind="text", help="Native YAML/JSON design specification used with the condition_spec input method."),
        field("input_files", "Input structure files (one path per line)", "", kind="text", help="For annotated_cif mode, list one annotated CIF path per design. Structure scoring can also take multiple inputs."),
        field("options_json", "Task options (JSON)", "{}", kind="text", help="Multiline JSON matching the RFD4 task options, including motif residues, active-site atoms and ligand interactions. The standalone RFD4 form exposes these controls directly."),
        COUNT, field("length_min", "Minimum length", 60, minimum=20, maximum=1024), field("length_max", "Maximum length", 120, minimum=20, maximum=1024),
        field("nsteps", "Sampling steps", 60, minimum=60, maximum=120), SEED], formats=("cif",)),
    "rfd4_fold": _spec("rfd4_fold", "RFD4-Proteina Structure Prediction", "Structure prediction", "Fold supplied protein sequences or connected generated structures using RFD4 itself.", accepts=["sequence","prediction","target","molecule_task"], output="prediction", workflow="refolding.rfd4_proteina", fields=[
        field("protein_sequences", "Protein sequence(s) / FASTA", "", kind="text"), field("target_structure", "Input structure", "", kind="text"),
        field("num_designs", "Predictions per input", 1, minimum=1, maximum=1000),field("nsteps", "Sampling steps",60,minimum=60,maximum=120),SEED], formats=("cif",)),
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
    "sabr": _spec("sabr", "SAbR · Antibody CDR Annotation", "Antibody analysis", "Apply SAbR numbering to each connected antibody design and annotate CDR1, CDR2, and CDR3. This step passes the same structures and candidate metadata downstream.", accepts=["backbone", "sequence", "prediction", "candidates"], output="sequence", color="#54866c", fields=[
        field("chain_type", "Antibody chain type", "auto", kind="select", options=["auto", "H", "K", "L"]),
        field("mode", "Alignment mode", "sabr", kind="select", options=["sabr", "softalign"]),
        field("noise_level", "Reference backbone noise (Å)", 0.0, minimum=0.0, maximum=2.0),
    ]),
    "abmpnn": _spec("abmpnn", "AbMPNN · CDR Sequence Redesign", "Sequence design", "Redesign selected antibody CDRs using the connected SAbR annotation; keep the target, framework, and unselected CDRs fixed.", accepts=["backbone", "sequence", "prediction", "candidates"], output="sequence", workflow="sequence.abmpnn", fields=[
        field("sequences", "Sequences per design", 1, minimum=1, maximum=10000, section="basic"),
        field("redesign_cdr1", "Redesign CDR1", True, kind="boolean", section="basic"),
        field("redesign_cdr2", "Redesign CDR2", True, kind="boolean", section="basic"),
        field("redesign_cdr3", "Redesign CDR3", True, kind="boolean", section="basic"),
        field("sampling_temp", "Sampling temperature", 0.1, minimum=0.000001, maximum=1.0),
        field("omit_aas", "Amino acids to omit", "C", kind="text"),
        field("seed", "Random seed (0 = random)", 0, minimum=0, maximum=4294967295),
    ], color="#a88845", formats=("pdb", "cif")),
    "rf3": _spec("rf3", "RosettaFold3", "Structure prediction", "Refold designed target–binder complexes with Foundry RF3.", accepts=["sequence", "prediction", "molecule_task"], output="prediction", workflow="refolding.rf3", fields=[MSA, TEMPLATE], color="#477e9b", formats=("cif",)),
    "boltz2": _spec("boltz2", "Boltz-2", "Structure prediction", "Refold target–binder complexes with configurable target MSA and optional target template conditioning.", accepts=["sequence", "prediction", "molecule_task"], output="prediction", workflow="refolding.optimized_prediction", fields=[
        field("samples", "Samples per candidate", 1, minimum=1, maximum=5),
        field("num_recycles", "Recycles", 3, minimum=1, maximum=100),
        field("sampling_steps", "Diffusion steps", 200, minimum=1, maximum=1000),
        field("seed_count", "Number of seeds", 1, minimum=1, maximum=16),
        SEED, field("target_use_msa", "Use cached target MSA; fetch if missing", True, kind="boolean"),
        TEMPLATE,
    ], color="#477e9b", formats=("cif",)),
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
    "af2": _spec("af2", "AlphaFold2 evaluation", "Structure prediction", "Choose conditioned pose self-consistency or a less-conditioned target-only template prediction.", accepts=["sequence", "prediction", "molecule_task"], output="prediction", workflow="refolding.af2_initial_guess", fields=[field("num_recycles", "Recycles", 3, minimum=1, maximum=30)], color="#477e9b"),
    "esmfold": _spec("esmfold", "ESMFold · binder alone", "Structure prediction", "Predict each designed binder without the target.", accepts=["sequence", "prediction"], output="prediction", workflow="refolding.esmfold_monomer", fields=[], color="#477e9b"),
    "filter": _spec("filter", "Filter candidates", "Analysis", "Apply a score threshold, require passing filters, or retain the best candidates.", accepts=["backbone", "sequence", "prediction", "candidates"], output="candidates", color="#9a7854", fields=[
        field("passing_only", "Require an explicit native filter pass", False, kind="boolean"), field("metric", "Threshold measurement", "", kind="text"), field("operator", "Comparison", ">=", kind="select", options=[">=", "<=", ">", "<"]), field("threshold", "Threshold (only used with a measurement)", 0.7), field("limit", "Shortlist per comparable protocol (0 = all)", 0, minimum=0, maximum=100000), field("higher_is_better", "Higher scores rank first", True, kind="boolean"),
        field("required_metrics", "Required finite measurements", [], kind="json"),
        field("rules", "Additional threshold rules", [], kind="json", help='List of {"field": "...", "operator": ">=", "value": ...}. All rules must pass.'),
        field("ranking_metric", "Ranking measurement (blank = threshold measurement)", "", kind="text"),
        field("manual_ids", "Manually selected candidate IDs (optional)", [], kind="json"),
        field("diversity_field", "Cluster measurement (optional)", "", kind="text"),
        field("max_per_cluster", "Maximum selected per cluster", 1, minimum=1, maximum=100000)]),
    "surface_analysis": _spec("surface_analysis", "Protein Surface Analysis", "Analysis", "Run PeSTo, ScanNet, Surf2Spot, and MaSIF-seed against one target and compare predicted PPI and hotspot residues.", accepts=["target"], output="analysis", workflow="detection.surface_analysis", color="#397b5b", formats=("pdb",), fields=[
        field("pesto_interface", "PeSTo interface class", "protein", kind="select", options=["protein", "dna_rna", "ion", "ligand", "lipid"]),
        field("scannet_mode", "ScanNet prediction mode", "interface", kind="select", options=["interface", "epitope", "idp"]),
        field("scannet_use_msa", "Use MSA with ScanNet", False, kind="boolean"),
        field("masif_chain_id", "MaSIF-seed chain", "", kind="text")]),
    "results": _spec("results", "Results & comparison", "Results", "Collect candidate structures, metrics, and surface-analysis runs from connected branches.", accepts=["backbone", "sequence", "prediction", "candidates", "analysis"], output=None, color="#5d7b6c"),
    "ipsae": _spec("ipsae", "ipSAE Interface Analysis", "Analysis", "Score predicted interfaces using the saved pairwise PAE matrix and complex structure.", accepts=["prediction"], output="prediction", workflow="analysis.ipsae", color="#9a7854", fields=[field("pae_cutoff", "PAE cutoff (Å)", 10., minimum=0., maximum=100.), field("dist_cutoff", "Interface distance cutoff (Å)", 10., minimum=0., maximum=100.)]),
    "rosetta_metrics": _spec("rosetta_metrics", "Rosetta Interface Analysis", "Analysis", "Calculate PyRosetta interface, packing, SASA, hydrogen-bond, and binder surface metrics.", accepts=["prediction"], output="prediction", workflow="analysis.rosetta", color="#9a7854", fields=[]),
}

# Active-learning docking is represented as ordinary editable canvas nodes.
# Its feedback edge is marked in the graph and executed by the dedicated
# controller; the generic DAG runner must not try to schedule these nodes.
_ACTIVE_STATE = "active_learning_state"
NODE_CATALOG.update({
    "active_library": _spec("active_library", "Compound library", "Active learning", "Use a registered external library's packed fingerprints for prediction and random exploration. The initial centroid sample has its own independent source. Compound rows remain in external files; only docked batches are materialized.", output=_ACTIVE_STATE, color="#3879a8", formats=(), fields=[
        field("library_id", "Registered library ID", "", kind="text", section="basic"),
        field("feature_manifest", "Fingerprint manifest", "", kind="text", section="basic"),
        # Retained only so older saved drafts can migrate this value to the seed input.
        field("centroid_manifest", "Legacy centroid manifest", "", kind="text", section="legacy"),
    ]),
    "active_target": _spec("active_target", "Receptor", "Active learning", "Select an imported receptor structure and define the docking box for this study.", output=_ACTIVE_STATE, color="#758f84", formats=("pdb",), fields=[
        field("target_pdb", "Saved receptor structure", "", kind="text", section="basic"),
        field("target_name", "Receptor name", "", kind="text", section="basic"),
        field("align_longest_axis_to_x", "Align reference ligand longest axis to X", False, kind="boolean", section="basic"),
        field("alignment_ligand_sdf", "Reference ligand SDF for axis alignment", "", kind="text", section="basic"),
        field("reference_ligand_sdf", "Reference ligand SDF for docking mode", "", kind="text", section="basic",
              help="Hybrid mode uses both the receptor and this bound/reference ligand. Ligand-based mode uses the ligand as its placement reference. The imported THRα target preselects its co-crystal ligand."),
        field("center_x", "Box center X (Å)", 0.0, minimum=-500, maximum=500, section="basic"),
        field("center_y", "Box center Y (Å)", 0.0, minimum=-500, maximum=500, section="basic"),
        field("center_z", "Box center Z (Å)", 0.0, minimum=-500, maximum=500, section="basic"),
        field("size_x", "Box size X (Å)", 25.0, minimum=1, maximum=60, section="basic"),
        field("size_y", "Box size Y (Å)", 25.0, minimum=1, maximum=60, section="basic"),
        field("size_z", "Box size Z (Å)", 25.0, minimum=1, maximum=60, section="basic"),
    ]),
    "active_seed": _spec("active_seed", "Random centroid seed", "Active learning", "Choose a registered centroid dataset and dock a reproducible random sample for the first training round. This source can differ from the full screening library; fingerprint settings must be compatible.", output=_ACTIVE_STATE, color="#3985a7", fields=[
        field("library_id", "Centroid source library ID", "", kind="text", section="basic"),
        field("centroid_manifest", "Centroid manifest", "", kind="text", section="basic"),
        field("sample_size", "Initial compounds", 10000, minimum=1, section="basic"),
        field("sampling_seed", "Sampling seed", 42, minimum=0, maximum=2147483647, section="basic"),
    ]),
    "active_docking": _spec("active_docking", "Docking", "Active learning", "Dock each selected round with the chosen engine, splitting the round across selected Studio computers. The initial centroid sample sets the first round size; the selection node sets later round sizes.", accepts=[_ACTIVE_STATE], output=_ACTIVE_STATE, color="#427f68", fields=[
        field("engine", "Docking engine", "UniDock Pro", kind="select", options=["UniDock Pro"], section="basic",
              help="The active-learning controller currently connects UniDock Pro. Other engines need their score and output adapters before they can run in this loop."),
        field("scrub_protocol", "Ligand-preparation protocol", DEFAULT_SCRUB_PROTOCOL, kind="select",
              options=list(SCRUB_PROTOCOLS), section="basic",
              help="Preparation path: SMILES → RDKit 3D starting structure → MolScrub pH-dependent acid/base preparation and ETKDG conformer/minimization → Meeko PDBQT → Uni-Dock Pro. Acid/base preparation stays enabled; pH defaults to 7.4 and tautomers are skipped. In the same 10,000-molecule preparation benchmark (docking excluded), 0.1.0 took 94.8 s (104.5 PDBQT/s), 0.3.0 defaults took 317.8 s (31.2/s), and 0.3.0 matched settings took 201.8 s (49.1/s). All used MMFF94s; the modes differ in version, conformer count, and force-field iterations."),
        field("scrub_ph", "Ligand preparation pH", 7.4, minimum=0.0, maximum=14.0, section="basic"),
        field("scrub_skip_tautomer", "Skip tautomer enumeration", True, kind="boolean", section="basic"),
        field("cleanup_preparation", "Remove batch preparation files after docking", True, kind="boolean", section="basic",
              help="Use temporary scratch for preparation and docking. Keep scores, preparation timing and failure records; never retry failed molecules individually."),
        field("cleanup_after_round", "Remove worker poses after collecting round winners", True, kind="boolean", section="basic",
              help="After each round, transfer and verify the cumulative best poses on the organizer, then remove bulk worker poses. Retention count is set in Accumulate results."),
        field("score_field", "Score field", "vina_score", kind="text", section="basic"),
        field("score_units", "Score units", "kcal/mol", kind="text", section="basic", help="Declare the units for the selected scoring function, or use unitless."),
        field("score_direction", "Score ranking", "Lower is better", kind="select", options=["Lower is better", "Higher is better"], section="basic"),
        field("batch_size", "Maximum compounds in later rounds (0 = no cap)", 0, minimum=0, section="basic",
              help="Optional cap for active-selection batches after the first round. By default, later rounds use the configured percentage of the initial centroid sample. Large rounds are split into UniDock Pro-sized launches on each selected computer."),
        field("worker_batch_size", "Compounds per computer launch", 2000, minimum=1, maximum=150000, section="basic",
              help="Each selected computer docks batches of this size in its scheduled UniDock Pro jobs. Preparation and docking run on that computer. Increase this based on available GPU memory and measured throughput; the Studio adapter supports up to 150,000 compounds per launch."),
        field("docking_mode", "UniDock Pro docking mode", "hybrid", kind="select",
              options=["classic", "ligand_based", "hybrid"], section="basic",
              help="Classic is structure-based (receptor only). Ligand-based uses the reference ligand. Hybrid uses both receptor and reference ligand. This is separate from Search effort."),
        field("search_mode", "UniDock Pro search effort", "balance", kind="select", options=["fast", "balance", "detail"], section="basic",
              help="Search effort is independent of docking mode. Balance is the selected effort for this study."),
        field("num_modes", "Docking poses retained per compound during scoring", 1, minimum=1, maximum=20, section="basic"),
        field("worker_nodes", "Compute computers", [], kind="json", section="basic", help="Select the organizer and available compute computers. Their UniDock Pro batches run in parallel through each computer's Studio scheduler."),
        field("ssh_user", "SSH user", "compute", kind="text", section="basic"),
    ]),
    "active_accumulate": _spec("active_accumulate", "Accumulate docking results", "Active learning", "Merge every completed docking round, keep the best configured score per source library and stable compound ID, and sort the cumulative table using the docking node's ranking direction. Source namespaces preserve distinct compounds when two libraries reuse an ID. Snapshot the top ranks for comparison with the previous cumulative ranking. Keep poses only for the best-ranked compounds.", accepts=[_ACTIVE_STATE], output=_ACTIVE_STATE, color="#668d75", fields=[
        field("retained_pose_count", "Best pose files retained externally", 1000, minimum=0, maximum=1000000, section="basic", help="Physical pose files kept in the campaign folder. This does not import compounds into Database; set that separately in Final results."),
    ]),
    "active_binder_rule": _spec("active_binder_rule", "Binder labeling rule", "Active learning", "Recalculate binder labels on the cumulative docked dataset before each training round. Change the binder percentage logarithmically from the starting percentage to the final percentage across the selected rounds; keep the final percentage afterwards. Label the best-scoring percentage using the ranking direction declared by the docking node.", accepts=[_ACTIVE_STATE], output=_ACTIVE_STATE, color="#668d75", fields=[
        field("binder_percent_start", "Starting binder percentage", 10.0, minimum=0.001, maximum=100, section="basic"),
        field("binder_percent_target", "Final binder percentage", 0.01, minimum=0.001, maximum=100, section="basic"),
        field("tightening_rounds", "Rounds to reach final percentage", 10, minimum=1, maximum=100, section="basic"),
    ]),
    "active_train": _spec("active_train", "Train OpenVS-style neural network", "Active learning", "Train a binary classifier on 1024-bit Morgan fingerprints using the OpenVS dense architecture: three hidden layers of 3000 nodes and dropout 0.5.", accepts=[_ACTIVE_STATE], output=_ACTIVE_STATE, color="#a88845", fields=[
        field("hidden_nodes", "Nodes per hidden layer", 3000, minimum=1, maximum=20000, section="basic"),
        field("dropout", "Dropout", 0.5, minimum=0, maximum=0.9, section="basic"),
        field("epochs", "Training epochs", 20, minimum=1, maximum=500, section="basic"),
        field("training_batch_size", "Training batch size", 5000, minimum=1, maximum=100000, section="basic"),
    ]),
    "active_predict": _spec("active_predict", "Predict compound library", "Active learning", "Stream through the registered fingerprint shards and write only compound ID plus model probability; do not copy the whole library into Studio's molecule database.", accepts=[_ACTIVE_STATE], output=_ACTIVE_STATE, color="#8268a2", fields=[
        field("prediction_batch_size", "Prediction batch size", 10000, minimum=1, maximum=100000, section="basic"),
    ]),
    "active_select": _spec("active_select", "Select predicted hits and random compounds", "Active learning", "Select high-probability undocked compounds plus a uniform random exploration sample from the remaining external library. The combined next-round batch is a percentage of the actual initial centroid sample.", accepts=[_ACTIVE_STATE], output=_ACTIVE_STATE, color="#8268a2", fields=[
        field("next_round_seed_fraction_percent", "Next-round batch (% of initial centroid sample)", 50.0, minimum=1, maximum=100, section="basic"),
        field("predicted_fraction_percent", "Next round from top predictions (%)", 80.0, minimum=0, maximum=100, section="basic"),
    ]),
    "active_decision": _spec("active_decision", "Stopping decision", "Active learning", "After accumulating and reordering all docking rounds, compare the top N cumulative scores with the previous round's cumulative top N. Compare the score at each rank using the score units declared by the docking node. Stop after the selected number of stable comparisons, when the library is exhausted, or at the optional iteration limit.", accepts=[_ACTIVE_STATE], output=_ACTIVE_STATE, color="#bd8d42", fields=[
        field("top_n", "Top ranks compared", 100, minimum=1, maximum=100000, section="basic"),
        field("score_tolerance", "Stable score change", 0.05, minimum=0, maximum=100, section="basic", help="Use the units declared by the docking node for its selected score field."),
        field("stable_rounds", "Stable rounds before stopping", 3, minimum=1, maximum=100, section="basic"),
        field("max_iterations", "Maximum iterations (0 = unlimited)", 0, minimum=0, maximum=100000, section="basic"),
    ]),
    "active_results": _spec("active_results", "Active learning results", "Results", "All compound scores, IDs, SMILES, round history, models and retained poses stay in external campaign files. Intermediate rounds never import compound records into Database. Optionally import a bounded selection from the final cumulative ranking after the loop stops.", output=None, accepts=[_ACTIVE_STATE], color="#427f68", formats=(), fields=[
        field("final_database_import_limit", "Final compounds imported into Database (0 = none)", 1000,
              minimum=0, maximum=1000, section="basic",
              help="Import only this many top-ranked final compounds, with scores, original IDs and available retained poses. All other results remain linked as external Parquet files, never as individual Database rows."),
    ]),
})
NODE_CATALOG["active_docking"]["inputs"][0]["multiple"] = True
NODE_CATALOG["active_predict"]["inputs"][0]["multiple"] = True
NODE_CATALOG["active_select"]["inputs"][0]["multiple"] = True
NODE_CATALOG["complexa_generator"]["inputs"][0]["optional"] = True

from mn_molecular_design_studio.quality_tools import TOOLS as QUALITY_TOOLS
for _key, _profile in QUALITY_TOOLS.items():
    _fields = []
    if _key in {"sequence_qc", "pep_patch"}:
        _fields.append(field("ph", "pH", 7.4, minimum=0., maximum=14.))
    if _key in {"biotite_interface", "interface_clustering"}:
        _fields.append(field("contact_distance", "Contact cutoff (Å)", 5., minimum=.1, maximum=20.))
    if _profile.get("cohort"):
        _fields.append(field("similarity_threshold", "Similarity threshold", .8 if _key == "sequence_clustering" else .5, minimum=.01, maximum=1.))
    if _key == "foldseek":
        _fields.append(field("coverage", "Alignment coverage", .8, minimum=.01, maximum=1.))
    NODE_CATALOG[_key] = _spec(_key, _profile["name"], "Analysis", "Analyze one collection as a batch; retain all candidates and their per-item outcomes.",
        accepts=["sequence", "prediction", "candidates"], output="prediction", workflow="analysis." + _key, fields=_fields, color="#9a7854")


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
    "boltz2": [field("generate_ipsae", "Calculate ipSAE", True, kind="boolean")],
    "af2": [SEED, field("af2_validation_mode", "AF2 validation mode", "conditioned_self_consistency", kind="select", section="basic",
        options=["conditioned_self_consistency", "target_only"], help="Conditioned: initialize with the designed complex. Target-only: retain target template/sequence/side chains and supply binder sequence without binder coordinates. Neither mode measures experimental binding probability."),
        field("model_count", "Number of models", 1, minimum=1, maximum=5, section="basic"),
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
    elif _name == "sequence_stage":
        _mpnn_fields.append(field(_name, "Sequence-design stage", _default, kind="select",
            options=["redesign", "initial", "soluble_optimization"], section="basic",
            help="Initial: ProteinMPNN for backbone-only candidates, preserve co-designed sequences. Soluble optimization: optional interface-fixed redesign, retaining parents. Trained on soluble proteins; no experimental solubility guarantee."))
    elif _name == "backend":
        _mpnn_fields.append(field(_name, "Execution backend", _default, kind="select", options=["auto", "shared_ligandmpnn", "foundry"]))
    else:
        _default = 0 if _name == "seed" else "" if _default is None or isinstance(_default, dict) else _default
        _mpnn_fields.append(field(
            _name,
            "Fix interface residues from valid generator sequences" if _name == "fix_interface_residues" else _name.replace("_", " ").capitalize(),
            _default,
            kind="boolean" if isinstance(_default, bool) else "number" if isinstance(_default, (float, int)) else "text",
            help="Fix binder residues with a target heavy-atom contact at the interface cutoff (Å). Requires explicit chain roles and matching sequence/structure; ambiguous mappings are rejected." if _name == "fix_interface_residues" else None,
        ))
for _kind in ("ligandmpnn", "foundry_mpnn"):
    _specification = NODE_CATALOG[_kind]
    _specification.update(label="ProteinMPNN", workflow="sequence.mpnn", formats=["pdb", "cif"], fields=copy.deepcopy(_mpnn_fields))
    _specification["inputs"][0]["accepts"] = ["backbone", "sequence", "prediction"]
    _specification["inputs"][0]["multiple"] = True
NODE_CATALOG["foundry_mpnn"]["hidden"] = True
NODE_CATALOG["foundry_mpnn"]["fields"].append(field("number_of_batches", "Number of batches", 1, minimum=1, maximum=1000))
for _kind in ("sabr", "abmpnn"):
    NODE_CATALOG[_kind]["inputs"][0]["multiple"] = True

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
NODE_CATALOG["rfdiffusion_beta_pairing"] = _spec(
    "rfdiffusion_beta_pairing",
    "RFdiffusion β-pairing",
    "Design engine",
    "Generate binder backbones with explicit target–binder β-strand SS/adjacency interface conditioning. Connect ProteinMPNN and AlphaFold2 Initial Guess after this step.",
    accepts=["target"],
    output="backbone",
    workflow="design.rfdiffusion_beta_pairing",
    formats=("pdb",),
    fields=[
        field("target_pdb", "Target PDB (optional when connected)", "", kind="text", section="basic"),
        field("target_chain", "Target chain (blank = connected chain)", "", kind="text", section="basic"),
        field("target_beta_strands", "Target β-strand residues", "", kind="text", section="basic", help="Use ranges such as B99-103,B120-124. Each continuous segment gets a corresponding designed binder β-strand."),
        field("binder_length", "Binder length", 65, minimum=30, maximum=500, section="basic"),
        field("binder_strand_length", "Designed β-strand length", 6, minimum=3, maximum=40, section="basic"),
        field("num_designs", "Number of backbones", 10, minimum=1, maximum=10000, section="basic"),
        field("timesteps", "Diffusion timesteps", 50, minimum=1, maximum=200),
        field("systematic_sampling", "Systematic register sampling", False, kind="boolean"),
        SEED,
    ],
    color="#27674e",
)
NODE_CATALOG["rfdiffusion_beta_pairing"]["inputs"][0]["optional"] = True
NODE_CATALOG["rfdiffusion_beta_pairing"]["inputs"][0]["label"] = "Target (optional)"
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

# Generic collection sockets use the existing normalized-candidate transport.
NODE_CATALOG.update({
    "candidate_merge": _spec("candidate_merge", "Candidate collection", "Collection", "Merge native outputs without replacing sequences or deleting scientific histories.",
        accepts=["backbone", "sequence", "prediction", "candidates"], output="candidates", fields=[
            field("deduplication", "Sequence deduplication", "mark_duplicates", kind="select", options=["keep_all", "mark_duplicates"], section="basic")]),
    "shared_filter": _spec("shared_filter", "Shared scientific filters", "Collection", "Keep every decision; route sequence-ready, eligible variants to prediction.",
        accepts=["candidates", "sequence", "prediction"], output="sequence", fields=[
            field("rules", "Scientific threshold rules", [], kind="json", section="basic"),
            field("exclude_native_rejected", "Exclude native-rejected designs from independent prediction", True, kind="boolean", section="basic")]),
    "candidate_analyzer": _spec("candidate_analyzer", "Unified candidate analyzer", "Collection", "Join independent predictor evidence to exact sequence variants, with explicit missing measurements.",
        accepts=["candidates", "sequence", "prediction"], output="candidates", fields=[
            field("calculate_self_consistency", "Calculate target-aligned binder pose RMSD when a reference exists", False, kind="boolean", section="basic")]),
    "candidate_rank": _spec("candidate_rank", "Ranking and diversity", "Collection", "Rank declared measurements within comparable protocols; retain unselected and unresolved records.",
        accepts=["candidates"], output="candidates", fields=[
            field("rules", "Scientific threshold rules", [], kind="json", section="basic"),
            field("ranking_rules", "Ordered ranking measurements", [], kind="json", section="basic"),
            field("limit", "Selected candidates per comparable protocol (0 = unlimited)", 0, minimum=0, maximum=10000, section="basic"),
            field("unique_sequences", "Select one representative per exact sequence and protocol", True, kind="boolean", section="basic")]),
})
for _kind in ("candidate_merge", "candidate_analyzer"):
    NODE_CATALOG[_kind]["inputs"][0]["multiple"] = True
for _kind, _specification in NODE_CATALOG.items():
    if _specification["category"] in {"Sequence design", "Structure prediction"}:
        for _socket in _specification["inputs"]:
            if "candidates" not in _socket["accepts"]:
                _socket["accepts"].append("candidates")
    if _specification['category'] == 'Sequence design' and not any(f['name'] == 'sequence_stage' for f in _specification['fields']):
        _specification['fields'].append(field('sequence_stage', 'Sequence-design stage', 'redesign', kind='select',
            options=['redesign', 'initial'], help='Initial designs backbone-only inputs and preserves existing native sequences. Redesign explicitly creates additional variants.'))

TEMPLATE_REFERENCES={
    'multi_generator_binder': ('https://www.anthropic.com/research/claude-scorpion-toxin-binder-design', 'General composition inspired by Figure 2: complete native generators, conditional initial ProteinMPNN, independent predictors and consolidated scientific evidence. No paper-specific target, campaign, thresholds or autonomous agent is reproduced.'),
    'peptide_trap': ('https://doi.org/10.1038/s41589-024-01578-5', 'Published peptide-trap scaffolds → burial-pattern matching and threading → Rosetta interface redesign → AlphaFold2 Initial Guess → Rosetta interface analysis. Uses the authors’ archived models and an adaptation of their deposited FastDesign/FastRelax protocol. C37 and the complete C104.1 model are available initially. Initial scaffold generation, hairpin retargeting and the paper’s complete selection procedure are not reproduced.'),
    'kejia_idr': ('https://github.com/drhicks/Kejia_peptide_binders', 'Kejia template threading and flexible peptide docking → ProteinMPNN → AlphaFold2 multimer with a whole-complex initial guess. This starter uses five sequences per backbone and five AF2 models. Native Kejia threading is reused; downstream Studio implementations replace the upstream cluster scripts. Iterative redesign, clustering and the modified diffusion refinement branch are not included automatically.'),
    'unidock_active_learning':('https://github.com/Novartis/OpenVS','Active-learning study template inspired by OpenVS: an independent centroid sample, configurable docking engine, cumulative score ranking, changing binder percentage, a 1024-bit/3000-node classifier, prediction across the compound library, and selection of predicted hits plus random exploration.'),
    'abmpnn_cdr_redesign':('https://app.tamarind.bio/tools/abmpnn','Select saved project designs, annotate antibody CDRs with SAbR, redesign selected loops with AbMPNN, then independently refold the redesigned sequences with Boltz-2 before sequence analysis and Results.'),
    'rfd4_protein_binder':('https://huggingface.co/spaces/nvidia/rfd4-proteina-zerogpu','RFD4-Proteina co-designs sequence and structure together. This starter sends the generated complex to Boltz-2, then calculates ipSAE and Rosetta interface metrics. Use a protein target structure; edit the RFD4 task options for exact condition controls.'),
    'rfd4_small_molecule':('https://github.com/RosettaCommons/RFD4-Proteina/blob/main/docs/examples/python_api_tutorial.ipynb','Small-molecule pocket design based on RFD4’s annotated-CIF Python tutorial. The RFD4 node accepts either guided ligand controls or an annotated CIF/native condition specification. The connected Boltz-2 step independently predicts the generated protein–ligand complex.'),
    'rfd4_enzyme':('https://github.com/RosettaCommons/RFD4-Proteina/blob/main/docs/examples/python_api_tutorial.ipynb','Enzyme motif scaffolding based on the annotated-CIF tutorial. Specify motif residues, exact active-site atoms and ligand interaction roles in RFD4 task options; Boltz-2 provides an independent structure-prediction step.'),
    'pxdesign':('https://github.com/bytedance/PXDesign#-quick-start','Binder Design Context feeds supported target-site annotations into PXDesign, followed by ProteinMPNN and independent AF2 Initial Guess and Protenix evaluation branches.'),
    'rfdiffusion':('https://github.com/RosettaCommons/RFdiffusion#sequence-design-for-binders','Binder Design Context supplies supported target-site annotations, then ProteinMPNN designs sequences followed by target-templated AF2 Initial Guess. The authors also permit ProteinMPNN without FastRelax; this template uses that option.'),
    'genie3':('https://github.com/aqlaboratory/genie3#application-3-binder-design','Binder Design Context supplies supported target-site annotations. One ProteinMPNN sequence per backbone is followed by AF2 target-template evaluation with five models and twenty recycles. Studio uses its AF2 Initial Guess adapter; the authors use ColabFold template mode.'),
    'protpardelle':('https://github.com/ProteinDesignLab/protpardelle-1c#sampling','Structure generation, eight ProteinMPNN sequences per backbone, then ESMFold self-consistency evaluation. The starter task is unconditional generation; switch to a conditioned task in the generator form.'),
    'esmfold2_binder':('https://github.com/Biohub/esm/tree/main/cookbook/tutorials','ESMFold2 Binder Design performs gradient-guided sequence optimization with its built-in ESMFold2 and ESMC scoring. Its candidate output can be sent through ProteinMPNN and then any compatible independent structure predictor; edit or replace these downstream starter steps on the canvas.'),
    'all_generators_boltz2':('https://github.com/RosettaCommons/RFdiffusion','Miniprotein binder workflow: branch from one Miniprotein Design Context into each available generator, use compact 60–100 residue designs where that engine exposes a binder-length range, use ProteinMPNN only for backbone-only outputs, preserve native co-designed sequences, optionally create interface-fixed SolubleMPNN variants alongside parents, then predict with Boltz-2 and calculate ipSAE and Rosetta interface metrics. SolubleMPNN is trained on soluble proteins and does not guarantee experimental solubility. Each generator receives only annotations it supports.'),
    'vhh_design':('https://github.com/Biohub/esm/tree/main/cookbook/tutorials','VHH-focused workflow using only engines with an explicit VHH or nanobody mode: BindCraft2 VHH scaffold with CDR edit regions, BoltzGen nanobody task, and ESMFold2 Binder Design in VHH mode with a selectable framework. Their original sequences pass directly to Boltz-2 and interface analysis. Use the SAbR/AbMPNN CDR-redesign workflow when explicit redesign is needed. Hotspot and coldspot annotations are passed only to compatible generator paths.'),
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


ACTIVE_LEARNING_TEMPLATE_VERSION = 7
ACTIVE_LEARNING_LAYOUT = {
    "active_seed": {"x": 80, "y": 40, "width": 320, "height": 150},
    "active_library": {"x": 760, "y": 40, "width": 320, "height": 150},
    "active_target": {"x": 455, "y": 215, "width": 230, "height": 110},
    "active_docking": {"x": 80, "y": 330, "width": 320, "height": 180},
    "active_select": {"x": 760, "y": 330, "width": 320, "height": 180},
    "active_accumulate": {"x": 80, "y": 610, "width": 320, "height": 180},
    "active_predict": {"x": 760, "y": 610, "width": 320, "height": 180},
    "active_decision": {"x": 80, "y": 870, "width": 320, "height": 240},
    "active_train": {"x": 760, "y": 900, "width": 320, "height": 180},
    "active_binder_rule": {"x": 455, "y": 930, "width": 230, "height": 120},
    "active_results": {"x": 80, "y": 1220, "width": 320, "height": 130},
}


def new_graph(template="foundry") -> dict[str, Any]:
    if template == "multi_generator_binder":
        from mn_molecular_design_studio.candidate_pipeline import new_pipeline
        return clean_graph(new_pipeline())
    from mn_molecular_design_studio.scientific_tasks import TASKS, new_task_graph
    if template in TASKS:
        return new_task_graph(template)
    if template == "unidock_active_learning":
        kinds = ["active_library", "active_seed", "active_docking", "active_accumulate",
                 "active_binder_rule", "active_train", "active_predict", "active_select",
                 "active_decision", "active_results", "active_target"]
        labels = {"active_library": "Compound library", "active_seed": "Centroid library",
                  "active_docking": "Docking", "active_accumulate": "Accumulate results",
                  "active_decision": "Stopping criteria reached?", "active_binder_rule": "Binder labels",
                  "active_train": "Train NN", "active_predict": "Predict compound library",
                  "active_select": "Select next docking batch", "active_target": "Receptor",
                  "active_results": "Final results"}
        nodes = []
        for index, kind in enumerate(kinds):
            position = ACTIVE_LEARNING_LAYOUT[kind]
            spec = NODE_CATALOG[kind]
            params = {item["name"]: copy.deepcopy(item["default"]) for item in spec["fields"]}
            nodes.append({"id": f"{kind}-{index}", "type": "studio", "position": {"x": position["x"], "y": position["y"]},
                          "data": {"kind": kind, "label": labels[kind], "params": params}})
        by_kind = {node["data"]["kind"]: node["id"] for node in nodes}
        pairs = [("active_seed", "active_docking"),
                 ("active_target", "active_docking"), ("active_docking", "active_accumulate"),
                 ("active_accumulate", "active_decision"), ("active_decision", "active_binder_rule"),
                 ("active_binder_rule", "active_train"),
                 ("active_train", "active_predict"), ("active_predict", "active_select"),
                 ("active_library", "active_predict"), ("active_library", "active_select"),
                 ("active_decision", "active_results")]
        edges = [{"id": f"edge-active-{index}", "source": by_kind[source], "sourceHandle": "out",
                  "target": by_kind[target], "targetHandle": "in"} for index, (source, target) in enumerate(pairs)]
        edges.append({"id": "edge-active-loop", "source": by_kind["active_select"], "sourceHandle": "out",
                      "target": by_kind["active_docking"], "targetHandle": "in", "loop_back": True})
        return clean_graph({"schema_version": SCHEMA, "id": uuid.uuid4().hex,
                            "name": "Active-learning docking screen",
                            "nodes": nodes, "edges": edges,
                            "viewport": {"x": 0, "y": 0, "zoom": 0.62},
                            "template_id": template, "template_version": ACTIVE_LEARNING_TEMPLATE_VERSION})
    if template == "abmpnn_cdr_redesign":
        kinds = ["project_design_source", "sabr", "abmpnn", "boltz2", "sequence_qc", "results"]
        positions = [(60, 190), (370, 190), (680, 190), (990, 190), (1300, 80), (1610, 190)]
        nodes = []
        for index, (kind, position) in enumerate(zip(kinds, positions)):
            spec = NODE_CATALOG[kind]
            params = {f["name"]: copy.deepcopy(f["default"]) for f in spec["fields"]}
            if kind == "project_design_source":
                params.update(source_run_dir="", selected_candidate_ids=[])
            nodes.append({
                "id": f"{kind}-{index}", "type": "studio",
                "position": {"x": position[0], "y": position[1]},
                "data": {"kind": kind, "label": spec["label"], "params": params},
            })
        edges = [
            {"id": "edge-sabr", "source": nodes[0]["id"], "sourceHandle": "out", "target": nodes[1]["id"], "targetHandle": "in"},
            {"id": "edge-abmpnn", "source": nodes[1]["id"], "sourceHandle": "out", "target": nodes[2]["id"], "targetHandle": "in"},
            {"id": "edge-refold", "source": nodes[2]["id"], "sourceHandle": "out", "target": nodes[3]["id"], "targetHandle": "in"},
            {"id": "edge-qc", "source": nodes[3]["id"], "sourceHandle": "out", "target": nodes[4]["id"], "targetHandle": "in"},
            {"id": "edge-results", "source": nodes[4]["id"], "sourceHandle": "out", "target": nodes[5]["id"], "targetHandle": "in"},
        ]
        return clean_graph({
            "schema_version": SCHEMA, "id": uuid.uuid4().hex,
            "name": "Antibody CDR redesign · SAbR → AbMPNN → sequence analysis",
            "nodes": nodes, "edges": edges,
            "viewport": {"x": 0, "y": 0, "zoom": 0.72},
            "template_id": template,
        })
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
    if template in {"rfd4_protein_binder", "rfd4_small_molecule", "rfd4_enzyme"}:
        task, label, description = {
            "rfd4_protein_binder": (
                "protein_binder", "RFD4 Protein Binder",
                "Co-design a protein binder, then run Boltz-2, ipSAE and Rosetta interface analysis.",
            ),
            "rfd4_small_molecule": (
                "small_molecule_binder", "RFD4 Small-Molecule Pocket Design",
                "Design a protein around a ligand in an mmCIF complex, then independently predict the protein–ligand complex with Boltz-2. Atom restraints are prefilled for the tutorial’s 5ENI/5QA example; edit them for another structure.",
            ),
            "rfd4_enzyme": (
                "enzyme", "RFD4 Enzyme Motif Scaffolding",
                "Scaffold the tutorial’s 5ENI/5QA active-site motif around its ligand, then independently predict the resulting complex with Boltz-2. Update residue and atom names if using another structure.",
            ),
        }[template]
        kinds = ["target"]
        if task == "protein_binder":
            kinds.append("binder_context")
        kinds.extend(["rfd4", "boltz2"])
        if task == "protein_binder":
            kinds.extend(["ipsae", "rosetta_metrics"])
        kinds.append("results")

        nodes = []
        for index, kind in enumerate(kinds):
            spec = NODE_CATALOG[kind]
            params = {f["name"]: copy.deepcopy(f["default"]) for f in spec["fields"]}
            if kind == "rfd4":
                params.update(task=task, input_mode="guided", num_designs=8)
                if task == "small_molecule_binder":
                    params.update(length_min=100, length_max=155)
                    params["options_json"] = json.dumps({
                        "ligand_pose": "generated",
                        "ligand_atoms_to_hold": "5QA:CL,5QA:CL1",
                        "ligand_interactions": [
                            {"interaction": "Hydrogen-bond donor", "ligand_atoms": "5QA:N"},
                            {"interaction": "Hydrogen-bond acceptor", "ligand_atoms": "5QA:O1,5QA:O"},
                        ],
                    }, indent=2)
                elif task == "enzyme":
                    params.update(length_min=150, length_max=200)
                    params["options_json"] = json.dumps({
                        "motif_residues": "A41,A36,A87",
                        "motif_atoms": "all",
                        "unindex_motif": True,
                        "ligand_pose": "generated",
                        "active_site_atoms": "A36:CG1 A87:OG1",
                        "ligand_atoms_to_hold": "5QA:CL,5QA:CL1",
                        "ligand_hbond_donors": "5QA:N",
                        "ligand_hbond_acceptors": "5QA:O1,5QA:O",
                    }, indent=2)
                else:
                    params["options_json"] = "{}"
            if kind == "boltz2":
                params["write_full_pae"] = True
            nodes.append({
                "id": f"{kind}-{index}",
                "type": "studio",
                "position": {"x": 0, "y": 0},
                "data": {"kind": kind, "label": spec["label"], "params": params},
            })

        node_by_kind = {node["data"]["kind"]: node for node in nodes}
        if task == "protein_binder":
            for kind, x in (("target", 40), ("binder_context", 340), ("rfd4", 700),
                            ("boltz2", 1030), ("ipsae", 1360), ("rosetta_metrics", 1690), ("results", 2020)):
                node_by_kind[kind]["position"] = {"x": x, "y": 180}
            node_by_kind["binder_context"]["data"].update(
                label="Binder Design Context",
                description="Set target chains and optional hotspots. RFD4 receives the selected structure and hotspot annotations.",
            )
        else:
            positions = [(40, 180), (390, 180), (740, 180), (1090, 180)]
            for node, position in zip(nodes, positions):
                node["position"] = {"x": position[0], "y": position[1]}

        pairs = list(zip(nodes, nodes[1:]))
        edges = [
            {"id": f"edge-rfd4-{index}", "source": source["id"], "sourceHandle": "out",
             "target": destination["id"], "targetHandle": "in"}
            for index, (source, destination) in enumerate(pairs)
        ]
        return clean_graph({
            "schema_version": SCHEMA,
            "id": uuid.uuid4().hex,
            "name": label + " · RFD4 → independent prediction" + (" + interface analysis" if task == "protein_binder" else ""),
            "nodes": nodes,
            "edges": edges,
            "viewport": {"x": 0, "y": 0, "zoom": 0.5 if task == "protein_binder" else 0.7},
            "template_id": template,
        })
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
                params["sequence_stage"] = "initial"
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
        downstream[0]["data"]["label"] = "Initial sequences · ProteinMPNN or native"
        soluble = copy.deepcopy(downstream[0])
        soluble["id"] = "optional-soluble-mpnn"
        soluble["data"].update(label="Optional SolubleMPNN optimization")
        soluble["data"]["params"].update(sequence_stage="soluble_optimization", model_type="soluble_mpnn",
            optimization_enabled=False, fix_interface_residues=True, backend="shared_ligandmpnn")
        soluble["position"]["x"] += 300
        nodes.insert(nodes.index(downstream[0]) + 1, soluble)
        for item in downstream[1:]:
            item["position"]["x"] += 330
        next(e for e in edges if e["source"] == downstream[0]["id"])["source"] = soluble["id"]
        edges.append({"id": "edge-initial-soluble", "source": downstream[0]["id"], "sourceHandle": "out", "target": soluble["id"], "targetHandle": "in"})
        result = {
            "schema_version": SCHEMA, "id": uuid.uuid4().hex,
            "name": "Generator comparison → initial sequences → optional SolubleMPNN → Boltz-2 + analysis",
            "nodes": nodes, "edges": edges,
            "viewport": {"x": 0, "y": 0, "zoom": 0.42},
            "template_id": template, "template_version": 4,
        }
        return clean_graph(result)
    if template == "vhh_design":
        generator_kinds = ["bindcraft2_generator", "boltzgen_generator", "esm_binder"]
        kinds = ["target", "binder_context", *generator_kinds, "sabr", "boltz2", "ipsae", "rosetta_metrics", "results"]
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
            edges.append({"id": f"edge-vhh-predict-{index}", "source": node["id"], "sourceHandle": "out", "target": downstream[0]["id"], "targetHandle": "in"})
        for index, (source, dest) in enumerate(zip(downstream, downstream[1:])):
            edges.append({"id": f"edge-vhh-analysis-{index}", "source": source["id"], "sourceHandle": "out", "target": dest["id"], "targetHandle": "in"})
        return clean_graph({
            "schema_version": SCHEMA, "id": uuid.uuid4().hex,
            "name": "VHH design · supported generators → Boltz-2 + analysis",
            "nodes": nodes, "edges": edges,
            "viewport": {"x": 0, "y": 0, "zoom": 0.62},
            "template_id": template, "template_version": 2,
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
    if engine == "peptide_trap":
        kinds = ["generation_source" if template.endswith("_generated") else "peptide_trap", "af2", "rosetta_metrics", "results"]
    elif engine == "kejia_idr":
        kinds = ["backbone_source" if template.endswith("_generated") else "kejia_idr", "ligandmpnn", "af2", "results"]
    elif engine == "rfdiffusion_beta_pairing":
        kinds = ["backbone_source" if template.endswith("_generated") else engine, "ligandmpnn", "af2", "results"]
    elif engine in generator_input.ENGINES:
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
    elif engine in {"peptide_trap", "kejia_idr", "rfdiffusion_beta_pairing"}:
        for node, x in zip(nodes, (50, 370, 690, 1010)):
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
    if engine == "peptide_trap":
        next(node for node in nodes if node["data"]["kind"] == "af2")["data"]["params"].update(
            model_count=1, num_recycles=3, multimer=True, use_binder_template=False,
            use_interface_template=False, msa_source="single_sequence")
    if engine == "kejia_idr":
        for node in nodes:
            if node["data"]["kind"] == "ligandmpnn":
                node["data"]["params"].update(sequences=5, sampling_temp=0.1, omit_aas="CX")
            elif node["data"]["kind"] == "af2":
                node["data"]["params"].update(model_count=5, num_recycles=3, multimer=True,
                    use_binder_template=False, use_interface_template=False,
                    msa_source="single_sequence")
    if engine == "rfdiffusion_beta_pairing":
        sequence = next(n for n in nodes if n["data"]["kind"] == "ligandmpnn")
        sequence["data"]["params"]["sequences"] = 10
        af2 = next((n for n in nodes if n["data"]["kind"] == "af2"), None)
        if af2:
            af2["data"]["params"].update(
                model_count=1,
                num_recycles=3,
                multimer=True,
                use_binder_template=False,
                use_interface_template=False,
                msa_source="single_sequence",
            )
    result={"schema_version": SCHEMA, "id": uuid.uuid4().hex, "name": {"foundry": "Foundry binder design", "bindcraft2": "BindCraft2 design", "boltzgen": "BoltzGen design", "complexa":"Proteina-Complexa custom pipeline", "complexa_generated":"Continue generated designs", "esmfold2_binder":"ESMFold2 Binder Design → native sequence → RosettaFold3", "rfdiffusion_beta_pairing":"RFdiffusion β-pairing → ProteinMPNN → AF2 Initial Guess", "rfdiffusion_beta_pairing_generated":"Continue β-pairing backbones → ProteinMPNN → AF2 Initial Guess"}.get(template, generator_input.ENGINES[engine]['title']+' · author-inspired pipeline' if engine in generator_input.ENGINES else "Untitled pipeline"), "nodes": nodes, "edges": edges, "viewport": {"x": 0, "y": 0, "zoom": 0.75}}
    if engine == "kejia_idr":
        result["name"] = "Kejia IDR → ProteinMPNN → AlphaFold2 Initial Guess"
    elif engine == "peptide_trap":
        result["name"] = "Peptide Trap → AlphaFold2 → interface analysis"
    for node in nodes:
        if node["data"]["kind"] in {"ligandmpnn", "foundry_mpnn"}:
            node["data"]["params"]["sequence_stage"] = "initial"
            node["data"]["label"] = "Initial sequences · ProteinMPNN or native"
    if template in TEMPLATE_REFERENCES:result['template_id']=template
    return result


def migrate_template_graph(graph: dict[str, Any]) -> dict[str, Any]:
    """Migrate supported starter drafts while retaining explicit user settings."""
    if graph.get("template_id") == "all_generators_boltz2" and graph.get("template_version") == 3:
        stages = [n for n in graph.get("nodes", []) if n.get("data", {}).get("kind") == "ligandmpnn"]
        defaults = {f["name"]: f["default"] for f in NODE_CATALOG["ligandmpnn"]["fields"]}
        defaults["sequences"] = 1
        if len(stages) == 1 and "sequence_stage" not in stages[0]["data"].get("params", {}) and all(
            value == defaults.get(key) for key, value in stages[0]["data"].get("params", {}).items()
        ):
            migrated = copy.deepcopy(graph)
            initial = next(n for n in migrated["nodes"] if n["id"] == stages[0]["id"])
            initial["data"]["params"]["sequence_stage"] = "initial"
            soluble = next(n for n in new_graph("all_generators_boltz2")["nodes"] if n["id"] == "optional-soluble-mpnn")
            soluble["id"] = "soluble-" + uuid.uuid4().hex[:8]
            soluble["position"] = {"x": initial["position"]["x"] + 300, "y": initial["position"]["y"]}
            for edge in migrated["edges"]:
                if edge["source"] == initial["id"]:
                    edge["source"] = soluble["id"]
            migrated["nodes"].append(soluble)
            migrated["edges"].append({"id": "initial-" + soluble["id"], "source": initial["id"], "target": soluble["id"], "sourceHandle": "out", "targetHandle": "in"})
            migrated["template_version"] = 4
            return migrated
        return graph
    if graph.get("template_id") == "vhh_design" and graph.get("template_version") in (None, 1):
        # Repair only the unmodified, implicit generic redesign in old starter drafts.
        # Explicitly configured redesign nodes remain the user's choice.
        candidates = [node for node in graph.get("nodes", []) if node.get("data", {}).get("kind") == "ligandmpnn"]
        if len(candidates) == 1:
            node = candidates[0]
            defaults = {item["name"]: item["default"] for item in NODE_CATALOG["ligandmpnn"]["fields"]}
            defaults["sequences"] = 1
            params = node["data"].get("params") or {}
            if all(value == defaults.get(key) for key, value in params.items()):
                migrated = copy.deepcopy(graph)
                replacement = next(item for item in migrated["nodes"] if item["id"] == node["id"])
                specification = NODE_CATALOG["sabr"]
                replacement["data"] = {"kind": "sabr", "label": specification["label"],
                    "params": {item["name"]: copy.deepcopy(item["default"]) for item in specification["fields"]}}
                migrated["template_version"] = 2
                return migrated
        return graph
    if graph.get("template_id") == "unidock_active_learning":
        if graph.get("template_version") == ACTIVE_LEARNING_TEMPLATE_VERSION:
            return graph
        canonical = new_graph("unidock_active_learning")
        existing = {node["data"]["kind"]: node for node in graph.get("nodes", [])}
        if set(existing) != set(ACTIVE_LEARNING_LAYOUT):
            return graph
        migrated = copy.deepcopy(graph)
        seed = next(node for node in migrated["nodes"] if node["data"]["kind"] == "active_seed")
        library = existing["active_library"]["data"]["params"]
        if not seed["data"]["params"].get("centroid_manifest") and library.get("centroid_manifest"):
            seed["data"]["params"]["centroid_manifest"] = library["centroid_manifest"]
            seed["data"]["params"]["library_id"] = library.get("library_id", "")
        docking_params = next(node for node in migrated["nodes"]
                              if node["data"]["kind"] == "active_docking")["data"]["params"]
        template_version = int(graph.get("template_version") or 1)
        # Older active-learning drafts stored one CPU count directly on the
        # docking node. Worker resources now belong to each selected computer;
        # retaining the old field makes strict graph validation reject an
        # otherwise runnable saved draft.
        docking_params.pop("cpu_cores", None)
        if template_version < 5 and docking_params.get("batch_size") == 10_000:
            # Before v5 this was the fixed batch size for every round.
            docking_params["batch_size"] = 0
        elif template_version < 6 and docking_params.get("batch_size") == 1_000_000:
            # v5 introduced a default cap. Keep user-selected smaller caps,
            # but remove the old default so rounds follow the centroid sample.
            docking_params["batch_size"] = 0
        ids = {}
        for reference in canonical["nodes"]:
            node = next(n for n in migrated["nodes"] if n["data"]["kind"] == reference["data"]["kind"])
            ids[reference["id"]] = node["id"]
            node["position"] = copy.deepcopy(reference["position"])
            old_defaults = {"Annotated library", "Annotated compound library", "THRα receptor", "THRα receptor and docking box", "UniDock Pro · normal mode", "Predict across annotated library", "Predict across external library"}
            if node["data"].get("label") == NODE_CATALOG[node["data"]["kind"]]["label"] or node["data"].get("label") in old_defaults:
                node["data"]["label"] = reference["data"]["label"]
        migrated["edges"] = [{**edge, "source": ids[edge["source"]], "target": ids[edge["target"]]}
                             for edge in canonical["edges"]]
        migrated["template_version"] = ACTIVE_LEARNING_TEMPLATE_VERSION
        if migrated.get("name") == "THRα · UniDock Pro active-learning screen":
            migrated["name"] = canonical["name"]
        return migrated
    if graph.get("template_id") != "all_generators_boltz2" or graph.get("template_version") in (3, 4):
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
    return clean_graph(migrate_template_graph(migrated))


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
        if kind == 'boltz2':
            from mn_molecular_design_studio.prediction_configuration import boltz_settings
            params = boltz_settings(params, saved_graph=True)
        if kind == 'af2' and 'af2_validation_mode' not in params:
            upstream_kinds = {n.get('id'): n.get('data', {}).get('kind') for n in nodes}
            direct = any(e.get('target') == node['id'] and upstream_kinds.get(e.get('source')) == 'cofold_input' for e in edges)
            params['af2_validation_mode'] = 'target_only' if direct else 'conditioned_self_consistency'
        allowed = {f["name"] for f in NODE_CATALOG[kind]["fields"]} | ({"selected_candidate_ids"} if kind in {"design_source","generation_source","backbone_source","project_design_source"} else set())
        if kind in {"bindcraft2", "rfd4"}:
            allowed.add("distributed_campaign")
        if kind == 'boltz2':
            allowed |= {'prediction_config_version', 'configuration_migration', 'write_full_pae'}
        # Migrate generation-only kinds from the previous small canvas forms.
        if kind=='rfdiffusion':params.pop('rfdiffusion_guidance_preset',None)
        if kind=='pxdesign':
            params.pop('use_deepspeed_evo_attention',None)
            params.pop('prepare_target_msa',None)
        if set(params) - allowed:
            raise ValueError(f"Unsupported settings for {NODE_CATALOG[kind]['label']}: {', '.join(sorted(set(params)-allowed))}.")
        if "distributed_campaign" in params:
            campaign = params["distributed_campaign"]
            if kind not in {"bindcraft2", "rfd4"} or not isinstance(campaign, dict) or campaign.get("enabled") is not True:
                raise ValueError("Distributed campaign settings are available only when enabled on a supported BindCraft2 or RFD4 tile.")
            from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import validate_campaign_config
            campaign = {**campaign, "generator_id": "bindcraft2" if kind == "bindcraft2" else "rfd4_proteina"}
            params["distributed_campaign"] = validate_campaign_config(campaign)
            if kind == "rfd4" and (params.get("task") != "protein_binder" or params.get("input_mode") not in {"guided", "guided_binder"}):
                raise ValueError("RFD4 distributed campaigns support guided protein-binder tiles only.")
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
        if "enabled" in data or "retained_run_dir" in data:
            if NODE_CATALOG[kind]["category"] not in {"Design engine", "Design pipeline", "Structure prediction"}:
                raise ValueError("Branch controls apply to generator and predictor nodes only.")
            if type(data.get("enabled", True)) is not bool or not isinstance(data.get("retained_run_dir", ""), str):
                raise ValueError("Branch enabled must be boolean and retained_run_dir must be a path string.")
            node_data.update(enabled=data.get("enabled", True), retained_run_dir=data.get("retained_run_dir", ""))
        if "execution" in data:
            from mn_molecular_design_studio.candidate_execution import policy
            node_data["execution"] = policy(data["execution"])
        if description and description.strip():
            node_data["description"] = description.strip()[:1000]
        cleaned.append({"id": node["id"], "type": "studio", "position": position, "data": node_data})
    viewport = graph.get("viewport") or {"x": 0, "y": 0, "zoom": 1}
    if not isinstance(viewport, dict) or any(isinstance(viewport.get(k), bool) or not isinstance(viewport.get(k), (float, int)) or not math.isfinite(viewport[k]) for k in ("x", "y", "zoom")) or not .05 <= viewport["zoom"] <= 5:
        raise ValueError("Invalid canvas viewport.")
    kinds_by_id = {node["id"]: node.get("data", {}).get("kind") for node in nodes}
    return {"schema_version": SCHEMA, "id": graph["id"], "name": str(graph.get("name") or "Untitled pipeline")[:160], "nodes": cleaned,
            "edges": [{k: edge[k] for k in ("id", "source", "target", "sourceHandle", "targetHandle") if k in edge}
                      | ({"loop_back": True} if edge.get("loop_back") is True or (
                          graph.get("template_id") == "unidock_active_learning"
                          and kinds_by_id.get(edge.get("source")) == "active_select"
                          and kinds_by_id.get(edge.get("target")) == "active_docking"
                      ) else {})
                      for edge in edges if isinstance(edge, dict)],
            **({"candidate_execution": __import__("mn_molecular_design_studio.candidate_execution", fromlist=["policy"]).policy(graph["candidate_execution"])} if "candidate_execution" in graph else {}),
            "viewport": viewport,**({'template_id':graph['template_id']} if graph.get('template_id') in SUPPORTED_TEMPLATE_IDS else {}),
            **({'task_configuration': copy.deepcopy(graph['task_configuration'])}
               if isinstance(graph.get('task_configuration'), dict) else {}),
            **({'template_version':graph['template_version']} if type(graph.get('template_version')) is int and (
                graph.get('template_id') == 'all_generators_boltz2' and 1 <= graph['template_version'] <= 4
                or graph.get('template_id') == 'vhh_design' and graph['template_version'] == 2
                or graph.get('template_id') == 'unidock_active_learning' and 2 <= graph['template_version'] <= ACTIVE_LEARNING_TEMPLATE_VERSION) else {})}


def validate_active_library_sources(screening: dict, seed: dict) -> tuple[dict, dict]:
    """Independent datasets share the model's fingerprint contract, not an ID namespace."""
    from mn_molecular_design_studio.external_libraries import read_feature_manifest
    if not screening.get("library_id") or not screening.get("feature_manifest"):
        raise ValueError("Choose the compound screening library and its fingerprint dataset.")
    if not seed.get("library_id") or not seed.get("centroid_manifest"):
        raise ValueError("Choose the initial centroid library and its centroid dataset independently.")
    try:
        _path, features = read_feature_manifest(screening["feature_manifest"])
        centroids = json.loads(Path(seed["centroid_manifest"]).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError("A selected fingerprint or centroid manifest is unavailable or invalid.") from exc
    if centroids.get("schema_version") != "studio.library-centroids.v2" or centroids.get("status") != "completed":
        raise ValueError("Choose a completed centroid manifest.")
    if features.get("library_id") != screening["library_id"]:
        raise ValueError("The fingerprint dataset must belong to its selected screening library.")
    if centroids.get("library_id") != seed["library_id"]:
        raise ValueError("The centroid dataset must belong to its selected centroid source library.")
    left, right = features.get("fingerprint_spec", {}), centroids.get("fingerprint_spec", {})
    required = ("algorithm", "bits", "radius", "use_chirality", "use_features")
    if any(key not in spec for spec in (left, right) for key in required):
        raise ValueError("Both datasets must declare their Morgan fingerprint settings.")
    if left["bits"] != 1024 or right["bits"] != 1024:
        raise ValueError("Both model inputs require 1024-bit Morgan fingerprints.")
    differences = [key for key in (*required, "hydrogen_policy", "standardization") if left.get(key) != right.get(key)]
    if differences:
        raise ValueError("Centroid and screening fingerprints must use compatible settings: " + ", ".join(differences) + ".")
    return features, centroids


def validate_graph(graph: dict, *, require_inputs=False) -> dict:
    graph = clean_graph(migrate_template_graph(graph))
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
        if edge.get("loop_back") is True:
            if graph.get("template_id") != "unidock_active_learning" or a["key"] != "active_select" or b["key"] != "active_docking":
                raise ValueError("Only next-batch selection can return compounds to docking.")
            continue
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
        from mn_molecular_design_studio.candidate_sets import KINDS, validate_settings
        if kind in KINDS:
            validate_settings(kind, nodes[node_id]["data"]["params"])
        if kind == "filter":
            from mn_molecular_design_studio.candidate_selection import selection_policy
            selection_policy(nodes[node_id]["data"]["params"])
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
        from mn_molecular_design_studio.candidate_execution import policy, supports
        execution = policy(node["data"].get("execution", graph.get("candidate_execution")))
        if (execution["mode"] == "sharded" and "execution" in node["data"]
                and not supports(kind, {"candidates_jsonl": "preflight"}) and kind not in KINDS
                and kind not in {"target", "binder_context", "results"}):
            raise ValueError(f"{node['data']['label']} has no collection-sharding adapter; use native execution.")
        if node["data"].get("enabled") is False:
            continue
        if NODE_CATALOG[kind]["inputs"] and not incoming[node_id] and not NODE_CATALOG[kind]["inputs"][0].get("optional"):
            raise ValueError(f"Connect an input to {node['data']['label']}.")
        if kind not in {"results", "active_results"} and not any(e["source"] == node_id for e in graph["edges"]):
            raise ValueError(f"Connect {node['data']['label']} to a downstream step or Results node.")
        if require_inputs and kind == "target":
            if not node["data"]["params"].get("target_pdb"):
                raise ValueError("Load a target structure before running the pipeline.")
        if require_inputs and kind == "kejia_idr":
            from mn_molecular_design_studio.kejia_idr_input import prepare_request
            prepare_request(node["data"]["params"])
        if require_inputs and kind == "peptide_trap":
            from mn_molecular_design_studio.peptide_trap_input import prepare_request
            prepare_request(node["data"]["params"])
        if require_inputs and kind == "active_library":
            seed = next((n for n in nodes.values() if n["data"]["kind"] == "active_seed"), None)
            if seed is not None:
                validate_active_library_sources(node["data"]["params"], seed["data"]["params"])
        if require_inputs and kind == "active_target":
            target_path = str(node["data"]["params"].get("target_pdb") or "").strip()
            if not target_path or not Path(target_path).expanduser().is_file():
                raise ValueError("Choose an existing imported receptor structure.")
        if require_inputs and kind == 'cofold_input':
            params = node['data']['params']
            task = params.get('task_type')
            if not params.get('target_protein_sequences', '').strip() and not params.get('sequence_file'):
                raise ValueError('Enter at least one target protein sequence in Cofolding Inputs.')
            target_count = len([part for part in re.split(r'[:\n\r]+', params.get('target_protein_sequences', '')) if part.strip()])
            partner_count = len([part for part in re.split(r'[:\n\r]+', params.get('partner_protein_sequences', '')) if part.strip()])
            if task == 'protein' and target_count != 1 and not params.get('sequence_file'):
                raise ValueError('The protein monomer task accepts exactly one protein sequence.')
            if task in {'protein_protein', 'protein_macrocyclic'} and not params.get('partner_protein_sequences', '').strip():
                raise ValueError('Enter at least one partner protein sequence for this cofolding task.')
            if task == 'protein_macrocyclic' and partner_count != 1:
                raise ValueError('The protein–macrocyclic peptide task accepts exactly one cyclic partner sequence.')
            if task == 'protein_ligand' and not (params.get('compound_set_file') or params.get('ligand_ccd_codes', '').strip() or params.get('ligand_smiles', '').strip()):
                raise ValueError('Enter a CCD code or SMILES string for the protein–ligand task.')
            ligand_format = params.get('ligand_input_type', 'smiles')
            if task == 'protein_ligand' and not params.get('compound_set_file') and ligand_format == 'ccd' and not params.get('ligand_ccd_codes', '').strip():
                raise ValueError('Enter a CCD code for the selected ligand notation.')
            if task == 'protein_ligand' and not params.get('compound_set_file') and ligand_format == 'smiles' and not params.get('ligand_smiles', '').strip():
                raise ValueError('Enter a SMILES string for the selected ligand notation.')
            if task == 'protein_dna' and not params.get('dna_sequences', '').strip():
                raise ValueError('Enter a DNA sequence for the protein–DNA task.')
            if task == 'protein_rna' and not params.get('rna_sequences', '').strip():
                raise ValueError('Enter an RNA sequence for the protein–RNA task.')
            if task == 'protein_macrocyclic' and not params.get('cyclic_partner'):
                raise ValueError('Mark the protein partner as cyclic for the macrocyclic task.')
            optional_target_path = str(params.get('target_pdb') or '').strip()
            for selection_key in ('compound_set_file', 'sequence_file'):
                if params.get(selection_key) and not Path(params[selection_key]).expanduser().is_file():
                    raise ValueError(f'The {selection_key} collection in Cofolding Inputs does not exist.')
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
        scaffold_mode = kind in {"bindcraft2", "bindcraft2_generator"} and node["data"]["params"].get("modality") in {"VHH", "scFv", "Fab", "ARP"}
        if isinstance(length, str) and not (length == "" and scaffold_mode):
            match = re.fullmatch(r"(\d+)(?:-(\d+))?", length)
            if not match or not 1 <= int(match[1]) <= int(match[2] or match[1]) <= 10000:
                raise ValueError("Binder length must be a positive length or range, such as 74 or 60-100.")
        if require_inputs and kind == "design_source" and not node["data"]["params"].get("source_run_dir"):
            raise ValueError("Select a successful design job from this project.")
        if require_inputs and kind == "molecule_source" and not node["data"]["params"].get("source_run_dir"):
            raise ValueError("Select molecules in Database and choose Use in pipeline.")
        if require_inputs and kind == "generation_source" and not node["data"]["params"].get("source_run_dir"):
            raise ValueError("Select a completed generator job from this project.")
        if require_inputs and kind == "project_design_source":
            if not node["data"]["params"].get("source_run_dir"):
                raise ValueError("Select a completed design job from this project.")
            if not node["data"]["params"].get("selected_candidate_ids"):
                raise ValueError("Select one or more design candidates from that job.")
        if require_inputs and kind=='backbone_source' and not node['data']['params'].get('source_run_dir'):
            raise ValueError('Select a completed generator job from this project.')
        if require_inputs and kind in {"sabr", "abmpnn"} and not incoming[node_id]:
            raise ValueError(f"Connect antibody designs to {NODE_CATALOG[kind]['label']}.")
        if require_inputs and kind == "abmpnn" and not any(node["data"]["params"].get(key) for key in ("redesign_cdr1", "redesign_cdr2", "redesign_cdr3")):
            raise ValueError("Select at least one CDR region for AbMPNN redesign.")
        if require_inputs and kind in generator_input.ENGINES and not incoming[node_id]:
            p=node['data']['params']
            if p.get('task') not in {'unconditional','symmetry','custom','fold'} and not p.get('target_pdb'):
                raise ValueError('Load a structure in '+node['data']['label']+' or connect a Target block.')
        if require_inputs and kind == "rfdiffusion_beta_pairing":
            p = node["data"]["params"]
            if not incoming[node_id] and not p.get("target_pdb"):
                raise ValueError("Load a PDB target in RFdiffusion β-pairing or connect a Target block.")
            if not str(p.get("target_beta_strands") or "").strip():
                raise ValueError("Select target β-strand residue ranges for RFdiffusion β-pairing.")
    if not any(n["data"]["kind"] in {"results", "active_results"} for n in nodes.values()):
        raise ValueError("Add a Results node to collect the pipeline outputs.")
    if graph.get("template_id") == "unidock_active_learning":
        required = {"active_library", "active_target", "active_seed", "active_docking", "active_accumulate",
                    "active_binder_rule", "active_train", "active_predict", "active_select", "active_decision", "active_results"}
        present = {node["data"]["kind"] for node in nodes.values()}
        if present != required:
            raise ValueError("The active-learning template needs its complete set of source, docking, learning, decision and Results nodes.")
        if sum(edge.get("loop_back") is True for edge in graph["edges"]) != 1:
            raise ValueError("The active-learning template needs exactly one selection-to-docking feedback connection.")
    return {"graph": graph, "order": order, "incoming": incoming, "outgoing": outgoing}


def fingerprint(kind: str, params: dict, upstream: list[str]) -> str:
    from mn_molecular_design_studio.cache_provenance import referenced_files
    return hashlib.sha256(json.dumps({"schema": SCHEMA, "cache_contract": 2, "kind": kind, "params": params,
        "referenced_contents": referenced_files(params), "upstream": sorted(upstream)}, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


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


NODE_CATALOG["rfd4"]["inputs"][0]["optional"] = True
NODE_CATALOG["rfd4_fold"]["inputs"][0]["optional"] = True

WORKFLOW_ROLES = {
    "design.rfd4_proteina": "Design engine",
    "refolding.rfd4_proteina": "Structure prediction engine",
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


# Shared by the canvas and REST/MCP template discovery.
LEGACY_STARTER_TEMPLATES = {"Active-learning docking": "unidock_active_learning",
                 "Foundry · RFD3 → MPNN → RF3": "foundry",
                 "Kejia IDR · Template threading → ProteinMPNN → AF2": "kejia_idr",
                 "Peptide Trap · Rosetta interface design → AF2 → analysis": "peptide_trap",
                 "Antibody CDRs · Existing designs → SAbR → AbMPNN → sequence analysis": "abmpnn_cdr_redesign",
                 "PXDesign · MPNN → AF2 + Protenix":"pxdesign","RFdiffusion · MPNN → AF2":"rfdiffusion",
                 "Genie3 · MPNN → AF2 template":"genie3","Protpardelle-1C · MPNN → ESMFold":"protpardelle",
                 "ESMFold2 Binder Design → ProteinMPNN → RF3":"esmfold2_binder",
                 "Generator comparison · initial sequences → optional SolubleMPNN → Boltz-2 + analysis":"all_generators_boltz2",
                 "VHH binders · BindCraft2 + BoltzGen + ESMFold2 → Boltz-2":"vhh_design",
                 "RFD4 Protein Binder → Boltz-2 → interface analysis":"rfd4_protein_binder",
                 "RFD4 Small-Molecule Pocket Design → Boltz-2":"rfd4_small_molecule",
                 "RFD4 Enzyme Motif Scaffolding → Boltz-2":"rfd4_enzyme",
                 "Cofolding · Protein monomer":"cofold_protein",
                 "Cofolding · Protein / protein":"cofold_protein_protein",
                 "Cofolding · Protein / ligand":"cofold_protein_ligand",
                 "Cofolding · Protein / DNA":"cofold_protein_dna",
                 "Cofolding · Protein / RNA":"cofold_protein_rna",
                 "Cofolding · Protein / macrocyclic peptide":"cofold_protein_macrocyclic",
                 "Blank canvas": "blank"}

# Navigation changes do not migrate or reinterpret saved graphs. All former IDs,
# hidden native recipes and generated-result handoffs remain constructible.
from mn_molecular_design_studio.scientific_tasks import TASKS
STARTER_TEMPLATES = {spec['label']: key for key, spec in TASKS.items()} | {"Multi-Generator Binder Design & Evaluation": "multi_generator_binder", "Blank canvas": "blank"}
SUPPORTED_TEMPLATE_IDS = (set(LEGACY_STARTER_TEMPLATES.values()) | set(STARTER_TEMPLATES.values())
                          | {'bindcraft2', 'boltzgen', 'complexa', 'rfdiffusion_beta_pairing'}
                          | {key + '_generated' for key in ('complexa', 'rfdiffusion', 'genie3', 'pxdesign',
                                                           'protpardelle', 'kejia_idr', 'peptide_trap', 'rfdiffusion_beta_pairing')})
