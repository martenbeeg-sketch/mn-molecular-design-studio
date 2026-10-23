"""Allowlisted, schema-checked CLI access to app workflows.

Workflow IDs are the public local-automation names. They map to existing app
workflow functions and job groups; the CLI never imports or executes a function
named by an untrusted request.
"""

from __future__ import annotations

import importlib
import inspect
import json
import math
import types
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Union, get_args, get_origin, get_type_hints

from mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler import submission_priority
from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import spawn_worker_for_run
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import enqueue_workflow_call
from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import (
    AutomationError,
    job_reference,
    resolve_job,
)


WORKFLOW_SCHEMA_VERSION = 1
_schema_only_paths = ContextVar("studio_schema_only_paths", default=False)


@contextmanager
def deferred_path_validation():
    """Check request types before remote files are staged; execution validates them."""
    token = _schema_only_paths.set(True)
    try:
        yield
    finally:
        _schema_only_paths.reset(token)


def _request_priority(payload: dict[str, Any]) -> int:
    priority = payload.get("priority", 0)
    if type(priority) is not int or not 0 <= priority <= 1000:
        raise AutomationError("priority must be an integer from 0 to 1000.")
    return priority


_INTERNAL_PARAMETERS = frozenset(
    {
        "existing_job",
        "_defer_finish",
        "progress_callback",
        "internal_parent_run_dir",
        "internal_parent_task_group",
        "internal_parent_role",
        "internal_parent_engine",
        "run_dir",
        "resume",
        "input_zip_bytes",
        "table_bytes",
        "task_group",
        "job_type",
        "tool_name",
        "direct_candidates",
        "direct_input",
    }
)
_COMMON_REFOLDING_FIELDS = frozenset(
    {
        "models", "max_records", "num_loops", "num_sampling_steps", "seed", "device", "docker_image",
        "run_pyrosetta_input_metrics", "pyrosetta_nprocs", "run_common_interface_metrics",
        "run_predicted_rosetta_metrics", "run_pymol_metrics", "run_esmfold2", "esmfold2_modes",
        "esmfold2_use_target_msa", "run_af2_initial_guess", "af2_num_recycles", "af2_multimer",
        "af2_use_initial_guess", "af2_use_binder_template", "af2_use_interface_template",
        "run_boltz2_initial_guess", "boltz2_use_target_template", "boltz2_use_target_msa",
        "boltz2_recycling_steps", "boltz2_sampling_steps", "boltz2_diffusion_samples",
        "boltz2_write_full_pae", "run_rf3", "rf3_checkpoint_path", "rf3_use_target_msa",
        "rf3_use_target_template", "rf3_recycles", "rf3_num_steps", "rf3_diffusion_batch_size",
        "rf3_seed", "run_openfold3", "openfold3_checkpoint_path", "openfold3_use_target_msa",
        "openfold3_num_diffusion_samples", "openfold3_num_model_seeds", "openfold3_num_recycles",
        "openfold3_use_msa_server", "run_protenix", "protenix_use_msa", "protenix_cycle",
        "protenix_diffusion_steps", "protenix_samples", "run_protenix_v1", "protenix_v1_model_name",
        "protenix_v1_use_msa", "protenix_v1_use_template", "protenix_v1_use_default_params",
        "protenix_v1_cycle", "protenix_v1_diffusion_steps", "protenix_v1_samples", "run_protenix_v2",
        "protenix_v2_model_name", "protenix_v2_use_msa", "protenix_v2_use_template",
        "protenix_v2_use_default_params", "protenix_v2_cycle", "protenix_v2_diffusion_steps",
        "protenix_v2_samples", "run_boltzgen_fold", "boltzgen_recycling_steps", "boltzgen_sampling_steps",
        "boltzgen_diffusion_samples", "run_colabfold", "colabfold_cache_dir", "colabfold_msa_source",
        "msa_repository_dir", "require_real_target_msa", "colabfold_num_recycles", "colabfold_num_models",
        "colabfold_use_target_templates", "colabfold_use_target_msa", "colabfold_max_template_hits",
        "colabfold_gpu_device", "run_alphafast_af3", "alphafast_db_dir", "alphafast_weights_dir",
        "alphafast_batch_size", "alphafast_num_recycles", "alphafast_use_target_templates",
        "alphafast_query_only_msa", "alphafast_gpu_device", "gpu_device", "target_override_pdb",
        "target_override_chains",
    }
)


@dataclass(frozen=True)
class WorkflowSpec:
    workflow_id: str
    module: str
    function: str
    task_group: str
    tool: str
    job_type: str
    description: str
    mode: str = "queued"


def _spec(
    workflow_id: str,
    module: str,
    function: str,
    group: str,
    tool: str,
    job_type: str,
    description: str,
    mode: str = "queued",
) -> WorkflowSpec:
    return WorkflowSpec(workflow_id, f"mn_molecular_design_studio.engine.workflow_engine.workflows.{module}", function, group, tool, job_type, description, mode)


_SPECS = (
    _spec("ligand.complex.prepare", "ligand_tools", "run_ligand_complex_prepare", "ligand-preparation", "studio_ligand_complex_preparation", "ligand_complex_preparation", "Import and repair a protein-ligand complex with optional C-terminal extension and chain-aware trimming."),
    _spec("ligand.docking.vina", "ligand_tools", "run_ligand_vina", "ligand-docking", "vina", "ligand_docking", "Dock a Database compound set against a target pocket with AutoDock Vina."),
    _spec("ligand.docking.autodock_gpu", "ligand_tools", "run_ligand_autodock_gpu", "ligand-docking", "autodock_gpu", "ligand_docking", "Dock a Database compound set with AutoDock-GPU and AutoGrid affinity maps."),
    _spec("ligand.docking.gnina", "ligand_tools", "run_ligand_gnina", "ligand-docking", "gnina", "ligand_docking", "Dock a Database compound set with GNINA and CNN pose rescoring."),
    _spec("ligand.docking.unidock_pro", "ligand_tools", "run_ligand_unidock_pro", "ligand-docking", "unidock_pro", "ligand_docking", "GPU-accelerated classic receptor docking with Uni-Dock Pro."),
    _spec("ligand.docking.rosetta_ligand", "ligand_tools", "run_ligand_rosetta_ligand", "ligand-docking", "rosetta_ligand", "ligand_docking", "RosettaLigand docking with VSH, VSX, or convergence protocols."),
    _spec("ligand.affinity.nesso", "ligand_tools", "run_ligand_nesso", "ligand-affinity", "nesso_affinity", "ligand_affinity", "Predict protein–ligand affinity with Nesso-1; this tile returns affinity scores and does not generate a complex structure."),
    _spec("ligand.generation.lddm", "ligand_tools", "run_ligand_lddm", "ligand-generation", "lddm_generation", "ligand_generation", "Generate pocket-conditioned ligands with the LDDM CD+BB checkpoint."),
    _spec("ligand.generation.pocketxmol", "ligand_tools", "run_ligand_pocketxmol", "ligand-generation", "pocketxmol_generation", "ligand_generation", "Generate de novo ligands in a protein pocket with PocketXMol."),
    _spec("ligand.generation.omtra", "ligand_tools", "run_ligand_omtra", "ligand-generation", "omtra_generation", "ligand_generation", "Generate ligands with OMTRA's multi-condition foundation model."),
    _spec("ligand.generation.flowr_root", "ligand_tools", "run_ligand_flowr_root", "ligand-generation", "flowr_root_generation", "ligand_generation", "Generate or redesign pocket-conditioned ligands with FLOWR.root."),
    _spec("ligand.generation.conditar", "ligand_tools", "run_ligand_conditar", "ligand-generation", "conditar_generation", "ligand_generation", "Generate pocket-conditioned ligands with conDitar."),
    _spec("ligand.generation.paopt", "ligand_tools", "run_ligand_paopt", "ligand-generation", "paopt_generation", "ligand_generation", "Steer conDitar ligand generation toward selected ADMET objectives with paOPT."),
    _spec("ligand.generation.drugrpg", "ligand_tools", "run_ligand_drugrpg", "ligand-generation", "drugrpg_generation", "ligand_generation", "Generate pocket-conditioned ligands with DrugRPG."),
    _spec("ligand.generation.pfm", "ligand_tools", "run_ligand_pfm", "ligand-generation", "pfm_generation", "ligand_generation", "Generate pocket-conditioned ligands with perturbed flow matching (PFM)."),
    _spec("ligand.generation.pocketflow", "ligand_tools", "run_ligand_pocketflow", "ligand-generation", "pocketflow_generation", "ligand_generation", "Generate ligands inside a protein pocket with PocketFlow."),
    _spec("ligand.generation.pgmg", "ligand_tools", "run_ligand_pgmg", "ligand-generation", "pgmg_generation", "ligand_generation", "Generate compounds from a PGMG-compatible pharmacophore hypothesis."),
    _spec("ligand.md.gromacs", "ligand_tools", "run_ligand_gromacs_md", "ligand-md", "gromacs_md", "ligand_md", "Prepare a bound protein–ligand complex and run explicit-solvent GROMACS dynamics."),
    _spec("ligand.md.openmm", "ligand_tools", "run_ligand_openmm_md", "ligand-md", "openmm_md", "ligand_md", "Prepare a bound protein–ligand complex and run explicit-solvent OpenMM dynamics."),
    _spec("ligand.md.analysis", "md_analysis", "run_md_analysis", "md-analysis", "md_analysis", "md_analysis", "Analyze a recorded MD segment with protein-aligned RMSD and residue contact occupancy."),
    _spec("target.prepare", "target_prep", "enqueue_target_preparation", "target-prep", "internal_pdb_cleaner", "target_preparation", "Clean, trim, and optionally prepare MSAs for a PDB target.", "target_prepare"),
    _spec("target.crop", "target_crop", "run_target_crop", "target-crop", "internal_pdb_cropper", "target_cropping", "Crop a prepared target by residue ranges or Mol* selection coordinates."),
    _spec("target.mask", "target_masking", "create_masked_target", "target-prep", "internal_target_masker", "target_masking", "Create a mutation/masked-target job with recorded mutation details.", "direct"),
    _spec("target.refolding_candidates", "target_prep", "create_target_refolding_candidate_set", "target-prep", "target_refolding_input_builder", "target_refolding_candidate_set", "Create a target-only normalized candidate set.", "direct"),
    _spec("target.refolding_evaluation", "target_prep", "enqueue_target_refolding_evaluation", "target-refolding", "target_refolding_evaluation_engines", "target_refolding_evaluation", "Build target-only candidates and run selected refolding engines.", "target_refolding_evaluation"),
    _spec("detection.scannet", "detection", "run_scannet", "detection", "scannet", "ppi_detection", "Run ScanNet PPI/hotspot detection."),
    _spec("detection.pesto", "detection", "run_pesto", "detection", "pesto", "ppi_detection", "Run PeSTo interface detection."),
    _spec("detection.surf2spot", "detection", "run_surf2spot", "detection", "surf2spot", "ppi_detection", "Run Surf2Spot hotspot detection."),
    _spec("detection.masif_seed", "detection", "run_masif_seed", "detection", "masif_seed", "ppi_detection", "Run MaSIF-seed surface/site detection."),
    _spec("detection.surface_analysis", "detection", "run_surface_analysis", "detection", "surface_analysis", "surface_analysis", "Run PeSTo, ScanNet, Surf2Spot, and MaSIF-seed together on one target."),
    _spec("design.rfdiffusion_classic", "design", "run_rfdiffusion_classic", "design", "rfdiffusion_classic", "design_campaign", "Run the RFdiffusion classic design workflow."),
    _spec("design.rfdiffusion_beta_pairing", "rfdiffusion_beta_pairing", "run_rfdiffusion_beta_pairing", "design", "rfdiffusion_beta_pairing", "design_generation", "Generate RFdiffusion binder backbones with explicit target–binder β-strand SS/adjacency interface conditioning."),
    _spec("design.bindcraft", "design", "run_bindcraft", "design", "bindcraft", "design_campaign", "Run the BindCraft 1.x vanilla workflow."),
    _spec("design.bindcraft2", "bindcraft2", "enqueue_bindcraft2_design", "design", "bindcraft2", "bindcraft2_design", "Run a native BindCraft 2 design job.", "direct"),
    _spec("design.bindcraft2_generator", "bindcraft2", "enqueue_bindcraft2_design", "design", "bindcraft2_generator", "bindcraft2_generator", "Run BindCraft2 gradient design and emit only trajectories that pass the pre-ProteinMPNN handoff gate.", "direct"),
    _spec("design.rfdiffusion3", "rfdiffusion3", "run_rfdiffusion3_design", "design", "rfdiffusion3", "rfdiffusion3_design", "All-atom protein, enzyme, nucleic-acid and small-molecule conditioned design with RFdiffusion3."),
    _spec("design.rfd4_proteina", "rfd4_proteina", "run_rfd4_proteina_design", "design", "rfd4_proteina", "rfd4_proteina_design", "Generate sequence and structure together with RFDiffusion4-Proteina using guided binder inputs, annotated CIFs, or native YAML/JSON conditions."),
    _spec("refolding.rfd4_proteina", "rfd4_proteina", "run_rfd4_structure_prediction", "refolding", "rfd4_proteina", "rfd4_structure_prediction", "Predict protein structure with RFD4-Proteina’s own sequence-conditioned model from FASTA or a structure."),
    _spec("design.rfantibody", "rfantibody", "run_rfantibody_design", "design", "rfantibody", "rfantibody_design", "Run RFdiffusion, ProteinMPNN and RF2 antibody or nanobody design."),
    _spec("design.rfantibody_mimosa", "rfantibody", "run_rfantibody_design", "design", "rfantibody_mimosa", "rfantibody_design", "Run paper-based MIMOSA geometry or sequence-fixed guidance with the RFantibody VHH workflow."),
    _spec("design.rfpeptides", "rfpeptides", "run_rfpeptides_design", "design", "rfpeptides", "rfpeptides_design", "Generate cyclic peptide binder backbones with the RFpeptides RFdiffusion protocol."),
    _spec("design.rfdiffusion3_foundry", "design", "run_rfdiffusion3_foundry", "design", "rfdiffusion3_foundry", "design_campaign", "Run RFdiffusion3 / Foundry design."),
    _spec("design.boltzgen", "design", "run_boltzgen", "design", "boltzgen", "design_campaign", "De novo protein binder design with BoltzGen's protein-anything protocol."),
    _spec("design.boltzgen_mimosa", "design", "run_boltzgen", "design", "boltzgen_mimosa", "design_campaign", "Run paper-based MIMOSA geometry or sequence-fixed guidance with BoltzGen nanobody templates."),
    _spec("design.boltzgen_generator", "design", "run_boltzgen", "design", "boltzgen_generator", "design_generation", "Generate BoltzGen backbones only; sequence redesign, refolding, and ranking are separate steps."),
    _spec("design.pxdesign", "design", "run_pxdesign", "design", "pxdesign", "design_campaign", "Run PXDesign."),
    _spec("design.genie3", "design", "run_genie3", "design", "genie3", "design_campaign", "Run Genie3."),
    _spec("design.esmfold2_binder", "esm_binder", "run_esmfold2_native_binder_design", "design", "esmfold2_binder_design", "design_campaign", "Run Biohub ESMFold2 gradient-guided binder design."),
    _spec("design.esmfold2_screening", "esm_binder", "run_esmfold2_binder_screening", "design", "esmfold2_binder", "esmfold2_binder_screening", "Run the app's ESMFold2 binder sequence screening workflow."),
    _spec("design.protpardelle_1c", "design", "run_protpardelle_1c", "design", "protpardelle_1c", "design_campaign", "Run Protpardelle-1c binder generation."),
    _spec("design.proteina_complexa", "design", "run_proteina_complexa", "design", "proteina_complexa", "design_campaign", "Full protein binder, ligand binder, and motif + ligand design pipelines with generation, sequence design, refolding and native filtering."),
    _spec("design.proteina_complexa_generate", "design", "run_proteina_complexa_generation", "design", "proteina_complexa_generator", "design_generation", "Generate protein binders, ligand binders, and motif scaffolds with backbone and sequence. Choose sequence redesign and structure prediction separately on the pipeline canvas."),
    _spec("design.pxdesign_generate", "generators", "run_pxdesign_generation", "design", "pxdesign_generator", "design_generation", "Generate protein binder backbones with PXDesign. Choose sequence redesign and evaluation on the pipeline canvas."),
    _spec("design.rfdiffusion_generate", "generators", "run_rfdiffusion_generation", "design", "rfdiffusion_generator", "design_generation", "Generate or diversify protein backbones with target, motif, fold, or symmetry conditioning."),
    _spec("design.genie3_generate", "generators", "run_genie3_generation", "design", "genie3_generator", "design_generation", "Generate binder, unconditional, or motif-conditioned structures with Genie3."),
    _spec("design.protpardelle_generate", "generators", "run_protpardelle_generation", "design", "protpardelle_generator", "design_generation", "Generate protein structures with sidechain, motif, partial diffusion, and multichain conditioning."),
    _spec("sequence.mpnn", "mpnn", "run_mpnn_sequence_design", "design", "ligandmpnn", "sequence_design", "Design sequences from structures or project designs using ProteinMPNN, LigandMPNN, SolubleMPNN, HyperMPNN, AbMPNN or membrane models."),
    _spec("sequence.abmpnn", "mpnn", "run_abmpnn_sequence_design", "design", "abmpnn", "sequence_design", "Redesign selected antibody CDRs with AbMPNN, using SAbR numbering while keeping the framework and target fixed."),
    _spec("analysis.sabr", "sabr", "run_sabr_annotation", "analysis", "sabr", "antibody_numbering", "Number antibody chains and annotate CDRs with SAbR; optionally save a renumbered structure copy."),
    _spec("sequence.esm_if1", "inverse_folding", "run_esm_if1", "design", "esm_if1", "sequence_design", "Design protein sequences from any structure with ESM-IF1, optionally using other chains as structural context."),
    _spec("sequence.caliby", "inverse_folding", "run_caliby", "design", "caliby", "sequence_design", "Design or score sequences with Caliby or SolubleCaliby, using single structures or structural ensembles."),
    _spec("sequence.boltzgen_inverse", "inverse_folding", "run_boltzgen_inverse", "design", "boltzgen_inverse", "sequence_design", "Use BoltzGen's inverse-folding model on existing structures from any generator or prediction engine."),
    _spec("sequence.ligandmpnn", "sequence_design", "run_ligandmpnn_sequence_design", "design", "ligandmpnn", "sequence_design", "Design sequences for normalized candidates with ProteinMPNN/LigandMPNN."),
    _spec("sequence.foundry_mpnn", "sequence_design", "run_foundry_mpnn_sequence_design", "design", "foundry_mpnn", "sequence_design", "Design sequences with the Foundry MPNN workflow."),
    _spec("sequence.pipeline", "sequence_design", "enqueue_sequence_design_pipeline", "design", "ligandmpnn", "sequence_design", "Run sequence design with an optional queued refolding continuation.", "sequence_pipeline"),
    _spec("refolding.monomer", "refolding", "run_monomer_refolding_contract", "refolding-validation", "af2_monomer", "monomer_refolding", "Run a configured monomer refolding backend."),
    _spec("refolding.esmfold2_monomer", "refolding", "run_esmfold2_monomer_refolding", "refolding-validation", "esmfold2_monomer", "monomer_refolding", "Refold candidate sequences with ESMFold2."),
    _spec("refolding.esmfold2_cofolding", "esmfold2_cofolding", "run_esmfold2_cofolding", "refolding-validation", "esmfold2_cofolding", "complex_refolding", "Cofold protein sequences with the optimized ESMFold2 image."),
    _spec("refolding.af3_jax_cofolding", "af3_jax_cofolding", "run_af3_jax_cofolding", "refolding-validation", "af3_jax_cofolding", "complex_refolding", "Cofold protein and biomolecular inputs with selectable AF3-family JAX models in MN Cofolding's optimized image."),
    _spec("refolding.optimized_prediction", "optimized_prediction", "run_optimized_prediction", "refolding-validation", "optimized_prediction", "complex_refolding", "Run cofolding engines through the shared molecule input and prediction result format."),
    _spec("simulation.bioemu2", "bioemu2", "run_bioemu2", "simulation", "bioemu2", "bioemu2_sampling", "Sample BioEmu-2 conformational ensembles from uploaded target–design complexes or sequence-defined protein systems."),
    _spec("refolding.esmfold_monomer", "refolding", "run_esmfold_monomer_refolding", "refolding-validation", "esmfold", "monomer_refolding", "Refold candidate sequences with ESMFold."),
    _spec("refolding.boltz2_monomer", "refolding", "run_boltz2_monomer_refolding", "refolding-validation", "boltz2_monomer", "monomer_refolding", "Refold candidate sequences with Boltz-2 monomer mode."),
    _spec("refolding.complex", "refolding", "run_complex_refolding_contract", "refolding-validation", "af2_initial_guess", "complex_refolding", "Run a configured complex refolding backend."),
    _spec("refolding.esmfold2_complex", "refolding", "run_esmfold2_complex_validation", "refolding-validation", "esmfold2_complex_validation", "complex_refolding", "Validate candidate complexes with ESMFold2."),
    _spec("refolding.af2_initial_guess", "refolding", "run_af2_initial_guess_complex_refolding", "refolding-validation", "af2_initial_guess", "complex_refolding", "Run AlphaFold2 initial-guess complex refolding."),
    _spec("refolding.af2_initial_guess_input", "refolding", "run_af2_initial_guess_from_inputs", "refolding-validation", "af2_initial_guess", "complex_refolding", "Run AlphaFold2 initial-guess prediction from PDB templates and optional partner sequences."),
    _spec("refolding.af2_design_validation", "design_validation", "run_af2_design_validation", "refolding-validation", "af2_initial_guess", "complex_refolding", "Refold selected project designs with target-only AlphaFold2 templates."),
    _spec("refolding.boltz2_complex", "refolding", "run_boltz2_complex_refolding", "refolding-validation", "boltz2_initial_guess", "complex_refolding", "Run Boltz-2 complex refolding."),
    _spec("refolding.openfold3", "refolding", "run_openfold3_complex_refolding", "refolding-validation", "openfold3", "complex_refolding", "Run OpenFold-3 complex refolding."),
    _spec("refolding.rf3", "refolding", "run_rf3_complex_refolding", "refolding-validation", "rf3", "complex_refolding", "Run RosettaFold3 complex refolding."),
    _spec("refolding.protenix", "refolding", "run_protenix_complex_refolding", "refolding-validation", "protenix", "complex_refolding", "Run Protenix complex refolding."),
    _spec("refolding.protenix_cli", "refolding", "run_protenix_cli_complex_refolding", "refolding-validation", "protenix_v1", "complex_refolding", "Run the Protenix CLI complex-refolding workflow."),
    _spec("refolding.boltzgen_fold", "refolding", "run_boltzgen_fold_complex_refolding", "refolding-validation", "boltzgen_fold", "complex_refolding", "Run BoltzGen Fold complex refolding."),
    _spec("analysis.ipsae", "analysis", "run_ipsae_analysis_contract", "analysis", "ipsae", "analysis", "Calculate standalone ipSAE scores from saved complex structures and pairwise PAE matrices; scalar PAE metrics are not a substitute."),
    _spec("analysis.rosetta", "analysis", "run_rosetta_metrics_contract", "analysis", "rosetta_metrics", "analysis", "Calculate PyRosetta interface and energy metrics from saved complex structures."),
    _spec("import.bindcraft", "candidate_import", "run_bindcraft_import", "candidate-import", "bindcraft_import", "candidate_import", "Import a BindCraft result folder as normalized candidates."),
    _spec("import.table", "candidate_import", "run_generic_table_import", "candidate-import", "generic_table_import", "candidate_import", "Import candidates from CSV/TSV/Excel and referenced structures."),
    _spec("benchmark.dataset", "benchmark", "run_de_novo_binder_scoring_dataset", "benchmark", "de_novo_binder_scoring_scripts", "de_novo_binder_scoring_dataset", "Run the de novo binder scoring dataset and configured structure metrics."),
    _spec("benchmark.refold_candidates", "benchmark", "enqueue_candidate_refolding_evaluation", "benchmark", "refolding_evaluation_engines", "refolding_evaluation", "Queue selected structure-prediction engines for an existing candidate set.", "candidate_refolding"),
    _spec("benchmark.precomputed_metrics", "benchmark", "run_precomputed_metric_benchmark", "benchmark", "de_novo_binder_scoring_metrics", "metric_dataset_benchmark", "Analyze precomputed benchmark metric tables."),
    _spec("benchmark.esmfold2", "benchmark", "run_esmfold2_binder_benchmark", "benchmark", "esmfold2_benchmark", "binder_benchmark", "Run the ESMFold2 binder benchmark."),
    _spec("campaign.validation_sequence", "campaigns", "run_validation_sequence", "refolding-validation", "full_validation_pipeline", "validation_sequence", "Run a configured multi-step refolding/validation lineage."),
    _spec("benchmark.refolding_capacity", "capacity_benchmark", "create_refolding_capacity_benchmark", "benchmark", "refolding_capacity_matrix", "refolding_capacity_benchmark", "Create and launch a refolding capacity matrix.", "direct"),
    _spec("benchmark.design_capacity", "capacity_benchmark", "create_design_capacity_benchmark", "benchmark", "design_generator_capacity_ladder", "design_generator_capacity_benchmark", "Create and launch a design-generator capacity ladder.", "direct"),
    _spec("benchmark.practical_capacity", "capacity_benchmark", "create_practical_capacity_benchmark_from_parent", "benchmark", "design_generator_capacity_ladder", "design_generator_capacity_benchmark", "Create a practical capacity benchmark from an existing parent run.", "direct"),
    _spec("benchmark.collection", "benchmark", "create_benchmark_collection", "benchmark", "benchmark_collection_merge", "benchmark_collection", "Merge metrics from selected completed benchmark runs.", "direct"),
    _spec("benchmark.matrix_workspace", "benchmark", "create_benchmark_matrix_workspace", "benchmark", "benchmark_matrix_workspace", "benchmark_matrix_workspace", "Create a metric-matrix workspace from selected benchmark runs.", "direct"),
)

from mn_molecular_design_studio.quality_tools import TOOLS as _QC_TOOLS
_SPECS = (*_SPECS, *(
    _spec("analysis." + key, "quality", "run_" + key, "analysis", key, "analysis",
        profile["name"] + ": analyze the selected collection in one batch with per-item status and preserved source lineage.")
    for key, profile in _QC_TOOLS.items()))
WORKFLOWS: dict[str, WorkflowSpec] = {spec.workflow_id: spec for spec in _SPECS}


def _get_function(spec: WorkflowSpec):
    module = importlib.import_module(spec.module)
    function = getattr(module, spec.function, None)
    if not callable(function):
        raise AutomationError(f"Workflow implementation is unavailable: {spec.module}.{spec.function}")
    return function


def _type_hints(function) -> dict[str, Any]:
    try:
        return get_type_hints(function)
    except Exception:
        return {}


def _parameters_for(spec: WorkflowSpec) -> dict[str, inspect.Parameter]:
    if spec.workflow_id in {"design.bindcraft2", "design.bindcraft2_generator", "design.boltzgen", "design.boltzgen_generator", "design.boltzgen_mimosa"}:
        function = _get_function(spec)
        if spec.workflow_id in {"design.bindcraft2", "design.bindcraft2_generator"}:
            hidden = {"trajectory_only"}
        elif spec.workflow_id == "design.boltzgen_generator":
            hidden = {"run_vanilla_pipeline", "generator_only_tile"}
        else:
            # The full BoltzGen workflow supports generation-only runs as well
            # as its bundled downstream stages. Keep the internal tile marker
            # private while exposing the actual pipeline toggle to API clients.
            hidden = {"generator_only_tile"}
        return {
            name: parameter
            for name, parameter in inspect.signature(function).parameters.items()
            if name not in hidden
            and name not in _INTERNAL_PARAMETERS
            and parameter.kind not in {inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD}
        }
    if spec.workflow_id == "campaign.validation_sequence":
        function = _get_function(spec)
        signature = inspect.signature(function)
        parameters = {
            name: signature.parameters[name]
            for name in ("campaign_name", "steps", "cpu_cores", "gpu_device")
        }
        parameters["initial_source"] = inspect.Parameter(
            "initial_source",
            inspect.Parameter.KEYWORD_ONLY,
            default=None,
            annotation=dict[str, Any] | None,
        )
        return parameters
    if spec.mode == "target_refolding_evaluation":
        function = _get_function(spec)
        signature = inspect.signature(function)
        return {
            "target_entries": signature.parameters["target_entries"],
            "evaluation_name": signature.parameters["evaluation_name"],
            "options": inspect.Parameter(
                "options",
                inspect.Parameter.KEYWORD_ONLY,
                default={},
                annotation=dict[str, Any],
            ),
        }
    if spec.mode == "sequence_pipeline":
        return {
            "backend": inspect.Parameter("backend", inspect.Parameter.KEYWORD_ONLY, annotation=str),
            "source_run_dir": inspect.Parameter("source_run_dir", inspect.Parameter.KEYWORD_ONLY, annotation=Path),
            "candidates_jsonl": inspect.Parameter("candidates_jsonl", inspect.Parameter.KEYWORD_ONLY, annotation=Path),
            "design_kwargs": inspect.Parameter("design_kwargs", inspect.Parameter.KEYWORD_ONLY, annotation=dict[str, Any]),
            "validation_kwargs": inspect.Parameter("validation_kwargs", inspect.Parameter.KEYWORD_ONLY, default=None, annotation=dict[str, Any] | None),
        }
    if spec.mode == "candidate_refolding":
        signature = inspect.signature(_get_function(spec))
        return {
            "source_run_dir": signature.parameters["source_run_dir"],
            "candidates_jsonl": signature.parameters["candidates_jsonl"],
            "max_candidates": signature.parameters["max_candidates"],
            "selected_candidate_ids": signature.parameters["selected_candidate_ids"],
            "evaluation_name": signature.parameters["evaluation_name"],
            "options": inspect.Parameter(
                "options",
                inspect.Parameter.KEYWORD_ONLY,
                default={},
                annotation=dict[str, Any],
            ),
        }
    function = _get_function(spec)
    return {
        name: parameter
        for name, parameter in inspect.signature(function).parameters.items()
        if name not in _INTERNAL_PARAMETERS
        and parameter.kind not in {inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD}
    }


def _type_name(annotation: Any) -> str:
    if annotation is inspect.Parameter.empty:
        return "any JSON value"
    if annotation is Any:
        return "any JSON value"
    if get_origin(annotation) is Literal:
        return "one of " + ", ".join(repr(item) for item in get_args(annotation))
    return str(annotation).replace("typing.", "").replace("<class '", "").replace("'>", "")


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool, list, dict)):
        return value
    return None


def workflow_catalog(*, category: str | None = None) -> list[dict[str, Any]]:
    """Return public workflow IDs and their request fields."""
    rows: list[dict[str, Any]] = []
    for spec in _SPECS:
        group = spec.workflow_id.split(".", 1)[0]
        if category and group != category:
            continue
        function = _get_function(spec)
        hints = _type_hints(function)
        parameters = []
        for name, parameter in _parameters_for(spec).items():
            annotation = hints.get(name, parameter.annotation)
            item = {
                "name": name,
                "type": _type_name(annotation),
                "required": parameter.default is inspect.Parameter.empty,
            }
            if parameter.default is not inspect.Parameter.empty:
                item["default"] = _json_default(parameter.default)
            parameters.append(item)
        rows.append(
            {
                "workflow": spec.workflow_id,
                "category": group,
                "description": spec.description,
                "parameters": parameters,
            }
        )
    return rows


def workflow_schema(workflow_id: str) -> dict[str, Any]:
    try:
        spec = WORKFLOWS[workflow_id]
    except KeyError as exc:
        raise AutomationError(f"Unknown workflow '{workflow_id}'. Use 'workflow list' to see supported IDs.") from exc
    row = next(item for item in workflow_catalog() if item["workflow"] == workflow_id)
    function = _get_function(spec)
    request_fields = ["schema_version", "workflow", "parameters", "resources"]
    if {"source_run_dir", "candidates_jsonl"}.issubset(inspect.signature(function).parameters) or workflow_id == "campaign.validation_sequence":
        request_fields.append("source_job")
    schema = {"schema_version": WORKFLOW_SCHEMA_VERSION, **row, "request_fields": request_fields}
    schema["collection_execution"] = {
        "automatic": True, "one_job_per_selection": True,
        "parameters_accepts": ["object", "array"] if spec.mode == "queued" else ["object"],
        "strategy": "native_batch_when_available",
    }
    if spec.mode in {"target_refolding_evaluation", "candidate_refolding"}:
        from mn_molecular_design_studio.engine.workflow_engine.workflows.benchmark import run_de_novo_binder_scoring_dataset

        option_fields = set(_COMMON_REFOLDING_FIELDS)
        if spec.mode == "target_refolding_evaluation":
            option_fields.update({"split_chain_breaks", "chain_break_mode"})
        signature = inspect.signature(run_de_novo_binder_scoring_dataset)
        hints = _type_hints(run_de_novo_binder_scoring_dataset)
        options_schema = []
        for name in sorted(option_fields):
            parameter = signature.parameters.get(name)
            if parameter is None:
                if name == "split_chain_breaks":
                    annotation, default = bool, False
                elif name == "chain_break_mode":
                    annotation, default = str, "preserve_original_chain"
                else:
                    continue
            else:
                annotation = hints.get(name, parameter.annotation)
                default = parameter.default
            options_schema.append(
                {
                    "name": name,
                    "type": _type_name(annotation),
                    "required": False,
                    "default": _json_default(default) if default is not inspect.Parameter.empty else None,
                }
            )
        schema["options_schema"] = options_schema
    return schema


def _resolve_path(value: Any, *, request_file: Path, field_name: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise AutomationError(f"{field_name} must be a non-empty path or managed path URI.")
    raw = value.strip()
    if _schema_only_paths.get():
        return Path(raw)
    if raw.startswith(("app://", "reference://", "runs://")):
        path = resolve_stored_path(raw, must_exist=True)
        if path is None:
            raise AutomationError(f"{field_name} does not resolve to an existing path: {raw}")
    else:
        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = request_file.parent / path
        path = path.resolve()
    if not path.exists():
        raise AutomationError(f"{field_name} does not exist: {path}")
    return path.resolve()


def _convert(value: Any, annotation: Any, *, request_file: Path, field_name: str) -> Any:
    if annotation is inspect.Parameter.empty or annotation is Any or annotation is object:
        return value
    if value is None:
        if annotation is type(None) or type(None) in get_args(annotation):
            return None
        raise AutomationError(f"{field_name} may not be null.")
    origin = get_origin(annotation)
    args = get_args(annotation)
    if origin in {Union, types.UnionType}:
        failures = []
        for choice in args:
            if choice is type(None):
                continue
            try:
                return _convert(value, choice, request_file=request_file, field_name=field_name)
            except AutomationError as exc:
                failures.append(str(exc))
        raise AutomationError(f"{field_name} has the wrong type for {_type_name(annotation)}.")
    if origin is Literal:
        if value not in args:
            raise AutomationError(f"{field_name} must be one of {', '.join(map(str, args))}.")
        return value
    if annotation is Path:
        return _resolve_path(value, request_file=request_file, field_name=field_name)
    if annotation is bool:
        if type(value) is not bool:
            raise AutomationError(f"{field_name} must be true or false.")
        return value
    if annotation is int:
        if type(value) is not int:
            raise AutomationError(f"{field_name} must be an integer.")
        return value
    if annotation is float:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise AutomationError(f"{field_name} must be a finite number.")
        return float(value)
    if annotation is str:
        if not isinstance(value, str):
            raise AutomationError(f"{field_name} must be a string.")
        return value
    if origin in {list, tuple, set}:
        if not isinstance(value, list):
            raise AutomationError(f"{field_name} must be a JSON list.")
        item_type = args[0] if args else Any
        converted = [_convert(item, item_type, request_file=request_file, field_name=f"{field_name}[]") for item in value]
        return tuple(converted) if origin is tuple else set(converted) if origin is set else converted
    if origin is dict:
        if not isinstance(value, dict):
            raise AutomationError(f"{field_name} must be a JSON object.")
        key_type, value_type = args if len(args) == 2 else (Any, Any)
        return {
            _convert(key, key_type, request_file=request_file, field_name=f"{field_name} key"): _convert(
                child,
                value_type,
                request_file=request_file,
                field_name=f"{field_name}.{key}",
            )
            for key, child in value.items()
        }
    if isinstance(annotation, type) and isinstance(value, annotation):
        return value
    return value


def _resolve_source_job(
    parameters: dict[str, Any],
    source_job: object,
    function,
    spec: WorkflowSpec,
) -> dict[str, Any]:
    if source_job is None:
        return parameters
    if not isinstance(source_job, str) or not source_job.strip():
        raise AutomationError("source_job must be a visible job ID or code.")
    row = resolve_job(source_job)
    run_dir = Path(row["run_dir"]).resolve()
    candidates = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    if spec.workflow_id == "campaign.validation_sequence" and "initial_source" not in parameters:
        if not candidates.is_file():
            raise AutomationError(f"Source job {row['run_id']} has no normalized candidates file: {candidates}")
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

        metadata = read_json(run_dir / "metadata.json")
        parameters["initial_source"] = {
            "task_group": row.get("task_group") or metadata.get("task_group"),
            "run_id": row.get("run_id") or run_dir.name,
            "run_dir": run_dir,
            "job_code": row.get("job_code") or metadata.get("job_code"),
            "tool": row.get("tool") or metadata.get("tool"),
            "candidates_jsonl": candidates,
            "candidate_count": row.get("candidate_count"),
            "stage_counts": row.get("stage_counts") or {},
        }
        return parameters
    signature = inspect.signature(function)
    if "source_run_dir" in signature.parameters:
        parameters.setdefault("source_run_dir", run_dir)
    if "candidates_jsonl" in signature.parameters:
        if not candidates.is_file():
            if signature.parameters["candidates_jsonl"].default is not None:
                raise AutomationError(f"Source job {row['run_id']} has no normalized candidates file: {candidates}")
        else:
            parameters.setdefault("candidates_jsonl", candidates)
    return parameters


def _resolve_campaign_initial_source(value: Any, request_file: Path) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise AutomationError("initial_source must be an object with run_dir and candidates_jsonl.")
    source = dict(value)
    for field in ("run_dir", "candidates_jsonl"):
        if field not in source:
            raise AutomationError(f"initial_source.{field} is required.")
        source[field] = str(_resolve_path(source[field], request_file=request_file, field_name=f"initial_source.{field}"))
    if not Path(source["run_dir"]).is_dir():
        raise AutomationError("initial_source.run_dir must be a directory.")
    if not Path(source["candidates_jsonl"]).is_file():
        raise AutomationError("initial_source.candidates_jsonl must be a file.")
    return source


def _validate_campaign_steps(value: Any) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not value:
        raise AutomationError("steps must be a non-empty JSON list.")
    supported = {
        ("sequence_design", "ligandmpnn"),
        ("sequence_design", "foundry_mpnn"),
        ("monomer_refolding", "af2_monomer"),
        ("monomer_refolding", "boltz2_monomer"),
        ("monomer_refolding", "esmfold"),
        ("complex_refolding", "af2_initial_guess"),
        ("complex_refolding", "boltz2_initial_guess"),
        ("analysis", "ranking"),
        ("analysis", "filters"),
        ("analysis", "reports"),
    }
    normalized = []
    for index, raw in enumerate(value):
        if not isinstance(raw, dict):
            raise AutomationError(f"steps[{index}] must be an object.")
        unknown = sorted(set(raw) - {"module", "tool", "params"})
        if unknown:
            raise AutomationError(f"Unsupported field(s) in steps[{index}]: {', '.join(unknown)}")
        module, tool = raw.get("module"), raw.get("tool")
        if (module, tool) not in supported:
            raise AutomationError(f"Unsupported campaign step at index {index}: {module}/{tool}.")
        params = raw.get("params") or {}
        if not isinstance(params, dict):
            raise AutomationError(f"steps[{index}].params must be a JSON object.")
        normalized.append({"module": module, "tool": tool, "params": params})
    return normalized


def _validate_generic_parameters(
    spec: WorkflowSpec,
    value: object,
    request_file: Path,
    *,
    source_job: object = None,
) -> tuple[dict[str, Any], Any]:
    if not isinstance(value, dict):
        raise AutomationError("parameters must be a JSON object.")
    function = _get_function(spec)
    if source_job is not None and (
        {"source_run_dir", "candidates_jsonl"}.intersection(value)
        or (spec.workflow_id == "campaign.validation_sequence" and "initial_source" in value)
    ):
        raise AutomationError("Use source_job or explicit source paths, not both.")
    if source_job is not None and not (
        {"source_run_dir", "candidates_jsonl"}.issubset(inspect.signature(function).parameters)
        or spec.workflow_id == "campaign.validation_sequence"
    ):
        raise AutomationError(f"source_job is not supported for {spec.workflow_id}.")
    public_parameters = _parameters_for(spec)
    unknown = sorted(set(value) - set(public_parameters))
    if unknown:
        raise AutomationError(f"Unsupported parameter(s) for {spec.workflow_id}: {', '.join(unknown)}")
    hints = _type_hints(function)
    normalized: dict[str, Any] = {}
    for name, raw in value.items():
        parameter = public_parameters[name]
        annotation = hints.get(name, parameter.annotation)
        normalized[name] = _convert(raw, annotation, request_file=request_file, field_name=name)
    if "target_entries" in normalized:
        normalized["target_entries"] = _resolve_target_entries(normalized["target_entries"], request_file)
    if spec.workflow_id == "campaign.validation_sequence" and "initial_source" in normalized and normalized["initial_source"] is not None:
        normalized["initial_source"] = _resolve_campaign_initial_source(normalized["initial_source"], request_file)
    if spec.workflow_id == "campaign.validation_sequence" and "steps" in normalized:
        normalized["steps"] = _validate_campaign_steps(normalized["steps"])
    if spec.workflow_id == "design.bindcraft2_generator":
        normalized["trajectory_only"] = True
    elif spec.workflow_id == "design.bindcraft2":
        normalized["trajectory_only"] = False
    elif spec.workflow_id == "design.boltzgen_generator":
        normalized["run_vanilla_pipeline"] = False
        normalized["generator_only_tile"] = True
    elif spec.workflow_id == "design.boltzgen":
        normalized.setdefault("run_vanilla_pipeline", True)
        normalized["generator_only_tile"] = False
    if spec.mode in {"target_refolding_evaluation", "candidate_refolding"}:
        _validate_refolding_options(
            normalized.get("options") or {},
            request_file,
            allow_chain_break_options=spec.mode == "target_refolding_evaluation",
        )
    missing = [
        name
        for name, parameter in public_parameters.items()
        if parameter.default is inspect.Parameter.empty
        and name not in normalized
        and not (source_job and name in {"source_run_dir", "candidates_jsonl"})
        and not (source_job and spec.workflow_id == "campaign.validation_sequence" and name == "initial_source")
    ]
    if spec.mode == "sequence_pipeline":
        missing.extend(name for name in ("backend", "source_run_dir", "candidates_jsonl", "design_kwargs") if name not in normalized and name not in missing)
    if spec.mode == "target_refolding_evaluation":
        missing.extend(name for name in ("target_entries",) if name not in normalized and name not in missing)
    if spec.mode == "candidate_refolding":
        missing.extend(name for name in ("source_run_dir", "candidates_jsonl") if name not in normalized and name not in missing)
    if spec.workflow_id == "campaign.validation_sequence" and "initial_source" not in normalized and not source_job:
        missing.append("initial_source or source_job")
    if missing:
        raise AutomationError(f"Missing required parameter(s) for {spec.workflow_id}: {', '.join(missing)}")
    return normalized, function


def _validate_resources(payload: object) -> tuple[int | None, str | None]:
    if payload is None:
        return None, None
    if not isinstance(payload, dict):
        raise AutomationError("resources must be a JSON object.")
    unknown = sorted(set(payload) - {"cpu_cores", "gpu_device"})
    if unknown:
        raise AutomationError(f"Unsupported resource field(s): {', '.join(unknown)}")
    cpu_cores = payload.get("cpu_cores")
    if cpu_cores is not None and (type(cpu_cores) is not int or cpu_cores < 1):
        raise AutomationError("resources.cpu_cores must be a positive integer.")
    gpu_device = payload.get("gpu_device")
    if gpu_device is not None:
        if isinstance(gpu_device, bool) or not isinstance(gpu_device, (str, int)) or not str(gpu_device).strip():
            raise AutomationError("resources.gpu_device must be a device ID such as '0', 'cpu', or 'all'.")
        gpu_device = str(gpu_device).strip()
    return cpu_cores, gpu_device


def _set_resource_parameter(parameters: dict[str, Any], name: str, value: Any, *, label: str) -> None:
    provided = parameters.get(name)
    if provided is not None and str(provided) != str(value):
        raise AutomationError(f"resources.{label} conflicts with parameters.{name}.")
    parameters[name] = value


def _apply_resource_parameters(
    spec: WorkflowSpec,
    function,
    parameters: dict[str, Any],
    *,
    cpu_cores: int | None,
    gpu_device: str | None,
) -> None:
    """Apply resource choices consistently to validated and submitted requests."""
    signature = inspect.signature(function)
    if gpu_device is not None and "gpu_device" in signature.parameters:
        _set_resource_parameter(parameters, "gpu_device", gpu_device, label="gpu_device")
    if gpu_device in {"cpu", "none", "off", "false", "0-gpu"} and "device" in signature.parameters:
        _set_resource_parameter(parameters, "device", "cpu", label="gpu_device")
    if cpu_cores is not None and "cpu_cores" in signature.parameters:
        _set_resource_parameter(parameters, "cpu_cores", cpu_cores, label="cpu_cores")

    if spec.mode in {"target_refolding_evaluation", "candidate_refolding"}:
        options = parameters.setdefault("options", {})
        if gpu_device is not None:
            _set_resource_parameter(options, "gpu_device", gpu_device, label="gpu_device")
            if gpu_device in {"cpu", "none", "off", "false", "0-gpu"}:
                _set_resource_parameter(options, "device", "cpu", label="gpu_device")
        if cpu_cores is not None and options.get("run_pyrosetta_input_metrics"):
            current = options.get("pyrosetta_nprocs")
            if current is None or current == 0:
                options["pyrosetta_nprocs"] = cpu_cores
            elif type(current) is int and current > cpu_cores:
                # Keep explicit PyRosetta parallelism within the scheduler allocation.
                options["pyrosetta_nprocs"] = cpu_cores

    if spec.mode == "sequence_pipeline":
        design_kwargs = parameters.setdefault("design_kwargs", {})
        if gpu_device is not None:
            _set_resource_parameter(design_kwargs, "gpu_device", gpu_device, label="gpu_device")


def _validate_workflow_inputs(spec: WorkflowSpec, parameters: dict[str, Any]) -> None:
    if spec.workflow_id.startswith("ligand."):
        from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import validate_ligand_tool_request

        try:
            validate_ligand_tool_request(
                spec.workflow_id,
                parameters,
                require_files=not _schema_only_paths.get(),
            )
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id.startswith("detection.") and not _schema_only_paths.get():
        from mn_molecular_design_studio.engine.workflow_engine.workflows.detection import validate_surface_target

        try:
            selected = parameters.get("chain_ids")
            if spec.workflow_id == "detection.surface_analysis" and not selected:
                raise ValueError("Select at least one protein chain for surface analysis.")
            if spec.workflow_id == "detection.masif_seed":
                selected = [parameters["chain_id"]]
            chains = validate_surface_target(parameters["target_pdb"], selected or None)
            if spec.workflow_id == "detection.surface_analysis" and parameters.get("masif_chain_id") not in (None, "", *chains):
                raise ValueError("The MaSIF-seed chain must be one of the selected target chains.")
        except (ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    from mn_molecular_design_studio import inverse_folding_input
    if spec.workflow_id in inverse_folding_input.WORKFLOWS.values():
        engine = next(k for k, v in inverse_folding_input.WORKFLOWS.items() if v == spec.workflow_id)
        try:
            inverse_folding_input.prepare(engine, parameters, require_input=True)
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "sequence.mpnn":
        from mn_molecular_design_studio.mpnn_input import prepare
        try:
            prepare(parameters, require_input=True)
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "sequence.abmpnn":
        from mn_molecular_design_studio.mpnn_input import prepare
        from mn_molecular_design_studio.sabr_runtime import CDR_NAMES
        try:
            regions = parameters.get("cdr_regions") or list(CDR_NAMES)
            if not isinstance(regions, list) or not regions or not set(regions) <= set(CDR_NAMES):
                raise ValueError("Choose one or more of CDR1, CDR2, and CDR3.")
            prepare(parameters | {"model_type": "ab_mpnn", "backend": "shared_ligandmpnn"}, require_input=True)
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "analysis.sabr":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.sabr import validate_sabr_request
        try:
            validate_sabr_request(parameters, require_files=not _schema_only_paths.get())
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    from mn_molecular_design_studio import generator_input
    if spec.workflow_id in generator_input.WORKFLOWS.values():
        engine=next(k for k,v in generator_input.WORKFLOWS.items() if v==spec.workflow_id)
        try:
            generator_input.prepare_request(engine,parameters)
        except (TypeError,ValueError,OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "design.rfdiffusion3":
        from mn_molecular_design_studio.rfdiffusion3_input import prepare_request
        try:
            prepare_request(parameters)
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "design.rfd4_proteina":
        from tempfile import TemporaryDirectory
        from mn_molecular_design_studio.rfd4_proteina_input import prepare_request, TASKS
        try:
            if _schema_only_paths.get():
                task = parameters.get("task", "protein_binder")
                mode = parameters.get("input_mode", "guided_binder")
                if task not in TASKS:
                    raise ValueError("Choose a supported RFDiffusion4-Proteina design task.")
                if mode not in {"guided", "guided_binder", "annotated_cif", "condition_spec", "unconditional"}:
                    raise ValueError("Choose a supported RFDiffusion4-Proteina input method.")
            else:
                with TemporaryDirectory(prefix="studio-rfd4-validate-") as temporary:
                    prepare_request(parameters, run_dir=Path(temporary))
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "refolding.rfd4_proteina":
        if not any(parameters.get(k) for k in ("target_structure","protein_sequences","input_files")):
            raise AutomationError("Provide protein sequence(s), a structure, or connected input structures.")
        if parameters.get("target_structure") and (parameters.get("protein_sequences") or parameters.get("input_files")):
            raise AutomationError("Choose one folding input source.")
    if spec.workflow_id == "design.rfantibody":
        from mn_molecular_design_studio.rfantibody_input import validate_request
        try:
            validate_request(parameters, require_files=not _schema_only_paths.get())
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "design.rfpeptides":
        from mn_molecular_design_studio.rfpeptides_input import validate_request
        try:
            validate_request(parameters, require_files=not _schema_only_paths.get())
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "design.rfdiffusion_beta_pairing":
        from mn_molecular_design_studio import rfdiffusion_beta_pairing_input
        try:
            rfdiffusion_beta_pairing_input.prepare_request(
                parameters, require_files=not _schema_only_paths.get()
            )
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id in {"design.proteina_complexa", "design.proteina_complexa_generate"}:
        from mn_molecular_design_studio.complexa_input import prepare_request
        try:
            prepare_request(parameters | ({"generator_only":True,"search_strategy":"single-pass","replicas":1} if spec.workflow_id.endswith("_generate") else {}))
        except (TypeError, ValueError, OSError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "refolding.af2_design_validation":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.design_validation import validate_af2_design_request
        try:
            validate_af2_design_request(parameters)
        except (TypeError, ValueError) as exc:
            raise AutomationError(str(exc)) from exc
    if spec.workflow_id == "refolding.optimized_prediction":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import validate_prediction_request

        try:
            validate_prediction_request(parameters)
        except (TypeError, ValueError) as exc:
            raise AutomationError(str(exc)) from exc


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(child) for key, child in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(child) for child in value]
    return value


def _resolve_target_entries(value: Any, request_file: Path) -> Any:
    if not isinstance(value, list):
        raise AutomationError("target_entries must be a JSON list.")
    normalized = []
    for index, raw_entry in enumerate(value):
        if not isinstance(raw_entry, dict):
            raise AutomationError(f"target_entries[{index}] must be an object.")
        entry = dict(raw_entry)
        if entry.get("target_pdb"):
            entry["target_pdb"] = str(
                _resolve_path(entry["target_pdb"], request_file=request_file, field_name=f"target_entries[{index}].target_pdb")
            )
        normalized.append(entry)
    return normalized


def _validate_refolding_options(
    options: Any,
    request_file: Path,
    *,
    allow_chain_break_options: bool,
) -> None:
    if not isinstance(options, dict):
        raise AutomationError("options must be a JSON object.")
    allowed = _COMMON_REFOLDING_FIELDS | ({"split_chain_breaks", "chain_break_mode"} if allow_chain_break_options else set())
    unknown = sorted(set(options) - allowed)
    if unknown:
        raise AutomationError(f"Unsupported refolding option(s): {', '.join(unknown)}")
    from mn_molecular_design_studio.engine.workflow_engine.workflows.benchmark import run_de_novo_binder_scoring_dataset

    hints = _type_hints(run_de_novo_binder_scoring_dataset)
    signature = inspect.signature(run_de_novo_binder_scoring_dataset)
    for name, value in options.items():
        if name == "split_chain_breaks":
            if type(value) is not bool:
                raise AutomationError("options.split_chain_breaks must be true or false.")
            continue
        if name == "chain_break_mode":
            if value not in {"preserve_original_chain", "split_fragments"}:
                raise AutomationError("options.chain_break_mode must be preserve_original_chain or split_fragments.")
            continue
        parameter = signature.parameters.get(name)
        if parameter is None:
            raise AutomationError(f"Unsupported refolding option: {name}")
        options[name] = _convert(
            value,
            hints.get(name, parameter.annotation),
            request_file=request_file,
            field_name=f"options.{name}",
        )


def _direct_result(spec: WorkflowSpec, function, parameters: dict[str, Any]) -> tuple[Path, dict[str, Any]]:
    if spec.mode == "target_prepare":
        run_dir = Path(function(**parameters))
        return run_dir, {}
    if spec.mode == "enqueue":
        run_dir = Path(function(**parameters))
        return run_dir, {}
    if spec.mode == "target_refolding_evaluation":
        options = dict(parameters.pop("options", {}) or {})
        source_run_dir, evaluation_run_dir = function(**parameters, **options)
        evaluation_run_dir = Path(evaluation_run_dir)
        return evaluation_run_dir, {"target_candidate_set_run_dir": str(source_run_dir)}
    if spec.mode == "sequence_pipeline":
        run_dir = Path(function(**parameters))
        return run_dir, {}
    if spec.mode == "candidate_refolding":
        options = dict(parameters.pop("options", {}) or {})
        run_dir = Path(function(**parameters, **options, task_group="benchmark"))
        return run_dir, {}
    result = function(**parameters)
    if isinstance(result, (tuple, list)):
        paths = [Path(item) for item in result if isinstance(item, (str, Path))]
        if not paths:
            raise AutomationError(f"Workflow {spec.workflow_id} did not return a run directory.")
        return paths[-1], {"related_run_dirs": [str(path) for path in paths[:-1]]}
    if not isinstance(result, (str, Path)):
        raise AutomationError(f"Workflow {spec.workflow_id} did not return a run directory.")
    return Path(result), {}


def _prepare_selection(spec, payload, request_path):
    """Validate every selected input, then prefer the tool's native batch API."""
    from mn_molecular_design_studio.batching import parameter_items, coalesce_prediction_parameters
    try:
        selected = parameter_items(payload.get("parameters", {}))
    except ValueError as exc:
        raise AutomationError(str(exc)) from exc
    if len(selected) > 1 and spec.mode != "queued":
        raise AutomationError(f"{spec.workflow_id} uses native selection inputs; pass its collection in one parameters object.")
    cpu, gpu = _validate_resources(payload.get("resources"))
    normalized = []
    for item in selected:
        parameters, function = _validate_generic_parameters(spec, item, request_path, source_job=payload.get("source_job"))
        parameters = _resolve_source_job(parameters, payload.get("source_job"), function, spec)
        _apply_resource_parameters(spec, function, parameters, cpu_cores=cpu, gpu_device=gpu)
        _validate_workflow_inputs(spec, parameters)
        normalized.append(parameters)
    merged = coalesce_prediction_parameters(spec.workflow_id, normalized)
    if merged is not None:
        _validate_workflow_inputs(spec, merged)
        return merged, function, cpu, gpu
    if len(normalized) == 1:
        return normalized[0], function, cpu, gpu
    signature = inspect.signature(function).parameters
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import _nested_cpu_thread_requests
    default_cpu_parameter = signature.get("cpu_cores") or signature.get("cpu_threads")
    default_cpu = default_cpu_parameter.default if default_cpu_parameter is not None else 4
    cores = [item.get("cpu_cores") or item.get("cpu_threads") or
             max(_nested_cpu_thread_requests(item), default=default_cpu) for item in normalized]
    cpu = cpu or max((int(value) for value in cores if value is not None and value is not inspect.Parameter.empty), default=4)
    if gpu is None:
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import _queue_resource_for
        default = signature["gpu_device"].default if "gpu_device" in signature else ("all" if _queue_resource_for(spec.tool, selected[0]) else "none")
        devices = {str(item.get("gpu_device", default)) for item in normalized}
        if len(devices) > 1:
            raise AutomationError("Selected inputs must share one GPU allocation; set resources.gpu_device for the batch.")
        gpu = devices.pop()
    return normalized, function, cpu, gpu


def submit_workflow_request(request_file: str | Path) -> dict[str, Any]:
    """Validate and submit one allowlisted workflow request."""
    request_path = Path(request_file).expanduser().resolve()
    try:
        payload = json.loads(request_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise AutomationError(f"Could not read workflow request JSON {request_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise AutomationError("Workflow request must be a JSON object.")
    allowed_fields = {"schema_version", "workflow", "parameters", "resources", "source_job", "priority"}
    unknown = sorted(set(payload) - allowed_fields)
    if unknown:
        raise AutomationError(f"Unsupported workflow request field(s): {', '.join(unknown)}")
    if type(payload.get("schema_version")) is not int or payload["schema_version"] != WORKFLOW_SCHEMA_VERSION:
        raise AutomationError(f"schema_version must be {WORKFLOW_SCHEMA_VERSION}.")
    priority = _request_priority(payload)
    workflow_id = payload.get("workflow")
    if not isinstance(workflow_id, str) or workflow_id not in WORKFLOWS:
        raise AutomationError(f"Unknown workflow '{workflow_id}'. Use 'workflow list' to see supported IDs.")
    spec = WORKFLOWS[workflow_id]
    parameters, function, cpu_cores, resource_gpu = _prepare_selection(spec, payload, request_path)

    with submission_priority(priority):
        if isinstance(parameters, list):
            from mn_molecular_design_studio.engine.workflow_engine.workflows.input_batch import run_input_batch
            run_dir = enqueue_workflow_call(
                run_input_batch, task_group=spec.task_group, tool=spec.tool, job_type=spec.job_type,
                kwargs={"workflow_id": workflow_id, "items": parameters, "cpu_cores": cpu_cores,
                        "gpu_device": resource_gpu, "job_name": parameters[0].get("job_name") or "Selected inputs"},
                cpu_cores=cpu_cores, gpu_device=resource_gpu,
            )
            related = {}
        elif spec.mode != "queued":
            run_dir, related = _direct_result(spec, function, parameters)
            if spec.mode in {"target_prepare", "target_refolding_evaluation", "sequence_pipeline", "candidate_refolding", "enqueue"}:
                _apply_custom_resource_request(run_dir, cpu_cores=cpu_cores, gpu_device=resource_gpu)
                spawn_worker_for_run(run_dir)
        else:
            run_dir = enqueue_workflow_call(
                function,
                task_group=spec.task_group,
                tool=spec.tool,
                job_type=spec.job_type,
                kwargs=parameters,
                cpu_cores=cpu_cores,
                gpu_device=resource_gpu,
            )
            related = {}
    return {"workflow": workflow_id, "job": job_reference(run_dir), **related}


def _apply_custom_resource_request(
    run_dir: Path,
    *,
    cpu_cores: int | None,
    gpu_device: str | None,
) -> None:
    if cpu_cores is None and gpu_device is None:
        return
    from mn_molecular_design_studio.engine.workflow_engine.core.gpu import gpu_queue_resource
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json

    metadata_path = run_dir / "metadata.json"
    metadata = read_json(metadata_path)
    if cpu_cores is not None:
        resource_request = metadata.get("resource_request") if isinstance(metadata.get("resource_request"), dict) else {}
        resource_request["cpu_cores"] = cpu_cores
        metadata["resource_request"] = resource_request
    if gpu_device is not None:
        queue_resource = gpu_queue_resource(gpu_device)
        if queue_resource:
            metadata["queue_resource"] = queue_resource
        else:
            metadata.pop("queue_resource", None)
    write_json(metadata_path, metadata)


def validate_workflow_request(request_file: str | Path) -> dict[str, Any]:
    """Validate a workflow request without making a job or starting a worker."""
    request_path = Path(request_file).expanduser().resolve()
    try:
        payload = json.loads(request_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise AutomationError(f"Could not read workflow request JSON {request_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise AutomationError("Workflow request must be a JSON object.")
    if type(payload.get("schema_version")) is not int or payload["schema_version"] != WORKFLOW_SCHEMA_VERSION:
        raise AutomationError(f"schema_version must be {WORKFLOW_SCHEMA_VERSION}.")
    workflow_id = payload.get("workflow")
    if not isinstance(workflow_id, str) or workflow_id not in WORKFLOWS:
        raise AutomationError(f"Unknown workflow '{workflow_id}'. Use 'workflow list' to see supported IDs.")
    spec = WORKFLOWS[workflow_id]
    allowed_fields = {"schema_version", "workflow", "parameters", "resources", "source_job", "priority"}
    unknown = sorted(set(payload) - allowed_fields)
    if unknown:
        raise AutomationError(f"Unsupported workflow request field(s): {', '.join(unknown)}")
    priority = _request_priority(payload)
    parameters, function, cpu_cores, gpu_device = _prepare_selection(spec, payload, request_path)
    return {
        "valid": True,
        "schema_version": WORKFLOW_SCHEMA_VERSION,
        "workflow": workflow_id,
        "priority": priority,
        "parameters": _jsonable(parameters),
        "resources": {key: value for key, value in (("cpu_cores", cpu_cores), ("gpu_device", gpu_device)) if value is not None},
    }
