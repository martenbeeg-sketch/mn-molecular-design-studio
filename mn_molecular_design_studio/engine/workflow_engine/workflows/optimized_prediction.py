from __future__ import annotations

import hashlib
import json
import math
import os
import re
import secrets
import shlex
import subprocess
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import yaml

from mn_molecular_design_studio import cofolding_msa
from mn_molecular_design_studio import design_validation, structure_templates

from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
    STAGE_COMPLEX_REFOLDING,
    STAGE_MONOMER_REFOLDING,
    read_candidates,
    write_candidates,
)
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import (
    docker_gpu_args,
    gpu_queue_resource,
    normalize_gpu_device,
)
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job,
    finish_job,
    update_status,
    write_json,
)
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root
from mn_molecular_design_studio.engine.workflow_engine.workflows import target_msa as target_msa_workflow


# Prediction engine profiles share the cofolding task payload and result format.
# Optimized engines use the V7 batch API; ESMFold runs through its native CLI.
# BoltzGen uses its fixed batch3 wrapper; the paused V7 campaign's batch2 wrapper
# repeated its full input set. Runtime-gated modes are narrower than some images'
# declared capability lists.
ESMFOLD2_PARAMETER_PRESETS: dict[str, dict[str, int]] = {
    "Fast": {"num_loops": 3, "num_sampling_steps": 50},
    "Standard": {"num_loops": 10, "num_sampling_steps": 68},
    "Careful": {"num_loops": 20, "num_sampling_steps": 68},
    "High diffusion": {"num_loops": 10, "num_sampling_steps": 200},
    "Design rank": {"num_loops": 3, "num_sampling_steps": 200},
}
ESMFOLD2_PARAMETER_PRESET_CHOICES = (
    ("Fast", "Fast · 3 loops / 50 steps"),
    ("Standard", "Standard · 10 loops / 68 steps"),
    ("Careful", "Careful · 20 loops / 68 steps"),
    ("High diffusion", "High diffusion · 10 loops / 200 steps"),
    ("Design rank", "Design rank · 3 loops / 200 steps"),
    ("Custom", "Custom"),
)
ESMFOLD2_PRESET_CONFIG = {
    "selector_key": "parameter_preset",
    "custom_label": "Custom",
    "values": ESMFOLD2_PARAMETER_PRESETS,
}
ESMFOLD2_PARAMETER_FIELDS = (
    {"key": "parameter_preset", "label": "ESMFold2 parameter preset", "type": "select",
     "default": "Standard", "choices": ESMFOLD2_PARAMETER_PRESET_CHOICES,
     "help": "Choose a starting combination; select Custom to enter your own loop and diffusion-step values."},
    {"key": "num_loops", "label": "Number of loops", "default": 10,
     "min": 1, "max": 64, "help": "ESMFold2 recurrent trunk refinement loops."},
    {"key": "num_sampling_steps", "label": "Sampling steps", "default": 68,
     "min": 1, "max": 1000, "step": 1, "help": "Number of ESMFold2 diffusion sampling steps."},
)

ENGINE_PROFILES: dict[str, dict[str, Any]] = {
    "alphafold3": {
        "label": "AlphaFold3 JAX · official weights",
        "code_version": "AlphaFold3 3.1.14",
        "optimized": True,
        "image": "mn-cofolding-af3-memory:3.1.14-cu13-af3-batch1",
        "api_engine": "alphafold3",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "af3_json",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"),
        "default_recycles": 10,
        "input_summary": "Protein, DNA, RNA, and CCD- or SMILES-defined ligands.",
        "terms": True,
    },
    "af3_torch": {
        "label": "AlphaFold3 PyTorch · OpenFold3 P2 weights",
        "code_version": "xfold AF3 22bdeed",
        "optimized": True,
        "image": "mn-af3-torch-opt:22bdeed-cu130-batch3",
        "api_engine": "af3_torch",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "af3_json",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"),
        "default_recycles": 10,
        "input_summary": "Protein, DNA, RNA, and CCD- or SMILES-defined ligands.",
        "settings": ({"key": "diffusion_steps", "label": "Diffusion steps", "default": 200,
                      "min": 1, "max": 1000, "step": 10,
                      "help": "AF3-Torch denoising steps. More steps increase runtime."},),
    },
    "openfold3_jax": {
        "label": "AlphaFold3 JAX · OpenFold3 P2 weights",
        "code_version": "AlphaFold3 JAX 3.1.14 · OpenFold3 P2",
        "optimized": True,
        "image": "mn-cofolding-af3-memory:3.1.14-cu13-openfold3-batch1",
        "api_engine": "openfold3",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "af3_json",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"),
        "default_recycles": 10,
        "input_summary": "Protein, DNA, RNA, and CCD- or SMILES-defined ligands.",
    },
    "boltz2": {
        "label": "Boltz-2",
        "code_version": "Boltz-2 2.2.1",
        "optimized": True,
        "image": "mn-boltz2-opt:2.2.1-cu130-batch2",
        "api_engine": "boltz2",
        "modes": ("off", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "boltz_yaml",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"),
        "default_recycles": 3,
        "input_summary": "Protein, DNA, RNA, and CCD- or SMILES-defined ligands.",
        "supports_cyclic": True,
        "supports_modifications": True,
        "supports_restraints": True,
        "supports_affinity": True,
        "settings": (
            {"key": "sampling_steps", "label": "Diffusion steps", "default": 200,
             "min": 1, "max": 1000, "step": 10,
             "help": "Boltz-2 sampling steps. More steps increase runtime."},
            {"key": "max_parallel_samples", "label": "Parallel samples", "default": 5,
             "min": 1, "max": 5, "step": 1,
             "help": "Limits how many requested diffusion samples Boltz-2 evaluates concurrently; lower values reduce peak memory."},
            {"key": "step_scale", "type": "float", "label": "Step scale", "default": 1.5,
             "min": 0.1, "max": 3.0, "step": 0.1,
             "help": "Lower values generally increase sample diversity; Boltz-2 recommends values around 1–2."},
            {"key": "output_format", "type": "select", "label": "Output format", "default": "pdb",
             "choices": ("pdb", "mmcif")},
            {"key": "use_msa_server", "type": "bool", "label": "Use MSA search", "default": True,
             "help": "Generate missing protein MSAs with Boltz's configured MSA server. Uploaded A3M files are used for their protein chains."},
            {"key": "use_msa", "type": "bool", "label": "Use multiple sequence alignment", "default": True},
            {"key": "generate_ipsae", "type": "bool", "label": "Generate ipSAE metrics", "default": True},
            {"key": "boltz_version", "type": "select", "label": "Boltz version", "default": "2.2.1",
             "choices": ("2.2.1",)},
            {"key": "max_msa_seqs", "label": "Maximum MSA sequences", "default": 8192,
             "min": 1, "max": 8192, "step": 128,
             "help": "Caps the number of sequences used from each protein MSA."},
            {"key": "seed_count", "label": "Number of seeds", "default": 1,
             "min": 1, "max": 16, "step": 1,
             "help": "Run the prediction for this many consecutive seeds, starting from the seed above."},
            {"key": "use_potentials", "type": "bool", "label": "Use potentials", "default": False,
             "help": "Apply Boltz inference-time potentials. Forced contact or pocket restraints enable potentials automatically."},
            {"key": "method", "type": "select", "label": "Experimental method", "default": "",
             "choices": ("", "md", "x-ray diffraction", "electron microscopy", "solution nmr",
                         "solid-state nmr", "neutron diffraction", "theoretical model", "afdb", "boltz-1", "other"),
             "help": "Optional experimental method conditioning."},
        ),
        "supports_msa_upload": True,
        "fixed_settings": ("Full PAE output is always enabled by the installed Boltz-2 engine and retained for ipSAE.",),
        "legacy_settings": ({"key": "write_full_pae", "type": "bool", "label": "Full PAE", "default": True},),
    },
    "protenix_v2": {
        "label": "Protenix V2",
        "code_version": "Protenix 2.0.0",
        "optimized": True,
        "image": "mn-protenix-opt:2.0.0-cu128-batch2",
        "api_engine": "protenix_v2_opt",
        "modes": ("off", "fast"),
        "gpu_count": {},
        "input_kind": "protenix_json",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"),
        "default_recycles": 10,
        "input_summary": "Protein, DNA, RNA, and CCD- or SMILES-defined ligands.",
        "settings": (
            {"key": "step", "label": "Diffusion steps", "default": 200,
             "min": 1, "max": 1000, "step": 10,
             "help": "Protenix V2 diffusion steps."},
            {"key": "use_msa", "type": "bool", "label": "Use MSA search", "default": False,
             "help": "Engine-level setting; this form cannot set MSA behavior per protein chain or upload per-chain A3M files."},
        ),
    },
    "rf3_foundry": {
        "label": "RosettaFold3",
        "code_version": "Foundry 0.2.1.dev13",
        "optimized": True,
        "image": "mn-rf3-opt:0.2.1-dev13-cu130-batch1",
        "api_engine": "rf3_foundry",
        "modes": ("off", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "rf3_json",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"),
        "default_recycles": 10,
        "input_summary": "Protein, DNA, RNA, and CCD- or SMILES-defined ligands.",
        "settings": ({"key": "num_steps", "label": "Diffusion steps", "default": 50,
                      "min": 1, "max": 1000, "step": 10,
                      "help": "RosettaFold3 sampling steps."},),
    },
    "esmfold": {
        "label": "ESMFold",
        "code_version": "fair-esm 2.0.1 · OpenFold 4b41059",
        "optimized": False,
        "image": "mn-esmfold:2.0.1-cu128-sm89",
        "api_engine": "esmfold",
        "modes": ("off",),
        "gpu_count": {},
        "input_kind": "multimer_fasta",
        "default_recycles": 4,
        "max_samples": 1,
        "supports_seed": False,
        "input_summary": "Protein chains entered in the shared cofolding editor. ESMFold can accept multiple chains, but was trained on monomers; treat complex predictions as exploratory.",
        "settings": (
            {"key": "max_tokens_per_batch", "label": "Maximum tokens per GPU batch", "default": 1024,
             "min": 0, "max": 16384, "step": 128,
             "help": "Lower values reduce peak GPU memory for batches of short inputs. Set 0 to disable batching."},
            {"key": "chunk_size", "type": "select", "label": "Attention chunk size", "default": "auto",
             "choices": (("auto", "Auto (default)"), ("128", "128"), ("64", "64"), ("32", "32")),
             "help": "Smaller chunks reduce attention memory use at the cost of speed. ESM recommends 128, 64, or 32 when memory is limited."},
            {"key": "cpu_offload", "type": "bool", "label": "Enable CPU offload", "default": False,
             "help": "Moves some model parameters to system memory to help with longer sequences; this can make inference slower."},
        ),
    },
    "esmfold2": {
        "label": "ESMFold2 Full",
        "code_version": "ESM 3.3.0",
        "optimized": True,
        "image": "mn-esmfold2-opt:0.1.0-cu130-batch1",
        "api_engine": "esmfold2",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "esmfold2_json",
        "default_recycles": 10,
        "input_summary": "Sequence-only protein chains; chain order assigns IDs but does not assign target/binder roles.",
        "settings_replace_num_recycles": True,
        "parameter_presets": ESMFOLD2_PRESET_CONFIG,
        "settings": ESMFOLD2_PARAMETER_FIELDS,
    },
    "esmfold2_fast": {
        "label": "ESMFold2 Fast",
        "code_version": "ESM 3.3.0",
        "optimized": True,
        "image": "mn-esmfold2-opt:0.1.0-cu130-fast-batch1",
        "api_engine": "esmfold2_fast",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "esmfold2_json",
        "default_recycles": 10,
        "input_summary": "Sequence-only protein chains; chain order assigns IDs but does not assign target/binder roles.",
        "settings_replace_num_recycles": True,
        "parameter_presets": ESMFOLD2_PRESET_CONFIG,
        "settings": ESMFOLD2_PARAMETER_FIELDS,
    },
    "chai1": {
        "label": "Chai-1",
        "code_version": "chai_lab 0.6.1",
        "optimized": True,
        "image": "mn-chai1-opt:0.6.1-cu130-batch1",
        "api_engine": "chai1",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {},
        "input_kind": "chai1_json",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_smiles"),
        "default_recycles": 3,
        "input_summary": "Protein, DNA, RNA, and SMILES-defined ligands; generic CCD ligand codes are not accepted by this FASTA adapter.",
        "settings": ({"key": "low_memory", "type": "bool", "label": "Low-memory mode", "default": True,
                      "help": "Reduces Chai-1 peak GPU memory. The pinned batch adapter fixes diffusion timesteps at 200."},),
        "fixed_settings": ("Diffusion timesteps: 200", "MSA and template servers: disabled"),
    },
    "atlasfold": {
        "label": "AtlasFold-M",
        "code_version": "AtlasFold 1.0.0",
        "optimized": True,
        "image": "mn-atlasfold-opt:1.0.0-cu128-batch2",
        "api_engine": "atlasfold",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {},
        "input_kind": "multimer_fasta",
        "default_recycles": 10,
        "input_summary": "Sequence-only protein chains; chain order assigns IDs but does not assign target/binder roles.",
        "settings": ({"key": "num_steps", "label": "Sampling steps", "default": 200,
                      "min": 1, "max": 1000, "step": 10,
                      "help": "AtlasFold prediction steps."},),
    },
    "openfold3_ob0": {
        "label": "OpenFold3 OpenBind-0",
        "code_version": "OpenFold3 0.5.0 · OpenBind-0",
        "optimized": True,
        "image": "mn-openfold3-native-opt:0.5.0-ob0-cu128-batch2",
        "api_engine": "openfold3_ob0",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "openfold3_ob0_json",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"),
        "default_recycles": 10,
        "input_summary": "Protein, DNA, RNA, and CCD- or SMILES-defined ligands.",
        "settings": (
            {"key": "checkpoint", "type": "select", "label": "Model checkpoint",
             "default": "openbind0_174k",
             "choices": (("openbind0_174k", "OpenBind-0 · openbind-2025-06-30-174k · 174k steps"),),
             "help": "This OpenFold3 0.5.0 runner supports OpenBind-0. OpenFold3 P2 (of3-p2-155k.pt) is a legacy checkpoint for OpenFold3 versions below 0.5.0 and cannot be selected for this runner."},
            {"key": "use_msa_server", "type": "bool", "label": "Use MSA server", "default": False,
             "help": "Engine-level setting; this form cannot set MSA behavior per protein chain or upload per-chain A3M files."},
            {"key": "use_templates", "type": "bool", "label": "Use templates", "default": False,
             "help": "Enable the OpenFold3 template search path when available in the image."},
        ),
    },
    "opendde_native": {
        "label": "OpenDDE",
        "code_version": "OpenDDE 1.1.1",
        "optimized": True,
        "image": "mn-opendde-opt:1.1.1-cu128-batch3",
        "api_engine": "opendde_native",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "opendde_json",
        "default_recycles": 10,
        "input_summary": "Sequence-only protein chains; chain order assigns IDs but does not assign target/binder roles.",
        "settings": (
            {"key": "step", "label": "Diffusion steps", "default": 200,
             "min": 1, "max": 1000, "step": 10,
             "help": "OpenDDE diffusion steps."},
            {"key": "use_msa", "type": "bool", "label": "Use MSA search", "default": False,
             "help": "Engine-level setting; this form cannot set MSA behavior per protein chain or upload per-chain A3M files."},
            {"key": "use_template", "type": "bool", "label": "Use templates", "default": False,
             "help": "Enable OpenDDE's configured template search; this may require a template database or network access."},
        ),
    },
    "protenix_v1": {
        "label": "Protenix V1",
        "code_version": "Protenix 1.1.0",
        "optimized": True,
        "image": "mn-protenix-v1-opt:1.1.0-cu130-batch2",
        "api_engine": "protenix_v1",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "protenix_json",
        "molecule_inputs": ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"),
        "default_recycles": 10,
        "input_summary": "Protein, DNA, RNA, and CCD- or SMILES-defined ligands.",
        "settings": (
            {"key": "step", "label": "Diffusion steps", "default": 200,
             "min": 1, "max": 1000, "step": 10,
             "help": "Protenix V1 diffusion steps."},
            {"key": "use_msa", "type": "bool", "label": "Use MSA search", "default": False,
             "help": "Engine-level setting; this form cannot set MSA behavior per protein chain or upload per-chain A3M files."},
        ),
    },
    "colabfold": {
        "label": "ColabFold",
        "code_version": "ColabFold 1.6.1",
        "optimized": True,
        "image": "mn-colabfold-opt:1.6.1-jax0.5.3-cu124-batch3",
        "api_engine": "colabfold_opt",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {"big": 2},
        "input_kind": "multimer_fasta",
        "default_recycles": 3,
        "input_summary": "Sequence-only protein chains; chain order assigns IDs but does not assign target/binder roles.",
        "settings": ({"key": "msa_mode", "type": "select", "label": "MSA mode", "default": "single_sequence",
                      "choices": ("single_sequence", "mmseqs2_uniref_env", "mmseqs2_uniref"),
                      "help": "Single-sequence is fastest; MMseqs2 modes search their configured sequence databases for the multimer input."},),
    },
    "boltzgen_fold": {
        "label": "BoltzGen Fold",
        "code_version": "BoltzGen 0.3.2",
        "optimized": True,
        "image": "mn-boltzgen-opt:a07fea1-cu130-batch3",
        "api_engine": "boltzgen_fold_opt_msa",
        "modes": ("off", "exact", "fast", "big"),
        "gpu_count": {},
        "input_kind": "boltzgen",
        "default_recycles": 3,
        "input_summary": "Target PDB chain and binder protein sequence; MSAs are reused from reference files or fetched from ColabFold.",
        "settings": (
            {"key": "sampling_steps", "label": "Diffusion steps", "default": 200,
             "min": 1, "max": 1000, "step": 10,
             "help": "BoltzGen Fold sampling steps."},
            {"key": "external_msa_max_seqs", "label": "Maximum MSA sequences", "default": 1024,
             "min": 1, "max": 8192, "step": 128,
             "help": "Maximum rows used per target or binder MSA. The full alignment is retained in the reference cache."},
        ),
    },
}


# Every MSA-capable native adapter receives job-local A3Ms prepared by Studio.
ENGINE_PROFILES["alphafold3"]["settings"] = (
    {"key": "model", "type": "select", "label": "AF3-family parameter set", "default": "alphafold3", "hidden": True,
     "choices": ("alphafold3", "openfold3", "openbind0", "boltz2", "protenix2", "rosettafold3", "intellifold2", "opendde")},
)
# Server flags remain accepted for older saved/API requests, but are hidden in
# the common editor: new jobs resolve the reference cache before ColabFold.
for _engine, _profile in ENGINE_PROFILES.items():
    _profile["supports_templates"] = _engine in structure_templates.TEMPLATE_ENGINES
    _profile["requires_target_template"] = _engine == "boltzgen_fold"
    _profile["supports_msa"] = _engine not in {"atlasfold", "esmfold", "esmfold2_fast"}
    if not _profile["supports_msa"]:
        msa_note = (
            "ESMFold uses the submitted protein sequences directly; MSA search is not available."
            if _engine == "esmfold"
            else "This engine uses single sequences and cannot consume an MSA."
        )
        _profile["fixed_settings"] = (*_profile.get("fixed_settings", ()), msa_note)
        continue
    _profile["legacy_settings"] = (*_profile.get("legacy_settings", ()), *(dict(field) for field in _profile.get("settings", ()) if field["key"] in {"use_msa_server", "msa_mode"}))
    _fields = [dict(field) for field in _profile.get("settings", ()) if field["key"] not in {"use_msa_server", "msa_mode"}]
    _fields = [field for field in _fields if field["key"] != "use_msa"]
    if _engine != "boltzgen_fold":
        _fields.append({"key": "use_msa", "type": "bool", "label": "Use multiple sequence alignment", "default": True,
                        "help": "Reuse an exact matching reference MSA first; fetch missing alignments from ColabFold. Save the alignment with the job and in reference_files/msa_cache."})
    if not any(field["key"] in {"max_msa_seqs", "external_msa_max_seqs"} for field in _fields):
        _fields.append({"key": "max_msa_seqs", "label": "Maximum MSA sequences", "default": 1024,
                        "min": 1, "max": 8192, "step": 128,
                        "help": "Limit the job input alignment; the full alignment remains in the reference cache."})
    _profile["settings"] = tuple(_fields)
    _profile["fixed_settings"] = tuple(value for value in _profile.get("fixed_settings", ()) if "MSA and template servers" not in value)
    if _engine == "esmfold2":
        _profile["input_summary"] = "Protein chains with optional per-chain MSAs."


PredictionEngine = Literal[
    "alphafold3", "af3_torch", "openfold3_jax", "boltz2", "protenix_v2",
    "rf3_foundry", "esmfold", "esmfold2", "esmfold2_fast", "chai1", "atlasfold",
    "openfold3_ob0", "opendde_native", "protenix_v1", "colabfold", "boltzgen_fold",
]
OptimizationMode = Literal["off", "exact", "fast", "big"]
_AA = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")
_DNA = frozenset("ACGTN")
_RNA = frozenset("ACGUN")
_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def _parse_proteins(text: str, *, required: bool = True) -> list[dict[str, str]]:
    rows = [part.strip().upper() for part in re.split(r"[:\n\r]+", str(text or "")) if part.strip()]
    if not rows:
        if required:
            raise ValueError("Enter at least one protein sequence. Separate chains with colons or new lines.")
        return []
    if len(rows) > len(_CHAIN_IDS):
        raise ValueError(f"At most {len(_CHAIN_IDS)} protein chains are supported.")
    parsed = []
    for index, sequence in enumerate(rows):
        sequence = "".join(sequence.split())
        invalid = sorted(set(sequence) - _AA)
        if not sequence or invalid:
            raise ValueError(f"Protein chain {index + 1} contains invalid amino-acid symbols: {''.join(invalid)}")
        parsed.append({"id": _CHAIN_IDS[index], "sequence": sequence})
    return parsed


def _parse_nucleic_acids(text: str, *, molecule: str) -> list[dict[str, str]]:
    alphabet = _DNA if molecule == "DNA" else _RNA
    rows = [part.strip().upper() for part in re.split(r"[:\n\r]+", str(text or "")) if part.strip()]
    parsed = []
    for index, value in enumerate(rows, start=1):
        sequence = "".join(value.split())
        invalid = sorted(set(sequence) - alphabet)
        if not sequence or invalid:
            raise ValueError(
                f"{molecule} chain {index} contains unsupported nucleotide symbols: {''.join(invalid)}."
            )
        parsed.append({"id": "", "sequence": sequence})
    return parsed


def _parse_ligands(ccd_text: str, smiles_text: str) -> tuple[list[str], list[str]]:
    codes = [part.strip() for part in re.split(r"[:\n\r]+", str(ccd_text or "")) if part.strip()]
    for code in codes:
        if not re.fullmatch(r"[A-Za-z0-9]{1,8}", code):
            raise ValueError(f"Invalid CCD code {code!r}; use a 1–8 character component ID.")
    smiles = [line.strip() for line in str(smiles_text or "").splitlines() if line.strip()]
    if any(len(value) > 50_000 for value in smiles):
        raise ValueError("A ligand SMILES string exceeds 50,000 characters.")
    return [code.upper() for code in codes], smiles


def _assign_molecule_chain_ids(
    proteins: list[dict[str, str]],
    dna: list[dict[str, str]],
    rna: list[dict[str, str]],
    ligand_ccd_codes: list[str],
    ligand_smiles: list[str],
    *,
    task_index: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    entity_count = len(proteins) + len(dna) + len(rna) + len(ligand_ccd_codes) + len(ligand_smiles)
    if entity_count > len(_CHAIN_IDS):
        raise ValueError(f"Task {task_index} exceeds the {len(_CHAIN_IDS)} supported chains and ligands.")
    next_id = len(proteins)
    for entity in (*dna, *rna):
        entity["id"] = _CHAIN_IDS[next_id]
        next_id += 1
    ccd_entities = []
    for code in ligand_ccd_codes:
        ccd_entities.append({"id": _CHAIN_IDS[next_id], "code": code})
        next_id += 1
    smiles_entities = []
    for value in ligand_smiles:
        smiles_entities.append({"id": _CHAIN_IDS[next_id], "smiles": value})
        next_id += 1
    return proteins, dna, rna, ccd_entities, smiles_entities


def _boltz_residue_reference(
    chain_id: object,
    value: object,
    *,
    chain_lengths: dict[str, int],
    ligand_ids: set[str],
    label: str,
) -> list[Any]:
    selected_chain = str(chain_id or "").strip()
    if selected_chain not in chain_lengths and selected_chain not in ligand_ids:
        raise ValueError(f"{label} refers to unknown chain {selected_chain!r}.")
    residue_text = str(value or "").strip()
    if selected_chain in ligand_ids:
        if not residue_text or residue_text.isdigit():
            raise ValueError(f"{label} for ligand chain {selected_chain} must use a ligand atom name.")
        return [selected_chain, residue_text]
    try:
        residue_index = int(residue_text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must use a one-based residue index.") from exc
    if not 1 <= residue_index <= chain_lengths[selected_chain]:
        raise ValueError(
            f"{label} residue {residue_index} is outside chain {selected_chain} "
            f"(1–{chain_lengths[selected_chain]})."
        )
    return [selected_chain, residue_index]


def _normalize_boltz_chain_options(
    raw_options: object,
    *,
    proteins: list[dict[str, str]],
    dna: list[dict[str, str]],
    rna: list[dict[str, str]],
    ligand_ccd_codes: list[dict[str, str]],
    ligand_smiles: list[dict[str, str]],
    supported: bool,
    task_index: int,
) -> dict[str, dict[str, Any]]:
    if raw_options in (None, []):
        return {}
    if not isinstance(raw_options, list):
        raise ValueError(f"Task {task_index} molecule_chain_options must be a list.")
    protein_by_id = {row["id"]: row for row in proteins}
    normalized: dict[str, dict[str, Any]] = {}
    for row in raw_options:
        if not isinstance(row, dict):
            raise ValueError(f"Task {task_index} molecule chain options must be objects.")
        chain_id = str(row.get("id") or "").strip()
        if chain_id not in protein_by_id:
            raise ValueError(f"Task {task_index} chain options refer to unknown protein chain {chain_id!r}.")
        if chain_id in normalized:
            raise ValueError(f"Task {task_index} has duplicate settings for protein chain {chain_id}.")
        cyclic = row.get("cyclic", False)
        if type(cyclic) is not bool:
            raise ValueError(f"Task {task_index} cyclic setting for chain {chain_id} must be true or false.")
        modifications_text = str(row.get("modifications_text") or "")
        if not supported and (cyclic or modifications_text.strip()):
            raise ValueError(f"Task {task_index} selected cyclic chains or modifications unsupported by this engine.")
        modifications = []
        seen_positions: set[int] = set()
        for line_number, raw_line in enumerate(modifications_text.splitlines(), start=1):
            line = raw_line.strip()
            if not line:
                continue
            match = re.fullmatch(r"(\d+)\s*[,=:]\s*([A-Za-z0-9]{1,8})", line)
            if not match:
                raise ValueError(
                    f"Task {task_index} chain {chain_id} modification line {line_number} must look like '3, MSE'."
                )
            position = int(match.group(1))
            if not 1 <= position <= len(protein_by_id[chain_id]["sequence"]):
                raise ValueError(
                    f"Task {task_index} chain {chain_id} modification position {position} is outside the sequence."
                )
            if position in seen_positions:
                raise ValueError(f"Task {task_index} chain {chain_id} modifies residue {position} more than once.")
            seen_positions.add(position)
            modifications.append({"position": position, "ccd": match.group(2).upper()})
        msa_text = str(row.get("msa_a3m_text") or "").strip()
        if msa_text:
            lines = [line.strip() for line in msa_text.replace("\r", "").splitlines() if line.strip()]
            if not lines or not lines[0].startswith(">"):
                raise ValueError(f"Task {task_index} custom MSA for chain {chain_id} must be an A3M file with a FASTA header.")
            query = []
            for line in lines:
                if line.startswith(">"):
                    if query:
                        break
                    continue
                query.append(line)
            query_sequence = re.sub(r"[a-z.-]", "", "".join(query)).upper()
            if query_sequence and query_sequence != protein_by_id[chain_id]["sequence"]:
                raise ValueError(f"Task {task_index} custom MSA query does not match protein chain {chain_id}.")
        normalized[chain_id] = {
            "cyclic": cyclic,
            "modifications": modifications,
            "msa_a3m_text": msa_text,
            "use_msa": row.get("use_msa", True),
        }
        if type(normalized[chain_id]["use_msa"]) is not bool:
            raise ValueError(f"Task {task_index} chain {chain_id} use_msa must be true or false.")
    return normalized


def _normalize_boltz_restraints(
    raw_restraints: object,
    *,
    proteins: list[dict[str, str]],
    dna: list[dict[str, str]],
    rna: list[dict[str, str]],
    ligand_ccd_codes: list[dict[str, str]],
    ligand_smiles: list[dict[str, str]],
    supported: bool,
    task_index: int,
) -> list[dict[str, Any]]:
    if raw_restraints in (None, []):
        return []
    if not supported:
        raise ValueError(f"Task {task_index} contains restraints unsupported by this engine.")
    if not isinstance(raw_restraints, list):
        raise ValueError(f"Task {task_index} restraints must be a list.")
    chain_lengths = {
        row["id"]: len(row["sequence"])
        for rows in (proteins, dna, rna) for row in rows
    }
    ligand_ids = {
        row["id"] for rows in (ligand_ccd_codes, ligand_smiles) for row in rows
    }

    def distance(value: object, label: str) -> float:
        try:
            result = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{label} must be a distance from 4 to 20 Å.") from exc
        if not math.isfinite(result) or not 4.0 <= result <= 20.0:
            raise ValueError(f"{label} must be a distance from 4 to 20 Å.")
        return result

    def force_value(value: object, label: str) -> bool:
        if type(value) is not bool:
            raise ValueError(f"{label} force setting must be true or false.")
        return value

    result: list[dict[str, Any]] = []
    for index, restraint in enumerate(raw_restraints, start=1):
        label = f"Task {task_index} restraint {index}"
        if not isinstance(restraint, dict) or len(restraint) != 1:
            raise ValueError(f"{label} must have exactly one type: contact, pocket, or bond.")
        kind, details = next(iter(restraint.items()))
        if not isinstance(details, dict):
            raise ValueError(f"{label} details must be an object.")
        if kind == "contact":
            left = details.get("token1")
            right = details.get("token2")
            if not isinstance(left, list) or len(left) != 2 or not isinstance(right, list) or len(right) != 2:
                raise ValueError(f"{label} contact needs two chain and residue references.")
            result.append({"contact": {
                "token1": _boltz_residue_reference(left[0], left[1], chain_lengths=chain_lengths, ligand_ids=ligand_ids, label=f"{label} token 1"),
                "token2": _boltz_residue_reference(right[0], right[1], chain_lengths=chain_lengths, ligand_ids=ligand_ids, label=f"{label} token 2"),
                "max_distance": distance(details.get("max_distance"), label),
                "force": force_value(details.get("force", False), label),
            }})
        elif kind == "pocket":
            binder = str(details.get("binder") or "").strip()
            if binder not in chain_lengths and binder not in ligand_ids:
                raise ValueError(f"{label} refers to unknown binder chain {binder!r}.")
            raw_contacts = details.get("contacts")
            if raw_contacts is None and "contacts_input" in details:
                raw_contacts = []
                for contact_line in str(details.get("contacts_input") or "").splitlines():
                    line = contact_line.strip()
                    if not line:
                        continue
                    parts = [part.strip() for part in line.split(",", maxsplit=1)]
                    if len(parts) != 2:
                        raise ValueError(f"{label} pocket contacts must look like 'A, 56', one per line.")
                    raw_contacts.append(parts)
            if not isinstance(raw_contacts, list) or not raw_contacts:
                raise ValueError(f"{label} pocket restraint needs at least one contact residue.")
            contacts = []
            for contact_index, contact in enumerate(raw_contacts, start=1):
                if not isinstance(contact, list) or len(contact) != 2:
                    raise ValueError(f"{label} pocket contact {contact_index} needs a chain and residue reference.")
                contacts.append(_boltz_residue_reference(
                    contact[0], contact[1], chain_lengths=chain_lengths, ligand_ids=ligand_ids,
                    label=f"{label} pocket contact {contact_index}",
                ))
            result.append({"pocket": {
                "binder": binder,
                "contacts": contacts,
                "max_distance": distance(details.get("max_distance"), label),
                "force": force_value(details.get("force", False), label),
            }})
        elif kind == "bond":
            endpoints = []
            for endpoint_name in ("atom1", "atom2"):
                endpoint = details.get(endpoint_name)
                if not isinstance(endpoint, list) or len(endpoint) != 3:
                    raise ValueError(f"{label} covalent bond needs chain, residue, and atom for {endpoint_name}.")
                endpoint_chain = str(endpoint[0] or "").strip()
                if endpoint_chain in ligand_ids:
                    if str(endpoint[1] or "").strip() != "1":
                        raise ValueError(f"{label} ligand bond endpoints use residue index 1.")
                    reference = [endpoint_chain, 1]
                else:
                    reference = _boltz_residue_reference(
                        endpoint_chain, endpoint[1], chain_lengths=chain_lengths, ligand_ids=ligand_ids,
                        label=f"{label} {endpoint_name}",
                    )
                atom_name = str(endpoint[2] or "").strip()
                if not atom_name:
                    raise ValueError(f"{label} {endpoint_name} needs an atom name.")
                endpoints.append([reference[0], reference[1], atom_name])
            result.append({"bond": {"atom1": endpoints[0], "atom2": endpoints[1]}})
        else:
            raise ValueError(f"{label} type {kind!r} is not supported.")
    return result


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _af3_json(
    name: str,
    proteins: list[dict[str, str]],
    dna: list[dict[str, str]],
    rna: list[dict[str, str]],
    ligand_ccd_codes: list[dict[str, str]],
    ligand_smiles: list[dict[str, str]],
    seed: int,
) -> dict[str, Any]:
    sequences = [
        {"protein": {
            "id": [chain["id"]],
            "sequence": chain["sequence"],
            "modifications": [],
            "templates": [],
            "unpairedMsa": f">query\n{chain['sequence']}\n",
            "pairedMsa": "",
        }}
        for chain in proteins
    ]
    sequences.extend(
        {molecule: {"id": [chain["id"]], "sequence": chain["sequence"]}}
        for molecule, chains in (("dna", dna), ("rna", rna))
        for chain in chains
    )
    sequences.extend(
        {"ligand": {"id": [ligand["id"]], "ccdCodes": [ligand["code"]]}}
        for ligand in ligand_ccd_codes
    )
    sequences.extend(
        {"ligand": {"id": [ligand["id"]], "smiles": ligand["smiles"]}}
        for ligand in ligand_smiles
    )
    return {
        "name": name,
        "sequences": sequences,
        "modelSeeds": [seed],
        "dialect": "alphafold3",
        "version": 1,
    }


def _write_boltzgen_adapter(
    run_dir: Path,
    *,
    item_id: str,
    target_pdb_text: str,
    target_chain_id: str,
    target_msa_a3m: str,
    binder_sequence: str,
    auto_msa: bool = False,
    msa_depth: int = 1024,
) -> tuple[Path, dict[str, Any]]:
    try:
        import numpy as np
        from Bio.PDB import Atom, Chain, Model, PDBIO, PDBParser, Residue, Structure
        from Bio.SeqUtils import seq1, seq3
    except ImportError as exc:
        raise RuntimeError("BoltzGen Fold setup needs Biopython and NumPy in the Studio environment.") from exc

    normalized_binder = "".join(str(binder_sequence).split()).upper()
    if not normalized_binder or set(normalized_binder) - _AA:
        raise ValueError("BoltzGen Fold needs a valid binder protein sequence.")
    if not str(target_pdb_text or "").strip():
        raise ValueError("BoltzGen Fold requires an uploaded target PDB structure.")
    if not re.fullmatch(r"[A-Za-z0-9]", str(target_chain_id or "")) or target_chain_id == "Z":
        raise ValueError("Choose one target chain ID other than the reserved binder chain Z.")

    raw_root = run_dir / "artifacts" / "raw" / "optimized_prediction"
    target_dir = raw_root / "target_template"
    target_dir.mkdir(parents=True, exist_ok=True)
    uploaded = target_dir / f"{item_id}_uploaded_target.pdb"
    uploaded.write_text(target_pdb_text.rstrip() + "\n", encoding="utf-8")
    parser = PDBParser(QUIET=True)
    try:
        source = next(parser.get_structure("uploaded-target", str(uploaded)).get_models())
    except Exception as exc:
        raise ValueError(f"Could not parse the target PDB: {exc}") from exc
    if target_chain_id not in source:
        raise ValueError(f"Target PDB has no chain {target_chain_id!r}.")

    template_structure = Structure.Structure("boltzgen-target")
    template_model = Model.Model(0)
    template_structure.add(template_model)
    target_copy = Chain.Chain(target_chain_id)
    target_sequence = ""
    for residue in source[target_chain_id]:
        if residue.id[0] != " ":
            continue
        letter = seq1(residue.resname, custom_map={"UNK": "X"}).upper()
        if letter not in _AA:
            continue
        target_copy.add(residue.copy())
        target_sequence += letter
    if not target_sequence or not list(target_copy.get_atoms()):
        raise ValueError(f"Target chain {target_chain_id!r} contains no standard protein residues.")
    template_model.add(target_copy)
    template_path = target_dir / f"{item_id}_target_only.pdb"
    writer = PDBIO()
    writer.set_structure(template_structure)
    writer.save(str(template_path))

    msa_text = str(target_msa_a3m or "").strip()
    if auto_msa:
        prepared = cofolding_msa.prepare_chain(run_dir, target_sequence, target_chain_id,
                                               item_id=item_id, depth=msa_depth, supplied=msa_text)
        msa_text = prepared.read_text()
    msa_lines = [line.strip() for line in msa_text.splitlines() if line.strip()]
    if not msa_lines or not msa_lines[0].startswith(">"):
        raise ValueError("BoltzGen Fold needs a target-chain A3M MSA beginning with a FASTA header.")
    query_lines: list[str] = []
    for line in msa_lines:
        if line.startswith(">"):
            if query_lines:
                break
            continue
        query_lines.append(line)
    msa_query = re.sub(r"[a-z.-]", "", "".join(query_lines)).upper()
    if msa_query.upper() != target_sequence:
        raise ValueError("The first BoltzGen target MSA sequence does not match the selected PDB chain.")

    data_dir = raw_root / "boltzgen_data"
    msa_dir = raw_root / "boltzgen_external_msas"
    spec_dir = raw_root / "boltzgen_specs"
    for directory in (data_dir, msa_dir, spec_dir):
        directory.mkdir(parents=True, exist_ok=True)
    try:
        model = next(parser.get_structure("target-template", str(template_path)).get_models())
        target_atoms = [atom.coord.copy() for atom in model.get_atoms()]
        center = np.asarray(target_atoms, dtype=np.float64).mean(axis=0)
        x_start = float(np.asarray(target_atoms)[:, 0].max() + 30.0)
        design_structure = Structure.Structure("boltzgen-input-adapter")
        design_model = Model.Model(0)
        design_structure.add(design_model)
        design_model.add(model[ target_chain_id ].copy())
        binder_chain = Chain.Chain("Z")
        for residue_index, letter in enumerate(normalized_binder, start=1):
            residue = Residue.Residue((" ", residue_index, " "), seq3(letter).upper(), " ")
            x = x_start + (residue_index - 1) * 3.8
            y, z = float(center[1]), float(center[2])
            atom_rows = [("N", (x - 1.2, y + 0.3, z), "N"), ("CA", (x, y, z), "C"),
                         ("C", (x + 1.25, y, z), "C"), ("O", (x + 1.75, y + 0.95, z), "O")]
            if letter != "G":
                atom_rows.append(("CB", (x, y - 0.8, z + 1.4), "C"))
            for atom_index, (atom_name, coord, element) in enumerate(atom_rows, start=1):
                atom = Atom.Atom(atom_name, np.asarray(coord, dtype=np.float64), 0.0, 1.0, " ",
                                 f"{atom_name:>4s}", (residue_index - 1) * 5 + atom_index, element=element)
                residue.add(atom)
            binder_chain.add(residue)
        design_model.add(binder_chain)
        staged_pdb = data_dir / f"{item_id}.pdb"
        writer = PDBIO()
        writer.set_structure(design_structure)
        writer.save(str(staged_pdb))
    except Exception as exc:
        raise ValueError(f"Could not stage the BoltzGen target and binder sequence: {exc}") from exc

    design_mask = np.asarray(
        [chain.id == "Z" for chain in design_model for residue in chain if residue.id[0] == " "], dtype=bool
    )
    np.savez(data_dir / f"{item_id}.npz", design_mask=design_mask)
    item_msa_dir = msa_dir / item_id
    item_msa_dir.mkdir(parents=True, exist_ok=True)
    (item_msa_dir / f"{target_chain_id}.a3m").write_text(msa_text.rstrip() + "\n", encoding="utf-8")

    spec_path = spec_dir / f"{item_id}.yaml"
    spec_path.write_text(yaml.safe_dump({"entities": [
        {"protein": {"id": "Z", "sequence": normalized_binder}},
        {"file": {"path": f"/work/{template_path.relative_to(run_dir).as_posix()}",
                  "include": [{"chain": {"id": target_chain_id}}], "structure_groups": "all"}},
    ]}, sort_keys=False), encoding="utf-8")
    return spec_path, {
        "target_sequence": target_sequence,
        "target_chain_id": target_chain_id,
        "target_only_template": str(template_path),
        "design_structure": str(staged_pdb),
        "target_msa": str(item_msa_dir / f"{target_chain_id}.a3m"),
    }


def _stage_uploaded_boltz_yaml(
    run_dir: Path,
    *,
    item_id: str,
    yaml_text: str,
    auxiliary_files: object,
    use_msa: bool = True,
) -> tuple[Path, dict[str, Any]]:
    try:
        payload = yaml.safe_load(str(yaml_text or ""))
    except yaml.YAMLError as exc:
        raise ValueError(f"Boltz YAML could not be parsed: {exc}") from exc
    if not isinstance(payload, dict) or not isinstance(payload.get("sequences"), list) or not payload["sequences"]:
        raise ValueError("Boltz YAML must be an object with a non-empty sequences list.")
    if payload.get("version", 1) not in (1, "1"):
        raise ValueError("This Boltz adapter accepts YAML input version 1.")

    raw_root = run_dir / "artifacts" / "raw" / "optimized_prediction"
    upload_dir = raw_root / "uploaded_inputs" / item_id
    upload_dir.mkdir(parents=True, exist_ok=True)
    staged_files: dict[str, str] = {}
    if auxiliary_files is None:
        auxiliary_files = []
    if not isinstance(auxiliary_files, list):
        raise ValueError("Boltz YAML auxiliary files must be a list of uploaded files.")
    for row in auxiliary_files:
        if not isinstance(row, dict):
            raise ValueError("Each Boltz YAML auxiliary upload must include a name and file text.")
        original_name = str(row.get("name") or "").replace("\\", "/")
        file_name = original_name.rsplit("/", 1)[-1]
        if not file_name or file_name in {".", ".."}:
            raise ValueError("Boltz YAML auxiliary uploads need a plain file name.")
        if file_name in staged_files:
            raise ValueError(f"Boltz YAML auxiliary file name {file_name!r} was uploaded more than once.")
        if Path(file_name).suffix.lower() not in {".a3m", ".txt", ".cif", ".mmcif", ".pdb", ".ent"}:
            raise ValueError(f"Boltz YAML auxiliary file {file_name!r} has an unsupported extension.")
        destination = upload_dir / file_name
        destination.write_text(str(row.get("text") or "").rstrip() + "\n", encoding="utf-8")
        staged_files[file_name] = f"/work/{destination.relative_to(run_dir).as_posix()}"

    referenced_files: set[str] = set()

    def rewrite_file_references(value: Any) -> Any:
        if isinstance(value, dict):
            rewritten = {}
            for key, nested in value.items():
                if key in {"msa", "path", "cif", "pdb"} and isinstance(nested, str) and nested.strip():
                    if key == "msa" and not use_msa:
                        rewritten[key] = "empty"
                        continue
                    reference = nested.strip()
                    if key == "msa" and reference.lower() in {"empty", "none"}:
                        rewritten[key] = reference
                        continue
                    file_name = reference.replace("\\", "/").rsplit("/", 1)[-1]
                    if file_name in staged_files:
                        rewritten[key] = staged_files[file_name]
                        referenced_files.add(file_name)
                    else:
                        raise ValueError(
                            f"Boltz YAML refers to {reference!r}, but no uploaded file with that name was provided."
                        )
                else:
                    rewritten[key] = rewrite_file_references(nested)
            return rewritten
        if isinstance(value, list):
            return [rewrite_file_references(nested) for nested in value]
        return value

    payload = rewrite_file_references(payload)
    path = raw_root / "inputs" / f"{item_id}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    forced = any(
        isinstance(row, dict)
        and any(isinstance(details, dict) and details.get("force") is True for details in row.values())
        for row in payload.get("constraints", [])
    ) if isinstance(payload.get("constraints", []), list) else False
    return path, {
        "raw_yaml": True,
        "uploaded_files": sorted(staged_files),
        "referenced_files": sorted(referenced_files),
        "force": forced,
    }


def _apply_boltz_input_options(
    run_dir: Path, input_path: Path, *, item_id: str, templates: object = None,
    template_distance_threshold: object = None, use_msa: bool = True,
    custom_msas: object = None, affinity_binder: str = "",
) -> dict[str, Any]:
    payload = yaml.safe_load(input_path.read_text(encoding="utf-8"))
    if templates is None:
        templates = []
    if not isinstance(templates, list) or len(templates) > 16:
        raise ValueError("Upload at most 16 Boltz template files.")
    threshold = None
    if template_distance_threshold not in (None, ""):
        threshold = float(template_distance_threshold)
        if not math.isfinite(threshold) or threshold <= 0:
            raise ValueError("Template distance threshold must be a positive number of Ångstroms.")
        if not templates and not payload.get("templates"):
            raise ValueError("A template distance threshold requires at least one template.")
    upload_dir = run_dir / "artifacts" / "raw" / "optimized_prediction" / "templates" / item_id
    seen_names = set()
    for template in templates:
        if not isinstance(template, dict):
            raise ValueError("Each Boltz template must include a file name and text.")
        name = str(template.get("name") or "").replace("\\", "/").rsplit("/", 1)[-1]
        extension = Path(name).suffix.lower()
        text = str(template.get("text") or "")
        if extension not in {".pdb", ".ent", ".cif", ".mmcif"} or not text.strip():
            raise ValueError("Boltz templates must be non-empty CIF or PDB files.")
        if name in seen_names:
            raise ValueError(f"Template file name {name!r} was uploaded more than once.")
        if len(text.encode("utf-8")) > 20 * 1024 * 1024:
            raise ValueError("Each Boltz template must be smaller than 20 MB.")
        seen_names.add(name)
        upload_dir.mkdir(parents=True, exist_ok=True)
        destination = upload_dir / name
        destination.write_text(text.rstrip() + "\n", encoding="utf-8")
        payload.setdefault("templates", []).append({
            "pdb" if extension in {".pdb", ".ent"} else "cif": f"/work/{destination.relative_to(run_dir).as_posix()}",
        })
    if threshold is not None:
        for template in payload.get("templates", []):
            template.update({"force": True, "threshold": threshold})
    proteins = [row["protein"] for row in payload["sequences"] if isinstance(row, dict) and isinstance(row.get("protein"), dict)]
    if custom_msas is None:
        custom_msas = []
    if not isinstance(custom_msas, list):
        raise ValueError("Custom A3M files must be a list.")
    used_msa_chains = set()
    for msa_file in custom_msas if use_msa else []:
        if not isinstance(msa_file, dict):
            raise ValueError("Each custom A3M must include its protein chain, file name, and text.")
        chain_id = str(msa_file.get("chain_id") or "").strip()
        if not re.fullmatch(r"[A-Za-z0-9_-]{1,32}", chain_id):
            raise ValueError("Custom A3M chain IDs must use 1–32 letters, digits, underscores, or hyphens.")
        text = str(msa_file.get("text") or "").strip()
        matching = [protein for protein in proteins if chain_id in (
            protein.get("id") if isinstance(protein.get("id"), list) else [protein.get("id")]
        )]
        if not chain_id or not matching or chain_id in used_msa_chains:
            raise ValueError("Assign each custom A3M to one distinct protein chain in the input.")
        if not text.startswith(">"):
            raise ValueError(f"Custom A3M for chain {chain_id} must begin with a FASTA header.")
        query = re.split(r"\n>", text, maxsplit=1)[0].splitlines()[1:]
        query_sequence = re.sub(r"[a-z.\-\s]", "", "".join(query)).upper()
        if any(query_sequence != str(protein.get("sequence") or "").upper() for protein in matching):
            raise ValueError(f"The custom A3M query does not match protein chain {chain_id}.")
        used_msa_chains.add(chain_id)
        msa_dir = run_dir / "artifacts" / "raw" / "optimized_prediction" / "msas"
        msa_dir.mkdir(parents=True, exist_ok=True)
        msa_path = msa_dir / f"{item_id}_custom_{chain_id}.a3m"
        msa_path.write_text(text.rstrip() + "\n", encoding="utf-8")
        for protein in matching:
            protein["msa"] = f"/work/{msa_path.relative_to(run_dir).as_posix()}"
    if not use_msa:
        for protein in proteins:
            protein["msa"] = "empty"
    if affinity_binder:
        ligand_ids = {
            chain_id for row in payload["sequences"] if isinstance(row, dict) and isinstance(row.get("ligand"), dict)
            for chain_id in (row["ligand"]["id"] if isinstance(row["ligand"].get("id"), list) else [row["ligand"].get("id")])
        }
        if affinity_binder not in ligand_ids:
            raise ValueError("Affinity prediction requires a ligand chain as the binder.")
        properties = payload.setdefault("properties", [])
        if not any(isinstance(row, dict) and isinstance(row.get("affinity"), dict)
                   and row["affinity"].get("binder") == affinity_binder for row in properties):
            properties.append({"affinity": {"binder": affinity_binder}})
    input_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return {
        "template_count": len(payload.get("templates", [])), "use_msa": use_msa,
        "affinity_binder": affinity_binder,
        "force": any(bool(row.get("force")) for row in payload.get("templates", [])),
    }


def _find_reference_boltz_msa(sequence: str, reference_directory: Path) -> tuple[Path | None, str]:
    """Find an exact-sequence A3M in the configured reference-file caches."""
    cleaned = "".join(str(sequence or "").split()).upper()
    if not cleaned:
        return None, "empty sequence"
    repositories = (
        reference_directory / "msa_cache",
        reference_directory / "boltz_models" / "msa_repository",
    )
    seen: set[Path] = set()
    for repository in repositories:
        try:
            resolved = repository.resolve()
        except OSError:
            resolved = repository
        if resolved in seen:
            continue
        seen.add(resolved)
        cached, source = target_msa_workflow.find_cached_msa_for_sequence(
            cleaned, msa_repository_dir=repository,
        )
        if cached is None or not cached.is_file():
            continue
        # The reference search also repairs older cache files and can crop a
        # longer prefixed query. Confirm the staged query itself is exact.
        try:
            content = cached.read_text(encoding="utf-8").rstrip("\x00").replace("\r\n", "\n")
        except (OSError, UnicodeDecodeError):
            continue
        query_lines: list[str] = []
        in_query = False
        for raw_line in content.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_query:
                    break
                in_query = True
            elif in_query:
                query_lines.extend(char for char in line if char.isupper())
        if "".join(query_lines).replace("-", "") == cleaned:
            return cached, source
    return None, "missing"


def _attach_reference_boltz_msas(
    run_dir: Path,
    input_path: Path,
    *,
    item_id: str,
    use_msa: bool,
) -> dict[str, Any]:
    """Attach reference A3Ms, leaving cache misses for Boltz's ColabFold fallback."""
    if not use_msa:
        return {"reference_msa_count": 0, "reference_msa_chains": [], "reference_msa_missing_chains": []}
    try:
        payload = yaml.safe_load(input_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, yaml.YAMLError) as exc:
        raise ValueError(f"Could not read staged Boltz input while looking up reference MSAs: {exc}") from exc
    if not isinstance(payload, dict) or not isinstance(payload.get("sequences"), list):
        raise ValueError("Staged Boltz input must contain a sequences list before reference MSA lookup.")

    reference_directory = Path(reference_root()).expanduser().resolve()
    msa_dir = run_dir / "artifacts" / "raw" / "optimized_prediction" / "msas"
    attached: list[dict[str, str]] = []
    missing: list[str] = []
    changed = False
    reference_path_by_sequence: dict[str, str] = {}
    protein_rows = [
        row.get("protein") for row in payload["sequences"]
        if isinstance(row, dict) and isinstance(row.get("protein"), dict)
    ]
    msa_by_sequence: dict[str, str] = {}
    for protein in protein_rows:
        sequence = "".join(str(protein.get("sequence") or "").split()).upper()
        current_msa = str(protein.get("msa") or "").strip()
        if not sequence or not current_msa or current_msa.lower() in {"empty", "none"}:
            continue
        previous = msa_by_sequence.get(sequence)
        if previous is not None and previous != current_msa:
            raise ValueError(
                "Boltz requires identical protein sequences to share one MSA file; "
                f"the input assigns different MSAs to sequence {sequence[:12]}…"
            )
        msa_by_sequence[sequence] = current_msa

    for index, row in enumerate(payload["sequences"], start=1):
        protein = row.get("protein") if isinstance(row, dict) else None
        if not isinstance(protein, dict):
            continue
        protein_ids = protein.get("id")
        chain_ids = [str(value) for value in protein_ids] if isinstance(protein_ids, list) else [str(protein_ids or "?")]
        chain_label = ",".join(chain_ids)
        current_msa = str(protein.get("msa") or "").strip()
        if current_msa and current_msa.lower() not in {"empty", "none"}:
            continue  # Keep an explicitly uploaded or YAML-provided A3M.
        sequence = "".join(str(protein.get("sequence") or "").split()).upper()
        if sequence in msa_by_sequence:
            protein["msa"] = msa_by_sequence[sequence]
            changed = True
            if sequence in reference_path_by_sequence:
                attached.append({
                    "chain_id": chain_label,
                    "source": "same_sequence_reference_cache",
                    "reference_path": reference_path_by_sequence[sequence],
                })
            continue
        cached, source = _find_reference_boltz_msa(sequence, reference_directory)
        if cached is None:
            # An explicit "empty" sentinel disables Boltz's auto-generation.
            # Remove it when MSA use is enabled so the configured ColabFold
            # server can fill this cache miss during prediction.
            if "msa" in protein:
                protein.pop("msa", None)
                changed = True
            missing.extend(chain_ids)
            continue
        msa_dir.mkdir(parents=True, exist_ok=True)
        sequence_key = hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]
        destination = msa_dir / f"{item_id}_reference_{sequence_key}.a3m"
        destination.write_bytes(cached.read_bytes())
        protein["msa"] = f"/work/{destination.relative_to(run_dir).as_posix()}"
        msa_by_sequence[sequence] = protein["msa"]
        reference_path_by_sequence[sequence] = str(cached)
        changed = True
        attached.append({"chain_id": chain_label, "source": source, "reference_path": str(cached)})

    if changed:
        input_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return {
        "reference_msa_count": len(attached),
        "reference_msa_chains": attached,
        "reference_msa_missing_chains": missing,
        "msa_fallback": "ColabFold server (https://api.colabfold.com)" if missing else "reference files",
    }


def _cache_generated_boltz_msas(run_dir: Path, output_dir: Path) -> dict[str, list[dict[str, str]]]:
    """Save exact-query Boltz-generated A3Ms in the shared reference cache."""
    input_root = run_dir / "artifacts" / "raw" / "optimized_prediction" / "inputs"
    generated_by_item: dict[str, list[dict[str, str]]] = {}
    reference_directory = Path(reference_root()).expanduser().resolve()
    cache_dir = reference_directory / "msa_cache"

    for input_path in sorted(input_root.glob("*.yaml")):
        item_id = input_path.stem
        if not re.fullmatch(r"[A-Za-z0-9_-]{1,100}", item_id):
            continue
        try:
            payload = yaml.safe_load(input_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, yaml.YAMLError):
            continue
        sequence_chains: dict[str, list[str]] = {}
        for row in payload.get("sequences", []) if isinstance(payload, dict) else []:
            protein = row.get("protein") if isinstance(row, dict) else None
            if not isinstance(protein, dict):
                continue
            sequence = "".join(str(protein.get("sequence") or "").split()).upper()
            raw_ids = protein.get("id")
            chain_ids = [str(value) for value in raw_ids] if isinstance(raw_ids, list) else [str(raw_ids or "Protein")]
            current_msa = str(protein.get("msa") or "").strip().casefold()
            # These chains already had an explicit, staged, or reference MSA.
            if current_msa and current_msa not in {"empty", "none"}:
                continue
            if sequence:
                sequence_chains.setdefault(sequence, []).extend(chain_ids)
        if not sequence_chains:
            continue

        generated_dirs = output_dir.glob(f"**/boltz_results_{item_id}/msa")
        candidates_by_sequence: dict[str, list[Path]] = {sequence: [] for sequence in sequence_chains}
        for msa_dir in generated_dirs:
            for msa_path in msa_dir.rglob("*.a3m"):
                valid, _reason = target_msa_workflow.validate_a3m_file(msa_path)
                if not valid:
                    continue
                try:
                    query_sequence = target_msa_workflow._first_a3m_query_sequence(msa_path)
                except (OSError, UnicodeError):
                    continue
                if query_sequence in candidates_by_sequence:
                    candidates_by_sequence[query_sequence].append(msa_path)

        for sequence, paths in candidates_by_sequence.items():
            if not paths:
                continue
            digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
            cache_path = cache_dir / f"{digest}.a3m"
            temporary_path: Path | None = None
            try:
                cache_dir.mkdir(parents=True, exist_ok=True)
                existing, _source = _find_reference_boltz_msa(sequence, reference_directory)
                if existing is not None:
                    continue
                # Prefer the fullest exact-query alignment when Boltz emits several
                # databases for one chain, while retaining the original job files.
                source_path = max(paths, key=lambda path: (path.stat().st_size, path.name))
                temporary_path = cache_dir / f".{digest}.{run_dir.name}.tmp"
                temporary_path.write_bytes(source_path.read_bytes())
                valid, _reason = target_msa_workflow.validate_a3m_file(temporary_path)
                query_sequence = target_msa_workflow._first_a3m_query_sequence(temporary_path) if valid else ""
                if not valid or query_sequence != sequence:
                    try:
                        temporary_path.unlink(missing_ok=True)
                    except OSError:
                        pass
                    continue
                os.replace(temporary_path, cache_path)
            except Exception as exc:
                if temporary_path is not None:
                    try:
                        temporary_path.unlink(missing_ok=True)
                    except OSError:
                        pass
                try:
                    with (run_dir / "stderr.log").open("a", encoding="utf-8") as stderr:
                        stderr.write(f"Could not cache generated Boltz MSA for chain(s) {', '.join(sequence_chains[sequence])}: {exc}\n")
                except OSError:
                    pass
                continue
            generated_by_item.setdefault(item_id, []).extend(
                {"chain_id": chain_id, "file_name": source_path.name, "cache_path": str(cache_path)}
                for chain_id in dict.fromkeys(sequence_chains[sequence])
            )

    if generated_by_item:
        mapping_path = run_dir / "artifacts" / "raw" / "optimized_prediction" / "input_mapping.json"
        try:
            mapping_payload = json.loads(mapping_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            mapping_payload = {}
        items = mapping_payload.get("items") if isinstance(mapping_payload, dict) else None
        for mapping_item in items if isinstance(items, list) else []:
            if not isinstance(mapping_item, dict):
                continue
            item_id = str(mapping_item.get("id") or "")
            cached_rows = generated_by_item.get(item_id)
            if not cached_rows:
                continue
            details = mapping_item.setdefault("input_details", {})
            if isinstance(details, dict):
                details["generated_msa_cached_chains"] = cached_rows
        if isinstance(mapping_payload, dict):
            write_json(mapping_path, mapping_payload)
    return generated_by_item


def _stage_native_input_without_msa(
    run_dir: Path,
    *,
    profile: dict[str, Any],
    item_id: str,
    job_name: str,
    proteins: list[dict[str, str]],
    dna: list[dict[str, str]],
    rna: list[dict[str, str]],
    ligand_ccd_codes: list[dict[str, str]],
    ligand_smiles: list[dict[str, str]],
    seed: int,
    samples: int,
    target_pdb_text: str,
    target_chain_id: str,
    target_msa_a3m: str,
    binder_sequence: str,
    num_recycles: int,
    molecule_chain_options: dict[str, dict[str, Any]] | None = None,
    restraints: list[dict[str, Any]] | None = None,
    boltz_yaml_text: str = "",
    boltz_yaml_files: list[dict[str, str]] | None = None,
    boltz_templates: list[dict[str, str]] | None = None,
    template_distance_threshold: object = None,
    use_msa: bool = True,
    boltz_custom_msas: list[dict[str, str]] | None = None,
    affinity_binder: str = "",
    auto_msa: bool = False,
    msa_depth: int = 1024,
) -> tuple[Path, dict[str, Any]]:
    raw_root = run_dir / "artifacts" / "raw" / "optimized_prediction"
    input_dir = raw_root / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    kind = profile["input_kind"]
    api_engine = profile["api_engine"]
    if kind == "boltzgen":
        path, details = _write_boltzgen_adapter(
            run_dir, item_id=item_id, target_pdb_text=target_pdb_text,
            target_chain_id=target_chain_id, target_msa_a3m=target_msa_a3m,
            binder_sequence=binder_sequence, auto_msa=auto_msa, msa_depth=msa_depth,
        )
        return path, details

    if kind == "af3_json":
        path = input_dir / f"{item_id}.json"
        _write_json(path, _af3_json(job_name, proteins, dna, rna, ligand_ccd_codes, ligand_smiles, seed))
        return path, {}
    if kind == "esmfold2_json":
        path = input_dir / f"{item_id}.json"
        _write_json(path, {
            "id": job_name,
            "sequences": [{"type": "protein", "id": chain["id"], "sequence": chain["sequence"]} for chain in proteins],
            "seeds": [seed],
            "num_diffusion_samples": samples,
        })
        return path, {}
    if kind == "chai1_json":
        fasta_path = input_dir / f"{item_id}.fasta"
        records = [part for chain in proteins for part in (f">protein|name={chain['id']}", chain["sequence"])]
        records.extend(part for chain in dna for part in (f">dna|name={chain['id']}", chain["sequence"]))
        records.extend(part for chain in rna for part in (f">rna|name={chain['id']}", chain["sequence"]))
        records.extend(part for ligand in ligand_smiles for part in (f">ligand|name={ligand['id']}", ligand["smiles"]))
        fasta_path.write_text("\n".join(records) + "\n", encoding="utf-8")
        path = input_dir / f"{item_id}.json"
        _write_json(path, {"items": [{"id": job_name, "fasta": fasta_path.name, "seeds": [seed]}]})
        return path, {}
    if kind == "multimer_fasta":
        path = input_dir / f"{item_id}.fasta"
        path.write_text(f">{job_name}\n{':'.join(chain['sequence'] for chain in proteins)}\n", encoding="utf-8")
        return path, {}
    if kind in {"protenix_json", "opendde_json"}:
        path = input_dir / f"{item_id}.json"
        chains = []
        for chain in proteins:
            protein_chain: dict[str, Any] = {"sequence": chain["sequence"], "count": 1, "id": [chain["id"]]}
            chains.append({"proteinChain": protein_chain})
        if kind == "protenix_json":
            for molecule, molecule_chains in (("dnaSequence", dna), ("rnaSequence", rna)):
                chains.extend({molecule: {"sequence": chain["sequence"], "count": 1, "id": [chain["id"]]}}
                              for chain in molecule_chains)
            chains.extend({"ligand": {"ligand": f"CCD_{ligand['code']}", "count": 1,
                                       "id": [ligand["id"]]}}
                          for ligand in ligand_ccd_codes)
            chains.extend({"ligand": {"ligand": ligand["smiles"], "count": 1,
                                       "id": [ligand["id"]]}}
                          for ligand in ligand_smiles)
        payload = [{"name": job_name, "modelSeeds": [seed], "sequences": chains}]
        _write_json(path, payload)
        return path, {}
    if kind == "rf3_json":
        path = input_dir / f"{item_id}.json"
        components = [
            {"seq": chain["sequence"], "chain_id": chain["id"]}
            for chains in (proteins, dna, rna) for chain in chains
        ]
        components.extend({"ccd_code": ligand["code"], "chain_id": ligand["id"]}
                          for ligand in ligand_ccd_codes)
        components.extend({"smiles": ligand["smiles"], "chain_id": ligand["id"]}
                          for ligand in ligand_smiles)
        _write_json(path, {"name": job_name, "components": components})
        return path, {}
    if kind == "boltz_yaml":
        if boltz_yaml_text.strip():
            path, details = _stage_uploaded_boltz_yaml(
                run_dir, item_id=item_id, yaml_text=boltz_yaml_text, auxiliary_files=boltz_yaml_files,
                use_msa=use_msa,
            )
            extras = _apply_boltz_input_options(
                run_dir, path, item_id=item_id, templates=boltz_templates,
                template_distance_threshold=template_distance_threshold, use_msa=use_msa,
                custom_msas=boltz_custom_msas, affinity_binder=affinity_binder,
            )
            reference_msas = _attach_reference_boltz_msas(
                run_dir, path, item_id=item_id, use_msa=use_msa,
            )
            return path, {
                **details, **extras, **reference_msas,
                "force": bool(details.get("force") or extras.get("force")),
            }
        path = input_dir / f"{item_id}.yaml"
        chain_options = molecule_chain_options or {}
        sequences = []
        for chain in proteins:
            protein = {"id": chain["id"], "sequence": chain["sequence"]}
            chain_option = chain_options.get(chain["id"], {})
            if chain_option.get("cyclic"):
                protein["cyclic"] = True
            if chain_option.get("modifications"):
                protein["modifications"] = chain_option["modifications"]
            msa_text = str(chain_option.get("msa_a3m_text") or "").strip()
            if msa_text:
                msa_dir = raw_root / "msas"
                msa_dir.mkdir(parents=True, exist_ok=True)
                msa_path = msa_dir / f"{item_id}_chain_{chain['id']}.a3m"
                msa_path.write_text(msa_text.rstrip() + "\n", encoding="utf-8")
                protein["msa"] = f"/work/{msa_path.relative_to(run_dir).as_posix()}"
            sequences.append({"protein": protein})
        sequences.extend({"dna": {"id": chain["id"], "sequence": chain["sequence"]}} for chain in dna)
        sequences.extend({"rna": {"id": chain["id"], "sequence": chain["sequence"]}} for chain in rna)
        sequences.extend({"ligand": {"id": ligand["id"], "ccd": ligand["code"]}}
                         for ligand in ligand_ccd_codes)
        sequences.extend({"ligand": {"id": ligand["id"], "smiles": ligand["smiles"]}}
                         for ligand in ligand_smiles)
        payload: dict[str, Any] = {"version": 1, "sequences": sequences}
        if restraints:
            payload["constraints"] = restraints
        path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
        extras = _apply_boltz_input_options(
            run_dir, path, item_id=item_id, templates=boltz_templates,
            template_distance_threshold=template_distance_threshold, use_msa=use_msa,
            affinity_binder=affinity_binder,
        )
        reference_msas = _attach_reference_boltz_msas(
            run_dir, path, item_id=item_id, use_msa=use_msa,
        )
        return path, {**extras, **reference_msas}
    if kind == "openfold3_ob0_json":
        chain_rows = [
            {"molecule_type": molecule, "chain_ids": chain["id"], "sequence": chain["sequence"]}
            for molecule, molecule_chains in (("protein", proteins), ("dna", dna), ("rna", rna))
            for chain in molecule_chains
        ]
        chain_rows.extend({"molecule_type": "ligand", "chain_ids": ligand["id"], "ccd_codes": [ligand["code"]]}
                          for ligand in ligand_ccd_codes)
        chain_rows.extend({"molecule_type": "ligand", "chain_ids": ligand["id"], "smiles": ligand["smiles"]}
                          for ligand in ligand_smiles)
        path = input_dir / f"{item_id}.json"
        _write_json(path, {"queries": {job_name: {
            "chains": chain_rows, "use_msas": False, "use_main_msas": False, "use_paired_msas": False,
        }}})
        runner_yaml = input_dir / f"{item_id}_runner.yaml"
        runner_yaml.write_text(
            "experiment_settings:\n"
            f"  seeds: [{seed}]\n"
            "  use_msa_server: false\n"
            "  use_templates: false\n"
            "model_update:\n"
            "  presets: [predict]\n"
            "  custom:\n"
            "    architecture:\n"
            "      shared:\n"
            f"        num_recycles: {int(num_recycles)}\n",
            encoding="utf-8",
        )
        return path, {"runner_yaml": str(runner_yaml)}
    raise ValueError(f"No input adapter is registered for {api_engine!r}.")



def _stage_native_input(run_dir: Path, *, prepare_msas: bool = False, msa_depth: int = 1024,
                        allow_msa_server: bool = True, structure_uploads: list[dict] | None = None,
                        binder_use_msa: bool = True, **kwargs: Any) -> tuple[Path, dict[str, Any]]:
    profile = kwargs["profile"]
    enabled = bool(kwargs.get("use_msa", True) and profile.get("supports_msa"))
    kind = profile["input_kind"]
    path, details = _stage_native_input_without_msa(
        run_dir, **kwargs, auto_msa=prepare_msas and enabled, msa_depth=msa_depth)
    details["msa_source"] = "reference_then_colabfold" if enabled else "single_sequence"
    if structure_uploads:
        details.update(structure_templates.stage_templates(
            run_dir, path, item_id=kwargs["item_id"], profile=profile,
            proteins=kwargs["proteins"], uploads=structure_uploads))
    if kind == "boltzgen" and not binder_use_msa:
        destination = Path(details["target_msa"]).parent / 'Z.a3m'
        destination.write_text(f'>query\n{kwargs["binder_sequence"]}\n')
        details["query_only_chains"] = ["Z"]
    if not prepare_msas or not enabled:
        return path, details
    item_id = kwargs["item_id"]
    proteins = kwargs["proteins"]
    if kind == "boltzgen":
        destination = Path(details["target_msa"]).parent / 'Z.a3m'
        if binder_use_msa:
            binder = cofolding_msa.prepare_chain(run_dir, kwargs["binder_sequence"], "Z", item_id=item_id, depth=msa_depth)
            destination.write_text(binder.read_text())
        else:
            destination.write_text(f'>query\n{kwargs["binder_sequence"]}\n')
            details["query_only_chains"] = ["Z"]
        return path, details
    payload = yaml.safe_load(path.read_text()) if kind == "boltz_yaml" else (
        json.loads(path.read_text()) if path.suffix == '.json' else None)
    chain_options = kwargs.get("molecule_chain_options") or {}
    supplied = {chain_id: option.get("msa_a3m_text", "") for chain_id, option in chain_options.items()}
    if kind == "boltz_yaml":
        proteins = []
        reference_ids = {chain for row in details.get("reference_msa_chains", []) for chain in row["chain_id"].split(',')}
        for row in payload["sequences"]:
            protein = row.get("protein")
            if not protein:
                continue
            ids = protein["id"] if isinstance(protein["id"], list) else [protein["id"]]
            for chain_id in ids:
                proteins.append({"id": chain_id, "sequence": protein["sequence"]})
                msa = str(protein.get("msa") or '')
                if msa.startswith('/work/') and chain_id not in reference_ids:
                    supplied[chain_id] = (run_dir / msa.removeprefix('/work/')).read_text()
    msas = {}
    for chain_index, chain in enumerate(proteins):
        staged_name = (
            f"chain_{chain_index}/colabfold_main.a3m"
            if kind == "openfold3_ob0_json" else None
        )
        if not chain_options.get(chain["id"], {}).get("use_msa", True):
            single = run_dir / "artifacts/raw/cofolding_msa" / item_id / (staged_name or f"chain_{chain['id']}.a3m")
            single.parent.mkdir(parents=True, exist_ok=True)
            single.write_text(f'>query\n{chain["sequence"]}\n')
            msas[chain["id"]] = single
            details.setdefault("query_only_chains", []).append(chain["id"])
            continue
        msas[chain["id"]] = cofolding_msa.prepare_chain(
            run_dir, chain["sequence"], chain["id"], item_id=item_id, depth=msa_depth,
            supplied=supplied.get(chain["id"], ''), allow_server=allow_msa_server,
            staged_name=staged_name,
        )
    native = {key: f'/work/{value.relative_to(run_dir).as_posix()}' for key, value in msas.items()}
    if kind == "af3_json":
        for row in payload["sequences"]:
            if "protein" in row:
                protein = row["protein"]
                ids = protein["id"] if isinstance(protein["id"], list) else [protein["id"]]
                protein["unpairedMsa"] = msas[ids[0]].read_text()
                protein["pairedMsa"] = ''
    elif kind == "esmfold2_json":
        for chain in payload["sequences"]:
            chain["msa"] = native[chain["id"]]
    elif kind == "chai1_json":
        payload["msa_dir"] = str(Path(next(iter(native.values()))).parent) if native else ''
    elif kind == "multimer_fasta":
        path = path.with_suffix('.a3m')
        path.write_text(cofolding_msa.complex_a3m(proteins, msas))
        payload = None
    elif kind in {"protenix_json", "opendde_json"}:
        for row in payload[0]["sequences"]:
            chain = row.get("proteinChain")
            if chain:
                key = chain["id"][0]
                chain["unpairedMsaPath"] = native[key]
                paired = msas[key].with_name(msas[key].stem + '_query_only_paired.a3m')
                paired.write_text(f'>query\n{chain["sequence"]}\n')
                chain["pairedMsaPath"] = f'/work/{paired.relative_to(run_dir).as_posix()}'
    elif kind == "rf3_json":
        payload["msa_paths"] = native
    elif kind == "openfold3_ob0_json":
        for query in payload["queries"].values():
            query.update(use_msas=bool(native), use_main_msas=bool(native), use_paired_msas=False)
            for chain in query["chains"]:
                if chain["molecule_type"] == 'protein':
                    chain["main_msa_file_paths"] = native[chain["chain_ids"]]
    elif kind == "boltz_yaml":
        by_sequence = {}
        for row in payload["sequences"]:
            protein = row.get('protein')
            if protein:
                key = protein['id'][0] if isinstance(protein['id'], list) else protein['id']
                # Boltz requires identical sequences to share one MSA path.
                protein['msa'] = by_sequence.setdefault(protein['sequence'], native[key])
        details["reference_msa_missing_chains"] = []
    if payload is not None:
        path.write_text(yaml.safe_dump(payload, sort_keys=False) if kind == 'boltz_yaml' else json.dumps(payload, indent=2) + '\n')
    details["msa_files_by_chain"] = {key: value.relative_to(run_dir).as_posix() for key, value in msas.items()}
    details["msa_manifest"] = 'artifacts/raw/cofolding_msa/manifest.json'
    return path, details


def _batch_options(
    profile: dict[str, Any], *, seed: int, num_recycles: int,
    engine_options: dict[str, Any] | None = None,
    use_potentials: bool = False,
) -> dict[str, Any]:
    engine = profile["api_engine"]
    recycle_count = int(num_recycles)
    settings = engine_options or {}
    use_msa = bool(settings.get("use_msa", True) and profile.get("supports_msa"))
    msa_depth = int(settings.get("max_msa_seqs", 1024)) if use_msa else 1
    if engine == "alphafold3":
        return {
            "model": settings.get("model", "alphafold3"),
            **({"model_dir": "/reference/alphafold3"} if settings.get("model", "alphafold3") == "alphafold3" else {}),
            "num_recycles": recycle_count,
            "flash_attention_implementation": "xla",
        }
    if engine == "openfold3":
        return {"model": "openfold3", "num_recycles": recycle_count, "flash_attention_implementation": "xla"}
    if engine == "af3_torch":
        return {"num_recycles": recycle_count, "diffusion_steps": int(settings.get("diffusion_steps", 200))}
    if engine in {"esmfold2", "esmfold2_fast"}:
        return {"variant": "fast" if engine == "esmfold2_fast" else "full_msa" if use_msa else "full_nomsa",
                "num_loops": int(settings.get("num_loops", recycle_count)),
                "num_sampling_steps": int(settings.get("num_sampling_steps", 200)), "msa_max_depth": msa_depth}
    if engine == "chai1":
        return {"num_trunk_recycles": recycle_count,
                "msa_depth": msa_depth, "low_memory": bool(settings.get("low_memory", True))}
    if engine == "atlasfold":
        return {"num_recycles": recycle_count, "num_steps": int(settings.get("num_steps", 200)), "seed": seed}
    if engine == "colabfold_opt":
        return {"model_type": "alphafold2_multimer_v3",
                "msa_mode": "mmseqs2_uniref_env" if use_msa else "single_sequence", "num_recycle": recycle_count}
    if engine == "openfold3_ob0":
        checkpoint_paths = {
            "openbind0_174k": "/reference/openfold3/of3-ob-2025-06-30-174k.pt",
        }
        selected_checkpoint = str(settings.get("checkpoint", "openbind0_174k"))
        if selected_checkpoint not in checkpoint_paths:
            raise ValueError(f"Unsupported OpenFold3 OpenBind-0 checkpoint: {selected_checkpoint!r}.")
        return {"checkpoint": checkpoint_paths[selected_checkpoint],
                "num_recycles": recycle_count,
                "use_msa_server": False,
                "use_templates": bool(settings.get("use_templates", False))}
    if engine == "opendde_native":
        return {"cycle": recycle_count, "step": int(settings.get("step", 200)), "dtype": "bf16", "model_name": "opendde_v1",
                "use_msa": use_msa,
                "use_template": bool(settings.get("use_template", False))}
    if engine == "protenix_v1":
        return {"cycle": recycle_count, "step": int(settings.get("step", 200)), "dtype": "bf16", "model_name": "protenix_base_default_v1.0.0",
                "use_msa": use_msa,
                "use_template": bool(settings.get("use_template", False))}
    if engine == "protenix_v2_opt":
        return {"model_name": "protenix-v2", "cycle": recycle_count,
                "step": int(settings.get("step", 200)), "dtype": "bf16",
                "use_msa": use_msa,
                "use_template": bool(settings.get("use_template", False)), "use_default_params": True}
    if engine == "boltz2":
        if settings.get("write_full_pae", True) is not True:
            raise ValueError("Full PAE is always emitted by the installed Boltz-2 engine; write_full_pae=False is unsupported.")
        options = {
            # The optimized Boltz image defaults to /models, so mount the
            # shared Boltz cache there instead of relying on its empty image
            # directory. This also gives the runner a writable location for
            # its persisted JIT cache.
            "cache": "/models",
            "recycling_steps": recycle_count,
            "sampling_steps": int(settings.get("sampling_steps", 200)),
            "write_full_pae": True,
            "output_format": str(settings.get("output_format", "mmcif")),
            "max_parallel_samples": int(settings.get("max_parallel_samples", 5)),
            "step_scale": float(settings.get("step_scale", 1.5)),
            "use_potentials": bool(use_potentials or settings.get("use_potentials", False)),
            "use_msa_server": bool(settings.get("use_msa", True) and settings.get("use_msa_server", True)),
            "max_msa_seqs": int(settings.get("max_msa_seqs", 8192)),
            "seed_count": int(settings.get("seed_count", 1)),
        }
        method = str(settings.get("method", "")).strip()
        if method:
            options["method"] = method
        return options
    if engine == "rf3_foundry":
        return {"checkpoint": "/reference/foundry/rf3_foundry_01_24_latest_remapped.ckpt",
                "n_recycles": recycle_count, "num_steps": int(settings.get("num_steps", 50))}
    if engine == "boltzgen_fold_opt_msa":
        raw_root = Path("/work/artifacts/raw/optimized_prediction")
        return {
            "data_dir": str(raw_root / "boltzgen_data"),
            "external_msa_dir": str(raw_root / "boltzgen_external_msas"),
            "external_msa_max_seqs": int(settings.get("external_msa_max_seqs", 1024)),
            "moldir": "/cache/datasets--boltzgen--inference-data/snapshots/c3d36fd276e9caf098c75d4113c6d5eb320b1a4c/mols.zip",
            "checkpoint": "/cache/models--boltzgen--boltzgen-1/snapshots/c1be29e1f82ffcc72264f64b993c43fb4e0d17f0/boltz2_conf_final.ckpt",
            "recycling_steps": recycle_count, "sampling_steps": int(settings.get("sampling_steps", 200)),
        }
    return {}


def _mean_plddt(path: Path) -> float | None:
    try:
        from Bio.PDB import MMCIFParser, PDBParser

        parser = PDBParser(QUIET=True) if path.suffix.lower() == ".pdb" else MMCIFParser(QUIET=True)
        structure = parser.get_structure("optimized_prediction", str(path))
        values = [float(atom.get_bfactor()) for atom in structure.get_atoms()
                  if atom.get_name().strip() == "CA" and float(atom.get_bfactor()) >= 0]
        return round(sum(values) / len(values), 3) if values else None
    except Exception:
        return None


def _confidence_sidecar(path: Path) -> dict[str, Any]:
    from mn_molecular_design_studio.cofolding_results import confidence_files, read_payload, summary_metrics, prediction_row, sample_payload

    scores: dict[str, Any] = prediction_row(path)
    for sidecar in confidence_files(path):
        if any(token in sidecar.name for token in (
            "summary_confidence", "confidences_aggregated", "confidence_", "_confidence.json",
            "_scores_", "scores.model_idx_", "affinity_",
        )) or 'fold_out_npz' in sidecar.parts:
            scores.update(summary_metrics(sample_payload(path, sidecar, read_payload(sidecar))))
    return scores


def _metric_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return round(result, 5) if result == result else None


def _add_boltz_interface_metrics(run_dir: Path, candidates: list[dict[str, Any]], workers: int) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _run_ipsae

    pair_jobs = []
    for candidate in candidates:
        raw = candidate.get("raw_metadata") or {}
        chains = list(dict.fromkeys(row["id"] for row in raw.get("protein_chains", [])))
        if len(chains) < 2:
            continue
        structure = Path(candidate["complex_pdb"])
        native_structure = Path(raw.get("native_prediction_structure") or structure)
        pae = next((path for path in (
            native_structure.with_name(f"pae_{native_structure.stem}.npz"),
            native_structure.with_name(f"{native_structure.stem}_pae.npz"),
        ) if path.is_file()), None)
        if pae is None:
            candidate["metrics"]["ipsae_error"] = "The matching pairwise PAE file is missing."
            continue
        for i, first in enumerate(chains):
            for second in chains[i + 1:]:
                pair_jobs.append({
                    **candidate, "candidate_id": f"{candidate['candidate_id']}__{first}_{second}",
                    "target_chains": [first], "binder_chains": [second], "metrics": {},
                    "raw_metadata": {**raw, "pae_path": str(pae),
                                     "ipsae_parent_id": candidate["candidate_id"],
                                     "ipsae_chain_pair": [first, second]},
                })
    if not pair_jobs:
        return
    try:
        scored = _run_ipsae(run_dir, run_dir, pair_jobs, workers=workers)
    except Exception as exc:
        affected = {row["raw_metadata"]["ipsae_parent_id"] for row in pair_jobs}
        for candidate in candidates:
            if candidate["candidate_id"] in affected:
                candidate["metrics"]["ipsae_error"] = str(exc)
        return
    for candidate in candidates:
        interfaces = [{"chains": row["raw_metadata"]["ipsae_chain_pair"], "metrics": row.get("metrics", {})}
                      for row in scored if row["raw_metadata"]["ipsae_parent_id"] == candidate["candidate_id"]]
        if not interfaces:
            continue
        candidate["raw_metadata"]["ipsae_interfaces"] = interfaces
        scores = [score for row in interfaces if (score := _metric_float(row["metrics"].get("ipsae"))) is not None]
        if scores:
            candidate["metrics"].update({"ipsae": min(scores), "ipsae_min": min(scores),
                                         "ipsae_max": max(scores), "ipsae_avg": round(sum(scores) / len(scores), 5)})
        errors = [row["metrics"]["ipsae_error"] for row in interfaces if row["metrics"].get("ipsae_error")]
        if errors:
            candidate["metrics"]["ipsae_error"] = "; ".join(dict.fromkeys(errors))


def _profile_gpu_count(profile: dict[str, Any], mode: str) -> int:
    return int(profile.get("gpu_count", {}).get(mode, 1))


def _validate_engine_options(profile: dict[str, Any], values: dict[str, Any] | None) -> dict[str, Any]:
    supplied = values or {}
    if not isinstance(supplied, dict):
        raise ValueError("Engine-specific settings must be an object.")
    definitions = {field["key"]: field for field in (*profile.get("settings", ()), *profile.get("legacy_settings", ()))}
    legacy_keys = {field["key"] for field in profile.get("legacy_settings", ())}
    unknown = sorted(set(supplied) - set(definitions))
    if unknown:
        raise ValueError(f"{profile['label']} does not support these editable settings: {', '.join(unknown)}.")
    normalized: dict[str, Any] = {}
    for key, field in definitions.items():
        if key in legacy_keys and key not in supplied:
            continue
        value = supplied.get(key, field["default"])
        field_type = field.get("type", "int")
        if field_type == "bool":
            if type(value) is not bool:
                raise ValueError(f"{field['label']} must be true or false.")
        elif field_type == "select":
            choices = [choice[0] if isinstance(choice, (tuple, list)) else choice for choice in field.get("choices", ())]
            if value not in choices:
                raise ValueError(f"{field['label']} must be one of {', '.join(str(choice) for choice in choices)}.")
        elif field_type == "float":
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise ValueError(f"{field['label']} must be a finite number.")
            if not float(field["min"]) <= float(value) <= float(field["max"]):
                raise ValueError(f"{field['label']} must be from {field['min']} to {field['max']}.")
            value = float(value)
        elif field_type == "int":
            if type(value) is not int or not int(field["min"]) <= value <= int(field["max"]):
                raise ValueError(
                    f"{field['label']} must be an integer from {field['min']} to {field['max']}."
                )
        else:
            raise ValueError(f"{field['label']} has an unsupported setting type.")
        normalized[key] = value
    if "msa_mode" in supplied:
        normalized["use_msa"] = supplied["msa_mode"] != "single_sequence"
    return normalized


@dataclass(frozen=True)
class _PreparedPrediction:
    profile: dict[str, Any]
    selected_mode: str
    recycle_count: int
    seed: int
    clean_name: str
    selected_options: dict[str, Any]
    normalized_tasks: list[dict[str, Any]]
    requested_gpu_count: int
    gpu_text: str


def _prepare_prediction(
    engine: PredictionEngine,
    protein_sequences: str = "",
    optimization_mode: OptimizationMode = "off",
    job_name: str = "optimized-prediction",
    seed: int = 0,
    samples: int = 3,
    num_recycles: int | None = None,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    accept_alphafold3_terms: bool = False,
    target_pdb_text: str = "",
    target_chain_id: str = "A",
    target_msa_a3m: str = "",
    binder_sequence: str = "",
    dna_sequences: str = "",
    rna_sequences: str = "",
    ligand_ccd_codes: str = "",
    ligand_smiles: str = "",
    compound_id: str = "",
    compound_name: str = "",
    tasks: list[dict[str, Any]] | None = None,
    engine_options: dict[str, Any] | None = None,
    validation_sources: list[dict[str, Any]] | None = None,
    validation_mode: Literal["complex", "binder"] | None = None,
    use_target_template: bool = True,
    target_use_msa: bool = True,
    compound_set_file: Path | None = None,
    sequence_file: Path | None = None,
    *,
    resolve_random_seed: bool = True,
) -> _PreparedPrediction:
    """Normalize and check prediction inputs without starting a job or fetching MSAs."""
    profile = ENGINE_PROFILES[str(engine)]
    selected_mode = str(optimization_mode)
    recycle_count = int(profile.get("default_recycles", 10) if num_recycles is None else num_recycles)
    if selected_mode not in profile["modes"]:
        raise ValueError(f"{profile['label']} supports these V7-validated modes: {', '.join(profile['modes'])}.")
    official_weights = str((engine_options or {}).get("model", "alphafold3")) == "alphafold3"
    if profile.get("terms") and official_weights and not accept_alphafold3_terms:
        raise ValueError("Review and accept AlphaFold3's published weight terms before starting this engine.")
    seed = int(seed)
    if not 0 <= seed <= 2**31 - 1:
        raise ValueError("Seed must be between 0 and 2,147,483,647; 0 means random.")
    if resolve_random_seed and seed == 0 and profile.get("supports_seed", True):
        seed = secrets.randbelow(2**31 - 1) + 1
    max_samples = int(profile.get("max_samples", 5))
    if not 1 <= int(samples) <= max_samples:
        raise ValueError(f"{profile['label']} supports between 1 and {max_samples} sample(s) per input.")
    if profile["api_engine"] == "boltzgen_fold_opt_msa" and selected_mode in {"exact", "fast"} and int(samples) < 2:
        raise ValueError(f"BoltzGen Fold {selected_mode} mode requires at least 2 samples.")
    if not 1 <= recycle_count <= 100:
        raise ValueError("Recycle count must be between 1 and 100.")
    if not 1 <= int(cpu_cores) <= 256:
        raise ValueError("CPU cores must be between 1 and 256.")

    is_boltzgen = profile["input_kind"] == "boltzgen"
    supports_boltz_features = str(engine) == "boltz2"
    clean_name = re.sub(r"[^A-Za-z0-9_-]+", "_", str(job_name or "optimized-prediction")).strip("_")[:80]
    clean_name = clean_name or "optimized_prediction"
    selected_options = _validate_engine_options(profile, engine_options)
    if str(engine) == "boltz2":
        if selected_options.get("write_full_pae", True) is not True:
            raise ValueError("Full PAE is always emitted by the installed Boltz-2 engine and retained for ipSAE; write_full_pae=False is unsupported.")
        selected_options["write_full_pae"] = True
    if str(engine) == "alphafold3":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.af3_jax_cofolding import MODEL_PROFILES
        model_profile = MODEL_PROFILES[selected_options.get("model", "alphafold3")]
        if selected_mode not in model_profile["modes"]:
            raise ValueError("This AF3-family parameter set does not support the selected optimization mode.")
        if not official_weights:
            profile = {**profile, "label": model_profile["label"],
                       "gpu_count": {"big": 2} if model_profile["two_gpu_big"] else {}}
    if validation_sources is not None:
        if tasks is not None or protein_sequences.strip() or binder_sequence.strip() or compound_set_file or sequence_file:
            raise ValueError("Use selected designs or new molecule inputs in one batch, not both.")
        tasks = design_validation.build_tasks(
            validation_sources, mode=validation_mode or "complex", profile=profile,
            target_use_msa=target_use_msa, use_target_template=use_target_template,
        )
        if profile.get("supports_msa"):
            selected_options["use_msa"] = bool(target_use_msa and validation_mode != "binder")
    elif validation_mode is not None:
        raise ValueError("Select design runs before choosing a validation mode.")
    if supports_boltz_features and int(seed) + int(selected_options.get("seed_count", 1)) - 1 > 2**31 - 1:
        raise ValueError("The selected Boltz-2 seed range must not exceed 2,147,483,647.")
    raw_tasks = tasks if tasks is not None else [{
        "name": clean_name,
        "protein_sequences": protein_sequences,
        "target_pdb_text": target_pdb_text,
        "target_chain_id": target_chain_id,
        "target_msa_a3m": target_msa_a3m,
        "binder_sequence": binder_sequence,
        "dna_sequences": dna_sequences,
        "rna_sequences": rna_sequences,
        "ligand_ccd_codes": ligand_ccd_codes,
        "ligand_smiles": ligand_smiles,
        "compound_id": compound_id,
        "compound_name": compound_name,
    }]
    from mn_molecular_design_studio.batching import MAX_PREDICTION_ITEMS, prediction_collection_tasks
    if compound_set_file or sequence_file:
        if tasks is not None:
            raise ValueError("Use a collection or explicit prediction tasks, not both.")
        raw_tasks = prediction_collection_tasks(raw_tasks[0], compound_set_file=compound_set_file, sequence_file=sequence_file)
    if supports_boltz_features and isinstance(raw_tasks, list) and any(task.get("source_compound_id") for task in raw_tasks if isinstance(task, dict)):
        selected_options["use_potentials"] = bool((engine_options or {}).get("use_potentials", True))
    task_limit = MAX_PREDICTION_ITEMS
    if not isinstance(raw_tasks, list) or not 1 <= len(raw_tasks) <= task_limit:
        raise ValueError(f"Submit between 1 and {task_limit} independent tasks per batch.")
    normalized_tasks: list[dict[str, Any]] = []
    supported_molecule_inputs = set(profile.get("molecule_inputs", ()))
    for index, raw_task in enumerate(raw_tasks, start=1):
        if not isinstance(raw_task, dict):
            raise ValueError(f"Task {index} must be an object.")
        if validation_sources is None and raw_task.get("validation_reference"):
            raise ValueError("Choose validation_sources to compare a prediction with an existing design.")
        task_job_name = clean_name if len(raw_tasks) == 1 else f"{clean_name}_{index:02d}"[:80]
        default_input_name = task_job_name if supports_boltz_features else f"task_{index:02d}"
        task_name = re.sub(
            r"[^A-Za-z0-9_-]+", "_", str(raw_task.get("name") or default_input_name),
        ).strip("_")[:60]
        task_name = task_name or default_input_name
        boltz_yaml_text = str(raw_task.get("boltz_yaml_text") or "")
        boltz_yaml_files = raw_task.get("boltz_yaml_files") or []
        raw_boltz_yaml: dict[str, Any] | None = None
        if boltz_yaml_text.strip():
            if not supports_boltz_features:
                raise ValueError(f"Task {index} uploaded Boltz YAML for an engine that does not accept it.")
            try:
                raw_boltz_yaml = yaml.safe_load(boltz_yaml_text)
            except yaml.YAMLError as exc:
                raise ValueError(f"Task {index} Boltz YAML could not be parsed: {exc}") from exc
            if (
                not isinstance(raw_boltz_yaml, dict)
                or not isinstance(raw_boltz_yaml.get("sequences"), list)
                or not raw_boltz_yaml["sequences"]
            ):
                raise ValueError(f"Task {index} Boltz YAML must contain a non-empty sequences list.")
            if raw_boltz_yaml.get("version", 1) not in (1, "1"):
                raise ValueError(f"Task {index} uses an unsupported Boltz YAML version; expected version 1.")
        elif boltz_yaml_files:
            raise ValueError(f"Task {index} has Boltz YAML auxiliary files but no YAML input file.")
        task_binder = "".join(str(raw_task.get("binder_sequence") or "").split()).upper()
        molecule_text = {
            "dna_sequences": str(raw_task.get("dna_sequences") or ""),
            "rna_sequences": str(raw_task.get("rna_sequences") or ""),
            "ligand_ccd_codes": str(raw_task.get("ligand_ccd_codes") or ""),
            "ligand_smiles": str(raw_task.get("ligand_smiles") or ""),
        }
        unsupported_molecules = sorted(
            key for key, value in molecule_text.items() if value.strip() and key not in supported_molecule_inputs
        )
        if unsupported_molecules:
            raise ValueError(
                f"{profile['label']} does not support these molecule inputs: {', '.join(unsupported_molecules)}."
            )
        dna_chains = _parse_nucleic_acids(molecule_text["dna_sequences"], molecule="DNA")
        rna_chains = _parse_nucleic_acids(molecule_text["rna_sequences"], molecule="RNA")
        task_ccd_codes, task_smiles = _parse_ligands(
            molecule_text["ligand_ccd_codes"], molecule_text["ligand_smiles"],
        )
        yaml_has_ligand = bool(raw_boltz_yaml and any(
            isinstance(row, dict) and isinstance(row.get("ligand"), dict)
            for row in raw_boltz_yaml.get("sequences", [])
        ))
        yaml_ligand_smiles = [
            str(row["ligand"].get("smiles") or "").strip()
            for row in (raw_boltz_yaml or {}).get("sequences", [])
            if isinstance(row, dict) and isinstance(row.get("ligand"), dict)
            and str(row["ligand"].get("smiles") or "").strip()
        ]
        has_ligand = bool(task_ccd_codes or task_smiles or yaml_has_ligand)
        task_compound_name = str(raw_task.get("compound_name") or "").strip()
        task_compound_id = str(raw_task.get("compound_id") or raw_task.get("source_compound_id") or "").strip()
        if has_ligand and not task_compound_name:
            raise ValueError(f"Task {index} needs a compound name for its ligand input.")
        if has_ligand and not task_compound_id:
            task_compound_id = task_compound_name
        target_ids: list[str] = []
        binder_ids: list[str] = []
        if is_boltzgen:
            task_proteins = []
        elif str(raw_task.get("target_protein_sequences") or "").strip() or task_binder:
            target_chains = _parse_proteins(str(raw_task.get("target_protein_sequences") or ""), required=False)
            binder_chains = _parse_proteins(task_binder) if task_binder else []
            if target_chains and not binder_chains:
                raise ValueError(f"Task {index} needs one binder/partner protein sequence after its target chains.")
            if len(target_chains) + len(binder_chains) > len(_CHAIN_IDS):
                raise ValueError(f"Task {index} exceeds the {len(_CHAIN_IDS)} supported protein chains.")
            for chain_index, chain in enumerate(binder_chains, start=len(target_chains)):
                chain["id"] = _CHAIN_IDS[chain_index]
            task_proteins = [*target_chains, *binder_chains]
            target_ids = [chain["id"] for chain in target_chains]
            binder_ids = [chain["id"] for chain in binder_chains]
        elif str(raw_task.get("protein_sequences") or "").strip():
            task_proteins = _parse_proteins(str(raw_task.get("protein_sequences") or ""))
        else:
            task_proteins = []
        if is_boltzgen and (not task_binder or set(task_binder) - _AA):
            raise ValueError(f"Task {index} needs a valid BoltzGen Fold binder protein sequence.")
        if not is_boltzgen:
            task_proteins, dna_chains, rna_chains, ccd_entities, smiles_entities = _assign_molecule_chain_ids(
                task_proteins, dna_chains, rna_chains, task_ccd_codes, task_smiles, task_index=index,
            )
            if not (task_proteins or dna_chains or rna_chains or ccd_entities or smiles_entities or raw_boltz_yaml):
                raise ValueError(f"Task {index} needs at least one protein, DNA, RNA, or ligand input.")
            if raw_boltz_yaml and (
                str(raw_task.get("protein_sequences") or "").strip()
                or str(raw_task.get("target_protein_sequences") or "").strip()
                or task_binder
                or any(value.strip() for value in molecule_text.values())
            ):
                raise ValueError(f"Task {index} cannot mix uploaded Boltz YAML with molecule-builder inputs.")
            task_ccd_codes, task_smiles = ccd_entities, smiles_entities
            if binder_ids:
                task_binder = ":".join(chain["sequence"] for chain in task_proteins if chain["id"] in binder_ids)
        else:
            ccd_entities, smiles_entities = [], []
        scientific_contract = raw_task.get("scientific_contract") or {}
        if scientific_contract:
            expected = {scientific_contract["chain_map"][row["chain"]]: row["sequence"]
                        for row in scientific_contract["components"] if row["type"] == "protein"}
            if expected != {row["id"]: row["sequence"] for row in task_proteins}:
                raise ValueError("Scientific contract protein sequence/chain mapping differs from task.")
            target_ids = list(scientific_contract["target_protein_ids"])
            binder_ids = list(scientific_contract["binder_protein_ids"])
            if set(target_ids) & set(binder_ids) or set(target_ids + binder_ids) - expected.keys():
                raise ValueError("Scientific contract has invalid target/binder roles.")
            task_binder = ":".join(row["sequence"] for row in task_proteins if row["id"] in binder_ids)
        molecule_chain_options = _normalize_boltz_chain_options(
            raw_task.get("molecule_chain_options", []),
            proteins=task_proteins, dna=dna_chains, rna=rna_chains,
            ligand_ccd_codes=ccd_entities, ligand_smiles=smiles_entities,
            supported=supports_boltz_features, task_index=index,
        )
        if selected_options.get("use_msa", True) and profile["input_kind"] in {"boltz_yaml", "multimer_fasta", "chai1_json"}:
            msa_policy_by_sequence = {}
            for chain in task_proteins:
                policy = molecule_chain_options.get(chain["id"], {}).get("use_msa", True)
                previous = msa_policy_by_sequence.setdefault(chain["sequence"], policy)
                if previous != policy:
                    raise ValueError(f"Task {index}: this engine shares MSAs between identical sequences. Turn target MSA off when the target and binder sequences are identical.")
        restraints = _normalize_boltz_restraints(
            raw_task.get("restraints", []),
            proteins=task_proteins, dna=dna_chains, rna=rna_chains,
            ligand_ccd_codes=ccd_entities, ligand_smiles=smiles_entities,
            supported=supports_boltz_features, task_index=index,
        )
        raw_yaml_proteins = [
            row.get("protein") for row in raw_boltz_yaml.get("sequences", [])
            if isinstance(row, dict) and isinstance(row.get("protein"), dict)
        ] if raw_boltz_yaml else []
        raw_yaml_force = any(
            isinstance(row, dict)
            and any(isinstance(details, dict) and details.get("force") is True for details in row.values())
            for row in raw_boltz_yaml.get("constraints", [])
        ) if raw_boltz_yaml and isinstance(raw_boltz_yaml.get("constraints", []), list) else False
        if raw_boltz_yaml and isinstance(raw_boltz_yaml.get("templates"), list):
            raw_yaml_force = raw_yaml_force or any(
                isinstance(row, dict) and row.get("force") is True for row in raw_boltz_yaml["templates"]
            )
        item_id = "i" + hashlib.sha256(
            f"{job_name}:{engine}:{index}".encode("utf-8")
        ).hexdigest()[:16]
        normalized_tasks.append({
            "index": index,
            "id": item_id,
            "name": task_name,
            "source_input_name": str(task_compound_name or raw_task.get("name") or default_input_name),
            "compound_id": task_compound_id,
            "compound_name": task_compound_name,
            "compound_smiles": [str(row.get("smiles") or "") for row in task_smiles] + yaml_ligand_smiles,
            "job_name": task_job_name,
            "proteins": task_proteins,
            "scientific_contract": scientific_contract,
            "target_protein_ids": target_ids,
            "binder_protein_ids": binder_ids,
            "dna": dna_chains,
            "rna": rna_chains,
            "ligand_ccd_codes": task_ccd_codes,
            "ligand_smiles": task_smiles,
            "target_pdb_text": str(raw_task.get("target_pdb_text") or ""),
            "target_chain_id": str(raw_task.get("target_chain_id") or "A"),
            "target_msa_a3m": str(raw_task.get("target_msa_a3m") or ""),
            "binder_sequence": task_binder,
            "molecule_chain_options": molecule_chain_options,
            "restraints": restraints,
            "boltz_yaml_text": boltz_yaml_text,
            "boltz_yaml_files": boltz_yaml_files,
            "boltz_templates": raw_task.get("boltz_templates") or [],
            "template_distance_threshold": raw_task.get("template_distance_threshold"),
            "boltz_custom_msas": raw_task.get("boltz_custom_msas") or [],
            "affinity_binder": str(raw_task.get("affinity_binder") or (smiles_entities[0]["id"] if supports_boltz_features and raw_task.get("source_compound_id") and smiles_entities else "")).strip(),
            "structure_templates": structure_templates.validate_uploads(
                raw_task.get("structure_templates"), supported=bool(profile.get("supports_templates"))),
            "validation_reference": raw_task.get("validation_reference"),
            "source_candidate_id": str(raw_task.get("source_candidate_id") or ""),
            "source_compound_id": str(raw_task.get("source_compound_id") or task_compound_id),
            "source_compound_name": task_compound_name,
            "source_run_dir": str(raw_task.get("source_run_dir") or ""),
            "binder_use_msa": bool(raw_task.get("binder_use_msa", True)),
            "raw_yaml_proteins": raw_yaml_proteins,
            "force_potentials": bool(raw_yaml_force or raw_task.get("template_distance_threshold")),
        })
    if supports_boltz_features and selected_options.get("use_msa", True) and not bool(selected_options.get("use_msa_server", True)):
        configured_reference = Path(reference_root()).expanduser().resolve()
        missing_msa_chains = [
            chain["id"]
            for task in normalized_tasks
            for chain in task["proteins"]
            if task["molecule_chain_options"].get(chain["id"], {}).get("use_msa", True)
            if not task["molecule_chain_options"].get(chain["id"], {}).get("msa_a3m_text")
            and _find_reference_boltz_msa(chain["sequence"], configured_reference)[0] is None
        ]
        missing_msa_chains.extend(
            f"Task {task['index']} protein {protein.get('id', '?')}"
            for task in normalized_tasks
            for protein in task["raw_yaml_proteins"]
            if not protein.get("msa") or str(protein.get("msa")).strip().lower() in {"empty", "none"}
            if _find_reference_boltz_msa(
                str(protein.get("sequence") or ""), configured_reference,
            )[0] is None
        )
        if missing_msa_chains:
            raise ValueError(
                "Boltz-2 MSA search is off, but protein chains are missing custom A3M files: "
                + ", ".join(missing_msa_chains)
                + ". Upload an A3M for each protein chain or turn MSA search on."
            )
    requested_gpu_count = _profile_gpu_count(profile, selected_mode)
    if any(task["structure_templates"] for task in normalized_tasks):
        if engine == "openfold3_ob0":
            selected_options["use_templates"] = True
        elif engine in {"protenix_v1", "protenix_v2"}:
            selected_options["use_template"] = True
    gpu_text = normalize_gpu_device(gpu_device)
    gpu_ids = [part.strip() for part in gpu_text.removeprefix("device=").split(",") if part.strip()]
    if any(not part.isdigit() for part in gpu_ids) or len(gpu_ids) < requested_gpu_count:
        raise ValueError(
            f"{profile['label']} {selected_mode} mode requires {requested_gpu_count} reserved GPU(s); "
            "set the GPU allocation before starting."
        )

    return _PreparedPrediction(
        profile=profile, selected_mode=selected_mode, recycle_count=recycle_count,
        seed=seed, clean_name=clean_name, selected_options=selected_options,
        normalized_tasks=normalized_tasks, requested_gpu_count=requested_gpu_count,
        gpu_text=gpu_text,
    )


def validate_prediction_request(parameters: dict[str, Any]) -> None:
    """Apply worker input checks before queueing, without jobs, searches, or seed generation."""
    _prepare_prediction(**parameters, resolve_random_seed=False)


def boltz2_prediction_parameters(
    tasks: list[dict[str, Any]], *, settings: dict[str, Any], job_name: str,
) -> dict[str, Any]:
    """Build the canonical optimized-prediction request used by the Boltz-2 tile."""
    from mn_molecular_design_studio.prediction_configuration import boltz_settings
    settings = boltz_settings(settings)
    mode = str(settings.get("optimization_mode") or "off")
    option_keys = (
        "sampling_steps", "max_parallel_samples", "step_scale", "output_format",
        "seed_count", "use_potentials", "method", "use_msa", "generate_ipsae", "boltz_version", "write_full_pae",
    )
    engine_options = {key: settings[key] for key in option_keys if key in settings}
    engine_options["max_msa_seqs"] = int(settings.get("max_msa_seqs") or 8192)
    # The tile exposes MSA as one user-facing switch; disabling it also disables
    # the server fallback. Callers may pass the derived value explicitly.
    engine_options["use_msa_server"] = bool(settings.get("use_msa_server", settings.get("use_msa", True)))
    return {
        "engine": "boltz2",
        "optimization_mode": mode,
        "job_name": str(job_name or "boltz2-prediction"),
        "seed": int(settings.get("seed") or 0),
        "samples": int(settings.get("samples", 1)),
        "num_recycles": int(settings.get("num_recycles", ENGINE_PROFILES["boltz2"]["default_recycles"])),
        "cpu_cores": int(settings.get("cpu_cores") or 4),
        "gpu_device": str(settings.get("gpu_device") or ("0,1" if mode == "big" else "0")),
        "tasks": [dict(task) for task in tasks],
        "engine_options": engine_options,
    }


def run_optimized_prediction(
    engine: PredictionEngine,
    protein_sequences: str = "",
    optimization_mode: OptimizationMode = "off",
    job_name: str = "optimized-prediction",
    seed: int = 0,
    samples: int = 3,
    num_recycles: int | None = None,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    accept_alphafold3_terms: bool = False,
    target_pdb_text: str = "",
    target_chain_id: str = "A",
    target_msa_a3m: str = "",
    binder_sequence: str = "",
    dna_sequences: str = "",
    rna_sequences: str = "",
    ligand_ccd_codes: str = "",
    ligand_smiles: str = "",
    compound_id: str = "",
    compound_name: str = "",
    tasks: list[dict[str, Any]] | None = None,
    engine_options: dict[str, Any] | None = None,
    validation_sources: list[dict[str, Any]] | None = None,
    validation_mode: Literal["complex", "binder"] | None = None,
    use_target_template: bool = True,
    target_use_msa: bool = True,
    compound_set_file: Path | None = None,
    sequence_file: Path | None = None,
) -> Path:
    """Run one engine using the same input checks as API preflight validation."""
    prepared = _prepare_prediction(
        engine=engine, protein_sequences=protein_sequences,
        optimization_mode=optimization_mode, job_name=job_name, seed=seed,
        samples=samples, num_recycles=num_recycles, cpu_cores=cpu_cores,
        gpu_device=gpu_device, accept_alphafold3_terms=accept_alphafold3_terms,
        target_pdb_text=target_pdb_text, target_chain_id=target_chain_id,
        target_msa_a3m=target_msa_a3m, binder_sequence=binder_sequence,
        dna_sequences=dna_sequences, rna_sequences=rna_sequences,
        ligand_ccd_codes=ligand_ccd_codes, ligand_smiles=ligand_smiles,
        compound_id=compound_id, compound_name=compound_name,
        tasks=tasks, engine_options=engine_options,
        validation_sources=validation_sources, validation_mode=validation_mode,
        use_target_template=use_target_template, target_use_msa=target_use_msa,
        compound_set_file=compound_set_file, sequence_file=sequence_file,
    )
    if engine == "boltz2":
        from mn_molecular_design_studio.scientific_contracts import validate_component_chemistry
        for task in prepared.normalized_tasks:
            validate_component_chemistry(task, reference_root())
    profile = prepared.profile
    selected_mode = prepared.selected_mode
    recycle_count = prepared.recycle_count
    seed = prepared.seed
    clean_name = prepared.clean_name
    selected_options = prepared.selected_options
    normalized_tasks = prepared.normalized_tasks
    requested_gpu_count = prepared.requested_gpu_count
    gpu_text = prepared.gpu_text
    is_boltzgen = profile["input_kind"] == "boltzgen"
    supports_boltz_features = str(engine) == "boltz2"

    job = create_job(
        "refolding-validation",
        job_type="monomer_refolding" if validation_mode == "binder" else "complex_refolding",
        tool="optimized_prediction_" + str(engine),
        inputs={
            "tasks": [{
                "name": task["name"],
                "source_input_name": task["source_input_name"],
                "protein_chains": [
                    {"id": row["id"], "sequence": row["sequence"],
                     "role": (
                         "binder" if row["id"] in task["binder_protein_ids"]
                         else "target" if row["id"] in task["target_protein_ids"]
                         else "unspecified"
                     )}
                    for row in task["proteins"]
                ],
                "dna_chains": [{"id": row["id"], "sequence": row["sequence"]} for row in task["dna"]],
                "rna_chains": [{"id": row["id"], "sequence": row["sequence"]} for row in task["rna"]],
                "ligand_ccd_codes": task["ligand_ccd_codes"],
                "ligand_smiles": task["ligand_smiles"],
                "compound_id": task["compound_id"], "compound_name": task["compound_name"],
                "restraints": task["restraints"],
                "target_chain_id": task["target_chain_id"] if is_boltzgen else None,
                "binder_sequence": task["binder_sequence"] if is_boltzgen else None,
                "has_target_structure": bool(task["target_pdb_text"]) if is_boltzgen else False,
            } for task in normalized_tasks],
            "engine": str(engine),
        },
        params={
            "engine": str(engine), "optimization_mode": selected_mode, "image": profile["image"],
            "api_engine": profile["api_engine"], "seed": int(seed), "samples": int(samples),
            "num_recycles": recycle_count, "engine_options": selected_options,
            "task_count": len(normalized_tasks), "cpu_cores": int(cpu_cores),
            "gpu_device": gpu_text, "gpu_count": requested_gpu_count,
            "queue_resource": gpu_queue_resource(gpu_text),
            "accept_alphafold3_terms": bool(accept_alphafold3_terms),
            "validation_mode": validation_mode,
            "validation_sources": validation_sources,
        },
    )
    # My Results displays the user-facing job name from campaign_name, matching
    # the naming convention used by BindCraft and other Studio workflows.
    update_status(
        job.run_dir, "preparing", campaign_name=clean_name,
        current_phase="Preparing prediction inputs",
    )
    raw_root = job.run_dir / "artifacts" / "raw" / "optimized_prediction"
    output_dir = raw_root / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        manifest_items = []
        input_mappings = []
        runner_yamls = {}
        chai_template_hits = []
        for task in normalized_tasks:
            if task["validation_reference"]:
                task["validation_reference"] = design_validation.stage_reference(
                    job.run_dir, task["id"], task["validation_reference"])
            if profile["api_engine"] == "esmfold":
                input_dir = raw_root / "inputs"
                input_dir.mkdir(parents=True, exist_ok=True)
                input_path = input_dir / f"{task['id']}.fasta"
                complex_sequence = ":".join(chain["sequence"] for chain in task["proteins"])
                input_path.write_text(f">{task['id']}\n{complex_sequence}\n", encoding="utf-8")
                input_details = {"msa_source": "single_sequence"}
            else:
                input_path, input_details = _stage_native_input(
                    job.run_dir, profile=profile, item_id=task["id"], job_name=task["job_name"],
                    proteins=task["proteins"], dna=task["dna"], rna=task["rna"],
                    ligand_ccd_codes=task["ligand_ccd_codes"], ligand_smiles=task["ligand_smiles"],
                    seed=int(seed), samples=int(samples),
                    target_pdb_text=task["target_pdb_text"], target_chain_id=task["target_chain_id"],
                    target_msa_a3m=task["target_msa_a3m"], binder_sequence=task["binder_sequence"],
                    num_recycles=recycle_count,
                    molecule_chain_options=task["molecule_chain_options"],
                    restraints=task["restraints"],
                    boltz_yaml_text=task["boltz_yaml_text"],
                    boltz_yaml_files=task["boltz_yaml_files"],
                    boltz_templates=task["boltz_templates"],
                    template_distance_threshold=task["template_distance_threshold"],
                    use_msa=bool(selected_options.get("use_msa", True)),
                    boltz_custom_msas=task["boltz_custom_msas"],
                    affinity_binder=task["affinity_binder"],
                    prepare_msas=True,
                    msa_depth=int(selected_options.get("max_msa_seqs", selected_options.get("external_msa_max_seqs", 1024))),
                    allow_msa_server=bool(selected_options.get("use_msa_server", True)) if engine == "boltz2" else True,
                    structure_uploads=task["structure_templates"], binder_use_msa=task["binder_use_msa"],
                )
            manifest_items.append({
                "id": task["id"],
                "input": f"/work/{input_path.relative_to(job.run_dir).as_posix()}",
                **({"template_path": input_details["batch_template_path"]} if input_details.get("batch_template_path") else {}),
            })
            if input_details.get("chai_template_hits"):
                chai_template_hits.append(input_details["chai_template_hits"])
            input_mappings.append({
                "id": task["id"], "name": task["name"], "source_input_name": task["source_input_name"],
                "compound_id": task["compound_id"], "compound_name": task["compound_name"],
                "protein_chains": [
                    {"id": row["id"], "length": len(row["sequence"]),
                     "role": (
                         "binder" if row["id"] in task["binder_protein_ids"]
                         else "target" if row["id"] in task["target_protein_ids"]
                         else "unspecified"
                     )}
                    for row in task["proteins"]
                ],
                "dna_chains": [{"id": row["id"], "length": len(row["sequence"])} for row in task["dna"]],
                "rna_chains": [{"id": row["id"], "length": len(row["sequence"])} for row in task["rna"]],
                "ligand_ccd_codes": task["ligand_ccd_codes"],
                "ligand_smiles": [row.get("smiles") for row in task["ligand_smiles"]],
                "ligand_smiles_count": len(task["ligand_smiles"]),
                "raw_yaml": bool(task["boltz_yaml_text"]),
                "input_file": input_path.relative_to(job.run_dir).as_posix(),
                "input_details": input_details,
            })
            if input_details.get("runner_yaml"):
                runner_path = Path(str(input_details["runner_yaml"]))
                runner_yamls[task["id"]] = f"/work/{runner_path.relative_to(job.run_dir).as_posix()}"
        container_name = "mn-studio-predict-" + hashlib.sha256(job.run_id.encode("utf-8")).hexdigest()[:12]
        command = [
            "docker", "run", "--name", container_name, "--rm", "--cpus", str(int(cpu_cores)),
            "--shm-size=64G" if profile["api_engine"] == "esmfold" else "--shm-size=32G",
            *docker_gpu_args(gpu_text),
            "--mount", f"type=bind,src={job.run_dir.resolve()},dst=/work",
        ]
        if profile["api_engine"] == "esmfold":
            esmfold_cache = Path(reference_root()).expanduser().resolve() / "esm_models" / "esmfold_v1" / "torch"
            esmfold_cache.mkdir(parents=True, exist_ok=True)
            command.extend([
                "--mount", f"type=bind,src={esmfold_cache},dst=/reference/esmfold_v1_torch",
                "-e", "TORCH_HOME=/reference/esmfold_v1_torch",
            ])
            combined_fasta = raw_root / "inputs" / "esmfold_batch.fasta"
            combined_fasta.write_text(
                "".join(
                    f">{task['id']}\n{':'.join(chain['sequence'] for chain in task['proteins'])}\n"
                    for task in normalized_tasks
                ),
                encoding="utf-8",
            )
            command.extend([
                "-w", "/work", str(profile["image"]), "esm-fold",
                "-i", f"/work/{combined_fasta.relative_to(job.run_dir).as_posix()}",
                "-o", f"/work/{output_dir.relative_to(job.run_dir).as_posix()}",
                "--num-recycles", str(recycle_count),
                "--max-tokens-per-batch", str(int(selected_options.get("max_tokens_per_batch", 1024))),
            ])
            chunk_size = str(selected_options.get("chunk_size", "auto"))
            if chunk_size != "auto":
                command.extend(("--chunk-size", chunk_size))
            if selected_options.get("cpu_offload", False):
                command.append("--cpu-offload")
        else:
            reference = Path(reference_root()).expanduser().resolve()
            cache = Path(os.getenv("MN_MOLECULAR_DESIGN_STUDIO_CACHE_DIR", str(reference / "mn-cofolding" / "cache"))).expanduser().resolve()
            boltz_cache = reference / "boltz_models"
            if profile["api_engine"] == "boltz2":
                required_boltz_files = (
                    boltz_cache / "boltz2_conf.ckpt",
                    boltz_cache / "boltz2_aff.ckpt",
                    boltz_cache / "ccd.pkl",
                )
                missing_boltz_files = [str(path) for path in required_boltz_files if not path.is_file()]
                if not boltz_cache.is_dir() or missing_boltz_files or not (
                    (boltz_cache / "mols").is_dir() or (boltz_cache / "mols.tar").is_file()
                ):
                    missing = missing_boltz_files or [str(boltz_cache / "mols"), str(boltz_cache / "mols.tar")]
                    raise FileNotFoundError(
                        "Boltz-2 reference cache is incomplete; missing: " + ", ".join(missing)
                    )
            if profile["api_engine"] == "boltzgen_fold_opt_msa":
                cache = reference / "boltzgen-cache"
                required_cache = (
                    cache / "datasets--boltzgen--inference-data/snapshots/c3d36fd276e9caf098c75d4113c6d5eb320b1a4c/mols.zip",
                    cache / "models--boltzgen--boltzgen-1/snapshots/c1be29e1f82ffcc72264f64b993c43fb4e0d17f0/boltz2_conf_final.ckpt",
                )
                missing = [str(path) for path in required_cache if not path.is_file()]
                if missing:
                    raise FileNotFoundError("BoltzGen Fold cache is missing: " + ", ".join(missing))
            cache.mkdir(parents=True, exist_ok=True)
            command.extend((
                "--mount", f"type=bind,src={cache},dst=/cache",
                "--mount", f"type=bind,src={reference},dst=/reference,readonly",
            ))
            if profile["api_engine"] == "af3_torch":
                # The image intentionally ships without model parameters. Reuse the
                # converted OpenFold3 P2 weights stored by MN Cofolding, matching its
                # deployment contract and keeping the reference mount read-only.
                command.extend((
                    "--env", "HOME=/cache/home",
                    "--env", "XDG_CACHE_HOME=/cache/xdg",
                    "--env", "AF3_TORCH_PARAMS_DIR=/reference/mn-cofolding/af3_torch",
                    "--env", "AF3_TORCH_GPU=RTX4090",
                    "--env", "AF3_TORCH_CACHE_ROOT=/cache/af3_torch/torch2.13.0-cu130-sm89",
                    "--env", "MODEL_OPT_JIT_ROOT=/cache/af3_torch/jit",
                    "--env", f"OMP_NUM_THREADS={int(cpu_cores)}",
                ))
            if profile["api_engine"] == "boltz2":
                command.extend(("--mount", f"type=bind,src={boltz_cache},dst=/models"))
            options = _batch_options(
                profile, seed=int(seed), num_recycles=recycle_count, engine_options=selected_options,
                use_potentials=any(
                    bool(details.get("force"))
                    for task in normalized_tasks
                    for restraint in task["restraints"]
                    for details in restraint.values()
                ) or any(task["force_potentials"] for task in normalized_tasks),
            )
            if runner_yamls:
                options["runner_yamls"] = runner_yamls
            if chai_template_hits:
                hits_path = raw_root / "inputs/all_template_hits.m8"
                hits_path.write_text("".join((job.run_dir / value.removeprefix("/work/")).read_text() for value in chai_template_hits))
                options["template_hits_path"] = f"/work/{hits_path.relative_to(job.run_dir).as_posix()}"
            manifest = {
                "api_version": "mn-predict/v1", "batch_id": job.run_id,
                "engine": profile["api_engine"], "mode": selected_mode,
                "gpu_count": requested_gpu_count, "samples_per_item": int(samples),
                "seed": int(seed), "items": manifest_items, "options": options,
            }
            manifest_path = raw_root / "mn_predict_manifest.json"
            _write_json(manifest_path, manifest)
            batch_arguments = ["batch", "--manifest", f"/work/{manifest_path.relative_to(job.run_dir).as_posix()}",
                               "--output-dir", f"/work/{output_dir.relative_to(job.run_dir).as_posix()}"]
            patches = []
            if chai_template_hits:
                patches.append(structure_templates.stage_chai_runtime_patch(job.run_dir))
            if profile["api_engine"] == "alphafold3" and selected_options.get("model", "alphafold3") not in {"alphafold3", "openfold3"}:
                from mn_molecular_design_studio.engine.workflow_engine.workflows.af3_jax_cofolding import stage_batch_model_patch
                patches.append(stage_batch_model_patch(job.run_dir))
            if patches:
                scripts = [shlex.join(["python3", f"/work/{patch.relative_to(job.run_dir).as_posix()}"]) for patch in patches]
                script = " && ".join([*scripts, "exec " + shlex.join(["/usr/local/bin/mn-predict", *batch_arguments])])
                command.extend(["--entrypoint", "/bin/bash", "-w", "/work", str(profile["image"]), "-lc", script])
            else:
                command.extend(["-w", "/work", str(profile["image"]), *batch_arguments])
        from mn_molecular_design_studio.prediction_provenance import image_identity, file_identity
        installed_engine = image_identity(profile["image"])
        prediction_inputs = [file_identity(path) for path in (raw_root / "inputs").rglob("*") if path.is_file()]
        checkpoint_provenance = [file_identity(boltz_cache / name) for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "ccd.pkl")] if profile["api_engine"] == "boltz2" else []
        write_json(raw_root / "execution_inputs.json", {**installed_engine, "input_files": prediction_inputs, "checkpoints": checkpoint_provenance})
        write_json(job.run_dir / "command.json", {
            "mode": "docker", "command": command, "image": profile["image"],
            "api_version": "esmfold-cli/v1" if profile["api_engine"] == "esmfold" else "mn-predict/v1",
            "engine": profile["api_engine"],
        })
        write_json(raw_root / "input_mapping.json", {
            "engine": engine, "api_engine": profile["api_engine"], "image": profile["image"],
            "optimization_mode": selected_mode, "gpu_count": requested_gpu_count,
            "items": input_mappings,
        })
        update_status(job.run_dir, "running", current_phase="Running optimized structure prediction",
                      progress_label=profile["label"], queue_resource=gpu_queue_resource(gpu_text))
        with (job.run_dir / "stdout.log").open("a", encoding="utf-8") as stdout, (job.run_dir / "stderr.log").open("a", encoding="utf-8") as stderr:
            stdout.write("$ " + " ".join(command) + "\n")
            stdout.flush()
            process = subprocess.run(command, cwd=job.run_dir, stdout=stdout, stderr=stderr, check=False)
        if profile["api_engine"] == "esmfold":
            esmfold_rows = []
            for task in normalized_tasks:
                structure_path = output_dir / f"{task['id']}.pdb"
                if not structure_path.is_file():
                    matches = sorted(output_dir.glob(f"{task['id']}*.pdb"))
                    structure_path = matches[0] if matches else structure_path
                completed = process.returncode == 0 and structure_path.is_file()
                esmfold_rows.append({
                    "id": task["id"],
                    "status": "completed" if completed else "failed",
                    "outputs": [structure_path.relative_to(job.run_dir).as_posix()] if completed else [],
                    **({} if completed else {"error": "ESMFold did not produce a PDB for this input."}),
                })
            result_path = output_dir / "mn-batch-result.json"
            write_json(result_path, {"api_version": "esmfold-cli/v1", "items": esmfold_rows})
        if profile["api_engine"] == "boltz2":
            _cache_generated_boltz_msas(job.run_dir, output_dir)
        result_path = output_dir / "mn-batch-result.json"
        batch_result = json.loads(result_path.read_text(encoding="utf-8")) if result_path.is_file() else {}
        rows = batch_result.get("items") if isinstance(batch_result, dict) else []
        rows_by_id = {
            str(value.get("id") or ""): value
            for value in rows or []
            if isinstance(value, dict)
        }
        candidates = []
        structures: list[Path] = []
        item_results = []
        failures = []
        source_collections = {}
        for task in normalized_tasks:
            row = rows_by_id.get(task["id"])
            task_structures: list[Path] = []
            task_errors = []
            if not row:
                task_errors.append("No result row was returned for this task.")
            elif str(row.get("status") or "") != "completed":
                task_errors.append(str(row.get("error") or f"Batch status: {row.get('status', 'missing')}"))
            if row:
                raw_task_outputs = row.get("outputs")
                if not isinstance(raw_task_outputs, list):
                    raw_task_outputs = []
                for raw_path in raw_task_outputs:
                    raw_path = str(raw_path)
                    container_path = Path(raw_path)
                    if container_path.is_absolute():
                        try:
                            container_path = container_path.relative_to("/work")
                        except ValueError as exc:
                            raise ValueError(f"Prediction output escaped the job mount: {raw_path}") from exc
                    host_path = (job.run_dir / container_path).resolve()
                    if job.run_dir.resolve() not in host_path.parents or not host_path.is_file():
                        raise ValueError(f"Prediction output is missing or outside the job folder: {raw_path}")
                    task_structures.append(host_path)
            if not task_structures and not task_errors:
                task_errors.append("The completed task returned no structure files.")
            if task_errors:
                failures.append({"item_id": task["id"], "task_name": task["name"], "errors": task_errors})
            target_chains = list(task["target_protein_ids"])
            binder_chains = list(task["binder_protein_ids"])
            for sample_index, structure_path in enumerate(task_structures, start=1):
                from mn_molecular_design_studio.scientific_contracts import verify_prediction_structure, restore_prediction_ccd_names
                native_structure_path = structure_path
                confidence = _confidence_sidecar(native_structure_path)
                structure_path = restore_prediction_ccd_names(native_structure_path, task)
                identity_check = verify_prediction_structure(structure_path, task)
                structures.append(structure_path)
                pae_path = next((value for value in (
                    native_structure_path.with_name(f"pae_{native_structure_path.stem}.npz"),
                    native_structure_path.with_name(f"{native_structure_path.stem}_pae.npz"),
                    native_structure_path.with_name(f"{native_structure_path.stem}_pae.json"),
                    native_structure_path.with_name(f"pae_{native_structure_path.stem}.json"),
                ) if value.is_file()), None)
                metrics = {
                    "plddt_mean": _mean_plddt(structure_path),
                    "ptm": _metric_float(confidence.get("ptm")),
                    "iptm": _metric_float(confidence.get("iptm")),
                    "ranking_score": _metric_float(confidence.get("ranking_score", confidence.get("confidence_score"))),
                    "affinity_pred_value": _metric_float(confidence.get("affinity_pred_value")),
                    "affinity_probability_binary": _metric_float(confidence.get("affinity_probability_binary")),
                    "prediction_sample": sample_index,
                    "optimization_mode": selected_mode,
                    "engine": engine,
                }
                metrics = {key: value for key, value in metrics.items() if value is not None}
                candidates.append({
                    "candidate_id": f"{job.run_id}_{engine}_{task['index']:02d}_{sample_index:03d}",
                    "stage": STAGE_MONOMER_REFOLDING if validation_mode == "binder" else STAGE_COMPLEX_REFOLDING,
                    "source_tool": "optimized_prediction_" + str(engine),
                    "complex_pdb": str(structure_path),
                    "binder_sequence": task["binder_sequence"] if is_boltzgen or binder_chains else None,
                    "target_chains": [task["target_chain_id"]] if is_boltzgen else target_chains,
                    "binder_chains": ["Z"] if is_boltzgen else binder_chains,
                    "metrics": metrics,
                    "parents": [task["source_candidate_id"]] if task.get("source_candidate_id") else [],
                    "raw_metadata": {
                        "scientific_contract": task["scientific_contract"],
                        "native_prediction_structure": str(native_structure_path),
                        "pae_path": str(pae_path) if pae_path else None,
                        "identity_check": identity_check,
                        "molecule_chain_options": task["molecule_chain_options"],
                        "engine": engine, "image": profile["image"], "api_engine": profile["api_engine"],
                        "optimization_mode": selected_mode, "seed": int(seed),
                        **installed_engine,
                        "input_provenance": str(raw_root / "execution_inputs.json"),
                        "checkpoint_provenance": checkpoint_provenance,
                        "engine_version": profile.get("code_version"),
                        "conditioning_mode": "target_template" if (task.get("structure_templates") or task.get("boltz_templates") or task.get("target_pdb_text")) else "sequence_only",
                        "effective_settings": {"samples_per_item": int(samples), "seed": int(seed), **options} if profile["api_engine"] == "boltz2" else {},
                        "execution_manifest": str(manifest_path) if profile["api_engine"] == "boltz2" else None,
                        "native_execution": batch_result.get("execution", {}) if profile["api_engine"] == "boltz2" else {},
                        "prediction_protocol": {"engine": engine, "image": profile["image"], "optimization_mode": selected_mode,
                            "samples": int(samples), "seed": int(seed),
                            "num_recycles": recycle_count, "engine_options": selected_options, "validation_mode": validation_mode,
                            "target_use_msa": bool(target_use_msa), "use_target_template": bool(use_target_template),
                            "molecule_chain_options": task["molecule_chain_options"], "restraints": task["restraints"]},
                        "target_sequences": [chain["sequence"] for chain in task["proteins"] if chain["id"] in target_chains],
                        "confidence_outputs": confidence, "batch_item_id": task["id"],
                        "input_name": task["source_input_name"], "source_compound_id": task["source_compound_id"],
                        "source_compound_name": task["source_compound_name"],
                        "compound_id": task["compound_id"], "compound_name": task["compound_name"],
                        "ligand_smiles": task["compound_smiles"],
                        "ligand_ccd_codes": task["ligand_ccd_codes"],
                        **({"source_run_dir": task["source_run_dir"]} if task.get("source_run_dir") else {}),
                        "protein_chains": [
                            {
                                "id": chain["id"],
                                "sequence": chain["sequence"],
                                "role": (
                                    "binder" if chain["id"] in binder_chains
                                    else "target" if chain["id"] in target_chains
                                    else "unspecified"
                                ),
                            }
                            for chain in (task["proteins"] or [
                                {"id": chain_id, "sequence": protein.get("sequence", "")}
                                for protein in task["raw_yaml_proteins"]
                                for chain_id in (protein["id"] if isinstance(protein.get("id"), list) else [protein.get("id")])
                                if chain_id
                            ])
                        ],
                    },
                })
                reference = task.get("validation_reference") or {}
                source_run = task.get("source_run_dir") or reference.get("source_run_dir")
                source_id = task.get("source_candidate_id") or reference.get("candidate_id")
                if source_run and source_id:
                    if source_run not in source_collections:
                        source_collections[source_run] = {item["candidate_id"]: item for item in read_candidates(Path(source_run))}
                    if source_id in source_collections[source_run]:
                        candidates[-1]["raw_metadata"]["source_candidate"] = source_collections[source_run][source_id]
                if task["validation_reference"]:
                    candidate = candidates[-1]
                    reference = task["validation_reference"]
                    candidate["parents"] = [reference["candidate_id"]]
                    candidate["raw_metadata"]["source_run_dir"] = reference["source_run_dir"]
                    candidate["raw_metadata"]["validation_reference"] = reference
                    try:
                        comparison, artifacts = design_validation.compare_prediction(
                            job.run_dir, task["id"], sample_index, structure_path, reference)
                        candidate["metrics"].update(comparison)
                        candidate["raw_metadata"].update(artifacts)
                    except (ValueError, KeyError, OSError) as exc:
                        candidate["metrics"].update(validation_status="failed", validation_error=str(exc))
                        failures.append({"scope": "validation", "task_name": task["name"], "error": str(exc)})
            item_results.append({
                "item_id": task["id"], "task_name": task["name"],
                "status": str((row or {}).get("status") or "missing"),
                "structure_count": len(task_structures),
                "error": "; ".join(task_errors) if task_errors else None,
            })
        if process.returncode != 0:
            failures.append({
                "scope": "container", "error": str(batch_result.get("error") or f"Container exited with code {process.returncode}.")
            })
        success = process.returncode == 0 and not failures and len(item_results) == len(normalized_tasks)
        if supports_boltz_features and selected_options.get("generate_ipsae", True) and candidates:
            update_status(job.run_dir, "running", current_phase="Calculating protein interface ipSAE metrics")
            _add_boltz_interface_metrics(job.run_dir, candidates, workers=int(cpu_cores))
        if candidates:
            write_candidates(job.run_dir, "optimized_prediction_" + str(engine), candidates)
        validation_report = design_validation.write_report(job.run_dir, candidates) if validation_sources is not None else None
        failure_text = "; ".join(
            f"{failure.get('task_name') or failure.get('scope')}: {'; '.join(failure.get('errors') or [failure.get('error', 'failed')])}"
            for failure in failures
        )
        finish_job(job.run_dir, success, {
            "outputs": {
                "structures": [str(path) for path in structures],
                "batch_manifest": str(result_path),
                "candidates": candidates,
                "items": item_results,
                "validation_report": str(validation_report) if validation_report else None,
            },
            "metrics": {
                "engine": engine, "image": profile["image"], "optimization_mode": selected_mode,
                "model_label": profile["label"],
                "gpu_count": requested_gpu_count, "task_count": len(normalized_tasks),
                "completed_task_count": sum(item["status"] == "completed" and item["structure_count"] > 0 for item in item_results),
                "failed_task_count": sum("item_id" in failure for failure in failures),
                "container_failure_count": sum(failure.get("scope") == "container" for failure in failures),
                "comparison_failure_count": sum(failure.get("scope") == "validation" for failure in failures),
                "sample_count": len(structures),
                "candidate_count": len(candidates), "engine_options": selected_options,
            },
            "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl"},
            **({"error": failure_text} if not success else {}),
        })
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a", encoding="utf-8", errors="replace") as stderr:
            stderr.write("\n" + traceback.format_exc() + "\n")
        finish_job(job.run_dir, False, {
            "metrics": {"engine": engine, "image": profile["image"], "optimization_mode": selected_mode,
                        "prediction_error": str(exc)},
            "error": f"{type(exc).__name__}: {exc}",
        })
    return job.run_dir
