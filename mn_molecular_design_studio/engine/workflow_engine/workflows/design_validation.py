"""Design-source adapter for AlphaFold2 initial-guess validation."""
from __future__ import annotations

import copy
import tempfile
from pathlib import Path

from Bio.PDB import Model, PDBIO, Structure

from mn_molecular_design_studio import design_validation as comparison
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, STAGE_COMPLEX_REFOLDING


def validate_af2_design_request(parameters: dict) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
    tasks = comparison.build_tasks(parameters.get("validation_sources"), mode="complex", profile=ENGINE_PROFILES["colabfold"],
                                   target_use_msa=parameters.get("target_use_msa", True), use_target_template=False)
    if any(len(task["validation_reference"]["binder_chains"]) != 1 for task in tasks):
        raise ValueError("AlphaFold2 Initial Guess currently supports one designed binder chain per complex; use another cofolding engine for multichain binders.")
    for key, default, maximum in (("model_count", 1, 5), ("num_recycles", 3, 100), ("msa_max_depth", 1024, 4096), ("cpu_cores", 4, 256)):
        if not 1 <= int(parameters.get(key, default)) <= maximum:
            raise ValueError(f"{key} must be between 1 and {maximum}.")


def run_af2_design_validation(validation_sources: list[dict], job_name: str = "af2-design-validation",
                              num_recycles: int = 3, model_count: int = 1, multimer: bool = True,
                              target_use_msa: bool = True, msa_max_depth: int = 1024,
                              cpu_cores: int = 4, gpu_device: str = "0") -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
    from mn_molecular_design_studio.engine.workflow_engine.workflows.refolding import run_af2_initial_guess_complex_refolding
    validate_af2_design_request({"validation_sources": validation_sources, "target_use_msa": target_use_msa,
                                 "model_count": model_count, "num_recycles": num_recycles,
                                 "msa_max_depth": msa_max_depth, "cpu_cores": cpu_cores})
    tasks = comparison.build_tasks(validation_sources, mode="complex", profile=ENGINE_PROFILES["colabfold"],
                                   target_use_msa=target_use_msa, use_target_template=False)
    with tempfile.TemporaryDirectory(prefix="mds-af2-design-validation-") as temporary:
        root = Path(temporary)
        candidates = []
        for index, task in enumerate(tasks):
            reference = task["validation_reference"]
            source_dir = Path(reference["source_run_dir"])
            source = next(row for row in read_candidates(source_dir) if row["candidate_id"] == reference["candidate_id"])
            structure = comparison.read_structure(Path(reference["reference_structure_path"]))
            original_model = next(structure.get_models())
            normalized = Structure.Structure("design"); model = Model.Model(0); normalized.add(model)
            chain_map = {}
            for i, original_id in enumerate([*reference["target_chains"], *reference["binder_chains"]]):
                chain = copy.deepcopy(original_model[original_id]); chain.detach_parent()
                chain.id = comparison.CHAIN_IDS[i]; chain_map[original_id] = chain.id
                model.add(chain)
            pdb = root / f"design_{index}.pdb"
            writer = PDBIO(); writer.set_structure(normalized); writer.save(str(pdb))
            candidates.append({**source, "candidate_id": f"selected_design_{index}", "stage": STAGE_COMPLEX_REFOLDING,
                               "complex_pdb": str(pdb), "target_pdb": str(pdb), "binder_pdb": str(pdb),
                               "binder_sequence": task["binder_sequence"],
                               "target_chains": [chain_map[c] for c in reference["target_chains"]],
                               "binder_chains": [chain_map[c] for c in reference["binder_chains"]],
                               "raw_metadata": {**source.get("raw_metadata", {}), "validation_reference": reference,
                                                "validation_binder_single_sequence": True}})
        return run_af2_initial_guess_complex_refolding(
            source_run_dir=root, candidates_jsonl=root / "selected.jsonl", require_monomer_success=False,
            num_recycles=num_recycles, model_count=model_count, multimer=multimer,
            use_initial_guess=False, use_binder_template=False, use_interface_template=False,
            msa_source="reference_then_colabfold" if target_use_msa else "single_sequence", msa_max_depth=msa_max_depth,
            gpu_device=gpu_device, direct_candidates=candidates,
            direct_input={"job_name": job_name, "validation_sources": validation_sources})
