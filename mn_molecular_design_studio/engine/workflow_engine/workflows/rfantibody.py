"""Queued, standalone RFantibody design using the pinned upstream CLI."""
from __future__ import annotations

import csv
import re
import secrets
import shlex
import shutil
import subprocess
import traceback
from pathlib import Path
from typing import Literal

from Bio.PDB import PDBIO, Select

from mn_molecular_design_studio import rfantibody_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import STAGE_COMPLEX_REFOLDING, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, update_status, write_json
from mn_molecular_design_studio.runtime import reference_root

STANDARD_RESIDUES = frozenset({"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"})
PDB_CHAIN_IDS = "".join(char for char in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" if char not in "HLT")


class _SelectedProtein(Select):
    def __init__(self, chain_ids: set[str]):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

    def accept_residue(self, residue):
        return residue.id[0] == " " and residue.resname in STANDARD_RESIDUES and "CA" in residue

    def accept_atom(self, atom):
        return atom.element != "H" and (not atom.is_disordered() or atom.get_altloc() in {"A", " "})


def _write_selected_structure(structure, destination: Path, chain_ids: set[str]) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(destination), _SelectedProtein(chain_ids))
    if not destination.is_file() or "ATOM" not in destination.read_text(errors="ignore"):
        raise ValueError(f"No standard protein atoms could be written to {destination.name}.")


def _prepare_target(source: Path, selected_chains: list[str], hotspots: list[str], destination: Path) -> dict:
    structure = read_structure(source)
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("The selected target structure has no coordinate model.")
    chains = protein_chains(structure)
    missing = [chain for chain in selected_chains if chain not in chains]
    if missing:
        raise ValueError(f"Selected antigen chain(s) are missing from the target: {', '.join(missing)}.")

    used: set[str] = set()
    chain_map: dict[str, str] = {}
    for source_chain in selected_chains:
        if len(source_chain) == 1 and source_chain in PDB_CHAIN_IDS and source_chain not in used:
            output_chain = source_chain
        else:
            output_chain = next((char for char in PDB_CHAIN_IDS if char not in used), "")
        if not output_chain:
            raise ValueError("The selected target uses more chains than the PDB format can represent.")
        used.add(output_chain)
        chain_map[source_chain] = output_chain

    mapped_hotspots = []
    for hotspot in hotspots:
        match = re.fullmatch(r"([A-Za-z0-9]+?)(-?\d+)", hotspot)
        source_chain, residue_number = match.groups()
        if source_chain not in chain_map:
            raise ValueError(f"Hotspot {hotspot} is not part of the selected antigen chains.")
        residues = chains[source_chain]
        if not any(int(residue.id[1]) == int(residue_number) and not str(residue.id[2]).strip() for residue in residues):
            raise ValueError(f"Hotspot {hotspot} is missing or has an insertion code, which RFantibody cannot address.")
        mapped_hotspots.append(f"{chain_map[source_chain]}{int(residue_number)}")

    for chain in list(model):
        if chain.id not in chain_map:
            model.detach_child(chain.id)
            continue
        original_id = chain.id
        chain.id = chain_map[original_id]
        for residue in list(chain):
            if residue.id[0] != " " or residue.resname not in STANDARD_RESIDUES or "CA" not in residue:
                chain.detach_child(residue.id)
    _write_selected_structure(structure, destination, set(chain_map.values()))
    return {"chain_map": chain_map, "target_chains": list(chain_map.values()), "hotspots": mapped_hotspots}


def _score_from_pdb(path: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    for line in path.read_text(errors="ignore").splitlines():
        match = re.match(r"SCORE\s+([A-Za-z0-9_]+):\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)", line)
        if match:
            scores[match.group(1)] = float(match.group(2))
    return scores


def _write_binder_structure(source: Path, destination: Path, chains: list[str]) -> None:
    structure = read_structure(source)
    model = next(structure.get_models(), None)
    if model is None:
        return
    allowed = set(chains)
    for chain in list(model):
        if chain.id not in allowed:
            model.detach_child(chain.id)
    _write_selected_structure(structure, destination, allowed)


def _collect_candidates(
    run_dir: Path,
    target_pdb: Path,
    target_chain_ids: list[str],
    hotspots: list[str],
    *,
    calculate_epitope_distance: bool,
) -> list[dict]:
    output_root = run_dir / "artifacts" / "raw" / "rfantibody"
    structures_dir = output_root / "structures"
    candidates = []
    pdb_paths = sorted(output_root.glob("rf2_models/model_*/*_best.pdb"))
    for index, path in enumerate(pdb_paths, 1):
        try:
            structure = read_structure(path)
            chains = protein_chains(structure)
            binder_chains = [chain for chain in ("H", "L") if chain in chains]
            target_chains = [chain for chain in chains if chain not in binder_chains]
            if not binder_chains or not target_chains:
                continue
            binder_sequence = "".join(sequence(chains[chain]) for chain in binder_chains)
            binder_path = structures_dir / f"{path.stem}_binder.pdb"
            _write_binder_structure(path, binder_path, binder_chains)
            metrics = _score_from_pdb(path)
            model_match = re.search(r"model_(\d+)", path.parent.name)
            raw_hotspots = hotspots if calculate_epitope_distance else []
            candidates.append({
                "candidate_id": path.stem,
                "source_tool": "rfantibody",
                "stage": STAGE_COMPLEX_REFOLDING,
                "target_pdb": str(target_pdb),
                "complex_pdb": str(path),
                "binder_pdb": str(binder_path),
                "binder_sequence": binder_sequence,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": raw_hotspots,
                "binder_length": len(binder_sequence),
                "metrics": {
                    **metrics,
                    "rf2_repeat": int(model_match.group(1)) if model_match else None,
                    "binder_length": len(binder_sequence),
                    "hotspot_distance_analysis": bool(calculate_epitope_distance),
                },
                "raw_metadata": {
                    "upstream_commit": inputs.UPSTREAM_COMMIT,
                    "input_target_chains": target_chain_ids,
                    "input_hotspots": hotspots,
                    "rf2_repeat": int(model_match.group(1)) if model_match else None,
                    "result_kind": "rf2_complex_prediction",
                },
            })
        except (OSError, ValueError, KeyError) as exc:
            (output_root / "collection_warnings.txt").open("a").write(f"{path.name}: {type(exc).__name__}: {exc}\n")
    return sorted(candidates, key=lambda candidate: (-float(candidate["metrics"].get("pred_lddt", -1)), float(candidate["metrics"].get("interaction_pae", 1e9))))


def run_rfantibody_design(
    task: Literal["antibody", "nanobody"] = "antibody",
    target_pdb: Path | None = None,
    target_chains: list[str] | None = None,
    hotspots: str = "",
    framework: str = "hu-4D5-8_Fv",
    framework_pdb: Path | None = None,
    custom_heavy_chain: str = "H",
    custom_light_chain: str = "L",
    regions: list[str] | None = None,
    loop_lengths: dict[str, str] | None = None,
    number_of_designs: int = 1,
    number_of_sequences: int = 1,
    use_abmpnn: bool = False,
    calculate_epitope_distance: bool = True,
    mpnn_temperature: float = 0.1,
    bias_aas: str = "",
    bias_by_res_json: str = "{}",
    omit_aas: str = "CX",
    augment_eps: float = 0.05,
    num_connections: int = 48,
    num_recycles: int = 3,
    rf2_repeats: int = 5,
    hotspot_show_prop: float = 0.1,
    save_trajectories: bool = False,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    job_name: str = "rfantibody-design",
) -> Path:
    parameters = dict(locals())
    normalized = inputs.validate_request(parameters, require_files=True)
    if normalized["target_pdb"] is None:
        raise ValueError("Load a target PDB or mmCIF structure first.")
    weights_dir = reference_root() / "rfantibody"
    required_weights = (
        "RFdiffusion_Ab.pt",
        "ProteinMPNN_v48_noise_0.2.pt",
        "RF2_ab.pt",
    )
    missing_weights = [name for name in required_weights if not (weights_dir / name).is_file()]
    if missing_weights:
        raise ValueError(
            f"RFantibody checkpoints are missing from {weights_dir}: "
            + ", ".join(missing_weights)
        )
    seed = secrets.randbelow(2_000_000_000) + 1
    gpu_device = normalize_gpu_device(gpu_device)
    parameters.update({
        "target_pdb": str(normalized["target_pdb"]),
        "framework_pdb": str(normalized["framework_pdb"]) if normalized["framework_pdb"] else None,
        "gpu_device": gpu_device,
        "cpu_cores": normalized["cpu_cores"],
        "seed": seed,
        "engine_version": inputs.CODE_VERSION,
        "image": inputs.IMAGE,
    })
    inputs_payload = {"target_pdb": str(normalized["target_pdb"]), "framework_pdb": str(normalized["framework_pdb"]) if normalized["framework_pdb"] else None}
    job = create_job("design", "rfantibody_design", "rfantibody", inputs_payload, parameters)
    update_status(job.run_dir, "preparing", campaign_name=job_name, current_phase="Preparing RFantibody target and framework")
    root = job.run_dir / "artifacts" / "raw" / "rfantibody"
    inputs_dir = root / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    target_path = inputs_dir / "target.pdb"
    try:
        mapping = _prepare_target(normalized["target_pdb"], normalized["target_chains"], normalized["hotspots"], target_path)
        shutil.copy2(normalized["target_pdb"], inputs_dir / f"source_target{normalized['target_pdb'].suffix.lower()}")
        if normalized["framework"] == "custom":
            shutil.copy2(normalized["framework_pdb"], inputs_dir / "custom_framework.pdb")

        design_loops = normalized["design_loops"]
        hotspot_text = ",".join(mapping["hotspots"])
        diff_qv = root / "stages" / "1_rfdiffusion.qv"
        mpnn_qv = root / "stages" / "2_proteinmpnn.qv"
        predictions_root = root / "rf2_models"
        output_root = root / "structures"
        (diff_qv.parent).mkdir(parents=True, exist_ok=True)
        predictions_root.mkdir(parents=True, exist_ok=True)
        output_root.mkdir(parents=True, exist_ok=True)
        target_in_container = "/work/artifacts/raw/rfantibody/inputs/target.pdb"
        diff_qv_container = "/work/artifacts/raw/rfantibody/stages/1_rfdiffusion.qv"
        mpnn_qv_container = "/work/artifacts/raw/rfantibody/stages/2_proteinmpnn.qv"
        if normalized["framework"] == "custom":
            framework_container = "/work/artifacts/raw/rfantibody/inputs/framework_hlt.pdb"
            framework_command = [
                "python", "scripts/util/chothia2HLT.py", "/work/artifacts/raw/rfantibody/inputs/custom_framework.pdb",
                "--heavy", normalized["custom_heavy_chain"], "--output", framework_container,
            ]
            if normalized["task"] == "antibody":
                framework_command.extend(["--light", normalized["custom_light_chain"]])
        else:
            framework_name = "h-NbBCII10" if normalized["framework"] == "h-NbBCII10" else "hu-4D5-8_Fv"
            framework_container = f"scripts/examples/example_inputs/{framework_name}.pdb"
            framework_command = []

        image_mounts = [
            "--mount", f"type=bind,src={job.run_dir.resolve()},dst=/work",
            "--mount", f"type=bind,src={weights_dir.resolve()},dst=/ref/rfantibody,readonly",
        ]
        mpnn_checkpoint_in_container = "/ref/rfantibody/ProteinMPNN_v48_noise_0.2.pt"
        if normalized["use_abmpnn"]:
            checkpoint = normalized["abmpnn_checkpoint"]
            image_mounts.extend(["--mount", f"type=bind,src={checkpoint.resolve()},dst=/mnt/studio-weights/abmpnn.pt,readonly"])
            mpnn_checkpoint_in_container = "/mnt/studio-weights/abmpnn.pt"

        command = [
            "docker", "run", "--rm", "--name", f"mn-rfantibody-{job.run_id[-8:]}",
            "--cpus", str(normalized["cpu_cores"]), "--shm-size=16G", *docker_gpu_args(gpu_device),
            *image_mounts, "-w", "/opt/rfantibody",
            "-e", f"OMP_NUM_THREADS={normalized['cpu_cores']}", "-e", "PYTHONUNBUFFERED=1",
            "-e", "RFANTIBODY_WEIGHTS=/ref/rfantibody",
            inputs.IMAGE, "bash", "-euo", "pipefail", "-c",
        ]
        steps: list[str] = []
        if framework_command:
            steps.append("echo '[RFantibody] Converting the custom antibody framework'; " + shlex.join(framework_command))
        diffusion = [
            "rfdiffusion", "--target", target_in_container, "--framework", framework_container,
            "--output-quiver", diff_qv_container, "--num-designs", str(normalized["number_of_designs"]),
            "--design-loops", design_loops, "--diffuser-t", "50",
        ]
        if hotspot_text:
            diffusion.extend(["--hotspots", hotspot_text])
        if not save_trajectories:
            diffusion.append("--no-trajectory")
        steps.extend(["echo '[RFantibody] Stage 1/3 · RFdiffusion'; " + shlex.join(diffusion)])
        mpnn = [
            "proteinmpnn", "--input-quiver", diff_qv_container, "--output-quiver", mpnn_qv_container,
            "--loops", ",".join(normalized["regions"]), "--seqs-per-struct", str(normalized["number_of_sequences"]),
            "--temperature", str(normalized["mpnn_temperature"]), "--augment-eps", str(normalized["augment_eps"]),
            "--num-connections", str(normalized["num_connections"]), "--omit-aas", normalized["omit_aas"],
            "--bias-aas", normalized["bias_aas"], "--bias-by-res-json", normalized["bias_by_res_json"],
        ]
        if use_abmpnn:
            mpnn.extend(["--weights", mpnn_checkpoint_in_container])
        steps.append("echo '[RFantibody] Stage 2/3 · ProteinMPNN'; " + shlex.join(mpnn))
        for repeat in range(1, normalized["rf2_repeats"] + 1):
            model_output = f"/work/artifacts/raw/rfantibody/rf2_models/model_{repeat:02d}"
            rf2 = [
                "rf2", "--input-quiver", mpnn_qv_container, "--output-dir", model_output,
                "--num-recycles", str(normalized["num_recycles"]), "--hotspot-show-prop", str(normalized["hotspot_show_prop"]),
                "--seed", str(seed + repeat - 1), "--no-cautious",
            ]
            steps.append(f"echo '[RFantibody] Stage 3/3 · RF2 prediction {repeat}/{normalized['rf2_repeats']}'; " + shlex.join(rf2))
        script = "\n".join(steps)
        write_json(job.run_dir / "command.json", {"mode": "docker", "image": inputs.IMAGE, "command": command + [script]})
        update_status(job.run_dir, "running", current_phase="RFdiffusion → ProteinMPNN → RF2")
        with (job.run_dir / "stdout.log").open("a") as stdout, (job.run_dir / "stderr.log").open("a") as stderr:
            process = subprocess.run(command + [script], stdout=stdout, stderr=stderr, cwd=job.run_dir, check=False)
        update_status(job.run_dir, "running", current_phase="Collecting RFantibody designs and RF2 scores")
        candidates = _collect_candidates(
            job.run_dir, target_path, mapping["target_chains"], mapping["hotspots"],
            calculate_epitope_distance=calculate_epitope_distance,
        )
        normalized_candidates = write_candidates(job.run_dir, "rfantibody", candidates) if candidates else []
        table_path = root / "designs.csv"
        with table_path.open("w", newline="") as handle:
            fields = ["rank", "candidate_id", "binder_sequence", "binder_chains", "pred_lddt", "interaction_pae", "rf2_repeat", "structure"]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for rank, candidate in enumerate(normalized_candidates, 1):
                metrics = candidate.get("metrics") or {}
                writer.writerow({"rank": rank, "candidate_id": candidate["candidate_id"], "binder_sequence": candidate.get("binder_sequence"),
                                 "binder_chains": ",".join(candidate.get("binder_chains") or []), "pred_lddt": metrics.get("pred_lddt"),
                                 "interaction_pae": metrics.get("interaction_pae"), "rf2_repeat": metrics.get("rf2_repeat"),
                                 "structure": candidate.get("complex_pdb")})
        artifacts = [{"path": str(path), "name": path.name, "type": "rfantibody_output"}
                     for path in root.rglob("*") if path.is_file()]
        success = process.returncode == 0 and bool(normalized_candidates)
        finish_job(job.run_dir, success, {
            "outputs": {"artifacts": artifacts, "candidates": normalized_candidates, "design_table": str(table_path)},
            "metrics": {"return_code": process.returncode, "candidate_count": len(normalized_candidates),
                        "requested_designs": normalized["number_of_designs"], "sequences_per_structure": normalized["number_of_sequences"],
                        "rf2_repeats": normalized["rf2_repeats"], "engine_version": inputs.CODE_VERSION,
                        "selected_target_chain_map": mapping["chain_map"], "selected_hotspots": mapping["hotspots"]},
            "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                                     "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
                                     "design_table": str(table_path.relative_to(job.run_dir))},
        })
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a") as handle:
            handle.write(traceback.format_exc())
        finish_job(job.run_dir, False, {"error": str(exc), "metrics": {"error": str(exc), "engine_version": inputs.CODE_VERSION}})
    return job.run_dir
