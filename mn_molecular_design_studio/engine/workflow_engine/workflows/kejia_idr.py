"""Kejia template threading as a CPU generator in the Studio scheduler."""
from __future__ import annotations

import copy
from pathlib import Path
import secrets
import subprocess
import traceback
from typing import Literal

from Bio.PDB import Chain, Model, PDBIO, Select, Structure

from mn_molecular_design_studio import kejia_idr_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import STAGE_GENERATION_BACKBONE, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, update_status, write_json


def stage_template(source: Path, destination: Path, binder_chain: str, peptide_chain: str) -> None:
    """Native Kejia assumes pose chain 1 is binder and chain 2 is peptide."""
    chains = protein_chains(read_structure(source))
    structure = Structure.Structure("kejia")
    model = Model.Model(0)
    structure.add(model)
    for selected, identifier in ((binder_chain, "A"), (peptide_chain, "B")):
        chain = Chain.Chain(identifier)
        for index, residue in enumerate(chains[selected], 1):
            residue = copy.deepcopy(residue)
            residue.detach_parent()
            residue.id = (" ", index, " ")
            chain.add(residue)
        model.add(chain)
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(destination))


class _TargetOnly(Select):
    def accept_chain(self, chain):
        return chain.id == "B"


def normalize_output(output: Path, target_sequence: str, template: dict) -> list[dict]:
    rows = []
    for index, record in enumerate(read_json(output / "threading.json").get("records", []), 1):
        path = output / record["pdb"]
        if path.parent.resolve() != output.resolve() or not path.is_file():
            raise ValueError("The threading result references an invalid structure path.")
        structure = read_structure(path)
        chains = protein_chains(structure)
        scores = record["scores"]
        start, end = int(scores["target_peptide_start"]), int(scores["target_peptide_end"])
        if not 1 <= start <= end <= len(target_sequence) or sequence(chains.get("B", [])) != target_sequence[start - 1:end]:
            raise ValueError("The threaded peptide sequence does not match its reported target fragment.")
        target = path.with_suffix(".target.pdb")
        writer = PDBIO(); writer.set_structure(structure); writer.save(str(target), _TargetOnly())
        rows.append({
            "candidate_id": f"kejia_{output.name}_{index:05d}", "stage": STAGE_GENERATION_BACKBONE,
            "source_tool": "kejia_idr", "complex_pdb": str(path), "target_pdb": str(target),
            "binder_chains": ["A"], "target_chains": ["B"], "binder_sequence": None,
            "metrics": {key: value for key, value in scores.items() if isinstance(value, (int, float))} |
                       {"sequence_design_required": True, "result_kind": "generated"},
            "raw_metadata": {"engine_version": inputs.CODE_VERSION, "source_revision": inputs.SOURCE_REVISION,
                             "template": template, "native_tag": record["tag"], "native_scores": scores,
                             "requested_target_sequence": target_sequence,
                             "threaded_target_sequence": record["target_sequence"],
                             "target_fragment_range": [start, end], "sequence_design_required": True},
        })
    return rows


def run_kejia_idr_threading(
    target_sequence: str,
    template_structure: Path | None = None,
    template_library: Literal["upload", "walle", "walle2", "mini", "strand"] = "upload",
    binder_chain: str = "A", peptide_chain: str = "B", max_templates: int = 10,
    flex_docking_repeats: int = 3, min_seq_identity: float = 0.0,
    min_blosum: float = -10.0, max_farep: float = 10.0,
    keep_bidentates: bool = False, seed: int = 0,
    job_name: str = "kejia-idr-threading", cpu_cores: int = 4,
) -> Path:
    parameters = dict(locals())
    parameters["template_structure"] = str(template_structure) if template_structure else None
    prepared = inputs.prepare_request(parameters)
    seed = seed or secrets.randbelow(2147483646) + 1
    parameters.update(target_sequence=prepared["target_sequence"], seed=seed, image=inputs.IMAGE,
                      engine_version=inputs.CODE_VERSION, gpu_device="none", generator_only=True)
    job = create_job("design", "kejia_idr_threading", "kejia_idr",
                     {"template_structure": parameters["template_structure"], "template_library": template_library}, parameters)
    root = job.run_dir / "artifacts/raw/kejia_idr"
    root.mkdir(parents=True, exist_ok=True)
    commands, candidates = [], []
    rejected = 0
    try:
        for index, template in enumerate(prepared["templates"], 1):
            identifier = f"template_{index:05d}"
            update_status(job.run_dir, "running", current_phase=f"Threading template {index}/{len(prepared['templates'])}")
            staged = root / "inputs" / f"{identifier}.pdb"
            inputs.template_chains(Path(template["path"]))
            stage_template(Path(template["path"]), staged, template["binder_chain"], template["peptide_chain"])
            output = root / identifier
            output.mkdir()
            request = {key: prepared[key] for key in ("target_sequence", "flex_docking_repeats", "min_seq_identity", "min_blosum", "max_farep")}
            request.update(template="/work/" + staged.relative_to(job.run_dir).as_posix(),
                           seed=(seed + index - 1) % 2147483647 or 1, keep_bidentates=keep_bidentates)
            request_path = root / "inputs" / f"{identifier}.json"
            write_json(request_path, request)
            command = ["docker", "run", "--rm", "--network", "none", "--name", f"mn-kejia-{job.run_id[-8:]}-{index}",
                       "--cpus", str(cpu_cores),
                       "--mount", f"type=bind,src={job.run_dir.resolve()},dst=/work", inputs.IMAGE,
                       "--request", "/work/" + request_path.relative_to(job.run_dir).as_posix(),
                       "--output", "/work/" + output.relative_to(job.run_dir).as_posix()]
            commands.append({"template": template, "command": command})
            write_json(job.run_dir / "command.json", {"mode": "docker", "image": inputs.IMAGE, "steps": commands})
            with (job.run_dir / "stdout.log").open("a") as stdout, (job.run_dir / "stderr.log").open("a") as stderr:
                result = subprocess.run(command, stdout=stdout, stderr=stderr, cwd=job.run_dir, check=False)
            if result.returncode or not (output / "threading.json").is_file():
                raise RuntimeError(f"Kejia threading failed on template {index}; see the job logs.")
            rejected += int(read_json(output / "threading.json").get("clash_rejected", 0))
            candidates.extend(normalize_output(output, prepared["target_sequence"], template))
        normalized = write_candidates(job.run_dir, "kejia_idr", candidates)
        finish_job(job.run_dir, True, {
            "outputs": {"candidates": normalized, "artifacts": [{"path": str(p), "name": p.name, "type": "kejia_output"} for p in root.rglob("*") if p.is_file()]},
            "metrics": {"candidate_count": len(candidates), "template_count": len(prepared["templates"]),
                        "clash_rejected": rejected, "engine_version": inputs.CODE_VERSION},
            "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl"},
        })
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a") as handle:
            handle.write(traceback.format_exc())
        finish_job(job.run_dir, False, {"error": str(exc)})
    return job.run_dir
