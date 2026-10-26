"""Paper-scaffold retargeting with native Rosetta interface design."""
from __future__ import annotations

from pathlib import Path
import secrets
import shutil
import subprocess
import traceback

from Bio.PDB import PDBIO, Select

from mn_molecular_design_studio import peptide_trap_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import STAGE_GENERATION_BACKBONE_SEQUENCE, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, update_status, write_json


def stage_scaffold(source: Path, destination: Path) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.workflows.kejia_idr import stage_template
    # The author's roles are peptide A, binder B. Retain that exact chain order.
    stage_template(source, destination, "A", "B")


class _PeptideOnly(Select):
    def accept_chain(self, chain):
        return chain.id == "A"


def normalize_output(output: Path, register: dict, target: str) -> list[dict]:
    rows = []
    for index, record in enumerate(read_json(output / "designs.json").get("records", []), 1):
        path = output / record["pdb"]
        if not path.resolve().is_relative_to(output.resolve()) or not path.is_file() or Path(record["pdb"]).name != record["pdb"]:
            raise ValueError("A peptide-trap result references an invalid structure path.")
        structure = read_structure(path)
        chains = protein_chains(structure)
        expected = target[register["start"] - 1:register["end"]]
        if set(chains) != {"A", "B"} or sequence(chains["A"]) != expected or record["target_sequence"] != expected:
            raise ValueError("The peptide-trap output did not preserve the requested target fragment.")
        if sequence(chains["B"]) != record["binder_sequence"]:
            raise ValueError("The designed binder sequence does not match its output structure.")
        if record["accepted"] is not True:
            continue
        peptide = path.with_suffix(".target.pdb")
        writer = PDBIO(); writer.set_structure(structure); writer.save(str(peptide), _PeptideOnly())
        rows.append({
            "candidate_id": f"peptide_trap_{output.name}_{index:04d}",
            "source_tool": "peptide_trap", "stage": STAGE_GENERATION_BACKBONE_SEQUENCE,
            "complex_pdb": str(path), "target_pdb": str(peptide),
            "binder_chains": ["B"], "target_chains": ["A"],
            "binder_sequence": record["binder_sequence"], "binder_length": str(len(chains["B"])),
            "metrics": dict(record["metrics"]) | {"sequence_design_required": False, "result_kind": "generated"},
            "raw_metadata": {"engine_version": inputs.CODE_VERSION, "source_doi": "10.5281/zenodo.10391229",
                             "scaffold": register["scaffold"], "burial_pattern": register["burial_pattern"],
                             "template_model_file": inputs.SCAFFOLDS[register["scaffold"]]["file"],
                             "target_fragment_range": [register["start"], register["end"]],
                             "threaded_target_sequence": expected, "requested_target_sequence": target,
                             "sequence_design_required": False,
                             "protocol": "Paper interface FastDesign/FastRelax adaptation; binding is unvalidated."},
        })
    return rows


def run_peptide_trap_design(
    target_sequence: str, scaffold_ids: list[str] | None = None,
    max_registers: int = 20, designs_per_register: int = 1,
    design_repeats: int = 2, relax_repeats: int = 2, max_backbone_rmsd: float = 2.0,
    seed: int = 0, cpu_cores: int = 4, job_name: str = "peptide-trap",
) -> Path:
    parameters = dict(locals())
    parameters["scaffold_ids"] = scaffold_ids if scaffold_ids is not None else list(inputs.SCAFFOLDS)
    prepared = inputs.prepare_request(parameters)
    seed = seed or secrets.randbelow(2147483646) + 1
    parameters.update(target_sequence=prepared["target_sequence"], seed=seed,
                      engine_version=inputs.CODE_VERSION, image=inputs.IMAGE, gpu_device="none", generator_only=True)
    job = create_job("design", "peptide_trap_design", "peptide_trap", {"scaffold_ids": parameters["scaffold_ids"]}, parameters)
    root = job.run_dir / "artifacts/raw/peptide_trap"
    staged = root / "inputs"
    staged.mkdir(parents=True, exist_ok=True)
    candidates, commands, attempted = [], [], 0
    try:
        native_protocol = staged / "design_after_seqthread.xml"
        shutil.copy2(inputs.library_root() / "source/design_after_seqthread.xml", native_protocol)
        write_json(root / "registers.json", {"matching_register_count": len(prepared["all_matches"]), "selected_registers": prepared["registers"]})
        for index, register in enumerate(prepared["registers"], 1):
            identifier = f"{register['scaffold'].replace('.', '_')}_{register['start']:05d}_{index:04d}"
            update_status(job.run_dir, "running", current_phase=f"Designing peptide-trap register {index}/{len(prepared['registers'])}")
            template = staged / f"{identifier}.pdb"
            stage_scaffold(inputs.template_path(register["scaffold"]), template)
            peptide = sequence(protein_chains(read_structure(template))["A"])
            request = {key: prepared[key] for key in ("designs_per_register", "design_repeats", "relax_repeats", "max_backbone_rmsd")}
            request.update(sequence=register["sequence"], template_peptide_sequence=peptide,
                           template="/work/" + template.relative_to(job.run_dir).as_posix(),
                           native_protocol="/work/" + native_protocol.relative_to(job.run_dir).as_posix(),
                           seed=(seed + index - 1) % 2147483647 or 1)
            request_path = staged / f"{identifier}.json"
            write_json(request_path, request)
            output = root / identifier
            output.mkdir()
            command = ["docker", "run", "--rm", "--network", "none", "--name", f"mn-trap-{job.run_id[-8:]}-{index}",
                       "--cpus", str(cpu_cores), "--mount", f"type=bind,src={job.run_dir.resolve()},dst=/work",
                       inputs.IMAGE, "--request", "/work/" + request_path.relative_to(job.run_dir).as_posix(),
                       "--output", "/work/" + output.relative_to(job.run_dir).as_posix()]
            commands.append({"register": register, "command": command})
            write_json(job.run_dir / "command.json", {"mode": "docker", "image": inputs.IMAGE, "steps": commands})
            with (job.run_dir / "stdout.log").open("a") as stdout, (job.run_dir / "stderr.log").open("a") as stderr:
                process = subprocess.run(command, stdout=stdout, stderr=stderr, cwd=job.run_dir, check=False)
            if process.returncode or not (output / "designs.json").is_file():
                raise RuntimeError(f"Peptide-trap interface design failed (exit {process.returncode}); see the job logs.")
            attempted += len(read_json(output / "designs.json").get("records", []))
            candidates.extend(normalize_output(output, register, prepared["target_sequence"]))
        normalized = write_candidates(job.run_dir, "peptide_trap", candidates)
        finish_job(job.run_dir, True, {
            "outputs": {"candidates": normalized, "artifacts": [{"path": str(p), "name": p.name, "type": "peptide_trap_output"} for p in root.rglob("*") if p.is_file()]},
            "metrics": {"candidate_count": len(candidates), "attempted_design_count": attempted,
                        "geometry_rejected": attempted - len(candidates), "matching_register_count": len(prepared["all_matches"]),
                        "selected_register_count": len(prepared["registers"]), "engine_version": inputs.CODE_VERSION},
            "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl"},
        })
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a") as handle:
            handle.write(traceback.format_exc())
        finish_job(job.run_dir, False, {"error": str(exc)})
    return job.run_dir
