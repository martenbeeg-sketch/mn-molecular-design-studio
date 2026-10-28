"""Independent inverse-folding jobs with a shared, non-destructive graph contract."""
from __future__ import annotations

import copy
import csv
import json
import secrets
import shutil
from pathlib import Path

from mn_molecular_design_studio import inverse_folding_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.runtime import reference_root

NATIVE = Path(__file__).with_name("inverse_folding_native.py")


def stage_structure(source, destination, *, protein_only=True):
    """Map author IDs to contiguous native IDs; keep the mapping for outputs."""
    from Bio.PDB import PDBIO, Select
    from Bio.SeqUtils import seq1
    from mn_molecular_design_studio.design_validation import read_structure
    model = next(read_structure(Path(source)).get_models())
    info = inputs.inventory(source)
    protein = {k for k, v in info["residues"].items() if v["kind"] == "protein"}
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    names = [c.id for c in model]
    if len(names) > len(alphabet):
        raise ValueError("This inverse-folding adapter supports at most 62 structure chains.")
    used = {c for c in names if c in alphabet}
    available = iter(c for c in alphabet if c not in used)
    chain_map = {c: c if c in alphabet else next(available) for c in names}
    rows = []
    chains = list(model)
    for chain in chains:
        model.detach_child(chain.id)
    for chain in chains:
        original = chain.id
        chain.id = chain_map[original]
        residue_objects = list(chain)
        for residue in residue_objects:
            chain.detach_child(residue.id)
        position = 0
        other = 0
        for residue in residue_objects:
            original_id = list(residue.id)
            label = f"{original}{residue.id[1]}{residue.id[2].strip()}"
            if label in protein:
                if not all(atom in residue for atom in ("N", "CA", "C")):
                    raise ValueError(f"Residue {label} needs N, CA and C backbone atoms.")
                aa = seq1(residue.resname, custom_map={"MSE": "M", "UNK": "X"}, undef_code="X")
                position += 1
                rows.append(dict(label=label, chain=original, original_id=original_id, aa=aa,
                                 engine_chain=chain.id, engine_label=f"{chain.id}{position}", position=position))
                if aa == "X":
                    # Backbone-only generator outputs often use UNK. Identity is
                    # a placeholder only at positions which will be redesigned.
                    residue.resname = "GLY"
                if residue.resname == "MSE":
                    residue.resname = "MET"
                residue.id = (" ", position, " ")
            elif protein_only or residue.id[0] == "W":
                continue
            else:
                other += 1
                residue.id = (residue.id[0], other, " ")
            for atom in residue.get_atoms():
                if atom.occupancy is None:
                    atom.occupancy = 1.
            chain.add(residue)
        if len(chain):
            model.add(chain)
    if not rows:
        raise ValueError("The structure has no protein backbone to design.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(str(destination))
    return rows, chain_map


def restore_backbone_sequence(destination, rows, reference_rows):
    """Label a topology-matched glycine/UNK backbone with its primary sequence."""
    if [(r["label"], r["chain"]) for r in rows] != [(r["label"], r["chain"]) for r in reference_rows]:
        raise ValueError("Provided conformers must match every primary chain and residue ID.")
    if [r["aa"] for r in rows] == [r["aa"] for r in reference_rows]:
        return
    if not rows or any(r["aa"] not in {"G", "X"} for r in rows):
        raise ValueError("Provided conformers must match every primary chain, residue ID and sequence.")
    from Bio.PDB import PDBIO
    from Bio.SeqUtils import seq3
    from mn_molecular_design_studio.design_validation import read_structure

    structure = read_structure(Path(destination))
    model = next(structure.get_models())
    for observed, reference in zip(rows, reference_rows):
        residue = model[observed["engine_chain"]][(" ", observed["position"], " ")]
        residue.resname = seq3(reference["aa"]).upper()
        observed["aa"] = reference["aa"]
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(str(destination))


def design_positions(source, candidate, p):
    info = inputs.inventory(source)
    proteins = [c for c, v in info["chains"].items() if v["kind"] == "protein"]
    requested = inputs.chain_list(p["design_chains"])
    inherited = list(candidate.get("binder_chains") or [])
    targets = set(candidate.get("target_chains") or [])
    chains = requested or inherited or [c for c in proteins if c not in targets]
    if not chains or not set(chains) <= set(proteins):
        raise ValueError("Select protein chains present in this input structure.")
    if targets & set(chains):
        raise ValueError("A connected candidate's fixed target chains cannot be redesigned.")
    raw = candidate.get("raw_metadata") or {}
    if raw.get("fixed_motif_required") and not (raw.get("motif_mapping_verified") and raw.get("fixed_residues")):
        raise ValueError("Sequence design needs verified fixed motif residue mappings.")
    fixed = set(inputs.residues(p["fixed_residues"], info)) | set(raw.get("fixed_residues") or [])
    if not fixed <= set(info["residues"]):
        raise ValueError("Fixed motif residues must exist in the source structure.")
    selected = set(inputs.residues(p["designed_residues"], info, chains=chains))
    eligible = {r["label"] for c in chains for r in info["chains"][c]["residues"]}
    mutable = (selected or eligible) - fixed
    if not mutable and p["task"] == "design":
        raise ValueError("No residues remain to design after applying fixed positions.")
    return chains, mutable


def map_constraints(path, source, candidate_id, mapping, mutable):
    """Translate Caliby's author-position CSV to the staged native numbering."""
    if not path:
        return {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "pdb_key" not in (reader.fieldnames or []):
            raise ValueError("Position constraints CSV needs a pdb_key column.")
        matching = [row for row in reader if row["pdb_key"] in {Path(source).stem, candidate_id, "*"}]
    if len(matching) > 1:
        raise ValueError("Use one position-constraints row per structure.")
    if not matching:
        return {}
    info = inputs.inventory(source)
    labels = {r["label"]: r["engine_label"] for r in mapping}
    result = {}
    for key, value in matching[0].items():
        if not value or key == "pdb_key":
            continue
        if key in {"fixed_pos_seq", "fixed_pos_scn"}:
            result[key] = ",".join(labels[v] for v in inputs.residues(value, info))
        elif key in {"fixed_pos_override_seq", "pos_restrict_aatype"}:
            entries = []
            for item in value.split(","):
                position, aa = item.strip().split(":", 1)
                if position not in labels or position not in mutable:
                    raise ValueError("Sequence overrides and amino-acid restrictions must refer to designable positions.")
                if not aa or any(v not in inputs.AA for v in aa) or key == "fixed_pos_override_seq" and len(aa) != 1:
                    raise ValueError("Invalid amino-acid identity in constraints CSV.")
                entries.append(labels[position] + ":" + aa)
            result[key] = ",".join(entries)
        elif key == "symmetry_pos":
            groups = []
            for group in value.split("|"):
                positions = inputs.residues(group, info)
                if any(v not in mutable for v in positions):
                    raise ValueError("Symmetry constraints may tie only designable residues.")
                groups.append(",".join(labels[v] for v in positions))
            result[key] = "|".join(groups)
        else:
            raise ValueError(f"Unsupported Caliby constraint column: {key}")
    return result


def source_rows(p):
    paths = p.get("structure_paths") or []
    if paths:
        return [dict(candidate_id=f"input_{i:04d}", complex_pdb=str(Path(v).resolve()), raw_metadata={}) for i, v in enumerate(paths)]
    selected = p.get("selected_candidate_ids")
    rows = read_candidates(Path(p["candidates_jsonl"]))
    if selected is not None:
        if not set(selected) <= {r["candidate_id"] for r in rows}:
            raise ValueError("A selected candidate is absent from the source job.")
        rows = [r for r in rows if r["candidate_id"] in selected]
    if not rows:
        raise ValueError("Select at least one structure candidate.")
    for row in rows:
        for key in ("complex_pdb", "binder_pdb", "target_pdb"):
            if row.get(key):
                resolved = resolve_stored_path(row[key], run_dir=Path(p["source_run_dir"]), must_exist=True)
                if resolved is None:
                    raise ValueError(f"The selected candidate's {key} is unavailable.")
                row[key] = str(resolved)
    return rows


def _cp(run, path):
    return "/work/" + str(Path(path).relative_to(run))


def prepare_items(engine, p, run):
    items = []
    for index, candidate in enumerate(source_rows(p)):
        source = Path(candidate.get("complex_pdb") or candidate.get("binder_pdb") or "")
        if not source.is_file():
            raise ValueError("Each candidate must contain a readable structure.")
        directory = run / "artifacts/inverse_folding/inputs" / f"input_{index:04d}"
        directory.mkdir(parents=True, exist_ok=True)
        original = directory / ("original" + source.suffix)
        shutil.copy2(source, original)
        staged = directory / "backbone.pdb"
        mapping, chains = stage_structure(source, staged, protein_only=engine != "boltzgen_inverse")
        designed, mutable = design_positions(source, candidate, p)
        if any(r["aa"] == "X" and r["label"] not in mutable for r in mapping):
            raise ValueError("Unknown amino-acid identities can only occur at redesigned positions.")
        constraints = map_constraints(p.get("constraints_csv"), source, candidate["candidate_id"], mapping, mutable)
        fixed_csv = set(inputs.residues(constraints.get("fixed_pos_seq", "")))
        mutable -= {r["label"] for r in mapping if r["engine_label"] in fixed_csv}
        if not mutable and p["task"] == "design":
            raise ValueError("No residues remain to design after CSV constraints.")
        entry = dict(key=f"input_{index:04d}", path=_cp(run, staged), original=_cp(run, original),
            mapping=mapping, chain_map=chains, design_chains=designed,
            mutable=sorted(mutable), candidate=candidate, constraints=constraints, ensemble_paths=[])
        for j, conformer in enumerate(p.get("ensemble_paths") or []):
            path = directory / f"conformer_{j:04d}.pdb"
            mapped, _ = stage_structure(conformer, path)
            restore_backbone_sequence(path, mapped, mapping)
            entry["ensemble_paths"].append(_cp(run, path))
        items.append(entry)
    return items


def designed_structure(original, destination, mapping, chain_sequences, mutable):
    """Keep coordinates and fixed context; remove obsolete mutated sidechains."""
    from Bio.PDB import MMCIFIO
    from Bio.SeqUtils import seq3
    from mn_molecular_design_studio.design_validation import read_structure
    structure = read_structure(original)
    model = next(structure.get_models())
    expected = {}
    for row in mapping:
        expected.setdefault(row["engine_chain"], []).append(row)
    if set(chain_sequences) != set(expected):
        raise ValueError("The model returned unexpected protein chain IDs.")
    for chain, rows in expected.items():
        sequence = chain_sequences[chain]
        if len(sequence) != len(rows) or any(v not in inputs.AA for v in sequence):
            raise ValueError("The model output must have one standard amino acid per input residue.")
        for row, aa in zip(rows, sequence):
            residue = model[row["chain"]][tuple(row["original_id"])]
            if row["label"] not in mutable:
                if aa != row["aa"]:
                    raise ValueError(f"The model changed fixed residue {row['label']}.")
                continue
            residue.resname = seq3(aa).upper()
            for atom in list(residue):
                if atom.id not in {"N", "CA", "C", "O", "OXT"}:
                    residue.detach_child(atom.id)
    writer = MMCIFIO()
    writer.set_structure(model)
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer.save(str(destination))
    return {r["chain"]: chain_sequences[r["engine_chain"]] for r in mapping}


def normalize_outputs(engine, p, run, items, *, allow_partial=False):
    records = json.loads((run / "artifacts/inverse_folding/native_results.json").read_text())
    by_key = {v["key"]: v for v in items}
    outputs, table, fasta = [], [], []
    counts = {key: 0 for key in by_key}
    for record in records:
        item = by_key[record["key"]]
        counts[item["key"]] += 1
        index = counts[item["key"]]
        cid = item["key"] + f"_seq_{index:04d}"
        path = run / "artifacts/inverse_folding/structures" / (cid + ".cif")
        original = run / item["original"].removeprefix("/work/")
        mutable = set(item["mutable"]) if p["task"] == "design" else set()
        sequences = designed_structure(original, path, item["mapping"], record["sequences"], mutable)
        designed = item["design_chains"]
        sequence = ":".join(sequences[c] for c in designed)
        parent = item["candidate"]
        metrics = {k: v for k, v in record.get("metrics", {}).items() if isinstance(v, (int, float))}
        # Prediction confidence metrics do not describe the redesigned sequence.
        raw = copy.deepcopy(parent.get("raw_metadata") or {})
        raw.update(operation="scoring" if p["task"] == "score" else "inverse_folding", structure_kind="sequence_on_input_backbone",
                   source_candidate_id=parent["candidate_id"], source_tool=parent.get("source_tool"), engine_version=inputs.ENGINES[engine]["version"],
                   fixed_residues=[r["label"] for r in item["mapping"] if r["label"] not in mutable],
                   native_metrics=record.get("metrics", {}), original_structure=str(original))
        outputs.append(dict(candidate_id=cid, stage="sequence_design", source_tool=engine, complex_pdb=str(path),
            binder_sequence=sequence, binder_chains=designed,
            target_chains=[c for c in sequences if c not in designed], target_pdb=parent.get("target_pdb"),
            hotspots=parent.get("hotspots") or [], parents=[parent["candidate_id"]], metrics=metrics, raw_metadata=raw))
        table.append(dict(candidate_id=cid, source_candidate_id=parent["candidate_id"], sequence=sequence, **metrics))
        fasta += [f">{cid}\n{sequence}\n"]
    required = 1 if p["task"] == "score" else p["num_sequences"]
    if any(n > required or not allow_partial and n != required for n in counts.values()):
        raise ValueError("The engine did not produce the requested number of sequences for every structure.")
    root = run / "artifacts/inverse_folding"
    (root / "sequences.fasta").write_text("".join(fasta))
    columns = list(dict.fromkeys(k for row in table for k in row))
    with (root / "sequences.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(table)
    return write_candidates(run, engine, outputs)


def recover(run):
    """Normalize durable native sequence receipts; never initialize a model."""
    root = Path(run) / 'artifacts/inverse_folding'
    if not (root / 'native_request.json').is_file() or not (root / 'native_results.json').is_file():
        return []
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    request = read_json(root / 'native_request.json')
    return normalize_outputs(request['engine'], request['parameters'], Path(run), request['items'], allow_partial=True)


def _run(engine, parameters, existing_job=None):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.sequence_design import _run_shell_steps
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import _json_value
    p = inputs.prepare(engine, parameters, require_input=True)
    ready = inputs.readiness(engine, p)
    if not ready["ready"]:
        raise ValueError("Missing reference weights: " + ", ".join(ready["missing"]))
    p["seed"] = p["seed"] or secrets.randbelow(2147483646) + 1
    job_inputs = _json_value({"structure_paths": p["structure_paths"],
        "source_run_dir": p.get("source_run_dir"), "candidates_jsonl": p.get("candidates_jsonl")})
    job_params = _json_value(p)
    job = existing_job or create_job("design", "sequence_design", engine, inputs=job_inputs, params=job_params)
    run = job.run_dir
    update_status(run, "preparing", current_phase="Preparing structure-to-sequence inputs", campaign_name=p.get("job_name"))
    try:
        items = prepare_items(engine, p, run)
        root = run / "artifacts/inverse_folding"
        if p.get("constraints_csv"):
            shutil.copy2(p["constraints_csv"], root / "position_constraints.csv")
        runner = root / "runner.py"
        shutil.copy2(NATIVE, runner)
        payload = dict(engine=engine, parameters=p, items=items, output=_cp(run, root / "native_results.json"),
                       checkpoint="/ref/" + str(inputs.weights(engine, p)[0].relative_to(reference_root())))
        # Only job-owned paths enter the container request; external input paths
        # stay in the host-side provenance and were copied above.
        for key in ("structure_path", "structure_paths", "source_run_dir", "candidates_jsonl", "constraints_csv", "ensemble_paths"):
            payload["parameters"].pop(key, None)
        request = root / "native_request.json"
        write_json(request, payload)
        cache = root / "cache"
        cache.mkdir(exist_ok=True)
        command = ["docker", "run", "--rm", *docker_gpu_args(p.get("gpu_device", "all")), "--shm-size=4G",
            "--entrypoint", "python", "-v", f"{run}:/work", "-v", f"{reference_root()}:/ref:ro", "-w", "/work",
            "-e", "HF_HUB_OFFLINE=1", "-e", "MODEL_PARAMS_DIR=/ref/caliby", "-e", "PDB_MIRROR_PATH=", "-e", "CCD_MIRROR_PATH=",
            "-e", "TORCH_HOME=/work/artifacts/inverse_folding/cache/torch", "-e", "XDG_CACHE_HOME=/work/artifacts/inverse_folding/cache",
            "-e", "WANDB_MODE=disabled", "-e", "BOLTZGEN_OPT=off", "-e", "OMP_NUM_THREADS=4",
            inputs.ENGINES[engine]["image"], _cp(run, runner), _cp(run, request)]
        rc = _run_shell_steps(run, [dict(name=inputs.ENGINES[engine]["title"], command=command,
            candidate_ids=[item['candidate']['candidate_id'] for item in items], native_invocations=1)])
        candidates = normalize_outputs(engine, p, run, items, allow_partial=rc != 0) if (root / 'native_results.json').is_file() else []
        from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import Artifact
        artifacts = [Artifact(path.stem, path, path.suffix.lstrip(".")).to_json(run)
                     for path in sorted(root.rglob("*")) if path.is_file() and "cache" not in path.relative_to(root).parts]
        finish_job(run, rc == 0 and bool(candidates), dict(exit_code=rc,
            outputs=dict(artifacts=artifacts, candidates=candidates, sequences_fasta=str(root / "sequences.fasta"), sequences_csv=str(root / "sequences.csv")),
            metrics=dict(sequence_count=len(candidates)), downstream_artifacts=dict(candidates_jsonl=str(run / "artifacts/normalized_candidates/candidates.jsonl"))))
    except Exception as exc:
        finish_job(run, False, dict(error=str(exc)))
        raise
    return run


def run_esm_if1(structure_path: Path | None = None, structure_paths: list[Path] | None = None,
    source_run_dir: Path | None = None, candidates_jsonl: Path | None = None, selected_candidate_ids: list[str] | None = None,
    design_chains: str = "", designed_residues: str = "", fixed_residues: str = "", num_sequences: int = 10,
    temperature: float = 1., multichain_backbone: bool = True, seed: int = 0, job_name: str = "esm-if1-design",
    gpu_device: object = "all", existing_job=None) -> Path:
    return _run("esm_if1", {k: v for k, v in locals().items() if k != "existing_job"}, existing_job)


def run_caliby(structure_path: Path | None = None, structure_paths: list[Path] | None = None,
    source_run_dir: Path | None = None, candidates_jsonl: Path | None = None, selected_candidate_ids: list[str] | None = None,
    design_chains: str = "", designed_residues: str = "", fixed_residues: str = "", num_sequences: int = 4,
    task: str = "design", model_type: str = "caliby", ensemble_mode: str = "generate", num_conformers: int = 32,
    ensemble_paths: list[Path] | None = None, constraints_csv: Path | None = None, omit_aas: str = "",
    temperature: float = .01, num_sweeps: int = 500, batch_size: int = 1, optimization_mode: str = "off",
    seed: int = 0, job_name: str = "caliby-design", gpu_device: object = "all", existing_job=None) -> Path:
    return _run("caliby", {k: v for k, v in locals().items() if k != "existing_job"}, existing_job)


def run_boltzgen_inverse(structure_path: Path | None = None, structure_paths: list[Path] | None = None,
    source_run_dir: Path | None = None, candidates_jsonl: Path | None = None, selected_candidate_ids: list[str] | None = None,
    design_chains: str = "", designed_residues: str = "", fixed_residues: str = "", num_sequences: int = 10,
    temperature: float = .1, omit_aas: str = "", seed: int = 0, job_name: str = "boltzgen-inverse-design",
    gpu_device: object = "all", existing_job=None) -> Path:
    return _run("boltzgen_inverse", {k: v for k, v in locals().items() if k != "existing_job"}, existing_job)
