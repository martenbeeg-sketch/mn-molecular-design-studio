"""One sequence-design entry point for structures, saved candidates and graphs."""
from __future__ import annotations

import json
from pathlib import Path

from mn_molecular_design_studio import mpnn_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import JobPaths, create_job, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates

STAGES = ["generation.backbone", "generation.backbone_sequence", "sequence_design", "complex_refolding"]


def verification_request(run_dir, gpu_device="0"):
    """Fold each designed chain separately, including paired antibody outputs."""
    from Bio.SeqUtils import seq1
    from mn_molecular_design_studio.design_validation import read_structure
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    tasks = []
    for candidate in read_candidates(run_dir / "artifacts/normalized_candidates/candidates.jsonl"):
        path = resolve_stored_path(candidate.get("complex_pdb") or candidate.get("binder_pdb"), run_dir=run_dir, must_exist=True)
        if not path:
            raise ValueError("Verification needs a saved sequence-design structure.")
        model = next(read_structure(path).get_models())
        for chain in candidate.get("binder_chains") or []:
            sequence = "".join(seq1(r.resname, undef_code="X") for r in model[chain] if "CA" in r)
            if sequence:
                tasks.append({"name": candidate["candidate_id"] + "-" + chain, "protein_sequences": sequence})
    if not tasks:
        raise ValueError("There are no designed protein chains to verify.")
    return {"schema_version": 1, "workflow": "refolding.optimized_prediction", "parameters": {
        "engine": "esmfold", "job_name": run_dir.name + "-verification", "tasks": tasks, "samples": 1},
        "resources": {"gpu_device": gpu_device}}


def _pdb_copy(source, destination):
    """Convert CIF author IDs without truncating colliding chain names."""
    from Bio.PDB import PDBIO
    from mn_molecular_design_studio.design_validation import read_structure
    structure = read_structure(Path(source))
    model = next(structure.get_models())
    names = list(dict.fromkeys(c.id for c in model))
    used = {c for c in names if len(c) == 1}
    available = iter(c for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" if c not in used)
    mapping = {name: name if len(name) == 1 else next(available, None) for name in names}
    if None in mapping.values():
        raise ValueError("This MPNN backend supports at most 62 structure chains.")
    # Detach first: swaps and shortened IDs must not collide in Bio.PDB's map.
    chains = list(model)
    for chain in chains:
        model.detach_child(chain.id)
    for chain in chains:
        chain.id = mapping[chain.id]
        for atom in chain.get_atoms():
            if atom.occupancy is None:
                atom.occupancy = 1.
        model.add(chain)
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(str(destination))
    return mapping


def residue_map(source_path, role_metadata):
    """The legacy adapter renumbers residues; retain a full author-ID mapping."""
    from mn_molecular_design_studio.design_validation import read_structure
    mapping = {}
    role_rows = (role_metadata.get("chain_roles") or {}).get("chain_map") or {}
    chains = {r["original_chain"]: r["engine_chain"] for section in ("binder", "targets") for r in role_rows.get(section, [])}
    model = next(read_structure(Path(source_path)).get_models())
    for chain in model:
        if chain.id in chains:
            for index, residue in enumerate(chain, 1):
                label = f"{chain.id}{residue.id[1]}{residue.id[2].strip()}"
                mapping[label] = f"{chains[chain.id]}{index}"
    return mapping


def candidate_options(options, source, design_chains, role_metadata=None, *, normalized_source=None):
    """Validate MPNN options against the staged, engine-safe structure.

    Generator backbones can contain CA-only ``UNK`` residues or long CIF chain
    IDs.  The staged structure repairs those representations and renumbers
    chains/residues for LigandMPNN, so validate selectors there after mapping
    any user-entered source residue labels into staged labels.
    """
    source = Path(source)
    structure = Path(normalized_source) if normalized_source else source
    source_chains = inputs.chain_list(design_chains)
    mapping = residue_map(source, role_metadata) if role_metadata else {}
    engine_chains = inputs.chain_list(
        (role_metadata or {}).get("binder_chains") or source_chains
    )

    def remap_labels(value, *, chains=None):
        labels = inputs.residues(value, info=None, chains=chains)
        try:
            return [mapping[label] for label in labels] if role_metadata else labels
        except KeyError as exc:
            raise ValueError(f"Residue {exc.args[0]} could not be mapped to the MPNN input.") from exc

    translated = dict(options)
    for key in ("designed_residues", "fixed_residues"):
        translated[key] = " ".join(remap_labels(translated.get(key), chains=source_chains))
    for key in ("transmembrane_buried", "transmembrane_interface"):
        translated[key] = " ".join(remap_labels(translated.get(key)))
    for key in ("bias_aa_per_residue", "omit_aa_per_residue"):
        raw = inputs._json_map(translated.get(key), key.replace("_", " "))
        translated[key] = {
            remap_labels(label, chains=source_chains)[0]: value
            for label, value in raw.items()
        }

    p = inputs.prepare(translated | {
        "structure_path": str(structure),
        "source_run_dir": None,
        "candidates_jsonl": None,
        "design_chains": ",".join(engine_chains),
    })
    info = inputs.inventory(structure)
    fixed = inputs.residues(p["fixed_residues"], info, chains=engine_chains)
    designed = inputs.residues(p["designed_residues"], info, chains=engine_chains)
    if designed:
        fixed += [r["label"] for c in engine_chains for r in info["chains"][c]["residues"] if r["label"] not in designed]
    p["fixed_residues"] = list(dict.fromkeys(fixed))
    p["symmetry_groups"] = []
    if p["homo_oligomer"] > 1:
        rows = [info["chains"][c]["residues"] for c in engine_chains]
        p["symmetry_groups"] = [[chain_rows[i]["label"] for chain_rows in rows] for i in range(len(rows[0]))]
    return p


def shared_arguments(p, directory, container_dir):
    args = []
    if p["bias_aa"]:
        args += ["--bias_AA", p["bias_aa"]]
    for key, flag in [("bias_aa_per_residue", "--bias_AA_per_residue"), ("omit_aa_per_residue", "--omit_AA_per_residue")]:
        if p[key]:
            name = key + ".json"
            write_json(directory / name, p[key])
            args += [flag, container_dir + "/" + name]
    if p["symmetry_groups"]:
        groups = p["symmetry_groups"]
        args += ["--symmetry_residues", "|".join(",".join(group) for group in groups),
                 "--symmetry_weights", "|".join(",".join(str(1 / len(group)) for _ in group) for group in groups)]
    if p["model_type"] == "global_label_membrane_mpnn":
        args += ["--global_transmembrane_label", str(int(p["global_transmembrane_label"]))]
    if p["model_type"] == "per_residue_label_membrane_mpnn":
        args += ["--transmembrane_buried", p["transmembrane_buried"], "--transmembrane_interface", p["transmembrane_interface"]]
    if p["use_side_chain_context"]:
        args += ["--ligand_mpnn_use_side_chain_context", "1"]
    if p["pack_side_chains"]:
        args += ["--pack_side_chains", "1", "--number_of_packs_per_design", str(p["number_of_packs_per_design"])]
    return args


def run_mpnn_sequence_design(
    source_run_dir: Path | None = None, candidates_jsonl: Path | None = None,
    structure_path: Path | None = None, model_type: str = "protein_mpnn", noise_level: float | None = None,
    design_chains: str = "", designed_residues: str = "", fixed_residues: str = "",
    fix_interface_residues: bool = False,
    num_seq_per_target: int = 2, sampling_temp: float = .1, omit_aas: str = "C", seed: int | None = None,
    bias_aa: str = "", bias_aa_per_residue: dict | None = None, omit_aa_per_residue: dict | None = None,
    homo_oligomer: int = 0, global_transmembrane_label: bool = True,
    transmembrane_buried: str = "", transmembrane_interface: str = "",
    use_side_chain_context: bool = False, pack_side_chains: bool = False, number_of_packs_per_design: int = 1,
    backend: str = "auto", selected_candidate_ids: list[str] | None = None,
    verify_sequences: bool = False, job_name: str = "proteinmpnn-design",
    require_backbone_hotspot_filter_pass: bool = False, gpu_device: object = "0",
    existing_job: JobPaths | None = None,
) -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design as native
    p = inputs.prepare({k: v for k, v in locals().items() if k not in {"existing_job", "native"}}, require_input=True)
    weights = inputs.checkpoint(model_type, noise_level)
    if not weights.is_file():
        raise ValueError(f"Missing {inputs.MODELS[model_type]['label']} checkpoint: {weights}")
    if structure_path:
        rows = []
    else:
        rows = [r for r in read_candidates(Path(candidates_jsonl)) if r.get("stage") in STAGES
                and (selected_candidate_ids is None or r.get("candidate_id") in selected_candidate_ids)
                and (not require_backbone_hotspot_filter_pass or r.get("metrics", {}).get("passes_backbone_hotspot_filter") is not False)]
        if not rows:
            raise ValueError("Select one or more available input designs.")
    job = existing_job or create_job("design", job_type="sequence_design", tool="ligandmpnn",
                                     inputs={k: str(v) for k, v in {"structure_path": structure_path, "source_run_dir": source_run_dir, "candidates_jsonl": candidates_jsonl}.items() if v}, params=p)
    if structure_path:
        chains = inputs.chain_list(design_chains)
        info = inputs.inventory(structure_path)
        rows = [{"candidate_id": "input", "stage": "generation.backbone", "source_tool": "structure_input",
                 "complex_pdb": str(Path(structure_path).resolve()), "binder_chains": chains,
                 "target_chains": [c for c in info["chains"] if c not in chains], "raw_metadata": {}}]
    source = job.run_dir / "artifacts/mpnn_source"
    source.mkdir(parents=True, exist_ok=True)
    # Save only the selected candidates; original jobs and structures remain immutable.
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    for row in rows:
        for key in ("target_pdb", "complex_pdb", "binder_pdb"):
            if row.get(key):
                path = resolve_stored_path(row[key], run_dir=Path(source_run_dir) if source_run_dir else None, must_exist=True)
                if path:
                    row[key] = str(path)
    candidate_file = source / "candidates.jsonl"
    can_foundry = model_type in {"protein_mpnn", "ligand_mpnn"} and not pack_side_chains
    all_cif = all(native._is_cif_path(Path(r.get("complex_pdb") or r.get("binder_pdb") or "")) for r in rows)
    chosen = "foundry" if backend == "foundry" or backend == "auto" and can_foundry and all_cif else "shared_ligandmpnn"
    if chosen == "foundry":
        from Bio.PDB import MMCIFIO
        from mn_molecular_design_studio.design_validation import read_structure
        for row in rows:
            path = Path(row.get("complex_pdb") or row.get("binder_pdb") or "")
            if not native._is_cif_path(path):
                converted = source / f"{row['candidate_id']}.cif"
                writer = MMCIFIO()
                writer.set_structure(read_structure(path))
                writer.save(str(converted))
                row["complex_pdb"] = str(converted)
    candidate_file.write_text("".join(json.dumps(r) + "\n" for r in rows))
    p["resolved_backend"] = chosen
    write_json(job.run_dir / "artifacts/mpnn_settings.json", p | {"checkpoint": str(weights)})
    if chosen == "foundry":
        result = native.run_foundry_mpnn_sequence_design(source, candidate_file, batch_size=num_seq_per_target,
            model_type=model_type, gpu_device=gpu_device, existing_job=job, direct_input=p)
    else:
        result = native.run_ligandmpnn_sequence_design(source, candidate_file, model_type=model_type,
            design_chains=design_chains, num_seq_per_target=num_seq_per_target, sampling_temp=sampling_temp,
            omit_aas=omit_aas, seed=p["seed"], accepted_stages=STAGES, gpu_device=gpu_device, existing_job=job, direct_input=p)
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    if verify_sequences and read_json(result / "result.json").get("success") is True:
        from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import submit_workflow_request
        from mn_molecular_design_studio.workspace import mark_studio_job
        request = result / "requests/verify_sequences.json"
        write_json(request, verification_request(result, gpu_device))
        verification = submit_workflow_request(request)
        marker = read_json(result / "molecular_design_studio.json")
        if marker.get("project_id") and marker.get("studio_owner_id"):
            mark_studio_job(verification, project_id=marker["project_id"], owner_id=marker["studio_owner_id"])
        metadata = read_json(result / "metadata.json")
        metadata["validation_run_dir"] = verification["job"]["run_dir"]
        write_json(result / "metadata.json", metadata)
    return result
