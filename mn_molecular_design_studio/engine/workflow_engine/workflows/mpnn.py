"""One sequence-design entry point for structures, saved candidates and graphs."""
from __future__ import annotations

import json
import copy
import secrets
import re
from pathlib import Path

from mn_molecular_design_studio import mpnn_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import JobPaths, create_job, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

STAGES = ["generation.backbone", "generation.backbone_sequence", "sequence_design", "complex_refolding"]


def inherited_fixed_residues(candidate, structure, binders):
    """Map inherited constraints by explicit chain order and residue position.

    Prediction may rename/renumber chains. The full source snapshot makes this
    mapping checkable; neither fuzzy sequence matching nor unrestricted fallback
    is allowed when a protected ancestor cannot be mapped.
    """
    from mn_molecular_design_studio import design_validation as dv
    from mn_molecular_design_studio.sabr_runtime import fixed_residues_from_annotation, CDR_NAMES
    current = dv.protein_chains(dv.read_structure(structure))
    current_rows = [(c, r) for c in binders for r in current[c]]
    current_labels = [f"{c}{r.id[1]}{r.id[2].strip()}" for c, r in current_rows]
    fixed = [label for label, (_, r) in zip(current_labels, current_rows) if r.resname == "CYS"]
    ancestor, seen = candidate, set()
    while isinstance(ancestor, dict) and ancestor:
        if id(ancestor) in seen or len(seen) >= 64:
            raise ValueError("Ambiguous or cyclic sequence-constraint lineage.")
        seen.add(id(ancestor))
        raw = ancestor.get("raw_metadata") or {}
        conditions = (raw.get("reference_constraints") or {}).get("condition_spec")
        if conditions and not (raw.get("fixed_motif_required") and raw.get("motif_mapping_verified")):
            import yaml
            spec = yaml.safe_load(conditions) if isinstance(conditions, str) else conditions
            def sequence_selectors(value, in_sequence=False):
                if isinstance(value, dict):
                    for key, item in value.items():
                        if in_sequence and key == "select":
                            yield item
                        else:
                            yield from sequence_selectors(item, in_sequence or key == "C_SEQ")
                elif isinstance(value, list):
                    for item in value:
                        yield from sequence_selectors(item, in_sequence)
                elif in_sequence and value:
                    yield None
            for selector in sequence_selectors(spec):
                match = re.fullmatch(r"sel\(['\"]([^'\"]+)['\"]\)", str(selector).strip())
                if not match or match[1].split('/')[0] not in set(ancestor.get("target_chains") or []):
                    raise ValueError("Native sequence conditions need a verified fixed-residue mapping before soluble redesign.")
        protected = list(raw.get("fixed_residues") or [])
        prompt, selection = raw.get("sampled_binder_prompt"), raw.get("sabr_cdr_selection")
        if protected or prompt or selection:
            source = resolve_stored_path(ancestor.get("complex_pdb") or ancestor.get("binder_pdb"), must_exist=True)
            if source is None:
                raise ValueError("Inherited fixed-residue constraints have no verifiable source structure.")
            proteins = dv.protein_chains(dv.read_structure(source))
            source_binders = list(ancestor.get("binder_chains") or [])
            if len(source_binders) != len(binders) or any(c not in proteins for c in source_binders) or [len(proteins[c]) for c in source_binders] != [len(current[c]) for c in binders]:
                raise ValueError("Inherited fixed-residue constraints have an ambiguous sequence/structure mapping.")
            rows = [(c, r) for c in source_binders for r in proteins[c]]
            labels = [f"{c}{r.id[1]}{r.id[2].strip()}" for c, r in rows]
            if prompt:
                if len(prompt) != len(labels):
                    raise ValueError("Native fixed framework prompt cannot be mapped to the candidate structure.")
                protected += [label for label, aa in zip(labels, prompt) if aa != "#"]
            if selection:
                antibody_fixed, _ = fixed_residues_from_annotation(source, source_binders,
                    selection.get("regions") or CDR_NAMES, selection["residues_by_region"])
                protected += antibody_fixed
            targets = set(ancestor.get("target_chains") or [])
            target_labels = {f"{c}{r.id[1]}{r.id[2].strip()}" for c in targets if c in proteins for r in proteins[c]}
            if set(protected) - set(labels) - target_labels:
                raise ValueError("Inherited fixed residues cannot be mapped to the sequence-design chains.")
            for index, (label, (_, residue)) in enumerate(zip(labels, rows)):
                if label in protected:
                    if residue.resname != current_rows[index][1].resname:
                        raise ValueError("Candidate sequence already differs from an inherited fixed-residue identity.")
                    fixed.append(current_labels[index])
        ancestor = raw.get("source_candidate")
    return list(dict.fromkeys(fixed))


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
    fix_interface_residues: bool = False, interface_cutoff: float = 4.0,
    sequence_stage: str = "redesign", optimization_enabled: bool = False,
    num_seq_per_target: int = 2, sampling_temp: float = .1, omit_aas: str = "C", seed: int | None = None,
    bias_aa: str = "", bias_aa_per_residue: dict | None = None, omit_aa_per_residue: dict | None = None,
    homo_oligomer: int = 0, global_transmembrane_label: bool = True,
    transmembrane_buried: str = "", transmembrane_interface: str = "",
    use_side_chain_context: bool = False, pack_side_chains: bool = False, number_of_packs_per_design: int = 1,
    backend: str = "auto", selected_candidate_ids: list[str] | None = None,
    verify_sequences: bool = False, job_name: str = "proteinmpnn-design",
    require_backbone_hotspot_filter_pass: bool = False, gpu_device: object = "0",
    auto_detect_cdrs: bool = False, cdr_regions: list[str] | None = None,
    sabr_chain_type: str = "auto", sabr_mode: str = "sabr", sabr_noise_level: float = 0.0,
    existing_job: JobPaths | None = None,
) -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design as native
    if sequence_stage == "initial":
        model_type = "protein_mpnn"
    elif sequence_stage == "soluble_optimization":
        model_type, fix_interface_residues, backend = "soluble_mpnn", True, "shared_ligandmpnn"
    p = inputs.prepare({k: v for k, v in locals().items() if k not in {"existing_job", "native"}}, require_input=True)
    p["seed"] = p["seed"] or secrets.randbelow(2147483646) + 1
    weights = inputs.checkpoint(model_type, noise_level)
    if not isinstance(auto_detect_cdrs, bool):
        raise ValueError("Automatic SAbR CDR selection must be enabled or disabled.")
    if structure_path:
        chains = inputs.chain_list(design_chains)
        info = inputs.inventory(Path(structure_path))
        rows = [{"candidate_id": "input", "stage": "generation.backbone", "source_tool": "structure_input",
                 "complex_pdb": str(Path(structure_path).resolve()), "binder_chains": chains,
                 "target_chains": [c for c in info["chains"] if c not in chains], "raw_metadata": {}}]
        if fix_interface_residues or sequence_stage == "soluble_optimization":
            from mn_molecular_design_studio.engine.workflow_engine.workflows.refolding import _sequences_by_chain
            sequences = _sequences_by_chain(Path(structure_path))
            rows[0].update(stage="generation.backbone_sequence", binder_sequence="".join(sequences.get(c, "") for c in chains))
    else:
        rows = [r for r in read_candidates(Path(candidates_jsonl)) if r.get("stage") in STAGES
                and (selected_candidate_ids is None or r.get("candidate_id") in selected_candidate_ids)
                and (not require_backbone_hotspot_filter_pass or r.get("metrics", {}).get("passes_backbone_hotspot_filter") is not False)]
    if not rows:
        raise ValueError("Select one or more available input designs.")
    if source_run_dir:
        from mn_molecular_design_studio.mimosa_motif import hydrate_candidate_fixed_motif

        rows = [hydrate_candidate_fixed_motif(row, source_run_dir) for row in rows]
    # Validate before partitioning: a co-designed sequence is never an implicit
    # request to replace it when its structure/roles are unavailable or ambiguous.
    for row in rows:
        for key in ("target_pdb", "complex_pdb", "binder_pdb"):
            if row.get(key):
                path = resolve_stored_path(row[key], run_dir=source_run_dir, must_exist=True)
                if path:
                    row[key] = str(path)
        if sequence_stage != "redesign" or fix_interface_residues:
            path = Path(row.get("complex_pdb") or row.get("binder_pdb") or "")
            binders = inputs.chain_list(row.get("binder_chains"))
            targets = inputs.chain_list(row.get("target_chains"))
            from mn_molecular_design_studio.design_validation import read_structure
            available = {c.id for c in next(read_structure(path).get_models())}
            if not binders or len(binders) != len(set(binders)) or len(targets) != len(set(targets)) or set(binders) & set(targets) or not set(binders + targets) <= available:
                raise ValueError("Sequence design requires unambiguous target and binder chain roles.")
            if design_chains and inputs.chain_list(design_chains) != binders:
                raise ValueError("This sequence stage must preserve the candidate's target and binder chain roles.")
            if row.get("stage") != "generation.backbone" or fix_interface_residues:
                if not native._candidate_sequence_matches_structure(row, path, binders):
                    raise ValueError("Candidate sequence does not match its structure; refusing unrestricted redesign.")
            if fix_interface_residues and not targets and (sequence_stage != "soluble_optimization" or optimization_enabled):
                raise ValueError("Interface fixing requires explicit target chains in the candidate complex.")
            if sequence_stage == "soluble_optimization" and optimization_enabled:
                raw = dict(row.get("raw_metadata") or {})
                raw["fixed_residues"] = inherited_fixed_residues(row, path, binders)
                row["raw_metadata"] = raw
    preserved = copy.deepcopy(rows) if sequence_stage == "soluble_optimization" else []
    if sequence_stage == "initial":
        preserved = [row for row in rows if row.get("stage") != "generation.backbone"]
        rows = [row for row in rows if row.get("stage") == "generation.backbone"]
    if sequence_stage == "soluble_optimization" and not optimization_enabled:
        rows = []
    if rows and not weights.is_file():
        raise ValueError(f"Missing {inputs.MODELS[model_type]['label']} checkpoint: {weights}")
    if not rows:
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import finish_job
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
        job = existing_job or create_job("design", job_type="sequence_design", tool="ligandmpnn", inputs={
            "source_run_dir": str(source_run_dir) if source_run_dir else None,
            "candidates_jsonl": str(candidates_jsonl) if candidates_jsonl else None,
            "structure_path": str(structure_path) if structure_path else None}, params=p)
        candidates = write_candidates(job.run_dir, model_type, preserved)
        finish_job(job.run_dir, success=True, result={"outputs": {"candidates": candidates}, "metrics": {"candidate_count": len(candidates)}})
        return job.run_dir
    cdr_fixed_by_candidate: dict[str, list[str]] = {}
    cdr_regions_by_candidate: dict[str, dict[str, list[str]]] = {}
    if auto_detect_cdrs:
        from mn_molecular_design_studio.sabr_runtime import (
            CDR_NAMES, cdr_fixed_residues, fixed_residues_from_annotation,
        )

        regions = list(cdr_regions or CDR_NAMES)
        if not set(regions) <= set(CDR_NAMES) or not regions:
            raise ValueError("Choose one or more of CDR1, CDR2, and CDR3.")
        autodetect_rows = rows
        if structure_path:
            info = inputs.inventory(Path(structure_path))
            binder_chains = inputs.chain_list(design_chains)
            autodetect_rows = [{
                "candidate_id": "input", "complex_pdb": str(Path(structure_path).resolve()),
                "binder_chains": binder_chains,
                "target_chains": [chain for chain in info["chains"] if chain not in binder_chains],
            }]
        for candidate in autodetect_rows:
            source_id = str(candidate.get("candidate_id") or "")
            source_path_text = candidate.get("complex_pdb") or candidate.get("binder_pdb")
            source_path = resolve_stored_path(
                source_path_text,
                run_dir=Path(source_run_dir) if source_run_dir else None,
                must_exist=True,
            )
            if source_path is None:
                raise ValueError(f"SAbR could not find the structure for candidate {source_id}.")
            chains = inputs.chain_list(design_chains) or inputs.chain_list(candidate.get("binder_chains"))
            if not chains:
                raise ValueError(f"SAbR could not identify an antibody chain for candidate {source_id}.")
            upstream_sabr = (candidate.get("raw_metadata") or {}).get("sabr_cdr_selection") or {}
            upstream_regions = upstream_sabr.get("residues_by_region")
            if upstream_regions and set(upstream_sabr.get("chain_ids") or ()) >= set(chains):
                fixed, selected = fixed_residues_from_annotation(source_path, chains, regions, upstream_regions)
            else:
                fixed, selected = cdr_fixed_residues(
                    source_path, chains, regions, chain_type=sabr_chain_type,
                    noise_level=float(sabr_noise_level), mode=sabr_mode,
                )
            cdr_fixed_by_candidate[source_id] = list(dict.fromkeys(fixed))
            cdr_regions_by_candidate[source_id] = selected
            raw = dict(candidate.get("raw_metadata") or {})
            raw["sabr_cdr_selection"] = {
                "chain_ids": chains, "regions": regions, "residues_by_region": selected,
                "chain_type": sabr_chain_type, "mode": sabr_mode,
                "noise_level": float(sabr_noise_level),
            }
            candidate["raw_metadata"] = raw
        if not cdr_fixed_by_candidate:
            raise ValueError("SAbR did not identify any CDRs to redesign.")
    job = existing_job or create_job("design", job_type="sequence_design", tool="ligandmpnn",
                                     inputs={k: str(v) for k, v in {"structure_path": structure_path, "source_run_dir": source_run_dir, "candidates_jsonl": candidates_jsonl}.items() if v}, params=p)
    source = job.run_dir / "artifacts/mpnn_source"
    source.mkdir(parents=True, exist_ok=True)
    # Save only the selected candidates; original jobs and structures remain immutable.
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
    from mn_molecular_design_studio.prediction_provenance import file_identity, image_identity
    from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
    runtime_manifest = load_manifest("foundry_mpnn" if chosen == "foundry" else "ligandmpnn")
    write_json(job.run_dir / "artifacts/mpnn_settings.json", p | {
        "checkpoint": file_identity(weights), "runtime": image_identity(runtime_manifest["image"])})
    if cdr_regions_by_candidate:
        write_json(job.run_dir / "artifacts/sabr_cdr_selection.json", cdr_regions_by_candidate)
    if chosen == "foundry":
        result = native.run_foundry_mpnn_sequence_design(source, candidate_file, batch_size=num_seq_per_target,
            model_type=model_type, gpu_device=gpu_device, existing_job=job, direct_input=p)
    else:
        direct_options = {k: v for k, v in p.items() if k not in {
            "auto_detect_cdrs", "cdr_regions", "sabr_chain_type", "sabr_mode", "sabr_noise_level",
        }}
        result = native.run_ligandmpnn_sequence_design(source, candidate_file, model_type=model_type,
            design_chains=design_chains, num_seq_per_target=num_seq_per_target, sampling_temp=sampling_temp,
            omit_aas=omit_aas, seed=p["seed"], accepted_stages=STAGES,
            fixed_residues_by_candidate=cdr_fixed_by_candidate,
            gpu_device=gpu_device, existing_job=job, direct_input=direct_options)
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    if read_json(result / "result.json").get("success") is True:
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
        redesigned = read_candidates(result)
        for child in redesigned:
            raw = child.setdefault("raw_metadata", {})
            raw["sequence_design"] = {"stage": sequence_stage, "model_type": model_type,
                "seed": p["seed"], "interface_contact": {"atoms": "heavy", "cutoff_angstrom": interface_cutoff},
                "parent_candidate_id": child["parents"][0]}
            # Component metadata must describe the child structure, while the
            # exact parent assembly and constraints remain in source_candidate.
            if raw.get("components"):
                from mn_molecular_design_studio.molecule_imports import parse_structure
                child_path = resolve_stored_path(child["complex_pdb"], run_dir=result, must_exist=True)
                components = parse_structure(child_path.read_bytes(), child_path.name, include_ligands=True)[0]["components"]
                original_components = raw["components"]
                from collections import Counter
                def context(rows):
                    return Counter((r["type"], r.get("sequence"), r.get("ccd"), r.get("smiles")) for r in rows if r["type"] != "protein")
                if context(components) != context(original_components):
                    raise ValueError("Sequence redesign changed the non-protein molecular context.")
                raw["components"] = components
                raw["studio_molecule_id"] = child["candidate_id"]
                raw["sequences_by_chain"] = native.refolding_workflow._sequences_by_chain(child_path)
        candidates = write_candidates(result, model_type, [*preserved, *redesigned])
        payload = read_json(result / "result.json")
        payload.setdefault("outputs", {})["candidates"] = candidates
        payload.setdefault("metrics", {})["candidate_count"] = len(candidates)
        write_json(result / "result.json", payload)
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


def run_abmpnn_sequence_design(
    source_run_dir: Path | None = None,
    candidates_jsonl: Path | None = None,
    structure_path: Path | None = None,
    design_chains: str = "",
    cdr_regions: list[str] | None = None,
    num_seq_per_target: int = 1,
    sampling_temp: float = 0.1,
    omit_aas: str = "C",
    seed: int | None = None,
    selected_candidate_ids: list[str] | None = None,
    sabr_chain_type: str = "auto",
    sabr_mode: str = "sabr",
    sabr_noise_level: float = 0.0,
    verify_sequences: bool = False,
    job_name: str = "abmpnn-cdr-redesign",
    gpu_device: object = "0",
    existing_job: JobPaths | None = None,
) -> Path:
    """Design selected antibody CDRs with the app's AbMPNN checkpoint."""
    if structure_path and (source_run_dir or candidates_jsonl):
        raise ValueError("Choose a structure file or existing design results, not both.")
    if bool(source_run_dir) != bool(candidates_jsonl):
        raise ValueError("Existing design results require both the source job and candidate collection.")
    return run_mpnn_sequence_design(
        source_run_dir=source_run_dir, candidates_jsonl=candidates_jsonl,
        structure_path=structure_path, model_type="ab_mpnn", design_chains=design_chains,
        num_seq_per_target=num_seq_per_target, sampling_temp=sampling_temp,
        omit_aas=omit_aas, seed=seed, selected_candidate_ids=selected_candidate_ids,
        verify_sequences=verify_sequences, job_name=job_name, gpu_device=gpu_device,
        auto_detect_cdrs=True, cdr_regions=cdr_regions,
        sabr_chain_type=sabr_chain_type, sabr_mode=sabr_mode,
        sabr_noise_level=sabr_noise_level, existing_job=existing_job,
    )
