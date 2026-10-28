"""Execute a saved visual graph through the normal scheduler job adapters."""
from __future__ import annotations

import copy
import hashlib
import inspect
import json
import re
import shutil
import signal
import tempfile
import time
import secrets
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.pipeline_graph import (
    NODE_CATALOG, OPTIMIZED_PREDICTION_NODE_ENGINES, fingerprint, scope_root, validate_graph,
)
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, update_status, utc_now, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.mimosa_motif import hydrate_candidate_fixed_motif


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _path(value, run_dir=None) -> Path:
    path = resolve_stored_path(value, run_dir=run_dir, must_exist=True)
    if path is None:
        raise ValueError(f"Pipeline input is unavailable: {value}")
    return path


def _absolute_candidates(source: dict) -> list[dict]:
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    run_dir = _path(source["run_dir"])
    rows = copy.deepcopy(read_candidates(_path(source["candidates_jsonl"], run_dir)))
    for row in rows:
        for key in ("target_pdb", "complex_pdb", "binder_pdb"):
            if row.get(key):
                resolved = resolve_stored_path(row[key], run_dir=run_dir)
                if resolved:
                    row[key] = str(resolved)
    return rows


def _selected_source_rows(kind: str, selected: list[dict], all_rows: list[dict]) -> list[dict]:
    """Keep source-loader annotations when materializing a pipeline input node.

    ``project_design_source`` enriches selected rows with verified MIMOSA
    fixed-residue identities. Re-reading the normalized candidate file is
    still useful for source kinds that only select rows, but must not replace
    those enriched project-design rows with their unhydrated originals.
    """
    if kind == "project_design_source":
        return copy.deepcopy(selected)
    selected_ids = {str(candidate.get("candidate_id") or "") for candidate in selected}
    return [candidate for candidate in all_rows
            if str(candidate.get("candidate_id") or "") in selected_ids]


def _local_source(parent: Path, node_id: str, candidates: list[dict], target: dict | None = None, analysis_runs: list[dict] | None = None, *, final_evidence=False) -> dict:
    """Write a new candidate collection, keeping upstream artifacts read-only."""
    destination = parent / "artifacts" / "pipeline_nodes" / node_id
    normalized = destination / "artifacts" / "normalized_candidates"
    normalized.mkdir(parents=True, exist_ok=True)
    candidate_path = normalized / "candidates.jsonl"
    candidate_path.write_text("".join(json.dumps(c) + "\n" for c in candidates))
    if final_evidence:
        import pandas as pd
        from mn_molecular_design_studio.output_vocabulary import evidence_rows
        pd.DataFrame(evidence_rows(candidates)).to_csv(normalized / "final_evidence.csv", index=False)
    write_json(destination / "metadata.json", {"status": "completed", "tool": "pipeline_collection", "run_id": node_id})
    write_json(destination / "input.json", {"inputs": {}, "params": {}})
    write_json(destination / "result.json", {"success": True, "outputs": {"candidates": candidates}})
    return {"run_dir": str(destination), "candidates_jsonl": str(candidate_path), "target": target or {}, "candidate_count": len(candidates), "analysis_runs": analysis_runs or []}


def source_designs(source_run_dir: str, project_id: str, selected_ids=None) -> tuple[Path, list[dict]]:
    from mn_molecular_design_studio.design_validation import design_candidates
    run_dir = _path(source_run_dir)
    marker = read_json(run_dir / "molecular_design_studio.json")
    if str(marker.get("project_id") or "") != project_id:
        raise ValueError("The selected design job must belong to this project.")
    candidates = design_candidates(run_dir)
    if selected_ids is not None:
        wanted = set(selected_ids)
        if not wanted.issubset({str(c.get("candidate_id")) for c in candidates}):
            raise ValueError("Only filter-passing designs from the selected successful job can be used.")
        candidates = [c for c in candidates if str(c.get("candidate_id")) in wanted]
    if not candidates:
        raise ValueError("This job has no successful, filter-passing designs.")
    return run_dir, candidates


def source_generations(source_run_dir: str, project_id: str, selected_ids=None, *, stage=None) -> tuple[Path, list[dict]]:
    from mn_molecular_design_studio.generation_io import generated_candidates
    run_dir=_path(source_run_dir)
    if str(read_json(run_dir/'molecular_design_studio.json').get('project_id') or '')!=project_id:
        raise ValueError('The selected generator job must belong to this project.')
    candidates=generated_candidates(run_dir,stage=stage)
    if selected_ids is not None:
        wanted=set(selected_ids)
        if not wanted.issubset({c['candidate_id'] for c in candidates}):
            raise ValueError('Select generated designs from the completed generator job.')
        candidates=[c for c in candidates if c['candidate_id'] in wanted]
    if not candidates:raise ValueError('This job has no completed generation outputs.')
    return run_dir,candidates


def source_project_designs(source_run_dir: str, project_id: str, selected_ids=None) -> tuple[Path, list[dict]]:
    """Load structure-bearing candidates from any completed design job in this project."""
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates

    run_dir = _path(source_run_dir)
    marker = read_json(run_dir / "molecular_design_studio.json")
    metadata = read_json(run_dir / "metadata.json")
    result = read_json(run_dir / "result.json")
    if str(marker.get("project_id") or "") != str(project_id):
        raise ValueError("The selected design job must belong to this project.")
    if metadata.get("task_group") != "design" or metadata.get("status") not in {"completed", "succeeded"} or result.get("success") is not True:
        raise ValueError("Choose a completed design job with successful outputs.")

    rows = []
    for candidate in read_candidates(run_dir / "artifacts/normalized_candidates/candidates.jsonl"):
        if candidate.get("stage") not in {"generation.backbone", "generation.backbone_sequence", "sequence_design", "complex_refolding"}:
            continue
        if not candidate.get("binder_chains"):
            continue
        structure = candidate.get("complex_pdb") or candidate.get("binder_pdb")
        if not structure:
            continue
        resolved = resolve_stored_path(structure, run_dir=run_dir, must_exist=True)
        if resolved is None or not resolved.is_file():
            continue
        row = copy.deepcopy(candidate)
        for key in ("target_pdb", "complex_pdb", "binder_pdb"):
            if row.get(key):
                path = resolve_stored_path(row[key], run_dir=run_dir, must_exist=True)
                if path:
                    row[key] = str(path)
        row = hydrate_candidate_fixed_motif(row, run_dir)
        rows.append(row)

    if selected_ids is not None:
        wanted = {str(value) for value in selected_ids}
        available = {str(row.get("candidate_id") or "") for row in rows}
        if not wanted or not wanted.issubset(available):
            raise ValueError("Select one or more structure-bearing candidates from this completed design job.")
        rows = [row for row in rows if str(row.get("candidate_id") or "") in wanted]
    if not rows:
        raise ValueError("This completed job has no structure-bearing design candidates.")
    return run_dir, rows


def source_molecule_selection(source_run_dir: str, project_id: str, selected_ids=None):
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    run_dir = _path(source_run_dir)
    marker = read_json(run_dir / "molecular_design_studio.json")
    metadata = read_json(run_dir / "metadata.json")
    snapshot = read_json(run_dir / "selection_snapshot.json")
    if (marker.get("project_id") != project_id or snapshot.get("project_id") != project_id
            or metadata.get("tool") != "database_selection" or metadata.get("status") != "completed"):
        raise ValueError("Use a completed Database selection snapshot from this project.")
    candidates = read_candidates(run_dir)
    expected = set(snapshot.get("molecule_ids") or [])
    if not expected or {row["candidate_id"] for row in candidates} != expected:
        raise ValueError("The Database selection snapshot is incomplete.")
    if selected_ids is not None:
        if not set(selected_ids).issubset(expected):
            raise ValueError("Selected molecules are not present in this frozen selection.")
        candidates = [row for row in candidates if row["candidate_id"] in selected_ids]
    return run_dir, candidates


def _generation_target(params, source=None):
    target=copy.deepcopy((source or {}).get('target') or {})
    own_input=bool(params.get('target_pdb'))
    target.update({k:v for k,v in params.items() if k in {'target_pdb','target_chains','hotspots','pipeline_task','ligand','motif_atom_spec'} and (own_input or v not in ('',None,[]))})
    if params.get('ligand_residue'):
        target['ligand_original_code']=params['ligand_residue']['name']
    chains=target.get('target_chains',[])
    if isinstance(chains,str):chains=[c for c in re.split(r'[,;\s]+',chains) if c]
    target['target_chains']=chains
    return target


def _optimized_prediction_tasks(engine: str, source: dict, params: dict) -> list[dict]:
    """Convert connected candidate complexes into shared prediction tasks."""
    import io
    from Bio.PDB import MMCIFIO, PDBIO, Select
    from mn_molecular_design_studio import design_validation
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES

    profile = ENGINE_PROFILES[engine]
    direct_input = source.get("molecule_task")
    if isinstance(direct_input, dict):
        task_input = copy.deepcopy(direct_input)
        task_type = str(task_input.get("task_type") or "")
        target_sequences = str(task_input.get("target_protein_sequences") or "")
        partner_sequences = str(task_input.get("partner_protein_sequences") or "")
        target_pdb = str(task_input.get("target_pdb") or "").strip()
        target_chain = str(task_input.get("target_chain_id") or "A").strip()
        if profile.get("input_kind") == "boltzgen":
            if not target_pdb or not partner_sequences:
                raise ValueError("BoltzGen Fold needs a target PDB path and a protein partner sequence.")
            if re.search(r"[:\n\r]", partner_sequences):
                raise ValueError("BoltzGen Fold accepts one protein partner sequence per task.")
            if len(target_chain) != 1 or target_chain.upper() == "Z":
                raise ValueError("BoltzGen Fold requires a one-character target chain ID other than Z.")
            target_path = _path(target_pdb)
            if target_path.suffix.lower() in {".cif", ".mmcif"}:
                from Bio.PDB import MMCIFParser
                structure = MMCIFParser(QUIET=True).get_structure("boltzgen-target", str(target_path))
            else:
                from Bio.PDB import PDBParser
                structure = PDBParser(QUIET=True).get_structure("boltzgen-target", str(target_path))
            model = next(structure.get_models())
            if target_chain not in model:
                raise ValueError(f"BoltzGen Fold cannot find target chain {target_chain!r} in {target_path.name}.")

            class TargetChain(Select):
                def accept_chain(self, candidate_chain):
                    return candidate_chain.id == target_chain

            stream = io.StringIO()
            writer = PDBIO()
            writer.set_structure(structure)
            writer.save(stream, TargetChain())
            return [{
                "name": "cofold_task", "target_pdb_text": stream.getvalue(),
                "target_chain_id": target_chain, "binder_sequence": partner_sequences,
            }]

        if task_type == "protein":
            task = {"name": "cofold_task", "protein_sequences": target_sequences}
            protein_count = len([part for part in re.split(r"[:\n\r]+", target_sequences) if part.strip()])
            target_chain_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"[:protein_count])
        elif task_type in {"protein_protein", "protein_macrocyclic"}:
            task = {
                "name": "cofold_task", "target_protein_sequences": target_sequences,
                "binder_sequence": partner_sequences,
            }
            target_count = len([part for part in re.split(r"[:\n\r]+", target_sequences) if part.strip()])
            partner_count = len([part for part in re.split(r"[:\n\r]+", partner_sequences) if part.strip()])
            chain_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
            target_chain_ids = chain_ids[:target_count]
            partner_chain_ids = chain_ids[target_count:target_count + partner_count]
        else:
            task = {"name": "cofold_task", "protein_sequences": target_sequences}
            protein_count = len([part for part in re.split(r"[:\n\r]+", target_sequences) if part.strip()])
            target_chain_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"[:protein_count])

        for molecule_key in ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles"):
            task[molecule_key] = str(task_input.get(molecule_key) or "")
        if profile.get("supports_msa"):
            task["molecule_chain_options"] = [
                {"id": chain_id, "use_msa": bool(params.get("target_use_msa", params.get("use_target_msa", True)))}
                for chain_id in target_chain_ids
            ]
            if task_type in {"protein_protein", "protein_macrocyclic"}:
                task["molecule_chain_options"].extend(
                    {"id": chain_id, "use_msa": False}
                    for chain_id in partner_chain_ids
                )
        if task_type == "protein_macrocyclic" and task_input.get("cyclic_partner"):
            partner_chain = partner_chain_ids[0]
            partner_option = next((option for option in task.setdefault("molecule_chain_options", [])
                                   if option.get("id") == partner_chain), None)
            if partner_option is None:
                partner_option = {"id": partner_chain, "use_msa": False}
                task["molecule_chain_options"].append(partner_option)
            partner_option["cyclic"] = True

        if target_pdb and profile.get("supports_templates") and params.get("use_target_template", True):
            target_path = _path(target_pdb)
            if target_path.suffix.lower() in {".cif", ".mmcif"}:
                from Bio.PDB import MMCIFParser
                structure = MMCIFParser(QUIET=True).get_structure("cofold-template", str(target_path))
            else:
                from Bio.PDB import PDBParser
                structure = PDBParser(QUIET=True).get_structure("cofold-template", str(target_path))
            model = next(structure.get_models())
            selected_chains = [part.strip() for part in re.split(r"[,;\s]+", target_chain) if part.strip()]
            templates = []
            for index, chain_id in enumerate(selected_chains):
                if chain_id not in model:
                    raise ValueError(f"Could not find target chain {chain_id!r} in {target_path.name}.")

                class TargetChain(Select):
                    def accept_chain(self, candidate_chain):
                        return candidate_chain.id == chain_id

                stream = io.StringIO()
                writer = MMCIFIO()
                writer.set_structure(structure)
                writer.save(stream, TargetChain())
                query_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
                if index >= len(target_chain_ids):
                    raise ValueError("The target structure contains more chains than the target protein sequences entered.")
                templates.append({
                    "name": f"target_{query_ids[index]}.cif", "text": stream.getvalue(),
                    "query_chain_id": query_ids[index], "template_chain_id": chain_id,
                })
            task["structure_templates"] = templates
        from mn_molecular_design_studio.batching import prediction_collection_tasks
        return prediction_collection_tasks(
            task, compound_set_file=task_input.get("compound_set_file"),
            sequence_file=task_input.get("sequence_file"),
        )

    rows = _absolute_candidates(source)
    if not rows:
        raise ValueError(f"{profile['label']} needs at least one connected candidate complex.")
    from mn_molecular_design_studio.batching import MAX_PREDICTION_ITEMS
    if len(rows) > MAX_PREDICTION_ITEMS:
        raise ValueError(f"{profile['label']} accepts up to {MAX_PREDICTION_ITEMS} candidate complexes per pipeline step; reduce the upstream selection.")
    if engine == "boltz2" and source.get("target", {}).get("pipeline_task") in {"ligand_binder", "motif_ligand"}:
        tasks = []
        for candidate in rows:
            code = (candidate.get("raw_metadata") or {}).get("ligand_original_code", "")
            if not re.fullmatch(r"[A-Z0-9]{1,5}", code):
                raise ValueError("Boltz-2 needs a CCD or SMILES ligand definition. This motif uses a custom ligand fragment; connect RosettaFold3 to preserve its input atom graph.")
            if not candidate.get("binder_sequence"):
                raise ValueError("The generated design has no binder sequence.")
            task = {"name": candidate["candidate_id"], "binder_sequence": candidate["binder_sequence"],
                    "ligand_ccd_codes": code, "compound_id": code, "compound_name": code,
                    "molecule_chain_options": [{"id": "A", "use_msa": False}],
                    "source_candidate_id": candidate["candidate_id"]}
            if source.get("run_dir"):
                task["source_run_dir"] = str(_path(source["run_dir"]))
            tasks.append(task)
        return tasks
    run_dir = _path(source["run_dir"])
    tasks = []
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    for candidate in rows:
        candidate_id = str(candidate.get("candidate_id") or "")
        if not candidate_id:
            raise ValueError(f"{profile['label']} cannot map a candidate without an ID back to its input.")
        library = candidate.get("raw_metadata") or {}
        if library.get("studio_molecule_id") and library.get("components"):
            from mn_molecular_design_studio.scientific_contracts import component_prediction_task
            task = component_prediction_task(
                candidate, run_dir, target_use_msa=bool(params.get("target_use_msa", True)),
            )
            tasks.append(task)
            continue
        info = design_validation.candidate_description(run_dir, candidate)
        if candidate.get("source_tool") in {"bindcraft2_generator", "boltzgen_generator", "esmfold2_binder_design"}:
            observed = ":".join(design_validation.sequence(info["chains"][chain]) for chain in info["binder_chains"])
            if observed.replace(":", "") != str(candidate.get("binder_sequence") or "").replace(":", ""):
                raise ValueError("Generator binder sequence differs from structure; cannot predict a mismatched VHH.")
        target_chains = list(info["target_chains"])
        binder_chains = list(info["binder_chains"])
        if not target_chains:
            raise ValueError(f"Candidate {candidate_id} has no identified target chain for {profile['label']} cofolding.")
        ordered = [*target_chains, *binder_chains]
        if len(ordered) > len(chain_ids):
            raise ValueError(f"Candidate {candidate_id} exceeds {profile['label']}'s protein-chain limit.")
        chain_map = {chain: chain_ids[index] for index, chain in enumerate(ordered)}
        if profile["api_engine"] == "boltzgen_fold_opt_msa":
            if len(target_chains) != 1:
                raise ValueError("BoltzGen Fold needs exactly one target chain; select a single target chain before connecting this node.")
            target_chain = target_chains[0]
            if len(target_chain) != 1 or target_chain.upper() == "Z":
                raise ValueError("BoltzGen Fold requires a one-character target chain ID other than Z.")
            structure = design_validation.read_structure(info["structure_path"])
            class TargetChain(Select):
                def accept_chain(self, candidate_chain):
                    return candidate_chain.id == target_chain
            stream = io.StringIO()
            writer = PDBIO()
            writer.set_structure(structure)
            writer.save(stream, TargetChain())
            tasks.append({
                "name": candidate_id,
                "target_pdb_text": stream.getvalue(),
                "target_chain_id": target_chain,
                "binder_sequence": info["binder_sequence"],
                "source_candidate_id": candidate_id,
                "source_run_dir": str(run_dir),
            })
            continue

        task = {
            "name": candidate_id,
            "target_protein_sequences": ":".join(design_validation.sequence(info["chains"][chain]) for chain in target_chains),
            "binder_sequence": info["binder_sequence"],
            "molecule_chain_options": [
                {"id": chain_map[chain], "use_msa": bool(profile.get("supports_msa") and params.get("target_use_msa", True) and chain in target_chains)}
                for chain in ordered
            ],
            "source_candidate_id": candidate_id,
            "source_run_dir": str(run_dir),
        }
        if params.get("use_target_template", True) and profile.get("supports_templates"):
            structure = design_validation.read_structure(info["structure_path"])
            templates = []
            for chain in target_chains:
                class TargetChain(Select):
                    def accept_chain(self, candidate_chain):
                        return candidate_chain.id == chain
                stream = io.StringIO()
                writer = MMCIFIO()
                writer.set_structure(structure)
                writer.save(stream, TargetChain())
                templates.append({
                    "name": f"target_{chain_map[chain]}.cif",
                    "text": stream.getvalue(),
                    "query_chain_id": chain_map[chain],
                    "template_chain_id": chain,
                })
            task["structure_templates"] = templates
        tasks.append(task)
    return tasks


def _af3_pipeline_tasks(source: dict, params: dict) -> list[dict]:
    """Convert connected candidate complexes into official AF3 tasks."""
    return _optimized_prediction_tasks("alphafold3", source, params)


def _optimized_prediction_request(engine: str, params: dict, source: dict, name: str) -> dict:
    """Build a shared prediction request for a direct molecule task or design candidates."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES

    if engine == "boltz2":
        from mn_molecular_design_studio.prediction_configuration import boltz_settings
        params = boltz_settings(params)
    profile = ENGINE_PROFILES[engine]
    option_keys = {setting["key"] for setting in profile.get("settings", ())}
    engine_options = {
        setting["key"]: params[setting["key"]]
        for setting in profile.get("settings", ())
        if setting["key"] in params and setting["key"] not in {"use_msa", "use_msa_server"}
    }
    if engine == "boltz2":
        engine_options["write_full_pae"] = True
    msa_setting = "target_use_msa" if "target_use_msa" in params else "use_target_msa"
    if "use_msa" in option_keys:
        engine_options["use_msa"] = bool(params.get(msa_setting, True))
    ligand_only = engine == "boltz2" and source.get("target", {}).get("pipeline_task") in {"ligand_binder", "motif_ligand"}
    if ligand_only:
        engine_options["use_msa"] = False
    if engine == "alphafold3":
        engine_options["model"] = params.get("model", "alphafold3")
    terms_accepted = bool(params.get("accept_alphafold3_terms", False))
    return {
        "schema_version": 1,
        "workflow": "refolding.optimized_prediction",
        "parameters": {
            "engine": engine,
            "optimization_mode": params.get("optimization_mode", "off"),
            "job_name": name,
            "seed": int(params.get("seed", 0)),
            "samples": int(params.get("samples", 1)),
            "num_recycles": int(params.get("num_recycles", profile.get("default_recycles", 10))),
            "cpu_cores": 4,
            "tasks": _optimized_prediction_tasks(engine, source, params),
            "engine_options": engine_options,
            "use_target_template": bool(params.get("use_target_template", profile.get("requires_target_template", False))),
            "target_use_msa": False if ligand_only else bool(params.get(msa_setting, True)),
            "accept_alphafold3_terms": terms_accepted,
        },
    }


def filter_candidates(rows: list[dict], settings: dict) -> list[dict]:
    from mn_molecular_design_studio.candidate_selection import evaluate_selection, selected_candidates
    outcomes = evaluate_selection(rows, settings)
    if outcomes and all(row["raw_metadata"]["selection"]["status"] == "not_evaluated" for row in outcomes):
        return outcomes[:settings["limit"]] if settings.get("limit") else outcomes
    return selected_candidates(outcomes)


def node_request(kind: str, params: dict, source: dict, name: str) -> dict:
    """Only allowlisted adapters can turn graph data into a workflow call."""
    workflow = NODE_CATALOG[kind]["workflow"]
    p = copy.deepcopy(params)
    if kind == "af2" and source.get("target", {}).get("pipeline_task") in {"ligand_binder", "motif_ligand"}:
        raise ValueError("AlphaFold2 does not predict ligand complexes. Use RosettaFold3, or ESMFold for the binder alone.")
    if kind == "peptide_trap":
        p["job_name"] = name
        return {"schema_version": 1, "workflow": workflow, "parameters": p}
    if kind == "kejia_idr":
        p["template_structure"] = str(_path(p["template_structure"])) if p.get("template_structure") else None
        p["job_name"] = name
        return {"schema_version": 1, "workflow": workflow, "parameters": p}
    if kind == "boltz2":
        # Historical saved graphs retain target MSA by default. New recipes may
        # explicitly select single-sequence conditioning; binder chains retain
        # their per-chain policy through _optimized_prediction_tasks.
        p.setdefault("target_use_msa", True)
    if kind in {"rfd4", "rfd4_fold"}:
        from mn_molecular_design_studio import design_validation
        target=source.get("target") or {}
        p["target_structure"]=p.get("target_structure") or target.get("target_pdb") or None
        p["job_name"]=name
        if kind=="rfd4_fold":
            if source.get("molecule_task"):
                task=source["molecule_task"]
                if any(task.get(k) for k in ("dna_sequences","rna_sequences","ligand_smiles","ligand_ccd_codes")):
                    raise ValueError("The RFD4 protein folding node accepts protein sequences; use the design tool for RNA or molecule-conditioned tasks.")
                p["protein_sequences"]=":".join(str(task.get(k) or "").strip() for k in ("target_protein_sequences","partner_protein_sequences") if task.get(k))
            if source.get("candidates_jsonl"):
                rows=_absolute_candidates(source)
                p["input_files"]=[str(design_validation.candidate_description(_path(source["run_dir"]),row)["structure_path"]) for row in rows]
                p.setdefault('options',{})['source_candidates']={Path(path).name:{key:row.get(key) for key in ('candidate_id','target_chains','binder_chains','complex_pdb','target_pdb')} for path,row in zip(p['input_files'],rows)}
                p["target_structure"]=None
        else:
            options=json.loads(p.pop("options_json","{}"))
            if not isinstance(options,dict): raise ValueError("RFD4 task options must be a JSON object.")
            p["options"]={**options,"refold":False}
            p["input_files"]=[x.strip() for x in str(p.get("input_files") or "").splitlines() if x.strip()]
            chains=target.get("target_chains") or []
            p["target_chains"]=[x for x in re.split(r"[,;\s]+",chains) if x] if isinstance(chains,str) else chains
            p["hotspots"]=target.get("hotspots","")
            campaign=p.pop("distributed_campaign",None)
            if campaign is not None:
                if p.get("task")!="protein_binder" or p.get("input_mode") not in {"guided","guided_binder"} or p["input_files"]:
                    raise ValueError("RFD4 distributed campaigns support one guided protein-binder target only.")
                campaign={**campaign,"generator_id":"rfd4_proteina"}
                campaign.setdefault("seed",int(p.get("seed") or 0))
                from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import validate_campaign_config
                p["distributed_campaign"]=validate_campaign_config(campaign)
        return {"schema_version":1,"workflow":workflow,"parameters":p}
    if isinstance(source.get("molecule_task"), dict):
        molecule_task = source["molecule_task"]
        if kind == "af2":
            if p.get("use_binder_template") or p.get("use_interface_template"):
                raise ValueError("Target-only evaluation cannot use binder or interface templates.")
            if p.get("af2_validation_mode") == "conditioned_self_consistency":
                raise ValueError("Conditioned self-consistency requires a designed complex. Choose target_only for a target plus binder sequence.")
            target_pdb = str(molecule_task.get("target_pdb") or "").strip()
            partner_sequence = str(molecule_task.get("partner_protein_sequences") or "").strip()
            if not target_pdb or not partner_sequence:
                raise ValueError("AlphaFold2 Initial Guess needs a target PDB path and partner protein sequence.")
            target_chains = [value for value in re.split(r"[,;\s]+", str(molecule_task.get("target_chain_id") or "A")) if value]
            return {
                "schema_version": 1,
                "workflow": "refolding.af2_initial_guess_input",
                "parameters": {
                    "input_mode": "target_template", "job_name": name,
                    "target_pdb": str(_path(target_pdb)), "target_chains": target_chains,
                    "partner_sequence": partner_sequence,
                    "num_recycles": int(p.get("num_recycles", 3)),
                    "model_count": int(p.get("model_count", 1)),
                    "multimer": bool(p.get("multimer", True)),
                    "msa_source": p.get("msa_source", "single_sequence"), "msa_max_depth": int(p.get("msa_max_depth", 256)),
                    "seed": int(p.get("seed", 0)),
                },
            }
        direct_engines = {"af3": "alphafold3", "boltz2": "boltz2", "rf3": "rf3_foundry"}
        engine = direct_engines.get(kind) or OPTIMIZED_PREDICTION_NODE_ENGINES.get(kind)
        if engine:
            return _optimized_prediction_request(engine, p, source, name)
    if kind == "boltz2":
        # Use the same MSA-aware Boltz-2 runner as the prediction tile. The
        # legacy initial-guess adapter only injected an MSA when a benchmark
        # CSV happened to be supplied, despite the graph tile promising one.
        p.setdefault("target_use_msa", True)
        return _optimized_prediction_request("boltz2", p, source, name)
    if kind == "surface_analysis":
        target = source.get("target") or {}
        target_pdb = target.get("target_pdb")
        chains = target.get("target_chains") or []
        chains = [value for value in re.split(r"[,;\s]+", chains) if value] if isinstance(chains, str) else list(chains)
        if not target_pdb or not chains:
            raise ValueError("Protein Surface Analysis needs a connected target with selected protein chains.")
        masif_chain_id = str(p.get("masif_chain_id") or chains[0])
        if masif_chain_id not in chains:
            raise ValueError("Choose a MaSIF-seed chain from the connected target chains.")
        return {
            "schema_version": 1,
            "workflow": workflow,
            "parameters": {
                "target_pdb": target_pdb,
                "chain_ids": chains,
                "pesto_interface": p.get("pesto_interface", "protein"),
                "scannet_mode": p.get("scannet_mode", "interface"),
                "scannet_use_msa": bool(p.get("scannet_use_msa", False)),
                "masif_chain_id": masif_chain_id,
            },
        }
    if kind == "rfdiffusion_beta_pairing":
        inherited = source.get("target") or {}
        p["target_pdb"] = p.get("target_pdb") or inherited.get("target_pdb") or None
        inherited_chains = inherited.get("target_chains") or []
        if isinstance(inherited_chains, str):
            inherited_chains = [value for value in re.split(r"[,;\s]+", inherited_chains) if value]
        p["target_chain"] = p.get("target_chain") or (inherited_chains[0] if len(inherited_chains) == 1 else "")
        p["job_name"] = name
        return {"schema_version": 1, "workflow": workflow, "parameters": p}
    from mn_molecular_design_studio import inverse_folding_input
    if kind in inverse_folding_input.ENGINES:
        if source.get("candidates_jsonl"):
            p.update(source_run_dir=source["run_dir"], candidates_jsonl=source["candidates_jsonl"])
        else:
            target = source.get("target") or {}
            p["structure_path"] = target.get("target_pdb")
            if not p.get("design_chains"):
                chains = target.get("target_chains") or []
                p["design_chains"] = chains if isinstance(chains, str) else ",".join(chains)
        p["job_name"] = name
        return {"schema_version": 1, "workflow": workflow, "parameters": p}
    from mn_molecular_design_studio import generator_input
    if kind in generator_input.ENGINES:
        inherited=source.get('target') or {}
        own=bool(p.get('target_pdb'))
        p['target_pdb']=p.get('target_pdb') or inherited.get('target_pdb') or None
        chains=p.get('target_chains') if own else inherited.get('target_chains',p.get('target_chains'))
        p['target_chains']=[v for v in re.split(r'[,;\s]+',chains) if v] if isinstance(chains,str) else list(chains or [])
        if not own and not p.get('hotspots'):p['hotspots']=inherited.get('hotspots','')
        if NODE_CATALOG[kind].get('hotspot_support') == 'unsupported':
            p['hotspots'] = ''
        if p.get('task') in {'unconditional','symmetry'}:p.update(target_pdb=None,target_chains=[],hotspots='')
        for key in ('scaffold_pdb','scaffold_ss_file','scaffold_adj_file','target_ss_file','target_adj_file'):
            if key in p:p[key]=p[key] or None
        p['campaign_name']=name
        return {'schema_version':1,'workflow':workflow,'parameters':p}
    if kind=='complexa_generator':
        target=_generation_target(p,source)
        p.update(target_pdb=target.get('target_pdb',''),target_chains=target['target_chains'],hotspots=target.get('hotspots',''),campaign_name=name)
        if p['pipeline_task']!='protein_binder':p.update(target_chains=[],hotspots='')
        return {"schema_version":1,"workflow":workflow,"parameters":p}
    if NODE_CATALOG[kind]["category"] in {"Design engine", "Design pipeline"}:
        target = source.get("target") or {}
        if not target.get("target_pdb"):
            raise ValueError("This design step needs a target structure.")
        p.update(target_pdb=target["target_pdb"], target_chains=target["target_chains"], hotspots=target.get("hotspots", ""), campaign_name=name)
        if NODE_CATALOG[kind].get("hotspot_support") == "unsupported":
            p["hotspots"] = ""
        if NODE_CATALOG[kind].get("coldspot_support") == "unsupported":
            p.pop("coldspots", None)
            p.pop("non_binding_residues", None)
        if kind in {"bindcraft2", "bindcraft2_generator"}:
            p["coldspots"] = p.get("coldspots") or target.get("coldspots", "")
        elif kind in {"boltzgen", "boltzgen_generator"}:
            p["non_binding_residues"] = p.get("non_binding_residues") or target.get("coldspots", "")
        if kind == "rfd3":
            length = str(p.pop("binder_length"))
            match = re.fullmatch(r"(\d+)(?:-(\d+))?", length)
            if not match or int(match[1]) < 1 or int(match[2] or match[1]) < int(match[1]):
                raise ValueError("Binder length must be a positive length or range, such as 74 or 60-100.")
            p.update(task="protein_binder", target_structure=p.pop("target_pdb"), length_min=int(match[1]), length_max=int(match[2] or match[1]), job_name=p.pop("campaign_name"))
        elif kind == "rfdiffusion":
            from mn_molecular_design_studio.engine.workflow_engine.workflows.design import default_target_contig
            p["contig"] = default_target_contig(Path(p["target_pdb"]), p["target_chains"], p["binder_length"])
            if not p["contig"]:
                raise ValueError("RFdiffusion requires a PDB target with readable selected chains.")
            p["run_vanilla_pipeline"] = False
        elif kind in {"genie3", "genie3_pipeline"}:
            p["run_mode"] = "generation_only" if kind == "genie3" else "full_vanilla_pipeline"
        elif kind in {"pxdesign", "pxdesign_pipeline"}:
            p["run_mode"] = "generation_only" if kind == "pxdesign" else "pipeline"
        elif kind == "bindcraft2":
            attempts = int(p.pop("num_designs")); p.update(number_of_final_designs=attempts, max_trajectories=attempts)
            campaign=p.get("distributed_campaign")
            if campaign is not None:
                campaign={**campaign,"generator_id":"bindcraft2"}
                campaign.setdefault("seed",int(p.get("campaign_seed") or 0))
                from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import validate_campaign_config
                p["distributed_campaign"]=validate_campaign_config(campaign)
        elif kind == "bindcraft2_generator":
            trajectories = int(p.pop("num_designs")); p.update(number_of_final_designs=1, max_trajectories=trajectories)
        elif kind == "boltzgen_generator":
            # The generator workflow pins the generation-only contract in the
            # workflow API adapter; this node exposes only generation controls.
            p.pop("budget", None)
        elif kind == "protpardelle_pipeline":
            p["num_mpnn_seqs"] = int(p.pop("sequences"))
        elif kind == "esm_binder":
            chains = p.pop("target_chains")
            if len(chains) != 1:
                raise ValueError("ESMFold2 Binder Design currently accepts one target chain.")
            p["target_chain"] = chains[0]
    elif kind == "af3":
        p = {
            "engine": "alphafold3",
            "optimization_mode": p["optimization_mode"],
            "job_name": name,
            "seed": 0,
            "samples": int(p["samples"]),
            "num_recycles": int(p["num_recycles"]),
            "cpu_cores": 4,
            "accept_alphafold3_terms": bool(p["accept_alphafold3_terms"]),
            "tasks": _af3_pipeline_tasks(source, p),
            "engine_options": {"model": p["model"], "use_msa": bool(p["target_use_msa"])},
            "use_target_template": bool(p["use_target_template"]),
            "target_use_msa": bool(p["target_use_msa"]),
        }
    elif kind in OPTIMIZED_PREDICTION_NODE_ENGINES:
        from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES

        engine = OPTIMIZED_PREDICTION_NODE_ENGINES[kind]
        profile = ENGINE_PROFILES[engine]
        option_keys = {setting["key"] for setting in profile.get("settings", ())}
        engine_options = {
            setting["key"]: p[setting["key"]]
            for setting in profile.get("settings", ())
            if setting["key"] in p and setting["key"] not in {"use_msa", "use_msa_server"}
        }
        if "use_msa" in option_keys:
            engine_options["use_msa"] = bool(p.get("target_use_msa", True))
        p = {
            "engine": engine,
            "optimization_mode": p.get("optimization_mode", "off"),
            "job_name": name,
            "seed": int(p.get("seed", 0)),
            "samples": int(p.get("samples", 1)),
            "num_recycles": int(p.get("num_recycles", profile.get("default_recycles", 10))),
            "cpu_cores": 4,
            "tasks": _optimized_prediction_tasks(engine, source, p),
            "engine_options": engine_options,
            "use_target_template": bool(p.get("use_target_template", profile.get("requires_target_template", False))),
            "target_use_msa": bool(p.get("target_use_msa", True)),
        }
    else:
        p.update(source_run_dir=source["run_dir"], candidates_jsonl=source["candidates_jsonl"])
        if kind == "abmpnn":
            from mn_molecular_design_studio import mpnn_input
            p["num_seq_per_target"] = int(p.pop("sequences", 1))
            p["cdr_regions"] = [
                region for region, enabled in (
                    ("CDR1", p.pop("redesign_cdr1", True)),
                    ("CDR2", p.pop("redesign_cdr2", True)),
                    ("CDR3", p.pop("redesign_cdr3", True)),
                ) if enabled
            ]
            if not p["cdr_regions"]:
                raise ValueError("Select at least one CDR region for AbMPNN redesign.")
            p["seed"] = int(p.get("seed") or 0) or None
            p["job_name"] = name
            # Validate common sequence-design controls; the dedicated workflow
            # pins the model to AbMPNN and derives fixed framework residues from SAbR.
            p = mpnn_input.prepare(p | {"model_type": "ab_mpnn", "backend": "shared_ligandmpnn"})
            from mn_molecular_design_studio.engine.workflow_engine.workflows.mpnn import run_abmpnn_sequence_design
            allowed = inspect.signature(run_abmpnn_sequence_design).parameters
            p = {key: value for key, value in p.items() if key in allowed}
        elif kind in {"foundry_mpnn", "ligandmpnn"}:
            from mn_molecular_design_studio import mpnn_input
            p["num_seq_per_target"] = int(p.pop("sequences")) * int(p.pop("number_of_batches", 1))
            p["noise_level"] = float(p["noise_level"]) if p.get("noise_level") not in (None, "") else None
            p = mpnn_input.prepare(p)
            from mn_molecular_design_studio.engine.workflow_engine.workflows.mpnn import run_mpnn_sequence_design
            allowed = inspect.signature(run_mpnn_sequence_design).parameters
            p = {k: v for k, v in p.items() if k in allowed}
        elif kind in {"rf3", "boltz2", "af2", "protenix"}:
            p["require_monomer_success"] = False
            if kind == "af2":
                from mn_molecular_design_studio.prediction_configuration import af2_mode
                p["af2_validation_mode"], p["use_initial_guess"] = af2_mode(
                    p.get("af2_validation_mode"), p.get("use_initial_guess"), legacy_initial_guess=True)
            if kind == "boltz2":
                p.pop("use_target_msa", None)  # This adapter always uses shared target MSA preparation.
        elif kind in {"openfold3_legacy", "protenix_cli", "boltzgen_fold_batch"}:
            # These standalone cofolding tiles consume normalized candidate
            # collections directly. Pipeline candidates may come from design,
            # generation, sequence redesign, or a preceding predictor, so do
            # not require the old monomer-refolding stage gate here.
            p["require_monomer_success"] = False
    return {"schema_version": 1, "workflow": workflow, "parameters": p}


def _submit(request: dict, parent: Path, node: dict, *, owner_id: str, project_id: str, gpu_device: str) -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import submit_workflow_request
    from mn_molecular_design_studio.workspace import mark_studio_job
    request = copy.deepcopy(request)
    from mn_molecular_design_studio.quality_tools import WORKFLOWS as quality_workflows
    workflow = str(request.get("workflow", ""))
    cpu_analysis = workflow.startswith("analysis.") and not quality_workflows.get(workflow, {}).get("gpu", False)
    request["resources"] = {**request.get("resources", {}), "gpu_device": "none" if cpu_analysis or workflow in {"design.kejia_idr", "design.peptide_trap"} else gpu_device, "cpu_cores": request.get("resources", {}).get("cpu_cores", 4)}
    if cpu_analysis:
        request["resources"]["min_vram_gb"] = 0
    distributed_campaign = request.get("parameters", {}).get("distributed_campaign")
    if (
        workflow in {"design.bindcraft2", "design.rfd4_proteina"}
        and isinstance(distributed_campaign, dict)
        and distributed_campaign.get("enabled")
    ):
        distributed_campaign.setdefault("gpu_selection", str(gpu_device))
        request["parameters"]["distributed_campaign"] = distributed_campaign
        request["resources"]["gpu_device"] = "none"
        request["resources"]["min_vram_gb"] = 0
    if workflow in {"design.kejia_idr", "design.peptide_trap"}:
        request["resources"]["cpu_cores"] = int(request["parameters"].get("cpu_cores", 4))
    if node.get('candidate_batch') and request['resources']['gpu_device'] not in {'none', 'cpu'}:
        from mn_molecular_design_studio.cluster.gpu_runtime import profiles, request_images, compatible
        from mn_molecular_design_studio.cluster.capabilities import gpu_details
        runtime_profiles = profiles(request_images(request))
        devices = gpu_details()
        eligible = [i for i, device in devices.items() if compatible(device, runtime_profiles)]
        selector = str(request['resources']['gpu_device'])
        if runtime_profiles and not eligible:
            raise ValueError('No local GPU satisfies the installed CUDA binary targets.')
        if runtime_profiles and selector not in {'all', 'auto', 'any'} and any(i not in eligible for i in selector.removeprefix('device=').split(',')):
            raise ValueError('The selected GPU does not satisfy the installed CUDA binary targets.')
        if selector in {'auto', 'any'} and runtime_profiles and set(eligible) != set(devices):
            request['resources']['gpu_device'] = eligible[int(node.get('batch_offset', 0)) % len(eligible)]
        elif selector == 'all' and runtime_profiles and set(eligible) != set(devices):
            raise ValueError('The requested GPU group includes an incompatible CUDA architecture.')
    path = parent / "requests" / f"{node['id']}.json"
    # Workflow submission resolves relative paths against the request file's
    # parent directory. `write_json` stores managed paths relative to the
    # enclosing run directory instead, which is a different base here because
    # requests live in a subdirectory. Keep these operational requests as JSON
    # with their absolute/portable input paths intact.
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(request, indent=2, sort_keys=True) + "\n")
    result = submit_workflow_request(path)
    mark_studio_job(result, owner_id=owner_id, project_id=project_id)
    child = _path(result["job"]["run_dir"])
    metadata = read_json(child / "metadata.json")
    from mn_molecular_design_studio.cluster import config as cluster_config
    metadata["pipeline_parent_node_id"] = cluster_config.identity()
    metadata.update(pipeline_parent_run_dir=str(parent), campaign_id=parent.name, campaign_name=read_json(parent / "input.json").get("params", {}).get("pipeline_name", "Pipeline"), graph_node_id=node["id"])
    write_json(child / "metadata.json", metadata)
    return child


def _child_output(child: Path, target: dict) -> dict:
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    candidate_path = child / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    if not candidate_path.is_file():
        raise ValueError(f"Job {child.name} completed without a candidate collection.")
    candidates = read_candidates(candidate_path)
    params=read_json(child/'input.json').get('params',{})
    if params.get('generator_only'):
        target=_generation_target(params)
    # Bundled tools retain intermediate stage records; the next node receives
    # only final outputs, rather than redesigning/refolding all three pools.
    final = [c for c in candidates if c.get("metrics", {}).get("candidate_pool_level") == "native_pipeline"]
    if final and len(final) != len(candidates):
        output = _local_source(child, "final_outputs", _absolute_candidates({"run_dir": str(child), "candidates_jsonl": str(candidate_path)}), target)
        rows = [c for c in _absolute_candidates(output) if c.get("metrics", {}).get("candidate_pool_level") == "native_pipeline"]
        return _local_source(child, "final_outputs", rows, target)
    return {"run_dir": str(child), "candidates_jsonl": str(candidate_path), "target": target, "candidate_count": len(candidates)}


def _attach_refold_input_references(output: dict, source: dict, parent: Path, node_id: str) -> dict:
    """Keep the exact input backbone alongside each refolded result for pose RMSD."""
    source_rows = _absolute_candidates(source) if source.get("candidates_jsonl") else []
    by_id = {str(row.get("candidate_id") or ""): row for row in source_rows}
    rows = _absolute_candidates(output)
    for row in rows:
        raw = dict(row.get("raw_metadata") or {})
        direct = raw.get("source_candidate") if isinstance(raw.get("source_candidate"), dict) else {}
        possible_ids = [str(value) for value in row.get("parents") or [] if str(value)]
        if direct.get("candidate_id"):
            possible_ids.insert(0, str(direct["candidate_id"]))
        reference = next((by_id[value] for value in possible_ids if value in by_id), None)
        if reference is None and direct:
            reference = direct
        if not reference:
            continue
        raw.setdefault("source_candidate", copy.deepcopy(reference))
        structure = reference.get("complex_pdb") or reference.get("binder_pdb")
        if not structure:
            continue
        raw["pipeline_refold_input"] = {
            "complex_pdb": str(structure),
            "binder_chains": list(reference.get("binder_chains") or []),
            "target_chains": list(reference.get("target_chains") or []),
            "source_candidate_id": str(reference.get("candidate_id") or ""),
        }
        row["raw_metadata"] = raw
    return _local_source(parent, f"{node_id}_pose_inputs", rows, output.get("target"))


def _add_target_aligned_pose_rmsd(candidate: dict, parent: Path) -> dict:
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    reference = raw.get("pipeline_refold_input") if isinstance(raw.get("pipeline_refold_input"), dict) else None
    if not reference:
        return candidate
    metrics = dict(candidate.get("metrics") or {})
    reference_path = Path(str(reference.get("complex_pdb") or ""))
    prediction_path = Path(str(candidate.get("complex_pdb") or ""))
    if not reference_path.is_file() or not prediction_path.is_file():
        metrics["target_aligned_binder_rmsd"] = None
        metrics["target_aligned_binder_rmsd_status"] = "missing reference or Boltz-2 complex"
        candidate["metrics"] = metrics
        return candidate
    try:
        from mn_molecular_design_studio.engine.workflow_engine.workflows import mpnn as mpnn_workflow
        from mn_molecular_design_studio.engine.workflow_engine.workflows.benchmark import _target_aligned_binder_pose_metrics

        with tempfile.TemporaryDirectory(prefix="pipeline-pose-rmsd-") as temporary:
            temp = Path(temporary)
            reference_pdb, prediction_pdb = temp / "reference.pdb", temp / "prediction.pdb"
            reference_map = mpnn_workflow._pdb_copy(reference_path, reference_pdb)
            prediction_map = mpnn_workflow._pdb_copy(prediction_path, prediction_pdb)
            ref_binder = [reference_map[chain] for chain in reference.get("binder_chains") or [] if chain in reference_map]
            ref_target = [reference_map[chain] for chain in reference.get("target_chains") or [] if chain in reference_map]
            pred_binder = [prediction_map[chain] for chain in candidate.get("binder_chains") or [] if chain in prediction_map]
            pred_target = [prediction_map[chain] for chain in candidate.get("target_chains") or [] if chain in prediction_map]
            result = _target_aligned_binder_pose_metrics(
                reference_pdb, prediction_pdb,
                reference_binder_chains=ref_binder,
                reference_target_chains=ref_target,
                prediction_binder_chains=pred_binder or ref_binder,
                prediction_target_chains=pred_target or ref_target,
            )
        metrics.update(result)
        metrics["target_aligned_binder_rmsd_status"] = result.get("pose_rmsd_status", "")
    except Exception as exc:
        metrics["target_aligned_binder_rmsd"] = None
        metrics["target_aligned_binder_rmsd_status"] = f"calculation failed: {type(exc).__name__}: {exc}"
    candidate["metrics"] = metrics
    return candidate


def _cache_valid(output: dict) -> bool:
    try:
        from mn_molecular_design_studio.cache_provenance import file_digest
        manifest = output.get("content_manifest")
        if not manifest or any(not Path(path).is_file() or file_digest(path) != digest for path, digest in manifest.items()):
            return False
        analysis_runs = list(output.get("analysis_runs") or [])
        if output.get("analysis_run"):
            analysis_runs.append(output["analysis_run"])
        for analysis_run in analysis_runs:
            if not analysis_run.get("success"):
                return False
            root = _path(analysis_run.get("run_dir"))
            result = read_json(root / "result.json")
            if read_json(root / "metadata.json").get("status") != "completed" or result.get("success") is not True:
                return False
            artifacts = result.get("outputs", {}).get("artifacts") or []
            if not artifacts or any(not _path(item["path"], root).is_file() for item in artifacts):
                return False
        if output.get("target", {}).get("target_pdb") and not _path(output["target"]["target_pdb"]).is_file():
            return False
        if output.get("analysis_run"):
            return True
        path = _path(output["candidates_jsonl"])
        if output.get("review_candidates_jsonl"):
            review_path = _path(output["review_candidates_jsonl"])
            if output.get("review_digest") != file_digest(review_path):
                return False
        root = _path(output["run_dir"])
        if read_json(root / "metadata.json").get("status") != "completed" or read_json(root / "result.json").get("success") is not True:
            return False
        for row in _absolute_candidates(output):
            structure = row.get("complex_pdb") or row.get("binder_pdb")
            if structure and not Path(structure).is_file():
                return False
        if output.get("target", {}).get("target_pdb") and not Path(output["target"]["target_pdb"]).is_file():
            return False
        return path.is_file()
    except (OSError, ValueError, KeyError, TypeError):
        return False


def _child_analysis_output(child: Path, target: dict) -> dict:
    metadata = read_json(child / "metadata.json")
    result = read_json(child / "result.json")
    return {
        "analysis_run": {
            "run_dir": str(child),
            "run_id": str(metadata.get("run_id") or child.name),
            "task_group": str(metadata.get("task_group") or "detection"),
            "job_code": str(metadata.get("job_code") or ""),
            "success": result.get("success") is True,
        },
        "target": copy.deepcopy(target),
        "candidate_count": 0,
    }


def cancel_graph_children(run_dir: Path) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
    from mn_molecular_design_studio.candidate_execution import cancel
    cancel(run_dir, read_json(run_dir / "graph_execution.json"))
    for node in read_json(run_dir / "graph_execution.json").get("nodes", {}).values():
        raw = node.get("job_run_dir")
        if not raw:
            continue
        child = Path(raw)
        metadata = read_json(child / "metadata.json")
        if metadata.get("pipeline_parent_run_dir") == str(run_dir) and metadata.get("status") in {"queued", "running", "pending"}:
            stop_job(child, reason="Parent pipeline was cancelled")


def enqueue_graph_pipeline(graph: dict, *, pipeline_name: str, owner_id: str, project_id: str, reuse_completed=True, gpu_device="all") -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import enqueue_workflow_call
    graph = validate_graph(graph, require_inputs=True)["graph"]
    if graph.get("template_id") == "unidock_active_learning":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.active_learning_controller import run_active_learning_pipeline
        from mn_molecular_design_studio.collection_policy import collection_scope
        with collection_scope("manual"):
            return enqueue_workflow_call(
                run_active_learning_pipeline, task_group="campaign", tool="visual_pipeline", job_type="visual_pipeline",
                kwargs={"graph": graph, "pipeline_name": pipeline_name, "owner_id": owner_id,
                        "project_id": project_id, "reuse_completed": False, "gpu_device": "none"},
                cpu_cores=1, gpu_device="none",
            )
    return enqueue_workflow_call(run_graph_pipeline, task_group="campaign", tool="visual_pipeline", job_type="visual_pipeline",
        kwargs={"graph": graph, "pipeline_name": pipeline_name, "owner_id": owner_id, "project_id": project_id, "reuse_completed": reuse_completed, "gpu_device": gpu_device}, cpu_cores=1, gpu_device="none")


def run_graph_pipeline(graph: dict, pipeline_name: str, owner_id: str, project_id: str, reuse_completed: bool = True, gpu_device: str = "all") -> Path:
    plan = validate_graph(graph, require_inputs=True)
    graph = plan["graph"]
    job = create_job("campaign", "visual_pipeline", "visual_pipeline", {}, {"pipeline_name": pipeline_name, "graph": graph, "reuse_completed": reuse_completed, "gpu_device": "none"})
    parent = job.run_dir
    from mn_molecular_design_studio.workspace import mark_studio_job
    mark_studio_job({"job": {"run_dir": str(parent)}}, owner_id=owner_id, project_id=project_id)
    snapshot = parent / "graph_snapshot.json"
    if snapshot.exists() and read_json(snapshot) != graph:
        raise ValueError("A submitted pipeline snapshot cannot be edited. Submit the changed draft as a new run.")
    if not snapshot.exists():
        with snapshot.open("x") as handle:
            json.dump(graph, handle, indent=2)
    state = read_json(parent / "graph_execution.json") or {"nodes": {}, "graph_id": graph["id"]}
    node_states, outputs = state["nodes"], {}
    nodes = {n["id"]: n for n in graph["nodes"]}
    cache_root = scope_root(owner_id, project_id) / "cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    cancelled = [False]
    runtime_identities = {}
    old_handler = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGTERM, lambda *_: cancelled.__setitem__(0, True))
    def persist():
        state["updated_at"] = utc_now(); write_json(parent / "graph_execution.json", state)
    def ancestors_of(node_id):
        pending, ancestors = list(plan["incoming"][node_id]), set()
        while pending:
            ancestor = pending.pop()
            if ancestor not in ancestors:
                ancestors.add(ancestor)
                pending.extend(plan["incoming"][ancestor])
        return ancestors
    def complete(node_id, output, key, status="completed"):
        from mn_molecular_design_studio.candidate_execution import inherit_partial_pool
        pool = inherit_partial_pool(output, [outputs[p] for p in plan['incoming'][node_id] if p in outputs],
                                    missing_branches=node_states.get(node_id, {}).get('missing_branches', []))
        if pool:
            output['partial_pool'] = pool
            node_states.setdefault(node_id, {})['partial'] = True
        from mn_molecular_design_studio.cache_provenance import referenced_files
        output["content_manifest"] = referenced_files({k: output.get(k) for k in
            ("candidates_jsonl", "review_candidates_jsonl", "target", "molecule_task", "analysis_runs", "analysis_run")}, base=output.get("run_dir"))
        if output.get("candidates_jsonl"):
            output["content_manifest"].update(referenced_files(_absolute_candidates(output), base=output.get("run_dir")))
        if output.get("review_candidates_jsonl"):
            output["review_digest"] = file_digest(Path(output["review_candidates_jsonl"]))
        output["fingerprint"] = fingerprint(key, {"candidate_digest": file_digest(Path(output["candidates_jsonl"])),
            "review_digest": output.get("review_digest"), "artifact_digests": sorted(output["content_manifest"].values()),
            **({'partial_pool': output['partial_pool']} if output.get('partial_pool') else {})}, []) if output.get("candidates_jsonl") else key
        outputs[node_id] = output
        node_states.setdefault(node_id, {}).update(status=status, output=output, fingerprint=key, candidate_count=output.get("candidate_count", 0), completed_at=utc_now())
        analyses = [*(output.get("analysis_runs") or []), *([output["analysis_run"]] if output.get("analysis_run") else [])]
        if (output.get("candidates_jsonl") or output.get("analysis_run", {}).get("success")) and all(a.get("success") for a in analyses) and status not in {"cached", "partial", "failed"} and not output.get('partial_pool') and not output.get('candidate_execution', {}).get('missing_batches'):
            write_json(cache_root / f"{key}.json", {"output": output})
        persist()
    try:
        update_status(parent, "running", campaign_name=pipeline_name, current_phase="Scheduling pipeline steps")
        while len(outputs) < len(nodes):
            if cancelled[0] or read_json(parent / "metadata.json").get("status") == "cancelled":
                cancel_graph_children(parent)
                update_status(parent, "cancelled", progress_label="Pipeline cancelled; completed results preserved")
                return parent
            progressed = False
            for node_id in plan["order"]:
                if node_id in outputs:
                    continue
                node, parents = nodes[node_id], plan["incoming"][node_id]
                current = node_states.setdefault(node_id, {"status": "pending"})
                from mn_molecular_design_studio import candidate_execution
                settings = candidate_execution.policy(node["data"].get("execution", graph.get("candidate_execution")))
                if current.get("status") in {"failed", "blocked"}:
                    continue
                if (settings['mode'] == 'sharded' and settings['incremental_cpu'] and len(parents) == 1
                        and candidate_execution.incremental_cpu_kind(node['data']['kind'])
                        and node_states.get(parents[0], {}).get('batches')):
                    upstream_state = node_states[parents[0]]
                    ready_batch = next((b for b in upstream_state['batches'] if b['status'] in {'completed', 'partial'}), None)
                    if ready_batch:
                        import sys
                        streaming_rows = json.loads(Path(ready_batch['rows_path']).read_text())
                        streaming_source = _local_source(parent, node_id + '_stream_probe', streaming_rows)
                        target = outputs.get(parents[0], {}).get('target') or next(
                            (outputs[a].get('target') for a in ancestors_of(parents[0]) if a in outputs and outputs[a].get('target')), {})
                        try:
                            streamed = candidate_execution.tick_incremental(parent, node, current, upstream_state,
                                {**streaming_source, 'target': target}, settings,
                                graph=sys.modules[__name__], owner_id=owner_id, project_id=project_id, persist=persist)
                            if streamed:
                                output, stream_key, stream_status = streamed
                                complete(node_id, output, stream_key, stream_status)
                                progressed = True
                        except Exception as exc:
                            candidate_execution.cancel(parent, {'nodes': {node_id: current}})
                            current.update(status='failed', error=str(exc)); persist(); progressed = True
                    if ready_batch or upstream_state.get('status') not in {'failed', 'blocked', 'cancelled'}:
                        continue
                failed_parents = [p for p in parents if node_states.get(p, {}).get("status") in {"failed", "blocked"}
                                  and (p not in outputs or node_states.get(p, {}).get("batches"))]
                partial_fan_in = bool(failed_parents and settings['failure_policy'] == 'partial')
                if failed_parents and not partial_fan_in:
                    current.update(status="blocked", error="An upstream step failed."); persist(); progressed=True; continue
                if not all(p in outputs or partial_fan_in and p in failed_parents for p in parents):
                    continue
                if partial_fan_in:
                    current['missing_branches'] = failed_parents
                    current['partial'] = True
                    parents = [p for p in parents if p in outputs and p not in failed_parents]
                    if not parents:
                        current.update(status='blocked', error='No successful upstream branches are available.')
                        persist(); progressed = True; continue
                try:
                    kind, params = node["data"]["kind"], node["data"]["params"]
                    upstream = [outputs[p] for p in parents]
                    source = upstream[0] if upstream else {}
                    if kind in {"ligandmpnn", "foundry_mpnn", "sabr", "abmpnn"} and len(upstream) > 1:
                        merged_rows = []
                        merged_target = next((item.get("target") for item in upstream if item.get("target")), {})
                        for source_id, output in zip(parents, upstream):
                            for row in _absolute_candidates(output):
                                candidate = copy.deepcopy(row)
                                original_id = str(candidate.get("candidate_id") or "candidate")
                                candidate["candidate_id"] = f"{source_id}:{original_id}"
                                raw = dict(candidate.get("raw_metadata") or {})
                                raw["pipeline_generator_node"] = source_id
                                raw["pipeline_original_candidate_id"] = original_id
                                candidate["raw_metadata"] = raw
                                merged_rows.append(candidate)
                        source = _local_source(parent, f"{node_id}_fan_in", merged_rows, merged_target)
                    key_params = copy.deepcopy(params)
                    if settings["mode"] != "native":
                        key_params["candidate_execution"] = settings
                    if kind in {"filter", "results", "ipsae", "rosetta_metrics"} or NODE_CATALOG[kind]["workflow"] and NODE_CATALOG[kind]["workflow"].startswith("analysis."):
                        from mn_molecular_design_studio import output_vocabulary, candidate_selection, measurement_status
                        key_params["scientific_results_revision"] = [file_digest(Path(module.__file__)) for module in
                            (output_vocabulary, candidate_selection, measurement_status)]
                    if kind == "target" or kind == 'cofold_input' and params.get('target_pdb') or kind=='complexa_generator' and params.get('target_pdb'):
                        target_path = _path(params["target_pdb"])
                        key_params["structure_digest"] = file_digest(target_path)
                    if kind in {"design_source","generation_source","backbone_source","molecule_source","project_design_source"}:
                        loader=(source_molecule_selection if kind == 'molecule_source' else source_project_designs if kind == 'project_design_source' else source_generations if kind in {'generation_source','backbone_source'} else source_designs)
                        extra={'stage':'generation.backbone' if kind=='backbone_source' else 'generation.backbone_sequence'} if kind in {'generation_source','backbone_source'} else {}
                        source_dir, _ = loader(params["source_run_dir"], project_id, params.get("selected_candidate_ids"),**extra)
                        key_params["source_digest"] = file_digest(source_dir / "artifacts/normalized_candidates/candidates.jsonl")
                        from mn_molecular_design_studio.cache_provenance import referenced_files
                        key_params["source_artifact_contents"] = referenced_files(_absolute_candidates({
                            "run_dir": str(source_dir), "candidates_jsonl": str(source_dir / "artifacts/normalized_candidates/candidates.jsonl")}), base=source_dir)
                        if kind == "project_design_source":
                            from mn_molecular_design_studio import mimosa_motif
                            key_params["mimosa_motif_adapter_digest"] = file_digest(Path(mimosa_motif.__file__))
                    if kind=='complexa_generator':
                        from mn_molecular_design_studio import complexa_input, complexa_container
                        key_params['input_adapter_digest']=file_digest(Path(complexa_input.__file__))
                        key_params['container_adapter_digest']=file_digest(Path(complexa_container.__file__))
                    from mn_molecular_design_studio import generator_input
                    from mn_molecular_design_studio import inverse_folding_input
                    if kind in inverse_folding_input.ENGINES:
                        from mn_molecular_design_studio.engine.workflow_engine.workflows import inverse_folding
                        key_params["input_adapter_digest"] = file_digest(Path(inverse_folding_input.__file__))
                        key_params["container_adapter_digest"] = file_digest(inverse_folding.NATIVE)
                    if kind in generator_input.ENGINES:
                        key_params['input_adapter_digest']=file_digest(Path(generator_input.__file__))
                        for key in ('target_pdb','scaffold_pdb','scaffold_ss_file','scaffold_adj_file','target_ss_file','target_adj_file'):
                            if params.get(key):key_params[key+'_digest']=file_digest(_path(params[key]))
                    elif kind == "peptide_trap":
                        from mn_molecular_design_studio import peptide_trap_input
                        from mn_molecular_design_studio.engine.workflow_engine.workflows import peptide_trap
                        prepared = peptide_trap_input.prepare_request(params)
                        key_params["input_adapter_digest"] = file_digest(Path(peptide_trap_input.__file__))
                        key_params["workflow_adapter_digest"] = file_digest(Path(peptide_trap.__file__))
                        key_params["template_digests"] = [file_digest(peptide_trap_input.template_path(key)) for key in prepared["scaffold_ids"]]
                        key_params["native_protocol_digest"] = file_digest(peptide_trap_input.library_root() / "source/design_after_seqthread.xml")
                        key_params["container_adapter_digest"] = file_digest(Path(__file__).parents[1] / "tools/peptide_trap_runtime.py")
                    elif kind == "kejia_idr":
                        from mn_molecular_design_studio import kejia_idr_input
                        key_params["input_adapter_digest"] = file_digest(Path(kejia_idr_input.__file__))
                        key_params["source_revision"] = kejia_idr_input.SOURCE_REVISION
                        prepared = kejia_idr_input.prepare_request(params)
                        key_params["template_digests"] = [file_digest(Path(row["path"])) for row in prepared["templates"]]
                        key_params["container_adapter_digest"] = file_digest(Path(__file__).parents[1] / "tools/kejia_idr_threading.py")
                    elif kind == "rfdiffusion_beta_pairing":
                        from mn_molecular_design_studio import rfdiffusion_beta_pairing_input
                        key_params["input_adapter_digest"] = file_digest(Path(rfdiffusion_beta_pairing_input.__file__))
                        if params.get("target_pdb"):
                            key_params["target_pdb_digest"] = file_digest(_path(params["target_pdb"]))
                    elif kind == "sabr":
                        from mn_molecular_design_studio import sabr_runtime
                        key_params["implementation_digest"] = file_digest(Path(sabr_runtime.__file__))
                    elif kind == "abmpnn":
                        from mn_molecular_design_studio import mpnn_input
                        from mn_molecular_design_studio import mimosa_motif
                        key_params["checkpoint_digest"] = file_digest(mpnn_input.checkpoint("ab_mpnn"))
                        key_params["mimosa_motif_adapter_digest"] = file_digest(Path(mimosa_motif.__file__))
                    if kind in {"af2", "boltz2"}:
                        from mn_molecular_design_studio import prediction_configuration
                        key_params["prediction_configuration_digest"] = file_digest(Path(prediction_configuration.__file__))
                        if kind == "af2":
                            from mn_molecular_design_studio.engine.workflow_engine.workflows.refolding import AF2_BINDER_EVAL
                            key_params["native_runner_digest"] = file_digest(AF2_BINDER_EVAL)
                    if NODE_CATALOG[kind]["workflow"]:
                        from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function
                        effective_request = node_request(kind, params, source, "cache-specification")
                        implementation = _get_function(WORKFLOWS[effective_request["workflow"]])
                        key_params["implementation_digest"] = file_digest(Path(inspect.getfile(implementation)))
                        from mn_molecular_design_studio.cache_provenance import execution_identity
                        effective = effective_request["parameters"]
                        identity_key = json.dumps([kind, params], sort_keys=True)
                        if identity_key not in runtime_identities:
                            runtime_identities[identity_key] = execution_identity(kind, effective, effective_request["workflow"])
                        key_params["runtime_identity"] = runtime_identities[identity_key]
                        # Cache the effective request, excluding run-local transport
                        # paths whose scientific contents are covered by upstream keys.
                        key_params["effective_parameters"] = {k: v for k, v in effective.items() if k not in
                            {"source_run_dir", "candidates_jsonl", "job_name", "campaign_name", "target_pdb", "target_structure"}}
                        if "seed" in effective:
                            seed_key = fingerprint(kind, key_params, [o["fingerprint"] for o in upstream])
                            seed_file = cache_root / (seed_key + ".seed.json")
                            resolved_seed = current.get("resolved_seed") or effective.get("seed")
                            if not resolved_seed:
                                resolved_seed = (read_json(seed_file).get("seed") if reuse_completed else None) or secrets.randbelow(2147483646) + 1
                            current["resolved_seed"] = resolved_seed
                            if reuse_completed:
                                write_json(seed_file, {"seed": resolved_seed})
                            params = {**params, "seed": resolved_seed}
                            key_params["effective_parameters"]["seed"] = resolved_seed
                        current["execution_specification"] = key_params
                    key = fingerprint(kind, key_params, [o["fingerprint"] for o in upstream])
                    if settings['mode'] == 'sharded' and candidate_execution.supports(kind, source):
                        cached = read_json(cache_root / f'{key}.json').get('output') if reuse_completed and not current.get('batches') else None
                        if cached and _cache_valid(cached):
                            complete(node_id, cached, key, 'cached')
                            progressed = True
                            continue
                        if current.get('batches') and current.get('fingerprint') != key:
                            raise ValueError('Candidate inputs changed during execution. Submit a new pipeline.')
                        import sys
                        result = candidate_execution.tick(parent, node, current, source, params, key, settings,
                            graph=sys.modules[__name__], owner_id=owner_id, project_id=project_id,
                            gpu_device=gpu_device, persist=persist)
                        if result:
                            output, batch_status = result
                            complete(node_id, output, key, batch_status)
                            progressed = True
                        continue
                    if current.get("job_run_dir") and (current.get("fingerprint") != key or current.get("output") and not _cache_valid(current["output"])):
                        previous = {k: current.get(k) for k in ("job_run_dir", "fingerprint", "status", "resolved_seed")}
                        current.setdefault("superseded_attempts", []).append(previous)
                        child = Path(current.pop("job_run_dir"))
                        if read_json(child / "metadata.json").get("status") not in {"completed", "failed", "cancelled"}:
                            from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
                            stop_job(child, "Pipeline scientific inputs changed; superseded by a fresh execution")
                        current.update(status="pending")
                        current.pop("output", None)
                    if current.get("job_run_dir"):
                        child = Path(current["job_run_dir"])
                        metadata = read_json(child / "metadata.json")
                        current["status"] = str(metadata.get("status") or "queued")
                        current["job_code"] = metadata.get("job_code")
                        child_result = read_json(child / "result.json")
                        is_analysis_node = bool(NODE_CATALOG[kind]["outputs"] and NODE_CATALOG[kind]["outputs"][0]["type"] == "analysis")
                        if is_analysis_node and current["status"] in {"completed", "failed"} and isinstance(child_result.get("outputs", {}).get("tool_status"), dict):
                            output = _child_analysis_output(child, source.get("target", {}))
                            step_status = "completed" if output["analysis_run"]["success"] else "failed"
                            complete(node_id, output, key, step_status)
                            if step_status == "failed":
                                node_states[node_id]["error"] = "One or more surface-analysis methods failed; partial results are available."
                                persist()
                            progressed=True
                        elif current["status"] == "completed" and child_result.get("success") is True:
                            if is_analysis_node:
                                current.update(status="failed", error="The surface-analysis job completed without reporting method results.")
                                persist();progressed=True
                            else:
                                output = _child_output(child, source.get("target", {}))
                                if NODE_CATALOG[kind]["category"] == "Structure prediction":
                                    output = _attach_refold_input_references(output, source, parent, node_id)
                                complete(node_id, output, key); progressed=True
                        elif current["status"] == "failed" and child_result.get("scientific_status") == "failed" and (child / "artifacts/normalized_candidates/candidates.jsonl").is_file():
                            output = _child_output(child, source.get("target", {}))
                            complete(node_id, output, key, "failed")
                            node_states[node_id]["error"] = "No valid required measurements; per-candidate causes are retained."
                            persist(); progressed = True
                        elif current["status"] in {"failed", "cancelled"}:
                            current.update(status="failed", error=f"Child job {current.get('job_code') or child.name} did not complete successfully.");persist();progressed=True
                        continue
                    cached = read_json(cache_root / f"{key}.json").get("output") if reuse_completed else None
                    if cached and _cache_valid(cached):
                        if NODE_CATALOG[kind]["category"] == "Structure prediction":
                            cached = _attach_refold_input_references(cached, source, parent, node_id)
                        complete(node_id, cached, key, "cached");progressed=True;continue
                    if kind == "target":
                        chains = [x for x in re.split(r"[\s,;]+", params["target_chains"]) if x]
                        if not chains:
                            raise ValueError("Select at least one target chain.")
                        staged_target = parent / "artifacts" / "pipeline_inputs" / node_id / ("target-" + key_params["structure_digest"][:16] + target_path.suffix)
                        staged_target.parent.mkdir(parents=True, exist_ok=True)
                        if not staged_target.exists():
                            shutil.copy2(target_path, staged_target)
                        output = {"target": {"target_pdb": str(staged_target), "target_chains": chains,
                                             "hotspots": params["hotspots"], "coldspots": params.get("coldspots", "")}}
                        complete(node_id, output, key);progressed=True;continue
                    if kind == 'cofold_input':
                        task = {key: copy.deepcopy(params.get(key, '')) for key in (
                            'task_type', 'target_protein_sequences', 'partner_protein_sequences',
                            'dna_sequences', 'rna_sequences', 'ligand_input_type', 'ligand_ccd_codes', 'ligand_smiles', 'compound_set_file', 'sequence_file',
                            'cyclic_partner', 'target_chain_id',
                        )}
                        target_pdb = str(params.get('target_pdb') or '').strip()
                        task['target_pdb'] = str(_path(target_pdb)) if target_pdb else ''
                        complete(node_id, {'molecule_task': task}, key);progressed=True;continue
                    if kind == "binder_context":
                        target = copy.deepcopy(source.get("target") or {})
                        if not target.get("target_pdb") or not target.get("target_chains"):
                            raise ValueError("Binder Design Context needs a connected Target Structure with selected chains.")
                        if params.get("clear_unselected"):
                            target["hotspots"] = params.get("hotspots", "")
                            target["coldspots"] = params.get("coldspots", "")
                        else:
                            target["hotspots"] = params.get("hotspots") or target.get("hotspots", "")
                            target["coldspots"] = params.get("coldspots") or target.get("coldspots", "")
                        complete(node_id, {"target": target}, key);progressed=True;continue
                    if kind in {"design_source","generation_source","backbone_source","molecule_source","project_design_source"}:
                        loader=(source_molecule_selection if kind == 'molecule_source' else source_project_designs if kind == 'project_design_source' else source_generations if kind in {'generation_source','backbone_source'} else source_designs)
                        extra={'stage':'generation.backbone' if kind=='backbone_source' else 'generation.backbone_sequence'} if kind in {'generation_source','backbone_source'} else {}
                        source_dir, selected = loader(params["source_run_dir"], project_id, params.get("selected_candidate_ids"),**extra)
                        all_rows = _absolute_candidates({"run_dir": str(source_dir), "candidates_jsonl": str(source_dir / "artifacts/normalized_candidates/candidates.jsonl")})
                        selected_rows = _selected_source_rows(kind, selected, all_rows)
                        if kind in {'generation_source','backbone_source'}:
                            source_input = read_json(source_dir/'input.json')
                            source_params = source_input.get('params', {}) if isinstance(source_input.get('params'), dict) else {}
                            source_inputs = source_input.get('inputs', {}) if isinstance(source_input.get('inputs'), dict) else {}
                            target_inputs = dict(source_inputs)
                            if not target_inputs.get('target_pdb'):
                                target_inputs['target_pdb'] = target_inputs.get('target_structure')
                            target = _generation_target({**target_inputs, **source_params})
                            # Candidate paths have already been resolved relative to the
                            # source run, so prefer them over external/original upload paths.
                            candidate_target = next((row.get('target_pdb') for row in all_rows if row.get('target_pdb')), None)
                            if candidate_target:
                                target['target_pdb'] = candidate_target
                            elif target.get('target_pdb'):
                                resolved_target = resolve_stored_path(target['target_pdb'], run_dir=source_dir)
                                if resolved_target:
                                    target['target_pdb'] = str(resolved_target)
                            if not target.get('target_chains'):
                                target['target_chains'] = next((row.get('target_chains') for row in all_rows if row.get('target_chains')), [])
                            if not target.get('hotspots'):
                                target['hotspots'] = next((row.get('hotspots') for row in all_rows if row.get('hotspots')), [])
                        elif kind == 'project_design_source':
                            first = selected_rows[0] if selected_rows else {}
                            target = {
                                "target_pdb": first.get("target_pdb", ""),
                                "target_chains": first.get("target_chains") or [],
                                "hotspots": first.get("hotspots", ""),
                            }
                        else:
                            target = {}
                        complete(node_id, _local_source(parent, node_id, selected_rows,target), key);progressed=True;continue
                    if kind in {"filter", "results"}:
                        rows = []
                        analysis_runs = []
                        for source_id, output in zip(parents, upstream):
                            analysis_runs.extend(output.get("analysis_runs") or [])
                            if output.get("analysis_run"):
                                analysis_runs.append(output["analysis_run"])
                            review_source = {**output, "candidates_jsonl": output["review_candidates_jsonl"]} if kind == "results" and output.get("review_candidates_jsonl") else output
                            for row in _absolute_candidates(review_source) if review_source.get("candidates_jsonl") else []:
                                if kind == "results":
                                    row = _add_target_aligned_pose_rmsd(row, parent)
                                    row.setdefault("raw_metadata", {}).update(pipeline_source_node=source_id, pipeline_original_candidate_id=row["candidate_id"])
                                    row["candidate_id"] = f"{source_id}:{row['candidate_id']}"
                                rows.append(row)
                        if kind == "filter":
                            from mn_molecular_design_studio.candidate_selection import evaluate_selection, selected_candidates, write_selection_artifacts
                            outcomes = evaluate_selection(rows, params)
                            write_selection_artifacts(parent / "artifacts/pipeline_nodes" / node_id / "artifacts/selection", outcomes)
                            rows = selected_candidates(outcomes)
                            if outcomes and all(row["raw_metadata"]["selection"]["status"] == "not_evaluated" for row in outcomes):
                                rows = outcomes[:params["limit"]] if params.get("limit") else outcomes
                        missing = sorted(set(current.get('missing_branches', [])) | {
                            missing for upstream_output in upstream for missing in upstream_output.get('partial_pool', {}).get('missing_branches', [])})
                        partial_batches = [upstream_output['candidate_execution'] for upstream_output in upstream
                                           if upstream_output.get('candidate_execution', {}).get('missing_batches')]
                        partial_batches.extend(record for upstream_output in upstream
                                               for record in upstream_output.get('partial_pool', {}).get('candidate_batches', []))
                        if missing or partial_batches:
                            for row in rows:
                                row.setdefault('raw_metadata', {})['partial_pool'] = dict(missing_branches=missing, candidate_batches=partial_batches)
                        output = _local_source(parent, node_id, rows, source.get("target"), analysis_runs if kind == "results" else [], final_evidence=kind == "results")
                        if missing or partial_batches:
                            output['partial_pool'] = dict(missing_branches=missing, candidate_batches=partial_batches)
                            current['partial'] = True
                        if kind == "filter":
                            output["review_candidates_jsonl"] = str(parent / "artifacts/pipeline_nodes" / node_id / "artifacts/selection/all_candidates.jsonl")
                        complete(node_id, output, key);progressed=True;continue
                    if kind == "sabr":
                        from mn_molecular_design_studio import sabr_runtime
                        from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
                        rows = _absolute_candidates(source) if source.get("candidates_jsonl") else []
                        if not rows:
                            raise ValueError("SAbR needs at least one connected antibody design with a structure.")
                        for row in rows:
                            path_text = row.get("complex_pdb") or row.get("binder_pdb")
                            path = resolve_stored_path(path_text, run_dir=_path(source["run_dir"]), must_exist=True)
                            if path is None:
                                raise ValueError(f"SAbR could not locate a structure for {row.get('candidate_id', 'candidate')}.")
                            chains = [str(chain) for chain in row.get("binder_chains") or []]
                            if not chains:
                                raise ValueError(f"SAbR could not identify the antibody chain for {row.get('candidate_id', 'candidate')}.")
                            annotation = sabr_runtime.annotate_file(
                                path, chains, chain_type=params.get("chain_type", "auto"),
                                noise_level=float(params.get("noise_level", 0.0)),
                                mode=params.get("mode", "sabr"), scheme="imgt",
                            )
                            residues_by_region = {region: [] for region in sabr_runtime.CDR_NAMES}
                            for chain_data in annotation["chains"].values():
                                for region, labels in chain_data["cdr_residues"].items():
                                    residues_by_region[region].extend(labels)
                            raw = dict(row.get("raw_metadata") or {})
                            raw["sabr_cdr_selection"] = {
                                "chain_ids": chains,
                                "regions": list(sabr_runtime.CDR_NAMES),
                                "residues_by_region": residues_by_region,
                                "chain_type": params.get("chain_type", "auto"),
                                "mode": params.get("mode", "sabr"),
                                "noise_level": float(params.get("noise_level", 0.0)),
                                "scheme": "imgt",
                            }
                            raw["sabr_annotation"] = annotation
                            row["raw_metadata"] = raw
                        complete(node_id, _local_source(parent, node_id, rows, source.get("target")), key)
                        progressed = True
                        continue
                    if source.get("candidates_jsonl") and not source.get("candidate_count"):
                        complete(node_id, _local_source(parent, node_id, [], source.get("target")), key, "skipped");progressed=True;continue
                    request = node_request(kind, params, source, f"{pipeline_name} · {node['data']['label']}")
                    child = _submit(request, parent, node, owner_id=owner_id, project_id=project_id, gpu_device=gpu_device)
                    metadata = read_json(child / "metadata.json")
                    current.update(status="queued", job_run_dir=str(child), job_code=metadata.get("job_code"), fingerprint=key, started_at=utc_now())
                    persist();progressed=True
                except Exception as exc:
                    if current.get('batches'):
                        candidate_execution.cancel(parent, {'nodes': {node_id: current}})
                    current.update(status="failed", error=str(exc));persist();progressed=True
            unresolved = [n for n in nodes if n not in outputs and node_states.get(n, {}).get("status") not in {"failed", "blocked"}]
            if not unresolved:
                break
            if not progressed:
                persist();time.sleep(2)
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
        candidates=[]
        surface_analyses=[]
        seen_surface_runs=set()
        for node_id, output in outputs.items():
            analyses = [*(output.get("analysis_runs") or []), *([output["analysis_run"]] if output.get("analysis_run") else [])]
            for analysis in analyses:
                identity = (str(analysis.get("task_group") or "detection"), str(analysis.get("run_id") or analysis.get("run_dir") or ""))
                if identity not in seen_surface_runs:
                    seen_surface_runs.add(identity)
                    surface_analyses.append(analysis)
            if nodes[node_id]["data"]["kind"] != "results":
                continue
            ancestors = ancestors_of(node_id)
            for candidate in _absolute_candidates(output) if output.get("candidates_jsonl") else []:
                candidate["candidate_id"] = f"{node_id}:{candidate['candidate_id']}"
                candidate.setdefault("raw_metadata", {})["pipeline_result_node"] = node_id
                if graph.get("task_configuration"):
                    candidate["raw_metadata"]["pipeline_task"] = graph["task_configuration"].get("task_id")
                    # Record actual graph ancestry, including edited recipes,
                    # rather than infer provenance from menu labels/settings.
                    candidate["raw_metadata"]["pipeline_method_nodes"] = [
                        {"node_id": ancestor, "kind": nodes[ancestor]["data"]["kind"],
                         "parameters": nodes[ancestor]["data"]["params"]}
                        for ancestor in sorted(ancestors)
                        if NODE_CATALOG[nodes[ancestor]["data"]["kind"]]["category"] in
                        {"Design engine", "Design pipeline", "Structure prediction", "Sequence design"}]
                candidates.append(candidate)
        write_candidates(parent, "visual_pipeline", candidates)
        from mn_molecular_design_studio.output_vocabulary import evidence_rows
        import pandas as pd
        pd.DataFrame(evidence_rows(candidates)).to_csv(parent / "artifacts/normalized_candidates/final_evidence.csv", index=False)
        failures = [v for v in node_states.values() if v.get("status") in {"failed", "blocked"}]
        # Only new, explicitly configured comparisons permit partial success.
        # Legacy graphs retain their original failure semantics. Each comparison
        # branch owns a Results node, so a failed branch has no fan-in barrier.
        comparison = graph.get("task_configuration", {}).get("comparison") is True
        successful_results = [node_id for node_id in outputs if nodes[node_id]["data"]["kind"] == "results"
                              and node_states.get(node_id, {}).get("status") in {"completed", "cached"}
                              and (node_states[node_id].get('partial') or not any(node_states.get(a, {}).get("status") in {"failed", "blocked"}
                                          for a in ancestors_of(node_id)))]
        explicit_partial_results = [n for n in successful_results if node_states[n].get('partial')]
        success = not failures or bool(comparison or explicit_partial_results) and bool(successful_results)
        comparison_status = "partial" if (failures or any(n.get('partial') for n in node_states.values())) and success else "failed" if failures else "completed"
        state["status"] = "completed" if success else "failed"
        state["comparison_status"] = comparison_status
        persist()
        finish_job(parent, success, {"outputs": {"candidates": candidates, "surface_analyses": surface_analyses, "pipeline_nodes": node_states,
                    "comparison_status": comparison_status, "successful_result_nodes": successful_results}, "metrics": {"candidate_count": len(candidates), "surface_analysis_count": len(surface_analyses), "node_count": len(nodes), "cached_steps": sum(v.get("status")=="cached" for v in node_states.values()), "failed_steps": len(failures)}, "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl", "graph_snapshot": "graph_snapshot.json", "graph_execution": "graph_execution.json"}})
        return parent
    except Exception:
        cancel_graph_children(parent)
        raise
    finally:
        signal.signal(signal.SIGTERM, old_handler)
