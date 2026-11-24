"""Canvas adapters for standalone tool workflows with typed inputs and results."""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlsplit


@dataclass(frozen=True)
class TileWorkflow:
    label: str
    workflow: str
    category: str
    output: str = "tool_result"
    inputs: tuple[tuple[str, str], ...] = ()
    cpu: bool = False


TARGET = (("target", "target"),)
CANDIDATES = (("candidates", "candidates"),)
COMPOUNDS = (("compounds", "compound_set"),)
REFERENCE_LIGAND = (("reference_ligand", "ligand_sdf"),)
LIGAND = (("ligand", "ligand_sdf"),)
TILES = {
    **{
        "tile_docking_" + engine: TileWorkflow(label, "ligand.docking." + engine, "Docking",
            "docking_result", TARGET + COMPOUNDS + (REFERENCE_LIGAND if engine == "rosetta_ligand" else ()),
            cpu=engine in {"vina", "rosetta_ligand"})
        for engine, label in (
            ("vina", "AutoDock Vina"), ("autodock_gpu", "AutoDock-GPU"), ("gnina", "GNINA"),
            ("rosetta_ligand", "RosettaLigand"), ("unidock_pro", "Uni-Dock Pro"),
        )
    },
    **{
        "tile_generation_" + engine: TileWorkflow(label, "ligand.generation." + engine,
            "Ligand generation", "compound_set", TARGET + (REFERENCE_LIGAND if engine == "lddm" else ()))
        for engine, label in (
            ("conditar", "conDitar"), ("drugrpg", "DrugRPG"), ("flowr_root", "FLOWR.root"),
            ("lddm", "LDDM"), ("omtra", "OMTRA"), ("paopt", "paOPT"), ("pfm", "PFM"),
            ("pgmg", "PGMG"), ("pocketflow", "PocketFlow"), ("pocketxmol", "PocketXMol"),
        )
    },
    **{
        "tile_detection_" + tool: TileWorkflow(label, "detection." + tool, "Analysis",
            "analysis", TARGET, cpu=tool in {"masif_seed", "scannet"})
        for tool, label in (
            ("masif_seed", "MaSIF-seed"), ("pesto", "PeSTo"), ("scannet", "ScanNet"), ("surf2spot", "Surf2Spot"),
        )
    },
    **{
        "tile_design_" + tool: TileWorkflow(label, "design." + tool, "Antibody and peptide design", "candidates", TARGET)
        for tool, label in (
            ("rfantibody", "RFantibody"), ("rfantibody_mimosa", "MIMOSA · RFantibody"),
            ("boltzgen_mimosa", "MIMOSA · BoltzGen"), ("rfpeptides", "RFpeptides"),
        )
    },
    "tile_affinity_nesso": TileWorkflow("Nesso-1", "ligand.affinity.nesso", "Docking", inputs=(*TARGET, *COMPOUNDS)),
    "tile_prepare_library": TileWorkflow("Prepare library", "library.prepare", "Ligand library",
                                        "fingerprint_manifest", cpu=True),
    "tile_library_centroids": TileWorkflow("Create centroids", "library.centroids", "Ligand library",
                                          "centroid_manifest", (("features", "fingerprint_manifest"),), cpu=True),
    "tile_complex_prepare": TileWorkflow("Complex Preparation", "ligand.complex.prepare", "Molecular dynamics",
                                         "prepared_complex", cpu=True),
    "tile_md_openmm": TileWorkflow("OpenMM", "ligand.md.openmm", "Molecular dynamics", "md_trajectory",
                                   (("complex", "prepared_complex"),) + LIGAND),
    "tile_md_gromacs": TileWorkflow("GROMACS", "ligand.md.gromacs", "Molecular dynamics", "md_trajectory",
                                    (("complex", "prepared_complex"),) + LIGAND),
    "tile_bioemu2": TileWorkflow("BioEmu-2", "simulation.bioemu2", "Molecular dynamics", inputs=CANDIDATES),
    "tile_af3_jax": TileWorkflow("ColabFold 2.0 preview", "refolding.af3_jax_cofolding", "Structure prediction",
                                 "prediction"),
    "tile_binderscope": TileWorkflow("BinderScope", "analysis.binderscope", "Analysis", cpu=True),
    "tile_validation_acceptance": TileWorkflow("Validation & Acceptance", "analysis.validation_acceptance",
                                               "Analysis", "candidates", CANDIDATES, cpu=True),
}
WORKFLOWS = {tile.workflow: tile for tile in TILES.values()}
ARTIFACT_OUTPUTS = {
    "compound_set": "generated_compound_set",
    "fingerprint_manifest": "feature_manifest",
    "centroid_manifest": "centroid_manifest",
    "prepared_complex": "prepared_complex",
}
PORT_LABELS = {
    "target": "Target structure", "candidates": "Candidate collection", "compound_set": "Compound set",
    "fingerprint_manifest": "Fingerprint manifest", "centroid_manifest": "Centroid manifest",
    "prepared_complex": "Prepared complex", "md_trajectory": "Trajectory + topology", "tool_result": "Saved Results",
    "analysis": "Analysis results", "prediction": "Predicted structures", "docking_result": "Docking results",
    "ligand_sdf": "Reference ligand (SDF)",
}
ARTIFACT_INPUT_FIELDS = {"compounds": "compound_set_file", "features": "feature_manifest",
                         "complex": "complex_pdb", "reference_ligand": "reference_ligand_sdf",
                         "ligand": "ligand_sdf"}


@lru_cache(maxsize=None)
def parameter_descriptors(workflow: str) -> tuple[dict, ...]:
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_catalog

    return tuple(workflow_catalog(workflow_id=workflow)[0]["parameters"])


@lru_cache(maxsize=None)
def parameters_for(workflow: str) -> tuple[tuple[str, bool], ...]:
    return tuple((p["name"], p["required"]) for p in parameter_descriptors(workflow))


def connected_parameter_names(tile: TileWorkflow, port: str) -> set[str]:
    """Fields supplied by a connection without a manual override."""
    names = {
        "target": {"target_pdb", "target_chains", "chain_ids"},
        "candidates": {"source_run_dir", "candidates_jsonl"},
        "features": {"feature_manifest", "library_id"},
    }.get(port, {ARTIFACT_INPUT_FIELDS.get(port)})
    return names.intersection(dict(parameters_for(tile.workflow)))


def validate_parameters(tile: TileWorkflow, settings: dict, connected_ports: set[str], *, runnable: bool) -> None:
    parameters = dict(parameters_for(tile.workflow))
    unknown = set(settings) - parameters.keys()
    if unknown:
        raise ValueError(f"{tile.label} has no settings named {', '.join(sorted(unknown))}.")
    supplied = {
        name for name, value in settings.items()
        if value is not None and (not isinstance(value, str) or value.strip())
    }
    for port in connected_ports:
        supplied.update(connected_parameter_names(tile, port))
    if "target" in connected_ports:
        supplied.add("chain_id")  # A single connected chain can fill this optional override.
    missing = [name for name, required in parameters.items() if required and name not in supplied]
    if runnable and missing:
        raise ValueError(f"{tile.label} needs {', '.join(missing)}.")
    if runnable and tile.workflow == "analysis.binderscope":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.binderscope import validate_request

        if not settings.get("native_run_dirs"):
            raise ValueError("Select native BindCraft campaigns and their input target for this canvas step. "
                             "Analyze saved database designs from Data Browser.")
        validate_request(settings, require_files=False)


def workflow_request(tile: TileWorkflow, settings: dict, inputs: dict, name: str) -> dict:
    parameters = dict(settings)
    accepted = dict(parameters_for(tile.workflow))
    for handle, source in inputs.items():
        if not source or not isinstance(source, dict):
            raise ValueError(f"{tile.label} received an empty input on {handle}.")
        if handle == "target":
            target = source.get("target") or {}
            if not target.get("target_pdb"):
                raise ValueError(f"{tile.label} requires a saved target structure.")
            parameters["target_pdb"] = target["target_pdb"]
            chains = target.get("target_chains") or []
            for key in ("target_chains", "chain_ids"):
                if key in accepted and chains:
                    parameters[key] = chains
            if "chain_id" in accepted and "chain_id" not in parameters and len(chains) == 1:
                parameters["chain_id"] = chains[0]
            if "hotspots" in accepted and target.get("hotspots") and not parameters.get("hotspots"):
                parameters["hotspots"] = target["hotspots"]
        elif handle == "candidates":
            if source.get("run_sources") or not source.get("candidates_jsonl"):
                raise ValueError("Connect a materialized, local candidate collection for this tool.")
            parameters.update(source_run_dir=source["run_dir"], candidates_jsonl=source["candidates_jsonl"])
        else:
            artifact = source.get("tile_artifact") or {}
            if not artifact.get("path"):
                raise ValueError(f"The upstream workflow produced no {handle} artifact.")
            key = ARTIFACT_INPUT_FIELDS[handle]
            parameters[key] = artifact["path"]
            if handle == "features" and artifact.get("library_id"):
                parameters["library_id"] = artifact["library_id"]
    if "job_name" in accepted:
        parameters.setdefault("job_name", name)
    validate_parameters(tile, parameters, set(), runnable=True)
    return {"schema_version": 1, "workflow": tile.workflow, "parameters": parameters}


def artifact_output(tile: TileWorkflow, child: Path, result: dict) -> dict:
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

    key = ARTIFACT_OUTPUTS.get(tile.output)
    if key is None:
        return {}
    downstream = result.get("downstream_artifacts") or {}
    outputs = result.get("outputs") or {}
    raw = downstream.get(key) or outputs.get(key)
    if not isinstance(raw, str) or not raw:
        raise ValueError(f"{tile.label} saved no {PORT_LABELS[tile.output].lower()}.")
    resolved = resolve_stored_path(raw, run_dir=child, must_exist=True)
    if resolved is None or not resolved.is_file():
        raise ValueError(f"The {PORT_LABELS[tile.output].lower()} from {child.name} is unavailable.")
    library_id = outputs.get("library_id") or downstream.get("library_id") or (result.get("metrics") or {}).get("library_id")
    return {"tile_artifact": {"path": str(resolved), **({"library_id": library_id} if library_id else {})}}


def input_digests(parameters: dict, run_dir: Path) -> dict[str, str]:
    """Hash explicit file inputs through the same keys as submission path validation."""
    from mn_molecular_design_studio.cache_provenance import file_digest
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import (
        is_user_input_path_key, resolve_stored_path,
    )

    digests = {}
    for key, value in parameters.items():
        if not is_user_input_path_key(key):
            continue
        for item in value if isinstance(value, list) else [value]:
            if not isinstance(item, str) or not item:
                continue
            if item.startswith("cluster://"):
                from mn_molecular_design_studio.cluster import config

                # The authorized Run snapshot identifies peer-owned files.
                # Only the source workstation can resolve them to local paths.
                if urlsplit(item).netloc != config.identity():
                    continue
            path = resolve_stored_path(item, run_dir=run_dir, must_exist=True)
            if path is not None and path.is_file():
                digests[str(path)] = file_digest(path)
    return digests
