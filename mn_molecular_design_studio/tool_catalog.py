"""Shared definitions for the All tools tiles and workflow discovery."""
from __future__ import annotations

from mn_molecular_design_studio import generator_input
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction
from mn_molecular_design_studio.pipeline_graph import WORKFLOW_ROLES
from mn_molecular_design_studio.quality_tools import WORKFLOWS as QUALITY_WORKFLOWS
from mn_molecular_design_studio.tool_references import references_for_tool


WORKFLOW_TITLES = {
    "Bioemu2": "BioEmu-2",
    "Esm If1": "ESM-IF1",
    "Caliby": "Caliby",
    "Boltzgen Inverse": "BoltzGen Inverse Folding",
    "Pxdesign Generate": "PXDesign",
    "Rfdiffusion Generate": "RFdiffusion",
    "Rfdiffusion Beta Pairing": "RFdiffusion β-pairing",
    "Kejia Idr": "Kejia IDR Binder",
    "Peptide Trap": "Peptide Trap",
    "Genie3 Generate": "Genie3",
    "Protpardelle Generate": "Protpardelle-1C",
    "Proteina Complexa": "Proteina-Complexa",
    "Proteina Complexa Generate": "Proteina-Complexa Generator",
    "Bindcraft2": "BindCraft2",
    "Bindcraft2 Generator": "BindCraft2 Generator",
    "Rfantibody": "RFantibody",
    "Rfantibody Mimosa": "MIMOSA · RFantibody",
    "Rfpeptides": "RFpeptides",
    "Bindcraft": "BindCraft",
    "Boltzgen": "BoltzGen",
    "Boltzgen Mimosa": "MIMOSA · BoltzGen",
    "Boltzgen Generator": "BoltzGen Generator",
    "Af2 Initial Guess": "AlphaFold2 Initial Guess",
    "Esmfold2": "ESMFold2",
    "Esmfold2 Cofolding": "ESMFold2 Cofolding",
    "Esmfold Monomer": "ESMFold",
    "Af3 Jax Cofolding": "ColabFold 2.0 preview",
    "Rf3": "RosettaFold3",
    "Optimized Prediction": "Structure Prediction",
    "Esmfold2 Binder": "ESMFold2 Binder Design",
    "Boltz2 Complex": "Boltz-2 Complex",
    "Boltz2 Monomer": "Boltz-2 Monomer",
    "Ligandmpnn": "LigandMPNN",
    "Mpnn": "ProteinMPNN",
    "Abmpnn": "AbMPNN",
    "Sabr": "SAbR Antibody Numbering",
    "Rfdiffusion Classic": "RFdiffusion",
    "Rfdiffusion3 Foundry": "RFdiffusion3",
    "Rfdiffusion3": "RFdiffusion3",
    "Rfd4 Proteina": "RFD4-Proteina",
    "Ipsae": "ipSAE & Structure Analysis",
    "Rosetta Metrics": "Rosetta Interface Metrics",
    "Pesto": "PeSTo",
    "Scannet": "ScanNet",
    "Masif Seed": "MaSIF-seed",
    "Surf2Spot": "Surf2Spot",
    "Surface Analysis": "Protein Surface Analysis",
    "Docking.Vina": "AutoDock Vina",
    "Docking.Autodock Gpu": "AutoDock-GPU",
    "Docking.Gnina": "GNINA",
    "Docking.Unidock Pro": "Uni-Dock Pro",
    "Docking.Rosetta Ligand": "RosettaLigand",
    "Affinity.Nesso": "Nesso-1",
    "Complex.Prepare": "Complex Preparation",
    "Generation.Lddm": "LDDM",
    "Generation.Pocketxmol": "PocketXMol",
    "Generation.Omtra": "OMTRA",
    "Generation.Flowr Root": "FLOWR.root",
    "Generation.Conditar": "conDitar",
    "Generation.Paopt": "paOPT",
    "Generation.Drugrpg": "DrugRPG",
    "Generation.Pfm": "PFM",
    "Generation.Pocketflow": "PocketFlow",
    "Generation.Pgmg": "PGMG",
    "Md.Gromacs": "GROMACS",
    "Md.Openmm": "OpenMM",
}

# Keep the Tools catalog focused on supported Studio tiles. Internal workflow
# adapters may remain available to the workflow API without a standalone tile.
TOOL_CATALOG_WORKFLOW_IDS = {
    "library.prepare",
    "library.centroids",
    "detection.scannet",
    "detection.pesto",
    "detection.surf2spot",
    "detection.masif_seed",
    "detection.surface_analysis",
    "design.rfdiffusion_generate",
    "design.rfdiffusion_beta_pairing",
    "design.kejia_idr",
    "design.peptide_trap",
    "design.bindcraft2",
    "design.bindcraft2_generator",
    "design.rfdiffusion3",
    "design.rfd4_proteina",
    "refolding.rfd4_proteina",
    "design.rfantibody",
    "design.rfantibody_mimosa",
    "design.rfpeptides",
    "design.boltzgen",
    "design.boltzgen_mimosa",
    "design.boltzgen_generator",
    "design.pxdesign_generate",
    "design.genie3_generate",
    "design.esmfold2_binder",
    "design.protpardelle_generate",
    "design.proteina_complexa",
    "design.proteina_complexa_generate",
    "sequence.mpnn",
    "sequence.abmpnn",
    "analysis.sabr",
    "sequence.esm_if1",
    "sequence.caliby",
    "sequence.boltzgen_inverse",
    "refolding.af3_jax_cofolding",
    "refolding.optimized_prediction",
    "refolding.af2_initial_guess",
    "refolding.openfold3",
    "refolding.rf3",
    "refolding.protenix",
    "refolding.protenix_cli",
    "refolding.boltzgen_fold",
    "analysis.ipsae",
    "analysis.binderscope",
    "analysis.rosetta",
    "analysis.validation_acceptance",
    "ligand.complex.prepare",
    "ligand.docking.vina",
    "ligand.docking.autodock_gpu",
    "ligand.docking.gnina",
    "ligand.docking.unidock_pro",
    "ligand.docking.rosetta_ligand",
    "ligand.affinity.nesso",
    "ligand.generation.lddm",
    "ligand.generation.pocketxmol",
    "ligand.generation.omtra",
    "ligand.generation.flowr_root",
    "ligand.generation.conditar",
    "ligand.generation.paopt",
    "ligand.generation.drugrpg",
    "ligand.generation.pfm",
    "ligand.generation.pocketflow",
    "ligand.generation.pgmg",
    "ligand.md.gromacs",
    "ligand.md.openmm",
    "ligand.md.analysis",
    "simulation.bioemu2",
}

# Upstream code versions for tool tiles. Values describe engine code or a pinned
# source revision, never the Docker image tag. Some legacy shared images have no
# pinned upstream release, which is shown explicitly.
TOOL_WORKFLOW_CODE_VERSIONS = {
    "library.prepare": ("Studio · library fingerprint preparation v1", False),
    "library.centroids": ("Studio · library centroid preparation v1", False),
    "simulation.bioemu2": ("BioEmu-2 · 0.1.0rc0 · Microsoft source f99e0d8", False),
    "ligand.complex.prepare": ("Complex Preparation · Studio bundled runtime", False),
    "ligand.docking.vina": ("AutoDock Vina · 1.2.7", False),
    "ligand.docking.autodock_gpu": ("AutoDock-GPU · 1.6", False),
    "ligand.docking.gnina": ("GNINA · installed docking-suite version", False),
    "ligand.docking.unidock_pro": ("Uni-Dock Pro · installed docking-suite version", False),
    "ligand.docking.rosetta_ligand": ("RosettaLigand · Studio bundled OpenVS runner", False),
    "ligand.affinity.nesso": ("Nesso-1 · 1.0.0", False),
    "ligand.generation.lddm": ("LDDM · f254fb4", False),
    "ligand.generation.pocketxmol": ("PocketXMol · 65488cf", False),
    "ligand.generation.omtra": ("OMTRA · Studio workflow adapter", False),
    "ligand.generation.flowr_root": ("FLOWR.root · Studio workflow adapter", False),
    "ligand.generation.conditar": ("conDitar · permissioned local image", False),
    "ligand.generation.paopt": ("paOPT · permissioned local image", False),
    "ligand.generation.drugrpg": ("DrugRPG · Studio workflow adapter", False),
    "ligand.generation.pfm": ("PFM · Studio workflow adapter", False),
    "ligand.generation.pocketflow": ("PocketFlow · Studio workflow adapter", False),
    "ligand.generation.pgmg": ("PGMG · Studio workflow adapter", False),
    "ligand.md.gromacs": ("GROMACS · Studio bundled runtime", False),
    "ligand.md.openmm": ("OpenMM · Studio bundled runtime", False),
    "ligand.md.analysis": ("MDTraj · Studio trajectory analysis", False),
    "sequence.esm_if1": ("ESM-IF1 · fair-esm 2.0.1 · 2b369911", False),
    "sequence.caliby": ("Caliby 0.1 · 41d31560", True),
    "sequence.boltzgen_inverse": ("BoltzGen 0.3.2 · a07fea1", False),
    "design.pxdesign_generate": ("PXDesign f788441", False),
    "design.rfdiffusion_generate": ("RFdiffusion · version not pinned", False),
    "design.rfdiffusion_beta_pairing": ("mn-rfdiffusion:latest · published SS/adjacency conditioning", False),
    "design.kejia_idr": ("Kejia 7f65c18 · native template threading", False),
    "design.peptide_trap": ("Sahtoe 2024 · published trap scaffolds", False),
    "design.genie3_generate": ("Genie3 d77ae5a", False),
    "design.protpardelle_generate": ("Protpardelle-1C ee37840", False),
    "detection.scannet": ("ScanNet a61623c", False),
    "detection.pesto": ("PeSTo ba651aa", False),
    "detection.surf2spot": ("Surf2Spot a38a7a7", False),
    "detection.masif_seed": ("MaSIF-seed 7715391", False),
    "detection.surface_analysis": ("PeSTo · ScanNet · Surf2Spot · MaSIF-seed", False),
    "design.bindcraft2": ("BindCraft 2 v1.0.4", False),
    "design.bindcraft2_generator": ("BindCraft 2 v1.0.4", False),
    "design.rfdiffusion_classic": ("RFdiffusion · version not pinned", False),
    "design.rfdiffusion3": ("RFdiffusion3 a36d29c · Foundry 0.1.9", False),
    "refolding.rfd4_proteina": ("RFD4-Proteina 0.0.3 · native folding", False),
    "design.rfd4_proteina": ("RFDiffusion4-Proteina 0.0.3 · 9d35ee4", False),
    "design.rfantibody": ("RFantibody 3e5a5bb · CUDA 11.8", False),
    "design.rfantibody_mimosa": ("MIMOSA-RFA · RFantibody 3e5a5bb", False),
    "design.rfpeptides": ("RFpeptides · RFdiffusion 86507b6 · CUDA 12.8", False),
    "design.boltzgen": ("BoltzGen 0.3.2", True),
    "design.boltzgen_mimosa": ("MIMOSA-BG · BoltzGen 0.3.2", True),
    "design.boltzgen_generator": ("BoltzGen 0.3.2", True),
    "design.pxdesign": ("PXDesign f788441 · Protenix 0.5.0+pxd", False),
    "design.genie3": ("Genie3 d77ae5a", False),
    "design.esmfold2_binder": ("Biohub ESM 3.4.1.post1", False),
    "design.protpardelle_1c": ("Protpardelle-1C ee37840", False),
    "design.proteina_complexa": ("Proteina-Complexa 1.0.0 · 32b71ae", False),
    "design.proteina_complexa_generate": ("Proteina-Complexa 1.0.0 · 32b71ae", False),
    "sequence.ligandmpnn": ("LigandMPNN · version not pinned", False),
    "sequence.foundry_mpnn": ("Foundry MPNN · version not pinned", False),
    "sequence.mpnn": ("LigandMPNN 26ec57a · Foundry 0.1.9", False),
    "refolding.af3_jax_cofolding": (
        "ColabFold AF3 preview f850f0f · AlphaFold3 3.1.14",
        True,
    ),
    "refolding.af2_initial_guess": ("AlphaFold2 · version not pinned", False),
    "refolding.openfold3": ("OpenFold3 0.5.0 · P2 weights", False),
    "refolding.rf3": ("RosettaFold3 · version not pinned", False),
    "refolding.protenix": ("Protenix 0.5.0+pxd", False),
    "refolding.protenix_cli": ("Protenix 2.0.0", False),
    "refolding.boltzgen_fold": ("BoltzGen 0.3.2", True),
    "analysis.ipsae": ("ipSAE · version not pinned", False),
    "analysis.binderscope": ("BinderScope c61ef9c · integrated source", False),
    "analysis.rosetta": ("PyRosetta · version not pinned", False),
}

# Related prediction modes share one Tools card. The existing workflows remain
# separate behind the card because their inputs and settings differ.
TOOL_CATALOG_GROUPS = {
    "boltzgen_fold": {
        "title": "BoltzGen Fold",
        "description": "Candidate-batch refolding or target PDB plus binder sequence.",
        "engine_version": "BoltzGen 0.3.2",
        "engine_optimized": True,
        "modes": {
            "Refold candidate batch": "refolding.boltzgen_fold",
            "Target PDB + binder sequence": "refolding.optimized_prediction",
        },
        "engine_profiles": {"Target PDB + binder sequence": "boltzgen_fold"},
    },
}

TOOL_CATALOG_WORKFLOW_IDS = (*TOOL_CATALOG_WORKFLOW_IDS, *QUALITY_WORKFLOWS)

# Each optimized prediction profile gets its own engine card. Version labels
# below refer to the upstream engine code, not the container or batch wrapper.
OPTIMIZED_IMAGE_TILE_IDS = tuple(optimized_prediction.ENGINE_PROFILES)
AF2_TILE_WORKFLOW = "refolding.af2_initial_guess_input"

def _workflow_discovery_tags(row: dict) -> dict[str, list[str]]:
    workflow = str(row.get("workflow", ""))
    if workflow in QUALITY_WORKFLOWS:
        profile = QUALITY_WORKFLOWS[workflow]
        function = "Diversity Analysis" if profile.get("cohort") else {
            "sequence": "Sequence Analysis", "interface": "Interface Analysis",
            "ligand": "Protein–Ligand Analysis", "structure": "Structure Quality",
        }[profile["kind"]]
        return {"functions": [function], "modalities": ["Protein", "Ligand"] if profile["kind"] == "ligand" else ["Protein", "Antibody", "Peptide"],
                "input_types": ["Existing Results", "Database Selection"]}
    prefix = workflow.split(".", 1)[0]
    function_by_prefix = {
        "library": "Library Preparation",
        "ligand": "Molecular Docking & Affinity",
        "simulation": "Biomolecular Simulation",
        "target": "Target Preparation",
        "detection": "Binding-site Detection",
        "design": "De Novo Binder Design",
        "sequence": "Protein Design",
        "refolding": "Structure Prediction",
        "analysis": "Structure Scoring",
        "import": "Molecule Import",
        "benchmark": "Benchmarking",
        "campaign": "Multi-step Validation",
    }
    function_by_workflow = {
        "simulation.bioemu2": "Conformational Ensemble Sampling",
        "ligand.complex.prepare": "Complex Preparation",
        "target.prepare": "Structure Preparation",
        "target.crop": "Structure Cropping",
        "target.mask": "Target Masking",
        "target.refolding_candidates": "Candidate Preparation",
        "target.refolding_evaluation": "Structure Prediction",
        "design.esmfold2_screening": "Binder Screening",
        "refolding.esmfold_monomer": "Cofolding & Structure Prediction",
        "refolding.esmfold2_cofolding": "Cofolding & Structure Prediction",
        "refolding.af2_initial_guess": "Cofolding & Structure Prediction",
        "refolding.af3_jax_cofolding": "Cofolding & Structure Prediction",
        "refolding.optimized_prediction": "Cofolding & Structure Prediction",
        "sequence.ligandmpnn": "Protein Sequence Design",
        "sequence.mpnn": "Protein Sequence Design",
        "sequence.abmpnn": "Antibody Sequence Design",
        "analysis.sabr": "Antibody Numbering",
        "sequence.esm_if1": "Inverse Folding",
        "sequence.caliby": "Inverse Folding",
        "sequence.boltzgen_inverse": "Inverse Folding",
        "sequence.foundry_mpnn": "Protein Sequence Design",
        "sequence.pipeline": "Protein Sequence Design",
        "analysis.ipsae": "ipSAE & Structure Analysis",
        "analysis.binderscope": "BinderScope",
        "analysis.rosetta": "Rosetta Structure Scoring",
        "ligand.docking.vina": "Molecular Docking",
        "ligand.docking.autodock_gpu": "Molecular Docking",
        "ligand.docking.gnina": "Molecular Docking",
        "ligand.docking.unidock_pro": "Molecular Docking",
        "ligand.docking.rosetta_ligand": "Molecular Docking",
        "ligand.affinity.nesso": "Affinity Prediction",
        "ligand.generation.lddm": "Ligand Generation",
        "ligand.generation.pocketxmol": "Ligand Generation",
        "ligand.generation.omtra": "Ligand Generation",
        "ligand.generation.flowr_root": "Ligand Generation",
        "ligand.generation.conditar": "Ligand Generation",
        "ligand.generation.paopt": "Ligand Generation",
        "ligand.generation.drugrpg": "Ligand Generation",
        "ligand.generation.pfm": "Ligand Generation",
        "ligand.generation.pocketflow": "Ligand Generation",
        "ligand.generation.pgmg": "Ligand Generation",
        "ligand.md.gromacs": "Molecular Dynamics",
        "ligand.md.openmm": "Molecular Dynamics",
        "ligand.md.analysis": "MD Analysis",
    }
    modalities_by_prefix = {
        "library": ["Ligand"],
        "simulation": ["Protein", "Ligand", "DNA", "RNA"],
        "ligand": ["Protein", "Ligand"],
        "target": ["Protein"],
        "detection": ["Protein"],
        "design": ["Protein"],
        "sequence": ["Protein", "Antibody", "Peptide"],
        "refolding": ["Protein", "Antibody", "Peptide"],
        "analysis": ["Protein", "Antibody", "Peptide"],
        "import": ["Protein", "Antibody", "Peptide"],
        "benchmark": ["Protein", "Antibody", "Peptide"],
        "campaign": ["Protein", "Antibody", "Peptide"],
    }
    input_types_by_prefix = {
        "library": ["ID / SMILES table", "Prepared fingerprint manifest"],
        "simulation": ["Protein Complex (PDB/mmCIF)", "Protein Sequences"],
        "ligand": ["Protein Structure (PDB)", "Compound Set"],
        "target": ["Structure (PDB)"],
        "detection": ["Structure (PDB)"],
        "design": ["Structure (PDB)"],
        "sequence": ["Existing Results"],
        "refolding": ["Existing Results"],
        "analysis": ["Existing Results"],
        "import": ["Structure (PDB)", "Sequence"],
        "benchmark": ["Existing Results"],
        "campaign": ["Existing Results"],
    }
    input_types_by_workflow = {
        "design.kejia_idr": ["FASTA / Sequence", "Structure (PDB)", "Structure (mmCIF)"],
        "design.peptide_trap": ["FASTA / Sequence", "Template Library"],
        "ligand.complex.prepare": ["Protein–Ligand Complex (PDB)"],
        "sequence.mpnn": ["Structure (PDB)", "Structure (mmCIF)", "Existing Results"],
        "sequence.abmpnn": ["Antibody Structure (PDB/mmCIF)", "Existing Antibody Designs"],
        "analysis.sabr": ["Antibody Structure (PDB/mmCIF)", "Existing Antibody Designs"],
        "sequence.esm_if1": ["Structure (PDB)", "Structure (mmCIF)", "Existing Results"],
        "sequence.caliby": ["Structure (PDB)", "Structure (mmCIF)", "Existing Results"],
        "sequence.boltzgen_inverse": ["Structure (PDB)", "Structure (mmCIF)", "Existing Results"],
        **{workflow:["Structure (PDB)","Structure (mmCIF)"] for workflow in generator_input.WORKFLOWS.values()},
        "design.boltzgen": ["Structure (PDB)", "Structure (mmCIF)"],
        "design.boltzgen_generator": ["Structure (PDB)", "Structure (mmCIF)"],
        "design.proteina_complexa": ["Structure (PDB)", "Structure (mmCIF)"],
        "design.proteina_complexa_generate": ["Structure (PDB)", "Structure (mmCIF)"],
        "design.rfdiffusion3": ["Structure (PDB)", "Structure (mmCIF)", "JSON"],
        "refolding.rfd4_proteina": ["FASTA / Sequence", "Structure (PDB)", "Structure (mmCIF)"],
        "design.rfd4_proteina": ["Structure (PDB)", "Structure (mmCIF)", "Annotated CIF", "YAML / JSON"],
        "design.rfpeptides": ["Structure (PDB)", "Structure (mmCIF)"],
        "refolding.esmfold2_cofolding": ["Protein Sequences"],
        "refolding.af3_jax_cofolding": ["Molecule Fields"],
        "refolding.optimized_prediction": ["Molecule Fields"],
        "analysis.ipsae": ["Existing Results"],
        "analysis.rosetta": ["Existing Results"],
        "ligand.docking.vina": ["Protein Structure (PDB)", "Database Compound Set"],
        "ligand.docking.autodock_gpu": ["Protein Structure (PDB)", "Database Compound Set"],
        "ligand.docking.gnina": ["Protein Structure (PDB)", "Database Compound Set"],
        "ligand.docking.unidock_pro": ["Protein Structure (PDB)", "Database Compound Set"],
        "ligand.docking.rosetta_ligand": ["Protein Structure (PDB)", "Database Compound Set", "Optional Reference Ligand (SDF)"],
        "ligand.affinity.nesso": ["Protein Structure (PDB)", "Database Compound Set"],
        "ligand.generation.lddm": ["Protein Structure (PDB)", "Reference Ligand (SDF)"],
        "ligand.generation.pocketxmol": ["Protein Structure (PDB)", "Pocket Structure (PDB)"],
        "ligand.generation.omtra": ["Protein Structure (PDB)", "Optional Pocket (PDB)", "Optional Reference Ligand (SDF)"],
        "ligand.generation.flowr_root": ["Protein Structure (PDB)", "Pocket Structure (PDB)", "Reference Ligand (SDF)"],
        "ligand.generation.conditar": ["Protein Structure (PDB)", "Pocket Structure (PDB)"],
        "ligand.generation.paopt": ["Protein Structure (PDB)", "Pocket Structure (PDB)", "Reference Ligand (SDF)", "ADMET Objectives"],
        "ligand.generation.drugrpg": ["Protein Structure (PDB)", "Pocket Structure (PDB)"],
        "ligand.generation.pfm": ["Protein Structure (PDB)", "Pocket Structure (PDB)"],
        "ligand.generation.pocketflow": ["Protein Structure (PDB)", "Pocket Structure (PDB)"],
        "ligand.generation.pgmg": ["Protein Structure (PDB)", "Pharmacophore Hypothesis (.posp)"],
        "ligand.md.gromacs": ["Protein–Ligand Complex (PDB)", "Bound Ligand (SDF)"],
        "ligand.md.openmm": ["Protein–Ligand Complex (PDB)", "Bound Ligand (SDF)"],
        "ligand.md.analysis": ["MD Trajectory", "Matching Topology"],
        "simulation.bioemu2": ["Target–Design Complexes (PDB/mmCIF)", "Protein Sequences"],
        "import.bindcraft": ["Result Folder"],
        "import.table": ["CSV / Table", "Structure (PDB)", "Sequence"],
        "benchmark.dataset": ["CSV / Table", "Result Folder", "Structure (PDB)", "Existing Results"],
        "benchmark.precomputed_metrics": ["CSV / Table"],
        "benchmark.esmfold2": ["CSV / Table", "Sequence"],
        "benchmark.refold_candidates": ["Existing Results"],
        "benchmark.refolding_capacity": ["Existing Results", "Structure (PDB)"],
        "benchmark.design_capacity": ["Structure (PDB)"],
        "benchmark.practical_capacity": ["Existing Results"],
        "benchmark.collection": ["Existing Results"],
        "benchmark.matrix_workspace": ["Existing Results"],
        "campaign.validation_sequence": ["Existing Results"],
    }
    functions = [function_by_workflow.get(workflow, function_by_prefix.get(prefix, "Molecular Design"))]
    if workflow in WORKFLOW_ROLES:
        functions.insert(0, WORKFLOW_ROLES[workflow])
    if workflow in {"design.bindcraft2", "design.bindcraft2_generator"}:
        functions.append("Antibody Design")
    return {
        "functions": functions,
        "modalities": (
            ["Protein", "Antibody", "Peptide"]
            if workflow in {"design.bindcraft2", "design.bindcraft2_generator"}
            else ["Protein", "DNA", "RNA", "Ligand"]
            if workflow in {"refolding.af3_jax_cofolding", "design.rfdiffusion3"}
            else ["Protein", "Peptide"] if workflow in {"design.rfpeptides", "design.kejia_idr", "design.peptide_trap"}
            else ["Protein", "Ligand"] if workflow in {"design.proteina_complexa", "design.proteina_complexa_generate"}
            else ["Protein"]
            if workflow == "refolding.esmfold2_cofolding"
            else modalities_by_prefix.get(prefix, ["Protein"])
        ),
        "input_types": input_types_by_workflow.get(
            workflow,
            input_types_by_prefix.get(prefix, ["Structure (PDB)"]),
        ),
    }


def catalog_rows(workflows: list[dict]) -> list[dict]:
    all_workflows = [
        row for row in workflows
        if row["workflow"] in TOOL_CATALOG_WORKFLOW_IDS
    ]
    workflow_by_id = {row["workflow"]: row for row in all_workflows}
    tags_by_workflow = {row["workflow"]: _workflow_discovery_tags(row) for row in all_workflows}
    grouped_workflow_ids: set[str] = set()
    tool_rows: list[dict] = []
    for group_id, group in TOOL_CATALOG_GROUPS.items():
        modes = {
            label: workflow_id
            for label, workflow_id in group["modes"].items()
            if workflow_id in workflow_by_id
        }
        if len(modes) < 2:
            continue
        member_rows = [workflow_by_id[workflow_id] for workflow_id in modes.values()]
        group_tags = {
            key: list(dict.fromkeys(tag for row in member_rows for tag in tags_by_workflow[row["workflow"]][key]))
            for key in ("functions", "modalities", "input_types")
        }
        tool_rows.append(
            {
                "key": f"group-{group_id}",
                "group_id": group_id,
                "title": group["title"],
                "description": group["description"],
                "workflow_ids": list(modes.values()),
                "workflow_rows": member_rows,
                "modes": modes,
                "tags": group_tags,
                "mode_count": len(modes),
                "engine_version": group.get("engine_version"),
                "engine_optimized": group.get("engine_optimized", False),
                "references": references_for_tool(
                    list(modes.values()),
                    optimized_engine=next(iter(group.get("engine_profiles", {}).values()), None),
                ),
            }
        )
        grouped_workflow_ids.update(modes.values())

    for row in all_workflows:
        workflow_id = row["workflow"]
        if workflow_id == "refolding.optimized_prediction":
            continue
        if workflow_id == "analysis.validation_acceptance":
            continue
        if workflow_id in grouped_workflow_ids:
            continue
        base_title = workflow_id.split(".", 1)[-1].replace("_", " ").title()
        tool_rows.append(
            {
                "key": workflow_id,
                "group_id": None,
                "workflow": workflow_id,
                "title": {"library.prepare": "Prepare library", "library.centroids": "Create centroids"}.get(workflow_id) or ("RFD4-Proteina Structure Prediction" if workflow_id == "refolding.rfd4_proteina" else (QUALITY_WORKFLOWS.get(workflow_id, {}).get("name") or WORKFLOW_TITLES.get(base_title, base_title))),
                "description": row["description"],
                "workflow_ids": [workflow_id],
                "workflow_rows": [row],
                "tags": tags_by_workflow[workflow_id],
                "engine_version": (
                    TOOL_WORKFLOW_CODE_VERSIONS.get(workflow_id, (None, False))[0]
                ),
                "engine_optimized": (
                    TOOL_WORKFLOW_CODE_VERSIONS.get(workflow_id, (None, False))[1]
                ),
                "references": references_for_tool([workflow_id]),
            }
        )

    optimized_row = workflow_by_id.get("refolding.optimized_prediction")
    if optimized_row:
        for engine_id in OPTIMIZED_IMAGE_TILE_IDS:
            if engine_id == "boltzgen_fold":
                continue
            profile = optimized_prediction.ENGINE_PROFILES[engine_id]
            supported = set(profile.get("molecule_inputs", ()))
            tool_rows.append({
                "key": f"optimized-image-{engine_id}",
                "group_id": None,
                "workflow": "refolding.optimized_prediction",
                "optimized_engine": engine_id,
                "title": profile["label"],
                "description": profile.get("input_summary", "Predict protein structures."),
                "engine_version": profile.get("code_version"),
                "engine_optimized": profile.get("optimized", False),
                "workflow_ids": ["refolding.optimized_prediction"],
                "workflow_rows": [optimized_row],
                "references": references_for_tool(
                    ["refolding.optimized_prediction"], optimized_engine=engine_id
                ),
                "tags": {
                    "functions": ["Cofolding & Structure Prediction"],
                    "modalities": ["Protein", *(["DNA"] if "dna_sequences" in supported else []),
                                   *(["RNA"] if "rna_sequences" in supported else []),
                                   *(["Ligand"] if supported.intersection({"ligand_ccd_codes", "ligand_smiles"}) else [])],
                    "input_types": ["Molecule Fields" if supported else "Protein Sequences"],
                },
            })

        # Acceptance is an evaluation-only tool over saved candidate evidence.
        validation_row = workflow_by_id.get("analysis.validation_acceptance")
        tool_rows.append({
            "key": "studio.validation-acceptance",
            "group_id": None,
            "workflow": "studio.validation_acceptance",
            "title": "Validation & Acceptance",
            "description": "Apply configurable, format-aware acceptance criteria to saved design measurements and structures. No prediction engine is launched.",
            "workflow_ids": ["studio.validation_acceptance"],
            "workflow_rows": [validation_row] if validation_row else [],
            "engine_version": None,
            "engine_optimized": False,
            "references": [],
            "tags": {
                "functions": ["Validation & Acceptance"],
                "modalities": ["Protein", "Antibody", "Peptide"],
                "input_types": ["Existing Results"],
            },
        })

    tool_rows.sort(key=lambda row: row["title"].casefold())
    return tool_rows
