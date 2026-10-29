"""Representative research references for the Studio tool catalog."""
from __future__ import annotations


def _doi(short: str, title: str, doi: str) -> dict[str, str]:
    return {"short": short, "title": title, "doi": doi, "url": f"https://doi.org/{doi}"}


def _url(short: str, title: str, url: str) -> dict[str, str]:
    return {"short": short, "title": title, "url": url}


PAPERS = {
    "morgan": _doi("Morgan fingerprints", "Extended-Connectivity Fingerprints", "10.1021/ci100050t"),
    "butina": _doi("Butina clustering", "Unsupervised data base clustering based on daylight's fingerprint and Tanimoto similarity: a fast and automated way to cluster small and large data sets", "10.1021/ci9803381"),
    "scannet": _doi("ScanNet", "ScanNet: an interpretable geometric deep learning model for protein-protein interface prediction", "10.1038/s41592-022-01490-7"),
    "pesto": _doi("PeSTo", "PeSTo: modeling protein-protein interfaces with geometric deep learning", "10.1038/s41467-023-37701-8"),
    "surf2spot": _doi("Surf2Spot", "Surf2Spot: protein surface hotspot prediction", "10.1093/bib/bbag542"),
    "masif": _doi("MaSIF", "MaSIF: a deep learning approach to protein surface interaction prediction", "10.1038/s41592-019-0666-6"),
    "rfdiffusion": _doi("RFdiffusion", "De novo design of protein structure and function with RFdiffusion", "10.1038/s41586-023-06415-8"),
    "kejia": _doi("Kejia IDR", "De novo design of intrinsically disordered region binders", "10.1126/science.adr8063"),
    "peptide_trap": _doi("Peptide Trap", "De novo design of peptide traps for targeting intrinsically disordered proteins", "10.1038/s41589-024-01578-5"),
    "bindcraft": _doi("BindCraft", "De novo design of high-affinity protein binders with BindCraft", "10.1038/s41586-025-09429-6"),
    "rfd3": _doi("RFdiffusion3", "De novo Design of All-atom Biomolecular Interactions with RFdiffusion3", "10.1101/2025.09.18.676967"),
    "rfd4": _url("RFD4-Proteina", "Biomolecular design through latent reasoning with RFD4-Proteina", "https://rfd4-proteina.org/assets/paper/RFD4-Proteina.pdf"),
    "rfantibody": _doi("RFantibody", "De novo design of high-affinity antibodies", "10.1038/s41586-025-09721-5"),
    "mimosa": _doi("MIMOSA", "MIMOSA: model-informed sequence optimization for antibody design", "10.64898/2026.09.24.753155"),
    "rfpeptides": _doi("RFpeptides", "De novo design of peptide binders with RFpeptides", "10.1038/s41589-025-01929-w"),
    "boltzgen": _doi("BoltzGen", "BoltzGen: generative modeling for biomolecular design", "10.1101/2025.11.20.689494"),
    "pxdesign": _doi("PXDesign", "PXDesign: Fast, Modular, and Accurate De Novo Design of Protein Binders", "10.1101/2025.08.15.670450"),
    "genie3": _doi("Genie3", "Genie 3: an all-atom generative model for protein design", "10.64898/2026.05.01.722168"),
    "esmfold2": _doi("ESMFold2", "Language modeling materializes a world model of protein biology", "10.64898/2026.06.03.729735"),
    "protpardelle": _doi("Protpardelle-1C", "Protpardelle-1C: a generative model for protein design", "10.1101/2025.08.18.670959"),
    "complexa": _url("Proteina-Complexa", "Scaling Atomistic Protein Binder Design with Generative Pretraining and Test-Time Compute", "https://arxiv.org/abs/2603.27950"),
    "proteinmpnn": _doi("ProteinMPNN", "Robust deep learning-based protein sequence design using ProteinMPNN", "10.1126/science.add2187"),
    "ligandmpnn": _doi("LigandMPNN", "Atomic context-conditioned protein sequence design using LigandMPNN", "10.1038/s41592-025-02626-1"),
    "abmpnn": _url("AbMPNN", "Antibody sequence design with AbMPNN", "https://arxiv.org/abs/2310.19513"),
    "anarci": _doi("ANARCI", "ANARCI: antigen receptor numbering and receptor classification", "10.1093/bioinformatics/btv552"),
    "esmif": _doi("ESM-IF1", "Robust deep learning-based protein sequence design using ESM-IF1", "10.1101/2022.04.10.487779"),
    "caliby": _doi("Caliby", "Caliby: sequence design and scoring with Potts models", "10.1101/2025.09.30.679633"),
    "alphafold3": _doi("AlphaFold 3", "Accurate structure prediction of biomolecular interactions with AlphaFold 3", "10.1038/s41586-024-07487-w"),
    "alphafold2": _doi("AlphaFold 2", "Highly accurate protein structure prediction with AlphaFold", "10.1038/s41586-021-03819-2"),
    "rf3": _doi("RosettaFold3", "RoseTTAFold All-Atom: enabling chemistry-aware modeling of biomolecular interactions", "10.1101/2025.08.14.670328"),
    "protenix1": _doi("Protenix V1", "Protenix: advancing structure prediction through a comprehensive AlphaFold 3 reproduction", "10.64898/2026.02.05.703733"),
    "protenix2": _doi("Protenix V2", "Protenix V2: scaling open all-atom structure prediction", "10.64898/2026.04.10.717613"),
    "boltz2": _doi("Boltz-2", "Boltz-2: towards accurate and efficient binding affinity prediction", "10.1101/2025.06.14.659707"),
    "esmfold": _doi("ESMFold", "Evolutionary-scale prediction of atomic-level protein structure with a language model", "10.1126/science.ade2574"),
    "chai1": _doi("Chai-1", "Chai-1: decoding the molecular interactions of life", "10.1101/2024.10.10.615955"),
    "atlasfold": _doi("AtlasFold-M", "AtlasFold: Protein structure prediction with metagenomic-scale language models", "10.64898/2026.09.04.749352"),
    "opendde": _doi("OpenDDE", "Folding, Reasoning, and Scaling with Open-source Drug Discovery Engine", "10.48550/arXiv.2607.03787"),
    "colabfold": _doi("ColabFold", "ColabFold: making protein folding accessible to all", "10.1038/s41592-022-01488-1"),
    "bioemu": _doi("BioEmu", "Scalable emulation of protein equilibrium ensembles with generative deep learning", "10.1126/science.adv9817"),
    "vina": _doi("AutoDock Vina", "AutoDock Vina: improving the speed and accuracy of docking", "10.1002/jcc.21334"),
    "autodock_gpu": _doi("AutoDock-GPU", "AutoDock-GPU: a graphics processing unit-accelerated docking tool", "10.1021/acs.jctc.0c01006"),
    "gnina": _doi("GNINA", "GNINA 1.0: molecular docking with deep learning", "10.1186/s13321-021-00522-2"),
    "unidock": _doi("Uni-Dock", "Uni-Dock: GPU-Accelerated Docking Enables Ultralarge Virtual Screening", "10.1021/acs.jctc.2c01145"),
    "rosetta_ligand": _doi("RosettaLigand", "RosettaLigand docking with full ligand and receptor flexibility", "10.1016/j.jmb.2008.11.010"),
    "nesso": _doi("Nesso-1", "Nesso-1: protein-ligand binding affinity prediction", "10.64898/2026.08.01.742196"),
    "lddm": _doi("LDDM", "LDDM: ligand design with diffusion models", "10.64898/2026.09.15.751537"),
    "pocketxmol": _doi("PocketXMol", "PocketXMol: efficient molecular generation for structure-based drug design", "10.1016/j.cell.2026.01.003"),
    "omtra": _url("OMTRA", "OMTRA: optimizing molecular generation with target-aware rewards", "https://arxiv.org/abs/2512.05080"),
    "flowr": _doi("FLOWR.root", "FLOWR: a flow-matching model for molecular generation", "10.1038/s41467-026-74130-9"),
    "conditar": _url("conDitar", "Generating Developable 3D Molecules via Pocket-Conditioned Diffusion and Property-Aware Optimization", "https://arxiv.org/abs/2607.12349"),
    "paopt": _url("paOPT", "Generating Developable 3D Molecules via Pocket-Conditioned Diffusion and Property-Aware Optimization", "https://arxiv.org/abs/2607.12349"),
    "drugrpg": _doi("DrugRPG", "Integrating chemical priors and physical laws to mitigate hallucinations in structure-based drug design", "10.1038/s42004-026-02115-2"),
    "pfm": _doi("PFM", "PFM: perturbed flow matching for structure-based drug design", "10.64898/2026.07.11.737913"),
    "pocketflow": _doi("PocketFlow", "PocketFlow: sampling 3D molecules for target protein pockets", "10.1038/s42256-024-00808-8"),
    "pgmg": _doi("PGMG", "Pharmacophore-guided molecular generation with a deep generative model", "10.1038/s41467-023-41454-9"),
    "gromacs": _doi("GROMACS", "GROMACS: high performance molecular simulations through multi-level parallelism from laptops to supercomputers", "10.1016/j.softx.2015.06.001"),
    "openmm": _doi("OpenMM", "OpenMM 7: rapid development of high performance algorithms for molecular dynamics", "10.1371/journal.pcbi.1005659"),
    "mdtraj": _doi("MDTraj", "MDTraj: a modern open library for the analysis of molecular dynamics trajectories", "10.1016/j.bpj.2015.08.015"),
    "biopython": _doi("Biopython", "Biopython: freely available Python tools for computational molecular biology and bioinformatics", "10.1093/bioinformatics/btp163"),
    "proteinsol": _doi("Protein-Sol", "Protein–Sol: a web tool for predicting protein solubility from sequence", "10.1093/bioinformatics/btx345"),
    "biotite": _doi("Biotite", "Biotite: a unifying open source computational biology framework in Python", "10.1186/s12859-018-2367-z"),
    "mmseqs": _doi("MMseqs2", "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets", "10.1038/nbt.3988"),
    "foldseek": _doi("Foldseek", "Fast and accurate protein structure search with Foldseek", "10.1038/s41587-023-01773-0"),
    "peppatch": _doi("PEP-Patch", "PEP-Patch: Electrostatics in Protein–Protein Recognition, Specificity, and Antibody Developability", "10.1021/acs.jcim.3c01490"),
    "esm1v": _doi("ESM-1v", "Language models enable zero-shot prediction of the effects of mutations on protein function", "10.1101/2021.07.09.450648"),
    "plip": _doi("PLIP", "PLIP: fully automated protein-ligand interaction profiler", "10.1093/nar/gkab294"),
    "pandamap": _doi("PandaMap", "PandaMap: protein-ligand interaction analysis with geometric deep learning", "10.64898/2026.08.06.743421"),
    "posebusters": _doi("PoseBusters", "PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel molecules", "10.1039/D3SC04185A"),
    "ipsae": _doi("ipSAE", "ipSAE: a score for confident satisfactory protein-protein interactions predicted by AlphaFold2 and AlphaFold3", "10.1101/2025.02.10.637595"),
    "rosetta": _doi("Rosetta", "Macromolecular modeling and design in Rosetta: recent methods and frameworks", "10.1038/s41592-020-0848-2"),
}


# Keys are workflow IDs or optimized-prediction engine IDs. A tile may cite
# several closely related papers when it bundles more than one method.
WORKFLOW_REFERENCES = {
    "library.prepare": ["morgan"],
    "library.centroids": ["butina"],
    "simulation.bioemu2": ["bioemu"],
    "detection.scannet": ["scannet"],
    "detection.pesto": ["pesto"],
    "detection.surf2spot": ["surf2spot"],
    "detection.masif_seed": ["masif"],
    "detection.surface_analysis": ["scannet", "pesto", "surf2spot", "masif"],
    "design.rfdiffusion_generate": ["rfdiffusion"],
    "design.rfdiffusion_beta_pairing": ["rfdiffusion"],
    "design.kejia_idr": ["kejia"],
    "design.peptide_trap": ["peptide_trap"],
    "design.bindcraft2": ["bindcraft"],
    "design.bindcraft2_generator": ["bindcraft"],
    "design.rfdiffusion3": ["rfd3"],
    "design.rfd4_proteina": ["rfd4"],
    "refolding.rfd4_proteina": ["rfd4"],
    "design.rfantibody": ["rfantibody"],
    "design.rfantibody_mimosa": ["rfantibody", "mimosa"],
    "design.rfpeptides": ["rfpeptides"],
    "design.boltzgen": ["boltzgen"],
    "design.boltzgen_mimosa": ["boltzgen", "mimosa"],
    "design.boltzgen_generator": ["boltzgen"],
    "design.pxdesign_generate": ["pxdesign"],
    "design.genie3_generate": ["genie3"],
    "design.esmfold2_binder": ["esmfold2"],
    "design.protpardelle_generate": ["protpardelle"],
    "design.proteina_complexa": ["complexa"],
    "design.proteina_complexa_generate": ["complexa"],
    "sequence.mpnn": ["proteinmpnn", "ligandmpnn"],
    "sequence.abmpnn": ["abmpnn"],
    "analysis.sabr": ["anarci"],
    "sequence.esm_if1": ["esmif"],
    "sequence.caliby": ["caliby"],
    "sequence.boltzgen_inverse": ["boltzgen"],
    "refolding.af3_jax_cofolding": ["alphafold3"],
    "refolding.af2_initial_guess": ["alphafold2"],
    "refolding.openfold3": ["alphafold3"],
    "refolding.rf3": ["rf3"],
    "refolding.protenix": ["protenix1"],
    "refolding.protenix_cli": ["protenix2"],
    "refolding.boltzgen_fold": ["boltzgen"],
    "analysis.ipsae": ["ipsae"],
    "analysis.rosetta": ["rosetta"],
    "ligand.complex.prepare": ["openmm"],
    "ligand.docking.vina": ["vina"],
    "ligand.docking.autodock_gpu": ["autodock_gpu"],
    "ligand.docking.gnina": ["gnina"],
    "ligand.docking.unidock_pro": ["unidock"],
    "ligand.docking.rosetta_ligand": ["rosetta_ligand"],
    "ligand.affinity.nesso": ["nesso"],
    "ligand.generation.lddm": ["lddm"],
    "ligand.generation.pocketxmol": ["pocketxmol"],
    "ligand.generation.omtra": ["omtra"],
    "ligand.generation.flowr_root": ["flowr"],
    "ligand.generation.conditar": ["conditar"],
    "ligand.generation.paopt": ["paopt"],
    "ligand.generation.drugrpg": ["drugrpg"],
    "ligand.generation.pfm": ["pfm"],
    "ligand.generation.pocketflow": ["pocketflow"],
    "ligand.generation.pgmg": ["pgmg"],
    "ligand.md.gromacs": ["gromacs"],
    "ligand.md.openmm": ["openmm"],
    "ligand.md.analysis": ["mdtraj"],
    "analysis.sequence_qc": ["biopython"],
    "analysis.proteinsol": ["proteinsol"],
    "analysis.biotite_interface": ["biotite"],
    "analysis.sequence_clustering": ["mmseqs"],
    "analysis.interface_clustering": ["foldseek"],
    "analysis.foldseek": ["foldseek"],
    "analysis.pep_patch": ["peppatch"],
    "analysis.esm1v_score": ["esm1v"],
    "analysis.esmif_score": ["esmif"],
    "analysis.plip": ["plip"],
    "analysis.pandamap": ["pandamap"],
    "analysis.posebusters": ["posebusters"],
}


ENGINE_REFERENCES = {
    "alphafold3": ["alphafold3"],
    "af3_torch": ["alphafold3"],
    "openfold3_jax": ["alphafold3"],
    "boltz2": ["boltz2"],
    "protenix_v2": ["protenix2"],
    "rf3_foundry": ["rf3"],
    "esmfold": ["esmfold"],
    "esmfold2": ["esmfold2"],
    "esmfold2_fast": ["esmfold2"],
    "chai1": ["chai1"],
    "atlasfold": ["atlasfold"],
    "openfold3_ob0": ["alphafold3"],
    "opendde_native": ["opendde"],
    "protenix_v1": ["protenix1"],
    "colabfold": ["colabfold"],
    "boltzgen_fold": ["boltzgen"],
}


def references_for_tool(
    workflow_ids: list[str] | tuple[str, ...] = (),
    *,
    optimized_engine: str | None = None,
) -> list[dict[str, str]]:
    """Return stable, de-duplicated paper links for a catalog tile."""
    keys: list[str] = []
    for workflow_id in workflow_ids:
        keys.extend(WORKFLOW_REFERENCES.get(workflow_id, ()))
    if optimized_engine:
        keys.extend(ENGINE_REFERENCES.get(optimized_engine, ()))
    references: list[dict[str, str]] = []
    seen: set[str] = set()
    for key in keys:
        paper = PAPERS[key]
        if paper["url"] not in seen:
            seen.add(paper["url"])
            references.append(paper)
    return references
