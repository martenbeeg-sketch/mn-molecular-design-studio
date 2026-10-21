"""Public QC tool contracts shared by tiles, scheduling, and pipeline nodes."""
TOOLS = {
    "sequence_qc": {"name": "Sequence Composition", "image": "mn-protein-qc:1.0.0", "kind": "sequence", "metrics": ["qc_sequence_length", "qc_sequence_pi", "qc_sequence_gravy"]},
    "proteinsol": {"name": "Protein-Sol", "image": "mn-proteinsol:1.0.0", "kind": "sequence", "metrics": ["proteinsol_scaled_solubility"]},
    "biotite_interface": {"name": "Biotite Interface", "image": "mn-protein-qc:1.0.0", "kind": "interface", "metrics": ["qc_interface_contact_pairs"]},
    "sequence_clustering": {"name": "Sequence Clustering", "image": "mn-protein-qc:1.0.0", "kind": "sequence", "cohort": True, "metrics": ["diversity_cluster"]},
    "interface_clustering": {"name": "Interface Clustering", "image": "mn-protein-qc:1.0.0", "kind": "interface", "cohort": True, "metrics": ["diversity_cluster"]},
    "foldseek": {"name": "Foldseek", "image": "mn-foldseek:10-941cd33-qc1", "kind": "structure", "cohort": True, "metrics": ["diversity_cluster"]},
    "pep_patch": {"name": "PEP-Patch", "image": "mn-pep-patch:65da35b-qc1", "kind": "structure", "metrics": ["peppatch_total_area"]},
    "esm1v_score": {"name": "ESM-1v Scoring", "image": "mn-esm-if1-opt:2.0.1-2b369911", "kind": "sequence", "gpu": True,
        "references": ["esm_models/esm1v_t33_650M_UR90S_" + str(i) + ".pt" for i in range(1, 6)], "metrics": ["esm1v_native_probability"]},
    "esmif_score": {"name": "ESM-IF1 Scoring", "image": "mn-esm-if1-opt:2.0.1-2b369911", "kind": "structure", "gpu": True,
        "references": ["esm_models/esm_if1_gvp4_t16_142M_UR50.pt"], "metrics": ["esmif_native_log_likelihood"]},
    "plip": {"name": "PLIP", "image": "mn-plip:latest", "kind": "ligand", "metrics": ["plip_interaction_count"]},
    "pandamap": {"name": "PandaMap", "image": "mn-pandamap:0007347", "kind": "ligand", "metrics": ["pandamap_interaction_count"]},
    "posebusters": {"name": "PoseBusters", "image": "mn-posebusters:1a5f26a", "kind": "ligand", "metrics": ["posebusters_passed_all"]},
}

WORKFLOWS = {"analysis." + key: value for key, value in TOOLS.items()}
