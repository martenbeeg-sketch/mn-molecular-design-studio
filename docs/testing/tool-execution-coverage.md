# Tool execution and Results coverage

Each row covers one visible tile mode. The live submissions were made in **Tool Validation** through MCP against release `8edf4e5`.

The project-wide evidence includes **90 completed jobs**, **29 failed jobs**, and **0 awaiting final Results verification**. Some jobs exercise the same tile or a workflow available only through the API.

**81 of 95 tiles** have a completed production job for every visible mode. **95 of 95** have completed evidence for every visible mode when patched native runs are included. Setup and schema checks are covered in the connection map.

Live completion requires a successful scientific result retrieved through MCP. Earlier saved results were opened through MCP without rerunning those jobs. A dash means the check has no evidence in this audit. Isolated replays use the patched checkout and the installed native runtime.

The live Interface Clustering run rejected the selected models because they lacked a verified shared target structure. An isolated native fixture with identical target coordinates and numbering clustered both inputs and returned two contact maps.

The patched end-to-end MCP retests completed for Sequence Composition, AutoDock-GPU and Surf2Spot. AutoDock-GPU had one opted-in SQL observation and a byte-identical saved pose. The same two-item Boltz-2 batch completed on the affected rootful Docker worker after the permissions fix. The combined Surface Analysis replay completed all four methods with thirteen artifacts on the other worker. The Results-only follow-up runs produced zero SQL result observations across the three workstation databases checked.

The next batch verified ESM-1v Scoring, validation from saved Boltz-2 models, three inverse-folding or sequence tools, BioEmu-2, several prediction engines, and two ligand generators. The native OpenMM trajectory was verified through patched result registration and a fresh MCP-to-MDTraj analysis. A patched GROMACS worker also completed preparation, equilibration, and production. MD Analysis read three GROMACS frames with the exact bonded preparation topology. The table shows which runs used the deployed release and which needed patched replays.

The latest batch exercised eight small-molecule generators and four prediction tiles. It verified the saved LDDM-to-Vina compound-set handoff using the patched MCP API. The native Mol-1 replay on PC004033 saved both its structure and a 20 × 20 PAE matrix. The final Results-only SQL checks are recorded by workstation in the matrix.

The next batch used the same saved Boltz-2 cohort for native OpenFold3, RosettaFold3 and both Protenix adapters. The isolated BoltzGen Fold replays repaired cache mounting and worker access to confidence archives. An authenticated BinderScope MCP test attached one evaluation to its exact SQL structure. A separate live ColabFold check acquired one alignment and the repeat reused it from disk.

The OpenFold3 P2 JAX image was updated on PC004032, PC004033 and the RTX 4090 coordinator. Both RTX 5090 workers completed single-sequence jobs submitted through production MCP and returned a verified model through Results. The earlier failures remain in the history.

The design batch submitted 15 jobs through production MCP. Thirteen completed and two required adapter fixes. The patched native replays recovered RFdiffusion fixed-length generation and Proteina-Complexa Generator result registration. The BoltzGen Generator replay also opened its saved archive as the worker on PC004033. Both downstream handoffs retained their parent candidate IDs.

The follow-up checked all ten remaining tiles and exercised two additional tile checks through live MCP. Eleven jobs were queued, including the MIMOSA RFantibody retry, the RFdiffusion range check, and Sequence Composition on one saved ProteinMPNN candidate. Seven completed and four failed. MIMOSA BoltzGen completed with a documented backbone integrity rejection and zero accepted designs. Isolated replays reached native inference on both BindCraft2 workers, registered the full Proteina-Complexa output, staged the MIMOSA RFantibody inputs across workstations, and ran MD Analysis through authenticated MCP. BinderScope attached one evaluation to the exact Tool Validation structure. No other test job produced a molecular-library observation on the three workstations. Nine headless Streamlit Results checks and 223 focused tests passed. [Official AlphaFold3 model-weight terms](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)

The final fifteen-tile pass completed official AlphaFold3 through production MCP with one 20-residue structure and a finite 20 × 20 PAE matrix. The other fourteen completed through authenticated MCP in isolated Studio workspaces using the patched checkout. Both BindCraft2 smoke tests completed and recorded rejected trajectories with zero accepted designs. 15 Results pages passed the headless Streamlit check. The production AlphaFold3 test and the isolated retests registered zero molecular-library result observations. The official term acknowledgment is stored per user for the current terms version. [Detailed follow-up](tool-connections.md#final-fifteen-tile-retests)

| Tile or mode | Workflow ID | Live execution and MCP Results | Streamlit Results | Earlier saved Results | Patched replay |
| --- | --- | --- | --- | --- | --- |
| AbMPNN | `sequence.abmpnn` | 1 completed | Verified | — | — |
| AlphaFold2 Initial Guess | `refolding.af2_initial_guess_input` | 1 completed, 1 failed | Verified | — | — |
| AlphaFold3 JAX · official weights | `refolding.optimized_prediction` | 1 completed | Verified | Retrieved | — |
| AlphaFold3 JAX · OpenFold3 P2 weights | `refolding.optimized_prediction` | 3 completed, 2 failed | Verified | — | — |
| AlphaFold3 PyTorch · OpenFold3 P2 weights | `refolding.optimized_prediction` | 1 completed | Verified | — | — |
| AtlasFold-M | `refolding.optimized_prediction` | 1 completed | Verified | — | — |
| AutoDock Vina | `ligand.docking.vina` | 2 completed | Verified | — | Completed |
| AutoDock-GPU | `ligand.docking.autodock_gpu` | 1 failed | Patched replay | — | Completed |
| BindCraft2 | `design.bindcraft2` | 1 failed | Patched replay | — | Completed |
| BindCraft2 Generator | `design.bindcraft2_generator` | 1 failed | Patched replay | — | Completed |
| Binderscope | `analysis.binderscope` | 1 completed | Verified | — | Completed |
| BioEmu-2 | `simulation.bioemu2` | 1 completed | — | — | — |
| Biotite Interface | `analysis.biotite_interface` | 1 completed | — | — | — |
| Boltz-2 | `refolding.optimized_prediction` | 3 completed, 2 failed | Verified | — | Completed |
| BoltzGen | `design.boltzgen` | 1 completed | Verified | Retrieved | — |
| BoltzGen Fold · Refold candidate batch | `refolding.boltzgen_fold` | 1 failed | Patched replay | — | Completed |
| BoltzGen Fold · Target PDB + binder sequence | `refolding.optimized_prediction` | 1 completed | — | — | Completed |
| BoltzGen Generator | `design.boltzgen_generator` | 1 completed | Verified | Retrieved | Completed |
| BoltzGen Inverse Folding | `sequence.boltzgen_inverse` | 1 completed | Verified | — | — |
| Caliby | `sequence.caliby` | 1 completed | Verified | — | — |
| Chai-1 | `refolding.optimized_prediction` | 1 completed | — | — | — |
| ColabFold | `refolding.optimized_prediction` | 1 completed, 1 failed | — | — | — |
| ColabFold 2.0 preview | `refolding.af3_jax_cofolding` | 1 completed | Verified | — | — |
| Complex Preparation | `ligand.complex.prepare` | 1 completed | — | — | — |
| conDitar | `ligand.generation.conditar` | 1 completed, 1 failed | Verified | — | — |
| Create centroids | `library.centroids` | 1 completed | — | — | — |
| DrugRPG | `ligand.generation.drugrpg` | 1 completed | Verified | — | — |
| ESM-1v Scoring | `analysis.esm1v_score` | 1 completed | Verified | — | — |
| ESM-IF1 | `sequence.esm_if1` | 1 completed | Verified | — | — |
| ESM-IF1 Scoring | `analysis.esmif_score` | 1 completed | — | — | — |
| ESMFold | `refolding.optimized_prediction` | 1 completed | — | — | — |
| ESMFold2 Binder Design | `design.esmfold2_binder` | 1 completed | Verified | Retrieved | — |
| ESMFold2 Fast | `refolding.optimized_prediction` | 1 completed | — | — | — |
| ESMFold2 Full | `refolding.optimized_prediction` | 1 completed, 1 failed | — | — | — |
| FLOWR.root | `ligand.generation.flowr_root` | 1 completed | Verified | — | — |
| Foldseek | `analysis.foldseek` | 1 completed | — | — | — |
| Genie3 | `design.genie3_generate` | 1 completed | — | — | — |
| GNINA | `ligand.docking.gnina` | 1 completed | — | — | — |
| GROMACS | `ligand.md.gromacs` | 1 failed | Patched replay | — | Completed |
| Interface Clustering | `analysis.interface_clustering` | 1 failed | Patched replay | — | Completed |
| ipSAE & Structure Analysis | `analysis.ipsae` | 1 completed | — | — | — |
| Kejia IDR Binder | `design.kejia_idr` | 1 completed | Verified | — | — |
| LDDM | `ligand.generation.lddm` | 1 completed | Verified | — | — |
| LocalFold | `refolding.optimized_prediction` | 1 completed | Verified | — | — |
| MaSIF-seed | `detection.masif_seed` | 1 completed | — | — | — |
| Md.Analysis | `ligand.md.analysis` | API validation failed | Patched replay | — | Completed |
| MIMOSA · BoltzGen | `design.boltzgen_mimosa` | 1 completed with rejection | Verified | — | — |
| MIMOSA · RFantibody | `design.rfantibody_mimosa` | 1 completed, 1 failed | Verified | — | Completed |
| Mol-1 | `refolding.optimized_prediction` | 1 failed | Patched replay | — | Completed |
| Nesso-1 | `ligand.affinity.nesso` | 1 completed | Verified | — | — |
| OMTRA | `ligand.generation.omtra` | 1 failed | Patched replay | — | Completed |
| OpenDDE | `refolding.optimized_prediction` | 1 completed | Verified | — | — |
| Openfold3 | `refolding.openfold3` | 1 completed | Verified | — | — |
| OpenFold3 OpenBind-0 | `refolding.optimized_prediction` | 1 completed | Verified | — | — |
| OpenMM | `ligand.md.openmm` | 1 failed | Patched replay | — | Completed |
| PandaMap | `analysis.pandamap` | 1 completed | — | — | — |
| paOPT | `ligand.generation.paopt` | 1 completed, 1 failed | Verified | — | — |
| PEP-Patch | `analysis.pep_patch` | 1 completed | — | — | — |
| Peptide Trap | `design.peptide_trap` | 1 completed | Verified | — | — |
| PeSTo | `detection.pesto` | 1 completed | Verified | — | — |
| PFM | `ligand.generation.pfm` | 1 completed | Verified | — | — |
| PGMG | `ligand.generation.pgmg` | 1 completed | Verified | — | — |
| PLIP | `analysis.plip` | 1 completed | — | — | — |
| PocketFlow | `ligand.generation.pocketflow` | 1 completed, 1 failed | Verified | — | — |
| PocketXMol | `ligand.generation.pocketxmol` | 1 completed, 1 failed | Verified | — | — |
| PoseBusters | `analysis.posebusters` | 1 completed | — | — | — |
| Prepare library | `library.prepare` | 1 completed | — | — | — |
| Protein Surface Analysis | `detection.surface_analysis` | 1 failed | Patched replay | — | Completed |
| Protein-Sol | `analysis.proteinsol` | 1 completed | — | — | — |
| Proteina-Complexa | `design.proteina_complexa` | 1 failed | Patched replay | — | Completed |
| Proteina-Complexa Generator | `design.proteina_complexa_generate` | 1 failed | Patched replay | — | Completed |
| ProteinMPNN | `sequence.mpnn` | 2 completed | Verified | — | — |
| Protenix | `refolding.protenix` | 1 completed | Verified | — | — |
| Protenix Cli | `refolding.protenix_cli` | 1 completed | Verified | — | — |
| Protenix V1 | `refolding.optimized_prediction` | 1 completed | — | Retrieved | — |
| Protenix V2 | `refolding.optimized_prediction` | 1 completed | — | — | — |
| Protpardelle-1C | `design.protpardelle_generate` | 1 completed | Verified | — | — |
| PXDesign | `design.pxdesign_generate` | 1 completed | Verified | — | — |
| RFantibody | `design.rfantibody` | 1 completed | Verified | Retrieved | — |
| RFD4-Proteina | `design.rfd4_proteina` | 1 completed | Verified | — | — |
| RFD4-Proteina Structure Prediction | `refolding.rfd4_proteina` | 1 completed, 1 failed | Verified | — | — |
| RFdiffusion | `design.rfdiffusion_generate` | 1 completed, 1 failed | Verified | — | Completed |
| RFdiffusion β-pairing | `design.rfdiffusion_beta_pairing` | 1 completed | Verified | — | — |
| RFdiffusion3 | `design.rfdiffusion3` | 1 completed | Verified | — | — |
| RFpeptides | `design.rfpeptides` | 1 completed | Verified | — | — |
| Rosetta | `analysis.rosetta` | 1 completed | — | — | — |
| RosettaFold3 | `refolding.rf3` | 1 completed, 1 failed | Verified | — | — |
| RosettaFold3 | `refolding.optimized_prediction` | 1 completed | — | — | — |
| RosettaLigand | `ligand.docking.rosetta_ligand` | 1 completed | Verified | — | — |
| SAbR Antibody Numbering | `analysis.sabr` | 1 completed | Verified | — | — |
| ScanNet | `detection.scannet` | 1 completed | Verified | — | — |
| Sequence Clustering | `analysis.sequence_clustering` | 2 completed | Verified | — | — |
| Sequence Composition | `analysis.sequence_qc` | 1 completed, 1 failed | Verified | Retrieved | Completed |
| Surf2Spot | `detection.surf2spot` | 1 failed | Patched replay | — | Completed |
| Uni-Dock Pro | `ligand.docking.unidock_pro` | 1 completed | — | — | — |
| Validation & Acceptance | `analysis.validation_acceptance` | 1 completed | Verified | — | — |

The machine-readable matrix contains the run references and result counts. [Connection matrix](tool-connection-matrix.json)
