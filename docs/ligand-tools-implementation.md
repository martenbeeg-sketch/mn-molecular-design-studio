# Ligand tool tiles and pipeline design

## Current tile slice

The Studio tool catalog has standalone workflows for AutoDock Vina, AutoDock-GPU, GNINA, Uni-Dock Pro, RosettaLigand, Nesso-1, LDDM, PocketXMol, OMTRA, FLOWR.root, conDitar, conDitar + paOPT, DrugRPG, PFM, PocketFlow, PGMG, and bound-ligand MD with GROMACS or OpenMM. Docking and Nesso accept a compound set selected from Database. Every generator produces a compound table and a versioned compound-set artifact; models with coordinate outputs also produce SDF poses. The Database importer accepts generated compound tables as new immutable sets. PGMG uses a PGMG-compatible `.posp` pharmacophore input; OMTRA can additionally consume Pharmit JSON; paOPT exposes its six supported ADMET endpoints.

Each submitted job snapshots its input compound subset under `artifacts/input/compound_set.json` and records the source set ID, source digest, consumed digest, and count in `artifact_contract.json`. Docking results distinguish scores from structures. GNINA's CNN score is explicitly higher-is-better; Vina, AutoDock-GPU, and Uni-Dock scores are lower-is-better. AutoDock-GPU prepares a receptor with Meeko, computes AutoGrid maps, then sends each ligand and the `.maps.fld` descriptor to the packaged AutoDock-GPU binary. Nesso emits affinity metrics only and does not claim a docked pose.

The MD tiles take a protein–ligand complex PDB, a selected bound-residue key, and matching ligand coordinates in SDF. They check the SDF's heavy-atom count and element order against the selected PDB residue. GROMACS uses the bundled AmberTools/GAFF2 preparation and Roe–Brooks equilibration/production runner, returning the prepared system, XTC trajectory, checkpoint, and trajectory analytics. OpenMM uses the bundled OpenFF/Amber preparation and returns its production trajectory and run summary. The OpenMM tile explicitly selects the `current_staged` preparation protocol; Roe–Brooks is selected by the GROMACS workflow. Both runtimes are included under `mn_molecular_design_studio/engine/workflow_engine/vendor/ligand_runtime`; neither requires a sibling source checkout. Engine container images, MODELLER, force fields, and model/reference data remain runtime prerequisites.

These are queueable container workflows. Container execution has been tested through mocked runners and Studio tests; scientific execution still depends on the relevant installed Docker image and reference data on the host.

## Shared docking input page

AutoDock Vina, GNINA, Uni-Dock Pro, and RosettaLigand share one form based on Tamarind's live Vina input page. It accepts PDB/CIF uploads and RCSB accessions, SDF or SMILES/CSV ligand uploads, and existing Database compound sets. Uploaded ligands are staged as the same immutable compound-set contract used by pipelines. The bundled offline example loads 1IEP and imatinib (PubChem CID 5291), including the reference 20 Å box; its job name follows the selected engine.

The compact Mol* viewer shows a magenta search box in the receptor coordinate frame. Manual box edits update the shape without reloading the protein or resetting its camera. Pocket modes include a bound ligand with padding, a chain-aware click-and-drag sequence picker, typed author-residue ranges, and Select all/Clear controls. Selection in the sequence and structure stays synchronized. Vina and GNINA use the displayed center and dimensions as their search box, and Uni-Dock Pro uses its center, dimensions, and selected GPU search mode. RosettaLigand uses the displayed center with its own Rosetta protocol padding; the box outline is a visual pocket aid for that engine.

Each engine keeps its own search controls inside Optional Settings: exhaustiveness and pose count for Vina/GNINA, fast/balanced/detailed search for Uni-Dock Pro, and protocol, placement mode, ligand preparation, replicas, and CPU workers for RosettaLigand. GNINA exposes its CNN scoring mode; results rank by CNN score except when CNN scoring is disabled, when docking energy is the ranking metric. Uni-Dock Pro requires a CUDA GPU. Rosetta scores are in REU and are intended for comparison within a protocol.

The Boltz pocket mode submits 1–10 independent probe-SMILES tasks through `refolding.optimized_prediction` using Boltz-2. Refreshing a completed prediction reads its normalized candidate structures, detects probe contacts, aligns predicted protein sequences back to the uploaded author residues, and derives the Vina box from those original receptor coordinates. This prevents a cofolded coordinate frame from being applied directly to an unrelated receptor frame. Probe runs still require the Studio Boltz image and reference resources; unit/UI checks validate request construction and contact mapping without running model inference.

## Database model

- A compound set is a versioned JSON record with normalized compound IDs, SMILES, retained source properties, a stable set ID, and a SHA-256 content digest.
- Database groups refer to compound-set records and compound IDs. A launch resolves that selection to an immutable snapshot; it never depends on later changes to the source group.
- A generated set can be imported from the generation tile's CSV/SDF outputs. Import creates a new Database record and does not mutate the job output.
- Future compound annotations should be separate versioned records keyed by `(compound_set_id, compound_id, annotation_type, source_job_id)`. This keeps docking scores, affinity estimates, and MD readouts attributable to the method and run that produced them.

## Pipeline direction

The current visual graph schema and executor are centered on protein structures and sequence candidates. Ligand tools have typed job artifact contracts, but no ligand graph nodes or artifact-edge resolution have been added yet. Add those through a ligand-specific typed port model rather than treating every output as an untyped result folder.

Proposed typed flow:

1. **Target and pocket**: a protein structure plus a pocket definition (co-crystal ligand coordinates, residue selection, or a validated center/box).
2. **Compound source**: a Database compound-set reference or an upstream generated compound set.
3. **Generation**: LDDM consumes a target and reference ligand; PocketXMol, conDitar, DrugRPG, PFM, and PocketFlow consume a target and pocket; FLOWR.root uses a pocket plus reference ligand; paOPT steers conDitar with selected ADMET objectives; OMTRA accepts a target with optional pocket/reference context; PGMG consumes a PGMG `.posp` hypothesis. These standalone tiles emit `compound_set/v1` and retain coordinate outputs when the model produces them.
4. **Docking**: Vina, AutoDock-GPU, GNINA, Uni-Dock Pro, and RosettaLigand consume a target/pocket and compound set. They emit score rows and pose collections, with engine-specific metrics and ranking direction or units in the artifact metadata.
5. **Selection and rescoring**: a filter node selects compounds/poses using an explicitly named metric and direction. Nesso can consume the selected target/SMILES set for affinity prediction, but its output is not a structure and cannot directly feed a pose-dependent step.
6. **Molecular dynamics**: the standalone OpenMM and GROMACS tiles accept a protein–ligand complex. A pipeline still needs a pose-to-complex adapter for docking outputs, followed by separate typed preparation, production-replica, and trajectory-analysis stages. That preserves the bundled runtime's preparation and restart boundaries.
7. **Promotion**: a user-selected generated or filtered molecule collection is imported into Database as a new immutable compound-set record.

## Recommended implementation sequence

1. Add a Database compound-set source node and ligand target/pocket node to the graph catalog; validate ports and project ownership before scheduling.
2. Add the standalone docking, affinity, generation, and MD tiles as graph nodes that reuse their workflow functions and artifact contracts, including the newly added generator adapters.
3. Add a score-aware selection node. Never compare values across engines without an explicit normalization or ranking policy.
4. Add an explicit pose-to-complex adapter so the current PDBQT pose outputs can become validated complex inputs for MD.
5. Expand the GROMACS tile into a multi-stage graph workflow for system preparation, production replicas, and analysis, while retaining the current standalone tile as a useful first pass.

This sequence keeps the standalone tiles useful immediately while making each pipeline edge auditable and reproducible from the exact Database record and job artifacts.
