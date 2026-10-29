# RFD4–Proteina implementation and testing handoff

Updated 22 October 2026.

## Current status

Ready for user acceptance testing. This is not a declaration that every task,
combination of controls, input chemistry, and folding engine has completed a GPU
run. Start with one design and inspect its native metrics before larger runs.

The design page exposes Tamarind's 13 task choices and task-dependent controls.
Nanobody frameworks appear before a structure is uploaded. Selecting a supplied
nanobody exposes design-chain and CDR selectors. Existing Structure switches
between diversification, inverse folding, structure prediction, and side-chain
packing. The form uses the reference's field order, rectangular task selectors,
green selection state, stacked fields, and optional-settings section within
Studio's own navigation. It is not a pixel-identical copy of Tamarind's application.

### Controls connected to native conditions

- All four bundled VHH frameworks and a supplied nanobody; per-loop CDR length
  ranges; framework backbone/all-atom conditioning; separate antibody/TCR chain
  docking; constant-domain trimming with ANARCII antibody or TCR numbering.
- Indexed and unindexed motif placement, per-gap scaffold lengths, enzyme total
  length bounds, and per-residue active-site atom overrides. Explicit total
  length bounds require an unindexed motif; incompatible combinations fail with
  an explanation.
- Ligand selection, pose and atom restraints, burial/exposure, hydrogen-bond
  donors/acceptors, pi interactions, salt bridges, hydrophobic contacts and metal
  coordination. Crystallization additives GOL, SO4, PEG, PGE, EDO and MPD are removed.
- DNA sequence/structure conditioning; RNA dot-bracket conditioning.
- Diversification by trajectory fraction or independent backbone/sequence times;
  fixed sequences/coordinates, a redesign shell and secondary-structure regions.
- Sampling effort, step count, ODE sampling, noise and TSR settings, seed zero,
  sequence preferences, reward search, optional refolding, pass thresholds,
  apo checks, interface energy, developability, trajectory output and clustering.
- Annotated CIFs, native condition files, and batch structure scoring. Author
  chain/residue IDs are mapped separately for each uploaded structure.
- The experimental IDR peptide-binder form accepts one IDR FASTA sequence and
  creates an ideal extended backbone CIF as a sequence/atom carrier. Its native
  condition spec fixes `C_SEQ`, adds optional residue-indexed `C_HOT`, and omits
  `C_CRD` and `C_CRDGRP`; the target coordinates remain free. This is a flexible
  co-generation hypothesis, not an IDR ensemble prediction or validated binder
  campaign. Hotspot positions are 1-based in the pasted FASTA.
- Structure scoring can select the matching protein binder, nanobody, antibody,
  TCR, motif, enzyme, small-molecule, DNA, or single-protein criteria. Enzyme
  scoring now selects RFD4's no-ligand pass criteria when all inputs lack a
  ligand; mixed ligand/no-ligand enzyme batches must be split. It retains the
  uploaded RFD4 annotations required for motif metrics and chain-of-interest
  detection.

RFD4's native sequence-conditioned protein folding is a separate tool tile and
pipeline node. Independent refolding uses a selectable compatible Studio engine.
The native OpenFold3 bridge uses NVIDIA's matching published runtime, not Studio's
incompatible general-purpose OpenFold3 container. Alternative engines use the
shared Studio prediction workflow and must return mappable confidence arrays;
missing confidence is an error, never fabricated.

The visual pipeline starter menu now has RFD4 protein-binder, small-molecule
pocket, and enzyme-motif workflows. They use the Python API tutorial's 5ENI/5QA
selections for the two ligand examples and continue to Boltz-2; the binder path
also adds ipSAE and Rosetta interface analysis. These graph starters have not yet
been submitted as pipeline jobs or exercised with new GPU runs. The Hugging Face
ZeroGPU implementation additionally uses AtomWorks Studio for direct annotated-
CIF editing and has a SecureDNA screening sidecar; those deployment features are
not currently part of Studio's local RFD4 runner. The Space's AOT build pins model
revision `388fcba` while Studio remains pinned to the newer local bundle
`a06e250`; the new templates borrow the Space's interface/workflow idea without
changing Studio's model pin.

## Evidence from this workstation

- 57 focused Python tests passed (request validation, registration, pipeline
  contracts and capabilities). One existing FastAPI/httpx deprecation warning.
- 17 native input cases passed, covering the task types, the four structure
  operations and both custom input methods. Unconditional generation was also
  checked on the GPU.
- 17 additional native control cases passed: four VHH frameworks × two atom
  selections with CDR ranges, supplied antibody/TCR loops, real antibody/TCR
  variable-domain trimming, indexed/unindexed scaffold gaps, diversification
  controls, ligand contacts and active-site overrides.
- 21 native configuration compositions passed, including reward-search variants.
- Browser checks exercised all task forms, all five nanobody framework choices,
  operation switches, independent noise controls and object editors.
- A 30-residue unconditional design plus independent native OpenFold3 refolding
  completed; native backbone RMSD 0.8439 Å and pLDDT 80.25. Run:
  `design/20261022-101323-51c6e425`.
- Native RFD4 protein structure prediction preserved the supplied 30-residue
  sequence. Run: `design/20261022-094449-2cf7215f`.
- Native RNA generation produced ten nucleotides (`GGCUGAAGCU`) and normalized as
  an RNA molecule. No independent RNA refolding claim is made.
- Structure scoring with refolding disabled completed and normalized its output.
  Run: `design/20261022-104702-99c53243`.
- A Boltz-2 independent folding call completed through the shared prediction
  bridge and returned a mapped structure, pLDDT and PAE.

Local detailed logs and native input/configuration reports are in
`/tmp/mds-rfd4-validation`. This temporary directory is not a durable audit store.
Actual workflow runs and their staged adapter/configuration/provenance files are
in Studio's configured results directory.

The score-type and enzyme no-ligand follow-up changes described below were made
after that validation pass; they have not yet been rerun through the focused test
suite or a GPU workflow.

## Limits requiring further testing

Full GPU generation remains to be exercised for binder, antibody/TCR, enzyme,
small-molecule, DNA and motif tasks, including combined advanced options. Large
batches, multi-structure specificity checks, clustering, energy scoring, search
and every alternative predictor are not comprehensively validated. A successful
input/configuration check does not establish biological validity or model quality.

The upstream agent skills also define campaign ledgers, immutable pass-criteria
versions, recoverable run reconciliation, campaign packaging, and reward-guided
wave experiments. Studio currently records per-run artifacts and pipeline
lineage, but does not implement that full campaign-ledger lifecycle or the
autoresearch/commitment-curve workflows. Reward search is limited to three fixed
presets and one terminal checkpoint; custom reward terms, gates, and checkpoint
schedules are not exposed. The score form applies one design type to its whole
batch; it does not yet expose RFD4's per-file task CSV mapping or arbitrary
`metrics_tags`/custom-metric configuration.

Use mmCIF for guided enzyme, ligand and DNA tasks. PDB inputs cannot retain custom
bond orders/charges and are rejected for those tasks. The general prediction bridge
requires known CCD ligands; use the native OpenFold3 bridge for custom chemistry.
Some predictor confidence tokenizations cannot be mapped safely for atomized
ligands and will fail explicitly. ANARCII must recognize a variable domain when
constant-domain removal is requested.

## Pinned runtime

Generator: package 0.0.3, upstream commit
`9d35ee4f6328fd9b17f5b84756ee5db99d6d89ac`, image
`mn-rfd4-proteina:0.0.3-cu130`.
Model bundle revision: `a06e2509c7033b010fc43614e9fd840303588502`.
Native oracle image: `mn-rfd4-oracle:598acab-cu13`, imported from NVIDIA's runtime
at revision `598acab93d5d0aeacc4362f6e88dc20e7db85ef6` after checking SHA256
`7bbd35478ed5e95131efa7756d03a2e754f3508e603ad08bb8f8db0bd89484a1`.
The published XLA flags are applied to avoid excessive first-run GEMM compilation.

To install that oracle image on another workstation, after building the generator:

```bash
python scripts/rfd4-proteina/setup_oracle.py --reference-dir /mnt/db/reference_files
```

The setup requires Docker, `unsquashfs`, sufficient disk space and network access.
The manual equivalent was performed here; the convenience script itself has not
been exercised from an empty machine. No host Python dependencies are installed.
