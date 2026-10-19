# Cofolding input editor

For project-scoped batches of filter-passing designs, target templates and result
comparisons, see [Refolding selected designs](design-validation-batches.md).

The cofolding pages share `cofolding_input.py` and the static frontend in
`mn_molecular_design_studio/components/cofolding_editor/`. The frontend uses
Streamlit's component messages and needs no JavaScript build step. It renders
only molecule types and features enabled by the selected engine profile.

One protein card is shown initially. Additional cards can be added, duplicated,
removed, or reordered. Copies, sequence formatting, FASTA import, cyclic
proteins, residue modifications, ligand CCD/SMILES input, and restraints are
converted to the existing workflow task schema in Python. Empty cards prevent
submission. Chain IDs use the same protein/DNA/RNA/CCD/SMILES order as the engine
adapters.

Boltz-2 uses a single complex editor with its optional settings in the same
component. Template files are staged into the job mount and referenced in native
YAML. A template distance threshold enables forced template restraints and
potentials. Disabling MSA explicitly writes `msa: empty`; custom A3Ms are assigned
and checked against their protein sequences. Uploaded YAML can reference
uploaded A3M, CIF, or PDB files by name. Ligand affinity is requested through
Boltz's native `properties` field. Protein interface ipSAE runs after prediction,
with the PAE file that belongs to each structure sample.

## Shared MSA preparation

All MSA-capable cofolding engines use `cofolding_msa.py`. With MSA enabled
(the default), Studio checks exact-query alignments in the configured
reference directory's `msa_cache` and legacy `boltz_models/msa_repository`.
Cropped target prefixes can reuse a longer reference alignment, cropped by
aligned columns. Discovery reads the query first; full validation is performed
only for matching files.

Missing alignments are fetched with the installed ColabFold image using
`--msa-only --af3-json`, without a GPU prediction. The full validated alignment
is saved atomically as `msa_cache/<sequence-sha256>.a3m`. Per-sequence locks
prevent concurrent jobs from fetching the same missing alignment. The job gets
its own copy capped to the selected depth; uploaded MSAs are also cached, and
a shallow upload never replaces a deeper existing reference alignment.
`artifacts/raw/cofolding_msa/manifest.json` records chain, complex, source,
sequence hash, cache location, and the exact job-local input.

Adapters supply embedded unpaired MSAs to AF3, native main-MSA paths to
OpenFold3, unpaired/paired paths to Protenix and OpenDDE, chain MSA paths to
RF3 and ESMFold2 Full, an A3M directory to Chai, and a multimer A3M to
ColabFold. Boltz uses explicit per-sequence YAML MSA paths. BoltzGen Fold
prepares MSAs for the coordinate-derived target sequence and binder; its target
A3M upload is optional. Per-chain searches provide unpaired homologs, not
inferred cross-chain evolutionary pairings. The adapters use query-only paired
rows where a native schema requires them.

ESMFold2 Fast and AtlasFold are sequence-only and do not fetch an MSA.
Explicit single-sequence settings remain available on capable engines. The AF3
and ESMFold2 family workflows share this cache and staging service, including
an explicitly selected local MMseqs fallback. Historical jobs retain their
original settings: a run submitted with `single_sequence` is reported as such;
viewing its results never fetches an MSA or implies one was used.

Seed 0 means random. The version selector reflects the installed 2.2.1 model.
Boltz-2 Server is omitted. The updated local wrapper is
`mn-boltz2-opt:2.2.1-cu130-batch2`; rebuild it from `mn-tool-containers` with
`./build.sh boltz2-opt-batch2`. All optimized engines, the AF3-family model selector, and ESMFold2 now use the
same chain cards, sample count, Optional settings, and Submit controls. Engine
settings are described by `prediction_settings_for()` and the installed engine
profile. Only supported molecule types and parameters are exposed. Boltz keeps
its native YAML, template, restraint and affinity controls. BoltzGen Fold adds
target PDB/A3M uploads and uses the shared protein card for the binder. The
optimized engine pages retain multiple independent complexes under Batch mode.

SMILES previews use RDKit from the app environment or an interpreter set through
`MN_MOLECULAR_DESIGN_STUDIO_CHEM_PYTHON`. Preview generation is optional and runs
locally. Input changes preserve the editor state across Streamlit reruns;
frontend asset changes remount the component using its saved model.

Run the focused input contract checks with:

```sh
python -m pytest -q
```


## Shared results

`cofolding_results.py` reads exact model-associated artifacts and normalizes the
small quality summaries for the shared results renderer in `run_app.py`. The
results order is quality/affinity cards, py2Dmol, colored metrics table, all-model
confidence figures, and the A3M viewer. The PAE axes are square and match the
pLDDT axes height; pLDDT uses the DeepMind palette and PAE uses dark green for low
error on a fixed 0–30 Å scale. Tables remain compact and confidence figures use
the existing scroll container.

Adapters cover AF3 JAX/Torch and RF3 JSONs, native OpenFold3 aggregated/full JSONs,
Protenix/OpenDDE per-sample summary/full-data JSONs, ColabFold rank-specific score
JSONs, Chai NPZ scores, ESMFold2 pred_rows.jsonl and NPZ arrays, AtlasFold
JSON/NPZ, and Boltz/BoltzGen NPZ confidence outputs. Missing arrays are reported
as unavailable. Arrays are never borrowed from another model. BoltzGen batch3
keeps all samples in a shared fold_out_npz archive; the copied CIF's sample
number selects its own scores and arrays from that archive.

Candidates are ranked within each complex using the engine's ranking score,
then available ipTM, pTM, or pLDDT. Separate complexes use a Complex selector.
Existing jobs are adapted at display time; new optimized jobs also save the
normalized scalar confidence metrics. A3M previews prefer the manifest's exact staged input, report its source and
row count, and show one alignment per protein chain within the selected complex.
Legacy jobs fall back to discovered staged/generated A3Ms.
