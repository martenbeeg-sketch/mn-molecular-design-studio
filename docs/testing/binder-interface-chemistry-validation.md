# PPI Interface Chemistry comparison with BinderScope

Validation performed on 2026-11-08 against the local BinderScope MIT source at
`/home/user/programs/git-projects/tools-to-implement/binderscope`, revision
`c61ef9c`. The comparison called BinderScope's `binder_review.interface`
implementation directly; it did not install or modify BinderScope and did not
change any Molecular Database or source structure.

## Definition comparison

| Evidence | Studio `studio.ppi-interface-chemistry.v1` | BinderScope | Interpretation |
|---|---|---|---|
| Buried surface | Bio.PDB Shrake–Rupley; 960 sphere points; 1.4 Å probe; sum of target and binder SASA loss is the two-side `interface_total_delta_sasa_angstrom2`; separately reports half that total | Same two-side dSASA convention, same 960-point/1.4 Å setup in the examined implementation | Same coordinates and selected chains produce the same value. The half-delta field is explicitly not the two-side total. |
| Buried unsatisfied polar atoms | N/O polar heavy atoms, burial ≥0.5 Å², complex SASA ≤0.1 Å², and no qualifying intracomplex H-bond | Same geometric burial/complex-SASA/H-bond rule | This is a geometric flag, not a desolvation energy; no protonation, solvent network, or energetic penalty is modeled. |
| H-bond evidence | Explicit H where present; otherwise constructed backbone amide H and a labeled heavy-atom distance/angle proxy; evidence basis counts are retained | Explicit/computed geometric donor-acceptor checks; reports a coarser “computed” basis | Counts agree on the checked structures. Studio exposes the evidence basis more explicitly. Neither value is an energy or proof of experimental bonding. |
| Hydrophobic burial | Buried carbon/sulfur SASA divided by two-side dSASA (`hydrophobic_cs_buried_fraction`) | `hyd_frac` uses carbon/sulfur buried SASA divided by dSASA | Same definition. BinderScope `Interface_Hydrophobicity` is a different statistic: hydrophobic-residue percentage in its 4.0 Å interface-residue set; Studio does not label its atom-area fraction as that statistic. |
| Residues | Author chain, residue number, insertion code and a viewer selector are saved per residue | Structure parsing retains insertion codes; some interface summary rows reduce position to an integer | Studio preserves author insertion codes in its residue evidence. The regression fixture checks `1A` separately from `1`. |
| Missing H | Uses explicit H when available; otherwise backbone geometry and side-chain heavy-atom proxy; records the split | Similar geometry fallback; summary basis is less detailed | Missing H does not silently turn H-bond count into an experimental or energy measurement. |
| Input scope | Selected protein chains only; rejects a selected component with no supported amino-acid residues; multi-model files must be explicitly resolved to one model | The cited interface routine is likewise for protein complexes | Protein–ligand cases are unsupported by this PPI evaluator and must use ligand-aware analyses. |

## Real saved structures

Two already evaluated staging structures were compared. No predictions or
analyses were rerun for this validation.

| Staging structure | Same-input two-side dSASA | Buried unsatisfied polar atoms | H-bonds | Contact residue pairs | C/S buried fraction |
|---|---:|---:|---:|---:|---:|
| RFD4 `mol-6d5d88c26498f1b6fd6a788792a4ae08`, `structure-f7db4f80ee4ec9f4db908f7c9d7ce355` | 2848.392 Å² | 2 | 17 | 83 | 0.678902 |
| BindCraft2 `mol-15c613d2ca9ad9f83fba2aedb7f42c56`, `peer-structure-726b3f36b76f240713dac46314712f7e` | 3177.605 Å² | 4 | 18 | 93 | 0.639529 |

For both structures, running Studio and BinderScope over the exact same
coordinate representation gave identical values for the listed quantities.
The RFD4 saved mmCIF also matched directly: dSASA 2848.392 Å², BUNS 2,
H-bonds 17, 83 contact residue pairs, and C/S fraction 0.678902.

The BindCraft2 original mmCIF omitted `_atom_site.occupancy`: the standard
Bio.PDB parser raised `KeyError`, and BinderScope's Gemmi loader reported no
models. Studio's parser accepts the optional field with the CIF default
occupancy. For a comparison only, the same 1,126 selected protein atoms were
written to a temporary PDB in `/tmp`; atom identity was unchanged and the
largest coordinate shift from 3-decimal PDB serialization was 0.000822 Å.
Comparing Studio's original-CIF saved result to BinderScope on the temporary
PDB gave ΔdSASA = +0.589 Å² (0.0185%), ΔC/S fraction = −0.000039, and identical
BUNS, H-bond and contact counts. Studio and BinderScope matched exactly when
both read the temporary normalized PDB. The small dSASA difference is therefore
coordinate serialization/format handling, not an observed algorithmic
definition discrepancy. The original CIF was left unchanged.

Both structures have the exact same recorded 31-residue target sequence
`VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQ` on chain A, with the selected binder on
chain B. The existing compatibility check passes. A target-only C-alpha fit
maps all 31 residues and reports 9.446 Å RMSD (`Bio.PDB.Superimposer Kabsch fit
on mapped C-alpha atoms`). That is a display alignment between these saved
models, not a native-state accuracy measure. Pooled contact evidence still
requires compatible saved Biotite settings and successful residue mapping.

## Focused fixtures

`tests/test_binder_interface_chemistry.py` now covers separated chains with no
contacts (empty evidence and `None`, not zero, for an undefined buried-area
fraction), absent hydrogens (explicit fallback basis counts), insertion-code
preservation, a ligand selected as a PPI partner (rejected), exact digest
validation and multi-model rejection. The result remains a versioned,
method-specific Evaluation with exact structure/model provenance.

## Staging cohort review and export

The existing Database structure loader and Evaluation row projection were used
for a bounded RFD4/BindCraft2 review. The two rows share the exact target above.
The RFD4 structure is from `design/20261030-070659-2044e3d6` on PC004033;
BindCraft2's recorded source is `analysis/20261031-090826-6bf41866` on
PC004032. The latter source record does not establish a final Studio acceptance
decision, so it remains a native-stage pass only.

RFD4's saved rows include PPI Interface Chemistry, Biotite Interface,
Protein-Sol, Sequence QC, and Validation & Acceptance. BindCraft2 has PPI
Interface Chemistry, Biotite Interface, Protein-Sol, Sequence QC, and a
BindCraft2 evaluation. The comprehensive profile separately reports missing
PEP-Patch, Sequence Liabilities, Rosetta Interface, ipSAE, and Tag/Fusion
Suitability evidence; the export does not mark those methods as failed or
measured zero. Each saved PPI Chemistry row exports `studio.ppi-interface-
chemistry.v1`, its exact method parameters, structure digest, scalar metrics,
and its per-residue annotations, contact pairs, and H-bond evidence.

For a workflow-path demonstration in the isolated staging database only, the
review annotation for RFD4 was set to “Shortlisted for experimental testing”
with a rationale citing its recorded final pass and saved assessments while
stating that this is not an affinity claim. BindCraft2 was set to “Not
shortlisted for experimental testing” for this demonstration because only a
native pass is present and there is no final Studio decision; its rationale
explicitly avoids interpreting that as evidence of poor design. These are
account-scoped human annotations, not acceptance decisions or Evaluations.
No production Database or scientific Evaluation was changed. The generated
staging exports are `/tmp/mds-phase4-staging-design-review.json`,
`/tmp/mds-phase4-staging-design-review.csv`, and
`/tmp/mds-phase4-staging-experimental-shortlist.csv`.
