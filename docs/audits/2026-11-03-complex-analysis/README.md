# Independent Complex Analysis audit — 3 November 2026

## Outcome and scope

Only previously undecided RFD4 candidate/run pairs received new acceptance decisions. Of 780 pairs from 192 runs, **35 passed, 727 failed and 18 remain unresolved** because required native confidence evidence is absent. The frozen pre-audit acceptance records have zero accepted/rejected overlap with this selection. No BindCraft2 acceptance decision was recalculated or changed. Additional Database evaluations on one accepted Tau BindCraft2 complex were scientific spot-checks requested in the original audit; subsequent clarification stopped further submissions.

The 35 passes are available on pc002975 as **35 canonical binder identities in 34 run-linked accepted collections**. Names begin `R3 (31 aa) · accepted RFD4` or `Tau214–283 (70 aa) · accepted RFD4`. There are 29 R3 and 6 Tau passes. The older source entries used whole-complex molecular identities. Promotion uses the binder identity, preserving the older complex identity as provenance and attaching the exact existing complex coordinates as a representation. Canonical sequence fingerprints prevent duplicate binder identities.

Every promoted record has an explicit existing target molecule ID, target sequence SHA-256, source node/run/candidate/complex-record IDs, exact complex digest and binder/target chains. Each collection preserves its source run, job reference and available campaign/parent metadata. Each filter Evaluation refers to the new structure representation. All 35 passed the real Database authorization/structure-selection boundary for Complex Analysis without launching another calculation. Sequence carriers are **not** asserted to be validated native structural references.

All 34 source runs lack a formal campaign ID in the recoverable metadata. Their recorded campaign names are preserved as labels; shared names were not used to invent parent-campaign membership. The exact run links are present. All 35 saved selections are protected from automatic run cleanup. Results navigation exposes the accepted collection for all 19 locally sourced passing candidate/run pairs, including when no automatic Database outputs remain. Remote Results require their source workstation; the newly curated collections remain local under the existing sharing configuration.

The 727 newly rejected source memberships and 248 pre-existing explicitly rejected report memberships were removed through owner/admin-checked Database operations. All 975 are absent from the source and both peers, including their group-scoped structures, observations and source records. Source-side application checked Results coordinate digests before/after and preserved every other membership. Unresolved entries remain pending. New accepted collections were materialized on pc002975; automatic bulk Molecular Database sharing remains configured off. Explicit removal events propagate independently of that setting.

## Saved policy and canvas

The inspected canvas draft is revision `92d93d…`, saved 2 November 2026, at `pipelines/761cb2224390f726b26c6230d09c2ca2ae2ede127e3a1ae78c594f2d35ef13b2/2b423c1d015f4ec7af82d2c43d763b4b/draft.json` under the Studio home. Its route is Target → RFD4 → Candidate collection → structure/ipSAE → shared acceptance filter → Campaign Stop → Results. The internal target-node name is misleading; the saved target sequence identifies PD-L1.

For the pending sequence-defined IDR candidates, the applicable saved miniprotein policy requires pTM ≥0.55, iPTM ≥0.70, native oracle iPAE/31 Å ≤0.35, at least seven binder residues contacting the target within 4 Å, and zero interchain Cα pairs within 2.5 Å. RFD4 percentage pTM/iPTM values are converted to fractions. These are the existing RFD4-compatible subset of the BindCraft2-style policy, not a claim of cross-engine calibration. An ipSAE calculation is a separate measurement; this saved policy contains no mandatory ipSAE threshold. No matrix or threshold was invented.

The registered structure-analysis workflow followed by acceptance was executed on existing outputs. RFD4's saved scalar OF3 confidence JSON does not contain a compatible pairwise PAE matrix, so ipSAE remains unavailable for those structures. Geometry is measured from the saved coordinates. Eighteen additional records outside the applicable saved-policy selection were left untouched; they were not assigned a fabricated policy.

## Independently verified campaign leads

| Source run | Target | Native production evidence | Normalized/Database ranked designs |
|---|---|---|---:|
| PC004032 `20261019-181522-21f9a96c` | R3, 31 aa | 5,000 trajectories; 9,475 refolds: 143 pass, 9,332 reject | 65 / 65 |
| PC004032 `20261021-184403-b0c849c0` | R3, 31 aa | 3,665 trajectories; 6,614 refolds: 89 pass, 6,525 reject | 35 / 35 |
| PC004033 `20261021-190947-7ea80083` | Tau214–283, 70 aa input | 3,873 trajectories; 4,799 refolds: 37 pass, 4,762 reject; user-requested early finish | 17 / 17 |
| PC004032 `20261022-224027-4c669b6e` | Tau214–283, 70 aa input | 1,636 trajectories; 1,541 refolds: 18 pass, 1,523 reject; early finish | 8 / 8 |
| pc002975 `20261102-081647-0ff0921a` | PD-L1, 115 aa | Campaign ledger generated 300; saved final policy accepted/materialized 162 and rejected/deferred 138 | 162 / 162 |

The passing refold variants have exactly 65, 35, 17 and 8 distinct trajectory hashes, respectively; every hash has a ranked result. Keeping selected redesigns per trajectory explains the difference. These were not missing Database transfers. The named Tau inputs really contain the 70-aa sequence, although saved target crops can be shorter.

For the PD-L1 run, native PPI evidence now has 118 passes and **44 failures**, with none missing. Its separate saved final policy explicitly accepts all 162 materialized candidates; it does not exclude native PPI failures. They were retained. The older 42-fail/two-missing lead was a transient snapshot, not the completed result.

R3-labeled RFD4 groups contain **636 distinct historical molecular identities: 286 with the 31-aa R3 target and 350 with the 70-aa Tau target**. The R3-named acceptance subset has 84 distinct accepted and 240 distinct rejected identities; membership counts are 84 accepted, 248 rejected and 12 unresolved. Unresolved records overlap other outcomes in different runs and must not be added as disjoint molecules.

Actual sequences:

- R3: `VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQ`
- Tau214–283: `SLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLD`

## Confirmed defects and corrections

| Severity | Evidence / correction | Commit |
|---|---|---|
| High | Replicas refused owner/admin removals; source-routed removal now retains durable tombstones and prevents forced ingestion from restoring membership. Results files and other groups are preserved. Nested explicit rejection metadata is honored during ingestion. | `9ae0500` |
| High | PAE filename discovery could associate model_10 with model_1 or select an ambiguous shared sidecar. Discovery now requires bounded/exact linkage. Multi-model Database files are refused unless an exact single-model representation exists. | `48ea325` |
| High | Native Campaign Stop silently discarded the visible ipSAE step. Unsupported native-stop/analysis graphs now fail validation rather than count premature acceptance. Full native-controller execution of this route remains unsupported. | `2a610d5` |
| Medium | Liability reuse mixed electrostatic PEP-Patch rows with hydrophobic evidence and dropped insertion identifiers. Method-specific selection and author-number/insertion parsing were corrected. | `48ea325` |
| Medium | Database staging discarded genuine source PAE. Exact authorized original-candidate/file-digest association and PAE digest identity are now retained. | `739baae` |
| High operational | A stale demo scheduler bridge claimed production jobs with the wrong Python and app home. It was stopped; queued calls now bind to the submitting executable/workspace. | `5f70100` |
| Medium | Acceptance after structure analysis lost native geometry aliases. Exact-file, exact-chain generator-parent geometry now remains available, with producing-tool attribution. | `1a1eb53` |
| Medium | Tau CIF omitted optional occupancy; ipSAE could report valid despite failed geometry. Prepared copies preserve atom order and record occupancy defaults/digests; incomplete geometry is reported separately. PEP-Patch received a writable Numba cache for peer UID mapping. | `3dda98a` |
| Medium operational | Workstation acceptance released a staging reservation without dispatching the next queued request, leaving CPU capacity idle until expensive catalog reconciliation. | `50cbc63` |
| Low | Initial collection pagination stayed disabled; numeric LocalFold choice defaults prevented canvas validation. | `bbcbd5a`, `affade4` |
| Medium | Database analysis defaulted to chain A even when the saved representation declared binder B/target A. Interface, PEP-Patch, liability and tag/fusion defaults now follow validated saved roles. Generated binders retain an exact modeled-target molecule relationship, and Database evaluations preserve that target identity for the matching chain selection. A modeled crop is not asserted to be the full campaign target or a native reference. | `cbeed24`, `8a2f52d` |
| Medium | Explicit run-linked imports into curated groups could be removed with automatic run records. Normal imports now mark these as saved selections; Results UI/API expose accessible saved collections separately from automatic outputs. Source-run and destination-collection ACLs both apply, including legitimate selections saved to another accessible project/private collection. | `cbeed24` |
| High scientific | Standalone Rosetta called total interface dSASA / bound binder SASA a binder percentage; the live Tau check produced values above 100%. That field is excluded from new normalized scores and historical score/filter projections, with raw evidence preserved. Valid Rosetta energy, area and residue-count fields remain available. Compatible cached projections upgrade without recalculating the whole history. BindCraft2 native records are unaffected. | `04983bd`, `e0a050b` |
| Medium operational | Rootful Docker created Rosetta output CSVs as root, preventing the host from finishing an otherwise successful calculation. Rosetta and ipSAE now use the existing rootful/rootless-aware UID mapping helper. A real rootful container check on PC004033 verified UID 1000 output and successful host rewriting. | `932909d` |

The four pre-existing modified files were inspected and preserved separately in `6966d17`. No Git history was reset. Code was fast-forwarded on both peers using local Git bundles under the user's later synchronization authorization; no GitHub push occurred.

## Verification and scientific limits

The broad audit suite passed **408 tests**, with 10 opt-in Docker tests skipped and 80 existing Biopython deprecation warnings. Canvas/campaign tests passed **156**. Scheduler/workflow queue tests passed **88**. The final peer staging suite passed **36**. These suites overlap; do not sum them as unique tests. Removal tests cover owner/admin, nonowner denial, unauthenticated access, disconnected source, durable retry, other memberships and Results preservation. Adversarial scientific tests include same-length wrong-model PAE and mixed-method/insertion-code patch data.

The product follow-up regression suite passed **208 tests** with 26 existing Biopython warnings; the final target-identity cases passed twice, including reversed chain roles. The projection-upgrade suite passed **60**. Both compute machines passed the **77-test** focused product suite and the final five scorer/projection tests. The final historical-chain-default suite passed **17** locally; all three existing representative Database structures reopened with the expected binder B/target A roles without performing calculations. Live verification checked all 35 exact saved structures, target links, default chain groups and cleanup-protection flags; 19 local-source Results links; and the new API `saved_collections` response. These counts also overlap.

`scientific-executions.csv` records 28 real registered Database jobs separately from tests, including failed attempts. Biotite, Rosetta, PEP-Patch, Sequence Liability and Tag/Fusion were exercised on selected passing PD-L1/R3 structures; Tau's Biotite, PEP-Patch, Liability, Tag/Fusion and ipSAE completed. The earlier Tau BindCraft2 Rosetta attempt failed because of output ownership and was not repeated after the user's scope clarification. The replacement used one newly accepted Tau RFD4 binder. Its final registered Rosetta job `analysis/20261103-190338-8acfa0a6` completed with valid scientific output, the original exact Database representation/target identity, and the mislabeled SASA ratio excluded. Tau BindCraft2 ipSAE was 0.384856, directional min/max 0.384856/0.675084, at 10 Å PAE/distance cutoffs with the exact source PAE. PD-L1/R3 RFD4 ipSAE remained unavailable. No design generation or invented reference was used.

Biotite contact distances are Å; surface areas are Å², with hydrogen-bond evidence unavailable when explicit hydrogens are absent ([Biotite](https://www.biotite-python.org/latest/apidoc/biotite.structure.hbond.html)). Rosetta interface energies retain Rosetta score units and protocol dependence, not experimental affinity ([InterfaceAnalyzer](https://docs.rosettacommons.org/docs/latest/application_documentation/analysis/interface-analyzer)). Hydrophobic and electrostatic surface patches remain separate [PEP-Patch](https://github.com/liedllab/surface_analyses) methods. [ipSAE](https://github.com/DunbrackLab/IPSAE/blob/main/ipsae.py) needs compatible pairwise PAE and chain/residue correspondence.

R3/Tau measurements describe a saved predicted conformation of sequence-defined/disordered targets. They do not establish a stable native fold or binding. Reference placement uses a checked source-to-reference rigid transform and exact residue mapping; environment overlap, tag/fusion clearance and sequence motifs remain advisory geometry/evidence. No genuine validated reference association was available for these selected IDR complexes, so reference-dependent live evaluations were not run.

Remaining limits: general replica metadata editing remains read-only outside the implemented removal path; automatic bulk molecular sharing remains disabled; optional-evidence arrivals do not automatically refresh every previously cached liability result; native Campaign Stop plus intervening structure analysis is guarded rather than fully implemented. Historical/smoke/failed runs with no recoverable sequence or generation ledger remain explicitly unverified in the inventory. Normalized output counts are not substituted for discarded/deferred generation counts.

## Evidence files

- `run-inventory.csv`: 2,410 deduplicated source runs, including historical/smoke runs, source policy, target sequence hashes, normalized/Database counts and distinct audit outcomes. Blank generated counts mean unavailable evidence.
- `pending-decisions.csv`: all 780 original candidate/run pairs, exact coordinate digests, registered filter jobs, measurements and reasons.
- `existing-acceptance-records.csv`: frozen pre-audit candidate/run decisions.
- `promoted-records.json`: target, binder, source complex, collection, structure, run/campaign and filter links for the 35 passes.
- `verification.json`: three-machine removal checks and all 35 Complex Analysis source-authorization checks.
- `scientific-executions.csv`: real Database workflow attempts and outcomes.

Raw snapshots, local backup, complete execution journals and reproduction scripts are retained in the ignored `.tmp/independent-audit-20261103/` directory. No direct SQL edits were used to mutate production membership or acceptance records.
