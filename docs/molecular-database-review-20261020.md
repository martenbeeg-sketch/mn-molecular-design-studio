# Molecular Database review — 20 October 2026

## Recommendation

Make Database the place to find, compare, select, and reuse molecular inputs and scientific outputs across Studio tools. Keep My results as the execution history, with queue state, logs, settings, and raw artifacts. A database row should follow a molecular entity through its calculations; it should not simply reproduce one candidate row from each job.

Automatically register molecular outputs and attach new evaluations to their inputs. Users should control the displayed project, groups, filters, shortlist, and exact selection sent into a tool. Group membership should not determine whether the underlying record exists.

## Evidence reviewed

- [Tamarind Database](https://app.tamarind.bio/database): inspected while signed in, including import-source filters, personal/organization scopes, project/user filters, File/Edit/View menus, and the import dialog. The inspected account views returned no populated groups, so populated-record behavior was assessed from its documentation and published preview rather than exercised against account data.
- [Database documentation](https://docs.tamarind.bio/tamarind/database).
- [Molecules API concepts](https://app.tamarind.bio/api-docs/molecules-concepts): inspected in the browser, including group semantics and computed columns.
- Studio's live Database page on the local development instance and its implementation in `run_app.py`, `workflows/modules.py`, `ligand_compounds.py`, and the ligand artifact contracts.

## Tamarind's documented model

Tamarind identifies a molecule by its chain composition, allows membership in several groups, and preserves accumulated scores on repeat import. Structures carry source labels. Design tools create output groups; prediction and scoring attach results to input groups. It distinguishes molecular ancestry from the originating job/user. Group schemas can specify fields and chain roles. These are documented semantics, rather than assumptions inferred from appearance.

The signed-in UI exposes separate filters for tool inputs, agent inputs, user uploads, and tool outputs. Its API concepts also describe computed columns for derived values and selection rules.

## Findings in the current Studio

| Finding | Evidence | Practical consequence |
| --- | --- | --- |
| Records mainly represent job candidates | `_render_database()` builds entries from `_studio_candidate_sources()` and legacy `(task_group, run_id, candidate_id)` references. Compounds use separate set/member references. | There is no common molecular identity that directly connects an imported molecule to all of its subsequent predictions and scores. |
| Automatic coverage is incomplete | `candidate_sources()` scans only `candidate-import`, `design`, `refolding-validation`, and `analysis`, requiring normalized candidate files. Ligand jobs use `ligand-docking`, `ligand-affinity`, `ligand-generation`, and MD-specific groups/contracts. | Docking, Nesso, and generated-compound artifacts are not automatically collected by this discovery path. Expanding the group list alone would still leave incompatible data formats. |
| Grouping can hide a source collection | `source_entries` excludes every source mentioned by any custom-group member. | Creating a shortlist from a job can remove the full original source entry from the sidebar. |
| Compound subsets can expand at tool launch | `_database_compound_tool_options()` filters selected members in its single-source branch but returns the original full-set path and count. | A subgroup can dock compounds that were excluded from the group. This needs an immediate correction. |
| A normal user's save can remove other users' groups | `_load_database_groups()` returns only that user's groups, while `_save_database_groups()` replaces the shared JSON file with the supplied list. | Creating/importing a group as a non-admin can drop groups belonging to other users. Ownership-preserving transactional writes are required. |
| Scores are difficult to compare | Table columns are fixed; sort/filter controls use candidate name, length, stage, and source. Metrics are displayed in the selected row's JSON detail. | Users cannot build a practical cross-tool scoring sheet or filter a shortlist by available scientific metrics. |
| Details are technical rather than visual | Selected-row detail renders sequence/SMILES, JSON, and structure downloads. | No integrated structure/pose browser or compound depiction is provided on this page. |
| Sidebar scales with job count | The live page lists many separate generation, prediction, analysis, and smoke jobs as groups. | Finding a campaign or useful shortlist becomes increasingly difficult. |
| Shared scope is unfinished | The Organization control explicitly reports that shared molecule groups are not enabled. Database controls have no project filter. | The UI does not yet provide a consistent project/team library experience. |

The subset finding was reproduced with the existing function in an isolated temporary directory: a two-compound source contained `selected` and `excluded`; a group referencing only `selected` resolved to a tool input containing both, with displayed count 2. No production records were changed.

The ownership finding was also reproduced with the existing load/save functions in a temporary directory. Starting with Alice's and Bob's groups, creating a new Alice group as a normal user left only Alice's two groups in the stored file; Bob's entry disappeared. This is a confirmed storage issue, not a hypothetical permission concern.

The existing implementation has useful foundations: immutable compound sets, integrity digests, exact job-input snapshots, protein candidate ancestry, typed ligand artifacts, and project access checks. Preserve those contracts and add a common library service above them.

## Proposed Studio behavior

### What to collect

- Imported sequences, proteins, compounds, and complexes, including inputs uploaded directly inside tiles.
- Generated candidates with explicit generation/validation/acceptance status.
- Predicted structures and docking poses as distinct representations attached to the relevant entity or complex.
- Score observations, measurements, and annotations, preserving the calculation that produced each value.
- Preparation and pocket definitions as linked context, including original receptor, chain/residue locators, repair history, pH, prepared ligand state, and coordinate transforms.

Keep raw output folders and large trajectories in their current artifact storage. Database should index them with source job, workstation, and retrievable artifact references. A preview on another workstation must resolve the source through the existing access-controlled API rather than assume a local absolute path exists.

### Identity and score context

Introduce stable molecule, representation, and observation IDs alongside the existing run/candidate IDs. Preserve complete chain composition, repeated copies, nucleic acids, stereochemistry, and explicit chemical-state distinctions. Link prepared protonation/tautomer states to their parent compound instead of silently merging them.

Store a docking score against the exact ligand state, receptor, pocket, pose, and run. Store a sequence-level property against the relevant sequence. Nesso observations need the protein–compound context. Repeated calculations should remain accessible with their settings, units, ranking direction, and model version.

For example, imported imatinib remains reusable as a compound; its Vina poses against 1IEP are related results; GNINA rescoring refers to the scored pose; and a later MD run refers to its selected prepared complex. These relationships should be inspectable from the compound and complex views.

Deduplication and shared lookup must respect project/account permissions. Cross-project identity matching must not expose inaccessible records.

### Page layout and selection

- Left: project/scope selector, All molecules, Inputs, Generated, Shortlisted, saved groups, and saved filters. Put job-created collections in a collapsible section.
- Center: a paginated molecular table with type, name, composition, tags, structure availability, chosen scores, and status. Support numerical filters and explicit score/run selection.
- Right: a detail panel with Mol*, a 2D compound depiction, chains and roles, structure/pose selection, scores, history, annotations, and source-job links.
- Selection actions: Add to group, Save shortlist, Export selected, and Use in tool/pipeline. Compatibility should determine the available tools and required inputs.

Groups should contain references. One record may be in several groups, and creating a group must leave the full source collection accessible. Distinguish a saved live filter from a frozen selection. Every submitted job consumes a frozen selection of exact entity/representation IDs and content digests.

Offer explicit policies such as a chosen run, latest successful run, or best pose within one specified protocol when displaying repeated scores. Missing data should remain missing. A generic cross-engine score should not be invented.

## Implementation order

1. **Correct selection and group writes.** Snapshot partial single-source compound groups correctly; preserve source collections when creating a group. Preserve other owners' groups and protect concurrent updates. Add focused regression coverage for exact selected membership and multiple users saving groups.
2. **Add the molecular library service.** Implement stable entities, components, representations, observations, groups/memberships, provenance, and project access. Keep durable library annotations separate from a rebuildable job index. Expose the service through the Studio API so tiles, pipelines, and UI share it.
3. **Backfill and connect job completion.** Import legacy candidate/group references and typed ligand artifacts idempotently. Attach evaluations to identified inputs. Record failed/unavailable artifacts accurately and support remote artifact resolution. Avoid rerunning scientific jobs during migration.
4. **Build the molecular browser.** Add the table/detail layout, scientific columns and filters, exact selection actions, project scope, and a common molecule picker for tiles. Use one implementation of selection resolution everywhere.
5. **Connect pipeline selections.** Allow a Database group or frozen selection to supply typed graph inputs. Write generated entities and evaluation observations back through the same service. Synchronize library metadata across workstations through the supported API/peer mechanism.
6. **Add advanced schemas and formulas.** Introduce reusable field/chain-role templates, experimental measurements, saved derived columns, and plotting after identity and result association work reliably.

Acceptance should include a round trip from import → selection → tool → recorded results → reuse, with exactly the selected compounds, preserved ancestry, multiple structures/runs, and correct project permissions. A meaningful cluster check should confirm that records remain browsable when artifacts live on another workstation.

This review changes no application behavior. Its proposed redesign is a Studio recommendation; it is not a claim to reproduce Tamarind's internal implementation.
