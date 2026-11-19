# Datasets and scientific results

## Separate purposes, one SQL model

- Datasets registers and browses reusable inputs, including structure-free
  compounds and sequences. The visible Import data action uses the existing
  reviewed import wizard. Project references are optional; new input groups use
  an empty owning project ID. Existing groups, ownership and project access are
  not automatically migrated or broadened.
  Only input members and an explicit project-reference summary are shown here;
  calculation scores, poses and attached analyses stay in Data browser. Connected
  result buttons navigate to that separate page rather than embedding results.
- Data browser explores scientific outcomes, their structures, native scores,
  parameters and attached analyses. The former Database route remains an alias.
  Input/reference placeholders are not counted as scientific outcome rows.
- Dataset members and connected results retain exact molecular identities.
  References do not copy molecules. Creating an input collection from results
  does not silently freeze every model or copy all evaluations; original saved
  results remain accessible through the member/result handoff. Input structure
  imports register their own managed, verified coordinate records.

## Project references and safety

The additive `dataset_project_links` SQL junction has an indexed project lookup.
Its portable source is the versioned dataset-registration metadata in each
collection origin; SQL triggers maintain the index after insert, update and
delete, including normal peer replay. Existing databases without this table can
still be read without mutation; the normal writer installs the additive schema.
Neither molecular identities nor scientific observation schemas change.

Only a dataset owner or administrator on the source workstation can change
references. Referenced projects must exist, be accessible to the actor and be
mutable. Reference forms compare the complete source revision under a write
transaction. Deletion, peer authority, old-writer origin edits, cleanup locks and
stale forms fail closed. References are organizational links, not ACL grants;
private datasets remain private, and existing project-shared datasets retain
their original project access. Legacy source-project deletion/collection cleanup
is blocked while another project references its dataset. Independent datasets
have no owning project for project cleanup to remove.

Inputs have an explicit registration role and per-membership input marker. Even
if imported annotations mention a generator or a pass, they are references,
not admitted design outcomes. Existing accepted designs and native filters are
unchanged. An input preview reads only its own membership, annotations, sources
and coordinates; it does not hydrate scientific observations/evaluations or
borrow another group's scores or model. Canonical dataset-entry/target names and
aliases remain searchable. Protein Library homonyms get display-only project
context without renaming source records.

## API and verification

Authenticated clients can list/create datasets at `/api/v1/molecules/datasets`,
edit project references with revision checks at
`/api/v1/molecules/datasets/{group_id}/projects`, and use the existing imports
endpoint with `reusable_dataset=true` and optional `project_references`.
Existing imports retain their previous defaults. The UI and API use the same
registration/input-marker/reference functions. Stale references return 409;
inaccessible datasets remain inaccessible even to members of a linked project.

Isolated tests cover reviewed imports, private/multi-project reuse, row selection,
structure-free pagination, connected-result links, read-only previews, input-only
exports, canonical name/alias search, SQL projection/restart/peer deletion,
stale/old-writer protection, authenticated HTTP access and project-cleanup safety.
No research analyses are started and no production credentials are reset.
The separate, previously reproduced naming-replay SQLite-lock failures described
in the workspace consolidation report are not modified by this change.
