# Canonical project and target names

Project and target labels are metadata, not scientific identity. Renaming must
never change a molecule ID, sequence, chemical identity, Evaluation selector,
admission decision, structure, job reference, project ID or filesystem path.

## Store and compatibility

`entity_names` in `database/molecular_library.sqlite3` stores the current label,
previous-name aliases, opaque edit revision, editing account and source node.
`entity_name_details` stores allowlisted descriptive properties using the same
edit revision and authority as the name. `entity_name_history` retains applied
name/property revisions. SQL labels are resolved
by stable IDs; historical scientific metadata is not rewritten.

- Project labels use the stable project ID. Existing organization JSON is seeded
  into SQL on first organization load/save. Before an explicit rename, legacy
  project-name edits remain compatible. After an explicit rename, SQL is
  authoritative: stale JSON cannot undo it. Organization JSON remains the
  compatibility projection and store for ownership and lifecycle. Once edited
  through Properties, the SQL project description is authoritative too.
- Target/ligand labels use the existing browser's exact component identity:
  molecular type plus sequence, stereospecific SMILES, or CCD identifier.
  Renaming does not assert a new target role or create a scientific connection.
  A label applies within a project, or within one owner's personal library.
  This scope prevents a private alias from being exposed through an otherwise
  identical public molecule. Different sequence crops and stereoisomers stay
  distinct. Conflicting names across accessible projects are shown together,
  not silently resolved alphabetically; selecting a project disambiguates them.
- Targets need not have standalone records: a receptor component in an imported
  complex can have a shared label. Generated binders are not promoted to targets.
- Original collection names remain unchanged and searchable. Canonical aliases
  are also searchable in the folder browser and SQL molecule search. Display
  resolution applies in both browsers and the existing record inspector.
- Read-only browsing of an older SQL snapshot without these tables still works
  and does not migrate it. Explicit rename actions use a writable connection.

## Editing

**Project properties** on Projects and **Folder properties** in the browser
open the same Properties dialog. Selecting a folder or row never starts editing.
The left tree and filters contain no editing forms.

**Result / molecule properties** opens an object selector distinguishing a saved
outcome, its target/compound identities and each accessible dataset membership.
The dialog explains the exact scope before **Edit metadata**. Only the owner or
an administrator may save, subject to the metadata authority, group source and
project cleanup locks. Downstream job ownership alone does not grant target
rename permission. Production sessions are revalidated during dialog actions.

Typed fields are deliberately limited:

- Project: name, aliases, description, tags, notes.
- Dataset/collection: name, aliases, description, source annotation, citation,
  tags, notes. Dataset metadata uses source workstation plus original collection
  ID, including on peer replicas; similarly named or locally equal-ID datasets
  on different machines do not share annotations.
- Target/compound: shared name/aliases, description, reference identifier
  annotations, tags, notes, scoped to one project or personal library.
- Dataset entry: local label/aliases, dataset identifier annotation, supplier,
  batch/sample, tags, notes. No other membership or compound identity is changed.
- Saved outcome: label/aliases, description, tags, notes. This does not rewrite
  any generating tool name, score, filter decision, model or Evaluation.

Root categories, activity/method branches and other virtual folders are not
editable objects. Compound and entry editors are separate even when the same
compound appears in several datasets. Scientific identifiers, sequences,
SMILES, coordinates, links, parameters and computed scores are not writable
through these dialogs. Existing collection membership, imports and stored
source metadata remain unchanged. The old generic scalar-metadata form and its
automatic first-membership selection have been removed.

Each form submits the revision it originally displayed. SQL checks it under a
write transaction. Another edit, a restart, or a stale form cannot silently
overwrite saved metadata. Source annotations are also checked for legacy-writer
changes before saving dataset entries or project details. **Reload current
values** discards the stale draft; **Cancel editing**, **Close**, and dismissal
write nothing. Opening or saving a dialog retains browser navigation and
selection; successful saves invalidate the SQL-backed browser cache.

## Workstation synchronization and rollout

`entity-names.v1` is an optional peer replication capability. Updated peers
exchange `entity_name` records, preserve their original authority, and do not
turn absent remote labels into local deletions. Name synchronization is disabled
with molecular-library sharing. Separate feature cursors allow a peer upgraded
later to catch up on labels that were excluded while it was older.

`entity-properties.v1` adds typed properties and dataset/entry/outcome labels.
Peers advertising only `entity-names.v1` receive compatible project/target name
records without the new properties field, not unsupported entity types. A later
upgrade uses a fresh feature cursor. Name-only replay cannot erase properties
already stored on an updated peer. All nodes still need upgrading for uniform
editing and display; this compatibility bridge is not permission to deploy a
mixed release indefinitely.

Each explicitly named entity has one source workstation. Subsequent edits must
be made there, including by administrators. A conflicting first claim from a
different, disconnected workstation is retained as a replication conflict;
it is not silently accepted as a replacement. This is a deliberately
single-writer metadata model, not distributed consensus. Establish shared names
while peers are connected before editing them independently. Conflict/authority
transfer requires operator review; no automatic takeover or name deletion is
provided.

Older peers do not receive unsupported record kinds. They can display the
project-name compatibility projection, but cannot enforce canonical target
labels or revision checks. Upgrade all UI/API/agent processes before promising
uniform naming. Do not deploy mixed source files into an old release. Back up
the SQL database and organization JSON before deployment. The migration is
additive and does not rewrite scientific records or existing research names.

## Verification

`tests/test_entity_properties.py`, `tests/test_entity_names.py` and
`tests/test_entity_name_replication.py` cover
SQL name/alias resolution, unchanged scientific rows, stale edits, restart and
reimport, permissions, private names, project scope, exact sequence crops,
cache invalidation, old read-only snapshots, project JSON compatibility,
cleanup guards, form revisions, three simulated workstations, relay ownership,
capability filtering and authenticated rename API conflicts. Existing identity
and browser suites cover molecular type and stereochemical separation.
