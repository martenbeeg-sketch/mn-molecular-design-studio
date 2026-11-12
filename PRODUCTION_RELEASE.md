# Production release — 12 November 2026

## Decision and access

**The corrected research application is deployed on all three workstations with
real instructor authentication and preserved research ownership.** This is an
engineering release acceptance, not a new independent audit or certification of
every advertised scientific engine. Untested predictors, generator modes and
specialist pipelines remain unqualified.

- Research: **https://pc002975.marionegri.it:18503/** or
  **https://172.18.24.170:18503/**. The resolvable short hostname `pc002975` also works.
- Hosting PC only: **http://127.0.0.1:8502/**.
- Login: **`research_instructor`**, using the current password chosen by the
  instructor. Original bootstrap credentials are obsolete. Account recovery
  issues a new private temporary password and requires another password change;
  do not share credential files or the instructor password with students.
- Existing research owner `00000000000000000000000000000001` is retained. No
  campaign, import, project, library or Evaluation ownership was reassigned.
- Production is **instructor-only** pending the unanswered choice about three
  separate student accounts. The isolated student pilot remains at port **18502**,
  separate from research data and GPUs. Do not share the instructor password.

This report supersedes earlier statements that research still runs `25e003c`
or that AUD-008 remains page-local. Historical audit/test reports and their
scientific limitations are retained; earlier failed probes are not recast as passes.

### Client trust and network boundary

Institute IT/the operator must trust the **research** CA on actual client PCs.
Do not bypass browser certificate warnings. This is a different root from the
teaching-pilot CA. Public certificate:
`/home/user/Downloads/MDS-Research-LAN-CA.crt`.

SHA-256 certificate fingerprint:
`E5:D8:88:E9:E4:12:C1:64:50:25:31:57:D5:37:81:FC:C0:2B:CF:50:C0:64:5D:89:36:17:65:68:22:89:12:F6`.

The dedicated root is name-constrained to the research hostname/short hostname
and `172.18.24.170`; it expires 11 November 2031. Server certificates last 90 days.
`mds-research-tls.timer` checks daily and renews within 21 days of expiry without
changing the root. Renewal uses a process lock, verifies hostname/IP and Caddy
configuration, and restores the previous certificate/key pair on failure. A
forced real renewal and normal no-op check pass. Failed systemd units require
operator attention; no external paging/monitoring service is configured.

The listener binds only `172.18.24.170:18503`, allowing the existing explicit
institute subnet `172.18.24.0/22`. Other source networks receive 403, including a
loopback caller forging `X-Forwarded-For`. Plain HTTP receives 400; LAN `/api/*`
receives 404. Research APIs remain loopback-only behind trusted SSH peer links.
Gateway uploads are capped at 256 MB, with application limits still applicable.
IT must confirm client addresses belong to this subnet. No public-Internet
listener or automatic client trust installation was added.

## Database layout correction — 12 November 2026

The current deployed application is **`7a82f70af11f91f109777bee639e5f32825ab198`**
on all three workstations, from
`~/programs/deployments/mn-molecular-design-studio/prod-7a82f70af11f/source`.
Archive SHA-256: `a6b7c0898b7a08b0f3ae9eb9a60b3e549b41b60a817e8ef0e6eb5b4f69ca4b3b`.
All agents report code fingerprint
`57ddf999e716986ddbbc17c98f3e840c0ef200a7054511fd533fef272d999001`.

This corrects the preceding grouped views' unwanted page redesign, without
removing their grouping behavior. All perspectives use the established
**left filters / center table / right molecule inspector** layout and proportions.
View switching, searching, grouping, ordering, access/project/source filters,
and original-record selection are on the left. Collections' molecule, numeric
score and evaluation-context filters are also on the left. The established
molecular sheet and native structure/evaluation inspector are reused.

The first current item is inspected automatically; no “Load molecule inspector”
toggle is required. A row click updates the inspector and left-side selection.
The inspector loads only the selected original molecule, preserving model,
structure-integrity and evaluation behavior. Grouped empty/search states do not
show a stale off-filter record. Grouped sheets use focus-only navigation while
Collections retain their existing multi-selection and editing/export actions.

Exact-source primary regression suite: **221 passed**; focused isolated peer
suites: **59 passed on each compute PC**. New checks assert the three-column
placement and automatic inspector in Molecules, Runs and Collections, plus
focus/filter changes and backward-compatible sheet selection. Browser checks
use synthetic molecules and confirm the native 3D viewer and row-click updates.
Evidence: `docs/audits/2026-11-12-database-viewer-layout/`.

All six active application databases per host have fresh private,
integrity-checked cold backups under
`~/.local/state/mn-molecular-design-studio/releases/viewer-7a82f70af11f/`.
The coordinated Studio-only restart began with no queued/running/claimed Studio
jobs or pending cluster commands. Post-restart paths, process environments,
production authentication, API/UI health and fingerprints are verified.
All seven research table counts and all account rows, including the current
password, are exactly unchanged on each host. No authentication reset, research
cleanup, ownership change or scientific migration was performed. The shared
scheduler, teaching pilot and TLS gateway remain unchanged. Refresh the browser
and sign in again with the current password if prompted.

## Earlier grouped molecule/run views — 12 November 2026

The earlier application was **`6064409bf7da583b9f742438513de7573389ced5`**
on PC002975, PC004032 and PC004033. It superseded the navigation-only update
below. Source: `~/programs/deployments/mn-molecular-design-studio/prod-6064409bf7da/source`.
Archive SHA-256: `d89cf74aec858431ebbff61cbf227722e53d5a1f2969ca09082087eeeb349032`.
All three agents report application-code fingerprint
`29398c0d25485deb2c5062c8de127ca33323454e9f1cb2573f2fd2d36ee32d21`.

- **Molecules** is the new default perspective. Whole-molecule grouping uses
  recorded component contents with chain multiplicity retained. Individual
  sequence grouping links exact polymer type/sequence occurrences across jobs
  and complexes. Molecular-context variants retain their distinct original IDs.
  An overview group is not a new scientific identity or an averaged result.
- **Runs** combines collections from the same node/task/run/owner/project and
  optionally groups explicitly linked campaign children. Shared evaluation
  containers are resolved through per-row observation/source provenance, not
  their stale first-run group origin. Different workstations, task groups and
  ownership contexts remain distinct even when run labels match.
- Campaign parent links/IDs must be explicit. Cycles, conflicting lineage and
  missing identities are not guessed. In particular, old accepted-RFD4 singleton
  registrations with blank parent/campaign IDs cannot safely be joined merely
  from their names, targets or registration batch. Verified lineage backfill
  remains separate future work; no production lineage was rewritten here.
- **Collections** retains the existing import/edit/export workflow. Grouped
  items list contributing jobs and permit drill-down to original collections.
  The original molecule inspector is loaded only on explicit request. Scores,
  structures, targets, evaluator versions and settings are not combined by
  these views. No group, record or scientific identity was merged or deleted.
- ACL-filtered metadata grouping and search cover the whole eligible scope
  before 50-item pagination. Related history is drawn only from accessible
  groups. Browsing does not start ingestion or load coordinate files. The
  metadata scan is linear, not an indexed constant-time materialized view.
- A read-only primary Library probe finds 2,341 molecular compositions, 1,527
  individual polymer sequences, and 1,506 recorded jobs represented by 1,481
  campaign/run groups. One composition links 135 jobs; one sequence links 206;
  one explicitly linked campaign groups 24. These are view-scope counts, not
  physical database totals. The probe completes in 3.07 seconds on the host;
  this is one observation, not a latency guarantee for larger databases.
- Exact-source isolated primary regression checks: **215 passed**. Focused
  read-only/container checks: **53 passed on each compute workstation**. Cases
  cover aliases, complex partners/stoichiometry, stereochemistry, polymer types,
  identity variants, incomplete legacy data, shared evaluation containers,
  homonymous runs, explicit/cyclic/conflicting lineage, ACLs and late-page
  search/navigation. Browser checks use labelled synthetic data.
- All three active application SQLite/configuration sets have fresh private,
  integrity-checked cold snapshots under
  `~/.local/state/mn-molecular-design-studio/releases/grouped-6064409bf7da/`.
  Jobs and commands were idle before the coordinated Studio-only restart.
  After restart, process paths/environments, production authentication, API/UI
  health and all agent code fingerprints are verified. All seven research table
  counts and all account rows (including the current password) are exactly
  unchanged per host. The scheduler, teaching pilot and TLS gateway are unchanged.

Refresh the research tab and sign in again if prompted, using the **current**
password, not an expired temporary password. Hot reload remains disabled.
Evidence is in `docs/audits/2026-11-12-grouped-database-views/`. This UI release
does not repair the previously recorded peer-sync discrepancies or extend the
scientific-tool qualification of the original acceptance.

## Earlier library navigation update — 12 November 2026

The earlier navigation-only application was **`304fc295416995f0e338c745582c2772053940fd`**
on PC002975, PC004032 and PC004033. It superseded the original application
source below and is now superseded by the grouped views above.

- Archive SHA-256: `10b747f4b50995112b238fb1226566f923799a3bef52ece73acbf9630898799a`.
- Agent code fingerprint on all three workstations:
  `5e040dcd4cf4d8b3097e35fb7d0848b4f83c896e97c749938894c177d0c662e8`.
- Source: `~/programs/deployments/mn-molecular-design-studio/prod-304fc2954169/source`.
- Collections now have prominent search, provenance-based categories, deterministic
  libraries-first/name/date/size sorting and 15-collection pages. Imports are
  easier to find among generated single-result collections. Access/source filters
  are collapsed, the selected collection has a clear heading and source badges,
  and an empty project filter explicitly means **All projects**.
- Categories/search narrow collection navigation only; **All molecules** retains
  the existing access-filtered scope. No groups are merged, renamed or hidden
  based on names or size. Direct Results-to-Database navigation resets browsing
  filters and opens the correct page; changing pages retains the current selection.
- Final isolated primary database/review suite: **200 passed**. The focused
  **38-test** suite also passes independently on each compute workstation in
  network-disabled containers without production storage, GPU or Docker access.
  Browser checks use explicitly labelled synthetic preview data. These are UI
  regression checks, not new qualification of scientific tools or pipelines.
- Studio jobs and cluster commands were idle before the coordinated restart.
  Six active application SQLite databases per host have fresh integrity-checked,
  private cold snapshots with systemd/configuration recovery copies under
  `~/.local/state/mn-molecular-design-studio/releases/library-304fc2954169/`.
  This supplements the earlier full recovery snapshots; it does not recopy
  scientific files or establish an off-site backup service.
- Runtime working directories, process environments, production authentication,
  API health, primary Streamlit health and all three agent fingerprints were
  checked after restart. Accounts, passwords and node IDs are unchanged. All
  seven checked research table counts are exactly unchanged on each host.
  Primary retains 21,473 molecules and 1,050 groups; each peer retains 21,470
  molecules and 1,050 groups. The previously reported peer-sync warnings/count
  differences are **not resolved by this UI-only update**.
- The shared scheduler, teaching pilot, TLS gateway and authentication policy
  are unchanged. Live hot reload remains disabled: refresh the research browser
  after deployment and sign in again if prompted.

Focused test XML and non-secret runtime receipts are retained in
`docs/audits/2026-11-12-library-navigation/`. The earlier evidence archive below
continues to describe the original production acceptance, not this follow-up.

## Original production source (superseded above)

- Application commit: `87cb209e3884844076434df601e80a6152bf3c18`.
- Archive SHA-256: `db8a705786b7e8bab5f41cfe79a1e51b1c66b203c1c7d8a29d04e73bf566353f`.
- Runtime application-code fingerprint, all three agents:
  `62bf990569d7d8ab64bab9ade6ddd194d02b0214192826e977c0f7831db3d487`.
- Source: `~/programs/deployments/mn-molecular-design-studio/prod-87cb209e3884/source`,
  under `user` on PC002975 and **`compute`** on PC004032/PC004033.
- Every active API, agent and the primary UI loads that source with production
  authentication. Fresh PIDs, executable, `/proc` environment, working directory
  and agent-reported fingerprint were checked, not just Git HEAD/API health.
- Peer UIs remain inactive as before, with updated pinned launch definitions.
  Primary hot-reload/file watching is disabled. Node UUIDs, peer credentials,
  cursors and scientific paths are retained. Compute development checkouts remain
  on their original `main`/`d5df4d66`; bundles import objects without resetting them.

The research gateway, renewal helper and corresponding tests are separate
operational artifacts delivered with this report. The original acceptance used
the pinned commit in this section; the navigation follow-up is identified above.

## Corrections closed

### AUD-008: whole-scope review

The new query scans all authorized scoped metadata in batches of at most 250
models inside one SQLite read snapshot. Saved-system membership, acceptance,
search and method/version/settings/image-qualified numeric predicates apply
**before** selecting a page. Counts, facets and measurement ranges no longer
come from the first page. Filter changes reset pagination; an empty page does
not discard off-page selected IDs. Actions still reauthorize exact selections.

Evidence linking is shared with exact-model review. No coordinate hydration,
name-only joins, failed-result numbers or arbitrary conflict resolution was
added. Conflicting repeated measurements remain missing; measured zero stays
zero. Different/unknown image identities no longer collapse into one series.
Tests cover 262 models, late-page system/search/acceptance/numeric matches,
version changes, conflicts, ACL denial, concurrent mutation and Streamlit journeys.

Cost is linear in scoped metadata. This is not an indexed predicate cache;
unusually large per-molecule evidence histories can still increase hydration
cost. Bounded batching must not be mistaken for a constant-time whole-campaign query.

### Safe production sign-in

Local `create-admin --preserve-development-owner` binds only the first admin to
the existing development owner and requires a password change. Default bootstrap
still creates a fresh ID. Existing-account and competing-bootstrap attempts
cannot overwrite accounts. No HTTP ownership-transfer endpoint was added.

The instructor was created once on the primary. Its one authoritative user event
was received through the normal replica receiver on each peer. No account or
replica database was substituted between hosts. Real sign-in works on all three
APIs; anonymous access and temporary-password workflow access are denied.
Acceptance sessions were revoked. The `legacy-unassigned` group is not silently
reassigned and no existing research was shared with students.

## Validation results

| Completed check | Result |
|---|---|
| Unique full CPU coverage, disjoint shards | **3,172 passed, 37 skipped** |
| Core final shard | 3,119 passed, 37 skipped, 103 warnings |
| Workflow-mode / bootstrap / TLS shards | 44 / 6 / 3 passed |
| Final focused review/catalog/auth/API/concurrency slice | 107 passed |
| Same focused slice on each compute host | 107 passed on each |
| Coordinated recovery and isolated application rehearsal | All three passed |
| Seven canonical scientific tables per host, before/after ordinary peer resumption | Digests unchanged |
| Real production sign-in over three loopback/SSH paths | Passed; acceptance sessions revoked |
| Browser UI | Real sign-in form; no development auto-login |
| TLS/network | CA-verified health 200; API 404; outside/forged source 403; HTTP 400 |

Final shards reconcile exactly to **3,207 collected case IDs plus two optional
Chaperone module skips**, with no duplicate cases counted. Optional native/GPU/
container skips are not passes. Earlier attempts exposed incorrect assumptions
in the new bootstrap test fixture about project metadata/API naming; the fixture
was corrected and affected final shards rerun. Full coverage is not claimed as
one successful monolithic pytest invocation.

Native scientific evidence remains in `PRODUCTION_TOOL_FIXES.md`: 12 unique GPU
cases, both ESM methods, seven CPU protein QC methods, native Genie3 on the primary
and a real 5090 generation, plus a CPU graph and authenticated self-contained ZIP.
Those runs used the earlier corrected scientific implementation, principally
`86804ff`, **not fresh inference of every engine at this commit**. Current changes
affect review/bootstrap, not those inference adapters. Prior ligand-tool results
remain historical evidence, not newly rerun results.

The missing verified CPU QC image was delivered to PC004032. All three hosts now
have `mn-protein-qc:1.0.0` executable-content fingerprint
`sha256:2a0cee2e4df83591fc2e2088167b16c500179b2f8609ae58930c024179af9088`.
The separately verified SM120-compatible Genie3 tag is preserved. No old image
tag, driver, model weight or shared runtime cache was blindly overwritten.

## Recovery and preservation

Studio writers and peer links were quiesced on all hosts after confirming zero
Studio queued/running/claimed scheduler jobs. Shared schedulers, unrelated
applications and the teaching pilot were not stopped. Private recovery sets:

- Primary: `/mnt/data/RESULTS/.mn-mds-production-recovery/20261112T093000Z`.
- PC004032: `/data/RESULTS/.mn-mds-production-recovery/20261112T093000Z`.
- PC004033: `/data/results/.mn-mds-production-recovery/20261112T093000Z`.

They include complete application/workdir trees, sealed inputs, settings, per-host
cluster state, consistent read-only-source SQLite backups, scheduler state and
systemd configuration. **42 SQLite files**, including retained older copies,
pass integrity checks; foreign-key results match their sources. Shared-scheduler
snapshots are online while unrelated work continues; the Studio queue is empty
and no pending Studio cluster command is introduced.

All non-SQLite entry sizes/mtimes/types and symlink targets are compared. Unchanged
files are hard-linked from that host's prior isolated backup, **never live files**;
changed files are copied. Do not edit either backup set. Previously approved
exact unreadable `campaign_settings.json` exceptions remain recorded, with
readable source provenance retained. No `.nv` traversal/permission change or
broad cache exclusion was added.

Separate writable application clones run the new app with networking disabled and
recovered workdirs read-only. Seven canonical scientific-table digests, node IDs,
three protected campaign groups and 60 retained campaign models per host remain
intact. Fifty structure references are checked per host: local byte/hash checks
cover 50/48/44 files respectively; the remaining 2/6 peer entries are explicit
`cluster-library://` remote references, not local file paths. They are not claimed
as failed canonical coordinates or as successful authenticated remote downloads.
Legacy-owner bootstrap and denied student access pass without changing science.

Initial peer rehearsals omitted an empty reference-directory mount and failed;
only private copies were affected, and corrected rehearsals pass. Rehearsal
accounts are disposable: **never restore mutated `isolated-restore` clones into
production**. Use authoritative snapshots and the reviewed recovery procedure.

This is not a new cryptographic reread of every historical trajectory byte, an
off-site disaster-recovery copy or a clean-machine/model-corpus reconstruction.
No campaign/import/test group was deleted or old job relabelled. Historical
queued/running/collecting catalog rows remain; no live lease/pending command makes
them active, and scheduler absence alone is not deletion authority.

### Rollback boundary

Each host also has a private `postcutover-identity-config` supplement containing
fresh auth/replica SQLite backups and current service configuration. The instructor
account payload matches its replica event in these paired snapshots. The primary
also retains its private bootstrap credentials and research TLS configuration/keys
there. These six additional online database backups are distinct from the 42-file
quiet-point capture; no secret supplement is included in public evidence. Future
research/password changes require newer backups; this is not an ongoing backup service.

Do not revert authentication to development, delete admission-v2 state, replace
another host's replica catalog, reset epochs/cursors or mix old blind review
writers. Stop/drain all writers and peer delivery, retain current paired
scientific/auth/replica state and trial recovery privately before rollback.
The staged `tools-86804ff29f50/source` is a nearer compatible investigation fallback
with production authentication retained. Former `p11-25e003c` is retained for
forensics, **not approved for unconditional code-only rollback after new v2 work**.
Do not substitute a rehearsal clone or pre-account backup for current state.

## Remaining operator actions and limitations

1. Choose your password, open important campaigns/imports and complete an intended
   real research workflow before relying on new experimental output.
2. IT/operator: trust the research CA on actual clients, confirm their subnet and
   check browser/WebSocket access. Local tests do not test students' laptops.
3. Confirm student admission/roster and compute/storage fairness. The isolated
   pilot remains the immediate student route; production access is not silently granted.
4. Many predictors, generator modes and MD/free-energy/QM/ADMET workflows lack
   fresh native end-to-end qualification. The 273-row historical matrix is
   overlapping inventory, not 273 independently tested tools. Computational
   affinity/confidence is not experimental binding or efficacy evidence.
5. Unknown historical links are not repaired. Atomic admission covers exact
   Database assessments, not every generic campaign/pipeline. Off-site backup,
   alerting and institutional resource policy still need an operator.

Evidence: `PRODUCTION_RELEASE_EVIDENCE.tar.gz`. It excludes account databases,
passwords, tokens, private identity-event payloads, private keys and scientific
database contents. Original independent audit evidence/worktrees are unchanged.

Evidence archive: 67 allowlisted files, 175,276 bytes; SHA-256
`783cb4fafc1415de037a17448fe1f88c681ea3305208b0ce1632bace23f29e80`.
The allowlist was checked against the actual private bootstrap password and
replicated password hash before publication; neither appears in the archive.
