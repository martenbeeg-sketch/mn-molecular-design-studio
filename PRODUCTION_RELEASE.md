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
- Login: **`research_instructor`**. The temporary password is in the local,
  owner-only `~/.config/mn-molecular-design-studio/production-credentials.json`.
  Do not send that file to students or commit it. Change the password at first
  sign-in; normal API/workflow access is blocked until this happens.
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

## Exact deployed source

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
operational artifacts delivered with this report. Scientific application code
remains the pinned commit above.

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
