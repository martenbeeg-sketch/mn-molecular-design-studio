# Scientific data organization v2 — Phase 1

Implemented and verified locally on 2026-11-09, starting from `2c30eb1`.
This is an additive first slice, not a production migration. No production
Database, Results, replication epoch, service or scientific file was changed.

## Implemented contract

`studio.scientific-lineage.v1` is a rebuildable relational projection inside the
existing Molecular Database SQLite file. It adds scientific jobs, output
occurrences, typed dependencies, campaign summaries and campaign membership.
It does not replace molecules, structures, observations, Evaluations or reviews.

- Jobs are qualified by stable source node, task group and run ID. Explicit
  parent metadata and distributed-ledger cluster-job receipts establish execution
  relationships. Names and filenames do not establish parentage. Project and
  owner boundaries prevent cross-scope campaign grouping.
- Output occurrences retain candidate/native IDs, model selectors, CandidateSet
  target/backbone/variant/prediction lineage, sequences, chain roles, raw source
  metadata, original measurements, existing normalized measurements and recorded
  acceptance. Equal coordinate hashes do not merge model occurrences. Model `0`
  is preserved; an unknown model is not assigned to a known model.
- Artifact metadata retains its original path, portable source reference, owner,
  digest and verification evidence. Absolute container paths are provenance,
  not permission to read outside a job. `owner_present` describes observed source
  availability; `owner_local_unverified` does not promise a local/offline copy.
- Retained references join existing structures by source, owner/project,
  candidate, exact recorded model and available digest. Unmapped historical
  representations remain explicit references with warnings. No molecule or
  structure is created by projection/backfill.
- Evaluation outputs retain the existing descriptor ID and observation ID,
  method-specific measurements, parameters, units and source context. Typed
  input edges represent exact representations and multiple inputs, including
  PAE/reference files or recorded MD trajectory snapshots. Missing historical
  input identity remains unresolved.
- MD records distinguish study, independent replica, continuation segment and
  trajectory-analysis input. GB and PB remain separate existing Evaluations.
  Docking records preserve compound/preparation/pose context; a compound is not
  equated with a pose.

The existing completion hook captures lineage independently of the existing
Molecular Database admission policy, including terminal failed/stopped jobs.
It neither promotes rejected candidates nor changes scientific completion.
Capture failures are logged and retriable by explicit backfill.

## Synchronization and locality

An optional `scientific_lineage` extension travels in the existing
owner-authoritative job replication record. Updated receivers validate its
author, job envelope and ACL scope before materializing the projection. Older
extension-free job records remain supported; older receivers retain the optional
JSON without needing new record kinds or epochs. Two isolated peers were tested
through the existing capture/receive path, replay and reopen.

Metadata synchronization transfers no coordinate files or trajectories. Retained
structures still use existing authenticated Database/artifact access. The
`analysis_placement` helper describes the owner of required Results inputs and
rejects unavailable or multi-owner inputs without explicit staging. This phase
does **not** introduce or change scheduler placement, dispatch or transport.

## Campaign UI and API

Complex Analysis **Select** now opens the indexed campaign browser:

1. Select project, scientific activity and recorded exact construct.
2. Open a campaign to see its distributed jobs and original IDs.
3. Inspect native producer counters with their original units/stages/scopes,
   separate indexed occurrence counts and retained representation counts.
4. Keep final computational pass as the default, or explicitly choose an earlier
   recorded passing stage. Native passing never becomes final acceptance.
5. Open exact retained models directly in the existing Analyze → Compare →
   Shortlist workspace, without another Collect operation.

Primary Results outputs define the default cohort. Controller/worker references
remain visible under an explicit diagnostic option, avoiding repeated counting
of a parent aggregate and its underlying executions. Ambiguous multiple output
sections do not silently become one final cohort. Saved cohorts and individual
Database selection remain available as an advanced path. Existing Results
inspection remains intact.

Authorized, paginated API:

- `GET /api/v1/molecules/scientific/campaigns`
- `GET /api/v1/molecules/scientific/campaigns/{campaign_ref}`
- `GET /api/v1/molecules/scientific/campaigns/{campaign_ref}/models`

Reads use indexed tables only: no historical directory discovery, remote SSH or
completion capture on page render. Full records/coordinates hydrate only when
models are explicitly opened. Source engine display uses an exact linked
scientific candidate, rather than calling a pipeline/controller a generator.

## Controlled historical validation

Read-only source connections and explicit selected files populated an isolated
subset at `/tmp/mds-scientific-v2-u51u83xz`. Networking, peer synchronization,
admissions and dispatch are disabled. Original canonical IDs were copied into
staging; none were reimported or regenerated. This is not a completeness audit.

| Historical evidence | Result and mapping boundary |
| --- | --- |
| BindCraft2 parent `campaign/20261029-094447-d4f389fb`, controller `design/20261029-094457-fb2fd80a`, PC004033 worker `design/20261029-094502-6ad69fc5` | One three-job campaign through metadata plus ledger receipt; one authoritative native-passing parent output. Supplemental execution/representation references remain separate. |
| RFD4 `design/20261028-121331-91d15b3f` and `design/20261028-121326-52a55d64` | Separate occurrences/campaigns with identical recorded 115-residue target sequence. Existing independent representation references retained; unknown prediction/model associations stay unknown. Cross-engine target equivalence is not asserted where only different target identifiers/digests are available. |
| PC004033 `design/20261102-081701-97f7d42e` | Job metadata preserved; no candidate/structure rows available in this bounded source replica. Explicit unavailable-source warning, no replacement model. |
| `refolding-validation/20261102-184847-e1c4c3e4` | Five authentic LocalFold samples and five existing Evaluations; absent explicit model IDs remain unknown. Supplemental representation references are not fabricated samples. |
| PC004033 `docking/mnligand-3f28268e02567ac4` | Three saved TRα / 2H79 / Medchem / UDP pose/score representation references and three existing Evaluations, including original chemistry/preparation context. Raw full docking output was not fetched. |
| MD `ligand-md/mnligand-md-99d094e268164eca`, endpoint `ligand-mmgbsa/mnligand-md-469ca9ed9a8e156c`, `md-analysis/mnligand-md-85995d463974aade` | One saved study `study-98bebb39c3ebe5f3`: three replicas, six segments, six production Evaluations, six GB plus six PB Evaluations and six trajectory analyses. Metadata reconstructed solely from saved observation contexts; exact recorded trajectory references/hashes preserved. Execution jobs without explicit parent evidence remain separate. |

Real browser, isolated port 8510:

- Opened Complex Analysis → BindCraft2 campaign: parent/controller/worker shown.
- Final cohort correctly empty; explicitly selected Native generator acceptance.
- One authoritative output/one exact retained structure offered.
- Opened Analyze: generator `bindcraft2`, one matching saved Evaluation:
  `evaluation-7a3bc3b55327dccbc1651776de6bec56`.
- Opened Scientific Inspector with molecule
  `mol-861abd93038523a290940671b4b62662`, structure
  `structure-d2772c6c80ed2c92a942d89081ab57ff`, source parent run and digest
  `a6ebc954057e53aa4907f277c13bec13e22fcf324a5272762a03b8669f3fea9b`.
  Returning preserved the selected model and native-stage cohort. No analysis ran.

Local evidence: `validation-report.json`, `source-metadata/`, `preview.json`,
`campaign-browser.jpg`, `analyze-browser.jpg`, `inspector-browser.jpg` in the staging
directory. These local validation files are not part of a production backfill.

The bounded SSH retry connected as `compute@PC004032`. PC004033 returned **No route
to host** before authentication; no live PC004033 artifact retrieval is claimed.

## Chaperone compatibility

The small versioned schema translator `studio.chaperone-lineage.v1` was compared
with the local native schema at Chaperone
`2d3ac1e7ce8078e7b99ce1511940fda61c7198ae` and the
[upstream syntax](https://rosettacommons.github.io/chaperone/syntax/).
It preserves design/parent IDs, sequence-only records, JSON sequences/metrics,
ordered structural history, step metadata and CatCIF/Parquet references. Ordered
history is not relabeled as independent predictions. A unique native parent in
the supplied batch resolves; otherwise parentage stays unresolved. This is a
fixture-tested schema translator, not a new Chaperone execution integration or
archive reader. The existing Chaperone adapter remains unchanged and optional;
no Python 3.12 dependency was introduced into Studio.

## Tests, responsiveness and limits

Focused command:

```sh
python -m pytest -q tests/test_scientific_lineage.py tests/test_logical_runs.py \
  tests/test_complex_analysis_design_collection_labels.py \
  tests/test_complex_analysis_database_hydration.py \
  tests/test_complex_analysis_collection.py tests/test_cluster.py \
  -k 'scientific_lineage or logical or design_collection or database_candidate or hydration or replica_import_promotes_canonical_job or heartbeat_promotes_legacy_fallback_job' \
  --disable-warnings --maxfail=1
```

**43 passed, 100 deselected**, one existing deprecation warning. Coverage includes
ACLs, source spoofing, duplicate coordinates, model zero/unknown, target
incompatibility, distributed grouping/counts, existing Evaluation IDs, MD inputs,
docking chemistry, two-peer replay, Results-only completion without promotion,
staging backfill idempotence and an index-only Streamlit test.
Three additional focused source-library/admission regressions passed, including
queued submission registration through the existing completion hook.

Ten metadata-query samples in the small staging subset: campaign page median
**1.07 ms**, selected detail median **1.61 ms**. Browser observations: campaign
heading visible in **674 ms** after navigation; selected campaign controls in
**347 ms**. These are staging observations, not a full production benchmark or a
claim of equivalent performance on a fully backfilled cluster.

Current bounds:

- Snapshot publication explicitly rejects metadata above 1.5 MB; it does not
  truncate or claim a complete index. Large docking/ensemble jobs need paged
  metadata publication before comprehensive historical backfill.
- Campaign summaries are computed at write time from compact indexed job facts;
  regrouping currently visits all indexed job metadata, not full Results. A
  future bounded improvement is affected-campaign-only regrouping.
- Historical target/model/proxy/trajectory associations absent from recorded
  evidence remain unresolved. Target names alone never establish equivalence.
- Results-only locality is tested as a routing contract; no new owner-local
  calculation or live remote transfer was run in this phase.
- Historical backfill CLI is staging-only and explicit-manifest based. A marked
  existing staging destination is required. Production adoption requires a
  separately reviewed bounded backfill; no automatic migration is enabled.

Next bounded step: add paged owner-published output detail and validate one full
distributed historical campaign manifest, including explicit execution proxies,
before enabling a controlled production backfill.
