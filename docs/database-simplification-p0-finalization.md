# P0 finalization receipt

Date: 2026-11-05. **Code completion: GO. Production deployment: NO-GO pending operational acceptance and explicit approval.**

The accepted independent [closure review](database-simplification-p0-closure-review.md) remains the baseline. No additional defect was found, no implementation correction or delegation was needed, and the approved code/test bytes were preserved. No reconstruction, scientific reimport, migration, deployment, service restart or workstation synchronization was performed.

## Git boundary

- Before finalization: branch `main`, HEAD `1d4c8fb7b74b0494700a60d200dc452a9c2b3417`, upstream `origin/main`; 16 local commits ahead.
- Final release identity: the full SHA of the single commit containing this receipt, reported separately after committing. Use that full SHA for all deployment receipts. Do not substitute the pre-finalization HEAD.
- Included paths: `mn_molecular_design_studio/engine/workflow_engine/workflows/input_batch.py`, `tests/test_api_integrations.py`, `docs/database-simplification-p0-closure-review.md`, this receipt, and `docs/database-simplification-p0-deployment.md`.
- Only the existing approved 16-line dispatch-lineage diff and 201 added regression-test lines change implementation/tests. The closure review is included unchanged. The other documents record final verification and deployment requirements.
- SHA-256 manifests verified all 76 initially untracked files unchanged and both approved implementation/test files unchanged before staging. The closure review is intentionally committed; the other 75 untracked audit files (three audit directories and four Phase 1 reports) remain untracked. No unrelated tracked changes were present.

Remote: `ssh://git@ssh.github.com:443/martenbeeg-sketch/mn-molecular-design-studio.git`. A read-only `git ls-remote` confirmed `refs/heads/main` at `d5df4d6616b624c69776645020d40307588ce6e3`, matching the local remote-tracking ref. The existing 16 unpublished commits include unrelated analysis/Phase 1 work. Publishing `main` would include them, so **no push is authorized under this task's push policy**. No branch switch, new branch, amend or history rewrite is needed.

No CI workflow files were found at the conventional GitHub/GitLab/Jenkins/Azure/CircleCI paths in either HEAD or the known remote tree, and local hooks contain only samples. Server-side branch protection, repository policies, external CI/webhooks and automatic deployment effects could not be verified (`gh` is unavailable and no connected GitHub policy API is configured). Absence of local workflow files does not establish that a push is safe. Leave this commit local; do not publish unrelated history to distribute the release.

## Final verification

Interpreter: `/home/user/mambaforge/envs/mn-molecular-design-studio/bin/python`. Suites ran once each, with `MDS_GPU_TESTS=0`, `MDS_QC_DOCKER_TESTS=0`, `-q -p no:cacheprovider --tb=short`, unique temporary state directories and JUnit output. No broad forensic audit or historical inventory was repeated.

| Selection | Final result |
|---|---:|
| Same 25-module combined selection documented in [F1–F5 corrections](database-simplification-p0-f1-f5-corrections.md#full-regression-and-state-boundary) | **607 passed in 152.39 s** |
| `tests/test_api_integrations.py -k 'vina_batch_worker_before_api_marker or vina_batch_parameters_cannot_forge'` | **21 passed, 24 deselected in 22.24 s** |
| Documented independent scenarios, legacy completion, previous final review and closure probes | **97 passed in 40.38 s** |
| `git diff --check` | Pass |

Independent modules: `/tmp/mds-p0-final-independent-review.qWvGDQ/{test_original_scenarios,test_legacy_completion,test_final_review}.py` and `/tmp/mds-p0-closure-review.5cZ0XR/test_closure.py`. Evidence for this final run: `/tmp/mds-p0-finalization.eyri92zn/{combined,focused,independent}.{log,xml}`; preservation manifest: `/tmp/mds-p0-finalization.eyri92zn/initial-state.json`. Temporary evidence paths are local and are not a permanent archive.

Verified boundaries:

- **F1:** omitted HTTP/MCP import selections create no scientific imports; explicit exact IDs and positive top-N remain available; invalid/conflicting selections are rejected and retries stay idempotent.
- **F2:** saved positive graph limits remain preserved but ordinary campaign completion publishes zero imports.
- **F3:** both exact owner/run SpiroChem shortlists remain excluded from Library on source and real-adapter replicas. Direct history, independent curated memberships, private ACLs and tombstones remain effective; development filtering stays restricted to the eight audited owner/run pairs.
- **F4:** Evaluation deduplication still requires exact historical observation/source/group/owner/run/compound identity and authoritative replica markers. Unrelated or private foreign observations cannot suppress a selected local Evaluation.
- **F5 including dispatch race:** authenticated `studio_owner_id` and string `studio_project_id` metadata tags authorize child lineage before the parent Studio marker exists. Generic ownership, nested/item parameters and `input.json` cannot authorize it. Real worker timing cases cover coalesced and separate/mixed docking choices, omitted and explicit policies, exact scores/source/settings/owner/project attribution, pose/CSV retention, canonical-ID reuse, parent stability and serialized-worker retry. Non-coalesced parents remain manual to avoid aggregate duplication.
- **Preservation/navigation:** intentional HY-L126 fixture import, binder memberships, historical Results/references, structures, accepted designs and owner-qualified replication regressions pass. Ordinary Database/Results navigation does not initiate scientific ingestion in the exercised paths; existing derived projection caches may update.

The reviewed `input_batch.py` SHA-256 is `dedbd6e6b376564110b8429fd78d64a36029e902bfbe44f33688552468d5368a`; `tests/test_api_integrations.py` is `8220c1c0c333b3512b786046d817889586986fe6b612c65c4eaef22fc45835eb`. Expected package source fingerprint using the existing `capabilities.code_revision()` algorithm: `a8c872bf2e2de4b541e20b111aef0304b1c6e5d1429a46af5ff2f3b8f39ac0db`. This is source evidence, not proof of any deployed process's loaded code.

## Deployment and limits

Follow the prepared [three-workstation deployment and rollback procedure](database-simplification-p0-deployment.md). Remote login on PC004032 and PC004033 is `compute`, as confirmed by the user; this does not change scientific owners. No remote connection was made for this finalization.

Before any approved rollout: drain Studio jobs across all three hosts, gate writers/peer traffic, verify pinned SHA, take independent WAL-consistent backups of databases and scientific files/configuration/replication state, prove isolated restoration and record old runtime environments. Install the same reviewed source and verify the code actually loaded by APIs and workers before reopening cross-host work. Use isolated canaries for every check that might write; inspect R3, other IDR/Tau, PD-L1, HY-L126, selected binders/Evaluations, shortlist visibility, ACLs, tombstones and ordinary navigation against a small existing production baseline.

Prefer code/environment rollback while retaining current scientific data. A Database restore requires explicit accounting for every legitimate post-backup record/file and coherent ownership/replica/tombstone state; never silently replace new data with an older snapshot.

Limitations remain those accepted in the independent review: temporary fixture databases and Docker output substitutes establish ingestion/selection behavior, not live computation, transport, per-workstation version parity, or a new production HY-L126/campaign inventory. No live peer synchronization or production preservation checksum verification was performed. **P0 is ready to close at code level; production deployment remains blocked on the operational evidence and explicit approval above.**
