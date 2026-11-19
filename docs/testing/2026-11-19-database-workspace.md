# One Database workspace

## Change and compatibility

The single Database navigation entry now renders the scientific folder browser.
The old Data Browser session route is an alias of that same destination. The
legacy screen is no longer dispatched by the application; its shared import,
collection, export, score-rendering and tool-input helpers remain available.

Library actions exposes imports, empty collection creation, collection inspection,
metadata columns and confirmed membership removal. Collection removal states its
actual scope: group-owned coordinates, observations and provenance are removed;
other collections and original Results folders are not. Whole-run/project cleanup
remains on Projects. Confirmation is bound to the exact selected member IDs.

Data actions operates only on explicitly selected table rows. Existing molecule
exports, comparison and tool-input handoffs remain available there. The scientific
table's own CSV export still preserves exact displayed evaluation context across
all matching pages. Molecule exports retain their existing broader accessible-
molecule semantics and say so explicitly.

Project links, Results links and import/group completion requests now select the
new tree. Requests clear stale searches, refinements and batch selection, preserve
the authorized exact source collection and focus the requested molecule, including
on a later table page. Focusing does not select designs for analysis. Explicit
source collections can be inspected without reintroducing every automatic batch
as a normal library folder. Inaccessible destinations produce a warning rather
than silently opening an unrelated collection.

Curated data remains the default. Historical indexed SQL is an explicit advanced
option; it does not register jobs, admit designs or imply that anything passed.
Read-only previews cannot import, create collections or open write dialogs.
Management writes use the same database root as the browser, never a separate
default database. Descriptive Properties editing and BinderScope submission keep
their existing access checks and scientific contracts.

## Verification scope

Streamlit AppTest fixtures exercise the real controls, including import review and
commit, empty collection creation, metadata columns, exact source navigation,
project navigation, row selection, selection-bound removal, molecule download and
tool reuse. SQL fixtures check that removing one membership preserves another.
Navigation tests cover the old session alias, one sidebar entry, both application
routes and no automatic analysis selection. Additional browser/layout, identity,
metadata, cache, scientific tree, Results and comparison regressions use isolated
homes, databases and cluster directories.

No research records, scientific scores, admission decisions, passwords or grants
need migration for this change. Tests do not submit research analyses. The real
in-app browser was signed out during verification; authenticated production mouse-
click verification is not claimed.

The focused and adjacent regression set passes 340 tests. A wider run reports
five failures in `test_entity_name_replication.py`: replay calls
`workspace.save_organization` while a molecular write transaction is active, and
the separate naming projection connection times out with `database is locked`.
The first identical failure also reproduces on the unchanged previous release
`026547e` in a fresh isolated home. This existing synchronization defect is not
corrected or hidden by the page consolidation; this change does not modify the
replication, naming projection or transaction implementations.
