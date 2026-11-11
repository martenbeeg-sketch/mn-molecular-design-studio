# Saved design review workbench

The local BinderScope viewer's board/query pages were inspected as a usability
reference: compact design rows, evidence plots and a focused decision panel.
Studio now uses that review pattern inside its existing Database workspace;
BinderScope code, calculations, identities and storage are not duplicated.

Opening a Protein design cohort now shows one review workbench rather than a
source-membership table followed by a molecule/target/study selector. Each row
is an exact retained Database representation. Short display labels disambiguate
repeated worker names with a source number; original IDs stay in details.

- Review up to 250 models per metadata page. Final computational passes are the
  default, with explicit earlier-stage and diagnostic cohort choices.
- Saved method/version/parameter-specific values appear as selectable columns
  and an interactive scatterplot. Row and point/box selections use one stable
  structure-ID selection. Checked design labels and plot symbols identify it.
- Missing, conflicting and failed measurements are not plotted as zeros or
  chosen arbitrarily. Different parameter sets remain separate series. Units
  and original definitions/provenance are preserved in focused details.
- Native confidence, scientific assessments, computational acceptance and human
  shortlisting remain distinct. No new ranking or scientific calculation runs.
- The focused design exposes saved sequences, chain roles, acceptance policies,
  all saved measurements/status, exact identities and existing Inspector,
  assessment, 3D and persistent shortlist/export actions.
- Coordinates load only for an explicit inspection/assessment action. Review
  and shortlist controls reuse existing Database/Evaluation/review contracts.
- Campaign totals/input provenance are available on explicit request; MD and
  other work retain their recorded system/study navigation. No schema changes,
  Results scans, remote discovery, scientific edits or deployment occurred.

## Evidence and limitations

The existing isolated 8510 Database snapshot's RFD4 10-design historical batch
group contains 174 source models. The read-only model/evidence projection took
0.181 s. All have saved evidence; four have explicitly linked final passes on
these exact representations. Earlier native records and downstream structures
are not assigned a final pass merely because a candidate name matches. Analyses
on a different coordinate digest remain attached to their actual representation.

Browser verification opened this real cohort, rendered its four final-passed
models with saved generator ipTM/contact values, selected a table row and used
plot box selection to select Design 2 (ipTM 0.82, 31 binder-contact residues).
The checked row, focused evidence and exact downstream actions synchronized.
No structure coordinates, new assessments or scientific records were required.
Shortlist/export integration and mixed RFD4/BindCraft2 rows are covered by small
Database fixtures; existing scientific calculation/export evidence is reused.

The metadata projection reads only existing authorized memberships, structures,
observations, saved Evaluation descriptors and user review annotations. It
requires an exact representation link or matching recorded digest with compatible
candidate/model identity. Source run, molecule or name alone cannot join evidence.
Multiple representations of a molecule remain separate; duplicate measurement
conflicts remain blank and all source values remain inspectable. This strict
linking can expose historical gaps instead of borrowing another model's evidence.

Focused regressions cover ACLs (including private same-molecule evidence), exact
digest aliases, model-specific acceptance, failed/conflicting/missing/zero values,
parameter separation, paging, readable labels, shared selection, no-coordinate
review, existing shortlist/export controls and surrounding navigation contracts.

Validation: 103 focused tests passed across the design workbench, system and
study navigation, saved evidence, collection labels/state, authorized catalog,
Database hydration and existing shortlist exports. No scientific calculations
were repeated. Production and unrelated cleanup files remain unchanged.
