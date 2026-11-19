"""Design-level navigation without combining scientific evaluation records."""
from __future__ import annotations

from mn_molecular_design_studio.database_browser import _component_identity, _key


def _design_identity(molecule, group_id=None):
    components = molecule["components"]
    roles = molecule.get("component_roles_by_group", {}).get(group_id, molecule.get("component_roles", {}))
    if not any(values & {"binder", "designed_binder", "design", "target", "receptor"} for values in roles.values()):
        roles = molecule.get("component_roles", {})
        if any(values & {"binder", "designed_binder", "design"} and values & {"target", "receptor"} for values in roles.values()):
            return None
    designed = [component for component in components
                if roles.get(_component_identity(component), set())
                & {"binder", "designed_binder", "design"}]
    if not designed:
        designed = [component for component in components if component["type"] == "small_molecule"]
    identities = [_component_identity(component) for component in designed or components]
    if not identities or any(identity is None for identity in identities):
        return None
    return ("design" if designed else "molecule", tuple(sorted(identities)))


def design_subject(index, result):
    molecule = index["molecules"].get(result["molecule_id"])
    identity = _design_identity(molecule, result.get("group_id")) if molecule else None
    if identity is None:
        return None
    if "design_cohorts" not in index:
        cohorts = {}
        for molecule_id, entry in index["molecules"].items():
            for group_id in entry["group_ids"]:
                contents = _design_identity(entry, group_id)
                if contents is None:
                    continue
                cohort = cohorts.setdefault(contents, {"molecule_ids": set(), "group_ids": set()})
                cohort["molecule_ids"].add(molecule_id)
                cohort["group_ids"].add(group_id)
        index["design_cohorts"] = cohorts
    cohort = index["design_cohorts"].get(identity)
    if not cohort:
        return None
    return {"id": _key(["design-evaluations", identity, sorted(cohort["molecule_ids"]), sorted(cohort["group_ids"])]),
            "molecule_ids": cohort["molecule_ids"], "direct_molecule_ids": cohort["molecule_ids"], "group_ids": cohort["group_ids"]}


def _target_identity(index, target_id):
    molecule = index["molecules"].get(target_id)
    if molecule:
        components = [_component_identity(component) for component in molecule["components"]]
        if components and all(components):
            return tuple(sorted(components))
    return target_id


def connected_evaluations(index, selected, results):
    """Return exact design evaluations in the selected target context, plus unresolved count."""
    anchor = selected.get("evaluation") or {}
    molecule = index["molecules"].get(selected["molecule_id"], {})
    target = anchor.get("target_id") or ""
    if not target:
        links = molecule.get("metadata_links_by_group", {}).get(selected["group_id"], set())
        if len(links) == 1:
            target = next(iter(links))
    if not target:
        targets = {row["evaluation"]["target_id"] for row in results if row.get("evaluation", {}).get("target_id")}
        if len(targets) == 1:
            target = next(iter(targets))
    target_identity = _target_identity(index, target) if target else None
    anchor_ids = set(selected["structure_ids"])
    if selected.get("source_structure_id"):
        anchor_ids.add(selected["source_structure_id"])
    anchor_digests = {value for value in (selected.get("structure_digest"), selected.get("source_structure_digest"),
                                        anchor.get("structure_digest"), anchor.get("source_structure_digest")) if value}
    selected_digests = {value for value in (selected.get("structure_digest"), anchor.get("structure_digest")) if value}
    connected = []
    unresolved = 0
    cohort = design_subject(index, selected)
    cohort_ids = cohort["molecule_ids"] if cohort else {selected["molecule_id"]}
    identity = _design_identity(molecule, selected.get("group_id")) if molecule else None
    for row in results:
        if selected.get("target_pdb_id") and row.get("target_pdb_id") and selected["target_pdb_id"] != row["target_pdb_id"]:
            continue
        if identity is not None and _design_identity(index["molecules"][row["molecule_id"]], row.get("group_id")) != identity:
            continue
        evaluation = row.get("evaluation") or {}
        if not row["evaluation_id"]:
            continue
        candidate_target = evaluation.get("target_id") or ""
        candidate_ids = set(row["structure_ids"])
        if row.get("source_structure_id"):
            candidate_ids.add(row["source_structure_id"])
        candidate_digests = {value for value in (evaluation.get("structure_digest"), evaluation.get("source_structure_digest"),
                                               row.get("source_structure_digest")) if value}
        shared_structure = bool(anchor_ids & candidate_ids or anchor_digests & candidate_digests)
        sequence_only = row.get("source_reference_kind") == "sequence"
        current = row["evaluation_id"] == selected.get("evaluation_id")
        if candidate_target and target_identity and _target_identity(index, candidate_target) != target_identity:
            continue
        source = row.get("source_molecule_id")
        reference_digests = {value for value in (row.get("source_structure_digest"), evaluation.get("source_structure_digest")) if value}
        if not current and ((source and source not in cohort_ids) or
                            (anchor_ids & candidate_ids and reference_digests and anchor_digests
                             and not reference_digests & anchor_digests)):
            unresolved += 1
            continue
        if not current and not (shared_structure or sequence_only or candidate_target and target_identity):
            unresolved += 1
            continue
        if shared_structure:
            calculation = evaluation.get("calculation_structure_digest")
            if current or calculation in selected_digests or not calculation and set(selected["structure_ids"]) & set(row["structure_ids"]):
                relation = "Selected structure"
            elif calculation and calculation in anchor_digests or selected.get("source_structure_id") in row["structure_ids"]:
                relation = "Source model (not selected coordinates)"
            elif reference_digests & selected_digests or row.get("source_structure_id") in selected["structure_ids"]:
                relation = "Derived from selected structure"
            else:
                relation = "Other structure of this design"
        elif sequence_only:
            relation = "Sequence only"
        elif not candidate_ids and not candidate_digests and not anchor_ids and not anchor_digests:
            relation = "Selected calculation (no coordinates)" if current else "Same molecule and target (no coordinates)"
        elif candidate_target and target_identity:
            relation = "Other structure of this design"
        else:
            relation = "Structure link unverified"
        connected.append({**row, "structure_relation": relation})
    connected.sort(key=lambda row: (row["evaluation_id"] != selected.get("evaluation_id"), row["tool"].casefold(),
                                    row["run_id"], row["model"], row["evaluation_id"]))
    return connected, unresolved


def evaluation_label(result):
    evaluation = result["evaluation"]
    version = evaluation.get("engine_version") or (evaluation.get("image_digest") or "")[:19] or "version not recorded"
    representation = evaluation.get("representation_id") or evaluation.get("structure_digest") or "no structure"
    return " · ".join(str(value) for value in (result["tool"], result["run_id"], result["model"], version,
                                              representation[:16], result["structure_relation"], result["evaluation_id"][-8:]) if value)


def analysis_results(results):
    return [row for row in results if row["evaluation"].get("stage") in {"analysis", "benchmark"}
            or any(marker in row["tool"].casefold() for marker in ("ipsae", "pep_patch", "peppatch", "biotite", "sequence_qc"))]
