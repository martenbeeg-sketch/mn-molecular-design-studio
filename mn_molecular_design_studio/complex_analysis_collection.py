"""Collection-level projections of saved Database structures and Evaluations.

This module only maps saved residue identities and formats existing records. It
does not run structural calculations or infer missing Evaluation values.
"""
from __future__ import annotations

from typing import Any


def _overlapping_offsets(sequence: str, fragment: str) -> list[int]:
    if not sequence or not fragment:
        return []
    result = []
    start = 0
    while True:
        index = sequence.find(fragment, start)
        if index < 0:
            return result
        result.append(index)
        start = index + 1


def sequence_correspondence(reference: str, candidate: str) -> dict[str, Any]:
    """Map candidate sequence positions onto a reference by unique exact match.

    Exact matches and unique exact subsequences are accepted. A reference that
    occurs once inside a longer candidate is also mapped over its shared span;
    candidate overhangs remain ``None``. Repeated segments and non-identical
    sequences are reported without attempting alignment or structural fitting.
    """
    reference = "".join(str(reference or "").split()).upper()
    candidate = "".join(str(candidate or "").split()).upper()
    if not reference or not candidate:
        return {"status": "unsupported", "compatible": False, "mapping": [],
                "reason": "A component sequence is unavailable."}
    if reference == candidate:
        return {"status": "exact", "compatible": True,
                "mapping": list(range(len(candidate))), "reason": "Exact sequence match."}

    candidate_offsets = _overlapping_offsets(reference, candidate)
    if len(candidate_offsets) == 1:
        offset = candidate_offsets[0]
        return {"status": "unique_exact_subsequence", "compatible": True,
                "mapping": list(range(offset, offset + len(candidate))),
                "reason": f"Unique exact subsequence at reference index {offset + 1}; partial coverage is expected."}
    if len(candidate_offsets) > 1:
        return {"status": "ambiguous", "compatible": False, "mapping": [],
                "reason": "The component sequence occurs more than once in the reference."}

    reference_offsets = _overlapping_offsets(candidate, reference)
    if len(reference_offsets) == 1:
        offset = reference_offsets[0]
        mapping: list[int | None] = [None] * len(candidate)
        for index in range(len(reference)):
            mapping[offset + index] = index
        return {"status": "reference_is_unique_subsequence", "compatible": True,
                "mapping": mapping,
                "reason": "The reference covers one unique exact segment; candidate overhangs remain unmapped."}
    if len(reference_offsets) > 1:
        return {"status": "ambiguous", "compatible": False, "mapping": [],
                "reason": "The reference sequence occurs more than once in the component."}
    return {"status": "incompatible", "compatible": False, "mapping": [],
            "reason": "The component sequence is not an exact match or exact subsequence."}


def component_kind(components: list[dict[str, Any]]) -> str:
    """Describe molecular components without assuming design provenance."""
    kinds = {str(row.get("type") or "").casefold() for row in components if isinstance(row, dict)}
    proteins = sum(1 for row in components if isinstance(row, dict)
                   and str(row.get("type") or "").casefold() == "protein")
    ligand = bool(kinds & {"small_molecule", "ligand"})
    if proteins >= 2 and ligand:
        return "protein-protein-and-ligand"
    if proteins >= 2:
        return "protein-protein"
    if proteins and ligand:
        return "protein-ligand"
    if proteins:
        return "protein"
    return "other"


def author_residue_index(track: dict[str, Any], reference: str, chain_id: str) -> int | None:
    """Resolve an evaluator author residue reference to one sequence index."""
    from mn_molecular_design_studio.complex_analysis_inspector import parse_author_residue_selector

    parsed = parse_author_residue_selector(reference)
    if not parsed or parsed.get("chain_id") != str(chain_id or ""):
        return None
    matches = []
    for row in track.get("residues") or []:
        if int(row.get("author_number", 0)) != parsed["number"]:
            continue
        if parsed["exact_insertion"] and str(row.get("insertion_code") or "") != str(parsed.get("insertion_code") or ""):
            continue
        matches.append(row)
    return int(matches[0]["sequence_index"]) if len(matches) == 1 else None


def mapped_contact_positions(
    contact_rows: list[dict[str, Any]], *, chain_id: str, track: dict[str, Any],
    candidate_to_reference: list[int | None],
) -> set[int]:
    """Transcribe saved Biotite residue references onto a compatible sequence."""
    positions: set[int] = set()
    for row in contact_rows:
        if not isinstance(row, dict):
            continue
        for field in ("binder_residue", "target_residue"):
            local_index = author_residue_index(track, str(row.get(field) or ""), chain_id)
            if local_index is None or local_index >= len(candidate_to_reference):
                continue
            reference_index = candidate_to_reference[local_index]
            if reference_index is not None:
                positions.add(int(reference_index))
    return positions


def aggregate_contact_footprints(
    members: list[dict[str, Any]], reference_track: dict[str, Any],
) -> list[dict[str, Any]]:
    """Aggregate only saved contact maps, keeping missing and zero distinct."""
    member_positions: dict[str, set[int]] = {}
    member_status: dict[str, str] = {}
    evaluable = 0
    for member in members:
        identity = str(member.get("identity") or member.get("structure_id") or member.get("label") or "")
        rows = member.get("contact_rows")
        status = str(member.get("contact_status") or ("saved" if isinstance(rows, list) else "missing evaluation"))
        member_status[identity] = status
        if isinstance(rows, list):
            evaluable += 1
            member_positions[identity] = mapped_contact_positions(
                rows, chain_id=str(member.get("tracked_chain") or ""),
                track=member.get("tracked_track") or {},
                candidate_to_reference=member.get("candidate_to_reference") or [],
            )
        else:
            member_positions[identity] = set()

    output = []
    for residue in reference_track.get("residues") or []:
        index = int(residue.get("sequence_index", -1))
        count = sum(index in positions for positions in member_positions.values())
        output.append({
            "reference_index": index,
            "amino_acid": str(residue.get("amino_acid") or ""),
            "author_number": int(residue.get("author_number", 0)),
            "insertion_code": str(residue.get("insertion_code") or ""),
            "selector": str(residue.get("selector") or ""),
            "contact_count": count,
            "evaluable_count": evaluable,
            "missing_count": sum(status == "missing evaluation" for status in member_status.values()),
            "unavailable_count": sum(status == "artifact unavailable" for status in member_status.values()),
            "failed_count": sum(status == "failed evaluation" for status in member_status.values()),
            "frequency": count / evaluable if evaluable else None,
            "member_status": {
                identity: ("contact" if index in member_positions[identity] else
                           "no contact" if member_status[identity] == "saved" else member_status[identity])
                for identity in member_positions
            },
        })
    return output


def heatmap_selection(
    cells: list[tuple[int, str]], member_ids: list[str], residue_columns: dict[str, int],
) -> dict[str, Any] | None:
    """Resolve a Streamlit dataframe single-cell selection to a saved member."""
    if not cells:
        return None
    row_index, column = cells[0]
    try:
        row_index = int(row_index)
        member_id = str(member_ids[row_index])
    except (TypeError, ValueError, IndexError):
        return None
    residue_index = residue_columns.get(str(column))
    if residue_index is None:
        return None
    return {"member_id": member_id, "reference_index": int(residue_index), "column": str(column)}


def prepare_inspector_navigation_state(
    state: Any, *, molecule_id: str, structure_id: str, evaluation_id: str = "", selector: str = "",
) -> None:
    """Carry an exact collection selection into the existing Database Inspector."""
    structure_id = str(structure_id)
    viewer_key = f"complex-analysis-db-viewer-{structure_id}"
    state[f"{viewer_key}-inspector-highlights"] = [str(selector)] if selector else []
    state[f"complex-analysis-structure-id-{molecule_id}"] = structure_id
    state[f"complex-analysis-evaluation-{molecule_id}"] = str(evaluation_id or "")


def evaluation_table_rows(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Project exact saved Evaluation descriptors, values and source context."""
    result = []
    confidence_tokens = ("plddt", "pae", "ptm", "iptm", "ipae", "confidence")
    for entry in entries:
        observation = entry.get("observation") if isinstance(entry.get("observation"), dict) else {}
        descriptor = entry.get("evaluation") if isinstance(entry.get("evaluation"), dict) else {}
        context = entry.get("context") if isinstance(entry.get("context"), dict) else {}
        candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        settings = context.get("settings") if isinstance(context.get("settings"), dict) else {}
        metrics = observation.get("metrics") if isinstance(observation.get("metrics"), dict) else {}
        evaluator = str(descriptor.get("evaluator") or observation.get("namespace") or "unknown")
        confidence = any(any(token in str(key).casefold() for token in confidence_tokens) for key in metrics)
        result.append({
            "molecule_id": str(entry.get("molecule_id") or observation.get("molecule_id") or ""),
            "structure_id": str(entry.get("structure_id") or descriptor.get("representation_id") or ""),
            "structure_label": str(entry.get("structure_label") or ""),
            "evaluation_id": str(descriptor.get("id") or ""),
            "observation_id": str(observation.get("id") or descriptor.get("observation_id") or ""),
            "method": evaluator,
            "category": "Prediction confidence" if confidence else "Scientific evaluation",
            "source_engine": str(descriptor.get("engine_version") or raw.get("engine_version") or settings.get("engine_version") or evaluator),
            "run_id": str(descriptor.get("run_id") or observation.get("run_id") or context.get("run_id") or ""),
            "candidate_id": str(descriptor.get("candidate_id") or context.get("candidate_id") or candidate.get("candidate_id") or ""),
            "structure_digest": str(context.get("source_structure_digest") or context.get("structure_digest") or ""),
            "parameters": settings,
            "provenance": context,
            "metrics": metrics,
            "created_at": str(descriptor.get("created_at") or observation.get("created_at") or ""),
            "status": str(descriptor.get("status") or context.get("status") or "completed"),
        })
    return result


def filter_evaluation_rows(
    rows: list[dict[str, Any]], *, method: str = "", source_engine: str = "",
    metric: str = "", minimum: float | None = None, maximum: float | None = None,
) -> list[dict[str, Any]]:
    """Filter saved evaluations by recorded method, engine and metric value."""
    output = []
    for row in rows:
        if method and row.get("method") != method:
            continue
        if source_engine and row.get("source_engine") != source_engine:
            continue
        if metric:
            values = row.get("metrics") if isinstance(row.get("metrics"), dict) else {}
            value = values.get(metric)
            if type(value) not in {int, float}:
                continue
            if minimum is not None and float(value) < float(minimum):
                continue
            if maximum is not None and float(value) > float(maximum):
                continue
        output.append(row)
    return output


def collection_page(items: list[Any], *, page: int, page_size: int) -> list[Any]:
    """Return one stable display page while keeping the full collection intact."""
    size = max(1, int(page_size))
    start = max(0, int(page)) * size
    return items[start:start + size]


def assignment_rows(evaluations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Expose assignments already saved by registered clustering workflows."""
    output = []
    def assignment_kind(method: str) -> tuple[tuple[str, ...], tuple[str, ...]] | None:
        normalized = method.casefold()
        if "interface_clustering" in normalized:
            return (("interface_cluster", "cluster_id", "cluster", "diversity_cluster", "diversity.cluster"),
                    ("is_representative", "diversity_representative", "diversity.representative"))
        if "foldseek" in normalized:
            return (("foldseek_cluster", "cluster_id", "cluster"), ("is_representative", "representative"))
        return None

    # Some saved collection evaluations store only the per-member representative
    # flag, not the representative candidate ID. Resolve that case only inside
    # the same method, run, and cluster; identical labels from separate runs are
    # independent assignments and must not be joined.
    group_representatives: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for candidate in evaluations:
        method = str(candidate.get("method") or "")
        keys = assignment_kind(method)
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        if not keys:
            continue
        cluster_keys, representative_keys = keys
        cluster = next((metrics[key] for key in cluster_keys if key in metrics), None)
        representative = next((metrics[key] for key in representative_keys if key in metrics), None)
        if cluster is None or representative is not True:
            continue
        group = (method.casefold(), str(candidate.get("run_id") or ""), str(cluster))
        group_representatives.setdefault(group, []).append(candidate)

    for row in evaluations:
        method = str(row.get("method") or "")
        normalized = method.casefold()
        keys = assignment_kind(method)
        if not keys:
            continue
        if "interface_clustering" in normalized:
            kind = "Interface Clustering · binding-mode assignment"
        elif "foldseek" in normalized:
            kind = "Foldseek · structural-similarity assignment"
        cluster_keys, representative_keys = keys
        metrics = row.get("metrics") if isinstance(row.get("metrics"), dict) else {}
        cluster = next((metrics[key] for key in cluster_keys if key in metrics), None)
        representative = next((metrics[key] for key in representative_keys if key in metrics), None)
        representative_candidate_id = str(metrics.get("diversity_representative_id") or metrics.get("representative_id") or "")
        candidate_matches = [candidate for candidate in evaluations
                             if representative_candidate_id
                             and str(candidate.get("candidate_id") or "") == representative_candidate_id
                             and str(candidate.get("method") or "").casefold() == normalized
                             and str(candidate.get("run_id") or "") == str(row.get("run_id") or "")]
        representative_entry = candidate_matches[0] if len(candidate_matches) == 1 else None
        representative_status = "identified by saved candidate ID" if representative_entry else "not identified"
        if representative is True and representative_entry is None:
            representative_entry = row
            representative_status = "saved representative member"
        elif representative_entry is None and cluster is not None:
            group = (normalized, str(row.get("run_id") or ""), str(cluster))
            candidates = group_representatives.get(group, [])
            if len(candidates) == 1:
                representative_entry = candidates[0]
                representative_status = "unique saved representative in this run and cluster"
            elif len(candidates) > 1:
                representative_status = "ambiguous: multiple saved representatives in this run and cluster"
        if cluster is None and representative is None:
            continue
        output.append({
            "Structure": row.get("structure_label") or row.get("structure_id"),
            "Molecule ID": row.get("molecule_id"),
            "Structure ID": row.get("structure_id"),
            "Candidate ID": row.get("candidate_id"),
            "Method": kind,
            "Assignment": cluster if cluster is not None else "Not recorded",
            "Representative": representative if representative is not None else "Not recorded",
            "Representative candidate ID": representative_candidate_id or "Not recorded",
            "Representative structure": (representative_entry.get("structure_label") if representative_entry
                                         else "Not uniquely identified in this run and cluster"),
            "Representative molecule ID": (representative_entry.get("molecule_id") if representative_entry
                                            else ""),
            "Representative structure ID": (representative_entry.get("structure_id") if representative_entry
                                             else ""),
            "Representative resolution": representative_status,
            "Evaluation ID": row.get("evaluation_id"),
            "Run ID": row.get("run_id"),
            "Parameters": row.get("parameters") or {},
        })
    return output
