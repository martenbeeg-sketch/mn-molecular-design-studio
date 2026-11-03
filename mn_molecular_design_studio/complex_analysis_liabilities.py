"""Post-hoc, residue-linked sequence and developability liability screening.

Definitions for sequence motifs and SG-distance disulfide inference are adapted
from BinderScope's ``binder_review.liabilities`` and ``disulfides`` modules
(MIT; see ``docs/third-party-notices.md``). Studio applies them to explicitly
selected Database chains and retains its own structure, Evaluation, access and
residue-numbering provenance. This is an advisory screen, not a stability,
modification, expression or binding prediction.
"""
from __future__ import annotations

import ast
import csv
import itertools
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np

METHOD = "Studio sequence and structure liability screen (BinderScope liabilities-inspired)"
METHOD_VERSION = "studio.sequence-liabilities.v1"
DISULFIDE_CLOSE_MAX_ANGSTROM = 3.0
DISULFIDE_LOOSE_MAX_ANGSTROM = 3.4
DISULFIDE_MATCH_TIE_TOLERANCE_ANGSTROM = 0.01
EXPOSED_RSA_THRESHOLD = 0.25
DIBASIC_RSA_THRESHOLD = 0.40
HYDROPHOBIC_STRETCH_MINIMUM = 6
_STANDARD = frozenset("ACDEFGHIKLMNPQRSTVWY")
_HYDROPHOBIC = frozenset("VILMFWY")
_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}
_MOTIF_DEFINITIONS = {
    "n_glycosylation_sequon": "N-X-S/T; X is not proline. Sequence motif only; occupancy is not measured.",
    "deamidation_ng": "Asn-Gly sequence motif; susceptibility is not an observed chemical modification.",
    "deamidation_ns": "Asn-Ser sequence motif; susceptibility is not an observed chemical modification.",
    "isomerization_dg": "Asp-Gly sequence motif; propensity is not an observed isomerization.",
    "isomerization_ds": "Asp-Ser sequence motif; propensity is not an observed isomerization.",
    "acid_cleavage_dp": "Asp-Pro sequence motif; it is not an observed cleavage.",
    "furin_consensus": "R-X-K/R-R sequence motif; it is not evidence of protease exposure or cleavage.",
    "dibasic_protease_motif": "Adjacent Lys/Arg residues; it is not evidence of protease exposure or cleavage.",
    "hydrophobic_stretch": "Six or more consecutive V/I/L/M/F/W/Y residues; exposure is evaluated separately when available.",
}


def _identity(chain_id: str, residue) -> dict[str, Any]:
    number, insertion = int(residue.id[1]), str(residue.id[2] or "").strip()
    from mn_molecular_design_studio.complex_analysis_inspector import author_residue_selector
    return {
        "chain_id": str(chain_id), "author_number": number, "insertion_code": insertion,
        "author_residue": f"{number}{insertion}",
        "residue_name": str(residue.resname).strip().upper(),
        "selector": author_residue_selector(str(chain_id), number, insertion),
    }


def _coordinate(atom) -> np.ndarray | None:
    try:
        value = np.asarray(atom.get_coord(), dtype=float)
    except (AttributeError, TypeError, ValueError):
        return None
    return value if value.shape == (3,) and np.isfinite(value).all() else None


def _selected_residues(model, chains: list[str]) -> dict[tuple[str, int, str], dict[str, Any]]:
    from Bio.PDB.Polypeptide import is_aa
    output = {}
    for chain_id in chains:
        if chain_id not in model:
            raise ValueError(f"Selected protein chain {chain_id!r} is absent from the selected model.")
        for residue in model[chain_id].get_residues():
            if not is_aa(residue, standard=False):
                continue
            letter = _AA3_TO_1.get(str(residue.resname).strip().upper(), "X")
            ident = _identity(chain_id, residue)
            ident.update({"amino_acid": letter, "residue_id": residue.id, "residue_object": residue})
            output[(chain_id, ident["author_number"], ident["insertion_code"])] = ident
    return output


def _track_for_chain(path: Path, chain_id: str, model) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Use file-declared polymer sequence when unambiguous, else observed residues."""
    from mn_molecular_design_studio.complex_analysis_reference_context import reference_component_track
    from mn_molecular_design_studio.complex_analysis_inspector import protein_sequence_tracks

    track = reference_component_track(path, chain_id, [])
    declared_rows = list(track.get("residues") or []) if track.get("status") == "available" else []
    observed_tracks = protein_sequence_tracks(path, [chain_id])
    observed_rows = list(observed_tracks[0].get("residues") or []) if observed_tracks else []
    if declared_rows:
        rows = declared_rows
        source = "file-declared polymer sequence with exact author mapping where unambiguous"
        status = "declared_sequence"
    else:
        rows = observed_rows
        source = "coordinate-observed residues only; unresolved construct sequence is unknown"
        status = "observed_sequence_only" if rows else "unavailable"
    if not rows:
        return [], {"status": "unavailable", "reason": track.get("reason") or "No supported protein sequence is available."}

    coord_by_key = {}
    if chain_id in model:
        for residue in model[chain_id]:
            if str(residue.resname).strip().upper() not in _AA3_TO_1:
                continue
            key = (int(residue.id[1]), str(residue.id[2] or "").strip())
            coord_by_key.setdefault(key, residue)
    for row in rows:
        row["selector"] = ""
        row["coordinates_present"] = False
        row["relative_sasa"] = None
        row["absolute_sasa_angstrom2"] = None
        row["secondary_structure"] = None
        row["ca_b_factor"] = None
        number = row.get("author_number")
        if number is None:
            continue
        key = (int(number), str(row.get("insertion_code") or "").strip())
        residue = coord_by_key.get(key)
        if residue is None or _AA3_TO_1.get(str(residue.resname).strip().upper(), "X") != row["amino_acid"]:
            continue
        from mn_molecular_design_studio.complex_analysis_inspector import author_residue_selector
        row["selector"] = author_residue_selector(chain_id, int(number), key[1])
        row["coordinates_present"] = True
        row["residue_id"] = residue.id
        row["residue_object"] = residue
        ca = residue["CA"] if "CA" in residue else None
        row["ca_b_factor"] = float(ca.get_bfactor()) if ca is not None and math.isfinite(float(ca.get_bfactor())) else None
    unresolved = sum(not row["coordinates_present"] for row in rows)
    return rows, {
        "status": status,
        "sequence_source": source,
        "sequence_length": len(rows),
        "coordinate_residue_count": len(rows) - unresolved,
        "unresolved_or_unmapped_residue_count": unresolved,
        "complete_coordinate_coverage": unresolved == 0,
        "reason": "" if unresolved == 0 else "Some declared sequence positions have no uniquely mapped coordinates; structural findings are unavailable for them.",
    }


def _recorded_pairs(path: Path, model) -> list[tuple[tuple[str, int, str], tuple[str, int, str], str]]:
    """Read explicit PDB SSBOND/LINK and mmCIF disulfide struct_conn records."""
    result: list[tuple[tuple[str, int, str], tuple[str, int, str], str]] = []
    suffix = path.name.casefold()
    if suffix.endswith((".cif", ".mmcif")):
        try:
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict
            data = MMCIF2Dict(str(path))
            def values(key):
                value = data.get(key, [])
                return list(value) if isinstance(value, (list, tuple)) else ([value] if value else [])
            conn_types = values("_struct_conn.conn_type_id")
            n = len(conn_types)
            for index, conn_type in enumerate(conn_types):
                if "disulf" not in str(conn_type).casefold():
                    continue
                endpoints = []
                for prefix in ("ptnr1", "ptnr2"):
                    chain_values = values(f"_struct_conn.{prefix}_auth_asym_id") or values(f"_struct_conn.{prefix}_label_asym_id")
                    number_values = values(f"_struct_conn.{prefix}_auth_seq_id") or values(f"_struct_conn.{prefix}_label_seq_id")
                    comp_values = values(f"_struct_conn.{prefix}_auth_comp_id") or values(f"_struct_conn.{prefix}_label_comp_id")
                    atom_values = values(f"_struct_conn.{prefix}_auth_atom_id") or values(f"_struct_conn.{prefix}_label_atom_id")
                    insert_values = values(f"_struct_conn.pdbx_{prefix}_PDB_ins_code")
                    chain = str(chain_values[index]).strip() if index < len(chain_values) else ""
                    number = int(str(number_values[index]).strip()) if index < len(number_values) and str(number_values[index]).strip() not in {".", "?"} else None
                    comp = str(comp_values[index]).strip().upper() if index < len(comp_values) else ""
                    atom = str(atom_values[index]).strip().upper() if index < len(atom_values) else ""
                    insertion = str(insert_values[index]).strip() if index < len(insert_values) and str(insert_values[index]).strip() not in {".", "?"} else ""
                    if number is None or comp != "CYS" or atom not in {"SG", "?", "."}:
                        endpoints = []
                        break
                    endpoints.append((chain, number, insertion))
                if len(endpoints) == 2:
                    result.append((endpoints[0], endpoints[1], "mmCIF _struct_conn disulfide record"))
        except (OSError, ValueError, KeyError, IndexError, ImportError):
            return []
        return result

    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    for line in lines:
        if line.startswith("SSBOND") and len(line) >= 36:
            try:
                first = (line[15:16].strip() or "_", int(line[17:21]), line[21:22].strip())
                second = (line[29:30].strip() or "_", int(line[31:35]), line[35:36].strip())
                result.append((first, second, "PDB SSBOND record"))
            except (ValueError, IndexError):
                continue
        elif line.startswith("LINK") and len(line) >= 57:
            if line[12:16].strip().upper() != "SG" or line[42:46].strip().upper() != "SG":
                continue
            if line[17:20].strip().upper() != "CYS" or line[47:50].strip().upper() != "CYS":
                continue
            try:
                first = (line[21:22].strip() or "_", int(line[22:26]), line[26:27].strip())
                second = (line[51:52].strip() or "_", int(line[52:56]), line[56:57].strip())
                result.append((first, second, "PDB LINK SG–SG record"))
            except (ValueError, IndexError):
                continue
    return result


def _optimal_matching(nodes: list[tuple[str, int, str]], edges: dict[tuple[tuple[str, int, str], tuple[str, int, str]], float]):
    """Return BinderScope-style maximum-margin matching and detect tied optima."""
    memo: dict[tuple[tuple[str, int, str], ...], tuple[float, int, list[tuple[tuple[str, int, str], tuple[str, int, str]]]]] = {}

    def solve(remaining):
        if not remaining:
            return 0.0, 1, []
        if remaining in memo:
            return memo[remaining]
        first, rest = remaining[0], remaining[1:]
        candidates = []
        score, count, pairs = solve(rest)
        candidates.append((score, count, pairs))  # leave first cysteine unpaired
        for other in rest:
            edge = (first, other) if first < other else (other, first)
            if edge not in edges:
                continue
            next_remaining = tuple(item for item in rest if item != other)
            child_score, child_count, child_pairs = solve(next_remaining)
            candidates.append((child_score + DISULFIDE_LOOSE_MAX_ANGSTROM - edges[edge], child_count, [edge, *child_pairs]))
        best_score = max(row[0] for row in candidates)
        best = [row for row in candidates if abs(row[0] - best_score) <= DISULFIDE_MATCH_TIE_TOLERANCE_ANGSTROM]
        result = (best_score, min(2, sum(row[1] for row in best)), best[0][2])
        memo[remaining] = result
        return result

    return solve(tuple(sorted(nodes)))


def _disulfide_assessment(path: Path, model, component_chains: list[str], residue_index: dict) -> dict[str, Any]:
    selected_chain_set = set(component_chains)
    selected_cysteines = {key: value for key, value in residue_index.items()
                          if key[0] in selected_chain_set and value["amino_acid"] == "C"}
    if not selected_cysteines:
        return {"status": "not_present", "basis": "none", "pairs": [], "findings": [],
                "reason": "No cysteine residues occur in the selected sequence component."}
    # Geometry must include cysteine partners on other modeled protein chains,
    # while findings remain scoped to the explicitly selected component.
    from Bio.PDB.Polypeptide import is_aa
    all_cysteines = dict(selected_cysteines)
    for chain in model:
        for residue in chain:
            if not is_aa(residue, standard=False) or str(residue.resname).strip().upper() != "CYS":
                continue
            key = (str(chain.id), int(residue.id[1]), str(residue.id[2] or "").strip())
            if key in all_cysteines:
                continue
            identity = _identity(str(chain.id), residue)
            identity.update({"amino_acid": "C", "residue_id": residue.id, "residue_object": residue})
            all_cysteines[key] = identity
    records = _recorded_pairs(path, model)
    relevant_records = []
    for first, second, record_method in records:
        # A recorded covalent partner outside the selected component still
        # resolves the selected cysteine's explicit pairing identity.
        if first in selected_cysteines or second in selected_cysteines:
            relevant_records.append((first, second, record_method))
    pairs, findings = [], []
    if relevant_records:
        paired = set()
        for first, second, record_method in relevant_records:
            paired.update((first, second))
            pairs.append({"first": first, "second": second, "distance_angstrom": _distance_sg(model, first, second),
                          "basis": "recorded", "record": record_method})
        for key, cys in selected_cysteines.items():
            if key not in paired:
                residue = cys.get("residue_object")
                atom = residue["SG"] if residue is not None and "SG" in residue else None
                if _coordinate(atom) is None:
                    findings.append(_finding("cysteine_sg_coordinates_missing", "info", "sequence_only", cys,
                        "Cysteine is in the selected sequence but has no finite SG coordinate; pairing status is unknown.",
                        evidence={"disulfide_basis": "explicit_record_present_but_no_pair_recorded"}))
                    continue
                findings.append(_finding(
                    "unpaired_cysteine", "warn", "structure_supported", cys,
                    "No disulfide involving this cysteine is present in the structure's explicit disulfide records. This is a file-record observation, not proof of redox state.",
                    evidence={"disulfide_basis": "recorded", "recorded_pair_count": len(relevant_records)},
                ))
        return {"status": "recorded", "basis": "explicit_structure_record", "pairs": pairs,
                "findings": findings, "record_count": len(relevant_records),
                "reason": "Explicit SSBOND, SG–SG LINK or mmCIF disulfide struct_conn records take precedence over geometry."}

    coordinate_sg = {}
    missing_sg = []
    for key, cys in all_cysteines.items():
        residue = cys.get("residue_object")
        if residue is None:
            if key in selected_cysteines:
                missing_sg.append((key, cys))
            continue
        atom = residue["SG"] if "SG" in residue else None
        xyz = _coordinate(atom) if atom is not None else None
        if xyz is None:
            if key in selected_cysteines:
                missing_sg.append((key, cys))
        else:
            coordinate_sg[key] = xyz
    if not coordinate_sg:
        missing_findings = [_finding("cysteine_sg_coordinates_missing", "info", "sequence_only", cys,
            "Cysteine is in the selected sequence but no finite SG coordinate is present; pairing status is unknown.",
            evidence={"disulfide_basis": "unavailable"}) for _key, cys in missing_sg]
        return {"status": "unavailable", "basis": "none", "pairs": [], "findings": missing_findings,
                "missing_sg": [item[0] for item in missing_sg],
                "reason": "No explicit disulfide records or usable Cys SG coordinates are present; unpaired status is unknown."}

    edges = {}
    adjacency = {key: set() for key in coordinate_sg}
    for first, second in itertools.combinations(sorted(coordinate_sg), 2):
        distance = float(np.linalg.norm(coordinate_sg[first] - coordinate_sg[second]))
        if distance <= DISULFIDE_LOOSE_MAX_ANGSTROM:
            edges[(first, second)] = distance
            adjacency[first].add(second)
            adjacency[second].add(first)
    components, seen = [], set()
    for node in sorted(adjacency):
        if node in seen or not adjacency[node]:
            continue
        stack, component = [node], []
        seen.add(node)
        while stack:
            current = stack.pop()
            component.append(current)
            for adjacent in adjacency[current]:
                if adjacent not in seen:
                    seen.add(adjacent)
                    stack.append(adjacent)
        components.append(sorted(component))

    unresolved = set()
    paired = set()
    for component in components:
        if len(component) > 14:
            unresolved.update(component)
            for key in component:
                if key not in selected_cysteines:
                    continue
                findings.append(_finding("ambiguous_disulfide_pairing", "info", "structure_supported", selected_cysteines[key],
                    "More than 14 cysteine sulfur atoms form one local candidate network; this implementation does not guess a matching.",
                    evidence={"candidate_network_cysteines": len(component), "disulfide_basis": "inferred_geometry"}))
            continue
        _score, count, selected = _optimal_matching(component, edges)
        if count > 1:
            unresolved.update(component)
            for key in component:
                if key not in selected_cysteines:
                    continue
                alternatives = sorted(round(edges[pair], 2) for pair in edges if key in pair and pair in selected)
                findings.append(_finding("ambiguous_disulfide_pairing", "info", "structure_supported", selected_cysteines[key],
                    "Multiple one-to-one SG-distance matchings are equally supported within the 0.01 Å matching tolerance; no disulfide or free-cysteine assignment is made.",
                    evidence={"candidate_network_cysteines": len(component), "optimal_matching_count": count,
                              "disulfide_basis": "inferred_geometry", "matched_candidate_distances_angstrom": alternatives}))
            continue
        for first, second in selected:
            paired.update((first, second))
            if first not in selected_cysteines and second not in selected_cysteines:
                continue
            distance = edges[(first, second)]
            pairs.append({"first": first, "second": second, "distance_angstrom": round(distance, 2),
                          "basis": "inferred_geometry",
                          "geometry_class": "close" if distance <= DISULFIDE_CLOSE_MAX_ANGSTROM else "loose"})
            if distance > DISULFIDE_CLOSE_MAX_ANGSTROM:
                for key in (first, second):
                    if key not in selected_cysteines:
                        continue
                    findings.append(_finding("loose_disulfide_geometry", "info", "structure_supported", selected_cysteines[key],
                        f"This cysteine is paired by the selected one-to-one SG geometry at {distance:.2f} Å (loose range {DISULFIDE_CLOSE_MAX_ANGSTROM:.1f}–{DISULFIDE_LOOSE_MAX_ANGSTROM:.1f} Å). This is model geometry, not an explicit covalent-bond record.",
                        evidence={"partner": _endpoint(second if key == first else first, residue_index),
                                  "sg_sg_distance_angstrom": round(distance, 2), "disulfide_basis": "inferred_geometry"}))
    for key, cys in selected_cysteines.items():
        if key not in coordinate_sg or key in paired or key in unresolved:
            continue
        nearest = min((float(np.linalg.norm(coordinate_sg[key] - xyz)) for other, xyz in coordinate_sg.items() if other != key), default=None)
        reason = "No SG partner lies within the 3.4 Å geometry-assignment limit among modeled protein chains. This is a model-based candidate for an unpaired cysteine, not an experimental redox measurement."
        findings.append(_finding("unpaired_cysteine", "warn", "structure_supported", cys, reason,
            evidence={"disulfide_basis": "inferred_geometry", "nearest_other_sg_distance_angstrom": round(nearest, 2) if nearest is not None else None,
                      "geometry_pairing_limit_angstrom": DISULFIDE_LOOSE_MAX_ANGSTROM}))
    for key, cys in missing_sg:
        findings.append(_finding("cysteine_sg_coordinates_missing", "info", "sequence_only", cys,
            "Cysteine is in the selected sequence but no finite SG coordinate is present; pairing status is unknown.",
            evidence={"disulfide_basis": "unavailable"}))
    return {"status": "inferred_geometry", "basis": "modeled_protein_sg_geometry", "pairs": pairs,
            "findings": findings, "missing_sg": [item[0] for item in missing_sg],
            "record_count": 0,
            "method_parameters": {"close_max_angstrom": DISULFIDE_CLOSE_MAX_ANGSTROM,
                                  "loose_max_angstrom": DISULFIDE_LOOSE_MAX_ANGSTROM,
                                  "matching_objective": "maximum summed loose-limit margin, one-to-one",
                                  "ambiguity_tolerance_angstrom": DISULFIDE_MATCH_TIE_TOLERANCE_ANGSTROM,
                                  "maximum_exact_network_size": 14}}


def _distance_sg(model, first, second) -> float | None:
    endpoints = []
    for chain_id, number, insertion in (first, second):
        if chain_id not in model:
            return None
        residue = next((res for res in model[chain_id]
                        if int(res.id[1]) == int(number) and str(res.id[2] or "").strip() == str(insertion or "")), None)
        if residue is None or "SG" not in residue:
            return None
        coord = _coordinate(residue["SG"])
        if coord is None:
            return None
        endpoints.append(coord)
    return round(float(np.linalg.norm(endpoints[0] - endpoints[1])), 2)


def _endpoint(key, residue_index):
    item = residue_index.get(key)
    if item is not None:
        return {name: item.get(name) for name in ("chain_id", "author_number", "insertion_code", "author_residue", "residue_name", "selector")}
    chain, number, insertion = key
    return {"chain_id": chain, "author_number": number, "insertion_code": insertion,
            "author_residue": f"{number}{insertion}", "residue_name": "CYS", "selector": ""}


def _finding(kind: str, severity: str, basis: str, residue: dict, message: str, *, evidence=None, method=None):
    return {
        "kind": kind, "severity": severity, "basis": basis, "chain_id": residue.get("chain_id"),
        "sequence_index": residue.get("sequence_index"), "author_number": residue.get("author_number"),
        "insertion_code": residue.get("insertion_code") or "", "author_residue": residue.get("author_residue"),
        "residue_name": residue.get("residue_name"), "amino_acid": residue.get("amino_acid"),
        "selector": residue.get("selector") or "", "motif": "", "message": message,
        "interface_status": "unknown", "evidence": evidence or {}, "method": method or METHOD,
    }


def _safe_protein_metrics(sequence: str, saved: dict | None = None) -> dict[str, Any]:
    if saved and str(saved.get("sequence") or "") == sequence:
        metrics = saved.get("metrics") if isinstance(saved.get("metrics"), dict) else {}
        required = ("qc_sequence_pi", "qc_sequence_charge", "qc_sequence_gravy")
        saved_ph = _number(metrics.get("qc_sequence_ph"))
        if (all(metrics.get(key) is not None for key in required)
                and saved_ph is not None and abs(saved_ph - 7.4) <= 1e-6):
            return {"status": "reused_saved_evaluation", "evaluator": "sequence_qc",
                    "evaluation_id": saved.get("evaluation_id"), "run_id": saved.get("run_id"),
                    "metrics": {key: metrics.get(key) for key in ("qc_sequence_length", *required)},
                    "method": "Studio Sequence Composition / Biopython ProteinAnalysis", "sequence": sequence}
    if not sequence or set(sequence) - _STANDARD:
        return {"status": "unavailable", "reason": "Sequence Composition requires a complete sequence of standard amino acids; unknown or non-standard residues are retained but not silently discarded.",
                "sequence": sequence}
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_runner import sequence_metrics
    metrics = sequence_metrics(sequence, 7.4)
    return {"status": "measured", "evaluator": "sequence_qc-compatible calculation",
            "evaluation_id": "", "run_id": "", "metrics": {key: metrics[key] for key in (
                "qc_sequence_length", "qc_sequence_pi", "qc_sequence_charge", "qc_sequence_gravy")},
            "method": "Studio Sequence Composition; Biopython ProteinAnalysis", "sequence": sequence}


def _residue_findings(sequence: str, rows: list[dict], exposure: dict[int, float | None], interface_indices: set[int],
                      *, exposure_status: str) -> list[dict]:
    findings = []

    def emit(kind, start, length, severity, basis, message, evidence=None):
        for offset in range(length):
            idx = start + offset
            residue = dict(rows[idx])
            residue["sequence_index"] = idx
            finding = _finding(kind, severity, basis, residue, message, evidence=evidence,
                               method=METHOD)
            finding["motif"] = sequence[start:start + length]
            finding["sequence_position"] = idx + 1
            finding["interface_status"] = "in_interface" if idx in interface_indices else "not_in_interface" if interface_indices is not None else "unknown"
            finding["evidence"].setdefault("sequence_position_start", start + 1)
            finding["evidence"].setdefault("sequence_position_end", start + length)
            findings.append(finding)

    n = len(sequence)
    # Scan tri- and dipeptide motifs separately so the final two residues are
    # still considered for NG/NS/DG/DS/DP alerts.
    for index in range(max(0, n - 2)):
        tri = sequence[index:index + 3]
        if tri[0] == "N" and tri[1] in (_STANDARD - {"P"}) and tri[2] in "ST":
            emit("n_glycosylation_sequon", index, 3, "info", "sequence_only",
                 "Potential N-X-S/T sequon. This does not establish glycan occupancy; no sugar coordinates or expression context are inferred.",
                 {"motif_definition": _MOTIF_DEFINITIONS["n_glycosylation_sequon"]})
    for index in range(max(0, n - 1)):
        di = sequence[index:index + 2]
        if di == "NG":
            emit("deamidation_ng", index, 2, "warn", "sequence_only",
                 "NG motif associated with higher sequence-based deamidation susceptibility; modification is not measured.",
                 {"motif_definition": _MOTIF_DEFINITIONS["deamidation_ng"]})
        if di == "NS":
            emit("deamidation_ns", index, 2, "info", "sequence_only",
                 "NS motif associated with sequence-based deamidation susceptibility; modification is not measured.",
                 {"motif_definition": _MOTIF_DEFINITIONS["deamidation_ns"]})
        if di == "DG":
            emit("isomerization_dg", index, 2, "warn", "sequence_only",
                 "DG motif associated with sequence-based Asp isomerization susceptibility; isomerization is not measured.",
                 {"motif_definition": _MOTIF_DEFINITIONS["isomerization_dg"]})
        if di == "DS":
            emit("isomerization_ds", index, 2, "info", "sequence_only",
                 "DS motif associated with sequence-based Asp isomerization susceptibility; isomerization is not measured.",
                 {"motif_definition": _MOTIF_DEFINITIONS["isomerization_ds"]})
        if di == "DP":
            emit("acid_cleavage_dp", index, 2, "info", "sequence_only",
                 "DP motif can be acid-sensitive in sequence screens; cleavage is not measured.",
                 {"motif_definition": _MOTIF_DEFINITIONS["acid_cleavage_dp"]})
        if all(letter in "KR" for letter in di):
            rsa = [exposure.get(index), exposure.get(index + 1)]
            if all(value is not None for value in rsa):
                if min(rsa) >= DIBASIC_RSA_THRESHOLD:
                    emit("dibasic_protease_motif", index, 2, "info", "sequence_and_structure",
                         "Dibasic sequence motif is solvent exposed in the supplied model; cleavage is not measured.",
                         {"motif_definition": _MOTIF_DEFINITIONS["dibasic_protease_motif"],
                          "relative_sasa": rsa, "relative_sasa_threshold": DIBASIC_RSA_THRESHOLD,
                          "accessibility_method": "DSSP Sander scale"})
                else:
                    emit("dibasic_protease_motif", index, 2, "info", "sequence_and_structure",
                         "Dibasic sequence motif is below the saved-model exposure threshold; this does not establish protease resistance or cleavage.",
                         {"motif_definition": _MOTIF_DEFINITIONS["dibasic_protease_motif"],
                          "relative_sasa": rsa, "relative_sasa_threshold": DIBASIC_RSA_THRESHOLD,
                          "exposure_status": "below_threshold", "accessibility_method": "DSSP Sander scale"})
            else:
                emit("dibasic_protease_motif", index, 2, "info", "sequence_only",
                     "Potential dibasic protease motif. Exposure is unavailable; cleavage is not measured.",
                     {"motif_definition": _MOTIF_DEFINITIONS["dibasic_protease_motif"],
                      "relative_sasa": rsa, "accessibility_status": exposure_status})
    for match in re.finditer(r"R[ACDEFGHIKLMNPQRSTVWY][KR]R", sequence):
        index = match.start()
        rsa = exposure.get(index)
        basis = "sequence_and_structure" if rsa is not None else "sequence_only"
        emit("furin_consensus", index, 4, "warn" if rsa is None or rsa >= EXPOSED_RSA_THRESHOLD else "info", basis,
             "Potential R-X-K/R-R sequence consensus. DSSP exposure is shown where available; cleavage is not measured.",
             {"motif_definition": _MOTIF_DEFINITIONS["furin_consensus"], "relative_sasa": rsa,
              "accessibility_method": "DSSP Sander scale" if rsa is not None else "unavailable"})
    for index, residue in enumerate(sequence):
        if residue not in {"M", "W"}:
            continue
        rsa = exposure.get(index)
        if rsa is None:
            message = f"{residue} has sequence-based oxidation susceptibility; solvent exposure is unknown in this analysis. Oxidation is not measured."
            severity, basis = "info", "sequence_only"
        elif rsa >= EXPOSED_RSA_THRESHOLD:
            message = f"{residue} is relatively exposed in this model (DSSP relative SASA {rsa:.2f}) and may be more oxidation-susceptible; oxidation is not measured."
            severity, basis = ("warn" if residue == "M" else "info"), "sequence_and_structure"
        else:
            message = f"{residue} is below the DSSP relative-SASA exposure threshold in this model ({rsa:.2f}); this is not evidence of zero oxidation risk."
            severity, basis = "info", "sequence_and_structure"
        emit("methionine_oxidation_susceptibility" if residue == "M" else "tryptophan_oxidation_susceptibility",
             index, 1, severity, basis, message,
             {"relative_sasa": rsa, "exposure_threshold": EXPOSED_RSA_THRESHOLD,
              "accessibility_method": "DSSP Sander scale" if rsa is not None else "unavailable"})
    index = 0
    while index < n:
        if sequence[index] not in _HYDROPHOBIC:
            index += 1
            continue
        end = index + 1
        while end < n and sequence[end] in _HYDROPHOBIC:
            end += 1
        if end - index >= HYDROPHOBIC_STRETCH_MINIMUM:
            known = [exposure.get(pos) for pos in range(index, end)]
            if all(value is not None for value in known):
                average = sum(known) / len(known)
                severity = "warn" if average >= EXPOSED_RSA_THRESHOLD else "info"
                basis = "sequence_and_structure"
                message = ("Long hydrophobic run is relatively exposed in this supplied model; this is a surface-risk observation, not an aggregation prediction."
                           if average >= EXPOSED_RSA_THRESHOLD else
                           "Long hydrophobic run is mostly below the DSSP exposure threshold in this supplied model; packing and alternate conformations are not assessed.")
            else:
                average, severity, basis = None, "info", "sequence_only"
                message = "Long hydrophobic run is a sequence-only alert because residue accessibility is incomplete or unavailable."
            emit("hydrophobic_stretch", index, end - index, severity, basis, message,
                 {"motif_definition": _MOTIF_DEFINITIONS["hydrophobic_stretch"],
                  "mean_relative_sasa": average, "exposure_threshold": EXPOSED_RSA_THRESHOLD,
                  "accessibility_method": "DSSP Sander scale" if average is not None else "unavailable"})
        index = end
    for index, residue in enumerate(sequence):
        if residue in {"K", "M"} and index in interface_indices:
            emit("interface_lysine" if residue == "K" else "interface_methionine", index, 1, "info", "structure_supported",
                 "This residue contacts the selected interface in the saved Biotite contact artifact. Chemical modification or its effect on binding is not measured.",
                 {"interface_evidence": "saved Biotite Interface contact artifact"})
    return findings


def _parse_csv(path: Path) -> list[dict[str, Any]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return [dict(row) for row in csv.DictReader(handle)]
    except (OSError, UnicodeError, csv.Error):
        return []


def _parse_nested(value):
    if isinstance(value, (dict, list)):
        return value
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return json.loads(text)
    except (json.JSONDecodeError, TypeError):
        try:
            return ast.literal_eval(text)
        except (ValueError, SyntaxError):
            return None


def _pep_patch_findings(reuse: list[dict], component_chains: list[str], chain_rows: dict[str, list[dict]],
                        *, interface_indices: dict[str, set[int]]) -> tuple[list[dict], list[dict]]:
    findings, summaries = [], []
    from mn_molecular_design_studio.complex_analysis_inspector import resolve_author_residue_selector
    for item in reuse:
        if item.get("evaluator") != "pep_patch":
            continue
        for artifact in item.get("artifacts") or []:
            name, path_text = str(artifact.get("name") or ""), str(artifact.get("path") or "")
            path = Path(path_text)
            if name == "patch_metrics.json":
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    continue
                summaries.append({"evaluation_id": item.get("evaluation_id"), "run_id": item.get("run_id"),
                                  "metrics": payload, "method": "PEP-Patch", "method_version": item.get("method_version") or "unrecorded"})
            elif name.startswith("residue_patches") and name.endswith(".csv"):
                for row in _parse_csv(path):
                    # The native artifact contains both surface methods. An
                    # electrostatic patch is not hydrophobic liability evidence.
                    if str(row.get("method") or "").strip().casefold() != "hydrophobic":
                        continue
                    nested = _parse_nested(row.get("per_residue_patches"))
                    if isinstance(nested, list):
                        patch_rows = nested
                    elif row.get("chain_id") and row.get("residue_id"):
                        patch_rows = [row]
                    else:
                        patch_rows = []
                    patch_config = _parse_nested(row.get("patch_config")) or {}
                    for patch in patch_rows:
                        patch = dict(patch)
                        chain_id = str(patch.get("chain_id") or row.get("chain_id") or "")
                        if chain_id not in component_chains:
                            continue
                        try:
                            residue_id = str(patch.get("residue_id", patch.get("author_number", ""))).strip()
                            identity = re.fullmatch(r"(-?\d+)([A-Za-z]?)", residue_id)
                            if identity is None:
                                continue
                            number = int(identity.group(1))
                        except (TypeError, ValueError):
                            continue
                        chain_track = chain_rows.get(chain_id) or []
                        candidates = [entry for entry in chain_track if entry.get("author_number") == number]
                        if identity.group(2):
                            candidates = [entry for entry in candidates
                                          if str(entry.get("insertion_code") or "") == identity.group(2)]
                        if len(candidates) != 1:
                            continue
                        residue = candidates[0]
                        patch_labels = _parse_nested(patch.get("patch_labels")) or {}
                        patch_id = patch.get("patch_id", patch.get("patch_label"))
                        patch_id_value = int(float(patch_id)) if str(patch_id or "").strip() else 0
                        if patch_id_value <= 0 and isinstance(patch_labels, dict):
                            patch_ids = [int(float(key)) for key in patch_labels if str(key).strip()]
                        else:
                            patch_ids = [patch_id_value] if patch_id_value > 0 else []
                        for patch_id_value in patch_ids:
                            area = patch_labels.get(str(patch_id_value), patch_labels.get(patch_id_value)) if isinstance(patch_labels, dict) else None
                            if area is None:
                                area = patch.get("area")
                            selector = resolve_author_residue_selector(
                                [{"chain_id": chain_id, "residues": [{"author_number": number,
                                  "insertion_code": residue.get("insertion_code") or "", "selector": residue.get("selector") or ""}]}],
                                chain_id, number, str(residue.get("insertion_code") or ""))
                            evidence = {
                                "patch_id": patch_id_value, "patch_configuration": patch_config,
                                "residue_patch_area_angstrom2": _number(area),
                                "residue_mean_hydrophobic_value": _number(patch.get("value")),
                                "saved_patch_metrics_evaluation_id": item.get("evaluation_id"),
                                "saved_patch_metrics_run_id": item.get("run_id"),
                                "source_evaluator": "PEP-Patch",
                                "limitation": "This is the existing PEP-Patch surface-patch assignment; it is not a universal aggregation score.",
                            }
                            finding = _finding("peppatch_hydrophobic_patch_member", "info", "structure_supported", residue,
                                "Residue belongs to a hydrophobic surface patch in a saved, matching PEP-Patch artifact.", evidence=evidence,
                                method="PEP-Patch saved result")
                            finding["selector"] = selector or finding["selector"]
                            finding["sequence_position"] = residue.get("sequence_index", 0) + 1
                            finding["interface_status"] = "in_interface" if finding["sequence_position"] - 1 in interface_indices.get(chain_id, set()) else "not_in_interface" if chain_id in interface_indices else "unknown"
                            findings.append(finding)
    return findings, summaries


def _number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _matching_saved_evidence(reused: list[dict], component_chains: list[str], partner_chains: list[str],
                             interface_ligand_selection: dict | None = None) -> tuple[set[tuple[str, int, str]], dict[str, Any], list[dict]]:
    from mn_molecular_design_studio.complex_analysis_inspector import parse_author_residue_selector
    interface = set()
    interface_provenance = {"status": "unavailable", "reason": "No exact, saved interface Evaluation was staged for this structure."}
    confidence = []
    for evaluation in reused:
        evaluator = str(evaluation.get("evaluator") or "")
        if evaluator in {"sequence_qc"}:
            continue
        for artifact in evaluation.get("artifacts") or []:
            if evaluator == "biotite_interface" and str(artifact.get("name") or "") == "contacts.json":
                try:
                    payload = json.loads(Path(str(artifact["path"])).read_text(encoding="utf-8"))
                except (OSError, KeyError, json.JSONDecodeError):
                    continue
                saved_a = list(map(str, payload.get("binder_chains") or []))
                saved_b = list(map(str, payload.get("target_chains") or []))
                if not partner_chains or not ((saved_a == component_chains and saved_b == partner_chains)
                                             or (saved_a == partner_chains and saved_b == component_chains)):
                    continue
                rows = payload.get("residue_pairs") or []
                for row in rows:
                    if not isinstance(row, dict):
                        continue
                    for field in ("binder_residue", "target_residue"):
                        parsed = parse_author_residue_selector(str(row.get(field) or ""))
                        if parsed and parsed["chain_id"] in component_chains:
                            interface.add((parsed["chain_id"], parsed["number"], str(parsed.get("insertion_code") or "")))
                interface_provenance = {
                    "status": "measured", "evaluator": "biotite_interface",
                    "evaluation_id": evaluation.get("evaluation_id"), "run_id": evaluation.get("run_id"),
                    "contact_distance_angstrom": payload.get("distance_angstrom"),
                    "component_chains": component_chains, "partner_chains": partner_chains,
                    "residue_pair_count": len(rows), "source_artifact": artifact.get("name"),
                    "method": "Saved Biotite Interface heavy-atom contact artifact",
                }
            if evaluator in {"plip", "pandamap"} and str(artifact.get("name") or "") == "interactions.csv":
                saved_ligand = evaluation.get("ligand_selection") if isinstance(evaluation.get("ligand_selection"), dict) else {}
                if not interface_ligand_selection or saved_ligand != interface_ligand_selection:
                    continue
                for row in _parse_csv(Path(str(artifact.get("path") or ""))):
                    chain_id = str(row.get("protein_chain") or row.get("chain_id") or "")
                    if chain_id not in component_chains:
                        continue
                    try:
                        number = int(row.get("protein_residue_number", row.get("residue_number", row.get("residue_no"))))
                    except (TypeError, ValueError):
                        continue
                    insertion = str(row.get("protein_insertion_code") or "")
                    interface.add((chain_id, number, insertion))
                if interface:
                    interface_provenance = {"status": "measured", "evaluator": evaluator,
                        "evaluation_id": evaluation.get("evaluation_id"), "run_id": evaluation.get("run_id"),
                        "ligand_selection": evaluation.get("ligand_selection") or {},
                        "method": f"Saved {evaluator} protein-residue interaction artifact"}
        metrics = evaluation.get("metrics") if isinstance(evaluation.get("metrics"), dict) else {}
        if any("plddt" in str(key).casefold() or "confidence" in str(key).casefold() for key in metrics):
            confidence.append({"evaluation_id": evaluation.get("evaluation_id"), "evaluator": evaluator,
                               "run_id": evaluation.get("run_id"), "metrics": metrics,
                               "scope": "saved Evaluation-wide values; not residue-local unless stated by that evaluator"})
    return interface, interface_provenance, confidence


def _findings_for_chain(sequence: str, rows: list[dict], chain_id: str, dssp: dict, interface_keys: set,
                        coverage: dict) -> list[dict]:
    exposure: dict[int, float | None] = {}
    for index, row in enumerate(rows):
        if not row.get("coordinates_present"):
            exposure[index] = None
            continue
        key = (chain_id, row.get("residue_id"))
        value = dssp.get(key) if isinstance(dssp, dict) else None
        exposure[index] = _number((value or {}).get("relative_accessibility"))
        row["relative_sasa"] = exposure[index]
        row["secondary_structure"] = (value or {}).get("secondary_structure")
    # Sequence-only evidence never jumps across an unknown or unmapped residue.
    normalized_chars = [char if char in _STANDARD else "X" for char in sequence]
    if coverage.get("status") == "observed_sequence_only":
        # Without a declared polymer sequence, do not create sequence motifs
        # across unresolved or non-peptide coordinate gaps.
        for index in range(1, len(rows)):
            previous, current = rows[index - 1], rows[index]
            first_residue, second_residue = previous.get("residue_object"), current.get("residue_object")
            if not previous.get("coordinates_present") or not current.get("coordinates_present"):
                normalized_chars[index] = "X"
                continue
            if first_residue is None or second_residue is None or "C" not in first_residue or "N" not in second_residue:
                normalized_chars[index] = "X"
                continue
            first_coord, second_coord = _coordinate(first_residue["C"]), _coordinate(second_residue["N"])
            if first_coord is None or second_coord is None or float(np.linalg.norm(first_coord - second_coord)) > 2.0:
                normalized_chars[index] = "X"
    normalized = "".join(normalized_chars)
    findings = _residue_findings(normalized, rows, exposure,
        {index for index, row in enumerate(rows) if (chain_id, int(row.get("author_number") or -999), str(row.get("insertion_code") or "")) in interface_keys},
        exposure_status="measured" if any(value is not None for value in exposure.values()) else "unavailable")
    for finding in findings:
        idx = int(finding.get("sequence_position") or 1) - 1
        if 0 <= idx < len(rows):
            finding["sequence_source"] = coverage.get("sequence_source")
            finding["coordinate_status"] = "present" if rows[idx].get("coordinates_present") else "unresolved_or_unmapped"
            finding["evidence"]["relative_sasa"] = exposure.get(idx)
            finding["evidence"]["secondary_structure"] = rows[idx].get("secondary_structure")
            finding["evidence"]["absolute_sasa_angstrom2"] = rows[idx].get("absolute_sasa_angstrom2")
            if rows[idx].get("declared_plddt") is not None:
                finding["evidence"]["declared_ca_plddt"] = rows[idx]["declared_plddt"]
                finding["evidence"]["confidence_method"] = "Cα B-factor explicitly declared as pLDDT"
            else:
                finding["evidence"].setdefault("confidence_status", "unavailable_or_b_factor_uninterpreted")
    return findings


def calculate_sequence_liabilities(
    structure_path: Path | str,
    *,
    component_chains: list[str],
    partner_chains: list[str],
    model_id: str,
    source_identity: dict[str, Any],
    reused_evidence: list[dict[str, Any]] | None = None,
    saved_confidence_context: dict[str, Any] | None = None,
    interface_ligand_selection: dict[str, Any] | None = None,
    component_format: str = "linear_protein",
    b_factor_interpretation: str = "uninterpreted",
) -> dict[str, Any]:
    """Assess explicit sequence motifs and coordinate-supported liabilities."""
    from mn_molecular_design_studio.complex_analysis_inspector import structure_digest
    from mn_molecular_design_studio.complex_analysis_tag_fusion import _dssp_data, _read_structure, _sasa_data

    path = Path(structure_path).expanduser().resolve(strict=True)
    source_digest = structure_digest(path)
    expected_digest = str(source_identity.get("structure_digest") or "")
    if not expected_digest or source_digest != expected_digest:
        raise ValueError("The selected structure digest does not match the authorized Database representation.")
    if str(model_id) != str(source_identity.get("model_id") or ""):
        raise ValueError("The selected model does not match the authorized Database representation.")
    component_chains = list(dict.fromkeys(str(value).strip() for value in component_chains if str(value).strip()))
    partner_chains = list(dict.fromkeys(str(value).strip() for value in partner_chains if str(value).strip()))
    if not component_chains or len(component_chains) != len(set(component_chains)):
        raise ValueError("Select at least one protein chain as the assessed component.")
    if set(component_chains) & set(partner_chains):
        raise ValueError("Assessed component and interface partner chain groups must be disjoint.")
    if component_format not in {"linear_protein", "cyclic_construct", "multichain_construct"}:
        raise ValueError("Choose a supported protein construct format.")
    if b_factor_interpretation not in {"uninterpreted", "plddt_0_100"}:
        raise ValueError("B-factors must remain uninterpreted unless the source explicitly confirms pLDDT values.")
    structure, model = _read_structure(path)
    selected_residues = _selected_residues(model, component_chains)
    for chain in partner_chains:
        if chain not in model:
            raise ValueError(f"Selected partner chain {chain!r} is absent from the selected model.")
    if not selected_residues:
        raise ValueError("The selected chains contain no supported amino-acid residues.")

    rows_by_chain, coverage_by_chain, sequence_properties = {}, {}, {}
    all_findings, all_residue_index = [], {}
    saved = list(reused_evidence or [])
    sequence_qc_by_chain = {str(item.get("chain_id") or ""): item for item in saved if item.get("evaluator") == "sequence_qc"}
    for chain_id in component_chains:
        rows, coverage = _track_for_chain(path, chain_id, model)
        if not rows:
            continue
        for index, row in enumerate(rows):
            row["sequence_index"] = index
            row["chain_id"] = chain_id
            row["residue_name"] = str(row.get("residue_name") or "")
            key = (chain_id, int(row["author_number"]), str(row.get("insertion_code") or "")) if row.get("author_number") is not None else None
            if key and key in selected_residues:
                row["residue_name"] = selected_residues[key]["residue_name"]
                row["amino_acid"] = selected_residues[key]["amino_acid"]
                row["residue_id"] = selected_residues[key]["residue_id"]
                row["residue_object"] = selected_residues[key]["residue_object"]
                row["coordinates_present"] = True
            elif row.get("coordinates_present"):
                row["coordinates_present"] = False
        rows_by_chain[chain_id] = rows
        coverage_by_chain[chain_id] = coverage
        sequence = "".join(str(row.get("amino_acid") or "X") for row in rows)
        saved_sequence = sequence_qc_by_chain.get(chain_id)
        sequence_properties[chain_id] = _safe_protein_metrics(sequence, saved_sequence)
        all_residue_index.update({(chain_id, row["author_number"], row.get("insertion_code") or ""): row
                                  for row in rows if row.get("author_number") is not None})

    if not rows_by_chain:
        raise ValueError("The selected protein chain(s) have no usable sequence track.")
    dssp, dssp_provenance = _dssp_data(path, model)
    sasa, sasa_provenance = _sasa_data(model)
    for chain_id, rows in rows_by_chain.items():
        for row in rows:
            residue_id = row.get("residue_id")
            if residue_id is None:
                continue
            row["absolute_sasa_angstrom2"] = _number(sasa.get((chain_id, residue_id)))
            dssp_row = dssp.get((chain_id, residue_id)) if isinstance(dssp, dict) else None
            row["relative_sasa"] = _number((dssp_row or {}).get("relative_accessibility"))
            row["secondary_structure"] = (dssp_row or {}).get("secondary_structure")
            if b_factor_interpretation == "plddt_0_100" and row.get("ca_b_factor") is not None:
                row["declared_plddt"] = row["ca_b_factor"]
            else:
                row["declared_plddt"] = None

    interface_keys, interface_provenance, confidence_evidence = _matching_saved_evidence(
        saved, component_chains, partner_chains, interface_ligand_selection=interface_ligand_selection,
    )
    declared_plddt = {
        chain_id: [{"author_number": row.get("author_number"), "insertion_code": row.get("insertion_code") or "",
                    "selector": row.get("selector") or "", "ca_plddt": row.get("declared_plddt"),
                    "interface_status": "in_interface" if (chain_id, int(row.get("author_number") or -999), str(row.get("insertion_code") or "")) in interface_keys
                                       else "not_in_interface" if interface_provenance.get("status") == "measured" else "unknown"}
                   for row in rows if row.get("declared_plddt") is not None]
        for chain_id, rows in rows_by_chain.items()
    }
    confidence_context = dict(saved_confidence_context or {})
    confidence_context_current = bool(
        confidence_context
        and str(confidence_context.get("structure_digest") or "") == expected_digest
        and str(confidence_context.get("model_id") or "") == str(model_id)
    )
    if confidence_context_current:
        metrics = confidence_context.get("native_scalar_metrics") if isinstance(confidence_context.get("native_scalar_metrics"), dict) else {}
        outputs = confidence_context.get("native_confidence_outputs") if isinstance(confidence_context.get("native_confidence_outputs"), dict) else {}
        selected_metrics = {key: value for key, value in metrics.items()
                            if any(token in str(key).casefold() for token in ("plddt", "confidence", "iptm", "ipae", "pae", "ptm"))}
        confidence_evidence.append({
            "status": "evaluation_context_only", "evaluation_id": confidence_context.get("evaluation_id"),
            "evaluator": confidence_context.get("evaluator"), "run_id": confidence_context.get("run_id"),
            "metrics": selected_metrics, "outputs": outputs,
            "scope": "Saved Evaluation-wide confidence context; not residue-local and not used to grade liabilities.",
            "structure_digest": expected_digest, "model_id": str(model_id),
        })
    elif confidence_context:
        confidence_evidence.append({"status": "not_reused_stale_or_mismatched",
                                    "reason": "Saved confidence context did not match the selected structure digest and model."})
    if any(declared_plddt.values()):
        confidence_evidence.append({
            "status": "residue_local_values", "basis": "Cα B-factors explicitly declared as pLDDT by the user",
            "method": "Source structure Cα B-factor field; no confidence engine was rerun or inferred.",
            "per_chain": declared_plddt,
            "scope": "Residue-local values are shown as stored; they are not converted into a liability grade.",
        })
    for chain_id, rows in rows_by_chain.items():
        sequence = "".join(str(row.get("amino_acid") or "X") for row in rows)
        chain_findings = _findings_for_chain(sequence, rows, chain_id, dssp, interface_keys, coverage_by_chain[chain_id])
        if component_format == "cyclic_construct":
            coverage_by_chain[chain_id]["terminal_assessment"] = {
                "status": "not_applicable", "reason": "The selected construct is declared cyclic; free N/C termini are not assumed."}
        elif component_format == "multichain_construct" or len(component_chains) > 1:
            coverage_by_chain[chain_id]["terminal_assessment"] = {
                "status": "limited", "reason": "The selected multichain component is reported chain by chain; biological chain connectivity and free termini are not inferred."}
        elif not coverage_by_chain[chain_id].get("complete_coordinate_coverage"):
            coverage_by_chain[chain_id]["terminal_assessment"] = {
                "status": "incomplete", "reason": "One or both declared construct termini or internal sequence positions are unresolved/unmapped; observed coordinate ends are not promoted to biological termini."}
        else:
            coverage_by_chain[chain_id]["terminal_assessment"] = {
                "status": "sequence_track_complete", "reason": "Sequence endpoints are present in the selected source track; this does not establish expression-construct processing."}
        if coverage_by_chain[chain_id].get("status") == "observed_sequence_only":
            coverage_by_chain[chain_id]["terminal_assessment"] = {
                "status": "unknown", "reason": "No declared polymer sequence is available; observed coordinate endpoints may be truncations and are not called biological termini."}
        all_findings.extend(chain_findings)

    disulfides = _disulfide_assessment(path, model, component_chains, all_residue_index)
    for chain_id, rows in rows_by_chain.items():
        for row in rows:
            if row.get("amino_acid") != "C" or row.get("coordinates_present") or row.get("author_number") is not None:
                continue
            # A declared polymer Cys can lack a unique author-number mapping
            # when unresolved. Preserve the sequence alert without inventing
            # an author selector or treating it as a geometrically free Cys.
            disulfides["findings"].append(_finding(
                "cysteine_sg_coordinates_missing", "info", "sequence_only", row,
                "Cysteine is present in the declared sequence but has no uniquely mapped coordinates; disulfide pairing is unknown.",
                evidence={"disulfide_basis": "unavailable", "coordinate_status": "unresolved_or_unmapped"},
            ))
            disulfides["findings"][-1]["sequence_index"] = row.get("sequence_index")
            disulfides["findings"][-1]["sequence_position"] = int(row.get("sequence_index") or 0) + 1
            if disulfides.get("status") == "not_present":
                disulfides["status"] = "unavailable"
                disulfides["reason"] = "A declared sequence cysteine is unresolved or unmapped; pairing status is unknown."
    all_findings.extend(disulfides["findings"])
    patch_findings, patch_summaries = _pep_patch_findings(saved, component_chains, rows_by_chain,
        interface_indices={chain: {idx for idx, row in enumerate(rows) if (chain, int(row.get("author_number") or -999), str(row.get("insertion_code") or "")) in interface_keys}
                           for chain, rows in rows_by_chain.items()})
    all_findings.extend(patch_findings)
    for finding in all_findings:
        chain_id = str(finding.get("chain_id") or "")
        residue_key = (chain_id, int(finding.get("author_number") or -999), str(finding.get("insertion_code") or ""))
        finding["interface_status"] = "in_interface" if residue_key in interface_keys else "not_in_interface" if interface_provenance.get("status") == "measured" else "unknown"
        if finding.get("selector"):
            finding["residue_ref"] = f"{chain_id}:{finding.get('author_residue')}"
        else:
            finding["residue_ref"] = f"{chain_id}:sequence#{int(finding.get('sequence_position') or 0)}"

    for index, finding in enumerate(all_findings):
        finding["finding_id"] = f"liability-{index + 1:05d}"

    missing_data = []
    for chain_id, properties in sequence_properties.items():
        if properties.get("status") == "unavailable":
            missing_data.append({"data": f"Sequence properties for chain {chain_id}", "status": "unavailable",
                                 "reason": properties.get("reason") or "No complete standard-amino-acid sequence was available for the existing Sequence Composition method."})
    for chain_id, coverage in coverage_by_chain.items():
        if not coverage.get("complete_coordinate_coverage"):
            missing_data.append({"data": f"Sequence-to-coordinate coverage for chain {chain_id}", "status": "partial",
                                 "reason": coverage.get("reason") or "Some declared sequence positions are unresolved or unmapped."})
    if disulfides.get("status") == "unavailable" and any(
        row.get("amino_acid") == "C" for rows in rows_by_chain.values() for row in rows
    ):
        missing_data.append({"data": "Cysteine disulfide pairing", "status": "unavailable",
                             "reason": disulfides.get("reason") or "No explicit bond records or usable sulfur geometry were available."})
    if dssp_provenance.get("status") != "measured":
        missing_data.append({"data": "DSSP relative accessibility / secondary structure",
                             "status": "unavailable", "reason": dssp_provenance.get("reason") or "DSSP returned no assignments."})
    if sasa_provenance.get("status") != "measured":
        missing_data.append({"data": "absolute residue SASA", "status": "unavailable",
                             "reason": sasa_provenance.get("reason") or "ShrakeRupley returned no assignments."})
    if interface_provenance.get("status") != "measured":
        missing_data.append({"data": "selected interface residue map", "status": "not_available",
                             "reason": "Run and save Biotite Interface with these exact chain groups, or save matching PLIP/PandaMap interaction analysis for this structure."})
    if not any(item.get("status") == "residue_local_values" for item in confidence_evidence):
        missing_data.append({"data": "residue-local structural confidence", "status": "not_available",
                             "reason": "No explicitly declared pLDDT-valued Cα B-factors or residue-local confidence Evaluation was available. Evaluation-wide confidence context, if shown, is not localized or used for liability grading."})
    if not patch_summaries and not patch_findings:
        missing_data.append({"data": "PEP-Patch hydrophobic surface patches", "status": "not_available",
                             "reason": "No exact-digest saved PEP-Patch Evaluation with residue-level patch artifacts was available; the evaluator did not duplicate PEP-Patch."})
    if any(not item.get("selector") for item in all_findings):
        missing_data.append({"data": "Mol* author-residue selection", "status": "partial",
                             "reason": "Some sequence findings map to unresolved or ambiguous residues and therefore have no coordinate selector."})
    if component_format == "cyclic_construct":
        missing_data.append({"data": "free termini", "status": "not_applicable", "reason": "The selected construct was explicitly declared cyclic."})
    if component_format == "multichain_construct" or len(component_chains) > 1:
        missing_data.append({"data": "construct-level chain connectivity and termini", "status": "unsupported", "reason": "The chains are screened separately; no covalent connectivity is inferred between them."})

    counts = {name: sum(finding.get("severity") == name for finding in all_findings) for name in ("warn", "info")}
    status = "partial" if missing_data else "measured"
    return {
        "schema": "studio.sequence_liabilities.v1", "method": METHOD, "method_version": METHOD_VERSION,
        "attribution": "Adapted from BinderScope binder_review.liabilities and binder_review.disulfides (MIT); see docs/third-party-notices.md.",
        "assessment_status": status,
        "source": {key: str(source_identity.get(key) or "") for key in
                   ("molecule_id", "structure_id", "group_id", "structure_digest", "model_id", "evaluation_id")},
        "interface_ligand_selection": dict(interface_ligand_selection or {}),
        "selected_component_chains": component_chains, "selected_partner_chains": partner_chains,
        "component_format": component_format,
        "parameters": {"ph": 7.4, "exposed_relative_sasa_threshold": EXPOSED_RSA_THRESHOLD,
                       "dibasic_relative_sasa_threshold": DIBASIC_RSA_THRESHOLD,
                       "hydrophobic_stretch_minimum_residues": HYDROPHOBIC_STRETCH_MINIMUM,
                       "b_factor_interpretation": b_factor_interpretation,
                       "interface_ligand_selection": dict(interface_ligand_selection or {}),
                       "disulfide_close_max_angstrom": DISULFIDE_CLOSE_MAX_ANGSTROM,
                       "disulfide_loose_max_angstrom": DISULFIDE_LOOSE_MAX_ANGSTROM,
                       "disulfide_matching_tie_tolerance_angstrom": DISULFIDE_MATCH_TIE_TOLERANCE_ANGSTROM},
        "method_provenance": {"sequence_properties": "Studio Sequence Composition / Biopython ProteinAnalysis",
                              "sasa": sasa_provenance, "dssp": dssp_provenance,
                              "interface": interface_provenance,
                              "reused_evaluations": [{key: item.get(key) for key in ("evaluator", "evaluation_id", "run_id", "source_digest", "chain_id", "partner_chains", "method_version")} for item in saved],
                              "pep_patch": [{"evaluation_id": item.get("evaluation_id"), "run_id": item.get("run_id"),
                                            "method": "PEP-Patch saved result"} for item in patch_summaries]},
        "sequence_properties": sequence_properties,
        "sequence_coverage": coverage_by_chain,
        "disulfides": disulfides,
        "interface": interface_provenance,
        "confidence_evidence": confidence_evidence,
        "hydrophobic_patch_summaries": patch_summaries,
        "findings": all_findings,
        "missing_data": missing_data,
        "measurement_counts": {"finding_count": len(all_findings), "warning_count": counts["warn"],
                               "information_count": counts["info"],
                               "unpaired_cysteine_count": sum(item.get("kind") == "unpaired_cysteine" for item in all_findings),
                               "interface_liability_count": sum(item.get("kind") in {"interface_lysine", "interface_methionine"} for item in all_findings),
                               "peppatch_patch_residue_count": len(patch_findings)},
        "source_structure_unchanged": True,
        "interpretation_limitations": [
            "Sequence motifs are alerts only; deamidation, isomerization, oxidation, glycosylation, cleavage and instability were not experimentally measured.",
            "A coordinate model does not establish redox state, conformational populations, expression, formulation behavior, or biological activity.",
            "SG distance-based disulfide assignment is model geometry when explicit file records are absent; ambiguous matching is reported without forcing a bond.",
            "Only a matching saved Biotite, PLIP, PandaMap, Sequence QC or PEP-Patch Evaluation is reused. No scientific tool output is regenerated inside this evaluator.",
            "Absolute SASA uses the existing Studio Bio.PDB ShrakeRupley implementation; relative accessibility and secondary structure use DSSP when available.",
            "Confidence is not converted into a developability grade. Saved Evaluation-wide confidence values are not treated as residue-local.",
        ],
    }
