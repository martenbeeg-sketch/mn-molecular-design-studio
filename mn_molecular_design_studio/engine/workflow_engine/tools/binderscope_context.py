"""Exact residue mapping and reference comparisons for the offline adapter."""
from __future__ import annotations

import re


def prepare_context(context, reference, mapping, sequences, canonical, reader):
    import gemmi
    import numpy as np

    declaration = context.get("hotspots") or {"status": "not_recorded", "requested": []}
    lookup, reference_rows = {}, []
    reference_status = {"status": "unavailable", "reason": (context.get("reference") or {}).get("reason") or "No verified original input target."}
    if reference:
        original = reader(reference["path"])
        if len(original) != 1 or [chain.name for chain in original[0]].count(reference["chain"]) != 1:
            raise ValueError("Original target reference must identify one exact model and chain.")
        chain = original[0][reference["chain"]]
        residues = [residue for residue in chain if gemmi.find_tabulated_residue(residue.name).is_amino_acid()
                    and any(atom.name == "CA" for atom in residue)]
        sequence = "".join(gemmi.find_tabulated_residue(residue.name).one_letter_code.upper() for residue in residues)
        positions = [index for index in range(len(sequence) - len(sequences["A"]) + 1)
                     if sequence[index:index + len(sequences["A"])] == sequences["A"]]
        if len(positions) == 1:
            offset = positions[0]
            for index, author in enumerate(mapping["A"]):
                residue = residues[offset + index]
                key = (reference["chain"], residue.seqid.num, residue.seqid.icode.strip())
                lookup.setdefault(key, []).append(author["canonical_number"])
                reference_rows.append(next(atom for atom in residue if atom.name == "CA").pos)
            reference_status = {"status": "assessed", "reference": context["reference"],
                "mapped_residues": len(mapping["A"]), "reference_residues": len(residues),
                "scope": "whole target" if len(residues) == len(mapping["A"]) else "unique contiguous target segment",
                "interpretation": "Coordinate conformation comparison, not experimental proof of fold preservation."}
        else:
            reference_status = {"status": "unavailable", "reason": "Target sequence does not map to one unique contiguous segment of the original reference."}
    if declaration.get("numbering") == "saved_target_author":
        lookup = {}
        for author in mapping["A"]:
            lookup.setdefault((author["chain"], author["number"], author["insertion_code"]), []).append(author["canonical_number"])
    rows = []
    for token in declaration.get("requested") or []:
        row = {"requested": token, "status": "unmapped"}
        match = re.fullmatch(r"(?:(?P<chain>[^:,\s]+):)?(?P<number>-?\d+)(?P<insertion>[A-Za-z]?)", token)
        if match is None:
            match = re.fullmatch(r"(?P<chain>[A-Za-z]+)(?P<number>-?\d+)(?P<insertion>[A-Za-z]?)", token)
        if declaration.get("status") == "conflicting":
            row["reason"] = "Conflicting original hotspot declarations."
        elif match is None:
            row["reason"] = "Unsupported hotspot syntax; chain and residue identity were not guessed."
        else:
            chain = match.group("chain")
            candidates = [number for (original_chain, original_number, insertion), numbers in lookup.items()
                          if (not chain or original_chain == chain) and original_number == int(match.group("number"))
                          and insertion == match.group("insertion") for number in numbers]
            if len(candidates) == 1:
                number = candidates[0]
                author = next(item for item in mapping["A"] if item["canonical_number"] == number)
                row.update(status="mapped", canonical_number=number,
                           saved_residue=f"{author['chain']}:{author['number']}{author['insertion_code']}")
            else:
                row["reason"] = "Residue is absent, ambiguous or cannot be mapped from the original target numbering."
        rows.append(row)
    if context.get("target_kind") == "sequence":
        reference_status = {"status": "not_applicable", "reason": "Sequence-only or disordered target: no fixed input fold is assumed."}
    elif reference_rows:
        selected = reader(canonical)
        mobile = [next(atom for atom in residue if atom.name == "CA").pos for residue in selected[0]["A"]]
        fixed_xyz = np.array([[point.x, point.y, point.z] for point in reference_rows])
        mobile_xyz = np.array([[point.x, point.y, point.z] for point in mobile])
        fixed_xyz -= fixed_xyz.mean(axis=0)
        mobile_xyz -= mobile_xyz.mean(axis=0)
        if len(mobile) >= 3 and min(np.linalg.matrix_rank(fixed_xyz), np.linalg.matrix_rank(mobile_xyz)) >= 2:
            left, singular, right = np.linalg.svd(mobile_xyz.T @ fixed_xyz)
            correction = np.diag([1, 1, np.linalg.det(left @ right)])
            difference = mobile_xyz @ (left @ correction @ right) - fixed_xyz
            reference_status["target_ca_rmsd_angstrom"] = float(np.sqrt(np.mean(np.sum(difference ** 2, axis=1))))
        else:
            reference_status.update(status="unavailable", reason="Insufficient non-collinear Cα coordinates for a rigid-body comparison.")
    return {"declaration": declaration, "rows": rows, "target_reference_comparison": reference_status,
            "interpretation": "Contacts in this saved complex; for flexible targets they are conformation-specific. Admission decisions are unchanged."}


def hotspot_result(prepared, analysis):
    from binder_review.model import HOTSPOT_TOUCH_CUTOFF

    rows = []
    for original in prepared["rows"]:
        row = dict(original)
        if row["status"] == "mapped":
            number = row["canonical_number"]
            distance = analysis.hotspot_dist.get(number)
            if distance is None:
                row.update(status="unmapped", reason="No hotspot distance was returned for this residue.")
            else:
                row.update(status="contacted" if number in analysis.hot else "missed", minimum_distance_angstrom=distance)
        rows.append(row)
    requested = len(rows)
    contacted = sum(row["status"] == "contacted" for row in rows)
    missed = sum(row["status"] == "missed" for row in rows)
    unmapped = sum(row["status"] == "unmapped" for row in rows)
    declaration_status = prepared["declaration"]["status"]
    status = declaration_status if declaration_status != "recorded" else "partial" if unmapped else "assessed"
    return {"status": status, "requested_count": requested, "contacted_count": contacted, "missed_count": missed,
            "unmapped_count": unmapped, "coverage_fraction": contacted / requested if requested and status == "assessed" else None,
            "cutoff_angstrom": HOTSPOT_TOUCH_CUTOFF, "distance_rule": "Minimum heavy-atom distance; BinderScope residue contacts.",
            "residues": rows, "evidence": prepared["declaration"].get("evidence") or [], "interpretation": prepared["interpretation"]}
