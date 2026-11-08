"""Residue-resolved protein interface surface and polar-interaction evidence.

This is a Studio-native Python 3.11 implementation of a defined protein–protein
interface assessment.  It follows BinderScope's documented interface choices
where noted, while preserving a separately named/versioned Studio Evaluation.
It is not Rosetta InterfaceAnalyzer and does not infer binding affinity.
"""
from __future__ import annotations

import copy
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

METHOD = "Studio protein–protein buried interface surface and polar chemistry"
METHOD_VERSION = "studio.ppi-interface-chemistry.v1"
SASA_POINTS = 960
SASA_PROBE_ANGSTROM = 1.4
CONTACT_CUTOFF_ANGSTROM = 4.5
HBOND_DA_CUTOFF_ANGSTROM = 3.5
HBOND_HA_CUTOFF_ANGSTROM = 2.5
HBOND_DHA_MIN_DEGREES = 120.0
HBOND_ACCEPTOR_ANGLE_MIN_DEGREES = 90.0
BURIED_POLAR_MIN_DELTA_SASA = 0.5
BURIED_POLAR_MAX_COMPLEX_SASA = 0.1

_DONORS = {("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ"), ("ASN", "ND2"),
           ("GLN", "NE2"), ("HIS", "ND1"), ("HIS", "NE2"), ("TRP", "NE1"), ("SER", "OG"),
           ("THR", "OG1"), ("TYR", "OH")}
_ACCEPTORS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"), ("ASN", "OD1"),
              ("GLN", "OE1"), ("HIS", "ND1"), ("HIS", "NE2"), ("SER", "OG"), ("THR", "OG1"),
              ("TYR", "OH")}
_ACCEPTOR_ANTECEDENT = {
    ("ASP", "OD1"): ("CG",), ("ASP", "OD2"): ("CG",),
    ("GLU", "OE1"): ("CD",), ("GLU", "OE2"): ("CD",),
    ("ASN", "OD1"): ("CG",), ("GLN", "OE1"): ("CD",),
    ("HIS", "ND1"): ("CG", "CE1"), ("HIS", "NE2"): ("CD2", "CE1"),
    ("SER", "OG"): ("CB",), ("THR", "OG1"): ("CB",), ("TYR", "OH"): ("CZ",),
}
_AMINO_ACIDS = frozenset({"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                          "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                          "MSE", "SEC", "PYL"})
_AA3_TO_1 = {"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C", "GLN":"Q", "GLU":"E",
             "GLY":"G", "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F",
             "PRO":"P", "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V", "MSE":"M",
             "SEC":"U", "PYL":"O"}


def _element(atom) -> str:
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if element:
        return element
    name = str(atom.get_name() or "").strip().upper()
    return "H" if name.startswith(("H", "D")) else (name[:1] or "X")


def _key(chain_id: str, residue, atom) -> tuple[str, int, str, str]:
    return (str(chain_id), int(residue.id[1]), str(residue.id[2] or "").strip(), str(atom.get_name()).strip())


def _copy_protein_model(model, chains: set[str], *, heavy_only: bool):
    from Bio.PDB.Polypeptide import is_aa

    selected = copy.deepcopy(model)
    for chain in list(selected.get_chains()):
        if str(chain.id) not in chains:
            selected.detach_child(chain.id)
            continue
        for residue in list(chain.get_residues()):
            if not is_aa(residue, standard=False):
                chain.detach_child(residue.id)
                continue
            for atom in list(residue.get_atoms()):
                if heavy_only and _element(atom) in {"H", "D"}:
                    residue.detach_child(atom.id)
    if not list(selected.get_atoms()):
        raise ValueError("The selected protein chains contain no atoms for interface assessment.")
    return selected


def _atom_sasa(model, chains: set[str]) -> dict[tuple[str, int, str, str], float]:
    from Bio.PDB.SASA import ShrakeRupley
    from mn_molecular_design_studio.complex_analysis_tag_fusion import _SASA_COMPAT_LOCK

    selected = _copy_protein_model(model, chains, heavy_only=True)
    # Biopython 1.79 references the removed np.int alias; the shared lock keeps
    # this narrowly scoped compatibility shim safe alongside Tag/Fusion SASA.
    with _SASA_COMPAT_LOCK:
        added_alias = not hasattr(np, "int")
        if added_alias:
            setattr(np, "int", int)
        try:
            ShrakeRupley(n_points=SASA_POINTS).compute(selected, level="A")
        finally:
            if added_alias:
                delattr(np, "int")
    result = {}
    for chain in selected:
        for residue in chain:
            for atom in residue.get_atoms():
                value = getattr(atom, "sasa", None)
                if value is not None and math.isfinite(float(value)):
                    result[_key(chain.id, residue, atom)] = float(value)
    return result


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float | None:
    u, v = a - b, c - b
    denominator = float(np.linalg.norm(u) * np.linalg.norm(v))
    if denominator < 1e-8:
        return None
    return math.degrees(math.acos(float(np.clip(np.dot(u, v) / denominator, -1.0, 1.0))))


def _polar_atom(atom, residue) -> tuple[bool, bool]:
    name, residue_name, element = str(atom.get_name()).strip().upper(), str(residue.resname).strip().upper(), _element(atom)
    if element in {"H", "D"}:
        return False, False
    donor = (name == "N" and residue_name != "PRO") or (residue_name, name) in _DONORS
    acceptor = name in {"O", "OXT"} or (residue_name, name) in _ACCEPTORS
    return donor, acceptor


def _hbond_evidence(model, binder_chains: set[str], target_chains: set[str]):
    """BinderScope-style geometry: explicit H, placed backbone amide H, or acceptor-angle proxy."""
    atoms = []
    donor_rows, acceptor_rows = [], []
    for chain in model:
        if str(chain.id) not in binder_chains | target_chains:
            continue
        residues = [residue for residue in chain if str(residue.resname).strip().upper() in _AMINO_ACIDS]
        for index, residue in enumerate(residues):
            residue_atoms = list(residue.get_atoms())
            heavy_atoms = [atom for atom in residue_atoms if _element(atom) not in {"H", "D"}]
            for atom in heavy_atoms:
                xyz = np.asarray(atom.get_coord(), dtype=float)
                if xyz.shape != (3,) or not np.isfinite(xyz).all():
                    continue
                row = {"atom": atom, "residue": residue, "chain": str(chain.id), "xyz": xyz,
                       "key": _key(chain.id, residue, atom), "index": index}
                atoms.append(row)
                donor, acceptor = _polar_atom(atom, residue)
                if donor:
                    explicit = [h for h in residue_atoms if _element(h) in {"H", "D"}
                                and float(np.linalg.norm(np.asarray(h.get_coord(), dtype=float) - xyz)) <= 1.25]
                    if explicit:
                        row["hydrogens"] = [np.asarray(h.get_coord(), dtype=float) for h in explicit]
                        row["h_basis"] = "explicit"
                    elif str(atom.get_name()).strip().upper() == "N" and str(residue.resname).upper() != "PRO" and index > 0:
                        previous = residues[index - 1]
                        if all(name in previous for name in ("C", "O")) and float(np.linalg.norm(previous["C"].coord - xyz)) <= 2.0:
                            vector = np.asarray(previous["C"].coord - previous["O"].coord, dtype=float)
                            norm = float(np.linalg.norm(vector))
                            if norm > 1e-8:
                                row["hydrogens"] = [xyz + vector / norm]
                                row["h_basis"] = "computed_backbone_amide"
                    donor_rows.append(row)
                if acceptor:
                    antecedents = ("C",) if str(atom.get_name()).strip().upper() in {"O", "OXT"} else _ACCEPTOR_ANTECEDENT.get((str(residue.resname).upper(), str(atom.get_name()).strip().upper()), ())
                    row["antecedent_xyz"] = [np.asarray(residue[name].coord, dtype=float) for name in antecedents if name in residue]
                    acceptor_rows.append(row)

    if not donor_rows or not acceptor_rows:
        return [], set(), {"explicit": 0, "computed_backbone_amide": 0, "heavy_atom_angle_proxy": 0}
    tree = cKDTree(np.asarray([row["xyz"] for row in acceptor_rows]))
    donor_xyz = np.asarray([row["xyz"] for row in donor_rows])
    evidence, satisfied, basis_counts = [], set(), {"explicit": 0, "computed_backbone_amide": 0, "heavy_atom_angle_proxy": 0}
    for donor_index, acceptor_indexes in enumerate(tree.query_ball_point(donor_xyz, HBOND_DA_CUTOFF_ANGSTROM)):
        donor = donor_rows[donor_index]
        for acceptor_index in acceptor_indexes:
            acceptor = acceptor_rows[acceptor_index]
            if donor["key"] == acceptor["key"] or donor["chain"] == acceptor["chain"] and donor["residue"].id == acceptor["residue"].id:
                continue
            if (donor["chain"] == acceptor["chain"] and donor["atom"].get_name() == "N"
                    and acceptor["atom"].get_name() == "O" and donor["index"] == acceptor["index"] + 1):
                continue
            distance = float(np.linalg.norm(donor["xyz"] - acceptor["xyz"]))
            h_distance, dha_angle, basis = None, None, "heavy_atom_angle_proxy"
            donor_h = donor.get("hydrogens") or []
            if donor_h:
                candidates = []
                for hydrogen in donor_h:
                    hdist = float(np.linalg.norm(hydrogen - acceptor["xyz"]))
                    angle = _angle(donor["xyz"], hydrogen, acceptor["xyz"])
                    if hdist <= HBOND_HA_CUTOFF_ANGSTROM and angle is not None and angle >= HBOND_DHA_MIN_DEGREES:
                        candidates.append((hdist, angle))
                if not candidates:
                    continue
                h_distance, dha_angle = min(candidates)
                basis = donor.get("h_basis", "explicit")
            else:
                antecedents = acceptor.get("antecedent_xyz") or []
                if antecedents:
                    antecedent = np.mean(np.asarray(antecedents), axis=0)
                    angle = _angle(antecedent, acceptor["xyz"], donor["xyz"])
                    if angle is not None and angle < HBOND_ACCEPTOR_ANGLE_MIN_DEGREES:
                        continue
            satisfied.add(donor["key"])
            satisfied.add(acceptor["key"])
            is_interfacial = ((donor["chain"] in binder_chains and acceptor["chain"] in target_chains)
                              or (donor["chain"] in target_chains and acceptor["chain"] in binder_chains))
            if is_interfacial:
                basis_counts[basis] += 1
            if is_interfacial:
                evidence.append({"donor": {"chain_id": donor["chain"], "author_number": int(donor["residue"].id[1]),
                                    "insertion_code": str(donor["residue"].id[2] or "").strip(),
                                    "residue_name": str(donor["residue"].resname), "atom": str(donor["atom"].get_name())},
                                 "acceptor": {"chain_id": acceptor["chain"], "author_number": int(acceptor["residue"].id[1]),
                                    "insertion_code": str(acceptor["residue"].id[2] or "").strip(),
                                    "residue_name": str(acceptor["residue"].resname), "atom": str(acceptor["atom"].get_name())},
                                 "donor_acceptor_distance_angstrom": round(distance, 3),
                                 "hydrogen_acceptor_distance_angstrom": round(h_distance, 3) if h_distance is not None else None,
                                 "donor_hydrogen_acceptor_angle_degrees": round(dha_angle, 2) if dha_angle is not None else None,
                                 "evidence_basis": basis})
    return evidence, satisfied, basis_counts


def _residue_totals(atom_sasa: dict, chain_ids: set[str]) -> dict[tuple[str, int, str], float]:
    result: dict[tuple[str, int, str], float] = {}
    for (chain, number, insertion, _atom), area in atom_sasa.items():
        if chain in chain_ids:
            key = (chain, number, insertion)
            result[key] = result.get(key, 0.0) + float(area)
    return result


def assess_protein_interface(
    structure_path: Path | str,
    *,
    binder_chains: list[str],
    target_chains: list[str],
    source_identity: dict[str, Any],
) -> dict[str, Any]:
    """Measure one exact, authorized single-model Database protein–protein representation."""
    from Bio.PDB.Polypeptide import is_aa
    from mn_molecular_design_studio.complex_analysis_comparison import _read_biopython_structure
    from mn_molecular_design_studio.complex_analysis_inspector import author_residue_selector, structure_digest

    path = Path(structure_path).expanduser().resolve(strict=True)
    digest = structure_digest(path)
    if not source_identity.get("structure_digest") or digest != str(source_identity["structure_digest"]):
        raise ValueError("The selected structure digest does not match the authorized Database representation.")
    if not source_identity.get("model_id"):
        raise ValueError("The authorized Database representation has no exact model identity.")
    binder = list(dict.fromkeys(str(value).strip() for value in binder_chains if str(value).strip()))
    target = list(dict.fromkeys(str(value).strip() for value in target_chains if str(value).strip()))
    if not binder or not target or set(binder) & set(target):
        raise ValueError("Choose explicit, nonempty and disjoint binder and target protein chain groups.")
    if len(binder) != len(binder_chains) or len(target) != len(target_chains):
        raise ValueError("Duplicate or empty chain identifiers are not supported.")
    structure = _read_biopython_structure(path)
    models = list(structure.get_models())
    if len(models) != 1:
        raise ValueError("Interface chemistry requires one exact single-model Database representation.")
    model = models[0]
    selected = set(binder + target)
    for chain_id in selected:
        if chain_id not in model:
            raise ValueError(f"Selected protein chain {chain_id!r} is absent from the exact model.")
        if not any(is_aa(residue, standard=False) for residue in model[chain_id]):
            raise ValueError(f"Selected chain {chain_id!r} contains no supported amino-acid residues.")

    sasa_complex = _atom_sasa(model, selected)
    sasa_binder = _atom_sasa(model, set(binder))
    sasa_target = _atom_sasa(model, set(target))
    target_atoms = [(chain, residue, atom) for chain in target for residue in model[chain]
                    if is_aa(residue, standard=False) for atom in residue.get_atoms() if _element(atom) not in {"H", "D"}]
    binder_atoms = [(chain, residue, atom) for chain in binder for residue in model[chain]
                    if is_aa(residue, standard=False) for atom in residue.get_atoms() if _element(atom) not in {"H", "D"}]
    if not target_atoms or not binder_atoms:
        raise ValueError("Selected chain groups must both contain protein heavy atoms.")

    txyz = np.asarray([atom.coord for _chain, _residue, atom in target_atoms], dtype=float)
    bxyz = np.asarray([atom.coord for _chain, _residue, atom in binder_atoms], dtype=float)
    contacts: dict[tuple[tuple[str, int, str], tuple[str, int, str]], dict[str, Any]] = {}
    atom_contact_counts: dict[tuple[str, int, str, str], int] = {}
    atom_tree = cKDTree(txyz)
    for binder_index, target_indexes in enumerate(atom_tree.query_ball_point(bxyz, CONTACT_CUTOFF_ANGSTROM)):
        bchain, bresidue, batom = binder_atoms[binder_index]
        bres = (str(bchain), int(bresidue.id[1]), str(bresidue.id[2] or "").strip())
        bkey = _key(bchain, bresidue, batom)
        for target_index in target_indexes:
            tchain, tresidue, tatom = target_atoms[target_index]
            tres = (str(tchain), int(tresidue.id[1]), str(tresidue.id[2] or "").strip())
            tkey = _key(tchain, tresidue, tatom)
            distance = float(np.linalg.norm(bxyz[binder_index] - txyz[target_index]))
            pair = (tres, bres)
            row = contacts.setdefault(pair, {"target": tres, "binder": bres, "minimum_distance_angstrom": distance,
                                             "heavy_atom_contact_count": 0})
            row["minimum_distance_angstrom"] = min(row["minimum_distance_angstrom"], distance)
            row["heavy_atom_contact_count"] += 1
            atom_contact_counts[tkey] = atom_contact_counts.get(tkey, 0) + 1
            atom_contact_counts[bkey] = atom_contact_counts.get(bkey, 0) + 1

    hbonds, satisfied, hbond_basis = _hbond_evidence(model, set(binder), set(target))
    atom_rows = []
    for chain, residue, atom in target_atoms + binder_atoms:
        chain_id, atom_key = str(chain), _key(chain, residue, atom)
        alone = (sasa_binder if chain_id in binder else sasa_target).get(atom_key, 0.0)
        bound = sasa_complex.get(atom_key, 0.0)
        buried = float(alone - bound)
        donor, acceptor = _polar_atom(atom, residue)
        atom_rows.append({"key": atom_key, "chain_id": chain_id, "residue": residue, "atom": atom,
                          "isolated_sasa": float(alone), "complex_sasa": float(bound), "buried_sasa": buried,
                          "polar": bool(donor or acceptor), "hydrophobic_element": _element(atom) in {"C", "S"},
                          "contact_count": int(atom_contact_counts.get(atom_key, 0)),
                          "satisfied": atom_key in satisfied})
    total_buried = float(sum(row["buried_sasa"] for row in atom_rows))
    hydrophobic_buried = float(sum(row["buried_sasa"] for row in atom_rows if row["hydrophobic_element"]))
    unsatisfied = [row for row in atom_rows if row["polar"] and row["buried_sasa"] >= BURIED_POLAR_MIN_DELTA_SASA
                   and row["complex_sasa"] <= BURIED_POLAR_MAX_COMPLEX_SASA and not row["satisfied"]]
    t_sasa = _residue_totals(sasa_target, set(target))
    b_sasa = _residue_totals(sasa_binder, set(binder))
    t_bound = _residue_totals(sasa_complex, set(target))
    b_bound = _residue_totals(sasa_complex, set(binder))
    residue_rows = []
    for side, chains, isolated, bound in (("target", target, t_sasa, t_bound), ("binder", binder, b_sasa, b_bound)):
        residue_keys = set(isolated) | set(bound)
        for chain_id, number, insertion in sorted(residue_keys):
            residue = next((residue for residue in model[chain_id] if int(residue.id[1]) == number
                            and str(residue.id[2] or "").strip() == insertion and is_aa(residue, standard=False)), None)
            if residue is None:
                continue
            key = (chain_id, number, insertion)
            buried = float(isolated.get(key, 0.0) - bound.get(key, 0.0))
            partner_keys = sorted({pair[1] if side == "target" else pair[0] for pair in contacts
                                   if (pair[0] if side == "target" else pair[1]) == key})
            atom_subset = [row for row in atom_rows if row["chain_id"] == chain_id
                           and int(row["residue"].id[1]) == number and str(row["residue"].id[2] or "").strip() == insertion]
            buns = [row for row in unsatisfied if row["key"] in {atom["key"] for atom in atom_subset}]
            relevant_hbonds = [row for row in hbonds if any(row[role]["chain_id"] == chain_id
                               and int(row[role]["author_number"]) == number
                               and str(row[role].get("insertion_code") or "") == insertion
                               for role in ("donor", "acceptor"))]
            residue_rows.append({
                "side": side, "chain_id": chain_id, "author_number": number, "insertion_code": insertion,
                "author_residue": f"{number}{insertion}", "selector": author_residue_selector(chain_id, number, insertion),
                "residue_name": str(residue.resname).strip().upper(),
                "amino_acid": _AA3_TO_1.get(str(residue.resname).strip().upper(), "X"),
                "isolated_sasa_angstrom2": round(float(isolated.get(key, 0.0)), 3),
                "complex_sasa_angstrom2": round(float(bound.get(key, 0.0)), 3),
                "buried_sasa_angstrom2": round(buried, 3),
                "buried_fraction_of_isolated_sasa": round(buried / isolated[key], 5) if isolated.get(key, 0.0) > 0 else None,
                "heavy_atom_contact_count": sum(int(row["contact_count"]) for row in atom_subset),
                "partner_residues": [f"{c}:{n}{i}" for c, n, i in partner_keys],
                "polar_atoms": sum(bool(row["polar"]) for row in atom_subset),
                "buried_unpaired_polar_atoms": [{"atom": str(row["atom"].get_name()),
                    "buried_sasa_angstrom2": round(row["buried_sasa"], 3),
                    "complex_sasa_angstrom2": round(row["complex_sasa"], 3)} for row in buns],
                "hydrogen_bond_evidence_count": len(relevant_hbonds),
            })

    source = {str(key): value for key, value in source_identity.items()}
    binder_delta = float(sum(row["buried_sasa"] for row in atom_rows if row["chain_id"] in binder))
    target_delta = float(sum(row["buried_sasa"] for row in atom_rows if row["chain_id"] in target))
    return {
        "schema": "studio.ppi-interface-chemistry.v1", "method": METHOD, "method_version": METHOD_VERSION,
        "source": source, "selected_chains": {"binder": binder, "target": target},
        "parameters": {"sasa_algorithm": "Bio.PDB ShrakeRupley", "sasa_sphere_points": SASA_POINTS,
            "sasa_probe_angstrom": SASA_PROBE_ANGSTROM, "heavy_atoms_only_for_sasa": True,
            "contact_cutoff_angstrom": CONTACT_CUTOFF_ANGSTROM,
            "hbond_donor_acceptor_cutoff_angstrom": HBOND_DA_CUTOFF_ANGSTROM,
            "hbond_hydrogen_acceptor_cutoff_angstrom": HBOND_HA_CUTOFF_ANGSTROM,
            "hbond_min_donor_hydrogen_acceptor_angle_degrees": HBOND_DHA_MIN_DEGREES,
            "hbond_min_acceptor_antecedent_angle_degrees_without_hydrogen": HBOND_ACCEPTOR_ANGLE_MIN_DEGREES,
            "buried_polar_minimum_delta_sasa_angstrom2": BURIED_POLAR_MIN_DELTA_SASA,
            "buried_polar_maximum_complex_sasa_angstrom2": BURIED_POLAR_MAX_COMPLEX_SASA,
            "heteroatoms_waters_and_ligands": "excluded from this protein–protein interface calculation"},
        "metrics": {"interface_total_delta_sasa_angstrom2": round(total_buried, 3),
            "interface_buried_surface_area_half_delta_sasa_angstrom2": round(total_buried / 2.0, 3),
            "binder_delta_sasa_angstrom2": round(binder_delta, 3),
            "target_delta_sasa_angstrom2": round(target_delta, 3),
            "hydrophobic_cs_buried_area_angstrom2": round(hydrophobic_buried, 3),
            "hydrophobic_cs_buried_fraction": round(hydrophobic_buried / total_buried, 6) if abs(total_buried) > 1e-9 else None,
            "buried_unpaired_polar_atom_count": len(unsatisfied),
            "interfacial_hydrogen_bond_count": len(hbonds),
            "interfacial_hydrogen_bond_evidence_basis_counts": hbond_basis,
            "contacting_heavy_atom_pair_count": sum(row["heavy_atom_contact_count"] for row in contacts.values()),
            "contacting_residue_pair_count": len(contacts),
            "target_interface_residue_count": len({pair[0] for pair in contacts}),
            "binder_interface_residue_count": len({pair[1] for pair in contacts}),
            "buried_polar_assessment_status": "measured"},
        "contacts": [{**row, "minimum_distance_angstrom": round(row["minimum_distance_angstrom"], 3)}
                     for _pair, row in sorted(contacts.items())],
        "hydrogen_bonds": hbonds, "residue_annotations": residue_rows,
        "methodological_assumptions": [
            "Buried surface uses solvent-accessible heavy-atom area lost on binding; the reported half-delta value is half of the two-side total. BinderScope's dSASA convention reports the two-side total.",
            "Hydrophobic buried area fraction uses buried carbon/sulfur atom SASA divided by total two-side delta SASA, matching BinderScope's hyd_frac definition. It is distinct from BinderScope's Interface_Hydrophobicity metric, the hydrophobic-residue percentage of its 4.0 Å interface-residue set. Numerical identity depends on the same coordinate representation and SASA implementation.",
            "Buried unpaired polar atoms follow the BinderScope geometric rule: N/O donor or acceptor heavy atom, burial >=0.5 Å², complex SASA <=0.1 Å², and no qualifying intracomplex hydrogen bond. This is geometry-based; no protonation, solvent network, or energetic penalty is modeled.",
            "Hydrogen-bond evidence uses explicit H atoms when present; otherwise backbone amide H positions are constructed from preceding C=O geometry, while side-chain donors without H use donor–acceptor distance and acceptor-angle geometry. Counts retain the evidence basis; they are not energies.",
            "SASA is Bio.PDB Shrake–Rupley with 960 sphere points and 1.4 Å probe. Only selected protein chains and heavy atoms contribute; waters, ligands, glycans, alternate protonation states, and unselected components are excluded.",
            "Contact residues use any selected binder–target heavy-atom distance <=4.5 Å. Residue identifiers retain chain, author numbering, and insertion code.",
        ],
        "attribution": {"reference": "BinderScope binder_review.interface.py / binder_review.structure.py", "license": "MIT",
                        "relationship": "Studio-native implementation; Bio.PDB Shrake-Rupley and documented BinderScope cutoffs/definitions, not an imported BindCraft-specific pipeline."},
    }
