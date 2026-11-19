"""Read-only preparation of the linked, locally bundled design inspector."""
from __future__ import annotations

from functools import lru_cache
import hashlib
import html
import json
import math
from pathlib import Path
import re

from mn_molecular_design_studio.binderscope_presentation import contact_rows, residue_label, result_coordinates


ASSETS = Path(__file__).parent / "components/binderscope_inspector"
CHARGED_ATOMS = {"ARG": ("NE", "NH1", "NH2", "CZ"), "LYS": ("NZ",), "ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}


def _structure(text):
    import gemmi

    if "_atom_site." in text:
        document = gemmi.cif.read_string(text)
        blocks = [block for block in document if len(block.find_values("_atom_site.id"))]
        if len(blocks) != 1:
            raise ValueError("One exact coordinate-bearing block is required for inspection.")
        category = blocks[0].get_mmcif_category("_atom_site.")
        category.setdefault("occupancy", ["1.0"] * len(category["id"]))
        category.setdefault("B_iso_or_equiv", ["0.0"] * len(category["id"]))
        blocks[0].set_mmcif_category("_atom_site.", category)
        structure = gemmi.make_structure_from_block(blocks[0])
    else:
        structure = gemmi.read_pdb_string(text)
    if len(structure) != 1:
        raise ValueError("The assessment must identify one saved model.")
    structure.remove_alternative_conformations()
    return structure


def _finite(value):
    return float(value) if type(value) in {int, float} and math.isfinite(value) else None


def build_view_payload(result, coordinates, *, title=""):
    import gemmi

    expected = result.get("source", {}).get("input_structure_digest")
    if not expected or hashlib.sha256(coordinates.encode()).hexdigest() != expected:
        raise ValueError("The inspection coordinates do not match the assessed source checksum.")
    try:
        structure = _structure(coordinates)
    except RuntimeError as exc:
        raise ValueError("Saved coordinates could not be parsed for the linked inspector.") from exc
    mapping = result.get("source", {}).get("residue_map") or {}
    if not all(mapping.get(side) for side in ("A", "B")):
        raise ValueError("Exact target and binder residue mappings are required.")
    analysis = result.get("analysis") or {}
    confidence = analysis.get("plddt_res") or []
    confidence_confirmed = (result.get("source", {}).get("assessment", {}).get("b_factor_interpretation") == "plddt_0_100"
        and len(confidence) == len(mapping["B"]) and all(_finite(value) is not None and 0 <= value <= 100 for value in confidence))
    display = gemmi.Structure()
    display.add_model(gemmi.Model("1"))
    residues, atoms, selected_chains = {}, {}, set()
    for side in ("A", "B"):
        identities = mapping[side]
        chains = {row["chain"] for row in identities}
        author_ids = {(row["number"], row["insertion_code"]) for row in identities}
        if (len(chains) != 1 or len(author_ids) != len(identities)
                or sorted(row["canonical_number"] for row in identities) != list(range(1, len(identities) + 1))):
            raise ValueError("The assessment's residue mapping is ambiguous.")
        chain_name = next(iter(chains))
        matching = [chain for chain in structure[0] if chain.name == chain_name]
        if len(matching) != 1 or chain_name in selected_chains:
            raise ValueError("The assessed target and binder must identify unique saved chains.")
        selected_chains.add(chain_name)
        chain = matching[0]
        copied = gemmi.Chain(side)
        by_author = {}
        for residue in chain:
            by_author.setdefault((residue.seqid.num, residue.seqid.icode.strip()), []).append(residue)
        identities = sorted(identities, key=lambda row: row["canonical_number"])
        residues[side], atoms[side] = [], {}
        for position, identity in enumerate(identities):
            matches = by_author.get((identity["number"], identity["insertion_code"])) or []
            if len(matches) != 1:
                raise ValueError("An assessed residue is absent or ambiguous in the saved coordinates.")
            residue = matches[0]
            number = identity["canonical_number"]
            info = gemmi.find_tabulated_residue(residue.name)
            ca = [atom for atom in residue if atom.name == "CA"]
            if not info.is_amino_acid() or len(ca) != 1:
                raise ValueError("An assessed protein residue has no unique Cα coordinate.")
            atom_map = {atom.name: [atom.pos.x, atom.pos.y, atom.pos.z] for atom in residue if not atom.element.is_hydrogen}
            if any(not math.isfinite(value) for point in atom_map.values() for value in point):
                raise ValueError("Saved coordinates contain non-finite atom positions.")
            atoms[side][number] = atom_map
            burial = analysis.get("dsasa_res") or []
            secondary = analysis.get("ss") or ""
            author = residue_label(result, side, number)
            aa = info.one_letter_code.upper()
            residues[side].append({"id": f"{side}:{number}", "side": side, "number": number, "aa": aa,
                "label": f"{aa}{identity['number']}{identity['insertion_code']}", "author": author,
                "chain": chain_name, "xyz": atom_map["CA"],
                "confidence": _finite(confidence[position]) if side == "B" and confidence_confirmed else None,
                "burial": _finite(burial[position]) if side == "B" and position < len(burial) else
                    _finite((analysis.get("t_dsasa") or {}).get(str(number), (analysis.get("t_dsasa") or {}).get(number))) if side == "A" else None,
                "ss": secondary[position] if side == "B" and position < len(secondary) else None,
                "interface": number in (analysis.get("tc" if side == "A" else "bc") or [])})
            saved = residue.clone()
            saved.seqid = gemmi.SeqId(number, " ")
            copied.add_residue(saved)
        if len(by_author) != len(identities):
            raise ValueError("Saved protein contents no longer match the assessed residue inventory.")
        display[0].add_chain(copied)
    extra_count = 0
    for chain in structure[0]:
        if chain.name in selected_chains:
            continue
        extra_count += 1
        if extra_count > 24:
            raise ValueError("Too many additional saved chains for this inspection view.")
        copied = chain.clone()
        copied.name = chr(ord("B") + extra_count)
        display[0].add_chain(copied)
    lookup = {row["id"]: row for rows in residues.values() for row in rows}
    contacts = []
    for row in contact_rows(result):
        target, binder = f"A:{row['canonical_target']}", f"B:{row['canonical_binder']}"
        if target not in lookup or binder not in lookup:
            raise ValueError("A contact references an unmapped residue.")
        contacts.append({"target": target, "binder": binder, "distance": row["Minimum distance (Å)"],
                         "pairs": row["Atom pairs"], "kind": row["Interaction"], "hotspot": row["Intended hotspot"]})
    interactions, missing = [], 0
    for row in analysis.get("hbonds") or []:
        target = atoms["A"].get(row["t"], {}).get(row.get("ta"))
        binder = atoms["B"].get(row["b"], {}).get(row.get("ba"))
        if target is None or binder is None:
            missing += 1
            continue
        interactions.append({"kind": "H-bond", "target": f"A:{row['t']}", "binder": f"B:{row['b']}",
            "start": target, "end": binder, "distance": row.get("d"), "basis": "Recorded donor/acceptor heavy atoms",
            "atoms": [row.get("ta"), row.get("ba")]})
    for row in analysis.get("salt") or []:
        endpoints = []
        for side, number, residue_name in (("A", row["t"], row.get("tn")), ("B", row["b"], row.get("bn"))):
            group = [point for atom, point in atoms[side].get(number, {}).items() if atom in CHARGED_ATOMS.get(residue_name, ())]
            endpoints.append([sum(point[axis] for point in group) / len(group) for axis in range(3)] if group else None)
        if any(point is None for point in endpoints):
            missing += 1
            continue
        interactions.append({"kind": "Salt bridge", "target": f"A:{row['t']}", "binder": f"B:{row['b']}",
            "start": endpoints[0], "end": endpoints[1], "distance": row.get("d"),
            "basis": "Opposite-charge side-chain group centroids, not minimum atom distance"})
    hotspots = [{**row, "id": f"A:{row['canonical_number']}" if row.get("canonical_number") is not None else None}
                for row in result.get("hotspot_compliance", {}).get("residues") or []]
    polar = [{**row, "id": f"{'A' if row['side'] == 't' else 'B'}:{row['resi']}"} for row in analysis.get("buns") or []]
    source = result.get("source") or {}
    return {"schema": "studio.binderscope-inspector.v1", "title": title or source.get("source", {}).get("record_name") or source.get("label") or "Saved design",
        "pdb": display.make_pdb_string(), "source_sha256": source.get("input_structure_digest"),
        "residues": residues, "contacts": contacts, "interactions": interactions, "unmapped_interactions": missing,
        "polar": polar, "hotspots": hotspots, "hotspot_status": result.get("hotspot_compliance", {}).get("status", "not_recorded"),
        "patches": analysis.get("patches") or [], "liabilities": analysis.get("liab") or [],
        "confidence_confirmed": confidence_confirmed,
        "target_label": result.get("_studio_target_label") or f"Saved target · {len(residues['A'])} residues",
        "reference": result.get("target_reference_comparison") or {"status": "unavailable", "reason": "No verified original target reference."},
        "hbond_basis": analysis.get("hbond_basis", "unavailable"), "extra_chains": extra_count,
        "target_caption": "Saved target numbering; biological offsets are not guessed from names.",
        "limitation": "Post-analysis of one saved conformation. Original admission remains unchanged. No specificity, full-protein context or PAE is inferred from missing inputs."}


@lru_cache(maxsize=4)
def _asset(name):
    return (ASSETS / name).read_text()


@lru_cache(maxsize=1)
def _viewer_library():
    return (Path(__file__).parent / "vendor/binderscope/binderscope/viewer/vendor/3Dmol-min.js").read_text()


def inspector_document(payload):
    data = json.dumps(payload, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    replacements = {"TITLE": html.escape(payload["title"]), "STYLE": _asset("inspector.css"),
        "VIEWER": _viewer_library().replace("</script", "<\\/script"), "DATA": data, "SCRIPT": _asset("inspector.js")}
    return re.sub(r"__STUDIO_(TITLE|STYLE|VIEWER|DATA|SCRIPT)__", lambda match: replacements[match[1]], _asset("index.html"))


def render_linked_inspector(result, database, user, *, key=""):
    if database is None or user is None:
        return False
    import streamlit as st
    import streamlit.components.v1 as components

    try:
        coordinates = result_coordinates(database, user, result)
        payload = build_view_payload(result, coordinates)
        document = inspector_document(payload)
    except (OSError, ValueError, PermissionError) as exc:
        st.caption("Linked inspector unavailable: " + str(exc))
        return False
    components.html(document, height=985, scrolling=True)
    suffix = "cif" if "_atom_site." in coordinates else "pdb"
    st.download_button("Download assessed complex", coordinates, f"binderscope-assessed-complex.{suffix}",
        "chemical/x-mmcif" if suffix == "cif" else "chemical/x-pdb", key=key + "-linked-complex-download")
    return True
