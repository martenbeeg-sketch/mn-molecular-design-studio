"""MD interaction diagrams with source-named ligand-atom endpoints."""
from __future__ import annotations

from dataclasses import dataclass
from html import escape
import math

from mn_molecular_design_studio.result_details_ui import render_result_details


INTERACTIONS = {
    "contact": ("Geometric contact", "contact_occupancy", "#64748b"),
    "hydrogen_bond": ("Hydrogen bond", "hydrogen_bond_occupancy", "#2563eb"),
    "hydrophobic": ("Hydrophobic", "hydrophobic_occupancy", "#b8791c"),
    "water_bridge": ("Water bridge", "water_bridge_occupancy", "#0891b2"),
    "salt_bridge": ("Salt bridge", "salt_bridge_occupancy", "#a855b5"),
}


def _occupancy(value):
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) and 0 <= number <= 1 else None


def interaction_residues(structural):
    contacts = structural.get("contacts")
    contacts = contacts if isinstance(contacts, dict) else {}
    network = structural.get("interaction_network")
    network = network if isinstance(network, dict) else {}
    source = contacts.get("residues")
    if not isinstance(source, list) or not source:
        source = network.get("nodes")
    if not isinstance(source, list):
        return []
    grouped = {}
    for row in source:
        if not isinstance(row, dict):
            continue
        residue = str(row.get("residue") or row.get("id") or "").strip()
        if not residue:
            continue
        values = {kind: value for kind, (_, field, _) in INTERACTIONS.items()
            if (value := _occupancy(row.get(field))) is not None}
        grouped.setdefault(residue, []).append({"residue": residue, "occupancies": values,
            "ligand_atom": str(row.get("top_ligand_atom") or "").strip(),
            "atom_endpoints": row.get("atom_endpoints") if isinstance(row.get("atom_endpoints"), dict) else {},
            "protein_atom": str(row.get("top_protein_atom") or "").strip()})
    return [rows[0] for _, rows in sorted(grouped.items()) if all(row == rows[0] for row in rows)]


@dataclass(frozen=True)
class LigandDrawing:
    svg: str
    atoms: dict[str, tuple[float, float]]
    atom_mapping_note: str


def _named_ligand(smiles, pdb_text):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    records = [line for line in pdb_text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    residues = {(line[17:20], line[21:22], line[22:27]) for line in records}
    if not records or len(records) > 1200 or len(residues) != 1:
        raise ValueError("The saved ligand block must identify one residue.")
    template = Chem.MolFromSmiles(smiles)
    source = Chem.MolFromPDBBlock(pdb_text, sanitize=False, removeHs=False)
    if template is None or source is None:
        raise ValueError("The saved ligand chemistry could not be read.")
    source = Chem.RemoveHs(source, sanitize=False)
    if not 1 <= source.GetNumAtoms() <= 300 or source.GetNumAtoms() != template.GetNumHeavyAtoms():
        raise ValueError("The saved ligand atom count differs from its linked chemical identity.")
    names = [atom.GetPDBResidueInfo().GetName().strip() if atom.GetPDBResidueInfo() else "" for atom in source.GetAtoms()]
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("Saved ligand atom names are missing or ambiguous.")
    molecule = AllChem.AssignBondOrdersFromTemplate(template, source)
    Chem.AssignStereochemistryFrom3D(molecule, replaceExistingTags=True)
    if Chem.MolToSmiles(molecule, isomericSmiles=False) != Chem.MolToSmiles(template, isomericSmiles=False) or not molecule.HasSubstructMatch(template, useChirality=True):
        raise ValueError("The saved ligand block disagrees with the linked chemical identity.")
    if [atom.GetPDBResidueInfo().GetName().strip() for atom in molecule.GetAtoms()] != names:
        raise ValueError("Ligand atom names changed while assigning chemical bonds.")
    return molecule, names


def ligand_drawing(smiles, *, ligand_pdb="", label_atoms=()):
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > 4096:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D

        names = []
        note = "A saved, matching ligand atom-name block is unavailable; atom endpoints are not guessed."
        molecule = None
        if isinstance(ligand_pdb, str) and ligand_pdb and len(ligand_pdb) <= 256 * 1024:
            try:
                molecule, names = _named_ligand(smiles, ligand_pdb)
                note = "Original atom names are retained from the source ligand block, matched to its linked chemical identity."
            except (RuntimeError, ValueError) as exc:
                note = str(exc) + " Atom endpoints are not guessed."
        if molecule is None:
            molecule = Chem.MolFromSmiles(smiles)
        if molecule is None or not 1 <= molecule.GetNumAtoms() <= 300:
            return None
        rdDepictor.Compute2DCoords(molecule)
        drawer = rdMolDraw2D.MolDraw2DSVG(440, 330)
        drawer.drawOptions().clearBackground = False
        drawer.drawOptions().padding = 0.08
        for atom_index, name in enumerate(names):
            if name in label_atoms:
                drawer.drawOptions().atomLabels[atom_index] = name
        drawer.DrawMolecule(molecule)
        drawer.FinishDrawing()
        positions = {name: (230 + drawer.GetDrawCoords(index).x, 200 + drawer.GetDrawCoords(index).y)
            for index, name in enumerate(names)}
        drawing = drawer.GetDrawingText()
        svg = drawing[drawing.index("<svg"):].replace("<svg", '<svg x="230" y="200"', 1)
        return LigandDrawing(svg, positions, note)
    except (ImportError, RuntimeError, ValueError):
        return None


@dataclass(frozen=True)
class InteractionDiagram:
    svg: str
    shown: int
    eligible: int
    chemical_depiction: bool
    mapped_residues: int
    atom_mapping_note: str
    atom_rows: tuple[dict, ...]


def _node_positions(rows, drawing):
    from mn_molecular_design_studio.md_interaction_layout import distribute_labels, reduce_crossings

    atoms = drawing.atoms if drawing else {}
    anchors = [(atoms[row["ligand_atom"]][0] - 450, atoms[row["ligand_atom"]][1] - 365)
        if row["ligand_atom"] in atoms else (0, 0) for row in rows]
    positions = distribute_labels(anchors, x_radius=350, y_radius=255)
    return reduce_crossings(anchors, positions) + (450, 365)


def build_interaction_diagram(structural, *, ligand_smiles="", ligand_pdb="", ligand_label="Ligand",
                              interaction="all", minimum_occupancy=0.1, residue_limit=12, scope="Selected production segment"):
    if interaction != "all" and interaction not in INTERACTIONS:
        raise ValueError("Unknown recorded interaction type.")
    threshold = _occupancy(minimum_occupancy)
    if threshold is None or isinstance(residue_limit, bool) or not isinstance(residue_limit, int) or not 1 <= residue_limit <= 16:
        raise ValueError("Invalid diagram occupancy threshold or residue limit.")
    selected = []
    for row in interaction_residues(structural):
        values = {kind: value for kind, value in row["occupancies"].items()
            if value > 0 and value >= threshold and (interaction == "all" or kind == interaction)}
        if interaction == "all" and any(kind != "contact" for kind in values):
            values.pop("contact", None)
        if values:
            selected.append({**row, "occupancies": values})
    selected.sort(key=lambda row: (-max(row["occupancies"].values()), row["residue"]))
    eligible = len(selected)
    selected = selected[:residue_limit]
    def ligand_atom(row, kind):
        endpoint = row["atom_endpoints"].get(kind) or {}
        return endpoint.get("ligand_atom") or row["ligand_atom"]

    atom_names = {ligand_atom(row, kind) for row in selected for kind in row["occupancies"]}
    drawing = ligand_drawing(ligand_smiles, ligand_pdb=ligand_pdb, label_atoms=atom_names) if selected else None
    pose = structural.get("input_type") == "pose"
    if pose:
        description = "Interactions in one static pose."
    elif structural.get("method_version"):
        description = "Percentages report residue-level trajectory occupancy. Lines use one representative atom pair per interaction type."
    else:
        description = ("Endpoints use recorded representative ligand atoms (top_ligand_atom). "
                       "Percentages are residue-level occupancies, not atom-pair occupancies or atom-specific bonds.")
    note = drawing.atom_mapping_note if drawing else "Chemical depiction and atom-name mapping are unavailable."
    atom_positions = drawing.atoms if drawing else {}
    mapped = sum(any(ligand_atom(row, kind) in atom_positions for kind in row["occupancies"]) for row in selected)
    positions = _node_positions(selected, drawing)
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="900" height="790" viewBox="0 0 900 790" role="img">',
        '<title>Recorded ligand atom–residue associations</title>',
        '<desc>' + escape(str(scope) + ". " + description + " Recorded analysis window: " +
            str(structural.get("analysis_window") or "not recorded")) + '</desc>',
        '<rect width="900" height="790" rx="24" fill="#f8fafc"/>',
        '<g font-family="Arial, sans-serif" fill="#1e293b">',
        '<text x="450" y="35" text-anchor="middle" font-size="23" font-weight="bold">Ligand atom–residue associations</text>',
        '<text x="450" y="61" text-anchor="middle" font-size="14" fill="#475569">' +
        ("One saved pose · named atoms from the interaction pairs" if pose else "Named atom = representative · % = residue occupancy") + '</text>']
    radius = 42 if len(selected) > 12 else 49
    for row, (node_x, node_y) in zip(selected, positions):
        for offset, (kind, occupancy) in enumerate(row["occupancies"].items()):
            name = ligand_atom(row, kind)
            atom = atom_positions.get(name)
            if atom is None:
                continue
            direction_x, direction_y = node_x - atom[0], node_y - atom[1]
            distance = math.hypot(direction_x, direction_y)
            if distance <= radius:
                continue
            end_x, end_y = node_x - radius * direction_x / distance, node_y - radius * direction_y / distance
            bend = (offset - (len(row["occupancies"]) - 1) / 2) * 24
            control_x = (atom[0] + end_x) / 2 - direction_y / distance * bend
            control_y = (atom[1] + end_y) / 2 + direction_x / distance * bend
            color = INTERACTIONS[kind][2]
            measure = "" if pose else f" {occupancy:.1%} residue occupancy"
            title = escape(f'{row["residue"]} · {INTERACTIONS[kind][0]}{measure} · ligand atom {name}')
            parts.append(f'<path data-ligand-atom="{escape(name)}" data-residue="{escape(row["residue"])}" '
                f'd="M {atom[0]:.2f},{atom[1]:.2f} Q {control_x:.2f},{control_y:.2f} {end_x:.2f},{end_y:.2f}" '
                f'stroke="{color}" stroke-width="{1.2 + 2.5 * occupancy:.2f}" stroke-dasharray="5 4" opacity="0.72" fill="none"><title>{title}</title></path>')
    if drawing:
        parts.append(drawing.svg)
    else:
        parts.extend(['<circle cx="450" cy="365" r="65" fill="#e2ece8"/>',
            '<text x="450" y="365" text-anchor="middle" font-size="20">Ligand</text>',
            '<text x="450" y="463" text-anchor="middle" font-size="14">Chemical depiction unavailable</text>'])
    label = str(ligand_label)
    parts.append(f'<text x="450" y="557" text-anchor="middle" font-size="14" font-weight="bold">{escape(label[:42] + ("…" if len(label) > 42 else ""))}</text>')
    for row, (node_x, node_y) in zip(selected, positions):
        values = row["occupancies"]
        dominant = max(values, key=values.get)
        residue, separator, chain = row["residue"].partition(" · ")
        mapped_names = sorted({ligand_atom(row, kind) for kind in values} & set(atom_positions))
        atom_name = ", ".join(mapped_names) or "unmapped"
        title = escape(row["residue"] + ": " + "; ".join(
            INTERACTIONS[kind][0] if pose else f"{INTERACTIONS[kind][0]} {value:.1%}" for kind, value in values.items()
        ))
        measure = "" if pose else f"{max(values.values()):.1%}"
        parts.extend([f'<g><title>{title}</title><circle cx="{node_x:.2f}" cy="{node_y:.2f}" r="{radius}" fill="white" stroke="{INTERACTIONS[dominant][2]}" stroke-width="2"/>',
            f'<text x="{node_x:.2f}" y="{node_y - 18:.2f}" text-anchor="middle" font-size="14" font-weight="bold">{escape(residue[:17])}</text>',
            f'<text x="{node_x:.2f}" y="{node_y - 1:.2f}" text-anchor="middle" font-size="11" fill="#64748b">{escape(chain[:18] if separator else "Residue")}</text>',
            f'<text x="{node_x:.2f}" y="{node_y + 16:.2f}" text-anchor="middle" font-size="14">{measure}</text>',
            f'<text x="{node_x:.2f}" y="{node_y + 33:.2f}" text-anchor="middle" font-size="11">atom: {escape(atom_name[:16])}</text></g>'])
    for atom_name in sorted(atom_names & set(atom_positions)):
        atom_x, atom_y = atom_positions[atom_name]
        parts.append(f'<circle data-atom-name="{escape(atom_name)}" cx="{atom_x:.2f}" cy="{atom_y:.2f}" r="4" fill="none" stroke="#0f766e" stroke-width="1.5"/>')
    used_kinds = {kind for row in selected for kind in row["occupancies"]}
    legend = [value for kind, value in INTERACTIONS.items() if kind in used_kinds]
    for position, (title, _, color) in enumerate(legend):
        legend_x, legend_y = 80 + (position % 3) * 275, 690 + (position // 3) * 28
        parts.append(f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 26}" y2="{legend_y}" stroke="{color}" stroke-width="3" stroke-dasharray="5 4"/>'
            f'<text x="{legend_x + 35}" y="{legend_y + 5}" font-size="14">{escape(title)}</text>')
    parts.extend([f'<text x="450" y="759" text-anchor="middle" font-size="13" fill="#64748b">{len(selected)} of {eligible} residues · {mapped} named-atom endpoints · no guessed atom links</text>', '</g></svg>'])
    atom_rows = tuple({"Residue": row["residue"], "Recorded representative ligand atom": row["ligand_atom"] or "not recorded",
        "Reported protein atom (separate statistic)": row["protein_atom"] or "not recorded",
        "Atom endpoint": "mapped" if row["ligand_atom"] in atom_positions else "unavailable / ambiguous",
        **({INTERACTIONS[kind][0]: ligand_atom(row, kind) for kind in row["occupancies"]} if pose else
           {INTERACTIONS[kind][0] + " residue occupancy (%)": 100 * value for kind, value in row["occupancies"].items()})}
        for row in selected)
    return InteractionDiagram("".join(parts), len(selected), eligible, drawing is not None, mapped, note, atom_rows)


def render_interaction_diagram(st, structural, *, ligand_smiles="", ligand_pdb="", atom_source_note="", ligand_label="Ligand", key_prefix="md", scope="Selected production segment"):
    rows = interaction_residues(structural)
    st.markdown("**Ligand atom–residue interactions**")
    if not rows or not any(value > 0 for row in rows for value in row["occupancies"].values()):
        st.info("No valid residue interaction occupancies are recorded for this segment. A contact map is not inferred from another structure or replica.")
        return
    with st.expander("Interaction diagram options"):
        interaction = st.selectbox("Interaction type", ["all", *INTERACTIONS],
            format_func=lambda value: "All recorded types" if value == "all" else INTERACTIONS[value][0], key=key_prefix + "-interaction-kind")
        threshold = 0 if structural.get("input_type") == "pose" else st.slider(
            "Minimum occupancy (%)", 0, 100, 10, key=key_prefix + "-interaction-threshold")
        limit = st.slider("Maximum residues shown", 4, 16, 12, key=key_prefix + "-interaction-limit")
    diagram = build_interaction_diagram(structural, ligand_smiles=ligand_smiles, ligand_pdb=ligand_pdb,
        ligand_label=ligand_label, interaction=interaction, minimum_occupancy=threshold / 100, residue_limit=limit, scope=scope)
    if not diagram.shown:
        st.info("No recorded interactions match these display settings. Lower the occupancy threshold or select another interaction type.")
        return
    if structural.get("input_type") == "pose":
        caption = scope + ". Dashed lines link the residue to a measured ligand atom. The pair table identifies distance-only evidence."
    elif structural.get("method_version"):
        caption = scope + ". Each line ends at the most frequent ligand atom for that interaction type. The pair table lists occupancies for individual atom pairs."
    else:
        caption = (scope + ". Dashed lines end at each residue's recorded representative ligand atom (top_ligand_atom). "
        "Percentages are residue-level trajectory occupancies, not atom-pair occupancies. "
        "The archive does not identify the exact ligand atom for each hydrogen bond; these are representative atom associations, not proven atom-specific bonds.")
    st.caption(caption)
    st.caption(" ".join(value for value in (diagram.atom_mapping_note, atom_source_note) if value))
    if diagram.mapped_residues < diagram.shown:
        st.warning(f"{diagram.shown - diagram.mapped_residues} displayed residue(s) have no verified atom endpoint. Their atom links are omitted, not placed on an arbitrary ligand atom.")
    st.image(diagram.svg, width="stretch")
    st.caption(f"Showing {diagram.shown} of {diagram.eligible} matching residues. The layout is schematic, not a spatial binding-site map or binding-energy estimate.")
    if st.button("Enlarge interaction diagram", key=key_prefix + "-interaction-enlarge"):
        @st.dialog("Ligand atom–residue interactions", width="large")
        def enlarged():
            with st.container(width=650):
                st.image(diagram.svg, width="stretch")
            st.caption(caption)
        enlarged()
    st.download_button("Download interaction diagram (SVG)", diagram.svg, file_name="md_ligand_atom_interactions.svg",
        mime="image/svg+xml", key=key_prefix + "-interaction-download")
    with st.expander("Atom assignments and recorded residue occupancies"):
        st.dataframe(list(diagram.atom_rows), hide_index=True, width="stretch")
        window = structural.get("analysis_window")
        if isinstance(window, dict) and window:
            render_result_details(st, window)
