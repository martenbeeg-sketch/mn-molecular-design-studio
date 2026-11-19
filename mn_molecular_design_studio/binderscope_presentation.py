"""Native, representation-bound views of the integrated BinderScope measurements."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import math


INTERACTION_COLORS = {
    "Contact": "#829594", "H-bond": "#148e91", "Salt bridge": "#935bb5",
    "H-bond + salt bridge": "#be7826",
}


def residue_identity(result, side, number):
    rows = result.get("source", {}).get("residue_map", {}).get(side) or []
    matches = [row for row in rows if row["canonical_number"] == number]
    if len(matches) != 1:
        return None
    return matches[0]


def residue_label(result, side, number):
    row = residue_identity(result, side, number)
    if row is None:
        return f"Unmapped {side}:{number}"
    return f"{row['chain'] or '_'}:{row['number']}{row['insertion_code']}"


def residue_selection(result, side, number):
    row = residue_identity(result, side, number)
    if row is None:
        return None
    return f"{row['chain'] or '_'}:{row['number']}^{row['insertion_code'] or '_'}"


def contact_rows(result):
    analysis = result.get("analysis") or {}
    hbonds = {(row["t"], row["b"]) for row in analysis.get("hbonds") or []}
    salts = {(row["t"], row["b"]) for row in analysis.get("salt") or []}
    hotspots = {row["canonical_number"] for row in result.get("hotspot_compliance", {}).get("residues") or []
                if row.get("canonical_number") is not None}
    rows = []
    for contact in analysis.get("contacts") or []:
        pair = (contact["t"], contact["b"])
        interaction = " + ".join(kind for kind, inventory in (("H-bond", hbonds), ("salt bridge", salts)) if pair in inventory)
        interaction = interaction or "Contact"
        if interaction == "salt bridge":
            interaction = "Salt bridge"
        rows.append({"Target residue": residue_label(result, "A", pair[0]),
            "Binder residue": residue_label(result, "B", pair[1]), "Interaction": interaction,
            "Minimum distance (Å)": contact.get("d"), "Atom pairs": contact.get("n"),
            "Intended hotspot": pair[0] in hotspots,
            "canonical_target": pair[0], "canonical_binder": pair[1]})
    return rows


def selected_contact(event, rows):
    points = (event or {}).get("selection", {}).get("points") or []
    if not points:
        return None
    identity = points[-1].get("customdata") or []
    if len(identity) < 2:
        return None
    return next((row for row in rows if [row["canonical_target"], row["canonical_binder"]] == list(identity[:2])), None)


def highlight_selections(result, layer, contact=None):
    analysis = result.get("analysis") or {}
    residues = []
    if layer == "Interface":
        residues = [("A", number) for number in analysis.get("tc") or []]
        residues += [("B", number) for number in analysis.get("bc") or []]
    elif layer == "Intended hotspots":
        residues = [("A", row["canonical_number"]) for row in result.get("hotspot_compliance", {}).get("residues") or []
                    if row.get("canonical_number") is not None]
    elif layer == "Buried polar atoms":
        residues = [("A" if row["side"] == "t" else "B", row["resi"]) for row in analysis.get("buns") or []]
    elif layer == "Selected contact" and contact:
        residues = [("A", contact["canonical_target"]), ("B", contact["canonical_binder"])]
    return list(dict.fromkeys(selection for side, number in residues
                             if (selection := residue_selection(result, side, number)) is not None))


def result_coordinates(database, user, result):
    source = result.get("source") or {}
    binding = result.get("_studio_sql_binding")
    identifier = binding.get("structure_id") if binding is not None else source.get("source", {}).get("structure_id")
    expected = source.get("input_structure_digest")
    if not identifier or not expected:
        raise ValueError("This assessment has no verifiable SQL structure identity.")
    if binding and binding.get("structure_digest") and binding["structure_digest"] != expected:
        raise ValueError("SQL attachment and assessment disagree about the exact model checksum.")
    path = database.structure_file(user, identifier)
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != expected:
        raise ValueError("Saved coordinates do not match the exact model used for this assessment.")
    return content.decode("utf-8")


def _base_layout(figure, *, height=470):
    figure.update_layout(template="plotly_white", height=height, margin=dict(l=50, r=25, t=20, b=65),
                         font=dict(size=12), legend=dict(orientation="h", y=1.14, x=0), clickmode="event+select")
    return figure


def contact_map(result, rows=None):
    import plotly.graph_objects as graph

    rows = contact_rows(result) if rows is None else rows
    target_order = sorted({row["canonical_target"] for row in rows})
    binder_order = sorted({row["canonical_binder"] for row in rows})
    figure = graph.Figure()
    for interaction, color in INTERACTION_COLORS.items():
        subset = [row for row in rows if row["Interaction"] == interaction]
        if not subset:
            continue
        figure.add_trace(graph.Scatter(x=[row["canonical_target"] for row in subset],
            y=[row["canonical_binder"] for row in subset], mode="markers", name=interaction,
            marker=dict(color=color, size=[8 + min(16, 2 * math.sqrt(row["Atom pairs"] or 0)) for row in subset],
                        line=dict(width=[2 if row["Intended hotspot"] else 0 for row in subset], color="#e7a331")),
            customdata=[[row["canonical_target"], row["canonical_binder"], row["Target residue"], row["Binder residue"],
                         row["Minimum distance (Å)"], row["Atom pairs"], row["Intended hotspot"]] for row in subset],
            hovertemplate="Target %{customdata[2]} · binder %{customdata[3]}<br>%{customdata[4]:.2f} Å · %{customdata[5]} atom pairs"
                          "<br>Intended hotspot: %{customdata[6]}<extra>%{fullData.name}</extra>"))
    figure.update_xaxes(title="Target · saved residue numbering", tickvals=target_order,
                        ticktext=[residue_label(result, "A", number) for number in target_order], tickangle=-55)
    figure.update_yaxes(title="Binder · saved residue numbering", tickvals=binder_order,
                        ticktext=[residue_label(result, "B", number) for number in binder_order])
    return _base_layout(figure)


def contact_network(result, rows=None):
    import plotly.graph_objects as graph

    rows = contact_rows(result) if rows is None else rows
    target = sorted({row["canonical_target"] for row in rows})
    binder = sorted({row["canonical_binder"] for row in rows})
    neighbors = {"A": defaultdict(list), "B": defaultdict(list)}
    for row in rows:
        neighbors["A"][row["canonical_target"]].append(row["canonical_binder"])
        neighbors["B"][row["canonical_binder"]].append(row["canonical_target"])
    for sweep in range(6):
        positions = {number: position for position, number in enumerate(binder)}
        target.sort(key=lambda number: (sum(positions[value] for value in neighbors["A"][number]) / len(neighbors["A"][number]), number))
        positions = {number: position for position, number in enumerate(target)}
        binder.sort(key=lambda number: (sum(positions[value] for value in neighbors["B"][number]) / len(neighbors["B"][number]), number))
    positions = {side: {number: position / max(1, len(numbers) - 1) for position, number in enumerate(numbers)}
                 for side, numbers in (("A", target), ("B", binder))}
    figure = graph.Figure()
    for interaction, color in INTERACTION_COLORS.items():
        subset = [row for row in rows if row["Interaction"] == interaction]
        if not subset:
            continue
        horizontal, vertical = [], []
        for row in subset:
            horizontal.extend([0, 1, None])
            vertical.extend([positions["A"][row["canonical_target"]], positions["B"][row["canonical_binder"]], None])
        figure.add_trace(graph.Scatter(x=horizontal, y=vertical, mode="lines", name=interaction,
                                     line=dict(color=color, width=1.5), hoverinfo="skip"))
        figure.add_trace(graph.Scatter(x=[0.5] * len(subset),
            y=[(positions["A"][row["canonical_target"]] + positions["B"][row["canonical_binder"]]) / 2 for row in subset],
            mode="markers", showlegend=False, marker=dict(color=color, size=9),
            customdata=[[row["canonical_target"], row["canonical_binder"]] for row in subset],
            text=[f"{row['Target residue']} ↔ {row['Binder residue']} · {row['Minimum distance (Å)']} Å · {row['Atom pairs']} atom pairs" for row in subset],
            hovertemplate="%{text}<extra>%{fullData.name}</extra>"))
    hotspots = {row["canonical_target"] for row in rows if row["Intended hotspot"]}
    for side, numbers, horizontal, color, text_position in (("A", target, 0, "#477aa9", "middle left"), ("B", binder, 1, "#318c70", "middle right")):
        figure.add_trace(graph.Scatter(x=[horizontal] * len(numbers), y=[positions[side][number] for number in numbers],
            mode="markers+text", text=[residue_label(result, side, number) for number in numbers], textposition=text_position,
            marker=dict(size=9, color=["#e7a331" if side == "A" and number in hotspots else color for number in numbers]),
            showlegend=False, hoverinfo="text"))
    figure.update_xaxes(visible=False, range=[-0.35, 1.35])
    figure.update_yaxes(visible=False, range=[1.04, -0.04])
    return _base_layout(figure, height=max(430, 22 * max(len(target), len(binder))))


def residue_rows(result):
    analysis = result.get("analysis") or {}
    residues = result.get("source", {}).get("residue_map", {}).get("B") or []
    tracks = analysis.get("tracks") or {}
    interface = set(analysis.get("bc") or [])
    rows = []
    for position, residue in enumerate(residues):
        number = residue["canonical_number"]
        def value(track):
            values = analysis.get(track) or []
            return values[position] if position < len(values) else None
        contacts = tracks.get("n_target_res") or []
        rows.append({"Residue": residue_label(result, "B", number), "canonical_number": number,
            "pLDDT": value("plddt_res"), "Buried area (Å²)": value("dsasa_res"), "Relative SASA": value("rsa_res"),
            "Secondary structure": value("ss"), "Interface": number in interface,
            "Contacted target residues": contacts[position] if position < len(contacts) else None,
            "Liabilities": "; ".join(row["kind"] for row in analysis.get("liab") or [] if row.get("pos") == number)})
    return rows


def residue_profile(result):
    import plotly.graph_objects as graph
    from plotly.subplots import make_subplots

    rows = residue_rows(result)
    figure = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.1,
                           subplot_titles=("Binding burial · binder side", "Per-residue confidence · binder only"))
    labels = [row["Residue"] for row in rows]
    figure.add_trace(graph.Bar(x=labels, y=[row["Buried area (Å²)"] for row in rows], name="Buried area (Å²)",
        marker_color=["#318c70" if row["Interface"] else "#b9c6c4" for row in rows]), row=1, col=1)
    if any(row["pLDDT"] is not None for row in rows):
        figure.add_trace(graph.Scatter(x=labels, y=[row["pLDDT"] for row in rows], name="pLDDT", mode="lines", line_color="#477aa9"), row=2, col=1)
    else:
        figure.add_annotation(text="No confirmed per-residue pLDDT; B-factors are not relabelled.", x=0.5, y=0.2, xref="paper", yref="paper", showarrow=False)
    figure.update_yaxes(title_text="Å²", row=1, col=1)
    figure.update_yaxes(title_text="pLDDT (0–100)", range=[0, 100], row=2, col=1)
    return _base_layout(figure, height=470)
