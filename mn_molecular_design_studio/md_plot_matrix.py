"""Compact, replica-aware plots reconstructed from saved MN-Ligand analyses."""
from __future__ import annotations

import io
import math

from mn_molecular_design_studio.md_atom_interactions import INTERACTIONS, _occupancy
from mn_molecular_design_studio.md_interaction_layout import distribute_labels, place_edge_percentages, reduce_crossings


COLORS = {"contact": "#777777", "hydrogen_bond": "#800000", "hydrophobic": "#62a854",
    "water_bridge": "#1e90ff", "salt_bridge": "#cc33aa"}
TIME_PLOTS = {
    "backbone_rmsd": ("Protein backbone RMSD", "RMSD (Å)", ("backbone_rmsd_angstrom",)),
    "ligand_rmsd": ("Ligand RMSD", "RMSD (Å)", ("ligand_rmsd_angstrom",)),
    "ligand_displacement": ("Ligand displacement", "Displacement (Å)", ("ligand_centroid_displacement_angstrom",)),
    "minimum_distance": ("Minimum protein–ligand distance", "Distance (Å)", ("minimum_protein_distance_angstrom",)),
    "radius_of_gyration": ("Radius of gyration", "Rg (Å)", ("protein_rg_angstrom", "ligand_rg_angstrom", "complex_rg_angstrom")),
    "secondary_structure": ("Secondary structure", "Fraction", ("helix_fraction", "sheet_fraction", "coil_fraction")),
}
SUMMARY_PLOTS = {
    "site_retention": ("Reference-site retention", "Fraction", "reference_site_retained_fraction"),
    "endpoint_energy": ("MM/GBSA binding energy", "kcal/mol", "delta_g_bind_kcal_mol"),
    "endpoint_energy_pb": ("MM/PBSA binding energy", "kcal/mol", "pb_delta_g_bind_kcal_mol"),
    "throughput": ("Simulation throughput", "ns/day", "performance_ns_per_day"),
}
PLOT_LABELS = {"interaction_network": "Ligand–residue interaction map",
    **{key: value[0] for key, value in TIME_PLOTS.items()},
    **{kind + "_occupancy": title + " occupancy" for kind, (title, _, _) in INTERACTIONS.items()},
    "protein_rmsf": "Protein Cα RMSF", "ligand_rmsf": "Ligand atom RMSF",
    "contact_matrix": "Contact occupancy by replica",
    **{key: value[0] for key, value in SUMMARY_PLOTS.items()}}


def _number(value):
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def _rows(value):
    return [row for row in value if isinstance(row, dict)] if isinstance(value, list) else []


def replica_ids(report):
    rows = _rows(report.get("replicas"))
    identifiers = [row.get("replica") for row in rows]
    if any(isinstance(value, bool) or not isinstance(value, int) or value <= 0 for value in identifiers) or len(set(identifiers)) != len(identifiers):
        raise ValueError("The stored analysis has missing or ambiguous replica identities.")
    return sorted(row["replica"] for row in rows if row.get("status") in {"completed", "success"} and row.get("success") is not False)


def validate_report_sources(report, productions):
    if not report.get("replicas"):
        raise ValueError("This report does not identify its source production replicas.")
    replica_ids(report)
    for row in _rows(report.get("replicas")):
        source = str(row.get("run_id") or "")
        if not source or source not in productions or str(productions[source]) != str(row["replica"]):
            raise ValueError("An analysis replica is not linked to an accessible production in this study.")


def _replica_values(report, field, replica):
    matrix = report.get("contact_matrix")
    if not isinstance(matrix, dict):
        return {}
    axis = matrix.get("replicas")
    residues = matrix.get("residues")
    values = matrix.get(field)
    if not isinstance(axis, list) or not isinstance(residues, list) or not isinstance(values, list):
        return {}
    if any(isinstance(value, bool) or not isinstance(value, int) for value in axis) or len(set(axis)) != len(axis) or axis.count(replica) != 1:
        return {}
    if any(not isinstance(value, str) for value in residues) or len(set(residues)) != len(residues):
        return {}
    column = axis.index(replica)
    return {residue: row[column] for residue, row in zip(residues, values)
        if isinstance(row, list) and column < len(row)}


def contact_rows(report, replica=None, *, limit=12):
    rows = _rows(report.get("contact_consensus"))
    counts = {}
    for row in rows:
        label = row.get("residue")
        if isinstance(label, str) and label:
            counts[label] = counts.get(label, 0) + 1
    rows = [row for row in rows if counts.get(row.get("residue"), 0) == 1]
    ranked = any(_number(row.get("mean_binding_importance_score")) is not None for row in rows)
    score = "mean_binding_importance_score" if ranked else "mean_contact_occupancy"
    rows.sort(key=lambda row: (-(_number(row.get(score)) or 0), row["residue"]))
    fields = [field for _, field, _ in INTERACTIONS.values()]
    mapped = {field: _replica_values(report, field, replica) for field in fields} if replica is not None else {}
    return [{"residue": row["residue"], "ligand_atom": str(row.get("top_ligand_atom") or ""),
        "occupancies": {kind: value for kind, (_, field, _) in INTERACTIONS.items()
            if (value := _occupancy(row.get("mean_" + field) if replica is None else mapped[field].get(row["residue"]))) is not None}}
        for row in rows[:limit]]


def _source_depiction(report):
    depiction = report.get("ligand_depiction")
    if not isinstance(depiction, dict):
        return {}, []
    atoms, bonds = _rows(depiction.get("atoms")), _rows(depiction.get("bonds"))
    if not 1 <= len(atoms) <= 300 or len(bonds) > 900:
        return {}, []
    names = [atom.get("name") for atom in atoms]
    if any(not isinstance(name, str) or not name for name in names) or len(set(names)) != len(names):
        return {}, []
    if any(_number(atom.get(field)) is None for atom in atoms for field in ("x", "y")):
        return {}, []
    horizontal = [float(atom["x"]) for atom in atoms]
    vertical = [float(atom["y"]) for atom in atoms]
    extent_x, extent_y = max(horizontal) - min(horizontal), max(vertical) - min(vertical)
    scale = min(2.35 / max(extent_x, 0.1), 1.65 / max(extent_y, 0.1))
    if not math.isfinite(scale):
        return {}, []
    center_x, center_y = min(horizontal) + extent_x / 2, min(vertical) + extent_y / 2
    positioned = {name: {**atom, "xy": ((float(atom["x"]) - center_x) * scale,
        (float(atom["y"]) - center_y) * scale + .16)} for name, atom in zip(names, atoms)}
    prepared = []
    for bond in bonds:
        begin, end = bond.get("begin"), bond.get("end")
        if any(isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(atoms) for index in (begin, end)):
            return {}, []
        order = _number(bond.get("order"))
        if begin == end or order not in (1, 1.5, 2, 3):
            return {}, []
        prepared.append((names[begin], names[end], order))
    return positioned, prepared


def interaction_figure(report, replica=None, *, residue_limit=12):
    from matplotlib.figure import Figure
    from matplotlib.patches import FancyArrowPatch

    figure = Figure(figsize=(4.2, 3.0))
    axis = figure.subplots()
    atoms, bonds = _source_depiction(report)
    rows = contact_rows(report, replica, limit=residue_limit)
    if not atoms:
        axis.text(.5, .5, "Named ligand depiction not recorded\nNo atom endpoints are guessed", transform=axis.transAxes,
            ha="center", va="center", fontsize=8)
        axis.axis("off")
        return figure
    for begin, end, order in bonds:
        start, finish = atoms[begin]["xy"], atoms[end]["xy"]
        axis.plot([start[0], finish[0]], [start[1], finish[1]], color="#333333", linewidth=.8 + .25 * (order - 1), zorder=3)
    atom_colors = {"O": "#d73027", "N": "#2166c2", "S": "#b18a00", "F": "#37a055", "Cl": "#37a055"}
    for atom in atoms.values():
        horizontal, vertical = atom["xy"]
        symbol = str(atom.get("symbol") or "C")
        if symbol != "C":
            axis.text(horizontal, vertical, symbol, color=atom_colors.get(symbol, "#444444"), fontsize=6,
                ha="center", va="center", bbox={"facecolor": "white", "edgecolor": "none", "pad": .1}, zorder=5)
    visible = []
    used = set()
    center = (0, .16)
    anchors = [(atoms[row["ligand_atom"]]["xy"][0], atoms[row["ligand_atom"]]["xy"][1] - center[1])
        if row["ligand_atom"] in atoms else (0, 0) for row in rows]
    positions = distribute_labels(anchors, x_radius=2.18, y_radius=1.40)
    positions = reduce_crossings(anchors, positions)
    positions += center
    edge_labels = []
    for row, end in zip(rows, positions):
        source = atoms.get(row["ligand_atom"])
        values = {kind: value for kind, value in row["occupancies"].items() if value > 0}
        if any(kind != "contact" for kind in values):
            values.pop("contact", None)
        if not values:
            continue
        if source:
            start = source["xy"]
            for style_index, (kind, occupancy) in enumerate(values.items()):
                curve = .09 * (style_index - (len(values) - 1) / 2)
                patch = FancyArrowPatch(start, end, connectionstyle=f"arc3,rad={curve}", arrowstyle="-",
                    linestyle=":" if kind == "contact" else "--", linewidth=.55 + 3 * occupancy,
                    color=COLORS[kind], alpha=.45 + .5 * occupancy, shrinkA=2, shrinkB=10, zorder=1)
                patch.set_gid(f"{row['residue']}|{row['ligand_atom']}|{kind}|{occupancy}")
                axis.add_patch(patch)
                text = axis.text(0, 0, f"{occupancy:.1%}", fontsize=5.4, color=COLORS[kind],
                    ha="center", va="center", bbox={"facecolor": "white", "edgecolor": "none", "alpha": .82, "pad": .3}, zorder=6)
                text.set_gid("occupancy|" + patch.get_gid())
                edge_labels.append((patch, text))
                used.add(kind)
        residue, _, chain = row["residue"].partition(" · ")
        text = residue + ("\n" + chain if chain else "") + ("\n(atom unavailable)" if not source else "")
        axis.scatter(*end, s=360, facecolor="#e9f1e6", edgecolor="#4c738c", linewidth=.7, zorder=2)
        axis.text(*end, text, fontsize=5.1, ha="center", va="center", zorder=4)
        visible.append(row)
    if not visible:
        axis.text(0, -1.3, "No valid positive interaction occupancies recorded", ha="center", fontsize=7)
    for name in sorted({row["ligand_atom"] for row in visible} & set(atoms)):
        axis.annotate(name, atoms[name]["xy"], xytext=(3, 4), textcoords="offset points",
            fontsize=4.2, color="#475569", zorder=5)
    for kind in INTERACTIONS:
        if kind in used:
            axis.plot([], [], color=COLORS[kind], linestyle="--", label=INTERACTIONS[kind][0], linewidth=1.5)
    if used:
        axis.legend(loc="lower center", bbox_to_anchor=(.5, -.1), ncol=3, fontsize=4.8, frameon=False)
    axis.set(xlim=(-2.58, 2.58), ylim=(-1.78, 1.9), aspect="equal")
    axis.axis("off")
    figure.subplots_adjust(left=.02, right=.98, top=.93, bottom=.14)
    molecule_positions = [atom["xy"] for atom in atoms.values()]
    molecule_bounds = [(min(point[0] for point in molecule_positions), min(point[1] for point in molecule_positions)),
        (max(point[0] for point in molecule_positions), max(point[1] for point in molecule_positions))]
    place_edge_percentages(figure, axis, edge_labels, residue_positions=positions, molecule_bounds=molecule_bounds)
    return figure


def _series(report, replica):
    rows = _rows(report.get("replica_series"))
    matches = [row for row in rows if row.get("replica") == replica]
    return matches[0] if len(matches) == 1 else {}


def available_plots(report):
    series = _rows(report.get("replica_series"))
    available = []
    if report.get("contact_consensus"):
        available.append("interaction_network")
    available += [key for key, (_, _, fields) in TIME_PLOTS.items()
        if any(isinstance(row.get(field), list) and any(_number(value) is not None for value in row[field]) for row in series for field in fields)]
    available += [kind + "_occupancy" for kind in INTERACTIONS
        if (kind != "salt_bridge" or report.get("salt_bridge_applicable") is not False)
        if any(kind in row["occupancies"] for row in contact_rows(report, limit=10000))]
    for key, field in (("protein_rmsf", "protein_rmsf_angstrom"), ("ligand_rmsf", "ligand_rmsf_angstrom")):
        if any(isinstance(row.get(field), list) and row[field] for row in series):
            available.append(key)
    if isinstance(report.get("contact_matrix"), dict) and report["contact_matrix"].get("contact_occupancy"):
        available.append("contact_matrix")
    available += [key for key, (_, _, field) in SUMMARY_PLOTS.items()
        if any(_number(row.get(field)) is not None for row in _rows(report.get("replicas")))]
    return available


def plot_figure(report, plot, replica=None, *, residue_limit=12, selected_replicas=None):
    from matplotlib.figure import Figure

    if plot == "interaction_network":
        return interaction_figure(report, replica, residue_limit=residue_limit)
    figure = Figure(figsize=(4.2, 3.0))
    axis = figure.subplots()
    plotted = False
    identifiers = [identifier for identifier in replica_ids(report) if selected_replicas is None or identifier in selected_replicas] if replica is None else [replica]
    if plot in TIME_PLOTS:
        _, units, fields = TIME_PLOTS[plot]
        for identifier in identifiers:
            row = _series(report, identifier)
            times = row.get("time_ns")
            if not isinstance(times, list):
                continue
            for field in fields:
                values = row.get(field)
                if not isinstance(values, list):
                    continue
                points = [(_number(time), _number(value)) for time, value in zip(times, values)]
                if not any(time is not None and value is not None for time, value in points):
                    continue
                axis.plot([time if time is not None else float("nan") for time, _ in points],
                    [value if value is not None else float("nan") for _, value in points], linewidth=.8,
                    label=(f"R{identifier}" if replica is None else field.replace("_angstrom", "").replace("_fraction", "").replace("_", " ")))
                plotted = True
        axis.set(xlabel="Time (ns)", ylabel=units)
        if plotted and (len(fields) > 1 or replica is None):
            axis.legend(fontsize=5, frameon=False)
    elif plot.endswith("_occupancy"):
        kind = plot.removesuffix("_occupancy")
        rows = contact_rows(report, replica, limit=residue_limit)
        pairs = [(row["residue"], row["occupancies"][kind] * 100) for row in rows if kind in row["occupancies"]]
        if pairs:
            axis.bar(range(len(pairs)), [pair[1] for pair in pairs], color=COLORS[kind], width=.65)
            axis.set_xticks(range(len(pairs)), [pair[0].replace(" · chain ", ":") for pair in pairs], rotation=65, ha="right", fontsize=5)
            axis.set(ylabel="Residue occupancy (%)", ylim=(0, 105))
            plotted = True
    elif plot in {"protein_rmsf", "ligand_rmsf"}:
        prefix = "protein" if plot == "protein_rmsf" else "ligand"
        for identifier in identifiers:
            row = _series(report, identifier)
            values = row.get(prefix + "_rmsf_angstrom")
            if isinstance(values, list) and any(_number(value) is not None for value in values):
                axis.plot([_number(value) if _number(value) is not None else float("nan") for value in values], linewidth=.8, label=f"R{identifier}")
                plotted = True
        axis.set(xlabel="Recorded residue index" if prefix == "protein" else "Recorded atom index", ylabel="RMSF (Å)")
        if replica is None and plotted:
            axis.legend(fontsize=5, frameon=False)
    elif plot in SUMMARY_PLOTS:
        _, units, field = SUMMARY_PLOTS[plot]
        rows = [row for row in _rows(report.get("replicas")) if row.get("replica") in identifiers and _number(row.get(field)) is not None
            and (not plot.startswith("endpoint_energy") or row.get("endpoint_status") in {"completed", "success"})]
        if rows:
            axis.bar([f"R{row['replica']}" for row in rows], [float(row[field]) for row in rows], color="#4c78a8", width=.5)
            axis.set_ylabel(units)
            plotted = True
    elif plot == "contact_matrix":
        import numpy as np

        rows = contact_rows(report, limit=residue_limit)
        matrix = [[_occupancy(_replica_values(report, "contact_occupancy", identifier).get(row["residue"])) for identifier in identifiers] for row in rows]
        if rows and identifiers:
            image = axis.imshow(np.array([[value * 100 if value is not None else float("nan") for value in row] for row in matrix]),
                aspect="auto", vmin=0, vmax=100, cmap="Blues")
            axis.set_xticks(range(len(identifiers)), [f"R{identifier}" for identifier in identifiers], fontsize=6)
            axis.set_yticks(range(len(rows)), [row["residue"].replace(" · chain ", ":") for row in rows], fontsize=5)
            figure.colorbar(image, ax=axis, label="Occupancy (%)", fraction=.06)
            plotted = True
    if not plotted:
        axis.text(.5, .5, "Not recorded for this replica", transform=axis.transAxes, ha="center", va="center", fontsize=8)
        axis.axis("off")
    axis.tick_params(labelsize=6)
    axis.xaxis.label.set_size(7)
    axis.yaxis.label.set_size(7)
    figure.tight_layout(pad=.7)
    return figure


def render_plot_matrix(st, report, *, key_prefix="md-overview"):
    available = available_plots(report)
    if not available or not report.get("replicas"):
        return False
    try:
        identifiers = replica_ids(report)
    except ValueError as exc:
        st.warning(str(exc))
        return False
    if report.get("success") is False or not identifiers:
        st.warning("The stored report has no confirmed completed replicas. Successful analysis is not asserted.")
        return False
    st.markdown("### Compact MD overview")
    st.caption("MN-Ligand-style plot matrix · saved analysis · no trajectory recalculation. Each column is one replica; compare the same plot across repeats.")
    with st.expander("Plot matrix options"):
        display = st.selectbox("Replica display", ["Separate", "Overlay / recorded consensus", "Both"], key=key_prefix + "-display")
        selected_replicas = st.multiselect("Replicas shown together", identifiers, default=identifiers, key=key_prefix + "-replicas")
        selected_replicas = selected_replicas if isinstance(selected_replicas, list) else identifiers
        defaults = [key for key in ("interaction_network", "backbone_rmsd", "ligand_rmsd", "contact_occupancy") if key in available]
        selected = st.multiselect("Plots", available, default=defaults or available[:4], format_func=lambda value: PLOT_LABELS[value], key=key_prefix + "-plots")
        selected = selected if isinstance(selected, list) else defaults
        limit = st.slider("Residues in interaction plots", 4, 16, 12, key=key_prefix + "-limit")
        st.caption("Interaction plots use the original saved consensus hotspot ranking (or contact occupancy when no hotspot score was saved). "
            "Every positive saved interaction for those residues is drawn, including occupancies below 10%. "
            "Map percentages are residue-level occupancies; named endpoints are the source report's representative ligand atoms, not measured atom-pair occupancies. "
            "Overlay shows separate recorded curves, not a newly calculated mean. Consensus bars/maps use the report's saved all-replica means, independent of the replica selection.")
    if not selected or not selected_replicas:
        st.info("Choose at least one plot and replica.")
        return True
    window = report.get("analysis_window")
    if isinstance(window, dict):
        st.caption("Recorded analysis window: " + " · ".join(f"{key.replace('_', ' ')}: {value}" for key, value in window.items() if not isinstance(value, (list, dict))))
    if st.button("Expand overview", key=key_prefix + "-expand-overview"):
        @st.dialog("MD overview · all selected replicas", width="large")
        def expanded_overview():
            st.caption("One saved analysis report · replica columns · interaction percentages on the lines. No recalculation.")
            _render_matrix_panels(st, report, selected, selected_replicas, display, limit, key_prefix + "-expanded", allow_enlarge=False)
        expanded_overview()
    _render_matrix_panels(st, report, selected, selected_replicas, display, limit, key_prefix)
    return True


def _render_matrix_panels(st, report, selected, selected_replicas, display, limit, key_prefix, *, allow_enlarge=True):
    for plot in selected:
        specs = ([None] if display != "Separate" else []) + (selected_replicas if display != "Overlay / recorded consensus" else [])
        if plot in SUMMARY_PLOTS or plot == "contact_matrix":
            specs = [None]
        st.markdown("**" + PLOT_LABELS[plot] + "**")
        columns = st.columns(min(3, len(specs)))
        for index, identifier in enumerate(specs):
            figure = plot_figure(report, plot, identifier, residue_limit=limit, selected_replicas=selected_replicas)
            with columns[index % len(columns)]:
                label = f"Replica {identifier}" if identifier is not None else "Recorded replicas / consensus"
                st.caption(label)
                st.pyplot(figure, clear_figure=False, width="stretch")
                panel_key = key_prefix + "-" + plot + "-" + str(identifier)
                if allow_enlarge and st.button("Enlarge", key=panel_key + "-enlarge"):
                    @st.dialog(PLOT_LABELS[plot] + " · " + label, width="large")
                    def enlarged():
                        st.pyplot(figure, clear_figure=False, width="stretch")
                    enlarged()
                with st.expander("Save figure"):
                    buffer = io.BytesIO()
                    figure.savefig(buffer, format="svg")
                    st.download_button("Download SVG", buffer.getvalue(), file_name=plot + "-" + str(identifier) + ".svg",
                        mime="image/svg+xml", key=panel_key + "-download")
