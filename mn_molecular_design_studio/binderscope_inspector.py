"""Studio-owned BinderScope inspector, not an embedded external report."""
from __future__ import annotations

import csv
import hashlib
import io
import json

import streamlit as st

from mn_molecular_design_studio.result_details_ui import render_result_details
from mn_molecular_design_studio.binderscope_view import render_linked_inspector
from mn_molecular_design_studio.binderscope_presentation import (
    contact_map, contact_network, contact_rows, highlight_selections, residue_label,
    residue_profile, residue_rows, result_coordinates, selected_contact,
)


def inspector_key(result, prefix=""):
    source = result.get("source") or {}
    identity = [prefix, result.get("candidate_id"), source.get("input_structure_digest"), source.get("assessment")]
    return "binderscope-inspector-" + hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]


def _csv(rows):
    output = io.StringIO()
    if rows:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output.getvalue()


def render_summary(result):
    analysis = result.get("analysis") or {}
    metrics = st.columns(4)
    for column, label, field, unit in zip(metrics, ("Buried area · both sides", "H-bonds", "Salt bridges", "Buried unpaired polar atoms"),
            ("dsasa", "n_hbonds", "n_salt", "n_unsat"), (" Å²", "", "", "")):
        value = analysis.get(field)
        column.metric(label, "Not available" if value is None else f"{value:,.0f}{unit}")
    st.caption("Post-analysis of this exact saved model, not a new filter decision or experimental binding evidence.")


def _render_context(result):
    hotspot = result.get("hotspot_compliance") or {}
    st.markdown("**Intended hotspot contacts**")
    if hotspot.get("requested_count"):
        st.write(f"{hotspot['contacted_count']} contacted · {hotspot['missed_count']} missed · {hotspot['unmapped_count']} unmapped "
                 f"of {hotspot['requested_count']} requested · {hotspot['cutoff_angstrom']} Å heavy-atom cutoff")
        st.dataframe([{"Requested residue": row["requested"], "Saved target residue": row.get("saved_residue", "—"),
                       "Result": row["status"], "Minimum distance (Å)": row.get("minimum_distance_angstrom"),
                       "Mapping note": row.get("reason", "")} for row in hotspot.get("residues") or []],
                     hide_index=True, use_container_width=True)
    else:
        st.caption("Not assessed: " + hotspot.get("status", "not recorded in this older assessment").replace("_", " ") + ". This is not a pass.")
    if hotspot.get("interpretation"):
        st.caption(hotspot["interpretation"])
    reference = result.get("target_reference_comparison") or {}
    if reference.get("status") == "assessed" and reference.get("target_ca_rmsd_angstrom") is not None:
        st.write(f"Original-target comparison: Cα RMSD {reference['target_ca_rmsd_angstrom']:.3f} Å "
                 f"over {reference['mapped_residues']} residues ({reference['scope']}).")
        st.caption(reference.get("interpretation", "Coordinate comparison, not experimental validation."))
    elif reference:
        st.caption("Original-target comparison: " + reference.get("reason", reference.get("status", "unavailable")))
    st.markdown("**Tags and fusions**")
    termini = result.get("termini") or {}
    if termini.get("status") == "unavailable" or termini.get("applies") is False:
        st.info(termini.get("reason") or termini.get("why") or "Free-terminal analysis does not apply to this topology.")
    else:
        rows = []
        for end in ("n", "c"):
            terminal = termini.get(end) or {}
            for tag, assessment in (terminal.get("g") or {}).items():
                rows.append({"End": end.upper(), "Tag / fusion": {"his6": "His6", "strep": "Twin-Strep", "domain": "GFP / SNAP / mCherry", "large": "Fc / MBP"}.get(tag, tag),
                    "Prediction": assessment.get("s"), "Linker residues": assessment.get("l"),
                    "Target distance (Å)": terminal.get("d"), "Confidence": terminal.get("conf"),
                    "Anchored through tail": terminal.get("via_tail")})
        if rows:
            st.dataframe(rows, hide_index=True, use_container_width=True)
        else:
            st.caption("No terminal assessment recorded.")
    st.caption("Tag/fusion assessments are geometric probes of one conformation, not simulations or guarantees. Missing environment context is not an all-clear.")
    for label, field in (("Binding-mode evidence", "binding_modes"), ("Specificity evidence", "specificity"), ("Scaffold / paratope evidence", "paratope")):
        with st.expander(label):
            render_result_details(st, result.get(field) or {"status": "Unavailable"})


def _render_structure(result, database, user, pair, key):
    source = result.get("source") or {}
    st.markdown("**Saved complex · 3D**")
    if database is None or user is None:
        st.info("Open this assessment from its saved design in Data Browser to view the verified complex.")
        return
    try:
        coordinates = result_coordinates(database, user, result)
    except (OSError, ValueError, PermissionError) as exc:
        st.warning("Structure unavailable: " + str(exc))
        return
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import ChainVisualization, StructureVisualization

    choices = ["Whole complex", "Interface", "Intended hotspots", "Buried polar atoms"]
    if pair:
        choices.insert(0, "Selected contact")
    layer = st.selectbox("Highlight in complex", choices, index=0 if pair else 1, key=key + "-layer-" + str(bool(pair)))
    color_choices = ["Chain role", "Hydrophobicity", "Secondary structure"]
    if result.get("analysis", {}).get("plddt_res") and source.get("assessment", {}).get("b_factor_interpretation") == "plddt_0_100":
        color_choices.insert(1, "pLDDT")
    binder_color = st.selectbox("Binder coloring", color_choices, key=key + "-binder-color")
    themes = {"Chain role": "uniform", "pLDDT": "plddt", "Hydrophobicity": "hydrophobicity", "Secondary structure": "secondary-structure"}
    chains = []
    for side, color, label in (("A", "0x477aa9", "Target"), ("B", "0x318c70", "Binder")):
        for chain in dict.fromkeys(row["chain"] for row in source.get("residue_map", {}).get(side) or []):
            chains.append(ChainVisualization(chain_id=chain, color=themes[binder_color] if side == "B" else "uniform",
                color_params={"value": color} if side == "A" or binder_color == "Chain role" else None,
                representation_type="cartoon", label=label))
    molstar_custom_component([StructureVisualization(pdb=coordinates, representation_type="cartoon", chains=chains or None,
        highlighted_selections=highlight_selections(result, layer, pair))], key=key + "-complex", height=510,
        show_controls=False, show_structure_tools=False)
    st.caption("Blue = target · binder uses the chosen coloring (green in chain-role mode). Highlighted residues use exact saved chain, author number and insertion code. Proteins remain cartoons.")
    if pair:
        st.caption(f"Selected pair: {pair['Target residue']} ↔ {pair['Binder residue']} · {pair['Interaction']} · {pair['Minimum distance (Å)']} Å")
    suffix = "cif" if "loop_" in coordinates else "pdb"
    st.download_button("Download assessed complex", coordinates, f"binderscope-assessed-complex.{suffix}",
                       "chemical/x-mmcif" if suffix == "cif" else "chemical/x-pdb", key=key + "-complex-download")


def _render_interface(result, database, user, key):
    rows = contact_rows(result)
    st.caption("Contacts: minimum heavy-atom distance ≤4.5 Å. H-bond geometry ≤3.5 Å; salt bridges use charged-group centroids ≤4.0 Å. Orange outlines mark mapped intended hotspots.")
    if not rows:
        st.info("No residue contacts were recorded for this model.")
        _render_structure(result, database, user, None, key)
        return
    diagram, structure = st.columns([1, 1.1], gap="large")
    with diagram:
        view = st.radio("2D interface", ["Contact map", "Residue network"], horizontal=True, key=key + "-diagram")
        interactions = st.multiselect("Interactions shown", list(dict.fromkeys(row["Interaction"] for row in rows)),
            default=list(dict.fromkeys(row["Interaction"] for row in rows)), key=key + "-interactions")
        shown = [row for row in rows if row["Interaction"] in interactions]
        if shown:
            figure = contact_map(result, shown) if view == "Contact map" else contact_network(result, shown)
            event = st.plotly_chart(figure, use_container_width=True, key=key + "-plot-" + view,
                on_select="rerun", selection_mode="points", config={"displaylogo": False, "scrollZoom": False})
            selected = selected_contact(event, shown)
        else:
            st.info("No interactions selected. The underlying assessment is unchanged.")
            selected = None
        st.caption(f"{len(shown)} of {len(rows)} residue pairs shown · marker size = atom-pair count. Click a marker or choose a pair below to highlight it in 3D.")
        if view == "Residue network":
            st.caption("Target left · binder right. Ordering reduces crossings; positions are not spatial coordinates. Click an edge's midpoint.")
        options = [None, *[(row["canonical_target"], row["canonical_binder"]) for row in shown]]
        lookup = {(row["canonical_target"], row["canonical_binder"]): row for row in shown}
        pair_key = key + "-selected-pair"
        if selected:
            identity = (selected["canonical_target"], selected["canonical_binder"])
            if st.session_state.get(key + "-plot-selected") != (view, identity):
                st.session_state[key + "-plot-selected"] = (view, identity)
                st.session_state[pair_key] = identity
        if st.session_state.get(pair_key) not in options:
            st.session_state[pair_key] = None
        chosen = st.selectbox("Inspect residue pair", options, key=pair_key, format_func=lambda identity:
            "No individual pair" if identity is None else lookup[identity]["Target residue"] + " ↔ " + lookup[identity]["Binder residue"])
    with structure:
        _render_structure(result, database, user, lookup.get(chosen), key)
    _render_contact_measurements(result, key)


def _render_contact_measurements(result, key):
    with st.expander("Contact measurements and CSV"):
        rows = contact_rows(result)
        public = [{field: value for field, value in row.items() if not field.startswith("canonical_")} for row in rows]
        st.dataframe(public, hide_index=True, use_container_width=True)
        st.download_button("Download all contacts (CSV)", _csv(public), "binderscope-contacts.csv", "text/csv", key=key + "-contacts-csv")


def _render_chemistry(result):
    analysis = result.get("analysis") or {}
    st.caption("Polar interactions are geometric assignments, not measured interaction energies. Hydrogen basis: " + str(analysis.get("hbond_basis", "unavailable")))
    for label, field in (("Hydrogen bonds", "hbonds"), ("Salt bridges", "salt")):
        st.markdown("**" + label + "**")
        rows = []
        for row in analysis.get(field) or []:
            rows.append({"Target": residue_label(result, "A", row["t"]), "Target amino acid": row.get("tn"), "Target atom": row.get("ta"),
                "Binder": residue_label(result, "B", row["b"]), "Binder amino acid": row.get("bn"), "Binder atom": row.get("ba"),
                "Distance (Å)": row.get("d"), "Direction": {"D>A": "Target donates", "A<D": "Target accepts"}.get(row.get("k"), "—")})
        if rows:
            st.dataframe(rows, hide_index=True, use_container_width=True)
        else:
            st.caption("None recorded." if field in analysis else "Not available.")
    st.markdown("**Buried unpaired polar atoms**")
    buried = [{"Side": "Target" if row["side"] == "t" else "Binder", "Residue": residue_label(result, "A" if row["side"] == "t" else "B", row["resi"]),
               "Amino acid": row.get("resn"), "Atom": row.get("atom")} for row in analysis.get("buns") or []]
    if buried:
        st.dataframe(buried, hide_index=True, use_container_width=True)
    else:
        st.caption("None recorded." if "buns" in analysis else "Not available.")
    st.markdown("**Exposed hydrophobic patches**")
    patches = [{"Patch": position + 1, "Exposed area (Å²)": row["area"],
                "Binder residues": ", ".join(residue_label(result, "B", number) for number in row["res"])}
               for position, row in enumerate(analysis.get("patches") or [])]
    if patches:
        st.dataframe(patches, hide_index=True, use_container_width=True)
    st.markdown("**Sequence liabilities**")
    liabilities = [{"Finding": row["kind"], "Binder residue": residue_label(result, "B", row["pos"]) if row.get("pos") else "Global",
                    "Motif": row.get("motif"), "Level": row.get("sev"), "Meaning": row.get("reason")}
                   for row in analysis.get("liab") or []]
    if liabilities:
        st.dataframe(liabilities, hide_index=True, use_container_width=True)
    else:
        st.caption("None recorded." if "liab" in analysis else "Not available.")
    if analysis.get("global_flags"):
        render_result_details(st, analysis["global_flags"])
    with st.expander("Sequence properties and secondary structure"):
        render_result_details(st, {field: analysis.get(field) for field in ("seq_props", "comp", "disulfides", "ssbond_basis", "segments")},
                              table_fields=("segments", "disulfides"))


def render_inspector(result, *, database=None, user=None, key=""):
    key = inspector_key(result, key)
    source = result.get("source") or {}
    label = source.get("source", {}).get("record_name") or source.get("label") or result.get("candidate_id") or "Saved design"
    st.markdown("#### " + label)
    structure_label = source.get("source", {}).get("structure_label") or "Exact assessed complex"
    binding = result.get("_studio_sql_binding") or {}
    if (structure_label == "Complex Pdb" and database is not None and user and binding.get("structure_id")
            and binding.get("structure_digest") and binding["structure_digest"] == source.get("input_structure_digest")):
        from mn_molecular_design_studio.structure_labels import ranked_bindcraft2_labels
        visible = database._visible_ids(user, scope="all", collection="all")
        with database.connection() as connection:
            labels = ranked_bindcraft2_labels(connection, [binding["structure_id"]], visible_groups=visible)
        structure_label = labels.get(binding["structure_id"], structure_label)
    st.caption(f"{structure_label} · saved author residue numbering")
    render_summary(result)
    analysis = result.get("analysis") or {}
    if analysis.get("confidence_status"):
        st.info(analysis["confidence_status"] + ". Geometry-based measurements remain available.")
    tabs = st.tabs(["Interface · 2D & 3D", "Residue profile", "Chemistry & liabilities", "Hotspots, target & tags", "Parameters & evidence"])
    with tabs[0]:
        if not render_linked_inspector(result, database, user, key=key):
            _render_interface(result, database, user, key)
        else:
            _render_contact_measurements(result, key)
    with tabs[1]:
        rows = residue_rows(result)
        if rows:
            st.plotly_chart(residue_profile(result), use_container_width=True, key=key + "-profile")
            public = [{field: value for field, value in row.items() if field != "canonical_number"} for row in rows]
            st.dataframe(public, hide_index=True, use_container_width=True)
            st.download_button("Download residue measurements (CSV)", _csv(public), "binderscope-residues.csv", "text/csv", key=key + "-residues-csv")
        else:
            st.caption("No mapped per-residue measurements are available.")
        st.caption("Confidence is shown only when confirmed for this exact binder model; target B-factors remain uninterpreted.")
    with tabs[2]:
        _render_chemistry(result)
    with tabs[3]:
        _render_context(result)
    with tabs[4]:
        render_result_details(st, source.get("assessment") or {"status": "No recorded per-design settings"})
        st.caption("Original admission, generator scores and model identities are retained separately in the design's score history. These views do not recompute or silently substitute another structure.")
        render_result_details(st, source)
        st.download_button("Download this assessment (JSON)", json.dumps(result, indent=2), "binderscope-assessment.json", "application/json", key=key + "-json")


@st.dialog("BinderScope · design inspector", width="large")
def inspection_dialog(result, database, user, key):
    render_inspector(result, database=database, user=user, key=key + "-dialog")
