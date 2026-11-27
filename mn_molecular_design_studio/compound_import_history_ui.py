"""Review difficult compound rows and inspect saved decisions in Datasets."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import streamlit as st


def review_with_pubchem(invalid, key):
    from mn_molecular_design_studio.compound_import_history import pubchem_review
    from mn_molecular_design_studio.engine.workflow_engine.vendor.ligand_runtime.mn_ligand.workflows import compound_pubchem

    approved = st.session_state.setdefault(key + "-approved", {})
    with st.expander("Find a replacement in PubChem"):
        number = st.selectbox("Source row", list(invalid), format_func=lambda value: f"{value} · {invalid[value]['name']}", key=key + "-row")
        row = invalid[number]
        st.code(row["raw_smiles"] or "Blank SMILES", language=None, wrap_lines=True)
        fields = {label.casefold().strip(): str(value).strip() for label, value in (row.get("source_data") or {}).items()}
        cas = fields.get("cas no.") or fields.get("cas number") or fields.get("cas") or ""
        query_type = st.radio("Search by", ["CAS number", "Compound name"], index=0 if cas else 1,
                              horizontal=True, key=f"{key}-mode-{number}")
        search = st.text_input("PubChem query", value=cas if query_type == "CAS number" else row["name"],
                               key=f"{key}-query-{number}-{query_type}").strip()
        result_key = f"{key}-results-{number}"
        if st.button("Search PubChem", key=key + "-search", disabled=not search):
            try:
                with st.spinner("Searching PubChem…"):
                    hits = compound_pubchem.search_pubchem_candidates(search,
                        query_type="cas" if query_type == "CAS number" else "name")
                st.session_state[result_key] = hits
                st.session_state.pop(result_key + "-candidate", None)
                if not hits:
                    st.info("PubChem returned no matching structures.")
            except (ValueError, OSError) as exc:
                st.error(str(exc))
        if hits := st.session_state.get(result_key):
            choice = st.selectbox("PubChem structure", range(len(hits)), index=None, key=result_key + "-candidate",
                                  format_func=lambda index: f"CID {hits[index]['cid']} · {hits[index]['molecular_formula']}")
            if choice is not None and choice < len(hits):
                candidate = hits[choice]
                st.code(candidate["smiles"], language=None, wrap_lines=True)
                vendor_formula = fields.get("molecular formula") or fields.get("formula")
                if vendor_formula:
                    comparison = compound_pubchem.compare_molecular_formulas(vendor_formula, candidate["molecular_formula"])
                    st.caption(comparison["message"])
                    candidate = {**candidate, "vendor_formula": vendor_formula, "formula_comparison": comparison}
                components = compound_pubchem.pubchem_component_options(candidate["smiles"])
                selected = 0
                if len(components) > 1 or components[0]["occurrences"] > 1:
                    selected = st.selectbox("Component to import", range(len(components)), index=None,
                        format_func=lambda index: f"{components[index]['formula']} · {components[index]['smiles']}",
                        key=f"{result_key}-component-{candidate['cid']}")
                if selected is not None:
                    parent = compound_pubchem.selected_parent_candidate(candidate, components[selected])
                    from mn_molecular_design_studio.cofolding_input import _smiles_preview

                    if picture := _smiles_preview(parent["smiles"]):
                        st.image(picture)
                    if st.button("Accept selected structure", key=key + "-accept"):
                        if len(components) > 1 or components[0]["occurrences"] > 1:
                            parent["component_selection"] = "user-selected component"
                        approved[number] = pubchem_review(parent, source="Studio", query=candidate.get("query"),
                            query_type=candidate.get("query_type"), reviewed_at=datetime.now(timezone.utc).isoformat())
        if number in approved:
            st.caption(f"Accepted PubChem CID {approved[number]['candidate']['cid']} for this row.")
            if st.button("Clear accepted structure", key=key + "-clear"):
                approved.pop(number)
    return approved


def render_import_history(database, user, group, *, prefix):
    from mn_molecular_design_studio.database_ui import _import_report_frame, _report_csv

    batches = database.import_batches(user, group["id"])
    if not batches:
        return
    if not st.toggle("Show import history", key=prefix + "-import-history-" + group["id"]):
        return
    lookup = {row["id"]: row for row in batches}
    batch_id = st.selectbox("Import batch", list(lookup), key=prefix + "-import-batch-" + group["id"],
                            format_func=lambda value: f"{lookup[value]['source_name']} · {lookup[value]['created_at'][:19]}")
    batch = lookup[batch_id]
    summary = batch["summary"]
    rows = database.import_batch_rows(user, batch_id)
    reviewed = sum(bool(row.get("review") or row["descriptors"].get("corrected_smiles")) for row in rows)
    st.caption(f"{len(rows):,} source rows · {summary.get('unique_structures', 0):,} unique structures · {reviewed:,} reviewed replacements")
    view = st.selectbox("Show", ["Reviewed replacements", "All rows"] if reviewed else ["All rows"],
                        key=prefix + "-import-filter-" + batch_id)
    visible = [row for row in rows if row.get("review") or row["descriptors"].get("corrected_smiles")] if view == "Reviewed replacements" else rows
    table = _import_report_frame(visible[:1000])
    if view == "Reviewed replacements":
        columns = ("Row", "Compound ID", "Input SMILES", "Replacement SMILES", "Review method", "PubChem CID",
                   "Component selection", "Status", "Review run")
        table = [{key: row[key] for key in columns} for row in table]
    st.dataframe(table, hide_index=True, width="stretch")
    if len(visible) > 1000:
        st.caption(f"Showing the first 1,000 of {len(visible):,} rows.")
    st.download_button("Download validation report", lambda: _report_csv(rows),
                       file_name=Path(batch["source_name"]).stem + "-validation.csv",
                       mime="text/csv", key=prefix + "-import-download-" + batch_id)
    try:
        original = database.import_source_file(user, batch_id)
    except OSError:
        st.caption("Original upload is unavailable on this workstation.")
    else:
        st.download_button("Download original upload", original.read_bytes, file_name=batch["source_name"],
                           key=prefix + "-import-upload-" + batch_id)


def render_member_reviews(record):
    reviews = [source for source in record["sources"] if source.get("compound_review") or source.get("corrected_smiles")]
    if not reviews:
        return
    with st.expander("Compound import history"):
        for source in reviews:
            review = source.get("compound_review") or {}
            candidate = review.get("candidate") or {}
            label = f"{source.get('source_id') or record['name']} · {review.get('method') or 'Manual replacement'}"
            if cid := candidate.get("cid"):
                label += f" · PubChem CID {cid}"
            st.write(label)
            st.caption(f"{source.get('source_name')} · row {source.get('source_row')}")
            st.caption("Original SMILES")
            st.code(source.get("raw_smiles") or "Blank", language=None, wrap_lines=True)
            st.caption("Accepted SMILES")
            st.code(source.get("canonical_smiles") or source.get("corrected_smiles"), language=None, wrap_lines=True)
            if selection := candidate.get("component_selection"):
                st.caption(selection)
            if run := review.get("review_run_id"):
                st.caption(f"MN-Ligand review {run}")
