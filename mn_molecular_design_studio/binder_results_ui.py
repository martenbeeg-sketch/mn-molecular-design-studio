"""Binder and IDR decisions in the ordinary Results view."""
from __future__ import annotations

import math
import sqlite3

import streamlit as st

from mn_molecular_design_studio.binder_result_evidence import binder_outcomes, binder_results


_DATABASE_LABELS = {
    "registered": "Registered",
    "not_registered": "No saved structure matched",
    "different_saved_representation": "Different saved representation",
    "structure_unverified": "Structure unverified",
}


def _label(status):
    return str(status or "not_evaluated").replace("_", " ").capitalize()


def _assessment_choices(row):
    roles = {"refolded_validation": "Refolded model", "generated_design": "Generated model"}
    choices = [
        {**item, "relation": roles.get(item["representation_role"], label)}
        for field, label in (("binderscope", "Results model"), ("binderscope_derived", "Saved model"))
        for item in row[field]
    ]
    choices.sort(key=lambda item: item["observation_id"] != row["preferred_binderscope_observation_id"])
    return choices


def render_binder_outcomes(job, candidates, user, *, inputs=None, digests=None, database=None, key="binder"):
    if not binder_results(job, candidates, inputs):
        return False
    with st.expander("Binder design outcomes", expanded=True):
        if not candidates:
            st.info("No candidate records were published for this run.")
            return True
        if database is None:
            from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
            database = ReadOnlyLibrary()
        st.caption(f"{len(candidates):,} recorded candidates")
        page_size = 50
        pages = math.ceil(len(candidates) / page_size)
        page = st.selectbox("Candidate page", range(1, pages + 1), key=key + "-page") if pages > 1 else 1
        subset = candidates[(page - 1) * page_size:page * page_size]
        try:
            outcomes = binder_outcomes(database, user, job, subset, digests=digests)
        except (OSError, ValueError, sqlite3.Error) as exc:
            st.warning(f"Could not verify saved structure links: {exc}")
            return False
        assessments_by_row = [_assessment_choices(row) for row in outcomes]
        table = [{
            "Design": row["candidate_id"],
            "Native filter": _label(row["native_filter"]["status"]),
            "Final acceptance": _label(row["final_acceptance"]["status"]),
            "Database": _DATABASE_LABELS[row["database"]["status"]],
            "BinderScope": assessments[0]["relation"] if assessments else "—",
        } for row, assessments in zip(outcomes, assessments_by_row)]
        selection = st.dataframe(table, hide_index=True, use_container_width=True,
                                 on_select="rerun", selection_mode="single-row", key=f"{key}-table-{page}")
        indexes = selection.selection.rows if selection else []
        if not indexes or indexes[0] >= len(outcomes):
            return True
        row = outcomes[indexes[0]]
        reasons = [row["native_filter"]["reason"], row["final_acceptance"]["basis"]]
        reasons += [str(reason) for reason in row["final_acceptance"]["reasons"]]
        if explanation := " · ".join(text for text in reasons if text):
            st.caption(explanation)
        assessments = assessments_by_row[indexes[0]]
        if not assessments:
            return True
        selected = assessments[0] if len(assessments) == 1 else st.selectbox(
            "BinderScope assessment", assessments,
            format_func=lambda item: f"{item['relation']} · {item['model_label']} · {item['run_id']}",
            key=key + "-assessment-" + str(row["candidate_id"]),
        )
        st.caption(f"{selected['relation']} · {selected['model_label']} · {selected['run_id']}")
        from mn_molecular_design_studio.binderscope import saved_result
        from mn_molecular_design_studio.binderscope_ui import render_design_result

        linked = saved_result(database, user, selected)
        binding = (linked or {}).get("_studio_sql_binding") or {}
        if (not linked or binding.get("structure_id") != selected["structure_id"]
                or binding.get("structure_digest") != selected["sha256"]
                or (linked.get("source") or {}).get("input_structure_digest") != selected["sha256"]):
            st.caption("The saved assessment link is no longer verifiable.")
        else:
            render_design_result(linked, database=database, user=user, key=key + "-" + selected["observation_id"], compact=True)
    return True
