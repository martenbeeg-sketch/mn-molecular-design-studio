"""One comparison view shared by the molecular library and job results."""
from __future__ import annotations

import hashlib
import json
import streamlit as st

from mn_molecular_design_studio.comparisons import compare, compatibility_warnings, prepare_rows
from mn_molecular_design_studio.output_vocabulary import PREFIX, DEFINITIONS, field_label


def render_comparison(rows, *, key, database=None, user=None, project_id="", selected_ids_key=None, settings=None):
    if not rows:
        st.info("No matching evaluations are available.")
        return
    settings = settings or {}
    result_count = len(rows)
    st.subheader("Compare candidates")
    st.info(
        "An evaluation is one tool's result for one candidate structure. Start by choosing the scores you care about, "
        "then set which direction is better. Add cutoffs only when you have a pass/fail requirement. "
        "Select rows in the table or points in the plot to keep or send candidates to another tool."
    )
    st.caption(
        f"This comparison currently has {result_count:,} evaluation rows. Each row is one selected tool result. "
        "A candidate gets a comparison rank only when all chosen scores are present and it meets every cutoff. "
        "With no cutoffs, every candidate says Pass; the rank still compares their scores."
    )
    with st.expander("Advanced: choose runs or combine results"):
        combine = st.checkbox("Combine evaluations of the same structure", value=settings.get("combine_same_subject", False), key=key + "-combine")
        runs = sorted({(row.get("evaluation") or row.get("selected_evaluation") or {}).get("evaluator", "") + " / " +
                       (row.get("evaluation") or row.get("selected_evaluation") or {}).get("run_id", "") for row in rows})
        saved_runs = settings.get("evaluator_runs")
        selected_runs = st.multiselect("Include tool runs", runs, default=runs if saved_runs is None else [run for run in saved_runs if run in runs], key=key + "-runs")
        st.caption("Use this when you want to compare a different run, or combine a prediction with later QC scores. Combining is allowed only when the same candidate, exact structure, target revision and chain roles were used. If two runs report conflicting values for one score, that score stays blank.")
    rows = prepare_rows(rows, combine_same_subject=combine, evaluator_runs=selected_runs)
    if not rows:
        st.info("Select at least one contributing evaluation run.")
        return
    fields = sorted({field for row in rows for field, value in row.get("scores", {}).items() if type(value) in {int, float}},
                    key=lambda field: (not field.startswith(PREFIX), field))
    st.caption(f"{len(rows):,} comparison rows. Every score retains its original evaluation provenance.")
    defaults = settings.get("fields") or [f for f in fields if f.startswith(PREFIX) and DEFINITIONS.get(f.removeprefix(PREFIX)) and DEFINITIONS[f.removeprefix(PREFIX)].direction][:2]
    selected = st.multiselect("Scores to compare", fields, default=[f for f in defaults if f in fields], format_func=field_label, key=key + "-fields")
    directions, weights, thresholds = {}, {}, []
    with st.expander("Ranking and acceptance rules"):
        methods = ["mean", "worst", "product", "pareto"]
        method = st.selectbox("How to combine the selected scores", methods, index=methods.index(settings.get("method", "mean")), key=key + "-method",
                              help="Mean gives each selected score equal influence after ranking. Worst uses the weakest score. Product penalizes a poor rank more strongly. Pareto shows candidates that cannot be beaten on every score at once.")
        for field in selected:
            definition = DEFINITIONS.get(field.removeprefix(PREFIX)) if field.startswith(PREFIX) else None
            help_text = f"{definition.scope}; {definition.method}; {definition.unit}" if definition else "Original tool value; choose its direction explicitly."
            columns = st.columns([3, 1.5, 1, 1.3, 1.3])
            columns[0].caption(field_label(field) + " · " + help_text)
            suffix = hashlib.sha1(field.encode()).hexdigest()[:10]
            direction = settings.get("directions", {}).get(field) or (definition.direction if definition else "higher")
            saved_test = next((t for t in settings.get("thresholds", []) if t["field"] == field), None)
            directions[field] = columns[1].selectbox("Better when", ["higher", "lower"], index=int(direction == "lower"), key=key + suffix + "-direction")
            weights[field] = columns[2].number_input("Weight", min_value=.01, value=float(settings.get("weights", {}).get(field, 1.)), key=key + suffix + "-weight")
            enabled = columns[3].checkbox("Threshold", value=bool(saved_test), key=key + suffix + "-enable")
            value = columns[4].number_input("Cutoff", value=float(saved_test["value"]) if saved_test else 0., key=key + suffix + "-cutoff")
            if enabled:
                operators = [">=", "<=", ">", "<", "="]
                operator = columns[3].selectbox("Rule", operators, index=operators.index(saved_test["operator"] if saved_test else ">=" if directions[field] == "higher" else "<="), key=key + suffix + "-operator")
                thresholds.append({"field": field, "operator": operator, "value": value})
    warnings = compatibility_warnings(rows, selected)
    for warning in warnings:
        st.warning(warning)
    ranking_fields = [] if warnings else selected
    try:
        frame = compare(rows, fields=ranking_fields, directions=directions, weights=weights, thresholds=thresholds, method=method)
    except ValueError as exc:
        st.error(str(exc))
        return
    st.caption(f"{int(frame['Pass'].sum()):,} / {len(frame):,} meet the cutoffs. With no cutoffs, all rows pass. Native ranks from the original tool are shown separately from this comparison rank.")
    visible = ["Compound ID", "Name", "Evaluator", "Target", "Run", "Native rank", "Comparison rank", "Pass", "Metric coverage", *selected, "Filter failures", "Missing metrics", "Conflicting metrics"]
    selection_version = hashlib.sha1(json.dumps([frame["Evaluation ID"].tolist(), selected, directions, weights, thresholds, method], sort_keys=True).encode()).hexdigest()[:12]
    event = st.dataframe(frame[visible].rename(columns={field: field_label(field) for field in selected}), hide_index=True,
                         on_select="rerun", selection_mode="multi-row", key=key + "-table-" + selection_version, use_container_width=True)
    selected_indices = event.selection.rows
    if len(selected) >= 2:
        import plotly.express as px
        xcol, ycol = st.columns(2)
        x = xcol.selectbox("X axis", selected, format_func=field_label, key=key + "-x")
        y = ycol.selectbox("Y axis", selected, index=1, format_func=field_label, key=key + "-y")
        color = st.selectbox("Color by", ["Evaluator", "Pass", "Run", "Target"], key=key + "-color")
        fig = px.scatter(frame, x=x, y=y, color=color, hover_name="Name", custom_data=["Molecule ID", "Evaluation ID"],
                         labels={x: field_label(x), y: field_label(y)})
        fig.update_layout(dragmode="lasso")
        plot = st.plotly_chart(fig, on_select="rerun", selection_mode=("points", "box", "lasso"), key=key + "-scatter-" + selection_version, use_container_width=True)
        points = plot.selection.points
        plot_pairs = {(p["customdata"][0], p["customdata"][1]) for p in points if p.get("customdata")}
        if plot_pairs:
            selected_indices = [i for i, row in frame.iterrows() if (row["Molecule ID"], row["Evaluation ID"]) in plot_pairs]
    with st.expander("Distributions and correlations"):
        if selected:
            import plotly.express as px
            field = st.selectbox("Distribution", selected, format_func=field_label, key=key + "-histfield")
            figure = px.histogram(frame, x=field, color="Evaluator", barmode="overlay", labels={field: field_label(field)})
            for threshold in thresholds:
                if threshold["field"] == field:
                    figure.add_vline(x=threshold["value"], line_dash="dash", annotation_text="Acceptance cutoff")
            st.plotly_chart(figure, key=key + "-hist", use_container_width=True)
            coverage = [{"Metric": field_label(f), "Available": int(frame[f].notna().sum()), "Missing": int(frame[f].isna().sum())} for f in selected]
            st.dataframe(coverage, hide_index=True)
            if thresholds:
                counts = [{"Rule": field_label(t["field"]) + " " + t["operator"] + " " + str(t["value"]),
                           "Passing": int(compare(rows, thresholds=[t])["Pass"].sum()), "Total": len(rows)} for t in thresholds]
                st.dataframe(counts, hide_index=True)
            if len(selected) > 1:
                numeric = frame[selected]
                st.plotly_chart(px.imshow(numeric.corr(method="spearman"), zmin=-1, zmax=1, color_continuous_scale="RdBu_r"), key=key + "-corr", use_container_width=True)
                st.caption("Spearman correlations use available pairs; missing values are excluded. Coverage is shown above.")
    selected_frame = frame.iloc[selected_indices] if selected_indices else frame.iloc[:0]
    if not selected_frame.empty:
        ids = list(dict.fromkeys(selected_frame["Molecule ID"]))
        st.caption(f"{len(ids)} molecules selected. Selection can be sent to another tile through Database.")
        if selected_ids_key and st.button("Use this selection", key=key + "-use"):
            st.session_state[selected_ids_key] = ids
            st.session_state["database-focus"] = ids[0]
            st.session_state["database-evaluation-" + ids[0]] = selected_frame.iloc[0]["Source evaluations"][0]
            st.rerun()
        if database and user:
            name = st.text_input("Save comparison as", value="Selected comparison", key=key + "-name")
            if st.button("Save comparison", key=key + "-save"):
                result = database.save_comparison(user, name, ids, project_id=project_id,
                    evaluation_filters={"evaluation_ids": sorted({identifier for ids in selected_frame["Source evaluations"] for identifier in ids})},
                    settings={"fields": selected, "directions": directions, "weights": weights, "thresholds": thresholds, "method": method,
                              "combine_same_subject": combine, "evaluator_runs": selected_runs})
                st.success("Saved the selected evaluations and comparison rules.")
                st.session_state["molecular-preselect-group"] = result["group_id"]
    st.download_button("Download comparison CSV", frame.to_csv(index=False).encode(), file_name="comparison.csv", mime="text/csv", key=key + "-csv")
