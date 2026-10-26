"""Campaign-wide results, independent of the underlying docking batch jobs."""
import math

import streamlit as st


def _finite(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def saturation_data(rounds):
    """Small, cumulative docking trajectories; never scan the compound table."""
    by_iteration = {}
    for row in rounds:
        number = _finite(row.get("round"))
        if number is None or number < 1 or not number.is_integer():
            continue
        scores = [_finite(item.get("score")) for item in row.get("top_scores", []) if isinstance(item, dict)]
        # A partial vector must not silently become a mean of a smaller top set.
        valid_vector = bool(scores) and all(score is not None for score in scores)
        by_iteration[int(number)] = {
            "iteration": int(number), "best_score": _finite(row.get("best_score")),
            "top_mean": math.fsum(scores) / len(scores) if valid_vector else None,
            "top_n": row.get("top_n") or len(scores),
            "maximum_change": _finite(row.get("maximum_rank_score_change")),
            "stable_comparisons": row.get("stable_round_count"),
        }
    return [by_iteration[number] for number in sorted(by_iteration)]


def render_saturation(payload, metrics):
    import altair as alt
    import pandas as pd

    rounds = payload.get("rounds") or []
    trajectory = saturation_data(rounds)
    summaries = {int(row["round"]): row.get("score_summaries") or [] for row in rounds
                 if _finite(row.get("round")) is not None}
    points, changes = [], []
    top_counts = set()
    for row in trajectory:
        if row["best_score"] is not None:
            points.append({"Iteration": row["iteration"], "Score": row["best_score"], "Series": "Best score", "Compounds": 1})
        groups = summaries.get(row["iteration"]) or [{"top_n": row["top_n"], "mean_score": row["top_mean"],
            "maximum_rank_score_change": row["maximum_change"], "stable_round_count": row["stable_comparisons"]}]
        for group in groups:
            count = int(group["top_n"])
            top_counts.add(count)
            label = f"Top {count:,}"
            mean = _finite(group.get("mean_score"))
            if mean is not None:
                points.append({"Iteration": row["iteration"], "Score": mean, "Series": label,
                               "Compounds": group.get("count", count)})
            change = _finite(group.get("maximum_rank_score_change"))
            if change is not None:
                changes.append({"Iteration": row["iteration"], "Maximum change": change, "Series": label,
                                "Stable comparisons": group.get("stable_round_count")})
    if not points:
        return
    units = payload.get("score_units") or metrics.get("score_units") or ""
    score_title = f"Docking score ({units})" if units else "Docking score"
    direction = payload.get("score_direction") or metrics.get("score_direction") or ""
    st.subheader("Score convergence")
    labels = [f"Top {number:,}" for number in sorted(top_counts)]
    colors = ["#d68116", "#38876b", "#8954a8", "#c54f64", "#606c38"]
    st.caption("Cumulative docking scores, ranked again after each iteration. "
               "Colored curves show the mean score of each top-ranked group; the blue curve shows the best compound. " + direction)
    chart = alt.Chart(pd.DataFrame(points)).mark_line(point=True).encode(
        x=alt.X("Iteration:Q", axis=alt.Axis(tickMinStep=1, format="d")),
        y=alt.Y("Score:Q", title=score_title, scale=alt.Scale(zero=False)),
        color=alt.Color("Series:N", title=None, scale=alt.Scale(
            domain=["Best score", *labels], range=["#1467ad", *colors[:len(labels)]])),
        tooltip=[alt.Tooltip("Iteration:Q", format="d"), "Series:N", alt.Tooltip("Score:Q", format=".5g"),
                 alt.Tooltip("Compounds:Q", format="d")],
    ).properties(height=280).interactive()
    st.altair_chart(chart, width="stretch")
    if not changes:
        return
    st.caption("Largest absolute score change at the same rank in the tracked top set, compared with the previous iteration. "
               "Stable comparisons also require unchanged compound order. "
               "A flat best-score line alone does not establish convergence of the top set.")
    change_title = f"Maximum rank score change ({units})" if units else "Maximum rank score change"
    change_chart = alt.Chart(pd.DataFrame(changes)).mark_line(point=True).encode(
        x=alt.X("Iteration:Q", axis=alt.Axis(tickMinStep=1, format="d")),
        y=alt.Y("Maximum change:Q", title=change_title, scale=alt.Scale(zero=True)),
        color=alt.Color("Series:N", title=None, scale=alt.Scale(domain=labels, range=colors[:len(labels)])),
        tooltip=[alt.Tooltip("Iteration:Q", format="d"), "Series:N", alt.Tooltip("Maximum change:Q", format=".5g"),
                 alt.Tooltip("Stable comparisons:Q", format="d")],
    )
    criteria = payload.get("stopping_criteria") or {}
    tolerance = _finite(criteria.get("score_tolerance"))
    if tolerance is not None:
        threshold = alt.Chart(pd.DataFrame([{"Tolerance": tolerance}])).mark_rule(
            strokeDash=[6, 4], color="#b04d3c").encode(
                y="Tolerance:Q", tooltip=[alt.Tooltip("Tolerance:Q", format=".5g")])
        change_chart = change_chart + threshold
    st.altair_chart(change_chart.properties(height=190).interactive(), width="stretch")
    stopping_top = int(criteria.get("top_n") or trajectory[-1]["top_n"])
    latest = max(row["Iteration"] for row in changes)
    last = next((row for row in changes if row["Iteration"] == latest and row["Series"] == f"Top {stopping_top:,}"), changes[-1])
    details = f"Run's stopping criterion: top {stopping_top:,} · Latest maximum change: {last['Maximum change']:g}" + (f" {units}" if units else "")
    if tolerance is not None:
        details += f" · Dashed line: tolerance {tolerance:g}"
    if last["Stable comparisons"] is not None:
        details += f" · Stable comparisons: {last['Stable comparisons']}"
        if criteria.get("stable_rounds"):
            details += f"/{criteria['stable_rounds']} required"
    st.caption(details)
    latest_groups = summaries.get(trajectory[-1]["iteration"]) or []
    if latest_groups:
        st.markdown("**Latest comparison by group**")
        st.dataframe([{"Group": f"Top {int(group['top_n']):,}", "Compounds": group.get("count"),
            "Mean score": group.get("mean_score"), "Score at cutoff": group.get("boundary_score"),
            "Maximum rank score change": group.get("maximum_rank_score_change"),
            "Same compound order": group.get("same_top_order"),
            "Stable comparisons": group.get("stable_round_count")} for group in latest_groups],
            hide_index=True, width="stretch")


def render_summary(payload, metrics, *, name, status, preview=True):
    terminal = str(status).casefold() in {"completed", "succeeded"}
    st.title("Final screening results" if terminal else "Screening results so far")
    st.caption(name)
    st.caption("Combined cumulative ranking across all iterations and all computers.")
    columns = st.columns(4)
    columns[0].metric("Iterations", metrics.get("round_count") or len((payload or {}).get("rounds", [])))
    columns[1].metric("Scored compounds", f"{int((payload or {}).get('compound_count') or metrics.get('unique_docked_compounds') or 0):,}")
    columns[2].metric("Failed compounds", f"{int(metrics.get('failed_compounds') or 0):,}")
    best = metrics.get("best_score")
    columns[3].metric("Best docking score", f"{float(best):g}" if best is not None else "—")
    if metrics.get("score_units"):
        st.caption(f"Score: {metrics.get('score_field') or 'docking score'} · {metrics['score_units']} · {metrics.get('score_direction') or ''}")
    reason = metrics.get("stop_reason") or (payload or {}).get("stop_reason")
    if reason:
        st.caption("Stopped: " + str(reason).replace("_", " "))
    if not payload:
        st.info("The combined score table has not been published yet. Iteration and batch jobs remain available under the campaign's Steps action.")
        return
    rounds = payload.get("rounds") or []
    if rounds:
        render_saturation(payload, metrics)
        st.subheader("Iteration history")
        st.dataframe([{"Iteration": row.get("round"), "Docked": row.get("successful_count", row.get("docked_count")),
            "Failed": row.get("failed_count"), "Cumulative scored": row.get("unique_tested_count"),
            "Best score": row.get("best_score"), "Binder %": row.get("binder_percent"),
            "Stable comparisons": row.get("stable_round_count")} for row in rounds], hide_index=True, width="stretch")
    if preview:
        st.subheader("Combined compound ranking")
        columns = {"rank", "compound_id", "compound_name", "smiles", "score", "score_units", "first_round", "best_round"}
        st.dataframe([{k: v for k, v in row.items() if k in columns} for row in payload.get("rows", [])],
                     hide_index=True, width="stretch")
        st.caption(f"Preview of {len(payload.get('rows', [])):,} compounds. All {int(payload['compound_count']):,} scored compounds remain in the external Parquet table.")
        st.code(payload["scores_path"], language=None)
        for name, path in (payload.get("files") or {}).items():
            st.caption(f"{name}: {path}")


def render_local_results(run_dir, user, result, *, name, status):
    from mn_molecular_design_studio.bulk_screening_results import external_results
    from mn_molecular_design_studio.molecular_database_results_ui import _render_external_screening_results
    payload = external_results(run_dir, user)
    render_summary(payload, result.get("metrics") or {}, name=name, status=status, preview=False)
    if payload:
        _render_external_screening_results(run_dir, user)
