"""Bounded result tables, charts and portable reports without Streamlit state."""
from __future__ import annotations

import csv
import html
import io
import math

from mn_molecular_design_studio.result_details_ui import result_field_label, result_value_text


def candidate_state_rows(candidates):
    from mn_molecular_design_studio.binder_result_evidence import binder_decisions
    from mn_molecular_design_studio.output_vocabulary import candidate_evidence, format_candidate_parents
    from mn_molecular_design_studio.workflow_modes import native_ranking
    rows = []
    for candidate in candidates:
        current = candidate_evidence(candidate)[0]
        raw = candidate.get("raw_metadata") or {}
        selection = raw.get("selection") or {}
        collection = raw.get("candidate_set")
        collection = collection if isinstance(collection, dict) else {}
        source = raw.get("source_candidate")
        source = source if isinstance(source, dict) else {}
        source_candidate_id = (collection.get("native_candidate_id") or collection.get("native_backbone_identity")
                               or raw.get("native_candidate_id") or raw.get("pipeline_original_candidate_id")
                               or source.get("candidate_id") or candidate.get("candidate_id"))
        native_acceptance = binder_decisions(candidate)[0]["status"]
        validation = raw.get("validation_acceptance")
        validation = validation if isinstance(validation, dict) else {}
        shared_filter = collection.get("shared_filter")
        shared_filter = shared_filter if isinstance(shared_filter, dict) else {}
        measurement_detail = current["error"] or (
            "No independent measurement was recorded for this candidate."
            if current["status"] == "not_measured" else "")
        rows.append({"Candidate": candidate.get("candidate_id"), "Source candidate": source_candidate_id,
                     "Native generator acceptance": native_acceptance or "unavailable",
                     "Studio filtering": shared_filter.get("status", "not_evaluated"),
                     "Filter decision": selection.get("status", "not_evaluated"),
                     "Manual selection": bool((raw.get("human_review") or {}).get("selected", False)),
                     "Validation": validation.get("status", "not_evaluated"),
                     "Measurement": current["status"], "Measurement detail": measurement_detail,
                     "Selection": selection.get("status", "not_evaluated"), "Selected": selection.get("selected", False),
                     "Filter route": selection.get("routed_to", ""),
                     "Filter backend": (selection.get("filter_backend") or {}).get("status", "not_used"),
                     "Rank within protocol": selection.get("rank"), "Cause": current["error"],
                     "Selection reasons": "; ".join(selection.get("reasons") or []),
                     "Sequence": candidate.get("binder_sequence"), "Structure": candidate.get("complex_pdb"),
                     "Parents": format_candidate_parents(candidate.get("parents"))})
        from mn_molecular_design_studio.candidate_sets import _campaign_acceptance
        campaign = _campaign_acceptance(candidate)
        if campaign:
            rows[-1].update({
                "Campaign acceptance": campaign["status"],
                "Campaign acceptance level": campaign["acceptance_level"],
                "Campaign decision": campaign["reason"],
                "Campaign attempt": campaign["attempt_id"],
            })
        if collection:
            native_execution = (candidate.get('raw_metadata') or {}).get('distributed_campaign', {})
            predictor_execution = [{'predictor': item['branch'],
                'execution': item['candidate'].get('raw_metadata', {}).get('candidate_execution') or
                    item['candidate'].get('raw_metadata', {}).get('execution_provenance') or
                    item['candidate'].get('raw_metadata', {}).get('candidate_set', {}).get('prediction_execution', {})}
                for item in collection.get('predictions', [])]
            rows[-1].update({'Generator': collection.get('generator'), 'Source run': collection.get('source_run'),
                'Native attempt': collection.get('native_attempt_id'), 'Backbone': collection.get('backbone_id'),
                'Sequence variant': collection.get('variant_id'), 'Parent variants': ', '.join(collection.get('parent_variant_ids', [])),
                'Original sequence': collection.get('original_sequence'), 'Original structure': collection.get('original_structure'),
                'Native measurements (original source)': result_value_text(collection.get('native_measurements', {}), summarize_fields=True),
                'Native ranking': result_value_text(native_ranking(candidate), summarize_fields=True),
                'Workflow mode': (candidate.get('raw_metadata') or {}).get('pipeline_workflow', {}).get('effective_label'),
                'Comparative ranking measurement': selection.get('ranking_metric'),
                'Native filter': collection.get('native_filter'), 'Shared filter': collection.get('shared_filter', {}).get('status', 'not_evaluated'),
                'Predictor outcomes': result_value_text(collection.get('predictor_matrix', {}), summarize_fields=True),
                'Measurements': result_value_text(collection.get('analysis', {}).get('measurements', {}), summarize_fields=True),
                'Predictions': len(collection.get('predictions', [])), 'Diversity': result_value_text(collection.get('deduplication', {}), summarize_fields=True),
                'Worker / GPU provenance': result_value_text(
                    {'native': native_execution or collection.get('source_provenance', {}).get('metadata', {}).get('resource_allocation', {}),
                     'predictions': predictor_execution}, summarize_fields=True)})
    return rows


def convergence_table(payload, top_counts):
    return [{"round": row["round"], **summary} for row in payload.get("rounds", [])
            for summary in row.get("score_summaries", []) if summary.get("top_n") in top_counts]


def render_report(rows, *, format, title, x="round", y="mean_score", series="top_n", units=""):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    if format == "CSV":
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        # Keep spreadsheet applications from interpreting labels as formulas.
        writer.writerows({key: "'" + value if isinstance(value, str) and value.startswith(("=", "+", "-", "@")) else value
                          for key, value in row.items()} for row in rows)
        return stream.getvalue().encode(), "text/csv", "csv"
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_svg import FigureCanvasSVG
    figure = Figure(figsize=(9, 5), layout="constrained")
    axis = figure.subplots()
    groups = list(dict.fromkeys(row.get(series, "Results") for row in rows))
    for group in groups:
        selected = [row for row in rows if row.get(series, "Results") == group
                    and isinstance(row.get(x), (float, int)) and isinstance(row.get(y), (float, int))
                    and math.isfinite(row[x]) and math.isfinite(row[y])]
        if selected:
            axis.plot([r[x] for r in selected], [r[y] for r in selected], marker="o",
                      label=f"Top {group}" if series == "top_n" else str(group))
    axis.set(title=title, xlabel=x.replace("_", " ").capitalize(),
             ylabel=y.replace("_", " ").capitalize() + (" (" + units + ")" if units else ""))
    if axis.lines:
        axis.legend()
    axis.grid(alpha=0.2)
    stream = io.BytesIO()
    FigureCanvasSVG(figure).print_svg(stream)
    svg = stream.getvalue()
    if format == "SVG":
        return svg, "image/svg+xml", "svg"
    if format != "HTML":
        raise ValueError("Choose CSV, SVG or HTML.")
    table = "<table><thead><tr>" + "".join("<th>" + html.escape(key) + "</th>" for key in fields) + "</tr></thead><tbody>"
    table += "".join("<tr>" + "".join("<td>" + html.escape(str(row.get(key, ""))) + "</td>" for key in fields) + "</tr>" for row in rows)
    content = ("<!doctype html><html><head><meta charset='utf-8'><title>" + html.escape(title) +
        "</title><style>body{font:16px sans-serif;max-width:1100px;margin:2em auto}svg{max-width:100%;height:auto}"
        "td,th{padding:.5em;border-bottom:1px solid #ddd;text-align:left}table{border-collapse:collapse}</style></head><body><h1>" +
        html.escape(title) + "</h1>" + "<svg" + svg.decode().split("<svg", 1)[1] + table + "</tbody></table></body></html>")
    return content.encode(), "text/html", "html"


def render_workflow_results_summary(run_dir):
    """Normal Results panel, including successful runs with zero candidates."""
    from pathlib import Path
    import streamlit as st
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    summary = read_json(Path(run_dir) / 'result.json').get('outputs', {}).get('workflow_summary')
    if not isinstance(summary, dict) or not summary.get('preset_mode'):
        return
    with st.expander('Scientific workflow outcomes', expanded=True):
        st.write(summary['effective_label'])
        st.caption('Starting preset: ' + summary['preset_label'] + ' · Candidates: ' + str(summary['candidate_count']))
        st.caption(summary['ranking'])
        st.dataframe(summary['generator_operations'], hide_index=True)
        steps = {row["node"]: row.get("tool") or row["node"] for row in summary["execution_plan"]}
        counts = [{"Step": steps.get(step, step), "Outcome": result_field_label(outcome), "Count": count}
                  for step, outcomes in summary["stage_counts"].items() for outcome, count in outcomes.items()]
        if counts:
            st.dataframe(counts, hide_index=True)
        st.dataframe(summary['execution_plan'], hide_index=True)
