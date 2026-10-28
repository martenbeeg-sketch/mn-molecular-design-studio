"""Bounded result tables, charts and portable reports without Streamlit state."""
from __future__ import annotations

import csv
import html
import io
import math


def candidate_state_rows(candidates):
    from mn_molecular_design_studio.output_vocabulary import candidate_evidence
    rows = []
    for candidate in candidates:
        current = candidate_evidence(candidate)[0]
        selection = (candidate.get("raw_metadata") or {}).get("selection") or {}
        rows.append({"Candidate": candidate.get("candidate_id"), "Measurement": current["status"],
                     "Selection": selection.get("status", "not_evaluated"), "Selected": selection.get("selected", False),
                     "Rank within protocol": selection.get("rank"), "Cause": current["error"],
                     "Selection reasons": "; ".join(selection.get("reasons") or []),
                     "Sequence": candidate.get("binder_sequence"), "Structure": candidate.get("complex_pdb"),
                     "Parents": ", ".join(candidate.get("parents") or [])})
        collection = (candidate.get('raw_metadata') or {}).get('candidate_set')
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
                'Native filter': collection.get('native_filter'), 'Shared filter': collection.get('shared_filter', {}).get('status', 'not_evaluated'),
                'Predictor outcomes': str(collection.get('predictor_matrix', {})),
                'Measurements': str(collection.get('analysis', {}).get('measurements', {})),
                'Predictions': len(collection.get('predictions', [])), 'Diversity': str(collection.get('deduplication', {})),
                'Worker / GPU provenance': str({'native': native_execution or collection.get('source_provenance', {}).get('metadata', {}).get('resource_allocation', {}),
                                               'predictions': predictor_execution})})
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
