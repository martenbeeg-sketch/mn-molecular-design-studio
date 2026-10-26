"""Bounded result tables, charts and portable reports without Streamlit state."""
from __future__ import annotations

import csv
import html
import io
import math


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
