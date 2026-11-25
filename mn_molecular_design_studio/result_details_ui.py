"""Bounded tables for supplementary run details with heterogeneous schemas."""

from collections.abc import Mapping
from itertools import islice
import math

_MAX_ROWS = 100
_MAX_COLUMNS = 20
_MAX_TEXT = 240
_METRIC_LABELS = {
    "pae": "PAE",
    "plddt": "pLDDT",
    "ptm": "pTM",
    "iptm": "ipTM",
    "ipsae": "ipSAE",
    "rmsd": "RMSD",
    "rmsf": "RMSF",
}


def result_field_label(field):
    text = str(field).replace("_", " ")
    words = text.split()
    if not words:
        return text
    words[0] = words[0][:1].upper() + words[0][1:]
    return " ".join(_METRIC_LABELS.get(word.casefold(), word) for word in words)


def result_value_text(value, *, summarize_fields=False):
    """Format a table cell, keeping large structures out of the browser payload."""
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, float):
        return str(value) if math.isfinite(value) else "—"
    if isinstance(value, Mapping):
        if summarize_fields and value:
            items = [f"{result_field_label(key)} = {result_value_text(item)}" for key, item in islice(value.items(), 6)]
            if len(value) > len(items):
                items.append(f"+{len(value) - len(items):,} more")
            return _short_text(" · ".join(items))
        return f"{len(value):,} fields" if value else "—"
    if isinstance(value, (list, tuple)):
        if len(value) <= 8 and all(isinstance(item, (str, int, float, bool)) for item in value):
            return _short_text(", ".join(result_value_text(item) for item in value)) or "—"
        return f"{len(value):,} entries"
    return _short_text(str(value))


def _short_text(text):
    return text if len(text) <= _MAX_TEXT else text[:_MAX_TEXT - 1] + "…"


def _field_rows(fields, path="", depth=0, omit=()):
    for field, value in fields.items():
        if field in omit or value is None or (isinstance(value, (str, list, tuple, Mapping)) and not value):
            continue
        label = result_field_label(field)
        label = f"{path} › {label}" if path else label
        if isinstance(value, Mapping) and depth < 3:
            yield from _field_rows(value, label, depth + 1)
        else:
            yield {"Field": label, "Value": result_value_text(value)}


def render_result_details(st, details, *, field="Value", table_fields=()):
    """Display scalar details and small record previews without serializing a payload."""
    if isinstance(details, Mapping):
        fields = _field_rows(details, omit=table_fields)
        rows = list(islice(fields, _MAX_ROWS + 1))
        if rows:
            st.dataframe(rows[:_MAX_ROWS], hide_index=True)
            if len(rows) > _MAX_ROWS:
                st.caption(f"Showing the first {_MAX_ROWS} fields.")
        for key in table_fields:
            if details.get(key):
                st.markdown(f"**{result_field_label(key)}**")
                render_result_details(st, details[key], field=result_field_label(key))
        return

    if isinstance(details, (list, tuple)):
        if not details:
            return
        visible = details[:_MAX_ROWS]
        if all(isinstance(row, Mapping) for row in visible):
            # Ordered union of the first available fields, bounded per record.
            columns = list(dict.fromkeys(
                key for row in visible for key in islice(row, _MAX_COLUMNS)
            ))
            st.dataframe(
                [{key: result_value_text(row.get(key)) for key in columns[:_MAX_COLUMNS]} for row in visible],
                hide_index=True,
                column_config={key: result_field_label(key) for key in columns[:_MAX_COLUMNS]},
            )
            if len(columns) > _MAX_COLUMNS or any(len(row) > _MAX_COLUMNS for row in visible):
                st.caption(f"Showing the first {_MAX_COLUMNS} fields per record.")
        elif len(details) <= 8 and all(isinstance(value, (str, int, float, bool)) for value in visible):
            st.dataframe([{"Field": field, "Value": result_value_text(details)}], hide_index=True)
        else:
            st.dataframe([{"Entry": index, "Value": result_value_text(value)}
                          for index, value in enumerate(visible, 1)], hide_index=True)
        if len(details) > _MAX_ROWS:
            st.caption(f"Showing {_MAX_ROWS} of {len(details):,} entries.")
    elif details is not None:
        st.dataframe([{"Field": field, "Value": result_value_text(details)}], hide_index=True)
