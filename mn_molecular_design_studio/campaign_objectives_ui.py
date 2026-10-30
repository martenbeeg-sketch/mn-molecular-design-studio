"""Common duration/deadline controls for existing campaign scopes."""
from __future__ import annotations

from datetime import datetime
from zoneinfo import ZoneInfo

from mn_molecular_design_studio.campaign_objectives import timestamp


def render_time_limits(st, current, *, key, required=False):
    duration = float(st.number_input(
        "Run duration (minutes; 0 = no duration limit)", min_value=0.0,
        value=float(current.get("max_runtime_seconds") or (3600 if required and not current.get("deadline_at") else 0)) / 60,
        key=key + "-duration-minutes",
        help="Measured from the campaign's original start; recovery does not reset this timer.",
    ))
    deadline = current.get("deadline_at")
    try:
        displayed = timestamp(deadline).astimezone(ZoneInfo("Europe/Rome")).isoformat(timespec="minutes") if deadline else ""
    except ValueError:
        displayed = str(deadline)  # An invalid draft must remain editable.
    value = st.text_input(
        "Stop submitting at (Rome date and time, optional)", value=displayed,
        placeholder="YYYY-MM-DD HH:MM", key=key + "-deadline",
        help="For example 2026-10-31 08:00. An explicit timezone offset is also accepted. Active jobs finish and retain Results.",
    ).strip()
    limits = {"max_runtime_seconds": duration * 60} if duration else {}
    error = None
    if value:
        try:
            parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
            if parsed.tzinfo is None:
                zone = ZoneInfo("Europe/Rome")
                local = parsed.replace(tzinfo=zone)
                if local.utcoffset() != local.replace(fold=1).utcoffset():
                    raise ValueError("This Rome time is ambiguous; supply its explicit +01:00 or +02:00 offset.")
                parsed = local
            limits["deadline_at"] = timestamp(parsed.isoformat()).isoformat()
        except ValueError as exc:
            limits["deadline_at"] = value  # Preserve invalid input so graph/API validation cannot silently drop it.
            error = f"Invalid campaign deadline: {exc}"
    if required and not limits and not error:
        error = "Set a duration or deadline for a timed campaign."
    if error:
        st.error(error)
    st.caption("The first configured target or limit reached stops new dispatch. Active work finishes; starting another batch after campaign completion requires researcher review.")
    return limits, error is None
