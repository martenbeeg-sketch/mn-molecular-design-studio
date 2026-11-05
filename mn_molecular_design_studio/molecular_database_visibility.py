"""Evidence-qualified presentation rules for historical Database collections."""
from __future__ import annotations


# Exact owner/run collections confirmed by the 2026-11-05 Database
# simplification audit. The saved automatic memberships remain in All indexed.
UNREVIEWED_SCREENING_SHORTLISTS = {
    "screening-final-a5ba730c962dc449fb40bbab7702a25b": ("pc002975", "20261025-220608-e91da528"),
    "screening-final-7483c3d8846ffad1263b3b8d78078d3d": ("pc002975", "20261026-170220-29a4a53e"),
}

# Only these eight executions are positively identified as bounded development
# or validation runs by the saved smoke reports and their original run IDs.
CONFIRMED_DEVELOPMENT_RUNS = {
    ("pc002975", "20261016-234215-9257efd9"),
    ("pc002975", "20261016-234216-3e724091"),
    ("pc002975", "20261017-000358-50e2305b"),
    ("pc002975", "20261017-000358-96b9e13c"),
    ("pc002975", "20261017-000359-1987191a"),
    ("pc002975", "20261017-000359-06035874"),
    ("pc002975", "20261017-000359-c16f32e2"),
    ("pc002975", "20261017-000359-e7fcc679"),
}


def _origin(group) -> dict:
    value = group.get("origin") if isinstance(group, dict) else None
    return value if isinstance(value, dict) else {}


def _owner(group) -> str:
    return str(group.get("owner_id") or "").strip().casefold() if isinstance(group, dict) else ""


def is_unreviewed_screening_shortlist(group) -> bool:
    """Match only an audited collection ID and its owner/run pair."""
    if not isinstance(group, dict):
        return False
    owner_run = UNREVIEWED_SCREENING_SHORTLISTS.get(str(group.get("id") or ""))
    return bool(owner_run and owner_run == (_owner(group), str(_origin(group).get("run_id") or "")))


def is_unreviewed_screening_shortlist_id(group_id) -> bool:
    """Return whether an exact audited collection ID needs historical navigation."""
    return str(group_id or "") in UNREVIEWED_SCREENING_SHORTLISTS


def is_confirmed_development_group(group) -> bool:
    """Match confirmed development executions; names and account flags are ignored."""
    if not isinstance(group, dict):
        return False
    return (_owner(group), str(_origin(group).get("run_id") or "")) in CONFIRMED_DEVELOPMENT_RUNS
