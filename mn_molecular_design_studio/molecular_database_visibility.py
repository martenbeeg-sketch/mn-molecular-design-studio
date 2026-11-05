"""Evidence-qualified presentation rules for historical Database collections."""
from __future__ import annotations


# Exact collection identities confirmed by the 2026-11-05 Database
# simplification audit. Workstation provenance and user ownership are separate
# fields; saved automatic memberships remain in All indexed.
UNREVIEWED_SCREENING_SHORTLISTS = {
    "screening-final-a5ba730c962dc449fb40bbab7702a25b": {
        "owner_id": "00000000000000000000000000000001",
        "run_id": "20261025-220608-e91da528",
        "source_node": "pc002975",
        "source_node_id": "50e484e4-0ad5-47c9-a8c4-b46d4000a59b",
    },
    "screening-final-7483c3d8846ffad1263b3b8d78078d3d": {
        "owner_id": "00000000000000000000000000000001",
        "run_id": "20261026-170220-29a4a53e",
        "source_node": "pc002975",
        "source_node_id": "50e484e4-0ad5-47c9-a8c4-b46d4000a59b",
    },
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
    """Match the audited source group or its owner-qualified peer replica."""
    if not isinstance(group, dict):
        return False
    origin = _origin(group)
    group_id = str(group.get("id") or "")
    replica_node = str(origin.get("_replica_node") or "").strip().casefold()
    if replica_node:
        # A peer marker is authoritative. Do not let a forged source_node
        # field turn a replica from another workstation into an audited match.
        original_group_id = str(origin.get("_replica_group") or "")
    else:
        # The source database records the workstation name separately from
        # both the user owner_id and the UUID used in peer replication.
        original_group_id = group_id
    identity = UNREVIEWED_SCREENING_SHORTLISTS.get(original_group_id)
    if not identity:
        return False
    if replica_node:
        if replica_node != identity["source_node_id"].casefold():
            return False
    elif str(origin.get("source_node") or "").strip().casefold() != identity["source_node"].casefold():
        return False
    return (
        _owner(group) == identity["owner_id"].casefold()
        and str(origin.get("run_id") or "") == identity["run_id"]
    )


def is_unreviewed_screening_shortlist_id(group_id) -> bool:
    """Return whether an exact audited collection ID needs historical navigation."""
    return str(group_id or "") in UNREVIEWED_SCREENING_SHORTLISTS


def is_confirmed_development_group(group) -> bool:
    """Match confirmed development executions; names and account flags are ignored."""
    if not isinstance(group, dict):
        return False
    return (_owner(group), str(_origin(group).get("run_id") or "")) in CONFIRMED_DEVELOPMENT_RUNS
