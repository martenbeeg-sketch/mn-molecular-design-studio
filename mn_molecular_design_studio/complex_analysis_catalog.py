"""ACL-revalidated, paged Database structure choices for Complex Analysis UI."""
from __future__ import annotations

from typing import Any, MutableMapping


def authorized_structure_catalog(
    database,
    user: dict,
    *,
    session_state: MutableMapping[str, Any],
    cache_key: str,
    query: str = "",
    limit: int = 100,
    advance: bool = False,
    refresh: bool = False,
) -> dict[str, Any]:
    """Return a small fresh authorized page, rechecking cached IDs on reruns.

    Broad discovery is performed once per search/page. Streamlit reruns reuse
    only the IDs from pages already discovered, then ask the canonical
    ``list_structure_records`` API to reapply current ACLs and return fresh
    labels/representations. No molecule names or source metadata are cached.
    """
    page_size = int(limit)
    if not 1 <= page_size <= 5000:
        raise ValueError("Structure picker page size must be between 1 and 5000.")
    user_id = str(user.get("id") or "")
    if not user_id:
        raise PermissionError("A signed-in user is required to list Database structures.")
    search = str(query or "")
    cached = session_state.get(cache_key)
    reusable = (
        not refresh and isinstance(cached, dict)
        and str(cached.get("user_id") or "") == user_id
        and str(cached.get("query") or "") == search
        and int(cached.get("page_size") or 0) == page_size
    )

    if not reusable:
        page = database.list_structure_records(
            user, collection="all", query=search, limit=page_size, offset=0,
        )
        ids = list(dict.fromkeys(str(row.get("id") or "") for row in page.get("items", []) if row.get("id")))
        cached = {
            "user_id": user_id, "query": search, "page_size": page_size,
            "ids": ids, "authorized_ids": ids, "next_offset": page.get("next_offset"),
        }
        session_state[cache_key] = cached
        rows = list(page.get("items") or [])
        return {"records": rows, "next_offset": cached["next_offset"], "discovered": True,
                "authorized_count": len(rows)}

    ids = list(dict.fromkeys(str(value) for value in cached.get("ids", []) if str(value)))
    next_offset = cached.get("next_offset")
    if advance and next_offset is not None:
        previously_authorized = set(str(value) for value in cached.get("authorized_ids", []))
        if previously_authorized != set(ids):
            # ACL changes can shift offset pagination. Rebuild the loaded prefix
            # from the canonical current listing. Use the previously loaded
            # ID count (including revoked entries) plus one page, so revocation
            # does not silently shorten the newly requested page.
            page = database.list_structure_records(
                user, collection="all", query=search,
                limit=min(5000, len(ids) + page_size), offset=0,
            )
            ids = [str(row.get("id") or "") for row in page.get("items", []) if row.get("id")]
        else:
            page = database.list_structure_records(
                user, collection="all", query=search, limit=page_size, offset=int(next_offset),
            )
            new_ids = [str(row.get("id") or "") for row in page.get("items", []) if row.get("id")]
            ids = list(dict.fromkeys([*ids, *new_ids]))
        next_offset = page.get("next_offset")
        cached = {**cached, "ids": ids, "next_offset": next_offset}
        session_state[cache_key] = cached

    # ACLs are re-evaluated on every rerun, including cached picker pages.
    # Chunk to the Database API's supported page size if a user loads many pages.
    rows_by_id: dict[str, dict] = {}
    for start in range(0, len(ids), 5000):
        chunk = ids[start:start + 5000]
        page = database.list_structure_records(
            user, collection="all", query=search, molecule_ids=chunk,
            limit=max(1, len(chunk)), offset=0,
        )
        rows_by_id.update({str(row.get("id") or ""): row for row in page.get("items", []) if row.get("id")})
    rows = [rows_by_id[value] for value in ids if value in rows_by_id]
    cached = session_state.get(cache_key)
    if isinstance(cached, dict):
        session_state[cache_key] = {**cached, "authorized_ids": [str(row.get("id") or "") for row in rows]}
    return {"records": rows, "next_offset": next_offset, "discovered": bool(advance and reusable),
            "authorized_count": len(rows)}
