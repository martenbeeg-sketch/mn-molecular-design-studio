"""ACL-revalidated, paged Database structure choices for Complex Analysis UI."""
from __future__ import annotations

from typing import Any, Callable, MutableMapping


def authorized_structure_catalog(
    database,
    user: dict,
    *,
    session_state: MutableMapping[str, Any],
    cache_key: str,
    query: str = "",
    project_id: str = "",
    collection: str = "all",
    group_id: str = "",
    group_ids: list[str] | None = None,
    limit: int = 100,
    advance: bool = False,
    refresh: bool = False,
    eligible: Callable[[dict[str, Any]], bool] | None = None,
) -> dict[str, Any]:
    """Return a small fresh authorized page, rechecking cached IDs on reruns.

    Broad discovery is performed once per search/page. Streamlit reruns reuse
    only the IDs from pages already discovered, then ask the canonical
    ``list_structure_records`` API to reapply current ACLs and return fresh
    labels/representations. When ``eligible`` is supplied, raw Database pages
    are scanned only until the requested number of eligible rows is found. This
    prevents a page of unrelated structures from appearing empty when later
    Database pages contain valid complexes. No molecule names or source
    metadata are cached.
    """
    page_size = int(limit)
    if not 1 <= page_size <= 5000:
        raise ValueError("Structure picker page size must be between 1 and 5000.")
    user_id = str(user.get("id") or "")
    if not user_id:
        raise PermissionError("A signed-in user is required to list Database structures.")
    search = str(query or "")
    project = str(project_id or "")
    collection_name = str(collection or "all")
    group = str(group_id or "")
    selected_groups = list(dict.fromkeys(str(value) for value in (group_ids or []) if str(value))) if group_ids is not None else None
    filters = {"project_id": project, "collection": collection_name, "group_id": group,
               "group_ids": selected_groups}
    group_kwargs = ({"group_ids": selected_groups, "restrict_structure_groups": True}
                    if selected_groups is not None else {"group_id": group})
    cached = session_state.get(cache_key)
    reusable = (
        not refresh and isinstance(cached, dict)
        and str(cached.get("user_id") or "") == user_id
        and str(cached.get("query") or "") == search
        and cached.get("filters") == filters
        and int(cached.get("page_size") or 0) == page_size
        and bool(cached.get("eligibility_filter")) == bool(eligible)
    )

    def scan_eligible(start_offset: int, wanted: int) -> tuple[list[dict[str, Any]], int | None]:
        found: list[dict[str, Any]] = []
        offset = int(start_offset)
        next_offset: int | None = offset
        while len(found) < wanted and next_offset is not None:
            page = database.list_structure_records(
                user, collection=collection_name, project_id=project, **group_kwargs, query=search,
                limit=min(page_size, wanted - len(found)), offset=next_offset,
            )
            page_rows = [row for row in page.get("items", []) if isinstance(row, dict)]
            found.extend(row for row in page_rows if eligible is None or eligible(row))
            returned_offset = page.get("next_offset")
            next_offset = int(returned_offset) if returned_offset is not None else None
            if not page_rows:
                next_offset = None
        return found[:wanted], next_offset

    if eligible is not None:
        if not reusable:
            rows, next_offset = scan_eligible(0, page_size)
            ids = list(dict.fromkeys(str(row.get("id") or "") for row in rows if row.get("id")))
            cached = {
                "user_id": user_id, "query": search, "page_size": page_size,
                "filters": filters,
                "ids": ids, "authorized_ids": ids, "next_offset": next_offset,
                "eligibility_filter": True,
            }
            session_state[cache_key] = cached
            return {"records": rows, "next_offset": next_offset, "discovered": True,
                    "authorized_count": len(rows)}

        ids = list(dict.fromkeys(str(value) for value in cached.get("ids", []) if str(value)))
        next_offset = cached.get("next_offset")
        authorized_rows: list[dict[str, Any]] = []
        if ids:
            for start in range(0, len(ids), 5000):
                chunk = ids[start:start + 5000]
                page = database.list_structure_records(
                    user, collection=collection_name, project_id=project, **group_kwargs,
                    query=search, molecule_ids=chunk,
                    limit=max(1, len(chunk)), offset=0,
                )
                authorized_rows.extend(row for row in page.get("items", [])
                                       if isinstance(row, dict) and eligible(row))
        authorized_by_id = {str(row.get("id") or ""): row for row in authorized_rows if row.get("id")}
        authorized_ids = [value for value in ids if value in authorized_by_id]
        access_changed = authorized_ids != list(cached.get("authorized_ids", []))

        if advance and next_offset is not None and not access_changed:
            new_rows, next_offset = scan_eligible(int(next_offset), page_size)
            new_ids = [str(row.get("id") or "") for row in new_rows if row.get("id")]
            for row in new_rows:
                if row.get("id"):
                    authorized_by_id[str(row["id"])] = row
            ids = list(dict.fromkeys([*ids, *new_ids]))
        elif access_changed:
            # Offset pages can shift when access or eligibility changes. Rebuild
            # the loaded eligible prefix from the canonical ACL-filtered order.
            wanted = len(ids) + (page_size if advance and next_offset is not None else 0)
            rebuilt, next_offset = scan_eligible(0, wanted)
            ids = list(dict.fromkeys(str(row.get("id") or "") for row in rebuilt if row.get("id")))
            authorized_by_id = {str(row.get("id") or ""): row for row in rebuilt if row.get("id")}

        records = [authorized_by_id[value] for value in ids if value in authorized_by_id]
        cached = {**cached, "ids": ids, "authorized_ids": [str(row["id"]) for row in records],
                  "next_offset": next_offset}
        session_state[cache_key] = cached
        return {"records": records, "next_offset": next_offset,
                "discovered": bool(advance), "authorized_count": len(records)}

    if not reusable:
        page = database.list_structure_records(
            user, collection=collection_name, project_id=project, **group_kwargs,
            query=search, limit=page_size, offset=0,
        )
        ids = list(dict.fromkeys(str(row.get("id") or "") for row in page.get("items", []) if row.get("id")))
        cached = {
            "user_id": user_id, "query": search, "page_size": page_size,
            "filters": filters,
            "ids": ids, "authorized_ids": ids, "next_offset": page.get("next_offset"),
            "eligibility_filter": False,
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
                user, collection=collection_name, project_id=project, **group_kwargs, query=search,
                limit=min(5000, len(ids) + page_size), offset=0,
            )
            ids = [str(row.get("id") or "") for row in page.get("items", []) if row.get("id")]
        else:
            page = database.list_structure_records(
                user, collection=collection_name, project_id=project, **group_kwargs,
                query=search, limit=page_size, offset=int(next_offset),
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
            user, collection=collection_name, project_id=project, **group_kwargs,
            query=search, molecule_ids=chunk,
            limit=max(1, len(chunk)), offset=0,
        )
        rows_by_id.update({str(row.get("id") or ""): row for row in page.get("items", []) if row.get("id")})
    rows = [rows_by_id[value] for value in ids if value in rows_by_id]
    cached = session_state.get(cache_key)
    if isinstance(cached, dict):
        session_state[cache_key] = {**cached, "authorized_ids": [str(row.get("id") or "") for row in rows]}
    return {"records": rows, "next_offset": next_offset, "discovered": bool(advance and reusable),
            "authorized_count": len(rows)}
