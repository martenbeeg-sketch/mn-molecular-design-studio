"""ACL-revalidated, paged Database structure choices for Complex Analysis UI."""
from __future__ import annotations

from typing import Any, Callable, MutableMapping


# Query keys are existing scientific references, not new system/study identities.
_RELATION_FIELDS = ("group_id", "molecule_id", "target_ref", "target_revision", "preparation_id", "study_id")


def _evidence_scope_query(database, user, group_ids, relationship=None, *, replica=None, include_history=True):
    import json

    ids = list(dict.fromkeys(map(str, group_ids)))
    database._authorize_group_ids_for_query(user, ids)
    if relationship is not None and (
        set(relationship) != set(_RELATION_FIELDS) or relationship["group_id"] not in ids
        or any(not isinstance(value, str) for value in relationship.values())
    ):
        raise ValueError("The selected relationship is outside these source memberships.")
    cte = """WITH evidence AS (
        SELECT o.*,
          COALESCE(NULLIF(json_extract(context,'$.target_molecule_id'),''),
                   NULLIF(json_extract(context,'$.target_id'),''),
                   CASE WHEN json_type(context,'$.target')='text' THEN json_extract(context,'$.target') END,'') AS target_ref,
          COALESCE(json_extract(context,'$.target_structure_digest'),'') AS target_revision,
          COALESCE(NULLIF(json_extract(context,'$.prepared_target_run_id'),''),
                   NULLIF(json_extract(context,'$.target_provenance.prepared_target_run_id'),''),
                   json_extract(context,'$.receptor_preparation_run_id'),'') AS preparation_id,
          COALESCE(json_extract(context,'$.study_id'),'') AS study_id,
          json_extract(context,'$.replica_index') AS replica,
          COALESCE(NULLIF(json_extract(context,'$.source_production_run_id'),''),
                   json_extract(context,'$.production_run_id'),'') AS production_run,
          COALESCE(json_extract(context,'$.continuation_of_run_id'),'') AS continuation_of,
          CASE WHEN json_type(context,'$.endpoint_superseded_by_run_id')='text'
               THEN COALESCE(json_extract(context,'$.endpoint_superseded_by_run_id'),'') ELSE '' END AS superseded_by
        FROM molecule_observations o
        WHERE group_id IN (SELECT value FROM json_each(?))
    ) """
    conditions, parameters = [], [json.dumps(ids)]
    if relationship is not None:
        for field in _RELATION_FIELDS:
            conditions.append(f"o.{field}=?")
            parameters.append(relationship[field])
    if replica is not None:
        conditions.append("o.replica=?")
        parameters.append(replica)
    if not include_history:
        conditions.append("o.superseded_by=''")
    return cte, " AND ".join(conditions) or "1", parameters


def saved_relationship_page(database, user: dict, group_ids: list[str], *, offset=0, limit=25, query="") -> dict:
    """Page molecule–target/preparation–study relationships from saved evidence.

    Membership scope prevents unrelated sources with missing study/target IDs
    from being merged. Neither labels nor shared normalized run IDs establish a
    relationship. Historical records contribute coverage, not new simulations.
    """
    if not 1 <= limit <= 100 or offset < 0:
        raise ValueError("Invalid relationship page.")
    cte, _, parameters = _evidence_scope_query(database, user, group_ids)
    fields = ",".join("o." + field for field in _RELATION_FIELDS)
    clause = "1"
    if query:
        clause = "(m.name LIKE ? OR o.target_ref LIKE ? OR o.study_id LIKE ? OR o.preparation_id LIKE ?)"
        parameters += ["%" + query + "%"] * 4
    grouped = ("SELECT " + fields + ",MAX(m.name) AS molecule_name, COUNT(*) AS record_count, "
        "COUNT(DISTINCT replica) AS replica_count, "
        "COUNT(DISTINCT NULLIF(production_run,'')) AS production_count, "
        "JSON_GROUP_ARRAY(DISTINCT namespace) AS methods, "
        "MAX(json_extract(o.context,'$.target_receptor')) AS target_name, "
        "MAX(json_extract(o.context,'$.target_pdb_id')) AS target_pdb, "
        "SUM(CASE WHEN superseded_by<>'' THEN 1 ELSE 0 END) AS historical_count "
        "FROM evidence o LEFT JOIN molecule_members m ON m.molecule_id=o.molecule_id AND m.group_id=o.group_id "
        "WHERE " + clause + " GROUP BY " + fields)
    with database.connection() as connection:
        totals = dict(connection.execute(cte + "SELECT COUNT(*) AS relationships, "
            "COUNT(DISTINCT molecule_id) AS molecules, "
            "COUNT(DISTINCT CASE WHEN target_ref<>'' OR target_revision<>'' OR preparation_id<>'' "
            "THEN json_array(target_ref,target_revision,preparation_id) END) AS target_contexts, "
            "COUNT(DISTINCT CASE WHEN study_id<>'' THEN json_array(group_id,study_id) END) AS studies, "
            "COALESCE(SUM(replica_count),0) AS replica_memberships FROM (" + grouped + ")", parameters).fetchone())
        rows = connection.execute(cte + grouped + " ORDER BY o.target_ref,o.study_id,o.molecule_id,o.group_id LIMIT ? OFFSET ?",
                                  [*parameters, limit, offset]).fetchall()
    import json
    items = [{**dict(row), "methods": json.loads(row["methods"]),
              "relationship": {key: str(row[key]) for key in _RELATION_FIELDS}} for row in rows]
    return {"items": items, "totals": totals, "next_offset": offset + limit if offset + limit < totals["relationships"] else None}


def relationship_structure_ids(database, user, group_ids, relationship):
    """Return only explicit representation/digest or exact target-context links.

    A shared molecule or normalized run ID alone cannot attach a model to a
    study. This does not fetch coordinates or infer receptor correspondence.
    """
    cte, clause, parameters = _evidence_scope_query(database, user, group_ids, relationship)
    with database.connection() as connection:
        rows = connection.execute(cte + "SELECT s.id FROM molecule_structures s "
            "WHERE s.group_id=? AND s.molecule_id=? AND (EXISTS (SELECT 1 FROM evidence o WHERE " + clause + " AND ("
            "NULLIF(json_extract(o.context,'$.representation_id'),'')=s.id OR "
            "COALESCE(NULLIF(json_extract(o.context,'$.source_structure_digest'),''),"
            "NULLIF(json_extract(o.context,'$.structure_digest'),''))="
            "COALESCE(NULLIF(json_extract(s.context,'$.portable_sha256'),''),"
            "NULLIF(json_extract(s.context,'$.content_sha256'),''),"
            "NULLIF(json_extract(s.context,'$.structure_digest'),'')))) OR "
            "(?='' AND ?<>'' AND "
            "COALESCE(NULLIF(json_extract(s.context,'$.target_molecule_id'),''),json_extract(s.context,'$.target_id'),'')=? "
            "AND COALESCE(json_extract(s.context,'$.target_structure_digest'),'')=? "
            "AND COALESCE(NULLIF(json_extract(s.context,'$.prepared_target_run_id'),''),"
            "json_extract(s.context,'$.target_provenance.prepared_target_run_id'),'')=?))",
            [parameters[0], relationship["group_id"], relationship["molecule_id"], *parameters[1:],
             relationship["study_id"], relationship["target_ref"], relationship["target_ref"],
             relationship["target_revision"], relationship["preparation_id"]]).fetchall()
    return {row["id"] for row in rows}


def saved_execution_groups(database, user: dict, group_ids: list[str], relationship: dict) -> list[dict]:
    """Show recorded replicas/segments and parent references without finalizing jobs."""
    cte, clause, parameters = _evidence_scope_query(database, user, group_ids, relationship)
    with database.connection() as connection:
        rows = connection.execute(cte + "SELECT replica,production_run, "
            "GROUP_CONCAT(DISTINCT NULLIF(continuation_of,'')) AS recorded_parents, "
            "COUNT(*) AS records, SUM(CASE WHEN superseded_by<>'' THEN 1 ELSE 0 END) AS historical_count "
            "FROM evidence o WHERE " + clause + " GROUP BY replica,production_run ORDER BY replica,production_run",
            parameters).fetchall()
    return [dict(row) for row in rows]


def retained_target_references(database, user: dict, group_ids: list[str]) -> dict[str, list[str]]:
    """Use existing Evaluation target links, without guessing biological identity."""
    import json

    ids = list(dict.fromkeys(map(str, group_ids)))
    database._authorize_group_ids_for_query(user, ids)
    with database.connection() as connection:
        if not connection.execute("SELECT 1 FROM sqlite_master WHERE name='molecule_evaluations' AND type='table'").fetchone():
            return {}
        rows = connection.execute(
            "SELECT DISTINCT group_id,target_id FROM molecule_evaluations "
            "WHERE group_id IN (SELECT value FROM json_each(?)) AND target_id<>'' ORDER BY target_id",
            (json.dumps(ids),),
        ).fetchall()
    result = {}
    for row in rows:
        result.setdefault(row["group_id"], []).append(row["target_id"])
    return result


def saved_evidence_summary(database, user: dict, group_ids: list[str], *, relationship=None,
                           replica=None, include_history=True) -> dict[str, dict]:
    """Count saved observation methods without hydrating models or reading Results."""
    import json

    ids = list(dict.fromkeys(map(str, group_ids)))
    if relationship is not None or replica is not None or not include_history:
        cte, clause, parameters = _evidence_scope_query(database, user, ids, relationship,
            replica=replica, include_history=include_history)
        sql = cte + "SELECT group_id,namespace,COUNT(*) AS n FROM evidence o WHERE " + clause + " GROUP BY group_id,namespace"
    else:
        database._authorize_group_ids_for_query(user, ids)
        sql = ("SELECT group_id,namespace,COUNT(*) AS n FROM molecule_observations "
               "WHERE group_id IN (SELECT value FROM json_each(?)) GROUP BY group_id,namespace")
        parameters = [json.dumps(ids)]
    with database.connection() as connection:
        rows = connection.execute(sql, parameters).fetchall()
    result = {key: {"observation_count": 0, "methods": {}} for key in ids}
    for row in rows:
        result[row["group_id"]]["methods"][row["namespace"]] = row["n"]
        result[row["group_id"]]["observation_count"] += row["n"]
    return result


def saved_evidence_page(database, user: dict, group_ids: list[str], *, namespace: str = "",
                        offset: int = 0, limit: int = 25, relationship=None,
                        replica=None, include_history=True) -> dict:
    """Page exact, authorized observations, including evidence without coordinates.

    No structure matching, metric flattening or identity creation occurs here.
    Evaluation descriptors use the same canonical projection as the Database.
    """
    import json
    from mn_molecular_design_studio.evaluations import evaluation_descriptor
    from mn_molecular_design_studio.output_vocabulary import normalize_metrics

    if not 1 <= limit <= 100 or offset < 0:
        raise ValueError("Invalid saved-evidence page.")
    cte, clause, parameters = _evidence_scope_query(database, user, group_ids, relationship,
        replica=replica, include_history=include_history)
    if namespace:
        clause += " AND o.namespace=?"
        parameters.append(namespace)
    with database.connection() as connection:
        has_evaluations = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE name='molecule_evaluations' AND type='table'").fetchone()
        total = connection.execute(cte + "SELECT COUNT(*) FROM evidence o WHERE " + clause,
                                   parameters).fetchone()[0]
        evaluation_select = ",e.descriptor AS saved_evaluation" if has_evaluations else ""
        evaluation_join = " LEFT JOIN molecule_evaluations e ON e.observation_id=o.id " if has_evaluations else " "
        rows = connection.execute(
            cte + "SELECT o.*,m.name AS molecule_name" + evaluation_select + " FROM evidence o "
            "LEFT JOIN molecule_members m ON m.group_id=o.group_id AND m.molecule_id=o.molecule_id "
            + evaluation_join + "WHERE " + clause + " ORDER BY o.created_at DESC,o.id LIMIT ? OFFSET ?",
            [*parameters, limit, offset],
        ).fetchall()
    items = []
    for row in rows:
        observation = dict(row)
        observation["metrics"] = json.loads(row["metrics"])
        observation["context"] = json.loads(row["context"])
        saved = observation.pop("saved_evaluation", None)
        observation["evaluation"] = json.loads(saved) if saved else evaluation_descriptor(observation)
        observation["normalized_output"] = normalize_metrics(observation["metrics"], row["namespace"], observation["context"])
        items.append(observation)
    return {"items": items, "total": total,
            "next_offset": offset + limit if offset + limit < total else None}


def collection_source_contexts(database, user: dict, groups: list[dict]) -> dict[str, dict]:
    """Read indexed source-run facts and saved input-target evidence, never Results.

    Native Evaluations may describe the original input target even when older
    group origins omitted it. Only the group's exact producing job is used;
    downstream analyses stored in that group cannot redefine its source target.
    """
    import json
    from mn_molecular_design_studio import scientific_lineage as lineage

    ids = [str(group["id"]) for group in groups if group.get("id")]
    if not ids:
        return {}
    database._authorize_group_ids_for_query(user, ids)
    result = {key: {"targets": [], "counts": []} for key in ids}
    sources = {}
    for group in groups:
        origin = group.get("origin") or {}
        source = origin.get("source_run") or origin
        sources[str(group["id"])] = {
            "run_id": str(source.get("run_id") or ""),
            "task_group": str(source.get("task_group") or ""),
            "node_id": str(origin.get("_replica_node") or source.get("source_node") or source.get("node_id") or ""),
            "workstation": str(source.get("source_workstation") or source.get("workstation") or ""),
            "project_id": str(group.get("project_id") or ""), "owner_id": str(group.get("owner_id") or ""),
        }
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT DISTINCT o.group_id, json_extract(o.context,'$.target_id') AS target_id, "
            "json_extract(o.context,'$.target_structure_digest') AS target_structure_digest, "
            "json_extract(o.context,'$.settings.target_name') AS target_name, "
            "json_extract(o.context,'$.settings.target_crop_lengths') AS crop_lengths, "
            "json_extract(o.context,'$.candidate.raw_metadata.target_sequences_by_chain') AS sequences, "
            "COALESCE(json_extract(o.context,'$._replica_node'),json_extract(o.context,'$.source_node')) AS node_id "
            "FROM molecule_observations o JOIN molecule_groups g ON g.id=o.group_id "
            "WHERE o.group_id IN (SELECT value FROM json_each(?)) "
            "AND o.run_id=COALESCE(json_extract(g.origin,'$.source_run.run_id'),json_extract(g.origin,'$.run_id')) "
            "AND o.namespace=COALESCE(json_extract(g.origin,'$.source_run.tool'),json_extract(g.origin,'$.tool'))",
            (json.dumps(ids),),
        ).fetchall()
        for row in rows:
            key = str(row["group_id"])
            source = sources[key]
            if source["node_id"] and row["node_id"] and source["node_id"] != row["node_id"]:
                continue
            if row["node_id"]:
                source["node_id"] = row["node_id"]
            target = {"label": row["target_name"] or "", "target_id": row["target_id"] or "",
                      "target_structure_digest": row["target_structure_digest"] or "",
                      "evidence": "saved native Evaluation input context"}
            if row["crop_lengths"]:
                target["crop_lengths"] = json.loads(row["crop_lengths"])
            if row["sequences"]:
                target["sequences_by_chain"] = json.loads(row["sequences"])
            if any(target.get(key) for key in ("label", "target_id", "target_structure_digest", "sequences_by_chain")):
                if target not in result[key]["targets"]:
                    result[key]["targets"].append(target)
        run_ids = sorted({source["run_id"] for source in sources.values() if source["run_id"]})
        clause, parameters = lineage._access(user, alias="j")
        jobs = connection.execute(
            "SELECT j.ref,j.job_metadata,j.facts,m.campaign_ref,c.name AS campaign_name,c.kind AS campaign_kind, "
            "(SELECT count(*) FROM scientific_campaign_jobs members WHERE members.campaign_ref=m.campaign_ref) AS campaign_jobs "
            "FROM scientific_jobs j LEFT JOIN scientific_campaign_jobs m ON m.job_ref=j.ref "
            "LEFT JOIN scientific_campaigns c ON c.ref=m.campaign_ref "
            "WHERE j.run_id IN (SELECT value FROM json_each(?)) AND " + clause,
            [json.dumps(run_ids), *parameters],
        ).fetchall() if run_ids else []
    jobs_by_scope = {}
    for row in jobs:
        job = json.loads(row["job_metadata"])
        scope = tuple(str(job.get(field) or "") for field in ("run_id", "task_group", "project_id", "owner_id"))
        jobs_by_scope.setdefault(scope, []).append((row, job))
    for key, source in sources.items():
        result[key]["source"] = source
        matching = []
        scope = tuple(source[field] for field in ("run_id", "task_group", "project_id", "owner_id"))
        for row, job in jobs_by_scope.get(scope, []):
            if source["node_id"]:
                if job.get("node_id") != source["node_id"]:
                    continue
            elif not source["workstation"] or str(job.get("workstation") or "").casefold() != source["workstation"].casefold():
                continue
            matching.append((row, job))
        if len(matching) != 1:
            continue  # Missing or ambiguous source facts remain unknown.
        row, job = matching[0]
        facts = json.loads(row["facts"])
        target = facts.get("target") or {}
        if target:
            result[key]["targets"] = target.get("recorded_constructs") or [target]
        # Every indexed standalone job also has a one-job campaign projection.
        # That projection must not erase Results' historical batch presentation.
        campaign_link = (row["campaign_ref"] or "") if row["campaign_jobs"] > 1 or any(
            job.get(field) for field in ("studio_run_id", "campaign_group_id", "distributed_campaign_id")) else ""
        result[key].update(job=job, job_ref=row["ref"], counts=facts.get("counts") or [],
                           acceptance=(facts.get("source_summary") or {}).get("acceptance_counts") or {},
                           campaign_ref=campaign_link, campaign_name=row["campaign_name"] or "",
                           campaign_kind=row["campaign_kind"] or "")
    return result


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
