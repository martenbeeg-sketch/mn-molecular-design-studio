"""Choose coordinate files from the reusable input catalogue."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sqlite3

import streamlit as st

from mn_molecular_design_studio import datasets
from mn_molecular_design_studio.molecular_database import MolecularDatabase

COORDINATE_FORMATS = ("pdb", "ent", "cif", "mmcif")


@dataclass(frozen=True)
class DatasetFile:
    name: str
    path: Path
    molecule_id: str
    structure_id: str

    def getvalue(self) -> bytes:
        return self.path.read_bytes()


def options(database, user, *, project_id="", formats=COORDINATE_FORMATS, search=""):
    groups = {group["id"]: group["name"] for group in database.groups(user, collection="library", scope="all")
              if datasets.is_dataset(group) and (not group["project_id"] or
                  not project_id or project_id in datasets.project_references(group))}
    formats = sorted({str(value).lower().lstrip(".") for value in formats})
    if not groups or not formats:
        return []
    needle = "%" + search.strip().casefold() + "%"
    with database.connection() as connection:
        rows = connection.execute("""
            SELECT s.id, s.molecule_id, s.group_id, s.format, s.context, s.label, mm.name
            FROM molecule_structures s JOIN molecule_members mm ON mm.molecule_id=s.molecule_id AND mm.group_id=s.group_id
            JOIN molecule_groups g ON g.id=s.group_id
            WHERE s.group_id IN (SELECT value FROM json_each(?))
              AND lower(s.format) IN (SELECT value FROM json_each(?))
              AND lower(mm.name || ' ' || g.name || ' ' || s.label || ' ' ||
                  coalesce(json_extract(s.context,'$.source_filename'), '')) LIKE ?
            ORDER BY mm.name COLLATE NOCASE, g.name, s.label, s.id LIMIT 200
            """, (json.dumps(list(groups)), json.dumps(formats), needle)).fetchall()
    answer = []
    for row in rows:
        context = json.loads(row["context"] or "{}")
        name = Path(context.get("source_filename") or (row["name"] + "." + row["format"])).name
        if Path(name).suffix.lower().lstrip(".") != row["format"]:
            name = Path(name).stem + "." + row["format"]
        answer.append({"id": row["id"], "molecule_id": row["molecule_id"], "name": name,
                       "label": f"{row['name']} · {row['label']} · {groups[row['group_id']]} · {row['format'].upper()}"})
    return answer


def select_dataset_files(label, *, key, formats=COORDINATE_FORMATS, multiple=False, project_id="", help=None):
    from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
    from mn_molecular_design_studio.database_actions_ui import read_only_workspace
    from mn_molecular_design_studio.navigation import request_dataset_import

    if st.button("Import coordinates in Datasets", key=key + "-datasets", disabled=read_only_workspace()):
        request_dataset_import(st.session_state, "structure", project_id=project_id or "",
                               return_to=st.session_state.get("page", ""))
        st.rerun(scope="app")
    user = st.session_state.get("_molecular_database_user")
    if not user:
        return [] if multiple else None
    search = st.text_input("Find a saved structure", placeholder="Name or dataset", key=key + "-find")
    try:
        root = os.environ.get("MDS_DATA_BROWSER_DATABASE_ROOT", "").strip()
        database = ReadOnlyLibrary(Path(root) if root else None)
        available = {row["id"]: row for row in options(database, user, project_id=project_id, formats=formats, search=search)
                     } if database.path.is_file() else {}
        choice_key = key + "-selection-" + hashlib.sha256(repr((search, tuple(formats))).encode()).hexdigest()[:10]
        if multiple:
            selected = st.multiselect(label, list(available), format_func=lambda value: available[value]["label"],
                                      key=choice_key, help=help)
        else:
            selected = st.selectbox(label, list(available), index=None, format_func=lambda value: available[value]["label"],
                                    key=choice_key, help=help)
            selected = [selected] if selected else []
        files = []
        for identifier in selected:
            row = available[identifier]
            path = database.structure_file(user, identifier)
            files.append(DatasetFile(row["name"], path, row["molecule_id"], identifier))
        if not available:
            st.caption("No saved structures match.")
        return files if multiple else next(iter(files), None)
    except (OSError, ValueError, PermissionError, sqlite3.Error) as exc:
        st.error(f"Could not read saved structures: {exc}")
        return [] if multiple else None


def register_example(path, owner_id):
    """Register an explicitly loaded example once, retaining its exact coordinate bytes."""
    from mn_molecular_design_studio.database_actions_ui import read_only_workspace
    from mn_molecular_design_studio.molecule_imports import parse_structure

    source = Path(path)
    user = st.session_state.get("_molecular_database_user")
    if not user or str(user["id"]) != str(owner_id) or read_only_workspace():
        return source
    content = source.read_bytes()
    database = MolecularDatabase()
    group_id = "inputs-example-" + hashlib.sha256(str(owner_id).encode()).hexdigest()[:24]
    source_key = group_id + ":" + source.name + ":" + hashlib.sha256(content).hexdigest()
    with database.connection() as connection:
        exists = connection.execute("SELECT 1 FROM molecule_sources WHERE source_key=?", (source_key,)).fetchone()
    if exists:
        return source
    origin = datasets.input_origin(user)
    origin["bundled_examples"] = True
    records = parse_structure(content, source.name)
    records[0]["source_key"] = source_key
    group_id = database.create_group(user, "Bundled structure examples", group_id=group_id, origin=origin)
    registered = database.group(group_id, user)
    database.import_records(user, datasets.input_records(records), group_id=group_id, origin=registered["origin"])
    return source
