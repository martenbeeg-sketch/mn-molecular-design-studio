"""Prepare external compound libraries without importing compounds into Database."""
from __future__ import annotations

from pathlib import Path
import secrets

import streamlit as st

from mn_molecular_design_studio.external_libraries import LibraryCatalog, source_columns
from mn_molecular_design_studio.engine.workflow_engine.vendor.ligand_runtime.mn_ligand.runtime import cpu_process_limit
from mn_molecular_design_studio.cluster import config as cluster_config


def _worker_choices():
    from mn_molecular_design_studio.cluster.store import Store

    organizer_id = cluster_config.identity()
    settings = cluster_config.load()
    local_limit = cpu_process_limit()
    by_id = {}
    if settings["mode"] != "disabled":
        try:
            for node in Store().nodes():
                snapshot = node.get("snapshot") or {}
                if node.get("disabled"):
                    continue
                resources = snapshot.get("resources") or {}
                total = int(resources.get("cpu_total") or 1)
                free = int(resources.get("cpu_free") or 0)
                node_id = str(node["id"])
                local = node_id == organizer_id
                # Peer snapshots can expire while metadata sync is running.
                # Keep registered enabled nodes selectable; the tile's submit
                # validator checks their live scheduler and free slots over SSH.
                # The organizer's selected CPU allocation also becomes the
                # local chunk worker's allocation, so expose its full pool;
                # remote workers can only lease slots currently free there.
                available = min(256, max(1, total if local else (free or total)))
                name = cluster_config.display_name(node.get("name"))
                by_id[node_id] = {"node_id": node_id, "name": name, "local": local,
                                  "cpu_available": available, "cpu_total": total,
                                  "availability": node.get("availability", "unknown")}
        except (OSError, ValueError, KeyError):
            pass
    if organizer_id not in by_id:
        name = cluster_config.display_name(settings.get("name") or "This workstation")
        by_id[organizer_id] = {"node_id": organizer_id, "name": name, "local": True,
                               "cpu_available": local_limit, "cpu_total": local_limit}
    return sorted(by_id.values(), key=lambda item: (not item["local"], item["name"].casefold()))


def render_setup(schema, *, user, project_selector, submit):
    st.title("Prepare library")
    st.write("Register a compound source and prepare Morgan fingerprints for the library. Similarity analysis and centroid generation are a separate tile.")
    st.caption("Input: a table with ID and SMILES (Name optional), or SDF. The Studio Database stores library metadata and file locations; compounds and fingerprints stay in external files.")
    project_id = project_selector("library-preparation")
    if user.get("role") != "admin":
        st.info("An administrator can register external library paths on this workstation.")
        return

    catalog = LibraryCatalog()
    records = catalog.list()
    choices = {row["id"]: row for row in records}
    selected = st.selectbox(
        "Library record", ["new", *choices],
        format_func=lambda key: "Register a new library" if key == "new" else choices[key]["name"] + " · " + key[-8:],
        key="library-prepare-choice",
    )
    record = choices.get(selected)
    prefix = "library-prep-" + selected
    nodes = _worker_choices()
    local_node = next(item for item in nodes if item["local"])
    if record and record.get("owner_node_id") and record["owner_node_id"] != local_node["node_id"]:
        st.info(f"This library record points to source files on workstation {record['owner_node_id']}. Select files that exist on this organizer, or open the tile on their owning workstation.")

    library_name = st.text_input("Library name", value=record["name"] if record else "", key=prefix + "-name")
    files_text = st.text_area(
        "Source file paths on this organizer",
        value="\n".join(item["path"] for item in record["sources"] if item.get("owner_node_id") in {None, local_node["node_id"]}) if record else "",
        help="One or more absolute paths on the organizer. Supported tables: CSV, TSV, Parquet, Feather, Arrow; SDF is also supported. Rows stream from disk.",
        key=prefix + "-files",
    )
    paths = [line.strip() for line in files_text.splitlines() if line.strip()]
    sdf_only = bool(paths) and all(Path(path).suffix.lower() == ".sdf" for path in paths)
    if st.button("Inspect source fields", disabled=not paths, key=prefix + "-inspect"):
        try:
            st.write({str(path): source_columns(Path(path).expanduser().resolve(strict=True)) for path in paths})
        except (OSError, ValueError) as exc:
            st.error(str(exc))

    fields = st.columns(3)
    id_column = fields[0].text_input(
        "ID column / SDF property", value=record["id_column"] if record else "id",
        help="Required for tables. For SDF, use an ID property when available, or leave blank to generate stable IDs from file and record order.",
        key=prefix + "-id",
    )
    name_column = fields[1].text_input(
        "Name column / SDF property (optional)", value=(record.get("name_column") or "") if record else "",
        help="Optional display name. Leave blank to use the compound ID as its label.",
        key=prefix + "-name-column",
    )
    smiles_column = fields[2].text_input(
        "SMILES column (blank for SDF)", value=(record.get("smiles_column") or "") if record else ("" if sdf_only else "smiles"),
        help="Required for tables. For SDF, leave blank and SMILES are derived from the molecular structure.",
        key=prefix + "-smiles",
    )
    output = st.text_input(
        "External project folder on this organizer",
        help="Choose the library project folder on the organizer. Studio writes the fingerprint Parquet table here and leaves the original library in place.",
        key=prefix + "-output",
    )

    default_node_ids = [local_node["node_id"]]
    chosen_ids = st.multiselect(
        "Worker computers",
        options=[item["node_id"] for item in nodes],
        default=default_node_ids,
        format_func=lambda node_id: next(item["name"] for item in nodes if item["node_id"] == node_id)
            + (" · organizer" if node_id == local_node["node_id"] else ""),
        help="One scheduled Studio job streams records to worker processes on these computers. Include the organizer because it owns the source library.",
        key=prefix + "-worker-nodes",
    )
    if local_node["node_id"] not in chosen_ids:
        st.warning("Include the organizer among the workers; it reads the source file and submits the single parent job.")
    worker_nodes = []
    local_workers = None
    selected_nodes = [item for item in nodes if item["node_id"] in chosen_ids]
    if selected_nodes:
        st.markdown("**CPU processes per computer**")
        columns = st.columns(min(3, len(selected_nodes)))
        for index, node in enumerate(selected_nodes):
            limit = node["cpu_available"]
            key = prefix + "-workers-v3-" + node["node_id"]
            count = int(columns[index % len(columns)].number_input(
                node["name"], min_value=1, max_value=limit, value=limit, step=1,
                help=(f"The scheduler reports {limit} available CPU slots. The default uses all of them; "
                      "lower this value to leave capacity for other work."), key=key,
            ))
            worker_nodes.append({"node_id": node["node_id"], "cpu_cores": count})
            if node["local"]:
                local_workers = count

    selection_mode = st.selectbox(
        "Record selection", ["first_rows", "random_sample"],
        format_func=lambda value: "First rows (0 = all rows)" if value == "first_rows"
        else "Random sample across all source rows",
        help="Random sampling scans the full source and selects the requested number of rows uniformly without replacement.",
        key=prefix + "-selection-mode",
    )
    max_records_text = st.text_input(
        "Random sample size" if selection_mode == "random_sample" else "Rows to process (0 = all rows)", value="0",
        help=("Choose the number of rows to sample from the entire source. The sample is reproducible with the seed below."
              if selection_mode == "random_sample" else
              "For a test, enter a smaller positive row count. Zero processes the full source with bounded memory and no fixed compound limit."),
        key=prefix + "-max-records",
    ).strip()
    try:
        max_records = int(max_records_text)
        max_records_valid = max_records >= 0
    except ValueError:
        max_records = 0
        max_records_valid = False
    if not max_records_valid:
        st.error("Enter a whole number of zero or greater.")
    sample_seed = int(st.number_input(
        "Random sample seed", min_value=0, max_value=2**31 - 1, value=42, step=1,
        help="The same seed and unchanged source select the same random rows.",
        key=prefix + "-sample-seed",
    )) if selection_mode == "random_sample" else 42
    if selection_mode == "random_sample" and max_records == 0:
        st.error("Set a positive random sample size.")
    ssh_user = st.text_input("SSH account on worker computers", value="compute", key=prefix + "-ssh-user")

    step_cols = st.columns(3)
    step_cols[0].markdown("**1 · Library source**  \nID + SMILES table or SDF")
    step_cols[1].markdown("**2 · Prepare fingerprints**  \nOpenVS Morgan · radius 2 · 1024 bits")
    step_cols[2].markdown("**3 · Compact library table**  \nID + optional Name + SMILES + packed fingerprint")
    st.caption("Valid records are written as Zstandard-compressed Parquet. Fingerprints use feature Morgan radius 2, 1024 bits, feature atom invariants and chirality. CPU processes stream bounded batches; invalid records are written separately. Limited runs retain their selected source rows as Parquet. Random mode makes a reproducible uniform sample across the full source without replacement.")
    st.info("One job is queued on the organizer's scheduler. Its worker processes run on the selected computers, and generated outputs are collected only after every worker finishes. The original library remains on the organizer.")

    ready = bool(project_id and library_name.strip() and paths and (id_column.strip() or sdf_only)
                 and output.strip() and local_workers and worker_nodes and ssh_user.strip() and max_records_valid
                 and not (selection_mode == "random_sample" and max_records == 0))
    if st.button("Queue fingerprint preparation", type="primary", disabled=not ready, key=prefix + "-submit"):
        parameters = {
            "source_files": paths, "library_name": library_name.strip(),
            "id_column": id_column.strip(), "name_column": name_column.strip(),
            "smiles_column": smiles_column.strip(), "output_directory": output.strip(),
            "cpu_cores": int(local_workers), "worker_nodes": worker_nodes,
            "max_records": max_records, "selection_mode": selection_mode,
            "sample_seed": sample_seed, "ssh_user": ssh_user.strip(),
        }
        try:
            submit(schema, "library.prepare", parameters, project_id,
                   {"_execution": {"placement": "local", "idempotency_key": secrets.token_urlsafe(24)}},
                   "Library fingerprint preparation")
        except (OSError, ValueError) as exc:
            st.error(str(exc))

    if record:
        with st.expander("Database record for this external library"):
            st.json(record)
