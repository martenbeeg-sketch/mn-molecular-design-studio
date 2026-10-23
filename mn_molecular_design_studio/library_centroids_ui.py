"""Studio tile for creating a separately cataloged centroid set."""
from __future__ import annotations

import secrets

import streamlit as st

from mn_molecular_design_studio.external_libraries import LibraryCatalog, FEATURE_SCHEMA
from mn_molecular_design_studio.cluster import config as cluster_config
from mn_molecular_design_studio.library_preparation_ui import _worker_choices
from mn_molecular_design_studio.centroid_chunks import balanced_chunk_sizes


def render_setup(schema, *, user, project_selector, submit):
    st.title("Create centroids")
    st.write("Cluster a prepared fingerprint library in random, balanced chunks. CPU workers on the selected computers select representative compounds; Studio stores the final compressed table externally and registers its metadata.")
    project_id = project_selector("library-centroids")
    if user.get("role") != "admin":
        st.info("An administrator can select external library paths and create centroid preparations.")
        return

    catalog = LibraryCatalog()
    records = catalog.list()
    choices = {item["id"]: item for item in records}
    if not choices:
        st.info("Register a library and prepare fingerprints first.")
        return
    library_id = st.selectbox(
        "External library", list(choices),
        format_func=lambda key: choices[key]["name"] + " · " + key[-8:], key="library-centroid-library",
    )
    record = choices[library_id]
    settings = cluster_config.load()
    placement = st.session_state.get("cluster-placement", settings["default_placement"])
    if not cluster_config.enabled():
        placement = "local"
    if placement == "auto":
        placement = "local"
        st.caption("Automatic placement stays on this workstation because the fingerprint files are stored on their owning workstation.")
    remote_target = placement not in {"local", None, ""}
    selected_node = placement if remote_target else cluster_config.identity()
    wrong_owner = bool(record.get("owner_node_id") and record["owner_node_id"] != selected_node)
    if wrong_owner:
        st.info(f"This library's files are registered on workstation {record['owner_node_id']}. Select that workstation in the sidebar to use its local paths.")
    prepared = [item for item in record.get("preparations", [])
                if item.get("schema_version") == FEATURE_SCHEMA
                and item.get("owner_node_id") in {None, selected_node}]
    if not prepared:
        st.info("This library has no registered fingerprint preparation yet.")
        return
    manifest_index = st.selectbox(
        "Prepared fingerprint table", range(len(prepared)),
        format_func=lambda index: f"{prepared[index].get('rows', 0):,} compounds · {prepared[index].get('created_at', 'unknown date')}",
        key="library-centroid-manifest",
    )
    selected_manifest = prepared[manifest_index]

    if remote_target:
        st.info("Centroid preparation reads the fingerprint files on the selected compute workstation. After completion, Studio copies the compressed centroid table and manifest to the organizer; the compute-side files remain in place.")

    threshold = st.number_input("Tanimoto threshold", min_value=0.01, max_value=1.0, value=0.6, step=0.01, key="library-centroid-threshold")
    methods = {
        "sampled": "Random chunks → sample representatives → final clustering",
        "hierarchical": "Random chunks → merge all representatives in further rounds",
        "exact": "Global exact clustering (serial reference method)",
    }
    method = st.selectbox("Clustering method", list(methods), format_func=methods.get,
                          key="library-centroid-method")
    st.caption("Fingerprints stay unchanged: feature Morgan, radius 2, 1024 bits, chirality enabled. Each chunk uses OpenVS-style best-first selection, ordered by molecular weight, with exact Tanimoto comparisons and no pairwise matrix.")
    row_count = int(selected_manifest.get("rows") or 0)
    chunk_size, chunk_count, candidate_count, merge_rounds, seed = 300_000, 0, 300_000, 3, 42
    worker_nodes, source_cpus = [], 1
    ssh_user = "compute"
    workers_ready = True
    if method != "exact":
        sizing = st.radio("Size the chunks by", ["Target compounds per chunk", "Number of chunks"],
                          horizontal=True, key="library-centroid-sizing")
        if sizing == "Target compounds per chunk":
            chunk_size = int(st.number_input("Target compounds per chunk", min_value=1, value=300_000,
                                            step=10_000, key="library-centroid-chunk-size"))
        else:
            default_count = max(1, (row_count + 299_999) // 300_000)
            chunk_count = int(st.number_input("Number of random chunks", min_value=1,
                                             max_value=max(1, row_count), value=default_count,
                                             key="library-centroid-chunk-count"))
            chunk_size = max(1, (row_count + chunk_count - 1) // chunk_count)
        if row_count:
            sizes = balanced_chunk_sizes(row_count, chunk_size, chunk_count)
            st.caption(f"{row_count:,} compounds → {len(sizes):,} random chunks of {min(sizes):,}–{max(sizes):,} compounds. Chunk sizes differ by at most one compound.")
        if method == "sampled":
            candidate_count = int(st.number_input(
                "Maximum representatives in the final comparison", min_value=1, value=300_000,
                step=10_000, key="library-centroid-final-size",
                help="Sample an equal quota from each chunk's representatives. Unused quota from small chunks is redistributed. If fewer representatives exist, use all of them.",
            ))
            st.info("This method samples local representatives before the final clustering. It produces a diverse subset; it does not guarantee a representative within the threshold for every original compound. Final cluster sizes count the sampled candidates.")
        else:
            merge_rounds = int(st.number_input(
                "Maximum additional merge rounds (0 = until stable)", min_value=0, value=3,
                key="library-centroid-merge-rounds",
            ))
            st.info("Every local representative continues into the next round. Cluster sizes retain the original compound populations, but the similarity threshold applies at each merge step. Merging does not guarantee the same threshold between the final representative and every original compound.")
        seed = int(st.number_input("Random seed", min_value=0, max_value=2**53 - 1, value=42,
                                   key="library-centroid-seed"))
        nodes = _worker_choices()
        source_available = any(node["node_id"] == selected_node for node in nodes)
        chosen = st.multiselect(
            "Worker computers", options=[node["node_id"] for node in nodes],
            default=[node["node_id"] for node in nodes],
            format_func=lambda node_id: next(node["name"] for node in nodes if node["node_id"] == node_id)
                + (" · source / coordinator" if node_id == selected_node else ""),
            key="library-centroid-workers-" + selected_node,
            help="One Studio job coordinates all workers. Each computer reserves its CPU slots with its scheduler. Include the source workstation.",
        )
        workers_ready = source_available and selected_node in chosen
        if not workers_ready:
            st.warning("Include the source workstation among the workers; it must be online and accepting jobs.")
        selected_nodes = [node for node in nodes if node["node_id"] in chosen]
        if selected_nodes:
            st.markdown("**CPU slots per computer**")
            columns = st.columns(min(3, len(selected_nodes)))
            for index, node in enumerate(selected_nodes):
                limit = node["cpu_available"]
                count = int(columns[index % len(columns)].number_input(
                    node["name"], min_value=1, max_value=limit, value=limit,
                    help=f"Default: all {limit} CPU slots currently available from this computer's scheduler.",
                    key="library-centroid-cpus-v1-" + node["node_id"],
                ))
                worker_nodes.append({"node_id": node["node_id"], "cpu_cores": count})
                if node["node_id"] == selected_node:
                    source_cpus = count
        ssh_user = st.text_input("SSH account on worker computers", value="compute", key="library-centroid-ssh")
        st.caption("Independent chunks run concurrently across computers and CPU processes. When fewer chunks are available, CPU threads work within each chunk. Completed outputs return to the coordinator after all workers in the round finish. Temporary chunks are removed; the final Parquet table and manifest are retained.")
    else:
        st.caption("This reference method runs on one CPU and uses a disk-backed index. Large diverse libraries can take substantial time.")
    output_directory = st.text_input(
        "New centroid output directory on the selected workstation",
        help="Use an absolute path outside the Studio workspace and choose a new folder.", key="library-centroid-output",
    )
    organizer_output = ""
    if remote_target:
        organizer_output = st.text_input(
            "New destination directory on this organizer",
            help="Generated centroid outputs are copied here after the job finishes. The original fingerprint files stay on the compute workstation.",
            key="library-centroid-organizer-output",
        )

    ready = bool(project_id and output_directory.strip() and selected_manifest.get("path") and not wrong_owner
                 and workers_ready and ssh_user.strip()
                 and (not remote_target or organizer_output.strip()))
    if st.button("Queue centroid creation", type="primary", disabled=not ready, key="library-centroid-submit"):
        parameters = {
            "library_id": library_id,
            "feature_manifest": selected_manifest["path"],
            "output_directory": output_directory.strip(),
            "similarity_threshold": float(threshold),
            "clustering_method": method, "cpu_cores": source_cpus,
            "chunk_size": chunk_size, "chunk_count": chunk_count,
            "final_candidate_count": candidate_count, "max_merge_rounds": merge_rounds,
            "random_seed": seed, "worker_nodes": worker_nodes, "ssh_user": ssh_user.strip(),
        }
        try:
            submit(schema, "library.centroids", parameters, project_id,
                   {"_execution": {"placement": placement, "idempotency_key": secrets.token_urlsafe(24),
                                   **({"organizer_output_directory": organizer_output.strip()} if remote_target else {})}},
                   "Library centroid creation")
        except (OSError, ValueError) as exc:
            st.error(str(exc))

    with st.expander("External library record"):
        st.json(record)
