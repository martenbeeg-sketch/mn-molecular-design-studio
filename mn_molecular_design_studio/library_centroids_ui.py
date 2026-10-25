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
    source_sample_size = 0
    if method != "exact":
        source_sample_size = int(st.number_input(
            "Random input sample for a pilot (0 = use the full library)",
            min_value=0, max_value=max(0, row_count), value=0, step=1_000_000,
            key="library-centroid-source-sample-size",
            help="Select a deterministic random subset from the prepared fingerprints. The centroid manifest records both the sample size and full source-library size.",
        ))
    active_row_count = source_sample_size or row_count
    if source_sample_size:
        st.info(f"Pilot input: {source_sample_size:,} randomly selected compounds from {row_count:,} prepared compounds.")
    chunk_size, chunk_count, candidate_count, merge_rounds, seed = 300_000, 0, 300_000, 0, 42
    target_reduction_factor, plateau_reduction_percent, plateau_consecutive_rounds = 10, 1.0, 3
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
            default_count = max(1, (active_row_count + 299_999) // 300_000)
            chunk_count = int(st.number_input("Number of random chunks", min_value=1,
                                             max_value=max(1, active_row_count), value=default_count,
                                             key="library-centroid-chunk-count"))
            chunk_size = max(1, (active_row_count + chunk_count - 1) // chunk_count)
        if active_row_count:
            sizes = balanced_chunk_sizes(active_row_count, chunk_size, chunk_count)
            st.caption(f"{active_row_count:,} compounds → {len(sizes):,} random chunks of {min(sizes):,}–{max(sizes):,} compounds. Chunk sizes differ by at most one compound.")
        if method == "sampled":
            candidate_count = int(st.number_input(
                "Maximum representatives in the final comparison", min_value=1, value=300_000,
                step=10_000, key="library-centroid-final-size",
                help="Sample an equal quota from each chunk's representatives. Unused quota from small chunks is redistributed. If fewer representatives exist, use all of them.",
            ))
            st.info("This method samples local representatives before the final clustering. It produces a diverse subset; it does not guarantee a representative within the threshold for every original compound. Final cluster sizes count the sampled candidates.")
        else:
            merge_rounds = int(st.number_input(
                "Maximum additional merge rounds (0 = stopping criteria or stable)", min_value=0, value=0,
                key="library-centroid-merge-rounds",
            ))
            stop_columns = st.columns(3)
            target_reduction_factor = int(stop_columns[0].number_input(
                "Stop at reduction factor (0 disables)", min_value=0, max_value=1000, value=10,
                key="library-centroid-stop-factor",
                help="10 stops when the representative count is at most one tenth of the input count.",
            ))
            plateau_reduction_percent = float(stop_columns[1].number_input(
                "Plateau threshold (%) (0 disables)", min_value=0.0, max_value=100.0, value=1.0,
                step=0.1, key="library-centroid-plateau-percent",
                help="A round is considered stalled when it reduces representatives by less than this percentage.",
            ))
            plateau_consecutive_rounds = int(stop_columns[2].number_input(
                "Consecutive plateau rounds", min_value=1, max_value=100, value=3,
                key="library-centroid-plateau-rounds",
                help="Require this many consecutive low-reduction rounds before stopping on the plateau rule.",
            ))
            st.info("Stop when either the reduction factor is reached or the plateau threshold holds for the configured consecutive rounds. Every local representative continues into the next round until then. Cluster sizes retain original compound populations, but the similarity threshold applies at each merge step; merging does not guarantee the same threshold between the final representative and every original compound.")
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
            help="One Studio job coordinates all workers. Enabled computers remain selectable if their peer heartbeat is stale; the job checks each scheduler and CPU allocation live before starting. Include the source workstation.",
        )
        workers_ready = source_available and selected_node in chosen
        if not workers_ready:
            st.warning("Include the source workstation among the workers.")
        selected_nodes = [node for node in nodes if node["node_id"] in chosen]
        if selected_nodes:
            st.markdown("**CPU slots and chunk share per computer**")
            columns = st.columns(min(3, len(selected_nodes)))
            for index, node in enumerate(selected_nodes):
                limit = node["cpu_available"]
                default_cpus = min(limit, 32) if node["node_id"] == selected_node else limit
                count = int(columns[index % len(columns)].number_input(
                    node["name"], min_value=1, max_value=limit, value=default_cpus,
                    help=(f"Default: {default_cpus} slots. The organizer uses its tested 32-slot setting; remote computers default to all {limit} slots currently available from their schedulers."),
                    key="library-centroid-cpus-v3-" + node["node_id"],
                ))
                max_parallel_default = 1 if node["node_id"] == selected_node else count
                max_parallel = int(columns[index % len(columns)].number_input(
                    f"Concurrent chunks · {node['name']}", min_value=1, max_value=count,
                    value=min(count, max_parallel_default),
                    help=("Organizer default is one chunk at a time: the measured large-round throughput was better than running many 2-core chunk tasks concurrently. "
                          "Remote defaults use all selected CPU slots, distributing them across simultaneous chunks."),
                    key="library-centroid-max-parallel-v1-" + node["node_id"],
                ))
                chunk_weight = int(columns[index % len(columns)].number_input(
                    f"Chunk share · {node['name']}", min_value=1, max_value=1000,
                    value=1 if node["node_id"] == selected_node else 4,
                    help="Relative number of chunks assigned to this computer. Defaults reflect the measured first-round aggregate throughput (organizer 1, each remote about 4); adjust after benchmarking.",
                    key="library-centroid-chunk-weight-v2-" + node["node_id"],
                ))
                worker_nodes.append({"node_id": node["node_id"], "cpu_cores": count,
                                     "chunk_weight": chunk_weight,
                                     "max_parallel_chunks": max_parallel})
                if node["node_id"] == selected_node:
                    source_cpus = count
        ssh_user = st.text_input("SSH account on worker computers", value="compute", key="library-centroid-ssh")
        st.caption("Independent chunks run across computers. CPU slots are divided across the configured concurrent chunks; chunk shares set each computer's relative amount of work. Defaults use the measured throughput of these three computers. The selected schedulers and free CPU slots are checked live. Completed outputs return to the coordinator after all workers in the round finish. Temporary chunks are removed; the final Parquet table and manifest are retained.")
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
            "source_sample_size": source_sample_size,
            "final_candidate_count": candidate_count, "max_merge_rounds": merge_rounds,
            "target_reduction_factor": target_reduction_factor,
            "plateau_reduction_percent": plateau_reduction_percent,
            "plateau_consecutive_rounds": plateau_consecutive_rounds,
            "random_seed": seed, "worker_nodes": worker_nodes, "ssh_user": ssh_user.strip(),
        }
        try:
            submit(schema, "library.centroids", parameters, project_id,
                   {"_execution": {"placement": placement, "idempotency_key": secrets.token_urlsafe(24),
                                   **({"organizer_output_directory": organizer_output.strip()} if remote_target else {})}},
                   (f"Library centroid pilot · {source_sample_size:,} input compounds"
                    if source_sample_size else "Library centroid creation"))
        except (OSError, ValueError) as exc:
            st.error(str(exc))

    with st.expander("External library record"):
        st.json(record)
