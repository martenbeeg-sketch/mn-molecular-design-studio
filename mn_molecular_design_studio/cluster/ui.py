"""Shared workstation controls; remote files are opened only on request."""
from __future__ import annotations

import hashlib
import json
import sys
import uuid
from pathlib import Path

import pandas as pd
import streamlit as st

from . import config, transport


_TOOL_LABELS = {
    "visual_pipeline": "Visual pipeline",
    "optimized_prediction": "Structure Prediction",
    "bindcraft": "BindCraft",
    "bindcraft2": "BindCraft2",
    "bindcraft2_generator": "BindCraft2 Generator",
    "boltzgen": "BoltzGen",
    "boltzgen_generator": "BoltzGen Generator",
    "boltzgen_fold": "BoltzGen Fold",
    "proteina_complexa": "Proteina-Complexa",
    "proteina_complexa_generator": "Proteina-Complexa Generator",
    "pxdesign": "PXDesign",
    "pxdesign_generator": "PXDesign",
    "rfdiffusion": "RFdiffusion",
    "rfdiffusion_classic": "RFdiffusion",
    "rfdiffusion_generator": "RFdiffusion Generator",
    "rfdiffusion3": "RFdiffusion3",
    "rfdiffusion3_foundry": "RFdiffusion3",
    "genie3": "Genie3",
    "genie3_generator": "Genie3 Generator",
    "protpardelle_1c": "Protpardelle-1C",
    "protpardelle_generator": "Protpardelle-1C Generator",
    "rfantibody": "RFantibody",
    "esmfold2_binder_design": "ESMFold2 Binder Design",
    "esmfold2_binder": "ESMFold2 Binder Design",
    "esmfold2_cofolding": "ESMFold2",
    "esmfold2_monomer": "ESMFold2",
    "esmfold2_complex_validation": "ESMFold2",
    "esmfold": "ESMFold",
    "af2_initial_guess": "AlphaFold2 Initial Guess",
    "alphafold2": "AlphaFold2",
    "af3_jax_cofolding": "ColabFold 2.0 preview",
    "alphafold3": "AlphaFold3 JAX",
    "af3_torch": "AlphaFold3 PyTorch",
    "openfold3": "OpenFold3",
    "rf3": "RosettaFold3",
    "protenix": "Protenix",
    "protenix_v1": "Protenix CLI",
    "boltz2": "Boltz-2",
    "boltz2_complex": "Boltz-2",
    "boltz2_monomer": "Boltz-2",
    "boltz2_initial_guess": "Boltz-2",
    "ligandmpnn": "ProteinMPNN",
    "foundry_mpnn": "ProteinMPNN (Foundry)",
    "esm_if1": "ESM-IF1",
    "caliby": "Caliby",
    "boltzgen_inverse": "BoltzGen Inverse Folding",
    "scannet": "ScanNet",
    "pesto": "PeSTo",
    "surf2spot": "Surf2Spot",
    "masif_seed": "MaSIF-seed",
    "ipsae": "ipSAE Scoring",
    "rosetta_metrics": "Rosetta Interface Metrics",
}


def display_tool_label(row):
    tool = str(row.get("tool") or row.get("job_type") or "").strip()
    optimized_prefix = "optimized_prediction_"
    if tool.startswith(optimized_prefix):
        engine = tool[len(optimized_prefix):]
        try:
            from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
            profile = ENGINE_PROFILES.get(engine)
            if profile:
                return str(profile.get("label") or engine)
        except (ImportError, AttributeError):
            pass
    if "." in tool:
        try:
            from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
            spec = WORKFLOWS.get(tool)
            if spec:
                tool = spec.tool
        except (ImportError, AttributeError):
            pass
    if tool in _TOOL_LABELS:
        return _TOOL_LABELS[tool]
    workflow_aliases = {
        "design.bindcraft2": "BindCraft2",
        "design.bindcraft2_generator": "BindCraft2 Generator",
        "design.boltzgen": "BoltzGen",
        "design.boltzgen_generator": "BoltzGen Generator",
        "design.proteina_complexa": "Proteina-Complexa",
        "design.proteina_complexa_generate": "Proteina-Complexa Generator",
        "design.pxdesign_generate": "PXDesign",
        "design.rfdiffusion_generate": "RFdiffusion Generator",
        "design.rfdiffusion3": "RFdiffusion3",
        "design.genie3_generate": "Genie3 Generator",
        "design.protpardelle_generate": "Protpardelle-1C Generator",
        "sequence.mpnn": "ProteinMPNN",
        "sequence.ligandmpnn": "ProteinMPNN",
        "sequence.foundry_mpnn": "ProteinMPNN (Foundry)",
        "sequence.esm_if1": "ESM-IF1",
        "sequence.caliby": "Caliby",
        "sequence.boltzgen_inverse": "BoltzGen Inverse Folding",
        "refolding.af3_jax_cofolding": "ColabFold 2.0 preview",
        "refolding.af2_initial_guess": "AlphaFold2 Initial Guess",
        "refolding.openfold3": "OpenFold3",
        "refolding.rf3": "RosettaFold3",
        "refolding.protenix": "Protenix",
        "refolding.boltzgen_fold": "BoltzGen Fold",
    }
    if tool in workflow_aliases:
        return workflow_aliases[tool]
    return tool.replace("_", " ").strip().title() or "Unknown tool"


def _request(path, token, payload=None, method="GET", **kwargs):
    return transport.request("/api/v1/cluster" + path, payload, token=token, method=method, **kwargs)


def _cache(user_id, name):
    safe = hashlib.sha256(str(user_id).encode()).hexdigest()[:24]
    return config.home() / (name + "-" + safe + ".json")


def catalog(user, token):
    path = _cache(user["id"], "catalog-cache")
    try:
        result = _request("/jobs?limit=1000", token)
        total = result["total"]
        while len(result["items"]) < total:
            page = _request(f"/jobs?limit=1000&offset={len(result['items'])}", token)
            if not page["items"]:
                break
            result["items"].extend(page["items"])
        config.atomic_json(path, result)
        return result, False
    except transport.ClusterUnavailable:
        try:
            return json.loads(path.read_text()), True
        except (FileNotFoundError, json.JSONDecodeError):
            return {"items": [], "total": 0}, True


def execution_selector(user, token):
    if not config.enabled():
        return
    settings = config.load()
    try:
        nodes = _request("/nodes", token)["nodes"]
        config.atomic_json(_cache(user["id"], "nodes-cache"), {"nodes": nodes})
    except (transport.ClusterUnavailable, ValueError):
        try:
            nodes = json.loads(_cache(user["id"], "nodes-cache").read_text())["nodes"]
        except (FileNotFoundError, json.JSONDecodeError):
            nodes = []
    names = {n["id"]: f"{config.display_name(n['name'])} · {n['availability']}" for n in nodes if n["id"] != config.identity()}
    options = ["auto", "local", *names]
    if st.session_state.get("cluster-placement") not in options:
        st.session_state["cluster-placement"] = settings["default_placement"]
    st.sidebar.selectbox("Run new jobs on", options, format_func=lambda value: {"auto": "Automatic · available workstation", "local": "This workstation"}.get(value, names.get(value, value)), key="cluster-placement")
    st.sidebar.caption(f"This workstation: {config.display_name(settings['name'])}")
    if nodes:
        with st.sidebar.expander("Workstation resources"):
            for node in nodes:
                resource = node["snapshot"].get("resources", {})
                st.write(f"**{config.display_name(node['name'])}** · {node['availability']}")
                st.caption(f"{len(resource.get('free_gpus', []))} free GPUs · {resource.get('cpu_free', 0)} free CPU threads · {resource.get('ram_available_gb', 0):.0f} GB RAM available")


def settings(user, token):
    st.subheader("Workstations")
    current = config.load()
    st.caption("Share the job catalog and send new work to available workstations. Full results stay on their execution machine.")
    if user.get("role") != "admin":
        st.write(f"{current['name']} · {current['mode']}")
        return
    modes = ["disabled", "peer", "coordinator", "worker"]
    mode = st.selectbox("Workstation connection", modes, index=modes.index(current["mode"]),
                        format_func=lambda v: {"disabled": "Standalone", "peer": "Independent workstations", "coordinator": "Coordinator (legacy)", "worker": "Connected workstation (legacy)"}[v], key="cluster-settings-mode")
    name_key = "cluster-settings-name"
    if st.session_state.get(name_key) == current["name"]:
        st.session_state[name_key] = config.display_name(current["name"])
    name = st.text_input("Workstation name", value=config.display_name(current["name"]), key=name_key)
    address = st.text_input("Local API address" if mode == "peer" else "Coordinator API address", value=current["coordinator_url"] or "http://127.0.0.1:8766", key="cluster-settings-url")
    accept = st.checkbox("Accept new jobs on this workstation", value=current["accept_jobs"], key="cluster-settings-accept")
    st.caption(f"Permanent workstation ID: {config.identity()}")
    if mode == "peer" and current["mode"] != "peer":
        st.info("Connect trusted peers with the cluster enable-peers command. Each installation keeps its accounts, scheduler and results available independently.")
    elif mode == "worker" and current["mode"] != "worker":
        st.info("Enroll this workstation on the coordinator, then import its private enrollment file with the cluster join command.")
        st.code("mn-molecular-design-studio cluster join workstation-enrollment.json --coordinator-url " + address, language="bash")
    elif st.button("Save workstation settings", key="cluster-settings-save"):
        try:
            changes = {"mode": mode, "name": name, "coordinator_url": address, "accept_jobs": accept}
            if mode == "coordinator" and not current["node_token"]:
                from .store import Store
                enrollment = Store().enroll(name, config.identity())
                changes.update(node_token=enrollment["node_token"], default_placement="auto")
            config.save(changes)
            st.success("Workstation settings saved. Start or restart the API and agent services to apply connection changes.")
        except (ValueError, OSError) as exc:
            st.error(str(exc))
    if mode != "disabled":
        st.code("mn-molecular-design-studio cluster install-services --enable", language="bash")
        st.caption("The background agent publishes resources and accepts work while the browser interface is closed.")
    if current["mode"] == "disabled":
        return
    try:
        nodes = _request("/nodes", token)["nodes"]
        st.dataframe(pd.DataFrame([{"Workstation": config.display_name(n["name"]), "Connection": n["availability"], "Last seen": n["last_seen_at"],
                                   "Free GPUs": len(n["snapshot"].get("resources", {}).get("free_gpus", [])),
                                   "Free CPU threads": n["snapshot"].get("resources", {}).get("cpu_free"),
                                   "Accepts jobs": n["snapshot"].get("accept_jobs", False)} for n in nodes]), hide_index=True, use_container_width=True)
    except (transport.ClusterUnavailable, ValueError) as exc:
        st.info(str(exc))
    if current["mode"] == "peer":
        st.caption("Accounts and projects are stored locally and synchronized with trusted workstations. Each machine remains usable when the others are offline.")
        from .replication import Replica
        with Replica().db() as db:
            conflicts = db.execute("SELECT key,error FROM conflicts ORDER BY key").fetchall()
        if conflicts:
            st.warning(f"{len(conflicts)} metadata changes need administrator attention. Existing account and job identities have been preserved.")
            st.dataframe([dict(r) for r in conflicts], hide_index=True)
    if current["mode"] == "coordinator":
        with st.expander("Enroll another workstation"):
            new_name = st.text_input("New workstation name", placeholder="pc004032", key="cluster-enroll-name")
            if st.button("Create enrollment file", key="cluster-enroll-create"):
                try:
                    enrollment = _request("/nodes", token, {"name": new_name}, method="POST")
                    st.session_state["cluster-enrollment"] = enrollment
                except (ValueError, transport.ClusterUnavailable) as exc:
                    st.error(str(exc))
            enrollment = st.session_state.get("cluster-enrollment")
            if enrollment:
                st.download_button("Download private enrollment file", json.dumps(enrollment, indent=2),
                                   file_name=f"studio-{enrollment['node_id']}-enrollment.json", mime="application/json")
                st.caption("This file grants access to one workstation agent. Copy it privately to that workstation.")


def render_runs(user, token, local_search):
    st.title("My results")
    page, stale = catalog(user, token)
    if stale:
        st.info("Coordinator unavailable. Showing the last synchronized catalog and current local jobs. Remote status is last known.")
    rows = list(page["items"])
    local_page = local_search(limit=1000)
    local = list(local_page["items"])
    while len(local) < local_page["total"]:
        more = local_search(limit=1000, offset=len(local))["items"]
        if not more:
            break
        local.extend(more)
    known = {(r.get("node_id"), r.get("task_group"), r.get("run_id")): i for i, r in enumerate(rows)}
    for row in local:
        key = (config.identity(), row["task_group"], row["run_id"])
        item = dict(row, node_id=config.identity(), node_name=config.load()["name"], availability="online")
        if key in known:
            rows[known[key]] = rows[known[key]] | item
        else:
            rows.append(item)
    search_col, machine_col, refresh_col = st.columns([0.5, 0.35, 0.15])
    query = search_col.text_input("Search jobs", placeholder="Name, engine, or run ID", key="cluster-results-query").casefold()
    names = {r.get("node_id"): config.display_name(r.get("node_name") or "Not assigned") for r in rows if r.get("node_id")}
    machine = machine_col.selectbox("Workstation", ["", *names], format_func=lambda v: names.get(v, "All workstations"), key="cluster-results-machine")
    if refresh_col.button("Refresh", key="cluster-results-refresh"):
        st.rerun()
    projects = sorted({r.get("project_id") for r in rows if r.get("project_id")})
    from mn_molecular_design_studio.workspace import load_organization
    project_names = {p["id"]: p["name"] for p in load_organization().get("projects", [])}
    project = st.selectbox("Project", ["", *projects], format_func=lambda v: project_names.get(v, "All projects"), key="cluster-results-project")
    chosen = [r for r in rows if (not machine or r.get("node_id") == machine) and (not project or r.get("project_id") == project)
              and (not query or query in " ".join(f"{r.get(k, '')} {config.display_text(r.get(k, ''))}" for k in ("campaign_name", "pipeline_name", "run_id", "tool", "job_code")).casefold())]
    chosen.sort(key=lambda r: r.get("created_at", ""), reverse=True)
    st.caption(f"{len(chosen)} jobs · Job status and file availability are tracked separately.")
    pages = max(1, (len(chosen) + 49) // 50)
    page_key = "cluster-results-page"
    if st.session_state.get(page_key, 1) > pages:
        st.session_state[page_key] = 1
    selected_page = int(st.number_input("Page", min_value=1, max_value=pages, value=1, step=1, key=page_key))
    header = st.columns([3, 1, 1, 1.3, 0.7])
    for column, label in zip(header, ("Job name · tool · pipeline", "Status", "Workstation", "File availability", "Open")):
        column.caption(label)
    for row in chosen[(selected_page - 1) * 50:selected_page * 50]:
        with st.container(border=True):
            title, state, machine, access, action = st.columns([3, 1, 1, 1.3, 0.7])
            title.write(config.display_text(row.get("campaign_name") or row.get("tool") or "Job"))
            pipeline_name = str(row.get("pipeline_name") or "").strip() or "—"
            title.caption(f"Tool: {display_tool_label(row)} · Pipeline: {config.display_text(pipeline_name)}")
            title.caption(project_names.get(row.get("project_id"), "Project") + " · " + str(row.get("job_code", "")))
            state.write(str(row.get("status", "unknown")).replace("_", " ").title())
            machine.write(config.display_name(row.get("node_name") or "Not assigned"))
            availability = "last known" if stale and row.get("node_id") != config.identity() else row.get("availability", "pending")
            access.write(availability.title())
            if row.get("last_seen_at"):
                access.caption(str(row["last_seen_at"]))
            ref = "cluster:" + row["cluster_job_id"] if row.get("cluster_job_id") else row["run_id"]
            if action.button("View", key="cluster-view-" + ref):
                if row.get("node_id") == config.identity() and row.get("run_id"):
                    ref = row["run_id"]
                st.session_state["selected-job-reference"] = ref
                st.session_state["page"] = "Job details"
                st.rerun()


def render_job(reference, user, token, submit):
    jid = reference.removeprefix("cluster:")
    try:
        row = _request("/jobs/" + jid, token)
    except transport.ClusterUnavailable:
        cached, _ = catalog(user, token)
        row = next((j for j in cached["items"] if j.get("cluster_job_id") == jid), None)
        if not row:
            st.info("The coordinator and this job's cached record are unavailable.")
            return
        row = dict(row, availability="last known")
    st.title(config.display_text(row.get("campaign_name") or row.get("tool") or "Job"))
    st.caption(f"Tool: {display_tool_label(row)} · Pipeline: {config.display_text(str(row.get('pipeline_name') or '').strip() or '—')}")
    st.caption(f"Files stored on {config.display_name(row.get('node_name') or 'Not assigned')} · {row.get('location') or 'Waiting for placement'}")
    a, b, c = st.columns(3)
    a.metric("Job status", row.get("status", "unknown").replace("_", " ").title())
    b.metric("Workstation", config.display_name(row.get("node_name") or "Not assigned"))
    c.metric("Data availability", row.get("availability", "unknown").title())
    if row.get("availability") != "online":
        st.info("The saved job record remains available. Full files and live analysis will be available when the owning workstation reconnects.")
    if row.get("wait_reason"):
        st.info(row["wait_reason"])
    if row.get("error"):
        st.error(row["error"])
    if row.get("status") not in {"completed", "succeeded", "failed", "cancelled"}:
        if st.button("Request cancellation", key="cluster-cancel-" + jid):
            response = _request("/jobs/" + jid + "/actions", token, {"action": "cancel"}, method="POST")
            st.success("Cancellation requested. The workstation will confirm when execution has stopped.")
    st.json(row.get("metrics", {}), expanded=False)
    if row.get("local_reference"):
        if st.button("Load full results from workstation", key="cluster-inspect-" + jid):
            operation = _request("/jobs/" + jid + "/actions", token, {"action": "inspect"}, method="POST")
            st.session_state["cluster-details-op-" + jid] = operation["command_id"]
        @st.fragment(run_every=2)
        def inspect_operation():
            command_id = st.session_state.get("cluster-details-op-" + jid)
            if not command_id:
                return
            try:
                op = _request("/commands/" + command_id, token)
            except (ValueError, transport.ClusterUnavailable) as exc:
                st.info(str(exc)); return
            if op["status"] == "failed":
                st.error(op["error"]); return
            if op["status"] != "completed":
                st.info("Waiting for the owning workstation to return the results."); return
            data = op["response"]
            candidates = data.get("results", {}).get("candidates", [])
            if candidates:
                st.dataframe(pd.json_normalize([{ "candidate_id": r.get("candidate_id"), **r.get("metrics", {})} for r in candidates]), hide_index=True)
            with st.expander("Run settings and inputs"):
                st.json(data.get("input", {}))
            with st.expander("Run summary"):
                st.json(data.get("result", {}))
        inspect_operation()
        artifacts = row.get("artifacts", [])
        if artifacts:
            selected = st.selectbox("Saved file", range(len(artifacts)), format_func=lambda i: artifacts[i].get("relative_path", artifacts[i].get("name", "File")), key="cluster-file-" + jid)
            if st.button("Prepare file for viewing or download", key="cluster-file-load-" + jid):
                path = artifacts[selected].get("relative_path")
                operation = _request("/jobs/" + jid + "/actions", token, {"action": "artifact", "parameters": {"path": path}}, method="POST")
                st.session_state["cluster-file-op-" + jid] = operation["command_id"]
            @st.fragment(run_every=2)
            def file_operation():
                command_id = st.session_state.get("cluster-file-op-" + jid)
                if not command_id:
                    return
                try:
                    operation = _request("/commands/" + command_id, token)
                    if operation["status"] == "failed":
                        st.error(operation["error"]); return
                    if operation["status"] != "completed":
                        st.info("Waiting for the owning workstation to serve this file."); return
                    key = "cluster-file-bytes-" + command_id
                    if key not in st.session_state:
                        st.session_state[key] = _request("/commands/" + command_id + "/file", token, raw=True, timeout=60)
                    content = st.session_state[key]
                    name = operation["response"].get("name", "result")
                    st.download_button("Download " + name, content, file_name=name, key="cluster-download-" + command_id)
                    if Path(name).suffix.lower() in {".pdb", ".cif", ".mmcif"} and len(content) < 50 * 1024 * 1024:
                        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
                        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization
                        molstar_custom_component([StructureVisualization(pdb=content.decode("utf-8", errors="replace"), color="plddt")], key="cluster-viewer-" + command_id, height=520, show_controls=False, show_structure_tools=False)
                except (ValueError, transport.ClusterUnavailable) as exc:
                    st.info(str(exc))
            file_operation()
        with st.expander("Continue on the workstation that owns these results"):
            workflow = st.selectbox("Next step", ["analysis.ipsae", "analysis.rosetta", "sequence.mpnn", "refolding.boltz2_complex"],
                                    format_func=lambda v: {"analysis.ipsae": "ipSAE interface confidence", "analysis.rosetta": "Rosetta interface analysis", "sequence.mpnn": "ProteinMPNN sequence redesign", "refolding.boltz2_complex": "Boltz-2 complex prediction"}[v], key="cluster-analysis-" + jid)
            if workflow == "analysis.ipsae":
                parameters = {"pae_cutoff": st.number_input("PAE cutoff", 1.0, 40.0, 10.0, key="cluster-pae-" + jid),
                              "dist_cutoff": st.number_input("Distance cutoff", 1.0, 40.0, 10.0, key="cluster-dist-" + jid)}
            elif workflow == "analysis.rosetta":
                parameters = {"nprocs": st.number_input("CPU processes", 1, 64, 4, key="cluster-nprocs-" + jid)}
            elif workflow == "sequence.mpnn":
                parameters = {"model_type": st.selectbox("MPNN model", ["protein_mpnn", "ligand_mpnn", "soluble_mpnn", "hyper_mpnn", "ab_mpnn"], key="cluster-mpnn-" + jid),
                              "num_seq_per_target": st.number_input("Sequences per design", 1, 1000, 2, key="cluster-sequences-" + jid),
                              "fix_interface_residues": st.checkbox("Keep interface residues when a valid source sequence is available", key="cluster-interface-" + jid)}
            else:
                parameters = {"require_monomer_success": False, "recycling_steps": st.number_input("Recycles", 1, 20, 3, key="cluster-recycles-" + jid),
                              "sampling_steps": st.number_input("Sampling steps", 1, 1000, 200, key="cluster-steps-" + jid),
                              "diffusion_samples": st.number_input("Samples", 1, 20, 1, key="cluster-samples-" + jid)}
            if st.button("Queue on " + config.display_name(row.get("node_name") or "owning workstation"), type="primary", key="cluster-next-" + jid):
                result = submit("/api/v1/workflows", {"schema_version": 1, "workflow": workflow, "parameters": parameters,
                                                     "project_id": row["project_id"], "source_job": "cluster:" + jid,
                                                     "execution": {"placement": "auto", "idempotency_key": str(uuid.uuid4())}})
                st.success("Queued with data locality preserved · " + result["job"].get("job_code", "submitted"))
