"""Shared workstation controls; remote files are opened only on request."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import sys
import uuid
import csv
import io
from pathlib import Path, PurePosixPath
from urllib.parse import quote

import pandas as pd
import streamlit as st

from . import config, progress as job_progress, transport


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
    "rfdiffusion_beta_pairing": "RFdiffusion β-pairing",
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
        "design.rfdiffusion_beta_pairing": "RFdiffusion β-pairing",
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


def _is_bindcraft2_job(row):
    identity = " ".join(str(row.get(key) or "") for key in ("tool", "job_type", "workflow")).casefold()
    return "bindcraft2" in identity.replace(" ", "") or "bindcraft 2" in identity


def _row_progress(row):
    progress = row.get("progress")
    if isinstance(progress, dict) and progress.get("engine") == "bindcraft2":
        return progress
    if row.get("node_id") != config.identity() or not _is_bindcraft2_job(row):
        return None
    run_dir = row.get("run_dir")
    if not isinstance(run_dir, str) or not run_dir:
        return None
    try:
        return job_progress.bindcraft2_progress(Path(run_dir))
    except OSError:
        return None


def _progress_summary(progress):
    if not isinstance(progress, dict):
        return ""
    completed = int(progress.get("trajectories_completed") or 0)
    active = int(progress.get("trajectories_in_progress") or 0)
    accepted = int(progress.get("designs_accepted") or 0)
    trajectory_limit = progress.get("trajectory_limit")
    design_target = progress.get("design_target")
    trajectory_label = f"{completed:,} / {int(trajectory_limit):,}" if trajectory_limit else f"{completed:,} completed"
    design_label = f"{accepted:,} / {int(design_target):,}" if design_target else f"{accepted:,} accepted"
    return f"Progress: {trajectory_label} trajectories · {active:,} running · {design_label} accepted"


def _render_bindcraft2_progress(row, job_id, token):
    terminal = {"completed", "succeeded", "failed", "cancelled"}
    run_every = 15 if row.get("status") not in terminal else None

    @st.fragment(run_every=run_every)
    def live_progress():
        current = row
        try:
            current = _request("/jobs/" + job_id, token)
        except (transport.ClusterUnavailable, ValueError):
            pass
        progress = _row_progress(current)
        st.subheader("BindCraft progress")
        if not progress:
            st.caption("Live campaign counters will appear when BindCraft starts writing progress.")
            return

        completed = int(progress.get("trajectories_completed") or 0)
        active = int(progress.get("trajectories_in_progress") or 0)
        accepted = int(progress.get("designs_accepted") or 0)
        trajectory_limit = progress.get("trajectory_limit")
        design_target = progress.get("design_target")
        attempts_col, accepted_col = st.columns(2)
        attempts_value = f"{completed:,}" + (f" / {int(trajectory_limit):,}" if trajectory_limit else "")
        accepted_value = f"{accepted:,}" + (f" / {int(design_target):,}" if design_target else "")
        attempts_col.metric("Completed trajectories", attempts_value, delta=f"{active:,} running" if active else None)
        accepted_col.metric("Accepted designs", accepted_value)
        if trajectory_limit:
            attempts_col.progress(min(1.0, completed / max(1, int(trajectory_limit))), text="Attempt cap")
        if design_target:
            accepted_col.progress(min(1.0, accepted / max(1, int(design_target))), text="Accepted-design goal")
        if progress.get("updated_at"):
            st.caption("Counters update while the campaign runs; the trajectory count is completed attempts, with active workers shown separately.")

    live_progress()


def _candidate_structure_path(candidate, data):
    """Return a run-relative candidate structure path safe for the artifact API."""
    if not isinstance(candidate, dict) or not isinstance(data, dict):
        return None
    run_dir = (data.get("job") or {}).get("run_dir")
    if not isinstance(run_dir, str) or not run_dir:
        return None
    root = PurePosixPath(run_dir)
    for key in ("complex_pdb", "binder_pdb"):
        value = candidate.get(key)
        if not isinstance(value, str) or not value.strip():
            continue
        path = PurePosixPath(value.strip())
        if path.is_absolute():
            try:
                path = path.relative_to(root)
            except ValueError:
                continue
        if not path.parts or any(part in {".", ".."} for part in path.parts):
            continue
        return str(path), key
    return None


_BINDCRAFT2_RESULT_COLUMNS = (
    "i_pDAE",
    "i_pTM",
    "pLDDT",
    "pTM",
    "i_pAE",
    "Unbound_Binder_pLDDT",
    "Target_pLDDT",
    "Interface_Residues",
    "Binder_RMSD",
    "Target_Crop_Length",
    "Backbone_Clashes",
    "All_Atom_Clashes",
    "Interface_BuriedArea",
    "Interface_BuriedArea_Fraction",
    "Surface_Hydrophobicity",
    "Interface_Hydrophobicity",
    "Binder_Helix_Fraction",
    "Binder_BetaSheet_Fraction",
    "Binder_Loop_Fraction",
    "SS_pLDDT",
    "Interface_Helix_Fraction",
    "Interface_BetaSheet_Fraction",
    "Interface_Loop_Fraction",
    "Binder_Disulfides",
    "Binder_Mass_kDa",
    "Binder_pI",
    "Binder_Net_Charge",
    "Binder_Cysteines",
    "Binder_Extinction",
    "Binder_Free_Cysteines",
    "Interface_Binder_Residues",
    "Interface_Target_Residues",
)


def _candidate_result_rows(candidates, tool=""):
    """Build a cross-tool result table, retaining BindCraft2's native metrics."""
    is_bindcraft2 = "bindcraft2" in str(tool).casefold() or any(
        isinstance(candidate, dict)
        and "bindcraft2" in str(candidate.get("source_tool") or candidate.get("tool") or "").casefold()
        for candidate in candidates
    )
    metric_columns = (
        ("pred_lddt", "pLDDT"),
        ("interaction_pae", "Interaction PAE"),
        ("pae", "PAE"),
        ("hotspot_atom_contact_fraction", "Hotspot atom contacts"),
        ("hotspot_ca_contact_fraction", "Hotspot CA contacts"),
        ("rf2_repeat", "RF2 repeat"),
        ("cyclic", "Cyclic"),
    )
    rows = []
    for rank, candidate in enumerate(candidates, 1):
        if not isinstance(candidate, dict):
            continue
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        raw_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        native_row = raw_metadata.get("bindcraft2_row") if isinstance(raw_metadata.get("bindcraft2_row"), dict) else {}
        row = {
            "Rank": metrics.get("bindcraft2_rank", native_row.get("rank", rank)) if is_bindcraft2 else rank,
            "Design": (
                native_row.get("design")
                or raw_metadata.get("bindcraft2_design")
                or candidate.get("candidate_id")
                or f"Design {rank}"
            ),
            "Sequence": candidate.get("binder_sequence") or candidate.get("sequence") or "—",
            "Length": candidate.get("binder_length") or metrics.get("binder_length"),
        }
        if is_bindcraft2:
            for key in _BINDCRAFT2_RESULT_COLUMNS:
                value = metrics.get(key, native_row.get(key))
                if value not in (None, ""):
                    row[key] = value
            if "native_pass_filters" in metrics:
                row["Passed filters"] = metrics["native_pass_filters"]
        else:
            for key, label in metric_columns:
                if key in metrics:
                    row[label] = metrics[key]
        rows.append(row)
    return rows


def _bindcraft2_run_parameter_rows(data):
    """Select human-readable BindCraft2 inputs and filter thresholds from inspect data."""
    input_payload = data.get("input") if isinstance(data, dict) else None
    if not isinstance(input_payload, dict):
        return [], []
    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    inputs = input_payload.get("inputs") if isinstance(input_payload.get("inputs"), dict) else {}

    def present(value):
        return value not in (None, "", [], {})

    rows = []
    values = (
        ("Campaign", params.get("campaign_name") or input_payload.get("campaign_name")),
        ("Target", params.get("target_name") or inputs.get("target_name")),
        ("Target type", params.get("target_type") or inputs.get("target_type")),
        ("Target chains", params.get("target_chains") or inputs.get("target_chains")),
        ("Design format", params.get("modality")),
        ("Binder length range", params.get("binder_length")),
        ("Designs requested", params.get("number_of_final_designs")),
        ("Maximum trajectories", params.get("max_trajectories")),
        ("Target crop count", params.get("target_crop_count")),
        ("Target crop lengths", params.get("target_crop_lengths")),
        ("Hotspots", params.get("hotspots")),
        ("Coldspots", params.get("coldspots")),
        ("Random seed", params.get("campaign_seed")),
        ("GPU device", params.get("gpu_device")),
    )
    for label, value in values:
        if present(value):
            if isinstance(value, (list, tuple)):
                value = ", ".join(map(str, value))
            rows.append({"Parameter": label, "Value": value})

    quality = params.get("quality_filters") if isinstance(params.get("quality_filters"), dict) else {}
    filters = quality.get("filters") if isinstance(quality.get("filters"), dict) else {}
    filter_rows = []
    for name, value in quality.items():
        if name == "filters" or value is None:
            continue
        try:
            value = f"{float(value):g}"
        except (TypeError, ValueError):
            value = str(value)
        key = str(name)
        direction = "≤" if key.startswith("max_") else "≥" if key.startswith("min_") else ""
        detail = key.removeprefix("min_").removeprefix("max_")
        parts = detail.split("_")
        if parts[0] in {"iptm", "plddt"} and len(parts) > 1:
            metric = "i_pTM" if parts[0] == "iptm" else "pLDDT"
            label = f"{parts[-1].title()} {metric} stage minimum"
        else:
            label = detail.replace("_", " ").title()
        filter_rows.append({"Filter": label, "Requirement": f"{direction} {value}".strip()})
    for name, config in filters.items():
        if not isinstance(config, dict) or config.get("threshold") is None:
            continue
        threshold = config["threshold"]
        try:
            threshold = float(threshold)
            threshold = f"{threshold:g}"
        except (TypeError, ValueError):
            threshold = str(threshold)
        direction = "≥" if config.get("higher") else "≤"
        filter_rows.append({"Filter": str(name).replace("_", " "), "Requirement": f"{direction} {threshold}"})
    filter_rows.sort(key=lambda item: item["Filter"].casefold())
    return rows, filter_rows


def _safe_run_artifact_path(value):
    if not isinstance(value, str) or not value.strip():
        return None
    path = PurePosixPath(value.strip())
    if path.is_absolute() or not path.parts or any(part in {".", ".."} for part in path.parts):
        return None
    return str(path)


def _bindcraft2_refolded_table_path(data):
    result = data.get("result") if isinstance(data, dict) else None
    outputs = result.get("outputs") if isinstance(result, dict) else None
    artifacts = outputs.get("artifacts") if isinstance(outputs, dict) else None
    for artifact in artifacts if isinstance(artifacts, list) else []:
        if isinstance(artifact, dict) and artifact.get("type") == "bindcraft2_refolded_table":
            path = _safe_run_artifact_path(artifact.get("path"))
            if path:
                return path
    return "artifacts/raw/bindcraft2/output/2_Refolded/!_Refolded.csv"


def _bindcraft2_monomer_path(refolded_csv, candidate):
    raw = candidate.get("raw_metadata") if isinstance(candidate, dict) else None
    raw = raw if isinstance(raw, dict) else {}
    direct_path = _safe_run_artifact_path(raw.get("bindcraft2_binder_monomer"))
    if direct_path:
        return direct_path
    if not refolded_csv:
        return None
    native = raw.get("bindcraft2_row") if isinstance(raw.get("bindcraft2_row"), dict) else {}
    wanted_sequence = "".join(str(candidate.get("binder_sequence") or native.get("Binder_Sequence") or "").split()).upper()
    matches = []
    if wanted_sequence:
        csv_text = str(refolded_csv or "").lstrip("\ufeff")
        for row in csv.DictReader(io.StringIO(csv_text)):
            sequence = "".join(str(row.get("Binder_Sequence") or "").split()).upper()
            if sequence != wanted_sequence:
                continue
            failed = str(row.get("failed_filters") or "").strip()
            outcome = str(row.get("outcome") or "").strip().casefold()
            matches.append((not failed and outcome in {"", "passed"}, row))
    if not matches:
        design = str(raw.get("bindcraft2_design") or native.get("design") or "").strip()
    else:
        _, matched = next((item for item in matches if item[0]), matches[0])
        design = str(matched.get("design") or "").strip()
    design = PurePosixPath(design).name
    if not design or design in {".", ".."}:
        return None
    if design.endswith("_monomer.cif"):
        return _safe_run_artifact_path(
            f"artifacts/raw/bindcraft2/output/2_Refolded/BinderMonomer/{design}"
        )
    if design.endswith((".cif", ".mmcif", ".pdb")):
        design = PurePosixPath(design).stem
    return _safe_run_artifact_path(
        f"artifacts/raw/bindcraft2/output/2_Refolded/BinderMonomer/{design}_monomer.cif"
    )


def _fetch_remote_artifact_text(jid, path, token, state_key, label="saved structure"):
    relative = _safe_run_artifact_path(path)
    if not relative:
        st.warning("The saved structure path is invalid.")
        return None
    command_id = st.session_state.get(state_key)
    if not command_id:
        try:
            operation = _request(
                "/jobs/" + jid + "/actions",
                token,
                {"action": "artifact", "parameters": {"path": relative}},
                method="POST",
            )
            command_id = operation["command_id"]
            st.session_state[state_key] = command_id
        except (ValueError, transport.ClusterUnavailable) as exc:
            st.warning(f"Could not request the saved structure: {exc}")
            return None
    try:
        status = _request("/commands/" + command_id, token)
        if status["status"] == "failed":
            error = status.get("error") or f"The workstation could not serve this {label}."
            st.warning(f"{error} Requested path: {relative}")
            return None
        if status["status"] != "completed":
            st.info(f"Loading the {label} from the owning workstation…")
            return None
        content_key = "cluster-artifact-bytes-" + command_id
        if content_key not in st.session_state:
            st.session_state[content_key] = _request(
                "/commands/" + command_id + "/file", token, raw=True, timeout=60
            )
        return st.session_state[content_key].decode("utf-8", errors="replace")
    except (ValueError, transport.ClusterUnavailable) as exc:
        st.info(str(exc))
        return None


def _align_bindcraft2_monomer_to_complex(complex_text, monomer_text, candidate):
    """Fit the binder-alone CIF onto its binder chain in the predicted complex."""
    from Bio.PDB import MMCIFIO, MMCIFParser, Superimposer
    from mn_molecular_design_studio.design_validation import protein_chains, residue_pairs, sequence
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import _mmcif_text_with_default_occupancy

    def parse(text, name):
        prepared = _mmcif_text_with_default_occupancy(str(text))
        return MMCIFParser(QUIET=True).get_structure(name, io.StringIO(prepared))

    complex_structure = parse(complex_text, "bindcraft2_complex")
    monomer_structure = parse(monomer_text, "bindcraft2_binder_monomer")
    complex_chains = protein_chains(complex_structure)
    monomer_chains = protein_chains(monomer_structure)
    expected = str(candidate.get("binder_sequence") or "").replace(" ", "").upper()
    if not expected:
        raise ValueError("The candidate sequence is missing, so its binder chain cannot be matched.")

    def best_chain(chains, preferred=()):
        ordered = [chain for chain in preferred if chain in chains]
        ordered.extend(chain for chain in chains if chain not in ordered)
        scored = []
        for chain in ordered:
            observed = sequence(chains[chain])
            pairs = residue_pairs(expected, observed)
            identity = sum(expected[i] == observed[j] for i, j in pairs) / max(1, len(pairs))
            scored.append((identity, chain, pairs))
        return max(scored, default=(0.0, None, []))

    binder_chains = candidate.get("binder_chains") or []
    complex_identity, complex_chain, _ = best_chain(complex_chains, binder_chains)
    monomer_identity, monomer_chain, _ = best_chain(monomer_chains)
    if complex_chain is None or monomer_chain is None or min(complex_identity, monomer_identity) < 0.5:
        raise ValueError("Could not identify the binder chain in both structures.")
    complex_residues, monomer_residues = complex_chains[complex_chain], monomer_chains[monomer_chain]
    complex_sequence, monomer_sequence = sequence(complex_residues), sequence(monomer_residues)
    pairs = residue_pairs(complex_sequence, monomer_sequence)
    pairs = [(i, j) for i, j in pairs if "CA" in complex_residues[i] and "CA" in monomer_residues[j]]
    if len(pairs) < 3:
        raise ValueError("Fewer than three matching binder C-alpha atoms are available for alignment.")
    fixed = [complex_residues[i]["CA"] for i, _ in pairs]
    moving = [monomer_residues[j]["CA"] for _, j in pairs]
    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    for atom in monomer_structure.get_atoms():
        atom.transform(*superimposer.rotran)
    output = io.StringIO()
    writer = MMCIFIO()
    writer.set_structure(monomer_structure)
    writer.save(output)
    return output.getvalue(), float(superimposer.rms), len(pairs)


def _candidate_structure_color(candidate, structure_text=""):
    """Use pLDDT colors only when per-atom confidence values are in the PDB."""
    metrics = candidate.get("metrics") if isinstance(candidate, dict) else None
    metrics = metrics if isinstance(metrics, dict) else {}
    has_confidence_metric = any(
        isinstance(metrics.get(key), (int, float))
        for key in ("pred_lddt", "plddt", "mean_plddt")
    )
    if has_confidence_metric:
        factors = []
        for line in str(structure_text).splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    factors.append(float(line[60:66]))
                except ValueError:
                    continue
        if len(factors) > 1 and max(factors) - min(factors) > 0.01:
            return "plddt"
    return "chain-id"


def _cache(user_id, name):
    safe = hashlib.sha256(str(user_id).encode()).hexdigest()[:24]
    return config.home() / (name + "-" + safe + ".json")


def catalog(user, token, *, include_hidden=False):
    path = _cache(user["id"], "catalog-cache-hidden" if include_hidden else "catalog-cache")
    try:
        query = "?limit=1000" + ("&include_hidden=true" if include_hidden else "")
        result = _request("/jobs" + query, token)
        total = result["total"]
        while len(result["items"]) < total:
            page = _request(f"/jobs?limit=1000&offset={len(result['items'])}" + ("&include_hidden=true" if include_hidden else ""), token)
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


def render_runs(user, token, local_search, *, elapsed_label=None, elapsed_since=None):
    run_every = 2 if st.session_state.get("cluster-results-purge-operations") else 15
    @st.fragment(run_every=run_every)
    def live_results():
        _render_runs_content(user, token, local_search, elapsed_label=elapsed_label, elapsed_since=elapsed_since)

    live_results()


def _render_runs_content(user, token, local_search, *, elapsed_label=None, elapsed_since=None):
    if st.session_state.pop("cluster-results-reset-selection", False):
        for key in list(st.session_state):
            if str(key).startswith("cluster-results-select-"):
                st.session_state.pop(key, None)
        st.session_state["cluster-results-selected"] = {}
    st.title("My results")
    if message := st.session_state.pop("cluster-results-flash", None):
        st.success(message)
    if failures := st.session_state.pop("cluster-results-purge-errors", None):
        for failure in failures:
            st.error(failure)
    purge_operations = list(st.session_state.get("cluster-results-purge-operations") or [])
    if purge_operations:
        remaining_operations = []
        completed_operations = 0
        for operation in purge_operations:
            try:
                result = _request("/commands/" + operation["command_id"], token)
            except (ValueError, transport.ClusterUnavailable) as exc:
                remaining_operations.append(operation)
                st.info(f"Waiting for purge on {operation.get('workstation', 'the owning workstation')}: {exc}")
                continue
            if result.get("status") == "completed":
                completed_operations += 1
            elif result.get("status") == "failed":
                st.error(
                    f"Could not purge {operation.get('name') or 'a selected run'} on "
                    f"{operation.get('workstation', 'its workstation')}: {result.get('error') or 'operation failed'}"
                )
            else:
                remaining_operations.append(operation)
        if remaining_operations:
            st.session_state["cluster-results-purge-operations"] = remaining_operations
            st.info(f"Physically removing {len(remaining_operations)} run folder(s) on the owning workstation(s)…")
        else:
            st.session_state.pop("cluster-results-purge-operations", None)
            if completed_operations:
                st.success(f"Physically removed {completed_operations} remote run folder(s) and their automatic Database records.")
    show_removed = st.toggle("Show removed jobs", key="cluster-results-show-removed")
    # Keep hidden rows in the merged view so local jobs do not reappear when
    # they are filtered out by the coordinator's default Results query.
    page, stale = catalog(user, token, include_hidden=True)
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
        if _is_bindcraft2_job(item):
            item["progress"] = _row_progress(item)
        if key in known:
            rows[known[key]] = rows[known[key]] | item
        else:
            rows.append(item)
    status_options = sorted({str(r.get("status") or "unknown").replace("_", " ").title() for r in rows}, key=str.casefold)
    status_key = "cluster-results-status-filter"
    if not set(st.session_state.get(status_key, [])).issubset(status_options):
        st.session_state[status_key] = []
    search_col, machine_col, status_col, refresh_col = st.columns([0.38, 0.26, 0.24, 0.12])
    query = search_col.text_input("Search jobs", placeholder="Name, engine, or run ID", key="cluster-results-query").casefold()
    names = {r.get("node_id"): config.display_name(r.get("node_name") or "Not assigned") for r in rows if r.get("node_id")}
    machine = machine_col.selectbox("Workstation", ["", *names], format_func=lambda v: names.get(v, "All workstations"), key="cluster-results-machine")
    statuses = status_col.multiselect("Status", status_options, key=status_key, placeholder="All statuses")
    if refresh_col.button("Refresh", key="cluster-results-refresh"):
        st.rerun()
    projects = sorted({r.get("project_id") for r in rows if r.get("project_id")})
    from mn_molecular_design_studio.workspace import load_organization
    project_names = {p["id"]: p["name"] for p in load_organization().get("projects", [])}
    project = st.selectbox("Project", ["", *projects], format_func=lambda v: project_names.get(v, "All projects"), key="cluster-results-project")
    chosen = [r for r in rows if (show_removed or not r.get("hidden"))
              and (not machine or r.get("node_id") == machine) and (not project or r.get("project_id") == project)
              and (not statuses or str(r.get("status") or "unknown").replace("_", " ").title() in statuses)
              and (not query or query in " ".join(f"{r.get(k, '')} {config.display_text(r.get(k, ''))}" for k in ("campaign_name", "pipeline_name", "run_id", "tool", "job_code")).casefold())]
    chosen.sort(key=lambda r: r.get("created_at", ""), reverse=True)
    st.caption(f"{len(chosen)} jobs · Job status and file availability are tracked separately.")
    pages = max(1, (len(chosen) + 49) // 50)
    page_key = "cluster-results-page"
    if st.session_state.get(page_key, 1) > pages:
        st.session_state[page_key] = 1
    selected_page = int(st.number_input("Page", min_value=1, max_value=pages, value=1, step=1, key=page_key))
    header = st.columns([0.35, 2.8, 1, 1.25, 1, 1.3, 0.7])
    header[0].caption("Select")
    for column, label in zip(
        header[1:],
        ("Job name · tool · pipeline", "Status", "Elapsed", "Workstation", "File availability", "Open"),
    ):
        column.caption(label)
    selected_jobs = st.session_state.setdefault("cluster-results-selected", {})
    for row in chosen[(selected_page - 1) * 50:selected_page * 50]:
        ref = str(row.get("cluster_job_id") or (str(row.get("task_group") or "") + "/" + str(row.get("run_id") or "")))
        checkbox_key = "cluster-results-select-" + hashlib.sha1(ref.encode()).hexdigest()[:12]
        with st.container(border=True):
            select_col, title, state, elapsed, machine, access, action = st.columns([0.35, 2.8, 1, 1.25, 1, 1.3, 0.7])
            if select_col.checkbox("Select job", value=ref in selected_jobs, key=checkbox_key, label_visibility="collapsed"):
                selected_jobs[ref] = row
            else:
                selected_jobs.pop(ref, None)
            title.write(config.display_text(row.get("campaign_name") or row.get("tool") or "Job"))
            pipeline_name = str(row.get("pipeline_name") or "").strip() or "—"
            title.caption(f"Tool: {display_tool_label(row)} · Pipeline: {config.display_text(pipeline_name)}")
            title.caption(project_names.get(row.get("project_id"), "Project") + " · " + str(row.get("job_code", "")))
            progress_label = _progress_summary(_row_progress(row))
            if progress_label:
                title.caption(progress_label)
            state.write(str(row.get("status", "unknown")).replace("_", " ").title())
            elapsed.write(elapsed_label(row) if elapsed_label else "—")
            since = elapsed_since(row) if elapsed_since else ""
            if since:
                elapsed.caption(since)
            machine.write(config.display_name(row.get("node_name") or "Not assigned"))
            availability = "last known" if stale and row.get("node_id") != config.identity() else row.get("availability", "pending")
            access.write(availability.title())
            if row.get("last_seen_at"):
                access.caption(str(row["last_seen_at"]))
            view_ref = "cluster:" + row["cluster_job_id"] if row.get("cluster_job_id") else row["run_id"]
            if row.get("node_id") == config.identity() and row.get("run_id"):
                view_ref = row["run_id"]
            view_href = "?job_ref=" + quote(str(view_ref), safe="")
            action.markdown(
                f'<a href="{view_href}" target="_self" '
                'style="display:inline-block;padding:0.35rem 0.65rem;border:1px solid #d5dbe1;'
                'border-radius:0.4rem;background:#fff;color:#26352d;text-decoration:none">View</a>',
                unsafe_allow_html=True,
            )
    st.session_state["cluster-results-selected"] = selected_jobs
    selected_rows = list(selected_jobs.values())
    if selected_rows:
        st.caption(f"{len(selected_rows)} job{'s' if len(selected_rows) != 1 else ''} selected")
        if show_removed:
            available_to_restore = [row for row in selected_rows if row.get("hidden")]
            if available_to_restore and st.button("Restore selected jobs", key="cluster-results-restore-selected"):
                failures = []
                restored = 0
                for row in available_to_restore:
                    try:
                        jid = str(row.get("cluster_job_id") or "")
                        if jid:
                            _request("/jobs/" + jid + "/actions", token, {"action": "restore"}, method="POST")
                        elif row.get("node_id") == config.identity() and row.get("task_group") and row.get("run_id"):
                            from mn_molecular_design_studio.auth import can_access_job
                            from mn_molecular_design_studio.job_catalog import set_job_hidden, get_job_access_record
                            access_record = get_job_access_record(str(row["task_group"]), str(row["run_id"]))
                            if not can_access_job(user, access_record):
                                raise PermissionError("Job not found.")
                            set_job_hidden(str(row["task_group"]), str(row["run_id"]), False)
                        else:
                            raise ValueError("This job cannot be updated while its workstation is unavailable.")
                        restored += 1
                    except (ValueError, PermissionError, transport.ClusterUnavailable) as exc:
                        failures.append(str(exc))
                st.session_state["cluster-results-selected"] = {}
                st.session_state["cluster-results-reset-selection"] = True
                if restored:
                    st.session_state["cluster-results-flash"] = f"Restored {restored} job{'s' if restored != 1 else ''}."
                for failure in failures:
                    st.error(failure)
                st.rerun()
        else:
            pending_key = "cluster-results-pending-remove-and-purge"
            if st.button("Remove and purge selected", key="cluster-results-remove-and-purge-request"):
                st.session_state[pending_key] = [dict(row) for row in selected_rows]
            pending = st.session_state.get(pending_key)
            if pending:
                st.warning("Confirming permanently deletes each selected finished run folder from its workstation's Results directory and removes its automatic Database records. Queued or running jobs cannot be purged. Molecules explicitly saved to user collections are kept. Remote runs are deleted when their workstation receives the request.")
                confirm, cancel = st.columns(2)
                if confirm.button("Confirm removal and purge", type="primary", key="cluster-results-remove-and-purge-confirm"):
                    removed_sources = removed_members = succeeded = 0
                    queued_operations = []
                    failures = []
                    for row in pending:
                        try:
                            jid = str(row.get("cluster_job_id") or "")
                            if jid:
                                result = _request("/jobs/" + jid + "/purge", token, method="DELETE")
                                if result.get("status") == "purge_queued":
                                    queued_operations.append({
                                        "command_id": result["command_id"],
                                        "workstation": config.display_name(row.get("node_name") or result.get("node_id")),
                                        "name": row.get("campaign_name") or row.get("run_id") or jid,
                                    })
                                else:
                                    database = result.get("database") or {}
                                    removed_sources += int(database.get("removed_sources", 0))
                                    removed_members += int(database.get("removed_members", 0))
                                    succeeded += 1
                            elif row.get("node_id") == config.identity() and row.get("task_group") and row.get("run_id"):
                                from mn_molecular_design_studio.auth import can_access_job
                                from mn_molecular_design_studio.job_catalog import get_job_access_record, reconcile
                                from mn_molecular_design_studio.molecular_database_results import remove_run_records
                                from mn_molecular_design_studio.engine.workflow_engine.core.jobs import delete_job_run, runs_root
                                access_record = get_job_access_record(str(row["task_group"]), str(row["run_id"]))
                                if not can_access_job(user, access_record):
                                    raise PermissionError("Job not found.")
                                if str(row.get("status") or "").casefold() not in {"completed", "succeeded", "cancelled", "failed"}:
                                    raise ValueError("Only finished jobs can be permanently removed.")
                                result = remove_run_records(
                                    row["task_group"], row["run_id"], user,
                                    owner_id=row.get("owner_id"), project_id=row.get("project_id"),
                                )
                                removed_sources += result["removed_sources"]
                                removed_members += result["removed_members"]
                                try:
                                    delete_job_run(row["task_group"], row["run_id"])
                                except FileNotFoundError:
                                    pass
                                reconcile(runs_root())
                                succeeded += 1
                            else:
                                raise ValueError("This job cannot be updated while its workstation is unavailable.")
                        except (ValueError, PermissionError, OSError, sqlite3.Error, transport.ClusterUnavailable) as exc:
                            failures.append(f"{row.get('campaign_name') or row.get('run_id') or 'Selected job'}: {exc}")
                    st.session_state.pop(pending_key, None)
                    st.session_state["cluster-results-selected"] = {}
                    st.session_state["cluster-results-reset-selection"] = True
                    if queued_operations:
                        st.session_state["cluster-results-purge-operations"] = queued_operations
                    messages = []
                    if succeeded:
                        messages.append(
                            f"Physically removed {succeeded} run folder{'s' if succeeded != 1 else ''} "
                            f"and their automatic Database records "
                            f"({removed_sources} source links, {removed_members} run records)."
                        )
                    if queued_operations:
                        messages.append(
                            f"Queued physical deletion of {len(queued_operations)} run folder"
                            f"{'s' if len(queued_operations) != 1 else ''} on the owning workstation(s)."
                        )
                    if messages:
                        st.session_state["cluster-results-flash"] = " ".join(messages)
                    if failures:
                        st.session_state["cluster-results-purge-errors"] = failures
                    st.rerun()
                if cancel.button("Cancel", key="cluster-results-remove-and-purge-cancel"):
                    st.session_state.pop(pending_key, None)
                    st.rerun()
        clear_col, _ = st.columns([1, 4])
        if clear_col.button("Clear selection", key="cluster-results-clear-selection"):
            st.session_state["cluster-results-selected"] = {}
            st.session_state["cluster-results-reset-selection"] = True
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
    if _is_bindcraft2_job(row):
        _render_bindcraft2_progress(row, str(row.get("cluster_job_id") or jid), token)
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
    if row.get("metrics"):
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
                tool = " ".join(
                    str(row.get(key) or "")
                    for key in ("tool", "job_type", "workflow", "pipeline_name", "campaign_name")
                ).casefold() + " " + display_tool_label(row).casefold()
                if "rfpeptides" in tool:
                    st.info("RFpeptides generated cyclic peptide backbones. It does not design amino-acid sequences; run ProteinMPNN or another sequence-design tool as a follow-up.")
                rows = _candidate_result_rows(candidates, tool=tool)
                if "bindcraft2" in tool:
                    st.caption("The table includes BindCraft2's ranked confidence, interface, clash, and developability metrics. Scroll horizontally to see all columns.")
                st.dataframe(
                    pd.DataFrame(rows),
                    hide_index=True,
                    use_container_width=True,
                    height=min(410, 38 * (len(rows) + 1) + 3),
                    column_config={"Sequence": st.column_config.TextColumn(width="large")},
                )
                parameter_rows, filter_rows = ([], [])
                if "bindcraft2" in tool:
                    parameter_rows, filter_rows = _bindcraft2_run_parameter_rows(data)
                selected = st.selectbox(
                    "Design to visualize",
                    range(len(candidates)),
                    format_func=lambda index: f"{index + 1}. {candidates[index].get('candidate_id') or 'Design'}",
                    key="cluster-candidate-" + jid,
                )
                candidate = candidates[selected]
                sequence = candidate.get("binder_sequence") or candidate.get("sequence")
                if sequence:
                    st.markdown("**Designed molecule sequence**")
                    st.code(sequence, language=None)
                elif "rfpeptides" in tool:
                    st.caption("Sequence: not applicable yet; this run contains a backbone structure only.")

                structure = _candidate_structure_path(candidate, data)
                if not structure:
                    st.info("No candidate structure file was recorded for this design.")
                elif row.get("availability") != "online":
                    st.info("The structure is saved on the owning workstation and will load when that workstation reconnects.")
                else:
                    relative_path, structure_field = structure
                    candidate_key = str(candidate.get("candidate_id") or selected)
                    file_key = f"cluster-candidate-file-{jid}-{candidate_key}"
                    command_id = st.session_state.get(file_key)
                    if not command_id:
                        try:
                            file_operation = _request(
                                "/jobs/" + jid + "/actions",
                                token,
                                {"action": "artifact", "parameters": {"path": relative_path}},
                                method="POST",
                            )
                            command_id = file_operation["command_id"]
                            st.session_state[file_key] = command_id
                        except (ValueError, transport.ClusterUnavailable) as exc:
                            st.warning(f"Could not request the representative structure: {exc}")
                    if command_id:
                        try:
                            file_status = _request("/commands/" + command_id, token)
                            if file_status["status"] == "failed":
                                st.warning(file_status.get("error") or "The workstation could not serve this structure.")
                            elif file_status["status"] != "completed":
                                st.info("Loading the representative structure from the owning workstation…")
                            else:
                                metadata = file_status.get("response") or {}
                                name = metadata.get("name") or Path(relative_path).name
                                if metadata.get("size_bytes", 0) > 50 * 1024 * 1024:
                                    st.info(f"{name} is too large for the inline structure viewer. Use the saved-file controls below to download it.")
                                else:
                                    content_key = "cluster-candidate-bytes-" + command_id
                                    if content_key not in st.session_state:
                                        st.session_state[content_key] = _request(
                                            "/commands/" + command_id + "/file",
                                            token,
                                            raw=True,
                                            timeout=60,
                                        )
                                    content = st.session_state[content_key]
                                    structure_text = content.decode("utf-8", errors="replace")
                                    st.subheader("Representative structure")
                                    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
                                    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization
                                    visualizations = [StructureVisualization(
                                        pdb=structure_text,
                                        color=_candidate_structure_color(candidate, structure_text),
                                        representation_type="cartoon",
                                    )]
                                    overlay_key = None
                                    if "bindcraft2" in tool:
                                        candidate_token = hashlib.sha256(candidate_key.encode()).hexdigest()[:12]
                                        overlay_key = f"cluster-bc2-overlay-{jid}-{candidate_token}"
                                        if st.checkbox(
                                            "Overlay binder-alone monomer",
                                            key=overlay_key + "-enabled",
                                            help="Align the saved binder-only prediction onto the binder chain in this complex.",
                                        ):
                                            monomer_path = _bindcraft2_monomer_path("", candidate)
                                            if not monomer_path:
                                                refolded_csv = _fetch_remote_artifact_text(
                                                    jid,
                                                    _bindcraft2_refolded_table_path(data),
                                                    token,
                                                    overlay_key + "-refolded-csv",
                                                    label="BindCraft2 refolded results table",
                                                )
                                                if refolded_csv:
                                                    monomer_path = _bindcraft2_monomer_path(refolded_csv, candidate)
                                            if monomer_path:
                                                monomer_token = hashlib.sha256(monomer_path.encode()).hexdigest()[:12]
                                                monomer_text = _fetch_remote_artifact_text(
                                                    jid,
                                                    monomer_path,
                                                    token,
                                                    overlay_key + "-monomer-" + monomer_token,
                                                    label="binder-alone monomer",
                                                )
                                                if monomer_text:
                                                    aligned_key = overlay_key + "-aligned-" + monomer_token
                                                    if aligned_key not in st.session_state:
                                                        try:
                                                            aligned_text, overlay_rmsd, matched_residues = _align_bindcraft2_monomer_to_complex(
                                                                structure_text, monomer_text, candidate
                                                            )
                                                            st.session_state[aligned_key] = {
                                                                "text": aligned_text,
                                                                "rmsd": overlay_rmsd,
                                                                "matched_residues": matched_residues,
                                                            }
                                                        except (ImportError, KeyError, ValueError, OSError) as exc:
                                                            st.warning(f"Could not align the binder-alone monomer: {exc}")
                                                    aligned = st.session_state.get(aligned_key)
                                                    if aligned:
                                                        visualizations.append(StructureVisualization(
                                                            pdb=aligned["text"],
                                                            color="uniform",
                                                            color_params={"value": "0xff7d45"},
                                                            representation_type="cartoon",
                                                        ))
                                                        st.caption(
                                                            "Orange: binder-alone monomer, fitted onto the complex binder. "
                                                            f"Cα RMSD: {aligned['rmsd']:.2f} Å over {aligned['matched_residues']} matched residues."
                                                        )
                                            else:
                                                st.warning("No matching binder-alone monomer file was found for this accepted sequence.")
                                    molstar_custom_component(
                                        visualizations,
                                        key=f"cluster-candidate-view-{jid}-{selected}-{structure_field}-overlay-{bool(overlay_key and overlay_key + '-enabled' in st.session_state and st.session_state[overlay_key + '-enabled'])}",
                                        height=520,
                                        show_controls=True,
                                        show_structure_tools=False,
                                    )
                                    st.download_button(
                                        "Download candidate structure",
                                        content,
                                        file_name=name,
                                        key=f"cluster-candidate-download-{command_id}",
                                    )
                        except (ValueError, transport.ClusterUnavailable) as exc:
                            st.info(str(exc))
                if parameter_rows:
                    st.markdown("**BindCraft2 run parameters**")
                    st.dataframe(
                        pd.DataFrame(parameter_rows),
                        hide_index=True,
                        use_container_width=True,
                        height=min(360, 35 * (len(parameter_rows) + 1) + 3),
                    )
                if filter_rows:
                    st.markdown("**Quality filters used**")
                    st.dataframe(
                        pd.DataFrame(filter_rows),
                        hide_index=True,
                        use_container_width=True,
                        height=min(260, 35 * (len(filter_rows) + 1) + 3),
                    )
            with st.expander("Run settings and inputs"):
                st.json(data.get("input", {}))
            with st.expander("Run summary"):
                st.json(data.get("result", {}))
            if candidates and "rfantibody" in tool:
                from mn_molecular_design_studio.rfantibody_ui import render_confidence_guide
                render_confidence_guide()
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
