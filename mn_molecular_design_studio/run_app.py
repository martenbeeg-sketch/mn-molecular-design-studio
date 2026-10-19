from __future__ import annotations

import csv
import hashlib
from html import escape
import io
import json
import importlib
import os
import re
import secrets
import sqlite3
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas as pd
import streamlit as st

from mn_molecular_design_studio.cofolding_input import render_molecule_editor, prediction_settings_for
from mn_molecular_design_studio import cofolding_results
from mn_molecular_design_studio.cluster import config as cluster_config
from mn_molecular_design_studio.cluster.ui import display_tool_label
from mn_molecular_design_studio import pipeline_graph, complexa_input, generator_input
importlib.reload(generator_input)
importlib.reload(pipeline_graph)
importlib.reload(complexa_input)
WORKFLOW_ROLES = pipeline_graph.WORKFLOW_ROLES
from mn_molecular_design_studio.design_validation_ui import render_batch_controls, render_comparison

from mn_molecular_design_studio.workspace import (
    ensure_studio_workspace,
    is_studio_job,
    load_organization,
    load_settings,
    mark_studio_job,
    save_settings,
    create_project,
    validate_project_id,
)
from mn_molecular_design_studio.auth import (
    accessible_project_ids,
    authenticate,
    can_access_job,
    change_password,
    create_user,
    development_user,
    grant_project_access,
    is_development_mode,
    list_users,
    password_min_length,
    project_member_ids,
    require_project_access,
    reset_user_password,
    resolve_session,
    revoke_project_access,
    revoke_session,
    set_user_disabled,
    user_input_root,
    user_count,
)

STUDIO_WORKSPACE = ensure_studio_workspace()

from mn_molecular_design_studio.runtime import app_home, reference_root, runs_root
from mn_molecular_design_studio.job_catalog import (
    catalog_path,
    get_job_access_record,
    job_artifacts,
    job_list_facets as _catalog_job_list_facets,
    search_jobs as _catalog_search_jobs,
)
from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import (
    ACTIVE_STATUSES,
    AutomationError,
    cancel_job,
    inspect_job,
    job_results,
    resolve_job,
    submit_design_campaign_request,
)
from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import build_run_zip
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import STAGE_COMPLEX_REFOLDING, read_candidates
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
# Refresh catalog additions in an already running Streamlit session.
importlib.reload(workflow_automation)
workflow_catalog = workflow_automation.workflow_catalog
workflow_schema = workflow_automation.workflow_schema
from mn_molecular_design_studio.engine.workflow_engine.workflows import (
    af3_jax_cofolding,
    bindcraft2,
    candidate_import,
    design,
    generators,
    esmfold2_cofolding,
    optimized_prediction,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows.design_campaigns import ENGINE_LABELS, ENGINE_ORDER
from mn_molecular_design_studio.engine.workflow_engine.workflows.modules import candidate_sources, load_source_candidates
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization
from mn_molecular_design_studio import target_editor
importlib.reload(target_editor)
from mn_molecular_design_studio.target_editor import (
    TargetDesignInput,
    TargetSource,
    _sequence_entries_by_chain,
    render_target_coldspots,
    render_target_editor,
    render_target_source,
)

# Reload the Studio-owned BindCraft2 workflow on Streamlit reruns so workflow
# edits take effect without dropping form state.
importlib.reload(af3_jax_cofolding)
importlib.reload(bindcraft2)
importlib.reload(design)
importlib.reload(generators)

API_BASE_URL = (os.environ.get("MDS_API_BASE_URL") or load_settings(STUDIO_WORKSPACE)["api_base_url"]).rstrip("/")

WORKFLOW_TITLES = {
    "Esm If1": "ESM-IF1",
    "Caliby": "Caliby",
    "Boltzgen Inverse": "BoltzGen Inverse Folding",
    "Pxdesign Generate": "PXDesign",
    "Rfdiffusion Generate": "RFdiffusion",
    "Genie3 Generate": "Genie3",
    "Protpardelle Generate": "Protpardelle-1C",
    "Proteina Complexa": "Proteina-Complexa",
    "Proteina Complexa Generate": "Proteina-Complexa Generator",
    "Bindcraft2": "BindCraft2",
    "Bindcraft2 Generator": "BindCraft2 Generator",
    "Rfantibody": "RFantibody",
    "Rfpeptides": "RFpeptides",
    "Bindcraft": "BindCraft",
    "Boltzgen": "BoltzGen",
    "Boltzgen Generator": "BoltzGen Generator",
    "Af2 Initial Guess": "AlphaFold2 Initial Guess",
    "Esmfold2": "ESMFold2",
    "Esmfold2 Cofolding": "ESMFold2 Cofolding",
    "Esmfold Monomer": "ESMFold",
    "Af3 Jax Cofolding": "ColabFold 2.0 preview",
    "Rf3": "RosettaFold3",
    "Optimized Prediction": "Structure Prediction",
    "Esmfold2 Binder": "ESMFold2 Binder Design",
    "Boltz2 Complex": "Boltz-2 Complex",
    "Boltz2 Monomer": "Boltz-2 Monomer",
    "Ligandmpnn": "LigandMPNN",
    "Mpnn": "ProteinMPNN",
    "Rfdiffusion Classic": "RFdiffusion",
    "Rfdiffusion3 Foundry": "RFdiffusion3",
    "Rfdiffusion3": "RFdiffusion3",
    "Ipsae": "ipSAE Scoring",
    "Rosetta Metrics": "Rosetta Interface Metrics",
    "Pesto": "PeSTo",
    "Scannet": "ScanNet",
    "Masif Seed": "MaSIF-seed",
    "Surf2Spot": "Surf2Spot",
    "Surface Analysis": "Protein Surface Analysis",
}

# Keep the Tools catalog focused on workflows that launch a specific tool or
# model. App-level preparation, import, ranking, benchmarking, and orchestration
# workflows remain available to the workflow API without appearing as tool cards.
TOOL_CATALOG_WORKFLOW_IDS = {
    "detection.scannet",
    "detection.pesto",
    "detection.surf2spot",
    "detection.masif_seed",
    "detection.surface_analysis",
    "design.rfdiffusion_generate",
    "design.bindcraft2",
    "design.bindcraft2_generator",
    "design.rfdiffusion3",
    "design.rfantibody",
    "design.rfpeptides",
    "design.boltzgen",
    "design.boltzgen_generator",
    "design.pxdesign_generate",
    "design.genie3_generate",
    "design.esmfold2_binder",
    "design.protpardelle_generate",
    "design.proteina_complexa",
    "design.proteina_complexa_generate",
    "sequence.mpnn",
    "sequence.esm_if1",
    "sequence.caliby",
    "sequence.boltzgen_inverse",
    "refolding.af3_jax_cofolding",
    "refolding.optimized_prediction",
    "refolding.af2_initial_guess",
    "refolding.openfold3",
    "refolding.rf3",
    "refolding.protenix",
    "refolding.protenix_cli",
    "refolding.boltzgen_fold",
    "analysis.ipsae",
    "analysis.rosetta",
}

# Upstream code versions for tool tiles. Values describe engine code or a pinned
# source revision, never the Docker image tag. Some legacy shared images have no
# pinned upstream release, which is shown explicitly.
TOOL_WORKFLOW_CODE_VERSIONS = {
    "sequence.esm_if1": ("ESM-IF1 · fair-esm 2.0.1 · 2b369911", False),
    "sequence.caliby": ("Caliby 0.1 · 41d31560", True),
    "sequence.boltzgen_inverse": ("BoltzGen 0.3.2 · a07fea1", False),
    "design.pxdesign_generate": ("PXDesign f788441", False),
    "design.rfdiffusion_generate": ("RFdiffusion · version not pinned", False),
    "design.genie3_generate": ("Genie3 d77ae5a", False),
    "design.protpardelle_generate": ("Protpardelle-1C ee37840", False),
    "detection.scannet": ("ScanNet a61623c", False),
    "detection.pesto": ("PeSTo ba651aa", False),
    "detection.surf2spot": ("Surf2Spot a38a7a7", False),
    "detection.masif_seed": ("MaSIF-seed 7715391", False),
    "detection.surface_analysis": ("PeSTo · ScanNet · Surf2Spot · MaSIF-seed", False),
    "design.bindcraft2": ("BindCraft 2 c0ee638", False),
    "design.bindcraft2_generator": ("BindCraft 2 c0ee638", False),
    "design.rfdiffusion_classic": ("RFdiffusion · version not pinned", False),
    "design.rfdiffusion3": ("RFdiffusion3 a36d29c · Foundry 0.1.9", False),
    "design.rfantibody": ("RFantibody 3e5a5bb · CUDA 11.8", False),
    "design.rfpeptides": ("RFpeptides · RFdiffusion 86507b6 · CUDA 12.8", False),
    "design.boltzgen": ("BoltzGen 0.3.2", True),
    "design.boltzgen_generator": ("BoltzGen 0.3.2", True),
    "design.pxdesign": ("PXDesign f788441 · Protenix 0.5.0+pxd", False),
    "design.genie3": ("Genie3 d77ae5a", False),
    "design.esmfold2_binder": ("Biohub ESM 3.4.1.post1", False),
    "design.protpardelle_1c": ("Protpardelle-1C ee37840", False),
    "design.proteina_complexa": ("Proteina-Complexa 1.0.0 · 32b71ae", False),
    "design.proteina_complexa_generate": ("Proteina-Complexa 1.0.0 · 32b71ae", False),
    "sequence.ligandmpnn": ("LigandMPNN · version not pinned", False),
    "sequence.foundry_mpnn": ("Foundry MPNN · version not pinned", False),
    "sequence.mpnn": ("LigandMPNN 26ec57a · Foundry 0.1.9", False),
    "refolding.af3_jax_cofolding": (
        "ColabFold AF3 preview f850f0f · AlphaFold3 3.1.14",
        True,
    ),
    "refolding.af2_initial_guess": ("AlphaFold2 · version not pinned", False),
    "refolding.openfold3": ("OpenFold3 0.5.0 · P2 weights", False),
    "refolding.rf3": ("RosettaFold3 · version not pinned", False),
    "refolding.protenix": ("Protenix 0.5.0+pxd", False),
    "refolding.protenix_cli": ("Protenix 2.0.0", False),
    "refolding.boltzgen_fold": ("BoltzGen 0.3.2", True),
    "analysis.ipsae": ("ipSAE · version not pinned", False),
    "analysis.rosetta": ("PyRosetta · version not pinned", False),
}

# Related prediction modes share one Tools card. The existing workflows remain
# separate behind the card because their inputs and settings differ.
TOOL_CATALOG_GROUPS = {
    "boltzgen_fold": {
        "title": "BoltzGen Fold",
        "description": "Candidate-batch refolding or target PDB plus binder sequence.",
        "engine_version": "BoltzGen 0.3.2",
        "engine_optimized": True,
        "modes": {
            "Refold candidate batch": "refolding.boltzgen_fold",
            "Target PDB + binder sequence": "refolding.optimized_prediction",
        },
        "engine_profiles": {"Target PDB + binder sequence": "boltzgen_fold"},
    },
}

ANALYSIS_SOURCE_WORKFLOWS = {"analysis.ipsae", "analysis.rosetta"}

# Each optimized prediction profile gets its own engine card. Version labels
# below refer to the upstream engine code, not the container or batch wrapper.
OPTIMIZED_IMAGE_TILE_IDS = tuple(optimized_prediction.ENGINE_PROFILES)

BINDCRAFT2_MODALITY_LABELS = {
    "binder": "De novo binder",
    "large_binder": "Large binder",
    "peptide": "Peptide",
    "cyclic_peptide": "Cyclic peptide",
    "VHH": "VHH (single-domain antibody)",
    "scFv": "scFv",
    "Fab": "Fab",
    "ARP": "ARP (ankyrin repeat)",
    "multidomain": "Multidomain",
    "homo_oligomer": "Homo-oligomer",
}
BINDCRAFT2_MODALITY_HELP = {
    "binder": "De novo miniprotein, with binder length sampled from the selected range.",
    "large_binder": "Long de novo binder (250–600 residues), with longer optimization stages.",
    "peptide": "Linear peptide binder (12–25 residues); complex pTM and free-binder pLDDT filters are disabled by the preset.",
    "cyclic_peptide": "Head-to-tail cyclic peptide (7–16 residues), with peptide-specific validation and closure.",
    "VHH": "VHH scaffold template; framework is held and CDR/edit spans are designed. Loop lengths come from the edit ranges.",
    "scFv": "Two-chain scFv scaffold; antibody loops are designed and variable loop lengths are sampled from the edit ranges.",
    "Fab": "Fab scaffold template; antibody loops are designed while framework residues are held.",
    "ARP": "Ankyrin-repeat scaffold; repeat loop lengths come from the edit ranges.",
    "multidomain": "One long binder chain assembled from multiple designed domains.",
    "homo_oligomer": "A binder designed as a symmetric assembly; copies sets the oligomer size.",
}
BINDCRAFT2_CONFORMATIONAL_OBJECTIVE_LABELS = {
    "none": "None",
    "induced_fit": "Induced Fit",
    "fold_switch": "Fold Switch",
}
BINDCRAFT2_CONFORMATIONAL_OBJECTIVE_HELP = {
    "none": "No extra objective; use the selected binder format's normal design objective.",
    "induced_fit": "The preset aims for 5 Å of interface movement between free and target-bound predictions.",
    "fold_switch": "The preset aims for a whole-fold change (TM-score ≤0.6) between free and target-bound predictions.",
}
BINDCRAFT2_SHORT_LABELS = {
    "binder": "Binder",
    "large_binder": "Large binder",
    "peptide": "Peptide",
    "cyclic_peptide": "Cyclic",
    "VHH": "VHH",
    "scFv": "scFv",
    "Fab": "Fab",
    "ARP": "ARP",
    "multidomain": "Multidomain",
    "homo_oligomer": "Homo-olig.",
}
BINDCRAFT2_PROPERTY_LABELS = {
    "humanize": "Humanize",
    "forced_targeting": "Forced Targeting",
    "termini_accessible": "Accessible Termini",
    "termini_together": "Termini Together",
    "disulfide_staple": "Disulfide Staple",
    "protease_stable": "Protease Stable",
    "mixed_topology": "Mixed Topology",
    "bigbang": "Big Bang Start",
    "initial_guess": "Initial Guess",
}
BINDCRAFT2_PROPERTY_HELP = {
    "humanize": "For therapeutic leads. Pushes the sequence toward human germline and away from MHC-II anchor motifs.",
    "forced_targeting": "Forces contact onto your hotspots. Needs hotspots set above, and a structured target.",
    "termini_accessible": "Keeps both chain ends pointing away from the target, so you can add tags or fusions.",
    "termini_together": "Brings the N and C ends close together, within about 10 angstroms.",
    "disulfide_staple": "Allows cysteine and requires at least one real disulfide in the design.",
    "protease_stable": "Selects against protease cleavage sites, exposed loops and exposed ends, so the binder is harder to cut.",
    "mixed_topology": "Pushes for beta sheet instead of an all-helical bundle. At most 50% helix, at least 20% sheet.",
    "bigbang": "Seeds the design from coordinates on hand. Only helps on large complexes, over about 600 residues.",
    "initial_guess": "Re-predicts each candidate from the pose it folded into. Helps on hard interfaces, but raises false positives.",
}
BINDCRAFT2_PROPERTY_ORDER = (
    "humanize",
    "forced_targeting",
    "termini_accessible",
    "termini_together",
    "disulfide_staple",
    "protease_stable",
    "mixed_topology",
    "bigbang",
    "initial_guess",
)
BINDCRAFT2_BASE_QUALITY_FILTERS = {
    "min_interface_iptm": 0.7,
    "min_complex_ptm": 0.55,
    "max_interface_ipae": 0.35,
    "min_unbound_binder_plddt": 0.8,
    "min_interface_residues": 7,
    "max_backbone_clashes": 0,
    "max_binder_rmsd": None,
    "min_plddt_screen": 0.6,
    "min_plddt_refine": 0.6,
    "min_plddt_anneal": 0.65,
    "min_plddt_harden": 0.65,
    "min_plddt_mutate": 0.6,
    "min_plddt_final": 0.7,
    "min_iptm_anneal": 0.5,
    "min_iptm_harden": 0.5,
    "min_iptm_mutate": 0.5,
    "max_interdomain_contact_final": None,
    "min_domain_separation_ratio_final": None,
    "min_scaffold_sequence_retained_final": None,
    "max_scaffold_framework_rmsd_final": None,
}
BINDCRAFT2_MODALITY_QUALITY_FILTERS = {
    "binder": {"min_unbound_binder_plddt": 0.7, "max_binder_rmsd": 3.5},
    "large_binder": {
        "min_unbound_binder_plddt": 0.6,
        "max_binder_rmsd": 3.5,
    },
    "peptide": {
        "min_complex_ptm": None,
        "max_interface_ipae": 0.33,
        "min_unbound_binder_plddt": None,
    },
    "cyclic_peptide": {
        "min_complex_ptm": None,
        "max_interface_ipae": 0.33,
        "min_unbound_binder_plddt": None,
    },
    "homo_oligomer": {"min_unbound_binder_plddt": 0.7, "max_binder_rmsd": 3.5},
    "multidomain": {
        "min_complex_ptm": 0.4,
        "max_interface_ipae": 0.45,
        "min_unbound_binder_plddt": 0.7,
        "min_plddt_final": 0.75,
        "max_interdomain_contact_final": 0.05,
        "min_domain_separation_ratio_final": 1.0,
    },
    "VHH": {
        "min_unbound_binder_plddt": 0.7,
        "min_scaffold_sequence_retained_final": 1.0,
    },
    "ARP": {
        "min_unbound_binder_plddt": 0.7,
        "min_scaffold_sequence_retained_final": 1.0,
        "max_scaffold_framework_rmsd_final": 3.0,
    },
    "scFv": {
        "min_unbound_binder_plddt": 0.7,
        "min_scaffold_sequence_retained_final": 1.0,
        "max_scaffold_framework_rmsd_final": 3.0,
    },
    "Fab": {
        "min_unbound_binder_plddt": 0.7,
        "min_scaffold_sequence_retained_final": 1.0,
        "max_scaffold_framework_rmsd_final": 3.0,
    },
}


def _bindcraft2_quality_filter_defaults(modality: str) -> dict[str, float | int | None]:
    defaults = dict(BINDCRAFT2_BASE_QUALITY_FILTERS)
    defaults.update(BINDCRAFT2_MODALITY_QUALITY_FILTERS.get(modality, {}))
    return defaults


st.set_page_config(
    page_title="Molecular Design Studio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


def _render_authentication_gate() -> dict:
    if is_development_mode() and not cluster_config.is_worker():
        return development_user()
    try:
        configured_users = user_count()
    except (OSError, sqlite3.Error) as exc:
        st.error(f"Could not open the Studio account database: {exc}")
        st.stop()

    if configured_users == 0:
        st.title("Set up Molecular Design Studio")
        st.warning("No admin account exists yet. Create your admin account from the computer running the Studio.")
        st.code("mn-molecular-design-studio create-admin", language="bash")
        st.caption("The command prompts for your username and a password without displaying the password as you type.")
        st.stop()

    token = st.session_state.get("studio-session-token")
    user = resolve_session(token)
    if user is None:
        if token:
            revoke_session(token)
            # Clear cached job rows and result data when a session expires.
            for key in list(st.session_state.keys()):
                del st.session_state[key]
        st.title("Sign in to Molecular Design Studio")
        with st.form("studio-login-form"):
            username = st.text_input("Username or email")
            password = st.text_input("Password", type="password")
            submitted = st.form_submit_button("Sign in", type="primary", use_container_width=True)
        if submitted:
            result = authenticate(username, password)
            if result is None:
                st.error("Sign-in failed. Check your credentials or ask the admin to restore your account.")
            else:
                user, token, _expires_at = result
                st.session_state["studio-session-token"] = token
                st.rerun()
        st.stop()

    if user.get("must_change_password"):
        st.title("Set a new password")
        st.info("The admin gave you a temporary password. Choose a new password before continuing.")
        with st.form("studio-forced-password-form"):
            new_password = st.text_input("New password", type="password")
            confirm_password = st.text_input("Confirm new password", type="password")
            submitted = st.form_submit_button("Save password", type="primary")
        if submitted:
            if new_password != confirm_password:
                st.error("The passwords do not match.")
            else:
                try:
                    change_password(str(user["id"]), new_password, forced=True)
                except (OSError, sqlite3.Error, TypeError, ValueError) as exc:
                    st.error(str(exc))
                else:
                    revoke_session(token)
                    st.session_state.clear()
                    st.success("Password changed. Sign in again with your new password.")
                    st.rerun()
        st.caption(f"Use at least {password_min_length()} characters.")
        st.stop()
    return user


CURRENT_USER = _render_authentication_gate()


def search_jobs(**kwargs):
    """Search only rows this signed-in user may access (admins can see all)."""
    return _catalog_search_jobs(
        access_owner_id=str(CURRENT_USER["id"]),
        accessible_project_ids=accessible_project_ids(CURRENT_USER),
        access_admin=CURRENT_USER.get("role") == "admin",
        **kwargs,
    )


def job_list_facets() -> dict[str, list[str]]:
    return _catalog_job_list_facets(
        access_owner_id=str(CURRENT_USER["id"]),
        accessible_project_ids=accessible_project_ids(CURRENT_USER),
        access_admin=CURRENT_USER.get("role") == "admin",
    )


def _require_current_user_job_access(job: dict) -> None:
    record = get_job_access_record(str(job.get("task_group") or ""), str(job.get("run_id") or ""))
    if not can_access_job(CURRENT_USER, record):
        raise PermissionError("Job not found.")


def _inspect_current_user_job(reference: str) -> dict:
    try:
        row = resolve_job(reference)
    except AutomationError as exc:
        if CURRENT_USER.get("role") != "admin":
            raise AutomationError("No accessible job matches that reference.") from exc
        raise
    _require_current_user_job_access(row)
    return inspect_job(reference)


def _cancel_current_user_job(reference: str, *, reason: str) -> dict:
    try:
        row = resolve_job(reference)
    except AutomationError as exc:
        if CURRENT_USER.get("role") != "admin":
            raise AutomationError("No accessible job matches that reference.") from exc
        raise
    _require_current_user_job_access(row)
    return cancel_job(reference, reason=reason)

if "page" not in st.session_state:
    st.session_state["page"] = "All tools"

# Same-tab View links in My results pass the job reference through the query
# string. Consume it once so the detail view is selected without opening a tab.
requested_job = st.query_params.get("job_ref")
if requested_job:
    st.session_state["selected-job-reference"] = str(requested_job)
    st.session_state["page"] = "Job details"
    del st.query_params["job_ref"]

with st.sidebar:
    st.markdown(
        """
        <style>
        [data-testid="stSidebar"] {
            background: #f7f8f7;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] button {
            min-height: 2.65rem;
            justify-content: flex-start;
            border: 0;
            border-radius: 0.65rem;
            padding: 0.4rem 0.65rem;
            text-align: left;
            box-shadow: none;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] button p {
            text-align: left;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] [data-testid="stBaseButton-tertiary"] {
            background: transparent;
            color: #46514b;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] [data-testid="stBaseButton-tertiary"]:hover {
            background: #edf1ee;
            color: #1f5d41;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] [data-testid="stBaseButton-primary"] {
            background: #e5ece7;
            color: #1c6342;
            font-weight: 600;
        }
        [data-testid="stSidebar"] [class*="st-key-nav-"] [data-testid="stBaseButton-primary"]:hover {
            background: #dce7df;
            color: #174f35;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.markdown("### 🧬 Molecular Design Studio")
    st.caption(f"{CURRENT_USER['display_name']} · {str(CURRENT_USER['role']).title()}")
    if is_development_mode() and not cluster_config.is_worker():
        st.warning("Development mode · sign-in disabled")
    elif st.button("Sign out", key="account-sign-out", icon=":material/logout:", use_container_width=True):
        revoke_session(st.session_state.get("studio-session-token"))
        st.session_state.clear()
        st.rerun()
    st.divider()

    if st.button(
        "All tools",
        key="nav-all-tools",
        icon=":material/apps:",
        type="primary" if st.session_state["page"] in {"All tools", "Workflow setup"} else "tertiary",
        use_container_width=True,
    ):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.session_state.pop("selected_workflow_group", None)
        st.rerun()

    st.divider()
    for label, icon, page in (
        ("My results", ":material/query_stats:", "My runs"),
        ("Database", ":material/storage:", "Database"),
        ("Pipelines", ":material/account_tree:", "Pipelines"),
        ("Projects", ":material/folder_open:", "Projects"),
        *(([("Users & access", ":material/manage_accounts:", "Users & access")] if CURRENT_USER.get("role") == "admin" and not is_development_mode() else [])),
        ("Account", ":material/person:", "Account"),
    ):
        if st.button(
            label,
            key=f"nav-{label.lower().replace(' ', '-')}",
            icon=icon,
            type="primary" if st.session_state["page"] == page or (page == "My runs" and st.session_state["page"] == "Job details") else "tertiary",
            use_container_width=True,
        ):
            st.session_state["page"] = page
            st.rerun()

    if CURRENT_USER.get("role") == "admin":
        st.divider()
        if st.button(
            "Settings",
            key="nav-settings",
            icon=":material/settings:",
            type="primary" if st.session_state["page"] == "Settings" else "tertiary",
            use_container_width=True,
        ):
            st.session_state["page"] = "Settings"
            st.rerun()

if cluster_config.enabled():
    from mn_molecular_design_studio.cluster.ui import execution_selector
    execution_selector(CURRENT_USER, str(st.session_state.get("studio-session-token") or "development-mode"))


def _workflow_discovery_tags(row: dict) -> dict[str, list[str]]:
    workflow = str(row.get("workflow", ""))
    prefix = workflow.split(".", 1)[0]
    function_by_prefix = {
        "target": "Target Preparation",
        "detection": "Binding-site Detection",
        "design": "De Novo Binder Design",
        "sequence": "Protein Design",
        "refolding": "Structure Prediction",
        "analysis": "Structure Scoring",
        "import": "Molecule Import",
        "benchmark": "Benchmarking",
        "campaign": "Multi-step Validation",
    }
    function_by_workflow = {
        "target.prepare": "Structure Preparation",
        "target.crop": "Structure Cropping",
        "target.mask": "Target Masking",
        "target.refolding_candidates": "Candidate Preparation",
        "target.refolding_evaluation": "Structure Prediction",
        "design.esmfold2_screening": "Binder Screening",
        "refolding.esmfold_monomer": "Cofolding & Structure Prediction",
        "refolding.esmfold2_cofolding": "Cofolding & Structure Prediction",
        "refolding.af2_initial_guess": "Cofolding & Structure Prediction",
        "refolding.af3_jax_cofolding": "Cofolding & Structure Prediction",
        "refolding.optimized_prediction": "Cofolding & Structure Prediction",
        "sequence.ligandmpnn": "Protein Sequence Design",
        "sequence.mpnn": "Protein Sequence Design",
        "sequence.esm_if1": "Inverse Folding",
        "sequence.caliby": "Inverse Folding",
        "sequence.boltzgen_inverse": "Inverse Folding",
        "sequence.foundry_mpnn": "Protein Sequence Design",
        "sequence.pipeline": "Protein Sequence Design",
        "analysis.ipsae": "Interface Confidence Scoring",
        "analysis.rosetta": "Rosetta Structure Scoring",
    }
    modalities_by_prefix = {
        "target": ["Protein"],
        "detection": ["Protein"],
        "design": ["Protein"],
        "sequence": ["Protein", "Antibody", "Peptide"],
        "refolding": ["Protein", "Antibody", "Peptide"],
        "analysis": ["Protein", "Antibody", "Peptide"],
        "import": ["Protein", "Antibody", "Peptide"],
        "benchmark": ["Protein", "Antibody", "Peptide"],
        "campaign": ["Protein", "Antibody", "Peptide"],
    }
    input_types_by_prefix = {
        "target": ["Structure (PDB)"],
        "detection": ["Structure (PDB)"],
        "design": ["Structure (PDB)"],
        "sequence": ["Existing Results"],
        "refolding": ["Existing Results"],
        "analysis": ["Existing Results"],
        "import": ["Structure (PDB)", "Sequence"],
        "benchmark": ["Existing Results"],
        "campaign": ["Existing Results"],
    }
    input_types_by_workflow = {
        "sequence.mpnn": ["Structure (PDB)", "Structure (mmCIF)", "Existing Results"],
        "sequence.esm_if1": ["Structure (PDB)", "Structure (mmCIF)", "Existing Results"],
        "sequence.caliby": ["Structure (PDB)", "Structure (mmCIF)", "Existing Results"],
        "sequence.boltzgen_inverse": ["Structure (PDB)", "Structure (mmCIF)", "Existing Results"],
        **{workflow:["Structure (PDB)","Structure (mmCIF)"] for workflow in generator_input.WORKFLOWS.values()},
        "design.boltzgen": ["Structure (PDB)", "Structure (mmCIF)"],
        "design.boltzgen_generator": ["Structure (PDB)", "Structure (mmCIF)"],
        "design.proteina_complexa": ["Structure (PDB)", "Structure (mmCIF)"],
        "design.proteina_complexa_generate": ["Structure (PDB)", "Structure (mmCIF)"],
        "design.rfdiffusion3": ["Structure (PDB)", "Structure (mmCIF)", "JSON"],
        "design.rfpeptides": ["Structure (PDB)", "Structure (mmCIF)"],
        "refolding.esmfold2_cofolding": ["Protein Sequences"],
        "refolding.af3_jax_cofolding": ["Molecule Fields"],
        "refolding.optimized_prediction": ["Molecule Fields"],
        "analysis.ipsae": ["Existing Results"],
        "analysis.rosetta": ["Existing Results"],
        "import.bindcraft": ["Result Folder"],
        "import.table": ["CSV / Table", "Structure (PDB)", "Sequence"],
        "benchmark.dataset": ["CSV / Table", "Result Folder", "Structure (PDB)", "Existing Results"],
        "benchmark.precomputed_metrics": ["CSV / Table"],
        "benchmark.esmfold2": ["CSV / Table", "Sequence"],
        "benchmark.refold_candidates": ["Existing Results"],
        "benchmark.refolding_capacity": ["Existing Results", "Structure (PDB)"],
        "benchmark.design_capacity": ["Structure (PDB)"],
        "benchmark.practical_capacity": ["Existing Results"],
        "benchmark.collection": ["Existing Results"],
        "benchmark.matrix_workspace": ["Existing Results"],
        "campaign.validation_sequence": ["Existing Results"],
    }
    functions = [function_by_workflow.get(workflow, function_by_prefix.get(prefix, "Molecular Design"))]
    if workflow in WORKFLOW_ROLES:
        functions.insert(0, WORKFLOW_ROLES[workflow])
    if workflow in {"design.bindcraft2", "design.bindcraft2_generator"}:
        functions.append("Antibody Design")
    return {
        "functions": functions,
        "modalities": (
            ["Protein", "Antibody", "Peptide"]
            if workflow in {"design.bindcraft2", "design.bindcraft2_generator"}
            else ["Protein", "DNA", "RNA", "Ligand"]
            if workflow in {"refolding.af3_jax_cofolding", "design.rfdiffusion3"}
            else ["Protein", "Peptide"] if workflow == "design.rfpeptides"
            else ["Protein", "Ligand"] if workflow in {"design.proteina_complexa", "design.proteina_complexa_generate"}
            else ["Protein"]
            if workflow == "refolding.esmfold2_cofolding"
            else modalities_by_prefix.get(prefix, ["Protein"])
        ),
        "input_types": input_types_by_workflow.get(
            workflow,
            input_types_by_prefix.get(prefix, ["Structure (PDB)"]),
        ),
    }


def _render_tool_catalog() -> None:
    st.title("All tools")
    st.caption("Explore tools by modality, function, or input type.")
    st.markdown(
        """
        <style>
        [class*="st-key-workflow-card-"] button {
            box-sizing: border-box;
            min-height: 236px !important;
            height: auto !important;
            padding: 1.05rem 1.15rem;
            border: 1px solid #e3e9e5;
            border-radius: 14px;
            background: #ffffff;
            box-shadow: 0 2px 8px rgba(26, 52, 39, 0.06);
            text-align: left;
            align-items: flex-start;
            justify-content: flex-start;
            white-space: pre-wrap !important;
            transition: border-color 120ms ease, box-shadow 120ms ease;
        }
        [class*="st-key-workflow-card-"] button:hover {
            border-color: #79a58b;
            box-shadow: 0 4px 12px rgba(26, 52, 39, 0.11);
        }
        [class*="st-key-workflow-card-"] button * {
            width: 100%;
            text-align: left;
            white-space: pre-wrap !important;
        }
        [class*="st-key-workflow-card-"] button p {
            margin: 0 0 0.35rem;
            line-height: 1.4;
        }
        [class*="st-key-workflow-card-"] button [data-testid="stMarkdownContainer"] p:first-of-type {
            color: #68766e;
            font-size: 0.76rem;
            font-weight: 600;
            letter-spacing: 0.035em;
            line-height: 1.2;
            width: fit-content;
            padding: 0.18rem 0.42rem;
            border-radius: 0.3rem;
            background: #edf3ef;
        }
        [class*="st-key-workflow-card-"] button [data-testid="stMarkdownContainer"] p:nth-of-type(2) {
            color: #1d3026;
            font-size: 1.2rem;
            font-weight: 700;
            line-height: 1.25;
        }
        [class*="st-key-workflow-card-"] button [data-testid="stMarkdownContainer"] p:nth-of-type(n+3) {
            font-size: 0.9rem;
        }
        [class*="st-key-workflow-card-engine-"] button [data-testid="stMarkdownContainer"] p:nth-of-type(3) {
            color: #77837c;
            font-size: 0.74rem;
            line-height: 1.2;
            margin-top: -0.16rem;
        }
        [class*="st-key-workflow-card-engine-"] button [data-testid="stMarkdownContainer"] p:nth-of-type(3) em {
            font-size: 0.62rem;
            font-style: normal;
            font-weight: 600;
        }
        [class*="st-key-workflow-card-"]:not([class*="st-key-workflow-card-engine-"]) button [data-testid="stMarkdownContainer"] p:nth-of-type(4),
        [class*="st-key-workflow-card-engine-"] button [data-testid="stMarkdownContainer"] p:nth-of-type(5) {
            color: #68766e;
            font-size: 0.68rem;
            font-weight: 600;
            letter-spacing: 0.02em;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

    all_workflows = [
        row for row in workflow_catalog()
        if row["workflow"] in TOOL_CATALOG_WORKFLOW_IDS
    ]
    workflow_by_id = {row["workflow"]: row for row in all_workflows}
    tags_by_workflow = {row["workflow"]: _workflow_discovery_tags(row) for row in all_workflows}
    grouped_workflow_ids: set[str] = set()
    tool_rows: list[dict] = []
    for group_id, group in TOOL_CATALOG_GROUPS.items():
        modes = {
            label: workflow_id
            for label, workflow_id in group["modes"].items()
            if workflow_id in workflow_by_id
        }
        if len(modes) < 2:
            continue
        member_rows = [workflow_by_id[workflow_id] for workflow_id in modes.values()]
        group_tags = {
            key: list(dict.fromkeys(tag for row in member_rows for tag in tags_by_workflow[row["workflow"]][key]))
            for key in ("functions", "modalities", "input_types")
        }
        tool_rows.append(
            {
                "key": f"group-{group_id}",
                "group_id": group_id,
                "title": group["title"],
                "description": group["description"],
                "workflow_ids": list(modes.values()),
                "workflow_rows": member_rows,
                "tags": group_tags,
                "mode_count": len(modes),
                "engine_version": group.get("engine_version"),
                "engine_optimized": group.get("engine_optimized", False),
            }
        )
        grouped_workflow_ids.update(modes.values())

    for row in all_workflows:
        workflow_id = row["workflow"]
        if workflow_id == "refolding.optimized_prediction":
            continue
        if workflow_id in grouped_workflow_ids:
            continue
        base_title = workflow_id.split(".", 1)[-1].replace("_", " ").title()
        tool_rows.append(
            {
                "key": workflow_id,
                "group_id": None,
                "workflow": workflow_id,
                "title": WORKFLOW_TITLES.get(base_title, base_title),
                "description": row["description"],
                "workflow_ids": [workflow_id],
                "workflow_rows": [row],
                "tags": tags_by_workflow[workflow_id],
                "engine_version": (
                    TOOL_WORKFLOW_CODE_VERSIONS.get(workflow_id, (None, False))[0]
                ),
                "engine_optimized": (
                    TOOL_WORKFLOW_CODE_VERSIONS.get(workflow_id, (None, False))[1]
                ),
            }
        )

    optimized_row = workflow_by_id.get("refolding.optimized_prediction")
    if optimized_row:
        for engine_id in OPTIMIZED_IMAGE_TILE_IDS:
            if engine_id == "boltzgen_fold":
                continue
            profile = optimized_prediction.ENGINE_PROFILES[engine_id]
            supported = set(profile.get("molecule_inputs", ()))
            tool_rows.append({
                "key": f"optimized-image-{engine_id}",
                "group_id": None,
                "workflow": "refolding.optimized_prediction",
                "optimized_engine": engine_id,
                "title": profile["label"],
                "description": profile.get("input_summary", "Predict protein structures."),
                "engine_version": profile.get("code_version"),
                "engine_optimized": profile.get("optimized", False),
                "workflow_ids": ["refolding.optimized_prediction"],
                "workflow_rows": [optimized_row],
                "tags": {
                    "functions": ["Cofolding & Structure Prediction"],
                    "modalities": ["Protein", *(["DNA"] if "dna_sequences" in supported else []),
                                   *(["RNA"] if "rna_sequences" in supported else []),
                                   *(["Ligand"] if supported.intersection({"ligand_ccd_codes", "ligand_smiles"}) else [])],
                    "input_types": ["Molecule Fields" if supported else "Protein Sequences"],
                },
            })

    tool_rows.sort(key=lambda row: row["title"].casefold())
    modalities = sorted({tag for row in tool_rows for tag in row["tags"]["modalities"]})
    functions = sorted({tag for row in tool_rows for tag in row["tags"]["functions"]})
    input_types = sorted({tag for row in tool_rows for tag in row["tags"]["input_types"]})
    search_col, modality_col, function_col, input_col = st.columns([2.2, 1, 1, 1])
    query = search_col.text_input("Search tools", placeholder="Try BindCraft2, RFdiffusion, or structure prediction")
    with modality_col.popover("Modality", use_container_width=True):
        selected_modalities = st.multiselect("Modality", modalities, key="tools-modality-filter", label_visibility="collapsed")
    with function_col.popover("Function", use_container_width=True):
        selected_functions = st.multiselect("Function", functions, key="tools-function-filter", label_visibility="collapsed")
    with input_col.popover("Input type", use_container_width=True):
        selected_input_types = st.multiselect("Input type", input_types, key="tools-input-type-filter", label_visibility="collapsed")

    rows = tool_rows
    if selected_modalities:
        rows = [row for row in rows if set(row["tags"]["modalities"]) & set(selected_modalities)]
    if selected_functions:
        rows = [row for row in rows if set(row["tags"]["functions"]) & set(selected_functions)]
    if selected_input_types:
        rows = [row for row in rows if set(row["tags"]["input_types"]) & set(selected_input_types)]
    if query.strip():
        normalized_query = query.strip().lower()
        rows = [
            row for row in rows
            if normalized_query in " ".join(
                [
                    row["key"],
                    row["title"],
                    row["description"],
                    *row["workflow_ids"],
                    *(member["description"] for member in row["workflow_rows"]),
                    *row["tags"]["functions"],
                    *row["tags"]["modalities"],
                    *row["tags"]["input_types"],
                ]
            ).lower()
        ]

    st.caption(f"{len(rows)} tools")
    column_count = 3
    for row_start in range(0, len(rows), column_count):
        row_columns = st.columns(column_count, gap="medium")
        for column, row in zip(row_columns, rows[row_start : row_start + column_count]):
            title = row["title"]
            description = row["description"]
            if len(description) > 118:
                description = description[:115].rsplit(" ", 1)[0] + "…"
            tags = row["tags"]
            function = tags["functions"][0].upper()
            modality_tags = " · ".join(tags["modalities"]).upper()
            input_tags = " · ".join(tags["input_types"]).upper()
            engine_version = row.get("engine_version")
            if engine_version and row.get("engine_optimized"):
                engine_version += " *(opt)*"
            version_line = f"**{engine_version}**\n\n" if engine_version else ""
            if row.get("optimized_engine"):
                setup_note = f"{len(optimized_prediction.ENGINE_PROFILES[row['optimized_engine']]['modes'])} validated modes · Set up →"
            elif row.get("workflow") in {"design.proteina_complexa", "design.proteina_complexa_generate"}:
                setup_note = "3 design tasks · Set up →" if row["workflow"].endswith("_generate") else "3 pipeline tasks · Set up →"
            elif row.get('workflow') in generator_input.WORKFLOWS.values():
                engine=next(k for k,v in generator_input.WORKFLOWS.items() if v==row['workflow'])
                count=len(generator_input.ENGINES[engine]['tasks'])
                setup_note=str(count)+' generation task'+('s' if count!=1 else '')+' · Set up →'
            elif row["group_id"]:
                setup_note = f"{row['mode_count']} prediction modes · Set up →"
            else:
                workflow_id = row["workflow"]
                hidden_analysis_inputs = {"source_run_dir", "candidates_jsonl", "tool"}
                hidden_tool_inputs = (
                    {"gpu_device"}
                    if workflow_id in {"refolding.esmfold2_cofolding", "refolding.af3_jax_cofolding"}
                    else set()
                )
                input_count = sum(
                    1
                    for parameter in row["workflow_rows"][0]["parameters"]
                    if parameter["name"] not in hidden_analysis_inputs | hidden_tool_inputs
                )
                setup_note = f"{input_count} configurable inputs · Set up →"
            label = (
                f"**{function}**\n\n"
                f"**{title}**\n\n"
                f"{version_line}"
                f"{description}\n\n"
                f"{modality_tags}  ·  INPUT: {input_tags}\n\n"
                f"{setup_note}"
            )
            with column:
                card_prefix = "workflow-card-engine-" if engine_version else "workflow-card-"
                card_key = card_prefix + row["key"].replace(".", "-").replace("_", "-").lower()
                with st.container(key=card_key):
                    if st.button(label, key=card_key + "-button", use_container_width=True):
                        if row.get("optimized_engine"):
                            st.session_state.pop("selected_workflow_group", None)
                            st.session_state["selected_workflow"] = row["workflow"]
                            st.session_state["selected_workflow_engine"] = row["optimized_engine"]
                        elif row["group_id"]:
                            st.session_state["selected_workflow_group"] = row["group_id"]
                            st.session_state.pop("selected_workflow", None)
                            st.session_state.pop("selected_workflow_engine", None)
                        else:
                            st.session_state.pop("selected_workflow_group", None)
                            st.session_state["selected_workflow"] = row["workflow"]
                            st.session_state.pop("selected_workflow_engine", None)
                        st.session_state["page"] = "Workflow setup"
                        st.rerun()


def _parameter_label(parameter: dict) -> str:
    return parameter["name"].replace("_", " ").strip().capitalize()


def _parameter_widget(parameter: dict, *, workflow_id: str, group: str):
    name = parameter["name"]
    annotation = parameter["type"].lower()
    default = parameter.get("default")
    required = parameter.get("required", False)
    label = _parameter_label(parameter) + (" *" if required else "")
    help_text = f"Input name: `{name}` · Type: `{parameter['type']}`"
    widget_key = f"setup-{workflow_id}-{group}-{name}"

    choices = parameter.get("choices")
    if isinstance(choices, list) and choices:
        index = choices.index(default) if default in choices else 0
        return st.selectbox(label, choices, index=index, help=help_text, key=widget_key)

    if "bool" in annotation:
        return st.checkbox(label, value=bool(default) if default is not None else False, help=help_text, key=widget_key)

    if any(type_name in annotation for type_name in ("list[", "tuple[", "set[")):
        if isinstance(default, (list, tuple, set)):
            initial = ", ".join(str(item) for item in default)
        else:
            initial = "" if default is None else str(default)
        return st.text_input(
            label,
            value=initial,
            placeholder="Separate values with commas",
            help=help_text,
            key=widget_key,
        )

    if "dict[" in annotation or annotation == "object":
        initial = "" if default is None else json.dumps(default, indent=2)
        return st.text_area(label, value=initial, height=110, help=help_text, key=widget_key)

    numeric = annotation in {"int", "float"} or annotation in {"int | none", "float | none"}
    if numeric and default is not None:
        value = float(default) if "float" in annotation else int(default)
        return st.number_input(
            label,
            value=value,
            step=0.1 if "float" in annotation else 1,
            help=help_text,
            key=widget_key,
        )

    if isinstance(default, str) and ("sequence" in name or "text" in name or "json" in name):
        return st.text_area(label, value=default, height=110, help=help_text, key=widget_key)

    return st.text_input(
        label,
        value="" if default is None else str(default),
        help=help_text,
        key=widget_key,
    )


def _parameter_value(raw_value, parameter: dict):
    annotation = parameter["type"].lower()
    default = parameter.get("default")
    required = parameter.get("required", False)

    if "bool" in annotation:
        return bool(raw_value)
    if isinstance(raw_value, (int, float)):
        return raw_value

    text = str(raw_value).strip()
    if text == "":
        if "none" in annotation:
            return None
        if required:
            return ""
        return default

    if any(type_name in annotation for type_name in ("list[", "tuple[", "set[")):
        values = [part.strip() for part in text.replace("\n", ",").split(",") if part.strip()]
        return values
    if "dict[" in annotation or annotation == "object":
        return json.loads(text)
    if annotation in {"int", "int | none"}:
        return int(text)
    if annotation in {"float", "float | none"}:
        return float(text)
    return text


def _parse_amino_acid_bias(text: str) -> dict[str, float]:
    bias: dict[str, float] = {}
    for part in re.split(r"[,;\n]+", text.strip()):
        if not part.strip():
            continue
        amino_acid, separator, raw_weight = part.partition(":")
        if not separator or len(amino_acid.strip()) != 1:
            raise ValueError("Write amino-acid biases as one-letter codes and weights, e.g. W:0.4, F:0.8.")
        code = amino_acid.strip().upper()
        try:
            weight = float(raw_weight.strip())
        except ValueError as exc:
            raise ValueError(f"Bias for {code} must be a number, for example {code}:0.8.") from exc
        if code in bias:
            raise ValueError(f"Amino-acid bias for {code} is listed more than once.")
        bias[code] = weight
    return bias


def _studio_projects() -> list[dict]:
    organization = load_organization(STUDIO_WORKSPACE)
    projects = organization.get("projects")
    if not isinstance(projects, list):
        return []
    visible_ids = set(accessible_project_ids(CURRENT_USER))
    return sorted(
        [
            project for project in projects
            if isinstance(project, dict)
            and project.get("id")
            and project.get("name")
            and (CURRENT_USER.get("role") == "admin" or str(project.get("id")) in visible_ids)
        ],
        key=lambda project: str(project.get("name", "")).casefold(),
    )


def _render_project_selector(
    key_prefix: str,
    *,
    label: str = "Project",
    allow_create: bool = True,
) -> str | None:
    projects = _studio_projects()
    project_ids = ["", *[str(project["id"]) for project in projects]]
    key = f"{key_prefix}-project-id"
    pending = st.session_state.get("project-pending-selection", {})
    if isinstance(pending, dict) and key in pending:
        st.session_state[key] = pending.pop(key)
        if not pending:
            st.session_state.pop("project-pending-selection", None)
    if st.session_state.get(key) not in project_ids:
        st.session_state[key] = ""

    if allow_create:
        project_col, add_col = st.columns([0.76, 0.24], vertical_alignment="bottom")
    else:
        project_col = st.container()
        add_col = None
    selected = project_col.selectbox(
        f"{label} *",
        project_ids,
        format_func=lambda value: "Select a project" if value == "" else next(
            (str(project["name"]) for project in projects if project["id"] == value), "Unknown project"
        ),
        key=key,
        help="Every job must belong to a project. Select an existing project or create one.",
    )
    if add_col is not None and add_col.button("New", icon=":material/add:", key=f"{key_prefix}-create-project", help="Create a project"):
        _render_project_dialog(key)
    if not projects:
        project_col.caption("Create a project before submitting or importing a job.")
    return selected or None


@st.dialog("Manage projects", width="medium")
def _render_project_dialog(select_key: str | None = None) -> None:
    st.markdown(
        """
        <style>
        div[role="dialog"] [data-testid="stFormSubmitButton"] button,
        div[role="dialog"] form button[type="submit"],
        [data-testid="stFormSubmitButton"] button {
            background: #2f6b50 !important;
            border-color: #2f6b50 !important;
            border-radius: 999px !important;
            color: #ffffff !important;
            font-weight: 600 !important;
            padding-left: 1.1rem !important;
            padding-right: 1.1rem !important;
        }
        div[role="dialog"] [data-testid="stFormSubmitButton"] button:hover,
        div[role="dialog"] form button[type="submit"]:hover,
        [data-testid="stFormSubmitButton"] button:hover {
            background: #24553e !important;
            border-color: #24553e !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.caption("Create a project to keep jobs from different tools together.")
    with st.form("studio-project-create-form"):
        name = st.text_input("Project name", key="studio-project-new-name").strip()
        description = st.text_area("Description (optional)", key="studio-project-new-description", height=90)
        color = st.color_picker("Color", value="#ef5350", key="studio-project-new-color")
        # Form inputs do not rerun the script as they change, so deriving the
        # button's disabled state from ``name`` leaves it disabled forever.
        submitted = st.form_submit_button("Create project", type="primary")
    if submitted:
        try:
            project = create_project(
                name,
                description,
                color,
                STUDIO_WORKSPACE,
                owner_id=str(CURRENT_USER["id"]),
            )
            if CURRENT_USER.get("role") != "admin":
                grant_project_access(
                    project["id"],
                    str(CURRENT_USER["id"]),
                    actor=CURRENT_USER,
                    allow_self_grant=True,
                )
        except (OSError, TypeError, ValueError) as exc:
            st.error(f"Could not create project: {exc}")
        else:
            if select_key:
                pending = st.session_state.get("project-pending-selection", {})
                pending = pending if isinstance(pending, dict) else {}
                pending[select_key] = project["id"]
                st.session_state["project-pending-selection"] = pending
            else:
                st.session_state["projects-last-created"] = project["name"]
            st.rerun()

    projects = _studio_projects()
    st.divider()
    st.markdown("**Existing projects**")
    if projects:
        for project in projects:
            color = str(project.get("color") or "#75817b")
            if not re.fullmatch(r"#[0-9A-Fa-f]{6}", color):
                color = "#75817b"
            description = str(project.get("description") or "")
            description_html = f'<br><span style="color:#75817b;font-size:.9em">{escape(description)}</span>' if description else ""
            st.markdown(
                f'<div style="padding:.35rem 0"><span style="color:{color}">●</span> '
                f'<strong>{escape(str(project["name"]))}</strong>{description_html}</div>',
                unsafe_allow_html=True,
            )
    else:
        st.caption("No projects yet.")


def _render_bindcraft2_setup(schema: dict, *, generator_only: bool = False) -> None:
    if st.button("← Back to All tools", key=f"back-from-bindcraft2-{int(generator_only)}"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()

    with st.container():
        st.markdown(
            """
            <style>
            .block-container { padding-top: 1.35rem; max-width: 1800px; }
            .bc2-tool-header { text-align: center; padding: .1rem 0 1rem; }
            .bc2-tool-header h1 { margin: 0 0 .3rem; font-size: 2rem; }
            .bc2-tool-header p { color: #66736d; margin: .4rem 0; }
            .bc2-tool-meta { color: #738077; font-size: .84rem; }
            div[data-testid="stButton"] button { border-radius: 999px; }
            div[data-testid="stButton"] button[kind="primary"] {
                background: #2f6b50 !important;
                border-color: #2f6b50 !important;
                color: #ffffff !important;
            }
            div[data-testid="stButton"] button[kind="primary"]:hover {
                background: #24553e !important;
                border-color: #24553e !important;
            }
            [role="radio"][aria-checked="true"] {
                background: #2f6b50 !important;
                border-color: #2f6b50 !important;
                color: #ffffff !important;
            }
            [data-baseweb="tag"] { background: #e5f0e9 !important; color: #28583f !important; }
            </style>
            <div class="bc2-tool-header">
              <h1>__TITLE__</h1>
              <div class="bc2-tool-meta">__META__</div>
              <p>Design binders for your target as miniproteins, antibodies, peptides or scaffolds.</p>
            </div>
            """.replace("__TITLE__", "BindCraft2 Generator" if generator_only else "BindCraft2")
            .replace("__META__", "Gradient generation only · ProteinMPNN sequence redesign is a separate step" if generator_only else "De novo binder design · Long-running workflow"),
            unsafe_allow_html=True,
        )
        with st.expander("About this workflow", expanded=False):
            st.write(schema["description"])
            st.markdown("[BindCraft2 source code](https://github.com/PacesaLab/BindCraft2)")

        job_name_col, project_col = st.columns([1.25, 1.25])
        with job_name_col:
            job_name = st.text_input(
                "Job name",
                value="bindcraft2",
                placeholder="For example, 7JZU binder design",
                key=f"bindcraft2-job-name-{int(generator_only)}",
            ).strip()
        with project_col:
            project_id = _render_project_selector(f"bindcraft2-{int(generator_only)}")

        st.subheader("Modality")
        modality = st.session_state.get("bindcraft2-selected-modality", "binder")
        for row_values in (list(BINDCRAFT2_MODALITY_LABELS)[:5], list(BINDCRAFT2_MODALITY_LABELS)[5:]):
            modality_columns = st.columns(5, gap="small")
            for column, value in zip(modality_columns, row_values):
                with column:
                    if st.button(
                        BINDCRAFT2_SHORT_LABELS[value],
                        key=f"bindcraft2-modality-{value}",
                        type="primary" if modality == value else "secondary",
                        use_container_width=True,
                    ):
                        st.session_state["bindcraft2-selected-modality"] = value
                        st.rerun()
        st.caption(BINDCRAFT2_MODALITY_HELP[modality])
        conformational_objective_key = "bindcraft2-conformational-objective"
        objective_options = ["none", "induced_fit", "fold_switch"]
        try:
            requested_sequence_crops = int(st.session_state.get("bindcraft2-idr-crop-count", 1) or 1)
        except (TypeError, ValueError):
            requested_sequence_crops = 1
        objective_target_type = st.session_state.get("bindcraft2-target-type", "structure")
        second_target_requested = objective_target_type == "structure" and bool(
            st.session_state.get("bindcraft2-use-second-target", False)
        )
        available_objective_options = ["none"]
        one_target_state = objective_target_type == "structure" or (
            objective_target_type == "sequence" and requested_sequence_crops == 1
        )
        if one_target_state and not second_target_requested:
            available_objective_options.append("induced_fit")
            if modality not in bindcraft2.BC2_SCAFFOLD_FORMATS:
                available_objective_options.append("fold_switch")
        if st.session_state.get(conformational_objective_key, "none") not in available_objective_options:
            st.session_state[conformational_objective_key] = "none"
        conformational_objective = st.session_state.get(conformational_objective_key, "none")
        st.markdown("**Optional conformational objective**")
        objective_columns = st.columns(3, gap="small")
        for column, objective in zip(objective_columns, objective_options):
            with column:
                if st.button(
                    BINDCRAFT2_CONFORMATIONAL_OBJECTIVE_LABELS[objective],
                    key=f"bindcraft2-objective-{objective}",
                    type="primary" if conformational_objective == objective else "secondary",
                    use_container_width=True,
                    disabled=(objective not in available_objective_options),
                ):
                    st.session_state[conformational_objective_key] = objective
                    st.rerun()
        st.caption("This adds a conformational objective to the selected binder format; it does not select a separate binder format.")
        st.caption(BINDCRAFT2_CONFORMATIONAL_OBJECTIVE_HELP[conformational_objective])
        if objective_target_type == "sequence" and requested_sequence_crops != 1:
            st.caption("Conformational objectives require one target state. Set sequence crops to one to use them.")
        elif second_target_requested:
            st.caption("Conformational objectives require one target. Turn off the second target to use one.")
        elif modality in bindcraft2.BC2_SCAFFOLD_FORMATS:
            st.caption("Fold Switch is unavailable for scaffold formats because the framework is held fixed.")

        current_target_type = st.session_state.get("bindcraft2-target-type", "structure")
        st.subheader(
            "Target structure" if current_target_type == "structure" else "Disordered target sequence"
        )
        target_type = st.selectbox(
            "Target type",
            options=["structure", "sequence"],
            format_func=lambda value: "Structure (PDB or mmCIF)" if value == "structure" else "Sequence (disordered / IDR)",
            key="bindcraft2-target-type",
            help="Use a folded target structure or an intrinsically disordered target sequence. The selected modality above controls the binder format; Binder is the general-purpose choice for IDR targets, while Peptide means a short linear binder.",
        )
        second_target_requested = target_type == "structure" and bool(
            st.session_state.get("bindcraft2-use-second-target", False)
        )
        target_path: Path | None = None
        target_name = ""
        target_sequence = ""
        target_fasta_text = ""
        target_error = ""
        target_crop_lengths: tuple[int, int] | None = None
        target_crop_count: int | None = None
        target_chains: list[str] = []
        hotspots = ""
        coldspots = ""
        design_target_path: Path | None = None
        valid_trim_ranges = True
        target_input = None
        if target_type == "structure":
            structure_format = st.radio(
                "Structure format",
                options=["pdb", "cif"],
                format_func=lambda value: "PDB" if value == "pdb" else "CIF / mmCIF",
                horizontal=True,
                key="bindcraft2-structure-format",
            )
            if structure_format == "pdb":
                st.caption("Fetch by PDB accession or upload a PDB file.")
                target_source = render_target_source(
                    key_prefix="bindcraft2",
                    label="Target",
                    example_code="7JZU",
                    upload_label="Drop a PDB file here or browse",
                    upload_help="Upload a PDB (.pdb/.ent) structure, or fetch it by PDB accession above.",
                    upload_file_types=["pdb", "ent"],
                    upload_expander_label="Upload a PDB file",
                    upload_in_expander=True,
                    owner_id=str(CURRENT_USER["id"]),
                )
            else:
                st.caption("CIF/mmCIF structures are provided by upload; accession fetch is available in PDB mode.")
                target_source = render_target_source(
                    key_prefix="bindcraft2-cif",
                    label="Target",
                    upload_label="Drop a CIF/mmCIF file here or browse",
                    upload_help="Upload a CIF (.cif) or mmCIF (.mmcif) structure file.",
                    allow_cif=True,
                    fetch_enabled=False,
                    upload_file_types=["cif", "mmcif"],
                    upload_in_expander=False,
                    owner_id=str(CURRENT_USER["id"]),
                )
            target_input = render_target_editor(
                target_source,
                key_prefix="bindcraft2",
                label="Target",
                prefer_chain_b=True,
                show_structure_tools=False,
            )
            target_path = target_input.source_path
            target_name = target_source.name
            target_chains = target_input.chains
            hotspots = target_input.hotspots
            design_target_path = target_input.design_path
            valid_trim_ranges = target_input.valid_trim_ranges

            st.markdown("#### Coldspot selection")
            use_coldspots = st.checkbox(
                "Enable coldspot selection",
                value=False,
                key="bindcraft2-enable-coldspots",
                help="Choose target residues the binder should avoid. Disabled by default.",
                disabled=target_input.source_path is None or not target_input.chains,
            )
            if use_coldspots:
                coldspots = render_target_coldspots(
                    target_input,
                    key_prefix="bindcraft2",
                    legacy_coldspots_key="bindcraft2-coldspots",
                    show_structure_tools=False,
                )
            else:
                target_input.coldspots = ""
        else:
            st.caption("BindCraft2 treats this as an intrinsically disordered region (IDR) and co-folds it with the binder. Upload one named FASTA record; its header becomes the target name. Only the 20 standard amino acids are accepted. The Binder modality is the general-purpose starting point; Peptide selects a short linear binder format.")
            sequence_file = st.file_uploader(
                "Upload target FASTA",
                type=["fasta", "fa", "faa"],
                key="bindcraft2-sequence-upload",
                help="Upload a FASTA file with one record. The text after > is used as the target name.",
            )
            if sequence_file is not None:
                try:
                    target_fasta_text = sequence_file.getvalue().decode("utf-8-sig")
                    target_name, target_sequence = bindcraft2.parse_disordered_target_fasta(target_fasta_text)
                    st.caption(f"Target name from FASTA: {target_name}")
                except UnicodeDecodeError:
                    target_error = "The FASTA file must be UTF-8 text."
                except ValueError as exc:
                    target_error = str(exc)
            st.markdown("**Disordered-target sampling**")
            crop_count_col, crop_min_col, crop_max_col = st.columns(3)
            with crop_count_col:
                target_crop_count = int(st.number_input(
                    "Sequence crops",
                    min_value=1,
                    max_value=100,
                    value=1,
                    step=1,
                    key="bindcraft2-idr-crop-count",
                    help="Number of target fragments sampled from the disordered sequence. The BindCraft2 default is one crop.",
                ))
            with crop_min_col:
                crop_min = int(st.number_input(
                    "Minimum crop length",
                    min_value=1,
                    max_value=10000,
                    value=10,
                    step=1,
                    key="bindcraft2-idr-crop-min",
                ))
            with crop_max_col:
                crop_max = int(st.number_input(
                    "Maximum crop length",
                    min_value=1,
                    max_value=10000,
                    value=40,
                    step=1,
                    key="bindcraft2-idr-crop-max",
                ))
            target_crop_lengths = (crop_min, crop_max)
            use_coldspots = False

        if bindcraft2.BC2_PRESET_LENGTHS.get(modality):
            default_min, default_max = bindcraft2.BC2_PRESET_LENGTHS[modality]
            st.subheader("Binder length")
            min_col, max_col = st.columns(2)
            with min_col:
                binder_min = st.number_input(
                    "Minimum length",
                    min_value=1,
                    max_value=2000,
                    value=default_min,
                    step=1,
                    key=f"bindcraft2-min-length-{modality}",
                )
            with max_col:
                binder_max = st.number_input(
                    "Maximum length",
                    min_value=1,
                    max_value=2000,
                    value=default_max,
                    step=1,
                    key=f"bindcraft2-max-length-{modality}",
                )
            binder_length = f"{binder_min}-{binder_max}"
        else:
            binder_length = None
            st.caption("The scaffold sets the framework length; variable loops can change the final length.")

        st.subheader("Generation" if generator_only else "Design goal")
        if generator_only:
            number_of_final_designs = 1
            max_trajectories = st.number_input(
                "Number of gradient trajectories",
                min_value=1,
                max_value=1000,
                value=10,
                step=1,
                key="bindcraft2-generator-max-trajectories",
                help="Each trajectory generates a target–binder pose before ProteinMPNN sequence redesign and downstream acceptance filters.",
            )
            st.caption(
                "Only trajectories that pass BindCraft2's internal gradient-stage gates are forwarded: "
                "the normal handoff rule is that the trajectory finishes with no terminated stage. "
                "Default final gates are i_pTM ≥ 0.70 and binder pLDDT ≥ 0.70; stage-specific confidence "
                "gates also apply. Later sequence, refolding, and campaign acceptance filters run after ProteinMPNN."
            )
        else:
            goal_col, attempt_col = st.columns(2)
            with goal_col:
                number_of_final_designs = st.number_input(
                    "Number of designs",
                    min_value=1,
                    max_value=100,
                    value=1,
                    step=1,
                    key="bindcraft2-final-designs",
                    help="Stop when this many accepted designs have been found, or the attempt limit is reached.",
                )
            with attempt_col:
                max_trajectories = st.number_input(
                    "Maximum design attempts",
                    min_value=1,
                    max_value=1000,
                    value=10,
                    step=1,
                    key="bindcraft2-max-trajectories",
                )

        second_target_path: Path | None = None
        second_target_design_path: Path | None = None
        second_target_chains: list[str] = []
        second_target_hotspots = ""
        second_target_coldspots = ""
        second_target_trim_valid = True
        second_target_goal = "bind"
        second_target_strength = 1.0
        redesign_interface = False
        binder_copies = 2
        number_of_domains = 2
        excluded_amino_acids: list[str] = []
        amino_acid_bias: dict[str, float] = {}
        quality_filter_defaults = _bindcraft2_quality_filter_defaults(modality)
        if conformational_objective in bindcraft2.BC2_CONFORMATIONAL_OBJECTIVES:
            quality_filter_defaults["max_binder_rmsd"] = None
        min_interface_iptm = quality_filter_defaults["min_interface_iptm"]
        max_interface_ipae = quality_filter_defaults["max_interface_ipae"]
        min_complex_ptm = quality_filter_defaults["min_complex_ptm"]
        min_unbound_binder_plddt = quality_filter_defaults["min_unbound_binder_plddt"]
        min_interface_residues = quality_filter_defaults["min_interface_residues"]
        max_backbone_clashes = quality_filter_defaults["max_backbone_clashes"]
        max_binder_rmsd = quality_filter_defaults["max_binder_rmsd"]
        min_plddt_screen = quality_filter_defaults["min_plddt_screen"]
        min_plddt_refine = quality_filter_defaults["min_plddt_refine"]
        min_plddt_anneal = quality_filter_defaults["min_plddt_anneal"]
        min_plddt_harden = quality_filter_defaults["min_plddt_harden"]
        min_plddt_mutate = quality_filter_defaults["min_plddt_mutate"]
        min_plddt_final = quality_filter_defaults["min_plddt_final"]
        min_iptm_anneal = quality_filter_defaults["min_iptm_anneal"]
        min_iptm_harden = quality_filter_defaults["min_iptm_harden"]
        min_iptm_mutate = quality_filter_defaults["min_iptm_mutate"]
        max_interdomain_contact_final = quality_filter_defaults["max_interdomain_contact_final"]
        min_domain_separation_ratio_final = quality_filter_defaults["min_domain_separation_ratio_final"]
        min_scaffold_sequence_retained_final = quality_filter_defaults["min_scaffold_sequence_retained_final"]
        max_scaffold_framework_rmsd_final = quality_filter_defaults["max_scaffold_framework_rmsd_final"]
        optional_error = ""
        engine_overrides: dict = {}
        with st.expander("Optional settings", expanded=False):
            st.markdown("**Target preferences**")
            if target_type == "structure":
                use_second_target = st.checkbox(
                    "Add a second target",
                    key="bindcraft2-use-second-target",
                    disabled=(conformational_objective != "none"),
                    help="Add another structure to bind as well or avoid for selectivity.",
                )
                if conformational_objective != "none":
                    st.caption("Turn off the conformational objective before adding a second target.")
            else:
                st.caption("A second target structure is available when the primary target type is Structure.")
                use_second_target = False
            if use_second_target:
                st.caption("Second target structure · PDB format")
                second_source = render_target_source(
                    key_prefix="bindcraft2-second",
                    label="Second target",
                    example_code="1ABC",
                    fetch_button_label="Fetch second target",
                    fetch_button_key="bindcraft2-fetch-second-pdb",
                    upload_label="Upload a second target PDB file",
                    upload_help="Choose a PDB file or fetch it by accession above.",
                    upload_in_expander=False,
                    fetched_path_key="bindcraft2-fetched-second-target",
                    fetched_name_key="bindcraft2-fetched-second-target-name",
                    owner_id=str(CURRENT_USER["id"]),
                )
                second_target = render_target_editor(
                    second_source,
                    key_prefix="bindcraft2-second",
                    label="Second target",
                    legacy_hotspots_key="bindcraft2-second-target-hotspots",
                )
                second_target_path = second_target.source_path
                second_target_design_path = second_target.design_path
                second_target_chains = second_target.chains
                second_target_hotspots = second_target.hotspots
                second_target_trim_valid = second_target.valid_trim_ranges
                if second_target_path is not None:
                    second_target_coldspots = render_target_coldspots(
                        second_target,
                        key_prefix="bindcraft2-second",
                        legacy_coldspots_key="bindcraft2-second-target-coldspots",
                    )
                    second_target_design_path = second_target.design_path
                else:
                    st.info("Upload a second target PDB or fetch one by accession to enable the additional target settings.")

                second_goal_col, second_strength_col = st.columns([1.2, 1])
                with second_goal_col:
                    second_target_goal = st.selectbox(
                        "Second target goal",
                        options=["bind", "avoid"],
                        format_func=lambda value: "Bind it" if value == "bind" else "Avoid it",
                        key="bindcraft2-second-target-goal",
                        help="Bind this target as well, or push the binder away from it. Avoiding one gives you selectivity.",
                    )
                with second_strength_col:
                    second_target_strength = st.number_input(
                        "Second target strength",
                        min_value=0.0,
                        value=1.0,
                        step=0.25,
                        key="bindcraft2-second-target-strength",
                        help="Relative weight for this target in the design objective.",
                    )

            if not generator_only:
                redesign_interface = st.checkbox(
                    "**Redesign Interface**",
                    key="bindcraft2-redesign-interface",
                )
                st.markdown(
                    '<div style="margin:-.65rem 0 .8rem 1.85rem;color:#66736d;font-size:.88rem;line-height:1.35">'
                    "Allow ProteinMPNN to change interface residues instead of holding the designed interface fixed.</div>",
                    unsafe_allow_html=True,
                )
            if modality == "homo_oligomer":
                binder_copies = st.number_input(
                    "Binder copies",
                    min_value=1,
                    max_value=24,
                    value=2,
                    step=1,
                    key="bindcraft2-binder-copies",
                    help="Number of binder copies (for example, 2 for a dimer or 3 for a trimer). This sets the oligomer size.",
                )
            if modality == "multidomain":
                number_of_domains = st.number_input(
                    "Number of domains",
                    min_value=2,
                    max_value=12,
                    value=2,
                    step=1,
                    key="bindcraft2-number-of-domains",
                    help="How many domains to build into the single binder chain. Each domain needs at least 50 residues.",
                )

            st.markdown("**Design options**")
            previous_design_properties = st.session_state.get("bindcraft2-design-properties", [])
            previous_design_properties = set(previous_design_properties) if isinstance(previous_design_properties, (list, tuple, set)) else set()
            design_properties = []
            for prop in BINDCRAFT2_PROPERTY_ORDER:
                property_key = f"bindcraft2-design-property-{prop}"
                if property_key not in st.session_state:
                    st.session_state[property_key] = prop in previous_design_properties
                if st.checkbox(
                    f"**{BINDCRAFT2_PROPERTY_LABELS[prop]}**",
                    key=property_key,
                    disabled=(prop == "mixed_topology" and modality in bindcraft2.BC2_SCAFFOLD_FORMATS),
                ):
                    design_properties.append(prop)
                st.markdown(
                    '<div style="margin:-.65rem 0 .65rem 1.85rem;color:#66736d;font-size:.88rem;line-height:1.35">'
                    f"{escape(BINDCRAFT2_PROPERTY_HELP[prop])}</div>",
                    unsafe_allow_html=True,
                )

            st.markdown("**Sequence composition**")
            excluded_amino_acids = st.multiselect(
                "Exclude amino acids",
                options=list(bindcraft2.BC2_AMINO_ACIDS),
                key="bindcraft2-excluded-amino-acids",
                help="Selected amino acids will not be used in the designed binder.",
            )
            amino_acid_bias_text = st.text_input(
                "Amino acid bias",
                placeholder="W:0.4, F:0.8",
                key="bindcraft2-amino-acid-bias",
                help="Optional relative propensities. 1.0 is neutral; values below 1 reduce and values above 1 increase use.",
            )
            try:
                amino_acid_bias = _parse_amino_acid_bias(amino_acid_bias_text)
            except ValueError as exc:
                optional_error = str(exc)

            if generator_only:
                quality_filters = {}
            else:
                st.markdown("**Candidate quality filters**")
                filter_cols = st.columns(3)
                with filter_cols[0]:
                    min_interface_iptm = st.number_input("Min interface i_pTM", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_interface_iptm"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-interface-iptm")
                    min_complex_ptm = st.number_input("Min complex pTM", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_complex_ptm"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-complex-ptm", placeholder="Disabled for this modality", help="Minimum confidence for the full target-binder complex. The upstream peptide presets disable this filter.")
                    min_interface_residues = st.number_input("Min interface residues", min_value=0, max_value=500, value=quality_filter_defaults["min_interface_residues"], step=1, key=f"bindcraft2-{modality}-min-interface-residues")
                with filter_cols[1]:
                    max_interface_ipae = st.number_input("Max interface i_pAE", min_value=0.0, max_value=1.0, value=quality_filter_defaults["max_interface_ipae"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-max-interface-ipae")
                    min_unbound_binder_plddt = st.number_input("Min unbound binder pLDDT", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_unbound_binder_plddt"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-unbound-plddt", placeholder="Disabled for this modality", help="Confidence of the binder on its own. The upstream peptide presets disable this filter.")
                    max_backbone_clashes = st.number_input("Max backbone clashes", min_value=0, max_value=10000, value=quality_filter_defaults["max_backbone_clashes"], step=1, key=f"bindcraft2-{modality}-max-backbone-clashes", help="The shared BindCraft2 default is zero backbone clashes.")
                with filter_cols[2]:
                    max_binder_rmsd = st.number_input(
                        "Max binder RMSD",
                        min_value=0.0,
                        max_value=100.0,
                        value=quality_filter_defaults["max_binder_rmsd"],
                        step=0.1,
                        format="%.1f",
                        key=f"bindcraft2-{modality}-{conformational_objective}-max-binder-rmsd",
                        placeholder="No limit",
                        help=(
                            "The conformational objective presets disable this filter by default. "
                            "Without an objective, binder, large-binder, and homo-oligomer presets default to 3.5 Å. Clear to disable."
                        ),
                    )

                with st.expander("Stage confidence filters", expanded=False):
                    plddt_cols = st.columns(3)
                    with plddt_cols[0]:
                        min_plddt_screen = st.number_input("Min pLDDT · screen", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_plddt_screen"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-plddt-screen")
                        min_plddt_refine = st.number_input("Min pLDDT · refine", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_plddt_refine"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-plddt-refine")
                    with plddt_cols[1]:
                        min_plddt_anneal = st.number_input("Min pLDDT · anneal", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_plddt_anneal"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-plddt-anneal")
                        min_plddt_harden = st.number_input("Min pLDDT · harden", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_plddt_harden"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-plddt-harden")
                    with plddt_cols[2]:
                        min_plddt_mutate = st.number_input("Min pLDDT · mutate", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_plddt_mutate"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-plddt-mutate")
                        min_plddt_final = st.number_input("Min pLDDT · final", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_plddt_final"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-plddt-final")
                    iptm_cols = st.columns(3)
                    with iptm_cols[0]:
                        min_iptm_anneal = st.number_input("Min i_pTM · anneal", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_iptm_anneal"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-iptm-anneal")
                    with iptm_cols[1]:
                        min_iptm_harden = st.number_input("Min i_pTM · harden", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_iptm_harden"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-iptm-harden")
                    with iptm_cols[2]:
                        min_iptm_mutate = st.number_input("Min i_pTM · mutate", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_iptm_mutate"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-iptm-mutate")

                if modality == "multidomain":
                    st.markdown("**Multidomain filters**")
                    domain_filter_cols = st.columns(2)
                    with domain_filter_cols[0]:
                        max_interdomain_contact_final = st.number_input("Max interdomain contacts", min_value=0.0, max_value=1.0, value=quality_filter_defaults["max_interdomain_contact_final"], step=0.01, format="%.2f", key=f"bindcraft2-{modality}-max-interdomain-contact", placeholder="No limit", help="Maximum fraction of contacts between binder domains. Clear to disable.")
                    with domain_filter_cols[1]:
                        min_domain_separation_ratio_final = st.number_input("Min domain separation ratio", min_value=0.0, max_value=10.0, value=quality_filter_defaults["min_domain_separation_ratio_final"], step=0.1, format="%.1f", key=f"bindcraft2-{modality}-min-domain-separation", placeholder="No limit", help="Minimum separation between the designed binder domains. Clear to disable.")
                elif modality in {"VHH", "ARP", "scFv", "Fab"}:
                    st.markdown("**Scaffold filters**")
                    scaffold_filter_cols = st.columns(2)
                    with scaffold_filter_cols[0]:
                        min_scaffold_sequence_retained_final = st.number_input("Min scaffold sequence retained", min_value=0.0, max_value=1.0, value=quality_filter_defaults["min_scaffold_sequence_retained_final"], step=0.05, format="%.2f", key=f"bindcraft2-{modality}-min-scaffold-sequence-retained", placeholder="No limit", help="Minimum fraction of the scaffold sequence retained in the designed binder. Clear to disable.")
                    with scaffold_filter_cols[1]:
                        max_scaffold_framework_rmsd_final = st.number_input("Max scaffold framework RMSD", min_value=0.0, max_value=20.0, value=quality_filter_defaults["max_scaffold_framework_rmsd_final"], step=0.1, format="%.1f", key=f"bindcraft2-{modality}-max-scaffold-framework-rmsd", placeholder="No limit", help="Maximum structural deviation of the scaffold framework.")

                quality_filters = {
                    "filters": {
                        "i_pTM": {"threshold": float(min_interface_iptm), "higher": True},
                        "pTM": {"threshold": float(min_complex_ptm) if min_complex_ptm is not None else None, "higher": True},
                        "i_pAE": {"threshold": float(max_interface_ipae), "higher": False},
                        "Unbound_Binder_pLDDT": {
                            "threshold": float(min_unbound_binder_plddt) if min_unbound_binder_plddt is not None else None,
                            "higher": True,
                        },
                        "Interface_Residues": {"threshold": int(min_interface_residues), "higher": True},
                        "Backbone_Clashes": {"threshold": int(max_backbone_clashes), "higher": False},
                        "Binder_RMSD": {
                            "threshold": float(max_binder_rmsd) if max_binder_rmsd is not None else None,
                            "higher": False,
                        },
                    },
                    "min_plddt_screen": float(min_plddt_screen),
                    "min_plddt_refine": float(min_plddt_refine),
                    "min_plddt_anneal": float(min_plddt_anneal),
                    "min_plddt_harden": float(min_plddt_harden),
                    "min_plddt_mutate": float(min_plddt_mutate),
                    "min_plddt_final": float(min_plddt_final),
                    "min_iptm_anneal": float(min_iptm_anneal),
                    "min_iptm_harden": float(min_iptm_harden),
                    "min_iptm_mutate": float(min_iptm_mutate),
                    "max_interdomain_contact_final": (
                        float(max_interdomain_contact_final)
                        if max_interdomain_contact_final is not None
                        else None
                    ),
                    "min_domain_separation_ratio_final": (
                        float(min_domain_separation_ratio_final)
                        if min_domain_separation_ratio_final is not None
                        else None
                    ),
                    "min_scaffold_sequence_retained_final": (
                        float(min_scaffold_sequence_retained_final)
                        if min_scaffold_sequence_retained_final is not None
                        else None
                    ),
                    "max_scaffold_framework_rmsd_final": (
                        float(max_scaffold_framework_rmsd_final)
                        if max_scaffold_framework_rmsd_final is not None
                        else None
                    ),
                }

            st.caption("Tamarind’s maximum run time is a service-side limit; this local BindCraft2 build does not expose a matching setting.")
            resource_cols = st.columns(3)
            with resource_cols[0]:
                gpu_device = st.text_input("GPU device", value="0", key="bindcraft2-gpu-device", help="Use a device ID such as 0, or enter cpu.")
            with resource_cols[1]:
                cpu_cores = st.number_input("CPU cores", min_value=1, max_value=256, value=4, step=1, key="bindcraft2-cpu-cores")
            with resource_cols[2]:
                workers_per_gpu = st.selectbox(
                    "Workers per GPU",
                    options=list(bindcraft2.WORKERS_PER_GPU_OPTIONS),
                    format_func=lambda value: "Automatic" if value == "auto" else str(value),
                    key="bindcraft2-workers-per-gpu",
                )
            campaign_seed = st.number_input(
                "Random seed",
                min_value=0,
                max_value=2_147_483_647,
                value=0,
                step=1,
                key="bindcraft2-campaign-seed",
            )

            with st.expander("Advanced BindCraft2 settings", expanded=False):
                st.caption(
                    "Optional JSON overrides for BindCraft2 campaign settings. Any setting accepted by the selected engine image can be supplied here; "
                    "these values override the visible controls and modality presets."
                )
                engine_overrides_text = st.text_area(
                    "Engine settings JSON",
                    placeholder='{"design_recycles": 2, "sequence_candidates": 20}',
                    height=120,
                    key="bindcraft2-engine-overrides",
                )
                if engine_overrides_text.strip():
                    try:
                        decoded_overrides = json.loads(engine_overrides_text)
                        if not isinstance(decoded_overrides, dict):
                            raise ValueError("Advanced BindCraft2 settings must be a JSON object.")
                        engine_overrides = decoded_overrides
                    except (json.JSONDecodeError, ValueError) as exc:
                        optional_error = f"Invalid advanced BindCraft2 settings: {exc}"

        settings_preview = None
        settings_error = optional_error or target_error
        target_is_ready = (
            target_path is not None and bool(target_chains) and design_target_path is not None and valid_trim_ranges
            if target_type == "structure"
            else bool(target_sequence)
        )
        if target_is_ready and not settings_error:
            try:
                settings_preview = bindcraft2.build_bindcraft2_settings(
                    target_chains=target_chains,
                    binder_length=binder_length,
                    target_type=target_type,
                    target_sequence=target_sequence,
                    target_name=target_name or "app_target",
                    target_crop_lengths=target_crop_lengths,
                    target_crop_count=target_crop_count,
                    target_filename=(
                        "target.fasta"
                        if target_type == "sequence"
                        else "target.cif"
                        if target_path is not None and target_path.suffix.lower() in {".cif", ".mmcif"}
                        else "target.pdb"
                    ),
                    hotspots=hotspots,
                    coldspots=coldspots,
                    campaign_name=job_name,
                    modality=[modality, conformational_objective]
                    if conformational_objective != "none"
                    else [modality],
                    design_properties=design_properties,
                    number_of_final_designs=int(number_of_final_designs),
                    max_trajectories=int(max_trajectories),
                    campaign_seed=int(campaign_seed),
                    workers_per_gpu=workers_per_gpu,
                    redesign_interface=redesign_interface,
                    binder_copies=binder_copies if modality == "homo_oligomer" else None,
                    number_of_domains=number_of_domains if modality == "multidomain" else None,
                    second_target_filename="second_target.pdb" if use_second_target and second_target_path else None,
                    second_target_chains=second_target_chains,
                    second_target_hotspots=second_target_hotspots,
                    second_target_coldspots=second_target_coldspots,
                    second_target_goal=second_target_goal,
                    second_target_strength=second_target_strength,
                    excluded_amino_acids=excluded_amino_acids,
                    amino_acid_bias=amino_acid_bias,
                    quality_filters=quality_filters,
                    engine_overrides=engine_overrides,
                    trajectory_only=generator_only,
                )
            except (TypeError, ValueError) as exc:
                settings_error = str(exc)
        if target_type == "structure" and use_second_target and (second_target_path is None or not second_target_chains):
            settings_error = "Load a second target structure and select at least one chain, or turn off Add a second target."
        if target_type == "structure" and use_second_target and not second_target_trim_valid:
            settings_error = "Choose a second-target trim range that contains at least one protein residue."

        if settings_error:
            st.error(settings_error)
        with st.expander("Review settings to send", expanded=False):
            if settings_preview is not None and not settings_error:
                st.json(settings_preview)
            else:
                st.caption("Choose valid target inputs to preview the BindCraft2 settings.")

        ready = target_is_ready and bool(job_name) and bool(project_id) and not settings_error
        if binder_length and binder_min >= binder_max:
            st.warning("Maximum binder length must be larger than minimum binder length.")
            ready = False
        if target_type == "structure" and not target_path:
            st.info("Fetch or upload a PDB/mmCIF target structure to continue.")
            if target_type == "structure" and target_path and not target_chains:
                st.info("Select at least one target chain to continue.")
        if target_type == "sequence" and not target_sequence and not target_error:
            st.info("Upload a named FASTA file to continue.")
        if not job_name:
            st.caption("Add a job name before submitting.")
        if not project_id:
            st.info("Select or create a project before submitting this design.")

        prioritize_job = st.checkbox(
            "Prioritize in scheduler queue",
            key="bindcraft2-prioritize-queue",
            help="Moves this design ahead of normal queued work when resources become available. Running jobs are never interrupted.",
        )

        if st.button("Submit trajectories" if generator_only else "Submit design", type="primary", use_container_width=True, disabled=not ready, key=f"bindcraft2-submit-{int(generator_only)}"):
            parameters = {
                "target_type": target_type,
                "target_fasta": None,
                "target_pdb": str(design_target_path) if target_type == "structure" else None,
                "target_chains": target_chains,
                "binder_length": binder_length,
                "target_crop_lengths": target_crop_lengths,
                "target_crop_count": target_crop_count,
                "hotspots": hotspots,
                "coldspots": coldspots,
                "campaign_name": job_name,
                "modality": [modality, conformational_objective]
                if conformational_objective != "none"
                else [modality],
                "design_properties": design_properties,
                "number_of_final_designs": int(number_of_final_designs),
                "max_trajectories": int(max_trajectories),
                "campaign_seed": int(campaign_seed),
                "workers_per_gpu": workers_per_gpu,
                "second_target_pdb": str(second_target_design_path) if second_target_design_path is not None and use_second_target else None,
                "second_target_chains": second_target_chains,
                "second_target_hotspots": second_target_hotspots,
                "second_target_coldspots": second_target_coldspots,
                "second_target_goal": second_target_goal,
                "second_target_strength": float(second_target_strength),
                "redesign_interface": redesign_interface,
                "binder_copies": int(binder_copies) if modality == "homo_oligomer" else None,
                "number_of_domains": int(number_of_domains) if modality == "multidomain" else None,
                "excluded_amino_acids": excluded_amino_acids,
                "amino_acid_bias": amino_acid_bias,
                "quality_filters": quality_filters,
                "engine_overrides": engine_overrides,
                "cpu_cores": int(cpu_cores),
                "gpu_device": gpu_device.strip() or "0",
            }
            # The workflow selected by this tile fixes generator versus full
            # campaign mode at the API boundary.
            target_fasta_staging: Path | None = None
            payload = {
                "schema_version": schema["schema_version"],
                "workflow": schema["workflow"],
                "parameters": parameters,
                "project_id": project_id,
                "priority": 1000 if prioritize_job else 0,
            }
            try:
                if target_type == "sequence":
                    with tempfile.NamedTemporaryFile(
                        mode="w",
                        encoding="utf-8",
                        suffix=".fasta",
                        prefix="bindcraft2-target-",
                        dir=STUDIO_WORKSPACE / "tmp",
                        delete=False,
                    ) as staged_fasta:
                        staged_fasta.write(f">{target_name}\n{target_sequence}\n")
                        target_fasta_staging = Path(staged_fasta.name)
                    parameters["target_fasta"] = str(target_fasta_staging)
                validation = _post_local_api("/api/v1/workflows/validate", payload)
                if not validation.get("valid"):
                    st.error("The workflow inputs did not pass validation.")
                else:
                    result = _post_local_api("/api/v1/workflows", payload)
                    mark_studio_job(result, project_id=project_id, owner_id=str(CURRENT_USER["id"]))
                    job_reference = result.get("job", {})
                    queue_label = "Generator queued" if generator_only else "Design queued"
                    st.success(f"{queue_label} · Job {job_reference.get('job_code') or job_reference.get('run_id') or 'submitted'}")
                    st.session_state["last_submitted_job"] = result["job"]
            except (ValueError, TypeError, json.JSONDecodeError) as exc:
                st.error(f"Could not start BindCraft2: {exc}")
            finally:
                if target_fasta_staging is not None:
                    target_fasta_staging.unlink(missing_ok=True)


def _render_boltzgen_setup(schema: dict, *, generator_only: bool = False) -> None:
    if st.session_state.pop("boltzgen-example-job-pending", False):
        st.session_state[f"boltzgen-job-name-{int(generator_only)}"] = "boltzgen-7jzu-example"
    if st.button("← Back to All tools", key=f"back-from-boltzgen-{int(generator_only)}"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()

    st.markdown(
        """
        <style>
        .block-container { padding-top: 1.35rem; max-width: 1800px; }
        .boltzgen-tool-header { text-align: center; padding: .1rem 0 1rem; }
        .boltzgen-tool-header h1 { margin: 0 0 .3rem; font-size: 2rem; }
        .boltzgen-tool-meta { color: #738077; font-size: .84rem; }
        .boltzgen-tool-header p { color: #66736d; margin: .4rem 0; }
        div[data-testid="stButton"] button { border-radius: 999px; }
        div[data-testid="stButton"] button[kind="primary"] {
            background: #2f6b50 !important;
            border-color: #2f6b50 !important;
            color: #ffffff !important;
        }
        div[data-testid="stButton"] button[kind="primary"]:hover {
            background: #24553e !important;
            border-color: #24553e !important;
        }
        div[data-testid="stButton"] button:disabled {
            opacity: .58;
        }
        [role="radio"][aria-checked="true"] {
            background: #2f6b50 !important;
            border-color: #2f6b50 !important;
            color: #ffffff !important;
        }
        .boltzgen-runtime-card {
            border: 1px solid #dce5df;
            background: #f8faf8;
            border-radius: 12px;
            padding: .9rem 1.1rem;
            margin: .5rem 0 1rem;
        }
        .boltzgen-runtime-title { font-size: 1rem; font-weight: 650; color: #25352d; }
        .boltzgen-runtime-note { color: #738077; font-size: .86rem; margin-top: .2rem; }
        [data-baseweb="tag"] { background: #e5f0e9 !important; color: #28583f !important; }
        </style>
        <div class="boltzgen-tool-header">
          <h1>__TITLE__</h1>
          <div class="boltzgen-tool-meta">__META__</div>
          <p>Run BoltzGen protein, nanobody, peptide, small-molecule, and redesign workflows.</p>
        </div>
        """.replace("__TITLE__", "BoltzGen Generator" if generator_only else "BoltzGen")
        .replace("__META__", "Backbone generation only · sequence redesign and evaluation are separate steps" if generator_only else "De novo binder design · Long-running workflow"),
        unsafe_allow_html=True,
    )
    with st.expander("About this workflow", expanded=False):
        st.write(schema["description"])
        st.caption("Each task maps to a native BoltzGen protocol or native design YAML feature.")
        if generator_only:
            st.info("This tile runs the design stage only. It does not run inverse folding, refolding, analysis, filtering, or ranking.")
        st.markdown(
            "[BoltzGen source code](https://github.com/HannesStark/boltzgen) · "
            "[BoltzGen paper](https://hannes-stark.com/assets/boltzgen.pdf)"
        )
        st.caption("Uses the same BoltzGen container image as the Studio Cofolding workflow.")

    name_col, project_col = st.columns([1.25, 1.25])
    with name_col:
        job_name = st.text_input(
            "Job name",
            value="boltzgen-generator" if generator_only else "boltzgen",
            placeholder="For example, targetX binder design",
            key=f"boltzgen-job-name-{int(generator_only)}",
        ).strip()
    with project_col:
        project_id = _render_project_selector(f"boltzgen-{int(generator_only)}")

    st.subheader("Task")
    boltzgen_tasks = [
        ("protein", "Protein"),
        ("nanobody", "Nanobody"),
        ("small_molecule", "Small Molecule"),
        ("peptide", "Peptide"),
        ("cyclotide", "Cyclotide"),
        ("redesign", "Protein Redesign"),
        ("yaml", "YAML Configuration"),
    ]
    available_tasks = {task_id for task_id, _ in boltzgen_tasks}
    selected_task = st.session_state.get("boltzgen-selected-task", "protein")
    if selected_task not in available_tasks:
        selected_task = "protein"
    task_columns = st.columns(4, gap="small")
    for index, (task_id, task_label) in enumerate(boltzgen_tasks):
        with task_columns[index % 4]:
            if st.button(
                task_label,
                key=f"boltzgen-task-{task_id}",
                type="primary" if selected_task == task_id else "secondary",
                use_container_width=True,
            ):
                st.session_state["boltzgen-selected-task"] = task_id
                st.rerun()
    task_descriptions = {
        "protein": "De novo protein design with BoltzGen's protein-anything protocol.",
        "nanobody": "Design a single-domain antibody using BoltzGen's bundled nanobody scaffolds and nanobody-anything protocol.",
        "small_molecule": "Design a protein for a CCD ligand or SMILES input with protein-small_molecule.",
        "peptide": "Design a peptide binder with peptide-anything.",
        "cyclotide": "Design a cyclic peptide using BoltzGen's native cyclic sequence and disulfide constraints.",
        "redesign": "Redesign selected residues in the uploaded protein with BoltzGen's native file.design feature.",
        "yaml": "Edit a native BoltzGen design YAML and choose one of its four protocols.",
    }
    st.caption(task_descriptions[selected_task])

    st.markdown("#### Example inputs")
    example_col, example_help_col = st.columns([1.5, 4.5], vertical_alignment="center")
    with example_col:
        load_example = st.button(
            "Load example inputs",
            key="boltzgen-load-example",
            type="secondary",
            use_container_width=True,
        )
    with example_help_col:
        st.caption("Loads the 7JZU target and starter campaign values. Fetching a new example target can take a few seconds.")
    if load_example:
        try:
            request = Request(
                "https://files.rcsb.org/download/7JZU.pdb",
                headers={"User-Agent": "MN-Molecular-Design-Studio/0.1"},
            )
            with urlopen(request, timeout=30) as response:
                content = response.read(50_000_001)
            if len(content) > 50_000_000:
                raise ValueError("The example structure is larger than the 50 MB import limit.")
            example_root = user_input_root(str(CURRENT_USER["id"])) / "structures"
            example_root.mkdir(parents=True, exist_ok=True)
            example_path = example_root / f"7JZU-{hashlib.sha256(content).hexdigest()[:12]}.pdb"
            if not example_path.exists():
                example_path.write_bytes(content)
            st.session_state["boltzgen-fetched-target"] = str(example_path)
            st.session_state["boltzgen-fetched-target-name"] = "7JZU example"
            st.session_state["boltzgen-pdb-code"] = "7JZU"
            st.session_state["boltzgen-example-job-pending"] = True
            st.session_state["boltzgen-binder-length-min"] = 60
            st.session_state["boltzgen-binder-length-max"] = 180
            st.session_state["boltzgen-num-designs"] = 10
            st.session_state["boltzgen-budget"] = 2
            st.session_state["boltzgen-target-source-view"] = "my_files"
            st.session_state["boltzgen-selected-task"] = "protein"
            st.rerun()
        except (HTTPError, URLError, OSError, ValueError) as exc:
            st.error(f"Could not load the BoltzGen example: {exc}")

    if selected_task == "small_molecule":
        target_source = TargetSource(owner_id=str(CURRENT_USER["id"]))
        target_input = TargetDesignInput(label="Ligand", key_prefix="boltzgen-ligand")
        target_key_prefix = "boltzgen-ligand"
        st.subheader("Ligand")
        ligand_type = st.radio(
            "Ligand format",
            ["ccd", "smiles"],
            format_func=lambda value: "CCD code" if value == "ccd" else "SMILES",
            horizontal=True,
            key="boltzgen-ligand-type",
        )
        ligand_value = st.text_input(
            "CCD code" if ligand_type == "ccd" else "SMILES string",
            value="TSA" if ligand_type == "ccd" else "",
            key="boltzgen-ligand-value",
            help="Provide a PDB Chemical Component Dictionary code or a valid SMILES string.",
        ).strip()
    else:
        ligand_type = "ccd"
        ligand_value = ""
        st.subheader("Protein to redesign" if selected_task == "redesign" else "Target structure")
        structure_format = st.radio(
            "Structure format",
            options=["pdb", "cif"],
            format_func=lambda value: "PDB" if value == "pdb" else "CIF / mmCIF",
            horizontal=True,
            key="boltzgen-structure-format",
        )
        target_source_view = st.session_state.get("boltzgen-target-source-view", "my_files")
        source_buttons = st.columns([1, 1, 5], gap="small", vertical_alignment="center")
        with source_buttons[0]:
            if st.button(
                "My files",
                key="boltzgen-source-my-files",
                type="primary" if target_source_view == "my_files" else "secondary",
                use_container_width=True,
            ):
                st.session_state["boltzgen-target-source-view"] = "my_files"
                st.rerun()
        with source_buttons[1]:
            if st.button(
                "My proteins",
                key="boltzgen-source-my-proteins",
                type="primary" if target_source_view == "my_proteins" else "secondary",
                use_container_width=True,
            ):
                st.session_state["boltzgen-target-source-view"] = "my_proteins"
                st.rerun()

        saved_proteins: list[Path] = []
        structures_root = user_input_root(str(CURRENT_USER["id"])) / "structures"
        try:
            resolved_structures_root = structures_root.resolve()
            saved_proteins = []
            for candidate_path in structures_root.iterdir():
                if not candidate_path.is_file() or candidate_path.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
                    continue
                resolved_path = candidate_path.resolve()
                try:
                    resolved_path.relative_to(resolved_structures_root)
                except ValueError:
                    continue
                saved_proteins.append(resolved_path)
            saved_proteins = sorted(
                saved_proteins,
                key=lambda path: path.stat().st_mtime,
                reverse=True,
            )
        except OSError:
            saved_proteins = []

        if structure_format == "pdb":
            if target_source_view == "my_proteins":
                saved_pdb = [path for path in saved_proteins if path.suffix.lower() in {".pdb", ".ent"}]
                if saved_pdb:
                    selected_saved = st.selectbox(
                        "Saved protein structures",
                        options=[str(path) for path in saved_pdb],
                        format_func=lambda raw: Path(raw).name,
                        key="boltzgen-saved-protein-pdb",
                        help="Structures you previously uploaded or fetched in this account.",
                    )
                    selected_path = Path(selected_saved)
                    target_source = TargetSource(selected_path, selected_path.stem, str(CURRENT_USER["id"]))
                else:
                    target_source = TargetSource(owner_id=str(CURRENT_USER["id"]))
                    st.info("No saved protein structures yet. Choose My files to fetch or upload a target.")
            else:
                st.caption("Fetch a structure by PDB accession or upload a PDB file.")
                target_source = render_target_source(
                    key_prefix="boltzgen",
                    label="Target",
                    example_code="7JZU",
                    upload_label="Drop a PDB file here or browse",
                    upload_help="Upload a PDB (.pdb/.ent) structure or fetch it by PDB accession.",
                    upload_file_types=["pdb", "ent"],
                    owner_id=str(CURRENT_USER["id"]),
                )
            target_key_prefix = "boltzgen"
        else:
            if target_source_view == "my_proteins":
                saved_cif = [path for path in saved_proteins if path.suffix.lower() in {".cif", ".mmcif"}]
                if saved_cif:
                    selected_saved = st.selectbox(
                        "Saved mmCIF protein structures",
                        options=[str(path) for path in saved_cif],
                        format_func=lambda raw: Path(raw).name,
                        key="boltzgen-saved-protein-cif",
                    )
                    selected_path = Path(selected_saved)
                    target_source = TargetSource(selected_path, selected_path.stem, str(CURRENT_USER["id"]))
                else:
                    target_source = TargetSource(owner_id=str(CURRENT_USER["id"]))
                    st.info("No saved mmCIF protein structures yet. Choose My files to upload a target.")
            else:
                st.caption("Upload a CIF/mmCIF structure. Chain and hotspot selections use the author IDs in the file.")
                target_source = render_target_source(
                    key_prefix="boltzgen-cif",
                    label="Target",
                    upload_label="Drop a CIF/mmCIF file here or browse",
                    upload_help="Upload a CIF (.cif) or mmCIF (.mmcif) protein structure.",
                    allow_cif=True,
                    fetch_enabled=False,
                    upload_file_types=["cif", "mmcif"],
                    upload_in_expander=False,
                    owner_id=str(CURRENT_USER["id"]),
                )
            target_key_prefix = "boltzgen-cif"

        target_input = render_target_editor(
            target_source,
            key_prefix=target_key_prefix,
            label="Protein" if selected_task == "redesign" else "Target",
            prefer_chain_b=False,
            residue_selection_title="Residues to redesign" if selected_task == "redesign" else "Binding-site residues",
            residue_selection_caption=(
                "Select protein positions BoltzGen should redesign."
                if selected_task == "redesign"
                else "Selected residues condition the designed binder's binding site. This is soft conditioning."
            ),
            residue_selection_input_label="Chosen redesign residues" if selected_task == "redesign" else "Chosen residues",
        )

    if selected_task == "nanobody":
        st.subheader("Nanobody scaffold")
        scaffold_name = st.selectbox(
            "Bundled BoltzGen scaffold",
            options=["7eow", "7xl0", "8coh", "8z8v"],
            format_func=lambda value: {"7eow": "7EOW · nanobody scaffold", "7xl0": "7XL0 · nanobody scaffold", "8coh": "8COH · nanobody scaffold", "8z8v": "8Z8V · nanobody scaffold"}[value],
            key="boltzgen-nanobody-scaffold",
            help="Scaffold YAML and CIF files are bundled from the BoltzGen repository examples.",
        )
    else:
        scaffold_name = "7eow"
    if selected_task == "cyclotide":
        st.subheader("Cyclotide sequence")
        cyclotide_sequence = "3C8C6C5C3C1C2"
        st.code(cyclotide_sequence)
        st.caption("BoltzGen's native cyclic sequence notation with the three disulfide constraints from its cyclotide example.")
    else:
        cyclotide_sequence = "3C8C6C5C3C1C2"
    if selected_task == "peptide":
        st.subheader("Peptide length")
        peptide_length_min_col, peptide_length_max_col = st.columns(2)
        with peptide_length_min_col:
            peptide_min = int(st.number_input("Minimum peptide length", min_value=1, max_value=2000, value=12, key="boltzgen-peptide-min"))
        with peptide_length_max_col:
            peptide_max = int(st.number_input("Maximum peptide length", min_value=1, max_value=2000, value=20, key="boltzgen-peptide-max"))
        peptide_sequence = f"{peptide_min}..{peptide_max}"
    else:
        peptide_min, peptide_max, peptide_sequence = 12, 20, "12..20"
    if selected_task == "yaml":
        st.subheader("BoltzGen YAML")
        yaml_protocol = st.selectbox(
            "Protocol",
            ["protein-anything", "peptide-anything", "protein-small_molecule", "nanobody-anything"],
            key="boltzgen-yaml-protocol",
        )
        yaml_default = "entities:\n  - protein:\n      id: {{DESIGN_CHAIN}}\n      sequence: 60..180\n  - file:\n      path: {{TARGET_FILE}}\n      include:\n        - chain:\n            id: {{TARGET_CHAIN}}\n{{BINDING_TYPES}}      structure_groups: all\n"
        yaml_state_key = "boltzgen-design-yaml-v2"
        if yaml_state_key not in st.session_state:
            previous_yaml = str(st.session_state.get("boltzgen-design-yaml", ""))
            is_legacy_template = (
                previous_yaml.startswith("entities:")
                and "id: B\n      sequence: 60..180" in previous_yaml
                and "{{TARGET_FILE}}" in previous_yaml
                and "{{TARGET_CHAIN}}" in previous_yaml
                and "{{BINDING_TYPES}}" in previous_yaml
            )
            if previous_yaml.strip() and not is_legacy_template:
                st.session_state[yaml_state_key] = previous_yaml
        design_yaml = st.text_area(
            "Design specification",
            value=st.session_state.get(yaml_state_key, yaml_default),
            height=220,
            key=yaml_state_key,
            help="TARGET_FILE, TARGET_CHAIN, DESIGN_CHAIN, and BINDING_TYPES are replaced from the selected inputs. BoltzGen validates the YAML before running.",
        )
    else:
        yaml_protocol = "protein-anything"
        design_yaml = ""

    if selected_task in {"protein", "small_molecule"}:
        st.subheader("Designed protein length")
        length_min_col, length_max_col = st.columns(2)
        with length_min_col:
            binder_min = int(
                st.number_input(
                    "Minimum length",
                    min_value=1,
                    max_value=2000,
                    value=60,
                    step=1,
                    key="boltzgen-binder-length-min",
                )
            )
        with length_max_col:
            binder_max = int(
                st.number_input(
                    "Maximum length",
                    min_value=1,
                    max_value=2000,
                    value=180,
                    step=1,
                    key="boltzgen-binder-length-max",
                )
            )
        binder_length = f"{binder_min}-{binder_max}"
    elif selected_task == "peptide":
        binder_min, binder_max = peptide_min, peptide_max
        binder_length = peptide_sequence
    elif selected_task == "cyclotide":
        binder_min = binder_max = 32
        binder_length = "32"
    elif selected_task == "nanobody":
        binder_min = binder_max = 1
        binder_length = "1"
    elif selected_task == "redesign":
        binder_min = binder_max = 1
        binder_length = "1"
    else:
        binder_min, binder_max = 60, 180
        binder_length = "60-180"

    st.subheader("Generation" if generator_only else "Design campaign")
    design_col, budget_col = st.columns(2) if not generator_only else (st.container(), st.container())
    with design_col:
        num_designs = int(
            st.number_input(
                "Number of designs",
                min_value=1,
                max_value=500000,
                value=10,
                step=1,
                key="boltzgen-num-designs",
                help="How many designs BoltzGen should generate.",
            )
        )
    if generator_only:
        budget = num_designs
        st.caption("BoltzGen will save generated structures without running inverse folding, refolding, filtering, or ranking.")
    else:
        with budget_col:
            budget = int(
                st.number_input(
                    "Budget",
                    min_value=1,
                    max_value=500000,
                    value=2,
                    step=1,
                    key="boltzgen-budget",
                    help="How many top-ranked designs to keep in the full pipeline.",
                )
            )

    with st.expander("Optional settings", expanded=False):
        run_vanilla_pipeline = not generator_only
        if generator_only:
            st.caption("Generation-only mode is fixed for this tile.")
        else:
            st.caption("Native Boltz2 refolding and ranking are enabled. The image does not include BoltzGen's optional ESMFold2 runtime.")
        if selected_task in {"protein", "nanobody", "peptide", "cyclotide", "yaml"}:
            st.markdown("#### Non-binding site")
            non_binding_residues = render_target_coldspots(
                target_input,
                key_prefix=target_key_prefix,
            )
        else:
            non_binding_residues = ""

        batch_col, omit_col = st.columns(2) if not generator_only else (st.container(), st.container())
        with batch_col:
            diffusion_batch_text = st.text_input(
                "Diffusion batch size",
                placeholder="Auto",
                key="boltzgen-diffusion-batch-size",
                help=(
                    "Number of diffusion samples generated per GPU pass. Leave blank to use BoltzGen's default: "
                    "1 for fewer than 100 designs and 10 for 100 or more. Designs in a batch share a sampled length."
                ),
            ).strip()
        if generator_only:
            omit_amino_acids = ""
        else:
            with omit_col:
                omit_amino_acids = st.text_input(
                    "Omit amino acids",
                    placeholder="For example, C",
                    key="boltzgen-omit-amino-acids",
                    help="One-letter amino-acid codes excluded during inverse folding. Leave blank to use BoltzGen's protocol default.",
                ).strip().upper()

        if run_vanilla_pipeline:
            can_filter_adherence = selected_task in {"protein", "nanobody", "peptide", "cyclotide", "yaml"}
            if can_filter_adherence:
                filter_binding_site_adherence = st.checkbox(
                    "Filter by binding-site adherence",
                    value=False,
                    key="boltzgen-filter-bindingsite",
                    disabled=not bool(target_input.hotspots),
                    help=(
                        "Enable BoltzGen's native filter, which keeps designs with a residue within 4 Å of a selected binding-site residue. "
                        "Requires binding-site residues and the full pipeline."
                    ),
                )
                if not target_input.hotspots:
                    filter_binding_site_adherence = False
                    st.caption("Select binding-site residues above to enable BoltzGen's native adherence filter.")
            else:
                filter_binding_site_adherence = False
            st.markdown("#### Additional filters")
            filter_rows_key = "boltzgen-additional-filter-rows"
            if filter_rows_key not in st.session_state:
                st.session_state[filter_rows_key] = [0]
            filter_rows = st.session_state[filter_rows_key]
            additional_filters: list[str] = []
            for row_number, row_id in enumerate(filter_rows, start=1):
                filter_col, remove_col = st.columns([6, 1])
                with filter_col:
                    filter_text = st.text_input(
                        f"Filter {row_number}",
                        placeholder="For example, design_ipsae_min>0.5",
                        key=f"boltzgen-additional-filter-{row_id}",
                        help="BoltzGen syntax: metric>threshold or metric<threshold. Every filter must pass.",
                    ).strip()
                if filter_text:
                    additional_filters.append(filter_text)
                with remove_col:
                    if len(filter_rows) > 1 and st.button(
                        "Remove",
                        key=f"boltzgen-remove-filter-{row_id}",
                        help="Remove this BoltzGen filter.",
                    ):
                        st.session_state[filter_rows_key] = [item for item in filter_rows if item != row_id]
                        st.rerun()
            if st.button("+ Add filter", key="boltzgen-add-filter"):
                st.session_state[filter_rows_key] = [*filter_rows, max(filter_rows, default=-1) + 1]
                st.rerun()
            st.caption(
                "BoltzGen computes ipSAE during refolding; this optimized image is configured to save full pairwise PAE matrices."
            )
        else:
            filter_binding_site_adherence = False
            additional_filters = []

        sampling_steps = int(
            st.number_input(
                "Sampling steps",
                min_value=1,
                max_value=1000,
                value=20,
                step=1,
                key="boltzgen-sampling-steps",
                help="Diffusion sampling steps used during design.",
            )
        )
        gpu_device = st.text_input(
            "GPU device",
            value="0",
            key="boltzgen-gpu-device",
            help="GPU device ID used by the scheduler and container, such as 0.",
        ).strip()

    settings_error = ""
    diffusion_batch_size: int | None = None
    if diffusion_batch_text:
        if not diffusion_batch_text.isdigit():
            settings_error = "Diffusion batch size must be a positive whole number."
        else:
            diffusion_batch_size = int(diffusion_batch_text)
            if diffusion_batch_size < 1:
                settings_error = "Diffusion batch size must be a positive whole number."
    if not settings_error:
        if binder_min > binder_max:
            settings_error = "Maximum binder length must be at least the minimum binder length."
        elif budget > num_designs:
            settings_error = "Final design budget cannot exceed the number of designs."
        elif not gpu_device:
            settings_error = "Enter a GPU device ID."
        elif omit_amino_acids and (
            not all(letter in "ACDEFGHIKLMNPQRSTVWY" for letter in omit_amino_acids)
            or len(set(omit_amino_acids)) != len(omit_amino_acids)
            or len(omit_amino_acids) == 20
        ):
            settings_error = "Omit amino acids must use unique standard codes and leave at least one amino acid available."
        elif filter_binding_site_adherence and not target_input.hotspots:
            settings_error = "Select at least one binding-site residue for BoltzGen's native adherence filter."
        elif additional_filters:
            for expression in additional_filters:
                operator = ">" if ">" in expression else "<" if "<" in expression else ""
                feature, _, threshold = expression.partition(operator) if operator else ("", "", "")
                try:
                    valid_filter = bool(feature.strip()) and float(threshold) == float(threshold)
                except ValueError:
                    valid_filter = False
                if not valid_filter:
                    settings_error = (
                        "Each BoltzGen filter must use metric>number or metric<number syntax."
                    )
                    break
    if (
        not settings_error
        and set(filter(None, target_input.hotspots.split(",")))
        & set(filter(None, non_binding_residues.split(",")))
    ):
        settings_error = "A residue cannot be selected as both binding and non-binding."

    target_ready = bool(ligand_value) if selected_task == "small_molecule" else target_input.ready
    if selected_task == "redesign" and not target_input.hotspots:
        target_ready = False
    if selected_task == "yaml" and not design_yaml.strip():
        target_ready = False
    ready = (
        selected_task in available_tasks
        and target_ready
        and bool(job_name)
        and bool(project_id)
        and not settings_error
    )
    if selected_task == "small_molecule" and not ligand_value:
        st.info("Enter a CCD code or SMILES string to continue.")
    elif selected_task != "small_molecule" and target_source.path is None:
        st.info("Fetch or upload a target structure to continue.")
    elif selected_task != "small_molecule" and not target_input.chains:
        st.info("Select at least one target chain to continue.")
    if selected_task == "redesign" and not target_input.hotspots:
        st.info("Select at least one residue to redesign.")
    if selected_task == "yaml" and not design_yaml.strip():
        st.info("Enter a BoltzGen YAML design specification to continue.")
    if not project_id:
        st.info("Select or create a project before submitting this design.")
    if settings_error:
        st.error(settings_error)

    try:
        from mn_molecular_design_studio.engine.workflow_engine.core.runtime_estimator import estimate_job_runtime

        runtime = estimate_job_runtime(
            {"tool": "boltzgen_generator" if generator_only else "boltzgen", "job_type": "design_generation" if generator_only else "design_campaign"},
            {"params": {"num_designs": num_designs}},
        )
        runtime_label = runtime.get("estimated_time", "n/a")
        runtime_basis = "Estimated from the current campaign size; actual time depends on GPU and queue load."
    except (ImportError, TypeError, ValueError, KeyError):
        runtime_label = "Varies"
        runtime_basis = "Actual time depends on GPU, campaign size, and queue load."
    st.markdown(
        f'<div class="boltzgen-runtime-card"><div class="boltzgen-runtime-title">Estimated runtime · {escape(str(runtime_label))}</div>'
        f'<div class="boltzgen-runtime-note">{escape(runtime_basis)}</div></div>',
        unsafe_allow_html=True,
    )

    parameters = {
        "task": selected_task,
        "target_chains": target_input.chains,
        "binder_length": binder_length,
        "hotspots": target_input.hotspots,
        "campaign_name": job_name,
        "num_designs": num_designs,
        "budget": budget,
        "sampling_steps": sampling_steps,
        "non_binding_residues": non_binding_residues,
        "diffusion_batch_size": diffusion_batch_size,
        "omit_amino_acids": omit_amino_acids,
        "additional_filters": additional_filters,
        "filter_binding_site_adherence": filter_binding_site_adherence,
        "scaffold_name": scaffold_name,
        "ligand_type": ligand_type,
        "ligand_value": ligand_value,
        "peptide_sequence": peptide_sequence,
        "cyclotide_sequence": cyclotide_sequence,
        "redesign_residues": target_input.hotspots if selected_task == "redesign" else "",
        "design_yaml": design_yaml,
        "yaml_protocol": yaml_protocol,
    }
    if target_input.design_path:
        parameters["target_pdb"] = str(target_input.design_path)
    with st.expander("Review settings to send", expanded=False):
        if ready:
            st.json(
                {
                    "workflow": schema["workflow"],
                    "parameters": parameters,
                    "resources": {"gpu_device": gpu_device},
                    "project": project_id,
                }
            )
        else:
            st.caption("Choose valid target, project, and campaign settings to preview the request.")

    prioritize_job = st.checkbox(
        "Prioritize in scheduler queue",
        key="boltzgen-prioritize-queue",
        help="Moves this design ahead of normal queued work when resources become available. Running jobs are never interrupted.",
    )
    if st.button("Submit generator job" if generator_only else "Submit pipeline job", type="primary", use_container_width=True, disabled=not ready, key=f"boltzgen-submit-{int(generator_only)}"):
        payload = {
            "schema_version": schema["schema_version"],
            "workflow": schema["workflow"],
            "parameters": parameters,
            "resources": {"gpu_device": gpu_device},
            "project_id": project_id,
            "priority": 1000 if prioritize_job else 0,
        }
        try:
            validation = _post_local_api("/api/v1/workflows/validate", payload)
            if not validation.get("valid"):
                st.error("The workflow inputs did not pass validation.")
            else:
                result = _post_local_api("/api/v1/workflows", payload)
                mark_studio_job(result, project_id=project_id, owner_id=str(CURRENT_USER["id"]))
                job_reference = result.get("job", {})
                queue_label = "Generator queued" if generator_only else "Design pipeline queued"
                st.success(f"{queue_label} · Job {job_reference.get('job_code') or job_reference.get('run_id') or 'submitted'}")
                st.session_state["last_submitted_job"] = result["job"]
        except (ValueError, TypeError, json.JSONDecodeError) as exc:
            st.error(f"Could not start BoltzGen: {exc}")


def _post_local_api(path: str, payload: dict) -> dict:
    token = str(st.session_state.get("studio-session-token") or ("development-mode" if is_development_mode() else ""))
    if not token:
        raise ValueError("Your sign-in session expired. Sign in again before submitting jobs.")
    pending_key = None
    if cluster_config.enabled() and path in {"/api/v1/workflows", "/api/v1/workflows/validate", "/api/v1/pipelines/run", "/api/v1/pipelines/validate"}:
        payload = dict(payload)
        fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()
        pending_key = "cluster-submission-" + fingerprint
        key = st.session_state.setdefault(pending_key, secrets.token_urlsafe(24))
        payload.setdefault("execution", {"placement": st.session_state.get("cluster-placement", cluster_config.load()["default_placement"]), "idempotency_key": key})
    request = Request(
        f"{API_BASE_URL}{path}",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "Authorization": f"Bearer {token}"},
        method="POST",
    )
    try:
        with urlopen(request, timeout=60) as response:
            result = json.loads(response.read().decode("utf-8"))
            if pending_key and result.get("job"):
                st.session_state.pop(pending_key, None)
            return result
    except HTTPError as exc:
        try:
            error = json.loads(exc.read().decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            error = {"detail": str(exc)}
        raise ValueError(error.get("detail", str(exc))) from exc
    except URLError as exc:
        raise ValueError(f"Could not reach the local workflow API at {API_BASE_URL}: {exc.reason}") from exc


def _render_parameter_fields(parameters: list[dict], *, workflow_id: str, group: str) -> dict:
    values = {}
    columns = st.columns(2, gap="large")
    for index, parameter in enumerate(parameters):
        with columns[index % 2]:
            raw_value = _parameter_widget(parameter, workflow_id=workflow_id, group=group)
        values[parameter["name"]] = raw_value
    return values


def _analysis_source_job_options(project_id: str) -> list[dict]:
    """List completed, accessible complex-result jobs from one selected project."""
    result = search_jobs(statuses=["completed", "succeeded"], project_id=project_id, limit=1000)
    options: list[dict] = []
    for row in result.get("items", []):
        run_dir = Path(str(row.get("run_dir") or ""))
        candidates_path = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
        if not candidates_path.is_file():
            continue
        candidates = read_candidates(candidates_path)
        complex_candidates = [candidate for candidate in candidates if candidate.get("complex_pdb")]
        if not complex_candidates:
            continue
        reference = f"{row.get('task_group', '')}/{row.get('run_id', '')}"
        name = str(row.get("campaign_name") or row.get("tool") or row.get("job_type") or "Results")
        options.append(
            {
                "reference": reference,
                "label": f"{name} · {row.get('job_code') or row.get('run_id')} · {len(complex_candidates)} complex candidate(s)",
                "row": row,
            }
        )
    return options


def _source_job_pae_availability(source_job: dict) -> tuple[int, int]:
    """Count source complex candidates with a saved pairwise PAE matrix file."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _find_pae_path, _resolve_path

    run_dir = Path(str(source_job.get("run_dir") or ""))
    candidates_path = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidates = [candidate for candidate in read_candidates(candidates_path) if candidate.get("complex_pdb")]
    available = 0
    for candidate in candidates:
        structure = _resolve_path(run_dir, str(candidate.get("complex_pdb") or ""))
        pae_path = _find_pae_path(run_dir, candidate, structure)
        if pae_path and pae_path.is_file():
            available += 1
    return available, len(candidates)


def _render_esmfold2_cofolding_setup(schema: dict) -> None:
    _render_family_cofolding_setup(schema, family="esmfold2")


def _af2_pdb_chain_sequences(source: TargetSource) -> dict[str, str]:
    if not source.path or not Path(source.path).is_file():
        return {}
    try:
        entries = _sequence_entries_by_chain(Path(source.path).read_text(errors="replace"))
    except OSError:
        return {}
    return {chain: "".join(residue for _number, residue in residues) for chain, residues in entries.items()}


def _render_af2_pdb_chain_selector(
    source: TargetSource,
    *,
    key: str,
    label: str,
    default_chains: list[str] | None = None,
) -> tuple[list[str], dict[str, str]]:
    sequences = _af2_pdb_chain_sequences(source)
    if not sequences:
        if source.path:
            st.error(f"No protein chains could be read from {label}. Upload a PDB containing protein ATOM records.")
        return [], {}
    digest = hashlib.sha1(str(source.path).encode()).hexdigest()[:10]
    selected = st.multiselect(
        label,
        list(sequences),
        default=[chain for chain in (default_chains or list(sequences)[:1]) if chain in sequences],
        format_func=lambda chain: f"Chain {chain} · {len(sequences[chain]):,} residues",
        key=f"{key}-{digest}",
        help="Select one or more protein chains from the PDB.",
    )
    return selected, sequences


def _render_af2_initial_guess_setup(schema: dict) -> None:
    title = "AlphaFold2 Initial Guess"
    if st.button("← Back to All tools", key="back-from-af2-initial-guess"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    _render_cofolding_header(
        title,
        "Predict protein structures from PDB templates and optional partner sequences.",
    )
    with st.expander("About this workflow", expanded=False):
        st.write(
            "Use one target PDB, separate target and partner PDBs, or a complex PDB. "
            "Batch inputs in a job must use the same template type and AF2 parameter family."
        )

    name_col, project_col = st.columns([1.25, 1.25])
    with name_col:
        job_name = st.text_input("Job name", value="alphafold2-initial-guess", key="af2-ig-job-name").strip()
    with project_col:
        project_id = _render_project_selector("af2-initial-guess")
    validation_profile = {
        "label": "AlphaFold2 Initial Guess", "modes": ("off",), "input_kind": "af2_initial_guess",
        "default_recycles": 3, "supports_templates": True, "requires_target_template": True, "supports_msa": True,
        "supports_seed": False, "default_samples": 1,
        "settings": ({"key": "multimer", "type": "bool", "label": "Multimer parameters", "default": True},
                     {"key": "max_msa_seqs", "label": "Maximum target MSA rows", "default": 1024, "min": 1, "max": 4096}),
    }
    def submit_af2_validation(_schema, _workflow, parameters, selected_project, settings, label):
        values = {
            "validation_sources": parameters["validation_sources"], "job_name": parameters["job_name"],
            "num_recycles": parameters["num_recycles"], "model_count": parameters["samples"],
            "multimer": parameters["engine_options"]["multimer"], "target_use_msa": parameters["target_use_msa"],
            "msa_max_depth": parameters["engine_options"]["max_msa_seqs"],
            "cpu_cores": parameters["cpu_cores"], "gpu_device": parameters["gpu_device"],
        }
        _submit_cofolding_job(workflow_schema("refolding.af2_design_validation"), "refolding.af2_design_validation",
                              values, selected_project, settings, label)
    handled, _ = render_batch_controls(
        "af2_initial_guess", validation_profile, project_id, job_name, schema,
        search_jobs=search_jobs, submit=submit_af2_validation, allow_new_batch=False)
    if handled:
        return

    mode_labels = {
        "Target PDB + optional partner sequence": "target_template",
        "Separate target and partner PDBs": "separate_templates",
        "One complex PDB": "complex_template",
    }
    mode_by_id = {value: label for label, value in mode_labels.items()}
    mode_state_key = "af2-ig-input-mode"
    st.markdown("**EXAMPLES**")
    example_cols = st.columns(2)
    with example_cols[0]:
        if st.button("7QF4 · MiniSOG monomer", key="af2-ig-example-minisog"):
            st.session_state[mode_state_key] = mode_by_id["target_template"]
            st.session_state["af2-ig-item-1-target-pdb-code"] = "7QF4"
            st.session_state["af2-ig-item-1-partner-sequence"] = ""
            st.rerun()
    with example_cols[1]:
        if st.button("8JYR · Antibody complex", key="af2-ig-example-antibody"):
            st.session_state[mode_state_key] = mode_by_id["complex_template"]
            st.session_state["af2-ig-item-1-complex-pdb-code"] = "8JYR"
            st.rerun()

    with st.expander("Multiple PDB inputs", expanded=False):
        batch_size = int(st.number_input(
            "Independent predictions", min_value=1, max_value=32, value=1, step=1,
            key="af2-ig-batch-size",
            help="Run up to 32 independent PDB inputs with the same AF2 settings in one job.",
        ))

    selected_mode = st.radio(
        "Template input",
        list(mode_labels),
        horizontal=True,
        key=mode_state_key,
        help="Choose the same input type for every prediction in this batch. Each prediction requires at least one PDB.",
    )
    input_mode = mode_labels[selected_mode]

    input_tabs = st.tabs([f"Prediction {index + 1}" for index in range(batch_size)]) if batch_size > 1 else [st.container()]
    batch_inputs: list[dict[str, Any]] = []
    input_validity: list[bool] = []
    partner_presence: list[bool] = []
    owner_id = str(CURRENT_USER["id"])
    for index, tab in enumerate(input_tabs, start=1):
        with tab:
            input_row: dict[str, Any] = {"input_mode": input_mode}
            valid_input = False
            if input_mode == "target_template":
                target_pdb = render_target_source(
                    key_prefix=f"af2-ig-item-{index}-target",
                    label="Target structure",
                    example_code="7QF4",
                    upload_file_types=["pdb"],
                    owner_id=owner_id,
                )
                selected_target_chains, _target_sequences = _render_af2_pdb_chain_selector(
                    target_pdb, key=f"af2-ig-item-{index}-target-chain", label="Target chain(s)"
                )
                partner_sequence = st.text_area(
                    "Partner protein sequence",
                    height=105,
                    key=f"af2-ig-item-{index}-partner-sequence",
                    placeholder="Paste a partner sequence, or leave blank to refold only the target.",
                    help="Use protein sequence only. Ligands and nucleic acids are not supported.",
                ).strip().upper()
                input_row.update(
                    target_pdb=str(target_pdb.path) if target_pdb.path else None,
                    target_chains=selected_target_chains,
                    partner_sequence=partner_sequence,
                )
                valid_input = bool(target_pdb.path and selected_target_chains)
                if partner_sequence and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", partner_sequence) is None:
                    st.error(f"Prediction {index}: partner sequence must use the 20 standard amino-acid letters.")
                    valid_input = False
                partner_presence.append(bool(partner_sequence))
            elif input_mode == "separate_templates":
                target_col, partner_col = st.columns(2, gap="large")
                with target_col:
                    target_pdb = render_target_source(
                        key_prefix=f"af2-ig-item-{index}-target",
                        label="Target structure",
                        example_code="7QF4",
                        upload_file_types=["pdb"],
                        owner_id=owner_id,
                    )
                    selected_target_chains, _target_sequences = _render_af2_pdb_chain_selector(
                        target_pdb, key=f"af2-ig-item-{index}-target-chain", label="Target chain(s)"
                    )
                with partner_col:
                    partner_pdb = render_target_source(
                        key_prefix=f"af2-ig-item-{index}-partner",
                        label="Partner structure",
                        example_code="8JYR",
                        upload_file_types=["pdb"],
                        owner_id=owner_id,
                    )
                    selected_partner_chains, _partner_sequences = _render_af2_pdb_chain_selector(
                        partner_pdb, key=f"af2-ig-item-{index}-partner-chain", label="Partner chain(s)"
                    )
                input_row.update(
                    target_pdb=str(target_pdb.path) if target_pdb.path else None,
                    target_chains=selected_target_chains,
                    partner_pdb=str(partner_pdb.path) if partner_pdb.path else None,
                    partner_chains=selected_partner_chains,
                )
                valid_input = bool(target_pdb.path and selected_target_chains and partner_pdb.path and selected_partner_chains)
                partner_presence.append(True)
            else:
                complex_pdb = render_target_source(
                    key_prefix=f"af2-ig-item-{index}-complex",
                    label="Complex structure",
                    example_code="8JYR",
                    upload_file_types=["pdb"],
                    owner_id=owner_id,
                )
                available = _af2_pdb_chain_sequences(complex_pdb)
                if available:
                    chain_options = list(available)
                    digest = hashlib.sha1(str(complex_pdb.path).encode()).hexdigest()[:10]
                    chain_col, partner_col = st.columns(2, gap="large")
                    with chain_col:
                        selected_complex_targets = st.multiselect(
                            "Target chain(s)", chain_options, default=chain_options[:1],
                            format_func=lambda chain: f"Chain {chain} · {len(available[chain]):,} residues",
                            key=f"af2-ig-item-{index}-complex-targets-{digest}",
                        )
                    with partner_col:
                        selected_complex_binders = st.multiselect(
                            "Partner chain(s)", chain_options, default=chain_options[1:2],
                            format_func=lambda chain: f"Chain {chain} · {len(available[chain]):,} residues",
                            key=f"af2-ig-item-{index}-complex-binders-{digest}",
                        )
                    overlap = sorted(set(selected_complex_targets) & set(selected_complex_binders))
                    if overlap:
                        st.error(f"Prediction {index}: chain(s) selected as both target and partner: {', '.join(overlap)}.")
                    valid_input = bool(complex_pdb.path and selected_complex_targets and selected_complex_binders and not overlap)
                    input_row.update(
                        complex_pdb=str(complex_pdb.path) if complex_pdb.path else None,
                        complex_target_chains=selected_complex_targets,
                        complex_binder_chains=selected_complex_binders,
                    )
                elif complex_pdb.path:
                    st.error("No protein chains could be read from the complex PDB.")
                partner_presence.append(True)
            batch_inputs.append(input_row)
            input_validity.append(valid_input)

    target_only_batch = input_mode == "target_template" and partner_presence and not any(partner_presence)
    mixed_monomer_complex = input_mode == "target_template" and any(partner_presence) and not all(partner_presence)
    if mixed_monomer_complex:
        st.error("For one batch, either enter a partner sequence for every prediction or leave every partner sequence blank.")
    if target_only_batch:
        model_family = "ptm"
        st.caption("Target-only input uses AF2 monomer (PTM) parameters.")
    else:
        model_family = st.selectbox(
            "AF2 parameter set",
            ["multimer_v3", "ptm"],
            format_func=lambda value: "AF2 Multimer v3" if value == "multimer_v3" else "AF2 monomer (PTM)",
            key=f"af2-ig-model-family-{input_mode}",
            help="Multimer v3 is the usual setting for complexes; PTM uses the AF2 monomer parameter family.",
        )

    model_count = int(st.number_input(
        "Number of models *", min_value=1, max_value=5, value=1, step=1,
        key="af2-ig-model-count",
    ))
    with st.expander("Optional settings", expanded=False):
        num_recycles = int(st.number_input(
            "Number of recycles", min_value=1, max_value=20, value=3, step=1,
            key="af2-ig-num-recycles",
        ))
        msa_col, depth_col = st.columns(2, gap="large")
        with msa_col:
            msa_default = "reference_then_colabfold" if input_mode == "target_template" and not target_only_batch else "single_sequence"
            msa_choices = ["reference_then_colabfold", "single_sequence"]
            msa_source = st.selectbox(
                "Protein MSA source",
                msa_choices,
                index=msa_choices.index(msa_default),
                format_func=lambda value: (
                    "Reference folder, then ColabFold 1.6.1 opt"
                    if value == "reference_then_colabfold"
                    else "Single sequence (no MSA)"
                ),
                key=f"af2-ig-msa-source-{input_mode}-{'target-only' if target_only_batch else 'complex'}",
                help="Exact-sequence reference MSAs are reused first. Missing alignments are generated and saved for later jobs.",
            )
        with depth_col:
            msa_max_depth = int(st.number_input(
                "Maximum MSA rows per chain", min_value=1, max_value=4096, value=1024, step=32,
                key="af2-ig-msa-depth",
                disabled=msa_source == "single_sequence",
            ))

    st.caption("AF2 initial guess supports protein chains only. Ligands, DNA, RNA, and affinity scoring are not part of this workflow.")
    submitted = st.button(
        "Run AlphaFold2 prediction",
        type="primary",
        disabled=not bool(project_id) or not all(input_validity) or mixed_monomer_complex,
        key="af2-ig-submit",
    )
    if not submitted:
        return
    if not all(input_validity):
        st.error("Complete the PDB and chain selections for every prediction in the batch.")
        return

    parameters: dict[str, Any] = {
        "input_mode": input_mode,
        "batch_inputs": batch_inputs,
        "num_recycles": num_recycles,
        "model_count": model_count,
        "multimer": model_family == "multimer_v3",
        "msa_source": msa_source,
        "msa_max_depth": msa_max_depth,
        "docker_image": "mn-colabdesign:latest",
        "gpu_device": "0",
    }
    parameters["job_name"] = job_name or "alphafold2-initial-guess"
    try:
        _submit_cofolding_job(schema, "refolding.af2_initial_guess_input", parameters, project_id, {}, title)
    except (ValueError, TypeError, json.JSONDecodeError) as exc:
        st.error(f"Could not start {title}: {exc}")


def _resolve_seed_values(seed_text: str, *, label: str, max_count: int = 32) -> list[int]:
    """Treat zero as a request for a fresh, distinct random seed."""
    try:
        requested = [int(value.strip()) for value in str(seed_text or "").split(",") if value.strip()]
    except ValueError as exc:
        raise ValueError(f"{label} seeds must be integers separated by commas.") from exc
    if not requested or len(requested) > max_count or any(value < 0 or value > 2**31 - 1 for value in requested):
        raise ValueError(f"Enter 1 to {max_count} non-negative {label} seeds no greater than 2,147,483,647.")
    fixed = [value for value in requested if value]
    if len(set(fixed)) != len(fixed):
        raise ValueError(f"Enter distinct fixed {label} seeds; use 0 for independently randomized seeds.")
    resolved = []
    used = set(fixed)
    for value in requested:
        if value == 0:
            value = secrets.randbelow(2**31 - 1) + 1
            while value in used:
                value = secrets.randbelow(2**31 - 1) + 1
            used.add(value)
        resolved.append(value)
    return resolved




def _render_af3_jax_cofolding_setup(schema: dict) -> None:
    _render_family_cofolding_setup(schema, family="af3_jax")


def _render_family_cofolding_setup(schema: dict, *, family: str) -> None:
    is_jax = family == "af3_jax"
    title = "ColabFold 2.0 preview" if is_jax else "ESMFold2 Cofolding"
    workflow_id = "refolding.af3_jax_cofolding" if is_jax else "refolding.esmfold2_cofolding"
    if st.button("← Back to All tools", key=f"back-from-{family}-cofolding"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    _render_cofolding_header(title, "Choose a model and enter molecular chains for structure prediction.")
    with st.expander("About this workflow", expanded=False):
        st.write(schema["description"])
    name_col, project_col = st.columns([1.25, 1.25])
    with name_col:
        job_name = st.text_input("Job name", value=f"{family.replace('_', '-')}-cofold", key=f"{family.replace('_', '-')}-cofold-job-name").strip()
        if is_jax:
            model = st.selectbox("Model", list(af3_jax_cofolding.MODEL_PROFILES),
                                 index=list(af3_jax_cofolding.MODEL_PROFILES).index("alphafold3"),
                                 format_func=lambda value: af3_jax_cofolding.MODEL_PROFILES[value]["label"], key="mn-cofolding-model")
            profile = {**af3_jax_cofolding.MODEL_PROFILES[model], "default_recycles": 10}
        else:
            model = st.selectbox("Model", ["full", "fast"], format_func=lambda value: "ESMFold2 Full · MSA capable" if value == "full" else "ESMFold2 Fast · sequence only", key="esmfold2-cofold-model")
            profile = {"modes": ("off", "exact", "fast", "big"), "default_recycles": 10}
    with project_col:
        project_id = _render_project_selector(f"{family}-cofold")
    batch_engine = "alphafold3" if is_jax else "esmfold2_fast" if model == "fast" else "esmfold2"
    batch_profile = {**optimized_prediction.ENGINE_PROFILES[batch_engine], "form_key": f"{family}-{model}"}
    if is_jax:
        batch_profile.update(label=title, modes=tuple(profile["modes"]), terms=bool(profile.get("terms")),
                             gpu_count={"big": 2} if profile["two_gpu_big"] else {})
        batch_profile["settings"] = tuple(
            {**field, "default": model} if field["key"] == "model" else field
            for field in batch_profile.get("settings", ()))
    handled, family_batch_count = render_batch_controls(
        batch_engine, batch_profile, project_id, job_name, workflow_schema("refolding.optimized_prediction"),
        search_jobs=search_jobs, submit=_submit_cofolding_job, allow_new_batch=False)
    if handled:
        return
    extra_fields = [
        {"key": "seed_count", "label": "Number of seeds", "default": 1, "min": 1, "max": 32,
         "help": "Run each seed with the requested number of samples. Zero as the starting seed randomizes every run."},
    ]
    if is_jax or model == "full":
        choices = [("reference_then_colabfold", "Reference MSA, then ColabFold server"),
                   ("local_mmseqs_gpu", "Reference MSA, then local MMseqs GPU"),
                   ("single_sequence", "Single sequence, no MSA")]
        extra_fields.extend([
            {"key": "msa_source", "type": "select", "label": "Protein MSA source", "choices": choices,
             "default": choices[0][0], "help": "Reuse matching reference alignments before generating missing per-chain MSAs."},
            {"key": "msa_max_depth", "label": "Maximum MSA rows per chain", "default": 1024, "min": 1, "max": 65536 if is_jax else 8192, "step": 128},
        ])
    else:
        profile["fixed_settings"] = ("ESMFold2 Fast uses single-sequence input without an MSA.",)
    values, settings, submitted = render_molecule_editor(
        f"{family}-cofold-{model}", ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles") if is_jax else (),
        supports_examples=True, supports_affinity=bool(profile.get("supports_affinity")),
        supports_templates=is_jax,
        modes=tuple(profile["modes"]), can_submit=bool(project_id),
        prediction_settings=prediction_settings_for(profile, extra_fields=tuple(extra_fields)),
    )
    if not submitted:
        return
    try:
        if values.get("input_error"):
            raise ValueError(values["input_error"])
        seed = int(settings.get("seed") or 0)
        count = int(settings.get("seed_count") or 1)
        seeds = _resolve_seed_values(",".join(str(seed + index if seed else 0) for index in range(count)), label=title)
        mode = settings["optimization_mode"]
        parameters = {
            "protein_sequences": values["protein_sequences"], "job_name": job_name or f"{family}-cofold",
            "optimization_mode": mode, "seeds": seeds, "num_diffusion_samples": int(settings["samples"]),
            "msa_source": settings.get("msa_source", "single_sequence"), "msa_max_depth": int(settings.get("msa_max_depth", 1024)),
            "gpu_device": "0,1" if mode == "big" and (not is_jax or profile["two_gpu_big"]) else "0",
        }
        if is_jax:
            parameters.update({key: values[key] for key in ("dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles")})
            parameters["structure_templates"] = values.get("structure_templates", [])
            parameters.update(model_name=model, num_recycles=int(settings["num_recycles"]),
                              cpu_cores=int(settings["cpu_cores"]), accept_alphafold3_terms=bool(settings.get("accept_terms")))
        else:
            parameters.update(model_variant=model, num_loops=int(settings["num_recycles"]), cpu_threads=int(settings["cpu_cores"]))
        _submit_cofolding_job(schema, workflow_id, parameters, project_id, settings, title)
    except (ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
        st.error(f"Could not start {title}: {exc}")


def _render_boltz_prediction_setup(profile: dict[str, Any], schema: dict[str, Any]) -> None:
    _render_cofolding_header("Boltz-2", "Predict protein, nucleic acid, and small molecule complexes · v2.2.1")
    with st.expander("About this workflow", expanded=False):
        st.write("Runs Boltz-2 locally through the Studio GPU scheduler.")
        st.code(profile["image"], language=None)
        st.caption("Available inference profiles: " + ", ".join(profile["modes"]))
    name_col, project_col = st.columns([1.25, 1.25])
    with name_col:
        job_name = st.text_input("Job name", value="boltz2-prediction", key="optimized-job-name-boltz2").strip()
    with project_col:
        project_id = _render_project_selector("optimized-prediction-boltz2")
    handled, task_count = render_batch_controls(
        "boltz2", profile, project_id, job_name, schema, search_jobs=search_jobs, submit=_submit_cofolding_job)
    if handled:
        return
    tasks, settings, submitted = [], {}, False
    tabs = st.tabs([f"Complex {i + 1}" for i in range(task_count)]) if task_count > 1 else [st.container()]
    for index, tab in enumerate(tabs, start=1):
        with tab:
            values, item_settings, item_submitted = render_molecule_editor(
                f"optimized-boltz2-{index}", tuple(profile["molecule_inputs"]), supports_cyclic=True,
                supports_modifications=True, supports_restraints=True, supports_msa_upload=True,
                supports_affinity=bool(profile.get("supports_affinity")),
                supports_affinity_control=bool(profile.get("supports_affinity")),
                boltz_settings=task_count == 1, supports_templates=task_count > 1,
                modes=tuple(profile["modes"]), can_submit=bool(project_id),
            )
            tasks.append(values)
            if task_count == 1:
                settings, submitted = item_settings, item_submitted
    if task_count > 1:
        _, settings, submitted = render_molecule_editor(
            "optimized-boltz2-batch-settings", (), boltz_settings=True, settings_only=True,
            modes=tuple(profile["modes"]), can_submit=bool(project_id))
    if not submitted:
        return
    try:
        for values in tasks:
            if values.pop("input_error", ""):
                raise ValueError("Complete every complex before submitting the batch.")
            values.pop("chain_catalog", None)
        engine_options = {
            key: settings[key] for key in (
                "sampling_steps", "max_parallel_samples", "step_scale", "output_format",
                "seed_count", "use_potentials", "method", "use_msa", "generate_ipsae", "boltz_version",
            )
        }
        engine_options["max_msa_seqs"] = int(settings.get("max_msa_seqs") or 8192)
        engine_options["use_msa_server"] = bool(settings["use_msa"])
        parameters = {
            "engine": "boltz2", "optimization_mode": settings["optimization_mode"],
            "job_name": job_name or "boltz2-prediction", "seed": int(settings.get("seed") or 0),
            "samples": int(settings["samples"]), "num_recycles": int(settings["num_recycles"]),
            "cpu_cores": 4, "gpu_device": "0,1" if settings["optimization_mode"] == "big" else "0",
            "tasks": tasks, "engine_options": engine_options,
        }
        payload = {
            "schema_version": schema["schema_version"], "workflow": "refolding.optimized_prediction",
            "parameters": parameters, "project_id": project_id,
            "priority": 1000 if settings.get("prioritize") else 0,
        }
        validation = _post_local_api("/api/v1/workflows/validate", payload)
        if not validation.get("valid"):
            st.error("The Boltz-2 inputs did not pass validation.")
            return
        result = _post_local_api("/api/v1/workflows", payload)
        mark_studio_job(result, project_id=project_id, owner_id=str(CURRENT_USER["id"]))
        st.session_state["last_submitted_job"] = result["job"]
        job = result.get("job", {})
        st.success(f"Boltz-2 job queued · Job {job.get('job_code') or job.get('run_id') or 'submitted'}")
    except (ValueError, TypeError, json.JSONDecodeError) as exc:
        st.error(f"Could not start the Boltz-2 prediction: {exc}")


def _render_cofolding_header(title: str, description: str) -> None:
    st.html(
        '<style>.block-container{padding-top:1.3rem;max-width:1500px}'
        '.cofolding-tool-header{text-align:center;margin:0 0 1.4rem}'
        '.cofolding-tool-header h1{font-size:2rem;margin:0 0 .4rem}'
        '.cofolding-tool-header p{color:#748087;font-size:.9rem;margin:.5rem 0}</style>'
        f'<div class="cofolding-tool-header"><h1>{escape(title)}</h1>'
        f'<p>{escape(description)}</p></div>'
    )


def _submit_cofolding_job(schema: dict, workflow_id: str, parameters: dict, project_id: str, settings: dict, label: str) -> None:
    payload = {
        "schema_version": schema["schema_version"], "workflow": workflow_id,
        "parameters": parameters, "project_id": project_id,
        "priority": 1000 if settings.get("prioritize") else 0,
    }
    validation = _post_local_api("/api/v1/workflows/validate", payload)
    if not validation.get("valid"):
        st.error(f"The {label} inputs did not pass validation.")
        return
    result = _post_local_api("/api/v1/workflows", payload)
    mark_studio_job(result, project_id=project_id, owner_id=str(CURRENT_USER["id"]))
    st.session_state["last_submitted_job"] = result["job"]
    job = result.get("job", {})
    st.success(f"{label} job queued · Job {job.get('job_code') or job.get('run_id') or 'submitted'}")


def _render_optimized_prediction_setup(engine_id: str) -> None:
    profile = optimized_prediction.ENGINE_PROFILES[engine_id]
    workflow_id = "refolding.optimized_prediction"
    schema = workflow_schema(workflow_id)
    if st.button("← Back to All tools", key="back-to-optimized-images"):
        st.session_state["page"] = "All tools"
        for key in ("selected_workflow", "selected_workflow_engine", "selected_workflow_group"):
            st.session_state.pop(key, None)
        st.rerun()
    if engine_id == "boltz2":
        _render_boltz_prediction_setup(profile, schema)
        return
    _render_cofolding_header(profile["label"], profile.get("input_summary", "Predict molecular structures."))
    with st.expander("About this workflow", expanded=False):
        st.write("Predict molecular structures through the Studio GPU scheduler.")
        st.code(profile["image"], language=None)
    name_col, project_col = st.columns([1.25, 1.25])
    with name_col:
        job_name = st.text_input("Job name", value=f"{engine_id}-prediction", key=f"optimized-job-name-{engine_id}").strip()
    with project_col:
        project_id = _render_project_selector(f"optimized-prediction-{engine_id}")
    handled, task_count = render_batch_controls(
        engine_id, profile, project_id, job_name, schema, search_jobs=search_jobs, submit=_submit_cofolding_job)
    if handled:
        return
    task_tabs = st.tabs([f"Complex {index + 1}" for index in range(task_count)]) if task_count > 1 else [st.container()]
    tasks = []
    settings = {}
    submitted = False
    config = prediction_settings_for(profile)
    for index, task_tab in enumerate(task_tabs, start=1):
        with task_tab:
            target_values = {}
            if profile["input_kind"] == "boltzgen":
                st.caption("Upload the target structure and enter one binder protein chain below. MSAs are reused from reference files or fetched from ColabFold; a target A3M upload is optional.")
                target_file = st.file_uploader("Target PDB structure", type=["pdb", "ent"], key=f"optimized-target-pdb-{engine_id}-{index}")
                target_chain_id = st.text_input("Target chain ID", value="A", max_chars=1, key=f"optimized-target-chain-{engine_id}-{index}")
                msa_file = st.file_uploader("Target chain A3M (optional)", type=["a3m", "txt"], key=f"optimized-target-msa-{engine_id}-{index}")
                target_values = {
                    "target_pdb_text": target_file.getvalue().decode("utf-8", errors="replace") if target_file else "",
                    "target_chain_id": target_chain_id,
                    "target_msa_a3m": msa_file.getvalue().decode("utf-8", errors="replace") if msa_file else "",
                }
            values, item_settings, item_submitted = render_molecule_editor(
                f"optimized-{engine_id}-{index}", tuple(profile.get("molecule_inputs", ())),
                supports_cyclic=bool(profile.get("supports_cyclic")),
                supports_modifications=bool(profile.get("supports_modifications")),
                supports_restraints=bool(profile.get("supports_restraints")),
                supports_msa_upload=bool(profile.get("supports_msa_upload")),
                supports_templates=bool(profile.get("supports_templates") and profile["input_kind"] != "boltzgen"),
                supports_examples=profile["input_kind"] != "boltzgen",
                supports_affinity=bool(profile.get("supports_affinity")),
                modes=tuple(profile["modes"]), can_submit=bool(project_id),
                prediction_settings=config if task_count == 1 else None,
            )
            if profile["input_kind"] == "boltzgen":
                values = {**target_values, "binder_sequence": values["protein_sequences"], "input_error": values.get("input_error", "")}
            tasks.append(values)
            if task_count == 1:
                settings, submitted = item_settings, item_submitted
    if task_count > 1:
        _, settings, submitted = render_molecule_editor(
            f"optimized-{engine_id}-batch-settings", (), modes=tuple(profile["modes"]), can_submit=bool(project_id),
            prediction_settings=config, settings_only=True,
        )
    if not submitted:
        return
    try:
        for task in tasks:
            if task.pop("input_error", ""):
                raise ValueError("Complete all molecule chains before submitting.")
            task.pop("chain_catalog", None)
        mode = str(settings.get("optimization_mode", profile["modes"][0]))
        parameters = {
            "engine": engine_id, "optimization_mode": mode, "job_name": job_name or f"{engine_id}-prediction",
            "seed": int(settings.get("seed") or 0), "samples": int(settings["samples"]),
            "num_recycles": int(settings.get("num_loops", settings.get("num_recycles", profile.get("default_recycles", 10)))),
            "cpu_cores": int(settings["cpu_cores"]),
            "gpu_device": "0,1" if profile.get("gpu_count", {}).get(mode, 1) == 2 else "0",
            "accept_alphafold3_terms": bool(settings.get("accept_terms")), "tasks": tasks,
            "engine_options": {field["key"]: settings[field["key"]] for field in profile.get("settings", ())},
        }
        _submit_cofolding_job(schema, workflow_id, parameters, project_id, settings, profile["label"])
    except (ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
        st.error(f"Could not start the prediction: {exc}")


def _render_workflow_setup() -> None:
    group_id = st.session_state.get("selected_workflow_group")
    if group_id in TOOL_CATALOG_GROUPS:
        group = TOOL_CATALOG_GROUPS[group_id]
        mode_label = st.selectbox(
            "Prediction mode",
            list(group["modes"]),
            key=f"tool-group-mode-{group_id}",
        )
        workflow_id = group["modes"][mode_label]
        st.session_state["selected_workflow"] = workflow_id
        engine_profile = group.get("engine_profiles", {}).get(mode_label)
        if engine_profile:
            st.session_state["selected_workflow_engine"] = engine_profile
        else:
            st.session_state.pop("selected_workflow_engine", None)
    else:
        st.session_state.pop("selected_workflow_group", None)
        workflow_id = st.session_state.get("selected_workflow")
    if not workflow_id:
        st.session_state["page"] = "All tools"
        st.rerun()

    # Keep previously selected ESMFold entries on the common cofolding form.
    # The old workflow remains available to queued/API consumers, but its
    # candidate-refolding input mask is no longer the user-facing ESMFold UI.
    if workflow_id == "refolding.esmfold_monomer":
        workflow_id = "refolding.optimized_prediction"
        st.session_state["selected_workflow"] = workflow_id
        st.session_state["selected_workflow_engine"] = "esmfold"

    # Keep the existing AlphaFold2 tile while routing its setup through the
    # direct PDB and chain-role form shared by the other cofolding workflows.
    if workflow_id == "refolding.af2_initial_guess":
        workflow_id = "refolding.af2_initial_guess_input"
        st.session_state["selected_workflow"] = workflow_id

    if workflow_id == "refolding.optimized_prediction":
        engine_id = st.session_state.get("selected_workflow_engine")
        if engine_id not in optimized_prediction.ENGINE_PROFILES:
            engine_id = st.selectbox(
                "Prediction engine",
                list(optimized_prediction.ENGINE_PROFILES),
                format_func=lambda key: optimized_prediction.ENGINE_PROFILES[key]["label"],
                key="optimized-prediction-engine-fallback",
            )
        _render_optimized_prediction_setup(str(engine_id))
        return

    schema = workflow_schema(workflow_id)
    from mn_molecular_design_studio import inverse_folding_input
    if workflow_id in inverse_folding_input.WORKFLOWS.values():
        from mn_molecular_design_studio import inverse_folding_ui
        importlib.reload(inverse_folding_ui)
        engine = next(k for k, v in inverse_folding_input.WORKFLOWS.items() if v == workflow_id)
        inverse_folding_ui.render_setup(engine, schema, owner_id=str(CURRENT_USER["id"]),
            project_selector=_render_project_selector, submit=_submit_cofolding_job, search_jobs=search_jobs)
        return
    if workflow_id in {"sequence.mpnn", "sequence.ligandmpnn", "sequence.foundry_mpnn"}:
        from mn_molecular_design_studio import mpnn_ui
        importlib.reload(mpnn_ui)
        mpnn_ui.render_setup(workflow_schema("sequence.mpnn"), owner_id=str(CURRENT_USER["id"]),
            project_selector=_render_project_selector, submit=_submit_cofolding_job, search_jobs=search_jobs)
        return
    if workflow_id in generator_input.WORKFLOWS.values():
        from mn_molecular_design_studio import generator_ui
        importlib.reload(generator_ui)
        generator_ui.render_setup(schema,owner_id=str(CURRENT_USER['id']),project_selector=_render_project_selector,submit=_submit_cofolding_job)
        return
    if workflow_id in {"design.rfdiffusion3", "design.rfdiffusion3_foundry"}:
        from mn_molecular_design_studio import rfdiffusion3_ui
        importlib.reload(rfdiffusion3_ui)
        rfdiffusion3_ui.render_setup(workflow_schema("design.rfdiffusion3"), owner_id=str(CURRENT_USER["id"]),
                     header=_render_cofolding_header, project_selector=_render_project_selector,
                     submit=_submit_cofolding_job)
        return
    if workflow_id == "design.rfantibody":
        from mn_molecular_design_studio import rfantibody_ui
        importlib.reload(rfantibody_ui)
        rfantibody_ui.render_setup(schema, owner_id=str(CURRENT_USER["id"]),
                                   project_selector=_render_project_selector, submit=_submit_cofolding_job)
        return
    if workflow_id == "design.rfpeptides":
        from mn_molecular_design_studio import rfpeptides_ui
        importlib.reload(rfpeptides_ui)
        rfpeptides_ui.render_setup(schema, owner_id=str(CURRENT_USER["id"]),
                                   project_selector=_render_project_selector, submit=_submit_cofolding_job)
        return
    if workflow_id in {"design.proteina_complexa", "design.proteina_complexa_generate"}:
        from mn_molecular_design_studio import complexa_ui
        importlib.reload(complexa_ui)
        complexa_ui.render_setup(schema,owner_id=str(CURRENT_USER["id"]),project_selector=_render_project_selector,submit=_submit_cofolding_job)
        return
    if workflow_id == "design.esmfold2_binder":
        from mn_molecular_design_studio import esmfold2_binder_ui
        importlib.reload(esmfold2_binder_ui)
        esmfold2_binder_ui.render_setup(schema, project_selector=_render_project_selector, submit=_submit_cofolding_job)
        return
    if workflow_id in {"design.bindcraft2", "design.bindcraft2_generator"}:
        _render_bindcraft2_setup(schema, generator_only=workflow_id.endswith("_generator"))
        return
    if workflow_id in {"design.boltzgen", "design.boltzgen_generator"}:
        _render_boltzgen_setup(schema, generator_only=workflow_id.endswith("_generator"))
        return
    if workflow_id == "refolding.esmfold2_cofolding":
        _render_esmfold2_cofolding_setup(schema)
        return
    if workflow_id == "refolding.af3_jax_cofolding":
        _render_af3_jax_cofolding_setup(schema)
        return
    if workflow_id == "refolding.af2_initial_guess_input":
        _render_af2_initial_guess_setup(schema)
        return

    if workflow_id in {
        "detection.pesto",
        "detection.scannet",
        "detection.masif_seed",
        "detection.surf2spot",
        "detection.surface_analysis",
    }:
        from mn_molecular_design_studio import surface_ui
        importlib.reload(surface_ui)
        surface_ui.render_setup(
            schema,
            owner_id=str(CURRENT_USER["id"]),
            project_selector=_render_project_selector,
            submit=_submit_cofolding_job,
        )
        return

    base_title = workflow_id.split(".", 1)[-1].replace("_", " ").title()
    title = WORKFLOW_TITLES.get(base_title, base_title)

    if st.button("← Back to All tools", key="back-to-all-tools"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.session_state.pop("selected_workflow_group", None)
        st.rerun()

    st.caption(schema["category"].replace("_", " ").upper())
    st.title(title)
    st.write(schema["description"])
    st.caption("This setup launches the existing workflow with the inputs you provide.")
    project_id = _render_project_selector(f"workflow-{workflow_id}")

    analysis_source_workflow = workflow_id in ANALYSIS_SOURCE_WORKFLOWS
    source_jobs: list[dict] = []
    selected_source_job: str | None = None
    selected_source_pae_count: int | None = None
    if analysis_source_workflow:
        st.subheader("Existing results")
        if not project_id:
            st.info("Select a project first. Only completed results from that project will be available here.")
        else:
            try:
                source_jobs = _analysis_source_job_options(project_id)
            except Exception as exc:
                st.error(f"Could not load existing candidate results: {exc}")
        if project_id and source_jobs:
            selected_source_job = st.selectbox(
                "Results to analyze",
                [item["reference"] for item in source_jobs],
                format_func=lambda reference: next(
                    (item["label"] for item in source_jobs if item["reference"] == reference), reference
                ),
                key=f"analysis-source-job-{workflow_id}-{project_id}",
                help="The selected job is read as input. This tool creates a separate analysis result and leaves its original ranking and files unchanged.",
            )
        elif project_id:
            st.info("No completed jobs with saved complex structures are available in this project yet.")

        if workflow_id == "analysis.ipsae" and selected_source_job:
            selected_source = next(
                (item for item in source_jobs if item["reference"] == selected_source_job),
                None,
            )
            if selected_source:
                pae_count, complex_count = _source_job_pae_availability(selected_source["row"])
                selected_source_pae_count = pae_count
                if pae_count:
                    st.caption(f"Pairwise PAE matrix files found for {pae_count} of {complex_count} complex candidate(s).")
                else:
                    st.warning(
                        "No pairwise PAE matrix file was found for this job. ipSAE needs the PAE matrix; "
                        "scalar i_pAE or i_pDAE metrics are not sufficient."
                    )

    parameters = schema.get("parameters", [])
    if analysis_source_workflow:
        parameters = [
            parameter
            for parameter in parameters
            if parameter.get("name") not in {"source_run_dir", "candidates_jsonl", "tool"}
        ]
    if workflow_id == "analysis.rosetta":
        parameters = [parameter for parameter in parameters if parameter.get("name") != "nprocs"]
    required = [parameter for parameter in parameters if parameter.get("required")]
    optional = [parameter for parameter in parameters if not parameter.get("required")]
    options = schema.get("options_schema", [])

    with st.form(key=f"workflow-setup-form-{workflow_id}"):
        st.subheader("Required inputs" if required else "Workflow setup")
        if required:
            raw_values = _render_parameter_fields(required, workflow_id=workflow_id, group="required")
        else:
            st.caption("This workflow can start with its default settings.")
            raw_values = {}

        if optional or options:
            with st.expander(f"Advanced settings · {len(optional) + len(options)} options"):
                if optional:
                    st.markdown("#### Workflow options")
                    raw_values.update(_render_parameter_fields(optional, workflow_id=workflow_id, group="advanced"))
                if options:
                    st.markdown("#### Analysis options")
                    raw_options = _render_parameter_fields(options, workflow_id=workflow_id, group="options")
                else:
                    raw_options = {}
        else:
            raw_options = {}

        prioritize_job = st.checkbox(
            "Prioritize in scheduler queue",
            key=f"workflow-prioritize-queue-{workflow_id}",
            help="Moves this job ahead of normal queued work when resources become available. Running jobs are never interrupted.",
        )

        submitted = st.form_submit_button(
            "Run workflow",
            type="primary",
            use_container_width=True,
            disabled=not project_id
            or (analysis_source_workflow and not selected_source_job)
            or (workflow_id == "analysis.ipsae" and selected_source_pae_count == 0),
        )

    if submitted:
        try:
            values = {
                parameter["name"]: _parameter_value(raw_values.get(parameter["name"]), parameter)
                for parameter in parameters
            }
            if options:
                values["options"] = {
                    parameter["name"]: _parameter_value(raw_options.get(parameter["name"]), parameter)
                    for parameter in options
                }
            payload = {
                "schema_version": schema["schema_version"],
                "workflow": workflow_id,
                "parameters": values,
                "project_id": project_id,
                "priority": 1000 if prioritize_job else 0,
            }
            if analysis_source_workflow:
                payload["source_job"] = selected_source_job
            validation = _post_local_api("/api/v1/workflows/validate", payload)
            if not validation.get("valid"):
                st.error("The workflow inputs did not pass validation.")
            else:
                result = _post_local_api("/api/v1/workflows", payload)
                mark_studio_job(result, project_id=project_id, owner_id=str(CURRENT_USER["id"]))
                job_reference = result.get("job", {})
                st.success(f"Workflow queued · Job {job_reference.get('job_code') or job_reference.get('run_id') or 'submitted'}")
                st.session_state["last_submitted_job"] = result["job"]
        except (ValueError, TypeError, json.JSONDecodeError) as exc:
            st.error(f"Could not start workflow: {exc}")


def _render_pipelines() -> None:
    st.title("Pipelines")
    view = st.segmented_control("Pipeline workspace", ["Node canvas", "Compare design engines"], default="Node canvas", key="pipeline-workspace-view", label_visibility="collapsed")
    if view == "Compare design engines":
        _render_multi_engine_campaign()
        return
    st.caption("Connect tools into a reusable experiment. Move blocks freely, branch outputs, and keep every run as a separate snapshot.")
    from mn_molecular_design_studio import pipeline_node_ui, pipeline_graph_ui
    importlib.reload(pipeline_node_ui)
    importlib.reload(pipeline_graph_ui)
    pipeline_graph_ui.render_pipeline_builder(owner_id=str(CURRENT_USER["id"]), project_selector=_render_project_selector, submit_api=_post_local_api, search_jobs=search_jobs)


def _render_multi_engine_campaign() -> None:
    st.caption("Send one target through several design engines and collect their candidates in one campaign.")
    flow_cols = st.columns([1, 0.2, 1, 0.2, 1])
    flow_cols[0].markdown("**1 · Target**")
    flow_cols[1].markdown("→")
    flow_cols[2].markdown("**2 · Design engines**")
    flow_cols[3].markdown("→")
    flow_cols[4].markdown("**3 · Shared results**")

    name_col, project_col = st.columns([1.25, 1.25])
    with name_col:
        job_name = st.text_input("Pipeline name", value="multi-engine-design", key="pipeline-job-name").strip()
    with project_col:
        project_id = _render_project_selector("pipeline")
    st.subheader("Target")
    st.caption("Upload a PDB file or fetch a structure by PDB accession. Choose the target chains and optional hotspots below.")
    target_source = render_target_source(
        key_prefix="pipeline",
        label="Target",
        example_code="7JZU",
        upload_label="Drop a PDB file here or browse",
        upload_help="PDB files only. CIF conversion is not part of the campaign workflow yet.",
        upload_in_expander=True,
        owner_id=str(CURRENT_USER["id"]),
    )
    target_input = render_target_editor(target_source, key_prefix="pipeline", label="Target")

    st.subheader("Design engines")
    engine_options = [engine for engine in ENGINE_ORDER if engine != "template_redesign"]
    engines_key = "pipeline-selected-engines"
    if engines_key not in st.session_state:
        st.session_state[engines_key] = [
            engine for engine in ("rfdiffusion_classic", "bindcraft", "rfdiffusion3_foundry")
            if engine in engine_options
        ]
    else:
        st.session_state[engines_key] = [
            engine for engine in st.session_state[engines_key] if engine in engine_options
        ]
    selected_engines = st.multiselect(
        "Run the target through",
        engine_options,
        format_func=lambda engine: ENGINE_LABELS[engine],
        key=engines_key,
        help="Each selected engine runs as part of one campaign, and its candidates are gathered in the campaign results.",
    )
    if selected_engines:
        st.caption("Selected: " + " · ".join(ENGINE_LABELS[engine] for engine in selected_engines))

    setting_cols = st.columns(4)
    with setting_cols[0]:
        binder_length = st.text_input("Binder length", value="60-100", key="pipeline-binder-length").strip()
    with setting_cols[1]:
        design_attempts = st.number_input("Design attempts per engine", min_value=1, max_value=1000, value=2, step=1, key="pipeline-design-attempts")
    with setting_cols[2]:
        sequences_per_backbone = st.number_input("Sequences per backbone", min_value=1, max_value=1000, value=2, step=1, key="pipeline-sequences-per-backbone")
    with setting_cols[3]:
        random_seed = st.number_input("Random seed", min_value=0, max_value=2_147_483_647, value=0, step=1, key="pipeline-random-seed")
    gpu_device = st.text_input("GPU device", value="0", key="pipeline-gpu-device", help="Enter a GPU ID such as 0, or cpu where supported.").strip() or "0"
    prioritize_pipeline = st.checkbox(
        "Prioritize in scheduler queue",
        key="pipeline-prioritize-queue",
        help="Moves this campaign ahead of normal queued work when resources become available. Running jobs are never interrupted.",
    )

    if target_input.hotspots:
        st.caption(f"Hotspots: {target_input.hotspots}")
    st.info("This shared campaign currently applies hotspots across engines. Coldspot selection is available in the standalone BindCraft2 workflow.")
    if not selected_engines:
        st.warning("Select at least one design engine.")
    if not project_id:
        st.info("Select or create a project before submitting this pipeline.")
    if not target_input.source_path:
        st.info("Load a target structure to configure the pipeline.")

    with st.expander("Review pipeline", expanded=False):
        st.json(
            {
                "target": target_input.source_path.name if target_input.source_path else None,
                "target_chains": target_input.chains,
                "hotspots": target_input.hotspots,
                "engines": [ENGINE_LABELS[engine] for engine in selected_engines],
                "binder_length": binder_length,
                "design_attempts_per_engine": int(design_attempts),
                "sequences_per_backbone": int(sequences_per_backbone),
                "gpu_device": gpu_device,
            }
        )

    ready = bool(
        job_name
        and project_id
        and target_input.ready
        and target_input.design_path
        and selected_engines
        and binder_length
        and target_input.valid_trim_ranges
    )
    if st.button("Run pipeline", type="primary", use_container_width=True, disabled=not ready, key="pipeline-submit"):
        attempts = int(design_attempts)
        engine_configs = {}
        if "rfdiffusion_classic" in selected_engines:
            engine_configs["rfdiffusion_classic"] = {
                "timesteps": 50,
                "model_weights": "Complex_base",
                "sequence_design_method": "protein_mpnn",
                "execution_backend": "docker",
                "analysis_keep_top_n": attempts,
            }
        if "bindcraft" in selected_engines:
            engine_configs["bindcraft"] = {
                "number_of_final_designs": min(2, attempts),
                "advanced_settings_file": "default_4stage_multimer_mpnn.json",
                "filter_settings": "no_filters.json",
            }
        if "bindcraft2" in selected_engines:
            engine_configs["bindcraft2"] = {
                "number_of_final_designs": attempts,
                "max_trajectories": attempts,
                "cpu_cores": 1,
            }
        if "rfdiffusion3_foundry" in selected_engines:
            engine_configs["rfdiffusion3_foundry"] = {"timesteps": 50, "mpnn_model_type": "protein_mpnn"}
        if "boltzgen" in selected_engines:
            engine_configs["boltzgen"] = {"budget": min(2, attempts), "sampling_steps": 20}
        request_payload = {
            "schema_version": 1,
            "campaign_name": job_name,
            "project_id": project_id,
            "target_pdb": str(target_input.design_path),
            "target_chains": target_input.chains,
            "binder_length": binder_length,
            "hotspots": target_input.hotspots,
            "design_attempts": attempts,
            "sequences_per_backbone": int(sequences_per_backbone),
            "random_seed": int(random_seed),
            "engines": selected_engines,
            "engine_configs": engine_configs,
            "workflow_recipe": "vanilla",
            "gpu_device": gpu_device,
            "priority": 1000 if prioritize_pipeline else 0,
            "survivors_per_engine": 1_000_000,
            "passing_only": False,
            "keep_best_failed": True,
            "continue_after_failure": True,
            "common_validation": {"enabled": False},
            "sequence_refinement": {"mode": "none"},
            "template_redesign": {"enabled": False},
            "evaluation": {"mode": "none"},
        }
        try:
            app_home().mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(prefix="pipeline-request-", dir=app_home()) as request_dir:
                request_path = Path(request_dir) / "campaign.json"
                request_path.write_text(json.dumps(request_payload, indent=2), encoding="utf-8")
                result = submit_design_campaign_request(request_path)
            mark_studio_job(result, project_id=project_id, owner_id=str(CURRENT_USER["id"]))
            job_reference = result.get("job_code") or result.get("run_id")
            st.session_state["pipeline-last-submission"] = job_reference
            st.success(f"Pipeline queued · Job {job_reference}")
            st.caption("Track engine progress and open the combined candidate set in My results and Database.")
        except (AutomationError, OSError, TypeError, ValueError) as exc:
            st.error(f"Could not start the pipeline: {exc}")


def _parse_job_datetime(value: object) -> datetime | None:
    if not isinstance(value, str) or not value.strip():
        return None
    try:
        parsed = datetime.fromisoformat(value.strip().replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed


def _job_submitted_label(value: object) -> str:
    parsed = _parse_job_datetime(value)
    return parsed.astimezone().strftime("%Y-%m-%d %I:%M %p") if parsed else "—"


def _job_runtime_timestamps(row: dict) -> dict[str, object]:
    timestamps = {
        key: row.get(key)
        for key in ("worker_started_at", "queue_wait_started_at", "completed_at")
    }
    if all(value for value in timestamps.values()):
        return timestamps

    raw_run_dir = row.get("run_dir")
    if not isinstance(raw_run_dir, str) or not raw_run_dir.strip():
        return timestamps
    try:
        root = runs_root().expanduser().resolve()
        run_dir = Path(raw_run_dir).expanduser().resolve()
        run_dir.relative_to(root)
        metadata_path = (run_dir / "metadata.json").resolve(strict=True)
        metadata_path.relative_to(root)
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return timestamps
    if isinstance(metadata, dict):
        for key in timestamps:
            timestamps[key] = timestamps[key] or metadata.get(key)
    return timestamps


def _duration_label(seconds: float) -> str:
    total = max(0, int(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes:02d}m"
    if minutes:
        return f"{minutes}m {seconds:02d}s"
    return f"{seconds}s"


def _job_elapsed_label(row: dict) -> str:
    status = str(row.get("status") or "").casefold()
    timestamps = _job_runtime_timestamps(row)
    now = datetime.now(timezone.utc)
    if status == "queued":
        queued_at = _parse_job_datetime(timestamps.get("queue_wait_started_at")) or _parse_job_datetime(
            row.get("created_at")
        )
        return f"Queued · {_duration_label((now - queued_at).total_seconds())}" if queued_at else "Queued"

    started_at = _parse_job_datetime(timestamps.get("worker_started_at"))
    if status in {"running", "preparing"}:
        if not started_at:
            return "Starting"
        return f"Running · {_duration_label((now - started_at).total_seconds())}"
    if not started_at:
        return "—"

    completed_at = _parse_job_datetime(timestamps.get("completed_at")) or _parse_job_datetime(row.get("updated_at"))
    if not completed_at:
        return "—"
    return _duration_label((completed_at - started_at).total_seconds())


def _jobs_csv(rows: list[dict]) -> bytes:
    stream = io.StringIO(newline="")
    fields = ("Job name", "Project", "Tool", "Pipeline", "Job code", "Status", "Submitted", "Time elapsed", "Stage", "Run ID")
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    project_names = {str(project["id"]): str(project["name"]) for project in _studio_projects()}
    for row in rows:
        writer.writerow(
            {
                "Job name": row.get("campaign_name") or row.get("tool") or row.get("job_type") or "Job",
                "Project": project_names.get(str(row.get("project_id") or ""), ""),
                "Tool": display_tool_label(row),
                "Pipeline": row.get("pipeline_name") or "—",
                "Job code": row.get("job_code", ""),
                "Status": row.get("status", "unknown").title(),
                "Submitted": _job_submitted_label(row.get("created_at")),
                "Time elapsed": _job_elapsed_label(row),
                "Stage": row.get("current_phase") or row.get("current_engine") or "—",
                "Run ID": row.get("run_id", ""),
            }
        )
    return stream.getvalue().encode("utf-8-sig")


def _render_runs() -> None:
    if cluster_config.enabled():
        from mn_molecular_design_studio.cluster.ui import render_runs
        render_runs(CURRENT_USER, str(st.session_state.get("studio-session-token") or "development-mode"), search_jobs)
        return
    st.title("My results")
    try:
        facets = job_list_facets()
    except Exception as exc:
        st.error(f"Could not load jobs: {exc}")
        return

    refresh_col, search_col, project_col, filter_col = st.columns([0.06, 0.42, 0.28, 0.24], vertical_alignment="bottom")
    if refresh_col.button("", icon=":material/refresh:", help="Refresh job list", key="results-refresh"):
        st.rerun()
    query = search_col.text_input(
        "Search jobs",
        placeholder="Search job name, workflow, or job ID…",
        label_visibility="collapsed",
        key="results-search",
    ).strip().casefold()
    projects = _studio_projects()
    project_options = ["", "__none__", *[str(project["id"]) for project in projects]]
    if st.session_state.get("results-project-filter") not in project_options:
        st.session_state["results-project-filter"] = ""
    selected_project = project_col.selectbox(
        "Project",
        project_options,
        format_func=lambda value: (
            "All projects" if value == "" else "No project" if value == "__none__" else next(
                (str(project["name"]) for project in projects if project["id"] == value), "Unknown project"
            )
        ),
        key="results-project-filter",
        label_visibility="collapsed",
    )
    with filter_col.popover("Add filter", icon=":material/filter_alt:", use_container_width=True):
        chosen_statuses = st.multiselect(
            "Status", facets["statuses"], key="results-status-filter"
        )
        chosen_tools = st.multiselect("Tool", facets["tools"], key="results-tool-filter")
    selected_project_id = selected_project if selected_project else None

    page_size = st.selectbox("Rows per page", [10, 25, 50], index=0, key="results-page-size")
    filters_key = (query, tuple(sorted(chosen_statuses)), tuple(sorted(chosen_tools)), selected_project_id, page_size)
    if st.session_state.get("results-filter-state") != filters_key:
        st.session_state["results-filter-state"] = filters_key
        st.session_state["results-page"] = 0
        st.session_state.pop("results-export-data", None)
        st.session_state.pop("results-export-count", None)
        st.session_state["results-selected-jobs"] = {}
        for state_key in list(st.session_state):
            if str(state_key).startswith("results-row-select-"):
                st.session_state.pop(state_key, None)
    page_number = max(0, int(st.session_state.get("results-page", 0)))
    offset = page_number * page_size
    try:
        result_page = search_jobs(
            query=query,
            statuses=chosen_statuses,
            tools=chosen_tools,
            project_id=selected_project_id,
            limit=page_size,
            offset=offset,
        )
    except Exception as exc:
        st.error(f"Could not search jobs: {exc}")
        return
    page_rows = result_page["items"]
    total_jobs = result_page["total"]
    page_count = max(1, (total_jobs + page_size - 1) // page_size)
    if page_number >= page_count:
        page_number = page_count - 1
        st.session_state["results-page"] = page_number
        offset = page_number * page_size
        result_page = search_jobs(
            query=query,
            statuses=chosen_statuses,
            tools=chosen_tools,
            project_id=selected_project_id,
            limit=page_size,
            offset=offset,
        )
        page_rows = result_page["items"]

    action_col, columns_col, count_col = st.columns([0.25, 0.25, 0.5], vertical_alignment="center")
    if action_col.button(
        "Prepare CSV export",
        icon=":material/download:",
        disabled=not total_jobs,
        use_container_width=True,
        key="results-export-prepare",
    ):
        export_rows: list[dict] = []
        export_offset = 0
        while export_offset < total_jobs:
            export_page = search_jobs(
                query=query,
                statuses=chosen_statuses,
                tools=chosen_tools,
                project_id=selected_project_id,
                limit=1000,
                offset=export_offset,
            )
            if not export_page["items"]:
                break
            export_rows.extend(export_page["items"])
            export_offset += len(export_page["items"])
        st.session_state["results-export-data"] = _jobs_csv(export_rows)
        st.session_state["results-export-count"] = len(export_rows)
    export_data = st.session_state.get("results-export-data")
    if export_data:
        action_col.download_button(
            f"Download CSV ({st.session_state.get('results-export-count', 0)})",
            data=export_data,
            file_name="molecular-design-jobs.csv",
            mime="text/csv",
            use_container_width=True,
            key="results-export-download",
        )
    all_columns = ["Job name", "Project", "Tool", "Pipeline", "#", "Status", "Submitted", "Time elapsed", "Stage", "Run ID"]
    default_columns = ["Job name", "Project", "Tool", "Pipeline", "#", "Status", "Submitted", "Time elapsed"]
    stored_columns = st.session_state.get("results-visible-columns")
    if isinstance(stored_columns, list) and "Type" in stored_columns:
        migrated_columns = ["Tool" if column == "Type" else column for column in stored_columns if column != "Type" or "Tool" not in stored_columns]
        if "Pipeline" not in migrated_columns:
            migrated_columns.insert(migrated_columns.index("Tool") + 1 if "Tool" in migrated_columns else 0, "Pipeline")
        st.session_state["results-visible-columns"] = [column for column in migrated_columns if column in all_columns]
    if not st.session_state.get("results-elapsed-column-default-applied"):
        existing_columns = st.session_state.get("results-visible-columns")
        if isinstance(existing_columns, list) and "Time elapsed" not in existing_columns:
            st.session_state["results-visible-columns"] = [*existing_columns, "Time elapsed"]
        st.session_state["results-elapsed-column-default-applied"] = True
    with columns_col.popover("Columns", icon=":material/view_column:", use_container_width=True):
        visible_columns = st.multiselect(
            "Choose columns",
            all_columns,
            default=default_columns,
            key="results-visible-columns",
            label_visibility="collapsed",
        )
    count_col.caption(f"{total_jobs} job{'s' if total_jobs != 1 else ''}")

    if not total_jobs:
        if facets["tools"]:
            st.info("No jobs match your search and filters.")
        else:
            st.info("No jobs yet. Choose a workflow from All tools to get started.")
        return
    if not visible_columns:
        st.warning("Choose at least one column to display jobs.")
        return

    start_index = offset
    display_rows = [
        {
            "Job name": row.get("campaign_name") or row.get("tool") or row.get("job_type") or "Job",
            "Project": next((str(project["name"]) for project in projects if project.get("id") == row.get("project_id")), "—"),
            "Tool": display_tool_label(row),
            "Pipeline": row.get("pipeline_name") or "—",
            "#": row.get("job_code", ""),
            "Status": str(row.get("status", "unknown")).title(),
            "Submitted": _job_submitted_label(row.get("created_at")),
            "Time elapsed": _job_elapsed_label(row),
            "Stage": row.get("current_phase") or row.get("current_engine") or "—",
            "Run ID": row.get("run_id", ""),
        }
        for row in page_rows
    ]
    row_weights = [0.24]
    row_weights.extend(
        {
            "Job name": 2.0,
            "Project": 1.1,
            "Tool": 1.5,
            "Pipeline": 1.4,
            "#": 0.45,
            "Status": 0.8,
            "Submitted": 1.2,
            "Time elapsed": 0.9,
            "Stage": 1.0,
            "Run ID": 1.45,
        }.get(column, 1.0)
        for column in visible_columns
    )
    row_weights.append(0.5)

    header_cols = st.columns(row_weights, vertical_alignment="center")
    header_cols[0].markdown("**Select**")
    for index, column_name in enumerate(visible_columns, start=1):
        header_cols[index].markdown(f"**{escape(column_name)}**")
    header_cols[-1].markdown("**View**")

    selected_job_map = st.session_state.setdefault("results-selected-jobs", {})
    for row, display_row in zip(page_rows, display_rows):
        reference = f"{row['task_group']}/{row['run_id']}"
        row_key = "results-row-select-" + hashlib.sha1(reference.encode("utf-8")).hexdigest()[:12]
        with st.container(border=True):
            row_cols = st.columns(row_weights, vertical_alignment="center")
            is_selected = row_cols[0].checkbox(
                "Select job",
                value=reference in selected_job_map,
                key=row_key,
                label_visibility="collapsed",
            )
            if is_selected:
                selected_job_map[reference] = row
            else:
                selected_job_map.pop(reference, None)

            for index, column_name in enumerate(visible_columns, start=1):
                value = display_row.get(column_name, "—")
                if column_name == "Status":
                    status_text = escape(str(value))
                    status_color = "#1d6848" if str(value).casefold() in {"completed", "succeeded"} else "#48554e"
                    row_cols[index].html(
                        f'<span style="color:{status_color};font-weight:600">{status_text}</span>'
                    )
                else:
                    row_cols[index].write(value)

            view_href = f"/?job_ref={quote(str(row['run_id']), safe='')}"
            row_cols[-1].html(
                f'<a href="{escape(view_href, quote=True)}" target="_self" '
                'style="display:inline-block;padding:0.25rem 0.55rem;border-radius:0.4rem;'
                'background:#e8f1ec;color:#1f6849;font-weight:600;text-decoration:none">View</a>'
            )
    st.session_state["results-selected-jobs"] = selected_job_map

    previous_col, range_col, next_col = st.columns([0.2, 0.6, 0.2], vertical_alignment="center")
    if previous_col.button("← Previous", disabled=page_number == 0, key="results-previous"):
        st.session_state["results-page"] = page_number - 1
        st.rerun()
    end_index = min(start_index + len(page_rows), total_jobs)
    range_col.markdown(
        f"<div style='text-align:center'>{start_index + 1}–{end_index} of {total_jobs}</div>",
        unsafe_allow_html=True,
    )
    if next_col.button("Next →", disabled=page_number + 1 >= page_count, key="results-next"):
        st.session_state["results-page"] = page_number + 1
        st.rerun()

    selected_rows = list(selected_job_map.values())
    active_selected = [row for row in selected_rows if str(row.get("status", "")).lower() in ACTIVE_STATUSES]
    if selected_rows:
        st.caption(f"{len(selected_rows)} job{'s' if len(selected_rows) != 1 else ''} selected")
        if active_selected:
            if st.button("Stop selected jobs", icon=":material/stop_circle:", key="results-stop-selected"):
                st.session_state["results-pending-stop"] = [
                    f"{row['task_group']}/{row['run_id']}" for row in active_selected
                ]
                st.rerun()
    pending_stop = st.session_state.get("results-pending-stop", [])
    if pending_stop:
        st.warning(f"Stop {len(pending_stop)} active job{'s' if len(pending_stop) != 1 else ''}? This will request cancellation.")
        confirm_col, cancel_col = st.columns(2)
        if confirm_col.button("Confirm stop", type="primary", key="results-confirm-stop"):
            errors = []
            stopped = 0
            for reference in pending_stop:
                try:
                    _cancel_current_user_job(reference, reason="Cancelled from Molecular Design Studio results")
                    stopped += 1
                except AutomationError as exc:
                    errors.append(str(exc))
            st.session_state.pop("results-pending-stop", None)
            if stopped:
                st.success(f"Cancellation requested for {stopped} job{'s' if stopped != 1 else ''}.")
            for error in errors:
                st.error(error)
            st.rerun()
        if cancel_col.button("Keep running", key="results-cancel-stop"):
            st.session_state.pop("results-pending-stop", None)
            st.rerun()

    if len(selected_rows) == 1:
        selected_row = selected_rows[0]
        open_col, stop_hint_col = st.columns([0.25, 0.75], vertical_alignment="center")
        if open_col.button("Open job results", icon=":material/open_in_new:", type="primary", key="results-open-job"):
            st.session_state["selected-job-reference"] = f"{selected_row['task_group']}/{selected_row['run_id']}"
            st.session_state["page"] = "Job details"
            st.rerun()
        stop_hint_col.caption("Open the saved outputs, files, inputs, and logs for this job.")


def _job_artifact_paths(details: dict) -> list[tuple[dict, Path]]:
    """Resolve catalog artifact links while keeping downloads inside this run."""
    _require_current_user_job_access(details["job"])
    run_dir = Path(details["job"]["run_dir"]).expanduser().resolve()
    root = runs_root().expanduser().resolve()
    artifacts: list[tuple[dict, Path]] = []
    try:
        catalog_rows = job_artifacts(str(details["job"]["task_group"]), str(details["job"]["run_id"]))
    except Exception:
        catalog_rows = []
    for item in catalog_rows:
        try:
            path = (root / str(item["relative_path"])).resolve(strict=True)
            path.relative_to(run_dir)
            if path.is_file():
                artifacts.append((item, path))
        except (KeyError, OSError, ValueError):
            continue
    return artifacts


_AA_TYPE_STYLE = {
    "hydrophobic": ("#f7d794", "#392b00"),
    "polar": ("#b7e4c7", "#17422b"),
    "basic": ("#b8d5f5", "#163257"),
    "acidic": ("#f3b7b7", "#5a1111"),
    "special": ("#d8c2f0", "#392257"),
    "other": ("#d9dee5", "#2d3748"),
}
_AA_TYPE_GROUPS = {
    "hydrophobic": set("AVILMFWY"),
    "polar": set("STNQ"),
    "basic": set("KRH"),
    "acidic": set("DE"),
    "special": set("CGP"),
}
def _sequence_column_name(columns) -> object | None:
    normalized = {
        re.sub(r"[^a-z0-9]+", "", str(column).casefold()): column for column in columns
    }
    return next(
        (normalized[key] for key in ("bindersequence", "binderseq", "sequence", "seq") if key in normalized),
        None,
    )


def _bindcraft2_resubmit_parameters(details: dict) -> tuple[dict | None, str | None]:
    job = details.get("job") if isinstance(details.get("job"), dict) else {}
    metadata = details.get("metadata") if isinstance(details.get("metadata"), dict) else {}
    input_payload = details.get("input") if isinstance(details.get("input"), dict) else {}
    inputs = input_payload.get("inputs") if isinstance(input_payload.get("inputs"), dict) else {}
    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    if metadata.get("demo_fixture"):
        return None, "This imported demo does not contain the original target FASTA or structure file."
    if not str(job.get("job_type") or input_payload.get("job_type") or "").startswith("bindcraft2"):
        return None, "Resubmission is available for reproducible BindCraft2 design jobs."

    run_dir = Path(str(job.get("run_dir") or "")).expanduser().resolve()
    if not run_dir.is_dir():
        return None, "The original job folder is unavailable."
    raw_dir = run_dir / "artifacts" / "raw" / "bindcraft2"
    try:
        settings_value = json.loads((raw_dir / "settings.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        settings_value = {}
    settings = settings_value if isinstance(settings_value, dict) else {}
    targets = settings.get("targets") if isinstance(settings.get("targets"), list) else []
    primary = targets[0] if targets and isinstance(targets[0], dict) else {}
    second = targets[1] if len(targets) > 1 and isinstance(targets[1], dict) else {}
    target_type = str(params.get("target_type") or inputs.get("target_type") or "structure").casefold()

    def stored_file(reference: object, fallback: str) -> Path | None:
        candidate = Path(str(reference or fallback)).expanduser()
        if not candidate.is_absolute():
            candidate = run_dir / candidate
        try:
            candidate = candidate.resolve(strict=True)
            candidate.relative_to(run_dir)
        except (OSError, ValueError):
            return None
        return candidate if candidate.is_file() else None

    if target_type == "sequence":
        fasta = stored_file(inputs.get("target_fasta"), "artifacts/raw/bindcraft2/target.fasta")
        if fasta is None:
            return None, "The original target FASTA is missing, so this job cannot be reproduced."
        target_pdb = None
    else:
        target_pdb = None
        for reference in (inputs.get("target_structure"), inputs.get("target_pdb")):
            target_pdb = stored_file(reference, "") if reference else None
            if target_pdb:
                break
        if target_pdb is None:
            for name in ("target.pdb", "target.cif", "target.mmcif"):
                target_pdb = stored_file(f"artifacts/raw/bindcraft2/{name}", "")
                if target_pdb:
                    break
        if target_pdb is None:
            return None, "The original target structure is missing, so this job cannot be reproduced."
        fasta = None

    def chains(value: object) -> list[str]:
        if isinstance(value, (list, tuple)):
            return [str(item) for item in value if str(item).strip()]
        return [item.strip() for item in str(value or "").split(",") if item.strip()]

    second_pdb = None
    if second:
        second_pdb = stored_file(inputs.get("second_target_pdb"), "artifacts/raw/bindcraft2/second_target.pdb")
        if second_pdb is None:
            return None, "The original second target is missing, so this job cannot be reproduced."
    second_weight = float(second.get("weight", 1.0) or 0.0)
    binder_length = params.get("binder_length")
    if not binder_length and isinstance(settings.get("binder_lengths"), (list, tuple)):
        binder_length = "-".join(str(value) for value in settings["binder_lengths"])
    quality = params.get("quality_filters")
    if not isinstance(quality, dict):
        quality = {key: settings[key] for key in (
            "min_plddt_screen", "min_plddt_refine", "min_plddt_anneal", "min_plddt_harden",
            "min_plddt_mutate", "min_plddt_final", "min_iptm_anneal", "min_iptm_harden",
            "min_iptm_mutate", "max_interdomain_contact_final", "min_domain_separation_ratio_final",
            "min_scaffold_sequence_retained_final", "max_scaffold_framework_rmsd_final",
        ) if key in settings}
        if isinstance(settings.get("filters"), dict):
            quality["filters"] = settings["filters"]
    aa_bias_setting = settings.get("aa_bias") if isinstance(settings.get("aa_bias"), dict) else {}
    excluded = params.get("excluded_amino_acids")
    aa_bias = params.get("amino_acid_bias")
    if not isinstance(excluded, list) or not isinstance(aa_bias, dict):
        excluded = [str(key) for key, value in aa_bias_setting.items() if float(value) == 0]
        aa_bias = {str(key): float(value) for key, value in aa_bias_setting.items() if float(value) != 0}
    crop_lengths = params.get("target_crop_lengths", settings.get("crop_fasta_sequence"))
    if isinstance(crop_lengths, tuple):
        crop_lengths = list(crop_lengths)

    return {
        "target_type": target_type,
        "target_fasta": str(fasta) if fasta else None,
        "target_pdb": str(target_pdb) if target_pdb else None,
        "target_chains": chains(params.get("target_chains", primary.get("chains"))),
        "binder_length": binder_length,
        "target_crop_lengths": crop_lengths,
        "target_crop_count": params.get("target_crop_count", settings.get("idr_crop_count")),
        "hotspots": params.get("hotspots", primary.get("hotspots", "")),
        "coldspots": params.get("coldspots", primary.get("coldspots", "")),
        "campaign_name": params.get("campaign_name") or settings.get("campaign_name") or job.get("campaign_name") or "bindcraft2",
        "modality": params.get("modality", settings.get("modality", "binder")),
        "design_properties": params.get("design_properties") or [key for key in bindcraft2.BC2_PROPERTY_PRESETS if settings.get(key)],
        "number_of_final_designs": int(params.get("number_of_final_designs", settings.get("number_of_final_designs", 1))),
        "max_trajectories": int(params.get("max_trajectories", settings.get("max_trajectories", 1))),
        "campaign_seed": int(params.get("campaign_seed", settings.get("campaign_seed", 0))),
        "workers_per_gpu": params.get("workers_per_gpu", settings.get("workers_per_gpu", "auto")),
        "redesign_interface": bool(params.get("redesign_interface", settings.get("redesign_interface", False))),
        "binder_copies": params.get("binder_copies", settings.get("copies")),
        "number_of_domains": params.get("number_of_domains", settings.get("n_domains")),
        "second_target_pdb": str(second_pdb) if second_pdb else None,
        "second_target_chains": chains(params.get("second_target_chains", second.get("chains"))),
        "second_target_hotspots": params.get("second_target_hotspots", second.get("hotspots", "")),
        "second_target_coldspots": params.get("second_target_coldspots", second.get("coldspots", "")),
        "second_target_goal": params.get("second_target_goal", "avoid" if second_weight < 0 else "bind"),
        "second_target_strength": float(params.get("second_target_strength", abs(second_weight))),
        "excluded_amino_acids": excluded,
        "amino_acid_bias": aa_bias,
        "quality_filters": quality,
        "scaffold_edits": params.get("scaffold_edits", settings.get("mutate_positions")),
        "engine_overrides": params.get("engine_overrides", {}),
        "trajectory_only": bool(params.get("trajectory_only") or settings.get("trajectory_only")),
        "cpu_cores": int(params.get("cpu_cores", (metadata.get("resource_request") or {}).get("cpu_cores", 4))),
        "gpu_device": str(params.get("gpu_device", (metadata.get("resource_allocation") or {}).get("gpu_device", "0")) or "0"),
    }, None


def _colored_amino_acid_sequence(value: object) -> str:
    sequence = re.sub(r"\s+", "", str(value or "")).upper()
    if sequence.casefold() in {"", "nan", "none"}:
        return ""
    residues = []
    for index, amino_acid in enumerate(sequence, start=1):
        group = next((name for name, members in _AA_TYPE_GROUPS.items() if amino_acid in members), "other")
        background, foreground = _AA_TYPE_STYLE[group]
        residues.append(
            f'<span title="{index}: {escape(amino_acid)} · {group}" '
            f'style="display:inline-block;min-width:1.05em;text-align:center;'
            f'margin:1px 0;padding:2px 0;border-radius:2px;background:{background};'
            f'color:{foreground};">{escape(amino_acid)}</span>'
        )
    return "".join(residues)


def _tool_output_column_width(column: object, sequence_column: object | None) -> int:
    """Choose starting widths; every column remains user-resizable."""
    if column == sequence_column:
        return 420
    normalized = re.sub(r"[^a-z0-9]+", "", str(column).casefold())
    if normalized in {"rank", "length", "binderlength"}:
        return 76
    if normalized in {"design", "candidate", "candidateid", "name", "modelid"}:
        return 240
    if "residue" in normalized or "position" in normalized or "setting" in normalized:
        return 180
    return max(104, min(164, len(str(column)) * 8 + 24))


def _render_tool_output_table(
    dataframe: pd.DataFrame,
    max_rows: int = 500,
    metric_color_rules: dict[str, tuple[str, float | None]] | None = None,
    scroll_height: int = 520,
) -> str:
    """Render an interactive tool-output grid with optional metric and sequence colors."""
    sequence_column = _sequence_column_name(dataframe.columns)
    shown = dataframe.head(max_rows)
    column_widths = [_tool_output_column_width(column, sequence_column) for column in dataframe.columns]
    table_width = max(sum(column_widths), 640)
    color_rules = metric_color_rules or {}
    numeric_values: dict[str, dict[int, float]] = {}
    color_ranges: dict[str, tuple[float, float]] = {}
    for column, (mode, _threshold) in color_rules.items():
        if column not in dataframe.columns:
            continue
        parsed: dict[int, float] = {}
        for row_index, value in dataframe[column].items():
            try:
                number = float(str(value).strip().removesuffix("%"))
            except (TypeError, ValueError):
                continue
            if pd.notna(number):
                parsed[row_index] = number
        numeric_values[column] = parsed
        if mode == "gradient_high" and parsed:
            color_ranges[column] = (min(parsed.values()), max(parsed.values()))

    def mix_color(start: tuple[int, int, int], end: tuple[int, int, int], fraction: float) -> str:
        fraction = min(1.0, max(0.0, fraction))
        channels = [round(left + (right - left) * fraction) for left, right in zip(start, end)]
        return "#" + "".join(f"{channel:02x}" for channel in channels)

    def metric_background(column: object, row_index: object) -> str | None:
        rule = color_rules.get(str(column))
        number = numeric_values.get(str(column), {}).get(row_index)
        if rule is None or number is None:
            return None
        mode, threshold = rule
        low = (249, 178, 178)  # soft red
        neutral = (216, 211, 194)  # warm neutral
        high = (139, 215, 182)  # mint green
        if mode == "threshold_high":
            return "#8bd7b6" if threshold is not None and number >= threshold else "#f9b2b2"
        if mode == "gradient_high":
            minimum, maximum = color_ranges.get(str(column), (number, number))
            fraction = (number - minimum) / (maximum - minimum) if maximum > minimum else 1.0
            return mix_color(low, neutral, fraction * 2) if fraction < 0.5 else mix_color(neutral, high, (fraction - 0.5) * 2)
        return None

    legend = "".join(
        f'<span class="mnds-aa-key"><i style="background:{background};"></i>{escape(name.title())}</span>'
        for name, (background, _) in _AA_TYPE_STYLE.items()
    )
    parts = ['<div class="mnds-results-shell">']
    if sequence_column is not None:
        parts.extend([
            '<div class="mnds-aa-legend"><span>AA type:</span>', legend, "</div>",
        ])
    parts.extend([
        '<div class="mnds-results-toolbar">',
        '<label>Find <input class="mnds-results-search" type="search" placeholder="Filter rows" aria-label="Filter table rows"></label>',
        '<label>Rows <select class="mnds-results-page-size" aria-label="Rows per page">',
        '<option>10</option><option selected>25</option><option>50</option><option>100</option><option>500</option>',
        '</select></label>',
        '<span class="mnds-results-count" aria-live="polite"></span>',
        '<button class="mnds-results-prev" type="button">Previous</button>',
        '<button class="mnds-results-next" type="button">Next</button>',
        '</div>',
        f'<div class="mnds-results-scroll" style="height:{max(120, int(scroll_height))}px"><table class="mnds-results-table" style="width:{table_width}px"><colgroup>',
    ])
    parts.extend(f'<col style="width:{width}px">' for width in column_widths)
    parts.append("</colgroup><thead><tr>")
    for index, column in enumerate(dataframe.columns):
        label = f"{column} · AA type" if column == sequence_column else str(column)
        parts.append(
            f'<th scope="col" data-column="{index}" title="Click to sort">'
            f'{escape(label)}<span class="mnds-column-resize" role="separator" tabindex="0" '
            f'aria-label="Resize {escape(label, quote=True)} column" title="Drag to resize"></span></th>'
        )
    parts.append("</tr></thead><tbody>")
    for row_index, row in shown.iterrows():
        parts.append("<tr>")
        for column in dataframe.columns:
            value = row[column]
            background = metric_background(column, row_index)
            cell_style = f' style="background:{background};"' if background else ""
            if pd.isna(value):
                cell = ""
            elif column == sequence_column:
                cell = (
                    f'<span class="mnds-sequence" title="{escape(str(value), quote=True)}">'
                    f'{_colored_amino_acid_sequence(value)}</span>'
                )
            else:
                cell = escape(str(value))
                if len(cell) > 500:
                    cell = f'<span title="{escape(str(value), quote=True)}">{cell[:497]}…</span>'
                elif cell:
                    cell = f'<span title="{escape(str(value), quote=True)}">{cell}</span>'
            parts.append(f"<td{cell_style}>{cell}</td>")
        parts.append("</tr>")
    parts.extend(("</tbody></table></div></div>", ""))
    if len(dataframe) > max_rows:
        parts.append(
            f'<div class="mnds-results-note">Showing the first {max_rows:,} of {len(dataframe):,} rows. '
            "Use the download button above to export the selected design set.</div>"
        )
    parts.append(
        """<script>
(() => {
  const root = document.currentScript.parentElement;
  const table = root.querySelector('.mnds-results-table');
  const body = table.tBodies[0];
  const columns = table.querySelectorAll('col');
  const search = root.querySelector('.mnds-results-search');
  const pageSizeSelect = root.querySelector('.mnds-results-page-size');
  const count = root.querySelector('.mnds-results-count');
  const previous = root.querySelector('.mnds-results-prev');
  const next = root.querySelector('.mnds-results-next');
  const rows = Array.from(body.rows);
  let page = 0;
  let sortColumn = -1;
  let sortDirection = 1;

  function filteredRows() {
    const term = search.value.trim().toLocaleLowerCase();
    return rows.filter(row => !term || row.textContent.toLocaleLowerCase().includes(term));
  }

  function renderPage() {
    const filtered = filteredRows();
    const pageSize = Number(pageSizeSelect.value) || 25;
    const pageCount = Math.max(1, Math.ceil(filtered.length / pageSize));
    page = Math.min(page, pageCount - 1);
    const start = page * pageSize;
    const visible = new Set(filtered.slice(start, start + pageSize));
    for (const row of rows) row.hidden = !visible.has(row);
    count.textContent = filtered.length
      ? `Rows ${start + 1}–${Math.min(start + pageSize, filtered.length)} of ${filtered.length}`
      : 'No matching rows';
    previous.disabled = page === 0;
    next.disabled = page >= pageCount - 1;
  }

  search.addEventListener('input', () => { page = 0; renderPage(); });
  pageSizeSelect.addEventListener('change', () => { page = 0; renderPage(); });
  previous.addEventListener('click', () => { page -= 1; renderPage(); });
  next.addEventListener('click', () => { page += 1; renderPage(); });

  for (const header of table.tHead.rows[0].cells) {
    const index = Number(header.dataset.column);
    header.addEventListener('click', event => {
      if (event.target.closest('.mnds-column-resize')) return;
      sortDirection = sortColumn === index ? -sortDirection : 1;
      sortColumn = index;
      const value = row => row.cells[index]?.textContent.trim() || '';
      rows.sort((left, right) => {
        const a = value(left);
        const b = value(right);
        const numberPattern = /^[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[-+]?\\d+)?%?$/i;
        if (numberPattern.test(a) && numberPattern.test(b)) {
          return (Number(a.replace(/%$/, '')) - Number(b.replace(/%$/, ''))) * sortDirection;
        }
        return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' }) * sortDirection;
      });
      rows.forEach(row => body.appendChild(row));
      table.querySelectorAll('th[aria-sort]').forEach(cell => cell.removeAttribute('aria-sort'));
      header.setAttribute('aria-sort', sortDirection === 1 ? 'ascending' : 'descending');
      page = 0;
      renderPage();
    });

    const handle = header.querySelector('.mnds-column-resize');
    handle.addEventListener('pointerdown', event => {
      event.preventDefault();
      event.stopPropagation();
      const startX = event.clientX;
      const startWidth = header.getBoundingClientRect().width;
      const startTableWidth = table.getBoundingClientRect().width;
      const onMove = moveEvent => {
        const width = Math.max(56, startWidth + moveEvent.clientX - startX);
        const delta = width - startWidth;
        columns[index].style.width = `${width}px`;
        table.style.width = `${Math.max(640, startTableWidth + delta)}px`;
      };
      const onUp = () => {
        document.removeEventListener('pointermove', onMove);
        document.removeEventListener('pointerup', onUp);
        document.removeEventListener('pointercancel', onUp);
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
      };
      document.body.style.cursor = 'col-resize';
      document.body.style.userSelect = 'none';
      document.addEventListener('pointermove', onMove);
      document.addEventListener('pointerup', onUp, { once: true });
      document.addEventListener('pointercancel', onUp, { once: true });
    });
    handle.addEventListener('keydown', event => {
      if (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight') return;
      event.preventDefault();
      const currentWidth = columns[index].getBoundingClientRect().width;
      const width = Math.max(56, currentWidth + (event.key === 'ArrowRight' ? 24 : -24));
      const appliedDelta = width - currentWidth;
      columns[index].style.width = `${width}px`;
      table.style.width = `${Math.max(640, table.getBoundingClientRect().width + appliedDelta)}px`;
    });
  }

  renderPage();
})();
</script>"""
    )
    return "".join(parts)


def _render_tool_output_table_css() -> str:
    return """
    <style>
    .mnds-results-shell { width:100%; color:#17212b; }
    .mnds-aa-legend { display:flex; align-items:center; gap:10px; flex-wrap:wrap; margin:0 0 8px; font:12px sans-serif; color:#52606d; }
    .mnds-aa-key { display:inline-flex; align-items:center; gap:4px; }
    .mnds-aa-key i { width:10px; height:10px; display:inline-block; border-radius:2px; }
    .mnds-results-toolbar { display:flex; align-items:center; gap:10px; flex-wrap:wrap; margin:0 0 8px; font:12px sans-serif; color:#52606d; }
    .mnds-results-toolbar label { display:inline-flex; align-items:center; gap:5px; }
    .mnds-results-toolbar input, .mnds-results-toolbar select, .mnds-results-toolbar button { border:1px solid #d5dbe1; border-radius:5px; padding:5px 8px; background:#fff; color:#26352d; font:inherit; }
    .mnds-results-toolbar input { min-width:160px; }
    .mnds-results-toolbar button { cursor:pointer; }
    .mnds-results-toolbar button:disabled { cursor:default; opacity:.45; }
    .mnds-results-count { margin-left:auto; }
    .mnds-results-scroll { width:100%; height:520px; overflow:auto; border:1px solid #dfe3e8; border-radius:7px; }
    .mnds-results-table { border-collapse:separate; border-spacing:0; min-width:100%; table-layout:fixed; font:13px sans-serif; }
    .mnds-results-table th { position:sticky; top:0; z-index:1; background:#f0f2f5; color:#4d5966; text-align:left; font-size:11px; text-transform:uppercase; letter-spacing:.02em; cursor:pointer; user-select:none; }
    .mnds-results-table th, .mnds-results-table td { position:relative; padding:8px 10px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
    .mnds-results-table th { padding-right:14px; }
    .mnds-column-resize { position:absolute; z-index:2; top:0; right:0; width:10px; height:100%; cursor:col-resize; touch-action:none; background:transparent; }
    .mnds-column-resize:hover, .mnds-column-resize:focus { background:#8daf9a; outline:none; }
    .mnds-results-table td { background:#fff; }
    .mnds-results-table tr:nth-child(even) td { background:#fafbfc; }
    .mnds-results-table .mnds-sequence { display:block; white-space:nowrap; overflow:hidden; text-overflow:clip; font:13px/1.3 ui-monospace,SFMono-Regular,Menlo,monospace; }
    .mnds-results-note { margin-top:6px; color:#687684; font-size:12px; }
    </style>
    """


def _structure_design_status(item: dict, path: Path) -> str:
    relative_path = str(item.get("relative_path", path)).casefold().replace("\\", "/")
    if "2_refolded" in relative_path or "not_accepted" in path.stem.casefold():
        return "Not accepted"
    if "3_ranked" in relative_path:
        return "Accepted"
    return "Other"


def _structure_design_name(path: Path) -> str:
    name = path.stem
    for suffix in ("_best_candidate_not_accepted", "_not_accepted"):
        if name.casefold().endswith(suffix):
            return name[: -len(suffix)]
    return name


@st.cache_data(show_spinner=False)
def _mmcif_binder_sequence(path_text: str, modified_ns: int) -> str:
    """Read the binder polymer sequence recorded in an mmCIF result file."""
    path = Path(path_text)
    try:
        import gemmi

        block = gemmi.cif.read_file(str(path)).sole_block()
        entity_ids = list(block.find_values("_entity.id"))
        descriptions = list(block.find_values("_entity.pdbx_description"))
        description_by_id = {
            entity_id: descriptions[index] if index < len(descriptions) else ""
            for index, entity_id in enumerate(entity_ids)
        }
        poly_ids = list(block.find_values("_entity_poly.entity_id"))
        sequences = list(block.find_values("_entity_poly.pdbx_seq_one_letter_code_can"))
        candidates = []
        for index, entity_id in enumerate(poly_ids):
            raw_sequence = sequences[index] if index < len(sequences) else ""
            sequence = re.sub(r"[^A-Z]", "", str(raw_sequence).upper())
            if sequence:
                candidates.append(("binder" in str(description_by_id.get(entity_id, "")).casefold(), sequence))
        if candidates:
            binder_sequences = [sequence for is_binder, sequence in candidates if is_binder]
            return max(binder_sequences or [sequence for _, sequence in candidates], key=len)
    except (ImportError, OSError, RuntimeError, ValueError, IndexError):
        pass

    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict

        values = MMCIF2Dict(str(path))
        entity_ids = values.get("_entity.id", [])
        descriptions = values.get("_entity.pdbx_description", [])
        if isinstance(entity_ids, str):
            entity_ids = [entity_ids]
        if isinstance(descriptions, str):
            descriptions = [descriptions]
        description_by_id = {
            str(entity_id): descriptions[index] if index < len(descriptions) else ""
            for index, entity_id in enumerate(entity_ids)
        }
        poly_ids = values.get("_entity_poly.entity_id", [])
        sequences = values.get("_entity_poly.pdbx_seq_one_letter_code_can", [])
        if isinstance(poly_ids, str):
            poly_ids = [poly_ids]
        if isinstance(sequences, str):
            sequences = [sequences]
        candidates = []
        for index, entity_id in enumerate(poly_ids):
            raw_sequence = sequences[index] if index < len(sequences) else ""
            sequence = re.sub(r"[^A-Z]", "", str(raw_sequence).upper())
            if sequence:
                candidates.append(("binder" in str(description_by_id.get(str(entity_id), "")).casefold(), sequence))
        if candidates:
            binder_sequences = [sequence for is_binder, sequence in candidates if is_binder]
            return max(binder_sequences or [sequence for _, sequence in candidates], key=len)
    except (ImportError, OSError, KeyError, ValueError, IndexError):
        pass
    return ""


def _structure_summary_rows(structures: list[tuple[dict, Path]], columns=()) -> pd.DataFrame:
    """Make minimal table rows for structures without a ranked CSV record."""
    columns = list(columns)
    normalized = {re.sub(r"[^a-z0-9]+", "", str(column).casefold()): column for column in columns}
    design_column = next(
        (normalized[key] for key in ("design", "candidateid", "candidate", "name", "modelid") if key in normalized),
        None,
    )
    sequence_column = _sequence_column_name(columns)
    length_column = normalized.get("length", normalized.get("binderlength"))
    if design_column is None:
        design_column = "Design"
        columns.append(design_column)
    if sequence_column is None:
        sequence_column = "Binder_Sequence"
        columns.append(sequence_column)
    if length_column is None:
        length_column = "Length"
        columns.append(length_column)
    if "Design Status" not in columns:
        columns.append("Design Status")

    rows = []
    for item, path in structures:
        row = {column: None for column in columns}
        row[design_column] = _structure_design_name(path)
        try:
            sequence = _mmcif_binder_sequence(str(path), path.stat().st_mtime_ns)
        except OSError:
            sequence = ""
        row[sequence_column] = sequence
        row[length_column] = len(sequence) if sequence else None
        row["Design Status"] = _structure_design_status(item, path)
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def _bindcraft2_csv_kind(path: Path, dataframe: pd.DataFrame) -> str | None:
    normalized = {
        re.sub(r"[^a-z0-9]+", "", str(column).casefold()) for column in dataframe.columns
    }
    if {"trajectory", "design", "length", "bindcraftversion", "terminated"}.issubset(normalized):
        return "trajectories"
    if {"campaign", "scope", "metric", "samples", "mean", "std", "min", "max"}.issubset(normalized):
        return "summary"
    if {"phase", "round"}.issubset(normalized) and any(
        str(column).casefold().startswith("app_target.") for column in dataframe.columns
    ):
        return "losses"
    return None


_IPSAE_RESULT_COLUMN_ORDER = (
    # Put the scores calculated by this analysis first so they are visible
    # immediately, before the inherited BindCraft2 metrics and sequence.
    "ipsae",
    "ipsae_min",
    "ipsae_max",
    "ipsae_avg",
    "ipsae_d0chn",
    "ipsae_d0dom",
    "lis",
    "ipae",
    "iptm_af",
    "iptm_d0chn",
    "pdockq_min",
    "pdockq_max",
    "pdockq2_min",
    "pdockq2_max",
    # Keep the original design ranking visible without changing it.
    "rank",
    "bindcraft2_rank",
    "candidate_id",
    "binder_sequence",
    "length",
    "binder_chain",
    "target_chain",
    # Original prediction confidence scores and complex quality measures.
    "i_pTM",
    "i_pAE",
    "i_pDAE",
    "pTM",
    "pLDDT",
    "Target_pLDDT",
    "Unbound_Binder_pLDDT",
    "SS_pLDDT",
    "Binder_RMSD",
    "Target_RMSD",
    "Interface_BuriedArea",
    "Interface_BuriedArea_Fraction",
    "Interface_Residues",
    "Epitope_Residues_Contacted",
    "Hotspot_Contact_Fraction",
    "Off_Epitope_Contact_Fraction",
    "All_Atom_Clashes",
    "Backbone_Clashes",
    # Compact developability / composition indicators.
    "Binder_Mass_kDa",
    "Binder_Net_Charge",
    "Surface_Hydrophobicity",
    "Binder_Cysteines",
    "Binder_Disulfides",
)


def _ipsae_display_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Return the score- and complex-quality-focused ipSAE view."""
    available = {str(column).casefold(): column for column in dataframe.columns}
    columns = [
        available[name.casefold()]
        for name in _IPSAE_RESULT_COLUMN_ORDER
        if name.casefold() in available
    ]
    return dataframe.loc[:, columns].copy() if columns else dataframe.copy()


def _bindcraft2_report_read_json(path: Path) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _bindcraft2_report_read_csv(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, low_memory=False)
    except (OSError, ValueError, pd.errors.ParserError, UnicodeError):
        return pd.DataFrame()


def _bindcraft2_report_counts(dataframe: pd.DataFrame) -> tuple[int, int]:
    """Return passed and rejected refold counts from BindCraft2's native table."""
    if dataframe.empty:
        return 0, 0
    columns = {
        re.sub(r"[^a-z0-9]+", "", str(column).casefold()): column
        for column in dataframe.columns
    }
    outcome_column = columns.get("outcome")
    failed_column = columns.get("failedfilters")
    if outcome_column is not None:
        outcomes = dataframe[outcome_column].fillna("").astype(str).str.strip().str.casefold()
        passed = int(outcomes.isin({"passed", "pass", "accepted"}).sum())
        rejected = int((~outcomes.isin({"passed", "pass", "accepted"}) & outcomes.ne("")).sum())
        if passed + rejected < len(dataframe) and failed_column is not None:
            failed = dataframe[failed_column].fillna("").astype(str).str.strip()
            unresolved = outcomes.eq("")
            passed += int((unresolved & failed.eq("")).sum())
            rejected += int((unresolved & failed.ne("")).sum())
        return passed, rejected
    if failed_column is not None:
        failed = dataframe[failed_column].fillna("").astype(str).str.strip()
        return int(failed.eq("").sum()), int(failed.ne("").sum())
    return 0, 0


def _bindcraft2_report_filter_rejections(dataframe: pd.DataFrame) -> dict[str, int]:
    if dataframe.empty:
        return {}
    columns = {
        re.sub(r"[^a-z0-9]+", "", str(column).casefold()): column
        for column in dataframe.columns
    }
    failed_column = columns.get("failedfilters")
    if failed_column is None:
        return {}
    counts: dict[str, int] = {}
    for value in dataframe[failed_column].fillna("").astype(str):
        reasons = [part.strip() for part in re.split(r"[;,|]", value) if part.strip()]
        for reason in reasons:
            counts[reason] = counts.get(reason, 0) + 1
    return counts


def _bindcraft2_report_closest_misses(dataframe: pd.DataFrame, filters: dict) -> list[dict[str, object]]:
    if dataframe.empty or not isinstance(filters, dict):
        return []
    normalized_columns = {
        re.sub(r"[^a-z0-9]+", "", str(column).casefold()): column
        for column in dataframe.columns
    }
    normalized_filters = {
        re.sub(r"[^a-z0-9]+", "", str(name).casefold()): (str(name), config)
        for name, config in filters.items()
        if isinstance(config, dict)
    }
    failed_column = normalized_columns.get("failedfilters")
    design_column = normalized_columns.get("design")
    if failed_column is None:
        return []

    misses: list[dict[str, object]] = []
    for _, row in dataframe.iterrows():
        reasons = [part.strip() for part in re.split(r"[;,|]", str(row.get(failed_column) or "")) if part.strip()]
        for reason in reasons:
            reason_key = re.sub(r"[^a-z0-9]+", "", reason.casefold())
            matched = normalized_filters.get(reason_key)
            if matched is None:
                matched = next(
                    (pair for key, pair in normalized_filters.items() if key and (key in reason_key or reason_key in key)),
                    None,
                )
            if matched is None:
                continue
            filter_name, config = matched
            metric_key = re.sub(r"[^a-z0-9]+", "", filter_name.casefold())
            metric_column = normalized_columns.get(metric_key)
            if metric_column is None:
                continue
            try:
                value = float(row.get(metric_column))
                threshold = float(config.get("threshold"))
            except (TypeError, ValueError):
                continue
            higher = bool(config.get("higher"))
            short_by = threshold - value if higher else value - threshold
            if short_by < 0:
                continue
            misses.append(
                {
                    "Design": str(row.get(design_column) or "—") if design_column else "—",
                    "Filter": filter_name.replace("_", " "),
                    "Got": value,
                    "Needed": f"{'≥' if higher else '≤'} {threshold:g}",
                    "Short by": short_by,
                }
            )
    misses.sort(key=lambda item: (float(item["Short by"]), str(item["Design"]).casefold()))
    return misses[:12]


def _bindcraft2_campaign_report(run_dir: Path, input_payload: dict, result: dict) -> dict | None:
    """Read BindCraft2's saved native outputs for the in-page campaign report."""
    raw_dir = run_dir / "artifacts" / "raw" / "bindcraft2"
    output_dir = raw_dir / "output"
    trajectory_csv = _bindcraft2_report_read_csv(output_dir / "1_Trajectories" / "!_Trajectories.csv")
    refolded_csv = _bindcraft2_report_read_csv(output_dir / "2_Refolded" / "!_Refolded.csv")
    ranked_csv = _bindcraft2_report_read_csv(output_dir / "3_Ranked" / "!_Ranked.csv")
    if not any(path.is_file() for path in (
        output_dir / "1_Trajectories" / "!_Trajectories.csv",
        output_dir / "2_Refolded" / "!_Refolded.csv",
        output_dir / "3_Ranked" / "!_Ranked.csv",
        output_dir / "campaign_metadata.json",
    )):
        return None

    app_settings = _bindcraft2_report_read_json(raw_dir / "settings.json")
    campaign = _bindcraft2_report_read_json(output_dir / "campaign_metadata.json")
    state = _bindcraft2_report_read_json(output_dir / ".campaign_state.json")
    settings = campaign.get("settings") if isinstance(campaign.get("settings"), dict) else app_settings
    settings = settings if isinstance(settings, dict) else {}
    inputs = input_payload.get("inputs") if isinstance(input_payload.get("inputs"), dict) else {}
    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    targets = settings.get("targets") if isinstance(settings.get("targets"), list) else []
    target = targets[0] if targets and isinstance(targets[0], dict) else {}
    resolved = campaign.get("resolved") if isinstance(campaign.get("resolved"), dict) else {}
    state_rejections = state.get("rejections") if isinstance(state.get("rejections"), dict) else {}

    passed_count, rejected_count = _bindcraft2_report_counts(refolded_csv)
    candidate_metrics = result.get("metrics") if isinstance(result.get("metrics"), dict) else {}
    candidate_count = int(candidate_metrics.get("accepted_design_count") or candidate_metrics.get("candidate_count") or 0)
    trajectory_count = len(trajectory_csv) or int(state.get("trajectories") or 0)
    refolded_count = len(refolded_csv) or int(state_rejections.get("candidates_scored") or 0)
    ranked_count = len(ranked_csv) or candidate_count or int(state.get("accepted") or 0)
    if refolded_csv.empty and not passed_count:
        passed_count = ranked_count
    if refolded_csv.empty:
        rejected_count = max(0, int(state_rejections.get("candidates_rejected") or 0))

    filters = settings.get("filters") if isinstance(settings.get("filters"), dict) else {}
    filter_rejections = state_rejections.get("failed_filters")
    if not isinstance(filter_rejections, dict) or not filter_rejections:
        filter_rejections = _bindcraft2_report_filter_rejections(refolded_csv)
    normalized_filter_rejections = {}
    for name, count in filter_rejections.items():
        try:
            count_value = int(float(count))
        except (TypeError, ValueError):
            continue
        if str(name).strip() and count_value > 0:
            normalized_filter_rejections[str(name)] = count_value
    filter_rejections = normalized_filter_rejections
    closest_misses = _bindcraft2_report_closest_misses(refolded_csv, filters)

    modality = settings.get("modality") or params.get("modality") or "binder"
    if isinstance(modality, (list, tuple)):
        modality = ", ".join(str(item) for item in modality if str(item).strip())
    modality_label = str(modality).replace("_", " ").strip().title() or "Binder"
    objective = params.get("conformational_objective") or params.get("objective")
    if objective:
        modality_label += f" · {str(objective).replace('_', ' ').title()} objective"

    raw_lengths = settings.get("binder_lengths") or params.get("binder_length")
    if isinstance(raw_lengths, (list, tuple)) and len(raw_lengths) >= 2:
        length_label = f"{raw_lengths[0]}–{raw_lengths[1]} aa"
    elif isinstance(raw_lengths, (list, tuple)) and raw_lengths:
        length_label = f"{raw_lengths[0]} aa"
    else:
        length_label = str(raw_lengths or "—").replace("-", "–")
        if length_label != "—" and "aa" not in length_label.casefold():
            length_label += " aa"

    target_type = str(inputs.get("target_type") or params.get("target_type") or "structure").casefold()
    target_type_label = "Sequence (disordered)" if target_type in {"sequence", "idr", "disordered"} else "Structure"
    target_chains = inputs.get("target_chains") or params.get("target_chains") or target.get("chains") or []
    if isinstance(target_chains, str):
        target_chains = [part.strip() for part in target_chains.split(",") if part.strip()]
    chain_label = ", ".join(str(chain) for chain in target_chains) or "—"
    hotspots = params.get("hotspots") or target.get("hotspots") or ""
    target_name = (
        inputs.get("target_name")
        or params.get("target_name")
        or target.get("name")
        or target.get("target_name")
    )
    target_label = f"{target_type_label} · chain {chain_label}"
    if target_name:
        target_label = f"{target_name} · {target_label}"
    if hotspots:
        target_label += f" · hotspots {hotspots}"

    stages = ("Screen", "Refine", "Anneal", "Harden", "Mutate", "Final")
    stage_values = {
        "Interface i_pTM": {
            "Anneal": settings.get("min_iptm_anneal"),
            "Harden": settings.get("min_iptm_harden"),
            "Mutate": settings.get("min_iptm_mutate"),
            "Final": settings.get("min_iptm_final"),
        },
        "Binder pLDDT": {
            "Screen": settings.get("min_plddt_screen"),
            "Refine": settings.get("min_plddt_refine"),
            "Anneal": settings.get("min_plddt_anneal"),
            "Harden": settings.get("min_plddt_harden"),
            "Mutate": settings.get("min_plddt_mutate"),
            "Final": settings.get("min_plddt_final"),
        },
        "Unbound binder pLDDT": {"Final": settings.get("min_monomer_plddt_final")},
    }
    stage_rows = []
    for metric, values in stage_values.items():
        row = {"Metric": metric}
        for stage in stages:
            value = values.get(stage)
            try:
                row[stage] = f"≥ {float(value):g}" if value is not None else "—"
            except (TypeError, ValueError):
                row[stage] = "—"
        if any(value != "—" for value in row.values() if value != metric):
            stage_rows.append(row)

    other_filter_rows = []
    for name, config in filters.items():
        if not isinstance(config, dict):
            continue
        threshold = config.get("threshold")
        if threshold in (None, "Infinity", "-Infinity", "inf", "-inf"):
            continue
        try:
            numeric_threshold = float(threshold)
        except (TypeError, ValueError):
            continue
        if numeric_threshold == 0 and not config.get("mandatory") and name != "Backbone_Clashes":
            continue
        direction = "≥" if config.get("higher") else "≤"
        other_filter_rows.append({"Filter": str(name).replace("_", " "), "Threshold": f"{direction} {numeric_threshold:g}"})
    other_filter_rows.sort(key=lambda row: row["Filter"].casefold())

    settings_rows = [
        ("Design format", modality_label),
        ("Target", target_label),
        ("Binder length range", length_label),
        ("Designs requested", settings.get("number_of_final_designs", params.get("number_of_final_designs", "—"))),
        ("Sequences tried per design", settings.get("sequence_candidates", "—")),
        ("Passing sequences sought per design", settings.get("enough_passing_sequences", "—")),
        ("Random seed", settings.get("campaign_seed", params.get("campaign_seed", "—"))),
        ("Maximum trajectories", settings.get("max_trajectories", params.get("max_trajectories", "—"))),
    ]
    termination = state_rejections.get("terminated") if isinstance(state_rejections.get("terminated"), dict) else {}
    if not termination and not trajectory_csv.empty:
        terminated_column = next((column for column in trajectory_csv.columns if str(column).casefold() == "terminated"), None)
        if terminated_column:
            termination = trajectory_csv[terminated_column].fillna("").astype(str).str.strip().value_counts().to_dict()
    file_guide = [
        ("artifacts/raw/bindcraft2/output/1_Trajectories/!_Trajectories.csv", "Design attempts and trajectory-level scores."),
        ("artifacts/raw/bindcraft2/output/2_Refolded/!_Refolded.csv", "Redesigned candidates, refold scores, and filter outcomes."),
        ("artifacts/raw/bindcraft2/output/3_Ranked/!_Ranked.csv", "Accepted candidates retained and ranked by BindCraft2."),
        ("artifacts/raw/bindcraft2/output/summary.csv", "Campaign-level aggregate metrics."),
        ("artifacts/raw/bindcraft2/output/campaign_metadata.json", "Resolved settings, stage steps, software version, and model metadata."),
        ("artifacts/raw/bindcraft2/settings.json", "Input settings requested for this BindCraft2 run."),
    ]
    file_guide = [
        {"File": relative, "Contents": description}
        for relative, description in file_guide
        if (run_dir / relative).is_file()
    ]

    return {
        "campaign_name": campaign.get("campaign_name") or settings.get("campaign_name") or params.get("campaign_name") or "BindCraft2 campaign",
        "version": campaign.get("version") or "",
        "modality": modality_label,
        "trajectories": trajectory_count,
        "refolded": refolded_count,
        "passed": passed_count,
        "rejected": rejected_count,
        "ranked": ranked_count,
        "settings_rows": settings_rows,
        "score_definitions": [
            {"Score": "i_pTM", "Meaning": "Predicted confidence in the binder–target interface; higher is better."},
            {"Score": "i_pAE", "Meaning": "Normalized predicted aligned error across the interface; lower is better."},
            {"Score": "i_pDAE", "Meaning": "BindCraft2 distance-masked interface confidence used for ranking; higher is better."},
            {"Score": "pLDDT", "Meaning": "Mean local structure confidence for the binder; higher is better."},
            {"Score": "pTM", "Meaning": "Predicted confidence in the overall complex structure; higher is better."},
            {"Score": "Unbound_Binder_pLDDT", "Meaning": "Predicted binder confidence when folded without the target; higher is better."},
            {"Score": "Interface_Residues", "Meaning": "Number of binder residues participating in the predicted interface."},
            {"Score": "Hotspot_Contact_Fraction", "Meaning": "Fraction of requested target hotspots contacted by the binder."},
            {"Score": "All_Atom_Clashes / Backbone_Clashes", "Meaning": "Detected atomic overlaps; lower is better, with zero ideal."},
        ],
        "stage_rows": stage_rows,
        "filters": other_filter_rows,
        "filter_rejections": filter_rejections,
        "closest_misses": closest_misses,
        "termination": termination,
        "has_refold_diagnostics": not refolded_csv.empty,
        "file_guide": file_guide,
    }


def _render_bindcraft2_campaign_report(run_dir: Path, input_payload: dict, result: dict) -> None:
    report = _bindcraft2_campaign_report(run_dir, input_payload, result)
    if report is None:
        return

    st.divider()
    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    if bool(params.get("trajectory_only")):
        st.subheader("BindCraft2 generator summary")
        st.caption("Only trajectories that passed BindCraft2's gradient-design stage gates are forwarded. These are eligible for ProteinMPNN sequence redesign; they have not passed the later sequence, refolding, or acceptance filters.")
        generated_count = int((result.get("metrics") or {}).get("generated_trajectory_count") or 0)
        metric_cols = st.columns(2)
        metric_cols[0].metric("Gradient trajectories run", report["trajectories"])
        metric_cols[1].metric("Eligible for ProteinMPNN", generated_count)
        st.caption("Continue these structures through ProteinMPNN, then choose a structure prediction workflow to evaluate redesigned sequences.")
        return

    st.subheader("BindCraft2 campaign summary")
    description = report["modality"]
    if report["version"]:
        description += f" · {report['version']}"
    st.caption(description)

    count_rows = [
        ("Trajectories run", report["trajectories"], "Design attempts started from a fresh initialization."),
        ("Candidates refolded", report["refolded"], "Redesigned sequences predicted again as target–binder complexes."),
        ("Passed every filter", report["passed"], "Refolded candidates that cleared the saved acceptance filters."),
        ("Designs ranked", report["ranked"], "Accepted designs retained and ranked by BindCraft2."),
    ]
    bar_max = max((count for _, count, _ in count_rows), default=0)
    summary_html = [
        '<div style="display:flex;flex-direction:column;gap:1rem;margin:.35rem 0 1.2rem">'
    ]
    for label, count, help_text in count_rows:
        width = max(0, min(100, round(100 * count / bar_max))) if bar_max else 0
        summary_html.append(
            '<div style="display:grid;grid-template-columns:minmax(180px,1.3fr) 3rem minmax(120px,1.5fr);'
            'column-gap:1rem;row-gap:.1rem;align-items:center">'
            f'<span style="font-size:1rem;color:#414b55">{escape(label)}</span>'
            f'<strong style="font-size:1.2rem;text-align:right;color:#111">{int(count)}</strong>'
            '<span style="height:.75rem;background:#edf1f5;border-radius:2px;overflow:hidden">'
            f'<span style="display:block;width:{width}%;height:100%;background:#507da8"></span></span>'
            f'<span style="grid-column:1 / -1;font-size:.82rem;color:#81878c">{escape(help_text)}</span>'
            '</div>'
        )
    summary_html.append('</div>')
    st.markdown("".join(summary_html), unsafe_allow_html=True)
    if report["refolded"] > report["trajectories"]:
        st.caption("Refolded candidates can outnumber trajectories because one trajectory can produce several redesigned sequences.")

    st.markdown("#### What this campaign asked for")
    settings_frame = pd.DataFrame(report["settings_rows"], columns=["Setting", "Value"])
    st.dataframe(settings_frame, hide_index=True, use_container_width=True, height=min(360, 38 + 36 * len(settings_frame)))

    if report["stage_rows"]:
        st.markdown("#### Acceptance thresholds by stage")
        st.caption("Thresholds are read from the effective BindCraft2 settings saved with this run.")
        st.dataframe(pd.DataFrame(report["stage_rows"]), hide_index=True, use_container_width=True)
    if report["filters"]:
        with st.expander("Other active acceptance filters", expanded=False):
            st.dataframe(pd.DataFrame(report["filters"]), hide_index=True, use_container_width=True)

    st.markdown("#### What the scores mean")
    st.caption("These score values are reported by BindCraft2 in its native CSV output; the descriptions below are explanatory labels from Studio.")
    st.dataframe(pd.DataFrame(report["score_definitions"]), hide_index=True, use_container_width=True)

    st.markdown("#### Why designs were rejected")
    if report["rejected"] <= 0:
        if report["has_refold_diagnostics"]:
            st.caption(f"No refolded candidates were rejected in this run: {report['passed']} cleared every filter.")
        else:
            st.caption("The run did not produce a refolded-candidate table, so filter rejection reasons are unavailable.")
    elif report["filter_rejections"]:
        rejection_rows = sorted(report["filter_rejections"].items(), key=lambda item: (-item[1], item[0].casefold()))
        rejection_frame = pd.DataFrame(rejection_rows, columns=["Filter", "Candidates rejected"])
        st.bar_chart(rejection_frame.set_index("Filter"), horizontal=True)
        st.caption("Counts come from BindCraft2's native failed_filters field; candidates may fail more than one filter.")
    else:
        st.caption(f"{report['rejected']} candidates were rejected, but the saved native outputs do not identify their filter reasons.")
    if report["closest_misses"]:
        st.markdown("**Closest to passing**")
        st.dataframe(pd.DataFrame(report["closest_misses"]), hide_index=True, use_container_width=True)
    termination = report["termination"]
    if termination:
        termination_text = ", ".join(f"{name or 'not recorded'}: {count}" for name, count in sorted(termination.items()))
        st.caption(f"Trajectory outcomes: {termination_text}.")
    st.caption("The original BindCraft2 CSVs, settings, and structures are available in the Files tab and Download job ZIP.")
    if report["file_guide"]:
        st.markdown("#### Where the files are")
        st.dataframe(pd.DataFrame(report["file_guide"]), hide_index=True, use_container_width=True)


def _render_boltzgen_campaign_report(run_dir: Path, result: dict) -> None:
    """Show native BoltzGen scores and a saved refold PAE matrix."""
    candidates = result.get("candidates") if isinstance(result.get("candidates"), list) else []
    if not candidates:
        candidates_path = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
        candidates = read_candidates(candidates_path) if candidates_path.is_file() else []
    if not candidates:
        return

    st.divider()
    st.subheader("BoltzGen campaign summary")
    pae_candidates = []
    metric_rows = []
    for index, candidate in enumerate(candidates):
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        raw_pae_path = raw.get("pae_path") or raw.get("pae_source_path")
        pae_path = Path(str(raw_pae_path)) if raw_pae_path else None
        if pae_path is not None and not pae_path.is_absolute():
            pae_path = run_dir / pae_path
        if pae_path is not None and pae_path.is_file():
            pae_candidates.append((index, candidate, metrics, pae_path))
        metric_rows.append(
            {
                "Design": candidate.get("candidate_id") or f"Design {index + 1}",
                "Stage": candidate.get("stage") or "—",
                "Full PAE saved": bool(pae_path and pae_path.is_file()),
                "Minimum design ipSAE": metrics.get("design_ipsae_min"),
                "Design→target ipSAE": metrics.get("design_to_target_ipsae"),
                "Target→design ipSAE": metrics.get("target_to_design_ipsae"),
                "Design ipTM": metrics.get("design_iptm"),
                "Target ipTM": metrics.get("target_iptm"),
                "pTM": metrics.get("ptm"),
            }
        )
    full_pae_count = sum(bool(row["Full PAE saved"]) for row in metric_rows)
    summary_cols = st.columns(3)
    summary_cols[0].metric("Design records", len(metric_rows))
    summary_cols[1].metric("Full PAE matrices", f"{full_pae_count} / {len(metric_rows)}")
    summary_cols[2].metric("Native ranked designs", sum(row["Stage"] == STAGE_COMPLEX_REFOLDING for row in metric_rows))
    st.dataframe(pd.DataFrame(metric_rows), hide_index=True, use_container_width=True)
    if not pae_candidates:
        st.caption("A full PAE matrix is written during the native refolding stage. Generation-only runs have no refold PAE matrix.")
        return

    selection = st.selectbox(
        "Refold PAE matrix",
        range(len(pae_candidates)),
        format_func=lambda item: str(pae_candidates[item][1].get("candidate_id") or f"Design {item + 1}"),
        key=f"boltzgen-pae-design-{run_dir.name}",
    )
    _candidate_index, candidate, _metrics, pae_path = pae_candidates[selection]
    try:
        import matplotlib.pyplot as plt
        import numpy as np

        with np.load(pae_path, allow_pickle=False) as archive:
            pae_key = next((key for key in archive.files if "pae" in key.casefold()), None)
            if pae_key is None:
                raise ValueError("This archive does not contain a PAE array.")
            matrix = np.asarray(archive[pae_key], dtype=float)
        sample_label = ""
        if matrix.ndim > 2:
            sample_count = int(np.prod(matrix.shape[:-2]))
            matrix = np.nanmedian(matrix.reshape(sample_count, *matrix.shape[-2:]), axis=0)
            sample_label = f" · median of {sample_count} refolds"
        if matrix.ndim != 2:
            raise ValueError(f"Expected a 2D PAE matrix, found shape {matrix.shape}.")
        fig, axis = plt.subplots(figsize=(7.2, 5.8), constrained_layout=True)
        image = axis.imshow(matrix, origin="lower", interpolation="nearest", cmap="magma", vmin=0, vmax=max(30.0, float(np.nanpercentile(matrix, 99))))
        axis.set_title(f"{candidate.get('candidate_id') or 'BoltzGen design'} · {matrix.shape[0]} residues{sample_label}")
        axis.set_xlabel("Residue index")
        axis.set_ylabel("Residue index")
        fig.colorbar(image, ax=axis, label="Predicted aligned error (Å)")
        st.pyplot(fig, clear_figure=True)
        st.download_button(
            "Download PAE matrix as CSV",
            data=pd.DataFrame(matrix).to_csv(index=False, header=False).encode("utf-8"),
            file_name=f"{candidate.get('candidate_id') or 'boltzgen'}_pae.csv",
            mime="text/csv",
            key=f"boltzgen-pae-download-{run_dir.name}-{selection}",
        )
        st.caption(f"PAE array `{pae_key}` loaded from `{pae_path.name}`{sample_label}.")
    except (ImportError, OSError, ValueError, KeyError) as exc:
        st.warning(f"Could not render the saved PAE matrix: {exc}")


def _bindcraft2_design_key(value: object) -> str:
    name = str(value or "").replace("\\", "/").rsplit("/", 1)[-1].strip()
    for suffix in (
        "_trajectory_monomer",
        "_best_candidate_not_accepted",
        "_not_accepted",
        "_trajectory",
    ):
        if name.casefold().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name.casefold()


def _bindcraft2_trajectory_table(
    dataframe: pd.DataFrame,
    structures: list[tuple[dict, Path]],
) -> pd.DataFrame:
    table = dataframe.copy()
    normalized = {
        re.sub(r"[^a-z0-9]+", "", str(column).casefold()): column
        for column in table.columns
    }
    design_column = normalized.get("design")
    terminated_column = normalized.get("terminated")
    if design_column is None:
        return table

    structures_by_design: dict[str, list[tuple[dict, Path]]] = {}
    for item, path in structures:
        if path.stem.casefold() == "target":
            continue
        design_key = _bindcraft2_design_key(path.stem)
        if design_key:
            structures_by_design.setdefault(design_key, []).append((item, path))

    sequence_column = _sequence_column_name(table.columns) or "Binder_Sequence"
    if sequence_column not in table.columns:
        table[sequence_column] = ""
    design_statuses: list[str] = []
    binder_sequences: list[str] = []
    for _, row in table.iterrows():
        design_key = _bindcraft2_design_key(row.get(design_column))
        matches = structures_by_design.get(design_key, [])
        matches = sorted(
            matches,
            key=lambda pair: (
                {"Accepted": 0, "Not accepted": 1, "Other": 2}.get(_structure_design_status(*pair), 2),
                "monomer" in pair[1].stem.casefold(),
            ),
        )
        design_status = next(
            (_structure_design_status(*pair) for pair in matches if _structure_design_status(*pair) != "Other"),
            "Other",
        )
        terminated = str(row.get(terminated_column) or "").strip().casefold() if terminated_column else ""
        if design_status == "Other" and terminated not in {"", "nan", "none", "-"}:
            design_status = "Not accepted"
        design_statuses.append(design_status)

        sequence = ""
        for _, path in matches:
            if path.suffix.casefold() not in {".cif", ".mmcif"}:
                continue
            try:
                sequence = _mmcif_binder_sequence(str(path), path.stat().st_mtime_ns)
            except OSError:
                sequence = ""
            if sequence:
                break
        binder_sequences.append(sequence)

    table[sequence_column] = binder_sequences
    table["Design Status"] = design_statuses
    return table


def _format_bytes(size: int) -> str:
    value = max(0, int(size))
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} TB"


def _citation_values(value: object, path: str = "") -> list[tuple[str, object]]:
    found: list[tuple[str, object]] = []
    if isinstance(value, dict):
        for key, nested in value.items():
            key_text = str(key)
            key_path = f"{path}.{key_text}" if path else key_text
            if any(marker in key_text.casefold() for marker in ("citation", "doi", "paper", "publication", "reference")):
                if nested not in (None, "", [], {}):
                    found.append((key_path, nested))
            else:
                found.extend(_citation_values(nested, key_path))
    elif isinstance(value, list):
        for index, nested in enumerate(value):
            found.extend(_citation_values(nested, f"{path}[{index}]"))
    return found


def _af3_local_artifact(run_dir: Path, raw_path: object) -> Path | None:
    """Resolve an AF3 artifact path without allowing it to escape the job folder."""
    if not raw_path:
        return None
    candidate = Path(str(raw_path)).expanduser()
    if not candidate.is_absolute():
        candidate = run_dir / candidate
    try:
        candidate = candidate.resolve(strict=True)
        candidate.relative_to(run_dir.resolve())
    except (OSError, ValueError):
        return None
    return candidate if candidate.is_file() else None


def _af3_metric_label(value: object, *, percent: bool = False) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "—"
    if not pd.notna(number):
        return "—"
    if percent and abs(number) <= 1.5:
        number *= 100
    return f"{number:.1f}%" if percent else f"{number:.3f}"


def _af3_optional_metric(value: object) -> float | None:
    if isinstance(value, dict):
        value = value.get("value")
    if isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if pd.notna(number) else None


def _read_a3m_preview(path: Path, max_records: int = 20) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []

    def save_record() -> None:
        if header is not None:
            aligned = "".join(sequence_parts)
            # Lowercase letters and dots encode insertions in A3M files. Omit
            # them for a readable aligned preview while retaining gap markers.
            aligned = "".join(char for char in aligned if not char.islower() and char != ".")
            records.append((header, aligned))

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    save_record()
                    if len(records) >= max_records:
                        break
                header = line[1:]
                sequence_parts = []
            elif header is not None:
                sequence_parts.append(line)
        else:
            save_record()
    return records


def _a3m_preview_html(records: list[tuple[str, str]]) -> str:
    residue_colors = {
        **{residue: "#315fba" for residue in "AVLIMFWY"},
        **{residue: "#24834a" for residue in "STNQ"},
        **{residue: "#bd4141" for residue in "DE"},
        **{residue: "#8b51ad" for residue in "KRH"},
        "C": "#be8a00",
        "G": "#606b66",
        "P": "#cc6b2c",
        "-": "#a0a8a3",
    }
    table_rows = []
    for header, sequence in records:
        letters = "".join(
            f'<span style="color:{residue_colors.get(letter.upper(), "#4d5b54")}">{escape(letter)}</span>'
            for letter in sequence
        )
        table_rows.append(
            "<tr>"
            f'<th title="{escape(header, quote=True)}">{escape(header[:72])}</th>'
            f'<td><code>{letters}</code></td>'
            "</tr>"
        )
    return (
        "<!doctype html><html><head><meta charset=\"utf-8\"><style>"
        "body{margin:0;font:12px system-ui,sans-serif;color:#24332b}"
        ".wrap{max-height:390px;overflow:auto;border:1px solid #dfe7e1;border-radius:6px}"
        "table{border-collapse:collapse;white-space:nowrap;width:max-content;min-width:100%}"
        "th{position:sticky;left:0;background:#f3f8f4;text-align:left;max-width:180px;"
        "overflow:hidden;text-overflow:ellipsis;padding:7px 9px;border-right:1px solid #dfe7e1;"
        "border-bottom:1px solid #e7ece8;font-weight:500}"
        "td{padding:7px 10px;border-bottom:1px solid #e7ece8}"
        "code{font:12px ui-monospace,SFMono-Regular,Consolas,monospace;letter-spacing:.06em}"
        "</style></head><body><div class=\"wrap\"><table><tbody>"
        + "".join(table_rows)
        + "</tbody></table></div></body></html>"
    )


def _plot_plddt_confidence_points(axis: Any, positions: Any, values: Any) -> None:
    """Plot pLDDT points using the same four confidence bands as Mol*."""
    import numpy as np

    scores = np.asarray(values, dtype=float)
    x_values = np.asarray(positions)
    valid = np.isfinite(scores)
    bands = (
        ("≤50", valid & (scores <= 50), "#ff7d45"),
        (">50–70", valid & (scores > 50) & (scores <= 70), "#ffdb13"),
        (">70–90", valid & (scores > 70) & (scores <= 90), "#65cbf3"),
        (">90", valid & (scores > 90), "#0053d6"),
    )
    for label, mask, color in bands:
        axis.scatter(x_values[mask], scores[mask], s=9, color=color, label=label, linewidths=0)
    axis.legend(
        title="pLDDT", loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=4,
        fontsize=8, title_fontsize=8, frameon=False,
    )


_PY2DMOL_PAE_GREEN_COLORS = ("#246b3c", "#78b889", "#c6e3cd", "#f4faf5")


def _pyplot_pae_green_colormap() -> Any:
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("studio_pae_green", _PY2DMOL_PAE_GREEN_COLORS)












@st.cache_data(show_spinner="Reading prediction confidence outputs…", ttl=3600, max_entries=16)
def _cached_cofolding_confidence_data(
    structure_file: str,
    structure_mtime_ns: int,
    structure_size: int,
    pae_file: str,
    pae_mtime_ns: int,
    pae_size: int,
) -> dict:
    """Read per-model confidence, selecting BoltzGen's exported refolding sample."""
    del structure_mtime_ns, structure_size, pae_mtime_ns, pae_size
    structure = Path(structure_file)
    files = cofolding_results.confidence_files(structure)
    if pae_file and Path(pae_file) not in files:
        files.insert(0, Path(pae_file))
    return cofolding_results.read_confidence(structure, tuple(files))


@st.cache_data(show_spinner="Preparing structure viewer…", ttl=3600, max_entries=12)
def _cached_cofolding_py2dmol_html(
    structure_file: str,
    structure_mtime_ns: int,
    structure_size: int,
    pae_file: str,
    pae_mtime_ns: int,
    pae_size: int,
) -> str:
    """Build the shared py2Dmol viewer with the saved PAE map when available."""
    from py2Dmol import view

    confidence = _cached_cofolding_confidence_data(
        structure_file, structure_mtime_ns, structure_size, pae_file, pae_mtime_ns, pae_size,
    )
    pae_preview = confidence["pae_preview"]
    viewer = view(
        size=(560, 560), controls=True, box=True, color="deepmind",
        style="richardson", bg="white", heatmap=pae_preview is not None, heatmap_size=560,
    )
    pae_maps = (
        {"pae": {
            "data": pae_preview.tolist(), "vmin": 0.0, "vmax": 30.0,
            "colors": _PY2DMOL_PAE_GREEN_COLORS,
        }}
        if pae_preview is not None else None
    )
    viewer.add_pdb(str(structure_file), name="prediction", maps=pae_maps)
    return viewer.to_html(title="Cofolding prediction viewer")


def _cofolding_candidate_pae_path(run_dir: Path, candidate: dict, structure_path: Path | None) -> Path | None:
    _, files = cofolding_results.candidate_files(run_dir, candidate)
    for path in files:
        name = path.name.casefold()
        if any(token in name for token in ("pae", "full_data_sample", "_confidences.json", "_confidence.npz", "scores.model_idx_", "_scores_rank_")):
            return path
    return next((path for path in files if path.suffix == ".npz"), None)


def _boltz2_msa_files(run_dir: Path, raw_metadata: dict) -> list[tuple[str, Path]]:
    """Find staged and Boltz-generated protein A3Ms for this batch item."""
    import yaml

    item_id = str(raw_metadata.get("batch_item_id") or "")
    if not item_id or not re.fullmatch(r"[A-Za-z0-9_-]{1,100}", item_id):
        return []
    input_path = _af3_local_artifact(
        run_dir, run_dir / "artifacts" / "raw" / "optimized_prediction" / "inputs" / f"{item_id}.yaml",
    )
    if input_path is None:
        return []
    try:
        payload = yaml.safe_load(input_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError):
        return []
    sequences = payload.get("sequences") if isinstance(payload, dict) else None
    if not isinstance(sequences, list):
        return []
    found: list[tuple[str, Path]] = []
    seen: set[tuple[str, Path]] = set()
    reference_chains: set[str] = set()
    staged_chains: set[str] = set()
    cached_generated_files: set[tuple[str, str]] = set()
    mapping_path = _af3_local_artifact(
        run_dir, run_dir / "artifacts" / "raw" / "optimized_prediction" / "input_mapping.json",
    )
    if mapping_path is not None:
        try:
            mapping_payload = json.loads(mapping_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            mapping_payload = {}
        mapping_items = mapping_payload.get("items") if isinstance(mapping_payload, dict) else None
        for mapping_item in mapping_items if isinstance(mapping_items, list) else []:
            if not isinstance(mapping_item, dict) or str(mapping_item.get("id") or "") != item_id:
                continue
            input_details = mapping_item.get("input_details") if isinstance(mapping_item.get("input_details"), dict) else {}
            reference_rows = input_details.get("reference_msa_chains")
            for reference_row in reference_rows if isinstance(reference_rows, list) else []:
                if isinstance(reference_row, dict) and reference_row.get("chain_id"):
                    reference_chains.update(str(reference_row["chain_id"]).split(","))
            generated_rows = input_details.get("generated_msa_cached_chains")
            for generated_row in generated_rows if isinstance(generated_rows, list) else []:
                if isinstance(generated_row, dict) and generated_row.get("chain_id") and generated_row.get("file_name"):
                    cached_generated_files.add((str(generated_row["chain_id"]), str(generated_row["file_name"])))
            break
    protein_sequences: list[tuple[str, str]] = []
    for row in sequences:
        protein = row.get("protein") if isinstance(row, dict) else None
        if not isinstance(protein, dict):
            continue
        chain_ids = protein.get("id") if isinstance(protein.get("id"), list) else [protein.get("id")]
        sequence = str(protein.get("sequence") or "")
        for chain_id in chain_ids:
            protein_sequences.append((str(chain_id or "Protein"), sequence))
        raw_msa = protein.get("msa")
        if not isinstance(raw_msa, str) or not raw_msa.strip() or raw_msa.casefold() in {"empty", "none"}:
            continue
        msa_value = raw_msa.strip()
        if msa_value.startswith("/work/"):
            msa_value = msa_value.removeprefix("/work/")
        path = _af3_local_artifact(run_dir, msa_value)
        if path is None:
            continue
        for chain_id in chain_ids:
            chain_label = str(chain_id or "Protein")
            staged_chains.add(chain_label)
            source_label = "reference folder" if chain_label in reference_chains else "staged input"
            key = (f"{chain_label} · {source_label}", path)
            if key not in seen:
                found.append(key)
                seen.add(key)

    # Boltz can generate MSAs during prediction (for example through the
    # ColabFold server) without recording an input path in the YAML. Those
    # alignments live beside the prediction outputs under boltz_results_<id>.
    output_root = run_dir / "artifacts" / "raw" / "optimized_prediction" / "output"
    generated_dirs = output_root.glob(f"**/boltz_results_{item_id}/msa")
    sequence_by_chain: dict[str, str] = {}
    for chain_id, sequence in protein_sequences:
        sequence_by_chain[chain_id] = "".join(char for char in sequence.upper() if char.isalpha())
    generated_by_chain: dict[str, list[Path]] = {}
    for msa_dir in generated_dirs:
        for path in sorted(msa_dir.rglob("*.a3m")):
            resolved_path = _af3_local_artifact(run_dir, path)
            if resolved_path is None:
                continue
            preview = _read_a3m_preview(resolved_path, max_records=1)
            if not preview:
                continue
            query_sequence = "".join(char for char in preview[0][1].upper() if char.isalpha())
            matching_chains = [
                chain_id for chain_id, sequence in sequence_by_chain.items()
                if sequence and sequence == query_sequence
            ]
            if not matching_chains and len(sequence_by_chain) == 1:
                matching_chains = list(sequence_by_chain)
            for chain_id in matching_chains:
                if chain_id not in staged_chains:
                    generated_by_chain.setdefault(chain_id, []).append(resolved_path)

    def generated_alignment_priority(chain_id: str, path: Path) -> tuple[int, str]:
        if (chain_id, path.name) in cached_generated_files:
            return -1, path.name.casefold()
        filename = path.name.casefold()
        if filename == "uniref.a3m":
            priority = 0
        elif filename == "pair.a3m":
            priority = 1
        elif filename.startswith("bfd."):
            priority = 2
        else:
            priority = 3
        return priority, filename

    for chain_id, paths in sorted(generated_by_chain.items()):
        path = min(paths, key=lambda candidate: generated_alignment_priority(chain_id, candidate))
        if path.name.casefold() == "uniref.a3m":
            source_label = "UniRef A3M"
        elif path.name.casefold() == "pair.a3m":
            source_label = "paired A3M"
        else:
            source_label = path.name.removesuffix(".a3m")
        cache_note = " · saved for reuse" if (chain_id, path.name) in cached_generated_files else ""
        found.append((f"{chain_id} · generated · {source_label}{cache_note}", path))
    return found


def _render_cofolding_prediction_results(
    run_id: str,
    run_dir: Path,
    candidates: list[dict],
    input_payload: dict,
    result: dict,
) -> None:
    """Shared Boltz-style results, with native artifact adapters for each engine."""
    import numpy as np

    params = input_payload.get("params") or input_payload.get("parameters") or {}
    run_metrics = result.get("metrics") or {}
    engine = str(params.get("engine") or run_metrics.get("engine") or "")
    model_label = str(run_metrics.get("model_label") or params.get("model_label") or "")
    if not model_label and any(str(item.get("source_tool") or "") == "af3_jax_cofolding" for item in candidates):
        model_name = str(params.get("model_name") or run_metrics.get("model_name") or "alphafold3")
        model_label = str(af3_jax_cofolding.MODEL_PROFILES.get(model_name, {}).get("label") or model_name)
    if not model_label:
        if any(str(item.get("source_tool") or "") == "af2_initial_guess" for item in candidates):
            model_label = "AlphaFold2 Initial Guess"
    if not model_label:
        profile = optimized_prediction.ENGINE_PROFILES.get(engine, {})
        model_label = str(profile.get("label") or ("ESMFold2" if "model_variant" in params else "ColabFold 2.0 preview")).removesuffix(" optimized")
    candidates = cofolding_results.ranked_candidates([
        cofolding_results.enrich_candidate(run_dir, candidate) for candidate in candidates if isinstance(candidate, dict)
    ])
    st.subheader(f"{model_label} predictions")
    comparison_report = cofolding_results.local_artifact(run_dir, (result.get("outputs") or {}).get("validation_report"))
    if comparison_report:
        st.download_button("Download design comparisons", comparison_report.read_bytes(),
                           file_name=f"{run_id}_design_comparisons.csv", mime="text/csv", key=f"validation-report-{run_id}")
    if not candidates:
        st.info("No predictions have been saved yet. Refresh after the job finishes, or check Logs for progress.")
    else:
        complex_ids = list(dict.fromkeys(str((item.get("raw_metadata") or {}).get("batch_item_id") or "") for item in candidates))
        if len(complex_ids) > 1:
            complex_labels = []
            for index, complex_id in enumerate(complex_ids):
                item = next(c for c in candidates if str((c.get("raw_metadata") or {}).get("batch_item_id") or "") == complex_id)
                reference = (item.get("raw_metadata") or {}).get("validation_reference") or {}
                complex_labels.append(f"{Path(reference['source_run_dir']).name} · {reference['candidate_id']}" if reference else f"Complex {index + 1}")
            complex_index = st.selectbox("Design" if comparison_report else "Complex", range(len(complex_ids)), format_func=lambda index: complex_labels[index], key=f"cofolding-complex-{run_id}")
            candidates = [item for item in candidates if str((item.get("raw_metadata") or {}).get("batch_item_id") or "") == complex_ids[complex_index]]

        def structure_path_for(candidate: dict) -> Path | None:
            return cofolding_results.structure_for(run_dir, candidate)

        def prediction_label(index: int) -> str:
            candidate = candidates[index]
            item_metrics = candidate.get("metrics") or {}
            raw = candidate.get("raw_metadata") or {}
            label = f"Model {cofolding_results.model_index(candidate, index)}"
            seeds = {str((item.get("metrics") or {}).get("seed") or (item.get("raw_metadata") or {}).get("seed") or "") for item in candidates}
            if len(seeds) > 1:
                label += f" · seed {item_metrics.get('seed') or raw.get('seed')}"
            if index == 0 and any(cofolding_results.number(item_metrics.get(key)) is not None for key in ("ranking_score", "confidence_score", "iptm", "ptm", "plddt_mean")):
                label += " (best)"
            return label

        selected_index = st.selectbox(
            "Model", range(len(candidates)), format_func=prediction_label,
            key=f"cofolding-prediction-{run_id}",
        )
        candidate = candidates[selected_index]
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        raw_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        confidence = raw_metadata.get("confidence_outputs") if isinstance(raw_metadata.get("confidence_outputs"), dict) else {}
        structure_path = structure_path_for(candidate)
        pae_path = _cofolding_candidate_pae_path(run_dir, candidate, structure_path)

        def confidence_metric(*names: str) -> object:
            normalized = {str(key).casefold(): value for key, value in confidence.items()}
            normalized_metrics = {str(key).casefold(): value for key, value in metrics.items()}
            for name in names:
                key = name.casefold()
                value = normalized_metrics.get(key)
                if value is None:
                    value = normalized.get(key)
                if value is not None:
                    return value
            return None

        top_cards: list[tuple[str, object, bool]] = [
            ("Mean pLDDT", metrics.get("plddt_mean"), True),
            ("Confidence score", confidence_metric("confidence_score", "ranking_score"), False),
            ("pTM", confidence_metric("ptm"), False),
            ("ipTM", confidence_metric("iptm"), False),
        ]
        affinity_value = confidence_metric("affinity_pred_value", "affinity_prediction", "pic50_binding_affinity")
        affinity_probability = confidence_metric("affinity_probability_binary", "binding_probability")
        if affinity_value is not None:
            top_cards.append(("Affinity prediction", affinity_value, False))
        if affinity_probability is not None:
            top_cards.append(("Binding probability", affinity_probability, True))
        top_cards = [card for card in top_cards if card[1] is not None]
        metric_columns = st.columns(max(1, len(top_cards)))
        for column, (label, value, percent) in zip(metric_columns, top_cards):
            column.metric(label, _af3_metric_label(value, percent=percent))

        affinity_rows = []
        for key, value in confidence.items():
            normalized_key = str(key).casefold()
            if not (normalized_key.startswith("affinity_") or normalized_key == "binding_probability"):
                continue
            number = _af3_optional_metric(value)
            if number is None:
                continue
            if normalized_key.startswith("affinity_pred_value"):
                label = "Affinity prediction" + normalized_key.removeprefix("affinity_pred_value")
            elif normalized_key.startswith("affinity_probability_binary"):
                label = "Binding probability" + normalized_key.removeprefix("affinity_probability_binary")
            else:
                label = normalized_key.replace("_", " ").title()
            affinity_rows.append({"Metric": label, "Value": _af3_metric_label(number, percent="probability" in normalized_key)})
        if len(affinity_rows) > 2:
            with st.expander("Additional affinity outputs", expanded=False):
                st.dataframe(pd.DataFrame(affinity_rows), hide_index=True, use_container_width=True)

        def confidence_data_for(item: dict) -> tuple[Path | None, Path | None, dict | None, str | None]:
            item_structure = structure_path_for(item)
            item_pae = _cofolding_candidate_pae_path(run_dir, item, item_structure)
            if item_structure is None:
                return item_structure, item_pae, None, "The structure file is unavailable."
            try:
                item_structure_stat = item_structure.stat()
                item_pae_stat = item_pae.stat() if item_pae is not None else None
                data = _cached_cofolding_confidence_data(
                    str(item_structure), item_structure_stat.st_mtime_ns, item_structure_stat.st_size,
                    str(item_pae) if item_pae else "", item_pae_stat.st_mtime_ns if item_pae_stat else 0,
                    item_pae_stat.st_size if item_pae_stat else 0,
                )
                return item_structure, item_pae, data, None
            except (ImportError, OSError, ValueError, KeyError) as exc:
                return item_structure, item_pae, None, str(exc)

        confidence_data = None
        confidence_error = None
        selected_confidence = confidence_data_for(candidate)
        if selected_confidence[0] == structure_path and selected_confidence[1] == pae_path:
            confidence_data = selected_confidence[2]
            confidence_error = selected_confidence[3]

        render_comparison(run_dir, candidate, key=f"validation-{run_id}-{candidate.get('candidate_id')}")
        st.subheader("Display structures")
        if structure_path is None:
            st.info("The selected prediction has no available PDB or mmCIF structure file.")
        elif structure_path.stat().st_size > 50 * 1024 * 1024:
            st.info("This structure is over 50 MB; download it from the Files tab to inspect it locally.")
        else:
            try:
                structure_stat = structure_path.stat()
                pae_stat = pae_path.stat() if pae_path is not None else None
                viewer_html = _cached_cofolding_py2dmol_html(
                    str(structure_path), structure_stat.st_mtime_ns, structure_stat.st_size,
                    str(pae_path) if pae_path else "", pae_stat.st_mtime_ns if pae_stat else 0,
                    pae_stat.st_size if pae_stat else 0,
                )
                import streamlit.components.v1 as components

                components.html(viewer_html, height=660, scrolling=False)
                st.caption(
                    "Interactive py2Dmol viewer with per-residue pLDDT coloring. "
                    + ("The PAE map is available in the viewer and confidence plots." if confidence_data is not None and confidence_data.get("pae") is not None else "No PAE array was saved for this prediction.")
                )
            except Exception as exc:
                st.error(f"Could not build the interactive structure viewer: {exc}")

        candidate_keys = {
            str(key).casefold()
            for item in candidates
            for source in (
                item.get("metrics") if isinstance(item.get("metrics"), dict) else {},
                (item.get("raw_metadata") or {}).get("confidence_outputs")
                if isinstance(item.get("raw_metadata"), dict) else {},
            )
            for key in source
        }
        display_keys = [
            ("confidence_score", "Confidence score", False), ("ranking_score", "Ranking score", False),
            ("ptm", "pTM", False), ("iptm", "ipTM", False),
            ("complex_plddt", "Complex pLDDT", True), ("ligand_iptm", "Ligand ipTM", False),
            ("protein_iptm", "Protein ipTM", False), ("complex_pde", "Complex PDE", False),
            ("complex_ipde", "Complex iPDE", False), ("ipsae", "ipSAE", False),
        ]
        known_metric_keys = {key for key, _label, _percent in display_keys} | {"plddt_mean", "prediction_sample"}
        extra_metric_keys = []
        for key in sorted(candidate_keys - known_metric_keys):
            if key.startswith("affinity_") or key == "binding_probability":
                continue
            for item in candidates:
                item_metrics = item.get("metrics") if isinstance(item.get("metrics"), dict) else {}
                item_raw = item.get("raw_metadata") if isinstance(item.get("raw_metadata"), dict) else {}
                item_confidence = item_raw.get("confidence_outputs") if isinstance(item_raw.get("confidence_outputs"), dict) else {}
                value = {str(name).casefold(): nested for name, nested in item_metrics.items()}.get(key)
                if value is None:
                    value = {str(name).casefold(): nested for name, nested in item_confidence.items()}.get(key)
                if _af3_optional_metric(value) is not None:
                    extra_metric_keys.append(key)
                    break
        affinity_table_keys = sorted(
            key for key in candidate_keys
            if key.startswith("affinity_") or key == "binding_probability"
        )
        boltz_params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
        boltz_tasks = boltz_params.get("tasks") if isinstance(boltz_params.get("tasks"), list) else []
        first_task = boltz_tasks[0] if boltz_tasks and isinstance(boltz_tasks[0], dict) else {}
        input_sequence = str(first_task.get("protein_sequences") or "").strip()
        table_rows = []
        for index, item in enumerate(candidates):
            item_metrics = item.get("metrics") if isinstance(item.get("metrics"), dict) else {}
            item_raw = item.get("raw_metadata") if isinstance(item.get("raw_metadata"), dict) else {}
            item_confidence = item_raw.get("confidence_outputs") if isinstance(item_raw.get("confidence_outputs"), dict) else {}
            row: dict[str, object] = {}
            item_proteins = item_raw.get("protein_chains") or item_raw.get("chain_sequences") or {}
            if isinstance(item_proteins, dict):
                sequence = ":".join(str(value) for value in item_proteins.values())
            else:
                sequence = ":".join(str(chain.get("sequence") or "") for chain in item_proteins if isinstance(chain, dict))
            if sequence or input_sequence:
                row["Sequence"] = sequence or input_sequence
            row.update({
                "Model": prediction_label(index),
                "Mean pLDDT": _af3_metric_label(item_metrics.get("plddt_mean"), percent=True),
            })
            sources = ({str(key).casefold(): value for key, value in item_metrics.items()},
                       {str(key).casefold(): value for key, value in item_confidence.items()})
            for key, label, percent in display_keys:
                if key not in candidate_keys or (key == "confidence_score" and "ranking_score" in candidate_keys):
                    continue
                value = next((source[key] for source in sources if key in source and source[key] is not None), None)
                row[label] = _af3_metric_label(value, percent=percent)
            for key in extra_metric_keys:
                value = next((source[key] for source in sources if key in source and source[key] is not None), None)
                label = key.replace("_", " ").title()
                row[label] = _af3_metric_label(value)
            for key in affinity_table_keys:
                if key.startswith("affinity_pred_value"):
                    label = "Affinity prediction" + key.removeprefix("affinity_pred_value")
                elif key.startswith("affinity_probability_binary"):
                    label = "Binding probability" + key.removeprefix("affinity_probability_binary")
                else:
                    label = key.replace("_", " ").title()
                value = next((source[key] for source in sources if key in source and source[key] is not None), None)
                row[label] = _af3_metric_label(value, percent="probability" in key)
            row["Structure"] = str(item.get("complex_pdb") or "—")
            table_rows.append(row)
        st.subheader("Quality and affinity metrics" if affinity_table_keys else "Quality metrics")
        table_frame = pd.DataFrame(table_rows)
        table_scroll_height = min(520, max(170, 34 + min(len(table_frame), 12) * 32))
        table_frame_height = table_scroll_height + (86 if "Sequence" in table_frame.columns else 66)
        download_column, _ = st.columns([1, 5])
        download_column.download_button(
            "Download metrics CSV",
            data=table_frame.to_csv(index=False).encode("utf-8"),
            file_name=f"{run_id}_prediction_metrics.csv",
            mime="text/csv",
            key=f"boltz2-metrics-download-{run_id}",
        )
        import streamlit.components.v1 as components

        components.html(
            _render_tool_output_table_css() + _render_tool_output_table(
                table_frame,
                metric_color_rules={
                    "Confidence score": ("gradient_high", None),
                    "Ranking score": ("gradient_high", None),
                    "pTM": ("threshold_high", 0.8),
                    "ipTM": ("threshold_high", 0.8),
                },
                scroll_height=table_scroll_height,
            ),
            height=table_frame_height,
            scrolling=False,
        )
        st.caption("Greener confidence scores indicate higher confidence. pTM and ipTM are green at or above 0.8.")
        with st.expander("Raw metrics table", expanded=False):
            st.dataframe(table_frame, hide_index=True, use_container_width=True)

        chain_details = (confidence_data or {}).get("details", {})
        if chain_details:
            with st.expander("Chain and interface metrics", expanded=False):
                for key, value in chain_details.items():
                    if key == "chain_ids":
                        continue
                    st.markdown(f"**{key.replace('_', ' ').title()}**")
                    if isinstance(value, list):
                        try:
                            array = np.asarray(value)
                            while array.ndim > 2 and array.shape[0] == 1:
                                array = array[0]
                            if array.ndim in {1, 2} and np.issubdtype(array.dtype, np.number):
                                frame = pd.DataFrame(array)
                                labels = chain_details.get("chain_ids") or []
                                if len(labels) == len(frame):
                                    frame.index = labels
                                    if len(labels) == len(frame.columns):
                                        frame.columns = labels
                                st.dataframe(frame, use_container_width=True)
                                continue
                        except (ValueError, TypeError):
                            pass
                    st.json(value, expanded=False)

        st.subheader("Confidence plots · all models")
        any_confidence_plot = False
        try:
            import matplotlib

            matplotlib.use("Agg", force=True)
            import matplotlib.pyplot as plt

            with st.container(height=760, border=False):
                for index, item in enumerate(candidates):
                    item_structure, _item_pae, item_confidence, item_error = (
                        selected_confidence if index == selected_index else confidence_data_for(item)
                    )
                    st.markdown(f"**{prediction_label(index)}**")
                    if item_structure is None or item_confidence is None:
                        st.caption(f"Confidence figures unavailable: {item_error or 'prediction files are missing.'}")
                        continue
                    plddt = item_confidence["plddt"]
                    pae_matrix = item_confidence["pae"]
                    figure = plt.figure(figsize=(16, 4.5), constrained_layout=True)
                    grid = figure.add_gridspec(1, 2, width_ratios=(1, 1))
                    plddt_axis = figure.add_subplot(grid[0, 0])
                    pae_axis = figure.add_subplot(grid[0, 1])
                    if plddt.size:
                        stride = max(1, int(np.ceil(plddt.size / 5000)))
                        positions = np.arange(0, plddt.size, stride)
                        _plot_plddt_confidence_points(plddt_axis, positions, plddt[::stride])
                        plddt_axis.set(xlabel=item_confidence.get("plddt_kind", "Atom"), ylabel="pLDDT", ylim=(0, 100), title="Predicted pLDDT")
                        for threshold in (50, 70, 90):
                            plddt_axis.axhline(threshold, linestyle="--", linewidth=0.7, color="#94a39a", alpha=0.7)
                        any_confidence_plot = True
                    else:
                        plddt_axis.axis("off")
                        plddt_axis.text(0.5, 0.5, "Per-atom pLDDT values were not found.", ha="center", va="center")

                    if pae_matrix is not None:
                        display_matrix = item_confidence["pae_preview"]
                        heatmap = pae_axis.imshow(
                            display_matrix, cmap=_pyplot_pae_green_colormap(), vmin=0, vmax=30,
                            interpolation="nearest", aspect="auto",
                        )
                        pae_axis.set_box_aspect(1)
                        colorbar = figure.colorbar(heatmap, ax=pae_axis, fraction=0.046, pad=0.04)
                        colorbar.set_label("PAE (Å)", fontsize=8)
                        sampled = item_confidence["pae_sampled"]
                        pae_axis.set_xlabel("Scored residue" + (" (sampled)" if sampled else ""))
                        pae_axis.set_ylabel("Aligned residue" + (" (sampled)" if sampled else ""))
                        pae_axis.set_title("Predicted Aligned Error (PAE)" + (" · sampled preview" if sampled else ""))
                        any_confidence_plot = True
                    else:
                        pae_axis.axis("off")
                        pae_axis.text(0.5, 0.5, "No PAE matrix was saved for this prediction.", ha="center", va="center")
                    st.pyplot(figure, clear_figure=True)
                    plt.close(figure)
                    st.caption("PAE matrices are square and height-matched to pLDDT; darker green indicates lower error on the fixed 0–30 Å scale.")
        except (ImportError, OSError, ValueError, TypeError) as exc:
            st.warning(f"Could not render confidence plots: {exc}")
        if not any_confidence_plot and not confidence_error:
            st.info("No pLDDT or PAE figures were available for the saved predictions.")

        msa_files = cofolding_results.manifest_alignment_files(run_dir, candidate) or (_boltz2_msa_files(run_dir, raw_metadata) if engine == "boltz2" else cofolding_results.alignment_files(run_dir, candidate))
        with st.expander("A3M Sequence Alignment Viewer", expanded=False):
            if not msa_files:
                if engine == "esmfold":
                    st.info("ESMFold predicts from the submitted protein sequences and does not use an MSA.")
                else:
                    st.info(cofolding_results.alignment_status(run_dir, candidate))
            else:
                if len(msa_files) == 1:
                    selected_msa = 0
                    st.caption(msa_files[0][0])
                else:
                    selected_msa = st.selectbox(
                        "Alignment", range(len(msa_files)),
                        format_func=lambda index: msa_files[index][0],
                        key=f"boltz2-msa-chain-{run_id}-{selected_index}",
                    )
                alignment_label, msa_path = msa_files[selected_msa]
                if " · generated · " in alignment_label:
                    if " · saved for reuse" in alignment_label:
                        st.caption("Generated during prediction and saved in the shared reference MSA cache for later jobs.")
                    else:
                        st.caption("Showing one A3M preview from the MSA server. Boltz also uses environmental search results and, for multichain inputs, paired results internally; these remain in the job outputs.")
                elif " · reference folder" in alignment_label:
                    st.caption("Loaded from the shared reference folder; a copy is saved with this job's inputs.")
                try:
                    alignment_rows = _read_a3m_preview(msa_path, max_records=20)
                    if alignment_rows:
                        import streamlit.components.v1 as components

                        components.html(_a3m_preview_html(alignment_rows), height=400, scrolling=False)
                        st.caption(f"{alignment_label} · showing up to {len(alignment_rows)} alignment rows.")
                        frequency_rows = _read_a3m_preview(msa_path, max_records=8192)
                        aligned_sequences = [sequence for _, sequence in frequency_rows if sequence]
                        if aligned_sequences:
                            import matplotlib

                            matplotlib.use("Agg", force=True)
                            import matplotlib.pyplot as plt
                            from matplotlib.patches import Patch
                            from collections import Counter

                            alignment_length = min(map(len, aligned_sequences))
                            residue_colors = {
                                **{residue: "#315fba" for residue in "AVLIMFWY"},
                                **{residue: "#24834a" for residue in "STNQ"},
                                **{residue: "#bd4141" for residue in "DE"},
                                **{residue: "#8b51ad" for residue in "KRH"},
                                "C": "#be8a00", "G": "#606b66", "P": "#cc6b2c",
                                "-": "#a0a8a3",
                            }
                            residue_order = list("ACDEFGHIKLMNPQRSTVWY-")
                            residue_frequencies = {residue: np.zeros(alignment_length) for residue in residue_order}
                            for position in range(alignment_length):
                                counts = Counter(sequence[position].upper() for sequence in aligned_sequences)
                                for residue, count in counts.items():
                                    if residue not in residue_frequencies:
                                        residue_frequencies[residue] = np.zeros(alignment_length)
                                        residue_order.append(residue)
                                        residue_colors.setdefault(residue, "#4d5b54")
                                    residue_frequencies[residue][position] = count / len(aligned_sequences)
                            figure, axis = plt.subplots(figsize=(10, 2.8), constrained_layout=True)
                            positions = np.arange(alignment_length)
                            bottom = np.zeros(alignment_length)
                            for residue in residue_order:
                                values = residue_frequencies[residue]
                                if not np.any(values):
                                    continue
                                axis.bar(
                                    positions, values, bottom=bottom,
                                    color=residue_colors[residue], width=0.95, linewidth=0,
                                )
                                bottom += values
                            axis.set(
                                xlim=(-0.5, max(0.5, alignment_length - 0.5)), ylim=(0, 1),
                                xlabel="Sequence position", ylabel="Fraction of rows",
                                title="Residue frequencies by position",
                            )
                            axis.set_xticks(np.arange(0, alignment_length, max(1, alignment_length // 20)))
                            axis.set_yticks((0, 0.25, 0.5, 0.75, 1), labels=("0%", "25%", "50%", "75%", "100%"))
                            axis.grid(axis="y", alpha=0.18)
                            axis.set_axisbelow(True)
                            legend_items = [
                                ("Hydrophobic", "#315fba"), ("Polar", "#24834a"),
                                ("Acidic", "#bd4141"), ("Basic", "#8b51ad"),
                                ("Cys", "#be8a00"), ("Gly", "#606b66"),
                                ("Pro", "#cc6b2c"), ("Gap", "#a0a8a3"),
                            ]
                            axis.legend(
                                handles=[Patch(facecolor=color, label=label) for label, color in legend_items],
                                loc="upper center", bbox_to_anchor=(0.5, 1.28), ncol=4,
                                frameon=False, fontsize=7,
                            )
                            st.pyplot(figure, clear_figure=True)
                            plt.close(figure)
                            st.caption(
                                f"Each stacked column shows residue frequencies across {len(aligned_sequences)} A3M rows "
                                "(query included, capped at 8,192); gray is a gap. Colors group amino-acid types."
                            )
                    else:
                        st.info("The staged A3M file contains no alignment rows.")
                except (OSError, UnicodeError, ValueError) as exc:
                    st.warning(f"Could not read this chain's A3M alignment: {exc}")

        if structure_path is not None:
            try:
                st.download_button(
                    "Download predicted structure", data=structure_path.read_bytes(),
                    file_name=structure_path.name,
                    mime="chemical/x-pdb" if structure_path.suffix.casefold() in {".pdb", ".ent"} else "chemical/x-mmcif",
                    key=f"boltz2-structure-download-{run_id}-{selected_index}",
                )
            except OSError as exc:
                st.caption(f"Could not prepare the structure download: {exc}")
        if confidence_data is not None and confidence_data.get("pae") is not None:
            pae_matrix = confidence_data["pae"]
            if max(pae_matrix.shape) <= 2048:
                st.download_button(
                    "Download full PAE CSV",
                    data=pd.DataFrame(pae_matrix).to_csv(index=False, header=False).encode("utf-8"),
                    file_name=f"{candidate.get('candidate_id') or 'boltz2_prediction'}_pae.csv",
                    mime="text/csv", key=f"boltz2-pae-download-{run_id}-{selected_index}",
                )
            else:
                st.caption("The full PAE archive is available in the Files tab; CSV export is limited to 2,048 residues.")

    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    settings = [
        ("Engine", params.get("engine") or (result.get("metrics") or {}).get("engine")),
        ("Optimization mode", params.get("optimization_mode")), ("Seed", params.get("seed") or params.get("seeds")),
        ("Diffusion samples", params.get("samples") or params.get("num_diffusion_samples")), ("Recycles", params.get("num_recycles") or params.get("num_loops")),
        ("Engine options", params.get("engine_options")),
    ]
    with st.expander("Run settings and inputs", expanded=False):
        for label, value in settings:
            if value not in (None, "", [], {}):
                st.markdown(f"**{label}**  \n{escape(str(value))}")
        input_path = next(iter(sorted(run_dir.glob("artifacts/raw/*/input_mapping.json"))), run_dir / "missing-input-mapping.json")
        try:
            mapping = json.loads(input_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            mapping = None
        if mapping:
            st.json(mapping, expanded=False)
    if result.get("metrics"):
        with st.expander("Run summary", expanded=False):
            st.json(result["metrics"], expanded=False)


def _render_job_details() -> None:
    reference = str(st.session_state.get("selected-job-reference") or "")
    if not reference:
        st.info("Choose a job from My results to open its saved outputs.")
        if st.button("Back to My results", key="job-details-empty-back"):
            st.session_state["page"] = "My runs"
            st.rerun()
        return

    back_col, refresh_col = st.columns([0.85, 0.15], vertical_alignment="center")
    if back_col.button("← Back to My results", key="job-details-back"):
        st.session_state["page"] = "My runs"
        st.rerun()
    if refresh_col.button("Refresh", icon=":material/refresh:", key="job-details-refresh"):
        st.rerun()

    if reference.startswith("cluster:"):
        from mn_molecular_design_studio.cluster.ui import render_job
        render_job(reference, CURRENT_USER, str(st.session_state.get("studio-session-token") or "development-mode"), _post_local_api)
        return

    try:
        details = _inspect_current_user_job(reference)
        result_payload = job_results(reference, include_candidates=True)
    except (AutomationError, OSError, ValueError) as exc:
        st.error(f"Could not open this job: {exc}")
        return

    job = details["job"]
    metadata = details.get("metadata") or {}
    input_payload = details.get("input") or {}
    result = details.get("result") or {}
    run_id = str(job.get("run_id") or "job")
    surface_titles = {"pesto": "PeSTo", "scannet": "ScanNet", "masif_seed": "MaSIF-seed", "surf2spot": "Surf2Spot", "surface_analysis": "Protein Surface Analysis"}
    raw_title = str(job.get("tool") or job.get("job_type") or run_id)
    title = str(job.get("campaign_name") or surface_titles.get(raw_title.casefold()) or raw_title)
    st.title(title)
    st.caption(f"Tool: {display_tool_label(job)} · Pipeline: {job.get('pipeline_name') or '—'}")
    if metadata.get("demo_fixture"):
        st.info("Imported demo fixture · These results were copied from an existing ZIP and were not run by Molecular Design Studio.")
    submitted = _job_submitted_label(job.get("created_at"))
    st.caption(f"{submitted} · Run ID: `{run_id}`")
    project_id = None
    try:
        catalog_matches = search_jobs(
            query=run_id,
            task_group=str(job.get("task_group") or "") or None,
            limit=10,
            include_hidden=True,
        )["items"]
        catalog_job = next((item for item in catalog_matches if str(item.get("run_id")) == run_id), None)
        if catalog_job:
            project_id = catalog_job.get("project_id")
    except (OSError, sqlite3.Error, ValueError):
        pass
    project_id = project_id or job.get("project_id") or metadata.get("project_id")
    if not project_id:
        try:
            studio_marker = json.loads((Path(job["run_dir"]) / "molecular_design_studio.json").read_text(encoding="utf-8"))
        except (KeyError, OSError, json.JSONDecodeError):
            studio_marker = {}
        if isinstance(studio_marker, dict):
            project_id = studio_marker.get("project_id")
    project_name = next(
        (str(project.get("name")) for project in _studio_projects() if project.get("id") == project_id),
        "No project",
    )
    output_summary = ("Candidates", str(details.get("candidate_count", 0)))
    job_metrics = result.get("metrics") or {}
    if raw_title == "surface_analysis":
        output_summary = ("Methods completed", f"{job_metrics.get('successful_tools', 0)} / {job_metrics.get('tool_count', 4)}")
    elif raw_title in surface_titles:
        output_summary = ("Output files", str(job_metrics.get("artifact_count", 0)))
    elif raw_title == "visual_pipeline" and job_metrics.get("surface_analysis_count") and not details.get("candidate_count"):
        output_summary = ("Surface analyses", str(job_metrics["surface_analysis_count"]))
    summary_values = (
        ("Status", str(job.get("status", "unknown")).replace("_", " ").title()),
        ("Tool", display_tool_label(job)),
        ("Pipeline", str(job.get("pipeline_name") or "—")),
        output_summary,
        ("Time elapsed", _job_elapsed_label(job)),
        ("Project", project_name),
    )
    summary_items = "".join(
        '<div style="display:flex;align-items:center;gap:.38rem;padding:.38rem .62rem;'
        'border:1px solid #e2e8e4;border-radius:.55rem;background:#f8faf9;white-space:nowrap">'
        f'<span style="font-size:.76rem;color:#68756e">{escape(label)}</span>'
        f'<span style="font-size:.88rem;font-weight:600;color:#26352d">{escape(value)}</span></div>'
        for label, value in summary_values
    )
    st.markdown(
        '<div style="display:flex;align-items:center;gap:.45rem;flex-wrap:wrap;margin:.4rem 0 .8rem">'
        f"{summary_items}</div>",
        unsafe_allow_html=True,
    )
    if str(job.get("tool")) == "visual_pipeline":
        pipeline_run = Path(job["run_dir"])
        try:
            saved_graph = json.loads((pipeline_run / "graph_snapshot.json").read_text())
            execution = json.loads((pipeline_run / "graph_execution.json").read_text())
        except (OSError, ValueError):
            saved_graph, execution = {}, {}
        with st.expander("Pipeline steps & saved graph", expanded=True):
            labels = {n["id"]: n["data"]["label"] for n in saved_graph.get("nodes", [])}
            st.dataframe([{"Step": labels.get(n, n), "Status": v.get("status"), "Job": v.get("job_code", ""), "Candidates": v.get("candidate_count"), "Detail": v.get("error", "")} for n, v in execution.get("nodes", {}).items()], hide_index=True, use_container_width=True)
            st.caption("This graph is the saved run snapshot. Editing it as a draft keeps these results intact.")
            st.download_button("Download run graph", json.dumps(saved_graph, indent=2), file_name=f"pipeline-{run_id}.json", mime="application/json", key=f"pipeline-graph-download-{run_id}")
            if saved_graph and project_id and st.button("Edit graph as a draft", key=f"pipeline-graph-edit-{run_id}"):
                prefix = f"visual-pipeline-{project_id}"
                st.session_state[prefix+"-graph"] = saved_graph
                st.session_state[prefix+"-reset"] = secrets.token_hex(16)
                st.session_state[prefix+"-name"] = saved_graph["name"]
                st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = project_id
                st.session_state["pipeline-workspace-view"] = "Node canvas"
                st.session_state["page"] = "Pipelines"
                st.rerun()
    if not project_id:
        st.warning("This existing job has no project assignment. Assign it before resubmitting or analyzing it.")
        assignment_project_id = _render_project_selector(f"job-details-{run_id}", label="Assign project")
        if st.button(
            "Assign job to project",
            type="primary",
            disabled=not assignment_project_id,
            key=f"job-assign-project-{run_id}",
        ):
            try:
                mark_studio_job(
                    {"run_dir": str(job.get("run_dir") or "")},
                    project_id=assignment_project_id,
                    owner_id=str(CURRENT_USER["id"]),
                )
            except (OSError, TypeError, ValueError) as exc:
                st.error(f"Could not assign this job to a project: {exc}")
            else:
                st.success("Job assigned to the selected project.")
                st.rerun()
    resubmit_notice = st.session_state.pop("job-resubmit-notice", None)
    if resubmit_notice:
        st.success(resubmit_notice)
    action_cols = st.columns([1.35, 1, 4.65], vertical_alignment="center")
    run_dir = Path(job["run_dir"]).expanduser().resolve()
    archive_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", title).strip("._-") or "job"
    action_cols[0].download_button(
        "Download job ZIP",
        data=lambda: build_run_zip(run_dir),
        file_name=f"{archive_stem}_{run_id}.zip",
        mime="application/zip",
        use_container_width=True,
        key=f"job-download-zip-{run_id}",
    )
    job_type = str(job.get("job_type") or input_payload.get("job_type") or "")
    is_bindcraft2 = job_type.startswith("bindcraft2") or bool(metadata.get("demo_fixture"))
    analysis_tool = str(metadata.get("tool") or job.get("tool") or input_payload.get("tool") or "").casefold()
    workflow_identity = " ".join(
        str(value or "")
        for value in (
            metadata.get("workflow"),
            input_payload.get("workflow"),
            metadata.get("tool"),
            job.get("tool"),
            input_payload.get("tool"),
        )
    ).casefold().replace("-", "_")
    surface_tool_names = {"pesto", "scannet", "masif_seed", "surf2spot", "surface_analysis"}
    surface_tool = analysis_tool if analysis_tool in surface_tool_names else next(
        (tool for tool in surface_tool_names if tool in workflow_identity),
        None,
    )
    raw_result_candidates = result_payload.get("candidates")
    result_candidates = raw_result_candidates if isinstance(raw_result_candidates, list) else []
    is_af3_jax_cofolding = "af3_jax_cofolding" in workflow_identity or any(
        isinstance(candidate, dict)
        and str(candidate.get("source_tool") or candidate.get("tool") or "").casefold() == "af3_jax_cofolding"
        for candidate in result_candidates
    )
    is_esmfold2_cofolding = "esmfold2_cofolding" in workflow_identity or any(
        isinstance(candidate, dict)
        and str(candidate.get("source_tool") or candidate.get("tool") or "").casefold() == "esmfold2_cofolding"
        for candidate in result_candidates
    )
    input_params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    result_metrics = result.get("metrics") if isinstance(result.get("metrics"), dict) else {}
    is_optimized_prediction = (
        "optimized_prediction" in workflow_identity
        or str(input_params.get("engine") or result_metrics.get("engine") or "") in optimized_prediction.ENGINE_PROFILES
        or any(isinstance(item, dict) and str(item.get("source_tool") or "").startswith("optimized_prediction_") for item in result_candidates)
    )
    is_design_validation = any(isinstance(item, dict) and (item.get("raw_metadata") or {}).get("validation_reference") for item in result_candidates)
    is_rfdiffusion3 = analysis_tool == "rfdiffusion3"
    is_rfantibody = analysis_tool == "rfantibody" or "rfantibody" in workflow_identity
    is_rfpeptides = analysis_tool == "rfpeptides" or "rfpeptides" in workflow_identity
    has_dedicated_prediction_results = is_rfantibody or is_rfpeptides or is_rfdiffusion3 or is_af3_jax_cofolding or is_esmfold2_cofolding or is_optimized_prediction or is_design_validation
    pipeline_execution = {}
    surface_nodes = []
    if str(job.get("tool")) == "visual_pipeline":
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
        pipeline_execution = read_json(run_dir / "graph_execution.json")
        surface_nodes = [
            (node_id, state["output"]["analysis_run"])
            for node_id, state in (pipeline_execution.get("nodes") or {}).items()
            if isinstance(state.get("output"), dict) and isinstance(state["output"].get("analysis_run"), dict)
        ]
    has_dedicated_surface_results = surface_tool is not None or bool(surface_nodes and not result_candidates)
    is_ipsae = analysis_tool in {"ipsae", "analysis.ipsae"} or analysis_tool.endswith(".ipsae")
    if is_bindcraft2:
        resubmit_parameters, resubmit_reason = _bindcraft2_resubmit_parameters(details)
        if not project_id:
            resubmit_reason = "Assign this job to a project before resubmitting it."
        if action_cols[1].button(
            "Resubmit",
            type="secondary",
            disabled=resubmit_reason is not None,
            help=resubmit_reason or "Queue another BindCraft2 run with the saved inputs and settings.",
            use_container_width=True,
            key=f"job-resubmit-{run_id}",
        ):
            request_payload = {
                "schema_version": 1,
                "workflow": "design.bindcraft2_generator" if resubmit_parameters.get("trajectory_only") else "design.bindcraft2",
                "parameters": {key: value for key, value in resubmit_parameters.items() if key != "trajectory_only"},
                "project_id": project_id,
                "priority": 0,
            }
            try:
                validation = _post_local_api("/api/v1/workflows/validate", request_payload)
                if not validation.get("valid"):
                    st.error("The saved BindCraft2 inputs did not pass validation.")
                else:
                    submission_result = _post_local_api("/api/v1/workflows", request_payload)
                    submitted_job = submission_result.get("job", {})
                    mark_studio_job(submission_result, project_id=project_id, owner_id=str(CURRENT_USER["id"]))
                    new_group = str(submitted_job.get("task_group") or "design")
                    new_run_id = str(submitted_job.get("run_id") or "")
                    if new_run_id:
                        st.session_state["selected-job-reference"] = f"{new_group}/{new_run_id}"
                        st.session_state["job-resubmit-notice"] = f"Resubmitted as job {submitted_job.get('job_code') or new_run_id}."
                        st.rerun()
                    st.success(f"Resubmitted job {submitted_job.get('job_code') or 'successfully'}.")
            except (ValueError, TypeError, json.JSONDecodeError) as exc:
                st.error(f"Could not resubmit BindCraft2: {exc}")
        if resubmit_reason:
            action_cols[2].caption(resubmit_reason)

    st.divider()
    results_tab, files_tab, settings_tab, logs_tab, citations_tab = st.tabs(
        ["Results", "Files", "Job settings", "Logs", "Citations"]
    )
    artifact_paths = _job_artifact_paths(details)
    is_generator_tool=str(job.get('tool')) in generator_input.TOOLS.values()

    with results_tab:
        if str(job.get("tool")) == "visual_pipeline":
            if surface_nodes:
                snapshot = read_json(run_dir / "graph_snapshot.json")
                node_labels = {node["id"]: node["data"].get("label", node["id"]) for node in snapshot.get("nodes", [])}
                st.subheader("Protein surface analysis")
                analysis_tabs = st.tabs([node_labels.get(node_id, node_id) for node_id, _analysis in surface_nodes])
                for analysis_tab, (node_id, analysis) in zip(analysis_tabs, surface_nodes):
                    with analysis_tab:
                        try:
                            child_reference = f"{analysis.get('task_group') or 'detection'}/{analysis.get('run_id') or ''}"
                            child_details = _inspect_current_user_job(child_reference)
                            from mn_molecular_design_studio import surface_ui
                            importlib.reload(surface_ui)
                            surface_ui.render_surface_analysis_results(
                                f"{run_id}-{node_id}",
                                Path(child_details["job"]["run_dir"]).resolve(),
                                _job_artifact_paths(child_details),
                                child_details.get("input") if isinstance(child_details.get("input"), dict) else {},
                            )
                        except (AutomationError, OSError, ValueError) as exc:
                            st.warning(f"Could not load the completed surface-analysis step: {exc}")
        if surface_tool is not None:
            from mn_molecular_design_studio import surface_ui
            importlib.reload(surface_ui)
            if surface_tool == "surface_analysis":
                surface_ui.render_surface_analysis_results(run_id, run_dir, artifact_paths, input_payload)
            else:
                surface_ui.render_results(surface_tool, run_id, run_dir, artifact_paths, input_payload)
        elif is_rfantibody:
            from mn_molecular_design_studio import rfantibody_ui
            rfantibody_ui.render_results(run_id, run_dir, result_candidates, input_params, result)
        elif is_rfpeptides:
            from mn_molecular_design_studio import rfpeptides_ui
            rfpeptides_ui.render_results(run_id, run_dir, result_candidates, input_params, result)
        elif is_rfdiffusion3:
            from mn_molecular_design_studio.rfdiffusion3_ui import render_results
            render_results(run_id, run_dir, result_candidates, input_params, result)
        elif str(job.get("tool")) in {"proteina_complexa", "proteina_complexa_generator"}:
            from mn_molecular_design_studio.complexa_ui import render_results
            render_results(run_id,run_dir,result_candidates,input_params,result)
        elif is_generator_tool:
            from mn_molecular_design_studio.generator_ui import render_results
            engine=next(k for k,v in generator_input.TOOLS.items() if v==str(job.get('tool')))
            render_results(run_id,run_dir,result_candidates,{**input_params,'generator_engine':engine},result)
        elif has_dedicated_prediction_results:
            _render_cofolding_prediction_results(run_id, run_dir, result_candidates, input_payload, result)

        structure_artifacts = (
            []
            if has_dedicated_prediction_results or has_dedicated_surface_results or is_generator_tool or str(job.get("tool")) in {"proteina_complexa", "proteina_complexa_generator"}
            else [
                (item, path)
                for item, path in artifact_paths
                if path.suffix.casefold() in {".pdb", ".ent", ".cif", ".mmcif"}
            ]
        )
        structure_artifacts.sort(
            key=lambda pair: (
                0 if "3_ranked" in str(pair[0].get("relative_path", "")).casefold() else 1,
                pair[1].name.casefold(),
            )
        )
        accepted_structures = [pair for pair in structure_artifacts if _structure_design_status(*pair) == "Accepted"]
        not_accepted_structures = [pair for pair in structure_artifacts if _structure_design_status(*pair) == "Not accepted"]
        other_structures = [pair for pair in structure_artifacts if _structure_design_status(*pair) == "Other"]
        if accepted_structures and not_accepted_structures:
            design_scope = st.radio(
                "Design set shown in the table and viewer",
                ["Accepted", "Not accepted", "All"],
                horizontal=True,
                format_func=lambda value: {
                    "Accepted": f"Accepted ({len(accepted_structures)})",
                    "Not accepted": f"Not accepted ({len(not_accepted_structures)})",
                    "All": f"All ({len(structure_artifacts)})",
                }[value],
                key=f"job-design-scope-{run_id}",
            )
        else:
            design_scope = "Accepted" if accepted_structures else "All"
        if design_scope == "Accepted":
            visible_structure_artifacts = accepted_structures
        elif design_scope == "Not accepted":
            visible_structure_artifacts = not_accepted_structures
        else:
            visible_structure_artifacts = structure_artifacts

        csv_artifacts = (
            []
            if has_dedicated_prediction_results or has_dedicated_surface_results
            else [(item, path) for item, path in artifact_paths if path.suffix.casefold() == ".csv"]
        )
        preferred_csv_names = (
            ("ipsae_metrics.csv",) if is_ipsae
            else ("!_ranked.csv", "!_refolded.csv", "!_trajectories.csv")
        )
        csv_artifacts.sort(
            key=lambda pair: (
                preferred_csv_names.index(pair[1].name.casefold())
                if pair[1].name.casefold() in preferred_csv_names
                else len(preferred_csv_names),
                pair[1].name.casefold(),
            )
        )
        csv_loaded = False
        csv_df = None
        if csv_artifacts:
            selected_csv = st.selectbox(
                "Tool output table",
                range(len(csv_artifacts)),
                format_func=lambda index: csv_artifacts[index][1].name,
                key=f"job-output-csv-v2-{run_id}",
            )
            csv_item, csv_path = csv_artifacts[selected_csv]
            try:
                csv_df = pd.read_csv(csv_path, low_memory=False, nrows=50_000)
                native_csv_kind = _bindcraft2_csv_kind(csv_path, csv_df) if is_bindcraft2 and not metadata.get("demo_fixture") else None
                table_df = csv_df.copy()
                if native_csv_kind == "trajectories":
                    table_df = _bindcraft2_trajectory_table(csv_df, structure_artifacts)
                    if design_scope == "Accepted":
                        table_df = table_df[table_df["Design Status"] == "Accepted"].copy()
                    elif design_scope == "Not accepted":
                        table_df = table_df[table_df["Design Status"] == "Not accepted"].copy()
                    st.caption("BindCraft2 trajectory log · one row per attempted design. Terminated records the last stage reached.")
                elif native_csv_kind == "summary":
                    st.caption("BindCraft2 campaign summary · aggregate counts and metrics, not one row per design.")
                    campaign_rows = csv_df[csv_df["scope"].astype(str).str.casefold() == "campaign"]
                    headline_metrics = [
                        (label, metric)
                        for label, metric in (
                            ("Trajectories", "trajectories"),
                            ("Accepted designs", "accepted_designs"),
                            ("Redesign candidates", "redesign_candidates"),
                        )
                        if (campaign_rows["metric"].astype(str) == metric).any()
                    ]
                    if headline_metrics:
                        metric_cols = st.columns(len(headline_metrics))
                        for metric_col, (label, metric) in zip(metric_cols, headline_metrics):
                            metric_row = campaign_rows[campaign_rows["metric"].astype(str) == metric].iloc[0]
                            metric_col.metric(label, str(metric_row.get("mean", "—")))
                elif native_csv_kind == "losses":
                    st.caption("BindCraft2 optimization trace · per-round metrics for one trajectory, rather than a candidate table.")
                    numeric_columns = [
                        column for column in csv_df.columns
                        if column not in {"phase", "round"}
                        and pd.to_numeric(csv_df[column], errors="coerce").notna().any()
                    ]
                    if numeric_columns:
                        metric_name = st.selectbox(
                            "Trace metric",
                            numeric_columns,
                            key=f"job-loss-metric-{run_id}-{selected_csv}",
                        )
                        trace = csv_df.copy()
                        trace["_value"] = pd.to_numeric(trace[metric_name], errors="coerce")
                        trace["_round"] = pd.to_numeric(trace["round"], errors="coerce")
                        phase_order = {name: index for index, name in enumerate(("screen", "refine", "anneal", "harden", "mutate"))}
                        trace["_phase_order"] = trace["phase"].astype(str).str.casefold().map(phase_order).fillna(99)
                        trace = trace.sort_values(["_round", "_phase_order"], kind="stable")
                        trace["Step"] = trace["phase"].astype(str) + " · " + trace["round"].astype(str)
                        st.line_chart(trace.set_index("Step")[["_value"]].rename(columns={"_value": metric_name}))

                if is_ipsae:
                    table_df = _ipsae_display_table(csv_df)
                    st.caption(
                        "ipSAE scores are listed first, followed by the original BindCraft2 ranking, "
                        "prediction confidence, and complex-quality metrics. This view does not rerank designs. "
                        "The full output CSV remains available below."
                    )

                csv_caption, csv_download = st.columns([0.78, 0.22], vertical_alignment="center")
                csv_caption.caption(f"Tool output: `{csv_item.get('relative_path', csv_path.name)}`")
                csv_download.download_button(
                    "Download displayed table",
                    data=table_df.to_csv(index=False).encode("utf-8"),
                    file_name=f"{csv_path.stem}_{design_scope.casefold().replace(' ', '_')}.csv",
                    mime="text/csv",
                    key=f"job-csv-download-{run_id}-{selected_csv}-{design_scope.casefold().replace(' ', '-')}",
                )
                if native_csv_kind == "losses":
                    with st.expander("Raw BindCraft2 trace table", expanded=False):
                        st.dataframe(csv_df, hide_index=True, use_container_width=True)
                elif _sequence_column_name(table_df.columns) is not None:
                    import streamlit.components.v1 as components

                    components.html(
                        _render_tool_output_table_css()
                        + _render_tool_output_table(table_df),
                        height=650,
                        scrolling=False,
                    )
                    if is_ipsae:
                        with st.expander("Full raw ipSAE CSV", expanded=False):
                            st.dataframe(csv_df, hide_index=True, use_container_width=True)
                    else:
                        with st.expander("Interactive raw table", expanded=False):
                            st.dataframe(table_df, hide_index=True, use_container_width=True)
                else:
                    st.dataframe(table_df, hide_index=True, use_container_width=True)
                if design_scope == "Not accepted" and native_csv_kind != "trajectories" and len(table_df):
                    st.caption("Not accepted structures have no ranked metric row; their design, binder sequence, and length are read from the structure file.")
                if len(csv_df) >= 50_000:
                    st.caption("Tool output preview limited to 50,000 rows. Download the original CSV from the Files tab.")
                csv_loaded = True
            except (OSError, ValueError, pd.errors.ParserError, UnicodeError) as exc:
                st.warning(f"Could not read the tool's CSV output: {exc}")

        candidates = [] if has_dedicated_prediction_results or has_dedicated_surface_results else result_payload.get("candidates") or []
        if candidates and not csv_loaded:
            candidate_rows = []
            for candidate in candidates:
                row = {
                    "Candidate": candidate.get("candidate_id") or "—",
                    "Stage": candidate.get("stage") or "—",
                    "Tool": candidate.get("source_tool") or candidate.get("tool") or "—",
                    "Length": candidate.get("binder_length") or "—",
                    "Binder sequence": candidate.get("binder_sequence") or "—",
                }
                metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
                row.update({str(key): value for key, value in metrics.items() if isinstance(value, (str, int, float, bool)) or value is None})
                candidate_rows.append(row)
            st.subheader("Normalized candidates")
            candidate_df = pd.DataFrame(candidate_rows)
            design_id_column = next(
                (column for column in candidate_df.columns if re.sub(r"[^a-z0-9]+", "", str(column).casefold()) == "candidate"),
                None,
            )
            if design_scope == "Accepted":
                if design_id_column is not None and accepted_structures:
                    candidate_ids = candidate_df[design_id_column].astype(str).str.casefold()
                    matching_rows = candidate_ids.isin({path.stem.casefold() for _, path in accepted_structures})
                    if matching_rows.any():
                        candidate_df = candidate_df[matching_rows].copy()
                candidate_df["Design Status"] = "Accepted"
            elif design_scope == "Not accepted":
                candidate_df = _structure_summary_rows(not_accepted_structures, candidate_df.columns)
            else:
                candidate_df["Design Status"] = "Accepted"
                extra_structures = [*not_accepted_structures, *other_structures]
                if extra_structures:
                    candidate_df = pd.concat(
                        [candidate_df, _structure_summary_rows(extra_structures, candidate_df.columns)],
                        ignore_index=True,
                    )
            if _sequence_column_name(candidate_df.columns) is not None:
                import streamlit.components.v1 as components

                components.html(
                    _render_tool_output_table_css()
                    + _render_tool_output_table(candidate_df),
                    height=650,
                    scrolling=False,
                )
                with st.expander("Interactive raw table", expanded=False):
                    st.dataframe(candidate_df, hide_index=True, use_container_width=True)
            else:
                st.dataframe(candidate_df, hide_index=True, use_container_width=True)
            candidate_labels = [str(item.get("candidate_id") or f"Candidate {index + 1}") for index, item in enumerate(candidates)]
            selected_candidate_index = st.selectbox(
                "Candidate details",
                range(len(candidates)),
                format_func=lambda index: candidate_labels[index],
                key=f"job-candidate-{run_id}",
            )
            candidate = candidates[selected_candidate_index]
            sequence = candidate.get("binder_sequence")
            if sequence:
                st.text_area("Binder sequence", str(sequence), height=110, disabled=True, key=f"job-sequence-{run_id}")
            metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
            if metrics:
                st.markdown("**Metrics**")
                st.json(metrics)
            if candidate.get("raw_metadata"):
                with st.expander("Candidate metadata"):
                    st.json(candidate["raw_metadata"])
        elif not candidates and not csv_loaded and not has_dedicated_prediction_results and not has_dedicated_surface_results:
            st.info("No normalized candidate records are available for this job yet.")
            if result:
                st.subheader("Run result")
                st.json(result)

        if visible_structure_artifacts and not has_dedicated_prediction_results and not has_dedicated_surface_results:
            st.subheader("Structure preview")
            structure_index = st.selectbox(
                "Structure file",
                range(len(visible_structure_artifacts)),
                format_func=lambda index: (
                    f"[{_structure_design_status(*visible_structure_artifacts[index])}] "
                    f"{visible_structure_artifacts[index][0].get('relative_path', visible_structure_artifacts[index][1].name)}"
                ),
                key=f"job-structure-{run_id}-{design_scope.casefold().replace(' ', '-')}",
            )
            _, structure_path = visible_structure_artifacts[structure_index]
            try:
                if structure_path.stat().st_size > 50 * 1024 * 1024:
                    st.info("This structure is over 50 MB; download it from the Files tab to inspect it locally.")
                else:
                    structure_text = structure_path.read_text(encoding="utf-8", errors="replace")
                    molstar_custom_component(
                        [StructureVisualization(pdb=structure_text, color="plddt", representation_type="cartoon")],
                        key=f"job-structure-viewer-compact-{run_id}-{structure_index}",
                        height=520,
                        show_controls=False,
                    )
                    st.caption("Structure colored by per-residue pLDDT confidence (stored in the structure B-factor field).")
            except OSError as exc:
                st.warning(f"Could not load this structure preview: {exc}")
        elif not has_dedicated_prediction_results and not has_dedicated_surface_results:
            st.info(f"No {design_scope.casefold()} structure files were found in this job's Results folder.")

        if str(metadata.get("tool") or "").casefold() == "boltzgen":
            _render_boltzgen_campaign_report(run_dir, result_payload)

        if is_bindcraft2:
            _render_bindcraft2_campaign_report(run_dir, input_payload, result)

        result_files = []
        for raw_path in (details.get("campaign_result_path"), str(Path(job["run_dir"]) / "result.json")):
            if raw_path:
                path = Path(raw_path)
                if path.is_file() and path not in result_files:
                    result_files.append(path)
        if result_files:
            st.markdown("**Download result data**")
            download_cols = st.columns(min(len(result_files), 3))
            for index, path in enumerate(result_files):
                try:
                    download_cols[index % len(download_cols)].download_button(
                        f"Download {path.name}",
                        data=path.read_bytes(),
                        file_name=path.name,
                        mime="application/json" if path.suffix.lower() == ".json" else "application/octet-stream",
                        key=f"job-result-download-{run_id}-{index}",
                    )
                except OSError as exc:
                    st.caption(f"Could not read {path.name}: {exc}")

    with files_tab:
        st.subheader("Job files")
        if not artifact_paths:
            st.info("No files are indexed for this job yet.")
        else:
            file_rows = [
                {
                    "File": str(item.get("relative_path", path.name)),
                    "Type": str(item.get("kind") or path.suffix.lstrip(".") or "file"),
                    "Size": _format_bytes(int(item.get("size_bytes") or path.stat().st_size)),
                    "Modified": str(item.get("modified_at") or "—"),
                }
                for item, path in artifact_paths
            ]
            st.dataframe(file_rows, hide_index=True, use_container_width=True)
            file_index = st.selectbox(
                "Choose a file to download",
                range(len(artifact_paths)),
                format_func=lambda index: str(artifact_paths[index][0].get("relative_path", artifact_paths[index][1].name)),
                key=f"job-file-{run_id}",
            )
            item, path = artifact_paths[file_index]
            try:
                st.download_button(
                    f"Download {path.name}",
                    data=path.read_bytes(),
                    file_name=path.name,
                    mime="application/octet-stream",
                    key=f"job-artifact-download-{run_id}",
                )
            except OSError as exc:
                st.error(f"Could not read this file: {exc}")

    with settings_tab:
        st.subheader("Inputs and recorded settings")
        st.json({"input": input_payload, "metadata": metadata})

    with logs_tab:
        logs = [(item, path) for item, path in artifact_paths if path.suffix.casefold() in {".log", ".out", ".err"} or "log" in path.name.casefold()]
        if not logs:
            st.info("No log files are indexed for this job.")
        else:
            log_index = st.selectbox(
                "Log file",
                range(len(logs)),
                format_func=lambda index: str(logs[index][0].get("relative_path", logs[index][1].name)),
                key=f"job-log-{run_id}",
            )
            log_path = logs[log_index][1]
            try:
                content = log_path.read_text(encoding="utf-8", errors="replace")
                st.code(content[-200_000:], language="text")
                if len(content) > 200_000:
                    st.caption("Showing the last 200,000 characters.")
            except OSError as exc:
                st.error(f"Could not read this log: {exc}")

    with citations_tab:
        citation_rows = _citation_values({"input": input_payload, "metadata": metadata, "result": result})
        if citation_rows:
            for label, citation in citation_rows:
                st.markdown(f"**{label}**")
                st.write(citation)
        else:
            st.info("This job does not have citation details recorded in its metadata or result files.")


def _database_groups_file() -> Path:
    return app_home() / "database" / "groups.json"


def _load_database_groups() -> list[dict]:
    path = _database_groups_file()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return []
    groups = [group for group in payload if isinstance(group, dict) and group.get("id") and group.get("name")]
    if CURRENT_USER.get("role") == "admin":
        return groups
    return [group for group in groups if str(group.get("owner_id") or "") == str(CURRENT_USER["id"])]


def _save_database_groups(groups: list[dict]) -> None:
    path = _database_groups_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(".json.tmp")
    temporary_path.write_text(json.dumps(groups, indent=2), encoding="utf-8")
    temporary_path.replace(path)


def _database_source_key(source: dict) -> str:
    return f"{source.get('task_group', '')}/{source.get('run_id', '')}"


def _database_source_name(source: dict) -> str:
    return str(source.get("campaign_name") or source.get("tool") or source.get("job_type") or source.get("run_id") or "Candidate set")


def _guess_database_column(columns: list[str], names: tuple[str, ...]) -> str:
    normalized = {column.lower().replace(" ", "_").replace("-", "_"): column for column in columns}
    for name in names:
        if name in normalized:
            return normalized[name]
    for column in columns:
        key = column.lower().replace(" ", "_").replace("-", "_")
        if any(name in key for name in names):
            return column
    return ""


def _database_column_picker(label: str, columns: list[str], guessed: str, key: str) -> str:
    options = ["", *columns]
    return st.selectbox(
        label,
        options,
        index=options.index(guessed) if guessed in options else 0,
        format_func=lambda value: value or "None",
        key=key,
    )


def _parse_fasta_sequences(text: str) -> list[dict[str, str]]:
    sequences: list[dict[str, str]] = []
    current_name: str | None = None
    current_sequence: list[str] = []
    fallback_lines: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                sequences.append({"candidate_id": current_name, "binder_sequence": "".join(current_sequence)})
            current_name = line[1:].strip().split(maxsplit=1)[0] if line[1:].strip() else f"sequence_{len(sequences) + 1:03d}"
            current_sequence = []
        elif current_name is not None:
            current_sequence.append(line)
        else:
            fallback_lines.append(line)
    if current_name is not None:
        sequences.append({"candidate_id": current_name, "binder_sequence": "".join(current_sequence)})
    elif fallback_lines:
        for index, sequence in enumerate(fallback_lines, start=1):
            sequences.append({"candidate_id": f"sequence_{index:03d}", "binder_sequence": sequence})
    return [row for row in sequences if row["binder_sequence"]]


def _database_table_bytes(rows: list[dict], columns: list[str]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _studio_candidate_sources(task_group: str | None = None) -> list[dict]:
    """Keep this studio's molecule database separate from other apps' jobs."""
    visible = []
    for source in candidate_sources(task_group):
        if not is_studio_job(source.get("run_dir", "")):
            continue
        record = get_job_access_record(str(source.get("task_group") or ""), str(source.get("run_id") or ""))
        if can_access_job(CURRENT_USER, record):
            visible.append(source)
    return visible


def _database_import_complete(run_dir: Path) -> None:
    project_id = validate_project_id(st.session_state.pop("database-import-project-id", None))
    if not project_id:
        raise ValueError("Select a project before importing molecules.")
    require_project_access(CURRENT_USER, project_id)
    st.session_state["database-last-import"] = run_dir.name
    mark_studio_job({"run_dir": str(run_dir)}, project_id=project_id, owner_id=str(CURRENT_USER["id"]))
    group_id = st.session_state.pop("database-import-group-id", None)
    group_name = st.session_state.pop("database-import-group-name", None)
    imported_refs: list[dict] = []
    try:
        source = next(
            (item for item in _studio_candidate_sources("candidate-import") if Path(item.get("run_dir", "")).resolve() == run_dir.resolve()),
            None,
        )
        if source:
            imported_refs = [
                {
                    "task_group": source.get("task_group"),
                    "run_id": source.get("run_id"),
                    "candidate_id": candidate.get("candidate_id"),
                }
                for candidate in load_source_candidates(source)
            ]
    except Exception:
        imported_refs = []
    if group_id and group_name and imported_refs:
        groups = _load_database_groups()
        group = next((item for item in groups if item.get("id") == group_id), None)
        if group is None:
            group = {
                "id": group_id,
                "name": group_name,
                "owner_id": str(CURRENT_USER["id"]),
                "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "members": [],
            }
            groups.append(group)
        existing_refs = {
            (item.get("task_group"), item.get("run_id"), item.get("candidate_id"))
            for item in group.get("members", [])
        }
        group["members"] = [
            *group.get("members", []),
            *[
                ref for ref in imported_refs
                if (ref.get("task_group"), ref.get("run_id"), ref.get("candidate_id")) not in existing_refs
            ],
        ]
        _save_database_groups(groups)
        st.session_state["database-preselect-group"] = f"custom:{group_id}"
    else:
        st.session_state["database-preselect-group"] = f"source:{run_dir.parent.name}/{run_dir.name}"
    st.session_state["database-import-dialog-requested"] = False
    st.session_state["database-import-step"] = 0
    st.rerun()


def _import_uploaded_structures(
    files: list[tuple[str, bytes]],
    *,
    import_name: str,
    binder_chains_text: str,
    target_chains_text: str,
    structures_are_complexes: bool,
) -> Path:
    valid_suffixes = (".pdb", ".ent", ".cif", ".mmcif", ".pdb.gz", ".cif.gz", ".mmcif.gz")
    selected_files = [(name, content) for name, content in files if name.lower().endswith(valid_suffixes)]
    if not selected_files:
        raise ValueError("No PDB or mmCIF structure files were found.")
    if len(selected_files) > 500:
        raise ValueError("A single import can contain at most 500 structures.")
    app_home().mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="molecular-import-", dir=app_home()) as temporary_dir:
        source_dir = Path(temporary_dir)
        rows = []
        for index, (name, content) in enumerate(selected_files, start=1):
            safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", Path(name).name).strip("._") or f"structure_{index}.pdb"
            if any((source_dir / row["structure_path"]).name == safe_name for row in rows):
                safe_name = f"{index:04d}_{safe_name}"
            structure_path = source_dir / safe_name
            structure_path.write_bytes(content)
            rows.append(
                {
                    "candidate_id": Path(safe_name).name.split(".")[0],
                    "structure_path": safe_name,
                    "binder_chains": binder_chains_text.strip(),
                    "target_chains": target_chains_text.strip(),
                }
            )
        table_bytes = _database_table_bytes(rows, ["candidate_id", "structure_path", "binder_chains", "target_chains"])
        kwargs = {
            "table_bytes": table_bytes,
            "table_filename": "structure_import.csv",
            "base_dir": source_dir,
            "id_column": "candidate_id",
            "binder_chains_column": "binder_chains" if binder_chains_text.strip() else "",
            "target_chains_column": "target_chains" if target_chains_text.strip() else "",
            "default_binder_chains": [part.strip() for part in re.split(r"[,; ]+", binder_chains_text) if part.strip()] or ["A"],
            "default_target_chains": [part.strip() for part in re.split(r"[,; ]+", target_chains_text) if part.strip()],
            "import_name": import_name,
            "source_tool": "structure_import",
        }
        if structures_are_complexes:
            kwargs["complex_path_column"] = "structure_path"
        else:
            kwargs["binder_path_column"] = "structure_path"
            kwargs["binder_chains_column"] = ""
            kwargs["target_chains_column"] = ""
        return candidate_import.run_generic_table_import(**kwargs)


@st.dialog("Import molecules", width="large")
def _render_database_import_dialog() -> None:
    if int(st.session_state.get("database-import-step", 0)) == 0:
        groups = _load_database_groups()
        choices = ["__new__", *[str(group["id"]) for group in groups]]
        group_by_id = {str(group["id"]): group for group in groups}
        if st.session_state.get("database-import-group-choice") not in choices:
            st.session_state["database-import-group-choice"] = "__new__"
        st.caption("Group  ›  What are you importing?")
        group_choice = st.selectbox(
            "Group",
            choices,
            format_func=lambda value: "New group" if value == "__new__" else group_by_id[value]["name"],
            key="database-import-group-choice",
        )
        if group_choice == "__new__":
            group_name = st.text_input("New group name", placeholder="e.g. EGFR binder candidates", key="database-import-new-group-name").strip()
        else:
            group_name = str(group_by_id[group_choice]["name"])
            st.session_state["database-import-new-group-name"] = group_name
        st.session_state["database-import-project-id"] = _render_project_selector(
            "database-import",
            allow_create=True,
        )
        back_col, next_col = st.columns([1, 1])
        if back_col.button("Cancel", key="database-import-cancel-from-group"):
            st.session_state["database-import-dialog-requested"] = False
            st.session_state["database-import-step"] = 0
            st.session_state.pop("database-import-project-id", None)
            st.rerun()
        if next_col.button("Next", type="primary", disabled=not group_name or not st.session_state.get("database-import-project-id"), key="database-import-next"):
            if group_choice == "__new__":
                group_id = hashlib.sha1(f"{datetime.now(timezone.utc).isoformat()}-{group_name}".encode("utf-8")).hexdigest()[:16]
            else:
                group_id = group_choice
            st.session_state["database-import-group-id"] = group_id
            st.session_state["database-import-group-name"] = group_name
            st.session_state["database-import-step"] = 1
            st.rerun()
        return

    import_group_name = str(st.session_state.get("database-import-group-name") or "Molecule import")
    st.caption(f"Group  ›  **{import_group_name}**  ›  What are you importing?")
    nav_cols = st.columns([1, 1, 5])
    if nav_cols[0].button("Cancel", key="database-import-cancel-from-method"):
        st.session_state["database-import-dialog-requested"] = False
        st.session_state["database-import-step"] = 0
        st.session_state.pop("database-import-project-id", None)
        st.rerun()
    if nav_cols[1].button("← Back", key="database-import-back-to-group"):
        st.session_state["database-import-step"] = 0
        st.rerun()
    method = st.radio(
        "What are you importing?",
        ["Table", "Structure file", "ZIP of structures", "FASTA", "Fetch by identifier", "Paste sequences"],
        horizontal=True,
        key="database-import-method",
    )

    if method == "Table":
        st.caption("Import candidate sequences or structure paths from CSV, TSV, or Excel.")
        uploaded = st.file_uploader("Candidate table", type=["csv", "tsv", "txt", "xlsx", "xls"], key="database-import-table")
        if not uploaded:
            return
        data = uploaded.getvalue()
        try:
            if uploaded.name.lower().endswith((".xlsx", ".xls")):
                preview = pd.read_excel(io.BytesIO(data))
            else:
                separator = "\t" if uploaded.name.lower().endswith((".tsv", ".txt")) else ","
                preview = pd.read_csv(io.BytesIO(data), sep=separator)
        except Exception as exc:
            st.error(f"Could not read this table: {exc}")
            return
        columns = [str(column) for column in preview.columns]
        st.dataframe(preview.head(10), hide_index=True, use_container_width=True)
        map_cols = st.columns(3)
        id_column = map_cols[0].selectbox("Candidate ID", ["", *columns], index=(["", *columns].index(_guess_database_column(columns, ("candidate_id", "binder_id", "design_id", "name", "id"))) if _guess_database_column(columns, ("candidate_id", "binder_id", "design_id", "name", "id")) in columns else 0), format_func=lambda value: value or "None", key="database-table-id-column")
        sequence_guess = _guess_database_column(columns, ("binder_sequence", "sequence", "seq", "binder_seq", "aa_sequence"))
        sequence_column = map_cols[1].selectbox("Binder sequence", ["", *columns], index=(["", *columns].index(sequence_guess) if sequence_guess in columns else 0), format_func=lambda value: value or "None", key="database-table-sequence-column")
        complex_guess = _guess_database_column(columns, ("complex_pdb", "complex_path", "structure_path", "pdb", "structure"))
        admin_path_import = CURRENT_USER.get("role") == "admin"
        structure_path_options = ["", *columns] if admin_path_import else [""]
        complex_column = map_cols[2].selectbox("Complex structure path", structure_path_options, index=(structure_path_options.index(complex_guess) if complex_guess in structure_path_options else 0), format_func=lambda value: value or "None", key="database-table-complex-column", disabled=not admin_path_import)
        map_cols_2 = st.columns(3)
        binder_guess = _guess_database_column(columns, ("binder_pdb", "binder_path", "binder_structure"))
        binder_column = map_cols_2[0].selectbox("Binder structure path", structure_path_options, index=(structure_path_options.index(binder_guess) if binder_guess in structure_path_options else 0), format_func=lambda value: value or "None", key="database-table-binder-column", disabled=not admin_path_import)
        binder_chain_guess = _guess_database_column(columns, ("binder_chains", "binder_chain", "designed_chain"))
        binder_chains_column = map_cols_2[1].selectbox("Binder chain(s) column", ["", *columns], index=(["", *columns].index(binder_chain_guess) if binder_chain_guess in columns else 0), format_func=lambda value: value or "None", key="database-table-binder-chains-column")
        target_chain_guess = _guess_database_column(columns, ("target_chains", "target_chain"))
        target_chains_column = map_cols_2[2].selectbox("Target chain(s) column", ["", *columns], index=(["", *columns].index(target_chain_guess) if target_chain_guess in columns else 0), format_func=lambda value: value or "None", key="database-table-target-chains-column")
        options = st.columns(3)
        base_dir = options[0].text_input("Path base folder", placeholder="/path/to/relative/structure/files", key="database-table-base-dir", disabled=not admin_path_import)
        default_binder_chains = options[1].text_input("Default binder chain(s)", value="A", key="database-table-default-binder-chains")
        default_target_chains = options[2].text_input("Default target chain(s)", value="", key="database-table-default-target-chains")
        target_structure = st.text_input("Optional target structure path", placeholder="/path/to/target.pdb", key="database-table-target-structure", disabled=not admin_path_import)
        import_name = import_group_name or Path(uploaded.name).stem
        import_ready = bool(sequence_column or (admin_path_import and (complex_column or binder_column)))
        if not import_ready:
            if admin_path_import:
                st.info("Select a sequence column or a structure path column to import candidates.")
            else:
                st.info("Select a sequence column to import candidates. Host filesystem structure paths are available to the admin only.")
        if st.button("Import molecules", type="primary", disabled=not import_ready or not st.session_state.get("database-import-project-id"), key="database-table-import-button"):
            try:
                if not admin_path_import and (complex_column or binder_column or base_dir.strip() or target_structure.strip()):
                    raise PermissionError("Host filesystem structure paths can only be imported by the admin. Upload structures or import sequences instead.")
                run_dir = candidate_import.run_generic_table_import(
                    table_bytes=data,
                    table_filename=uploaded.name,
                    base_dir=Path(base_dir).expanduser() if base_dir.strip() else None,
                    id_column=id_column,
                    sequence_column=sequence_column,
                    complex_path_column=complex_column,
                    binder_path_column=binder_column,
                    binder_chains_column=binder_chains_column,
                    target_chains_column=target_chains_column,
                    target_structure=Path(target_structure).expanduser() if target_structure.strip() else None,
                    default_binder_chains=[part.strip() for part in re.split(r"[,; ]+", default_binder_chains) if part.strip()] or ["A"],
                    default_target_chains=[part.strip() for part in re.split(r"[,; ]+", default_target_chains) if part.strip()],
                    import_name=import_name.strip() or Path(uploaded.name).stem,
                    source_tool="table_import",
                )
                _database_import_complete(run_dir)
            except Exception as exc:
                st.error(f"Import failed: {exc}")
        return

    if method in {"FASTA", "Paste sequences"}:
        st.caption("One candidate per FASTA record or non-empty sequence line.")
        if method == "FASTA":
            fasta_file = st.file_uploader("FASTA file", type=["fasta", "fa", "faa", "txt"], key="database-import-fasta")
            fasta_text = fasta_file.getvalue().decode("utf-8", errors="replace") if fasta_file else ""
        else:
            fasta_text = st.text_area("Protein sequences", placeholder=">binder_1\nM...\n>binder_2\nM...", height=180, key="database-import-sequences")
        candidates = _parse_fasta_sequences(fasta_text)
        if candidates:
            st.dataframe(candidates[:10], hide_index=True, use_container_width=True)
        import_name = import_group_name or "Sequence import"
        if st.button("Import molecules", type="primary", disabled=not candidates or not st.session_state.get("database-import-project-id"), key="database-sequence-import-button"):
            try:
                table_bytes = _database_table_bytes(candidates, ["candidate_id", "binder_sequence"])
                run_dir = candidate_import.run_generic_table_import(
                    table_bytes=table_bytes,
                    table_filename="sequences.csv",
                    id_column="candidate_id",
                    sequence_column="binder_sequence",
                    import_name=import_name.strip() or "Sequence import",
                    source_tool="sequence_import",
                )
                _database_import_complete(run_dir)
            except Exception as exc:
                st.error(f"Import failed: {exc}")
        return

    if method in {"Structure file", "ZIP of structures"}:
        complex_mode = st.checkbox("These are binder–target complexes", value=False, key="database-structure-complex-mode")
        structure_cols = st.columns(2)
        binder_chains = structure_cols[0].text_input("Binder chain(s)", value="A", key="database-structure-binder-chains")
        target_chains = structure_cols[1].text_input("Target chain(s)", value="", key="database-structure-target-chains")
        if method == "Structure file":
            uploads = st.file_uploader(
                "Structure files",
                type=["pdb", "ent", "cif", "mmcif", "gz"],
                accept_multiple_files=True,
                key="database-import-structures",
            )
            files = [(upload.name, upload.getvalue()) for upload in uploads]
        else:
            archive = st.file_uploader("ZIP of PDB/mmCIF structures", type=["zip"], key="database-import-zip")
            files = []
            if archive:
                try:
                    with zipfile.ZipFile(io.BytesIO(archive.getvalue())) as bundle:
                        archive_root = Path("molecules").resolve()
                        members = [item for item in bundle.infolist() if not item.is_dir()]
                        if len(members) > 500:
                            raise ValueError("A single import can contain at most 500 structures.")
                        for item in members:
                            member_path = (archive_root / item.filename).resolve()
                            if archive_root not in member_path.parents:
                                raise ValueError("The ZIP contains an unsafe path.")
                            if item.filename.lower().endswith((".pdb", ".ent", ".cif", ".mmcif", ".pdb.gz", ".cif.gz", ".mmcif.gz")):
                                files.append((Path(item.filename).name, bundle.read(item)))
                except (OSError, zipfile.BadZipFile, RuntimeError, ValueError) as exc:
                    st.error(f"Could not read this ZIP: {exc}")
                    return
        import_name = import_group_name or "Structure import"
        if st.button("Import molecules", type="primary", disabled=not files or not st.session_state.get("database-import-project-id"), key="database-structure-import-button"):
            try:
                run_dir = _import_uploaded_structures(
                    files,
                    import_name=import_name.strip() or "Structure import",
                    binder_chains_text=binder_chains,
                    target_chains_text=target_chains,
                    structures_are_complexes=complex_mode,
                )
                _database_import_complete(run_dir)
            except Exception as exc:
                st.error(f"Import failed: {exc}")
        return

    if method == "Fetch by identifier":
        st.caption("Fetch a PDB structure from RCSB using its four-character accession.")
        code_col, fetch_col = st.columns([3, 1], vertical_alignment="bottom")
        pdb_code = code_col.text_input("PDB accession", placeholder="e.g. 7JZU", key="database-import-pdb-code").strip().upper()
        if fetch_col.button("Fetch PDB", key="database-import-fetch-pdb"):
            if not re.fullmatch(r"[A-Z0-9]{4}", pdb_code):
                st.error("Enter a four-character PDB accession.")
            else:
                try:
                    request = Request(
                        f"https://files.rcsb.org/download/{pdb_code}.pdb",
                        headers={"User-Agent": "MN-Molecular-Design-Studio/0.1"},
                    )
                    with urlopen(request, timeout=30) as response:
                        content = response.read(50_000_001)
                    if len(content) > 50_000_000:
                        raise ValueError("The PDB file is larger than the 50 MB import limit.")
                    st.session_state["database-fetched-pdb"] = {"name": f"{pdb_code}.pdb", "content": content, "code": pdb_code}
                    st.rerun()
                except (HTTPError, URLError, OSError, ValueError) as exc:
                    st.error(f"Could not fetch {pdb_code}: {exc}")
        fetched = st.session_state.get("database-fetched-pdb")
        if fetched and fetched.get("code") == pdb_code:
            st.success(f"Fetched {pdb_code}.pdb")
            complex_mode = st.checkbox("This is a binder–target complex", value=True, key="database-fetched-complex-mode")
            chain_cols = st.columns(2)
            binder_chains = chain_cols[0].text_input("Binder chain(s)", value="B", key="database-fetched-binder-chains")
            target_chains = chain_cols[1].text_input("Target chain(s)", value="", key="database-fetched-target-chains")
            import_name = import_group_name or pdb_code
            if st.button("Import molecule", type="primary", disabled=not st.session_state.get("database-import-project-id"), key="database-fetched-import-button"):
                try:
                    run_dir = _import_uploaded_structures(
                        [(fetched["name"], fetched["content"])],
                        import_name=import_name.strip() or pdb_code,
                        binder_chains_text=binder_chains,
                        target_chains_text=target_chains,
                        structures_are_complexes=complex_mode,
                    )
                    _database_import_complete(run_dir)
                except Exception as exc:
                    st.error(f"Import failed: {exc}")


@st.dialog("Create molecule group", width="small")
def _render_database_group_dialog(filtered_refs: list[dict], selected_refs: list[dict]) -> None:
    st.caption(f"{len(filtered_refs)} molecule{'s' if len(filtered_refs) != 1 else ''} match the current view.")
    use_selected = False
    if selected_refs:
        use_selected = st.checkbox("Use selected molecules only", value=True, key="database-group-selected-only")
    group_refs = selected_refs if use_selected else filtered_refs
    group_name = st.text_input("Group name", key="database-new-group-name")
    if st.button("Cancel", key="database-create-group-cancel"):
        st.session_state["database-create-group-requested"] = False
        st.rerun()
    if st.button("Create group", type="primary", disabled=not group_refs or not group_name.strip(), key="database-create-group-button"):
        groups = _load_database_groups()
        group = {
            "id": hashlib.sha1(f"{datetime.now(timezone.utc).isoformat()}-{group_name}".encode("utf-8")).hexdigest()[:16],
            "name": group_name.strip(),
            "owner_id": str(CURRENT_USER["id"]),
            "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "members": group_refs,
        }
        groups.append(group)
        _save_database_groups(groups)
        st.session_state["database-preselect-group"] = f"custom:{group['id']}"
        st.session_state["database-create-group-requested"] = False
        st.session_state.pop("database-create-group-refs", None)
        st.rerun()


def _render_database() -> None:
    st.title("Database")
    if st.session_state.pop("database-last-import", None):
        st.success("Molecules imported. The new candidate set is available in Groups.")
    try:
        sources = _studio_candidate_sources()
    except Exception as exc:
        st.error(f"Could not load candidate sets: {exc}")
        sources = []
    groups = _load_database_groups()
    source_by_key = {_database_source_key(source): source for source in sources}
    groups = [
        {
            **group,
            "members": [
                member for member in group.get("members", [])
                if isinstance(member, dict)
                and f"{member.get('task_group', '')}/{member.get('run_id', '')}" in source_by_key
            ],
        }
        for group in groups
    ]
    grouped_source_keys = {
        f"{member.get('task_group', '')}/{member.get('run_id', '')}"
        for group in groups
        for member in group.get("members", [])
        if isinstance(member, dict)
    }
    source_entries = [
        {
            "id": f"source:{_database_source_key(source)}",
            "name": _database_source_name(source),
            "count": source.get("candidate_count", 0),
            "method": source.get("tool") or source.get("job_type") or "Workflow",
            "source": source,
            "members": None,
        }
        for source in sources
        if _database_source_key(source) not in grouped_source_keys
    ]
    custom_entries = [
        {"id": f"custom:{group['id']}", "name": group["name"], "count": len(group.get("members", [])), "method": "Custom group", "source": None, "members": group.get("members", [])}
        for group in groups
    ]
    entries = [*custom_entries, *source_entries]

    side_col, main_col = st.columns([0.23, 0.77], gap="medium")
    with side_col:
        title_col, actions_col = st.columns([0.7, 0.3], vertical_alignment="center")
        title_col.markdown("#### Groups")
        if actions_col.button("+ New", key="database-new-group", use_container_width=True, help="Create a group from the current filtered molecules"):
            st.session_state["database-create-group-requested"] = True
            st.rerun()
        if st.button("Refresh groups", icon=":material/refresh:", key="database-refresh-groups", use_container_width=True):
            st.rerun()
        group_query = st.text_input("Search groups or jobs", placeholder="Search groups or jobs…", label_visibility="collapsed", key="database-group-search").strip().casefold()
        group_scope = st.radio("Group scope", ["My groups", "Organization"], horizontal=True, label_visibility="collapsed", key="database-group-scope")
        method_options = ["All methods", *sorted({entry["method"] for entry in entries})]
        if st.session_state.get("database-method-filter") not in method_options:
            st.session_state["database-method-filter"] = "All methods"
        method_filter = st.selectbox("Import method", method_options, key="database-method-filter")
        filtered_entries = [
            entry for entry in entries
            if (not group_query or group_query in entry["name"].casefold() or group_query in entry["method"].casefold())
            and (method_filter == "All methods" or entry["method"] == method_filter)
        ]
        if group_scope == "Organization":
            organization = load_organization(STUDIO_WORKSPACE)
            st.info(f"Shared molecule groups for {organization.get('name', 'the organization')} are not enabled yet.")
            filtered_entries = []
        if filtered_entries:
            preselected_group = st.session_state.pop("database-preselect-group", None)
            available_ids = [entry["id"] for entry in filtered_entries]
            if preselected_group in available_ids:
                st.session_state["database-selected-group"] = preselected_group
            elif st.session_state.get("database-selected-group") not in available_ids:
                st.session_state["database-selected-group"] = available_ids[0]
            selected_group_id = st.radio(
                "Candidate sets",
                available_ids,
                format_func=lambda entry_id: next(
                    f"{entry['name']}  ·  {entry['count']}" for entry in filtered_entries if entry["id"] == entry_id
                ),
                key="database-selected-group",
                label_visibility="collapsed",
            )
        else:
            selected_group_id = ""
            if group_scope == "My groups":
                st.caption("No molecule groups found.")

    selected_entry = next((entry for entry in entries if entry["id"] == selected_group_id), None)
    all_candidates: list[dict] = []
    if selected_entry:
        if selected_entry["source"]:
            try:
                source = selected_entry["source"]
                task_group, _, run_id = _database_source_key(source).partition("/")
                all_candidates = [
                    {
                        **candidate,
                        "_database_ref": {
                            "task_group": task_group,
                            "run_id": run_id,
                            "candidate_id": candidate.get("candidate_id"),
                        },
                    }
                    for candidate in load_source_candidates(source)
                ]
            except Exception as exc:
                st.error(f"Could not open this group: {exc}")
        else:
            for member in selected_entry["members"]:
                source = source_by_key.get(f"{member.get('task_group', '')}/{member.get('run_id', '')}")
                if not source:
                    continue
                try:
                    matched = next(
                        (candidate for candidate in load_source_candidates(source) if candidate.get("candidate_id") == member.get("candidate_id")),
                        None,
                    )
                except Exception:
                    matched = None
                if matched:
                    all_candidates.append({**matched, "_database_ref": member})

    with main_col:
        heading_col, menu_col, export_col, columns_col = st.columns([0.44, 0.18, 0.18, 0.2], vertical_alignment="center")
        heading_col.subheader("Molecules")
        with menu_col.popover("File", use_container_width=True):
            if st.button("Upload molecules…", icon=":material/upload_file:", key="database-open-import", use_container_width=True):
                st.session_state["database-import-dialog-requested"] = True
                st.session_state["database-import-step"] = 0
                st.rerun()
            group_button = st.button(
                "Group from filter",
                icon=":material/create_new_folder:",
                disabled=not all_candidates,
                key="database-group-from-filter",
                use_container_width=True,
                help="Save the filtered molecules as a group.",
            )
            if group_button:
                st.session_state["database-create-group-requested"] = True
        with export_col:
            export_placeholder = st.empty()
        all_columns = ["Candidate", "Stage", "Source", "Binder sequence", "Length", "Binder chains", "Target chains"]
        with columns_col.popover("Columns", icon=":material/view_column:", use_container_width=True):
            visible_columns = st.multiselect(
                "Choose columns",
                all_columns,
                default=all_columns[:5],
                key="database-visible-columns",
                label_visibility="collapsed",
            )
        if not selected_entry:
            st.info("Choose a candidate group or import molecules to get started.")
            filtered_candidates: list[dict] = []
            selected_candidates: list[dict] = []
        else:
            search_col, sort_col, order_col, filter_col = st.columns([0.48, 0.2, 0.13, 0.19], vertical_alignment="center")
            candidate_query = search_col.text_input("Search within this group", placeholder="Search within this group…", key="database-candidate-search").strip().casefold()
            sort_by = sort_col.selectbox("Sort", ["Added", "Candidate", "Binder length"], key="database-sort")
            descending = order_col.toggle("Desc", value=True, key="database-sort-desc")
            with filter_col.popover("Filter", icon=":material/filter_alt:", use_container_width=True):
                stages = sorted({str(candidate.get("stage") or "Unknown") for candidate in all_candidates})
                if not set(st.session_state.get("database-stage-filter", [])).issubset(stages):
                    st.session_state["database-stage-filter"] = []
                selected_stages = st.multiselect("Stage", stages, key="database-stage-filter")
                tools = sorted({str(candidate.get("source_tool") or "Unknown") for candidate in all_candidates})
                if not set(st.session_state.get("database-source-filter", [])).issubset(tools):
                    st.session_state["database-source-filter"] = []
                selected_tools = st.multiselect("Source", tools, key="database-source-filter")
            filtered_candidates = []
            for candidate in all_candidates:
                sequence = str(candidate.get("binder_sequence") or "")
                searchable = " ".join((str(candidate.get("candidate_id") or ""), str(candidate.get("source_tool") or ""), sequence)).casefold()
                if candidate_query and candidate_query not in searchable:
                    continue
                if selected_stages and str(candidate.get("stage") or "Unknown") not in selected_stages:
                    continue
                if selected_tools and str(candidate.get("source_tool") or "Unknown") not in selected_tools:
                    continue
                filtered_candidates.append(candidate)
            if sort_by == "Candidate":
                filtered_candidates.sort(key=lambda candidate: str(candidate.get("candidate_id") or "").casefold(), reverse=descending)
            elif sort_by == "Binder length":
                filtered_candidates.sort(key=lambda candidate: len(str(candidate.get("binder_sequence") or "")), reverse=descending)
            else:
                filtered_candidates.sort(key=lambda candidate: int((candidate.get("metrics") or {}).get("import_index", 0) or 0), reverse=descending)
            table_rows = [
                {
                    "Candidate": candidate.get("candidate_id", ""),
                    "Stage": str(candidate.get("stage") or "").replace("_", " ").title(),
                    "Source": candidate.get("source_tool") or "—",
                    "Binder sequence": candidate.get("binder_sequence") or "",
                    "Length": len(str(candidate.get("binder_sequence") or "")),
                    "Binder chains": ", ".join(candidate.get("binder_chains") or []),
                    "Target chains": ", ".join(candidate.get("target_chains") or []),
                }
                for candidate in filtered_candidates
            ]
            selected_cols = visible_columns or all_columns[:1]
            export_rows = [{key: row[key] for key in all_columns} for row in table_rows]
            export_placeholder.download_button(
                "Export",
                data=_database_table_bytes(export_rows, all_columns),
                file_name=f"{selected_entry['name'][:60].replace(' ', '-')}-molecules.csv",
                mime="text/csv",
                icon=":material/download:",
                disabled=not table_rows,
                use_container_width=True,
            )
            if selected_cols:
                table_key = "database-molecules-" + hashlib.sha1(
                    "|".join(str(row.get("candidate_id") or "") for row in filtered_candidates).encode("utf-8")
                ).hexdigest()[:12]
                selection = st.dataframe(
                    table_rows,
                    hide_index=True,
                    use_container_width=True,
                    column_order=selected_cols,
                    on_select="rerun",
                    selection_mode="multi-row",
                    key=table_key,
                )
                selected_indices = list(selection.selection.rows) if selection and selection.selection else []
                selected_candidates = [filtered_candidates[index] for index in selected_indices if 0 <= index < len(filtered_candidates)]
            else:
                st.warning("Choose at least one column to display molecules.")
                selected_candidates = []
            st.caption(f"{len(filtered_candidates)} molecules in this group")
            if len(selected_candidates) == 1:
                candidate = selected_candidates[0]
                with st.expander(f"{candidate.get('candidate_id', 'Molecule')} details", expanded=True):
                    st.code(candidate.get("binder_sequence") or "No binder sequence recorded", language="text")
                    st.json({key: candidate.get(key) for key in ("stage", "source_tool", "target_chains", "binder_chains", "hotspots", "metrics", "raw_metadata")})
                    for path_key in ("binder_pdb", "complex_pdb", "target_pdb"):
                        path_text = candidate.get(path_key)
                        path = Path(str(path_text)) if path_text else None
                        if path and path.is_file():
                            st.download_button(
                                f"Download {path_key.replace('_', ' ')}",
                                data=path.read_bytes(),
                                file_name=path.name,
                                mime="chemical/x-pdb" if path.suffix.lower() == ".pdb" else "application/octet-stream",
                                key=f"database-download-{path_key}-{candidate.get('candidate_id')}",
                            )

    if st.session_state.get("database-import-dialog-requested", False):
        _render_database_import_dialog()
    elif st.session_state.get("database-create-group-requested", False):
        refs = [candidate["_database_ref"] for candidate in filtered_candidates if candidate.get("_database_ref")]
        selection_refs = [candidate["_database_ref"] for candidate in selected_candidates if candidate.get("_database_ref")]
        _render_database_group_dialog(refs, selection_refs)


def _render_account() -> None:
    st.title("Account")
    st.write(f"Signed in as **{CURRENT_USER['display_name']}** (`{CURRENT_USER['username']}`)")
    if is_development_mode():
        st.info("Development mode uses a temporary local admin identity. Switch to production mode to test real accounts and password changes.")
        return
    st.caption(f"Role: {str(CURRENT_USER['role']).title()} · Sessions expire after 12 hours.")
    with st.form("account-change-password-form"):
        current_password = st.text_input("Current password", type="password")
        new_password = st.text_input("New password", type="password")
        confirm_password = st.text_input("Confirm new password", type="password")
        submitted = st.form_submit_button("Change password", type="primary")
    if submitted:
        if new_password != confirm_password:
            st.error("The new passwords do not match.")
        else:
            try:
                change_password(str(CURRENT_USER["id"]), new_password, current_password=current_password)
            except (OSError, sqlite3.Error, TypeError, ValueError) as exc:
                st.error(str(exc))
            else:
                revoke_session(st.session_state.get("studio-session-token"))
                st.session_state.clear()
                st.rerun()
    st.caption(f"Passwords require at least {password_min_length()} characters.")


def _render_user_access() -> None:
    st.title("Users & access")
    if cluster_config.is_worker():
        st.info("Accounts and project access are managed on the coordinator so the same user has the same identity on every workstation.")
        return
    if is_development_mode():
        st.info("User management is disabled in development mode. Switch to production mode and create the real admin account to test account grants.")
        return
    st.caption("Only admins can create accounts, reset passwords, disable accounts, and grant project access.")

    with st.form("admin-create-user-form", clear_on_submit=True):
        st.subheader("Add a user")
        username_col, name_col = st.columns(2)
        username = username_col.text_input("Username or email", key="admin-new-username")
        display_name = name_col.text_input("Display name", key="admin-new-display-name")
        password_col, confirm_col = st.columns(2)
        temporary_password = password_col.text_input("Temporary password", type="password", key="admin-new-password")
        confirm_password = confirm_col.text_input("Confirm temporary password", type="password", key="admin-new-password-confirm")
        create_submitted = st.form_submit_button("Create account", type="primary")
    if create_submitted:
        if temporary_password != confirm_password:
            st.error("The temporary passwords do not match.")
        else:
            try:
                created = create_user(username, display_name, temporary_password, actor=CURRENT_USER)
            except (OSError, sqlite3.Error, TypeError, ValueError) as exc:
                st.error(f"Could not create account: {exc}")
            else:
                st.success(f"Account created for {created['username']}. Share the temporary password with them securely; they must change it at first sign-in.")

    users = list_users()
    projects = _studio_projects()
    project_options = {str(project["id"]): str(project["name"]) for project in projects}
    st.divider()
    st.subheader("Accounts")
    for user in users:
        with st.container(border=True):
            title_col, status_col = st.columns([0.7, 0.3], vertical_alignment="center")
            title_col.markdown(f"**{escape(str(user['display_name']))}** · `{escape(str(user['username']))}`")
            if user["role"] == "admin":
                status_col.caption("Admin · all data")
            else:
                status_col.caption("Disabled" if user["disabled"] else "User")

            if user["role"] == "admin":
                continue

            access_key = f"admin-project-access-{user['id']}"
            current_access = [project_id for project_id in project_member_ids(str(user["id"])) if project_id in project_options]
            selected_access = st.multiselect(
                "Projects this user can access",
                options=list(project_options),
                default=current_access,
                format_func=lambda value: project_options.get(value, "Unknown project"),
                key=access_key,
            )
            access_col, reset_col, account_col = st.columns([0.45, 0.3, 0.25], vertical_alignment="center")
            if access_col.button("Save project access", key=f"admin-save-access-{user['id']}"):
                try:
                    for project_id in selected_access:
                        grant_project_access(project_id, str(user["id"]), actor=CURRENT_USER)
                    for project_id in set(current_access) - set(selected_access):
                        revoke_project_access(project_id, str(user["id"]), actor=CURRENT_USER)
                except (OSError, sqlite3.Error, ValueError) as exc:
                    st.error(f"Could not update project access: {exc}")
                else:
                    st.success(f"Project access updated for {user['username']}.")
                    st.rerun()

            with reset_col.popover("Reset password", use_container_width=True):
                with st.form(f"admin-reset-password-form-{user['id']}"):
                    reset_password_value = st.text_input("New temporary password", type="password", key=f"admin-reset-password-{user['id']}")
                    reset_password_confirm = st.text_input("Confirm temporary password", type="password", key=f"admin-reset-confirm-{user['id']}")
                    reset_submitted = st.form_submit_button("Reset and sign out")
                if reset_submitted:
                    if reset_password_value != reset_password_confirm:
                        st.error("The temporary passwords do not match.")
                    else:
                        try:
                            reset_user_password(str(user["id"]), reset_password_value, actor=CURRENT_USER)
                        except (OSError, sqlite3.Error, ValueError) as exc:
                            st.error(f"Could not reset the password: {exc}")
                        else:
                            st.success(f"Password reset for {user['username']}; their current sessions were signed out.")
                            st.rerun()

            action_label = "Restore account" if user["disabled"] else "Disable account"
            if account_col.button(action_label, key=f"admin-toggle-account-{user['id']}"):
                try:
                    set_user_disabled(str(user["id"]), not bool(user["disabled"]), actor=CURRENT_USER)
                except (OSError, sqlite3.Error, ValueError) as exc:
                    st.error(str(exc))
                else:
                    st.rerun()

    st.caption(f"Passwords require at least {password_min_length()} characters and are stored as one-way hashes.")


def _render_projects() -> None:
    st.title("Projects")
    organization = load_organization(STUDIO_WORKSPACE)
    st.subheader(organization.get("name", "Molecular Design Studio"))
    st.caption("Projects group jobs from different tools around a shared target or research question.")
    last_created = st.session_state.pop("projects-last-created", None)
    if last_created:
        st.success(f"Project created: {last_created}")
    if st.button("New project", icon=":material/add:", type="primary", key="projects-new"):
        _render_project_dialog()

    projects = _studio_projects()
    if not projects:
        st.info("No projects yet. Create a project to group jobs from multiple tools.")
        return

    for start in range(0, len(projects), 3):
        columns = st.columns(3, gap="medium")
        for column, project in zip(columns, projects[start : start + 3]):
            with column:
                color = str(project.get("color") or "#75817b")
                if not re.fullmatch(r"#[0-9A-Fa-f]{6}", color):
                    color = "#75817b"
                st.markdown(f'<span style="color:{color};font-size:1.4rem">●</span>', unsafe_allow_html=True)
                st.markdown(f"**{project['name']}**")
                st.caption(project.get("description") or "No description")
                try:
                    count = search_jobs(project_id=str(project["id"]), limit=1)["total"]
                except Exception:
                    count = "—"
                st.caption(f"{count} job{'s' if count != 1 else ''} · Created {_job_submitted_label(project.get('created_at'))}")
                st.divider()
    st.caption(f"Project metadata is stored in the Studio workspace at `{STUDIO_WORKSPACE}`. Job files remain in the Results folder configured in Settings.")


def _render_settings() -> None:
    st.title("Settings")
    st.subheader("Runtime folders")
    settings = load_settings(STUDIO_WORKSPACE)
    st.caption(
        f"Job catalog: `{catalog_path()}`. It stores searchable job details and links to files; run files remain in the Results folder."
    )
    if st.session_state.pop("settings-storage-saved", False):
        st.success("Storage settings saved.")
    results_directory = st.text_input(
        "Results folder",
        value=settings["results_directory"],
        help="Molecular Design Studio reads and stores job runs in this folder. The default is a new empty folder on the same data drive as the previous app.",
        key="settings-results-directory",
    )
    reference_directory = st.text_input(
        "Model reference files folder",
        value=settings["reference_directory"],
        help="Shared, mostly read-only model weights and reference files. Point this at an existing reference folder to reuse large checkpoints without copying them into the Studio workspace.",
        key="settings-reference-directory",
    )
    if st.button("Save storage settings", type="primary", key="settings-save-storage"):
        try:
            save_settings(
                {
                    "results_directory": results_directory,
                    "reference_directory": reference_directory,
                },
                STUDIO_WORKSPACE,
            )
            ensure_studio_workspace()
            st.session_state["settings-storage-saved"] = True
            st.rerun()
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            st.error(f"Could not save results settings: {exc}")
    st.code(
        f"App home: {app_home()}\nRuns: {runs_root()}\nReferences: {reference_root()}",
        language="text",
    )
    st.subheader("Local API")
    st.markdown(f"API address: `{API_BASE_URL}` · [Interactive API docs]({API_BASE_URL}/docs)")
    api_address = st.text_input("Workflow API address", value=settings["api_base_url"], key="settings-api-address", disabled=bool(os.environ.get("MDS_API_BASE_URL")))
    if os.environ.get("MDS_API_BASE_URL"):
        st.caption("The API address is supplied by the launch environment.")
    elif st.button("Save API address", key="settings-save-api"):
        try:
            save_settings({**settings, "api_base_url": api_address}, STUDIO_WORKSPACE)
            st.rerun()
        except (OSError, ValueError) as exc:
            st.error(str(exc))
    st.caption("The API currently accepts requests from this computer only.")

    st.subheader("MN Compute Scheduler")
    st.caption(
        "Molecular Design Studio sends new workflow runs to the shared scheduler. "
        "This check is manual and read-only; it does not submit jobs or change scheduler state."
    )
    if st.button("Check scheduler connection", key="settings-check-scheduler"):
        try:
            from mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler import scheduler_status

            scheduler = scheduler_status()
            st.success(f"Scheduler is online at {scheduler['url']}.")
            st.caption(f"Molecular Design Studio workers online: {len(scheduler['workers'])}.")
        except Exception as exc:
            st.error(f"Could not connect to the scheduler: {exc}")

    from mn_molecular_design_studio.cluster.ui import settings as cluster_settings
    cluster_settings(CURRENT_USER, str(st.session_state.get("studio-session-token") or "development-mode"))


page = st.session_state["page"]
if page == "All tools":
    _render_tool_catalog()
elif page == "Workflow setup":
    _render_workflow_setup()
elif page == "My runs":
    _render_runs()
elif page == "Job details":
    _render_job_details()
elif page == "Database":
    _render_database()
elif page == "Pipelines":
    _render_pipelines()
elif page == "Projects":
    _render_projects()
elif page == "Users & access" and CURRENT_USER.get("role") == "admin":
    _render_user_access()
elif page == "Account":
    _render_account()
elif page == "Settings" and CURRENT_USER.get("role") == "admin":
    _render_settings()
else:
    st.session_state["page"] = "All tools"
    _render_tool_catalog()
