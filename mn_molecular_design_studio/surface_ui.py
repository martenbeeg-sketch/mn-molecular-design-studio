"""Input and result views for the standalone protein-surface tools."""

from __future__ import annotations

import hashlib
from html import escape
import json
import math
from pathlib import Path
import re

import altair as alt
import pandas as pd
import streamlit as st

from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
    StructureVisualization,
    molstar_custom_component,
)
from mn_molecular_design_studio.engine.workflow_engine.core.structures import STANDARD_PROTEIN_RESIDUES, filter_pdb_text, residue_inventory
from mn_molecular_design_studio.target_editor import (
    _sequence_entries_by_chain,
    render_target_source,
)


SURFACE_TOOLS = {
    "detection.pesto": {
        "tool": "pesto",
        "title": "PeSTo",
        "description": "Predict protein interface scores for protein, DNA/RNA, ion, ligand, and lipid interactions.",
        "version": "PeSTo · ba651aa",
        "selection": "multiple",
        "about": "PeSTo assigns residue-level interface probabilities for protein, DNA/RNA, ion, ligand, or lipid interactions. Select the protein chains and the interaction class to score.",
    },
    "detection.scannet": {
        "tool": "scannet",
        "title": "ScanNet",
        "description": "Identify likely protein interaction and hotspot residues.",
        "version": "ScanNet · a61623c",
        "selection": "multiple",
        "about": "ScanNet predicts binding-site probabilities from a protein structure. Select one or more chains, then choose interface, epitope, or disordered-region prediction. MSA-assisted prediction is optional.",
    },
    "detection.masif_seed": {
        "tool": "masif_seed",
        "title": "MaSIF-seed",
        "description": "Predict a protein interaction site on one molecular surface.",
        "version": "MaSIF-seed · 7715391",
        "selection": "single",
        "about": "MaSIF-seed searches a single prepared chain for a likely interaction site. The result includes a predicted surface mesh and per-vertex scores, which are projected onto nearby residues for the structure view.",
    },
    "detection.surf2spot": {
        "tool": "surf2spot",
        "title": "Surf2Spot",
        "description": "Map predicted protein hotspot regions onto a molecular surface.",
        "version": "Surf2Spot · a38a7a7",
        "selection": "multiple",
        "about": "Surf2Spot predicts hotspot regions from the selected protein chains. Results include hotspot tables and surface files, plus the tool's PyMOL session when available.",
    },
    "detection.surface_analysis": {
        "tool": "surface_analysis",
        "title": "Protein Surface Analysis",
        "description": "Run four protein-surface predictors together and compare their PPI and hotspot evidence.",
        "version": "PeSTo · ScanNet · Surf2Spot · MaSIF-seed",
        "selection": "multiple",
        "about": "Runs PeSTo, ScanNet, Surf2Spot, and MaSIF-seed against one selected target. The results begin with a cross-tool residue comparison, followed by each method's full tables and visualizations.",
    },
}


_STYLE = """<style>
.block-container{max-width:1240px!important;padding-top:2.3rem!important}
.st-key-surface-tool-surface{--surface-green:#27674e;--surface-mint:#edf6f0}
.st-key-surface-tool-surface .surface-heading{text-align:center;padding:.3rem 0 1.15rem}
.st-key-surface-tool-surface .surface-eyebrow{font-size:.68rem;letter-spacing:.12em;font-weight:650;color:#698575}
.st-key-surface-tool-surface h1{font-size:2.05rem;color:#173d31;letter-spacing:-.035em;margin:.45rem 0}
.st-key-surface-tool-surface .surface-heading p{font-size:.92rem;color:#718478;margin:.55rem auto;max-width:620px}
.st-key-surface-tool-surface .surface-version{color:#7c9284;font-size:.75rem}
.st-key-surface-tool-surface [data-testid="stVerticalBlockBorderWrapper"]{border-color:#d7e6dc!important;border-radius:12px!important}
.st-key-surface-tool-surface .surface-section{font-size:.72rem;font-weight:650;letter-spacing:.09em;color:#5b7b67;margin:.1rem 0 .8rem}
.st-key-surface-tool-surface [data-baseweb="input"],.st-key-surface-tool-surface [data-baseweb="textarea"],.st-key-surface-tool-surface [data-baseweb="select"]>div{background:#f4f8f5!important;border-color:#d7e6dc!important;border-radius:7px!important}
.st-key-surface-tool-surface [data-testid="stExpander"]{border-color:#d7e6dc!important;border-radius:8px!important;background:#fcfefc}
.st-key-surface-tool-surface [data-testid="stFileUploaderDropzone"]{background:#f3f8f5!important;border:1px dashed #b9d2c2!important;border-radius:8px}
.st-key-surface-tool-surface [data-testid="stButton"] button{border-radius:999px!important;border-color:#bad2c3;color:#27674e}
.st-key-surface-tool-surface [data-testid="stButton"] button[kind="primary"]{background:#2f6b50!important;border-color:#2f6b50!important;color:white!important}
.st-key-surface-tool-surface [data-testid="stButton"] button[kind="primary"]:hover{background:#24553e!important;border-color:#24553e!important}
.st-key-surface-tool-surface [data-testid="stButton"] button:disabled{opacity:.58}
.st-key-surface-tool-surface [role="radio"][aria-checked="true"]{background:#2f6b50!important;border-color:#2f6b50!important;color:white!important}
.st-key-surface-tool-surface [data-baseweb="tag"]{background:#e0eee5!important;color:#27674e!important}
</style>"""

_SCORE_COLOR_PARAMS = {
    "domain": [0, 100],
    "list": {"colors": [0xD72323, 0xF5F5F5, 0x1952E6]},
}


def _slug(workflow_id: str) -> str:
    return re.sub(r"[^a-z0-9-]+", "-", workflow_id.casefold()).strip("-")


def _filtered_target_path(source_path: Path, chains: list[str], *, owner_id: str, label: str) -> Path:
    """Write a user-owned, chain-filtered PDB when a tool needs a single target file."""
    text = source_path.read_text(encoding="utf-8", errors="replace")
    filtered = filter_pdb_text(text, keep_chains=set(chains), remove_waters=True, remove_hetero=True)
    digest = hashlib.sha256(filtered.encode("utf-8")).hexdigest()[:12]
    safe_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("._-") or "target"
    path = user_input_root(owner_id) / "structures" / f"{safe_label}-{digest}.pdb"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(filtered, encoding="utf-8")
    return path


def render_setup(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    workflow_id = str(schema.get("workflow") or "")
    config = SURFACE_TOOLS.get(workflow_id)
    if config is None:
        st.error("This surface-analysis workflow is not configured in the Studio UI.")
        return
    slug = _slug(workflow_id)

    st.html(_STYLE)
    with st.container(key="surface-tool-surface"):
        if st.button("← Back to All tools", key=f"{slug}-back"):
            st.session_state["page"] = "All tools"
            st.session_state.pop("selected_workflow", None)
            st.rerun()
        st.html(
            '<div class="surface-heading"><div class="surface-eyebrow">PROTEIN SURFACE ANALYSIS</div>'
            f'<h1>{escape(config["title"])}</h1>'
            f'<p>{escape(config["description"])}</p>'
            f'<span class="surface-version">{escape(config["version"])}</span></div>'
        )
        with st.expander("About this tool", expanded=False):
            st.write(config["about"])

        st.markdown('<div class="surface-section">TARGET STRUCTURE</div>', unsafe_allow_html=True)
        st.caption("Fetch by PDB accession or upload a PDB file.")
        source = render_target_source(
            key_prefix=f"{slug}-target",
            label="Target",
            example_code="7JZU",
            allow_cif=False,
            owner_id=owner_id,
            upload_label="Drop a PDB file here or browse",
            upload_help="Upload a PDB (.pdb/.ent) structure, or fetch it by PDB accession above.",
            upload_expander_label="Upload a PDB file",
        )

        selected_chains: list[str] = []
        target_path = source.path
        if source.path is not None and source.path.is_file():
            try:
                structure_text = source.path.read_text(encoding="utf-8", errors="replace")
                sequence_by_chain = _sequence_entries_by_chain(structure_text)
                chain_options = sorted(sequence_by_chain)
            except (OSError, ValueError) as exc:
                st.error(f"Could not read the target chains: {exc}")
                chain_options = []
                structure_text = ""

            if chain_options:
                chain_key = f"{slug}-chains-{hashlib.sha1(str(source.path).encode()).hexdigest()[:10]}"
                if config["selection"] == "single":
                    default_chain = chain_options[0]
                    if st.session_state.get(chain_key) not in chain_options:
                        st.session_state[chain_key] = default_chain
                    selected = st.selectbox("Target chain", chain_options, key=chain_key)
                    selected_chains = [str(selected)] if selected else []
                else:
                    if st.session_state.get(chain_key) is None:
                        st.session_state[chain_key] = chain_options
                    st.session_state[chain_key] = [
                        chain for chain in st.session_state[chain_key] if chain in chain_options
                    ]
                    selected_chains = st.multiselect("Target chains", chain_options, key=chain_key)
                residues = sum(len(sequence_by_chain.get(chain, [])) for chain in selected_chains)
                st.caption(f"{len(selected_chains)} chain(s) selected · {residues:,} residues")
                if selected_chains:
                    view_text = filter_pdb_text(structure_text, keep_chains=set(selected_chains))
                    st.markdown("#### MolView")
                    molstar_custom_component(
                        [StructureVisualization(
                            pdb=view_text,
                            color="chain-id",
                            representation_type="molecular-surface",
                        )],
                        key=f"{slug}-input-view-{hashlib.sha1((str(source.path)+','.join(selected_chains)).encode()).hexdigest()[:12]}",
                        height=450,
                        show_controls=True,
                        show_structure_tools=False,
                    )

                if workflow_id in {"detection.surf2spot", "detection.masif_seed"} and selected_chains:
                    target_path = _filtered_target_path(
                        source.path,
                        selected_chains,
                        owner_id=owner_id,
                        label=f"{slug}-{source.name or source.path.stem}",
                    )

        if workflow_id == "detection.surf2spot" and selected_chains:
            st.caption("Surf2Spot will analyze the selected chains from the target structure.")
        if workflow_id == "detection.masif_seed":
            st.caption("MaSIF-seed analyzes the selected chain as a single target surface.")
        if workflow_id == "detection.pesto":
            st.caption("PeSTo scores the selected protein chains. The interface class selects the interaction type to predict; partner molecules are not required.")
        elif workflow_id == "detection.scannet":
            st.caption("Only the selected chains are included in the prediction.")
        elif workflow_id == "detection.surface_analysis":
            st.caption("All methods use the selected protein chains. PeSTo's interface class selects the interaction type to predict; partner molecules are not required.")

        pesto_interface = "protein"
        if workflow_id in {"detection.pesto", "detection.surface_analysis"}:
            st.markdown("**Interface class**")
            interface_labels = {
                "protein": "Protein–protein",
                "dna_rna": "DNA / RNA",
                "ion": "Ion",
                "ligand": "Ligand",
                "lipid": "Lipid",
            }
            interface_key = f"{slug}-interface"
            current_interface = str(st.session_state.get(interface_key, "protein"))
            if current_interface not in interface_labels:
                current_interface = "protein"
                st.session_state[interface_key] = current_interface
            interface_columns = st.columns(len(interface_labels), gap="small")
            for column, interface_value in zip(interface_columns, interface_labels):
                with column:
                    if st.button(
                        interface_labels[interface_value],
                        key=f"{slug}-interface-button-{interface_value}",
                        type="primary" if current_interface == interface_value else "secondary",
                        use_container_width=True,
                    ):
                        st.session_state[interface_key] = interface_value
                        st.rerun()
            pesto_interface = str(st.session_state.get(interface_key, "protein"))
            st.caption(f"Score protein residues for {interface_labels[pesto_interface].casefold()} interaction propensity.")

        scan_mode = "interface"
        scannet_use_msa = False
        if workflow_id in {"detection.scannet", "detection.surface_analysis"}:
            st.markdown("**Prediction mode**")
            mode_labels = {"interface": "Interface", "epitope": "Epitope", "idp": "Disordered region"}
            mode_key = f"{slug}-mode"
            current_mode = str(st.session_state.get(mode_key, "interface"))
            if current_mode not in mode_labels:
                current_mode = "interface"
                st.session_state[mode_key] = current_mode
            mode_columns = st.columns(3, gap="small")
            for column, mode_value in zip(mode_columns, mode_labels):
                with column:
                    if st.button(
                        mode_labels[mode_value],
                        key=f"{slug}-mode-button-{mode_value}",
                        type="primary" if current_mode == mode_value else "secondary",
                        use_container_width=True,
                    ):
                        st.session_state[mode_key] = mode_value
                        st.rerun()
            scan_mode = str(st.session_state.get(mode_key, "interface"))
            st.caption({
                "interface": "Predict residue-level protein interaction interface scores.",
                "epitope": "Predict epitope residues on the selected protein chains.",
                "idp": "Predict interaction regions in intrinsically disordered proteins.",
            }[scan_mode])
            with st.expander("Advanced settings", expanded=False):
                use_msa = st.checkbox(
                    "Use MSA",
                    value=False,
                    key=f"{slug}-use-msa",
                    help="Slower and requires HHblits/database paths inside the ScanNet container.",
                )
                scannet_use_msa = use_msa
        else:
            use_msa = False

        masif_chain_id = selected_chains[0] if selected_chains else ""
        if workflow_id == "detection.surface_analysis" and selected_chains:
            masif_chain_id = st.selectbox(
                "MaSIF-seed chain",
                selected_chains,
                index=selected_chains.index(masif_chain_id) if masif_chain_id in selected_chains else 0,
                key=f"{slug}-masif-chain",
                help="MaSIF-seed searches one chain surface; the other methods use all selected chains.",
            )

        st.markdown('<div class="surface-section">SAVE AND RUN</div>', unsafe_allow_html=True)
        project_id = project_selector(f"surface-{slug}")
        prioritize = st.checkbox(
            "Prioritize in scheduler queue",
            key=f"{slug}-prioritize",
            help="Moves this run ahead of normal queued work when resources become available.",
        )

        ready = bool(target_path and selected_chains and project_id)
        submitted = st.button(
            f"Run {config['title']}",
            type="primary",
            use_container_width=True,
            disabled=not ready,
            key=f"{slug}-run",
        )
        if submitted and target_path and project_id:
            parameters: dict[str, object] = {"target_pdb": str(target_path)}
            if workflow_id in {"detection.pesto", "detection.scannet", "detection.surface_analysis"}:
                parameters["chain_ids"] = selected_chains
            if workflow_id == "detection.surface_analysis":
                parameters.update(
                    pesto_interface=pesto_interface,
                    scannet_mode=scan_mode,
                    scannet_use_msa=scannet_use_msa,
                    masif_chain_id=masif_chain_id,
                )
            if workflow_id == "detection.pesto":
                parameters["interface"] = pesto_interface
            elif workflow_id == "detection.masif_seed":
                parameters["chain_id"] = selected_chains[0]
            if workflow_id == "detection.scannet":
                parameters["mode"] = scan_mode
                parameters["use_msa"] = use_msa
            submit(schema, workflow_id, parameters, project_id, {"prioritize": prioritize}, config["title"])


def _artifact_paths(artifacts: list[tuple[dict, Path]], suffixes: set[str]) -> list[tuple[dict, Path]]:
    return [(item, path) for item, path in artifacts if path.suffix.casefold() in suffixes]


def _target_structure(tool: str, run_dir: Path, input_payload: dict) -> Path | None:
    candidates = {
        "surface_analysis": [run_dir / "artifacts/surface_analysis/input/target.pdb"],
        "pesto": [run_dir / "artifacts/pesto/input/target.pdb"],
        "scannet": [run_dir / "artifacts/input/target.pdb"],
        "surf2spot": [run_dir / "artifacts/surf2spot/source/target.pdb", run_dir / "artifacts/surf2spot/input/target.pdb"],
        "masif_seed": [run_dir / "artifacts/masif_seed/input/target.pdb"],
    }.get(tool, [])
    raw_parameters = input_payload.get("parameters") if isinstance(input_payload.get("parameters"), dict) else {}
    raw = raw_parameters.get("target_pdb") or input_payload.get("target_pdb")
    if raw:
        candidates.append(Path(str(raw)).expanduser())
    return next((path.resolve() for path in candidates if path.is_file()), None)


def _residue_token(value: object) -> str:
    text = str(value).strip()
    match = re.fullmatch(r"(-?\d+)(?:\.0+)?([A-Za-z]?)", text)
    return f"{int(match.group(1))}{match.group(2)}" if match else text


def _residue_number(value: object) -> int | None:
    match = re.match(r"^\s*(-?\d+)", str(value))
    return int(match.group(1)) if match else None


def _normalized_score_table(raw: pd.DataFrame, tool: str, target_pdb: Path | None, *, source_path: Path | None = None) -> pd.DataFrame | None:
    if raw.empty:
        return None
    normalized = {re.sub(r"[^a-z0-9]+", "_", str(column).casefold()).strip("_"): column for column in raw.columns}
    chain_column = next((normalized[name] for name in ("chain", "chain_id", "chainid", "chain_name") if name in normalized), None)
    residue_column = next((normalized[name] for name in ("residue", "residue_index", "residue_id", "resnum", "residue_number", "aa_id") if name in normalized), None)
    amino_column = next((normalized[name] for name in ("amino_acid", "sequence", "aa", "resname") if name in normalized), None)
    score_column = next((normalized[name] for name in (
        "score", "binding_site_probability", "probability", "interface_probability",
        "ppi_probability", "binding_probability", "prediction_score", "hotspot_score",
    ) if name in normalized), None)
    if score_column is None:
        score_column = next((column for key, column in normalized.items() if any(word in key for word in ("score", "probability", "prob", "confidence"))), None)
    if tool == "surf2spot" and residue_column is None and "aa_id" in normalized:
        residue_column = normalized["aa_id"]
    if residue_column is None or score_column is None:
        return None
    if chain_column is None:
        if tool != "surf2spot" or target_pdb is None:
            return None
        try:
            chain_options = sorted(_sequence_entries_by_chain(target_pdb.read_text(encoding="utf-8", errors="replace")))
        except (OSError, ValueError):
            chain_options = []
        # Surf2Spot emits one table per chain, named target_<chain>[...].csv.
        matching_chains = [chain for chain in chain_options if source_path and re.match(rf"^target_{re.escape(chain)}(?:_|$)", source_path.stem)]
        if len(matching_chains) == 1:
            chain_values = matching_chains[0]
        elif len(chain_options) == 1:
            chain_values = chain_options[0]
        else:
            return None
    else:
        chain_values = raw[chain_column].fillna("_").astype(str).str.strip().replace("", "_")
    residue_labels = raw[residue_column].map(_residue_token)
    if tool == "surf2spot" and residue_column == normalized.get("aa_id"):
        if target_pdb is None:
            return None
        # aa_id is a one-based index in Surf2Spot's renumbered chain, not an
        # author residue ID. Restore author numbering, gaps, and insertion codes.
        by_chain: dict[str, list[str]] = {}
        for row in residue_inventory(target_pdb.read_text(encoding="utf-8", errors="replace")):
            if row["resname"] in STANDARD_PROTEIN_RESIDUES:
                by_chain.setdefault(row["chain"], []).append(f"{row['residue_number']}{row['insertion_code']}")
        row_chains = [chain_values] * len(raw) if isinstance(chain_values, str) else list(chain_values)
        mapped = []
        for chain, token in zip(row_chains, residue_labels):
            residues = by_chain.get(chain, [])
            index = int(token) if re.fullmatch(r"\d+", token) else 0
            mapped.append(residues[index - 1] if 1 <= index <= len(residues) else "")
        residue_labels = pd.Series(mapped, index=raw.index)
    table = pd.DataFrame({
        "Chain": chain_values,
        "Residue": residue_labels,
        "Residue Number": pd.to_numeric(residue_labels.map(_residue_number), errors="coerce"),
        "Amino acid": raw[amino_column].astype(str) if amino_column else "",
        "Score": pd.to_numeric(raw[score_column], errors="coerce"),
    }).dropna(subset=["Residue Number", "Score"])
    table = table[table["Score"].map(math.isfinite)]
    if table.empty:
        return None
    table["Residue Number"] = table["Residue Number"].astype(int)
    table["Residue ID"] = table["Chain"] + table["Residue"].astype(str)
    table = table.sort_values("Score", ascending=False).drop_duplicates(["Chain", "Residue"], keep="first")
    return table.drop(columns=["Residue Number"])


def _combined_surface_consensus(tables: list[pd.DataFrame]) -> pd.DataFrame:
    if not tables:
        return pd.DataFrame()
    combined = pd.concat(tables, ignore_index=True)
    # Tools can emit multiple tables for a target. Count each method once per
    # residue and retain its strongest score before computing within-tool ranks.
    combined = combined.groupby(["Method", "Chain", "Residue", "Residue ID"], as_index=False)["Score"].max()
    combined["Within-method percentile"] = combined.groupby("Method")["Score"].rank(pct=True)
    return (
        combined.groupby(["Chain", "Residue", "Residue ID"], as_index=False)
        .agg(**{
            "Methods": ("Method", "nunique"),
            "Scoring methods": ("Method", lambda values: ", ".join(sorted(set(values)))),
            "Mean within-method percentile": ("Within-method percentile", "mean"),
        })
        .sort_values(["Methods", "Mean within-method percentile", "Chain", "Residue"], ascending=[False, False, True, True])
    )


def _masif_score_table(run_dir: Path, target_pdb: Path | None) -> pd.DataFrame | None:
    if target_pdb is None:
        return None
    surfaces = sorted((run_dir / "artifacts/masif_seed/pred_surfaces").glob("*.ply"))
    predictions = sorted((run_dir / "artifacts/masif_seed/pred_data").glob("*.npy"))
    if not surfaces or not predictions:
        return None
    try:
        from mn_molecular_design_studio.engine.workflow_engine.core.detection_annotations import _masif_residue_scores

        table = _masif_residue_scores(surfaces[0], predictions[0], target_pdb)
    except Exception:
        return None
    if table.empty:
        return None
    table = table.rename(columns={"chain": "Chain", "residue": "Residue", "amino_acid": "Amino acid", "score": "Score"})
    table["Chain"] = table["Chain"].astype(str).str.strip()
    table["Residue"] = table["Residue"].map(_residue_token)
    table["Residue ID"] = table["Chain"] + table["Residue"].astype(str)
    return table.sort_values("Score", ascending=False)


def _pdb_with_residue_scores(pdb_text: str, table: pd.DataFrame, *, maximum: float) -> str:
    scores = {
        (str(row.Chain).strip(), _residue_token(row.Residue)): float(row.Score)
        for row in table.itertuples()
        if pd.notna(row.Residue) and pd.notna(row.Score)
    }
    scale = 100.0 if maximum <= 1.0 else 100.0 / maximum
    lines: list[str] = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            lines.append(line)
            continue
        chain = line[21].strip() or "_"
        residue = _residue_token(line[22:26].strip() + line[26].strip())
        score = max(0.0, min(100.0, scores.get((chain, residue), 0.0) * scale))
        padded = line if len(line) >= 66 else line.ljust(66)
        lines.append(f"{padded[:60]}{score:6.2f}{padded[66:]}")
    return "\n".join(lines) + "\n"


def _ordered_score_table(table: pd.DataFrame) -> pd.DataFrame:
    return table.assign(_number=table["Residue"].map(_residue_number)).sort_values(["Chain", "_number", "Residue"]).drop(columns="_number")


def _render_score_profile(table: pd.DataFrame) -> None:
    ordered = _ordered_score_table(table)
    chart = alt.Chart(ordered).mark_line(color="#39785b").encode(
        x=alt.X("Residue ID:N", sort=ordered["Residue ID"].tolist(), title="Residue"),
        y=alt.Y("Score:Q"),
        detail="Chain:N",
        tooltip=[alt.Tooltip("Residue ID:N"), alt.Tooltip("Score:Q", format=".3f")],
    ).properties(height=220)
    st.altair_chart(chart, use_container_width=True)


def _masif_point_cloud(run_dir: Path, *, limit: int = 6000) -> str | None:
    surface_paths = sorted((run_dir / "artifacts/masif_seed/pred_surfaces").glob("*.ply"))
    prediction_paths = sorted((run_dir / "artifacts/masif_seed/pred_data").glob("*.npy"))
    if not surface_paths:
        return None
    try:
        lines = surface_paths[0].read_text(encoding="utf-8", errors="replace").splitlines()
        header_end = lines.index("end_header")
        vertex_count = next(int(line.split()[-1]) for line in lines[:header_end] if line.startswith("element vertex "))
        property_names: list[str] = []
        in_vertices = False
        for line in lines[:header_end]:
            if line.startswith("element vertex "):
                in_vertices = True
            elif line.startswith("element "):
                in_vertices = False
            elif in_vertices and line.startswith("property "):
                property_names.append(line.split()[-1])
        property_index = {name: index for index, name in enumerate(property_names)}
        xyz_indices = [property_index[name] for name in ("x", "y", "z")]
        vertices = [[float(value) for value in lines[header_end + 1 + index].split()] for index in range(vertex_count)]
        import numpy as np

        if prediction_paths:
            values = np.asarray(np.load(prediction_paths[0]), dtype=float).reshape(-1)
        elif "iface" in property_index:
            values = np.asarray([row[property_index["iface"]] for row in vertices], dtype=float)
        else:
            return None
        if len(vertices) != len(values) or not len(values):
            return None
        max_score = float(np.nanmax(values))
        if not np.isfinite(max_score):
            return None
        scale = 100.0 if max_score <= 1.0 else 100.0 / max_score
        stride = max(1, (len(vertices) + limit - 1) // limit)
        point_lines: list[str] = []
        for index in range(0, len(vertices), stride):
            row = vertices[index]
            x, y, z = (row[position] for position in xyz_indices)
            score = max(0.0, min(100.0, float(values[index]) * scale))
            serial = len(point_lines) + 1
            point_lines.append(
                f"HETATM{serial:5d}  C   PNT A{serial:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00{score:6.2f}           C"
            )
        return "\n".join([*point_lines, "END", ""])
    except Exception:
        return None


def render_results(
    tool: str,
    run_id: str,
    run_dir: Path,
    artifacts: list[tuple[dict, Path]],
    input_payload: dict,
) -> None:
    """Render tool output tables, predicted residues, structures, and surface downloads."""
    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    if tool == "pesto":
        interface_labels = {
            "protein": "Protein–protein",
            "dna_rna": "DNA / RNA",
            "ion": "Ion",
            "ligand": "Ligand",
            "lipid": "Lipid",
        }
        interface = str(params.get("interface") or params.get("mode") or "protein")
        st.caption(f"PeSTo interface class · {interface_labels.get(interface, interface)}")
    elif tool == "scannet":
        mode = str(params.get("mode") or "interface")
        st.caption(f"ScanNet prediction mode · {mode.replace('_', ' ').title()}")
    target = _target_structure(tool, run_dir, input_payload)
    csv_artifacts = _artifact_paths(artifacts, {".csv"})
    pdb_artifacts = _artifact_paths(artifacts, {".pdb", ".ent"})
    surface_artifacts = _artifact_paths(artifacts, {".ply"})
    session_artifacts = [pair for pair in artifacts if pair[1].suffix.casefold() == ".pse"]

    score_table: pd.DataFrame | None = None
    selected_csv_path: Path | None = None
    raw_table: pd.DataFrame | None = None
    if csv_artifacts:
        selected_index = st.selectbox(
            "Prediction table",
            range(len(csv_artifacts)),
            format_func=lambda index: str(csv_artifacts[index][0].get("relative_path") or csv_artifacts[index][1].name),
            key=f"surface-result-csv-{run_id}",
        )
        selected_item, selected_csv_path = csv_artifacts[selected_index]
        try:
            raw_table = pd.read_csv(selected_csv_path, low_memory=False, nrows=50_000)
            score_table = _normalized_score_table(raw_table, tool, target, source_path=selected_csv_path)
            st.caption(f"Output table: `{selected_item.get('relative_path') or selected_csv_path.name}`")
            if score_table is not None:
                summary_cols = st.columns(3)
                summary_cols[0].metric("Residues scored", f"{len(score_table):,}")
                summary_cols[1].metric("Top residue", str(score_table.iloc[0]["Residue ID"]))
                summary_cols[2].metric("Top score", f"{float(score_table.iloc[0]['Score']):.3f}")
                _render_score_profile(score_table)
                with st.expander("Top predicted residues", expanded=True):
                    st.dataframe(score_table.head(20), hide_index=True, use_container_width=True)
            st.dataframe(raw_table, hide_index=True, use_container_width=True)
            st.download_button(
                "Download prediction table",
                data=selected_csv_path.read_bytes(),
                file_name=selected_csv_path.name,
                mime="text/csv",
                key=f"surface-result-csv-download-{run_id}-{selected_index}",
            )
        except (OSError, ValueError, pd.errors.ParserError, UnicodeError) as exc:
            st.warning(f"Could not read this prediction table: {exc}")

    if tool == "masif_seed" and score_table is None:
        score_table = _masif_score_table(run_dir, target)
        if score_table is not None:
            summary_cols = st.columns(3)
            summary_cols[0].metric("Residues with surface scores", f"{len(score_table):,}")
            summary_cols[1].metric("Top residue", str(score_table.iloc[0]["Residue ID"]))
            summary_cols[2].metric("Top score", f"{float(score_table.iloc[0]['Score']):.3f}")
            _render_score_profile(score_table)

    highlighted: list[str] = []
    if score_table is not None:
        finite_scores = pd.to_numeric(score_table["Score"], errors="coerce")
        finite_scores = finite_scores[finite_scores.map(pd.notna)]
        score_maximum = float(finite_scores.max()) if not finite_scores.empty else 1.0
        score_maximum = score_maximum if score_maximum > 0 else 1.0
        score_default = min(0.5, score_maximum)
        threshold = st.slider(
            "Minimum score to highlight",
            min_value=0.0,
            max_value=score_maximum,
            value=score_default,
            step=max(score_maximum / 100.0, 0.001),
            key=f"surface-result-threshold-{run_id}",
            help="Filters the residue list and controls which predicted sites are highlighted in the 3D view.",
        )
        visible_scores = score_table[pd.to_numeric(score_table["Score"], errors="coerce") >= threshold].copy()
        st.caption(f"{len(visible_scores):,} residue(s) meet the selected score threshold.")
        if tool == "surf2spot":
            st.caption("Surf2Spot's native aa_id indices are mapped back to the input structure's chain and residue IDs in the residue scores and MolView.")
        st.download_button(
            "Download mapped residue scores",
            score_table.to_csv(index=False).encode("utf-8"),
            file_name=f"{tool}-{run_id}-residue-scores.csv",
            mime="text/csv",
            key=f"surface-result-mapped-scores-{run_id}",
        )
        if raw_table is None:
            st.dataframe(visible_scores.head(100), hide_index=True, use_container_width=True)
        insertion_positions = {
            (str(row.Chain), _residue_number(row.Residue)) for row in score_table.itertuples()
            if re.search(r"[A-Za-z]$", str(row.Residue))
        }
        # Mol*'s selection helper accepts integer IDs. Exact insertion-code
        # scores remain visible through the B-factor coloring without selecting
        # another residue sharing the same integer ID.
        highlighted = [
            f"{row.Chain}{_residue_number(row.Residue)}"
            for row in visible_scores.head(300).itertuples()
            if str(row.Chain).strip() and _residue_number(row.Residue) is not None
            and (str(row.Chain), _residue_number(row.Residue)) not in insertion_positions
        ]

    preview_structure = next(
        (path for _item, path in pdb_artifacts if path != target and "scored" in path.name.casefold()),
        None,
    )
    if preview_structure is None:
        preview_structure = next((path for _item, path in pdb_artifacts if path != target), None)
    if score_table is not None and target is not None:
        preview_structure = target
    if preview_structure is None:
        preview_structure = target

    if preview_structure is not None and preview_structure.is_file():
        st.subheader("MolView")
        caption = "Predicted high-scoring residues are highlighted on the molecular surface." if highlighted else "Molecular surface of the analyzed target."
        show_surface = st.checkbox(
            "Show molecular surface",
            value=True,
            key=f"surface-result-show-surface-{run_id}",
            help="Toggle surface rendering. Residues and surface points are colored by their prediction score.",
        )
        representation = "molecular-surface+ball-and-stick" if show_surface else "ball-and-stick"
        st.caption(
            "Mol* reads the residue prediction scores from the PDB B-factor field. Red marks stronger scores."
            if score_table is not None else caption
        )
        if score_table is not None:
            st.markdown(
                '<div style="display:flex;align-items:center;gap:10px;margin:4px 0 14px;max-width:520px;">'
                '<span style="font-size:12px;color:#59665f;">lower score</span>'
                '<div style="height:12px;flex:1;border-radius:6px;border:1px solid #d1d5db;'
                'background:linear-gradient(90deg,rgb(25,82,230),rgb(245,245,245),rgb(215,35,35));"></div>'
                '<span style="font-size:12px;color:#59665f;">higher score</span></div>',
                unsafe_allow_html=True,
            )
        try:
            if preview_structure.stat().st_size > 50 * 1024 * 1024:
                st.info("This structure exceeds the 50 MB preview limit. Download it from the Files tab to inspect it locally.")
            else:
                pdb_text = preview_structure.read_text(encoding="utf-8", errors="replace")
                if score_table is not None:
                    finite_scores = pd.to_numeric(score_table["Score"], errors="coerce")
                    finite_scores = finite_scores[finite_scores.map(pd.notna)]
                    score_maximum = float(finite_scores.max()) if not finite_scores.empty else 1.0
                    score_maximum = score_maximum if score_maximum > 0 else 1.0
                    pdb_text = _pdb_with_residue_scores(pdb_text, score_table, maximum=score_maximum)
                visualizations = [StructureVisualization(
                    pdb=pdb_text,
                    color="uncertainty" if score_table is not None else "chain-id",  # type: ignore[arg-type]
                    color_params=_SCORE_COLOR_PARAMS if score_table is not None else None,
                    representation_type=representation,  # type: ignore[arg-type]
                    highlighted_selections=highlighted,
                )]
                if tool == "masif_seed":
                    point_cloud = _masif_point_cloud(run_dir)
                    if point_cloud:
                        visualizations.append(StructureVisualization(
                            pdb=point_cloud,
                            color="uncertainty",  # type: ignore[arg-type]
                            color_params=_SCORE_COLOR_PARAMS,
                            representation_type="ball-and-stick",
                        ))
                molstar_custom_component(
                    visualizations,
                    key=f"surface-result-view-{run_id}-{preview_structure.name}",
                    height=550,
                    show_controls=True,
                    show_structure_tools=False,
                )
        except OSError as exc:
            st.warning(f"Could not load the structure preview: {exc}")

    if surface_artifacts:
        st.subheader("Predicted surface files")
        st.caption("Download the native PLY surface output for inspection in a molecular-surface viewer.")
        columns = st.columns(min(3, len(surface_artifacts)))
        for index, (item, path) in enumerate(surface_artifacts):
            columns[index % len(columns)].download_button(
                f"Download {path.name}",
                data=path.read_bytes(),
                file_name=path.name,
                mime="application/octet-stream",
                key=f"surface-result-ply-{run_id}-{index}",
            )
    for index, (item, path) in enumerate(session_artifacts):
        st.download_button(
            f"Download PyMOL session · {path.name}",
            data=path.read_bytes(),
            file_name=path.name,
            mime="application/octet-stream",
            key=f"surface-result-pse-{run_id}-{index}",
        )

    if not csv_artifacts and score_table is None and preview_structure is None and not surface_artifacts and not session_artifacts:
        st.info("No prediction output is available yet. Refresh the job after it finishes to load its results.")


def render_surface_analysis_results(
    run_id: str,
    run_dir: Path,
    artifacts: list[tuple[dict, Path]],
    input_payload: dict,
) -> None:
    """Show cross-method residue evidence followed by each standalone renderer."""
    try:
        result_payload = json.loads((run_dir / "result.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        result_payload = {}
    outputs = result_payload.get("outputs") if isinstance(result_payload.get("outputs"), dict) else {}
    statuses = outputs.get("tool_status") if isinstance(outputs.get("tool_status"), dict) else {}
    if not statuses:
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
        statuses = read_json(run_dir / "artifacts/surface_analysis/tool_status.json")
    tools = ("pesto", "scannet", "surf2spot", "masif_seed")
    labels = {"pesto": "PeSTo", "scannet": "ScanNet", "surf2spot": "Surf2Spot", "masif_seed": "MaSIF-seed"}
    target = _target_structure("surface_analysis", run_dir, input_payload)

    normalized_tables: list[pd.DataFrame] = []
    by_tool: dict[str, list[tuple[dict, Path]]] = {}
    for tool in tools:
        root = (run_dir / "artifacts" / tool).resolve()
        tool_artifacts = [(item, path) for item, path in artifacts if path.is_relative_to(root)]
        by_tool[tool] = tool_artifacts
        if tool == "masif_seed":
            # The combined job may select only one chain for MaSIF while the
            # other methods use several. Project surface vertices back onto
            # MaSIF's filtered input to avoid assigning nearby partner-chain
            # residues the selected chain's surface scores.
            masif_target = _target_structure("masif_seed", run_dir, input_payload) or target
            table = _masif_score_table(run_dir, masif_target)
            if table is not None:
                normalized_tables.append(table.assign(Method=labels[tool]))
            continue
        for item, path in _artifact_paths(tool_artifacts, {".csv"}):
            try:
                raw = pd.read_csv(path, low_memory=False, nrows=50_000)
                table = _normalized_score_table(raw, tool, target, source_path=path)
            except (OSError, ValueError, pd.errors.ParserError, UnicodeError):
                table = None
            if table is not None:
                normalized_tables.append(table.assign(Method=labels[tool]))

    st.subheader("Combined surface analysis")
    st.caption("Compare PeSTo, ScanNet, Surf2Spot, and MaSIF-seed on the selected target. Individual method outputs remain below.")
    status_columns = st.columns(len(tools))
    for column, tool in zip(status_columns, tools):
        status = statuses.get(tool) if isinstance(statuses.get(tool), dict) else {}
        is_success = bool(status.get("success"))
        artifact_count = len(by_tool.get(tool, []))
        column.metric(labels[tool], "Complete" if is_success else "Partial / failed" if status else "Pending", f"{artifact_count} output file(s)")
    errors = [f"{labels[tool]}: {status.get('error')}" for tool in tools if isinstance((status := statuses.get(tool)), dict) and status.get("error")]
    if errors:
        st.warning("Some methods did not complete cleanly; any outputs they produced are still shown below. " + " ".join(errors))

    if normalized_tables:
        consensus = _combined_surface_consensus(normalized_tables)
        shared_count = int((consensus["Methods"] > 1).sum())
        top = str(consensus.iloc[0]["Residue ID"]) if not consensus.empty else "—"
        overview = st.columns(3)
        overview[0].metric("Residues compared", f"{len(consensus):,}")
        overview[1].metric("Residues scored by multiple methods", f"{shared_count:,}")
        overview[2].metric("Top combined rank", top)
        st.caption("Methods counts score availability, including low scores. It does not measure agreement on a hotspot. Residues are compared by coverage, then by their mean rank within each method. Raw scores from different methods and prediction modes are not calibrated against each other.")
        st.dataframe(consensus.head(100), hide_index=True, use_container_width=True)
        st.download_button(
            "Download combined residue table",
            consensus.to_csv(index=False).encode("utf-8"),
            file_name=f"surface-analysis-{run_id}-residue-consensus.csv",
            mime="text/csv",
            key=f"surface-analysis-consensus-download-{run_id}",
        )
    else:
        st.info("The methods completed without a shared residue-score table. Open each method below to inspect its native output files.")

    method_tabs = st.tabs([labels[tool] for tool in tools])
    raw_params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    for tab, tool in zip(method_tabs, tools):
        with tab:
            method_params = (
                {"interface": raw_params.get("pesto_interface", "protein")}
                if tool == "pesto" else
                {"mode": raw_params.get("scannet_mode", "interface"), "use_msa": raw_params.get("scannet_use_msa", False)}
                if tool == "scannet" else
                raw_params
            )
            method_payload = {**input_payload, "params": method_params}
            render_results(tool, f"{run_id}-{tool}", run_dir, by_tool[tool], method_payload)
