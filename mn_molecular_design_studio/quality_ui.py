"""Progress and compact artifact views for every QC batch."""
from pathlib import Path
import json
import streamlit as st


def _read(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return {}


def render_progress(run_dir, status, key):
    terminal = {"completed", "succeeded", "failed", "cancelled"}
    @st.fragment(run_every=5 if status not in terminal else None)
    def progress():
        data = _read(Path(run_dir) / "artifacts/analysis/batch_status.json")
        if data:
            total, finished = int(data.get("total", 0)), int(data.get("finished", 0))
            failed = sum(item.get("status") == "failed" for item in data.get("items", []))
            st.progress(min(1., finished / max(1, total)), text=f"{finished:,} / {total:,} evaluated · {failed:,} failed")
            if data.get("phase"):
                st.caption(data["phase"])
            if failed:
                with st.expander("Per-item failures"):
                    st.dataframe([{ "Candidate": item["id"], "Reason": item.get("error", "")} for item in data.get("items", []) if item.get("status") == "failed"], hide_index=True)
        elif status not in terminal:
            st.caption("Batch counters appear when the analysis container starts.")
        if status not in terminal and _read(Path(run_dir) / "metadata.json").get("status") in terminal:
            st.rerun(scope="app")
    progress()


def render_native_ligand_poses(run_dir, key, *, candidate_id=None):
    from mn_molecular_design_studio.md_results_ui import render_classified_interactions

    root = Path(run_dir)
    batch = _read(root / "artifacts/analysis/batch_manifest.json")
    paths = {}
    for item in batch.get("items", []):
        item_id = str(item.get("id") or "")
        if not item_id.startswith("item_") or not item_id[5:].isdigit():
            continue
        if candidate_id is not None and str(item.get("candidate_id")) != str(candidate_id):
            continue
        path = root / "artifacts/analysis" / item_id / "ligand_interactions.json"
        if path.is_file():
            paths[item_id] = (str(item.get("candidate_id") or item_id), path)
    if not paths:
        return
    item_id = next(iter(paths)) if len(paths) == 1 else st.selectbox(
        "Analyzed pose", list(paths), format_func=lambda value: paths[value][0], key=key + "-native-ligand-pose")
    _label, path = paths[item_id]
    report = _read(path)
    if report:
        pair_path = path.with_name("ligand_atom_pairs.csv")
        ligand_path = root / "artifacts/input" / item_id / "selected_ligand.pdb"
        try:
            pairs = pair_path.read_text(encoding="utf-8") if pair_path.is_file() else ""
            ligand = ligand_path.read_text(encoding="utf-8") if ligand_path.is_file() else ""
        except (OSError, UnicodeError):
            pairs, ligand = "", ""
        render_classified_interactions(st, report, pairs_csv=pairs, ligand_pdb=ligand,
                                       key_prefix=f"{key}-native-ligand-{item_id}")


def render_contact_maps(run_dir, key):
    import pandas as pd
    paths = sorted((Path(run_dir) / "artifacts/analysis").glob("item_*/contacts.json"))
    if not paths:
        return
    with st.expander("Interface residue maps"):
        manifest = _read(Path(run_dir) / "artifacts/analysis/batch_manifest.json")
        labels = {row["id"]: row["candidate_id"] for row in manifest.get("items", [])}
        path = st.selectbox("Complex", paths, format_func=lambda p: labels.get(p.parent.name, p.parent.name), key=key + "-contacts")
        data = _read(path)
        rows = data.get("residue_pairs", [])
        st.caption(f"{len(rows):,} residue pairs within {data.get('distance_angstrom')} Å; source chain and residue numbering is retained.")
        if rows:
            import plotly.express as px
            st.plotly_chart(px.scatter(pd.DataFrame(rows), x="target_residue", y="binder_residue"), key=key + "-contact-map", use_container_width=True)
        st.download_button("Download contact map JSON", path.read_bytes(), file_name=labels.get(path.parent.name, path.parent.name) + "_contacts.json", key=key + "-contacts-download")
