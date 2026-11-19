"""Reader for BindCraft 1 (v1.5-style) run folders.

Expected folder (case-insensitive folder names; users prune, so every part may be missing but the CSVs):

    <run>/trajectory_stats.csv  mpnn_design_stats.csv  final_design_stats.csv  failure_csv.csv
    <run>/accepted/*.pdb  accepted/Ranked/<rank>_<design>_model<k>.pdb  (+ plots, animations)
    <run>/rejected/*.pdb
    <run>/<settings>.json  (optional; its basename is the TargetSettings column)
    <run>/*.xlsx           (optional copy of final_design_stats; may carry cell highlights)

Joins (verified on real runs, docs/data_inventory.md):
  accepted: final_design_stats.Design (no _model suffix) <-> accepted/<Design>_model<k>.pdb
  rejected: mpnn_design_stats rows that are not in final_design_stats <-> rejected/<Design>_model<k>.pdb
  backbone: Design without _mpnn<n> <-> trajectory_stats.Design
Structures are only referenced and lightly checked here (chains, residue numbers, sequence, B-factor
scale). Parsing for analysis lives in structure.py.
"""

from __future__ import annotations

import ast
import hashlib
import math
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

import pandas as pd

from . import ReaderError
from ..config import (
    ACCEPTED_DIR, CSV_FAILURE, CSV_FINAL, CSV_MPNN, CSV_TRAJECTORY, DEFAULT_FILTERS, MODEL_NUMBERS,
    REJECTED_DIR, TRAJECTORY_COUNTERS, check_filters,
)
from ..model import (
    METRIC_STEMS, SINGLE_VALUE_STEMS, BackboneRecord, DesignRecord, FunnelCounts, InputFile, RunData,
    RunSettings, TargetConstruct,
)

NAME = "bindcraft1"

NAME_RE = re.compile(r"^(?P<binder>.+)_l(?P<length>\d+)_s(?P<seed>\d+)(?:_mpnn(?P<mpnn>\d+))?(?:_model(?P<model>\d+))?$")
_TIME_RE = re.compile(r"(?:(\d+)\s*hours?)?[, ]*(?:(\d+)\s*minutes?)?[, ]*(?:(\d+)\s*seconds?)?")
_ID_IN_TEXT = re.compile(r"\b[A-Za-z0-9]+_l\d+_s\d+(?:_mpnn\d+)?\b")

# Columns BindCraft writes that are not metrics.
_STR_COLS = {"Design": str, "Protocol": str, "Target_Hotspot": str, "Sequence": str, "InterfaceResidues": str,
             "Notes": str, "TargetSettings": str, "Filters": str, "AdvancedSettings": str,
             "TrajectoryTime": str, "DesignTime": str}


# ---------------------------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------------------------


def parse_name(stem: str) -> Optional[dict[str, Any]]:
    """Split a BindCraft 1 design or file stem. Returns None if it does not follow the scheme.

    {"binder": "binder01", "length": 77, "seed": 305024, "mpnn": 4 or None, "model": 2 or None,
     "design_id": "binder01_l77_s305024_mpnn4", "backbone_id": "binder01_l77_s305024"}
    """
    m = NAME_RE.match(stem)
    if not m:
        return None
    binder, length, seed = m["binder"], int(m["length"]), int(m["seed"])
    mpnn = int(m["mpnn"]) if m["mpnn"] is not None else None
    model = int(m["model"]) if m["model"] is not None else None
    backbone = f"{binder}_l{length}_s{seed}"
    design = backbone + (f"_mpnn{mpnn}" if mpnn is not None else "")
    return {"binder": binder, "length": length, "seed": seed, "mpnn": mpnn, "model": model,
            "design_id": design, "backbone_id": backbone}


def _find_child(root: Path, name: str) -> Optional[Path]:
    """Case-insensitive child lookup (macOS is case-insensitive, Linux is not)."""
    if not root.is_dir():
        return None
    exact = root / name
    if exact.exists():
        return exact
    low = name.lower()
    for p in root.iterdir():
        if p.name.lower() == low:
            return p
    return None


def detect(path: str | Path) -> bool:
    """True when the folder looks like a BindCraft 1 run (the three stats CSVs are present)."""
    root = Path(path)
    if not root.is_dir():
        return False
    return all(_find_child(root, n) is not None for n in (CSV_FINAL, CSV_MPNN, CSV_TRAJECTORY))


def _f(x: Any) -> Optional[float]:
    if x is None:
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def _s(x: Any) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return str(x)


def parse_duration(text: Any) -> Optional[float]:
    """'0 hours, 9 minutes, 2 seconds' -> 542.0"""
    t = _s(text)
    if not t:
        return None
    m = _TIME_RE.fullmatch(t.strip())
    if not m or not any(m.groups()):
        return None
    h, mi, s = (int(g) if g else 0 for g in m.groups())
    return float(h * 3600 + mi * 60 + s)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=_STR_COLS)


def _records(df: pd.DataFrame) -> list[dict[str, Any]]:
    return df.to_dict("records")


def _parse_aas(value: Any) -> dict[str, float]:
    t = _s(value)
    if not t:
        return {}
    try:
        d = ast.literal_eval(t)
    except (ValueError, SyntaxError):
        return {}
    return {str(k): float(v) for k, v in d.items()} if isinstance(d, dict) else {}


def _metric_columns(row: Mapping[str, Any]) -> tuple[dict[str, Optional[float]], dict[int, dict[str, Optional[float]]]]:
    avg: dict[str, Optional[float]] = {}
    models: dict[int, dict[str, Optional[float]]] = {}
    for stem in METRIC_STEMS:
        if stem in SINGLE_VALUE_STEMS:
            if stem in row:
                avg[stem] = _f(row[stem])
            continue
        key = f"Average_{stem}"
        if key in row:
            avg[stem] = _f(row[key])
        for k in MODEL_NUMBERS:
            col = f"{k}_{stem}"
            if col in row:
                v = _f(row[col])
                if v is not None:
                    models.setdefault(k, {})[stem] = v
    return avg, models


def _parse_interface_residues(value: Any) -> list[str]:
    t = _s(value)
    return [x.strip() for x in t.split(",") if x.strip()] if t else []


def parse_hotspots(text: Any) -> list[int]:
    """'26,122,125' -> [26, 122, 125]; also accepts chain letters ('A26') and ranges ('30-33')."""
    out: list[int] = []
    for tok in _s(text).replace(";", ",").split(","):
        tok = re.sub(r"^[A-Za-z]+", "", tok.strip())
        if not tok:
            continue
        if "-" in tok:
            a, _, b = tok.partition("-")
            if a.isdigit() and b.isdigit():
                out.extend(range(int(a), int(b) + 1))
        elif tok.isdigit():
            out.append(int(tok))
    return out


def _aggregate(items: list[str], what: str, limit: int = 4) -> str:
    shown = ", ".join(items[:limit]) + (f", ... ({len(items)} in all)" if len(items) > limit else "")
    return f"{what}: {shown}"


# ---------------------------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------------------------


def _find_settings_json(root: Path, target_settings_name: Optional[str]) -> Optional[Path]:
    import json

    candidates: list[Path] = []
    if target_settings_name:
        c = _find_child(root, f"{target_settings_name}.json")
        if c:
            candidates.append(c)
    candidates += sorted(p for p in root.glob("*.json") if p not in candidates)
    for p in candidates:
        try:
            data = json.loads(p.read_text())
        except (ValueError, OSError):
            continue
        if isinstance(data, dict) and "binder_name" in data and ("starting_pdb" in data or "lengths" in data):
            return p
    return None


def _read_settings(
    run_id: str, label: str, root: Path, settings_path: Optional[Path], traj: pd.DataFrame, mpnn: pd.DataFrame,
) -> tuple[RunSettings, Optional[Path]]:
    import json

    st = RunSettings(run_id=run_id, label=label)
    first = _records(mpnn.head(1))[0] if len(mpnn) else (_records(traj.head(1))[0] if len(traj) else {})
    st.protocol = _s(first.get("Protocol")) or None
    st.advanced_settings_name = _s(first.get("AdvancedSettings")) or None
    st.filters_name = _s(first.get("Filters")) or None
    st.target_settings_name = _s(first.get("TargetSettings")) or None
    st.helicity = _f(first.get("Helicity"))
    st.hotspot_string = _s(first.get("Target_Hotspot")) or None
    if st.hotspot_string:
        st.hotspots = parse_hotspots(st.hotspot_string)
        st.hotspot_source = "csv"

    path = settings_path or _find_settings_json(root, st.target_settings_name)
    if path is not None and Path(path).is_file():
        data = json.loads(Path(path).read_text())
        st.settings_path = str(Path(path).resolve())
        st.settings_raw = data
        st.binder_name = data.get("binder_name")
        st.design_path_recorded = data.get("design_path")
        sp = data.get("starting_pdb")
        st.starting_pdb_recorded = sp
        st.starting_pdb_name = Path(sp).name if sp else None
        st.chains = data.get("chains")
        hs = data.get("target_hotspot_residues")
        if hs not in (None, ""):
            st.hotspots = parse_hotspots(hs)
            st.hotspot_string = str(hs)
            st.hotspot_source = "settings_json"
            if first and _s(first.get("Target_Hotspot")) and parse_hotspots(first["Target_Hotspot"]) != st.hotspots:
                st.notes.append("Target_Hotspot column differs from the settings file; the settings file is used.")
        lengths = data.get("lengths")
        if isinstance(lengths, (list, tuple)) and len(lengths) == 2:
            st.length_min, st.length_max = int(lengths[0]), int(lengths[1])
            st.length_source = "settings_json"
        if data.get("number_of_final_designs") is not None:
            st.final_designs_target = int(data["number_of_final_designs"])
            st.final_designs_target_source = "settings_json"
    else:
        path = None
        st.notes.append(
            f"No settings file in the run folder (TargetSettings = {st.target_settings_name!r}): hotspots come "
            "from the Target_Hotspot column, lengths from the observed Length range."
        )
        info = parse_name(_s(first.get("Design"))) if first else None
        if info:
            st.binder_name = info["binder"]

    if st.length_min is None and len(traj):
        st.length_min, st.length_max = int(traj["Length"].min()), int(traj["Length"].max())
        st.length_source = "inferred: min and max of trajectory_stats Length"
    return st, (Path(path) if path else None)


# ---------------------------------------------------------------------------------------------
# Structures
# ---------------------------------------------------------------------------------------------


def _one_letter(resname: str) -> str:
    import gemmi

    info = gemmi.find_tabulated_residue(resname)
    return info.one_letter_code.upper() if info and info.one_letter_code else "X"


def inspect_structure(path: Path) -> dict[str, Any]:
    """Light check of one BindCraft 1 PDB: chains, target residue numbers, sequences, B-factor scale."""
    import gemmi

    st = gemmi.read_structure(str(path))
    model = st[0]
    chains = [c.name for c in model]
    out: dict[str, Any] = {"chains": chains}
    if "A" in chains:
        a = model["A"]
        out["target_numbers"] = tuple(r.seqid.num for r in a)
        out["target_seq"] = "".join(_one_letter(r.name) for r in a)
        out["hydrogens"] = any(at.element.name == "H" for r in a for at in r)
    if "B" in chains:
        b = model["B"]
        out["binder_seq"] = "".join(_one_letter(r.name) for r in b)
        bs = [at.b_iso for r in b for at in r]
        out["binder_mean_b"] = sum(bs) / len(bs) if bs else None
    return out


def _structure_files(directory: Optional[Path]) -> dict[str, tuple[Path, Optional[int]]]:
    """design id -> (path, model k) for the *.pdb files directly inside a folder."""
    out: dict[str, tuple[Path, Optional[int]]] = {}
    if directory is None or not directory.is_dir():
        return out
    for p in sorted(directory.glob("*.pdb")):
        info = parse_name(p.stem)
        key = info["design_id"] if info else p.stem
        out[key] = (p.resolve(), info["model"] if info else None)
    return out


def _ranked_files(ranked_dir: Optional[Path]) -> dict[int, Path]:
    out: dict[int, Path] = {}
    if ranked_dir is None or not ranked_dir.is_dir():
        return out
    for p in sorted(ranked_dir.glob("*.pdb")):
        m = re.match(r"^(\d+)_", p.name)
        if m:
            out[int(m.group(1))] = p.resolve()
    return out


def _gaps(numbers: list[int]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for a, b in zip(numbers, numbers[1:]):
        if b > a + 1:
            out.append((a + 1, b - 1))
    return out


def _missing_ranges(missing: list[int]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for n in sorted(missing):
        if out and n == out[-1][1] + 1:
            out[-1] = (out[-1][0], n)
        else:
            out.append((n, n))
    return out


def _reference_construct(path: Path) -> Optional[tuple[list[int], str]]:
    try:
        info = inspect_structure(path)
    except Exception:  # unreadable reference: the caller notes it
        return None
    if "target_numbers" not in info:
        return None
    return list(info["target_numbers"]), info["target_seq"]


# ---------------------------------------------------------------------------------------------
# Annotations found next to the data
# ---------------------------------------------------------------------------------------------


def _rtf_to_text(raw: str) -> str:
    text = re.sub(r"\\'[0-9a-fA-F]{2}", "", raw)
    text = re.sub(r"\\[a-zA-Z]+-?\d* ?", " ", text)
    text = re.sub(r"[{}]", " ", text)
    return re.sub(r"[ \t]+", " ", text).strip()


def _read_note_files(accepted_dir: Optional[Path], root: Path) -> list[dict[str, Any]]:
    notes: list[dict[str, Any]] = []
    for base in [d for d in (accepted_dir, root) if d is not None]:
        for p in sorted(base.glob("*")):
            if p.suffix.lower() not in (".rtf", ".txt", ".md") or not p.is_file():
                continue
            raw = p.read_text(errors="replace")
            text = _rtf_to_text(raw) if p.suffix.lower() == ".rtf" else raw
            notes.append({
                "path": str(p.resolve()), "text": text[:600],
                "ids": sorted(set(_ID_IN_TEXT.findall(text))),
            })
    return notes


def read_xlsx_highlights(path: Path) -> Optional[dict[str, Any]]:
    """Designs whose row carries a solid cell fill in a final_design_stats-style workbook.

    Returns {"file", "sheet", "colours": {rgb: n_rows}, "designs": [ids]} or None. Rows where every
    filled cell has the same single colour count as highlighted. The meaning of the highlight is not
    recorded in the file (in some runs it matches the designs that were sent to the lab, but that is an inference
    from names).
    """
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=False, data_only=True)
    for ws in wb.worksheets:
        header = [c.value for c in ws[1]]
        if "Design" not in header:
            continue
        di = header.index("Design")
        colours: dict[str, int] = {}
        designs: list[str] = []
        for row in ws.iter_rows(min_row=2):
            name = row[di].value
            fills = {c.fill.fgColor.rgb for c in row if c.fill is not None and c.fill.fill_type == "solid"}
            fills = {f for f in fills if isinstance(f, str) and f not in ("FFFFFFFF", "00000000")}
            if name and fills:
                designs.append(str(name))
                for f in fills:
                    colours[f] = colours.get(f, 0) + 1
        return {"file": str(path.resolve()), "sheet": ws.title, "colours": colours, "designs": designs}
    return None


# ---------------------------------------------------------------------------------------------
# read_run
# ---------------------------------------------------------------------------------------------


def read_run(
    path: str | Path,
    run_id: str,
    label: Optional[str] = None,
    settings_path: Optional[str | Path] = None,
    construct_path: Optional[str | Path] = None,
    check_structures: bool = True,
    hash_structures: bool = True,
) -> RunData:
    """Read one BindCraft 1 run folder into a RunData.

    construct_path     the PDB that was given to BindCraft, when it is known (the numbers present in the
                       structures are compared with it and deleted regions are reported against it).
    check_structures   open every structure once with gemmi (target numbers, chains, sequences).
    hash_structures    include the structure files in RunData.files with their sha256.
    """
    root = Path(path).resolve()
    if not detect(root):
        raise ReaderError(f"{root} does not look like a BindCraft 1 run folder (missing the stats CSVs)")
    label = label or root.name
    warnings: list[str] = []
    notes: list[str] = []
    files: list[InputFile] = []

    def add_file(role: str, p: Path) -> None:
        files.append(InputFile(role=role, path=str(p.resolve()), sha256=sha256_file(p), size=p.stat().st_size))

    csv_final = _find_child(root, CSV_FINAL)
    csv_mpnn = _find_child(root, CSV_MPNN)
    csv_traj = _find_child(root, CSV_TRAJECTORY)
    csv_fail = _find_child(root, CSV_FAILURE)
    assert csv_final and csv_mpnn and csv_traj
    final_df, mpnn_df, traj_df = _read_csv(csv_final), _read_csv(csv_mpnn), _read_csv(csv_traj)
    for p in (csv_final, csv_mpnn, csv_traj, csv_fail):
        if p is not None:
            add_file("csv", p)

    failure_counters: dict[str, int] = {}
    if csv_fail is not None:
        fdf = pd.read_csv(csv_fail)
        if len(fdf):
            failure_counters = {str(k): int(v) for k, v in fdf.iloc[0].items() if _f(v) is not None}
    else:
        notes.append("failure_csv.csv is missing: failure counters are unknown.")

    settings, settings_file = _read_settings(
        run_id, label, root, Path(settings_path) if settings_path else None, traj_df, mpnn_df)
    if settings_file is not None:
        add_file("settings", settings_file)
    notes.extend(settings.notes)
    ver = settings.design_path_recorded or ""
    mver = re.search(r"BindCraft[_-]?v?(\d+)[_.](\d+)", ver)
    tool_version = f"v{mver.group(1)}.{mver.group(2)} (from design_path in the settings file)" if mver else None

    acc_dir = _find_child(root, ACCEPTED_DIR)
    rej_dir = _find_child(root, REJECTED_DIR)
    ranked_dir = _find_child(acc_dir, "Ranked") if acc_dir else None
    acc_files = _structure_files(acc_dir)
    rej_files = _structure_files(rej_dir)
    ranked = _ranked_files(ranked_dir)

    # -- designs -----------------------------------------------------------------------------
    final_rows = {_s(r["Design"]): r for r in _records(final_df)}
    mpnn_rows = _records(mpnn_df)
    mpnn_ids = {_s(r["Design"]) for r in mpnn_rows}
    designs: list[DesignRecord] = []
    unparsed: list[str] = []

    def build(row: Mapping[str, Any], status: str) -> Optional[DesignRecord]:
        did = _s(row["Design"])
        info = parse_name(did)
        if info is None or info["mpnn"] is None:
            unparsed.append(did)
            return None
        avg, models = _metric_columns(row)
        aas = {"avg": _parse_aas(row.get("Average_InterfaceAAs"))}
        for k in MODEL_NUMBERS:
            if f"{k}_InterfaceAAs" in row:
                d = _parse_aas(row[f"{k}_InterfaceAAs"])
                if d:
                    aas[str(k)] = d
        rec = DesignRecord(
            id=info["design_id"], run_id=run_id, backbone_id=info["backbone_id"], status=status,
            binder_name=info["binder"], length=int(row.get("Length", info["length"])), seed=int(row.get("Seed", info["seed"])),
            mpnn_index=info["mpnn"], sequence=_s(row.get("Sequence")), protocol=_s(row.get("Protocol")) or None,
            target_hotspot=_s(row.get("Target_Hotspot")) or None, metrics=avg, models=models,
            interface_residues=_parse_interface_residues(row.get("InterfaceResidues")), interface_aas=aas,
            design_time_s=parse_duration(row.get("DesignTime")), notes=_s(row.get("Notes")),
        )
        if status == "accepted" and _f(row.get("Rank")) is not None:
            rec.rank = int(row["Rank"])
        rec.failed_filters = check_filters(avg, models, aas, DEFAULT_FILTERS)
        return rec

    for did, row in final_rows.items():
        rec = build(row, "accepted")
        if rec is not None:
            designs.append(rec)
    accepted_ids = {d.id for d in designs}
    if set(final_rows) - mpnn_ids:
        missing = sorted(set(final_rows) - mpnn_ids)
        warnings.append(_aggregate(missing, "accepted designs without a row in mpnn_design_stats.csv"))
    for row in mpnn_rows:
        if _s(row["Design"]) in final_rows:
            continue
        rec = build(row, "rejected")
        if rec is not None:
            designs.append(rec)
    if unparsed:
        warnings.append(_aggregate(unparsed, "design names that do not follow <name>_l<len>_s<seed>_mpnn<n>, skipped"))
    if len({d.id for d in designs}) != len(designs):
        warnings.append("duplicate design ids across the CSVs; the first row of each was kept")
        seen: set[str] = set()
        designs = [d for d in designs if not (d.id in seen or seen.add(d.id))]

    # structures
    no_struct_acc, no_struct_rej = [], []
    for d in designs:
        pool = acc_files if d.status == "accepted" else rej_files
        hit = pool.get(d.id)
        if hit is not None:
            d.structure_path, d.structure_model = str(hit[0]), hit[1]
        elif d.status == "accepted":
            no_struct_acc.append(d.id)
        elif rej_dir is not None:
            no_struct_rej.append(d.id)
        if d.status == "accepted" and d.rank is not None and d.rank in ranked:
            d.ranked_path = str(ranked[d.rank])
    if no_struct_acc:
        warnings.append(_aggregate(no_struct_acc, f"{len(no_struct_acc)} accepted designs have a CSV row but no PDB in accepted/"))
    if no_struct_rej:
        warnings.append(_aggregate(no_struct_rej, f"{len(no_struct_rej)} not-accepted designs have a CSV row but no PDB in rejected/"))
    known = {d.id for d in designs}
    orphan_acc = sorted(set(acc_files) - known)
    orphan_rej = sorted(set(rej_files) - known)
    if orphan_acc:
        warnings.append(_aggregate(orphan_acc, f"{len(orphan_acc)} PDBs in accepted/ have no CSV row"))
    if orphan_rej:
        warnings.append(_aggregate(orphan_rej, f"{len(orphan_rej)} PDBs in rejected/ have no CSV row"))
    if rej_dir is None:
        n_rej = sum(1 for d in designs if d.status == "rejected")
        notes.append(f"No rejected/ folder: {n_rej} not-accepted sequences are metric-only (no structure, no near-miss analysis).")
    stray_ranked = [r for r in ranked if r not in {d.rank for d in designs if d.status == "accepted"}]
    if stray_ranked:
        warnings.append(f"Ranked/ holds rank numbers with no Rank row: {sorted(stray_ranked)[:5]}")

    # -- backbones ---------------------------------------------------------------------------
    traj_rows = _records(traj_df)
    backbones: dict[str, BackboneRecord] = {}
    for i, row in enumerate(traj_rows):
        did = _s(row["Design"])
        info = parse_name(did)
        if info is None:
            continue
        tm: dict[str, Optional[float]] = {s: _f(row[s]) for s in METRIC_STEMS if s in row and s not in SINGLE_VALUE_STEMS}
        backbones[did] = BackboneRecord(
            id=did, run_id=run_id, binder_name=info["binder"], length=int(row.get("Length", info["length"])),
            seed=int(row.get("Seed", info["seed"])), order=i, trajectory_metrics=tm,
            trajectory_time_s=parse_duration(row.get("TrajectoryTime")), start_sequence=_s(row.get("Sequence")) or None,
            interface_residues=_parse_interface_residues(row.get("InterfaceResidues")), notes=_s(row.get("Notes")),
        )
    if len(backbones) != len(traj_rows):
        warnings.append("trajectory_stats.csv has duplicate or unparseable Design names")
    orphan_bb: list[str] = []
    for d in designs:
        bb = backbones.get(d.backbone_id)
        if bb is None:
            orphan_bb.append(d.backbone_id)
            bb = backbones[d.backbone_id] = BackboneRecord(
                id=d.backbone_id, run_id=run_id, binder_name=d.binder_name, length=d.length, seed=d.seed, order=-1)
        bb.design_ids.append(d.id)
        if d.status == "accepted":
            bb.accepted_ids.append(d.id)
    if orphan_bb:
        warnings.append(_aggregate(sorted(set(orphan_bb)), "backbones with designs but no row in trajectory_stats.csv"))
    by_mpnn = {d.id: d.mpnn_index or 0 for d in designs}
    for bb in backbones.values():
        bb.design_ids.sort(key=lambda i: by_mpnn[i])
        bb.accepted_ids.sort(key=lambda i: by_mpnn[i])
    backbone_list = sorted(backbones.values(), key=lambda b: (b.order < 0, b.order, b.id))

    # -- structure check, target construct ---------------------------------------------------
    target = TargetConstruct(run_id=run_id)
    if check_structures:
        _check_structures(root, designs, settings, target, warnings, construct_path)
    for n in settings.hotspots:
        letter = target.residue_at(n)
        if letter is not None:
            target.hotspot_residues[n] = letter
    if target.numbers and settings.hotspots:
        absent = [n for n in settings.hotspots if n not in target.numbers]
        if absent:
            warnings.append(f"hotspot residues {absent} are not present in the target of this run")

    # -- consistency checks on the data ------------------------------------------------------
    hs = {d.target_hotspot for d in designs if d.target_hotspot}
    if len(hs) > 1:
        warnings.append(f"Target_Hotspot differs between rows: {sorted(hs)[:4]}; the run is not homogeneous")
    acc = [d for d in designs if d.status == "accepted"]
    if acc and settings.filters_name and settings.filters_name != "default_filters":
        notes.append(f"The run used the filter file {settings.filters_name!r}, not default_filters; DEFAULT_FILTERS may not apply.")
    bad_acc = [d.id for d in acc if d.failed_filters]
    if bad_acc:
        warnings.append(_aggregate(bad_acc, f"{len(bad_acc)} accepted designs fail DEFAULT_FILTERS (the run probably used other thresholds)"))
    unexplained = [d for d in designs if d.status == "rejected" and not d.failed_filters]
    if unexplained:
        notes.append(f"{len(unexplained)} not-accepted sequences pass every inferred filter (rejected by a rule DEFAULT_FILTERS does not cover).")

    # -- inputs, hashes, annotations ---------------------------------------------------------
    for x in sorted(root.glob("*.xlsx")):
        add_file("xlsx", x)
    if hash_structures:
        for d in designs:
            if d.structure_path:
                add_file(f"{d.status}_structure", Path(d.structure_path))
    annotations: dict[str, Any] = {}
    for x in sorted(root.glob("*.xlsx")):
        if x.name.startswith("~$"):
            continue
        try:
            hl = read_xlsx_highlights(x)
        except Exception as exc:  # corrupt or locked workbook
            notes.append(f"xlsx {x.name} could not be read: {exc}")
            continue
        if hl and hl["designs"]:
            annotations.setdefault("xlsx_highlight", []).append(hl)
            marked = set(hl["designs"])
            for d in designs:
                if d.id in marked:
                    d.xlsx_marked = True
    note_files = _read_note_files(acc_dir, root)
    if note_files:
        annotations["note_files"] = note_files
        for nf in note_files:
            add_file("note", Path(nf["path"]))

    funnel = _funnel(designs, backbone_list, failure_counters, settings, mpnn_rows)
    return RunData(
        run_id=run_id, label=label, reader=NAME, source_dir=str(root), settings=settings, target=target,
        designs=designs, backbones=backbone_list, funnel=funnel, failure_counters=failure_counters, files=files,
        warnings=warnings, notes=notes, annotations=annotations, tool_version=tool_version,
    )


def _check_structures(
    root: Path, designs: list[DesignRecord], settings: RunSettings, target: TargetConstruct,
    warnings: list[str], construct_path: Optional[str | Path],
) -> None:
    from collections import Counter

    sets: Counter[tuple[int, ...]] = Counter()
    seqs: dict[tuple[int, ...], str] = {}
    bad_chain: list[str] = []
    bad_len: list[str] = []
    bad_seq: list[str] = []
    bad_b: list[str] = []
    set_of: dict[str, tuple[int, ...]] = {}
    n = 0
    for d in designs:
        if not d.structure_path:
            continue
        try:
            info = inspect_structure(Path(d.structure_path))
        except Exception as exc:
            warnings.append(f"{d.id}: structure unreadable ({exc})")
            continue
        n += 1
        if info["chains"] != ["A", "B"]:
            bad_chain.append(d.id)
            continue
        nums = tuple(info["target_numbers"])
        sets[nums] += 1
        seqs.setdefault(nums, info["target_seq"])
        set_of[d.id] = nums
        if len(info["binder_seq"]) != d.length:
            bad_len.append(d.id)
        if d.sequence and info["binder_seq"] != d.sequence:
            bad_seq.append(d.id)
        mb = info.get("binder_mean_b")
        want = d.pose_metric("pLDDT")
        if mb is not None and want is not None and abs(mb / 100.0 - want) > 0.02:
            bad_b.append(d.id)
    target.n_structures_checked = n
    if bad_chain:
        warnings.append(_aggregate(bad_chain, f"{len(bad_chain)} structures do not have exactly chains A and B"))
    if bad_len:
        warnings.append(_aggregate(bad_len, f"{len(bad_len)} structures whose binder chain length differs from the CSV Length"))
    if bad_seq:
        warnings.append(_aggregate(bad_seq, f"{len(bad_seq)} structures whose binder sequence differs from the CSV Sequence"))
    if bad_b:
        warnings.append(_aggregate(bad_b, f"{len(bad_b)} structures whose binder mean B-factor / 100 differs from the pLDDT of the saved model by more than 0.02"))
    if not sets:
        return
    majority, _ = sets.most_common(1)[0]
    target.numbers = list(majority)
    target.sequence = seqs[majority]
    target.first, target.last = majority[0], majority[-1]
    target.consistent = len(sets) == 1
    if not target.consistent:
        odd = [i for i, s in set_of.items() if s != majority]
        warnings.append(_aggregate(odd, f"{len(odd)} structures have a different target residue set than the majority ({len(majority)} residues)"))

    ref_path = Path(construct_path) if construct_path else None
    if ref_path is None and settings.starting_pdb_name:
        for base in (root, root.parent):
            c = _find_child(base, settings.starting_pdb_name)
            if c is not None and c.is_file():
                ref_path = c
                break
    ref = _reference_construct(ref_path) if ref_path and ref_path.is_file() else None
    if ref is not None:
        ref_numbers, ref_seq = ref
        target.reference_path = str(ref_path.resolve())  # type: ignore[union-attr]
        target.reference_length = len(ref_numbers)
        by_num = dict(zip(ref_numbers, ref_seq))
        mism = [nn for nn, aa in zip(target.numbers, target.sequence) if by_num.get(nn) != aa]
        target.sequence_matches_reference = not mism
        if mism:
            warnings.append(f"{len(mism)} target residues differ from the construct file {Path(ref_path).name} (first at {mism[0]})")  # type: ignore[arg-type]
        target.deleted_ranges = _missing_ranges([x for x in ref_numbers if x not in set(target.numbers)])
    else:
        target.deleted_ranges = _gaps(target.numbers)


def _funnel(
    designs: list[DesignRecord], backbones: list[BackboneRecord], counters: Mapping[str, int],
    settings: RunSettings, mpnn_rows: list[dict[str, Any]],
) -> FunnelCounts:
    f = FunnelCounts()
    st, nt = f.status, f.notes
    traj_rows = sum(1 for b in backbones if b.order >= 0)
    f.trajectory_rows = traj_rows
    st["trajectory_rows"] = "observed"
    nt["trajectory_rows"] = "Rows in trajectory_stats.csv: trajectories that finished and were scored."
    f.mpnn_rows = len(designs)
    st["mpnn_rows"] = "observed"
    nt["mpnn_rows"] = ("Sequences with a row in mpnn_design_stats.csv. Sequences that failed the AlphaFold2 pre-filters "
                       "(pLDDT, i_pTM, i_pAE) get no row and are only counted in failure_csv.csv.")
    f.accepted_designs = sum(1 for d in designs if d.status == "accepted")
    f.mpnn_rows_rejected = len(designs) - f.accepted_designs
    f.backbones_with_mpnn_rows = sum(1 for b in backbones if b.design_ids)
    f.accepted_backbones = sum(1 for b in backbones if b.accepted_ids)
    for k in ("mpnn_rows_rejected", "backbones_with_mpnn_rows", "accepted_designs", "accepted_backbones"):
        st[k] = "observed"
    nt["accepted_backbones"] = "A family is one backbone (trajectory seed); counts of families, never of sequences."

    f.trajectories_attempted = None
    st["trajectories_attempted"] = "unknown"
    nt["trajectories_attempted"] = ("Unknown: the Trajectory/ folders that would show every attempt, including the ones "
                                    "terminated early, are not in the run folder.")
    present = [c for c in TRAJECTORY_COUNTERS if c in counters]
    if present:
        f.trajectories_terminated = sum(counters[c] for c in present)
        f.trajectories_attempted_estimate = traj_rows + f.trajectories_terminated
        st["trajectories_terminated"] = st["trajectories_attempted_estimate"] = "derived"
        nt["trajectories_terminated"] = "Sum of the Trajectory_* failure counters (" + ", ".join(f"{c}={counters[c]}" for c in present) + ")."
        nt["trajectories_attempted_estimate"] = (
            "Derived: trajectory rows plus the Trajectory_* counters, assuming a terminated trajectory is counted "
            "once (BindCraft stops a trajectory at its first failed trajectory check).")
    idx = [int(parse_name(_s(r["Design"]))["mpnn"]) for r in mpnn_rows if parse_name(_s(r["Design"])) and parse_name(_s(r["Design"]))["mpnn"] is not None]  # type: ignore[index]
    if idx:
        f.mpnn_sequences_per_trajectory = max(idx)
        st["mpnn_sequences_per_trajectory"] = "observed"
        f.mpnn_sequences_scored_estimate = max(idx) * traj_rows
        st["mpnn_sequences_scored_estimate"] = "derived"
        nt["mpnn_sequences_scored_estimate"] = (
            "Upper estimate: the highest mpnn index seen times the trajectory rows, assuming every written "
            "trajectory entered sequence design. Consistent with the i_pAE counter plus the kept rows.")
    f.final_designs_target = settings.final_designs_target
    if f.final_designs_target is None:
        _infer_final_target(backbones, f, settings)
    else:
        st["final_designs_target"] = "observed"
    return f


def _infer_final_target(backbones: list[BackboneRecord], f: FunnelCounts, settings: RunSettings) -> None:
    """BindCraft stops after the trajectory that takes the accepted count to number_of_final_designs."""
    ordered = [b for b in backbones if b.order >= 0]
    if not ordered or not f.accepted_designs:
        return
    last_row = max(ordered, key=lambda b: b.order)
    if not last_row.accepted_ids:
        return
    before = f.accepted_designs - len(last_row.accepted_ids)
    lo, hi = before + 1, f.accepted_designs
    if lo == hi:
        f.final_designs_target = hi
        f.status["final_designs_target"] = "derived"
        f.notes["final_designs_target"] = "Inferred: the last trajectory took the accepted count from %d to %d, which ended the run." % (before, hi)
        settings.final_designs_target = hi
        settings.final_designs_target_source = "inferred: the run ended on the trajectory that reached this count"
    else:
        f.notes["final_designs_target"] = f"Not recoverable: between {lo} and {hi}."
