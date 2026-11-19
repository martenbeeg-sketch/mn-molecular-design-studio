"""Project -> the blobs of docs/data_contract.md: campaign (RAW), inspector (INS), homolog (HOM), profile and manifest.

The exporter is the only module that writes blobs. It converts the construct numbering of ``RunData`` and of the interface
analysis to UniProt numbering once (``TargetConstruct.to_uniprot``), copies the scores BindCraft wrote, recomputes the rest from
the structures (the interface, liabilities and context modules), assembles the design records, the per-run blocks, the target
block and the manifest, and validates every output against docs/schema (``scripts/validate_schema.py``). A violation raises
``ExportError``: the build fails rather than writing a bundle the pages cannot read. Binder sequences are written into the
campaign and inspector blobs only (the contract requires ``seq``), never into the manifest, a doc or a log.

Public entry points: ``build_blobs(project)`` -> ``Blobs`` and ``write_export(blobs, out_dir)``. docs/export.md has the layout,
the numbering rules, the metric map and the validation steps.
"""
from __future__ import annotations

import datetime as _dt
import importlib.util
import json
import math
import platform
import re
import statistics
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import numpy as np

from . import __version__, config, families, liabilities, model, paratope, termini
from . import refdata as rd
from . import structure as S
from . import target_frame as TF
from .model import DesignRecord, RunData
from .pipeline import Project, ProjectKey, sha256_file, sha256_json, site_residues, utc_now

EXPORT_VERSION = 1
GENERATOR = "binderscope export.py"
INSPECTOR_MODALITY = "minibinder"

#: The 38 capability flags (profile_spec.md section 4). The exporter writes every one of them.
CAP_FLAGS = (
    "outcomes", "weights", "olig", "rest_ctx", "flag_region", "footprints", "has_pae", "has_traj", "traj_total", "has_fold",
    "families", "has_boltz", "models", "has_seeds", "ddg_scan", "species_panel", "paralog_panel", "decoy_cofold", "consurf",
    "glycan", "membrane", "multi_run", "rejected", "reviews", "ipsae", "min_ipae", "complex_rmsd", "dg", "sc", "iptm_b",
    "ipae_norm", "surf_hyd", "hotspot_rmsd", "binder_rmsd", "seed_spread", "n_agree", "esm_plddt", "ciplddt",
)
#: capability flag -> the key of design.m that has to carry data for the flag to be true
METRIC_FLAGS = {k: k for k in ("ipsae", "min_ipae", "complex_rmsd", "dg", "sc", "iptm_b", "ipae_norm", "surf_hyd", "hotspot_rmsd",
                               "binder_rmsd", "seed_spread", "n_agree", "esm_plddt", "ciplddt")}

#: all keys of design.m the contract lists as required (null when the data has no value)
M_CORE = ("plddt", "ptm", "iptm", "iptm_b", "ipae_norm", "ipae_mean", "min_ipae", "ipsae", "complex_rmsd", "binder_rmsd",
          "hotspot_rmsd", "target_rmsd", "dsasa", "sc", "dg", "n_hbonds", "n_salt", "n_unsat", "hyd_frac", "hphob_patch",
          "surf_hyd", "loop_frac", "helix_frac", "strand_frac", "iface_helix", "pI", "charge", "gravy", "n_cys", "n_met",
          "n_trp", "epi_prec", "epi_rec", "mem_excess", "nterm_d", "cterm_d", "seed_spread", "n_agree")

#: The checks that compare a design with something outside the construct, as capability flags, with the words each uses for what it would have
#: tested. Which of them are off and why is written to ``profile.caps_off`` (see ``context_off``).
CONTEXT_FLAGS = {
    "olig": "the other subunits", "rest_ctx": "the rest of the protein and the cut ends", "glycan": "the glycan sites", "membrane": "the membrane plane",
}
NO_REFERENCE_NOTE = ("No reference structure was given: the construct is the frame, so this run's designs were not checked against the other subunits, "
                     "the rest of the protein, cut ends, glycan sites or a membrane. A design without a finding there has not been found clear of "
                     "them; it was not tested.")

#: a run without any hotspot residue (BindCraft can be run that way): what is not available for it
NO_HOTSPOT_NOTE = ["No hotspots recorded for this run: its designs have no hotspots touched and no intended epitope, so the epitope checks "
                   "and the epitope numbers (epi_prec, epi_rec) are left empty."]

FAILURE_LABELS = {
    "Trajectory_logits_pLDDT": "Trajectory: logits pLDDT", "Trajectory_softmax_pLDDT": "Trajectory: softmax pLDDT",
    "Trajectory_one-hot_pLDDT": "Trajectory: one-hot pLDDT", "Trajectory_final_pLDDT": "Trajectory: final pLDDT",
    "Trajectory_Contacts": "Trajectory: contacts", "Trajectory_Clashes": "Trajectory: clashes",
    "Trajectory_WrongHotspot": "Trajectory: wrong hotspot", "pLDDT": "pLDDT", "pTM": "pTM", "i_pTM": "i_pTM", "pAE": "pAE",
    "i_pAE": "i_pAE", "ShapeComplementarity": "Shape complementarity", "Surface_Hydrophobicity": "Surface hydrophobicity",
    "n_InterfaceHbonds": "Interface H-bonds", "n_InterfaceUnsatHbonds": "Unsatisfied H-bonds at the interface",
    "InterfaceAAs_K": "Interface lysines", "InterfaceAAs_M": "Interface methionines", "Binder_pLDDT": "Binder pLDDT (alone)",
    "Binder_RMSD": "Binder RMSD", "HotspotRMSD": "Hotspot RMSD", "Hotspot_RMSD": "Hotspot RMSD", "Target_RMSD": "Target RMSD",
    "n_InterfaceResidues": "Interface residues", "Binder_pTM": "Binder pTM", "Binder_pAE": "Binder pAE",
}


#: short phrases for sentences ("12 on unsatisfied H-bonds")
FAILURE_SHORT = {
    "n_InterfaceUnsatHbonds": "unsatisfied H-bonds", "ShapeComplementarity": "shape complementarity", "Surface_Hydrophobicity": "surface hydrophobicity",
    "InterfaceAAs_M": "interface methionines", "InterfaceAAs_K": "interface lysines", "n_InterfaceHbonds": "too few H-bonds",
    "Binder_RMSD": "binder RMSD", "Binder_pLDDT": "binder pLDDT", "pLDDT": "pLDDT", "i_pTM": "i_pTM", "i_pAE": "i_pAE", "pTM": "pTM",
    # names of a campaign's own filters (BindCraft 2) that have no BindCraft 1 stem
    "Unbound_Binder_pLDDT": "unbound binder pLDDT", "Interface_Residues": "too few interface residues", "Target_pLDDT": "target pLDDT",
}


class ExportError(RuntimeError):
    """The blobs do not validate against docs/schema (or a cross-check that must hold does not)."""


# ------------------------------------------------------------------------------------------------------------
# result containers
# ------------------------------------------------------------------------------------------------------------
@dataclass
class Blobs:
    """The validated outputs of one export. ``structures`` maps a file name to PDB text (binder plus target, reference frame)."""

    campaign: dict[str, Any]
    inspect: dict[str, Any]
    homolog: dict[str, Any]
    profile: dict[str, Any]
    manifest: dict[str, Any]
    structures: dict[str, str]
    warnings: list[str] = field(default_factory=list)
    validation: Any = None


# ------------------------------------------------------------------------------------------------------------
# small helpers
# ------------------------------------------------------------------------------------------------------------
def _num(x: Any) -> Any:
    if isinstance(x, (np.floating,)):
        x = float(x)
    if isinstance(x, (np.integer,)):
        x = int(x)
    return x


def _r(x: Any, nd: int) -> Optional[float]:
    """Round a number; None and NaN stay None."""
    x = _num(x)
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    return round(float(x), nd)


def clean(obj: Any) -> Any:
    """Plain JSON types: numpy scalars and arrays converted, NaN and infinity turned into null, tuples into lists."""
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return clean(obj.tolist())
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return None if math.isnan(f) or math.isinf(f) else f
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def dumps(obj: Any, indent: Optional[int] = None) -> str:
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":") if indent is None else None, indent=indent, allow_nan=False)


# BindCraft's own columns copied into the inspector records (design.src.average): the page reads these
INSPECTOR_SRC_STEMS = ("n_InterfaceUnsatHbonds", "n_InterfaceHbonds", "dSASA")


def structure_filename(pid: str) -> str:
    """File name of a design's complex in structures/: the project id, with '/' written as '__'."""
    return pid.replace("/", "__") + ".pdb"


def _load_validator():
    """scripts/validate_schema.py loaded as a module (it is a script, not a package module)."""
    root = Path(__file__).resolve().parents[2]
    path = root / "scripts" / "validate_schema.py"
    if not path.is_file():
        return None
    name = "validate_schema"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _pct(a: int, b: int) -> str:
    return f"{a} of {b}"


# ------------------------------------------------------------------------------------------------------------
# per-run blocks
# ------------------------------------------------------------------------------------------------------------
def _source_kind(text: Optional[str]) -> str:
    """settings_json, csv, inferred or unknown, from the reader's free-text source label."""
    t = (text or "unknown").lower()
    for k in ("settings_json", "csv", "inferred"):
        if t.startswith(k):
            return k
    return "unknown"


def _tool_text(run: RunData) -> str:
    """The tool as the pages name it: BindCraft with its major version ('BindCraft 1'), the same for every run of a tool whether or not
    its settings file records the exact version (that stays in run.tool_version and the manifest)."""
    m = re.match(r"v?(\d+)", (run.tool_version or "").strip())
    return f"BindCraft {m.group(1)}" if m else "BindCraft 1"


def failure_rows(run: RunData) -> list[list[Any]]:
    """[[stem, label, count], ...] of the nonzero counters of failure_csv.csv, most frequent first (the counters overlap)."""
    rows = [[k, FAILURE_LABELS.get(k, k.replace("_", " ")), int(v)] for k, v in run.failure_counters.items() if v]
    rows.sort(key=lambda r: (-r[2], r[0]))
    return rows


def filter_rule_dict(rule: config.FilterRule) -> dict[str, Any]:
    """One entry of runs[].filters. A reconstructed BindCraft 1 default keeps exactly the shape of model.filter_rule_to_dict. A rule the
    campaign recorded itself (provenance "campaign") also says where it was read (`source`), whether BindCraft counts it as mandatory,
    whether BindCraft judges it on each validation model (`judged_per_model`) and the stem its values are stored under (`value_stem`)."""
    out = model.filter_rule_to_dict(rule)
    if rule.from_campaign:
        out["metric"] = model.page_metric(rule.value_stem)
        out.update({"source": rule.source, "mandatory": rule.mandatory, "judged_per_model": rule.judged_per_model,
                    "value_stem": rule.value_stem})
    return out


def filter_set_dict(run: RunData, fs: config.FilterSet) -> dict[str, Any]:
    """runs[].filter_set, written only for a run judged by the thresholds it recorded itself (a run on BindCraft 1's reconstructed
    defaults has no such key: its rules say "counter-verified", "bounded" or "memory", and no key means exactly that)."""
    return {
        "source": fs.source, "label": fs.label, "filters_name": run.settings.filters_name, "n_rules": len(fs.rules),
        "unevaluated": [{"name": u.name, "op": u.op, "threshold": u.threshold, "mandatory": u.mandatory, "source": u.source,
                         "reason": u.reason} for u in fs.unevaluated],
        "unlimited": list(fs.unlimited),
        "note": ("These are the thresholds this campaign recorded in its own metadata, not BindCraft 1's reconstructed defaults. "
                 "They are screens, not verdicts: a design BindCraft rejected stays rejected, and where its CSV values pass a "
                 "threshold that BindCraft judged on each validation model the rejection still stands."
                 + (f" {len(fs.unevaluated)} recorded threshold{'s' if len(fs.unevaluated) != 1 else ''} on metrics Binderscope does not carry "
                    f"({', '.join(u.name for u in fs.unevaluated)}) cannot be re-checked here; BindCraft's own verdict covers {'them' if len(fs.unevaluated) != 1 else 'it'}."
                    if fs.unevaluated else "")),
    }


def project_filter_rules(project: Project) -> list[dict[str, Any]]:
    """thresholds.filters of the manifest. For a project whose runs are all judged by the default filters it is the default list, as it
    always was. Otherwise each run contributes its own rules, every entry tagged with its run."""
    sets = {rid: config.filters_for_run(run) for rid, run in project.runs.items()}
    if not any(fs.from_campaign for fs in sets.values()):
        return [model.filter_rule_to_dict(r) for r in config.DEFAULT_FILTERS]
    return [{**filter_rule_dict(r), "run": rid} for rid, fs in sets.items() for r in fs.rules]


def run_paratopes(project: Project) -> dict[str, tuple[Optional[paratope.RunParatope], Optional[str]]]:
    """{run id: (analysis, reason it is off)} for every run: (None, None) for a run that is not scaffolded (docs/paratope.md)."""
    return {rid: paratope.analyse_run(run) for rid, run in project.runs.items()}


def paratope_items(project: Project, rid: str, para: tuple[Optional[paratope.RunParatope], Optional[str]]) -> tuple[dict[str, Any], list[str]]:
    """The scaffold blocks of one campaign runs[] entry and the notes that go with them. Nothing at all for a run that is not scaffolded; a
    scaffolded run whose analysis is off carries no paratope block, a note with the reason, and a diversity_paratope block whose fold is off
    and whose paratope and loop-length axes are unknown for that reason (the epitope axis does not depend on the scaffold)."""
    res, why = para
    run = project.runs[rid]
    if res is None and why is None:
        return {}, []
    epitope = project.families["diversity"]["runs"][rid].get("epitope")
    if res is None:
        unknown = {"status": "unknown", "why": why}
        return ({"diversity_paratope": {"fold": {"status": "off", "why": families.FOLD_OFF_WHY}, "paratope": dict(unknown),
                                        "loop_length": dict(unknown), "epitope": epitope or dict(unknown)}},
                [f"Paratope analysis is off for this run: {why}. Fold diversity is off too, because every design is built on the same scaffold."])
    out: dict[str, Any] = {}
    acc, rej = res.block("accepted"), res.block("rejected")
    if acc:
        out["paratope"] = acc
    if rej:
        out["paratope_rejected"] = rej
    fams = {f["id"]: f["accepted"] for f in project.families["families"].get(rid, [])}
    out["diversity_paratope"] = families.diversity_paratope(res, fams, epitope)
    notes = list(res.layout.notes) + list(res.notes)
    if res.unaligned:
        notes.append(f"{len(res.unaligned)} of {len(res.designs) + len(res.unaligned)} designs could not be placed on the scaffold and are left out of the paratope figures.")
    return out, notes


def run_item(project: Project, rid: str, para: Optional[tuple[Optional[paratope.RunParatope], Optional[str]]] = None) -> dict[str, Any]:
    """One entry of campaign runs[] (data_contract.md 3.3). ``para`` is this run's entry of ``run_paratopes``."""
    run, frame = project.runs[rid], project.frame
    extra_blocks, extra_notes = paratope_items(project, rid, para if para is not None else paratope.analyse_run(run))
    st = run.settings
    div = project.families["diversity"]["runs"][rid]
    n_acc, n_rej, n_rej_s = len(run.accepted()), len(run.rejected()), len(run.rejected_with_structure())
    fn = run.funnel
    lo, hi = st.length_min, st.length_max
    settings = {
        "binder_name": st.binder_name, "protocol": st.protocol, "advanced_settings": st.advanced_settings_name,
        "filters_name": st.filters_name, "target_settings": st.target_settings_name, "helicity": st.helicity, "chains": st.chains,
        "length_range": [lo, hi] if lo is not None and hi is not None else None,
        "final_designs_target": st.final_designs_target,
        "hotspots": list(frame.hotspots(rid)), "hotspots_input": [int(h) for h in st.hotspots],
        "sources": {"hotspots": _source_kind(st.hotspot_source), "length_range": _source_kind(st.length_source),
                    "final_designs_target": _source_kind(st.final_designs_target_source)},
    }
    funnel = [{"key": r["field"], "label": r["label"], "n": r["value"], "status": r["status"], "note": r["note"]}
              for r in fn.stage_rows()]
    fset = config.filters_for_run(run)
    item = {
        "id": rid, "label": run.label, "reader": run.reader, "source": _tool_text(run), "tool_version": run.tool_version,
        "experimental": bool(run.experimental), "folder": Path(run.source_dir).name,
        "n_designs": {"accepted": n_acc, "rejected": n_rej, "rejected_with_structure": n_rej_s, "rows": len(run.designs),
                      "families": len(run.accepted_backbones())},
        "settings": settings,
        "target": frame.run_blocks([rid])[rid],
        "filters": [filter_rule_dict(r) for r in fset.rules],
        **({"filter_set": filter_set_dict(run, fset)} if fset.from_campaign else {}),
        "funnel": funnel, "failure": failure_rows(run),
        "n_traj": fn.trajectories_attempted, "n_traj_pass": fn.trajectory_rows, "n_scored": fn.mpnn_rows,
        "n_accepted": fn.accepted_designs,
        "traj_basis": div.get("order", {}).get("traj_basis") or families.traj_info(run)[0],
        "epoch": {"hash": div["epoch"]["hash"], "fields": div["epoch"]["components"]},
        "warnings": list(run.warnings),
        "notes": list(dict.fromkeys([*run.notes, *st.notes,                  # readers copy settings notes into run.notes
                                     *([] if frame.hotspots(rid) else NO_HOTSPOT_NOTE), *([] if frame.reference_given else [NO_REFERENCE_NOTE]),
                                     *_placement_notes(project, rid),
                                     *_filter_notes(fset), *extra_notes])),
    }
    item.update(extra_blocks)
    item["term"] = term_run_block(project, rid, termini.format_of(run))
    return item


def _filter_notes(fs: config.FilterSet) -> list[str]:
    """One note for a run judged by the thresholds it recorded itself: which claim the filter list on the pages makes."""
    if not fs.from_campaign:
        return []
    skipped = (f"; {len(fs.unevaluated)} recorded threshold{'s' if len(fs.unevaluated) != 1 else ''} on metrics Binderscope does not carry "
               f"cannot be re-checked here ({', '.join(u.name for u in fs.unevaluated)})" if fs.unevaluated else "")
    return [f"Filters: this run is judged by the {len(fs.rules)} thresholds the campaign recorded itself, not by BindCraft 1's reconstructed "
            f"defaults{skipped}. A design BindCraft rejected stays rejected even where its CSV values pass a threshold, because BindCraft judges "
            "confidence filters on each validation model and the CSV holds the mean."]


def _placement_notes(project: Project, rid: str) -> list[str]:
    """A note when some designs of the run could not be placed by the interface fit and were placed by the whole-target fit instead."""
    mine = [p["placement"] for (r, _), p in project.poses.items() if r == rid]
    fell = [p for p in mine if p.get("fit_note")]
    if not fell:
        return []
    need = int(project.frame.config["geometry"].get("iface_min_n", 12))
    return [f"{len(fell)} of {len(mine)} designs were placed by the whole-target fit: fewer than {need} target residues lie near the binder, so there is "
            "nothing to weight the fit around (prov.fit and prov.fit_note of each design say which)"]


def run_summary_block(project: Project, runs: Sequence[dict[str, Any]], profile: Mapping[str, Any]) -> dict[str, Any]:
    """The legacy project-level ``run`` block: sums over the runs, with unknown counts staying unknown."""
    n_runs = len(runs)
    pooled: Counter = Counter()
    for r in runs:
        for stem, _label, n in r["failure"]:
            pooled[stem] += n
    labels = {stem: label for r in runs for stem, label, _ in r["failure"]}
    short = (profile.get("target") or {}).get("short") or project.spec.name
    sources = sorted({r["source"] for r in runs})
    return {
        "source": " and ".join(sources), "name": f"{short}: {n_runs} BindCraft {'run' if n_runs == 1 else 'runs'} ({', '.join(r['id'] for r in runs)})",
        "n_traj": None if any(r["n_traj"] is None for r in runs) else sum(r["n_traj"] for r in runs),
        "n_traj_pass": sum(r["n_traj_pass"] or 0 for r in runs), "n_scored": sum(r["n_scored"] or 0 for r in runs),
        "n_accepted": sum(r["n_accepted"] or 0 for r in runs), "n_runs": n_runs,
        "failure": [[s, labels[s], n] for s, n in sorted(pooled.items(), key=lambda kv: (-kv[1], kv[0]))],
    }


# ------------------------------------------------------------------------------------------------------------
# designs
# ------------------------------------------------------------------------------------------------------------
def _identity_rf() -> list[float]:
    return [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]


def _terminus_env(project: Project, rid: str, pose: Mapping[str, Any]):
    """How close points of a design's own frame come to what lies outside its construct, as the target frame sees it: the other subunits, the
    rest of the chain and the glycans. The placed pose is what makes that question answerable, so the points are carried onto the reference
    chain by the design's placement first. None when no reference was given: the check is then off, never a clear result (a design with
    nothing outside its construct and a design that was never compared with it must not look alike)."""
    frame = project.frame
    if not frame.reference_given:
        return None
    R = np.asarray(pose["placement"]["R"], dtype=float).reshape(3, 3)
    t = np.asarray(pose["placement"]["t"], dtype=float)

    def env(points: Any) -> dict[str, Optional[float]]:
        res = frame.context(np.asarray(points, dtype=float).reshape(-1, 3) @ R.T + t, rid)
        m = res.metrics
        return {"olig": m.get("olig_min"), "rest": m.get("rest_min"), "glycan": (res.glycan_min or {}).get("d")}
    return env


def _terminus_block(project: Project, rid: str, d: DesignRecord, a: Any, pose: Mapping[str, Any], fmt: termini.Format) -> dict[str, Any]:
    """designs[].term and inspect.designs[].term (docs/termini.md): which end can carry a tag, and what it would cost."""
    if not fmt.applies:
        return termini.assess(a, None, None, fmt)
    cx = S.load_complex(d.structure_path, n_binder_chains=model.binder_chain_count(d))
    block = termini.assess(a, cx.binder, cx.target, fmt, _terminus_env(project, rid, pose))
    # the exit directions are in the design's own frame; the inspector draws them on the placed binder
    R = np.asarray(pose["placement"]["R"], dtype=float).reshape(3, 3)
    for end in ("n", "c"):
        if (block.get(end) or {}).get("dir"):
            block[end]["dir"] = [round(float(x), 3) for x in R @ np.asarray(block[end]["dir"], dtype=float)]
    return block


def term_run_block(project: Project, rid: str, fmt: termini.Format) -> dict[str, Any]:
    """runs[].term: what kind of molecule the run designed, whether BindCraft was asked to keep its ends free, and the tag sizes and
    constants the per-design blocks were judged with."""
    run = project.runs[rid]
    return {"v": termini.SCHEMA, "format": fmt.to_dict(), "designed_for": termini.designed_for(run), "tags": termini.tag_table(),
            "const": termini.constants(), "env": bool(project.frame.reference_given),
            "partner": bool(project.frame.partner_chains)}


def _metrics(d: DesignRecord, a: Any, pose: Mapping[str, Any], run: RunData, epitope: set[int], tc: list[int]) -> dict[str, Any]:
    """design.m: the 38 core keys (null when the data has no value) and the context keys."""
    M = d.metrics
    sp = a.seq_props or {}
    pm = pose["metrics"]
    ep = set(epitope)
    inter = len(set(tc) & ep)
    ipae = M.get("i_pAE")
    m: dict[str, Any] = {
        "plddt": M.get("pLDDT"), "ptm": M.get("pTM"), "iptm": M.get("i_pTM"), "iptm_b": M.get("i_pTM"),
        "ipae_norm": ipae, "ipae_mean": _r(ipae * 31, 1) if ipae is not None else None,
        "min_ipae": None, "ipsae": None, "ipsae_max": None, "complex_rmsd": None,
        "binder_rmsd": M.get("Binder_RMSD"), "hotspot_rmsd": M.get("Hotspot_RMSD"), "target_rmsd": M.get("Target_RMSD"),
        "dsasa": _r(a.dsasa, 1), "sc": M.get("ShapeComplementarity"), "dg": M.get("dG"),
        "n_hbonds": int(a.n_hbonds), "n_salt": int(a.n_salt), "n_unsat": int(a.n_unsat),
        "hyd_frac": _r(a.hyd_frac, 3), "hphob_patch": _r(a.hphob_patch, 1), "surf_hyd": M.get("Surface_Hydrophobicity"),
        "loop_frac": _r(a.loop_frac, 3), "helix_frac": _r(a.helix_frac, 3), "strand_frac": _r(a.strand_frac, 3),
        "iface_helix": _r(a.iface_helix, 3),
        "pI": _r(sp.get("pI"), 2), "charge": _r(sp.get("charge"), 1), "gravy": _r(sp.get("gravy"), 3),
        "n_cys": int(sp.get("n_cys", 0)), "n_met": int(sp.get("n_met", 0)), "n_trp": int(sp.get("n_trp", 0)),
        # null when the run records no hotspot (there is no intended epitope to be precise about or to recall)
        "epi_prec": ((_r(inter / len(tc), 3) if tc else 0.0) if ep else None), "epi_rec": _r(inter / len(ep), 3) if ep else None,
        "mem_excess": _r(pm["mem_excess"], 1), "nterm_d": _r(a.nterm_d, 1), "cterm_d": _r(a.cterm_d, 1),
        "seed_spread": None, "n_agree": None, "ciplddt": M.get("i_pLDDT"), "cplddt": None, "esm_plddt": None,
        "pde": None, "lis": None,
        "olig_min": _r(pm["olig_min"], 1), "olig_n3": pm["olig_n3"], "olig_n4": pm["olig_n4"], "olig_n": pm["olig_n"],
        "stalk_n": None,
        "rest_min": _r(pm["rest_min"], 1), "rest_n3": pm["rest_n3"], "rest_n4": pm["rest_n4"], "cut_n": pm["cut_n"],
        "cov_n": pm["cov_n"], "place_rmsd": _r(pose["placement"]["rmsd"], 2),
        # how well the placement fits where the binder sits (the part of the target within about 13 A of it) and what the whole-target fit would give
        "place_iface_rmsd": _r(pose["placement"].get("iface_rmsd"), 2), "place_iface_n": pose["placement"].get("iface_n"),
        "place_global_rmsd": _r(pose["placement"].get("global_rmsd"), 2),
    }
    for k in M_CORE:
        m.setdefault(k, None)
    return m


def _src(d: DesignRecord, run: RunData, names: Optional[Mapping[str, str]] = None) -> dict[str, Any]:
    """design.src: BindCraft's own numbers exactly as the CSV wrote them (CSV scales).

    `names` is {campaign filter name: stem the value is stored under} for a campaign rule whose name is not the stem (BindCraft 2's
    Interface_Residues is stored as n_InterfaceResidues). The value is repeated under the campaign's own name, which is what its CSV
    calls the column, so that the page can find the value of a filter by the filter's name."""
    average = {k: v for k, v in d.metrics.items() if v is not None}
    for name, stem in (names or {}).items():
        if stem in average:
            average.setdefault(name, average[stem])
    return {
        "tool": _tool_text(run),
        "average": average,
        "models": {str(k): {s: v for s, v in vals.items() if v is not None} for k, vals in sorted(d.models.items())},
        "interface_aas": {k: dict(v) for k, v in d.interface_aas.items()},
    }


def _prov(d: DesignRecord, run: RunData, pose: Optional[Mapping[str, Any]] = None) -> dict[str, Any]:
    rel = None
    if d.structure_path:
        src = d.extras.get("source_structure") or d.structure_path          # a one-based copy of a BindCraft 2 file names its original
        try:
            rel = str(Path(src).resolve().relative_to(Path(run.source_dir).resolve()))
        except ValueError:
            rel = Path(src).name
    out = {
        "run": run.run_id, "binder_name": d.binder_name, "backbone": d.backbone_id, "length": d.length, "seed": d.seed,
        "mpnn": d.mpnn_index, "protocol": d.protocol, "model": d.structure_model, "model_file": rel, "rank": d.rank,
        "csv": "final_design_stats.csv" if d.status == "accepted" else "mpnn_design_stats.csv",
        "design_time_s": d.design_time_s, "notes": d.notes, "structure": "present" if d.structure_path else "missing",
        "marked": bool(d.xlsx_marked),
    }
    pl = (pose or {}).get("placement") or {}
    if pl.get("fit"):
        out["fit"] = pl["fit"]                    # how the design was placed on the reference chain: interface, all or core
        if pl.get("fit_note"):
            out["fit_note"] = pl["fit_note"]
    return out


@dataclass
class _Ctx:
    """Everything the per-design builders share."""

    project: Project
    profile: dict[str, Any]
    mismatches: list[str]
    near: dict[ProjectKey, bool]
    fold_label: dict[str, str]
    q: dict[ProjectKey, float]


def _design_records(project: Project, profile: dict[str, Any], mismatches: list[str],
                    para: Optional[Mapping[str, tuple[Optional[paratope.RunParatope], Optional[str]]]] = None) -> tuple[list[dict[str, Any]], dict[ProjectKey, bool]]:
    """Campaign designs[] for every accepted design and every rejected design with a structure. ``para`` (``run_paratopes``) adds ``para`` to
    the designs of a scaffolded run that were placed on the scaffold."""
    runs, frame, pids = project.runs, project.frame, project.project_ids
    assign = project.modes["assign"]
    labels = project.families["fold_labels"]
    rej_fold = project.families["rejected_fold"]
    out: list[dict[str, Any]] = []
    inspected: dict[ProjectKey, bool] = {}
    for rid, run in runs.items():
        fmt = termini.format_of(run)
        off = run.target.uniprot_offset
        names = {r.stem: r.metric for r in config.filters_for_run(run).rules if r.metric}
        hs = list(frame.hotspots(rid))
        epi = set(frame.epitope(rid))
        qs = families.quality_scores(run)
        placed = ((para or {}).get(rid) or (None, None))[0]
        para_blobs = placed.design_blobs() if placed is not None else {}
        for d in run.accepted() + run.rejected_with_structure():
            key = (rid, d.id)
            pid = pids[key]
            pose, a = project.poses[key], project.analyses[key]
            tc = sorted(int(c) + off for c in a.tc)
            if tc != sorted(pose["contacted"]):
                mismatches.append(f"{pid}: contact set of the interface module differs from the target frame's "
                                  f"({len(tc)} against {len(pose['contacted'])} residues)")
            hot = sorted(set(tc) & set(hs))
            if hot != sorted(int(h) + off for h in a.hot) or hot != sorted(pose["hot"]):
                mismatches.append(f"{pid}: hotspots touched differ between the interface module and the target frame")
            rec_fail = project.rescue.get(key)
            near = bool(d.status == "rejected" and rec_fail and rec_fail["nm"] and rec_fail["nm"]["near"])
            inspected[key] = d.status == "accepted" or near
            bb = run.backbone(d.backbone_id)
            m = _metrics(d, a, pose, run, epi, tc)
            sc = labels.get(pid) if d.status == "accepted" else (rej_fold.get(pid) or {}).get("label", "unclustered")
            asg = assign.get(pid) or {}
            rec: dict[str, Any] = {
                "id": pid, "name": d.id, "run": rid, "status": d.status, "rank": d.rank, "bb": d.backbone_id, "fam": d.backbone_id,
                "s": d.mpnn_index if d.mpnn_index is not None else 1, "seed": d.seed, "traj": bb.traj if bb else 0,
                "sc": sc or "unclustered", "site": asg.get("site", "A1"), "mode": asg["mode"], "len": d.length, "seq": d.sequence,
                "flip": False, "ca": [_r(x, 1) for x in pose["ca"]], "rf": _identity_rf(),
                "view": pose["view"], "vc": [_r(x, 1) for x in pose["ic"]], "emb": [_r(x, 4) for x in asg["emb"]],
                "tc": tc, "bc": [int(b) - 1 for b in a.bc], "nc": int(a.nc), "hot": hot, "hot_of": len(hs), "m": m,
                "glycan": {k: _r(v, 1) for k, v in pose["glycan"].items()},
                "models": [{"model": "AF2 (BindCraft)", "ipsae": None, "iptm": m["iptm"], "rmsd": m["hotspot_rmsd"]}],
                "liab": liabilities.to_campaign(a.liab),
                "inspect": pid if inspected[key] else None,
                "pl": [int(round(x)) if x is not None else 0 for x in a.plddt_res],
                "prov": _prov(d, run, pose), "src": _src(d, run, names),
                "fail": [] if d.status == "accepted" else (rec_fail["fail"] if rec_fail else []),
            }
            if d.status == "accepted":
                rec["q"] = qs[d.id]
            elif rec_fail and rec_fail["nm"]:
                rec["nm"] = rec_fail["nm"]
            if d.id in para_blobs:
                rec["para"] = para_blobs[d.id]
            if inspected[key]:      # the designs a person can reach from the board or the inspector; the rest of the rejected designs carry none, not an empty one
                rec["term"] = clean(_terminus_block(project, rid, d, a, pose, fmt))
            out.append(rec)
    return out, inspected


def _near_binder_rmsd(pl: Mapping[str, Any]) -> float:
    """The CA RMSD that matters for the context checks: over the residues near the binder (the whole target when a record has none)."""
    v = pl.get("iface_rmsd")
    return float(pl["rmsd"] if v is None else v)


def placement_warnings(project: Project, designs: Sequence[Mapping[str, Any]]) -> list[str]:
    """One warning per run whose designs fit the crystal chain worse than ``place_warn_rmsd`` near the binder (CONTRACTS open item 5): the
    design's own target differs more from the reference chain at the interface, so its context distances are less certain. The trimmed (core)
    RMSD says whether the rigid part of the target still fits."""
    geo = project.frame.config["geometry"]
    thr = float(geo["place_warn_rmsd"])
    out = []
    for rid in project.runs:
        mine = [d for d in designs if d["run"] == rid]
        parts = []
        for status in ("accepted", "rejected"):
            ds = [d for d in mine if d["status"] == status]
            bad = [d for d in ds if _near_binder_rmsd(project.poses[(rid, d["name"])]["placement"]) > thr]
            if bad:
                parts.append(f"{len(bad)} of {len(ds)} {status}")
        if parts:
            core = max(project.poses[(rid, d["name"])]["placement"]["core_rmsd"] for d in mine)
            out.append(f"{rid}: placement RMSD near the binder above {thr:g} \u00c5 for {' and '.join(parts)} designs: their own target differs more from "
                       f"the crystal chain at the interface, so distances to the other subunit and the rest of the protein are less certain (the rigid "
                       f"core of every design fits within {core:.2f} \u00c5)")
    return out


def _sibs(project: Project, rid: str, d: DesignRecord) -> list[dict[str, Any]]:
    """Sequences on the same backbone (the design included): accepted ones by rank, then the rejected ones by index."""
    run = project.runs[rid]
    pids = project.project_ids
    sibs = [x for x in run.designs if x.backbone_id == d.backbone_id]
    sibs.sort(key=lambda x: (x.status != "accepted", x.rank if x.rank is not None else 10 ** 6, x.mpnn_index or 0, x.id))
    first = sibs[0].sequence
    off = run.target.uniprot_offset
    hs = set(project.frame.hotspots(rid))
    out = []
    for x in sibs:
        key = (rid, x.id)
        a = project.analyses.get(key)
        pose = project.poses.get(key)
        diff = [i + 1 for i, (p, q) in enumerate(zip(x.sequence, first)) if p != q]
        M = x.metrics
        out.append({
            "id": pids[key], "seq": x.sequence, "iptm": M.get("i_pTM"), "ipae_norm": M.get("i_pAE"), "plddt": M.get("pLDDT"),
            "unsat": int(a.n_unsat) if a else None, "hb": int(a.n_hbonds) if a else None, "dg": M.get("dG"),
            "sc": M.get("ShapeComplementarity"),
            "hot": len(hs & {int(c) + off for c in a.tc}) if a else None,
            "ndiff": len(diff), "diff": diff, "status": x.status,
        })
    return out


def _pdb_binder(ss: str, seq: str, placed_b: S.Atoms) -> str:
    """pdb_b: HELIX and SHEET records first, then the binder's heavy atoms (B-factor = pLDDT) in the reference frame.

    The inspector indexes residues 1..n against ``ss`` and ``seq``. A binder of
    one chain is already numbered that way. A binder of several (an scFv's VH
    and VL) has each chain numbered from 1, so the numbers repeat: the residues
    are renumbered continuously across the chains, and each chain keeps its own
    letter so the viewer still draws two chains and a TER between them.
    """
    nums = [int(x) for x in placed_b.residues.number]
    chains = [str(c) for c in placed_b.residues.chain]
    runs: list[list[int]] = []                       # residue index ranges, one per chain in file order
    for i, c in enumerate(chains):
        if runs and chains[runs[-1][0]] == c:
            runs[-1][1] = i + 1
        else:
            runs.append([i, i + 1])
    for a, b in runs:
        if nums[a:b] != list(range(1, b - a + 1)):
            where = f" of chain {chains[a]}" if len(runs) > 1 else ""
            raise ExportError(f"binder residues{where} are not numbered 1..n; the inspector expects them to be")
    names = [S.one_to_three(c) for c in seq]
    lines = TF.ss_records(ss, list(range(1, len(nums) + 1)), names, "B" if len(runs) == 1 else chains[0])
    serial = 1
    for a, b in runs:
        lo, hi = int(placed_b.residues.start[a]), int(placed_b.residues.stop[b - 1])
        atoms = placed_b.select(np.arange(lo, hi))
        lines += S.to_pdb_lines(atoms, serial, chain="B" if len(runs) == 1 else chains[a], resnum_offset=a)
        lines.append("TER")
        serial += len(atoms)
    if len(runs) == 1:
        lines = lines[:-1] + ["TER"]
    return "\n".join(lines) + "\n"


def is_pdb_id(name: Optional[str]) -> bool:
    """A four-character PDB entry id (digit, then three letters or digits)."""
    return bool(name) and len(str(name)) == 4 and str(name)[0].isdigit() and str(name).isalnum()


def reference_name(frame: TF.TargetFrame) -> str:
    """The reference structure as labels name it: ``PDB 1ABC`` for a PDB entry, the id as given for any other file."""
    ref = frame.meta["reference"]["pdb"]
    return f"PDB {ref}" if is_pdb_id(ref) else str(ref)


def _placement_remark(pl: Mapping[str, Any]) -> str:
    """The REMARK line that says how the structure was placed on the reference chain and how well the target fits (CA RMSD, residues)."""
    whole = f"{pl['rmsd']} A over {pl['n']}"
    near = f"{pl['iface_rmsd']} A over {pl['iface_n']} residues" if pl.get("iface_rmsd") is not None else None
    fit = pl.get("fit")
    if fit == "interface" and near:
        return f"REMARK 250 placement: fitted around the interface, CA RMSD {near} (whole target {whole})"
    if fit == "core":
        return (f"REMARK 250 placement: fitted on the rigid core of the target, CA RMSD {pl['core_rmsd']} A over {pl['core_n']} residues (whole target {whole}"
                + (f"; near the interface {near}" if near else "") + ")")
    if fit in ("all", "interface"):
        return (f"REMARK 250 placement: whole-target fit, CA RMSD {whole} residues" + (f" (near the interface {near})" if near else "")
                + (f"; {pl['fit_note']}" if pl.get("fit_note") else ""))
    return f"REMARK 250 placement: CA RMSD {whole} shared residues (rigid fit of the target on the reference chain)"


def _pdb_complex(d: DesignRecord, rid: str, pid: str, ss: str, seq: str, placed_t: S.Atoms, placed_b: S.Atoms, off: int,
                 pose: Mapping[str, Any], frame: TF.TargetFrame) -> str:
    """structures/<id>.pdb: the design's own target (chain A, UniProt numbers) and the binder (chain B), reference frame."""
    pl = pose["placement"]
    numbering = f"UniProt {frame.accession} numbering" if frame.accession else "numbering of the reference structure"
    head = [
        "REMARK 250 Binderscope structure: a BindCraft design placed on the reference chain",
        f"REMARK 250 design {pid}  run {rid}  AlphaFold2 model {d.structure_model}",
        f"REMARK 250 frame: {reference_name(frame)} chain {frame.meta['chain']} ({numbering}); chain A is this design's own target (the residues of its run),",
        "REMARK 250 chain B is the binder (residues 1 to n, B-factor = pLDDT). Hydrogens are not written.",
        _placement_remark(pl),
        "REMARK 250 A predicted complex is a screen, not a measurement.",
    ]
    nums_b = [int(x) for x in placed_b.residues.number]
    names_b = [S.one_to_three(c) for c in seq]
    ss_lines = TF.ss_records(ss, nums_b, names_b, "B")
    lines = head + ss_lines
    lines += S.to_pdb_lines(placed_t, 1, chain="A", resnum_offset=off)
    lines.append("TER")
    lines += S.to_pdb_lines(placed_b, len(placed_t) + 2, chain="B")
    lines += ["TER", "END"]
    return "\n".join(lines) + "\n"


def _inspector_designs(project: Project, profile: dict[str, Any], camp_designs: list[dict[str, Any]], inspected: dict[ProjectKey, bool]
                       ) -> tuple[dict[str, dict[str, Any]], list[str], dict[str, str]]:
    """inspect.designs (tag = project id), the display order and the structure files."""
    runs, frame, pids = project.runs, project.frame, project.project_ids
    by_pid = {r["id"]: r for r in camp_designs}
    modes = project.modes["modes"]
    ins: dict[str, dict[str, Any]] = {}
    structures: dict[str, str] = {}
    order_acc: list[tuple] = []
    order_rej: list[tuple] = []
    for rid, run in runs.items():
        off = run.target.uniprot_offset
        for d in run.accepted() + run.rejected_with_structure():
            key = (rid, d.id)
            if not inspected.get(key):
                continue
            pid = pids[key]
            rec = by_pid[pid]
            pose, a = project.poses[key], project.analyses[key]
            cx = S.load_complex(d.structure_path, n_binder_chains=model.binder_chain_count(d))
            R = np.asarray(pose["placement"]["R"], dtype=float).reshape(3, 3)
            t = np.asarray(pose["placement"]["t"], dtype=float)
            placed_b, placed_t = cx.binder.transformed(R, t), cx.target.transformed(R, t)
            seq = d.sequence
            # The CSV separates the chains of a multi-chain binder with "/", the
            # structure does not, so they are compared chain by chain: that also
            # catches chains read in the wrong order, which comparing one joined
            # string would not.
            want = seq.split("/")
            got = [cx.binder.select(cx.binder.chain == c).sequence for c in cx.binder_chains] or [cx.binder.sequence]
            if len(want) != len(got) or any(a != b for a, b in zip(want, got)):
                where = next((i for i, (a, b) in enumerate(zip(want, got)) if a != b), None)
                detail = (f" (chain {cx.binder_chains[where]}: {len(got[where])} residues in the structure, "
                          f"{len(want[where])} in the CSV)" if where is not None and where < len(cx.binder_chains)
                          else f" ({len(got)} chains in the structure, {len(want)} in the CSV)")
                raise ExportError(f"{pid}: the binder sequence of the structure differs from the CSV{detail}")
            m = rec["m"]
            tspan = lambda mp, fill: frame.union_array({int(k): v for k, v in mp.items()}, fill=fill)  # noqa: E731 - construct numbers in, union order out
            octx = clean(pose["octx"])
            ctx = clean(pose["ctx"]) or None       # a design that was not placed against anything outside its construct has no context block (null, not {})
            tag = pid
            ins[tag] = {
                "tag": tag, "id": pid, "label": pid, "run": rid, "status": d.status, "rank": d.rank, "prov": rec["prov"], "bb": d.backbone_id,
                "mode": rec["mode"], "mode_name": modes[rec["mode"]]["name"], "scaffold": rec["sc"], "site": rec["site"],
                "modality": INSPECTOR_MODALITY, "seq": seq, "n": d.length, "ss": a.ss, "pdb_b": _pdb_binder(a.ss, seq, placed_b),
                "plddt_res": [_r(x, 1) for x in a.plddt_res],
                "plddt_t": [None if v is None else _r(v, 1) for v in tspan(a.plddt_t, None)],
                "rsa_res": [_r(x, 3) for x in a.rsa_res], "dsasa_res": [_r(x, 1) for x in a.dsasa_res],
                "ddg_res": None, "t_dsasa": [_r(v, 1) for v in tspan(a.t_dsasa, 0.0)], "t_ddg": None,
                "contacts": [{"t": int(c["t"]) + off, "b": int(c["b"]), "d": _r(c["d"], 2), "n": int(c["n"])} for c in a.contacts],
                "hbonds": [{"t": int(h["t"]) + off, "tn": h["tn"], "ta": h["ta"], "b": int(h["b"]), "bn": h["bn"], "ba": h["ba"],
                            "d": _r(h["d"], 2), "k": h["k"]} for h in a.hbonds],
                "salt": [{"t": int(s["t"]) + off, "tn": s["tn"], "b": int(s["b"]), "bn": s["bn"], "d": _r(s["d"], 2)} for s in a.salt],
                "buns": [{"side": b["side"], "resi": int(b["resi"]) + (off if b["side"] == "t" else 0), "resn": b["resn"], "atom": b["atom"]}
                         for b in a.buns],
                "patches": [{"area": _r(p["area"], 1), "res": [int(x) for x in p["res"]]} for p in a.patches],
                "voids": [], "dsasa": m["dsasa"], "hyd_frac": m["hyd_frac"],
                "t_contact": list(rec["tc"]), "b_contact": [int(b) for b in a.bc], "hot": list(rec["hot"]),
                "ssbonds": [[int(i), int(j)] for i, j in a.ssbonds], "ssbond_basis": a.ssbond_basis,
                "disulfides": [{"a": int(p["a"]), "b": int(p["b"]), "d": p["d"], "band": p["band"]} for p in a.disulfides],
                "nterm_d": m["nterm_d"], "cterm_d": m["cterm_d"],
                "mem_margin": _r(-m["mem_excess"], 1) if m["mem_excess"] is not None else None,
                "glycan": dict(rec["glycan"]), "pae": None, "pae_n": None,
                "poses": [{"model": "AF2 (BindCraft)", "ipsae": None, "iptm": m["iptm"], "rmsd": m["hotspot_rmsd"],
                           "R": [1, 0, 0, 0, 1, 0, 0, 0, 1], "t": [0, 0, 0]}],
                "seeds": [], "sibs": _sibs(project, rid, d), "metrics": {**m, "length": d.length},
                "liab": liabilities.to_inspector(a.liab), "flip": False, "center": [_r(x, 2) for x in pose["ic"]],
                "ctx": ctx, "octx": octx,
            }
            # the filters a rejected design fails (copy of the campaign record): the inspector names them without loading the campaign blob
            ins[tag]["fail"] = list(rec.get("fail") or [])
            if rec.get("nm"):
                ins[tag]["nm"] = rec["nm"]
            # BindCraft's own counts of the quantities the inspector recomputes (Rosetta's unsatisfied hydrogen bonds differ from a geometric count):
            # shown next to ours and used for the unsatisfied-polars meter
            own = {k: v for k, v in (rec.get("src") or {}).get("average", {}).items() if k in INSPECTOR_SRC_STEMS}
            if own:
                ins[tag]["src"] = {"average": own}
            # which scaffold segment each interface residue falls in, for a design built
            # on a scaffold: copied from the campaign record for the same reason as the
            # filters above, since the inspector never loads the campaign blob. Absent
            # for an unscaffolded design, and the page draws nothing rather than zeroes.
            if rec.get("para"):
                ins[tag]["para"] = rec["para"]
            if rec.get("term"):
                ins[tag]["term"] = rec["term"]
            structures[structure_filename(pid)] = _pdb_complex(d, rid, pid, a.ss, seq, placed_t, placed_b, off, pose, frame)
            row = (list(runs).index(rid), d.rank if d.rank is not None else 10 ** 6, d.id)
            if d.status == "accepted":
                order_acc.append((row, tag))
            else:
                nm = rec.get("nm") or {}
                order_rej.append(((-(nm.get("score") or 0.0), row), tag))
    order = [t for _, t in sorted(order_acc)] + [t for _, t in sorted(order_rej)]
    return ins, order, structures


# ------------------------------------------------------------------------------------------------------------
# profile
# ------------------------------------------------------------------------------------------------------------
def completed_profile(project: Project) -> dict[str, Any]:
    """A copy of the project's profile with what the exporter reads filled in when the profile lacks it (id, kind, target name, the
    limits of checks.context from the frame's defaults), the target's accession taken from the frame and, for the built-in default
    profile (whose target is a placeholder), the name, short name, gene and organism taken from the UniProt record when an accession
    was given. A profile that names its target (any other id) keeps its own words."""
    base = json.loads(json.dumps(project.profile))
    frame = project.frame
    base.setdefault("schema", "binderscope.profile/1")
    base.setdefault("id", "custom")
    base.setdefault("kind", "run")
    tg = base.setdefault("target", {})
    tg.setdefault("short", tg.get("gene") or project.spec.name)
    tg.setdefault("name", tg["short"])
    ctx = base.setdefault("checks", {}).setdefault("context", {})        # a profile without context limits gets the frame's defaults
    for part, vals in TF.DEFAULT_CONFIG["checks"].items():
        have = ctx.setdefault(part, {})
        for k, v in vals.items():
            have.setdefault(k, v)
    if frame.accession:
        tg["uniprot"] = frame.accession
    if base.get("id") == "default" and frame.accession:
        ann = frame.meta["annotations"]
        if ann.get("name"):
            tg["name"] = ann["name"]
        short = ann.get("gene") or frame.accession
        tg["short"] = short
        for k, v in (("gene", ann.get("gene")), ("organism", ann.get("organism"))):
            if v:
                tg[k] = v
        base.setdefault("specificity", {})["subject"] = short
    tg.setdefault("uniprot", "")
    return base


def _name_segments(target: dict[str, Any], profile_frame: Mapping[str, Any] | None) -> None:
    """The profile names its own frame segments (its words win over the CATH names of the UniProt annotation): the labels and phrases of the
    target block's frame strip and its per-residue domain labels (``dom``) take the profile's names for the segments with the same range."""
    segs = (profile_frame or {}).get("segments") or []
    if not segs or not (target.get("frame") or {}).get("segments"):
        return
    by = {(s["start"], s["end"], s["kind"]): s for s in segs}
    for d in target["frame"]["segments"]:
        m = by.get((d["start"], d["end"], d["kind"]))
        if m:
            d["label"] = m["label"]
            if "phrase" in m:
                d["phrase"] = m["phrase"]
    target["dom"] = [next((s["label"] for s in target["frame"]["segments"] if s["start"] <= u <= s["end"]), "") for u in target["resi"]]


def _same_frame(given: Mapping[str, Any], derived: Mapping[str, Any]) -> bool:
    """A profile's frame strip agrees with the one derived from UniProt when ranges, kinds and sites match; the profile names its own segments."""
    g, d = clean(given), clean(derived)
    seg = lambda f: [[x["start"], x["end"], x["kind"]] for x in f.get("segments") or []]  # noqa: E731
    return all(g.get(k) == d.get(k) for k in ("native_range", "glycan_sites", "disulfides")) and seg(g) == seg(d)


def _olig_off_note(frame: Any) -> str:
    """Why the other-subunit check is off, from what the reference actually holds.

    A reference with no second chain, one whose other chains were deliberately not
    counted, and one the user asked to ignore are three different situations. The
    first is a fact about the protein; the others are choices, and reporting them as
    "this protein has no other subunit" would be a claim the file does not support.
    """
    comp = getattr(frame, "composition", None)
    if not comp:                                        # a frame saved before compositions were kept
        return "the reference structure has no other chain of this target, so there is no other subunit to overlap"
    partners = comp.get("partners") or {}
    if partners.get("mode") == "none":
        return "partner chains were switched off for this review, so no other subunit was checked against"
    # Chains that exist but were not counted: a different molecule, or one outside the assembly.
    skipped = [r for r in (comp.get("chains") or [])
               if r.get("role") in ("other_molecule", "outside_assembly") and not r.get("used")]
    if skipped:
        names = ", ".join(str(r.get("chain")) for r in skipped[:4])
        return (f"the reference holds other chains ({names}) but none of them is another subunit of this target, "
                "so there was nothing to overlap; name them as partners if they belong to the complex")
    return "the reference structure has no other chain of this target, so there is no other subunit to overlap"


def context_off(project: Project, caps: Mapping[str, bool], present: Mapping[str, bool], wanted: Mapping[str, Any]) -> dict[str, dict[str, str]]:
    """``profile.caps_off``: for every context flag that is false, why. Three situations look the same in ``caps`` and are not the same:

    * ``no_reference``: no full-length reference structure was given, so nothing outside the construct was ever known. The checks did not run.
    * ``reference_lacks``: a reference was given but it holds nothing for the check (no other subunit, nothing outside the construct, no sugar at a
      site inside it, no transmembrane region). The check ran on nothing, which is not the same as passing.
    * ``profile``: the profile switched the flag off.

    ``note`` is the specific reason, in words the Overview prints under its chip."""
    frame = project.frame
    out: dict[str, dict[str, str]] = {}
    for flag in CONTEXT_FLAGS:
        if caps.get(flag):
            continue
        if wanted.get(flag, True) is False and present.get(flag):
            out[flag] = {"why": "profile", "note": "switched off by the profile"}
        # (see _olig_off_note: "the reference has no other chain" is only one of the
        #  reasons the other-subunit check can be off, and saying it for the others
        #  would describe the file wrongly)
        elif not frame.reference_given:
            out[flag] = {"why": "no_reference", "note": "needs a full-length reference structure, and none was given"}
        else:
            out[flag] = {"why": "reference_lacks", "note": {
                "olig": _olig_off_note(frame),
                "rest_ctx": "the reference chain holds nothing outside the construct, and no run deleted a stretch from it",
                "glycan": "no N-glycan site lies inside the construct (none in the reference structure, none in UniProt)",
                "membrane": "no transmembrane region is known: the UniProt topology is read only when an accession is given",
            }[flag]}
    return out


def resolve_context_hard(project: Project, base: Mapping[str, Any]) -> dict[str, Any]:
    """Does the Context check decide a verdict (hard) or only inform it (advisory)?

    Advisory by default, because the check places the binder against the reference as one static model. Against a predicted reference that is
    a prediction of what lies outside the construct, and a region AlphaFold cannot place (a flexible linker) lands at an arbitrary position:
    on a real campaign every design collided with such a region and the check alone turned the BindCraft-accepted designs to NO PASS, which also
    hid the findings that were real. It becomes hard when the reference is an experimental structure (the geometry is then measured), or when the
    user says so. The setting is ``--context-check`` or ``profile.checks.context.hard`` (``auto``, ``hard``, ``advisory``; true and false are
    read as hard and advisory). ``cannot tell`` is treated as predicted."""
    raw = project.spec.context_hard
    if raw in (None, ""):
        raw = ((base.get("checks") or {}).get("context") or {}).get("hard", "auto")
    setting = {True: "hard", False: "advisory"}.get(raw, str(raw).lower())
    if setting not in ("auto", "hard", "advisory"):
        raise ExportError(f"checks.context.hard / --context-check must be auto, hard or advisory, not {raw!r}")
    cls = project.frame.reference_class
    if setting == "auto" and cls["kind"] == "none":
        return {"hard": False, "setting": setting, "basis": "the Context check was not run: no reference structure was given"}
    where = "--context-check" if project.spec.context_hard else "the profile"
    if setting == "hard":
        return {"hard": True, "setting": setting, "basis": f"set to hard by {where}"}
    if setting == "advisory":
        return {"hard": False, "setting": setting, "basis": f"set to advisory by {where}"}
    if cls["treated_as"] == "experimental":
        return {"hard": True, "setting": setting, "basis": "the reference is an experimental structure, so the geometry it gives is measured"}
    what = {"predicted": "the reference is a predicted model", "unknown": "the file does not say whether the reference is experimental or predicted, "
            "so it is read as predicted"}[cls["kind"]]
    return {"hard": False, "setting": setting, "basis": f"{what}, so the placement of the rest of the protein is a prediction, not a measurement"}


def build_profile(project: Project, campaign_parts: Mapping[str, Any], warnings: list[str]) -> dict[str, Any]:
    """The project profile: the profile completed from the data. ``caps`` is written from the data (a flag is true exactly
    when its data is present; a user profile can only switch one off), the target construct, length, reference and numbering come
    from the frame, and a frame strip that differs from the one derived from UniProt and the reference is reported."""
    base = completed_profile(project)
    frame = project.frame
    designs, target, homolog, fold, runs = (campaign_parts[k] for k in ("designs", "target", "homolog", "fold", "runs"))
    have = lambda key: any((d["m"].get(key) is not None) for d in designs)  # noqa: E731
    bbs = Counter(d["bb"] for d in designs)
    present = {
        "outcomes": False, "weights": False,
        "olig": bool((target.get("olig") or {}).get("neigh")) and have("olig_min"),
        "rest_ctx": have("rest_min") and bool(((campaign_parts.get("rest") or target.get("rest")) or {}).get("traces")),
        "flag_region": False, "footprints": False, "has_pae": False,
        "has_traj": any(d["traj"] > 0 for d in designs), "traj_total": all(r["n_traj"] is not None for r in runs),
        "has_fold": bool(fold), "families": any(n > 1 for n in bbs.values()), "has_boltz": False, "models": False,
        "has_seeds": False, "ddg_scan": False, "species_panel": bool(homolog.get("orth")), "paralog_panel": bool(homolog.get("par")),
        "decoy_cofold": False, "consurf": False, "glycan": bool(target["glycan"]), "membrane": have("mem_excess"),
        "multi_run": len(runs) > 1, "rejected": any(d["status"] == "rejected" for d in designs), "reviews": True,
    }
    for flag, key in METRIC_FLAGS.items():
        present[flag] = have(key)
    wanted = base.get("caps") or {}
    caps = {}
    for flag in CAP_FLAGS:
        caps[flag] = bool(present.get(flag)) and wanted.get(flag, True) is not False
    if project.spec.consurf and homolog.get("consurf"):
        caps["consurf"] = True
    base["caps"] = caps
    ctxres = resolve_context_hard(project, base)
    cc = base.setdefault("checks", {}).setdefault("context", {})
    cc["hard"], cc["hard_setting"], cc["hard_basis"] = ctxres["hard"], ctxres["setting"], ctxres["basis"]
    off = context_off(project, caps, present, wanted)
    if off:
        base["caps_off"] = off
    else:
        base.pop("caps_off", None)
    tg = base.setdefault("target", {})
    first, last = target["resi"][0], target["resi"][-1]
    short = tg.get("short") or tg.get("gene") or project.spec.name
    tg.setdefault("construct", {})
    tg["construct"] = {**tg["construct"], "first": first, "last": last}
    tg["construct"].setdefault("label", "construct")
    tg["length"] = frame.length
    given = frame.reference_given
    # Without a reference there is no structure to name: the pages print "Reference X" from pdb, so it stays empty rather than naming the construct.
    tg["reference"] = {**(tg.get("reference") or {}), "pdb": frame.meta["reference"]["pdb"] if given else "", "chain": frame.meta["chain"],
                       "assembly": frame.meta["reference"]["assembly"], "given": given}
    cls = frame.reference_class
    tg["reference"].update({"kind": cls["kind"], "treated_as": cls["treated_as"], "evidence": list(cls["evidence"]),
                            "method": cls.get("method"), "resolution": cls.get("resolution")})
    if given:
        tg["reference"].setdefault("label", frame.meta["reference"]["label"])
    else:
        tg["reference"]["label"] = frame.meta["reference"]["label"]
    tg["numbering"] = {**(tg.get("numbering") or {}), "input_offset": -frame.offset}
    tg["numbering"].setdefault("input_label", "BindCraft input")
    tg["label"] = f"{short} {first}-{last}"
    tg["span_label"] = f"{short} {first}–{last}"
    derived = frame.profile_frame()
    if not base.get("frame"):
        base["frame"] = derived
    elif frame.accession and not _same_frame(base["frame"], derived):
        warnings.append("profile.frame differs from the frame strip derived from UniProt and the reference structure; the profile's own strip is used")
    base["version"] = base.get("version", 1)
    return base


# ------------------------------------------------------------------------------------------------------------
# counts against the CSVs
# ------------------------------------------------------------------------------------------------------------
def _child(root: Path, name: str) -> Optional[Path]:
    for p in root.iterdir():
        if p.name.lower() == name.lower():
            return p
    return None


def _read_csv(path: Optional[Path]):
    import pandas as pd
    if path is None:
        return None
    try:
        return pd.read_csv(path)
    except Exception:                                            # an unreadable table is "unknown", not a mismatch
        return None


def _counts_bindcraft1(run: Any, exp: Sequence[Mapping[str, Any]]) -> dict[str, tuple]:
    root = Path(run.source_dir)
    final = _read_csv(_child(root, config.CSV_FINAL))
    mpnn = _read_csv(_child(root, config.CSV_MPNN))
    traj = _read_csv(_child(root, config.CSV_TRAJECTORY))
    acc_dir, rej_dir = _child(root, config.ACCEPTED_DIR), _child(root, config.REJECTED_DIR)
    n_acc_pdb = len([p for p in acc_dir.glob("*.pdb")]) if acc_dir else 0
    n_rej_pdb = len([p for p in rej_dir.glob("*.pdb")]) if rej_dir else 0
    n_acc, n_rej = sum(d["status"] == "accepted" for d in exp), sum(d["status"] == "rejected" for d in exp)
    fams = {str(x).rsplit("_mpnn", 1)[0] for x in final["Design"]} if final is not None else set()
    funnel = {r["field"]: r["value"] for r in run.funnel.stage_rows()}
    return {
        "accepted designs": (len(final) if final is not None else None, n_acc, "final_design_stats.csv rows against exported accepted designs"),
        "accepted structures": (n_acc_pdb, n_acc, "accepted/*.pdb files against exported accepted designs"),
        "rejected structures": (n_rej_pdb, n_rej, "rejected/*.pdb files against exported rejected designs"),
        "sequence rows": (len(mpnn) if mpnn is not None else None, funnel.get("mpnn_rows"), "mpnn_design_stats.csv rows against the funnel"),
        "trajectory rows": (len(traj) if traj is not None else None, funnel.get("trajectory_rows"), "trajectory_stats.csv rows against the funnel"),
        "families": (len(fams), len(run.accepted_backbones()), "distinct backbones of final_design_stats.csv against the families"),
    }


def _counts_bindcraft2(run: Any, exp: Sequence[Mapping[str, Any]]) -> dict[str, tuple]:
    """The same comparison for a BindCraft 2 campaign folder.

    Read here from the tables themselves, with BindCraft's own rule for passing
    (`outcome` is exactly `passed`), rather than through the reader: a check
    that went through the reader would agree with it by construction, which is
    the one thing it must not do. It is what lets a campaign of any modality
    say whether its numbers were read correctly, instead of being called
    unverified on the grounds that nobody has looked.
    """
    # The folder and table names come from the reader, so a rename cannot leave
    # this checking files that are no longer there and reporting None for both.
    from .readers import bindcraft2 as bc2

    root = Path(run.source_dir)

    def stage(folder: str):
        here = _child(root, folder)
        return _read_csv(_child(here, bc2.TABLES[folder]) if here is not None else None)

    refolded = stage(bc2.REFOLDED_DIR)
    ranked = stage(bc2.RANKED_DIR)
    traj = stage(bc2.TRAJECTORY_DIR)

    passed = rejected = families = None
    if refolded is not None and "outcome" in refolded.columns:
        outcome = refolded["outcome"].fillna("").astype(str).str.strip()
        is_passed = outcome == "passed"
        passed = int(is_passed.sum())
        # A blank outcome is a row that was not scored, which is neither.
        rejected = int((~is_passed & (outcome != "")).sum())
        if "hash" in refolded.columns:
            families = int(refolded.loc[is_passed, "hash"].astype(str).nunique())

    n_acc = sum(d["status"] == "accepted" for d in exp)
    n_rej = sum(d["status"] == "rejected" for d in exp)
    n_ranked = sum(d.get("rank") is not None for d in exp)
    funnel = {r["field"]: r["value"] for r in run.funnel.stage_rows()}
    return {
        "accepted designs": (passed, n_acc, "!_Refolded.csv rows with outcome `passed` against exported accepted designs"),
        "rejected designs": (rejected, n_rej, "!_Refolded.csv rows scored and not passed against exported rejected designs"),
        "sequence rows": (len(refolded) if refolded is not None else None, funnel.get("mpnn_rows"),
                          "!_Refolded.csv rows against the funnel"),
        "trajectory rows": (len(traj) if traj is not None else None, funnel.get("trajectories_attempted"),
                            "!_Trajectories.csv rows against the funnel"),
        "families": (families, len(run.accepted_backbones()),
                     "distinct `hash` among passed rows against the families"),
        "ranked designs": (len(ranked) if ranked is not None else None, n_ranked,
                           "!_Ranked.csv rows against exported designs carrying a rank"),
    }


def check_counts(project: Project, camp_designs: Sequence[Mapping[str, Any]]) -> tuple[list[str], dict[str, Any]]:
    """Counts read independently from the CSVs and the structure folders against what was exported. Returns (mismatches, table).

    Both BindCraft layouts are compared. A run read by any other reader has an empty table entry and is counted neither
    way. A quantity a layout does not record reads as None on both sides and so does not count as a mismatch."""
    mismatches: list[str] = []
    table: dict[str, Any] = {}
    for rid, run in project.runs.items():
        exp = [d for d in camp_designs if d["run"] == rid]
        if run.reader == "bindcraft1":
            rows = _counts_bindcraft1(run, exp)
        elif run.reader == "bindcraft2":
            rows = _counts_bindcraft2(run, exp)
        else:
            table[rid] = {}
            continue
        table[rid] = {k: {"csv": a, "exported": b, "match": a == b} for k, (a, b, _) in rows.items()}
        for k, (a, b, why) in rows.items():
            if a != b:
                mismatches.append(f"{rid}: {k}: {a} in the files, {b} exported ({why})")
    return mismatches, table


def metric_agreement(project: Project) -> dict[str, Any]:
    """Recomputed interface metrics against BindCraft's own columns for every accepted design (docs/interface_validation.md).

    Keyed by the page metric (``dsasa``) or, when the page has none, the BindCraft stem; compared in the page scale (percent
    columns divided by 100) against the column of the AlphaFold2 model whose pose is in the file."""
    from scipy.stats import spearmanr

    out: dict[str, Any] = {}
    stems = list(next(iter(project.analyses.values())).vs_bindcraft.keys()) if project.analyses else []
    acc = [(k, a) for k, a in project.analyses.items() if project.runs[k[0]].design(k[1]).status == "accepted"]
    for stem in stems:
        if stem == "InterfaceResidues":
            continue
        x, y = [], []
        scale = 0.01 if stem.endswith("%") else 1.0
        for _, a in acc:
            v = a.vs_bindcraft.get(stem)
            if v and v.get("bindcraft") is not None and v.get("recomputed") is not None and v["recomputed"] == v["recomputed"]:
                x.append(float(v["recomputed"]) * scale)
                y.append(float(v["bindcraft"]) * scale)
        if len(x) < 3:
            continue
        xa, ya = np.array(x), np.array(y)
        const = xa.std() == 0 or ya.std() == 0
        err = np.abs(xa - ya)
        name = model.PAGE_METRIC.get(stem) or stem
        out[name] = {
            "stem": stem, "n": len(x), "mean_abs_diff": _r(err.mean(), 4), "max_abs_diff": _r(err.max(), 4),
            "mean_diff": _r((xa - ya).mean(), 4), "r": None if const else _r(np.corrcoef(xa, ya)[0, 1], 3),
            "rho": None if const else _r(spearmanr(xa, ya)[0], 3), "mean_bindcraft": _r(ya.mean(), 4),
            "same_definition": bool(next(iter(project.analyses.values())).vs_bindcraft[stem].get("same_definition")),
            "scale": "percent columns divided by 100" if scale != 1.0 else "CSV scale",
        }
    # the contact residues BindCraft lists against the recomputed 4.0 A set (exact for _model2 files)
    ex, tot, jac = 0, 0, []
    for _, a in acc:
        v = a.vs_bindcraft.get("InterfaceResidues")
        if v:
            tot += 1
            ex += 1 if v.get("exact") else 0
            jac.append(v.get("jaccard") or 0.0)
    if tot:
        out["InterfaceResidues"] = {"stem": "InterfaceResidues", "n": tot, "exact": ex, "mean_jaccard": _r(float(np.mean(jac)), 3),
                                    "note": "BindCraft writes the column from model 2 whatever file is saved; exact for _model2 files"}
    return out


# ------------------------------------------------------------------------------------------------------------
# manifest
# ------------------------------------------------------------------------------------------------------------
_ABS_PATH = re.compile(r"(?<![\w:./])(?:/[\w.@+~-]+){2,}/?")


def portable_paths(obj: Any) -> Any:
    """A copy of ``obj`` in which every absolute path inside a string is cut to its last component. A bundle is meant to be handed
    on, and the paths of the machine that made it (a user name, a cluster's directories: the settings file of a BindCraft run records
    where it ran) say nothing a reader needs; the file name is what identifies a source. The original settings stay in the run folder."""
    if isinstance(obj, str):
        return _ABS_PATH.sub(lambda m: m.group(0).rstrip("/").rsplit("/", 1)[-1], obj)
    if isinstance(obj, dict):
        return {k: portable_paths(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [portable_paths(v) for v in obj]
    return obj


def _git_commit() -> Optional[str]:
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[2], capture_output=True, text=True, timeout=10)
        return r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else None
    except (OSError, subprocess.SubprocessError):
        return None


def _libraries() -> dict[str, str]:
    from importlib import metadata

    out = {}
    for name in ("numpy", "scipy", "pandas", "scikit-learn", "biopython", "gemmi", "tmtools", "Jinja2", "openpyxl"):
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    return out


def _sources(project: Project) -> list[dict[str, Any]]:
    """Every input file with its hash. Structure folders are summarised, one entry per run and role."""
    spec = project.spec
    out: list[dict[str, Any]] = []
    for rid, run in project.runs.items():
        groups: dict[str, list] = {}
        for f in run.files:
            if f.role in ("accepted_structure", "rejected_structure"):
                groups.setdefault(f.role, []).append(f)
            else:
                out.append({"role": f.role, "path": f.path, "sha256": f.sha256, "size": f.size, "run": rid})
        for role, fs in groups.items():
            fs.sort(key=lambda f: f.path)
            folder = str(Path(fs[0].path).parent)
            out.append({"role": role, "path": folder, "sha256": sha256_json([[Path(f.path).name, f.sha256] for f in fs]),
                        "size": sum(f.size for f in fs), "run": rid, "count": len(fs)})
    for role, p in (("target_structure", spec.target_pdb), ("reference_structure", spec.reference_cif)):
        if p is None:
            continue                         # no reference structure was given: there is no such input, and the manifest does not invent one
        out.append({"role": role, "path": str(p), "sha256": sha256_file(p), "size": Path(p).stat().st_size})
    return out


def _reference_entries(project: Project, include_consurf: bool) -> list[dict[str, Any]]:
    """The public reference data the build used, from data/reference/manifest.json (only accessions and ids were sent). Without an
    accession no record was read, so the list is empty."""
    if not project.frame.accession:
        return []
    try:
        entries = rd.manifest().get("entries", {})
    except Exception:  # noqa: BLE001
        return []
    wanted = {w for w in (project.frame.accession, str(project.frame.meta["reference"]["pdb"]).lower() if project.frame.reference_given else "") if w}
    for o in project.homolog.get("orth", []) + project.homolog.get("par", []):
        wanted.add(o["acc"])
    out = []
    for rel, e in sorted(entries.items()):
        if not e.get("file"):
            continue
        low = rel.lower()
        is_consurf = low.startswith("consurf")
        if is_consurf and not include_consurf:
            continue
        if not any(w.lower() in low for w in wanted):
            continue
        host = (e.get("url") or "").split("/")[2] if e.get("url") else ""
        src = {"rest.uniprot.org": "UniProt", "www.ebi.ac.uk": "PDBe SIFTS", "alphafold.ebi.ac.uk": "AlphaFold DB",
               "consurfdb.tau.ac.il": "ConSurf-DB"}.get(host, host or "reference")
        item = {"id": rel, "source": src, "url": e.get("url", ""), "retrieved": e.get("retrieved"), "sha256": e.get("sha256"),
                "bytes": int(e.get("bytes") or 0)}
        item = {k: v for k, v in item.items() if v not in (None, "")}
        if is_consurf:
            item["optional"] = True
            item["licence"] = "academic, non-commercial use"
        out.append(item)
    return out


def _licences(include_consurf: bool, frame: TF.TargetFrame) -> list[dict[str, Any]]:
    out = [
        {"item": "BindCraft output", "status": "included", "terms": "the user's own run output; BindCraft is released under MIT",
         "url": "https://github.com/martinpacesa/BindCraft"},
        {"item": "3Dmol.js", "status": "included", "terms": "BSD-3-Clause, copy in assets/3Dmol-LICENSE.txt", "url": "https://3dmol.csb.pitt.edu"},
    ]
    if frame.accession:     # the UniProt row appears only when records were read for an accession
        out.append({"item": "UniProt sequences and annotation", "status": "included", "terms": "CC BY 4.0", "url": "https://www.uniprot.org"})
    ref = frame.meta["reference"]["pdb"]
    if is_pdb_id(ref):
        out.append({"item": f"PDB {ref}", "status": "included", "terms": "CC0 (wwPDB)", "url": f"https://www.rcsb.org/structure/{ref}"})
    else:
        out.append({"item": f"Reference structure {ref}", "status": "included",
                    "terms": "supplied with the project; its terms are the user's to check before the bundle is shared"})
    out.append({"item": "ConSurf-DB", "status": "optional-on" if include_consurf else "optional-off",
                "terms": "academic, non-commercial use; off by default", "url": "https://consurfdb.tau.ac.il"})
    return out


# ------------------------------------------------------------------------------------------------------------
# highlights
# ------------------------------------------------------------------------------------------------------------
def _join(items: Sequence[str]) -> str:
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def _ranges_text(rs: Sequence[Sequence[int]]) -> str:
    return " and ".join(f"{a}-{b}" for a, b in rs)


def _shared_sites_text(camp: Mapping[str, Any], rids: Sequence[str]) -> str:
    """Whether the runs of a project bind the same site, from the pooled site clusters, or "" when it cannot be said.

    Families are clustered by the target residues they contact, so a cluster
    holding families of more than one run means those runs land on the same
    patch of the target. Said only for a project of several runs, and only when
    the pooled block actually carries site clusters.
    """
    pooled = ((camp.get("diversity") or {}).get("pooled") or {})
    clusters = ((pooled.get("clusters") or {}).get(families.SITE_BASIS) or
                (pooled.get("clusters") or {}).get("site") or [])
    if len(rids) < 2 or not clusters:
        return ""
    run_of = {d["fam"]: d["run"] for d in camp.get("designs", []) if d.get("status") == "accepted" and d.get("fam")}
    shared, alone = [], []
    for c in clusters:
        runs = {run_of.get(f) for f in (c.get("fams") or [])}
        runs.discard(None)
        (shared if len(runs) > 1 else alone).append((c, runs))
    n = len(clusters)
    if not shared:
        where = "; ".join(f"{_join(sorted(r))} in {len(c.get('fams') or [])}" for c, r in alone[:4])
        return (f"The runs do not share a binding site: each of the {n} site "
                f"{'group' if n == 1 else 'groups'} holds families of one run only ({where}). "
                "They are designed against the same target but land on different patches of it, "
                "so a decision in one says nothing about the other.")
    biggest = max(shared, key=lambda pair: len(pair[0].get("fams") or []))
    cov = sum(len(c.get("fams") or []) for c, _ in shared)
    total = sum(len(c.get("fams") or []) for c in clusters)
    if len(shared) == n and n == 1:
        lead = (f"Every family of {_join(sorted(biggest[1]))} falls in one site group: the runs converge on the same "
                "patch of the target")
    else:
        lead = (f"{len(shared)} of the {n} site groups hold families of more than one run "
                f"({cov} of {total} families), the largest shared by {_join(sorted(biggest[1]))}")
    return (lead + ". Binding site is compared across runs even when nothing else can be: it is measured on the "
            "target's residues, not on the binder, so a difference of modality or length does not touch it.")


def compute_highlights(project: Project, profile: Mapping[str, Any], camp: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Four to eight data-derived statements for the overview page: every number is counted from the exported data.

    Each item: ``id``, ``text``, ``severity`` (info, warn or bad) and an optional ``link`` ({page, run, design, filter})."""
    out: list[dict[str, Any]] = []
    designs = camp["designs"]
    acc = [d for d in designs if d["status"] == "accepted"]
    rids = [r["id"] for r in camp["runs"]]
    runs = {r["id"]: r for r in camp["runs"]}
    frame = project.frame
    chk = (profile.get("checks") or {}).get("context") or {}
    by_run = {rid: [d for d in acc if d["run"] == rid] for rid in rids}
    lab = lambda hs: [frame.label(h) for h in hs]  # noqa: E731

    # 0. the context checks did not run: say so before anything that could be read as "nothing found"
    off = profile.get("caps_off") or {}
    if "olig" in off and "rest_ctx" in off and {off["olig"]["why"], off["rest_ctx"]["why"]} <= {"no_reference", "reference_lacks"}:
        what = _join([CONTEXT_FLAGS[f] for f in CONTEXT_FLAGS if f in off and off[f]["why"] in ("no_reference", "reference_lacks")])
        if not frame.reference_given:
            text = (f"No reference structure was given, so no design was checked against {what}. A design without a finding there has not been "
                    "found clear of them; it was not tested. Add the full-length structure (--reference-cif) to run these checks.")
        else:
            text = (f"The reference structure ({frame.meta['reference']['pdb']}) holds nothing outside the construct and no other subunit, so there was nothing "
                    f"to test a design against: {what} are off, and a design without a finding there has not been found clear of them.")
        out.append({"id": "context_off", "text": text, "severity": "warn"})
    # 1. overlap with the rest of the protein
    if all(d["m"].get("rest_n3") is not None for d in acc) and acc:
        thr = (chk.get("rest") or {}).get("fail_n3", 10)
        bad = [d for d in acc if d["m"]["rest_n3"] >= thr]
        clash = [d for d in acc if (project.poses[(d["run"], d["name"])]["diag"].get("rest_clash") or 0) > 0]
        per = [f"{sum(1 for d in by_run[r] if d['m']['rest_n3'] >= thr)} of {len(by_run[r])} in {r}" for r in rids if by_run[r]]
        if bad:
            text = (f"{len(bad)} of {len(acc)} accepted designs fail the rest-of-protein check ({_join(per)}): {thr} or more binder atoms sit within "
                    f"3 \u00c5 of residues outside the construct that BindCraft never saw, as the reference structure places them. {len(clash)} of {len(acc)} have at least one atom closer than 2.2 \u00c5.")
        else:
            text = (f"None of the {len(acc)} accepted designs fails the rest-of-protein check (fewer than {thr} binder atoms within 3 \u00c5 of residues "
                    f"outside the construct, in the reference structure); {len(clash)} have at least one atom closer than 2.2 \u00c5.")
        worst = max(acc, key=lambda d: d["m"]["rest_n3"])
        sev = "bad" if bad and len(bad) / len(acc) > 0.25 else ("warn" if bad else "info")
        out.append({"id": "rest_overlap", "text": text, "severity": sev,
                    "link": {"page": "inspector", "design": worst["id"], "run": worst["run"],
                             "filter": f"largest overlap: {worst['m']['rest_n3']} atoms within 3 \u00c5 of the rest of the protein"}})
        notes = _notes_clause(project, profile)
        if notes:
            out.append({"id": "run_notes", "text": notes, "severity": "info", "link": {"page": "board", "filter": "families named in the run notes"}})
    # 2. the other subunit
    if acc and all(d["m"].get("olig_min") is not None for d in acc):
        thr = (chk.get("partner") or {}).get("fail_n3", 10)
        over = [d for d in acc if d["m"]["olig_n3"] > 0]
        fail = [d for d in acc if d["m"]["olig_n3"] >= thr]
        near = min(acc, key=lambda d: d["m"]["olig_min"])
        med = statistics.median(d["m"]["olig_min"] for d in acc)
        if not over:
            text = (f"No accepted design reaches the other subunit in the reference assembly: the closest binder atom is {near['m']['olig_min']:.1f} Å from it "
                    f"(median {med:.0f} Å), so all {len(acc)} pass the dimer-overlap check.")
            sev = "info"
        else:
            text = (f"{len(over)} of {len(acc)} accepted designs put binder atoms within 3 Å of the other subunit in the reference assembly ({len(fail)} fail the "
                    f"dimer-overlap check at {thr} or more atoms); the closest is {near['m']['olig_min']:.1f} Å.")
            sev = "bad" if fail else "warn"
        out.append({"id": "partner_overlap", "text": text, "severity": sev,
                    "link": {"page": "inspector", "design": near["id"], "run": near["run"], "filter": "closest approach to the other subunit"}})
    # 3. hotspots per run
    hs_parts, low_share, any_none = [], 1.0, False
    for rid in rids:
        ds = by_run[rid]
        hs = runs[rid]["settings"]["hotspots"]
        if not ds or not hs:
            continue
        touched = {h: sum(1 for d in ds if h in d["hot"]) for h in hs}
        allh = sum(1 for d in ds if len(d["hot"]) == len(hs))
        none = sum(1 for d in ds if not d["hot"])
        order = sorted(hs, key=lambda h: -touched[h])
        first = f"{frame.label(order[0])} is touched by {touched[order[0]]} of {len(ds)}"
        rest = [f"{frame.label(h)} by {touched[h]}" for h in order[1:]]
        hs_parts.append(f"{rid} ({_join(lab(hs))}): {first}" + (f", {_join(rest)}" if rest else "")
                        + f"; {allh} {'touches' if allh == 1 else 'touch'} all {len(hs)}"
                        + (f", {none} {'touches' if none == 1 else 'touch'} none" if none else ""))
        low_share = min(low_share, min(touched.values()) / len(ds))
    none_rec = [rid for rid in rids if not runs[rid]["settings"]["hotspots"]]
    if hs_parts or none_rec:
        text = ("Hotspots touched by the accepted designs. " + ". ".join(hs_parts) + "." if hs_parts else "")
        if none_rec:
            text += (" " if text else "") + (f"No hotspots recorded for {_join(none_rec)}: its designs have no hotspots touched and no intended epitope" if len(none_rec) == 1 else
                                              f"No hotspots recorded for {_join(none_rec)}: their designs have no hotspots touched and no intended epitope")
            text += ("; the binding-site view (site A1) is built from the residues the designs contact most." if len(none_rec) == len(rids) else ".")
        out.append({"id": "hotspots", "text": text, "severity": "warn" if (hs_parts and low_share < 0.25) else "info",
                    "link": {"page": "map", "filter": "hotspot coverage per run"}})
    # 4. what differs between the runs
    if len(rids) > 1:
        diffs = []
        hot_txt = " against ".join(_join(lab(runs[r]["settings"]["hotspots"])) or "none recorded" for r in rids)
        if len({tuple(runs[r]["settings"]["hotspots"]) for r in rids}) > 1:
            diffs.append(f"hotspots {hot_txt}")
        nres = [runs[r]["target"]["n_resi"] for r in rids]
        if len(set(nres)) > 1:
            removed = {r: runs[r]["target"]["removed"] for r in rids if runs[r]["target"]["removed"]}
            extra = "".join(f" ({r} lacks {_ranges_text(rem)})" for r, rem in removed.items())
            diffs.append(f"construct {' against '.join(str(n) for n in nres)} residues{extra}")
        lens = [runs[r]["settings"]["length_range"] for r in rids]
        if len({tuple(x) if x else None for x in lens}) > 1:
            diffs.append("binder length " + " against ".join(f"{x[0]}-{x[1]}" if x else "unknown" for x in lens))
        if diffs:
            out.append({"id": "run_differences", "severity": "info", "link": {"page": "board", "filter": "per-run diversity"},
                        "text": f"The runs differ in {_join(diffs)}. Diversity and the need for more designs are therefore estimated per run and never pooled."})
        # 4b. where the runs land on the target. This is the one comparison that
        # survives a difference of modality or length: a footprint is on the
        # target, so it does not matter what bound it. Without this the review
        # lists what cannot be compared across runs and says nothing about what
        # can -- and two campaigns converging on one epitope, or missing each
        # other entirely, is the thing a project review exists to show.
        site_text = _shared_sites_text(camp, rids)
        if site_text:
            out.append({"id": "shared_sites", "severity": "info", "link": {"page": "map", "filter": "sites"},
                        "text": site_text})
    # 5. families and coverage
    div = camp["diversity"]["runs"]
    parts, sev = [], "info"
    basis = camp["diversity"]["basis_default"]
    fold_off: list[str] = []
    why_off: dict[str, str] = {}
    for rid in rids:
        blk = div.get(rid) or {}
        off = (blk.get("fold") or {}).get("status") == "off"
        if off:
            fold_off.append(rid)
            why_off[rid] = str((blk.get("fold") or {}).get("why") or "")
        adv = (blk.get("advice") or {}).get(blk.get("basis_default") if off else basis)
        n = blk.get("n_families")
        if n is None:
            continue
        if adv:
            parts.append(f"{rid} has {n} (coverage {adv['coverage']:.2f}, level {adv['level']})")
            if adv["level"] == "keep":
                sev = "warn"
        else:
            parts.append(f"{rid} has {n}")
    if parts:
        low = [rid for rid in rids if (div.get(rid) or {}).get("n_families", 0) < families.LOW_CONFIDENCE_BELOW]
        unk = [rid for rid in rids if runs[rid]["n_traj"] is None]
        tail = []
        if low:
            tail.append(f"The estimates for {_join(low)} are low-confidence (fewer than {families.LOW_CONFIDENCE_BELOW} families).")
        if unk:
            tail.append(f"Attempted trajectories are unknown for {_join(unk)}, so yield is per passed trajectory.")
        # Why fold is off differs per run -- a shared scaffold, binders too short
        # for TM-score to be calibrated, a cyclic backbone -- so the reason is
        # taken from the block that made the decision. Saying "its designs share
        # a scaffold" about a peptide campaign is simply untrue.
        reasons = sorted({why_off[rid] for rid in fold_off if why_off.get(rid)})
        because = "; and ".join(reasons) if reasons else "fold distance does not apply to these designs"
        if not fold_off:
            how = "measured on fold and site together"
        elif len(fold_off) == len(rids):
            how = f"measured on the binding site alone: {because}"
        else:
            how = f"measured on fold and site together, except for {_join(fold_off)} (binding site alone: {because})"
        out.append({"id": "families", "severity": sev, "link": {"page": "board", "filter": "diversity panel"},
                    "text": f"Families, not sequences: {_join(parts)}, {how}." + "".join(" " + t for t in tail)})
    # 5b. what a scaffolded run's designs actually bind with (docs/paratope.md)
    for rid in rids:
        text = _paratope_highlight(rid, runs[rid])
        if text:
            out.append({"id": f"paratope_{rid}" if len(rids) > 1 else "paratope", "severity": "info", "text": text,
                        "link": {"page": "board", "filter": "paratope panel"}})
    # 5c. which end of a family can carry a tag (docs/termini.md)
    for rid in rids:
        text = _terminus_highlight(runs[rid], by_run[rid])
        if text:
            out.append({"id": f"termini_{rid}" if len(rids) > 1 else "termini", "severity": "info", "text": text,
                        "link": {"page": "board", "filter": "tag panel"}})
    # 6. near-miss rescues
    rej = [d for d in designs if d["status"] == "rejected"]
    if rej:
        near = [d for d in rej if (d.get("nm") or {}).get("near")]
        stems = Counter((d["nm"]["stems"] or ["?"])[0] for d in near)
        where = [r for r in rids if any(d["run"] == r for d in rej)]
        scope = f"all in {where[0]}" if len(where) == 1 else ", ".join(f"{sum(1 for d in rej if d['run'] == r)} in {r}" for r in where)
        detail = _join([f"{n} on {FAILURE_SHORT.get(s, s)}" for s, n in stems.most_common()])
        # a threshold the campaign recorded on a metric Binderscope does not carry cannot be re-checked: a near miss is a miss of one
        # filter among those that can be checked
        blind = sorted({u.name for rid in where for u in config.filters_for_run(project.runs[rid]).unevaluated if u.mandatory})
        caveat = (f" Thresholds on {_join(blind)} are the campaign's own and were not re-checked here, so a design may also miss one of them."
                  if blind and near else "")
        out.append({"id": "near_misses", "severity": "info", "link": {"page": "query", "filter": "near-miss rescue"},
                    "text": f"{len(near)} of the {len(rej)} rejected designs that have a structure ({scope}) are near-miss rescues: a single filter "
                            f"missed by at most {int(model.NEAR_MISS_MAX_REL * 100)}% of its threshold" + (f" ({detail})." if detail else ".") + caveat})
    # 7. specificity at the hotspot residues
    spec = _specificity_highlight(project, profile, camp)
    if spec:
        out.append(spec)
    return out


def _notes_clause(project: Project, profile: Mapping[str, Any]) -> str:
    """A sentence on the backbones a run's own notes list, when the checks say something about them (data from the note file)."""
    chk = (profile.get("checks") or {}).get("context") or {}
    sentences = []
    for rid, run in project.runs.items():
        for nf in (run.annotations or {}).get("note_files", []) or []:
            ids = [i for i in nf.get("ids", []) if run.backbone(i) is not None]
            if len(ids) < 2:
                continue
            worst = {}
            for bb in ids:
                lv = "pass"
                for d in run.accepted():
                    if d.backbone_id != bb:
                        continue
                    pose = project.poses.get((rid, d.id))
                    if not pose:
                        continue
                    v = TF.verdicts(pose["metrics"], pose["diag"], chk)
                    for part in ("partner", "rest"):
                        if part in v and v[part]["level"] == "fail":
                            lv = "fail"
                    if lv != "fail" and (v.get("cut", {}).get("level") == "warn" or v.get("rest", {}).get("level") == "warn"):
                        lv = "warn"
                worst[bb] = lv
            n_fail = sum(1 for v in worst.values() if v == "fail")
            n_warn = sum(1 for v in worst.values() if v == "warn")
            if n_fail or n_warn:
                sentences.append(f"The notes file of {rid} lists {len(ids)} backbones; the context checks fail "
                                 f"{n_fail} of them (rest-of-protein overlap)" + (f" and warn on {n_warn} more." if n_warn else "."))
    return " ".join(sentences)


def _specificity_highlight(project: Project, profile: Mapping[str, Any], camp: Mapping[str, Any]) -> Optional[dict[str, Any]]:
    hom = project.homolog
    if not hom.get("orth"):
        return None
    sp = (profile.get("specificity") or {})
    hot = sorted({h for r in camp["runs"] for h in r["settings"]["hotspots"]})
    if not hot:
        return None
    idx = {r: i for i, r in enumerate(hom["resis"])}
    hot = [h for h in hot if h in idx]
    names = {e["key"]: (e.get("phrase") or (e.get("label") or e["key"]).lower()) for e in sp.get("species", [])}
    same, differ = [], {}
    for o in hom["orth"]:
        key = o.get("key") or o["acc"]
        diff = {h: o["map"][idx[h]] for h in hot if o["map"][idx[h]] != hom["seq"][idx[h]]}
        (differ.__setitem__(key, diff) if diff else same.append(key))
    label = lambda k: names.get(k, k)  # noqa: E731
    parts = []
    if same:
        parts.append(f"{_join([label(k) for k in same])} keep all {len(hot)}")
    if differ:
        by_res: dict[int, dict[str, list[str]]] = {}
        for k, diff in differ.items():
            for h, aa in diff.items():
                by_res.setdefault(h, {}).setdefault(aa, []).append(label(k))
        for h in sorted(by_res, key=lambda h: -sum(len(v) for v in by_res[h].values())):
            var = _join([f"{aa} in {_join(ks)}" for aa, ks in sorted(by_res[h].items())])
            parts.append(f"{frame_label(project, h)} varies ({var})")
    foc = ((sp.get("paralogs") or {}).get("focus") or {}).get("key")
    par = next((p for p in hom.get("par", []) if p.get("key") == foc or p.get("gene") == foc), None)
    if par:
        n_same = sum(1 for h in hot if par["map"][idx[h]] == hom["seq"][idx[h]])
        parts.append(f"{foc} matches {n_same} of the {len(hot)}")
    headline = set(sp.get("headline") or [])
    sev = "warn" if headline & set(differ) else "info"
    return {"id": "specificity", "severity": sev, "link": {"page": "map", "filter": "species and paralog specificity"},
            "text": f"Sequence screen at the {len(hot)} hotspot residues of the runs: " + "; ".join(parts) + "."}


def _terminus_highlight(run: Mapping[str, Any], accepted: Sequence[Mapping[str, Any]]) -> Optional[str]:
    """One statement for a run: of its families, how many can take a His6 and a Twin-Strep at which end, and what BindCraft was asked about the ends.
    Counted from the exported blocks, in families and not in sequences. None when the question has no answer for the run or nothing was judged."""
    rt = run.get("term")
    if not rt or not (rt.get("format") or {}).get("applies", True):
        return None
    fams: dict[str, list[Mapping[str, Any]]] = {}
    for d in accepted:
        if (d.get("term") or {}).get("pick"):
            fams.setdefault(d["bb"], []).append(d)
    if not fams:
        return None
    pk = {"his6": {}, "strep": {}}
    for bb, ds in fams.items():
        for k in pk:
            pk[k][bb] = termini.family_pick([d["term"]["pick"][k] for d in ds])
    n = len(fams)
    cnt = {k: Counter(v.values()) for k, v in pk.items()}
    fmt = rt["format"]
    word = "family" if n == 1 else "families"

    def line(k: str, label: str) -> str:
        c = cnt[k]
        return f"{label}: {c['C']} at the C-terminus, {c['N']} at the N-terminus, {c['both']} at either end, {c['none']} at neither"
    text = (f"Where a tag can go ({n} {word} of {run['id']}, predicted from one conformation). {line('his6', 'His6')}. {line('strep', 'Twin-Strep')}.")
    if fmt.get("convention"):
        conv = fmt["convention"]
        text += (f" A {fmt['kind']} is conventionally tagged at the {conv}-terminus; the geometry alone puts {cnt['his6'][conv]} of {n} {word} there "
                 f"for a His6 and {cnt['strep'][conv]} of {n} for a Twin-Strep.")
    df = rt.get("designed_for") or {}
    if df.get("known"):
        if df.get("accessible"):
            text += " BindCraft was asked to angle both ends away from the target, so this confirms that setting more than it discovers anything."
        elif df.get("buried"):
            text += " BindCraft was asked to bury the ends, the opposite of what a tag wants, so tight ends are expected."
        else:
            text += " BindCraft was not asked to keep the ends free, so the spread is unconstrained."
    return text


def _paratope_highlight(rid: str, run: Mapping[str, Any]) -> Optional[str]:
    """One sentence on the paratope of a scaffolded run, counted from its paratope block (None for a run without one)."""
    blk, rej = run.get("paratope"), run.get("paratope_rejected")
    if not blk or not blk.get("share"):
        return None
    segs = [s["name"] for s in blk["segments"]]
    shares = "; ".join(f"{n} {blk['share'][n]:.1f}%" for n in segs if n in blk["share"])
    full = next((u for u in blk["loop_use"] if u["loops"] == segs), None)
    n_full = full["n"] if full else 0
    where = f"{shares}; framework {blk['framework_share']:.1f}%"
    rfull = next((u["share"] for u in (rej or {}).get("loop_use", []) if u["loops"] == segs), 0.0)
    contrast = (f" Among the rejected designs framework carries {rej['framework_share']:.1f}% and {rfull:.1f}% use all {len(segs)}." if rej else "")
    return (f"{rid} is built on one scaffold ({blk['scaffold'] or 'recorded in the campaign'}), so fold diversity is switched off. Interface residues of its "
            f"{blk['n']} accepted designs fall in {where}; {n_full} of {blk['n']} use all {len(segs)} designed segments.{contrast}")


def frame_label(project: Project, h: int) -> str:
    return project.frame.label(h)


# ------------------------------------------------------------------------------------------------------------
# assembly
# ------------------------------------------------------------------------------------------------------------
def build_blobs(project: Project, *, created: Optional[str] = None, validate: bool = True, strict: bool = True) -> Blobs:
    """Assemble and validate the blobs of a Project. ``strict`` raises ExportError on any schema or semantic error."""
    warnings: list[str] = []
    mismatches: list[str] = []
    frame, runs = project.frame, project.runs
    created = created or utc_now()
    today = created[:10]
    profile0 = completed_profile(project)

    para = run_paratopes(project)
    run_items = [run_item(project, rid, para[rid]) for rid in runs]
    site = site_residues(runs, project.poses)           # where the designs sit: site A1 of a project without any hotspot
    target = clean(frame.target_block(list(runs), name=(profile0.get("target") or {}).get("name"), site_residues=site))
    _name_segments(target, profile0.get("frame"))
    designs, inspected = _design_records(project, profile0, mismatches, para)
    warnings += placement_warnings(project, designs)
    ins_designs, order, structures = _inspector_designs(project, profile0, designs, inspected)

    fold = project.families["fold"]
    parts = {"designs": designs, "target": target, "homolog": project.homolog, "fold": fold, "runs": run_items,
             "rest": target.get("rest")}
    profile = build_profile(project, parts, warnings)
    generated = {"tool": "binderscope", "version": __version__, "date": today}
    modes = {k: dict(v) for k, v in project.modes["modes"].items()}
    camp: dict[str, Any] = {
        "schema": "binderscope.campaign/1", "profile": profile, "target": target, "designs": designs, "modes": modes,
        "run": run_summary_block(project, run_items, profile), "runs": run_items, "fold": fold,
        "diversity": project.families["diversity"], "generated": generated,
    }
    camp["highlights"] = compute_highlights(project, profile, camp)

    ib = clean(frame.inspector_block(list(runs), site_residues=site))
    # The paratope blocks travel with the inspector's own runs too: it labels a design's
    # segments from the run's scaffold layout and never reads the campaign blob. Only the
    # layout is needed here, not the shares, so paratope_rejected stays where it is.
    _para_by_run = {r["id"]: r.get("paratope") for r in camp["runs"] if r.get("paratope")}
    for _r in ib.get("runs") or []:
        _p = _para_by_run.get(_r.get("id"))
        if _p:
            _r["paratope"] = _p
    _term_by_run = {r["id"]: r.get("term") for r in camp["runs"] if r.get("term")}
    for _r in ib.get("runs") or []:
        if _term_by_run.get(_r.get("id")):
            _r["term"] = _term_by_run[_r["id"]]
    olig = ib.get("olig")
    if olig and "buried" in olig:
        # the whole interface on the reference chain (olig.buried holds only the part inside the union construct): the inspector's frame strip draws it
        olig = {**olig, "buried_all": {str(u): round(a, 1) for u, a in sorted(frame.partner_buried().items())}}
    inspect_blob: dict[str, Any] = {
        "schema": "binderscope.inspect/1", "profile": profile, "target_pdb": ib["target_pdb"], "target_seq": ib["target_seq"],
        "target_resi": ib["target_resi"], "t_ss": ib["t_ss"], "designs": ins_designs, "order": order, "hotspots": ib["hotspots"],
        "epitope": ib["epitope"], "glycan": ib["glycan"], "glycan_xyz": ib["glycan_xyz"], "glycan_atoms": ib["glycan_atoms"],
        "membrane": ib["membrane"], "runs": ib["runs"], "olig": olig, "rest": ib["rest"], "generated": generated,
    }
    homolog = json.loads(json.dumps(project.homolog))
    homolog.setdefault("generated", {}).update({"date": today})
    count_mismatches, count_table = check_counts(project, designs)
    agreement = metric_agreement(project)
    manifest = build_manifest(project, profile, camp, inspect_blob, agreement, count_mismatches, mismatches, count_table, created, warnings)
    manifest = portable_paths(manifest)  # file names only: no user name or directory of this machine or of the cluster that ran BindCraft
    camp, inspect_blob, homolog, manifest = (clean(x) for x in (camp, inspect_blob, homolog, manifest))
    profile = camp["profile"]
    blobs = Blobs(camp, inspect_blob, homolog, profile, manifest, structures, warnings)
    if validate:
        validate_blobs(blobs, strict=strict)
    return blobs


def build_manifest(project: Project, profile: dict[str, Any], camp: Mapping[str, Any], ins: Mapping[str, Any],
                   agreement: Mapping[str, Any], count_mismatches: list[str], other_mismatches: list[str], count_table: Mapping[str, Any],
                   created: str, warnings: list[str]) -> dict[str, Any]:
    spec = project.spec
    sources = _sources(project)
    profile_sha = sha256_json(profile)
    bundle_id = sha256_json(sorted(s["sha256"] for s in sources) + [profile_sha])[:16]
    runs = []
    for r in camp["runs"]:
        run = project.runs[r["id"]]
        runs.append({
            "id": r["id"], "label": r["label"], "reader": r["reader"], "experimental": r["experimental"],
            "tool_version": r["tool_version"], "source_dir": str(run.source_dir), "settings": r["settings"],
            "settings_raw": run.settings.settings_raw, "target": r["target"], "filters": r["filters"],
            **({"filter_set": r["filter_set"]} if "filter_set" in r else {}),
            "counts": {"n_designs": r["n_designs"], "funnel": r["funnel"], "failure": r["failure"],
                       "csv": count_table.get(r["id"], {})},
            "epoch": r["epoch"], "warnings": r["warnings"], "notes": r["notes"],
            "annotations": _annotations(run),
        })
    designs = camp["designs"]
    n_acc = sum(d["status"] == "accepted" for d in designs)
    families_n = len({d["bb"] for d in designs if d["status"] == "accepted"})
    tgt = profile["target"]
    all_warn = list(warnings) + list(project.warnings)
    return {
        "schema": "binderscope.manifest/1",
        "tool": {"name": "binderscope", "package": "binder-review", "version": __version__, "git_commit": _git_commit(),
                 "python": platform.python_version(), "platform": platform.platform(), "command": spec.command,
                 "libraries": _libraries()},
        "created": created,
        "bundle": {"id": bundle_id, "name": spec.name, "title": f"{tgt.get('short', spec.name)} · {len(camp['runs'])} BindCraft run{'' if len(camp['runs']) == 1 else 's'}",
                   "pages": {k: f"{k}.html" for k in ("index", "board", "inspector", "map", "query")},
                   "data": {k: f"data/{k}.json" for k in ("campaign", "inspect", "homolog", "profile", "manifest")}, "files": []},
        "profile": {"id": profile["id"], "sha256": profile_sha},
        "target": {"name": tgt["name"], "uniprot": tgt.get("uniprot", ""), "reference_pdb": tgt["reference"]["pdb"], "reference_given": project.frame.reference_given,
                   "reference_class": project.frame.reference_class,
                   "context_check": {k: profile["checks"]["context"][k] for k in ("hard", "hard_setting", "hard_basis")},
                   "organism": tgt.get("organism", "")},
        "sources": sources, "reference": _reference_entries(project, spec.consurf), "runs": runs, "caps": dict(profile["caps"]),
        **({"caps_off": dict(profile["caps_off"])} if profile.get("caps_off") else {}),
        "thresholds": {
            "checks": profile["checks"], "context": profile["checks"]["context"],
            "filters": project_filter_rules(project),
            "diversity": {**project.families["diversity"]["distance_default"], "fold_weight": project.families["diversity"]["fold_weight"],
                          "basis_default": project.families["diversity"]["basis_default"]},
            "near_miss": {"max_failed_stems": model.NEAR_MISS_MAX_STEMS, "max_rel": model.NEAR_MISS_MAX_REL},
            "interface": {"contact": model.CONTACT_CUTOFF, "bindcraft_interface": model.BINDCRAFT_INTERFACE_CUTOFF,
                          "hbond": model.HBOND_MAX_DIST, "salt_bridge": model.SALT_BRIDGE_MAX_DIST, "partner_contact": 3.0,
                          "partner_near": 4.0, "clash": 2.2, "sasa_probe": 1.4, "sasa_points": 960,
                          "place_warn_rmsd": project.frame.config["geometry"]["place_warn_rmsd"]},
            "placement": _placement_thresholds(project),
        },
        "counts": {"designs": len(designs), "accepted": n_acc, "rejected": len(designs) - n_acc, "families": families_n,
                   "runs": len(camp["runs"]), "inspected": len(ins["designs"])},
        "validation": {"counts_match_csv": (not count_mismatches) if any(count_table.values()) else None,    # null: no run could be compared
                       "mismatches": list(count_mismatches) + list(other_mismatches),
                       "metric_agreement": agreement},
        "licences": _licences(spec.consurf, project.frame),
        "privacy": {"network_requests": project.network_requests, "binder_sequences_sent": False, "design_ids_sent": False},
        "warnings": all_warn,
        "notes": [
            "Scores are screens, not verdicts: a passing design is a predicted binder, not a measured one.",
            _placement_note(project),
            *project.frame.meta.get("notes", []),
        ],
        "highlights": camp.get("highlights", []),
    }


def _placement_note(project: Project) -> str:
    """The manifest note on how designs were placed and how far the distances can be trusted. With a reference the distances are estimates
    against one static model; without one nothing outside the construct was measured."""
    fit = {"interface": ", weighted around the binder (the residues near it count most)",
           "core": " over the residues that a trimmed fit keeps"}.get(project.frame.config["geometry"].get("fit"), " over all shared residues")
    if project.frame.reference_given:
        return (f"The placement of each design on {project.frame.meta['reference']['pdb']} is a rigid fit of its own target on the crystal chain{fit}; "
                "distances to the other subunit and to the rest of the protein are estimates to a few angstrom, because the reference is one static model.")
    return (f"Each design's own target is superposed on the construct file{fit}, so that all designs share one frame; the construct is the file given "
            "to BindCraft, not a structure of the full-length protein, and nothing outside it was measured.")


def _placement_thresholds(project: Project) -> dict[str, Any]:
    """manifest.thresholds.placement: how the designs were placed on the reference chain and what the interface fit weighs."""
    geo = {**TF.DEFAULT_CONFIG["geometry"], **project.frame.config["geometry"]}
    return {"mode": getattr(project.spec, "placement", "interface"), "fit": geo["fit"], "iface_full": geo["iface_full"], "iface_decay": geo["iface_decay"],
            "iface_floor": geo["iface_floor"], "iface_cutoff": geo["iface_cutoff"], "iface_min_n": geo["iface_min_n"],
            "near_weight": TF.IFACE_NEAR_WEIGHT, "report_weight": TF.IFACE_REPORT_WEIGHT, "warn_rmsd": geo["place_warn_rmsd"]}


def _annotations(run: RunData) -> dict[str, Any]:
    ann = run.annotations or {}
    out: dict[str, Any] = {}
    # What kind of campaign this is, as the reader worked it out. Carried here
    # because the bundle otherwise only has it buried in a note, and a page (or
    # a person reading the manifest) that wants to say "this is a peptide
    # campaign" should not have to parse prose for it. "mixed" is a real
    # answer: a run whose designs do not agree.
    if ann.get("modality"):
        out["modality"] = ann["modality"]
        if ann.get("modalities"):
            out["modalities"] = dict(ann["modalities"])
    if ann.get("xlsx_highlight"):
        out["xlsx_highlight"] = [{"file": Path(h["file"]).name, "sheet": h.get("sheet"), "colours": h.get("colours"), "n_designs": len(h.get("designs", [])),
                                  "meaning": "not recorded in the file; shown only behind the outcomes switch"} for h in ann["xlsx_highlight"]]
    if ann.get("note_files"):
        out["note_files"] = [{"file": Path(n["path"]).name, "ids": n.get("ids", [])} for n in ann["note_files"]]
    sc = ann.get("scaffold")
    if sc and sc.get("scaffold_edits"):
        out["scaffold"] = {"file": Path(sc["scaffold"]).name if sc.get("scaffold") else None, "scaffold_edits": sc["scaffold_edits"],
                           "retained_fraction": {"threshold": sc.get("threshold"), "mandatory": sc.get("mandatory")},
                           "meaning": "the scaffold positions the campaign was allowed to change; the paratope analysis divides the scaffold by it"}
    return out


# ------------------------------------------------------------------------------------------------------------
# validation and writing
# ------------------------------------------------------------------------------------------------------------
def validate_blobs(blobs: Blobs, *, strict: bool = True) -> Any:
    """Every blob against its JSON Schema plus the semantic and cross rules (scripts/validate_schema.py). Errors raise
    ``ExportError`` (strict); warnings are kept in ``blobs.warnings``."""
    try:
        vs = _load_validator()
    except ImportError as ex:           # jsonschema is a dependency of the package; a source tree without it gets a clear message
        raise ExportError(f"the blobs cannot be validated: {ex}") from ex
    if vs is None:
        blobs.warnings.append("scripts/validate_schema.py not found: the blobs were not validated")
        if strict:
            raise ExportError("scripts/validate_schema.py not found; cannot validate (pass strict=False to skip)")
        return None
    # the validators mutate nothing that stays, but work on copies to be safe about transient helper keys
    rep = vs.check_blobs(blobs.campaign, blobs.inspect, blobs.homolog, "export")
    rep.merge(vs.check_manifest(blobs.manifest), "manifest ")
    if blobs.manifest["caps"] != blobs.campaign["profile"]["caps"]:
        rep.error("manifest.caps", "differs from campaign.profile.caps")
    blobs.validation = rep
    blobs.warnings += [f"validator: {w}" for w in rep.warnings]
    if rep.errors and strict:
        shown = "\n  ".join(rep.errors[:25])
        raise ExportError(f"{len(rep.errors)} schema or consistency errors in the exported blobs:\n  {shown}")
    return rep


def write_export(blobs: Blobs, out_dir: str | Path) -> list[Path]:
    """work/<project>/export: the blobs as JSON, the structure files and the profile (the input of ``bundle.write_bundle``)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "structures").mkdir(exist_ok=True)
    for old in (out / "structures").glob("*.pdb"):
        old.unlink()
    paths = []
    for name, obj in (("campaign", blobs.campaign), ("inspect", blobs.inspect), ("homolog", blobs.homolog), ("profile", blobs.profile),
                      ("manifest", blobs.manifest)):
        p = out / f"{name}.json"
        p.write_text(dumps(obj) + "\n", encoding="utf-8")
        paths.append(p)
    for fn, text in blobs.structures.items():
        p = out / "structures" / fn
        p.write_text(text, encoding="utf-8")
        paths.append(p)
    return paths
