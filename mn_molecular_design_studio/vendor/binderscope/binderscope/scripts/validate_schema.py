#!/usr/bin/env python3
"""Validate Binderscope data against docs/schema and against the rules JSON Schema cannot express.

What it checks
  * every blob against its JSON Schema (draft 2020-12, cross-file $ref through a referencing registry):
    campaign (RAW), inspect (INS), homolog (HOM), manifest, profile, review event;
  * semantic rules per blob (array lengths, index ranges, numbering, run references, capability flags against the data);
  * cross rules between blobs (INS against RAW, HOM against RAW and the profile);
  * a whole bundle folder (data files, hashes in the manifest, the review files).

Command line (use the project venv)
  scripts/validate_schema.py --examples                 profiles and example blobs in docs/schema/examples
  scripts/validate_schema.py --bundle path/to/bundle    a bundle folder (data/*.json, reviews/*.jsonl)
  scripts/validate_schema.py --kind campaign FILE.json  one file against one schema
Exit status 0 = no errors (warnings may remain), 1 = errors.

Library use (tests, the exporter):
  import importlib.util, pathlib, sys
  spec = importlib.util.spec_from_file_location('validate_schema', ROOT / 'scripts' / 'validate_schema.py')
  vs = importlib.util.module_from_spec(spec); sys.modules['validate_schema'] = vs; spec.loader.exec_module(vs)
  rep = vs.check_bundle(bundle_dir)      # Report with .errors and .warnings (lists of strings), .ok
Messages read '<json path>: <what is wrong>'. Errors mean a page would break or show a wrong thing; warnings mean the
data is odd but renders. docs/data_contract.md says which rule comes from which field.
"""
from __future__ import annotations

import argparse
import base64
import binascii
import copy
import functools
import hashlib
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
try:
    from jsonschema import Draft202012Validator
    from referencing import Registry, Resource
    from referencing.jsonschema import DRAFT202012
except ImportError as _ex:     # the exporter loads this file as a module: say what is missing instead of a bare traceback
    _msg = "validate_schema.py needs the jsonschema package, version 4.18 or later (it brings referencing): pip install 'jsonschema>=4.26'"
    if __name__ == "__main__":
        sys.exit(_msg)
    raise ImportError(_msg) from _ex

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA_DIR = ROOT / "docs" / "schema"
EXAMPLE_DIR = SCHEMA_DIR / "examples"

KINDS = ("campaign", "inspect", "homolog", "manifest", "profile", "review_event")

# capability flags that promise a metric column (docs/profile_spec.md section 4)
METRIC_CAPS = (
    "ipsae", "min_ipae", "complex_rmsd", "dg", "sc", "iptm_b", "ipae_norm", "surf_hyd",
    "hotspot_rmsd", "binder_rmsd", "seed_spread", "n_agree", "esm_plddt", "ciplddt",
)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}

TOKEN = re.compile(r"\{([^{}|]+)(?:\|([a-z]+))?\}")


# --------------------------------------------------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------------------------------------------------
@dataclass
class Report:
    name: str = ""
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)  # what was looked at (designs, residues, ...), so tests can see that checks ran

    def count(self, key: str, n: int = 1) -> None:
        self.stats[key] = self.stats.get(key, 0) + n

    def error(self, where: str, msg: str) -> None:
        self.errors.append(f"{where}: {msg}")

    def warn(self, where: str, msg: str) -> None:
        self.warnings.append(f"{where}: {msg}")

    def merge(self, other: "Report", prefix: str = "") -> "Report":
        self.errors += [prefix + e for e in other.errors]
        self.warnings += [prefix + w for w in other.warnings]
        for k, v in other.stats.items():
            self.count(k, v)
        return self

    @property
    def ok(self) -> bool:
        return not self.errors

    def format(self, limit: int = 40, warnings: bool = True) -> str:
        seen = ", ".join(f"{k} {v}" for k, v in sorted(self.stats.items()))
        lines = [f"== {self.name}: {len(self.errors)} error(s), {len(self.warnings)} warning(s)" + (f" [{seen}]" if seen else "")]
        for e in self.errors[:limit]:
            lines.append("  ERROR " + e)
        if len(self.errors) > limit:
            lines.append(f"  ... {len(self.errors) - limit} more errors")
        if warnings:
            for w in self.warnings[:limit]:
                lines.append("  warn  " + w)
            if len(self.warnings) > limit:
                lines.append(f"  ... {len(self.warnings) - limit} more warnings")
        return "\n".join(lines)


# --------------------------------------------------------------------------------------------------------------------
# JSON Schema layer
# --------------------------------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=1)
def registry() -> Registry:
    resources = []
    for p in sorted(SCHEMA_DIR.glob("*.schema.json")):
        doc = json.loads(p.read_text(encoding="utf-8"))
        resources.append((doc["$id"], Resource.from_contents(doc, default_specification=DRAFT202012)))
    return Registry().with_resources(resources)


@functools.lru_cache(maxsize=None)
def validator(kind: str) -> Draft202012Validator:
    schema = json.loads((SCHEMA_DIR / f"{kind}.schema.json").read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema, registry=registry())


def _path(parts) -> str:
    out = "$"
    for p in parts:
        out += f"[{p}]" if isinstance(p, int) else f".{p}"
    return out


def schema_errors(kind: str, instance, limit: int = 40) -> list[str]:
    """First `limit` schema violations of `instance` against the schema of `kind`, as 'path: message' strings."""
    out = []
    for e in validator(kind).iter_errors(instance):
        msg = e.message if len(e.message) <= 220 else e.message[:217] + "..."
        out.append(f"{_path(e.absolute_path)}: {msg}")
        if len(out) >= limit:
            out.append("... (more schema errors not listed)")
            break
    return out


def _schema_into(rep: Report, kind: str, instance, prefix: str = "") -> bool:
    errs = schema_errors(kind, instance)
    for e in errs:
        rep.errors.append(f"schema {prefix}{e}")
    return not errs


# --------------------------------------------------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------------------------------------------------
def _asc(xs) -> bool:
    return all(a < b for a, b in zip(xs, xs[1:]))


def _expand(ranges) -> set:
    s = set()
    for a, b in ranges or []:
        s.update(range(a, b + 1))
    return s


def _unit(v, tol=0.02) -> bool:
    return abs(float(np.linalg.norm(np.asarray(v, dtype=float))) - 1.0) <= tol


def decode_b64(s):
    try:
        return base64.b64decode(s, validate=True)
    except (binascii.Error, ValueError, TypeError):
        return None


def parse_pdb_atoms(text: str) -> list[dict]:
    """Minimal fixed-column PDB reader (ATOM and HETATM): enough to check chains, numbering, CA and B-factors."""
    atoms = []
    for line in text.splitlines():
        if line[:6] in ("ATOM  ", "HETATM") and len(line) >= 54:
            try:
                atoms.append({
                    "rec": line[:6].strip(), "name": line[12:16].strip(), "resn": line[17:20].strip(), "chain": line[21],
                    "resi": int(line[22:26]), "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
                    "b": float(line[60:66]) if len(line) >= 66 and line[60:66].strip() else None,
                })
            except ValueError:
                continue
    return atoms


def effective_caps(blob: dict, profile: dict | None = None) -> dict:
    """Capability flags in force: legacy blob caps underneath, profile.caps on top (the profile wins)."""
    caps = dict((blob or {}).get("caps") or {})
    if "traj" in caps and "has_traj" not in caps:
        caps["has_traj"] = caps["traj"]
    prof = profile if profile is not None else (blob or {}).get("profile")
    caps.update(((prof or {}).get("caps")) or {})
    return caps


def _lookup(obj, path: str):
    cur = obj
    for k in path.split("."):
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        else:
            return None
    return cur


@functools.lru_cache(maxsize=1)
def key_registry() -> dict:
    p = SCHEMA_DIR / "profile.keys.json"
    if not p.is_file():
        return {}
    return {e["key"]: e for e in json.loads(p.read_text(encoding="utf-8"))["keys"]}


@functools.lru_cache(maxsize=1)
def _key_patterns() -> tuple:
    """(regex, entry) for the registry keys that stand for a family of keys: ``<entry>`` matches one dotted part (finding.spec_orth.<entry>.some)."""
    out = []
    for k, e in key_registry().items():
        if "<" in k:
            out.append((re.compile("^" + re.sub(r"<[^<>.]+>", "[^.]+", re.escape(k).replace("\\<", "<").replace("\\>", ">")) + "$"), e))
    return tuple(out)


def registry_entry(reg: dict, key: str):
    """The registry entry of a text key: the key itself, else the pattern key it belongs to (None when no page reads it)."""
    e = reg.get(key)
    if e is not None:
        return e
    for pat, entry in _key_patterns():
        if pat.match(key):
            return entry
    return None


# --------------------------------------------------------------------------------------------------------------------
# profile
# --------------------------------------------------------------------------------------------------------------------
def check_profile(p: dict, name: str = "profile") -> Report:
    rep = Report(name)
    if not _schema_into(rep, "profile", p):
        return rep
    caps = p["caps"]
    rep.count("profile.text_keys", len(p.get("text") or {}))

    # text templates: every token resolves, every key is a known key
    reg = key_registry()
    for key, tpl in (p.get("text") or {}).items():
        if tpl.count("{") != tpl.count("}"):
            rep.error(f"text.{key}", "unbalanced braces")
            continue
        entry = registry_entry(reg, key)
        if reg and entry is None:
            rep.warn(f"text.{key}", "key is not in profile.keys.json (no page reads it)")
        for m in TOKEN.finditer(tpl):
            tok = m.group(1)
            if tok.startswith("@"):
                if _lookup(p, tok[1:]) is None:
                    rep.error(f"text.{key}", f"lookup {{{tok}}} does not resolve inside the profile")
            elif entry is not None and tok not in entry["tokens"]:
                rep.error(f"text.{key}", f"token {{{tok}}} is not one of {sorted(entry['tokens'])}")

    # frame strip
    fr = p.get("frame") or {}
    segs = fr.get("segments") or []
    prev = None
    for i, s in enumerate(segs):
        if s["start"] > s["end"]:
            rep.error(f"frame.segments[{i}]", "start after end")
        if prev is not None and s["start"] <= prev:
            rep.error(f"frame.segments[{i}]", "segments overlap or are not ascending")
        prev = s["end"]
    nat = fr.get("native_range")
    if nat and nat[0] > nat[1]:
        rep.error("frame.native_range", "first after last")
    if nat:
        for g in fr.get("glycan_sites") or []:
            if not nat[0] <= g <= nat[1]:
                rep.error("frame.glycan_sites", f"{g} outside native_range")
    tc = (p.get("target") or {}).get("construct") or {}
    if tc and nat and not (nat[0] <= tc["first"] <= tc["last"] <= nat[1]):
        rep.warn("target.construct", "construct lies outside frame.native_range")

    # specificity
    sp = p.get("specificity") or {}
    keys = [s["key"] for s in sp.get("species") or []]
    if len(set(keys)) != len(keys):
        rep.error("specificity.species", "duplicate keys")
    kset = set(keys)
    for s in sp.get("species") or []:
        if s.get("fallback") and s["fallback"] not in kset:
            rep.error("specificity.species", f"fallback {s['fallback']} is not a species key")
    for name_ in ("show", "headline", "decoy_species"):
        for k in sp.get(name_) or []:
            if k not in kset:
                rep.error(f"specificity.{name_}", f"{k} is not a species key")
    for d in sp.get("detail") or []:
        for k in ([d] if isinstance(d, str) else d.get("keys", [])):
            if k not in kset:
                rep.error("specificity.detail", f"{k} is not a species key")
        entry = d if isinstance(d, str) else d.get("group", "")
        if caps.get("species_panel"):
            for tail in ("some", "none"):
                if f"finding.spec_orth.{entry}.{tail}" not in (p.get("text") or {}):
                    rep.error("specificity.detail", f"text key finding.spec_orth.{entry}.{tail} is missing")
    for r in (sp.get("rules") or []) + ((p.get("rules") or {}).get("extra") or []):
        expr = r.get("expr") or {}
        for k in ([expr["species"]] if "species" in expr else []) + list(expr.get("min_species") or []):
            if k not in kset:
                rep.error(f"rules {r['id']}", f"expression uses {k}, not a species key")

    # checks
    ck = p.get("checks") or {}
    pd = ck.get("preset_default")
    if pd and pd not in ("lit", "relaxed", "strict", "bindcraft") and pd not in (ck.get("presets") or {}):
        rep.error("checks.preset_default", f"{pd} is neither a base preset nor in checks.presets")

    # capability promises against the sections that carry the data
    need = [
        ("olig", ("checks", "context", "partner"), "checks.context.partner"),
        ("rest_ctx", ("checks", "context", "rest"), "checks.context.rest"),
        ("species_panel", ("specificity", "species"), "specificity.species"),
        ("flag_region", ("flag_region",), "flag_region"),
        ("footprints", ("footprints",), "footprints"),
        ("glycan", ("frame", "glycan_sites"), "frame.glycan_sites"),
    ]
    for cap, path_, label in need:
        cur = p
        for k in path_:
            cur = cur.get(k) if isinstance(cur, dict) else None
        if caps.get(cap) and not cur:
            rep.error(f"caps.{cap}", f"is true but {label} is empty or missing")
    if caps.get("paralog_panel") and not (sp.get("paralogs") or {}).get("focus"):
        rep.error("caps.paralog_panel", "is true but specificity.paralogs.focus is missing")
    if caps.get("real") is not None:
        rep.warn("caps.real", "deprecated flag; new exports must not write it")
    for k in METRIC_CAPS + ("outcomes", "olig", "rest_ctx", "has_pae", "has_traj", "has_boltz", "species_panel") if caps else ():
        if k not in caps:    # an empty caps object is a source profile: the exporter writes every flag from the data
            rep.warn(f"caps.{k}", "flag missing (read as false); exporters write every flag")
    return rep


# --------------------------------------------------------------------------------------------------------------------
# campaign (RAW)
# --------------------------------------------------------------------------------------------------------------------
def _check_target(t: dict, rep: Report, caps: dict) -> None:
    w = "target"
    n = len(t["resi"])
    if not _asc(t["resi"]):
        rep.error(f"{w}.resi", "not strictly ascending")
    for k in ("resn", "rsa", "dom", "orig_resi"):
        if k in t and len(t[k]) != n:
            rep.error(f"{w}.{k}", f"{len(t[k])} entries, resi has {n}")
    if len(t["ss"]) != n:
        rep.error(f"{w}.ss", f"{len(t['ss'])} letters, resi has {n}")
    if len(t["ca"]) != 3 * n:
        rep.error(f"{w}.ca", f"{len(t['ca'])} numbers, expected {3 * n}")
    a = t["atoms"]
    na = len(a["res"])
    if len(a["x"]) != 3 * na:
        rep.error(f"{w}.atoms.x", f"{len(a['x'])} numbers for {na} atoms")
    for k in ("el", "name", "sasa"):
        if len(a[k]) != na:
            rep.error(f"{w}.atoms.{k}", f"{len(a[k])} entries for {na} atoms")
    if na and (max(a["res"]) >= n):
        rep.error(f"{w}.atoms.res", "index beyond resi")
    rset = set(t["resi"])
    for k in ("hotspots", "epitope", "glycan"):
        bad = [r for r in t[k] if r not in rset]
        if bad:
            rep.error(f"{w}.{k}", f"residues not in resi: {bad[:8]}")
    miss = [g for g in t["glycan"] if str(g) not in t["glycan_xyz"]]
    if miss:
        rep.error(f"{w}.glycan_xyz", f"no anchor for sites {miss}")
    if "A1" in t["sites"]:
        s = t["sites"]["A1"]
        if not _unit(s["normal"]):
            rep.error(f"{w}.sites.A1.normal", "not a unit vector")
        bad = [r for r in s["residues"] if r not in rset]
        if bad:
            rep.error(f"{w}.sites.A1.residues", f"residues not in resi: {bad[:8]}")
    if not _unit(t["membrane"]["axis"]):
        rep.error(f"{w}.membrane.axis", "not a unit vector")
    num = t.get("numbering") or {}
    if "orig_resi" in t and "input_offset" in num:
        bad = [(r, o) for r, o in zip(t["resi"], t["orig_resi"]) if o is not None and o != r + num["input_offset"]]
        if bad:
            rep.error(f"{w}.orig_resi", f"orig != uniprot + input_offset ({num['input_offset']}) for {bad[:4]}")
    if "unresolved" in t and _expand(t["unresolved"]) & rset:
        rep.error(f"{w}.unresolved", "ranges overlap residues of resi")
    if "rest" in t and t["rest"].get("ranges") and _expand(t["rest"]["ranges"]) & rset:
        rep.error(f"{w}.rest.ranges", "ranges overlap residues of resi")
    for k, lst in (t.get("hotspots_by_run") or {}).items():
        bad = [r for r in lst if r not in rset]
        if bad:
            rep.error(f"{w}.hotspots_by_run.{k}", f"residues not in resi: {bad[:8]}")
    if t.get("hotspots_by_run"):
        union = sorted(set().union(*map(set, t["hotspots_by_run"].values())))
        if union != sorted(t["hotspots"]):
            rep.error(f"{w}.hotspots", f"is not the union of hotspots_by_run ({union})")
    if t.get("epitope_by_run"):
        union = sorted(set().union(*map(set, t["epitope_by_run"].values())))
        if union != sorted(t["epitope"]):
            rep.warn(f"{w}.epitope", f"is not the union of epitope_by_run ({union[:8]}...)")
    if "glycan_atoms" in t:
        for k, g in t["glycan_atoms"].items():
            if len(g["xyz"]) != 3 * len(g["el"]):
                rep.error(f"{w}.glycan_atoms.{k}", "xyz and el have different lengths")
            if int(k) not in t["glycan"]:
                rep.warn(f"{w}.glycan_atoms.{k}", "site is not listed in glycan")
    if caps.get("olig") and not (t.get("olig") or {}).get("neigh"):
        rep.error(f"{w}.olig", "caps.olig is true but target.olig.neigh is missing")
    for k, xs in ((t.get("olig") or {}).get("neigh") or {}).items():
        if len(xs) % 3:
            rep.error(f"{w}.olig.neigh.{k}", "length is not a multiple of 3")
    if caps.get("glycan") and not t["glycan"]:
        rep.error(f"{w}.glycan", "caps.glycan is true but the list is empty")
    if not caps.get("glycan", True) and t["glycan"]:
        rep.warn(f"{w}.glycan", "caps.glycan is false but sites are listed")


def _run_present(run: dict) -> set:
    return _expand((run.get("target") or {}).get("resi_ranges"))


def _check_runs(raw: dict, rep: Report, caps: dict) -> dict:
    """Checks runs[] and returns run id -> set of present residues."""
    runs = raw.get("runs") or []
    present = {}
    ids = [r["id"] for r in runs]
    if len(set(ids)) != len(ids):
        rep.error("runs", "duplicate run ids")
    rset = set(raw["target"]["resi"])
    union = set()
    for i, r in enumerate(runs):
        w = f"runs[{i}] ({r['id']})"
        t = r.get("target")
        if not t:
            rep.warn(w, "no target block: the viewer cannot draw this run's construct")
            continue
        rr = t["resi_ranges"]
        if not all(a <= b for a, b in rr) or not _asc([a for a, _ in rr]):
            rep.error(f"{w}.target.resi_ranges", "ranges must be ascending, each first <= last")
        if any(rr[k][0] <= rr[k - 1][1] for k in range(1, len(rr))):
            rep.error(f"{w}.target.resi_ranges", "ranges overlap")
        pres = _expand(rr)
        present[r["id"]] = pres
        union |= pres
        gaps = [[rr[k - 1][1] + 1, rr[k][0] - 1] for k in range(1, len(rr)) if rr[k][0] - rr[k - 1][1] > 1]
        if "removed" in t and t["removed"] != gaps:
            rep.error(f"{w}.target.removed", f"{t['removed']} differs from the gaps between resi_ranges {gaps}")
        if "first" in t and rr and t["first"] != rr[0][0]:
            rep.error(f"{w}.target.first", "differs from resi_ranges")
        if "last" in t and rr and t["last"] != rr[-1][1]:
            rep.error(f"{w}.target.last", "differs from resi_ranges")
        if "n_resi" in t and t["n_resi"] != len(pres):
            rep.error(f"{w}.target.n_resi", f"{t['n_resi']} but {len(pres)} residues are present")
        cps = t.get("cut_points")
        if cps is not None:
            if not _asc(cps):
                rep.error(f"{w}.target.cut_points", "not ascending")
            if set(cps) - pres:
                rep.error(f"{w}.target.cut_points", f"cut points not present in the construct: {sorted(set(cps) - pres)[:6]}")
            for a, b in gaps:
                for flank in (a - 1, b + 1):
                    if flank not in cps:
                        rep.error(f"{w}.target.cut_points", f"flank {flank} of removed stretch [{a}, {b}] is missing")
        for k in ("hotspots", "epitope"):
            bad = [x for x in t.get(k, []) if x not in pres]
            if bad:
                rep.error(f"{w}.target.{k}", f"residues not in this run's construct: {bad}")
        # funnel
        seen = set()
        for j, f in enumerate(r.get("funnel") or []):
            if f["key"] in seen:
                rep.error(f"{w}.funnel[{j}]", f"duplicate key {f['key']}")
            seen.add(f["key"])
            if (f["n"] is None) != (f["status"] == "unknown"):
                rep.error(f"{w}.funnel[{j}]", "n is null exactly when status is 'unknown'")
        if r.get("traj_basis") == "attempted" and r.get("n_traj") is None:
            rep.error(f"{w}.traj_basis", "'attempted' but n_traj is unknown (null)")
        if caps.get("traj_total") and r.get("n_traj") is None:
            rep.warn(f"{w}.n_traj", "caps.traj_total is true but this run's n_traj is null")
    if runs and union != rset:
        extra, missing = sorted(rset - union), sorted(union - rset)
        rep.error("target.resi", f"is not the union of the runs' constructs (in resi only: {extra[:6]}, in runs only: {missing[:6]})")
    if caps.get("multi_run") and len(runs) < 2:
        rep.error("runs", "caps.multi_run is true but fewer than two runs are listed")
    if not caps.get("multi_run") and len(runs) > 1:
        rep.warn("runs", "several runs but caps.multi_run is false: the run filter stays hidden")
    return present


def _check_design(d: dict, i: int, raw: dict, present: dict, run_hot: dict, runs_ids: set, caps: dict, rep: Report) -> None:
    w = f"designs[{i}] ({d['id']})"
    n = d["len"]
    t = raw["target"]
    rset = t["_rset"]
    # A multi-chain binder's seq separates its chains with "/", which is not a
    # residue: len counts residues, so the separators come out of the count.
    residues = len(d["seq"].replace("/", ""))
    if residues != n:
        rep.error(f"{w}.seq", f"{residues} residues, len says {n}")
    if len(d["ca"]) != 3 * n:
        rep.error(f"{w}.ca", f"{len(d['ca'])} numbers, expected {3 * n}")
    if d["bc"] and (max(d["bc"]) >= n):
        rep.error(f"{w}.bc", "index beyond len")
    if "pl" in d and len(d["pl"]) != n:
        rep.error(f"{w}.pl", f"{len(d['pl'])} values, len says {n}")
    for k in ("tc", "hot"):
        if not _asc(d[k]):
            rep.error(f"{w}.{k}", "not strictly ascending")
    bad = [r for r in d["tc"] if r not in rset]
    if bad:
        rep.error(f"{w}.tc", f"residues not in target.resi: {bad[:6]}")
    own = present.get(d.get("run"))
    if own is not None:
        bad = [r for r in d["tc"] if r not in own]
        if bad:
            rep.error(f"{w}.tc", f"contacts on residues the design's own run lacks: {bad[:6]}")
    if set(d["hot"]) - set(d["tc"]):
        rep.error(f"{w}.hot", "hotspots that are not in tc")
    hs = run_hot.get(d.get("run"), t["hotspots"])
    bad = [r for r in d["hot"] if r not in hs]
    if bad:
        rep.error(f"{w}.hot", f"not hotspots of the design's run: {bad}")
    if "hot_of" in d and d["hot_of"] != len(hs):
        rep.error(f"{w}.hot_of", f"{d['hot_of']} but the run has {len(hs)} hotspots")
    if d["mode"] not in raw["modes"]:
        rep.error(f"{w}.mode", f"{d['mode']} is not defined in modes")
    if runs_ids and d.get("run") not in runs_ids:
        rep.error(f"{w}.run", f"{d.get('run')!r} is not a run id")
    bad = [g for g in d["glycan"] if int(g) not in t["glycan"]]
    if bad:
        rep.error(f"{w}.glycan", f"sites not in target.glycan: {bad}")
    missing = [g for g in t["glycan"] if str(g) not in d["glycan"]]
    if missing and caps.get("glycan"):
        rep.warn(f"{w}.glycan", f"no distance for sites {missing}")
    if any(p["p"] > n for p in d["liab"]):
        rep.error(f"{w}.liab", "position beyond len")
    status = d.get("status", "accepted")
    if status == "rejected":
        if "nm" not in d:
            rep.warn(f"{w}.nm", "rejected design without a near-miss summary")
        if not d.get("fail"):
            rep.warn(f"{w}.fail", "rejected design without failed filters")
    else:
        if d.get("fail"):
            rep.warn(f"{w}.fail", "accepted design lists failed filters (per-model failures are allowed, so this is a warning)")
        if "nm" in d:
            rep.error(f"{w}.nm", "near-miss summary on an accepted design")
    if "outcome" in d and not caps.get("outcomes"):
        rep.error(f"{w}.outcome", "present but caps.outcomes is false (the reveal switch would leak it)")
    if d.get("traj", 0) < 0:
        rep.error(f"{w}.traj", "negative")


def check_campaign(raw: dict, name: str = "campaign") -> Report:
    rep = Report(name)
    if not _schema_into(rep, "campaign", raw):
        return rep
    prof = raw["profile"]
    rep.merge(check_profile(prof, "profile"), "profile: ")
    caps = effective_caps(raw)
    t = raw["target"]
    rep.count("campaign.designs", len(raw["designs"]))
    rep.count("campaign.target_residues", len(t["resi"]))
    rep.count("campaign.runs", len(raw.get("runs") or []))
    t["_rset"] = set(t["resi"])  # transient helper, removed below
    try:
        _check_target(t, rep, caps)
        present = _check_runs(raw, rep, caps)
        runs = raw.get("runs") or []
        run_hot = {r["id"]: r["target"].get("hotspots", t["hotspots"]) for r in runs if r.get("target")}
        run_ids = {r["id"] for r in runs}
        ids = [d["id"] for d in raw["designs"]]
        if len(set(ids)) != len(ids):
            dup = sorted({x for x in ids if ids.count(x) > 1})
            rep.error("designs", f"duplicate ids {dup[:5]}")
        for i, d in enumerate(raw["designs"]):
            _check_design(d, i, raw, present, run_hot, run_ids, caps, rep)
        bb_runs: dict = {}
        for d in raw["designs"]:
            bb_runs.setdefault(d["bb"], set()).add(d.get("run"))
        shared = sorted(b for b, r in bb_runs.items() if len(r) > 1)
        if shared:
            rep.error("designs.bb", f"backbone ids used by more than one run (a family belongs to one run): {shared[:4]}")
        fam_bad = sorted({d["fam"] for d in raw["designs"] if d["fam"] != d["bb"]})
        if fam_bad:
            rep.warn("designs.fam", f"differs from bb for {fam_bad[:3]} (every data set so far has fam == bb)")
        # metrics promised by capability flags
        ds = raw["designs"]
        for k in METRIC_CAPS:
            have = any(d["m"].get(k) is not None for d in ds)
            if caps.get(k) and not have:
                rep.error(f"caps.{k}", f"is true but no design carries m.{k}")
            if k in caps and not caps[k] and have:
                rep.warn(f"caps.{k}", f"is false but designs carry m.{k}")
        for cap, key in (("olig", "olig_min"), ("rest_ctx", "rest_min"), ("flag_region", "stalk_n")):
            have = any(d["m"].get(key) is not None for d in ds)
            if caps.get(cap) and not have:
                rep.error(f"caps.{cap}", f"is true but no design carries m.{key}")
            if cap in caps and not caps[cap] and have:
                rep.warn(f"caps.{cap}", f"is false but designs carry m.{key}")
        has_out = any("outcome" in d for d in ds)
        if caps.get("outcomes") and not has_out:
            rep.error("caps.outcomes", "is true but no design carries an outcome")
        has_rej = any(d.get("status") == "rejected" for d in ds)
        if caps.get("rejected") and not has_rej:
            rep.error("caps.rejected", "is true but no design is rejected")
        if has_rej and "rejected" in caps and not caps["rejected"]:
            rep.error("caps.rejected", "is false but rejected designs are present (they would show as accepted)")
        if caps.get("weights") and not any("wt" in d for d in ds):
            rep.error("caps.weights", "is true but no design carries wt")
        if caps.get("has_fold") and not raw.get("fold"):
            rep.error("caps.has_fold", "is true but fold is missing")
        if raw.get("fold") and "has_fold" in caps and not caps["has_fold"]:
            rep.warn("fold", "present but caps.has_fold is false")
        if caps.get("has_traj"):
            if not any(d["traj"] > 0 for d in ds):
                rep.error("caps.has_traj", "is true but every design has traj 0")
        # modes
        for k, m in raw["modes"].items():
            bad = [r for r in m.get("top", []) if r not in t["_rset"]] + [r for r in m.get("enr", []) if r not in t["_rset"]]
            if bad:
                rep.error(f"modes.{k}", f"residues not in target.resi: {bad[:6]}")
            if "n" in m:
                cnt = sum(1 for d in ds if d["mode"] == k and d.get("status", "accepted") == "accepted")
                if cnt != m["n"]:
                    rep.warn(f"modes.{k}.n", f"{m['n']} but {cnt} accepted designs carry this mode")
            if "site" in m and m["site"] not in t["sites"]:
                rep.error(f"modes.{k}.site", f"{m['site']} is not in target.sites")
        # per-run counts
        for r in runs:
            if "n_designs" in r and "accepted" in r["n_designs"]:
                cnt = sum(1 for d in ds if d.get("run") == r["id"] and d.get("status", "accepted") == "accepted")
                if cnt != r["n_designs"]["accepted"]:
                    rep.error(f"runs[{r['id']}].n_designs.accepted", f"{r['n_designs']['accepted']} but {cnt} accepted designs reference the run")
            if "n_designs" in r and "rejected" in r["n_designs"] and "rejected_with_structure" in r["n_designs"]:
                cnt = sum(1 for d in ds if d.get("run") == r["id"] and d.get("status") == "rejected")
                if cnt != r["n_designs"]["rejected_with_structure"]:
                    rep.error(f"runs[{r['id']}].n_designs.rejected_with_structure", f"{r['n_designs']['rejected_with_structure']} but {cnt} rejected designs with a pose reference the run")
            if r.get("traj_basis") == "none":
                bad = [d["id"] for d in ds if d.get("run") == r["id"] and d["traj"] != 0]
                if bad:
                    rep.error(f"runs[{r['id']}].traj_basis", f"'none' but designs carry traj != 0 ({bad[:3]})")
        _check_fold(raw, rep)
        _check_diversity(raw, rep)
    finally:
        t.pop("_rset", None)
    return rep


def _check_fold(raw: dict, rep: Report) -> None:
    f = raw.get("fold")
    if not f:
        return
    n = f["n"]
    if len(f["ids"]) != n:
        rep.error("fold.ids", f"{len(f['ids'])} ids for n = {n}")
    if len(set(f["ids"])) != len(f["ids"]):
        rep.error("fold.ids", "duplicate ids")
    dids = {d["id"] for d in raw["designs"]}
    bad = [x for x in f["ids"] if x not in dids]
    if bad:
        rep.error("fold.ids", f"ids that are not designs: {bad[:4]}")
    buf = decode_b64(f["tm"])
    if buf is None:
        rep.error("fold.tm", "not valid base64")
        return
    if len(buf) != n * n:
        rep.error("fold.tm", f"{len(buf)} bytes, expected {n * n}")
        return
    m = np.frombuffer(buf, dtype=np.uint8).reshape(n, n)
    if not np.array_equal(m, m.T):
        rep.error("fold.tm", "matrix is not symmetric")
    if not np.all(np.diag(m) == 100):
        rep.warn("fold.tm", "diagonal is not 100")


def _check_diversity(raw: dict, rep: Report) -> None:
    dv = raw.get("diversity")
    if not dv:
        return
    run_ids = {r["id"] for r in raw.get("runs") or []}
    bb_by_run: dict = {}
    for d in raw["designs"]:
        bb_by_run.setdefault(d.get("run"), set()).add(d["bb"])
    for rid, blk in (dv.get("runs") or {}).items():
        if run_ids and rid not in run_ids:
            rep.error(f"diversity.runs.{rid}", "not a run id")
        _check_div_block(blk, f"diversity.runs.{rid}", rep, bb_by_run.get(rid))
    if dv.get("pooled"):
        _check_div_block(dv["pooled"], "diversity.pooled", rep, set().union(*bb_by_run.values()) if bb_by_run else None)
        runs = raw.get("runs") or []
        epochs = {(r.get("epoch") or {}).get("hash") for r in runs}
        if len(epochs) > 1 and dv["pooled"].get("estimates_valid"):
            rep.error("diversity.pooled.estimates_valid", "runs have different epochs: pooled coverage and forecasts must not be given")
        if not dv["pooled"].get("estimates_valid", False) and dv["pooled"].get("advice"):
            rep.error("diversity.pooled.advice", "advice present although estimates_valid is false")


def _check_div_block(blk: dict, w: str, rep: Report, bbs) -> None:
    for basis, cl in (blk.get("clusters") or {}).items():
        sizes = sum(c["size"] for c in cl)
        fams = [f for c in cl for f in c["fams"]]
        if len(fams) != len(set(fams)):
            rep.error(f"{w}.clusters.{basis}", "a family sits in two clusters")
        if sizes != len(fams):
            rep.error(f"{w}.clusters.{basis}", "cluster size differs from the number of listed families")
        if bbs is not None and set(fams) - bbs:
            rep.error(f"{w}.clusters.{basis}", f"families that are not backbones of this run: {sorted(set(fams) - bbs)[:3]}")
        c = (blk.get("counts") or {}).get(basis)
        if c:
            if c.get("S") is not None and c["S"] != len(cl):
                rep.error(f"{w}.counts.{basis}.S", f"{c['S']} but {len(cl)} clusters are listed")
            if c.get("n") is not None and c["n"] != len(fams):
                rep.error(f"{w}.counts.{basis}.n", f"{c['n']} but {len(fams)} families are listed")
            if c.get("f1") is not None and c.get("f2") is not None and c["f1"] + c["f2"] > c.get("S", 10 ** 9):
                rep.error(f"{w}.counts.{basis}", "f1 + f2 exceeds S")
    if blk.get("estimates_valid") is False and blk.get("advice"):
        rep.error(f"{w}.advice", "advice present although estimates_valid is false")


# --------------------------------------------------------------------------------------------------------------------
# inspect (INS)
# --------------------------------------------------------------------------------------------------------------------
def _check_ins_design(tag: str, d: dict, ins: dict, caps: dict, rep: Report) -> None:
    w = f"designs.{tag}"
    n = d["n"]
    if d["tag"] != tag:
        rep.error(f"{w}.tag", f"{d['tag']!r} differs from the key")
    for k in ("seq", "ss", "plddt_res", "rsa_res", "dsasa_res"):
        # seq separates the chains of a multi-chain binder with "/", which is not
        # a residue; the per-residue arrays beside it have no separator.
        got = len(d[k].replace("/", "")) if k == "seq" and isinstance(d[k], str) else len(d[k])
        if got != n:
            rep.error(f"{w}.{k}", f"{got} entries, n is {n}")
    if d.get("ddg_res") is not None and len(d["ddg_res"]) != n:
        rep.error(f"{w}.ddg_res", f"{len(d['ddg_res'])} entries, n is {n}")
    tres = d.get("t_resi") or ins["target_resi"]
    tset = set(tres)
    for k in ("plddt_t", "t_dsasa", "t_ddg"):
        if d.get(k) is not None and len(d[k]) != len(tres):
            rep.error(f"{w}.{k}", f"{len(d[k])} entries, the target has {len(tres)} residues")
    for j, c in enumerate(d["contacts"]):
        if c["t"] not in tset or not 1 <= c["b"] <= n:
            rep.error(f"{w}.contacts[{j}]", f"residue out of range (t {c['t']}, b {c['b']}, n {n})")
            break
        if c["d"] > 4.55:
            rep.warn(f"{w}.contacts[{j}]", f"distance {c['d']} A above the 4.5 A contact cut-off")
            break
    ct = sorted({c["t"] for c in d["contacts"]})
    cb = sorted({c["b"] for c in d["contacts"]})
    if ct != sorted(d["t_contact"]):
        rep.warn(f"{w}.t_contact", "differs from the target residues of contacts")
    if cb != sorted(d["b_contact"]):
        rep.warn(f"{w}.b_contact", "differs from the binder residues of contacts")
    if set(d["hot"]) - set(d["t_contact"]):
        rep.warn(f"{w}.hot", "hotspots that are not contacted")
    for j, h in enumerate(d["hbonds"]):
        if h["t"] not in tset or not 1 <= h["b"] <= n:
            rep.error(f"{w}.hbonds[{j}]", "residue out of range")
            break
    for j, h in enumerate(d["salt"]):
        if h["t"] not in tset or not 1 <= h["b"] <= n:
            rep.error(f"{w}.salt[{j}]", "residue out of range")
            break
    for j, b in enumerate(d["buns"]):
        if (b["side"] == "b" and not 1 <= b["resi"] <= n) or (b["side"] == "t" and b["resi"] not in tset):
            rep.error(f"{w}.buns[{j}]", "residue out of range")
            break
    for pa in d["patches"]:
        if any(not 1 <= r <= n for r in pa["res"]):
            rep.error(f"{w}.patches", "residue out of range")
    for sb in d["ssbonds"]:
        if any(r > n for r in sb):
            rep.error(f"{w}.ssbonds", "residue out of range")
    # binder structure text
    # Every atom of the block: it holds the binder and nothing else, and a binder
    # of several chains keeps its own chain letters, so filtering on "B" used to
    # drop all but the first chain and then report the rest as missing.
    atoms = parse_pdb_atoms(d["pdb_b"])
    if not atoms:
        rep.error(f"{w}.pdb_b", "no binder atoms")
    else:
        resis = sorted({a["resi"] for a in atoms})
        if resis != list(range(1, n + 1)):
            rep.error(f"{w}.pdb_b", f"binder residues are not 1..{n}")
        cas = [a for a in atoms if a["name"] == "CA"]
        if len(cas) != n:
            rep.error(f"{w}.pdb_b", f"{len(cas)} CA atoms for {n} residues")
        bs = [a["b"] for a in cas if a["b"] is not None]
        # synthetic data keeps the source structure's B-factors here; real data carries pLDDT (0..100)
        strict_b = (ins["profile"].get("kind") != "synthetic") and not d.get("source_note")
        flag = rep.error if strict_b else (lambda *a: None)
        if bs and (min(bs) < 0 or max(bs) > 100.0):
            flag(f"{w}.pdb_b", "B-factors (pLDDT) outside 0..100")
        letters = "".join(THREE_TO_ONE.get(a["resn"], "X") for a in sorted(cas, key=lambda a: a["resi"]))
        # The PDB block carries no separator, and a multi-chain binder's seq does.
        if letters != d["seq"].upper().replace("/", ""):
            rep.error(f"{w}.pdb_b", "residue names differ from seq")
        if bs and len(bs) == len(d["plddt_res"]):
            diff = float(np.max(np.abs(np.asarray(bs) - np.asarray(d["plddt_res"], dtype=float))))
            if diff > 0.51:
                flag(f"{w}.plddt_res", f"differs from the B-factors of pdb_b by up to {diff:.1f}")
    # PAE
    if caps.get("has_pae") and d.get("pae") is None:
        rep.error(f"{w}.pae", "caps.has_pae is true but pae is null")
    if "has_pae" in caps and not caps["has_pae"] and d.get("pae") is not None:
        rep.error(f"{w}.pae", "caps.has_pae is false but a PAE matrix is present (the page would show a score the profile hides)")
    if d.get("pae") is not None:
        pn = d.get("pae_n")
        buf = decode_b64(d["pae"])
        if pn is None:
            rep.error(f"{w}.pae_n", "pae without pae_n")
        elif buf is None:
            rep.error(f"{w}.pae", "not valid base64")
        else:
            tot = pn[0] + pn[1]
            if len(buf) != tot * tot:
                rep.error(f"{w}.pae", f"{len(buf)} bytes, expected {tot * tot}")
            if pn[0] != n:
                rep.error(f"{w}.pae_n", f"binder residues {pn[0]} differ from n {n}")
            expect = len(d.get("pae_tres") or tres)
            if pn[1] != expect:
                rep.error(f"{w}.pae_n", f"target residues {pn[1]} differ from {expect}")
    # poses and seeds
    if not d["poses"]:
        rep.error(f"{w}.poses", "no entries (the first pose is the design itself)")
    for j, p in enumerate(d["poses"]):
        R = np.asarray(p["R"], dtype=float).reshape(3, 3)
        if abs(np.linalg.det(R) - 1.0) > 0.02 or np.max(np.abs(R @ R.T - np.eye(3))) > 0.02:
            rep.error(f"{w}.poses[{j}].R", "not a proper rotation")
    for j, p in enumerate(d["seeds"]):
        R = np.asarray(p["R"], dtype=float).reshape(3, 3)
        if abs(np.linalg.det(R) - 1.0) > 0.02:
            rep.error(f"{w}.seeds[{j}].R", "not a proper rotation")
    if d["seeds"] and "has_seeds" in caps and not caps["has_seeds"]:
        rep.warn(f"{w}.seeds", "present but caps.has_seeds is false")
    # siblings
    if not any(s["id"] == d["id"] for s in d["sibs"]):
        rep.warn(f"{w}.sibs", "does not list the design itself")
    for s in d["sibs"]:
        if any(x > n for x in s.get("diff", [])):
            rep.error(f"{w}.sibs", "diff position beyond n")
        if "ndiff" in s and "diff" in s and s["ndiff"] != len(s["diff"]):
            rep.warn(f"{w}.sibs", f"ndiff {s['ndiff']} but {len(s['diff'])} positions listed")
    # metrics echo
    m = d["metrics"]
    if d["mem_margin"] is not None and m.get("mem_excess") is not None and abs(d["mem_margin"] + m["mem_excess"]) > 0.1:
        rep.warn(f"{w}.mem_margin", "is not -metrics.mem_excess")
    if any(p["pos"] > n for p in d["liab"]):
        rep.error(f"{w}.liab", "position beyond n")
    if d.get("ctx"):
        c = d["ctx"]
        if c.get("pair") and not 1 <= c["pair"][0] <= n:
            rep.error(f"{w}.ctx.pair", "binder residue out of range")
    if caps.get("rest_ctx") and not d.get("ctx"):
        rep.error(f"{w}.ctx", "caps.rest_ctx is true but the design has no ctx block")
    if caps.get("olig") and not d.get("octx"):
        rep.error(f"{w}.octx", "caps.olig is true but the design has no octx block")
    if d.get("octx") and d["octx"].get("olig_buried"):
        bad = [r for r in d["octx"]["olig_buried"] if r not in d["t_contact"]]
        if bad:
            rep.warn(f"{w}.octx.olig_buried", "residues that are not contacted")
    if d.get("status") == "rejected" and "rejected" in caps and not caps["rejected"]:
        rep.error(f"{w}.status", "rejected design but caps.rejected is false")


def check_inspect(ins: dict, raw: dict | None = None, name: str = "inspect") -> Report:
    rep = Report(name)
    if not _schema_into(rep, "inspect", ins):
        return rep
    prof = ins["profile"]
    rep.merge(check_profile(prof, "profile"), "profile: ")
    caps = effective_caps(ins)
    rep.count("inspect.designs", len(ins["designs"]))
    if raw is not None:
        caps = {**effective_caps(raw), **{k: v for k, v in effective_caps(ins).items()}}
    n = len(ins["target_resi"])
    if len(ins["target_seq"]) != n:
        rep.error("target_seq", f"{len(ins['target_seq'])} letters for {n} residues")
    if len(ins["t_ss"]) != n:
        rep.error("t_ss", f"{len(ins['t_ss'])} letters for {n} residues")
    if not _asc(ins["target_resi"]):
        rep.error("target_resi", "not strictly ascending")
    ca = [a for a in parse_pdb_atoms(ins["target_pdb"]) if a["chain"] == "A" and a["rec"] == "ATOM" and a["name"] == "CA"]
    if [a["resi"] for a in ca] != ins["target_resi"]:
        rep.error("target_pdb", "CA records of chain A differ from target_resi")
    order = ins["order"]
    if len(set(order)) != len(order):
        rep.error("order", "duplicate tags")
    if set(order) - set(ins["designs"]):
        rep.error("order", f"tags that are not designs: {sorted(set(order) - set(ins['designs']))[:4]}")
    if set(ins["designs"]) - set(order):
        rep.warn("order", "designs that are not in order (they cannot be selected)")
    for tag, d in ins["designs"].items():
        _check_ins_design(tag, d, ins, caps, rep)
    for k in ("hotspots", "epitope", "glycan"):
        bad = [r for r in ins[k] if r not in set(ins["target_resi"])]
        if bad:
            rep.error(k, f"residues not in target_resi: {bad[:6]}")
    if caps.get("olig") and not (ins.get("olig") or {}).get("neigh"):
        rep.error("olig", "caps.olig is true but olig.neigh is missing")
    if caps.get("rest_ctx") and not (ins.get("rest") or {}).get("traces"):
        rep.warn("rest", "caps.rest_ctx is true but rest.traces is empty (no removed region is drawn)")
    for r in ins.get("runs") or []:
        if r["target"] and not r["target"].get("resi_ranges"):
            rep.error(f"runs[{r['id']}]", "target without resi_ranges")
    if raw is not None:
        rep.merge(_cross_raw_ins(raw, ins))
    return rep


def _cross_raw_ins(raw: dict, ins: dict) -> Report:
    rep = Report("cross raw/inspect")
    t = raw["target"]
    if ins["target_resi"] != t["resi"]:
        rep.error("inspect.target_resi", "differs from campaign target.resi")
    for k in ("hotspots", "epitope", "glycan"):
        if sorted(ins[k]) != sorted(t[k]):
            rep.error(f"inspect.{k}", f"differs from campaign target.{k}")
    by_id = {d["id"]: d for d in raw["designs"]}
    by_tag = {d.get("inspect"): d for d in raw["designs"] if d.get("inspect")}
    for tag, d in ins["designs"].items():
        r = by_tag.get(tag) or by_id.get(d["id"]) or by_id.get(tag)
        if r is None:
            if not d.get("source_note"):
                rep.warn(f"inspect.designs.{tag}", "no campaign design with this tag or id")
            continue
        if r["seq"] != d["seq"]:
            rep.error(f"inspect.designs.{tag}.seq", "differs from the campaign design")
        if r.get("run") != d.get("run") and d.get("run") is not None:
            rep.error(f"inspect.designs.{tag}.run", "differs from the campaign design")
        if "tc" in r and sorted(r["tc"]) != sorted(d["t_contact"]) and d.get("t_resi") is None:
            rep.warn(f"inspect.designs.{tag}.t_contact", "differs from the campaign tc")
        if abs(r["m"]["dsasa"] - d["dsasa"]) > 1.0 and d.get("t_resi") is None:
            rep.warn(f"inspect.designs.{tag}.dsasa", f"{d['dsasa']} vs campaign m.dsasa {r['m']['dsasa']}")
        a = [x for x in parse_pdb_atoms(d["pdb_b"]) if x["name"] == "CA"]    # every binder chain, not just the first
        if a and len(r["ca"]) == 3 * len(a):
            xyz = np.array([[x["x"], x["y"], x["z"]] for x in sorted(a, key=lambda x: x["resi"])])
            diff = float(np.max(np.abs(xyz.reshape(-1) - np.asarray(r["ca"], dtype=float))))
            if diff > 0.2:
                rep.warn(f"inspect.designs.{tag}.pdb_b", f"CA positions differ from campaign ca by up to {diff:.2f} A")
    return rep


# --------------------------------------------------------------------------------------------------------------------
# homolog (HOM)
# --------------------------------------------------------------------------------------------------------------------
def species_entries(profile: dict) -> list:
    """The species the profile asks the pages to show: profile.specificity.species."""
    return ((profile or {}).get("specificity") or {}).get("species") or []


def resolve_species(hom: dict, entry: dict):
    """The orth entry a profile species entry points at (organism name or 'acc:<accession>'), or None."""
    m = entry["match"]
    for o in hom["orth"]:
        if m.startswith("acc:"):
            if o["acc"] == m[4:]:
                return o
        elif o["org"] == m:
            return o
    return None


def check_homolog(hom: dict, raw: dict | None = None, profile: dict | None = None, name: str = "homolog") -> Report:
    rep = Report(name)
    if not _schema_into(rep, "homolog", hom):
        return rep
    n = len(hom["resis"])
    rep.count("homolog.residues", n)
    rep.count("homolog.orthologs", len(hom["orth"]))
    rep.count("homolog.paralogs", len(hom["par"]))
    if not _asc(hom["resis"]):
        rep.error("resis", "not strictly ascending")
    for k in ("seq", "cons", "ident", "grade", "cover"):
        if len(hom[k]) != n:
            rep.error(k, f"{len(hom[k])} entries for {n} residues")
    for i, o in enumerate(hom["orth"]):
        if len(o["map"]) != n:
            rep.error(f"orth[{i}] ({o['acc']}).map", f"{len(o['map'])} letters for {n} residues")
    for i, o in enumerate(hom["par"]):
        if len(o["map"]) != n:
            rep.error(f"par[{i}] ({o['gene']}).map", f"{len(o['map'])} letters for {n} residues")
    core = [o for o in hom["orth"] if not o.get("extra")]
    if hom["n_orth"] != len(core):
        rep.warn("n_orth", f"{hom['n_orth']} but {len(core)} orthologs are not marked extra")
    if hom.get("consurf") and len(hom["consurf"].get("grade", [])) != n:
        rep.error("consurf.grade", f"{len(hom['consurf'].get('grade', []))} entries for {n} residues")
    prof = profile or (raw or {}).get("profile")
    if raw is not None and hom["resis"] != raw["target"]["resi"]:
        rep.error("resis", "differs from campaign target.resi")
    if prof is not None:
        caps = effective_caps({}, prof)
        if caps.get("species_panel"):
            for e in species_entries(prof):
                if resolve_species(hom, e) is None and not e.get("optional") and not e.get("fallback"):
                    rep.error(f"profile.specificity.species[{e['key']}]", f"no ortholog matches {e['match']!r}")
        if caps.get("paralog_panel"):
            focus = ((prof.get("specificity") or {}).get("paralogs") or {}).get("focus")
            genes = {p["gene"] for p in hom["par"]}
            if not hom["par"]:
                rep.error("par", "caps.paralog_panel is true but par is empty")
            elif focus and focus["key"] not in genes:
                rep.error("par", f"focus paralog {focus['key']} is not in par")
            order = ((prof.get("specificity") or {}).get("paralogs") or {}).get("order") or []
            miss = [g for g in order if g not in genes]
            if miss:
                rep.warn("par", f"profile paralog order lists genes not in par: {miss}")
        if caps.get("consurf") and not hom.get("consurf"):
            rep.error("consurf", "caps.consurf is true but the overlay is missing")
        if "consurf" in caps and not caps["consurf"] and hom.get("consurf"):
            rep.warn("consurf", "overlay present but caps.consurf is false (it stays hidden)")
    return rep


# --------------------------------------------------------------------------------------------------------------------
# manifest and review events
# --------------------------------------------------------------------------------------------------------------------
def check_manifest(m: dict, name: str = "manifest") -> Report:
    rep = Report(name)
    if not _schema_into(rep, "manifest", m):
        return rep
    ids = [r["id"] for r in m["runs"]]
    if len(set(ids)) != len(ids):
        rep.error("runs", "duplicate run ids")
    paths = [f["path"] for f in (m["bundle"].get("files") or [])]
    if len(set(paths)) != len(paths):
        rep.error("bundle.files", "duplicate paths")
    for i, s in enumerate(m["sources"]):
        if s.get("run") and s["run"] not in ids:
            rep.error(f"sources[{i}].run", f"{s['run']!r} is not a run id")
    c = m.get("counts") or {}
    if "runs" in c and c["runs"] != len(ids):
        rep.warn("counts.runs", f"{c['runs']} but {len(ids)} runs are listed")
    return rep


def check_review_events(events: list[dict], name: str = "reviews") -> Report:
    rep = Report(name)
    seen = set()
    for i, e in enumerate(events):
        if not _schema_into(rep, "review_event", e, f"[{i}] "):
            continue
        if e["id"] in seen:
            rep.error(f"[{i}].id", f"duplicate event id {e['id']}")
        seen.add(e["id"])
        if e["decision"] == "clear" and e.get("sequences"):
            rep.error(f"[{i}].sequences", "a clear event carries no sequences")
        if e["decision"] != "keep" and e.get("sequences"):
            rep.warn(f"[{i}].sequences", "only keep decisions put sequences on the order list")
        if e["unit"] == "design" and e["decision"] == "keep" and e.get("sequences") and e["sequences"] != [e["design"]]:
            rep.warn(f"[{i}].sequences", "a keep on one design lists only that design")
    return rep


def read_jsonl(path: Path) -> list[dict]:
    out = []
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines()):
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError as ex:
                out.append({"__invalid__": f"line {i + 1}: {ex}"})
    return out


# --------------------------------------------------------------------------------------------------------------------
# whole bundle and legacy blobs
# --------------------------------------------------------------------------------------------------------------------
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_blobs(campaign: dict | None, inspect: dict | None, homolog: dict | None, name: str = "blobs") -> Report:
    """All blobs of one project with the cross rules."""
    rep = Report(name)
    if campaign is not None:
        rep.merge(check_campaign(campaign, "campaign"), "campaign ")
    if inspect is not None:
        rep.merge(check_inspect(inspect, campaign, "inspect"), "inspect ")
    if homolog is not None:
        rep.merge(check_homolog(homolog, campaign, None, "homolog"), "homolog ")
    if campaign is not None and inspect is not None and campaign["profile"] != inspect["profile"]:
        rep.error("profile", "campaign.profile and inspect.profile differ")
    return rep


def check_bundle(bundle_dir, name: str | None = None, verify_hashes: bool = True) -> Report:
    """Validate data/*.json and reviews/*.jsonl of a bundle folder, plus the manifest hashes and the file links."""
    root = Path(bundle_dir)
    rep = Report(name or str(root))
    data = root / "data"

    def load(fn):
        p = data / fn
        return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else None

    camp, ins, hom, prof, man = (load(f) for f in ("campaign.json", "inspect.json", "homolog.json", "profile.json", "manifest.json"))
    if camp is None:
        rep.error("data/campaign.json", "missing")
    rep.merge(check_blobs(camp, ins, hom))
    if prof is not None:
        rep.merge(check_profile(prof, "profile"), "data/profile.json ")
        if camp is not None and camp["profile"] != prof:
            rep.error("data/profile.json", "differs from campaign.profile")
    if man is None:
        rep.error("data/manifest.json", "missing")
    else:
        rep.merge(check_manifest(man), "manifest ")
        if camp is not None:
            if man["caps"] != camp["profile"]["caps"]:
                rep.error("manifest.caps", "differs from campaign.profile.caps")
            ids = {r["id"] for r in camp.get("runs") or []}
            if ids and ids != {r["id"] for r in man["runs"]}:
                rep.error("manifest.runs", "run ids differ from campaign.runs")
        if verify_hashes:
            for f in man["bundle"].get("files") or []:
                p = root / f["path"]
                if not p.is_file():
                    rep.error(f"manifest.bundle.files {f['path']}", "file is missing")
                elif p.stat().st_size != f["bytes"] or _sha256(p) != f["sha256"]:
                    rep.error(f"manifest.bundle.files {f['path']}", "size or sha256 differs from the file")
            for k, rel in (man["bundle"].get("pages") or {}).items():
                if not (root / rel).is_file():
                    rep.error(f"manifest.bundle.pages.{k}", f"{rel} is missing")
            for k, rel in (man["bundle"].get("data") or {}).items():
                if not (root / rel).is_file():
                    rep.error(f"manifest.bundle.data.{k}", f"{rel} is missing")
    for p in sorted((root / "reviews").glob("*.jsonl")) if (root / "reviews").is_dir() else []:
        evs = read_jsonl(p)
        bad = [e for e in evs if "__invalid__" in e]
        for b in bad:
            rep.error(f"reviews/{p.name}", b["__invalid__"])
        rep.merge(check_review_events([e for e in evs if "__invalid__" not in e]), f"reviews/{p.name} ")
    if camp is not None and (root / "structures").is_dir():
        # the file naming of structures/ is the bundle writer's choice: <id>.pdb with '/' written as '__' (or the part after '/') is accepted anywhere below it
        names = {p.name for p in (root / "structures").rglob("*.pdb")}
        for d in camp["designs"]:
            cand = {f"{d['id']}.pdb", f"{d['id'].replace('/', '__')}.pdb", f"{d['id'].split('/')[-1]}.pdb"}
            if d.get("inspect") and not (cand & names):
                rep.warn(f"structures/{d['id']}.pdb", "no structure file for an inspected design")
    return rep


def load_example_profile(pid: str) -> dict:
    return json.loads((EXAMPLE_DIR / f"{pid}.profile.json").read_text(encoding="utf-8"))


def inject_profile(blob: dict, profile: dict) -> dict:
    """Older blobs predate profiles: a shallow copy with a profile added under 'profile'."""
    out = dict(blob)
    out["profile"] = profile
    return out


def check_examples(example_dir: Path | None = None) -> Report:
    """Every ``*.profile.json`` of the examples folder, and every example project ``<id>.{campaign,inspect,homolog,manifest}.json`` with its
    ``<id>.reviews.jsonl``."""
    d = Path(example_dir) if example_dir else EXAMPLE_DIR
    rep = Report("examples")
    for p in sorted(d.glob("*.profile.json")):
        rep.merge(check_profile(json.loads(p.read_text(encoding="utf-8")), p.name), f"{p.name} ")
    for cp in sorted(d.glob("*.campaign.json")):
        pid = cp.name[: -len(".campaign.json")]
        blobs = {}
        for kind in ("campaign", "inspect", "homolog", "manifest"):
            p = d / f"{pid}.{kind}.json"
            if p.is_file():
                blobs[kind] = json.loads(p.read_text(encoding="utf-8"))
        rep.merge(check_blobs(blobs.get("campaign"), blobs.get("inspect"), blobs.get("homolog"), f"{pid} example"), f"{pid} example ")
        if "manifest" in blobs:
            rep.merge(check_manifest(blobs["manifest"]), f"{pid} manifest ")
        p = d / f"{pid}.reviews.jsonl"
        if p.is_file():
            evs = read_jsonl(p)
            for e in evs:
                if "__invalid__" in e:
                    rep.error(p.name, e["__invalid__"])
            rep.merge(check_review_events([e for e in evs if "__invalid__" not in e]), f"{p.name} ")
    return rep


# --------------------------------------------------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("file", nargs="?", help="one JSON (or JSONL for review_event) file; needs --kind")
    ap.add_argument("--kind", choices=KINDS)
    ap.add_argument("--examples", action="store_true", help="example profiles and example blobs in docs/schema/examples")
    ap.add_argument("--bundle", help="bundle folder")
    ap.add_argument("--no-hashes", action="store_true", help="skip the sha256 check of bundle files")
    ap.add_argument("-q", "--quiet", action="store_true", help="hide warnings")
    a = ap.parse_args(argv)
    reports = []
    if a.examples:
        reports.append(check_examples())
    if a.bundle:
        reports.append(check_bundle(a.bundle, verify_hashes=not a.no_hashes))
    if a.file:
        if not a.kind:
            ap.error("--kind is required with a file")
        path = Path(a.file)
        rep = Report(path.name)
        if a.kind == "review_event":
            evs = read_jsonl(path)
            rep.merge(check_review_events([e for e in evs if "__invalid__" not in e]))
            for e in evs:
                if "__invalid__" in e:
                    rep.error(path.name, e["__invalid__"])
        else:
            obj = json.loads(path.read_text(encoding="utf-8"))
            fn = {"campaign": check_campaign, "inspect": check_inspect, "homolog": check_homolog, "manifest": check_manifest, "profile": check_profile}[a.kind]
            rep.merge(fn(obj, name=path.name))
        reports.append(rep)
    if not reports:
        ap.print_help()
        return 2
    bad = 0
    for r in reports:
        print(r.format(warnings=not a.quiet))
        bad += len(r.errors)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
