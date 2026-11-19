"""The whole analysis as cached functions: readers, target frame, poses and context checks, interface analysis, families
and diversity, modes and near-miss records, specificity panel, and the assembled in-memory ``Project``.

Every stage writes its result under ``work/<project>/`` keyed by the hashes of its inputs and of the code it ran, so a
rerun with a warm cache only reads JSON back (docs/export.md has the layout and the timings). Nothing here writes a blob
for the viewer: ``export.py`` turns a ``Project`` into the validated blobs and ``bundle.py`` into the static site.

Numbering follows docs/CONTRACTS.md: ``RunData`` and ``DesignAnalysis`` use construct numbers, everything stored in a
pose record is already numbered as the reference chain (UniProt numbering) and placed in its frame. The exporter converts the
construct-numbered interface values once.

Nothing here is specific to a target: the UniProt accession, the reference structure id, its chain and partner chains and the
profile are inputs (``ProjectSpec``). Without an accession no UniProt, SIFTS or AlphaFold record is read and the specificity
panel is empty; the frame then comes from the reference structure and the construct alone.
"""
from __future__ import annotations

import concurrent.futures as cf
import contextlib
import dataclasses
import datetime as _dt
import functools
import hashlib
import inspect as _inspect
import json
import multiprocessing as mp
import os
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Mapping, Optional, Sequence

import numpy as np

from . import __version__, config, families, interface, model, modes, specificity
from . import refdata as rd
from . import structure as S
from . import target_frame as TF
from .model import DesignRecord, RunData, binder_chain_count
from .readers import detect_reader, get_reader

PIPELINE_VERSION = 1          # bump when a stage's output changes shape; part of every cache key
POSE_VERSION = 2
#: How designs are placed on the reference chain (ProjectSpec.placement, ``--placement``) and the fit of target_frame each one stands for.
PLACEMENTS = {"interface": "interface", "whole": "all", "core": "core"}
DEFAULT_PROFILE = "default"
BUILTIN_PROFILES = ("default",)

ProjectKey = tuple[str, str]  # (run id, design id as the run wrote it)


class PipelineError(RuntimeError):
    """The inputs do not fit together or a stage failed; the message says what to look at."""


# ------------------------------------------------------------------------------------------------------------
# log with timings
# ------------------------------------------------------------------------------------------------------------
class Log:
    """Progress lines with elapsed seconds: ``[ 12.3s] stage ... done in 4.5 s``. ``quiet`` keeps only the timings."""

    def __init__(self, stream=None, quiet: bool = False):
        self.stream = stream if stream is not None else sys.stderr
        self.quiet = quiet
        self.t0 = time.perf_counter()
        self.timings: dict[str, float] = {}
        self._last_pct: dict[str, int] = {}

    def __call__(self, msg: str) -> None:
        if not self.quiet:
            print(f"[{time.perf_counter() - self.t0:6.1f}s] {msg}", file=self.stream, flush=True)

    @contextlib.contextmanager
    def stage(self, name: str) -> Iterator[None]:
        self(f"{name} ...")
        t = time.perf_counter()
        try:
            yield
        finally:
            dt = time.perf_counter() - t
            self.timings[name] = round(self.timings.get(name, 0.0) + dt, 2)
            self(f"{name} done in {dt:.1f} s")

    def progress(self, name: str, step: int = 25) -> Callable[[int, int], None]:
        """A ``progress(done, total)`` callback that prints at most one line per ``step`` percent."""

        def cb(done: int, total: int) -> None:
            pct = int(100 * done / total) if total else 100
            if pct // step > self._last_pct.get(name, -1) // step or done == total:
                if self._last_pct.get(name) != pct:
                    self(f"  {name}: {done}/{total}")
                self._last_pct[name] = pct

        self._last_pct[name] = -1
        return cb


# ------------------------------------------------------------------------------------------------------------
# hashing
# ------------------------------------------------------------------------------------------------------------
def sha256_file(path: str | os.PathLike) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_json(obj: Any) -> str:
    return sha256_text(json.dumps(obj, sort_keys=True, separators=(",", ":"), default=_json_default))


def _json_default(o: Any) -> Any:
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    if isinstance(o, Path):
        return str(o)
    if dataclasses.is_dataclass(o) and not isinstance(o, type):
        return dataclasses.asdict(o)
    raise TypeError(f"not JSON serialisable: {type(o)}")


@functools.lru_cache(maxsize=256)
def code_hash(*modules: Any) -> str:
    """Hash of the source files of modules (or of functions), so a code change can never serve a stale cache.
    Memoised per call signature: the sources do not change while a process runs."""
    h = hashlib.sha256()
    for m in modules:
        try:
            src = _inspect.getsource(m)
        except (OSError, TypeError):
            src = repr(m)
        h.update(src.encode("utf-8"))
    return h.hexdigest()[:16]


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def _read_cached(path: Path, key: str) -> Optional[Any]:
    """Payload of a cache file ``{"key", "data"}`` when the key matches, else None (a missing or corrupt file is a miss)."""
    try:
        blob = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return blob.get("data") if isinstance(blob, dict) and blob.get("key") == key else None


def _write_cached(path: Path, key: str, data: Any) -> Any:
    """Store ``data`` and return it as it will read back (plain JSON types), so a cold and a warm run see identical objects."""
    text = json.dumps({"key": key, "data": data}, default=_json_default, separators=(",", ":"))
    _atomic_write(path, text)
    return json.loads(text)["data"]


# ------------------------------------------------------------------------------------------------------------
# specification and project
# ------------------------------------------------------------------------------------------------------------
@dataclass
class RunSpec:
    """One run folder: the path, the short id (``run01``), an optional label and an optional settings file."""

    path: Path
    run_id: str
    label: Optional[str] = None
    settings: Optional[Path] = None


@dataclass
class ProjectSpec:
    """Everything ``run_pipeline`` needs: the runs, the target files and the options. Paths are kept as given."""

    name: str
    runs: list[RunSpec]
    target_pdb: Path
    reference_cif: Optional[Path] = None   # None: no full-length reference; the construct is the frame and every check that needs the rest of the protein is off
    uniprot: Optional[str] = None       # UniProt accession of the target; None: the profile's target.uniprot, else none (no UniProt records)
    profile: Any = DEFAULT_PROFILE      # a built-in id, a path to a profile JSON, or a dict
    consurf: bool = False               # include the optional ConSurf-DB overlay (academic licence), off by default
    online: bool = False                # allow refdata to fetch missing public reference data by accession
    jobs: Optional[int] = None
    use_cache: bool = True
    command: str = ""                   # the command line that produced the project (manifest tool.command)
    reference_id: Optional[str] = None  # name of the reference structure in labels (default: the file stem of reference_cif)
    reference_chain: Optional[str] = None   # chain of the target in the reference structure; None = the chain that carries the construct's sequence
    partner_chains: Optional[Sequence[str] | str] = None   # None or "auto" = the other subunits (same protein, in the target's assembly), () or "none" = none, "all" = every other chain of the assembly, or the names of chains
    assembly: Optional[str] = None      # the assembly of the reference the subunits come from; None = the first one that holds the target chain
    context_hard: Optional[str] = None  # auto, hard or advisory: whether the Context check decides a verdict (None: the profile's checks.context.hard, else auto)
    placement: str = "interface"        # how designs are placed on the reference chain: interface (weighted around the binder), whole (every shared CA) or core


@dataclass
class Project:
    """The assembled analysis, in memory. Plain dicts from the cached stages, objects only where the exporter needs them."""

    spec: ProjectSpec
    work: Path
    runs: dict[str, RunData]
    frame: TF.TargetFrame
    profile: dict[str, Any]
    poses: dict[ProjectKey, dict[str, Any]]
    analyses: dict[ProjectKey, interface.DesignAnalysis]
    families: dict[str, Any]
    modes: dict[str, Any]
    rescue: dict[ProjectKey, dict[str, Any]]
    homolog: dict[str, Any]
    project_ids: dict[ProjectKey, str]
    digests: dict[str, Any]
    timings: dict[str, float] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    network_requests: int = 0

    @property
    def run_ids(self) -> list[str]:
        return list(self.runs)

    def designs(self, run_id: Optional[str] = None, *, status: Optional[str] = None, with_structure: bool = True) -> list[DesignRecord]:
        out = []
        for rid, run in self.runs.items():
            if run_id is not None and rid != run_id:
                continue
            for d in run.accepted() + run.rejected_with_structure():
                if status is not None and d.status != status:
                    continue
                if with_structure and not d.structure_path:
                    continue
                out.append(d)
        return out

    def summary(self) -> dict[str, Any]:
        runs = {}
        for rid, run in self.runs.items():
            runs[rid] = {"accepted": len(run.accepted()), "rejected": len(run.rejected()),
                         "rejected_with_structure": len(run.rejected_with_structure()),
                         "families": len(run.accepted_backbones()), "rows": len(run.designs)}
        return {"name": self.spec.name, "runs": runs, "designs_exported": len(self.poses),
                "timings": dict(self.timings), "warnings": list(self.warnings)}


# ------------------------------------------------------------------------------------------------------------
# profile
# ------------------------------------------------------------------------------------------------------------
def profiles_dir() -> Path:
    return Path(__file__).resolve().parent / "profiles"


def load_profile(profile: Any = DEFAULT_PROFILE) -> dict[str, Any]:
    """A built-in profile by id (``default``), a profile JSON by path, or a dict (copied)."""
    if isinstance(profile, Mapping):
        return json.loads(json.dumps(profile))
    if isinstance(profile, str) and profile in BUILTIN_PROFILES:
        path = profiles_dir() / f"{profile}.profile.json"
    else:
        path = Path(profile)
    if not path.is_file():
        raise PipelineError(f"profile {profile!r} not found (built-in profiles: {', '.join(BUILTIN_PROFILES)}; or give the path of a profile JSON)")
    return json.loads(path.read_text(encoding="utf-8"))


def effective_uniprot(spec: ProjectSpec, profile: Mapping[str, Any]) -> Optional[str]:
    """The accession the stages use: the project's own (``--uniprot``), else the profile's ``target.uniprot``, else None."""
    acc = (spec.uniprot or (profile.get("target") or {}).get("uniprot") or "").strip()
    return acc or None


# ------------------------------------------------------------------------------------------------------------
# stage 1: read the runs
# ------------------------------------------------------------------------------------------------------------
def read_runs(spec: ProjectSpec, log: Optional[Log] = None) -> dict[str, RunData]:
    """Read every run folder with the reader that detects it. BindCraft 1 runs are compared with the construct file."""
    log = log or Log(quiet=True)
    runs: dict[str, RunData] = {}
    for rs in spec.runs:
        if rs.run_id in runs:
            raise PipelineError(f"run id {rs.run_id!r} is used twice")
        path = Path(rs.path)
        if not path.is_dir():
            raise PipelineError(f"{path} is not a folder")
        reader = detect_reader(path)
        if reader is None:
            raise PipelineError(f"{path} is not a run folder of a known reader")
        opts: dict[str, Any] = {}
        if reader.NAME == "bindcraft1":
            opts["construct_path"] = spec.target_pdb
        t = time.perf_counter()
        run = reader.read_run(path, rs.run_id, label=rs.label, settings_path=rs.settings, **opts)
        log(f"  {rs.run_id}: {len(run.accepted())} accepted, {len(run.rejected())} rejected rows "
            f"({len(run.rejected_with_structure())} with a structure), {len(run.accepted_backbones())} families "
            f"[{reader.NAME}, {time.perf_counter() - t:.1f} s]")
        for w in run.warnings:
            log(f"  WARNING {rs.run_id}: {w}")
        runs[rs.run_id] = run
    return runs


def one_based_structures(runs: Mapping[str, RunData], work: Path, log: Optional[Log] = None) -> int:
    """Point the designs of a zero-based BindCraft 2 run at one-based copies of their structures; returns how many were copied.

    BindCraft 2 versions 1.0.0 and 1.0.1 write the residue numbers of the binder chains zero-based, while the exporter, the interface
    analysis and the inspector expect 1 to n (docs/bindcraft2_reader.md). The reader records which numbering a campaign uses
    (``annotations["binder_numbering"]``); for a zero-based run every design structure is rewritten by ``bindcraft2.one_based_copy`` under
    ``work/onebased/<run>/`` (named by a prefix of the original's sha256 and a tag of the shift and of the code that writes the copy, so a
    copy is reused only while all of them are the same) and the design's ``structure_path`` becomes the copy. The original files are
    inputs and are never touched; their path stays in ``design.extras["source_structure"]``. Runs of any other reader, and one-based
    campaigns, are left as they are."""
    from .readers import bindcraft2

    log = log or Log(quiet=True)
    total = 0
    for rid, run in runs.items():
        shift = int((run.annotations.get("binder_numbering") or {}).get("offset_to_one_based", 0))
        if run.reader != "bindcraft2" or not shift:
            continue
        sha = {f.path: f.sha256 for f in run.files if f.role.endswith("_structure")}
        tag = sha256_json({"shift": shift, "code": code_hash(bindcraft2.one_based_copy)})[:6]
        out_dir = work / "onebased" / rid
        out_dir.mkdir(parents=True, exist_ok=True)
        n = 0
        for d in run.designs:
            if not d.structure_path or "source_structure" in d.extras:        # no structure, or already pointed at its copy
                continue
            src = str(Path(d.structure_path).resolve())
            dst = out_dir / f"{(sha.get(src) or sha256_file(src))[:12]}_{tag}_{Path(src).stem}.cif"
            if not dst.is_file():
                tmp = dst.with_name(dst.name + ".tmp")
                bindcraft2.one_based_copy(run, d, tmp)
                tmp.replace(dst)
            d.extras["source_structure"] = src
            d.structure_path = str(dst)
            n += 1
        if n:
            log(f"  {rid}: binder numbering is zero-based, {n} structures renumbered to one-based copies")
            note = ("The binder residue numbers of this campaign start at 0; the analysis and the exported structures use copies of the structure files "
                    "renumbered from 1 (the binder chains only), and provenance names the original files, which are not modified.")
            if note not in run.notes:
                run.notes.append(note)
        total += n
    return total


def structure_hashes(runs: Mapping[str, RunData]) -> dict[str, str]:
    """{resolved structure path: sha256} for every structure the readers hashed."""
    out = {}
    for run in runs.values():
        for f in run.files:
            if f.role in ("accepted_structure", "rejected_structure"):
                out[f.path] = f.sha256
    return out


def _partner_key(partners: Any) -> Any:
    """The partner request as it enters a cache key: None and "auto" are the same request, a string such as "all" stays one word."""
    if partners is None or (isinstance(partners, str) and partners.strip().lower() in ("", "auto")):
        return None
    return partners if isinstance(partners, str) else list(partners)


def input_digests(spec: ProjectSpec, runs: Mapping[str, RunData]) -> dict[str, Any]:
    """Hashes of every input and of the run definitions: the root of every cache key."""
    per_run = {}
    for rid, run in runs.items():
        files = sorted((f.role, Path(f.path).name, f.sha256) for f in run.files)
        per_run[rid] = {
            "files": sha256_json(files),
            "definition": sha256_json({
                "hotspots": run.settings.hotspots, "numbers": run.target.numbers, "label": run.label,
                "length": [run.settings.length_min, run.settings.length_max], "reader": run.reader}),
        }
    return {
        "pipeline": PIPELINE_VERSION, "version": __version__,
        "reference_cif": sha256_file(spec.reference_cif) if spec.reference_cif else None, "target_pdb": sha256_file(spec.target_pdb),
        "uniprot": spec.uniprot, "reference": [spec.reference_id, spec.reference_chain, _partner_key(spec.partner_chains), spec.assembly],
        "runs": per_run,
    }


# ------------------------------------------------------------------------------------------------------------
# stage 2: target frame
# ------------------------------------------------------------------------------------------------------------
def frame_key(spec: ProjectSpec, runs: Mapping[str, RunData], digests: Mapping[str, Any]) -> str:
    return sha256_json({
        "ref": digests["reference_cif"], "construct": digests["target_pdb"], "uniprot": spec.uniprot,
        "reference": digests["reference"], "placement": spec.placement,
        "runs": {rid: {"hotspots": r.settings.hotspots, "numbers": r.target.numbers, "label": r.label or rid,
                       "sequence": r.target.sequence, "consistent": r.target.consistent} for rid, r in runs.items()},
        "code": code_hash(TF, S), "refdata": _refdata_digest(), "v": PIPELINE_VERSION,
    })


def _refdata_digest() -> str:
    try:
        entries = rd.manifest().get("entries", {})
    except Exception:  # noqa: BLE001 - a missing manifest only widens the key
        return "none"
    return sha256_json({k: v.get("sha256") for k, v in sorted(entries.items())})


def build_target_frame(spec: ProjectSpec, runs: Mapping[str, RunData], digests: Mapping[str, Any], work: Path,
                       log: Optional[Log] = None) -> TF.TargetFrame:
    """The target frame, loaded from ``work/target`` when its key matches, built and saved otherwise. Always assigns the
    construct-to-UniProt offset to every run (the readers leave it empty)."""
    log = log or Log(quiet=True)
    d = work / "target"
    key = frame_key(spec, runs, digests)
    frame: Optional[TF.TargetFrame] = None
    kf = d / "frame.key"
    if spec.use_cache and kf.is_file() and kf.read_text().strip() == key:
        try:
            frame = TF.TargetFrame.load(d)
            log("  target frame: cache hit")
        except TF.FrameError:
            frame = None
    if frame is None:
        frame = TF.build_frame(spec.reference_cif, spec.target_pdb, dict(runs), spec.uniprot, chain=spec.reference_chain,
                               partner_chains=spec.partner_chains, assembly_id=spec.assembly, reference_id=spec.reference_id,
                               config={"geometry": {"fit": PLACEMENTS[spec.placement]}})
        frame.save(d)
        _atomic_write(kf, key + "\n")
        for w in frame.meta.get("warnings", []):
            log(f"  WARNING target frame: {w}")
        frame = TF.TargetFrame.load(d)      # a cold and a warm build read the same file, so their blobs have the same key order
    frame.assign(runs)
    if frame.composition:               # which chain is the target, which are partners and why: said in the log, never left to be guessed
        for line in TF.format_reference_chains(frame.composition):
            log(f"  reference chains: {line}")
    return frame


# ------------------------------------------------------------------------------------------------------------
# stage 3: pose, context and thumbnail view of every design
# ------------------------------------------------------------------------------------------------------------
def canonical_view(tgt_xyz: np.ndarray, bnd_xyz: np.ndarray, cutoff: float = model.CONTACT_CUTOFF) -> tuple[np.ndarray, np.ndarray]:
    """Thumbnail view of a placed design: rows are the view x, y and z axes in the reference frame, and the interface centre.

    The mockup's recipe (mockups/build/real/rdesign.py): the centre is the mean midpoint of the target-binder atom pairs
    within the contact cut-off, x points from it to the binder centroid, y is the first principal axis of the two chains
    projected normal to x, z completes the frame."""
    from scipy.spatial import cKDTree

    tt, bt = cKDTree(tgt_xyz), cKDTree(bnd_xyz)
    pairs = tt.query_ball_tree(bt, cutoff)
    ia = np.fromiter((i for i, lst in enumerate(pairs) for _ in lst), dtype=np.int64)
    jb = np.fromiter((j for lst in pairs for j in lst), dtype=np.int64)
    if len(ia):
        ic = (0.5 * (tgt_xyz[ia] + bnd_xyz[jb])).mean(axis=0)
    else:
        ic = tgt_xyz.mean(axis=0) + 0.5 * (bnd_xyz.mean(axis=0) - tgt_xyz.mean(axis=0))
    a = bnd_xyz.mean(axis=0) - ic
    if np.linalg.norm(a) < 1e-6:
        a = bnd_xyz.mean(axis=0) - tgt_xyz.mean(axis=0)
    a = a / np.linalg.norm(a)
    pts = np.vstack([tgt_xyz[::4], bnd_xyz[::4]]) - ic
    pp = pts - np.outer(pts @ a, a)
    _, _, vt = np.linalg.svd(pp - pp.mean(axis=0), full_matrices=False)
    y = vt[0] - np.dot(vt[0], a) * a
    y = y / np.linalg.norm(y)
    return np.stack([a, y, np.cross(a, y)]), ic


def pose_one(frame: TF.TargetFrame, cx: S.Complex, run_id: str, design_id: str, status: str, rank: Optional[int],
             model_k: Optional[int]) -> dict[str, Any]:
    """Placement on the reference chain frame, the context checks and the thumbnail view of one loaded design.

    Everything is in the reference frame and UniProt numbered; the record is small plain JSON (docs/export.md).
    ``metrics`` holds full-precision values, the exporter rounds."""
    chk = TF.check_design(frame, cx, run_id, design_id=design_id, status=status, rank=rank, model_k=model_k)
    pl, res = chk.placement, chk.context
    placed_b = pl.apply(cx.binder)
    placed_t = pl.apply(cx.target)
    ca = placed_b.ca_xyz()
    if np.isnan(ca).any():
        raise PipelineError(f"{design_id}: a binder residue has no CA atom")
    view, ic = canonical_view(placed_t.heavy().xyz, placed_b.xyz)
    diag = res.diag
    return {
        "run": run_id, "id": design_id, "n_binder": int(len(ca)),
        "placement": {"R": [round(float(x), 6) for x in pl.R.reshape(-1)], "t": [round(float(x), 4) for x in pl.t],
                      "rmsd": round(float(pl.rmsd), 3), "n": int(pl.n), "core_rmsd": round(float(pl.core_rmsd), 3),
                      "core_n": int(pl.core_n), "dropped": [int(x) for x in pl.dropped_numbers()],
                      "max_dev": round(float(pl.max_deviation), 2), "n_mismatch": int(pl.n_mismatch),
                      "warnings": list(pl.warnings), **pl.fit_fields()},
        "metrics": {k: (None if v is None else (int(v) if isinstance(v, (int, np.integer)) else float(v)))
                    for k, v in res.metrics.items()},
        "glycan": {k: float(v) for k, v in res.glycan.items()}, "glycan_min": res.glycan_min,
        "ctx": json.loads(json.dumps(res.ctx, default=_json_default)), "octx": json.loads(json.dumps(res.octx, default=_json_default)),
        "diag": {k: json.loads(json.dumps(diag.get(k), default=_json_default))
                 for k in ("partner_pair", "rest_pair", "partner_clash", "rest_clash", "rest_n_by_region", "dimer_buried_contacts",
                           "cut_contacts", "covered", "buried_in_assembly", "partner_role", "partner_roles") if k in diag},
        "contacted": [int(x) for x in chk.contacted], "hot": [int(x) for x in chk.hotspots_touched],
        "hot_distance": {str(k): round(float(v), 2) for k, v in chk.hotspot_distance.items()},
        "contact_dev_max": None if chk.contact_deviation_max is None else round(float(chk.contact_deviation_max), 2),
        "termini": {k: float(v) for k, v in chk.termini.items()},
        "ca": [round(float(x), 2) for x in ca.reshape(-1)],
        "view": [round(float(x), 4) for x in view.reshape(-1)], "ic": [round(float(x), 2) for x in ic],
    }


_FRAMES: dict[tuple[str, str], TF.TargetFrame] = {}


def _frame_for(frame_dir: str, fkey: str = "") -> TF.TargetFrame:
    """The frame saved in ``frame_dir``, read once per process and frame key (a folder that is built again with other settings is read again)."""
    f = _FRAMES.get((frame_dir, fkey))
    if f is None:
        f = _FRAMES[(frame_dir, fkey)] = TF.TargetFrame.load(frame_dir)
    return f


def _pose_worker(args: tuple) -> tuple[str, dict[str, Any]]:
    frame_dir, run_id, did, path, status, rank, model_k, fkey, n_chains = args
    cx = S.load_complex(path, n_binder_chains=n_chains)
    return did, pose_one(_frame_for(frame_dir, fkey), cx, run_id, did, status, rank, model_k)


def pose_key(structure_sha: str, fkey: str, n_binder_chains: int = 1) -> str:
    # The chain count is in the key: the same file read with one binder chain and
    # with two gives different poses, and the file's hash cannot say which.
    return sha256_json({"s": structure_sha, "frame": fkey, "nb": int(n_binder_chains),
                        "code": code_hash(TF, S, pose_one, canonical_view), "v": POSE_VERSION})


def compute_poses(runs: Mapping[str, RunData], frame: TF.TargetFrame, fkey: str, work: Path, *, jobs: Optional[int] = None,
                  use_cache: bool = True, log: Optional[Log] = None) -> dict[ProjectKey, dict[str, Any]]:
    """Pose records of every accepted design and every rejected design with a structure, cached per design in
    ``work/poses/<run>/<id>.json``. The frame is read back from ``work/target`` by the workers."""
    log = log or Log(quiet=True)
    hashes = structure_hashes(runs)
    out: dict[ProjectKey, dict[str, Any]] = {}
    todo: list[tuple] = []
    keys: dict[ProjectKey, str] = {}
    for rid, run in runs.items():
        for d in run.accepted() + run.rejected_with_structure():
            sp = str(Path(d.structure_path).resolve())
            sha = hashes.get(sp) or hashes.get(d.structure_path) or sha256_file(d.structure_path)
            k = pose_key(sha, fkey, binder_chain_count(d))
            keys[(rid, d.id)] = k
            hit = _read_cached(work / "poses" / rid / f"{d.id}.json", k) if use_cache else None
            if hit is not None:
                out[(rid, d.id)] = hit
            else:
                todo.append((str(work / "target"), rid, d.id, d.structure_path, d.status, d.rank, d.structure_model, fkey,
                             binder_chain_count(d)))
    total = len(keys)
    cb = log.progress("poses", 25)
    if out:
        cb(len(out), total)
    workers = jobs if jobs is not None else _default_jobs()

    def store(rid: str, did: str, rec: dict[str, Any]) -> None:
        out[(rid, did)] = _write_cached(work / "poses" / rid / f"{did}.json", keys[(rid, did)], rec)
        cb(len(out), total)

    if todo and (workers <= 1 or len(todo) <= 2):
        for args in todo:
            did, rec = _pose_worker(args)
            store(args[1], did, rec)
    elif todo:
        ctx = mp.get_context("spawn")
        with cf.ProcessPoolExecutor(max_workers=min(workers, len(todo)), mp_context=ctx) as ex:
            futs = {ex.submit(_pose_worker, a): a for a in todo}
            for fut in cf.as_completed(futs):
                a = futs[fut]
                did, rec = fut.result()
                store(a[1], did, rec)
    return out


def _default_jobs() -> int:
    return max(1, min(8, (os.cpu_count() or 2) - 1))


# ------------------------------------------------------------------------------------------------------------
# stage 4: interface analysis and liabilities (interface.analyze_run caches per design on its own)
# ------------------------------------------------------------------------------------------------------------
def compute_interface(runs: Mapping[str, RunData], work: Path, *, jobs: Optional[int] = None, use_cache: bool = True,
                      log: Optional[Log] = None) -> dict[ProjectKey, interface.DesignAnalysis]:
    """Interface analysis of every accepted design and every rejected design with a structure (construct numbering,
    liabilities included). Delegates the process pool and the disk cache (``work/analysis``) to ``interface.analyze_run``."""
    log = log or Log(quiet=True)
    designs = [d for run in runs.values() for d in run.accepted() + run.rejected_with_structure()]
    ids = [d.id for d in designs]
    clash = len(set(ids)) != len(ids)
    hotspots = {rid: list(run.settings.hotspots) for rid, run in runs.items()}
    res = interface.analyze_run(designs, n_jobs=jobs if jobs is not None else _default_jobs(),
                                cache_dir=(work / "analysis") if use_cache else None, hotspots=hotspots,
                                progress=log.progress("interface", 25))
    out: dict[ProjectKey, interface.DesignAnalysis] = {}
    for d in designs:
        out[(d.run_id, d.id)] = res[f"{d.run_id}/{d.id}" if clash else d.id]
    return out


# ------------------------------------------------------------------------------------------------------------
# stage 5: families, diversity and the fold matrix
# ------------------------------------------------------------------------------------------------------------
def compute_families(runs: Mapping[str, RunData], digests: Mapping[str, Any], work: Path, *, jobs: Optional[int] = None,
                     use_cache: bool = True, log: Optional[Log] = None) -> dict[str, Any]:
    """Families, the pooled fold matrix and labels, per-run diversity and the fold label of every rejected structure,
    as one plain dict (``families.analyze_families`` serialised); cached in ``work/stages/families.json``."""
    log = log or Log(quiet=True)
    key = sha256_json({"d": digests, "code": code_hash(families, S), "v": PIPELINE_VERSION,
                       "offsets": {rid: r.target.uniprot_offset for rid, r in runs.items()}})
    path = work / "stages" / "families.json"
    hit = _read_cached(path, key) if use_cache else None
    if hit is not None:
        log("  families: cache hit")
        return hit
    fa = families.analyze_families(list(runs.values()), n_jobs=jobs if jobs is not None else _default_jobs(), rejected_labels=True)
    data = {
        "schema": families.SCHEMA, "runs": fa.runs, "fold": fa.fold, "fold_labels": fa.fold_labels, "diversity": fa.diversity,
        "rejected_fold": fa.rejected_fold, "chronology": fa.chronology, "warnings": fa.warnings, "seconds": fa.seconds,
        "families": {rid: [f.to_dict() for f in fams] for rid, fams in fa.families.items()},
        "design_ids": [[r, i, pid] for (r, i), pid in sorted(fa.design_ids.items())],
    }
    return _write_cached(path, key, data)


# ------------------------------------------------------------------------------------------------------------
# stage 6: modes, the map embedding and near-miss records
# ------------------------------------------------------------------------------------------------------------
def compute_modes(runs: Mapping[str, RunData], frame: TF.TargetFrame, poses: Mapping[ProjectKey, Mapping[str, Any]],
                  pids: Mapping[ProjectKey, str], digests: Mapping[str, Any], work: Path, *, use_cache: bool = True,
                  log: Optional[Log] = None) -> tuple[dict[str, Any], dict[ProjectKey, dict[str, Any]]]:
    """Binding modes, per-design mode, site and map position (``modes.build_modes``) and the near-miss records of the
    rejected designs (``modes.rescue_records``). Footprints are the contact sets of the pose stage (UniProt numbers, 4.5 A heavy
    atom, the page definition), so no structure is read again. Cached in ``work/stages/modes.json``."""
    log = log or Log(quiet=True)
    key = sha256_json({"d": digests, "code": code_hash(modes, config, model), "v": PIPELINE_VERSION,
                       "poses": sha256_json({f"{k[0]}/{k[1]}": v["contacted"] for k, v in sorted(poses.items())}),
                       "frame": _frame_signature(frame, site_residues(runs, poses))})
    path = work / "stages" / "modes.json"
    hit = _read_cached(path, key) if use_cache else None
    if hit is None:
        foot: dict[str, list[int]] = {}
        infos: dict[str, modes.RunInfo] = {}
        for rid, run in runs.items():
            info = modes.RunInfo(rid, hotspots=run.hotspots_uniprot(), residues=set(run.target.uniprot_numbers()), label=run.label)
            for d in run.accepted() + run.rejected_with_structure():
                pid = pids[(rid, d.id)]
                foot[pid] = list(poses[(rid, d.id)]["contacted"])
                info.designs[pid] = d.status
            infos[rid] = info
        tb = frame.target_block(list(runs), site_residues=site_residues(runs, poses))
        res = modes.build_modes(foot, None, infos, residue_names=modes.residue_names_from_runs(runs), ref_ca=frame.ref_ca,
                                sites=tb["sites"])
        rescue = {pids[(rid, did)]: rec for rid, run in runs.items() for did, rec in modes.rescue_records(run).items()}
        summary = {rid: modes.rescue_summary(run) for rid, run in runs.items()}
        hit = _write_cached(path, key, {"modes": res.to_dict(), "rescue": rescue, "rescue_summary": summary})
    else:
        log("  modes: cache hit")
    by_key = {(rid, did): hit["rescue"][pids[(rid, did)]] for rid, run in runs.items()
              for did in (d.id for d in run.rejected_with_structure())}
    return hit["modes"], by_key


def site_residues(runs: Mapping[str, RunData], poses: Mapping[ProjectKey, Mapping[str, Any]], share: float = 0.5) -> list[int]:
    """The residues the designs contact most, UniProt numbers: the contacted residues (pose records) that at least ``share`` of the
    most contacted residue's count of the accepted designs contact (of all designs with a pose when none was accepted). Site A1 of the
    target block is built from them when no run records a hotspot; with a hotspot anywhere the result is not used."""
    keys = [k for k in poses if runs[k[0]].design(k[1]).status == "accepted"] or list(poses)
    counts: dict[int, int] = {}
    for k in keys:
        for r in poses[k]["contacted"]:
            counts[int(r)] = counts.get(int(r), 0) + 1
    if not counts:
        return []
    top = max(counts.values())
    return sorted(r for r, c in counts.items() if c >= share * top)


def _frame_signature(frame: TF.TargetFrame, site: Sequence[int] = ()) -> str:
    return sha256_json({"offset": frame.offset, "runs": {r: frame.run_info(r)["hotspots"] for r in frame.run_ids},
                        "sites": frame.target_block(site_residues=site)["sites"]})


# ------------------------------------------------------------------------------------------------------------
# stage 7: specificity panel
# ------------------------------------------------------------------------------------------------------------
def compute_specificity(runs: Mapping[str, RunData], profile: Mapping[str, Any], work: Path, *, include_consurf: bool = False,
                        use_cache: bool = True, log: Optional[Log] = None, accession: Optional[str] = None,
                        frame: Optional[TF.TargetFrame] = None) -> dict[str, Any]:
    """The HOM blob over the union of the runs' constructs (``specificity.build_panel``), cached in ``work/stages/homolog.json``.

    A panel is built when the target's UniProt accession is known and the profile names species or paralogs; otherwise (no accession,
    or nothing to compare with) the blob is empty (no orthologs, no paralogs) and the pages hide the specificity features. Only
    accessions and gene names are ever used; the ConSurf-DB overlay is included only when asked for."""
    log = log or Log(quiet=True)
    resis = model.union_uniprot_residues([r.target for r in runs.values()])
    want = specificity.panel_requested(profile)
    ref = ({"pdb": frame.meta["reference"]["pdb"], "chain": frame.meta["chain"],
            "observed": [frame.meta["numbering"]["reference_first"], frame.meta["numbering"]["reference_last"]]} if frame is not None else {})
    key = sha256_json({"resis": resis, "spec": profile.get("specificity"), "acc": accession, "ref": ref, "consurf": include_consurf,
                       "refdata": _refdata_digest(), "code": code_hash(specificity), "v": PIPELINE_VERSION})
    path = work / "stages" / "homolog.json"
    hit = _read_cached(path, key) if use_cache else None
    if hit is not None:
        log("  specificity: cache hit")
        return hit
    if not (accession and want):
        why = ("no UniProt accession was given" if want and not accession else
               "the profile names no species and no paralogs")
        log(f"  specificity: skipped ({why})")
        letters = {int(u): frame.residue(u) for u in resis} if frame is not None else modes.residue_names_from_runs(runs)
        sp = profile.get("specificity") or {}
        blob = specificity.empty_homolog_json(resis, letters, subject=sp.get("subject") or "", kind=sp.get("kind") or "host",
                                              note=f"No specificity panel: {why}.")
        return _write_cached(path, key, blob)
    panel = specificity.build_panel(resis, profile, include_consurf=include_consurf, accession=accession,
                                    observed=ref.get("observed"), pdb=ref.get("pdb"), chain=ref.get("chain", "A"))
    return _write_cached(path, key, specificity.to_homolog_json(panel))


# ------------------------------------------------------------------------------------------------------------
# the whole analysis
# ------------------------------------------------------------------------------------------------------------
def design_ids(runs: Mapping[str, RunData]) -> dict[ProjectKey, str]:
    """Project design ids: the id itself when it is unique across the project, ``<run>/<id>`` otherwise."""
    return families.project_design_ids(list(runs.values()))


@contextlib.contextmanager
def _offline(enforce: bool) -> Iterator[None]:
    """Keep the reference data off the network (the public data are cached by accession) unless the build was told it may fetch."""
    old = os.environ.get(rd.OFFLINE_ENV)
    if enforce:
        os.environ[rd.OFFLINE_ENV] = "1"
    try:
        yield
    finally:
        if enforce:
            if old is None:
                os.environ.pop(rd.OFFLINE_ENV, None)
            else:
                os.environ[rd.OFFLINE_ENV] = old


def run_pipeline(spec: ProjectSpec, work: str | os.PathLike, *, log: Optional[Log] = None) -> Project:
    """Run every stage, with the disk cache under ``work``, and return the assembled Project.

    Order (docs/CONTRACTS.md 5.8): read runs, target frame, poses and context checks, interface analysis and
    liabilities, families and diversity, modes and near-miss records, specificity panel. A rerun with a warm cache reads
    JSON back and takes seconds (the readers and the families stage dominate it)."""
    log = log or Log()
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    for p in (spec.target_pdb, spec.reference_cif):
        if p is not None and not Path(p).is_file():
            raise PipelineError(f"{p} not found")
    if spec.placement not in PLACEMENTS:
        raise PipelineError(f"placement must be one of {', '.join(PLACEMENTS)}, not {spec.placement!r}")
    if not spec.use_cache:                      # recompute everything: drop the stage caches, they are rewritten below
        drop_caches(work)
    n_before = _refdata_entries()
    with _offline(not spec.online):
        with log.stage("read runs"):
            runs = read_runs(spec, log)
            one_based_structures(runs, work, log)
        profile = load_profile(spec.profile)
        spec = dataclasses.replace(spec, uniprot=effective_uniprot(spec, profile))   # the stages see the accession the profile may give
        digests = input_digests(spec, runs)
        with log.stage("target frame"):
            frame = build_target_frame(spec, runs, digests, work, log)
        fkey = frame_key(spec, runs, digests)
        with log.stage("poses and context checks"):
            poses = compute_poses(runs, frame, fkey, work, jobs=spec.jobs, use_cache=True, log=log)
        with log.stage("interface analysis"):
            analyses = compute_interface(runs, work, jobs=spec.jobs, use_cache=True, log=log)
        with log.stage("families and diversity"):
            fams = compute_families(runs, digests, work, jobs=spec.jobs, use_cache=True, log=log)
        pids = design_ids(runs)
        with log.stage("modes and near-miss records"):
            modes_block, rescue = compute_modes(runs, frame, poses, pids, digests, work, use_cache=True, log=log)
        with log.stage("specificity panel"):
            homolog = compute_specificity(runs, profile, work, include_consurf=spec.consurf, use_cache=True, log=log,
                                          accession=spec.uniprot, frame=frame)
    warnings: list[str] = []
    for rid, run in runs.items():
        warnings += [f"{rid}: {w}" for w in run.warnings]
    warnings += [f"target frame: {w}" for w in frame.meta.get("warnings", [])]
    warnings += [f"families: {w}" for w in fams.get("warnings", [])]
    if specificity.panel_requested(profile) and not spec.uniprot:
        warnings.append("specificity: the profile names species or paralogs but no UniProt accession was given (--uniprot, or target.uniprot "
                        "in the profile), so no panel was built")
    project = Project(spec=spec, work=work, runs=runs, frame=frame, profile=profile, poses=poses, analyses=analyses,
                      families=fams, modes=modes_block, rescue=rescue, homolog=homolog, project_ids=pids, digests=digests,
                      timings=dict(log.timings), warnings=warnings, network_requests=max(0, _refdata_entries() - n_before))
    _write_runrecords(project)
    return project


def drop_caches(work: Path) -> int:
    """Delete the stage cache files (``*.json`` under analysis/, poses/ and stages/, the frame key) of a project folder and return how
    many were removed. Only files of that shape are touched, never a whole folder, so a project folder that holds anything else is safe."""
    n = 0
    for sub in ("analysis", "poses", "stages"):
        d = work / sub
        if not d.is_dir():
            continue
        for p in list(d.rglob("*.json")) + list(d.rglob("*.tmp")):
            p.unlink(missing_ok=True)
            n += 1
        for q in sorted((x for x in d.rglob("*") if x.is_dir()), reverse=True):
            try:
                q.rmdir()                       # only empty folders go
            except OSError:
                pass
        try:
            d.rmdir()
        except OSError:
            pass
    kf = work / "target" / "frame.key"
    if kf.is_file():
        kf.unlink()
        n += 1
    return n


def _refdata_entries() -> int:
    try:
        return len(rd.manifest().get("entries", {}))
    except Exception:  # noqa: BLE001
        return 0


def _write_runrecords(project: Project) -> None:
    """work/<project>/runs/<run>/runrecord.json: ``RunData.to_dict()`` for audit (no sequences, CONTRACTS section 6)."""
    given = project.frame.reference_given
    for rid, run in project.runs.items():
        rec = run.to_dict()
        # which situation the context checks were in, so an audit of the run does not have to infer it from absent values
        rec["target_frame"] = {"reference_given": given, "reference": Path(project.spec.reference_cif).name if given else None,
                               "context_checks": "run against the reference structure" if given else
                               "not run: no reference structure was given, so nothing outside the construct is known"}
        _atomic_write(project.work / "runs" / rid / "runrecord.json", json.dumps(rec, indent=1, default=_json_default, sort_keys=True) + "\n")


# ------------------------------------------------------------------------------------------------------------
# helpers shared with the exporter
# ------------------------------------------------------------------------------------------------------------
def utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def default_run_id(folder: Path) -> str:
    """Short run id from a folder name: campaign_run01 -> run01, otherwise the folder name."""
    import re

    m = re.search(r"(run\d+)", Path(folder).name, re.IGNORECASE)
    return m.group(1).lower() if m else Path(folder).name
