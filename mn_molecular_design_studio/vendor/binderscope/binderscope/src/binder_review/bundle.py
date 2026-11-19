"""Write the static bundle: the pages of viewer/ (through ``viewer_build.render_site``), the data blobs as JSON and as same-folder
scripts, the structures, the vendored 3Dmol.js and the review folder. The bundle opens from disk (file://) with no server and
works offline and has no network reference (system fonts; ``fonts="google"`` adds the optional Google Fonts link).

Layout (docs/export.md):

    index.html board.html inspector.html map.html query.html [guide.html]
    assets/css/shared.css  assets/js/*.js  assets/3Dmol-min.js  assets/3Dmol-LICENSE.txt  assets/binderscope-mark.svg  [assets/guide/*]
    data/{campaign,inspect,homolog,profile,manifest}.json   canonical copies for tools and tests
    data/{campaign,inspect,homolog,profile,manifest}.js     the same objects as scripts the pages load (window.RAW, INS, HOM, PROFILE, MANIFEST);
                                                            data/inspect.js leaves the binder structures (pdb_b) out: the inspector loads
                                                            structures/<id>.js on demand (``embed_structures`` keeps them in the script)
    structures/<design id>.pdb                              the design's own target and the binder, in the reference chain frame
    structures/<design id>.js                               the same text as a lazy-loaded script (window.BINDERSCOPE.structures[<id>])
    reviews/                                                where the user keeps saved review files

``write_bundle`` reads the export folder written by ``export.write_export`` (``work/<project>/export``), so ``build-bundle``
is fast and can be repeated whenever the viewer templates change.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Optional

from . import viewer_build
from .export import dumps, _load_validator

BLOB_NAMES = ("campaign", "inspect", "homolog", "profile", "manifest")
README = """Binderscope bundle
==================

Open index.html in a browser (Chrome, Edge, Firefox or Safari). Everything runs from this folder, offline, with no server:
the pages load their data from data/*.js, the structures are in structures/, and 3Dmol.js is in assets/.

  index.html      overview of the runs, what stands out, the build record
  board.html      family board: tiles, verdicts, diversity, decisions, order list
  inspector.html  one design in depth, findings linked to the 3D view
  map.html        binding modes, hotspot coverage, species and paralog specificity
  query.html      rule stack with a funnel, near-miss rescue, plain-words query
  data/           campaign, inspect, homolog, profile and manifest (JSON for tools, JS for the pages)
  structures/     the design models placed on the reference chain (PDB)
  reviews/        put saved review files here (the board's "Save review file" downloads one)

Scores are screens, not verdicts. data/manifest.json records the tool version, the input hashes, the thresholds and the checks
made at build time. It records file and folder names only (the absolute paths of the machine that built the bundle are cut off) but copies the
settings of each run: read it before sharing the bundle outside your lab.
"""


class BundleError(RuntimeError):
    """The bundle cannot be written (missing export, a page that points at a missing file)."""


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_export(export_dir: str | Path) -> dict[str, Any]:
    """The blobs of ``work/<project>/export`` as objects, plus the list of structure files."""
    d = Path(export_dir)
    missing = [n for n in ("campaign", "inspect", "homolog", "profile", "manifest") if not (d / f"{n}.json").is_file()]
    if missing:
        raise BundleError(f"{d} is not an export folder: {', '.join(m + '.json' for m in missing)} missing "
                          "(run 'binder-review import' first)")
    out: dict[str, Any] = {n: json.loads((d / f"{n}.json").read_text(encoding="utf-8")) for n in BLOB_NAMES}
    out["structures"] = sorted((d / "structures").glob("*.pdb")) if (d / "structures").is_dir() else []
    return out


def structure_script(design_id: str, text: str) -> str:
    """A lazy-loaded structure: ``<script src="structures/<id>.js">`` defines ``window.BINDERSCOPE.structures[<design id>]`` as the PDB text.
    A JSON string literal is a valid JavaScript string, so the text needs no further escaping."""
    return ("window.BINDERSCOPE = window.BINDERSCOPE || {};\n(window.BINDERSCOPE.structures = window.BINDERSCOPE.structures || {})["
            + json.dumps(design_id) + "] = " + json.dumps(text) + ";\n")


def design_id_of(file_name: str) -> str:
    """The design id of a structure file name (the exporter writes '/' of run-prefixed ids as '__')."""
    stem = file_name[:-4] if file_name.endswith((".pdb", ".js")) and file_name[-4] == "." else Path(file_name).stem
    return stem.replace("__", "/")


def _list_files(root: Path, skip: set[str]) -> list[dict[str, Any]]:
    files = []
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root).as_posix()
        if rel in skip or p.name == ".DS_Store":
            continue
        files.append({"path": rel, "bytes": p.stat().st_size, "sha256": sha256_file(p)})
    return files


def _slim_inspect_script(out: Path) -> int:
    """Rewrite data/inspect.js without ``pdb_b`` of any design; returns how many were left out. The inspector reads a design's
    binder from structures/<id>.js when its record has none (the file the bundle already holds for every design), so the page
    opens faster and the script stays small however many designs the bundle has. data/inspect.json is left as it is."""
    ins = json.loads((out / "data" / "inspect.json").read_text(encoding="utf-8"))
    n = 0
    for d in (ins.get("designs") or {}).values():
        if d.pop("pdb_b", None) is not None:
            n += 1
    (out / "data" / "inspect.js").write_text(viewer_build.js_wrapper("inspect", dumps(ins)), encoding="utf-8")
    return n


def check_clean_target(out: Path) -> None:
    """Refuse to wipe a folder that is not an earlier Binderscope bundle: ``write_bundle`` cleans its output folder, and a
    mistyped ``--out`` must never delete someone's files. An empty or missing folder and a previous bundle (it has index.html and
    data/manifest.json) are fine."""
    if not out.exists():
        return
    if not out.is_dir():
        raise BundleError(f"{out} exists and is not a folder")
    entries = [p for p in out.iterdir() if p.name != ".DS_Store"]
    if entries and not ((out / "index.html").is_file() and (out / "data" / "manifest.json").is_file()):
        raise BundleError(f"{out} exists, is not empty and is not a Binderscope bundle: refusing to clean it "
                          "(choose another output folder, or keep its files with --keep)")


def write_bundle(export_dir: str | Path, out_dir: str | Path, *, name: Optional[str] = None, fonts: str = "system", csp: bool = False,
                 guide: Optional[bool] = None, clean: bool = True, validate: bool = True, structure_scripts: bool = True,
                 embed_structures: bool = False, log=None) -> Path:
    """Write the bundle into ``out_dir`` and return ``out_dir/index.html``.

    export_dir  the folder written by ``export.write_export`` (campaign, inspect, homolog, profile, manifest, structures/)
    name        bundle name recorded in the manifest (default: the name of ``out_dir``)
    fonts       'system' (default: no network reference at all, so opening the bundle contacts no server) or 'google' (the optional Google Fonts link, system fallback)
    csp         add a Content-Security-Policy meta tag to every page
    guide       include guide.html and its assets: None (default) when viewer/pages/guide.html exists, else skipped silently
    validate    run the bundle validator at the end (file hashes, links, schemas); errors raise BundleError
    structure_scripts  also write structures/<id>.js, the lazy-loadable form of each structure (default yes)
    embed_structures   keep the binder structures (pdb_b) inside data/inspect.js. Default no: with structure scripts the script is written
                       without them (7 MB instead of 1.5 MB for 97 designs) and the inspector loads structures/<id>.js when a design opens;
                       data/inspect.json, the canonical copy, always keeps them
    """
    log = log or (lambda msg: None)
    t0 = time.perf_counter()
    out = Path(out_dir)
    if clean:
        check_clean_target(out)
    blobs = load_export(export_dir)
    manifest = blobs["manifest"]
    profile = blobs["profile"]
    manifest["bundle"]["name"] = name or out.name
    data_files = {k: blobs[k] for k in ("campaign", "inspect", "homolog")}

    def render(with_guide: Optional[bool]) -> Path:
        return viewer_build.render_site(None, out, profile=profile, data_files=data_files, manifest=manifest, guide=with_guide,
                                        fonts=fonts, csp=csp, created=manifest.get("created"), name=out.name, clean=clean)

    try:
        index = render(guide)
    except viewer_build.BuildError as ex:
        if guide is False or "guide" not in str(ex):
            raise BundleError(str(ex)) from ex
        log(f"  guide skipped: {ex}")
        index = render(False)
    log(f"  pages and data written ({time.perf_counter() - t0:.1f} s)")

    (out / "structures").mkdir(exist_ok=True)
    for src in blobs["structures"]:
        shutil.copyfile(src, out / "structures" / src.name)
        if structure_scripts:
            (out / "structures" / (src.stem + ".js")).write_text(structure_script(design_id_of(src.name), src.read_text(encoding="utf-8")), encoding="utf-8")
    if structure_scripts and not embed_structures and blobs["structures"]:
        n_slim = _slim_inspect_script(out)
        log(f"  data/inspect.js without the {n_slim} embedded structures (loaded on demand from structures/)")
    (out / "README.txt").write_text(README, encoding="utf-8")
    (out / "reviews").mkdir(exist_ok=True)

    pages = {p.stem: p.name for p in sorted(out.glob("*.html"))}
    manifest["bundle"]["pages"] = {k: pages[k] for k in ("index", "board", "inspector", "map", "query", "guide") if k in pages}
    manifest["bundle"]["data"] = {n: f"data/{n}.json" for n in BLOB_NAMES}
    manifest["bundle"]["structures"] = {"format": "pdb", "folder": "structures", "count": len(blobs["structures"]), "scripts": bool(structure_scripts),
                                        "naming": "<design id>.pdb, and <design id>.js when scripts is true; '/' of run-prefixed ids is written as '__'"}
    # every file of the bundle with its hash, except the two manifest files, which cannot hash themselves
    manifest["bundle"]["files"] = _list_files(out, {"data/manifest.json", "data/manifest.js"})
    text = dumps(manifest)
    (out / "data" / "manifest.json").write_text(text + "\n", encoding="utf-8")
    (out / "data" / "manifest.js").write_text(viewer_build.js_wrapper("manifest", text), encoding="utf-8")
    log(f"  manifest with {len(manifest['bundle']['files'])} file hashes")

    if validate:
        try:
            vs = _load_validator()
        except ImportError as ex:          # jsonschema missing: the message says what to install
            raise BundleError(f"the bundle cannot be validated: {ex} (or build with --no-validate)") from ex
        if vs is not None:
            rep = vs.check_bundle(out)
            if rep.errors:
                raise BundleError(f"{len(rep.errors)} errors in the bundle:\n  " + "\n  ".join(rep.errors[:20]))
            for w in rep.warnings:
                log(f"  validator warning: {w}")
    log(f"bundle written to {out} ({sum(f['bytes'] for f in manifest['bundle']['files']) / 1e6:.1f} MB in {len(manifest['bundle']['files'])} files, {time.perf_counter() - t0:.1f} s)")
    return index


def bundle_summary(bundle_dir: str | Path) -> dict[str, Any]:
    """Facts about a written bundle, for ``binder-review inspect``."""
    root = Path(bundle_dir)
    man = json.loads((root / "data" / "manifest.json").read_text(encoding="utf-8"))
    return {"id": man["bundle"]["id"], "name": man["bundle"]["name"], "created": man["created"], "tool": man["tool"],
            "counts": man.get("counts", {}), "pages": man["bundle"].get("pages", {}), "n_files": len(man["bundle"].get("files", [])),
            "bytes": sum(f["bytes"] for f in man["bundle"].get("files", [])), "validation": man.get("validation", {}),
            "highlights": man.get("highlights", []), "warnings": man.get("warnings", [])}


def run_build_bundle(args: argparse.Namespace) -> int:
    """CLI: ``binder-review build-bundle PROJECT_DIR --out BUNDLE_DIR``."""
    project = Path(args.project)
    export = project / "export" if (project / "export").is_dir() else project
    out = Path(args.out)
    try:
        index = write_bundle(export, out, name=getattr(args, "name", None), fonts=args.fonts, csp=args.csp,
                             guide=False if args.no_guide else None, clean=not args.keep, validate=not args.no_validate,
                             structure_scripts=not getattr(args, "no_structure_scripts", False),
                             embed_structures=getattr(args, "embed_structures", False),
                             log=lambda m: print(m, file=sys.stderr, flush=True))
    except BundleError as ex:
        print(f"binder-review build-bundle: {ex}", file=sys.stderr)
        return 1
    print(index)
    return 0
