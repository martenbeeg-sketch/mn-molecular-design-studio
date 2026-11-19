"""Assemble the Binderscope viewer: page templates + profile + data blobs -> a static site that opens at file://.

The templates live in ``viewer/pages`` (index, board, inspector, map, query), shared code in ``viewer/js`` and ``viewer/css``,
the vendored 3Dmol.js in ``viewer/vendor``. A template is plain HTML with four placeholders, replaced here:

    <!--@head-->    charset, viewport, title, theme script, fonts link, shared stylesheet, icon (and the CSP meta when asked)
    <!--@brand-->   the wordmark: inline mark, product name, target label
    <!--@badge-->   the data badge (profile.badge) with its tooltip
    <!--@nav-->     the run filter slot and the page links

Shared code is linked, never inlined and never loaded from a CDN: ``assets/css/shared.css``, ``assets/js/*.js``,
``assets/3Dmol-min.js``. Data is loaded as same-folder scripts (``data/campaign.js`` sets ``window.RAW`` and so on) because
``fetch()`` of a local file is blocked at file://. The only network reference is the optional Google Fonts link (a system
font stack is the fallback; ``fonts="system"`` leaves it out).

``render_site`` is the function the exporter calls; ``render_page`` renders one page to a string (tests, previews).
See docs/viewer.md.
"""
from __future__ import annotations

import hashlib
import html
import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

from . import __version__

PAGES = ("index", "board", "inspector", "map", "query")
NAV_ORDER = ("index", "board", "inspector", "map", "query")
NAV_LABEL = {"index": "Overview", "board": "Board", "inspector": "Inspector", "map": "Map", "query": "Query"}
DEFAULT_TITLES = {"index": "Overview", "board": "Binder Family Board", "inspector": "Binder Design Inspector", "map": "Binder Campaign Map", "query": "Binder Structure Query"}
BLOBS = ("campaign", "inspect", "homolog", "profile", "manifest")
GLOBAL_OF = {"campaign": "RAW", "inspect": "INS", "homolog": "HOM", "profile": "PROFILE", "manifest": "MANIFEST"}
# data blobs each page loads (the same list is written into its template; checked after rendering)
NEEDS = {"index": ("campaign", "manifest"), "board": ("campaign", "manifest"), "map": ("campaign", "homolog", "manifest"),
         "query": ("campaign", "homolog", "manifest"), "inspector": ("inspect", "homolog", "manifest")}
PLACEHOLDERS = ("head", "brand", "badge", "nav")

FONTS_URL = "https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&display=swap"
THEME_SCRIPT = ('<script>if(location.hash=="#light")document.documentElement.setAttribute("data-theme","light");'
                'if(location.hash=="#dark")document.documentElement.setAttribute("data-theme","dark");</script>')
CSP_STRICT = ("default-src 'none'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; "
              "font-src https://fonts.gstatic.com data:; img-src 'self' data: blob:; connect-src 'none'; worker-src blob:; base-uri 'none'; form-action 'none'")
CSP_NO_FONTS = CSP_STRICT.replace(" https://fonts.googleapis.com", "").replace("font-src https://fonts.gstatic.com data:", "font-src data:")


class BuildError(Exception):
    """The site cannot be assembled (missing data, a link that points nowhere, a profile that disagrees with a blob)."""


def viewer_dir() -> Path:
    """The viewer sources: $BINDERSCOPE_VIEWER, else ``viewer/`` of the source checkout, else ``viewer/`` inside the package."""
    env = os.environ.get("BINDERSCOPE_VIEWER")
    if env:  # an explicit override is never replaced silently by another folder
        c = Path(env)
        if (c / "pages").is_dir() and (c / "js" / "bs.js").is_file():
            return c
        raise BuildError(f"BINDERSCOPE_VIEWER={env} has no pages/ and js/bs.js")
    here = Path(__file__).resolve()
    cands = [here.parents[2] / "viewer", here.parent / "viewer"]
    for c in cands:
        if (c / "pages").is_dir() and (c / "js" / "bs.js").is_file():
            return c
    raise BuildError("viewer sources not found (expected viewer/pages and viewer/js/bs.js next to src/); set BINDERSCOPE_VIEWER")


# --------------------------------------------------------------------------------------------------------------------
# text templates (profile_spec.md section 3): {token}, {token|n}, {@profile.path}
# --------------------------------------------------------------------------------------------------------------------
_TOKEN = re.compile(r"\{([^{}|]+)(?:\|([a-z]+))?\}")


def render_template(tpl: str, profile: Mapping[str, Any], tokens: Mapping[str, Any] | None = None) -> str:
    """Python twin of ``BS.render``: an unknown token stays visible, an unresolved {@path} renders empty."""
    tokens = tokens or {}

    def sub(m: re.Match) -> str:
        name, fmt = m.group(1), m.group(2)
        if name.startswith("@"):
            cur: Any = profile
            for k in name[1:].split("."):
                cur = cur[k] if isinstance(cur, Mapping) and k in cur else ""
            v = cur
        else:
            v = tokens[name] if name in tokens else m.group(0)
        if fmt == "n" and isinstance(v, (int, float)) and not isinstance(v, bool):
            v = f"{v:,}"
        return str(v)

    return _TOKEN.sub(sub, tpl)


# --------------------------------------------------------------------------------------------------------------------
# pieces of a page
# --------------------------------------------------------------------------------------------------------------------
def _esc(s: Any) -> str:
    return html.escape(str(s), quote=True)


def mark_svg(vdir: Path | None = None) -> str:
    """The inline logo mark (viewer/assets/binderscope-mark.svg) as a decorative element; the ring follows currentColor."""
    svg = ((vdir or viewer_dir()) / "assets" / "binderscope-mark.svg").read_text(encoding="utf-8").strip()
    svg = re.sub(r'\s+role="img"\s+aria-label="[^"]*"', "", svg)
    return svg.replace("<svg ", '<svg class="mark" aria-hidden="true" focusable="false" ', 1)


def favicon_svg(vdir: Path | None = None) -> str:
    """The mark with explicit colours (currentColor means black in a tab strip, invisible on a dark one)."""
    svg = ((vdir or viewer_dir()) / "assets" / "binderscope-mark.svg").read_text(encoding="utf-8").strip()
    return svg.replace("currentColor", "#4a3aa7")


def _profile(profile: Mapping[str, Any] | None) -> Mapping[str, Any]:
    return profile or {}


def render_brand(profile: Mapping[str, Any] | None, vdir: Path | None = None) -> str:
    p = _profile(profile)
    name = (p.get("brand") or {}).get("name") or "Binderscope"
    short = (p.get("target") or {}).get("short") or ""
    proj = f'<span class="bproj" title="{_esc((p.get("target") or {}).get("name") or short)}">{_esc(short)}</span>' if short else ""
    return (f'<a class="brand" href="index.html" title="Project overview">{mark_svg(vdir)}<span class="bname">{_esc(name)}</span>{proj}</a>')


def render_badge(profile: Mapping[str, Any] | None, page: str, tokens: Mapping[str, Any] | None = None) -> str:
    p = _profile(profile)
    b = p.get("badge")
    if not b or not b.get("text"):
        return '<span class="tag-mock" id="tagMock" hidden></span>'
    cls = {"mock": "tag-mock", "real": "tag-real", "run": "tag-run"}.get(b.get("class"), "tag-mock")
    tpl = (p.get("text") or {}).get(f"badge.title.{page}")
    title = render_template(tpl, p, tokens or {}) if tpl else ""
    t = f' title="{_esc(title)}"' if title else ""
    return f'<span class="{cls}" id="tagMock"{t}>{_esc(b["text"])}</span>'


def render_nav(page: str, available: Iterable[str], guide: bool = False) -> str:
    have = set(available)
    links = []
    for name in NAV_ORDER:
        if name not in have:
            continue
        cur = ' aria-current="page"' if name == page else ""
        links.append(f'<a href="{name}.html" data-nav="{name}"{cur}>{NAV_LABEL[name]}</a>')
    if guide:
        links.append('<a class="extra" href="guide.html" data-nav="guide">Guide</a>')
    return '<div class="runslot" id="bsRunSlot" hidden></div>\n    <nav class="views" aria-label="Views">' + "".join(links) + "</nav>"


def render_head(page: str, profile: Mapping[str, Any] | None, *, fonts: str = "google", csp: bool | str = False) -> str:
    p = _profile(profile)
    title = ((p.get("brand") or {}).get("titles") or {}).get(page)
    if not title:
        title = DEFAULT_TITLES[page]
        if page == "index":  # the other pages read 'Family board . <target> . Binderscope': the tab of the first page says the same
            title = " · ".join(x for x in (title, (p.get("target") or {}).get("short"), (p.get("brand") or {}).get("name")) if x)
    parts = ['<meta charset="utf-8">', '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">']
    if csp:
        policy = csp if isinstance(csp, str) else (CSP_STRICT if fonts == "google" else CSP_NO_FONTS)
        parts.append(f'<meta http-equiv="Content-Security-Policy" content="{html.escape(policy, quote=False).replace(chr(34), "&quot;")}">')
    parts += [f"<title>{_esc(title)}</title>", THEME_SCRIPT, '<link rel="icon" type="image/svg+xml" href="assets/favicon.svg">']
    if fonts == "google":
        parts.append(f'<link rel="stylesheet" href="{FONTS_URL}">')
    parts.append('<link rel="stylesheet" href="assets/css/shared.css">')
    return "\n".join(parts)


def render_page(name: str, profile: Mapping[str, Any] | None = None, *,
                run: Mapping[str, Any] | None = None, available: Iterable[str] | None = None, guide: bool = False,
                fonts: str = "google", csp: bool | str = False, viewer: Path | None = None) -> str:
    """One page as an HTML string. ``profile`` drives the wordmark, titles, badge and tooltip; ``run`` (campaign.run) supplies the
    tokens of the badge tooltip ({run_license}, {run_credit}); ``available`` lists the pages that exist (the nav shows only those)."""
    if name not in PAGES:
        raise BuildError(f"unknown page {name!r}; pages are {', '.join(PAGES)}")
    if fonts not in ("google", "system"):
        raise BuildError("fonts must be 'google' or 'system'")
    vdir = viewer or viewer_dir()
    src = (vdir / "pages" / f"{name}.html").read_text(encoding="utf-8")
    # without a run block the tokens stay visible in the static tooltip; bs.js fills them in from the data when the page opens
    tokens = {} if run is None else {"run_license": run.get("license", ""), "run_credit": run.get("credit", ""),
                                     "run_source": run.get("source", ""), "run_name": run.get("name", "")}
    avail = list(available) if available is not None else list(PAGES)
    rep = {
        "head": render_head(name, profile, fonts=fonts, csp=csp),
        "brand": render_brand(profile, vdir),
        "badge": render_badge(profile, name, tokens),
        "nav": render_nav(name, avail, guide),
    }
    for key, val in rep.items():
        tag = f"<!--@{key}-->"
        if src.count(tag) != 1:
            raise BuildError(f"template {name}.html must contain {tag} exactly once (found {src.count(tag)})")
        src = src.replace(tag, val)
    if "<body>" not in src:
        raise BuildError(f"template {name}.html has no <body> tag")
    src = src.replace("<body>", f'<body data-page="{name}">', 1)
    left = re.findall(r"<!--@\w+-->", src)
    if left:
        raise BuildError(f"template {name}.html has unknown placeholders {left}")
    return src


# --------------------------------------------------------------------------------------------------------------------
# data blobs
# --------------------------------------------------------------------------------------------------------------------
def js_wrapper(blob: str, json_text: str) -> str:
    """data/<blob>.js: the same object as data/<blob>.json for file:// (fetch cannot read local files)."""
    return f"window.BINDERSCOPE = window.BINDERSCOPE || {{}};\nwindow.BINDERSCOPE.{blob} = window.{GLOBAL_OF[blob]} = {json_text.strip()};\n"


def _load_profile(profile: Any) -> dict | None:
    if profile is None:
        return None
    if isinstance(profile, Mapping):
        return dict(profile)
    return json.loads(Path(profile).read_text(encoding="utf-8"))


def _dump(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def collect_blobs(project_or_data_dir: Any = None, data_files: Mapping[str, Any] | None = None) -> dict[str, str]:
    """JSON text of every blob found: ``<dir>/<name>.json`` or ``<dir>/data/<name>.json``, overridden by ``data_files`` (a path or an object per blob)."""
    out: dict[str, str] = {}
    if project_or_data_dir is not None:
        d = Path(project_or_data_dir)
        for b in BLOBS:
            for cand in (d / f"{b}.json", d / "data" / f"{b}.json"):
                if cand.is_file():
                    out[b] = cand.read_text(encoding="utf-8")
                    break
    for b, v in (data_files or {}).items():
        if b not in BLOBS:
            raise BuildError(f"unknown data blob {b!r}; blobs are {', '.join(BLOBS)}")
        out[b] = _dump(v) if isinstance(v, (Mapping, list)) else Path(v).read_text(encoding="utf-8")
    return out


def synthesize_manifest(blobs: Mapping[str, str], profile: Mapping[str, Any] | None, name: str, pages: Iterable[str],
                        created: str | None = None, title: str | None = None) -> dict:
    """A minimal manifest for sites built without one (previews, the regression sites). A real build supplies its own."""
    h = hashlib.sha256()
    for b in sorted(blobs):
        h.update(hashlib.sha256(blobs[b].encode("utf-8")).digest())
    p = profile or {}
    tgt = p.get("target") or {}
    camp = json.loads(blobs["campaign"]) if "campaign" in blobs else {}
    runs = camp.get("runs")
    if not runs:
        r = camp.get("run") or {}
        runs = [{"id": "run", "label": r.get("name") or tgt.get("short") or "run", "reader": "legacy", "source": r.get("source", ""), "experimental": False}]
    n_acc = sum(1 for d in camp.get("designs", []) if d.get("status", "accepted") == "accepted")
    return {
        "schema": "binderscope.manifest/1",
        "tool": {"name": "binderscope-viewer", "package": "binder-review", "version": __version__, "command": "viewer_build.render_site"},
        "created": created or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "bundle": {"id": h.hexdigest()[:16], "name": name, "title": title or tgt.get("name") or name,
                   "pages": {k: f"{k}.html" for k in pages}, "data": {b: f"data/{b}.json" for b in sorted(set(blobs) | {"manifest"})}, "files": []},
        "profile": {"id": p.get("id", "default")},
        "target": {"name": tgt.get("name", ""), "uniprot": tgt.get("uniprot", ""), "reference_pdb": (tgt.get("reference") or {}).get("pdb", ""), "organism": tgt.get("organism", "")},
        "sources": [],
        "runs": runs,
        "caps": dict(p.get("caps") or {}),
        "counts": {"designs": len(camp.get("designs", [])), "accepted": n_acc},
        "notes": ["Manifest synthesised by viewer_build.render_site: the site was built from data blobs without an export manifest."],
    }


_REF = re.compile(r"""(?:src|href)=["']([^"'#]+)(?:#[^"']*)?["']""")


def check_links(out_dir: Path, pages: Iterable[str]) -> list[str]:
    """Relative src and href targets of the written pages that do not exist in the site folder (external URLs are ignored)."""
    missing = []
    for p in pages:
        text = (out_dir / f"{p}.html").read_text(encoding="utf-8")
        text = re.sub(r"<script(?![^>]*\bsrc=)[^>]*>.*?</script>", "<script></script>", text, flags=re.S)  # code strings are not links
        text = re.sub(r"<style[^>]*>.*?</style>", "", text, flags=re.S)
        for ref in _REF.findall(text):
            if re.match(r"^(https?:|data:|mailto:|javascript:|#)", ref) or ref == "":
                continue
            if not (out_dir / ref).exists():
                missing.append(f"{p}.html -> {ref}")
    return sorted(set(missing))


def render_site(project_or_data_dir: Any, out_dir: Any, profile: Any = None, data_files: Mapping[str, Any] | None = None, *,
                manifest: Mapping[str, Any] | None = None, pages: Iterable[str] | None = None, guide: bool | None = None,
                fonts: str = "google", csp: bool | str = False, created: str | None = None, name: str | None = None,
                clean: bool = False, inject_profile: bool = True) -> Path:
    """Write the viewer site into ``out_dir`` and return ``out_dir/index.html``.

    project_or_data_dir  folder holding campaign.json, inspect.json, homolog.json, profile.json, manifest.json (here or in data/),
                         plus optional structures/; None when ``data_files`` gives everything.
    profile              the project profile (dict or path); default: data/profile.json, else the profile inside the campaign blob.
    data_files           {blob: path or object} overriding or adding blobs.
    manifest             dict (or path); default data/manifest.json, else a minimal one is written (see ``synthesize_manifest``).
    pages                pages to build (default: every page whose data is present; index always).
    guide                copy guide.html and its assets and add the nav link (default: yes when viewer/pages/guide.html exists).
    fonts                'google' (link with system-font fallback) or 'system' (no network reference at all).
    csp                  True or a policy string: add a Content-Security-Policy meta tag to every page.
    inject_profile       a blob without ``profile`` gets the profile embedded (the three copies must be identical).
    """
    vdir = viewer_dir()
    out = Path(out_dir)
    if clean and out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    blobs = collect_blobs(project_or_data_dir, data_files)
    if "campaign" not in blobs and "inspect" not in blobs:
        raise BuildError("no campaign.json or inspect.json: nothing to show")

    prof = _load_profile(profile)
    if prof is None and "profile" in blobs:
        prof = json.loads(blobs["profile"])
    if prof is None:
        for b in ("campaign", "inspect"):
            if b in blobs and json.loads(blobs[b]).get("profile"):
                prof = json.loads(blobs[b])["profile"]
                break
    # the profile travels inside the blobs the pages load (campaign, inspect); every copy must be the same object
    if prof is not None:
        for b in ("campaign", "inspect"):
            if b not in blobs:
                continue
            obj = json.loads(blobs[b])
            if "profile" not in obj:
                if not inject_profile:
                    raise BuildError(f"{b} blob has no profile and inject_profile is off")
                obj["profile"] = prof
                blobs[b] = _dump(obj)
            elif obj["profile"] != prof:
                raise BuildError(f"the profile embedded in the {b} blob differs from the profile given")
        blobs["profile"] = _dump(prof)

    wanted = list(pages) if pages is not None else [p for p in PAGES if p == "index" or all(b in blobs or b == "manifest" for b in NEEDS[p])]
    if "index" not in wanted:
        wanted.insert(0, "index")
    for p in wanted:
        if p not in PAGES:
            raise BuildError(f"unknown page {p!r}")
        missing = [b for b in NEEDS[p] if b not in blobs and b != "manifest"]
        if missing:
            raise BuildError(f"page {p} needs the {', '.join(missing)} blob(s)")
    if guide is None:
        guide = (vdir / "pages" / "guide.html").is_file()

    man: dict | None = None
    if manifest is not None:
        man = dict(manifest) if isinstance(manifest, Mapping) else json.loads(Path(manifest).read_text(encoding="utf-8"))
    elif "manifest" in blobs:
        man = json.loads(blobs["manifest"])
    if man is None:
        man = synthesize_manifest(blobs, prof, name or out.name, wanted, created)
    blobs["manifest"] = _dump(man)

    run = (json.loads(blobs["campaign"]).get("run") if "campaign" in blobs else None) or (json.loads(blobs["inspect"]).get("run") if "inspect" in blobs else None)
    for p in wanted:
        text = render_page(p, prof, run=run, available=wanted, guide=bool(guide), fonts=fonts, csp=csp, viewer=vdir)
        (out / f"{p}.html").write_text(text, encoding="utf-8")

    # shared files, copied (never linked to the source folder): the bundle must stand alone
    for sub in ("assets/css", "assets/js"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    shutil.copyfile(vdir / "css" / "shared.css", out / "assets" / "css" / "shared.css")
    for js in sorted((vdir / "js").glob("*.js")):
        shutil.copyfile(js, out / "assets" / "js" / js.name)
    shutil.copyfile(vdir / "vendor" / "3Dmol-min.js", out / "assets" / "3Dmol-min.js")
    shutil.copyfile(vdir / "vendor" / "3Dmol-LICENSE.txt", out / "assets" / "3Dmol-LICENSE.txt")
    shutil.copyfile(vdir / "assets" / "binderscope-mark.svg", out / "assets" / "binderscope-mark.svg")
    (out / "assets" / "favicon.svg").write_text(favicon_svg(vdir) + "\n", encoding="utf-8")
    if guide:
        ttext = (vdir / "pages" / "guide.html").read_text(encoding="utf-8")
        if fonts != "google":  # the page source links the web fonts; a system-fonts bundle has no network reference
            ttext = re.sub(r'<link rel="stylesheet" href="https://fonts\.googleapis\.com[^>]*>\n?', "", ttext)
        (out / "guide.html").write_text(ttext, encoding="utf-8")
        tdir = vdir / "assets" / "guide"
        if tdir.is_dir():
            shutil.copytree(tdir, out / "assets" / "guide", dirs_exist_ok=True)

    (out / "data").mkdir(exist_ok=True)
    for b, text in blobs.items():
        (out / "data" / f"{b}.json").write_text(text.strip() + "\n", encoding="utf-8")
        (out / "data" / f"{b}.js").write_text(js_wrapper(b, text), encoding="utf-8")

    if project_or_data_dir is not None:
        sdir = Path(project_or_data_dir)
        for cand in (sdir / "structures", sdir / "data" / "structures"):
            if cand.is_dir():
                shutil.copytree(cand, out / "structures", dirs_exist_ok=True)
                break
    rdir = out / "reviews"
    rdir.mkdir(exist_ok=True)
    (rdir / "README.txt").write_text("Review files saved from the board (Save review file) belong here, one JSONL file per reviewer.\n"
                                     "Browsers cannot write into this folder themselves: move the downloaded file here to keep it with the bundle.\n", encoding="utf-8")

    bad = check_links(out, wanted + (["guide"] if guide else []))
    if bad:
        raise BuildError("links to files that are not in the site: " + "; ".join(bad[:10]))
    return out / "index.html"


def main(argv: list[str] | None = None) -> int:
    """python -m binder_review.viewer_build DATA_DIR OUT_DIR: write the viewer site for a folder of data blobs."""
    import argparse

    ap = argparse.ArgumentParser(prog="python -m binder_review.viewer_build", description=main.__doc__)
    ap.add_argument("data_dir", help="folder with campaign.json, inspect.json, homolog.json, profile.json, manifest.json (here or in data/)")
    ap.add_argument("out_dir")
    ap.add_argument("--profile", help="profile JSON (default: data/profile.json or the profile inside the campaign blob)")
    ap.add_argument("--fonts", choices=("google", "system"), default="google")
    ap.add_argument("--csp", action="store_true", help="add a Content-Security-Policy meta tag")
    ap.add_argument("--no-guide", dest="no_guide", action="store_true")
    ap.add_argument("--no-tutorial", dest="no_guide", action="store_true", help=argparse.SUPPRESS)   # old name of --no-guide
    ap.add_argument("--clean", action="store_true", help="remove OUT_DIR first")
    ap.add_argument("--name")
    a = ap.parse_args(argv)
    try:
        idx = render_site(a.data_dir, a.out_dir, profile=a.profile, fonts=a.fonts, csp=a.csp, guide=False if a.no_guide else None, clean=a.clean, name=a.name)
    except BuildError as e:
        print(f"error: {e}")
        return 1
    print(idx)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
