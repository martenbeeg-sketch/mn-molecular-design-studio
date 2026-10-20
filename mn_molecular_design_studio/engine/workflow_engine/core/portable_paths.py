from __future__ import annotations

from pathlib import Path, PurePosixPath
from typing import Final
from urllib.parse import urlsplit

from mn_molecular_design_studio.engine.workflow_engine.runtime import app_home, reference_root, runs_root


RUNS_SCHEME: Final = "runs"
REFERENCE_SCHEME: Final = "reference"
APP_SCHEME: Final = "app"
PORTABLE_SCHEMES: Final = frozenset({RUNS_SCHEME, REFERENCE_SCHEME, APP_SCHEME, "cluster"})
_PATH_KEY_PARTS: Final = (
    "cache", "checkpoint", "database", "dataset", "directory", "file", "folder", "input",
    "model_path", "output", "path", "repository", "root", "run_dir", "source", "structure",
    "target", "weights",
)
_PATH_KEY_SUFFIXES: Final = (
    "_a3m", "_cif", "_csv", "_dir", "_fasta", "_fa", "_file", "_folder", "_json", "_jsonl",
    "_mmcif", "_msa", "_pdb", "_path", "_root", "_structure", "_table", "_tsv", "_weights", "_yaml", "_yml",
)
_PROVENANCE_KEY_PARTS: Final = ("command", "argv", "executable", "queued_command")
_INLINE_CONTENT_KEYS: Final = frozenset({
    "bias_by_res_json", "design_yaml", "target_msa_a3m", "text", "sequence", "sequences", "smiles",
})
_INLINE_CONTENT_SUFFIXES: Final = ("_text", "_sequence", "_sequences", "_smiles")
_PATH_VALUE_EXTENSIONS: Final = frozenset({
    ".a3m", ".cif", ".csv", ".fa", ".fasta", ".json", ".jsonl", ".mmcif",
    ".msa", ".pdb", ".tsv", ".yaml", ".yml",
})
_RELATIVE_PATH_KEY_SUFFIXES: Final = (
    "_a3m", "_cif", "_csv", "_dir", "_directory", "_fa", "_fasta", "_file", "_folder",
    "_json", "_jsonl", "_mmcif", "_msa", "_pdb", "_path", "_repository", "_root",
    "_structure", "_table", "_tsv", "_weights", "_yaml", "_yml",
)


def _safe_relative(value: str) -> Path | None:
    relative = PurePosixPath(value.lstrip("/"))
    if not relative.parts or ".." in relative.parts:
        return None
    return Path(*relative.parts)


def _root_for(scheme: str) -> Path:
    return {
        RUNS_SCHEME: runs_root(),
        REFERENCE_SCHEME: reference_root(),
        APP_SCHEME: app_home(),
    }[scheme]


def _uri_candidate(raw: str) -> Path | None:
    parsed = urlsplit(raw)
    if parsed.scheme == "cluster":
        from mn_molecular_design_studio.cluster.bundles import resolve_artifact
        return resolve_artifact(raw)
    if parsed.scheme not in PORTABLE_SCHEMES:
        return None
    joined = "/".join(part for part in (parsed.netloc, parsed.path) if part)
    relative = _safe_relative(joined)
    if relative is None:
        return None
    return _root_for(parsed.scheme) / relative


def is_portable_path(value: object) -> bool:
    return urlsplit(str(value or "").strip()).scheme in PORTABLE_SCHEMES


def is_operational_path_key(key: str) -> bool:
    lowered = key.lower()
    if _is_inline_content_key(lowered):
        return False
    # Format names such as "pdb" and "mmcif" are values, even when the key
    # contains "input" or "output"; they must not be resolved as filesystem paths.
    if lowered.endswith("_format"):
        return False
    return any(part in lowered for part in _PATH_KEY_PARTS) or lowered.endswith(_PATH_KEY_SUFFIXES)


def _is_inline_content_key(key: str) -> bool:
    """Distinguish pasted engine inputs from file references with similar names."""
    lowered = key.lower()
    return lowered in _INLINE_CONTENT_KEYS or lowered.endswith(_INLINE_CONTENT_SUFFIXES)


def is_provenance_key(key: str) -> bool:
    lowered = key.lower()
    return any(part in lowered for part in _PROVENANCE_KEY_PARTS)


def _relative_path_key(key: str) -> bool:
    lowered = key.lower()
    return lowered == "path" or lowered.endswith(_RELATIVE_PATH_KEY_SUFFIXES)


def _is_named_resource_key(key: str) -> bool:
    """Return whether a file-like setting value names a bundled preset, not a path."""
    return key.lower() in {"advanced_settings_file", "model_weights"}


def store_managed_paths(value: object, *, run_dir: Path | None = None, key: str = "", provenance: bool = False) -> object:
    """Encode managed absolute paths in JSON payloads without touching the caller's data."""
    provenance = provenance or is_provenance_key(key)
    if isinstance(value, dict):
        return {
            child_key: store_managed_paths(child, run_dir=run_dir, key=str(child_key), provenance=provenance)
            for child_key, child in value.items()
        }
    if isinstance(value, list):
        return [store_managed_paths(child, run_dir=run_dir, key=key, provenance=provenance) for child in value]
    if not isinstance(value, str) or provenance or not is_operational_path_key(key) or not value.startswith("/"):
        return value
    stored = portable_path(value, run_dir=run_dir)
    if run_dir is not None and stored != value and not is_portable_path(stored) and not _relative_path_key(key):
        try:
            relative_to_runs = Path(value).expanduser().resolve().relative_to(runs_root().expanduser().resolve())
        except ValueError:
            pass
        else:
            return f"runs:///{relative_to_runs.as_posix()}"
    return stored if stored != value else value


def resolve_managed_paths(value: object, *, run_dir: Path | None = None, key: str = "", provenance: bool = False) -> object:
    """Resolve stored managed path references in a JSON payload for current use."""
    provenance = provenance or is_provenance_key(key)
    if isinstance(value, dict):
        return {
            child_key: resolve_managed_paths(child, run_dir=run_dir, key=str(child_key), provenance=provenance)
            for child_key, child in value.items()
        }
    if isinstance(value, list):
        return [resolve_managed_paths(child, run_dir=run_dir, key=key, provenance=provenance) for child in value]
    if not isinstance(value, str) or provenance or _is_inline_content_key(key):
        return value
    if is_portable_path(value):
        resolved = resolve_stored_path(value, run_dir=run_dir)
        return str(resolved) if resolved is not None else value
    if _is_named_resource_key(key) and not value.startswith((".", "/")) and "/" not in value and "\\" not in value:
        return value
    path_value = value.startswith(".") or Path(value).suffix.lower() in _PATH_VALUE_EXTENSIONS
    path_key = _relative_path_key(key)
    if is_operational_path_key(key) and not Path(value).is_absolute() and (path_value or path_key):
        resolved = resolve_stored_path(value, run_dir=run_dir)
        return str(resolved) if resolved is not None else value
    return value


def _legacy_candidate(raw: str) -> Path | None:
    normalized = raw.replace("\\", "/")
    mappings = (
        ("/workdir/runs/", runs_root()),
        ("/reference_files/", reference_root()),
        ("/workdir/targets/", app_home() / "workdir" / "targets"),
    )
    for marker, root in mappings:
        if marker not in normalized:
            continue
        relative = _safe_relative(normalized.split(marker, 1)[1])
        if relative is not None:
            return root / relative
    return None


def resolve_stored_path(
    value: str | Path | None,
    *,
    run_dir: Path | None = None,
    must_exist: bool = False,
) -> Path | None:
    """Resolve a portable URI, run-relative path, or legacy absolute path.

    A still-existing recorded absolute path wins. Legacy relocation is tried
    only when that path is gone, so this remains compatible on the source
    computer while allowing managed app data to move to a new root.
    """
    if value is None:
        return None
    raw = str(value).strip()
    if not raw:
        return None

    candidate = _uri_candidate(raw)
    if candidate is None:
        recorded = Path(raw).expanduser()
        if recorded.exists():
            return recorded
        if not recorded.is_absolute():
            base = Path(run_dir).expanduser() if run_dir is not None else runs_root()
            candidate = base / recorded
        else:
            candidate = _legacy_candidate(raw)
            if candidate is None:
                return None if must_exist else recorded

    if must_exist and not candidate.exists():
        return None
    return candidate


def portable_path(value: str | Path, *, run_dir: Path | None = None) -> str:
    """Encode paths under managed roots without changing external paths."""
    raw = str(value).strip()
    if urlsplit(raw).scheme in PORTABLE_SCHEMES:
        return raw
    path = Path(value).expanduser()
    if not path.is_absolute() and run_dir is not None:
        path = Path(run_dir).expanduser() / path
    path = path.resolve()
    if run_dir is not None:
        try:
            return path.relative_to(Path(run_dir).expanduser().resolve()).as_posix()
        except ValueError:
            pass
    for scheme, root in (
        (RUNS_SCHEME, runs_root()),
        (REFERENCE_SCHEME, reference_root()),
        (APP_SCHEME, app_home()),
    ):
        try:
            relative = path.relative_to(root.resolve()).as_posix()
        except ValueError:
            continue
        return f"{scheme}:///{relative}"
    return str(path)
