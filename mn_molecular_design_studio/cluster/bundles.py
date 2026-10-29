"""Input staging and machine-qualified artifact references."""
from __future__ import annotations

import copy
import hashlib
import json
import os
import stat
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

from . import config

MAX_BYTES = 512 * 1024 * 1024


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def blob_path(digest, *, root=None):
    if len(str(digest)) != 64 or any(c not in "0123456789abcdef" for c in str(digest)):
        raise ValueError("Invalid content digest.")
    root = Path(root) if root else config.home()
    path = root / "transfer-cache" / digest
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    return path


def resolve_artifact(uri, *, node_id=None, root=None):
    from mn_molecular_design_studio.runtime import runs_root
    parsed = urlsplit(str(uri))
    if parsed.scheme != "cluster" or parsed.netloc != (node_id or config.identity()):
        raise ValueError("This artifact belongs to another workstation. Open it through the cluster catalog.")
    relative = PurePosixPath(parsed.path.lstrip("/"))
    if len(relative.parts) < 2 or ".." in relative.parts or parsed.query or parsed.fragment:
        raise ValueError("Invalid machine-qualified artifact path.")
    base = Path(root) if root else runs_root()
    path = (base / str(relative)).resolve()
    if not path.is_relative_to(base.resolve()):
        raise ValueError("Artifact path escapes this workstation's Results folder.")
    return path


def referenced_nodes(value):
    nodes = set()
    if isinstance(value, dict):
        for child in value.values():
            nodes.update(referenced_nodes(child))
    elif isinstance(value, list):
        for child in value:
            nodes.update(referenced_nodes(child))
    elif isinstance(value, str) and value.startswith("cluster://"):
        parsed = urlsplit(value)
        uuid_value = parsed.netloc
        import uuid
        nodes.add(str(uuid.UUID(uuid_value)))
    return nodes


def pack(payload, *, base_dir, node_id, authorize_path=lambda p: None, preserve_path_keys=()):
    """Copy only submitted inputs; existing runs and reference files stay local."""
    from mn_molecular_design_studio.runtime import runs_root, reference_root
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import is_operational_path_key, resolve_stored_path
    run_root, references = runs_root().resolve(), reference_root().resolve()
    files, required_references = {}, set()
    size = 0
    preserve_path_keys = {str(key).casefold() for key in preserve_path_keys}

    def transform(value, key=""):
        nonlocal size
        if isinstance(value, dict):
            return {k: transform(v, str(k)) for k, v in value.items()}
        if isinstance(value, list):
            return [transform(v, key) for v in value]
        if not isinstance(value, str) or not value or not is_operational_path_key(key):
            return value
        if key.casefold() in preserve_path_keys:
            return value
        if value.startswith("cluster://") or value.startswith("cluster-input:///"):
            return value
        path = resolve_stored_path(value, run_dir=Path(base_dir), must_exist=True)
        if path is None:
            # Preserve a reference URI even if this origin has no local weights.
            if value.startswith("reference:///"):
                required_references.add(value.removeprefix("reference:///"))
            return value
        path = path.resolve()
        authorize_path(path)
        staged_reference_input = path.is_file() and path.suffix.lower() in {".a3m", ".pdb", ".cif", ".mmcif", ".fasta", ".fa"} and any(part in key.casefold() for part in ("msa", "template", "target", "structure"))
        if path.is_relative_to(references) and not staged_reference_input:
            relative = path.relative_to(references).as_posix()
            required_references.add(relative)
            return "reference:///" + relative
        if path.is_relative_to(run_root):
            relative = path.relative_to(run_root)
            if len(relative.parts) < 2:
                raise ValueError("Select a specific source run.")
            return "cluster://" + node_id + "/" + relative.as_posix()
        if path.is_dir():
            raise ValueError(f"Remote input directories must be a managed source run or reference folder: {key}")
        if not path.is_file():
            return value
        if path not in files:
            size += path.stat().st_size
            if size > MAX_BYTES:
                raise ValueError("Submitted input files exceed the 512 MB staging limit. Use a workstation that already owns these inputs.")
            files[path] = f"inputs/{len(files):05d}/{path.name}"
        return "cluster-input:///" + files[path]

    prepared = transform(copy.deepcopy(payload))
    if not files:
        return prepared, None, sorted(required_references)
    with tempfile.NamedTemporaryFile(dir=config.home(), suffix=".zip", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path, relative in files.items():
                archive.write(path, relative)
        digest = digest_file(temporary)
        destination = blob_path(digest)
        if not destination.exists():
            temporary.replace(destination)
            destination.chmod(0o600)
        return prepared, digest, sorted(required_references)
    finally:
        temporary.unlink(missing_ok=True)


def unpack(path, destination, *, expected_digest=None):
    path, destination = Path(path), Path(destination)
    if expected_digest and digest_file(path) != expected_digest:
        raise ValueError("The input transfer failed its SHA-256 integrity check.")
    destination.mkdir(parents=True, exist_ok=True, mode=0o700)
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        if len(entries) > 10000 or sum(i.file_size for i in entries) > MAX_BYTES:
            raise ValueError("Input archive exceeds the staging limit.")
        for entry in entries:
            relative = PurePosixPath(entry.filename)
            target = (destination / str(relative)).resolve()
            if relative.is_absolute() or ".." in relative.parts or not target.is_relative_to(destination.resolve()) or stat.S_ISLNK(entry.external_attr >> 16):
                raise ValueError("Unsafe input archive member.")
        for entry in entries:
            archive.extract(entry, destination)
    manifest = destination / 'input-manifest.json'
    if manifest.is_file():
        records = json.loads(manifest.read_text()).get('sha256', {})
        for relative, checksum in records.items():
            file = (destination / relative).resolve()
            if not file.is_relative_to(destination.resolve()) or not file.is_file() or digest_file(file) != checksum:
                raise ValueError('A candidate input failed its SHA-256 integrity check.')

    # Candidate collections contain file references inside the transferred file.
    # Verify native bytes first, then derive a worker-local collection.
    if manifest.is_file():
        for relative in records:
            path = destination / relative
            if path.suffix == '.jsonl':
                rows = [restore(json.loads(line), staging=destination, node_id=config.identity())
                        for line in path.read_text().splitlines() if line.strip()]
                path.write_text(''.join(json.dumps(row) + '\n' for row in rows))


def restore(payload, *, staging, node_id):
    if isinstance(payload, dict):
        # Historical evidence may contain origin-only paths and cluster URIs.
        # It describes prior work and is never an input to restore on this node.
        return {k: (copy.deepcopy(v) if k in {'metrics', 'native_measurements', 'source_provenance'}
                    else restore(v, staging=staging, node_id=node_id))
                for k, v in payload.items()}
    if isinstance(payload, list):
        return [restore(v, staging=staging, node_id=node_id) for v in payload]
    if isinstance(payload, str) and payload.startswith("cluster://"):
        return str(resolve_artifact(payload, node_id=node_id))
    if isinstance(payload, str) and payload.startswith("cluster-input:///"):
        relative = PurePosixPath(payload.removeprefix("cluster-input:///"))
        path = (Path(staging) / str(relative)).resolve()
        if relative.is_absolute() or ".." in relative.parts or not path.is_relative_to(Path(staging).resolve()) or not path.exists():
            raise ValueError("A staged input file or directory is unavailable.")
        return str(path)
    return payload
