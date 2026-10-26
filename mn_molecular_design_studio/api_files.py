"""Private file handles shared by REST clients and workflow submission.

Files remain on disk. Only small descriptors are stored beside them, and no
upload is implicitly imported into the molecular database.
"""
from __future__ import annotations

import hashlib
import json
import re
import uuid
from pathlib import Path

from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.external_libraries import atomic_json

PREFIX = "studio-file:///"


def storage_name(name):
    suffixes = Path(str(name)).suffixes[-2:]
    return "content" + "".join(s for s in suffixes if re.fullmatch(r"\.[A-Za-z0-9_-]{1,16}", s))


def directory(user, file_id):
    if not re.fullmatch(r"[a-f0-9]{32}", str(file_id)):
        raise FileNotFoundError("File not found.")
    root = user_input_root(str(user["id"])) / "api_files"
    folder = root / file_id
    if folder.is_symlink() or not folder.resolve().is_relative_to(root.resolve()):
        raise FileNotFoundError("File not found.")
    return folder


def describe(user, file_id):
    folder = directory(user, file_id)
    try:
        record = json.loads((folder / "file.json").read_text())
        stored = record.get("storage_name", "content")
        if not re.fullmatch(r"content(?:\.[A-Za-z0-9_-]{1,16}){0,2}", stored):
            raise FileNotFoundError("File not found.")
        path = folder / stored
        if record["id"] != file_id or not path.is_file() or path.is_symlink():
            raise FileNotFoundError("File not found.")
        return record
    except (OSError, ValueError, KeyError) as exc:
        raise FileNotFoundError("File not found.") from exc


def resolve(user, file_id):
    record = describe(user, file_id)
    return directory(user, file_id) / record.get("storage_name", "content")


def finish(user, file_id, *, name, size, sha256, media_type):
    # Names are response metadata, never filesystem paths or HTTP headers.
    clean_name = Path(str(name).replace("\\", "/")).name
    clean_name = "".join(c for c in clean_name if c.isprintable())[:200] or "file"
    record = {"id": file_id, "name": clean_name, "size_bytes": size,
              "sha256": sha256, "media_type": media_type or "application/octet-stream",
              "storage_name": storage_name(name),
              "reference": PREFIX + file_id, "download_url": f"/api/v1/files/{file_id}"}
    atomic_json(directory(user, file_id) / "file.json", record)
    return record


def store_bytes(user, content, *, name, media_type):
    file_id = uuid.uuid4().hex
    folder = directory(user, file_id)
    folder.mkdir(parents=True, mode=0o700)
    path = folder / storage_name(name)
    with path.open("xb") as stream:
        path.chmod(0o600)
        stream.write(content)
    return finish(user, file_id, name=name, size=len(content),
                  sha256=hashlib.sha256(content).hexdigest(), media_type=media_type)


def resolve_references(value, user):
    if isinstance(value, dict):
        return {k: resolve_references(v, user) for k, v in value.items()}
    if isinstance(value, list):
        return [resolve_references(v, user) for v in value]
    if isinstance(value, str) and value.startswith(PREFIX):
        return str(resolve(user, value.removeprefix(PREFIX)))
    return value
