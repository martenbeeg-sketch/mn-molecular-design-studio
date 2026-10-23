"""On-demand Results rendering for remote MD jobs."""
from __future__ import annotations

import hashlib
import json
from pathlib import PurePosixPath
from typing import Any, Callable

from mn_molecular_design_studio.md_results_ui import (
    _render_md_analysis_summary,
    _safe_relative_parts,
    render_md_results,
)


def _decode_object(text: str | None, label: str) -> dict[str, Any] | None:
    if text is None:
        return None
    try:
        value = json.loads(text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not valid JSON.") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} JSON is not an object.")
    return value


def render_remote_md_results(
    st,
    jid: str,
    row: dict[str, Any],
    token: str,
    fetch_text: Callable[..., str | None],
) -> None:
    """Render MD summaries fetched from the job's owning workstation.

    ``fetch_text`` is the authenticated artifact callback, normally
    ``cluster.ui._fetch_remote_artifact_text``. Only safe paths recorded in the
    manifest or the fixed MDTraj summary filenames are requested.
    """
    job_id = str(jid or "").strip()
    if not job_id:
        st.warning("Remote MD result is missing its job ID.")
        return
    prefix = f"remote-md-{job_id}"
    refresh_key = f"{prefix}-refresh"
    if st.button("Refresh MD view", key=refresh_key):
        session_state = getattr(st, "session_state", None)
        if session_state is not None:
            for key in list(session_state):
                if str(key).startswith(f"{prefix}-fetch-"):
                    session_state.pop(key, None)

    def load_text(relative_path: str, label: str) -> str | None:
        parts = _safe_relative_parts(relative_path)
        relative = PurePosixPath(*parts).as_posix()
        digest = hashlib.sha256(relative.encode("utf-8")).hexdigest()[:16]
        return fetch_text(job_id, relative, token, f"{prefix}-fetch-{digest}", label)

    def load_json(relative_path: str, label: str) -> dict[str, Any] | None:
        return _decode_object(load_text(relative_path, label), label)

    native_analysis = (
        str(row.get("tool") or "").casefold() == "md_analysis"
        or str(row.get("workflow") or "").casefold() == "ligand.md.analysis"
    )
    if native_analysis:
        try:
            summary = load_json("artifacts/md_analysis/summary.json", "MD analysis summary")
        except ValueError as exc:
            st.warning(str(exc))
            return
        if summary is None:
            st.info("MD analysis results are still loading. Use Refresh MD view to check again.")
            return
        if not summary.get("source_study_id"):
            _render_md_analysis_summary(st, summary, text_loader=load_text, key_prefix=prefix)
            return
        try:
            manifest = load_json("artifacts/md/study.json", "MD study manifest")
        except ValueError:
            manifest = None
        if manifest is None:
            _render_md_analysis_summary(st, summary, text_loader=load_text, key_prefix=prefix)
            return
        if manifest.get("schema_version") != 1 or not isinstance(manifest.get("studies"), list):
            st.warning("The remote MD study manifest has an unsupported schema.")
            return
        render_md_results(
            st,
            result={},
            manifest=manifest,
            json_loader=load_json,
            text_loader=load_text,
            key_prefix=prefix,
            remote=True,
        )
        return

    try:
        manifest = load_json("artifacts/md/study.json", "MD study manifest")
    except ValueError as exc:
        st.warning(str(exc))
        return
    if manifest is None:
        # This may be a standalone trajectory-analysis job with no linked study.
        try:
            summary = load_json("artifacts/md_analysis/summary.json", "MD analysis summary")
        except ValueError as exc:
            st.warning(str(exc))
            return
        if summary is None:
            st.info("MD results are still loading. Use Refresh MD view to check again.")
            return
        _render_md_analysis_summary(st, summary, text_loader=load_text, key_prefix=prefix)
        return
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("studies"), list):
        st.warning("The remote MD study manifest has an unsupported schema.")
        return
    render_md_results(
        st,
        result={},
        manifest=manifest,
        json_loader=load_json,
        text_loader=load_text,
        key_prefix=prefix,
        remote=True,
    )
