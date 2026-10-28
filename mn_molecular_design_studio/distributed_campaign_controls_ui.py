"""Safe pause, resume and cancellation controls for local campaign jobs."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable


_ACTIVE = {"queued", "pending", "starting", "preparing", "running", "cancelling", "cancel_requested", "scheduled", "waiting", "waiting_for_workstation"}
_PAUSED = {"paused", "stopped"}


def render_campaign_lifecycle_controls(
    st: Any,
    *,
    reference: str,
    run_dir: str | Path,
    status: str,
    context: dict[str, Any] | None,
    cancel: Callable[[str], Any],
) -> None:
    """Expose campaign-specific controls while respecting resumable job states."""
    if not context:
        return
    folder = Path(run_dir)
    resumable = (folder / "worker_request.json").is_file()
    normalized = str(status or "").casefold()
    if normalized in _ACTIVE:
        columns = st.columns(2 if resumable else 1)
        if resumable and columns[0].button("Pause campaign controller", key="campaign-pause-" + reference):
            try:
                from mn_molecular_design_studio.engine.workflow_engine.core.jobs import pause_job

                pause_job(folder, reason="Campaign controller paused from Results")
            except (OSError, RuntimeError, ValueError) as exc:
                st.error(f"Could not pause the campaign controller: {exc}")
            else:
                st.success("Campaign controller paused. Active attempts may finish and will be reconciled when resumed.")
                st.rerun()
        cancel_column = columns[1] if resumable else columns[0]
        if cancel_column.button("Cancel campaign", key="campaign-cancel-" + reference):
            try:
                cancel(reference)
            except (OSError, RuntimeError, ValueError) as exc:
                st.error(f"Could not cancel the campaign: {exc}")
            else:
                st.success("Campaign cancellation requested. Active attempt and validation jobs are being stopped.")
                st.rerun()
        if resumable:
            st.caption("Pausing stops the controller from authorizing attempts; cancel terminates the campaign and cannot be resumed.")
    elif normalized in _PAUSED:
        if resumable and st.button("Resume campaign", key="campaign-resume-" + reference):
            try:
                from mn_molecular_design_studio.engine.workflow_engine.core.jobs import resume_job

                resume_job(folder)
            except (FileNotFoundError, OSError, RuntimeError, ValueError) as exc:
                st.error(f"Could not resume the campaign: {exc}")
            else:
                st.success("Campaign controller resumed from its saved ledger and worker receipts.")
                st.rerun()
        elif not resumable:
            st.info("This campaign has no saved worker request and cannot be resumed in place.")
