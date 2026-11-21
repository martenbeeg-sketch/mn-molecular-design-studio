"""Shared prediction selection and confidence plots for local and peer Results."""
from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import streamlit as st

from mn_molecular_design_studio import cofolding_results as results


def select_models(candidates, *, key, input_labels=None, input_label="Input complex"):
    by_input = {}
    for candidate in candidates:
        by_input.setdefault(results.complex_group_key(candidate), []).append(candidate)
    input_ids = list(by_input)
    input_index = 0
    if len(by_input) > 1:
        labels = [results.complex_label(items[0], f"Complex {index + 1}", input_labels)
                  for index, items in enumerate(by_input.values())]
        input_index = st.selectbox(
            input_label, range(len(labels)), format_func=labels.__getitem__, key=f"{key}-input",
        )
    models = by_input[input_ids[input_index]]
    model_labels = results.model_labels(models)
    selected = st.selectbox(
        "Model", range(len(models)), format_func=model_labels.__getitem__, key=f"{key}-model-{input_index}",
    )
    selection = (models[selected].get("raw_metadata") or {}).get("legacy_model_selection") or {}
    if selection.get("rank_status") == "rounded_tie":
        st.caption("Native ranking is unavailable for this seed. The rounded ranking scores tie at first.")
    return models, model_labels, selected


def _plot_plddt_confidence_points(axis: Any, positions: Any, values: Any) -> None:
    """Plot pLDDT points using the same four confidence bands as Mol*."""
    import numpy as np

    scores = np.asarray(values, dtype=float)
    x_values = np.asarray(positions)
    valid = np.isfinite(scores)
    bands = (
        ("≤50", valid & (scores <= 50), "#ff7d45"),
        (">50–70", valid & (scores > 50) & (scores <= 70), "#ffdb13"),
        (">70–90", valid & (scores > 70) & (scores <= 90), "#65cbf3"),
        (">90", valid & (scores > 90), "#0053d6"),
    )
    for label, mask, color in bands:
        axis.scatter(x_values[mask], scores[mask], s=9, color=color, label=label, linewidths=0)
    axis.legend(
        title="pLDDT", loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=4,
        fontsize=8, title_fontsize=8, frameon=False,
    )


PAE_GREEN_COLORS = ("#246b3c", "#78b889", "#c6e3cd", "#f4faf5")


def _pyplot_pae_green_colormap() -> Any:
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("studio_pae_green", PAE_GREEN_COLORS)


def prediction_viewer_html(structure: Path, confidence: dict | None) -> str:
    from py2Dmol import view

    pae_preview = confidence.get("pae_preview") if confidence else None
    # Confidence arrays may be saved separately without atom B-factors.
    has_plddt = bool(confidence and confidence["plddt"].size and confidence["structure_plddt"])
    viewer = view(size=(560, 560), controls=True, box=True, color="deepmind" if has_plddt else "chain",
                  style="richardson", bg="white", heatmap=pae_preview is not None, heatmap_size=560)
    maps = ({"pae": {"data": pae_preview.tolist(), "vmin": 0.0, "vmax": 30.0, "colors": PAE_GREEN_COLORS}}
            if pae_preview is not None else None)
    with tempfile.TemporaryDirectory(prefix="studio-prediction-viewer-") as scratch:
        display_path = results.viewer_structure(structure, Path(scratch))
        viewer.add_pdb(str(display_path), name="prediction", maps=maps, use_biounit=False, filter_additives=False)
    return viewer.to_html(title="Prediction structure")


def render_confidence_plots(confidence: dict | None, *, unavailable: str = "") -> None:
    st.subheader("Confidence plots · selected model")
    if confidence is None:
        st.info(unavailable or "Confidence files are unavailable for the selected model.")
        return
    plddt = confidence["plddt"]
    pae = confidence["pae_preview"]
    if not plddt.size and pae is None:
        st.info("No pLDDT values or PAE matrix are available for the selected model.")
        return

    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(12, 4), constrained_layout=True)
    try:
        plddt_axis, pae_axis = figure.subplots(1, 2)
        if plddt.size:
            stride = max(1, int(np.ceil(plddt.size / 5000)))
            positions = np.arange(0, plddt.size, stride)
            _plot_plddt_confidence_points(plddt_axis, positions + 1, plddt[::stride])
            plddt_axis.set(xlabel=confidence.get("plddt_kind", "Atom"), ylabel="pLDDT",
                           ylim=(0, 100), title="Predicted pLDDT")
            for threshold in (50, 70, 90):
                plddt_axis.axhline(threshold, linestyle="--", linewidth=0.7, color="#94a39a", alpha=0.7)
        else:
            plddt_axis.axis("off")
            plddt_axis.text(.5, .5, "No pLDDT values available", ha="center", va="center")

        if pae is not None:
            heatmap = pae_axis.imshow(pae, cmap=_pyplot_pae_green_colormap(), vmin=0, vmax=30,
                                      interpolation="nearest", aspect="auto")
            pae_axis.set_box_aspect(1)
            figure.colorbar(heatmap, ax=pae_axis, fraction=0.046, pad=0.04, label="PAE (Å)")
            sampled = " (sampled)" if confidence["pae_sampled"] else ""
            pae_axis.set(xlabel="Scored position" + sampled, ylabel="Aligned position" + sampled,
                         title="Predicted aligned error")
        else:
            pae_axis.axis("off")
            pae_axis.text(.5, .5, "No PAE matrix saved", ha="center", va="center")
        st.pyplot(figure)
    finally:
        plt.close(figure)
