from __future__ import annotations

from pathlib import Path

_component = None


def render_hotspot_sequence(
    *,
    chain: str,
    sequence: list[tuple[int, str]],
    selected_residues: set[int],
    key: str,
):
    """Render the click-and-drag, 20-residue sequence selector."""
    global _component
    if _component is None:
        import streamlit.components.v1 as components

        component_path = Path(__file__).parent / "components" / "hotspot_sequence"
        _component = components.declare_component(
            "mn_molecular_design_studio_hotspot_sequence",
            path=str(component_path),
        )
    return _component(
        chainId=chain,
        sequence=[{"number": number, "letter": amino_acid} for number, amino_acid in sequence],
        selectedResidues=sorted(selected_residues),
        key=key,
        default=None,
    )
