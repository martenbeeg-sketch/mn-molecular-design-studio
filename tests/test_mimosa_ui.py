from __future__ import annotations

from pathlib import Path

import pytest

from mn_molecular_design_studio.mimosa_input import normalize_mimosa_source_positions
from mn_molecular_design_studio import mimosa_ui
from mn_molecular_design_studio.mimosa_ui import (
    contacting_partner_residues,
    fetch_paper_mimosa_example,
    source_complex_chains,
)


def _atom_line(serial: int, atom: str, residue: str, chain: str, number: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    element = atom[0]
    return (
        f"ATOM  {serial:5d} {atom:>4} {residue:>3} {chain}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2}\n"
    )


def _write_complex(path: Path) -> None:
    lines = []
    serial = 1
    target_coordinates = [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 4.0, 0.0)]
    target_residues = ["GLY", "SER", "ALA"]
    for chain, residues, coordinates in (
        ("A", target_residues, target_coordinates),
        ("B", ["PRO", "TYR"], [(1.0, 1.0, 3.0), (4.0, 1.0, 3.0)]),
        ("C", ["LYS"], [(2.0, 2.0, 4.0)]),
    ):
        for number, (residue, ca) in enumerate(zip(residues, coordinates), start=1):
            x, y, z = ca
            for atom, atom_xyz in (
                ("N", (x - 1.0, y, z)),
                ("CA", ca),
                ("C", (x + 1.0, y, z)),
            ):
                number_offset = {"A": 0, "B": 20, "C": 30}[chain]
                lines.append(_atom_line(serial, atom, residue, chain, number + number_offset, atom_xyz))
                serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def test_normalize_mimosa_source_positions_accepts_author_ids() -> None:
    assert normalize_mimosa_source_positions("21, 22; 26a") == ["21", "22", "26A"]
    assert normalize_mimosa_source_positions(["-2", "105B"]) == ["-2", "105B"]
    assert normalize_mimosa_source_positions(["A:33", "b:94"]) == ["A:33", "B:94"]


@pytest.mark.parametrize("value", ["21,21", "bad", "21:22:3", ["22", "22"]])
def test_normalize_mimosa_source_positions_rejects_invalid_or_duplicate_ids(value) -> None:
    with pytest.raises(ValueError):
        normalize_mimosa_source_positions(value)


def test_mimosa_residue_picker_lists_contacting_partner_amino_acids(tmp_path: Path) -> None:
    complex_path = tmp_path / "complex.pdb"
    _write_complex(complex_path)

    assert source_complex_chains(complex_path) == ["A", "B", "C"]
    residues = contacting_partner_residues(complex_path, "B", "A", contact_cutoff=5.0)
    assert [row["position"] for row in residues] == ["21", "22"]
    assert [row["label"] for row in residues] == ["P21", "Y22"]


def test_mimosa_residue_picker_respects_contact_cutoff(tmp_path: Path) -> None:
    complex_path = tmp_path / "complex.pdb"
    _write_complex(complex_path)

    residues = contacting_partner_residues(complex_path, "B", "A", contact_cutoff=1.0)
    assert residues == []


def test_mimosa_residue_picker_qualifies_positions_across_partner_chains(tmp_path: Path) -> None:
    complex_path = tmp_path / "complex.pdb"
    _write_complex(complex_path)

    residues = contacting_partner_residues(complex_path, ["B", "C"], ["A"], contact_cutoff=5.0)

    assert [row["position"] for row in residues] == ["B:21", "B:22", "C:31"]
    assert [row["label"] for row in residues] == ["B:P21", "B:Y22", "C:K31"]


def test_keap1_nrf2_paper_example_loads_target_and_source_complex(tmp_path: Path, monkeypatch) -> None:
    fetched = []

    def fetch_structure(accession: str, extension: str, destination: str | Path) -> Path:
        fetched.append((accession, extension, Path(destination)))
        return Path(destination)

    def fetch_interaction(template_id: str, destination: str | Path) -> Path:
        fetched.append((template_id, "cif", Path(destination)))
        return Path(destination)

    monkeypatch.setattr(mimosa_ui, "_fetch_rcsb_structure", fetch_structure)
    monkeypatch.setattr(mimosa_ui, "fetch_paper_interaction_template", fetch_interaction)
    example = fetch_paper_mimosa_example("keap1_nrf2", tmp_path)

    assert example["target_id"] == "2FLU"
    assert example["target_chain"] == example["source_target_chain"] == "X"
    assert example["source_partner_chain"] == "P"
    assert example["condition"] == "sequence_fixed"
    assert [(accession, extension) for accession, extension, _ in fetched] == [
        ("2FLU", "pdb"),
        ("2FLU", "cif"),
    ]


def test_il8_paper_examples_use_scattered_multicdr_motifs() -> None:
    cxcr1 = mimosa_ui.PAPER_MIMOSA_EXAMPLES["il8_cxcr1"]
    antibody = mimosa_ui.PAPER_MIMOSA_EXAMPLES["il8_ly3041658"]

    assert (cxcr1["source_template"], cxcr1["source_partner_chain"], cxcr1["source_target_chain"]) == ("8IC0", "A", "F")
    assert cxcr1["motif_positions"] == ["21", "22", "23", "26", "27", "28", "29"]
    assert cxcr1["design_loop"] == "all_cdrs"
    assert (antibody["source_template"], antibody["source_partner_chain"], antibody["source_target_chain"]) == ("6WZM", "A,B", "E")
    assert antibody["motif_positions"] == ["A:33", "A:99", "A:102", "A:104", "A:105", "A:106", "A:110", "B:94"]
    assert antibody["design_loop"] == "all_cdrs"
