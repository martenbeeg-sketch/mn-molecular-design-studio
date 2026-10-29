from pathlib import Path

import pytest

from mn_molecular_design_studio import rfd4_proteina_input as inputs


def _parameters(target: Path):
    return {
        "task": "protein_binder",
        "input_mode": "condition_spec",
        "target_structure": str(target),
        "target_chains": ["A"],
        "input_files": [str(target)],
        "condition_spec": """idr_binder:
  input: R3.cif
  contig: "70-110, /0, A1-31"
  conditions:
    C_SEQ:
      true:
        - select: "sel('A')"
""",
        "num_designs": 1,
        "batch_size": 1,
        "nsteps": 60,
        "length_min": 70,
        "length_max": 110,
    }


def test_rfd4_campaign_accepts_one_native_condition_spec(tmp_path):
    target = tmp_path / "R3.cif"
    target.touch()
    assert inputs.validate_campaign_inputs(_parameters(target)) == 1


def test_rfd4_campaign_rejects_multi_design_condition_spec(tmp_path):
    target = tmp_path / "R3.cif"
    target.touch()
    parameters = _parameters(target)
    parameters["condition_spec"] += """
second_design:
  input: R3.cif
  contig: "20-30, /0, A1-31"
"""
    with pytest.raises(ValueError, match="exactly one named"):
        inputs.validate_campaign_inputs(parameters)
