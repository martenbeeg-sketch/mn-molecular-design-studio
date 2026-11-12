import builtins
from contextlib import nullcontext
import json
import math
import sys
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows import quality_runner


@pytest.fixture
def native_esm(monkeypatch):
    loaded, sequences = [], []
    control = {"probability": None}

    class Tensor:
        def __init__(self, values):
            self.values = values

        def __getitem__(self, key):
            return self

        def to(self, device):
            return self

        def softmax(self, axis):
            return self

        def gather(self, axis, indices):
            return self

        def squeeze(self, axis):
            return self

        def cpu(self):
            return self

        def tolist(self):
            return self.values

    class Model:
        def eval(self):
            return self

        def to(self, device):
            return self

        def __call__(self, tokens):
            value = control["probability"]
            if value is None:
                value = len(loaded) / 10
            return {"logits": Tensor([value] * len(tokens.values))}

    def convert(batch):
        sequence = batch[0][1]
        sequences.append(sequence)
        return None, None, Tensor(sequence)

    def load(path):
        loaded.append(path)
        return Model(), SimpleNamespace(get_batch_converter=lambda: convert, get_idx=lambda amino_acid: 0)

    monkeypatch.setitem(sys.modules, "esm", SimpleNamespace(pretrained=SimpleNamespace(load_model_and_alphabet_local=load)))
    monkeypatch.setitem(sys.modules, "torch", SimpleNamespace(manual_seed=lambda seed: None,
        cuda=SimpleNamespace(is_available=lambda: False), inference_mode=nullcontext,
        tensor=lambda values, **kwargs: Tensor(values)))
    original_import = builtins.__import__

    def without_biopython(name, *args, **kwargs):
        if name == "Bio" or name.startswith("Bio."):
            raise ModuleNotFoundError("No module named 'Bio'")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_biopython)
    return loaded, sequences, control


def test_five_model_native_scoring_does_not_require_biopython(tmp_path, native_esm):
    loaded, sequences, _ = native_esm
    item = {"id": "valid", "sequence": "acdef"}
    results = quality_runner.esm_scores([item], "esm1v_score", tmp_path, lambda *args, **kwargs: None)
    assert results[0]["status"] == "completed"
    assert results[0]["metrics"]["esm1v_native_probability"] == pytest.approx(0.3)
    assert [path.rsplit("_", 1)[1] for path in loaded] == [f"{number}.pt" for number in range(1, 6)]
    assert sequences == ["ACDEF"] * 5
    saved = json.loads((tmp_path / "valid/native_probabilities.json").read_text())
    assert saved["probabilities"] == pytest.approx([0.3] * 5)
    assert saved["sequence"] == "ACDEF"
    assert saved["protocol"] == "five-model unmasked wild-type marginals"
    assert item["sequence"] == "acdef"


@pytest.mark.parametrize("sequence,error", [("", "nonempty"), ("ACDX", "unknown residues"),
                                              ("A" * 1023, "1,022")])
def test_invalid_sequences_fail_before_loading_models(tmp_path, native_esm, sequence, error):
    results = quality_runner.esm_scores([{"id": "bad", "sequence": sequence}], "esm1v_score", tmp_path, lambda *args, **kwargs: None)
    assert results[0]["status"] == "failed"
    assert error in results[0]["error"]
    assert not native_esm[0]


@pytest.mark.parametrize("probability", [math.nan, math.inf, -0.1, 1.1])
def test_nonfinite_or_out_of_range_model_probabilities_cannot_be_success(tmp_path, native_esm, probability):
    native_esm[2]["probability"] = probability
    results = quality_runner.esm_scores([{"id": "bad", "sequence": "ACDE"}], "esm1v_score", tmp_path, lambda *args, **kwargs: None)
    assert results[0]["status"] == "failed"
    assert results[0]["metrics"] == {}
    assert not (tmp_path / "bad/native_probabilities.json").exists()


def test_one_invalid_item_does_not_block_valid_items(tmp_path, native_esm):
    results = quality_runner.esm_scores([{"id": "bad", "sequence": "X"}, {"id": "valid", "sequence": "ACDE"}],
                                      "esm1v_score", tmp_path, lambda *args, **kwargs: None)
    assert {row["id"]: row["status"] for row in results} == {"bad": "failed", "valid": "completed"}


def test_input_error_is_not_scored(tmp_path, native_esm):
    results = quality_runner.esm_scores([{"id": "bad", "sequence": "ACDE", "input_error": "unverified structure"}],
                                      "esm1v_score", tmp_path, lambda *args, **kwargs: None)
    assert results[0]["status"] == "failed"
    assert "unverified structure" in results[0]["error"]
    assert not native_esm[0]
