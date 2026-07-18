"""Offline seam tests for the ProSST structure method (dna_decode/forward/prosst_scorer.py).

The real ProSST forward needs `transformers` + `prosst` (+ torch_geometric) and a GPU — NOT exercised in CI;
the real run is `scripts/prosst_lift.py` on Kaggle (Step 6 of the plan). These pin the seam: the
unavailable-signal, the predict_effect('prosst') branch on a stub table, the tier thresholds, the reused
AlphaFold URL, the hybrid composition, and the prosst_lift pure analysis helpers.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import prosst_tier, prosst_variant_table, quantize_structure  # noqa: E402  (exports)
from dna_decode.forward.prosst_scorer import (  # noqa: E402
    StructureMethodUnavailable, alphafold_pdb_url, prosst_model_name,
)
from dna_decode.forward.variant_effect import predict_effect, rank_average_hybrid  # noqa: E402


# ---- the unavailable-signal (this host has no prosst / torch_geometric) -----------------------------------

def test_quantize_structure_raises_when_prosst_absent():
    with pytest.raises(StructureMethodUnavailable, match="quantizer unavailable"):
        quantize_structure("nonexistent.pdb")


def test_load_prosst_raises_when_transformers_absent(monkeypatch):
    monkeypatch.setitem(sys.modules, "transformers", None)   # force `from transformers import ...` to fail
    from dna_decode.forward import prosst_scorer
    prosst_scorer._BUNDLE.clear()
    with pytest.raises(StructureMethodUnavailable, match="transformers"):
        prosst_scorer._load_prosst(2048)


def test_prosst_variant_table_needs_tokens_or_pdb():
    # reaches the input guard before any model load
    with pytest.raises((ValueError, StructureMethodUnavailable)):
        prosst_variant_table("MK", ["M1L"])


# ---- tiers + reused helpers -------------------------------------------------------------------------------

def test_prosst_tier_thresholds():
    assert prosst_tier(0.0) == "preserved"       # >= -1.0
    assert prosst_tier(-5.0) == "damaging"        # <= -3.0
    assert prosst_tier(-2.0) == "uncertain"       # in between


def test_alphafold_url_and_model_name():
    assert alphafold_pdb_url("P12345").endswith("AF-P12345-F1-model_v4.pdb")
    assert prosst_model_name(2048) == "AI4Protein/ProSST-2048"
    assert prosst_model_name(512) == "AI4Protein/ProSST-512"


# ---- the predict_effect seam (stub table -> no model) -----------------------------------------------------

def test_predict_effect_prosst_on_a_stub_table():
    p = predict_effect("M", "M1W", method="prosst", prosst_table={"M1W": -4.0})
    assert p.predicted_effect == "damaging" and p.regime == "B_molecular" and p.method == "prosst"
    q = predict_effect("M", "M1L", method="prosst", prosst_table={"M1L": 0.2})
    assert q.predicted_effect == "preserved"


def test_predict_effect_prosst_refuses_missing_table_or_variant():
    with pytest.raises(ValueError, match="prosst_table"):
        predict_effect("M", "M1W", method="prosst")                       # no table
    with pytest.raises(ValueError, match="prosst_table"):
        predict_effect("M", "M1W", method="prosst", prosst_table={"M1L": 0.1})   # variant absent


def test_unknown_method_message_names_prosst():
    with pytest.raises(NotImplementedError, match="prosst"):
        predict_effect("M", "M1L", method="magic")


def test_prosst_table_composes_into_the_hybrid():
    esm = {"M1L": -1.0, "M1W": -3.0, "M1V": 0.0}
    prosst = {"M1L": -0.5, "M1W": -4.0, "M1V": 0.1}
    combined = rank_average_hybrid([esm, prosst])
    assert set(combined) == {"M1L", "M1W", "M1V"}
    assert combined["M1V"] > combined["M1L"] > combined["M1W"]     # concordant order preserved


# ---- prosst_variant_table with a FAKE HF bundle (offline; never loads/downloads the real ProSST) ----------
#
# The real forward pass needs transformers+prosst+GPU (Step-6 Kaggle). A fake (model, tokenizer) bundle
# mimics the HF interface with real torch tensors, so the model_bundle reuse path, the token-length guard,
# the mutant-skip branches, and the log-ratio math are all exercised WITHOUT touching the network or the
# 50s model download. Fake logits are `[0,1,...,19]` at every position, so log_softmax cancels the constant
# and the score is EXACTLY idx(alt) − idx(wt) over "ACDEFGHIKLMNPQRSTVWY" — a deterministic, checkable sign.

_VOCAB_AA = "ACDEFGHIKLMNPQRSTVWY"


class _FakeOut:
    def __init__(self, logits):
        self.logits = logits


class _FakeModel:
    def __call__(self, input_ids, attention_mask=None, ss_input_ids=None):
        import torch
        L = input_ids.shape[1]
        row = torch.arange(len(_VOCAB_AA), dtype=torch.float)   # logit(aa) = its index
        return _FakeOut(row.repeat(L, 1).unsqueeze(0))          # (1, L, 20)


class _FakeTokenizer:
    def __call__(self, seq, return_tensors="pt"):
        import torch
        L = len(seq) + 2                                        # CLS + residues + EOS
        return {"input_ids": torch.zeros((1, L), dtype=torch.long),
                "attention_mask": torch.ones((1, L), dtype=torch.long)}

    def get_vocab(self):
        return {aa: i for i, aa in enumerate(_VOCAB_AA)}


def _fake_bundle():
    return (_FakeModel(), _FakeTokenizer())


def test_prosst_variant_table_length_mismatch_guard():
    # tokens (3) != sequence length (2): fires AFTER bundle binding but BEFORE any forward pass
    with pytest.raises(ValueError, match="structure tokens"):
        prosst_variant_table("WA", ["W1A"], structure_tokens=[1, 2, 3], model_bundle=_fake_bundle())


def test_prosst_variant_table_scores_and_skips_via_fake_bundle():
    from dna_decode.forward import prosst_scorer
    prosst_scorer._BUNDLE.clear()
    seq = "WA"                                                  # W=idx18, A=idx0
    mutants = [
        "W1A",          # valid, damaging direction: idx(A)-idx(W) = 0-18 = -18
        "A2W",          # valid, preserved direction: idx(W)-idx(A) = 18-0 = +18
        "A2A",          # valid, self-sub: 0
        "W1A:A2W",      # multi-mutant (has ':') -> skipped
        "XX",           # too short -> skipped
        "W1Z",          # Z not a standard AA -> skipped
        "K1A",          # WT mismatch (seq[0]='W' != 'K') -> skipped
        "W9A",          # position past length -> skipped
    ]
    table = prosst_variant_table(seq, mutants, structure_tokens=[3, 5], model_bundle=_fake_bundle())
    assert set(table) == {"W1A", "A2W", "A2A"}                  # only the valid singles survive
    assert table["W1A"] == -18.0 and table["A2W"] == 18.0 and table["A2A"] == 0.0
    assert prosst_tier(table["W1A"]) == "damaging" and prosst_tier(table["A2W"]) == "preserved"
    # model_bundle reuse path: the real loader was never invoked, so nothing got cached/downloaded
    assert prosst_scorer._BUNDLE == {}


def test_prosst_variant_table_pdb_path_quantizes_and_raises_when_absent():
    # structure_tokens omitted + pdb_path given: reaches quantize_structure (no prosst stack) -> unavailable
    with pytest.raises(StructureMethodUnavailable, match="quantizer unavailable"):
        prosst_variant_table("WA", ["W1A"], pdb_path="nonexistent.pdb", model_bundle=_fake_bundle())


# ---- prosst_lift.py pure analysis helpers -----------------------------------------------------------------

def test_prosst_lift_pure_helpers():
    from scripts.prosst_lift import paired_in_subset, by_category, reproduction_median
    recs = [
        {"dms": "A_Activity", "status": "OK", "hybrid_minus_esm": 0.05, "prosst_vs_pg_repro": 0.9},
        {"dms": "B_Stability", "status": "OK", "hybrid_minus_esm": 0.07, "prosst_vs_pg_repro": 0.8},
        {"dms": "C_Stability", "status": "OK", "hybrid_minus_esm": -0.01, "prosst_vs_pg_repro": 0.85},
        {"dms": "D_bad", "status": "ERROR"},
    ]
    cats = {"A_Activity": "Activity", "B_Stability": "Stability", "C_Stability": "Stability"}
    p = paired_in_subset(recs)
    assert p["n"] == 3 and p["win"] == "2/3"
    bc = by_category(recs, cats)
    assert bc["Stability"]["n"] == 2 and bc["Activity"]["n"] == 1
    assert reproduction_median(recs) == 0.85


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
