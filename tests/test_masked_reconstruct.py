"""Offline tests for the masked-genome reconstruction harness (F1 engine).

Pure: uses MockFoundationModel (deterministic 6-mer masked-LM) + the Markov baseline. No torch, no
network, no weights. Pins the Markov math, the mock reconstruction, and the DELTA wiring (the headline
metric is LM-minus-Markov, never raw accuracy). Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.interp.markov_baseline import MarkovModel  # noqa: E402
from dna_decode.interp.masked_reconstruct import (  # noqa: E402
    _lm_base_counts,
    _masked_base_indices,
    score_reconstruction,
)
from dna_decode.models.foundation import MockFoundationModel, TokenPrediction  # noqa: E402


# ----- Markov baseline -------------------------------------------------------------------------
def test_markov_perfect_on_periodic():
    # In "ACAC...AC" order-1: after A always C, after C always A -> perfect on interior positions.
    seq = "AC" * 30
    m = MarkovModel.fit([seq], k=1)
    idx = list(range(2, len(seq)))  # skip the first period (no left context of order 1 at pos 0)
    assert m.accuracy_on_masked(seq, idx) == 1.0


def test_markov_backoff_and_tiebreak():
    m = MarkovModel.fit(["AAAAAA"], k=3)
    assert m.predict_base("AA") == "A"          # seen high-order context
    assert m.predict_base("GGGGG") == "A"       # unseen context -> backs off to mono (only A seen)
    empty = MarkovModel.fit([""], k=3)
    assert empty.predict_base("ACG") == "A"     # nothing learned -> deterministic "A"


def test_markov_skips_non_acgt():
    m = MarkovModel.fit(["ACGTN"], k=1)
    # index 4 is 'N' -> skipped; scored set is empty here
    assert m.accuracy_on_masked("ACGTN", [4]) == 0.0


# ----- mock masked-LM --------------------------------------------------------------------------
def test_mock_predictions_even_correct_odd_corrupt():
    model = MockFoundationModel()
    seq = "ACGTACGTACGTACGTACGTACGT"  # 24 bp -> four 6-mers
    preds = model.masked_token_predictions(seq)
    assert len(preds) == 4
    assert preds[0].pred_kmer == preds[0].true_kmer            # even -> correct
    assert preds[1].pred_kmer != preds[1].true_kmer            # odd -> corrupted
    assert preds[1].pred_kmer[1:] == preds[1].true_kmer[1:]    # only first base changed
    assert isinstance(preds[0], TokenPrediction)


def test_lm_base_counts_and_indices():
    model = MockFoundationModel()
    seq = "ACGTACGTACGTACGTACGTACGT"
    preds = model.masked_token_predictions(seq)
    correct, total = _lm_base_counts(preds)
    assert total == 24 and correct == 22          # 2 correct tokens (12) + 2 tokens missing 1 base (10)
    assert _masked_base_indices(preds) == list(range(24))


# ----- delta wiring ----------------------------------------------------------------------------
def test_score_reconstruction_delta_wiring():
    model = MockFoundationModel()
    seq = "ACGTACGTACGTACGTACGTACGT"
    score = score_reconstruction(model, seq, markov_k=3, region_label="unit")
    assert score.n_tokens_masked == 4
    assert score.n_bases_scored == 24
    assert score.lm_base_accuracy == pytest.approx(22 / 24)
    assert score.lm_token_accuracy == pytest.approx(0.5)
    assert score.lm_mean_true_prob == pytest.approx(0.5)      # (0.9+0.1+0.9+0.1)/4
    # delta is exactly LM minus the Markov baseline over the same base set
    m = MarkovModel.fit([seq], k=3)
    expected_markov = m.accuracy_on_masked(seq, list(range(24)))
    assert score.markov_base_accuracy == pytest.approx(expected_markov)
    assert score.delta == pytest.approx(22 / 24 - expected_markov)
    d = score.as_dict()
    assert d["headline_delta_lm_minus_markov"] == round(score.delta, 4)
    assert "delta" in d["note"] and "NOT the claim" in d["note"]


def test_non_mlm_model_raises():
    from dna_decode.models.foundation import FoundationModel, ModelMetadata

    class _Emb(FoundationModel):
        def _load_weights(self):
            pass

        def _embed_window(self, sequence):
            import numpy as np
            return np.zeros(4, dtype="float32")

    m = _Emb(ModelMetadata(name="emb", huggingface_id="x", embedding_dim=4, max_context=64))
    with pytest.raises(ValueError):
        score_reconstruction(m, "ACGTACGTACGT")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
