"""Offline tests for the masked-genome reconstruction harness (F1/F1' engine).

Pure: MockFoundationModel (deterministic 6-mer masked-LM with per-base marginals) + the Markov
baseline. No torch, no network, no weights. Pins the Markov math (incl. leave-one-out + smoothed
per-base NLL), the mock reconstruction, and the DELTA wiring (PRIMARY = per-base NLL delta, NT
marginals vs disjoint/LOO Markov). Runnable via pytest OR standalone.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.interp.markov_baseline import MarkovModel  # noqa: E402
from dna_decode.interp.masked_reconstruct import (  # noqa: E402
    _masked_base_indices,
    _nt_per_base,
    score_reconstruction,
)
from dna_decode.models.foundation import MockFoundationModel, TokenPrediction  # noqa: E402


# ----- Markov baseline -------------------------------------------------------------------------
def test_markov_perfect_on_periodic():
    seq = "AC" * 30
    m = MarkovModel.fit([seq], k=1)
    idx = list(range(2, len(seq)))
    assert m.accuracy_on_masked(seq, idx) == 1.0


def test_markov_base_distribution_normalized_and_smoothed():
    m = MarkovModel.fit(["ACGTACGTAC"], k=2)
    d = m.base_distribution("AC", alpha=1.0)
    assert abs(sum(d.values()) - 1.0) < 1e-9
    assert all(v > 0 for v in d.values())          # Laplace smoothing -> no zeros
    u = m.base_distribution("ZZZ")                  # unseen context -> uniform-ish (smoothed)
    assert abs(sum(u.values()) - 1.0) < 1e-9


def test_markov_leave_one_out_excludes_self_count():
    # context "AA" is followed by G exactly once. Without LOO, predict_base("AA") == 'G'.
    # With LOO excluding G, that lone count is removed -> backs off / no longer picks G.
    m = MarkovModel.fit(["AAG"], k=2)
    assert m.predict_base("AA") == "G"
    assert m.predict_base("AA", exclude_base="G") != "G"
    # LOO also raises the NLL of the leaked target
    seq = "AAG"
    nll_leak = m.nll_on_masked(seq, [2], leave_one_out=False)
    nll_loo = m.nll_on_masked(seq, [2], leave_one_out=True)
    assert nll_loo > nll_leak


def test_markov_nll_finite_and_positive():
    m = MarkovModel.fit(["ACGT" * 20], k=3)
    nll = m.nll_on_masked("ACGT" * 20, list(range(10, 40)), leave_one_out=True)
    assert 0.0 < nll < 10.0


# ----- mock masked-LM (argmax + per-base marginals) --------------------------------------------
def test_mock_predictions_and_marginals():
    model = MockFoundationModel()
    seq = "ACGTACGTACGTACGTACGTACGT"  # four 6-mers
    preds = model.masked_token_predictions(seq)
    assert len(preds) == 4 and isinstance(preds[0], TokenPrediction)
    assert preds[0].pred_kmer == preds[0].true_kmer                 # even -> correct token
    assert preds[1].pred_kmer[0] != preds[1].true_kmer[0]           # odd -> first base wrong
    for p in preds:                                                  # marginals normalized
        for marg in p.base_marginals:
            assert abs(sum(marg) - 1.0) < 1e-9
    assert _masked_base_indices(preds) == list(range(24))


def test_nt_per_base_nll_and_accuracy_from_mock():
    model = MockFoundationModel()
    seq = "ACGTACGTACGTACGTACGTACGT"
    preds = model.masked_token_predictions(seq)
    nll, acc, n = _nt_per_base(preds)
    assert n == 24
    assert acc == pytest.approx(22 / 24)                            # 2 even (12) + 2 odd (10)
    # 22 correct bases at -log(0.9) + 2 wrong (odd offset0) at -log(0.1)
    expected = (22 * -math.log(0.9) + 2 * -math.log(0.1)) / 24
    assert nll == pytest.approx(expected)


# ----- delta wiring (PRIMARY = NLL delta) ------------------------------------------------------
def test_score_reconstruction_nll_delta_wiring_loo_default():
    model = MockFoundationModel()
    seq = "ACGTACGTACGTACGTACGTACGT"
    score = score_reconstruction(model, seq, markov_k=3, region_label="unit")
    assert score.markov_mode == "leave-one-out"                    # no fit seqs -> LOO default
    assert score.n_bases_scored == 24
    m = MarkovModel.fit([seq], k=3)
    exp_markov_nll = m.nll_on_masked(seq, list(range(24)), leave_one_out=True)
    exp_nt_nll = _nt_per_base(model.masked_token_predictions(seq))[0]
    assert score.markov_per_base_nll == pytest.approx(exp_markov_nll)
    assert score.nt_per_base_nll == pytest.approx(exp_nt_nll)
    assert score.nll_delta == pytest.approx(exp_markov_nll - exp_nt_nll)
    d = score.as_dict()
    assert d["PRIMARY_nll_delta_markov_minus_nt"] == round(score.nll_delta, 4)
    assert "identical masked bases" in d["note"]


def test_score_reconstruction_disjoint_mode():
    model = MockFoundationModel()
    seq = "ACGTACGTACGTACGTACGTACGT"
    train = "GGGGCCCCAAAATTTT" * 4                                  # disjoint training region
    score = score_reconstruction(model, seq, markov_k=3, markov_fit_sequences=[train],
                                 region_label="unit")
    assert score.markov_mode == "disjoint-fit"


def test_strict_raises_on_dropped_position():
    from dna_decode.models.foundation import FoundationModelError
    model = MockFoundationModel()
    seq = "ACGTACGTACGT"  # two 6-mer tokens (indices 0,1)
    with pytest.raises(FoundationModelError):
        model.masked_token_predictions(seq, positions=[0, 99], strict=True)  # 99 out of range


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
