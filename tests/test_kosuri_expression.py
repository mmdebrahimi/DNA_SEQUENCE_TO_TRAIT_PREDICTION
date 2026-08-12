"""Track B expression validation — pure logic (wheel-only) + the real-data gate when sd03 is present."""
from __future__ import annotations

import os

import numpy as np
import pytest

from scripts.kosuri_expression_validate import PREREGISTERED_BAR, additive_predict, r2, verdict

_SD03 = os.environ.get("KOSURI_SD03", "D:/PythonProjects/DNA_AI_Decoder/sd03.xls")


# ---- pure: the scoring + baseline ----

def test_r2_is_1_for_a_perfect_fit_and_0_for_the_mean():
    y = np.array([1.0, 2.0, 3.0, 4.0])
    assert r2(y, y) == pytest.approx(1.0)
    assert r2(y, np.full_like(y, y.mean())) == pytest.approx(0.0)


def test_r2_goes_negative_for_a_predictor_worse_than_the_mean():
    """Load-bearing: the headline finding is a NEGATIVE R2 (-0.014 on held-out promoters). If r2 clipped
    at 0 that result would silently read as 'no signal' instead of 'worse than predicting the mean'."""
    y = np.array([1.0, 2.0, 3.0, 4.0])
    assert r2(y, np.array([4.0, 3.0, 2.0, 1.0])) < 0


def test_r2_ignores_non_finite_pairs():
    y = np.array([1.0, 2.0, np.nan, 4.0])
    assert np.isfinite(r2(y, np.array([1.0, 2.0, 3.0, 4.0])))


def _frame(rows):
    pd = pytest.importorskip("pandas")
    return pd.DataFrame(rows)


def test_additive_baseline_recovers_a_clean_additive_signal():
    """y = mu + promoter effect + RBS effect must be recovered exactly when that IS the generating model."""
    rows = [{"p_code": p, "r_code": r, "y": 10.0 + p - r} for p in range(3) for r in range(3)]
    d = _frame(rows)
    assert r2(d.y.values, additive_predict(d, d)) == pytest.approx(1.0, abs=1e-9)


def test_unseen_element_falls_back_rather_than_erroring():
    """An unseen element contributes 0, so the baseline degrades to (other element + grand mean) instead
    of crashing. That fallback is exactly why it still scores 0.26-0.50 where an identity model scores ~0."""
    tr = _frame([{"p_code": p, "r_code": r, "y": 10.0 + p - r} for p in range(3) for r in range(3)])
    te = _frame([{"p_code": 99, "r_code": 0, "y": 12.0}])          # promoter never seen in training
    pred = additive_predict(tr, te)
    assert np.isfinite(pred).all()


# ---- pure: the verdict contract ----

def test_verdict_reports_both_halves_and_they_can_disagree():
    """The finding IS the disagreement: PASS on combinations, FAIL on elements. A verdict function that
    collapsed to one boolean would erase it."""
    v = verdict({
        "held_out_combination": {"additive_baseline": 0.795, "gbm_identity": 0.893,
                                 "gbm_identity_plus_deltaG": 0.919},
        "held_out_promoter": {"additive_baseline": 0.263, "gbm_identity": -0.014,
                              "gbm_identity_plus_deltaG": 0.144},
        "held_out_rbs": {"additive_baseline": 0.499, "gbm_identity": 0.268,
                         "gbm_identity_plus_deltaG": 0.327},
    })
    assert v["combination_split_verdict"] == "PASS"
    assert v["element_split_verdict"] == "FAIL"
    assert v["preregistered_bar"] == PREREGISTERED_BAR


def test_a_combination_result_below_the_bar_is_reported_as_FAIL():
    v = verdict({
        "held_out_combination": {"additive_baseline": 0.79, "gbm_identity": 0.80,
                                 "gbm_identity_plus_deltaG": 0.81},
        "held_out_promoter": {"additive_baseline": 0.2, "gbm_identity": 0.0,
                              "gbm_identity_plus_deltaG": 0.1},
        "held_out_rbs": {"additive_baseline": 0.4, "gbm_identity": 0.2,
                         "gbm_identity_plus_deltaG": 0.3},
    })
    assert v["combination_split_verdict"] == "FAIL"


# ---- pure: the RBS sequence features (the design-question arm) ----

def test_sd_motif_scores_highest_on_a_perfect_shine_dalgarno():
    from scripts.kosuri_expression_validate import rbs_features

    perfect = rbs_features("TTTTAGGAGGTTTTATG")[-2]     # best SD match
    none = rbs_features("TTTTTTTTTTTTTTTATG")[-2]
    assert perfect == 6 and none < 6


def test_sd_spacing_is_measured_from_the_motif_to_the_end():
    from scripts.kosuri_expression_validate import rbs_features

    near = rbs_features("TTTTTTAGGAGGATG")[-1]
    far = rbs_features("AGGAGGTTTTTTTTTTTTATG")[-1]
    assert far > near > 0


def test_features_use_only_the_sequence_and_are_fixed_width():
    """No identity, no measured strength -- a held-out RBS must be scored from its letters alone."""
    from scripts.kosuri_expression_validate import rbs_features

    a, b = rbs_features("ACGTACGTACGTAAGGAGGATG"), rbs_features("TTTT")
    assert len(a) == len(b) == 4 + 16 + 64 + 4
    assert rbs_features("ACGT") == rbs_features("acgt")     # case-insensitive


def _split(head, base=0.4991, ident=0.2678, other=0.4991, seqonly=0.1429, ridge=0.0681, oracle=0.8068):
    def s(m): return {"mean": m, "std": 0.03, "p5": m - 0.05, "p95": m + 0.05}
    return {"additive_baseline": s(base), "identity": s(ident), "other_element_only": s(other),
            "sequence_only": s(seqonly), "other_plus_sequence": s(head),
            "ridge_other_plus_sequence": s(ridge),
            "other_plus_sequence_plus_deltaG_ORACLE": s(oracle), "n_splits": 25}


def test_verdict_headlines_the_NO_deltaG_arm():
    """Load-bearing. deltaG spans promoter TSS -> +30 GFP, so it is not design-time recomputable; it
    must never be the headline. A previous version headlined it (0.781) and overstated the result."""
    from scripts.kosuri_expression_validate import sequence_verdict

    v = sequence_verdict(_split(0.7762), {"r2": 0.6123, "n_elements": 111}, "rbs")
    assert v["headline_from_sequence"] == 0.7762          # the no-deltaG arm
    assert v["deltaG_oracle_upper_bound"] == 0.8068       # reported, but separately
    assert v["headline_from_sequence"] != v["deltaG_oracle_upper_bound"]
    assert "not design-time recomputable" in v["oracle_note"].lower() or "oracle" in v["oracle_note"].lower()


def test_verdict_requires_a_real_margin_over_the_baseline():
    """A sequence model that merely ties the baseline is NOT generalisation."""
    from scripts.kosuri_expression_validate import sequence_verdict

    real = sequence_verdict(_split(0.7762), {"r2": 0.6123, "n_elements": 111}, "rbs")
    assert real["generalises_from_sequence"] is True
    assert real["vs_additive_baseline"] == pytest.approx(0.2771, abs=1e-3)

    tie = sequence_verdict(_split(0.51), {"r2": 0.2, "n_elements": 111}, "rbs")
    assert tie["generalises_from_sequence"] is False


def test_verdict_carries_the_conditional_scope_and_both_numbers():
    """The per-construct and per-element numbers answer different questions; both must survive."""
    from scripts.kosuri_expression_validate import sequence_verdict

    v = sequence_verdict(_split(0.7762), {"r2": 0.6123, "n_elements": 111}, "rbs")
    assert v["per_element_mean_r2"] == 0.6123 and v["per_element_n"] == 111
    assert "characterised" in v["scope"].lower()


def test_promoter_features_score_the_sigma70_boxes():
    from scripts.kosuri_expression_validate import promoter_features

    f = promoter_features("AAATTGACAGGGGGGGGGGGGGGGGGTATAATAAA")
    assert f[-4] == 6 and f[-3] == 6      # perfect -35 (TTGACA) and -10 (TATAAT) matches
    weak = promoter_features("A" * 40)
    assert weak[-4] < 6 and weak[-3] < 6


def test_promoter_features_exclude_measured_TSS_and_are_fixed_width():
    """TSS.best is MEASURED by RNA-seq in this dataset -- including it would repeat the deltaG mistake."""
    from scripts.kosuri_expression_validate import promoter_features

    a, b = promoter_features("ACGTTTGACAAAGGTATAATGG"), promoter_features("TTTT")
    assert len(a) == len(b) == 4 + 16 + 64 + 6
    assert promoter_features("acgt") == promoter_features("ACGT")


# ---- real data (slow, skipped when the uncommitted supplementary is absent) ----

@pytest.mark.slow
@pytest.mark.skipif(not os.path.exists(_SD03), reason="Kosuri sd03.xls not present (not committed)")
def test_reproduces_the_papers_own_published_numbers():
    """Gate before trusting anything downstream: recompute the paper's results from its own columns.

    Also pins the units trap -- `model.prot.simple` is log2 while `prot` is raw; comparing in the wrong
    space returns R2 = -15.
    """
    pytest.importorskip("xlrd")
    from scripts.kosuri_expression_validate import reproduce_published

    rep = reproduce_published(_SD03)
    assert rep["rna_simple_log10"] == pytest.approx(0.92, abs=0.01)
    assert rep["rna_full_log10"] == pytest.approx(0.96, abs=0.01)
    assert rep["protein_simple_log2"] == pytest.approx(0.76, abs=0.02)


# ---- pure: library provenance + the OOD stress test ----

def test_library_is_recovered_from_the_name_prefix():
    """Leave-library-out is only possible because provenance survives in the part names."""
    from scripts.kosuri_expression_validate import library_of

    assert library_of("apFAB871") == "BIOFAB"
    assert library_of("BBa_J61133") == "BioBrick/Anderson"
    assert library_of("J23101") == "BioBrick/Anderson"      # Anderson promoter series
    assert library_of("salis-3-11") == "Salis"
    assert library_of("lacUV5") == "vector/other"
    assert library_of("DeadRBS") == "vector/other"


def test_library_classifier_is_total_and_never_raises():
    from scripts.kosuri_expression_validate import library_of

    for odd in ("", "   ", "x", "APFAB871"):     # case-sensitive by design; must still return a bucket
        assert library_of(odd) in {"BIOFAB", "BioBrick/Anderson", "Salis", "vector/other"}


def test_control_percentile_is_empirical_and_assumes_no_distribution_shape():
    """Replaces a `mean - 2*sd` heuristic that compared ONE structured point to a control spread while
    assuming normality and correcting for nothing. The percentile assumes no shape at all."""
    from scripts.kosuri_expression_validate import control_percentile

    ctl = [0.1 * i for i in range(10)]            # 0.0 .. 0.9
    assert control_percentile(-1.0, ctl) == 0.0   # worse than every random split
    assert control_percentile(1.0, ctl) == 1.0    # better than every random split
    assert control_percentile(0.35, ctl) == pytest.approx(0.4)


def test_control_percentile_survives_non_finite_controls_and_never_raises():
    """A degenerate same-size split can return nan; it must be dropped, not poison the percentile."""
    from scripts.kosuri_expression_validate import control_percentile

    assert control_percentile(0.5, [0.1, float("nan"), 0.9]) == pytest.approx(0.5)
    assert np.isnan(control_percentile(0.5, [float("nan")]))
    assert np.isnan(control_percentile(float("nan"), [0.1, 0.2]))


def test_within_group_r2_is_zero_for_a_model_that_only_knows_the_group():
    """THE load-bearing metric. A model predicting each library's mean scores well against the global mean
    while having no part-level ranking at all -- exactly the failure the global-mean denominator hid
    (RBS-Salis: 0.625 offset-inclusive, 0.100 within-library)."""
    from scripts.kosuri_expression_validate import _within_group_r2, r2

    y = np.array([1.0, 2.0, 3.0, 11.0, 12.0, 13.0])
    groups = np.array(["a", "a", "a", "b", "b", "b"])
    group_mean_only = np.array([2.0, 2.0, 2.0, 12.0, 12.0, 12.0])

    assert r2(y, group_mean_only) > 0.85                            # looks strong globally
    assert _within_group_r2(y, group_mean_only, groups) == pytest.approx(0.0)   # ranks nothing


def test_within_group_r2_goes_negative_when_ranking_is_inverted():
    """Two of three promoter libraries came out NEGATIVE within-library; the metric must be able to say so."""
    from scripts.kosuri_expression_validate import _within_group_r2

    y = np.array([1.0, 2.0, 3.0])
    assert _within_group_r2(y, np.array([3.0, 2.0, 1.0]), np.zeros(3)) < 0


def test_rmse_and_spearman_are_denominator_free_and_rank_only():
    from scripts.kosuri_expression_validate import _rmse, _spearman

    y = np.array([1.0, 2.0, 3.0, 4.0])
    assert _rmse(y, y) == pytest.approx(0.0)
    assert _rmse(y, y + 2.0) == pytest.approx(2.0)
    assert _spearman(y, 10 * y + 5) == pytest.approx(1.0)      # monotone => rank-identical
    assert _spearman(y, -y) == pytest.approx(-1.0)


def test_sign_agreement_is_what_makes_a_holdout_verdict_robust():
    """The rule that caught a wrong published number. Promoter-BIOFAB scored -0.6547 under the default
    regressor (which degenerates to a constant at n_train=22) and POSITIVE under three others -- the signs
    disagree, so the verdict was the model, not the data."""
    agree = lambda v: bool(all(x > 0 for x in v) or all(x < 0 for x in v))  # noqa: E731

    assert agree([-0.6547, 0.1220, 0.0952, 0.0144]) is False   # promoter BIOFAB: model-dependent
    assert agree([0.2524, 0.3344, 0.4709, 0.3750]) is True     # RBS BIOFAB: robust


def test_degenerate_constant_prediction_is_detectable():
    """A regressor that cannot split emits ONE value; its within-library R2 is then 0.0 BY CONSTRUCTION and
    means nothing. That must be visible as degeneracy rather than read as 'no signal'."""
    from scripts.kosuri_expression_validate import _within_group_r2

    y = np.array([1.0, 2.0, 3.0, 4.0])
    constant = np.full(4, 13.1199)
    assert len(np.unique(np.round(constant, 9))) == 1
    assert _within_group_r2(y, constant, np.zeros(4)) == pytest.approx(0.0)


def test_leave_library_out_defaults_to_enough_controls_for_a_percentile():
    """20 controls cannot resolve a percentile below 0.05; the two load-bearing verdicts sit at 0.005/0.000."""
    import inspect

    from scripts.kosuri_expression_validate import leave_library_out_with_size_control as f

    assert inspect.signature(f).parameters["n_control"].default >= 200
