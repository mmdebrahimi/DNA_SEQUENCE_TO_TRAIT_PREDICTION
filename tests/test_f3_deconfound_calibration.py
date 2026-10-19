"""F3 Step 4: the calibration harness's frozen verdict machinery (offline).

Every verdict state has a test that REACHES it -- a state with no such test is decoration.
The ordering is pinned, because a case satisfying two conditions must return the earlier one.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import numpy as np

from scripts.f3_deconfound_calibration import (
    MAX_TRAITS,
    MIN_TRAITS_FOR_SIGN_TEST,
    N_PERM,
    VERDICT_MEANING,
    VERDICTS,
    _phenotype_traits,
    _select_controls,
    compute_verdict,
    sign_test,
)
from scripts.f3_substrate_census import PANELS

REPO = Path(__file__).resolve().parent.parent


def _arm(n_traits=8, median_drop=0.01, grouped_beats=True):
    return {"n_traits": n_traits, "median_drop": median_drop,
            "grouped_beats_own_null": grouped_beats}


def _args(**kw):
    base = dict(n_intact=4, has_negative=True, has_positive=True,
                positive_fold_compatible=True, neg=_arm(), pos=_arm(median_drop=0.10))
    base.update(kw)
    return base


# --------------------------------------------------------------------------- every state reachable

def test_blocked_no_intact_panel():
    assert compute_verdict(**_args(n_intact=0)) == "blocked:no_intact_panel"


def test_blocked_no_negative_control():
    assert compute_verdict(**_args(has_negative=False)) == "blocked:no_negative_control"


def test_blocked_no_positive_control():
    assert compute_verdict(**_args(has_positive=False)) == "blocked:no_positive_control"


def test_blocked_positive_control_not_fold_compatible():
    assert (compute_verdict(**_args(positive_fold_compatible=False))
            == "blocked:positive_control_not_fold_compatible")


def test_indeterminate_underpowered_on_either_control():
    assert compute_verdict(**_args(neg=_arm(n_traits=MIN_TRAITS_FOR_SIGN_TEST - 1))) \
        == "INDETERMINATE_UNDERPOWERED"
    assert compute_verdict(**_args(pos=_arm(n_traits=MIN_TRAITS_FOR_SIGN_TEST - 1))) \
        == "INDETERMINATE_UNDERPOWERED"


def test_indeterminate_when_an_arm_is_missing_entirely():
    assert compute_verdict(**_args(neg=None)) == "INDETERMINATE_UNDERPOWERED"


def test_invalid_over_aggressive_when_the_negative_control_loses_its_null():
    """The de-confounder destroys signal where there is no structure to remove."""
    assert compute_verdict(**_args(neg=_arm(grouped_beats=False))) == "INVALID_OVER_AGGRESSIVE"


def test_invalid_insensitive_when_the_positive_drop_does_not_exceed_the_negative():
    assert compute_verdict(**_args(pos=_arm(median_drop=0.01), neg=_arm(median_drop=0.01))) \
        == "INVALID_INSENSITIVE"
    assert compute_verdict(**_args(pos=_arm(median_drop=-0.05), neg=_arm(median_drop=0.01))) \
        == "INVALID_INSENSITIVE"


def test_valid_differential_response():
    assert compute_verdict(**_args()) == "VALID_DIFFERENTIAL_RESPONSE"


def test_every_declared_verdict_is_reachable_and_documented():
    """No decoration: each state is produced by some input above, and each has a meaning."""
    produced = {
        compute_verdict(**_args(n_intact=0)),
        compute_verdict(**_args(has_negative=False)),
        compute_verdict(**_args(has_positive=False)),
        compute_verdict(**_args(positive_fold_compatible=False)),
        compute_verdict(**_args(neg=_arm(n_traits=1))),
        compute_verdict(**_args(neg=_arm(grouped_beats=False))),
        compute_verdict(**_args(pos=_arm(median_drop=0.0))),
        compute_verdict(**_args()),
    }
    assert produced == set(VERDICTS)
    for v in VERDICTS:
        assert v in VERDICT_MEANING and len(VERDICT_MEANING[v]) > 20


# --------------------------------------------------------------------------- ordering

def test_blocked_outranks_underpowered():
    """A blocked run must not be reported as merely underpowered."""
    assert compute_verdict(**_args(positive_fold_compatible=False, neg=_arm(n_traits=1))) \
        == "blocked:positive_control_not_fold_compatible"


def test_over_aggression_outranks_insensitivity():
    """Both failures at once reports the one that says the instrument is actively harmful."""
    assert compute_verdict(**_args(neg=_arm(grouped_beats=False, median_drop=0.5),
                                   pos=_arm(median_drop=0.0))) == "INVALID_OVER_AGGRESSIVE"


def test_no_intact_panel_outranks_everything():
    assert compute_verdict(**_args(n_intact=0, has_negative=False, has_positive=False,
                                   positive_fold_compatible=False)) == "blocked:no_intact_panel"


# --------------------------------------------------------------------------- sign test

def test_sign_test_is_paired_and_two_sided():
    assert sign_test([0.1] * 6)["p"] == pytest.approx(2 / 64)
    assert sign_test([0.1, -0.1])["p"] == pytest.approx(1.0)


def test_sign_test_ignores_exact_zeros_in_the_denominator():
    out = sign_test([0.0, 0.0, 0.1, 0.1, 0.1])
    assert out["n_nonzero"] == 3 and out["n_positive"] == 3


def test_sign_test_on_no_nonzero_drops_returns_no_p():
    out = sign_test([0.0, 0.0])
    assert out["n_nonzero"] == 0 and out["p"] is None


# --------------------------------------------------------------------------- control selection

def test_controls_are_selected_by_design_not_by_a_measured_statistic():
    neg, pos = _select_controls(PANELS)
    assert neg is not None and neg.design == "constructed_biparental"
    assert pos is not None and pos.design == "natural_population"


def test_the_selection_does_not_consult_pc1():
    """pc1 was MEASURED unfit (it ranked the NEGATIVE control first), so the selector must not
    read it. Pinned by source inspection: a regression would be silent otherwise."""
    src = (REPO / "scripts" / "f3_deconfound_calibration.py").read_text(encoding="utf-8")
    body = src.split("def _select_controls")[1].split("def ")[0]
    assert "pc1" not in body


# --------------------------------------------------------------------------- phenotype shapes

# `_phenotype_traits` is what decides WHICH VECTOR IS A TRAIT, and it does so by matching an
# axis length against the individual count. An orientation error there does not raise: it hands
# the calibration a transposed matrix whose "traits" are individuals, and every downstream number
# stays plausible. So each accepted shape is pinned by VALUE, not just by key count.

def test_a_trait_dict_keeps_only_vectors_whose_length_matches_the_rows():
    """A wrong-length vector is DROPPED, not padded or recycled -- a recycled vector would be
    scored against rows it never described."""
    out = _phenotype_traits({"a": [1.0, 2.0, 3.0], "b": [4.0, 5.0, 6.0], "short": [1.0, 2.0]},
                            3, 10)
    assert set(out) == {"a", "b"}
    assert out["a"].tolist() == [1.0, 2.0, 3.0]


def test_a_trait_dict_is_bounded_by_max_traits():
    out = _phenotype_traits({f"t{i}": [1.0, 2.0, 3.0] for i in range(6)}, 3, 2)
    assert len(out) == 2


def test_a_traits_by_rows_matrix_is_TRANSPOSED_so_a_trait_stays_a_trait():
    """The orientation branch, pinned by value. M is 3 traits x 4 individuals; each returned
    vector must be one ROW of M (that trait across individuals), length 4 -- not one column."""
    M = np.arange(12, dtype=float).reshape(3, 4)
    out = _phenotype_traits(M, 4, 10)
    assert len(out) == 3
    assert out["trait0"].tolist() == [0.0, 1.0, 2.0, 3.0]
    assert out["trait2"].tolist() == [8.0, 9.0, 10.0, 11.0]
    assert all(v.shape[0] == 4 for v in out.values())


def test_a_rows_by_traits_matrix_is_left_alone():
    """The same function must NOT transpose a matrix that is already (individuals, traits) --
    otherwise the fix for one orientation breaks the other."""
    M = np.arange(12, dtype=float).reshape(4, 3)
    out = _phenotype_traits(M, 4, 10)
    assert out["trait0"].tolist() == [0.0, 3.0, 6.0, 9.0]
    assert all(v.shape[0] == 4 for v in out.values())


def test_a_names_plus_matrix_tuple_keeps_the_real_trait_names():
    out = _phenotype_traits((["p", "q"], np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])), 3, 10)
    assert set(out) == {"p", "q"}
    assert out["q"].tolist() == [2.0, 4.0, 6.0]


def test_an_unusable_phenotype_shape_yields_no_traits_rather_than_a_guess():
    """A tuple with no matrix, and a bare 1-D vector, both return {}. The arm then scores
    nothing, which surfaces as INDETERMINATE_UNDERPOWERED rather than as a silent result."""
    assert _phenotype_traits((["a", "b"],), 3, 10) == {}
    assert _phenotype_traits(np.array([1.0, 2.0, 3.0]), 3, 10) == {}


# --------------------------------------------------------------------------- the artifact

def _artifact():
    hits = sorted((REPO / "wiki").glob("f3_deconfound_calibration_*.json"))
    if not hits:
        pytest.skip("calibration artifact not generated in this checkout")
    return json.loads(hits[-1].read_text(encoding="utf-8"))


def test_the_committed_artifact_parses_and_declares_a_known_verdict():
    a = _artifact()
    assert a["schema"] == "f3-deconfound-calibration-v1"
    assert a["verdict"] in VERDICTS


def test_a_blocked_artifact_emits_NO_metrics():
    """Pinned by key ABSENCE, not emptiness: a blocked run must not read as a partial result."""
    a = _artifact()
    if not a["verdict"].startswith("blocked:"):
        pytest.skip("this run is not blocked")
    assert "negative_control_result" not in a
    assert "positive_control_result" not in a
    assert "metrics_withheld" in a


def test_the_artifact_names_its_own_limits_including_the_grouping_alternative():
    lim = " ".join(_artifact()["honest_limits"]).lower()
    assert "two-sample" in lim                       # between-control comparison is not paired
    assert "by construction" in lim                  # the evidence asymmetry
    assert "one grouping" in lim                     # INVALID_INSENSITIVE's alternative cause
    assert str(N_PERM) in " ".join(_artifact()["honest_limits"])


def test_the_artifact_records_the_control_selection_rule_and_its_basis():
    sel = _artifact()["control_selection"]
    assert "BY CONSTRUCTION" in sel["rule"]
    assert "NOT by pc1" in sel["rule"]
    assert sel["negative_basis"] and sel["positive_basis"]


def test_the_asserted_constants_are_read_from_the_module_not_restated():
    """A shipped default is a promise: restating a constant in a test lets the two drift."""
    assert MIN_TRAITS_FOR_SIGN_TEST > 0
    assert MAX_TRAITS > 0
    assert N_PERM > 0
