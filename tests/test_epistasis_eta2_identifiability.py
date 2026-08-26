"""The pure math of the eta^2 identifiability probe + regressions on the claims it withdrew.

WHY THIS EXISTS. `wiki/forward_epistasis_h2_confirmed_2026-08-25.md` claimed H2 CONFIRMED by correlating the
additive score's pooling gain against eta^2 of FITNESS. A pooling gain needs BOTH the label and the predictor
to separate by group, and the score side was never measured. Measured, the two eta^2 are near-collinear and
the two well-powered proteins point in OPPOSITE directions when partialled -- so the contributions are not
separably identifiable, and the memo credited the wrong quantity.

The real-data half lives on the gitignored D: cache, so these tests pin (a) the pure helpers on synthetic
data where the answer is known by construction, and (b) that the withdrawn wording stays withdrawn.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from epistasis_eta2_identifiability import (  # noqa: E402
    eta2_of, fitness_groups, mean_score_by_order, partial_spearman, score_groups,
)

H2_MEMO = ROOT / "wiki" / "forward_epistasis_h2_confirmed_2026-08-25.md"
CORRECTION = ROOT / "wiki" / "forward_epistasis_pooling_correction_2026-08-25.md"
NEW_MEMO = ROOT / "wiki" / "forward_epistasis_eta2_identifiability_2026-08-25.md"


# ---------------------------------------------------------------------------- pure helpers

def test_eta2_is_one_when_groups_are_constant_and_differ():
    """All variance between groups, none within -> eta^2 == 1 exactly."""
    assert eta2_of([[1.0, 1.0, 1.0], [5.0, 5.0, 5.0]]) == pytest.approx(1.0)


def test_eta2_is_zero_when_groups_share_a_mean():
    """Identical groups -> no between-group variance at all."""
    assert eta2_of([[1.0, 3.0], [1.0, 3.0]]) == pytest.approx(0.0, abs=1e-12)


def test_score_and_fitness_groups_read_the_two_axes_apart():
    """The (fitness, score) pair ordering is load-bearing -- swapping them silently inverts the probe."""
    by_order = {2: [(0.1, -1.0), (0.2, -2.0)], 3: [(0.5, -9.0)]}
    assert fitness_groups(by_order, (2, 3)) == [[0.1, 0.2], [0.5]]
    assert score_groups(by_order, (2, 3)) == [[-1.0, -2.0], [-9.0]]


def test_partial_spearman_zeroes_out_a_fully_mediated_relationship():
    """If x~y is exactly what z~y and x~z imply, the partial is 0 -- the GFP case in miniature."""
    assert partial_spearman(r_xy=0.9 * 0.9, r_zy=0.9, r_xz=0.9) == pytest.approx(0.0, abs=1e-12)


def test_partial_spearman_returns_nan_under_perfect_collinearity_rather_than_a_fallback():
    """Under-determined is the HONEST answer. A fallback number here would fabricate identifiability."""
    assert partial_spearman(r_xy=0.9, r_zy=0.9, r_xz=1.0) != partial_spearman(r_xy=0.9, r_zy=0.9, r_xz=1.0)


def test_partial_spearman_can_go_negative_which_is_the_his7_finding():
    """HIS7's score-side partial is -0.135. A helper that clamped at 0 would have hidden the
    opposite-directions result that IS the evidence for non-identifiability."""
    assert partial_spearman(r_xy=0.80, r_zy=0.95, r_xz=0.90) < 0


def test_mean_score_by_order_exposes_the_sum_of_k_scaling():
    """The additive score is a SUM of k terms, so its mean must fall with k. That structural fact is the
    mechanism behind the whole finding, so the probe reports it directly rather than asserting it."""
    by_order = {2: [(0.0, -2.0), (0.0, -2.0)], 3: [(0.0, -3.0)], 4: [(0.0, -4.0)]}
    means = mean_score_by_order(by_order)
    assert means == {"2": -2.0, "3": -3.0, "4": -4.0}
    vals = [means[k] for k in sorted(means, key=int)]
    assert all(b < a for a, b in zip(vals, vals[1:])), "monotone decline in k is the claim under test"


# ------------------------------------------------------------------- withdrawn-claim regressions

_WITHDRAWAL_MARKERS = ("WITHDRAWN", "withdrawn", "AMENDED", "~~", "NOT inferential")


def _every_occurrence_is_retracted(text: str, phrase: str, window: int = 400) -> tuple[bool, int]:
    """Is every appearance of `phrase` accompanied by a withdrawal marker nearby?

    A bare `phrase not in text` CANNOT distinguish an assertion from a quotation inside its own retraction
    -- and a corrected memo necessarily quotes what it is correcting (struck through, or as "X is
    withdrawn"). That naive form failed on this very file. So the contract mirrors `_DISCLOSED_ABSENT` in
    `test_claude_md_citations.py`: the phrase MAY appear, but never unqualified.
    """
    i, seen = text.find(phrase), 0
    while i != -1:
        seen += 1
        around = text[max(0, i - window):i + len(phrase) + window]
        if not any(m in around for m in _WITHDRAWAL_MARKERS):
            return False, seen
        i = text.find(phrase, i + 1)
    return True, seen


def test_the_retraction_guard_is_not_vacuous():
    """A guard satisfied by MENTIONING the phrase is a mute button. Prove it still fails the bad case."""
    asserted = "The controls show the gain is **100%** between-order. Deployment rule follows."
    ok, seen = _every_occurrence_is_retracted(asserted, "the gain is **100%** between-order")
    assert not ok and seen == 1, "the guard would pass an un-retracted assertion"

    retracted = "~~the gain is **100%** between-order~~ - WITHDRAWN, the collapse is definitional."
    ok, seen = _every_occurrence_is_retracted(retracted, "the gain is **100%** between-order")
    assert ok and seen == 1

    # and one retracted occurrence must not launder a second, un-retracted one elsewhere
    mixed = retracted + ("\n\n" + "filler. " * 120) + asserted
    ok, seen = _every_occurrence_is_retracted(mixed, "the gain is **100%** between-order")
    assert not ok and seen == 2, "a nearby retraction must not cover a distant assertion"


@pytest.mark.skipif(not H2_MEMO.exists(), reason="H2 memo absent")
def test_the_h2_memo_no_longer_asserts_causality_or_a_100_percent_share():
    """Control A re-deals the SAME group sizes from one shuffled pool, so its collapse to 0.000 is
    definitional -- it checks the arithmetic, not a mechanism. 'causal' and '100%' were withdrawn."""
    text = H2_MEMO.read_text(encoding="utf-8", errors="replace")
    assert "AMENDED" in text, "the banner is the load-bearing part -- readers arrive via existing citations"
    for phrase in ("the gain is **100%** between-order",
                   "the gain is untouched"):
        ok, seen = _every_occurrence_is_retracted(text, phrase)
        assert ok, f"{phrase!r} appears asserted, not retracted"
        assert seen, f"{phrase!r} vanished entirely -- if intentional, drop it from this guard deliberately"


@pytest.mark.skipif(not H2_MEMO.exists(), reason="H2 memo absent")
def test_the_h2_p_values_are_marked_non_inferential():
    """26 nested subsets sharing most of their variants are not 26 independent observations; n_eff is the
    number of ORDERS. The p-values may stay on the page, but never unqualified."""
    text = H2_MEMO.read_text(encoding="utf-8", errors="replace")
    for p in ("7.4e-17", "9.3e-26"):
        i = text.find(p)
        assert i != -1, f"{p} vanished -- if the table was rewritten, update this guard deliberately"
        assert "NOT inferential" in text[i:i + 120], f"{p} is quoted without its non-inferential marker"


@pytest.mark.skipif(not CORRECTION.exists(), reason="correction memo absent")
def test_the_correction_memo_does_not_still_import_the_causal_h2_claim():
    """REGRESSION: the memo that was RIGHT was carrying the part that was wrong -- it cited H2 as
    'SINCE CONFIRMED ... two negative controls make it causal'."""
    text = CORRECTION.read_text(encoding="utf-8", errors="replace")
    ok, seen = _every_occurrence_is_retracted(text, "make it causal")
    assert ok, "the correction memo still asserts the causal H2 claim it imported"
    assert seen, "phrase vanished -- if the paragraph was rewritten, update this guard deliberately"
    assert "REFINED AND PARTLY SUPERSEDED" in text


@pytest.mark.skipif(not NEW_MEMO.exists(), reason="identifiability memo absent")
def test_the_new_memo_does_not_swap_in_the_score_side_as_the_next_overclaim():
    """The tempting correction -- 'eta^2 of the SCORE governs it' -- is refuted by HIS7's -0.135. The memo
    must land on the joint condition, not crown a new winner."""
    text = NEW_MEMO.read_text(encoding="utf-8", errors="replace")
    assert "not separably identifiable" in text
    assert "next overclaim" in text


@pytest.mark.skipif(not (ROOT / "dna_decode" / "forward" / "variant_effect.py").exists(), reason="absent")
def test_predict_multi_effect_docstring_no_longer_claims_the_joint_test_is_deferred():
    """STALE-DEFERRAL regression. The docstring said the joint epistasis test was 'deferred to a GPU run'
    while the sweep had scored joint on 5 proteins. Fifth instance of this class in this repo, and the one
    surface a docstring reader hits first."""
    text = (ROOT / "dna_decode" / "forward" / "variant_effect.py").read_text(encoding="utf-8",
                                                                            errors="replace")
    assert "deferred to a GPU run" not in text
    assert "HAS BEEN RUN" in text
    assert (ROOT / "wiki" / "forward_epistasis_sweep_2026-07-27.json").exists()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
