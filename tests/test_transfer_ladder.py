"""Guards for the transfer-ladder design — the metric swap and the derived ordering.

The plan this implements replaced AUROC with `coverage_lift` on a measurement, not a preference, so the
load-bearing tests are (a) that the base-rate-robustness CLAIM actually holds, and (b) that the ladder
ordering is DERIVED from sourced lineages rather than asserted from familiarity.
"""
from __future__ import annotations

import pytest

from dna_decode.eval import transfer_ladder as tl


# --------------------------------------------------------------------------------------------------
# coverage / coverage_lift
# --------------------------------------------------------------------------------------------------
def test_coverage_counts_nonzero_not_magnitude():
    assert tl.coverage([0, 0, 0, 0]) == 0.0
    assert tl.coverage([1, 2, 3, 4]) == 1.0
    assert tl.coverage([0, 5, 0, 5]) == 0.5
    # magnitude is irrelevant -- "did it fire" is the question
    assert tl.coverage([0.01, 0, 99.0, 0]) == 0.5


def test_coverage_refuses_an_empty_class():
    """An absent class is a missing MEASUREMENT, not 0.0 coverage. Returning 0.0 would let an empty
    label set produce a clean-looking lift on no data."""
    with pytest.raises(ValueError, match="empty"):
        tl.coverage([])


def test_coverage_lift_reproduces_both_measured_rungs():
    """The two real rungs, from the plan's Pre-Change Baseline. If these drift, the metric changed."""
    # E. coli: 351 essential at 0.3789 coverage, 3432 non-essential at 0.0664
    ess = [1] * 133 + [0] * (351 - 133)          # 133/351 = 0.3789
    non = [1] * 228 + [0] * (3432 - 228)         # 228/3432 = 0.0664
    assert round(tl.coverage_lift(ess, non), 4) == 0.3125
    # human: 115/681 = 0.1689, 5/899 = 0.0056
    ess_h = [1] * 115 + [0] * (681 - 115)
    non_h = [1] * 5 + [0] * (899 - 5)
    assert round(tl.coverage_lift(ess_h, non_h), 4) == 0.1633


def test_coverage_lift_is_BASE_RATE_ROBUST_the_whole_reason_it_replaced_auroc():
    """THE load-bearing claim. Each term is class-conditioned, so reweighting a class cannot move the
    statistic. AUROC does NOT have this property, which is why it is not comparable across rungs whose
    prevalence differs ~5x (E. coli 0.0928 genome-wide vs BAGEL's curated extremes 0.431)."""
    ess = [1, 1, 0, 0]          # coverage 0.5
    non = [1, 0, 0, 0]          # coverage 0.25
    base = tl.coverage_lift(ess, non)
    assert tl.coverage_lift(ess, non * 50) == base       # thin/fatten the negative class
    assert tl.coverage_lift(ess * 37, non) == base       # and the positive class
    assert tl.coverage_lift(ess * 7, non * 13) == base   # both at once


# --------------------------------------------------------------------------------------------------
# ordering, derived
# --------------------------------------------------------------------------------------------------
def test_shared_rank_depth_on_synthetic_lineages():
    assert tl.shared_rank_depth("a; b; c", "a; b; c") == 3
    assert tl.shared_rank_depth("a; b; c", "a; b; z") == 2
    assert tl.shared_rank_depth("a; b", "x; y") == 0
    assert tl.shared_rank_depth("a; b; c", "a; b") == 2      # prefix, shorter is the bound


def test_the_derived_order_puts_paeruginosa_NEARER_than_saureus():
    """The counter-intuitive case, and the reason ordering is derived rather than asserted: E. coli and
    P. aeruginosa are both Gammaproteobacteria (5 shared ranks) while S. aureus is a different phylum
    (2). A human ordering by familiarity would very likely get this backwards."""
    d = {r.key: tl.depth_from_reference(r) for r in tl.RUNGS}
    assert d["paeruginosa"] == 5, d
    assert d["saureus"] == 2, d
    assert d["paeruginosa"] > d["saureus"]


def test_the_two_eukaryotes_TIE_so_the_ladder_has_four_distance_classes():
    """Found by doing the derivation. Both share only 'cellular organisms' with E. coli, so they are
    ONE distance class -- which the design must report rather than invent an ordering for."""
    d = {r.key: tl.depth_from_reference(r) for r in tl.RUNGS}
    assert d["scerevisiae"] == 1 and d["human"] == 1, d
    classes = tl.distance_classes()
    assert len(classes) == 4, classes
    tied = [keys for depth, keys in classes if len(keys) > 1]
    assert tied == [["human", "scerevisiae"]], classes


def test_every_rung_carries_a_sourced_lineage_and_a_url():
    for r in tl.RUNGS:
        assert r.lineage.startswith("cellular organisms"), r.key
        assert r.source_url.startswith("https://eutils.ncbi.nlm.nih.gov/"), r.key
        assert str(r.taxid) in r.source_url, r.key


def test_a_rung_without_provenance_RAISES_rather_than_defaulting():
    """The ecoff_catalog gate, applied here: an unsourced lineage would let a remembered taxonomy decide
    the ladder's ordered axis, and a wrong ordering produces a fully self-consistent wrong verdict."""
    ok = dict(key="x", taxid=1, scientific_name="X", label_source_id="s",
              technology="transposon_insertion", lineage="cellular organisms; Bacteria",
              source_url="https://eutils.ncbi.nlm.nih.gov/x?id=1", annotation_kind="k", key_space="ks")
    tl.Rung(**ok)                                          # control: valid
    for bad in ({"lineage": ""}, {"lineage": "   "}, {"source_url": ""}, {"source_url": "nope"},
                {"technology": "vibes"}):
        with pytest.raises(tl.RungProvenanceError):
            tl.Rung(**{**ok, **bad})


# --------------------------------------------------------------------------------------------------
# verdict
# --------------------------------------------------------------------------------------------------
def _r(key, depth, lift, null=0.02, adj=None, scored=True):
    d = {"key": key, "depth": depth, "coverage_lift": lift, "null_p95": null, "scored": scored}
    if adj is not None:
        d["coverage_lift_adjusted"] = adj
    return d


def test_powering_is_checked_BEFORE_any_decay_arithmetic():
    """A run with too few rungs necessarily has no decay to find, so checking decay first would publish
    an unpowered run as a finding."""
    v = tl.classify_ladder([_r("ecoli", 8, 0.31), _r("human", 1, 0.16)],
                           cliff_drop=0.10, plateau_tol=0.05, min_rungs=3)
    assert v.verdict == tl.INDETERMINATE_INSUFFICIENT_RUNGS
    assert v.n_scored == 2


def test_an_unscored_rung_does_not_count_toward_powering():
    rows = [_r("ecoli", 8, 0.31), _r("human", 1, 0.16), _r("saureus", 2, 0.0, scored=False)]
    v = tl.classify_ladder(rows, cliff_drop=0.10, plateau_tol=0.05, min_rungs=3)
    assert v.verdict == tl.INDETERMINATE_INSUFFICIENT_RUNGS


def test_a_lift_inside_its_own_null_blocks_a_verdict():
    rows = [_r("ecoli", 8, 0.31), _r("paeruginosa", 5, 0.29), _r("human", 1, 0.01, null=0.05)]
    v = tl.classify_ladder(rows, cliff_drop=0.10, plateau_tol=0.05, min_rungs=3)
    assert v.verdict == tl.INDETERMINATE_NULL_NOT_CLEARED
    assert v.detail["not_cleared"] == ["human"]


def test_bar_sensitivity_refuses_a_verdict_the_vocabulary_lesson_encoded():
    """If subtracting the phrasing floor moves a rung into a different band, the instrument is not
    deciding. Encodes the retraction recorded in wiki/essentiality_missed_vocabulary_2026-10-05."""
    rows = [_r("ecoli", 8, 0.31), _r("paeruginosa", 5, 0.30),
            _r("human", 1, 0.16, adj=0.29)]      # raw says low, adjusted says high
    v = tl.classify_ladder(rows, cliff_drop=0.10, plateau_tol=0.05, min_rungs=3)
    assert v.verdict == tl.INDETERMINATE_BAR_SENSITIVE
    assert v.detail["bar_sensitive"] == ["human"]


def test_no_decay_when_every_rung_sits_together():
    rows = [_r("ecoli", 8, 0.31), _r("paeruginosa", 5, 0.30), _r("saureus", 2, 0.29)]
    v = tl.classify_ladder(rows, cliff_drop=0.10, plateau_tol=0.05, min_rungs=3)
    assert v.verdict == tl.NO_DECAY_DISTANCE_INSENSITIVE


def test_domain_cliff_when_bacteria_hold_and_eukaryotes_fall_away():
    rows = [_r("ecoli", 8, 0.31), _r("paeruginosa", 5, 0.30), _r("saureus", 2, 0.29),
            _r("scerevisiae", 1, 0.15), _r("human", 1, 0.16)]
    v = tl.classify_ladder(rows, cliff_drop=0.10, plateau_tol=0.05, min_rungs=3)
    assert v.verdict == tl.DOMAIN_CLIFF


def test_smooth_decay_when_the_drop_is_graded():
    rows = [_r("ecoli", 8, 0.31), _r("paeruginosa", 5, 0.24), _r("saureus", 2, 0.17),
            _r("human", 1, 0.10)]
    v = tl.classify_ladder(rows, cliff_drop=0.10, plateau_tol=0.02, min_rungs=3)
    assert v.verdict == tl.SMOOTH_DECAY_WITH_DISTANCE


def test_classify_ladder_takes_no_default_thresholds():
    """Thresholds are the runner's to freeze AFTER the distribution analysis. A default here would be a
    bar frozen at design time, which is the failure this repo has recorded twice."""
    with pytest.raises(TypeError):
        tl.classify_ladder([_r("ecoli", 8, 0.31)])          # type: ignore[call-arg]


def test_preregistered_names_every_verdict_classify_ladder_can_return():
    """A bar silent on an outcome cannot adjudicate it."""
    assert set(tl.PREREGISTERED["verdicts"]) == set(tl.VERDICTS)
    assert tl.PREREGISTERED["primary_test"] == tl.PRIMARY_TEST
    assert len(tl.VERDICTS) == 6


def test_preregistered_records_why_auroc_was_demoted_and_the_eukaryote_tie():
    """Both are findings a later reader would otherwise have to re-derive."""
    assert "base-rate robust" in tl.PREREGISTERED["primary_metric"]
    assert "tie-mass" in tl.PREREGISTERED["secondary_metric"]
    assert "equal lineage depth" in tl.PREREGISTERED["eukaryote_tie"]
