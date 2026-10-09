"""Guards for the G-A encoder gate.

This gate can turn Step 3's PREDICTED negative into a MEASURED one, so the thing most worth pinning is
that it cannot do so dishonestly: a dirty null, a drifted protocol, or a failed reconcile must each REFUSE
to grade rather than return a band. All six branches are unit-tested on synthetic aggregates, with no
torch and no data.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path("scripts/glm_encoder_gate.py")


def _mod():
    spec = importlib.util.spec_from_file_location("eg", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _arms(gc, w20, multi, null, *, gc_seed0=None, protocol="registered"):
    """Synthetic per-arm aggregate in the shape `verdict` consumes."""
    eg = _mod()
    anchor = eg.PREREGISTERED["reconcile_anchor"]["expected"]
    s0 = anchor if gc_seed0 is None else gc_seed0
    reg = eg.REGISTERED_PROTOCOL.as_dict()
    p = reg if protocol == "registered" else {**reg, "lr": 0.5}

    def mk(med, seed0=None, proto=None):
        return {"median": med, "per_seed": [seed0 if seed0 is not None else med] + [med] * 9,
                "protocol": proto}

    return {
        "gc": mk(gc, s0, None),
        "encoder-w20": mk(w20, None, p),
        "encoder-multi": mk(multi, None, p),
        "encoder-shuffled": mk(null, None, p),
    }


def _seeds():
    return list(range(10))


# ---------------------------------------------------------------------------------------------------
# provenance and validity come BEFORE the margin
# ---------------------------------------------------------------------------------------------------
def test_a_seed_list_that_differs_from_the_registered_one_FAILS_CLOSED():
    """A median over a different seed set is not the registered gate, however good it looks. Give it a
    spectacular margin AND the wrong seeds -- it must still refuse."""
    eg = _mod()
    v, d, flags = eg.verdict(_arms(0.30, 0.90, 0.95, 0.00), seeds_used=[0, 1, 2])
    assert v == "INDETERMINATE_PROTOCOL_DRIFT"
    assert flags["registered_protocol_used"] is False
    assert "differ from the registered list" in d


def test_an_arm_running_under_a_DRIFTED_protocol_refuses_to_grade():
    """The artifact must not be able to claim the registered protocol was used when it was not -- that
    flag is what Step 9's MVP predicate reads."""
    eg = _mod()
    v, _d, flags = eg.verdict(_arms(0.30, 0.90, 0.95, 0.00, protocol="drifted"), seeds_used=_seeds())
    assert v == "INDETERMINATE_PROTOCOL_DRIFT"
    assert flags["registered_protocol_used"] is False


def test_a_MISSING_protocol_fails_closed_rather_than_attesting_from_absence():
    """The guard was `is not None and != registered`, which PASSES an arm reporting no protocol at all --
    so `registered_protocol_used=True` would be asserted from absence of evidence, and that is the exact
    flag Step 9's MVP predicate reads. Only the ridge comparator is legitimately protocol-free."""
    eg = _mod()
    arms = _arms(0.30, 0.90, 0.95, 0.00)
    arms["encoder-multi"]["protocol"] = None
    v, d, flags = eg.verdict(arms, seeds_used=_seeds())
    assert v == "INDETERMINATE_PROTOCOL_DRIFT"
    assert flags["registered_protocol_used"] is False
    assert "NO protocol" in d


def test_the_protocol_FREE_comparator_does_not_trip_the_guard():
    """Non-vacuity for the test above: gc is a ridge with no training protocol, and it must not be treated
    as drift -- otherwise the gate could never grade anything."""
    eg = _mod()
    arms = _arms(0.30, 0.33, 0.36, 0.00)
    assert arms["gc"]["protocol"] is None
    v, _d, flags = eg.verdict(arms, seeds_used=_seeds())
    assert v == "GENERALISES"
    assert flags["registered_protocol_used"] is True


def test_the_auroc_framing_records_its_own_DENOMINATOR():
    """A median over 3 measurable seeds must not sit beside a 10-seed primary with nothing saying so."""
    eg = _mod()
    runs = [{"spearman": 0.3, "auroc_active": 0.6, "protocol": None},
            {"spearman": 0.3, "auroc_active": None, "protocol": None},
            {"spearman": 0.3, "auroc_active": 0.7, "protocol": None}]
    a = eg._agg(runs)
    assert a["auroc_n_seeds_measurable"] == 2
    assert a["auroc_covers_all_seeds"] is False
    full = eg._agg([{"spearman": 0.3, "auroc_active": 0.6, "protocol": None}] * 3)
    assert full["auroc_covers_all_seeds"] is True


def test_a_failed_reconcile_blocks_EVERYTHING_even_with_a_huge_margin():
    """A script reading a different split than the published baseline cannot grade anything."""
    eg = _mod()
    v, d, flags = eg.verdict(_arms(0.30, 0.90, 0.95, 0.00, gc_seed0=0.11), seeds_used=_seeds())
    assert v == "INDETERMINATE_RECONCILE_FAILED"
    assert flags["reconcile_ok"] is False
    assert "ot reading the same split" in d  # sentence-initial capital in the message


def test_a_dirty_null_refuses_to_grade_rather_than_returning_a_band():
    """The validity gate. A shuffled-label arm that scores means the pipeline fabricates signal, so no
    band may be reported -- not even a cautious one."""
    eg = _mod()
    v, d, flags = eg.verdict(_arms(0.30, 0.40, 0.42, 0.20), seeds_used=_seeds())
    assert v == "INDETERMINATE_NULL_NOT_CLEAN"
    assert flags["null_clean"] is False
    assert "fabricates signal" in d
    assert "WEAK" not in v and "GENERALISES" not in v


def test_a_NEGATIVE_null_is_also_dirty_because_the_bar_is_on_MAGNITUDE():
    """A shuffled arm at -0.20 is just as much a pipeline artifact as +0.20; an unsigned bar would pass it."""
    eg = _mod()
    v, _d, _f = eg.verdict(_arms(0.30, 0.40, 0.42, -0.20), seeds_used=_seeds())
    assert v == "INDETERMINATE_NULL_NOT_CLEAN"


# ---------------------------------------------------------------------------------------------------
# the three substantive branches
# ---------------------------------------------------------------------------------------------------
def test_GENERALISES_needs_the_full_margin():
    eg = _mod()
    v, d, flags = eg.verdict(_arms(0.30, 0.33, 0.36, 0.00), seeds_used=_seeds())
    assert v == "GENERALISES", d
    assert flags["null_clean"] and flags["reconcile_ok"] and flags["registered_protocol_used"]


def test_a_small_positive_margin_is_WEAK_not_a_pass():
    """The distinction that keeps a 0.01 edge from being reported as the encoder earning its keep."""
    eg = _mod()
    v, _d, _f = eg.verdict(_arms(0.30, 0.31, 0.315, 0.00), seeds_used=_seeds())
    assert v == "WEAK"


def test_DOES_NOT_GENERALISE_when_the_encoder_fails_to_beat_one_GC_number():
    """The outcome Step 3 predicted. This is what upgrades its predicted-negative to a measured one."""
    eg = _mod()
    v, d, _f = eg.verdict(_arms(0.3000, 0.2926, 0.2880, 0.00), seeds_used=_seeds())
    assert v == "DOES_NOT_GENERALISE"
    assert "MARGIN" in d


def test_an_exact_tie_is_NOT_a_pass():
    """margin == 0 must fall to DOES_NOT_GENERALISE; a boundary that leaked upward would let a model that
    merely matched a single GC number be reported as beating it."""
    eg = _mod()
    v, _d, _f = eg.verdict(_arms(0.30, 0.30, 0.30, 0.00), seeds_used=_seeds())
    assert v == "DOES_NOT_GENERALISE"


def test_the_margin_is_taken_against_the_BEST_encoder_arm():
    """If the single-width arm wins, it must get the credit -- otherwise the multi-width arm's weakness
    would mask a real positive from a simpler model."""
    eg = _mod()
    v, d, _f = eg.verdict(_arms(0.30, 0.40, 0.31, 0.00), seeds_used=_seeds())
    assert v == "GENERALISES"
    assert "encoder-w20" in d


# ---------------------------------------------------------------------------------------------------
# the second framing must be able to measure something
# ---------------------------------------------------------------------------------------------------
def test_auroc_is_rank_based_and_handles_the_obvious_cases():
    eg = _mod()
    assert eg.auroc([0.1, 0.2, 0.3, 0.4], [0, 0, 1, 1]) == pytest.approx(1.0)
    assert eg.auroc([0.4, 0.3, 0.2, 0.1], [0, 0, 1, 1]) == pytest.approx(0.0)
    # [0,1,0,1] is NOT a coin flip -- the positives sit at ranks 2 and 4, i.e. better than chance.
    assert eg.auroc([0.1, 0.2, 0.3, 0.4], [0, 1, 0, 1]) == pytest.approx(0.75)
    # the genuine 0.5 case: one positive above, one below
    assert eg.auroc([0.1, 0.2, 0.3, 0.4], [0, 1, 1, 0]) == pytest.approx(0.5)


def test_auroc_with_ALL_TIED_scores_is_exactly_one_half_not_an_accident_of_sort_order():
    """Ties are handled by average ranks. Sort-order tie-breaking has silently shifted a median in this
    project before, so a fully-tied input is pinned explicitly."""
    eg = _mod()
    assert eg.auroc([0.5] * 6, [1, 0, 1, 0, 1, 0]) == pytest.approx(0.5)


def test_auroc_returns_None_when_a_CLASS_IS_ABSENT_rather_than_0_point_5():
    """0.5 would read as a measured coin-flip; None says the cell is unmeasurable. The difference matters
    because a seed whose held-out peaks happen to be all-inactive must not contribute a fake 0.5."""
    eg = _mod()
    assert eg.auroc([0.1, 0.2, 0.3], [0, 0, 0]) is None
    assert eg.auroc([0.1, 0.2, 0.3], [1, 1, 1]) is None


def test_the_incumbent_arm_ALSO_gets_an_auroc_so_the_second_framing_can_compare():
    """Scoring only the encoder arms on AUROC would leave the second framing unable to answer the question
    it exists for. `ridge_gc` must therefore return predictions, not just a correlation.

    This test previously ended in `... or True`, which can never fail -- a vacuous assertion of exactly
    the kind this repo keeps finding. It now checks the real contract: a two-element return whose second
    element is a usable prediction vector.
    """
    eg = _mod()

    class _T:
        def __init__(self, seq, expression):
            self.seq, self.expression = seq, expression

    import random
    r = random.Random(0)
    tr = [_T("".join(r.choice("ACGT") for _ in range(60)), r.random()) for _ in range(80)]
    te = [_T("".join(r.choice("ACGT") for _ in range(60)), r.random()) for _ in range(40)]
    out = eg.ridge_gc(tr, te)
    assert isinstance(out, tuple) and len(out) == 2, "ridge_gc must return (spearman, predictions)"
    sp, pred = out
    assert isinstance(sp, float)
    assert len(pred) == len(te), "one prediction per test row, or AUROC cannot be computed"
    assert all(isinstance(v, float) and v == v for v in pred)
    # and those predictions must actually drive an AUROC
    assert eg.auroc(pred, [1, 0] * 20) is not None


def test_fit_encoder_returns_predictions_so_the_gate_can_compute_a_second_framing():
    """Pins the seam. Without it the AUROC column would be silently all-None -- a check that measures
    nothing while appearing to."""
    from dna_decode.glm import encoder
    import inspect
    assert '"predictions": pred' in inspect.getsource(encoder.fit_encoder)


# ---------------------------------------------------------------------------------------------------
# the paired analysis -- the instrument the frozen median rule is NOT
# ---------------------------------------------------------------------------------------------------
def _ps(gc, other):
    return {"gc": {"per_seed": gc, "median": 0.0},
            "encoder-multi": {"per_seed": other, "median": 0.0}}


def test_a_positive_median_with_sign_flipping_margins_is_flagged_as_including_zero():
    """THE case the first real run produced: median margin +0.0122, but 6/10 wins with margins flipping
    sign, so the paired CI spans zero. Reporting the median alone would overstate it."""
    eg = _mod()
    gc = [0.3000, 0.3151, 0.3073, 0.2917, 0.3191, 0.3347, 0.3091, 0.3191, 0.3133, 0.2747]
    enc = [0.3235, 0.3526, 0.3293, 0.3233, 0.2872, 0.3126, 0.2896, 0.2995, 0.3312, 0.3600]
    pa = eg.paired_margins(_ps(gc, enc))["encoder-multi"]
    assert pa["wins"] == 6 and pa["n_seeds"] == 10
    assert pa["sign_test_p_two_sided"] > 0.5
    assert pa["ci95_includes_zero"] is True
    assert pa["mean_margin"] > 0, "the mean IS positive -- that is exactly why the CI matters"


def test_a_CONSISTENT_win_on_every_seed_does_NOT_include_zero():
    """Non-vacuity: the flag must be able to come out False, or it says nothing."""
    eg = _mod()
    gc = [0.30] * 10
    enc = [0.40, 0.41, 0.39, 0.42, 0.40, 0.38, 0.41, 0.40, 0.39, 0.41]
    pa = eg.paired_margins(_ps(gc, enc))["encoder-multi"]
    assert pa["wins"] == 10
    assert pa["sign_test_p_two_sided"] < 0.01
    assert pa["ci95_includes_zero"] is False


def test_the_margins_are_PAIRED_not_a_difference_of_aggregates():
    """If it subtracted medians instead of pairing per seed, two arms with identical medians but opposite
    per-seed ordering would look identical. They must not."""
    eg = _mod()
    a = eg.paired_margins(_ps([0.1, 0.2], [0.2, 0.1]))["encoder-multi"]
    b = eg.paired_margins(_ps([0.1, 0.2], [0.1, 0.2]))["encoder-multi"]
    assert a["per_seed_margin"] == [0.1, -0.1]
    assert b["per_seed_margin"] == [0.0, 0.0]
    assert a["sd_margin"] > b["sd_margin"]


def test_an_exact_tie_counts_as_a_NON_win_which_is_the_conservative_direction():
    eg = _mod()
    pa = eg.paired_margins(_ps([0.3] * 4, [0.3, 0.3, 0.4, 0.4]))["encoder-multi"]
    assert pa["wins"] == 2, "ties must not be credited as wins"


def test_the_comparator_is_excluded_from_its_own_paired_table():
    eg = _mod()
    assert "gc" not in eg.paired_margins(_ps([0.3] * 3, [0.3] * 3))


def test_the_paired_block_is_namespace_separate_from_the_verdict():
    """It must QUALIFY the frozen verdict, never silently alter it -- re-sizing a pre-registered bar after
    seeing the result is an authority call."""
    eg = _mod()
    src = SCRIPT.read_text(encoding="utf-8")
    assert '"paired_analysis": paired_margins(arms)' in src
    # the verdict function must not consult the paired analysis
    import inspect
    assert "paired" not in inspect.getsource(eg.verdict)


# ---------------------------------------------------------------------------------------------------
# pre-registration hygiene
# ---------------------------------------------------------------------------------------------------
def test_the_registered_seed_list_is_at_least_ten_seeds():
    eg = _mod()
    assert len(eg.PREREGISTERED["registered_seeds"]) >= 10


def test_every_verdict_the_function_can_return_is_documented():
    """A verdict string with no entry in the frozen table is an undocumented outcome."""
    eg = _mod()
    documented = set(eg.PREREGISTERED["verdicts"])
    for case in (_arms(0.30, 0.33, 0.36, 0.00), _arms(0.30, 0.305, 0.306, 0.00),
                 _arms(0.30, 0.29, 0.28, 0.00), _arms(0.30, 0.40, 0.42, 0.20),
                 _arms(0.30, 0.40, 0.42, 0.00, gc_seed0=0.11),
                 _arms(0.30, 0.40, 0.42, 0.00, protocol="drifted")):
        v, _d, _f = eg.verdict(case, seeds_used=_seeds())
        assert v in documented, v


def test_the_reconcile_anchor_points_at_a_committed_artifact_that_exists():
    eg = _mod()
    src = Path(eg.PREREGISTERED["reconcile_anchor"]["source"])
    assert src.exists(), f"reconcile anchor cites {src}, which does not exist"
    import json
    json.loads(src.read_text(encoding="utf-8"))


def test_the_gate_records_that_it_was_built_against_a_predicted_negative():
    """Provenance that matters for reading the result: the encoder was built on a user decision AGAINST
    its own stop gate, and the artifact must say so rather than presenting itself as neutral."""
    eg = _mod()
    src = SCRIPT.read_text(encoding="utf-8")
    assert "built_against_a_predicted_negative" in src
    assert "glm_tile_headroom_2026-10-08.json" in src
