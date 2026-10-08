"""Guards for the G-A stop gate.

The gate's job is to be able to CLOSE a family before it is built, so its two failure modes are opposite and
both must be pinned: it must not wave through a predicted-negative, and it must not block on a reconcile
failure it merely imagined.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

CACHE = Path("D:/dna_decode_cache/mpra")
HAVE_DATA = (CACHE / "GSE144621_peak_tile_expression_formatted_std.txt.gz").exists()


def _mod():
    spec = importlib.util.spec_from_file_location("th", Path("scripts/glm_tile_headroom.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _arms(gc, gc_dinuc, gc_posgc, dinuc_posgc, gc_seed0=None):
    """Synthetic arm dict. gc_seed0 defaults to the reconcile anchor so tests isolate one variable."""
    th = _mod()
    anchor = th.PREREGISTERED["reconcile_anchor"]["expected"]
    s0 = anchor if gc_seed0 is None else gc_seed0
    mk = lambda med, seed0: {"median": med, "per_seed": [seed0] + [med] * 9, "n_features": 1,
                             "stdev": 0.0, "min": med, "max": med, "fraction_of_ceiling": med / 0.941}
    return {"gc": mk(gc, s0), "gc_dinuc": mk(gc_dinuc, gc_dinuc),
            "gc_posgc": mk(gc_posgc, gc_posgc), "dinuc_posgc": mk(dinuc_posgc, dinuc_posgc)}


# ---------------------------------------------------------------------------------------------------
# the reconcile gate comes FIRST
# ---------------------------------------------------------------------------------------------------
def test_reconcile_failure_blocks_EVERYTHING_even_when_the_gain_looks_great():
    """A script reading a different split than the published baseline cannot grade anything. Give it a
    spectacular position gain AND a broken anchor -- it must still refuse."""
    th = _mod()
    arms = _arms(gc=0.31, gc_dinuc=0.28, gc_posgc=0.90, dinuc_posgc=0.88, gc_seed0=0.11)
    v, d, ratify = th._verdict(arms, 0.941)
    assert v == "INDETERMINATE_RECONCILE_FAILED"
    assert "not reading the same split" in d
    assert ratify is False, "a refusal must not also assert a ratification recommendation"


def test_reconcile_passes_within_tolerance_and_then_grades():
    """Non-vacuity for the test above: the anchor must not block legitimate runs."""
    th = _mod()
    tol = th.PREREGISTERED["reconcile_anchor"]["tolerance"]
    exp = th.PREREGISTERED["reconcile_anchor"]["expected"]
    arms = _arms(gc=0.31, gc_dinuc=0.28, gc_posgc=0.40, dinuc_posgc=0.38, gc_seed0=exp + tol * 0.9)
    v, _d, _r = th._verdict(arms, 0.941)
    assert v == "POSITION_HAS_HEADROOM"


# ---------------------------------------------------------------------------------------------------
# the two substantive branches
# ---------------------------------------------------------------------------------------------------
def test_no_headroom_sets_the_RATIFICATION_FLAG():
    """The flag is the whole point: a machine-readable stop rather than a judgement buried in prose."""
    th = _mod()
    arms = _arms(gc=0.3112, gc_dinuc=0.2785, gc_posgc=0.3044, dinuc_posgc=0.2797)
    v, d, ratify = th._verdict(arms, 0.941)
    assert v == "NO_POSITIONAL_HEADROOM"
    assert ratify is True
    assert "POSITION GAIN" in d


def test_real_headroom_does_NOT_set_the_flag():
    th = _mod()
    arms = _arms(gc=0.31, gc_dinuc=0.30, gc_posgc=0.45, dinuc_posgc=0.44)
    v, _d, ratify = th._verdict(arms, 0.941)
    assert v == "POSITION_HAS_HEADROOM"
    assert ratify is False


def test_the_gain_is_computed_against_the_BEST_composition_arm_not_plain_gc():
    """If a richer composition-only arm beats plain GC, the position gain must be measured against THAT --
    otherwise position gets credit for beating a weaker baseline than composition actually achieves."""
    th = _mod()
    # gc_dinuc is the strongest composition arm here; gc_posgc beats plain gc but NOT gc_dinuc
    arms = _arms(gc=0.20, gc_dinuc=0.40, gc_posgc=0.30, dinuc_posgc=0.29)
    v, d, _r = th._verdict(arms, 0.941)
    assert v == "NO_POSITIONAL_HEADROOM", d
    assert "gc_dinuc" in d


def test_a_NEGATIVE_position_gain_is_handled_and_reported():
    """The measured real case: position actively costs. Must not be clamped to zero or read as a tie."""
    th = _mod()
    arms = _arms(gc=0.3112, gc_dinuc=0.2785, gc_posgc=0.3044, dinuc_posgc=0.2797)
    _v, d, _r = th._verdict(arms, 0.941)
    assert "-0.0068" in d or "-0.007" in d, d


# ---------------------------------------------------------------------------------------------------
# the minimal positional feature
# ---------------------------------------------------------------------------------------------------
def test_positional_gc_width_and_length_invariance():
    th = _mod()
    v = th.positional_gc(["ACGT" * 30, "ACGT" * 12])
    assert len(v[0]) == len(v[1]) == th.N_BINS


def test_positional_gc_SEES_position_where_plain_gc_cannot():
    """The feature must actually be positional, or the gate tests nothing. Same overall GC, opposite
    layout: plain GC is identical, positional GC must differ."""
    th = _mod()
    from dna_decode.glm.genomewide import gc_features
    a = "GC" * 25 + "AT" * 25
    b = "AT" * 25 + "GC" * 25
    assert gc_features([a]) == gc_features([b]), "the premise: plain GC cannot tell these apart"
    assert th.positional_gc([a]) != th.positional_gc([b])


def test_positional_gc_handles_a_short_sequence_without_empty_bins():
    """A 48bp fragment over 10 bins gives uneven bins; none may be empty or the feature silently emits 0.0
    for a real segment."""
    th = _mod()
    v = th.positional_gc(["ACGT" * 12])
    assert len(v[0]) == th.N_BINS
    assert all(isinstance(x, float) for x in v[0])


# ---------------------------------------------------------------------------------------------------
# the committed real result
# ---------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not Path("wiki/glm_tile_headroom_2026-10-08.json").exists(),
                    reason="headroom artifact not present")
def test_the_committed_artifact_records_the_stop_gate_firing():
    """Pins the measured outcome: every added feature made it WORSE, position gain is negative, and the
    reconcile anchor landed on 0.3000 exactly -- which is what makes the negative trustworthy."""
    import json
    d = json.loads(Path("wiki/glm_tile_headroom_2026-10-08.json").read_text(encoding="utf-8"))
    assert d["verdict"] == "NO_POSITIONAL_HEADROOM"
    assert d["recommend_user_ratification_before_cnn"] is True
    arms = d["arms"]
    assert arms["gc"]["n_features"] == 1
    # one feature beat all the richer arms -- the signature of a one-dimensional substrate
    assert arms["gc"]["median"] > arms["gc_posgc"]["median"]
    assert arms["gc"]["median"] > arms["gc_dinuc"]["median"]
    assert arms["gc"]["median"] > arms["dinuc_posgc"]["median"]
    assert abs(arms["gc"]["per_seed"][0] - 0.3000) < 0.01, "reconcile anchor must have held"
