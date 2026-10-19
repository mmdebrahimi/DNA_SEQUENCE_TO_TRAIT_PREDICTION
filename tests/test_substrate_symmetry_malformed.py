"""Malformed-census shapes: what the validator does when the census is not the shape it expects.

Separate from `tests/test_substrate_symmetry.py`, which pins the intended ladder on well-formed
input. The question here is the FAIL DIRECTION, because this module has exactly one verdict it must
never produce by accident -- `exhausted` -- and exactly one reading it must never make available --
a rung silently absent, which reads as "nothing to say here".

Two of these pin KNOWN RESIDUALS rather than desired behaviour. They are labelled, they assert the
current shape exactly, and they say what a fix would look like, so a future fix fails them loudly
and the author updates the docstring deliberately instead of a silent drift.
"""
from __future__ import annotations

import pytest

from dna_decode.eval.substrate_symmetry import (
    SYMMETRY_REQUIRED,
    audit,
    ceiling_guard,
    check_symmetry,
    utilization,
)

FULL = {f: 1 for f in SYMMETRY_REQUIRED}


def _c(rungs) -> dict:
    return {"schema": "test", "rungs": rungs}


# ==================================================================================================
# the `rungs` container is not a mapping
# ==================================================================================================
@pytest.mark.parametrize("rungs", [
    [{"peak_tiles": {}}],          # a LIST of rungs -- a plausible schema-v2 shape
    "peak_tiles",                  # a bare string
    None,                          # explicit null
    42,
])
def test_a_non_mapping_rungs_container_is_indeterminate_never_ok(rungs):
    """`_rungs`'s isinstance guard is the only thing standing between a schema change and a vacuous
    green: without it, `.items()` would raise and -- in the renderer, which swallows everything --
    print UNAVAILABLE, while a direct caller of `check_symmetry` would crash. Returning
    `indeterminate` means no caller can mistake a shape change for a passing census."""
    r = check_symmetry(_c(rungs))
    assert r.verdict == "indeterminate:no-rungs"
    assert not r.ok
    assert r.n_rungs == 0
    assert audit(_c(rungs))["utilization"] == []


# ==================================================================================================
# row counts that are present but unusable
# ==================================================================================================
@pytest.mark.parametrize("bad", [0, -5, -1])
def test_zero_or_negative_row_count_is_unknown_never_exhausted(bad):
    """A nonsense row count must not become a denominator. `used/0` would raise and a negative
    denominator would invert the fraction, so either could manufacture `exhausted` from garbage."""
    u = utilization("x", {"n": bad, "largest_n_measured": 30})
    assert u.verdict == "unknown:no row count"
    assert u.exhausted is False
    assert u.fraction is None


def test_a_zero_in_a_higher_precedence_key_FALLS_THROUGH_instead_of_poisoning():
    """A post-filter count of 0 beside a real raw count must not read as "this rung has no data" --
    that is the module's own target error class (a missing/unusable key reading as absent data). The
    precedence loop keeps scanning rather than returning on the first key it finds."""
    u = utilization("peak_tiles", {"rows_after_published_filters": 0, "rows_in_file": 46713,
                                   "largest_n_measured": 30000})
    assert u.rows_available == 46713, "the loop returned on the first key instead of scanning on"
    assert not u.verdict.startswith("unknown:"), "a usable denominator exists; refusing it is the bug"
    assert u.verdict == "partial:64.2%"


@pytest.mark.parametrize("bad", [0, -5, 30000.0, "30000", None, [30000]])
def test_a_non_positive_int_largest_n_measured_is_unknown_never_exhausted(bad):
    """Includes the float: JSON `30000.0` parses to a Python float and `isinstance(x, int)` is False,
    so it degrades to `unknown` -- the SAFE direction. Pinned because the safe direction here is a
    consequence of a strict isinstance check that a future "be lenient about numbers" edit would
    remove, and leniency would have to be applied to the row count too, where `True` is already a
    live hazard (below)."""
    u = utilization("x", {"n": 100, "largest_n_measured": bad})
    assert u.verdict.startswith("unknown:"), u.verdict
    assert u.exhausted is False


# ==================================================================================================
# FORMERLY KNOWN RESIDUALS -- all three FIXED 2026-10-13, these now pin the fix
# ==================================================================================================
# Each of these was written as a characterization test asserting defect-shaped CURRENT behaviour, with
# a docstring naming the fix and instructing the fixer to rewrite it. All three failed toward
# `exhausted` or toward silence -- the two outcomes this module exists to prevent -- so leaving them
# would have made the module's own documented fail-safe an over-claim. Fixed, and these now assert the
# fixed shape. The defect each one describes is kept in the docstring so the fix cannot be casually
# reverted without reading what it cost.
def test_a_non_dict_rung_is_REPORTED_not_silently_dropped():
    """WAS: a rung whose value is not a dict was skipped by every `isinstance(r, dict)` filter, so
    symmetry returned `ok` and `audit()`'s table omitted the rung entirely -- this module's own stated
    failure class ("a rung absent reads as nothing to say here") reproduced inside the validator built
    to prevent it. `n_rungs` counted it while nothing scored it.

    NOW: a `malformed:<rung>` verdict ahead of the asymmetry ladder, and the rung gets a ROW in the
    table rather than being filtered out of it."""
    census = _c({"a": dict(FULL), "b": dict(FULL), "malformed": "not a dict"})
    rep = check_symmetry(census)
    assert rep.verdict == "malformed:malformed", rep.detail
    assert not rep.ok
    assert "only the 2 scoreable of 3" in rep.detail
    scored = [u["rung"] for u in audit(census)["utilization"]]
    assert scored == ["a", "b", "malformed"], "the malformed rung must APPEAR, flagged"
    bad = [u for u in audit(census)["utilization"] if u["rung"] == "malformed"][0]
    assert bad["verdict"].startswith("malformed:not-a-mapping")


def test_largest_n_measured_ABOVE_rows_available_is_SURFACED():
    """WAS: a clean `exhausted`, with an impossible `fraction > 1.0` as the only trace. `exhausted` is
    the one verdict that licenses "stop looking for more data", so hiding a census arithmetic error
    behind it is the costly direction. The most plausible cause is a stale `largest_n_measured` left
    when a row count is corrected DOWNWARD -- the edit this census just underwent.

    NOW: an explicit `inconsistent:` verdict, and it is NOT `exhausted`."""
    u = utilization("x", {"n": 100, "largest_n_measured": 500})
    assert u.verdict.startswith("inconsistent:"), u.verdict
    assert u.exhausted is False, "an impossible count must never license 'stop looking'"
    assert u.fraction == 5.0, "the impossible ratio is still reported, not hidden"


def test_a_bool_in_a_row_count_key_is_REFUSED():
    """WAS the worst fail direction available: `bool` is an `int` subclass, so `n: true` passed
    `isinstance(v, int) and v > 0`, `rows_available` became 1, and ANY measured N produced
    `exhausted` -- the module's forbidden claim, out of a type error. The live census is full of
    booleans (`clears_project_bar`, `trained_on_this_substrate`) and `n` is a generic enough key to
    collide.

    NOW: `type(v) is int` excludes bool without excluding any real count."""
    u = utilization("x", {"n": True, "largest_n_measured": 50})
    assert u.rows_available is None, "a bool must not be read as a row count"
    assert u.verdict.startswith("unknown:"), u.verdict
    assert u.exhausted is False
    # and a bool in largest_n_measured is equally refused
    assert utilization("x", {"n": 100, "largest_n_measured": True}).verdict.startswith("unknown:")


def test_a_real_int_row_count_still_works_after_the_bool_fix():
    """Non-vacuity on the fix itself: `type(v) is int` must not have excluded ordinary counts."""
    u = utilization("x", {"n": 1000, "largest_n_measured": 1000})
    assert u.rows_available == 1000
    assert u.verdict == "exhausted"


# ==================================================================================================
# RESTORED 2026-10-13. These two were lost when the residuals section above was rewritten by
# truncating the file at a marker without checking what followed it -- a 3-collected-case deletion
# that `6496 -> 6493` on the full suite is what surfaced. Recorded rather than quietly re-added,
# because "I only read to line 142 and assumed the file ended" is the same not-looking failure the
# module under test exists to catch.
def test_an_explicit_null_counts_as_PRESENT_not_absent():
    """LOAD-BEARING: it is WHY the live census passes. `sigma70_designed` declares
    `ceiling_for_a_sequence_model: null` together with a `ceiling_reason`, and that must read as
    "we know this is unknowable and said so" -- not as the silent omission the validator refuses.

    The detectable defect is the MISSING KEY, never the unknown VALUE. If this inverted, the honest
    way to record an unknowable ceiling would itself become a violation, and the only way to pass
    would be to invent a number."""
    declared = _c({"a": dict(FULL), "b": {**FULL, "ceiling_for_a_sequence_model": None}})
    assert check_symmetry(declared).verdict == "ok", "a declared null must be PRESENT"
    omitted = dict(FULL)
    del omitted["ceiling_for_a_sequence_model"]
    rep = check_symmetry(_c({"a": dict(FULL), "b": omitted}))
    assert rep.verdict == "asymmetric:ceiling_for_a_sequence_model", rep.detail


def test_every_asymmetric_field_is_in_the_dict_though_the_verdict_names_one():
    """The verdict string carries only the alphabetically-first offender, so a consumer reading the
    verdict alone would fix one field and believe it was done. `asymmetric_fields` carries them all,
    and the live census had TWO (`size_lever` AND `best_measured_arm`, both 1-of-4) -- so a one-field
    reading would have left half the defect in place."""
    stripped = {k: v for k, v in FULL.items() if k not in ("size_lever", "best_measured_arm")}
    rep = check_symmetry(_c({"a": dict(FULL), "b": stripped}))
    assert rep.verdict == "asymmetric:best_measured_arm", "alphabetically first"
    assert set(rep.asymmetric_fields) == {"size_lever", "best_measured_arm"}, rep.asymmetric_fields
    assert rep.asymmetric_fields["size_lever"] == ["b"]
