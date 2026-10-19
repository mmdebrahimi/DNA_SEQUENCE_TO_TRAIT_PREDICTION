"""Tests for the substrate-symmetry validator.

The load-bearing one is `test_refuses_the_REAL_live_census`: a validator that passes the file we
already know is asymmetric has not been shown to detect anything. It skips (never silently passes)
when the census is absent, so a fresh clone cannot turn it into a vacuous green.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_decode.eval.substrate_symmetry import (
    EXHAUSTED_AT,
    ROW_COUNT_KEYS,
    SYMMETRY_REQUIRED,
    UNDER_USED_BELOW,
    audit,
    ceiling_guard,
    check_symmetry,
    utilization,
)

CENSUS = Path(__file__).resolve().parent.parent / "wiki" / "mpra_substrate_census_2026-10-12.json"


def _census(rungs: dict) -> dict:
    return {"schema": "test", "rungs": rungs}


# --------------------------------------------------------------------------------------------------
# the ladder -- every rung reachable
# --------------------------------------------------------------------------------------------------
def test_no_rungs_is_indeterminate_not_ok():
    r = check_symmetry({})
    assert r.verdict == "indeterminate:no-rungs"
    assert not r.ok


def test_single_rung_REFUSES_rather_than_passing():
    """Symmetry is undefined over one rung. Returning ok would mean the check silently starts failing
    when a second rung is added, for a reason the author cannot connect to their edit."""
    r = check_symmetry(_census({"only": {f: 1 for f in SYMMETRY_REQUIRED}}))
    assert r.verdict == "indeterminate:single-rung"
    assert not r.ok


def test_symmetric_census_is_ok():
    rungs = {n: {f: 1 for f in SYMMETRY_REQUIRED} for n in ("a", "b", "c")}
    r = check_symmetry(_census(rungs))
    assert r.verdict == "ok", r.detail
    assert r.ok
    assert r.n_rungs == 3


def test_field_on_one_of_three_rungs_is_ASYMMETRIC_and_names_the_missing():
    """Non-vacuity: the detector must actually fire on the shape it exists for, and say which rungs."""
    rungs = {
        "has_it": {f: 1 for f in SYMMETRY_REQUIRED},
        "missing_a": {f: 1 for f in SYMMETRY_REQUIRED if f != "size_lever"},
        "missing_b": {f: 1 for f in SYMMETRY_REQUIRED if f != "size_lever"},
    }
    r = check_symmetry(_census(rungs))
    assert r.verdict == "asymmetric:size_lever", r.detail
    assert r.asymmetric_fields["size_lever"] == ["missing_a", "missing_b"]


def test_field_absent_from_EVERY_rung_is_not_a_violation():
    """Uniform absence is a real choice (nothing to say on any rung). Only a SPLIT is the defect."""
    rungs = {n: {f: 1 for f in SYMMETRY_REQUIRED if f != "best_measured_arm"} for n in ("a", "b")}
    assert check_symmetry(_census(rungs)).verdict == "ok"


# --------------------------------------------------------------------------------------------------
# utilization -- and the fail-safe direction
# --------------------------------------------------------------------------------------------------
def test_under_used_reports_the_real_fragment_shape():
    u = utilization("fragments_LB", {"n": 321123, "largest_n_measured": 30000})
    assert u.verdict.startswith("under-used:"), u.verdict
    assert "9.3" in u.verdict
    assert u.fraction == pytest.approx(0.0934, abs=1e-3)


def test_exhausted_reports_the_real_tile_shape():
    u = utilization("peak_tiles", {"rows_after_published_filters": 44106,
                                   "rows_in_file": 46713, "largest_n_measured": 44106})
    assert u.verdict == "exhausted"
    assert u.exhausted


def test_unknown_NEVER_counts_as_exhausted():
    """The one failure direction that matters: an unmeasured rung must not read as spent data."""
    for rung in ({"n": 321123}, {"largest_n_measured": 30000}, {}):
        u = utilization("x", rung)
        assert u.verdict.startswith("unknown:"), u.verdict
        assert u.exhausted is False


def test_partial_sits_between_the_two_bands():
    u = utilization("x", {"n": 1000, "largest_n_measured": 700})
    assert u.verdict.startswith("partial:")
    assert UNDER_USED_BELOW <= u.fraction < EXHAUSTED_AT


@pytest.mark.parametrize("key", list(ROW_COUNT_KEYS))
def test_every_declared_row_count_key_is_actually_read(key):
    """A key missing from ROW_COUNT_KEYS does not read as "the list is incomplete" -- it reads as
    "that rung has no data", which is a silent `unknown`.

    MEASURED 2026-10-13: the first version of this module listed only the tile and fragment keys, so
    `sigma70_designed` reported `unknown:no row count` while recording 10,898 rows under THREE names
    (`designed_rows`, `usable_joined`, `n_distinct_sequences`). The error class this whole module
    targets -- claiming data is absent without reading the fields -- was reproduced inside the fix
    for it. Parametrising over the real tuple means a key that is declared but not honoured fails."""
    u = utilization("x", {key: 500, "largest_n_measured": 500})
    assert u.rows_available == 500, f"{key} is declared in ROW_COUNT_KEYS but not read"
    assert u.verdict == "exhausted"


def test_row_count_precedence_prefers_post_filter_over_raw():
    """A model sees the POST-filter rows, so utilization must be measured against those. Reading the
    raw file count instead would understate utilization and could turn `exhausted` into `partial`."""
    u = utilization("peak_tiles", {"rows_in_file": 46713, "rows_after_published_filters": 44106,
                                   "largest_n_measured": 44106})
    assert u.rows_available == 44106
    assert u.verdict == "exhausted"


# --------------------------------------------------------------------------------------------------
# the ceiling guard -- enforcing in code a rule the census only stated in prose
# --------------------------------------------------------------------------------------------------
def test_ceiling_guard_fires_on_the_exact_number_that_was_misquoted():
    w = ceiling_guard("fragments_LB", {"single_replicate_agreement": 0.4583,
                                       "ceiling_for_a_sequence_model": 0.7928})
    assert len(w) == 1
    assert "0.4583" in w[0] and "0.7928" in w[0]


def test_ceiling_guard_silent_when_only_one_value_exists():
    assert ceiling_guard("x", {"ceiling_for_a_sequence_model": 0.941}) == []
    assert ceiling_guard("x", {"single_replicate_agreement": 0.79}) == []


def test_ceiling_guard_silent_when_the_two_coincide():
    """If a rung genuinely has one replicate, the step-up is identity and there is nothing to warn."""
    assert ceiling_guard("x", {"single_replicate_agreement": 0.5,
                               "ceiling_for_a_sequence_model": 0.5}) == []


# --------------------------------------------------------------------------------------------------
# the rule cannot be silently disabled
# --------------------------------------------------------------------------------------------------
def test_symmetry_required_cannot_be_emptied():
    """An empty required-set makes check_symmetry return ok on everything -- a green that means
    nothing. Pinning non-emptiness makes deleting the rule a visible test failure."""
    assert len(SYMMETRY_REQUIRED) >= 3
    assert "size_lever" in SYMMETRY_REQUIRED
    rungs = {n: {} for n in ("a", "b")}
    assert check_symmetry(_census(rungs), required=()).verdict == "ok"  # documents the degenerate case


# --------------------------------------------------------------------------------------------------
# THE REAL FILE
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not CENSUS.exists(), reason="census artifact absent (gitignored / fresh clone)")
def test_the_live_census_is_symmetric_and_fully_populated():
    """After the Step-4 amendment every rung carries every required field. This is the regression
    guard that stops a FUTURE rung being added silently asymmetric."""
    census = json.loads(CENSUS.read_text(encoding="utf-8"))
    rep = check_symmetry(census)
    assert rep.ok, f"{rep.verdict}: {rep.detail}"
    for name, rung in census["rungs"].items():
        for f in SYMMETRY_REQUIRED:
            assert f in rung, f"rung {name} is missing {f}"


@pytest.mark.skipif(not CENSUS.exists(), reason="census artifact absent (gitignored / fresh clone)")
def test_live_census_utilization_is_known_for_every_rung():
    """A rung whose utilization is `unknown` is one nobody can reason about -- which is how
    'no rows available' happened. Every rung must carry a measurable figure or an explicit reason."""
    census = json.loads(CENSUS.read_text(encoding="utf-8"))
    rep = audit(census)
    unknown = [u["rung"] for u in rep["utilization"] if u["verdict"].startswith("unknown:")]
    assert not unknown, f"utilization unknown for {unknown}"
