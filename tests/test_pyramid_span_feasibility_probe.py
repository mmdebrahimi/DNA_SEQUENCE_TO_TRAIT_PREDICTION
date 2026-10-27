"""Guards for the pyramid-span feasibility probe.

This probe is an ARM'S EXECUTED EVIDENCE, so the tests that matter are the ones proving it could
have come out the other way. A probe that can only ever print INFEASIBLE is not a measurement, it is
the prose claim it replaced wearing a JSON costume.

Three failure directions are pinned:
  1. NON-VACUITY -- a rung long enough DOES come out FEASIBLE (`test_..._CAN_return_feasible`).
  2. "did not look" is NOT "infeasible" -- zero measurable rungs returns INDETERMINATE, which is the
     exact distinction the 2026-10-26 prose version collapsed.
  3. the honest floor is IMPORTED, so a 1 bp branch cannot be smuggled in to make an arm reachable.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.encoder import PYRAMID_WIDTHS  # noqa: E402
from scripts.pyramid_span_feasibility_probe import (  # noqa: E402
    FINEST_HONEST_WIDTH,
    REQUIRED_DECADES,
    ProbeFault,
    ProbeReport,
    RungSpan,
    _rung_from_stats,
    honest_span_decades,
    render,
    required_bp,
)

ARTIFACT = REPO / "wiki" / "pyramid_span_feasibility_2026-10-27.json"


class _Stats:
    """The fields `_rung_from_stats` reads off a LoadStats. Mirrors the real attribute names."""

    def __init__(self, max_len, n_distinct=1, path="synthetic.txt.gz"):
        self.length_min = min(48, max_len)
        self.length_median = max_len
        self.length_max = max_len
        self.n_distinct_lengths = n_distinct
        self.path = path


# --------------------------------------------------------------------- the arithmetic

def test_the_span_formula_reproduces_the_census_field_exactly():
    """log10(499/3) must reproduce `rungs.fragments_LB.max_span_decades_at_3bp_floor` = 2.22."""
    assert round(honest_span_decades(499), 2) == 2.22


def test_the_span_formula_reproduces_the_measured_sibling_arms_own_name():
    """The sibling is named `pyramid-4scale-1.70decades-150bp`; 150 bp must give 1.70."""
    assert round(honest_span_decades(150), 2) == 1.70


def test_four_decades_above_the_honest_floor_needs_30000_bp():
    assert required_bp() == 30000
    assert required_bp(decades=4.0, finest=1) == 10000, \
        "10,000 bp is the 1 bp-floor figure -- the degenerate branch the encoder rejects"


def test_the_finest_width_is_IMPORTED_from_the_encoder_not_restated():
    """A literal here would let a 1 bp floor be introduced locally and inflate every span."""
    assert FINEST_HONEST_WIDTH == PYRAMID_WIDTHS[0]
    assert FINEST_HONEST_WIDTH == 3


def test_a_sequence_shorter_than_the_finest_branch_REFUSES():
    with pytest.raises(ProbeFault, match="below the finest"):
        honest_span_decades(2)
    with pytest.raises(ProbeFault):
        honest_span_decades(0)


# --------------------------------------------------------------------- non-vacuity

def test_the_probe_CAN_return_feasible_on_a_long_enough_rung():
    """THE non-vacuity proof. Without this the INFEASIBLE verdict proves nothing."""
    rep = ProbeReport()
    rep.rungs.append(_rung_from_stats("synthetic_long", 10, _Stats(required_bp())))
    assert rep.verdict() == "FEASIBLE"
    assert rep.rungs[0].feasible_for_4_decades is True
    assert rep.rungs[0].shortfall_factor is None, "a feasible rung has no shortfall"
    assert rep.rungs[0].honest_span_decades >= REQUIRED_DECADES


def test_one_bp_short_of_the_requirement_is_still_infeasible():
    """The boundary bites, so FEASIBLE is not reachable by rounding."""
    rep = ProbeReport()
    rep.rungs.append(_rung_from_stats("synthetic", 10, _Stats(required_bp() - 1)))
    assert rep.verdict() == "INFEASIBLE_SUBSTRATE_GATED"
    assert rep.rungs[0].shortfall_factor is not None


def test_a_mixed_report_is_FEASIBLE_if_ANY_rung_reaches_it():
    rep = ProbeReport()
    rep.rungs.append(_rung_from_stats("short", 10, _Stats(499)))
    rep.rungs.append(_rung_from_stats("long", 10, _Stats(required_bp())))
    assert rep.verdict() == "FEASIBLE"


# --------------------------------------------------------------------- "did not look" != infeasible

def test_zero_measurable_rungs_is_INDETERMINATE_and_never_infeasible():
    """The 2026-10-26 prose claim collapsed exactly this distinction."""
    rep = ProbeReport()
    rep.unreachable.append({"rung": "fragments_LB", "reason": "OSError: drive detached"})
    assert rep.verdict() == "INDETERMINATE_NO_RUNG_MEASURED"
    assert "INFEASIBLE" not in rep.verdict()


def test_an_unreachable_rung_is_RECORDED_not_silently_dropped():
    rep = ProbeReport()
    rep.rungs.append(_rung_from_stats("ok", 10, _Stats(499)))
    rep.unreachable.append({"rung": "gone", "reason": "OSError: detached"})
    d = rep.to_dict()
    assert d["rungs_unreachable"] and d["rungs_unreachable"][0]["rung"] == "gone"
    assert "what_was_searched" in d


# --------------------------------------------------------------------- the row-limit guard

def test_a_row_limited_render_REFUSES_to_state_the_verdict_as_established():
    """Measured: at --limit 2000 the fragment max reads 324 bp against a true 499.

    The first version printed INFEASIBLE to stdout while the artifact said
    INDETERMINATE_ROW_LIMITED -- two surfaces of one run disagreeing, the louder one wrong.
    """
    rep = ProbeReport()
    rep.rungs.append(_rung_from_stats("fragments_LB", 2000, _Stats(324)))
    capped = render(rep, limit=2000)
    assert "INDETERMINATE_ROW_LIMITED" in capped
    assert "cannot ground a verdict" in capped
    uncapped = render(rep, limit=None)
    assert "INDETERMINATE_ROW_LIMITED" not in uncapped
    assert "INFEASIBLE_SUBSTRATE_GATED" in uncapped


# --------------------------------------------------------------------- the committed artifact

@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe artifact not present")
def test_the_committed_artifact_records_a_REAL_uncapped_measurement():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["schema"] == "pyramid-span-feasibility-v1"
    assert "row_limit_applied" not in d, "a capped run must not back the committed verdict"
    assert d["verdict"] == "INFEASIBLE_SUBSTRATE_GATED"
    assert d["required_sequence_length_bp"] == 30000
    assert d["finest_honest_width_bp"] == PYRAMID_WIDTHS[0]

    by_rung = {r["rung"]: r for r in d["rungs_measured"]}
    assert set(by_rung) == {"fragments_LB", "fragments_M9", "peak_tiles"}, \
        "all three canonical rungs must be covered, so the search is bounded and stated"
    # Full-file row counts -- a capped run would show fewer and a lower max.
    assert by_rung["fragments_LB"]["n_sequences"] == 321123
    assert by_rung["fragments_LB"]["length_max"] == 499
    assert by_rung["fragments_M9"]["length_max"] == 499
    assert by_rung["peak_tiles"]["length_max"] == 150
    assert all(not r["feasible_for_4_decades"] for r in d["rungs_measured"])


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe artifact not present")
def test_the_artifact_names_the_REACHABLE_span_not_only_the_shortfall():
    """A verdict that only says 'no' leaves the next session to re-derive the re-scope target."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    best = d["best_achievable"]
    assert best["rung"] == "fragments_LB"
    assert best["max_sequence_length_bp"] == 499
    assert round(best["honest_span_decades"], 2) == 2.22
    assert math.isclose(best["honest_span_decades"], honest_span_decades(499), abs_tol=1e-4)


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe artifact not present")
def test_the_artifact_carries_R7_SCOPE_and_a_REVIVAL_condition():
    """R7: a kill must say what was tested and what re-opens it."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    scope = d["scope"].lower()
    assert "this host" in scope, "the scope must be the cached substrate, not the world"
    assert "not a claim that no such" in scope
    assert "merit" in scope, "it must disclaim saying anything about the arm's predictive merit"
    assert "30000" in d["revival_condition"] or "30,000" in d["revival_condition"]


def test_a_rung_span_round_trips_through_its_dict():
    r = RungSpan("x", 1, 48, 200, 499, 300, honest_span_decades(499), False, 60.1, "f.gz")
    d = r.to_dict()
    assert d["length_max"] == 499 and round(d["honest_span_decades"], 2) == 2.22
