"""Pin the axis dynamic-range summariser -- the predictor of knockout-ratio flatness.

Pure; no cobrapy solve, no DB, no network.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from fba_axis_dynamic_range import OBSERVED_FLATNESS, _summarise  # noqa: E402


def test_identical_growths_collapse_to_one_distinct_value():
    """The nitrogen case that motivated all of this: conditions the model literally cannot tell apart."""
    s = _summarise("x", [0.92593] * 6)
    assert s["n_distinct_growths"] == 1 and s["distinct_fraction"] == round(1 / 6, 4)
    assert s["cv"] == 0.0


def test_distinct_counting_is_tolerance_rounded_not_exact_float():
    """Same discipline as the within-gene metric: LP noise in the 12th decimal is not a distinct medium."""
    s = _summarise("x", [0.5, 0.5 + 1e-13, 0.5 - 1e-13, 0.9])
    assert s["n_distinct_growths"] == 2


def test_cv_is_zero_for_a_flat_axis_and_positive_for_a_spread_one():
    assert _summarise("flat", [1.0, 1.0, 1.0])["cv"] == 0.0
    assert _summarise("spread", [0.1, 0.5, 1.0])["cv"] > 0.5


def test_a_dead_axis_reports_none_cv_rather_than_dividing_by_zero():
    s = _summarise("dead", [0.0, 0.0])
    assert s["cv"] is None and s["n_distinct_growths"] == 1


def test_the_two_summaries_are_not_redundant():
    """distinct_fraction asks 'can the model tell these apart at all', cv asks 'by how much'. An axis can
    score 1.0 on the first while barely spreading -- which is why both are reported."""
    a = _summarise("all_distinct_but_narrow", [1.000, 1.001, 1.002, 1.003])
    assert a["distinct_fraction"] == 1.0 and a["cv"] < 0.01


def test_observed_flatness_covers_every_axis_the_probe_reports():
    """A new axis added to the probe without its measured flatness would silently sort as 0 and could
    invert the monotonicity verdict."""
    assert set(OBSERVED_FLATNESS) == {"media4", "carbon", "nitrogen"}
    assert all(0.0 < v < 1.0 for v in OBSERVED_FLATNESS.values())
