"""Tests for scripts/cipro_mic_audit.py — pure-logic functions only.

Pins the tier-classification + MIC-parsing logic that drives the EP1 audit
infrastructure verdict (HIGH_R / HIGH_S / DECISIVE / BORDERLINE / AMBIGUOUS /
CONFLICT / NO_MIC). The cohort-load + CSV-read + packet-write paths are
orchestration (skipped). Tests cover the breakpoint thresholds, conflict
detection, ambiguous (CLSI vs EUCAST disagreement), and the decisive subset
output.
"""
from __future__ import annotations

import pytest

from scripts.cipro_mic_audit import (
    CLSI_R,
    CLSI_S,
    EUCAST_R,
    EUCAST_S,
    _confidence_tier,
    _parse_mic,
)


# ---- _parse_mic --------------------------------------------------------------


def test_parse_mic_plain_numeric():
    assert _parse_mic("8.0", "") == 8.0


def test_parse_mic_strips_leading_operator():
    assert _parse_mic("<=2", "") == 2.0
    assert _parse_mic(">32", "") == 32.0
    assert _parse_mic("=0.5", "") == 0.5


def test_parse_mic_empty_returns_none():
    assert _parse_mic("", "") is None


def test_parse_mic_na_variants_return_none():
    for marker in ("NA", "N/A", "NULL", "NONE", "-", "n/a", "na"):
        assert _parse_mic(marker, "") is None, marker


def test_parse_mic_unparseable_returns_none():
    assert _parse_mic("not a number", "") is None


# ---- _confidence_tier --------------------------------------------------------


def test_tier_high_r_at_extreme_mic():
    # MIC >= 4 * CLSI_R = 8.0 -> HIGH_R
    tier, detail = _confidence_tier([8.0, 16.0], ["R", "Resistant"])
    assert tier == "HIGH_R"
    assert detail["median_mic"] == 12.0
    assert detail["clsi_call"] == "R"
    assert detail["eucast_call"] == "R"


def test_tier_high_s_at_extreme_mic():
    # MIC <= CLSI_S / 4 = 0.125 -> HIGH_S
    tier, detail = _confidence_tier([0.06, 0.125], ["S", "S"])
    assert tier == "HIGH_S"
    assert detail["clsi_call"] == "S"
    assert detail["eucast_call"] == "S"


def test_tier_borderline_in_narrow_band():
    # MIC in [CLSI_S/2, 2*CLSI_R] = [0.25, 4.0] -> BORDERLINE (when CLSI/EUCAST agree)
    # 1.0 puts us at CLSI_I (between 0.5 and 2.0) but EUCAST_R (>= 1.0); they disagree
    # So we pick a value where they agree on R: e.g., 2.5 -> CLSI_R AND EUCAST_R
    tier, detail = _confidence_tier([2.5], ["R"])
    assert tier == "BORDERLINE"
    assert detail["clsi_call"] == "R"
    assert detail["eucast_call"] == "R"


def test_tier_ambiguous_clsi_eucast_disagree():
    # MIC = 0.75 -> CLSI_S (<=0.5? no, 0.75 > 0.5; CLSI_I since 0.5<0.75<2.0)
    # EUCAST_S (0.75 <= 0.25? no) / EUCAST_R (0.75 >= 1.0? no) / EUCAST_I
    # Need a value where CLSI != EUCAST: 1.0 -> CLSI_I (0.5<1.0<2.0); EUCAST_R (>=1.0)
    tier, detail = _confidence_tier([1.0], ["R"])
    assert tier == "AMBIGUOUS"
    assert detail["clsi_call"] != detail["eucast_call"]


def test_tier_no_mic_when_empty():
    tier, detail = _confidence_tier([], ["R"])
    assert tier == "NO_MIC"
    assert detail == {}


def test_tier_conflict_when_calls_disagree():
    # R + S calls in same strain -> CONFLICT regardless of MIC values
    tier, detail = _confidence_tier([8.0], ["R", "S"])
    assert tier == "CONFLICT"
    assert "S" in detail["distinct_calls"] or "SUSCEPTIBLE" in detail["distinct_calls"]


def test_tier_conflict_takes_precedence_over_high_r():
    # Even with HIGH_R-level MICs, conflicting calls win
    tier, _ = _confidence_tier([8.0, 16.0], ["Resistant", "Susceptible"])
    assert tier == "CONFLICT"


def test_tier_detail_includes_distance_metrics():
    tier, detail = _confidence_tier([8.0], ["R"])
    assert tier == "HIGH_R"
    assert detail["distance_to_clsi_r"] == 8.0 / CLSI_R  # 4.0
    assert detail["distance_to_clsi_s"] == 8.0 / CLSI_S  # 16.0


def test_tier_uses_median_for_multiple_mics():
    # 3 MICs -> median = 4.0 -> in BORDERLINE band [0.25, 4.0]
    tier, detail = _confidence_tier([1.0, 4.0, 16.0], ["R", "R", "R"])
    assert detail["median_mic"] == 4.0


def test_tier_decisive_outside_borderline_band():
    # MIC > 2*CLSI_R but < 4*CLSI_R -> DECISIVE (between borderline and HIGH_R)
    # 2*CLSI_R = 4.0; 4*CLSI_R = 8.0; so MIC in (4.0, 8.0) -> DECISIVE
    tier, detail = _confidence_tier([5.0], ["R"])
    assert tier == "DECISIVE"
    assert detail["clsi_call"] == "R"
    assert detail["eucast_call"] == "R"


def test_breakpoint_constants_match_documented_values():
    # Lock the breakpoints as part of the regression-guard contract
    assert CLSI_R == 2.0
    assert CLSI_S == 0.5
    assert EUCAST_R == 1.0
    assert EUCAST_S == 0.25
