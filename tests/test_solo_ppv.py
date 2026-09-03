"""Pins for SOLO-PPV grading: the Wilson interval, the WHO bar, and the co-carriage split."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "wiki" / "solo_ppv_2026-09-03.json"


def _mod():
    spec = importlib.util.spec_from_file_location("solo_ppv", ROOT / "scripts" / "solo_ppv.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_wilson_does_not_assert_certainty_on_a_perfect_run():
    """A Wald interval on 99/99 is [1.0, 1.0] -- certainty the data cannot support. Wilson is why the
    headline reads 'CI lower 0.963' rather than 'PPV 1.0, no interval'."""
    lo, hi = _mod().wilson(99, 99)
    assert lo < 1.0 and hi == pytest.approx(1.0)
    assert 0.95 < lo < 0.98


def test_wilson_is_wide_on_tiny_n_and_narrows_with_evidence():
    m = _mod()
    small_lo, _ = m.wilson(3, 3)
    big_lo, _ = m.wilson(300, 300)
    assert small_lo < big_lo


def test_the_who_bar_needs_both_enough_solos_and_a_clearing_interval():
    m = _mod()
    assert m.score(4, 0)["meets_who_grade1_bar"] is False      # only 4 solo occurrences
    assert m.score(5, 0)["meets_who_grade1_bar"] is True
    assert m.score(2, 30)["meets_who_grade1_bar"] is False     # plenty of n, PPV far too low


def test_an_empty_stratum_reports_no_ppv_rather_than_a_number():
    s = _mod().score(0, 0)
    assert s["n"] == 0 and s["ppv"] is None


@pytest.mark.skipif(not ART.is_file(), reason="solo artifact not present")
def test_solo_denominator_is_strictly_smaller_than_pooled():
    """The whole point: solo DISCARDS co-carriage isolates. If it did not, it would not be controlling
    for anything."""
    d = json.loads(ART.read_text(encoding="utf-8"))["result"]["strata"]
    assert d["solo_no_aac3"]["n"] < d["pooled_all_carriers"]["n"]


@pytest.mark.skipif(not ART.is_file(), reason="solo artifact not present")
def test_the_strata_reconcile_with_the_pooled_total():
    """Arithmetic check -- a mis-split would silently move isolates between strata."""
    d = json.loads(ART.read_text(encoding="utf-8"))["result"]["strata"]
    assert (d["solo_no_aac3"]["resistant"] + d["co_carriage_with_aac3"]["resistant"]
            == d["pooled_all_carriers"]["resistant"])
    assert (d["solo_no_aac3"]["susceptible"] + d["co_carriage_with_aac3"]["susceptible"]
            == d["pooled_all_carriers"]["susceptible"])


@pytest.mark.skipif(not ART.is_file(), reason="solo artifact not present")
def test_the_artifact_names_solo_as_the_conservative_estimator():
    """Solo is not universally better -- penalised regression beats it on TB sensitivity. Claiming it
    as strictly superior would be the overclaim."""
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert any("CONSERVATIVE" in lim for lim in d["honest_limits"])
    assert any("aac(3)" in lim and "NOT a full co-determinant screen" in lim
               for lim in d["honest_limits"])
