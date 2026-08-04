"""FBA cell tests: pure logic (wheel-only, no cobra) + real-iML1515 smoke (needs cobra, marked slow)."""
from __future__ import annotations

import pytest

from dna_decode.fba.keio import confusion, metrics_from_confusion, parse_keio_fitness
from dna_decode.fba.model import call_essential


# ---- pure logic: essential-call threshold (no cobra import) ----

def test_call_essential_zero_growth_is_essential():
    assert call_essential(0.0, 0.877) is True
    assert call_essential(1e-9, 0.877) is True


def test_call_essential_full_growth_is_not_essential():
    assert call_essential(0.877, 0.877) is False
    assert call_essential(0.5, 0.877) is False


def test_call_essential_below_one_percent_is_essential():
    # 0.5% of WT -> essential; 5% -> viable (default frac=0.01)
    assert call_essential(0.004, 0.877) is True
    assert call_essential(0.05, 0.877) is False


def test_call_essential_none_or_nan_is_essential():
    assert call_essential(None, 0.877) is True
    assert call_essential(float("nan"), 0.877) is True


def test_call_essential_dead_wildtype_edge():
    # WT can't grow -> only ~0 KO growth counts as essential
    assert call_essential(0.0, 0.0) is True
    assert call_essential(0.1, 0.0) is False


# ---- pure logic: Keio fitness parsing ----

_TSV = (
    "orgId\tlocusId\tsysName\tgeneName\tdesc\tset1 D-Glucose (C)\tset2 D-Glucose (C)\tset3 L-Arabinose (C)\n"
    "Keio\t1\tb0001\tthrA\tx\t-0.1\t0.1\t-0.2\n"          # glucose mean 0.0 -> non-essential
    "Keio\t2\tb0720\tgltA\tcitrate synthase\t-3.0\t-3.1\t0.4\n"  # glucose mean -3.05 -> essential
    "Keio\t3\tb9999\tzzz\ty\tNA\t\t-1.0\n"                # no usable glucose value -> dropped
)


def test_parse_keio_fitness_glucose():
    labels = parse_keio_fitness(_TSV, carbon="D-Glucose", threshold=-2.0)
    assert set(labels) == {"b0001", "b0720"}          # b9999 dropped (no usable glucose value)
    assert labels["b0001"]["essential"] is False
    assert labels["b0720"]["essential"] is True
    assert labels["b0720"]["fitness"] == pytest.approx(-3.05, abs=1e-6)


def test_parse_keio_fitness_carbon_switch():
    labels = parse_keio_fitness(_TSV, carbon="L-Arabinose", threshold=-2.0)
    assert set(labels) == {"b0001", "b0720", "b9999"}  # all have an arabinose value
    assert all(v["essential"] is False for v in labels.values())  # none below -2 on arabinose


# ---- pure logic: confusion + metrics ----

def test_confusion_and_metrics():
    exp = {"a": True, "b": True, "c": False, "d": False, "e": True}
    pred = {"a": True, "b": False, "c": False, "d": True, "e": True}  # e also
    cm = confusion(exp, pred)
    assert cm == {"tp": 2, "fp": 1, "tn": 1, "fn": 1, "n": 5}
    m = metrics_from_confusion(cm)
    assert m["accuracy"] == pytest.approx(0.6)
    assert m["precision"] == pytest.approx(2 / 3)
    assert m["recall"] == pytest.approx(2 / 3)
    assert -1.0 <= m["mcc"] <= 1.0


def test_confusion_intersection_only():
    # only shared keys are scored
    cm = confusion({"a": True, "x": True}, {"a": True, "y": False})
    assert cm["n"] == 1 and cm["tp"] == 1


# ---- real-model smoke (needs cobra + iML1515; slow) ----

@pytest.mark.slow
def test_iml1515_wildtype_and_known_essential():
    cobra = pytest.importorskip("cobra")  # noqa: F841
    from dna_decode.fba.model import knockout_growth, load_model, wildtype_growth

    m = load_model()
    wt = wildtype_growth(m)
    assert 0.5 < wt < 1.2  # iML1515 glucose-M9 growth ~0.877 /h
    # gltA (b0720, citrate synthase) is essential on glucose minimal
    g = knockout_growth(m, "b0720")
    assert call_essential(g, wt) is True
    # a transporter/peripheral gene should stay viable (pick one that is non-essential): pgi b4025
    g2 = knockout_growth(m, "b4025")
    assert call_essential(g2, wt) is False
