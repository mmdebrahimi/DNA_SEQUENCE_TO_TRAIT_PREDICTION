"""VME/ME — the confusion matrix in the direction that can hurt someone.

These pin the two things that make the layer trustworthy rather than decorative: the arithmetic identity
(so nobody can later "improve" VME into a different quantity while keeping the name), and the separate
denominators (so the asymmetry the convention exists to expose cannot be pooled away).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.eval.error_rates import vme_me, worst_vme_first  # noqa: E402

CARD = ROOT / "wiki" / "decoder_validation_report_card.json"


def test_vme_is_exactly_one_minus_sensitivity():
    """The identity IS the claim. If a future edit makes VME something else, the name becomes a lie."""
    tp, fp, tn, fn = 14, 3, 27, 16
    r = vme_me(tp, fp, tn, fn)
    sens, spec = tp / (tp + fn), tn / (tn + fp)
    assert r.vme == pytest.approx(1 - sens)
    assert r.me == pytest.approx(1 - spec)


def test_each_rate_carries_its_own_denominator():
    """Pooling them into one 'error rate' would hide exactly the asymmetry the convention exposes -- a
    cell can be near-perfect in one direction and dangerous in the other."""
    r = vme_me(tp=10, fp=1, tn=99, fn=10)
    assert r.vme_n == 20 and r.me_n == 100
    assert r.vme == pytest.approx(0.5) and r.me == pytest.approx(0.01)


def test_a_zero_denominator_is_none_not_zero():
    """A single-class cohort has no rate in the missing direction. Reporting 0.0 would read as a perfect
    score on a question nobody asked."""
    r = vme_me(tp=5, fp=0, tn=0, fn=0)
    assert r.vme == pytest.approx(0.0) and r.vme_n == 5
    assert r.me is None and r.me_n == 0 and r.me_ci is None


def test_absent_counts_return_none_rather_than_guessing():
    assert vme_me(1, 2, 3, None) is None
    assert vme_me(None, None, None, None) is None


def test_ordering_puts_the_dangerous_cell_first():
    """Alphabetical order would put a 0.533 wherever the alphabet happened to put it."""
    cells = [{"organism": "a", "drug": "x", "vme": 0.0},
             {"organism": "z", "drug": "y", "vme": 0.53},
             {"organism": "m", "drug": "w", "vme": 0.2}]
    assert [c["vme"] for c in worst_vme_first(cells)] == [0.53, 0.2, 0.0]


def test_a_cell_without_vme_is_dropped_from_the_ranking_not_sorted_as_zero():
    """Treating a missing rate as 0.0 would rank an UNMEASURED cell as the safest on the card."""
    assert worst_vme_first([{"organism": "a", "drug": "x", "vme": None}]) == []


# --- the committed card ------------------------------------------------------------------------

@pytest.mark.skipif(not CARD.exists(), reason="report card not generated")
def test_the_card_carries_error_rates_on_every_scored_cell_and_only_those():
    d = json.loads(CARD.read_text(encoding="utf-8"))
    scored = [c for c in d["cells"] if c.get("state") == "SCORED"]
    assert scored, "no SCORED cells; the check would be vacuous"
    assert all(c.get("error_rates") for c in scored)
    others = [c for c in d["cells"] if c.get("state") != "SCORED"]
    assert not any(c.get("error_rates") for c in others), (
        "a non-SCORED cell has error rates -- it has no confusion matrix to derive them from")


@pytest.mark.skipif(not CARD.exists(), reason="report card not generated")
def test_the_cards_own_rates_reconcile_with_its_own_sens_and_spec():
    """Guards the wiring, not the maths: a block computed from the wrong cell's counts would still be
    internally consistent, and only reconciling against THAT cell's published sens/spec catches it.

    Tolerance is 1e-3 because the card PUBLISHES sens/spec rounded to 3dp while error_rates keeps full
    precision (1 - 0.967 = 0.033 vs an exact 1/30 = 0.0333). Loose enough for the rounding, still far
    tighter than the gap between any two cells on the card, which is what a mis-wiring would look like.
    """
    d = json.loads(CARD.read_text(encoding="utf-8"))
    for c in d["cells"]:
        er = c.get("error_rates")
        if not er:
            continue
        if c.get("sens") is not None and er["vme"] is not None:
            assert er["vme"] == pytest.approx(1 - c["sens"], abs=1e-3), (c["organism"], c["drug"])
        if c.get("spec") is not None and er["me"] is not None:
            assert er["me"] == pytest.approx(1 - c["spec"], abs=1e-3), (c["organism"], c["drug"])


@pytest.mark.skipif(not CARD.exists(), reason="report card not generated")
def test_no_acceptance_bar_is_asserted():
    """Regulatory VME/ME ceilings exist; their exact values are not verifiable from this repo. Writing a
    remembered number beside a real measurement is the fabrication hazard, so the card must not."""
    md = (ROOT / "wiki" / "decoder_validation_report_card.md").read_text(encoding="utf-8")
    assert "No acceptance bar is asserted" in md
    import re
    section = md.split("## Clinical error rates", 1)[1].split("\n## ", 1)[0]
    banned = re.findall(r"(?:VME|ME)\s*(?:<=|<|≤)\s*\d", section)
    assert not banned, f"the section asserts a numeric acceptance threshold: {banned}"
