"""Guards for the deep-research harvest + its two-screen driver.

The load-bearing tests are the ones that pin **what cannot be called promising**. A literature campaign's
failure mode is not a wrong number, it is an unscreened lead written up as a finding — which is exactly the
pattern this project recorded when most `/innovate` survivors fell as soon as they were acted on, and when a
run proposed re-scoring gLM2 in a regime cell already measured as `LOSES_TO_CATALOG`.

So: unscreened is never passed, a method with no regime triple is never promising, and a locator the
harvester could not confirm is never promising.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.deep_research_screen import (
    LEAD_ONLY,
    PROMISING,
    REJECTED,
    UNSCREENABLE,
    screen_harvest,
    screen_one,
)

REPO = Path(__file__).resolve().parents[1]
HARVEST = REPO / "wiki" / "deep_research_candidates_2026-09-28.json"


# --- classifier: the three ways a candidate must FAIL to be promising --------------------------

def test_a_dataset_with_no_label_screen_is_lead_only_never_promising():
    """Unscreened is not passed. This is the whole point of the class."""
    out = screen_one({"name": "x", "kind": "dataset", "locator": "https://doi.org/10.1000/x"})
    assert out["label_verdict"] is None
    assert out["screen_class"] == LEAD_ONLY


def test_a_method_with_no_regime_triple_is_unscreenable():
    """The gLM2 hole: the regime screen is precisely what catches a re-run of a recorded negative, so a
    method that was never put through it cannot be promising."""
    out = screen_one({"name": "m", "kind": "method", "locator": "https://doi.org/10.1000/m"})
    assert out["screen_class"] == UNSCREENABLE


@pytest.mark.parametrize("locator", ["", "   ", "UNVERIFIED possibly PRJNA999999", "TODO find doi"])
def test_a_missing_or_unverified_locator_is_unscreenable(locator):
    """An unverified accession is the fabrication hazard this project guards hardest. It is preserved and
    it blocks promotion — never silently dropped to make a candidate look citable."""
    out = screen_one({"name": "x", "kind": "dataset", "locator": locator,
                      "regime_screen": {"population": "constructed", "endpoint": "molecular",
                                        "method": "supervised"}})
    assert out["screen_class"] == UNSCREENABLE, out
    assert out["locator_ok"] is False


def test_a_recorded_negative_regime_is_rejected_even_with_a_good_locator():
    """natural x molecular x zero_shot is LOSES_TO_CATALOG — measured, not asserted."""
    out = screen_one({"name": "another likelihood scorer", "kind": "method",
                      "locator": "https://doi.org/10.1000/plm",
                      "regime_screen": {"population": "natural", "endpoint": "molecular",
                                        "method": "zero_shot"}})
    assert out["regime_verdict"] == "LOSES_TO_CATALOG"
    assert out["screen_class"] == REJECTED


def test_the_natural_supervised_cell_is_NOT_treated_as_closed():
    """Scope rail: the 0-for-5 negative is ZERO-SHOT-scoped. Compressing it to 'natural populations are
    closed' has hidden a live direction three separate times in this project, so a supervised
    natural-population proposal must come back open-with-conditions, never rejected."""
    out = screen_one({"name": "supervised within-lineage design", "kind": "method",
                      "locator": "https://doi.org/10.1000/abc",
                      "regime_screen": {"population": "natural", "endpoint": "organism",
                                        "method": "supervised"}})
    assert out["regime_verdict"] == "REQUIRES_DECONFOUNDING"
    assert out["screen_class"] == PROMISING


# --- the driver reproduces the three COMMITTED hand verdicts ------------------------------------

def test_the_driver_reproduces_the_committed_candidate_verdicts():
    """Wiring check against known answers: PEAR clears, Oxford and HBV are rejected. If this drifts, the
    driver is wrong — not the harvest."""
    import scripts.screen_candidate_gates as scg

    for name, expected in (("PEAR", "CLEARS"), ("OXFORD", "REJECTED"), ("HBV", "REJECTED")):
        spec = getattr(scg, name)
        out = screen_one({"name": spec["candidate"], "kind": "dataset",
                          "locator": "https://doi.org/10.1000/committed",
                          "label_screen": {"intended_layer": spec["intended_layer"],
                                           "evidence": spec["evidence"]}})
        assert out["label_verdict"] == expected, (name, out["label_verdict"], out["label_reason"])


# --- the harvest artifact itself ----------------------------------------------------------------

@pytest.mark.skipif(not HARVEST.exists(), reason="harvest not produced yet")
def test_harvest_is_wellformed_and_every_candidate_is_accounted_for():
    doc = json.loads(HARVEST.read_text(encoding="utf-8"))
    cands = doc.get("candidates")
    assert isinstance(cands, list) and cands, "harvest carries no candidates"

    ids = [c.get("id") for c in cands]
    assert len(ids) == len(set(ids)), "duplicate candidate ids"
    for c in cands:
        for req in ("id", "family", "kind", "name", "locator"):
            assert c.get(req) not in (None, ""), f"{c.get('id')}: missing {req}"
        assert c["family"] in {"F1", "F2", "F3", "F4"}, c["family"]
        assert c["kind"] in {"dataset", "method", "thesis"}, c["kind"]

    rep = screen_harvest(HARVEST)
    # Every candidate must land in exactly one class, and the counts must reconcile with the list --
    # a summary that does not sum to its own rows is how a dropped candidate hides.
    assert sum(rep["class_counts"].values()) == rep["n_candidates"] == len(cands)
    for s in rep["candidates"]:
        assert s["screen_class"] in {PROMISING, LEAD_ONLY, REJECTED, UNSCREENABLE}


@pytest.mark.skipif(not HARVEST.exists(), reason="harvest not produced yet")
def test_no_candidate_is_promising_without_a_screen_that_actually_ran():
    """The invariant the whole campaign rests on."""
    rep = screen_harvest(HARVEST)
    for s in rep["candidates"]:
        if s["screen_class"] != PROMISING:
            continue
        assert s["locator_ok"], s["name"]
        ran = (s["label_verdict"] is not None) or (s["regime_verdict"] is not None)
        assert ran, f"{s['name']} is PROMISING with neither screen having run"
