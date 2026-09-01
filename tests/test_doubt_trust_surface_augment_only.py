"""The L2 doubt layer is AUGMENT-ONLY on the trust surface: it must not move any cell's tier.

"Does this cell have a known catalog-completeness gap?" and "how well is this cell validated?" are
DIFFERENT questions. Folding one into the other is the shared-key silent-overwrite trap this project
has hit before, so the guard below compares the badge computed WITH the doubt layer against the same
badge computed with it disabled, and requires every pre-existing key to be byte-identical.

A guard like this passes trivially if nothing is ever attached, so the non-vacuity test is not
optional decoration -- it is what makes the other tests mean anything.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.data import trust_surface as ts  # noqa: E402
from dna_decode.eval.amr_rules import DRUG_RULE  # noqa: E402

_TARGET_SITE_DRUGS = ["efavirenz", "nevirapine", "etravirine", "rilpivirine", "doravirine",
                      "nirmatrelvir", "fluconazole", "voriconazole", "lamivudine"]
_ALL = sorted(DRUG_RULE) + _TARGET_SITE_DRUGS


def _badge_without_doubt(drug, monkeypatch):
    monkeypatch.setattr(ts, "doubt_layer_for", lambda _d: None)
    return ts.trust_block(drug)


def test_attaching_the_doubt_layer_changes_no_other_field(monkeypatch):
    """The decisive comparison: same badge, doubt layer on vs off."""
    with_doubt = {d: ts.trust_block(d) for d in _ALL}
    monkeypatch.setattr(ts, "doubt_layer_for", lambda _d: None)
    without = {d: ts.trust_block(d) for d in _ALL}

    for d in _ALL:
        a, b = dict(with_doubt[d]), dict(without[d])
        a.pop("doubt_layer", None)
        assert a == b, f"{d}: the doubt layer altered a pre-existing badge field"


def test_the_guard_is_not_vacuous_some_cell_actually_gets_a_doubt_layer():
    """If nothing ever attaches, the augment-only guard above proves nothing."""
    attached = [d for d in _ALL if "doubt_layer" in ts.trust_block(d)]
    assert attached, "no cell carries a doubt_layer -- the augment-only guard would be vacuous"


def test_the_gentamicin_doubt_block_tracks_the_deployed_rule():
    """Gentamicin was the one cell with an independently-confirmed completeness gap.

    This asserted `known_gap_recovered is True` with an `rmt` family listed until 2026-08-31. The user
    then authorized the v2 lock; the rule now counts `rmt*`/`npmA`, so the gap is CLOSED and the block
    must stop reporting it. Written against the deployed rule so it cannot drift either way.
    """
    import pytest

    from dna_decode.eval.amr_rules import rule_for
    dl = ts.doubt_layer_for("gentamicin")
    if dl is None:
        pytest.skip("doubt-layer artifact not generated on this checkout")
    assert dl["arm"] == "determinant_completeness"
    if rule_for("gentamicin").get("symbol_rescue"):
        assert not any(s.lower().startswith(("rmt", "npma")) for s in dl["strong_families"]), (
            f"the rescue ships but the doubt block still flags {dl['strong_families']}")
        assert dl["known_gap_recovered"] is False, "the closed gap must stop reading as recovered"
    else:
        assert dl["known_gap_recovered"] is True
        assert any(s.lower().startswith("rmt") for s in dl["strong_families"]), dl["strong_families"]


def test_a_cell_with_no_measurement_returns_none_not_a_clean_bill():
    """Never-measured must be distinguishable from measured-and-clean."""
    assert ts.doubt_layer_for("a-drug-that-does-not-exist") is None


def test_tier_is_never_derived_from_the_doubt_layer():
    """A completeness gap is not a validation tier and must never be read as one."""
    for d in _ALL:
        badge = ts.trust_block(d)
        assert badge["tier"] == ts.lookup_trust(d)["tier"], f"{d}: tier moved"
