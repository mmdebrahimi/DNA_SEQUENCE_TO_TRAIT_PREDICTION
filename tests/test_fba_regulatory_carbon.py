"""Carbon-panel regulatory arm — pure verdict logic (no solver, no feba.db)."""
from __future__ import annotations

from scripts.fba_regulatory_carbon_test import verdict_for


def test_the_STRONG_null_is_the_binding_test_not_the_weak_one():
    """The whole point of the re-test. A result that only beats the rate-matched null is named as a
    weak-null artifact rather than reported as a lift."""
    assert verdict_for(0.75, 0.0, 0.30) == "REGULATORY_LIFT_IS_A_WEAK_NULL_ARTIFACT"


def test_clearing_the_strong_null_comfortably_confirms_it():
    assert verdict_for(0.75, 0.0, 0.004) == "REGULATORY_LIFT_CONFIRMED_ON_WIDE_PANEL"


def test_a_marginal_pass_is_named_marginal():
    assert verdict_for(0.75, 0.0, 0.03) == "REGULATORY_LIFT_SURVIVES_MARGINALLY"


def test_failing_both_nulls_is_plain_not_supported():
    assert verdict_for(0.75, 0.40, 0.60) == "REGULATORY_LIFT_NOT_SUPPORTED"


def test_a_missing_strong_null_refuses_to_render_a_verdict():
    """Without the binding test there is no verdict -- reporting the weak-null number alone is exactly
    the mistake this script exists to correct."""
    assert verdict_for(0.75, 0.0, None) == "INDETERMINATE"
    assert verdict_for(None, 0.0, 0.01) == "INDETERMINATE"


def test_an_intervention_that_beats_the_null_but_LOSES_to_the_baseline_is_not_a_lift():
    """The bug this function shipped with, caught by verify-in-batch on the real run. Beating a null
    built from the restricted arm's OWN margins says the calls are well-placed GIVEN how many were made
    -- fully compatible with making more calls being a bad idea. The real run: baseline 0.7368 ->
    restricted 0.6839 with strong-null p 0.0, and the first version printed CONFIRMED."""
    assert verdict_for(0.6839, 0.0, 0.0, baseline=0.7368) == "REGULATORY_RESTRICTION_MAKES_IT_WORSE"


def test_the_baseline_check_precedes_the_null_check():
    """Order matters: a worse-than-baseline arm must never reach the null branch at all."""
    assert verdict_for(0.50, 0.0, 0.0, baseline=0.60) == "REGULATORY_RESTRICTION_MAKES_IT_WORSE"
    assert verdict_for(0.70, 0.0, 0.0, baseline=0.60) == "REGULATORY_LIFT_CONFIRMED_ON_WIDE_PANEL"


def test_omitting_the_baseline_preserves_the_old_null_only_semantics():
    assert verdict_for(0.6839, 0.0, 0.0) == "REGULATORY_LIFT_CONFIRMED_ON_WIDE_PANEL"
