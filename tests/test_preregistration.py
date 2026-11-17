"""Tests for the pre-registration validator.

Each test names the project failure it guards. Several use the REAL numbers from the failure, so a
regression is visible as a changed headline rather than as an abstract assertion.
"""
from __future__ import annotations

import math

import pytest

from dna_decode.eval.preregistration import (
    DEFAULT_MAX_CLEAN_FAILURE,
    MIN_UNITS,
    Bar,
    IndependenceUnit,
    PreregistrationRefused,
    Statistic,
    assert_commensurable,
    assert_non_vacuous,
    assert_not_a_function_of_the_design,
    assert_range_reported_with_aggregate,
    clean_failure_rate,
    derive_bar,
    validate_bar,
)

# The measured active frame, 2026-11-17. Reproduced to the row: t.active == "active" selects
# exactly the 14,913 tiles the gate recorded.
ACTIVE_FRAME_ICC = 0.6931
ACTIVE_FRAME_ROWS = 3822
ACTIVE_FRAME_UNITS = 462


def active_frame_unit() -> IndependenceUnit:
    return IndependenceUnit.from_icc(
        "peak",
        n_rows=ACTIVE_FRAME_ROWS,
        n_units=ACTIVE_FRAME_UNITS,
        icc=ACTIVE_FRAME_ICC,
        measured_from="one-way ANOVA ICC of rep1 within peak, 14,913-tile active frame",
    )


# --------------------------------------------------------------------------------------
# F3 -- the independence unit must be declared and MEASURED
# --------------------------------------------------------------------------------------
def test_design_effect_has_no_default_so_it_cannot_be_silently_one() -> None:
    """The whole of F3 is that design_effect=1.0 was assumed, never stated. Constructing a unit
    without it must be a TypeError, not a quiet 1.0."""
    with pytest.raises(TypeError):
        IndependenceUnit("peak", 3822, 462)  # type: ignore[call-arg]


def test_design_effect_below_one_is_refused() -> None:
    with pytest.raises(PreregistrationRefused, match="cannot ADD independent"):
        IndependenceUnit("peak", 3822, 462, 0.5, "somewhere")


def test_unmeasured_design_effect_is_refused_even_when_it_is_one() -> None:
    """A caller claiming independent rows must still say where that was established."""
    with pytest.raises(PreregistrationRefused, match="where it was MEASURED"):
        IndependenceUnit("row", 100, 100, 1.0, "")


def test_unit_coarser_than_the_row_is_refused() -> None:
    with pytest.raises(PreregistrationRefused, match="exceeds n_rows"):
        IndependenceUnit("peak", 462, 3822, 1.0, "x")


def test_from_icc_reproduces_the_measured_active_frame_design_effect() -> None:
    unit = active_frame_unit()
    m = ACTIVE_FRAME_ROWS / ACTIVE_FRAME_UNITS
    assert unit.design_effect == pytest.approx(1.0 + (m - 1.0) * ACTIVE_FRAME_ICC)
    # the headline: 3,822 rows are a few hundred independent units, not a few thousand
    assert 600 < unit.n_effective < 700
    assert unit.se_spearman_null == pytest.approx(
        1.0 / math.sqrt(unit.n_effective - 1.0)
    )


def test_icc_outside_the_unit_interval_is_refused() -> None:
    with pytest.raises(PreregistrationRefused, match="outside"):
        IndependenceUnit.from_icc("peak", n_rows=10, n_units=5, icc=1.5, measured_from="x")


# --------------------------------------------------------------------------------------
# the load-bearing number: a clean pipeline's failure rate
# --------------------------------------------------------------------------------------
def test_the_frozen_null_bar_fails_a_clean_pipeline_about_ninety_percent_of_the_time() -> None:
    """THE finding this module exists to make visible. The frozen 0.05 bar is not a strict test --
    it is one a clean pipeline almost cannot pass, so INDETERMINATE carried no information."""
    stat = Statistic("encoder_shuffled", "spearman", "max-of-k", 10, "two-sided")
    cfr = clean_failure_rate(stat, active_frame_unit(), 0.05)
    assert cfr > 0.85, f"expected a near-certain clean failure, got {cfr}"


def test_a_derived_bar_has_the_clean_failure_rate_it_was_asked_for() -> None:
    stat = Statistic("encoder_shuffled", "spearman", "max-of-k", 10, "two-sided")
    unit = active_frame_unit()
    bar = derive_bar(stat, unit, max_clean_failure=DEFAULT_MAX_CLEAN_FAILURE)
    assert clean_failure_rate(stat, unit, bar) == pytest.approx(
        DEFAULT_MAX_CLEAN_FAILURE, abs=0.01
    )


# --------------------------------------------------------------------------------------
# F4 -- the bar must be calibrated for the AGGREGATE actually computed
# --------------------------------------------------------------------------------------
def test_a_max_of_k_bar_is_strictly_looser_than_a_single_draw_bar() -> None:
    """Bounding a max-of-10 with a single-draw bar is F4. The two differ by construction, so a bar
    that ignores the aggregation is testing something else."""
    unit = active_frame_unit()
    single = derive_bar(Statistic("s", "spearman", "single", 1, "two-sided"), unit)
    max10 = derive_bar(Statistic("s", "spearman", "max-of-k", 10, "two-sided"), unit)
    assert max10 > single


def test_single_aggregation_requires_k_of_one_and_vice_versa() -> None:
    with pytest.raises(PreregistrationRefused, match="requires k == 1"):
        Statistic("s", "spearman", "single", 10, "upward")
    with pytest.raises(PreregistrationRefused, match="requires k >= 2"):
        Statistic("s", "spearman", "max-of-k", 1, "upward")


# --------------------------------------------------------------------------------------
# F5 -- sidedness must match the hazard
# --------------------------------------------------------------------------------------
def test_a_one_sided_hazard_gives_a_tighter_bar_than_a_two_sided_one() -> None:
    unit = active_frame_unit()
    one = derive_bar(Statistic("s", "spearman", "max-of-k", 10, "upward"), unit)
    two = derive_bar(Statistic("s", "spearman", "max-of-k", 10, "two-sided"), unit)
    assert one < two


def test_validate_flags_a_two_sided_bar_for_confirmation() -> None:
    stat = Statistic("encoder_shuffled", "spearman", "max-of-k", 10, "two-sided")
    v = validate_bar(stat, active_frame_unit(), Bar(0.5, derived_from="simulation"))
    assert any("two-sided" in r for r in v.reasons)


# --------------------------------------------------------------------------------------
# F2 -- an inherited bar must declare its origin design
# --------------------------------------------------------------------------------------
def test_an_inherited_bar_is_refused_even_when_its_value_would_pass() -> None:
    """A generous inherited bar still hides a design mismatch, so inheritance is flagged
    independently of the clean-failure rate."""
    stat = Statistic("m", "difference-of-spearman", "single", 1, "upward")
    v = validate_bar(
        stat,
        active_frame_unit(),
        Bar(0.9, derived_from="simulation", inherited_from="the pooled 44,106-row run"),
    )
    assert v.clean_failure_rate == pytest.approx(0.0, abs=1e-6)
    assert not v.ok
    assert any("INHERITED" in r for r in v.reasons)


def test_a_well_specified_bar_passes_so_the_validator_is_not_vacuous() -> None:
    """A validator that refuses everything is as useless as one that refuses nothing."""
    stat = Statistic("encoder_shuffled", "spearman", "max-of-k", 10, "upward")
    unit = active_frame_unit()
    good = Bar(derive_bar(stat, unit), derived_from="derive_bar at the measured design")
    v = validate_bar(stat, unit, good)
    assert v.ok, v.reasons


def test_a_bar_without_provenance_is_refused() -> None:
    stat = Statistic("s", "spearman", "single", 1, "upward")
    unit = active_frame_unit()
    v = validate_bar(stat, unit, Bar(derive_bar(stat, unit), derived_from=""))
    assert any("derived_from" in r for r in v.reasons)


# --------------------------------------------------------------------------------------
# refusal rather than extrapolation at small n
# --------------------------------------------------------------------------------------
def test_derive_bar_refuses_below_the_minimum_unit_count() -> None:
    tiny = IndependenceUnit("peak", 100, MIN_UNITS - 1, 2.0, "measured")
    with pytest.raises(PreregistrationRefused, match="permutation null"):
        derive_bar(Statistic("s", "spearman", "single", 1, "upward"), tiny)


# --------------------------------------------------------------------------------------
# F6 -- the retracted statistic, used as the test case
# --------------------------------------------------------------------------------------
def test_the_actually_retracted_statistic_is_caught() -> None:
    """rho = sqrt(3k(n-k)/(n^2-1)) is a closed form in the split proportion alone. It was
    published from this project and retracted. Fed same-design / different-data inputs it returns
    one value, which is exactly what this check refuses."""

    def retracted(dataset: object) -> float:
        n, k, _values = dataset  # type: ignore[misc]
        return math.sqrt(3.0 * k * (n - k) / (n * n - 1.0))

    same_design_different_data = [
        (44106, 14913, [1.0, 2.0, 3.0]),
        (44106, 14913, [9.0, -4.0, 0.5]),
    ]
    with pytest.raises(PreregistrationRefused, match="function of the DESIGN"):
        assert_not_a_function_of_the_design(retracted, same_design_different_data)


def test_a_statistic_that_does_read_the_data_passes() -> None:
    def reads_data(dataset: object) -> float:
        _n, _k, values = dataset  # type: ignore[misc]
        return sum(values) / len(values)

    assert_not_a_function_of_the_design(
        reads_data, [(10, 5, [1.0, 2.0]), (10, 5, [3.0, 4.0])]
    )


def test_one_dataset_cannot_demonstrate_variation() -> None:
    with pytest.raises(PreregistrationRefused, match="at least 2"):
        assert_not_a_function_of_the_design(lambda d: 1.0, [(1, 1, [1.0])])


# --------------------------------------------------------------------------------------
# F7 -- a control that removes nothing
# --------------------------------------------------------------------------------------
def test_a_filter_that_removes_nothing_is_refused() -> None:
    with pytest.raises(PreregistrationRefused, match="removed 0 of"):
        assert_non_vacuous(1818, 1818, what="the POINT-row filter")


def test_a_filter_that_removes_something_passes() -> None:
    assert_non_vacuous(4951, 1527, what="the position restriction")


# --------------------------------------------------------------------------------------
# F9 -- incommensurable metrics
# --------------------------------------------------------------------------------------
def test_a_fraction_beside_a_spearman_is_refused() -> None:
    with pytest.raises(PreregistrationRefused, match="not commensurable"):
        assert_commensurable(
            {
                "fragment-rung": "fraction-of-ceiling",
                "asgiven-position-resolved": "spearman",
            }
        )


def test_like_for_like_metrics_pass() -> None:
    assert_commensurable({"a": "spearman", "b": "spearman"})


# --------------------------------------------------------------------------------------
# F8 -- an aggregate quoted as a scalar
# --------------------------------------------------------------------------------------
def test_the_real_pooled_ceiling_is_refused_as_a_scalar() -> None:
    """0.941 pooled over deciles whose own ceilings run 0.5013 to 0.9823."""
    with pytest.raises(PreregistrationRefused, match="spans strata"):
        assert_range_reported_with_aggregate(
            0.941, [0.5013, 0.62, 0.71, 0.80, 0.86, 0.90, 0.93, 0.95, 0.97, 0.9823],
            what="the replicate noise ceiling",
        )


def test_a_tight_aggregate_may_stand_alone() -> None:
    assert_range_reported_with_aggregate(
        0.50, [0.49, 0.50, 0.51], what="a homogeneous aggregate"
    )


def test_an_aggregate_needs_its_strata_to_be_judged() -> None:
    with pytest.raises(PreregistrationRefused, match="per-stratum"):
        assert_range_reported_with_aggregate(0.941, [0.941], what="one stratum")
