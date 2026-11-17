"""Refuse a statistical bar that cannot be meaningful, BEFORE the run that it grades.

WHY THIS EXISTS
---------------
Every check here is scoped to a failure this project has actually made and logged. It is not a
generic statistics checklist and it should not grow into one: a check that cannot name a real past
failure is speculation, and this project's own record says a filter that removes nothing is not a
control.

The failures, each mapped to the check that catches it:

F1  THRESHOLD ASSERTED, NOT DERIVED. ``glm_active_frame_gate.PREREGISTERED_ACTIVE_FRAME`` records
    its own thresholds as "ASSERTED". -> ``derive_bar`` computes one from the design instead.
F2  BAR INHERITED ACROSS DESIGNS. ``null_max_abs = 0.05`` was taken verbatim from the pooled
    44,106-row gate and applied to the 14,913-row active frame. -> ``Bar.inherited_from`` must be
    None, or name the design it was derived under so a mismatch is visible.
F3  INDEPENDENCE UNIT UNACCOUNTED. The active frame's 3,822-tile test set is about 648 independent
    units: tiles within a peak share ICC 0.6931, a design effect of 5.900. Nothing in the project
    computed this before 2026-11-17. -> ``IndependenceUnit`` REQUIRES a measured design effect.
    There is deliberately no default of 1.0, because that default is the bug.
F4  BAR CALIBRATED FOR THE WRONG STATISTIC. The null criterion bounds max|rho| over 10 seeds, but
    its bar was set as though for a single draw. -> ``Statistic.aggregation`` plus ``k``, and the
    simulated null is of the AGGREGATE rather than of one draw.
F5  TWO-SIDED BAR ON A ONE-SIDED HAZARD. Leakage can only inflate a null upward; ``abs()`` made the
    test two-sided, and all three breaches were negative. -> ``Statistic.hazard`` drives sidedness.
F6  A STATISTIC THAT IS A FUNCTION OF THE DESIGN. rho = sqrt(3k(n-k)/(n^2-1)) is a closed form in
    the split proportion alone. It was published and retracted. ->
    ``assert_not_a_function_of_the_design``.
F7  A CONTROL THAT REMOVES NOTHING. The HIV fold-censoring filter and the ResFinder POINT-row
    filter each excluded zero rows while being presented as controls. -> ``assert_non_vacuous``.
F8  AN AGGREGATE QUOTED AS A SCALAR. 0.941 is a pooled average over deciles whose own ceilings run
    0.5013 to 0.9823. -> ``assert_range_reported_with_aggregate``.
F9  INCOMMENSURABLE METRICS IN ONE COMPARISON. The fragment arm's 0.185 is a FRACTION sitting in
    the same max() as bacterial Spearman values. -> ``assert_commensurable``.

THE LOAD-BEARING NUMBER IS ``clean_failure_rate``: the probability that a perfectly clean pipeline
FAILS the bar. Nothing in this project computed it before, and it is what separates a test from a
formality. The active frame's frozen null bar has a clean-failure rate of about 0.90, so its
INDETERMINATE verdict carried almost no information about leakage. Crucially this is knowable
BEFORE the run, from the design alone, which is what keeps re-deriving a bar from being
bar-shopping: the argument never mentions the observed result.

HONESTY (this project's OT1 class). This module owns the SIMULATION, the DERIVATION and the
REFUSALS, which are deterministic and data-free. It does NOT own whether the declared unit, hazard
or design effect are TRUE -- a caller can declare ``design_effect=1.0`` and be believed. So the
honest claim is "a code-owned refusal over a model-declared design", never "a bar cannot be
mis-specified here". The null model is asymptotic: Spearman under H0 is treated as
N(0, 1/(n_eff - 1)), standard for large n_eff and unreliable below roughly 20 units, where
``derive_bar`` refuses rather than extrapolating.

TWO NAMED MODELLING GAPS, recorded because each loosens a derived bar rather than tightening it,
and a gap that errs toward passing is the dangerous direction:

G1  ``from_icc`` USES THE LABEL'S ICC ALONE. For a CORRELATION the design effect is closer to
    1 + (m-1)*rho_x*rho_y, where rho_y is the label's intra-cluster correlation and rho_x is the
    PREDICTION's. Passing only rho_y implicitly sets rho_x = 1, which OVERSTATES the design
    effect, understates the effective n, and therefore LOOSENS the derived bar. Treat a bar
    derived this way as an upper bound on the threshold, and prefer a verdict that also survives
    ``design_effect=1.0`` -- the tightest bar available -- so the conclusion does not rest on an
    ICC modelling choice at all.
G2  ``kind`` IS NOW LOAD-BEARING, because it was not. A DIFFERENCE of two correlations measured on
    the same test set is not one Spearman: the two are dependent, so its null variance is not
    1/(n_eff - 1). ``VARIANCE_MODELS`` lists the kinds this module can actually simulate, and
    ``derive_bar`` REFUSES anything else unless the caller passes ``allow_proxy_variance_model``
    and records that it did. For a difference, prefer an EMPIRICAL null -- re-measure the same
    difference with the labels shuffled, on the same splits -- which needs no variance model and
    no design effect.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

__all__ = [
    "PreregistrationRefused",
    "IndependenceUnit",
    "Statistic",
    "Bar",
    "BarVerdict",
    "clean_failure_rate",
    "derive_bar",
    "validate_bar",
    "assert_not_a_function_of_the_design",
    "assert_non_vacuous",
    "assert_commensurable",
    "assert_range_reported_with_aggregate",
]

#: Below this many independent units the normal approximation to Spearman under H0 is not
#: trustworthy, and ``derive_bar`` refuses instead of extrapolating.
MIN_UNITS = 20

#: A bar may fail a clean pipeline this often and no more. 0.05 is itself a convention and is
#: recorded as one rather than derived -- but it is a convention about the BAR's own error rate,
#: which is a far weaker commitment than asserting a threshold on the data.
DEFAULT_MAX_CLEAN_FAILURE = 0.05

_HAZARDS = ("upward", "downward", "two-sided")
_AGGREGATIONS = ("single", "max-of-k", "median-of-k", "min-of-k")

#: Statistic kinds this module can actually simulate a null for, and the model it uses. A kind
#: absent from here has NO variance model, and ``derive_bar`` refuses it rather than silently
#: borrowing Spearman's. See gap G2 in the module docstring: a DIFFERENCE of two correlations on
#: one test set is not one correlation, because the two are dependent.
VARIANCE_MODELS = {
    "spearman": "asymptotic normal, SD = 1/sqrt(n_eff - 1)",
    "pearson": "asymptotic normal, SD = 1/sqrt(n_eff - 1)",
}


class PreregistrationRefused(Exception):
    """A bar or a claim was refused. The message names the check and what would discharge it."""


@dataclass(frozen=True)
class IndependenceUnit:
    """What one independent observation IS, and how many there really are.

    ``design_effect`` has NO default on purpose. Defaulting it to 1.0 is exactly F3: that
    assumption silently made a 3,822-row test set look six times more powerful than it is. A
    caller whose rows genuinely are independent passes ``design_effect=1.0`` and says so in
    ``measured_from``, which makes the assumption auditable instead of invisible.
    """

    name: str
    n_rows: int
    n_units: int
    design_effect: float
    measured_from: str

    def __post_init__(self) -> None:
        if not self.name:
            raise PreregistrationRefused("the independence unit needs a NAME (what one row is)")
        if self.n_rows <= 0 or self.n_units <= 0:
            raise PreregistrationRefused("n_rows and n_units must both be positive")
        if self.n_units > self.n_rows:
            raise PreregistrationRefused(
                f"n_units ({self.n_units}) exceeds n_rows ({self.n_rows}); the unit must be at "
                f"least as coarse as the row"
            )
        if self.design_effect < 1.0:
            raise PreregistrationRefused(
                f"design_effect {self.design_effect} < 1.0. Clustering cannot ADD independent "
                f"information, so a value below 1 means the estimate is wrong, not that the "
                f"design is efficient"
            )
        if not self.measured_from:
            raise PreregistrationRefused(
                "design_effect must say where it was MEASURED (artifact path or command). An "
                "unmeasured design effect is the failure this field exists to prevent; if the "
                "rows really are independent, pass design_effect=1.0 and record that here"
            )

    @property
    def n_effective(self) -> float:
        """Rows divided by the design effect -- the n that actually governs sampling noise."""
        return self.n_rows / self.design_effect

    @property
    def se_spearman_null(self) -> float:
        """SE of a Spearman correlation under H0 at the EFFECTIVE n."""
        if self.n_effective <= 1:
            raise PreregistrationRefused("effective n must exceed 1")
        return 1.0 / math.sqrt(self.n_effective - 1.0)

    @classmethod
    def from_icc(
        cls, name: str, *, n_rows: int, n_units: int, icc: float, measured_from: str
    ) -> "IndependenceUnit":
        """Build from a measured intra-cluster correlation, which is the usual route."""
        if not 0.0 <= icc <= 1.0:
            raise PreregistrationRefused(f"icc {icc} outside [0, 1]")
        m = n_rows / n_units
        return cls(
            name=name,
            n_rows=n_rows,
            n_units=n_units,
            design_effect=1.0 + (m - 1.0) * icc,
            measured_from=measured_from,
        )


@dataclass(frozen=True)
class Statistic:
    """The quantity a bar bounds, including how it is AGGREGATED and which way the hazard runs.

    ``aggregation``/``k`` exist because of F4: bounding a max-of-10 with a single-draw bar is a
    different test from the one intended. ``hazard`` exists because of F5.
    """

    name: str
    kind: str
    aggregation: str
    k: int
    hazard: str

    def __post_init__(self) -> None:
        if self.aggregation not in _AGGREGATIONS:
            raise PreregistrationRefused(
                f"aggregation {self.aggregation!r} not one of {_AGGREGATIONS}"
            )
        if self.hazard not in _HAZARDS:
            raise PreregistrationRefused(f"hazard {self.hazard!r} not one of {_HAZARDS}")
        if self.aggregation == "single" and self.k != 1:
            raise PreregistrationRefused("aggregation 'single' requires k == 1")
        if self.aggregation != "single" and self.k < 2:
            raise PreregistrationRefused(f"aggregation {self.aggregation!r} requires k >= 2")

    @property
    def two_sided(self) -> bool:
        return self.hazard == "two-sided"


@dataclass(frozen=True)
class Bar:
    """A threshold plus where it came from. ``inherited_from`` is what makes F2 visible."""

    value: float
    derived_from: str
    inherited_from: str | None = None


@dataclass(frozen=True)
class BarVerdict:
    ok: bool
    clean_failure_rate: float
    derived_value: float
    declared_value: float
    reasons: tuple[str, ...]

    def __str__(self) -> str:  # pragma: no cover - presentation only
        head = "OK" if self.ok else "REFUSED"
        body = "".join(f"\n  - {r}" for r in self.reasons) if self.reasons else ""
        return (
            f"{head}  declared={self.declared_value:.4f} derived={self.derived_value:.4f} "
            f"clean_failure_rate={self.clean_failure_rate:.3f}{body}"
        )


def _require_variance_model(stat: Statistic, allow_proxy: bool) -> None:
    """Refuse a kind this module cannot simulate (gap G2), unless the caller opts in explicitly."""
    if stat.kind in VARIANCE_MODELS or allow_proxy:
        return
    raise PreregistrationRefused(
        f"no variance model for kind {stat.kind!r}; this module can simulate only "
        f"{sorted(VARIANCE_MODELS)}. Borrowing Spearman's SD for a DIFFERENCE of two correlations "
        f"measured on one test set is wrong, because the two are dependent. Prefer an EMPIRICAL "
        f"null -- re-measure the same difference with the labels shuffled, on the same splits, "
        f"which needs no variance model and no design effect. To proceed anyway with Spearman's SD "
        f"as an acknowledged PROXY, pass allow_proxy_variance_model=True and record that you did"
    )


def _simulate(stat: Statistic, unit: IndependenceUnit, *, trials: int, seed: int) -> np.ndarray:
    """Sampling distribution of the AGGREGATE statistic under H0."""
    rng = np.random.default_rng(seed)
    draws = rng.normal(0.0, unit.se_spearman_null, size=(trials, stat.k))
    if stat.two_sided:
        signed = np.abs(draws)
    elif stat.hazard == "upward":
        signed = draws
    else:
        signed = -draws
    if stat.aggregation == "single":
        return signed[:, 0]
    if stat.aggregation == "max-of-k":
        return signed.max(axis=1)
    if stat.aggregation == "min-of-k":
        return signed.min(axis=1)
    return np.median(signed, axis=1)


def clean_failure_rate(
    stat: Statistic,
    unit: IndependenceUnit,
    bar_value: float,
    *,
    trials: int = 40_000,
    seed: int = 1117,
    allow_proxy_variance_model: bool = False,
) -> float:
    """P(a perfectly clean pipeline FAILS this bar).

    The number that separates a test from a formality: near 0.90 the verdict says nothing about
    the hazard it was written to detect.
    """
    _require_variance_model(stat, allow_proxy_variance_model)
    return float((_simulate(stat, unit, trials=trials, seed=seed) >= bar_value).mean())


def derive_bar(
    stat: Statistic,
    unit: IndependenceUnit,
    *,
    max_clean_failure: float = DEFAULT_MAX_CLEAN_FAILURE,
    trials: int = 40_000,
    seed: int = 1117,
    allow_proxy_variance_model: bool = False,
) -> float:
    """The smallest bar a clean pipeline clears with probability >= 1 - max_clean_failure.

    Derived from the DESIGN -- effective n, aggregation, sidedness -- and never from the observed
    result. That is what distinguishes this from bar-shopping.
    """
    _require_variance_model(stat, allow_proxy_variance_model)
    if unit.n_units < MIN_UNITS:
        raise PreregistrationRefused(
            f"{unit.n_units} independent units is below MIN_UNITS={MIN_UNITS}; the normal "
            f"approximation to Spearman under H0 is not trustworthy there. Use an exact "
            f"permutation null on the real grouping instead of this derivation"
        )
    if not 0.0 < max_clean_failure < 1.0:
        raise PreregistrationRefused("max_clean_failure must lie in (0, 1)")
    sim = _simulate(stat, unit, trials=trials, seed=seed)
    return float(np.quantile(sim, 1.0 - max_clean_failure))


def validate_bar(
    stat: Statistic,
    unit: IndependenceUnit,
    bar: Bar,
    *,
    max_clean_failure: float = DEFAULT_MAX_CLEAN_FAILURE,
    trials: int = 40_000,
    seed: int = 1117,
    allow_proxy_variance_model: bool = False,
) -> BarVerdict:
    """Validate a DECLARED bar against the design.

    Returns a verdict rather than raising, because the reasons are the product.
    """
    derived = derive_bar(
        stat, unit, max_clean_failure=max_clean_failure, trials=trials, seed=seed,
        allow_proxy_variance_model=allow_proxy_variance_model,
    )
    cfr = clean_failure_rate(
        stat, unit, bar.value, trials=trials, seed=seed,
        allow_proxy_variance_model=allow_proxy_variance_model,
    )
    reasons: list[str] = []

    if cfr > max_clean_failure:
        reasons.append(
            f"CLEAN-FAILURE RATE {cfr:.3f} exceeds {max_clean_failure:.3f}: a clean pipeline fails "
            f"this bar {cfr:.0%} of the time, so its verdict carries little information about the "
            f"hazard. The design-derived bar is {derived:.4f} against a declared {bar.value:.4f}"
        )
    if bar.inherited_from:
        reasons.append(
            f"INHERITED from {bar.inherited_from!r}. A bar transfers between designs only when the "
            f"effective n, the aggregation and the sidedness all match; state why they do, or "
            f"re-derive it at this design (effective n {unit.n_effective:.0f})"
        )
    if not bar.derived_from:
        reasons.append("no `derived_from`: a bar must record what produced it")
    if stat.two_sided:
        reasons.append(
            "HAZARD is 'two-sided' -- confirm that is deliberate. A leakage or contamination "
            "hazard runs in ONE direction, and taking abs() turns a one-sided test into a "
            "two-sided one, which is what made three negative null draws read as a failure"
        )
    return BarVerdict(
        ok=not reasons,
        clean_failure_rate=cfr,
        derived_value=derived,
        declared_value=bar.value,
        reasons=tuple(reasons),
    )


def assert_not_a_function_of_the_design(
    stat_fn: Callable[[object], float], datasets: Sequence[object], *, tol: float = 1e-9
) -> None:
    """Refuse a statistic that returns the same value on same-design, different-data inputs.

    F6: rho = sqrt(3k(n-k)/(n^2-1)) was published before anyone checked that it does not read the
    data at all. Pass two or more datasets that share the design and differ in their values.
    """
    if len(datasets) < 2:
        raise PreregistrationRefused(
            "need at least 2 same-design / different-data datasets; one cannot show variation"
        )
    vals = [float(stat_fn(d)) for d in datasets]
    if max(vals) - min(vals) <= tol:
        raise PreregistrationRefused(
            f"the statistic returned {vals[0]!r} on all {len(vals)} datasets, which share a design "
            f"and differ in their data. It is a function of the DESIGN, not a measurement of the "
            f"data -- the closed-form-in-the-split-proportion failure"
        )


def assert_non_vacuous(n_before: int, n_after: int, *, what: str) -> None:
    """Refuse a filter or control that removed nothing (F7)."""
    if n_before == n_after:
        raise PreregistrationRefused(
            f"{what} removed 0 of {n_before} rows. A filter that removes nothing is not a control; "
            f"either drop the claim that it is one, or state explicitly that its removal count is "
            f"zero on this data and that it guards a different input shape"
        )


def assert_commensurable(named_values: dict[str, str]) -> None:
    """Refuse a comparison that mixes metric kinds (F9), such as a fraction beside a Spearman."""
    kinds = set(named_values.values())
    if len(kinds) > 1:
        detail = ", ".join(f"{n}={k}" for n, k in sorted(named_values.items()))
        raise PreregistrationRefused(
            f"metrics are not commensurable: {detail}. A max() or a ranking over mixed kinds "
            f"silently promotes whichever happens to be numerically largest"
        )


def assert_range_reported_with_aggregate(
    aggregate: float, strata: Sequence[float], *, what: str, spread_tol: float = 0.05
) -> None:
    """Refuse quoting a stratified aggregate as a scalar when its strata disagree (F8).

    0.941 is a pooled average over deciles whose own ceilings run 0.5013 to 0.9823. Quoting the
    average alone invites comparing a score against the mean of a range it sits inside.
    """
    if len(strata) < 2:
        raise PreregistrationRefused(
            f"{what}: need the per-stratum values to judge whether the aggregate stands alone"
        )
    lo, hi = min(strata), max(strata)
    if hi - lo > spread_tol:
        raise PreregistrationRefused(
            f"{what}: the aggregate {aggregate:.4f} spans strata from {lo:.4f} to {hi:.4f} "
            f"(spread {hi - lo:.4f} > {spread_tol}). Report the RANGE with it; a scalar invites "
            f"comparing a score against the mean of a range it sits inside"
        )
