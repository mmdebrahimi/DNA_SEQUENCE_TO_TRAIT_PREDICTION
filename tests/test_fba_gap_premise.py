"""Track C premise check -- pure logic (wheel-only). The cobra/network path is exercised by the script."""
from __future__ import annotations

import pytest

from scripts.fba_gap_premise_check import fisher_exact_two_sided


def test_fisher_recovers_the_real_yeast_enrichment():
    """The observed table: FN 65/92 gap-adjacent vs TN 307/703. This is the number the verdict rests on."""
    p = fisher_exact_two_sided(65, 92 - 65, 307, 703 - 307)
    assert p == pytest.approx(1.224e-06, rel=0.05)
    assert p < 0.05          # premise SUPPORTED


def test_fisher_returns_1_for_identical_rates():
    """Equal FN and TN gap-adjacency is the falsifying outcome and must not be called significant."""
    assert fisher_exact_two_sided(50, 50, 100, 100) == pytest.approx(1.0)


def test_fisher_is_two_sided_and_symmetric_under_row_swap():
    """Two-sided means a DEPLETION is as detectable as an enrichment -- if the FNs were LESS gap-adjacent
    than the TNs that would also be a finding (and would still falsify the premise as written)."""
    assert fisher_exact_two_sided(65, 27, 307, 396) == pytest.approx(
        fisher_exact_two_sided(307, 396, 65, 27))


def test_fisher_handles_degenerate_margins_without_raising():
    """An empty cell class (e.g. zero false negatives) must return nan, not blow up the whole run."""
    import math
    assert math.isnan(fisher_exact_two_sided(0, 0, 10, 10))
    assert math.isnan(fisher_exact_two_sided(0, 0, 0, 0))


def test_the_verdict_rule_requires_direction_AND_significance():
    """A significant result in the WRONG direction (FNs less gap-adjacent) must not read as support."""
    def supported(fn_rate, tn_rate, p):
        return bool(fn_rate > tn_rate and p < 0.05)

    assert supported(0.707, 0.437, 1.2e-06) is True
    assert supported(0.300, 0.437, 1.2e-06) is False     # significant, wrong direction
    assert supported(0.707, 0.437, 0.400) is False       # right direction, not significant
