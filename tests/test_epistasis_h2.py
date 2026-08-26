"""The H2 within-protein test + its negative controls (scripts/epistasis_pooling_h2_test.py).

H2: the additive score's pooling gain is governed by eta^2(k) -- the share of fitness variance sitting
BETWEEN mutation orders. It was left "consistent but underpowered" (n=3 proteins) because eta^2 was
perfectly CONFOUNDED with protein identity. These pin the pure math and the two interventions that turn
the correlation into a causal claim.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.epistasis_pooling_h2_test import (  # noqa: E402
    additive_scores, controls, eta2, subset_stats,
)


def test_additive_score_sums_log_ratios_and_refuses_unparseable_tokens():
    table = {"1": {"A": -1.0, "M": 0.0}, "2": {"C": -3.0, "G": -1.0}}
    assert additive_scores(["M1A"], table) == [-1.0]              # -1.0 - 0.0
    assert additive_scores(["M1A:G2C"], table) == [-3.0]          # -1.0 + (-3.0 - -1.0)
    # a position outside the table, a malformed token, and an unknown residue must all yield None --
    # never a silently-truncated partial sum
    assert additive_scores(["M9A"], table) == [None]
    assert additive_scores(["notamutation"], table) == [None]
    assert additive_scores(["M1Z"], table) == [None]


def test_eta2_is_zero_when_groups_share_a_mean_and_one_when_they_do_not_overlap():
    assert eta2([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]) == pytest.approx(0.0, abs=1e-9)
    assert eta2([[5.0, 5.0], [9.0, 9.0]]) == pytest.approx(1.0)   # all variance is between groups


def test_pooling_gain_appears_only_when_groups_separate_in_BOTH_score_and_fitness():
    """The mechanism: pooling manufactures correlation when the groups differ on both axes."""
    # groups separated in fitness AND score -> pooling invents a positive rho the subgroups lack
    sep = {2: [(0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0)],
           3: [(10.0, 10.0), (11.0, 11.0), (10.0, 11.0), (11.0, 10.0)]}
    st = subset_stats(sep, (2, 3))
    assert st["eta2"] > 0.8 and st["pooling_gain"] > 0.5

    # same fitness distribution in both groups -> nothing to harvest
    flat = {2: [(0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0)],
            3: [(0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0)]}
    st2 = subset_stats(flat, (2, 3))
    assert st2["eta2"] == pytest.approx(0.0, abs=1e-9)
    assert abs(st2["pooling_gain"]) < 1e-9


def test_the_two_controls_isolate_between_order_structure():
    """A (shuffle order labels) must KILL the gain; B (shuffle fitness within order) must NOT.

    This is what makes H2 causal rather than correlational -- and B is the sharper of the two: it proves
    the gain was never made of within-order ranking skill in the first place.
    """
    # The fixture must leave HEADROOM: score and fitness are UNCORRELATED within each group (so
    # within_rho ~ 0) while the groups are separated on BOTH axes. A first version used score == fitness
    # exactly, which pins within_rho at +1.000 and makes the real gain 0 by construction -- the fixture
    # then tests nothing. Real assays sit near this shape, not that one.
    def _group(offset: float, n: int = 40):
        # a fixed permutation decorrelates score from fitness inside the group
        return [(offset + i, offset + ((i * 17) % n)) for i in range(n)]

    by_order = {2: _group(0.0), 3: _group(100.0)}
    c = controls(by_order, seed=0)
    assert c["real"]["pooling_gain"] > 0.1
    assert abs(c["ctl_a_order_labels_shuffled"]["pooling_gain"]) < 0.05      # collapses
    assert c["ctl_a_order_labels_shuffled"]["eta2"] == pytest.approx(0.0, abs=0.05)
    assert c["ctl_b_fitness_shuffled_within_order"]["pooling_gain"] > 0.1    # survives
    assert c["ctl_b_fitness_shuffled_within_order"]["eta2"] == pytest.approx(c["real"]["eta2"], abs=1e-6)


def test_the_committed_h2_artifact_shows_the_confirmed_pattern():
    p = ROOT / "wiki" / "forward_epistasis_h2_within_protein_2026-08-25.json"
    if not p.exists():
        pytest.skip("H2 artifact absent (needs the D: assay + ESM caches to regenerate)")
    d = json.loads(p.read_text(encoding="utf-8"))
    powered = {k: v for k, v in d["proteins"].items() if v["n_subsets"] >= 20}
    assert len(powered) >= 2, "expected 2 proteins with enough order-subsets to be powered"
    for name, v in powered.items():
        assert v["spearman_eta2_vs_gain"] > 0.9 and v["p"] < 1e-10, name
    # every protein, powered or not, must show the intervention pattern
    for name, v in d["proteins"].items():
        c = v["controls"]
        assert abs(c["ctl_a_order_labels_shuffled"]["pooling_gain"]) < 0.02, f"{name}: ctlA did not kill it"
        assert c["ctl_b_fitness_shuffled_within_order"]["pooling_gain"] > 0.5 * c["real"]["pooling_gain"], \
            f"{name}: ctlB should have left the gain intact"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
