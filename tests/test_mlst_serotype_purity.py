"""The MLST caller's sequence types are serotype-pure well beyond a shape-matched shuffle null.

The design point these tests protect: purity ALONE is gameable in both directions (one ST per genome
scores 1.000; one ST for everything scores the modal frequency), so the null must hold the observed
partition FIXED and shuffle only the labels. A test suite that pinned the observed number without
pinning that property would let the null be replaced by a weaker one and never notice.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from mlst_serotype_purity import (  # noqa: E402
    MIN_ST_SIZE, serotype_of, shuffle_null, weighted_purity,
)

ARTIFACT = ROOT / "wiki" / "mlst_serotype_purity_2026-09-05.json"


# --- the purity statistic -------------------------------------------------------------------------

def test_a_perfectly_pure_partition_scores_one():
    groups = {"11": ["O157:H7"] * 5, "21": ["O26:H11"] * 4}
    assert weighted_purity(groups)[0] == pytest.approx(1.0)


def test_singleton_sts_are_excluded_because_they_are_pure_by_construction():
    """An ST of size 1 is 100% pure no matter what the caller did; counting it inflates the result."""
    groups = {"a": ["O1:H1"], "b": ["O2:H2"], "c": ["O3:H3"]}
    purity, n_scored, n_st = weighted_purity(groups)
    assert (n_scored, n_st) == (0, 0)
    assert MIN_ST_SIZE >= 2


def test_top2_captures_a_bimodal_lineage_that_top1_penalizes():
    """A real lineage can carry two serotypes; top-1 cannot distinguish that from mis-grouping."""
    groups = {"131": ["O25:H4"] * 11 + ["O16:H5"] * 9}
    assert weighted_purity(groups, top_k=1)[0] == pytest.approx(11 / 20)
    assert weighted_purity(groups, top_k=2)[0] == pytest.approx(1.0)


def test_purity_is_weighted_by_genomes_not_by_st():
    """A big impure ST must not be outvoted by several tiny pure ones."""
    groups = {"big": ["O1:H1"] * 5 + ["O2:H2"] * 5, "s1": ["O3:H3"] * 3, "s2": ["O4:H4"] * 3}
    assert weighted_purity(groups)[0] == pytest.approx((5 + 3 + 3) / 16)


# --- the null -------------------------------------------------------------------------------------

def test_the_null_preserves_the_partition_shape():
    """Group sizes must be identical under the null, or it is not a matched comparison."""
    groups = {"a": ["O1:H1"] * 5, "b": ["O2:H2"] * 3, "c": ["O3:H3"] * 4}
    sizes = sorted(len(v) for v in groups.values())
    import mlst_serotype_purity as m
    captured = []
    real = m.weighted_purity

    def spy(g, top_k=1):
        captured.append(sorted(len(v) for v in g.values()))
        return real(g, top_k=top_k)

    m.weighted_purity = spy
    try:
        shuffle_null(groups, 5)
    finally:
        m.weighted_purity = real
    assert captured and all(c == sizes for c in captured)


def test_the_null_is_high_when_the_cohort_is_serotype_poor():
    """Sanity on the null itself: with few distinct serotypes, a random partition is already pure,
    which is exactly the confound the shuffle is there to absorb."""
    groups = {"a": ["O1:H1"] * 5, "b": ["O1:H1"] * 5}
    assert min(shuffle_null(groups, 20)) == pytest.approx(1.0)


def test_an_uninformative_partition_does_not_beat_its_own_null():
    """The decisive negative control: labels spread evenly across STs carry no lineage signal."""
    groups = {str(i): ["O1:H1", "O2:H2", "O3:H3"] for i in range(8)}
    obs = weighted_purity(groups)[0]
    assert obs <= max(shuffle_null(groups, 200))


def test_the_null_is_deterministic():
    groups = {"a": ["O1:H1"] * 4 + ["O2:H2"] * 2, "b": ["O3:H3"] * 5}
    assert shuffle_null(groups, 20) == shuffle_null(groups, 20)


# --- parsing --------------------------------------------------------------------------------------

def test_serotype_requires_both_axes():
    assert serotype_of({"O": "O157", "H": "H7"}) == "O157:H7"
    assert serotype_of({"O": "O157", "H": ""}) is None
    assert serotype_of({"O": "", "H": "H7"}) is None


# --- the committed artifact -----------------------------------------------------------------------

@pytest.fixture(scope="module")
def art() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("mlst purity artifact absent")
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_artifact_verdict(art):
    assert art["verdict"] == "MLST_RECOVERS_REAL_LINEAGES"


def test_every_statistic_beats_its_OWN_null_maximum(art):
    """Not just the mean, and not against a borrowed null."""
    assert art["observed_purity"] > art["null_max"]
    for m in art["additional_statistics"]:
        assert m["exceeds_null_max"], m["statistic"]
        assert m["observed"] > m["null_max"], m["statistic"]


def test_artifact_is_non_vacuous(art):
    """A cohort with one serotype, or no qualifying ST, would make every number meaningless."""
    assert art["n_distinct_serotypes"] > 1
    assert art["n_st_scored"] > 0 and art["n_genomes_scored"] > 0


def test_the_low_purity_sts_are_recorded_not_explained_away(art):
    """Their composition ships, so a reader can check the biology claim rather than take it."""
    comp = art["st_composition_for_low_purity_sts"]
    assert comp, "expected at least one ST below the 0.70 purity line"
    assert "131" in comp and len(comp["131"]) >= 2


def test_the_tier_does_not_move_and_the_limit_says_why(art):
    joined = " ".join(art["honest_limits"]).lower()
    assert "coherence check" in joined
    assert "faithful_to_tool" in joined
