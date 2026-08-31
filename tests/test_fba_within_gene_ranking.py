"""Pin the within-gene switch metric, and the four tie-handling defects found while building it.

Every one of these was a real bug in the first working version, and every one of them made the result
look BETTER than it was. They share a root: this data is full of exact ties (a flat gene emits one growth
ratio for all four media), and at every tie some arbitrary choice can masquerade as a result.

Offline; no cobrapy solve, no network.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from fba_within_gene_ranking import FLAT_EPS, score, within_gene_auroc  # noqa: E402


class _Rec:
    """Minimal GeneRecord stand-in: the scorer only reads `gene_id` and `experimental`."""

    def __init__(self, gene_id, experimental):
        self.gene_id = gene_id
        self.experimental = experimental


def _conds():
    from dna_decode.fba.conditional_essentiality import CONDITIONS

    return sorted(CONDITIONS)


def _mk(ess_by_cond, ratio_by_cond, gene_id="g1"):
    keys = _conds()
    rec = _Rec(gene_id, {c: e for c, e in zip(keys, ess_by_cond)})
    ratios = {c: {gene_id: r} for c, r in zip(keys, ratio_by_cond)}
    return ratios, [rec]


def test_perfect_within_gene_ranking_scores_one():
    assert within_gene_auroc([True, False, False, False], [0.0, 0.9, 0.9, 1.0]) == 1.0


def test_inverted_ranking_scores_zero():
    assert within_gene_auroc([True, False, False, False], [1.0, 0.1, 0.2, 0.3]) == 0.0


def test_one_sided_gene_is_undefined_not_a_half():
    """All-essential or all-dispensable: the switch question has no meaning. None, never 0.5."""
    assert within_gene_auroc([True] * 4, [0.1, 0.2, 0.3, 0.4]) is None
    assert within_gene_auroc([False] * 4, [0.1, 0.2, 0.3, 0.4]) is None


def test_flat_gene_scores_exactly_half_via_the_tie_tolerance():
    """DEFECT 1. Flatness used a 1e-9 tolerance while the tie test used exact float equality, so a gene
    whose ratios differ in the 15th decimal was scored as if its ordering meant something -- 36 of 41
    flat genes returned something other than 0.5, purely from LP noise."""
    assert within_gene_auroc([True, False, True, False], [0.5] * 4) == 0.5
    noise = [0.5, 0.5 + 1e-15, 0.5 - 1e-15, 0.5 + 2e-16]
    assert within_gene_auroc([True, False, True, False], noise) == 0.5


def test_flat_gene_is_never_an_oracle_hit():
    """DEFECT 2. With every ratio equal, `sorted` falls back to index order, so 'the k lowest conditions'
    means 'the first k conditions'. 12 of 41 flat genes matched their true pattern on condition ordering
    alone -- inflating the reported ceiling from 11 to 23."""
    keys = _conds()
    ess = [True, False, False, False]          # k=1, and the essential one is FIRST in sort order
    ratios, recs = _mk(ess, [0.7, 0.7, 0.7, 0.7])
    out = score(ratios, recs)
    assert out["n_flat"] == 1
    assert out["oracle_relative_exact_set"] == 0, "a flat gene has no ordering and cannot be a hit"
    assert keys  # the fixture really did use the live condition keys


def test_ambiguous_top_k_boundary_is_never_an_oracle_hit():
    """DEFECT 3. The partial version: when the k-th and (k+1)-th ratios tie, which condition lands inside
    the selection is decided by stable-sort order. Two genes at AUROC 0.833 were hitting that way."""
    # k=1; the essential condition ties with a dispensable one at the selection boundary.
    ratios, recs = _mk([True, False, False, False], [0.2, 0.2, 0.9, 1.0])
    out = score(ratios, recs)
    assert out["n_flat"] == 0, "this gene genuinely varies -- it is not the flat case"
    assert out["oracle_relative_exact_set"] == 0


def test_oracle_ceiling_equals_the_count_of_perfectly_ranked_genes():
    """The invariant the three fixes above enforce: top-k selection is exact precisely when every
    essential condition ranks strictly below every dispensable one, i.e. AUROC == 1.0."""
    keys = _conds()
    recs, ratios = [], {c: {} for c in keys}
    specs = [("perfect", [True, False, False, False], [0.0, 0.6, 0.7, 0.8]),
             ("partial", [True, True, False, False], [0.1, 0.7, 0.4, 0.9]),
             ("flat", [True, False, False, False], [0.5, 0.5, 0.5, 0.5]),
             ("inverted", [True, False, False, False], [0.9, 0.1, 0.2, 0.3])]
    for gid, ess, rat in specs:
        recs.append(_Rec(gid, {c: e for c, e in zip(keys, ess)}))
        for c, r in zip(keys, rat):
            ratios[c][gid] = r
    out = score(ratios, recs)
    n_perfect = sum(1 for g in out["per_gene"] if g["auroc"] == 1.0)
    assert out["oracle_relative_exact_set"] == n_perfect == 1


def test_flat_genes_are_reported_not_silently_dropped():
    """A flat gene is the FINDING (the model emits one number for four media), so it must stay in the
    denominator. Dropping it would inflate the headline into 'the model switches well'."""
    keys = _conds()
    recs, ratios = [], {c: {} for c in keys}
    for gid, rat in (("flat", [0.5] * 4), ("varying", [0.0, 0.6, 0.7, 0.8])):
        recs.append(_Rec(gid, {c: e for c, e in zip(keys, [True, False, False, False])}))
        for c, r in zip(keys, rat):
            ratios[c][gid] = r
    out = score(ratios, recs)
    assert out["n_genes_scored"] == 2 and out["n_flat"] == 1 and out["n_nonflat"] == 1
    assert out["flat_fraction"] == 0.5
    # the all-genes mean is dragged toward chance by the flat gene, exactly as intended
    assert out["mean_auroc_all"] == 0.75 and out["mean_auroc_nonflat"] == 1.0


def test_flat_eps_is_shared_by_both_tests():
    """DEFECT 4 was two tolerances disagreeing. Keep one constant; this pins that it is used, not shadowed."""
    assert FLAT_EPS > 0
    just_inside = [0.5, 0.5 + FLAT_EPS / 2, 0.5, 0.5]
    assert within_gene_auroc([True, False, True, False], just_inside) == 0.5


# --- the deployed-comparator lookup: an axis must be compared against its OWN baseline ---

def test_each_axis_reads_its_own_deployed_exact_set():
    """Quoting the 4-media 3/67 beside a 25-condition carbon result would compare a ceiling against the
    wrong baseline and make the ranking lever look ~4x more valuable than it is. Each axis reads its own
    committed artifact."""
    from fba_within_gene_ranking import deployed_exact_set

    m4, carbon = deployed_exact_set("media4"), deployed_exact_set("carbon")
    for got in (m4, carbon):
        if got is None:
            continue                                    # artifact not committed on this checkout
        assert got["match"] <= got["n"] and got["n"] > 0
    if m4 and carbon:
        assert (m4["match"], m4["n"]) != (carbon["match"], carbon["n"]), \
            "the two axes must not resolve to the same baseline"


def test_deployed_lookup_tolerates_both_artifact_schemas():
    """Schema drift between two generations of the same producer: the carbon artifact carries
    `n_scored_exact_set` (added with `exclude_cells`), the older 4-media one only
    `n_conditionally_essential`. Reading one key alone silently reported 'unknown'."""
    from fba_within_gene_ranking import deployed_exact_set

    for axis, expect_n in (("media4", 67), ("carbon", 217)):
        got = deployed_exact_set(axis)
        if got is not None:
            assert got["n"] == expect_n, f"{axis} resolved n={got['n']}, expected {expect_n}"


def test_unknown_axis_returns_none_rather_than_a_wrong_baseline():
    from fba_within_gene_ranking import deployed_exact_set

    assert deployed_exact_set("nitrogen_not_wired_yet") is None


# --- the three-stratum anatomy: a constant call can never match a two-sided truth ---

def _score_one(ess, rat):
    ratios, recs = _mk(ess, rat)
    return score(ratios, recs)


def test_flat_gene_lands_in_the_flat_stratum_and_never_hits():
    out = _score_one([True, False, False, False], [0.5, 0.5, 0.5, 0.5])
    assert out["strata"]["flat"]["n_genes"] == 1
    assert out["deployed_exact_set_recomputed"] == 0


def test_varying_ratio_that_never_crosses_the_cutoff_is_its_own_stratum():
    """The stratum a ranking rule could rescue: the ordering moves, the CALL does not."""
    out = _score_one([True, False, False, False], [0.30, 0.60, 0.70, 0.80])
    assert out["strata"]["varies_subthr"]["n_genes"] == 1
    assert out["strata"]["commits"]["n_genes"] == 0
    assert out["deployed_exact_set_recomputed"] == 0, "nothing crosses 1%, so the call is constant"


def test_a_committing_gene_that_is_right_is_the_only_way_to_score():
    out = _score_one([True, False, False, False], [0.001, 0.60, 0.70, 0.80])
    assert out["strata"]["commits"]["n_genes"] == 1
    assert out["strata"]["commits"]["deployed_exact"] == 1


def test_strata_partition_every_scored_gene_exactly_once():
    """If these stop summing, a gene is being double-counted or dropped and every rate is wrong."""
    keys = _conds()
    recs, ratios = [], {c: {} for c in keys}
    specs = [("flat", [0.5] * 4), ("subthr", [0.3, 0.6, 0.7, 0.8]), ("commits", [0.001, 0.6, 0.7, 0.8])]
    for gid, rat in specs:
        recs.append(_Rec(gid, {c: e for c, e in zip(keys, [True, False, False, False])}))
        for c, r in zip(keys, rat):
            ratios[c][gid] = r
    out = score(ratios, recs)
    total = sum(out["strata"][k]["n_genes"] for k in ("flat", "varies_subthr", "commits"))
    assert total == out["n_genes_scored"] == 3


def test_chance_null_ignores_genes_whose_predicted_count_differs_from_truth():
    """A gene predicting 2 essential conditions against a truth of 1 cannot match at ANY placement, so it
    contributes 0 to the chance expectation -- part of why the observed hit count is hard to reach."""
    out = _score_one([True, False, False, False], [0.001, 0.002, 0.7, 0.8])   # predicts 2, truth is 1
    assert out["strata"]["commits"]["n_genes"] == 1
    assert out["chance_exact_hits_among_committing"] == 0.0


def test_chance_null_is_one_over_n_choose_k_when_the_counts_agree():
    out = _score_one([True, False, False, False], [0.001, 0.6, 0.7, 0.8])     # predicts 1, truth is 1
    assert out["chance_exact_hits_among_committing"] == round(1 / 4, 4)


# --- the margin-preserving null: conditions are NOT interchangeable ---

def test_marginal_null_uses_the_condition_keyed_orientation():
    """The nulls module consumes {condition: {gene: bool}} -- the orientation every FBA caller uses.
    Building it gene-keyed raised KeyError on the first shuffle. This pins the wiring end-to-end."""
    out = _score_one([True, False, False, False], [0.001, 0.6, 0.7, 0.8])
    mp = out["marginal_preserving_null"]
    assert mp is not None and mp["n_draws"] > 0
    assert mp["mean"] is not None


def test_marginal_null_is_absent_when_nothing_commits():
    """With no committing gene there is no exact-set question, so the null is None rather than a
    meaningless 0.0 that would read as 'the model beat chance'."""
    out = _score_one([True, False, False, False], [0.5, 0.5, 0.5, 0.5])   # flat -> constant call
    assert out["strata"]["commits"]["n_genes"] == 0
    assert out["marginal_preserving_null"] is None


def test_the_two_nulls_are_reported_separately():
    """They answer different questions -- interchangeable-conditions vs both-margins-preserved -- and the
    stricter one changed a verdict (media4's 3 observed sits at the margin-preserving null's own max)."""
    out = _score_one([True, False, False, False], [0.001, 0.6, 0.7, 0.8])
    assert "chance_exact_hits_among_committing" in out
    assert "marginal_preserving_null" in out
