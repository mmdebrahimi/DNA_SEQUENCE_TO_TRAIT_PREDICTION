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
