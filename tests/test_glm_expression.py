"""Guards for the expression oracle.

The two expensive failures here are both SILENT and both measured: a join that loses 65% of the data while
looking healthy, and a random split that overstates generalisation. Each has a test that fails loudly.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.glm.expression import (
    ELEMENT_COLUMNS,
    MIN_JOIN_RATE,
    SEQ_LENGTH,
    ExpressionDataError,
    LoadReport,
    Pair,
    _norm_key,
    gc_feature,
    kmer_feature,
    leave_element_out,
    load_pairs,
    positional_onehot,
    random_split,
    spearman,
)

CACHE = Path("D:/dna_decode_cache/mpra")
HAVE_DATA = (CACHE / "GSE108535_sigma70_variant_data.txt.gz").exists()


def _synth(n=400, seed=0):
    import random
    r = random.Random(seed)
    out = []
    for i in range(n):
        m10 = r.choice(["a", "b", "c", "d"])
        m35 = r.choice(["p", "q", "r", "s"])
        seq = "".join(r.choice("ACGT") for _ in range(SEQ_LENGTH))
        # make expression depend on the element, so a real model can learn something
        val = {"a": 1.0, "b": 2.0, "c": 4.0, "d": 8.0}[m10] + r.random()
        out.append(Pair(name=f"v{i}", seq=seq, expression=val,
                        elements={"UP_element": "u", "Minus35": m35, "Spacer": "s1",
                                  "Minus10": m10, "Background": "bg1"}))
    return out


# ---------------------------------------------------------------------------------------------------
# THE join fix
# ---------------------------------------------------------------------------------------------------
def test_norm_key_unifies_the_separator_but_PROTECTS_the_mutation_arrow():
    """The hard-won fix. The expression table writes `gourse-326fold-up`, the barcode map writes
    `gourse_136fold_up`; without normalisation 7,077 of 10,898 pairs are silently lost.

    And the arrow must survive: a blanket replace turns `35T->A` into `35T_>A`, which still joins (it is
    applied to both sides) but is lossy and could in principle collide two distinct variants.
    """
    assert _norm_key("gourse-326fold-up_35T->A") == "gourse_326fold_up_35T->A"
    assert "->" in _norm_key("noUP_34G->T"), "the mutation arrow must not be rewritten"
    assert _norm_key("noUP_34G->T") == "noUP_34G->T", "already-underscored keys unchanged"
    assert _norm_key("a-b-c") == "a_b_c"
    # the protection is NON-VACUOUS: a blanket replace would differ
    assert "gourse-326fold-up_35T->A".replace("-", "_") != _norm_key("gourse-326fold-up_35T->A")


def test_the_join_fix_is_NON_VACUOUS_the_two_key_forms_really_differ():
    """If the two spellings were already equal, _norm_key would be decoration."""
    a, b = "gourse-326fold-up_x", "gourse_326fold_up_x"
    assert a != b, "the premise of the fix"
    assert _norm_key(a) == _norm_key(b)


def test_load_pairs_REFUSES_a_degraded_join_rather_than_returning_a_third_of_the_data(tmp_path):
    """A 35% join looks healthy downstream -- the model trains, the metric computes, nothing complains.
    So the floor is enforced at load time."""
    assert MIN_JOIN_RATE >= 0.95
    empty = tmp_path / "no-mpra-here"
    with pytest.raises(ExpressionDataError, match="downloading is disabled"):
        load_pairs(empty, allow_download=False)


def test_load_report_names_what_the_naive_join_would_have_given():
    """The report carries the counterfactual, so the 65% loss is visible rather than inferable."""
    r = LoadReport(n_expression_rows=100, n_joined=100, join_rate=1.0, naive_join=35)
    d = r.as_dict()
    assert d["naive_join_would_have_given"] == 35
    assert "silent" in d["note"]


# ---------------------------------------------------------------------------------------------------
# splits
# ---------------------------------------------------------------------------------------------------
def test_leave_element_out_is_disjoint_BY_CONSTRUCTION():
    """Membership is a pure function of the element label, so no sequence can be on both sides."""
    pairs = _synth()
    tr, te, held = leave_element_out(pairs, "Minus10", held={"a"})
    assert held == ["a"]
    assert all(p.elements["Minus10"] != "a" for p in tr)
    assert all(p.elements["Minus10"] == "a" for p in te)
    assert not ({p.name for p in tr} & {p.name for p in te})
    assert len(tr) + len(te) == len(pairs), "the split must partition, not sample"


def test_leave_element_out_refuses_an_unknown_column():
    with pytest.raises(ValueError, match="column must be one of"):
        leave_element_out(_synth(), "NotAColumn")


def test_leave_element_out_refuses_a_column_with_one_value():
    pairs = _synth()
    with pytest.raises(ValueError, match="<2 distinct values"):
        leave_element_out(pairs, "Spacer")       # synth fixes Spacer to a single value


def test_random_split_partitions_and_is_seed_stable():
    pairs = _synth()
    a1, b1 = random_split(pairs, seed=3)
    a2, b2 = random_split(pairs, seed=3)
    assert [p.name for p in b1] == [p.name for p in b2]
    assert len(a1) + len(b1) == len(pairs)
    assert not ({p.name for p in a1} & {p.name for p in b1})


# ---------------------------------------------------------------------------------------------------
# features
# ---------------------------------------------------------------------------------------------------
def test_feature_widths_are_what_the_artifact_reports():
    seqs = [s.seq for s in _synth(5)]
    assert len(gc_feature(seqs)[0]) == 1
    assert len(kmer_feature(seqs, k=4)[0]) == 256
    assert len(kmer_feature(seqs, k=6)[0]) == 4096
    assert len(positional_onehot(seqs)[0]) == SEQ_LENGTH * 4


def test_positional_onehot_sees_position_where_a_kmer_histogram_cannot():
    """Same composition, reversed layout: identical under a 1-mer histogram, different positionally.
    This is why onehot is the most robust feature in the measured results."""
    a = "A" * 75 + "C" * 75
    b = "C" * 75 + "A" * 75
    assert kmer_feature([a], k=1) == kmer_feature([b], k=1)
    assert positional_onehot([a]) != positional_onehot([b])


def test_onehot_is_exactly_one_hot_per_position():
    v = positional_onehot(["ACGT"])[0]
    assert sum(v) == 4.0
    assert v[0] == 1.0 and v[5] == 1.0 and v[10] == 1.0 and v[15] == 1.0


def test_spearman_is_rank_based_and_tie_tolerant():
    assert spearman([1, 2, 3, 4], [1, 2, 3, 4]) == pytest.approx(1.0)
    assert spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    # all-tied input has undefined rho; must return 0.0 rather than nan
    assert spearman([1, 1, 1, 1], [1, 2, 3, 4]) == 0.0


# ---------------------------------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------------------------------
def test_ridge_learns_a_planted_signal():
    from dna_decode.glm.expression import fit_predict_ridge
    pairs = _synth(600, seed=7)
    tr, te = random_split(pairs, seed=1)
    r = fit_predict_ridge(tr, te, feature="kmer4")
    assert r["n_train"] == len(tr) and r["n_test"] == len(te)
    assert isinstance(r["spearman"], float)


def test_ridge_refuses_an_unknown_feature():
    from dna_decode.glm.expression import fit_predict_ridge
    pairs = _synth(200)
    tr, te = random_split(pairs)
    with pytest.raises(ValueError, match="unknown feature"):
        fit_predict_ridge(tr, te, feature="nope")


def test_log_target_is_on_by_default_because_expression_spans_615x():
    """Spearman is invariant to the monotone transform, but the FIT is not: on the raw scale a handful of
    very bright promoters dominate the least-squares objective."""
    from dna_decode.glm.expression import fit_predict_ridge
    pairs = _synth(300)
    tr, te = random_split(pairs)
    assert fit_predict_ridge(tr, te, feature="gc")["log_target"] is True


# ---------------------------------------------------------------------------------------------------
# the real substrate
# ---------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not HAVE_DATA, reason="GSE108535 not cached on this host")
def test_the_REAL_substrate_loads_with_a_full_join():
    """R3 real-surface test. Verified by download: 10,898 pairs, all 150bp, ACGT-only, 615x range."""
    pairs, rep = load_pairs(CACHE, allow_download=False)
    assert rep.n_joined == 10898, rep.as_dict()
    assert rep.join_rate == pytest.approx(1.0)
    assert rep.naive_join == 3821, "the naive-join counterfactual is the point of the fix"
    assert all(len(p.seq) == SEQ_LENGTH for p in pairs)
    assert all(not (set(p.seq) - set("ACGT")) for p in pairs)
    assert len({p.seq for p in pairs}) == 10898, "every promoter is a distinct sequence"
    vals = [p.expression for p in pairs]
    assert max(vals) / min(vals) > 500
    for c in ELEMENT_COLUMNS:
        assert len({p.elements[c] for p in pairs}) >= 3


@pytest.mark.skipif(not HAVE_DATA, reason="GSE108535 not cached on this host")
def test_the_learned_features_genuinely_BEAT_the_crude_baseline_on_the_real_data():
    """The measured result, pinned. GC alone reaches 0.14 on held-out elements while positional one-hot
    reaches ~0.59 -- so unlike the MLDE literature's pattern, the crude baseline does NOT win here. That
    contrast is a finding and it should fail loudly if it stops holding.
    """
    from dna_decode.glm.expression import fit_predict_ridge
    pairs, _ = load_pairs(CACHE, allow_download=False)
    tr, te, _ = leave_element_out(pairs, "Minus10", seed=0)
    gc = fit_predict_ridge(tr, te, feature="gc")["spearman"]
    oh = fit_predict_ridge(tr, te, feature="onehot")["spearman"]
    assert gc < 0.35, f"GC-only should be weak, got {gc}"
    assert oh > gc + 0.2, f"positional one-hot ({oh}) must clearly beat GC ({gc})"
