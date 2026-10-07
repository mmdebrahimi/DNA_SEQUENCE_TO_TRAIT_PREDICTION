"""Guards for the baselines and the falsifier.

The load-bearing tests are the ones that pin the falsifier's REFUSALS. A discriminator AUROC of 0.52 reads
as "indistinguishable, excellent" and is equally consistent with "the classifier learned nothing" — so the
instrument's own controls must be checked, and the verdict must refuse when they fail.
"""
from __future__ import annotations

import random

import pytest

from dna_decode.glm.baselines import MarkovGenerator, UniformGenerator
from dna_decode.glm.falsifier import (
    BEATS_NULL,
    INSTRUMENT_UNDERPOWERED,
    INSUFFICIENT,
    MIN_PER_CLASS,
    VERDICTS,
    WORSE,
    heldout_auroc,
    kmer_features,
    run_falsifier,
    self_control_ceiling,
)


def _at_rich(n, length, seed=0):
    """A synthetic 'natural-like' corpus: AT-rich, so it is separable from uniform."""
    r = random.Random(seed)
    return ["".join(r.choices("ACGT", weights=[35, 15, 15, 35], k=length)) for _ in range(n)]


# ---------------------------------------------------------------------------------------------------
# baselines
# ---------------------------------------------------------------------------------------------------
def test_uniform_generator_shape_and_determinism():
    a = UniformGenerator(seed=3).generate(5, 20)
    b = UniformGenerator(seed=3).generate(5, 20)
    assert a == b, "same seed must reproduce"
    assert len(a) == 5 and all(len(s) == 20 for s in a)
    assert all(not (set(s) - set("ACGT")) for s in a)


def test_markov_fit_then_generate_is_deterministic_and_well_formed():
    train = _at_rich(50, 60, seed=1)
    g = MarkovGenerator(order=3, seed=5).fit(train)
    a = g.generate(10, 50)
    b = MarkovGenerator(order=3, seed=5).fit(train).generate(10, 50)
    assert a == b
    assert all(len(s) == 50 for s in a)
    assert all(not (set(s) - set("ACGT")) for s in a)
    assert g.name == "markov-k3"


def test_markov_refuses_to_generate_before_fit():
    with pytest.raises(RuntimeError, match="before fit"):
        MarkovGenerator().generate(1, 50)


def test_markov_refuses_a_length_not_exceeding_its_order():
    g = MarkovGenerator(order=5, seed=0).fit(_at_rich(20, 60))
    with pytest.raises(ValueError, match="must exceed order"):
        g.generate(1, 5)


def test_markov_refuses_an_all_unusable_training_set_rather_than_emitting_uniform():
    """Fitted on nothing, add-one smoothing would make it silently emit uniform noise."""
    with pytest.raises(ValueError, match="no usable training sequences"):
        MarkovGenerator(order=3).fit(["AC", "NNNNNNN", ""])


def test_markov_order_must_be_positive():
    with pytest.raises(ValueError, match="order must be"):
        MarkovGenerator(order=0)


def test_markov_reproduces_the_composition_of_its_training_set():
    """This is WHY it is the right null: it matches local composition with no learning."""
    from dna_decode.glm.corpus import gc_fraction
    train = _at_rich(200, 100, seed=2)
    gen = MarkovGenerator(order=3, seed=0).fit(train).generate(200, 100)
    gc_t = sum(map(gc_fraction, train)) / len(train)
    gc_g = sum(map(gc_fraction, gen)) / len(gen)
    assert abs(gc_t - gc_g) < 0.05, f"null should match GC: {gc_t:.3f} vs {gc_g:.3f}"


# ---------------------------------------------------------------------------------------------------
# the discriminator
# ---------------------------------------------------------------------------------------------------
def test_kmer_features_are_normalised_frequencies_of_the_right_width():
    rows = kmer_features(["ACGTACGT", "AAAAAAAA"], k=2)
    assert len(rows) == 2 and len(rows[0]) == 16
    assert abs(sum(rows[0]) - 1.0) < 1e-9
    assert abs(sum(rows[1]) - 1.0) < 1e-9


def test_kmer_features_refuses_k_below_one():
    with pytest.raises(ValueError):
        kmer_features(["ACGT"], k=0)


def test_discriminator_HAS_power_on_an_obviously_separable_pair():
    at = _at_rich(120, 100, seed=1)
    gc = ["".join(random.Random(i).choices("GC", k=100)) for i in range(120)]
    auc = heldout_auroc(at, gc, k=3, seed=0)
    assert auc is not None and auc > 0.95, f"AT-rich vs GC-only must separate, got {auc}"


def test_discriminator_is_UNBIASED_on_two_halves_of_one_distribution():
    """The instrument's own null. Identical distributions must land near chance."""
    pool = _at_rich(300, 100, seed=4)
    auc = heldout_auroc(pool[:150], pool[150:], k=3, seed=0)
    assert auc is not None and 0.35 < auc < 0.65, f"same distribution should be ~0.5, got {auc}"


def test_discriminator_returns_None_rather_than_a_number_on_too_few_sequences():
    assert heldout_auroc(["ACGT"], ["TGCA"], k=2) is None


def test_self_control_ceiling_is_near_chance_and_returns_every_value():
    ceiling, vals = self_control_ceiling(_at_rich(300, 80, seed=9), k=3)
    assert ceiling is not None and len(vals) == 5
    assert ceiling == max(vals)
    assert 0.35 < ceiling < 0.70, f"a same-distribution ceiling should sit near chance, got {ceiling}"


# ---------------------------------------------------------------------------------------------------
# the verdict, and its refusals
# ---------------------------------------------------------------------------------------------------
def test_uniform_output_is_judged_WORSE_than_the_markov_null():
    nat = _at_rich(200, 100, seed=11)
    uni = UniformGenerator(seed=1).generate(200, 100)
    null = MarkovGenerator(order=3, seed=1).fit(nat).generate(200, 100)
    r = run_falsifier(nat, uni, candidate_name="uniform", null_sequences=null,
                      uniform_sequences=uni, k=3, length=100)
    assert r.verdict == WORSE, r.as_dict()
    assert r.distinguishability > r.null_distinguishability


def test_a_candidate_that_IS_the_natural_distribution_beats_the_null():
    """Non-vacuity for BEATS_NULL: a held-out slice of the real distribution must win."""
    pool = _at_rich(600, 100, seed=13)
    nat, twin = pool[:300], pool[300:]
    null = UniformGenerator(seed=2).generate(300, 100)   # a deliberately bad 'null'
    r = run_falsifier(nat, twin, candidate_name="twin", null_sequences=null,
                      uniform_sequences=null, k=3, length=100)
    assert r.verdict == BEATS_NULL, r.as_dict()


def test_the_verdict_REFUSES_when_the_power_control_fails():
    """THE guard. If the discriminator cannot separate uniform noise from 'natural', every other number
    it produces is unsupported — so it must refuse rather than report a flattering 0.5."""
    nat = _at_rich(200, 100, seed=17)
    # Hand the power control a set IDENTICAL in distribution to natural, so it cannot separate.
    fake_uniform = _at_rich(200, 100, seed=18)
    r = run_falsifier(nat, _at_rich(200, 100, seed=19), candidate_name="x",
                      null_sequences=_at_rich(200, 100, seed=20),
                      uniform_sequences=fake_uniform, k=3, length=100)
    assert r.verdict == INSTRUMENT_UNDERPOWERED, r.as_dict()
    assert r.distinguishability is None, "no measurement may be reported once the instrument is void"
    assert "cannot separate" in r.reason


def test_too_few_sequences_yields_INSUFFICIENT_not_a_number():
    r = run_falsifier(_at_rich(10, 50), _at_rich(10, 50, seed=2), candidate_name="tiny", k=3, length=50)
    assert r.verdict == INSUFFICIENT
    assert r.distinguishability is None
    assert str(MIN_PER_CLASS) in r.reason


def test_every_verdict_string_is_in_the_declared_set():
    nat = _at_rich(200, 100, seed=21)
    r = run_falsifier(nat, UniformGenerator(seed=1).generate(200, 100), candidate_name="u",
                      null_sequences=MarkovGenerator(seed=1).fit(nat).generate(200, 100),
                      uniform_sequences=UniformGenerator(seed=1).generate(200, 100), k=3, length=100)
    assert r.verdict in VERDICTS


def test_the_constraint_scope_is_declared_INAPPLICABLE_rather_than_reporting_a_vacuous_pass():
    """A filter that cannot fail is not a control; this repo has retracted that claim twice."""
    r = run_falsifier(_at_rich(5, 50), _at_rich(5, 50, seed=3), candidate_name="x", k=3, length=50)
    scope = r.checks["constraint_scope"]
    assert "INAPPLICABLE" in scope and "CDS-shaped" in scope
    assert "alphabet_valid_rate" in r.checks and "exact_length_rate" in r.checks


def test_alphabet_and_length_checks_discriminate_INDEPENDENTLY():
    """The two checks are separate fields because they catch different defects, and this pins that.

    `"ACG"` is alphabet-VALID but the wrong LENGTH; `"ACGN"` is the right length but alphabet-INVALID. A
    single merged "valid" rate would conflate a truncated generation with a model emitting ambiguity
    codes, which have different causes and different fixes.
    """
    from dna_decode.glm.falsifier import alphabet_and_length_ok
    good = alphabet_and_length_ok(["ACGT", "TGCA"], 4)
    assert good["alphabet_valid_rate"] == 1.0 and good["exact_length_rate"] == 1.0

    mixed = alphabet_and_length_ok(["ACGN", "ACG"], 4)
    assert mixed["alphabet_valid_rate"] == 0.5, "'ACG' is alphabet-valid; only 'ACGN' is not"
    assert mixed["exact_length_rate"] == 0.5, "'ACGN' is length-4; only 'ACG' is not"

    # each axis can fail alone -- which is the point of keeping them apart
    assert alphabet_and_length_ok(["ACGN"], 4)["alphabet_valid_rate"] == 0.0
    assert alphabet_and_length_ok(["ACGN"], 4)["exact_length_rate"] == 1.0
    assert alphabet_and_length_ok(["ACG"], 4)["alphabet_valid_rate"] == 1.0
    assert alphabet_and_length_ok(["ACG"], 4)["exact_length_rate"] == 0.0


# ---------------------------------------------------------------------------------------------------
# positional features -- the measured fix for the k-mer ceiling effect
# ---------------------------------------------------------------------------------------------------
def test_positional_features_have_the_right_width_and_are_per_bin_normalised():
    from dna_decode.glm.falsifier import positional_features
    rows = positional_features(["ACGT" * 15], bins=5, k=1)
    assert len(rows[0]) == 5 * 4
    for b in range(5):
        seg = rows[0][b * 4:(b + 1) * 4]
        assert abs(sum(seg) - 1.0) < 1e-9, "each bin normalises independently"


def test_positional_features_SEE_position_where_kmer_features_are_blind():
    """The whole justification. Two sequences with identical composition but reversed layout are
    IDENTICAL under a k-mer histogram and DIFFERENT under positional features."""
    from dna_decode.glm.falsifier import kmer_features, positional_features
    a = "A" * 40 + "C" * 40
    b = "C" * 40 + "A" * 40
    ka, kb = kmer_features([a, b], k=1)
    assert ka == kb, "a 1-mer histogram cannot tell these apart -- that is the blindness"
    pa, pb = positional_features([a, b], bins=4, k=1)
    assert pa != pb, "positional features must distinguish them"


def test_positional_features_validate_their_arguments():
    from dna_decode.glm.falsifier import positional_features
    with pytest.raises(ValueError, match="bins must be"):
        positional_features(["ACGT"], bins=0)
    with pytest.raises(ValueError, match="k must be"):
        positional_features(["ACGT"], bins=2, k=0)


def test_build_features_modes_compose_and_refuse_an_unknown_mode():
    from dna_decode.glm.falsifier import build_features
    seqs = ["ACGT" * 20, "TGCA" * 20]
    nk = len(build_features(seqs, mode="kmer", k=3)[0])
    npos = len(build_features(seqs, mode="positional", k=2, bins=5)[0])
    nboth = len(build_features(seqs, mode="both", k=3, bins=5)[0])
    assert nboth == nk + npos, "'both' must be the concatenation, not a replacement"
    with pytest.raises(ValueError, match="mode must be"):
        build_features(seqs, mode="nonsense")


def test_the_self_control_ceiling_is_RE_DERIVED_per_feature_mode():
    """A richer feature set has more capacity to fit split noise, so reusing a k-mer ceiling to judge a
    positional measurement would understate the noise and manufacture a finding. Measured on the real
    corpus: the ceiling rose 0.5012 -> 0.5292 moving from kmer to both."""
    pool = _at_rich(400, 120, seed=31)
    c_k, _ = self_control_ceiling(pool, k=3, mode="kmer")
    c_b, _ = self_control_ceiling(pool, k=3, mode="both", bins=6)
    assert c_k is not None and c_b is not None
    # both must stay near chance on one distribution; the point is that each is computed, not shared
    assert 0.3 < c_k < 0.7 and 0.3 < c_b < 0.7


def test_the_feature_mode_is_recorded_in_the_artifact():
    nat = _at_rich(200, 120, seed=33)
    r = run_falsifier(nat, UniformGenerator(seed=1).generate(200, 120), candidate_name="u",
                      null_sequences=MarkovGenerator(seed=1).fit(nat).generate(200, 120),
                      uniform_sequences=UniformGenerator(seed=1).generate(200, 120),
                      k=3, length=120, mode="both", bins=6)
    d = r.as_dict()
    assert d["feature_mode"] == "both" and d["positional_bins"] == 6


def test_bins_are_recorded_as_None_when_the_mode_does_not_use_them():
    """Reporting bins=10 for a kmer-only run would imply positional features were involved."""
    nat = _at_rich(200, 120, seed=35)
    r = run_falsifier(nat, _at_rich(200, 120, seed=36), candidate_name="x",
                      uniform_sequences=UniformGenerator(seed=1).generate(200, 120),
                      k=3, length=120, mode="kmer", bins=10)
    assert r.as_dict()["positional_bins"] is None


def test_as_dict_labels_the_metric_direction_because_lower_is_better_is_counterintuitive():
    r = run_falsifier(_at_rich(5, 50), _at_rich(5, 50, seed=4), candidate_name="x", k=3, length=50)
    assert "LOWER_IS_BETTER" in "".join(r.as_dict().keys())
