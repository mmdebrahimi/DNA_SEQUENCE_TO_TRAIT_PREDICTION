"""THREE k-mer implementations exist in `dna_decode/glm/`, and their agreement is COINCIDENTAL.

Found while surveying reusable code for the G-A build:

- `falsifier.kmer_features`   -- used by the naturalness discriminator AND the reference-geometry triage
- `genomewide.kmer_features`  -- used by the genome-wide oracle arms
- `expression.kmer_feature`   -- delegates to the falsifier one

Nothing enforced that the first two return the same thing. They do today (verified: identical output, width
256 at k=4), but that is an accident of two independent implementations happening to normalise the same way.

**Why it matters rather than being tidiness.** G-A's gate compares a learned encoder against k-mer arms drawn
from `genomewide`, while the reference-geometry result it is interpreted against (grid vs tile AUROC 1.0000)
was computed with `falsifier`. If either drifts -- a different normaliser, a different handling of ambiguous
bases, a different k-mer ordering -- the two numbers silently stop being comparable, and nothing anywhere
would fail. That is the same class as the two GC baselines this project already had to reconcile, caught
before it happened instead of after.
"""
from __future__ import annotations

import random

import pytest

from dna_decode.glm.expression import kmer_feature as expr_kmer
from dna_decode.glm.falsifier import kmer_features as fals_kmer
from dna_decode.glm.genomewide import kmer_features as geno_kmer


def _seqs(n=12, lengths=(48, 150, 244, 499), seed=0):
    r = random.Random(seed)
    return ["".join(r.choice("ACGT") for _ in range(r.choice(lengths))) for _ in range(n)]


# ---------------------------------------------------------------------------------------------------
# NON-VACUITY FIRST -- the guard is worthless if these are the same object
# ---------------------------------------------------------------------------------------------------
def test_the_implementations_are_GENUINELY_SEPARATE_objects():
    """If `genomewide.kmer_features` were just a re-export of the falsifier one, every parity assertion
    below would pass trivially and guard nothing. Assert they are distinct code before trusting agreement.
    """
    assert fals_kmer is not geno_kmer, "parity guard is vacuous if these are the same function"
    assert fals_kmer.__code__ is not geno_kmer.__code__, "distinct names but shared code object"
    assert fals_kmer.__module__ != geno_kmer.__module__


def test_expression_DOES_delegate_and_that_is_recorded_not_assumed():
    """`expression.kmer_feature` is a thin wrapper over the falsifier implementation. Pinning the
    delegation means a future refactor that gives it its own body shows up here rather than silently
    creating a FOURTH implementation."""
    s = _seqs(5, seed=7)
    assert expr_kmer(s, 4) == fals_kmer(s, k=4)


# ---------------------------------------------------------------------------------------------------
# parity across k and across length
# ---------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("k", [1, 2, 4, 6])
def test_falsifier_and_genomewide_agree_exactly_for_each_k(k):
    """Exact equality, not approximate: both normalise to frequencies over the same alphabet ordering, so
    any difference at all is a drift signal rather than float noise."""
    s = _seqs(16, seed=k)
    a, b = fals_kmer(s, k=k), geno_kmer(s, k)
    assert len(a) == len(b) == len(s)
    for i, (x, y) in enumerate(zip(a, b)):
        assert x == y, f"k={k} diverged on sequence {i} (len {len(s[i])})"


@pytest.mark.parametrize("k", [1, 4, 6])
def test_widths_are_4_to_the_k_in_both(k):
    s = _seqs(3, seed=99)
    assert len(fals_kmer(s, k=k)[0]) == 4 ** k
    assert len(geno_kmer(s, k)[0]) == 4 ** k


def test_both_normalise_to_a_unit_sum():
    """A frequency vector, not a count vector. If one switched to raw counts the other's ridge arms would
    be on a different scale and the comparison would be meaningless rather than merely different."""
    s = _seqs(6, seed=3)
    for impl in (lambda q: fals_kmer(q, k=4), lambda q: geno_kmer(q, 4)):
        for v in impl(s):
            assert sum(v) == pytest.approx(1.0)


def test_parity_holds_on_the_SHORTEST_real_fragment_length():
    """The real fragment library goes down to 48bp. At k=6 a 48bp sequence yields only 43 k-mers, so any
    off-by-one in the window loop shows up at the short end first."""
    s = ["".join(random.Random(11).choice("ACGT") for _ in range(48))]
    assert fals_kmer(s, k=6) == geno_kmer(s, 6)


def test_parity_holds_when_k_EXCEEDS_the_sequence_length():
    """A degenerate input both must handle the same way. Whatever they return, they must return the SAME
    thing -- otherwise one silently contributes a zero vector where the other contributes something else.
    """
    s = ["ACG"]
    assert fals_kmer(s, k=6) == geno_kmer(s, 6)


def test_both_are_deterministic_and_order_preserving():
    s = _seqs(8, seed=5)
    assert geno_kmer(s, 4) == geno_kmer(s, 4)
    assert fals_kmer(s, k=4) == fals_kmer(s, k=4)
    # reversing the input must permute the output identically, never change the values
    assert geno_kmer(s, 4)[::-1] == geno_kmer(s[::-1], 4)
