"""Baseline sequence generators — the bar a learned generator has to clear to have earned anything.

Two of them, and they play OPPOSITE roles, which is why they live together:

* **`MarkovGenerator` — the NULL.** An order-k Markov chain fitted on the natural corpus reproduces local
  composition (GC, k-mer frequencies, codon-ish structure) with no learning and no model. Published work in
  this area repeatedly finds that a model which merely matches composition is matched by a baseline that
  *only* matches composition, so **a learned generator that does not beat this has demonstrated nothing**.
  This repo already has the general form of that lesson measured several times over: a single integer
  (mutation count) recovers 75% of a protein LM's headline, and a 1992 substitution matrix ties ESM2 on
  resistance ranking.
* **`UniformGenerator` — the POWER control.** I.i.d. uniform bases are obviously not DNA, so a
  discriminator that cannot separate *this* from the natural corpus has no power and its verdict on a real
  generator is void. It exists to make the falsifier's own instrument checkable.

Both are pure-Python + `random.Random(seed)`, so they are deterministic, need no GPU, no weights and no
network — which is what lets the whole falsifier be validated before any 4 GiB download.
"""
from __future__ import annotations

import random
from collections import defaultdict

BASES = "ACGT"


class UniformGenerator:
    """I.i.d. uniform bases. The discriminator's POWER control, not a serious generator."""

    name = "uniform"

    def __init__(self, seed: int = 0):
        self.rng = random.Random(seed)

    def generate(self, n: int, length: int) -> list[str]:
        return ["".join(self.rng.choice(BASES) for _ in range(length)) for _ in range(n)]


class MarkovGenerator:
    """Order-k Markov chain fitted on real sequence. THE NULL a learned generator must beat.

    Fitting is counts of (context -> next base) over the training sequences, with add-one smoothing so an
    unseen context cannot dead-end generation. The chain is seeded from a real k-mer drawn from the
    training set's observed start contexts, so generated sequences begin the way real ones do rather than
    from a uniform k-mer that the model was never fitted on.
    """

    def __init__(self, order: int = 3, seed: int = 0):
        if order < 1:
            raise ValueError("order must be >= 1")
        self.order = order
        self.rng = random.Random(seed)
        self._counts: dict[str, list[int]] = {}
        self._starts: list[str] = []
        self._fitted = False

    @property
    def name(self) -> str:
        return f"markov-k{self.order}"

    def fit(self, seqs: list[str]) -> MarkovGenerator:
        k = self.order
        counts: dict[str, list[int]] = defaultdict(lambda: [1, 1, 1, 1])  # add-one smoothing
        starts: list[str] = []
        usable = 0
        for s in seqs:
            s = s.upper()
            if len(s) <= k or set(s) - set(BASES):
                continue
            usable += 1
            starts.append(s[:k])
            for i in range(len(s) - k):
                ctx = s[i:i + k]
                nxt = s[i + k]
                counts[ctx][BASES.index(nxt)] += 1
        if usable == 0:
            raise ValueError(
                f"no usable training sequences (need length > order={k} and ACGT-only); "
                "a chain fitted on nothing would silently emit smoothing-uniform output"
            )
        self._counts = dict(counts)
        self._starts = starts
        self._fitted = True
        return self

    def _next(self, ctx: str) -> str:
        w = self._counts.get(ctx)
        if w is None:
            w = [1, 1, 1, 1]
        return self.rng.choices(BASES, weights=w, k=1)[0]

    def generate(self, n: int, length: int) -> list[str]:
        if not self._fitted:
            raise RuntimeError("MarkovGenerator.generate called before fit()")
        if length <= self.order:
            raise ValueError(f"length ({length}) must exceed order ({self.order})")
        out = []
        for _ in range(n):
            seq = list(self.rng.choice(self._starts))
            while len(seq) < length:
                seq.append(self._next("".join(seq[-self.order:])))
            out.append("".join(seq))
        return out
