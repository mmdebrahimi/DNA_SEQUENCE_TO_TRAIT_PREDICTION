"""Order-k Markov nucleotide baseline for masked-reconstruction DELTA scoring.

The world-model reconstruction bar is NEVER raw accuracy -- genomes are repetitive and a cheap
order-k Markov model already reconstructs a large fraction of masked bases with no "understanding".
The honest headline for the dog world-model probe is therefore

    delta = (LM per-base reconstruction accuracy) - (order-k Markov per-base accuracy)

This module is the cheap baseline. It is teacher-forced on the TRUE preceding bases (the strongest
honest form of a causal baseline), so the delta is conservative -- the LM has to beat a baseline that
is handed the real left context. Pure python, offline, no torch.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

_BASES = ("A", "C", "G", "T")


@dataclass
class MarkovModel:
    """Order-k Markov model over {A,C,G,T} with full backoff (k -> 0)."""

    k: int
    _tables: list  # _tables[order]: {context(str) -> {base: count}}; order 0 context is ""

    @classmethod
    def fit(cls, sequences, k: int) -> "MarkovModel":
        if k < 0:
            raise ValueError("k must be >= 0")
        tables = [defaultdict(lambda: defaultdict(int)) for _ in range(k + 1)]
        for seq in sequences:
            s = seq.upper()
            for i, b in enumerate(s):
                if b not in _BASES:
                    continue
                for order in range(0, k + 1):
                    if i < order:
                        continue
                    context = s[i - order:i]
                    if all(c in _BASES for c in context):
                        tables[order][context][b] += 1
        return cls(k=k, _tables=[{c: dict(d) for c, d in t.items()} for t in tables])

    def predict_base(self, left_context: str) -> str:
        """Argmax next base given the bases immediately to the left, backing off k -> 0.

        Deterministic tie-break: highest count, then A<C<G<T. Falls to "A" only if the model
        saw no bases at all (empty fit).
        """
        left = left_context.upper()
        top_order = min(self.k, len(left))
        for order in range(top_order, -1, -1):
            context = left[-order:] if order > 0 else ""
            dist = self._tables[order].get(context)
            if dist:
                return max(_BASES, key=lambda base: (dist.get(base, 0), -_BASES.index(base)))
        return "A"

    def accuracy_on_masked(self, full_sequence: str, masked_indices) -> float:
        """Per-base accuracy predicting each masked absolute index from its TRUE left context.

        Teacher-forced: the left context is the real sequence (not the model's own fills), so this
        is the STRONGEST form of the causal baseline. `masked_indices` are 0-based positions into
        `full_sequence`; positions whose true base is not ACGT are skipped.
        """
        s = full_sequence.upper()
        correct = 0
        scored = 0
        for j in masked_indices:
            true_base = s[j]
            if true_base not in _BASES:
                continue
            left = s[max(0, j - self.k):j]
            scored += 1
            if self.predict_base(left) == true_base:
                correct += 1
        return correct / scored if scored else 0.0
