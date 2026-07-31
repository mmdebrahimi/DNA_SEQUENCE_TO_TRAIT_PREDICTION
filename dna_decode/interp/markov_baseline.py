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

    def base_distribution(self, left_context: str, exclude_base: str | None = None,
                          alpha: float = 1.0) -> dict:
        """Laplace-smoothed predictive distribution over {A,C,G,T} given the left context.

        Backs off k -> 0 to the highest-order context with any counts. `exclude_base` implements
        LEAVE-ONE-OUT: decrement that base's count by 1 before normalizing, so a target position is
        never predicted using its own observation (removes the same-region transductive leakage).
        `alpha` is add-alpha smoothing (also the only mass when the context is unseen -> uniform).
        """
        left = left_context.upper()
        counts = {b: 0 for b in _BASES}
        top_order = min(self.k, len(left))
        for order in range(top_order, -1, -1):
            context = left[-order:] if order > 0 else ""
            dist = self._tables[order].get(context)
            if dist:
                counts = {b: dist.get(b, 0) for b in _BASES}
                break
        if exclude_base in _BASES and counts[exclude_base] > 0:
            counts[exclude_base] -= 1
        smoothed = {b: counts[b] + alpha for b in _BASES}
        total = sum(smoothed.values())
        if total <= 0:  # alpha=0 and all counts excluded/zero -> uniform
            return {b: 0.25 for b in _BASES}
        return {b: smoothed[b] / total for b in _BASES}

    def predict_base(self, left_context: str, exclude_base: str | None = None) -> str:
        """Argmax next base (deterministic tie-break count desc, then A<C<G<T). `exclude_base` = LOO."""
        dist = self.base_distribution(left_context, exclude_base=exclude_base, alpha=0.0
                                      if any(self._tables[0].values()) else 1.0)
        return max(_BASES, key=lambda base: (dist[base], -_BASES.index(base)))

    def accuracy_on_masked(self, full_sequence: str, masked_indices,
                           leave_one_out: bool = False) -> float:
        """Per-base argmax accuracy predicting each masked index from its TRUE left context.

        `leave_one_out=True` excludes each target's own count (use when the Markov was fit on the
        SAME sequence being scored; a no-op harm when fit on a disjoint region). `masked_indices`
        are 0-based; non-ACGT targets are skipped.
        """
        s = full_sequence.upper()
        correct = scored = 0
        for j in masked_indices:
            true_base = s[j]
            if true_base not in _BASES:
                continue
            left = s[max(0, j - self.k):j]
            scored += 1
            if self.predict_base(left, exclude_base=true_base if leave_one_out else None) == true_base:
                correct += 1
        return correct / scored if scored else 0.0

    def nll_on_masked(self, full_sequence: str, masked_indices,
                      leave_one_out: bool = False, alpha: float = 1.0) -> float:
        """Mean per-base negative log-likelihood (nats) of the true base under the Markov predictive
        distribution. The native surface for a probability model; the primary F1' endpoint on the
        Markov side. `leave_one_out` excludes each target's own count (same-region-fit leakage fix)."""
        import math

        s = full_sequence.upper()
        total_nll = scored = 0
        for j in masked_indices:
            true_base = s[j]
            if true_base not in _BASES:
                continue
            left = s[max(0, j - self.k):j]
            dist = self.base_distribution(left, exclude_base=true_base if leave_one_out else None,
                                          alpha=alpha)
            total_nll += -math.log(max(dist[true_base], 1e-12))
            scored += 1
        return total_nll / scored if scored else float("inf")
