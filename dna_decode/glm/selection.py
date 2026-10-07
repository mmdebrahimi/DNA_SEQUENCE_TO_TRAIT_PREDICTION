"""Does oracle-guided selection actually DELIVER? The downstream-task measurement.

**Why this exists.** The expression oracle reports Spearman rho ~0.59 on held-out promoter elements. That is
an UPSTREAM metric. The generative loop does not consume a correlation — it consumes a *choice*: propose k
candidates, score them, keep the best. Whether rho 0.59 converts into a useful choice is **not** derivable
from rho: it depends on the shape of the tail, on how ties fall, and on how much selection value is even
available in a pool of k. A correlation of 0.59 is compatible with capturing nearly all of the available
headroom or almost none of it.

This is the same move AlphaGenome's own authors make with self-distillation: the distillation metric proves
nothing, so the claim is validated DOWNSTREAM on real measured variant-effect benchmarks. The analogue here
is to stop quoting rho and measure the selection the loop will actually perform, against real measured
expression, on sequences the oracle has never seen.

**The design, and what each control is for.**

- `picked`   -- the scorer's argmax over k candidates. What the loop does.
- `random`   -- a uniformly random member of the SAME k candidates. The null. **Paired, not independent**:
                sharing the candidate set is what makes the comparison a measurement of the SCORER rather
                than of two different draws (this repo has a standing lesson on paired comparison).
- `ceiling`  -- the true best of those same k by MEASURED expression. The most any scorer could achieve.
                Without it a lift is uninterpretable, because the available headroom varies with k.

`headroom_captured = (picked - random) / (ceiling - random)` is therefore the informative number: the
fraction of the achievable selection value the scorer actually realises.

**Two correctness properties hold BY CONSTRUCTION and are asserted rather than hoped for.**

1. **At k = 1 the three arms are identical.** With one candidate, argmax *is* the random pick *is* the best.
   Any lift at k = 1 means the simulation is broken, so k = 1 is a free built-in self-test and
   `headroom_captured` is `None` there (a 0/0, never a fabricated number).
2. **Ties break RANDOMLY, not by index order.** `random.sample` returns its subset in random order and
   `max` keeps the first maximum, so tied predictions resolve uniformly. This matters: this repo has
   already measured a case where `sorted()`-order tie-breaking silently shifted a median, and a scorer that
   emits many ties (GC content over 150 bp emits few distinct values) would otherwise be credited with
   whatever the index order happened to give it.

Percentiles use **mid-ranks** for the same reason — the assay carries repeated expression values.

Pure functions only: no I/O, no model, no data loading, so every property above is testable offline.
"""
from __future__ import annotations

import random
from dataclasses import dataclass


def midrank_percentiles(values: list[float]) -> list[float]:
    """Mid-rank percentile in [0, 1] of each value within `values`; tied values share one percentile.

    Mid-ranks, not ordinal ranks: with repeated measurements an ordinal rank hands the whole tie-block's
    spread to whichever element sorted first, which is exactly the trap that shifted a published median
    elsewhere in this repo.
    """
    n = len(values)
    if n == 0:
        return []
    if n == 1:
        return [0.5]
    order = sorted(range(n), key=lambda i: values[i])
    pct = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        mid = (i + j) / 2.0
        p = mid / (n - 1)
        for t in range(i, j + 1):
            pct[order[t]] = p
        i = j + 1
    return pct


@dataclass(frozen=True)
class SelectionResult:
    """One (scorer, k) cell of the selection measurement."""

    scorer: str
    k: int
    n_trials: int
    n_pool: int
    mean_pct_picked: float
    mean_pct_random: float
    mean_pct_ceiling: float
    median_measured_picked: float
    median_measured_random: float
    median_measured_ceiling: float
    fold_over_random: float
    headroom_captured: float | None
    win_rate_vs_random: float
    tie_rate_vs_random: float

    def as_dict(self) -> dict:
        return {
            "scorer": self.scorer,
            "k": self.k,
            "n_trials": self.n_trials,
            "n_pool": self.n_pool,
            "mean_pct_picked": round(self.mean_pct_picked, 4),
            "mean_pct_random": round(self.mean_pct_random, 4),
            "mean_pct_ceiling": round(self.mean_pct_ceiling, 4),
            "median_measured_picked": round(self.median_measured_picked, 4),
            "median_measured_random": round(self.median_measured_random, 4),
            "median_measured_ceiling": round(self.median_measured_ceiling, 4),
            "fold_over_random": round(self.fold_over_random, 4),
            "headroom_captured": (None if self.headroom_captured is None
                                  else round(self.headroom_captured, 4)),
            "win_rate_vs_random": round(self.win_rate_vs_random, 4),
            "tie_rate_vs_random": round(self.tie_rate_vs_random, 4),
        }


def _median(xs: list[float]) -> float:
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return 0.0
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def simulate_selection(predicted: list[float], measured: list[float], *, k: int,
                       scorer: str = "scorer", trials: int = 2000, seed: int = 0,
                       percentiles: list[float] | None = None) -> SelectionResult:
    """Simulate propose-k / score-k / keep-best and measure what the pick is actually worth.

    `predicted` is the scorer's output, `measured` the real assay value, index-aligned. Candidate sets are
    drawn WITHOUT replacement from the pool, and the random and ceiling arms use the SAME set as the pick.
    """
    if len(predicted) != len(measured):
        raise ValueError(f"predicted ({len(predicted)}) and measured ({len(measured)}) must align")
    n = len(measured)
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if k > n:
        raise ValueError(f"k={k} exceeds pool size {n}; cannot draw k candidates without replacement")
    if trials < 1:
        raise ValueError(f"trials must be >= 1, got {trials}")

    pct = midrank_percentiles(measured) if percentiles is None else percentiles
    if len(pct) != n:
        raise ValueError("percentiles must align with measured")

    rng = random.Random(seed)
    p_pick, p_rand, p_ceil = [], [], []
    m_pick, m_rand, m_ceil = [], [], []
    wins = ties = 0
    for _ in range(trials):
        idx = rng.sample(range(n), k)
        # idx is in random order, so `max` breaking ties on the first maximum breaks them UNIFORMLY.
        pick = max(idx, key=lambda i: predicted[i])
        rand = idx[0]           # uniform over the pool, and paired with this candidate set
        ceil = max(idx, key=lambda i: measured[i])
        p_pick.append(pct[pick]); p_rand.append(pct[rand]); p_ceil.append(pct[ceil])
        m_pick.append(measured[pick]); m_rand.append(measured[rand]); m_ceil.append(measured[ceil])
        if measured[pick] > measured[rand]:
            wins += 1
        elif measured[pick] == measured[rand]:
            ties += 1

    mp, mr, mc = (sum(p_pick) / trials, sum(p_rand) / trials, sum(p_ceil) / trials)
    head = mc - mr
    med_r = _median(m_rand)
    return SelectionResult(
        scorer=scorer, k=k, n_trials=trials, n_pool=n,
        mean_pct_picked=mp, mean_pct_random=mr, mean_pct_ceiling=mc,
        median_measured_picked=_median(m_pick),
        median_measured_random=med_r,
        median_measured_ceiling=_median(m_ceil),
        fold_over_random=(_median(m_pick) / med_r if med_r else float("nan")),
        # 0/0 at k=1, where all three arms are identical by construction. None, never a fabricated number.
        headroom_captured=(None if abs(head) < 1e-12 else (mp - mr) / head),
        win_rate_vs_random=wins / trials,
        tie_rate_vs_random=ties / trials,
    )


def permutation_null(measured: list[float], *, k: int, trials: int = 400, n_perms: int = 40,
                     seed: int = 0) -> dict:
    """The null band for `headroom_captured` at this pool and this k. **Load-bearing, and measured.**

    A single cell's headroom carries SCORER-REALISATION uncertainty that `trials` cannot reduce. Measured on
    synthetic data: 60 independent useless scorers over a 300-element pool at k=10 gave mean headroom
    **0.0006** (correct -- a useless scorer captures none of the headroom) with **stdev 0.096**, and 31 of 60
    were negative. Raising `trials` shrinks the Monte-Carlo error in estimating ONE scorer's headroom and
    says nothing about how far that scorer sits from zero.

    So a reported headroom is only interpretable beside this band. A shuffled scorer is a uniform random
    permutation whatever it was shuffled FROM, so the null is a property of (pool, k) alone and is computed
    ONCE per cell rather than once per scorer.
    """
    n = len(measured)
    pct = midrank_percentiles(measured)
    rng = random.Random(seed)
    base = list(range(n))
    hs: list[float] = []
    for p in range(n_perms):
        shuffled = base[:]
        rng.shuffle(shuffled)
        fake = [float(v) for v in shuffled]       # a permutation IS the shuffled-scorer null
        r = simulate_selection(fake, measured, k=k, trials=trials, seed=seed + 1000 + p,
                              percentiles=pct)
        if r.headroom_captured is not None:
            hs.append(r.headroom_captured)
    if not hs:
        return {"n_perms": 0, "note": "degenerate at this k (headroom is 0/0); no null band"}
    hs.sort()
    mean = sum(hs) / len(hs)
    var = sum((h - mean) ** 2 for h in hs) / (len(hs) - 1) if len(hs) > 1 else 0.0
    return {
        "n_perms": len(hs),
        "mean": round(mean, 4),
        "stdev": round(var ** 0.5, 4),
        "min": round(hs[0], 4),
        "max": round(hs[-1], 4),
        "p95": round(hs[min(len(hs) - 1, int(0.95 * len(hs)))], 4),
    }


def assert_k1_is_degenerate(result: SelectionResult) -> None:
    """k=1 self-test: with one candidate the pick, the null and the ceiling must coincide exactly.

    This is the simulation's own correctness check — a lift at k=1 is impossible in a correct
    implementation, so it localises a bug to the harness rather than to the scorer.
    """
    if result.k != 1:
        raise ValueError(f"only meaningful at k=1, got k={result.k}")
    if not (abs(result.mean_pct_picked - result.mean_pct_random) < 1e-12
            and abs(result.mean_pct_picked - result.mean_pct_ceiling) < 1e-12):
        raise AssertionError(
            f"k=1 must be degenerate but arms differ: picked={result.mean_pct_picked} "
            f"random={result.mean_pct_random} ceiling={result.mean_pct_ceiling} -- harness bug"
        )
    if result.headroom_captured is not None:
        raise AssertionError("k=1 headroom is 0/0 and must be None, not a number")
