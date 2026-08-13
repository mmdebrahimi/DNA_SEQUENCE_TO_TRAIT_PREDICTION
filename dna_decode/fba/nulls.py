"""Nulls for the conditional-essentiality metric — including the margin-preserving one.

The rate-matched null (`scripts/fba_regulatory_conditional_test.py:rate_matched_null`) samples cells
INDEPENDENTLY: it draws k of the N gene x condition cells uniformly at random. That controls for the
BASE RATE and nothing else, and it was already recorded as the pFBA result's named weakness:

> "The rate-matched null samples cells INDEPENDENTLY, while real gene patterns are correlated across
>  conditions. A null preserving gene-level and condition-level margins would be a stronger test."

Why the weakness bites. Real essentiality is structured on BOTH axes:

  * **rows (genes)** -- a gene essential on one carbon source tends to be essential on chemically
    related ones; most genes are essential in 1-2 of 25 conditions, a few in many.
  * **columns (conditions)** -- a hard carbon source has more essential genes than an easy one
    (glycolate wild-type growth 0.153 vs maltose 1.780).

An independent-sampling null destroys both structures, so it is EASY to beat: some of the observed
agreement may come from merely reproducing the marginal shape rather than from getting the pattern right.
A null that preserves both margins asks the sharper question: **given that the predictor commits this
often per gene and this often per condition, is its PLACEMENT better than chance?**

`curveball_shuffle` implements swap randomization (Strona et al. 2014, *Nat Commun* 5:4114 -- the
"Curveball" algorithm), which samples uniformly from the set of binary matrices with the SAME row and
column sums. Each trade picks two rows, finds the columns where exactly one of them has a 1, and swaps a
random half of those -- an operation that provably leaves both margins invariant.

Pure: no solver, no model, no I/O. Deterministic given a seed.
"""
from __future__ import annotations

import random


def matrix_from_calls(gene_ids: list[str], conditions: tuple[str, ...],
                      calls: dict[str, dict[str, bool]]) -> list[list[int]]:
    """{condition: {gene: bool}} -> a genes x conditions 0/1 matrix, row order = `gene_ids`."""
    return [[1 if calls.get(c, {}).get(g, False) else 0 for c in conditions] for g in gene_ids]


def calls_from_matrix(gene_ids: list[str], conditions: tuple[str, ...],
                      mat: list[list[int]]) -> dict[str, dict[str, bool]]:
    """Inverse of `matrix_from_calls`, back into the shape `switch_accuracy` consumes."""
    return {c: {g: bool(mat[i][j]) for i, g in enumerate(gene_ids)}
            for j, c in enumerate(conditions)}


def margins(mat: list[list[int]]) -> tuple[list[int], list[int]]:
    """(row sums, column sums) — the two invariants a valid shuffle must preserve exactly."""
    rows = [sum(r) for r in mat]
    cols = [sum(mat[i][j] for i in range(len(mat))) for j in range(len(mat[0]))] if mat else []
    return rows, cols


def curveball_shuffle(mat: list[list[int]], n_trades: int | None = None,
                      seed: int = 0) -> list[list[int]]:
    """Swap-randomize a binary matrix, preserving BOTH row and column sums exactly.

    One trade: pick two rows; the columns where exactly one row has a 1 are the "swappable" set; deal a
    random half of them to each row. Cells where both rows agree are untouched, so every column sum is
    invariant by construction, and each row keeps its own count because it gives up exactly as many 1s as
    it receives.

    `n_trades` defaults to 5x the number of rows, the mixing heuristic in the original paper.
    """
    n_rows = len(mat)
    if n_rows < 2:
        return [row[:] for row in mat]
    out = [row[:] for row in mat]
    rng = random.Random(seed)
    trades = n_trades if n_trades is not None else 5 * n_rows

    for _ in range(trades):
        i, j = rng.randrange(n_rows), rng.randrange(n_rows)
        if i == j:
            continue
        a, b = out[i], out[j]
        # columns held by exactly one of the two rows -- the only ones a trade may move
        only_a = [k for k in range(len(a)) if a[k] and not b[k]]
        only_b = [k for k in range(len(a)) if b[k] and not a[k]]
        if not only_a or not only_b:
            continue
        pool = only_a + only_b
        rng.shuffle(pool)
        # a keeps exactly len(only_a) of the pool; b takes the rest -> both row sums preserved
        to_a = set(pool[:len(only_a)])
        for k in pool:
            a[k], b[k] = (1, 0) if k in to_a else (0, 1)
    return out


def margin_preserving_null(gene_ids: list[str], conditions: tuple[str, ...],
                           calls: dict[str, dict[str, bool]], score_fn,
                           n_draws: int = 200, seed0: int = 0) -> dict:
    """Score `n_draws` margin-preserving shuffles of `calls` through `score_fn(shuffled_calls)`.

    `score_fn` takes a calls-dict and returns a float (or None). Every draw preserves each gene's number
    of essential conditions AND each condition's number of essential genes, so a predictor only beats
    this null by placing its calls on the RIGHT cells -- not by matching the marginal shape.
    """
    mat = matrix_from_calls(gene_ids, conditions, calls)
    rows0, cols0 = margins(mat)
    scores: list[float] = []
    for s in range(n_draws):
        shuffled = curveball_shuffle(mat, seed=seed0 + s)
        r, c = margins(shuffled)
        if r != rows0 or c != cols0:                     # must never happen; fail loudly if it does
            raise AssertionError("curveball_shuffle broke a margin -- the null would be invalid")
        v = score_fn(calls_from_matrix(gene_ids, conditions, shuffled))
        if v is not None:
            scores.append(v)
    if not scores:
        return {"n_draws": 0, "mean": None, "max": None}
    scores.sort()
    return {
        "n_draws": len(scores),
        "mean": round(sum(scores) / len(scores), 4),
        "sd": round((sum((x - sum(scores) / len(scores)) ** 2 for x in scores)
                     / len(scores)) ** 0.5, 4),
        "max": round(scores[-1], 4),
        "p95": round(scores[int(0.95 * (len(scores) - 1))], 4),
        "scores": scores,
    }
