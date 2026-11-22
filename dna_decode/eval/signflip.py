"""Family-wise sign-flip permutation over whole replicate rows.

RULE_2 of wiki/human_replication_and_multiplicity_protocol_2026-11-22.json, relayed from the P2
tokenization-axis session: "Reporting whichever candidate cleared a per-candidate 95% CI among N
examined prices nothing. A sign-flip permutation over WHOLE replicate rows keeps every candidate's
value for that row together, so the cross-candidate correlation is carried into the null instead of
assumed away, and yields one family-wise p with no independence assumption."

WHAT A ROW IS, AND WHY IT MUST BE A WHOLE ROW. A row is ONE replicate -- for the human arm, one
TRAINED RUN -- holding that replicate's value for every candidate in the family (e.g. the four model
sizes s/m/l/xl). Flipping the sign of an entire row keeps those values together, so whatever makes
run 7 good for size s and also good for size l stays in the null. Flipping cells independently would
destroy exactly that structure and manufacture significance, which is the error this exists to
prevent.

WHAT THE VALUES MUST BE. Each cell is a PAIRED DIFFERENCE (candidate minus its baseline for that
replicate), because sign-flipping tests the null that the values are symmetric about ZERO. Handing
it raw scores tests whether those scores are symmetric about zero, which is a different and usually
meaningless question. `flag_uncentred` NOTES, and deliberately does not refuse, a matrix whose every
column sits far from zero relative to its own spread: that is the signature of raw scores, and
equally the signature of a strong consistent low-variance real effect, which is the most valuable
input this function can receive. They are not separable from the numbers, so the caller is warned
rather than overruled.

WHERE THIS APPLIES. Per the protocol: EXEMPT for the pre-specified primary verdict (M = 1) and for
the two per-cell reports; EXPOSED for the four-size curve and for any adaptive follow-on. Deciding
which case you are in is NOT this module's job -- it computes a family-wise p for a family you
declare. It does, however, refuse a one-candidate "family", because a family-wise correction over a
single candidate is a per-candidate test wearing the wrong name.

Also note what this does NOT do: it corrects for SELECTION across the family. It cannot repair a
wrong unit of replication. Rows must already be independent replicates -- for an architecture or
size claim that means separate TRAINED RUNS, not re-scorings of one generated set and not extra
sampling seeds from one trained model. That is RULE_1's job and it must be settled first.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

__all__ = [
    "SignFlipRefusal",
    "SignFlipResult",
    "family_wise_signflip",
    "bonferroni_p",
]

#: Above this many rows, exact enumeration of 2**n sign vectors is replaced by Monte Carlo.
#: 2**16 = 65536 complete enumerations is cheap; 2**24 is not.
EXACT_MAX_ROWS = 16


class SignFlipRefusal(ValueError):
    """The input cannot support a family-wise sign-flip test.

    Raised rather than returning a p-value, because every refusal here describes a case where a
    number WOULD be produced and would not mean what a reader takes it to mean.
    """


@dataclass(frozen=True)
class SignFlipResult:
    n_rows: int
    n_candidates: int
    statistic: str
    observed: list[float]
    observed_family: float
    family_wise_p: float
    per_candidate_p: list[float]
    n_draws: int
    exact: bool
    seed: int | None
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "test": "family-wise-signflip-over-whole-replicate-rows-v1",
            "n_rows_replicates": self.n_rows,
            "n_candidates_in_family": self.n_candidates,
            "statistic": self.statistic,
            "observed_per_candidate": [round(v, 6) for v in self.observed],
            "observed_family_statistic": round(self.observed_family, 6),
            "family_wise_p": round(self.family_wise_p, 6),
            "per_candidate_p_UNCORRECTED": [round(v, 6) for v in self.per_candidate_p],
            "n_draws": self.n_draws,
            "enumeration": "exact" if self.exact else "monte-carlo",
            "seed": self.seed,
            "reading": (
                "family_wise_p is ONE p for the whole family and assumes NO independence across "
                "candidates -- the cross-candidate correlation is carried in the null by flipping "
                "whole rows. per_candidate_p is UNCORRECTED and must not be reported as the "
                "survival of a candidate selected from this family."),
            "what_it_does_not_cover": (
                "selection across the family only. It cannot repair a wrong unit of replication: "
                "rows must already be independent replicates (separate trained runs for a size or "
                "architecture claim), which is RULE_1's question and must be settled first."),
            "notes": self.notes,
        }


def _validate(rows) -> list[list[float]]:
    if rows is None:
        raise SignFlipRefusal("rows is None")
    mat = [list(r) for r in rows]
    if len(mat) < 2:
        raise SignFlipRefusal(
            f"{len(mat)} replicate row(s): a sign-flip null over fewer than 2 rows has at most 2 "
            f"distinct draws and cannot produce a usable p. This is the UNIT-OF-REPLICATION "
            f"question (RULE_1), not a tuning parameter -- collapse scoring splits per generated "
            f"set, and use separate TRAINED RUNS for an architecture or size claim.")
    widths = {len(r) for r in mat}
    if len(widths) != 1:
        raise SignFlipRefusal(
            f"ragged rows: candidate counts {sorted(widths)}. Every replicate must carry a value "
            f"for every candidate, or the rows being flipped are not comparable vectors.")
    k = widths.pop()
    if k == 0:
        raise SignFlipRefusal("rows have zero candidates")
    if k == 1:
        raise SignFlipRefusal(
            "a one-candidate 'family' is a per-candidate test. A family-wise correction over a "
            "single candidate computes the same p and invites it to be reported as though a "
            "selection had been priced. Call this only on the family you actually chose from.")
    for i, r in enumerate(mat):
        for j, v in enumerate(r):
            fv = float(v)
            if math.isnan(fv) or math.isinf(fv):
                raise SignFlipRefusal(
                    f"row {i}, candidate {j} is {v!r}. A non-finite cell silently propagates "
                    f"through every draw; drop the replicate deliberately or explain the cell.")
            mat[i][j] = fv
    return mat


def _t_stats(mat: list[list[float]], statistic: str) -> list[float]:
    n = len(mat)
    out = []
    for j in range(len(mat[0])):
        col = [row[j] for row in mat]
        mean = sum(col) / n
        if statistic == "mean":
            out.append(mean)
            continue
        var = sum((v - mean) ** 2 for v in col) / (n - 1)
        # An EXACT `var <= 0` test does not fire on a constant column: [0.05, 0.05, 0.05] sums to
        # 0.15000000000000002, so the mean is off by ~1e-17, var lands near 3e-34, and t comes out
        # around 5e15 -- an enormous finite number where the honest answer is "undefined". Checked
        # against the column's own span AND a relative tolerance, because a guard that misses the
        # case it exists for is not a guard.
        span = max(col) - min(col)
        scale = max(abs(mean), max(abs(v) for v in col), 1e-300)
        if span == 0.0 or var <= 0.0 or math.sqrt(var) <= 1e-12 * scale:
            raise SignFlipRefusal(
                f"candidate {j} has ZERO (or numerically indistinguishable from zero) variance "
                f"across replicates (every row = {col[0]!r}), so its t statistic is undefined "
                f"rather than merely large. That is a real state, not a bug: pass "
                f"statistic='mean' to test it, and report that the column was constant.")
        out.append(mean / math.sqrt(var / n))
    return out


def _sign_vectors(n: int, n_draws: int, seed: int | None):
    """Yield (signs, exact, total). Exact enumeration includes the observed all-plus vector."""
    if n <= EXACT_MAX_ROWS:
        total = 2 ** n
        return (list(s) for s in itertools.product((1.0, -1.0), repeat=n)), True, total
    import random
    rng = random.Random(seed)

    def gen():
        # The observed configuration is included deliberately, so p can never be 0.
        yield [1.0] * n
        for _ in range(n_draws - 1):
            yield [1.0 if rng.random() < 0.5 else -1.0 for _ in range(n)]

    return gen(), False, n_draws


def family_wise_signflip(
    rows,
    *,
    statistic: str = "max_abs_t",
    n_draws: int = 20000,
    seed: int | None = 0,
    flag_uncentred: bool = True,
) -> SignFlipResult:
    """One family-wise p for a family of candidates measured on shared replicate rows.

    `rows[i][j]` is replicate i's PAIRED DIFFERENCE for candidate j. Signs are flipped per ROW, so
    the within-row correlation across candidates is preserved in the null.
    """
    if statistic not in ("max_abs_t", "mean"):
        raise SignFlipRefusal(f"unknown statistic {statistic!r}; use 'max_abs_t' or 'mean'")
    mat = _validate(rows)
    n, k = len(mat), len(mat[0])
    notes: list[str] = []

    if flag_uncentred:
        col_means = [sum(r[j] for r in mat) / n for j in range(k)]
        spread = max(abs(m) for m in col_means)
        pooled = [v for r in mat for v in r]
        scale = (max(pooled) - min(pooled)) or 1.0
        if spread / scale > 5.0:
            # A NOTE, not a refusal. This pattern -- every column tightly clustered far from zero
            # -- is the signature of raw scores, and ALSO the signature of a strong, consistent,
            # low-variance real effect, which is the single most valuable input this function can
            # receive. Refusing it (as a first version did) blocks the best legitimate case to
            # catch a caller error that the docstring already states. Undecidable from the numbers,
            # so the caller is told and not overruled.
            notes.append(
                "UNCENTRED: every column sits far from zero relative to the data's own spread. If "
                "these are RAW SCORES rather than paired differences, the p below answers nothing, "
                "because sign-flipping tests symmetry about ZERO. If they are differences with a "
                "large consistent effect and small variance, this is expected and the p is valid. "
                "This is a HEURISTIC and cannot tell the two apart.")

    stat_kind = "mean" if statistic == "mean" else "max_abs_t"
    obs = _t_stats(mat, stat_kind)
    obs_family = max(abs(v) for v in obs)

    signs_iter, exact, total = _sign_vectors(n, n_draws, seed)
    ge_family = 0
    ge_cand = [0] * k
    used = 0
    for signs in signs_iter:
        flipped = [[s * v for v in row] for s, row in zip(signs, mat)]
        t = _t_stats(flipped, stat_kind)
        used += 1
        if max(abs(v) for v in t) >= obs_family - 1e-12:
            ge_family += 1
        for j in range(k):
            if abs(t[j]) >= abs(obs[j]) - 1e-12:
                ge_cand[j] += 1

    if exact:
        p_family = ge_family / used
        p_cand = [c / used for c in ge_cand]
        notes.append(f"exact enumeration of all 2**{n} = {total} sign vectors; no sampling error")
    else:
        p_family = ge_family / used
        p_cand = [c / used for c in ge_cand]
        notes.append(f"monte-carlo over {used} draws including the observed configuration, so p is "
                     f"bounded below by {1 / used:.2e} and is never 0")

    if n < 5:
        # TWO-SIDED: |t| and |mean| are unchanged by flipping EVERY row, so the all-minus vector
        # always attains the observed statistic alongside the all-plus one. The floor is therefore
        # 2/2**n, not 1/2**n -- stating 1/2**n would understate the smallest resolvable p by 2x.
        floor = (2 / total) if exact else (1 / used)
        notes.append(
            f"only {n} replicate rows: the smallest attainable family-wise p is {floor:.4f} "
            f"(two-sided, so the all-minus sign vector always ties the observed statistic), and an "
            f"effect of ANY size cannot be resolved below it. A non-significant result here is "
            f"UNDERPOWERED, not evidence of absence.")

    return SignFlipResult(
        n_rows=n, n_candidates=k, statistic=stat_kind, observed=obs,
        observed_family=obs_family, family_wise_p=p_family, per_candidate_p=p_cand,
        n_draws=used, exact=exact, seed=None if exact else seed, notes=notes,
    )


def bonferroni_p(per_candidate_p: list[float]) -> float:
    """min(p) * k, capped at 1 -- the independence-assuming comparator.

    Provided ONLY so the two can be shown side by side: Bonferroni assumes the candidates are
    independent, so on a correlated family it over-corrects and throws away real power. It is the
    thing the sign-flip test replaces, not a check on it.
    """
    if not per_candidate_p:
        raise SignFlipRefusal("no per-candidate p values")
    return min(1.0, min(per_candidate_p) * len(per_candidate_p))
