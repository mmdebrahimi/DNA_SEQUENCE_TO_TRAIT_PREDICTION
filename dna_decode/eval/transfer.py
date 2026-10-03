"""Held-out-GROUP transfer harness: folds, k-shot bookkeeping, power check, leakage audit.

WHAT THIS IS NOT. It contains no de-confounding mathematics. Every statistic comes from
`dna_decode.deconfound` -- the canonical, exported, tested package -- and this module only decides
WHICH rows go where and reports what the primitives silently drop. A second implementation of
`within_group_r2` would create two truth surfaces for the metric the falsifier's replay arm is
compared against, which is the drift this project has paid for repeatedly.

IMPORT POSTURE IS A CONSTRAINT, MEASURED. `dna_decode/eval/__init__.py` is docstring-only and
`eval/genomic_prediction.py` imports sklearn lazily INSIDE its functions (0.002 s module import).
`deconfound/deconfound.py` imports scipy + sklearn at MODULE scope. So every `deconfound` import here
is lazy, inside the function that uses it; a module-scope import would give this module an eager
scipy+sklearn load and slow every test collection that touches `eval/`.

TWO FLOORS, AND THE TIGHTER ONE IS NOT THE OBVIOUS ONE. `MIN_SCORED` (10) is how many points must
remain to score a cell at all. `WITHIN_GROUP_MIN_N` (30) mirrors `within_group_r2`'s own `min_n`: that
function SKIPS any group smaller than it and returns `(nan, 0)`. Measured on 6 groups x 20 members,
`within_group_r2` returns `nan` with zero groups used while `cv_r2` returns a confident-looking
-0.4230 on the identical data. Since the de-confounded cell is the one a PASS verdict reads, 30 is the
binding floor -- which is why `derive_k_grid` uses `max()` of the two rather than `MIN_SCORED` alone.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

#: Minimum scored rows for a k-shot cell to produce a number at all.
MIN_SCORED = 10

#: Mirrors `deconfound.within_group_r2`'s `min_n` default. Passed EXPLICITLY at every call site rather
#: than inherited, so a change on either side cannot silently move this plan's PASS condition.
WITHIN_GROUP_MIN_N = 30

UNDERPOWERED = "UNDERPOWERED_METHODS_NEVER_DIFFERED"
POWERED = "POWERED_METHODS_DIFFER"


class UnassignedMemberError(ValueError):
    """An id carries no group. Silent bucketing warps fold structure -- one giant fold for everything
    unassigned -- so this raises by default, matching `eval/cv.py`'s convention."""


class KShotGridError(ValueError):
    """The group is too small to give even one shot and still satisfy the binding floor."""


@dataclass(frozen=True)
class TransferFold:
    held_out_group: str
    train_ids: tuple[str, ...]
    test_ids: tuple[str, ...]

    @property
    def n_train(self) -> int:
        return len(self.train_ids)

    @property
    def n_test(self) -> int:
        return len(self.test_ids)

    def as_dict(self) -> dict:
        return {"held_out_group": self.held_out_group, "n_train": self.n_train,
                "n_test": self.n_test}


@dataclass
class KShotCell:
    held_out_group: str
    k: int
    score: float | None
    n_scored: int
    status: str                       # "scored" | "unscorable"
    shot_ids: tuple[str, ...] = ()
    scored_ids: tuple[str, ...] = ()
    predictions: tuple[float, ...] = ()

    def as_dict(self) -> dict:
        return {"held_out_group": self.held_out_group, "k": self.k,
                "score": self.score, "n_scored": self.n_scored, "status": self.status,
                "n_shots": len(self.shot_ids)}


@dataclass
class KShotCurve:
    ks: tuple[int, ...]
    cells: list[KShotCell] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {"ks": list(self.ks), "cells": [c.as_dict() for c in self.cells],
                "n_unscorable": sum(1 for c in self.cells if c.status == "unscorable")}

    def cell(self, group: str, k: int) -> KShotCell | None:
        for c in self.cells:
            if c.held_out_group == group and c.k == k:
                return c
        return None


# --- folds ----------------------------------------------------------------------------------------

def held_out_group_folds(group_of: Mapping[str, str], *,
                         allow_unassigned: bool = False) -> list[TransferFold]:
    """One fold per distinct group: that group is held out, everything else trains.

    Disjointness is BY CONSTRUCTION -- membership is a pure function of the group key, so no code path
    can place an id on both sides. Checking two lists for overlap afterwards would only report that it
    did not happen on that run.
    """
    unassigned = [i for i, g in group_of.items() if g is None or str(g).strip() == ""]
    if unassigned and not allow_unassigned:
        raise UnassignedMemberError(
            f"{len(unassigned)} id(s) carry no group (e.g. {unassigned[:3]}). Silent bucketing would "
            f"put them all in one giant fold; pass allow_unassigned=True only after auditing them.")

    resolved = {i: (str(g) if g is not None and str(g).strip() else "__unassigned__")
                for i, g in group_of.items()}
    groups = sorted(set(resolved.values()))
    folds = []
    for g in groups:
        test = tuple(sorted(i for i, gg in resolved.items() if gg == g))
        train = tuple(sorted(i for i, gg in resolved.items() if gg != g))
        folds.append(TransferFold(g, train, test))
    return folds


def fold_report(folds: Sequence[TransferFold]) -> dict:
    """Cluster STRUCTURE, not just a count.

    The recorded TB lesson: 43 clusters hid 97% of the isolates, so a count alone cannot show that a
    partition is dominated by one member. Largest-group and singleton fractions can.
    """
    if not folds:
        return {"n_folds": 0, "n_members": 0, "largest_group_fraction": None,
                "singleton_fraction": None}
    sizes = [f.n_test for f in folds]
    total = sum(sizes)
    return {
        "n_folds": len(folds),
        "n_members": total,
        "largest_group_fraction": round(max(sizes) / total, 4),
        "singleton_fraction": round(sum(1 for s in sizes if s == 1) / len(folds), 4),
        "min_group_size": min(sizes),
        "max_group_size": max(sizes),
    }


def skipped_groups(folds: Sequence[TransferFold], *, min_train: int = 5,
                   min_test: int = 1) -> tuple[str, ...]:
    """Groups that `deconfound.cv_r2` will silently `continue` past (deconfound.py:30-31).

    `cv_r2` drops a group whose train side is under 5 rows or whose test side is empty, leaving NaN
    and quietly scoring a SUBSET of the groups. Reporting which ones is the difference between a
    leave-one-organism-out number and a leave-SOME-organisms-out number wearing the same name.

    REPORTS, never repairs. This duplicates a condition that lives inside another module, so the test
    asserts the prediction against `cv_r2`'s observed NaNs rather than against the literal constant.
    """
    return tuple(f.held_out_group for f in folds
                 if f.n_train < min_train or f.n_test < min_test)


# --- the derived k grid ----------------------------------------------------------------------------

def derive_k_grid(min_group_size: int, *, min_scored: int = MIN_SCORED,
                  within_group_min_n: int = WITHIN_GROUP_MIN_N) -> tuple[int, ...]:
    """`[0] + powers of two` capped by the BINDING floor, which is the larger of the two.

    Raises rather than returning an empty or degenerate grid: a group too small to give one shot and
    still satisfy the within-group floor cannot carry a k-shot curve, and silently returning `(0,)`
    would be the fail-open shape this harness exists to avoid.
    """
    cap = int(min_group_size) - max(int(min_scored), int(within_group_min_n))
    if cap < 1:
        raise KShotGridError(
            f"smallest group has {min_group_size} members; the binding floor is "
            f"max(min_scored={min_scored}, within_group_min_n={within_group_min_n}) so the usable k "
            f"ceiling is {cap} (<1). A k-shot curve is not computable here -- either enlarge the "
            f"groups or state explicitly that the within-group cell will be unscorable.")
    ks = [0]
    p = 1
    while p <= cap:
        ks.append(p)
        p *= 2
    return tuple(ks)


# --- the k-shot curve ------------------------------------------------------------------------------

def _shot_ids(test_ids: Sequence[str], k: int, seed: int) -> tuple[str, ...]:
    """Deterministic k-of-n draw. Keyed on (seed, group) only -- NOT on the method -- so every method
    compared sees the identical shots and the comparison is PAIRED. A difference of medians over
    different items is not a lift."""
    if k == 0:
        return ()
    rng = np.random.default_rng(abs(hash((seed, tuple(test_ids)))) % (2 ** 32))
    idx = rng.permutation(len(test_ids))[:k]
    return tuple(sorted(str(test_ids[i]) for i in idx))


def k_shot_curve(X, y, ids: Sequence[str], group_of: Mapping[str, str], *,
                 fit_predict, ks: tuple[int, ...] | None = None, seed: int = 0,
                 allow_unassigned: bool = False) -> KShotCurve:
    """Score each (held-out group, k) cell. `ids[i]` names row `i` of `X`/`y`.

    The k shots are drawn FROM the held-out group, added to train, and EXCLUDED from the scored set by
    construction -- so a shot can never be scored as if it were held out.

    `ks=None` derives the grid from the smallest group (see `derive_k_grid`).
    """
    from dna_decode.deconfound import r2          # lazy: keeps scipy/sklearn off the import path

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    ids = [str(i) for i in ids]
    if not (len(ids) == len(y) == X.shape[0]):
        raise ValueError(f"ids/X/y length mismatch: {len(ids)}, {X.shape[0]}, {len(y)}")
    pos = {i: n for n, i in enumerate(ids)}

    folds = held_out_group_folds({i: group_of[i] for i in ids},
                                 allow_unassigned=allow_unassigned)
    if ks is None:
        ks = derive_k_grid(min(f.n_test for f in folds))

    curve = KShotCurve(ks=tuple(ks))
    for f in folds:
        for k in ks:
            if k >= f.n_test:
                raise KShotGridError(
                    f"k={k} >= group {f.held_out_group!r} size {f.n_test}: no rows would remain to "
                    f"score. Use derive_k_grid() or pass a smaller ks.")
            shots = _shot_ids(f.test_ids, k, seed)
            scored = tuple(i for i in f.test_ids if i not in set(shots))
            assert not (set(shots) & set(scored)), "shot leaked into the scored set"

            if len(scored) < MIN_SCORED:
                curve.cells.append(KShotCell(f.held_out_group, k, None, len(scored), "unscorable",
                                             shots, scored))
                continue

            tr = [pos[i] for i in f.train_ids] + [pos[i] for i in shots]
            te = [pos[i] for i in scored]
            yhat = np.asarray(fit_predict(X[tr], y[tr], X[te]), dtype=float)
            curve.cells.append(KShotCell(
                f.held_out_group, k, float(r2(y[te], yhat)), len(scored), "scored",
                shots, scored, tuple(float(v) for v in yhat)))
    return curve


def power_check(curve_a: KShotCurve, curve_b: KShotCurve) -> str:
    """Did the two methods ever actually produce different predictions?

    MUST be consulted BEFORE any effect size. A run where the methods never differ necessarily has a
    zero gain, so testing the gain first would publish an unpowered run as a refutation -- the exact
    trap the serotype lineage-disjoint run had to guard against.
    """
    for a in curve_a.cells:
        b = curve_b.cell(a.held_out_group, a.k)
        if b is None:
            continue
        if a.status != b.status:
            return POWERED
        if a.status == "scored" and a.predictions != b.predictions:
            return POWERED
    return UNDERPOWERED
