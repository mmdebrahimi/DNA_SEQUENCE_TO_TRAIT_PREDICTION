"""Frame logic for the bacterial peak-tile substrate: one definition, shared by three scorers.

WHY THIS MODULE EXISTS, AND WHAT IT DELIBERATELY DOES NOT CONTAIN
-----------------------------------------------------------------
On 2026-11-16 a scratch measurement reported that a single binary split of the peak-tile
label recovers 0.8194 of the 0.941 replicate-derived ceiling, and inferred from that the
pooled Spearman is "overwhelmingly a detection metric". **That statistic was VACUOUS and the
inference is WITHDRAWN.** Spearman is Pearson on RANKS, and a threshold indicator is 1 on
exactly the top *k* ranks, so the joint distribution of (indicator, rank) depends only on
*k* and *n*::

    rho = sqrt( 3k(n-k) / (n^2 - 1) )

At ``n=44,106`` / ``k=14,913`` that closed form equals **0.8194 exactly**, and uniform,
Gaussian, heavy-tailed Pareto and deliberately-bimodal samples all return 0.8194 to four
decimal places. The maximum is 0.866 at ANY median split of ANY continuous label. The number
measured the split proportion and nothing else.

So **there is no best-single-threshold-versus-label statistic in this module, and adding one
would be a regression.** ``shape_summary`` keeps only the two distributional quantities that
are NOT determined by the design (modal-bin share, distinct-value count), which are the same
quantities the project's shipped ``assay_degeneracy`` screen uses.

Its replacement is ``reproducibility_curve``, and the replacement is selection-independent
BY CONSTRUCTION rather than by care. ``GSE144621`` publishes a NESTED replicate design --
four sub-replicate columns ``RNA_exp_sum_{1,2}_{1,2}``, two repeats inside each of two groups,
each group's published summary being exactly the mean of its own two sub-columns. So a curve
can BIN on group 1 and SCORE agreement inside group 2: the binning variable and the scored
pair then share no measurement, and the act of bucketing cannot mechanically force the scored
pair to agree or disagree. Binning on one group and scoring the mirror group is a free
replication.

Even that is not perfectly clean -- binning on group 1 restricts group 1, which correlates
with group 2 -- so ``attenuation_null`` generates synthetic data with the DESIGN HELD FIXED
(same n, same bin count, same global pair correlation, same binning-to-pair correlation, and
reliability CONSTANT across the range by construction) and reports how much trend the
plumbing alone produces. A trend inside that band is range restriction, not a finding.

``dna_decode/glm/genomewide.py`` is NOT edited by this module. The Spearman-Brown arithmetic
and the tie-handling correlation are IMPORTED from it, and ``spearman_brown_ceiling`` is
pinned by test to agree with ``noise_ceiling`` to the decimal rather than forked from it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from dna_decode.glm.genomewide import (
    KEY_COLUMNS,
    SUB_REPLICATE_COLUMNS,
    GenomeWideDataError,
    SubReplicateReadStats,
    noise_ceiling,
    spearman,
    sub_replicate_key,
)
from dna_decode.glm import genomewide as _owner

__all__ = [
    "ActiveFrameError",
    "SubReplicatePair",
    "SubReplicateReadStats",
    "attenuation_null",
    "frame_decomposition",
    "join_sub_replicates",
    "read_sub_replicates",
    "read_sub_replicates_with_stats",
    "reproducibility_curve",
    "restrict_predictions",
    "shape_summary",
    "spearman_brown_ceiling",
    "sub_replicate_key",
    "threshold_audit",
    "tile_key",
]

#: The four sub-replicate columns, as published. Order is load-bearing: group 1 first.
#: ``SUB_REPLICATE_COLUMNS`` / ``KEY_COLUMNS`` / ``sub_replicate_key`` /
#: ``SubReplicateReadStats`` are RE-EXPORTED from ``genomewide``, which owns this
#: substrate. They are imported rather than redefined so there is exactly one spelling
#: of the four column names and one definition of the join key.

#: Below this many rows a per-bin correlation is not worth reporting at all.
MIN_ROWS_PER_BIN = 3

#: Below this many rows a per-bin correlation is reported but flagged as thin.
THIN_ROWS_PER_BIN = 30


class ActiveFrameError(RuntimeError):
    """Raised where a silent answer would be worse than no answer."""


# --------------------------------------------------------------------------------------
# the ceiling arithmetic, ONE definition
# --------------------------------------------------------------------------------------

def spearman_brown_ceiling(r: float) -> dict:
    """The Spearman-Brown step-up, from a BARE correlation.

    ``genomewide.noise_ceiling`` accepts only objects carrying ``.rep1`` / ``.rep2``, so it
    cannot be pointed at an arbitrary column pair -- which is exactly what a curve over the
    four sub-replicate columns needs. This is the same arithmetic reached through a different
    door, and ``tests/test_glm_active_frame.py`` pins the two to agree so the project does not
    end up carrying two drifting copies of its own ceiling formula.

    Pass a RAW correlation. ``noise_ceiling`` computes the step-up from the unrounded value and
    rounds once, so handing it back its own rounded ``single_replicate_agreement`` rounds twice
    and can move the last digit by 1e-4. That is a rounding artifact rather than a disagreement,
    but it is the kind of 1e-4 that gets mistaken for drift.
    """
    rel = 2 * r / (1 + r) if r > -1 else 0.0
    return {
        "single_replicate_agreement": round(float(r), 4),
        "reliability_of_2rep_mean_spearman_brown": round(float(rel), 4),
        "ceiling_for_a_sequence_model": round(math.sqrt(rel), 4) if rel > 0 else 0.0,
        "conditions_on_the_outcome": False,
    }


# --------------------------------------------------------------------------------------
# the narrow reader + the join
# --------------------------------------------------------------------------------------

def tile_key(tile) -> tuple[str, str]:
    """The same key, read off an already-loaded ``Tile``.

    ``Tile.peak`` is built as ``f"{peak_start}_{peak_end}_{strand}"`` and ``Tile.seq`` is
    upper-cased by the loader, so ``(seq, peak)`` reconstructs the verified-unique tuple with
    no re-parse and no edit to ``genomewide``.
    """
    return (tile.seq, tile.peak)


def read_sub_replicates_with_stats(
    path: str | Path, *, category: str = "tile", length: int = 150,
) -> tuple[dict[tuple[str, str], tuple[float, float, float, float]], SubReplicateReadStats]:
    """Call the OWNER of this substrate; translate its error class to this module's.

    The parsing lives in ``genomewide`` because that module owns ``GSE144621`` -- it is where
    ``sniff_delimiter`` solves the space-vs-tab trap, and where ``load_peak_tiles`` fixes the
    row-admission predicate this reader's key set must match EXACTLY for
    ``join_sub_replicates`` to succeed. A second reader here would have to keep two predicates
    in lock-step with nothing enforcing it, which is the failure
    ``tests/test_substrate_loader_ownership.py`` exists to catch.

    The error type is translated so this module's callers keep one exception to catch; the
    refusals themselves (duplicate key, zero rows, missing columns) are the owner's.
    """
    try:
        return _owner.read_sub_replicates_with_stats(path, category=category, length=length)
    except GenomeWideDataError as exc:
        raise ActiveFrameError(str(exc)) from exc


def read_sub_replicates(
    path: str | Path, *, category: str = "tile", length: int = 150,
) -> dict[tuple[str, str], tuple[float, float, float, float]]:
    """``read_sub_replicates_with_stats`` without the accounting. Same refusals."""
    table, _st = read_sub_replicates_with_stats(path, category=category, length=length)
    return table


@dataclass(frozen=True)
class SubReplicatePair:
    """One tile's four sub-replicate readings, plus the two group summaries they build."""

    key: tuple[str, str]
    g1_1: float
    g1_2: float
    g2_1: float
    g2_2: float

    @property
    def g1_mean(self) -> float:
        return (self.g1_1 + self.g1_2) / 2.0

    @property
    def g2_mean(self) -> float:
        return (self.g2_1 + self.g2_2) / 2.0


def join_sub_replicates(tiles, subs) -> list[SubReplicatePair]:
    """Join the sub-replicate table onto loaded tiles, in TILE order, or refuse.

    Raises on a key mismatch in EITHER direction and on a duplicated key on the tile side.
    A partial join is the dangerous failure here: it would drop or mis-pair rows and still
    return a curve that looks fine, which is the shape of defect this whole module exists to
    stop being shipped.
    """
    tile_keys = [tile_key(t) for t in tiles]
    counts: dict[tuple[str, str], int] = {}
    for k in tile_keys:
        counts[k] = counts.get(k, 0) + 1
    dup = [k for k, c in counts.items() if c > 1]
    if dup:
        raise ActiveFrameError(
            f"{len(dup)} tile(s) share a (seq, peak) key, e.g. {dup[0]}. The key is measured "
            f"unique on this substrate; joining on an ambiguous key mis-pairs replicates"
        )

    tset, sset = set(tile_keys), set(subs)
    only_tiles, only_subs = tset - sset, sset - tset
    if only_tiles or only_subs:
        raise ActiveFrameError(
            f"join key sets DIFFER: {len(only_tiles)} key(s) only in tiles, "
            f"{len(only_subs)} only in the sub-replicate table "
            f"(tiles {len(tset)}, subs {len(sset)}). Both sides must be read under the SAME "
            f"filter predicate; a partial join would score a mis-aligned vector silently. "
            f"First tiles-only key: {next(iter(only_tiles), None)}; "
            f"first subs-only key: {next(iter(only_subs), None)}"
        )

    pairs = []
    for k in tile_keys:
        a, b, c, d = subs[k]
        pairs.append(SubReplicatePair(key=k, g1_1=a, g1_2=b, g2_1=c, g2_2=d))
    return pairs


# --------------------------------------------------------------------------------------
# the curve, and the null that makes it readable
# --------------------------------------------------------------------------------------

def _rank_bins(values, n_bins: int) -> list[np.ndarray]:
    """Equal-sized bins by RANK of ``values``, as equal as an integer split allows.

    Binning by rank rather than by value holds the bin SIZES fixed, so the curve's shape
    cannot be driven by the marginal distribution of the binning variable.
    """
    arr = np.asarray(values, dtype=float)
    n = arr.size
    if n_bins < 2:
        raise ActiveFrameError(f"n_bins must be at least 2, got {n_bins}")
    if n < MIN_ROWS_PER_BIN * n_bins:
        raise ActiveFrameError(
            f"{n} rows cannot support {n_bins} bins at {MIN_ROWS_PER_BIN} rows per bin; a "
            f"correlation over fewer rows reports noise as structure"
        )
    order = np.argsort(arr, kind="stable")
    return list(np.array_split(order, n_bins))


def _curve_rows(bin_by, pair_a, pair_b, n_bins: int) -> list[dict]:
    a = np.asarray(pair_a, dtype=float)
    b = np.asarray(pair_b, dtype=float)
    v = np.asarray(bin_by, dtype=float)
    rows = []
    for i, idx in enumerate(_rank_bins(v, n_bins)):
        agree = spearman(a[idx].tolist(), b[idx].tolist())
        ceil = spearman_brown_ceiling(agree)
        rows.append({
            "bin": i,
            "n": int(idx.size),
            "agreement": round(float(agree), 4),
            "ceiling": ceil["ceiling_for_a_sequence_model"],
            "bin_by_range": [round(float(v[idx].min()), 4), round(float(v[idx].max()), 4)],
            "thin": bool(idx.size < THIN_ROWS_PER_BIN),
        })
    return rows


def _trend(rows: list[dict]) -> float:
    """Spearman of bin INDEX against bin agreement: does reliability move across the range?"""
    return float(spearman([r["bin"] for r in rows], [r["agreement"] for r in rows]))


def reproducibility_curve(
    bin_by, pair_a, pair_b, *, n_bins: int = 10, binning_is_independent: bool,
) -> dict:
    """How reliability moves across the label's range, binned on a DIFFERENT measurement.

    ``binning_is_independent`` has no default on purpose: the caller must state whether the
    binning variable shares a measurement with the scored pair, because that single fact
    decides whether the result is about the data or about the binning. Passing ``False``
    returns the per-bin table with ``trend`` set to ``None``, so a confounded curve can be
    inspected but its slope cannot be quoted.

    Refuses outright when the binning array IS one of the scored arrays -- the same object,
    or element-wise identical -- since that is the rigged case, not a declaration error.
    """
    if bin_by is pair_a or bin_by is pair_b:
        raise ActiveFrameError(
            "bin_by is the SAME OBJECT as one of the scored arrays. Sorting by the very "
            "measurement whose agreement is then scored forces the two to be near mirror "
            "images inside each bin, so the curve would be arithmetic, not data"
        )
    v = np.asarray(bin_by, dtype=float)
    a = np.asarray(pair_a, dtype=float)
    b = np.asarray(pair_b, dtype=float)
    if not (v.size == a.size == b.size):
        raise ActiveFrameError(
            f"length mismatch: bin_by {v.size}, pair_a {a.size}, pair_b {b.size}"
        )
    for name, other in (("pair_a", a), ("pair_b", b)):
        if v.size == other.size and np.array_equal(v, other):
            raise ActiveFrameError(
                f"bin_by is ELEMENT-WISE IDENTICAL to {name}. Same rigging as the "
                f"same-object case, reached by a different route"
            )

    rows = _curve_rows(v, a, b, n_bins)
    out = {
        "n": int(v.size),
        "n_bins": n_bins,
        "bins": rows,
        "binning_is_independent": bool(binning_is_independent),
        "pooled_agreement": round(float(spearman(a.tolist(), b.tolist())), 4),
        "conditions_on_the_outcome": not bool(binning_is_independent),
    }
    if binning_is_independent:
        out["trend"] = round(_trend(rows), 4)
        out["reads_as"] = (
            "trend is the Spearman of bin index against within-bin agreement. It is "
            "interpretable ONLY against attenuation_null: binning restricts the binning "
            "variable, which correlates with the scored pair, so some slope arrives from "
            "range restriction alone"
        )
    else:
        out["trend"] = None
        out["reads_as"] = (
            "trend WITHHELD: the caller declared the binning variable is not independent of "
            "the scored pair, so any slope is partly or wholly an artifact of the binning. "
            "The per-bin table is shown for inspection and must not be quoted as a finding"
        )
    return out


def _nearest_psd(c: np.ndarray) -> np.ndarray:
    """Clip negative eigenvalues and restore a unit diagonal."""
    w, v = np.linalg.eigh((c + c.T) / 2.0)
    w = np.clip(w, 1e-8, None)
    out = v @ np.diag(w) @ v.T
    d = np.sqrt(np.clip(np.diag(out), 1e-12, None))
    return out / np.outer(d, d)


def _spearman_to_gaussian_r(rho: float) -> float:
    """Latent Pearson correlation of a Gaussian copula with the given Spearman."""
    return float(2.0 * math.sin(math.pi * max(-1.0, min(1.0, rho)) / 6.0))


def attenuation_null(
    bin_by, pair_a, pair_b, *, n_bins: int = 10, n_draws: int = 200, seed: int = 0,
) -> dict:
    """How much trend the DESIGN alone produces, with reliability constant by construction.

    This is the control the 2026-11-16 retraction demands, and the demand is specific: a
    statistic is only evidence about data if it could have come out differently. So this
    draws synthetic samples that hold the design FIXED -- same n, same bin count, the observed
    marginals exactly, the same global ``pair_a``-``pair_b`` rank correlation and the same
    binning-to-pair correlations -- while making within-bin reliability CONSTANT across the
    range. Any trend the null produces is therefore range restriction and nothing else.

    Reports mean / p95 / MAX over ``n_draws``, in this project's established style: a real
    curve earns its keep by beating the null's MAXIMUM, not its p-value.
    """
    v = np.asarray(bin_by, dtype=float)
    a = np.asarray(pair_a, dtype=float)
    b = np.asarray(pair_b, dtype=float)
    if not (v.size == a.size == b.size):
        raise ActiveFrameError(
            f"length mismatch: bin_by {v.size}, pair_a {a.size}, pair_b {b.size}"
        )
    n = int(v.size)

    rho_ab = float(spearman(a.tolist(), b.tolist()))
    rho_va = float(spearman(v.tolist(), a.tolist()))
    rho_vb = float(spearman(v.tolist(), b.tolist()))

    r_ab = _spearman_to_gaussian_r(rho_ab)
    r_va = _spearman_to_gaussian_r(rho_va)
    r_vb = _spearman_to_gaussian_r(rho_vb)
    corr = _nearest_psd(np.array([
        [1.0, r_va, r_vb],
        [r_va, 1.0, r_ab],
        [r_vb, r_ab, 1.0],
    ]))
    chol = np.linalg.cholesky(corr)

    sorted_v, sorted_a, sorted_b = np.sort(v), np.sort(a), np.sort(b)
    rng = np.random.default_rng(seed)
    trends = []
    for _ in range(int(n_draws)):
        z = rng.standard_normal((n, 3)) @ chol.T
        cols = []
        for j, marg in enumerate((sorted_v, sorted_a, sorted_b)):
            ranks = np.argsort(np.argsort(z[:, j], kind="stable"), kind="stable")
            cols.append(marg[ranks])
        trends.append(_trend(_curve_rows(cols[0], cols[1], cols[2], n_bins)))

    t = np.asarray(trends, dtype=float)
    return {
        "n": n,
        "n_bins": n_bins,
        "n_draws": int(n_draws),
        "seed": int(seed),
        "design_held_fixed": {
            "pooled_spearman_pair": round(rho_ab, 4),
            "bin_by_to_pair_a": round(rho_va, 4),
            "bin_by_to_pair_b": round(rho_vb, 4),
            "marginals": "observed, preserved exactly by rank-mapping each synthetic column",
            "reliability_across_the_range": "CONSTANT by construction -- the copula has one "
                                            "correlation, so the null contains no real trend",
        },
        "trend_null": {
            "mean": round(float(t.mean()), 4),
            "sd": round(float(t.std(ddof=1)) if t.size > 1 else 0.0, 4),
            "p95": round(float(np.percentile(t, 95)), 4),
            "max": round(float(t.max()), 4),
            "min": round(float(t.min()), 4),
        },
        "reads_as": (
            "an observed trend INSIDE this band measures range restriction, not the assay. "
            "Only a trend beyond the null MAX over n_draws is evidence about the data"
        ),
        "conditions_on_the_outcome": False,
    }


# --------------------------------------------------------------------------------------
# distributional shape -- what survives of the retracted helper
# --------------------------------------------------------------------------------------

def shape_summary(values, *, n_bins: int = 20) -> dict:
    """Distributional shape, using ONLY quantities the design does not determine.

    Deliberately has NO best-single-threshold-versus-label field: that quantity is a closed
    form in the split proportion (see the module docstring) and was retracted. Modal-bin share
    and distinct-value count are the project's shipped ``assay_degeneracy`` quantities and do
    carry information about the measurement.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    n = int(arr.size)
    if n == 0:
        raise ActiveFrameError("shape_summary received no finite values")
    lo, hi = float(arr.min()), float(arr.max())
    if hi > lo:
        counts, _edges = np.histogram(arr, bins=n_bins, range=(lo, hi))
        modal_share = float(counts.max()) / n
    else:
        modal_share = 1.0
    distinct = int(np.unique(arr).size)
    return {
        "n": n,
        "n_bins": n_bins,
        "min": round(lo, 4),
        "max": round(hi, 4),
        "median": round(float(np.median(arr)), 4),
        "modal_bin_share": round(modal_share, 4),
        "n_distinct_values": distinct,
        "distinct_value_share": round(distinct / n, 4),
        "conditions_on_the_outcome": False,
        "deliberately_absent": (
            "no best-single-threshold-versus-label Spearman and no share-of-ceiling: that "
            "statistic is rho = sqrt(3k(n-k)/(n^2-1)), a closed form in the split proportion "
            "alone, retracted 2026-11-16"
        ),
    }


# --------------------------------------------------------------------------------------
# the active flag, measured rather than assumed
# --------------------------------------------------------------------------------------

def threshold_audit(tiles, *, active_value: str = "active") -> dict:
    """Does the ``active`` flag COINCIDE with a clean cut on ``expression``? Measure, do not assume.

    Reports the implied cut, the count on each side, and ``overlap_either_side``. On non-zero
    overlap the flag is a LABEL rather than a cut, and callers must treat it as one.
    """
    act = [float(t.expression) for t in tiles if t.active == active_value]
    ina = [float(t.expression) for t in tiles if t.active != active_value]
    if not act or not ina:
        raise ActiveFrameError(
            f"threshold_audit needs both classes; got {len(act)} active / {len(ina)} inactive "
            f"under active_value={active_value!r}"
        )
    lo_active, hi_inactive = min(act), max(ina)
    implied_cut = lo_active
    overlap = sum(1 for x in ina if x > implied_cut) + sum(1 for x in act if x < implied_cut)
    clean = overlap == 0
    return {
        "n": len(act) + len(ina),
        "n_active": len(act),
        "n_inactive": len(ina),
        "implied_cut": round(float(implied_cut), 4),
        "active_range": [round(float(min(act)), 4), round(float(max(act)), 4)],
        "inactive_range": [round(float(min(ina)), 4), round(float(hi_inactive), 4)],
        "active_median": round(float(np.median(act)), 4),
        "inactive_median": round(float(np.median(ina)), 4),
        "overlap_either_side": int(overlap),
        "flag_is_a_clean_cut_on_the_label": bool(clean),
        "conditions_on_the_outcome": False,
        "reads_as": (
            "a clean cut means every active-frame number is RANGE-RESTRICTED on the outcome"
            if clean else
            "overlap is non-zero, so the flag is a LABEL rather than a threshold and callers "
            "must not describe active-frame numbers as a pure range restriction"
        ),
    }


def frame_decomposition(tiles, *, active_value: str = "active") -> dict:
    """Replicate agreement and its implied ceiling, pooled and within each flag class.

    The pooled block is unconditioned. The active and inactive blocks are defined by a
    threshold ON THE LABEL, so each is stamped ``selection_attenuated`` and points at the
    curve: their within-group agreement is pushed apart by the cut itself, which is why the
    0.8849-versus-0.3927 contrast is suggestive rather than decisive and why
    ``reproducibility_curve`` exists.
    """
    act = [t for t in tiles if t.active == active_value]
    ina = [t for t in tiles if t.active != active_value]

    def block(group, conditioned: bool) -> dict:
        nc = dict(noise_ceiling(group))
        nc["conditions_on_the_outcome"] = conditioned
        if conditioned:
            nc["selection_attenuated"] = True
            nc["see_instead"] = (
                "reproducibility_curve(bin_by=group-1 summary, pair=group-2 sub-replicates) "
                "with attenuation_null -- this subset is defined by a cut on the label, so "
                "range restriction attenuates its within-group agreement by construction"
            )
        return nc

    return {
        "pooled": block(list(tiles), conditioned=False),
        "active": block(act, conditioned=True) if act else None,
        "inactive": block(ina, conditioned=True) if ina else None,
        "conditions_on_the_outcome": False,
        "reads_as": (
            "compare the POOLED block across runs. The active and inactive blocks are "
            "selection-attenuated and are not comparable to it, nor to each other"
        ),
    }


# --------------------------------------------------------------------------------------
# mapping stored predictions onto a tile subset
# --------------------------------------------------------------------------------------

def restrict_predictions(tile_indices, preds, tiles, keep) -> tuple[list, list]:
    """Map a stored per-tile prediction vector onto a subset, or fail loudly.

    ``tile_indices`` and ``preds`` are the stored parallel arrays; ``keep`` is a predicate on
    a tile. Raises on a length mismatch and on any index outside ``tiles``, because a shifted
    prediction vector scores perfectly happily and silently measures the wrong thing.
    """
    idx = list(tile_indices)
    pv = list(preds)
    if len(idx) != len(pv):
        raise ActiveFrameError(
            f"stored arrays disagree: {len(idx)} tile_indices vs {len(pv)} predictions. A "
            f"zipped pair of unequal length silently truncates to the shorter one"
        )
    n = len(tiles)
    bad = [i for i in idx if not isinstance(i, (int, np.integer)) or i < 0 or i >= n]
    if bad:
        raise ActiveFrameError(
            f"{len(bad)} stored tile index(es) fall outside the {n} loaded tiles, first {bad[0]}. "
            f"This is a prediction file written against a DIFFERENT tile load; scoring it would "
            f"align predictions to the wrong rows"
        )
    out_y, out_p = [], []
    for i, p in zip(idx, pv):
        t = tiles[int(i)]
        if keep(t):
            out_y.append(float(t.expression))
            out_p.append(float(p))
    return out_y, out_p
