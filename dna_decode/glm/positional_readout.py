"""PURE position-resolved reductions of a `(B, C, L)` latent. No model import, no I/O, no scoring.

The measured sibling arm takes PromoterAtlas's final latent and masked-mean-pools it over positions
to 128 features, reaching median Spearman 0.5076 -- the best arm on the tile-chain rung. Its own
recorded limit says a mean-pool is the weakest reasonable use of the encoder, because averaging over
positions discards WHERE in the window a signal sits. This module holds the reductions that keep
that information, so the arm is a different REDUCTION over an identical forward pass.

**Why every function here is pure and shape-agnostic.** The answer for a synthetic tensor is known
by construction, so each reduction is testable without the 686 MB checkpoint, without D:, and
without a 17-minute forward pass. The model-touching code lives in the scripts; nothing here
imports torch.

**Two invariants this module exists to protect.**

1. **The parity gate.** 150 positions divide evenly by {2, 5, 10}, so for EQUAL-size bins the global
   mean is exactly recoverable from the bin means. `global_mean_from_bin_means` is that inverse, and
   it is what lets a scoring run assert the STORED features reproduce the committed
   `tile_embeddings_native150_n44106.npy` without a second forward pass. A ragged final bin would
   break the recovery quietly, which is why `bin_means`/`bin_maxes` REFUSE a non-dividing
   `n_bins` instead of truncating or padding.

2. **A broken reduction must not read as a clean null.** `genomewide.spearman` maps NaN to 0.0, so a
   reduction that emitted NaN would be scored as "no signal" rather than as an error -- the exact
   failure mode the sibling arm's log1p trap produced on lentiMPRA. Every entry point therefore
   refuses a non-finite input up front.

**What the position axis here does NOT mean.** `Tile` exposes `seq, expression, peak, strand,
category, active, rep1, rep2` -- there is no TSS, no promoter register and no peak-relative offset.
Position is relative to an arbitrary 150 bp window, not to a -35/-10 frame, so a bin index is not a
promoter coordinate. `channel_autocorrelation` exists because of this: spacing is
translation-invariant, so an anchor-free lag reduction targets periodic structure where bins cannot,
and it tests a DIFFERENT hypothesis rather than taking a second shot at the same one.
"""

from __future__ import annotations

import numpy as np

#: The only strand tokens `genomewide.Tile.strand` carries on this substrate.
PLUS, MINUS = "+", "-"


class PositionalReadoutError(ValueError):
    """A reduction was asked for something it refuses to guess at."""


def _check_latent(z, name: str = "z") -> np.ndarray:
    """Coerce to float64 `(B, C, L)` and refuse anything non-finite or mis-shaped."""
    a = np.asarray(z, dtype=np.float64)
    if a.ndim != 3:
        raise PositionalReadoutError(
            "%s must be 3-D (B, C, L); got shape %r with %d dims" % (name, a.shape, a.ndim)
        )
    if 0 in a.shape:
        raise PositionalReadoutError("%s has an empty axis: shape %r" % (name, a.shape))
    if not np.all(np.isfinite(a)):
        n_bad = int((~np.isfinite(a)).sum())
        raise PositionalReadoutError(
            "%s holds %d non-finite value(s); refusing rather than propagating, because spearman "
            "maps NaN to 0.0 and a broken reduction would read as a clean null" % (name, n_bad)
        )
    return a


def _check_bins(length: int, n_bins: int) -> int:
    n = int(n_bins)
    if n < 1:
        raise PositionalReadoutError("n_bins must be >= 1, got %r" % (n_bins,))
    if n > length:
        raise PositionalReadoutError(
            "n_bins %d exceeds the position axis (L=%d); there is no such partition" % (n, length)
        )
    if length % n:
        raise PositionalReadoutError(
            "n_bins %d does not divide L=%d (remainder %d). This function REFUSES a ragged final "
            "bin rather than truncating or padding: unequal bins make the global mean "
            "non-recoverable from the bin means, which would break the parity gate silently."
            % (n, length, length % n)
        )
    return n


def bin_means(z, n_bins: int) -> np.ndarray:
    """Mean over equal-size position bins -> `(B, C, n_bins)`.

    Equal bins are the whole point: `global_mean_from_bin_means` is their exact inverse, so the
    stored features can be checked against the committed mean-pool embedding directly.
    """
    a = _check_latent(z)
    n = _check_bins(a.shape[2], n_bins)
    return a.reshape(a.shape[0], a.shape[1], n, a.shape[2] // n).mean(axis=3)


def bin_maxes(z, n_bins: int) -> np.ndarray:
    """Max over equal-size position bins -> `(B, C, n_bins)`.

    A declared SECONDARY. Unlike the mean it has no inverse to the global mean, so it is never the
    rung the parity gate runs on.
    """
    a = _check_latent(z)
    n = _check_bins(a.shape[2], n_bins)
    return a.reshape(a.shape[0], a.shape[1], n, a.shape[2] // n).max(axis=3)


def global_mean_from_bin_means(bm) -> np.ndarray:
    """The exact inverse for EQUAL bins: `(B, C, n_bins)` -> `(B, C)`.

    This is the parity gate's arithmetic. For equal-size bins the mean of the bin means is the mean
    over all positions, so this must reproduce `z.mean(axis=2)` to float tolerance. It takes the
    STORED bin features, not a model, so the gate tests what will actually be scored.
    """
    a = _check_latent(bm, name="bin_means")
    return a.mean(axis=2)


#: Fewest overlapping positions a lag may leave. Pearson over TWO points is exactly +-1 by
#: construction regardless of the data, so a 2-point overlap is a guaranteed-saturated column
#: carrying no information about the sequence. Three is the principled floor.
MIN_OVERLAP = 3


def channel_autocorrelation(z, lags) -> np.ndarray:
    """Per-channel Pearson autocorrelation across positions, at each lag -> `(B, C * n_lags)`.

    Anchor-free and translation-invariant, which is why it can speak to SPACING on a substrate that
    exposes no promoter register. For lag `k` the two overlapping sub-series are correlated as an
    ordinary Pearson coefficient -- each centred on ITS OWN window, then

        r_k = sum_t u'[t] v'[t] / sqrt(sum_t u'[t]**2 * sum_t v'[t]**2),
        u = x[0 : L-k],  v = x[k : L]

    so `|r_k| <= 1` exactly, and a pure sinusoid of period p gives `r_p = 1.0` and `r_{p/2} = -1.0`.

    **Why per-lag normalisation rather than the biased sample ACF.** The biased form divides every
    lag by the lag-0 energy, which attenuates `r_k` by `(L-k)/L` -- a shrinkage that depends only on
    the lag, not on the sequence. Measured on this geometry (L=150): a pure period-15 cosine gets
    `r_15 = 0.9000` against `r_1 = 0.9014`, so the length bias alone moves the peak off the true
    period. As a FEATURE SET that is worse than cosmetic: it makes a lag-20 column systematically
    smaller than a lag-1 column for reasons unrelated to the tile, which is exactly the kind of
    scale confound a ridge design matrix should not carry. The per-lag form has no such bias and
    every column means the same thing.

    Consequence worth knowing: because there is no attenuation, the period AND ALL ITS HARMONICS in
    range saturate at ~1.0, so `argmax` over lags identifies a MULTIPLE of the period, not the
    period itself. This readout is a feature block for a ridge, not a period estimator.

    **Lag 0 is REFUSED, not silently dropped.** `r_0` is identically 1.0, so including it would add
    one constant column per channel -- features that carry zero information and only inflate the
    design matrix. A caller asking for it has mis-specified the lag set.

    **A channel with no variance on a lag's window yields 0.0 there**, not NaN: the correlation is
    undefined, and 0.0 is the honest "no periodic structure detected" encoding. NaN would be mapped
    to 0.0 by `genomewide.spearman` anyway -- silently, and for the wrong reason.

    Column order is channel-major: `[c0@lag0, c0@lag1, ..., c1@lag0, ...]`, matching
    `flatten_bins` so the two feature blocks concatenate predictably.
    """
    a = _check_latent(z)
    L = a.shape[2]

    ks = [int(k) for k in lags]
    if not ks:
        raise PositionalReadoutError("lags is empty; there is nothing to compute")
    if len(set(ks)) != len(ks):
        dupes = sorted({k for k in ks if ks.count(k) > 1})
        raise PositionalReadoutError(
            "duplicate lag(s) %r would emit identical columns; de-duplicate the lag set" % dupes
        )
    if any(k <= 0 for k in ks):
        bad = sorted(k for k in ks if k <= 0)
        raise PositionalReadoutError(
            "non-positive lag(s) %r refused: r_0 is identically 1.0 (a constant column, zero "
            "information, pure design-matrix inflation) and a negative lag is the same as its "
            "absolute value for an autocorrelation" % bad
        )
    if any(k >= L for k in ks):
        bad = sorted(k for k in ks if k >= L)
        raise PositionalReadoutError(
            "lag(s) %r are >= the position axis (L=%d); no pairs exist at that separation"
            % (bad, L)
        )
    if any(L - k < MIN_OVERLAP for k in ks):
        bad = sorted(k for k in ks if L - k < MIN_OVERLAP)
        raise PositionalReadoutError(
            "lag(s) %r leave fewer than MIN_OVERLAP=%d overlapping positions at L=%d; a Pearson "
            "coefficient over 2 points is exactly +-1 whatever the data, so such a column is "
            "saturated by construction and carries no information" % (bad, MIN_OVERLAP, L)
        )

    out = np.zeros((a.shape[0], a.shape[1], len(ks)), dtype=np.float64)
    for j, k in enumerate(ks):
        u = a[:, :, : L - k]
        v = a[:, :, k:]
        u = u - u.mean(axis=2, keepdims=True)
        v = v - v.mean(axis=2, keepdims=True)
        den = np.sqrt((u * u).sum(axis=2) * (v * v).sum(axis=2))
        live = den > 0
        out[:, :, j] = np.where(live, (u * v).sum(axis=2) / np.where(live, den, 1.0), 0.0)
    return out.reshape(a.shape[0], -1)


def orient(features_by_bin, strands) -> np.ndarray:
    """Reverse the BIN axis for `-` strand rows -> same shape in, same shape out.

    Operates on STORED bin features, so putting both orientations on the ladder costs no second
    forward pass.

    **This is a declared empirical variant, not established biology.** Only 576 of 22,285 minus
    sequences (2.6%) have their reverse complement among the plus sequences, and ZERO of those share
    a peak, so plus and minus tiles are DISTINCT genomic windows rather than two orientations of one
    element. The substrate therefore does not settle whether `seq` is already given 5'->3' on the
    annotated strand. Reversing the bin axis puts both cohorts in a common window orientation IF it
    is not; if it already is, this reduction is the wrong one. Which holds is inferred from a
    performance difference between two rungs -- a 1-bit inference, and flagged as such.

    It does NOT complement the channels, only reverse position order: the channels are learned
    latent dimensions, and there is no known base-complement permutation over them.
    """
    a = _check_latent(features_by_bin, name="features_by_bin")
    s = list(strands)
    if len(s) != a.shape[0]:
        raise PositionalReadoutError(
            "got %d strands for %d rows; a positional misalignment here would reverse the WRONG "
            "tiles and still produce a plausible number" % (len(s), a.shape[0])
        )
    unknown = sorted({t for t in s if t not in (PLUS, MINUS)})
    if unknown:
        raise PositionalReadoutError(
            "unrecognised strand token(s) %r; expected only %r or %r. Defaulting an unknown strand "
            "to either orientation would silently mix the cohorts." % (unknown, PLUS, MINUS)
        )
    out = a.copy()
    neg = np.array([t == MINUS for t in s], dtype=bool)
    if neg.any():
        out[neg] = out[neg][:, :, ::-1]
    return out


def flatten_bins(features_by_bin) -> np.ndarray:
    """`(B, C, n_bins)` -> `(B, C * n_bins)`, channel-major.

    Separate from `bin_means` so `orient` can act on the bin axis before it is folded away.
    """
    a = _check_latent(features_by_bin, name="features_by_bin")
    return a.reshape(a.shape[0], -1)
