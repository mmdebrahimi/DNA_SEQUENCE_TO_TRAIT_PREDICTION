"""A multi-scale conv pyramid that pools over REAL positions only, for variable-length substrates.

WHY THIS IS A SEPARATE MODULE and not a flag on `encoder.fit_encoder`.

`encoder.py` is the tile-rung encoder and it REFUSES ragged input on purpose:

    A ragged batch has exactly two silent resolutions -- pad or truncate -- and both invent sequence
    that was never measured, then return a plausible number.

That refusal is correct and is left alone. The fragment rung is 48-499 bp over 359 distinct lengths,
so it needs a different contract rather than a loosened one: pad for tensor shape, then make the pad
UNABLE to influence any output. Adding a `lengths=` flag to the shipped encoder would have put a
silent pad path inside the module whose whole job is to refuse one, and every tile-rung caller would
have inherited it.

THE MECHANISM, and it is arithmetic rather than a hope.

A width-`w` convolution with no padding over a length-`L` input emits `L - w + 1` positions. Output
position `p` reads input positions `p .. p+w-1`. So for a sequence whose TRUE length is `t`, padded
to `L`, the outputs that touch only real bases are exactly `p <= t - w`, i.e. the first
`valid = t - w + 1` of them (and none at all when `t < w`). Averaging over just those is therefore
independent of `L`:

    pooled(seq padded to 499) == pooled(seq padded to 600) == pooled(seq unpadded)

That identity is the module's load-bearing claim, so it is a TEST
(`test_padding_cannot_change_a_score`), not a comment. If it ever fails, any arm built on this module
is invalid and must be withdrawn rather than explained -- which is why the arm's frozen bar carries
`INVALID_PAD_LEAKED` as a verdict it can actually reach.

A sequence shorter than a branch's width contributes NOTHING to that branch. That is reported, not
hidden: `pool_coverage` returns the fraction of rows with zero valid positions per branch, so a
width that is dead on most of the substrate is visible instead of silently averaging into zeros.
"""
from __future__ import annotations

from .encoder import BASES, PYRAMID_WIDTHS, EncoderDataError

__all__ = [
    "valid_positions",
    "pool_coverage",
    "masked_mean_pool",
    "one_hot_padded",
    "build_masked_pyramid",
]


def valid_positions(true_len: int, width: int) -> int:
    """How many conv outputs touch ONLY real bases: `max(0, true_len - width + 1)`."""
    if width < 1:
        raise EncoderDataError("kernel width must be >= 1, got %d" % width)
    if true_len < 0:
        raise EncoderDataError("negative true length %d" % true_len)
    return max(0, true_len - width + 1)


def pool_coverage(lengths: list[int], widths=PYRAMID_WIDTHS) -> dict:
    """Per width, the fraction of rows that contribute NOTHING to that branch.

    A branch dead on most of the substrate is a real finding about the width set -- on fragments the
    48 bp and 150 bp branches are the ones at risk. Reported so it cannot pass as signal.
    """
    if not lengths:
        raise EncoderDataError("no lengths given")
    out = {}
    for w in widths:
        dead = sum(1 for t in lengths if valid_positions(t, w) == 0)
        out[int(w)] = {
            "n_dead_rows": dead,
            "fraction_dead": round(dead / len(lengths), 4),
            "min_true_len_needed": int(w),
        }
    return out


def masked_mean_pool(feats, lengths, width: int):
    """Mean over the valid conv outputs only. `feats` is `(N, C, L_out)`; returns `(N, C)`.

    A row with zero valid positions returns zeros for that branch -- the same "no evidence" encoding
    `one_hot` uses for an unknown base, rather than a NaN that would poison the batch or a 0.25-style
    value that would assert evidence nobody measured.
    """
    import torch

    if feats.dim() != 3:
        raise EncoderDataError("expected (N, C, L_out), got shape %s" % (tuple(feats.shape),))
    n, _, l_out = feats.shape
    if len(lengths) != n:
        raise EncoderDataError("got %d lengths for %d rows" % (len(lengths), n))

    valid = torch.tensor([valid_positions(int(t), width) for t in lengths],
                         dtype=torch.long, device=feats.device)
    valid = torch.clamp(valid, max=l_out)
    ar = torch.arange(l_out, device=feats.device).unsqueeze(0)       # (1, L_out)
    mask = (ar < valid.unsqueeze(1)).to(feats.dtype).unsqueeze(1)    # (N, 1, L_out)
    summed = (feats * mask).sum(dim=2)
    denom = mask.sum(dim=2).clamp(min=1.0)
    pooled = summed / denom
    # Rows with no valid position: force exact zeros rather than whatever /1 produced.
    pooled = pooled * (valid > 0).to(feats.dtype).unsqueeze(1)
    return pooled


def one_hot_padded(seqs: list[str], pad_to: int):
    """`(N, 4, pad_to)` float32. Unknown/pad bases are ALL-ZERO columns.

    Deliberately does NOT call `encoder.require_uniform_length`: ragged input is the expected case
    here, and the pad is made inert by `masked_mean_pool` rather than by pretending it is sequence.
    """
    import numpy as np

    if not seqs:
        raise EncoderDataError("no sequences")
    idx = {b: i for i, b in enumerate(BASES)}
    out = np.zeros((len(seqs), len(BASES), pad_to), dtype=np.float32)
    for n, s in enumerate(seqs):
        if len(s) > pad_to:
            raise EncoderDataError(
                "sequence of %d bp exceeds pad_to=%d; truncating would discard measured sequence"
                % (len(s), pad_to))
        for p, ch in enumerate(s.upper()):
            i = idx.get(ch)
            if i is not None:
                out[n, i, p] = 1.0
    return out


def build_masked_pyramid(widths=PYRAMID_WIDTHS, *, channels: int = 32):
    """One conv branch per width, each masked-mean-pooled, concatenated into a linear head.

    Same shape as `encoder.build_encoder` so parameter counts stay comparable across rungs; the only
    difference is that pooling is length-aware.
    """
    import torch
    from torch import nn

    w = tuple(int(x) for x in widths)
    if not w:
        raise EncoderDataError("no kernel widths given")
    if len(set(w)) != len(w):
        raise EncoderDataError("duplicate kernel widths %s" % (w,))
    if any(x < 1 for x in w):
        raise EncoderDataError("kernel widths must be >= 1, got %s" % (w,))

    class MaskedPyramid(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.widths = w
            self.branches = nn.ModuleList(
                [nn.Conv1d(len(BASES), channels, kernel_size=k) for k in w])
            self.head = nn.Linear(channels * len(w), 1)

        def forward(self, x, lengths):
            parts = []
            for conv, k in zip(self.branches, self.widths):
                h = torch.relu(conv(x))
                parts.append(masked_mean_pool(h, lengths, k))
            return self.head(torch.cat(parts, dim=1)).squeeze(-1)

    return MaskedPyramid()
