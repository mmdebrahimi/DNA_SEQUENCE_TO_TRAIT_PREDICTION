"""Guards for the masked-pooling pyramid.

The module's whole justification is one identity: PADDING CANNOT CHANGE A SCORE. If that fails, the
fragment arm built on it is invalid, so the identity is tested in both directions --
`test_padding_cannot_change_a_score` asserts it holds under masked pooling, and
`test_an_UNMASKED_pool_WOULD_be_changed_by_padding` proves the first test is not vacuously true of
any pooling scheme. Without the second, a bug that silently disabled masking would leave the first
one green.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.encoder import PYRAMID_WIDTHS, EncoderDataError  # noqa: E402
from dna_decode.glm.masked_pyramid import (  # noqa: E402
    build_masked_pyramid,
    masked_mean_pool,
    one_hot_padded,
    pool_coverage,
    valid_positions,
)

torch = pytest.importorskip("torch", reason="torch not installed")


# ------------------------------------------------------------------ the arithmetic

@pytest.mark.parametrize("true_len,width,expected", [
    (499, 3, 497),      # the longest fragment against the finest branch
    (499, 150, 350),    # the longest fragment against the widest branch
    (150, 150, 1),      # exactly one valid position
    (149, 150, 0),      # one short -> the branch is dead for this row
    (48, 150, 0),       # the shortest fragment cannot feed the widest branch
    (48, 48, 1),
])
def test_valid_positions_is_true_len_minus_width_plus_one(true_len, width, expected):
    assert valid_positions(true_len, width) == expected


def test_valid_positions_refuses_nonsense():
    with pytest.raises(EncoderDataError):
        valid_positions(100, 0)
    with pytest.raises(EncoderDataError):
        valid_positions(-1, 3)


def test_pool_coverage_surfaces_a_branch_that_is_DEAD_on_most_of_the_substrate():
    """On real fragments the 150 bp branch is dead for short rows. That must be visible."""
    lengths = [48, 60, 100, 149, 150, 300, 499]
    cov = pool_coverage(lengths, PYRAMID_WIDTHS)
    assert cov[3]["n_dead_rows"] == 0, "a 3 bp branch is alive on every real fragment"
    assert cov[150]["n_dead_rows"] == 4, "48/60/100/149 cannot feed a 150 bp kernel"
    assert cov[150]["fraction_dead"] == pytest.approx(4 / 7, abs=1e-4)
    assert cov[150]["min_true_len_needed"] == 150


def test_pool_coverage_refuses_an_empty_substrate():
    with pytest.raises(EncoderDataError):
        pool_coverage([])


# ------------------------------------------------------------------ THE load-bearing identity

def _pooled_for_pad(seq: str, pad_to: int, width: int, conv) -> "torch.Tensor":
    x = torch.from_numpy(one_hot_padded([seq], pad_to))
    h = torch.relu(conv(x))
    return masked_mean_pool(h, [len(seq)], width)


@pytest.mark.parametrize("width", list(PYRAMID_WIDTHS))
def test_padding_cannot_change_a_score(width):
    """THE claim. Identical pooled features whatever the sequence is padded to."""
    torch.manual_seed(0)
    conv = torch.nn.Conv1d(4, 8, kernel_size=width)
    seq = "ACGT" * 50  # 200 bp -- long enough to feed every pyramid branch
    tight = _pooled_for_pad(seq, len(seq), width, conv)
    padded = _pooled_for_pad(seq, 499, width, conv)
    over = _pooled_for_pad(seq, 600, width, conv)
    assert torch.allclose(tight, padded, atol=1e-6), "padding to 499 moved the pooled features"
    assert torch.allclose(tight, over, atol=1e-6), "padding to 600 moved the pooled features"


def test_an_UNMASKED_pool_WOULD_be_changed_by_padding():
    """Non-vacuity. Proves the invariance above comes from MASKING, not from the pooling being
    trivially insensitive. Without this, masking could be silently disabled and the test above
    would still pass.
    """
    torch.manual_seed(0)
    conv = torch.nn.Conv1d(4, 8, kernel_size=3)
    seq = "ACGT" * 50
    tight = torch.relu(conv(torch.from_numpy(one_hot_padded([seq], len(seq))))).mean(dim=2)
    padded = torch.relu(conv(torch.from_numpy(one_hot_padded([seq], 499)))).mean(dim=2)
    assert not torch.allclose(tight, padded, atol=1e-6), \
        "an unmasked mean was unaffected by padding, so the invariance test proves nothing"


def test_a_row_too_short_for_a_branch_pools_to_EXACT_zeros():
    """Not NaN (which poisons a batch) and not a fabricated non-zero."""
    torch.manual_seed(0)
    conv = torch.nn.Conv1d(4, 8, kernel_size=150)
    short = "ACGT" * 10  # 40 bp, below the 150 bp kernel
    pooled = _pooled_for_pad(short, 499, 150, conv)
    assert torch.all(pooled == 0.0)
    assert not torch.isnan(pooled).any()


def test_two_rows_pool_independently_so_a_short_row_cannot_contaminate_a_long_one():
    torch.manual_seed(0)
    conv = torch.nn.Conv1d(4, 8, kernel_size=12)
    long_seq, short_seq = "ACGT" * 50, "ACGT" * 5
    x = torch.from_numpy(one_hot_padded([long_seq, short_seq], 499))
    h = torch.relu(conv(x))
    both = masked_mean_pool(h, [len(long_seq), len(short_seq)], 12)
    alone = _pooled_for_pad(long_seq, 499, 12, conv)
    assert torch.allclose(both[0:1], alone, atol=1e-6)


# ------------------------------------------------------------------ refusals

def test_masked_mean_pool_refuses_a_length_count_mismatch():
    feats = torch.zeros(3, 8, 100)
    with pytest.raises(EncoderDataError, match="lengths"):
        masked_mean_pool(feats, [10, 20], 3)


def test_masked_mean_pool_refuses_a_wrong_rank_tensor():
    with pytest.raises(EncoderDataError, match=r"\(N, C, L_out\)"):
        masked_mean_pool(torch.zeros(3, 8), [1, 2, 3], 3)


def test_one_hot_padded_refuses_to_truncate():
    with pytest.raises(EncoderDataError, match="discard measured sequence"):
        one_hot_padded(["A" * 600], 499)


def test_one_hot_padded_encodes_pad_and_N_identically_as_all_zero():
    """The pad must carry no new semantics -- all-zero already means 'no evidence' here."""
    a = one_hot_padded(["ACGT"], 8)
    b = one_hot_padded(["ACGTNNNN"], 8)
    assert (a == b).all()


def test_one_hot_padded_refuses_an_empty_batch():
    with pytest.raises(EncoderDataError):
        one_hot_padded([], 10)


# ------------------------------------------------------------------ the model

def test_the_model_runs_and_is_pad_invariant_END_TO_END():
    """The identity must survive the whole forward pass, not just one branch."""
    torch.manual_seed(0)
    model = build_masked_pyramid(PYRAMID_WIDTHS, channels=4)
    model.eval()
    seqs = ["ACGT" * 50, "ACGTT" * 40, "GGCCAA" * 30]
    lens = [len(s) for s in seqs]
    with torch.no_grad():
        a = model(torch.from_numpy(one_hot_padded(seqs, 499)), lens)
        b = model(torch.from_numpy(one_hot_padded(seqs, 600)), lens)
    assert a.shape == (3,)
    assert torch.allclose(a, b, atol=1e-6)


def test_the_model_refuses_a_duplicate_or_empty_width_set():
    with pytest.raises(EncoderDataError, match="duplicate"):
        build_masked_pyramid((3, 3, 12))
    with pytest.raises(EncoderDataError, match="no kernel widths"):
        build_masked_pyramid(())
