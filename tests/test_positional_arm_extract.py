"""Tests for the one-pass extraction runner.

Every test here uses a STUB encoder, so the file is green with no checkpoint, no D: and no
substrate. The stub is deterministic per sequence, which is what makes the resume test meaningful:
an uninterrupted pass and a resumed pass must produce the SAME BYTES, and that can only be asserted
against a reproducible encoder.

The two properties worth the most: a resumed run is byte-identical to an uninterrupted one (this
project has a measured precedent of a mid-run crash leaving a partial cache that looked complete to
its consumer), and `load_block` REFUSES a PARTIAL block (a smoke block reshapes and scores exactly
like a whole one, so the status check is the only thing between a `--limit` run and a published
number).
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

from dna_decode.glm.positional_readout import bin_means, channel_autocorrelation
from scripts.positional_arm_extract import (
    BLOCK_NAME,
    SIDECAR_NAME,
    ExtractError,
    as_bins,
    extract,
    feature_layout,
    load_block,
    rebin,
    reduce_batch,
    slice_block,
)

LAGS = list(range(1, 21))
CHANNELS = 128
LENGTH = 150


def _seq(i: int, n: int = LENGTH) -> str:
    rng = np.random.default_rng(1000 + i)
    return "".join(rng.choice(list("ACGT"), size=n))


def stub_encoder(calls: list | None = None, fail_after_batches: int | None = None):
    """Deterministic per-sequence latent, so a resumed pass is byte-comparable."""
    state = {"batches": 0}

    def encode(seqs: list[str]) -> np.ndarray:
        if fail_after_batches is not None and state["batches"] >= fail_after_batches:
            raise RuntimeError("simulated mid-pass crash")
        state["batches"] += 1
        if calls is not None:
            calls.append(len(seqs))
        out = np.zeros((len(seqs), CHANNELS, len(seqs[0])), dtype=np.float32)
        for b, s in enumerate(seqs):
            seed = int(hashlib.sha256(s.encode()).hexdigest()[:8], 16)
            out[b] = np.random.default_rng(seed).normal(size=(CHANNELS, len(s)))
        return out

    return encode


# --------------------------------------------------------------------------------------------------
# layout arithmetic
# --------------------------------------------------------------------------------------------------
def test_feature_layout_offsets_are_contiguous_and_total_5120():
    lay = feature_layout(20)
    assert lay["bins10_mean"] == {"start": 0, "stop": 1280, "width": 1280}
    assert lay["bins10_max"] == {"start": 1280, "stop": 2560, "width": 1280}
    assert lay["autocorr"] == {"start": 2560, "stop": 5120, "width": 2560}
    assert lay["_total"] == 5120
    # contiguous, no gaps, no overlap
    blocks = [lay[k] for k in ("bins10_mean", "bins10_max", "autocorr")]
    for a, b in zip(blocks, blocks[1:]):
        assert a["stop"] == b["start"]
    assert blocks[-1]["stop"] == lay["_total"]


def test_feature_layout_tracks_the_lag_set():
    assert feature_layout(9)["_total"] == 1280 + 1280 + 128 * 9
    assert feature_layout(9)["autocorr"]["width"] == 1152


def test_feature_layout_refuses_degenerate_inputs():
    for kwargs in ({"n_lags": 0}, {"n_lags": 5, "channels": 0}, {"n_lags": 5, "n_bins": 0}):
        with pytest.raises(ValueError):
            feature_layout(**kwargs)


def test_as_bins_round_trips_and_is_channel_major():
    lay = feature_layout(20)
    rng = np.random.default_rng(4)
    bins = rng.normal(size=(6, CHANNELS, 10))
    flat = bins.reshape(6, -1)
    assert np.array_equal(as_bins(flat, lay), bins)


# --------------------------------------------------------------------------------------------------
# slice_block
# --------------------------------------------------------------------------------------------------
def test_slice_block_returns_the_named_block():
    lay = feature_layout(20)
    block = np.arange(2 * 5120, dtype=np.float32).reshape(2, 5120)
    assert np.array_equal(slice_block(block, lay, "bins10_mean"), block[:, 0:1280])
    assert np.array_equal(slice_block(block, lay, "autocorr"), block[:, 2560:5120])


def test_slice_block_refuses_an_unknown_or_private_name():
    lay = feature_layout(20)
    block = np.zeros((2, 5120), dtype=np.float32)
    with pytest.raises(ExtractError, match="no stored block named"):
        slice_block(block, lay, "bins5_mean")
    with pytest.raises(ExtractError, match="no stored block named"):
        slice_block(block, lay, "_total")


def test_slice_block_refuses_a_block_whose_width_disagrees_with_the_layout():
    lay = feature_layout(20)
    with pytest.raises(ExtractError, match="refusing to slice a mismatched block"):
        slice_block(np.zeros((2, 999)), lay, "bins10_mean")


# --------------------------------------------------------------------------------------------------
# rebin -- why bins2 and bins5 need no stored block of their own
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("target", [1, 2, 5, 10])
def test_rebinning_stored_10_bin_means_equals_binning_the_latent_directly(target):
    """The claim the plan rests the nested ladder on, asserted rather than assumed."""
    rng = np.random.default_rng(31)
    z = rng.normal(size=(5, 8, LENGTH))
    via_stored = rebin(bin_means(z, 10), target)
    direct = bin_means(z, target)
    assert np.allclose(via_stored, direct, rtol=0, atol=1e-12)


def test_rebin_refuses_an_unequal_grouping():
    bins = np.zeros((2, 3, 10))
    for bad in (3, 4, 7, 11, 0):
        with pytest.raises(ExtractError, match="cannot coarsen"):
            rebin(bins, bad)


def test_rebin_refuses_a_wrong_rank_input():
    with pytest.raises(ExtractError, match=r"rebin expects"):
        rebin(np.zeros((2, 10)), 5)


# --------------------------------------------------------------------------------------------------
# the autocorr block's orientation invariance -- a real argument, so it is pinned
# --------------------------------------------------------------------------------------------------
def test_the_autocorr_block_is_invariant_under_reversal_while_the_bin_blocks_are_not():
    """This is WHY the autocorr co-primary is immune to the unresolved orientation convention.

    The set of position pairs separated by k is unchanged when the series is read backwards, so the
    autocorrelation cannot distinguish a tile from its reverse. The bin blocks obviously can, and
    the test asserts both halves -- an invariance claim with no non-invariant control beside it
    would be satisfied by a block of zeros.
    """
    rng = np.random.default_rng(77)
    z = rng.normal(size=(4, 16, LENGTH))
    zr = z[:, :, ::-1]

    lay = feature_layout(20, channels=16)
    f = reduce_batch(z, LAGS)
    fr = reduce_batch(zr, LAGS)

    ac = f[:, lay["autocorr"]["start"]:lay["autocorr"]["stop"]]
    acr = fr[:, lay["autocorr"]["start"]:lay["autocorr"]["stop"]]
    assert np.allclose(ac, acr, atol=1e-5), "autocorr must be reversal-invariant"

    bm = f[:, lay["bins10_mean"]["start"]:lay["bins10_mean"]["stop"]]
    bmr = fr[:, lay["bins10_mean"]["start"]:lay["bins10_mean"]["stop"]]
    assert not np.allclose(bm, bmr, atol=1e-5), "the bin block MUST distinguish the reversal"


def test_the_invariance_holds_for_the_pure_reduction_too():
    rng = np.random.default_rng(78)
    z = rng.normal(size=(3, 5, 80))
    assert np.allclose(channel_autocorrelation(z, [1, 4, 9]),
                       channel_autocorrelation(z[:, :, ::-1], [1, 4, 9]), atol=1e-9)


# --------------------------------------------------------------------------------------------------
# reduce_batch
# --------------------------------------------------------------------------------------------------
def test_reduce_batch_emits_the_declared_width_and_dtype():
    rng = np.random.default_rng(9)
    z = rng.normal(size=(7, CHANNELS, LENGTH))
    out = reduce_batch(z, LAGS)
    assert out.shape == (7, 5120)
    assert out.dtype == np.float32


def test_reduce_batch_blocks_agree_with_the_pure_reductions():
    rng = np.random.default_rng(10)
    z = rng.normal(size=(3, CHANNELS, LENGTH))
    out = reduce_batch(z, LAGS)
    lay = feature_layout(20)
    assert np.allclose(as_bins(out[:, 0:1280], lay), bin_means(z, 10), atol=1e-5)
    assert np.allclose(out[:, 2560:5120], channel_autocorrelation(z, LAGS), atol=1e-5)


# --------------------------------------------------------------------------------------------------
# the run: shape, dtype, sidecar
# --------------------------------------------------------------------------------------------------
def test_a_limit_run_writes_the_declared_shape_dtype_and_sidecar_fields(tmp_path):
    seqs = [_seq(i) for i in range(300)]
    side = extract(encode=stub_encoder(), seqs=seqs, lags=LAGS, block_dir=tmp_path,
                   batch=64, partial=True, checkpoint_sha="deadbeef")
    arr = np.load(tmp_path / BLOCK_NAME, mmap_mode="r")
    assert arr.shape == (300, 5120)
    assert arr.dtype == np.float32
    assert np.all(np.isfinite(arr)), "no row may be left at its zero fill"

    assert side["schema"] == "positional-block-sidecar-v1"
    assert side["geometry"] == {
        "n_rows": 300, "n_features": 5120, "dtype": "float32", "n_bins": 10,
        "lags": LAGS, "channels": CHANNELS, "sequence_length": LENGTH,
    }
    assert side["n_rows_done"] == 300
    assert side["complete"] is True
    assert side["status"] == "PARTIAL"
    assert side["partial_limit_run"] is True
    assert side["checkpoint_sha256"] == "deadbeef"
    assert side["autocorr_orientation_invariant"] is True
    assert side["stored_unoriented"] is True
    assert set(side["layout"]) == {"bins10_mean", "bins10_max", "autocorr"}


def test_a_full_run_is_not_stamped_partial(tmp_path):
    seqs = [_seq(i) for i in range(128)]
    side = extract(encode=stub_encoder(), seqs=seqs, lags=LAGS, block_dir=tmp_path,
                   batch=64, partial=False)
    assert side["status"] == "COMPLETE"
    assert side["partial_limit_run"] is False
    assert side["complete"] is True


def test_the_sidecar_is_rewritten_after_every_batch_so_a_crash_costs_one_batch(tmp_path):
    calls: list[int] = []
    seqs = [_seq(i) for i in range(200)]
    extract(encode=stub_encoder(calls), seqs=seqs, lags=LAGS, block_dir=tmp_path, batch=64)
    assert calls == [64, 64, 64, 8], "batches must be issued in order and sized as asked"
    assert json.loads((tmp_path / SIDECAR_NAME).read_text())["n_rows_done"] == 200


# --------------------------------------------------------------------------------------------------
# resume -- byte-identical to an uninterrupted pass
# --------------------------------------------------------------------------------------------------
def test_resuming_after_an_interruption_reproduces_the_uninterrupted_bytes(tmp_path):
    seqs = [_seq(i) for i in range(320)]

    whole = tmp_path / "whole"
    extract(encode=stub_encoder(), seqs=seqs, lags=LAGS, block_dir=whole, batch=64)
    want = np.array(np.load(whole / BLOCK_NAME, mmap_mode="r"))

    crashed = tmp_path / "crashed"
    with pytest.raises(RuntimeError, match="simulated mid-pass crash"):
        extract(encode=stub_encoder(fail_after_batches=2), seqs=seqs, lags=LAGS,
                block_dir=crashed, batch=64)
    mid = json.loads((crashed / SIDECAR_NAME).read_text())
    assert mid["n_rows_done"] == 128, "two batches should have landed before the crash"
    assert mid["complete"] is False
    assert mid["status"] == "IN_PROGRESS"

    side = extract(encode=stub_encoder(), seqs=seqs, lags=LAGS, block_dir=crashed, batch=64)
    assert side["complete"] is True
    got = np.array(np.load(crashed / BLOCK_NAME, mmap_mode="r"))
    assert np.array_equal(got, want), "a resumed pass must be byte-identical to an uninterrupted one"


def test_a_resumed_run_only_encodes_the_remaining_rows(tmp_path):
    """Non-vacuity: a resume that silently redid the whole pass would also be byte-identical."""
    seqs = [_seq(i) for i in range(320)]
    with pytest.raises(RuntimeError):
        extract(encode=stub_encoder(fail_after_batches=2), seqs=seqs, lags=LAGS,
                block_dir=tmp_path, batch=64)
    calls: list[int] = []
    extract(encode=stub_encoder(calls), seqs=seqs, lags=LAGS, block_dir=tmp_path, batch=64)
    assert sum(calls) == 320 - 128, "only the unfinished rows should be re-encoded"


def test_a_completed_block_is_a_no_op_on_re_run(tmp_path):
    seqs = [_seq(i) for i in range(128)]
    extract(encode=stub_encoder(), seqs=seqs, lags=LAGS, block_dir=tmp_path, batch=64)
    calls: list[int] = []
    extract(encode=stub_encoder(calls), seqs=seqs, lags=LAGS, block_dir=tmp_path, batch=64)
    assert calls == [], "a complete block must not be re-encoded"


# --------------------------------------------------------------------------------------------------
# refusals
# --------------------------------------------------------------------------------------------------
def test_extract_refuses_a_ragged_cohort_rather_than_padding(tmp_path):
    seqs = [_seq(0, 150), _seq(1, 140)]
    with pytest.raises(ExtractError, match="ragged cohort"):
        extract(encode=stub_encoder(), seqs=seqs, lags=LAGS, block_dir=tmp_path)


def test_extract_refuses_zero_sequences(tmp_path):
    with pytest.raises(ExtractError, match="zero sequences"):
        extract(encode=stub_encoder(), seqs=[], lags=LAGS, block_dir=tmp_path)


def test_extract_refuses_to_overwrite_a_block_with_a_different_geometry(tmp_path):
    extract(encode=stub_encoder(), seqs=[_seq(i) for i in range(64)], lags=LAGS,
            block_dir=tmp_path, batch=64)
    # a different lag set changes the width
    with pytest.raises(ExtractError, match="Refusing to overwrite"):
        extract(encode=stub_encoder(), seqs=[_seq(i) for i in range(64)], lags=[1, 2, 3],
                block_dir=tmp_path, batch=64)
    # a different row count changes the height
    with pytest.raises(ExtractError, match="Refusing to overwrite"):
        extract(encode=stub_encoder(), seqs=[_seq(i) for i in range(32)], lags=LAGS,
                block_dir=tmp_path, batch=64)


def test_extract_refuses_an_encoder_that_returns_the_wrong_shape(tmp_path):
    def bad(seqs):
        return np.zeros((len(seqs), CHANNELS))  # 2-D, not (b, C, L)

    with pytest.raises(ExtractError, match="expected"):
        extract(encode=bad, seqs=[_seq(i) for i in range(8)], lags=LAGS, block_dir=tmp_path,
                batch=8)


def test_a_corrupt_sidecar_is_refused_rather_than_guessed_at(tmp_path):
    extract(encode=stub_encoder(), seqs=[_seq(i) for i in range(64)], lags=LAGS,
            block_dir=tmp_path, batch=64)
    (tmp_path / SIDECAR_NAME).write_text("{not json", encoding="utf-8")
    with pytest.raises(ExtractError, match="not parseable JSON"):
        load_block(tmp_path)


# --------------------------------------------------------------------------------------------------
# load_block -- the refusal that makes --limit safe
# --------------------------------------------------------------------------------------------------
def test_load_block_refuses_a_partial_block(tmp_path):
    extract(encode=stub_encoder(), seqs=[_seq(i) for i in range(128)], lags=LAGS,
            block_dir=tmp_path, batch=64, partial=True)
    with pytest.raises(ExtractError, match="stamped PARTIAL"):
        load_block(tmp_path)
    arr, side, lay = load_block(tmp_path, allow_partial=True)
    assert arr.shape == (128, 5120) and side["partial_limit_run"] is True


def test_load_block_refuses_an_incomplete_block(tmp_path):
    with pytest.raises(RuntimeError, match="simulated"):
        extract(encode=stub_encoder(fail_after_batches=1), seqs=[_seq(i) for i in range(320)],
                lags=LAGS, block_dir=tmp_path, batch=64)
    with pytest.raises(ExtractError, match="incomplete"):
        load_block(tmp_path)


def test_load_block_refuses_an_absent_block(tmp_path):
    with pytest.raises(ExtractError, match="no extracted block"):
        load_block(tmp_path / "nothing_here")


def test_load_block_accepts_a_complete_block_and_returns_a_usable_layout(tmp_path):
    seqs = [_seq(i) for i in range(192)]
    extract(encode=stub_encoder(), seqs=seqs, lags=LAGS, block_dir=tmp_path, batch=64)
    arr, side, lay = load_block(tmp_path)
    assert arr.shape == (192, 5120)
    assert lay["_total"] == 5120
    assert side["status"] == "COMPLETE"
    assert slice_block(arr, lay, "bins10_mean").shape == (192, 1280)
    assert as_bins(slice_block(arr, lay, "bins10_mean"), lay).shape == (192, CHANNELS, 10)


def test_load_block_refuses_when_the_file_disagrees_with_its_sidecar(tmp_path):
    seqs = [_seq(i) for i in range(128)]
    extract(encode=stub_encoder(), seqs=seqs, lags=LAGS, block_dir=tmp_path, batch=64)
    # shrink the block on disk, leaving the sidecar claiming 128 rows
    np.save(tmp_path / BLOCK_NAME, np.zeros((64, 5120), dtype=np.float32))
    with pytest.raises(ExtractError, match="disagrees with its sidecar|but its sidecar says"):
        load_block(tmp_path)
