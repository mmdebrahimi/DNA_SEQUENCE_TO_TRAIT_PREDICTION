"""Guards for the G-A conv encoder.

The encoder is being built against a stop gate that PREDICTED it would fail, so the most valuable test
here is the one that separates the two ways it can fail: **the encoder is too weak to learn position at
all** versus **position is not there to be learned on the real substrate**. The planted-signal test is
that separator — a motif at a fixed offset, which the gate's positional-GC proxy could not see and a
linear ridge could not represent. If the encoder cannot learn a planted motif, a later negative on real
tiles says nothing about the substrate and the build is simply broken.

Pure helpers must pass with torch ABSENT (the lazy-import convention); only the training tests
`importorskip` it.
"""
from __future__ import annotations

import random

import pytest

from dna_decode.glm.encoder import (
    BASES,
    DEFAULT_WIDTHS,
    REGISTERED_PROTOCOL,
    EncoderDataError,
    TrainingProtocol,
    conv_output_length,
    count_parameters,
    encoder_feature_width,
    one_hot,
    peak_blocked_val_split,
    require_uniform_length,
    resolve_device,
    validate_widths,
)


class _T:
    """Minimal tile stand-in: the fields `fit_encoder` actually reads."""

    def __init__(self, seq, expression, peak):
        self.seq, self.expression, self.peak = seq, expression, peak


def _planted(n=600, L=150, motif="TTGACA", offset=40, seed=0, noise=0.0):
    """Half the tiles carry a motif at a FIXED offset; expression tracks its presence.

    This is deliberately the signal shape the gate's proxy is blind to: positional GC in 10 bins averages
    over a 15 bp window, and a 6 bp motif at one offset barely moves that average, while a linear model
    over composition cannot express "this pattern, here".

    **`noise=0.0` is load-bearing, and the first version of this fixture got it wrong.** With
    within-class noise, truth has n distinct ranks while any binary detector emits two values, so a
    PERFECT detector is capped by tie structure at Spearman **0.8661** — and the bar had been set at 0.8,
    demanding 92% of an unreachable maximum. The model's 0.71 then read as a failure when it was in fact
    82% of the achievable ceiling. With exactly two expression levels, truth is tied the same way the
    prediction is and a perfect detector scores 1.0, so a 0.8 bar means what it appears to mean.
    `test_the_planted_bar_is_reachable` pins that.
    """
    r = random.Random(seed)
    out = []
    for i in range(n):
        s = [r.choice(BASES) for _ in range(L)]
        has = i % 2 == 0
        if has:
            s[offset:offset + len(motif)] = list(motif)
        e = (3.0 if has else 0.5) + (r.random() * noise if noise else 0.0)
        # peaks are groups of 10 so the inner peak-blocked split has something to block on
        out.append(_T("".join(s), e, f"p{i // 10}"))
    return out


# ---------------------------------------------------------------------------------------------------
# pure helpers -- must work with torch absent
# ---------------------------------------------------------------------------------------------------
def test_conv_output_length_is_the_textbook_formula():
    assert conv_output_length(150, 6) == 145
    assert conv_output_length(150, 75) == 76
    assert conv_output_length(150, 150) == 1
    assert conv_output_length(150, 151) == 0


def test_feature_width_accounts_for_BOTH_pools_per_branch():
    """Two pools (max, mean) per channel per branch. Getting this wrong mis-sizes the head silently."""
    assert encoder_feature_width((6,), channels=32) == 64
    assert encoder_feature_width(DEFAULT_WIDTHS, channels=32) == 192


def test_parameter_count_is_small_by_design():
    """Capacity already HURT once on this substrate (6-mers, 4,096 features, scored worst at 0.2022), so
    the encoder being small is a design commitment and worth pinning rather than drifting."""
    n = count_parameters(DEFAULT_WIDTHS, channels=32)
    # 32*4*(6+20+75) + 3*32 biases + 193 head
    assert n == 32 * 4 * (6 + 20 + 75) + 3 * 32 + 193
    assert n < 15_000, f"encoder grew to {n} parameters; it is supposed to be small"


def test_one_hot_shape_is_N_C_L_which_is_what_Conv1d_requires():
    a = one_hot(["ACGT", "TGCA"])
    assert a.shape == (2, 4, 4)


def test_one_hot_is_a_proper_indicator_in_the_pinned_base_order():
    a = one_hot(["ACGT"])
    assert [a[0, :, p].argmax() for p in range(4)] == [0, 1, 2, 3]
    assert a.sum() == 4.0


def test_an_unknown_base_becomes_an_ALL_ZERO_column_not_a_uniform_prior():
    """An N means 'no evidence of any base'. A 0.25 column would assert equal evidence for all four, and
    defaulting to 'A' would fabricate a base that was never sequenced."""
    a = one_hot(["ANGT"])
    assert a[0, :, 1].sum() == 0.0
    assert a[0, :, 0].sum() == 1.0


def test_lowercase_soft_mask_is_read_not_dropped():
    assert (one_hot(["acgt"]) == one_hot(["ACGT"])).all()


def test_ragged_input_is_REFUSED_and_the_error_names_the_lengths():
    """A silent pad or truncate invents sequence and then returns a plausible number."""
    with pytest.raises(EncoderDataError) as e:
        require_uniform_length(["ACGT", "ACG"])
    assert "ragged" in str(e.value)
    assert "Refusing rather than padding" in str(e.value)


def test_empty_and_zero_length_inputs_are_refused():
    with pytest.raises(EncoderDataError):
        require_uniform_length([])
    with pytest.raises(EncoderDataError):
        require_uniform_length([""])


def test_a_kernel_wider_than_the_sequence_is_refused_not_silently_dropped():
    """A branch with no valid output position would contribute a zero vector, so the model would appear to
    train while one stated resolution was absent."""
    with pytest.raises(EncoderDataError) as e:
        validate_widths((6, 200), 150)
    assert "200" in str(e.value)
    validate_widths((6, 150), 150)  # exactly fits -- allowed


def test_duplicate_and_empty_width_sets_are_refused():
    with pytest.raises(EncoderDataError):
        validate_widths((6, 6), 150)
    with pytest.raises(EncoderDataError):
        validate_widths((), 150)


def test_the_three_registered_widths_cover_three_DIFFERENT_scales():
    """Not a sweep: motif box / motif+spacer / whole-promoter. If two collapsed to the same scale the
    multi-branch arm would be a single-resolution arm wearing three hats."""
    assert DEFAULT_WIDTHS == (6, 20, 75)
    assert len(set(DEFAULT_WIDTHS)) == 3


# ---------------------------------------------------------------------------------------------------
# the inner validation split must not leak
# ---------------------------------------------------------------------------------------------------
def test_the_inner_val_split_is_peak_blocked_and_a_partition():
    peaks = [f"p{i // 10}" for i in range(400)]
    tr, va = peak_blocked_val_split(peaks, val_fraction=0.2, seed=0)
    assert set(tr) | set(va) == set(range(400))
    assert set(tr) & set(va) == set()
    # no PEAK appears on both sides -- the whole point
    assert {peaks[i] for i in tr}.isdisjoint({peaks[i] for i in va})


def test_a_single_peak_cannot_be_split_and_says_so():
    """Tiles inside one peak overlap heavily; splitting them would early-stop against near-duplicates."""
    with pytest.raises(EncoderDataError):
        peak_blocked_val_split(["p0"] * 50, val_fraction=0.2, seed=0)


def test_a_record_with_no_peak_is_REFUSED_rather_than_silently_grouped_together():
    """The grouping variable must be NAMED per substrate. A silent `""` fallback would put every row in
    one group, and a later 'convenience' fallback to a random inner split would leak on a substrate that
    needs a blocked one -- which is exactly the bug class this whole project keeps finding."""
    from dna_decode.glm.encoder import tile_peak_key

    class _NoPeak:
        seq, expression = "ACGT", 1.0

    with pytest.raises(EncoderDataError) as e:
        tile_peak_key(_NoPeak())
    assert "group_key" in str(e.value)
    assert tile_peak_key(_T("ACGT", 1.0, "p3")) == "p3"


def test_the_grid_can_block_on_ELEMENT_identity_instead_of_peak():
    """The designed grid has no peaks; its leakage unit is the element identity. The group_key seam is
    what lets one encoder serve both substrates without either inheriting the other's split."""
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder

    class _P:
        def __init__(self, seq, expression, elem):
            self.seq, self.expression, self.elements = seq, expression, {"e": elem}

    r = random.Random(0)
    data = [_P("".join(r.choice(BASES) for _ in range(60)), r.random(), f"e{i // 8}")
            for i in range(240)]
    tr = [p for p in data if int(p.elements["e"][1:]) < 24]
    te = [p for p in data if int(p.elements["e"][1:]) >= 24]
    out = fit_encoder(tr, te, widths=(6,), seed=0, protocol=TrainingProtocol(max_epochs=2),
                      group_key=lambda p: p.elements["e"])
    assert out["n_train_groups"] > 1
    # and the default grouping must REFUSE this substrate rather than guess
    with pytest.raises(EncoderDataError):
        fit_encoder(tr, te, widths=(6,), seed=0, protocol=TrainingProtocol(max_epochs=2))


def test_the_inner_split_never_consumes_every_peak():
    """val_fraction=1.0 would leave no training rows; it must clamp rather than produce a degenerate fit."""
    tr, va = peak_blocked_val_split([f"p{i}" for i in range(5)], val_fraction=1.0, seed=0)
    assert tr and va


# ---------------------------------------------------------------------------------------------------
# protocol discipline
# ---------------------------------------------------------------------------------------------------
def test_the_protocol_is_FROZEN_so_a_tweak_cannot_be_silent():
    with pytest.raises(Exception):
        REGISTERED_PROTOCOL.lr = 0.5  # type: ignore[misc]


def test_the_registered_protocol_is_recorded_field_by_field():
    d = REGISTERED_PROTOCOL.as_dict()
    assert set(d) == {"optimizer", "lr", "weight_decay", "batch_size", "max_epochs",
                      "early_stop_patience", "channels", "val_fraction"}
    assert d["optimizer"] == "adam"


def test_resolve_device_refuses_an_unknown_string():
    with pytest.raises(EncoderDataError):
        resolve_device("tpu")
    assert resolve_device("cpu") == "cpu"


# ---------------------------------------------------------------------------------------------------
# training -- torch required
# ---------------------------------------------------------------------------------------------------
def test_the_planted_bar_is_reachable():
    """NON-VACUITY FOR THE BAR ITSELF, and the reason it exists is that the bar was mis-specified once.

    A PERFECT detector must score 1.0 on this fixture, or a 0.8 bar is secretly asking for a fraction of
    something unreachable. With within-class noise the same perfect detector caps at 0.8661 -- that is
    measured here too, so the fixture's `noise` parameter cannot be reintroduced without this failing.
    """
    from dna_decode.glm.genomewide import spearman
    data = _planted(n=120, seed=0)
    perfect = [1.0 if i % 2 == 0 else 0.0 for i in range(len(data))]
    assert spearman(perfect, [t.expression for t in data]) == pytest.approx(1.0)

    noisy = _planted(n=120, seed=0, noise=0.2)
    capped = spearman(perfect, [t.expression for t in noisy])
    assert capped < 0.87, "with within-class noise a perfect detector is tie-capped -- the old defect"


def test_the_encoder_CAN_learn_a_motif_at_a_fixed_offset():
    """THE load-bearing test. It separates 'the encoder is too weak' from 'position is not there on the
    real substrate' -- without it, a negative on real tiles is uninterpretable.

    This signal is specifically one the stop gate's proxy is blind to: a 6 bp motif at one offset barely
    moves GC-per-10-bins, and a linear model over composition cannot express 'this pattern, here'.

    Measured: ~0.85 at these settings. Not 1.0, because a 6 bp motif has random partial matches in 150 bp
    of noise and the filter is learned from a few hundred rows -- so the bar sits at 0.8, below what the
    architecture achieves and well above the 0.306 an under-trained run produced.
    """
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    data = _planted(n=1200, seed=1)
    tr = [t for t in data if int(t.peak[1:]) < 96]
    te = [t for t in data if int(t.peak[1:]) >= 96]
    r = fit_encoder(tr, te, widths=(6, 20), seed=0,
                    protocol=TrainingProtocol(max_epochs=80, batch_size=128,
                                              early_stop_patience=1000))
    assert r["spearman"] > 0.8, f"encoder failed to learn a planted motif: {r}"


def test_shuffled_training_labels_collapse_to_about_zero():
    """The validity null. A model whose labels are destroyed must not score; if it does, the pipeline
    fabricates signal and no arm may be graded."""
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    data = _planted(n=600, seed=2)
    tr = [t for t in data if int(t.peak[1:]) < 48]
    te = [t for t in data if int(t.peak[1:]) >= 48]
    r = fit_encoder(tr, te, widths=(6, 20), seed=0, shuffle_labels=True,
                    protocol=TrainingProtocol(max_epochs=25, batch_size=128))
    assert abs(r["spearman"]) < 0.25, f"shuffled-label arm scored {r['spearman']}"
    assert r["shuffle_labels"] is True


def test_the_null_does_not_touch_the_TEST_labels():
    """Non-vacuity for the null: it must permute training targets only, or it measures something else."""
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    data = _planted(n=400, seed=3)
    tr = [t for t in data if int(t.peak[1:]) < 32]
    te = [t for t in data if int(t.peak[1:]) >= 32]
    before = [t.expression for t in te]
    fit_encoder(tr, te, widths=(6,), seed=0, shuffle_labels=True,
                protocol=TrainingProtocol(max_epochs=3))
    assert [t.expression for t in te] == before


def test_CPU_training_is_DETERMINISTIC_for_a_fixed_seed():
    """Same seed, identical predictions. CPU is the determinism reference; CUDA convs are documented as
    potentially nondeterministic and are deliberately NOT claimed."""
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    data = _planted(n=300, seed=4)
    tr = [t for t in data if int(t.peak[1:]) < 24]
    te = [t for t in data if int(t.peak[1:]) >= 24]
    kw = dict(widths=(6,), device="cpu", protocol=TrainingProtocol(max_epochs=4))
    a = fit_encoder(tr, te, seed=7, **kw)
    b = fit_encoder(tr, te, seed=7, **kw)
    assert a["spearman"] == b["spearman"]
    assert a["best_inner_val_mse"] == b["best_inner_val_mse"]


def test_a_different_seed_actually_changes_the_fit():
    """Non-vacuity for the determinism test: if every seed gave the same answer, determinism would be
    trivially true and the multi-seed gate would be averaging one number."""
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    data = _planted(n=300, seed=5)
    tr = [t for t in data if int(t.peak[1:]) < 24]
    te = [t for t in data if int(t.peak[1:]) >= 24]
    kw = dict(widths=(6,), device="cpu", protocol=TrainingProtocol(max_epochs=4))
    a = fit_encoder(tr, te, seed=1, **kw)
    b = fit_encoder(tr, te, seed=2, **kw)
    assert a["best_inner_val_mse"] != b["best_inner_val_mse"]


def test_determinism_flag_is_RESTORED_after_a_fit():
    """A fit that leaves `use_deterministic_algorithms` on would change every later test in the process."""
    torch = pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    before = torch.are_deterministic_algorithms_enabled()
    data = _planted(n=200, seed=6)
    tr = [t for t in data if int(t.peak[1:]) < 16]
    te = [t for t in data if int(t.peak[1:]) >= 16]
    fit_encoder(tr, te, widths=(6,), seed=0, protocol=TrainingProtocol(max_epochs=2))
    assert torch.are_deterministic_algorithms_enabled() == before


def test_empty_train_and_empty_test_are_refused():
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    data = _planted(n=100, seed=7)
    with pytest.raises(EncoderDataError):
        fit_encoder([], data, widths=(6,))
    with pytest.raises(EncoderDataError):
        fit_encoder(data, [], widths=(6,))


def test_train_test_length_mismatch_is_refused():
    """Mixing the 150 bp tiles with the variable-length frag library is a real bug, not a shape to pad."""
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    tr = _planted(n=100, L=150, seed=8)
    te = [_T("A" * 120, 1.0, "q0") for _ in range(20)]
    with pytest.raises(EncoderDataError):
        fit_encoder(tr, te, widths=(6,))


def test_the_returned_dict_records_the_protocol_that_actually_ran():
    """Step 7's `registered_protocol_used` must read a machine-readable field, not trust a human claim."""
    pytest.importorskip("torch")
    from dna_decode.glm.encoder import fit_encoder
    data = _planted(n=200, seed=9)
    tr = [t for t in data if int(t.peak[1:]) < 16]
    te = [t for t in data if int(t.peak[1:]) >= 16]
    p = TrainingProtocol(max_epochs=2, channels=8)
    r = fit_encoder(tr, te, widths=(6, 20), seed=0, protocol=p)
    assert r["protocol"] == p.as_dict()
    assert r["widths"] == [6, 20]
    assert r["n_parameters"] == count_parameters((6, 20), channels=8)
    assert r["device"] == "cpu"
    assert r["n_test"] == len(te)
