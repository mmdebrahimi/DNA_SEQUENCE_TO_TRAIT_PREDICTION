"""Tests for the shipped-expression-head module.

Three of these exist to pin traps rather than behaviour, and they are the ones worth reading:

* `padding='same'` on an EVEN kernel is ASYMMETRIC. A symmetric-looking split produces a tensor of
  the right shape and wrong values, so nothing raises -- it just moves the arm's number.
* The two windowing routes are NOT the same function. A test asserting they agree everywhere would
  be asserting a bug, so the assertion here is the real one: identical in a window's interior,
  different at its edges.
* The base's state-dictionary ELEMENT count (1,847,492) is not its PARAMETER count (1,847,428). The
  repo enforces the second as a hard refusal, so a probe that sums the state dictionary reports a
  checkpoint mismatch that does not exist.

Everything except the clone-touching tests runs on synthetic tensors: no checkpoint, no D:, no
forward pass.
"""

from __future__ import annotations

import pathlib

import pytest

torch = pytest.importorskip("torch", reason="atlas_heads' forward routes need torch")

from dna_decode.glm import atlas_heads as AH
from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO, EXPECTED_PARAMS

EXPR_DIR = pathlib.Path(DEFAULT_PA_REPO) / "trained_weights/expression_prediction"
BASE_CK = pathlib.Path(DEFAULT_PA_REPO) / "trained_weights/base_model/promoteratlas-base.pt"
_CLONE = EXPR_DIR.is_dir() and BASE_CK.is_file()
needs_clone = pytest.mark.skipif(not _CLONE, reason="PromoterAtlas clone absent (volatile temp dir)")


# --------------------------------------------------------------------------------------------------
# synthetic head weights
# --------------------------------------------------------------------------------------------------
def _head_sd(*, seed: int = 11) -> dict:
    """A head-shaped state dict with small random weights. Shapes are the real ones."""
    g = torch.Generator().manual_seed(seed)
    c = AH.LATENT_CHANNELS
    return {
        "conv1.weight": torch.randn(c, c, AH.CONV_KERNEL, generator=g) * 0.05,
        "conv1.bias": torch.randn(c, generator=g) * 0.05,
        "conv2.weight": torch.randn(c, c, AH.CONV_KERNEL, generator=g) * 0.05,
        "conv2.bias": torch.randn(c, generator=g) * 0.05,
        "fc1.weight": torch.randn(*AH.FC1_SHAPE, generator=g) * 0.01,
        "fc1.bias": torch.randn(AH.FC1_SHAPE[0], generator=g) * 0.01,
        "fc2.weight": torch.randn(1, AH.FC1_SHAPE[0], generator=g) * 0.01,
        "fc2.bias": torch.randn(1, generator=g) * 0.01,
    }


def _latent(b: int = 2, length: int = 150, *, seed: int = 5):
    g = torch.Generator().manual_seed(seed)
    return torch.randn(b, AH.LATENT_CHANNELS, length, generator=g)


# --------------------------------------------------------------------------------------------------
# TRAP 1 -- the asymmetric 'same' pad
# --------------------------------------------------------------------------------------------------
def test_the_same_pad_split_is_asymmetric_and_the_pinned_one_is_the_right_way_round():
    """If this ever flips, every head output changes silently and nothing raises."""
    import torch.nn.functional as F

    x = torch.zeros(1, 1, 9)
    x[0, 0, 0] = 1.0                                  # impulse, so taps are distinguishable
    w = torch.arange(1.0, 1.0 + AH.CONV_KERNEL).reshape(1, 1, AH.CONV_KERNEL)
    same = F.conv1d(x, w, padding="same")

    pinned = F.conv1d(F.pad(x, (AH.SAME_PAD_LEFT, AH.SAME_PAD_RIGHT)), w)
    flipped = F.conv1d(F.pad(x, (AH.SAME_PAD_RIGHT, AH.SAME_PAD_LEFT)), w)

    assert torch.equal(pinned, same), "the pinned (left, right) split must equal padding='same'"
    assert not torch.equal(flipped, same), (
        "the reversed split must NOT equal padding='same' -- if it did, this test could not detect "
        "getting the asymmetry backwards"
    )
    assert (AH.SAME_PAD_LEFT, AH.SAME_PAD_RIGHT) == (3, 4)
    assert AH.SAME_PAD_LEFT + AH.SAME_PAD_RIGHT == AH.CONV_KERNEL - 1


def test_the_trunk_preserves_length_so_a_window_slice_means_what_it_says():
    sd = _head_sd()
    for length in (86, 120, 150):
        out = AH._trunk(_latent(1, length), sd)
        assert out.shape[-1] == length, "trunk must be length-preserving at L=%d" % length


# --------------------------------------------------------------------------------------------------
# window arithmetic
# --------------------------------------------------------------------------------------------------
def test_window_starts_on_the_real_geometry_is_65_windows_0_through_64():
    starts = AH.window_starts(150, AH.HEAD_SEQ_LENGTH)
    assert len(starts) == 65
    assert starts[0] == 0 and starts[-1] == 64
    assert AH.n_windows(150) == 65


def test_the_geometry_arithmetic_is_internally_consistent():
    assert AH.FC1_SHAPE == (64, 11008)
    assert AH.FC1_SHAPE[1] == AH.LATENT_CHANNELS * AH.HEAD_SEQ_LENGTH == 11008
    assert 150 - AH.HEAD_SEQ_LENGTH + 1 == 65
    # the 150-position latent would flatten to this, which is why it cannot be fed directly
    assert AH.LATENT_CHANNELS * 150 == 19200 != AH.FC1_SHAPE[1]


def test_a_window_wider_than_the_latent_is_REFUSED_not_padded():
    with pytest.raises(AH.AtlasHeadsError, match="does not fit"):
        AH.window_starts(50, AH.HEAD_SEQ_LENGTH)


def test_an_exactly_fitting_window_gives_one_window():
    assert AH.window_starts(86, 86) == (0,)


@pytest.mark.parametrize("length,width", [(0, 86), (150, 0), (-1, 86), (150, -3)])
def test_nonpositive_geometry_is_refused(length, width):
    with pytest.raises(AH.AtlasHeadsError, match="positive"):
        AH.window_starts(length, width)


def test_non_integer_geometry_is_refused_rather_than_coerced():
    with pytest.raises(AH.AtlasHeadsError, match="must be ints"):
        AH.window_starts(150.0, 86)


# --------------------------------------------------------------------------------------------------
# the two routes -- THE load-bearing test
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("s", [0, 30, 64])
def test_the_differing_columns_are_EXACTLY_the_closed_form_prediction(s):
    """The `padding='same'` consequence as a mechanical claim, not an adjective.

    Run in float64 deliberately. In float32 the two routes differ at every column by order 1e-7 --
    `F.conv1d` blocks a 150-position input differently from an 86-position one and `MaxPool1d`
    amplifies it at argmax ties -- so float32 cannot distinguish the real edge effect from
    arithmetic noise. In float64 the differing set collapses to the derived boundary columns.
    """
    sd = {k: v.double() for k, v in _head_sd().items()}
    z = _latent(1, 150).double()
    w = AH.HEAD_SEQ_LENGTH

    cheap = AH._trunk(z, sd)[:, :, s:s + w]
    true_ = AH._trunk(z[:, :, s:s + w], sd)
    assert cheap.shape == true_.shape

    per_col = (cheap - true_).abs().amax(dim=(0, 1))
    observed = {int(c) for c in (per_col > 1e-9).nonzero().flatten().tolist()}
    predicted = set(AH.contaminated_columns(s, 150, w))

    assert observed == predicted, (
        "extra %r / missing %r -- the closed form in `contaminated_columns` no longer describes "
        "what the functions do" % (sorted(observed - predicted)[:8], sorted(predicted - observed)[:8])
    )
    assert observed, "the routes MUST differ somewhere, or the comparison is vacuous"


def test_the_boundary_widths_follow_from_the_pad_split_and_are_NOT_equal():
    """The right edge is wider than the left precisely because the `same` split is asymmetric."""
    assert AH.EDGE_COLS_LEFT == 2 * AH.SAME_PAD_LEFT + AH.POOL_PADDING == 9
    assert AH.EDGE_COLS_RIGHT == 2 * AH.SAME_PAD_RIGHT + AH.POOL_PADDING == 11
    assert AH.EDGE_COLS_LEFT != AH.EDGE_COLS_RIGHT


def test_a_SHARED_boundary_contaminates_nothing_which_is_what_makes_the_offsets_differ():
    w = AH.HEAD_SEQ_LENGTH
    first = AH.contaminated_columns(0, 150, w)        # left edge shared with the latent's own
    last = AH.contaminated_columns(64, 150, w)        # right edge shared
    middle = AH.contaminated_columns(30, 150, w)      # neither shared

    assert len(first) == AH.EDGE_COLS_RIGHT and min(first) == w - AH.EDGE_COLS_RIGHT
    assert len(last) == AH.EDGE_COLS_LEFT and max(last) == AH.EDGE_COLS_LEFT - 1
    assert len(middle) == AH.EDGE_COLS_LEFT + AH.EDGE_COLS_RIGHT
    assert first | last == middle


def test_contaminated_columns_is_EMPTY_when_the_window_is_the_whole_latent():
    """Both boundaries shared, so the routes must coincide -- the degenerate case, predicted."""
    assert AH.contaminated_columns(0, AH.HEAD_SEQ_LENGTH, AH.HEAD_SEQ_LENGTH) == frozenset()


def test_contaminated_columns_refuses_an_offset_that_is_not_a_window():
    with pytest.raises(AH.AtlasHeadsError, match="not a valid window offset"):
        AH.contaminated_columns(65, 150, AH.HEAD_SEQ_LENGTH)


def test_in_float32_the_routes_differ_EVERYWHERE_which_is_the_noise_floor_the_bar_must_clear():
    """TRAP 4, pinned: a tolerance set below ~1e-7 would be measuring arithmetic, not geometry."""
    sd = _head_sd()
    z = _latent(1, 150)
    cheap = AH._trunk(z, sd)[:, :, 30:30 + AH.HEAD_SEQ_LENGTH]
    true_ = AH._trunk(z[:, :, 30:30 + AH.HEAD_SEQ_LENGTH], sd)
    per_col = (cheap - true_).abs().amax(dim=(0, 1))

    n_diff = int((per_col > 1e-9).sum())
    assert n_diff > len(AH.contaminated_columns(30, 150, AH.HEAD_SEQ_LENGTH)), (
        "float32 is expected to differ at MORE columns than the closed form predicts; if it no "
        "longer does, the float64 run is no longer necessary and this test should be revisited"
    )
    # the noise is small relative to the edge effect, which is why float64 separates them cleanly
    interior = sorted(set(range(AH.HEAD_SEQ_LENGTH))
                      - set(AH.contaminated_columns(30, 150, AH.HEAD_SEQ_LENGTH)))
    assert float(per_col[interior].max()) < 1e-3


def test_both_routes_return_one_value_per_window_and_they_are_not_identical():
    sd = _head_sd()
    z = _latent(2, 150)
    cheap = AH.head_forward_full(z, sd)
    true_ = AH.head_forward_per_window(z, sd)
    assert cheap.shape == true_.shape == (2, 65)
    assert not torch.allclose(cheap, true_), (
        "if the two routes produced identical window scores there would be nothing for the bar to "
        "choose between, and the probe's measurement would be meaningless"
    )


def test_on_a_latent_exactly_one_window_wide_the_routes_COINCIDE():
    """The degenerate case that proves the divergence is about slicing, not about the code paths."""
    sd = _head_sd()
    z = _latent(1, AH.HEAD_SEQ_LENGTH)
    assert torch.allclose(AH.head_forward_full(z, sd), AH.head_forward_per_window(z, sd),
                          atol=0.0, rtol=0.0)


def test_a_trunk_evaluated_on_a_slice_is_what_the_true_route_uses():
    """Pins the true route to the 'head handed an 86 bp sequence' semantics it claims."""
    sd = _head_sd()
    z = _latent(1, 150)
    s = 7
    direct = AH._classify(AH._trunk(z[:, :, s:s + AH.HEAD_SEQ_LENGTH], sd), sd)
    assert torch.allclose(AH.head_forward_per_window(z, sd)[:, s], direct, atol=0.0, rtol=0.0)


# --------------------------------------------------------------------------------------------------
# TRAP 3 -- non-finite input must refuse, never score as a clean null
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_a_non_finite_latent_is_refused_by_both_routes(bad):
    sd = _head_sd()
    z = _latent(1, 150)
    z[0, 0, 4] = bad
    for fn in (AH.head_forward_full, AH.head_forward_per_window):
        with pytest.raises(AH.AtlasHeadsError, match="non-finite"):
            fn(z, sd)


def test_reduce_windows_refuses_non_finite_values():
    v = torch.zeros(2, 5)
    v[1, 2] = float("nan")
    with pytest.raises(AH.AtlasHeadsError, match="non-finite"):
        AH.reduce_windows(v, "mean")


# --------------------------------------------------------------------------------------------------
# reductions
# --------------------------------------------------------------------------------------------------
def test_reduce_windows_computes_the_declared_reductions():
    v = torch.tensor([[1.0, 2.0, 3.0], [-1.0, 0.0, 5.0]])
    assert torch.allclose(AH.reduce_windows(v, "mean"), torch.tensor([2.0, 4.0 / 3.0]))
    assert torch.allclose(AH.reduce_windows(v, "max"), torch.tensor([3.0, 5.0]))


def test_an_undeclared_reduction_is_refused_rather_than_defaulted():
    """A silent default would decide a DECLARED PRIMARY by accident."""
    with pytest.raises(AH.AtlasHeadsError, match="unknown reduction"):
        AH.reduce_windows(torch.zeros(1, 3), "median")


def test_reduce_windows_refuses_the_wrong_rank():
    with pytest.raises(AH.AtlasHeadsError, match=r"expected \(B, n_windows\)"):
        AH.reduce_windows(torch.zeros(2, 3, 4), "mean")


# --------------------------------------------------------------------------------------------------
# TRAP 2 -- element counts versus parameter counts
# --------------------------------------------------------------------------------------------------
def test_the_pinned_element_counts_reconcile_with_each_other():
    assert AH.EXPECTED_HEAD_ONLY + AH.BASE_STATE_DICT_ELEMENTS == AH.EXPECTED_HEAD_ELEMENTS


def test_the_element_count_differs_from_the_PARAMETER_count_by_exactly_the_rotary_buffers():
    """8 blocks x an 8-element `rotary_emb.inv_freq` buffer = 64 elements `parameters()` omits."""
    assert AH.BASE_STATE_DICT_ELEMENTS - EXPECTED_PARAMS == 64
    assert EXPECTED_PARAMS == 1_847_428 and AH.BASE_STATE_DICT_ELEMENTS == 1_847_492


# --------------------------------------------------------------------------------------------------
# endpoint derivation
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("name,endpoint", [
    ("promoteratlas-trspred-hossain2020.pt", AH.TRANSCRIPTION),
    ("promoteratlas-trspred-lafleur2022.pt", AH.TRANSCRIPTION),
    ("promoteratlas-trspred-urtecho2018.pt", AH.TRANSCRIPTION),
    ("promoteratlas-trspred-yu2021.pt", AH.TRANSCRIPTION),
    ("promoteratlas-tslpred-kosuri2013.pt", AH.TRANSLATION),
])
def test_the_endpoint_is_derived_from_the_clones_own_filename(name, endpoint):
    assert AH.head_endpoint(name) == endpoint


def test_the_kosuri_head_is_the_only_TRANSLATION_endpoint():
    """Its endpoint does not match a promoter-activity rung; the bar records that as a limit."""
    names = ["promoteratlas-trspred-%s.pt" % s
             for s in ("hossain2020", "lafleur2022", "urtecho2018", "yu2021")]
    names.append("promoteratlas-tslpred-kosuri2013.pt")
    tsl = [n for n in names if AH.head_endpoint(n) == AH.TRANSLATION]
    assert tsl == ["promoteratlas-tslpred-kosuri2013.pt"]


def test_an_unrecognised_endpoint_token_is_refused_not_guessed():
    with pytest.raises(AH.AtlasHeadsError, match="cannot derive an endpoint"):
        AH.head_endpoint("promoteratlas-somethingelse-foo2030.pt")


# --------------------------------------------------------------------------------------------------
# load_head refusals -- driven with real torch.save'd checkpoints so the static audit runs too
# --------------------------------------------------------------------------------------------------
def _base_sd(*, seed: int = 3) -> dict:
    g = torch.Generator().manual_seed(seed)
    return {"blocks.0.feed_forward_up.weight": torch.randn(4, 4, generator=g),
            "to_nucleotides.bias": torch.randn(4, generator=g)}


def _write_head(tmp_path: pathlib.Path, base: dict, *, name="promoteratlas-trspred-x2020.pt",
                mutate=None, norm=None) -> pathlib.Path:
    sd = _head_sd()
    sd.update({"base_model." + k: v.clone() for k, v in base.items()})
    if mutate:
        mutate(sd)
    p = tmp_path / name
    torch.save({"model_state_dict": sd, "normalization_params": norm or {"mean": 0.0, "std": 1.0}}, p)
    return p


def test_load_head_accepts_a_wellformed_synthetic_head(tmp_path):
    base = _base_sd()
    h = AH.load_head(_write_head(tmp_path, base), base_state_dict=base)
    assert h.endpoint == AH.TRANSCRIPTION
    assert h.predicts_transcription and not h.log_transform
    assert h.normalization == {"mean": 0.0, "std": 1.0}


def test_load_head_REFUSES_a_head_whose_base_weights_differ_by_ONE_tensor(tmp_path):
    """Drift makes this a FINE-TUNED-encoder head, a different claim and a different arm."""
    base = _base_sd()

    def drift(sd):
        sd["base_model.to_nucleotides.bias"] = sd["base_model.to_nucleotides.bias"] + 1e-3

    p = _write_head(tmp_path, base, mutate=drift)
    with pytest.raises(AH.AtlasHeadsError, match="NOT byte-identical"):
        AH.load_head(p, base_state_dict=base)


def test_the_drift_refusal_is_non_vacuous_the_same_head_loads_without_the_drift(tmp_path):
    base = _base_sd()
    assert AH.load_head(_write_head(tmp_path, base), base_state_dict=base).name


def test_load_head_REFUSES_a_wrong_fc1_shape(tmp_path):
    base = _base_sd()

    def shrink(sd):
        sd["fc1.weight"] = torch.zeros(64, AH.LATENT_CHANNELS * 150)

    with pytest.raises(AH.AtlasHeadsError, match="requires"):
        AH.load_head(_write_head(tmp_path, base, mutate=shrink), base_state_dict=base)


def test_load_head_REFUSES_an_unexpected_head_only_key_set(tmp_path):
    base = _base_sd()

    def extra(sd):
        sd["fc3.weight"] = torch.zeros(1, 1)

    with pytest.raises(AH.AtlasHeadsError, match="head-only keys"):
        AH.load_head(_write_head(tmp_path, base, mutate=extra), base_state_dict=base)


def test_load_head_REFUSES_a_base_tensor_the_frozen_base_does_not_have(tmp_path):
    base = _base_sd()
    extra_base = dict(base)
    extra_base["blocks.9.invented"] = torch.zeros(2)
    p = _write_head(tmp_path, extra_base)
    with pytest.raises(AH.AtlasHeadsError, match="NOT byte-identical"):
        AH.load_head(p, base_state_dict=base)


def test_load_head_reports_a_log_transform_only_when_the_checkpoint_says_so(tmp_path):
    base = _base_sd()
    p = _write_head(tmp_path, base, name="promoteratlas-tslpred-k2013.pt",
                    norm={"mean": 8.0, "std": 1.2, "log_transform": True})
    h = AH.load_head(p, base_state_dict=base)
    assert h.endpoint == AH.TRANSLATION and h.log_transform is True


def test_load_head_refuses_an_absent_path(tmp_path):
    with pytest.raises(AH.AtlasHeadsError, match="absent"):
        AH.load_head(tmp_path / "promoteratlas-trspred-nope.pt", base_state_dict={})


def test_head_paths_refuses_an_absent_or_empty_directory(tmp_path):
    with pytest.raises(AH.AtlasHeadsError, match="directory absent"):
        AH.head_paths(tmp_path / "nothing")
    (tmp_path / "empty").mkdir()
    with pytest.raises(AH.AtlasHeadsError, match="no .pt checkpoints"):
        AH.head_paths(tmp_path / "empty")


# --------------------------------------------------------------------------------------------------
# route_deviation
# --------------------------------------------------------------------------------------------------
def test_route_deviation_on_identical_inputs_is_zero_and_declines_a_degenerate_correlation():
    v = torch.zeros(2, 5)
    d = AH.route_deviation(v, v)
    assert d["max_abs_dev"] == 0.0 and d["median_abs_dev"] == 0.0 and d["n"] == 10
    assert d["pearson"] is None, "a constant vector has no correlation to report"


def test_route_deviation_reports_a_real_spread_and_refuses_a_shape_mismatch():
    a = torch.tensor([[0.0, 1.0, 2.0, 3.0]])
    b = torch.tensor([[0.0, 1.0, 2.0, 9.0]])
    d = AH.route_deviation(a, b)
    assert d["max_abs_dev"] == pytest.approx(6.0)
    assert d["pearson"] is not None and d["pearson"] > 0.5
    with pytest.raises(AH.AtlasHeadsError, match="shape mismatch"):
        AH.route_deviation(a, torch.zeros(1, 3))


# --------------------------------------------------------------------------------------------------
# clone-touching -- skipped cleanly when the volatile temp dir has gone
# --------------------------------------------------------------------------------------------------
@needs_clone
def test_the_real_clone_ships_FIVE_heads_all_of_them_geometry_locked_to_86bp():
    base = AH.load_base_state_dict()
    assert len(base) == AH.BASE_STATE_DICT_TENSORS
    assert sum(v.numel() for v in base.values()) == AH.BASE_STATE_DICT_ELEMENTS

    heads = AH.load_heads(EXPR_DIR, base_state_dict=base)
    assert len(heads) == 5, "the arm's multiplicity divisor is the number of heads"
    for h in heads:
        assert tuple(h.state_dict["fc1.weight"].shape) == AH.FC1_SHAPE
        assert sum(v.numel() for v in h.state_dict.values()) == AH.EXPECTED_HEAD_ELEMENTS
        ident = AH.base_identity(h.state_dict, base)
        assert ident["n_identical"] == ident["n_checked"] == AH.BASE_STATE_DICT_TENSORS
        assert ident["max_abs_dev"] == 0.0

    # 'trspred' sorts BEFORE 'tslpred' ('r' < 's'), so four transcription heads come first
    assert sorted(h.endpoint for h in heads) == [AH.TRANSCRIPTION] * 4 + [AH.TRANSLATION]
    assert [h.name for h in heads if h.log_transform] == ["promoteratlas-tslpred-kosuri2013"]
