"""The diffusion generator must actually learn, must never emit a leftover mask, and must refuse junk.

The load-bearing test here is `test_the_model_can_memorise_a_single_sequence`. A training loop can run,
report a falling loss and produce well-formed output while learning nothing -- a mis-wired loss, a
detached graph, an optimizer stepping on the wrong parameters all look exactly like slow convergence.
Overfitting a single sequence is the cheap decisive check: if the model cannot reproduce one sequence it
has seen a thousand times, nothing else in the file means anything.

The second is `test_sampling_leaves_no_mask`. The absorbing state is token id 4 and `decode` would turn
a leftover into whatever base sits at that index -- except there is no such base, so it would raise, and
if the vocabulary ever grew it would silently become a real base scored as real sequence. The sampler
asserts completion rather than trusting its own schedule.
"""

from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from dna_decode.glm.diffusion import (  # noqa: E402
    BASES,
    MASK_ID,
    DiffusionConfig,
    assert_no_mask_remains,
    build_model,
    corrupt,
    decode,
    encode,
    loss_fn,
    sample,
)

SEQ = "ACGT" * 8                      # 32 bp, long enough to be non-trivial and fast to overfit


def tiny_cfg(**kw) -> DiffusionConfig:
    base = {"length": len(SEQ), "d_model": 32, "n_layers": 2, "n_heads": 2, "d_ff": 64,
            "dropout": 0.0, "seed": 0}
    base.update(kw)
    return DiffusionConfig(**base)


# ------------------------------------------------------------------ encoding

def test_encode_decode_roundtrip():
    x = encode([SEQ, SEQ[::-1]], len(SEQ))
    assert x.shape == (2, len(SEQ))
    assert decode(x) == [SEQ, SEQ[::-1]]


def test_encode_refuses_a_non_acgt_character():
    """An N coerced to a base would teach the model that ambiguity codes are real sequence."""
    with pytest.raises(ValueError, match="non-ACGT"):
        encode(["ACGN" * 8], len(SEQ))


def test_encode_refuses_a_wrong_length():
    """A short window silently padded would teach the model that windows end early."""
    with pytest.raises(ValueError, match="length"):
        encode(["ACGT"], len(SEQ))


def test_encode_accepts_lowercase():
    assert decode(encode([SEQ.lower()], len(SEQ))) == [SEQ]


# ------------------------------------------------------------------ the forward process

def test_corruption_rate_tracks_t():
    """t is the masking probability; if it did not track, the 1/t weighting would be meaningless."""
    x0 = encode([SEQ] * 400, len(SEQ))
    for t_val in (0.1, 0.5, 0.9):
        t = torch.full((x0.shape[0],), t_val)
        xt, mask = corrupt(x0, t)
        rate = float(mask.float().mean())
        assert abs(rate - t_val) < 0.05, (t_val, rate)
        assert bool((xt[mask] == MASK_ID).all())


def test_unmasked_positions_are_left_alone():
    x0 = encode([SEQ] * 50, len(SEQ))
    xt, mask = corrupt(x0, torch.full((50,), 0.5))
    assert bool((xt[~mask] == x0[~mask]).all())


def test_t_equal_to_one_masks_everything():
    x0 = encode([SEQ] * 20, len(SEQ))
    xt, mask = corrupt(x0, torch.ones(20))
    assert bool(mask.all()) and bool((xt == MASK_ID).all())


# ------------------------------------------------------------------ the loss

def test_the_loss_is_finite_and_positive():
    cfg = tiny_cfg()
    m = build_model(cfg)
    loss = loss_fn(m, encode([SEQ] * 8, cfg.length), generator=torch.Generator().manual_seed(0))
    assert torch.isfinite(loss) and float(loss) > 0


def test_the_loss_weights_by_inverse_t_not_uniformly():
    """The 1/t weight is what makes this a likelihood bound rather than an arbitrary denoising loss.

    Pinned on the source because the alternative -- a plain mean over masked positions -- trains a model
    that fills one gap well and cannot generate from scratch, which would not show up as an error.
    """
    from pathlib import Path

    src = (Path(__file__).resolve().parents[1] / "dna_decode" / "glm"
           / "diffusion.py").read_text(encoding="utf-8")
    assert "1.0 / t" in src
    assert "likelihood" in src


# ------------------------------------------------------------------ THE learning check

def test_the_model_can_memorise_a_single_sequence():
    """If it cannot reproduce one sequence it has seen a thousand times, the loop is broken.

    A falling loss alone proves nothing: a detached graph, a mis-wired target or an optimizer stepping
    on the wrong parameters all look like slow convergence. This asserts the model recovers the actual
    bases from a FULLY masked input, which is the hardest case and exactly what sampling starts from.
    """
    cfg = tiny_cfg()
    m = build_model(cfg)
    opt = torch.optim.AdamW(m.parameters(), lr=3e-3)
    x = encode([SEQ] * 16, cfg.length)
    g = torch.Generator().manual_seed(0)
    m.train()
    for _ in range(400):
        loss = loss_fn(m, x, generator=g)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
    m.eval()
    with torch.no_grad():
        fully_masked = torch.full((1, cfg.length), MASK_ID, dtype=torch.long)
        pred = m(fully_masked).argmax(-1)
    recovered = decode(pred)[0]
    matches = sum(a == b for a, b in zip(recovered, SEQ))
    assert matches >= int(0.9 * cfg.length), (matches, cfg.length, recovered)


def test_an_untrained_model_does_NOT_recover_the_sequence():
    """Control for the test above: the recovery must come from training, not from the architecture."""
    cfg = tiny_cfg(seed=3)
    m = build_model(cfg)
    m.eval()
    with torch.no_grad():
        pred = m(torch.full((1, cfg.length), MASK_ID, dtype=torch.long)).argmax(-1)
    matches = sum(a == b for a, b in zip(decode(pred)[0], SEQ))
    assert matches < int(0.9 * cfg.length), matches


# ------------------------------------------------------------------ sampling

def test_sampling_leaves_no_mask():
    """A leftover absorbing state would be decoded and scored as real sequence."""
    cfg = tiny_cfg()
    m = build_model(cfg)
    out = sample(m, 4, cfg, steps=16, seed=0)
    assert len(out) == 4
    for s in out:
        assert len(s) == cfg.length
        assert set(s) <= set(BASES)


def test_sampling_is_deterministic_under_a_seed():
    cfg = tiny_cfg()
    m = build_model(cfg)
    assert sample(m, 3, cfg, steps=8, seed=11) == sample(m, 3, cfg, steps=8, seed=11)


def test_different_seeds_give_different_samples():
    """Determinism must not be degeneracy: a model emitting one sequence regardless would pass above."""
    cfg = tiny_cfg()
    m = build_model(cfg)
    assert sample(m, 3, cfg, steps=8, seed=1) != sample(m, 3, cfg, steps=8, seed=2)


def test_the_completion_assert_actually_fires():
    """Non-vacuity for the guard: it must raise on a tensor that still holds the absorbing state."""
    x = torch.full((2, 8), MASK_ID, dtype=torch.long)
    with pytest.raises(RuntimeError, match="still MASK"):
        assert_no_mask_remains(x)
    assert_no_mask_remains(torch.zeros((2, 8), dtype=torch.long))


def test_one_sampling_step_still_completes():
    """The final step must force every remaining mask to resolve, even with no gradual schedule."""
    cfg = tiny_cfg()
    m = build_model(cfg)
    out = sample(m, 2, cfg, steps=1, seed=0)
    assert all(set(s) <= set(BASES) and len(s) == cfg.length for s in out)


# ------------------------------------------------------------------ config

def test_the_parameter_estimate_is_close_to_the_real_count():
    """The estimate is used to size a run against the corpus BEFORE anything is built."""
    cfg = tiny_cfg()
    real = sum(q.numel() for q in build_model(cfg).parameters())
    assert 0.9 <= cfg.n_params_estimate() / real <= 1.1


def test_the_config_records_the_vocabulary_and_mask_id():
    d = tiny_cfg().as_dict()
    assert d["vocab"] == len(BASES) + 1
    assert d["mask_id"] == MASK_ID


def test_positional_embeddings_exist_because_position_is_the_point():
    """Today's measurement says positional architecture is what the generator must learn; a model
    without a positional table could not represent where a motif sits."""
    m = build_model(tiny_cfg())
    assert hasattr(m, "pos")
    assert m.pos.num_embeddings == tiny_cfg().length
