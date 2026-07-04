"""Offline mechanism tests for the genome-JEPA prototype (CPU, fast, deterministic).

These pin the LOAD-BEARING JEPA properties: stop-gradient EMA target, no representation collapse, real
local-dependency capture (structured beats shuffled). They do NOT claim genome value (that's the
workhorse-gated falsifier) -- they verify the mechanism is correct.
"""
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import genome_jepa_prototype as gj  # noqa: E402


def test_shapes():
    x = gj.make_markov_dna(8, 32, seed=0)
    assert x.shape == (8, 32) and x.dtype == torch.long and int(x.max()) < gj.VOCAB
    m = gj.GenomeJEPA(d_model=16, block=8)
    enc = m.context_encoder(x)
    assert enc.shape == (8, 32, 16)
    loss, tvar = m(x, start=4)
    assert loss.ndim == 0 and tvar.ndim == 0


def test_target_encoder_is_frozen_and_stop_grad():
    m = gj.GenomeJEPA(d_model=16, block=8)
    assert all(not p.requires_grad for p in m.target_encoder.parameters())
    x = gj.make_markov_dna(8, 32, seed=1)
    loss, _ = m(x, start=0)
    loss.backward()
    # target-encoder params receive NO gradient (stop-gradient); context-encoder params DO.
    assert all(p.grad is None for p in m.target_encoder.parameters())
    assert any(p.grad is not None for p in m.context_encoder.parameters())


def test_ema_moves_target_toward_context():
    m = gj.GenomeJEPA(d_model=16, block=8, ema=0.9)
    # perturb the online encoder so context != target
    with torch.no_grad():
        for p in m.context_encoder.parameters():
            p.add_(torch.randn_like(p))
    before = [tp.clone() for tp in m.target_encoder.parameters()]
    m.update_target()
    moved = False
    for tp0, tp1, cp in zip(before, m.target_encoder.parameters(), m.context_encoder.parameters()):
        assert torch.all(torch.abs(tp1 - cp) <= torch.abs(tp0 - cp) + 1e-6)  # moved toward context
        if not torch.allclose(tp0, tp1):
            moved = True
    assert moved


def test_loss_decreases_on_structured_data():
    x = gj.make_markov_dna(64, 32, order=2, seed=0)
    _, hist = gj.train_prototype(x, steps=60, d_model=16, block=8, seed=0)
    start = sum(hist["loss"][:10]) / 10
    end = sum(hist["loss"][-10:]) / 10
    assert end < start, f"loss did not decrease: {start:.4f} -> {end:.4f}"


def test_no_representation_collapse():
    x = gj.make_markov_dna(64, 32, order=2, seed=0)
    _, hist = gj.train_prototype(x, steps=60, d_model=16, block=8, seed=0)
    final_var = sum(hist["tvar"][-10:]) / 10
    assert final_var > 1e-3, f"representation collapsed: target variance {final_var:.5f}"


def test_structured_beats_shuffled_negative_control():
    # Dependency-capture is a SCALE-dependent emergent property: it only manifests reliably at the
    # decisive config (measured empirically — it is noise at small scale). This test is the slow one
    # (~40s CPU) but it is the honest, deterministic check that JEPA learns real local dependency, not
    # a trivial constant. (The fast mechanism tests above cover collapse/EMA/loss-decrease.)
    struct = gj.make_markov_dna(256, 64, order=2, seed=0)
    shuf = gj.shuffle_positions(struct, seed=1)
    _, hs = gj.train_prototype(struct, steps=300, d_model=32, block=8, seed=0)
    _, hz = gj.train_prototype(shuf, steps=300, d_model=32, block=8, seed=0)
    s = sum(hs["loss"][-20:]) / 20
    z = sum(hz["loss"][-20:]) / 20
    assert z - s > 0.05, f"no real dependency captured: structured {s:.4f} vs shuffled {z:.4f} (gap {z-s:+.4f})"
