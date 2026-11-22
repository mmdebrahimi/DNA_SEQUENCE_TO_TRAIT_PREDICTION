"""Pin that `diffusion.sample` is device-correct, and that fixing it did not move CPU output.

WHY THIS EXISTS. `sample` used to build `torch.Generator(device="cpu")` and then call
`torch.rand(x.shape, generator=g)` with no device, so with a CUDA model it raised a device mismatch:
the sampler was CPU-only by construction. The failure surfaced at SAMPLING -- after a full training
spend -- which made it the most expensive failure the GLM arm had available.

TWO KINDS OF CHECK, because this host has no usable CUDA and a test that can only skip is not a
guard. The behavioural CPU pins run everywhere. The structural pins read the source and are what
actually fail if someone reintroduces a hardcoded CPU generator on a machine that cannot exercise
the CUDA path. The CUDA test runs where there is a GPU and skips honestly where there is not --
it is the strongest check and also the one least often executed, which is exactly why the other
two tiers are here.

CPU EQUIVALENCE IS LOAD-BEARING, not cosmetic: `wiki/glm_diffusion_samples_2026-11-21.txt` was
produced by the pre-fix code. A changed CPU RNG stream would silently invalidate it.
"""
from __future__ import annotations

import hashlib
import inspect
from pathlib import Path

import pytest

from dna_decode.glm import diffusion
from dna_decode.glm.diffusion import MASK_ID, DiffusionConfig, build_model, sample

torch = pytest.importorskip("torch")

#: Digests captured from the PRE-FIX implementation. Do not regenerate them to make a failure go
#: away -- a mismatch means CPU sampling moved, which is a real finding about committed artifacts.
CPU_DIGESTS = {
    "seed0_steps8":  "3792655cd9fdeb41251aa830cf91bbcf",
    "seed0_steps64": "e0465e8af088b4a8be1dbe8b817e17ac",
    "seed1_steps8":  "c891d213230d4c2281336487e77ceba0",
    "seed1_steps64": "ecc1c5e0df201150af957c76d368a4e9",
    "seed7_steps8":  "2405c19225b46deaa30999f94df4ab8b",
    "seed7_steps64": "3a5e9a62ed503241ab72443658a6b688",
}


def _tiny():
    cfg = DiffusionConfig(length=150, d_model=64, n_layers=2, n_heads=4, d_ff=128, seed=7)
    model = build_model(cfg)
    model.eval()
    return cfg, model


def _digest(seqs) -> str:
    return hashlib.sha256("\n".join(seqs).encode()).hexdigest()[:32]


@pytest.mark.parametrize("seed,steps", [(0, 8), (0, 64), (1, 8), (1, 64), (7, 8), (7, 64)])
def test_cpu_sampling_is_byte_identical_to_the_prefix_implementation(seed, steps):
    cfg, model = _tiny()
    got = _digest(sample(model, 6, cfg, steps=steps, seed=seed))
    want = CPU_DIGESTS["seed%d_steps%d" % (seed, steps)]
    assert got == want, (
        "CPU sampling moved (%s != %s). wiki/glm_diffusion_samples_2026-11-21.txt was produced by "
        "this code path, so a changed stream invalidates it. Do not update the digest to pass."
        % (got, want))


def test_sampling_is_deterministic_per_seed_and_differs_across_seeds():
    cfg, model = _tiny()
    a = sample(model, 4, cfg, steps=16, seed=3)
    assert a == sample(model, 4, cfg, steps=16, seed=3), "same seed must reproduce"
    assert a != sample(model, 4, cfg, steps=16, seed=4), (
        "different seeds produced identical output -- the seed is not reaching the sampler, which "
        "would make a seed-spread tolerance vacuous")


def test_output_is_complete_clean_sequence():
    cfg, model = _tiny()
    for s in sample(model, 8, cfg, steps=12, seed=0):
        assert len(s) == cfg.length
        assert set(s) <= set("ACGT"), "a leftover MASK would decode as a base and score as real"


def test_mask_check_runs_before_decode():
    """Order matters: decoding first would turn a leftover absorbing state into a base silently."""
    src = inspect.getsource(sample)
    assert "assert_no_mask_remains" in src
    assert src.index("assert_no_mask_remains") < src.index("return decode("), (
        "assert_no_mask_remains must precede decode, or an incomplete reverse process is reported "
        "as clean sequence")


# ---------------------------------------------------------------- structural pins (no GPU needed)

def test_sampler_does_not_hardcode_a_cpu_generator():
    src = inspect.getsource(sample)
    assert 'Generator(device="cpu")' not in src, (
        "a hardcoded CPU generator is the original defect: it raises against a CUDA model, and the "
        "raise lands at sampling time, after the whole training spend")
    assert "next(model.parameters()).device" in src, (
        "the sampler must take its device from the MODEL, not assume one")


def test_every_random_draw_in_the_sampler_is_device_placed():
    """`torch.rand` and `torch.full` without `device=` default to CPU and break a CUDA run."""
    src = inspect.getsource(sample)
    for call in ("torch.rand(", "torch.full("):
        for line in src.splitlines():
            if call in line:
                assert "device=" in line, (
                    "%r has no device= in: %s" % (call, line.strip()))


def test_the_shipped_module_is_the_one_the_standalone_mirrors():
    """The Server1 job embeds a copy of this formulation; drift makes the two runs incomparable.

    Pins the handful of constants the two implementations MUST agree on. It deliberately does not
    diff the files -- the standalone adds conditioning and is expected to differ in body.
    """
    standalone = (Path(__file__).resolve().parents[1]
                  / "scripts" / "server1" / "glm_human_diffusion_server1.py")
    if not standalone.is_file():          # the job may be retired; that is not a failure here
        pytest.skip("standalone Server1 job not present")
    text = standalone.read_text(encoding="utf-8")
    assert "MASK_ID = %d" % MASK_ID in text, "absorbing-state id must match the shipped module"
    assert "VOCAB = %d" % diffusion.VOCAB in text, "vocab must match"
    assert 'BASES = "%s"' % diffusion.BASES in text, "base order must match or encodings diverge"
    assert "device=dev, generator=g" in text, "the standalone must carry the device fix too"


# ---------------------------------------------------------------------------------- the real check

@pytest.mark.skipif(not torch.cuda.is_available(), reason="no CUDA on this host")
def test_sampling_runs_on_cuda_and_returns_clean_sequence():
    cfg, model = _tiny()
    model = model.to("cuda")
    seqs = sample(model, 4, cfg, steps=12, seed=0)
    assert len(seqs) == 4
    for s in seqs:
        assert len(s) == cfg.length
        assert set(s) <= set("ACGT")
    assert all(MASK_ID != 0 for _ in seqs)      # sanity: MASK is not a base id


@pytest.mark.skipif(not torch.cuda.is_available(), reason="no CUDA on this host")
def test_cuda_and_cpu_streams_differ_and_that_is_documented():
    """Recorded as a property, not a bug: a seed reproduces on the SAME device only."""
    cfg, model = _tiny()
    on_cpu = sample(model, 4, cfg, steps=12, seed=0)
    on_gpu = sample(model.to("cuda"), 4, cfg, steps=12, seed=0)
    assert len(on_cpu) == len(on_gpu)
    if on_cpu == on_gpu:
        pytest.skip("streams happened to agree; no claim either way from one tiny model")
