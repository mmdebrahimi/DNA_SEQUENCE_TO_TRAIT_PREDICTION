"""Guards for the generator seam — chiefly the Pascal fp16 trap, encoded in code rather than prose.

`HFGenerator` is constructed but never loaded here: the whole point of the seam is that the harness is
testable on a host with no torch and no weights.
"""
from __future__ import annotations

import pytest

from dna_decode.glm.baselines import MarkovGenerator, UniformGenerator
from dna_decode.glm.generate import (
    DEFAULT_MODEL,
    FAST_FP16_MIN_CC,
    Generator,
    HFGenerator,
    device_report,
    recommended_dtype,
)


def test_the_baselines_satisfy_the_generator_protocol():
    """The seam's whole value: the falsifier cannot tell which generator it is scoring."""
    assert isinstance(UniformGenerator(), Generator)
    assert isinstance(MarkovGenerator(), Generator)
    assert isinstance(HFGenerator(), Generator)


def test_fp16_is_only_recommended_above_compute_capability_7():
    """THE trap. Consumer Pascal (CC 6.1) runs fp16 at ~1/64 of fp32, so 'half precision to save VRAM'
    makes a 1070 unusable rather than faster; Maxwell (5.0) has no fp16 units at all."""
    assert FAST_FP16_MIN_CC == (7, 0)
    assert recommended_dtype("cpu") == "float32", "fp16 on CPU is emulated and slower"


def test_recommended_dtype_matches_this_hosts_actual_capability():
    """Derived from the live device, not asserted — and fp32 is the answer on anything below 7.0."""
    rep = device_report()
    if not rep["cuda"]:
        assert rep["recommended_dtype"] == "float32"
        return
    major = int(rep["compute_capability"].split(".")[0])
    minor = int(rep["compute_capability"].split(".")[1])
    expected = "float16" if (major, minor) >= FAST_FP16_MIN_CC else "float32"
    assert rep["recommended_dtype"] == expected
    assert rep["fp16_fast"] is ((major, minor) >= FAST_FP16_MIN_CC)


def test_a_slow_fp16_device_gets_an_explicit_NOTE_not_silent_fp32():
    """Silently downgrading would hide the reason; an operator needs to know why it is slow."""
    rep = device_report()
    if rep["cuda"] and not rep["fp16_fast"]:
        joined = " ".join(rep["notes"])
        assert "1/64" in joined or "NOT fast" in joined
        assert "fp32" in joined.lower()


def test_device_report_is_total_and_never_raises():
    rep = device_report()
    for key in ("torch", "cuda", "recommended_dtype", "fp16_fast", "notes"):
        assert key in rep
    assert rep["recommended_dtype"] in ("float32", "float16")
    assert isinstance(rep["notes"], list)


def test_hf_generator_constructs_without_torch_weights_or_network():
    """Construction must be free; only `load()`/`generate()` may touch the heavy path."""
    g = HFGenerator("some/model", device="cpu")
    assert g.name == "hf:model"
    assert g.dtype == "float32"
    assert g._model is None and g._tok is None


def test_the_default_model_is_the_verified_mit_prokaryote_checkpoint():
    """Pinned because the v1 id does not exist (401) and only the v2 repo does -- verified against the
    HF API, not taken from a survey."""
    assert DEFAULT_MODEL == "GenerTeam/GENERator-v2-prokaryote-1.2b-base"
    assert "v2" in DEFAULT_MODEL and "prokaryote" in DEFAULT_MODEL


def test_an_explicit_dtype_overrides_the_recommendation():
    g = HFGenerator("m", device="cpu", dtype="float16")
    assert g.dtype == "float16", "an operator may override; the default must not be silently forced"


def test_the_vram_guard_estimates_from_the_config_and_scales_with_dtype():
    """Guard-before-the-heavy-load, and it was earned: the 1.2B checkpoint OOMed on the local 4 GiB card
    with "tried to allocate 44.00 MiB ... 0 bytes free" AFTER paying for the full 4.4 GB download and
    load. device_report() already knew the VRAM, so the crash was avoidable.
    """
    class _Cfg:
        hidden_size = 2048
        num_hidden_layers = 24
        vocab_size = 4096
        num_parameters = 1_200_000_000

    class _Auto:
        @staticmethod
        def from_pretrained(*a, **k):
            return _Cfg()

    g32 = HFGenerator("m", device="cpu", dtype="float32")
    g16 = HFGenerator("m", device="cpu", dtype="float16")
    n32 = g32._estimated_vram_gib(_Auto)
    n16 = g16._estimated_vram_gib(_Auto)
    assert n32 is not None and abs(n32 - 1.2e9 * 4 / 2**30) < 0.01
    assert abs(n32 / n16 - 2.0) < 1e-6, "fp32 must estimate exactly twice fp16"


def test_the_vram_guard_falls_back_rather_than_refusing_when_it_CANNOT_estimate():
    """An un-inspectable config must not block a load — refusing on a number we could not compute would
    be worse than trying. None means 'do not fire the guard'."""
    class _Bad:
        @staticmethod
        def from_pretrained(*a, **k):
            raise RuntimeError("no config")

    assert HFGenerator("m", device="cpu")._estimated_vram_gib(_Bad) is None

    class _Sparse:
        @staticmethod
        def from_pretrained(*a, **k):
            return type("C", (), {})()          # no num_parameters, no dims

    assert HFGenerator("m", device="cpu")._estimated_vram_gib(_Sparse) is None


def test_an_explicit_device_is_never_overridden_by_the_guard():
    """The guard only runs under device='auto'. An operator asking for cuda gets cuda (and its OOM)."""
    g = HFGenerator("m", device="cpu")
    assert g.device == "cpu" and g.fallback_reason is None


def test_a_fresh_generator_has_no_fallback_reason_so_the_field_cannot_read_as_stale():
    assert HFGenerator("m").fallback_reason is None


def test_generate_before_load_raises_rather_than_returning_empty():
    """An empty list would read as 'the generator produced nothing', not 'torch is missing'."""
    g = HFGenerator("definitely/not-a-real-model-xyz", device="cpu")
    with pytest.raises(Exception):
        g.generate(1, 60)
