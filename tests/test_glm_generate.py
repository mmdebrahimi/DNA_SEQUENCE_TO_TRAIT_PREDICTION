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


def test_generate_before_load_raises_rather_than_returning_empty():
    """An empty list would read as 'the generator produced nothing', not 'torch is missing'."""
    g = HFGenerator("definitely/not-a-real-model-xyz", device="cpu")
    with pytest.raises(Exception):
        g.generate(1, 60)
