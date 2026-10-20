"""The GLM's `temperature` / `top_k` knobs were SHIPPED DEAD. These pin the fix.

MEASURED 2026-10-20: `generate(temperature=T)` had no effect at all -- T=0.7 vs T=1.5 at the
same seed produced 6/6 BYTE-IDENTICAL sequences. The mechanism is an ordering inversion, and
the upstream docstring asserts the opposite of what the code does:

  * the remote `_BPLogitsProcessor` marginalises k-mer logits to per-base probabilities, samples,
    and returns a ONE-HOT (`-inf` everywhere, `0.0` at the chosen token);
  * HF appends `TemperatureLogitsWarper` AFTER the merged user/custom processors
    (`transformers/generation/utils.py` -- merge at 1257, temperature at 1276, read on 5.14.1);
  * dividing a one-hot by any temperature is identity.

Sampling itself was never broken (sequences within one call differ) -- only the knobs.
"""
from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")
from transformers import TemperatureLogitsWarper, TopKLogitsWarper  # noqa: E402

from dna_decode.glm.generate import HFGenerator  # noqa: E402


def test_temperature_on_a_one_hot_is_IDENTITY_which_is_why_the_knob_was_dead():
    """The keystone. No model needed -- this is the whole bug in four lines.

    If this ever fails, the premise of the fix is gone and the fix should be re-derived.
    """
    one_hot = torch.full((1, 5), float("-inf"))
    one_hot[0, 2] = 0.0
    warped = TemperatureLogitsWarper(0.5)(torch.zeros((1, 1), dtype=torch.long), one_hot.clone())
    assert torch.equal(warped, one_hot), "temperature changed a one-hot; the bug analysis is wrong"


def test_the_same_warper_DOES_bite_on_real_logits_so_the_test_above_is_not_vacuous():
    """A warper that changed nothing anywhere would make the keystone meaningless."""
    lg = torch.tensor([[1.0, 2.0, 3.0]])
    warped = TemperatureLogitsWarper(0.5)(torch.zeros((1, 1), dtype=torch.long), lg.clone())
    assert not torch.equal(warped, lg)
    assert warped.tolist() == [[2.0, 4.0, 6.0]]


def _gen(**kw) -> HFGenerator:
    return HFGenerator("stub/model", **kw)


def test_a_non_default_temperature_installs_a_warper():
    procs = _gen(temperature=0.7)._sampling_processors()
    assert any(isinstance(p, TemperatureLogitsWarper) for p in procs)


def test_the_default_temperature_installs_NOTHING():
    """T=1.0 is identity, so adding a warper would be noise in the processor chain."""
    assert _gen(temperature=1.0)._sampling_processors() == []


def test_top_k_installs_a_warper_and_zero_means_disabled():
    # 40 is above the measured reachability floor; see knob_contract.TOP_K_KMER_FLOOR.
    assert any(isinstance(p, TopKLogitsWarper) for p in _gen(top_k=40)._sampling_processors())
    assert _gen(top_k=0)._sampling_processors() == []


def test_the_warpers_reach_generate_as_logits_processor_and_the_INERT_kwarg_is_GONE():
    """The regression pin, and the only one that can catch a silent revert.

    `temperature=` is the kwarg that was dead: HF appends it after the one-hot collapse. The
    fix works ONLY because the remote generate builds `existing + [bp_proc]`, so anything in
    `logits_processor` runs BEFORE the collapse. Passing `temperature=` again would be inert,
    and nothing else in the suite would notice.
    """
    seen: dict = {}

    class _Tok:
        bos_token_id = 1
        pad_token_id = 0
        eos_token_id = 2

        def __call__(self, *_a, **_k):
            class R:
                input_ids = torch.ones((1, 2), dtype=torch.long)
            return R()

        def decode(self, *_a, **_k):
            return "ACGTACGT"

    class _Model:
        def generate(self, _ids, **kwargs):
            seen.update(kwargs)
            return torch.ones((1, 3), dtype=torch.long)

    g = _gen(temperature=0.7, top_k=40)   # >= the measured top_k floor
    g.load = lambda: None            # type: ignore[method-assign]
    g._tok, g._model, g._dev = _Tok(), _Model(), "cpu"
    g.generate(1, 12)

    assert "logits_processor" in seen, "warpers never reached generate()"
    kinds = [type(p).__name__ for p in seen["logits_processor"]]
    assert "TemperatureLogitsWarper" in kinds and "TopKLogitsWarper" in kinds
    assert "temperature" not in seen, (
        "`temperature=` is back on the generate() call. HF applies it AFTER the bp processor's "
        "one-hot collapse, so it is inert -- the warper must go through logits_processor instead.")
    assert "top_k" not in seen, "`top_k=` is back on the generate() call; same inert path"
