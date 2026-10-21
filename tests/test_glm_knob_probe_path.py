"""The knob verification must exercise the path the sweep actually uses.

WHY (2026-10-21): the probe called the UNCONDITIONAL `generate` while every sweep point runs
`generate_conditional`. Both paths do install the sampling warpers (`generate.py:232` and
`:277`), so the practical risk was low -- but the entire subject of this episode is a knob that
LOOKED wired and was inert, and proving liveness on a path no grid point takes is a proxy for
the thing under test. That is the shape of evidence that failed the first time.
"""
from __future__ import annotations

import ast
import importlib.util
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent


def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "_glm_falsifier_probe", str(ROOT / "scripts" / "glm_generator_falsifier.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _StubGen:
    """Records which generation path was taken, and at what temperature."""

    def __init__(self):
        self.calls = []
        self.temperature = 1.0

    def generate_conditional(self, prompts, length, **kw):
        self.calls.append(("generate_conditional", tuple(prompts), length, self.temperature))
        return [f"COND{self.temperature}"]

    def generate(self, n, length, **kw):
        self.calls.append(("generate", n, length, self.temperature))
        return [f"UNCOND{self.temperature}"]


# ------------------------------------------------------------------ the fix, behaviourally
def test_the_probe_uses_the_CONDITIONAL_path_when_prompts_exist():
    mod = _load_runner()
    g = _StubGen()
    mod.make_knob_probe(g, ["ACGTACGT"], 150)(0.7)
    assert [c[0] for c in g.calls] == ["generate_conditional"]


def test_the_probe_passes_exactly_ONE_prompt():
    """A probe is a 2-call diff, not a run: generating 120 sequences twice would cost 44 min."""
    mod = _load_runner()
    g = _StubGen()
    mod.make_knob_probe(g, ["P1", "P2", "P3"], 150)(0.7)
    assert g.calls[0][1] == ("P1",)


def test_the_probe_sets_the_temperature_it_was_asked_for():
    """Non-vacuity: without this, a probe that ignored its argument would still pass above."""
    mod = _load_runner()
    g = _StubGen()
    probe = mod.make_knob_probe(g, ["ACGT"], 150)
    probe(0.7)
    probe(1.5)
    assert [c[3] for c in g.calls] == [0.7, 1.5]


def test_two_temperatures_produce_DIFFERENT_output_through_the_probe():
    """The property `verify_knob_live` actually tests, end to end through this helper."""
    mod = _load_runner()
    g = _StubGen()
    probe = mod.make_knob_probe(g, ["ACGT"], 150)
    assert probe(0.7) != probe(1.5)


def test_the_UNCONDITIONAL_fallback_still_works_when_there_are_no_prompts():
    """Unreachable today (promoter+conditional is the only valid pairing) but kept correct."""
    mod = _load_runner()
    for empty in (None, []):
        g = _StubGen()
        mod.make_knob_probe(g, empty, 150)(0.7)
        assert [c[0] for c in g.calls] == ["generate"]
        assert g.calls[0][1] == 1, "the probe generates one sequence, not n"


# ------------------------------------------- ordering + restoration, pinned on the source
# These two are STRUCTURAL and say so. The prompt hoist and the try/finally live inside `main`,
# which cannot be driven here without a real model, a real genome and ~22 minutes.
def _main_src() -> str:
    tree = ast.parse((ROOT / "scripts" / "glm_generator_falsifier.py").read_text(encoding="utf-8"))
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
    return ast.unparse(fn)


def test_prompts_are_built_BEFORE_the_knob_verification():
    """If the hoist is reverted, the probe has no prompts and silently falls back."""
    src = _main_src()
    i_prompts = src.find("window_context(")
    # The CALL, not the import. `verify_knob_live` is imported at the top of the same block, so
    # a bare-name search finds the import and reports the ordering backwards -- which is what
    # the first version of this test did.
    i_verify = src.find("verify_knob_live(")
    assert i_prompts != -1 and i_verify != -1
    assert i_prompts < i_verify, "prompt construction must precede the knob verification"


def test_the_genome_is_loaded_ONCE():
    """The hoist moved a load_fasta call; a careless hoist loads the 4.6 MB genome twice."""
    full = (ROOT / "scripts" / "glm_generator_falsifier.py").read_text(encoding="utf-8")
    assert full.count("_lf(REFSEQ") == 1, "genome loaded more than once after the hoist"


def test_the_temperature_is_restored_in_a_FINALLY():
    """An exception in the probe must not leave the generator at the probe's temperature."""
    tree = ast.parse((ROOT / "scripts" / "glm_generator_falsifier.py").read_text(encoding="utf-8"))
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
    tries = [t for t in ast.walk(fn)
             if isinstance(t, ast.Try) and "verify_knob_live" in ast.unparse(t)]
    assert tries, "verify_knob_live must be wrapped in a try"
    assert any("temperature" in ast.unparse(h) for t in tries for h in t.finalbody), \
        "the finally must restore gen.temperature"


def test_the_artifact_records_WHICH_path_was_verified():
    """A verification whose path is unrecorded cannot be audited later."""
    full = (ROOT / "scripts" / "glm_generator_falsifier.py").read_text(encoding="utf-8")
    assert "verified_on" in full
