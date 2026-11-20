"""Pins the three defects that cost three Kaggle pushes on 2026-10-31.

Each one was cheap to catch locally and expensive to catch remotely, because Kaggle returns ZERO
output files for a kernel whose status is ERROR -- so a failed push yields no log, no traceback and
an EMPTY failureMessage, which is indistinguishable from an infrastructure problem. The guards are
in the builder; these tests keep them there.

  1. A `.py`-ONLY PAYLOAD. The split script verifies the rebuilt corpus fingerprint against
     `wiki/glm_generator_falsifier_REAL_2026-10-07.json` and REFUSES (rc=3) without it. That
     refusal is correct -- a fingerprint it cannot compare is not one it has checked -- so the
     payload must carry the DATA input, not just the import closure.
  2. A BUILDER THAT VALIDATED ITSELF AND NOT ITS ARTIFACT. The kernel body is a raw triple-quoted
     string, so a literal newline inside a `print("...")` is legal in the builder and an
     unterminated string literal in the emitted file. Parse the ARTIFACT.
  3. A CHILD PROCESS WITH NO PATH. `sys.path.insert` in the kernel process does not reach the
     subprocess that runs the split script, which died on `ModuleNotFoundError: dna_decode`. The
     step-1 import check could not catch it, because the PARENT imported everything fine.

Pure static/structural checks: no torch, no network, no GPU, no Kaggle.
"""
from __future__ import annotations

import ast
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "kaggle" / "build_glm_finetune_kernel.py"
BASELINE = "wiki/glm_generator_falsifier_REAL_2026-10-07.json"


def _builder():
    spec = importlib.util.spec_from_file_location("_bgf", BUILDER)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_payload_ships_the_baseline_artifact_not_only_python():
    """Defect 1. A .py-only closure makes the split script refuse with rc=3."""
    b = _builder()
    assert BASELINE in b.FILES, (
        "the payload must ship the baseline artifact; without it the split script cannot check the "
        "corpus fingerprint and correctly refuses")
    assert (ROOT / BASELINE).is_file(), "the baseline named in FILES does not exist on disk"


def test_every_payload_file_exists():
    b = _builder()
    missing = [f for f in b.FILES if not (ROOT / f).is_file()]
    assert not missing, "payload names files that do not exist: %s" % missing


def test_the_emitted_kernel_parses():
    """Defect 2, the outcome: whatever the builder would write must be valid Python."""
    b = _builder()
    body = b.build_body()
    ast.parse(body)


def test_the_diagnostic_variant_also_parses():
    b = _builder()
    body = b.build_body(diagnostic=True)
    ast.parse(body)


def test_builder_refuses_to_write_a_kernel_that_does_not_parse(tmp_path, monkeypatch):
    """Defect 2, the guard. Proven by breaking the template, not by reading the code."""
    b = _builder()
    good = b.KERNEL
    b.KERNEL = good.replace("def sh(cmd):", "def sh(cmd:", 1)
    assert b.KERNEL != good, "could not inject a syntax error; the anchor moved"
    try:
        with pytest.raises(SystemExit) as e:
            b.main([])
        assert "does not parse" in str(e.value), str(e.value)
    finally:
        b.KERNEL = good


def test_the_split_subprocess_gets_pythonpath():
    """Defect 3. The child must be given the path; the parent's sys.path does not carry over."""
    b = _builder()
    body = b.KERNEL
    i = body.find("glm_finetune_window_split.py\"]")
    assert i > 0, "the split-script subprocess invocation moved; re-anchor this test"
    window = body[max(0, i - 1200):i + 400]
    assert "PYTHONPATH" in window, (
        "the subprocess that runs the split script must receive PYTHONPATH, or it dies on "
        "ModuleNotFoundError: dna_decode while the parent imports fine")
    assert "env=env" in window, "the env must actually be passed to subprocess.run"


def test_both_kernels_always_report_their_own_failure():
    """Kaggle exposes no output for an ERROR kernel, so a silent death is unreadable."""
    b = _builder()
    body = b.KERNEL
    assert "glm_finetune_failure.json" in body, (
        "the real run must write a failure report; without it a crash is invisible")
    assert "glm_finetune_diagnostic.json" in body
    assert body.count("raise SystemExit(0)") >= 2, (
        "both the diagnostic and the real run must exit 0 so their artifacts are retrievable")


def test_the_failure_report_is_never_the_verdict_file():
    """A dead run must not be mistakable for a measured one."""
    b = _builder()
    body = b.KERNEL
    i = body.find("glm-generator-finetune-FAILURE-v1")
    assert i > 0, "the failure record moved"
    block = body[i:i + 1200]
    assert "glm_generator_finetune_verdict.json" not in block, (
        "the failure path must not write the verdict file")
    assert '"verdict": None' in block or "'verdict': None" in block


def test_the_frozen_bar_matches_the_committed_acceptance_bar():
    """The kernel's thresholds must not drift from the pre-registration."""
    import json
    b = _builder()
    bar = json.loads((ROOT / "wiki" / "glm_generator_finetune_acceptance_bar.json")
                     .read_text(encoding="utf-8"))
    assert bar["primary_bar"]["threshold"] == 0.6679
    assert "BAR_BOTH = 0.6679" in b.KERNEL
    assert "BAR_KMER = 0.5990" in b.KERNEL
    m = bar["method_frozen"]
    assert '"r": 16' in b.KERNEL or "'r': 16" in b.KERNEL or '"r": 16' in str(m)
    for k in ("r", "lora_alpha", "lora_dropout"):
        assert str(m[k]) in b.KERNEL, "frozen %s=%s is not in the kernel" % (k, m[k])


def test_exactly_at_the_bar_is_not_a_pass():
    """The bar is STRICTLY below. A score equal to the null's own value is not beating it."""
    b = _builder()
    ns = {"__name__": "_v"}
    exec(compile(b.build_body(), "k.py", "exec"), ns)
    v, _why = ns["verdict"]({"both": 0.6679}, {"both": 0.6679}, True)
    assert v == "C_NO_IMPROVEMENT", v
    v, _why = ns["verdict"]({"both": 0.77}, {"both": 0.6679}, True)
    assert v == "B_IMPROVES_BUT_BAR_UNMET", v
    v, _why = ns["verdict"]({"both": 0.77}, {"both": 0.6678}, True)
    assert v == "A_BEATS_NULL", v
