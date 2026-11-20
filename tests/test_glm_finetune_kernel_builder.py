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
import re
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

# ---------------------------------------------------------------------------------------------
# Step 2 + Step 3, 2026-11-20: the paired-delta variant.
#
# THE NO-WRITE INVARIANT below is the load-bearing one. The builder's file-writing entry point
# defaults to scripts/kaggle/glm_finetune_kernel.py -- the provenance record of the two
# committed arms -- so a test that builds with no flag and then compares would OVERWRITE the
# artifact it is checking. Note that it is NOT "never call main()": the pre-existing
# test_builder_refuses_to_write_a_kernel_that_does_not_parse calls main([]) deliberately and
# safely, because the AST guard raises before write_text. The invariant is on the FILE, not the
# call, which is why it is a sha256 snapshot.
# ---------------------------------------------------------------------------------------------

FROZEN_KERNEL = ROOT / "scripts" / "kaggle" / "glm_finetune_kernel.py"
PAIRED_KERNEL = ROOT / "scripts" / "kaggle" / "glm_paired_sweep_kernel.py"


def _sha(p: Path) -> str:
    import hashlib
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.fixture(autouse=True, scope="module")
def _frozen_kernel_is_not_rewritten():
    """No test in this module may rewrite the frozen two-arm artifact."""
    before = _sha(FROZEN_KERNEL)
    yield
    after = _sha(FROZEN_KERNEL)
    assert before == after, (
        "a test in this module rewrote %s. That file is the provenance record of the two "
        "committed arms (0.7787 / 0.8058); regenerating it is a deliberate act, never a side "
        "effect of a test." % FROZEN_KERNEL.name)


def test_the_seed_set_is_defined_once_and_contains_the_committed_arms_seed():
    b = _builder()
    seeds = b.SEED_SET
    assert len(seeds) == len(set(seeds)), "duplicate seeds: %r" % (seeds,)
    assert len(seeds) >= 3, "fewer than 3 seeds cannot bound a spread"
    assert 7 in seeds, (
        "7 is the falsifier CLI default and therefore the seed both committed arms ran at; "
        "without it the anchors against 0.7787 / 0.8058 are not checkable")
    # Defined ONCE: the template must not carry a second literal copy.
    assert "SEEDS = __SEEDS__" in b.KERNEL, "the template must take SEEDS by substitution"
    assert "SEED_SET = (" in BUILDER.read_text(encoding="utf-8")


def test_seed_is_appended_only_under_a_conditional():
    """Step 2. The frozen two-arm argv must be byte-for-byte what it was, so --seed may never
    be appended unconditionally."""
    body = _builder().build_body()
    assert "if seed is not None:" in body
    # Every occurrence of the flag must be the guarded append.
    assert body.count('"--seed"') == body.count('argv += ["--seed", str(seed)]')


def test_run_arm_keeps_seed_keyword_only():
    """A positional seed would silently shift an existing caller's argument."""
    body = _builder().build_body()
    assert "def run_arm(cli, label, model_path, device, outfile, *, seed=None):" in body


def test_build_body_refuses_an_unsubstituted_placeholder():
    """Proven by injecting one, not by reading the guard."""
    b = _builder()
    good = b.KERNEL
    b.KERNEL = good + "\nNEW = __UNSUBSTITUTED__\n"
    try:
        with pytest.raises(SystemExit) as e:
            b.build_body()
        assert "__UNSUBSTITUTED__" in str(e.value), str(e.value)
    finally:
        b.KERNEL = good


def test_diagnostic_and_paired_are_mutually_exclusive():
    """The diagnostic stops before the model download, so it can never produce a grid."""
    b = _builder()
    with pytest.raises(SystemExit) as e:
        b.build_body(diagnostic=True, paired=True)
    assert "mutually exclusive" in str(e.value)


def test_the_paired_variant_parses():
    import ast
    ast.parse(_builder().build_body(paired=True))


def test_the_paired_variant_saves_one_adapter_per_epoch():
    body = _builder().build_body(paired=True)
    assert 'ep_dir = ADAPTER + "/ep%d" % (ep + 1)' in body
    assert "peft_model.save_pretrained(ep_dir)" in body
    # And it must re-merge each of them in a loop over the epochs.
    assert "for ep in range(1, EPOCHS + 1):" in body


def test_the_paired_variant_rewrites_the_grid_after_every_cell():
    """A session cut mid-grid must leave a readable partial grid, not an average over an
    unknown denominator."""
    body = _builder().build_body(paired=True)
    sm = body.split("def score_model_at_every_seed(")[1].split("\ndef ")[0]
    assert "glm_paired_sweep_grid.json" in sm
    assert "for s in SEEDS:" in sm


def test_the_paired_variant_merges_through_step3_model_for_the_tokenizer_wiring():
    """A bare from_pretrained would skip setup_tokenizer, and the remote modeling code only
    PRINTS a warning when its own re-fetch fails -- generate() then raises far from the cause."""
    body = _builder().build_body(paired=True)
    ps = body.split("def paired_sweep(")[1].split("\ndef ")[0]
    assert "step3_model()" in ps
    assert "AutoModelForCausalLM.from_pretrained" not in ps


def test_the_paired_path_never_calls_the_frozen_verdict():
    body = _builder().build_body(paired=True)
    ps = body.split("def paired_sweep(")[1].split("\ndef ")[0]
    assert "verdict(" not in ps, "the paired path must use paired_reading, not verdict()"
    assert "paired_reading(grid, null_ok)" in ps


def test_the_paired_reading_is_scale_free():
    """The delta must be judged against its OWN measured spread, never a frozen constant.
    Deriving the band from one arm's spread is the documented error in
    feedback_tolerance_band_from_the_tested_quantity, whose first instance was this generator."""
    body = _builder().build_body(paired=True)
    pr = body.split("def paired_reading(")[1].split("\ndef ")[0]
    assert "spread = max(vals) - min(vals)" in pr
    assert "mean > spread" in pr
    assert "DELTA_ROBUST" in pr and "DELTA_NOT_ESTABLISHED" in pr
    # The committed delta may appear to fix the SIGN, but must not be used as a threshold.
    assert "mean > 0.0271" not in pr and "mean >= 0.0271" not in pr


def test_the_paired_variant_brands_itself_and_disclaims_the_frozen_arm():
    body = _builder().build_body(paired=True)
    assert "DIAGNOSTIC_METHOD_VARIANT" in body
    assert "method_frozen.epochs = 3" in body
    assert "TRAINING variance" in body, "the training-noise limit must ship in the artifact"


def test_the_frozen_method_constants_are_unchanged_by_the_variant():
    """The variant must change WHAT is scored, never the method being scored."""
    b = _builder()
    for frozen in ("EPOCHS = 3", "LR = 1e-4", "BATCH = 8", 'DTYPE = "float32"',
                   '"r": 16', '"lora_alpha": 32', '"lora_dropout": 0.05',
                   "BAR_BOTH = 0.6679", "BAR_KMER = 0.5990"):
        assert frozen in b.KERNEL, "frozen constant drifted: %s" % frozen
    for recorded in ('"n": 120', '"limit_natural": 1200', '"bins": 6', '"length": 150',
                     '"context": 300'):
        assert recorded in b.KERNEL, "RECORDED geometry drifted: %s" % recorded


def test_the_emitted_paired_kernel_on_disk_matches_a_fresh_build_of_its_code():
    """The committed variant artifact must be the one this builder produces. Compared with the
    base64 payload normalised out: the payload embeds the current source tree, so a full-body
    comparison would fail for reasons unrelated to the code."""
    if not PAIRED_KERNEL.is_file():
        pytest.skip("the paired variant has not been built yet")
    norm = lambda s: re.sub(r"[A-Za-z0-9+/=]{500,}", "<PAYLOAD>", s)
    # The fresh build needs a payload long enough for the normaliser to catch, or this side
    # keeps build_body's short default while the on-disk side gets normalised, and the two
    # differ at exactly that spot for a reason that has nothing to do with the code.
    fresh = _builder().build_body(paired=True, payload="A" * 600)
    assert norm(fresh) == norm(PAIRED_KERNEL.read_text(encoding="utf-8"))
