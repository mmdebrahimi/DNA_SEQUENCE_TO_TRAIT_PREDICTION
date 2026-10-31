"""The sampling sweep must ATTEST its swept variable, in the builder and in the emitted kernel.

WHY THIS EXISTS. The 2026-10-29 sweep produced three usable-looking grid points of which only one
was plan-eligible: `scripts/glm_generator_falsifier.py` runs `verify_knob_live` only when
`--verify-knob` is present AND the cell is non-default, and `run_cell`'s argv never passed the flag.
Both T = 0.8 and T = 0.6 therefore scored `excluded:requested-but-unverified` -- two of three points
discarded on a knob that demonstrably DID move the output (0.6391..0.7355 with the null pinned).
A MISSING ATTESTATION, not an inert knob, and nothing failed: no test covered this builder at all.

The predecessor defect is why the attestation matters in the first place: the 2026-10-19 sweep
reported three identical numbers from three grid points as an INVARIANCE FINDING, when the knob was
simply inert. An unattested sweep cannot tell those two apart.

WHAT THIS GUARDS, and the second assertion is the load-bearing one. The emitted kernel is a tracked
generated file, so the builder and the artifact can drift apart -- and a stale COPY of this kernel
is a recorded failure mode here (the kernel's own refusal text tells the operator that the notebook
may hold an older copy than the repo). Checking only the builder would pass while the committed
kernel that actually runs still omits the flag.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "kaggle" / "build_glm_sampling_sweep_kernel.py"
EMITTED = ROOT / "scripts" / "kaggle" / "glm_sampling_sweep_kernel.py"

FLAG = "--verify-knob"


def _run_cell_argv(source: str, where: str) -> list[str]:
    """The string literals of `run_cell`'s `argv` list -- parsed, never regexed.

    A substring search for the flag anywhere in the file would also match a comment or the
    docstring, so it could pass on a file that merely TALKS about the attestation. The flag has to
    be in the list that is handed to the CLI.
    """
    tree = ast.parse(source)
    fns = [n for n in ast.walk(tree)
           if isinstance(n, ast.FunctionDef) and n.name == "run_cell"]
    assert len(fns) == 1, "%s: expected exactly one run_cell, found %d" % (where, len(fns))

    for node in ast.walk(fns[0]):
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "argv" for t in node.targets):
            assert isinstance(node.value, ast.List), "%s: argv is not a list literal" % where
            return [e.value for e in node.value.elts
                    if isinstance(e, ast.Constant) and isinstance(e.value, str)]
    raise AssertionError("%s: run_cell has no `argv = [...]` assignment" % where)


def _kernel_template(builder_source: str) -> str:
    """The KERNEL raw-string template the builder emits, so `run_cell` can be parsed from it."""
    m = re.search(r"^KERNEL = r'''(.*?)'''", builder_source, re.S | re.M)
    assert m, "KERNEL = r'''...''' template not found in the builder"
    return m.group(1)


def test_builder_template_passes_the_verify_knob_flag():
    argv = _run_cell_argv(_kernel_template(BUILDER.read_text(encoding="utf-8")), "builder template")
    assert FLAG in argv, (
        "run_cell's argv omits %s, so every non-default grid point will score "
        "`excluded:requested-but-unverified` and be discarded from the verdict. "
        "argv was: %s" % (FLAG, argv))


def test_the_emitted_kernel_that_actually_runs_also_passes_it():
    assert EMITTED.is_file(), "%s is absent -- rebuild it with the builder" % EMITTED
    argv = _run_cell_argv(EMITTED.read_text(encoding="utf-8"), "emitted kernel")
    assert FLAG in argv, (
        "the COMMITTED kernel omits %s even if the builder carries it -- the emitted file is a "
        "tracked artifact and has drifted. Re-run "
        "`uv run python scripts/kaggle/build_glm_sampling_sweep_kernel.py`. argv was: %s"
        % (FLAG, argv))


def test_the_flag_is_unconditional_so_a_new_grid_point_cannot_ship_unattested():
    """It must sit in the argv literal, not behind an `if temperature != 1.0` style guard.

    The falsifier already self-gates on `requested_non_default`, so the default (1.0, 0) self-check
    cell skips the probe at no cost. Passing it conditionally here would buy nothing and would let a
    future grid point be added un-attested by forgetting to extend the condition -- which is the
    shape of the defect this file exists for.
    """
    argv = _run_cell_argv(_kernel_template(BUILDER.read_text(encoding="utf-8")), "builder template")
    assert argv.count(FLAG) == 1, "expected exactly one %s in argv, got %d" % (FLAG, argv.count(FLAG))
    # Present in the flat list literal at all => not guarded by a conditional, since a conditional
    # would have to live outside the literal (an append, or an `if` around the list).
    assert argv.index(FLAG) < argv.index("--out"), (
        "%s must precede --out in argv, the same ordering the falsifier's CLI is exercised with"
        % FLAG)
