"""README's "Quickstart (verified output)" block must actually be the output.

WHY THIS EXISTS
That heading is a claim, and on 2026-08-23 it was false in every one of its three examples:

  * the banner read `dna-decode 0.5.0` while the package was at **0.12.1** (seven minor versions stale)
  * the `amr` call omitted `organism: Escherichia` and BOTH trailing lines, including
    `validation: INDEPENDENT_MEASURED -- acc 0.919 (N=8769)`
  * the fungal call likewise omitted its `validation: NO_FREE_PHENOTYPE_SOURCE` line

So the README under-sold the tool's own honesty rails: a reader saw bare `CALL: R` with no trust badge,
which is precisely the impression the project works hardest to avoid giving.

`scripts/verify_quickstart.py` did not catch it — it hardcodes its own STEPS list and reads neither
QUICKSTART.md nor README.md, so its docstring's "executes the exact commands QUICKSTART.md documents" is
a hand-maintained correspondence, not a checked one. Nothing read the README block at all.

WHAT IS ASSERTED: every non-blank line the README SHOWS must appear in the real output, in order. Subset
rather than equality, because the `list` excerpt is deliberately abridged (2 of 44 traits) — an
abridgement is honest, a WRONG line is not. This still catches every defect above, since each was a line
the README displayed that the tool does not produce.

Offline + in-process: the commands are pure-Python paths over committed fixtures (no Docker, no network).
"""
from __future__ import annotations

import contextlib
import io
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

README = REPO / "README.md"
_HEADING = "## Quickstart (verified output)"


def _quickstart_examples() -> list[tuple[str, list[str]]]:
    """-> [(command, [expected output lines])] parsed from the README's verified-output block."""
    text = README.read_text(encoding="utf-8")
    after = text.split(_HEADING, 1)
    if len(after) < 2:
        return []
    m = re.search(r"```text\n(.*?)```", after[1], re.S)
    if not m:
        return []
    examples, cmd, out = [], None, []
    for line in m.group(1).splitlines():
        if line.startswith("$ "):
            if cmd:
                examples.append((cmd, out))
            cmd, out = line[2:].strip(), []
        elif cmd is not None and line.strip():
            out.append(line.rstrip())
    if cmd:
        examples.append((cmd, out))
    return examples


def _run(command: str) -> str:
    """Run a `uv run <console-script> ...` command IN-PROCESS and return stdout."""
    import shlex
    import tomllib

    tokens = shlex.split(command.split("#")[0].strip())      # drop a trailing shell comment
    while tokens and tokens[0] in ("uv", "run"):
        tokens.pop(0)
    route, argv = tokens[0], tokens[1:]

    scripts = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["scripts"]
    spec = scripts.get(route)
    assert spec, f"README quickstart uses {route!r}, which is not a console script"
    import importlib
    mod, _, fn = spec.partition(":")
    main = getattr(importlib.import_module(mod), fn)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = main(argv)
    assert rc in (0, None), f"README quickstart command exited {rc}: {command}"
    return buf.getvalue()


_EXAMPLES = _quickstart_examples()


def test_the_quickstart_block_was_actually_found():
    """Non-vacuity: a renamed heading or fence would empty the parametrized set and pass silently."""
    assert len(_EXAMPLES) >= 3, f"parsed {len(_EXAMPLES)} quickstart examples from README"
    assert all(out for _, out in _EXAMPLES), "an example shows a command but no output"


@pytest.mark.parametrize("command,expected", _EXAMPLES, ids=lambda v: str(v)[:48])
def test_readme_quickstart_lines_are_really_produced(command, expected):
    actual = _run(command).splitlines()
    pos, missing = 0, []
    for want in expected:
        for i in range(pos, len(actual)):
            if actual[i].rstrip() == want:
                pos = i + 1
                break
        else:
            missing.append(want)
    assert not missing, (
        f"README's 'Quickstart (verified output)' shows lines that `{command}` does not produce "
        f"(in order):\n  " + "\n  ".join(missing[:4]) +
        "\nRe-run the command and paste its real output into README.md."
    )


def test_the_version_banner_is_not_stale():
    """The defect that started this: a hardcoded `dna-decode 0.5.0` banner surviving to 0.12.1.

    Called out separately because it is the one line a reader uses to judge whether the whole block is
    current, and a subset match would still flag it -- but not say why.
    """
    from dna_decode.cli import _version
    text = README.read_text(encoding="utf-8")
    banners = set(re.findall(r"^dna-decode (\d+\.\d+\.\d+) - deterministic", text, re.M))
    assert banners <= {_version()}, (
        f"README shows version banner(s) {sorted(banners)} but the package is {_version()}.")


def test_verify_quickstart_covers_every_command_quickstart_md_documents():
    """`scripts/verify_quickstart.py` says it "executes the exact commands QUICKSTART.md documents".

    That correspondence is real today (5 documented commands, 5 STEPS) but it is HAND-MAINTAINED: adding
    a command to QUICKSTART.md does not add a STEP, and the script would keep making the claim while
    silently not verifying it. This pins the claim rather than the current coincidence.

    Identity = (route, the --drug value, or the subcommand). Deliberately coarse: argv cannot be compared
    verbatim because the doc uses `--sample-id demo` where the verifier uses `q`, and that difference is
    cosmetic.
    """
    qs = REPO / "QUICKSTART.md"
    if not qs.exists():
        pytest.skip("QUICKSTART.md absent")

    def identity(tokens: list[str]) -> str | None:
        toks = [t for t in tokens if t not in ("uv", "run")]
        if not toks or not toks[0].startswith("dna-"):
            return None
        if "--drug" in toks:
            return f"{toks[0]}:{toks[toks.index('--drug') + 1]}"
        sub = toks[1] if len(toks) > 1 and not toks[1].startswith("-") else ""
        if sub in ("list",):                     # a pure listing, not a decode step
            return None
        return f"{toks[0]}:{sub}"

    documented = set()
    for line in qs.read_text(encoding="utf-8").splitlines():
        s = line.strip().rstrip("\\").strip()
        if not s.startswith("uv run dna-"):
            continue
        ident = identity(s.split("#")[0].split())
        if ident:
            documented.add(ident)
    assert documented, "parsed no documented commands from QUICKSTART.md"

    import scripts.verify_quickstart as vq
    covered = set()
    for step in vq.STEPS:
        argv = step[2]
        route = "dna-amr" if step[1] is vq.amr_main else "dna-decode"
        if "--drug" in argv:
            covered.add(f"{route}:{argv[argv.index('--drug') + 1]}")
        else:
            covered.add(f"{route}:profile")      # the profile_main steps

    missing = sorted(documented - covered)
    assert not missing, (
        f"QUICKSTART.md documents commands that scripts/verify_quickstart.py does not run: {missing}. "
        f"Add a STEPS entry, or the script's 'executes the exact commands QUICKSTART.md documents' claim "
        f"is false.")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
