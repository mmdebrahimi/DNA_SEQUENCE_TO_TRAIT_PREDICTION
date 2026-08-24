"""Every command the decode router ADVERTISES must actually be runnable.

WHY THIS EXISTS
`dna_decode/decode_router.py` is the "point at your DNA -> what can this tool decode?" coherence layer: it
hands the user an exact command to run for each applicable decoder. That makes each `example` string a
PROMISE, and on 2026-08-23 four of them were false:

  * `dna-clinvar --vcf sample.vcf`                  -> unrecognized arguments (the VCF is POSITIONAL)
  * `dna-hla --vcf sample.vcf`                      -> unrecognized arguments (same)
  * `dna-pgx --gene CYP2C19 --vcf sample.vcf`       -> invalid choice (the choices are LOWERCASE)
  * `dna-forward --mutation ... --genome-fasta ...` -> unrecognized arguments; `dna-forward` had NO
                                                       nucleotide input path at all

The first three were syntax drift. The fourth was worse in kind: it advertised a CAPABILITY that shipped
as a library (`forward/genome_edit.predict_genome_edit`, tested and validated) but had never been wired to
a CLI. It is now real (`--cds-fasta` / `--cds-seq`), which is why this file asserts *validity*, not the
absence of a flag -- the right fix for an advertisement of something real is to build it.

Nothing caught these because the router's own tests only checked that the table renders. A rendered string
is not a runnable command.

HOW IT WORKS (and why it is safe + fast)
`argparse.ArgumentParser.parse_args` is patched to raise a sentinel the instant a parse SUCCEEDS. So each
command is resolved through its real parser and then stops -- no file is read, no network, no Docker, no
model load. A parse FAILURE surfaces as argparse's own `SystemExit(2)`, which is the assertion. This tests
exactly the argument contract and nothing beyond it.

SCOPE: this proves the commands PARSE. It does not prove they produce a good answer -- that is each cell's
own evidence contract (`cell_registry`) and validation artifact.
"""
from __future__ import annotations

import argparse
import contextlib
import importlib
import io
import re
import shlex
import sys
import tomllib
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode import decode_router as router  # noqa: E402

_PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"


class _ParsedOK(Exception):
    """Raised the moment a parser accepts the argv -- stops before any real work happens."""


def _console_scripts() -> dict[str, str]:
    return tomllib.loads(_PYPROJECT.read_text(encoding="utf-8"))["project"]["scripts"]


def _resolve(route: str):
    """route ('dna-mlst') -> its console-entry callable. Fails if the route is not a real script."""
    spec = _console_scripts().get(route)
    assert spec, (f"decode_router advertises route {route!r}, which is NOT a console script in "
                  f"pyproject [project.scripts] -- a user could not run it at all.")
    mod, _, fn = spec.partition(":")
    return getattr(importlib.import_module(mod), fn)


def _assert_parses(command: str) -> None:
    """Assert `command` is accepted by its real parser. Stops at the parse; runs nothing."""
    tokens = shlex.split(command)
    route, argv = tokens[0], tokens[1:]

    if route == "dna-decode" and argv and not argv[0].startswith("-"):
        # resolve past the unified dispatcher to the SUB-parser, otherwise the outer parser would accept
        # anything after the trait and this guard would pass vacuously.
        from dna_decode.cli import TRAITS, _delegate
        trait, rest = argv[0], argv[1:]
        if trait not in TRAITS:                        # an ANALYSES command or `list` -> the unified main
            call = lambda: _resolve("dna-decode")(argv)    # noqa: E731
        else:
            call = lambda: _delegate(trait, rest)      # noqa: E731
    elif route == "dna-decode":                        # a bare flag, e.g. `dna-decode --version`
        call = lambda: _resolve("dna-decode")(argv)    # noqa: E731
    else:
        main = _resolve(route)
        call = lambda: main(argv)                      # noqa: E731

    orig_parse = argparse.ArgumentParser.parse_args

    def _stop(self, args=None, namespace=None):
        orig_parse(self, args, namespace)     # the REAL parse_args -- full validation, then stop
        raise _ParsedOK

    # ONLY parse_args is patched, deliberately. `parse_args` IS `parse_known_args` plus the
    # leftover-argument check, so patching `parse_known_args` to raise-on-success would skip exactly the
    # "unrecognized arguments" error that 3 of the 4 original bugs produced -- the guard passed all three
    # until this was removed. No CLI in the package calls `parse_known_args` directly.
    argparse.ArgumentParser.parse_args = _stop
    try:
        call()
    except _ParsedOK:
        return                                          # parsed cleanly -- the promise holds
    except SystemExit as e:
        if e.code not in (0, None):
            raise AssertionError(
                f"decode_router advertises a command that its own CLI REJECTS:\n    {command}\n"
                f"argparse exited {e.code} (see the captured stderr for 'unrecognized arguments' / "
                f"'invalid choice' / 'required'). Fix the router example, or build the capability it "
                f"promises."
            ) from None
    except Exception:
        return          # got PAST parsing and failed on real work (missing file etc.) -- contract is fine
    finally:
        argparse.ArgumentParser.parse_args = orig_parse


def _advertised() -> list[tuple[str, str]]:
    return [(kind, d.example) for kind, decs in router.DECODERS.items() for d in decs]


def test_the_router_advertises_something_for_every_input_kind():
    """Guards against a vacuous pass: an empty table would satisfy every other test in this file."""
    assert set(router.DECODERS) == {"vcf_human", "protein_fasta", "nucleotide_fasta"}
    for kind, decs in router.DECODERS.items():
        assert decs, f"no decoders advertised for {kind}"
    assert len(_advertised()) >= 18


@pytest.mark.parametrize("kind,command", _advertised(), ids=lambda v: v.split()[0] if " " in str(v) else v)
def test_every_advertised_command_is_accepted_by_its_own_cli(kind, command):
    _assert_parses(command)


def test_every_advertised_route_is_a_real_console_script():
    scripts = _console_scripts()
    missing = sorted({d.route for decs in router.DECODERS.values() for d in decs} - set(scripts))
    assert not missing, f"decode_router advertises routes with no console entry: {missing}"


def test_commands_printed_by_the_run_path_are_also_valid(tmp_path, capsys):
    """`run_decode_plan` PRINTS extra commands that are not in the DECODERS table (the human-VCF hints).

    They are the same kind of promise and drifted the same way -- these three carried the `--vcf` and
    `CYP2C19` bugs, plus stale text calling clinvar/hla 'standalone script' after both became routable
    `dna-decode` traits. Parametrized tests over a static table would never have seen them.
    """
    vcf = tmp_path / "sample.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")
    router.run_decode_plan(vcf, runner=lambda *a, **k: 0)
    out = capsys.readouterr().out

    cmds = [ln.split("- ", 1)[1].strip() for ln in out.splitlines()
            if ln.strip().startswith("- dna-")]
    assert len(cmds) == 3, f"expected the 3 human-VCF hints, got {cmds}"
    for c in cmds:
        # strip a trailing parenthetical annotation before parsing
        _assert_parses(c.split("   (")[0].strip())
    assert "standalone script" not in out, (
        "clinvar/hla are routable `dna-decode` traits now -- the 'standalone script' text is stale.")


# --------------------------------------------------------------------------------------------------
# The SAME promise, one level down: every CLI module's own docstring / --help examples.
#
# The router is not the only place that hands a user a command. Each decoder's module docstring opens with
# a worked example block, and those reach the user through `--help` too. Sweeping them found two more:
#   * `dna-decode inverse --protein-fasta X --cds-fasta Y`  -- omits the REQUIRED --target-percentile
#   * `dna-hla cohort.vcf --allele b5801 ...`               -- b5801 is not a shipped allele
# The second was the worse one: B*58:01/allopurinol and A*31:01/carbamazepine were MEASURED against 1000G
# HLA truth and DEMOTED (sens 0.61 / PPV 0.18, and sens 0.0 -- see catalog._UNVALIDATED_TAGS), yet the
# docstring example AND the argparse description still advertised those two screens. An over-claim in
# --help is a trust-surface falsehood in the most-read place there is.
# --------------------------------------------------------------------------------------------------

# A docstring line is treated as a literal command only if it has no placeholder/prose markers. This is
# deliberately conservative -- `<...>`, `...`, quotes, pipes and `$` mean "fill this in", and a template is
# not a promise. The count assertion below stops the filter from quietly excluding everything.
_PLACEHOLDER = ("...", "<", "|", "`", "'", '"', "$")
_CMD_LINE = re.compile(r"^\s*(dna-[a-z0-9-]+ [^\n#]*?)\s*(?:#.*)?$")


def _docstring_commands() -> list[tuple[str, str]]:
    root = Path(__file__).resolve().parent.parent
    out, seen = [], set()
    for spec in sorted(set(_console_scripts().values())):
        mod = spec.split(":")[0]
        path = root / (mod.replace(".", "/") + ".py")
        if not path.exists():
            continue
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            m = _CMD_LINE.match(line)
            if not m:
                continue
            cmd = m.group(1).strip()
            if any(t in cmd for t in _PLACEHOLDER) or cmd in seen:
                continue
            seen.add(cmd)
            out.append((mod, cmd))
    return out


_DOCSTRING_COMMANDS = _docstring_commands()


def test_the_docstring_sweep_actually_finds_commands():
    """Non-vacuity: if the placeholder filter ever swallowed the whole surface, every parametrized case
    below would vanish and the sweep would silently pass while checking nothing."""
    assert len(_DOCSTRING_COMMANDS) >= 60, f"docstring sweep collapsed to {len(_DOCSTRING_COMMANDS)} commands"
    assert len({m for m, _ in _DOCSTRING_COMMANDS}) >= 15


@pytest.mark.parametrize("mod,command", _DOCSTRING_COMMANDS, ids=lambda v: str(v)[:60])
def test_every_command_in_a_cli_docstring_is_accepted_by_its_own_cli(mod, command):
    _assert_parses(command)


def test_every_shipped_help_text_survives_a_legacy_console_encoding():
    """`--help` must not CRASH on a Windows console. Measured 2026-08-23, not assumed.

    On this host `sys.stdout` is **cp1252 with errors=surrogateescape**, so a plain `print()` of a
    character outside cp1252 raises UnicodeEncodeError -- it does not degrade to mojibake. Empirically
    `-> (U+2192)`, `>= (U+2265)`, `<= (U+2264)`, Greek letters, `~= (U+2248)` and `!= (U+2260)` all
    CRASH, while em-dash, en-dash, curly quotes, bullet and `x (U+00D7)` are all IN cp1252 and are fine.
    So this is not "avoid Unicode" -- it is a specific, small, checkable set.

    The shipped console surface was swept and is CLEAN (0 of 48 broken), which is exactly why this guard
    is worth having: it pins a property that currently holds. 24 internal `scripts/*.py` modules DO crash
    on `--help` (an arrow in the module docstring that argparse uses as its description); those are
    developer analysis scripts, not the shipped product, and are deliberately out of scope here.

    Checked by rendering each parser's real help text and encoding it -- no subprocess, no console needed.
    """
    scripts = _console_scripts()
    broken = []
    for name, spec in sorted(scripts.items()):
        main = _resolve(name)
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
                main(["--help"])
        except SystemExit:
            pass
        except Exception:                       # a main that cannot even render help is another test's job
            continue
        text = buf.getvalue()
        if not text:
            continue
        try:
            text.encode("cp1252")
        except UnicodeEncodeError as e:
            offender = text[e.start:e.end]
            broken.append(f"{name} (U+{ord(offender[0]):04X})")
    assert not broken, (
        "these shipped console scripts print a `--help` that CRASHES on a cp1252 console: "
        f"{broken}. Replace the character with its ASCII form (-> for an arrow, >= for U+2265); "
        "em/en-dashes and curly quotes are safe."
    )


def test_hla_help_does_not_advertise_the_demoted_tags():
    """The demotion must reach the user-facing surface, not only the catalog's internal comment.

    `catalog._UNVALIDATED_TAGS` records that rs9263726 misses 39% of B*58:01 carriers ("unsafe for an
    SJS/TEN screen") and rs1061235 cannot call A*31:01 at all -- while `dna-hla --help` claimed
    "abacavir/allopurinol/carbamazepine". Only the abacavir screen ships.
    """
    from dna_decode.hla import HLA_ALLELES
    from dna_decode.hla import cli as hla_cli
    assert set(HLA_ALLELES) == {"b5701"}

    err = io.StringIO()
    with contextlib.redirect_stdout(err), pytest.raises(SystemExit):
        hla_cli.main(["--help"])
    helptext = err.getvalue().lower()
    for demoted in ("b5801", "a3101"):
        assert demoted not in helptext, f"--help advertises the DEMOTED tag {demoted}"
    # the drug names may appear only while saying they are NOT shipped
    for drug in ("allopurinol", "carbamazepine"):
        if drug in helptext:
            assert "failed" in helptext or "not shipped" in helptext, (
                f"--help mentions {drug} without saying its tag failed validation")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))


# --------------------------------------------------------------------------------------------------
# Flags named in PROSE (not in a full command line).
#
# The command-line guards above parse whole invocations, so they are blind to a flag mentioned in
# running text -- which is how `docs/quickstart.md` and `examples/README.md` both came to advertise
# `--json` for `dna-amr` when that CLI exposes only `--json-only` (2026-08-24). The union check below
# then found a 4th instance: README advertised `--legacy-6class` for `dna-pathotype`, a flag that was
# never implemented anywhere in the repo.
#
# SCOPE, honestly: this asserts every backticked flag in the docs is DECLARED SOMEWHERE in the repo.
# It cannot tell that a flag real on one tool is wrong for the tool being described (`--json` IS real
# on `dna-forward`). That attribution problem is deliberately NOT solved by a proximity heuristic --
# a noisy guard gets disabled. The `--json` case is pinned explicitly below instead.
# --------------------------------------------------------------------------------------------------

_DOC_FILES = ("README.md", "QUICKSTART.md", "CLAUDE.md", "docs/quickstart.md", "examples/README.md")
# flags that are not ours: argparse builtins + flags of external tools the docs legitimately show
_FOREIGN_FLAGS = {"--help", "--version", "--rm", "--entrypoint", "--output", "--type", "-v"}


def _declared_flags() -> set[str]:
    """Every long flag declared anywhere in the repo, by static scan.

    Static rather than by import: importing every analysis script is heavy and side-effecty. The
    trade-off is that a flag handled by manual argv inspection instead of `add_argument` is invisible
    here -- `dna-decode decode --run` is exactly that -- so those are allowlisted explicitly.
    """
    import re
    root = Path(__file__).resolve().parent.parent
    pat = re.compile(r"add_argument\(\s*[\"'](--[a-z][a-z0-9-]*)")
    out: set[str] = set()
    for sub in ("scripts", "dna_decode"):
        for p in (root / sub).rglob("*.py"):
            out |= set(pat.findall(p.read_text(encoding="utf-8", errors="replace")))
    out.add("--run")        # dna-decode decode --run: parsed by hand in cli.py, not via add_argument
    return out


def test_every_flag_named_in_the_docs_is_declared_somewhere():
    import re
    root = Path(__file__).resolve().parent.parent
    declared = _declared_flags() | _FOREIGN_FLAGS
    pat = re.compile(r"`(--[a-z][a-z0-9-]*)`")
    undeclared = {}
    for doc in _DOC_FILES:
        p = root / doc
        if not p.exists():
            continue
        named = set(pat.findall(p.read_text(encoding="utf-8", errors="replace")))
        missing = sorted(named - declared)
        if missing:
            undeclared[doc] = missing
    assert not undeclared, (
        f"docs name flags that are declared NOWHERE in the repo: {undeclared}. Either implement the "
        f"flag or stop advertising it -- an advertised flag is a promise.")


def test_the_docs_do_not_advertise_json_for_the_amr_cli():
    """REGRESSION: `dna-amr` exposes `--json-only`, NOT `--json`.

    `--json` IS a real flag on `dna-forward`, so the union check above cannot catch this -- the flag
    exists, just not on the tool these files describe. Pinned explicitly because it shipped wrong in
    two files at once.
    """
    root = Path(__file__).resolve().parent.parent
    from dna_decode.amr import cli as amr_cli
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), pytest.raises(SystemExit):
        amr_cli.main(["--help"])
    helptext = buf.getvalue()
    assert "--json-only" in helptext and "--json " not in helptext.replace("--json-only", "")

    for doc in ("docs/quickstart.md", "examples/README.md"):
        text = (root / doc).read_text(encoding="utf-8", errors="replace")
        assert "`--json`" not in text, (
            f"{doc} advertises `--json`; the amr CLI exposes `--json-only`.")
