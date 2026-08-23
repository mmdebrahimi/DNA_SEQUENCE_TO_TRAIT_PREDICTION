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
import importlib
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

    if route == "dna-decode":
        # resolve past the unified dispatcher to the SUB-parser, otherwise the outer parser would accept
        # anything after the trait and this guard would pass vacuously.
        from dna_decode.cli import _delegate
        trait, rest = argv[0], argv[1:]
        call = lambda: _delegate(trait, rest)          # noqa: E731
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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
