"""A typing CLI's `--identity` / `--coverage` default must BE the module constant, not a copy of it.

THE LIVE DEFECT THIS CLOSES. `dna-salmserovar` shipped `--coverage default=80.0` from the cell's first
commit. On 2026-09-04 `SEROVAR_COVERAGE_THRESHOLD` was lowered 80 -> 40 against a pre-registered bar,
measured at +35 net correct calls -- and the CLI default was not updated. Because `main()` passes
`args.coverage` explicitly, the shipped entry point kept overriding the validated value with the old
one. Nothing caught it because every validation (the threshold sweep, the SeqSero2 comparison, the O
port) called `call_serovar` directly, which uses the constant. The seam between the validated function
and the shipped command had no test.

Same class as `tests/test_advertised_commands.py`: what ships and what was measured must be the same
thing. That guard resolves advertised commands through their real parser; this one resolves shipped
DEFAULTS through their real parser.

SCOPE HONESTY: salmserovar was the ONLY one of nine typing cells that had drifted -- every sibling's
restated default happened to still match, because their constants were never revised. The fix is
structural (derive, don't restate) so a future constant change cannot reintroduce it in any of them.
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# The pre-fix hardcoded values, recorded so the six no-behaviour-change edits are PROVABLE rather than
# asserted: any of these that no longer matches its constant means that edit changed behaviour.
PRE_FIX_DEFAULTS = {
    "serotype": (85.0, 60.0),
    "pneumoserotype": (90.0, 70.0),
    "ktype": (90.0, 80.0),
    "plasmid": (95.0, 60.0),
    "resfinder": (90.0, 60.0),
    "disinfinder": (90.0, 60.0),
}
# salmserovar is deliberately EXCLUDED above: its coverage default genuinely changed 80.0 -> 40.0.
SALMSEROVAR_INTENDED_CHANGE = {"identity": 90.0, "coverage": 40.0}


def _typing_clis() -> list[str]:
    """DISCOVERED, not hand-listed -- a hand-enumerated list beside the data that defines it drifts,
    which is the very failure this file exists for."""
    found = []
    for pkg in sorted(p.name for p in (ROOT / "dna_decode").iterdir() if p.is_dir()):
        cli = ROOT / "dna_decode" / pkg / "cli.py"
        if not cli.exists():
            continue
        src = cli.read_text(encoding="utf-8")
        if '"--identity"' in src or '"--coverage"' in src:
            found.append(pkg)
    return found


def _defaults(pkg: str) -> dict[str, float]:
    """Resolve the defaults through the REAL parser, so this cannot pass on a source-text match while
    the shipped command does something else."""
    mod = importlib.import_module(f"dna_decode.{pkg}.cli")
    captured: dict[str, float] = {}
    real_parse = argparse.ArgumentParser.parse_args

    def spy(self, *a, **k):
        for act in self._actions:
            for opt in act.option_strings:
                if opt in ("--identity", "--coverage") and act.default is not None:
                    captured[opt.lstrip("-")] = float(act.default)
        raise SystemExit(0)      # stop before the command does any work

    argparse.ArgumentParser.parse_args = spy
    try:
        try:
            mod.main(["dummy.fna"])
        except SystemExit:
            pass
    finally:
        argparse.ArgumentParser.parse_args = real_parse
    return captured


def _constants(pkg: str) -> dict[str, float]:
    runner = importlib.import_module(f"dna_decode.{pkg}.runner")
    out: dict[str, float] = {}
    for name in dir(runner):
        val = getattr(runner, name)
        if not isinstance(val, float):
            continue
        if "IDENTITY" in name:
            out["identity"] = val
        elif "COVERAGE" in name:
            out["coverage"] = val
    return out


def test_discovery_is_not_vacuous():
    """A zero-length discovery would make every test below pass while checking nothing."""
    pkgs = _typing_clis()
    assert len(pkgs) >= 6, pkgs
    assert "salmserovar" in pkgs


@pytest.mark.parametrize("pkg", _typing_clis())
def test_cli_default_equals_the_module_constant(pkg: str):
    """The seam. A shipped default that merely COPIES a constant is one revision away from being a
    lie about what was validated."""
    got, want = _defaults(pkg), _constants(pkg)
    assert got, f"{pkg}: no --identity/--coverage default resolved through the real parser"
    for flag, value in got.items():
        assert flag in want, f"{pkg}: --{flag} has no corresponding runner constant"
        assert value == want[flag], (
            f"{pkg}: CLI --{flag} default {value} != runner constant {want[flag]}. The shipped command "
            f"is not running the validated threshold.")


def test_salmserovar_now_ships_the_validated_coverage_threshold():
    """The live defect, named explicitly rather than left to the generic check: 80.0 was the pre-2026-09-04
    value and it reached users for five days after the change that replaced it."""
    got = _defaults("salmserovar")
    assert got["coverage"] == SALMSEROVAR_INTENDED_CHANGE["coverage"] == 40.0
    assert got["identity"] == SALMSEROVAR_INTENDED_CHANGE["identity"] == 90.0
    assert got["coverage"] != 80.0


@pytest.mark.parametrize("pkg", sorted(PRE_FIX_DEFAULTS))
def test_the_other_cells_kept_their_exact_previous_defaults(pkg: str):
    """Six CLIs were rewritten to derive their defaults instead of restating them. That edit must be a
    pure refactor -- if any value moved, it changed a shipped caller's behaviour without a bar."""
    got = _defaults(pkg)
    want_id, want_cov = PRE_FIX_DEFAULTS[pkg]
    assert got["identity"] == want_id, f"{pkg}: identity default MOVED {want_id} -> {got['identity']}"
    assert got["coverage"] == want_cov, f"{pkg}: coverage default MOVED {want_cov} -> {got['coverage']}"


@pytest.mark.parametrize("pkg", _typing_clis())
def test_the_default_is_not_a_hardcoded_literal(pkg: str):
    """Equality today is not enough -- a restated literal passes the check above and drifts tomorrow.
    The source must reference the constant, which makes the property true by construction."""
    src = (ROOT / "dna_decode" / pkg / "cli.py").read_text(encoding="utf-8")
    for flag in ("identity", "coverage"):
        marker = f'"--{flag}", type=float, default='
        if marker not in src:
            continue
        tail = src.split(marker, 1)[1].lstrip()
        assert not tail[0].isdigit(), (
            f"{pkg}: --{flag} default is a hardcoded literal; reference the runner constant instead")
