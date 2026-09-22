"""A shipped CLI's `--help` must render on this host's DEFAULT console, not only under PYTHONUTF8=1.

FOUND BY USING IT (2026-09-22): `uv run python -m scripts.fetch_prospective_cohort --help` died with
`UnicodeEncodeError: 'charmap' codec can't encode character '\\u2192'` -- argparse wrote the help text to
a cp1252 console and a `->` arrow in the module docstring killed it. A traceback instead of help is a
real usability defect: it hits a user at the exact moment they are trying to discover how to run the
thing. Measured across the tree: 24 of 336 argparse scripts crashed this way.

WHY THE EXISTING GUARD DID NOT CATCH IT. `tests/test_advertised_commands.py` resolves every advertised
command through its REAL parser, but it does so by patching `parse_args` to raise the instant a parse
SUCCEEDS -- so it never lets argparse RENDER anything. `--help` is precisely the path that renders. The
two guards are complements: that one checks a command PARSES, this one checks its help PRINTS.

SCOPE, STATED HONESTLY. This is a STATIC approximation: it checks the text argparse would print that we
control -- the module docstring (when it is passed as `description=__doc__`) and every `help=` string --
is encodable in the legacy console codepage. It does NOT drive argparse itself, so it does not cover
`prog`, `metavar` or `choices` values. That is a deliberate trade: importing 336 scripts to render help
for real would pull in torch and every heavy dependency in the tree. All 24 real failures lived in the
text this DOES cover.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"

# The console codepage this host actually defaults to. The point is NOT that cp1252 is good -- it is that
# a CLI must not assume a UTF-8 console, because on Windows it does not get one by default.
LEGACY_CODEPAGE = "cp1252"


def _argparse_scripts() -> list[Path]:
    """DISCOVERED from the tree, never hand-listed -- a hand-maintained list drifts the moment a script
    is added, which is the documented `hardcoded_exclusion_list_undercovers` failure (5 instances)."""
    return [p for p in sorted(SCRIPTS.glob("*.py"))
            if "argparse" in p.read_text(encoding="utf-8", errors="replace")]


def _help_reachable_text(path: Path) -> str:
    """The help text we control: the module docstring IF it is used as the description, plus every
    `help=` string. Reading the docstring unconditionally would over-report -- a script whose docstring
    is never passed to argparse can hold whatever it likes."""
    txt = path.read_text(encoding="utf-8", errors="replace")
    out = []
    if "description=__doc__" in txt:
        m = re.match(r'\s*"""(.*?)"""', txt, re.S)
        if m:
            out.append(m.group(1))
    out.extend(re.findall(r'help=["\'](.*?)["\']', txt, re.S))
    return "\n".join(out)


def test_there_are_argparse_scripts_to_check():
    """NON-VACUITY. If the discovery predicate silently matched nothing, every assertion below would
    pass while checking zero files -- a green test that guards nothing."""
    assert len(_argparse_scripts()) > 100


def test_every_cli_help_text_renders_on_the_legacy_console():
    """THE guard. A `--help` that raises is worse than a bad `--help`."""
    broken = []
    for p in _argparse_scripts():
        try:
            _help_reachable_text(p).encode(LEGACY_CODEPAGE)
        except UnicodeEncodeError as e:
            bad = _help_reachable_text(p)[e.start]
            broken.append(f"{p.name}: {bad!r} (U+{ord(bad):04X})")
    assert not broken, (
        "these scripts would raise UnicodeEncodeError when argparse prints --help on a "
        f"{LEGACY_CODEPAGE} console:\n  " + "\n  ".join(broken))


def test_the_guard_actually_detects_the_character_that_caused_this():
    """NON-VACUITY for the check itself: the exact character that broke fetch_prospective_cohort must
    still be rejected by the encoder this test relies on. If cp1252 ever gained it, the guard would be
    silently inert."""
    with pytest.raises(UnicodeEncodeError):
        "→".encode(LEGACY_CODEPAGE)


def test_ascii_replacements_are_what_the_fixed_scripts_now_use():
    """The 24 scripts were fixed by substituting ASCII equivalents in help-reachable text, not by
    suppressing the error. Spot-check the one that surfaced the defect."""
    txt = (SCRIPTS / "fetch_prospective_cohort.py").read_text(encoding="utf-8")
    m = re.match(r'\s*"""(.*?)"""', txt, re.S)
    assert m and "→" not in m.group(1)
