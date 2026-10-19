"""An env var a doc advertises must be one something actually reads.

Why this exists
---------------
`README.md` advertised `export DNA_DECODE_CACHE_ROOT=E:/dna_decode_cache` as a working
knob and **no code had ever read it**, so exporting it did nothing while looking like the
supported way to relocate ~25 GB of caches. `config/datasources.yaml` was honest about the
same variable in the same breath ("would be cleaner -- Phase 3 work"), so the two documents
contradicted each other and the README carried the false one.

This is the `--help` over-claim class that `tests/test_advertised_commands.py` already
guards for CLI surfaces, one level out: **an over-claim in a README is a trust-surface
falsehood, and no report-card, registry or advertised-commands test looks there.**

Scope, deliberately narrow
--------------------------
Only an ASSIGNMENT counts as advertising -- `export VAR=...` or `$env:VAR = ...` inside a
doc. A prose MENTION is exempt by construction, which is what lets the README now discuss
the variable's absence without tripping its own guard (the same structural exemption
`prior_art.py` uses for docstring literals). A bare ALL-CAPS token is NOT matched: this
repo's docs are dense with constants and verdict strings (`DRUG_RULE`, `FAITHFUL_TO_TOOL`,
`INDETERMINATE_INSUFFICIENT_RUNGS`), and a token-level rule over-fires on hundreds of them.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

DOCS = ("README.md", "CLAUDE.md")

#: Three assignment forms, never a mention:
#:   `export VAR=...`            POSIX
#:   `$env:VAR = ...`            PowerShell
#:   `VAR=value <command>`       inline one-shot prefix, at line start
#:
#: The inline form was MISSING from the first version of this guard, and it is the form this
#: repo uses most -- `CLAUDE.md` carries `HF_HOME=D:/hf_cache uv run python ...` three times
#: plus `MSYS_NO_PATHCONV=1 docker ...`. A detector built only for the shape you happen to be
#: fixing under-covers exactly the shapes already in the tree. Measured on adoption: the
#: inline arm adds HF_HOME and MSYS_NO_PATHCONV and nothing else (zero false positives).
ADVERTISED = re.compile(
    r"(?:export\s+([A-Z][A-Z0-9_]{2,})\s*=)"
    r"|(?:\$env:([A-Z][A-Z0-9_]{2,})\s*=)"
    r"|(?m:^[ 	]*\$?[ 	]*([A-Z][A-Z0-9_]{2,})=\S+[ 	]+[a-z])"
)

#: Read by the code under test.
READS = re.compile(
    r"os\.environ\.get\(\s*[\"']([A-Z][A-Z0-9_]+)[\"']"
    r"|os\.environ\[\s*[\"']([A-Z][A-Z0-9_]+)[\"']"
    r"|os\.getenv\(\s*[\"']([A-Z][A-Z0-9_]+)[\"']"
)

#: Variables owned by a third party or the OS, which this repo must NOT implement.
#: Each entry names its owner so the list cannot quietly become an excuse list.
EXTERNAL = {
    "HF_HOME": "huggingface_hub / transformers",
    "UV_CACHE_DIR": "astral-sh/uv",
    "UV_PROJECT_ENVIRONMENT": "astral-sh/uv",
    "MSYS_NO_PATHCONV": "MSYS2 / Git-Bash",
    "PATH": "the OS",
    "PYTHONPATH": "CPython",
}


def _env_vars_read_by_code() -> set[str]:
    found: set[str] = set()
    for sub in ("dna_decode", "scripts", "tests", "tools"):
        base = ROOT / sub
        if not base.is_dir():
            continue
        for py in base.rglob("*.py"):
            # Kaggle kernels embed verbatim copies of this repo's docs as test
            # fixtures, so they are documents wearing a .py extension.
            if "kaggle" in py.parts:
                continue
            try:
                text = py.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            for m in READS.finditer(text):
                found.add(next(g for g in m.groups() if g))
    return found


def _advertised() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for name in DOCS:
        p = ROOT / name
        if not p.exists():
            continue
        for m in ADVERTISED.finditer(p.read_text(encoding="utf-8", errors="ignore")):
            var = m.group(1) or m.group(2) or m.group(3)
            out.setdefault(var, []).append(name)
    return out


def test_every_advertised_env_var_is_read_by_code_or_owned_externally() -> None:
    """A read must EXIST. This cannot prove the variable does anything.

    An unused `os.environ.get("FOO")` sitting in any module satisfies this guard, so it
    catches the advertised-but-never-implemented case (the one that actually happened) and
    NOT the advertised-but-inert case. Scope is also two documents, README.md and CLAUDE.md
    -- not `wiki/`, not `plans/`.
    """
    advertised = _advertised()
    read = _env_vars_read_by_code()
    unimplemented = {
        var: docs
        for var, docs in advertised.items()
        if var not in read and var not in EXTERNAL
    }
    assert not unimplemented, (
        "these docs advertise env vars nothing reads, so setting them does nothing: "
        + "; ".join(f"{v} (in {', '.join(d)})" for v, d in sorted(unimplemented.items()))
        + ". Either implement the read, or correct the doc -- an over-claim in a README "
          "is a trust-surface falsehood."
    )


def test_the_scanner_is_not_vacuous() -> None:
    """A guard that matches nothing would pass no matter how wrong the docs got."""
    advertised = _advertised()
    assert advertised, "the assignment scanner found no advertised env var at all"
    assert "HF_HOME" in advertised, (
        "HF_HOME is advertised in README.md's USB-cache section; if the scanner cannot "
        "see it, the pattern has drifted and the guard above is vacuous"
    )


def test_the_scanner_would_have_caught_the_real_defect() -> None:
    """Replay the exact pre-2026-10-18 README line; the guard must reject it."""
    regression = "export DNA_DECODE_CACHE_ROOT=E:/dna_decode_cache"
    hits = [m.group(1) or m.group(2) for m in ADVERTISED.finditer(regression)]
    assert hits == ["DNA_DECODE_CACHE_ROOT"], hits
    assert "DNA_DECODE_CACHE_ROOT" not in _env_vars_read_by_code(), (
        "if something now reads DNA_DECODE_CACHE_ROOT, this regression fixture is stale "
        "and the README should advertise the variable again"
    )


def test_a_prose_mention_is_exempt_so_the_fix_does_not_trip_its_own_guard() -> None:
    """A MENTION must not count as an advertisement.

    Asserted against a synthetic document rather than the live README, so that
    reintroducing a real assignment fails the guard above with ITS message instead of
    this one -- a misleading failure message costs the next session the time this
    guard exists to save.
    """
    prose = (
        "**There is no `DNA_DECODE_CACHE_ROOT`** -- this README advertised one "
        "until 2026-10-18 and no code has ever read it. "
        "See also DNA_DECODE_CACHE_ROOT in the paragraph below."
    )
    assert not list(ADVERTISED.finditer(prose)), (
        "a prose mention was scored as an advertisement; the exemption is structural "
        "(assignments only) and has broken"
    )
    readme = (ROOT / "README.md").read_text(encoding="utf-8", errors="ignore")
    assert "DNA_DECODE_CACHE_ROOT" in readme, "the README should still explain its absence"


def test_the_narrow_fba_variable_is_not_mistaken_for_a_cache_root() -> None:
    """`DNA_DECODE_CACHE` is real but relocates only the FBA model download."""
    read = _env_vars_read_by_code()
    assert "DNA_DECODE_CACHE" in read, "DNA_DECODE_CACHE should still be read somewhere"
    owners = [
        py
        for py in (ROOT / "dna_decode").rglob("*.py")
        if "DNA_DECODE_CACHE" in py.read_text(encoding="utf-8", errors="ignore")
    ]
    assert [p.name for p in owners] == ["model.py"], (
        "DNA_DECODE_CACHE gained a second reader, so it may no longer be FBA-scoped; "
        f"README.md says it is. Readers: {[str(p) for p in owners]}"
    )


@pytest.mark.parametrize("var,owner", sorted(EXTERNAL.items()))
def test_the_external_allowlist_names_an_owner_for_every_entry(var: str, owner: str) -> None:
    """The allowlist must stay a list of third-party vars, not an excuse list."""
    assert owner and not owner.isspace(), var
    assert var not in _env_vars_read_by_code(), (
        f"{var} is on the EXTERNAL allowlist but this repo reads it; if we own the read, "
        "it belongs in the implemented set, not the allowlist"
    )
