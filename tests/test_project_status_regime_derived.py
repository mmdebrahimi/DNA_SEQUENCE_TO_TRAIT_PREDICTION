"""`project_status.py` must DERIVE the regime map, never restate it.

FOUND 2026-09-22. CLAUDE.md's opening line tells every session to run `scripts/project_status.py`
BEFORE making any claim about this project's scope -- it exists specifically to stop stale prose being
repeated. It was itself printing a hardcoded line:

    constructed variation -> molecular phenotype ..... WORKS (TEM-1 genome-edit, rho 0.761)

while `dna_decode/eval/regime.py`'s note for that exact regime says, verbatim, *"do not quote 0.761 as
the path's general strength"* -- the PEAR replication on a second beta-lactamase came in at 0.352, so
the honest range is 0.35-0.76. The orientation tool was making the one claim its own data layer forbids.

The hardcoded summary was also INCOMPLETE: it listed five regimes where REGIMES holds six, silently
omitting `natural + molecular (zero-shot) -> LOSES_TO_CATALOG` (HIV NNRTI, ESM2 below chance).

Both failures share one cause -- a hand-written summary of a structured source drifts from it -- which
is why the fix is derivation plus this guard, not a corrected literal.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "project_status.py"
sys.path.insert(0, str(ROOT))


def _run() -> str:
    r = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", cwd=ROOT)
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout


def test_the_source_carries_no_hardcoded_regime_figure():
    """THE guard. A literal here is a stale claim waiting to happen -- equality today, drift tomorrow."""
    src = SCRIPT.read_text(encoding="utf-8")
    body = src.split("THE REGIME MAP")[1]
    # COMMENTS ARE EXEMPT, deliberately. The first version of this guard flagged the comments that
    # EXPLAIN the defect ("used to hardcode ... rho 0.761", and the example of the truncation cutting
    # "do not quote 0.761 as "). Those are correct and worth keeping; a rule that fires on correct
    # code gets suppressed, which is how a guard stops being one. What must not reappear is the
    # figure in EXECUTABLE text.
    code = "\n".join(ln for ln in body.splitlines() if not ln.strip().startswith("#"))
    assert "0.761" not in code, "a regime figure is hardcoded again; derive it from REGIMES"
    assert "REGIMES" in code, "the regime map no longer derives from its structured source"


def test_every_regime_in_the_structured_source_is_printed():
    """The hardcoded block listed FIVE regimes where REGIMES holds six, omitting the HIV
    below-chance result entirely. An orientation map that silently drops a regime is worse than none."""
    from dna_decode.eval.regime import REGIMES

    out = _run()
    assert len(REGIMES) >= 6
    for r in REGIMES:
        assert r.verdict in out, r.key
        assert r.population in out and r.endpoint in out, r.key


def test_the_load_bearing_warning_is_not_truncated():
    """A hard [:150] cut this exact caveat mid-sentence at 'do not quote 0.761 as ' -- silently
    amputating the load-bearing half of a warning is the same failure as omitting it."""
    out = _run()
    assert "do not quote 0.761 as the path's general strength" in " ".join(out.split())


def test_it_refuses_rather_than_falling_back_to_a_remembered_summary():
    """If the structured source cannot be read, printing a remembered map is precisely the failure
    this script exists to prevent. The except branch must say UNAVAILABLE, not recite one."""
    src = SCRIPT.read_text(encoding="utf-8")
    body = src.split("THE REGIME MAP")[1]
    assert "UNAVAILABLE" in body
    assert "do NOT quote" in body or "do not quote" in body
