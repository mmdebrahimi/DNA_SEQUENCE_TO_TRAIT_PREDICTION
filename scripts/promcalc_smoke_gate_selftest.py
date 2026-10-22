# -*- coding: utf-8 -*-
"""Self-test for stage B's bounded smoke gate: does it actually CATCH the version breaks?

    python scripts/promcalc_smoke_gate_selftest.py        # exit 0 = gate is non-vacuous

**Why this is a script and not a pytest test.** It needs the real `promotercalculator` in the isolated
stage-B interpreter and spends ~3 minutes of CPU doing real decompositions. A unit test cannot
establish that a third-party library still behaves -- that is the whole point of a real-surface check
-- and putting it in the suite would make the suite depend on an environment the project env does not
have.

**Why it exists at all.** The first write-up of the smoke gate CITED an earlier run rather than
verifying the claim, and re-running it found a figure wrong: "caught both breaks in 24 s" held only for
the IMPORT break. The container-shape break fails only AFTER the expensive solve, so on the full 44,106
tiles it would have burned the entire ~12.7 h run and returned an empty table. A gate that has only
ever been observed PASSING proves nothing; this makes the non-vacuity re-runnable instead of remembered.

Safe to run while the real job is live: everything happens on throwaway copies in a temp directory, and
neither the live script nor the live run output directory is touched. It does compete for CPU, so it
uses 4 workers rather than 7.

The claim under test is not "the smoke run passes" but "the smoke gate CATCHES the two version
breaks". A gate that passes proves nothing on its own -- it has to fail on the broken input, or it is
decoration. So this runs:

  A. current code            -> expect PASS (73.0 sites/tile, 0 errors)
  B. break 1 re-introduced   -> expect every tile to fail (ImportError)
  C. break 2 re-introduced   -> expect every tile to fail (AttributeError)

Everything runs on COPIES in a temp dir. The live script and the live run's output directory are
never touched.
"""
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(r"C:\Users\Laptop Portable\pythonprojects\dna_decode")
PYISO = Path(r"C:\Users\Laptop Portable\tools\promcalc_env\Scripts\python.exe")
STAGE_B = REPO / "scripts" / "promcalc_decompose.py"
WORKERS = "4"          # the box is busy with the real 7-worker run; keep the extra load modest
GUARD = 'if __name__ == "__main__":'

# Break 1: exactly the pre-fix import -- the top-level path only, no submodule fallback.
BREAK1 = '''

def _load_calculator():
    from promoter_calculator import Promoter_Calculator
    return Promoter_Calculator

'''

# Break 2: exactly the pre-fix field access -- .get on an object that does not define it.
BREAK2 = '''

def _site_field(d, name, default=None):
    return d.get(name, default)

'''


def run_case(label: str, work: Path, src: Path, expect: str) -> dict:
    out = work / "out"
    out.mkdir(exist_ok=True)
    shutil.copy(REPO / "data" / "promcalc" / "tiles_smoke.tsv", out / "tiles_smoke.tsv")
    t0 = time.time()
    p = subprocess.run([str(PYISO), str(src), "--in-dir", str(out), "--smoke", "--workers", WORKERS],
                       capture_output=True, text=True)
    el = time.time() - t0
    txt = p.stdout
    rec = {}
    if "{" in txt:
        try:
            rec = json.loads(txt[txt.index("{"): txt.rindex("}") + 1])
        except Exception:
            rec = {}
    ok, err = rec.get("ok"), rec.get("err")
    sample = (rec.get("errors_sample") or [None])[0]
    kind = ""
    if sample:
        m = re.search(r"(\w*Error)", str(sample))
        kind = m.group(1) if m else ""
    print(f"\n--- {label} ---")
    print(f"  expect              : {expect}")
    print(f"  elapsed             : {el:.1f}s")
    print(f"  tiles ok            : {ok}")
    print(f"  tiles errored       : {err}")
    print(f"  sites/tile mean     : {rec.get('sites_per_tile_mean')}")
    print(f"  first error kind    : {kind or '(none)'}")
    if sample:
        print(f"  first error text    : {str(sample)[:150]}")
    return {"ok": ok, "err": err, "spt": rec.get("sites_per_tile_mean"), "kind": kind}


results = {}
with tempfile.TemporaryDirectory() as td:
    base = Path(td)
    src_text = STAGE_B.read_text(encoding="utf-8")
    assert src_text.count(GUARD) == 1

    # A -- current code, unmodified copy
    wa = base / "A"; wa.mkdir()
    sa = wa / "stage_b_current.py"
    sa.write_text(src_text, encoding="utf-8")
    results["A_current"] = run_case("A. CURRENT CODE", wa, sa, "PASS: ~73.0 sites/tile, 0 errors")

    # B -- break 1 re-introduced
    wb = base / "B"; wb.mkdir()
    sb = wb / "stage_b_break1.py"
    sb.write_text(src_text.replace(GUARD, BREAK1 + GUARD), encoding="utf-8")
    results["B_break1_import"] = run_case("B. BREAK 1 (top-level class import only)", wb, sb,
                                          "FAIL: every tile ImportError")

    # C -- break 2 re-introduced
    wc = base / "C"; wc.mkdir()
    sc = wc / "stage_b_break2.py"
    sc.write_text(src_text.replace(GUARD, BREAK2 + GUARD), encoding="utf-8")
    results["C_break2_dotget"] = run_case("C. BREAK 2 (d.get on PromoCalcResults)", wc, sc,
                                          "FAIL: every tile AttributeError")

print("\n================ VERDICT ================")
a, b, c = results["A_current"], results["B_break1_import"], results["C_break2_dotget"]
checks = [
    ("current code decomposes tiles",        (a["ok"] or 0) > 0),
    ("current code has ZERO errors",         a["err"] == 0),
    ("current code gives ~73 sites/tile",    a["spt"] is not None and abs(a["spt"] - 73.0) < 1.0),
    ("break 1 fails EVERY tile",             (b["ok"] or 0) == 0 and (b["err"] or 0) > 0),
    ("break 1 fails as ImportError",         b["kind"] == "ImportError"),
    ("break 2 fails EVERY tile",             (c["ok"] or 0) == 0 and (c["err"] or 0) > 0),
    ("break 2 fails as AttributeError",      c["kind"] == "AttributeError"),
]
for name, passed in checks:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
allp = all(p for _, p in checks)
print(f"\nGATE IS {'NON-VACUOUS -- it passes clean code and catches BOTH breaks' if allp else 'NOT VERIFIED'}")
sys.exit(0 if allp else 1)
