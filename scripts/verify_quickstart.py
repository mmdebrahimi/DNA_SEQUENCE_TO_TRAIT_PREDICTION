"""Verified quickstart -- runs the documented PURE-PYTHON decoder paths end-to-end and asserts they work.

The point of the productization pass: a default `uv sync` (no [ml] extra, no Docker, no external DBs)
must still deliver a working deterministic decoder. This script IS that proof -- it executes the exact
commands QUICKSTART.md documents and checks each one's call + inline trust badge. It needs NO
torch/transformers, NO Docker, NO network: only the committed fixtures + the pure-Python observed-
substitution paths + the offline-degradation path. Exit 0 = every quickstart step works.

SCOPE (honest): this verifies the SOURCE-TREE / EDITABLE install by calling each CLI's main(argv)
in-process. The complementary ARTIFACT-boundary check (build the wheel -> fresh-env install -> badges
served from the PACKAGED cards, no repo wiki/) is scripts/verify_wheel_ships_cards.py. The packaging gate
LANDED 2026-06-24: a built wheel now ships the 4 trust cards as package data (dna_decode/report_cards/ via
the pyproject force-include), so a wheel install no longer silently degrades the validation badges.

Run: `uv run python scripts/verify_quickstart.py`  (also pinned by tests/test_verify_quickstart.py).
"""
from __future__ import annotations

import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.amr.cli import main as amr_main          # noqa: E402
from dna_decode.profile.cli import main as profile_main  # noqa: E402

_FASTA = str(REPO / "tests" / "fixtures" / "ecoli_mini" / "genome.fna")
_AMR_RUN = str(REPO / "tests" / "fixtures" / "amr_mini")

# (label, main, argv, expected stdout substrings, expected exit code)
STEPS = [
    ("HIV observed -> R + free-wetlab badge (pure-Python)",
     amr_main, ["--drug", "efavirenz", "--observed", "RT:K103N", "--sample-id", "q"],
     ["CALL: R", "INDEPENDENT_WETLAB"], 0),
    ("fungal observed -> R + no-free-source badge (pure-Python)",
     amr_main, ["--drug", "fluconazole", "--observed", "ERG11:Y132F", "--sample-id", "q"],
     ["CALL: R", "NO_FREE_PHENOTYPE_SOURCE"], 0),
    ("SARS-CoV-2 observed -> R + in-distribution badge (pure-Python)",
     amr_main, ["--drug", "nirmatrelvir", "--observed", "Mpro:E166V", "--sample-id", "q"],
     ["CALL: R", "IN_DISTRIBUTION"], 0),
    ("unified profile (cached AMR + offline typing degradation)",
     profile_main, [_FASTA, "--amrfinder-run", _AMR_RUN],
     ["amr (Escherichia)", "INDEPENDENT_MEASURED", "ASSUMED"], 0),
    ("unified profile offline (AMR section degrades, exit 0)",
     profile_main, [_FASTA],
     ["amr: [unavailable]"], 0),
]


def run_step(main, argv) -> tuple[int, str]:
    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            rc = main(argv)
    except SystemExit as e:  # a CLI that calls sys.exit instead of returning
        rc = int(e.code or 0)
    return rc, buf.getvalue()


def main(argv=None) -> int:
    failures = []
    for label, fn, args, expect, exp_rc in STEPS:
        rc, out = run_step(fn, args)
        missing = [s for s in expect if s not in out]
        ok = (rc == exp_rc) and not missing
        print(f"[{'PASS' if ok else 'FAIL'}] {label}  (rc={rc})")
        if not ok:
            if rc != exp_rc:
                print(f"        expected rc={exp_rc}, got {rc}")
            for s in missing:
                print(f"        missing in output: {s!r}")
            failures.append(label)
    n = len(STEPS)
    print(f"\nquickstart: {n - len(failures)}/{n} steps passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
