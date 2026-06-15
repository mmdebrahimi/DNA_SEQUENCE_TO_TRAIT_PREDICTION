"""One-command Oxford external re-validation driver.

Chains the whole arm with run-scoping + gating as the PRIMARY invariant:
  W0 probe (advisory) -> build_oxford_labels -> external_cohort_preflight
  --cohort-manifest (abort != PASS unless --allow-degraded) -> external_cohort_revalidate
  per drug (propagate exit 3 hard-fail / 1 degraded) -> build_external_validation_report
  ONLY IF every required drug run is acceptable.

Pure helpers (worst_exit / roll_up_allowed / preflight_blocks / plan_steps) are split
from the subprocess execution so the control logic is unit-testable offline.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Driver exit severity: hard-fail (3) > unaccepted-degraded (1) > gate-refusal (2) > clean (0).
_EXIT_RANK = {3: 3, 1: 2, 2: 1, 0: 0}


def worst_exit(codes) -> int:
    """Return the worst child exit per the 3 > 1 > 2 > 0 severity order."""
    codes = list(codes) or [0]
    return max(codes, key=lambda c: _EXIT_RANK.get(c, 1))


def roll_up_allowed(drug_exits, allow_degraded: bool) -> bool:
    """Roll up ONLY if every drug run is acceptable: exit 0, or exit 1 when degraded
    is explicitly allowed. Any hard-fail (3) or gate-refusal (2) blocks publication."""
    for c in drug_exits:
        if c == 0:
            continue
        if c == 1 and allow_degraded:
            continue
        return False
    return bool(list(drug_exits))   # nothing ran -> nothing to publish


def preflight_blocks(verdict: str, allow_degraded: bool) -> bool:
    """Whether the preflight verdict aborts the run."""
    return verdict != "PASS" and not allow_degraded


def plan_steps(drugs) -> list[str]:
    """Ordered step labels (preview + test): probe, labels, preflight, per-drug, rollup."""
    return ["w0_probe", "build_labels", "preflight",
            *[f"revalidate:{d}" for d in drugs], "rollup"]


def _run(argv) -> int:
    print(f"  $ {' '.join(str(x) for x in argv)}")
    return subprocess.run([sys.executable, *map(str, argv)]).returncode


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", required=True)
    ap.add_argument("--mic-table", required=True)
    ap.add_argument("--key-col", required=True)
    ap.add_argument("--drug-col", action="append", default=[], metavar="COL=ALIAS")
    ap.add_argument("--drugs", nargs="+", required=True, help="canonical drugs to score")
    ap.add_argument("--run-id", default=f"oxford_{_date.today().isoformat()}")
    ap.add_argument("--allow-degraded", action="store_true")
    ap.add_argument("--skip-smoke", action="store_true")
    a = ap.parse_args()

    run_id = a.run_id
    manifest = f"wiki/cohort_manifest_external_{run_id}.json"
    deg = ["--allow-degraded"] if a.allow_degraded else []

    print("[1/5] W0 probe (advisory)")
    _run(["scripts/oxford_w0_probe.py", "--project", a.project, "--mic-table", a.mic_table,
          "--key-col", a.key_col, *sum([["--drug-col", c] for c in a.drug_col], [])])

    print("[2/5] build labels + cohort manifest")
    rc = _run(["scripts/build_oxford_labels.py", "--project", a.project, "--mic-table", a.mic_table,
               "--key-col", a.key_col, "--run-id", run_id,
               *sum([["--drug-col", c] for c in a.drug_col], [])])
    if rc != 0:
        print("build_oxford_labels failed (e.g. crosswalk conflict) — aborting"); return rc

    print("[3/5] exact-set preflight")
    pf_rc = _run(["scripts/external_cohort_preflight.py", "--project", a.project,
                  "--cohort-name", "oxford", "--cohort-manifest", manifest, "--mic-open", *deg])
    if pf_rc != 0 and not a.allow_degraded:
        print("preflight did not PASS — aborting (pass --allow-degraded to override)"); return pf_rc

    print("[4/5] score per drug")
    drug_exits = []
    for drug in a.drugs:
        rc = _run(["scripts/external_cohort_revalidate.py", "--cohort", "oxford", "--drug", drug,
                   "--labels-dir", f"data/raw/oxford_extval_{drug}", "--cohort-manifest", manifest,
                   "--run-id", run_id, "--preflight-json",
                   f"wiki/external_preflight_oxford_{_date.today().isoformat()}.json",
                   *(["--skip-smoke"] if a.skip_smoke else []), *deg])
        drug_exits.append(rc)

    print("[5/5] roll-up gate")
    if roll_up_allowed(drug_exits, a.allow_degraded):
        _run(["scripts/build_external_validation_report.py", "--run-id", run_id, *deg])
    else:
        print(f"roll-up BLOCKED — drug exits {drug_exits} not all acceptable "
              f"(hard-fail/gate-refusal, or degraded without --allow-degraded)")

    return worst_exit([pf_rc, *drug_exits])


if __name__ == "__main__":
    raise SystemExit(main())
