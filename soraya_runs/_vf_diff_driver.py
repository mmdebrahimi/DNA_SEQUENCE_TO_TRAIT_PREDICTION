"""Canonical minister driver for ep4-v01-vf-diff (round 1: in_place gate park).

Drives the mission through scripts/run_minister.py per the v0.3.1 canonical entry point.
Round 1 of an in_place mission PARKS on the interrogation gate (no receipt yet) → blocked:user-only.
The model-driven seams are stubbed for round 1: in_place does not generate, and execution only fires
AFTER the interrogation receipt lands, so none of them are reachable this round.
"""
import os
import sys
from pathlib import Path

SORAYA = os.path.expanduser("~/.claude/skills/soraya/scripts")
sys.path.insert(0, SORAYA)
import run_minister as rm  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent          # dna_decode/
BIG_IDEA = PROJECT_ROOT / "project_state" / "ep4-v01-vf-diff" / "big-idea.md"
MISSION_ID = "ep4-v01-vf-diff"
RUN_ID = "2026-06-04-0008-ep4-v01-vf-diff"
RUN_DIR = PROJECT_ROOT / "soraya_runs" / RUN_ID


def _unreachable_seam(*a, **k):
    raise AssertionError("seam reached in round 1 of an in_place gate-park — unexpected")


def main():
    os.chdir(PROJECT_ROOT)   # ledger paths in big-idea.md resolve relative to cwd

    ok, problems, advisories = rm.preflight(
        money_hook_hardened=True,     # deployed real shim, matcher '*', no '|| true' (verified)
        agent_present=True,           # Agent tool present in this session
        pat_rotated=True,             # advisory only
        big_idea_path=str(BIG_IDEA),
        cwd=str(PROJECT_ROOT),
    )
    print("PREFLIGHT ok =", ok)
    if problems:
        print("PROBLEMS:", problems)
        sys.exit(2)
    for a in advisories:
        print("ADVISORY:", a)

    # Round 2: interrogation receipts are recorded (gate satisfied) and the real --until-mvp work
    # (vf_runner + vf_diff + cli wiring + test + AC9 row) was done by the model-driven seam already,
    # so the injected runner reports the terminal outcome. The driver re-checks endpoints live.
    driver = rm.build_driver(
        str(BIG_IDEA), str(RUN_DIR),
        generate_candidate=_unreachable_seam,
        until_mvp_runner=lambda pid: "retire",
        project_init_runner=_unreachable_seam,
        mission_id=MISSION_ID,
        # ship-bound: the family ships a user-facing artifact (vf_diff in the CLI) → high-stakes.
        high_stakes=lambda fid: True,
    )

    verdict, rounds, report = rm.run_to_stop(
        driver, arm_lease=True, run_id=RUN_ID, mission_id=MISSION_ID,
    )
    print("VERDICT =", verdict)
    print("ROUNDS  =", rounds)
    print("PHI     =", report["mission_potential"])
    print("FRONTIER=", report["frontier"])
    print("PARKED  =", report["parked"])


if __name__ == "__main__":
    main()
