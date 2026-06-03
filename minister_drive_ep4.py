"""Soraya minister maiden real-run driver — mission ep4-v01-expec-recall.

Resumable orchestration wrapper around scripts/run_minister.py. Reusing the SAME run_id
across sessions reuses the on-disk lifecycle-journal + gate-receipts, so:
  - round-1 GENERATE is idempotent (journal one-active-proposal guard),
  - interrogation receipts recorded between sessions clear the PARK,
  - promotion + --until-mvp run on a later invocation.

Usage:
  SORAYA_RUN_ID=<rid> python minister_drive_ep4.py          # drive to next stop
This is a run artifact (not project source); safe to delete after the mission closes.
"""
import os
import sys
from datetime import datetime
from pathlib import Path

SKILL = Path(os.path.expanduser("~/.claude/skills/soraya/scripts"))
sys.path.insert(0, str(SKILL))
REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))

import run_minister as rm          # noqa: E402
import run_dir as rd              # noqa: E402
import gate_receipts as gr        # noqa: E402
import mvp_criteria as mvp        # noqa: E402

BASE = REPO
MISSION_ID = "ep4-v01-expec-recall"
BIG_IDEA = REPO / "project_state/ep4-v01-expec-recall/big-idea.md"
GAP_ID = "improve-expec-recall"
FAMILY_ID = "fam-per-gene-expec-scoring"
EXPEC_TEST = "tests/test_pathotype_expec_recall.py"

FAMILY_GOAL = (
    "GOAL: Lift v0 pathotype resolver ExPEC recall from 0.75 to >=0.85 on the 24-genome "
    "H4 cohort WITHOUT regressing confident-supported-call precision (1.0) — via per-gene "
    "ExPEC marker scoring (split SIDEROPHORES/CAPSULE_SERUM clusters into per-gene presence "
    "and re-tune the resolver against the cached coverage in data/pathotype_cov_cache/). "
    "Retires gap improve-expec-recall: " + EXPEC_TEST + " passes + AC9 ledger row recorded."
)


def _generate_candidate(gap_id):
    return FAMILY_ID if gap_id == GAP_ID else None


def _until_mvp_runner(family_id):
    # Endpoints are re-checked by the driver (mission_met). The actual engineering work
    # (per-gene scoring + recall test + ledger row) is done by the model BEFORE the round
    # that executes the frontier. Here we report retire iff both gap endpoints are live-MET.
    import big_idea as bi
    big = bi.validate(bi.parse_big_idea(BIG_IDEA.read_text(encoding="utf-8")))
    gap = next(g for g in big.gaps if g.id == GAP_ID)
    deps = dict(exists=os.path.exists, read_ledger=rm._default_read_ledger,
                gated_run=rm._default_gated_run())
    return "retire" if all(mvp.evaluate(c, **deps) for c in gap.endpoint_criteria) else "no-path"


class _ProjectInitRunner:
    """F6 effect runner: materialize project_state/<family>.md (mirrors /project-init).
    expected_content drives the digest guard so a replay is a no-op (not a collision)."""
    def expected_content(self, family_id):
        return _ledger_text(family_id)

    def __call__(self, family_id, ledger_path, expected_content):
        p = REPO / ledger_path
        if not p.exists():
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(expected_content, encoding="utf-8")


def _ledger_text(family_id):
    return FAMILY_GOAL


def main():
    rid = os.environ.get("SORAYA_RUN_ID") or rd.make_run_id(MISSION_ID, now=datetime.now(), base=BASE)
    run_dir = BASE / rid
    fresh = not run_dir.exists()
    if fresh:
        rd.open_run(BASE, rid)
    print(f"[minister] run_id={rid} fresh={fresh}", flush=True)

    init_runner = _ProjectInitRunner()
    with rd.run_lock(BASE, rid):
        driver = rm.build_driver(
            BIG_IDEA, run_dir,
            generate_candidate=_generate_candidate,
            until_mvp_runner=_until_mvp_runner,
            project_init_runner=init_runner,
            mission_id=MISSION_ID,
            exists=os.path.exists,
            read_ledger=rm._default_read_ledger,
            high_stakes=lambda fid: False,   # local code + test + ledger row; not ship/global/migration/auth
        )
        verdict, rounds, report = rm.run_to_stop(
            driver, arm_lease=True, run_id=rid, mission_id=MISSION_ID)

    print(f"[minister] verdict={verdict} rounds={rounds} Phi={report['mission_potential']}", flush=True)
    print(f"[minister] frontier={report['frontier']} parked={report['parked']}", flush=True)
    # what each parked attempt still needs
    rdir = run_dir / "gate-receipts"
    for fam in {FAMILY_ID}:
        missing = gr.missing_for_promotion(rdir, fam)
        have = gr._count(rdir, fam, gr.INTERROGATION)
        print(f"[minister] family={fam} interrogation_receipts={have} missing={missing}", flush=True)
    print(f"[minister] RUN_ID_FOR_RESUME={rid}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
