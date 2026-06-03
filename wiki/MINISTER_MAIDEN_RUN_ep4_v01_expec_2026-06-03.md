# Minister maiden real run — EP-4 v0.1 ExPEC recall — launch handoff (2026-06-03)

> Run THIS in a dna_decode session (cwd = C:\Users\Farshad\PythonProjects\dna_decode). This is the FIRST
> real `/soraya minister` run driven to mission-mvp-reached (prior 4 dogfoods parked at blocked:user-only
> or were disposable). Mission is BOUNDED (per the open-world-boundary lesson). Soraya is v0.3.2 (291+ tests);
> the money hook is LIVE-deployed + the lease now auto-arms via run_minister.py.

## #1 NEXT ACTION — launch from a session ROOTED in dna_decode (the cwd is decisive)
The first attempt (2026-06-03) preflight-PASSED but Soraya HELD because the session cwd was
`rca_engine/articles`, not `dna_decode` → `/project-state` path-refuses → the minister would park on its own
AC9 ledger writes. (Soraya v0.3.3 now HARD-refuses this at preflight instead of parking mid-run.) So:
```
# exit the current session, then:
cd C:\Users\Farshad\PythonProjects\dna_decode
claude              # FRESH session rooted here (don't resume a session whose cwd is elsewhere)
/soraya minister ep4-v01-expec-recall
```
Also: let the background NT bake-off finish first (GPU contention) before launching.

## The mission (already drafted + validated)
`project_state/ep4-v01-expec-recall/big-idea.md` — ONE gap `improve-expec-recall` (strategy-budget 2),
endpoints:
1. `test-exit-0 python -m pytest tests/test_pathotype_expec_recall.py` — a NEW test asserting ExPEC recall
   ≥ 0.85 on the 24-genome H4 cohort AND confident-supported-call precision still 1.0 (no regression).
2. `project-state-row project_state/ecoli-pathotype-prediction-cli-2026-05-26.md:ExPEC recall hardened` —
   an AC9 ledger decision row recording the achieved metric (provenance).
The gap retires only when BOTH are live-MET.

## What the run will do (expect)
1. preflight (money hook hardened ✓, Agent present ✓; PAT advisory) + bounded-mission check.
2. Freeze the gap, Φ = 2.
3. Ascend: GENERATE one candidate family for the gap (e.g. fam-per-gene-expec-scoring: split the
   SIDEROPHORES/CAPSULE marker clusters into per-gene presence in `dna_decode/pathotype/{markers,resolve}.py`,
   re-tune against the cached coverage in `data/pathotype_cov_cache/` — instant reruns, no re-detect).
4. EMIT `/interrogate-me` ×2 for THAT family (you run them — Soraya can't self-invoke). Not high-stakes
   (local code + a test + a ledger row; no shipped skill/global-state/migration/auth) → only the 2
   interrogation receipts are needed, NOT /technical-plan + /brainstorm.
5. On receipts → self-invoke `/project-init` to seed the family ledger → atomic promotion → family joins
   the frontier → real `--until-mvp` drives it: implement per-gene scoring + write the recall test until
   both endpoints live-MET → gap retires → **mission-mvp-reached**.

## Caveats / watch-fors (from the dogfood corpus)
- **test-runner: VERIFIED CLEAR (2026-06-03).** `python -m pytest tests/test_pathotype_resolver.py` from
  the dna_decode repo root with the SYSTEM python (Python311) resolves pytest 9.0.3 and passes 20/20 — the
  pure-Python pathotype package imports from cwd, NO venv activation needed. So the gap endpoint
  `test-exit-0 python -m pytest tests/test_pathotype_expec_recall.py` evaluates correctly through the gated
  runner (won't fail-closed on a missing runner). Just launch the minister from the repo root.
- **budget 2** gives one retry if the first candidate is rejected at interrogation (per the budget≥2 rule).
- **No external deps / no money** — the whole mission is pure-Python CPU-only against the coverage cache,
  so the money/destructive gates should never fire (clean path to mission-mvp-reached).
- **Don't drift to the apex** — if the model proposes broadening beyond ExPEC recall, that's scope creep;
  the mission is frozen to this one gap.
- **dna_decode is a git repo** → per-step commits + rollback available (unlike the Skill_Development repo).

## After the run
Capture observations to `dogfood-observations/minister-real-run-ep4-v01-<date>.md` (objective per-round:
Φ, verdict, emitted commands, promote vs park, --until-mvp outcome, did it reach mission-mvp-reached).
This is the validation that closes Soraya's "#3 — one complete real bounded minister run" item.
