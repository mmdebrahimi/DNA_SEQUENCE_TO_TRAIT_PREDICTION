# Return Decision Tree — 2026-05-16 (post-departure setup)

Quick reference when you return: which artifacts to check, in what order, and what to do with each combined outcome. Written 2026-05-16 14:50 before user departed.

## What was launched

| job | mechanism | expected output | ETA from launch (14:43) |
|---|---|---|---|
| Stage 1b (mean+max + scaled NT-logreg) | detached via `run_stage1b_detached.bat` (PID 4736 at launch) | `wiki/stage1_n40_cipro_mean-plus-max_<date>.md` + `stage1b_detached.log` | ~5 hr (verdict packet only written at end) |
| Preflight v2 (cipro attribution) | Bash bg job `bx13oua9k` | `wiki/cipro_attribution_preflight_<date>.md` + `.per_strain.json` sidecar | ~5-10 min |
| Heartbeat sidecar (Stage 1b monitor) | `heartbeat_stage1b.bat 4736` (broken sleep — loops fast, ~2 lines/sec; harmless but noisy) | `stage1b_heartbeat.log` (timestamped ALIVE/DEAD lines) | self-terminates when Stage 1b dies OR verdict lands |

## Step 1 — Triage what survived

```bash
# Check Stage 1b liveness
tasklist | grep python.exe         # is PID 4736 (or any heavy python.exe) still alive?
ls -la wiki/stage1_n40_cipro_mean-plus-max_*.md   # verdict packet present?
tail -20 stage1b_heartbeat.log     # last line: ALIVE / DEAD / VERDICT_LANDED

# Check preflight v2
ls -la wiki/cipro_attribution_preflight_*.md      # verdict packet present?
ls -la wiki/cipro_attribution_preflight_*.per_strain.json
```

## Step 2 — Decision tree by combined state

### A. Both verdict packets exist
- Read preflight verdict first (cipro_attribution_preflight_*.md):
  - **STRONG_POSITIVE** = NT has biology; Stage 1's FAIL was likely a head/pooling/LOSO-noise issue
  - **WEAK_POSITIVE** = suggestive but not load-bearing
  - **INCONCLUSIVE_MISS** = NOT damning; mean-pool dilution + symbol coverage may be confounders
  - **DAMNING_MISS** = requires v2 + mean+max preflight both miss (not run yet)
- Then read Stage 1b verdict packet:
  - **`stage2_action == BURST_STAGE_2`** → CLEAN PASS; Stage 2 N=150 cipro is unblocked
  - **`stage2_action == HOLD_STAGE_2_CI_DEGENERATE`** → NOISY PASS; hold burst budget
  - **`stage2_action == ALTERNATIVE_POOLING_RERUN`** → still FAIL under mean+max + scaling fix; escalate to PIVOT_TO_BAKTA
  - **`stage2_action == PIVOT_TO_BAKTA`** → architecture mismatch; Bakta re-annotation + gene-presence comparator pathway
- Combined call:
  - `BURST_STAGE_2` + STRONG/WEAK positive → fire Stage 2 Databricks setup (resolve Pending Decisions row 6); EP2 cef/tet smoke as parallel falsification probe
  - `BURST_STAGE_2` + INCONCLUSIVE_MISS → fire Stage 2 BUT note the attribution gap; require gyrA/parC/parE attribution at Stage 2 verdict
  - `HOLD_STAGE_2_CI_DEGENERATE` → AMRFinderPlus POINT* parser as classical comparator tightener (Pending Decisions row 7) before re-running
  - `ALTERNATIVE_POOLING_RERUN` outcome → PIVOT_TO_BAKTA: Bakta re-annotate N=38 cohort, rebuild gene-presence variant, re-run Stage 1c
  - `PIVOT_TO_BAKTA` outcome → end-of-line for NT-XGBoost cipro architecture; demote NT track or pivot to PATRIC-style curated-feature head

### B. Stage 1b packet missing, preflight verdict exists
- Check heartbeat: last ALIVE timestamp vs current time
  - heartbeat says ALIVE within last ~5 min → Stage 1b still running; wait
  - heartbeat says DEAD or stale > 30 min → Stage 1b crashed silently
    - inspect `stage1b_detached.log` for partial output (may need to flush by killing process explicitly with `taskkill /F /PID <pid>`)
    - relaunch with `python -u` + heartbeat wrapper (the prior wrapper is broken; rewrite with PowerShell `Start-Sleep`)
- Apply preflight verdict for partial direction:
  - STRONG_POSITIVE → architecture has signal; Stage 1b worth retrying
  - INCONCLUSIVE_MISS → ambiguous; await Stage 1b retry

### C. Preflight packet missing, Stage 1b packet exists
- Preflight likely crashed; check `C:/Users/Farshad/AppData/Local/Temp/preflight_v2.log` for traceback
- Re-run: `uv run python scripts/cipro_attribution_preflight.py 2>&1 | tee preflight_v2_retry.log`
- Apply Stage 1b verdict alone:
  - `BURST_STAGE_2` → Stage 2 burst is gate-supported; preflight is informational not load-bearing
  - Anything else → escalate per Section A's Stage 1b table

### D. Both packets missing
- This is the worst-case Codex flagged. Inspect:
  - `tasklist | grep python.exe` — anything alive?
  - `stage1b_heartbeat.log` last line — when did Stage 1b stop?
  - `stage1b_detached.log` + `C:/Users/Farshad/AppData/Local/Temp/preflight_v2.log` — partial output? traceback?
- Recovery: relaunch BOTH with observability:
  - Stage 1b: `set PYTHONUNBUFFERED=1` + .bat wrapper + heartbeat with PowerShell sleep
  - Preflight: foreground, NOT bg, NOT piped through `| tail`

## Step 3 — Cleanup / ledger updates regardless of outcome

1. `git status --short` — capture untracked artifacts produced by the runs (heartbeat log, smoke output, etc.)
2. `/project-state --append-action --class run-tests` with the Stage 1b outcome row (e.g., `Stage 1b verdict packet at wiki/stage1_n40_cipro_mean-plus-max_<date>.md: <verdict>; stage2_action <action>; aggregation mean+max; LR scaled`).
3. `/project-state --update-hypothesis H17 --status <falsified|confirmed|under-investigation> --note "<H17 evidence per EP2 status>"` ONLY if EP2 cef + tet smoke was run (it wasn't pre-departure; defer until EP2 fires).
4. `/project-state --refresh-frame --current-state "<post-Stage-1b state>"` once verdicts land.
5. Ledger commit batch (per session pattern).
6. Push.

## Step 4 — If user wants to fire EP2 cef + tet smoke next

Pre-requisite per `plans/EP2_Cef_Tet_Smoke_Design_Plan.md` Step 2: rename + generalize `scripts/smoke_gate_12strain_cipro.py` → `scripts/smoke_gate_12strain.py` (output paths cipro-hardcoded at lines 264, 272, 333-341). Cef + tet mini-cohorts already built at `data/processed/gate_b_mini_{cef,tet}_cohort.parquet`. EP2 smoke per drug uses existing NT cache only if cef/tet strains overlap cipro N=12 mini cohort cache — likely they DON'T fully overlap, so cef/tet NT populates would need to fire first (~2 hr each on GTX 860M).

## Quick reference — branch state at departure

- main: 31 commits ahead of origin pre-push, **0 ahead post-push** (push completed 14:50 just before departure)
- Last commit: `37af5f6 chore: detached-launch helper for long-running Stage 1 jobs`
- Untracked: `heartbeat_stage1b.bat` (new, runaway), `project_state/dna-decode-2026-05-11-scratch.md` (durable scratch), `wiki/smoke_gate_12strain_cipro_2026-05-15.md` (provenance unclear)

## Open questions surfacing on return

1. Did `start /B` detachment actually survive past 10 min? (heartbeat will tell — if log shows ALIVE lines after 15:00 EDT, yes)
2. Stage 1b verdict packet exit-code semantics — was the bug fixed or is it still "exit 0 ≠ BURST_STAGE_2" (filed in TODOS)
3. Whether the cef + tet mini-cohorts have downloadable assemblies (the prior cipro N=40 had a 38/40 GFF3-availability bottleneck; cef/tet may have different gaps)
