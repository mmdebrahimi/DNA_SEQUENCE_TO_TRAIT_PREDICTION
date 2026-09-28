# Result — advance: prospective scoring (Campylobacter x ciprofloxacin)

terminal_reason: **blocked:external — host memory pressure**. The scoring job was STOPPED BY THE
HARNESS because the system ran critically low on memory (measured after: 1.5 GB free of 15.9 GB,
91% used). This was not a command failure and says nothing about the command. Per the harness
instruction it was NOT restarted.

## What landed

- Accrual sweep re-run: 936,325 PD rows -> 14,320 with measured AST + assembly -> **442 eligible
  post-lock isolates** (2026-09-03..09-21), `overall_status=OK`. SEVEN TIMES the first accrual's 63.
- Per-cell powering derived: exactly ONE cell clears >=10-per-class --
  **Campylobacter x ciprofloxacin 149R/277S**. Klebsiella x meropenem (7R/9S) is UNDERPOWERED.
- `verify_lock OK` against the 2026-08-31 v2 lock, 5 files pinned, before any scoring.
- Scoring ran to **184 of 426 isolates (72R / 112S)**, all preserved on disk in
  `data/amrfinder_runs/` + `data/refseq_cache/` (the junction to D:). A restart resumes at 242
  remaining, not from scratch -- `ensure_run` is per-GCA cache-aware.
- Spillover (committed `7edc6c5`): property test closing the hidden-number extraction class, with a
  non-vacuity control that reconstructs the old guard and requires failure. It caught a dead branch
  in its own attribution function. tests 72 -> 75.

## NO prospective number was produced

The scorer writes its artifact at the END of the loop, so no
`prospective_lock_validation_Campylobacter_ciprofloxacin_*.json` exists. sens/spec are UNKNOWN and
are not estimated here.

## The one thing worth knowing for next time

The 184 completed isolates are **already powered** (72R/112S, vs a >=10-per-class bar and a frozen
cell of n=40). A CACHE-ONLY scoring pass over just those -- no downloads, no AMRFinder -- would
yield a real, honestly-labelled partial prospective number. It is cheap. It was NOT run, because
the harness said not to restart without being asked and memory is still short.

If taken, it must be labelled a PARTIAL cohort (184 of 426, selection = whatever the sweep ordering
reached first), which is a weaker claim than the full 426 and must not be presented as the latter.
