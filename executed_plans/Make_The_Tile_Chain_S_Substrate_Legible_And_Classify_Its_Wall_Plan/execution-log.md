# Execution Log — Make_The_Tile_Chain_S_Substrate_Legible_And_Classify_Its_Wall_Plan
Date: 2026-10-21
Waves: 3 (max parallelism 2) — run SEQUENTIALLY; see the mode note below
Files changed: dna_decode/glm/genomewide.py, scripts/glm_tile_predictions.py, scripts/promcalc_toolchain_probe.py, tests/test_glm_tile_substrate.py, tests/test_promcalc_toolchain_probe.py, tests/test_tile_chain_substrate_wall.py, tests/test_tile_chain_constraint_registered.py, tests/test_tile_chain_next_md_record.py, wiki/tile_chain_substrate_wall_2026-10-21.md, wiki/tile_chain_substrate_wall_2026-10-21.json, wiki/promcalc_toolchain_probe_2026-10-21.json, soraya_runs/constraint-expiry.json, NEXT.md, wiki/plans-index.md, CHANGELOG.md, LESSONS_LEARNED.md, TODOS.md, FUTURE_FEATURES.md, README.md, .claude/testable-modules.md, wiki/decisions-log.md
Sentrux verdict: n/a — sentrux not installed on this host
Commit: 7df1f6b4 → 0e659e89 (five commits on `main`), archived in the commit that carries this log

## Mode

Parallel was ELIGIBLE (no toolkit errors, two 2-step waves, remote configured) and
**sequential was chosen** — the same override as the preceding plan execution, for the same
measured reason. `.gitignore:4` is `/data/*` with 59 files force-added, so a worktree checkout
gets those 59 and none of the ~4 GB of untracked data beneath them. Its absent-data failure set
therefore differs sharply from the parent's (the recorded 101 → 22 story), which makes per-step
baseline reconciliation inside a worktree meaningless. Recorded rather than inferred.

## Test reconciliation — every moved count named

| run | result |
|---|---|
| baseline (pre-change) | 14 failed / 6834 passed / 186 skipped |
| final (post Steps 1–5) | 14 failed / **6905** passed / 186 skipped |
| after `/test-epilogue` added 8 cases | 14 failed / **6913** passed / 186 skipped |

+71 from baseline to final = exactly the 71 tests this plan added. +8 more = exactly the cases
`/test-epilogue` added to `tests/test_glm_tile_substrate.py` (20 → 28). **Failure set identical at
every stage, zero regressions.** All 14 are absent-data: 7 cyp4f2/GeT-RM, 5 essentiality ladder,
1 pneumococcal cps DB, 1 README quickstart.

## Verification signals

1. **Clean** — no-arg run prints the `external_root` report and exits 1 with zero `WinError` or
   `Traceback` in its output (grep-confirmed 0).
2. **Clean, and stronger than specified** — the plan asked only that `--mpra-cache E:/…` get PAST
   the tile load. It printed `tiles 44106 | ceiling 0.9410`, both the published anchors reproduced
   exactly, so the `E:` copy is *the* published substrate rather than a file of the right name.
   This is the measurement the whole wall classification rests on.
3. **Clean** — a reachable `--in-dir` passes the new gate and falls through to the pre-existing
   "no decomposition shards" refusal, so the gate is not an always-refuse.
4. **Clean** — probe verdict `ENV_ABSENT_DEPS_SATISFIABLE_HERE`, 6/6 checks resolved, no problems.
5. **Clean** — memo names both dependencies with distinct classes; quotes none of the frozen bar's
   headline figures (enforced by test, with a non-vacuity partner asserting it DOES carry the two
   numbers this work measured).
6. **Clean** — 12 constraints, `file-exists` expiry, the original 11 field-identical by diff.
7. **Clean** — frozen five byte-unchanged; prospective lock `verify_lock ok = True` at 2026-08-31;
   arms ledger sha256 identical before and after.

## What did not happen, deliberately

* **The frozen acceptance bar was not edited** despite its `no uv` clause being false. A test
  asserts it still carries that clause, so the audit cannot quietly become an edit.
* **The arms ledger was not edited.** Per R6 a declared arm cannot evaporate; demoting a
  substrate-blocked arm would destroy the information that it is still owed.
* **`CLAUDE.md` was not edited** — 1 word of headroom against the enforced ceiling, which was not
  raised.
* **Two pre-existing ruff findings in `genomewide.py`** were confirmed present at HEAD and left
  alone rather than tag-along-fixed.
* **Nothing was installed.** The probe reports that an isolated environment is creatable; that is a
  different claim from having created one, and the install is a user decision.

## Incidents worth carrying

* A non-vacuity mutation **failed to apply** (multi-line anchor against a CRLF file) and the harness
  reported `MUTATION-DID-NOT-APPLY` rather than letting a correctly-green test read as vacuous —
  the trap this discipline exists for, firing as designed.
* A test of mine **asserted my model of a schema** rather than the schema: the arms ledger maps arm
  NAME → record, not `criteria[*]["arms"]` as a list. It failed with the arm plainly present. Fixed
  by reading the structure.
* `/test-epilogue` caught its own near-miss: a `write_text` restore silently converted a source file
  LF → CRLF while a `read_text` round-trip still compared equal. Found by byte comparison and
  restored exactly; verified here independently.
