# Execution Log — Contract_Number_Provenance_Sound_Decoy_Control_Typed_Provenance_Kinds_Plan
Date: 2026-09-27
Waves: 5 (max parallelism 1)
Files changed: scripts/contract_number_audit.py, tests/test_contract_number_audit.py, dna_decode/data/cell_registry.py, wiki/contract_number_audit_2026-09-27.json, wiki/contract_number_provenance_triage_2026-09-27.md, CLAUDE.md, project_state/evidence-surface-2026-08-31.md, LESSONS_LEARNED.md, TODOS.md, wiki/decisions-log.md, wiki/plans-index.md, .claude/testable-modules.md
Sentrux verdict: n/a — sentrux not installed
Commit: 6f4b420 (last of 11 on main; steps: 4656f41, 92b1035, c938833, e90e8e7, bf8d100+b327ff2)

## As-built deviations from the plan

All 5 steps completed. Three deviations, each forced by a measurement taken during execution and each
recorded rather than absorbed:

1. **`HIGH` = rate exactly 0 was unreachable, so the band is now < 5%.** The plan carried "HIGH when the
   rate is 0" over from the 10-sample rule without re-deriving it for full-pool resolution. No cell
   reaches zero — not even `typing:Salmonella:salmserovar`, whose 30 figures still co-occur in 1 of 119
   haystacks. A grade no cell can earn is not a grade, and it collapsed MODERATE across 0.17%–41%.
2. **`essentiality:any:essentiality` grades MODERATE, not the LOW the plan predicted.** 44/600 = 7.3%
   uncited, 39/199 = 19.6% as now cited. The test asserts "not HIGH", which is the substantive claim.
3. **The citation rule admitted a cell the plan expected it to reject (`essentiality`) and rejected two
   it expected to admit** (`pneumoserotype`, and the single-file `hla_b5701` starting point, which is the
   wrong artifact — the prose spans three alleles). The rule was applied as fixed.

## Found during execution, NOT specified in the plan

- **The control was circular.** `main()` writes its dated artifact into the pool it scans, and that
  artifact lists the number tokens checked. Caught because `finder:any:forward` moved 74/300 → 75/300
  between two consecutive runs. Fixed (pool excludes `contract_number_audit_*.json`).
- **The artifact-side matcher accepts a substring inside a longer number** (`2.1` inside `72.1`) while the
  prose side is boundary-guarded. **4 of the 157 then-clean numbers rest on nothing else.** REPORTED, NOT
  FIXED — the fix flips the standing verdict to `ADJUDICATION_REQUIRED`, a user-authority call. It is also
  what caught two candidate citations that would otherwise have been added.
- **`DECIMAL_RE`'s `(?![\d.])` lookahead hides numbers before a sentence-final period** — 4 currently, a
  5th (`essentiality` 0.580) exposed incidentally by a Step 4 citation. REPORTED, NOT FIXED.

Both unfixed defects are tracked in `TODOS.md` and carried beside the headline verdict in `CLAUDE.md`.

## Outcome

- Standing audit: `NOTHING_TO_ADJUDICATE 0/176` across 19 cells (was `0/157` over 16).
- Unverifiable numbers: 34 → 9 in 3 cells. LOW-discrimination cells: 2 → 6.
- Tests: 4977 → 5014 passed, 11 skipped, 0 failed, 0 regressions.
- Frozen AMR surface byte-unchanged (`55c0920` / `c688ea7`); prospective lock re-verifies;
  `tier_evidence_audit` unchanged at under-claim 1 / over-claim 0 across 115 cells.
