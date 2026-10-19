# Execution Log — Tile_Substrate_Selection_Value_And_The_Untried_Hybrid_Arm_Plan
Date: 2026-10-14
Waves: 3 (max parallelism 2) — executed SEQUENTIALLY by deliberate override; see below
Files changed: scripts/glm_tile_predictions.py, scripts/glm_tile_selection.py, scripts/glm_tile_hybrid_arm.py, tests/test_glm_tile_arms_declared.py, tests/test_glm_tile_predictions.py, tests/test_glm_tile_selection.py, tests/test_glm_tile_hybrid_arm.py, soraya_runs/arms-ledger.json, wiki/glm_tile_selection_2026-10-14.{md,json}, wiki/glm_tile_hybrid_arm_2026-10-14.{md,json}, NEXT.md, LESSONS_LEARNED.md, README.md, TODOS.md, wiki/decisions-log.md, wiki/plans-index.md, .claude/testable-modules.md
Sentrux verdict: n/a — binary resolves (v0.5.7) but no `.sentrux/rules.toml` exists in this repo, so the gate has nothing to enforce
Commit: cff96c3 (chain: c12c1a1, 4464c68, f552fc0, 2b091fb, 4e46670, cff96c3 — committed directly to `soraya/f3-calibration-and-cabbage-cell`, no PR)

## Mode: sequential, against the toolkit's mechanical selection

The toolkit returned `errors: []` with a 2-step wave and a configured remote + authenticated `gh`,
which selects parallel. Sequential was chosen anyway, for three reasons recorded before execution:

1. Parallelism buys ~nothing — the critical path is one long step (Step 2), and Steps 1/3/4/5 are
   seconds each.
2. Steps 1 and 5 write `soraya_runs/arms-ledger.json`, a cross-session shared file, while **3 other
   lane worktrees were live** (`lane/refseq-cache-integrity`, `lane/fba-structure`,
   `lane/nt-transformers-fix`).
3. Fresh `claude -p` worktree agents would have re-derived the subtle facts this session had pinned —
   the `{id(t): position}` join, the `log1p`-train/raw-test target asymmetry, `REGISTERED_PROTOCOL`,
   and the published reconcile anchors. Getting the join wrong yields a plausible number.

## Results

| | verdict | headline |
|---|---|---|
| Step 3 selection value | `SELECTION_IS_THE_BOTTLENECK` | best tile scorer realises **0.4807** median of achievable selection value at k=10; **2 of 10** seeds clear the frozen 0.50 bar |
| Step 4 hybrid arm | `HYBRID_ADDS` | paired **+0.0339** over the max_tx scalar, **10/10** seeds, CI95 **[0.0238, 0.0409]** excluding zero |

**The Step 3 result refutes the recommendation that motivated the plan.** The proposal was that the
scoring half was probably already adequate (citing the designed grid's 0.7937 at k=10) and that the
prediction criterion should likely be retired for generator work. It does not transfer: 0.485 here
against 0.794 there. The plan's own Risk Flag 6 named this outcome in advance, which is the only reason
the run could settle it rather than confirm the proposer.

Qualifications that must travel with it: it is a **near-miss**, not a collapse (0.4807 vs an ASSERTED
0.50 bar; the k=10 null band is stdev 0.0485, so 0.48 sits ~10 stdevs above chance; bar B passes
10/10 at 0.485 vs gc 0.295), and the **k-trend is monotone** — 0.388 (k=2) → 0.485 (k=10) → **0.537
(k=20)** — so the current oracle clears the bar at k=20 and the answer is conditional on how many
candidates the loop can afford to assay.

Step 4 is the **first learned component to measurably contribute on this substrate** and reproduces the
published GEMME pattern on DNA: conv ALONE loses to max_tx (−0.0767, 0/10) and still lifts it. Control
decisive in the right direction — max_tx with a *shuffled* conv scores −0.1232, 0/10. It does **not**
overturn the G-A encoder-gate negative (conv vs one GC number is still CI-includes-zero; a different
question).

## Verification

- Reconcile far exceeded its gate: the 2-anchor requirement was met and the table independently
  reproduced **all four** published medians to <3e-5 (gc 0.3112, conv 0.3234, max_tx 0.3968,
  decomposition 0.3925), which also confirms `fit_encoder` under `REGISTERED_PROTOCOL` reproduces the
  published encoder gate.
- Tests 6495 → 6566 passed, 14 skipped, **0 failed**; delta +71 exactly accounts for the new tests.
  Epilogue added 35 more (71 → 106 in the four files); 111 passed including the CLAUDE.md ceiling guard.
- All five sha256-pinned frozen files byte-unchanged (checked on the REAL paths —
  `dna_decode/data/amr_rules.py` does not exist and ad-hoc guards naming it have passed vacuously).
- Active 2026-08-31 prospective lock re-verified `ok=True, drifted=[], missing=[], incomplete_pin=[]`.
- Arms ledger: 10 arms, 6 measured / 4 planned; the hybrid is now the best measured arm (0.4254) and
  `ceiling_verdict` is still **`blocked:untried-arm`** — R5 working as designed.
- The three scripts verified byte-identical to their committed blobs by `git hash-object` after the
  epilogue's mutation testing (size alone cannot see a length-preserving change).

## Defects found during execution

1. **The plan's own staged workflow would have silently destroyed its intermediate results.** Each
   `--scorers` invocation rewrote the seed file from scratch, so the second run deleted the first's
   columns — and the damaged table still validated, round-tripped and reconciled, reading as "those
   scorers were never run" rather than as data loss. Fixed with merge-or-refuse. The reconcile gate
   could not see it: it checks that columns PRESENT are faithful, never that columns EXPECTED are
   present.
2. **A cost estimate in the plan was wrong by an order of magnitude** — the conv arm was sized as "the
   long pole, hours" and measured ~75 s/seed, ~13 min for ten.
3. Two mutation checks masked each other when applied together (epilogue finding), and one of the
   original tests was weaker than it read (ascending input made its correct answer sorted).

## Carried forward, deliberately not done

- **Re-run the selection measurement with the HYBRID as primary scorer.** It is now the best measured
  arm and Step 3's bar was frozen with `primary_scorer: max_tx` before any number was read. Asking
  whether the new best scorer clears the *same unchanged* bar is legitimate and cheap, but needs its
  own pre-registered run rather than a retrofit into a committed artifact.
- **Authority fork:** does a learned component contributing only IN COMBINATION satisfy the
  self-distillation ban's precondition (*"forbidden until a learned representation lands"*)?
- **Unresolved conflict, surfaced not overridden:** `wiki/decisions-log.md` 2026-10-13 found the
  fragment rung attractive at the corrected 0.79 ceiling; this session measured only 3.30% of its
  321,123 rows above 2× median. Different quantities; both can hold. No step here touched that arm.
- `CLAUDE.md` is at its enforced word ceiling, so these findings live in the two `wiki/` artifacts plus
  a `NEXT.md` entry. Whether to compress it is a user decision.
