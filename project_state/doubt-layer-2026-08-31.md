# F-A · Doubt Layer (L2)
<!-- project-schema: 0.1 -->

> Initialized 2026-08-31. Project ID: doubt-layer-2026-08-31. Originating goal (verbatim user input): "Build the L2 DOUBT layer: generalize and surface the catalog-completeness signals so every decoder call can carry a machine-readable 'this call may be incomplete, and here is why' block that never emits a competing resistance call".

## Project Context
- **Project ID:** doubt-layer-2026-08-31
- **Project root:** C:/Users/Farshad/PythonProjects/dna_decode
- **Captured:** 2026-08-31
- **Originating goal:** Build the L2 DOUBT layer: generalize and surface the catalog-completeness signals so every decoder call can carry a machine-readable "this call may be incomplete, and here is why" block that never emits a competing resistance call
- **Refined goal (if 3c produced one):** Ship a pure, cell-agnostic doubt layer whose signals are functions over already-computed determinant calls, wired into the decoder record behind a guard test that the block can never contain a resistance call, and measured per cell against the position-novelty flag's 0.604 incumbent.
- **Horizon (months):** 3
- **Schema:** project-schema 0.1
- **Family of:** dna-decode-2026-05-11 (umbrella) · plan `plans/Hybrid_Decoder_Architecture_Plan.md`

## Empirical Concerns
- **Verdict:** N-A
- **Check status:** not-applicable
- **Provisional:** NO
- **Findings:** (The goal text is imperative and contains no factual-shape claim about the external world. Its one empirical presupposition — that catalog-completeness signals exist and that completeness is the measured failure mode — is a claim about THIS repository's own committed artifacts, not about the literature, so a web check cannot adjudicate it. It is instead grounded directly in Evidence rows E1-E4 below, each citing a committed artifact. Recording N-A rather than the `--no-web-check` FLAG because the flag would assert unresolved factual uncertainty where the premise is in fact directly verified.)

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Horizon 3 months, under the 12-month bar. The deliverable is bounded and named: a module, a record field, a guard test, and a per-cell measurement artifact. No unbounded verbatim signal ("decode all", "understand any") is present. The success criterion is checkable by test exit status and file existence.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  - **C1 (selected)** — A doubt-layer module exists whose signals are pure functions over already-computed determinant calls (no model, no network, no structures), with tests. *Falsifier:* `uv run pytest tests/test_doubt_layer.py` exits non-zero, or the module imports torch / requests / a structure parser.
  - **C2 (selected)** — The decoder record carries a `doubt` block, and a guard test asserts the block can never contain a resistance call. *Falsifier:* a record is produced whose `doubt` block contains an R/S prediction and the suite stays green.
  - **C3 (selected)** — An artifact reports doubt-layer sensitivity PER CELL against the position-novelty flag's measured 0.604 on the EFV blind spot as incumbent. *Falsifier:* the artifact reports a pooled number, or reports no comparator.
  - **C4 (deferred)** — Doubt signals are visible on the trust surface without changing any cell's evidence tier. *Falsifier:* a cell's `tier` differs before vs after registration.
  - **C5 (rejected)** — The doubt layer improves resistance prediction accuracy. REJECTED at definition time: this is the failed-predictor framing (0-for-5, de-confounded) and adopting it would put L2 back in the regime that keeps dying. L2 must never compete with L1.

## Goal Hierarchy
### Long-term (12+ months tier)
Every decoder call in the tool carries an honest, machine-readable statement of where its own catalog is least trustworthy.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Completeness signal generalized from the two known gaps | a screen detects the shared shape from cached determinant calls and rediscovers `rmtE1` 36R/0S blind | DONE 2026-08-31 |
| 2 | Signal measured per cell against a named incumbent | artifact reports per-cell sensitivity vs the flag's 0.604; never pooled | 1 month |
| 3 | `doubt` block in the decoder record, guarded | guard test asserts the block can never carry a call | 1 month |
| 4 | Doubt visible on the trust surface, augment-only | no cell's evidence tier changes | 2 months |
| 5 | Doubt layer covers the target-site cells as well as AMR | both signal kinds route through one vocabulary | 3 months |

### Short-term (≤1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Full-index completeness run (prior run capped at 220 genomes/drug) | run-tests | Soraya | days |
| 2 | Per-cell doubt measurement artifact vs the 0.604 incumbent | edit-local-code | Soraya | days |
| 3 | `dna_decode/eval/doubt.py` — pure signals over computed calls | edit-local-code | Soraya | days |
| 4 | Guard test: the `doubt` block can never contain a call | run-tests | Soraya | days |
| 5 | Register augment-only on the trust surface | edit-local-code | Soraya | weeks |

## State Snapshot
### Assumptions
- The catalog's dominant failure mode is completeness, not accuracy — **high** (measured twice, independently).
- A doubt signal that never competes with L1 stays out of the failed-predictor regime — **high** (this is a design constraint, enforceable by test).
- Cheap deterministic signals are competitive with model-based ones here — **medium** (one measurement: the flag's 0.604 at zero tool cost).
- Per-cell reporting is required; pooled reporting would hide exactly the variance that matters — **high** (repo-wide precedent).

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| E1 | Gentamicin rule cannot represent `rmt*`; sens 0.523 vs 0.893 | wiki/gentamicin_rmt_disjoint_validation_2026-08-28.md | high | 2026-08-28 |
| E2 | HIV NNRTI catalog misses 53 resistant isolates at uncatalogued positions | wiki/hiv_esm_vs_catalog_2026-07-09.md | high | 2026-07-09 |
| E3 | Position-novelty flag recovers 60.4% of the EFV blind spot, lift 4.69, zero tools | wiki/hiv_blindspot_position_novelty_2026-07-11.json | high | 2026-07-11 |
| E4 | A general completeness screen rediscovers `rmtE1` 36R/0S blind, ranked first | wiki/determinant_completeness_screen_2026-08-31.md | high | 2026-08-31 |
| E5 | Learned predictors in the natural-population regime are 0-for-5 de-confounded | wiki/organism_gp_regime_correction_2026-08-29.md | high | 2026-08-29 |
<!-- project-state:end:evidence -->

### Unknowns
- Whether the completeness screen's `rmt_like` heuristic (>=3 R carriers, 0 S) generalizes beyond the one family it was tuned to recognise.
- What the doubt layer's false-positive rate is on cells with no known gap — unmeasured.
- Whether the AMR-side signal (determinant unrepresentable by the rule) and the target-site signal (novel substitution at a catalogued position) belong in one vocabulary or two.
- Whether a doubt block changes user behaviour at all, or is ignored like every other caveat.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | The two known gaps share one detectable shape | confirmed | 2026-08-31 |
| H2 | A deterministic doubt signal beats a model-based one on cost-adjusted value | under-investigation | 2026-08-31 |
| H3 | Doubt can be surfaced without changing any cell's evidence tier | open | (untested) |
| H4 | Cells with no known gap produce few doubt flags (low FP rate) | open | (untested) |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| L2 never emits a competing call | 2026-08-31 | the constraint that keeps it out of the 0-for-5 regime; enforced by guard test, not convention |
| Probe the deployed rule, never reimplement it | 2026-08-31 | a verbatim one-row probe cannot drift from DRUG_RULE |
| Rank by signature purity, not raw volume | 2026-08-31 | volume ranking buried the known gap 5th beneath correct exclusions |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Does L2 become a headline product claim? | Soraya | user authority | changes what the tool IS; currently an internal diagnostic |
| One doubt vocabulary or two (AMR vs target-site)? | Soraya | measurement | resolvable by executor once both are measured |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Step 1 shipped and validated (the screen rediscovers the known gap blind); steps 2-4 specified and carrying no authority.

### Target state / terminal condition
C1 + C2 + C3 all met: a pure doubt module with tests, a guarded `doubt` block in the record, and a per-cell measurement artifact naming the 0.604 incumbent.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` count + `gates-passed` count (raw counts, unweighted). 2026-08-31: 3 / 3 MVP criteria met (module+tests / guarded record block / per-cell artifact); 1 unknown retired (the `rmt_like` heuristic's false-positive rate is now measured: 4 of 5 raw hits are noise).
- **v0.2+:** weighted combination of unknowns-retired, gates-passed, evidence-confidence-improved, hypotheses-falsified (TBD via v0.2 design)

### Candidate next actions
Steps 1-5 all completed 2026-08-31 (see Action Log rows 2-7); the table below is the FOLLOW-ON set.
A stale candidate table is not cosmetic — `advance_ranker` reads the FIRST row as the family's next
action, so leaving completed work here made a finished family rank first on the portfolio frontier.

| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Measure the doubt layer's FALSE-POSITIVE rate on cells with no known gap | run-tests | med | high | high | 1-2hr |
| 2 | Extend the completeness screen to the target-site catalogs (one vocabulary or two?) | edit-local-code | med | high | med | 1-2hr |
| 3 | Re-screen when any NEW independent label set lands (both known gaps needed one) | research | low | high | high | ongoing |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** re-run `/project-state` after any action class fires (auto-append to Action Log triggers stale-state check)
- **Manual override:** user invokes `/project-state doubt-layer-2026-08-31` at any time
- **Family-specific:** re-evaluate whenever a NEW independent label set arrives — both known gaps were invisible until one did.

## MVP Criteria
- `test-exit-0 uv run pytest tests/test_doubt_layer.py -q`
- `test-exit-0 uv run pytest tests/test_doubt_record_guard.py -q`
- `file-exists wiki/doubt_layer_per_cell_2026-08-31.json`

Attempt budget: 3.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` — auto
- `research` — auto
- `write-plan` — auto
- `edit-local-code` — REQUIRES per-action human approval
- `run-tests` — auto if local + sandboxed
- `ask-user` — auto
- `stop` — auto

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-08-31 | propose | /project-init invoked | ledger created |
| 2 | 2026-08-31 | edit-local-code | scripts/determinant_completeness_screen.py + 9 tests (step 1) | shipped; rediscovers rmtE1 36R/0S blind, ranked first |
| 3 | 2026-08-31 | run-tests | step 2: full-index completeness run (1818 genomes, 6 drugs) | raw signature fires on 5 families, only 1 is real -- powering needed |
| 4 | 2026-08-31 | edit-local-code | step 2: dna_decode/eval/doubt.py + per-cell artifact | 1 of 1279 families STRONG (rmtE1 p=4.11e-12); enrichment null rejected as wrong |
| 5 | 2026-08-31 | edit-local-code | step 3: doubt block wired at the _target_site_record seam | verified on the real CLI; guard raises rather than emit a call |
| 6 | 2026-08-31 | edit-local-code | step 4: trust_surface.doubt_layer_for, augment-only | badge-with-vs-without diff green + non-vacuity pinned |
| 7 | 2026-08-31 | run-tests | full suite | 4029 passed, 0 failed; frozen surface byte-unchanged; lock re-verified |
<!-- project-state:end:action-log -->

## Open Questions for User
- Whether L2 becomes a headline product claim or stays an internal diagnostic. It changes what the tool is, and no executor step can settle it.
- Whether a doubt block should ever be able to request abstention from L1 (currently: no — it may qualify and explain, never overrule).

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** 1 of 3 MVP criteria addressed (step 1 shipped); criteria themselves not yet met.
