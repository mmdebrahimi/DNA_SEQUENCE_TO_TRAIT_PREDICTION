# GLM G-A — Learned representation for the expression oracle
<!-- project-schema: 0.1 -->

> Initialized 2026-10-07. Project ID: glm-learned-representation-2026-10-07. Originating goal (verbatim user input): "Replace the expression oracle's fixed sequence features (GC / k-mer / positional one-hot fed to ridge) with a LEARNED multi-resolution encoder plus a trivial head, and measure whether it transfers to real genomic sequence."

## Project Context
- **Project ID:** glm-learned-representation-2026-10-07
- **Project root:** C:\Users\Farshad\PythonProjects\dna_decode
- **Captured:** 2026-10-07
- **Originating goal:** Replace the expression oracle's fixed sequence features with a learned multi-resolution encoder plus a trivial head, and measure whether it transfers to real genomic sequence.
- **Refined goal (3c):** Determine, with a committed artifact either way, whether a learned multi-resolution sequence encoder with a trivial head beats GC content (Spearman 0.3000) by >=0.05 on real *E. coli* genomic promoter tiles under leave-peak-out, reported as a fraction of the replicate-derived 0.9410 ceiling.
- **Horizon (months):** 12 (nominal; real expected duration ~1-2 weeks)
- **Schema:** project-schema 0.1
- **Portfolio:** GLM umbrella, family **G-A**. Namespace `G-*` chosen because `F-A`..`F-E` and `F1`..`F4` are taken by live work. See `plans/GLM_Portfolio_Decompose_2026-10-07.md`.

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** attempted
- **Provisional:** NO
- **Findings:** The goal carries factual-shape numeric claims, all INTERNAL measurements with a committed artifact, so the authoritative check is a read of that artifact rather than WebSearch (a web search cannot confirm a number produced in this repo today). Verified against `wiki/glm_genomewide_oracle_2026-10-07.json`: GC 0.3000 (fraction_of_ceiling 0.3188) · one-hot 0.2691 · 4-mer 0.2769 · 6-mer 0.2022 · tile ceiling 0.9410 (single-replicate 0.7945, Spearman-Brown reliability 0.8855) · n_kept 44,106 · verdict `DOES_NOT_GENERALISE` — ALL CONFIRMED. The designed-grid comparator 0.5914 is confirmed against `wiki/glm_expression_oracle_result_2026-10-07.md`. The distillation-degeneracy claim (student-teacher agreement 1.0000) is confirmed against `wiki/glm_selfdistill_probe_2026-10-07.json`. **One discrepancy found, in the ARTIFACT not in this goal:** that artifact's `honest_limits` states fragments are "48-475 bp with 280 distinct lengths", which is a 40,000-row SAMPLE figure, while its own `load_stats` records the full file as **48-499 bp with 359 distinct lengths**. It is a hardcoded prose string contradicting a derived field in the same file. It does not touch any claim in this goal (which concerns tiles, not fragments) but it is a real self-contradiction and is logged as a short-term action to fix by deriving the string from `load_stats`. **One EXTERNAL interpretive claim** — that AlphaGenome's transferable design principle is capacity-in-the-representation with a deliberately trivial head — rests on primary sources already read in full (the DeepMind technical talk `rdBtxtcS4nM` and the Kohli interview `h5zkzon0gM4`), recorded in `wiki/alphagenome_integration_analysis_2026-10-07.md`. It is an interpretation of those sources, not a quoted claim, and is flagged as such rather than asserted as the authors' wording.

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Bounded on every axis that matters. Named substrate (GSE144621 peak tiles, 44,106 at 150 bp, already downloaded), named split (`leave_peak_out`, already implemented), named baseline to beat (0.3000), named margin (>=0.05), named normaliser (0.9410, derived not asserted), and an explicit falsifier that closes the family rather than extending it. The nominal 12-month horizon is the skill default; the real expected duration is 1-2 weeks and nothing in the goal is open-ended.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  1. **(ADOPTED as the refined goal)** Determine with a committed artifact either way whether a learned multi-resolution encoder + trivial head beats GC's 0.3000 by >=0.05 on genomic peak tiles under leave-peak-out, reported against the 0.9410 ceiling.
  2. Narrower: hold the encoder architecture fixed at a single kernel width and test ONLY whether any learned encoder beats GC — isolating "learned vs fixed" from "multi-resolution vs single-resolution". Cheaper and answers the prior question first.
  3. Wider: train on the designed grid and TEST on genomic tiles (cross-substrate transfer), which answers a different and harder question — whether grid-learned representations transfer — rather than whether a genomic-trained encoder works.
  4. Diagnostic-first: before training anything, measure how much of the tile signal is linearly available at all by fitting GC + dinucleotide + positional-GC, to establish whether 0.3000 is near a composition ceiling. If it is, a learned encoder is predicted to fail and the family can close cheaply.
  5. Rejected: report the random-split number as the headline. Named only to exclude it — random splits inflate ~3x in this repo (measured four times) and tiles within a peak overlap heavily.

## Goal Hierarchy
### Long-term (12+ months tier)
An expression oracle whose representation is learned rather than hand-specified, validated on the same sequence distribution the generative loop will query it on.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Linear-headroom diagnostic | Composition-only ceiling on tiles measured; predicts whether a learned encoder can help | 1 week |
| 2 | Single-resolution encoder | A CNN trained under leave-peak-out, scored against GC 0.3000 and the 0.9410 ceiling | 2 weeks |
| 3 | Multi-resolution parameterisation | Kernel widths at motif/spacer/promoter scale compared against single-resolution, same split | 3 weeks |
| 4 | Verdict + artifact | Committed artifact stating PASS or the closed negative, with the falsifier honoured either way | 1 month |
| 5 | Unblock G-B or close it | Either a higher-capacity teacher exists for self-distillation, or G-B closes on a measured ceiling | 1 month |

### Short-term (<=1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Fix the 48-475/280 vs 48-499/359 self-contradiction by deriving the honest-limits string from `load_stats` | edit-local-code | Soraya | 1 day |
| 2 | Linear-headroom diagnostic on tiles (GC + dinucleotide + positional-GC, leave-peak-out) | run-tests | Soraya | 2 days |
| 3 | Implement `dna_decode/glm/encoder.py` — conv encoder + global pool + linear head, CPU-runnable | edit-local-code | Soraya | 3 days |
| 4 | Gate script `scripts/glm_encoder_gate.py` emitting a verdict artifact either way | edit-local-code | Soraya | 1 day |
| 5 | Run on tiles AND grid; commit the artifact; record the verdict in this ledger | run-tests | Soraya | 1 week |

## MVP Criteria

Attempt budget: 3 per criterion. **The MVP is an ANSWERED QUESTION, not a won bet** — a committed negative
satisfies this bar exactly as a positive does, because the family's deliverable is a verdict with a
falsifier honoured. A bar that only a success could clear would make the recorded negative unreachable.

| # | Criterion | Kind | Predicate |
|---|---|---|---|
| 1 | The encoder module exists | file-exists | `dna_decode/glm/encoder.py` |
| 2 | Its tests pass | test-exit-0 | `uv run pytest tests/test_glm_encoder.py -q` |
| 3 | The gate ran and produced a verdict artifact | file-exists | `wiki/glm_learned_representation_2026-10-07.json` |
| 4 | The verdict is recorded in this ledger | project-state-row | Action Log row whose outcome names the verdict |

## State Snapshot
### Assumptions
- A conv encoder has capacity a ridge-on-fixed-features model lacks — **high** confidence (architectural, not empirical).
- 44,106 tiles at 150 bp is enough data to train a small CNN without the overfitting that sank 6-mers (4,096 features reached only 0.2022) — **medium**.
- The 0.9410 ceiling is the right normaliser for the tile arm — **high** (derived from that subset's own replicates; the error of using the fragment ceiling was already caught and fixed).
- Promoter strength on real genomic tiles is predictable from sequence ALONE at better than composition level — **low**. This is the actual bet and it may be false.
- CPU is sufficient — **medium**; a GPU would shorten the loop but the data is small.

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | Fixed features fail off the designed grid: GC 0.3000 beats one-hot 0.2691 / 4-mer 0.2769 / 6-mer 0.2022 on genomic tiles | wiki/glm_genomewide_oracle_2026-10-07.json | high | 2026-10-07 |
| 2 | Tile noise ceiling 0.9410 (single-replicate 0.7945) — a CLEAN assay, so a weak model is a model finding | same artifact | high | 2026-10-07 |
| 3 | On the designed grid positional one-hot reaches 0.5914 and GC only 0.1414 — the ranking INVERTS between substrates | wiki/glm_expression_oracle_result_2026-10-07.md | high | 2026-10-07 |
| 4 | Self-distillation is degenerate on a linear teacher (student-teacher agreement 1.0000) so G-B is blocked on this family | wiki/glm_selfdistill_probe_2026-10-07.json | high | 2026-10-07 |
| 5 | GC's signal on tiles is NEGATIVE (-0.3253): AT-rich is more active, consistent with the -10 box TATAAT | wiki/glm_genomewide_oracle_result_2026-10-07.md | high | 2026-10-07 |
| 6 | 6-mers (4,096 features) score WORST on tiles, so capacity without inductive bias actively hurts here | wiki/glm_genomewide_oracle_2026-10-07.json | high | 2026-10-07 |
<!-- project-state:end:evidence -->

### Unknowns
- Whether 0.3000 is near a composition-only ceiling, which would predict failure before any training.
- Whether the tile signal is dominated by a few strong promoters (a tail problem) or spread across the set.
- Whether multi-resolution helps over single-resolution at 150 bp, where the whole window is one promoter.
- Whether a grid-trained encoder transfers to tiles (a different and harder question, candidate 3).
- Whether the 2/3-inactive composition of the tile set (15,171 active / 31,542 inactive) makes regression the wrong framing versus classification.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| GA-H1 | A learned encoder beats GC's 0.3000 by >=0.05 on genomic tiles under leave-peak-out | open | never |
| GA-H2 | Multi-resolution kernels beat a single resolution at 150 bp | open | never |
| GA-H3 | 0.3000 is near a composition-only ceiling, so GA-H1 fails | open | never |
| GA-H4 | The tile task is better framed as active/inactive classification than as strength regression | open | never |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Multi-resolution is a PARAMETERISATION of this family, not a separate family | 2026-10-07 | AlphaGenome item 1 at bacterial scale is kernel widths, not a 10 kb pyramid |
| MVP = an answered question, not a won bet | 2026-10-07 | Otherwise a recorded negative could never reach the bar |
| The tile arm is normalised by the TILE ceiling (0.9410), never the fragment ceiling (0.7928) | 2026-10-07 | The shared-normaliser error was caught and is pinned by this decision |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Ratify the drafted idea-anchor for this family | Soraya | user | `features/glm-learned-representation/idea-anchor.md`; idea-anchor is user-confirmed by standing directive |
| Whether a negative here also closes G-B or only defers it | Soraya | GA-H1 result | A measured capacity ceiling would close G-B; a merely-weak encoder would only defer it |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Fixed features are measured to fail on real genomic promoter tiles (GC 0.3000 beats every learned feature; best learned reaches 29.4% of the 0.9410 ceiling); no learned encoder has been built or tested.

### Target state / terminal condition
A committed artifact stating whether a learned multi-resolution encoder beats GC by >=0.05 on genomic peak tiles under leave-peak-out — with the falsifier honoured if it does not, which closes the learned-representation lever for this substrate rather than extending it.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` count + `gates-passed` count (raw counts, unweighted)
- 5 unknowns and 4 hypotheses at init; 0 retired.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Fix the 48-475/280 vs 48-499/359 artifact self-contradiction | edit-local-code | low | low | low | 15 min |
| 2 | Linear-headroom diagnostic (GC + dinucleotide + positional-GC on tiles) | run-tests | medium | **high** | low | 2 hr |
| 3 | Implement `glm/encoder.py` (conv + pool + linear head) | edit-local-code | high | medium | medium | 1 day |
| 4 | Gate script emitting a verdict artifact either way | edit-local-code | medium | low | low | 3 hr |
| 5 | Run tiles + grid, commit artifact, record verdict | run-tests | **high** | **high** | medium | 1 day |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** re-run `/project-state` after any action class fires.
- Immediate re-evaluation if action 2 shows 0.3000 is at a composition ceiling — that would confirm GA-H3 and make actions 3-5 a predicted-negative, which is a scope decision for the user.

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
| 1 | 2026-10-07 | propose | /project-init protocol executed for GLM family G-A | ledger created; 3a PASS (numbers verified against committed artifacts, one artifact self-contradiction found), 3b PASS project, 3c PASS |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Ratify or redirect the drafted idea-anchor** at `features/glm-learned-representation/idea-anchor.md`. Idea-anchor output is user-confirmed by standing directive, so it is drafted and parked rather than treated as anchored.
- **Is a recorded NEGATIVE an acceptable terminal for this family?** The MVP bar is written as "the question is answered", so a measured negative reaches MVP. If you want the family to stay open until a learned encoder actually wins, say so — that is an acceptance-bar redefinition and therefore yours, not mine.
- **Scope check:** candidate 3 (train on grid, test on tiles) is a genuinely different and harder question than the adopted goal. It is deliberately NOT in scope. Flag if you want it folded in.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-10-07
- **Progress signal:** (none yet — init only)
