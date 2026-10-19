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
| 3 | The gate ran and produced a verdict artifact | file-exists | `wiki/glm_learned_representation_2026-10-08.json` |
| 4 | The verdict is recorded in this ledger | project-state-row | Action Log row whose outcome names the verdict |
| 5 | The pre-registered control arm came out FLAT (validity) | project-state-row | Action Log row recording `null_clean=true` from the gate artifact |
| 6 | The frozen seed list + training protocol were used (validity) | project-state-row | Action Log row recording `registered_protocol_used=true` from the gate artifact |

**D1 (2026-10-07) — the two validity predicates above are LOAD-BEARING, not ceremony.** The original four
(module exists / tests exit 0 / artifact exists / ledger row) are ALL satisfied by a BROKEN RUN that writes a
null-result artifact. So the bar as first written made a *fake* negative reachable, which is worse than either
"require a win" (which would make a genuine recorded negative unreachable) or "accept any artifact". The
adopted rule is: **require the experiment to be provably VALID, then accept either outcome.** Both predicates
read a machine-readable field the gate script stamps into its own artifact — not a human claim.

**Hypothesis resolutions (2026-10-08), with the reasoning that is not obvious from the status word:**

- **GA-H1 FALSIFIED.** Margin +0.0122 against its own >=0.05 clause, and the paired CI includes zero, so it
  fails on both the frozen bar and the stronger instrument.
- **GA-H2 stays UNDER-INVESTIGATION rather than confirmed, because two instruments disagree.** Paired
  multi-vs-w20: 7/10 wins, mean +0.0114, normal-approx CI **[+0.0007, +0.0222]** which *barely* excludes
  zero — while the distribution-free sign test is **p=0.344**. A lower bound of +0.0007 on 10 paired
  observations with sign-flipping values is inside the fragility this project has already been bitten by
  (tie-handling and sort order silently shifting a median). Suggestive, underpowered, not established.
- **GA-H3 CONFIRMED**, and this is the load-bearing one: a conv encoder escapes both limits Step 3's proxy
  had (a precise-offset motif; a composition x position interaction) and still cannot reliably beat one GC
  number, so 0.3000 really does sit near a composition-only ceiling on this substrate.
- **GA-H4 stays OPEN.** The run shows the conclusion is framing-ROBUST (AUROC gives the same ordering and
  the same small edge), which is not the same as showing classification is the *better* framing. Spearman
  and AUROC are different scales and this run was not designed to rank them; asserting GA-H4 either way
  off it would be reading a question the design did not ask.

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
| GA-H1 | A learned encoder beats GC's 0.3000 by >=0.05 on genomic tiles under leave-peak-out | falsified | 2026-10-08 |
| GA-H2 | Multi-resolution kernels beat a single resolution at 150 bp | under-investigation | 2026-10-08 |
| GA-H3 | 0.3000 is near a composition-only ceiling, so GA-H1 fails | confirmed | 2026-10-08 |
| GA-H4 | The tile task is better framed as active/inactive classification than as strength regression | open | 2026-10-08 |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Multi-resolution is a PARAMETERISATION of this family, not a separate family | 2026-10-07 | AlphaGenome item 1 at bacterial scale is kernel widths, not a 10 kb pyramid |
| MVP = an answered question, not a won bet | 2026-10-07 | Otherwise a recorded negative could never reach the bar |
| The tile arm is normalised by the TILE ceiling (0.9410), never the fragment ceiling (0.7928) | 2026-10-07 | The shared-normaliser error was caught and is pinned by this decision |
| D1: MVP = an answered question, but ONLY behind validity predicates (`null_clean`, `registered_protocol_used`) | 2026-10-07 | The original four predicates are all satisfied by a broken run writing a null-result artifact; see wiki/glm_authority_decisions_2026-10-07.md |
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
| 2 | 2026-10-08 | run-tests | G-A encoder gate ran on the registered 10-seed list: gc / encoder-w20 / encoder-multi (6,20,75) / encoder-shuffled, on 44,106 peak tiles under leave-peak-out, CPU (the determinism reference) | **VERDICT WEAK by the frozen median rule, and the PAIRED analysis says that WEAK is NOT distinguishable from DOES_NOT_GENERALISE.** gc +0.3112 / encoder-multi +0.3234 = median margin **+0.0122** against a 0.05 bar; paired on identical splits, encoder-multi wins **6/10 seeds**, mean +0.0125, sign-test **p=0.754**, 95% CI **[−0.0099, +0.0349] INCLUDES ZERO**; encoder-w20 is a coin flip (5/10, p=1.000). The per-seed margins flip sign (+.024 +.038 +.022 +.032 −.032 −.022 −.020 −.020 +.018 +.085), so the positive median is carried by seed variation, not a reliable edge. `null_clean=true` (shuffled arm −0.0240, loses 10/10, its own CI **excludes** zero at p=0.002 — so the paired test CAN resolve a real difference on 10 seeds and simply does not resolve one here), `registered_protocol_used=true`, `reconcile_ok=true` with the anchor landing on **0.3000 exactly**. BOTH framings agree: AUROC on the active column 0.6774 vs 0.6672, the same small unreliable edge — so this is not an artifact of framing regression over 29,193 inactive tiles. Artifact `wiki/glm_learned_representation_2026-10-08.json`, memo `wiki/glm_encoder_gate_result_2026-10-08.md`. **Built on the user's explicit decision AGAINST Step 3's predicted-negative, and that was the right call**: Step 3's proxy was positional GC in 10 bins under a LINEAR ridge, so it could neither see a motif at a precise offset nor represent a composition×position interaction; an encoder escapes both limits and STILL does not win, which upgrades a predicted-negative into a measured one. The encoder is not broken — it learns a planted fixed-offset motif to Spearman 0.84 against a computed ceiling of 1.0, its null collapses, and it tracks the incumbent closely rather than failing to train. Honest scope: THIS encoder at THIS size (32 ch/branch, 13,217 params) on THIS substrate; capacity has already hurt here (6-mers scored worst at 0.2022), so a larger encoder is a different experiment. |
| 3 | 2026-10-11 | run-tests | R5/R6 CEILING arm executed: 4-scale multi-resolution pyramid (widths 3/12/48/150 = 1.70 decades) vs the 1.10-decade encoder-multi, 10 REGISTERED seeds, leave-peak-out, CPU | **NO GAIN** - pyramid median 0.3167 / PAIRED +0.0069 vs GC (6/10 wins, CI95 [-0.0124,0.0262] INCLUDES ZERO) against encoder-multi 0.3234 / +0.0125 at 2.1x the parameters (27,649 vs 13,217): widening the receptive field 1.10->1.70 decades and doubling capacity came out slightly WORSE. VALID negative - own shuffled null clean (median -0.0198, paired -0.3368, 0/10 wins, sign-p 0.002), registered_protocol_used=True, null_clean=True, reconcile_ok=True. Honest reading is 'no arm separates from GC' (all three encoder CIs include zero), not 'multi beats pyramid'. Corroborates the already-recorded capacity-hurts pattern (6-mers/4,096 feats worst at 0.2022) rather than being independent. SECOND AND LARGER FINDING: the declared arm `multiresolution-pyramid-4-decades-UNRUN` is SUBSTRATE-GATED, not expensive - 4 decades needs ~10,000bp and a census of the real substrate gives tile max 150 / neg_control 150 / random 150 / pos_control 48-75, so the absolute span ceiling is log10(150)=2.18 decades even with a degenerate 1bp branch. Its prerequisite is the OTHER planned arm (substrate-change), so the two planned arms are ORDERED, not parallel - a dependency the arms ledger does not encode. Arm left PLANNED, NOT killed (R6: an arm dies only by its own executed evidence; what was measured here is the substrate). New arm pyramid-4scale-1.70decades-150bp recorded measured=0.3167; `best` correctly stays conv-3width at 0.3234. Artifacts: wiki/glm_pyramid_arm_2026-10-11.{json,md}. --pyramid is ADDITIVE: without it the three published arms and the reconcile gate are byte-identical (28 tests pass). |
| 4 | 2026-10-12 | edit-local-code | Corrected two over-claims in wiki/glm_pyramid_arm_2026-10-11.md and pinned the substrate gate with tests. Retired phrase 1 VERBATIM: "This closes one of them with a number and reclassifies the other" -- it closed NEITHER declared arm; it measured a THIRD arm, contradicting the memo's own "the arm measured here is not either declared arm". Retired phrase 2 VERBATIM: "more scale is not the missing lever" -- a lever-level negative from one encoder family at one capacity on one substrate, over-reaching the memo's own honest_limits. Both now appear ONLY inside a dated ## Correction 2026-10-12 section, never as live claims. Added 3 tests to tests/test_glm_encoder.py: PYRAMID_WIDTHS pinned at (3,12,48,150) with its 1.70-decade span asserted; a FOUR-decade pyramid shown substrate-gated by arithmetic (>=4 decades needs a >=10,000bp coarsest branch vs a 150bp tile) rather than merely expensive; and the 4-branch feature width. PYRAMID_WIDTHS previously had ZERO test references. Plan Step 3 (record the scoped negative) was SKIPPED as already-satisfied: neg-5ad65fa5 recorded 2026-10-09 already carries both the unrunnable scope and the revival condition, so a second record would duplicate it. Plan Step 5 (arms-ledger resolution) NOT RUN -- authority-gated, awaiting user decision. | 64 passed in tests/test_glm_encoder.py + test_glm_encoder_gate.py (baseline 61, +3 new). Both new guards proven NON-VACUOUS by temporary mutation: widening PYRAMID_WIDTHS to include 1000 fails the pin; disabling the too_wide refusal fails the 4-decade test. encoder.py reverted and verified byte-identical to HEAD (19882 bytes). Frozen five byte-unchanged; prospective lock (2026-08-31 v2) re-verified ok=True, zero drift. Criterion still blocked:untried-arm -- unchanged, because resolving it is the authority-gated step. |
| 5 | 2026-10-12 | run-tests | Epilogue committed (63a270e, 6 files, unpushed). +2 tests pinning the pyramid's 27,649 params (the memo's '2.1x the parameters' figure, previously untested -- and the existing small-by-design guard covers DEFAULT_WIDTHS only, so it never covered this arm) and the 1.10/1.70/2.18 decade ladder (1.70 was pinned earlier; 1.10, the baseline it is compared against, was trusted). CLAUDE.md stale suite count 6365/13 -> 6407/14 corrected. NEXT.md item 8 added = the authority decision for the ceiling criterion, re-derived from the run JSON not the memo prose (which fixed two narrative errors: the 4-decade arm is left PLANNED not killed per R6, and is ORDERED behind the substrate-change arm, a dependency the ledger does not encode). A CLAUDE.md pointer for this arm was REFUSED, not forgotten: 8 words of headroom, the ceiling was raised to fit and then exceeded, so both raise and addition were reverted; WORD_CEILING byte-identical to HEAD. Paying by compression is reserved to the user (NEXT.md item 6). | 71 passed across test_glm_encoder / test_glm_encoder_gate / test_claude_md_growth_ceiling. Both new tests proven non-vacuous by mutation then reverted. Frozen five byte-unchanged. encoder.py identity claimed by hash fe71bd4 after finding that `git checkout --` rewrote line endings 19882->20306 bytes while `git status` still read clean. Concurrent sessions' pgx/capstone/feature edits left unstaged. Step 5 of the plan still open on the user's A/B/C decision. |
| 6 | 2026-10-12 | edit-local-code | Corrected two UNDER-claims in the auto-loaded orientation files (1229b74). F3: the phylogeny de-confounder was calibrated 2026-10-10 on 1001 Genomes and eval/regime.py records it correctly, but CLAUDE.md carried 'BUILT AND NOT CALIBRATED' as a bullet headline citing the superseded 10-09 artifact, and NEXT.md contradicted ITSELF (F3 row said NOT CALIBRATED while a row 12 lines above credited the 1001G calibration). Three surfaces, three answers, auto-loaded one wrong. Both rewrites now DEFER to regime.py::DECONFOUNDER_CALIBRATION instead of restating status -- a doc that restates a status goes stale silently, and this one did for two days. Also fixed a figure THIS branch invalidated: CLAUDE.md said 6407 passed, true until 63a270e added 2 tests. | Full clean suite 6409 passed / 14 skipped / 0 failed, exit 0 (23:36) = predicted 6407+2. CLAUDE.md edit is word-NEGATIVE: bullet 219->218, file 35452->35451, ceiling headroom 8->9 words. ENFORCED not asserted -- an inline word-budget assert refused three drafts at +47/+12/+3 before one fit. Nothing compressed to make room (that stays the user's call, NEXT.md item 6). Doc+citation guards 106 passed. Frozen five + regime.py byte-unchanged. NOT DONE, both authority: plan Step 5 arms ledger A/B/C, and the CLAUDE.md layered-memory restructure (43 bullets / 24,677 words compressible, 5 orphans / 2,035 words need writing out first). |
| 7 | 2026-10-13 | edit-local-code | **Reconciled two sigma70 loaders; the duplicate reproduced the exact failure the original exists to prevent.** `dna_decode/glm/constructed.py` + its 12 tests DELETED (neither ever committed, no importers outside its own test). `dna_decode/glm/expression.py` was already canonical: already carried `_norm_key` (hyphen/underscore, with the mutation arrow `->` PROTECTED), already enforced `MIN_JOIN_RATE = 0.95`, and `tests/test_glm_expression.py` already pinned `rep.naive_join == 3821` labelled *"the naive-join counterfactual is the point of the fix"* PLUS a live assertion that every element column has >=3 levels. So the correct figure, the wrong figure, the cause, the fix and two green tests asserting both were ALL in the repo when the wrong number shipped -- not missing information. **FOUR corrected numbers:** sigma70 usable 3,821 -> **10,898 (join rate 1.00)**; the 65% gap is a separator mismatch NOT absent data (open question CLOSED, nothing was missing); design 8x8x8x8 with UP_element collapsed to 1 -> **3x8x8x8x8 = 12,288 cells / 10,898 filled / UP_element has 3** (5-way crossed); and fragment ceilings 0.4583/0.4731 -> **0.7928/0.8014** (those were `single_replicate_agreement`, not `ceiling_for_a_sequence_model`). **The ceiling error REVERSED a recommendation:** the frag rung trades 7.3x the data for ~85% of headroom, not ~49%, so the wrong number argued against the right next move. Peak tiles measured 0.7945 -> 0.9410 on the 44,106 filtered rows. The duplicate's ONE original contribution (`irrelevant_feature_control`, non-vacuous by design) moved into `expression.py` and GENERALISED from Background-only to any of the 5 element columns. **The ambiguous JSON key `replicate_ceiling` was DELETED from every rung, not re-valued** -- one key standing for two quantities is what made the conflation possible. Plan Step 5 marked VOID with the body left unedited as the record of what was planned on a false premise. **Arms ledger corrected:** the measured best was `conv-3width` 0.3234 while the biophysical model's **0.3968** from 10-12 was never registered; now recorded with full evidence, and `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` DECLARED (planned/medium) so it cannot evaporate. **Filing the fix as an ARM was rejected on inspection:** `measure_arm` requires `metric: float`, so a bug-fix cannot be an arm without faking a number, and it would have left the predictor-ranking ledger `blocked:untried-arm` on a maintenance ticket. Multi-branch fan-out also rejected -- five artifacts all restating one figure, strictly sequential. | 61 passed on the affected subset (test_glm_expression 23 + test_glm_encoder 38); `test_glm_expression.py` 23/23 inside the full run with **zero F/E markers** through 60%, full suite still in flight at commit time. Collected **6,421** = 6,415 committed baseline + 6 new (the 12 deleted were never committed). **Duplicate-loader audit RUN and the honest-limit CLOSED same day:** exactly ONE loader per GEO series (`expression.py` owns GSE108535, `genomewide.py` owns GSE144621, each with its own `GEO_BASE`); only residue is `DEFAULT_CACHE` defined identically in both -- recorded, not fixed, since consolidating would couple two independent substrate loaders. **A silent no-op was caught:** `measure_arm` mutates memory only and `save()` is separate, so the first registration exited unsaved and the verdict read unchanged -- indistinguishable from success. Criterion correctly STAYS `blocked:untried-arm` (three declared arms unrun). Artifacts `wiki/sigma70_loader_reconciliation_2026-10-13.{md,json}`. Concurrent sessions' pgx/capstone edits left unstaged. |
| 8 | 2026-10-19 | edit-local-code | **Closed the 48-475/280 vs 48-499/359 self-contradiction, and the ranker's 'already-done' advisory was a FALSE POSITIVE worth reading.** `advance_ranker` flagged this candidate at a 100% term match, but the match is LEXICAL against the ledger prose, not evidence of completion -- read per v2.1 rather than auto-skipped, and the contradiction was still live. What WAS already done is the PRODUCER: `scripts/glm_genomewide_oracle.py:205-210` already derives the string from `load_stats` and carries a comment saying so. Only the COMMITTED ARTIFACT was stale, because regenerating it needs the absent `D:/dna_decode_cache/mpra` substrate -- so the fix was applied to the artifact from **its own** `load_stats.LB` (48/499/359), using no new data and re-running nothing. **The edit is RECORDED, not silent:** a `post_hoc_corrections` entry carries the old string, the new one, and why, because silently rewriting a committed run record is worse than the contradiction it repairs -- a reader could no longer tell a measurement from a later repair. **The real gap was that NOTHING PINNED IT:** the producer fix was unguarded, so a future edit could re-hardcode the prose silently. Added 3 artifact-level guards that catch the defect CLASS (prose disagreeing with the derived data beside it) rather than one string's spelling: the numbers in `honest_limits` must equal `load_stats.LB`; a non-vacuity check that the regex really extracts three comparable numbers AND that the superseded (48,475,280) sample figure has not returned; and a check that any post-hoc edit is recorded with its old value. [done:52c5f50e] | 23 passed / 2 skipped in `tests/test_glm_genomewide.py` (25 collected, was 22); the 2 skips are the real-substrate tests, correctly skipped with GSE144621 absent. Non-vacuity is asserted in-test, not claimed. Artifact numbers UNTOUCHED -- only a prose restatement changed, and `post_hoc_corrections.numbers_unchanged` says so. Frozen five byte-unchanged. Also surfaced from this run, NOT fixed: (a) `glm-generator-finetune` is `blocked_by ['sampling sweep + 8 gb gpu']`, free text rather than a family slug, so the ranker fail-closes it as 'not an accepted family' permanently -- a ledger-shape defect on the one family carrying the `high` ambition bucket; (b) R15 warns AMBITION CONCENTRATED in a single family; (c) 686/686 wiki docs carry no `doc-status` header, so nothing separates a live claim from a historical record. |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Ratify or redirect the drafted idea-anchor** at `features/glm-learned-representation/idea-anchor.md`. Idea-anchor output is user-confirmed by standing directive, so it is drafted and parked rather than treated as anchored.
- **Is a recorded NEGATIVE an acceptable terminal for this family?** The MVP bar is written as "the question is answered", so a measured negative reaches MVP. If you want the family to stay open until a learned encoder actually wins, say so — that is an acceptance-bar redefinition and therefore yours, not mine.
- **Scope check:** candidate 3 (train on grid, test on tiles) is a genuinely different and harder question than the adopted goal. It is deliberately NOT in scope. Flag if you want it folded in.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-10-07
- **Progress signal:** (none yet — init only)