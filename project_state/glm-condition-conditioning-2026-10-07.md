# GLM G-C — Condition-conditioned expression oracle
<!-- project-schema: 0.1 -->

> Initialized 2026-10-07. Project ID: glm-condition-conditioning-2026-10-07. Originating goal (verbatim user input): "Build the bacterial analogue of AlphaGenome's cell-state conditioning: one oracle conditioned on growth medium, tested against two separate per-medium oracles."

## Project Context
- **Project ID:** glm-condition-conditioning-2026-10-07
- **Project root:** C:\Users\Farshad\PythonProjects\dna_decode
- **Captured:** 2026-10-07
- **Originating goal:** One expression oracle conditioned on growth medium, versus two separate per-medium oracles.
- **Refined goal (3c):** Determine whether a single model taking (sequence, medium) beats two independent per-medium models on held-out position-blocked data, by >=0.02 Spearman, with each arm scored against its own replicate-derived ceiling.
- **Horizon (months):** 12 (nominal; real expected duration ~1 week)
- **Schema:** project-schema 0.1
- **Portfolio:** GLM umbrella, family **G-C**. Independent of G-A by construction — G-A changes the REPRESENTATION, G-C changes the CONDITIONING. See `plans/GLM_Portfolio_Decompose_2026-10-07.md`.

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** attempted
- **Provisional:** NO
- **Findings:** Verified against `wiki/glm_genomewide_oracle_2026-10-07.json`: within-condition single-replicate agreement LB **0.4578** / M9 **0.4774** (mean 0.4676); cross-condition single-vs-single mean **0.4062** over the four pairings; noise-matched gap **-0.0614**; n_shared **297,599**; unmatched average-vs-average **0.5025**; ceiling-if-condition-were-irrelevant **0.6371**; disattenuated true cross-condition **0.7887**; verdict `CONDITION_CARRIES_SIGNAL` — ALL CONFIRMED. The derived claim "~38% of real variance is condition-specific" is 1 - 0.7887^2 = 0.378, arithmetic on a confirmed figure, and is labelled derived rather than measured. No WebSearch performed: every claim is an internal measurement with a committed artifact, which a web search cannot adjudicate.

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Two media, one already-downloaded substrate (297,599 shared fragments), a named comparator (two separate models), a named margin (>=0.02), and a named split (`position_blocked_split`, already implemented). The motivating effect is already measured positive, which makes this the portfolio's lowest-risk family rather than an open-ended inquiry.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  1. **(ADOPTED)** One (sequence, medium) model vs two per-medium models, >=0.02 on held-out position-blocked data, each arm against its own ceiling.
  2. Cheaper first cut: a shared ridge with a medium indicator feature. Answers "does sharing help at all" before any architecture is involved, and is hours not days.
  3. Harder: predict the DELTA between media directly (sequence -> LB minus M9), which targets the condition-specific variance rather than the shared component. Higher value, much lower signal-to-noise.
  4. Rejected: pool both media into one training set with no medium input. Named only to exclude it — that discards the condition signal entirely and would look like a fair baseline while answering nothing.

## Goal Hierarchy
### Long-term (12+ months tier)
An oracle that answers "how strongly does this sequence express" conditioned on the environment, which is the prerequisite for a trait that depends on environment.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Shared-ridge + medium indicator | Beats two separate ridges, or does not, on held-out data | 3 days |
| 2 | Per-arm ceilings computed | LB and M9 ceilings derived separately; no shared normaliser | 1 day |
| 3 | Verdict + artifact | Committed artifact stating whether conditioning helps, either way | 1 week |
| 4 | Optional merge with G-A | If G-A produces an encoder, re-test conditioning on top of it | 1 month |

### Short-term (<=1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Compute LB and M9 ceilings separately and pin them | run-tests | Soraya | 2 hr |
| 2 | Implement the two-vs-one comparison with a shared position-blocked split | edit-local-code | Soraya | 1 day |
| 3 | Add a medium-shuffled null (shuffle the medium label) so a spurious gain is detectable | edit-local-code | Soraya | 3 hr |
| 4 | Run, commit artifact, record verdict | run-tests | Soraya | 1 day |

## MVP Criteria

Attempt budget: 3 per criterion. **MVP = an answered question**, not a won bet.

| # | Criterion | Kind | Predicate |
|---|---|---|---|
| 1 | The comparison script exists | file-exists | `scripts/glm_condition_conditioning.py` |
| 2 | Its tests pass | test-exit-0 | `uv run pytest tests/test_glm_condition_conditioning.py -q` |
| 3 | A verdict artifact exists | file-exists | `wiki/glm_condition_conditioning_2026-10-08.json` |
| 4 | The verdict is recorded here | project-state-row | Action Log row whose outcome names the verdict |
| 5 | The pre-registered control arm came out FLAT (validity) | project-state-row | Action Log row recording `null_clean=true` from the gate artifact |
| 6 | The frozen seed list + training protocol were used (validity) | project-state-row | Action Log row recording `registered_protocol_used=true` from the gate artifact |

**D1 (2026-10-07) — the two validity predicates above are LOAD-BEARING, not ceremony.** The original four
(module exists / tests exit 0 / artifact exists / ledger row) are ALL satisfied by a BROKEN RUN that writes a
null-result artifact. So the bar as first written made a *fake* negative reachable, which is worse than either
"require a win" (which would make a genuine recorded negative unreachable) or "accept any artifact". The
adopted rule is: **require the experiment to be provably VALID, then accept either outcome.** Both predicates
read a machine-readable field the gate script stamps into its own artifact — not a human claim.

**Path correction (2026-10-08).** Criterion 2 named `tests/test_glm_condition.py`; the test file written for
this family is `tests/test_glm_condition_conditioning.py` (the project's `test_<script>` convention, as with
`test_glm_tile_headroom.py`). As written the predicate was permanently unsatisfiable — **the second such
mis-specified path in this family's bar**, after criterion 3's `_2026-10-07.json`. It did fail CLOSED
(pytest exits 4 on a missing file, verified, so it would have blocked the MVP rather than passing vacuously),
which is the safe direction but still a bar that could never be met. Corrected to the real path; the
criterion itself ("its tests pass") is unchanged and the bar is not relaxed.

## State Snapshot
### Assumptions
- Two conditions is enough to show sharing helps, if it helps — **medium**. Two is the minimum; a null result at n=2 conditions is weak evidence against conditioning in general.
- The shared component dominates, so a shared trunk should win — **medium**; disattenuated cross-condition 0.789 says most variance IS shared.
- Position-blocked splitting is as necessary here as elsewhere — **high** (91.5% of consecutive fragments overlap, measured).
- The fragment substrate is usable despite its low ceiling (0.793) — **medium**; it is the only two-condition substrate available.

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | Condition carries signal beyond noise: cross 0.4062 below within 0.4676, gap -0.0614, n=297,599 | wiki/glm_genomewide_oracle_2026-10-07.json | high | 2026-10-07 |
| 2 | Disattenuated true cross-condition correlation 0.7887, so most variance is SHARED and ~38% is condition-specific (derived) | same artifact | high | 2026-10-07 |
| 3 | The unmatched average-vs-average comparison reads 0.5025 and INVERTS the answer | same artifact | high | 2026-10-07 |
| 4 | Fragment ceilings differ by medium: LB 0.7928, M9 0.8014 — so each arm needs its own | same artifact | high | 2026-10-07 |
<!-- project-state:end:evidence -->

### Unknowns
- Whether a medium indicator is enough, or whether conditioning needs to modulate the representation.
- Whether the ~38% condition-specific variance is concentrated in a subset of fragments (e.g. catabolic promoters) rather than spread.
- Whether two conditions can distinguish "sharing helps" from "more training data helps" — the shared model sees 2x the rows.
- Whether the peak-tile subset also has two media (it does not appear to), which would otherwise give a cleaner substrate.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| GC-H1 | A (sequence, medium) model beats two separate per-medium models by >=0.02 | open | never |
| GC-H2 | The gain, if any, survives a medium-shuffled null | open | never |
| GC-H3 | Any gain is explained by doubled training data rather than by conditioning | open | never |
| GC-H4 | Condition-specific variance is concentrated in catabolic/regulated promoters | open | never |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Each arm is scored against its OWN replicate-derived ceiling | 2026-10-07 | LB 0.7928 and M9 0.8014 differ; a shared normaliser rescales one arm |
| A medium-SHUFFLED null is mandatory, not optional | 2026-10-07 | GC-H3 is the live confound: the shared model sees 2x the data, so a gain is ambiguous without it |
| All comparisons are noise-matched | 2026-10-07 | The unmatched form inverted this family's own founding measurement |
| D1: MVP = an answered question, but ONLY behind validity predicates (`null_clean`, `registered_protocol_used`) | 2026-10-07 | The original four predicates are all satisfied by a broken run writing a null-result artifact; see wiki/glm_authority_decisions_2026-10-07.md |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Ratify the drafted idea-anchor | Soraya | user | `features/glm-condition-conditioning/idea-anchor.md` |
| Whether a 2-condition null result closes the family or defers it pending more conditions | Soraya | GC-H1 result | Two conditions is the minimum; a null is weak evidence against conditioning generally |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Condition is measured to carry real signal (noise-matched gap -0.0614 on 297,599 shared fragments; disattenuated true cross-condition 0.7887), and no conditioned model has been built.

### Target state / terminal condition
A committed artifact stating whether one (sequence, medium) model beats two separate per-medium models on held-out position-blocked data, with a medium-shuffled null ruling out the doubled-data explanation.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. 4 unknowns, 4 hypotheses at init; 0 retired.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Pin LB and M9 ceilings separately | run-tests | low | low | low | 2 hr |
| 2 | Shared-ridge + medium indicator vs two ridges | edit-local-code | **high** | **high** | low | 1 day |
| 3 | Medium-shuffled null | edit-local-code | medium | **high** | low | 3 hr |
| 4 | Run, commit, record verdict | run-tests | high | high | low | 1 day |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- **Default:** re-run `/project-state` after any action class fires.
- Immediate if action 3 shows the gain survives shuffling by less than the gain itself — that confirms GC-H3 and the headline must be withdrawn.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — REQUIRES per-action human approval

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-10-07 | propose | project-init protocol executed by hand for GLM family G-C | ledger created; 3a PASS (all figures verified against the committed artifact), 3b PASS project, 3c PASS |
| 2 | 2026-10-08 | run-tests | G-C 10-seed sweep executed: partial_pooling (PRIMARY) vs two_model on 297,868 coordinate-keyed LB/M9 fragments, position-blocked splits, pooling strength selected on an inner TRAIN-only split | **VERDICT NO_GAIN** — partial_pooling +0.2033 vs two_model +0.2032, gain **+0.0001 against a bar of 0.02** (two orders of magnitude short, not marginal). `null_clean=true` (shuffled gain −0.0006) and `registered_protocol_used=true`, both stamped by the script into `wiki/glm_condition_conditioning_2026-10-08.json`, so D1's two validity predicates hold on a machine-readable field rather than a human claim. Memo: `wiki/glm_condition_conditioning_result_2026-10-08.md`. THREE defects found en route, two of them mine: (a) the full-interaction arm is DEGENERATE with the comparator by construction (corr 0.9999998 on synthetic truth with medium-specific coefficients; equal to 4dp on real data at the median); (b) the pooling scale was a silent NO-OP — applied BEFORE StandardScaler, which divides out a constant column factor, so all three conditioned arms were secretly ONE model (r=0.9999970, same to 7dp); fixed by moving the penalty AFTER standardisation (block penalised by alpha/c²) and pinned by a test that asserts a pre-scaler scale is inert and a post-scaler one is not; (c) MVP criterion 2 named a test file that could never exist — failed closed (pytest exit 4), the second mis-specified path in this bar. The pooling curve is FLAT (0.0001 spread; argmax wanders 0.03/0.1/1.0 across seeds) so the selected scale is NOISE — pooling strength is not a lever here. Partial mechanistic account: CENTERING the unit-sum 3-mer block makes the design exactly rank 63 of 64 (cond 9.3e14), so ridge alpha is already load-bearing. Also measured: the pooled metric's excess over per-medium-mean IS the medium offset (shuffling collapses pooled 0.3085 → 0.2027 ≈ its own pmm 0.2026), upgrading an asserted caveat to a demonstrated one. Engineering: features were 97% of cost (181s vs 4s ridge); caching one pass per medium cut the sweep ~3.5h → ~35min, bit-faithful (two_model +0.2013/+0.1916 on seeds 0/1, reproduced on 3 separate runs) |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Ratify or redirect the drafted idea-anchor** at `features/glm-condition-conditioning/idea-anchor.md`.
- **Is two conditions enough to count?** If conditioning fails at n=2 media, that is weak evidence against conditioning in general. Whether that closes the family or defers it is an acceptance-bar call and therefore yours.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-10-07
- **Progress signal:** (none yet — init only)