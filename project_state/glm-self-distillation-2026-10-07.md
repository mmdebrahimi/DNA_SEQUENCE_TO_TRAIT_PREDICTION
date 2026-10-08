# GLM G-B — Self-distillation (BLOCKED by G-A, measured)
<!-- project-schema: 0.1 -->

> Initialized 2026-10-07. Project ID: glm-self-distillation-2026-10-07. Originating goal (verbatim user input): "Adopt AlphaGenome's phase-2 self-distillation — perturb the input and train on the model's own outputs as labels — to manufacture a variant-effect training signal where only reference-sequence labels exist."

## Project Context
- **Project ID:** glm-self-distillation-2026-10-07
- **Project root:** C:\Users\Farshad\PythonProjects\dna_decode
- **Captured:** 2026-10-07
- **Originating goal:** Adopt AlphaGenome's self-distillation to manufacture a variant-effect signal from perturbed inputs labelled by the model's own outputs.
- **Refined goal (3c):** AFTER G-A delivers a higher-capacity oracle, determine whether distilling it on perturbed sequences beats the teacher on held-out REAL MEASURED data by >=0.02 — never on distillation agreement, which a degenerate student also achieves.
- **Horizon (months):** 12 (nominal; real duration ~3 days ONCE UNBLOCKED)
- **Schema:** project-schema 0.1
- **Portfolio:** GLM umbrella, family **G-B**. **blocked_by: `glm-learned-representation-2026-10-07` (G-A) — HARD and MEASURED.**

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** attempted
- **Provisional:** NO
- **Findings:** Verified against `wiki/glm_selfdistill_probe_2026-10-07.json`: a linear student on the same feature basis reproduces its teacher at Spearman agreement **1.0000** with delta_rho **-0.0** at mutation rates 0.02 / 0.05 / 0.10; a properly-trained non-linear student lands **0.5933-0.5992** against the teacher's **0.5963** (straddles, never exceeds); the GC-only non-vacuity control collapses **0.5963 -> 0.1736**; verdict `DEGENERATE_DISTILLATION_BLOCKED_UNTIL_LEARNED_REPRESENTATION` — ALL CONFIRMED. **The EXTERNAL method claim** (AlphaGenome's phase 2 perturbs inputs and trains on the pre-trained model's own outputs, lifting the held-out-region restriction) comes from the DeepMind technical talk `rdBtxtcS4nM`, read in full, and is recorded in `wiki/alphagenome_integration_analysis_2026-10-07.md`. It is a primary-source paraphrase, not a quotation, and is flagged as such.

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Bounded and SMALL — once a higher-capacity oracle exists this is a training-loop variant, not a new system. The blocker is a predecessor, not unbounded scope.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  1. **(ADOPTED)** Distil G-A's encoder on perturbed sequences; headline is held-out performance on REAL measured data, >=0.02 over the teacher.
  2. Narrower: use distillation only to SMOOTH the oracle off-distribution, measured as reduced variance on generated sequence, with no claim of improved accuracy. Honest and weaker; also the only version testable without new labels.
  3. Rejected: report distillation AGREEMENT as success. Named to exclude — agreement 1.0000 is exactly what the measured degenerate case achieves.
  4. Rejected: distil without a real-data validation half. Named to exclude — at DeepMind the distilled model is validated on real measured benchmarks, and adopting the first half alone is circular by construction.

## Goal Hierarchy
### Long-term (12+ months tier)
An oracle that is self-consistent under perturbation, so the generative loop can query it off the measured distribution without the score silently degrading.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | G-A delivers a higher-capacity oracle | A learned encoder exists whose capacity exceeds its own feature span | blocked |
| 2 | Perturb-and-relabel loop implemented | Student trains on teacher-labelled perturbed sequence | +2 days |
| 3 | Real-data validation half | Student scored on held-out REAL measured data, never on agreement | +1 day |
| 4 | Verdict + artifact | Committed artifact stating whether distillation beats the teacher, either way | +3 days |

### Short-term (<=1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | WAIT on G-A. No code. The degeneracy is measured, so early work is predicted-void | stop | Soraya | blocked |
| 2 | (on unblock) Reuse `scripts/glm_selfdistill_probe.py` with the G-A encoder as teacher+student | edit-local-code | Soraya | +2 days |
| 3 | (on unblock) Keep the GC-only non-vacuity control — a probe that cannot detect an impoverished student detects nothing | edit-local-code | Soraya | +3 hr |

## MVP Criteria

Attempt budget: 3 per criterion. **MVP = an answered question.** Criterion 1 is the explicit block: it
cannot be met before G-A lands, which is the correct behaviour rather than a defect.

| # | Criterion | Kind | Predicate |
|---|---|---|---|
| 1 | The predecessor's encoder exists (THE BLOCK) | file-exists | `dna_decode/glm/encoder.py` |
| 2 | The distillation probe still passes with the new teacher | test-exit-0 | `uv run pytest tests/test_glm_selfdistill.py -q` |
| 3 | A verdict artifact exists | file-exists | `wiki/glm_self_distillation_verdict.json` |
| 4 | The verdict is recorded here | project-state-row | Action Log row whose outcome names the verdict |

## State Snapshot
### Assumptions
- A conv encoder has capacity its own feature span lacks, so distillation is non-degenerate for it — **medium**. This is the premise the block rests on and it is UNTESTED for the CNN.
- Distillation cannot create information the teacher lacks — **high** (measured: the non-linear student straddled the teacher at 0.5933-0.5992 and never exceeded it).
- Uniform random point substitution is a reasonable perturbation — **low**. A generator's candidates are not uniform, so this tests the method's mechanics rather than the real off-distribution gap.

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | Linear student reproduces its teacher at agreement 1.0000 at every mutation rate | wiki/glm_selfdistill_probe_2026-10-07.json | high | 2026-10-07 |
| 2 | A trained non-linear student straddles the teacher (0.5933-0.5992 vs 0.5963) and cannot exceed it | same artifact | high | 2026-10-07 |
| 3 | The GC-only control collapses 0.5963 -> 0.1736, so the probe is non-vacuous | same artifact | high | 2026-10-07 |
| 4 | The method's purpose is to manufacture a variant-effect signal where only reference labels exist | wiki/alphagenome_integration_analysis_2026-10-07.md (primary-source paraphrase) | medium | 2026-10-07 |
<!-- project-state:end:evidence -->

### Unknowns
- Whether a conv encoder's capacity makes distillation non-degenerate, or whether the information bound dominates regardless of capacity.
- What perturbation distribution matches the generator's actual output — unknowable until G-E characterises it.
- Whether the real gain, if any, is accuracy or only off-distribution smoothness (candidate 2).

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| GB-H1 | Distillation is degenerate on a LINEAR teacher with a fixed basis | **confirmed** | 2026-10-07 |
| GB-H2 | Distillation is non-degenerate once the student has capacity the teacher lacks | open | never |
| GB-H3 | A distilled student beats its teacher on held-out REAL data by >=0.02 | open | never |
| GB-H4 | The information bound holds regardless of capacity, so GB-H3 fails for any student | open | never |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| G-B is blocked on G-A by MEASUREMENT, not preference | 2026-10-07 | Student-teacher agreement 1.0000; early work is predicted-void |
| The headline may never be distillation agreement | 2026-10-07 | Agreement is what the degenerate case achieves |
| The real-data validation half ships WITH the distillation half or not at all | 2026-10-07 | Adopting the first half alone is circular |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Ratify the drafted idea-anchor | Soraya | user | `features/glm-self-distillation/idea-anchor.md` |
| Whether a G-A negative also closes G-B | Soraya | G-A result | A measured capacity ceiling closes it; a merely-weak encoder only defers it |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
BLOCKED by measurement: self-distillation is degenerate on the current linear oracle (student-teacher agreement 1.0000), so no work proceeds until G-A delivers a higher-capacity teacher.

### Target state / terminal condition
A committed artifact stating whether distilling a higher-capacity oracle on perturbed sequence beats the teacher on held-out real measured data — with the information bound honoured if it does not.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. 3 unknowns, 4 hypotheses at init; **1 hypothesis already confirmed (GB-H1)** before any work, which is what justifies the block.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | WAIT on G-A (no code) | stop | none | none | none | 0 |
| 2 | (on unblock) distil the G-A encoder, real-data headline | edit-local-code | **high** | **high** | high | 2 days |
| 3 | (on unblock) keep the GC non-vacuity control | edit-local-code | low | medium | low | 3 hr |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- When `dna_decode/glm/encoder.py` exists (criterion 1 flips), re-evaluate immediately.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — REQUIRES per-action human approval

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-10-07 | propose | project-init protocol executed by hand for GLM family G-B | ledger created; 3a PASS (degeneracy figures confirmed), 3b PASS project, 3c PASS; family opens BLOCKED with GB-H1 already confirmed |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Ratify or redirect the drafted idea-anchor** at `features/glm-self-distillation/idea-anchor.md`.
- **Confirm the block is acceptable.** This family deliberately does nothing until G-A lands. If you want it attempted anyway, say so — but the degeneracy is measured, so the expected outcome is a null that costs real time.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-10-07
- **Progress signal:** (none — opens blocked by design)
