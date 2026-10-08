# GLM G-E — Distribution alignment (generator ↔ oracle)
<!-- project-schema: 0.1 -->

> Initialized 2026-10-07. Project ID: glm-distribution-alignment-2026-10-07. Originating goal (verbatim user input): "Make sure the generator and the oracle are validated on the SAME sequence distribution, because right now they are not — the oracle works on grid-like sequence and fails on genome-like, and nobody has characterised which one the generator emits."

## Project Context
- **Project ID:** glm-distribution-alignment-2026-10-07
- **Project root:** C:\Users\Farshad\PythonProjects\dna_decode
- **Captured:** 2026-10-07
- **Originating goal:** Ensure generator and oracle are validated on the same sequence distribution.
- **Refined goal (3c):** Measure WHICH distribution generated sequence falls in — designed-grid-like, genome-like, or neither — using a stated discriminator, and report an oracle number on THAT distribution. If the answer is "neither", say so, because that invalidates both existing oracle validations for loop purposes.
- **Horizon (months):** 12 (nominal; real duration ~1 week once unblocked)
- **Schema:** project-schema 0.1
- **Portfolio:** GLM umbrella, family **G-E**. **blocked_by: `glm-learned-representation-2026-10-07` (G-A) + `glm-generator-finetune-2026-10-07` (G-D).** This is the INTEGRATION family and the critical-path endpoint.

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** attempted
- **Provisional:** NO
- **Findings:** The goal's load-bearing claim is a CONJUNCTION of two confirmed measurements rather than a new claim. Oracle side, verified against `wiki/glm_genomewide_oracle_2026-10-07.json`: on genomic peak tiles GC 0.3000 beats one-hot 0.2691 / 4-mer 0.2769 / 6-mer 0.2022, best learned at 29.4% of the 0.9410 ceiling, verdict `DOES_NOT_GENERALISE`; and on the designed grid (`wiki/glm_expression_oracle_result_2026-10-07.md`) one-hot 0.5914 with GC only 0.1414 — so the ranking INVERTS between substrates. Generator side, verified against `wiki/glm_generator_falsifier_REAL_2026-10-07.json`: distinguishability 0.7706 vs the Markov null 0.6679, `WORSE_THAN_MARKOV_NULL`. **The claim that nobody has characterised the generator's output distribution is an ABSENCE claim**, checked by grepping the wiki for any artifact scoring generated sequence against either substrate's distribution: none exists. An absence verified by one search route is weaker than a positive measurement and is labelled as such.

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Bounded: one question (which distribution), one stated instrument (a discriminator already built and self-validating), one deliverable (an oracle number on the matching distribution). The predecessors are families, not open scope.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  1. **(ADOPTED)** Classify generated sequence against both reference distributions with a stated discriminator; report an oracle number on the matching one; if neither matches, report that.
  2. Cheaper precursor, doable NOW with no predecessors: characterise the two REFERENCE distributions against each other — how distinguishable is designed-grid sequence from genomic tile sequence? If they are barely distinguishable, the whole alignment worry is smaller than it looks and the family can shrink.
  3. Stronger: make the generator CONDITIONAL on the target distribution (prompt with grid-like or genome-like context) so alignment is achieved by construction rather than measured after the fact.
  4. Rejected: assume generated sequence is genome-like because the generator was pre-trained on genomes. Named to exclude — pre-training corpus is not the same as sampling distribution, and that inference is exactly what this family exists to replace with a measurement.

## Goal Hierarchy
### Long-term (12+ months tier)
A generative loop whose proposal half and scoring half are validated on the same sequence distribution, so a ranked edit list means what it says.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Reference-vs-reference distinguishability (candidate 2) | A number for how far apart the two substrates are — DOABLE NOW | 3 days |
| 2 | G-A and G-D land | An oracle validated on genomic sequence + a generator worth characterising | blocked |
| 3 | Generated-sequence classification | Which reference distribution generated output falls in, with a stated discriminator | +3 days |
| 4 | Aligned oracle number | An oracle score on the matching distribution, or an explicit "neither" | +1 week |

### Short-term (<=1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Candidate 2: measure grid-vs-tile distinguishability with the existing falsifier. NO predecessors needed | run-tests | Soraya | 3 days |
| 2 | WAIT on G-A and G-D for the rest | stop | Soraya | blocked |
| 3 | (on unblock) classify generated sequence against both references | run-tests | Soraya | +3 days |

## MVP Criteria

Attempt budget: 3 per criterion. **MVP = an answered question.** Criterion 1 is reachable NOW via
candidate 2 and deliberately does not depend on the predecessors, so this family is not fully idle.

| # | Criterion | Kind | Predicate |
|---|---|---|---|
| 1 | Reference-vs-reference distinguishability measured (no predecessors) | file-exists | `wiki/glm_reference_distribution_distance_2026-10-07.json` |
| 2 | The alignment script exists | file-exists | `scripts/glm_distribution_alignment.py` |
| 3 | Its tests pass | test-exit-0 | `uv run pytest tests/test_glm_distribution_alignment.py -q` |
| 4 | The verdict is recorded here | project-state-row | Action Log row whose outcome names the verdict |

## State Snapshot
### Assumptions
- The existing naturalness discriminator can also separate grid-like from genome-like sequence — **medium**. It was built to separate generated from natural, which is a different axis.
- "Which distribution" is a well-posed question with a binary-ish answer — **low**. Generated sequence may sit somewhere between, or outside both, which is why "neither" is an explicit allowed outcome.
- An oracle validated on the matching distribution will then be usable in the loop — **medium**. It also has to be ACCURATE there, which is G-A's job, not this family's.

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | The oracle's feature ranking INVERTS between substrates (grid one-hot 0.5914 / GC 0.1414; tiles GC 0.3000 / one-hot 0.2691) | wiki/glm_genomewide_oracle_result_2026-10-07.md | high | 2026-10-07 |
| 2 | The generator is below a 3-mer Markov chain at default sampling (0.7706 vs 0.6679) | wiki/glm_generator_falsifier_REAL_2026-10-07.json | high | 2026-10-07 |
| 3 | A self-validating discriminator already exists and REFUSES when its own controls fail | dna_decode/glm/falsifier.py | high | 2026-10-07 |
| 4 | No artifact characterises generated sequence against either substrate's distribution (ABSENCE, one search route) | wiki/ grep | medium | 2026-10-07 |
<!-- project-state:end:evidence -->

### Unknowns
- How far apart the two reference distributions actually are — candidate 2 answers this with no predecessors.
- Whether generated sequence falls in either, between, or outside both.
- Whether "distribution" is better operationalised as k-mer composition, positional structure, or discriminator-separability — they can disagree.
- Whether conditional generation (candidate 3) sidesteps the whole measurement.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| GE-H1 | Designed-grid and genomic-tile sequence are strongly distinguishable from each other | open | never |
| GE-H2 | Generated sequence falls closer to the genomic distribution than the grid distribution | open | never |
| GE-H3 | Generated sequence falls in NEITHER, invalidating both oracle validations for loop use | open | never |
| GE-H4 | Conditional prompting aligns the generator to a chosen distribution by construction | open | never |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| "Neither" is an explicit allowed outcome | 2026-10-07 | Forcing a binary answer would hide the most consequential result |
| Candidate 2 (reference-vs-reference) is unblocked and runs first | 2026-10-07 | It needs no predecessor and could shrink the whole family |
| Pre-training corpus is NOT evidence of sampling distribution | 2026-10-07 | That inference is what this family replaces with a measurement |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Ratify the drafted idea-anchor | Soraya | user | `features/glm-distribution-alignment/idea-anchor.md` |
| If generated sequence is "neither", does the loop get a NEW oracle substrate or does the generator get constrained? | Soraya | GE-H3 result | A strategic fork: new substrate is data work, constraining the generator is modelling work |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
The oracle is validated on grid-like sequence and measured to fail on genome-like sequence; the generator's output distribution is uncharacterised; so the loop's two halves are validated on different distributions and nobody has measured which one matters.

### Target state / terminal condition
A committed artifact naming which reference distribution generated sequence falls in, by a stated instrument, plus an oracle number on that distribution — or an explicit "neither", which is itself the deliverable.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. 4 unknowns, 4 hypotheses at init; 0 retired. **One unknown is retirable without any predecessor** (candidate 2).

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Reference-vs-reference distinguishability (NO predecessors) | run-tests | medium | **high** | low | 3 days |
| 2 | WAIT on G-A + G-D | stop | none | none | n/a | 0 |
| 3 | (on unblock) classify generated sequence | run-tests | **high** | **high** | high | 3 days |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- After candidate 2 reports. If the two reference distributions are barely distinguishable, GE-H1 is falsified and this family shrinks to a footnote — which would be a cheap and valuable outcome.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — REQUIRES per-action human approval

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-10-07 | propose | project-init protocol executed by hand for GLM family G-E | ledger created; 3a PASS (conjunction of confirmed measurements; the absence claim labelled one-route), 3b PASS project, 3c PASS |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Ratify or redirect the drafted idea-anchor** at `features/glm-distribution-alignment/idea-anchor.md`.
- **The strategic fork, if GE-H3 fires:** if generated sequence matches NEITHER reference distribution, the choice is a new oracle substrate (data work) or constraining the generator (modelling work). That is a direction call and therefore yours, not mine.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-10-07
- **Progress signal:** (none yet — init only)
