# GLM G-D — Generator fine-tuning (the un-pulled lever)
<!-- project-schema: 0.1 -->

> Initialized 2026-10-07. Project ID: glm-generator-finetune-2026-10-07. Originating goal (verbatim user input): "Fine-tune GENERator-v2-prokaryote-1.2b on real E. coli promoter windows rather than consuming it zero-shot, and measure whether it gets below the 3-mer Markov null."

## Project Context
- **Project ID:** glm-generator-finetune-2026-10-07
- **Project root:** C:\Users\Farshad\PythonProjects\dna_decode
- **Captured:** 2026-10-07
- **Originating goal:** Fine-tune GENERator-1.2B on real *E. coli* promoter windows rather than consuming it zero-shot.
- **Refined goal (3c):** Determine whether LoRA fine-tuning on MG1655 upstream-of-CDS windows brings GENERator's distinguishability BELOW the 3-mer Markov null's 0.6679 under feature mode `both`, on the frozen natural comparison set, with the corpus fingerprint stamped so the comparison is auditable.
- **Horizon (months):** 12 (nominal; real duration ~1 week once hardware returns)
- **Schema:** project-schema 0.1
- **Portfolio:** GLM umbrella, family **G-D**. **blocked_by: the sampling sweep reporting first** (so sampling and fine-tuning stay separable) **+ 8 GB GPU availability** (the server PC is mid-move).

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** attempted
- **Provisional:** NO
- **Findings:** Verified against `wiki/glm_generator_falsifier_REAL_2026-10-07.json`: GENERator-v2-prokaryote-1.2b at default sampling scores distinguishability **0.7706** (mode `both`) and **0.7237** (mode `kmer`) against the 3-mer Markov null's **0.6679** / **0.5990** — verdict `WORSE_THAN_MARKOV_NULL` on both; window_length_bp **150**, context_bp **300**, markov_order **3**, n_generated 120, fit/eval 600/600; corpus fingerprint n_cds_total **4340**, n_windows **4340**, strand **+2118/-2222**, n_overlapping_cds **2708** — ALL CONFIRMED. Two mechanism hypotheses were tested and REFUTED and must not be restated as fact: mode collapse (pairwise identity 0.249 against a ~0.25 random expectation, 6/6 distinct) and coding-like output (stop density 4.83 vs natural 4.67; ATG rate 1.75 vs 1.75). **The EXTERNAL claim** — that the AlphaGenome author's own advice for a non-model organism is to fine-tune on your own data rather than consume predictions as-is — is a primary-source paraphrase from the technical talk `rdBtxtcS4nM`, and its scope is narrower than it looks: the advice was given to a PLANT researcher, and plants are eukaryotes, so the transferable part is the principle and not an endorsement for bacteria. Flagged rather than asserted. Hardware facts verified by measurement: cu118 `torch.cuda.get_arch_list()` includes **sm_61**, so a GTX 1070 (CC 6.1) is covered; consumer Pascal runs fp16 at ~1/64 of fp32, so fp32 is mandatory and `recommended_dtype()` already enforces it.

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** One model, one corpus already on disk, one named bar (below 0.6679), one named method (LoRA at fp32). The blockers are a predecessor measurement and a physical machine — neither is unbounded scope.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  1. **(ADOPTED)** LoRA fine-tune on MG1655 upstream windows; headline is distinguishability below 0.6679 at mode `both`, corpus fingerprint stamped.
  2. Cheaper precursor: continue pre-training on the *E. coli* genome at large (not just promoters), which is more data but less targeted.
  3. Narrower: fine-tune and evaluate at mode `kmer` only (bar 0.5990). Weaker, because the positional features are what broke the k-mer ceiling.
  4. Rejected: fine-tune AND sweep sampling in the same run. Named to exclude — it confounds the two levers, and the sweep is already queued.
  5. Rejected: train on the same windows used as the natural comparison set. Named to exclude — it would make the discriminator score the generator on its own training data.

## Goal Hierarchy
### Long-term (12+ months tier)
A generator that proposes bacterial regulatory sequence indistinguishable from natural sequence by any cheap discriminator, as the proposal half of the generative loop.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Sampling sweep reports | Distinguishability per temperature/top-k setting, so the levers stay separable | blocked (server) |
| 2 | Train/eval split of the promoter corpus | Fine-tuning windows DISJOINT from the discriminator's natural set | +1 day |
| 3 | LoRA fine-tune at fp32 on the 8 GB card | A trained adapter + loss curve | +3 days |
| 4 | Verdict + artifact | Distinguishability vs 0.6679, fingerprint stamped, either way | +1 week |

### Short-term (<=1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | WAIT on the sampling sweep + hardware. Do not start training | stop | Soraya | blocked |
| 2 | (now, CPU) Build the disjoint train/eval window split so fine-tuning cannot train on the comparison set | edit-local-code | Soraya | 1 day |
| 3 | (on unblock) LoRA fine-tune at fp32; never fp16 on Pascal | edit-local-code | server | +3 days |
| 4 | (on unblock) Re-run the falsifier with the fingerprint stamped; commit artifact | run-tests | server | +1 day |

## MVP Criteria

Attempt budget: 3 per criterion. **MVP = an answered question.**

| # | Criterion | Kind | Predicate |
|---|---|---|---|
| 1 | The disjoint fine-tuning split exists (CPU-doable NOW) | file-exists | `dna_decode/glm/finetune_corpus.py` |
| 2 | Its tests pass | test-exit-0 | `uv run pytest tests/test_glm_finetune_corpus.py -q` |
| 3 | A verdict artifact exists | file-exists | `wiki/glm_generator_finetune_verdict.json` |
| 4 | The verdict is recorded here | project-state-row | Action Log row whose outcome names the verdict |

## State Snapshot
### Assumptions
- The generator failure is a ZERO-SHOT failure that fine-tuning can move — **medium**. This is the bet; the alternative is that the architecture or the tokenizer is the limit.
- LoRA at fp32 fits 8 GB for a 1.2B model — **medium** (weights ~4.9 GB, adapters and their optimizer state small, 150-300 bp sequences at small batch). Untested on that card.
- ~4,340 promoter windows is enough to fine-tune on — **low**. It is a small corpus for a 1.2B model, which is why candidate 2 exists.
- Distinguishability is the right target — **medium**. It is necessary, not sufficient: a k-mer discriminator measures LOCAL composition, and published audits name long-range organisation as the real failure.

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | GENERator-1.2B at default sampling: 0.7706 (both) / 0.7237 (kmer) vs Markov null 0.6679 / 0.5990 — WORSE_THAN_MARKOV_NULL | wiki/glm_generator_falsifier_REAL_2026-10-07.json | high | 2026-10-07 |
| 2 | Mode collapse REFUTED: pairwise identity 0.249 vs ~0.25 random, 6/6 distinct | same artifact + result memo | high | 2026-10-07 |
| 3 | Coding-like output REFUTED: stop density 4.83 vs 4.67, ATG 1.75 vs 1.75 | same | high | 2026-10-07 |
| 4 | cu118 arch list includes sm_61, so the 1070 is covered; fp32 mandatory on consumer Pascal | measured via torch.cuda.get_arch_list() | high | 2026-10-07 |
| 5 | The model is natively `llama`; the HF repo ships a custom `tokenizer.py`, so trust_remote_code is load-bearing for the TOKENIZER | cached snapshot inspection | high | 2026-10-07 |
<!-- project-state:end:evidence -->

### Unknowns
- Whether the sampling sweep already closes the gap, which would change this family's premise before it starts.
- Whether 4,340 windows is enough, or whether whole-genome continued pre-training is needed first.
- Whether the mechanism is the tokenizer (custom remote code, nucleotide-level vs BPE) rather than the weights.
- Whether a k-mer+positional discriminator is even the right target, given that long-range organisation is the published failure mode.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| GD-H1 | LoRA fine-tuning brings distinguishability below 0.6679 at mode `both` | open | never |
| GD-H2 | Mode collapse explains the zero-shot failure | **falsified** | 2026-10-07 |
| GD-H3 | Coding-like output explains the zero-shot failure | **falsified** | 2026-10-07 |
| GD-H4 | Sampling settings alone close the gap, making fine-tuning unnecessary | under-investigation | queued on server |
| GD-H5 | The tokenizer, not the weights, is the limiting factor | open | never |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| Fine-tuning windows must be DISJOINT from the discriminator's natural set | 2026-10-07 | Otherwise the generator is scored on its own training data |
| fp32 only on the 1070 | 2026-10-07 | Consumer Pascal runs fp16 at ~1/64 of fp32; recommended_dtype() enforces it |
| Sampling sweep reports BEFORE fine-tuning starts | 2026-10-07 | Two generator levers in one run cannot be attributed |
| Every number carries the corpus fingerprint | 2026-10-07 | NCBI re-annotates under the same accession; two runs can look comparable and not be |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Ratify the drafted idea-anchor | Soraya | user | `features/glm-generator-finetune/idea-anchor.md` |
| Whether to accept whole-genome continued pre-training (candidate 2) if 4,340 windows proves too few | Soraya | GD-H1 result | Materially larger compute; a scope call |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
The generator is measured WORSE than a 3-mer Markov chain at default sampling (0.7706 vs 0.6679), two mechanism hypotheses are refuted, the sampling sweep is queued on an offline machine, and fine-tuning — the author's own prescription — has never been attempted.

### Target state / terminal condition
A committed artifact stating whether LoRA fine-tuning brings distinguishability below 0.6679 at mode `both`, with the corpus fingerprint stamped and the sampling lever already reported separately.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` + `gates-passed`. 4 unknowns, 5 hypotheses at init; **2 already falsified** (GD-H2, GD-H3).

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Build the disjoint train/eval window split (CPU, NOW) | edit-local-code | medium | low | low | 1 day |
| 2 | WAIT on the sampling sweep | stop | none | **high** | n/a | 0 |
| 3 | (on unblock) LoRA fine-tune at fp32 | edit-local-code | **high** | **high** | high | 3 days |
| 4 | (on unblock) re-run falsifier, stamp fingerprint, commit | run-tests | high | high | medium | 1 day |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- When the sampling sweep reports. If a sampling setting already clears 0.6679, GD-H4 is confirmed and this family's premise weakens materially — re-rank before any training.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — REQUIRES per-action human approval

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-10-07 | propose | project-init protocol executed by hand for GLM family G-D | ledger created; 3a PASS (all figures + hardware facts confirmed; 2 mechanism hypotheses recorded as falsified), 3b PASS project, 3c PASS |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Ratify or redirect the drafted idea-anchor** at `features/glm-generator-finetune/idea-anchor.md`.
- **Hardware authority:** this family needs the server PC, which is being physically moved. Action 1 is CPU-doable now; actions 3-4 are not. No decision needed unless you want a paid GPU instead — which would be a money gate and is therefore yours.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-10-07
- **Progress signal:** (none yet — init only)
