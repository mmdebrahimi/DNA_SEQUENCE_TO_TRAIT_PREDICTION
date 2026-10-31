# GLM G-D — Generator fine-tuning (the un-pulled lever)
<!-- project-schema: 0.1 -->

> Initialized 2026-10-07. Project ID: glm-generator-finetune-2026-10-07. Originating goal (verbatim user input): "Fine-tune GENERator-v2-prokaryote-1.2b on real E. coli promoter windows rather than consuming it zero-shot, and measure whether it gets below the 3-mer Markov null."

## Project Context
- **Project ID:** glm-generator-finetune-2026-10-07
- **Project root:** C:\Users\Laptop Portable\pythonprojects\dna_decode
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
| GD-H4 | Sampling settings alone close the gap, making fine-tuning unnecessary | falsified | 2026-10-31 |
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
| 2 | 2026-10-29 | run-tests | **The sampling sweep RAN and REPORTED -- the first half of this family's unblock trigger has an artifact, and the trigger is still NOT discharged.** Kaggle T4, the frozen Step-4 grid (temperature 1.0/0.8/0.6, `--bins 6`, `top_k` unswept), 3.2 + 2.4 + 2.4 min. Artifacts `wiki/glm_sampling_sweep_2026-10-29_*.json`, commit 1733d8ae. | Nulls reproduce the recorded 0.5990 (`kmer`) / 0.6679 (`both`) EXACTLY in both modes at all three points -- `sweep_gates.null_reproduced` ok=true with delta_vs_anchor AND delta_vs_recorded both 0.0, so temperature touches only the GLM's sampling and never the Markov null (measured, not assumed). Anchor self-check CLOSE BUT NOT EXACT: kmer 0.7355 vs recorded 0.7237 (+0.0118), both 0.7787 vs 0.7706 (+0.0081) -- the recorded run's seed is unknown per the plan's own limits, so do NOT report this as a byte-exact reproduction. RESULT (lower distinguishability = better): the kmer gap over the null falls 0.1365 -> 0.0401 -> 0.0733 across T=1.0 -> 0.8 -> 0.6, so **temperature is a real lever worth ~3.4x on the gap, best at T=0.8** -- and **NO point beats the null in either mode** (every cell `WORSE_THAN` or `NO_BETTER_THAN_MARKOV_NULL`; `both` = 0.7787 / 0.7097 / 0.7156, none clearing 0.6679). **So this family's `:5` note does NOT fire: the premise is NOT weakened and no re-rank is due.** TRIGGER NOT DISCHARGED, and the reason is a defect in the kernel rather than in the result: `--verify-knob` occurs ZERO times in `scripts/kaggle/build_glm_sampling_sweep_kernel.py` or the emitted kernel, so `sweep_gates.verdict_eligibility` returns `excluded:requested-but-unverified` for T=0.8 and T=0.6 and only `eligible:default-sampling` for T=1.0 -- 2 of 3 points are not plan-eligible. The knob demonstrably moved the output (0.6391..0.7355 with the null pinned), so this is a missing ATTESTATION, not an inert knob. `wiki/glm_generator_sampling_sweep.json` -- the canonical name the `gd-finetune-gated-on-sampling-sweep-not-compute` expiry test watches -- was **deliberately NOT created**: writing it would fire a LoRA fine-tune's unblock trigger on a sweep this project's own gate excludes. **`blocked_by` left UNTOUCHED per this plan's Step 6** (*do NOT change G-D's `blocked_by` edge for the 8 GB GPU*); an earlier edit of mine that rewrote that cell, the `### Current state` line and the Evidence row was REVERTED to byte-identical-to-HEAD, because the first is forbidden by the plan and the other two hard-route through their own ops under AC9. Whether the hardware half is now discharged -- server1 is a GTX 1070 with exactly 8 GiB and the Kaggle T4 has 16 GiB -- is an AUTHORITY call, surfaced not taken. Frozen five untouched. |
| 3 | 2026-10-31 | run-tests | **The sampling-sweep half of this family's unblock trigger is DISCHARGED, and the measurement did not move.** Re-ran the sweep with `--verify-knob` wired (Kaggle T4, kernel `glm-sampling-sweep-vk`, COMPLETE 14:00:38Z), pulled the grid, and emitted `wiki/glm_generator_sampling_sweep.json` under a code-asserted discharge condition. Artifacts `wiki/glm_sampling_sweep_2026-10-31_vk_*.json`; commit f9c47d7e. | Row 2 recorded the trigger as NOT discharged for a precise reason, and that exact reason is now closed: `sweep_gates.verdict_eligibility` returns **`eligible:verified` at BOTH non-default points** with `knob_contract.verification.verdict = TEMPERATURE_IS_LIVE` (low 0.7 vs high 1.5, identical 0 of n=1), and `eligible:default-sampling` at T=1.0, which requested no knob and therefore has NOTHING to verify rather than being unverified -- **3 of 3 plan-eligible, up from 1 of 3.** **THE NUMBERS ARE IDENTICAL TO 2026-10-29** -- all 6 falsifier rows at all 3 points, diffed field by field on distinguishability / markov_null / gc_candidate / verdict -- so the attestation added a self-check WITHOUT perturbing the result, which is the only honest outcome for a self-check; a moved number would have meant the knob was not what we thought. Still NOT a byte-exact reproduction of the RECORDED 2026-10-07 baseline (anchor kmer 0.7355 vs 0.7237, both 0.7787 vs 0.7706; that run's seed is unknown), exactly as row 2 warned. **THE RESULT REMAINS NEGATIVE AND THIS DISCHARGE CHANGES NOTHING ABOUT IT:** no point beats the Markov null in either mode (every cell `WORSE_THAN` or `NO_BETTER_THAN_MARKOV_NULL`; `both` = 0.7787 / 0.7097 / 0.7156, none clearing 0.6679), temperature is a real lever on the kmer gap over the null (0.1365 -> 0.0401 -> 0.0733, best at T=0.8) but NO setting clears it, and `top_k` was NOT swept. **So this family's `:5` note still does NOT fire: the premise is NOT weakened and no re-rank is due.** The trigger file is emitted BY CODE from the real artifacts via `scripts/emit_glm_sweep_unblock_artifact.py`, never by hand -- by-hand being precisely what CLAUDE.md and the handoff forbid. The emitter ASSERTS the discharge condition and refuses otherwise, with two refusals both pinned by `tests/test_glm_sweep_unblock_emitter.py`: (a) any unverified non-default point refuses -- run against the 10-29 artifacts BEFORE emitting, where it names T=0.8 and T=0.6 and writes nothing; (b) a grid carrying NO verified non-default point refuses **even though every point is eligible** -- the load-bearing guard, because `eligible:default-sampling` is a genuine eligibility and a naive all-eligible check would discharge this trigger off a grid that never swept, re-entering the same nothing-to-verify-vs-unverified collapse `sweep_gates` exists to prevent. **G-D IS NOT UNBLOCKED BY THIS and `blocked_by` is UNTOUCHED per this plan's Step 6.** Do NOT read the edge's `8 gb gpu` clause as unmet hardware: the constraint record's own analysis is that the compute half was satisfied from 2026-07-09 (free Kaggle T4 = 16 GB) and the SWEEP was the binding EPISTEMIC half, because sampling and weights are two generator levers and running them together means neither can be attributed. G-D stays blocked because the user HELD that edge on 2026-10-31; lifting it is an AUTHORITY call, surfaced and not taken. `constraint_expiry` now reports `gd-finetune-gated-on-sampling-sweep-not-compute` as EXPIRED -> RE-OPEN, which is R9 working as designed. Frozen five + `uv.lock` byte-unchanged; active 2026-08-31 prospective lock re-verified `ok=True`, no drift. 32 tests pass across the sweep gate, the kernel attestation and the new emitter. |
| 4 | 2026-10-31 | edit-local-code | **THIS FAMILY IS RE-OPENED AND THE FINE-TUNE IS RUNNING.** User authority, verbatim: *"re-open the 3-day fine-tune now"* — overriding the HOLD the same user recorded earlier that day and the frozen plan's Step 6. Built `scripts/kaggle/build_glm_finetune_kernel.py` + the emitted kernel, froze the bar in `wiki/glm_generator_finetune_acceptance_bar.json`, lifted the umbrella edge (umbrella row 21), and launched on a Kaggle T4. Commits `500acc1b`, `52103409`, `bd61f949`, `3113b9e2`. | **GD-H4 IS FALSIFIED by the sweep** — no temperature setting closes the gap, so sampling alone does NOT make fine-tuning unnecessary, and this family's premise survives the re-evaluation trigger its own ledger defines. GD-H1 remains OPEN; this run is its test. **THE BAR WAS FROZEN BEFORE ANY NUMBER EXISTED:** mode-`both` distinguishability strictly below **0.6679**, which is not a chosen threshold but the 3-mer Markov null's own measured value, with all four outcomes named in advance and an explicit prohibition on re-specifying it. LoRA hyperparameters are frozen too (r=16, alpha=32, lr 1e-4, 3 epochs, fp32), so a miss cannot be tuned into a pass — a different setting is a NEW arm, not a re-read. **THE COMPARATOR IS RE-MEASURED IN-RUN, not read off the 2026-10-07 artifact:** that 0.7706 came from CPU and this is CUDA, which the project has already measured yields different sampled sequences for this generator; the artifact also stamps neither its shuffle seed nor `limit_natural` (both recoverable from the launcher — no `--seed`, so default 7 — which is weaker than stamped). So zero-shot and fine-tuned run in the same session, same device, dtype, eval windows, seed and CLI, with the weights as the only difference. **DISJOINTNESS VERIFIED ON THE REAL SUBSTRATE, IN THE REAL ENVIRONMENT:** the Kaggle diagnostic reached `step2-complete` with 4340 windows / 600 fit / 600 eval / 3140 unused / **3737 train / 0 overlap** and fingerprint +2118/-2222, 2708 overlapping — matching the committed manifest exactly. **COST CORRECTION: the "+3 days" estimate was a GTX-1070-at-fp32 figure** (consumer Pascal runs fp16 at ~1/64 of fp32); the corpus is ~560 kb and the T4 replaces the ~12.4 s/sequence CPU rate, so this is a sub-session job. That revision is an estimate from corpus size and device, not a measurement; the run reports its own wall time. **THREE PUSHES WERE SPENT ON MY OWN DEFECTS, all now guarded:** a `.py`-only payload omitted the baseline artifact the split script needs to check the fingerprint (its rc=3 refusal was CORRECT); the builder validated its own syntax but not the emitted kernel's, so a collapsed escape shipped an unterminated string literal twice (now `ast.parse(body)` before the write, proven non-vacuous); and the split subprocess had no PYTHONPATH, which step 1's import check could not catch because the PARENT imported fine. **I also misread "no log" as evidence of infrastructure failure** — Kaggle returns zero output files for ANY ERROR kernel, so that absence distinguishes nothing; both kernels now exit 0 and report their own traceback, the real one to a DISTINCT `glm_finetune_failure.json` that is never a verdict. Frozen five + `uv.lock` untouched. |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Ratify or redirect the drafted idea-anchor** at `features/glm-generator-finetune/idea-anchor.md`.
- **Hardware authority:** this family needs the server PC, which is being physically moved. Action 1 is CPU-doable now; actions 3-4 are not. No decision needed unless you want a paid GPU instead — which would be a money gate and is therefore yours.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-10-07
- **Progress signal:** (none yet — init only)
