# Score The Checkpoint On Both Eval Sets, And Measure What The Orthology Filter Costs

> Run the committed step-6000 checkpoint against the sequence-withheld set and the orthology-free set, with the Markov null and uniform on both, and report the admissibility verdict per set.

## Lens status
Inputs: the recommended next step from this conversation's review of the "Make The Generator Test Discriminate" plan, grounded in `wiki/glm_beyond_kmer_2026-11-22.json`, `wiki/glm_training_corpus_manifest_2026-11-21.json` and `scripts/glm_beyond_kmer.py` read directly.
Degradations: repo-index absent (`repo-index-status.sh` does not exist on this host), so this plan is source-only — every path below was verified by direct `ls`/`grep`, not from an index. `python3` on PATH is the Microsoft Store stub, so the toolkit and every step run under `python`. `codex` and `lmstudio-auditor.sh` are absent, so the prior review of this subject was single-model.

## Problem Statement

The generator's only structural test returned `EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM`: on **6 windows × 8 positions**, the model scored **1.4098 nats** and the strong order-3 null **1.4291**, both against uniform **1.3863** (= ln 4). Two sub-chance predictors cannot be ranked, so the runner correctly refused a verdict.

**Two corrections to the review that produced this plan, both found by reading the code rather than reasoning about it, and both narrowing what this plan may claim.**

**First, the admissibility gate is NOT new.** The review called it "genuinely new, nothing in the repo does this." It is already implemented inline at `scripts/glm_beyond_kmer.py:228-245` — `null_beats_uniform`, `model_beats_uniform`, and the refusing verdict. What does not exist is a **standalone, model-independent gate callable BEFORE scoring and against an arbitrary candidate set**. That is a real but far smaller delta, and this plan extracts rather than invents.

**Second, the runner already records a remedy, and it contradicts the review's recommendation.** Verbatim from `scripts/glm_beyond_kmer.py:246-250`:

> *"filtering orthology by gene necessarily selects strain-specific (prophage / plasmid / HGT) loci, so no eval-side shortcut can be both orthology-free and representative. The representative re-pretrain -- excluding draft E. coli so that ALL MG1655 windows are orthology-free without any gene filtering -- is the only route to an informative measurement."*

The review proposed that the 3,349-window sequence-withheld set might simply replace the 6-window set. That is **refuted**: the withhold is applied *by sequence*, so a test window whose **95%-identical ortholog** sits in the draft *E. coli* the model trained on is NOT withheld. Sequence-withholding does not remove near-ortholog leakage; the gene-level orthology filter does, at the cost of selecting atypical loci.

**So this plan does NOT claim to produce an informative beyond-k-mer verdict.** It measures something narrower and still decision-relevant: **how large the near-ortholog leak actually is**, by scoring both sets and comparing. If the two agree, the leak is small and the orthology filter was costing sample for little; if they diverge, the filter is load-bearing and the recorded re-pretrain is the only route. Either outcome replaces an assumption with a number, and the second is the one the code predicts.

## Codebase Context

**The scorer** is `scripts/glm_beyond_kmer.py` (16,522 bytes). It builds three things and the eval set is the only one this plan changes:
- a **pooled** order-3 chain and an **E. coli-fit** order-3 chain (`chain_ec`, the designated STRONG null), the latter refusing below 500 fitting windows;
- the eval set, at lines 137-148: `only = symbols(MG1655/annotations.gff3) - symbols(DRAFT_GFF)`, then `held = [w for w in wins if w.seq.upper() not in cset]`, then `of = [w for w in held if (m.get(w.gene_id) or "") in only]`, with `REFUSING` (exit 3) when `len(of) < --n-windows`;
- per-window scoring where **the same sampled positions serve every arm and both nulls** — that identity is the basis of the pairing and must not be broken.

Defaults are `--n-windows 120 --positions 30`, so the committed run was invoked at **6 / 8** — the request was lowered to whatever the filter left. `len(of)` is never recorded, which is why Step 2 exists.

**Three candidate sets exist and they are NOT nested the way the review assumed:**
- **A — the designed withhold.** `dna_decode/glm/training_corpus.py:137 falsifier_eval_windows(reference_dir, *, length, limit, seeds, kind) -> set[str]`, a pure deterministic union over seeds, delegating to the scorer's own `build_natural_set` + `fit_eval_halves`. At the committed geometry (`REF` = MG1655, `LENGTH` 150, `SCORING_LIMIT` 1200, `KIND` "promoter", `seeds=range(10)`) it is the **3,349** of the manifest. **It is NOT persisted** — the manifest records only `withheld_union_size`, so it must be recomputed, not loaded.
- **B — `held`.** MG1655 promoter windows absent from the corpus TEXT. A different predicate from A: it also catches windows dropped as exact duplicates or never collected. Size unknown.
- **C — `of`.** B ∩ MG1655-only gene symbols = **6** at the geometry used.

**Verified present:** `data/cache/refseq/GCF_000005845.2` (MG1655 — note this is NOT `data/refseq_cache/`, where it is absent), `data/refseq_cache/GCA_059754755.1` (draft, with `annotations.gff3`), `wiki/glm_training_corpus_2026-11-21.txt`, and the checkpoint `data/glm_diffusion_ckpt/latest.pt` (9,821,617 bytes, 2026-11-22 01:19).

**The checkpoint is the project's single copy of ~3 hours of compute, is gitignored, and `--tag` does not cover it.** A backup exists at `data/glm_diffusion_ckpt_pooled_backup/latest.pt`. Every step here is read-only with respect to it.

### Reusable-Code Survey
- `scripts/glm_beyond_kmer.py:228-245` — the uniform sanity gate, already written. Step 1 extracts this logic verbatim rather than restating it; the extraction must be behaviour-identical on the committed numbers.
- `dna_decode/glm/training_corpus.py::falsifier_eval_windows` — set A, pure and deterministic. Do not write a second withhold computation; its docstring is explicit that the set is "defined by the scorer rather than by a second description of it."
- `dna_decode/glm/corpus.py::build_natural_set`, `fit_eval_halves` — the shared window builders both the scorer and the withhold function already route through.
- `scripts/glm_orthology_free_baseline.py::id_to_symbol`, `symbols` — the gene-symbol maps set C needs; already imported by the scorer.
- `dna_decode/glm/diffusion.py::MASK_ID, DiffusionConfig, build_model, encode` — the checkpoint loading and encoding surface.
- `dna_decode/glm/baselines.py::MarkovGenerator` — both nulls.
- Searched for an existing two-set or set-inventory probe: `grep -rn "orthology" --include=*.py scripts/ dna_decode/ tests/` returns only the scorer and `glm_orthology_free_baseline.py`; no inventory script and no persisted withhold list exist. `graphify-out/GRAPH_REPORT.md` absent.

## Pre-Change Baseline

From `wiki/glm_beyond_kmer_2026-11-22.json`, the run this plan re-measures:

- eval set: 6 windows × 8 positions = **48 predictions**; `set` = "orthology-free held-out MG1655 promoter windows"; `checkpoint_step` 6000.
- primary cell `left_only|markov_ecoli`: diffusion **1.4098**, null **1.4291**, paired delta **+0.0192**, CI **[-0.0107, 0.0494]**, `resolved: false`, diffusion better on **4 of 6** windows.
- `left_only|uniform`: delta **-0.0235**, CI **[-0.038, -0.0047]**, `resolved: true` — the model is resolvably WORSE than guessing.
- `left_only|markov_pooled`: delta **+0.0349**, CI **[0.0079, 0.0731]**, `resolved: true` — but the pooled chain is the WEAK null and "barely beats uniform."
- verdict: **`EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM`**; runtime 1.0 s for the scoring loop.
- corpus side: 63,925 windows / 9,652,674 chars / `tokens_per_param` 11.55 / `withheld_union_size` 3,349.

## Verification Signal

**The load-bearing check is that the extracted gate reproduces the committed refusal.** Fed the baseline numbers above, `eval_admissibility` must return the refusing state — if it does not, the extraction changed behaviour and nothing downstream can be trusted. Proven non-vacuous by feeding numbers that clear uniform and asserting it admits.

Decision-relevant outputs, in order:
1. **|A|, |B|, |C| at the committed geometry** — the first time the orthology filter's cost is a number rather than an inference. |A| must recompute to exactly **3,349** or the geometry has drifted and the run refuses.
2. **The gate verdict per set, evaluated BEFORE any delta is read.** The expectation, stated so it cannot be claimed as a discovery afterwards: set A is expected to ADMIT (it is ordinary core-promoter sequence a corpus-fit 4-mer chain should predict above uniform) and set C is expected to REFUSE (reproducing the baseline). A result where A also refuses would falsify the whole eval-side framing and is the most informative single outcome available.
3. **The paired delta on set A, reported only if the gate admits it**, with its CI and n — against the same `markov_ecoli` strong null.
4. **The A-vs-C divergence**, which is the leak-size estimate this plan exists to produce.

Not a signal: any improvement in the model. The checkpoint does not change, so a difference between runs is attributable to the instrument by construction.

## Implementation Steps

### Step 1: Extract the uniform sanity gate into a standalone model-independent predicate
Files: dna_decode/glm/eval_admissibility.py, tests/test_eval_admissibility.py
Depends on: none

**What changes:**
- New pure module exposing `admissibility(null_ce_nats, model_ce_nats=None, *, uniform=math.log(4)) -> Admissibility` returning one of `ADMISSIBLE` / `REFUSED_NULL_BELOW_UNIFORM` / `REFUSED_MODEL_BELOW_UNIFORM`, plus the margins.
- `model_ce_nats` is **optional and defaults to None**, which is the entire point of extracting: the gate becomes answerable from the NULL alone, before any model is loaded or any checkpoint is read. A candidate eval set can then be screened for free.
- Logic copied from `scripts/glm_beyond_kmer.py:228-245`, not re-derived. No new threshold: `uniform` stays `math.log(4)` and is a keyword argument so a different alphabet size is expressible without editing the module.
- A fourth state `GATE_UNSATISFIABLE_ON_THIS_COMPARTMENT` is deliberately **NOT** added here. It cannot be decided from one set's numbers, and inventing a state nothing can emit would be the vacuous-guard shape this repo has paid for. It belongs to Step 4's reading of the two results.

**Test strategy:**
- Pin the committed baseline: `admissibility(null_ce_nats=1.4291, model_ce_nats=1.4098)` must refuse, reproducing `EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM`'s condition.
- Non-vacuity: a null at 1.30 and model at 1.25 must ADMIT, so the predicate is not constant.
- Null-only call (`model_ce_nats=None`) with a null below uniform must ADMIT, and above uniform must refuse — the pre-scoring path.
- Boundary: a null exactly at `math.log(4)` refuses (strict `<` on the margin), pinned so a future float change cannot silently flip it.

### Step 2: Measure the three candidate eval sets, with no model involved
Files: scripts/glm_eval_set_inventory.py, tests/test_glm_eval_set_inventory.py, wiki/glm_eval_set_inventory_2026-11-23.json
Depends on: none

**What changes:**
- New read-only probe recomputing all three sets at the committed geometry (MG1655, length 150, limit 1200, `kind="promoter"`, `seeds=range(10)`) via the existing `falsifier_eval_windows`, `build_natural_set` and the symbol maps — importing them, never reimplementing.
- Reports `|A|`, `|B|`, `|C|`, the pairwise intersections, and `|C| / |A|` as the orthology filter's measured cost.
- **REFUSES (non-zero exit) if `|A| != 3349`**, because a drifted geometry would silently make every downstream number describe a different cohort than the manifest does.
- Fits the `markov_ecoli` strong null exactly as the scorer does and runs **Step 1's null-only gate on each set**, so each set carries an admissibility verdict with **no checkpoint loaded**. This is the cheap half of the whole plan and it can run before any GPU or model work.
- Writes the artifact; emits no verdict about the generator.

**Test strategy:**
- Offline unit tests on the set algebra with synthetic window lists: C ⊆ B, the symbol-difference filter, and that an empty `only` set yields `|C| == 0` rather than silently falling back to B.
- Assert the probe REFUSES on a stubbed `|A|` of 3348 — the drift guard must fire, proven by making it fire.
- Real-data test marked skip-if-absent on `data/cache/refseq/GCF_000005845.2`, so the suite stays green on a host without the cache.

### Step 3: Score the checkpoint on both sets through one parameterised scorer
Files: scripts/glm_beyond_kmer.py, tests/test_glm_beyond_kmer_eval_sets.py, wiki/glm_beyond_kmer_two_sets_2026-11-23.json
Depends on: Step 1, Step 2

**What changes:**
- Add `--eval-set {orthology_free,sequence_withheld}` (default `orthology_free`, so the existing invocation is byte-unchanged) selecting set C or set A. One scorer, one position-sampling path, one pair of nulls — the alternative of a second script would let the two sets diverge in how they are scored, which is the comparison this plan rests on.
- Replace the inline gate block with a call to Step 1's `admissibility`, so the gate has exactly one definition. The committed verdict string and `verdict_detail` are preserved verbatim for the refusing case.
- Run both sets against the **same** checkpoint (`data/glm_diffusion_ckpt/latest.pt`, step 6000), opened **read-only**; the plan performs no training and writes nothing under `data/glm_diffusion_ckpt/`.
- On set A, request `--n-windows 120 --positions 30` (the scorer's own defaults, i.e. the geometry the committed run could not afford) and record what was actually available.
- Emit one artifact holding both sets' cells, both gate verdicts, and the A-vs-C delta.

**Test strategy:**
- **The extraction must be behaviour-preserving:** re-running the default (`orthology_free`) path must reproduce the committed artifact's primary-cell numbers — diffusion 1.4098, null 1.4291, delta +0.0192, verdict `EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM`. Any drift means the refactor changed the measurement and the run is void.
- Assert the same sampled positions are used across arms and nulls within a window (the pairing identity), by stubbing the samplers and comparing the position lists.
- Assert `--eval-set sequence_withheld` draws from `falsifier_eval_windows` and NOT from the symbol-filtered path, by patching each and checking which is called.
- Assert the checkpoint file's mtime and size are unchanged after a run — a cheap guard on the single irreplaceable artifact.

### Step 4: Read the two results, and declare the re-pretrain arm so it cannot evaporate
Files: wiki/glm_beyond_kmer_two_sets_2026-11-23.json, soraya_runs/arms-ledger.json
Depends on: Step 3

**What changes:**
- Add to the artifact an `honest_limits` block stating, at minimum: this measures the **near-ortholog leak size and the gate verdict per set**, NOT a beyond-k-mer verdict; sequence-withholding does not remove near-ortholog leakage, so an admitting gate on set A does **not** make set A confound-free; and the runner's own recorded remedy (a re-pretrain excluding draft *E. coli*) remains the only route the code claims is informative.
- Record the reading as one of three pre-named outcomes: `LEAK_SMALL_FILTER_OVERPRICED` (A admits, deltas agree within CI) · `LEAK_LARGE_FILTER_LOAD_BEARING` (A admits, deltas diverge beyond CI) · `GATE_UNSATISFIABLE_ON_THIS_COMPARTMENT` (A also refuses — the eval-side framing is dead and the re-pretrain is forced).
- **`declare_arm`** the representative re-pretrain (excluding draft *E. coli*) at **declaration** time under R6, with `cost_tier` large, so a cheap eval-side result cannot retire it as redundant. It is declared, not run.

**Test strategy:**
- Artifact parses as JSON (this repo has shipped an unparseable hand-edited artifact before).
- Every number in the `honest_limits` prose appears in a structured field (the contract-number-audit discipline).
- Assert the arms ledger gains exactly one `planned` arm and that `ArmsLedger.blocking()` then names its criterion — the declaration must actually bind, not merely be written.

## Execution Preview

```
Wave 0 (2 parallel):  Step 1 — extract the gate, Step 2 — measure the three sets
Wave 1 (1 parallel):  Step 3 — parameterise and score both
Wave 2 (1 parallel):  Step 4 — read the results, declare the re-pretrain arm
```

Critical path: Step 1 → Step 3 → Step 4 (3 waves)
Max parallelism: 2 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **This plan cannot deliver the verdict the question wants, and the code says so.** `scripts/glm_beyond_kmer.py:246-250` records that no eval-side shortcut can be both orthology-free and representative. The plan is scoped to measuring the leak instead, and Step 4 forces that limit into the artifact. The failure mode to guard against is reading an admitting gate on set A as licence to publish a beyond-k-mer result.
- Severity: high — **Set A is not confound-free.** Withholding by sequence leaves a 95%-identical ortholog in training. So a good delta on set A is partly recall, and the plan must not call it structure.
- Severity: medium — **Refactoring the gate out of a scorer that produced a committed artifact.** Step 3's first test is the behaviour-preservation pin; without it the extraction could silently change the measurement that motivated the work.
- Severity: medium — **The single checkpoint.** `data/glm_diffusion_ckpt/latest.pt` is ~3 h of compute, gitignored, one copy, and `--tag` does not cover it. Every step here is read-only on it and Step 3 asserts size and mtime unchanged, but no step may be reordered after anything that trains.
- Severity: medium — **`|C|` is unknown and may be smaller than 6.** The committed run requested 6 and the scorer refuses below the request, so 6 is a lower bound on `|C|`, not its value. If `|C|` is exactly 6 the filter is even more punishing than assumed; Step 2 measures it rather than inferring.
- Severity: low — **`held` (set B) has no recorded size and is a different predicate from A.** Reported for completeness; nothing in the plan's conclusions rests on it.
- Severity: low — **Geometry coupling.** `SCORING_LIMIT` 1200 and `seeds=10` define set A. A future scorer at a different geometry needs its own withhold set, as `glm_training_corpus_manifest_2026-11-21.json`'s own honest limits already note. Step 2's `|A| == 3349` refusal is the tripwire.
- Severity: low — **Single-model review.** `codex` and `lmstudio-auditor.sh` are absent on this host, so the review that produced this plan had no independent reviewer and the plan's framing carries that anchoring risk.

## Open Questions

- Whether to spend the representative re-pretrain at all. It is the only route the code calls informative, it costs a fresh ~3 h training run plus a corpus rebuild, and it changes every trained number. Step 4 declares the arm; running it is a user call.
- Whether `tokens_per_param` 11.55 against a corpus the manifest calls "SMALL by construction" means the generation axis should pause for corpus expansion (more species) before any further eval work. That is a scope decision this plan deliberately does not take.
- Whether the committed `EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM` artifact should be superseded by the two-set artifact or kept beside it. The frozen-artifact convention here says keep both; confirming is cheap.

## Architecture Impact

**Expected structural changes:**
- Modules added: `dna_decode/glm/eval_admissibility.py` (a leaf — imports only `math`). Scripts added: `scripts/glm_eval_set_inventory.py`. One existing script modified: `scripts/glm_beyond_kmer.py`.
- Dependency direction: the new module is imported BY the scorer and the probe; it imports nothing from the project, so it cannot create a cycle. Nothing under `dna_decode/` imports from `scripts/`.
- The gate gains exactly one definition where it previously had one inline copy; the number of definitions does not rise.

**Risk areas:**
- `scripts/glm_beyond_kmer.py` produced a committed artifact, so it is the one file where a behaviour change is a measurement change. Pinned by Step 3's reproduction test.
- `dna_decode/glm/training_corpus.py` and `corpus.py` are imported read-only; no step modifies them.

**Sentrux gates:**
- Baseline before Wave 0: yes (`sentrux` resolves at `~/.local/bin/sentrux`).
- Compare after Wave 1: yes.
- Acceptable delta: a minor score change for one added leaf module. A new cycle or a forbidden-layer violation is not acceptable, and neither is any delta implying `training_corpus.py` or `corpus.py` changed.
- Regression-response stance: investigate. A cycle here would mean the leaf module grew a project import, which is a design error to revert.

**Anti-drift checks:**
- `git diff --stat` shows zero changed lines in `dna_decode/glm/training_corpus.py`, `dna_decode/glm/corpus.py` and `dna_decode/glm/diffusion.py`.
- The frozen five remain byte-identical and the prospective lock still verifies (this plan touches none of them, so a change is a defect).
- `data/glm_diffusion_ckpt/latest.pt` size is 9,821,617 bytes before and after.

## Verification

1. Re-measure the targeted test baseline, then run `uv run pytest tests/test_eval_admissibility.py tests/test_glm_eval_set_inventory.py tests/test_glm_beyond_kmer_eval_sets.py -q` and compare the failure SET, not the count, against the standing absent-data set.
2. Confirm the extracted gate reproduces the committed refusal on the baseline numbers, and admits on numbers that clear uniform.
3. Confirm `|A|` recomputes to exactly 3,349, and that the probe refuses when it does not.
4. Confirm the default `orthology_free` path reproduces diffusion 1.4098 / null 1.4291 / delta +0.0192 and the original verdict string.
5. Confirm both artifacts parse as JSON and every prose number appears in a structured field.
6. Confirm the checkpoint's size and mtime are unchanged, and that no file under `data/glm_diffusion_ckpt/` was written.
7. Confirm the arms ledger gained exactly one `planned` arm and that `blocking()` names its criterion.

## Save-time amendments

Captured at: 2026-11-23
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads ONLY `## Implementation Steps` for executable work, so amendments are NOT executable instructions. If an amendment changes a Step contract (file lists, dependencies, structure), re-validation will fire and `/technical-plan` must be re-run before `/execute-plan`.

- scoped to leak-size measurement, not a beyond-kmer verdict
- gate extracted not invented
- re-pretrain declared as an arm

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,severity-high -->
