# GLM — portfolio umbrella
<!-- project-schema: 0.1 -->

> Initialized 2026-10-07. Project ID: glm-umbrella-2026-10-07. Originating goal (verbatim user input): "name a trait, get the genome edits" — the project north star, with AMR as substrate #1 and never the subject.

## Project Families

The GLM portfolio. **Namespace `G-*` deliberately, because `F-A`..`F-E` belong to the honest-decoder
portfolio in `project_state/dna-decode-2026-05-11.md` and `F1`..`F4` to the deep-research campaign in
`NEXT.md`.** This is a SEPARATE umbrella, not a second section of that one, because that umbrella's Mission
Terminal Condition explicitly excludes a learned predictor and one umbrella carries one terminal condition.

| # | Family | Family slug | Status | blocked_by | Scope |
|---|---|---|---|---|---|
| G-A | Learned representation | glm-learned-representation-2026-10-07 | ACCEPTED | none | learned multi-resolution encoder + trivial head; does it transfer to genomic sequence |
| G-B | Self-distillation | glm-self-distillation-2026-10-07 | ACCEPTED | glm-learned-representation-2026-10-07 | perturb-and-relabel; BLOCKED by a measured degeneracy, not a preference |
| G-C | Condition conditioning | glm-condition-conditioning-2026-10-07 | ACCEPTED | none | one (sequence, medium) model vs two per-medium models |
| G-D | Generator fine-tuning | glm-generator-finetune-2026-10-07 | ACCEPTED | sampling sweep + 8 GB GPU | the un-pulled lever on a measured ZERO-SHOT failure |
| G-E | Distribution alignment | glm-distribution-alignment-2026-10-07 | ACCEPTED | glm-learned-representation-2026-10-07, glm-generator-finetune-2026-10-07 | do the loop's two halves share a distribution |

**Critical path: G-A → G-E.** G-D joins from the side; G-B hangs off G-A; G-C is independent.
**Running now: G-A + G-C (both CPU, substrate on disk, logically independent).**
**Hardware-blocked: G-D.** **Downstream: G-B, G-E** — though G-E's candidate 2 is unblocked today.

## Requirements Flow-down

- **G-A → G-B is HARD and MEASURED.** A ridge teacher's predictions lie in the span of its own features, so
  a linear student reproduces it at agreement **1.0000** (mutation rates 0.02/0.05/0.10; GC-only control
  collapses 0.5963 → 0.1736, proving the probe non-vacuous). There is nothing for perturbation to
  regularise until the student has capacity the teacher lacks. **Do not plan G-B work as startable early.**
- **G-A → G-E.** G-E must score generated sequence with an oracle validated on the matching distribution.
  None exists for genome-like sequence today (GC 0.3000 beats every learned feature there), so G-E cannot
  produce a sound number until G-A succeeds or closes.
- **G-D → G-E.** Characterising the output of a generator that is below a trivial baseline measures the
  distribution of a known-broken sampler.
- **G-D is gated on the sampling sweep REPORTING FIRST.** Sampling and fine-tuning are two generator
  levers; run together, neither can be attributed.
- **G-C must stay independent.** It is the only family whose motivating measurement is already positive, so
  it is the portfolio's hedge against a negative G-A.
- **Nothing downstream may be planned as if G-A will succeed.** Fixed features failing to transfer is
  measured; a learned encoder transferring is a hypothesis with a stated falsifier.

## Mission Terminal Condition

A user names a target molecular phenotype and receives a **ranked set of concrete genome edits**, each
scored by an oracle **validated on the same sequence distribution the edits are drawn from**, with an honest
statement of what the score does not support.

**Not terminal:** an oracle that scores well on a designed grid while the generator emits genome-like
sequence, or vice versa — both halves can look "working" and the loop still be unsound. That is the measured
2026-10-07 configuration and the reason G-E exists.

**Not terminal either:** beating the deterministic catalog. A lookup table cannot generate, so a head-to-head
AUROC measures the wrong quantity (orientation error 0). Comparators are other GENERATIVE methods or a
random-proposal null.

**Portfolio falsifier:** if G-A closes negative AND G-D closes negative, the sequence-model lever space for
this substrate is exhausted and the live question becomes the LABEL/SUBSTRATE question — which routes to the
existing `F-E` label-acquisition family, not to a sixth GLM family. Named in advance so a sixth family
cannot be invented to avoid the conclusion.

## Project Context
- **Project ID:** glm-umbrella-2026-10-07
- **Project root:** C:\Users\Laptop Portable\pythonprojects\dna_decode
- **Captured:** 2026-10-07
- **Originating goal:** name a trait, get the genome edits
- **Horizon (months):** 12
- **Schema:** project-schema 0.1
- **Decomposition:** `plans/GLM_Portfolio_Decompose_2026-10-07.md`

## Empirical Concerns
- **Verdict:** PASS
- **Check status:** attempted
- **Provisional:** NO
- **Findings:** The portfolio-level claims are the conjunction of the five families' verified measurements, each confirmed against its committed artifact at family init (see each family ledger's own Empirical Concerns). No new factual-shape claim is introduced at this level.

## Project vs Research-Program
- **Verdict:** FAIL
- **Provisional:** NO
- **Classification:** research-program
- **Rationale:** "Name a trait, get the genome edits" is unbounded at the top level, and saying so honestly is the point of this gate — it is WHY the work decomposes into five bounded families each carrying its own falsifier. The umbrella is deliberately a research-program; the families are projects. The actionable goal at this level is portfolio sequencing, not a deliverable.

## Refinement Candidates
- **Verdict:** FAIL
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:** the five families G-A..G-E above, each independently bounded with its own success criterion and falsifier. That decomposition IS the refinement.

## Goal Hierarchy
### Long-term (12+ months tier)
A generative decoder that turns a named molecular phenotype into ranked, distribution-matched genome edits.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Oracle question answered | G-A commits a verdict either way | 1 month |
| 2 | Conditioning question answered | G-C commits a verdict either way | 1 month |
| 3 | Generator question answered | G-D commits a verdict once hardware returns | 2 months |
| 4 | Loop soundness established or refuted | G-E reports which distribution generated sequence falls in | 3 months |
| 5 | Portfolio falsifier evaluated | If G-A and G-D both close negative, route to F-E rather than invent G-F | 3 months |

### Short-term (<=1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | Execute G-A (learned representation) | edit-local-code | Soraya | 2 weeks |
| 2 | Execute G-C (condition conditioning) in parallel | edit-local-code | Soraya | 1 week |
| 3 | G-E candidate 2 (reference-vs-reference distance), unblocked today | run-tests | Soraya | 3 days |
| 4 | G-D action 1 (disjoint fine-tuning split), CPU-doable now | edit-local-code | Soraya | 1 day |
| 5 | Collect user ratification of the five drafted idea-anchors | ask-user | user | 1 day |

## State Snapshot
### Assumptions
- The five families are the right cut of the problem — **medium**; derived from a measured triage, but a triage is a judgement.
- Two parallel families is the right concurrency for one attended session — **high**.
- A negative G-A does not strand the portfolio because G-C hedges it — **high**.

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| 1 | The oracle works on grid-like and fails on genome-like sequence | wiki/glm_genomewide_oracle_result_2026-10-07.md | high | 2026-10-07 |
| 2 | The generator is below a 3-mer Markov chain at default sampling | wiki/glm_generator_real_result_2026-10-07.md | high | 2026-10-07 |
| 3 | Condition carries real signal (noise-matched gap -0.0614, n=297,599) | wiki/glm_genomewide_oracle_2026-10-07.json | high | 2026-10-07 |
| 4 | Self-distillation is degenerate on the current oracle (agreement 1.0000) | wiki/glm_selfdistill_probe_2026-10-07.json | high | 2026-10-07 |
| 5 | The AlphaGenome architecture does NOT transfer (1 Mb context / TADs / CTCF are absent in bacteria); the training METHOD does | wiki/alphagenome_lessons_actionable_2026-10-07.md | medium | 2026-10-07 |
<!-- project-state:end:evidence -->

### Unknowns
- Whether any learned representation beats composition on real bacterial promoter sequence.
- Whether the generator's failure is weights, sampling, tokenizer, or architecture.
- Whether the loop's two halves can be brought onto one distribution at all.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| GU-H1 | The GLM bottleneck is the GENERATOR, not the scorer, on grid-like sequence | confirmed | 2026-10-07 |
| GU-H2 | The GLM bottleneck is the SCORER on genome-like sequence | confirmed | 2026-10-07 |
| GU-H3 | Both bottlenecks are closable with sequence-model levers on this substrate | open | never |
| GU-H4 | If GU-H3 fails, the binding constraint is labels/substrate, not architecture | open | never |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| The GLM gets its OWN umbrella, not a second section of the honest-decoder umbrella | 2026-10-07 | One umbrella carries one terminal condition; that one explicitly excludes a learned predictor |
| Namespace is G-* | 2026-10-07 | F-A..F-E and F1..F4 are both taken by live work |
| Two families run in parallel (G-A, G-C); three have real predecessors | 2026-10-07 | Five in parallel would produce five half-finished families |
| Each family's MVP is an ANSWERED QUESTION, not a won bet | 2026-10-07 | A bar only a success could clear makes a recorded negative unreachable |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Ratify the five drafted idea-anchors | Soraya | user | `features/glm-*/idea-anchor.md`; idea-anchor is user-confirmed by standing directive |
| Confirm the parallel/sequential sequencing | Soraya | user | 2 running / 1 hardware-blocked / 2 downstream — my judgement, yours to redirect |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Both halves of the generative loop are measured and both are deficient in different regimes; five bounded families now exist to close or close-out each deficiency, with two startable immediately.

### Target state / terminal condition
See `## Mission Terminal Condition`.

### Progress proxy
- **v0.1 metric:** families whose verdict is committed (0 of 5), plus `unknowns-retired` across the portfolio.

### Candidate next actions
| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | Execute G-A | edit-local-code | **high** | **high** | high | 2 weeks |
| 2 | Execute G-C in parallel | edit-local-code | high | high | medium | 1 week |
| 3 | G-E candidate 2 (unblocked) | run-tests | medium | **high** | low | 3 days |
| 4 | G-D action 1 (CPU split builder) | edit-local-code | low | low | low | 1 day |
| 5 | Collect idea-anchor ratifications | ask-user | low | medium | low | 1 day |
<!-- project-state:end:candidate-actions -->

### Re-evaluation trigger
- After any family commits a verdict, re-rank the portfolio; the flow-down edges change what is eligible.

## Allowed Action Classes (v0.2 placeholder — not enforced in v0.1)
- `propose` / `research` / `write-plan` / `run-tests` / `ask-user` / `stop` — auto; `edit-local-code` — REQUIRES per-action human approval

## Action Log
| # | Date | Action class | Description | Outcome |
|---|---|---|---|---|
| 1 | 2026-10-07 | propose | GLM portfolio decomposed into 5 families; umbrella + 5 family ledgers created | decompose at plans/GLM_Portfolio_Decompose_2026-10-07.md; namespace G-* chosen to avoid two live collisions |
| 2 | 2026-10-09 | propose | Portfolio-level VOI re-assessment after G-A + G-C both reported negative; read `eval/regime.py`, the Phase 5 decomposition, and the disk | **FINDING: TWO UNRECONCILED DECOMPOSITIONS, and the active one does not test the thesis the user corrected toward.** On **2026-10-01** `wiki/glm_phase5_decomposition_2026-10-01.md` decomposed the user's corrected goal — *"a global genome LM that learns trends from abundant sequence and adapts to a new organism from the few sequences we can get"* — into **F1–F9** with critical path `F1 → F3 → (F6 ‖ F2) → F9 → F4 → F5 → F7 → F8`. On **2026-10-07** this umbrella instead seeded **G-A..G-E**: a single-organism, MOLECULAR-endpoint portfolio on bacterial promoter tiles. **F1–F9 were never seeded as ledgers** (verified: `project_state/` holds 15 ledgers, none F-named). So 10-08/10-09 drove G-A and G-C to negatives inside a frame that contains no transfer axis — while the Phase 5 doc's own words are *"the central claim of the Phase 5 thesis has never been measured here"*. **STATUS OF THE TWO FREE CRITICAL-PATH FAMILIES:** **F1 PARTIAL** — the transfer axis DID land (`regime.py` carries `organism_transfer` with levels unmeasured/held_out_organism/held_out_clade, and `wiki/essentiality_transfer_ladder_2026-10-06.md` put the first `held_out_organism` row in), but F1's k-shot protocol, leakage audit and baseline gauntlet are absent. **F3 NOT BUILT AS SPECIFIED — but see the CORRECTION in row 3, which found the technique and four panel pipelines already committed under different vocabulary.** No *calibrated protocol* exists, despite the decomposition naming it **"THE central technical risk"** on the grounds that in a model trained across all life *phylogeny is the dominant axis of variation*. **AND F3 IS FULLY UNBLOCKED — $0, laptop, no GPU:** both of its KNOWN-ANSWER calibration substrates are already on disk (`D:/dna_decode_cache/bloom` = the yeast segregants that WORK at r 0.46–0.80; `data/arabidopsis` = the panel that FAILS, within-group r² −0.13 vs structure-only 0.48), plus three further constructed panels not counted in the plan (`bxd` mouse, `arabmagic`, `dgrp` fly). Corroborating asset found while checking: `wiki/bxd_generality_verdict_2026-08-02.md` already measured a SECOND constructed panel — mouse BXD brain weight r=0.574, 6/12 traits beating a permutation null — establishing that genomic prediction on confound-free crosses **generalises fungi → mammal**, which the G-* frame does not use. **RANKED VOI: (1) F3 phylogeny de-confounding** — critical path, $0, both calibration substrates cached, and without it every downstream global-gLM number is uninterpretable because it would reproduce the Arabidopsis outcome at scale; **(2) finish F1's remaining deliverables**; **(3) G-D fine-tune** — still the un-pulled lever but SINGLE-ORGANISM, so it cannot test the transfer thesis. **AUTHORITY FORK SURFACED, NOT TAKEN:** choosing F3 over G-D decides which decomposition is the real portfolio, which is a scope call for the user. |
| 3 | 2026-10-09 | research | CORRECTION to row 2, found by a background grep that outlived the claim: I searched the PLAN's vocabulary ("phylogeny de-confound") and missed the IMPLEMENTATION's ("lineage deconfound", "gp_arm", "fm_prep") | **Row 2 UNDERCOUNTED the built infrastructure, and one of today's filings was mis-specified as a result.** (a) The de-confounding TECHNIQUE exists and is worked: `scripts/crossaxis_lineage_deconfound.py` runs clade-grouped `GroupKFold` against naive KFold with a collapse-toward-0.5 reading — it is not F3 (it de-confounds one domain question, resistance→virulence in E. coli) but F3 would generalise it rather than invent it. (b) FOUR constructed-panel pipelines are committed, not one: `yeast_bloom_gp_arm.py`, `bxd_gp_arm.py`, `arabmagic_gp_arm.py`, `dgrp_learned_decoder.py`, plus `build_cross_kingdom_summary.py` and a committed `wiki/cross_kingdom_validation_summary.{md,json}`. (c) F1 has a committed falsifier (`scripts/transfer_benchmark_falsifier.py`). **So F3 is genuinely unbuilt AS SPECIFIED (no protocol calibrated on the yeast-works/Arabidopsis-fails pair) but much CHEAPER than row 2 implied — which STRENGTHENS it as the top VOI while correcting the claim.** **(d) THE MATERIAL FINDING, and it corrects a revival condition filed EARLIER TODAY:** `scripts/yeast_bloom_fm_prep.py` states the prior explicitly — *"in a bi-parental cross the markers are a sufficient statistic, so the FM is expected to tie"*. So **"test a learned representation on a CONSTRUCTED eukaryotic panel" is satisfiable UNINFORMATIVELY**: on a bi-parental cross the ±1 marker already encodes allele identity and a sequence model ties BY CONSTRUCTION. Both affected revival conditions in `soraya_runs/negative-results.json` (the Arabidopsis record and the 0-for-5 family) are narrowed to *a constructed panel whose markers are NOT a sufficient statistic* — mutagenesis / multi-parent / novel-allele, NOT a bi-parental cross. **A revival condition that can be met without being informative is not a revival condition.** Ranking unchanged (F3 top), rationale sharpened. |
| 4 | 2026-10-22 | edit-local-code | Closed C1 -- the last bare `load_peak_tiles()` call -- after finding that the AST guard NAMED for that defect read ONE file; both tile-chain guards now DERIVE their file set from the call sites. Fixed the file-as-cache-dir WinError 183 in the OWNER so both consumers inherit it. Recorded sweep death #4 and the drive's physical departure. Re-fetched the ~100 MB GSE144621 substrate locally, which CLOSES the MPRA half of the tile chain's two-dependency wall. | Commits cb4ff8a9 + 24aff506, PUSHED -- origin/main 24aff506, 0 unpushed. CORRECTION: the four commits before these were reported as pushed in an earlier session summary and were NOT; origin/main had been sitting at e15d66a0. Suite 14 failed / 6920 passed / 186 skipped -- the SAME 14 pre-existing data-absent failures (cyp4f2 7, essentiality-ladder 5, pneumo cps 1, README-quickstart 1), no regression; +16 tests, 80 passed across the 3 affected files. Frozen five byte-unchanged; active prospective lock (2026-08-31) re-verified. KEY STATE CHANGE: the tile chain now gets PAST the tile load and fails on ONLY promcalc (`tiles 44106 \| ceiling 0.9410`) -- the measurement the plan's step ordering was built for and which was unshowable until one dependency was actually resolved. Restoration validated by THREE figures reproduced from an independent re-download: 46,713 rows, 44,106 tiles (the pinned EXPECTED_TILES), 0.9410 noise ceiling. No checksum against the original -- it is on the detached volume. The promcalc probe had reported INPUTS_ABSENT while the substrate sat ON DISK (external_root enumerates DRIVE roots; the copy is repo-local) -- a false absence, now fixed; verdict moved to ENV_ABSENT_DEPS_SATISFIABLE_HERE, which still over-claims and is not a green light. `uv pip install --dry-run promotercalculator` resolves to ONE package with sklearn 1.8.0 untouched, confirming amendment 2's refuted isolation premise FIRST-HAND here. ALSO FOUND: all 15 ledgers recorded a project root from a different machine (C:\Users\Farshad\...), so /project-state's Step 4 path gate was blocking EVERY mutation op project-wide -- which is why prior sessions fell back to raw Edits. Corrected to the measured git toplevel; this row is the first append through the sanctioned op since. NOT DONE, surfaced as the user's calls: regenerate promcalc (parked as an open question by the saved plan; ~44k tiles); stage the 4.33 GiB GENERator weights off the detached volume (contradicts the D: rule); server1 (all three of its sessions read offline in ListAgents). Sweep arm unchanged -- wall class external-wall, 4 deaths / 1 completed point. |
| 5 | 2026-10-22 | run-tests | Launched the promcalc stage-B regeneration -- the single remaining blocker on the tile chain -- detached, resume-safe and supervised, after the bounded smoke gate caught the plan's one flagged [unverified] risk in 24 s. Fixed two independent promotercalculator-1.2.4 breaks, gave stage B its first-ever test coverage, and measured the real per-tile cost against the recorded figure. | Commits 15e081c2 + ff57530f + b2e92186, PUSHED. THE VERSION-DRIFT RISK WAS REAL. The plan flagged that the promoter_calculator version used for the 2026-10-12 run is recorded NOWHERE; on 1.2.4 two independent breaks each killed every tile. (1) `from promoter_calculator import Promoter_Calculator` raises ImportError -- 1.2.4 re-exports only the lowercase FUNCTION; the CLASS moved to a submodule. (2) each per-TSS value is now a PromoCalcResults object with __getitem__ but NO .get, so the worker's five d.get() calls each raised AttributeError and -- because decompose_one catches broadly -- turned every tile into a FAILED tile rather than a visibly wrong one. Caught in 24 s / 40 tiles by the smoke gate instead of after ~11 h. THE CLASS IS UNCHANGED, which is the load-bearing part: 1.2.4's own wrapper.py drives it with the same ctor and output key, so these are an import-PATH and container-SHAPE fix, NOT a model swap -- which keeps stage C's frozen PREREGISTERED reconcile anchors the right acceptance test. COST CORRECTION, and it misled my own plan this morning: the wiki records 0.26 s/tile and 785 tiles/min on 7 workers (~74 min); measured here one tile costs 8.7 CORE-SECONDS pinned to one BLAS thread. 785 tiles/min would need ~0.5 core-s/tile, a ~17x gap no thread setting explains -- left UNRESOLVED rather than argued away. The plan had budgeted '~1 h' and called it the one real cost. THREAD PINNING IS WORTH 15 HOURS. The solve is BLAS-parallel, so my single-process 3.0 s/tile was quietly using ~3 cores and 7 unpinned workers oversubscribed 8 cores into context switching. Pinning OMP/OPENBLAS/MKL/NUMEXPR to 1 per worker: 28 -> 65 tiles/min (2.34x), ~26 h -> ~11 h, measured as a true A/B on the real workload. The restart cost nothing -- --resume recovered 495 tiles from the killed run's truncated gzip, exercising the documented drop-the-last-tile rule. REUSABLE: when a worker pool underperforms its single-process timing, suspect that timing was ALREADY parallel; the per-WORKER rate, not the aggregate, is the diagnostic. SELF-CORRECTION: I had called the sklearn-1.9.1 isolation premise REFUTED. Too strong -- in a FRESH env 1.2.4 DOES pull sklearn 1.9.1; the project dry-run stayed at 1.8.0 only because the existing pin satisfied the resolver. The separate venv is right for BOTH reasons. STAGE A reproduced the published substrate exactly (44,106 tiles / 46,713 rows / ceiling 0.941, both of its own pinned asserts passing). SMOKE after the fix: 2,920 rows / 40 tiles = exactly 73.0 sites/tile, 0 NaN, 0 hex35 parse failures, spacer_len 15-19, header == SITE_FIELDS. NEW COVERAGE: tests/test_promcalc_decompose.py, 12 tests, the first this stage has ever had -- its docstring told readers to keep the dna_decode boundary and nothing enforced it. Proven non-vacuous: all four version-break assertions fire against pre-fix HEAD. SUPERVISION: scripts/promcalc_watchdog.ps1 + promcalc_coverage.py, because 4 long runs died unattended this week and resume-safety is worth nothing if nothing re-launches. Three states, a retry CAP so an instantly-failing job cannot hot-loop, and a coverage check so a COMPLETED run cannot be restarted over its own output; undetermined coverage resumes rather than declaring done. Both detectors tested against the live job before being trusted. Three launch defects found and recorded: an INDENTED PowerShell here-string terminator (parse error -- the first watchdog silently never started), an unquoted -File path splitting on the space in 'Laptop Portable', and uv venv's trampoline python.exe making multiprocessing spawn the BASE interpreter so a venv-path liveness grep misses the workers. IN FLIGHT at write time: ~1,533 of 44,106 tiles, ETA ~03:00 on 2026-10-23. NO completion, NO 44,106/44,106 coverage and NO stage-C verdict is claimed. NOT DONE: stage C and its reconcile gate (needs the full table); stage B's --in-dir parameterisation (deliberately deferred -- editing a live module while its Pool is running); the full suite (7 of 8 cores saturated, and competing with the sweep is the recorded 'observing the sweep slowed the sweep' trap). USER AUTHORITY if the reconcile gate fails: report INDETERMINATE_RECONCILE_FAILED; do NOT re-size the frozen bar. Frozen five byte-unchanged; uv.lock and pyproject.toml untouched. |
<!-- project-state:end:action-log -->

## Open Questions for User
- **Ratify the five drafted idea-anchors** (one per family, under `features/glm-*/idea-anchor.md`). Idea-anchor output is user-confirmed by standing directive, so all five are drafted and parked rather than treated as anchored.
- **Confirm or redirect the sequencing:** G-A + G-C in parallel now, G-D on hardware return, G-B and G-E downstream.
- **Confirm the MVP shape:** each family's bar is "the question is answered", so a measured negative reaches MVP. If you want a family to stay open until it WINS, that is an acceptance-bar redefinition and therefore yours.

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-10-07
- **Progress signal:** (none yet — portfolio just created)
