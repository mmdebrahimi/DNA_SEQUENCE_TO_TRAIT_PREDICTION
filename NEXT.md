## D2 AUTHORISED 2026-10-19 — sourmash is an admissible G8 distance source (user)

The AUTHORITY gate is cleared; the EVIDENCE gate is not, and Step 6 must not fire on the
first alone. `plans/.../Mission_Reconciliation_A_Native_Clonality_Route_Plan` pre-committed
Step 6 (amend the CABBAGE G8 bar's distance source) to **Step 5 exit 0**. Step 5 returned
**exit 2 `INSUFFICIENT_DATA`** — all 10 cohorts blocked, `data/raw/` empty, D: offline — so
**zero agreement numbers exist**. Amending the bar now would record sourmash as validated
against a bar it has never been measured against: bar-shopping in the honest direction,
which this plan already declined once (see the `clause_logic` annotation).

**To fire it:** reattach D: -> `uv run python -m scripts.sourmash_vs_mash_crosscheck`
(needs `uv sync --extra lineage`). Exit 0 = all three clauses met -> Step 6 is unblocked
and pre-authorised, no further decision needed. Exit 1 = sourmash does NOT substitute for
Mash, and that is a real finding, not a failure to retry. Exit 2 = still drive-blocked.

The scorer is now EOL-safe (`verify_derived_from`), so the frozen bar no longer refuses to
load on an LF checkout — previously it would have exited 4 claiming the derived-from table
had changed. Nothing about the bar's three requirements moved.

---

# Open threads — what a fresh session should pick up
## TILE-RUNG SELECTION VALUE + THE HYBRID ARM — both ran 2026-10-14

Two measurements, and the first one **refutes the recommendation I had just given**.

**1. `SELECTION_IS_THE_BOTTLENECK`** (`wiki/glm_tile_selection_2026-10-14.md`). I argued the scoring
half was probably adequate and the prediction criterion should likely be retired for generator work,
citing the designed grid's 0.7937 of achievable selection value at k=10. On the tile rung the best
scorer realises **0.4807 median (2 of 10 seeds clear the frozen 0.50 bar)** — so the 0.7937 does not
transfer and prediction work on this rung IS justified. **It is a near-miss, not a collapse:** the k=10
null band is stdev 0.0485, so 0.48 sits ~10 stdevs above chance, and bar B passes 10/10 (0.485 vs gc
0.295). **The k-trend is the actionable part and it is monotone — 0.388 (k=2) → 0.485 (k=10) → 0.537
(k=20) — so the current oracle DOES clear the bar at k=20.** The answer is conditional on how many
candidates the loop can afford to assay, not a flat verdict on the scorer.

**2. `HYBRID_ADDS`** (`wiki/glm_tile_hybrid_arm_2026-10-14.md`) — **the first learned component to
measurably contribute on this substrate.** conv ⊕ max_tx rank-average: paired **+0.0339 over the
scalar, 10/10 seeds, CI95 [0.0238, 0.0409]** excluding zero; 0.4254 vs 0.3968, moving the rung from
42.2% to 45.2% of its 0.941 ceiling. The published GEMME pattern reproduced on DNA: **conv alone LOSES
to max_tx (−0.0767, 0/10) and still lifts it.** Control clean and decisive — max_tx with a *shuffled*
conv scores −0.1232, 0/10, so a content-free vector hurts and the gain is the conv's content.

**This does NOT overturn the G-A encoder-gate negative.** That gate asked "does the conv beat one GC
number?" — still no, CI includes zero. This asked whether the conv carries anything the *biophysical
scalar* lacks — yes. Both hold.

### The authority fork these two create

**Does a learned component contributing only IN COMBINATION satisfy the self-distillation ban's
precondition** (*"forbidden until a learned representation lands"*)? Not mine to declare. Note the two
results pull in opposite directions on the criterion question: #1 says keep improving prediction, #2
says a learned component now helps. Both support continued work on the *prediction* side, which is the
opposite of the pivot I proposed.

### Named follow-up, deliberately NOT run

**Re-run the selection measurement with the HYBRID as primary scorer.** It is now the best measured arm
(0.4254) and Step 3's bar was frozen with `primary_scorer: max_tx` before any number was read. Asking
whether the new best scorer clears the *same unchanged* bar is legitimate and cheap — but it needs its
own pre-registered run, not a retrofit into a committed artifact, and adding it mid-execution is the
scope drift the plan's non-goals exist to prevent.

Ledger: 10 arms, 6 measured / 4 planned; `ceiling_verdict` still **`blocked:untried-arm`** — a cheap arm
winning does not settle the ceiling question (R5 holding as designed).

## THE MISSION IS NOW DEFINED — `project_state/trait-to-edit/big-idea.md` (authored 2026-10-09)

A validated mission definition for the north star, authored at the user's request from everything
established to date. **It parses and validates against the real `big_idea.py`**, and
`run_minister.preflight(...)` returns `(True, [], [])` — so the `minister` lane is genuinely available
for it, which it was not before (there was no mission big-idea.md; only two stale EP4 ones).

`mission-shape: greenfield` · 4 gaps · 10 checkable endpoint predicates · **Φ = 12** (finite by
construction, so the mission cannot generate unboundedly).

| gap | budget | what closes it |
|---|---|---|
| `amr-label-substrate-validated` | 3 | CABBAGE's G8 + BioProject-level diversity + one leakage-gated scored cell |
| `deconfounder-calibrated` | 2 | the F3 calibration actually run on 1001G + the IBS cross-check |
| `generator-beats-its-baseline-panel` | 4 | a generator beating the mandatory baseline panel on biologically-motivated splits |
| `trait-to-edit-end-to-end` | 3 | one bacterial trait → a ranked edit set, oracle-scored, rails audited |

It also names what must **not** be re-proposed (zero-shot embeddings on natural populations; ESM2 scale;
absolute cross-assay magnitude; the ETD route) so a generative lane cannot re-raise a closed negative.
**Step 4 (climb to animals) is deliberately absent** — premature until the generator gap has a result.

**Launching a minister mission on it is a separate decision and was NOT taken.** Authoring the definition
is not the same act as running it.



Short by design. **Transient** state only: what is in flight, what is waiting on the user, what was just
learned that isn't durable yet. Durable findings belong in `CLAUDE.md` / `wiki/`; scope facts belong in
`scripts/project_status.py` (derived, never written).

Prune aggressively. A stale entry here is worse than an empty file.

_Last updated: 2026-10-09._

> **Pruned 2026-10-05.** This file had sat at 2026-09-28 / 428 lines, and its own "known-stale" section
> had itself gone stale in the opposite direction (it warned that 46 traits was really 44; the live count
> is **46** again after `dna-identify` and `dna-tb` shipped). Everything removed was already durable in
> `CLAUDE.md` (gentamicin v2 lock, HCMV contracting, the declined NNRTI curation, PEAR's reclassification,
> the doubt-layer firing rate) or in `project_state/` (the five project families, F-C/F-D). Checked before
> deleting, not assumed.

---

## IN FLIGHT — deep-research campaign F1 / F3 / F4

Plan `plans/Deep_Research_Missing_Part_2026-09-28_Plan.md`. User ratified full breadth (F1–F4).
Screening spine is built and green (`scripts/deep_research_screen.py`, 13 tests; runs both
`rejection_gates` G1–G10 and `regime.screen_proposal`, and refuses to call anything promising until they
have run).

| family | state |
|---|---|
| **F2** condition coverage | **ANSWERED NEGATIVELY.** "Public K-12" adds samples, not condition breadth — glucose share 60% → **82.7%** of annotated rows, **42%** carry no carbon annotation, only **3** carbon sources reach 10 samples. Computed from the repo's own metadata; the "contact us" wall covers only `log_tpm`. |
| **F1** label substrates | **SEARCHED 2026-10-09 and it produced the campaign's strongest candidate.** `CABBAGE` / EMBL-EBI AMR Portal (release 2026-07) — 1.7M genome-phenotype pairs, lab-measured, CC-BY 4.0, **8.3 MB parquet, MD5-verified on the real FTP surface**. **9 of 10 gates pass**; only **G8 (`min_effective_lineages`)** is open. Through the repo's OWN `tier_for_isolate` it clears N=150 per-class at the STRICT bar on **four cells** — against the 2026-05-18 census verdict that **no** drug clears it. `wiki/f1_cabbage_label_substrate_2026-10-09.md`. Still unsearched: Zenodo/Figshare/Dryad + ETD aggregators; ATLAS is Vivli-request-gated. |
| **F3** de-confounded population designs | **The INSTRUMENT is CALIBRATED ON ONE PANEL PAIR — `VALID_DIFFERENTIAL_RESPONSE`** (2026-10-10, `wiki/f3_deconfound_calibration_1001g_2026-10-10.md`), superseding the 10-09 `blocked:positive_control_not_fold_compatible` whose substrate blocker was RESOLVED by using 1001 Genomes (1,135 natural Arabidopsis accessions + AraPheno phenotypes). It responds differentially to real population structure (**+0.0441** natural vs **0.0000** bi-parental cross). **What is NOT established is TRANSFER** — one positive and one negative panel only; it bounds the INSTRUMENT and validates no downstream organism-level prediction. `eval/regime.py::DECONFOUNDER_CALIBRATION` is the live status — read it rather than this row, which was stale for two days and contradicted this file's own `deconfounder-calibrated` row above. **The literature SEARCH half is still unrun.** Population DESIGN, not organism complexity, remains the measured barrier to the north star's step 4. |
| **F4** cross-assay dosage transfer | **NEVER SEARCHED.** |

**Relaunch SEQUENTIALLY (or 2 at a time) and forbid sub-agent spawning.** The first attempt spawned 12
concurrent agents; all died on a session rate limit with **zero** final reports. First-party search worked
fine immediately after — the fan-out was the constraint, not the session.

## IN FLIGHT — the cross-organism transfer ladder (drafted 2026-10-05, NOT saved)

A 9-step plan to turn the project's **single** step-4 data point (E. coli-tuned essentiality decoder applied
UNCHANGED to human, AUROC 0.5805) into a 5-rung ladder, to separate **(A)** decay-with-distance from
**(B)** a eukaryote-boundary cliff. Drafted + adversarially reviewed in conversation; **never written to
`plans/`** — deliberately, because the review mutated it (the save-once rule).

**Three findings that change it, all measured, none in the plan as drafted:**

1. **AUROC is the WRONG primary metric — it is tie-dominated.** 83.1% of human essentials score EXACTLY 0
   (E. coli 62.1%), and 11 of 13 patterns have P(hit|non-essential)=0.0000 — the decoder is **silent**, not
   wrong. Coverage and precision-where-it-fires move in OPPOSITE directions across the two rungs
   (coverage 0.3789→0.1689; precision 0.3684→**0.9583**) and AUROC blends them. **Use coverage.**
2. **The miss is ~2.5× host-specific vs phrasing, and the ladder survives** —
   `wiki/essentiality_missed_vocabulary_2026-10-05.md`. Mechanism (iii) floor 141/566 (24.9%) vs
   mechanism (ii) floor 57/566 (10.1%). So a cheap pre-test does **not** pre-empt the ladder, but ~10% is a
   catalogue PHRASING gap unrelated to phylogeny and must be separated out before any rung is scored.
3. **The FBA walls do not apply.** `fba/essentiality_labels.py` marks S. aureus `LABEL_WALLED` and
   P. aeruginosa `MODEL_WALLED`, but both are **model-gene-join** walls; this decoder joins to a genome
   annotation, so two "blocked" rungs are reachable. Reachability is still a HYPOTHESIS until fetched.

**Open AUTHORITY forks:** re-run `/technical-plan` around coverage, or `/save-plan` as-is? (The review also
wants a W0 schema probe before any parser, and the functional-class step demoted to a provenance-gated
audit — it is otherwise a second adjudication system.)

## DONE — GLM G-A + G-C v2 build plan (all 9 steps, 2026-10-08)

`plans/GLM_G_A_And_G_C_V2_Build_Plan/technical-plan.md`. Commits `e2aaf54` → `f015cff`.
**Both families answered, both negative, both VALID negatives** (clean nulls + registered protocols,
stamped machine-readably so a broken run could not satisfy the bar). Derivations in the two memos; this
entry is the pointer.

| step | outcome |
|---|---|
| 3 — G-A stop gate | **`NO_POSITIONAL_HEADROOM`.** Position gain **−0.0068** against +0.02. One GC feature (0.3112) beat the 11-, 17- and 26-feature arms. Reconcile anchor **0.3000 exactly**. `wiki/glm_tile_headroom_2026-10-08.json`. |
| 8 — G-C verdict | **`NO_GAIN`.** partial_pooling +0.2033 vs two_model +0.2032 = **+0.0001** against 0.02. `wiki/glm_condition_conditioning_result_2026-10-08.md`. |
| 9 — G-A gate | **`WEAK` by the frozen median rule, NOT distinguishable from zero when paired.** `wiki/glm_encoder_gate_result_2026-10-08.md`. |

**The user resolved the ratification fork as BUILD ANYWAY, and that was the right call even though the
gate's prediction held.** A real conv encoder (13,217 params, widths 6/20/75) reached median +0.3234 vs
gc +0.3112 — margin +0.0122 against a 0.05 bar — but **paired on identical splits it wins only 6/10
seeds, sign-test p=0.754, CI95 [−0.0099, +0.0349] INCLUDING ZERO**; the single-width arm is a literal
coin flip (5/10, p=1.000). **Do not quote the +0.0122 median as the encoder beating the baseline** — the
first two seeds looked like a clean win (+0.024, +0.038) and ten seeds plus a paired test are what
resolved it. The null's own CI *excludes* zero (p=0.002), so the instrument can resolve a real difference
and simply does not resolve one here.

**Why building it was still worth it:** Step 3's proxy was positional GC in 10 bins under a **linear**
ridge, so it could neither see a precise-offset motif nor represent a composition×position interaction.
An encoder escapes both limits and still does not win — which converts a predicted-negative from a weak
proxy into a **measured** one from the model class the proxy stood in for. The encoder is demonstrably
not broken: it learns a planted fixed-offset motif to 0.84 against a computed ceiling of 1.0, its null
collapses 10/10, and both framings agree (AUROC 0.6774 vs 0.6672).

**Scope, do not over-read:** this is THIS encoder at THIS size on THIS substrate. Capacity has already
hurt here (6-mers scored worst at 0.2022), so a larger encoder is a *different* experiment — and it says
nothing about G-B. GA-H1 falsified · GA-H3 confirmed · **GA-H2 under-investigation because two
instruments disagree** (multi-vs-w20 7/10 wins, normal-approx CI barely excludes zero at +0.0007, sign
test p=0.344) · GA-H4 open (the run shows the conclusion is framing-ROBUST, which is not the same as
showing classification is the better framing).

**Two of the G-C arms were degenerate and the tests, not the run, caught it.** A full `seq × medium`
interaction over a BINARY medium spans the same hypothesis space as two per-medium models (corr
0.9999998), and the partial-pooling scale that was supposed to fix that was a **silent no-op** — applied
before `StandardScaler`, which divides out a constant column factor (r 0.9999970, same to 7 decimals).
All three conditioned arms were secretly one model. Fixed by applying the penalty AFTER standardisation
(`alpha/c²`); pinned by a test asserting a pre-scaler scale is inert and a post-scaler one is not.
**Reusable:** a per-feature weight placed before a standardiser does nothing, and the failure is
invisible — every arm still runs and returns a plausible number.

**Do not re-propose "just scale the interaction block" or a wider pooling grid.** The pooling curve is
FLAT (0.0001 spread; argmax wanders 0.03/0.1/1.0 across seeds), so the selected scale is NOISE. Partial
mechanistic account, measured: CENTERING the unit-sum 3-mer block makes the design exactly **rank 63 of
64** (cond 9.3e14), so ridge alpha is already load-bearing and a per-block penalty has little room.

**Cheap and already banked:** 3-mer extraction was **97% of the sweep's cost** (181 s/arm vs 4 s for the
ridge fit) and the features are a pure function of the sequence — identical across arms and seeds. One
cached pass per medium took the 10-seed × 7-arm sweep from ~3.5 h to ~35 min, bit-faithfully
(`two_model` +0.2013/+0.1916 on seeds 0/1, reproduced on 3 runs). Any future arm on this substrate should
reuse `feature_matrix`, not recompute.

## The system design — drafted, awaiting ratification

`plans/Hybrid_Decoder_Architecture_Plan.md`. The hybrid is **not** "catalog + ML predictor" (that framing
scored 0 survivors). Measured shape: **CALL / DOUBT / EVIDENCE**.

| layer | status |
|---|---|
| **L1 CALL** — deterministic curated rules | shipped, 128 cells |
| **L2 DOUBT** — "this call may be incomplete, and why"; never a competing call | **shipped and firing** (`dna_decode/eval/doubt.py`); nobody in the field ships this |
| **L3 EVIDENCE** — de-confounding, nulls, denominators, leakage, provenance | built; **four** cohort-evidence layers now render (lineage, source-concentration, prospective, species-composition) |
| **L4 LEARNED** — forward/inverse, orthogonal modalities; molecular + constructed ONLY | shipped, bounded |

## Waiting on the user — authority calls, not executor tasks

1. **Re-score `Klebsiella × meropenem` on *K. pneumoniae* only?** (new 2026-10-05)
   Measured: the cohort is **38% *K. aerogenes*** and **all 16 false negatives are those off-species
   isolates** — on *K. pneumoniae* there are zero misses (sens 12/12). The pneumoniae-only matrix is
   **sens 1.000 / spec 0.880 / acc 0.919, N=37** against the published 0.467 / 0.900 / 0.683.
   Publishing it edits a **frozen unit** of the reproducibility freeze and changes a SCORED number on the
   trust surface. Current answer is **disclose** (`species_composition` layer, augment-only).
   `wiki/meropenem_fn16_diagnosis_2026-10-05.md`.

2. **Drop lone-porin counting from the meropenem rule?** *(evidence added 2026-10-05, pointing the
   other way — see the closed lever below)*
   All 3 false positives on the AR Bank cohort are porin-only, and the rule's `threshold=1` on a lone
   `ompK35`/`ompK36` truncation over-calls. But it edits sha256-pinned `amr_rules.py` → retires the
   active 2026-08-31 v2 lock and restarts the prospective clock for every cell, motivated by **3
   isolates** on a cohort whose own powering gate HARD-FAILED (n=5 susceptibles).

3. **Whether a single-source cell warrants more than disclosure.** 3 of 10 SCORED AMR cells rest on one
   BioProject. Current answer is *disclose*, namespace-separate. Demoting them is a scope decision.

4. **Do the two valid negatives CLOSE G-A and G-C, or defer them?** (new 2026-10-08; replaces the
   now-RESOLVED ratification fork, which the user answered *build the conv encoder anyway*)
   Both families now have a valid recorded negative, so the remaining question is purely an
   acceptance-bar one and it is the same shape for both. **G-C:** a `NO_GAIN` at n=2 media rules out
   conditioning *as tested on two growth media with 3-mer interactions* — not conditioning in general,
   and specifically not a model whose **representation** is modulated by condition. **G-A:** `WEAK`-but-
   not-distinguishable-from-zero rules out *this* encoder at *this* size on *this* substrate; a larger
   encoder is a different experiment and the result says nothing about G-B. Closing either is a scope
   decision; the executor default is to leave both as recorded negatives and start nothing new.

5. **~~Acquire a natural-population genotype panel so the de-confounder's positive control can run?~~
   ANSWERED BY MEASUREMENT 2026-10-09 — the substrate is resolved and it cost 330 MB.**
   `wiki/f3_substrate_resolved_2026-10-09.md`. The 1001-Genomes imputed SNP matrix is **0.33 GB**
   (HEAD-measured; a literature estimate said "a few GB"), HTTP 200, no auth, and **both MD5s match the
   publisher's own `md5sum.txt`**. It ships a **pre-computed IBS kinship matrix (10.3 MB)**, so the
   decisive question was answerable without touching the 919 MB genotypes. Through the SAME pipeline and
   guards the census ran on `dgrp`: **a fold-compatible grouping EXISTS at p5 — 108 clusters, largest
   9.5%, singletons 30.6%** — where `dgrp` passed at no rung. Phenotypes were **already on disk**: FT10
   overlaps **1,003 accessions (86.3%)**. And FT10 is the right control, not merely an available one —
   it is exactly the trait/panel of the closed Arabidopsis negative, where structure is MEASURED to
   dominate. **What remains is executor work, not an authority call:** run the calibration (needs the
   919 MB matrix for the ridge), and cross-check the grouping against this repo's own
   `ibs_distance_matrix` rather than the project's converted kinship — a disagreement there is a finding
   to diagnose, not a number to choose between. Caveats: genotypes are IMPUTED; p5 is a narrow window
   (p2 and p10 fail in opposite directions); 159 phenotyped accessions have no genotype.

6. **Whether to compress the other long CLAUDE.md bullets — NOW BINDING, not optional (escalated
   2026-10-08).** The growth-ceiling guard fired on an ordinary 60-word finding, because the file had
   arrived at **12 words of headroom**. The ceiling was raised 35,400 → 35,460 on the record (which is
   what the guard's docstring prescribes), but it is now constraining real findings rather than drift, and
   every future addition hits the same wall. `scripts/claude_md_weight.py` reports **43 bullets / 24,892
   words — 70% of the file — PROVABLY stored elsewhere in `wiki/`**, so the diet is available and would
   buy years of headroom. Compressing other sessions' bullets is a scope decision, which is why it is
   still here rather than done. Original note follows. The file loads every session; several long
   bullets cite a resolvable memo, so their derivations could become pointers. Measure first:
   `uv run python scripts/claude_md_weight.py`. Two long bullets have no external store and must stay
   whole — the tool already protects them.

7. **~~Does the CABBAGE label substrate reopen the AMR track?~~ RESOLVED BY THE USER 2026-10-09 — YES,
   it reopens.** *"CABBAGE should reopen the banked AMR track --> i think it should. use your best
   judgement."* The scope of that decision is narrow and was stated when it was made: what CABBAGE
   contradicts is the premise that **no usable label substrate exists** — it does **NOT** reopen the
   embedding negatives (closed on de-confounded signal-vs-structure grounds, which more labels cannot
   address, and independently corroborated at 24k-genome scale by Yu/Barquist 2025) and it does **NOT**
   change the generative north star, whose measured bottleneck is the generator.
   **Now executor work, no further authority needed:** close G8 (`min_effective_lineages`, needs the
   assemblies clustered) and measure BioProject-level source diversity (the `database` column is an
   overlapping aggregator label, and the screen flags its own G2 as "weak").
   `wiki/f1_cabbage_label_substrate_2026-10-09.md`.

8. **How should the GLM umbrella's `ceiling:` criterion be resolved?** (new 2026-10-12)
   Two arms were declared 2026-10-09. The 2026-10-11 run measured a **THIRD** (`encoder-pyramid`, 1.70
   decades) and closed **NEITHER** declared arm, so the criterion is still `blocked:untried-arm`. The
   accounting is now correct and the remaining question is an authority call, not a measurement: one
   declared arm is **SUBSTRATE-GATED** -- >=4 decades of span needs a >=10,000 bp coarsest branch against a
   150 bp tile, so **no budget reaches it** (pinned by arithmetic in `tests/test_glm_encoder.py`) -- which
   means resolving the criterion requires either changing the substrate or retiring the arm, and the
   arms-ledger resolution step was deliberately NOT run. `wiki/glm_pyramid_arm_2026-10-11.md`;
   negative already filed as `neg-5ad65fa5` with its revival condition.
   **Substrate fact for this fork, measured 2026-10-13 and now printed by `scripts/project_status.py`:**
   the "changing the substrate" branch has **cached** candidates — two fragment rungs sit at ~9% used
   (30,000 of 321,123 LB / 318,457 M9 rows, 48–499 bp), while `peak_tiles` and `sigma70_designed` are
   genuinely spent. **This does NOT un-gate the span arm:** 499 bp still cannot host a 10,000 bp coarsest
   branch, so that arm stays substrate-gated and the fork remains an authority call. What changed is only
   that "no data for longer sequences" was false — it came from `size_lever` being absent on three of the
   four rungs. `wiki/substrate_symmetry_2026-10-13.json`; `## The GLM` in `CLAUDE.md`.

## Cheap untried levers (executor work, no authority needed)

- ~~**`dna-identify` cross-check on the 16 meropenem FN.**~~ **DONE 2026-10-05** — ran it, and it
  CORROBORATES: 16/16 FN abstain (`above_max_distance`, 0.1129–0.1339, inside the 0.10849–0.13386 band
  `thresholds.py` records for K. aerogenes), 0 called `klebsiella_pneumoniae`, controls 12/14 called at
  0.0011–0.0090. **Partition agreement with the GenBank labels: 30/30, zero disagreements.** Scope: a
  CLOSED-SET router establishes NOT-pneumoniae, not *aerogenes* specifically.
  `wiki/meropenem_fn_identify_crosscheck_2026-10-05.json`.
- ~~**Re-run AMRFinder under the correct `-O` for the 23 off-species isolates.**~~ **ILL-POSED — CLOSED
  2026-10-05 without running it.** Measured from the pinned image (`amrfinder -l`): the organism list
  offers `Enterobacter_asburiae` / `Enterobacter_cloacae` / `Klebsiella_oxytoca` /
  `Klebsiella_pneumoniae` and **no `Klebsiella_aerogenes`**. So there is no "correct flag" to re-run
  under, and the species contamination in that cohort is **not remediable by re-running** — a structural
  limit of the tool's supported set, not a to-do. (The nearest available options sit under *Enterobacter*,
  which is where this organism's former name sits; choosing one would be a judgement about an unsupported
  species, not a correction.) This is also why the 23 isolates cannot simply be re-scored "properly".
- ~~**The 2 *K. aerogenes* true positives are unexplained.**~~ **DONE 2026-10-05, and they are NOT one
  class.** `GCA_003951185.1` is called R by `ompK35_G41TfsTer32` — a **lone porin truncation**, its only
  counted carbapenem determinant. `GCA_003951605.1` is called R by `blaNMC-A`, a real class-A
  carbapenemase. **This bears on authority fork 2 and points the other way:** a lone-porin call that is
  RIGHT, so dropping lone-porin counting converts a true positive into a false negative. n=1 against the
  3 AR Bank FPs — not settling, but the fix now has a measured cost on the same arm.
- **The 3 Klebsiella cohorts with incomplete cached AMRFinder runs** (gentamicin 3/60, tetracycline 33/60,
  ceftriaxone 54/60) have their species cross-tab withheld for that reason. Completing them needs Docker.

## The FBA switch cell, as it now stands

Direction is fine, **silence** is the problem: 61–76% of genes emit one identical ratio for every
condition. Do **not** build the relative ranking rule (oracle ceiling +4 genes, and a deployable version
must infer k = the original problem restated). **Do not quote a raw exact-set rate** — a conditionally
essential gene is two-sided by construction, so a constant call can never match one. The measured
bottleneck is fixed in both datasets: PRECISE-1K ∩ Keio carbon = **11 of 28**, 621 of 1,035 samples
glucose.

## THE ANCHOR REPRODUCED — the GLM sampling sweep, 3 of 6 steps shipped (updated 2026-10-21)

`plans/Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/`. Steps 1-2 landed: the sampling
knobs are live (`wiki/glm_temperature_knob_2026-10-20.md`), the comparison mode is now a **contract
refusal** rather than two comments nobody enforced (`knob_contract.validate_comparison_mode`), and
`wiki/glm_sweep_invalidation_2026-10-20.md` records that the 2026-10-19 sweep was **doubly invalid**
(inert knob AND `conditional=False` against the aligned promoter set).

**Step 3 then RAN — `E:` was reattached and the anchor reproduced the recorded 2026-10-07 baseline
exactly**, all six rows bit-identical at `positional_bins=6` (`wiki/glm_sweep_anchor_2026-10-20.json`,
**untracked — commit it**): `kmer` 0.7237 vs null 0.5990, `both` 0.7706 vs 0.6679, both
`WORSE_THAN_MARKOV_NULL`. So the harness, the device (CPU) and the corpus build are all vindicated;
what is NOT done is the sweep.

**G-D's unblock trigger is STILL UNDISCHARGED, and it is now ONE run away** — the anchor is the T=1.0
comparator, not a swept point. Three facts worth carrying: the sweep is **CPU inference** (the 8 GB GPU
blocker gates fine-tuning, not this); `E:` attachment is **per-session and never inherited**, so the
4.33 GiB weight file under `E:/hf_cache` can vanish between sessions; and launching must go through
`run_glm_anchor_detached.bat` / `glm_run_lock.py run`, since two concurrent torch jobs is this host's
access-violation hazard and `acquire`+`release` does NOT hold the lock from a launcher (measured — the
second job starts). Pass `--bins 6` explicitly: the runner default is 10 and one anchor attempt was
already voided by it (`..._bins10_UNCOMPARABLE_...json`, which also shows the `both` gap shrinking
0.1027 → 0.0612 from bins 6→10, so prefer `kmer` as the primary comparator). Full resume conditions +
the carried honest limits are in `TODOS.md`.

**Do not quote the 2026-10-19 headline** (GLM 0.8141 vs null 0.3876). Valid baseline is the committed
`wiki/glm_generator_falsifier_REAL_2026-10-07.json`: 0.7237 vs null 0.5990 (`kmer`),
`WORSE_THAN_MARKOV_NULL` — which is the project's **recorded premise**, not a new finding.

## User-authority calls open 2026-10-20

- **Install BLAST+ on this host?** Measured absent (`find_blastn()` -> `None`;
  `wiki/host_toolchain_and_crash_2026-10-20.md`). CLAUDE.md asserted it was installed for months —
  true of the host it was written on, never synced here. Until it is installed, every blastn path takes
  its offline branch: `vf_diff` / virulence overlay, fungal ERG11, HIV RT/PR/IN/CA and SARS-CoV-2 Mpro
  genome-FASTA modes, and the whole typing suite. **Not claimed: that any committed typing validation is
  invalid** — only that none can be RE-RUN here. `winget install --id NCBI.BLAST`, or the NCBI win64
  installer into `C:/Users/Laptop Portable/ncbi-blast/bin` (satisfies the resolver's 3rd branch, no env var).
- **Mash on server2?** Two blockers, neither executor-closable: Docker Desktop will not start from
  session 0 over SSH, and ~5.37 GB free RAM against a ~8 GB WSL2 default. Unblocks via a `.wslconfig`
  cap (`[wsl2] memory=2GB processors=2` — profile scope, so the user's file) or closing BlueStacks
  (~3.3 GB of the pagefile). Reversible route once capped: `wsl --install -d Ubuntu --no-launch` then
  `apt-get install -y mash`; undo with `wsl --unregister Ubuntu`.
- **The lineage metrics stay blocked regardless of both.** `amrfinder_runs/` is absent in 10 of 10
  cohorts and the committed artifacts carry zero per-isolate keys, so the chain ends at AMRFinderPlus ->
  Docker. Sourmash solved distances and does not touch this.

## Known-stale / do not trust without re-deriving

- Any **cell count or trait count written in prose**, including in this file.
  `uv run python scripts/project_status.py` is the authority — it is derived, never written.
- `wiki/project_distillation_2026-08-29.md` says 46 CLI traits. That was wrong when written (it was 44)
  and is now right again by coincidence, which is the cleanest possible illustration of why a written
  count is not evidence.
