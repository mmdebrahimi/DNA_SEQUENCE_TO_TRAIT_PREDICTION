# F-A · Doubt Layer (L2)
<!-- project-schema: 0.1 -->

> Initialized 2026-08-31. Project ID: doubt-layer-2026-08-31. Originating goal (verbatim user input): "Build the L2 DOUBT layer: generalize and surface the catalog-completeness signals so every decoder call can carry a machine-readable 'this call may be incomplete, and here is why' block that never emits a competing resistance call".

## Project Context
- **Project ID:** doubt-layer-2026-08-31
- **Project root:** C:/Users/Farshad/PythonProjects/dna_decode
- **Captured:** 2026-08-31
- **Originating goal:** Build the L2 DOUBT layer: generalize and surface the catalog-completeness signals so every decoder call can carry a machine-readable "this call may be incomplete, and here is why" block that never emits a competing resistance call
- **Refined goal (if 3c produced one):** Ship a pure, cell-agnostic doubt layer whose signals are functions over already-computed determinant calls, wired into the decoder record behind a guard test that the block can never contain a resistance call, and measured per cell against the position-novelty flag's 0.604 incumbent.
- **Horizon (months):** 3
- **Schema:** project-schema 0.1
- **Family of:** dna-decode-2026-05-11 (umbrella) · plan `plans/Hybrid_Decoder_Architecture_Plan.md`

## Empirical Concerns
- **Verdict:** N-A
- **Check status:** not-applicable
- **Provisional:** NO
- **Findings:** (The goal text is imperative and contains no factual-shape claim about the external world. Its one empirical presupposition — that catalog-completeness signals exist and that completeness is the measured failure mode — is a claim about THIS repository's own committed artifacts, not about the literature, so a web check cannot adjudicate it. It is instead grounded directly in Evidence rows E1-E4 below, each citing a committed artifact. Recording N-A rather than the `--no-web-check` FLAG because the flag would assert unresolved factual uncertainty where the premise is in fact directly verified.)

## Project vs Research-Program
- **Verdict:** PASS
- **Provisional:** NO
- **Classification:** project
- **Rationale:** Horizon 3 months, under the 12-month bar. The deliverable is bounded and named: a module, a record field, a guard test, and a per-cell measurement artifact. No unbounded verbatim signal ("decode all", "understand any") is present. The success criterion is checkable by test exit status and file existence.

## Refinement Candidates
- **Verdict:** PASS
- **Provisional:** NO
- **Refined-from:** originating-goal
- **Candidates:**
  - **C1 (selected)** — A doubt-layer module exists whose signals are pure functions over already-computed determinant calls (no model, no network, no structures), with tests. *Falsifier:* `uv run pytest tests/test_doubt_layer.py` exits non-zero, or the module imports torch / requests / a structure parser.
  - **C2 (selected)** — The decoder record carries a `doubt` block, and a guard test asserts the block can never contain a resistance call. *Falsifier:* a record is produced whose `doubt` block contains an R/S prediction and the suite stays green.
  - **C3 (selected)** — An artifact reports doubt-layer sensitivity PER CELL against the position-novelty flag's measured 0.604 on the EFV blind spot as incumbent. *Falsifier:* the artifact reports a pooled number, or reports no comparator.
  - **C4 (deferred)** — Doubt signals are visible on the trust surface without changing any cell's evidence tier. *Falsifier:* a cell's `tier` differs before vs after registration.
  - **C5 (rejected)** — The doubt layer improves resistance prediction accuracy. REJECTED at definition time: this is the failed-predictor framing (0-for-5, de-confounded) and adopting it would put L2 back in the regime that keeps dying. L2 must never compete with L1.

## Goal Hierarchy
### Long-term (12+ months tier)
Every decoder call in the tool carries an honest, machine-readable statement of where its own catalog is least trustworthy.

### Mid-term (3-12 months)
| # | Milestone | Success Criterion | Horizon |
|---|---|---|---|
| 1 | Completeness signal generalized from the two known gaps | a screen detects the shared shape from cached determinant calls and rediscovers `rmtE1` 36R/0S blind | DONE 2026-08-31 |
| 2 | Signal measured per cell against a named incumbent | artifact reports per-cell sensitivity vs the flag's 0.604; never pooled | 1 month |
| 3 | `doubt` block in the decoder record, guarded | guard test asserts the block can never carry a call | 1 month |
| 4 | Doubt visible on the trust surface, augment-only | no cell's evidence tier changes | 2 months |
| 5 | Doubt layer covers the target-site cells as well as AMR | both signal kinds route through one vocabulary | 3 months |

### Short-term (≤1 month)
| # | Action | Class | Owner | Horizon |
|---|---|---|---|---|
| 1 | DONE 2026-08-31 (Action Log row 3) -- full-index completeness run. The committed artifact records n_genomes_scanned=1818 for all six drugs with NO cap, over the full 200-genome labelled purity denominator. The old "(prior run capped at 220 genomes/drug)" text described a superseded state and read as pending work: it misdirected a 2026-09-27 --advance run into re-verifying a finished run. | run-tests | Soraya | done |
| 2 | Per-cell doubt measurement artifact vs the 0.604 incumbent | edit-local-code | Soraya | days |
| 3 | `dna_decode/eval/doubt.py` — pure signals over computed calls | edit-local-code | Soraya | days |
| 4 | Guard test: the `doubt` block can never contain a call | run-tests | Soraya | days |
| 5 | Register augment-only on the trust surface | edit-local-code | Soraya | weeks |

## State Snapshot
### Assumptions
- The catalog's dominant failure mode is completeness, not accuracy — **high** (measured twice, independently).
- A doubt signal that never competes with L1 stays out of the failed-predictor regime — **high** (this is a design constraint, enforceable by test).
- Cheap deterministic signals are competitive with model-based ones here — **medium** (one measurement: the flag's 0.604 at zero tool cost).
- Per-cell reporting is required; pooled reporting would hide exactly the variance that matters — **high** (repo-wide precedent).

### Evidence
| # | Claim | Source | Confidence | Captured |
|---|---|---|---|---|
| E1 | Gentamicin rule cannot represent `rmt*`; sens 0.523 vs 0.893 | wiki/gentamicin_rmt_disjoint_validation_2026-08-28.md | high | 2026-08-28 |
| E2 | HIV NNRTI catalog misses 53 resistant isolates at uncatalogued positions | wiki/hiv_esm_vs_catalog_2026-07-09.md | high | 2026-07-09 |
| E3 | Position-novelty flag recovers 60.4% of the EFV blind spot, lift 4.69, zero tools | wiki/hiv_blindspot_position_novelty_2026-07-11.json | high | 2026-07-11 |
| E4 | A general completeness screen rediscovers `rmtE1` 36R/0S blind, ranked first | wiki/determinant_completeness_screen_2026-08-31.md | high | 2026-08-31 |
| E5 | Learned predictors in the natural-population regime are 0-for-5 de-confounded | wiki/organism_gp_regime_correction_2026-08-29.md | high | 2026-08-29 |
<!-- project-state:end:evidence -->

### Unknowns
- Whether the completeness screen's `rmt_like` heuristic (>=3 R carriers, 0 S) generalizes beyond the one family it was tuned to recognise.
- RESOLVED 2026-08-31: the doubt layer's false-positive rate. Family-wise (Bonferroni) correction over families screened per drug drops 4 of 5 raw-signature hits; the single survivor is the confirmed `rmtE1` gap. Residual: 1 confirmed gap is a single case, NOT a rate.
- Whether the AMR-side signal (determinant unrepresentable by the rule) and the target-site signal (novel substitution at a catalogued position) belong in one vocabulary or two.
- Whether a doubt block changes user behaviour at all, or is ignored like every other caveat.

### Hypotheses (Active)
| ID | Statement | Status (open/under-investigation/falsified/confirmed) | Last-tested |
|---|---|---|---|
| H1 | The two known gaps share one detectable shape | confirmed | 2026-08-31 |
| H2 | A deterministic doubt signal beats a model-based one on cost-adjusted value | under-investigation | 2026-08-31 |
| H3 | Doubt can be surfaced without changing any cell's evidence tier | open | (untested) |
| H4 | Cells with no known gap produce few doubt flags (low FP rate) | open | (untested) |
<!-- project-state:end:hypotheses -->

### Decisions Made
| Decision | Date | Notes |
|---|---|---|
| L2 never emits a competing call | 2026-08-31 | the constraint that keeps it out of the 0-for-5 regime; enforced by guard test, not convention |
| Probe the deployed rule, never reimplement it | 2026-08-31 | a verbatim one-row probe cannot drift from DRUG_RULE |
| Rank by signature purity, not raw volume | 2026-08-31 | volume ranking buried the known gap 5th beneath correct exclusions |
<!-- project-state:end:decisions-made -->

### Pending Decisions
| Decision | Proposer | Blocker | Notes |
|---|---|---|---|
| Does L2 become a headline product claim? | Soraya | user authority | changes what the tool IS; currently an internal diagnostic |
| One doubt vocabulary or two (AMR vs target-site)? | Soraya | measurement | resolvable by executor once both are measured |
<!-- project-state:end:pending-decisions -->

## Bellman-Inspired Decision Frame

### Current state (one-line summary)
Step 1 shipped and validated (the screen rediscovers the known gap blind); steps 2-4 specified and carrying no authority.

### Target state / terminal condition
C1 + C2 + C3 all met: a pure doubt module with tests, a guarded `doubt` block in the record, and a per-cell measurement artifact naming the 0.604 incumbent.

### Progress proxy
- **v0.1 metric:** `unknowns-retired` count + `gates-passed` count (raw counts, unweighted). 2026-08-31: 3 / 3 MVP criteria met (module+tests / guarded record block / per-cell artifact); 1 unknown retired (the `rmt_like` heuristic's false-positive rate is now measured: 4 of 5 raw hits are noise).
- **v0.2+:** weighted combination of unknowns-retired, gates-passed, evidence-confidence-improved, hypotheses-falsified (TBD via v0.2 design)

### Candidate next actions
Steps 1-5 all completed 2026-08-31 (see Action Log rows 2-7); the table below is the FOLLOW-ON set.
A stale candidate table is not cosmetic — `advance_ranker` reads the FIRST row as the family's next
action, so leaving completed work here made a finished family rank first on the portfolio frontier.

| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 2c | Measure the screen on the OTHER mutant-level cells (sarscov2-mpro is TN-starved 37R/5S; fungal has no free phenotype source) -- currently declared UNMEASURED, which is honest but incomplete | research | low | med | high | days |
| 3 | Re-screen when any NEW independent label set lands (both known gaps needed one) | research | low | high | high | ongoing |
<!-- project-state:end:candidate-actions -->
### Retired candidates (record)

Moved out of the ranked table 2026-09-28: these were completed/answered but still sat INSIDE
the marker region, and `advance_ranker` takes row 1 as `next_action` unconditionally -- so a
post-compaction `--advance` re-picked finished work (measured elsewhere: one candidate redone
five times). Kept verbatim; the verdicts are the record.

| # | Action | Class | Expected progress | Expected info gain | Uncertainty | Cost |
|---|---|---|---|---|---|---|
| 1 | DONE 2026-08-31 -- FP rate measured: family-wise correction drops 4 of 5 raw-signature hits; 1 survivor is the confirmed `rmtE1` gap | run-tests | med | high | resolved | -- |
| 2 | ANSWERED 2026-09-02 -- ONE vocabulary: the purity signature is well-formed on HIV NNRTI (the best case) and fires (V179F, 15 carriers all R, p=8.8e-06 over 638 units) | run-tests | med | high | resolved | -- |
| 2b | DONE 2026-09-02 -- `dna_decode/data/target_site_completeness.py` + a 2nd signal in `target_site_doubt`; V179F now surfaces STRONG doubt on the real CLI. Augment-only VERIFIED by diff (6 cases, 8 fields differ, 0 non-doubt, calls identical) | edit-local-code | med | med | resolved | -- |
### Re-evaluation trigger
- **Default:** re-run `/project-state` after any action class fires (auto-append to Action Log triggers stale-state check)
- **Manual override:** user invokes `/project-state doubt-layer-2026-08-31` at any time
- **Family-specific:** re-evaluate whenever a NEW independent label set arrives — both known gaps were invisible until one did.

## MVP Criteria
- `test-exit-0 uv run pytest tests/test_doubt_layer.py -q`
- `test-exit-0 uv run pytest tests/test_doubt_record_guard.py -q`
- `file-exists wiki/doubt_layer_per_cell_2026-08-31.json`

Attempt budget: 3.

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
| 1 | 2026-08-31 | propose | /project-init invoked | ledger created |
| 2 | 2026-08-31 | edit-local-code | scripts/determinant_completeness_screen.py + 9 tests (step 1) | shipped; rediscovers rmtE1 36R/0S blind, ranked first |
| 3 | 2026-08-31 | run-tests | step 2: full-index completeness run (1818 genomes, 6 drugs) | raw signature fires on 5 families, only 1 is real -- powering needed |
| 4 | 2026-08-31 | edit-local-code | step 2: dna_decode/eval/doubt.py + per-cell artifact | 1 of 1279 families STRONG (rmtE1 p=4.11e-12); enrichment null rejected as wrong |
| 5 | 2026-08-31 | edit-local-code | step 3: doubt block wired at the _target_site_record seam | verified on the real CLI; guard raises rather than emit a call |
| 6 | 2026-08-31 | edit-local-code | step 4: trust_surface.doubt_layer_for, augment-only | badge-with-vs-without diff green + non-vacuity pinned |
| 7 | 2026-08-31 | run-tests | full suite | 4029 passed, 0 failed; frozen surface byte-unchanged; lock re-verified |
| 8 | 2026-08-31 | propose | /project-state --update-last-evaluation | 3 of 3 MVP criteria MET (doubt module + tests; guarded doubt block in the record; per-cell artifact vs the 0.604 incumbent). 1 unknown retired: the rmt_like heuristic's false-positive rate is measured -- 4 of 5 raw-signature hits are noise. Corrects a contradiction: this section said 1-of-3-not-yet-met while Progress proxy said 3/3. |
| 9 | 2026-09-02 | run-tests | target-site denominator probe (HIV NNRTI, best case) | ONE_VOCABULARY: shape holds (1170 S / 638 units / 12 pure) and fires -- V179F p=8.8e-06; V179D+V179E independently named by the 09-01 OLS curation |
| 10 | 2026-09-02 | edit-local-code | wire the target-site completeness screen into the shipped doubt block | V179F surfaces STRONG on the real CLI (call stays S); augment-only verified by diff; fixed doubt_one_line printing signals[0]'s reason under another signal's tier |
| 11 | 2026-09-27 | edit-local-code | Pin dna_decode/data/target_site_completeness.py to the artifact it names (commit 9bdb3ba) | The index the CLI PRINTS carried 5 run-derived numbers and a header saying every field is traceable to the probe artifact, with nothing enforcing it: one test read the artifact, another read the module, neither compared them. Now pinned via each unit's OWN declared artifact path (self-extending, no hand-maintained map). Verified: carriers 15, p 8.832e-06 through the module's 3-sig-fig rounding, n_units_tested 638 == n_units_with_min_carriers, base_s_rate == round(base_susceptible_rate,4). 7 in-memory mutations each caught by the intended test; dropping the cell exposed a hole in my first draft (grouping keyed off what the module cites looped over nothing), closed by discovering probes from disk. |
| 12 | 2026-09-27 | research | CORRECTION to commit 9bdb3ba's ledger claim | That message said retiring the '(prior run capped at 220 genomes/drug)' row stops advance_ranker re-picking it as table[0]. WRONG: the ranker reads the marker-anchored Candidate next actions table, whose row 1 is already marked DONE 2026-08-31, so this family had no re-pick problem. The stale text is Short-term Goal Hierarchy row 1 (a table the ranker does not read, and which has no /project-state mutation op -- the documented v0.5 API gap). Substance unaffected: Action Log row 3 already recorded the 1818-genome 6-drug run, and the committed artifact confirms n_genomes_scanned=1818 with no cap over the full 200-genome labelled purity denominator. Genuinely pending here remains candidate 2c (screen the other mutant-level cells) and 3 (re-screen when a new independent label set lands). |
| 13 | 2026-09-30 | edit-local-code | Wired the SHIPPED-but-unreachable supervised blind-spot complement into the L2 doubt layer as a third signal | dna_decode/data/hiv_supervised_complement.py shipped 2026-07-12 with leave-one-STUDY-out blind-spot AUROC NNRTI 0.81 / PI 0.89 / INSTI 0.89 and verdicts SUPERVISED_RESCUES_BLINDSPOT + GENERAL_RESCUE (8/11 drugs) + DEPLOYABLE_HOLDS_OOD -- and was reachable from NOWHERE: imported only by its own builder, 0 mentions in cell_registry, no CLI route, absent from the doubt layer whose entire job is flagging catalog incompleteness. THIRD 'shipped but unreachable' instance after HCMV contracts and clinvar/hla routing. WHY A THIRD SIGNAL IS JUSTIFIED (structural, not preference): position_novelty fires only at CATALOGUED positions and the completeness screen needs CARRIERS, so both are STRUCTURALLY silent on a novel substitution at an un-catalogued position with few carriers -- exactly the complement's case. Appended never merged; target_site_doubt gained an optional `call` kwarg and the CLI passes call.prediction so an already-R call reports that the blind-spot question does not arise. VERIFIED ON A REAL ISOLATE: K103S alone -- catalog says S (it carries 103N not 103S), measured wet-lab fold log10 +0.833 = ~6.8x GENUINELY RESISTANT, now DOUBT [strong] at risk 0.717 with the CALL UNCHANGED. NON-VACUITY MEASURED after THREE broken attempts: blind spot = 1,109 catalog-S isolates / 52 truly resistant (base rate 0.047); the shipped 0.5 threshold flags 48 at precision 0.729, recall 0.673, 15.55x ENRICHMENT -- so the threshold is GOOD. THREE ERRORS OF MINE, all mid-run: (a) echoed the call VALUE into the evidence and assert_no_call refused it (it checks values not just keys) -- CLI broken until fixed, evidence now records only WHETHER a call exists; (b) documented the token normalisation as LOAD-BEARING when measured both forms score 0.9585 (the module normalises internally) -- kept as defensive redundancy with a test pinning agreement; (c) TWO successive broken scans each yielding a confident wrong design conclusion -- scan 1's `except Exception: continue` swallowed an AttributeError on all 2,168 rows and I concluded the 0.5 threshold never fires and that I had misused the module, BOTH UNFOUNDED; scan 2 used 103N where the catalog needs K103N so everything became blind-spot and a 0.996 AUROC was the full cohort. The fix was never the format -- it was a one-line PARTITION CHECK (both buckets populated), since both failures put everything in one bucket. HONEST LIMITS: IN-DISTRIBUTION not independent; a RANKING complement NOT a rule (binary fold-in tested and REJECTED at -0.006 balacc) so the tier uses the module's own is_flagged/DEFAULT_THRESHOLD; a quiet complement is NOT reassurance (scores 0.180 on 179F, the one confirmed gap, which the purity screen flags at 15/15 p=8.8e-06 -- they fail on DIFFERENT cases); NRTI+CAI not covered, reported not-measured; NO new registry cell (it makes no independent call, so it belongs in the 26 existing HIV cells' doubt block and the evidence-tier surface is unchanged). 8 new tests + 41 existing doubt/CLI tests pass. Memo wiki/supervised_complement_wired_2026-09-30.md |
| 14 | 2026-09-30 | edit-local-code | Fixed a caller-divergence defect in the new supervised-complement doubt signal, found by an EXISTING test in the full suite | The signal took an optional `call` kwarg and treated UNSUPPLIED as 'assume susceptible', which FIRES. The CLI passes call.prediction and stayed correctly quiet on K103N; a library caller passing nothing got a STRONG blind-spot doubt on the SAME genotype -- one the deployed catalog itself calls RESISTANT. One function, two answers, decided by what the caller happened to know. Caught by tests/test_doubt_not_measured_rendering.py::test_the_only_honest_silence_is_assessed_and_quiet, whose docstring already said 'if this ever returns a line the renderer has become noise'. ALL FOUR of my own new tests passed throughout because every one supplies an explicit call, so none exercised the path the rest of the suite does -- a new signal's own tests agreeing with it says little. FIX: target_site_doubt resolves the call from the same deployed catalog the call came from (new `_deployed_call` -> call_hiv_observed), so both paths agree at the source; it READS the deployed call to decide whether its own question applies and never makes one, and the value never enters the evidence (assert_no_call checks VALUES). TWO deliberate properties: an UNRESOLVABLE call is NOT-ASSESSABLE rather than a fire (the fallback must not reinstate the bug in a new place), and a SUPPLIED call outranks a DERIVED one (a caller asserting S on a catalog-R genotype IS a blind-spot claim, so it still fires). Motivating case unaffected: K103S fires STRONG from either path; V179F still fires from the PURITY screen not the complement (0.180), so 'they fail on different cases' holds. 2 new regression guards VERIFIED NON-VACUOUS by reconstructing the pre-fix behaviour: stubbing the resolver makes the divergence test fail, restoring the assume-susceptible branch makes the already-R case fire again. 85 doubt/CLI tests pass; both CLI cases verified live. Memo updated: wiki/supervised_complement_wired_2026-09-30.md |
| 15 | 2026-09-30 | run-tests | Measured the supervised complement's PI + INSTI arms (the wiring memo's named next step) -- the two classes are NOT equal and PI's arm can never fire | METHOD: partition isolates by the DEPLOYED call_hiv_observed, keep the catalog-SUSCEPTIBLE ones (that set IS the blind spot), ask whether the shipped 0.5 threshold separates the truly-resistant (fold >= the builder's own cutoff_fold 3). Genotype extraction lifted from the complement's OWN builder so tokens match the trained features. RESULT (validation datasets): INSTI NON-VACUOUS -- 369 blind-spot / 19 truly R / 17 flagged / precision 0.824 / recall 0.737 / ENRICHMENT 15.99x, holding at 12.34x on the training set. PI THRESHOLD_NEVER_FIRES on BOTH datasets -- 2 truly-resistant isolates out of 614 (base rate 0.003), max risk 0.388 < 0.5. MECHANISM is the INCUMBENT's, not the complement's: PI v0 is POSITION-based over 12 major positions on a 99-aa protease heavily mutated in treatment-experienced isolates, so nearly every resistant isolate trips one and lands in the catalog-R bucket -- catalog-R fraction 0.660 (PI) vs 0.503 (INSTI) vs 0.488 (NNRTI). PI's catalog looks COMPLETE here for the WRONG reason: complete on sensitivity at the cost of specificity, exactly what the shipped PI v0.1 specificity work addressed. INSTI is position-based too but over a 288-aa integrase with a lower over-call rate, which is why it keeps a real blind spot. CONSEQUENCE STATED NOT HIDDEN: do NOT claim the PI complement adds call-time value; NO code change and NO threshold lowering (that is tuning a doubt layer to produce doubt, on 2 positives that cannot support a threshold either way); 'measured and quiet' is the truthful PI state and is what the signal already reports. SECOND FINDING: the NNRTI headline moved 2x between two files differing only by a '.Full' suffix -- 15.54x on *_DataSet.txt vs 7.57x on *_DataSet.Full.txt, because .Full IS the training set (n_train 4222 matches exactly) with a 2x higher blind-spot base rate (0.107 vs 0.047). Neither is out-of-sample; the deployable claim stays the leave-one-STUDY-out AUROC, and the script now names its in-sample AUROC in the KEY ITSELF (blindspot_auroc_IN_SAMPLE_not_performance) because it is 0.94-0.99 against a deployable 0.81-0.89. CONFIRMATION: this independently-written script REPRODUCES the wiring memo's NNRTI row exactly (1108 vs 1109 isolates, 52 R, 48 flagged, 0.729/0.673, 15.54x vs 15.55x). Promoted from scratch to scripts/supervised_complement_nonvacuity.py (--dataset validation|training|both, --self-check; RAISES ScanRefused rather than returning a number on a degenerate partition / any errored row / a zero base rate; exit 2 when gitignored data absent). 9 tests incl. a pin that the PI cell carries a verdict and NO precision/recall/enrichment. Memo wiki/supervised_complement_per_class_nonvacuity_2026-09-30.md |
| 16 | 2026-09-30 | run-tests | GLM go/no-go: the SUPERVISED genotype-token direction CLEARS the same pre-registered bar the ZERO-SHOT arm failed, on the identical blind spot | The 2026-07-09 CLOSED NEGATIVE ('a learned variant scorer does NOT fill the catalog's blind spot') was ESM2-650M ZERO-SHOT at AUROC 0.4485 vs a pre-registered '>=0.65 AND > burden AND null<0.55'. Put both modes on the IDENTICAL isolate set with the IDENTICAL bar: zero-shot 0.4485/0.4750 FAIL; SUPERVISED over position-resolved substitution tokens 0.8142 under LEAVE-ONE-STUDY-OUT (patient-grouped 0.8223, ungrouped 0.8446) = PASS, delta +0.339. PIPELINE ANCHOR: our ESM2 reproduced the committed 0.4485 EXACTLY on the non-Full subset, which is what licenses the comparison; a mismatch exits 3 and publishes nothing. IT CANNOT BE REDISCOVERING THE CATALOG -- the subset is catalog-NEGATIVE by construction, so every feature is a substitution the deployed catalog does not carry. WHAT IT LEARNED CORROBORATES THREE INDEPENDENT IN-REPO FINDINGS: top feature is 103S (coef 2.69) = exactly the isolate hand-verified today (catalog says S, wet-lab ~6.8x resistant), and 179D (2.08) = position 179, flagged separately by the doubt purity screen (V179F p=8.8e-06) and the 2026-09-01 curation OLS pass (V179D/E). Four and three unrelated routes respectively, no shared code path. RECONCILES the June 5-drug bacterial functional-alphabet negative rather than contradicting it: that probe tested whole-genome k-mers vs curated determinants and found the catalog IS the signal (0.93-1.00 vs 0.54-0.72); HIV NNRTI has 1,111/2,168 catalog-negative isolates with 53 truly resistant, so the blind spot is big enough to learn in. RULE: the learned layer's headroom EQUALS the catalog's incompleteness and varies by arm. TWO DEFECTS OF MINE, both caught by my own tests: (a) I PARAPHRASED the pre-registered bar ('esm'->'score') -- now rule_verbatim (byte-equal to source, pinned by test) beside rule_as_applied (generalisation stated explicitly); (b) a string-substitution edit broke the test file's line continuations. HONEST SCOPE: a LINEAR model over one-hot tokens is the WEAKEST supervised family member, so a pass is a FLOOR not a ceiling and says NOTHING about whether attention/context helps; in-distribution to Stanford; ~96% subtype B so subtype generalisation untested; one drug/gene/cutoff; the pass rests on the absolute >=0.65 floor since burden was anti-predictive (0.3946). It does NOT beat the catalog (0.926 full-cohort) -- it works where the catalog is SILENT, which is why it belongs in L2 and not L1. 9 tests incl. an anti-laundering pin that the old zero-shot number must still fail its own bar. Memo wiki/glm_alphabet_headroom_2026-09-30.md |
<!-- project-state:end:action-log -->

## Open Questions for User
- Whether L2 becomes a headline product claim or stays an internal diagnostic. It changes what the tool is, and no executor step can settle it.
- Whether a doubt block should ever be able to request abstention from L1 (currently: no — it may qualify and explain, never overrule).

## Last Evaluation (v0.2 placeholder — not enforced in v0.1)
- **Date:** 2026-08-31
- **Progress signal:** 3 of 3 MVP criteria MET (doubt module + tests; guarded doubt block in the record; per-cell artifact vs the 0.604 incumbent). 1 unknown retired: the rmt_like heuristic's false-positive rate is measured -- 4 of 5 raw-signature hits are noise. Corrects a contradiction: this section said 1-of-3-not-yet-met while Progress proxy said 3/3.
