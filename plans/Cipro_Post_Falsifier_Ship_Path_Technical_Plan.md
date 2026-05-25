# Cipro Post-Falsifier Ship-Path — Technical Plan

> Verdict-conditional implementation plan for handling Codex's bounded-falsifier results once they arrive from the Precision 7780. Locks the response-per-verdict BEFORE the result lands to prevent outcome-biased decisions (per the 2026-05-14 verdict-vs-budget-decision LESSON). All 4 verdict branches must be designed pre-result.
>
> **Supersedes/extends `wiki/cipro_bounded_falsifier_coordination_plan_2026-05-22.md`.** On conflict between the two documents, the coordination plan wins for PRE-result coordination contract; this plan wins for POST-result execution. (Brainstorm B8.)

---

## Problem Statement

Codex is executing `scripts/cipro_bounded_falsifier.py` (or its equivalent — Codex owns runner mechanics) against the trained N=147 cipro NT-XGBoost classifier on the Precision 7780. Results land at `Downloads/cipro_bounded_falsifier_results_2026-05-22.{md,json}` and may resolve to one of 4 verdicts:

| Verdict | Trigger condition (matches runner's `compute_verdict`) |
|---|---|
| **PASS** | Bucket A (ERS control) ≥ 3/4 in top-10; Bucket B (ELX failure) ≥ 2/4 in top-10 AND median rank shift ≥ 100x; Bucket C (all-negative) recovered into top-50 OR all-indeterminate-with-saturation |
| **FAIL** | Bucket A passes but Bucket B does not improve under positive-only Δ ranking |
| **RUNNER_REGRESSION** | Bucket A fails (working strains broke under the method change) — halt + debug |
| **INDETERMINATE_BUCKET_C** (formerly "REVERT") | Buckets A + B pass but Bucket C is neither recovered nor all-indeterminate-with-saturation — Bucket C's status is genuinely ambiguous. Ship as FAIL with scope-limit naming the ambiguity. *(Review M1 + L1: renamed for clarity + alignment with `compute_verdict`'s actual REVERT branch.)* |

Three independent gates also fire alongside the verdict:

1. **Leakage gate (BLOCKING):** `scripts/leakage_check_dup_accession.py` returns `loso_leakage_present` — if True, the LOSO AUROC is inflated by same-genome train/test leakage on `GCA_025200635.1` (= `562.109860` AND `562.111036`); the falsifier verdict is uninterpretable until the cohort is dedup'd + retrained.
2. **Saturation gate:** if `baseline_proba_R ≥ 0.95 AND max_abs_delta_all_genes < 0.01` on > 50% of Bucket B strains, classifier saturation is a load-bearing co-cause for the ELX failure mode (alongside or instead of lineage confound — they're not mutually exclusive; brainstorm B2). Routes to a SHIP-WITH-CALIBRATION-NOTE branch under PASS, OR to dual-cause naming in the scope-limit doc under FAIL.
3. **Coverage gate:** if any subset strain has `MISSING_FROM_CACHE` or `ISM_RETURNED_EMPTY` indeterminate reason, the verdict is partial — record + proceed only if buckets have ≥ 3 valid strains each.
4. **Verdict-reconciliation gate:** Step 0 sub-step 2.5 (added) re-computes the verdict locally from `per_known_locus` raw data — if local verdict ≠ Codex's emitted verdict, escalate (brainstorm B6).

Deliverable: a written verdict-conditional ship plan that maps each combination of (verdict × leakage × saturation × coverage) to a concrete sequence of code/doc/commit changes — executable without further design work when the result lands.

Non-goal: this plan does NOT execute any code paths. It locks the decision tree.

---

## Design Decisions

### D1: Verdict-conditional branches are PRE-COMMITTED before results land

**Decision:** This plan defines exactly what happens for each of the 4 verdicts × 3 gate states. Once results arrive, execution chooses the matching branch with zero additional design work.

**Rationale:** Outcome-biased decision is the failure mode the 2026-05-14 verdict-vs-budget-decision LESSON pinned. If we wait to design until after seeing the result, we'll motivated-reason toward whichever branch matches our prior. Pre-committing the response — including "what counts as PASS evidence" — turns the result into an automatic dispatch.

**Trade-off:** Some branches may turn out unnecessary (e.g., REVERT is only triggered if positive deltas vanish on working strains, which is unlikely given the bucket A baseline shows them in the audit JSON). Planning effort partially wasted. Acceptable because the cost is small (this document) vs the cost of a biased post-result decision.

---

### D2: PASS and FAIL paths both ship v0 — only the spec differs

**Decision:** PASS ships `v0` with `attribution_scope_confidence = HIGH` baseline. FAIL ships `v0` with `attribution_scope_confidence` field populated as `HIGH`/`PARTIAL`/`INDETERMINATE` per the rule in the coordination plan's scope-limit template (Section 7) + a referenced `wiki/cipro_v0_attribution_scope_limit_<DATE>.md` doc. Both are tagged `v0.0-cipro` in git.

**Rationale:** North star = "AI DNA decoder tool, not papers; failure-tolerant iteration." A v0 that's honest about where attribution works and where it doesn't IS a success outcome. Treating FAIL as a non-ship state would let the project stall on the attribution sub-problem when the predictive AUROC (≥ 0.70 LOSO) already meets the v0 Predictive criterion.

**Trade-off:** Two ship configurations to maintain. Mitigated by the `attribution_scope_confidence` field being a single enum on the JSON output — no parallel code paths in `pipeline.py predict`.

---

### D3: Saturation + lineage-confound are NOT mutually exclusive; report both, decide ship action on primary signal

**Decision (revised per brainstorm B2):** Saturation does NOT "supersede" lineage-confound diagnosis. Both can be present simultaneously (a near-clonal block can ALSO saturate the classifier). When the saturation gate fires, the scope-limit doc reports BOTH causes side-by-side. Ship action depends on the PRIMARY signal:
- **Saturation > 50% of Bucket B + verdict is PASS:** route to PASS-WITH-CALIBRATION-NOTE sub-branch (Step P, sub-step 4a). Spec change is gated on calibration fix landing first.
- **Saturation > 50% of Bucket B + verdict is FAIL:** report dual-cause in scope-limit doc; queue BOTH `Cipro_Calibration_Fix_Plan.md` AND a lineage-confound investigation. Don't pick.
- **Saturation ≤ 50% + verdict is FAIL:** lineage-confound framing is primary; saturation is secondary.

**Rationale:** Probability-scale ISM deltas can't discriminate them as separate causes — only as co-occurring patterns. Forcing one to "supersede" was over-simplifying. The scope-limit doc reader gets full evidence + makes the call.

**Trade-off:** Slightly longer scope-limit doc (one extra section). Acceptable — the doc is the artifact future Claude sessions will read; honest both-cause framing prevents premature lock-in.

### D3a: PASS-WITH-CALIBRATION-NOTE sub-branch (added per H1)

**Decision (added per review H1):** When verdict = PASS AND saturation gate fires on > 50% of Bucket B, run the PASS branch (Step P) BUT defer the v0 spec tightening (Step P sub-step 4) until calibration fix lands. The decoder still tags `v0.0-cipro`, but the v0 spec's interpretability criterion stays at the current floor + the release notes flag the calibration follow-up.

**Rationale:** PASS means the method change is safe (doesn't break working strains). Saturation means cohort-wide interpretability might still be unreliable. Shipping is fine; tightening the spec on PASS + saturation would over-promise.

**Trade-off:** Two PASS sub-branches to maintain. Acceptable — the divergence is at exactly one sub-step (spec update).

---

### D4: Mash-cluster N=147 fires ONLY on PASS (not on FAIL or REVERT)

**Decision:** PASS path runs `scripts/mash_cluster_n147.py` (new) → emits `wiki/cipro_mash_clades_n147_<DATE>.json` → re-stratified LOSO per-clade AUROC. FAIL and REVERT paths skip Mash entirely.

**Rationale:** Mash-cluster's value is in stratifying interpretability reliability per clade ("HIGH for diverse clades, PARTIAL for the ELX-family near-clonal block"). Under FAIL, the method change didn't fix the ELX failure → per-clade stratification doesn't help, and Mash compute is wasted. Under PASS, the ranking change works and we have a real interpretability signal to stratify across clades.

**Trade-off:** PASS path is heavier (~30 min Mash + LOSO recompute + per-clade report). Worth it under PASS; explicitly skipped under FAIL.

---

### D5: Plan emits 4 distinct git-ready commit-message templates, one per verdict

**Decision:** Each verdict branch ends with a pre-written commit message template that includes the verdict + leakage status + saturation status as one-liners in the commit message body (per the coordination plan's pre-commit safety check).

**Rationale:** Commit messages are the durable audit trail across sessions. Embedding the verdict + gate states there means a future `git log --grep=PASS` (or FAIL/REVERT) surfaces the entire decision history without grepping wikis.

**Trade-off:** Adds ~10 lines per branch. Cheap relative to the value of grep-able history.

---

### D6: Path A keeps the existing `cipro_bounded_falsifier.py` runner contract

**Decision:** The post-falsifier plan does NOT modify the runner draft. It consumes the JSON output schema as-is (per StrainResult dataclass defined in `scripts/cipro_bounded_falsifier.py`). If Codex's runtime version diverges from the draft schema, this plan's `--audit-merge-json` style consumers break — and Codex must surface that divergence at file-transfer time.

**Rationale:** Schema is the coordination contract — both sides depend on it. Locking the consumer side BEFORE the producer-side runs catches schema drift at transfer time, not at consumer-runtime.

**Trade-off:** If Codex's runner emits a richer/leaner schema, this plan needs a quick schema-mapping update step (under 30 min) before execution. Documented as a known risk in §Risk Flags.

---

## Codebase Context

Code modules this plan reads/touches in execution (not in this document):

| Path | Read or Touch | Why |
|---|---|---|
| `scripts/pipeline.py` | TOUCH | `cmd_predict` adds optional `attribution_scope_confidence` field (FAIL branch only); helper `_classify_attribution_scope(strain_id, mash_clades, falsifier_passes)` added |
| `wiki/decoder_v0_ux_and_success_criterion.md` | TOUCH | Interpretability success criterion text updated per PASS or FAIL branch (PASS: tighten target ≥ 50% top-10 recovery; FAIL: note scope-limit + reference) |
| `dna_decode/eval/phylogeny.py` | READ | Reuse `compute_mash_distances(..., use_docker=True)` for PASS-path Mash-cluster |
| `scripts/mash_cluster_n147.py` | CREATE (PASS only) | Orchestration: sketch + dist + agglomerative cluster + per-clade LOSO AUROC + JSON+MD sidecar |
| `wiki/cipro_v0_attribution_scope_limit_2026-05-22.md` | CREATE (FAIL only) | Fill from coordination-plan §7 template; lists the 4 ELX-family Bucket B strains + 4 negative-Δ Bucket C strains by name |
| `LESSONS_LEARNED.md` | TOUCH | One-line entry per verdict (PASS / FAIL / RUNNER_REGRESSION / REVERT) |
| `project_state/dna-decode-2026-05-11.md` | TOUCH | Action Log row + Pending Decisions row resolved + Bellman frame refresh |
| `wiki/decisions-log.md` | TOUCH | HIGH-salience decision entry (post-falsifier verdict + ship action) |
| `tests/test_pipeline_predict_v0.py` | TOUCH (FAIL only) | 2 new tests pinning `attribution_scope_confidence` field emission |
| `wiki/cipro_bounded_falsifier_results_2026-05-22.{md,json}` | CREATE (commit from Codex) | Transferred from Precision 7780 |
| `wiki/cipro_leakage_check_dup_accession_2026-05-22.{md,json}` | CREATE (commit from Codex) | Transferred from Precision 7780 |
| `README.md` | TOUCH | One sentence under "Current state" pointing to the appropriate verdict artifact |

---

## Implementation Plan

### Step 0: Pre-execution gates (always fire FIRST)

**Action:** Before branching on verdict, evaluate the 3 gates in this exact order. ANY gate firing halts dispatch.

1. **Read leakage check JSON** (`Downloads/cipro_leakage_check_dup_accession_2026-05-22.json` or equivalent landing path).
   - If `loso_leakage_present == True`: HALT. Branch into Step L (Leakage-Recovery Sub-Plan).
   - If JSON missing: BLOCK. Request Codex re-run the leakage check before proceeding.

2. **Read falsifier results JSON** for completeness:
   - Confirm 12 strain results present (4 per bucket).
   - For each bucket, confirm ≥ 3 strains valid (i.e., not `MISSING_FROM_CACHE` or `ISM_RETURNED_EMPTY`). If any bucket falls below 3 valid: HALT + escalate to Codex for re-run with cache populate of the missing strains.
   - Stale-artifact check: `run_date` field in JSON should match today within ±1 day. If older: confirm with user that the older run is intended (brainstorm L3 / B4-extension).

2.5. **Verdict-reconciliation gate (added per brainstorm B6):**
   - Re-compute verdict locally from the JSON's `per_known_locus` data using `compute_verdict()` from `scripts/cipro_bounded_falsifier.py`.
   - If local-computed verdict ≠ Codex's emitted `verdict` field: HALT + escalate. Either Codex's runner has a verdict-criteria divergence, or the JSON was corrupted in transfer.
   - Spot-check: manually grep cohort parquet for both `562.109860` + `562.111036`; confirm count matches leakage-check claim (brainstorm B9).

3. **Compute saturation gate:**
   - Count Bucket B strains with `saturation_flag == True`.
   - If count > 2 (i.e., > 50% of Bucket B): set `saturation_supersedes = True` → FAIL branch will name saturation as cause.
   - Otherwise: `saturation_supersedes = False` → lineage-confound framing remains the FAIL hypothesis.

4. **Copy artifacts** from `Downloads/` into the repo:
   - `wiki/cipro_bounded_falsifier_results_2026-05-22.{md,json}`
   - `wiki/cipro_leakage_check_dup_accession_2026-05-22.{md,json}`
   - Git add + stage; commit at end of Step 0 with message:
     ```
     artifact(cipro-falsifier): import results + leakage check from Precision 7780

     Verdict: <PASS|FAIL|RUNNER_REGRESSION|REVERT>
     loso_leakage_present: <true|false>
     saturation_supersedes: <true|false>
     ```

**Verification:** After Step 0, `git log -1` shows the import commit with verdict + gate states in the body. `wiki/cipro_bounded_falsifier_results_2026-05-22.json` is parseable; `verdict` field matches expected enum.

---

### Step P: PASS branch

Triggered if `verdict == "PASS"`. Sub-branches on `saturation_supersedes` (renamed `saturation_co_occurs` per D3 revision).

**Sub-steps:**

1. **Create Mash-cluster orchestration script** at `scripts/mash_cluster_n147.py`:
   - Args: `--cohort data/processed/stage2_n150_cipro_cohort.parquet --refseq-cache <path> --output-prefix wiki/cipro_mash_clades_n147_<DATE>`.
   - Body: extract per-strain FASTA paths; call `compute_mash_distances(fasta_paths, use_docker=True)`; agglomerative-cluster the distance matrix; persist clade assignments to JSON.
   - **Threshold sweep (added per brainstorm B5):** before committing clade assignments, plot the distance distribution + sweep thresholds in `[0.02, 0.04, 0.05, 0.07, 0.10]`. Pick the threshold that produces ≥ 3 clades AND no single clade > 60% of the cohort AND maximizes intra-clade-vs-inter-clade variance ratio. If no threshold satisfies all 3 → fall back to 0.05 + flag in JSON sidecar.

2. **Run Mash-cluster** on Precision 7780 (Claude on this laptop drafts the script + Codex runs it because Docker Desktop is needed for the Mash image route).

3. **Re-stratify LOSO** per clade: load existing N=147 NT cache + cohort, run `leave_one_*_out_cv` with the new clade groups, report per-clade AUROC.

4. **Update v0 spec** at `wiki/decoder_v0_ux_and_success_criterion.md` — **branch on `saturation_co_occurs`:**
   - **4a. saturation_co_occurs == False:** Tighten interpretability criterion cautiously (review M2 + brainstorm B1). Falsifier certifies 12 strains, not the whole cohort. Update spec to: "≥ 12/12 falsifier-subset strains exhibit working attribution under positive-only Δ ranking; cohort-wide top-10 QRDR recovery rate is reported BUT not gating until measured ≥ 50% on N=67 in a follow-up pass." Add Mash-cluster path to provenance.
   - **4b. saturation_co_occurs == True (D3a sub-branch):** DEFER spec tightening; spec stays at current floor. Release notes flag the calibration follow-up. `plans/Cipro_Calibration_Fix_Plan.md` queued to maintenance mode.

5. **Re-fire `pipeline.py predict`** on 3 sample strains (one ERS, one ELX, one random) to confirm v0 JSON includes the new attribution. Tag the smoke output `wiki/cipro_predict_smoke_postfalsifier_<DATE>.md`.

5.5. **User-confirmation gate (added per brainstorm B3):** before tagging v0, surface a 1-line summary to the user with verdict + gate states + spec-change decision. Wait for explicit "ship" confirmation. User can defer ship even under PASS if other criteria don't land.

6. **Tag git release** `v0.0-cipro` with body referencing `wiki/cipro_bounded_falsifier_results_<DATE>.md`. **Local-only** until user confirms push.

7. **Ledger updates** — Action Log row, Decisions Made row, LESSONS_LEARNED bullet, Bellman frame refresh.

**Verification:** `pipeline.py predict` returns `attribution_scope_confidence = HIGH` for strains in PASS-passing clades; LOSO per-clade AUROC ≥ 0.70 for ≥ 75% of clades; git tag `v0.0-cipro` points at the post-step-6 commit.

**Commit message template:**
```
ship(cipro-v0): tag v0.0-cipro post-falsifier PASS verdict

Falsifier verdict: PASS
loso_leakage_present: false
saturation_supersedes: false
n_clades_with_AUROC_ge_0.70: <X/Y>
top10_QRDR_recovery_rate: <Z>%

Method change applied: positive-only Δ ranking on N=67 cipro-R strains.
Mash-cluster N=147 stratifies attribution_scope_confidence per clade.
v0 spec interpretability criterion tightened.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
```

---

### Step F: FAIL branch

Triggered if `verdict == "FAIL"`. Sub-branches on `saturation_supersedes`.

**Sub-steps:**

1. **Fill scope-limit doc** at `wiki/cipro_v0_attribution_scope_limit_<DATE>.md`:
   - Use coordination plan §7 template verbatim.
   - Populate `<FILL>` placeholders:
     - LOSO AUROC: from `data/processed/models/ciprofloxacin_nucleotide_transformer.pkl` provenance (or recompute if missing).
     - X of 4 Bucket A recovery: from falsifier JSON's Bucket A results.
     - Median rank under failed-falsifier ranking: from Bucket B `best_known_locus_rank_pos_delta`.
   - If `saturation_supersedes == True`: insert "Cause analysis" subsection naming classifier saturation; reference `baseline_proba_R` + `max_abs_delta_all_genes` evidence per Bucket B strain. Queue `plans/Cipro_Calibration_Fix_Plan.md` to maintenance mode.
   - If `saturation_supersedes == False`: keep "batch/clade-associated failure mode" framing; cause unresolved.

2. **Add `attribution_scope_confidence` field** to `scripts/pipeline.py` (review H2 — clarified scope):
   - New helper `_classify_attribution_scope(strain_id, falsifier_subset, audit_locus_tag_prefix=None) -> str` returning `"HIGH"|"PARTIAL"|"INDETERMINATE"`.
   - Without Mash-cluster (FAIL path), the field uses the cheaper **locus-tag-prefix proxy** from the existing audit JSON's `known_loci_hits` (ERS / ELX / ELY / ELV / ELU / ELT / OLZ / AM groupings) — already a meaningful clade proxy without Mash. Rules:
     - `"INDETERMINATE"` if `saturation_flag == True` for the strain (computed on the fly via baseline proba check), OR strain matches the 4 Bucket C all-negative-Δ exemplars.
     - `"PARTIAL"` if strain's locus-tag prefix matches ELX-family (failure case at falsifier-tested fidelity).
     - `"HIGH"` if strain's locus-tag prefix is ERS (control case where attribution works) OR matches another prefix family not in the falsifier's failure-case set.
   - Emit in JSON output + markdown sidecar.
   - **Caveat documented inline:** the locus-tag-prefix proxy is a Bakta/RefSeq batch heuristic, NOT a phylogeny. FAIL-path field has lower fidelity than PASS-path Mash-clade-based field. Documented in the v0 spec under "Honest output."

3. **Pre-commit existing 24 tests, then add 2 new ones** (review M4):
   - `tests/test_pipeline_predict_v0.py` is currently untracked (un-committed). Commit the existing 24 tests as a standalone commit BEFORE the FAIL-branch edits.
   - Then add to the SAME file:
     - `test_predict_emits_attribution_scope_confidence_field` — pins field present + valid enum.
     - `test_predict_attribution_scope_indeterminate_when_saturated` — synthetic high-proba fixture.
     - `test_predict_attribution_scope_partial_when_elx_family` — synthetic fixture with locus-tag prefix `ELX_*`.

4. **Update v0 spec** at `wiki/decoder_v0_ux_and_success_criterion.md`:
   - Interpretability criterion: change from current to "ships at PARTIAL tier with documented scope-limit; `attribution_scope_confidence` field present + tested."
   - Add reference to scope-limit doc in §5 ("How v0 satisfies the criteria").

5. **Update README** with 1-2 sentences under Current State pointing at the scope-limit doc.

6. **Tag git release** `v0.0-cipro` (FAIL-with-scope-limit is still v0; release notes name the scope-limit doc).

7. **Ledger updates** — same fields as PASS, different verdict text.

7.5. **Append "Falsifier resolution" coda to `reports/cipro_interp_audit_analysis_2026-05-22.md`** (review M3). The memo's hypotheses (lineage confound / saturation / etc.) get resolved or refined by the falsifier. Coda is a section labeled `## Falsifier Resolution (<DATE>)` with: actual verdict, gate states, which memo hypothesis was confirmed / refined / falsified.

**Verification:** `uv run pytest tests/test_pipeline_predict_v0.py -v` passes including 2 new tests; `pipeline.py predict` JSON output includes `attribution_scope_confidence`; scope-limit doc populated with real numbers (no `<FILL>` placeholders left); git tag `v0.0-cipro` points at the post-step-7 commit.

**Commit message template:**
```
ship(cipro-v0): tag v0.0-cipro post-falsifier FAIL with documented scope-limit

Falsifier verdict: FAIL
loso_leakage_present: false
saturation_supersedes: <true|false>
n_bucket_B_recovered: <0-1>/4
cause_named: <classifier_saturation|batch_clade_failure_mode_unresolved>

Scope-limit doc: wiki/cipro_v0_attribution_scope_limit_<DATE>.md
attribution_scope_confidence field added to pipeline.py predict + 2 tests.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
```

---

### Step R: RUNNER_REGRESSION branch

Triggered if `verdict == "RUNNER_REGRESSION"` (Bucket A fails).

**Sub-steps:**

1. **Halt ship work entirely.** Do NOT tag v0.

2. **Diff Codex's runner** (whichever Codex used in production) against `scripts/cipro_bounded_falsifier.py` (Claude draft). Look for divergence in:
   - Ranking logic (positive-only vs abs-Δ-with-positive-filter — these differ on ties).
   - Baseline classifier-load idiom (model object structure).
   - Cache `bulk_get` arg shape (tuples vs lists).
   - Per-strain ISM exception handling (swallow vs re-raise — same trap as the smoke-runner LESSON 2026-05-17).

3. **Re-run falsifier** with corrected runner against the same 12-strain subset. If still RUNNER_REGRESSION on second pass: escalate — the audit baseline itself may have a bug.

4. **Ledger update:** Action Log row only (no LESSONS entry until root cause is identified). Bellman frame stays in "bounded falsifier in flight" state.

**Verification:** Second-pass falsifier returns PASS or FAIL (not RUNNER_REGRESSION); diff between Codex runner and Claude draft documented at `reports/cipro_runner_diff_<DATE>.md`.

**Commit message template:**
```
halt(cipro-falsifier): RUNNER_REGRESSION on first pass — diffing runners

Falsifier verdict: RUNNER_REGRESSION
loso_leakage_present: <true|false>
bucket_A_failure_count: <0-4>

Halting ship work. Diff Codex runner vs Claude draft at
reports/cipro_runner_diff_<DATE>.md before re-running.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
```

---

### Step V: INDETERMINATE_BUCKET_C branch (formerly "REVERT")

Triggered if `verdict == "INDETERMINATE_BUCKET_C"` per the runner's `compute_verdict` — Buckets A + B pass, but Bucket C (all-negative-Δ strains) is neither recovered nor all-indeterminate-with-saturation (review M1 + L1).

**Sub-steps:**

1. **Treat as a special case of FAIL** for ship purposes: run Step F sub-steps 1, 2, 3, 4, 5, 5.5, 6, 7 — but the scope-limit doc names "Bucket C ambiguous: 4 strains with all-negative known-locus Δ AND no saturation diagnostic agreement" in cause analysis.

2. **Per-Bucket-C-strain analysis:** for each of the 4 strains, report whether `baseline_proba_R` is high (suggesting saturation) or low (suggesting genuine non-signal). The scope-limit doc surfaces this as a maintenance-mode follow-up.

3. **Queue follow-up:** maintenance-mode plan for understanding the all-negative-Δ pattern — `plans/Cipro_All_Negative_Delta_Investigation_Plan.md`.

**Verification:** Same as Step F; cause section in scope-limit doc explicitly says "Bucket C status indeterminate at falsifier-tested fidelity; investigation queued."

**Commit message template:**
```
ship(cipro-v0): tag v0.0-cipro post-falsifier INDETERMINATE_BUCKET_C

Falsifier verdict: INDETERMINATE_BUCKET_C
loso_leakage_present: false
bucket_C_status: ambiguous (not all-recovered, not all-saturated)

Shipped v0 with current attribution + scope-limit doc.
Bucket C investigation queued.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
```

---

### Step L: Leakage-Recovery Sub-Plan (fires if loso_leakage_present == True)

Triggered if Step 0 gate 1 returns `loso_leakage_present == True`. Falsifier verdict is uninterpretable until cohort is dedup'd + retrained.

**Sub-steps:**

1. **Dedup cohort:** keep one of `(562.109860, 562.111036)`. Default: keep the lower strain_id (`562.109860`) to match alphabetical convention. Update `data/processed/stage2_n150_cipro_cohort.parquet` → new file `stage2_n147_cipro_cohort_dedup_<DATE>.parquet`.

2. **Retrain XGBoost classifier** at `data/processed/models/ciprofloxacin_nucleotide_transformer_dedup_<DATE>.pkl` using `scripts/pipeline.py train`. Document new LOSO AUROC.

3. **Re-run falsifier** with retrained classifier + dedup'd cohort against the SAME 12-strain subset (still in the cohort; `562.109860` survives).

4. **Branch on the NEW verdict** per Step P/F/R/V.

5. **Update cohort-build script** `scripts/build_stage2_n150_cohort.py` to add `assembly_accession` uniqueness assertion BEFORE persisting the parquet — prevents this leak class on future cohort builds (per the 2026-05-22 LESSON on duplicate-accession leakage by construction).

6. **Add regression test** at `tests/test_cohort_build_dedup.py` pinning that `build_cohort` rejects a candidate set with duplicate `assembly_accession`.

**Verification:** Old AUROC vs new AUROC delta is small (typically < 1 pp on N=147, because LOSO same-genome leakage upper bound is ~0.5/N pp rank-based, NOT 2/N pp probability-mass — review H3). Soft sanity check, not theoretical bound; cohort-build test passes; falsifier results on retrained model interpretable.

**Commit message template:**
```
fix(cohort): dedup GCA_025200635.1 + retrain + assert accession-uniqueness

Trigger: leakage check returned loso_leakage_present=true
Dropped duplicate strain: 562.111036 (kept 562.109860)
Old LOSO AUROC: <X>; new LOSO AUROC: <Y>; delta: <Z> pp

Added assembly_accession uniqueness assertion in build_stage2_n150_cohort.py
+ tests/test_cohort_build_dedup.py regression guard.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
```

---

## Execution Preview

| Wave | Steps | Parallelizable | Cost |
|---|---|---|---|
| 0 | Step 0 (pre-execution gates) | no | < 10 min |
| 1 | Branch dispatch on verdict | n/a | n/a |
| 2a | Step P sub-steps 1-2 (PASS path: Mash-cluster) | Mash runs on Precision 7780 | ~30 min |
| 2b | Step F sub-steps 1-2 (FAIL path: scope-limit doc + predict field) | doc + code edits independent | ~45 min |
| 2c | Step R (RUNNER_REGRESSION: halt + diff) | n/a | ~60 min for diff + investigation |
| 2d | Step V (REVERT: same as Step F sub-set) | doc + code edits | ~45 min |
| 2L | Step L (Leakage recovery, if gate fired) | retrain blocks falsifier rerun | ~2-4 hr (retrain wallclock) |
| 3 | Tests + verification (PASS or FAIL only) | parallel test runs | < 15 min |
| 4 | Tag v0.0-cipro + ledger + commit | sequential | < 20 min |

Worst-case path: Step L (leakage gate fires) → retrain → falsifier rerun → Step F (FAIL with scope-limit). Total ~4-5 hr wallclock dominated by retrain.

Best-case path: Step 0 → Step P (PASS) → tag v0. Total ~1-1.5 hr wallclock dominated by Mash run.

---

## Risk Flags

- **Schema drift between Codex's runner output and the Claude draft's StrainResult dataclass.** If Codex emits a richer or leaner JSON than `scripts/cipro_bounded_falsifier.py`, Step 0 sub-step 2 (completeness check) breaks. Mitigation: Step 0 reads the JSON defensively (`.get(...)` with explicit None handling). Escalation: schema-mapping update before resuming.
- **Mash Docker image unavailable on Precision 7780.** PASS path needs `quay.io/biocontainers/mash:2.3--hb105d93_10`. Mitigation: verify image is pulled BEFORE running `mash_cluster_n147.py`; pre-pull as part of Step 0 gate if Codex hasn't run Mash since the install artifact 2026-05-15.
- **Retrain wallclock under Step L may be longer on Precision 7780 with full N=147 cache than the GTX 860M.** Acceptable cost; cache is already populated.
- **Saturation gate threshold (50% of Bucket B) is heuristic.** With only 4 strains in Bucket B, the threshold is 3 strains for "supersedes" — coarse. Mitigation: report all 4 strains' saturation flags in the FAIL scope-limit doc regardless of threshold; the doc reader can re-evaluate. Threshold only changes the cause-naming verb, not the ship action.
- **`v0.0-cipro` git tag commits to a tag name across all 4 verdict branches.** Reverting a tag requires force-push to the tag, which the project's commit-history pattern hasn't required before. Mitigation: do NOT push the tag to origin until the verdict is locked + ledger fully reflects it. Tag is local-only until the post-step-N commit.
- **Plan assumes Codex's runner produces results in 24 hr.** If Codex hits a runtime regression or compute issue, the falsifier may not land. Mitigation: this plan's branches are stable + persist; re-trigger when results arrive.

---

## Verification (after execution, regardless of branch)

1. `git log --oneline 1e2ad95..HEAD` shows: Step 0 import commit + ≥ 1 branch-specific commit + ledger commit.
2. `git log --grep="Falsifier verdict"` returns ≥ 1 hit with the actual verdict.
3. `uv run pytest tests/ -v` passes — including any FAIL-branch-added tests.
4. `wiki/cipro_bounded_falsifier_results_<DATE>.json` exists in repo + readable.
5. Ledger: Action Log has a row; Pending Decisions row for "Bounded falsifier verdict" marked RESOLVED with the verdict.
6. If PASS: `wiki/cipro_mash_clades_n147_<DATE>.json` exists; per-clade AUROC table in MD sidecar.
7. If FAIL or REVERT: `wiki/cipro_v0_attribution_scope_limit_<DATE>.md` exists; no `<FILL>` placeholders.
8. If Leakage path fired: `data/processed/stage2_n147_cipro_cohort_dedup_<DATE>.parquet` exists; `scripts/build_stage2_n150_cohort.py` has accession-uniqueness assertion; regression test passes.
9. PASS or FAIL: `git tag --list 'v0.0-cipro'` returns a tag.
10. README "Current state" reflects shipped status.

---

## What this plan deliberately does NOT cover

- **Phase 2 strategic decisions** (4th-drug substrate, per-gene NT windows, multimodal). Those need their own `/idea-anchor` + `/project-init` cycle per Pending Decisions row 2026-05-17.
- **Cef v0.1 + tet pivot.** Cef v0.1 ships in a future cycle after a cef-specific cohort + audit; tet remains scope-out per the cross-drug architectural finding.
- **Maintenance-mode follow-ups queued from this plan:** `plans/Cipro_Calibration_Fix_Plan.md` (if saturation_supersedes) + `plans/Cipro_Alternative_Attribution_Method_Plan.md` (if REVERT). Both deferred until v0 ships.
- **Pushing the v0.0-cipro git tag to origin.** Tag is local-only until user confirms ship + ledger fully reflects it.
- **External publication / blog / arXiv.** Per north star — AI DNA decoder tool, not papers.

---

## Post-execution status (2026-05-24)

**Outcome:** Codex on Precision 7780 ran the falsifier; verdict = **FAIL** (Bucket A passed, Bucket B did not improve under positive-only Δ ranking). Per Step F, v0 shipped with documented scope-limit as `reports/cipro_v0_scope_limit_decision_2026-05-23.md`. Release packet at `reports/dna_decoder_v0_release_candidate_2026-05-24.md`. The v0 spec was RELOCKED at `wiki/decoder_v0_ux_and_success_criterion.md` to match the implemented cached-strain predictor (not the original genome-input decoder concept).

Codex also caught the duplicate-accession leakage independently + retrained with `leave_one_accession_out` CV (matches this plan's Step L); leakage-safe CV AUROC = 0.8697.

**Step P (PASS branch) did NOT fire** — Mash-cluster N=147 didn't run as part of the post-falsifier dispatch. HOWEVER, `scripts/mash_cluster_n147.py` is now reclassified as **v0.1 infrastructure** (useful for either v0.1 path: real-genome-input needs phylogenetic context for new strains; cef-cached expansion needs clade-balanced cohort). See `plans/v0.1_Ingestion_Contract_Plan.md`.

**Status of cross-machine sync as of 2026-05-24:** Codex's pipeline.py changes + relocked spec + release packet + reports/ artifacts are on Precision 7780 only — NOT pushed to origin. Diagnostic tool at `scripts/cross_machine_sync_check.py` reports the drift. Next step: get Codex to push, OR transfer the divergent files manually.

## Coordination protocol with Codex (carries forward from coordination plan)

- Codex transfers `cipro_bounded_falsifier_results_<DATE>.{md,json}` + `cipro_leakage_check_dup_accession_<DATE>.{md,json}` into `Downloads/`.
- Claude executes Step 0 + branch dispatch on this laptop.
- PASS path Mash-cluster step requires Codex on Precision 7780 (Docker Desktop + Mash image). Claude drafts `scripts/mash_cluster_n147.py`; Codex runs it.
- Retrain (Step L) runs on Precision 7780 — Codex re-invokes `scripts/pipeline.py train`.
- All other steps (doc edits, code edits, test additions, ledger updates, git operations) run on this laptop.

**Escalation triggers** (per coordination plan Section 1):
- Schema drift detected at Step 0 → halt, request Codex re-emit results with explicit schema.
- Mash image pull fails → halt PASS path, request Codex docker-pull first.
- Test failure during Step P/F → halt commit, debug + report.
- Force-push needed → halt, ask user.
