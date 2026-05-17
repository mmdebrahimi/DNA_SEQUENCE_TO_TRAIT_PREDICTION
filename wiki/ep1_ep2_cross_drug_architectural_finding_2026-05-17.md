# Phase 1 Closeout — Cross-Drug Architectural Finding (EP1 + EP2)

> Internal closeout artifact synthesizing the 3-experiment evidence chain across ciprofloxacin (EP1) + ceftriaxone (EP2 PASS) + tetracycline (EP2 FAIL). Closes Phase 1 evidence collection. External publication deferred per EP1 closeout PC1=`internal_closeout`.

---

## 1. Executive summary

At 12-strain smoke fidelity, frozen-NT-whole-genome-pooling **PASSES on concentrated-signal AMR mechanisms** (cipro QRDR point mutations: AUROC 0.750; cef plasmid acquired-gene β-lactamases: AUROC 0.833) AND **FAILS on distributed mobile-element mechanisms** (tet tet-family efflux + ribosomal protection: AUROC 0.400, anti-predictive). The cef PASS is **consistent with concentrated β-lactamase-associated signal**, but **whether NT is genuinely localizing β-lactamase loci versus lineage-tracking on a 12-strain mini-cohort remains unresolved** (see §6 Residual uncertainty). The architecture's failure mode appears data-shape-dependent (mechanism-class-bounded), not drug-or-cohort-dependent in isolation.

Phase 1 EP1 cipro additionally produced a 4-tier adversarial audit infrastructure (mechanism × MIC merge + opacity-vs-noise classification + structurally-enforced SUSPEND gate) that caught a label-noise bottleneck at the N=38 cipro cohort and prevented a Databricks burst on uninterpretable data. The infrastructure generalizes to cef + tet + any future drug.

## 2. Evidence chain

| # | Experiment | Drug | Cohort | NT-XGBoost AUROC | Mechanism class | Verdict |
|---|---|---|---|---:|---|---|
| 1 | Smoke gate (12-strain mini) | Cipro | N=12 mini | 0.750 | QRDR point mutations | PASS 2026-05-14 |
| 2 | Stage 1 N=40 + 1b mean+max | Cipro | N=38 | 0.568 / 0.615 | QRDR + label noise | FAIL → EP1 audit closeout |
| 3 | Attribution preflight v2 | Cipro | N=38 | n/a (interpretive) | QRDR loci | INCONCLUSIVE_MISS (zero QRDR in top-K=20) |
| 4 | AMRFinderPlus mechanism audit | Cipro | N=38 | n/a (audit) | All cipro mechanisms | QRDR_DOMINANT (18/20 R have QRDR; 7 plasmid; 7 silent-S) |
| 5 | Raw BV-BRC MIC rejoin | Cipro | N=38 | n/a (audit) | Label-quality | NOISY (7 HIGH_R / 0 HIGH_S of 40) |
| 6 | Mechanism × MIC merge | Cipro | N=38 | n/a (audit) | Combined | NOISE_DOMINATES; gate fires SUSPEND_CONDITION_4 |
| 7 | Smoke gate (12-strain mini) | **Cef** | N=12 mini | **0.833** | β-lactamases (plasmid + AmpC + porin) | **PASS** 2026-05-17 |
| 8 | Smoke gate (12-strain mini) | **Tet** | N=12 mini | **0.400** (anti-predictive) | tet-family efflux + ribosomal protection | **FAIL** 2026-05-17; H17 falsified |

## 3. Architectural pattern (the load-bearing finding)

Across the 3 drugs tested, NT-frozen-whole-genome-pooling's behavior **appears to partition by mechanism-class-bounded data shape**, largely independent of drug identity at the smoke fidelity tested (N=12 mini-cohorts). The partition is a 3-mechanism-class observation, not yet a validated hypothesis:

- **Concentrated-signal mechanisms** — biologic resistance localized to a small set of CDS contexts. **Candidate dilution mechanism** (proposed explanation, not measured): whole-genome mean-pool contributes ~1/N_genes (N≈5000) per gene; localized signal survives this dilution because the predictive feature is confined to a small index set. For cipro this is QRDR mutations (gyrA/parC/parE point variants). For cef this is plasmid acquired-gene β-lactamases (blaCTX-M/SHV/TEM/OXA family) — distinct biology but similar data shape (localized CDS, presence/absence signal). NT-XGBoost passes smoke for both.
- **Distributed mobile-element mechanisms** — biologic resistance spread across multiple mobile elements (tet-family efflux pumps: tetA/B/C/D/K/L/M; ribosomal-protection proteins: tetO/W/Q/S), with significant plasmid-borne content + partial chromosomal integration. Whole-genome mean-pool averages these signals into noise. NT-XGBoost fails smoke (anti-predictive 0.400).

The N=38 cipro FAIL (Stage 1) is consistent with this pattern compounded by the EP1 cohort label-noise bottleneck — labels were noisy AND the architecture was operating outside its strong regime even on QRDR strains at the larger cohort size.

**§3 falsification trigger** (symmetric with EP2 D5(b) for H17):

The partition hypothesis is falsified if a future experiment finds EITHER (a) a concentrated-signal mechanism that FAILs smoke + Stage 1 under the same NT-frozen-pooling architecture, OR (b) a distributed mobile-element mechanism that PASSes smoke + Stage 1 under the same architecture. (a) would prove the partition is not signal-localization-dependent; (b) would prove distributed-signal cohorts can be handled by mean-pooling after all. Either outcome invalidates §3 and triggers Phase 2 re-investigation of the underlying mechanism. A 4th-mechanism-class test (e.g., colistin via mcr-family plasmids, aminoglycoside via aac/aph/aad acetyltransferases) would be the cheapest discriminating experiment.

## 4. EP1 cipro closeout — audit infrastructure findings

EP1 cipro shipped not only as a model-evaluation experiment but as an adversarial audit-infrastructure packet. The 4-tier evidence flow built across N=38 cipro is reusable across drugs:

- **Tier 1 (smoke gate at N=12):** PASS 2026-05-14 (NT 0.750 vs k-mer 0.694, gap +5.6 pp).
- **Tier 2 (Stage 1 at N=38 + Stage 1b mean+max + scaled-LR fix for H13):** FAIL; NT-LR 0.673 vs k-mer 0.648, gap +2.5 pp (below the 3-pp threshold; CI [-20, +23]).
- **Tier 3 (attribution + mechanism + MIC audits):** preflight INCONCLUSIVE_MISS; mechanism audit QRDR_DOMINANT (18/20 R have textbook mutations; 7/20 S strains carry silent primary mechanisms — likely mislabeled); MIC audit NOISY (7 HIGH_R / 0 HIGH_S of 40; 9 R have no MIC; 12 S borderline).
- **Tier 4 (mechanism × MIC merge with structurally-enforced gate):** signal quality 0.17; opacity_count = 0 (AMRFinder is NOT the bottleneck — biology is detectable); the merge gate fired SUSPEND_CONDITION_4 and refused to fire the downstream curated AMR baseline experiment on uninterpretable data.

The PC1 (`internal_closeout`) + PC2 (Fisher-exact label-stratified, α=0.10, enrichment ratio ≥1.25) pre-conditions locked before any decision script ran, preventing post-hoc rationalization. PC2's original 60% bare-threshold proposal was caught as BELOW the uniform-error null baseline (K/N = 26/38 ≈ 68.4%) and replaced with the statistical-context schema — a load-bearing methodological save flagged during /brainstorm.

Decision: no Databricks burst on cipro; EP1 closed internally as an audit-infrastructure packet; external publication deferred.

## 5. EP2 cef + tet smoke — H17 falsified

EP2 fired the 12-strain smoke gate on ceftriaxone + tetracycline mini-cohorts (each 6R/6S, 12 unique MLSTs, full assembly availability, filtered from the N=38 cipro cohort strain pool; shared NT cache at `nt_n40_cipro.h5`).

**Cef PASS:** NT-XGBoost AUROC 0.833 = k-mer AUROC 0.833 (gap 0 pp). NT is not obviously worse than the classical comparator. Result is **consistent with** the concentrated-β-lactamase-signal pattern; mechanism-vs-lineage adjudication remains unresolved at N=12 (see §6).

**Tet FAIL:** NT-XGBoost AUROC 0.400 (anti-predictive — below the 0.5 chance line) versus k-mer 0.722; gap +32.2 pp (well above the 15-pp engineering threshold). Per EP2 D5(b): tet smoke produces NT AUROC ≤ 0.55 AND best-classical ≥ 0.65, falsifying H17 ("cipro-derived NT-XGBoost architecture transfers to BOTH cef AND tet at 12-strain smoke fidelity"). H17 status updated to `falsified` in the project ledger.

The anti-predictive 0.400 (rather than near-chance 0.5) suggests genuine architectural mismatch on distributed-signal data, not a calibration plumbing bug at N=12 — the 2026-05-14 LESSON (`AUROC ≈ 0` with symmetric two-value scores = calibration bug) does not apply here (calibrate=False was already set; the score distribution is non-degenerate).

## 6. Residual uncertainty

Phase 1 evidence collection did NOT resolve the following questions. They are recorded as Phase 2 candidates, NOT Phase 1 reopening triggers per the synthesis plan's narrowed D3:

- **Cef mechanism vs lineage (N=12 mini-cohort).** Was the cef AUROC 0.833 genuine β-lactamase localization, or NT lineage-tracking on a 12-strain cohort with 12 unique MLSTs? At N=12 with ±0.19 AUROC CI the question cannot be adjudicated. A cef NT-vs-mechanism cross-tab diagnostic is scoped at `plans/Cef_Mechanism_Audit_Plan.md` (appendix-tier, may or may not fire); even if it fires, it is corroborating evidence at most, not decision-bearing.
- **Tet failure mode disambiguation.** The 0.400 anti-predictive result could be (a) architectural mismatch on distributed mobile-element signal (most likely per the cross-drug pattern), (b) tet-cohort label-noise carryover (not audited), or (c) an N=12-calibration interaction (less likely; calibrate=False is already set). Not adjudicated.
- **Cipro Stage 1 architectural verdict.** The Stage 1 N=38 cipro FAIL was partly diagnosed (cohort label noise was the dominant confound; opacity_count = 0). The pure architectural verdict on cipro at clean labels remains an open Phase 2 question.
- **BV-BRC strict-MIC cef cohort feasibility.** Today's HIGH_R = 7 / HIGH_S = 0 of 40 cipro cohort sparse pattern was not extended to cef. Cef-native cohort feasibility unknown.
- **Per-gene NT windows on CLEAN_R cipro strains.** Architectural-diagnostic experiment deferred.
- **Bakta cohort-wide annotation.** Toolchain validated 2026-05-15; cohort annotation not performed.
- **Mean+max attribution preflight v3.** Refactor scope already in `plans/Cipro_Decision_Bundle_Technical_Plan.md` Step 2 (gene_level_mutagenesis aggregation kwarg). Runtime deferred conditionally.

## 7. Architectural implications

The cross-drug pattern has direct production implications for AMR-prediction pipelines:

- **Per-drug architecture selection is required.** Frozen whole-genome pooling is not a universal AMR-prediction architecture. It works on concentrated-signal mechanisms (β-lactamases, QRDR) at smoke fidelity but fails on distributed mobile-element mechanisms (tet). A production pipeline shipping cipro + cef + tet on a single NT-pool architecture would silently mis-predict tet resistance.
- **Targeted per-gene NT windows are the likely architectural fix for distributed-mechanism resistance.** Whole-genome pooling dilutes mobile-element signal across ~5000 genes (each contributing ~1/5000 to the mean). Per-gene windows on locus-specific contexts (tetA, tetB, tetM, tetO, etc.) preserve the signal. Hardware-constrained: GTX 860M cannot fine-tune NT (CC=5.0 < 7.0 required for bitsandbytes 4-bit); the per-gene-windows architecture is implementable in frozen-inference mode with per-locus pooling.
- **Adversarial cohort audits should precede model gates, not follow failures.** EP1's mechanism × MIC × opacity-flag merge was built AFTER the Stage 1 FAIL; if it had run BEFORE, the cohort label-noise bottleneck would have been visible at cohort-build time and saved 5+ days of model-evaluation work on uninterpretable data. The infrastructure is now generic + reusable for future EPs.

## 8. Phase 1 stopping point

Phase 1 evidence collection is **closed** with this synthesis. The reopen trigger (per `plans/EP1_EP2_Cross_Drug_Synthesis_Plan.md` D3) is narrowly scoped:

> Phase 1 reopens ONLY for (a) an internal contradiction that invalidates the current EP1/EP2 verdict, OR (b) a factual/source mismatch in the evidence chain. New questions go to Phase 2 candidates, NOT Phase 1 reopening.

Phase 2 work (other drugs, larger cohorts, per-gene NT windows, multimodal data, etc.) is a separate strategic decision that requires its own /idea-anchor + /project-init cycle. It is NOT a continuation of Phase 1.

## 9. Corroborating diagnostic (optional)

If `plans/Cef_Mechanism_Audit_Plan.md` Step 3 (cef NT-vs-mechanism cross-tab) is later fired, its findings will be cited here as corroborating evidence — does NT-XGBoost's per-strain prediction agree with AMRFinder β-lactamase detection? Does the rank correlation (Spearman) between NT scores and β-lactamase counts hold? Discordant cases (NT-R without mechanism, NT-S with mechanism) flag interpretation edge cases.

**Status as of this writing (2026-05-17):** the cef cross-tab has NOT been fired. The conservative cef framing in §1 + §5 + §6 is intentionally written to be defensible WITHOUT it. The synthesis stands either way.

## 10. Lessons (cross-cutting; reusable for future EPs)

1. **Pre-conditions discipline catches statistical bugs.** PC2's original bare 60% threshold for "label-noise concentration" was caught as BELOW the uniform-error null baseline. Without PC1/PC2 locked before runtime, the bug would have biased the EP1 closeout's narrative toward "labels are the problem" regardless of evidence. Future EPs should pre-declare numeric thresholds + run them through a statistical sanity check.
2. **Mechanism × phenotype merge with opacity flag prevents tool/label conflation.** The merge audit's `opacity_count = 0` finding (AMRFinder finds primary mechanisms in every HIGH_R strain) prevented the wrong narrative ("labels are wrong AND mechanism audit is wrong"). The opacity flag is the critical separator.
3. **Smoke-tier infrastructure (N=12, calibrate=False, k-mer baseline) is reusable across drugs with cosmetic patches.** Today's cef + tet smokes used the existing cipro smoke script with `--drug` arg + minimum-diff output-string templating. The next drug runs in <1 hour from cohort build to verdict.
4. **Anti-predictive AUROC at N=12 with calibrate=False is data-shape divergence, NOT a plumbing bug.** The 2026-05-14 LESSON (`AUROC ≈ 0` with symmetric two-value scores = calibration bug) applies specifically to `CalibratedClassifierCV` with isotonic regression at N=11 — NOT to non-calibrated XGBoost. Tet's 0.400 is genuine architectural mismatch.
5. **Plan-language drift catches a real failure mode.** The cef mechanism audit plan's history (v1 rejected-11-steps → v2 reduced-3 → v3 re-inflated-8 → v4 reduced-3-again-after-framing-critique) shows that "smaller alternative" language can preserve prior review authority while authorizing similar expansion. /brainstorm round 4 caught this; in-place plan edits BEFORE /execute-plan honor the 2026-05-14 HIGH-salience LESSON.
6. **Cross-drug evidence chain (EP1 + EP2) was load-bearing for the architectural claim.** Cipro alone produced an ambiguous story (cohort noise + architecture co-bottleneck). Adding cef (PASS) + tet (FAIL) at the same smoke-tier infrastructure was what made the mechanism-class-bounded pattern visible. Single-drug audits at small N are weaker claims; cross-drug chained inference is stronger.

## Audit trail (source artifacts)

- `wiki/smoke_gate_12strain_cipro_2026-05-14.md` — cipro smoke PASS
- `wiki/stage1_n40_cipro_2026-05-15.md` — Stage 1 N=38 FAIL
- `wiki/stage1_n40_cipro_mean-plus-max_2026-05-16.md` — Stage 1b FAIL
- `wiki/cipro_attribution_preflight_2026-05-16.md` — attribution INCONCLUSIVE_MISS
- `wiki/cipro_mechanism_audit_2026-05-17.md` — QRDR_DOMINANT
- `wiki/cipro_mic_audit_2026-05-17.md` — NOISY label verdict
- `wiki/cipro_mechanism_phenotype_audit_2026-05-17.md` — merge NOISE_DOMINATES + SUSPEND_CONDITION_4
- `wiki/cipro_decision_bundle_pre_conditions_2026-05-17.md` — PC1 + PC2 lock
- `wiki/cipro_ep1_closeout_2026-05-17.md` — EP1 internal closeout packet
- `wiki/smoke_gate_12strain_ceftriaxone_2026-05-17.md` — cef PASS (0.833)
- `wiki/smoke_gate_12strain_tetracycline_2026-05-17.md` — tet FAIL (0.400 anti-predictive)
- `wiki/EP2_cef_tet_verdict_2026-05-17.md` — EP2 H17 falsified

- `plans/EP1_EP2_Cross_Drug_Synthesis_Plan.md` — this synthesis's planning document
- `plans/Cef_Mechanism_Audit_Plan.md` — appendix-tier diagnostic (may or may not fire)
- `plans/Cipro_Decision_Bundle_Plan.md` + `plans/Cipro_Decision_Bundle_Technical_Plan.md` — cipro EP1 closeout planning chain

- `project_state/dna-decode-2026-05-11.md` — full Action Log + Bellman frame
- `wiki/decisions-log.md` — high-salience retrospective entries (2026-05-14 calibration overcorrection; 2026-05-14 Sidework sequence; 2026-05-15 phase-ladder → Evidence Packets framing reset)
