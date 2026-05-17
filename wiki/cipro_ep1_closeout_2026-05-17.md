# EP1 Cipro — Internal Closeout (2026-05-17)

> EP1 (cipro Phase 1 deliverable) is closed as an adversarial audit/infrastructure packet. Verdict: the N=38 BV-BRC cipro cohort is not a defensible model gate, and frozen whole-genome NT pooling did not rescue it. No Databricks burst. External publication deferred. Next EP: EP2 cef + tet smoke gate.

---

## Summary

After 5+ experiments and a 4-tier audit pipeline, EP1 cipro is closed internally. Evidence is consistent with TWO bottlenecks operating simultaneously:

1. **Cohort bottleneck (confirmed):** the N=38 BV-BRC cipro cohort has structurally noisy labels. Only 7 of 20 R-strains have decisive MICs; 0 of 20 S-strains have decisive MICs; 6 S-strains carry primary cipro mechanisms (likely mislabeled).
2. **Architecture bottleneck (suggested):** NT-frozen-whole-genome pooling fails to localize the QRDR mutations that AMRFinderPlus correctly identifies in 18/20 R strains. NT-LR errors do NOT concentrate on the noisy strains specifically — they're roughly uniform across clean and noisy strains, suggesting the pooling architecture is failing on BOTH.

Neither bottleneck is fixable at this cohort size with this hardware. The infrastructure built today (verify_complete cache integrity gate, mechanism-phenotype merge with opacity flag, structurally-enforced SUSPEND_CONDITION_4 gate) is preserved for EP2 + future EPs.

## Evidence chain

| # | Experiment | Verdict | Key finding |
|---|---|---|---|
| 1 | Stage 1 (mean-pool NT-XGB) | FAIL | NT-XGB 0.568 vs k-mer 0.648; gap -8.0 pp |
| 2 | Stage 1b (mean+max + scaled NT-LR) | FAIL | NT-LR 0.673 vs k-mer 0.648; gap +2.5 pp (below 3 pp threshold) |
| 3 | Attribution preflight v2 | INCONCLUSIVE_MISS | Zero QRDR + expanded cipro loci in top-K=20 across 19 R strains |
| 4 | AMRFinderPlus mechanism audit | QRDR_DOMINANT | 18/20 R have textbook gyrA/parC/parE; 7/20 R have plasmid quinolone protection; 7/20 S carry primary mechanism (silent rate 35%) |
| 5 | Raw BV-BRC AST/MIC rejoin | NOISY | 7 HIGH_R + 0 HIGH_S of 40; 9 R have no MIC; 12 S borderline |
| 6 | Mechanism × MIC merge | NOISE_DOMINATES | clean_count=7 (7R/0S), opacity_count=0, signal quality 0.17, gate fires SUSPEND_CONDITION_4 |

## Pre-declared error-audit interpretation (per `wiki/cipro_decision_bundle_pre_conditions_2026-05-17.md`)

Based on existing Stage 1b NT-LR per-strain predictions:
- ~15 NT-LR errors observed total.
- Strict NOISY_* errors ≈ 9/15 (60%) vs null baseline ≈ 66% → enrichment ratio ≈ 0.91 (below 1.25 threshold).
- +SUSPECT_S errors ≈ 13/15 (87%) vs null baseline ≈ 82% → enrichment ratio ≈ 1.06 (below 1.25 threshold).
- **Expected verdict:** Fisher exact label-stratified test fails to reject uniform-error null under both noise-class definitions. **Label-noise concentration NOT established.**

## What this means scientifically

The clean version of the EP1 narrative is:

> The N=38 cipro cohort is structurally noisy AND the NT model errors do not concentrate specifically on the noisy strains. This is inconsistent with a single-bottleneck story (labels-only or architecture-only) and consistent with both bottlenecks operating simultaneously. The cohort cannot be used as a model gate at this size, AND the NT-frozen-pooling architecture appears unable to localize the QRDR mutations that the AMRFinderPlus curated catalog correctly identifies in the same genomes.

The strongest current claim is therefore neither "NT doesn't work for cipro" nor "BV-BRC labels are unusable" but:

> **Frozen-NT-whole-genome-pooling on N=38 BV-BRC cipro is not a research question with a clean answer; the right next experiment is EP2 (cef + tet smoke) to determine whether the architecture transfers to drugs with different resistance biology.**

## Decisions locked

1. **No Databricks burst on cipro.** N=150 strict-MIC expansion is feasibility-uncertain (7 HIGH_R + 0 HIGH_S of 38 suggests the BV-BRC universe is sparse at strict thresholds) AND doesn't fix the architecture bottleneck.
2. **No cipro Stage 2 retraining on this cohort.** SUSPEND_CONDITION_4 stands.
3. **External publication (arXiv / blog) deferred to post-EP2.** The publishable angle depends on whether EP2 shows architecture transfers.
4. **Infrastructure preserved + reused.** Today's audit pipeline (mechanism × MIC merge + opacity flag + structurally-enforced gate) is generalizable to cef + tet for EP2 and to other drugs in Phase 2.

## What's NOT closed

The following remain open + are taken up by EP2 + later EPs:

- **H17 (cipro-NT architecture transfers to cef + tet at smoke fidelity).** Falsified or confirmed by EP2.
- **Per-gene NT windows on the 7 CLEAN_R cipro strains.** Architectural-diagnostic experiment; not a Phase 1 ship gate; queue for Phase 2 if EP2 shows architecture-wide failure.
- **Bakta cohort-wide annotation.** Tier 2 in `plans/Cipro_Decision_Bundle_Plan.md`; lifted out. Queue for Phase 2 if Path A becomes alive.
- **Mean+max attribution preflight v3.** Refactor scope already in `plans/Cipro_Decision_Bundle_Technical_Plan.md` Step 2; runtime use deferred per the technical plan's decision matrix.

## What ships as the EP1 deliverable

- This packet (`wiki/cipro_ep1_closeout_2026-05-17.md`).
- The pre-conditions artifact (`wiki/cipro_decision_bundle_pre_conditions_2026-05-17.md`).
- The 6 prior audit packets (Stage 1, Stage 1b, preflight v2, mechanism audit, MIC audit, merge).
- The infrastructure code (mechanism × phenotype audit scripts, AMRFinder mechanism audit script, MIC rejoin script, decision-cell wrapper design, runtime artifact schema). Production code waves (technical plan Steps 1-10.8) can ship as the EP1-supporting infrastructure deliverable when convenient.
- Project ledger entry (`project_state/dna-decode-2026-05-11.md` Action Log + Bellman frame refresh).

## Next step

**EP2 cef + tet smoke gate.** Plan exists at `plans/EP2_Cef_Tet_Smoke_Design_Plan.md`. Minimum-viable scope: NT-XGB-only smoke (no mechanism-aware baseline) for cef + tet at the existing N=12 mini-cohort. Wall-time: ~hours. Expand to mechanism-aware baseline only if NT shows signal on either drug. EP2 result determines whether cipro's failure is drug-specific (QRDR-mechanism-specific) or architecture-wide.

## Lessons (for future EPs)

1. **Adversarial cohort audits before model gates.** The mechanism × phenotype × MIC merge caught what neither the model verdict nor the attribution preflight alone would have caught: the cohort itself isn't decidable at this size. Future EPs should fire this audit BEFORE the model gates, not after the model has failed.
2. **Two-bottleneck reasoning.** Don't assume failure modes are single-cause. Today's evidence is consistent with cohort noise AND architecture issues simultaneously; insisting on one or the other would have been wrong.
3. **PC1/PC2 pre-condition discipline.** Locking the framing + numeric threshold BEFORE the runtime decision fires prevents post-hoc rationalization. Codex's Q2 statistical correction (60% threshold was below the uniform-error null) would have biased the EP1 outcome toward "labels were the bottleneck" without the pre-condition discipline.
4. **K/N noise prevalence ≠ error concentration.** A cohort being 68% noisy is not the same as a model's errors being 68% concentrated on the noisy strains. The two statistics measure different things; conflating them is a common analytical mistake.
