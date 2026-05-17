# Cipro Decision Bundle — Pre-Conditions (2026-05-17)

> Structurally-enforced inputs for `plans/Cipro_Decision_Bundle_Technical_Plan.md` Wave 7+. Pre-declared per the technical plan's PC1/PC2 requirements. `scripts/cipro_decision_cell.py` (Step 10.7) refuses to write the runtime artifact unless this file exists.

---

## PC1: Phase 1 EP1 deliverable framing

**Value:** `internal_closeout`

**Definition:** EP1 cipro is closed as an adversarial audit/infrastructure packet showing the N=38 BV-BRC cipro cohort is not a defensible model gate, and that frozen whole-genome NT pooling did not rescue it. The deliverable is a wiki packet + project_state ledger entry — NOT a Databricks burst, NOT an external publication (arXiv / blog) at this time.

**Decision rationale (from `/brainstorm` 2026-05-17 Round 2):**
- Solo personal project; no external deadline; no team or stakeholder demanding a classifier.
- Databricks burst on N=150 strict-MIC cohort doesn't fix the architectural bottleneck (GTX 860M can't fine-tune; whole-genome pooling dilutes QRDR signal across ~5000 genes).
- AMRFinderPlus correctly identifies cipro resistance in the SAME genomes NT misses — strong evidence the biology is detectable but NT-frozen-pooling is the wrong tool. More data doesn't fix tool mismatch.
- The publishable object is the integrated audit stack (MIC rejoin + AMRFinder mechanism audit + mechanism × phenotype merge + structurally-enforced SUSPEND gate), not the failed model. External publication of this audit stack is **deferred** until EP2 (cef + tet smoke) tells us whether the architecture transfers — that result determines whether the cipro failure is drug-specific or architecture-wide.

**Allowed downstream actions under `internal_closeout`:**
- Write `wiki/cipro_ep1_closeout_<date>.md` summarizing the evidence chain.
- Update `project_state/dna-decode-2026-05-11.md` Bellman frame: EP1 closed; EP2 elevated to primary EP.
- Fire EP2 cef + tet smoke (NT-XGB minimum, mechanism-aware baseline expansion if interesting).
- Revisit publication framing post-EP2.

**Disallowed under `internal_closeout`:**
- Databricks burst for N=150 cipro cohort expansion.
- Cipro Stage 2 retraining.
- arXiv preprint or blog post (deferred, not forbidden).

## PC2: Numeric threshold for D4 estimand (error concentration on noisy strains)

**Value:** **Not a bare fraction.** Use stored statistical context per the schema below.

**Schema:**
```json
{
  "pc2_test": "fisher_exact_label_stratified",
  "pc2_alpha": 0.10,
  "pc2_enrichment_ratio_min": 1.25,
  "pc2_noise_class_definitions": [
    "strict_NOISY_only",
    "include_SUSPECT_S"
  ]
}
```

The error-audit script (`scripts/cipro_error_audit.py`, Step 10.5) must:
1. Compute observed error fraction under BOTH `strict_NOISY_only` and `include_SUSPECT_S` definitions; report each separately.
2. Compute the null fraction (= K/N where K = |noise_class subset| in the effective Stage 1b cohort, N = effective Stage 1b cohort size) for each definition.
3. Run a label-stratified Fisher exact (or Freeman-Halton) test conditional on the observed true-label distribution.
4. Emit observed_fraction, null_fraction, p_value, and enrichment_ratio (= observed_fraction / null_fraction).
5. Decision rule: **label-noise concentration is supported** iff (p_value < pc2_alpha) AND (enrichment_ratio >= pc2_enrichment_ratio_min) for AT LEAST ONE definition.

**Why not a bare fraction (from `/brainstorm` 2026-05-17 Round 1 statistical analysis):**
- The originally-proposed 60% threshold is BELOW the uniform-error null baseline (K/N = 26/38 = 68.4%). It would fire on uniformly-distributed errors — biased decision toward "label noise."
- A pure α=0.05 binomial threshold (k=14/15 = 93.3%) is severely underpowered at n=15 errors total. Even strong signal (e.g., 13/15) wouldn't reject.
- The label-stratified Fisher exact correctly handles small samples + sparse classes + the structural fact that all clean strains in the cohort are R (no clean S to balance).
- Reporting BOTH `strict_NOISY_only` and `include_SUSPECT_S` definitions handles estimand instability — the "what counts as noisy" choice changes the answer; reporting both makes the decision auditable.

**Pre-declared expected outcome based on existing Stage 1b per-strain table (`wiki/stage1_n40_cipro_mean-plus-max_2026-05-16.md`):**
- ~15 NT-LR errors observed.
- Strict NOISY_* errors ≈ 9/15 (60%) vs null ≈ 66% → expected to FAIL the test (p > 0.10; enrichment_ratio ≈ 0.91).
- +SUSPECT_S errors ≈ 13/15 (87%) vs null ≈ 82% → expected to FAIL the test (p > 0.10; enrichment_ratio ≈ 1.06).

**Pre-declared interpretation:** if the error audit FAILS to reject under both definitions (as expected), the EP1 closeout narrative is:
> The cohort is structurally noisy (per merge audit), AND the NT model errors do not specifically concentrate on the noisy strains — suggesting that frozen-NT-pool is failing on both clean and noisy strains. This is consistent with an ARCHITECTURE bottleneck on top of (not just instead of) a label bottleneck.

This is a more honest narrative than "labels were the bottleneck."

## Decision matrix (Step 11 reference)

| `decision_cell` (from census + error audit) | PC1 = `internal_closeout` | recommended_next_step |
|---|---|---|
| True_high_threshold | (unreachable under PC1) | n/a |
| True_low_threshold | (unreachable under PC1) | n/a |
| False_low_threshold | EP1 closed; cef + tet smoke next; defer publication | Close EP1; fire EP2 smoke. |
| False_high_threshold | EP1 closed; cef + tet smoke next; defer publication | Close EP1; fire EP2 smoke. |
| AMBIGUOUS_* | Same as False_*_threshold | Close EP1; fire EP2 smoke. |

Under PC1 = `internal_closeout`, no decision_cell value authorizes Databricks burst or Stage 2 expansion. The cell value influences the closeout narrative's tone (label-noise-dominant vs architecture-dominant) but does not unlock burst spend.

## Ordering (locked)

1. **Cipro micro-close (now):** write `wiki/cipro_ep1_closeout_2026-05-17.md`; update project ledger; commit.
2. **EP2 cef + tet smoke (next):** ~1-2 days; minimum viable NT-XGB-only smoke first; expand to mechanism-aware baselines only if NT shows signal on either drug.
3. **Revisit publication framing (post-EP2):** if EP2 confirms architecture-wide failure, the publishable angle is "small AMR cohorts need adversarial cohort audits"; if EP2 shows NT works on cef/tet, the angle becomes "cipro QRDR mechanism resists whole-genome pooling specifically."

## Audit trail

- /brainstorm 2026-05-17 Rounds 1 + 2 produced the Q1 / Q2 analysis (transcripts in conversation).
- /review 2026-05-17 (technical plan) produced the structurally-enforced pre-conditions design.
- Source data: today's audit JSONs (`wiki/cipro_*_audit_2026-05-17.{md,json}`) + Stage 1b verdict (`wiki/stage1_n40_cipro_mean-plus-max_2026-05-16.md`).
