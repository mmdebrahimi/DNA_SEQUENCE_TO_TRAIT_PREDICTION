# TB RIF full-cohort v1b — first real lineage-collapsed number (2026-06-17)

**Status: `TB_SUBSET_PLUMBING`** — the FULL per-drug prevalence-preserving RIF cohort (8955 HIGH-quality
isolates, 3448R/5507S) but **callability pending the regeno fetch** (callability_assessed=False), so the
C3 gate honestly withholds the `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE` label. Artifact:
`wiki/tb_rif_cryptic_results_2026-06-17.json`.

## Result
| metric | value |
|---|---|
| raw sens / spec | 0.907 / 0.976 (TP 3129 / FP 133 / TN 5374 / FN 319) |
| **lineage-collapsed sens / spec** | **0.696 [Wilson 0.595–0.78] / 0.987 [0.954–0.996]** |
| effective lineages | **R 3448 → 92**, S 5507 → 154 |
| discordant / mixed-prediction clusters | 59 / 57 |
| unassigned lineage | 225 |

## Honest read
- **Well-powered, NOT demote-to-smoke-test.** 92 effective R-lineages (the 300-convenience PoC's 2 was a
  sampling artifact). Clonal inflation at full scale is moderate, not extreme (raw sens 0.907 vs lineage
  0.696 ≈ 21 pp; raw overstates because the determinant rule nails the big clones).
- **Spec is excellent** (0.987 lineage) — the rule rarely false-calls R.
- **Sens is limited by the determinant MATCHER, not biology.** Lineage sens 0.696 means ~30% of distinct
  R-lineages carry an rpoB determinant the EXACT `(pos,ref,alt)` match misses — the MNV / multi-encoding of
  the same AA change (S450L has several coord encodings) + out-of-RRDR determinants (rpoB I491F). The
  deferred Step-1 **codon-level / AA-equivalent normalization** is the clear next fix and would raise sens
  toward the rule's true performance.

## Caveats (do not over-read)
- PLUMBING, not BASELINE: callability not assessed (needs the ~1.6 TB regeno fetch -> D:, the multi-day job).
- The matcher recall gap means this 0.696 is a LOWER bound on the WHO-catalogue rule's lineage sens.
- Knowledge-baseline framing unchanged: the WHO catalogue was built partly from CRyPTIC -> in-distribution.

## Next (TB thread)
1. Codon-level determinant normalization (Step 1 enhancement) -> re-run -> expect higher lineage sens.
2. Regeno fetch (`scripts/populate_tb_regeno_detached.bat`) -> callability -> the true BASELINE label.
3. Hand-curated post-2023 independent gold set -> the only path to a SCIENTIFIC (non-baseline) number.
