# Independent TB — full-cohort lineage-collapsed number (2026-07-02)

The clonality-corrected headline for the **genuinely-independent** (out-of-CRyPTIC-build) TB validation.
Computed over the full assembly-available provenance-disjoint cohort; the raw number it supersedes was
clone-inflated. Artifact: `wiki/tb_independent_amr_portal_lineage_collapsed.json`.

## Cohort
- **N = 2,845** assembly-available disjoint isolates (of 26,941 not-leaked phenotyped; the 2,845 carry a
  downloadable GCA assembly = the fetchable/callable subset). `leaked=0` vs the CRyPTIC leakset; biosample
  cross-archive overlap 0/30 probed (`wiki/tb_independence_biosample_check.json`).
- Genotype: NCBI GCA assembly → minimap2 `asm5` + `paftools.js call` vs H37Rv → PASS-masked VCF (called
  2026-06-23, Docker; the 2,848 VCFs live on `D:/dna_decode_cache/tb_indep/vcf/`).
- Lineage: the pinned **Napier/tbdb barcode** read from each isolate's VCF (NO Mash, NO re-fetch).
- Scorer: the FROZEN `tb_amr.score_drug` + `clonality.cluster_weighted_confusion` (WHO catalogue UNCHANGED).

## Numbers (raw → clonality-corrected)

| drug | raw sens | raw spec | **lineage sens** | **lineage spec** | eff. lineages (R/S) | discordant |
|---|---|---|---|---|---|---|
| **RIF** | 0.920 | 0.955 | **0.444** `[0.246–0.663]` | **0.979** `[0.889–0.996]` | 20 / 47 | 43 |
| **INH** | 0.879 | 0.962 | **0.321** `[0.179–0.507]` | **0.972** `[0.858–0.995]` | 30 / 36 | 44 |

(95% Wilson CIs. `TB_SUBSET_PLUMBING` status = assembly-available subset, not a prevalence-preserving full
cohort — honest; the cohort IS independent, the earlier "in-distribution" label was a shared-fn copy artifact,
corrected in `tb_independent_lineage_collapse.py`.)

## Honest interpretation
1. **Raw sensitivity (0.92 / 0.88) is clonality inflation.** 2,845 isolates collapse to ~67 barcode-lineages;
   one over-sampled resistant clone carried the raw metric. The lineage-collapsed sens (0.44 RIF / 0.32 INH)
   is the honest per-lineage figure.
2. **Specificity is high and robust** (~0.97–0.98 lineage-level): determinant-absence → S generalizes well
   out-of-distribution.
3. **The independent lineage number ≈ the in-distribution CRyPTIC lineage number** (CLAUDE.md: CRyPTIC-parquet
   lineage-collapsed RIF 0.41 / INH 0.349). So at the lineage level the WHO-catalogue rule does **not degrade
   out-of-distribution** — it is modestly sensitive in BOTH regimes. That is the real positive: generalization
   holds; the modest sensitivity is a property of the catalogue+calling, not of independence.
4. **This is a conservative (lower-bound) sensitivity**, for two method reasons, both NAMED not hidden:
   (a) the assembly-`asm5` confident-difference VCF misses determinants a regeno/masked pipeline would call
   (callability unassessed → non-match counts as S); (b) the **Napier barcode is coarse** (~67 sub-lineages) —
   43–44 lineages are *discordant* (mixed R/S predictions within one barcode-lineage) and excluded, which
   depresses sens. A finer Mash-distance collapse would separate those and likely lift sens; barcode-collapse
   is the conservative view actually computed.
5. Wide CIs (RIF sens 0.25–0.66) reflect only ~20–30 effective R-lineages — read the interval, not the point.

## Status
Deliverable-b (independent TB validation) is **COMPLETE at both levels**: raw (`cefa97f`, verified 2026-07-02)
+ clonality-corrected (this run). Frozen AMR surface byte-unchanged (leak guard 9/9). Remaining refinement
(NOT blocking): a finer Mash-distance collapse + a regeno-VCF callability pass would tighten the sensitivity
bound — both deferred, both named.
