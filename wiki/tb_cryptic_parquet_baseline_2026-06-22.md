# M. tuberculosis RIF + INH baseline on CRyPTIC (in-distribution) — 2026-06-22

**The TB cell ran on real data for the first time.** The cell shipped 2026-06-17 but was data-blocked on
the ~1.6 TB regeno cohort. The CRyPTIC Zenodo dump (now on D:) ships the same determinant information as
`VARIANTS.parquet` (per-isolate genomic-nucleotide variants, e.g. `761155c>t`), which matches the WHO
catalogue's genomic `(pos,ref,alt)` determinants directly — so `scripts/score_tb_cryptic_parquet.py` is a
thin parquet→calls adapter over the FROZEN scorer (`tb_amr.score_drug` + `tb_lineage.lineage_clusters` +
`score_tb_cryptic.score_cohort`, all reused unchanged).

## Result (SNV + exact indel-determinant match — updated 2026-06-22)

| Drug | n (measured HIGH) | RAW sens / spec | lineage-collapsed sens / spec | R-lineages / S-lineages |
|---|---|---|---|---|
| Rifampicin (rpoB) | 8,955 | **0.918** / 0.974 | 0.41 / 0.991 | 156 / 230 |
| Isoniazid (katG+inhA/...) | 9,518 | 0.889 / 0.989 | 0.349 / 0.995 | 212 / 203 |

Status: `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`. Artifacts:
`wiki/tb_{rif,inh}_cryptic_parquet_baseline_2026-06-22.json` (carry `indel_matching` + `raw_snv_only` deltas).

**Indel normalization is now implemented (was the deferred lower-bound caveat).** The adapter exact-matches
WHO indel/`katG_LoF` determinants against CRyPTIC `pos_del_x` / `pos_ins_x` rows (conversion VERIFIED against
the data, `scripts/_tb_indel_probe.py` retired into tests). Effect (before → after):
- **RIF 0.916 → 0.918** raw sens: 8 isolates flipped S→R, **all 8 true positives, 0 false positives** (spec
  unchanged) — recovered real rpoB Phe433/Thr444-dup + RRDR-deletion FNs.
- **INH 0.889 → 0.889** (no change): the 924 mapped katG-LoF determinant strings are correct (7/9 spot-checked
  present in the data), but their carrier isolates fall outside the 9,518 HIGH-quality-INH-labelled cohort (or
  in the small repeat-region left-align residual below) — the measured **labels** aren't there to score, even
  though the genotype match works. Bounded at ≤12/9,518 ≤0.13 pp.

The match is EXACT (no positional tolerance) → **zero false-positive risk** (RIF FP=0 confirms it). The SNV-only
lower bound was already nearly tight; this removes the caveat and makes the number complete, not a floor.

## What this is — and is NOT (the honesty rails)

1. **In-distribution, NOT independent validation.** The WHO catalogue v2 was built partly FROM CRyPTIC, so
   scoring the catalogue rule on CRyPTIC is a knowledge baseline. A truly-independent TB number still needs
   a hand-curated post-2023 gold set (deliverable-b, `organism_rules/tb_goldset.py`,
   `INDEPENDENT_VALIDATION_BLOCKED_NO_GOLDSET`).
2. **Raw isolate-level is CLONALITY-INFLATED.** Raw sens ~0.9 collapses to ~0.35–0.41 at the sublineage
   level (Napier barcode): the catalogue catches the dominant R clones but a minority of distinct R
   sublineages. This is the project's lineage-disclosure discipline applied to TB; the raw number is NOT
   the honest headline.
3. **Indel/delins/`katG_LoF` now matched (was the lower-bound caveat); a tiny residual remains, NAMED.**
   The adapter exact-matches indel determinants (SNV + per-base-MNV + indel). Residual, bounded + documented:
   (a) **6 complex delins** per drug where alt/ref share no clean prefix (need reference left-alignment) are
   skipped; (b) a few **repeat-region insertions** left-align ±1 between the CRyPTIC and WHO-coords
   conventions (exact match misses them; ±1 fuzz was deliberately NOT added — it would trade the zero-FP
   guarantee for ~zero gain). Net effect on the headline is ≤0.2 pp (RIF) / ≤0.13 pp (INH).
4. **Callability unassessed** (no regeno): a non-match is S, never ABSTAIN.

## Forward
- Independent TB number: hand-curate a post-2023 gold set (the gated deliverable-b) — BLOCKED:external on a
  measured per-isolate label (see `wiki/tb_goldset_public_source_exhaustion_2026-06-22.md`; author-request
  emails drafted at `wiki/tb_goldset_author_emails_2026-06-22.md`).
- Indel normalization: DONE (this update). Remaining residual (complex delins + repeat left-align) is
  reference-based normalization for ≤0.2 pp — not worth the complexity unless an independent gold set with
  high per-isolate indel burden lands.
