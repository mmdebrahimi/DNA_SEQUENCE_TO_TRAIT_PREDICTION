# M. tuberculosis RIF + INH baseline on CRyPTIC (in-distribution) — 2026-06-22

**The TB cell ran on real data for the first time.** The cell shipped 2026-06-17 but was data-blocked on
the ~1.6 TB regeno cohort. The CRyPTIC Zenodo dump (now on D:) ships the same determinant information as
`VARIANTS.parquet` (per-isolate genomic-nucleotide variants, e.g. `761155c>t`), which matches the WHO
catalogue's genomic `(pos,ref,alt)` determinants directly — so `scripts/score_tb_cryptic_parquet.py` is a
thin parquet→calls adapter over the FROZEN scorer (`tb_amr.score_drug` + `tb_lineage.lineage_clusters` +
`score_tb_cryptic.score_cohort`, all reused unchanged).

## Result

| Drug | n (measured HIGH) | RAW sens / spec | lineage-collapsed sens / spec | R-lineages / S-lineages |
|---|---|---|---|---|
| Rifampicin (rpoB) | 8,955 | 0.916 / 0.974 | 0.41 / 0.991 | 156 / 230 |
| Isoniazid (katG+inhA/...) | 9,518 | 0.889 / 0.989 | 0.349 / 0.995 | 212 / 203 |

Status: `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`. Artifacts:
`wiki/tb_{rif,inh}_cryptic_parquet_baseline_2026-06-22.json`.

## What this is — and is NOT (the honesty rails)

1. **In-distribution, NOT independent validation.** The WHO catalogue v2 was built partly FROM CRyPTIC, so
   scoring the catalogue rule on CRyPTIC is a knowledge baseline. A truly-independent TB number still needs
   a hand-curated post-2023 gold set (deliverable-b, `organism_rules/tb_goldset.py`,
   `INDEPENDENT_VALIDATION_BLOCKED_NO_GOLDSET`).
2. **Raw isolate-level is CLONALITY-INFLATED.** Raw sens ~0.9 collapses to ~0.35–0.41 at the sublineage
   level (Napier barcode): the catalogue catches the dominant R clones but a minority of distinct R
   sublineages. This is the project's lineage-disclosure discipline applied to TB; the raw number is NOT
   the honest headline.
3. **Lineage-collapsed sens is a LOWER BOUND (parser scope).** The adapter parses SNVs + (via per-base
   `snv_components`) codon-level MNVs, but NOT true indels / delins / `katG_LoF` frameshifts. Isolates whose
   only resistance determinant is an indel are called S (FN) — RIF is barely affected (RRDR is SNV/MNV);
   INH loses the katG-LoF minority (raw 0.889 is a few pp below the catalogue's true INH sens). Closing this
   needs indel normalization across the CRyPTIC↔WHO-coords representations (code-closable; deferred).
4. **Callability unassessed** (no regeno): a non-match is S, never ABSTAIN.

## Forward
- Independent TB number: hand-curate a post-2023 gold set (the gated deliverable-b).
- Tighten this baseline: add indel/delins/LoF normalization to the parquet adapter (code-closable).
