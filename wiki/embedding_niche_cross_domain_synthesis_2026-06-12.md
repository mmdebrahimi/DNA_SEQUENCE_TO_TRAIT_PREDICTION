# Embedding-niche thesis — cross-domain synthesis (2026-06-12)

**Headline:** Genomic foundation-model embeddings capture **lineage / population structure, not
mechanism / causal signal** — consistently, across the bacterial→eukaryotic kingdom boundary, including in
the domain engineered to give embeddings their best shot. The embedding-niche hypothesis is **not
supported** after three independent de-confounded tests.

This is a **clean negative result, not an execution failure.** The project built the fairest possible test
for the embedding hypothesis and got a consistent answer. It sharpens strategy: stop hunting for the
embedding niche; invest in the deterministic decoder that demonstrably works.

## The three tests

| # | Domain | Model | Surface metric | De-confounded metric | Verdict |
|---|---|---|---|---|---|
| 1 | E. coli cipro AMR | Nucleotide Transformer | AUROC 0.914 (beat k-mer) | within-lineage ≈ chance; lost to QRDR knowledge baseline 0.943 | embedding learned **lineage** |
| 2 | AMR suite (10 cells) | deterministic rules + clonality audit | raw sens/spec inflated | cluster-weighted: cipro R-classes collapse (kleb 0.967→0.5, 2 lineages) | chromosomal-mutation resistance **tracks lineage** |
| 3 | Arabidopsis flowering time | PlantCaduceus_l32 | global r² −0.035 (beat structure −0.449) | within-group r² **−0.13** (3 seeds); spearman 0.22 < structure 0.48 | embedding learned **population structure** |

Each test isolates the same confound by a different de-confounder (held-out lineage / Mash-cluster weight /
kinship-PC group). Every time, the embedding's apparent win evaporates once population structure is removed.

## Why test #3 is decisive

Per the project's three-part embedding-niche criterion (`memory: feedback_embedding_niche_two_half_test`),
an embedding decoder needs: (a) a sampling-independent label, (b) no curated catalog to beat, (c)
organism-specific depth (≥100 same-organism strains). AMR fails (b) — AMRFinder is the catalog. Arabidopsis
flowering time satisfies **all three** (common-garden label; no catalog; 1003 accessions). It was the
designed best-case. The de-confounded result is negative. There is no remaining "but we didn't test it
fairly" escape.

## Strategic consequences

1. **The deterministic decoder suite is the product.** 10 lineage-disclosed SCORED AMR cells (E. coli +
   Klebsiella + Campylobacter × cipro/cef/gent/tet/mero), mostly lineage-robust on acquired-gene
   mechanisms. This is the "AI DNA decoder" north star — not the embedding bet.
2. **Acquired-gene resistance generalizes; chromosomal-mutation resistance is lineage-confounded.**
   β-lactamases, tet efflux, aminoglycoside-modifying genes survive lineage-correction; cipro QRDR point
   mutations do not (they are vertically inherited → over-sampled clones carry the raw metric).
3. **Do not spend money scaling embeddings.** A bigger GPU / larger window will not convert a negative
   de-confounded metric to positive. The question is answered at the structure-vs-signal level.
4. **The embedding research arm is closed**, honestly, with a recorded negative — not abandoned.

## Provenance

- Test 1: `memory: feedback_embedding_vs_knowledge_baseline_and_within_lineage`; cipro NT vs QRDR-POINT.
- Test 2: this session — `wiki/decoder_validation_report_card.md` + `wiki/provdisjoint_lineage_metrics.json`
  (clonality-disclosure layer, `dna_decode/eval/clonality.py`).
- Test 3: `wiki/phase2_arabidopsis_result_2026-06-12.md` ← `wiki/pathb_databricks_handoff_2026-06-12.md`
  (Databricks PlantCaduceus_l32, FT10 n=1003, 3 seeds).
