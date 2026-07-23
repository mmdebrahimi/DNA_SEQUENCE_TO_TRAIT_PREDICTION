# Clinical-gene landscape census — the data hunt that doubled the R2 decoder substrate (2026-07-22)

**Status:** ✅ a data-hunt deliverable. User directive "hunt for more data." The R2 molecular regime (DMS +
ClinVar) clears all 8 rejection gates in `wiki/negative_results_map_2026-06-13.md` by construction
(wet-lab/clinical labels, sampling-independent, non-circular, no clonality axis) — so the richest untapped
free data is **more MaveDB clinical genes**. This census maps that landscape and DOUBLES the AUROC-viable
clinical-decoder substrate. Frozen AMR surface byte-unchanged (READ-only).

## The hunt

MaveDB (CC0) holds **2,798 score sets → 560 distinct human protein-coding genes → 490 with a UniProt id**.
The join to ClinVar/AlphaMissense was previously thought blocked by a numbering OFFSET (MaveDB assays are
often domain/construct-numbered). **The offset was never a wall — MaveDB ships it in the metadata**
(`targetGenes[].externalIdentifiers[] {dbName: UniProt, offset: N}`), so `UniProt_pos = mavedb_pos + offset`
is AUTO-derivable (MLH1 offset 486, PSD95 309, p53 101, …). Applying it took MLH1's ClinVar join from **0 → 29**.

## Result (`scripts/clinical_gene_landscape_census.py`; real MaveDB + ClinVar + AlphaMissense on D:)

Censused 28 clinically-actionable genes → **4 are AUROC-viable (up from 2)**:

| Gene | AlphaMissense AUROC | BLOSUM floor | n (path/benign) | clinical relevance |
|---|---|---|---|---|
| **TP53** | 0.986 | 0.714 | 262 (173/89) | Li-Fraumeni / pan-cancer |
| **LDLR** | 0.975 | 0.748 | 145 (115/30) | **familial hypercholesterolemia (NEW)** |
| **F9** | 0.949 | 0.717 | 219 (197/22) | **hemophilia B (NEW)** |
| **MSH2** | 0.936 | 0.832 | 297 (64/233) | Lynch syndrome |

**LDLR and F9 are newly decoded** — unlocked by the offset crack + a curated gene→UniProt fallback (MaveDB's
per-assay UniProt cross-ref is inconsistent; TP53's highest-variant assay omits it, which had mis-marked TP53
`NO_UNIPROT_ID` in the first pass). Every viable gene shows the R2 signature: AlphaMissense near the ceiling
(0.936–0.986), a large gap over the deterministic BLOSUM floor (0.71–0.83) = the learned-decoder headroom.

## Full census state distribution (28 genes)

| State | Count | Meaning |
|---|---|---|
| AUROC_VIABLE | 4 | TP53, MSH2, LDLR, F9 — both classes ≥15 in the DMS-covered region |
| SINGLE_CLASS | 18 | joins fine (offset applied) but ClinVar-missense is one-sided in the covered region |
| NO_MAVEDB_HUMAN_DMS | 5 | MSH6, PMS2, PALB2, VHL, RAD51C — no human DMS in MaveDB |
| NO_JOIN | 1 | SCN5A (the 252-variant assay carries offset 1620 that doesn't map — honest edge case) |

**The dominant non-viable state is SINGLE_CLASS, and it is a real finding about the DATA, not a decoder
failure:** clinically-curated ClinVar missense for most tumor-suppressor / disease genes is
**pathogenic-dominated** (KRAS 40/0, GCK 294/5, CBS 107/5, PTEN 199/1) — benign missense are rarely submitted.
AUROC needs both classes, so a gene is viable only where the DMS-covered region also carries ≥15 curated
benign missense (TP53/MSH2/LDLR/F9 do; most don't). This is the same gene-dependent balance seen earlier for
BRCA1/PTEN, now mapped across 28 genes.

## What this unlocks + honest scope

- **The offset is now auto-applied**, so any of the 490 UniProt-bearing human MaveDB genes can be censused by
  adding it to the gene list — the census is the reusable substrate map, not a one-off.
- **4 AUROC-viable genes** is the honest expandable clinical-validation set today (2× the prior TP53/MSH2).
  Widening it further is DATA-gated (needs genes with both a DMS AND curated-benign ClinVar in the covered
  region), NOT code-gated — a genuine label-availability wall, consistent with the project's north-star
  finding that labels, not models, are the binding constraint.
- Tier = in-distribution clinical (AM saw these proteins; ClinVar labels independent of the DMS-fitness tuning).
- SCN5A's off=1620 NO_JOIN shows the auto-offset is not infallible for tiny/odd assays — reported, not hidden.

## Reproduce

```
uv run python scripts/clinical_gene_landscape_census.py            # 28-gene clinical census (real network)
uv run python scripts/clinical_gene_landscape_census.py --genes LDLR,F9,TP53
```

MaveDB landscape cached at `D:/dna_decode_cache/mavedb/human_landscape.json`; ClinVar E-utilities per gene
cached at `D:/dna_decode_cache/clinvar/eutils/`; AM filter at `D:/dna_decode_cache/alphamissense/`. The
ClinVar fetch now has 429-backoff (the rate-limit that killed the first run) + honors `NCBI_API_KEY`.
Artifact JSON: `wiki/clinical_gene_landscape_census_2026-07-22.json`. Builds on
`scripts/clinical_am_hybrid_auroc.py` + `scripts/clinical_variant_effect_validate.py`. 4 offline tests
`tests/test_clinical_gene_landscape_census.py`.
