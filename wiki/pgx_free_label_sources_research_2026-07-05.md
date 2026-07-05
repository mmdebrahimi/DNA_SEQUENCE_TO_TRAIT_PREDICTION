# PGx free-label sources beyond the CYP2C cluster — research memo (2026-07-05)

**Question (from the prior session's open frontier):** which *free, independent* label source extends the
validatable PGx decoder surface beyond the CYP2C cluster (CYP2C8/9/19, all GeT-RM-validated) + VKORC1?

## Answer — the full GeT-RM CDC reference-material consensus

The ursaPGx benchmark truth I used for CYP2C8/9/19 (`star-allele-comparison_common.tsv`) is a **narrow
4-gene** set (CYP2C8/9/19 + CYP2D6). The `_merged.tsv` has more *samples* (3261 rows) but the SAME 4 genes.
The broader free source is the **CDC-based Genetic Testing Reference Materials (GeT-RM) program** consensus:

- GeT-RM characterized **137 Coriell samples across 28 genes** (consolidated tool: **363 samples, 34
  loci**), each by **≥2 labs / multiple methods** → a per-sample **consensus diplotype**. Per-gene panels
  published in *J Mol Diagn*: TPMT/NUDT15 (2022), CYP2C8/9/19 (2022), **CYP3A4/CYP3A5 (2023)**, DPYD, etc.
- **59–88 of these samples overlap the public 1000 Genomes panel** (the Coriell NA/HG IDs are 1000G IDs),
  so the consensus diplotype is directly joinable to a 1000G-region VCF — the SAME validation shape as the
  CYP2C8 cell (independent caller vs the accepted consensus truth set).
- **Per-sample consensus tables are downloadable** from the CDC GeT-RM site as Excel/Word
  (`cdc.gov/lab-quality/php/get-rm/reference-materials.html`) — e.g. `CYP3A4-CYP3A5-website-table.xlsx`
  (fetched this session). Also the *J Mol Diagn* paper supplements + PharmCAT's Sangkuhl 2020 Table S1.

## The load-bearing HONESTY nuance — three tiers, do not conflate

Not every gene extension is an *independent* validation. Screen each candidate:

1. **Real GeT-RM consensus (wet-lab/multi-method) → NEAR_INDEPENDENT** (the CYP2C tier). Genes with a
   genuine multi-lab consensus on ≥1000G-overlap samples: CYP2C8/9/19, CYP2D6, CYP3A4, **CYP3A5**, CYP4F2,
   TPMT, NUDT15, DPYD, CYP2B6, UGT1A1. PharmCAT-vs-GeT-RM concordance (Sangkuhl 2020): CYP3A5 **59/59**,
   TPMT 58/59 — real, high.
2. **1KGP-PharmCAT-derived truth → FAITHFUL_TO_TOOL, NOT independent.** For SLCO1B1, G6PD, and (post
   haplotype-def updates) DPYD/NUDT15, the field generates the "truth" by running PharmCAT/pypgx on the
   1000G VCFs — because GeT-RM has too few samples or its calls are outdated. Validating our caller against
   *that* proves we match the tool, NOT an independent label. Report it as faithful-to-tool (like the
   typing/finder cells), never as an independent number.
3. **Single-SNP star systems → calling is near-tautological.** For CYP3A5*3 (rs776746), SLCO1B1*5
   (rs4149056), CYP4F2*3 (rs2108622), NUDT15*3 (rs116855232), the star call is a deterministic function of
   ONE genotype — so "GeT-RM concordance" mostly re-checks the same SNP the truth used. The real
   independent checks for these are **trio-Mendelian consistency** (we have `scripts/pgx_trio_concordance.py`)
   + **allele-frequency sanity** vs gnomAD. The independent-*calling* value is highest where an allele is
   MULTI-variant or aliases another (CYP2D6 structural; TPMT *3A = *3B+*3C; CYP2B6 *6 = two SNPs; the
   CYP3A5 *6/*7 no-function alleles that a *3-only caller mis-reads as *1/expressor).

## Prioritized build order (judgment)

| Gene | Drug | Defining | Truth tier | Verdict |
|---|---|---|---|---|
| **CYP3A5** | tacrolimus (transplant) | *3 rs776746 / *6 rs10264272 / *7 rs41303343(ins) | **real GeT-RM consensus** | **BUILD FIRST** — high value; *6/*7 give real multi-allele content (expressor vs non-expressor); CPIC phenotype exists |
| SLCO1B1 | statins (myopathy) | *5 rs4149056 (single SNP) | 1KGP-PharmCAT / tautological | build, tier honestly (faithful-to-tool / genotype+trio), NOT independent |
| CYP2B6 | efavirenz | *6 = rs3745274+rs2279343 | real GeT-RM consensus | strong *calling* validation (2-SNP) — good later target |
| TPMT | thiopurines | *3A = *3B(rs1800460)+*3C(rs1142345) | real GeT-RM (58/59) | strong *calling* validation (compound) — good later target |
| CYP2D6 | many | structural/hybrid/CNV | real GeT-RM | HARD (SV/CNV) — needs a Cyrius-class caller, not the SNP framework |

**This run:** build **CYP3A5** as the flagship extension (real GeT-RM consensus, obtained this session:
`CYP3A4-CYP3A5-website-table.xlsx` → 9 CYP3A5 samples with *1/*3/*6/*7 consensus, 1000G-overlapping →
UNDERPOWERED but genuine). SLCO1B1 as a secondary, honestly tiered. CYP2B6/TPMT/CYP2D6 = named future
targets.

**Sources:** GeT-RM CDC (`cdc.gov/lab-quality/php/get-rm`) · Sangkuhl 2020 PharmCAT (*Clin Pharmacol Ther*)
· J Mol Diagn GeT-RM characterization papers (TPMT/NUDT15 2022; CYP2C 2022; CYP3A4/5 2023) · the 2025
open-access star-allele-caller benchmark (88 1000G/GeT-RM samples, 9 VIP genes).
