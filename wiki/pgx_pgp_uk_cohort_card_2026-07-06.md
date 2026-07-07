# PGx real-people COHORT card — 5 PGP-UK individuals (2026-07-06)

_deployment/robustness demonstration on real independent-cohort individuals; tiny-N observed allele counts (NOT a population estimate); no GeT-RM truth here (accuracy lives in the concordance cells). NOT a clinical tool._

**Cohort:** PGP-UK (Personal Genome Project UK), ENA PRJEB17529, open-consent GRCh37 · **N = 5** · samples: FR07961000, FR07961003, FR07961006, FR07961007, FR07961009

## Per-individual calls

| sample | CYP2C19 | CYP2C9 | CYP2D6 | CYP3A5 | TPMT | CYP2B6 | VKORC1 | SLCO1B1 |
|---|---|---|---|---|---|---|---|---|
| FR07961000 | *1/*2 IM | *1/*1 NM | *1/*10 NM | *1/*6 IM | *1/*1 NM | *1/*6 IM | G/G | T/T |
| FR07961003 | *2/*17 IM | *1/*1 NM | *1/*1 NM | *1/*1 NM | *1/*1 NM | *1/*6 IM | G/A | T/T |
| FR07961006 | *1/*2 IM | *1/*2 IM | *1/*1 NM | *1/*1 NM | *1/*1 NM | *1/*6 IM | A/A | T/T |
| FR07961007 | *1/*17 RM | *1/*1 NM | *1/*41 NM | *1/*1 NM | *1/*1 NM | *1/*6 IM | A/A | T/T |
| FR07961009 | *1/*1 NM | *2/*2 IM | *1/*1 NM | *1/*1 NM | *1/*1 NM | *1/*1 NM | G/G | T/T |

## Per-gene distribution (across the cohort)

- **CYP2C19** — phenotypes: {'IM': 3, 'RM': 1, 'NM': 1}; diplotypes: {'*1/*2': 2, '*2/*17': 1, '*1/*17': 1, '*1/*1': 1}; observed allele counts: {'*1': 5, '*2': 3, '*17': 2}
- **CYP2C9** — phenotypes: {'NM': 3, 'IM': 2}; diplotypes: {'*1/*1': 3, '*1/*2': 1, '*2/*2': 1}; observed allele counts: {'*1': 7, '*2': 3}
- **CYP2D6** — phenotypes: {'NM': 5}; diplotypes: {'*1/*1': 3, '*1/*10': 1, '*1/*41': 1}; observed allele counts: {'*1': 8, '*10': 1, '*41': 1}
- **CYP3A5** — phenotypes: {'NM': 4, 'IM': 1}; diplotypes: {'*1/*1': 4, '*1/*6': 1}; observed allele counts: {'*1': 9, '*6': 1}
- **TPMT** — phenotypes: {'NM': 5}; diplotypes: {'*1/*1': 5}; observed allele counts: {'*1': 10}
- **CYP2B6** — phenotypes: {'IM': 4, 'NM': 1}; diplotypes: {'*1/*6': 4, '*1/*1': 1}; observed allele counts: {'*1': 6, '*6': 4}
- **VKORC1** — genotypes: {'G/G': 2, 'A/A': 2, 'G/A': 1}
- **SLCO1B1** — genotypes: {'T/T': 5}

_**CYP2D6 caveat (load-bearing honesty):** the CYP2D6 call here is a SNP-proxy diplotype from a called VCF — it CANNOT see the structural alleles (*5 deletion / *xN duplication / *13/*36/*68 hybrids), so every CYP2D6 cell carries `cnv_hybrid_unassessed`. The copy-number half is resolvable only from a BAM/CRAM (dna_decode.pgx.cyp2d6_structural, 26/26 on 1000G CRAMs); PGP-UK ships VCFs, not reads._

_Accuracy-vs-truth lives in the GeT-RM concordance cells (`wiki/pgx_report_card.md`); this card is real-world deployment coverage, not a new accuracy number._
