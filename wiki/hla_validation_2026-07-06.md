# HLA tag-SNP validation vs the free 1000G HLA truth (2026-07-06)

**The wrapper-vs-truth number that gates the `dna-hla` cell.** The tag SNP → HLA-allele call is an LD PROXY,
so its validity had to be measured against a REAL HLA truth set — not asserted from the literature LD. Found
the free 1000G HLA truth (`ftp.1000genomes.ebi.ac.uk/.../20140725_hla_genotypes/20140702_hla_diversity.txt`,
1277 samples, 4-digit A/B/C/DRB1/DQB1) and joined it to the rs-tag genotypes on the 1000G panel.

## Result — only 1 of 3 provisional tags is deployable

| allele / drug | tag SNP | n scored | sens | spec | PPV | verdict |
|---|---|---|---|---|---|---|
| **HLA-B\*57:01 / abacavir** | rs2395029 | 1103 | **0.979** | **0.992** | **0.855** | **VALIDATED — deployed** (matches the clinical screen) |
| HLA-B\*58:01 / allopurinol | rs9263726 | 1103 | 0.609 | 0.824 | 0.176 | WEAK proxy (mixed-population LD; misses 39% of carriers) — **demoted** |
| HLA-A\*31:01 / carbamazepine | rs1061235 | 1103 | 0.000 | 1.000 | — | tag NOT paneled on 1000G (0 TP / 74 FN) — **demoted** |

## The lesson (honest, load-bearing)

A single-SNP LD proxy that is near-perfect for **B\*57:01** (rs2395029, the deployed abacavir screen: 47 TP /
1 FN / 8 FP) does **NOT generalize** to other HLA drug-risk alleles:
- **B\*58:01** — rs9263726's LD to B\*58:01 is strong only in East Asians; on the mixed 1000G panel it is a
  weak proxy (sens 0.61, PPV 0.18), unsafe for an SJS/TEN screen where a missed carrier is severe.
- **A\*31:01** — rs1061235 is not even present on the 1000G panel at chr6:29945521; A\*31:01 has no clean
  single-SNP tag → it needs sequence-based typing (HLA\*LA / arcasHLA / OptiType).

So the shipped `dna-hla` cell is **B\*57:01/abacavir only** — the one tag that cleared real validation. The
other two are kept as a documented negative (`dna_decode.hla.catalog._UNVALIDATED_TAGS`), NOT routable cells
(a clinical screen that silently mis-calls would be worse than none). This is the
validate-the-wrapper-against-the-tool-on-independent-data discipline: in-cohort tag frequency looked fine for
all three, but only the sample-level concordance vs real HLA types exposed which tags actually work.

## Provenance
Truth: 1000G `20140702_hla_diversity` (committed `tests/data/pgx_getrm/1000g_hla_diversity_20140702.txt`).
Per-allele carrier truth TSVs + per-allele result JSONs (`wiki/hla_{b5701,b5801,a3101}_validation_2026-07-06.json`).
Harness: `scripts/hla_concordance.py --allele <k> --vcf <region> --truth <tsv>`. NOT a clinical tool.
