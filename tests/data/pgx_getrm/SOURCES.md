# GeT-RM CYP2C19 consensus truth set

`star-allele-comparison_common.tsv` is vendored from the ursaPGx benchmark
(coriell-research/ursaPGx, scripts/results/), column `CYP2C19_getrm_ngs` = the
GeT-RM NGS consensus diplotype (consensus of Astrolabe + Stargazer + Aldy;
Gaedigk et al. 2022, J Mol Diagn) for the ~87 samples overlapping the 1000 Genomes
30x panel by Coriell ID. Used as the INDEPENDENT truth set to validate
`dna_decode/pgx/` star-allele CALLING (our caller is independent of those 3 tools).

- Source repo: https://github.com/coriell-research/ursaPGx
- Used by: scripts/pgx_getrm_concordance.py + tests/test_pgx_getrm.py
- Genotypes: the public 1000 Genomes 30x phased VCF (fetched via Docker bcftools; gitignored, re-fetchable).
