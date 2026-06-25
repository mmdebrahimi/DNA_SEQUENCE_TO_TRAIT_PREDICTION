# CYP2C19 validation fixtures — PharmCAT test VCFs

These VCF files are vendored verbatim from the PharmCAT project's test resources
(`src/test/resources/org/pharmgkb/pharmcat/haplotype/cyp2c19/`), used here as an
INDEPENDENT (we-did-not-author) cohort to validate `dna_decode/pgx/` against the
reference tool's own expected diplotypes (the filename encodes the expected call,
e.g. `s1s2.vcf` -> *1/*2). This is a FAITHFUL-TO-PHARMCAT, in-distribution number —
NOT the GeT-RM independent panel (that needs a 1000 Genomes VCF fetch via tabix).

- Source: https://github.com/PharmGKB/PharmCAT  (License: Mozilla Public License 2.0)
- Fetched: 2026-06-25 from `main` via raw.githubusercontent.com
- Used by: scripts/pgx_cyp2c19_validate.py + tests/test_pgx_validate.py
