# CYP2C19 caller on REAL 1000 Genomes (2026-06-25)

**Cohort:** 1000 Genomes 30x phased panel (chr10 CYP2C19 region, fetched via Docker bcftools)  (n=3202)

**Honesty tier:** REAL-1000G-DATA population characterization + 1 grounded GeT-RM data point (NA19122). NOT the full GeT-RM consensus concordance % (labels are in paper supplements -- a data-access step, not a tooling wall).

## Phenotype distribution (v0 caller, 3202 real genomes)

- Normal Metabolizer (NM): 1225  (38.3%)
- Intermediate Metabolizer (IM): 1038  (32.4%)
- Rapid Metabolizer (RM): 625  (19.5%)
- Poor Metabolizer (PM): 213  (6.7%)
- Ultrarapid Metabolizer (UM): 101  (3.2%)

## Diplotype distribution

- *1/*1: 1225
- *1/*2: 804
- *1/*17: 625
- *2/*17: 185
- *2/*2: 179
- *17/*17: 101
- *1/*3: 41
- *2/*3: 33
- *3/*17: 8
- *3/*3: 1

## Real-world blind-spot exposure (the v0.1 sentinel gap, quantified)

- *4-family carriers (rs28399504 ALT, invisible to v0): **5** (0.16%)
- of those, *4b mis-called as *17 (also carries the *17 SNP): **0**
- *35 mis-called as *1 (rs12769205 ALT without rs4244285): **44** (1.37%)

_These are the real-population rates at which the v0 core-SNP proxy would mis-call vs a full PharmVar caller -- the quantified cost of the v0.1 sentinel gap the brainstorm flagged. The v0.1 sentinel layer (rs28399504 + rs12769205) converts these silent mis-calls into withholds._

## Grounded GeT-RM check -- NA19122 (consensus *2/*35)

- v0 caller: ***1/*2** (Intermediate Metabolizer)
- GeT-RM consensus: ***2/*35**
- genotypes: rs4244285=1|0, rs12769205=1|1, rs28399504=0|0
- GeT-RM consensus *2/*35; v0 sees rs4244285 het -> calls *1/*2; the *35 haplotype (rs12769205 without rs4244285) is invisible to the v0 core set -> CONFIRMED *35 blind spot on the real genome.
