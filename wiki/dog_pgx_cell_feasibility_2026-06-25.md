# Dog cells (MDR1/ABCB1 + coat-colour) — feasibility probe + WALL (2026-06-25)

Probed recs #2 (dog MDR1/ABCB1) + #3 (dog coat-colour) before building. **Verdict: BLOCKED — external/
data-access wall.** Building an unvalidated dog cell would violate the project's core discipline (validate
against a free independent per-sample label) + the "validate-wrapper-vs-tool" lesson, so it is NOT built.
This is the "labels-not-models" wall reasserting itself across the kingdom boundary.

## Why the human PGx pattern does NOT transfer cleanly to dog
The CYP2C19/CYP2C9 cells worked because of THREE free things: PharmVar/CPIC defs + 1000G genotypes +
**a free per-sample consensus panel (GeT-RM)** to score against. The dog cells lose the third — and hit two
more obstacles:

1. **The causal variant is an INDEL with ambiguous coordinate frames.** ABCB1-1Δ = nt230(del4) =
   c.230_233del (also written c.227_230 / c.296_299 / c.259_262 across groups). The SNP-based pgx caller
   does not handle indels; an authoritative canFam4 genomic coordinate is not cleanly exposed (OMIA CSV /
   Ensembl VEP needed; the NCBI CanFam3→ROS/canFam4 remap service is RETIRED). Genome-mode would need a dog
   reference + a BLAST/indel caller (a fungal/HIV-style build), not a SNP lookup.
2. **Free genotype source ships SNVs; the indel's coverage is unverified.** Dog10K's prominent posted file
   is the SNV VCF; indels had no truth set + weaker filtering ("verify whether the nt230del4 indel survived
   the non-SNV filters"). Darwin's Ark is low-coverage IMPUTED (indels imputed unreliably; no curated
   per-dog MDR1 table).
3. **No free per-dog genotype↔MEASURED-phenotype panel (the GeT-RM-equivalent does not exist for dog).**
   Dog phenotype is BREED-associated, not individually measured-and-free. The best achievable free
   validation is a breed-carrier-frequency sanity check — a strictly weaker tier than the human cells'
   per-sample independent concordance (CYP2C19 72/72, CYP2C9 73/73). The data-source research already
   flagged this risk ("per-genome phenotype labels are thinner; the validation denominator is the gate");
   the probe confirms it for MDR1 specifically.

Coat-colour (#3) is the same shape but worse: multiple genes (MC1R/ASIP/TYRP1/MFSD12), several non-SNV, and
the phenotype (coat colour) is owner-reported / breed-level, not a free per-dog measured label.

## Wall classification
- **External / data-access** (not a quick code win). To unblock, a USER decision on dependency-acceptance:
  (a) build a dog-reference + BLAST/indel caller (fungal/HIV pattern) — code-closable but a real build; AND
  (b) accept a WEAK validation tier (breed-carrier-frequency sanity) OR self-curate a per-dog
  genotype+phenotype set (e.g. extract chr14 ABCB1 from Darwin's Ark Dryad + join the owner survey — real
  effort, imputed-indel-cautious), since no free GeT-RM-equivalent exists for dog.
- **Honest tier ceiling:** even fully built, a dog OMIA cell is "faithful-to-OMIA-catalog + breed-frequency
  sanity" — NOT the independent per-sample concordance the human PGx cells achieve. That gap is the point.

## Decision: STOP at the wall (do not manufacture an unvalidated cell)
Recs #2/#3 are recorded as BLOCKED:external. The human-PGx expansion (recs #1, done: warfarin pair shipped
v0.6.2, CYP2C9 73/73) is the validated win. Further human PGx genes (NUDT15/TPMT/DPYD/SLCO1B1/VKORC1-done...)
remain the clean, free, GeT-RM-validatable path; the dog kingdom needs the user dependency-acceptance above.

_Sources: OMIA:001402-9615 (omia.variant:469); Mealey 2001; Dog10K (NAR 2025; Zenodo 8084059) SNV-VCF +
non-SNV filter note; Darwin's Ark (Dryad doi:10.5061/dryad.g4f4qrfr0, imputed autosomal). Probe 2026-06-25._
