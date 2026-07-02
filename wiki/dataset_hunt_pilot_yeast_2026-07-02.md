# F4 pilot — 1002 Yeast Genomes (real-surface-first go/no-go, 2026-07-02)

The #1 shortlisted substrate, pilot-fetched + join-confirmed + signal-checked on REAL data. **Verdict: GO.**

## What ran
- Fetched (free, http, no DUA) from `1002genomes.u-strasbg.fr/files/`:
  - `phenoMatrix_35ConditionsNormalizedByYPD.tab.gz` (176 KB) → **971 isolates × 35 lab growth conditions**
    (YPD-normalized ratios; e.g. YPDFLUCONAZOLE, YPGALACTOSE, YPDNACL1M …).
  - `genesMatrix_PresenceAbsence.tab.gz` (408 KB) → **1011 isolates × 7,796 genes** (0/1/NA).
- **Join confirmed:** both keyed by the same 3-letter isolate code → **971 isolates present in BOTH** matrices.
  Per-unit genotype ⋈ per-unit phenotype works (G4 satisfied on real data, not just claimed).

## Signal-vs-null (honest, permutation-controlled)
Condition = YPDFLUCONAZOLE (antifungal growth). For each variable gene (≥20 isolates per presence class), a
two-sample |t| of growth (present vs absent). Real vs 20 label-shuffled permutations:

| | max |t| | #genes |t|>4 | #genes |t|>6 |
|---|---|---|---|
| **REAL** | **7.63** | 90 | 23 |
| **NULL** (20 shuffles) | mean 3.53 / ceiling 5.50 | — | — |

**Real max|t| 7.63 > null ceiling 5.50 → genotype→phenotype signal is above the permutation null.** The
substrate carries real signal; the join is real; the data is free and deep. GO for a full build.

## Honest caveats (load-bearing — the project's own hard lessons apply)
1. **This is NOT a de-confounded test.** The marginal gene-presence association can partly be POPULATION
   STRUCTURE (yeast clades), exactly the failure mode that killed the embedding arm (cipro within-lineage,
   Arabidopsis flowering-time). A real learned/deterministic decoder here MUST de-confound by lineage
   (the 1011-genome clade structure / the provided SNP distance matrix) and check within-clade concordance
   before any "it works" claim. The pilot proves *a substrate exists*, not *that a decoder generalizes*.
2. Pilot used gene presence/absence (light genotype); the full SNP matrix (`1011Matrix.gvcf.gz`) is the real
   genotype surface for a build.
3. Single condition scanned; 34 others available. Fluconazole chosen for antifungal-relevance continuity
   with the C. auris cell.

## Where this leaves the hunt
- **Top substrate = 1002 Yeast Genomes** (PASS, learned-niche, depth 1011) — a NEW kingdom (fungi, whole-
  organism growth phenotype) that clears all 8 gates on verified real data. First non-AMR, non-pathogen
  substrate to pass the bar since the freeze.
- Runners-up (WEB-VERIFIED PASS): **DGRP2** (fly, 205 lines). VERIFY tier: CaeNDR (worm), ClinVar
  (human, deterministic), Mouse Phenome DB, Rice 3000 — each with one named gate to confirm.
- Arabidopsis flowering-time correctly REJECTED (closed embedding negative — do not re-run).

## Recommended next step (a real decoder build, gated on user GO)
Build a lineage-DE-CONFOUNDED yeast growth decoder on the 1002 substrate: full SNP matrix → a growth-trait
model → the domain-knowledge + within-clade concordance test (the embedding-niche 3-part bar). If it beats a
clade-only baseline within-clade, it's the first learned-decoder win; if not, it joins the honest negative map
as the 5th de-confounded embedding failure — either way a real result on a clean substrate.
