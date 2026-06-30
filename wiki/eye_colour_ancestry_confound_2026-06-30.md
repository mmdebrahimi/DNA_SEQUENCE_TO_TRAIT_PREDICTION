# Eye-colour M2 — rs12913832 ancestry-confound, quantified (D:-free first cut, 2026-06-30)

The eye-colour v0 (acc 0.993) + v0.1 cells carry a standing caveat: "rs12913832 is European-calibrated →
ancestry-confounded." With the D: drive disconnected (the OpenSNP genotypes are offline), the *within-cohort*
re-score is blocked — but the confound's **structural magnitude** can be quantified from public 1000 Genomes
population frequencies (Ensembl REST, no D:, no big download). This is the M2 increment that was achievable
while blocked.

## The number

Blue allele (G) frequency by 1000G phase-3 super-population (sourced: Ensembl REST `variation/human/
rs12913832?pops=1`, fetched 2026-06-30; committed in `dna_decode/data/eye_colour_ancestry.py`):

| Super-population | Blue allele (G) freq |
|---|---|
| **EUR** (European) | **0.636** |
| AMR (admixed American) | 0.202 |
| SAS (South Asian) | 0.071 |
| AFR (African) | 0.028 |
| **EAS** (East Asian) | **0.002** |
| ALL | 0.177 |

**EUR / EAS ratio = 318×.** The blue allele is overwhelmingly European-concentrated; it is essentially
absent in East Asians and rare in Africans.

## What this does and does NOT establish (the honesty rail)

**Does establish:** rs12913832 is a strong **ancestry-informative** marker. OpenSNP is a Europe-/US-majority
DTC cohort, so a user's rs12913832 genotype is *partly a tag for European ancestry* — and European ancestry
itself correlates with light (blue) eyes. So a portion of the v0 0.993 accuracy *could* be the SNP tracking
ancestry rather than the eye-colour mechanism. This bounds the confound's plausibility: it is real and
structurally large.

**Does NOT establish:** that the 0.993 is inflated. rs12913832 is **also a known CAUSAL variant** — a HERC2
intronic regulatory element controlling OCA2 expression, the dominant eye-colour locus. High within-European
predictive accuracy is the established biology, not merely an ancestry artifact. Ancestry-informativeness and
causality co-exist here.

**The disentangler (deferred, D:-gated):** the only clean resolution is the **within-European re-score** —
does rs12913832 predict eye colour *among European-ancestry users only*? If accuracy stays high within EUR,
the signal is mechanistic, not an ancestry tag. That requires the per-user OpenSNP genotypes + an ancestry
axis (self-reported ancestry column and/or genome-wide ancestry inference), both on the disconnected D:
drive. Queued for when D: is reconnected.

## Status
- Artifact: committed frequencies + `confound_summary()` helper + 3 tests (`tests/test_eye_colour_ancestry.py`).
- This is **M2 partial** (structural magnitude). M2 full (within-ancestry re-score), M1 real IrisPlex number,
  M3 (PGP), M4 (Mendelian/PGx) remain blocked on the D: reconnect.
- FROZEN AMR surface byte-unchanged. Reproduce: `confound_summary()` (offline) or re-fetch via the Ensembl
  REST URL above.
