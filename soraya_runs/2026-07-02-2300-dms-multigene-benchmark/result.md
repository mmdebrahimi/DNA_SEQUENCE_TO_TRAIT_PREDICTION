# Soraya --until-mvp (overnight): DMS variant-effect benchmark + AlphaMissense contrast

Verdict: mvp-reached (both deliverables shipped). 2 commits: d61a86c + 59b9083.

## What shipped tonight (the "many gene samples" + push-forward)
1. **24-assay multi-gene DMS BLOSUM benchmark** (17 human genes; abundance/function/binding; ~180k missense).
   Deterministic BLOSUM62 floor = 0.219 median polarity-corrected Spearman; modality-dependent
   (abundance 0.243 > function 0.163). scripts/dms_variant_effect_benchmark.py.
2. **AlphaMissense measured contrast** (13 assays / 11 genes; offset-corrected join, match_rate 1.0):
   AlphaMissense median 0.515 vs BLOSUM 0.240 = +0.265 (~2.1x). The strong precomputed predictor (free, no
   GPU) ~doubles the deterministic floor -- strongest where BLOSUM is weakest (function assays).
   scripts/dms_alphamissense_benchmark.py.

## The result (honest framing)
Together these give the project its FIRST hard-numbers map of the molecular-phenotype G2P frontier + identify
AlphaMissense as the deployable strong tool. NOT a "deterministic breakthrough" (AM is a learned/AF-distilled
predictor) -- but a genuine, systematic, MEASURED contribution across 11-17 human genes that sharpens the
deterministic-vs-learned boundary with real numbers. The 3-regime law is now measured, not just argued.

## Discipline
verify-in-batch caught 4 real bugs unattended (min_n test threshold; AM sorted-by-uniprot not genome-pos ->
15min stream died in the A's; CRLF in the grep pattern file -> 0 matches; polarity sign alignment). 5 offline
tests across both cells; frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9 throughout).

## Banked here (not further) -- why
Unattended overnight + a coherent measured result reached; each of the 4 caught bugs needed live attention, so
further cells risk quality for diminishing returns (plateau + Care). Queued for a fresh attended pass: a
deterministic conservation feature (needs per-protein MSAs) that could CLOSE the gap = the real deterministic
breakthrough; PGx v0.2 using AM as a functional prior for withheld non-core alleles; the parked PGx-gene cell.

No-resume: bounded attempts per active session only.
