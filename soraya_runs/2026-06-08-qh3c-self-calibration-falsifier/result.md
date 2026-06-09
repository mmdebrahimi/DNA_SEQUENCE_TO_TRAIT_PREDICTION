# Soraya --advance — result

- **Run:** 2026-06-08-qh3c-self-calibration-falsifier
- **Family:** eukaryotic-trait-decoding-cycle-2026-06-07
- **Stop reason:** batch complete (5/5 actions) · commit 46a716f
- **Plateau check:** NOT a plateau — produced two new conclusions (one strand confirmed, one
  killed-then-repaired with a concrete design constraint).

## Ran
Dual falsifier for the two cheapest HIGH `/hypothesise` strands, on the 8 cached AMRFinder cohorts
already on disk (pure logic; no Docker / downloads / GPU / money).

## Results
- **#1 auto-threshold — SURVIVES.** LOO threshold sweep recovered Campylobacter→1 and Klebsiella→2
  (2/2 ground-truth agreement, unanimous picks, LOO acc 1.000 both).
- **#2 auto-intrinsic-flag — KILLED as specified, REPAIRED at gene-family granularity.** Per-symbol
  flags nothing (OXA-51-family fragments across alleles); family-level flags blaOXA-51-family + blaADC
  (both genuine intrinsics), spares all acquired carbapenemases.

## Forecast calibration
- Predicted: both strands cheap-and-informative. Outcome: #1 hit (survives), #2 partial (killed as
  literally specified, repaired). Net = strong VOI, no wasted motion.

## Next-VOI (not auto-run)
1. **Build `calibrate_organism(cohort)`** — wire auto-threshold + family-level auto-intrinsic into a real
   routine, with regression pins that auto-recovered values match hand-validated E. coli / Klebsiella /
   Campylobacter before overriding a default. (HIGH; weeks; the confirmed strand.)
2. **Strand #3 IS-element-upstream** — the biggest blind-spot strand; needs ISAba1 BLAST vs cached
   Acinetobacter assemblies (Docker). (HIGH; weeks.)
3. **Strands #5/#6 kingdom breadth** (antimalarial kelch13 / viral DRM) — engine transfer, free data.
4. Standing user-gated: Path B G2 GPU run on the workhorse.
