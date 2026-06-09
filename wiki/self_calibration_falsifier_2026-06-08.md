# Self-calibrating-tool falsifiers — 2026-06-08

> Two grey-zone strands from `/hypothesise` (2026-06-08), falsified on cached AMRFinder runs already on
> disk. Pure logic — no Docker, no downloads, no money, no GPU. Strand = "the tool ports itself to a new
> organism by auto-calibrating from a small labeled cohort, instead of hand-set thresholds/curation."
> Falsifier code + raw JSON: `soraya_runs/2026-06-08-qh3c-self-calibration-falsifier/`.

## Strand #1 — AUTO-THRESHOLD — **SURVIVES**

**Claim:** leave-one-out balanced-accuracy threshold selection over cipro QRDR-point counts recovers the
biologically-correct per-organism threshold with no hand-tuning.

**Test:** LOO threshold sweep (grid {1,2,3}) on three cipro cohorts with direct cached `main.tsv`,
processed through the canonical `qrdr_point_determinants`. Ground truth: Campylobacter→1 (single gyrA
T86I, proven 2026-06-08), Klebsiella→2 (E. coli-style double-mutant, historically acc 1.0).

| cohort | N | LOO modal threshold | LOO picks | LOO acc | ground truth | verdict |
|---|---:|---:|---|---:|---:|---|
| campylobacter_ciprofloxacin | 30 | **1** | {1: 30} | 1.000 | 1 | **AGREE** |
| klebsiella_cipro | 30 | **2** | {2: 30} | 1.000 | 2 | **AGREE** |
| pseudomonas_aeruginosa_ciprofloxacin | 30 | 1 | {1: 30} | 0.900 | (none) | exploratory |

**Verdict: SURVIVES** — 2/2 ground-truth agreement, **unanimous** per-strain LOO picks (zero instability),
LOO acc 1.000 on both. The kill condition (disagree on ≥2 GT organisms) is nowhere near triggered. The
meta-rule recovers BOTH threshold regimes that motivated the hypothesis — the single-mutation (Campylobacter)
and double-mutation (Klebsiella) cases — from the counts alone.
*(Pseudomonas LOO acc 0.900 at threshold 1 is exploratory — no validated ground truth; flagged for a
dedicated look, not part of the kill test.)*

## Strand #2 — AUTO-INTRINSIC-FLAG — **KILLED as specified, REPAIRED at family level**

**Claim:** a determinant present in ≥90% of BOTH R and S strains can be auto-flagged as intrinsic
(reproducing the manual OXA-51-family exclusion) without literature curation.

**As specified (per AMRFinder SYMBOL): KILLED.** Nothing was flagged. Cause diagnosed by the falsifier:
*A. baumannii*'s intrinsic blaOXA-51-family is universal at the **family** level but fragments across many
distinct allele symbols (blaOXA-68 / -69 / -66 / -100 / -312 / ADC alleles…) — each strain carries a
*different* allele, so no single symbol reaches 90% prevalence. The per-symbol heuristic structurally
cannot see a polymorphic intrinsic.

**Repaired (collapse alleles to gene FAMILY): SURVIVES.** Re-running with `blaOXA-51-family` grouping +
`-<number>` allele-suffix stripping:

| family | prev_R | prev_S | intrinsic_flag |
|---|---:|---:|---|
| blaADC | 1.000 | 1.000 | **True** |
| blaOXA-51-family | 0.933 | 1.000 | **True** |
| blaOXA (acquired OXA, mixed) | 0.800 | 0.600 | False |
| (all strong acquired: OXA-23/24/72, NDM, IMP, VIM, KPC) | — | — | False |

Both flags are TRUE positives: blaOXA-51-family is the known intrinsic carbapenemase, and **blaADC is the
intrinsic *A. baumannii* AmpC cephalosporinase** (the heuristic surfaced a second genuine intrinsic the
hand-curation hadn't explicitly listed). Zero acquired carbapenemases falsely flagged. Kill condition
(intrinsic missed OR acquired falsely flagged) not triggered.

**Verdict: the strand is viable IFF intrinsic detection operates at gene-FAMILY granularity, not
per-allele symbols.** This is a concrete, cheap design constraint for the build — not a dead end.

## Combined takeaway

The "self-calibrating tool" thesis is **empirically supported on data already in the repo**: auto-threshold
works as-is; auto-intrinsic-flag works at family granularity. Together they would let the decoder absorb a
new organism from a ≥15R/15S labeled cohort — auto-picking the count threshold AND auto-excluding intrinsic
families — turning today's per-organism hand-curation (the Acinetobacter strength-tier, the Campylobacter
threshold) into an automated calibration step.

**NOT yet wired into the deployed rule.** This is falsifier evidence, not a shipped capability. Building it
= a `calibrate_organism(cohort)` routine + regression pins that the auto-recovered values match the
hand-validated ones on E. coli / Klebsiella / Campylobacter before it can override a default. Honesty
caveats: ground truth exists for only 2 organisms (auto-threshold); the family-grouping map is hand-built
(small, but it is curation); both tests are in-cohort on N=30. An independent per-organism cohort remains
the bar before production.
