# The 3-way ESM2+ProSST+GEMME — does adding evolution beat the 2-way? (mostly no; phenotype-conditional)

**Date:** 2026-07-18 · **Components:** OUR ESM2-650M tables + OUR ProSST-2048 (validated, reproduces the
column at 1.0) + ProteinGym's precomputed **GEMME** column · **Script:** `scripts/three_way_lift.py` ·
**Data:** `data/processed/three_way_lift_checkpoint.jsonl` · **N=56, LOCAL CPU**

## Question + the framing reframe (R2)

The modality-hybrid sweep's TOP combination was `ESM2+GEMME+ProSST` (0.547). The directive was "move forward
with the 3-way — needs the GEMME install." **An R2 pre-build check flipped that:** the 3-way VALIDATION does
NOT need the (Windows-hostile JET2/R/Java) GEMME install, because (a) `rank_average_hybrid` is already N-ary,
and (b) **GEMME has ZERO learned parameters** — a deterministic evolutionary-conservation model — so
ProteinGym's precomputed GEMME column IS canonical GEMME output, exactly like the pre-quantized ProSST
structures used to reproduce the ProSST column at 1.0. The install is needed ONLY to run GEMME on a NOVEL
protein (deployment), not for this validation. So the 3-way was validated with our-ESM2 ⊕ our-ProSST ⊕
canonical-GEMME — no install.

## Result — the 3-way does NOT beat the 2-way overall; the extra lift is phenotype-conditional

| comparison | median Δ | win-rate | sign-p |
|---|---:|---:|---:|
| 3-way vs **2-way (ESM2+ProSST)** | **+0.0035** | **31/56 (55%)** | **0.50 (n.s.)** |
| 3-way vs ESM2 baseline | +0.0709 | 53/56 (95%) | 8e-13 |

median |Spearman|: 3-way **0.599** vs 2-way **0.594** vs ESM2 0.515 (GEMME-alone 0.482). **Adding GEMME on
top of ESM2+ProSST buys essentially nothing on average (+0.005 Spearman, a coin-flip win-rate).** The
already-validated 2-way structure hybrid captures the lift; the third (evolution) modality is largely
redundant with it.

**But the aggregate ≈0 is a category-mix CANCELLATION, not "GEMME is useless"** — per phenotype (3-way − 2-way):

| phenotype | n | median Δ | win-rate |
|---|---:|---:|---:|
| Activity | 9 | **+0.0180** | 7/9 |
| OrganismalFitness | 12 | **+0.0145** | 9/12 |
| Expression | 9 | +0.0023 | 5/9 |
| Binding | 4 | −0.0033 | 2/4 |
| **Stability** | 22 | **−0.0109** | 8/22 |

GEMME (evolution) **helps** exactly the phenotypes where evolution is the informative modality
(Activity/function +0.018, OrganismalFitness +0.015) and **hurts** Stability (−0.011, only 8/22) where
structure already dominates and evolution is redundant. This is the SAME per-category modality rule the 2-way
work found (evolution→activity, structure→stability), now one layer up: a *third* modality only pays off
where its signal is the missing one, and Stability is the largest bucket (22/56), so it drags the aggregate
to ~0.

## Deployable conclusion

- **The 2-way `ESM2+ProSST` is the sweet spot for most cells** — it is the validated, powered lift (+0.067,
  93%, `wiki/prosst_lift_2026-07-18.md`); adding GEMME does not improve it overall.
- **Add GEMME (go 3-way) ONLY for evolution-favorable phenotypes** — activity/function/organismal-fitness
  cells (+0.015–0.018), never stability cells (where it hurts). Phenotype-conditional routing, again.
- The sweep's +0.017 (0.547 vs 0.530) 3-way-over-2-way margin shrinks to +0.005 in our pipeline (our own
  ESM2 vs ProteinGym's ESM2 column); the direction (tiny positive) matches, the magnitude does not clear
  significance here.

## Honest scope

- The MVP bar "3-way beats 2-way, median Δ > 0" is technically met (+0.0035 > 0) but the honest headline is
  **not-significant-overall + phenotype-conditional** — reported as such rather than as a clean win.
- The GEMME component is ProteinGym's PRECOMPUTED column (canonical, since GEMME is deterministic), NOT our
  own GEMME forward. For a novel protein the 3-way needs a real GEMME run (`gemme_scorer.run_gemme`, MSA from
  `msa_fetch` + the JET2/R/Java toolchain — Windows-hostile, deferred, the same class as ProSST's quantizer).
- N=56 (structure ∩ ESM2-table ∩ GEMME-column ∩ seq_len≤400).

## Shipped

- `dna_decode/forward/gemme_scorer.py` — `gemme_table_from_column` (the deployable-now column adapter) +
  `run_gemme` (deferred, raises `GemmeUnavailable` without the toolchain) + `gemme_tier`.
- `predict_effect(method="gemme", gemme_table=…)` + package exports; the N-ary `rank_average_hybrid` needed
  no change (3-way already composed).
- `scripts/three_way_lift.py` — the 3-way-vs-2-way harness.
- Frozen decoder surface byte-unchanged (`verify_lock OK`).
