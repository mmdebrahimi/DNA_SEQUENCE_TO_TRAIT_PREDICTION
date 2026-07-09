# A protein language model does not help where the HIV DRM catalog is blind

**Date:** 2026-07-09 · **Model:** ESM2-650M (the peak ESM2 checkpoint) · **CPU, free**
**Labels:** Stanford PhenoSense fold-change — the project's one free, independent, isolate-level,
**continuous** wet-lab genotype-phenotype label.
**Comparator:** the **deployed** rule `hiv_amr.call_hiv_observed` — verified to agree with a hand
re-implementation on **2168/2168** EFV isolates before any scoring was done (not a strawman).
**Scripts:** `scripts/hiv_esm_vs_catalog.py` (EFV + controls) · `scripts/hiv_esm_vs_catalog_allrt.py`
(all 11 RT drugs) · **Artifacts:** `wiki/hiv_esm_vs_catalog_2026-07-09.json`,
`wiki/hiv_esm_vs_catalog_allrt_2026-07-09.json`

## The question

The curated DRM catalog is deterministic and frozen, and it has a structural blind spot: variants it
does not list. On EFV at a 3× fold-change cutoff it reaches sensitivity 0.947 / specificity 0.904, but
**53 resistant isolates carry no catalog DRM at all** (19 of them ≥10× fold-change). A learned
variant-effect scorer is the obvious candidate to fill exactly that gap — and the project's recorded
regime boundary says *molecular-property → learned wins*.

So: **on the catalog-negative subset, does ESM2-650M separate resistant from susceptible?**

Pre-registered bar: **ESM AUROC ≥ 0.65 AND > mutation-burden baseline AND shuffled-null < 0.55.**

## Answer: no — and the positive controls say why

| | EFV, catalog-negative subset (n=1111; 53 R / 1058 S) |
|---|---|
| **ESM2-650M mean-damage AUROC** | **0.449** |
| mutation-burden baseline | 0.451 |
| shuffled-label null | 0.497 |

**VERDICT: FAIL.** ESM is indistinguishable from chance, and from counting mutations.

The two positive controls turn this from a shrug into a diagnosis:

- **Control A — full cohort (n=2168):** ESM AUROC **0.454**, *below chance*, while the deployed
  catalog scores **0.926**. ESM does not see NNRTI resistance anywhere, not just in the blind spot.
- **Control B — do the 16 real DRMs look damaging?** Median damage percentile at their own position:
  **0.29**. ESM rates the true resistance mutations as **less** damaging than a random substitution
  there. Only Y181C (0.95), Y188H (0.84) and K103N (0.63) rank high.

**Mechanism.** ESM's masked-marginal is *evolutionary likelihood*, and DRMs are frequently observed:
they are viable, circulating, and drug-selected, and UniRef contains many drug-experienced HIV RT
sequences. So the model rates them **likely**, not damaging. The signal is not merely absent — it is
**anti-aligned** with resistance. Without Control A, the 0.449 would have been misread as "the blind
spot is irreducible" instead of "the instrument points the wrong way."

## It generalizes across 11 RT drugs, and the exception is instructive

Both NNRTI and NRTI target reverse transcriptase, so the cached masked-marginal matrix scores all of
them for free.

| drug | class | catalog AUROC | ESM full-cohort | ESM catalog-negative | burden |
|---|---|---|---|---|---|
| EFV | NNRTI | 0.926 | 0.454 | 0.449 | 0.451 |
| NVP | NNRTI | 0.948 | 0.464 | 0.463 | 0.533 |
| ETR | NNRTI | 0.729 | 0.491 | 0.387 | 0.341 |
| RPV | NNRTI | 0.683 | 0.430 | 0.439 | 0.400 |
| DOR | NNRTI | 0.549 | 0.572 | 0.796 ⚠ | 0.783 |
| 3TC | NRTI | 0.865 | **0.766** | 0.592 | 0.494 |
| ABC | NRTI | 0.809 | **0.734** | — underpowered | — |
| AZT | NRTI | 0.745 | **0.661** | — underpowered | — |
| D4T | NRTI | 0.668 | 0.598 | — underpowered | — |
| DDI | NRTI | 0.650 | 0.585 | — underpowered | — |
| TDF | NRTI | 0.656 | 0.567 | — underpowered | — |

**Genuine passes: 0 / 6.** ⚠ DOR clears 0.65 *nominally only*: n=37 with 21 positives, a burden
baseline of 0.783 (ESM's lift is **+0.013**), against a catalog that is itself near-chance (0.549).
A raw `≥0.65` count would have recorded that as a win; the machine-readable JSON now carries
`n_drugs_passing_bar_genuine: 0` alongside the nominal count, plus the caveat.

**ESM's sensitivity tracks the fitness cost of the resistance mechanism, not resistance.** On NNRTI
it is at or below chance (0.43–0.49): K103N and Y181C are cheap, high-frequency, well-tolerated. On
NRTI it is genuinely sensitive (0.57–0.77, 3TC best): M184V, K65R and the TAMs carry real replicative
cost, so evolutionary likelihood does see them. **But even where ESM is sensitive it loses to the
catalog** (3TC 0.766 vs 0.865; AZT 0.661 vs 0.745), and it adds nothing on the blind spot.

## What this settles

The decoder should **not** add a zero-shot LM variant scorer to fill the DRM catalog's blind spot.
The deterministic catalog stays the right product. This is a clean, free, independent-label negative.

It also **refines the recorded regime boundary** (`feedback_g2p_decoder_regime_boundary`). "Molecular
property → learned wins" was measured on DMS *fitness* (ProteinGym; ESM2-650M median Spearman 0.490).
Drug resistance is a different kind of molecular property: it is **antagonistically selected**. The
resistant variant is the one evolution has *kept*, so an evolutionary-likelihood model scores it as
ordinary. Learned scoring wins where the phenotype aligns with fitness and loses — even inverts —
where the phenotype is defined against a drug the training distribution already adapted to.

## Caveats

- **Uniform illustrative cutoff (fold ≥ 3×)** for every drug; no per-drug clinical cutoff exists
  in-repo for these classes. This is the precedent set by the PI v0.1 work: report the **comparison**,
  not an absolute-calibration claim.
- ESM damage is **length-normalized** (mean over mutations), so the score cannot grow with mutation
  count by construction; the sum variant is reported separately and is worse (0.432).
- 4 RT positions (122, 214, 272, 277) are excluded — HXB2 differs from consensus B there, so the
  wild-type would be wrong. Derived empirically: 314/318 positions have zero pure HXB2 calls.
- Only masked-marginal zero-shot scoring was tested. A **supervised** model trained on fold-change
  (the VespaG shape) is untested and is a different claim — but it would need these same labels to
  fit, so it cannot then be validated on them without leakage.
- AlphaMissense is not a comparator here: it covers the human proteome, not viral RT.
