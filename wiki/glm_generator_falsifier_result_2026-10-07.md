# The GLM falsifier: first real numbers, and a measured defect in my own instrument

**Date:** 2026-10-07 · **Artifact:** `wiki/glm_generator_falsifier_2026-10-07.json`
**Reproduce:** `uv run python scripts/glm_generator_falsifier.py --n 500 --limit-natural 2400 --modes kmer,both --bins 6`
(no network, no GPU, no labels — the committed MG1655 reference is the only input)

## What was built

The measurement apparatus for GLM Family 2, validated against **real** *E. coli* regulatory sequence
before any model weights were downloaded. The natural corpus is the 150-bp window 5-prime of every
annotated CDS in `GCF_000005845.2` — thousands of genuine promoters, free, already on disk.

**The split that keeps it honest:** the corpus is halved. The Markov null is fitted on the FIT half; every
discrimination runs on the EVAL half. No generator is ever scored against sequence it was trained on.

## The three numbers, and why the first two come first

`distinguishability` is a held-out-AUROC discriminator between natural and generated sequence, and
**LOWER IS BETTER** — 0.5 means indistinguishable from real DNA, 1.0 means trivially detectable. It is
only trustworthy behind two controls:

| | natural-vs-natural (unbiased?) | natural-vs-uniform (powered?) | **Markov null** |
|---|---|---|---|
| **k-mer** features (k=4) | 0.5060 *(from 0.5053/0.4909/0.5048/0.5060/0.5004)* | 0.9578 | **0.5714** |
| **k-mer + positional** (6 bins) | 0.5442 | 0.9694 | **0.6737** |

Both controls behave as they must: the self-control centres on chance, so the instrument is not reacting
to the split; the uniform control is high, so it has power. **The verdict function refuses outright when
either fails**, and that refusal fired for real during this work (below) rather than only in a test.

## Finding 1 — an order-3 Markov chain is nearly indistinguishable from real promoters

Under k-mer features the null scores **0.5714 against a derived noise floor of 0.5060**. Its GC matches to
three decimals (0.4521 vs 0.4528). A chain with no learning, no weights and no model reproduces the local
composition of real *E. coli* promoters almost exactly.

**That is a ceiling effect in the instrument, not a compliment to the null.** It leaves ~0.065 of AUROC
between "as good as a trivial baseline" and "indistinguishable from reality" — far too little room to rank
a real generator in. Had I run GENERator against this and got 0.55, the number would have been
uninterpretable.

This is the same shape as results already measured in this repo: a single integer recovers 75% of a
protein LM's headline; a 1992 substitution matrix ties ESM2 on resistance ranking. **A metric a dumb
baseline saturates cannot grade a smart model.**

## Finding 2 — position-aware features break the ceiling, and the prediction was registered first

Pre-registered hypothesis: a real promoter has *architecture* (−35 and −10 boxes at characteristic
distances from the start), while a Markov chain is **translation-invariant** and has no notion of position
— so positional features should separate the null decisively, with the self-control staying near chance.

**Confirmed, and replicated at two sample sizes:**

| run | k-mer | k-mer+positional | Δ null | Δ ceiling | net |
|---|---|---|---|---|---|
| n_nat=1600 | 0.5902 | 0.7142 | **+0.124** | +0.028 | ≈ +0.096 |
| n_nat=2400 | 0.5714 | 0.6737 | **+0.102** | +0.038 | ≈ +0.064 |

The ceiling is **re-derived per feature mode**, because a richer feature set has more capacity to fit split
noise — reusing the k-mer ceiling would have understated the noise and inflated the finding. Subtracting
it, the real gain is ≈ +0.06 to +0.10 of AUROC.

**Reading:** the null reproduces *composition* and fails at *architecture*. Usable headroom goes from
**~0.065 → ~0.13** of AUROC, so the instrument can now actually rank a generator.

This independently reproduces, on this project's own data, the published audits that separate generated
from natural sequence while naming **long-range organisation** as the failure mode. Those audits were the
reason the design's original bar ("constraint-pass above a Markov baseline") was rejected as too weak; this
measures *why* it was too weak.

## Finding 3 — the refusal gate fired on real data, not just in a test

A positional-**only** run (per-bin dinucleotide, 160 features) reached **0.8893** on the power control,
below the 0.90 bar, and the falsifier returned `INSTRUMENT_UNDERPOWERED` with **no distinguishability
number at all**. Binning dilutes counts (~25 bp per bin), so it lost the k-mer power it needed.

That is the gate working as designed. The alternative would have been reporting a flattering number from
an instrument that could not separate i.i.d. noise from real DNA — and **a discriminator AUROC of 0.52
reads as "excellent, indistinguishable" while being equally consistent with "the classifier learned
nothing."** My hypothesis was *untested* at that point, not refuted; `mode="both"` is what tested it.

## Honest limits

- **An upstream window is a PROXY for a promoter**, not a mapped transcription start site. 2708 of 4340
  CDS windows overlap another CDS — real *E. coli* compactness and operons, measured and reported, not
  silently dropped.
- **One organism.** Nothing here transfers to another host without being re-measured.
- **A k-mer/positional-bin discriminator still only sees coarse structure.** Low distinguishability is
  **necessary, not sufficient** — the published audits' failure axis is longer-range than 6 bins.
- **`dna_decode/constraints/` is CDS-shaped** (frame, internal stop, codon reachability) and therefore
  **INAPPLICABLE** to non-coding sequence. It is not run, and no vacuous 100%-pass is reported; alphabet
  validity and exact length are the checks that do apply.
- **No learned generator has been scored yet.** GENERator-v2-prokaryote (MIT, 4.33 GiB, verified against
  the HF API) is wired behind `--model auto` but not downloaded. **The bar it must beat is now a real,
  replicated number: 0.6737 under `both`, against a floor of 0.5442.**
- The self-control ceiling is **sample-size dependent** (0.5447 at n=1200 → 0.5060 at n=2400 under
  k-mer), which is why it is derived per run rather than pinned as a constant.

## What this changes

The next GLM measurement has a defensible instrument and a real bar. Before this, "does the generator
produce DNA-like sequence?" had no answer that a 3-mer Markov chain could not also claim.
