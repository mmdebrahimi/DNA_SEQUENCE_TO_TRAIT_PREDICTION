# GENERator-1.2B, measured: WORSE than a 3-mer Markov chain at producing promoter-like sequence

**Date:** 2026-10-07 · **Artifact:** `wiki/glm_generator_falsifier_REAL_2026-10-07.json`
**Reproduce:** `uv run python scripts/glm_generator_falsifier.py --model auto --device cpu --conditional
--context 300 --n 120 --limit-natural 1200 --modes kmer,both --bins 6 --natural-set promoter`
(~35 min on CPU at ~17 s/seq; weights 4.33 GiB from HF, MIT)

**This is the first learned genome generator this project has ever measured against real sequence**, and the
bar it was scored against was fixed and published **before** the run (`glm_generator_falsifier_result_2026-10-07.md`).

## The number

`distinguishability` = held-out AUROC of a classifier separating generated from real *E. coli* promoters.
**LOWER IS BETTER** — 0.5 is indistinguishable from reality.

| candidate | k-mer (k=4) | k-mer + positional (6 bins) |
|---|---|---|
| **Markov-k3 null** | **0.5990** | **0.6679** |
| **GENERator-v2-prokaryote-1.2b, conditional** | **0.7237** | **0.7706** |
| uniform (power control) | 0.9408 | 0.9641 |
| self-control ceiling (noise floor) | 0.5447 | 0.5865 |
| comparison tolerance (derived) | 0.0447 | 0.0865 |

**Verdict `WORSE_THAN_MARKOV_NULL` under BOTH feature modes**, by 0.125 and 0.103 — well past the derived
tolerance in each case. A 1.2-billion-parameter genome language model, prompted with 300 bp of real
upstream context, produces promoter-like sequence that is **more** detectable than an order-3 Markov chain's.

Both controls held: the self-control sat near chance (0.5447 / 0.5865) so the instrument was unbiased, and
the uniform control was high (0.9408 / 0.9641) so it had power.

## What this is NOT — two mechanisms proposed, tested, and REFUTED

Neither of my explanations survived, so the mechanism is recorded as **undiagnosed** rather than invented.

**1. Mode collapse — REFUTED.** Two samples shared conspicuous sequence and I suspected a collapsed
distribution. Measured pairwise identity over 6 independent samples: **mean 0.249, min 0.167, max 0.333**,
against a ~0.25 random expectation for independent DNA; 6 of 6 distinct. The real cause of the apparent
similarity is that `manual_seed` is set once per run, so the FIRST sample of two same-seed runs shares an
RNG state. The generator is diverse.

**2. Coding-like output — REFUTED.** The first unprompted sample began `ATGCTGCAA…`, and since a prokaryotic
genome is mostly coding, "the model emits ORFs when asked for a promoter" was the obvious story. It is wrong:

| | stop codons /100 bp (3 frames) | longest stop-free run | ATG /100 bp |
|---|---|---|---|
| natural promoters | 4.67 | 117.4 bp | 1.75 |
| GENERator output | 4.83 | 133.9 bp | 1.75 |

Coding sequence would show *markedly lower* stop density and *much longer* stop-free runs. These are nearly
identical, ATG rate matches to two decimals, and **0 of 8 sequences start with ATG in either set** — the
earlier observation was a coincidence of the unprompted samples, not a property of conditional output.

**So the model matches every coarse statistic measured — GC (0.4494 vs 0.4552), stop density, ATG rate, and
roughly stop-free run length — and is still separable at 0.72 on 4-mer frequencies.** The detectable
difference lives in higher-order k-mer structure, and this run does not identify it.

## The untested lever, named rather than claimed

**Sampling was left at defaults: `temperature=1.0`, no top-k, no top-p.** For an autoregressive model that is
a real and well-known determinant of output statistics, and nothing here tests it. So the supportable claim
is narrow:

> **GENERator-1.2B, AS SAMPLED WITH DEFAULT SETTINGS and prompted with 300 bp of context, is more
> distinguishable from real E. coli promoters than an order-3 Markov chain is.**

It is **not** "the model is bad" and **not** "genome LMs cannot generate promoters". A temperature/top-p sweep
is cheap (each point ~35 min on CPU, far less on a T4) and is the obvious next measurement.

## Honest limits

- **n=120 generated vs 600 natural.** The Markov null was generated at the same n, so the head-to-head is
  matched — but 120 is small and no confidence interval is computed on the AUROC difference.
- **One model, one organism, one window length (150 bp), one context length (300 bp), one seed.**
- **Default sampling** (above) — the single largest untested factor.
- A k-mer/positional-bin discriminator sees coarse structure only; the published audits' failure axis is
  longer-range still. This bar is a floor, not a ceiling.
- The generated sequences were **not saved** in the artifact, so the diagnostic above required a fresh
  8-sequence sample. Saving them is a cheap fix and would make any future diagnosis free.
- CPU fp32 (the local 4 GiB card cannot hold 1.2B at fp32 — the VRAM guard now catches that before loading).

## Why this is progress, not a setback

The project's standing instruction is explicit: build the GLM, and a poor number is not a reason to stop.
What this run bought is a **working, controlled, pre-registered measurement apparatus** and the first real
reading from it. Before today, "does the generator produce DNA-like sequence?" had no answer that a 3-mer
Markov chain could not also claim. Now there is a number, a null that beats it, two refuted explanations, and
a named next lever.
