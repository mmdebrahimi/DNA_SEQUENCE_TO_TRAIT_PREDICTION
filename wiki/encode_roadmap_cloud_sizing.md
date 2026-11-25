# ENCODE / Roadmap regulatory goal — cloud sizing memo

**Date:** 2026-11-24
**Discharges:** constraint `encode-roadmap-regulatory-goal-unsized-on-cloud` (`soraya_runs/constraint-expiry.json`)
**Inputs:** `wiki/encode_sizing_t4_2026-11-24.json` (measured) · `wiki/encode_roadmap_training_budget_sources.json` (sourced) · `wiki/encode_substrate_census_2026-11-24.json` (measured)

---

## DECISION NEEDED

**The compute cut from 2026-05-12 no longer binds, and a different constraint has taken its place.
Whether to act on that is a user authority call, not this memo's.**

The goal abandoned on 2026-05-12 was a tissue-specific gene-regulatory-element activity predictor
built on a transformer-based DNA model trained on ENCODE / Roadmap multi-tissue data. It was given
up because the local GPU was a 4 GiB GTX 860M. Measured on a free Kaggle T4, **the reference
architecture's full 196,608 bp context fits, at 1,024 channels and 134.63M parameters, in 10.86 of
14.56 GiB — 3.70 GiB spare.** Compute is no longer the reason this cannot be attempted.

What binds now is **data volume: an estimated 64,046.7 GiB** across the four assay families such a
model needs, against 21.3 GB free on this host's C: drive (measured 2026-11-24; an earlier draft said "about 26 GB", carried over from the plan rather than measured). That is the finding Step 5 registers as a
new constraint with its own expiry test, rather than leaving it in prose.

Three decisions follow, and all three are yours:

1. **Is a reduced-context or reduced-width variant still "the goal"?** The context matches exactly;
   the width does not (1,024 measured against 1,536 sourced). Whether that counts as the original
   goal or a different one is a scope question.
2. **Is a subset substrate acceptable?** Histone ChIP-seq alone carries 58,757.2 of the 64,046.7
   GiB. DNase-seq plus ATAC-seq together is **4,672.1 GiB** — still far past this laptop, but a
   different proposition. (Stated in GiB because that is the census's own unit: the same quantity is
   4.563 TiB or 5.02 TB depending on convention, and an earlier draft of this memo said "about
   4.7 TB", which matches neither. Unit conversion is where this figure goes wrong.)
3. **Should anything be built at all?** This memo prices the goal. It does not recommend pursuing it.

---

## CORRECTIONS

**Two claims in the plan that commissioned this memo are not supported, and are corrected here
rather than repeated.**

**1. The reference parameter count cannot be compared.** The plan asked this memo to state
"134.63M params against the reference's larger count". Step 2 found **no published total parameter
count** for the main reference model, so that comparison has no sourced right-hand side. It is
replaced by one that is sourced on both sides: **1,024 channels measured against 1,536 channels
sourced.** The conclusion the plan wanted — that this reaches the context and not the scale — still
holds, on the channel-width evidence instead of a parameter count.

**2. The AMP speedup is not one number.** The plan cited roughly fourteen-fold from the 384-channel
pair. Measured, the ratio **varies strongly with shape**: 14.0× at 384 channels (2683.7 → 191.7 ms)
but **4.2× at 768** (2960.2 → 702.0 ms). The largest shape has no fp32 timing to compare against at
all, because that configuration is the one that does not fit. Quote the range, not the 14×.

---

## The measured envelope

Tesla T4, sm_75, 14.56 GiB total VRAM, torch 2.11.0+cu130, 128 bp output bins, 5,313 output tracks,
a 9.0-hour session assumed. 45 configurations swept; **44 fit.**

| configuration | mode | fits | ms/step | peak VRAM | steps per 9 h |
|---|---|---|---|---|---|
| 196,608 bp · 1024 ch · 11 L · 134.63M | fp32 | **No — OOM** | — | — | — |
| 196,608 bp · 1024 ch · 11 L · 134.63M | AMP | Yes | 1028.3 | 10.86 GiB | 31,509 |
| 196,608 bp · 1024 ch · 11 L · 134.63M | AMP + ckpt | Yes | 1067.5 | 9.53 GiB | 30,351 |
| 196,608 bp · 768 ch · 8 L | fp32 | Yes | 2960.2 | 7.64 GiB | 10,945 |
| 196,608 bp · 768 ch · 8 L | AMP | Yes | 702.0 | 5.18 GiB | 46,151 |
| 196,608 bp · 384 ch · 4 L · 11.97M | fp32 | Yes | 2683.7 | 3.51 GiB | 12,073 |
| 196,608 bp · 384 ch · 4 L · 11.97M | AMP | Yes | 191.7 | 2.4 GiB | 169,033 |

**The single non-fitting configuration is fp32 at the largest shape.** That is the load-bearing
detail: half precision is the **enabling lever here, not an optimisation**. Without it the largest
shape does not run at all; with it the same shape leaves 3.70 GiB unused. A non-fitting row carries
no timing fields whatsoever — only a `reason` — so an OOM can never be mistaken for a slow result.

### The two levers, recorded separately

**AMP — enabling, and a throughput lever of shape-dependent size.** It converts the largest shape
from OOM to fitting. Where both modes ran, it is 14.0× faster at 384 channels and 4.2× at 768,
consistent with the sm_75 tensor cores doing the work. One number would misrepresent it.

**Gradient checkpointing — near-neutral on this card.** At the largest shape it trades 1.33 GiB
(10.86 → 9.53) for 3.8% less throughput (1028.3 → 1067.5 ms). At 384 channels it *increased* peak
memory slightly (2.4 → 2.43 GiB). Once AMP was on, **memory stopped binding**, so the lever that
exists to buy memory has nothing left to buy.

---

## The measured capacity beside the sourced budget

| quantity | reference (sourced) | here (measured) | reading |
|---|---|---|---|
| input context | 196,608 bp | 196,608 bp | **matches exactly** |
| model width | 1,536 channels | 1,024 channels | reaches the context, not the width |
| training hardware | 64 TPU v3 cores | 1 T4 | different hardware class |
| examples per optimiser step | 64 | 1 | **64× fewer** |
| total parameters | *unsourced* | 134.63M | **not comparable** |
| main-run step count | *unsourced* | 31,509 per 9 h session | **not comparable** |
| wall-clock training time | *unsourced* | 9 h per session | **not comparable** |

**The ratio the plan asked for does not exist.** "Reference training compute divided by one measured
T4 session" needs the reference run's duration or its step count, and **neither is published** — a
scan for any duration-shaped expression across all 78,591 readable characters of the paper returned
zero matches. Step 2 records it `unsourced` with no number rather than approximating it, because a
fabricated budget produces a fully self-consistent feasibility verdict that is entirely wrong.

So "sufficient" is **a reading, and a partial one**: the architecture's *context* is reachable on
free hardware, and *how much of its training* one could complete is unknown and not knowable from
the published record. The honest statement is the first half only.

**Do not compare the step counts.** 31,509 measured steps carry one example each; the reference's
carried 64. A step here is not a step there.

---

## What this does NOT conclude

1. **No accuracy claim exists.** Every configuration was driven with random tensors
   (`torch.randn`). Nothing was learned, so nothing here speaks to whether such a model would work.
2. **Nothing was trained to convergence.** These are per-step timings and peak-memory figures
   extrapolated to a session, not a training run.
3. **This reaches the CONTEXT, not the scale.** 1,024 channels against the sourced 1,536. The
   parameter comparison is unavailable, not favourable.
4. **One free T4, one session, no scaling.** No paid tier, no multi-GPU, no claim about what a
   larger budget would buy.

A fifth, from the census: **the substrate was counted, not acquired.** Nothing has been downloaded.

---

## Why this memo's existence is what expires the constraint

The constraint `encode-roadmap-regulatory-goal-unsized-on-cloud` carries
`expiry_kind=file-exists` with `expiry_arg` pointing at this file. **The test fires on this file
existing — not on its verdict.** A pessimistic memo would expire the constraint exactly as a
favourable one does, and that is deliberate: the constraint records that the goal was *never sized
on cloud hardware*, so the thing that discharges it is **the question being answered**, whichever
way it comes out.

This is stated plainly because it is the half of `file-exists` semantics a later reader is most
likely to misread. The register will show this constraint expired; that means *it was measured*, and
nothing more. The verdict is in this document, and it is mixed: **compute cleared, data did not.**
