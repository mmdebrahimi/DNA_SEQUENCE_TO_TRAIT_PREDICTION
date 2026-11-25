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

**What binds now is MEASURED, and it is neither compute nor storage.** A first draft of this memo
asserted a data wall of 64,046.7 GiB against 21.3 GB of free local disk; both sides were the wrong
quantity, and the section below retracts it. Measured 2026-11-25, the real limits are **19.52 GiB of
persistent output** and **30 GPU hours per week** — and the open question is **read BANDWIDTH**, not
capacity. See *The data question*.

Two decisions follow, and both are yours:

1. **Is a reduced-width variant still "the goal"? DRAFTED: yes — same goal, not a new one.**
   The 2026-05-11 goal named an architecture CLASS and a dataset ("a transformer-based DNA model
   trained on ENCODE / Roadmap Epigenomics multi-tissue data") and a target (tissue-specific
   regulatory-element activity). It named no parameter count and no channel width. A 1,024-channel
   model at the reference's **full** 196,608 bp context is that same class, on that same data, for
   that same target — so this is a `same-premise-instantiation`, not a superseding goal. Recorded
   for ratification rather than left as a blocking question; redirect if you read the goal as
   width-bearing.
2. **Should anything be built at all?** This memo prices what can be priced and does not recommend
   pursuing the goal. **The blocker is no longer a lookup** — storage was measured on 2026-11-25 and
   does not bind. What is now untested is **read bandwidth** over a gcsfuse mount for a ~302 GiB
   tensor per epoch, which is an empirical question a short probe could settle.

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

**3. RETRACTED — "data is the wall", as this memo first stated it, compared the wrong numerator to
the wrong denominator.** The retraction is kept visible rather than edited away, because the claim
was published in commit `615c5467` and reported aloud before it was checked.

- **Wrong numerator.** 64,046.7 GiB is the census's estimate of *every released bigWig in those four
  assay families*. That is not a training set. It is unfiltered by organism (it includes mouse),
  by assembly, by output type and by replicate, and it is raw signal rather than 128 bp bins. The
  sourced track composition of one reference example is 2,131 TF ChIP-seq, 1,860 histone, 684
  DNase/ATAC and 638 CAGE tracks — so the census counted **20.7× the histone tracks and 16.4× the
  DNase/ATAC tracks the reference actually used.**
- **Wrong denominator.** Training would run on the cloud card this memo measures, so the relevant
  limit is Kaggle's storage, not this laptop's C: drive. Comparing a cloud-trained model's data
  requirement to local free space answers a question nobody asked.
- **Computed properly, from figures already sourced**, the training-set-shaped target tensor is
  34,021 sequences × bins × 5,313 tracks: **301.7 GiB** at the reference's 896 predicted bins in
  float16, 517.1 GiB at the 1,536 bins measured here, 603.3 GiB at 896 bins in float32. The input
  sequence adds 6.2 GiB and is negligible beside it. **The retracted figure was a 212×
  overstatement** of the quantity that matters.

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

**One measured difference runs in feasibility's favour, so it is stated rather than banked.** The
reference predicts the **central 896 bins (114,688 bp)** of its 196,608 bp input — sourced. The sweep
measured here used `bins = context / 128 = 1,536`, predicting the **full** context, so it carried a
**larger output head than the reference's** across all 5,313 tracks. The VRAM and throughput figures
above are therefore conservative with respect to the reference's actual output shape, not flattering
to it.

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
nothing more. The verdict is in this document, and it is partial: **compute cleared; the data
question is RESOLVED as of 2026-11-25 — storage does not bind; read bandwidth is the open
empirical question.**

---

## The data question

**Status: RESOLVED 2026-11-25 — storage does not block this goal.** Measured, not looked up: see
`wiki/kaggle_storage_limits.md`.

**The requirement, computed from sourced figures.** 34,021 training sequences × bins × 5,313 tracks:
**301.7 GiB** at the reference's 896 predicted bins in float16, 517.1 GiB at the 1,536 bins measured
here, 603.3 GiB at 896 bins in float32. Plus 6.2 GiB of input sequence. This is a derived figure, and
its inputs are the sequence count and track count sourced in Step 2 plus a bin count and dtype chosen
here — so it is an estimate of a *tensor*, not a measurement of a download.

**The limit, MEASURED.** The docs route stayed dead — every Kaggle-owned page returns only its HTML
title to a fetch, across seven attempts and five URLs — so the question was answered by measuring a
running session instead, which is better evidence than a published policy anyway. A CPU-only kernel
reports **8.0 TiB of ephemeral scratch with 2.5 TiB free**, and an attached dataset mounts separately
through **gcsfuse reporting 1.0 PiB, read-only, consuming no local disk at all**. `kaggle quota`
supplies the accelerator half from the API: **30.00 h GPU and 20.00 h TPU per week**.

I had framed this as a documentation problem. It was a measurement problem, and the fix was to stop
trying to read the number and go take it.

**So the status is:** a 301.7–603.3 GiB tensor fits ephemeral scratch roughly **four to eight times
over**, and a dataset large enough to hold it need not fit the session filesystem at all, because it
streams from cloud storage. **Storage is not the binding constraint.**

**What binds instead**, both measured: **persistent output at 19.52 GiB** — which limits what a run
may KEEP, not what it may read (a checkpoint of this 134.63M-param model is roughly 0.5–1.6 GiB, so
that holds a dozen) — and **30 GPU hours per week**, which at the measured 31,509 steps per 9-hour
session is about **3.3 sessions per week**.

**The untested question is BANDWIDTH, not capacity.** A mount reporting 1.0 PiB says nothing about
read speed, and a 302 GiB tensor streamed per epoch over gcsfuse is plausibly the real limit. That is
now the open empirical question, and it is a different one from the question this memo was asked.

**Not addressed here either way:** nothing has been downloaded, the per-file byte figures in the
census are themselves median-sampled estimates, and the pipeline that would turn released bigWigs
into a binned training tensor does not exist in this repo.
