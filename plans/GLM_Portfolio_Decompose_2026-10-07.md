# GLM portfolio — decomposition, flow-down, critical path

**Date:** 2026-10-07 · **Mission:** the GLM (name a trait, get the genome edits) · **Namespace: `G-A`…`G-E`**

## Namespace note — read this before adding a family

`F-A`…`F-E` are **taken** by the honest-decoder portfolio in `project_state/dna-decode-2026-05-11.md`
(doubt layer / catalog curation / evidence surface / learned-narrow / label acquisition), and `F1`…`F4` are
**taken** by the in-flight deep-research campaign in `NEXT.md`. Reusing either would collide with live work.
The GLM families are therefore `G-*`, and they live under their **own umbrella** — not that one — because
the existing umbrella's Mission Terminal Condition explicitly **excludes** a learned predictor
("*Not terminal: a learned predictor that beats the curated catalog in the natural-population regime*").
One umbrella carries one terminal condition; two missions need two umbrellas.

## Mission terminal condition (GLM)

A user names a target molecular phenotype and receives a **ranked set of concrete genome edits**, each
scored by an oracle that has been **validated on the same sequence distribution the edits are drawn from**,
with an honest statement of what the score does not support.

**Not terminal:** an oracle that scores well on a designed grid while the generator emits genome-like
sequence, or vice versa. That configuration can have both halves "working" and still be unsound — which is
precisely the gap measured on 2026-10-07 and the reason `G-E` exists.

## Where the mission actually stands (measured, 2026-10-06 → 10-07)

| component | state |
|---|---|
| edit primitive | **SHIPPED** — `GenomeEdit`/`EditSet`, one span-replacement primitive, 46 tests |
| discriminator (naturalness) | **SHIPPED** — self-validating controls; refuses when its own controls fail |
| oracle (scorer) | **WORKS ON A DESIGNED GRID ONLY** — ρ 0.59 leave-element-out; captures 79% of achievable selection value at k=10. On real genomic peak tiles a single GC number (0.300) **beats** every learned feature (one-hot 0.269), at 29% of the derived 0.941 ceiling |
| generator | **LOSES TO A 3-MER MARKOV CHAIN** — distinguishability 0.7706 vs the null's 0.6679 (lower is better), at default sampling |
| condition signal | **REAL** — noise-matched cross-condition 0.4062 below within-condition 0.4676; disattenuated true cross-condition 0.789, so ~38% of variance is condition-specific |
| self-distillation | **DEGENERATE on the current oracle** — linear student reproduces its teacher at agreement **1.0000** |

**The binding conclusion:** the scorer is not the bottleneck on grid-like sequence and *is* on genome-like
sequence; the generator is below a trivial baseline; and the two halves are validated on different
distributions. That is a three-way problem, which is why this decomposes into five families and not one.

## The five families

| # | Family | Slug | blocked_by | Hardware | Scope |
|---|---|---|---|---|---|
| **G-A** | Learned representation for the oracle | `glm-learned-representation-2026-10-07` | none | CPU ok, GPU faster | replace fixed features with a learned multi-resolution encoder; trivial head retained |
| **G-B** | Self-distillation | `glm-self-distillation-2026-10-07` | **G-A (hard, measured)** | CPU | manufacture a variant-effect signal by perturb-and-relabel; validate on real data or not at all |
| **G-C** | Condition-conditioned oracle | `glm-condition-conditioning-2026-10-07` | none | CPU | one model conditioned on growth medium vs two separate models |
| **G-D** | Generator fine-tuning | `glm-generator-finetune-2026-10-07` | sampling sweep reports first | **8 GB GPU** | the un-pulled lever on a measured ZERO-SHOT failure |
| **G-E** | Distribution alignment | `glm-distribution-alignment-2026-10-07` | **G-A + G-D** | CPU | characterise what the generator emits; validate an oracle on *that* distribution |

**Critical path: G-A → G-E.** `G-D` joins it from the side; `G-B` hangs off `G-A`; `G-C` runs independently.

## Requirements flow-down

- **G-A → G-B is HARD and measured, not a preference.** A ridge teacher's predictions lie in the span of
  its own features, so a linear student recovers the teacher exactly — agreement **1.0000** at mutation
  rates 0.02/0.05/0.10, with a GC-only control collapsing 0.596 → 0.174 to prove the probe is non-vacuous.
  There is nothing for the perturbation to regularise until the student has capacity the teacher lacks.
  **Do not plan G-B work as if it can start early.**
- **G-A → G-E.** `G-E` must score generated sequence with an oracle validated on the matching distribution.
  Today no such oracle exists for genome-like sequence (GC beats the learned features there), so `G-E`
  cannot produce a sound number until `G-A` either succeeds or closes.
- **G-D → G-E.** Characterising generator output is only worth doing on a generator that is not below a
  trivial baseline; otherwise `G-E` measures the distribution of a known-broken sampler.
- **G-D is gated on the sampling sweep reporting first.** Sampling and fine-tuning are two generator levers;
  running them together confounds which one moved the number. The sweep is already queued on the server.
- **G-C is independent and must stay independent.** It is the only family whose motivating measurement is
  already positive, so it is the portfolio's hedge: if `G-A` closes negative, `G-C` still has a live result.
- **Nothing downstream may be planned as if `G-A` will succeed.** Fixed features failing to transfer is a
  measured fact; a learned encoder transferring is a hypothesis with a stated falsifier.

## Parallel vs sequential — the judgement the user asked for

**Start two in parallel, now: G-A and G-C.** Both are CPU-feasible, both have their substrate already
downloaded and loaders already written, and they are logically independent — G-A changes the *representation*,
G-C changes the *conditioning*, and they touch different axes of the same model. Running them together also
hedges the portfolio: G-A is the higher-ceiling bet with a real chance of closing negative, while G-C's
motivating effect is already measured, so a negative G-A does not leave the quarter empty.

**G-D waits on hardware, not on logic.** The 1.2B model needs the 8 GB card and that machine is in a
basement mid-move. It also must wait for the sampling sweep so the two levers stay separable.

**G-B waits on G-A by construction** — the degeneracy is measured, not suspected. It is also the smallest
family: once a learned encoder exists, the distillation step is a training-loop variant, not a new system.

**G-E is last and is the integration.** It is the family that decides whether the loop is sound, and it
cannot be faked earlier: its inputs are an oracle validated on genomic sequence and a generator worth
characterising.

So: **2 running · 1 hardware-blocked · 2 downstream.** Not five in parallel — three of the five have real
predecessors, and pretending otherwise would produce five half-finished families instead of two finished ones.

## What each family must NOT do (anti-scope, carried into each ledger)

- **G-A** must not report a random-split number as its headline. Random splits inflate ~3× in this repo,
  measured four times; and on peak tiles the honest split is leave-peak-out because tiles within a peak
  overlap heavily.
- **G-B** must not report distillation *agreement* as success. Agreement with the teacher is what a
  degenerate student also achieves. The only acceptable headline is held-out performance on real measured
  data.
- **G-C** must not compare a 2-replicate mean against a single replicate. That unmatched comparison reads
  0.5025 and inverts the condition answer; matched, it is 0.4062 vs 0.4676.
- **G-D** must not quote a distinguishability number without its corpus fingerprint (`n_cds_total`,
  `n_windows`, `n_skipped_edge`, `n_overlapping_cds`, `strand_counts`). NCBI re-annotates assemblies under
  the same accession, so two runs can look comparable and not be.
- **G-E** must not normalise one arm by another arm's noise ceiling. Measured: the fragment ceiling (0.793)
  and the tile ceiling (0.941) differ enough that a shared normaliser silently rescaled an arm by 5 points.

## Portfolio-level falsifier

If **G-A closes negative** (a learned encoder cannot beat a single GC number on real genomic promoter
tiles) **and G-D closes negative** (fine-tuning does not get the generator below a 3-mer Markov chain), then
the sequence-model lever space for this substrate is exhausted, and the live question becomes the **label and
substrate** question rather than the architecture question — which routes back to the existing `F-E` label
acquisition family rather than to a sixth GLM family. Naming that in advance is what stops a sixth family
from being invented to avoid the conclusion.
