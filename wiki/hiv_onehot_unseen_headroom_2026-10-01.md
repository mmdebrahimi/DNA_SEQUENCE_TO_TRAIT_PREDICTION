# The modality channel is not worth building — one-hot's blindness is real and absorbed

**2026-10-01.** The gLM plan's T2 was the one remaining measurement-supported build: add an **orthogonal
modality** (structure/conservation) to the supervised blind-spot complement, rather than more parameters.
This probe was the cheap pre-check, and it **kills the build before it was paid for.**

## The hypothesis was sharp, and half of it is confirmed

One-hot encodes each `(position, amino-acid)` pair as an **independent column**, so a substitution with no
training signal in its own fold gets zero weight. Positions are islands — one-hot knows nothing about 181
being near 188, or valine resembling isoleucine. A position-level prior *can* score such a substitution;
that was the orthogonality argument.

**The blindness is real, and measured rather than assumed: the maximum |coefficient| over all
zero-training-signal columns is exactly `0.0`, on both genes.** The mechanism is exactly as hypothesised.

## But the room isn't there

| gene | blind spot | carry ≥1 unseen sub | exposed AUROC | unexposed AUROC | gap |
|---|---|---|---|---|---|
| RT (EFV) | n=1,883, 201 R | **7.5%** (R 12.4% / S 7.0%) | **0.8229** | **0.8163** | **−0.0066** |
| IN (RAL) | n=844, 59 R | 15.0% | *withheld — 8 R* | 0.9168 | — |

**Verdict by the frozen rule: `NO_HEADROOM_ONEHOT_ALREADY_COVERS_IT`.** Exposure is 7.5% against a 20%
bar, and — decisively — **the model scores slightly *better* where it is structurally blind**. The gap is
in the wrong direction. Either condition alone would have been enough.

The rule was frozen before the run and `self_check` drives all three verdicts, so it could have returned
`HEADROOM_EXISTS`. The floor model reproduces the committed 0.8102 / 0.8923 **exactly**, so the strata are
measured against the right model.

IN's exposed stratum holds **8 resistant isolates** — below the 10-per-class floor — so its AUROC is
**withheld** and the cell reads `UNDERPOWERED` rather than contributing a noisy number.

## Why — measured, not asserted

Per-isolate **redundancy** absorbs it. These are clinical isolates carrying many substitutions each:

- RT: **12.44** active columns per isolate on average; among exposed isolates the **median blind fraction
  is 0.0833 — one column in twelve.** The eleven seen substitutions outvote the blind one.
- IN: 7.32 columns, median blind fraction 0.1111.

**The honest exception, kept visible:** RT's **maximum** blind fraction is **1.0** — exactly **one** isolate
has *every* substitution unseen. One isolate cannot move an AUROC, but it is the single case where the
blindness is total, and it is reported rather than smoothed away.

This causal claim was measured because this repo has had **three asserted causes refuted by measurement**
in one arc; the descriptive statistic is the test.

## It reconciles the two modality results

The forward cell's orthogonal-modality lift (ESM2+GEMME+ProSST, +0.056, 90.5%) was measured on **single
variants** — exactly the case where a model must score one substitution with no other signal, so a prior is
all there is. Our blind spot is **multi-substitution isolates**, where the model already holds redundant
signal. The same mechanism predicts the blindness *would* bite on single-mutant prediction, which is
**consistent with** the 2026-07-03 leave-one-position-out head reaching only 0.691 and losing to the
catalog — a different slice and split, so consistent, not proven.

> **The reusable rule: a modality adds where the incumbent representation is starved of signal, not where
> it is merely blind in principle.** Count the incumbent's surviving signal per case before buying a
> channel.

## What this licenses

**Do not build the T2 channel.** With this, every *technical* lever on the gLM is measured and none pays:

| lever | verdict |
|---|---|
| zero-shot likelihood | fails (0.4485; BLOSUM62 matches it) |
| scale | regresses (650M → 3B → 15B) |
| capacity over the genotype | loses (4/4, 2 significant) |
| pretrained embeddings + head | PARTIAL — loses to catalog |
| **orthogonal modality** | **no room — this probe** |

The remaining lever is **labels / arm selection**, which is a user decision and explicitly not an executor
task. The gLM question is now closed on the technical side with evidence.

## Honest scope

- IN-DISTRIBUTION to Stanford, ~96% subtype B; two genes of one virus.
- A zero-training-signal column still **exists** (featurization is global); what is absent is *signal*,
  which is why the coefficient magnitude is reported directly.
- Exposure is per **fold**, against that fold's own training rows — the only honest denominator.
- This bounds the **room** for a position-level prior. It does not show such a prior would capture room if
  room existed.

## Reproduce

```
uv run python scripts/hiv_onehot_unseen_position_headroom.py --self-check   # rule reachability, no data
uv run python scripts/hiv_onehot_unseen_position_headroom.py              # 2 genes
```

Artifact `wiki/hiv_onehot_unseen_headroom_2026-10-01.json`; 10 tests
`tests/test_hiv_onehot_unseen_headroom.py`. Exits 3 if the floor does not reproduce. Frozen AMR surface
byte-unchanged (read-only).
