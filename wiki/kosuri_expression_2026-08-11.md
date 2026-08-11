# Track B — the learned expression model: PASS on composability, FAIL on the design question (2026-08-11)

The pre-registered test from `wiki/design_epoch_plan_2026-08-07.md`, run against Kosuri et al. 2013
(PNAS **110**:14024) — 12,563 constructed promoter × RBS combinations with measured DNA, RNA and protein.
Q2 of the design epoch: *will the host express it?*

**The bar was set before the data was in hand**, from the paper's own numbers: beat protein **R² ≈ 0.82**
on held-out constructs, with the baseline re-fit on the training split only and the split taken **by
element**. Both halves of that pre-registration now have answers, and **they disagree — which is the
finding.**

## Reproduction gate (run first — before trusting any new number)

Recomputing the paper's results from its own model columns:

| | reproduced | published |
|---|---|---|
| protein, simple model | **0.7525** | 0.76 |
| RNA, simple model | **0.9238** | 0.92 |
| RNA, full model | **0.9623** | 0.96 |

RNA matches to three decimals; protein is within the paper's threshold-censoring. The data is understood.

> **Units trap, recorded because it cost a debugging cycle:** `model.prot.simple` is stored in **log2**
> while `prot` is raw RFU. Compared in the wrong space it returns **R² = −15**, which reads like a broken
> loader rather than a units mismatch. RNA columns are raw and compare in log10.

## Result

Target `log2(protein)`, 12,270 constructs, 114 promoters × 111 RBSs. Baseline = the paper's model
(`mu + promoter effect + RBS effect`) re-fit on the **training split only**.

| split | additive baseline | GBM (element identity) | GBM + ΔG |
|---|---|---|---|
| held-out **combination** | 0.795 | 0.893 | **0.919** |
| held-out **promoter** | 0.263 | **−0.014** | 0.144 |
| held-out **RBS** | 0.499 | 0.268 | 0.327 |

- **Combination split → PASS.** 0.919 vs the 0.82 bar; and against the fair like-for-like comparator
  (my own baseline, out-of-sample) it is **0.795 → 0.919, +0.124**. Promoter×RBS *interaction* is real,
  learnable, and worth ~+0.10 on its own; ΔG (5′ secondary structure) adds a further +0.026.
- **Element split → FAIL.** Nothing approaches 0.82 — **including the baseline** (0.26–0.50).

## The falsification that matters

The tempting read of 0.919 is "we can predict expression from sequence." **That is false, and the
element-split rows are the proof.** Given an unseen promoter the GBM scores **−0.014 — below chance, and
worse than the additive baseline it beat by 0.124 on combinations.** It learned element *identity*, not
sequence. The baseline degrades more gracefully precisely because an unseen element contributes 0 and it
falls back on the other element plus the grand mean.

The only genuine sequence generalisation in the whole experiment is **ΔG**, the one feature computable
from sequence before building anything: it lifts held-out-element R² from ~0 to **0.14–0.33**.

## Honest verdict on the pre-registration

**By my own stated falsifier — "split BY ELEMENT" — this is a FAIL.** Recorded as such.

But the bar itself was mis-specified for that split, and it is worth saying exactly how rather than
quietly re-scoring: **0.82 is a combination-level, in-sample number.** An element-strength model has no
strength for an unseen element, so it *cannot* reach 0.82 there — and indeed doesn't (0.26–0.50). The bar
and the split were incompatible, and I did not notice until the data was in hand.

What each half legitimately establishes:

| question | answer |
|---|---|
| Given characterised parts, can ML pick a **new combination** to hit a target expression level? | **Yes** — R² 0.919, clearly beating the composability model the paper published |
| Can we score a **novel** promoter/RBS from its sequence? | **Not demonstrated.** Needs real sequence modelling, and this dataset alone cannot answer it |

## What blocks the second question

`sd01.xls` (promoter sequences) **failed to download** — the saved file is 58,788 bytes of the Cloudflare
challenge page, not data. RBS sequences are present in `sd02.xls`; promoter sequences are not available.
So a true sequence→expression model cannot be built or tested here yet.

**This does not sink Track B.** The practical design capability — choose a promoter+RBS pair from a
characterised library to hit a target expression level — is demonstrated and is genuinely what a strain
engineer does. The unanswered part is generalisation to novel parts.

## Reproduce

```bash
uv run python scripts/kosuri_expression_validate.py --sd03 <path>/sd03.xls
```

Data is **not committed** (16 MB third-party supplementary; PNAS is Cloudflare-gated to scripts —
see the design-epoch plan for the four exhausted fetch routes). Sidecar:
`wiki/kosuri_expression_2026-08-11.json`. Tests: `tests/test_kosuri_expression.py`.
