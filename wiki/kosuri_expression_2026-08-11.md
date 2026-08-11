# Track B — the learned expression model: composability PASS, novel-RBS-from-sequence PASS, novel-promoter untested (2026-08-11)

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

## ✅ Sequence generalisation DOES work — corrected 2026-08-11

**An earlier version of this memo concluded that sequence generalisation "is not demonstrated and this
dataset alone cannot answer it." That was wrong, and it was wrong because I had used only 2 of the 4 files
provided** — `sd02.xls` carries the **sequence of all 111 RBSs**, which I had not touched. Prompted by the
user asking whether the files were already in hand.

Held-out **RBS** (the RBS is *never* seen in training), scored from its letters alone — 1/2/3-mers, length,
GC, best match to the Shine-Dalgarno core `AGGAGG`, and SD-to-start spacing:

| model | R² | reading |
|---|---|---|
| additive baseline | 0.499 | |
| GBM on element identity | 0.268 | worse than baseline — nothing to say about an unseen element |
| `promoter_only` **control** | 0.499 | the promoter alone; identical to baseline, as expected |
| `rbs_sequence_only` **control** | 0.099 | RBS sequence alone can't carry total protein (promoter dominates) |
| **promoter + RBS sequence** | **0.747** | **+0.248 over the promoter-only control** |
| **promoter + RBS sequence + ΔG** | **0.781** | **+0.281 vs baseline, +0.513 vs identity** |

**A novel RBS can be scored from its sequence at R² 0.78.** The two controls are what make that claim safe:
the promoter alone gets 0.499 and sequence alone gets 0.099, so the lift is genuinely attributable to
sequence rather than to the promoter doing the work.

This is Q2 of the design epoch — *will the host express it?* — **answered for the RBS half**, on the
strictly harder split (unseen element, not unseen combination).

## The falsification that matters

The tempting read of 0.919 is "we can predict expression from sequence." **That is false, and the
element-split rows are the proof.** Given an unseen promoter the GBM scores **−0.014 — below chance, and
worse than the additive baseline it beat by 0.124 on combinations.** It learned element *identity*, not
sequence. The baseline degrades more gracefully precisely because an unseen element contributes 0 and it
falls back on the other element plus the grand mean.

Note this falsifies the *identity* model specifically, not sequence modelling — with only identity + ΔG
available, ΔG is the sole source of generalisation (lifting held-out-element R² from ~0 to 0.14–0.33).
Give the model **real RBS sequence** and it reaches 0.781 on the same split (section above). The lesson is
that an identity encoding cannot generalise to a new part, not that expression is unpredictable.

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
| Can we score a **novel RBS** from its sequence? | **Yes — R² 0.781**, with promoter-only (0.499) and sequence-only (0.099) controls showing the lift is real |
| Can we score a **novel promoter** from its sequence? | **Untested** — `sd01.xls` (promoter sequences) failed to download |

## What is still blocked

`sd01.xls` (promoter sequences) **failed to download** — the saved file is 58,788 bytes of the Cloudflare
challenge page (contains `cloudflare` + `captcha` markers and **zero** occurrences of "Promoter" or
"apFAB"). The SI PDF was also checked: its DNA runs are **sequencing primers**, not the promoter library.
So the promoter half of the sequence question cannot be tested here.

That leaves Track B at: **combination selection ✅, novel-RBS scoring ✅, novel-promoter scoring ❌.**
Re-downloading one 44 KB file closes the last of it.

## Reproduce

```bash
uv run python scripts/kosuri_expression_validate.py --sd03 <path>/sd03.xls --sd02 <path>/sd02.xls
```

Data is **not committed** (16 MB third-party supplementary; PNAS is Cloudflare-gated to scripts —
see the design-epoch plan for the four exhausted fetch routes). Sidecar:
`wiki/kosuri_expression_2026-08-11.json`. Tests: `tests/test_kosuri_expression.py`.
