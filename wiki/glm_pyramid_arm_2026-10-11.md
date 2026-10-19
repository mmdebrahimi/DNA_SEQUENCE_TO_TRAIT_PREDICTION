# The multi-scale pyramid arm: measured, no gain — and the declared 4-decade arm is substrate-gated, not expensive

**2026-10-11.** `--advance` picked the GLM umbrella's `ceiling:` criterion because R5/R6 had it
`blocked:untried-arm` and R15 reported that *every ambitious candidate in the portfolio lives in one
family*. Two arms had been declared on 2026-10-09 and never run. **This run closes NEITHER of them.** It
measures a THIRD arm that is neither declared arm (`encoder-pyramid`, 1.70 decades), and reclassifies one
of the declared two (the 4-decade arm) from *unrun* to *substrate-gated*. Artifact:
`wiki/glm_pyramid_arm_2026-10-11.json`. Reproduce:
`uv run python scripts/glm_encoder_gate.py --seeds 10 --pyramid --no-download`.

## The measurement

Substrate GSE144621 peak tiles (n=44,106, len 150, leave-peak-out), replicate-derived ceiling **0.941**,
10 registered seeds, CPU as the determinism reference.

| arm | span | median ρ | params | paired mean vs GC | wins | CI95 |
|---|---|---|---|---|---|---|
| `gc` | — | 0.3112 | — | — | — | — |
| `encoder-w20` | 1 scale | 0.3163 | 2,657 | +0.0010 | 5/10 | [−0.0202, 0.0222] |
| `encoder-multi` | 1.10 dec | **0.3234** | 13,217 | **+0.0125** | 6/10 | [−0.0099, 0.0349] |
| **`encoder-pyramid`** | **1.70 dec** | **0.3167** | **27,649** | **+0.0069** | 6/10 | [−0.0124, 0.0262] |
| `encoder-shuffled` | null | −0.0240 | 13,217 | −0.3288 | 0/10 | [−0.3597, −0.2980] |
| `encoder-pyramid-shuffled` | null | −0.0198 | 27,649 | −0.3368 | 0/10 | [−0.3651, −0.3086] |

**Widening the receptive field from 1.10 to 1.70 decades and 2.1× the parameters did not help — it came
out slightly worse** (+0.0125 → +0.0069 paired), barely separating from the single-scale arm. The
pyramid's own shuffled null is clean (0/10 wins, paired −0.3368, sign-p 0.002) and
`registered_protocol_used` / `null_clean` / `reconcile_ok` are all True, so this is a VALID negative
rather than a broken run.

**The honest reading is "no arm separates from GC", not "multi beats pyramid."** All three encoder CIs
include zero; the 0.3234-vs-0.3167 ordering is within noise. What the pyramid *does* settle is narrower
than a lever-level claim: **within 150 bp tiles, widening the receptive field from 1.10 to 1.70 decades
under the registered small-CNN protocol did not help.** That is one encoder family at one capacity on one
substrate — it is NOT evidence that scale is unhelpful in general, and this memo's own `honest_limits`
block says so. It is corroboration, not an independent finding, because this
project already recorded that capacity HURTS on this substrate (6-mers at 4,096 features scored worst of
every arm at 0.2022).

## The finding that matters more: the 4-decade arm cannot run here at any price

**Measured, not estimated.** A census of the substrate through its own loader:

| category | n | min | max |
|---|---|---|---|
| `tile` | 45,201 | 65 | **150** |
| `neg_control` | 470 | 150 | 150 |
| `random` | 936 | 150 | 150 |
| `pos_control` | 106 | 48 | 75 |

Four decades of scale needs ~**10,000 bp** of context. Every category here is ≤150 bp, so the absolute
span ceiling is log10(150/1) = **2.18 decades** — and that only with a degenerate 1 bp finest branch
that cannot see a motif at all. So `multiresolution-pyramid-4-decades-UNRUN` is **substrate-gated, not
expensive**, and its prerequisite is the *other* planned arm,
`substrate-change-constructed-eukaryotic-or-human-panel-UNRUN`.

**The two planned arms are ORDERED, not parallel** — a dependency the arms ledger does not encode, and
the reason the ceiling criterion cannot be closed by compute alone.

**The arm is left PLANNED, not killed.** Per R6 an arm dies only by its own executed evidence, and what
was measured here is the SUBSTRATE, not the arm. Recording "infeasible on the current substrate" is a
precondition; retiring the arm would be the bookkeeping-side-effect kill R6 exists to prevent. If the
substrate arm ever lands, the 4-decade pyramid becomes runnable and should run.

**Why the span tested is 1.70 and not 2.18.** Widths `(1,10,100,150)` would score the higher number, but
a width-1 kernel is a per-base recode with no motif visibility — buying span by shrinking the finest
branch inflates the decade metric without testing what the metric stands for. `(3,12,48,150)` keeps every
branch a scale someone would choose on its own: sub-motif, motif, element, whole tile.

## What this does and does not change

- A **new arm is measured** (`pyramid-4scale-1.70decades-150bp`, metric 0.3167) and recorded in
  `soraya_runs/arms-ledger.json` with its evidence. `best` correctly remains `conv-3width-13k-params-150bp`
  at 0.3234 — the pyramid did not win, and the code computed that rather than being told it.
- The ceiling criterion stays **`blocked:untried-arm`**. That is correct: the arm measured here is not
  either declared arm. Nothing was retired to make a verdict look settled.
- The gap to the ceiling is unchanged and large: best measured 0.3234 against 0.941.
- `scripts/glm_encoder_gate.py` gains `--pyramid` as an **additive** flag; without it the three published
  arms and the reconcile gate are byte-identical, so the 2026-10-08 result is untouched.

## Correction 2026-10-12

Two statements in the 2026-10-11 version of this memo were wrong and are corrected above. Recording them
rather than editing silently, because what a memo CLAIMED is part of its evidence trail:

- The opening said this run *"closes one of them with a number and reclassifies the other."* It closed
  NEITHER declared arm. It measured a THIRD arm that is neither of them. The sentence contradicted this
  memo's own `## What this does and does not change`, which stated correctly that "the arm measured here
  is not either declared arm" and that the criterion stays `blocked:untried-arm`. Both cannot hold.
- The reading said *"more scale is not the missing lever"* — a lever-level negative drawn from one encoder
  family at one capacity on one substrate, which over-reached this memo's own `honest_limits` block. The
  supportable claim is scoped to 150 bp tiles and the registered small-CNN protocol.

The MEASUREMENT was never in question: the results table, the paired statistics, the shuffled nulls and
the substrate census are unchanged. Only the two over-claims above were edited.
