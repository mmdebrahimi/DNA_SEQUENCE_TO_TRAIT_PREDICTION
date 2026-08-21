# FVA requirement classes — and why this run only *partly* answers the question

**Date:** 2026-08-21 · **Artifact:** `wiki/fba_fva_requirement_2026-08-21.json`
**Script:** `scripts/fba_fva_requirement_class.py` · Model-only (no `feba.db`; D: still disconnected)

The last run ended with a redirect: *"class A is dominated by unblocked-but-idle genes, so the deficit is
an OBJECTIVE problem."* This run tested that with FVA at 100 % of optimum, which is basis-independent
where the earlier single-solution read was not.

## Result — 131 committed essential genes, evaluated in glucose

| class | n | share | meaning |
|---|---|---|---|
| REQUIRED | 55 | 42.0 % | 0 ∉ [min,max] — every optimal solution routes flux here |
| CAPABLE_BUT_IDLE | 18 | **13.7 %** | can carry flux, doesn't have to → objective doesn't demand it |
| INACTIVE_IN_CONDITION | 58 | 44.3 % | pinned at zero *in this medium* |

**The redirect is WEAKENED, not confirmed.** CAPABLE_BUT_IDLE — the class that would make this an
objective-composition problem — is only **13.7 %**, not the dominant story I expected.

## The confound, stated plainly

**I evaluated all 131 genes in glucose, but they are conditionally essential across an 11-condition
panel.** A gene essential on xylose is *correctly* inactive on glucose. So `INACTIVE_IN_CONDITION` at
44.3 % is substantially an artefact of the evaluation condition, not a finding about the model.

The clean version evaluates each gene **in the condition where it is experimentally essential** — which
needs the per-gene × per-condition experimental labels, i.e. `feba.db` on the disconnected `D:`. The
committed ratios artifact carries the model's ratios but not the experimental labels.

So: **this run does not settle the objective hypothesis.** It bounds it — CAPABLE_BUT_IDLE is at least
13.7 % in glucose — and names exactly what finishes it.

## What *is* robust here

**1. `INACTIVE_IN_CONDITION` is medium-induced, not structural.** Cross-checked against
`fba_structural_blindspot`: of the 58, **57 carry flux fine once other exchanges are opened**. Exactly 1
is a true reconstruction dead-end — consistent with yesterday's 1/131. These genes need a *condition the
model was never given*, not model repair.

*(I first labelled this class `BLOCKED`, which contradicted my own 1/131 structural finding. Caught by
cross-checking the two artifacts against each other, and renamed.)*

**2. FVA and deletion disagree on 29 genes, and isozymes explain about half.** All 29 are REQUIRED by FVA
yet deletion leaves growth unchanged; **16 of 29 carry an isozyme `or` in their GPR** — the *reaction* is
required, but a paralog covers the *gene*. `b0077`/`b0078` (ilvIH) are the archetype.

That is a real, separable mechanism: FVA reasons about reactions, essentiality is about genes, and the GPR
sits between them. The remaining 13 need another explanation (likely multi-reaction genes with an
alternative route) and are **not** explained here.

## Honest limits

1. **The condition confound above is the dominant one.** Do not quote the 44.3 % as a property of the model.
2. **FVA is per-reaction.** Gene-level CAPABLE_BUT_IDLE means each reaction can *individually* reach zero;
   it does not prove they can be zero *simultaneously*. The deletion cross-check is the gene-level truth.
3. **131-gene, 11-condition panel**, not the full 217-gene set (needs D:).
4. **Isozymes explain 16 of 29 disagreements — 55 %, not all.**

## Next

The finishing move is now precisely specified and cheap once `D:` returns: re-run this classification with
each gene evaluated **in its own experimentally-essential condition**. That removes the confound and gives
the real CAPABLE_BUT_IDLE fraction — the number that decides whether the objective hypothesis survives.

Until then the objective hypothesis stands as **bounded below at 13.7 %, unresolved above** — not
confirmed, and explicitly not the "dominant" story the previous run suggested.
