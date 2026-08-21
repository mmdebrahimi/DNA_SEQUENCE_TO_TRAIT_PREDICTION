# The structural blind spot is real — and it does NOT explain the deficit

**Date:** 2026-08-21 · **Verdict:** my own hypothesis **FALSIFIED**, with a redirect
**Artifact:** `wiki/fba_structural_blindspot_2026-08-21.json` · **Script:** `scripts/fba_structural_blindspot.py`

The `/innovate` run named this as the best next move: make the class-A "no flux" finding
**basis-independent** instead of a lower bound, and split *blocked in the reconstruction* from *merely
unused*. Both imply different fixes. It is model-only, so it ran with `D:` disconnected.

## What was measured

Blocked reactions computed with **every exchange open**, so "blocked" means a dead-end in the
reconstruction itself, not starvation by the current medium:

| | count | share |
|---|---|---|
| structurally blocked reactions | 260 / 2,712 | 9.6 % |
| genes lying **entirely** on blocked reactions | **128 / 1,516** | **8.4 %** |

Those 128 genes **can never be called essential** by any FBA variant, in any condition, under any
constraint layer. That is a permanent, provable blind spot in the model — not a tuning problem.

They are also a coherent *kind* of gene: `copA`, `cusA/B/C/F` (copper efflux), `ybdG` (mechanosensitive
channel), `nfsB` (nitroreductase), `entH` (enterobactin hydrolase). Efflux, detox and stress functions —
exactly what a growth-maximising stoichiometric model has no demand for. This is the same boundary the
stress axis hit from the other side.

## The falsification

I expected this set to explain a large share of the conditional-essentiality deficit. It does not.

Intersecting the 128 blocked genes with the **131 experimentally conditionally-essential genes** from the
committed E-Flux panel:

> **1 gene. 0.8 % of the essential panel.**

**Blocked reactions are not the explanation.** The appealing story — "the model misses these genes because
they are dead-ends in the reconstruction" — is wrong by two orders of magnitude.

## What that redirects to

Class A (zero flux at the optimum, 25.1 % of true-essential cells) is therefore dominated by genes that
are **unblocked but idle**: the reaction *can* carry flux, and simply doesn't have to at the optimum.

That points squarely at the **objective**, not the reconstruction. If a reaction is capable of flux but
carries none while its gene is experimentally essential, then biomass is not demanding the product. The
fix is a condition-specific demand/maintenance term, not a repaired dead-end and not another constraint
layer — consistent with gap-fill, threshold retuning, pFBA and E-Flux all failing identically.

## Honest limits

1. **The essential set is the 131-gene, 11-condition E-Flux panel**, not the full 25-condition/217-gene
   set. The larger set is behind `feba.db` on the disconnected `D:`. A 1/131 overlap is decisive enough
   that a larger denominator will not rescue the hypothesis, but the exact figure will move.
2. **"Entirely on blocked reactions" is strict** by design. A gene with one blocked and one live reaction
   is excluded — correctly, since the live reaction can still make it matter.
3. **This did not deliver the basis-independent class A.** Blocked-ness is a *stronger* condition than
   "zero flux is attainable at the optimum"; the full FVA-at-optimum version is still not done, and it is
   what would convert class A from a lower bound into an exact set.
4. **The efflux/detox reading of the 128 is an observation on a sample of names**, not a systematic
   functional enrichment test.

## Next

The unblocked-but-idle finding makes the objective the live suspect. The cheap model-only follow-up is
FVA-at-optimum over the class-A reactions (is zero flux *attainable*, not merely *observed*), which
finishes the job this run started. The label-side candidates from the `/innovate` ledger unblock when
`D:` returns.
