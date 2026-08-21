# Route redundancy + orphan reactions — why a whole pathway can be permanently undeletable

**Date:** 2026-08-21 · Model-only (D: disconnected) · **Supersedes the explanation in**
`wiki/fba_bypass_closure_2026-08-21.md` (its INCONCLUSIVE verdict stands; its stated cause did not)

## What I got wrong first

The bypass-closure memo explained the iron result as: *every ferri-siderophore exchange is export-only,
so no loaded siderophore can be imported, so the machinery is structurally unable to deliver iron.*

True about the exchange bounds, **wrong as the cause**. The loaded siderophore never needs importing
from the medium:

| reaction | stoichiometry | genes |
|---|---|---|
| `FEENTERexs` | `enter_e + fe3_e → feenter_e` | spontaneous |
| `FEENTERtonex` | `feenter_e + h_p → feenter_p` | **fepA, tonB, exbB, exbD** |
| `FEENTERabcpp` | `feenter_p → feenter_c` | **fepB, fepC, fepD, fepG** |
| `FEENTERR1–4` / `FEENTERES` | `feenter_c → fe2_c` | fes, yqjH, spontaneous |

The cell secretes enterobactin, it chelates extracellular Fe(III) **outside**, and the complex comes back
in through the fep/ton machinery. The route is present, gene-associated, and traversable.

The real reason closing `EX_fe2_e`/`EX_fe3_e` starves the model is simpler: **`fe3_e` is the shared input
to both routes.** Closing the exchange kills the alternative along with the shortcut. My "close the tap"
design could never have worked for iron.

## What is actually going on

Iron is drawn at **0.014085 mmol/gDW/h** — a tiny demand — and iML1515 offers many ways to satisfy it.
Periplasm→cytosol alone:

| reaction | genes |
|---|---|
| `FE2tpp` | zupT |
| `FE2abcpp` | feoB |
| `FE2t2pp` | mntH |
| **`FE3abcpp`** | **none — empty GPR** |

Measured, blocking these in sequence (glucose minimal, everything else default):

| blocked | growth | siderophore flux (`FEENTERtonex`) |
|---|---|---|
| nothing | 0.876997 | 0.000000 |
| the 3 gene-associated Fe2 routes | 0.876920 | 0.000000 |
| those 3 **plus** the orphan `FE3abcpp` | **0.876880** | **0.000000** |

Blocking all four costs **0.013 % of growth**, and the siderophore machinery **still never activates** —
so at least a fifth route exists that I did not map. I stopped mapping; the point was already made.

## The mechanism, stated generally

Two properties combine, and the second is the one worth carrying:

1. **Route redundancy.** A cheap metabolite with many parallel entry routes cannot be made limiting by
   deleting any one of them.
2. **Orphan reactions.** `FE3abcpp` has an **empty GPR**. A reaction with no gene association can never
   be disabled by *any* gene deletion, in any condition, under any GPR-aware method. **iML1515 has 113
   such non-exchange reactions (4.2 % of 2,712).**

> **Any pathway whose function is duplicated by an orphan reaction is unconditionally protected from
> gene-deletion analysis.** No demand term, no medium change, and no constraint layer alters that — the
> knockout simply has nothing to switch off.

That is a *third* explanation class, alongside "medium realism" and "objective incompleteness", and it is
the one that fits the 12 iron/zinc acquisition genes.

## Why this matters beyond iron

It gives a **model-only, label-free screen** for genes that FBA can never call essential: a gene is
unconditionally safe from deletion if every reaction it uniquely enables has a parallel route that is
either orphan or carried by a different gene. That is computable from the reconstruction alone — no
`feba.db`, no experimental labels — and it would let the decoder *declare* the blind spot up front rather
than silently mispredicting it.

It also reframes the earlier free-iron observation for the third time, and this version is the measured
one:

- v1: "the model gets free iron, so the siderophore system idles."
- v2 *(wrong)*: "...and could not use that system even if it didn't."
- **v3 (measured): the system is fully usable — it simply never has to be used, because iron is cheap,
  demand is tiny, and at least one alternative route answers to no gene at all.**

## Honest limits

1. **I did not map every iron route.** Growth survives blocking four, so a fifth exists. The claim is
   "redundancy is deep", not "there are exactly N routes".
2. **113 orphan reactions is a count, not an impact estimate.** How many *genes* are protected by an
   orphan alternative is not computed here — that is the screen described above, and it is not built.
3. **One condition, one nutrient class.** Zinc was only tested jointly with iron.
4. **Model-structure claim, not a recovered true positive.** It explains why these genes can never be
   predicted essential; it does not improve any score.

## Next

Build the **orphan-protection screen**: for each gene, does every reaction it uniquely enables have a
parallel route that no deletion can close? Output is a declared, label-free list of genes the model
cannot call essential — which is exactly the abstention architecture this repo already ships elsewhere.
Model-only, so `D:` does not block it.
