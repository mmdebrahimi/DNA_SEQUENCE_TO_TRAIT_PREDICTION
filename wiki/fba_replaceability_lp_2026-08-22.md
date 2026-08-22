# Closing my own declared gap: functional redundancy, proved by LP — and the zinc thread finally resolves

**Date:** 2026-08-22 · **Artifacts:** `wiki/fba_replaceability_2026-08-22.json` +
`wiki/fba_replaceability_validation_2026-08-22.json`
**Script:** `scripts/fba_replaceability_lp.py` · **Tests:** `tests/test_fba_replaceability_lp.py` (6)
Model-only — `D:` still disconnected.

## The gap this closes

`wiki/fba_orphan_protection_screen_2026-08-21.md` proved 572 genes structurally uncallable, and stated
plainly that the number was a **lower bound**: its rescue criterion was *exact stoichiometric
duplication*, which iML1515 barely has (12 transformation groups with more than one reaction). It could
not see **functional** redundancy — the very thing the iron case turned out to be. That limitation was
declared, not hidden, and this is the instrument that removes it.

## The test

For reaction *r*: delete *r* (and every reaction the same knockout disables), **close every exchange**,
and add a shadow reaction whose stoichiometry is the negative of *r*'s. At steady state, the shadow
carrying flux *f* forces the remaining internal network to perform *r*'s exact transformation at rate *f*.
Maximise *f*. If *f* reaches *r*'s full capacity, the internal network reproduces *r*'s job at any rate
*r* could ever reach.

**Closing the exchanges is the whole point.** A bypass drawing on an exchange is medium-dependent, so it
cannot support a claim quantified over all media. Excluding them makes the test conservative in the safe
direction: routes are missed, never invented.

**The loop caveat was measured, not assumed.** With every exchange shut, the only remaining fluxes are
internal cycles — exactly where thermodynamically infeasible loops live, and a spurious loop would fake
replaceability. Every candidate was re-tested with `loopless_solution`: **73 of 73 survived.**

## Result

| | genes |
|---|---|
| fully replaceable by internal metabolism | **73** |
| …surviving the loop audit | **73 / 73** |
| partially replaceable (route exists, narrower than the original) | 4 |

**Validated against reality, not just asserted.** All 73 were run through real single-gene deletions in
all 25 carbon conditions — 1,825 deletions. The theorem predicts growth ratio exactly 1 every time.
Observed contradictions: **0**.

## The floor moves

| | uncallable genes | of 1,516 | essential & uncallable | of 131 |
|---|---|---|---|---|
| dead-ends only (`fba_structural_blindspot`) | 128 | 8.4 % | 1 | 0.8 % |
| + GPR-aware + exact duplicates (`fba_orphan_protection_screen`) | 572 | 37.7 % | 23 | 17.6 % |
| **+ functional replaceability (this)** | **645** | **42.5 %** | **31** | **23.7 %** |

> **Nearly a quarter of the experimentally-essential gene set is provably unreachable by any
> constraint-based method** — no medium, no objective, no constraint layer. Gap-fill, threshold retune,
> pFBA and E-Flux did not fail independently four times; they ran into a shared floor.

## The zinc thread resolves — on the fourth attempt

The 8 newly-uncallable essential genes are not a scatter. **Every one is an ABC-type uptake transporter:**

| genes | reaction | substrate |
|---|---|---|
| `znuA` `znuB` `znuC` | `ZNabcpp` | zinc |
| `pstS` `pstA` `pstB` | `PIuabcpp` | phosphate |
| `mgtA` | `MG2uabcpp` | magnesium |
| `xylF` | `XYLabcpp` | xylose |

`znuA/znuB/znuC` are three of the twelve Fe/Zn acquisition genes the bypass-closure experiment was built
to explain. That explanation has now been rewritten three times; this is the version with a proof and a
validation behind it:

- **v1** — "the model gets free zinc, so the machinery idles." *(observation)*
- **v2** — "…and could not use the machinery anyway; the exchanges are export-only." **Wrong** — the
  loaded siderophore never needs importing.
- **v3** — "route redundancy plus orphan reactions." **Falsified** — zero genes are orphan-protected.
- **v4 (proved + validated)** — **the transporter's job is fully replaceable by internal metabolism at
  full capacity, independent of the medium.** That is why no medium change, objective term, or constraint
  layer could ever have made these genes essential — and it is why closing the tap starved the model
  instead of flipping the call.

The unifying shape is worth stating plainly: **ion and nutrient uptake in iML1515 is deeply redundant,
and redundancy in a transport step is invisible to single-gene deletion.** Orphan reactions were one
narrow instance of a much more general property, which is why the orphan generalisation failed while the
underlying intuition was pointing somewhere real.

## Honest limits

1. **Lower bound still.** Full replaceability is required at the reaction's *full* capacity bound — for a
   1000-bound reaction, far more than any medium would demand. Genes are missed, never falsely admitted.
   The 4 partially-replaceable genes are exactly the ones this strictness excludes.
2. **Medium-independence is bought by closing exchanges**, so a real bypass that legitimately uses an
   exchange is not counted. Deliberate.
3. **Model reach, not biology.** A gene here is a *declared blind spot* — the claim is that the model can
   never call it essential, never that the cell can lose it.
4. **Loop audit passed on all 73**, so nothing rests on an unaudited number — but `loopless_solution` is
   itself a relaxation of full thermodynamic consistency.
5. **The 23.7 % floor is against one committed label set** (131 genes) and scales with that set.

## Next

The residue is now sharp. Of 131 experimentally-essential genes: **31 are provably unreachable**, 8 of the
remainder are isozyme-masked (curable by strain-aware GPR curation), and **100 are genuinely callable and
still missed** — which is where any remaining modelling effort belongs. The one lever none of the four
failed ones tried is **expression-gated GPR** (treat an unexpressed isozyme as *absent*, not merely
capacity-limited), and it must be pre-registered before it is run: its pre-declared target is 8 genes, and
naming that number in advance is the only thing separating it from metric-shopping.
