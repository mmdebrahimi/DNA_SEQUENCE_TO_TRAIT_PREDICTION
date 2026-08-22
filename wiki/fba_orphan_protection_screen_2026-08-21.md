# The orphan-protection screen: my own hypothesis does not survive it — but the screen finds a harder floor

**Date:** 2026-08-21 · **Artifact:** `wiki/fba_orphan_protection_2026-08-21.json`
**Script:** `scripts/fba_orphan_protection_screen.py` · **Tests:** `tests/test_fba_orphan_protection.py` (22)
Model-only — `D:` still disconnected. **Supersedes the generalisation in**
`wiki/fba_orphan_redundancy_2026-08-21.md`.

## The hypothesis, and the verdict

Yesterday's memo generalised the iron finding to: *"any pathway whose function is duplicated by an orphan
reaction is unconditionally protected from gene-deletion analysis."* It named the screen that would turn
that count into an impact estimate. The screen is now built, and it says:

> **Zero of 1,516 genes in iML1515 are protected by an orphan duplicate. Not a small number — zero.**

The statement remains *true* (an empty-GPR reaction genuinely cannot be deleted) and *vacuous* at gene
level. The iron case was **functional** redundancy — `FE3abcpp` moves Fe(III) with ATP where `FE2tpp`
symports Fe(II), a different transformation entirely — not **stoichiometric** duplication. One real
instance did not generalise into a mechanism, and I should not have written it as though it had.

## What the screen actually is

For each gene, take the reactions its knockout **actually** disables (full-GPR aware, so an isozyme `or`
keeps its reaction alive). The gene is **structurally uncallable** if every one of those reactions is:

- **(B) universally blocked** — no flux even with all exchanges open. Any medium is a *subset* of
  all-open, so it is blocked in every medium; deleting it removes nothing.
- **(C) duplicate-rescued** — a surviving, unblocked reaction performs the same transformation
  (identical stoichiometry up to a scalar, reversal included) with a covering flux interval.

Then the knockout leaves the feasible set over the retained reactions unchanged, so growth is unchanged
**for any objective, in any medium** — and the model can never call the gene essential. The empty case
(nothing disabled — fully isozyme-buffered) is the degenerate instance.

**This is a theorem, and it was tested against the thing it predicts.** A full single-gene deletion ran
over all 1,516 genes in all 25 carbon conditions — 37,900 deletions. Predicted violations: 0. Observed:
**0**. And the screen is not vacuous: in every condition, 196–209 of the 944 `CALLABLE` genes *are* called
essential, so the split separates something real.

## Results

| class | genes | share |
|---|---|---|
| CALLABLE | 944 | 62.3 % |
| NO_KO_EFFECT (fully isozyme-buffered) | 451 | 29.7 % |
| ALL_DISABLED_BLOCKED | 120 | 7.9 % |
| DUPLICATE_RESCUED | **1** | 0.1 % |
| **structurally uncallable** | **572** | **37.7 %** |

The single duplicate-rescued gene is `entC`, rescued by `ICHORS_copy2` — a *gene-carried* copy, not an
orphan. iML1515 has only **12** transformation groups containing more than one reaction; exact
stoichiometric duplication is essentially absent, which is why (C) contributes almost nothing.

*(Count correction: yesterday's "113 orphan non-exchange reactions" excluded the 2 biomass
pseudo-reactions without saying so. Precisely: **446** orphan reactions, **115** non-exchange, **113**
excluding biomass, **107** excluding biomass and the 6 demand reactions.)*

## The number that matters: a hard floor on every constraint-based method

Cross-tabbed against the committed experimental calls (`fba_eflux_bridge_2026-08-17_ratios.json`, 131
conditionally-essential genes):

> **23 of 131 (17.6 %) experimentally-essential genes are structurally uncallable.** No medium, no
> objective, no constraint layer — gap-fill, threshold retune, pFBA, E-Flux, or anything else — can ever
> move them.

That is a large piece of the shared-failure puzzle. It is also **23× more than the previous structural
account**: `fba_structural_blindspot.py` found only **1 of 131** by looking for dead-end reactions. The
difference is entirely GPR-awareness — 17 of the 23 are hidden by gene-level redundancy, which a
reaction-level analysis cannot see.

| mechanism | genes |
|---|---|
| NO_KO_EFFECT (GPR redundancy) | 17 |
| ALL_DISABLED_BLOCKED (dead-end) | 6 — `aspC` `dsbB` `cysA` `cysW` `cysU` `xylB` |

## Splitting the 17 — and rejecting the tidy story

The tempting conclusion is "the model carries isozymes the strain doesn't have; curate the GPRs." Before
publishing that I ran the unmask test: delete the gene **and every gene sharing its reactions**, then
re-optimise. If the function is genuinely essential, growth collapses.

**It holds for only 8 of the 17.**

| verdict | genes |
|---|---|
| **isozyme-masked** — function IS essential, only redundancy hides it | `ilvI` `ilvH` `ilvB` `ilvN` `aroF` `tktA` `trxA` `ompC` |
| **not masked** — function dispensable in the model regardless | `pgm` `sapD` `ndk` `trkA` `gntT` `glpD` `rbsD` `trkH` `fbp` |

The masked set is biologically legible: `ilvI/ilvH` (AHAS III) and `ilvB/ilvN` (AHAS I) are isozyme pairs
behind the same two reactions; `aroF` sits behind `aroG`/`aroH`; `tktA` behind `tktB`; `ompC` behind
`ompF` and seven other porins. The other 9 are not masked at all: their function is dispensable in the
model no matter what, which is an objective/network problem that GPR work will not touch.

> **Correction (2026-08-22).** An earlier version of this paragraph said the third AHAS isozyme
> (`ilvG`/`ilvM`) is frameshifted in K-12 and that "the model does not know that". **Both halves are
> wrong.** `ilvG` and `ilvM` are **absent from iML1515 entirely** — the model already excludes them — and
> the measured GPR of `ACHBS`/`ACLS` is `(b3670 and b3671) or (b0077 and b0078)`, i.e. AHAS I *or*
> AHAS III, two isozymes that are both genuinely functional in K-12. I asserted that from memory instead
> of reading the model.
>
> The correction *strengthens* rather than weakens the follow-on: if the redundancy is between two real
> isozymes, then an experimentally-essential `ilvI`/`ilvH` implies its partner is not **expressed** in the
> tested condition — which is exactly what an expression-gated GPR would capture, and not something a
> sequence-level curation fix would address.

**Method note that changed the answer.** My first unmask test zeroed the gene's *reactions* instead of
deleting the masking *genes*. That mis-scored both directions — it wrongly cleared `trxA` (10 reactions,
partly non-redundant) and wrongly credited `ompC`, whose 285 mapped diffusion reactions make
reaction-zeroing a test of "remove all outer-membrane transport", not of redundancy. Deleting the mask is
the only faithful version, and it is what the script does. `tests/test_fba_orphan_protection.py` pins the
distinction.

## Where the 131 now stand

| | genes | what would move it |
|---|---|---|
| structurally uncallable — isozyme-masked | 8 (6.1 %) | strain-aware GPR curation |
| structurally uncallable — dead-end reaction | 6 (4.6 %) | reconstruction repair |
| structurally uncallable — function dispensable | 9 (6.9 %) | objective / network realism |
| callable, still missed | 108 (82.4 %) | open — the model *could* call these |

## Honest limits

1. **Model reach, not biology.** A gene in the uncallable set is a *declared blind spot*, not a claim that
   it is dispensable in the cell.
2. **Medium means exchange bounds.** Changing internal bounds is outside the theorem.
3. **The rescue criterion is exact stoichiometry** — conservative by construction. It misses *functional*
   redundancy (the iron case), so 572 is a **lower bound** on the true uncallable set. Proving functional
   redundancy would need an LP-based test, not a stoichiometric one; I have not built one and do not claim
   the gap is small.
4. **The unmask test is single-condition** (glucose) and its 85/451 population figure is not a claim about
   the other 24.
5. **The 17.6 % floor is against one committed label set.** It scales with that set, not with all of
   *E. coli* essentiality.

## Next

The screen turned a hypothesis into a number by killing the hypothesis, and the residue is sharper than
what it replaced. The natural follow-on is the one lever none of the four failed ones tried:
**expression-gated GPR** — treat an isozyme as *absent* (GPR boolean false), not merely capacity-limited,
when its expression is below threshold. This is precisely why E-Flux could not have worked here: E-Flux
scales a reaction's *bounds*, but an `or` of two isozymes keeps the reaction functional no matter how the
bounds move, so a single-gene deletion is untouched by it. Gating the boolean is a different intervention
with a pre-declared target of 8 genes — small, but for the first time *named in advance* rather than
hoped for. It must be pre-registered before it is run.
