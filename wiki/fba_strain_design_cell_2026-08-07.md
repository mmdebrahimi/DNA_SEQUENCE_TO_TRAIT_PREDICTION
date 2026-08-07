# The design cell — from "what does this edit do?" to "which edits give me this product?" (2026-08-07)

Track A of the design epoch (`wiki/design_epoch_plan_2026-08-07.md`). The FBA cell could already answer
*edit → trait*. This adds the **inverse**, which is what strain engineering actually needs:

> given a PRODUCT, which knockouts make producing it **necessary for the cell to grow**?

That property is **growth coupling**. It is the thing that makes an engineered strain stable: if the
product can only be made while growing, then selection for growth becomes selection for production
instead of fighting it. Ships as `dna_decode/fba/design.py` + `scripts/fba_strain_design.py` +
`dna-decode fba --design-target`.

## The test

At a fixed growth floor, minimise **and** maximise the product exchange:

- `min_flux > 0` → **OBLIGATORY**: *every* flux distribution that grows also secretes the product. A design.
- `min_flux ≈ 0 < max_flux` → **POSSIBLE**: the cell can make it and can avoid it. The un-engineered case.
- `max_flux ≈ 0` → **INFEASIBLE**.

## Validation: it recovers a published design

**End-to-end gate (`tests/test_fba_design.py::test_reproduces_the_literature_anaerobic_succinate_design`):**
the search independently finds the OptKnock-lineage anaerobic succinate design — remove the competing
fermentation routes and succinate becomes the obligatory redox sink:

| knockouts | guaranteed succinate | wild-type floor | gain | growth |
|---|---|---|---|---|
| `PFL` + `LDH_D` + `ALCD2x` | **9.26** | 0.047 | **+9.22** | 0.082 /h |

That is a ~196× increase in the *guaranteed* floor, not the ceiling.

## Four defects the verification caught — each would have shipped a broken tool that looked fine

This is the whole reason the run took as long as it did. Every one was found by inspecting real output,
not by a failing test.

**1. "Is OBLIGATORY" ≠ "is a design" — reported 2096 of 2096 knockouts as designs.**
Anaerobic *E. coli* already secretes a little succinate obligatorily, so the **wild type is OBLIGATORY
before any edit**. Counting the verdict rather than the *improvement* made every single evaluated
knockout a "growth-coupled design." Fixed by `improves_on_baseline` — a design must **beat** the
wild-type floor. The artifact now reports `wildtype_already_coupled` and
`n_inheriting_baseline_only` so this can never hide again.

**2. Gene-level knockouts are blunted by GPR isozymes.**
`PFL`'s rule is `(b0902 and b0903) or (b0902 and b3114) or (b3951 and b3952) or ...`, so deleting `pflB`
leaves the reaction fully open — verified: bounds stayed `(0, 1000)` and anaerobic growth was unchanged
at 0.15754 /h across *four* different knockout sets. Default is now **reaction-level** (`--level`).

**3. The growth floor was relative to the wild type, and far too low.**
Coupling is defined at **near-optimal growth of the strain itself**. A wild-type-relative floor is
unreachable for any knockout that slows the cell, and a low floor leaves so much slack that `min_flux`
collapses to ~0 for everything. Same strain, same model, only the floor differs:
`PFL+LDH_D+ALCD2x` guarantees **0.0027 at 10%** of the mutant optimum and **12.28 at 99%**.
Default is now 90% of each strain's own maximum.

**4. The top-ranked design was a pseudo-reaction.**
The first full search returned `ATPM` — the non-growth-associated **ATP maintenance** requirement — as
its best succinate design. You cannot delete a maintenance requirement with a knockout; there is no gene
for it. Fixed by requiring a non-empty GPR, which excludes in one principled rule every non-constructible
candidate: exchanges, demands, sinks, spontaneous reactions, and maintenance.

Defects 2, 3 and 4 all failed in the *same* direction — they made the tool return **zero or garbage
designs while appearing to work**. A green test suite would not have caught any of them.

## Honest limits

- **The pair/triple search is a bounded heuristic, not exhaustive.** Pairs and triples are drawn from the
  best-ranked single knockouts, and the real design's three members are individually unremarkable — so
  the *unrestricted* search can miss a design that the tool finds immediately when the candidate set is
  scoped to a pathway. Reported in every artifact's `search` field. A MILP (OptKnock proper) would fix
  this; `cobra.design` was removed from cobrapy, so that is a future build.
- **Every design is a HYPOTHESIS FOR THE BENCH, never a validated strain.** FBA sees stoichiometry. It
  does not see regulation, enzyme kinetics, toxicity, metabolic burden, or whether the knockout strain is
  constructible. This is the same wall that separates the decoders from clinical claims, and it is
  stamped into every emitted record's `scope`.
- **Aerobic glucose yields nothing for succinate** (0 designs over 1320 singles + 780 pairs), which is
  correct: with respiration available the cell can fully oxidise substrate and is never forced to secrete.
  Coupling lives in the anaerobic/micro-aerobic regime.

## Run

```bash
uv run python scripts/fba_strain_design.py --target succ --anaerobic --max-knockouts 3
dna-decode fba --design-target succ --anaerobic --json
```
