# Selector preflight: the TRN is available, and two of my own claims are wrong

**Date:** 2026-08-22 · Preflight only — **feasibility and wildtype behaviour, never the endpoint.**
Follows the adversarial review that retracted the closure claim
(`wiki/fba_biomass_completion_result_2026-08-22.md` → RETRACTION).

The review named regulatory-network-constrained FBA as the untried in-family method. Before
pre-registering anything, three things needed checking. All three came back informative, and two of them
contradict claims I published earlier today.

## 1 · Feasibility: GO — the data exists and covers the targets

| check | result |
|---|---|
| RegulonDB direct | **unreachable** (HTTP 000) |
| SBRG PRECISE TRN (RegulonDB-derived, GitHub) | **fetchable** — `data/raw/ecoli_trn/TRN.csv` |
| size | **8,325 interactions / 2,797 genes / 237 regulators**, keyed by b-number |
| coverage of the 131 conditionally-essential genes | **121 (92 %)** |
| **the 8 isozyme-masked targets whose partner isozyme is in the TRN** | **8 of 8** |

Independence: the TRN is curated transcriptional regulation; the labels are RB-TnSeq fitness. No shared
provenance, so a TRN-based selector is not circular against the thing being predicted.

## 2 · My published diagnosis of the wildtype collapse was WRONG

`wiki/fba_expression_gated_gpr_result_2026-08-22.md` attributes the collapse to volume: *"a gate at the
Nth percentile marks ~N % of measured genes absent BY CONSTRUCTION"*. The obvious fix would then be a
**per-gene** threshold (is this gene unusually low *for itself*, against its own distribution across all
1,035 compendium samples) instead of a genome-wide one.

That fix does cut the volume sharply — and **the wildtype still collapses in 6 of 11 conditions**:

| condition | genes marked off (per-gene p10) | wildtype |
|---|---|---|
| D-Galactose | **7** | **0.000 — collapsed** |
| Potassium acetate | 8 | collapsed |
| Sodium pyruvate | 9 | survived |
| D-Ribose | 230 | collapsed |
| D-Glucose | 0 | survived |

**Seven genes collapse galactose.** So the collapse was never primarily about how many genes are gated —
it is about *which*. Some genes the model routes flux through are lowly expressed in that condition, and
gating them removes capability the model needs even though the real cell demonstrably grows there.

That is a sharper statement than the one I published, and it changes what a working selector must do:
not "gate fewer genes" but "never gate away capability."

## 3 · "Safe by construction" was also wrong — single-gene safety does not compose

The orphan-protection theorem says a gene whose knockout disables **no** reaction cannot change the
feasible set. 451 of 1,516 genes qualify. Restricting the gate to those genes should therefore leave the
wildtype untouched *by construction*.

It does — in 9 of 11 conditions, exactly as predicted. And it fails in two:

| condition | low & "safe" genes gated | wildtype |
|---|---|---|
| D-Gluconate | 29 | **collapsed** |
| D-Ribose | 60 | **collapsed** |
| the other 9 | 0–27 | **unchanged to 5 decimal places** |

The reason is a genuine composition error in my reasoning, not a bug: `gpr_disabled_reactions(g)` is a
**single-gene** property. Gate one member of an isozyme pair and the reaction survives; gate **both**
members simultaneously and it does not. Redundancy is pairwise, and safety proved one-gene-at-a-time does
not compose over a set.

> **Reusable:** a per-element safety property does not license applying it to a whole set at once. The
> composition has to be checked jointly — or constructed to hold jointly.

## What follows

The requirement for a working selector is now precise and constructive, rather than "find a better
threshold":

1. select candidate genes by a **per-gene** expression criterion (not genome-wide);
2. restrict to genes whose knockout disables nothing (necessary, not sufficient);
3. **verify jointly** that the selected SET disables no reaction, and resolve collisions (an isozyme pair
   both selected) before applying;
4. the wildtype must then be provably unchanged — as a *constructed invariant*, not an assumption.

Step 3 is what both failures above were missing. The pre-registration for the run is at
`wiki/fba_composed_selector_prereg_2026-08-22.md`.

## Honest limits

1. **Preflight only.** Nothing here measured recovery of any target gene, and no endpoint was computed.
2. **11 of 25 conditions** have matched expression; the collapse counts are over those 11.
3. **The per-gene reference is the full 1,035-sample compendium**, which mixes many stresses and strains;
   a condition-matched reference would be tighter and is not what was used.
4. **Low mRNA still does not prove absent protein** — unchanged caveat, and it is a candidate explanation
   for finding 2 that this preflight cannot separate from "the model routes through the wrong gene".
