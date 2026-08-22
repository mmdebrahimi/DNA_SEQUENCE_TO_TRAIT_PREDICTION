# Expression-gated GPR: the primary endpoint was MET, and the guardrail is what saved it from being published

**Date:** 2026-08-22 · **Verdict:** `FAILURE_GUARDRAIL_BREACHED` (pre-committed, §6 row 4)
**Pre-registration:** `wiki/fba_expression_gated_gpr_prereg_2026-08-22.md` (frozen at `2e95b7d`, + Amendment 1)
**Artifact:** `wiki/fba_expression_gated_gpr_2026-08-22.json` · **Script:** `scripts/fba_expression_gated_gpr.py`
**Tests:** `tests/test_fba_expression_gated_gpr.py` (9) · Ran on real PRECISE-1K + `feba.db` (`D:` reattached)

## The headline

| | |
|---|---|
| Pre-registered primary (≥4 of 8) | **5 of 8 — MET** |
| Recall over the 131 | **0.338 → 0.648, nearly doubled** |
| False-positive guardrail (≤ +20 %) | **+334.6 % — BREACHED** |
| Honest primary, growing conditions only | **0 of 8** |
| Determinism (§5) | **passed**, 0 differing cells, both runs, all 4 arms |

**If §4 had carried only the recovery count, this would have shipped as a success.** It is not one.

## What actually happened

A genome-wide percentile gate makes the **wildtype itself infeasible**. At the primary 20th percentile the
model stops growing in 4 of 11 conditions — and every scoring path in this repo codes a zero-growth
wildtype as *all genes essential*. That is where the "recoveries" came from.

The correspondence is exact, with no exceptions:

| gate | genes marked absent | wildtype collapsed in | recoveries occurred in |
|---|---|---|---|
| p10 | 151 | **0 of 11** | **0 conditions** |
| **p20 (primary)** | 302 | D-Galactose, D-Ribose, Glycerol, K-acetate | **the same 4** |
| p30 | 453 | those + D-Fructose, D-Sorbitol | **the same 6** |

Every single recovery event — 8 at p20, 12 at p30 — is in a collapsed condition. **Restricted to
conditions where the model can still grow, the recovery count is 0 of 8 at every gate level.**

So both pre-committed rows point the same way: the guardrail says `FAILURE`, and the growth-feasible
restriction independently lands on `H1 FALSIFIED`.

## The operator works. The selector does not.

This is the part worth carrying, and it is not "boolean gating can't work".

A positive control run before the experiment shows the mechanism is real and clean:

| gated off | `ilvI` (b0077) | `ilvB` (b3671) | wildtype |
|---|---|---|---|
| nothing | 1.000 | 1.000 | 0.876997 |
| AHAS I (`ilvN`+`ilvB`) | **0.000** | — | **0.876997 (unchanged)** |
| AHAS III (`ilvI`+`ilvH`) | — | **0.000** | **0.876997 (unchanged)** |
| `tktB` | `tktA` → **0.000** | | unchanged |
| `aroG`+`aroH` | `aroF` → **0.000** | | unchanged |

Gating the *right* partner makes the target essential and costs the wildtype **nothing**. H1's mechanism
is correct. What fails is the **selector**:

- at **p10** the gate is **inert** — false positives move by +0.0 %, nothing recovers, nothing collapses;
- at **p20/p30** it is **lethal** — it removes essential metabolism along with redundant isozymes.

**There is no window in between.** A genome-wide expression percentile cannot separate "this isozyme is
off" from "this pathway is needed", because ~20 % of genes being below the 20th percentile is true *by
construction* whether or not any of them is genuinely absent.

And the way out is circular: to gate only the redundant isozymes you would have to already know which
isozyme is off in that condition — which is the thing the experiment was trying to infer.

## Why the pre-registration earned its keep

This project retracted a positive result on 2026-08-17 for choosing an endpoint after the primary
disappointed. Here the opposite happened, and the frozen document did the work:

1. **The guardrail decided the experiment, not the primary.** It was written before any data was
   reachable, precisely because gating expression off trivially makes more things look essential.
2. **The target set was frozen** to 8 named genes from a committed artifact — `frozen_target_set()`
   *aborts* if that artifact ever changes size, so a regenerated upstream file cannot silently re-score a
   different endpoint.
3. **The sensitivity range {10, 20, 30} was declared in advance** and all three are reported. p30 gives
   the best-looking recall (0.755) and is not the headline, because reporting best-of was forbidden in
   writing.
4. **Determinism was gated before interpretation** — two full pipeline runs, 0 differing cells across all
   4 arms.
5. The one post-hoc addition (growth-feasible restriction) is recorded as **Amendment 1** and moves the
   result **5/8 → 0/8** — it strengthens the negative. That direction is what makes an unplanned analysis
   safe.

## Honest limits

1. **11 of 25 conditions** have matched PRECISE-1K expression; a target cannot be recovered where there is
   no expression data.
2. **Strain mismatch** — PRECISE-1K is K-12 MG1655, the labels are `Keio` (BW25113 parent). Inherited
   from the E-Flux bridge, unchanged.
3. **Low mRNA does not prove no protein.** The claim is about a gating rule, never about the cell.
4. **This falsifies the percentile selector, not every boolean-gating scheme.** A targeted selector —
   operon-level, regulon-informed, or one that protects the wildtype by construction — is untested, and
   would need its own pre-registration.

## Where this leaves the deficit

Of 131 experimentally-essential genes:

| | genes | status |
|---|---|---|
| provably uncallable by any constraint-based method | 31 | closed — proved, validated |
| …of which isozyme-masked | 8 | **expression-gating does not recover them (this memo)** |
| callable, still missed | 100 | the remaining frontier |

The 8 masked genes are now the *best-characterised* failures in the set: the mechanism that hides them is
understood exactly, the intervention that would unmask them is known to work in principle, and the only
available selector for it has been measured and rejected.
