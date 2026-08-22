# The hypothesis was exactly right, the mechanism check was exactly right, and it still fails — here is why that matters

**Date:** 2026-08-22 · **Verdict:** `FAILURE_GUARDRAIL_BREACHED` (pre-committed, §6)
**Pre-registration:** `wiki/fba_biomass_completion_prereg_2026-08-22.md` (frozen at `abc8040`, before scoring)
**Artifact:** `wiki/fba_biomass_completion_2026-08-22.json` · **Script:** `scripts/fba_biomass_completion.py`

## Result

| endpoint | pre-registered bar | measured |
|---|---|---|
| **Primary** — predicted-to-flip genes recovered | ≥ 8 of 10 | **10 of 10** |
| **Mechanism check** — iron genes that must NOT flip | 0 of 7 | **0 of 7** |
| Recall over the 131 | reported | 0.394 → **0.451** (TP 721 → 826) |
| **Guardrail** — false positives | ≤ +20 % | **+53.6 % — BREACHED** |
| Determinism | identical calls, 2 runs | **passed**, 0 differing cells |

Both scientific predictions landed **perfectly**. All nine LPS-core genes plus `cyoE` flipped; not one of
the seven iron-uptake genes did — exactly as §3 predicted, because `enter_c` is *unloaded* enterobactin
while those genes are the sole route to the *loaded* complex.

And the pre-committed rule still fails it. That is not a disappointment; it is the finding.

## Why it fails — and this generalises well past this experiment

Every one of the 10 recovered genes is now predicted essential in **25 of 25** conditions, while
experimentally essential in only 5–17:

| gene | predicted essential | experimentally essential | new false positives |
|---|---|---|---|
| `gmhA` | 25 | 16 | 9 |
| `gmhB` | 25 | 17 | 8 |
| `hldE` | 25 | 13 | 12 |
| `hldD` | 25 | 14 | 11 |
| `waaC` | 25 | 13 | 12 |
| `waaF` | 25 | 7 | 18 |
| `waaP` | 25 | 7 | 18 |
| `waaG` | 25 | 5 | 20 |
| `galU` | 25 | 5 | 20 |
| `cyoE` | 25 | 6 | 19 |

**Those 10 genes alone introduce 147 of the 170 new false positives — 86 % of the guardrail breach.** The
cost is not diffuse damage elsewhere; it is concentrated in precisely the genes the intervention
"recovered."

> **A biomass coefficient is condition-independent.** A gene that is the sole route to a demanded
> metabolite *in every condition* therefore becomes **constitutively** essential — never conditionally
> essential. The demand trades a false negative in a few conditions for a false positive in all the rest.

`waaG` is the cleanest case: experimentally essential in 5 of 25 conditions, predicted essential in 25.
Five true positives bought, twenty false positives created. The increment's precision is 105/275 ≈ **38 %**,
well below the baseline model's 721/1038 ≈ **69 %** — so the addition is worse than what it was added to.

**Precise scope of the claim:** this is not "a biomass demand can never yield conditional essentiality."
It could, *if the sole-route property were itself condition-dependent*. For these genes it is not — the
biosynthetic pathway is the only route in every condition. Which means the load-bearing caveat is the one
the diagnostic memo already flagged: **it was run in a single condition**, and single-condition
sole-route is exactly what produces constitutive essentiality.

## What this settles

This is the third structural ceiling measured in this arc, and together they are close to a complete
account of why four independent levers failed identically:

| finding | status |
|---|---|
| 31 of 131 genes unreachable by **any** constraint-based method | proved, validated over ~40k deletions |
| `MIS_CONDITIONED` = 0 — the condition-modelling hypothesis | empty |
| Expression-gated GPR — operator works, **selector** does not | measured, no window exists |
| **Objective completion — right genes, wrong *shape*** | **measured here** |

The residual is not waiting for a better constraint lever, a better threshold, or a richer objective.
Each of those has now been tried and its ceiling measured. **Conditional essentiality needs a
condition-dependent mechanism**, and every intervention tested so far — bounds, thresholds, GPR gating,
biomass demands — is condition-*independent* by construction. That is the through-line.

## Honest limits

1. **Model-side only.** "Recovered" means the model now calls the gene essential. Nothing here validates
   the biomass edit against measured cell composition — and the false-positive result is a reason to
   doubt the edit as stated, not merely its scoring.
2. **No free parameter, deliberately.** Every coefficient was transplanted from iML1515's own
   `BIOMASS_Ec_iML1515_WT_75p37M` and verified against it at run time. Tuning them to pass the guardrail
   was available and is exactly what the pre-registration forbade.
3. **Growth goes slightly *up*** (−0.46 % cost) under the modification, so no recovery is a
   growth-burden artifact.
4. **The wholesale alternative is unavailable** — preflight found the full WT objective infeasible in
   25/25 conditions on minimal medium (it also demands glycogen, cardiolipins, spermidine, adocbl and 27
   others). That was caught *before* freezing, which is why this design is a targeted transplant.
5. **ECA was excluded by design** (`wecE`, `wzxE`) — the WT biomass carries no ECA component, so there
   was no model-supplied coefficient, and inventing one is what this design refuses.

## Next

The one intervention shape not yet tested is a **condition-dependent** one. Everything tried so far is
constitutive. That is a genuinely different design and, on this arc's record, it needs a pre-registration
before it is built — not another lever applied to the same data.
