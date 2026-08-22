# PRE-REGISTRATION — biomass completion from iML1515's own WT equation

**Written:** 2026-08-22, **before the scoring run**, and **frozen on commit**. Amendments must be dated
sections appended here, never silent edits.

**What has already been looked at, stated plainly:** a preflight measured **wildtype feasibility only** —
whether the modified objective still permits growth. It did **not** touch the endpoint (per-gene
essentiality recovery, or the false-positive count). That boundary is the whole point: this morning's
expression-gating run was invalidated by wildtype collapse, so feasibility is checked *before* freezing,
and nothing else is.

Prior context: `wiki/fba_demand_completion_2026-08-22.md` (the diagnostic that named the targets),
`wiki/fba_expression_gated_gpr_result_2026-08-22.md` (the failure this design is shaped by).

## 1 · Hypothesis

34 of the 131 conditionally-essential genes have **no growth effect anywhere**. The diagnostic named,
for 21 of them, a metabolite they are the **sole route** to, and verified 21/21 flip to essential once
that metabolite is demanded.

> **H1.** iML1515's default objective `BIOMASS_Ec_iML1515_core_75p37M` demands a **truncated** LPS
> precursor (`kdo2lipid4_e`) and omits heme O and enterobactin entirely. Restoring those demands **from
> the model's own WT biomass equation** should make the genes that build them essential — without
> inventing any coefficient.

## 2 · Intervention (frozen) — every number READ from the model, none invented

iML1515 ships a second biomass reaction, `BIOMASS_Ec_iML1515_WT_75p37M`, which demands the completed core
LPS, heme O and enterobactin. **Wholesale objective-switching is not available** — preflight found the WT
objective infeasible in **25 of 25** conditions on minimal medium (it additionally demands glycogen,
cardiolipins, spermidine, putrescine, adenosylcobalamin and 27 other components). So the intervention is a
**targeted transplant** of three demands, at the WT equation's own coefficients:

| change | metabolite | coefficient | source |
|---|---|---|---|
| **drop** | `kdo2lipid4_e` | −0.019456 | core biomass (the truncated precursor) |
| **add** | `colipa_e` | −0.008151 | **WT biomass** (completed core LPS) |
| **add** | `hemeO_c` | −0.000223 | **WT biomass** |
| **add** | `enter_c` | −0.000223 | **WT biomass** |

`kdo2lipid4` is dropped because `colipa` is downstream of it; demanding both double-counts the same lipid.

**Preflight result (feasibility only):** wildtype survives in **25 of 25** conditions. Growth changes by
**−0.46 %** (i.e. slightly *up*), because the completed-core coefficient is smaller than the truncated
one. There is no wildtype-collapse artifact available to manufacture recoveries.

## 3 · Per-pathway predictions (frozen) — the sharp part

A bare recovery count would be weak. Each demand targets a specific gene set, so each gets its own
pre-declared prediction:

| demand | genes predicted to FLIP | n |
|---|---|---|
| `colipa_e` | `gmhA` `gmhB` `hldE` `hldD` `waaC` `waaF` `waaP` `waaG` `galU` | 9 |
| `hemeO_c` | `cyoE` | 1 |

| demand | genes predicted **NOT** to flip | n | why |
|---|---|---|---|
| `enter_c` | `fes` `fepB` `fepC` `fepD` `fepG` `tonB` `exbD` | 7 | `enter_c` is **unloaded** enterobactin. These genes are the sole route to the **loaded** complex (`feenter_c`, `fe3dhbzs3`), which no biomass equation demands. Demanding enterobactin should make its *biosynthesis* essential, not its *uptake*. |

**If the 7 iron genes DO flip, the stated mechanism is wrong** even though the headline number would look
better. Recording that now is what stops a bigger number from being read as a better result.

Out of scope: `wecE` / `wzxE` (ECA) — the WT biomass has no ECA component, so there is no
model-supplied coefficient and inventing one is exactly what this design refuses to do.

## 4 · Endpoints (frozen, ordered)

| | endpoint | success |
|---|---|---|
| **Primary** | of the **10** genes predicted to flip, how many become predicted-essential in ≥1 condition where they are experimentally essential | **≥ 8 of 10** |
| **Mechanism check** | of the 7 iron genes predicted NOT to flip, how many flip | **0 expected**; ≥1 ⇒ mechanism disconfirmed |
| Secondary | net change in recall over the 131 | reported, not a bar |
| **Guardrail** | false positives across the full gene × condition grid | **must not rise more than 20 % relative** |

**A primary hit with the guardrail breached is a FAILURE**, not a partial success — the identical rule
that decided this morning's run. Enriching biomass trivially makes more things essential; the guardrail
is what makes the primary mean anything.

## 5 · Determinism gate (frozen)

`processes=1`; the whole pipeline runs **twice**; the two runs must produce **identical** essentiality
calls. Disagreement ⇒ `INDETERMINATE`, never the better of the two.

## 6 · Pre-committed verdicts

| outcome | verdict |
|---|---|
| ≥8/10 flip, 0/7 iron flip, guardrail held, determinism passed | **H1 SUPPORTED** |
| ≥8/10 flip, but ≥1 iron gene flips | **H1 SUPPORTED, MECHANISM DISCONFIRMED** — report both |
| 1–7/10 flip, guardrail held | **H1 PARTIAL** — report; do **not** tune coefficients to reach 8 |
| 0/10 flip | **H1 FALSIFIED** |
| guardrail breached | **FAILURE** regardless of the primary |
| runs disagree | **INDETERMINATE** |

## 7 · Known limits, recorded in advance

1. **This is still a model-side change.** Recovering a gene means the *model* now calls it essential; it
   does not validate the biomass edit against measured biomass composition.
2. **Coefficients are transplanted, not fitted.** That is deliberate — a fitted coefficient would be a
   free parameter, and this design has none.
3. **Growth goes slightly UP** under the modification, so any recovery cannot be attributed to a general
   growth burden.
4. **The 13 no-effect genes with no named demand are untouched** and are expected to stay missed.
5. Recovery is only possible in conditions with labels; the panel is 25 sole-carbon-source conditions.
