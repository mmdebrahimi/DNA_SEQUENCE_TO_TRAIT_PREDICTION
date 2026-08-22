# "Objective incompleteness" was a label — here is the measurement, per gene

**Date:** 2026-08-22 · **Artifact:** `wiki/fba_demand_completion_2026-08-22.json`
**Script:** `scripts/fba_demand_completion_probe.py` · **Tests:** `tests/test_fba_demand_completion.py` (6)
Model-only. **DIAGNOSTIC ONLY — no recall claim is made here** (see *Scope* below).

## What was asked

The partition left **34 genes with no growth effect anywhere** — deleting them changes growth by nothing,
in all 25 labelled conditions, and they are not in the proved-unreachable set. I called that "objective
incompleteness." That was a name, not a measurement. This asks the constructive question instead:

> **Which missing biomass demand would make each of these genes essential?**

For gene *g* and each non-currency metabolite *m* its reactions touch: maximise a demand on *m* subject to
biomass ≥ 10 % of wildtype, **with and without *g***. If the model can make *m* normally but cannot make
it at all without *g*, then *g* is the sole route to *m*.

## Result: 21 of 34 have a named missing demand — and all 21 verify

**"Sole route" only *implies* "essential once biomass demands it."** So the script forces each named
demand and re-runs the deletion rather than asserting the implication:

> **21 of 21 flip to essential under their own demand.** And the demands are *cheap*: they cost between
> **0.11 % and 2.8 %** of wildtype growth (0.000965–0.0247 on wt 0.876997).

The 21 are not scattered. They are five coherent pathways, recovered from the reconstruction with no
prior knowledge supplied:

| pathway | genes | the missing metabolite(s) |
|---|---|---|
| **LPS inner core** | 9 — `gmhA` `gmhB` `hldE` `hldD` `waaC` `waaF` `waaP` `waaG` `galU` | `gmhep7p` → `gmhep1p` → `adphep_DD` → `adphep_LD` → `hlipa` → `hhlipa` → `phhlipa` → `gicolipa` |
| **Enterobactin / iron uptake** | 7 — `fes` `fepB` `fepC` `fepD` `fepG` `tonB` `exbD` | `feenter_c`, `fe3dhbzs3`, `fe3dcit_p` |
| **Enterobacterial common antigen** | 2 — `wecE` `wzxE` | `dtdp4addg`, `unagamuf_p` |
| **Heme O** | 1 — `cyoE` | `hemeO` |
| **Regulatory second messengers** | 2 — `relA` `cyaA` | `gdptp` (ppGpp precursor), `camp` |

The LPS column is the whole ADP-heptose → inner-core ladder, in order, recovered gene by gene. **The
model's biomass demands a truncated lipid-A/KDO species and never asks for the completed core**, so every
gene that builds it is invisible to deletion analysis. That is the single largest identified cause in the
no-effect class.

The iron cluster closes a thread this arc has now visited four times: the model gets free Fe²⁺/Fe³⁺, so
the loaded siderophore is never demanded — and here the missing demand is named explicitly (`feenter_c`).

## The two that are named but should NOT be "fixed"

`relA` → ppGpp and `cyaA` → cAMP flip like the others, and that is misleading. **ppGpp and cAMP are
regulatory second messengers, not biomass components.** A biomass equation that demanded them would be
biologically wrong; these genes are missed because FBA has **no regulatory layer at all**, which no
objective term repairs.

So the constructive reading applies to **19 of 34**, not 21 — and the 2 regulatory genes are a
categorically different miss that this method can name but cannot address.

## The 13 with no named demand

`lipB` `trxB` `fruK` `pta` `ppk` `gshB` `gltD` `zntA` `avtA` `spoT` `pfkA` `sthA` `dgkA`

For these the model has an alternative route to **every** metabolite they touch, so no demand can flip
them — 4 of them (`trxB`, `ppk`, `gltD`, `sthA`) touch nothing but currency metabolites at all. A missing
biomass component is **not** the explanation here, and adding one would not help.

## Where the no-effect class now stands

| | genes | cause | addressable by a biomass demand? |
|---|---|---|---|
| cell-envelope + cofactor biosynthesis | 19 | biomass demands a truncated product | **yes — named and verified** |
| regulatory second messengers | 2 | FBA has no regulatory layer | no |
| redundant to their own product | 13 | alternative route exists | no |

## Scope — load-bearing, not boilerplate

**This is diagnostic. It does not add demands to the objective and re-score, and makes no recall claim.**
Doing that would be a **new endpoint on the same data**, and this project retracted a positive result once
for exactly that move and had a second one caught by a guardrail earlier today. If the scoring run
happens, it gets its own pre-registration first — including what the false-positive cost of enriching
biomass is, which this probe deliberately does not measure.

## Honest limits

1. **One condition** (glucose). A route available only on another carbon source would read as sole-route
   here.
2. **The currency-exclusion list is a curated choice.** An over-broad entry silently hides a real answer;
   `tests/test_fba_demand_completion.py` checks against the artifact that no reported metabolite is in
   the list, which is the failure mode that would matter.
3. **"Sole route to *m*" identifies a candidate, it does not settle whether biomass *should* demand *m*.**
   That is a modelling judgment a human makes — and the ppGpp/cAMP pair is the proof that the probe alone
   cannot make it.
4. **The 10 % growth floor is a choice**, reported so it can be varied; without it, routes available only
   to a non-growing cell would count.
5. **Labels, not viability.** The experimental calls are RB-TnSeq fitness (`fit < −2`) in these
   conditions. Nothing here claims an LPS-core mutant is inviable — only that the labels call these genes
   essential and the model cannot.

## Next

The honest next step is **not** "add 19 demands and report the new recall." It is to decide — as a
modelling judgment, with the false-positive cost measured — whether iML1515's biomass *should* demand a
completed LPS core. That is a scoped, pre-registrable question, and for the first time in this arc the
target is named, verified, and small.
