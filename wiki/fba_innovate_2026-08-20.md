# /innovate on the 4-lever wall — 0 engine-survived, 2 hand-verified, 6 deferred

**Date:** 2026-08-20 · **Ledger:** `wiki/fba_innovate_ledger_2026-08-20.json`
**Trigger:** (b) ≥2 approaches failed the SAME way — gap-fill, threshold retune, pFBA restriction and
E-Flux all died against "the deletion changed nothing."

## Engine verdict (do not skip this line)

    0 survived · 0 killed · 8 unfalsified

**No candidate is presented as `survived`.** The two strongest have real, executable falsifiers that I
ran by hand and that PASS — but the gated runner could not execute them, so calling them survivors would
be exactly the theater this pipeline exists to prevent.

Why the engine could not run them, and why I did not route around it:
- `uv run pytest …` is not on the runner's known-safe allowlist → fail-closed → `unfalsified`.
- The allowlisted `python -m pytest` **would have produced a FALSE survival**: the ambient python has no
  `cobra`, so `pytest.importorskip("cobra")` skips the tests and pytest exits 0. Verified directly:
  `collected 0 items / 1 skipped`. An exit-0 from a skipped test is indistinguishable from a pass.

So the honest state is: **falsifier defined, committed, hand-executed, engine-unverified.**

## What was measured (this is the real output of the run)

### M1 — Zero-flux cells are unreachable by capacity-tightening *(theorem, then tested)*

> **SCOPE CORRECTED 2026-08-21.** This section originally said *"No constraint-based method can move such
> a cell."* That is too broad and is **refuted by a later measurement in this same repo**:
> `wiki/fba_structural_blindspot_2026-08-21.md` found that of 58 genes pinned at zero flux in glucose,
> **57 carry flux fine once other exchanges are opened** — and opening an exchange IS a constraint change.
> The corrected claim is below; the numbers in this section are unaffected.

If some optimal solution carries zero flux through every reaction of gene *g*, deleting *g* leaves that
solution feasible and optimal, so the growth ratio is **exactly 1.0**. **Holding the objective and the
medium FIXED**, no intervention that merely tightens or removes capacity on reactions carrying no flux
can move such a cell — that class of intervention only reshapes flux among reactions that already carry
it. Changing the medium, changing the objective, or adding a demand term *can* remove the zero-flux
optimum, and are not covered by this claim.

Measured on the E. coli 25-carbon panel, each cell **in its own condition** (1,832 true-essential cells):

| class | count | share |
|---|---|---|
| A — reaction carries no flux | 459 | 25.1 % |
| B — carries flux + isozyme `or` | 512 | 27.9 % |
| C — carries flux, no `or` | 861 | 47.0 % |

Class A is **necessarily missed** → 459 of the panel's 1,083 misses (**42.4 %**) are unreachable by any
of the four levers. That is a mechanism for their identical failure, not four coincidences.

**Class A is a LOWER BOUND, not an estimate** — the proof needs only *some* zero-flux optimum, so a cell
in B or C may also admit one. FVA-at-optimum would give the true set; not done.

### M2 — iML1515's internal bounds are NOT binding *(this explains the E-Flux exact null)*

Scaling every internal, GPR-carrying reaction bound:

| scale | growth | |
|---|---|---|
| 0.5 / 0.2 / 0.1 | 0.876997 | **bit-identical — inert** |
| 0.05 | 0.714826 | bites |
| 0.002 and below | infeasible | |

**Positive control:** halving *carbon uptake* does bite. So this is "internal bounds are slack, uptake is
binding" — not "the model is numb to everything."

E-Flux scales internal bounds and exempts exchanges **by design**. It was therefore constraining capacity
that was never limiting. That is why 0 of 1,441 calls moved.

**Found by a vacuity guard, not by design.** My first version of the zero-flux test emulated E-Flux at
0.5× and its own "did the perturbation actually bite?" assertion failed with bit-identical growth. The
guard's failure *is* the finding.

### M3 — A concrete instance of class A: the model gets free iron

`EX_fe2_e`, `EX_fe3_e`, `EX_zn2_e` all have bounds `(-1000, 1000)` and sit in the default medium. All
**12/12** named Fe/Zn acquisition genes (`fepA fes fepB fepC fepD fepG tonB exbB exbD znuA znuB znuC`)
carry **zero flux** at the optimum. The whole siderophore system is idle because iron is a free tap —
while in a real cell it is essential. *(Surfaced by the G5 operator; verified here independently.)*

## The 8 candidates

**Unfalsified — falsifier hand-executed and passing, engine-unverified (see above):**

| id | claim | evidence |
|---|---|---|
| S1 | zero-flux cells unreachable by capacity-tightening (fixed objective + medium) | `tests/test_fba_zero_flux_theorem.py` (2 pass) |
| S2 | capacity slack explains the E-Flux null | `tests/test_fba_capacity_slack.py` (3 pass, incl. positive control) |

**Unfalsified — deferred, falsifier defined but blocked:** `D: is disconnected`, so `feba.db` (7.4 GB,
uncommitted) and PRECISE-1K are unreachable. Every label-side test is unrunnable.

| id | operator | claim |
|---|---|---|
| C1 | G8 | insert expression at the **GPR layer**, not the flux-bound layer (an `A or B` rule is only redundant if B is transcribed) |
| C2 | G8 | the stress NO-GO is over-broad for antibiotics whose target is **metabolic** — the dissociation vs non-metabolic drugs is the claim |
| C3 | G8 | cofitness propagation from the FBA-essential seed, **recomputed leave-condition-out** (the precomputed `Cofit` table is circular here) |
| R1 | G1 | stop predicting genes; **audit the model** — flat cells implicate a repeating reaction set = a curation-target list |
| R2 | G1 | unit becomes the **OR-rule × condition** → a *conditional* GPR annotation |
| R3 | G1 | flatness is predictable from the **model alone** → a declared coverage map, the abstention architecture this repo already ships |

## Corrections made during the run

1. **My own class-C framing was wrong.** I described class C as "compensated by rerouting." The
   denominator is *all* true-essential cells, and class A is 100 % missed, so **every one of the model's
   true positives lives in B∪C** — class C mixes hits with misses and is not a failure class. Class-C
   misses lie in [112, 624] of 1,083; the split is unmeasured. *(Caught by the G9 operator.)* Corrected
   in the script docstring.
2. **The glucose-only shortcut overstated class A 2×** (50.2 % → 25.1 % per-condition). Caught by re-running
   properly before publishing.
3. **A sub-agent premise was partly overstated.** G8 reported D-cycloserine's negative tail as "entirely
   the target pathway"; the strongest negatives are actually `fur`/`dps`/`ycfM` (oxidative stress), with
   `mrcB` 4th. Its fosfomycin premise (`glpT` +11.31 at rank 1) is **exactly right**.

## Next

The two hand-verified findings are model-only and need no D:. The highest-value follow-up is also
model-only: **FVA-at-optimum over the class-A reactions**, which would (a) make class A basis-independent
instead of a lower bound and (b) split "blocked/dead-end in the reconstruction" from "merely unused" —
two different fixes. The six deferred candidates unblock the moment D: is reconnected.
