# Finding the model's missing parts — gap-finding + repair (2026-08-07)

Track C of the design epoch (`wiki/design_epoch_plan_2026-08-07.md`), and the honest form of the user's
"predict the missing parts". Not nucleotide infilling — **the missing biochemistry**, which is where a
genome-scale model is actually wrong and where being wrong silently corrupts every prediction built on it.

Ships as `dna_decode/fba/gapfill.py` + `dna-decode fba --dead-ends` / `--gapfill-target`.

## Two capabilities, deliberately separated by evidential weight

| | what it is | needs |
|---|---|---|
| **`--dead-ends`** | a **structural fact** about the model: metabolites produced but never consumed (or the reverse) cannot carry steady-state flux | nothing — no labels, no donor, no network |
| **`--gapfill-target`** | a **hypothesis**: which donor reactions would restore a trait the model wrongly predicts as absent | a donor model + independent evidence the organism has the trait |

Keeping these apart is the point. The first is checkable arithmetic; the second is a claim about biology.

## The worked case — a measured false negative, found and repaired

**The gap is real, verified at source.** iML1515 predicts **no growth on sucrose**, but *E. coli* BW25113
has a **Sucrose** carbon-source experiment among the 28 in the Wetmore/Keio RB-TnSeq set, and that assay
only runs sources the organism grows on. (Fetched and confirmed directly, not taken from the summary.)

**The diagnostic localises it unprompted.** Over iML1515: **138 dead-end metabolites, 42 of them fed by a
transport reaction** — the structurally suspicious class, where the model imports something nothing then
consumes. `suc6p_c` is in that list, produced by `SUCptspp` and consumed by nothing:

```
suc6p_c   no_consumer   from SUCptspp
```

So the model carries a sucrose PTS transporter that leads nowhere. **No sucrose knowledge was used to find
this** — it falls out of the stoichiometry.

**The repair is a single reaction, and it is measured.** Gap-filling against *Salmonella* iYS1720 (1125
candidate reactions the model lacks) proposes exactly one:

```
FFSD      h2o_c + suc6p_c --> fru_c + g6p_c
```

| | before | after |
|---|---|---|
| growth on sucrose | **0.000 /h** | **1.7798 /h** |
| glucose / xylose / cellobiose | — | **unchanged** |

The rate is internally consistent: a disaccharide yielding two hexoses should give roughly twice the
glucose rate (0.877 → 1.78). And the specificity check matters as much as the target — a "repair" that also
lifts growth on unrelated carbon sources has made the model permissive rather than correct, and looks
identical to success if you only watch the target.

## Honesty rails (in the code, the CLI output, and the tests)

- **A gap is not automatically a defect.** Organisms genuinely lack capabilities; "repairing" a correct
  model fabricates biology. Every proposal ships stamped `HYPOTHESIS — a gap may be correct biology; verify
  the organism truly has the trait`. This case earns its repair only because the trait is independently
  measured. Worth stating plainly: *E. coli* K-12 sucrose catabolism is strain- and stock-dependent, so the
  evidence here is the RB-TnSeq experiment's existence, not a universal claim about K-12.
- **Reversible reactions count as both producer and consumer.** Scoring them one-directionally invents dead
  ends that are not dead — and a diagnostic is worthless if its hits aren't real. Pinned by a test.
- **`demand_reactions` and `exchange_reactions` are both OFF** in the search. Letting the solver invent a
  demand or an exchange lets it satisfy the objective by inventing a sink instead of finding the missing
  biochemistry. That is not a repair.
- **Boundary reactions are excluded** from the dead-end scan: exchanges exist precisely to be one-sided, so
  counting them buries the real hits under hundreds of uninteresting ones.

## Limits

- The donor is one related organism's reconstruction, so the reachable repairs are bounded by what the donor
  knows. A true universal reaction database (BiGG universal / ModelSEED) would widen it.
- The remaining 41 transport-fed dead ends are **unexamined** — each is a candidate gap, not a known defect,
  and most may be curation choices rather than errors.
- Nothing here is validated at the bench. A repaired model is a better hypothesis, not a measured organism.

## Run

```bash
dna-decode fba --dead-ends                     # structural: where are the missing parts?
dna-decode fba --gapfill-target EX_sucr_e      # hypothesis + measured verification
dna-decode fba --gapfill-target EX_sucr_e --donor iYS1720 --json
```

Tests: `tests/test_fba_gapfill.py` (8 pure + 1 real-model gate asserting the diagnostic surfaces `suc6p_c`
unprompted, the proposal contains `FFSD`, and the repair does not leak onto other carbon sources).
