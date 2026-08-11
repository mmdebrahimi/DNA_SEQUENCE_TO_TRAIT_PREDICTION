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

## Cross-organism census — and a NAMED HYPOTHESIS (not a finding)

The diagnostic is cheap enough to run over every model the cell supports:

| model | organism | reactions | metabolites | dead ends | % of metabolites | transport-fed |
|---|---|---|---|---|---|---|
| iML1515 | *E. coli* K-12 | 2712 | 1877 | 138 | **7.4%** | 42 |
| iYS1720 | *Salmonella* (pan) | 3357 | 2436 | 312 | 12.8% | 186 |
| iYS854 | *S. aureus* USA300 | 1455 | 1335 | 196 | 14.7% | 29 |
| iMM904 | *S. cerevisiae* | 1577 | 1226 | 225 | **18.4%** | 71 |

**The hypothesis (labelled `unfalsified` — no kill-test run):** the weak transfer of FBA gene-essentiality
away from *E. coli* (`wiki/fba_per_organism_essentiality_2026-08-03.md`: E. coli MCC **0.652** strong →
yeast MCC **0.252** weak) may be driven by **model incompleteness rather than by FBA or by yeast biology**.
The two organisms with measured essentiality sit at the extremes of this table — E. coli's model is the most
complete at 7.4%, yeast's the least at 18.4%, ~2.5× the dead-end rate.

**Why this is a hypothesis and not a result:** it rests on **n = 2** organisms that have both numbers, and
the obvious confound is uncontrolled — iML1515 is among the most intensively curated reconstructions in
existence, so curation effort plausibly drives both its low dead-end rate *and* its strong essentiality
score without one causing the other. Two points are a coincidence, not a trend.

**The decisive test:** gap-fill iMM904's dead ends against a donor, re-run the SGD inviable-null validation,
and see whether MCC moves. If repairing the model lifts essentiality accuracy, incompleteness was the
mechanism; if it doesn't, the weak transfer is about something else.

> **Correction (same day, before this shipped):** an earlier draft of this section called that test
> "runnable". **It is not, with what is on hand** — it needs a *yeast-compatible* donor, and every donor
> here is bacterial (cross-kingdom gap-filling would mostly fail on compartment and identifier mismatch).
> BiGG's only other *S. cerevisiae* model is **iND750**, a smaller **predecessor** of iMM904, so it would
> contribute almost nothing. A real test needs an external reconstruction (Yeast8 / Yeast-GEM, GitHub,
> different identifier namespace) — a **data-acquisition step, not a config change**. Checked rather than
> assumed, because the whole point of this module is not to assert what hasn't been verified.

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
