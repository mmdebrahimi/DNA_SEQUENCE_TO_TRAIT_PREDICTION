# The Klebsiella × meropenem specificity question was already answered — and three later attempts were weaker, not rival (2026-10-04)

**Headline: there is no specificity conflict for this cell.** The provenance-disjoint reading
(**spec 0.900, N=60**) is the only one of four that clears every bar this project applies, and it has
done so since 2026-06-10. The three readings treated as rivals to it are each *weaker measurements of
the same thing*, not disagreeing evidence.

**This memo exists because I got the framing wrong, not the numbers.** A plan written this session set
out to "produce a fourth number to adjudicate 0.733 / 0.900 / 0.000" and built a 129-isolate external
cohort to do it. Every individual figure in that framing was correct. What was never checked is whether
one of the three was already decisive — and it was. This is
[[feedback_verify_the_framing_not_just_the_numbers]] applied to my own plan.

## The four readings, each re-read from its artifact

| route | N | sens | spec | status against this project's own bars |
|---|---|---|---|---|
| frozen rule's `validated` string (calibration) | 30 | 1.000 | 0.733 | **in-sample.** Not a validation number. |
| **provenance-disjoint** (`provenance_disjoint_validation_klebsiella_merop_2026-06-10.json`) | **60** | **0.467** | **0.900** | **CLEARS EVERY BAR** — see below |
| prospective-lock (`prospective_lock_validation_Klebsiella_meropenem_2026-10-03.json`) | 16 (7R/9S) | 1.000 | 0.000 | **`UNDERPOWERED`** by its own gate (`need >= 10/class; have R=7 S=9`); one submission batch |
| source-diverse arm (`source_diverse_validation_klebsiella_meropenem.json`) | 44 cand. | — | — | **REFUSED** — and on a *different pool* (see below) |
| AR Bank external (in flight, this session) | 129 lab. | — | — | **G5 attrition** — scorable S = 5 |

Confusion matrix for the provdisjoint cell, verbatim: `tp=14, fp=3, tn=27, fn=16, abstain=0,
acc=0.683, sens=0.467, spec=0.900`.

## Why the provdisjoint 0.900 clears the bar — the check that inverted the conclusion

From `wiki/provdisjoint_source_concentration.json`:

```
organism: klebsiella   drug: meropenem   n_cohort: 60   spec: 0.9   sens: 0.467
bioproject: distinct 15, largest PRJNA504784 = 22, largest_share 0.367
```

The project's source-diversity bar (`scripts/source_diverse_validate.py`) is **>=5 distinct sources AND
largest share <=0.60**. This cell is at **15 and 0.367** — it passes with margin. It is also
provenance-disjoint by the fail-closed accession manifest (799 accessions excluded across 31 prior
cohorts, `manifest_complete=True`, `manifest_degraded=False`) and has the ecosystem-domination
exclusion list applied.

So the specificity question has had a source-diverse, provenance-disjoint, powered answer for four
months.

## The source-diverse refusal was about a DIFFERENT cohort — this is the trap worth keeping

The source-diverse arm's refusal (`only 1 BioProject(s), bar is 5`) reads naturally as *"this cell is
too concentrated to measure"*, and I read it that way. It is not that. That arm assembles a **fresh
never-before-scored candidate pool** and refuses when the pool is concentrated:

```
n_candidates: 44   distinct: 1   largest_share: 1.0   dominant: PRJNA717739
```

`PRJNA717739` is not the provdisjoint cohort's dominant source (`PRJNA504784`, 22 of 60). The two pools
are disjoint in composition. **The arm refused its own feedstock; it never evaluated the 60-isolate
cohort.** A refusal is a statement about the thing refused — check which thing that was before letting
it overwrite an existing result.

## The prospective 0.000 is self-labelled underpowered, and the direction matters

9 susceptible isolates from one submission batch, and the artifact's own `powering` block says
`UNDERPOWERED`. It was already diagnosed: all 9 false positives carry a genuine carbapenemase
(blaOXA-48 ×5, blaNDM-1 ×2, blaOXA-181, blaKPC-3), i.e. **determinant-present vs measured-MIC**, which is
[[feedback_determinant_presence_is_not_phenotype]] — not a detection failure. A spec of 0.000 on 9
isolates from one batch cannot displace 0.900 on 60 across 15 BioProjects, and was never framed as
able to; the error was in treating it as a third rival reading rather than as a diagnosed outlier.

## What this does to the in-flight AR Bank run (R1) — it flips its purpose

R1 was launched to adjudicate specificity. It cannot (scorable S=5, per the `/brainstorm` finding and
the G5 gate) and it **does not need to** (already answered). But it should still finish, because it can
measure the axis that is genuinely weak:

**Sensitivity is this cell's failure mode, and the frozen rule says so itself.** Its `validated` string:

> `Excludes ESBL/AmpC; blind to porin-loss-mediated R (expected FN mode).`

Provdisjoint confirms the prediction: **sens 0.467 with fn=16** — sixteen resistant isolates missed, the
largest error term in the cell. The AR Bank cohort's resistant class is **43R across 5 panels, largest
share 0.417**, which clears the diversity bar numerically.

So R1's honest deliverable is an **independent external check on the weak sensitivity axis and the
declared porin-loss blind spot**, with specificity explicitly withheld. That is a more useful number
than the one it was launched to get.

**Caveat on the independence unit, carried forward:** AR Bank cohorts carry no BioProject field, so
"panel" is a proxy and the existing AR Bank source-concentration artifact flags
`grouping_is_a_proxy`. If panels are one collection process, "5 panels" overstates independence. Report
it as a panel-level count, never as a BioProject-equivalent.

## Reusable

1. **Before building a tie-breaker, check whether one candidate already wins.** Four readings that
   disagree are not automatically four comparable measurements. Rank them by the bars the project
   already enforces (powering, source diversity, leakage posture) *before* spending labour on a fifth.
2. **A refusal names a specific object.** `status: source_concentrated` on a fresh candidate pool says
   nothing about a different, older, diverse cohort for the same cell. Check the pool identity.
3. **When a cell's own rule declares a blind spot, that is where measurement is worth buying.** The
   meropenem rule declares porin-loss blindness and the data shows fn=16 — so sensitivity, not
   specificity, was always the open question for this cell.

## Scope limits

- This memo re-reads existing artifacts; **it measures nothing new** and changes no number.
- The frozen surface is untouched; the active 2026-08-31 prospective lock is unaffected.
- `spec 0.900` carries the provdisjoint tier's own standing caveats: isolate-level (not
  clonality-corrected here), and `sra_center` is `NULL` for all 60 so BioProject is the only resolvable
  source axis.
- It does **not** claim the prospective 0.000 is wrong about anything — a determinant-present /
  MIC-susceptible gap is real and is what the L2 doubt layer is for. It claims only that 9 isolates from
  one batch cannot set the cell's specificity.
