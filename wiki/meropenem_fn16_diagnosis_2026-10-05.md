# The 16 meropenem false negatives are all a different species (2026-10-05)

## The question

`Klebsiella × meropenem` is a SCORED cell on the standing trust surface: **acc 0.683 / sens 0.467 /
spec 0.900, tp=14 fp=3 tn=27 fn=16, N=60**. Its 16 false negatives had one recorded explanation — that the
frozen rule is *"blind to porin-loss-mediated R (expected FN mode)"*, stated verbatim in
`amr_rules.DRUG_RULE['meropenem']['validated']`.

**That explanation is false**, measured 2026-10-04: AMRFinder files `ompK35`/`ompK36` truncations under
`Subclass=CARBAPENEM`, and the rule is `threshold=1, subclass_any={'CARBAPENEM'}` — so it *counts* porin
loss (10/10 porin-only isolates are called R). The cell's dominant error mode therefore had **no**
explanation. This is what it is.

## The answer

**All 16 false negatives are *Klebsiella aerogenes*. Not one is *K. pneumoniae*.**

The cohort is **38% a different species**: 37 *K. pneumoniae* + **23 *K. aerogenes*** out of 60, scored
throughout with AMRFinder `-O Klebsiella_pneumoniae` under a rule whose own `validated` string records it
as calibrated on *"Klebsiella N=30"*.

| species | n | tp | fp | tn | **fn** |
|---|---|---|---|---|---|
| *Klebsiella pneumoniae* | 37 | 12 | 3 | 22 | **0** |
| *Klebsiella aerogenes* | 23 | 2 | 0 | 5 | **16** |
| **total** | **60** | **14** | **3** | **27** | **16** |

**On the species the rule was validated on there are zero misses — sensitivity 12/12 = 1.000.** Every
miss in the cell comes from the 23 off-species isolates.

So the published **sens 0.467 is substantially a measurement of cohort composition**, not of the rule.

## Why this is believable

1. **The committed confusion matrix reproduces EXACTLY** — tp=14 fp=3 tn=27 fn=16, 0 abstain — by
   re-running the deployed `call_resistance(main_tsv, 'meropenem', organism='Klebsiella')` over the 60
   cached AMRFinder runs. Nothing is attributed until that gate passes
   (`compute_lineage_metrics.reconcile_raw_metrics`, which raises rather than returns on disagreement).
2. **The cross-tab reconciles with the matrix it decomposes** — the per-species cells sum to
   14/3/27/16, enforced by `cohort_species.compose(..., confusion=...)` which raises `CrosstabMismatch`
   otherwise.
3. **The determinant fingerprint agrees independently of the species labels.** All 16 FN carry only
   `ampC-Kaer` (13) / `ampC_Kaer-1` (3) at `Subclass=CEPHALOSPORIN` — `Kaer` being *K. aerogenes*' own
   intrinsic AmpC — and **zero `ompK*` rows**. Two different routes to the same conclusion.
4. **The mechanism is coherent.** *K. aerogenes* reaches carbapenem resistance largely through AmpC
   hyperproduction plus porin loss, not through the acquired carbapenemases (`blaKPC`/`NDM`/`OXA-48`) the
   rule counts. A rule targeting acquired carbapenemases is *expected* to miss it.

**The porin story fails twice over:** the rule can see porins, and these isolates have no porin call for
it to see.

## It is not systemic — and the contrast is the useful part

All 10 provenance-disjoint cohorts were audited (`wiki/provdisjoint_species_audit_2026-10-05.{md,json}`):

| cohort | scored as | composition | off-species | FN off-species |
|---|---|---|---|---|
| campylobacter cipro | `Campylobacter` | MATCHES_SCORED_ORGANISM | 0 | 0/0 |
| e.coli ceftriaxone / cipro / gentamicin / tetracycline | `Escherichia` | MATCHES_SCORED_ORGANISM | 0 each | 0/1, 0/2, 0/3, 0/2 |
| klebsiella ceftriaxone | `Klebsiella_pneumoniae` | MIXED_SPECIES | 2 | *withheld* |
| klebsiella ciprofloxacin | `Klebsiella_pneumoniae` | MIXED_SPECIES | 2 | 0/1 |
| klebsiella gentamicin | `Klebsiella_pneumoniae` | MIXED_SPECIES | 1 | *withheld* |
| **klebsiella meropenem** | `Klebsiella_pneumoniae` | **MIXED_SPECIES** | **23** | **16/16** |
| klebsiella tetracycline | `Klebsiella_pneumoniae` | MIXED_SPECIES | 2 | *withheld* |

**Meropenem is a ~10x outlier.** Every Klebsiella cohort carries *some* *K. aerogenes* (1–2 isolates,
2–3%), which is ordinary cohort noise; meropenem carries 23 (38%). And meropenem is the **only** cohort
where off-species isolates account for the misses — Klebsiella cipro has 2 off-species and 0 of its 1 FN
is one of them.

**Three Klebsiella cohorts have their cross-tab WITHHELD, not computed and hidden.** Their cached
AMRFinder runs are incomplete (gentamicin 3/60, tetracycline 33/60, ceftriaxone 54/60), and attributing
outcomes on a partial run set would compute a rate over a silently-shrunken denominator. Composition is
still reported for them because it needs only the assemblies, which are complete everywhere.

## Two findings that correct the plan's own expectations

- **The E. coli cohorts contain no *Shigella* at all** — 60/60 *Escherichia coli* in all four. The
  directory name `escherichia_coli_shigella` is the NCBI-PD taxgroup label, not a statement about
  content. The plan listed this as an open question; it is answered negatively.
- **`Campylobacter` is two species by design** — 31 *C. jejuni* + 9 *C. coli*, scored with the
  **genus-level** `-O Campylobacter`. It legitimately matches what it was scored as while not being
  monomorphic, which is why the verdict is named `MATCHES_SCORED_ORGANISM`: a first draft called it
  `SINGLE_SPECIES_AS_EXPECTED` and this cohort proved that name false.

## What this does NOT do

- **It does not re-score anything.** The 10 `provenance_disjoint_validation_*.json` artifacts are named
  frozen units of the reproducibility freeze and are byte-unchanged. The *K. pneumoniae*-only matrix
  (tp=12 fp=3 tn=22 fn=0 → sens **1.000**, spec **0.880**, acc **0.919**, N=37) is stated here as a
  measurement, **not** published as the cell's number. Replacing the cell's metrics is a **user authority
  call**.
- **It does not change the meropenem rule.** The lone-porin over-call (all 3 FPs on the separate AR Bank
  cohort are porin-only) is a different, still-open authority fork.
- **It is not an independent species call.** Species come from each assembly's own GenBank `ORGANISM`
  line — the **submitter's** assertion, the same evidence class the project's G1 circular-label gate
  covers. It is strong evidence about what the cohort *contains*; it is not wet-lab identification.
  `dna-identify` would be a second signal of the same broad kind (sequence comparison), and
  *K. aerogenes* is outside its 14 supported organisms, so it would **abstain** — that cross-check was
  **not run** here (it needs Docker for Mash).

## Honest limits

- `-O Klebsiella_pneumoniae` on a *K. aerogenes* genome also drives species-specific point-mutation
  screening, so determinant calls for those 23 isolates may be wrong in ways beyond this one rule. Not
  investigated; AMRFinder was not re-run under a corrected flag.
- Whether 23 *K. aerogenes* entered this cohort by a selection artifact or because the source labelled
  them as *Klebsiella* at genus level is **not** established here.
- The two *K. aerogenes* true positives (2 tp) are unexplained — they were called R, so they carry some
  CARBAPENEM-subclass determinant. Not chased.
- `MAX_UNRESOLVED_FRACTION = 0.10` is an asserted bar, not a derived one. It never bound: every one of
  the 580 assemblies across all 10 cohorts resolved an `ORGANISM` line.

## Reproduce

```bash
uv run python scripts/provdisjoint_species_audit.py                # all 10, read-only, offline, exit 0
uv run python scripts/provdisjoint_species_audit.py --cohort klebsiella_provdisjoint_meropenem
```
