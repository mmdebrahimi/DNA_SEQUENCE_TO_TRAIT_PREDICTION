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

---

## Follow-up the same day: an independent cross-check, and the 2 unexplained true positives

Both items this memo left open are now closed. Docker was unavailable when the memo was written and is
up now, which is the only reason these were deferred.

### 1. `dna-identify` corroborates, by a genuinely unrelated route

The memo's species evidence (GenBank `ORGANISM` + the `ampC-Kaer` fingerprint) shares one annotation
provenance. `dna-identify` is different: Mash sketch k-mer distance against a 14-organism closed-set
reference that **does not contain *K. aerogenes***, whose `thresholds.py` records K. aerogenes landing at
**0.10849–0.13386**, above `MAX_DISTANCE = 0.1036`.

Prediction registered before running: the FN should ABSTAIN; the *K. pneumoniae* true positives should be
called. **A false negative returning `klebsiella_pneumoniae` would have killed the claim.**

| group | n | called `klebsiella_pneumoniae` | abstained | distance |
|---|---|---|---|---|
| false negatives | 16 | **0** | **16** | 0.1129–0.1339 |
| true positives | 14 | 12 | 2 | called: 0.0011–0.0090 |

**Partition agreement with the GenBank labels: 30/30, zero disagreements.** Every isolate the GenBank
line calls non-*pneumoniae* is exactly an isolate the router abstains on, and the abstaining distances sit
inside the band `thresholds.py` records for K. aerogenes — an order of magnitude further than the called
controls (0.11–0.13 vs 0.001–0.009).

**HONEST SCOPE, and it matters:** `dna-identify` is a CLOSED-SET router over 14 organisms. It can only
say `klebsiella_pneumoniae` or ABSTAIN, so what it independently establishes is **not-*pneumoniae***, not
*aerogenes* specifically. The distances falling in the recorded aerogenes band are suggestive, not a
species call. It is also a *sequence-comparison* method, so it is a second signal of the same broad kind
as a submitter's sequence-derived label — independent of the annotation, not of the evidence class.

### 2. The 2 *K. aerogenes* true positives are explained — and they are NOT one class

| isolate | sole counted carbapenem determinant | reading |
|---|---|---|
| `GCA_003951185.1` | `ompK35_G41TfsTer32` (porin truncation, `Subclass=CARBAPENEM`) | a **lone-porin call that is CORRECT** |
| `GCA_003951605.1` | `blaNMC-A` (class-A carbapenemase) | a legitimate acquired-carbapenemase call |

Both also carry `ampC-Kaer`/`ampC_Kaer-1` at `Subclass=CEPHALOSPORIN`, which the rule does not count —
the same fingerprint as all 16 false negatives. So within the 23 *K. aerogenes* isolates the rule fires
only where a CARBAPENEM-subclass determinant happens to be called, and the AmpC that actually drives
carbapenem resistance in this species is invisible to it.

**This bears on the open lone-porin authority fork, and it points the other way.** CLAUDE.md records all
3 false positives on the separate AR Bank cohort as porin-only, which argues for requiring a carbapenemase
rather than counting a lone porin truncation. `GCA_003951185.1` is the counter-case: a lone-porin call
that is **right**. Dropping lone-porin counting would convert it from a true positive into a false
negative. **n=1 against n=3 — this does not settle the question and is not offered as settling it**; it
means the fix has a measurable cost on the same arm, so the decision should be made against both numbers
rather than the FP count alone.

**Reproduce:** `scratchpad/identify_crosscheck.py` (Docker Mash; registers its prediction before reading
output and reports `INDETERMINATE_CONTROLS_FAILED` if the controls do not resolve, so a uniformly
abstaining router cannot pass by accident).
