# The gentamicin `rmt` gap: the fix is motivated, in-distribution-safe as far as measurable, and its risk is UNTESTED

The first prospective accrual located a real catalog gap — E. coli × gentamicin sens 0.429, with 24 of 28
false negatives carrying a 16S rRNA methyltransferase and 0 carrying `aac(3)`. This is the check that
gates any revision: **would adding `rmt` break the cells that currently work?**

Runs entirely on cached AMRFinder output. No Docker, no network. **The frozen surface is not touched** —
the candidate is applied scorer-locally, mirroring the `experimental_drug_rules.py` overlay pattern, and a
test asserts `amr_rules.py` and `calibrated_amr_rules.json` contain no trace of it.

## Result, over 150 labelled isolates with cached determinant calls

| | frozen | candidate |
|---|---:|---:|
| accuracy | 0.927 | **0.933** |
| sensitivity | 0.907 | **0.920** |
| specificity | 0.947 | 0.947 |
| tp / fp / tn / fn | 68 / 4 / 71 / 7 | 69 / 4 / 71 / **6** |

**Control first:** the local re-implementation of the frozen rule agrees with the real `call_resistance`
on **150/150** isolates. Without that, none of the above means anything.

**One call changes.** `GCA_022316245.1`, label **R**, carries `rmtB1` — a genuine false negative rescued.
It is the *only* methyltransferase-carrying false negative in the data, so the candidate rescues 1 of 1.

## The specificity result is VACUOUS, and that is the important part

"Specificity unchanged at 0.947" reads like the candidate was tested for over-calling. **It was not.**

Only an **S-labelled carrier** can become a new false positive. There are **zero** of them: all 4
methyltransferase carriers in the labelled gentamicin data are R. The candidate therefore *cannot* produce
a false positive here — the check is true and carries no information.

**Its over-calling risk is UNTESTED, not zero.** The script now prints exactly that instead of the
reassuring-sounding verdict it printed first, and a test pins the wording.

## Two things I got wrong, both caught by checking

1. **I expected a second blind spot that does not exist.** Compound subclasses that name gentamicin
   (`GENTAMICIN/KANAMYCIN/TOBRAMYCIN`, `APRAMYCIN/GENTAMICIN/TOBRAMYCIN`, 48 hits) looked like they would
   be missed by `subclass_any={"GENTAMICIN"}`. The rule matches by **token**, so they are already counted —
   verified empirically on a real genome before the claim went anywhere. Pinned by test.
2. **My first run scored 63 of ~148 isolates** and reported the rest as "missing AMRFinder", which reads
   as a data gap. It was a lookup bug: AMRFinder output is per-**genome**, not per-drug, so cohorts share
   runs (`klebsiella_gentamicin`'s live under `klebsiella_cipro/amrfinder_runs/`). A global accession
   index took the control from 63/63 to **150/150**. A "missing data" number should always be suspected of
   being a path bug first.

## Status and what would settle it

The candidate is **well-motivated** (mechanism exact: 16S methyltransferases confer high-level
4,6-deoxystreptamine resistance including gentamicin; AMRFinder files them under generic `AMINOGLYCOSIDE`,
which carries no `GENTAMICIN` token) and **does no in-distribution harm as far as the data can show**.

It is **not deployed**, and should not be: changing the frozen surface invalidates both the prospective
lock and the reproducibility freeze. A v2 lock is a **user authority call**, not an executor one.

What would actually settle it, in order of value:

1. **S-labelled `rmt` carriers.** Without them the over-calling risk stays unmeasured. This is the binding
   constraint, and it is a *label* problem — the familiar wall.
2. **The prospective cohort's per-isolate determinant calls** (not on disk; the committed artifact carries
   only a confusion matrix), which would measure the rescue rather than assume it from the memo.
