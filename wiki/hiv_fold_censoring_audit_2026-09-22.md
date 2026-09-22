# The one free continuous wet-lab label is right-censored at 100, and the guard against that is vacuous

Run: `uv run python scripts/hiv_fold_censoring_audit.py` (offline, read-only, seconds).
Artifact: `wiki/hiv_fold_censoring_audit_2026-09-22.json`.

**This is a methodological audit of published numbers. It is NOT a decoder change and it does NOT show
any published figure is wrong.** Read the scope rails before quoting it.

## What was measured

The Stanford PhenoSense fold-change tables (`data/raw/hiv/*_DataSet.txt`) are this project's one free,
independent, isolate-level, **continuous** wet-lab label — the substrate that made HIV the first cell to
beat the label wall. Across all four datasets, **24 fold columns, ~40,000 values**:

| | |
|---|---|
| upper bound, shared across all four datasets | **100.0** |
| values anywhere strictly above it | **0** |
| columns reaching it | 22 |
| columns staying below it | 2 (TDF 51.0, BIC 66.0) |

| dataset | worst column | mode-share at 100 | n at ceiling |
|---|---|---|---|
| NRTI | **3TC** | **0.452** | 831 |
| NNRTI | **NVP** | **0.330** | 678 |
| INI | EVG | 0.166 | 125 |
| PI | LPV | 0.144 | 260 |

**3TC and NVP fail this repo's own shipped degeneracy bar** (`assay_degeneracy`, mode-share > 25%) —
the bar imported here rather than restated, so it cannot drift from the rule enforced elsewhere. That
screen exists because CcdB posted the forward/inverse sweep's best number purely because 79.3% of its
variants sat at an assay ceiling. **It had never been pointed at the HIV arm**, which lives in a
different part of the tree.

A hard cap is the parsimonious reading: values immediately below 100 are sparse and irregular (NVP: 87,
89, 90, 91, 91.8, 93 … each n=1–2) and then 678 sit at exactly 100.

## The existing guard is vacuous — its own finding

`hiv_nnrti_mutant_catalog.py:234` and `hiv_epistasis.py:27` both state that censored folds (`<` / `>`
prefixed) are kept at their numeric bound. **Measured: zero operator-prefixed values in any column of
any of the four datasets.** The guard protects against a censoring form that does not occur here, while
the censoring that does occur — a bare `100` with nothing above it — passes straight through.

Same shape as the self-corrected vacuous ResFinder `POINT`-row filter: **a filter that removes nothing
is not a control**, and its presence makes the hazard look handled.

## Why it matters

The NRTI (+0.06..+0.14 on 5/6), PI (+0.056 mean, 8/8) and INSTI (+0.087 mean, 5/5) v0.1 catalogs were
each selected by thresholding a **multivariate-OLS log10-fold coefficient at ≥ log10(1.5)**
(`RESIST_COEF_MIN`), fit on a response right-censored at 100 on drugs where up to **45%** of
observations sit at the ceiling.

Right-censoring attenuates coefficients toward zero **for the strongest effects specifically** — exactly
the majors such a catalog exists to find — so a coefficient threshold fit on a censored response can
systematically exclude true majors. Those three gain figures are published headlines in `CLAUDE.md`.

## What this does NOT claim

- **Not a decoder change.** No v0.1 mutant catalog is deployed: `dna_decode/data/hiv_amr.py` routes
  NRTI/PI/INSTI through **position-based** classes (`NRTI_MAJOR_POSITIONS`, `PI_CLASS`, `INSTI_CLASS`).
- **Does not reopen the declined NNRTI curation.** That verdict was scored at
  `ILLUSTRATIVE_FOLD_CUTOFF = 3.0`, and the ceiling is 100 ≫ 3, so **every censored observation is R
  either way** — the R/S labels, and the sens/spec behind the verdict, are unaffected. Censoring bites
  only in the coefficient-estimation step.
- **Does not re-explain the NNRTI OLS dropping Y181C.** `hiv_amr.py` ships Y181C/I/A, so that is a fact
  about a derived catalog that was never deployed; and the metric that drove the verdict is blind-spot
  recovery, defined over isolates carrying *no* catalogued DRM, which Y181C carriers are not. Censoring
  and the co-occurrence account are **not rivals** — if two co-occurring DRMs both pin fold at the
  ceiling, censoring *removes the residual variance OLS needs to separate them*, which **amplifies** the
  co-occurrence story rather than replacing it.
- **Does not show any published v0.1 gain figure is wrong.** It shows the response they were fit on is
  censored. Whether the numbers *move* under a censoring-aware (Tobit) refit is a separate decidable
  question this audit does not answer.

## Two corrections I made to my own statistic mid-run

**"Every fold column shares one ceiling" was the wrong statistic and reported `None`.** TDF (51) and BIC
(66) never reach the cap in this cohort. Demanding a shared *maximum* made a true finding read as a
refutation. The correct claim is a shared **upper bound that nothing exceeds**; a drug whose folds stay
low is not counter-evidence. Pinned by test so it is not "fixed" back.

**An identifier column was masquerading as a measurement.** `SeqID` is all-integer, near-unique
(2,272 distinct in 2,272 rows), max **789,388** — admitted as a fold column it *became* the dataset
maximum, so the bound was being computed against an accession number. Excluded by a derived rule
(all-integer **and** near-unique), never a hand-listed drug name, which is the documented
`hardcoded_exclusion_list_undercovers` failure.

I also overstated once in the artifact text — "a pile-up in every column that reaches it" — when the
minimum pile-up among reaching columns is **0.0011**. Corrected in place: what is shared is the *bound*;
the heavy pile-ups are where the degeneracy bar fails.

## Honest limits

- The ceiling is **inferred** from a pile-up with an empty tail above, not read off assay documentation.
  Consistent with a reporting cap; this does not prove the assay's dynamic range ends there.
- Attenuation is a property of fitting OLS to a censored response. This audit does **not** measure how
  much any particular coefficient moved.
- The HIV datasets are gitignored, so this cannot run in a clean checkout — it exits 2 rather than
  reporting a vacuously clean result.
- Frozen AMR surface untouched; `hiv_amr.py` unmodified.

## Provenance

Surfaced by the `/innovate` framing sweep (`wiki/innovate_sweep_2026-09-21.md`, survivor `QS-HIVCEIL`)
under the winning framing — *replace the blocked question with one the data on hand can decide*. The
candidate's originally stated target (Y181C) was **killed** by the adversarial pass and replaced; what
survived is the audit above.
