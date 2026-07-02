# Off-pathogen track — BANKED (2026-07-01)

> **UPDATE 2026-07-01 (post-banking, user-requested deep research):** the horse-coat cell's data wall was
> CLEARED — a deep-research pass found the Sarcidano 2022 open-access study (PMC9558981, Table 3), yielding a
> functional joint genotype × visually-observed colour dataset. Horse rule now **SCORED: concordance 1.000
> (62/62 base-colour horses), non-circular** (`data/horse/sarcidano_2022.tsv`,
> `wiki/horse_coat_validation_2026-07-01.json`). This CLOSES the horse cell with a real number rather than a
> data-wall; the banking of the broader track still stands (it remains deployed-rule integration on a single
> isolated breed — confirming, not new signal). The horse row below is updated accordingly.

Closeout for the off-pathogen expansion of the deterministic gene→trait decoder. This calls the plateau: the
track is banked at its honest terminal. Further autonomous increments here are motion, not signal.

## What the track established (durable)

The deterministic curated-rule + free-independent-label architecture **generalizes beyond pathogens** — shown
across two trait classes and two cohorts, all with honest tiers + frozen AMR surface intact throughout:

| Cell | Rule | Validation | Tier |
|---|---|---|---|
| Eye colour v0 (rs12913832) | HERC2 single-SNP, strand-agnostic | OpenSNP acc 0.993 (decisive homozygotes; abstains on hets) | PILOT/DEMO (self-report) |
| Eye colour v0.1 (IrisPlex 6-SNP) | published Walsh/HIrisPlex model form | OpenSNP deployed-0.7 acc 0.997 / N=346; argmax 0.985 / N=391 (+49% coverage) | PILOT/DEMO (self-report) |
| Ancestry confound | 1000G structural (EUR/EAS 318×) + within-EUR re-score 45/45 | CONFOUND_REDUCED_NOT_RESOLVED | — |
| Eye colour PGP replication | same callers, 2nd cohort | v0 11/11 + v0.1 15/15, zero errors (underpowered, ~5% linkage yield) | PILOT/DEMO |
| ABO blood type (rs8176719 O-status) | c.261delG deletion → O | non-O spec 1.0, n=6 pilot | PILOT/DEMO (serological) |
| Horse coat colour (MC1R-E × ASIP-A) | Rieder/VGL Mendelian epistasis | **SCORED: concordance 1.000 (62/62) on Sarcidano-2022 real non-circular data** (update below) | rule-integration, scored |

## Why the track is banked (the demonstrated wall)

The label-first admission gate (`off_pathogen_cell_admission_gate_2026-07-01.md`) requires a **measured/observed,
free, no-DUA, per-individual, non-circular** label. The label-first scan + the horse build **demonstrated**
(not assumed) that even the *cleanest* off-pathogen Mendelian rule (horse coat colour) has no such open dataset:
across 5 avenues (Dryad = circular + auth-gated; 709-study = marginal + non-functional 2-SNP coding;
Rieder/Noma = PDF-paywalled; Figshare/Zenodo = none; Arabidopsis = quantitative), none clears the bar. This is
the project's **LABELS wall** at its finest grain. And every clean off-pathogen rule found (IrisPlex, ABO,
coat colour) is a *deployed* published tool — validating it is integration, not new signal.

## Terminal state

- **Frozen validated core** (AMR bacterial/viral, HIV, SARS-CoV-2, TB, PGx) — unchanged, reproducibility-locked.
- **Off-pathogen** — demonstrated at pilot/demo strength; horse rule built but data-walled. **Banked.**
- Everything committed; frozen AMR surface byte-unchanged throughout (leak guard 9/9).

## The only genuine forward paths — both are USER-authority decisions, NOT executor tasks

1. **Acquire a non-public / wet-lab / institutional measured-label source** (clears the label wall by
   construction; e.g. a UK Biobank-class DUA, or a lab-measured cohort). This is the reproducibility-freeze's
   forward-path #1 — a user acquisition decision.
2. **Prospective-lock accrual** — score the frozen decoder against post-lock measured data as it appears
   (`scripts/prospective_lock_validate.py`); passive, no new build.

Absent one of those, the honest recommendation is to **rest the project at this terminal**. The horse cell's
own unblock (a user-supplied joint functional-genotype × observed-colour TSV → `horse_coat_validate.py --data`)
remains available if the user ever wants to close that one cell, but it is low-signal (deployed rule).
