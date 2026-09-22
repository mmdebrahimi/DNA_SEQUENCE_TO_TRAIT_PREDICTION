# /innovate sweep — 2026-09-21

Framing sweep (H4) over 5 framings / 11 candidates, run through
`skills/soraya/scripts/framing_sweep.py`. Validator-clean, exit 0.
Ledger: `scratchpad/innovate_framings_2026-09-21.json` (session-local).

**Read the caveats section before quoting any survivor.** The engine killed nothing. Every real kill in
this run came from generation-time checking or from the adversarial pass, not from the engine.

## Why a framing sweep rather than a single framing

Self-invoke trigger (b): **≥4 independent approaches failed the SAME way.** BV-BRC continuous MIC
(91% of labels model-imputed), HBV (every free resource an interpretation rule), pathotype (label is the
isolation site), the 19-cell colour family (62% of loci record no causal variant). The standing summary
is "the binding constraint is LABELS, not models" — that is the framing the sweep puts on trial.

## The winning framing: FR2 — question substitution

G9's abduction, and its discriminator is the sharpest artifact of this run. If the project's wins came
from finding **independent witnesses**, an independent witness should MOVE the evidence tier. It doesn't:

| phrase | occurrences in `CLAUDE.md` |
|---|---|
| `TIER DOES NOT MOVE` | 3 |
| `bounds AGREEMENT, never correctness` | 1 |
| `COHERENCE check, NOT a correctness check` | 1 |

**At least 3 of the 5 wins deliberately refused the tier move.** That is inexplicable under a
witness-independence account and is the signature of a different mechanism: the wins **replaced the
blocked question with a weaker one the data on hand can decide**, generating the contrast internally
instead of importing truth. The failures kept the original question and went shopping for the label it
demanded. Artifact names corroborate — failures are named for acquisition
(`bvbrc_strict_mic_4drug_census`, `unscored_genome_label_census`, `fetch_prospective_cohort`), wins for
contrast construction (`tb_implied_boundary`, `mlst_serotype_purity`, `*_concordance`).

Runner-up framings, kept on the record: FR1 witness-independence, FR3 determinant-as-unit-of-evidence,
FR4 perishability/decay-rate, FR5 refusal-as-product.

## Killed

Both killed by the adversarial pass, on verified premise arithmetic.

**`QS-IMPLIED`** — generalize the TB cut-off recovery into a reusable module + a scanner over the
committed datasets. **Killed on the denominator.** `git ls-files data/` returns **58**, `data/raw` **4**,
`data/processed` **0** — the "~248 datasets" figure is the *inventory* (which counts gitignored trees on
D:), not what is tracked. Candidate files with both a continuous measurement and a derived categorical:
**2, and they are the same C. auris cohort.** A scanner over 2 overlapping files is ceremony, and
convergence across datasets is impossible at n=2. Killed independently on coherence: an ECOFF is a
**wild-type** cut-off and **no dataset here ships a WT/NWT column** — every cohort ships (numeric MIC) +
(clinical R/S), so the scanner could only ever re-derive `mic_tiers.py`, which already ships.

**`QS-CONTESTED`** — a "contested-case index" counting which cached genomes would flip under a stated
perturbation. **Already built, and strictly better.** `scripts/determinant_completeness_screen.py` +
`wiki/determinant_completeness_screen_2026-08-31.json` already cover 1,818 genomes × 1,279
unrepresentable families × 6 drugs with per-family label splits, and add a family-wise correction the
proposal never mentions — the correction that dropped 4 of 5 raw hits and kept `rmtE1`. The motivating
example also inverts: **Oxford was never a member of the 1,818** (no `amrfinder_runs/` tree; it ships one
combined file keyed by study guuid), so an index over the 1,818 could not have predicted Oxford's zero —
the probe's whole value was going outside the indexed corpus.

## Survivors — both retargeted, neither as stated

### 1. `QS-TBQUAL` → the 2026-09-21 TB commit is mis-powered under this repo's own standard

**The stated candidate (label-quality monotonicity) is dead** and the adversarial pass consumed the test
to kill it: discordance is **not monotone** (INH MEDIUM 0.0511 sits *below* HIGH 0.0582, and MEDIUM
sensitivity beats HIGH on both drugs); the proposed null is prevalence-confounded (R-rate ranges
0.244→0.638 across tiers, and discordance tracks prevalence mechanically when sens≠spec); and the flag is
**circular for this question** — fraction of isolates within one dilution of the cut-off is RIF HIGH 0.050
/ MED 0.080 / **LOW 0.267**, so quality tracks distance-from-breakpoint, which is the very thing
`tb_implied_boundary` measures.

**What survives is outcome-determinative and verified independently (numbers re-derived here, not taken
from the agent):**

| | RIF | INH |
|---|---|---|
| boundary rung n, as shipped | 188 | 275 |
| boundary rung n, HIGH quality only | **84** | 145 |
| frozen `MIN_BOUNDARY_N` floor | 100 | 100 |
| verdict under HIGH-only | **INDETERMINATE** | PASS |

`scripts/score_tb_cryptic.py:189` and `scripts/score_tb_cryptic_parquet.py:150` both hard-filter to
`PHENOTYPE_QUALITY == "HIGH"`. `scripts/tb_implied_boundary.py` does not filter on quality at all (zero
occurrences of the string). **The frozen bar requires both drugs, so applying the repo's own TB quality
standard downgrades yesterday's `SUPPORTED` to `INDETERMINATE`.**

The two compared strata also differ compositionally, and the artifact does not disclose it: the boundary
rung is **43.6% LOW-quality vs 17.6%** in the lower-S stratum it is compared against (INH 38.5% vs 8.7%)
— a 2.5–4.4× enrichment.

**The effect itself is NOT an artifact — it replicates within HIGH alone:**

| drug | stratum | boundary carriage | lower-S carriage | gap |
|---|---|---|---|---|
| RIF | all quality (as shipped) | 0.378 | 0.0285 | 0.3492 |
| RIF | **HIGH only** | 0.357 | 0.0208 | **0.3363** |
| INH | all quality (as shipped) | 0.160 | 0.0109 | 0.1491 |
| INH | **HIGH only** | 0.152 | 0.0071 | **0.1446** |

So the honest statement is: **the finding holds, the verdict label does not.** ~16× within HIGH for RIF.
This is a powering-gate failure and an undisclosed composition difference, not a refutation.

### 2. `QS-HIVCEIL` → the three published v0.1 OLS gain figures were fit on a right-censored response

**The ceiling is confirmed hard and it is systemic, not NVP-specific.** Across all four Stanford datasets
(~30 drug columns, ~40,000 values): **max = 100.0 everywhere, zero values above 100**, with sparse
irregular values just below (NVP: 87, 89, 90, 91, 91.8, 93 … each n=1–2) then 678 at exactly 100. That is
a right-censor, not a natural mode.

| dataset | worst drug | mode-share at 100 | n at ceiling |
|---|---|---|---|
| NRTI | **3TC** | **0.452** | 831 |
| NNRTI | NVP | 0.330 | 678 |
| INI | EVG | 0.166 | 125 |
| PI | LPV | 0.144 | 260 |

`assay_degeneracy()`'s >25% bar is failed by 3TC and NVP. **The existing censoring guard is vacuous** —
`hiv_nnrti_mutant_catalog.py:234` and `hiv_epistasis.py:27` state that operator-prefixed folds are kept
at their numeric bound, but there are **zero** operator-prefixed values in any column. The guard protects
against something that does not occur while the real censoring, bare `100`, passes straight through. Same
shape as the self-corrected vacuous ResFinder `POINT`-row filter.

**The stated Y181C target is wrong and is dropped:** `hiv_amr.py:60` ships `Y181C/I/A`, so "the OLS
dropped Y181C" is a fact about a derived catalog that was never deployed; and the curation verdict was
driven by blind-spot recovery, defined over isolates carrying *no* catalogued DRM, which Y181C carriers
are not. Censoring also cannot reopen that verdict — `ILLUSTRATIVE_FOLD_CUTOFF = 3.0` ≪ 100, so every
censored observation is R either way and the R/S labels are unaffected.

**Retargeted claim:** the NRTI (+0.06..+0.14 on 5/6), PI (+0.056 mean, 8/8) and INSTI (+0.087 mean, 5/5)
v0.1 gains were each selected by thresholding a multivariate-OLS coefficient at `≥ log10(1.5)`, fit on a
response right-censored at 100 on drugs where up to **45%** of observations sit at the ceiling.
Right-censoring attenuates coefficients for the *strongest* DRMs specifically. Those three figures are
published headlines in `CLAUDE.md`. A Tobit refit is cheap, offline and decidable.

**Scope rail:** no v0.1 catalog is deployed (`hiv_amr.py` routes NRTI/PI/INSTI through position-based
classes). This is a **methodological audit of published numbers, not a decoder change**, and must be
stated that way.

## Two defects found in passing

**`recover_cutoff` accepts a degenerate all-susceptible ladder.** Reproduced:
`recover_cutoff([-2,-1,0], ["S","S","S"]) -> 0.0`. With no resistance anywhere there is no threshold, but
`out` is set to the top rung and the `above` branch is skipped because the slice is empty, so it returns a
vacuously "pure" cut-off from an information-free column. Harmless for the committed run (both classes
present); it would have poisoned the killed scanner. Third defect in that script.

**`CLAUDE.md:392` is falsified by the repo's own artifact.** It states the ECOFF values are "every value
UNSOURCED" and that EUCAST publishes "no stable per-result URL". `dna_decode/data/ecoff_catalog.py` now
carries **four sourced ECOFFs**, each with `source_url = https://mic.eucast.org/search/diagram/<id>` and a
verbatim quote. Stable per-result URLs exist. The fabrication rail was never breached — the values
arrived *through* it, with locator and quote. What was wrong was declaring the fetch impossible.

## Outcome — both survivors were acted on the same day (2026-09-22)

| survivor | shipped as | result |
|---|---|---|
| `QS-TBQUAL` | `275d5bd` | the 2026-09-21 TB artifact now carries a namespace-separate `quality_gate`: `SUPPORTED` under the frozen bar, `INDETERMINATE` under HIGH-only labels, with the within-HIGH replication beside it. Plus a third defect in `recover_cutoff` (vacuous cut-off on an all-susceptible ladder). 6 new tests. |
| `QS-HIVCEIL` | `d6c83b5` | `scripts/hiv_fold_censoring_audit.py` + memo: one upper bound of 100.0 across all four datasets, nothing above it anywhere, 3TC 0.452 / NVP 0.330 failing the repo's own bar, and the existing operator-prefix guard measured **vacuous** (zero such values). 11 tests. |

Both were **retargeted before shipping** — neither shipped in the form the engine marked `survived`.
That gap between "survived a cheap kill-test" and "survived an attack on its substance" is the whole
reason the expensive pass is not optional.

## Honest caveats on this sweep

- **The engine killed nothing: 11 survived, 0 killed.** Every kill-test was a "does this artifact already
  exist" predicate that I had pre-run by hand before writing the ledger, so the engine's execution was
  confirmatory, not discriminating. **These kill-tests establish novelty-in-repo, never soundness.** The
  survivor list out of the engine was worth roughly nothing on its own; the adversarial pass then killed
  2 of the 4 top survivors on substance and forced both remaining ones to change target.
- **`test-exit-0` kill-tests cannot carry discrimination controls** (H1 is read-only-kinds only), so each
  grep was paired by hand with a positive control on the same search path proving the path was live — a
  typo'd path would otherwise have faked a survival on every one of them.
- **Four candidates were killed at generation time, before reaching the ledger**, by checking git: cut-off
  recovery × CRyPTIC (`56da283`, done yesterday), `assay_degeneracy` × MIC (`28189ac`), diversity bar ×
  AR Bank (`scripts/ar_bank_source_concentration.py`), lineage/doubt layers × non-AMR routes (`acaaf01`,
  a closed negative). VME/ME was likewise killed — `dna_decode/eval/error_rates.py` already ships it,
  which retires evidence-surface candidate action #4.
- **`QS-TBQUAL` is the powered version of a lead the 2026-09-10 sweep raised at n=1** and labelled "a
  lead, not a finding" (`ERRORS-isolate-clustered`: error mass may be data quality, not catalog gaps).
  Independent route, different axis (shipped label quality vs inferred genome quality), 3,144 isolates
  instead of one.
- **Deeper adversaries not run.** `/brainstorm` and `/idea-validation-council` are user-only. The two
  survivors are `unfalsified-until-council` on their *retargeted* forms — the adversarial pass attacked
  the forms as originally stated.
- **FR2's discriminator is a grep over prose**, not an experiment. It is strong evidence about how the
  project describes its own wins; it is not proof of the mechanism.
