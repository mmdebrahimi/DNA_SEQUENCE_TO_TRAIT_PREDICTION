# Plan: using the ChatGPT genome→phenotype catalog (2026-07-02)

Screened through the F1 8-gate scorecard (`wiki/chatgpt_catalog_screen_2026-07-02.md`). **Verdict: the catalog
is useful, but most of it is user-gated or circular-as-label.** Exactly one new source is free + measured +
joinable + deep + non-circular: **DepMap/CCLE**. The catalog's own headline warning (label leakage /
circularity) is this project's central wall, independently confirmed (ClinVar type-1 circularity). This plan
routes the useful sources into 4 concrete tracks by usage class.

## Track A — DepMap/CCLE pilot + decoder (the executor-actionable win; the cancer analog of yeast)
DepMap is the human-cell-line twin of the 1002 Yeast Genomes substrate: cell-line **genotype** (CCLE
mutations/expression, >1,000 lines) × **measured phenotype** (PRISM drug sensitivity 4,518 compounds; CRISPR
dependency), joined by `depmap_id`, all free CSVs from depmap.org. **Do it exactly like the yeast run:**
1. **Pilot-fetch** (F4-style): PRISM `Repurposing_Public_*_Data_Matrix.csv` + `OmicsSomaticMutations`/
   `CCLE_expression` → confirm the `depmap_id` join + a signal-vs-null on one drug (like fluconazole on yeast).
2. **De-confounded decoder** (the yeast pattern, the ONLY honest test): genotype → drug-response Ridge/XGB,
   scored by naive / **lineage-de-confounded within-lineage r²** (lineage = tissue-of-origin / cell-line
   lineage, the DepMap analog of yeast clades) + permutation null + finer-lineage robustness. Beat a
   lineage-only baseline within-lineage or it's the 5th de-confounded negative (honest either way).
3. **Attribution** (learn from the yeast capstone): for a known pharmacogenomic pair (e.g. BRAF-mut ×
   vemurafenib, EGFR × erlotinib), check the top within-lineage genes ARE the known biomarker — this substrate
   HAS canonical gene-level (not copy-number) mechanisms, so attribution should succeed where yeast's didn't.
Gate: free/network only; no money; frozen AMR surface untouched.

## Track B — deterministic-catalog EXTENSION (ClinVar / ClinGen / CIViC) — the validated paradigm, done right
These are curated determinant→phenotype catalogs — the HUMAN analog of the WHO-TB / CARD catalogs the AMR
decoder already uses. **They are NOT labels** (type-1 circularity: assertions derived from the in-silico
predictors a decoder competes with). Use them the way the project uses WHO-TB:
1. Ingest as a deterministic RULE catalog (variant → pathogenicity/actionability), branded a
   `*_KNOWLEDGE_BASELINE` exactly like the TB CRyPTIC number.
2. **Validate against a TEMPORAL HOLDOUT** — score only on variants deposited AFTER the catalog/predictor
   training cutoff (the prospective-lock discipline the project already built), never on the catalog's own
   assertions (the wrapper-vs-tool + prospective-lock lessons). CIViC is CC0 (cleanest license).
Gate: free; the honesty rail (KNOWLEDGE_BASELINE + temporal holdout) is mandatory, not optional.

## Track C — USER-AUTHORITY acquisition (All of Us / UK Biobank / FinnGen / dbGaP / EGA) — freeze forward-path #1
These fail G1 by design (controlled access: researcher agreement / DAC approval / institutional affiliation /
possible cost). They are the highest supervised-training value AND they are **the "acquire a non-public label
source" path the reproducibility freeze already named as forward-path #1** — a USER decision, NOT an executor
task. Draft asks for the user to action:
- **All of Us Researcher Workbench** — registered tier is the lowest-friction (US-affiliation + training);
  gives WGS + EHR + surveys. Best first application.
- **UK Biobank** — approved-researcher application + data-use fee; the default common-disease corpus.
- **FinnGen** — public summary stats are free NOW (usable without access); individual-level is controlled.
- **dbGaP / EGA** — per-study controlled access; treat as "dataset multipliers" AFTER a first model.
Soraya can draft the applications; the user submits/authorizes (identity + institutional + any fee = authority).

## Track D — human-variant FEATURE sources + benchmarks (only if the project pivots to human variant-effect)
GTEx / ENCODE / SCREEN / eQTLGen / 4DN / FANTOM5 / Factorbook / VEP / CADD are FEATURE generators for a human
noncoding-variant-effect model, not genotype→organismal-phenotype label substrates — out of scope unless the
project deliberately opens a human-variant-effect arm. If it does, adopt GIAB + CAGI + precisionFDA as the
anti-circularity evaluation harness (aligns with the project's existing leakage discipline).

## Recommended sequencing
1. **Track A (DepMap pilot + de-confounded decoder)** — the one free, executor-actionable, high-VOI move;
   directly transfers the yeast machinery; a plausible FIRST attribution win. Do next.
2. **Track C (draft the All of Us / UK Biobank access asks)** — surface for the USER; the real unlock for
   deep supervised human phenotype learning, but their authority.
3. **Track B (ClinVar/CIViC as validated catalogs)** — medium-VOI; only with temporal-holdout honesty.
4. **Track D** — parked unless a human-variant-effect arm is deliberately opened.

## Non-goals / honest bounds
- Do NOT treat ClinVar/CIViC/CADD as ground-truth labels (circular). Features/catalogs only, with independent
  validation.
- Do NOT force feature sources (GTEx/ENCODE) into a phenotype-label abstraction.
- The frozen AMR surface stays byte-unchanged across all tracks.
