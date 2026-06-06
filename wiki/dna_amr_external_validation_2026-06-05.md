# dna-amr external validation — held-out cohort — 2026-06-05

> Recommendation-#1 deliverable: confirm the deterministic AMR decoder's op-chars hold OUT-OF-COHORT
> (the in-cohort 0.939 was tuned on N=147 — threshold=2 picked there). Turns "a benchmark" into "a
> trustworthy tool."

## Result

| Set | N | accuracy | sensitivity (R) | specificity (S) |
|---|---|---|---|---|
| in-cohort (tuning) N=147 | 147 | 0.939 | 0.931 | 0.947 |
| **held-out (external) N=29** | 29 | **0.862** | **0.882** | **0.833** |

Confusion (held-out, threshold=2): TP 15 · FP 2 · TN 10 · FN 2 (0 indeterminate).

**Threshold choice validated on independent data:** held-out threshold=1 → acc 0.793 / spec 0.583
(over-calls); threshold=2 → acc 0.862 / spec 0.833. The ≥2-determinant rule is better out-of-cohort too,
not just on the tuning set.

## Interpretation
- **The tool generalizes.** Held-out accuracy 0.862 (sens 0.882 / spec 0.833) is well above chance and
  close to in-cohort — the deterministic `≥2 curated quinolone determinants → R` rule is not overfit to
  the N=147 cohort.
- **The in-cohort 0.939 was mildly optimistic (~8 pp).** Expected: threshold=2 was picked on N=147, and
  N=29 carries wide CIs. The honest headline op-chars for the shipped tool are now the held-out numbers.
- **Sensitivity stays high (0.882)** — the clinically safer error direction (few resistant strains called
  susceptible: 2/17 FN).

## Held-out set + honest caveat
29 cipro-labeled E. coli strains (17R/12S) from the `gate_b_cohort` (a separate, cef-focused BV-BRC AST
cohort), excluding any strain in the N=147 cipro tuning cohort. Genuinely **out-of-cohort** (not used to
pick threshold=2). **Caveat:** same source database (BV-BRC) + same broth-microdilution label methodology,
so this is an out-of-cohort holdout, NOT a cross-source / different-lab validation. A fully independent
check (different lab / NARMS / EUCAST-sourced labels) is a future step. N=29 is modest (wide CIs).

## Reproduce
```
# held-out = gate_b_cohort cipro-labeled strains not in stage2_n150_cipro_cohort; AMRFinder-cached in data/amrfinder_runs
uv run python -c "from dna_decode.data.cohort import load_cohort; from dna_decode.eval.amr_rules import evaluate_cohort; \
n147={s.assembly_accession for s in load_cohort('data/processed/stage2_n150_cipro_cohort.parquet').strains}; \
gb=load_cohort('data/processed/gate_b_cohort.parquet'); \
held=[s for s in gb.strains if 'ciprofloxacin' in s.ast_labels and s.assembly_accession not in n147]; \
print(evaluate_cohort('data/amrfinder_runs',[(s.assembly_accession,int(s.ast_labels['ciprofloxacin'])) for s in held],'ciprofloxacin'))"
```

## Provenance
`dna_decode/eval/amr_rules.py` (the caller), `data/amrfinder_runs/` (per-strain AMRFinder outputs, gitignored).
Per /soraya recommendation #1 (external-validate the shipped decoder before further scope).
