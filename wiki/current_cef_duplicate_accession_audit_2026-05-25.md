# Cipro duplicate accession audit

- Cohort file: `data\processed\gate_b_cohort.parquet`
- Drug: `ceftriaxone`
- Pool size: `50`
- Duplicate accessions in pool: `0`
- Strains participating in duplicate accessions: `0`
- LOSO leakage present: `False`

## Model bundle

- Model path: `data\processed\models\ceftriaxone_nucleotide_transformer.pkl`
- Drug: `ceftriaxone`
- Model name: `nucleotide_transformer`
- CV strategy: `loso`
- CV grouping: `strain_id`
- Training cohort: `gate_b_cohort`
- Trained on: `2026-05-26`
- Primary CV AUROC: `0.895`
- Legacy auroc_loso field: `0.895`
- n_strains: `49`

## Duplicate accessions

(none)

## Verdict

PASS: no duplicated non-empty assembly_accession values found in the cipro pool.
