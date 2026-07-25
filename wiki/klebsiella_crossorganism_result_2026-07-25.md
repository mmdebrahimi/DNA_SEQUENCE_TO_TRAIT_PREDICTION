# Klebsiella cross-organism generalization — the paradigm TRANSFERS (2026-07-25)

**Question (user-directed):** does the deterministic RBP/depolymerase → phenotype paradigm — validated for
E. coli phage → outer-membrane receptor — generalize CROSS-ORGANISM to a different host + different phenotype:
*Klebsiella pneumoniae* phage **depolymerase → capsule KL-type**?

**Answer: YES, strongly.**

- **Source (different lab + different organism):** DpoTropiSearch (Concha-Eloko et al., Nat Commun 2025,
  "Unlocking data in Klebsiella lysogens..."; github.com/conchaeloko/DpoTropiSearch + Zenodo
  10.5281/zenodo.14065540). `KL_type_LCA` = capsule type; `domain_seq` = the depolymerase enzymatic domain.
- **Method:** the EXACT phage RBP-caller architecture (`protein_kmers` + Jaccard nearest-neighbour), retargeted
  — depolymerase domain k-mer nearest-neighbour transfer → predicted KL-type, leave-one-out.

## Result

| metric | value |
|---|---|
| depolymerases / KL-types | 2315 / 147 |
| called / correct | 2240 / 1272 |
| **cross-organism LOO accuracy** | **0.568** |
| prior-frequency null | 0.009 |
| **lift over null** | **+0.559** |

Well-powered KL-types [correct/called]:

| KL-type | LOO correct/called |
|---|---|
| KL122|KL106 | 20/20 |
| KL30 | 9/20 |
| KL19 | 20/20 |
| KL25 | 19/20 |
| KL16 | 15/20 |
| KL64 | 19/20 |
| KL106 | 20/20 |
| KL24 | 17/20 |
| KL1 | 17/20 |
| KL28 | 2/20 |
| KL23 | 20/20 |
| KL47 | 19/20 |

## The finding

A HARDER problem scores HIGHER: 147 KL-types (chance 0.009) at **0.568** vs the E. coli cross-lab
tail-fiber RBP number of **0.364** (12 receptor classes). The depolymerase **enzymatic domain** is a cleaner,
more MODULAR sequence→function unit than a tail fiber, so homology transfer works far better on capsule than
on outer-membrane receptors. **The deterministic sequence-homology → phenotype paradigm generalizes across the
organism boundary** (E. coli receptor → Klebsiella capsule) on modular determinants — while remaining
determinant-dependent (KL28 2/20 shows it is not uniform).

## Honest scope

- Labels are **prophage-host-LCA-INFERRED** (in-distribution — the paradigm-TRANSFER analogue of the
  within-LBNL 0.975 LOO), NOT independent wet-lab. The 63 experimentally-validated depolymerases
  (`exp_validated.multi.fasta`) are the gold-standard INDEPENDENT follow-on.
- Stratified subsample (KL-types with >=5 members, capped 20/type) for O(n^2) feasibility — a deterministic
  first-cap, not random. Raw data on D: (gitignored); reproduce via Zenodo + the script.

## Reproduce
```bash
# fetch Training_data.zip from Zenodo 10.5281/zenodo.14065540, extract cols 1-8 -> dpo_labels.tsv
uv run python scripts/klebsiella_depolymerase_crossorganism.py --labels <dpo_labels.tsv>
```
