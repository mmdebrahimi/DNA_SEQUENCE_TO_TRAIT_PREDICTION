# Functional-alphabet probe (2026-06-27) -- shared_lineage_gentamicin_cohort.parquet

**Verdict: FAILS**

- functional within-lineage concordance: 0.5944444444444444
- k-mer within-lineage concordance:      0.6666666666666666
- gap (functional - k-mer):              -0.07222222222222219
- paired in-MLST permutation gap p:      0.865
- gap null 95% band:                     [-0.1305555555555556, 0.13895833333333296]
- powering: n_shared_lineages=20, n_pairs=180
- admitted N=120 (R=60/S=60); dropped: 0 no-AMRFinder, 0 no-FASTA
- vocab: functional=136 tokens, k-mer=10000

_functional within-lineage concordance EXCEEDS k-mer with paired in-MLST label-permutation gap p < 0.05; UNDERPOWERED when n_pairs < 10. Any bucket is a valid honest completion._

_Non-neural, CPU-only. A within-lineage win means 'beats k-mer within-lineage', NOT 'mechanism proven'._