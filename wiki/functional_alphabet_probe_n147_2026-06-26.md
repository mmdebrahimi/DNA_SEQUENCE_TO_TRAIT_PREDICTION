# Functional-alphabet probe (2026-06-26) -- stage2_n150_cipro_cohort.parquet

**Verdict: TIES**

- functional within-lineage concordance: 1.0
- k-mer within-lineage concordance:      0.7209302325581395
- gap (functional - k-mer):              0.2790697674418605
- paired in-MLST permutation gap p:      0.0565
- gap null 95% band:                     [-0.34883720930232553, 0.313953488372093]
- powering: n_shared_lineages=6, n_pairs=43
- admitted N=147 (R=72/S=75); dropped: 0 no-AMRFinder, 0 no-FASTA
- vocab: functional=243 tokens, k-mer=10000

_functional within-lineage concordance EXCEEDS k-mer with paired in-MLST label-permutation gap p < 0.05; UNDERPOWERED when n_pairs < 10. Any bucket is a valid honest completion._

_Non-neural, CPU-only. A within-lineage win means 'beats k-mer within-lineage', NOT 'mechanism proven'._