# Functional-alphabet probe (2026-06-27) -- shared_lineage_ceftriaxone_cohort.parquet

**Verdict: BEATS_KMER**

- functional within-lineage concordance: 0.9388888888888889
- k-mer within-lineage concordance:      0.6055555555555555
- gap (functional - k-mer):              0.33333333333333337
- paired in-MLST permutation gap p:      0.0
- gap null 95% band:                     [-0.1388888888888889, 0.14166666666666666]
- powering: n_shared_lineages=20, n_pairs=180
- admitted N=120 (R=60/S=60); dropped: 0 no-AMRFinder, 0 no-FASTA
- vocab: functional=202 tokens, k-mer=10000

_functional within-lineage concordance EXCEEDS k-mer with paired in-MLST label-permutation gap p < 0.05; UNDERPOWERED when n_pairs < 10. Any bucket is a valid honest completion._

_Non-neural, CPU-only. A within-lineage win means 'beats k-mer within-lineage', NOT 'mechanism proven'._