# Functional-alphabet probe (2026-06-26) -- shared_lineage_tetracycline_cohort.parquet

**Verdict: BEATS_KMER**

- functional within-lineage concordance: 0.9626436781609196
- k-mer within-lineage concordance:      0.5402298850574713
- gap (functional - k-mer):              0.4224137931034483
- paired in-MLST permutation gap p:      0.0
- gap null 95% band:                     [-0.14367816091954017, 0.15237068965517203]
- powering: n_shared_lineages=20, n_pairs=174
- admitted N=118 (R=58/S=60); dropped: 0 no-AMRFinder, 0 no-FASTA
- vocab: functional=104 tokens, k-mer=10000

_functional within-lineage concordance EXCEEDS k-mer with paired in-MLST label-permutation gap p < 0.05; UNDERPOWERED when n_pairs < 10. Any bucket is a valid honest completion._

_Non-neural, CPU-only. A within-lineage win means 'beats k-mer within-lineage', NOT 'mechanism proven'._