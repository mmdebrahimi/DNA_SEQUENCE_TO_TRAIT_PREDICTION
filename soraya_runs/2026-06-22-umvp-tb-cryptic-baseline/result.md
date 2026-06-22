# result — tb-cryptic-baseline  (verdict: mvp-reached)
Computed the TB cell's in-distribution baseline (RIF + INH) from the CRyPTIC Zenodo parquet on D:, via a
thin parquet→calls adapter over the FROZEN scorer (scripts/score_tb_cryptic_parquet.py). The TB cell ran
on real data for the first time (was data-blocked on ~1.6 TB regeno).
- RIF: raw 0.916/0.974, lineage-collapsed 0.41/0.991 (n=8955).
- INH: raw 0.889/0.989, lineage-collapsed 0.349/0.995 (n=9518).
Status WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE. 5 offline tests; full suite 1593; frozen surface +
TB leak-guard byte-unchanged. Recovery rounds: 1 (diagnosed the SNV-only parser vs codon-MNV via
snv_components — confirmed codon determinants ARE matched; raw sens 0.916 proves it).
Verify-in-batch: numbers are sane (RIF catalogue ~0.92 sens is the known ballpark); the raw→lineage drop is
the documented clonality inflation.
