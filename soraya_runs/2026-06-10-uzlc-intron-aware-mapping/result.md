# Soraya --advance #1 — result (intron-aware multi-HSP codon mapping)
- Generalized observed_substitutions to multi-HSP stitching (query-position-keyed) -> intron-aware.
- Improves every genome-mode caller (fungal/K13) for multi-exon/split genes; unblocks pfcrt-class targets.
- Validated on real K13 CDS split into exons (exon1 + deep-exon2 across intron; WT boundary = no false call).
- 2 tests; full suite 967 passed, 0 regressions (single-exon path unchanged).
- pfcrt genome-mode flip pends a committed pfcrt CDS reference (intron blocker removed; ref is the gate).
## Next
1. Commit a pfcrt CDS reference (PlasmoDB/3D7 exon-stitch) -> flip pfcrt genome mode + real genomic-allele validation.
2. pfmdr1 (directional catalog) / antiviral kingdom.
3. G2 (user-owned, Databricks).
