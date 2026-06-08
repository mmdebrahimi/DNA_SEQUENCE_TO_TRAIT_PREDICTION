"""Chromosomal AMR point-mutation decoder (PointFinder-style, deterministic).

BLASTs each PointFinder reference gene CDS (gyrA/parC/gyrB/parE for E. coli QRDR) vs a genome, maps the
subject amino acid at each catalogued codon position (gap-aware, via typing.codon_map), and calls a
resistance mutation when the subject AA matches a Res_codon in resistens-overview.txt. An INDEPENDENT
chromosomal point-mutation caller — complements dna-amr (AMRFinder POINT) + dna-resfinder (acquired only,
no point mutations). v0 scope: E. coli fluoroquinolone QRDR. Offline-safe.
"""
