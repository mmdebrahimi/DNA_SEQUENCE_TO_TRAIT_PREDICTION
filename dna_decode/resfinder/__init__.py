"""Acquired-AMR-gene detector (ResFinder-style, deterministic).

A new trait decoder on the shared curated-DB blastn engine (dna_decode.typing): BLASTs the ResFinder
acquired-resistance allele DB vs a genome and reports the acquired AMR genes present, grouped by antibiotic
class. Deliberately a SECOND, INDEPENDENT caller vs the AMRFinder-based dna-amr engine (different curated DB)
— a cross-tool concordance check the AMR decoder otherwise lacked. Offline-safe.
"""
