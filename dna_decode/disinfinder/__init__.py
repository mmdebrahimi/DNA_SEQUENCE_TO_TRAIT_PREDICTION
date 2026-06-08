"""Biocide / disinfectant resistance gene detector (DisinFinder-style, deterministic).

A new trait decoder on the shared curated-DB blastn engine (dna_decode.typing): BLASTs the DisinFinder
allele DB (qacA/qacB/formA/... quaternary-ammonium + formaldehyde resistance genes) vs a genome and reports
the acquired biocide-resistance genes present. Hospital infection-control relevant; qac genes frequently
ride the same plasmids as AMR (composes with dna-coloc). Sibling of dna-resfinder. Offline-safe.
"""
