"""E. coli O:H serotype decoder (SerotypeFinder-style, deterministic).

A new trait decoder on the shared curated-DB blastn engine (dna_decode.typing): BLASTs the SerotypeFinder
O-antigen (wzx/wzy/wzm/wzt) + H-antigen (fliC/...) allele DB vs a genome and reports the O:H serotype.
Sibling of dna-amr / dna-pathotype / dna-plasmid. Offline-safe.
"""
