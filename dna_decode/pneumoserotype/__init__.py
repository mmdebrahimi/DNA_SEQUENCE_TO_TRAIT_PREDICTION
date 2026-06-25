"""Streptococcus pneumoniae capsular serotype decoder (PneumoCaT/SeroBA-style, deterministic).

A trait decoder on the shared curated-DB blastn engine (dna_decode.typing): BLASTs a curated per-serotype
capsular (cps) reference DB vs a genome and calls the serotype of the best-matching reference. Sibling of
dna-serotype (E. coli O:H) and dna-ktype (Klebsiella capsule). Faithful to the cps-reference typing method;
NOT an independent baseline. Offline-safe.
"""
