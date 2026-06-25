"""Salmonella enterica serovar decoder (SeqSero2-style, deterministic).

A trait decoder on the shared curated-DB blastn engine (dna_decode.typing): BLASTs the Salmonella
antigen allele DB (O-antigen group genes + H1 = fliC + H2 = fljB) vs a genome, assembles the
White-Kauffmann-Le Minor antigenic formula (O:H1:H2), and looks up the serovar. Sibling of dna-serotype
(O:H) and dna-ktype (capsule). Faithful to the SeqSero2/Kauffmann-White method; NOT an independent
baseline. Offline-safe.
"""
