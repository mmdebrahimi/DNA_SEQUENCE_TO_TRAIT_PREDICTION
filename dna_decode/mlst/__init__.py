"""MLST (multi-locus sequence typing) decoder — deterministic, PubMLST-backed.

7-gene sequence typing: for each housekeeping locus, find the EXACT allele in the genome (blastn 100%
identity + 100% coverage via the shared engine), assemble the allele-number profile, and look up the
Sequence Type (ST) in the PubMLST profiles table. v0 scope: E. coli Achtman scheme (adk/fumC/gyrB/icd/mdh/
purA/recA). Differs from the presence/codon decoders: exact-allele match + profile->ST lookup. Offline-safe.
"""
