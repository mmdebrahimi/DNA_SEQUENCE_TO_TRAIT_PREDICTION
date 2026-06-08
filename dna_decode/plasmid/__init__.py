"""Plasmid replicon typing decoder (PlasmidFinder-style, deterministic).

A new trait decoder sibling to `dna_decode.pathotype` and `dna_decode.amr` — BLASTs the PlasmidFinder
replicon allele DB against a genome assembly and reports the plasmid incompatibility (Inc) replicons present.
Composes with the AMR decoder: AMR says WHAT resistance; plasmid typing says whether it likely rides a
known mobile element. Deterministic curated-DB caller (same method class as VirulenceFinder pathotype), NOT
embeddings. Reuses the pathotype blastn machinery. Offline-safe (degrades to status=unavailable).
"""
